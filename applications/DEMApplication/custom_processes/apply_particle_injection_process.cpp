//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt

// Application includes
#include "apply_particle_injection_process.hpp"
#include "utilities/parallel_utilities.h"
#include "utilities/function_parser_utility.h"

namespace Kratos
{


    /// Constructor
    ApplyParticleInjectionProcess::ApplyParticleInjectionProcess(
        ModelPart& rModelPart,
        Parameters rParameters
        ) : Process(Flags()) , mrModelPart(rModelPart), mParameters(rParameters), mInterval(rParameters)
    {
        KRATOS_TRY

        Parameters default_parameters( R"(
            {
                "model_part_name":"MODEL_PART_NAME",
                "injected_particle_radius" : 1.0,
                "injected_particle_velocity": [10.0, "3*t", "x+y"],
                "max_rand_deviation_angle": 1.0e-5,
                "injection_start" : 0,
                "injection_stop" : 0.0,
                "injector_element_type" : 'SphericParticle',
                "injected_element_type" : 'SphericParticle',
                "contains_clusters" : "false",
                "injected_number_of_particles": 1000,
                "imposed_mass_flow_option" : false,
                "mass_flow": 100.0,
                "probability_distribution": normal,
                "standard_deviation": 0.0,
                "dense_inlet_option" : false,
            }  )" );



        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);
        mParameters = rParameters;

        // This list of member variables could be avoided by integrating them all into mParameters
        // mParameters["injected_particle_radius"].GetDouble()

        // mInjectedRadius = rParameters["injected_particle_radius"].GetDouble();
        // mInjectedVelocity = rParameters["injected_particle_velocity"].GetVector();
        // mMaxRandDevAngle = rParameters["max_rand_deviation_angle"].GetDouble();
        // mInjectorType = rParameters["injector_element_type"].GetString();
        // mInjectedType = rParameters["injected_element_type"].GetString();
        // mContainsClusters = rParameters["contains_clusters"].GetBool();
        // mInjectedNumberOfParticles = rParameters["injected_number_of_particles"].GetInt();
        // mImposedMassFlowOpt = rParameters["imposed_mass_flow_option"].GetBool();
        // mMassFlowValue = rParameters["mass_flow"].GetDouble();
        // mProbDistrib = rParameters["probability_distribution"].GetString();
        // mStandardDev = rParameters["standard_deviation"].GetDouble();
        // mDenseInletOption = rParameters["dense_inlet_option"].GetBool();


        KRATOS_CATCH("");
    }


///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    //********************* moved out of inlet.cpp
    DEM_Inlet::DEM_Inlet(ModelPart& inlet_modelpart, const Parameters& r_inlet_settings, const int seed):
     mInletModelPart(inlet_modelpart), mInletsSettings(Parameters(r_inlet_settings))
        {
        const int number_of_submodelparts = inlet_modelpart.NumberOfSubModelParts();
        mPartialParticleToInsert.resize(number_of_submodelparts);
        mLastInjectionTimes.resize(number_of_submodelparts);
        //mTotalNumberOfDetachedParticles.resize(number_of_submodelparts, false);
        mLayerRemoved.resize(number_of_submodelparts);
        mNumberOfParticlesInjected.resize(number_of_submodelparts);
        mMassInjected.resize(number_of_submodelparts);

        std::mt19937 gen(seed);
        mGenerator = gen;

        int smp_iterator_number = 0;
        for (ModelPart::SubModelPartsContainerType::iterator sub_model_part = inlet_modelpart.SubModelPartsBegin(); sub_model_part != inlet_modelpart.SubModelPartsEnd(); ++sub_model_part) {
            mPartialParticleToInsert[smp_iterator_number] = 0.0;
            mLastInjectionTimes[smp_iterator_number] = 0.0;
            mLayerRemoved[smp_iterator_number] = false;
            //mTotalNumberOfDetachedParticles[smp_iterator_number] = 0.0;
            mNumberOfParticlesInjected[smp_iterator_number] = 0;
            mMassInjected[smp_iterator_number] = 0.0;
            smp_iterator_number++;
        }

        mFirstInjectionIsDone = false;
        mBallsModelPartHasSphericity = false;
        mBallsModelPartHasRotation   = false;
        mTotalNumberOfParticlesInjected = 0;
        mTotalMassInjected = 0.0;
        SetNormalizedMaxIndentationForRelease(0.0);
        SetNormalizedMaxIndentationForNewParticleCreation(0.0);

        mWarningTooSmallInlet = false;
        mWarningTooSmallInletForMassFlow = false;
    }


    double DEM_Inlet::SetDistributionMeanRadius(ModelPart& mp) {
        // to process - migrate to lambda
        double mean_radius = 0.0;

        if (mp[PROBABILITY_DISTRIBUTION] == "piecewise_linear" || mp[PROBABILITY_DISTRIBUTION] == "discrete"){
            mean_radius = mInletsRandomVariables[mp.Name()]->GetMean();
        }

        else {
            mean_radius = mp[RADIUS];
        }

        return mean_radius;
    }

    double DEM_Inlet::SetMaxDistributionRadius(ModelPart& mp) {
        // to process - migrate to lambda
        return 1.5 * GetMaxRadius(mp);
    }

    double DEM_Inlet::GetMaxRadius(ModelPart& mp){
        // to process - migrate to lambda - reconsider for tests
        double max_radius = 0.0;

        if (mp[PROBABILITY_DISTRIBUTION] == "piecewise_linear" || mp[PROBABILITY_DISTRIBUTION] == "discrete"){
            const array_1d<double, 2>& support = mInletsRandomVariables[mp.Name()]->GetSupport();
            max_radius = support[1];
        }

        else {
            max_radius = mp[RADIUS];
        }

        return max_radius;
    }

    void DEM_Inlet::InitializeDEM_Inlet(ModelPart& r_modelpart, ParticleCreatorDestructor& creator, const bool using_strategy_for_continuum) {
        // to process

        mStrategyForContinuum = using_strategy_for_continuum;
        unsigned int& max_Id=creator.mMaxNodeId;
        //CreatePropertiesProxies(mFastProperties, mInletModelPart);
        mFastProperties = PropertiesProxiesManager().GetPropertiesProxies(r_modelpart);
        VariablesList r_modelpart_nodal_variables_list = r_modelpart.GetNodalSolutionStepVariablesList();

        if (r_modelpart_nodal_variables_list.Has(PARTICLE_SPHERICITY)) mBallsModelPartHasSphericity = true;

        if (r_modelpart.GetProcessInfo()[ROTATION_OPTION]) {
            mBallsModelPartHasRotation = true;
            mInletModelPart.GetProcessInfo()[ROTATION_OPTION] = true;
        }
        else {
            mInletModelPart.GetProcessInfo()[ROTATION_OPTION] = false;
        }

        int smp_number = 0;

        for (ModelPart::SubModelPartsContainerType::iterator smp_it = mInletModelPart.SubModelPartsBegin(); smp_it != mInletModelPart.SubModelPartsEnd(); ++smp_it) {
            ModelPart* mp = &*smp_it;
            mListOfSubModelParts.push_back(mp);
        }
        std::sort(mListOfSubModelParts.begin(), mListOfSubModelParts.end(), SortSubModelPartsByName);

        for(int i=0; i<(int)mListOfSubModelParts.size(); i++) {

            ModelPart& mp = *mListOfSubModelParts[i];

            CheckSubModelPart(mp);

            int mesh_size = mp.NumberOfNodes();
            if (!mesh_size) continue;
            ModelPart::NodesContainerType::ContainerType& all_nodes = mp.NodesArray();
            mp[INLET_INITIAL_VELOCITY] = mp[LINEAR_VELOCITY];    //This is the velocity of the moving injector of particles
            mp[INLET_INITIAL_PARTICLES_VELOCITY] = mp[VELOCITY]; //This is the initial velocity vector of the injected particles

            array_1d<double, 3>& inlet_velocity = mp[VELOCITY];

            KRATOS_ERROR_IF((inlet_velocity[0] == 0.0) && (inlet_velocity[1] == 0.0) && (inlet_velocity[2] == 0.0)) << "The inlet velocity cannot be zero for inlet: " << mp.Name() << std::endl;

            double max_rand_dev_angle = mp[MAX_RAND_DEVIATION_ANGLE];

            KRATOS_ERROR_IF(max_rand_dev_angle < 0.0 || max_rand_dev_angle > 89.5) << "The velocity deviation angle must be between 0 and 89.5 degrees for inlet: "<< mp.Name() << std::endl;

            Properties::Pointer p_properties = r_modelpart.pGetProperties(mp[PROPERTIES_ID]);
            int general_properties_id = p_properties->Id();

            PropertiesProxy* p_fast_properties = NULL;

            for (unsigned int i = 0; i < mFastProperties.size(); i++) {
                int fast_properties_id = mFastProperties[i].GetId();
                if (fast_properties_id == general_properties_id) {
                    p_fast_properties = &(mFastProperties[i]);
                    break;
                }
                mLastInjectionTimes[smp_number] = mp[INLET_START_TIME];
            }

            if (mp[PROBABILITY_DISTRIBUTION] == "piecewise_linear" || mp[PROBABILITY_DISTRIBUTION] == "discrete"){
                if (!mInletsSettings.Has(mp.Name())){
                    KRATOS_ERROR << "dem_inlet_settings does not contain settings for the inlet" << mp.Name() << ". Please, provide them.";
                }
                const Parameters& inlet_settings = mInletsSettings[mp.Name()];
                mInletsRandomSettings.emplace(mp.Name(), inlet_settings["random_variable_settings"]);
                const Parameters& rv_settings = mInletsRandomSettings[mp.Name()];
                int seed = rv_settings["seed"].GetInt();
                if (!rv_settings["do_use_seed"].GetBool()){
                    seed = std::random_device{}();
                }
                if (mp[PROBABILITY_DISTRIBUTION] == "piecewise_linear"){

                    mInletsRandomVariables[mp.Name()] = std::unique_ptr<PiecewiseLinearRandomVariable>(new PiecewiseLinearRandomVariable(rv_settings, seed));
                }

                else if (mp[PROBABILITY_DISTRIBUTION] == "discrete"){
                    mInletsRandomVariables[mp.Name()] = std::unique_ptr<DiscreteRandomVariable>(new DiscreteRandomVariable(rv_settings, seed));
                }

                else {
                    KRATOS_ERROR << "Unknown DEM inlet random variable: " << mp[PROBABILITY_DISTRIBUTION] << ".";
                }
            }

            double max_radius = SetMaxDistributionRadius(mp);

            if (!mp[MINIMUM_RADIUS]) {
                mp[MINIMUM_RADIUS] = 0.5 * mp[RADIUS];
            }
            if (!mp[MAXIMUM_RADIUS]) {
                mp[MAXIMUM_RADIUS] = max_radius;
            }

            Element::Pointer dummy_element_pointer;
            std::string& ElementNameString = mp[INJECTOR_ELEMENT_TYPE];
            const Element& r_reference_element = KratosComponents<Element>::Get(ElementNameString);

            for (int i = 0; i < mesh_size; i++) {
                Element* p_element = creator.ElementCreatorWithPhysicalParameters(r_modelpart,
                                                                                max_Id+1,
                                                                                all_nodes[i],
                                                                                dummy_element_pointer,
                                                                                p_properties,
                                                                                mp,
                                                                                mInletsRandomVariables,
                                                                                r_reference_element,
                                                                                p_fast_properties,
                                                                                mBallsModelPartHasSphericity,
                                                                                mBallsModelPartHasRotation,
                                                                                true,
                                                                                mp.Elements());

                FixInjectorConditions(p_element);
                max_Id++;
                /*if(mStrategyForContinuum){
                    SphericContinuumParticle* p_continuum_spheric_particle = dynamic_cast<SphericContinuumParticle*>(p_element);
                    p_continuum_spheric_particle->mContinuumInitialNeighborsSize=0;
                    p_continuum_spheric_particle->mInitialNeighborsSize=0;
                }*/
            }
            smp_number++;
        } //for smp_it
    } //InitializeDEM_Inlet

    void DEM_Inlet::DettachElements(ModelPart& r_modelpart, unsigned int& max_Id) {
        // to process

        ProcessInfo& r_process_info = r_modelpart.GetProcessInfo();

        ///DIMENSION
        int dimension = r_process_info[DOMAIN_SIZE];

        typedef ElementsArrayType::iterator ElementIterator;
        // This vector collects the ids of the particles that have been dettached
        // so that their id can be removed from the mOriginInletSubmodelPartIndexes map
        std::vector<int> ids_to_remove;

        #pragma omp parallel
        {
        std::vector<int> ids_to_remove_partial;
        #pragma omp for
        for (int k = 0; k < (int)r_modelpart.GetCommunicator().LocalMesh().Elements().size(); k++) {
            ElementIterator elem_it = r_modelpart.GetCommunicator().LocalMesh().Elements().ptr_begin() + k;
            if (elem_it->IsNot(NEW_ENTITY)) continue;
            if (elem_it->Is(DEMFlags::BELONGS_TO_A_CLUSTER)) continue;

            SphericParticle& spheric_particle = dynamic_cast<SphericParticle&>(*elem_it);
            Node<3>& r_node = spheric_particle.GetGeometry()[0];

            bool have_just_stopped_touching = true;

            for (unsigned int i = 0; i < spheric_particle.mNeighbourElements.size(); i++) {
                SphericParticle* p_neighbour_particle = spheric_particle.mNeighbourElements[i];
                if(p_neighbour_particle == NULL) continue;

                Node<3>& neighbour_node = p_neighbour_particle->GetGeometry()[0];

                const double indentation = CalculateNormalizedIndentation(spheric_particle, *p_neighbour_particle);
                const bool indentation_is_significant_for_release = indentation > mNormalizedMaxIndentationForRelease*spheric_particle.GetInteractionRadius();
                const bool indentation_is_significant_for_injection = indentation > mNormalizedMaxIndentationForNewParticleCreation*spheric_particle.GetInteractionRadius();
                const bool i_am_injected_he_is_injector = r_node.IsNot(BLOCKED) && neighbour_node.Is(BLOCKED);
                const bool i_am_injector_he_is_injected = r_node.Is(BLOCKED) && neighbour_node.IsNot(BLOCKED);

                if (i_am_injected_he_is_injector && indentation_is_significant_for_release) {
                    have_just_stopped_touching = false;
                    UpdateInjectedParticleVelocity(spheric_particle, *p_neighbour_particle);
                    break;
                }

                if (i_am_injector_he_is_injected && indentation_is_significant_for_injection) {
                    have_just_stopped_touching = false;
                    break;
                }
            }

            if (have_just_stopped_touching) {
                if (r_node.IsNot(BLOCKED)) {//The ball must be freed
                    RemoveInjectionConditions(spheric_particle, dimension);
                    ids_to_remove_partial.push_back(spheric_particle.Id());
                    UpdateTotalThroughput(spheric_particle);
                }
                else {
                    //Inlet BLOCKED nodes are ACTIVE when injecting, so when they cease to be in contact with other balls, ACTIVE is set to 'false', as they become available for injecting new elements.
                    r_node.Set(ACTIVE, false);
                    elem_it->Set(ACTIVE, false);
                }
            }
        }

        // removing dettached particle ids from map
        #pragma omp critical
        {
            ids_to_remove.insert(ids_to_remove.end(), ids_to_remove_partial.begin(), ids_to_remove_partial.end());

            for (unsigned int i = 0; i < ids_to_remove.size(); ++i){
                mOriginInletSubmodelPartIndexes.erase(ids_to_remove[i]);
            }
        }
        }

    } //Dettach





    //************************ Typical process struct
    /// Execute method is used to execute the Injector algorithms.
    void ApplyParticleInjectionProcess::Execute() {}

    void ApplyParticleInjectionProcess::ExecuteInitialize() {}

    void ApplyParticleInjectionProcess::ExecuteInitializeSolutionStep()
    {
        KRATOS_TRY;

        // Check if it is injecting inside the defined interval.
        const double time = mrModelPart.GetProcessInfo()[TIME];
        if(!mInterval.IsInInterval(time)) return;

        KRATOS_CATCH("");
    }

    void ApplyParticleInjectionProcess::ExecuteFinalizeSolutionStep() {}


    std::string ApplyParticleInjectionProcess::Info() const
    {
        return "ApplyParticleInjectionProcess";
    }

    void ApplyParticleInjectionProcess::PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ApplyParticleInjectionProcess";
    }

    void ApplyParticleInjectionProcess::PrintData(std::ostream& rOStream) const {}

    /// output stream function
    inline std::ostream& operator << (
        std::ostream& rOStream,
        const ApplyParticleInjectionProcess& rThis)
    {
        rThis.PrintData(rOStream);
        return rOStream;
    }
}

