//
// Author: Miguel AngelCeligueta, maceli@cimne.upc.edu
//

#include "explicit_solver_continuum.h"
#include "utilities/parallel_utilities.h"
#include "utilities/atomic_utilities.h"

namespace Kratos {

    void ContinuumExplicitSolverStrategy::Initialize() {

        KRATOS_TRY

        ModelPart& r_model_part = GetModelPart();
        ModelPart& fem_model_part = GetFemModelPart();
        ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
        SendProcessInfoToClustersModelPart();

        if (r_model_part.GetCommunicator().MyPID() == 0) {
            KRATOS_INFO("DEM") << "------------------CONTINUUM SOLVER STRATEGY---------------------" << "\n" << std::endl;
        }

        mNumberOfThreads = OpenMPUtils::GetNumThreads();
        DisplayThreadInfo();

        RebuildListOfSphericParticles <SphericContinuumParticle> (r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericContinuumParticles);
        RebuildListOfSphericParticles <SphericContinuumParticle> (r_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostSphericContinuumParticles);
        RebuildListOfSphericParticles <SphericParticle> (r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericParticles);
        RebuildListOfSphericParticles <SphericParticle> (r_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostSphericParticles);

        mSearchControlVector.resize(mNumberOfThreads);
        for (int i = 0; i < mNumberOfThreads; i++) mSearchControlVector[i] = 0;

        PropertiesProxiesManager().CreatePropertiesProxies(r_model_part, *mpInlet_model_part, *mpCluster_model_part);

        RepairPointersToNormalProperties(mListOfSphericParticles);
        RepairPointersToNormalProperties(mListOfGhostSphericParticles);

        RebuildPropertiesProxyPointers(mListOfSphericParticles);
        RebuildPropertiesProxyPointers(mListOfGhostSphericParticles);

        GetSearchControl() = r_process_info[SEARCH_CONTROL];

        InitializeDEMElements();
        InitializeFEMElements();
        UpdateMaxIdOfCreatorDestructor();
        InitializeClusters(); // This adds elements to the balls modelpart

        RebuildListOfSphericParticles <SphericContinuumParticle> (r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericContinuumParticles);
        RebuildListOfSphericParticles <SphericParticle> (r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericParticles);
        RebuildListOfSphericParticles <SphericContinuumParticle> (r_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostSphericContinuumParticles);
        RebuildListOfSphericParticles <SphericParticle> (r_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostSphericParticles);

        InitializeSolutionStep();
        ApplyInitialConditions();

        // Search Neighbors with tolerance (after first repartition process)
        SetSearchRadiiOnAllParticles(r_model_part, r_process_info[SEARCH_RADIUS_INCREMENT], 1.0);
        SearchNeighbours();
        MeshRepairOperations();
        SearchNeighbours();

        const bool automatic_skin_computation = r_process_info[AUTOMATIC_SKIN_COMPUTATION];
        const double factor_radius = r_process_info[SKIN_FACTOR_RADIUS];

        if (automatic_skin_computation) {
            ResetSkinParticles(r_model_part);
            ComputeSkin(r_model_part, factor_radius);
        }

        if (GetDeltaOption() == 2) {
            SetCoordinationNumber(r_model_part);
            if (automatic_skin_computation) {
                ComputeSkin(r_model_part, factor_radius);
                SetCoordinationNumber(r_model_part);
            }
        }

        RebuildListOfSphericParticles <SphericContinuumParticle> (r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericContinuumParticles);
        RebuildListOfSphericParticles <SphericParticle> (r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericParticles);
        RebuildListOfSphericParticles <SphericContinuumParticle> (r_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostSphericContinuumParticles);
        RebuildListOfSphericParticles <SphericParticle> (r_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostSphericParticles);

        bool has_mpi = false;
        Check_MPI(has_mpi);

        if (has_mpi) {
            RepairPointersToNormalProperties(mListOfSphericParticles);
            RepairPointersToNormalProperties(mListOfGhostSphericParticles);
        }

        RebuildPropertiesProxyPointers(mListOfSphericParticles);
        RebuildPropertiesProxyPointers(mListOfGhostSphericParticles);

        if (has_mpi) {
            //RebuildListsOfPointersOfEachParticle(); //Serialized pointers are lost, so we rebuild them using Id's
        }

        // Set Initial Contacts
        if (r_process_info[CASE_OPTION] != 0) {
            SetInitialDemContacts();
        }

        if (r_process_info[CRITICAL_TIME_OPTION]) {
            //InitialTimeStepCalculation();   //obsolete call
            CalculateMaxTimeStep();
        }

        ComputeNewNeighboursHistoricalData();

        if (fem_model_part.Nodes().size() > 0) {
            SetSearchRadiiWithFemOnAllParticles(r_model_part, mpDem_model_part->GetProcessInfo()[SEARCH_RADIUS_INCREMENT_FOR_WALLS], 1.0);
            SearchRigidFaceNeighbours();
            SetInitialFemContacts();
            ComputeNewRigidFaceNeighboursHistoricalData();
        }

        if (mRemoveBallsInitiallyTouchingWallsOption) {
            MarkToDeleteAllSpheresInitiallyIndentedWithFEM(*mpDem_model_part);
            mpParticleCreatorDestructor->DestroyParticles(r_model_part);
            RebuildListOfSphericParticles<SphericParticle>(r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericParticles);
            RebuildListOfSphericParticles<SphericParticle>(r_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostSphericParticles);

            // Search Neighbours and related operations
            SetSearchRadiiOnAllParticles(*mpDem_model_part, mpDem_model_part->GetProcessInfo()[SEARCH_RADIUS_INCREMENT], 1.0);
            SearchNeighbours();
            ComputeNewNeighboursHistoricalData();

            SetSearchRadiiOnAllParticles(*mpDem_model_part, mpDem_model_part->GetProcessInfo()[SEARCH_RADIUS_INCREMENT_FOR_WALLS], 1.0);
            SearchRigidFaceNeighbours(); //initial search is performed with hierarchical method in any case MSI
            ComputeNewRigidFaceNeighboursHistoricalData();
        }

        AttachSpheresToStickyWalls();

        if (r_process_info[CONTACT_MESH_OPTION] == 1) {
            CreateContactElements();
            InitializeContactElements();
        }

        r_model_part.GetCommunicator().SynchronizeElementalNonHistoricalVariable(NEIGHBOUR_IDS);
        r_model_part.GetCommunicator().SynchronizeElementalNonHistoricalVariable(NEIGHBOURS_CONTACT_AREAS);

        if (r_process_info[CASE_OPTION] != 0) {
            CalculateMeanContactArea();
            CalculateMaxSearchDistance();
        }
        ComputeNodalArea();

        KRATOS_CATCH("")
    }// Initialize()

    double ContinuumExplicitSolverStrategy::SolveSolutionStep() {

        KRATOS_TRY

        ModelPart& r_model_part = GetModelPart();

        bool has_mpi = false;
        VariablesList r_modelpart_nodal_variables_list = r_model_part.GetNodalSolutionStepVariablesList();
        if (r_modelpart_nodal_variables_list.Has(PARTITION_INDEX)) has_mpi = true;

        SearchDEMOperations(r_model_part, has_mpi);
        SearchFEMOperations(r_model_part, has_mpi);
        ForceOperations(r_model_part);
        PerformTimeIntegrationOfMotion();

        KRATOS_CATCH("")

        return 0.0;

    }//SolveSolutionStep()

    void ContinuumExplicitSolverStrategy::SearchFEMOperations(ModelPart& r_model_part, bool has_mpi) {
        ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
        const int time_step = r_process_info[TIME_STEPS];
        const bool is_time_to_search_neighbours = (time_step + 1) % GetNStepSearch() == 0 && (time_step > 0); //Neighboring search. Every N times.

        if (is_time_to_search_neighbours) {
            SetSearchRadiiWithFemOnAllParticles(r_model_part, mpDem_model_part->GetProcessInfo()[SEARCH_RADIUS_INCREMENT_FOR_WALLS], 1.0);
            SearchRigidFaceNeighbours();
            ComputeNewRigidFaceNeighboursHistoricalData();
        }
    }

    void ContinuumExplicitSolverStrategy::SearchDEMOperations(ModelPart& r_model_part, bool has_mpi) {

        ProcessInfo& r_process_info = r_model_part.GetProcessInfo();

        if (r_process_info[SEARCH_CONTROL] == 0) {

            ElementsArrayType& rElements = r_model_part.GetCommunicator().LocalMesh().Elements();
            int some_bond_is_broken = 0;

            block_for_each(rElements, [&](ModelPart::ElementType& rElement) {

                SphericContinuumParticle& r_sphere = dynamic_cast<SphericContinuumParticle&>(rElement);

                for (int j=0; j<(int) r_sphere.mContinuumInitialNeighborsSize; j++) {
                    if (r_sphere.mIniNeighbourFailureId[j] != 0) {
                        AtomicAdd(some_bond_is_broken, 1);
                        break;
                    }
                }

            });

            if (some_bond_is_broken > 0) {
                r_process_info[SEARCH_CONTROL] = 1;
                KRATOS_WARNING("DEM") << "From now on, the search is activated because some failure occurred " << std::endl;
            }

        }

        const int time_step = r_process_info[TIME_STEPS];
        const double time = r_process_info[TIME];
        const bool is_time_to_search_neighbours = (time_step + 1) % GetNStepSearch() == 0 && (time_step > 0); //Neighboring search. Every N times.

        if (r_process_info[SEARCH_CONTROL] > 0) {

            if (is_time_to_search_neighbours) {

                if (r_process_info[BOUNDING_BOX_OPTION] && time >= r_process_info[BOUNDING_BOX_START_TIME] && time <= r_process_info[BOUNDING_BOX_STOP_TIME]) {
                    BoundingBoxUtility();
                } else {
                    GetParticleCreatorDestructor()->DestroyParticles(r_model_part);
                    GetParticleCreatorDestructor()->DestroyContactElements(*mpContact_model_part);
                }

                RebuildListOfSphericParticles <SphericContinuumParticle> (r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericContinuumParticles); //These lists are necessary for the loop in SearchNeighbours
                RebuildListOfSphericParticles <SphericParticle> (r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericParticles);

                SetSearchRadiiOnAllParticles(r_model_part, r_process_info[SEARCH_RADIUS_INCREMENT], r_process_info[CONTINUUM_SEARCH_RADIUS_AMPLIFICATION_FACTOR]);

                SearchNeighbours(); //the amplification factor has been modified after the first search.

                RebuildListOfSphericParticles <SphericContinuumParticle> (r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericContinuumParticles); //These lists are necessary because the elements in this partition might have changed.
                RebuildListOfSphericParticles <SphericParticle> (r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericParticles);
                RebuildListOfSphericParticles <SphericContinuumParticle> (r_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostSphericContinuumParticles);
                RebuildListOfSphericParticles <SphericParticle> (r_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostSphericParticles);

                if (has_mpi) {
                    RepairPointersToNormalProperties(mListOfSphericParticles);
                    RepairPointersToNormalProperties(mListOfGhostSphericParticles);
                }

                RebuildPropertiesProxyPointers(mListOfSphericParticles);
                RebuildPropertiesProxyPointers(mListOfGhostSphericParticles);

                ComputeNewNeighboursHistoricalData();

                MarkNewSkinParticles();

                r_process_info[SEARCH_CONTROL] = 2;
            } else {
                r_process_info[SEARCH_CONTROL] = 1;
            }

            if (r_process_info[CONTACT_MESH_OPTION]) {
                CreateContactElements();
                InitializeContactElements();
            }
        }
        //Synch this var.
        r_process_info[SEARCH_CONTROL] = r_model_part.GetCommunicator().GetDataCommunicator().MaxAll(r_process_info[SEARCH_CONTROL]);
    }

    void ContinuumExplicitSolverStrategy::MarkNewSkinParticles() {

        KRATOS_TRY

        #pragma omp parallel
        {
            const int number_of_particles = (int)mListOfSphericContinuumParticles.size();

            #pragma omp for
            for (int i = 0; i < number_of_particles; i++) {
                mListOfSphericContinuumParticles[i]->MarkNewSkinParticlesDueToBreakage();
            }
        }

        KRATOS_CATCH("")
    }

    void ContinuumExplicitSolverStrategy::ResetSkinParticles(ModelPart& r_model_part) {
        auto& pNodes = r_model_part.GetCommunicator().LocalMesh().Nodes();
        #pragma omp parallel for
        for (int k = 0; k < (int)pNodes.size(); k++) {
            auto it = pNodes.begin() + k;
            it->FastGetSolutionStepValue(SKIN_SPHERE) = 0.0;
        }
    }

    void ContinuumExplicitSolverStrategy::ComputeSkin(ModelPart& rSpheresModelPart, const double factor_radius) {

        ElementsArrayType& pElements = rSpheresModelPart.GetCommunicator().LocalMesh().Elements();
        ProcessInfo& r_process_info = rSpheresModelPart.GetProcessInfo();
        const unsigned int problem_dimension = r_process_info[DOMAIN_SIZE];
        unsigned int minimum_bonds_to_check_if_a_particle_is_skin;

        if (problem_dimension == 2) {
            minimum_bonds_to_check_if_a_particle_is_skin = 4;
        } else return; // TODO: IMPLEMENT THIS ALSO IN 3D

        #pragma omp parallel for
        for (int k = 0; k < (int)pElements.size(); k++) {

            ElementsArrayType::iterator it = pElements.ptr_begin() + k;
            Element* p_element = &(*it);
            SphericContinuumParticle* p_sphere = dynamic_cast<SphericContinuumParticle*>(p_element);

            const array_1d<double, 3> element_center = p_sphere->GetGeometry()[0].Coordinates();
            const double element_radius              = p_sphere->GetGeometry()[0].FastGetSolutionStepValue(RADIUS);
            array_1d<double, 3> vector_from_element_to_neighbour = ZeroVector(3);
            array_1d<double, 3> sum_of_vectors_from_element_to_neighbours = ZeroVector(3);
            unsigned const int number_of_neighbors = p_sphere->mNeighbourElements.size();
            double modulus_of_vector_from_element_to_neighbour = 0.0;

            for (unsigned int i = 0; i < number_of_neighbors; i++) {
                SphericContinuumParticle* neighbour_iterator = dynamic_cast<SphericContinuumParticle*>(p_sphere->mNeighbourElements[i]);
                if (!neighbour_iterator) continue;
                const array_1d<double, 3> neigbour_center = neighbour_iterator->GetGeometry()[0].Coordinates();
                vector_from_element_to_neighbour = neigbour_center - element_center;
                modulus_of_vector_from_element_to_neighbour = DEM_MODULUS_3(vector_from_element_to_neighbour);
                DEM_MULTIPLY_BY_SCALAR_3(vector_from_element_to_neighbour, element_radius / modulus_of_vector_from_element_to_neighbour);
                sum_of_vectors_from_element_to_neighbours += vector_from_element_to_neighbour;
            }
            const double sum_modulus = DEM_MODULUS_3(sum_of_vectors_from_element_to_neighbours);

            if ((number_of_neighbors < minimum_bonds_to_check_if_a_particle_is_skin) || (sum_modulus > factor_radius * element_radius)) {
                p_sphere->GetGeometry()[0].FastGetSolutionStepValue(SKIN_SPHERE) = 1.0;
            }
        }
    }

    void ContinuumExplicitSolverStrategy::ComputeNewNeighboursHistoricalData() {
        KRATOS_TRY
        ProcessInfo& r_process_info = GetModelPart().GetProcessInfo();

        #pragma omp parallel
        {
            DenseVector<int> temp_neighbours_ids; //We are passing all these temporal vectors as arguments because creating them inside the function is slower (memory allocation and deallocation)
            std::vector<array_1d<double, 3> > temp_neighbour_elastic_contact_forces;
            std::vector<SphericParticle*> temp_neighbour_elements;

            const int number_of_particles = (int) mListOfSphericContinuumParticles.size();

            #pragma omp for
            for (int i = 0; i < number_of_particles; i++) {
                mListOfSphericContinuumParticles[i]->ReorderAndRecoverInitialPositionsAndFilter(temp_neighbour_elements);
                mListOfSphericContinuumParticles[i]->UpdateContinuumNeighboursVector(r_process_info);
                mListOfSphericContinuumParticles[i]->ComputeNewNeighboursHistoricalData(temp_neighbours_ids, temp_neighbour_elastic_contact_forces);
            }
        }

        KRATOS_CATCH("")
    }

    void ContinuumExplicitSolverStrategy::CreateContactElements() {
        KRATOS_TRY

        std::string ElementName;
        ElementName = std::string("ParticleContactElement");
        const Element& rReferenceElement = KratosComponents<Element>::Get(ElementName);

        //Here we are going to create contact elements when we are on a target particle and we see a neighbor whose id is higher than ours.
        //We create also a pointer from the node to the element, after creating it.
        //When our particle has a higher ID than the neighbor we also create a pointer to the (previously) created contact element.
        //We proceed in this way because we want to have the pointers to contact elements in a list in the same order as the initial elements order.

        const int number_of_particles = (int) mListOfSphericContinuumParticles.size();
        int used_bonds_counter = 0;
        #pragma omp parallel
        {
            #pragma omp for
            for (int i = 0; i < number_of_particles; i++) {
                unsigned int continuous_initial_neighbors_size = mListOfSphericContinuumParticles[i]->mContinuumInitialNeighborsSize;
                mListOfSphericContinuumParticles[i]->mBondElements.resize(continuous_initial_neighbors_size);
                for (unsigned int j = 0; j < mListOfSphericContinuumParticles[i]->mBondElements.size(); j++) {
                    mListOfSphericContinuumParticles[i]->mBondElements[j] = NULL;
                }
            }

            int private_counter = 0;
            Element::Pointer p_new_contact_element;
            #pragma omp for
            for (int i = 0; i < number_of_particles; i++) {
                bool add_new_bond = true;
                std::vector<SphericParticle*>& neighbour_elements = mListOfSphericContinuumParticles[i]->mNeighbourElements;
                unsigned int continuous_initial_neighbors_size = mListOfSphericContinuumParticles[i]->mContinuumInitialNeighborsSize;

                for (unsigned int j = 0; j < continuous_initial_neighbors_size; j++) {
                    SphericContinuumParticle* neighbour_element = dynamic_cast<SphericContinuumParticle*> (neighbour_elements[j]);
                    if (neighbour_element == NULL) continue; //The initial neighbor was deleted at some point in time!!
                    if (mListOfSphericContinuumParticles[i]->Id() > neighbour_element->Id()) continue;

                    #pragma omp critical
                    {
                        if (used_bonds_counter < (int) (*mpContact_model_part).Elements().size()) {
                            add_new_bond = false;
                            private_counter = used_bonds_counter;
                            used_bonds_counter++;
                        }
                    }
                    if (!add_new_bond) {
                        Element::Pointer& p_old_contact_element = (*mpContact_model_part).Elements().GetContainer()[private_counter];
                        p_old_contact_element->GetGeometry()(0) = mListOfSphericContinuumParticles[i]->GetGeometry()(0);
                        p_old_contact_element->GetGeometry()(1) = neighbour_element->GetGeometry()(0);
                        p_old_contact_element->SetId(used_bonds_counter);
                        p_old_contact_element->SetProperties(mListOfSphericContinuumParticles[i]->pGetProperties());
                        ParticleContactElement* p_bond = dynamic_cast<ParticleContactElement*> (p_old_contact_element.get());
                        mListOfSphericContinuumParticles[i]->mBondElements[j] = p_bond;
                    } else {
                        Geometry<Node<3> >::PointsArrayType NodeArray(2);
                        NodeArray.GetContainer()[0] = mListOfSphericContinuumParticles[i]->GetGeometry()(0);
                        NodeArray.GetContainer()[1] = neighbour_element->GetGeometry()(0);
                        const Properties::Pointer& properties = mListOfSphericContinuumParticles[i]->pGetProperties();
                        p_new_contact_element = rReferenceElement.Create(used_bonds_counter + 1, NodeArray, properties);

                        #pragma omp critical
                        {
                            (*mpContact_model_part).Elements().push_back(p_new_contact_element);
                            used_bonds_counter++;
                        }
                        ParticleContactElement* p_bond = dynamic_cast<ParticleContactElement*> (p_new_contact_element.get());
                        mListOfSphericContinuumParticles[i]->mBondElements[j] = p_bond;
                    }

                }
            }

            #pragma omp single
            {
                if ((int) (*mpContact_model_part).Elements().size() > used_bonds_counter) {
                    (*mpContact_model_part).Elements().erase((*mpContact_model_part).Elements().ptr_begin() + used_bonds_counter, (*mpContact_model_part).Elements().ptr_end());
                }
            }

            #pragma omp for
            for (int i = 0; i < number_of_particles; i++) {
                std::vector<SphericParticle*>& neighbour_elements = mListOfSphericContinuumParticles[i]->mNeighbourElements;
                unsigned int continuous_initial_neighbors_size = mListOfSphericContinuumParticles[i]->mContinuumInitialNeighborsSize;

                for (unsigned int j = 0; j < continuous_initial_neighbors_size; j++) {
                    SphericContinuumParticle* neighbour_element = dynamic_cast<SphericContinuumParticle*> (neighbour_elements[j]);
                    if (neighbour_element == NULL) continue; //The initial neighbor was deleted at some point in time!!
                    //ATTENTION: Ghost nodes do not have mContinuumIniNeighbourElements in general, so this bond will remain as NULL!!
                    if (mListOfSphericContinuumParticles[i]->Id() < neighbour_element->Id()) continue;
                    //In all functions using mBondElements we must check that this bond is not used.

                    for (unsigned int k = 0; k < neighbour_element->mContinuumInitialNeighborsSize; k++) {
                        //ATTENTION: Ghost nodes do not have mContinuumIniNeighbourElements in general, so this bond will remain as NULL!!
                        //In all functions using mBondElements we must check that this bond is not used.
                        if (neighbour_element->mNeighbourElements[k] == NULL) continue; //The initial neighbor was deleted at some point in time!!
                        if (neighbour_element->mNeighbourElements[k]->Id() == mListOfSphericContinuumParticles[i]->Id()) {
                            ParticleContactElement* bond = neighbour_element->mBondElements[k];
                            mListOfSphericContinuumParticles[i]->mBondElements[j] = bond;
                            break;
                        }
                    }
                }
            }

            //Renumbering the Id's of the bonds to make them unique and consecutive (otherwise the Id's are repeated)
            #pragma omp for
            for(int i=0; i<(int)(*mpContact_model_part).Elements().size(); i++) {
                (*mpContact_model_part).Elements().GetContainer()[i]->SetId(i+1);
            }

        } //#pragma omp parallel
        KRATOS_CATCH("")
    } //CreateContactElements

    void ContinuumExplicitSolverStrategy::SetCoordinationNumber(ModelPart& r_model_part) {

        ProcessInfo& r_process_info = r_model_part.GetProcessInfo();

        const double desired_coordination_number = r_process_info[COORDINATION_NUMBER];
        double standard_dev = 0;
        double current_coordination_number = ComputeCoordinationNumber(standard_dev);
        int iteration = 0;
        int maxiteration = 300;
        double& added_search_distance = r_process_info[SEARCH_RADIUS_INCREMENT];
        const bool local_coordination_option = r_process_info[LOCAL_COORDINATION_NUMBER_OPTION];
        const bool global_coordination_option = r_process_info[GLOBAL_COORDINATION_NUMBER_OPTION];
        added_search_distance = 0.0;

        KRATOS_INFO("DEM") << "Setting up Coordination Number (input = "<<desired_coordination_number<<") by increasing or decreasing the search radius. ";
        KRATOS_INFO_IF("", local_coordination_option) << "Local extension activated. ";
        KRATOS_INFO_IF("", global_coordination_option) << "Global extension activated. ";
        KRATOS_INFO("") << std::endl;

        KRATOS_ERROR_IF(desired_coordination_number <= 0.0) << "The specified Coordination Number is less or equal to zero, N.C. = " << desired_coordination_number << std::endl;

        //STAGE 1, Local Coordination Number
        double tolerance = 0.5;
        double max_factor_between_iterations = 1.02;
        double relative_error = 1000;
        if(local_coordination_option) {
            KRATOS_INFO("DEM")<<"Now iterating for local coordination number..."<<std::endl;
            while (std::abs(relative_error) > tolerance) {
                if (iteration >= maxiteration) break;
                iteration++;
                if (current_coordination_number == 0.0) {
                    KRATOS_WARNING("DEM") << "Coordination Number method not supported in this case" << "\n" << std::endl;
                    KRATOS_ERROR << "The specified tangency method is not supported for this problem, please use absolute value instead" << std::endl;
                    break;
                }

                std::vector<double> total_error;
                mNumberOfThreads = OpenMPUtils::GetNumThreads();
                total_error.resize(mNumberOfThreads);

                #pragma omp parallel for
                for (int i = 0; i < static_cast<int>(mListOfSphericContinuumParticles.size()); ++i) {
                    const std::size_t neighbour_elements_size = mListOfSphericContinuumParticles[i]->mNeighbourElements.size();
                    const double old_amplification = mListOfSphericContinuumParticles[i]->mLocalRadiusAmplificationFactor;
                    double adapted_to_skin_or_not_skin_desired_cn = desired_coordination_number;
                    if(mListOfSphericContinuumParticles[i]->IsSkin()) {
                        adapted_to_skin_or_not_skin_desired_cn = 0.6*desired_coordination_number;
                    }

                    if(neighbour_elements_size) {
                        if (neighbour_elements_size != std::round(adapted_to_skin_or_not_skin_desired_cn)) {
                            mListOfSphericContinuumParticles[i]->mLocalRadiusAmplificationFactor *= std::sqrt(adapted_to_skin_or_not_skin_desired_cn / (double)(neighbour_elements_size));
                        }
                    }
                    else {
                        mListOfSphericContinuumParticles[i]->mLocalRadiusAmplificationFactor *= max_factor_between_iterations;
                    }
                    if(mListOfSphericContinuumParticles[i]->mLocalRadiusAmplificationFactor > max_factor_between_iterations * old_amplification) mListOfSphericContinuumParticles[i]->mLocalRadiusAmplificationFactor = max_factor_between_iterations* old_amplification;
                    if(mListOfSphericContinuumParticles[i]->mLocalRadiusAmplificationFactor < old_amplification / max_factor_between_iterations) mListOfSphericContinuumParticles[i]->mLocalRadiusAmplificationFactor = old_amplification / max_factor_between_iterations;

                    const int error_of_this_particle = std::abs(std::round(adapted_to_skin_or_not_skin_desired_cn) - (neighbour_elements_size));
                    total_error[OpenMPUtils::ThisThread()]  += (double)(error_of_this_particle);
                }
                double total_absolute_error = 0.0;
                for (size_t i=0; i<total_error.size();i++) total_absolute_error += total_error[i];
                relative_error = total_absolute_error/mListOfSphericContinuumParticles.size();

                SetSearchRadiiOnAllParticles(r_model_part, added_search_distance, 1.0);
                SearchNeighbours();
            }//while
            current_coordination_number = ComputeCoordinationNumber(standard_dev);
            KRATOS_INFO("DEM")<<"Coordination number reached after local operations = "<<current_coordination_number<<std::endl;
        }

        if(global_coordination_option) {
            KRATOS_INFO("DEM")<<"Now iterating for global coordination number..."<<std::endl;
            //STAGE 2, Global Coordination Number
            tolerance = 1e-4;
            max_factor_between_iterations = 1.1;
            double amplification = 1.0;
            double old_old_amplification = 1.0;
            double old_amplification = 1.0;

            while (std::abs(current_coordination_number / desired_coordination_number - 1.0) > tolerance) {
                if (iteration >= maxiteration) break;
                iteration++;
                if (current_coordination_number == 0.0) {
                    KRATOS_WARNING("DEM") << "Coordination Number method not supported in this case" << "\n" << std::endl;
                    KRATOS_ERROR << "The specified tangency method is not supported for this problem, please use absolute value instead" << std::endl;
                    break;
                }
                old_old_amplification = old_amplification;
                old_amplification = amplification;
                amplification *= std::pow(desired_coordination_number / current_coordination_number, 1.0/3.0);

                if(amplification > max_factor_between_iterations * old_amplification) amplification = max_factor_between_iterations* old_amplification;
                if(amplification < old_amplification / max_factor_between_iterations) amplification = old_amplification / max_factor_between_iterations;
                if(std::abs(amplification - old_amplification) >= std::abs(old_amplification - old_old_amplification) - std::numeric_limits<double>::epsilon()) {
                    amplification = 0.5 * (amplification + old_amplification);
                }
                if ( amplification < 1.0 && !local_coordination_option) {
                    iteration = maxiteration;
                    break;
                }
                SetSearchRadiiOnAllParticles(r_model_part, added_search_distance, amplification);
                SearchNeighbours();
                current_coordination_number = ComputeCoordinationNumber(standard_dev);
            }//while

                if (iteration < maxiteration){
                KRATOS_INFO("DEM") << "The iterative procedure converged after " << iteration << " iterations, to value \e[1m" << current_coordination_number << "\e[0m using a global amplification of radius of " << amplification << ". " << "\n" << std::endl;
                KRATOS_INFO("DEM") << "Standard deviation for achieved coordination number is " << standard_dev << ". " << "\n" << std::endl;
                //KRATOS_INFO("DEM") << "This means that most particles (about 68% of the total particles, assuming a normal distribution) have a coordination number within " <<  standard_dev << " contacts of the mean (" << current_coordination_number-standard_dev << "–" << current_coordination_number+standard_dev << " contacts). " << "\n" << std::endl;
                r_process_info[CONTINUUM_SEARCH_RADIUS_AMPLIFICATION_FACTOR] = amplification;
            }

            else {
                KRATOS_WARNING("DEM") << "Coordination Number iterative procedure did NOT converge after " << iteration << " iterations. Coordination number reached is " << current_coordination_number << ". Desired number was " <<desired_coordination_number << "\n" << std::endl;
                KRATOS_ERROR << "Please use a Absolute tolerance instead " << std::endl;
                //NOTE: if it doesn't converge, problems occur with contact mesh and rigid face contact.
            }
        }

        KRATOS_INFO("DEM") <<std::endl;
    } //SetCoordinationNumber

    double ContinuumExplicitSolverStrategy::ComputeCoordinationNumber(double& standard_dev) {
        KRATOS_TRY

        ModelPart& r_model_part = GetModelPart();

        int total_contacts = 0;
        const int number_of_particles = (int) mListOfSphericParticles.size();
        std::vector<int> neighbour_counter;
        std::vector<int> sum;
        std::vector<int> number_of_non_skin_particles;
        double total_sum = 0.0;
        int total_non_skin_particles = 0;

        mNumberOfThreads = OpenMPUtils::GetNumThreads();
        neighbour_counter.resize(mNumberOfThreads);
        sum.resize(mNumberOfThreads);
        number_of_non_skin_particles.resize(mNumberOfThreads);

        for (int i = 0; i < mNumberOfThreads; i++) {
            total_contacts = 0;
            total_sum = 0;
            neighbour_counter[i] = 0;
            sum[i] = 0;
            number_of_non_skin_particles[i] = 0;
        }

        #pragma omp parallel for
        for (int i = 0; i < number_of_particles; i++) {
            SphericParticle* element = mListOfSphericParticles[i];
            SphericContinuumParticle* continuous_element = dynamic_cast<SphericContinuumParticle*> (element);
            if (continuous_element->IsSkin()) {
                continue;
            }
            neighbour_counter[OpenMPUtils::ThisThread()] += mListOfSphericParticles[i]->mNeighbourElements.size();
            sum[OpenMPUtils::ThisThread()] += (mListOfSphericParticles[i]->mNeighbourElements.size() - 10.0 )*(mListOfSphericParticles[i]->mNeighbourElements.size() - 10.0 );
            number_of_non_skin_particles[OpenMPUtils::ThisThread()] += 1;
        }

        for (int i = 0; i < mNumberOfThreads; i++) {
            total_contacts += neighbour_counter[i];
            total_sum += sum[i];
            total_non_skin_particles += number_of_non_skin_particles[i];
        }

        int global_total_contacts = r_model_part.GetCommunicator().GetDataCommunicator().SumAll(total_contacts);
        int global_total_non_skin_particles = r_model_part.GetCommunicator().GetDataCommunicator().SumAll(total_non_skin_particles);

        double coord_number = double(global_total_contacts) / double(global_total_non_skin_particles);

        standard_dev = sqrt(total_sum / global_total_non_skin_particles);

        return coord_number;

        KRATOS_CATCH("")
    }

    void ContinuumExplicitSolverStrategy::SetSearchRadiiOnAllParticles(ModelPart& r_model_part, const double added_search_distance, const double amplification) {
        KRATOS_TRY
        const int number_of_elements = r_model_part.GetCommunicator().LocalMesh().NumberOfElements();
        #pragma omp parallel for
        for (int i = 0; i < number_of_elements; i++) {
            mListOfSphericContinuumParticles[i]->SetSearchRadius(amplification * mListOfSphericContinuumParticles[i]->mLocalRadiusAmplificationFactor * (added_search_distance + mListOfSphericContinuumParticles[i]->GetRadius()));
        }
        KRATOS_CATCH("")
    }

    void ContinuumExplicitSolverStrategy::BoundingBoxUtility(bool is_time_to_mark_and_remove) {
        KRATOS_TRY

        ModelPart& r_model_part = GetModelPart();
        ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
        ParticleCreatorDestructor::Pointer& p_creator_destructor = GetParticleCreatorDestructor();

        p_creator_destructor->MarkDistantParticlesForErasing(r_model_part);

        if (r_process_info[IS_TIME_TO_PRINT] && r_process_info[CONTACT_MESH_OPTION] == 1) {
            p_creator_destructor->MarkContactElementsForErasing(r_model_part, *mpContact_model_part);
            p_creator_destructor->DestroyContactElements(*mpContact_model_part);
        }
        p_creator_destructor->DestroyParticles(r_model_part);

        KRATOS_CATCH("")
    }

    void ContinuumExplicitSolverStrategy::Check_MPI(bool& has_mpi) {
        VariablesList r_modelpart_nodal_variables_list = GetModelPart().GetNodalSolutionStepVariablesList();
        if (r_modelpart_nodal_variables_list.Has(PARTITION_INDEX)) has_mpi = true;
    }

    void ContinuumExplicitSolverStrategy::CalculateMaxSearchDistance() {

        KRATOS_TRY

        ModelPart& r_model_part = GetModelPart();
        ProcessInfo& r_process_info = r_model_part.GetProcessInfo();

        bool has_mpi = false;
        Check_MPI(has_mpi);

        std::vector<double> thread_maxima(OpenMPUtils::GetNumThreads(), 0.0);
        const int number_of_particles = (int) mListOfSphericContinuumParticles.size();

        #pragma omp parallel for
        for (int i = 0; i < number_of_particles; i++) {
            double max_search_distance_for_this_sphere = mListOfSphericContinuumParticles[i]->CalculateMaxSearchDistance(has_mpi, r_process_info);
            const double ratio_search_distance_versus_radius = max_search_distance_for_this_sphere / mListOfSphericContinuumParticles[i]->GetRadius();
            if (ratio_search_distance_versus_radius > thread_maxima[OpenMPUtils::ThisThread()]) thread_maxima[OpenMPUtils::ThisThread()] = ratio_search_distance_versus_radius;
        }

        double maximum_across_threads = 0.0;
        for (int i = 0; i < OpenMPUtils::GetNumThreads(); i++) {
            if (thread_maxima[i] > maximum_across_threads) maximum_across_threads = thread_maxima[i];
        }

        double& ratio = r_process_info[CONTINUUM_SEARCH_RADIUS_AMPLIFICATION_FACTOR];

        if (maximum_across_threads > ratio) {
            ratio = maximum_across_threads;
        }

        const double max_ratio = r_process_info[MAX_AMPLIFICATION_RATIO_OF_THE_SEARCH_RADIUS];

        static unsigned int counter = 0;
        unsigned int maximum_number_of_prints = 5;

        if ((ratio > max_ratio) && (counter <= maximum_number_of_prints)) {
            KRATOS_INFO("DEM") <<std::endl;
            KRATOS_WARNING("DEM") <<"************************************************************************"<<std::endl;
            KRATOS_WARNING("DEM") <<"WARNING! The automatic extension of the search radius, based on mechanical"<<std::endl;
            KRATOS_WARNING("DEM") <<"reasons, is trying to extend more than "<<max_ratio<<" times the "<<std::endl;
            KRATOS_WARNING("DEM") <<"original particle radius!"<<std::endl;
            KRATOS_WARNING("DEM") <<"Some bonds might break for search reasons instead of mechanical reasons from now on"<<std::endl;
            KRATOS_WARNING("DEM") <<"because the ratio is limited to that value ("<<max_ratio<<" times as introduced by "<<std::endl;
            KRATOS_WARNING("DEM") <<"the input variable 'MaxAmplificationRatioOfSearchRadius')"<<std::endl;
            KRATOS_WARNING("DEM") <<"************************************************************************"<<std::endl;
            ratio = max_ratio;
        }

        ++counter;

        KRATOS_CATCH("")
    }


    void ContinuumExplicitSolverStrategy::MeshRepairOperations() {

        KRATOS_TRY

        const int number_of_particles = (int) mListOfSphericContinuumParticles.size();
        int particle_counter = 0.0;

        #pragma omp parallel for
        for (int i = 0; i < number_of_particles; i++) {
            bool result = mListOfSphericContinuumParticles[i]->OverlappedParticleRemoval();

            if (result == true) {particle_counter += 1;}
        }

        GetModelPart().GetCommunicator().SynchronizeElementalFlags();

        DestroyMarkedParticlesRebuildLists();

        //KRATOS_WARNING("DEM") << "Mesh repair complete. In MPI node " <<GetModelPart().GetCommunicator().MyPID()<<". "<< particle_counter << " particles were removed. " << "\n" << std::endl;
        int total_spheres_removed = GetModelPart().GetCommunicator().GetDataCommunicator().SumAll(particle_counter);

        if(GetModelPart().GetCommunicator().MyPID() == 0 && total_spheres_removed) {
            KRATOS_WARNING("DEM") << "A total of "<<total_spheres_removed<<" spheres were removed due to excessive overlapping." << std::endl;
        }

        KRATOS_CATCH("")
    }


    void ContinuumExplicitSolverStrategy::DestroyMarkedParticlesRebuildLists() {

        KRATOS_TRY

        ModelPart& r_model_part = GetModelPart();

        GetParticleCreatorDestructor()->DestroyParticles(r_model_part);

        RebuildListOfSphericParticles <SphericContinuumParticle> (r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericContinuumParticles); //These lists are necessary because the elements in this partition might have changed.
        RebuildListOfSphericParticles <SphericParticle> (r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericParticles);
        RebuildListOfSphericParticles <SphericContinuumParticle> (r_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostSphericContinuumParticles);
        RebuildListOfSphericParticles <SphericParticle> (r_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostSphericParticles);

        KRATOS_CATCH("")
    }

    void ContinuumExplicitSolverStrategy::CalculateMeanContactArea() {

        KRATOS_TRY

        ModelPart& r_model_part = GetModelPart();
        ProcessInfo& r_process_info = r_model_part.GetProcessInfo();

        bool has_mpi = false;
        Check_MPI(has_mpi);

        const int number_of_particles = (int) mListOfSphericContinuumParticles.size();

        #pragma omp parallel for
        for (int i = 0; i < number_of_particles; i++) { //Do not do this for the ghost particles!
            mListOfSphericContinuumParticles[i]->CalculateMeanContactArea(has_mpi, r_process_info);
        }

        KRATOS_CATCH("")
    }

    void ContinuumExplicitSolverStrategy::BreakAlmostBrokenSpheres() {

        KRATOS_TRY

        const int maximum_allowed_number_of_intact_bonds = (int) GetModelPart().GetProcessInfo()[MAX_NUMBER_OF_INTACT_BONDS_TO_CONSIDER_A_SPHERE_BROKEN];

        #pragma omp parallel for
        for (int i = 0; i < (int) mListOfSphericContinuumParticles.size(); i++) {

            int number_of_intact_bonds = 0;

            for (int j = 0; j < (int) mListOfSphericContinuumParticles[i]->mContinuumInitialNeighborsSize; j++) {

                if (!mListOfSphericContinuumParticles[i]->mIniNeighbourFailureId[j]) number_of_intact_bonds++;
                if (number_of_intact_bonds > maximum_allowed_number_of_intact_bonds) break;
            }

            if (number_of_intact_bonds <= maximum_allowed_number_of_intact_bonds) {

                for (int j = 0; j < (int) mListOfSphericContinuumParticles[i]->mContinuumInitialNeighborsSize; j++) {

                    if (!mListOfSphericContinuumParticles[i]->mIniNeighbourFailureId[j]) mListOfSphericContinuumParticles[i]->mIniNeighbourFailureId[j] = 8;
                }
            }
        }

        KRATOS_CATCH("")
    }

    void ContinuumExplicitSolverStrategy::SetInitialDemContacts() {
        KRATOS_TRY

        ProcessInfo& r_process_info = GetModelPart().GetProcessInfo();
        const int number_of_particles = (int) mListOfSphericContinuumParticles.size();

        #pragma omp parallel
        {
            #pragma omp for
            for (int i = 0; i < number_of_particles; i++) {
                mListOfSphericContinuumParticles[i]->SetInitialSphereContacts(r_process_info);
                mListOfSphericContinuumParticles[i]->CreateContinuumConstitutiveLaws();
            }

            #pragma omp for
            for (int i = 0; i < number_of_particles; i++) {
                mListOfSphericContinuumParticles[i]->ContactAreaWeighting();
            }
        }

        KRATOS_CATCH("")
    }

    void ContinuumExplicitSolverStrategy::SetInitialFemContacts() {
        KRATOS_TRY

        const int number_of_particles = (int) mListOfSphericContinuumParticles.size();

        #pragma omp parallel for
        for (int i = 0; i < number_of_particles; i++) {
            mListOfSphericContinuumParticles[i]->SetInitialFemContacts();
        }

        KRATOS_CATCH("")
    }

    void ContinuumExplicitSolverStrategy::FinalizeSolutionStep() {
        BaseType::FinalizeSolutionStep();
        FinalizeSolutionStepFEM();

        ProcessInfo& r_process_info = GetModelPart().GetProcessInfo();
        if (r_process_info[COMPUTE_STRESS_TENSOR_OPTION]) {
            const int number_of_particles = (int) mListOfSphericContinuumParticles.size();
            #pragma omp parallel
            {
                #pragma omp for
                for (int i = 0; i < number_of_particles; i++) {
                    mListOfSphericContinuumParticles[i]->GetStressTensorFromNeighbourStep1();
                }
                #pragma omp for
                for (int i = 0; i < number_of_particles; i++) {
                    mListOfSphericContinuumParticles[i]->GetStressTensorFromNeighbourStep2();
                }
                #pragma omp for
                for (int i = 0; i < number_of_particles; i++) {
                    mListOfSphericContinuumParticles[i]->GetStressTensorFromNeighbourStep3();
                }
            }
        }

        BreakAlmostBrokenSpheres();
    }

    void ContinuumExplicitSolverStrategy::FinalizeSolutionStepFEM() {
        KRATOS_TRY

        ConditionsArrayType& pConditions = GetFemModelPart().GetCommunicator().LocalMesh().Conditions();
        ProcessInfo& r_process_info = GetFemModelPart().GetProcessInfo();
        Vector rhs_cond;
        std::vector<unsigned int> condition_partition;
        OpenMPUtils::CreatePartition(mNumberOfThreads, pConditions.size(), condition_partition);

        #pragma omp parallel for private (rhs_cond)
        for (int k = 0; k < mNumberOfThreads; k++) {
            ConditionsArrayType::iterator it_begin = pConditions.ptr_begin() + condition_partition[k];
            ConditionsArrayType::iterator it_end = pConditions.ptr_begin() + condition_partition[k + 1];

            for (ConditionsArrayType::iterator it = it_begin; it != it_end; ++it) {
                it->FinalizeSolutionStep(r_process_info);
            }
        }

        KRATOS_CATCH("")
    }

}
