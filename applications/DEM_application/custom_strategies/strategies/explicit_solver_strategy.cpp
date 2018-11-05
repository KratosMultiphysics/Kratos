//
// Author: Miguel AngelCeligueta, maceli@cimne.upc.edu
//

#include "explicit_solver_strategy.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "geometries/point_3d.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

namespace Kratos {

    void ExplicitSolverStrategy::RebuildPropertiesProxyPointers(std::vector<SphericParticle*>& rCustomListOfSphericParticles) {
        //This function is called for the local mesh and the ghost mesh, so mListOfSphericElements must not be used here.
        KRATOS_TRY

        const int number_of_particles = (int) rCustomListOfSphericParticles.size();
        std::vector<PropertiesProxy>& vector_of_properties_proxies = PropertiesProxiesManager().GetPropertiesProxies(*mpDem_model_part);

        #pragma omp parallel for
        for (int i = 0; i < number_of_particles; i++) {
          rCustomListOfSphericParticles[i]->SetFastProperties(vector_of_properties_proxies);
        }
        return;
        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::SendProcessInfoToClustersModelPart() {

        KRATOS_TRY

        ProcessInfo& r_process_info = mpDem_model_part->GetProcessInfo();
        ProcessInfo& rClusters_process_info = mpCluster_model_part->GetProcessInfo();

        r_process_info[CONTAINS_CLUSTERS] = false;
        rClusters_process_info[CONTAINS_CLUSTERS] = true;

        rClusters_process_info[GRAVITY] = r_process_info[GRAVITY];
        rClusters_process_info[ROTATION_OPTION] = r_process_info[ROTATION_OPTION];
        rClusters_process_info[DELTA_TIME] = r_process_info[DELTA_TIME];
        rClusters_process_info[VIRTUAL_MASS_OPTION] = r_process_info[VIRTUAL_MASS_OPTION];
        rClusters_process_info[TRIHEDRON_OPTION] = r_process_info[TRIHEDRON_OPTION];
        rClusters_process_info[NODAL_MASS_COEFF] = r_process_info[NODAL_MASS_COEFF];

        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::UpdateMaxIdOfCreatorDestructor() {

        KRATOS_TRY

        int max_Id = mpParticleCreatorDestructor->GetCurrentMaxNodeId();
        ModelPart& r_model_part = GetModelPart();
        int max_DEM_Id = mpParticleCreatorDestructor->FindMaxNodeIdInModelPart(r_model_part);
        int max_FEM_Id = mpParticleCreatorDestructor->FindMaxNodeIdInModelPart(*mpFem_model_part);
        int max_cluster_Id = mpParticleCreatorDestructor->FindMaxNodeIdInModelPart(*mpCluster_model_part);

        max_Id = std::max(max_Id, max_DEM_Id);
        max_Id = std::max(max_Id, max_FEM_Id);
        max_Id = std::max(max_Id, max_cluster_Id);
        mpParticleCreatorDestructor->SetMaxNodeId(max_Id);

        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::RepairPointersToNormalProperties(std::vector<SphericParticle*>& rCustomListOfSphericParticles) {

        KRATOS_TRY

        bool found = false;
        const int number_of_particles = (int) rCustomListOfSphericParticles.size();
        #pragma omp parallel for
        for (int i = 0; i < number_of_particles; i++) {

            int own_properties_id = rCustomListOfSphericParticles[i]->GetProperties().Id();

            for (PropertiesIterator props_it = mpDem_model_part->GetMesh(0).PropertiesBegin(); props_it != mpDem_model_part->GetMesh(0).PropertiesEnd(); props_it++) {
                int model_part_id = props_it->GetId();
                if (own_properties_id == model_part_id) {
                    rCustomListOfSphericParticles[i]->SetProperties(*(props_it.base()));
                    found = true;
                    break;
                }
            }

            if (found) continue;

            for (PropertiesIterator props_it = mpInlet_model_part->GetMesh(0).PropertiesBegin(); props_it != mpInlet_model_part->GetMesh(0).PropertiesEnd(); props_it++) {
                int model_part_id = props_it->GetId();
                if (own_properties_id == model_part_id) {
                    rCustomListOfSphericParticles[i]->SetProperties(*(props_it.base()));
                    found = true;
                    break;
                }
            }

            if (found) continue;

            for (PropertiesIterator props_it = mpCluster_model_part->GetMesh(0).PropertiesBegin(); props_it != mpCluster_model_part->GetMesh(0).PropertiesEnd(); props_it++) {
                int model_part_id = props_it->GetId();
                if (own_properties_id == model_part_id) {
                    rCustomListOfSphericParticles[i]->SetProperties(*(props_it.base()));
                    found = true;
                    break;
                }
            }

            if (!found) KRATOS_THROW_ERROR(std::logic_error, "This particle could not find its properties!!", "");
        }

        KRATOS_CATCH("")
    }


    void ExplicitSolverStrategy::DisplayThreadInfo() {

        ModelPart& r_model_part = GetModelPart();
        KRATOS_INFO("DEM") << "          **************************************************" << std::endl;
        KRATOS_INFO("DEM") << "            Parallelism Info:  MPI number of nodes: " << r_model_part.GetCommunicator().TotalProcesses() << std::endl;
        if (r_model_part.GetCommunicator().TotalProcesses() > 1)
            KRATOS_INFO("DEM") << "            Parallelism Info:  MPI node Id: " << r_model_part.GetCommunicator().MyPID() << std::endl;
        KRATOS_INFO("DEM") << "            Parallelism Info:  OMP number of processors: " << mNumberOfThreads << std::endl;
        KRATOS_INFO("DEM") << "          **************************************************" << std::endl;
        KRATOS_INFO("DEM") << std::endl;

    }

    void ExplicitSolverStrategy::Initialize() {

        KRATOS_TRY

        ModelPart& r_model_part = GetModelPart();

        ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
        SendProcessInfoToClustersModelPart();

        mNumberOfThreads = OpenMPUtils::GetNumThreads();
        DisplayThreadInfo();

        RebuildListOfSphericParticles<SphericParticle>(r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericParticles);
        RebuildListOfSphericParticles<SphericParticle>(r_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostSphericParticles);

        PropertiesProxiesManager().CreatePropertiesProxies(*mpDem_model_part, *mpInlet_model_part, *mpCluster_model_part);

        RepairPointersToNormalProperties(mListOfSphericParticles); // The particles sent to this partition have their own copy of the Kratos properties they were using in the previous partition!!
        RepairPointersToNormalProperties(mListOfGhostSphericParticles);

        RebuildPropertiesProxyPointers(mListOfSphericParticles);
        RebuildPropertiesProxyPointers(mListOfGhostSphericParticles);

        GetSearchControl() = r_process_info[SEARCH_CONTROL];

        InitializeDEMElements();
        InitializeFEMElements();
        UpdateMaxIdOfCreatorDestructor();
        InitializeClusters(); // This adds elements to the balls modelpart

        RebuildListOfSphericParticles<SphericParticle>(r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericParticles);
        RebuildListOfSphericParticles<SphericParticle>(r_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostSphericParticles);

        InitializeSolutionStep();

        ApplyInitialConditions();

        // Search Neighbours and related operations
        SetSearchRadiiOnAllParticles(*mpDem_model_part, mpDem_model_part->GetProcessInfo()[SEARCH_RADIUS_INCREMENT], 1.0);
        SearchNeighbours();
        ComputeNewNeighboursHistoricalData();

        SetSearchRadiiOnAllParticles(*mpDem_model_part, mpDem_model_part->GetProcessInfo()[SEARCH_RADIUS_INCREMENT_FOR_WALLS], 1.0);
        SearchRigidFaceNeighbours(); //initial search is performed with hierarchical method in any case MSI
        ComputeNewRigidFaceNeighboursHistoricalData();

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

        //set flag to 2 (search performed this timestep)
        mSearchControl = 2;

        // Finding overlapping of initial configurations
        if (r_process_info[CLEAN_INDENT_OPTION]) {
            for (int i = 0; i < 10; i++) CalculateInitialMaxIndentations(r_process_info);
        }

        if (r_process_info[CRITICAL_TIME_OPTION]) {
            //InitialTimeStepCalculation();   //obsolete call
            CalculateMaxTimeStep();
        }

        r_process_info[PARTICLE_INELASTIC_FRICTIONAL_ENERGY] = 0.0;

        // 5. Finalize Solution Step.
        //FinalizeSolutionStep();
        ComputeNodalArea();

        KRATOS_CATCH("")
    }// Initialize()

    void ExplicitSolverStrategy::MarkToDeleteAllSpheresInitiallyIndentedWithFEM(ModelPart& rSpheresModelPart) {

        ElementsArrayType& pElements = rSpheresModelPart.GetCommunicator().LocalMesh().Elements();

        #pragma omp parallel for
        for (int k = 0; k < (int)pElements.size(); k++) {
            ElementsArrayType::iterator it = pElements.ptr_begin() + k;
            Element* p_element = &(*it);
            SphericParticle* p_sphere = dynamic_cast<SphericParticle*>(p_element);

            if (p_sphere->mNeighbourRigidFaces.size()) {
                p_sphere->Set(TO_ERASE);
                p_sphere->GetGeometry()[0].Set(TO_ERASE);
            }
        }
    }

    void ExplicitSolverStrategy::ComputeNodalArea() {

        KRATOS_TRY

        ModelPart& fem_model_part = GetFemModelPart();
        NodesArrayType& pNodes = fem_model_part.Nodes();
        NodesArrayType::iterator i_begin = pNodes.ptr_begin();
        NodesArrayType::iterator i_end = pNodes.ptr_end();

        for (ModelPart::NodeIterator i = i_begin; i != i_end; ++i) {
            double& node_area = i->GetSolutionStepValue(DEM_NODAL_AREA);
            node_area = 0.0;
        }

        ConditionsArrayType& pConditions = GetFemModelPart().GetCommunicator().LocalMesh().Conditions();
        ConditionsArrayType::iterator it_begin = pConditions.ptr_begin();
        ConditionsArrayType::iterator it_end = pConditions.ptr_end();

        for (ConditionsArrayType::iterator it = it_begin; it != it_end; ++it) { //each iteration refers to a different triangle or quadrilateral

            Condition::GeometryType& geometry = it->GetGeometry();
            double Element_Area = geometry.Area();

            for (unsigned int i = 0; i < geometry.size(); i++) { //talking about each of the three nodes of the condition
                double& node_area = geometry[i].FastGetSolutionStepValue(DEM_NODAL_AREA);
                node_area += 0.333333333333333 * Element_Area; //TODO: ONLY FOR TRIANGLE... Generalize for 3 or 4 nodes
            }
        }

        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::CalculateMaxTimeStep() {
        KRATOS_TRY
        ModelPart& r_model_part = GetModelPart();
        ProcessInfo& r_process_info = r_model_part.GetProcessInfo();

        bool has_mpi = false; //check MPI not available in this strategy. refer to continuum strategy
        //          Check_MPI(has_mpi);

        std::vector<double> thread_maxima(OpenMPUtils::GetNumThreads(), 0.0);
        const int number_of_particles = (int) mListOfSphericParticles.size();

        #pragma omp parallel for
        for (int i = 0; i < number_of_particles; i++) {
            double max_sqr_period = mListOfSphericParticles[i]->CalculateLocalMaxPeriod(has_mpi, r_process_info);
            if (max_sqr_period > thread_maxima[OpenMPUtils::ThisThread()]) thread_maxima[OpenMPUtils::ThisThread()] = max_sqr_period;
        }

        double max_across_threads = 0.0;
        for (int i = 0; i < OpenMPUtils::GetNumThreads(); i++) {
            if (thread_maxima[i] > max_across_threads) max_across_threads = thread_maxima[i];
        }

        double critical_period = sqrt(max_across_threads);
        double beta = 0.03;
        double critical_timestep = beta * Globals::Pi / critical_period;

        double t = CalculateMaxInletTimeStep();
        if (t<critical_timestep && t>0.0){critical_timestep = t;}

        r_process_info[DELTA_TIME] = critical_timestep;
        KRATOS_INFO("DEM") << " (Critical) Timestep set to " << critical_timestep << ". " << "\n" << std::endl;


       //PropertiesContainerType pprop1 = *mpInlet_model_part->pProperties();
       //double young = (*mpInlet_model_part)[YOUNG_MODULUS];  // no funciona pq no forma part de modelpart sino de properties
       //PropertiesContainerType pprop2 = mpInlet_model_part->PropertiesArray(0);
       //long unsigned int pprop4 = mpInlet_model_part->NumberOfSubModelParts();
       KRATOS_CATCH("")
    }

    double ExplicitSolverStrategy::CalculateMaxInletTimeStep() {
        for (PropertiesIterator props_it = mpInlet_model_part->GetMesh(0).PropertiesBegin(); props_it != mpInlet_model_part->GetMesh(0).PropertiesEnd(); props_it++) {
            if ((*props_it).Has(PARTICLE_DENSITY)) {
                int inlet_prop_id = props_it->GetId();
                double young = (*props_it)[YOUNG_MODULUS];
                double density = (*props_it)[PARTICLE_DENSITY];
                double poisson = (*props_it)[POISSON_RATIO];

                for (ModelPart::SubModelPartsContainerType::iterator sub_model_part = mpInlet_model_part->SubModelPartsBegin(); sub_model_part != mpInlet_model_part->SubModelPartsEnd(); ++sub_model_part) {
                    int smp_prop_id = (*sub_model_part)[PROPERTIES_ID];
                    if (smp_prop_id == inlet_prop_id) {
                        double radius = (*sub_model_part)[RADIUS];
                        double shear_modulus = young/(2.0*(1.0+poisson));
                        double t = (Globals::Pi*radius*sqrt(density/shear_modulus))/(0.1630*poisson+0.8766);
                        return t;
                    }
                }
            }
        }
        return 0.0;
    }

    void ExplicitSolverStrategy::InitializeClusters() {
        KRATOS_TRY
        ElementsArrayType& pElements = mpCluster_model_part->GetCommunicator().LocalMesh().Elements();
        const int number_of_clusters = pElements.size();
        ProcessInfo& r_process_info = GetModelPart().GetProcessInfo();
        bool continuum_strategy = r_process_info[CONTINUUM_OPTION];
        std::vector<PropertiesProxy>& vector_of_properties_proxies = PropertiesProxiesManager().GetPropertiesProxies(*mpDem_model_part);

        //mpParticleCreatorDestructor->FindAndSaveMaxNodeIdInModelPart(*mpDem_model_part); //This has been moved to python main script and checks both dem model part and walls model part (also important!)

        #pragma omp parallel for schedule(dynamic, 100)
        for (int k = 0; k < number_of_clusters; k++) {

            ElementsArrayType::iterator it = pElements.ptr_begin() + k;
            Cluster3D& cluster_element = dynamic_cast<Kratos::Cluster3D&> (*it);

            cluster_element.Initialize(r_process_info);

            PropertiesProxy* p_fast_properties = NULL;
            int general_properties_id = cluster_element.GetProperties().Id();
            for (unsigned int i = 0; i < vector_of_properties_proxies.size(); i++) {
                int fast_properties_id = vector_of_properties_proxies[i].GetId();
                if (fast_properties_id == general_properties_id) {
                    p_fast_properties = &(vector_of_properties_proxies[i]);
                    break;
                }
            }
            cluster_element.CreateParticles(mpParticleCreatorDestructor.get(), *mpDem_model_part, p_fast_properties, continuum_strategy);
        }
        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::GetClustersForce() {
        KRATOS_TRY
        ProcessInfo& r_process_info = GetModelPart().GetProcessInfo(); //Getting the Process Info of the Balls ModelPart!
        const array_1d<double, 3>& gravity = r_process_info[GRAVITY];

        ElementsArrayType& pElements = mpCluster_model_part->GetCommunicator().LocalMesh().Elements();
        const int number_of_clusters = pElements.size();

        #pragma omp parallel for schedule(dynamic, 50)
        for (int k = 0; k < number_of_clusters; k++) {

            ElementsArrayType::iterator it = pElements.ptr_begin() + k;
            Cluster3D& cluster_element = dynamic_cast<Kratos::Cluster3D&> (*it);

            cluster_element.GetGeometry()[0].FastGetSolutionStepValue(TOTAL_FORCES).clear();
            cluster_element.GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT).clear();

            cluster_element.GetClustersForce(gravity);
        } // loop over clusters
        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::GetRigidBodyElementsForce() {
        KRATOS_TRY
        CalculateConditionsRHSAndAdd();
        ProcessInfo& r_process_info = GetModelPart().GetProcessInfo(); //Getting the Process Info of the Balls ModelPart!
        const array_1d<double, 3>& gravity = r_process_info[GRAVITY];
        ModelPart& fem_model_part = GetFemModelPart();
        ElementsArrayType& pElements = fem_model_part.GetCommunicator().LocalMesh().Elements();
        const int number_of_rigid_body_elements = pElements.size();

        //DO NOT PARALLELIZE THIS LOOP, IT IS PARALLELIZED INSIDE
        for (int k = 0; k < number_of_rigid_body_elements; k++) {

            ElementsArrayType::iterator it = pElements.ptr_begin() + k;
            RigidBodyElement3D& rigid_body_element = dynamic_cast<Kratos::RigidBodyElement3D&> (*it);
            rigid_body_element.GetGeometry()[0].FastGetSolutionStepValue(TOTAL_FORCES).clear();
            rigid_body_element.GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT).clear();
            rigid_body_element.GetRigidBodyElementsForce(gravity);

        } // loop over rigid body elements

        KRATOS_CATCH("")
    }

    double ExplicitSolverStrategy::Solve() {
        KRATOS_TRY
        ModelPart& r_model_part = GetModelPart();

        InitializeSolutionStep();
        SearchDEMOperations(r_model_part);
        SearchFEMOperations(r_model_part);
        ForceOperations(r_model_part);
        PerformTimeIntegrationOfMotion();
        FinalizeSolutionStep();

        return 0.00;

        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::SearchDEMOperations(ModelPart& r_model_part, bool has_mpi) {
        KRATOS_TRY
        ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
        int time_step = r_process_info[TIME_STEPS];
        const double time = r_process_info[TIME];
        const bool is_time_to_search_neighbours = (time_step + 1) % mNStepSearch == 0 && (time_step > 0); //Neighboring search. Every N times.
        const bool is_time_to_mark_and_remove = is_time_to_search_neighbours && (r_process_info[BOUNDING_BOX_OPTION] && time >= r_process_info[BOUNDING_BOX_START_TIME] && time <= r_process_info[BOUNDING_BOX_STOP_TIME]);
        BoundingBoxUtility(is_time_to_mark_and_remove);

        if (is_time_to_search_neighbours) {
            if (!is_time_to_mark_and_remove) { //Just in case that some entities were marked as TO_ERASE without a bounding box (manual removal)
                mpParticleCreatorDestructor->DestroyParticles(*mpCluster_model_part);
                mpParticleCreatorDestructor->DestroyParticles(r_model_part);
            }

            RebuildListOfSphericParticles<SphericParticle>(r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericParticles);
            RebuildListOfSphericParticles<SphericParticle>(r_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostSphericParticles);

            SearchNeighbours();

            RebuildListOfSphericParticles <SphericParticle> (r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericParticles);
            RebuildListOfSphericParticles <SphericParticle> (r_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostSphericParticles);
            RepairPointersToNormalProperties(mListOfSphericParticles);
            RepairPointersToNormalProperties(mListOfGhostSphericParticles);
            RebuildPropertiesProxyPointers(mListOfSphericParticles);
            RebuildPropertiesProxyPointers(mListOfGhostSphericParticles);

            ComputeNewNeighboursHistoricalData();

            mSearchControl = 2; // Search is active and has been performed during this time step
            //ReorderParticles();
        } else {
            mSearchControl = 1; // Search is active but no search has been done this time step;
        }

        //RebuildPropertiesProxyPointers(mListOfSphericParticles);
        //RebuildPropertiesProxyPointers(mListOfGhostSphericParticles);
        KRATOS_CATCH("")
    }//SearchDEMOperations;

    void ExplicitSolverStrategy::SearchFEMOperations(ModelPart& r_model_part, bool has_mpi) {
        KRATOS_TRY
        ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
        int time_step = r_process_info[TIME_STEPS];
        const bool is_time_to_search_neighbours = (time_step + 1) % mNStepSearch == 0 && (time_step > 0); //Neighboring search. Every N times.

        if (is_time_to_search_neighbours) { // for the moment it is always true, until all issues have been solved
            SearchRigidFaceNeighbours();
            ComputeNewRigidFaceNeighboursHistoricalData();
            mSearchControl = 2; // Search is active and has been performed during this time step
        }

        else {
            ConditionsArrayType& pTConditions = mpFem_model_part->GetCommunicator().LocalMesh().Conditions();
            const int number_of_conditions = (int) pTConditions.size();

            if (number_of_conditions > 0) {
                CheckHierarchyWithCurrentNeighbours();
                ComputeNewRigidFaceNeighboursHistoricalData();
                mSearchControl = 1; // Search is active but no search has been done this time step;
            }
        }
        KRATOS_CATCH("")
    }//SearchFEMOperations

    void ExplicitSolverStrategy::ForceOperations(ModelPart& r_model_part) {

        KRATOS_TRY
        // 3. Get and Calculate the forces
        CleanEnergies();
        GetForce(); // Basically only calls CalculateRightHandSide()
        //FastGetForce();
        GetClustersForce();
        GetRigidBodyElementsForce();

        if (r_model_part.GetProcessInfo()[COMPUTE_FEM_RESULTS_OPTION]) {
            CalculateNodalPressuresAndStressesOnWalls();
        }

        // 4. Synchronize (should be just FORCE and TORQUE)
        SynchronizeRHS(r_model_part);

        KRATOS_CATCH("")
    }//ForceOperations;

    void ExplicitSolverStrategy::InitialTimeStepCalculation() // obsoleta delete
    {
        KRATOS_TRY
        ModelPart& r_model_part = GetModelPart();
        ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
        ElementsArrayType& pElements = r_model_part.GetCommunicator().LocalMesh().Elements();

        ElementsIterator it_begin = pElements.ptr_begin();
        ElementsIterator it_end = pElements.ptr_end();

        double& process_info_delta_time = r_process_info[DELTA_TIME];
        process_info_delta_time = mMaxTimeStep;
        double temp_time_step = std::numeric_limits<double>::infinity();
        double elem_critical_time_step = temp_time_step;

        for (ElementsIterator it = it_begin; it != it_end; it++) {
            it->Calculate(DELTA_TIME, elem_critical_time_step, r_process_info);

            if (elem_critical_time_step < temp_time_step) {
                temp_time_step = elem_critical_time_step;
            }

        }

        temp_time_step /= mSafetyFactor;

        if (temp_time_step < mMaxTimeStep) process_info_delta_time = temp_time_step;

        KRATOS_INFO("DEM") << std::scientific;
        KRATOS_INFO("DEM") << std::setprecision(3) << "************* Using " << process_info_delta_time << " time step. (Critical: "
                  << temp_time_step << " with a diving factor: " << mSafetyFactor << " ) *************" << "\n" << std::endl;
        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::GetForce() {

        KRATOS_TRY

        ProcessInfo& r_process_info = GetModelPart().GetProcessInfo();
        double dt = r_process_info[DELTA_TIME];
        const array_1d<double, 3>& gravity = r_process_info[GRAVITY];

        const int number_of_particles = (int) mListOfSphericParticles.size();

        #pragma omp parallel for schedule(dynamic, 100)
        for (int i = 0; i < number_of_particles; i++) {
            mListOfSphericParticles[i]->CalculateRightHandSide(r_process_info, dt, gravity, mSearchControl);
        }

        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::FastGetForce() {
        KRATOS_TRY
        ProcessInfo& r_process_info = GetModelPart().GetProcessInfo();
        double dt = r_process_info[DELTA_TIME];
        const array_1d<double, 3>& gravity = r_process_info[GRAVITY];
        const int number_of_particles = (int) mListOfSphericParticles.size();

        #pragma omp parallel
        {
            #pragma omp for
            for (int i = 0; i < number_of_particles; i++) {
                mListOfSphericParticles[i]->FirstCalculateRightHandSide(r_process_info, dt, mSearchControl);
            }
            #pragma omp for
            for (int i = 0; i < number_of_particles; i++) {
                mListOfSphericParticles[i]->CollectCalculateRightHandSide(r_process_info);
            }
            #pragma omp for
            for (int i = 0; i < number_of_particles; i++) {
                mListOfSphericParticles[i]->FinalCalculateRightHandSide(r_process_info, dt, gravity);
            }
        }

        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::PerformTimeIntegrationOfMotion(int StepFlag) {

        KRATOS_TRY

        ProcessInfo& r_process_info = GetModelPart().GetProcessInfo();
        double delta_t = r_process_info[DELTA_TIME];
        double virtual_mass_coeff = r_process_info[NODAL_MASS_COEFF]; //TODO: change the name of this variable to FORCE_REDUCTION_FACTOR
        bool virtual_mass_option = (bool) r_process_info[VIRTUAL_MASS_OPTION];
        double force_reduction_factor = 1.0;
        if (virtual_mass_option) {
            force_reduction_factor = virtual_mass_coeff;
            if ((force_reduction_factor > 1.0) || (force_reduction_factor < 0.0)) {
                KRATOS_THROW_ERROR(std::runtime_error, "The force reduction factor is either larger than 1 or negative: FORCE_REDUCTION_FACTOR= ", virtual_mass_coeff)
            }
        }

        bool rotation_option = r_process_info[ROTATION_OPTION];

        const int number_of_particles       = (int) mListOfSphericParticles.size();
        const int number_of_ghost_particles = (int) mListOfGhostSphericParticles.size();

        ModelPart& r_clusters_model_part  = *mpCluster_model_part;
        ElementsArrayType& pLocalClusters = r_clusters_model_part.GetCommunicator().LocalMesh().Elements();
        ElementsArrayType& pGhostClusters = r_clusters_model_part.GetCommunicator().GhostMesh().Elements();
        ModelPart& r_fem_model_part  = *mpFem_model_part;
        ElementsArrayType& pFemElements = r_fem_model_part.GetCommunicator().LocalMesh().Elements();

        #pragma omp parallel
        {
            #pragma omp for nowait
            for (int i = 0; i < number_of_particles; i++) {
                mListOfSphericParticles[i]->Move(delta_t, rotation_option, force_reduction_factor, StepFlag);
            }

            #pragma omp for nowait
            for (int i = 0; i < number_of_ghost_particles; i++) {
                mListOfGhostSphericParticles[i]->Move(delta_t, rotation_option, force_reduction_factor, StepFlag);
            }

            #pragma omp for nowait
            for (int k = 0; k < (int) pLocalClusters.size(); k++) {
                ElementsArrayType::iterator it = pLocalClusters.ptr_begin() + k;
                Cluster3D& cluster_element = dynamic_cast<Kratos::Cluster3D&> (*it);
                cluster_element.RigidBodyElement3D::Move(delta_t, rotation_option, force_reduction_factor, StepFlag);
            }

            #pragma omp for nowait
            for (int k = 0; k < (int) pGhostClusters.size(); k++) {
                 ElementsArrayType::iterator it = pGhostClusters.ptr_begin() + k;
                Cluster3D& cluster_element = dynamic_cast<Kratos::Cluster3D&> (*it);
                cluster_element.RigidBodyElement3D::Move(delta_t, rotation_option, force_reduction_factor, StepFlag);
            }

            #pragma omp for nowait
            for (int k = 0; k < (int) pFemElements.size(); k++) {
                ElementsArrayType::iterator it = pFemElements.ptr_begin() + k;
                RigidBodyElement3D& rigid_body_element = dynamic_cast<Kratos::RigidBodyElement3D&> (*it);
                rigid_body_element.Move(delta_t, rotation_option, force_reduction_factor, StepFlag);
            }
        }

        //GetScheme()->Calculate(GetModelPart(), StepFlag);
        //GetScheme()->Calculate(*mpCluster_model_part, StepFlag);
        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::InitializeSolutionStep() {
        KRATOS_TRY

        ModelPart& r_model_part = GetModelPart();
        ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
        ElementsArrayType& pElements = r_model_part.GetCommunicator().LocalMesh().Elements();

        ModelPart& r_fem_model_part = GetFemModelPart();
        ProcessInfo& r_fem_process_info = r_fem_model_part.GetProcessInfo();
        ConditionsArrayType& pConditions = r_fem_model_part.GetCommunicator().LocalMesh().Conditions();

        #pragma omp parallel
        {
            #pragma omp for nowait
            for (int k = 0; k < (int) pElements.size(); k++) {
                ElementsArrayType::iterator it = pElements.ptr_begin() + k;
                (it)->InitializeSolutionStep(r_process_info);
            }

            #pragma omp for nowait
            for (int k = 0; k < (int) pConditions.size(); k++) {
                ConditionsArrayType::iterator it = pConditions.ptr_begin() + k;
                (it)->InitializeSolutionStep(r_fem_process_info);
            }
        }

        ApplyPrescribedBoundaryConditions();
        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::BoundingBoxUtility(bool is_time_to_mark_and_remove) {
        KRATOS_TRY
        ModelPart& r_model_part = GetModelPart();
        ProcessInfo& r_process_info = r_model_part.GetProcessInfo();

        if (r_process_info[DOMAIN_IS_PERIODIC]) {
            mpParticleCreatorDestructor->MoveParticlesOutsideBoundingBoxBackInside(r_model_part);
        } else if (is_time_to_mark_and_remove) {
            mpParticleCreatorDestructor->DestroyParticlesOutsideBoundingBox(*mpCluster_model_part);
            mpParticleCreatorDestructor->DestroyParticlesOutsideBoundingBox(r_model_part);
        }
        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::FinalizeSolutionStep() {

        KRATOS_TRY
        ModelPart& r_model_part = GetModelPart();
        ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
        ElementsArrayType& pElements = r_model_part.GetCommunicator().LocalMesh().Elements();
        OpenMPUtils::CreatePartition(mNumberOfThreads, pElements.size(), this->GetElementPartition());

         #pragma omp parallel for if(pElements.size()>100)
        for (int k = 0; k < (int)pElements.size(); k++) {
            ElementsArrayType::iterator it = pElements.ptr_begin() + k;
            (it)->FinalizeSolutionStep(r_process_info); //we use this function to call the set initial contacts and the add continuum contacts
        } //loop over particles

        //if (true) AuxiliaryFunctions::ComputeReactionOnTopAndBottomSpheres(r_model_part);


        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::InitializeElements() {
        KRATOS_TRY
        ModelPart& r_model_part = GetModelPart();
        ElementsArrayType& pElements = r_model_part.GetCommunicator().LocalMesh().Elements();

        OpenMPUtils::CreatePartition(mNumberOfThreads, pElements.size(), this->GetElementPartition());

        #pragma omp parallel for
        for (int k = 0; k < mNumberOfThreads; k++) {
            ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
            ElementsArrayType::iterator it_end = pElements.ptr_begin() + this->GetElementPartition()[k + 1];

            for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it) {
                (it)->Initialize();
            }
        }
        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::InitializeDEMElements() {

        KRATOS_TRY

        ModelPart& r_model_part = GetModelPart();
        ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
        const int number_of_particles = (int) mListOfSphericParticles.size();
        double total_mass = 0.0;

        #pragma omp parallel for
        for (int i = 0; i < number_of_particles; i++) {
            mListOfSphericParticles[i]->Initialize(r_process_info);
            total_mass += mListOfSphericParticles[i]->GetMass();
        }

        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::InitializeFEMElements() {

        KRATOS_TRY

        ConditionsArrayType& pTConditions = mpFem_model_part->GetCommunicator().LocalMesh().Conditions();
        ModelPart& fem_model_part = GetFemModelPart();
        ProcessInfo& r_process_info = GetModelPart().GetProcessInfo();

        if (fem_model_part.NumberOfSubModelParts()) {
            for (ModelPart::SubModelPartsContainerType::iterator sub_model_part = fem_model_part.SubModelPartsBegin(); sub_model_part != fem_model_part.SubModelPartsEnd(); ++sub_model_part) {

                ModelPart& submp = *sub_model_part;
                NodesArrayType& pNodes = sub_model_part->Nodes();

                #pragma omp parallel for
                for (int i=0; i<(int)pTConditions.size(); i++) {
                    ConditionsArrayType::iterator it = pTConditions.ptr_begin() + i;
                    (it)->Initialize();
                }

                // Central Node
                Node<3>::Pointer central_node;
                Geometry<Node<3> >::PointsArrayType central_node_list;

                array_1d<double, 3> reference_coordinates = ZeroVector(3);

                if (submp.Has(RIGID_BODY_CENTER_OF_MASS)) {
                    reference_coordinates[0] = submp[RIGID_BODY_CENTER_OF_MASS][0];
                    reference_coordinates[1] = submp[RIGID_BODY_CENTER_OF_MASS][1];
                    reference_coordinates[2] = submp[RIGID_BODY_CENTER_OF_MASS][2];
                }

                int Node_Id_1 = mpParticleCreatorDestructor->FindMaxNodeIdInModelPart(fem_model_part);

                mpParticleCreatorDestructor->CentroidCreatorForRigidBodyElements(fem_model_part, central_node, Node_Id_1 + 1, reference_coordinates);

                central_node_list.push_back(central_node);

                int Element_Id_1 = mpParticleCreatorDestructor->FindMaxElementIdInModelPart(fem_model_part);

                Properties::Pointer properties;
                if (submp.Has(PROPERTIES_ID)) properties = fem_model_part.GetMesh().pGetProperties(submp[PROPERTIES_ID]);
                else properties = fem_model_part.GetMesh().pGetProperties(fem_model_part.GetMesh(0).PropertiesBegin()->GetId()); // JIG: Backward compatibility, it should be removed in the future

                std::string ElementNameString = "RigidBodyElement3D";

                if (submp.Has(FLOATING_OPTION)) {
                    if (submp[FLOATING_OPTION]) {
                        ElementNameString = "ShipElement3D";
                    }
                }

                const Element& r_reference_element = KratosComponents<Element>::Get(ElementNameString);
                Element::Pointer RigidBodyElement3D_Kratos = r_reference_element.Create(Element_Id_1 + 1, central_node_list, properties);
                RigidBodyElement3D* rigid_body_element = dynamic_cast<RigidBodyElement3D*>(RigidBodyElement3D_Kratos.get());

                fem_model_part.AddElement(RigidBodyElement3D_Kratos); //, Element_Id + 1);

                std::size_t element_id = Element_Id_1 + 1;
                std::vector<std::size_t> ElementIds;
                ElementIds.push_back(element_id);

                if (submp.Has(FREE_BODY_MOTION)) { // JIG: Backward compatibility, it should be removed in the future
                    if (submp[FREE_BODY_MOTION]) {

                        std::vector<std::vector<Node<3>::Pointer> > thread_vectors_of_node_pointers;
                        thread_vectors_of_node_pointers.resize(mNumberOfThreads);
                        std::vector<std::vector<array_1d<double, 3> > > thread_vectors_of_coordinates;
                        thread_vectors_of_coordinates.resize(mNumberOfThreads);

                        #pragma omp parallel for
                        for (int k = 0; k < (int)pNodes.size(); k++) {
                            ModelPart::NodeIterator i = pNodes.ptr_begin() + k;
                            thread_vectors_of_node_pointers[OpenMPUtils::ThisThread()].push_back(*(i.base())); //TODO: this could be raw pointers. It would be a lot faster here (same speed when reading later on)
                            thread_vectors_of_coordinates[OpenMPUtils::ThisThread()].push_back(i->Coordinates() - reference_coordinates);
                        }
                        for (int i = 0; i < mNumberOfThreads; i++) {
                            rigid_body_element->mListOfNodes.insert(rigid_body_element->mListOfNodes.end(), thread_vectors_of_node_pointers[i].begin(), thread_vectors_of_node_pointers[i].end());
                            rigid_body_element->mListOfCoordinates.insert(rigid_body_element->mListOfCoordinates.end(), thread_vectors_of_coordinates[i].begin(), thread_vectors_of_coordinates[i].end());
                        }

                        std::vector<std::vector<RigidFace3D*> > thread_vectors_of_rigid_faces;
                        thread_vectors_of_rigid_faces.resize(mNumberOfThreads);

                        #pragma omp parallel for
                        for (int k = 0; k < (int)pTConditions.size(); k++) {
                            ConditionsArrayType::iterator it = pTConditions.ptr_begin() + k;
                            RigidFace3D* it_face = dynamic_cast<RigidFace3D*>(&(*it));
                            thread_vectors_of_rigid_faces[OpenMPUtils::ThisThread()].push_back(it_face);
                        }
                        for (int i = 0; i < mNumberOfThreads; i++) {
                            rigid_body_element->mListOfRigidFaces.insert(rigid_body_element->mListOfRigidFaces.end(), thread_vectors_of_rigid_faces[i].begin(), thread_vectors_of_rigid_faces[i].end());
                        }
                    }
                }

                rigid_body_element->Initialize(r_process_info);
                rigid_body_element->CustomInitialize(submp);
            }
        }

        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::CalculateConditionsRHSAndAdd() {

        KRATOS_TRY
        ClearFEMForces();
        ConditionsArrayType& pConditions = GetFemModelPart().GetCommunicator().LocalMesh().Conditions();
        ProcessInfo& r_process_info = GetFemModelPart().GetProcessInfo();

        Vector rhs_cond;
        Vector rhs_cond_elas;
        DenseVector<unsigned int> condition_partition;
        OpenMPUtils::CreatePartition(mNumberOfThreads, pConditions.size(), condition_partition);
        unsigned int index;

        #pragma omp parallel for private (index, rhs_cond, rhs_cond_elas)
        for (int k = 0; k < mNumberOfThreads; k++) {

            ConditionsArrayType::iterator it_begin = pConditions.ptr_begin() + condition_partition[k];
            ConditionsArrayType::iterator it_end = pConditions.ptr_begin() + condition_partition[k + 1];

            for (ConditionsArrayType::iterator it = it_begin; it != it_end; ++it) { //each iteration refers to a different triangle or quadrilateral

                Condition::GeometryType& geom = it->GetGeometry();

                //double Element_Area = geom.Area();

                it->CalculateRightHandSide(rhs_cond, r_process_info);
                DEMWall* p_wall = dynamic_cast<DEMWall*> (&(*it));
                p_wall->CalculateElasticForces(rhs_cond_elas, r_process_info);
                array_1d<double, 3> Normal_to_Element = ZeroVector(3);

                if (geom.size()>2) p_wall->CalculateNormal(Normal_to_Element);

                const unsigned int& dim = geom.WorkingSpaceDimension();

                for (unsigned int i = 0; i < geom.size(); i++) { //talking about each of the three nodes of the condition
                    //we are studying a certain condition here
                    index = i * dim; //*2;

                    array_1d<double, 3>& node_rhs = geom[i].FastGetSolutionStepValue(CONTACT_FORCES);
                    array_1d<double, 3>& node_rhs_elas = geom[i].FastGetSolutionStepValue(ELASTIC_FORCES);
                    array_1d<double, 3>& node_rhs_tang = geom[i].FastGetSolutionStepValue(TANGENTIAL_ELASTIC_FORCES);
                    double& node_pressure = geom[i].FastGetSolutionStepValue(DEM_PRESSURE);
                    array_1d<double, 3> rhs_cond_comp;

                    geom[i].SetLock();

                    for (unsigned int j = 0; j < dim; j++) { //talking about each coordinate x, y and z, loop on them
                        node_rhs[j] += rhs_cond[index + j];
                        node_rhs_elas[j] += rhs_cond_elas[index + j];
                        rhs_cond_comp[j] = rhs_cond[index + j];
                    }
                    //node_area += 0.333333333333333 * Element_Area; //TODO: ONLY FOR TRIANGLE... Generalize for 3 or 4 nodes.
                    //node_pressure actually refers to normal force. Pressure is actually computed later in function Calculate_Nodal_Pressures_and_Stresses()
                    node_pressure += MathUtils<double>::Abs(GeometryFunctions::DotProduct(rhs_cond_comp, Normal_to_Element));
                    noalias(node_rhs_tang) += rhs_cond_comp - GeometryFunctions::DotProduct(rhs_cond_comp, Normal_to_Element) * Normal_to_Element;

                    geom[i].UnSetLock();
                }
            }
        }
        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::ClearFEMForces() {

        KRATOS_TRY
        ModelPart& fem_model_part = GetFemModelPart();
        NodesArrayType& pNodes = fem_model_part.Nodes();

        #pragma omp parallel for
        for (int k = 0; k < (int)pNodes.size(); k++) {
            ModelPart::NodeIterator i = pNodes.begin() + k;

            array_1d<double, 3>& node_rhs = i->FastGetSolutionStepValue(CONTACT_FORCES);
            array_1d<double, 3>& node_rhs_elas = i->FastGetSolutionStepValue(ELASTIC_FORCES);
            array_1d<double, 3>& node_rhs_tang = i->FastGetSolutionStepValue(TANGENTIAL_ELASTIC_FORCES);
            double& node_pressure = i->GetSolutionStepValue(DEM_PRESSURE);
            //double& node_area = i->GetSolutionStepValue(DEM_NODAL_AREA);
            double& shear_stress = i->FastGetSolutionStepValue(SHEAR_STRESS);

            noalias(node_rhs) = ZeroVector(3);
            noalias(node_rhs_elas) = ZeroVector(3);
            noalias(node_rhs_tang) = ZeroVector(3);
            node_pressure = 0.0;
            //node_area = 0.0;
            shear_stress = 0.0;
        }
        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::CalculateNodalPressuresAndStressesOnWalls() {
        KRATOS_TRY

        ModelPart& fem_model_part = GetFemModelPart();
        NodesArrayType& pNodes = fem_model_part.Nodes();

        #pragma omp parallel for
        for (int k = 0; k < (int)pNodes.size(); k++) {
            ModelPart::NodeIterator i = pNodes.begin() + k;

            double& node_pressure = i->FastGetSolutionStepValue(DEM_PRESSURE);
            double node_area = i->FastGetSolutionStepValue(DEM_NODAL_AREA);
            double& shear_stress = i->FastGetSolutionStepValue(SHEAR_STRESS);
            array_1d<double, 3>& node_rhs_tang = i->FastGetSolutionStepValue(TANGENTIAL_ELASTIC_FORCES);

            node_pressure = node_pressure / node_area;
            shear_stress = GeometryFunctions::module(node_rhs_tang) / node_area;
        }
        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::SetFlagAndVariableToNodes(const Kratos::Flags& r_flag_name, ComponentOf3ComponentsVariableType& r_variable_to_set, const double value, NodesArrayType& r_nodes_array) {
        KRATOS_TRY
        #pragma omp parallel for
        for (int i = 0; i < (int) r_nodes_array.size(); i++) {
            NodesArrayType::iterator node_i = r_nodes_array.ptr_begin() + i;
            node_i->FastGetSolutionStepValue(r_variable_to_set) = value;
            node_i->Set(r_flag_name, true);
        }
        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::SetVariableToNodes(ComponentOf3ComponentsVariableType& r_variable_to_set, const double value, NodesArrayType& r_nodes_array) {
        KRATOS_TRY
        #pragma omp parallel for
        for (int i = 0; i < (int) r_nodes_array.size(); i++) {
            NodesArrayType::iterator node_i = r_nodes_array.ptr_begin() + i;
            node_i->FastGetSolutionStepValue(r_variable_to_set) = value;
        }
        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::ResetPrescribedMotionFlagsRespectingImposedDofs() {
        KRATOS_TRY
        ModelPart& r_model_part = GetModelPart();

        NodesArrayType& r_model_part_nodes = r_model_part.Nodes();

        if (!r_model_part_nodes.size()) return;

        const unsigned int vel_x_dof_position = (r_model_part.NodesBegin())->GetDofPosition(VELOCITY_X);
        const unsigned int ang_vel_x_dof_position = (r_model_part.NodesBegin())->GetDofPosition(ANGULAR_VELOCITY_X);

        #pragma omp parallel for
        for (int i = 0; i < (int) r_model_part_nodes.size(); i++) {
            ModelPart::NodesContainerType::iterator node_i = r_model_part.NodesBegin() + i;
            if (node_i->Is(BLOCKED)) continue;
            Node<3>& node = *node_i;

            if (node.GetDof(VELOCITY_X, vel_x_dof_position).IsFixed()) {
                node.Set(DEMFlags::FIXED_VEL_X, true);
            } else {
                node.Set(DEMFlags::FIXED_VEL_X, false);
            }
            if (node.GetDof(VELOCITY_Y, vel_x_dof_position + 1).IsFixed()) {
                node.Set(DEMFlags::FIXED_VEL_Y, true);
            } else {
                node.Set(DEMFlags::FIXED_VEL_Y, false);
            }
            if (node.GetDof(VELOCITY_Z, vel_x_dof_position + 2).IsFixed()) {
                node.Set(DEMFlags::FIXED_VEL_Z, true);
            } else {
                node.Set(DEMFlags::FIXED_VEL_Z, false);
            }
            if (node.GetDof(ANGULAR_VELOCITY_X, ang_vel_x_dof_position).IsFixed()) {
                node.Set(DEMFlags::FIXED_ANG_VEL_X, true);
            } else {
                node.Set(DEMFlags::FIXED_ANG_VEL_X, false);
            }
            if (node.GetDof(ANGULAR_VELOCITY_Y, ang_vel_x_dof_position + 1).IsFixed()) {
                node.Set(DEMFlags::FIXED_ANG_VEL_Y, true);
            } else {
                node.Set(DEMFlags::FIXED_ANG_VEL_Y, false);
            }
            if (node.GetDof(ANGULAR_VELOCITY_Z, ang_vel_x_dof_position + 2).IsFixed()) {
                node.Set(DEMFlags::FIXED_ANG_VEL_Z, true);
            } else {
                node.Set(DEMFlags::FIXED_ANG_VEL_Z, false);
            }
        }
        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::ApplyPrescribedBoundaryConditions() {

        KRATOS_TRY

        ModelPart& r_model_part = GetModelPart();
        const ProcessInfo& r_process_info = GetModelPart().GetProcessInfo();
        const double time = r_process_info[TIME];

        for (ModelPart::SubModelPartsContainerType::iterator sub_model_part = r_model_part.SubModelPartsBegin(); sub_model_part != r_model_part.SubModelPartsEnd(); ++sub_model_part) {

            double vel_start = 0.0, vel_stop = std::numeric_limits<double>::max();
            if ((*sub_model_part).Has(VELOCITY_START_TIME)) {
                vel_start = (*sub_model_part)[VELOCITY_START_TIME];
            }
            if ((*sub_model_part).Has(VELOCITY_STOP_TIME)) {
                vel_stop = (*sub_model_part)[VELOCITY_STOP_TIME];
            }

            if (time < vel_start || time > vel_stop) continue;

            NodesArrayType& pNodes = sub_model_part->Nodes();

            if ((*sub_model_part).Has(IMPOSED_VELOCITY_X_VALUE)) {
                SetFlagAndVariableToNodes(DEMFlags::FIXED_VEL_X, VELOCITY_X, (*sub_model_part)[IMPOSED_VELOCITY_X_VALUE], pNodes);
            }
            if ((*sub_model_part).Has(IMPOSED_VELOCITY_Y_VALUE)) {
                SetFlagAndVariableToNodes(DEMFlags::FIXED_VEL_Y, VELOCITY_Y, (*sub_model_part)[IMPOSED_VELOCITY_Y_VALUE], pNodes);
            }
            if ((*sub_model_part).Has(IMPOSED_VELOCITY_Z_VALUE)) {
                SetFlagAndVariableToNodes(DEMFlags::FIXED_VEL_Z, VELOCITY_Z, (*sub_model_part)[IMPOSED_VELOCITY_Z_VALUE], pNodes);
            }
            if ((*sub_model_part).Has(IMPOSED_ANGULAR_VELOCITY_X_VALUE)) {
                SetFlagAndVariableToNodes(DEMFlags::FIXED_ANG_VEL_X, ANGULAR_VELOCITY_X, (*sub_model_part)[IMPOSED_ANGULAR_VELOCITY_X_VALUE], pNodes);
            }
            if ((*sub_model_part).Has(IMPOSED_ANGULAR_VELOCITY_Y_VALUE)) {
                SetFlagAndVariableToNodes(DEMFlags::FIXED_ANG_VEL_Y, ANGULAR_VELOCITY_Y, (*sub_model_part)[IMPOSED_ANGULAR_VELOCITY_Y_VALUE], pNodes);
            }
            if ((*sub_model_part).Has(IMPOSED_ANGULAR_VELOCITY_Z_VALUE)) {
                SetFlagAndVariableToNodes(DEMFlags::FIXED_ANG_VEL_Z, ANGULAR_VELOCITY_Z, (*sub_model_part)[IMPOSED_ANGULAR_VELOCITY_Z_VALUE], pNodes);
            }
        } // for each mesh

        ModelPart& fem_model_part = GetFemModelPart();

        unsigned int rigid_body_elements_counter = 0;

        for (ModelPart::SubModelPartsContainerType::iterator sub_model_part = fem_model_part.SubModelPartsBegin(); sub_model_part != fem_model_part.SubModelPartsEnd(); ++sub_model_part) {

            ModelPart& submp = *sub_model_part;

            ElementsArrayType& pElements = mpFem_model_part->Elements();
            ElementsArrayType::iterator it = pElements.ptr_begin() + rigid_body_elements_counter;
            RigidBodyElement3D& rigid_body_element = dynamic_cast<Kratos::RigidBodyElement3D&> (*it);

            rigid_body_elements_counter++;

            if (submp.Has(FREE_BODY_MOTION)) { // JIG: Backward compatibility, it should be removed in the future
                if (submp[FREE_BODY_MOTION]) {
                    double vel_start = 0.0, vel_stop = std::numeric_limits<double>::max();
                    if (submp.Has(VELOCITY_START_TIME)) vel_start = submp[VELOCITY_START_TIME];
                    if (submp.Has(VELOCITY_STOP_TIME)) vel_stop = submp[VELOCITY_STOP_TIME];

                    if (time > vel_start && time < vel_stop) {

                        rigid_body_element.GetGeometry()[0].Set(DEMFlags::FIXED_VEL_X, false);
                        rigid_body_element.GetGeometry()[0].Set(DEMFlags::FIXED_VEL_Y, false);
                        rigid_body_element.GetGeometry()[0].Set(DEMFlags::FIXED_VEL_Z, false);
                        rigid_body_element.GetGeometry()[0].Set(DEMFlags::FIXED_ANG_VEL_X, false);
                        rigid_body_element.GetGeometry()[0].Set(DEMFlags::FIXED_ANG_VEL_Y, false);
                        rigid_body_element.GetGeometry()[0].Set(DEMFlags::FIXED_ANG_VEL_Z, false);

                        if (submp.Has(IMPOSED_VELOCITY_X_VALUE)) {
                            rigid_body_element.GetGeometry()[0].FastGetSolutionStepValue(VELOCITY)[0] = submp[IMPOSED_VELOCITY_X_VALUE];
                            rigid_body_element.GetGeometry()[0].Set(DEMFlags::FIXED_VEL_X, true);
                        }
                        if (submp.Has(IMPOSED_VELOCITY_Y_VALUE)) {
                            rigid_body_element.GetGeometry()[0].FastGetSolutionStepValue(VELOCITY)[1] = submp[IMPOSED_VELOCITY_Y_VALUE];
                            rigid_body_element.GetGeometry()[0].Set(DEMFlags::FIXED_VEL_Y, true);
                        }
                        if (submp.Has(IMPOSED_VELOCITY_Z_VALUE)) {
                            rigid_body_element.GetGeometry()[0].FastGetSolutionStepValue(VELOCITY)[2] = submp[IMPOSED_VELOCITY_Z_VALUE];
                            rigid_body_element.GetGeometry()[0].Set(DEMFlags::FIXED_VEL_Z, true);
                        }
                        if (submp.Has(IMPOSED_ANGULAR_VELOCITY_X_VALUE)) {
                            rigid_body_element.GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY)[0] = submp[IMPOSED_ANGULAR_VELOCITY_X_VALUE];
                            rigid_body_element.GetGeometry()[0].Set(DEMFlags::FIXED_ANG_VEL_X, true);
                        }
                        if (submp.Has(IMPOSED_ANGULAR_VELOCITY_Y_VALUE)) {
                            rigid_body_element.GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY)[1] = submp[IMPOSED_ANGULAR_VELOCITY_Y_VALUE];
                            rigid_body_element.GetGeometry()[0].Set(DEMFlags::FIXED_ANG_VEL_Y, true);
                        }
                        if (submp.Has(IMPOSED_ANGULAR_VELOCITY_Z_VALUE)) {
                            rigid_body_element.GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY)[2] = submp[IMPOSED_ANGULAR_VELOCITY_Z_VALUE];
                            rigid_body_element.GetGeometry()[0].Set(DEMFlags::FIXED_ANG_VEL_Z, true);
                        }

                        if (submp.Has(TABLE_NUMBER_VELOCITY)) { // JIG: Backward compatibility, it should be removed in the future
                            if (submp[TABLE_NUMBER_VELOCITY][0] != 0) {
                                const int table_number = submp[TABLE_NUMBER_VELOCITY_X];
                                rigid_body_element.GetGeometry()[0].FastGetSolutionStepValue(VELOCITY)[0] = submp.GetTable(table_number).GetValue(time);
                                rigid_body_element.GetGeometry()[0].Set(DEMFlags::FIXED_VEL_X, true);
                            }
                            if (submp[TABLE_NUMBER_VELOCITY][1] != 0) {
                                const int table_number = submp[TABLE_NUMBER_VELOCITY_Y];
                                rigid_body_element.GetGeometry()[0].FastGetSolutionStepValue(VELOCITY)[1] = submp.GetTable(table_number).GetValue(time);
                                rigid_body_element.GetGeometry()[0].Set(DEMFlags::FIXED_VEL_Y, true);
                            }
                            if (submp[TABLE_NUMBER_VELOCITY][2] != 0) {
                                const int table_number = submp[TABLE_NUMBER_VELOCITY_Z];
                                rigid_body_element.GetGeometry()[0].FastGetSolutionStepValue(VELOCITY)[2] = submp.GetTable(table_number).GetValue(time);
                                rigid_body_element.GetGeometry()[0].Set(DEMFlags::FIXED_VEL_Z, true);
                            }
                        }
                        if (submp.Has(TABLE_NUMBER_ANGULAR_VELOCITY)) { // JIG: Backward compatibility, it should be removed in the future
                            if (submp[TABLE_NUMBER_ANGULAR_VELOCITY][0] != 0) {
                                const int table_number = submp[TABLE_NUMBER_ANGULAR_VELOCITY_X];
                                rigid_body_element.GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY)[0] = submp.GetTable(table_number).GetValue(time);
                                rigid_body_element.GetGeometry()[0].Set(DEMFlags::FIXED_ANG_VEL_X, true);
                            }
                            if (submp[TABLE_NUMBER_ANGULAR_VELOCITY][1] != 0) {
                                const int table_number = submp[TABLE_NUMBER_ANGULAR_VELOCITY_Y];
                                rigid_body_element.GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY)[1] = submp.GetTable(table_number).GetValue(time);
                                rigid_body_element.GetGeometry()[0].Set(DEMFlags::FIXED_ANG_VEL_Y, true);
                            }
                            if (submp[TABLE_NUMBER_ANGULAR_VELOCITY][2] != 0) {
                                const int table_number = submp[TABLE_NUMBER_ANGULAR_VELOCITY_Z];
                                rigid_body_element.GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY)[2] = submp.GetTable(table_number).GetValue(time);
                                rigid_body_element.GetGeometry()[0].Set(DEMFlags::FIXED_ANG_VEL_Z, true);
                            }
                        }
                    }

                    else {
                        rigid_body_element.GetGeometry()[0].FastGetSolutionStepValue(VELOCITY)[0] = 0.0;
                        rigid_body_element.GetGeometry()[0].FastGetSolutionStepValue(VELOCITY)[1] = 0.0;
                        rigid_body_element.GetGeometry()[0].FastGetSolutionStepValue(VELOCITY)[2] = 0.0;
                        rigid_body_element.GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY)[0] = 0.0;
                        rigid_body_element.GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY)[1] = 0.0;
                        rigid_body_element.GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY)[2] = 0.0;
                        rigid_body_element.GetGeometry()[0].Set(DEMFlags::FIXED_VEL_X, true);
                        rigid_body_element.GetGeometry()[0].Set(DEMFlags::FIXED_VEL_Y, true);
                        rigid_body_element.GetGeometry()[0].Set(DEMFlags::FIXED_VEL_Z, true);
                        rigid_body_element.GetGeometry()[0].Set(DEMFlags::FIXED_ANG_VEL_X, true);
                        rigid_body_element.GetGeometry()[0].Set(DEMFlags::FIXED_ANG_VEL_Y, true);
                        rigid_body_element.GetGeometry()[0].Set(DEMFlags::FIXED_ANG_VEL_Z, true);
                    }

                    if (submp.Has(TABLE_NUMBER_FORCE)) { // JIG: Backward compatibility, it should be removed in the future
                        if (submp[TABLE_NUMBER_FORCE][0] != 0) {
                            const int table_number = submp[TABLE_NUMBER_FORCE][0];
                            rigid_body_element.GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_APPLIED_FORCE)[0] = submp.GetTable(table_number).GetValue(time);
                        }
                        if (submp[TABLE_NUMBER_FORCE][1] != 0) {
                            const int table_number = submp[TABLE_NUMBER_FORCE][1];
                            rigid_body_element.GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_APPLIED_FORCE)[1] = submp.GetTable(table_number).GetValue(time);
                        }
                        if (submp[TABLE_NUMBER_FORCE][2] != 0) {
                            const int table_number = submp[TABLE_NUMBER_FORCE][2];
                            rigid_body_element.GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_APPLIED_FORCE)[2] = submp.GetTable(table_number).GetValue(time);
                        }
                    }

                    if (submp.Has(TABLE_NUMBER_MOMENT)) { // JIG: Backward compatibility, it should be removed in the future
                        if (submp[TABLE_NUMBER_MOMENT][0] != 0) {
                            const int table_number = submp[TABLE_NUMBER_MOMENT][0];
                            rigid_body_element.GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_APPLIED_MOMENT)[0] = submp.GetTable(table_number).GetValue(time);
                        }
                        if (submp[TABLE_NUMBER_MOMENT][1] != 0) {
                            const int table_number = submp[TABLE_NUMBER_MOMENT][1];
                            rigid_body_element.GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_APPLIED_MOMENT)[1] = submp.GetTable(table_number).GetValue(time);
                        }
                        if (submp[TABLE_NUMBER_MOMENT][2] != 0) {
                            const int table_number = submp[TABLE_NUMBER_MOMENT][2];
                            rigid_body_element.GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_APPLIED_MOMENT)[2] = submp.GetTable(table_number).GetValue(time);
                        }
                    }
                }
            }
        }
        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::ApplyInitialConditions() {

        KRATOS_TRY
        ModelPart& r_model_part = GetModelPart();

        for (ModelPart::SubModelPartsContainerType::iterator sub_model_part = r_model_part.SubModelPartsBegin(); sub_model_part != r_model_part.SubModelPartsEnd(); ++sub_model_part) {

            NodesArrayType& pNodes = sub_model_part->Nodes();

            if ((*sub_model_part).Has(INITIAL_VELOCITY_X_VALUE)) {
                SetVariableToNodes(VELOCITY_X, (*sub_model_part)[INITIAL_VELOCITY_X_VALUE], pNodes);
            }
            if ((*sub_model_part).Has(INITIAL_VELOCITY_Y_VALUE)) {
                SetVariableToNodes(VELOCITY_Y, (*sub_model_part)[INITIAL_VELOCITY_Y_VALUE], pNodes);
            }
            if ((*sub_model_part).Has(INITIAL_VELOCITY_Z_VALUE)) {
                SetVariableToNodes(VELOCITY_Z, (*sub_model_part)[INITIAL_VELOCITY_Z_VALUE], pNodes);
            }
            if ((*sub_model_part).Has(INITIAL_ANGULAR_VELOCITY_X_VALUE)) {
                SetVariableToNodes(ANGULAR_VELOCITY_X, (*sub_model_part)[INITIAL_ANGULAR_VELOCITY_X_VALUE], pNodes);
            }
            if ((*sub_model_part).Has(INITIAL_ANGULAR_VELOCITY_Y_VALUE)) {
                SetVariableToNodes(ANGULAR_VELOCITY_Y, (*sub_model_part)[INITIAL_ANGULAR_VELOCITY_Y_VALUE], pNodes);
            }
            if ((*sub_model_part).Has(INITIAL_ANGULAR_VELOCITY_Z_VALUE)) {
                SetVariableToNodes(ANGULAR_VELOCITY_Z, (*sub_model_part)[INITIAL_ANGULAR_VELOCITY_Z_VALUE], pNodes);
            }
        } // for each mesh

        ModelPart& fem_model_part = GetFemModelPart();

        unsigned int rigid_body_elements_counter = 0;

        for (ModelPart::SubModelPartsContainerType::iterator sub_model_part = fem_model_part.SubModelPartsBegin(); sub_model_part != fem_model_part.SubModelPartsEnd(); ++sub_model_part) {

            ElementsArrayType& pElements = mpFem_model_part->Elements();
            ElementsArrayType::iterator it = pElements.ptr_begin() + rigid_body_elements_counter;
            RigidBodyElement3D& rigid_body_element = dynamic_cast<Kratos::RigidBodyElement3D&> (*it);

            if ((*sub_model_part).Has(INITIAL_VELOCITY_X_VALUE)) {
                rigid_body_element.GetGeometry()[0].FastGetSolutionStepValue(VELOCITY)[0] = (*sub_model_part)[INITIAL_VELOCITY_X_VALUE];
            }
            if ((*sub_model_part).Has(INITIAL_VELOCITY_Y_VALUE)) {
                rigid_body_element.GetGeometry()[0].FastGetSolutionStepValue(VELOCITY)[1] = (*sub_model_part)[INITIAL_VELOCITY_Y_VALUE];
            }
            if ((*sub_model_part).Has(INITIAL_VELOCITY_Z_VALUE)) {
                rigid_body_element.GetGeometry()[0].FastGetSolutionStepValue(VELOCITY)[2] = (*sub_model_part)[INITIAL_VELOCITY_Z_VALUE];
            }
            if ((*sub_model_part).Has(INITIAL_ANGULAR_VELOCITY_X_VALUE)) {
                rigid_body_element.GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY)[0] = (*sub_model_part)[INITIAL_ANGULAR_VELOCITY_X_VALUE];
            }
            if ((*sub_model_part).Has(INITIAL_ANGULAR_VELOCITY_Y_VALUE)) {
                rigid_body_element.GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY)[1] = (*sub_model_part)[INITIAL_ANGULAR_VELOCITY_Y_VALUE];
            }
            if ((*sub_model_part).Has(INITIAL_ANGULAR_VELOCITY_Z_VALUE)) {
                rigid_body_element.GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY)[2] = (*sub_model_part)[INITIAL_ANGULAR_VELOCITY_Z_VALUE];
            }

            rigid_body_elements_counter++;
        }

        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::SetSearchRadiiOnAllParticles(ModelPart& r_model_part, const double added_search_distance, const double amplification) {
        KRATOS_TRY
        int number_of_elements = r_model_part.GetCommunicator().LocalMesh().ElementsArray().end() - r_model_part.GetCommunicator().LocalMesh().ElementsArray().begin();
        #pragma omp parallel for
        for (int i = 0; i < number_of_elements; i++) {
            mListOfSphericParticles[i]->SetSearchRadius(amplification * (added_search_distance + mListOfSphericParticles[i]->GetRadius()));
        }
        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::SetNormalRadiiOnAllParticles(ModelPart& r_model_part) {
        KRATOS_TRY
        int number_of_elements = r_model_part.GetCommunicator().LocalMesh().ElementsArray().end() - r_model_part.GetCommunicator().LocalMesh().ElementsArray().begin();
        #pragma omp parallel for
        for (int i = 0; i < number_of_elements; i++) {
            mListOfSphericParticles[i]->SetRadius();
        }
        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::SetSearchRadiiWithFemOnAllParticles(ModelPart& r_model_part, const double added_search_distance, const double amplification) {
        KRATOS_TRY
        int number_of_elements = r_model_part.GetCommunicator().LocalMesh().ElementsArray().end() - r_model_part.GetCommunicator().LocalMesh().ElementsArray().begin();
        #pragma omp parallel for
        for (int i = 0; i < number_of_elements; i++) {
            mListOfSphericParticles[i]->SetSearchRadius(amplification * (added_search_distance + mListOfSphericParticles[i]->GetRadius()));
        }
        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::SearchNeighbours() {

        KRATOS_TRY

        if (!mDoSearchNeighbourElements) {
            return;
        }

        ModelPart& r_model_part = GetModelPart();

        int number_of_elements = r_model_part.GetCommunicator().LocalMesh().ElementsArray().end() - r_model_part.GetCommunicator().LocalMesh().ElementsArray().begin();
        if (!number_of_elements) return;

        GetResults().resize(number_of_elements);
        GetResultsDistances().resize(number_of_elements);

        //SetSearchRadiiOnAllParticles(r_model_part, r_model_part.GetProcessInfo()[SEARCH_RADIUS_INCREMENT], 1.0);
        mpSpSearch->SearchElementsInRadiusExclusive(r_model_part, this->GetArrayOfAmplifiedRadii(), this->GetResults(), this->GetResultsDistances());

        const int number_of_particles = (int) mListOfSphericParticles.size();

        #pragma omp parallel for schedule(dynamic, 100)
        for (int i = 0; i < number_of_particles; i++) {
            mListOfSphericParticles[i]->mNeighbourElements.clear();
            for (SpatialSearch::ResultElementsContainerType::iterator neighbour_it = this->GetResults()[i].begin(); neighbour_it != this->GetResults()[i].end(); ++neighbour_it) {
                Element* p_neighbour_element = (*neighbour_it).get();
                SphericParticle* p_spheric_neighbour_particle = dynamic_cast<SphericParticle*> (p_neighbour_element);
                if (mListOfSphericParticles[i]->Is(DEMFlags::BELONGS_TO_A_CLUSTER) && (mListOfSphericParticles[i]->GetClusterId() == p_spheric_neighbour_particle->GetClusterId())) continue;
                mListOfSphericParticles[i]->mNeighbourElements.push_back(p_spheric_neighbour_particle);
            }
            this->GetResults()[i].clear();
            this->GetResultsDistances()[i].clear();
        }
        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::ComputeNewNeighboursHistoricalData() {

        KRATOS_TRY
        const int number_of_particles = (int) mListOfSphericParticles.size();

        #pragma omp parallel
        {
            DenseVector<int> mTempNeighboursIds;
            std::vector<array_1d<double, 3> > mTempNeighbourElasticContactForces;
            std::vector<array_1d<double, 3> > mTempNeighbourTotalContactForces;

            #pragma omp for
            for (int i = 0; i < number_of_particles; i++) {
                mListOfSphericParticles[i]->ComputeNewNeighboursHistoricalData(mTempNeighboursIds, mTempNeighbourElasticContactForces);
            }
        }

        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::ComputeNewRigidFaceNeighboursHistoricalData() {
        KRATOS_TRY
        const int number_of_particles = (int) mListOfSphericParticles.size();

        #pragma omp parallel for
        for (int i = 0; i < number_of_particles; i++) {
            mListOfSphericParticles[i]->ComputeNewRigidFaceNeighboursHistoricalData();
        }

        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::SearchRigidFaceNeighbours() {
        KRATOS_TRY

        ElementsArrayType& pElements = mpDem_model_part->GetCommunicator().LocalMesh().Elements();
        ConditionsArrayType& pTConditions = mpFem_model_part->GetCommunicator().LocalMesh().Conditions();

        if (pTConditions.size() > 0) {
            const int number_of_particles = (int) mListOfSphericParticles.size();

            this->GetRigidFaceResults().resize(number_of_particles);
            this->GetRigidFaceResultsDistances().resize(number_of_particles);

            //Fast Bins Search
            //SetSearchRadiiOnAllParticles(*mpDem_model_part, mpDem_model_part->GetProcessInfo()[SEARCH_RADIUS_INCREMENT_FOR_WALLS], 1.0);
            mpDemFemSearch->SearchRigidFaceForDEMInRadiusExclusiveImplementation(pElements, pTConditions, this->GetRigidFaceResults(), this->GetRigidFaceResultsDistances());


            #pragma omp parallel for schedule(dynamic, 100)
            for (int i = 0; i < number_of_particles; i++) {
                mListOfSphericParticles[i]->mNeighbourPotentialRigidFaces.clear();
                for (ResultConditionsContainerType::iterator neighbour_it = this->GetRigidFaceResults()[i].begin(); neighbour_it != this->GetRigidFaceResults()[i].end(); ++neighbour_it) {
                    Condition* p_neighbour_condition = (*neighbour_it).get();
                    DEMWall* p_wall = dynamic_cast<DEMWall*> (p_neighbour_condition);
                    mListOfSphericParticles[i]->mNeighbourPotentialRigidFaces.push_back(p_wall);
                }//for results iterator
                this->GetRigidFaceResults()[i].clear();
                this->GetRigidFaceResultsDistances()[i].clear();
            }

            CheckHierarchyWithCurrentNeighbours();
            //DoubleHierarchyMethod();

            //typedef WeakPointerVector<Condition >::iterator ConditionWeakIteratorType;
            const int number_of_conditions = (int) pTConditions.size();

            #pragma omp parallel
            {
                #pragma omp for
                for (int i = 0; i < number_of_conditions; i++) {
                    ConditionsArrayType::iterator ic = pTConditions.begin() + i;
                    DEMWall* wall = dynamic_cast<Kratos::DEMWall*> (&(*ic));
                    wall->mNeighbourSphericParticles.resize(0);
                }

                #pragma omp for
                for (int i = 0; i < number_of_particles; i++) {
                    for (unsigned int j = 0; j < mListOfSphericParticles[i]->mNeighbourRigidFaces.size(); j++) {
                        DEMWall* p_wall = mListOfSphericParticles[i]->mNeighbourRigidFaces[j];
                        #pragma omp critical
                        {
                            p_wall->mNeighbourSphericParticles.push_back(mListOfSphericParticles[i]);
                        }
                    }
                }
            }//#pragma omp parallel
        }
        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::DoubleHierarchyMethod() {
        KRATOS_TRY
        const int number_of_particles = (int) mListOfSphericParticles.size();

        #pragma omp parallel
        {
        std::vector< double > Distance_Array; //MACELI: reserve.. or take it out of the loop and have one for every thread
        std::vector< array_1d<double, 3> > Normal_Array;
        std::vector< array_1d<double, 4> > Weight_Array;
        std::vector< int > Id_Array;
        std::vector< int > ContactType_Array;
        std::vector<DEMWall*> temporal_neigh;
        std::vector< array_1d<double, 4> > temporal_contact_weights;
        std::vector< int > temporal_contact_types;

        #pragma omp for
        for (int i = 0; i < number_of_particles; i++) {
            SphericParticle* p_sphere_i = mListOfSphericParticles[i];
            p_sphere_i->mNeighbourRigidFaces.resize(0);
            p_sphere_i->mNeighbourPotentialRigidFaces.resize(0);
            p_sphere_i->mContactConditionWeights.resize(0);

            Distance_Array.clear();
            Normal_Array.clear();
            Weight_Array.clear();
            Id_Array.clear();
            ContactType_Array.clear();

            std::vector<DEMWall*>& potential_neighbour_rigid_faces = p_sphere_i->mNeighbourPotentialRigidFaces;

            for (ResultConditionsContainerType::iterator neighbour_it = this->GetRigidFaceResults()[i].begin(); neighbour_it != this->GetRigidFaceResults()[i].end(); ++neighbour_it) {

                Condition* p_neighbour_condition = (*neighbour_it).get();
                DEMWall* p_wall = dynamic_cast<DEMWall*> (p_neighbour_condition);
                RigidFaceGeometricalConfigureType::DoubleHierarchyMethod(p_sphere_i,
                        p_wall,
                        Distance_Array,
                        Normal_Array,
                        Weight_Array,
                        Id_Array,
                        ContactType_Array
                        );
                potential_neighbour_rigid_faces.push_back(p_wall);

            }//for results iterator

            std::vector<DEMWall*>& neighbour_rigid_faces = p_sphere_i->mNeighbourRigidFaces;
            std::vector< array_1d<double, 4> >& neighbour_weights = p_sphere_i->mContactConditionWeights;
            std::vector< int >& neighbor_contact_types = p_sphere_i->mContactConditionContactTypes;

            size_t neigh_size = neighbour_rigid_faces.size();

            temporal_neigh.clear();
            temporal_contact_weights.clear();
            temporal_contact_types.clear();

            for (unsigned int n = 0; n < neigh_size; n++) {

                if (ContactType_Array[n] != -1) //if(it is not a -1 contact neighbour, we copy it)
                {
                    temporal_neigh.push_back(neighbour_rigid_faces[n]);
                    temporal_contact_weights.push_back(Weight_Array[n]);
                    temporal_contact_types.push_back(ContactType_Array[n]);

                }//if(it is not a -1 contact neighbour, we copy it)

            }//loop over temporal neighbours

            //swap

            temporal_neigh.swap(neighbour_rigid_faces);
            temporal_contact_weights.swap(neighbour_weights);
            temporal_contact_types.swap(neighbor_contact_types);

            this->GetRigidFaceResults()[i].clear();
            this->GetRigidFaceResultsDistances()[i].clear();
        }//for particles
        } //parallel region
        KRATOS_CATCH("")
    }//DoubleHierarchyMethod

    void ExplicitSolverStrategy::CheckHierarchyWithCurrentNeighbours()
        {
        KRATOS_TRY
        const int number_of_particles = (int) mListOfSphericParticles.size();

        #pragma omp parallel
        {
            std::vector< double > Distance_Array; //MACELI: reserve.. or take it out of the loop and have one for every thread
            std::vector< array_1d<double, 3> > Normal_Array;
            std::vector< array_1d<double, 4> > Weight_Array;
            std::vector< int > Id_Array;
            std::vector< int > ContactType_Array;

            #pragma omp for schedule(dynamic, 100)
            for (int i = 0; i < number_of_particles; i++) {
                SphericParticle* p_sphere_i = mListOfSphericParticles[i];
                p_sphere_i->mNeighbourRigidFaces.resize(0);
                p_sphere_i->mContactConditionWeights.resize(0);

                Distance_Array.clear();
                Normal_Array.clear();
                Weight_Array.clear();
                Id_Array.clear();
                ContactType_Array.clear();
                std::vector<DEMWall*>& potential_neighbour_rigid_faces = p_sphere_i->mNeighbourPotentialRigidFaces;

                for (unsigned int n = 0; n < potential_neighbour_rigid_faces.size(); ++n) {
                    Condition* p_neighbour_condition = potential_neighbour_rigid_faces[n];
                    DEMWall* p_wall = dynamic_cast<DEMWall*> (p_neighbour_condition);
                    RigidFaceGeometricalConfigureType::DoubleHierarchyMethod(p_sphere_i,
                            p_wall,
                            Distance_Array,
                            Normal_Array,
                            Weight_Array,
                            Id_Array,
                            ContactType_Array
                            );

                }//loop over temporal neighbours

                std::vector<DEMWall*>& neighbour_rigid_faces = p_sphere_i->mNeighbourRigidFaces;
                std::vector< array_1d<double, 4> >& neighbour_weights = p_sphere_i->mContactConditionWeights;
                std::vector< int >& neighbor_contact_types = p_sphere_i->mContactConditionContactTypes;

                size_t neigh_size = neighbour_rigid_faces.size();

                std::vector<DEMWall*> temporal_neigh(0);
                std::vector< array_1d<double, 4> > temporal_contact_weights;
                std::vector< int > temporal_contact_types;

                for (unsigned int n = 0; n < neigh_size; n++) {

                    if (ContactType_Array[n] != -1) //if(it is not a -1 contact neighbour, we copy it)
                    {
                        temporal_neigh.push_back(neighbour_rigid_faces[n]);
                        temporal_contact_weights.push_back(Weight_Array[n]);
                        temporal_contact_types.push_back(ContactType_Array[n]);

                    }//if(it is not a -1 contact neighbour, we copy it)

                }//loop over temporal neighbours

                //swap

                temporal_neigh.swap(neighbour_rigid_faces);
                temporal_contact_weights.swap(neighbour_weights);
                temporal_contact_types.swap(neighbor_contact_types);


            }//for particles
        }

        KRATOS_CATCH("")
        }//CheckHierarchyWithCurrentNeighbours

    void ExplicitSolverStrategy::CalculateInitialMaxIndentations(ProcessInfo& r_process_info) {
        KRATOS_TRY
        std::vector<double> indentations_list, indentations_list_ghost;
        indentations_list.resize(mListOfSphericParticles.size());
        indentations_list_ghost.resize(mListOfGhostSphericParticles.size());

        const int number_of_particles = (int) mListOfSphericParticles.size();

        #pragma omp parallel
        {
            #pragma omp for
            for (int i = 0; i < number_of_particles; i++) {
                double indentation;
                mListOfSphericParticles[i]->CalculateMaxBallToBallIndentation(indentation, r_process_info);
                double max_indentation = std::max(0.0, 0.5 * indentation); // reducing the radius by half the indentation is enough

                mListOfSphericParticles[i]->CalculateMaxBallToFaceIndentation(indentation);
                max_indentation = std::max(max_indentation, indentation);
                indentations_list[i] = max_indentation;
            }

            #pragma omp for //THESE TWO LOOPS CANNOT BE JOINED, BECAUSE THE RADII ARE CHANGING.
            for (int i = 0; i < number_of_particles; i++) {
                mListOfSphericParticles[i]->SetInteractionRadius(mListOfSphericParticles[i]->GetInteractionRadius() - indentations_list[i]);
            }

            #pragma omp single
            {
                SynchronizeHistoricalVariables(GetModelPart());
            }
            const int number_of_ghost_particles = (int) mListOfGhostSphericParticles.size();

            #pragma omp for //THESE TWO LOOPS CANNOT BE JOINED, BECAUSE THE RADII ARE CHANGING.
            for (int i = 0; i < number_of_ghost_particles; i++) {
                mListOfGhostSphericParticles[i]->SetInteractionRadius(mListOfGhostSphericParticles[i]->GetInteractionRadius() - indentations_list_ghost[i]);
            }

            #pragma omp for
            for (int i = 0; i < number_of_particles; i++) {
                double indentation;
                mListOfSphericParticles[i]->CalculateMaxBallToBallIndentation(indentation, r_process_info);
            }
        } //#pragma omp parallel

        KRATOS_CATCH("")
    } // CalculateInitialMaxIndentations()

    void ExplicitSolverStrategy::PrepareContactModelPart(ModelPart& r_model_part, ModelPart& mcontacts_model_part) {
        mcontacts_model_part.GetCommunicator().SetNumberOfColors(r_model_part.GetCommunicator().GetNumberOfColors());
        mcontacts_model_part.GetCommunicator().NeighbourIndices() = r_model_part.GetCommunicator().NeighbourIndices();
    }

    void ExplicitSolverStrategy::PrepareElementsForPrinting() {
        KRATOS_TRY
        ProcessInfo& r_process_info = (*mpDem_model_part).GetProcessInfo();
        ElementsArrayType& pElements = (*mpDem_model_part).GetCommunicator().LocalMesh().Elements();

        DenseVector<unsigned int> element_partition;

        OpenMPUtils::CreatePartition(mNumberOfThreads, pElements.size(), element_partition);

        #pragma omp parallel for
        for (int k = 0; k < (int) pElements.size(); k++) {
            ElementsArrayType::iterator it = pElements.ptr_begin() + k;
            Element* raw_p_element = &(*it);
            SphericParticle* p_sphere = dynamic_cast<SphericParticle*> (raw_p_element);
            p_sphere->PrepareForPrinting(r_process_info);
        }
        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::SynchronizeHistoricalVariables(ModelPart& r_model_part) {
        r_model_part.GetCommunicator().SynchronizeNodalSolutionStepsData();
    }

    void ExplicitSolverStrategy::SynchronizeRHS(ModelPart& r_model_part) {
        r_model_part.GetCommunicator().SynchronizeVariable(TOTAL_FORCES);
        r_model_part.GetCommunicator().SynchronizeVariable(PARTICLE_MOMENT);
    }

    void ExplicitSolverStrategy::CleanEnergies() {
        KRATOS_TRY

        ProcessInfo& r_process_info = GetModelPart().GetProcessInfo();
        double& total_elastic_energy = r_process_info[PARTICLE_ELASTIC_ENERGY];
        total_elastic_energy = 0.0;
        double& total_inelastic_frictional_energy = r_process_info[PARTICLE_INELASTIC_FRICTIONAL_ENERGY];
        total_inelastic_frictional_energy  = 0.0;
        double& total_inelastic_viscodamping_energy = r_process_info[PARTICLE_INELASTIC_VISCODAMPING_ENERGY];
        total_inelastic_viscodamping_energy  = 0.0;

        KRATOS_CATCH("")
    }
}
