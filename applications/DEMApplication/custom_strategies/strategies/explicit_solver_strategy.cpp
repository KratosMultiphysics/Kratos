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


        std::vector<PropertiesProxy>& vector_of_properties_proxies = PropertiesProxiesManager().GetPropertiesProxies(*mpDem_model_part);

        IndexPartition<unsigned int>(rCustomListOfSphericParticles.size()).for_each([&](unsigned int i){
            rCustomListOfSphericParticles[i]->SetFastProperties(vector_of_properties_proxies);
        });

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
        // Using IndexPartition should be fine since 'break' affects the internal for loops while the replaced continues only has an effect on the for_each loop.
        IndexPartition<unsigned int>(rCustomListOfSphericParticles.size()).for_each([&](unsigned int i){

            int own_properties_id = rCustomListOfSphericParticles[i]->GetProperties().Id();
            for (PropertiesIterator props_it = mpDem_model_part->GetMesh(0).PropertiesBegin(); props_it != mpDem_model_part->GetMesh(0).PropertiesEnd(); props_it++) {
                int model_part_id = props_it->GetId();
                if (own_properties_id == model_part_id) {
                    rCustomListOfSphericParticles[i]->SetProperties(*(props_it.base()));
                    found = true;
                    break;
                }
            }
            if (found) return;

            for (PropertiesIterator props_it = mpInlet_model_part->GetMesh(0).PropertiesBegin(); props_it != mpInlet_model_part->GetMesh(0).PropertiesEnd(); props_it++) {
                int model_part_id = props_it->GetId();
                if (own_properties_id == model_part_id) {
                    rCustomListOfSphericParticles[i]->SetProperties(*(props_it.base()));
                    found = true;
                    break;
                }
            }
            if (found) return;

            for (PropertiesIterator props_it = mpCluster_model_part->GetMesh(0).PropertiesBegin(); props_it != mpCluster_model_part->GetMesh(0).PropertiesEnd(); props_it++) {
                int model_part_id = props_it->GetId();
                if (own_properties_id == model_part_id) {
                    rCustomListOfSphericParticles[i]->SetProperties(*(props_it.base()));
                    found = true;
                    break;
                }
            }

            KRATOS_ERROR_IF_NOT(found) << "This particle could not find its properties!!" << std::endl;
        });

        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::DisplayThreadInfo() {
        KRATOS_TRY
        ModelPart& r_model_part = GetModelPart();
        KRATOS_INFO("DEM") << "          **************************************************" << std::endl;
        KRATOS_INFO("DEM") << "            Parallelism Info:  MPI number of nodes: " << r_model_part.GetCommunicator().TotalProcesses() << std::endl;
        if (r_model_part.GetCommunicator().TotalProcesses() > 1)
            KRATOS_INFO("DEM") << "            Parallelism Info:  MPI node Id: " << r_model_part.GetCommunicator().MyPID() << std::endl;
        KRATOS_INFO("DEM") << "            Parallelism Info:  OMP number of processors: " << mNumberOfThreads << std::endl;
        KRATOS_INFO("DEM") << "          **************************************************" << std::endl;
        KRATOS_INFO("DEM") << std::endl;
        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::Initialize() {

        KRATOS_TRY

        ModelPart& r_model_part = GetModelPart();

        ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
        SendProcessInfoToClustersModelPart();

        if (r_model_part.GetCommunicator().MyPID() == 0) {
            KRATOS_INFO("DEM") << "------------------DISCONTINUUM SOLVER STRATEGY---------------------" << "\n" << std::endl;
        }

        mNumberOfThreads = ParallelUtilities::GetNumThreads();
        DisplayThreadInfo();

        RebuildListOfSphericParticles<SphericParticle>(r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericParticles);
        RebuildListOfSphericParticles<SphericParticle>(r_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostSphericParticles);

        PropertiesProxiesManager().CreatePropertiesProxies(*mpDem_model_part, *mpInlet_model_part, *mpCluster_model_part);

        bool has_mpi = false;
        Check_MPI(has_mpi);

        if (has_mpi) {
            RepairPointersToNormalProperties(mListOfSphericParticles); // The particles sent to this partition have their own copy of the Kratos properties they were using in the previous partition!!
            RepairPointersToNormalProperties(mListOfGhostSphericParticles);
        }

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
            mpParticleCreatorDestructor->DestroyParticles<SphericParticle>(r_model_part);
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

        //set flag to 2 (search performed this timestep)
        mSearchControl = 2;

        // Finding overlapping of initial configurations
        if (r_process_info[CLEAN_INDENT_OPTION]) {
            for (int i = 0; i < 10; i++) CalculateInitialMaxIndentations(r_process_info);
        }

        //FinalizeSolutionStep();

        ComputeNodalArea();

        RVEInitialize();

        KRATOS_CATCH("")
    } // Initialize()

    void ExplicitSolverStrategy::AttachSpheresToStickyWalls() {
        KRATOS_TRY
        for (ModelPart::SubModelPartsContainerType::iterator sub_model_part = GetFemModelPart().SubModelPartsBegin(); sub_model_part != GetFemModelPart().SubModelPartsEnd(); ++sub_model_part) {

            ModelPart& submp = *sub_model_part;
            if(!submp[IS_STICKY]) continue;

            ConditionsArrayType& rConditions = submp.GetCommunicator().LocalMesh().Conditions();

            block_for_each(rConditions, [&](ModelPart::ConditionType& rCondition){
                rCondition.Set(DEMFlags::STICKY, true);
            });
        }

        const int number_of_particles = (int) mListOfSphericParticles.size();

        #pragma omp parallel for schedule(dynamic, 100)
        for (int i = 0; i < number_of_particles; i++) {
            std::vector<DEMWall*>& neighbour_walls_vector = mListOfSphericParticles[i]->mNeighbourPotentialRigidFaces;
            for (int j = 0; j<(int)neighbour_walls_vector.size(); j++) {
                if( neighbour_walls_vector[j]->Is(DEMFlags::STICKY) ) {
                    const bool is_inside = mListOfSphericParticles[i]->SwapIntegrationSchemeToGluedToWall(neighbour_walls_vector[j]);
                    if(is_inside) {
                        #pragma omp critical
                        {
                            neighbour_walls_vector[j]->GetVectorOfGluedParticles().push_back(mListOfSphericParticles[i]);
                        }
                        mListOfSphericParticles[i]->Set(DEMFlags::STICKY, true);
                        break;
                    }
                }
            }
        }
        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::MarkToDeleteAllSpheresInitiallyIndentedWithFEM(ModelPart& rSpheresModelPart) {
        KRATOS_TRY
        ElementsArrayType& rElements = rSpheresModelPart.GetCommunicator().LocalMesh().Elements();

        block_for_each(rElements, [&](ModelPart::ElementType& rElement) {
            Element* p_element = &(rElement);
            SphericParticle* p_sphere = dynamic_cast<SphericParticle*>(p_element);

            if (p_sphere->mNeighbourRigidFaces.size()) {
                p_sphere->Set(TO_ERASE);
                p_sphere->GetGeometry()[0].Set(TO_ERASE);
            }
        });

        KRATOS_CATCH("")
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
            const double inv_geometry_size = 1.0 / geometry.size();
            for (unsigned int i = 0; i < geometry.size(); i++) {
                double& node_area = geometry[i].FastGetSolutionStepValue(DEM_NODAL_AREA);
                node_area += inv_geometry_size * Element_Area;
            }

        }

        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::Check_MPI(bool& has_mpi) {
        VariablesList r_modelpart_nodal_variables_list = GetModelPart().GetNodalSolutionStepVariablesList();
        if (r_modelpart_nodal_variables_list.Has(PARTITION_INDEX)) has_mpi = true;
    }

    double ExplicitSolverStrategy::CalculateMaxInletTimeStep() {
        for (PropertiesIterator props_it = mpInlet_model_part->GetMesh(0).PropertiesBegin(); props_it != mpInlet_model_part->GetMesh(0).PropertiesEnd(); props_it++) {
            if ((*props_it).Has(PARTICLE_DENSITY)) {
                int inlet_prop_id = props_it->GetId();
                double young = (*props_it)[YOUNG_MODULUS];
                double density = (*props_it)[PARTICLE_DENSITY];
                double poisson = (*props_it)[POISSON_RATIO];

                for (ModelPart::SubModelPartsContainerType::iterator sub_model_part = mpInlet_model_part->SubModelPartsBegin(); sub_model_part != mpInlet_model_part->SubModelPartsEnd(); ++sub_model_part) {
                    KRATOS_ERROR_IF(!(*sub_model_part).Has(PROPERTIES_ID))<<"PROPERTIES_ID is not set for SubModelPart "<<(*sub_model_part).Name()<<" . Make sure the Materials file contains material assignation for this SubModelPart"<<std::endl;
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

    double ExplicitSolverStrategy::SolveSolutionStep() {
        KRATOS_TRY
        ModelPart& r_model_part = GetModelPart();

        SearchDEMOperations(r_model_part);
        SearchFEMOperations(r_model_part);
        ForceOperations(r_model_part);
        PerformTimeIntegrationOfMotion();

        return 0.00;

        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::SearchDEMOperations(ModelPart& r_model_part, bool has_mpi) {
        KRATOS_TRY
        ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
        int time_step = r_process_info[TIME_STEPS];
        const double time = r_process_info[TIME];
        const bool is_time_to_search_neighbours = (time_step + 1) % mNStepSearch == 0 && (time_step > 0); //Neighboring search. Every N times.
        const bool is_time_to_print_results = r_process_info[IS_TIME_TO_PRINT];
        const bool is_time_to_mark_and_remove = is_time_to_search_neighbours && (r_process_info[BOUNDING_BOX_OPTION] && time >= r_process_info[BOUNDING_BOX_START_TIME] && time <= r_process_info[BOUNDING_BOX_STOP_TIME]);
        BoundingBoxUtility(is_time_to_mark_and_remove);

        if (is_time_to_search_neighbours) {
            if (!is_time_to_mark_and_remove) { //Just in case that some entities were marked as TO_ERASE without a bounding box (manual removal)
                mpParticleCreatorDestructor->DestroyParticles<Cluster3D>(*mpCluster_model_part);
                mpParticleCreatorDestructor->DestroyParticles<SphericParticle>(r_model_part);
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

        if (is_time_to_print_results && r_process_info[CONTACT_MESH_OPTION] == 1) {
            CreateContactElements();
            InitializeContactElements();
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
            SetSearchRadiiWithFemOnAllParticles(r_model_part, mpDem_model_part->GetProcessInfo()[SEARCH_RADIUS_INCREMENT_FOR_WALLS], 1.0);
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

        GetForce(); // Basically only calls CalculateRightHandSide()
        //FastGetForce();
        GetClustersForce();
        GetRigidBodyElementsForce();

        if (r_model_part.GetProcessInfo()[COMPUTE_FEM_RESULTS_OPTION]) {
            CalculateNodalPressuresAndStressesOnWalls();
        }

        // Synchronize (should be just FORCE and TORQUE)
        SynchronizeRHS(r_model_part);

        KRATOS_CATCH("")
    }//ForceOperations;

    void ExplicitSolverStrategy::GetForce() {

        KRATOS_TRY

        ProcessInfo& r_process_info = GetModelPart().GetProcessInfo();
        double dt = r_process_info[DELTA_TIME];
        const array_1d<double, 3>& gravity = r_process_info[GRAVITY];

        const int number_of_particles = (int) mListOfSphericParticles.size();

        #pragma omp parallel for schedule(dynamic, 100)
        for (int i = 0; i < number_of_particles; i++) {
            RVEExecuteParticlePre(mListOfSphericParticles[i]);
            mListOfSphericParticles[i]->CalculateRightHandSide(r_process_info, dt, gravity);
            #pragma omp critical
            {
              RVEExecuteParticlePos(mListOfSphericParticles[i]);
            }
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
                mListOfSphericParticles[i]->FirstCalculateRightHandSide(r_process_info, dt);
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
            KRATOS_ERROR_IF((force_reduction_factor > 1.0) || (force_reduction_factor < 0.0)) << "The force reduction factor is either larger than 1 or negative: FORCE_REDUCTION_FACTOR= "<< virtual_mass_coeff << std::endl;
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
              if (mListOfSphericParticles[i]->mMoving)
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
        const ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
        ElementsArrayType& pElements = r_model_part.GetCommunicator().LocalMesh().Elements();

        ModelPart& r_fem_model_part = GetFemModelPart();
        const ProcessInfo& r_fem_process_info = r_fem_model_part.GetProcessInfo();
        ConditionsArrayType& pConditions = r_fem_model_part.GetCommunicator().LocalMesh().Conditions();

        RebuildListOfSphericParticles<SphericParticle>(r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericParticles);

        SetNormalRadiiOnAllParticles(*mpDem_model_part);

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

        RVEInitializeSolutionStep();

        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::BoundingBoxUtility(bool is_time_to_mark_and_remove) {
        KRATOS_TRY
        ModelPart& r_model_part = GetModelPart();
        ProcessInfo& r_process_info = r_model_part.GetProcessInfo();

        if (r_process_info[DOMAIN_IS_PERIODIC]) {
            mpParticleCreatorDestructor->MoveParticlesOutsideBoundingBoxBackInside(r_model_part);
        } else if (is_time_to_mark_and_remove) {
            mpParticleCreatorDestructor->DestroyParticlesOutsideBoundingBox<Cluster3D>(*mpCluster_model_part);
            mpParticleCreatorDestructor->DestroyParticlesOutsideBoundingBox<SphericParticle>(r_model_part);
        }
        if (r_process_info[CONTACT_MESH_OPTION] == 1) {
            mpParticleCreatorDestructor->MarkContactElementsForErasing(r_model_part, *mpContact_model_part);
            mpParticleCreatorDestructor->DestroyContactElements(*mpContact_model_part);
        }
        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::FinalizeSolutionStep() {

        KRATOS_TRY
        ModelPart& r_model_part = GetModelPart();
        const ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
        ElementsArrayType& rElements = r_model_part.GetCommunicator().LocalMesh().Elements();

        block_for_each(rElements, [&r_process_info](ModelPart::ElementType& rElement) {
            rElement.FinalizeSolutionStep(r_process_info);
        });

        RVEFinalizeSolutionStep();
        RVEWriteCoords();
        RVEWriteForceParticles();
        RVEWriteForceContacts();

        //if (true) AuxiliaryFunctions::ComputeReactionOnTopAndBottomSpheres(r_model_part);
        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::InitializeElements() {
        KRATOS_TRY
        ModelPart& r_model_part = GetModelPart();
        const ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
        ElementsArrayType& rElements = r_model_part.GetCommunicator().LocalMesh().Elements();

        block_for_each(rElements, [&r_process_info](ModelPart::ElementType& rElement) {
            rElement.Initialize(r_process_info);
        });

        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::InitializeDEMElements() {

        KRATOS_TRY

        ModelPart& r_model_part = GetModelPart();
        ProcessInfo& r_process_info = r_model_part.GetProcessInfo();

        double total_mass = 0.0;
        IndexPartition<unsigned int>(mListOfSphericParticles.size()).for_each([&](unsigned int i){
            mListOfSphericParticles[i]->ComputeNewRigidFaceNeighboursHistoricalData();
            mListOfSphericParticles[i]->Initialize(r_process_info);
            total_mass += mListOfSphericParticles[i]->GetMass();
        });

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

                if (submp.Has(RIGID_BODY_OPTION)) {
                    if (submp[RIGID_BODY_OPTION] == false) {
                        continue;
                    }
                }

                block_for_each(pTConditions, [&](ModelPart::ConditionType& rTCondition){
                    rTCondition.Initialize(r_process_info);
                });

                if (!r_process_info[IS_RESTARTED]){
                // Central Node
                Node::Pointer central_node;
                Geometry<Node >::PointsArrayType central_node_list;

                array_1d<double, 3> reference_coordinates = ZeroVector(3);

                if (submp.Has(RIGID_BODY_CENTER_OF_MASS)) {
                    reference_coordinates[0] = submp[RIGID_BODY_CENTER_OF_MASS][0];
                    reference_coordinates[1] = submp[RIGID_BODY_CENTER_OF_MASS][1];
                    reference_coordinates[2] = submp[RIGID_BODY_CENTER_OF_MASS][2];
                }

                int max_fem_node_id = mpParticleCreatorDestructor->FindMaxNodeIdInModelPart(fem_model_part);
                int max_dem_node_id = mpParticleCreatorDestructor->FindMaxNodeIdInModelPart(GetModelPart());
                int max_id_across_mps = std::max(max_fem_node_id, max_dem_node_id);
                int max_cluster_node_id = mpParticleCreatorDestructor->FindMaxNodeIdInModelPart(GetClusterModelPart());
                max_id_across_mps = std::max(max_id_across_mps, max_cluster_node_id);

                mpParticleCreatorDestructor->CentroidCreatorForRigidBodyElements(fem_model_part, central_node, max_id_across_mps + 1, reference_coordinates);

                central_node_list.push_back(central_node);

                int Element_Id_1 = mpParticleCreatorDestructor->FindMaxElementIdInModelPart(fem_model_part);

                Properties::Pointer properties;
                KRATOS_ERROR_IF(!submp.Has(PROPERTIES_ID))<<"PROPERTIES_ID is not set for SubModelPart "<<submp.Name()<<" . Make sure the Materials file contains material assignation for this SubModelPart"<<std::endl;
                properties = GetModelPart().GetMesh().pGetProperties(submp[PROPERTIES_ID]);

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
                submp.AddElement(RigidBodyElement3D_Kratos); //, Element_Id + 1);

                std::size_t element_id = Element_Id_1 + 1;
                std::vector<std::size_t> ElementIds;
                ElementIds.push_back(element_id);

                if (submp.Has(FREE_BODY_MOTION)) { // JIG: Backward compatibility, it should be removed in the future
                    if (submp[FREE_BODY_MOTION]) {

                        std::vector<std::vector<Node::Pointer> > thread_vectors_of_node_pointers;
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
                else {

                    // There is no need to create the rigid body elements, they already there
                    // But they need to be initialized
                    ElementsArrayType& pFemElements = fem_model_part.GetCommunicator().LocalMesh().Elements();

                    for (int k = 0; k < (int) pFemElements.size(); k++) {
                        ElementsArrayType::iterator it = pFemElements.ptr_begin() + k;
                        RigidBodyElement3D& rigid_body_element = dynamic_cast<Kratos::RigidBodyElement3D&> (*it);
                        rigid_body_element.Initialize(r_process_info);
                        rigid_body_element.CustomInitialize(submp);
                    }
                }
            }
        }

        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::CalculateConditionsRHSAndAdd() {

        KRATOS_TRY
        ClearFEMForces();
        ConditionsArrayType& rConditions = GetFemModelPart().GetCommunicator().LocalMesh().Conditions();
        ProcessInfo& r_process_info = GetFemModelPart().GetProcessInfo();
        const ProcessInfo& r_const_process_info = GetFemModelPart().GetProcessInfo();


        struct my_tls {
            Vector rhs_cond;
            Vector rhs_cond_elas;
        };

        //here the my_tls is constructed in place, which is the equivalent of "private" in OpenMP
        block_for_each(rConditions, my_tls(), [&](Condition& rCondition, my_tls& rTLS){
            Condition::GeometryType& geom = rCondition.GetGeometry();
            rCondition.CalculateRightHandSide(rTLS.rhs_cond, r_const_process_info);
            DEMWall* p_wall = dynamic_cast<DEMWall*> (&(rCondition));
            p_wall->CalculateElasticForces(rTLS.rhs_cond_elas, r_process_info);

            array_1d<double, 3> Normal_to_Element = ZeroVector(3);
            const unsigned int& dim = geom.WorkingSpaceDimension();

            if (geom.size()>2 || dim==2) p_wall->CalculateNormal(Normal_to_Element);

            for (unsigned int i = 0; i < geom.size(); i++) { //talking about each of the three nodes of the condition
                //we are studying a certain condition here
                unsigned int index = i * dim; //*2;

                array_1d<double, 3>& node_rhs = geom[i].FastGetSolutionStepValue(CONTACT_FORCES);
                array_1d<double, 3>& node_rhs_elas = geom[i].FastGetSolutionStepValue(ELASTIC_FORCES);
                array_1d<double, 3>& node_rhs_tang = geom[i].FastGetSolutionStepValue(TANGENTIAL_ELASTIC_FORCES);
                double& node_pressure = geom[i].FastGetSolutionStepValue(DEM_PRESSURE);
                array_1d<double, 3> rhs_cond_comp;
                noalias(rhs_cond_comp) = ZeroVector(3);

                geom[i].SetLock();

                for (unsigned int j = 0; j < dim; j++) { //talking about each coordinate x, y and z, loop on them
                    node_rhs[j] += rTLS.rhs_cond[index + j];
                    node_rhs_elas[j] += rTLS.rhs_cond_elas[index + j];
                    rhs_cond_comp[j] = rTLS.rhs_cond[index + j];
                }
                //node_area += 0.333333333333333 * Element_Area; //TODO: ONLY FOR TRIANGLE... Generalize for 3 or 4 nodes.
                //node_pressure actually refers to normal force. Pressure is actually computed later in function Calculate_Nodal_Pressures_and_Stresses()
                node_pressure += std::abs(GeometryFunctions::DotProduct(rhs_cond_comp, Normal_to_Element));
                noalias(node_rhs_tang) += rhs_cond_comp - GeometryFunctions::DotProduct(rhs_cond_comp, Normal_to_Element) * Normal_to_Element;

                geom[i].UnSetLock();
            }
        });

        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::ClearFEMForces() {

        KRATOS_TRY
        ModelPart& fem_model_part = GetFemModelPart();
        NodesArrayType& rNodes = fem_model_part.Nodes();

        block_for_each(rNodes, [&](ModelPart::NodeType& rNode) {

            array_1d<double, 3>& node_rhs = rNode.FastGetSolutionStepValue(CONTACT_FORCES);
            array_1d<double, 3>& node_rhs_elas = rNode.FastGetSolutionStepValue(ELASTIC_FORCES);
            array_1d<double, 3>& node_rhs_tang = rNode.FastGetSolutionStepValue(TANGENTIAL_ELASTIC_FORCES);
            double& node_pressure = rNode.GetSolutionStepValue(DEM_PRESSURE);
            //double& node_area = rNode.GetSolutionStepValue(DEM_NODAL_AREA);
            double& shear_stress = rNode.FastGetSolutionStepValue(SHEAR_STRESS);

            noalias(node_rhs) = ZeroVector(3);
            noalias(node_rhs_elas) = ZeroVector(3);
            noalias(node_rhs_tang) = ZeroVector(3);
            node_pressure = 0.0;
            //node_area = 0.0;
            shear_stress = 0.0;
        });
        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::CalculateNodalPressuresAndStressesOnWalls() {
        KRATOS_TRY

        ModelPart& fem_model_part = GetFemModelPart();
        NodesArrayType& rNodes = fem_model_part.Nodes();

        block_for_each(rNodes, [&](ModelPart::NodeType& rNode) {

            double& node_pressure = rNode.FastGetSolutionStepValue(DEM_PRESSURE);
            double node_area = rNode.FastGetSolutionStepValue(DEM_NODAL_AREA);
            double& shear_stress = rNode.FastGetSolutionStepValue(SHEAR_STRESS);
            array_1d<double, 3>& node_rhs_tang = rNode.FastGetSolutionStepValue(TANGENTIAL_ELASTIC_FORCES);
            
            if (node_area > 0.0){
                node_pressure = node_pressure / node_area;
                shear_stress = GeometryFunctions::module(node_rhs_tang) / node_area;
            }

        });
        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::SetFlagAndVariableToNodes(const Kratos::Flags& r_flag_name, ComponentOf3ComponentsVariableType& r_variable_to_set, const double value, NodesArrayType& r_nodes_array) {
        KRATOS_TRY

        block_for_each(r_nodes_array, [&](ModelPart::NodeType& rNode) {
            rNode.FastGetSolutionStepValue(r_variable_to_set) = value;
            rNode.Set(r_flag_name, true);
        });
        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::SetVariableToNodes(ComponentOf3ComponentsVariableType& r_variable_to_set, const double value, NodesArrayType& r_nodes_array) {
        KRATOS_TRY
        block_for_each(r_nodes_array, [&](ModelPart::NodeType& rNode) {
            rNode.FastGetSolutionStepValue(r_variable_to_set) = value;
        });
        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::ResetPrescribedMotionFlagsRespectingImposedDofs() {
        KRATOS_TRY
        ModelPart& r_model_part = GetModelPart();

        NodesArrayType& r_model_part_nodes = r_model_part.Nodes();

        if (!r_model_part_nodes.size()) return;

        const unsigned int vel_x_dof_position = (r_model_part.NodesBegin())->GetDofPosition(VELOCITY_X);
        const unsigned int ang_vel_x_dof_position = (r_model_part.NodesBegin())->GetDofPosition(ANGULAR_VELOCITY_X);


        block_for_each(r_model_part_nodes, [&](ModelPart::NodeType& rNode) {

            if (rNode.Is(BLOCKED)) return;
            Node& node = rNode;

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
        });
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

            if (submp.Has(RIGID_BODY_OPTION)) {
                if (submp[RIGID_BODY_OPTION] == false) {
                    continue;
                }
            }

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

                    if (submp.Has(EXTERNAL_APPLIED_FORCE)) { // JIG: Backward compatibility, it should be removed in the future
                        rigid_body_element.GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_APPLIED_FORCE)[0] = submp[EXTERNAL_APPLIED_FORCE][0];
                        rigid_body_element.GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_APPLIED_FORCE)[1] = submp[EXTERNAL_APPLIED_FORCE][1];
                        rigid_body_element.GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_APPLIED_FORCE)[2] = submp[EXTERNAL_APPLIED_FORCE][2];
                    }

                    if (submp.Has(EXTERNAL_APPLIED_MOMENT)) { // JIG: Backward compatibility, it should be removed in the future
                        rigid_body_element.GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_APPLIED_MOMENT)[0] = submp[EXTERNAL_APPLIED_MOMENT][0];
                        rigid_body_element.GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_APPLIED_MOMENT)[1] = submp[EXTERNAL_APPLIED_MOMENT][1];
                        rigid_body_element.GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_APPLIED_MOMENT)[2] = submp[EXTERNAL_APPLIED_MOMENT][2];
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

            if ((*sub_model_part).Has(RIGID_BODY_OPTION)) {
                if ((*sub_model_part)[RIGID_BODY_OPTION] == false) {
                    continue;
                }
            }

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
        IndexPartition<unsigned int>(number_of_elements).for_each([&](unsigned int i) {
            mListOfSphericParticles[i]->SetSearchRadius(amplification * (added_search_distance + mListOfSphericParticles[i]->GetRadius()));
        });

        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::SetNormalRadiiOnAllParticles(ModelPart& r_model_part) {
        KRATOS_TRY
        int number_of_elements = r_model_part.GetCommunicator().LocalMesh().ElementsArray().end() - r_model_part.GetCommunicator().LocalMesh().ElementsArray().begin();

        IndexPartition<unsigned int>(number_of_elements).for_each([&](unsigned int i){
            mListOfSphericParticles[i]->SetRadius();
        });

        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::SetSearchRadiiWithFemOnAllParticles(ModelPart& r_model_part, const double added_search_distance, const double amplification) {
        KRATOS_TRY
        int number_of_elements = r_model_part.GetCommunicator().LocalMesh().ElementsArray().end() - r_model_part.GetCommunicator().LocalMesh().ElementsArray().begin();

        IndexPartition<unsigned int>(number_of_elements).for_each([&](unsigned int i){
            mListOfSphericParticles[i]->SetSearchRadius(amplification * (added_search_distance + mListOfSphericParticles[i]->GetRadius()));
        });

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

        mpSpSearch->SearchElementsInRadiusExclusive(r_model_part, this->GetArrayOfAmplifiedRadii(), this->GetResults(), this->GetResultsDistances());

        const int number_of_particles = (int) mListOfSphericParticles.size();

        typedef std::map<SphericParticle*,std::vector<SphericParticle*>> ConnectivitiesMap;
        std::vector<ConnectivitiesMap> thread_maps_of_connectivities;
        thread_maps_of_connectivities.resize(ParallelUtilities::GetNumThreads());

        #pragma omp parallel for schedule(dynamic, 100)
        for (int i = 0; i < number_of_particles; i++) {
            mListOfSphericParticles[i]->mNeighbourElements.clear();
            for (SpatialSearch::ResultElementsContainerType::iterator neighbour_it = this->GetResults()[i].begin(); neighbour_it != this->GetResults()[i].end(); ++neighbour_it) {
                Element* p_neighbour_element = (*neighbour_it).get();
                SphericParticle* p_spheric_neighbour_particle = dynamic_cast<SphericParticle*> (p_neighbour_element);
                if (mListOfSphericParticles[i]->Is(DEMFlags::BELONGS_TO_A_CLUSTER) && (mListOfSphericParticles[i]->GetClusterId() == p_spheric_neighbour_particle->GetClusterId())) continue;
                if (mListOfSphericParticles[i]->Is(DEMFlags::POLYHEDRON_SKIN)) continue;
                mListOfSphericParticles[i]->mNeighbourElements.push_back(p_spheric_neighbour_particle);
                std::vector<SphericParticle*>& neighbours_of_this_neighbour_for_this_thread = thread_maps_of_connectivities[OpenMPUtils::ThisThread()][p_spheric_neighbour_particle];
                neighbours_of_this_neighbour_for_this_thread.push_back(mListOfSphericParticles[i]);
            }
            this->GetResults()[i].clear();
            this->GetResultsDistances()[i].clear();
        }

        // the next loop ensures consistency in neighbourhood (if A is neighbour of B, B must be neighbour of A)
        #pragma omp parallel for schedule(dynamic, 100)
        for (int i = 0; i < number_of_particles; i++) {
            auto& current_neighbours = mListOfSphericParticles[i]->mNeighbourElements;
            std::vector<SphericParticle*> neighbours_to_add;
            for (size_t k = 0; k < thread_maps_of_connectivities.size(); k++){
                ConnectivitiesMap::iterator it = thread_maps_of_connectivities[k].find(mListOfSphericParticles[i]);
                if (it != thread_maps_of_connectivities[k].end()) {
                    neighbours_to_add.insert(neighbours_to_add.end(), it->second.begin(), it->second.end());
                }
            }
            for (size_t l = 0; l < neighbours_to_add.size(); l++) {
                bool found = false;
                for (size_t m = 0; m < current_neighbours.size(); m++){
                    if (neighbours_to_add[l] == current_neighbours[m]) {
                        found = true;
                        break;
                    }
                }
                if ( found == false ) {
                    current_neighbours.push_back(neighbours_to_add[l]);
                }
            }
        }
        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::ComputeNewNeighboursHistoricalData() {

        KRATOS_TRY
        const int number_of_particles = (int) mListOfSphericParticles.size();

        #pragma omp parallel
        {
            DenseVector<int> temp_neighbours_ids;
            std::vector<array_1d<double, 3> > temp_neighbour_elastic_contact_forces;

            #pragma omp for
            for (int i = 0; i < number_of_particles; i++) {
                mListOfSphericParticles[i]->ComputeNewNeighboursHistoricalData(temp_neighbours_ids, temp_neighbour_elastic_contact_forces);
            }
        }

        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::CreateContactElements() {
        KRATOS_TRY

        std::string ElementName;
        ElementName = std::string("ParticleContactElement");
        const Element& rReferenceElement = KratosComponents<Element>::Get(ElementName);

        //Here we are going to create contact elements when we are on a target particle and we see a neighbor whose id is higher than ours.
        //We create also a pointer from the node to the element, after creating it.
        //When our particle has a higher ID than the neighbor we also create a pointer to the (previously) created contact element.
        //We proceed in this way because we want to have the pointers to contact elements in a list in the same order as the initial elements order.

        const int number_of_particles = (int) mListOfSphericParticles.size();
        int used_bonds_counter = 0;

        #pragma omp parallel
        {
            #pragma omp for
            for (int i = 0; i < number_of_particles; i++) {
                unsigned int neighbors_size = mListOfSphericParticles[i]->mNeighbourElements.size();
                mListOfSphericParticles[i]->mBondElements.resize(neighbors_size);
                for (unsigned int j = 0; j < mListOfSphericParticles[i]->mBondElements.size(); j++) {
                    mListOfSphericParticles[i]->mBondElements[j] = NULL;
                }
            }

            int private_counter = 0;
            Element::Pointer p_new_contact_element;
            #pragma omp for
            for (int i = 0; i < number_of_particles; i++) {
                bool add_new_bond = true;
                std::vector<SphericParticle*>& neighbour_elements = mListOfSphericParticles[i]->mNeighbourElements;
                unsigned int neighbors_size = mListOfSphericParticles[i]->mNeighbourElements.size();

                for (unsigned int j = 0; j < neighbors_size; j++) {
                    SphericParticle* neighbour_element = dynamic_cast<SphericParticle*> (neighbour_elements[j]);
                    if (neighbour_element == NULL) continue; //The initial neighbor was deleted at some point in time!!
                    if (mListOfSphericParticles[i]->Id() > neighbour_element->Id()) continue;

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
                        p_old_contact_element->GetGeometry()(0) = mListOfSphericParticles[i]->GetGeometry()(0);
                        p_old_contact_element->GetGeometry()(1) = neighbour_element->GetGeometry()(0);
                        p_old_contact_element->SetId(used_bonds_counter);
                        p_old_contact_element->SetProperties(mListOfSphericParticles[i]->pGetProperties());
                        ParticleContactElement* p_bond = dynamic_cast<ParticleContactElement*> (p_old_contact_element.get());
                        mListOfSphericParticles[i]->mBondElements[j] = p_bond;
                    } else {
                        Geometry<Node >::PointsArrayType NodeArray(2);
                        NodeArray.GetContainer()[0] = mListOfSphericParticles[i]->GetGeometry()(0);
                        NodeArray.GetContainer()[1] = neighbour_element->GetGeometry()(0);
                        const Properties::Pointer& properties = mListOfSphericParticles[i]->pGetProperties();
                        p_new_contact_element = rReferenceElement.Create(used_bonds_counter + 1, NodeArray, properties);

                        #pragma omp critical
                        {
                            (*mpContact_model_part).Elements().push_back(p_new_contact_element);
                            used_bonds_counter++;
                        }
                        ParticleContactElement* p_bond = dynamic_cast<ParticleContactElement*> (p_new_contact_element.get());
                        mListOfSphericParticles[i]->mBondElements[j] = p_bond;
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
                std::vector<SphericParticle*>& neighbour_elements = mListOfSphericParticles[i]->mNeighbourElements;
                unsigned int neighbors_size = mListOfSphericParticles[i]->mNeighbourElements.size();

                for (unsigned int j = 0; j < neighbors_size; j++) {
                    SphericParticle* neighbour_element = dynamic_cast<SphericParticle*> (neighbour_elements[j]);
                    if (neighbour_element == NULL) continue; //The initial neighbor was deleted at some point in time!!
                    //ATTENTION: Ghost nodes do not have mContinuumIniNeighbourElements in general, so this bond will remain as NULL!!
                    if (mListOfSphericParticles[i]->Id() < neighbour_element->Id()) continue;
                    //In all functions using mBondElements we must check that this bond is not used.

                    for (unsigned int k = 0; k < neighbour_element->mNeighbourElements.size(); k++) {
                        //ATTENTION: Ghost nodes do not have mContinuumIniNeighbourElements in general, so this bond will remain as NULL!!
                        //In all functions using mBondElements we must check that this bond is not used.
                        if (neighbour_element->mNeighbourElements[k] == NULL) continue; //The initial neighbor was deleted at some point in time!!
                        if (neighbour_element->mNeighbourElements[k]->Id() == mListOfSphericParticles[i]->Id()) {
                            ParticleContactElement* bond = neighbour_element->mBondElements[k];
                            mListOfSphericParticles[i]->mBondElements[j] = bond;
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

    void ExplicitSolverStrategy::InitializeContactElements() {

        KRATOS_TRY

        //CONTACT MODEL PART
        ElementsArrayType& pContactElements = GetAllElements(*mpContact_model_part);
        const ProcessInfo& r_process_info = GetModelPart().GetProcessInfo();

        block_for_each(pContactElements, [&r_process_info](ModelPart::ElementType& rContactElement) {
            rContactElement.Initialize(r_process_info);
        });

        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::PrepareContactElementsForPrinting() {

        ElementsArrayType& pContactElements = GetAllElements(*mpContact_model_part);

        block_for_each(pContactElements, [&](ModelPart::ElementType& rContactElement) {
            Element* raw_p_contact_element = &(rContactElement);
            ParticleContactElement* p_bond = dynamic_cast<ParticleContactElement*> (raw_p_contact_element);
            p_bond->PrepareForPrinting();
        });
    }

    void ExplicitSolverStrategy::ComputeNewRigidFaceNeighboursHistoricalData() {
        KRATOS_TRY

        IndexPartition<unsigned int>(mListOfSphericParticles.size()).for_each([&](unsigned int i){
            mListOfSphericParticles[i]->ComputeNewRigidFaceNeighboursHistoricalData();
        });

        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::SearchRigidFaceNeighbours() {
        KRATOS_TRY

        if (!mDoSearchNeighbourFEMElements) {
            return;
        }

        ElementsArrayType& pElements = mpDem_model_part->GetCommunicator().LocalMesh().Elements();
        ConditionsArrayType& pTConditions = mpFem_model_part->GetCommunicator().LocalMesh().Conditions();

        if (pTConditions.size() > 0) {
            const int number_of_particles = (int) mListOfSphericParticles.size();

            this->GetRigidFaceResults().resize(number_of_particles);
            this->GetRigidFaceResultsDistances().resize(number_of_particles);

            //Fast Bins Search
            mpDemFemSearch->SearchRigidFaceForDEMInRadiusExclusiveImplementation(pElements, pTConditions, this->GetRigidFaceResults(), this->GetRigidFaceResultsDistances());

            #pragma omp parallel for schedule(dynamic, 100)
            for (int i = 0; i < number_of_particles; i++) {
                mListOfSphericParticles[i]->mNeighbourPotentialRigidFaces.clear();
                for (ResultConditionsContainerType::iterator neighbour_it = this->GetRigidFaceResults()[i].begin(); neighbour_it != this->GetRigidFaceResults()[i].end(); ++neighbour_it) {
                    Condition* p_neighbour_condition = (*neighbour_it).get();
                    DEMWall* p_wall = dynamic_cast<DEMWall*> (p_neighbour_condition);
                    if (mListOfSphericParticles[i]->Is(DEMFlags::POLYHEDRON_SKIN)) {
                        bool must_skip_this_one = false;
                        auto& geom = p_wall->GetGeometry();
                        const unsigned int number_of_nodes = geom.size();
                        const array_1d<double, 3>& sphere_center = mListOfSphericParticles[i]->GetGeometry()[0];
                        const double epsilon = std::numeric_limits<double>::epsilon();
                        for(unsigned int k = 0; k < number_of_nodes; k++) {
                            const double distance_x = std::abs(geom[k][0] - sphere_center[0]);
                            const double distance_y = std::abs(geom[k][1] - sphere_center[1]);
                            const double distance_z = std::abs(geom[k][2] - sphere_center[2]);
                            if(distance_x < epsilon && distance_y < epsilon && distance_z < epsilon) {
                                must_skip_this_one= true;
                                break;
                            }
                        }
                        if (must_skip_this_one) continue;
                    }
                    mListOfSphericParticles[i]->mNeighbourPotentialRigidFaces.push_back(p_wall);
                }//for results iterator
                this->GetRigidFaceResults()[i].clear();
                this->GetRigidFaceResultsDistances()[i].clear();
            }

            CheckHierarchyWithCurrentNeighbours();

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

    void ExplicitSolverStrategy::CheckHierarchyWithCurrentNeighbours()
        {
        KRATOS_TRY
        const int number_of_particles = (int) mListOfSphericParticles.size();

        #pragma omp parallel
        {
            std::vector< double > Distance_Array;
            std::vector< array_1d<double, 3> > Normal_Array;
            std::vector< array_1d<double, 4> > Weight_Array;
            std::vector< int > Id_Array;
            std::vector< int > ContactType_Array;

            #pragma omp for schedule(dynamic, 100)
            for (int i = 0; i < number_of_particles; i++) {
                SphericParticle* p_sphere_i = mListOfSphericParticles[i];
                p_sphere_i->mNeighbourRigidFaces.resize(0);
                p_sphere_i->mNeighbourNonContactRigidFaces.resize(0);
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

    void ExplicitSolverStrategy::CalculateInitialMaxIndentations(const ProcessInfo& r_process_info) {
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
        ElementsArrayType& rElements = (*mpDem_model_part).GetCommunicator().LocalMesh().Elements();

        block_for_each(rElements, [&](ModelPart::ElementType& rElement) {
            Element* raw_p_element = &(rElement);
            SphericParticle* p_sphere = dynamic_cast<SphericParticle*> (raw_p_element);
            p_sphere->PrepareForPrinting(r_process_info);
        });

        KRATOS_CATCH("")
    }

    void ExplicitSolverStrategy::SynchronizeHistoricalVariables(ModelPart& r_model_part) {
        r_model_part.GetCommunicator().SynchronizeNodalSolutionStepsData();
    }

    void ExplicitSolverStrategy::SynchronizeRHS(ModelPart& r_model_part) {
        r_model_part.GetCommunicator().SynchronizeVariable(TOTAL_FORCES);
        r_model_part.GetCommunicator().SynchronizeVariable(PARTICLE_MOMENT);
    }

    double ExplicitSolverStrategy::ComputeCoordinationNumber(double& standard_dev) {
        
        KRATOS_TRY

        return 0.0;

        KRATOS_CATCH("")
    }

    //==========================================================================================================================================
    // HIERARCHICAL MULTISCALE RVE - START
    //==========================================================================================================================================

    //-----------------------------------------------------------------------------------------------------------------------------------------
    void ExplicitSolverStrategy::RVEInitialize(void) {
      ModelPart::ConditionsContainerType& r_conditions = GetFemModelPart().GetCommunicator().LocalMesh().Conditions();
      ModelPart& r_dem_model_part = GetModelPart();
      ProcessInfo& r_process_info = r_dem_model_part.GetProcessInfo();

      // Initialize properties
      mRVE_MeanRadius  = 0.0;
      mRVE_FlatWalls   = r_conditions.size() > 0;
      mRVE_Compress    = true;
      mRVE_Equilibrium = false;
      mRVE_Dimension   = r_process_info[DOMAIN_SIZE];
      mRVE_EqSteps     = 0;
      
      // Assemble vectors of wall elements
      RVEAssembleWallVectors();

      // Number of wall particles
      if (mRVE_FlatWalls) mRVE_NumParticlesWalls = 0;
      else                mRVE_NumParticlesWalls = mRVE_WallParticleXMin.size() + mRVE_WallParticleXMax.size() + mRVE_WallParticleYMin.size() + mRVE_WallParticleYMax.size() + mRVE_WallParticleZMin.size() + mRVE_WallParticleZMax.size();
      mRVE_NumParticles = mListOfSphericParticles.size() - mRVE_NumParticlesWalls;
      
      // Loop over particles
      for (int i = 0; i < mListOfSphericParticles.size(); i++) {
        // Particles movement
        mListOfSphericParticles[i]->mMoving = true;

        // Mean radius
        mRVE_MeanRadius += mListOfSphericParticles[i]->GetRadius();
      }
      mRVE_MeanRadius /= mListOfSphericParticles.size();

      // Open files
      RVEOpenFiles();

      // Read old contact forces
      RVEReadOldForces();
    }

    //-----------------------------------------------------------------------------------------------------------------------------------------
    void ExplicitSolverStrategy::RVEInitializeSolutionStep(void) {
      ModelPart&   r_dem_model_part = GetModelPart();
      ProcessInfo& r_process_info   = r_dem_model_part.GetProcessInfo();

      // Set flag for evaluating RVE in current step
      const int time_step = r_process_info[TIME_STEPS];
      const int eval_freq = r_process_info[RVE_EVAL_FREQ];
      mRVE_Solve = (time_step > 0 && time_step % eval_freq == 0.0);

      // Initialize variables
      if (mRVE_Solve && time_step > 0) {
        const int dim  = mRVE_Dimension;
        const int dim2 = mRVE_Dimension * mRVE_Dimension;

        mRVE_NumContacts             = 0;
        mRVE_NumContactsInner        = 0;
        mRVE_NumParticlesInner       = 0;
        mRVE_AvgCoordNum             = 0.0;
        mRVE_AvgCoordNumInner        = 0.0;
        mRVE_VolSolid                = 0.0;
        mRVE_WallForces              = 0.0;
        mRVE_RoseDiagram             = ZeroMatrix(2,40);
        mRVE_RoseDiagramInner        = ZeroMatrix(2,40);
        mRVE_FabricTensor            = ZeroMatrix(dim,dim);
        mRVE_FabricTensorInner       = ZeroMatrix(dim,dim);
        mRVE_CauchyTensor            = ZeroMatrix(dim,dim);
        mRVE_CauchyTensorInner       = ZeroMatrix(dim,dim);
        mRVE_TangentTensor           = ZeroMatrix(dim2,dim2);
        mRVE_TangentTensorInner      = ZeroMatrix(dim2,dim2);
        mRVE_ConductivityTensor      = ZeroMatrix(dim,dim);
        mRVE_ConductivityTensorInner = ZeroMatrix(dim,dim);
        mRVE_InnerVolParticles.clear();
        mRVE_ForceChain.clear();
      }
    }

    //-----------------------------------------------------------------------------------------------------------------------------------------
    void ExplicitSolverStrategy::RVEExecuteParticlePre(SphericParticle* p_particle) {
      p_particle->mRVESolve = mRVE_Solve;
      if (!mRVE_Solve)
        return;

      const int dim  = mRVE_Dimension;
      const int dim2 = mRVE_Dimension * mRVE_Dimension;

      if (p_particle->mWall == 0) {
        p_particle->mInner                   = (p_particle->mNeighbourRigidFaces.size() == 0);
        p_particle->mSkin                    = false;
        p_particle->mNumContacts             = 0;
        p_particle->mNumContactsInner        = 0;
        p_particle->mCoordNum                = 0;
        p_particle->mVolOverlap              = 0.0;
        p_particle->mWallForces              = 0.0;
        p_particle->mRoseDiagram             = ZeroMatrix(2,40);
        p_particle->mFabricTensor            = ZeroMatrix(dim,dim);
        p_particle->mFabricTensorInner       = ZeroMatrix(dim,dim);
        p_particle->mCauchyTensor            = ZeroMatrix(dim,dim);
        p_particle->mCauchyTensorInner       = ZeroMatrix(dim,dim);
        p_particle->mTangentTensor           = ZeroMatrix(dim2,dim2);
        p_particle->mTangentTensorInner      = ZeroMatrix(dim2,dim2);
        p_particle->mConductivityTensor      = ZeroMatrix(dim,dim);
        p_particle->mConductivityTensorInner = ZeroMatrix(dim,dim);
        p_particle->mForceChain.clear();
      }
    }

    //-----------------------------------------------------------------------------------------------------------------------------------------
    void ExplicitSolverStrategy::RVEExecuteParticlePos(SphericParticle* p_particle) {
      if (!mRVE_Solve) return;

      if (p_particle->mWall == 0) {
        mRVE_NumContacts        += p_particle->mNumContacts;
        mRVE_AvgCoordNum        += p_particle->mCoordNum;
        mRVE_VolSolid           += RVEComputeParticleVolume(p_particle) - p_particle->mVolOverlap;
        mRVE_WallForces         += p_particle->mWallForces;
        mRVE_RoseDiagram        += p_particle->mRoseDiagram;
        mRVE_FabricTensor       += p_particle->mFabricTensor;
        mRVE_CauchyTensor       += p_particle->mCauchyTensor;
        mRVE_TangentTensor      += p_particle->mTangentTensor;
        mRVE_ConductivityTensor += p_particle->mConductivityTensor;
        mRVE_ForceChain.insert(mRVE_ForceChain.end(), p_particle->mForceChain.begin(), p_particle->mForceChain.end());

        if (p_particle->mInner) {
          mRVE_NumParticlesInner++;
          mRVE_NumContactsInner        += p_particle->mNumContactsInner;
          mRVE_AvgCoordNumInner        += p_particle->mCoordNum;
          mRVE_RoseDiagramInner        += p_particle->mRoseDiagram;
          mRVE_FabricTensorInner       += p_particle->mFabricTensorInner;
          mRVE_CauchyTensorInner       += p_particle->mCauchyTensorInner;
          mRVE_TangentTensorInner      += p_particle->mTangentTensorInner;
          mRVE_ConductivityTensorInner += p_particle->mConductivityTensorInner;
          if (p_particle->mNeighbourElements.size() > 0) {
            mRVE_InnerVolParticles.push_back(p_particle);
          }
        }
        else if (p_particle->mSkin) {
          mRVE_InnerVolParticles.push_back(p_particle);
        }
      }
      else {
        mRVE_VolSolid += 0.5 * RVEComputeParticleVolume(p_particle);
      }
    }

    //-----------------------------------------------------------------------------------------------------------------------------------------
    void ExplicitSolverStrategy::RVEFinalizeSolutionStep(void) {
      if (!mRVE_Solve) return;

      // Average coordination number
      mRVE_AvgCoordNum      /= mRVE_NumParticles;
      mRVE_AvgCoordNumInner /= mRVE_NumParticlesInner;

      // Compute volume
      RVEComputeCorners();
      mRVE_VolTotal = RVEComputeTotalVolume();
      mRVE_VolInner = RVEComputeInnerVolume();

      // Compute porosity and void ratio
      RVEComputePorosity();

      // Compute stress applied by walls
      if (mRVE_FlatWalls)
        mRVE_WallStress = mRVE_WallForces / RVEComputeTotalSurface();
      else
        mRVE_WallStress = 0.0;  // TODO: Not computed for particle walls

      // Save previous values
      double prev_effect_stress = mRVE_EffectStress;
      double prev_dev_stress    = mRVE_DevStress;

      // Compute homogenized parameters
      RVEHomogenization();

      // Compute uniformity of rose diagram
      //RVEComputeRoseUniformity();

      // Write files
      RVEWriteFiles();

      // Stop compression
      RVEStopCompression();

      // Check equilibrium
      double tol = 0.0000000001;
      double ratio_eff = std::abs((mRVE_EffectStress-prev_effect_stress) / mRVE_EffectStress);
      double ratio_dev = std::abs((mRVE_DevStress-prev_dev_stress) / mRVE_DevStress);

      if (!mRVE_Compress && ratio_eff < tol && ratio_dev < tol)
        mRVE_EqSteps++;
      else
        mRVE_EqSteps = 0;

      if (mRVE_EqSteps >= 10)
        mRVE_Equilibrium = true;
    }

    //-----------------------------------------------------------------------------------------------------------------------------------------
    void ExplicitSolverStrategy::Finalize(void) {
      RVEFinalize();
    }

    //-----------------------------------------------------------------------------------------------------------------------------------------
    void ExplicitSolverStrategy::RVEFinalize(void) {
      // Close files
      RVECloseFiles();
    }

    //-----------------------------------------------------------------------------------------------------------------------------------------
    void ExplicitSolverStrategy::RVEAssembleWallVectors(void) {
      if      (mRVE_Dimension == 2 && mRVE_FlatWalls)  RVEAssembleWallVectors2D_Flat();
      else if (mRVE_Dimension == 3 && mRVE_FlatWalls)  RVEAssembleWallVectors3D_Flat();
      else if (mRVE_Dimension == 2 && !mRVE_FlatWalls) RVEAssembleWallVectors2D_Particles();
      else if (mRVE_Dimension == 3 && !mRVE_FlatWalls) RVEAssembleWallVectors3D_Particles();
    }

    //-----------------------------------------------------------------------------------------------------------------------------------------
    // TODO: Generalize
    void ExplicitSolverStrategy::RVEAssembleWallVectors2D_Flat(void) {
      ModelPart::ConditionsContainerType& r_conditions = GetFemModelPart().GetCommunicator().LocalMesh().Conditions();
      const double eps = std::numeric_limits<double>::epsilon();

      // Set all particles as non-wall
      for (int i = 0; i < (int)mListOfSphericParticles.size(); i++) {
        mListOfSphericParticles[i]->mWall = 0;
      }

      // Determine min/max wall slopes
      double slope_min =  DBL_MAX;
      double slope_max = -DBL_MAX;

      for (unsigned int i = 0; i < r_conditions.size(); i++) {
        ModelPart::ConditionsContainerType::iterator it = r_conditions.ptr_begin() + i;
        DEMWall* p_wall = dynamic_cast<DEMWall*> (&(*it));
        Condition::GeometryType& geom = p_wall->GetGeometry();
        const double x1 = geom[0][0];
        const double y1 = geom[0][1];
        const double x2 = geom[1][0];
        const double y2 = geom[1][1];

        double slope;
        if (x1==x2) slope = DBL_MAX;
        else        slope = (y2-y1)/(x2-x1);

        if (std::abs(slope) < slope_min) slope_min = std::abs(slope);
        if (std::abs(slope) > slope_max) slope_max = std::abs(slope);
      }

      // Determine X & Y walls based on slopes, and determine their min/max intersections with X & Y axes
      std::vector<DEMWall*> wall_elems_x;
      std::vector<DEMWall*> wall_elems_y;
      double intersect_xmin =  DBL_MAX;
      double intersect_xmax = -DBL_MAX;
      double intersect_ymin =  DBL_MAX;
      double intersect_ymax = -DBL_MAX;

      for (unsigned int i = 0; i < r_conditions.size(); i++) {
        ModelPart::ConditionsContainerType::iterator it = r_conditions.ptr_begin() + i;
        DEMWall* p_wall = dynamic_cast<DEMWall*> (&(*it));
        Condition::GeometryType& geom = p_wall->GetGeometry();
        const double x1 = geom[0][0];
        const double y1 = geom[0][1];
        const double x2 = geom[1][0];
        const double y2 = geom[1][1];
        
        double slope;
        if (x1==x2) slope = DBL_MAX;
        else        slope = (y2-y1)/(x2-x1);

        if (std::abs(std::abs(slope)-slope_max) < eps) {
          wall_elems_x.push_back(p_wall);
          const double intersect_x = x1 - y1 / slope;
          if (intersect_x < intersect_xmin) intersect_xmin = intersect_x;
          if (intersect_x > intersect_xmax) intersect_xmax = intersect_x;
        }
        else {
          wall_elems_y.push_back(p_wall);
          const double intersect_y = y1 - slope * x1;
          if (intersect_y < intersect_ymin) intersect_ymin = intersect_y;
          if (intersect_y > intersect_ymax) intersect_ymax = intersect_y;
        }
      }

      // Determine side (min or max) of X walls based on intersections with X axis
      for (unsigned int i = 0; i < wall_elems_x.size(); i++) {
        Condition::GeometryType& geom = wall_elems_x[i]->GetGeometry();
        const double x1 = geom[0][0];
        const double y1 = geom[0][1];
        const double x2 = geom[1][0];
        const double y2 = geom[1][1];

        double slope;
        if (x1==x2) slope = DBL_MAX;
        else        slope = (y2-y1)/(x2-x1);
        const double intersect_x = x1 - y1 / slope;

        if (std::abs(intersect_x-intersect_xmin) < eps)
          mRVE_WallXMin.push_back(wall_elems_x[i]);
        else
          mRVE_WallXMax.push_back(wall_elems_x[i]);
      }

      // Determine side (min or max) of Y walls based on intersections with Y axis
      for (unsigned int i = 0; i < wall_elems_y.size(); i++) {
        Condition::GeometryType& geom = wall_elems_y[i]->GetGeometry();
        const double x1 = geom[0][0];
        const double y1 = geom[0][1];
        const double x2 = geom[1][0];
        const double y2 = geom[1][1];

        double slope;
        if (x1==x2) slope = DBL_MAX;
        else        slope = (y2-y1)/(x2-x1);
        const double intersect_y = y1 - slope * x1;

        if (std::abs(intersect_y-intersect_ymin) < eps)
          mRVE_WallYMin.push_back(wall_elems_y[i]);
        else
          mRVE_WallYMax.push_back(wall_elems_y[i]);
      }
    }

    //-----------------------------------------------------------------------------------------------------------------------------------------
    // TODO: Adapt it to consider non square RVEs, e.g. RVEs with shear applied.
    void ExplicitSolverStrategy::RVEAssembleWallVectors3D_Flat(void) {
      const double eps = std::numeric_limits<double>::epsilon();

      std::vector<DEMWall*> wall_elems_x;
      std::vector<DEMWall*> wall_elems_y;
      std::vector<DEMWall*> wall_elems_z;

      double xmin =  DBL_MAX;
      double xmax = -DBL_MAX;
      double ymin =  DBL_MAX;
      double ymax = -DBL_MAX;
      double zmin =  DBL_MAX;
      double zmax = -DBL_MAX;

      ModelPart::ConditionsContainerType& r_conditions = GetFemModelPart().GetCommunicator().LocalMesh().Conditions();
      for (unsigned int i = 0; i < r_conditions.size(); i++) {
        ModelPart::ConditionsContainerType::iterator it = r_conditions.ptr_begin() + i;
        DEMWall* p_wall = dynamic_cast<DEMWall*> (&(*it));

        Condition::GeometryType& geom = p_wall->GetGeometry();
        double coords1[3] = {geom[0][0], geom[0][1], geom[0][2]};
        double coords2[3] = {geom[1][0], geom[1][1], geom[1][2]};
        double coords3[3] = {geom[2][0], geom[2][1], geom[2][2]};

        if (std::abs(coords1[0]-coords2[0]) < eps && std::abs(coords1[0]-coords3[0]) < eps) {
          wall_elems_x.push_back(p_wall);
          if (coords1[0] < xmin) xmin = coords1[0];
          if (coords1[0] > xmax) xmax = coords1[0];
        }
        else if (std::abs(coords1[1]-coords2[1]) < eps && std::abs(coords1[1]-coords3[1]) < eps) {
          wall_elems_y.push_back(p_wall);
          if (coords1[1] < ymin) ymin = coords1[1];
          if (coords1[1] > ymax) ymax = coords1[1];
        }
        else if (std::abs(coords1[2]-coords2[2]) < eps && std::abs(coords1[2]-coords3[2]) < eps) {
          wall_elems_z.push_back(p_wall);
          if (coords1[2] < zmin) zmin = coords1[2];
          if (coords1[2] > zmax) zmax = coords1[2];
        }
      }

      for (unsigned int i = 0; i < wall_elems_x.size(); i++) {
        const double x = wall_elems_x[i]->GetGeometry()[0][0];
        if      (std::abs(x-xmin) < eps) mRVE_WallXMin.push_back(wall_elems_x[i]);
        else if (std::abs(x-xmax) < eps) mRVE_WallXMax.push_back(wall_elems_x[i]);
      }
      for (unsigned int i = 0; i < wall_elems_y.size(); i++) {
        const double y = wall_elems_y[i]->GetGeometry()[0][1];
        if      (std::abs(y-ymin) < eps) mRVE_WallYMin.push_back(wall_elems_y[i]);
        else if (std::abs(y-ymax) < eps) mRVE_WallYMax.push_back(wall_elems_y[i]);
      }
      for (unsigned int i = 0; i < wall_elems_z.size(); i++) {
        const double z = wall_elems_z[i]->GetGeometry()[0][2];
        if      (std::abs(z-zmin) < eps) mRVE_WallZMin.push_back(wall_elems_z[i]);
        else if (std::abs(z-zmax) < eps) mRVE_WallZMax.push_back(wall_elems_z[i]);
      }

      for (int i = 0; i < (int)mListOfSphericParticles.size(); i++) {
        mListOfSphericParticles[i]->mWall = 0;
      }
    }

    //-----------------------------------------------------------------------------------------------------------------------------------------
    // TODO: Adapt it to consider non square RVEs, e.g. RVEs with shear applied.
    void ExplicitSolverStrategy::RVEAssembleWallVectors2D_Particles(void) {
      const double eps = std::numeric_limits<double>::epsilon();
      double xmin =  DBL_MAX;
      double xmax = -DBL_MAX;
      double ymin =  DBL_MAX;
      double ymax = -DBL_MAX;

      const int number_of_particles = (int)mListOfSphericParticles.size();
      for (int i = 0; i < number_of_particles; i++) {
        const double coord_x = mListOfSphericParticles[i]->GetGeometry()[0][0];
        const double coord_y = mListOfSphericParticles[i]->GetGeometry()[0][1];
        if (coord_x < xmin) xmin = coord_x;
        if (coord_x > xmax) xmax = coord_x;
        if (coord_y < ymin) ymin = coord_y;
        if (coord_y > ymax) ymax = coord_y;
      }

      for (int i = 0; i < number_of_particles; i++) {
        const double coord_x = mListOfSphericParticles[i]->GetGeometry()[0][0];
        const double coord_y = mListOfSphericParticles[i]->GetGeometry()[0][1];
        if (std::abs(coord_x - xmin) < eps) {
          mRVE_WallParticleXMin.push_back(mListOfSphericParticles[i]);
          mListOfSphericParticles[i]->mWall = 1;
        }
        else if (std::abs(coord_x - xmax) < eps) {
          mRVE_WallParticleXMax.push_back(mListOfSphericParticles[i]);
          mListOfSphericParticles[i]->mWall = 2;
        }
        else if (std::abs(coord_y - ymin) < eps) {
          mRVE_WallParticleYMin.push_back(mListOfSphericParticles[i]);
          mListOfSphericParticles[i]->mWall = 3;
        }
        else if (std::abs(coord_y - ymax) < eps) {
          mRVE_WallParticleYMax.push_back(mListOfSphericParticles[i]);
          mListOfSphericParticles[i]->mWall = 4;
        }
        else {
          mListOfSphericParticles[i]->mWall = 0;
        }
      }
    }

    //-----------------------------------------------------------------------------------------------------------------------------------------
    // TODO: Adapt it to consider non square RVEs, e.g. RVEs with shear applied.
    void ExplicitSolverStrategy::RVEAssembleWallVectors3D_Particles(void) {
      const double eps = std::numeric_limits<double>::epsilon();
      double xmin =  DBL_MAX;
      double xmax = -DBL_MAX;
      double ymin =  DBL_MAX;
      double ymax = -DBL_MAX;
      double zmin =  DBL_MAX;
      double zmax = -DBL_MAX;

      const int number_of_particles = (int)mListOfSphericParticles.size();
      for (int i = 0; i < number_of_particles; i++) {
        const double coord_x = mListOfSphericParticles[i]->GetGeometry()[0][0];
        const double coord_y = mListOfSphericParticles[i]->GetGeometry()[0][1];
        const double coord_z = mListOfSphericParticles[i]->GetGeometry()[0][2];
        if (coord_x < xmin) xmin = coord_x;
        if (coord_x > xmax) xmax = coord_x;
        if (coord_y < ymin) ymin = coord_y;
        if (coord_y > ymax) ymax = coord_y;
        if (coord_z < zmin) zmin = coord_z;
        if (coord_z > zmax) zmax = coord_z;
      }

      for (int i = 0; i < number_of_particles; i++) {
        const double coord_x = mListOfSphericParticles[i]->GetGeometry()[0][0];
        const double coord_y = mListOfSphericParticles[i]->GetGeometry()[0][1];
        const double coord_z = mListOfSphericParticles[i]->GetGeometry()[0][2];
        if (std::abs(coord_x - xmin) < eps) {
          mRVE_WallParticleXMin.push_back(mListOfSphericParticles[i]);
          mListOfSphericParticles[i]->mWall = 1;
        }
        else if (std::abs(coord_x - xmax) < eps) {
          mRVE_WallParticleXMax.push_back(mListOfSphericParticles[i]);
          mListOfSphericParticles[i]->mWall = 2;
        }
        else if (std::abs(coord_y - ymin) < eps) {
          mRVE_WallParticleYMin.push_back(mListOfSphericParticles[i]);
          mListOfSphericParticles[i]->mWall = 3;
        }
        else if (std::abs(coord_y - ymax) < eps) {
          mRVE_WallParticleYMax.push_back(mListOfSphericParticles[i]);
          mListOfSphericParticles[i]->mWall = 4;
        }
        else if (std::abs(coord_z - zmin) < eps) {
          mRVE_WallParticleZMin.push_back(mListOfSphericParticles[i]);
          mListOfSphericParticles[i]->mWall = 5;
        }
        else if (std::abs(coord_z - zmax) < eps) {
          mRVE_WallParticleZMax.push_back(mListOfSphericParticles[i]);
          mListOfSphericParticles[i]->mWall = 6;
        }
        else {
          mListOfSphericParticles[i]->mWall = 0;
        }
      }
    }

    //-----------------------------------------------------------------------------------------------------------------------------------------
    // TODO: Generalize - Specific for shear tests
    void ExplicitSolverStrategy::RVEComputeCorners(void) {
      const double eps = std::numeric_limits<double>::epsilon();
      const double xmin_x1 = mRVE_WallXMin[0]->GetGeometry()[0][0]; const double xmin_x2 = mRVE_WallXMin[0]->GetGeometry()[1][0];
      const double xmin_y1 = mRVE_WallXMin[0]->GetGeometry()[0][1]; const double xmin_y2 = mRVE_WallXMin[0]->GetGeometry()[1][1];
      const double xmax_x1 = mRVE_WallXMax[0]->GetGeometry()[0][0]; const double xmax_x2 = mRVE_WallXMax[0]->GetGeometry()[1][0];
      const double ymin_y1 = mRVE_WallYMin[0]->GetGeometry()[0][1]; const double ymin_y2 = mRVE_WallYMin[0]->GetGeometry()[1][1];
      const double ymax_y1 = mRVE_WallYMax[0]->GetGeometry()[0][1]; const double ymax_y2 = mRVE_WallYMax[0]->GetGeometry()[1][1];

      bool is_square = (std::abs(xmin_x1-xmin_x2)<eps && std::abs(xmax_x1-xmax_x2)<eps && std::abs(ymin_y1-ymin_y2)<eps && std::abs(ymax_y1-ymax_y2)<eps);

      if (1) { //(is_square) {
        mRVE_CornerCoordsX = {xmin_x1, xmax_x1, xmax_x1, xmin_x1};
        mRVE_CornerCoordsY = {ymin_y1, ymin_y1, ymax_y1, ymax_y1};
      }
      else { // sheared RVE
        double ymin_x_min =  DBL_MAX;
        double ymin_x_max = -DBL_MAX;
        for (unsigned int i = 0; i < mRVE_WallYMin.size(); i++) {
          Condition::GeometryType& geom = mRVE_WallYMin[i]->GetGeometry();
          const double ymin_x1 = geom[0][0];
          const double ymin_x2 = geom[1][0];
          if (ymin_x1 < ymin_x_min) ymin_x_min = ymin_x1;
          if (ymin_x1 > ymin_x_max) ymin_x_max = ymin_x1;
          if (ymin_x2 < ymin_x_min) ymin_x_min = ymin_x2;
          if (ymin_x2 > ymin_x_max) ymin_x_max = ymin_x2;
        }

        const double slope = (xmin_x2-xmin_x1) / (xmin_y2-xmin_y1);
        const double height = ymax_y1 - ymin_y1;
        const double ymax_x_min = ymin_x_min + slope * height;
        const double ymax_x_max = ymin_x_max + slope * height;

        mRVE_CornerCoordsX = {ymin_x_min, ymin_x_max, ymax_x_min, ymax_x_max};
        mRVE_CornerCoordsY = {ymin_y1,    ymin_y1,    ymax_y1,    ymax_y1};
      }
    }

    //-----------------------------------------------------------------------------------------------------------------------------------------
    // TODO: Missing 3D; 2D Do not work with extended walls
    double ExplicitSolverStrategy::RVEComputeTotalSurface(void) {
      double surface = 0.0;

      if (mRVE_WallXMin.size() > 0 || mRVE_WallXMax.size() > 0 ||
          mRVE_WallYMin.size() > 0 || mRVE_WallYMax.size() > 0 ||
          mRVE_WallZMin.size() > 0 || mRVE_WallZMax.size() > 0)
      {
        if (mRVE_Dimension == 2) { // Perimeter
          const double x1 = mRVE_CornerCoordsX[0]; const double y1 = mRVE_CornerCoordsY[0];
          const double x2 = mRVE_CornerCoordsX[1]; const double y2 = mRVE_CornerCoordsY[1];
          const double x3 = mRVE_CornerCoordsX[2]; const double y3 = mRVE_CornerCoordsY[2];
          const double x4 = mRVE_CornerCoordsX[3]; const double y4 = mRVE_CornerCoordsY[3];
          
          const double line1 = sqrt(pow(x2-x1,2)+pow(y2-y1,2));
          const double line2 = sqrt(pow(x3-x2,2)+pow(y3-y2,2));
          const double line3 = sqrt(pow(x4-x3,2)+pow(y4-y3,2));
          const double line4 = sqrt(pow(x1-x4,2)+pow(y1-y4,2));
          
          surface = line1 + line2 + line3 + line4;
        }
        else if (mRVE_Dimension == 3) {
          //double dX = std::abs(mRVE_WallXMin[0]->GetGeometry()[0][0] - mRVE_WallXMax[0]->GetGeometry()[0][0]);
          //double dY = std::abs(mRVE_WallYMin[0]->GetGeometry()[0][1] - mRVE_WallYMax[0]->GetGeometry()[0][1]);
          //double dZ = std::abs(mRVE_WallZMin[0]->GetGeometry()[0][2] - mRVE_WallZMax[0]->GetGeometry()[0][2]);
          //surface = 2.0 * (dX * dY + dX * dZ + dY * dZ);
        }
      }

      return surface;
    }

    //-----------------------------------------------------------------------------------------------------------------------------------------
    // TODO: Missing 3D; 2D Do not work with extended walls
    double ExplicitSolverStrategy::RVEComputeTotalVolume(void) {
      double vol = 0.0;

      if (mRVE_FlatWalls) {
        const double x1 = mRVE_CornerCoordsX[0]; const double y1 = mRVE_CornerCoordsY[0];
        const double x2 = mRVE_CornerCoordsX[1]; const double y2 = mRVE_CornerCoordsY[1];
        const double x3 = mRVE_CornerCoordsX[2]; const double y3 = mRVE_CornerCoordsY[2];
        const double x4 = mRVE_CornerCoordsX[3]; const double y4 = mRVE_CornerCoordsY[3];

        vol = 0.5 * std::abs(x1*y2 + x2*y3 + x3*y4 + x4*y1 - x2*y1 - x3*y2 - x4*y3 - x1*y4); // shoelace formula
      }
      else {
        //if (mRVE_WallParticleXMin.size() > 0 && mRVE_WallParticleXMax.size() > 0)
        //  dX = std::abs(mRVE_WallParticleXMin[0]->GetGeometry()[0][0] - mRVE_WallParticleXMax[0]->GetGeometry()[0][0]);
        //if (mRVE_WallParticleYMin.size() > 0 && mRVE_WallParticleYMax.size() > 0)
        //  dY = std::abs(mRVE_WallParticleYMin[0]->GetGeometry()[0][1] - mRVE_WallParticleYMax[0]->GetGeometry()[0][1]);
        //if (mRVE_WallParticleZMin.size() > 0 && mRVE_WallParticleZMax.size() > 0)
        //  dZ = std::abs(mRVE_WallParticleZMin[0]->GetGeometry()[0][2] - mRVE_WallParticleZMax[0]->GetGeometry()[0][2]);
      }

      return vol;
    }

    //-----------------------------------------------------------------------------------------------------------------------------------------
    double ExplicitSolverStrategy::RVEComputeInnerVolume(void) {
      // Used to compute inner volume with convex hull,
      // but since this is not a reliable apporach in Kratos anymore,
      // the inner volume is simply being computed with the porosity offset to save time

      ModelPart& r_dem_model_part = GetModelPart();
      ProcessInfo& r_process_info = r_dem_model_part.GetProcessInfo();
      const double os   = r_process_info[INNER_CONSOLIDATION_POROSITY_OFFSET] * mRVE_MeanRadius;
      const double xmin = mRVE_WallXMin[0]->GetGeometry()[0][0] + os;
      const double xmax = mRVE_WallXMax[0]->GetGeometry()[0][0] - os;
      const double ymin = mRVE_WallYMin[0]->GetGeometry()[0][1] + os;
      const double ymax = mRVE_WallYMax[0]->GetGeometry()[0][1] - os;
      return (xmax - xmin) * (ymax - ymin);
      
      //const int num_particles = mRVE_InnerVolParticles.size();
      //if (!mRVE_FlatWalls || num_particles < 3)
      //  return mRVE_VolTotal;
      //
      //// Create and clear IO
      //std::string switches = "PQ";
      //struct triangulateio in, out, vorout;
      //ClearTriangle(in);
      //ClearTriangle(out);
      //ClearTriangle(vorout);

      //// Build input
      //in.numberofpoints = num_particles;
      //in.pointlist = (double*)malloc(sizeof(double) * in.numberofpoints * 2);

      //for (int i = 0; i < in.numberofpoints; i++) {
      //  array_1d<double,3> coords = mRVE_InnerVolParticles[i]->GetGeometry()[0].Coordinates();
      //  in.pointlist[2 * i + 0] = coords[0];
      //  in.pointlist[2 * i + 1] = coords[1];
      //}

      //// Perform triangulation
      //int fail = 0;
      //try {
      //  triangulate(&switches[0], &in, &out, &vorout);
      //}
      //catch (int error_code) {
      //  fail = error_code;
      //}

      //if (fail || out.numberoftriangles == 0 || in.numberofpoints != out.numberofpoints) {
      //  KRATOS_ERROR << "Fail to generate triangulation!" << std::endl;
      //  FreeTriangle(in);
      //  FreeTriangle(out);
      //  FreeTriangle(vorout);
      //  return 0.0;
      //}

      //// Compute volume of convex hull (area in 2D)
      //double total_volume = 0.0;
      //for (int i = 0; i < out.numberoftriangles; i++) {
      //  // Get vertices IDs
      //  const int v1 = out.trianglelist[3 * i + 0] - 1;
      //  const int v2 = out.trianglelist[3 * i + 1] - 1;
      //  const int v3 = out.trianglelist[3 * i + 2] - 1;

      //  // Get vertices coordinates
      //  const double x1 = out.pointlist[2 * v1 + 0];
      //  const double y1 = out.pointlist[2 * v1 + 1];
      //  const double x2 = out.pointlist[2 * v2 + 0];
      //  const double y2 = out.pointlist[2 * v2 + 1];
      //  const double x3 = out.pointlist[2 * v3 + 0];
      //  const double y3 = out.pointlist[2 * v3 + 1];

      //  // Add triangle area
      //  total_volume += fabs(0.5 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)));
      //}

      //// Free memory
      //FreeTriangle(in);
      //FreeTriangle(out);
      //FreeTriangle(vorout);

      //ClearTriangle(in);
      //ClearTriangle(out);
      //ClearTriangle(vorout);

      //return total_volume;
    }

    //-----------------------------------------------------------------------------------------------------------------------------------------
    double ExplicitSolverStrategy::RVEComputeParticleVolume(SphericParticle* p_particle) {
      if      (mRVE_Dimension == 2) return Globals::Pi * p_particle->GetRadius() * p_particle->GetRadius(); // Area
      else if (mRVE_Dimension == 3) return p_particle->CalculateVolume();
      else return 0.0;
    }

    //-----------------------------------------------------------------------------------------------------------------------------------------
    void ExplicitSolverStrategy::RVEComputePorosity(void) {
      mRVE_Porosity = 1.0 - mRVE_VolSolid / mRVE_VolTotal;
      mRVE_VoidRatio = mRVE_Porosity / (1.0 - mRVE_Porosity);
      RVEComputeInnerPorosity();
    }

    //-----------------------------------------------------------------------------------------------------------------------------------------
    // TODO: Missing 3D; DO NOT WORK FOR NON-SQUARE SHAPES;
    void ExplicitSolverStrategy::RVEComputeInnerPorosity(void) {
      ModelPart& r_dem_model_part = GetModelPart();
      ProcessInfo& r_process_info = r_dem_model_part.GetProcessInfo();
      const double eps = std::numeric_limits<double>::epsilon();

      const double x1 = mRVE_CornerCoordsX[0]; const double y1 = mRVE_CornerCoordsY[0];
      const double x2 = mRVE_CornerCoordsX[1]; const double y2 = mRVE_CornerCoordsY[1];
      const double x3 = mRVE_CornerCoordsX[2]; const double y3 = mRVE_CornerCoordsY[2];
      const double x4 = mRVE_CornerCoordsX[3]; const double y4 = mRVE_CornerCoordsY[3];

      bool is_square = (std::abs(y1-y2)<eps && std::abs(y3-y4)<eps && std::abs(x1-x4)<eps && std::abs(x2-x3)<eps);

      if (!r_process_info[INNER_CONSOLIDATION_POROSITY] || !is_square) {
        mRVE_VolSolidInner  = 0.0;
        mRVE_PorosityInner  = 0.0;
        mRVE_VoidRatioInner = 0.0;
        return;
      }
      
      // Inner boundaries and corners
      const double offset = r_process_info[INNER_CONSOLIDATION_POROSITY_OFFSET] * mRVE_MeanRadius;

      const double xmin = mRVE_WallXMin[0]->GetGeometry()[0][0] + offset;
      const double xmax = mRVE_WallXMax[0]->GetGeometry()[0][0] - offset;
      const double ymin = mRVE_WallYMin[0]->GetGeometry()[0][1] + offset;
      const double ymax = mRVE_WallYMax[0]->GetGeometry()[0][1] - offset;

      const double corner_1[3] = { xmin, ymin, 0.0 };
      const double corner_2[3] = { xmax, ymin, 0.0 };
      const double corner_3[3] = { xmin, ymax, 0.0 };
      const double corner_4[3] = { xmax, ymax, 0.0 };

      // Inner volume
      const double inner_vol = (xmax - xmin) * (ymax - ymin);

      // Compute solid volume
      mRVE_VolSolidInner = 0.0;

      // Loop over all particles
      for (int i = 0; i < mListOfSphericParticles.size(); i++) {
        const double r = mListOfSphericParticles[i]->GetRadius();
        auto& central_node = mListOfSphericParticles[i]->GetGeometry()[0];
        const double x = central_node[0];
        const double y = central_node[1];
        const double z = central_node[2];
        const double coords[3] = { x, y, z };

        // Distances from particle center to corners
        const double dist_corner_1 = GeometryFunctions::DistanceOfTwoPoint(coords, corner_1);
        const double dist_corner_2 = GeometryFunctions::DistanceOfTwoPoint(coords, corner_2);
        const double dist_corner_3 = GeometryFunctions::DistanceOfTwoPoint(coords, corner_3);
        const double dist_corner_4 = GeometryFunctions::DistanceOfTwoPoint(coords, corner_4);

        // Area inside inner region
        double area_inside = 0.0;
        
        // 1 - Particle inside region
        if (x + r <= xmax && x - r >= xmin && y + r <= ymax && y - r >= ymin) {
          area_inside = Globals::Pi * r * r;
        }

        // 2A - Cell corner (bottom-left) inside circle
        else if (dist_corner_1 < r) {
          const double dx = x - xmin;
          const double dy = y - ymin;
          const double sx = sqrt(r * r - dy * dy);
          const double sy = sqrt(r * r - dx * dx);
          const double lx = sx + dx;
          const double ly = sy + dy;
          const double h = sqrt(lx * lx + ly * ly);
          const double t = 2 * asin(h / (2 * r));
          const double area_tri = lx * ly / 2;
          const double area_seg = 0.5 * r * r * (t - sin(t));
          area_inside = area_tri + area_seg;
        }

        // 2B - Cell corner (bottom-right) inside circle
        else if (dist_corner_2 < r) {
          const double dx = x - xmax;
          const double dy = y - ymin;
          const double sx = sqrt(r * r - dy * dy);
          const double sy = sqrt(r * r - dx * dx);
          const double lx = sx - dx;
          const double ly = sy + dy;
          const double h = sqrt(lx * lx + ly * ly);
          const double t = 2 * asin(h / (2 * r));
          const double area_tri = lx * ly / 2;
          const double area_seg = 0.5 * r * r * (t - sin(t));
          area_inside = area_tri + area_seg;
        }

        // 2C - Cell corner (top-left) inside circle
        else if (dist_corner_3 < r) {
          const double dx = x - xmin;
          const double dy = y - ymax;
          const double sx = sqrt(r * r - dy * dy);
          const double sy = sqrt(r * r - dx * dx);
          const double lx = sx + dx;
          const double ly = sy - dy;
          const double h = sqrt(lx * lx + ly * ly);
          const double t = 2 * asin(h / (2 * r));
          const double area_tri = lx * ly / 2;
          const double area_seg = 0.5 * r * r * (t - sin(t));
          area_inside = area_tri + area_seg;
        }

        // 2D - Cell corner (top-right) inside circle
        else if (dist_corner_4 < r) {
          const double dx = x - xmax;
          const double dy = y - ymax;
          const double sx = sqrt(r * r - dy * dy);
          const double sy = sqrt(r * r - dx * dx);
          const double lx = sx - dx;
          const double ly = sy - dy;
          const double h = sqrt(lx * lx + ly * ly);
          const double t = 2 * asin(h / (2 * r));
          const double area_tri = lx * ly / 2;
          const double area_seg = 0.5 * r * r * (t - sin(t));
          area_inside = area_tri + area_seg;
        }

        // Particle crossed by edges
        else {
          bool is_inside = false;

          // Left edge
          if (std::abs(x - xmin) < r && y > ymin && y < ymax) {
            const double ident = r - std::abs(x - xmin);
            const double ovlp  = r * r * acos((r - ident) / r) - (r - ident) * sqrt(2 * r * ident - ident * ident);
            if (x > xmin) {
              is_inside = 1;
              area_inside = area_inside - ovlp;
            }
            else {
              area_inside = ovlp;
            }
          }

          // Right edge
          if (std::abs(x - xmax) < r && y > ymin && y < ymax) {
            const double ident = r - std::abs(x - xmax);
            const double ovlp  = r * r * acos((r - ident) / r) - (r - ident) * sqrt(2 * r * ident - ident * ident);
            if (x < xmax) {
              is_inside = 1;
              area_inside = area_inside - ovlp;
            }
            else {
              area_inside = ovlp;
            }
          }

          // Lower edge
          if (std::abs(y - ymin) < r && x > xmin && x < xmax) {
            const double ident = r - std::abs(y - ymin);
            const double ovlp = r * r * acos((r - ident) / r) - (r - ident) * sqrt(2 * r * ident - ident * ident);
            if (y > ymin) {
              is_inside = 1;
              area_inside = area_inside - ovlp;
            }
            else {
              area_inside = ovlp;
            }
          }

          // Top edge
          if (std::abs(y - ymax) < r && x > xmin && x < xmax) {
            const double ident = r - std::abs(y - ymax);
            const double ovlp = r * r * acos((r - ident) / r) - (r - ident) * sqrt(2 * r * ident - ident * ident);
            if (y < ymax) {
              is_inside = 1;
              area_inside = area_inside - ovlp;
            }
            else {
              area_inside = ovlp;
            }
          }

          // Add particle area if it is inside cell
          if (is_inside)
            area_inside = area_inside + Globals::Pi * r * r;
        }

        // Remove overlaps: loop over possible particle neighbors
        for (int j = i+1; j < mListOfSphericParticles.size(); j++) {
          const double rn = mListOfSphericParticles[j]->GetRadius();
          const double xn = mListOfSphericParticles[j]->GetGeometry()[0][0];
          const double yn = mListOfSphericParticles[j]->GetGeometry()[0][1];
          const double zn = mListOfSphericParticles[j]->GetGeometry()[0][2];

          // Identation
          auto& neighbour_node = mListOfSphericParticles[j]->GetGeometry()[0];
          array_1d<double, 3> v;
          noalias(v) = central_node.Coordinates() - neighbour_node.Coordinates();
          const double d = DEM_MODULUS_3(v);
          const double ident = r + rn - d;

          if (ident <= 0.0)
            continue;

          // Overlap area
          double area_ovlp = r*r * std::acos((d*d + r*r - rn*rn) / (2.0 * d * r)) + rn*rn * std::acos((d*d + rn*rn - r*r) / (2.0 * d * rn)) - 0.5 * sqrt((d + r + rn) * (-d + r + rn) * (d + r - rn) * (d - r + rn));

          // Unit vector
          array_1d<double, 3> vu = v;
          GeometryFunctions::normalize(vu);

          // Contact radius
          const double rc = sqrt(r*r - pow((r*r - rn* rn + d*d)/(2*d),2.0));

          // Coordinates of overlap center
          const double d_1O = sqrt(r*r - rc*rc);
          array_1d<double,3> vu_d1O = vu;
          DEM_MULTIPLY_BY_SCALAR_3(vu_d1O, d_1O);
          array_1d<double, 3> pO;
          noalias(pO) = central_node.Coordinates() + vu_d1O;

          // Coordiantes of overlap ends
          array_1d<double, 3> vu_OA;
          array_1d<double, 3> vu_OB;
          array_1d<double, 3> positive_z = ZeroVector(3);
          positive_z[2] = 1.0;
          GeometryFunctions::CrossProduct(positive_z, vu, vu_OA);
          GeometryFunctions::CrossProduct(vu, positive_z, vu_OB);

          array_1d<double, 3> v_OA = vu_OA;
          array_1d<double, 3> v_OB = vu_OB;
          DEM_MULTIPLY_BY_SCALAR_3(v_OA, rc);
          DEM_MULTIPLY_BY_SCALAR_3(v_OB, rc);

          array_1d<double, 3> pA;
          array_1d<double, 3> pB;
          noalias(pA) = pO + v_OA;
          noalias(pB) = pO + v_OB;

          // Check if overlap area is inside region or cut by edges
          const double x1 = GeometryFunctions::min(pA[0], pB[0]);
          const double x2 = GeometryFunctions::max(pA[0], pB[0]);
          const double y1 = GeometryFunctions::min(pA[1], pB[1]);
          const double y2 = GeometryFunctions::max(pA[1], pB[1]);
          double ratio_in = 0.0;

          if      (x1 > xmin && x2 < xmax && y1 > ymin && y2 < ymax) ratio_in = 1.0;                     // Inside
          else if (x1 < xmin && x2 > xmin && y1 > ymin && y2 < ymax) ratio_in = (x2 - xmin) / (x2 - x1); // Left edge
          else if (x1 < xmax && x2 > xmax && y1 > ymin && y2 < ymax) ratio_in = (xmax - x1) / (x2 - x1); // Right edge
          else if (y1 < ymin && y2 > ymin && x1 > xmin && x2 < xmax) ratio_in = (y2 - ymin) / (y2 - y1); // Bottom edge
          else if (y1 < ymax && y2 > ymax && x1 > xmin && x2 < xmax) ratio_in = (ymax - y1) / (y2 - y1); // Top edge

          area_ovlp *= ratio_in;

          // Discount overlap
          area_inside -= area_ovlp;
        }

        // Add particle area to solid volume
        mRVE_VolSolidInner += area_inside;
      }

      // Compute porosity
      mRVE_PorosityInner = 1.0 - mRVE_VolSolidInner / inner_vol;
      mRVE_VoidRatioInner = mRVE_PorosityInner / (1.0 - mRVE_PorosityInner);
    }

    //-----------------------------------------------------------------------------------------------------------------------------------------
    void ExplicitSolverStrategy::RVEHomogenization(void) {
      if (mRVE_NumContactsInner == 0.0) {
        const int dim  = mRVE_Dimension;
        const int dim2 = mRVE_Dimension * mRVE_Dimension;

        mRVE_FabricTensor            = ZeroMatrix(dim,dim);
        mRVE_FabricTensorInner       = ZeroMatrix(dim,dim);
        mRVE_CauchyTensor            = ZeroMatrix(dim,dim);
        mRVE_CauchyTensorInner       = ZeroMatrix(dim,dim);
        mRVE_TangentTensor           = ZeroMatrix(dim2,dim2);
        mRVE_TangentTensorInner      = ZeroMatrix(dim2,dim2);
        mRVE_ConductivityTensor      = ZeroMatrix(dim,dim);
        mRVE_ConductivityTensorInner = ZeroMatrix(dim,dim);
        mRVE_Anisotropy              = 0.0;
        mRVE_AnisotropyInner         = 0.0;
        mRVE_EffectStress            = 0.0;
        mRVE_EffectStressInner       = 0.0;
        mRVE_DevStress               = 0.0;
        mRVE_DevStressInner          = 0.0;
      }

      else {
        double deviatoric_fabric;
        double deviatoric_fabric_inner;
        double deviatoric_stress;
        double deviatoric_stress_inner;
        double cauchy_trace             = 0.0;
        double cauchy_trace_inner       = 0.0;
        double double_dot_product       = 0.0;
        double double_dot_product_inner = 0.0;

        for (int i = 0; i < mRVE_Dimension; i++) {
          for (int j = 0; j < mRVE_Dimension; j++) {
            for (int k = 0; k < mRVE_Dimension; k++) {
              for (int l = 0; l < mRVE_Dimension; l++) {
                const int idx_1 = 2 * i + j;
                const int idx_2 = 2 * k + l;
                mRVE_TangentTensor(idx_1,idx_2)      /= mRVE_VolTotal;
                mRVE_TangentTensorInner(idx_1,idx_2) /= mRVE_VolInner;
              }
            }
            mRVE_FabricTensor(i,j)      /= mRVE_NumContacts;
            mRVE_FabricTensorInner(i,j) /= mRVE_NumContactsInner;

            mRVE_CauchyTensor(i,j)      /= mRVE_VolTotal;
            mRVE_CauchyTensorInner(i,j) /= mRVE_VolInner;

            mRVE_ConductivityTensor(i,j)      /= mRVE_VolTotal;
            mRVE_ConductivityTensorInner(i,j) /= mRVE_VolInner;

            if (i == j) {
              deviatoric_fabric       = 4.0 * (mRVE_FabricTensor(i,j)      - (1.0 / mRVE_Dimension));
              deviatoric_fabric_inner = 4.0 * (mRVE_FabricTensorInner(i,j) - (1.0 / mRVE_Dimension));

              cauchy_trace       += mRVE_CauchyTensor(i,j);
              cauchy_trace_inner += mRVE_CauchyTensorInner(i,j);
            }
            else {
              deviatoric_fabric       = 4.0 * mRVE_FabricTensor(i,j);
              deviatoric_fabric_inner = 4.0 * mRVE_FabricTensor(i,j);
            }

            double_dot_product       += 0.5 * deviatoric_fabric       * deviatoric_fabric;
            double_dot_product_inner += 0.5 * deviatoric_fabric_inner * deviatoric_fabric_inner;
          }
        }
        mRVE_Anisotropy        = sqrt(double_dot_product);
        mRVE_AnisotropyInner   = sqrt(double_dot_product_inner);

        mRVE_EffectStress      = cauchy_trace       / mRVE_Dimension;
        mRVE_EffectStressInner = cauchy_trace_inner / mRVE_Dimension;

        // Second loop to compute deviatoric stress
        double_dot_product       = 0.0;
        double_dot_product_inner = 0.0;

        for (int i = 0; i < mRVE_Dimension; i++) {
          for (int j = 0; j < mRVE_Dimension; j++) {
            deviatoric_stress       = (i==j) ? mRVE_CauchyTensor(i,j)      - mRVE_EffectStress      : mRVE_CauchyTensor(i,j);
            deviatoric_stress_inner = (i==j) ? mRVE_CauchyTensorInner(i,j) - mRVE_EffectStressInner : mRVE_CauchyTensorInner(i,j);

            double_dot_product       += 0.5 * deviatoric_stress       * deviatoric_stress;
            double_dot_product_inner += 0.5 * deviatoric_stress_inner * deviatoric_stress_inner;
          }
        }
        mRVE_DevStress      = sqrt(double_dot_product);
        mRVE_DevStressInner = sqrt(double_dot_product_inner);
      }
    }

    //-----------------------------------------------------------------------------------------------------------------------------------------
    // TODO: Missing 3D
    void ExplicitSolverStrategy::RVEComputeRoseUniformity(void) {
      int len_rose_all = mRVE_RoseDiagram.size2();
      int len_rose_inn = mRVE_RoseDiagramInner.size2();

      // Summation
      double sum_rose_xy_all = 0.0;
      double sum_rose_xy_inn = 0.0;

      for (int i = 0; i < len_rose_all; i++)
        sum_rose_xy_all += mRVE_RoseDiagram(0,i);

      for (int i = 0; i < len_rose_inn; i++)
        sum_rose_xy_inn += mRVE_RoseDiagramInner(0,i);

      // Normalize
      std::vector<double> rose_xy_all_norm(len_rose_all, 0.0);
      std::vector<double> rose_xy_inn_norm(len_rose_inn, 0.0);
      
      if (sum_rose_xy_all != 0.0)
        for (int i = 0; i < len_rose_all; i++)
          rose_xy_all_norm[i] = mRVE_RoseDiagram(0,i) / sum_rose_xy_all;
      
      if (sum_rose_xy_inn != 0.0)
        for (int i = 0; i < len_rose_inn; i++)
          rose_xy_inn_norm[i] = mRVE_RoseDiagramInner(0,i) / sum_rose_xy_inn;

      // Mean
      double mean_rose_xy_all = 0.0;
      double mean_rose_xy_inn = 0.0;
      
      if (sum_rose_xy_all != 0.0)
        mean_rose_xy_all = sum_rose_xy_all / len_rose_all;

      if (sum_rose_xy_inn != 0.0)
        mean_rose_xy_inn = sum_rose_xy_inn / len_rose_inn;

      // Variance
      double var_rose_xy_all = 0.0;
      double var_rose_xy_inn = 0.0;
      
      for (int i = 0; i < len_rose_all; i++)
        var_rose_xy_all += pow(mRVE_RoseDiagram(0,i) - mean_rose_xy_all, 2.0);
      var_rose_xy_all /= len_rose_all;

      for (int i = 0; i < len_rose_inn; i++)
        var_rose_xy_inn += pow(mRVE_RoseDiagramInner(0,i) - mean_rose_xy_inn, 2.0);
      var_rose_xy_inn /= len_rose_inn;

      // Standard deviation
      mRVE_StdDevRoseXYAll = sqrt(var_rose_xy_all);
      mRVE_StdDevRoseXYInn = sqrt(var_rose_xy_inn);
      mRVE_StdDevRoseAzAll = 0.0;
      mRVE_StdDevRoseAzInn = 0.0;
    }

    //-----------------------------------------------------------------------------------------------------------------------------------------
    void ExplicitSolverStrategy::RVEStopCompression(void) {
      if (!mRVE_Compress)
        return;

      ModelPart& r_dem_model_part = GetModelPart();
      ProcessInfo& r_process_info = r_dem_model_part.GetProcessInfo();
      std::string criterion       = r_process_info[RVE_CONSOLIDATION_STOP_CRITERION];
      bool check = true;

      if (criterion.compare("time") == 0) {
        const double limit = r_process_info[LIMIT_CONSOLIDATION_TIME];
        const double time  = r_process_info[TIME];
        check = time >= limit;
      }
      else if (criterion.compare("stress") == 0) {
        const double limit = r_process_info[LIMIT_CONSOLIDATION_STRESS];
        check = std::abs(mRVE_EffectStressInner) >= limit;
      }
      else if (criterion.compare("porosity") == 0) {
        const double limit = r_process_info[LIMIT_CONSOLIDATION_POROSITY];
        const bool inner_porosity = r_process_info[INNER_CONSOLIDATION_POROSITY];
        if (inner_porosity)
          check = mRVE_PorosityInner <= limit;
        else
          check = mRVE_Porosity <= limit;
      }

      if (check) {
        mRVE_Compress = false;

        if (mRVE_FlatWalls) {
          ModelPart& fem_model_part = GetFemModelPart();
          ModelPart::ConditionsContainerType& r_conditions = fem_model_part.GetCommunicator().LocalMesh().Conditions();

          for (ModelPart::SubModelPartsContainerType::iterator sub_model_part = fem_model_part.SubModelPartsBegin(); sub_model_part != fem_model_part.SubModelPartsEnd(); ++sub_model_part) {
            ModelPart& submp = *sub_model_part;
            array_1d<double, 3>& linear_velocity = submp[LINEAR_VELOCITY];
            linear_velocity[0] = 0.0;
            linear_velocity[1] = 0.0;
            linear_velocity[2] = 0.0;
          }

          for (unsigned int i = 0; i < r_conditions.size(); i++) {
            ModelPart::ConditionsContainerType::iterator it = r_conditions.ptr_begin() + i;
            DEMWall* p_wall = dynamic_cast<DEMWall*> (&(*it));
            for (unsigned int inode = 0; inode < p_wall->GetGeometry().size(); inode++) {
              array_1d<double, 3>& wall_velocity = p_wall->GetGeometry()[inode].FastGetSolutionStepValue(VELOCITY);
              noalias(wall_velocity) = ZeroVector(3);
            }
          }
        }

        else {
          for (ModelPart::SubModelPartsContainerType::iterator sub_model_part = r_dem_model_part.SubModelPartsBegin(); sub_model_part != r_dem_model_part.SubModelPartsEnd(); ++sub_model_part) {
            ModelPart& submp = *sub_model_part;
            array_1d<double, 3>& linear_velocity = submp[LINEAR_VELOCITY];
            if (linear_velocity[0] != 0.0 || linear_velocity[1] != 0.0 || linear_velocity[2] != 0.0) {
              linear_velocity[0] = 0.0;
              linear_velocity[1] = 0.0;
              linear_velocity[2] = 0.0;
            }
          }
          for (int i = 0; i < (int)mListOfSphericParticles.size(); i++) {
            if (mListOfSphericParticles[i]->mWall)
              mListOfSphericParticles[i]->mMoving = false;
          }
        }
      }
    }

    //-----------------------------------------------------------------------------------------------------------------------------------------
    void ExplicitSolverStrategy::RVEReadOldForces(void) {
      std::fstream f_old_elastic_forces;
      std::fstream f_old_elastic_forces_walls;

      f_old_elastic_forces.open("OLD_FORCES.txt", std::ios::in);
      f_old_elastic_forces_walls.open("OLD_FORCES_WALLS.txt", std::ios::in);

      if (f_old_elastic_forces) {
        for (int i = 0; i < mListOfSphericParticles.size(); i++) {
          int particle_id, n_neighbors;
          f_old_elastic_forces >> particle_id >> n_neighbors;
          for (int j = 0; j < n_neighbors; j++) {
            double fx, fy, fz;
            f_old_elastic_forces >> fx >> fy >> fz;
            mListOfSphericParticles[i]->mNeighbourElasticContactForces[j][0] = fx;
            mListOfSphericParticles[i]->mNeighbourElasticContactForces[j][1] = fy;
            mListOfSphericParticles[i]->mNeighbourElasticContactForces[j][2] = fz;
          }
        }
        f_old_elastic_forces.close();
      }

      if (f_old_elastic_forces_walls) {
        for (int i = 0; i < mListOfSphericParticles.size(); i++) {
          int particle_id, n_neighbors;
          f_old_elastic_forces_walls >> particle_id >> n_neighbors;
          for (int j = 0; j < n_neighbors; j++) {
            double fx, fy, fz;
            f_old_elastic_forces_walls >> fx >> fy >> fz;
            mListOfSphericParticles[i]->mNeighbourRigidFacesElasticContactForce[j][0] = fx;
            mListOfSphericParticles[i]->mNeighbourRigidFacesElasticContactForce[j][1] = fy;
            mListOfSphericParticles[i]->mNeighbourRigidFacesElasticContactForce[j][2] = fz;
          }
        }
        f_old_elastic_forces_walls.close();
      }
    }

    //-----------------------------------------------------------------------------------------------------------------------------------------
    void ExplicitSolverStrategy::RVEWriteFiles(void) {
      ModelPart&   r_dem_model_part = GetModelPart();
      ProcessInfo& r_process_info   = r_dem_model_part.GetProcessInfo();
      const int    time_step        = r_process_info[TIME_STEPS];
      const double time             = r_process_info[TIME];

      const int eval_freq = r_process_info[RVE_WRITE_FREQ];
      if (time_step % eval_freq != 0.0)
        return;

      const int write_coords_freq = r_process_info[RVE_WRITE_COORDINATES_FREQ];
      if (mRVE_FileCoordinatesHistory.is_open() && time_step % write_coords_freq == 0.0) {
        const double x1 = mRVE_CornerCoordsX[0]; const double y1 = mRVE_CornerCoordsY[0];
        const double x2 = mRVE_CornerCoordsX[1]; const double y2 = mRVE_CornerCoordsY[1];
        const double x3 = mRVE_CornerCoordsX[2]; const double y3 = mRVE_CornerCoordsY[2];
        const double x4 = mRVE_CornerCoordsX[3]; const double y4 = mRVE_CornerCoordsY[3];

        mRVE_FileCoordinatesHistory << std::defaultfloat << time_step << " " << time << " ";
        mRVE_FileCoordinatesHistory << std::fixed << std::setprecision(12) << x1 << " " << y1 << " " << x2 << " " << y2 << " " << x3 << " " << y3 << " " << x4 << " " << y4 << " ";

        const int number_of_particles = (int)mListOfSphericParticles.size();
        for (int i = 0; i < number_of_particles; i++) {
          const double r = mListOfSphericParticles[i]->GetRadius();
          const double x = mListOfSphericParticles[i]->GetGeometry()[0][0];
          const double y = mListOfSphericParticles[i]->GetGeometry()[0][1];
          mRVE_FileCoordinatesHistory << std::fixed << std::setprecision(12) << r << " " << x << " " << y << " ";
        }
        mRVE_FileCoordinatesHistory << std::endl;
      }

      if (true) { //(r_process_info[POST_WRITE_COORDINATES_LAST]) {
        //RVEWriteCoords();
      }

      if (mRVE_FilePorosity.is_open()) {
        mRVE_FilePorosity << time_step                    << " "
                          << time                         << " "
                          << mRVE_VolInner                << " "
                          << mRVE_VolTotal                << " "
                          << mRVE_VolSolid                << " "
                          << mRVE_VolTotal-mRVE_VolSolid  << " "
                          << mRVE_Porosity                << " "
                          << mRVE_PorosityInner           << " "
                          << mRVE_VoidRatio               << " "
                          << mRVE_VoidRatioInner
                          << std::endl;
      }

      if (mRVE_FileContactNumber.is_open()) {
        mRVE_FileContactNumber << time_step << " " << time << " ";
        for (int i = 0; i < mListOfSphericParticles.size(); i++) {
          if (mListOfSphericParticles[i]->mWall == 0) {
            const int contacts = mListOfSphericParticles[i]->mCoordNum;
            mRVE_FileContactNumber << contacts << " ";
          }
        }
        mRVE_FileContactNumber << std::endl;
      }

      if (mRVE_FileCoordNumber.is_open()) {
        mRVE_FileCoordNumber << time << " " << mRVE_AvgCoordNum << std::endl;
      }

      if (mRVE_FileInnerVolumeParticles.is_open()) {
        mRVE_FileInnerVolumeParticles << time_step << " " << time << " ";
        mRVE_FileInnerVolumeParticles << mRVE_InnerVolParticles.size() << " ";
        for (int i = 0; i < mRVE_InnerVolParticles.size(); i++) {
          array_1d<double,3> coords = mRVE_InnerVolParticles[i]->GetGeometry()[0].Coordinates();
          const double radius       = mRVE_InnerVolParticles[i]->GetRadius();
          mRVE_FileInnerVolumeParticles << coords[0] << " " << coords[1] << " " << coords[2] << " " << radius << " ";
        }
        mRVE_FileInnerVolumeParticles << std::endl;
      }

      if (mRVE_FileForceChain.is_open()) {
        mRVE_FileForceChain << time_step << " " << time << " ";
        for (int i = 0; i < mRVE_ForceChain.size(); i++) mRVE_FileForceChain << mRVE_ForceChain[i] << " ";
        mRVE_FileForceChain << std::endl;
      }

      if (true) { //(mRVE_FileElasticContactForces.is_open()) {
        //RVEWriteForceParticles();
      }

      if (true) {
        //RVEWriteForceContacts();
      }

      if (mRVE_FileRoseDiagram.is_open()) {
        mRVE_FileRoseDiagram << time_step << " " << time << " ";
        
        mRVE_FileRoseDiagram << "[ ";
        for (unsigned int i = 0; i < mRVE_RoseDiagram.size2(); i++) mRVE_FileRoseDiagram << mRVE_RoseDiagram(0,i) << " ";
        mRVE_FileRoseDiagram << "] ";

        mRVE_FileRoseDiagram << "[ ";
        for (unsigned int i = 0; i < mRVE_RoseDiagram.size2(); i++) mRVE_FileRoseDiagram << mRVE_RoseDiagram(1,i) << " ";
        mRVE_FileRoseDiagram << "]";

        mRVE_FileRoseDiagram << std::endl;
      }

      if (mRVE_FileRoseDiagramInner.is_open()) {
        mRVE_FileRoseDiagramInner << time_step << " " << time << " ";

        mRVE_FileRoseDiagramInner << "[ ";
        for (unsigned int i = 0; i < mRVE_RoseDiagramInner.size2(); i++) mRVE_FileRoseDiagramInner << mRVE_RoseDiagramInner(0,i) << " ";
        mRVE_FileRoseDiagramInner << "] ";

        mRVE_FileRoseDiagramInner << "[ ";
        for (unsigned int i = 0; i < mRVE_RoseDiagramInner.size2(); i++) mRVE_FileRoseDiagramInner << mRVE_RoseDiagramInner(1,i) << " ";
        mRVE_FileRoseDiagramInner << "]";

        mRVE_FileRoseDiagramInner << std::endl;
      }

      if (mRVE_FileRoseDiagramUniformity.is_open()) {
        mRVE_FileRoseDiagramUniformity << time_step << " "
                                       << time      << " "
                                       << mRVE_StdDevRoseXYAll << " "
                                       << mRVE_StdDevRoseXYInn
                                       << std::endl;
      }

      if (mRVE_FileAnisotropy.is_open()) {
        mRVE_FileAnisotropy << time_step       << " "
                            << time            << " "
                            << mRVE_Anisotropy << " "
                            << mRVE_AnisotropyInner
                            << std::endl;
      }

      if (mRVE_FileFabricTensor.is_open()) {
        if (mRVE_Dimension == 2)
          mRVE_FileFabricTensor << time_step << " " << time << " "
                                << "[[" << mRVE_FabricTensor(0,0) << "],[" << mRVE_FabricTensor(0,1) << "]]" << " "
                                << "[[" << mRVE_FabricTensor(1,0) << "],[" << mRVE_FabricTensor(1,1) << "]]"
                                << std::endl;
        else if (mRVE_Dimension == 3)
          mRVE_FileFabricTensor << time_step << " " << time << " "
                                << "[[" << mRVE_FabricTensor(0,0) << "],[" << mRVE_FabricTensor(0,1) << "],[" << mRVE_FabricTensor(0,2) << "]]" << " "
                                << "[[" << mRVE_FabricTensor(1,0) << "],[" << mRVE_FabricTensor(1,1) << "],[" << mRVE_FabricTensor(1,2) << "]]" << " "
                                << "[[" << mRVE_FabricTensor(2,0) << "],[" << mRVE_FabricTensor(2,1) << "],[" << mRVE_FabricTensor(2,2) << "]]"
                                << std::endl;
      }

      if (mRVE_FileFabricTensorInner.is_open()) {
        if (mRVE_Dimension == 2)
          mRVE_FileFabricTensorInner << time_step << " " << time << " "
                                     << "[[" << mRVE_FabricTensorInner(0,0) << "],[" << mRVE_FabricTensorInner(0,1) << "]]" << " "
                                     << "[[" << mRVE_FabricTensorInner(1,0) << "],[" << mRVE_FabricTensorInner(1,1) << "]]"
                                     << std::endl;
        else if (mRVE_Dimension == 3)
          mRVE_FileFabricTensorInner << time_step << " " << time << " "
                                     << "[[" << mRVE_FabricTensorInner(0,0) << "],[" << mRVE_FabricTensorInner(0,1) << "],[" << mRVE_FabricTensorInner(0,2) << "]]" << " "
                                     << "[[" << mRVE_FabricTensorInner(1,0) << "],[" << mRVE_FabricTensorInner(1,1) << "],[" << mRVE_FabricTensorInner(1,2) << "]]" << " "
                                     << "[[" << mRVE_FabricTensorInner(2,0) << "],[" << mRVE_FabricTensorInner(2,1) << "],[" << mRVE_FabricTensorInner(2,2) << "]]"
                                     << std::endl;
      }

      if (mRVE_FileStress.is_open()) {
        mRVE_FileStress << time_step              << " "
                        << time                   << " "
                        << mRVE_WallStress        << " "
                        << mRVE_EffectStress      << " "
                        << mRVE_DevStress         << " "
                        << mRVE_EffectStressInner << " "
                        << mRVE_DevStressInner
                        << std::endl;
      }

      if (mRVE_FileCauchyTensor.is_open()) {
        if (mRVE_Dimension == 2)
          mRVE_FileCauchyTensor << time_step << " " << time << " "
                                << "[[" << mRVE_CauchyTensor(0,0) << "],[" << mRVE_CauchyTensor(0,1) << "]]" << " "
                                << "[[" << mRVE_CauchyTensor(1,0) << "],[" << mRVE_CauchyTensor(1,1) << "]]"
                                << std::endl;
        else if (mRVE_Dimension == 3)
          mRVE_FileCauchyTensor << time_step << " " << time << " "
                                << "[[" << mRVE_CauchyTensor(0,0) << "],[" << mRVE_CauchyTensor(0,1) << "],[" << mRVE_CauchyTensor(0,2) << "]]" << " "
                                << "[[" << mRVE_CauchyTensor(1,0) << "],[" << mRVE_CauchyTensor(1,1) << "],[" << mRVE_CauchyTensor(1,2) << "]]" << " "
                                << "[[" << mRVE_CauchyTensor(2,0) << "],[" << mRVE_CauchyTensor(2,1) << "],[" << mRVE_CauchyTensor(2,2) << "]]"
                                << std::endl;
      }

      if (mRVE_FileCauchyTensorInner.is_open()) {
        if (mRVE_Dimension == 2)
          mRVE_FileCauchyTensorInner << time_step << " " << time << " "
                                     << "[[" << mRVE_CauchyTensorInner(0,0) << "],[" << mRVE_CauchyTensorInner(0,1) << "]]" << " "
                                     << "[[" << mRVE_CauchyTensorInner(1,0) << "],[" << mRVE_CauchyTensorInner(1,1) << "]]"
                                     << std::endl;
        else if (mRVE_Dimension == 3)
          mRVE_FileCauchyTensorInner << time_step << " " << time << " "
                                     << "[[" << mRVE_CauchyTensorInner(0,0) << "],[" << mRVE_CauchyTensorInner(0,1) << "],[" << mRVE_CauchyTensorInner(0,2) << "]]" << " "
                                     << "[[" << mRVE_CauchyTensorInner(1,0) << "],[" << mRVE_CauchyTensorInner(1,1) << "],[" << mRVE_CauchyTensorInner(1,2) << "]]" << " "
                                     << "[[" << mRVE_CauchyTensorInner(2,0) << "],[" << mRVE_CauchyTensorInner(2,1) << "],[" << mRVE_CauchyTensorInner(2,2) << "]]"
                                     << std::endl;
      }

      if (mRVE_FileTangentTensor.is_open()) {
          if (mRVE_Dimension == 2)
            mRVE_FileTangentTensor << time_step << " " << time << " "
                                   << "[[" << mRVE_TangentTensor(0,0) << "],[" << mRVE_TangentTensor(0,1) << "],[" << mRVE_TangentTensor(0,2) << "],[" << mRVE_TangentTensor(0,3) << "]]" << " "
                                   << "[[" << mRVE_TangentTensor(1,0) << "],[" << mRVE_TangentTensor(1,1) << "],[" << mRVE_TangentTensor(1,2) << "],[" << mRVE_TangentTensor(1,3) << "]]" << " "
                                   << "[[" << mRVE_TangentTensor(2,0) << "],[" << mRVE_TangentTensor(2,1) << "],[" << mRVE_TangentTensor(2,2) << "],[" << mRVE_TangentTensor(2,3) << "]]" << " "
                                   << "[[" << mRVE_TangentTensor(3,0) << "],[" << mRVE_TangentTensor(3,1) << "],[" << mRVE_TangentTensor(3,2) << "],[" << mRVE_TangentTensor(3,3) << "]]"
                                   << std::endl;
          else if (mRVE_Dimension == 3)
            mRVE_FileTangentTensor << time_step << " " << time << " "
                                   << "[[" << mRVE_TangentTensor(0,0) << "],[" << mRVE_TangentTensor(0,1) << "],[" << mRVE_TangentTensor(0,2) << "],[" << mRVE_TangentTensor(0,3) << "],[" << mRVE_TangentTensor(0,4) << "],[" << mRVE_TangentTensor(0,5) << "],[" << mRVE_TangentTensor(0,6) << "],[" << mRVE_TangentTensor(0,7) << "],[" << mRVE_TangentTensor(0,8) << "]]" << " "
                                   << "[[" << mRVE_TangentTensor(1,0) << "],[" << mRVE_TangentTensor(1,1) << "],[" << mRVE_TangentTensor(1,2) << "],[" << mRVE_TangentTensor(1,3) << "],[" << mRVE_TangentTensor(1,4) << "],[" << mRVE_TangentTensor(1,5) << "],[" << mRVE_TangentTensor(1,6) << "],[" << mRVE_TangentTensor(1,7) << "],[" << mRVE_TangentTensor(1,8) << "]]" << " "
                                   << "[[" << mRVE_TangentTensor(2,0) << "],[" << mRVE_TangentTensor(2,1) << "],[" << mRVE_TangentTensor(2,2) << "],[" << mRVE_TangentTensor(2,3) << "],[" << mRVE_TangentTensor(2,4) << "],[" << mRVE_TangentTensor(2,5) << "],[" << mRVE_TangentTensor(2,6) << "],[" << mRVE_TangentTensor(2,7) << "],[" << mRVE_TangentTensor(2,8) << "]]" << " "
                                   << "[[" << mRVE_TangentTensor(3,0) << "],[" << mRVE_TangentTensor(3,1) << "],[" << mRVE_TangentTensor(3,2) << "],[" << mRVE_TangentTensor(3,3) << "],[" << mRVE_TangentTensor(3,4) << "],[" << mRVE_TangentTensor(3,5) << "],[" << mRVE_TangentTensor(3,6) << "],[" << mRVE_TangentTensor(3,7) << "],[" << mRVE_TangentTensor(3,8) << "]]" << " "
                                   << "[[" << mRVE_TangentTensor(4,0) << "],[" << mRVE_TangentTensor(4,1) << "],[" << mRVE_TangentTensor(4,2) << "],[" << mRVE_TangentTensor(4,3) << "],[" << mRVE_TangentTensor(4,4) << "],[" << mRVE_TangentTensor(4,5) << "],[" << mRVE_TangentTensor(4,6) << "],[" << mRVE_TangentTensor(4,7) << "],[" << mRVE_TangentTensor(4,8) << "]]" << " "
                                   << "[[" << mRVE_TangentTensor(5,0) << "],[" << mRVE_TangentTensor(5,1) << "],[" << mRVE_TangentTensor(5,2) << "],[" << mRVE_TangentTensor(5,3) << "],[" << mRVE_TangentTensor(5,4) << "],[" << mRVE_TangentTensor(5,5) << "],[" << mRVE_TangentTensor(5,6) << "],[" << mRVE_TangentTensor(5,7) << "],[" << mRVE_TangentTensor(5,8) << "]]" << " "
                                   << "[[" << mRVE_TangentTensor(6,0) << "],[" << mRVE_TangentTensor(6,1) << "],[" << mRVE_TangentTensor(6,2) << "],[" << mRVE_TangentTensor(6,3) << "],[" << mRVE_TangentTensor(6,4) << "],[" << mRVE_TangentTensor(6,5) << "],[" << mRVE_TangentTensor(6,6) << "],[" << mRVE_TangentTensor(6,7) << "],[" << mRVE_TangentTensor(6,8) << "]]" << " "
                                   << "[[" << mRVE_TangentTensor(7,0) << "],[" << mRVE_TangentTensor(7,1) << "],[" << mRVE_TangentTensor(7,2) << "],[" << mRVE_TangentTensor(7,3) << "],[" << mRVE_TangentTensor(7,4) << "],[" << mRVE_TangentTensor(7,5) << "],[" << mRVE_TangentTensor(7,6) << "],[" << mRVE_TangentTensor(7,7) << "],[" << mRVE_TangentTensor(7,8) << "]]" << " "
                                   << "[[" << mRVE_TangentTensor(8,0) << "],[" << mRVE_TangentTensor(8,1) << "],[" << mRVE_TangentTensor(8,2) << "],[" << mRVE_TangentTensor(8,3) << "],[" << mRVE_TangentTensor(8,4) << "],[" << mRVE_TangentTensor(8,5) << "],[" << mRVE_TangentTensor(8,6) << "],[" << mRVE_TangentTensor(8,7) << "],[" << mRVE_TangentTensor(8,8) << "]]"
                                   << std::endl;
      }

      if (mRVE_FileTangentTensorInner.is_open()) {
          if (mRVE_Dimension == 2)
            mRVE_FileTangentTensorInner << time_step << " " << time << " "
                                        << "[[" << mRVE_TangentTensorInner(0,0) << "],[" << mRVE_TangentTensorInner(0,1) << "],[" << mRVE_TangentTensorInner(0,2) << "],[" << mRVE_TangentTensorInner(0,3) << "]]" << " "
                                        << "[[" << mRVE_TangentTensorInner(1,0) << "],[" << mRVE_TangentTensorInner(1,1) << "],[" << mRVE_TangentTensorInner(1,2) << "],[" << mRVE_TangentTensorInner(1,3) << "]]" << " "
                                        << "[[" << mRVE_TangentTensorInner(2,0) << "],[" << mRVE_TangentTensorInner(2,1) << "],[" << mRVE_TangentTensorInner(2,2) << "],[" << mRVE_TangentTensorInner(2,3) << "]]" << " "
                                        << "[[" << mRVE_TangentTensorInner(3,0) << "],[" << mRVE_TangentTensorInner(3,1) << "],[" << mRVE_TangentTensorInner(3,2) << "],[" << mRVE_TangentTensorInner(3,3) << "]]"
                                        << std::endl;
          else if (mRVE_Dimension == 3)
            mRVE_FileTangentTensorInner << time_step << " " << time << " "
                                        << "[[" << mRVE_TangentTensorInner(0,0) << "],[" << mRVE_TangentTensorInner(0,1) << "],[" << mRVE_TangentTensorInner(0,2) << "],[" << mRVE_TangentTensorInner(0,3) << "],[" << mRVE_TangentTensorInner(0,4) << "],[" << mRVE_TangentTensorInner(0,5) << "],[" << mRVE_TangentTensorInner(0,6) << "],[" << mRVE_TangentTensorInner(0,7) << "],[" << mRVE_TangentTensorInner(0,8) << "]]" << " "
                                        << "[[" << mRVE_TangentTensorInner(1,0) << "],[" << mRVE_TangentTensorInner(1,1) << "],[" << mRVE_TangentTensorInner(1,2) << "],[" << mRVE_TangentTensorInner(1,3) << "],[" << mRVE_TangentTensorInner(1,4) << "],[" << mRVE_TangentTensorInner(1,5) << "],[" << mRVE_TangentTensorInner(1,6) << "],[" << mRVE_TangentTensorInner(1,7) << "],[" << mRVE_TangentTensorInner(1,8) << "]]" << " "
                                        << "[[" << mRVE_TangentTensorInner(2,0) << "],[" << mRVE_TangentTensorInner(2,1) << "],[" << mRVE_TangentTensorInner(2,2) << "],[" << mRVE_TangentTensorInner(2,3) << "],[" << mRVE_TangentTensorInner(2,4) << "],[" << mRVE_TangentTensorInner(2,5) << "],[" << mRVE_TangentTensorInner(2,6) << "],[" << mRVE_TangentTensorInner(2,7) << "],[" << mRVE_TangentTensorInner(2,8) << "]]" << " "
                                        << "[[" << mRVE_TangentTensorInner(3,0) << "],[" << mRVE_TangentTensorInner(3,1) << "],[" << mRVE_TangentTensorInner(3,2) << "],[" << mRVE_TangentTensorInner(3,3) << "],[" << mRVE_TangentTensorInner(3,4) << "],[" << mRVE_TangentTensorInner(3,5) << "],[" << mRVE_TangentTensorInner(3,6) << "],[" << mRVE_TangentTensorInner(3,7) << "],[" << mRVE_TangentTensorInner(3,8) << "]]" << " "
                                        << "[[" << mRVE_TangentTensorInner(4,0) << "],[" << mRVE_TangentTensorInner(4,1) << "],[" << mRVE_TangentTensorInner(4,2) << "],[" << mRVE_TangentTensorInner(4,3) << "],[" << mRVE_TangentTensorInner(4,4) << "],[" << mRVE_TangentTensorInner(4,5) << "],[" << mRVE_TangentTensorInner(4,6) << "],[" << mRVE_TangentTensorInner(4,7) << "],[" << mRVE_TangentTensorInner(4,8) << "]]" << " "
                                        << "[[" << mRVE_TangentTensorInner(5,0) << "],[" << mRVE_TangentTensorInner(5,1) << "],[" << mRVE_TangentTensorInner(5,2) << "],[" << mRVE_TangentTensorInner(5,3) << "],[" << mRVE_TangentTensorInner(5,4) << "],[" << mRVE_TangentTensorInner(5,5) << "],[" << mRVE_TangentTensorInner(5,6) << "],[" << mRVE_TangentTensorInner(5,7) << "],[" << mRVE_TangentTensorInner(5,8) << "]]" << " "
                                        << "[[" << mRVE_TangentTensorInner(6,0) << "],[" << mRVE_TangentTensorInner(6,1) << "],[" << mRVE_TangentTensorInner(6,2) << "],[" << mRVE_TangentTensorInner(6,3) << "],[" << mRVE_TangentTensorInner(6,4) << "],[" << mRVE_TangentTensorInner(6,5) << "],[" << mRVE_TangentTensorInner(6,6) << "],[" << mRVE_TangentTensorInner(6,7) << "],[" << mRVE_TangentTensorInner(6,8) << "]]" << " "
                                        << "[[" << mRVE_TangentTensorInner(7,0) << "],[" << mRVE_TangentTensorInner(7,1) << "],[" << mRVE_TangentTensorInner(7,2) << "],[" << mRVE_TangentTensorInner(7,3) << "],[" << mRVE_TangentTensorInner(7,4) << "],[" << mRVE_TangentTensorInner(7,5) << "],[" << mRVE_TangentTensorInner(7,6) << "],[" << mRVE_TangentTensorInner(7,7) << "],[" << mRVE_TangentTensorInner(7,8) << "]]" << " "
                                        << "[[" << mRVE_TangentTensorInner(8,0) << "],[" << mRVE_TangentTensorInner(8,1) << "],[" << mRVE_TangentTensorInner(8,2) << "],[" << mRVE_TangentTensorInner(8,3) << "],[" << mRVE_TangentTensorInner(8,4) << "],[" << mRVE_TangentTensorInner(8,5) << "],[" << mRVE_TangentTensorInner(8,6) << "],[" << mRVE_TangentTensorInner(8,7) << "],[" << mRVE_TangentTensorInner(8,8) << "]]"
                                        << std::endl;
      }

      if (mRVE_FileConductivityTensor.is_open()) {
        if (mRVE_Dimension == 2)
          mRVE_FileConductivityTensor << time_step << " " << time << " "
                                      << "[[" << mRVE_ConductivityTensor(0,0) << "],[" << mRVE_ConductivityTensor(0,1) << "]]" << " "
                                      << "[[" << mRVE_ConductivityTensor(1,0) << "],[" << mRVE_ConductivityTensor(1,1) << "]]"
                                      << std::endl;
        else if (mRVE_Dimension == 3)
          mRVE_FileConductivityTensor << time_step << " " << time << " "
                                      << "[[" << mRVE_ConductivityTensor(0,0) << "],[" << mRVE_ConductivityTensor(0,1) << "],[" << mRVE_ConductivityTensor(0,2) << "]]" << " "
                                      << "[[" << mRVE_ConductivityTensor(1,0) << "],[" << mRVE_ConductivityTensor(1,1) << "],[" << mRVE_ConductivityTensor(1,2) << "]]" << " "
                                      << "[[" << mRVE_ConductivityTensor(2,0) << "],[" << mRVE_ConductivityTensor(2,1) << "],[" << mRVE_ConductivityTensor(2,2) << "]]"
                                      << std::endl;
      }

      if (mRVE_FileConductivityTensorInner.is_open()) {
        if (mRVE_Dimension == 2)
          mRVE_FileConductivityTensorInner << time_step << " " << time << " "
                                           << "[[" << mRVE_ConductivityTensorInner(0,0) << "],[" << mRVE_ConductivityTensorInner(0,1) << "]]" << " "
                                           << "[[" << mRVE_ConductivityTensorInner(1,0) << "],[" << mRVE_ConductivityTensorInner(1,1) << "]]"
                                           << std::endl;
        else if (mRVE_Dimension == 3)
          mRVE_FileConductivityTensorInner << time_step << " " << time << " "
                                           << "[[" << mRVE_ConductivityTensorInner(0,0) << "],[" << mRVE_ConductivityTensorInner(0,1) << "],[" << mRVE_ConductivityTensorInner(0,2) << "]]" << " "
                                           << "[[" << mRVE_ConductivityTensorInner(1,0) << "],[" << mRVE_ConductivityTensorInner(1,1) << "],[" << mRVE_ConductivityTensorInner(1,2) << "]]" << " "
                                           << "[[" << mRVE_ConductivityTensorInner(2,0) << "],[" << mRVE_ConductivityTensorInner(2,1) << "],[" << mRVE_ConductivityTensorInner(2,2) << "]]"
                                           << std::endl;
      }

      if (mRVE_FileFKS.is_open()) {
        mRVE_FileFKS << time_step << " " << time << " "
                     << mRVE_FabricTensorInner(1,0) << " "
                     << mRVE_FabricTensorInner(0,0)-mRVE_FabricTensorInner(1,1) << " "
                     << mRVE_ConductivityTensorInner(0,0) << " "
                     << mRVE_ConductivityTensorInner(1,1) << " "
                     << mRVE_ConductivityTensorInner(1,0)
                     << std::endl;
      }
    }

    //-----------------------------------------------------------------------------------------------------------------------------------------
    void ExplicitSolverStrategy::RVEWriteCoords(void) {
      ModelPart& r_dem_model_part = GetModelPart();
      ProcessInfo& r_process_info = r_dem_model_part.GetProcessInfo();
      const int time_step  = r_process_info[TIME_STEPS];
      const int write_freq = 1000;
      if (time_step % write_freq != 0.0)
        return;

      mRVE_FileCoordinatesLast.open("rve_coordinates_last.txt", std::ios::out | std::ios::trunc);
      mRVE_FileCoordinatesLast << "LINE 1: WALL_X_LOW_LEFT WALL_Y_LOW_LEFT WALL_X_LOW_RIGHT WALL_Y_LOW_RIGHT WALL_X_UP_RIGHT WALL_Y_UP_RIGHT WALL_X_UP_LEFT WALL_Y_UP_LEFT | ";
      mRVE_FileCoordinatesLast << "NEXT LINES: R X Y";
      mRVE_FileCoordinatesLast << std::endl;

      if (mRVE_CornerCoordsX.empty() || mRVE_CornerCoordsY.empty()) {
        mRVE_FileCoordinatesLast << std::fixed << std::setprecision(16) << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << std::endl;
      }
      else {
        const double x1 = mRVE_CornerCoordsX[0]; const double y1 = mRVE_CornerCoordsY[0];
        const double x2 = mRVE_CornerCoordsX[1]; const double y2 = mRVE_CornerCoordsY[1];
        const double x3 = mRVE_CornerCoordsX[2]; const double y3 = mRVE_CornerCoordsY[2];
        const double x4 = mRVE_CornerCoordsX[3]; const double y4 = mRVE_CornerCoordsY[3];
        mRVE_FileCoordinatesLast << std::fixed << std::setprecision(16) << x1 << " " << y1 << " " << x2 << " " << y2 << " " << x3 << " " << y3 << " " << x4 << " " << y4 << std::endl;
      }

      const int number_of_particles = (int)mListOfSphericParticles.size();
      for (int i = 0; i < number_of_particles; i++) {
        const double r = mListOfSphericParticles[i]->GetRadius();
        const double x = mListOfSphericParticles[i]->GetGeometry()[0][0];
        const double y = mListOfSphericParticles[i]->GetGeometry()[0][1];
        mRVE_FileCoordinatesLast << std::fixed << std::setprecision(16) << r << " " << x << " " << y << std::endl;
      }
      mRVE_FileCoordinatesLast.close();
    }

    //-----------------------------------------------------------------------------------------------------------------------------------------
    void ExplicitSolverStrategy::RVEWriteForceParticles(void) {
      ModelPart& r_dem_model_part = GetModelPart();
      ProcessInfo& r_process_info = r_dem_model_part.GetProcessInfo();
      const int time_step  = r_process_info[TIME_STEPS];
      const int write_freq = 1000;
      if (time_step % write_freq != 0.0)
        return;

      mRVE_FileElasticContactForces.open("rve_force_particles.txt", std::ios::out | std::ios::trunc);
      mRVE_FileElasticContactForcesWalls.open("rve_force_particles_walls.txt", std::ios::out | std::ios::trunc);

      for (int i = 0; i < mListOfSphericParticles.size(); i++) {
        const int id = mListOfSphericParticles[i]->Id();
        const int n_neighbors_p = mListOfSphericParticles[i]->mNeighbourElements.size();
        const int n_neighbors_w = mListOfSphericParticles[i]->mNeighbourRigidFaces.size();

        mRVE_FileElasticContactForces << std::defaultfloat << id << " " << n_neighbors_p << " ";
        for (int j = 0; j < n_neighbors_p; j++) {
          array_1d<double, 3> force = mListOfSphericParticles[i]->mNeighbourElasticContactForces[j];
          mRVE_FileElasticContactForces << std::fixed << std::setprecision(16) << force[0] << " " << force[1] << " " << force[2] << " ";
        }
        mRVE_FileElasticContactForces << std::endl;

        mRVE_FileElasticContactForcesWalls << std::defaultfloat << id << " " << n_neighbors_w << " ";
        for (int j = 0; j < n_neighbors_w; j++) {
          array_1d<double, 3> force = mListOfSphericParticles[i]->mNeighbourRigidFacesElasticContactForce[j];
          mRVE_FileElasticContactForcesWalls << std::fixed << std::setprecision(16) << force[0] << " " << force[1] << " " << force[2] << " ";
        }
        mRVE_FileElasticContactForcesWalls << std::endl;
      }

      mRVE_FileElasticContactForces.close();
      mRVE_FileElasticContactForcesWalls.close();
    }

    //-----------------------------------------------------------------------------------------------------------------------------------------
    void ExplicitSolverStrategy::RVEWriteForceContacts(void) {
      ModelPart& r_dem_model_part = GetModelPart();
      ProcessInfo& r_process_info = r_dem_model_part.GetProcessInfo();
      const double time    = r_process_info[TIME];
      const int time_step  = r_process_info[TIME_STEPS];
      const int write_freq = 1000;
      if (time_step % write_freq != 0.0)
        return;

      std::ofstream file;
      file.open("rve_force_chain.txt", std::ios::out | std::ios::trunc);

      file << "ID1 R1 X1 Y1 ID2 R2 X2 Y2 F1 F2 (WALLS: ID2=0 R2=0 X2,Y2=CONTACT_POINT)" << std::endl;
      file << "Time: ";
      file << std::defaultfloat << time << std::endl;
      
      for (int i = 0; i < mListOfSphericParticles.size(); i++) {
        std::size_t id1 = mListOfSphericParticles[i]->Id();
        const double r1 = mListOfSphericParticles[i]->GetRadius();
        const double x1 = mListOfSphericParticles[i]->GetGeometry()[0][0];
        const double y1 = mListOfSphericParticles[i]->GetGeometry()[0][1];

        // Particle-particle
        for (int j = 0; j < mListOfSphericParticles[i]->mNeighbourElements.size(); j++) {
          std::size_t id2 = mListOfSphericParticles[i]->mNeighbourElements[j]->Id();
          if (id2 <= id1) continue;

          const double r2 = mListOfSphericParticles[i]->mNeighbourElements[j]->GetRadius();
          const double x2 = mListOfSphericParticles[i]->mNeighbourElements[j]->GetGeometry()[0][0];
          const double y2 = mListOfSphericParticles[i]->mNeighbourElements[j]->GetGeometry()[0][1];
          array_1d<double, 3> f = mListOfSphericParticles[i]->mNeighbourElasticContactForces[j];

          file << std::defaultfloat << id1 << " ";
          file << std::fixed << std::setprecision(16) << r1 << " " << x1 << " " << y1 << " ";
          file << std::defaultfloat << id2 << " ";
          file << std::fixed << std::setprecision(16) << r2 << " " << x2 << " " << y2 << " ";
          file << std::fixed << std::setprecision(16) << f[0] << " " << f[1] << std::endl;
        }

        // Particle-wall
        for (int j = 0; j < mListOfSphericParticles[i]->mNeighbourRigidFaces.size(); j++) {
          DEMWall* p_wall = mListOfSphericParticles[i]->mNeighbourRigidFaces[j];
          Condition::GeometryType& geom = p_wall->GetGeometry();
          const double xw1 = geom[0][0];
          const double yw1 = geom[0][1];
          const double xw2 = geom[1][0];
          const double yw2 = geom[1][1];
          double x2, y2;

          if (std::abs(xw2-xw1) > std::abs(yw2-yw1)) {
            x2 = x1;
            y2 = yw1;
          }
          else {
            x2 = xw1;
            y2 = y1;
          }
          array_1d<double, 3> f = mListOfSphericParticles[i]->mNeighbourRigidFacesTotalContactForce[j];

          file << std::defaultfloat << id1 << " ";
          file << std::fixed << std::setprecision(16) << r1 << " " << x1 << " " << y1 << " ";
          file << std::defaultfloat << 0 << " ";
          file << std::fixed << std::setprecision(16) << 0.0 << " " << x2 << " " << y2 << " ";
          file << std::fixed << std::setprecision(16) << f[0] << " " << f[1] << std::endl;
        }
      }
      file.close();
    }

    //-----------------------------------------------------------------------------------------------------------------------------------------
    void ExplicitSolverStrategy::RVEOpenFiles(void) {
      ModelPart&   r_dem_model_part = GetModelPart();
      ProcessInfo& r_process_info   = r_dem_model_part.GetProcessInfo();
      
      //if (r_process_info[POST_WRITE_COORDINATES_HISTORY]) {
      //  mRVE_FileCoordinatesHistory.open("rve_coordinates_history.txt", std::ios::out);
      //  KRATOS_ERROR_IF_NOT(mRVE_FileCoordinatesHistory) << "Could not open file rve_coordinates_history.txt!" << std::endl;
      //  mRVE_FileCoordinatesHistory << "1 - STEP | ";
      //  mRVE_FileCoordinatesHistory << "2 - TIME | ";
      //  mRVE_FileCoordinatesHistory << "3 - WALL_X_LOW_LEFT WALL_Y_LOW_LEFT WALL_X_LOW_RIGHT WALL_Y_LOW_RIGHT WALL_X_UP_RIGHT WALL_Y_UP_RIGHT WALL_X_UP_LEFT WALL_Y_UP_LEFT | ";
      //  mRVE_FileCoordinatesHistory << "4 - R X Y of all particles";
      //  mRVE_FileCoordinatesHistory << std::endl;
      //}

      //if (r_process_info[POST_WRITE_POROSITY]) {
      //  mRVE_FilePorosity.open("rve_porosity.txt", std::ios::out);
      //  KRATOS_ERROR_IF_NOT(mRVE_FilePorosity) << "Could not open file rve_porosity.txt!" << std::endl;
      //  mRVE_FilePorosity << "1 - STEP | ";
      //  mRVE_FilePorosity << "2 - TIME | ";
      //  mRVE_FilePorosity << "3 - INNER VOLUME | ";
      //  mRVE_FilePorosity << "4 - TOTAL VOLUME | ";
      //  mRVE_FilePorosity << "5 - TOTAL SOLID VOLUME | ";
      //  mRVE_FilePorosity << "6 - TOTAL VOID VOLUME | ";
      //  mRVE_FilePorosity << "7 - POROSITY | ";
      //  mRVE_FilePorosity << "8 - POROSITY INNER | ";
      //  mRVE_FilePorosity << "9 - VOID RATIO | ";
      //  mRVE_FilePorosity << "10 - VOID RATIO INNER";
      //  mRVE_FilePorosity << std::endl;
      //}

      //if (r_process_info[POST_WRITE_CONTACT_NUMBER]) {
      //  mRVE_FileContactNumber.open("rve_contact_number.txt", std::ios::out);
      //  KRATOS_ERROR_IF_NOT(mRVE_FileContactNumber) << "Could not open file rve_conact_number.txt!" << std::endl;
      //  mRVE_FileContactNumber << "1 - STEP | ";
      //  mRVE_FileContactNumber << "2 - TIME | ";
      //  mRVE_FileContactNumber << "3 - NUMBER OF CONTACTS OF ALL PARTICLES";
      //  mRVE_FileContactNumber << std::endl;
      //}

      if (true) {
        mRVE_FileCoordNumber.open("rve_coordination_number.txt", std::ios::out);
        KRATOS_ERROR_IF_NOT(mRVE_FileCoordNumber) << "Could not open file rve_coordination_number.txt!" << std::endl;
        mRVE_FileCoordNumber << "1 - TIME | ";
        mRVE_FileCoordNumber << "2 - AVG COORDINATION NUMBER - ALL";
        mRVE_FileCoordNumber << std::endl;
      }

      //if (r_process_info[POST_WRITE_INNER_VOLUME_PARTICLES]) {
      //  mRVE_FileInnerVolumeParticles.open("rve_inner_volume_particles.txt", std::ios::out);
      //  KRATOS_ERROR_IF_NOT(mRVE_FileInnerVolumeParticles) << "Could not open file rve_inner_volume_particles.txt!" << std::endl;
      //  mRVE_FileInnerVolumeParticles << "1 - STEP | ";
      //  mRVE_FileInnerVolumeParticles << "2 - TIME | ";
      //  mRVE_FileInnerVolumeParticles << "3 - Number of particles | ";
      //  mRVE_FileInnerVolumeParticles << "4 - [X Y Z R] of each particles";
      //  mRVE_FileInnerVolumeParticles << std::endl;
      //}

      //if (r_process_info[POST_WRITE_FORCE_CHAIN]) {
      //  mRVE_FileForceChain.open("rve_force_chain.txt", std::ios::out);
      //  KRATOS_ERROR_IF_NOT(mRVE_FileForceChain) << "Could not open file rve_force_chain.txt!" << std::endl;
      //  mRVE_FileForceChain << "1 - STEP | ";
      //  mRVE_FileForceChain << "2 - TIME | ";
      //  mRVE_FileForceChain << "3 - [X1 Y1 Z1 X2 Y2 Z2 F] of each contact";
      //  mRVE_FileForceChain << std::endl;
      //}

      //if (r_process_info[POST_WRITE_ROSE_DIAGRAM]) {
      //  mRVE_FileRoseDiagram.open("rve_rose_diagram.txt", std::ios::out);
      //  KRATOS_ERROR_IF_NOT(mRVE_FileRoseDiagram) << "Could not open file rve_rose_diagram.txt!" << std::endl;
      //  mRVE_FileRoseDiagram << "1 - STEP | ";
      //  mRVE_FileRoseDiagram << "2 - TIME | ";
      //  mRVE_FileRoseDiagram << "3 - [ARRAY OF ANGLES IN XY PLANE] | ";
      //  mRVE_FileRoseDiagram << "4 - [ARRAY OF AZIMUTH ANGLES]";
      //  mRVE_FileRoseDiagram << std::endl;
      //}

      //if (r_process_info[POST_WRITE_ROSE_DIAGRAM_INNER]) {
      //  mRVE_FileRoseDiagramInner.open("rve_rose_diagram_inner.txt", std::ios::out);
      //  KRATOS_ERROR_IF_NOT(mRVE_FileRoseDiagramInner) << "Could not open file rve_rose_diagram_inner.txt!" << std::endl;
      //  mRVE_FileRoseDiagramInner << "1 - STEP | ";
      //  mRVE_FileRoseDiagramInner << "2 - TIME | ";
      //  mRVE_FileRoseDiagramInner << "3 - [ARRAY OF ANGLES IN XY PLANE] | ";
      //  mRVE_FileRoseDiagramInner << "4 - [ARRAY OF AZIMUTH ANGLES]";
      //  mRVE_FileRoseDiagramInner << std::endl;
      //}

      //if (r_process_info[POST_WRITE_ROSE_DIAGRAM_UNIFORMITY]) {
      //  mRVE_FileRoseDiagramUniformity.open("rve_rose_diagram_uniformity.txt", std::ios::out);
      //  KRATOS_ERROR_IF_NOT(mRVE_FileRoseDiagramUniformity) << "Could not open file rve_rose_diagram_uniformity.txt!" << std::endl;
      //  mRVE_FileRoseDiagramUniformity << "1 - STEP | ";
      //  mRVE_FileRoseDiagramUniformity << "2 - TIME | ";
      //  mRVE_FileRoseDiagramUniformity << "3 - STD DEV XY - ALL | ";
      //  mRVE_FileRoseDiagramUniformity << "4 - STD DEV XY - INNER";
      //  mRVE_FileRoseDiagramUniformity << std::endl;
      //}

      //if (r_process_info[POST_WRITE_ANISOTROPY]) {
      //  mRVE_FileAnisotropy.open("rve_anisotropy.txt", std::ios::out);
      //  KRATOS_ERROR_IF_NOT(mRVE_FileAnisotropy) << "Could not open file rve_anisotropy.txt!" << std::endl;
      //  mRVE_FileAnisotropy << "1 - STEP | ";
      //  mRVE_FileAnisotropy << "2 - TIME | ";
      //  mRVE_FileAnisotropy << "3 - ANISOTROPY - ALL | ";
      //  mRVE_FileAnisotropy << "4 - ANISOTROPY - INNER";
      //  mRVE_FileAnisotropy << std::endl;
      //}

      //if (r_process_info[POST_WRITE_FABRIC_TENSOR]) {
      //  mRVE_FileFabricTensor.open("rve_fabric_tensor.txt", std::ios::out);
      //  KRATOS_ERROR_IF_NOT(mRVE_FileFabricTensor) << "Could not open file rve_fabric_tensor.txt!" << std::endl;
      //  mRVE_FileFabricTensor << "1 - STEP | ";
      //  mRVE_FileFabricTensor << "2 - TIME | ";
      //  mRVE_FileFabricTensor << "3 - [[1,1][1,2][1,3]] | ";
      //  mRVE_FileFabricTensor << "4 - [[2,1][2,2][2,3]] | ";
      //  mRVE_FileFabricTensor << "5 - [[3,1][3,2][3,3]]";
      //  mRVE_FileFabricTensor << std::endl;
      //}

      //if (r_process_info[POST_WRITE_FABRIC_TENSOR_INNER]) {
      //  mRVE_FileFabricTensorInner.open("rve_fabric_tensor_inner.txt", std::ios::out);
      //  KRATOS_ERROR_IF_NOT(mRVE_FileFabricTensorInner) << "Could not open file rve_fabric_tensor_inner.txt!" << std::endl;
      //  mRVE_FileFabricTensorInner << "1 - STEP | ";
      //  mRVE_FileFabricTensorInner << "2 - TIME | ";
      //  mRVE_FileFabricTensorInner << "3 - [[1,1][1,2][1,3]] | ";
      //  mRVE_FileFabricTensorInner << "4 - [[2,1][2,2][2,3]] | ";
      //  mRVE_FileFabricTensorInner << "5 - [[3,1][3,2][3,3]]";
      //  mRVE_FileFabricTensorInner << std::endl;
      //}

      //if (r_process_info[POST_WRITE_STRESS]) {
      //  mRVE_FileStress.open("rve_stresses.txt", std::ios::out);
      //  KRATOS_ERROR_IF_NOT(mRVE_FileStress) << "Could not open file rve_stresses.txt!" << std::endl;
      //  mRVE_FileStress << "1 - STEP | ";
      //  mRVE_FileStress << "2 - TIME | ";
      //  mRVE_FileStress << "3 - WALL STRESS | ";
      //  mRVE_FileStress << "4 - MEAN EFFECTIVE STRESS - ALL | ";
      //  mRVE_FileStress << "5 - DEVIATORIC STRESS - ALL | ";
      //  mRVE_FileStress << "6 - MEAN EFFECTIVE STRESS - INNER | ";
      //  mRVE_FileStress << "7 - DEVIATORIC STRESS - INNER";
      //  mRVE_FileStress << std::endl;
      //}

      //if (r_process_info[POST_WRITE_CAUCHY_TENSOR]) {
      //  mRVE_FileCauchyTensor.open("rve_cauchy_tensor.txt", std::ios::out);
      //  KRATOS_ERROR_IF_NOT(mRVE_FileCauchyTensor) << "Could not open file rve_cauchy_tensor.txt!" << std::endl;
      //  mRVE_FileCauchyTensor << "1 - STEP | ";
      //  mRVE_FileCauchyTensor << "2 - TIME | ";
      //  mRVE_FileCauchyTensor << "3 - [[1,1][1,2][1,3]] | ";
      //  mRVE_FileCauchyTensor << "4 - [[2,1][2,2][2,3]] | ";
      //  mRVE_FileCauchyTensor << "5 - [[3,1][3,2][3,3]]";
      //  mRVE_FileCauchyTensor << std::endl;
      //}

      //if (r_process_info[POST_WRITE_CAUCHY_TENSOR_INNER]) {
      //  mRVE_FileCauchyTensorInner.open("rve_cauchy_tensor_inner.txt", std::ios::out);
      //  KRATOS_ERROR_IF_NOT(mRVE_FileCauchyTensorInner) << "Could not open file rve_cauchy_tensor_inner.txt!" << std::endl;
      //  mRVE_FileCauchyTensorInner << "1 - STEP | ";
      //  mRVE_FileCauchyTensorInner << "2 - TIME | ";
      //  mRVE_FileCauchyTensorInner << "3 - [[1,1][1,2][1,3]] | ";
      //  mRVE_FileCauchyTensorInner << "4 - [[2,1][2,2][2,3]] | ";
      //  mRVE_FileCauchyTensorInner << "5 - [[3,1][3,2][3,3]]";
      //  mRVE_FileCauchyTensorInner << std::endl;
      //}

      //if (r_process_info[POST_WRITE_TANGENT_TENSOR]) {
      //  mRVE_FileTangentTensor.open("rve_tangent_tensor.txt", std::ios::out);
      //  KRATOS_ERROR_IF_NOT(mRVE_FileTangentTensor) << "Could not open file rve_tangent_tensor.txt!" << std::endl;
      //  mRVE_FileTangentTensor << "1 - STEP | ";
      //  mRVE_FileTangentTensor << "2 - TIME | ";
      //  mRVE_FileTangentTensor << "ROW1: [[D1111][D1112][D1113][D1121][D1122][D1123][D1131][D1132][D1133]] | ";
      //  mRVE_FileTangentTensor << "ROW2: [[D1211][D1212][D1213][D1221][D1222][D1223][D1231][D1232][D1233]] | ";
      //  mRVE_FileTangentTensor << "ROW3: [[D1311][D1312][D1313][D1321][D1322][D1323][D1331][D1332][D1333]] | ";
      //  mRVE_FileTangentTensor << "ROW4: [[D2111][D2112][D2113][D2121][D2122][D2123][D2131][D2132][D2133]] | ";
      //  mRVE_FileTangentTensor << "ROW5: [[D2211][D2212][D2213][D2221][D2222][D2223][D2231][D2232][D2233]] | ";
      //  mRVE_FileTangentTensor << "ROW6: [[D2311][D2312][D2313][D2321][D2322][D2323][D2331][D2332][D2333]] | ";
      //  mRVE_FileTangentTensor << "ROW7: [[D3111][D3112][D3113][D3121][D3122][D3123][D3131][D3132][D3133]] | ";
      //  mRVE_FileTangentTensor << "ROW8: [[D3211][D3212][D3213][D3221][D3222][D3223][D3231][D3232][D3233]] | ";
      //  mRVE_FileTangentTensor << "ROW9: [[D3311][D3312][D3313][D3321][D3322][D3323][D3331][D3332][D3333]]";
      //  mRVE_FileTangentTensor << std::endl;
      //}

      //if (r_process_info[POST_WRITE_TANGENT_TENSOR_INNER]) {
      //  mRVE_FileTangentTensorInner.open("rve_tangent_tensor_inner.txt", std::ios::out);
      //  KRATOS_ERROR_IF_NOT(mRVE_FileTangentTensorInner) << "Could not open file rve_tangent_tensor_inner.txt!" << std::endl;
      //  mRVE_FileTangentTensorInner << "1 - STEP | ";
      //  mRVE_FileTangentTensorInner << "2 - TIME | ";
      //  mRVE_FileTangentTensorInner << "ROW1: [[D1111][D1112][D1113][D1121][D1122][D1123][D1131][D1132][D1133]] | ";
      //  mRVE_FileTangentTensorInner << "ROW2: [[D1211][D1212][D1213][D1221][D1222][D1223][D1231][D1232][D1233]] | ";
      //  mRVE_FileTangentTensorInner << "ROW3: [[D1311][D1312][D1313][D1321][D1322][D1323][D1331][D1332][D1333]] | ";
      //  mRVE_FileTangentTensorInner << "ROW4: [[D2111][D2112][D2113][D2121][D2122][D2123][D2131][D2132][D2133]] | ";
      //  mRVE_FileTangentTensorInner << "ROW5: [[D2211][D2212][D2213][D2221][D2222][D2223][D2231][D2232][D2233]] | ";
      //  mRVE_FileTangentTensorInner << "ROW6: [[D2311][D2312][D2313][D2321][D2322][D2323][D2331][D2332][D2333]] | ";
      //  mRVE_FileTangentTensorInner << "ROW7: [[D3111][D3112][D3113][D3121][D3122][D3123][D3131][D3132][D3133]] | ";
      //  mRVE_FileTangentTensorInner << "ROW8: [[D3211][D3212][D3213][D3221][D3222][D3223][D3231][D3232][D3233]] | ";
      //  mRVE_FileTangentTensorInner << "ROW9: [[D3311][D3312][D3313][D3321][D3322][D3323][D3331][D3332][D3333]]";
      //  mRVE_FileTangentTensorInner << std::endl;
      //}

      //if (r_process_info[POST_WRITE_CONDUCTIVITY_TENSOR]) {
      //  mRVE_FileConductivityTensor.open("rve_conductivity_tensor.txt", std::ios::out);
      //  KRATOS_ERROR_IF_NOT(mRVE_FileConductivityTensor) << "Could not open file rve_conductivity_tensor.txt!" << std::endl;
      //  mRVE_FileConductivityTensor << "1 - STEP | ";
      //  mRVE_FileConductivityTensor << "2 - TIME | ";
      //  mRVE_FileConductivityTensor << "3 - [[1,1][1,2][1,3]] | ";
      //  mRVE_FileConductivityTensor << "4 - [[2,1][2,2][2,3]] | ";
      //  mRVE_FileConductivityTensor << "5 - [[3,1][3,2][3,3]]";
      //  mRVE_FileConductivityTensor << std::endl;
      //}

      //if (r_process_info[POST_WRITE_CONDUCTIVITY_TENSOR_INNER]) {
      //  mRVE_FileConductivityTensorInner.open("rve_conductivity_tensor_inner.txt", std::ios::out);
      //  KRATOS_ERROR_IF_NOT(mRVE_FileConductivityTensorInner) << "Could not open file rve_conductivity_tensor_inner.txt!" << std::endl;
      //  mRVE_FileConductivityTensorInner << "1 - STEP | ";
      //  mRVE_FileConductivityTensorInner << "2 - TIME | ";
      //  mRVE_FileConductivityTensorInner << "3 - [[1,1][1,2][1,3]] | ";
      //  mRVE_FileConductivityTensorInner << "4 - [[2,1][2,2][2,3]] | ";
      //  mRVE_FileConductivityTensorInner << "5 - [[3,1][3,2][3,3]]";
      //  mRVE_FileConductivityTensorInner << std::endl;
      //}

      //if (r_process_info[POST_WRITE_FKS]) {
      //  mRVE_FileFKS.open("rve_KFS.txt", std::ios::out);
      //  KRATOS_ERROR_IF_NOT(mRVE_FileFKS) << "Could not open file rve_KFS.txt!" << std::endl;
      //  mRVE_FileFKS << "1 - STEP | ";
      //  mRVE_FileFKS << "2 - TIME | ";
      //  mRVE_FileFKS << "3 - Fxy | ";
      //  mRVE_FileFKS << "4 - Fxx-Fyy | ";
      //  mRVE_FileFKS << "5 - Kxx Kyy Kxy";
      //  mRVE_FileFKS << std::endl;
      //}
    }

    //-----------------------------------------------------------------------------------------------------------------------------------------
    void ExplicitSolverStrategy::RVECloseFiles(void) {
      if (mRVE_FileCoordinatesHistory.is_open())      mRVE_FileCoordinatesHistory.close();
      if (mRVE_FilePorosity.is_open())                mRVE_FilePorosity.close();
      if (mRVE_FileContactNumber.is_open())           mRVE_FileContactNumber.close();
      if (mRVE_FileCoordNumber.is_open())             mRVE_FileCoordNumber.close();
      if (mRVE_FileInnerVolumeParticles.is_open())    mRVE_FileInnerVolumeParticles.close();
      if (mRVE_FileForceChain.is_open())              mRVE_FileForceChain.close();
      if (mRVE_FileRoseDiagram.is_open())             mRVE_FileRoseDiagram.close();
      if (mRVE_FileRoseDiagramInner.is_open())        mRVE_FileRoseDiagramInner.close();
      if (mRVE_FileRoseDiagramUniformity.is_open())   mRVE_FileRoseDiagramUniformity.close();
      if (mRVE_FileAnisotropy.is_open())              mRVE_FileAnisotropy.close();
      if (mRVE_FileFabricTensor.is_open())            mRVE_FileFabricTensor.close();
      if (mRVE_FileFabricTensorInner.is_open())       mRVE_FileFabricTensorInner.close();
      if (mRVE_FileStress.is_open())                  mRVE_FileStress.close();
      if (mRVE_FileCauchyTensor.is_open())            mRVE_FileCauchyTensor.close();
      if (mRVE_FileCauchyTensorInner.is_open())       mRVE_FileCauchyTensorInner.close();
      if (mRVE_FileTangentTensor.is_open())           mRVE_FileTangentTensor.close();
      if (mRVE_FileTangentTensorInner.is_open())      mRVE_FileTangentTensorInner.close();
      if (mRVE_FileConductivityTensor.is_open())      mRVE_FileConductivityTensor.close();
      if (mRVE_FileConductivityTensorInner.is_open()) mRVE_FileConductivityTensorInner.close();
      if (mRVE_FileFKS.is_open())                     mRVE_FileFKS.close();
    }

    //-----------------------------------------------------------------------------------------------------------------------
    void ExplicitSolverStrategy::ClearTriangle(struct triangulateio& rTr)
    {
      KRATOS_TRY

      rTr.pointlist                  = (REAL*)NULL;
      rTr.pointattributelist         = (REAL*)NULL;
      rTr.pointmarkerlist            = (int*)NULL;
      rTr.numberofpoints             = 0;
      rTr.numberofpointattributes    = 0;

      rTr.trianglelist               = (int*)NULL;
      rTr.triangleattributelist      = (REAL*)NULL;
      rTr.trianglearealist           = (REAL*)NULL;
      rTr.neighborlist               = (int*)NULL;
      rTr.numberoftriangles          = 0;
      rTr.numberofcorners            = 3; //for three node triangles
      rTr.numberoftriangleattributes = 0;

      rTr.segmentlist                = (int*)NULL;
      rTr.segmentmarkerlist          = (int*)NULL;
      rTr.numberofsegments           = 0;

      rTr.holelist                   = (REAL*)NULL;
      rTr.numberofholes              = 0;

      rTr.regionlist                 = (REAL*)NULL;
      rTr.numberofregions            = 0;

      rTr.edgelist                   = (int*)NULL;
      rTr.edgemarkerlist             = (int*)NULL;
      rTr.normlist                   = (REAL*)NULL;
      rTr.numberofedges              = 0;

      KRATOS_CATCH("")
    }

    //-----------------------------------------------------------------------------------------------------------------------
    void ExplicitSolverStrategy::FreeTriangle(struct triangulateio& rTr)
    {
      KRATOS_TRY

      if (rTr.numberoftriangles) {
        if (rTr.trianglelist)          trifree(rTr.trianglelist);
        if (rTr.triangleattributelist) trifree(rTr.triangleattributelist);
        if (rTr.trianglearealist)      trifree(rTr.trianglearealist);
        if (rTr.neighborlist)          trifree(rTr.neighborlist);
      }
      if (rTr.segmentlist)       trifree(rTr.segmentlist);
      if (rTr.segmentmarkerlist) trifree(rTr.segmentmarkerlist);
      if (rTr.holelist) {
        delete[] rTr.holelist;
        rTr.numberofholes = 0;
      }
      if (rTr.regionlist) {
        delete[] rTr.regionlist;
        rTr.numberofregions = 0;
      }
      if (rTr.edgelist)       trifree(rTr.edgelist);
      if (rTr.edgemarkerlist) trifree(rTr.edgemarkerlist);
      if (rTr.normlist)       trifree(rTr.normlist);
      if (rTr.numberofpoints) {
        if (rTr.pointlist)          trifree(rTr.pointlist);
        if (rTr.pointattributelist) trifree(rTr.pointattributelist);
        if (rTr.pointmarkerlist)    trifree(rTr.pointmarkerlist);
      }

      KRATOS_CATCH("")
    }

    //==========================================================================================================================================
    // HIERARCHICAL MULTISCALE RVE - FINISH
    //==========================================================================================================================================
}
