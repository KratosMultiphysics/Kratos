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
            mListOfSphericParticles[i]->CalculateRightHandSide(r_process_info, dt, gravity);
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
}
