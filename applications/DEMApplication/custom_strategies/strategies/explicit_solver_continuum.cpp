//
// Author: Miguel AngelCeligueta, maceli@cimne.upc.edu
//

#include "explicit_solver_continuum.h"

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

        r_process_info[SEARCH_CONTROL_VECTOR].resize(mNumberOfThreads);
        for (int i = 0; i < mNumberOfThreads; i++) r_process_info[SEARCH_CONTROL_VECTOR][i] = 0;

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

        if (GetDeltaOption() == 2) {
            SetCoordinationNumber(r_model_part);
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
            for (int i = 0; i < mNumberOfThreads; i++) {
                if (r_process_info[SEARCH_CONTROL_VECTOR][i] == 1) {
                    r_process_info[SEARCH_CONTROL] = 1;
                    if(r_model_part.GetCommunicator().MyPID() == 0) {
                        KRATOS_WARNING("DEM") << "From now on, the search is activated because some failure occurred " << std::endl;
                    }
                    break;
                }
            }
        }

        const int time_step = r_process_info[TIME_STEPS];
        const double time = r_process_info[TIME];
        const bool is_time_to_search_neighbours = (time_step + 1) % GetNStepSearch() == 0 && (time_step > 0); //Neighboring search. Every N times.

        if (r_process_info[SEARCH_CONTROL] > 0) {

            if (is_time_to_search_neighbours) {

                CalculateMaxSearchDistance(); //Modifies r_process_info[AMPLIFIED_CONTINUUM_SEARCH_RADIUS_EXTENSION] // Must be called before the bounding box or it uses non-existent elements

	        if (r_process_info[BOUNDING_BOX_OPTION] && time >= r_process_info[BOUNDING_BOX_START_TIME] && time <= r_process_info[BOUNDING_BOX_STOP_TIME]) {

                    BoundingBoxUtility();
                } else {
                    GetParticleCreatorDestructor()->DestroyParticles(r_model_part);
                    GetParticleCreatorDestructor()->DestroyContactElements(*mpContact_model_part);
                }

                RebuildListOfSphericParticles <SphericContinuumParticle> (r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericContinuumParticles); //These lists are necessary for the loop in SearchNeighbours
                RebuildListOfSphericParticles <SphericParticle> (r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericParticles);

                SetSearchRadiiOnAllParticles(r_model_part, r_process_info[SEARCH_RADIUS_INCREMENT] + r_process_info[AMPLIFIED_CONTINUUM_SEARCH_RADIUS_EXTENSION], 1.0);

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

            //if (r_process_info[BOUNDING_BOX_OPTION] == 1 && has_mpi) {  //This block rebuilds all the bonds between continuum particles

            if (r_process_info[CONTACT_MESH_OPTION]) {
                CreateContactElements();
                InitializeContactElements();
            }
            //}

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

        const double in_coordination_number = r_process_info[COORDINATION_NUMBER];
        double standard_dev = 0;
        double out_coordination_number = ComputeCoordinationNumber(standard_dev);
        int iteration = 0;
        int maxiteration = 100;
        double& added_search_distance = r_process_info[SEARCH_RADIUS_INCREMENT];

        if(r_model_part.GetCommunicator().MyPID() == 0) {
            KRATOS_WARNING("DEM") << "Setting up Coordination Number by increasing or decreasing the search radius... " << std::endl;
        }

        if (in_coordination_number <= 0.0) {
            KRATOS_THROW_ERROR(std::runtime_error, "The specified Coordination Number is less or equal to zero, N.C. = ", in_coordination_number)
        }

        while (fabs(out_coordination_number / in_coordination_number - 1.0) > 1e-3) {
            if (iteration >= maxiteration) break;
            iteration++;
            if(r_model_part.GetCommunicator().MyPID() == 0) { KRATOS_INFO("DEM") <<" * "<<std::flush; }
            if (out_coordination_number == 0.0) {
                KRATOS_WARNING("DEM") << "Coordination Number method not supported in this case" << "\n" << std::endl;
                KRATOS_THROW_ERROR(std::runtime_error, "The specified tangency method is not supported for this problem, please use absolute value instead", " ")
                break;
            }
            added_search_distance *= in_coordination_number / out_coordination_number;
            SetSearchRadiiOnAllParticles(r_model_part, added_search_distance, 1.0);
            SearchNeighbours(); //r_process_info[SEARCH_RADIUS_INCREMENT] will be used inside this function, and it's the variable we are updating in this while
            out_coordination_number = ComputeCoordinationNumber(standard_dev);
        }//while

        if(r_model_part.GetCommunicator().MyPID() == 0) { KRATOS_INFO("DEM") <<std::endl;}

        if (iteration < maxiteration){
            if(r_model_part.GetCommunicator().MyPID() == 0) {
                KRATOS_WARNING("DEM") << "Coordination Number iterative procedure converged after " << iteration << " iterations, to value " << out_coordination_number << " using an extension of " << added_search_distance << ". " << "\n" << std::endl;
                KRATOS_WARNING("DEM") << "Standard deviation for achieved coordination number is " << standard_dev << ". " << "\n" << std::endl;
                KRATOS_WARNING("DEM") << "This means that most particles (about 68% of the total particles, assuming a normal distribution) have a coordination number within " <<  standard_dev << " contacts of the mean (" << out_coordination_number-standard_dev << "–" << out_coordination_number+standard_dev << " contacts). " << "\n" << std::endl;
            }
        }

        else {
            if(r_model_part.GetCommunicator().MyPID() == 0) {
                KRATOS_WARNING("DEM") << "Coordination Number iterative procedure did NOT converge after " << iteration << " iterations. Coordination number reached is " << out_coordination_number << ". " << "\n" << std::endl;
                KRATOS_THROW_ERROR(std::runtime_error, "Please use a Absolute tolerance instead ", " ")
                    //NOTE: if it doesn't converge, problems occur with contact mesh and rigid face contact.
            }
        }

    } //SetCoordinationNumber

    double ContinuumExplicitSolverStrategy::ComputeCoordinationNumber(double& standard_dev) {
        KRATOS_TRY

        ModelPart& r_model_part = GetModelPart();
        ElementsArrayType& pElements = r_model_part.GetCommunicator().LocalMesh().Elements();

        int total_contacts = 0;
        const int number_of_particles = (int) mListOfSphericParticles.size();
        std::vector<int> neighbour_counter;
        std::vector<int> sum;
        double total_sum = 0.0;

        mNumberOfThreads = OpenMPUtils::GetNumThreads();
        neighbour_counter.resize(mNumberOfThreads);
        sum.resize(mNumberOfThreads);

        #pragma omp parallel
        {
            neighbour_counter[OpenMPUtils::ThisThread()] = 0;
            sum[OpenMPUtils::ThisThread()] = 0;
            #pragma omp for
            for (int i = 0; i < number_of_particles; i++) {
                neighbour_counter[OpenMPUtils::ThisThread()] += mListOfSphericParticles[i]->mNeighbourElements.size();
                sum[OpenMPUtils::ThisThread()] += (mListOfSphericParticles[i]->mNeighbourElements.size() - 10.0 )*(mListOfSphericParticles[i]->mNeighbourElements.size() - 10.0 );
            }
        }
        for (int i = 0; i < mNumberOfThreads; i++) {
            total_contacts += neighbour_counter[i];
            total_sum += sum[i];
        }

        int global_total_contacts = r_model_part.GetCommunicator().GetDataCommunicator().SumAll(total_contacts);
        int global_number_of_elements = r_model_part.GetCommunicator().GetDataCommunicator().SumAll((int) pElements.size());

        double coord_number = double(global_total_contacts) / double(global_number_of_elements);

        standard_dev = sqrt(total_sum / number_of_particles);

        return coord_number;

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
            double max_sphere = mListOfSphericContinuumParticles[i]->CalculateMaxSearchDistance(has_mpi, r_process_info);
            if (max_sphere > thread_maxima[OpenMPUtils::ThisThread()]) thread_maxima[OpenMPUtils::ThisThread()] = max_sphere;
        }

        double maximum_across_threads = 0.0;
        for (int i = 0; i < OpenMPUtils::GetNumThreads(); i++) {
            if (thread_maxima[i] > maximum_across_threads) maximum_across_threads = thread_maxima[i];
        }

        r_process_info[AMPLIFIED_CONTINUUM_SEARCH_RADIUS_EXTENSION] = maximum_across_threads;

        const double ratio = r_process_info[AMPLIFIED_CONTINUUM_SEARCH_RADIUS_EXTENSION] / r_process_info[SEARCH_RADIUS_INCREMENT];
        const double max_ratio = r_process_info[MAX_AMPLIFICATION_RATIO_OF_THE_SEARCH_RADIUS];

        static unsigned int counter = 0;
        unsigned int maximum_number_of_prints = 500;

        if ((ratio > max_ratio) && (counter <= maximum_number_of_prints)) {
            KRATOS_INFO("DEM") <<std::endl;
            KRATOS_WARNING("DEM") <<"************************************************************************"<<std::endl;
            KRATOS_WARNING("DEM") <<"WARNING! The automatic extension of the search radius, based on mechanical"<<std::endl;
            KRATOS_WARNING("DEM") <<"reasons, is trying to extend more than "<<max_ratio<<" times the "<<std::endl;
            KRATOS_WARNING("DEM") <<"previously set search radius!"<<std::endl;
            KRATOS_WARNING("DEM") <<"Some bonds might break for search reasons instead of mechanical reasons."<<std::endl;
            KRATOS_WARNING("DEM") <<"The ratio is limited to that value ("<<max_ratio<<" times by the input "<<std::endl;
            KRATOS_WARNING("DEM") <<"variable 'MaxAmplificationRatioOfSearchRadius'"<<std::endl;
            KRATOS_WARNING("DEM") <<"************************************************************************"<<std::endl;
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
