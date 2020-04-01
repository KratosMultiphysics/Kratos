//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Peter Wilson (thanks Klaus Sautter)
//
//

#if !defined(KRATOS_MPM_EXPLICIT_SCHEME)
#define KRATOS_MPM_EXPLICIT_SCHEME

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "includes/element.h"
#include "containers/array_1d.h"
#include "solving_strategies/schemes/scheme.h"
#include "custom_utilities/mpm_boundary_rotation_utility.h"
#include "custom_utilities/mpm_explicit_utilities.h"

namespace Kratos {

    ///@name Kratos Globals
    ///@{

    ///@}
    ///@name Type Definitions
    ///@{

    ///@}

    ///@name  Enum's
    ///@{

    ///@}
    ///@name  Functions
    ///@{

    ///@}
    ///@name Kratos Classes
    ///@{

    /**
     * @class MPMExplicitScheme
     * @ingroup KratosParticle
     * @brief A MPM explicit scheme
     * @details Scheme options include Forward Euler or Central Difference.
     * Stress update options include Update Stress First (USF), Update Stress Last (USL) and Modified Update Stress Last (MUSL).
     *
     * @author Peter
     */
    template <class TSparseSpace,
        class TDenseSpace //= DenseSpace<double>
    >
        class MPMExplicitScheme
        : public Scheme<TSparseSpace, TDenseSpace> {

        public:
            /**@name Type Definitions */

            /*@{ */
            KRATOS_CLASS_POINTER_DEFINITION(MPMExplicitScheme);

            typedef Scheme<TSparseSpace, TDenseSpace>                      BaseType;

            typedef typename BaseType::TDataType                         TDataType;

            typedef typename BaseType::DofsArrayType                 DofsArrayType;

            typedef typename Element::DofsVectorType                DofsVectorType;

            typedef typename BaseType::TSystemMatrixType         TSystemMatrixType;

            typedef typename BaseType::TSystemVectorType         TSystemVectorType;

            typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

            typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

            typedef ModelPart::ElementsContainerType             ElementsArrayType;

            typedef ModelPart::ConditionsContainerType         ConditionsArrayType;

            typedef typename BaseType::Pointer                     BaseTypePointer;

            /// The arrays of elements and nodes
            typedef ModelPart::NodesContainerType NodesArrayType;

            /// Definition for the node iterator
            typedef typename ModelPart::NodeIterator NodeIterator;

            /// Definition of the size type
            typedef std::size_t SizeType;

            /// Definition of the index type
            typedef std::size_t IndexType;

            /// The definition of the numerical limit
            static constexpr double numerical_limit = std::numeric_limits<double>::epsilon();

            ///@}
            ///@name Life Cycle
            ///@{

            /**
             * @brief Default constructor.
             * @details The MPMExplicitScheme method
             * @param MaximumDeltaTime The maximum delta time to be considered
             * @param DeltaTimeFraction The delta ttime fraction
             * @param DeltaTimePredictionLevel The prediction level
             */
            MPMExplicitScheme(
                ModelPart& grid_model_part,
                const int StressUpdateOption,
                const bool isCentralDifference
                )
                : Scheme<TSparseSpace, TDenseSpace>(),
                mr_grid_model_part(grid_model_part),
                mStressUpdateOption(StressUpdateOption),
                mIsCentralDifference(isCentralDifference)
            {
                //mStressUpdateOption = StressUpdateOption; // 0 = USF, 1 = USL, 2 = MUSL
                //mIsCentralDifference = isCentralDifference;

                std::cout << "\n\n =========================== USING MPM EXPLICIT ========================== \n\n" << std::endl;
                if (mIsCentralDifference)
                {
                    std::cout << "\n\n =========================== CENTRAL DIFFERENCE ========================== \n\n" << std::endl;
                }
                else
                {
                    if (mStressUpdateOption < 3 && mStressUpdateOption > -1)
                    {
                        std::cout << "\n\n =========================== FORWARD EULER ========================== \n\n" << std::endl;
                    }
                    else
                    {
                        KRATOS_ERROR << "Invalid MPM explicit scheme constructed." << std::endl;
                    }
                }
            }

            /** Destructor.
            */
            virtual ~MPMExplicitScheme() {}

            ///@}
            ///@name Operators
            ///@{

            /**
             * Clone
             */
            BaseTypePointer Clone() override
            {
                return BaseTypePointer(new MPMExplicitScheme(*this));
            }

            /**
             * @brief This is the place to initialize the Scheme. This is intended to be called just once when the strategy is initialized
             * @param rModelPart The model of the problem to solve
             */
            void Initialize(ModelPart& rModelPart) override
            {
                KRATOS_TRY

                    ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

                // Preparing the time values for the first step (where time = initial_time +
                // dt)
                mTime.Current = r_current_process_info[TIME] + r_current_process_info[DELTA_TIME];
                mTime.Delta = r_current_process_info[DELTA_TIME];
                mTime.Middle = mTime.Current - 0.5 * mTime.Delta;
                mTime.Previous = mTime.Current - mTime.Delta;
                mTime.PreviousMiddle = mTime.Current - 1.5 * mTime.Delta;

                /// Working in 2D/3D (the definition of DOMAIN_SIZE is check in the Check method)
                const SizeType dim = r_current_process_info[DOMAIN_SIZE];

                // Initialize scheme
                if (!BaseType::SchemeIsInitialized())
                    InitializeExplicitScheme(rModelPart, dim);
                else
                    SchemeCustomInitialization(rModelPart, dim);

                BaseType::SetSchemeIsInitialized();

                KRATOS_CATCH("")
            }

            void InitializeExplicitScheme(
                ModelPart& rModelPart,
                const SizeType DomainSize = 3
                )
            {
                KRATOS_TRY

                    /// The array of ndoes
                    NodesArrayType& r_nodes = rModelPart.Nodes();

                // The first iterator of the array of nodes
                const auto it_node_begin = rModelPart.NodesBegin();

#pragma omp parallel for schedule(guided,512)
                for (int i = 0; i < static_cast<int>(r_nodes.size()); ++i) {
                    auto it_node = (it_node_begin + i);


                    array_1d<double, 3>& r_current_residual = it_node->FastGetSolutionStepValue(FORCE_RESIDUAL);
                    for (IndexType j = 0; j < DomainSize; j++) {
                        r_current_residual[j] = 0.0;
                    }
                }

                KRATOS_CATCH("")
            }

            //***************************************************************************
            //***************************************************************************

            /**
             * Performing the update of the solution
             * Incremental update within newton iteration. It updates the state variables at the end of the time step: u_{n+1}^{k+1}= u_{n+1}^{k}+ \Delta u
             * @param r_model_part
             * @param rDofSet set of all primary variables
             * @param A	LHS matrix
             * @param Dx incremental update of primary variables
             * @param b RHS Vector
             */
            void Update(
                ModelPart& r_model_part,
                DofsArrayType& rDofSet,
                TSystemMatrixType& A,
                TSystemVectorType& Dx,
                TSystemVectorType& b) override
            {
                KRATOS_TRY
                    // The current process info
                    ProcessInfo& r_current_process_info = r_model_part.GetProcessInfo();

                // The array of nodes
                NodesArrayType& r_nodes = r_model_part.Nodes();

                /// Working in 2D/3D (the definition of DOMAIN_SIZE is check in the Check method)
                const SizeType dim = r_current_process_info[DOMAIN_SIZE];

                // Step Update
                // The first step is time =  initial_time ( 0.0) + delta time
                mTime.Current = r_current_process_info[TIME];
                mTime.Delta = r_current_process_info[DELTA_TIME];

                mTime.Middle = mTime.Current - 0.50 * mTime.Delta;
                mTime.Previous = mTime.Current - 1.00 * mTime.Delta;
                mTime.PreviousMiddle = mTime.Middle - 1.00 * mTime.Delta;

                if (mTime.Previous < 0.0) mTime.Previous = 0.00;
                if (mTime.PreviousMiddle < 0.0) mTime.PreviousMiddle = 0.00;
                // The iterator of the first node
                const auto it_node_begin = r_model_part.NodesBegin();

                // Getting dof position
                const IndexType disppos = it_node_begin->GetDofPosition(DISPLACEMENT_X);

                // TODO enable parallel again
                //#pragma omp parallel for schedule(guided,512)
                for (int i = 0; i < static_cast<int>(r_nodes.size()); ++i) {
                    // Current step information "N+1" (before step update).
                    this->UpdateTranslationalDegreesOfFreedom(it_node_begin + i, disppos, dim);
                } // for Node parallel
                KRATOS_CATCH("")
            }

            //***************************************************************************
            //***************************************************************************

            void UpdateTranslationalDegreesOfFreedom(
                NodeIterator itCurrentNode,
                const IndexType DisplacementPosition,
                const SizeType DomainSize = 3
                )
            {
                std::array<bool, 3> fix_displacements = { false, false, false };
                fix_displacements[0] = (itCurrentNode->GetDof(DISPLACEMENT_X, DisplacementPosition).IsFixed());
                fix_displacements[1] = (itCurrentNode->GetDof(DISPLACEMENT_Y, DisplacementPosition + 1).IsFixed());
                if (DomainSize == 3)
                    fix_displacements[2] = (itCurrentNode->GetDof(DISPLACEMENT_Z, DisplacementPosition + 2).IsFixed());

                array_1d<double, 3>& r_nodal_momenta = itCurrentNode->FastGetSolutionStepValue(NODAL_MOMENTUM);
                array_1d<double, 3>& r_current_residual = itCurrentNode->FastGetSolutionStepValue(FORCE_RESIDUAL);

                double alpha = 1.0;
                if (mIsCentralDifference)
                {
                    alpha = 0.5; // factor since we are only adding the corrector here
                }

                for (IndexType j = 0; j < DomainSize; j++) {
                    if (fix_displacements[j]) {
                        r_nodal_momenta[j] = 0.0;
                        r_current_residual[j] = 0.0;
                    }
                    else
                    {
                        r_nodal_momenta[j] += alpha * mTime.Delta * r_current_residual[j];
                    }

                } // for DomainSize

                // We need to set updated grid velocity here if we are using the USL formulation
                if (mStressUpdateOption == 1)
                {
                    array_1d<double, 3>& r_current_velocity = itCurrentNode->FastGetSolutionStepValue(VELOCITY);
                    r_current_velocity.clear();
                    const double nodal_mass = itCurrentNode->FastGetSolutionStepValue(NODAL_MASS);
                    if (nodal_mass > numerical_limit)
                    {
                        for (IndexType j = 0; j < DomainSize; j++)
                        {
                            r_current_velocity[j] = r_nodal_momenta[j] / nodal_mass;
                        } // for DomainSize
                    }
                }
            }


            /**
            This is the place to initialize the elements.
            This is intended to be called just once when the strategy is initialized
             */
            void InitializeElements(ModelPart& rModelPart) override
            {
                KRATOS_TRY

                    int num_threads = OpenMPUtils::GetNumThreads();
                OpenMPUtils::PartitionVector element_partition;
                OpenMPUtils::DivideInPartitions(rModelPart.Elements().size(), num_threads, element_partition);

#pragma omp parallel
                {
                    int k = OpenMPUtils::ThisThread();
                    ElementsArrayType::iterator element_begin = rModelPart.Elements().begin() + element_partition[k];
                    ElementsArrayType::iterator element_end = rModelPart.Elements().begin() + element_partition[k + 1];

                    for (ElementsArrayType::iterator itElem = element_begin; itElem != element_end; itElem++)
                    {
                        itElem->Initialize(); // function to initialize the element
                    }
                }

                this->mElementsAreInitialized = true;

                KRATOS_CATCH("")
            }

            //***************************************************************************
            //***************************************************************************

            /**
            This is the place to initialize the conditions.
            This is intended to be called just once when the strategy is initialized
            */
            void InitializeConditions(ModelPart& rModelPart) override
            {
                KRATOS_TRY

                    KRATOS_ERROR_IF(this->mElementsAreInitialized == false) << "Before initilizing Conditions, initialize Elements FIRST" << std::endl;

                int num_threads = OpenMPUtils::GetNumThreads();
                OpenMPUtils::PartitionVector condition_partition;
                OpenMPUtils::DivideInPartitions(rModelPart.Conditions().size(), num_threads, condition_partition);

#pragma omp parallel
                {
                    int k = OpenMPUtils::ThisThread();
                    ConditionsArrayType::iterator condition_begin = rModelPart.Conditions().begin() + condition_partition[k];
                    ConditionsArrayType::iterator condition_end = rModelPart.Conditions().begin() + condition_partition[k + 1];

                    for (ConditionsArrayType::iterator itCond = condition_begin; itCond != condition_end; itCond++)
                    {
                        itCond->Initialize(); // Function to initialize the condition
                    }
                }

                this->mConditionsAreInitialized = true;
                KRATOS_CATCH("")
            }

            //***************************************************************************
            //***************************************************************************

            /**
             * initializes time step solution
             * only for reasons if the time step solution is restarted
             * @param r_model_part
             * @param A	LHS matrix
             * @param Dx incremental update of primary variables
             * @param b RHS Vector
             */
            void InitializeSolutionStep(
                ModelPart& r_model_part,
                TSystemMatrixType& A,
                TSystemVectorType& Dx,
                TSystemVectorType& b) override
            {
                KRATOS_TRY

                    ProcessInfo CurrentProcessInfo = r_model_part.GetProcessInfo();
                BaseType::InitializeSolutionStep(r_model_part, A, Dx, b);
                // LOOP OVER THE GRID NODES PERFORMED FOR CLEAR ALL NODAL INFORMATION
#pragma omp parallel for
                for (int iter = 0; iter < static_cast<int>(mr_grid_model_part.Nodes().size()); ++iter)
                {
                    auto i = mr_grid_model_part.NodesBegin() + iter;

                    // Variables to be cleaned
                    double& nodal_mass = (i)->FastGetSolutionStepValue(NODAL_MASS);
                    array_1d<double, 3 >& nodal_momentum = (i)->FastGetSolutionStepValue(NODAL_MOMENTUM);
                    array_1d<double, 3 >& nodal_inertia = (i)->FastGetSolutionStepValue(NODAL_INERTIA);
                    array_1d<double, 3 >& nodal_force = (i)->FastGetSolutionStepValue(FORCE_RESIDUAL);
                    array_1d<double, 3 >& nodal_displacement = (i)->FastGetSolutionStepValue(DISPLACEMENT);
                    array_1d<double, 3 >& nodal_velocity = (i)->FastGetSolutionStepValue(VELOCITY);
                    array_1d<double, 3 >& nodal_acceleration = (i)->FastGetSolutionStepValue(ACCELERATION);


                    double& nodal_old_pressure = (i)->FastGetSolutionStepValue(PRESSURE, 1);
                    double& nodal_pressure = (i)->FastGetSolutionStepValue(PRESSURE);
                    if (i->SolutionStepsDataHas(NODAL_MPRESSURE)) {
                        double& nodal_mpressure = (i)->FastGetSolutionStepValue(NODAL_MPRESSURE);
                        nodal_mpressure = 0.0;
                    }

                    // Clear
                    nodal_mass = 0.0;
                    nodal_momentum.clear();
                    nodal_inertia.clear();
                    nodal_force.clear();

                    nodal_displacement.clear();
                    nodal_velocity.clear();
                    nodal_acceleration.clear();
                    nodal_old_pressure = 0.0;
                    nodal_pressure = 0.0;
                }

                // Extrapolate from Material Point Elements and Conditions
                Scheme<TSparseSpace, TDenseSpace>::InitializeSolutionStep(r_model_part, A, Dx, b);

                // If we are updating stress first (USF), calculate nodal velocities from momenta and apply BCs
                if (mStressUpdateOption == 0 || mIsCentralDifference)
                {
                    const IndexType DisplacementPosition = mr_grid_model_part.NodesBegin()->GetDofPosition(DISPLACEMENT_X);

#pragma omp parallel for
                    for (int iter = 0; iter < static_cast<int>(mr_grid_model_part.Nodes().size()); ++iter)
                    {
                        auto i = mr_grid_model_part.NodesBegin() + iter;
                        const SizeType DomainSize = CurrentProcessInfo[DOMAIN_SIZE];
                        double& nodal_mass = (i)->FastGetSolutionStepValue(NODAL_MASS);
                        array_1d<double, 3 >& nodal_momentum = (i)->FastGetSolutionStepValue(NODAL_MOMENTUM);
                        array_1d<double, 3 >& nodal_velocity = (i)->FastGetSolutionStepValue(VELOCITY);

                        std::array<bool, 3> fix_displacements = { false, false, false };
                        fix_displacements[0] = (i->GetDof(DISPLACEMENT_X, DisplacementPosition).IsFixed());
                        fix_displacements[1] = (i->GetDof(DISPLACEMENT_Y, DisplacementPosition + 1).IsFixed());
                        if (DomainSize == 3)
                            fix_displacements[2] = (i->GetDof(DISPLACEMENT_Z, DisplacementPosition + 2).IsFixed());

                        if (nodal_mass > numerical_limit)
                        {
                            for (IndexType j = 0; j < DomainSize; j++)
                            {
                                if (fix_displacements[j])
                                {
                                    nodal_velocity[j] = 0.0;
                                }
                                else
                                {
                                    nodal_velocity[j] = nodal_momentum[j] / nodal_mass;
                                }
                            }
                        }
                    }
                }
                KRATOS_CATCH("")
            }


            //*******************************************************FinalizeNonLinIteration********************
            //***************************************************************************
            /**
            Function called once at the end of a solution step, after convergence is reached if
            an iterative process is needed
             */
            void FinalizeSolutionStep(
                ModelPart& rModelPart,
                TSystemMatrixType& A,
                TSystemVectorType& Dx,
                TSystemVectorType& b) override
            {
                KRATOS_TRY

                    ElementsArrayType& rElements = rModelPart.Elements();
                const ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

                if (mStressUpdateOption == 2)
                {
                    // TODO maybe move to strategy so we can use the base scheme function
                    PerformModifiedUpdateStressLastMapping(rCurrentProcessInfo, rModelPart, rElements);
                }

                // Definition of the first element iterator
                const auto it_elem_begin = rModelPart.ElementsBegin();

                // Finalizes solution step for all of the elements
#pragma omp parallel for
                for (int i = 0; i < static_cast<int>(rElements.size()); ++i) {
                    auto it_elem = it_elem_begin + i;
                    it_elem->FinalizeSolutionStep(rCurrentProcessInfo);
                }

                // Definition of the first condition iterator
                const auto it_cond_begin = rModelPart.ConditionsBegin();

                // Finalizes solution step for all of the conditions
#pragma omp parallel for
                for (int i = 0; i < static_cast<int>(rModelPart.Conditions().size()); ++i) {
                    auto it_cond = it_cond_begin + i;
                    it_cond->FinalizeSolutionStep(rCurrentProcessInfo);
                }

                KRATOS_CATCH("")
            }

            //***************************************************************************
            //***************************************************************************
            void PerformModifiedUpdateStressLastMapping(const ProcessInfo& rCurrentProcessInfo, ModelPart& rModelPart, ElementsArrayType& rElements)
            {
                // MUSL stress update. This works by projecting the updated particle
                // velocity back to the nodes. The nodal velocity field is then
                // used for stress computations.

                // We need to call 'FinalizeSolutionStep' twice. 
                // First, to update the particles and then aggregate the new particle velocities on the grid.
                // Second, to calculate the stresses from the grid velocity
                for (int iter = 0; iter < static_cast<int>(mr_grid_model_part.Nodes().size()); ++iter)
                {
                    auto i = mr_grid_model_part.NodesBegin() + iter;
                    (i)->SetValue(MUSL_VELOCITY_FIELD_IS_COMPUTED, false);
                }


                // Call each particle and aggregate the nodal velocity field
                for (ElementsArrayType::iterator it = rElements.begin(); it != rElements.end(); ++it)
                {
                    (it)->FinalizeSolutionStep(rCurrentProcessInfo);
                }


                // Reapply dirichlet BCs to MUSL velocity field
                const SizeType DomainSize = rCurrentProcessInfo[DOMAIN_SIZE];
                NodesArrayType& r_nodes = rModelPart.Nodes();
                const auto it_node_begin = rModelPart.NodesBegin();
                const IndexType DisplacementPosition = it_node_begin->GetDofPosition(DISPLACEMENT_X);

                for (int i = 0; i < static_cast<int>(r_nodes.size()); ++i)
                {
                    NodeIterator itCurrentNode = it_node_begin + i;

                    std::array<bool, 3> fix_displacements = { false, false, false };
                    fix_displacements[0] = (itCurrentNode->GetDof(DISPLACEMENT_X, DisplacementPosition).IsFixed());
                    fix_displacements[1] = (itCurrentNode->GetDof(DISPLACEMENT_Y, DisplacementPosition + 1).IsFixed());
                    if (DomainSize == 3)
                        fix_displacements[2] = (itCurrentNode->GetDof(DISPLACEMENT_Z, DisplacementPosition + 2).IsFixed());

                    array_1d<double, 3>& r_current_velocity = itCurrentNode->FastGetSolutionStepValue(VELOCITY);

                    for (IndexType j = 0; j < DomainSize; j++)
                    {
                        if (fix_displacements[j])
                        {
                            r_current_velocity[j] = 0.0;
                        }
                    }
                }
            }



            void InitializeNonLinIteration(ModelPart& r_model_part,
                TSystemMatrixType& A,
                TSystemVectorType& Dx,
                TSystemVectorType& b) override
            {
                KRATOS_TRY

                    // This calculates the stresses before momenta update.
                    // Used for USF and Central Difference method
                    if (mStressUpdateOption == 0 || mIsCentralDifference)
                    {
                        ElementsArrayType& pElements = r_model_part.Elements();
                        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

                        for (ElementsArrayType::iterator it = pElements.begin(); it != pElements.end(); ++it)
                        {
                            (it)->InitializeNonLinearIteration(CurrentProcessInfo);
                        }

                        ConditionsArrayType& pConditions = r_model_part.Conditions();
                        for (ConditionsArrayType::iterator it = pConditions.begin(); it != pConditions.end(); ++it)
                        {
                            (it)->InitializeNonLinearIteration(CurrentProcessInfo);
                        }
                    }

                KRATOS_CATCH("")
            }

            //***************************************************************************
            //***************************************************************************

            void InitializeNonLinearIteration(Condition::Pointer rCurrentCondition,
                ProcessInfo& CurrentProcessInfo) override
            {
                (rCurrentCondition)->InitializeNonLinearIteration(CurrentProcessInfo);
            }


            //***************************************************************************
            //***************************************************************************

            void InitializeNonLinearIteration(Element::Pointer rCurrentElement,
                ProcessInfo& CurrentProcessInfo) override
            {
                (rCurrentElement)->InitializeNonLinearIteration(CurrentProcessInfo);
            }

            //***************************************************************************
            //***************************************************************************

            //***************************************************************************
            //***************************************************************************


            //***************************************************************************
            //***************************************************************************

            /** Function that returns the list of Degrees of freedom to be
            assembled in the system for a Given Element
             */
            void GetElementalDofList(
                Element::Pointer rCurrentElement,
                Element::DofsVectorType& ElementalDofList,
                ProcessInfo& CurrentProcessInfo) override
            {
                rCurrentElement->GetDofList(ElementalDofList, CurrentProcessInfo);
            }

            //***************************************************************************
            //***************************************************************************

            /** Function that returns the list of Degrees of freedom to be
            assembled in the system for a Given Element
             */
            void GetConditionDofList(
                Condition::Pointer rCurrentCondition,
                Element::DofsVectorType& ConditionDofList,
                ProcessInfo& CurrentProcessInfo) override
            {
                rCurrentCondition->GetDofList(ConditionDofList, CurrentProcessInfo);
            }

            //***************************************************************************
            //***************************************************************************

            /**
             * This function is designed to be called once to perform all the checks needed
             * on the input provided. Checks can be "expensive" as the function is designed
             * to catch user's errors.
             * @param r_model_part
             * @return 0 all ok
             */
            int Check(ModelPart& r_model_part) override
            {
                KRATOS_TRY

                    int err = Scheme<TSparseSpace, TDenseSpace>::Check(r_model_part);
                if (err != 0) return err;

                //check that the variables are correctly initialized
                KRATOS_ERROR_IF(DISPLACEMENT.Key() == 0) << "DISPLACEMENT has Key zero! (check if the application is correctly registered" << std::endl;
                KRATOS_ERROR_IF(VELOCITY.Key() == 0) << "VELOCITY has Key zero! (check if the application is correctly registered" << std::endl;
                KRATOS_ERROR_IF(ACCELERATION.Key() == 0) << "ACCELERATION has Key zero! (check if the application is correctly registered" << std::endl;

                //check that variables are correctly allocated
                for (ModelPart::NodesContainerType::iterator it = r_model_part.NodesBegin();
                    it != r_model_part.NodesEnd(); it++)
                {
                    KRATOS_ERROR_IF(it->SolutionStepsDataHas(DISPLACEMENT) == false) << "DISPLACEMENT variable is not allocated for node " << it->Id() << std::endl;
                    KRATOS_ERROR_IF(it->SolutionStepsDataHas(VELOCITY) == false) << "VELOCITY variable is not allocated for node " << it->Id() << std::endl;
                    KRATOS_ERROR_IF(it->SolutionStepsDataHas(ACCELERATION) == false) << "ACCELERATION variable is not allocated for node " << it->Id() << std::endl;
                }

                //check that dofs exist
                for (ModelPart::NodesContainerType::iterator it = r_model_part.NodesBegin();
                    it != r_model_part.NodesEnd(); it++)
                {
                    KRATOS_ERROR_IF(it->HasDofFor(DISPLACEMENT_X) == false) << "Missing DISPLACEMENT_X dof on node " << it->Id() << std::endl;
                    KRATOS_ERROR_IF(it->HasDofFor(DISPLACEMENT_Y) == false) << "Missing DISPLACEMENT_Y dof on node " << it->Id() << std::endl;
                    KRATOS_ERROR_IF(it->HasDofFor(DISPLACEMENT_Z) == false) << "Missing DISPLACEMENT_Z dof on node " << it->Id() << std::endl;
                }

                //check for minimum value of the buffer index
                KRATOS_ERROR_IF(r_model_part.GetBufferSize() < 2) << "Insufficient buffer size. Buffer size should be greater than 2. Current size is" << r_model_part.GetBufferSize() << std::endl;

                return 0;
                KRATOS_CATCH("")
            }

            virtual void SchemeCustomInitialization(
                ModelPart& rModelPart,
                const SizeType DomainSize = 3
                )
            {
                KRATOS_TRY

                    // The array containing the nodes
                    NodesArrayType& r_nodes = rModelPart.Nodes();

                // The fisrt node interator
                const auto it_node_begin = rModelPart.NodesBegin();

                // Auxiliar zero array
                const array_1d<double, 3> zero_array = ZeroVector(3);

                // Getting dof position
                const IndexType disppos = it_node_begin->GetDofPosition(DISPLACEMENT_X);

#pragma omp parallel for schedule(guided,512)
                for (int i = 0; i < static_cast<int>(r_nodes.size()); ++i) {
                    // Current step information "N+1" (before step update).
                    auto it_node = it_node_begin + i;

                    const double nodal_mass = it_node->GetValue(NODAL_MASS);

                    const array_1d<double, 3>& r_current_residual = it_node->FastGetSolutionStepValue(FORCE_RESIDUAL);

                    array_1d<double, 3>& r_current_velocity = it_node->FastGetSolutionStepValue(VELOCITY);
                    //             array_1d<double,3>& r_current_displacement = it_node->FastGetSolutionStepValue(DISPLACEMENT);
                    array_1d<double, 3>& r_current_acceleration = it_node->FastGetSolutionStepValue(ACCELERATION);

                    // Solution of the explicit equation:
                    if (nodal_mass > numerical_limit) {
                        r_current_acceleration = r_current_residual / nodal_mass;
                    }
                    else {
                        r_current_acceleration = zero_array;
                    }

                    std::array<bool, 3> fix_displacements = { false, false, false };

                    fix_displacements[0] = (it_node->GetDof(DISPLACEMENT_X, disppos).IsFixed());
                    fix_displacements[1] = (it_node->GetDof(DISPLACEMENT_Y, disppos + 1).IsFixed());
                    if (DomainSize == 3)
                        fix_displacements[2] = (it_node->GetDof(DISPLACEMENT_Z, disppos + 2).IsFixed());

                    for (IndexType j = 0; j < DomainSize; j++) {
                        if (fix_displacements[j]) {
                            r_current_acceleration[j] = 0.0;
                        }
                        r_current_velocity[j] += (mTime.Previous - mTime.PreviousMiddle) * r_current_acceleration[j];
                        // r_current_displacement[j]  = 0.0;
                    } // for DomainSize

                }     // for node parallel

                mTime.Previous = mTime.Current;
                mTime.PreviousMiddle = mTime.Middle;
                KRATOS_CATCH("")
            }

            void Calculate_RHS_Contribution(
                Element::Pointer pCurrentElement,
                LocalSystemVectorType& RHS_Contribution,
                Element::EquationIdVectorType& EquationId,
                ProcessInfo& rCurrentProcessInfo
                ) override
            {
                KRATOS_TRY

                    this->TCalculate_RHS_Contribution(pCurrentElement, RHS_Contribution, rCurrentProcessInfo);
                KRATOS_CATCH("")
            }

            /**
             * @brief Functions that calculates the RHS of a "condition" object
             * @param pCondition The condition to compute
             * @param RHS_Contribution The RHS vector contribution
             * @param EquationId The ID's of the condition degrees of freedom
             * @param rCurrentProcessInfo The current process info instance
             */
            void Condition_Calculate_RHS_Contribution(
                Condition::Pointer pCurrentCondition,
                LocalSystemVectorType& RHS_Contribution,
                Element::EquationIdVectorType& EquationId,
                ProcessInfo& rCurrentProcessInfo
                ) override
            {
                KRATOS_TRY

                    this->TCalculate_RHS_Contribution(pCurrentCondition, RHS_Contribution, rCurrentProcessInfo);

                KRATOS_CATCH("")
            }

            template <typename TObjectType>
            void TCalculate_RHS_Contribution(
                TObjectType pCurrentEntity,
                LocalSystemVectorType& RHS_Contribution,
                ProcessInfo& rCurrentProcessInfo
                )
            {
                KRATOS_TRY

                pCurrentEntity->CalculateRightHandSide(RHS_Contribution, rCurrentProcessInfo);
                pCurrentEntity->AddExplicitContribution(RHS_Contribution, RESIDUAL_VECTOR, FORCE_RESIDUAL, rCurrentProcessInfo);
                //pCurrentEntity->AddExplicitContribution(RHS_Contribution, RESIDUAL_VECTOR, MOMENT_RESIDUAL, rCurrentProcessInfo);

                KRATOS_CATCH("")
            }


            /*@} */
            /**@name Operations */
            /*@{ */
            /*@} */
            /**@name Access */
            /*@{ */
            /*@} */
            /**@name Inquiry */
            /*@{ */
            /*@} */
            /**@name Friends */
            /*@{ */

        protected:
            /**@name Static Member Variables */
            /*@{ */
            /*@} */
            /**@name Member Variables */
            /*@{ */

            struct DeltaTimeParameters {
                double PredictionLevel; // 0, 1, 2 // NOTE: Should be a integer?
                double Maximum;         // Maximum delta time
                double Fraction;        // Fraction of the delta time
            };

            /**
             * @brief This struct contains the details of the time variables
             */
            struct TimeVariables {
                double PreviousMiddle; // n-1/2
                double Previous;       // n
                double Middle;         // n+1/2
                double Current;        // n+1

                double Delta;          // Time step
            };

            ///@name Protected static Member Variables
            ///@{

            TimeVariables mTime;            /// This struct contains the details of the time variables

            ModelPart& mr_grid_model_part;

            const int mStressUpdateOption; // 0 = USF, 1 = USL, 2 = MUSL
            const bool mIsCentralDifference;

            /*@} */
            /**@name Protected Operators*/
            /*@{ */

            /*@} */
            /**@name Protected Operations*/
            /*@{ */
            /*@} */
            /**@name Protected  Access */
            /*@{ */
            /*@} */
            /**@name Protected Inquiry */
            /*@{ */
            /*@} */
            /**@name Protected LifeCycle */
            /*@{ */
        private:
            /**@name Static Member Variables */
            /*@{ */
            /*@} */
            /**@name Member Variables */
            /*@{ */
            /*@} */
            /**@name Private Operators*/
            /*@{ */
            /*@} */
            /**@name Private Operations*/
            /*@{ */
            /*@} */
            /**@name Private  Access */
            /*@{ */
            /*@} */
            /**@name Private Inquiry */
            /*@{ */
            /*@} */
            /**@name Unaccessible methods */
            /*@{ */
    }; /* Class MPMExplicitScheme */
}  /* namespace Kratos.*/

#endif /* KRATOS_MPM_EXPLICIT_SCHEME defined */