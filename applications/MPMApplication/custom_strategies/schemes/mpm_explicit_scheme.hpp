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
#include "solving_strategies/schemes/scheme.h"
#include "custom_utilities/mpm_boundary_rotation_utility.h"
#include "custom_utilities/mpm_explicit_utilities.h"

namespace Kratos {
    /**
     * @class MPMExplicitScheme
     * @ingroup KratosMPM
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

            /**
             * @brief Default constructor.
             * @details The MPMExplicitScheme method
             * @param MaximumDeltaTime The maximum delta time to be considered
             * @param DeltaTimeFraction The delta ttime fraction
             * @param DeltaTimePredictionLevel The prediction level
             */
            explicit MPMExplicitScheme(
                ModelPart& grid_model_part
                )
                : Scheme<TSparseSpace, TDenseSpace>(),
                mr_grid_model_part(grid_model_part)
            {
            }

            /// Destructor.
            virtual ~MPMExplicitScheme() {}

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
                // Initialize scheme
                BaseType::SetSchemeIsInitialized();
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
                const ProcessInfo& r_current_process_info = r_model_part.GetProcessInfo();

                /// Working in 2D/3D (the definition of DOMAIN_SIZE is check in the Check method)
                const SizeType dim = r_current_process_info[DOMAIN_SIZE];
                const double delta_time = r_current_process_info[DELTA_TIME];

                // The iterator of the first node
                const auto it_node_begin = r_model_part.NodesBegin();

                // Getting dof position
                const IndexType disppos = it_node_begin->GetDofPosition(DISPLACEMENT_X);

                #pragma omp parallel for schedule(guided,512)
                for (int i = 0; i < static_cast<int>(r_model_part.Nodes().size()); ++i) {
                    auto it_node = it_node_begin + i;
                    if ((it_node)->Is(ACTIVE))
                    {
                        this->UpdateTranslationalDegreesOfFreedom(r_current_process_info, it_node, disppos, delta_time, dim);
                    }
                } // for Node parallel

                KRATOS_CATCH("")
            }

            //***************************************************************************
            //***************************************************************************

            void UpdateTranslationalDegreesOfFreedom(
                const ProcessInfo& r_current_process_info,
                NodeIterator itCurrentNode,
                const IndexType DisplacementPosition,
                const double delta_time,
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

                const double gamma = (r_current_process_info.GetValue(IS_EXPLICIT_CENTRAL_DIFFERENCE))
                    ? 0.5
                    : 1.0;
                    // we are only adding the central difference corrector here

                for (IndexType j = 0; j < DomainSize; j++)
                {
                    if (fix_displacements[j])
                    {
                        r_nodal_momenta[j] = 0.0;
                        r_current_residual[j] = 0.0;
                    }
                    else {
                        r_nodal_momenta[j] += gamma * delta_time * r_current_residual[j];
                    }
                } // for DomainSize
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

                const ProcessInfo& rCurrentProcessInfo = r_model_part.GetProcessInfo();
                BaseType::InitializeSolutionStep(r_model_part, A, Dx, b);
                #pragma omp parallel for
                for (int iter = 0; iter < static_cast<int>(mr_grid_model_part.Nodes().size()); ++iter)
                {
                    auto i = mr_grid_model_part.NodesBegin() + iter;
                    
                    if(i->Is(ACTIVE))
                    {
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
                }

                // Extrapolate from Material Point Elements and Conditions
                Scheme<TSparseSpace, TDenseSpace>::InitializeSolutionStep(r_model_part, A, Dx, b);

                // If we are updating stress before momenta update (USF and central difference),
                if (rCurrentProcessInfo.GetValue(EXPLICIT_STRESS_UPDATE_OPTION) == 0 ||
                    rCurrentProcessInfo.GetValue(IS_EXPLICIT_CENTRAL_DIFFERENCE))
                {
                    // calculate nodal velocities from momenta and apply BCs
                    calculateGridVelocityAndApplyDirichletBC(rCurrentProcessInfo,true);

                    // calculate stresses
                    const auto it_elem_begin = r_model_part.ElementsBegin();
                    #pragma omp parallel for
                    for (int i = 0; i < static_cast<int>(r_model_part.Elements().size()); ++i) {
                        auto it_elem = it_elem_begin + i;
                        std::vector<bool> dummy;
                        it_elem->CalculateOnIntegrationPoints(CALCULATE_EXPLICIT_MP_STRESS, dummy, rCurrentProcessInfo);
                    }
                }
                KRATOS_CATCH("")
            }

            /// Apply Dirichlet BCs to nodal velocity field
            void calculateGridVelocityAndApplyDirichletBC(
                const ProcessInfo rCurrentProcessInfo,
                bool calculateVelocityFromMomenta = false)
            {
                KRATOS_TRY

                const IndexType DisplacementPosition = mr_grid_model_part.NodesBegin()->GetDofPosition(DISPLACEMENT_X);
                const SizeType DomainSize = rCurrentProcessInfo[DOMAIN_SIZE];

                #pragma omp parallel for
                for (int iter = 0; iter < static_cast<int>(mr_grid_model_part.Nodes().size()); ++iter)
                {
                    NodeIterator i = mr_grid_model_part.NodesBegin() + iter;

                    if ((i)->Is(ACTIVE))
                    {
                        double& nodal_mass = (i)->FastGetSolutionStepValue(NODAL_MASS);
                        array_1d<double, 3 >& nodal_momentum = (i)->FastGetSolutionStepValue(NODAL_MOMENTUM);
                        array_1d<double, 3 >& nodal_velocity = (i)->FastGetSolutionStepValue(VELOCITY);

                        std::array<bool, 3> fix_displacements = { false, false, false };
                        fix_displacements[0] = (i->GetDof(DISPLACEMENT_X, DisplacementPosition).IsFixed());
                        fix_displacements[1] = (i->GetDof(DISPLACEMENT_Y, DisplacementPosition + 1).IsFixed());
                        if (DomainSize == 3)
                            fix_displacements[2] = (i->GetDof(DISPLACEMENT_Z, DisplacementPosition + 2).IsFixed());

                        for (IndexType j = 0; j < DomainSize; j++)
                        {
                            if (fix_displacements[j])
                            {
                                nodal_velocity[j] = 0.0;
                            }
                            else if (calculateVelocityFromMomenta && nodal_mass > numerical_limit)
                            {
                                nodal_velocity[j] = nodal_momentum[j] / nodal_mass;
                            }
                        }
                    }
                }

                KRATOS_CATCH("")
            }


            /** Function called once at the end of a solution step, after convergence is reached if
            an iterative process is needed */
            void FinalizeSolutionStep(
                ModelPart& rModelPart,
                TSystemMatrixType& A,
                TSystemVectorType& Dx,
                TSystemVectorType& b) override
            {
                KRATOS_TRY

                ElementsArrayType& rElements = rModelPart.Elements();
                const ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();
                const auto it_elem_begin = rModelPart.ElementsBegin();


                // map grid to MPs
                #pragma omp parallel for
                for (int i = 0; i < static_cast<int>(rElements.size()); ++i)
                {
                    auto it_elem = it_elem_begin + i;
                    std::vector<bool> dummy;
                    it_elem->CalculateOnIntegrationPoints(EXPLICIT_MAP_GRID_TO_MP, dummy, rCurrentProcessInfo);
                }

                //update stress after momenta update for USL(1) and MUSL(2)
                if (rCurrentProcessInfo.GetValue(EXPLICIT_STRESS_UPDATE_OPTION) > 0)
                {
                    this->CalculateUpdatedGridVelocityField(rCurrentProcessInfo, rModelPart);

                    #pragma omp parallel for
                    for (int i = 0; i < static_cast<int>(rElements.size()); ++i)
                    {
                        auto it_elem = it_elem_begin + i;
                        std::vector<bool> dummy;
                        it_elem->CalculateOnIntegrationPoints(CALCULATE_EXPLICIT_MP_STRESS, dummy, rCurrentProcessInfo);
                    }
                }

                // Finalizes solution step for all of the conditions
                const auto it_cond_begin = rModelPart.ConditionsBegin();
                #pragma omp parallel for
                for (int i = 0; i < static_cast<int>(rModelPart.Conditions().size()); ++i) {
                    auto it_cond = it_cond_begin + i;
                    it_cond->FinalizeSolutionStep(rCurrentProcessInfo);
                }

                KRATOS_CATCH("")
            }
            //***************************************************************************
            //***************************************************************************

            void CalculateUpdatedGridVelocityField(const ProcessInfo& rCurrentProcessInfo, ModelPart& rModelPart)
            {
                if (rCurrentProcessInfo.GetValue(EXPLICIT_STRESS_UPDATE_OPTION) == 1)
                {
                    // USL
                    const SizeType DomainSize = rCurrentProcessInfo[DOMAIN_SIZE];
                    #pragma omp parallel for
                    for (int iter = 0; iter < static_cast<int>(rModelPart.Nodes().size()); ++iter)
                    {
                        NodeIterator i = rModelPart.NodesBegin() + iter;
                        if ((i)->Is(ACTIVE))
                        {
                            array_1d<double, 3 >& r_nodal_momenta = (i)->FastGetSolutionStepValue(NODAL_MOMENTUM);
                            array_1d<double, 3>& r_current_velocity = i->FastGetSolutionStepValue(VELOCITY);
                            r_current_velocity.clear();
                            const double nodal_mass = i->FastGetSolutionStepValue(NODAL_MASS);
                            if (nodal_mass > numerical_limit)
                            {
                                for (IndexType j = 0; j < DomainSize; j++)
                                {
                                    r_current_velocity[j] = r_nodal_momenta[j] / nodal_mass;
                                } // for DomainSize
                            }
                        }
                    }
                }
                else if (rCurrentProcessInfo.GetValue(EXPLICIT_STRESS_UPDATE_OPTION) == 2)
                {
                    // MUSL stress update. This works by projecting the updated material point
                    // velocity back to the nodes. The nodal velocity field is then
                    // used for stress computations.

                    // Reset grid velocities
                    #pragma omp parallel for
                    for (int iter = 0; iter < static_cast<int>(rModelPart.Nodes().size()); ++iter)
                    {
                        auto i = rModelPart.NodesBegin() + iter;
                        array_1d<double, 3 >& nodal_velocity = (i)->FastGetSolutionStepValue(VELOCITY);
                        nodal_velocity.clear();
                    }

                    // Map updated MP velocities back to grid
                    const auto it_elem_begin = rModelPart.ElementsBegin();
                    #pragma omp parallel for
                    for (int i = 0; i < static_cast<int>(rModelPart.Elements().size()); ++i) {
                        auto it_elem = it_elem_begin + i;
                        std::vector<bool> dummy;
                        it_elem->CalculateOnIntegrationPoints(CALCULATE_MUSL_VELOCITY_FIELD, dummy, rCurrentProcessInfo);
                    }

                    // Reapply dirichlet BCs to MUSL velocity field
                    calculateGridVelocityAndApplyDirichletBC(rCurrentProcessInfo);
                }
            }

            //***************************************************************************
            //***************************************************************************

            /** Function that returns the list of Degrees of freedom to be
            assembled in the system for a Given Element
             */
            void GetDofList(
                const Element& rCurrentElement,
                Element::DofsVectorType& ElementalDofList,
                const ProcessInfo& CurrentProcessInfo) override
            {
                rCurrentElement.GetDofList(ElementalDofList, CurrentProcessInfo);
            }

            //***************************************************************************
            //***************************************************************************

            /** Function that returns the list of Degrees of freedom to be
            assembled in the system for a Given Element
             */
            void GetDofList(
                const Condition& rCurrentCondition,
                Element::DofsVectorType& rConditionDofList,
                const ProcessInfo& rCurrentProcessInfo) override
            {
                rCurrentCondition.GetDofList(rConditionDofList, rCurrentProcessInfo);
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
            int Check(const ModelPart& rModelPart) const override
            {
                KRATOS_TRY

                int err = Scheme<TSparseSpace, TDenseSpace>::Check(rModelPart);
                if (err != 0) return err;

                //check that variables are correctly allocated
                for (auto it = rModelPart.NodesBegin();
                    it != rModelPart.NodesEnd(); ++it)
                {
                    KRATOS_ERROR_IF(it->SolutionStepsDataHas(DISPLACEMENT) == false) << "DISPLACEMENT variable is not allocated for node " << it->Id() << std::endl;
                    KRATOS_ERROR_IF(it->SolutionStepsDataHas(VELOCITY) == false) << "VELOCITY variable is not allocated for node " << it->Id() << std::endl;
                    KRATOS_ERROR_IF(it->SolutionStepsDataHas(ACCELERATION) == false) << "ACCELERATION variable is not allocated for node " << it->Id() << std::endl;
                }

                //check that dofs exist
                for (auto it = rModelPart.NodesBegin();
                    it != rModelPart.NodesEnd(); ++it)
                {
                    KRATOS_ERROR_IF(it->HasDofFor(DISPLACEMENT_X) == false) << "Missing DISPLACEMENT_X dof on node " << it->Id() << std::endl;
                    KRATOS_ERROR_IF(it->HasDofFor(DISPLACEMENT_Y) == false) << "Missing DISPLACEMENT_Y dof on node " << it->Id() << std::endl;
                    KRATOS_ERROR_IF(it->HasDofFor(DISPLACEMENT_Z) == false) << "Missing DISPLACEMENT_Z dof on node " << it->Id() << std::endl;
                }

                //check for minimum value of the buffer index
                KRATOS_ERROR_IF(rModelPart.GetBufferSize() < 2) << "Insufficient buffer size. Buffer size should be greater than 2. Current size is" << rModelPart.GetBufferSize() << std::endl;
                KRATOS_CATCH("")
                return 0;
            }

            void CalculateRHSContribution(
                Element& rCurrentElement,
                LocalSystemVectorType& RHS_Contribution,
                Element::EquationIdVectorType& EquationId,
                const ProcessInfo& rCurrentProcessInfo
                ) override
            {
                KRATOS_TRY
                    rCurrentElement.CalculateRightHandSide(RHS_Contribution, rCurrentProcessInfo);
                    rCurrentElement.AddExplicitContribution(RHS_Contribution, RESIDUAL_VECTOR, FORCE_RESIDUAL, rCurrentProcessInfo);
                KRATOS_CATCH("")
            }

            /**
             * @brief Functions that calculates the RHS of a "condition" object
             * @param pCondition The condition to compute
             * @param RHS_Contribution The RHS vector contribution
             * @param EquationId The ID's of the condition degrees of freedom
             * @param rCurrentProcessInfo The current process info instance
             */
            void CalculateRHSContribution(
                Condition& rCurrentCondition,
                LocalSystemVectorType& RHS_Contribution,
                Element::EquationIdVectorType& EquationId,
                const ProcessInfo& rCurrentProcessInfo
                ) override
            {
                KRATOS_TRY
                rCurrentCondition.CalculateRightHandSide(RHS_Contribution, rCurrentProcessInfo);
                rCurrentCondition.AddExplicitContribution(RHS_Contribution, RESIDUAL_VECTOR, FORCE_RESIDUAL, rCurrentProcessInfo);
                KRATOS_CATCH("")
            }

        protected:
            /// @name Member Variables
            ModelPart& mr_grid_model_part;

        private:
    }; /* Class MPMExplicitScheme */
}  /* namespace Kratos.*/

#endif /* KRATOS_MPM_EXPLICIT_SCHEME defined */
