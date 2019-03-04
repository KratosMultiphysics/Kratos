//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//


#if !defined(KRATOS_RESIDUALBASED_PREDICTOR_CORRECTOR_VELOCITY_BOSSAK_TURBULENT_SCHEME )
#define  KRATOS_RESIDUALBASED_PREDICTOR_CORRECTOR_VELOCITY_BOSSAK_TURBULENT_SCHEME


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/deprecated_variables.h"
#include "solving_strategies/schemes/scheme.h"
#include "includes/variables.h"
#include "includes/cfd_variables.h"
#include "containers/array_1d.h"
#include "utilities/openmp_utils.h"
#include "utilities/dof_updater.h"
#include "utilities/coordinate_transformation_utilities.h"
#include "processes/process.h"
#include "solving_strategies/schemes/residual_based_bossak_velocity_scheme.h"
#include "utilities/derivatives_extension.h"

namespace Kratos {

    /**@name Kratos Globals */
    /*@{ */


    /*@} */
    /**@name Type Definitions */
    /*@{ */

    /*@} */


    /**@name  Enum's */
    /*@{ */


    /*@} */
    /**@name  Functions */
    /*@{ */



    /*@} */
    /**@name Kratos Classes */
    /*@{ */

    /// Bossak time scheme for the incompressible flow problem.
    /** This class provides a second order time scheme of the generalized-alpha Newmark
        family of methods. It also includes code required to implement slip conditions
        on the incompressible flow problem and provides the possibility of using a RANS
        model by passing a turbulence model as an argument to the constructor.
        This time scheme is intended to be used in combination with elements of type
        ASGS2D, ASGS3D, VMS or derived classes.
        To use the slip condition, assign IS_STRUCTURE != 0.0 to the non-historic database
        of the relevant boundary nodes (that is, use SetValue(IS_STRUCTURE,...)). To use
        a wall law in combination with the slip condition, use MonolithicWallCondition to
        mesh the boundary
        @see ASGS2D, ASGS3D, VMS, MonolithicWallConditon
     */
    template<class TSparseSpace,
    class TDenseSpace //= DenseSpace<double>
    >
    class ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent : public ResidualBasedBossakVelocityScheme<TSparseSpace, TDenseSpace> {
        class ElementDerivativesExtension : public DerivativesExtension
        {
            Element* mpElement;
        public:
            explicit ElementDerivativesExtension(Element* pElement): mpElement(pElement)
            {}

            void GetZeroDerivativesVector(std::size_t NodeId,
                                          std::vector<IndirectScalar<double>>& rVector,
                                          std::size_t Step,
                                          ProcessInfo& rCurrentProcessInfo) override
            {
                rVector.resize(3);
                std::size_t index = 0;
                Node<3>& r_node = mpElement->GetGeometry()[NodeId];
                rVector[index++] = MakeIndirectScalar(r_node, DISPLACEMENT_X, Step);
                rVector[index++] = MakeIndirectScalar(r_node, DISPLACEMENT_Y, Step);
                rVector[index] = MakeIndirectScalar(r_node, DISPLACEMENT_Z, Step);
            }

            void GetFirstDerivativesVector(std::size_t NodeId,
                                           std::vector<IndirectScalar<double>>& rVector,
                                           std::size_t Step,
                                           ProcessInfo& rCurrentProcessInfo) override
            {
                rVector.resize(4);
                std::size_t index = 0;
                Node<3>& r_node = mpElement->GetGeometry()[NodeId];
                rVector[index++] = MakeIndirectScalar(r_node, VELOCITY_X, Step);
                rVector[index++] = MakeIndirectScalar(r_node, VELOCITY_Y, Step);
                rVector[index++] = MakeIndirectScalar(r_node, VELOCITY_Z, Step);
                rVector[index] = MakeIndirectScalar(r_node, PRESSURE, Step);
            }

            void GetSecondDerivativesVector(std::size_t NodeId,
                                            std::vector<IndirectScalar<double>>& rVector,
                                            std::size_t Step,
                                            ProcessInfo& rCurrentProcessInfo) override
            {
                rVector.resize(3);
                std::size_t index = 0;
                Node<3>& r_node = mpElement->GetGeometry()[NodeId];
                rVector[index++] = MakeIndirectScalar(r_node, ACCELERATION_X, Step);
                rVector[index++] = MakeIndirectScalar(r_node, ACCELERATION_Y, Step);
                rVector[index] = MakeIndirectScalar(r_node, ACCELERATION_Z, Step);
            }

            void GetFirstDerivativesDofsVector(std::size_t NodeId,
                                               std::vector<Dof<double>::Pointer>& rVector,
                                               ProcessInfo& rCurrentProcessInfo) override
            {
                rVector.resize(4);
                std::size_t index = 0;
                Node<3>& r_node = mpElement->GetGeometry()[NodeId];
                rVector[index++] = r_node.pGetDof(VELOCITY_X);
                rVector[index++] = r_node.pGetDof(VELOCITY_Y);
                rVector[index++] = r_node.pGetDof(VELOCITY_Z);
                rVector[index++] = r_node.pGetDof(PRESSURE);
            }

            void GetZeroDerivativesVariables(std::vector<VariableData const*>& rVariables,
                                             ProcessInfo& rCurrentProcessInfo) const override
            {
                rVariables.resize(1);
                rVariables[0] = &DISPLACEMENT;
            }
            void GetFirstDerivativesVariables(std::vector<VariableData const*>& rVariables,
                                              ProcessInfo& rCurrentProcessInfo) const override
            {
                rVariables.resize(2);
                rVariables[0] = &VELOCITY;
                rVariables[1] = &PRESSURE;
            }

            void GetSecondDerivativesVariables(std::vector<VariableData const*>& rVariables,
                                               ProcessInfo& rCurrentProcessInfo) const override
            {
                rVariables.resize(1);
                rVariables[0] = &ACCELERATION;
            }
        };

    public:
        /**@name Type Definitions */
        /*@{ */

        KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent);

        typedef ResidualBasedBossakVelocityScheme<TSparseSpace, TDenseSpace> BaseType;

        typedef typename BaseType::TDataType TDataType;

        typedef typename BaseType::DofsArrayType DofsArrayType;

        typedef typename Element::DofsVectorType DofsVectorType;

        typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

        typedef typename BaseType::TSystemVectorType TSystemVectorType;

        typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

        typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

        typedef Element::GeometryType  GeometryType;


        /*@} */
        /**@name Life Cycle
         */
        /*@{ */

        /** Constructor without a turbulence model
         */
        ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent(
            double NewAlphaBossak,
            double MoveMeshStrategy,
            unsigned int DomainSize)
        :
          ResidualBasedBossakVelocityScheme<TSparseSpace, TDenseSpace>(NewAlphaBossak),
          mRotationTool(DomainSize,DomainSize+1,IS_STRUCTURE,0.0), // Second argument is number of matrix rows per node: monolithic elements have velocity and pressure dofs.
          mrPeriodicIdVar(Kratos::Variable<int>::StaticObject())
          {
            //default values for the Newmark Scheme
            mAlphaBossak = NewAlphaBossak;
            mMeshVelocity = MoveMeshStrategy;
        }


        /** Constructor without a turbulence model with periodic conditions
         */
        ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent(
            double NewAlphaBossak,
            unsigned int DomainSize,
            const Variable<int>& rPeriodicIdVar)
        :
          ResidualBasedBossakVelocityScheme<TSparseSpace, TDenseSpace>(NewAlphaBossak),
          mRotationTool(DomainSize,DomainSize+1,IS_STRUCTURE,0.0), // Second argument is number of matrix rows per node: monolithic elements have velocity and pressure dofs.
          mrPeriodicIdVar(rPeriodicIdVar)
          {
            //default values for the Newmark Scheme
            mAlphaBossak = NewAlphaBossak;
            mMeshVelocity = 0.0;
        }


        /** Constructor without a turbulence model
         */
        ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent(
            double NewAlphaBossak,
            double MoveMeshStrategy,
            unsigned int DomainSize,
            Variable<double>& rSlipVar)
        :
          ResidualBasedBossakVelocityScheme<TSparseSpace, TDenseSpace>(NewAlphaBossak),
          mRotationTool(DomainSize,DomainSize+1,rSlipVar,0.0), // Second argument is number of matrix rows per node: monolithic elements have velocity and pressure dofs.
          mrPeriodicIdVar(Kratos::Variable<int>::StaticObject())
          {
            //default values for the Newmark Scheme
            mAlphaBossak = NewAlphaBossak;
            mMeshVelocity = MoveMeshStrategy;
        }

        /** Constructor with a turbulence model
         */
        ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent(
            double NewAlphaBossak,
            double MoveMeshStrategy,
            unsigned int DomainSize,
            Process::Pointer pTurbulenceModel)
        :
          ResidualBasedBossakVelocityScheme<TSparseSpace, TDenseSpace>(NewAlphaBossak),
          mRotationTool(DomainSize,DomainSize+1,IS_STRUCTURE,0.0), // Second argument is number of matrix rows per node: monolithic elements have velocity and pressure dofs
          mrPeriodicIdVar(Kratos::Variable<int>::StaticObject()),
          mpTurbulenceModel(pTurbulenceModel)
          {
            //default values for the Newmark Scheme
            mAlphaBossak = NewAlphaBossak;
            mMeshVelocity = MoveMeshStrategy;
        }

        /** Destructor.
         */
        ~ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent() override {
        }


        /*@} */
        /**@name Operators
         */
        /*@{ */

        void Initialize(ModelPart& rModelPart) override
        {
            KRATOS_TRY;

            BaseType::Initialize(rModelPart);

            const int number_of_elements = rModelPart.NumberOfElements();

#pragma omp parallel for
            for (int i = 0; i < number_of_elements; i++)
            {
                Element& r_element = *(rModelPart.ElementsBegin() + i);
                r_element.SetValue(DERIVATIVES_EXTENSION,
                                        Kratos::make_shared<ElementDerivativesExtension>(&r_element));
            }

            KRATOS_CATCH("");
        }

        /**
                Performing the update of the solution.
         */
        //***************************************************************************

        void Update(ModelPart& r_model_part,
                            DofsArrayType& rDofSet,
                            TSystemMatrixType& A,
                            TSystemVectorType& Dv,
                            TSystemVectorType& b) override
        {
            KRATOS_TRY;

            mRotationTool.RotateVelocities(r_model_part);

            mpDofUpdater->UpdateDofs(rDofSet,Dv);

            mRotationTool.RecoverVelocities(r_model_part);

            BaseType::mUpdateDisplacement = (mMeshVelocity == 2);

            this->UpdateTimeSchemeVariables(r_model_part);

            AdditionalUpdateOperations(r_model_part, rDofSet, A, Dv, b);

            KRATOS_CATCH("")
        }

        //***************************************************************************

        void AdditionalUpdateOperations(ModelPart& rModelPart,
                                        DofsArrayType& rDofSet,
                                        TSystemMatrixType& A,
                                        TSystemVectorType& Dv,
                                        TSystemVectorType& b)
        {
            KRATOS_TRY

            int NumThreads = OpenMPUtils::GetNumThreads();
            OpenMPUtils::PartitionVector NodePartition;
            OpenMPUtils::DivideInPartitions(rModelPart.Nodes().size(), NumThreads, NodePartition);

            //updating time derivatives (nodally for efficiency)
            #pragma omp parallel
            {
                array_1d<double, 3 > DeltaVel;

                int k = OpenMPUtils::ThisThread();

                ModelPart::NodeIterator NodesBegin = rModelPart.NodesBegin() + NodePartition[k];
                ModelPart::NodeIterator NodesEnd = rModelPart.NodesBegin() + NodePartition[k + 1];

                for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; itNode++) {
                    if (mMeshVelocity == 2)//Lagrangian
                    {
                        if((itNode)->FastGetSolutionStepValue(IS_LAGRANGIAN_INLET) >= 1e-15)
                        {
                            noalias(itNode->FastGetSolutionStepValue(MESH_VELOCITY)) = ZeroVector(3);
                            noalias(itNode->FastGetSolutionStepValue(DISPLACEMENT)) = ZeroVector(3);
                        }
                        else
                        {
                            noalias(itNode->FastGetSolutionStepValue(MESH_VELOCITY)) = itNode->FastGetSolutionStepValue(VELOCITY);
                        }

                    }
                }
            }

            KRATOS_CATCH("")

        }


        //***************************************************************************
        //predicts the solution at the current step as
        // v = vold
        void Predict(ModelPart& rModelPart,
                     DofsArrayType& rDofSet,
                     TSystemMatrixType& A,
                     TSystemVectorType& Dv,
                     TSystemVectorType& b) override
        {
            BaseType::mUpdateDisplacement = (mMeshVelocity == 2);

            BaseType::Predict(rModelPart, rDofSet, A, Dv, b);

            int NumThreads = OpenMPUtils::GetNumThreads();
            OpenMPUtils::PartitionVector NodePartition;
            OpenMPUtils::DivideInPartitions(rModelPart.Nodes().size(), NumThreads, NodePartition);
#pragma omp parallel
            {
                int k = OpenMPUtils::ThisThread();

                ModelPart::NodeIterator NodesBegin =
                    rModelPart.NodesBegin() + NodePartition[k];
                ModelPart::NodeIterator NodesEnd =
                    rModelPart.NodesBegin() + NodePartition[k + 1];

                for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; itNode++)
                {
                    if (mMeshVelocity == 2) // Lagrangian
                    {
                        if ((itNode)->FastGetSolutionStepValue(IS_LAGRANGIAN_INLET) < 1e-15)
                        {
                            noalias(itNode->FastGetSolutionStepValue(MESH_VELOCITY)) =
                                itNode->FastGetSolutionStepValue(VELOCITY);
                        }
                        else
                        {
                            itNode->FastGetSolutionStepValue(MESH_VELOCITY_X) = 0.0;
                            itNode->FastGetSolutionStepValue(MESH_VELOCITY_Y) = 0.0;
                            itNode->FastGetSolutionStepValue(DISPLACEMENT_X) = 0.0;
                            itNode->FastGetSolutionStepValue(DISPLACEMENT_Y) = 0.0;
                        }
                    }
                }
            }

        }

        //***************************************************************************

        /** this function is designed to be called in the builder and solver
        to introduce
        the selected time integration scheme. It "asks" the matrix needed to
        the element and
        performs the operations needed to introduce the seected time
        integration scheme.

          this function calculates at the same time the contribution to the
        LHS and to the RHS
          of the system
         */
        void CalculateSystemContributions(Element::Pointer rCurrentElement,
                                          LocalSystemMatrixType& LHS_Contribution,
                                          LocalSystemVectorType& RHS_Contribution,
                                          Element::EquationIdVectorType& EquationId,
                                          ProcessInfo& CurrentProcessInfo) override
        {
            KRATOS_TRY

            BaseType::CalculateSystemContributions(rCurrentElement, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

            // If there is a slip condition, apply it on a rotated system of coordinates
            mRotationTool.Rotate(LHS_Contribution,RHS_Contribution,rCurrentElement->GetGeometry());
            mRotationTool.ApplySlipCondition(LHS_Contribution,RHS_Contribution,rCurrentElement->GetGeometry());

            KRATOS_CATCH("")
        }

        void Calculate_RHS_Contribution(Element::Pointer rCurrentElement,
                                        LocalSystemVectorType& RHS_Contribution,
                                        Element::EquationIdVectorType& EquationId,
                                        ProcessInfo& CurrentProcessInfo) override
        {
            BaseType::Calculate_RHS_Contribution(rCurrentElement, RHS_Contribution, EquationId, CurrentProcessInfo);

            // If there is a slip condition, apply it on a rotated system of coordinates
            mRotationTool.Rotate(RHS_Contribution,rCurrentElement->GetGeometry());
            mRotationTool.ApplySlipCondition(RHS_Contribution,rCurrentElement->GetGeometry());
        }

        /** functions totally analogous to the precedent but applied to
        the "condition" objects
         */
        void Condition_CalculateSystemContributions(Condition::Pointer rCurrentCondition,
                                                            LocalSystemMatrixType& LHS_Contribution,
                                                            LocalSystemVectorType& RHS_Contribution,
                                                            Element::EquationIdVectorType& EquationId,
                                                            ProcessInfo& CurrentProcessInfo) override
        {
            KRATOS_TRY

            BaseType::Condition_CalculateSystemContributions(rCurrentCondition, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

            // Rotate contributions (to match coordinates for slip conditions)
            mRotationTool.Rotate(LHS_Contribution,RHS_Contribution,rCurrentCondition->GetGeometry());
            mRotationTool.ApplySlipCondition(LHS_Contribution,RHS_Contribution,rCurrentCondition->GetGeometry());

            KRATOS_CATCH("")
        }

        void Condition_Calculate_RHS_Contribution(Condition::Pointer rCurrentCondition,
                                                          LocalSystemVectorType& RHS_Contribution,
                                                          Element::EquationIdVectorType& EquationId,
                                                          ProcessInfo& rCurrentProcessInfo) override
        {
            KRATOS_TRY;

            BaseType::Condition_Calculate_RHS_Contribution(rCurrentCondition, RHS_Contribution, EquationId, rCurrentProcessInfo);

            // Rotate contributions (to match coordinates for slip conditions)
            mRotationTool.Rotate(RHS_Contribution,rCurrentCondition->GetGeometry());
            mRotationTool.ApplySlipCondition(RHS_Contribution,rCurrentCondition->GetGeometry());

            KRATOS_CATCH("");
        }
        //*************************************************************************************
        //*************************************************************************************

        void InitializeSolutionStep(ModelPart& r_model_part,
                                            TSystemMatrixType& A,
                                            TSystemVectorType& Dx,
                                            TSystemVectorType& b) override
        {
            BaseType::InitializeSolutionStep(r_model_part, A, Dx, b);
        }
        //*************************************************************************************
        //*************************************************************************************

        void InitializeNonLinIteration(ModelPart& r_model_part,
                                               TSystemMatrixType& A,
                                               TSystemVectorType& Dx,
                                               TSystemVectorType& b) override
        {

        }

        void FinalizeNonLinIteration(ModelPart &rModelPart, TSystemMatrixType &A, TSystemVectorType &Dx, TSystemVectorType &b) override
        {

            KRATOS_TRY

            ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

            //if orthogonal subscales are computed
            if (CurrentProcessInfo[OSS_SWITCH] == 1.0) {

                KRATOS_INFO_IF("ResidualBasedSimpleSteadyScheme", rModelPart.GetCommunicator().MyPID() == 0)
                    << "Computing OSS projections" << std::endl;

                const int nnodes = static_cast<int>(rModelPart.Nodes().size());
                auto nbegin = rModelPart.NodesBegin();
                #pragma omp parallel for firstprivate(nbegin,nnodes)
                for(int i=0; i<nnodes; ++i)
                {
                    auto ind = nbegin + i;
                    noalias(ind->FastGetSolutionStepValue(ADVPROJ)) = ZeroVector(3);

                    ind->FastGetSolutionStepValue(DIVPROJ) = 0.0;

                    ind->FastGetSolutionStepValue(NODAL_AREA) = 0.0;


                }//end of loop over nodes

                //loop on nodes to compute ADVPROJ   CONVPROJ NODALAREA
                array_1d<double, 3 > output = ZeroVector(3);

                const int nel = static_cast<int>(rModelPart.Elements().size());
                auto elbegin = rModelPart.ElementsBegin();
                #pragma omp parallel for firstprivate(elbegin,nel,output)
                for(int i=0; i<nel; ++i)
                {
                    auto elem = elbegin + i;
                    elem->Calculate(ADVPROJ, output, CurrentProcessInfo);
                }

                rModelPart.GetCommunicator().AssembleCurrentData(NODAL_AREA);
                rModelPart.GetCommunicator().AssembleCurrentData(DIVPROJ);
                rModelPart.GetCommunicator().AssembleCurrentData(ADVPROJ);

                // Correction for periodic conditions
                this->PeriodicConditionProjectionCorrection(rModelPart);

                #pragma omp parallel for firstprivate(nbegin,nnodes)
                for(int i=0; i<nnodes; ++i)
                {
                    auto ind = nbegin + i;
                    if (ind->FastGetSolutionStepValue(NODAL_AREA) == 0.0)
                    {
                        ind->FastGetSolutionStepValue(NODAL_AREA) = 1.0;
                        //KRATOS_WATCH("*********ATTENTION: NODAL AREA IS ZERRROOOO************");
                    }
                    const double Area = ind->FastGetSolutionStepValue(NODAL_AREA);
                    ind->FastGetSolutionStepValue(ADVPROJ) /= Area;
                    ind->FastGetSolutionStepValue(DIVPROJ) /= Area;
                }
            }

            if (mpTurbulenceModel != 0) // If not null
                mpTurbulenceModel->Execute();

            KRATOS_CATCH("")

        }

        void FinalizeSolutionStep(ModelPart &rModelPart, TSystemMatrixType &A, TSystemVectorType &Dx, TSystemVectorType &b) override
        {
            Element::EquationIdVectorType EquationId;
            LocalSystemVectorType RHS_Contribution;
            LocalSystemMatrixType LHS_Contribution;
            ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

            //for (ModelPart::NodeIterator itNode = rModelPart.NodesBegin(); itNode != rModelPart.NodesEnd(); ++itNode)

            #pragma omp parallel for
            for(int k = 0; k<static_cast<int>(rModelPart.Nodes().size()); k++)
            {
                auto itNode = rModelPart.NodesBegin() + k;
                (itNode->FastGetSolutionStepValue(REACTION)).clear();

                // calculating relaxed acceleration
                const array_1d<double, 3 > & CurrentAcceleration = (itNode)->FastGetSolutionStepValue(ACCELERATION, 0);
                const array_1d<double, 3 > & OldAcceleration = (itNode)->FastGetSolutionStepValue(ACCELERATION, 1);
                const array_1d<double, 3> relaxed_acceleration = (1 - mAlphaBossak) * CurrentAcceleration
                                                                    + mAlphaBossak * OldAcceleration;
                (itNode)->SetValue(RELAXED_ACCELERATION, relaxed_acceleration);
            }

            //for (ModelPart::ElementsContainerType::ptr_iterator itElem = rModelPart.Elements().ptr_begin(); itElem != rModelPart.Elements().ptr_end(); ++itElem)

            #pragma omp parallel for firstprivate(EquationId,RHS_Contribution,LHS_Contribution)
            for(int k = 0; k<static_cast<int>(rModelPart.Elements().size()); k++)
            {
                auto itElem = rModelPart.Elements().ptr_begin()+k;
                int thread_id = OpenMPUtils::ThisThread();

                (*itElem)->InitializeNonLinearIteration(CurrentProcessInfo);
                //KRATOS_WATCH(LHS_Contribution);
                //basic operations for the element considered
                (*itElem)->CalculateLocalSystem(LHS_Contribution, RHS_Contribution, CurrentProcessInfo);

                //std::cout << rCurrentElement->Id() << " RHS = " << RHS_Contribution << std::endl;
                (*itElem)->CalculateMassMatrix(this->mMassMatrix[thread_id], CurrentProcessInfo);
                (*itElem)->CalculateLocalVelocityContribution(this->mDampingMatrix[thread_id], RHS_Contribution, CurrentProcessInfo);

                (*itElem)->EquationIdVector(EquationId, CurrentProcessInfo);

                //adding the dynamic contributions (statics is already included)
                this->AddDynamicsToLHS(LHS_Contribution, this->mDampingMatrix[thread_id], this->mMassMatrix[thread_id], CurrentProcessInfo);
                this->AddDynamicsToRHS((*itElem), RHS_Contribution, this->mDampingMatrix[thread_id], this->mMassMatrix[thread_id], CurrentProcessInfo);

                GeometryType& rGeom = (*itElem)->GetGeometry();
                unsigned int NumNodes = rGeom.PointsNumber();
                unsigned int Dimension = rGeom.WorkingSpaceDimension();

                unsigned int index = 0;
                for (unsigned int i = 0; i < NumNodes; i++)
                {
                    auto& reaction = rGeom[i].FastGetSolutionStepValue(REACTION);

                    double& target_value0 = reaction[0];
                    const double& origin_value0 = RHS_Contribution[index++];
                    #pragma omp atomic
                    target_value0 -= origin_value0;

                    double& target_value1 = reaction[1];
                    const double& origin_value1 = RHS_Contribution[index++];
                    #pragma omp atomic
                    target_value1 -= origin_value1;

                    if (Dimension == 3)
                    {
                      double& target_value2 = reaction[2];
                      const double& origin_value2 = RHS_Contribution[index++];
                      #pragma omp atomic
                      target_value2 -= origin_value2;
                    }
            //        rGeom[i].FastGetSolutionStepValue(REACTION_X,0) -= RHS_Contribution[index++];
             //          rGeom[i].FastGetSolutionStepValue(REACTION_Y,0) -= RHS_Contribution[index++];
            //        if (Dimension == 3) rGeom[i].FastGetSolutionStepValue(REACTION_Z,0) -= RHS_Contribution[index++];
                    index++; // skip pressure dof
                }
            }

            rModelPart.GetCommunicator().AssembleCurrentData(REACTION);

            // Base scheme calls FinalizeSolutionStep method of elements and conditions
            Scheme<TSparseSpace, TDenseSpace>::FinalizeSolutionStep(rModelPart, A, Dx, b);
        }

        //************************************************************************************************
        //************************************************************************************************

        /// Free memory allocated by this object.
        void Clear() override
        {
            this->mpDofUpdater->Clear();
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


        /*@} */

    protected:
        /**@name Protected static Member Variables */
        /*@{ */


        /*@} */
        /**@name Protected member Variables */
        /*@{ */
        double mAlphaBossak;
        double mMeshVelocity;


        /*@} */
        /**@name Protected Operators*/
        /*@{ */

        /** On periodic boundaries, the nodal area and the values to project need to take into account contributions from elements on
         * both sides of the boundary. This is done using the conditions and the non-historical nodal data containers as follows:\n
         * 1- The partition that owns the PeriodicCondition adds the values on both nodes to their non-historical containers.\n
         * 2- The non-historical containers are added across processes, communicating the right value from the condition owner to all partitions.\n
         * 3- The value on all periodic nodes is replaced by the one received in step 2.
         */
        void PeriodicConditionProjectionCorrection(ModelPart& rModelPart)
        {
            const int num_nodes = rModelPart.NumberOfNodes();
            const int num_conditions = rModelPart.NumberOfConditions();

            #pragma omp parallel for
            for (int i = 0; i < num_nodes; i++) {
                auto it_node = rModelPart.NodesBegin() + i;

                it_node->SetValue(NODAL_AREA,0.0);
                it_node->SetValue(ADVPROJ,ZeroVector(3));
                it_node->SetValue(DIVPROJ,0.0);
            }

            #pragma omp parallel for
            for (int i = 0; i < num_conditions; i++) {
                auto it_cond = rModelPart.ConditionsBegin() + i;

                if(it_cond->Is(PERIODIC)) {
                    this->AssemblePeriodicContributionToProjections(it_cond->GetGeometry());
                }
            }

            rModelPart.GetCommunicator().AssembleNonHistoricalData(NODAL_AREA);
            rModelPart.GetCommunicator().AssembleNonHistoricalData(ADVPROJ);
            rModelPart.GetCommunicator().AssembleNonHistoricalData(DIVPROJ);

            #pragma omp parallel for
            for (int i = 0; i < num_nodes; i++) {
                auto it_node = rModelPart.NodesBegin() + i;
                this->CorrectContributionsOnPeriodicNode(*it_node);
            }
        }

        void AssemblePeriodicContributionToProjections(Geometry< Node<3> >& rGeometry)
        {
            unsigned int nodes_in_cond = rGeometry.PointsNumber();

            double nodal_area = 0.0;
            array_1d<double,3> momentum_projection = ZeroVector(3);
            double mass_projection = 0.0;
            for ( unsigned int i = 0; i < nodes_in_cond; i++ )
            {
                auto& r_node = rGeometry[i];
                nodal_area += r_node.FastGetSolutionStepValue(NODAL_AREA);
                noalias(momentum_projection) += r_node.FastGetSolutionStepValue(ADVPROJ);
                mass_projection += r_node.FastGetSolutionStepValue(DIVPROJ);
            }

            for ( unsigned int i = 0; i < nodes_in_cond; i++ )
            {
                auto& r_node = rGeometry[i];
                /* Note that this loop is expected to be threadsafe in normal conditions,
                * since each node should belong to a single periodic link. However, I am
                * setting the locks for openmp in case that we try more complicated things
                * in the future (like having different periodic conditions for different
                * coordinate directions).
                */
                r_node.SetLock();
                r_node.GetValue(NODAL_AREA) = nodal_area;
                noalias(r_node.GetValue(ADVPROJ)) = momentum_projection;
                r_node.GetValue(DIVPROJ) = mass_projection;
                r_node.UnSetLock();
            }
        }

        void CorrectContributionsOnPeriodicNode(Node<3>& rNode)
        {
            if (rNode.GetValue(NODAL_AREA) != 0.0) // Only periodic nodes will have a non-historical NODAL_AREA set.
            {
                rNode.FastGetSolutionStepValue(NODAL_AREA) = rNode.GetValue(NODAL_AREA);
                noalias(rNode.FastGetSolutionStepValue(ADVPROJ)) = rNode.GetValue(ADVPROJ);
                rNode.FastGetSolutionStepValue(DIVPROJ) = rNode.GetValue(DIVPROJ);
            }
        }

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



        /*@} */

    private:
        /**@name Static Member Variables */
        /*@{ */


        /*@} */
        /**@name Member Variables */
        /*@{ */

        CoordinateTransformationUtils<LocalSystemMatrixType,LocalSystemVectorType,double> mRotationTool;

        const Variable<int>& mrPeriodicIdVar;

        Process::Pointer mpTurbulenceModel;

        typename TSparseSpace::DofUpdaterPointerType mpDofUpdater = TSparseSpace::CreateDofUpdater();

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
        /**@name Un accessible methods */
        /*@{ */


        /*@} */

    }; /* Class Scheme */

    /*@} */

    /**@name Type Definitions */
    /*@{ */


    /*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_RESIDUALBASED_PREDICTOR_CORRECTOR_BOSSAK_SCHEME  defined */
