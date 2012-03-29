/*
==============================================================================
KratosFluidDynamicsApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel 
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
 */

#if !defined(KRATOS_RESIDUALBASED_PREDICTOR_CORRECTOR_VELOCITY_BOSSAK_TURBULENT_SCHEME )
#define  KRATOS_RESIDUALBASED_PREDICTOR_CORRECTOR_VELOCITY_BOSSAK_TURBULENT_SCHEME


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"
#include "boost/numeric/ublas/matrix_proxy.hpp"

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "includes/variables.h"
#include "containers/array_1d.h"
#include "utilities/openmp_utils.h"
#include "processes/process.h"

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
    class ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent : public Scheme<TSparseSpace, TDenseSpace> {
    public:
        /**@name Type Definitions */
        /*@{ */

        //typedef boost::shared_ptr< ResidualBasedPredictorCorrectorBossakScheme<TSparseSpace,TDenseSpace> > Pointer;

        KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent);

        typedef Scheme<TSparseSpace, TDenseSpace> BaseType;

        typedef typename BaseType::TDataType TDataType;

        typedef typename BaseType::DofsArrayType DofsArrayType;

        typedef typename Element::DofsVectorType DofsVectorType;

        typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

        typedef typename BaseType::TSystemVectorType TSystemVectorType;

        typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

        typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

        typedef boost::numeric::ublas::matrix_row<LocalSystemMatrixType>  LocalRowType;

        typedef boost::numeric::ublas::matrix_range<LocalSystemMatrixType> MatrixBlockType;

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
          Scheme<TSparseSpace, TDenseSpace>(),
          mDomainSize(DomainSize)
        {
            //default values for the Newmark Scheme
            mAlphaBossak = NewAlphaBossak;
            mBetaNewmark = 0.25 * pow((1.00 - mAlphaBossak), 2);
            mGammaNewmark = 0.5 - mAlphaBossak;
            mMeshVelocity = MoveMeshStrategy;


            //Allocate auxiliary memory
            int NumThreads = OpenMPUtils::GetNumThreads();
            mMass.resize(NumThreads);
            mDamp.resize(NumThreads);
            mvel.resize(NumThreads);
            macc.resize(NumThreads);
            maccold.resize(NumThreads);

            std::cout << "using the velocity Bossak Time Integration Scheme" << std::endl;
        }

        /** Constructor with a turbulence model
         */
        ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent(
            double NewAlphaBossak,
            double MoveMeshStrategy,
            unsigned int DomainSize,
            Process::Pointer pTurbulenceModel)
        :
          Scheme<TSparseSpace, TDenseSpace>(),
          mDomainSize(DomainSize),
          mpTurbulenceModel(pTurbulenceModel)
        {
            //default values for the Newmark Scheme
            mAlphaBossak = NewAlphaBossak;
            mBetaNewmark = 0.25 * pow((1.00 - mAlphaBossak), 2);
            mGammaNewmark = 0.5 - mAlphaBossak;
            mMeshVelocity = MoveMeshStrategy;


            //Allocate auxiliary memory
            int NumThreads = OpenMPUtils::GetNumThreads();
            mMass.resize(NumThreads);
            mDamp.resize(NumThreads);
            mvel.resize(NumThreads);
            macc.resize(NumThreads);
            maccold.resize(NumThreads);

            std::cout << "using the velocity Bossak Time Integration Scheme" << std::endl;
        }

        /** Destructor.
         */
        virtual ~ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent() {
        }


        /*@} */
        /**@name Operators
         */
        /*@{ */

        /**
                Performing the update of the solution.
         */
        //***************************************************************************

        virtual void Update(ModelPart& r_model_part,
                            DofsArrayType& rDofSet,
                            TSystemMatrixType& A,
                            TSystemVectorType& Dv,
                            TSystemVectorType& b)
        {
            KRATOS_TRY;

            this->RotateVelocities(r_model_part);

            BasicUpdateOperations(r_model_part, rDofSet, A, Dv, b);

            this->RecoverVelocities(r_model_part);

            AdditionalUpdateOperations(r_model_part, rDofSet, A, Dv, b);

            KRATOS_CATCH("")
        }

        //***************************************************************************

        virtual void BasicUpdateOperations(ModelPart& rModelPart,
                                           DofsArrayType& rDofSet,
                                           TSystemMatrixType& A,
                                           TSystemVectorType& Dv,
                                           TSystemVectorType& b)
        {
            KRATOS_TRY

            int NumThreads = OpenMPUtils::GetNumThreads();
            OpenMPUtils::PartitionVector DofSetPartition;
            OpenMPUtils::DivideInPartitions(rDofSet.size(), NumThreads, DofSetPartition);

            //update of velocity (by DOF)
            #pragma omp parallel
            {
                int k = OpenMPUtils::ThisThread();

                typename DofsArrayType::iterator DofSetBegin = rDofSet.begin() + DofSetPartition[k];
                typename DofsArrayType::iterator DofSetEnd = rDofSet.begin() + DofSetPartition[k + 1];

                for (typename DofsArrayType::iterator itDof = DofSetBegin; itDof != DofSetEnd; itDof++) {
                    if (itDof->IsFree()) {
                        itDof->GetSolutionStepValue() += TSparseSpace::GetValue(Dv, itDof->EquationId());
                    }
                }

            }

            KRATOS_CATCH("")
        }

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
                    noalias(DeltaVel) = (itNode)->FastGetSolutionStepValue(VELOCITY) - (itNode)->FastGetSolutionStepValue(VELOCITY, 1);

                    array_1d<double, 3 > & CurrentAcceleration = (itNode)->FastGetSolutionStepValue(ACCELERATION, 0);
                    array_1d<double, 3 > & OldAcceleration = (itNode)->FastGetSolutionStepValue(ACCELERATION, 1);

                    UpdateAcceleration(CurrentAcceleration, DeltaVel, OldAcceleration);

                    if (mMeshVelocity == 2)//Lagrangian
                    {
                        array_1d<double, 3 > & CurrentDisplacement = (itNode)->FastGetSolutionStepValue(DISPLACEMENT, 0);
                        array_1d<double, 3 > & OldDisplacement = (itNode)->FastGetSolutionStepValue(DISPLACEMENT, 1);

                        array_1d<double, 3 > & OldVelocity = (itNode)->FastGetSolutionStepValue(VELOCITY, 1);

                        noalias(itNode->FastGetSolutionStepValue(MESH_VELOCITY)) = itNode->FastGetSolutionStepValue(VELOCITY);
                        UpdateDisplacement(CurrentDisplacement, OldDisplacement, OldVelocity, OldAcceleration, CurrentAcceleration);
                    }
                }
            }

            KRATOS_CATCH("")

        }


        //***************************************************************************
        //predicts the solution at the current step as
        // v = vold
        virtual void Predict(ModelPart& rModelPart,
                             DofsArrayType& rDofSet,
                             TSystemMatrixType& A,
                             TSystemVectorType& Dv,
                             TSystemVectorType& b)
        {
            std::cout << "prediction" << std::endl;

            int NumThreads = OpenMPUtils::GetNumThreads();
            OpenMPUtils::PartitionVector NodePartition;
            OpenMPUtils::DivideInPartitions(rModelPart.Nodes().size(), NumThreads, NodePartition);

            #pragma omp parallel
            {
                //array_1d<double, 3 > DeltaDisp;

                int k = OpenMPUtils::ThisThread();

                ModelPart::NodeIterator NodesBegin = rModelPart.NodesBegin() + NodePartition[k];
                ModelPart::NodeIterator NodesEnd = rModelPart.NodesBegin() + NodePartition[k + 1];

                for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; itNode++) {
                    array_1d<double, 3 > & OldVelocity = (itNode)->FastGetSolutionStepValue(VELOCITY, 1);
                    double& OldPressure = (itNode)->FastGetSolutionStepValue(PRESSURE, 1);
                    double& OldAirPressure = (itNode)->FastGetSolutionStepValue(AIR_PRESSURE, 1);

                    //predicting velocity
                    //ATTENTION::: the prediction is performed only on free nodes
                    array_1d<double, 3 > & CurrentVelocity = (itNode)->FastGetSolutionStepValue(VELOCITY);
                    double& CurrentPressure = (itNode)->FastGetSolutionStepValue(PRESSURE);
                    double& CurrentAirPressure = (itNode)->FastGetSolutionStepValue(AIR_PRESSURE);

                    if ((itNode->pGetDof(VELOCITY_X))->IsFree())
                        (CurrentVelocity[0]) = OldVelocity[0];
                    if (itNode->pGetDof(VELOCITY_Y)->IsFree())
                        (CurrentVelocity[1]) = OldVelocity[1];
                    if (itNode->HasDofFor(VELOCITY_Z))
                        if (itNode->pGetDof(VELOCITY_Z)->IsFree())
                            (CurrentVelocity[2]) = OldVelocity[2];

                    if (itNode->pGetDof(PRESSURE)->IsFree())
                        CurrentPressure = OldPressure;
                    if (itNode->HasDofFor(AIR_PRESSURE))
                        if (itNode->pGetDof(AIR_PRESSURE)->IsFree())
                            CurrentAirPressure = OldAirPressure;

                    // updating time derivatives ::: please note that displacements and
                    // their time derivatives can not be consistently fixed separately
                    array_1d<double, 3 > DeltaVel;
                    noalias(DeltaVel) = CurrentVelocity - OldVelocity;
                    array_1d<double, 3 > & OldAcceleration = (itNode)->FastGetSolutionStepValue(ACCELERATION, 1);
                    array_1d<double, 3 > & CurrentAcceleration = (itNode)->FastGetSolutionStepValue(ACCELERATION);

                    UpdateAcceleration(CurrentAcceleration, DeltaVel, OldAcceleration);

                    if (mMeshVelocity == 2) //Lagrangian
                    {
                        array_1d<double, 3 > & OldDisplacement = (itNode)->FastGetSolutionStepValue(DISPLACEMENT, 1);
                        array_1d<double, 3 > & CurrentDisplacement = (itNode)->FastGetSolutionStepValue(DISPLACEMENT, 0);

                        noalias(itNode->FastGetSolutionStepValue(MESH_VELOCITY)) = itNode->FastGetSolutionStepValue(VELOCITY);
                        UpdateDisplacement(CurrentDisplacement, OldDisplacement, OldVelocity, OldAcceleration, CurrentAcceleration);
                    }
                }
            }

            std::cout << "end of prediction" << std::endl;

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
                                          ProcessInfo& CurrentProcessInfo)
        {
            KRATOS_TRY
            int k = OpenMPUtils::ThisThread();

            //Initializing the non linear iteration for the current element
            (rCurrentElement) -> InitializeNonLinearIteration(CurrentProcessInfo);
            //KRATOS_WATCH(LHS_Contribution);
            //basic operations for the element considered
            (rCurrentElement)->CalculateLocalSystem(LHS_Contribution, RHS_Contribution, CurrentProcessInfo);

            //std::cout << rCurrentElement->Id() << " RHS = " << RHS_Contribution << std::endl;
            (rCurrentElement)->MassMatrix(mMass[k], CurrentProcessInfo);
            (rCurrentElement)->CalculateLocalVelocityContribution(mDamp[k], RHS_Contribution, CurrentProcessInfo);

            (rCurrentElement)->EquationIdVector(EquationId, CurrentProcessInfo);

            //adding the dynamic contributions (statics is already included)

            AddDynamicsToLHS(LHS_Contribution, mDamp[k], mMass[k], CurrentProcessInfo);
            AddDynamicsToRHS(rCurrentElement, RHS_Contribution, mDamp[k], mMass[k], CurrentProcessInfo);

            // If there is a slip condition, apply it on a rotated system of coordinates
            this->Rotate(LHS_Contribution,RHS_Contribution,rCurrentElement->GetGeometry());
            this->ApplySlipCondition(LHS_Contribution,RHS_Contribution,rCurrentElement->GetGeometry());

            KRATOS_CATCH("")
        }

        void Calculate_RHS_Contribution(Element::Pointer rCurrentElement,
                                        LocalSystemVectorType& RHS_Contribution,
                                        Element::EquationIdVectorType& EquationId,
                                        ProcessInfo& CurrentProcessInfo)
        {
            int k = OpenMPUtils::ThisThread();

            //Initializing the non linear iteration for the current element
            (rCurrentElement) -> InitializeNonLinearIteration(CurrentProcessInfo);

            //basic operations for the element considered
            (rCurrentElement)->CalculateRightHandSide(RHS_Contribution, CurrentProcessInfo);
            (rCurrentElement)->MassMatrix(mMass[k], CurrentProcessInfo);

            (rCurrentElement)->CalculateLocalVelocityContribution(mDamp[k], RHS_Contribution, CurrentProcessInfo);

            (rCurrentElement)->EquationIdVector(EquationId, CurrentProcessInfo);

            //adding the dynamic contributions (static is already included)

            AddDynamicsToRHS(rCurrentElement, RHS_Contribution, mDamp[k], mMass[k], CurrentProcessInfo);

            // If there is a slip condition, apply it on a rotated system of coordinates
            this->Rotate(RHS_Contribution,rCurrentElement->GetGeometry());
            this->ApplySlipCondition(RHS_Contribution,rCurrentElement->GetGeometry());
        }

        /** functions totally analogous to the precedent but applied to
        the "condition" objects
         */
        virtual void Condition_CalculateSystemContributions(Condition::Pointer rCurrentCondition,
                                                            LocalSystemMatrixType& LHS_Contribution,
                                                            LocalSystemVectorType& RHS_Contribution,
                                                            Element::EquationIdVectorType& EquationId,
                                                            ProcessInfo& CurrentProcessInfo)
        {
            KRATOS_TRY
            int k = OpenMPUtils::ThisThread();

            //KRATOS_WATCH("CONDITION LOCALVELOCITYCONTRIBUTION IS NOT DEFINED");
            (rCurrentCondition) -> InitializeNonLinearIteration(CurrentProcessInfo);
            (rCurrentCondition)->CalculateLocalSystem(LHS_Contribution, RHS_Contribution, CurrentProcessInfo);
            (rCurrentCondition)->MassMatrix(mMass[k], CurrentProcessInfo);
            //(rCurrentCondition)->DampMatrix(VelocityBossakAuxiliaries::mDamp,CurrentProcessInfo);
            (rCurrentCondition)->CalculateLocalVelocityContribution(mDamp[k], RHS_Contribution, CurrentProcessInfo);
            (rCurrentCondition)->EquationIdVector(EquationId, CurrentProcessInfo);


            AddDynamicsToLHS(LHS_Contribution, mDamp[k], mMass[k], CurrentProcessInfo);

            AddDynamicsToRHS(rCurrentCondition, RHS_Contribution, mDamp[k], mMass[k], CurrentProcessInfo);

            // Rotate contributions (to match coordinates for slip conditions)
            this->Rotate(LHS_Contribution,RHS_Contribution,rCurrentCondition->GetGeometry());
            this->ApplySlipCondition(LHS_Contribution,RHS_Contribution,rCurrentCondition->GetGeometry());

            KRATOS_CATCH("")
        }

        virtual void Condition_Calculate_RHS_Contribution(Condition::Pointer rCurrentCondition,
                                                          LocalSystemVectorType& RHS_Contribution,
                                                          Element::EquationIdVectorType& EquationId,
                                                          ProcessInfo& rCurrentProcessInfo)
        {
            KRATOS_TRY;

            int k = OpenMPUtils::ThisThread();

            //KRATOS_WATCH("CONDITION LOCALVELOCITYCONTRIBUTION IS NOT DEFINED");
            //Initializing the non linear iteration for the current condition
            (rCurrentCondition) -> InitializeNonLinearIteration(rCurrentProcessInfo);

            //basic operations for the element considered
            (rCurrentCondition)->CalculateRightHandSide(RHS_Contribution,rCurrentProcessInfo);
            (rCurrentCondition)->MassMatrix(mMass[k],rCurrentProcessInfo);
            //(rCurrentCondition)->DampMatrix(VelocityBossakAuxiliaries::mDamp,CurrentProcessInfo);
            (rCurrentCondition)->CalculateLocalVelocityContribution(mDamp[k], RHS_Contribution,rCurrentProcessInfo);
            (rCurrentCondition)->EquationIdVector(EquationId,rCurrentProcessInfo);

            //adding the dynamic contributions (static is already included)
            AddDynamicsToRHS(rCurrentCondition, RHS_Contribution, mDamp[k], mMass[k],rCurrentProcessInfo);

            // Rotate contributions (to match coordinates for slip conditions)
            this->Rotate(RHS_Contribution,rCurrentCondition->GetGeometry());
            this->ApplySlipCondition(RHS_Contribution,rCurrentCondition->GetGeometry());

            KRATOS_CATCH("");
        }
        //*************************************************************************************
        //*************************************************************************************

        void InitializeSolutionStep(ModelPart& r_model_part,
                                    TSystemMatrixType& A,
                                    TSystemVectorType& Dx,
                                    TSystemVectorType& b)
        {
            ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

            Scheme<TSparseSpace, TDenseSpace>::InitializeSolutionStep(r_model_part, A, Dx, b);

            double DeltaTime = CurrentProcessInfo[DELTA_TIME];

            if (DeltaTime == 0)
                KRATOS_ERROR(std::logic_error, "detected delta_time = 0 in the Bossak Scheme ... check if the time step is created correctly for the current model part", "");

            //initializing constants
            ma0 = 1.0 / (mGammaNewmark * DeltaTime);
            ma1 = DeltaTime * mBetaNewmark / mGammaNewmark;
            ma2 = (-1 + mGammaNewmark) / mGammaNewmark;
            ma3 = DeltaTime;
            ma4 = pow(DeltaTime, 2)*(-2.0 * mBetaNewmark + 1.0) / 2.0;
            ma5 = pow(DeltaTime, 2) * mBetaNewmark;
            mam = 1.0 / (mGammaNewmark * DeltaTime);
        }
        //*************************************************************************************
        //*************************************************************************************

        virtual void InitializeNonLinIteration(ModelPart& r_model_part,
                                               TSystemMatrixType& A,
                                               TSystemVectorType& Dx,
                                               TSystemVectorType& b)
        {
            KRATOS_TRY

            ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

            if (mpTurbulenceModel != 0) // If not null
                mpTurbulenceModel->Execute();

            //if orthogonal subscales are computed
            if (CurrentProcessInfo[OSS_SWITCH] == 1.0) {
                std::cout << ">>>>>>>>>>>>>>>Using OSS<<<<<<<<<<<<<<<<<<<" << std::endl;
                for (typename ModelPart::NodesContainerType::iterator ind = r_model_part.NodesBegin(); ind != r_model_part.NodesEnd(); ind++) {

                    noalias(ind->FastGetSolutionStepValue(ADVPROJ)) = ZeroVector(3);

                    ind->FastGetSolutionStepValue(DIVPROJ) = 0.0;

                    ind->FastGetSolutionStepValue(NODAL_AREA) = 0.0;


                }//end of loop over nodes

                //loop on nodes to compute ADVPROJ   CONVPROJ NODALAREA
                array_1d<double, 3 > output;


                for (typename ModelPart::ElementsContainerType::iterator elem = r_model_part.ElementsBegin(); elem != r_model_part.ElementsEnd(); elem++)
                {
                    elem->Calculate(ADVPROJ, output, CurrentProcessInfo);
                }

                r_model_part.GetCommunicator().AssembleCurrentData(NODAL_AREA);
                r_model_part.GetCommunicator().AssembleCurrentData(DIVPROJ);
                r_model_part.GetCommunicator().AssembleCurrentData(ADVPROJ);


                for (typename ModelPart::NodesContainerType::iterator ind = r_model_part.NodesBegin(); ind != r_model_part.NodesEnd(); ind++)
                {
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

            KRATOS_CATCH("")
        }

        void FinalizeSolutionStep(ModelPart &rModelPart, TSystemMatrixType &A, TSystemVectorType &Dx, TSystemVectorType &b)
        {
            Element::EquationIdVectorType EquationId;
            LocalSystemVectorType RHS_Contribution;
            LocalSystemMatrixType LHS_Contribution;
            ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

            for (ModelPart::NodeIterator itNode = rModelPart.NodesBegin(); itNode != rModelPart.NodesEnd(); ++itNode)
            {
                itNode->FastGetSolutionStepValue(REACTION_X,0) = 0.0;
                itNode->FastGetSolutionStepValue(REACTION_Y,0) = 0.0;
                itNode->FastGetSolutionStepValue(REACTION_Z,0) = 0.0;
            }

            for (ModelPart::ElementsContainerType::ptr_iterator itElem = rModelPart.Elements().ptr_begin(); itElem != rModelPart.Elements().ptr_end(); ++itElem)
            {
                (*itElem)->InitializeNonLinearIteration(CurrentProcessInfo);
                //KRATOS_WATCH(LHS_Contribution);
                //basic operations for the element considered
                (*itElem)->CalculateLocalSystem(LHS_Contribution, RHS_Contribution, CurrentProcessInfo);

                //std::cout << rCurrentElement->Id() << " RHS = " << RHS_Contribution << std::endl;
                (*itElem)->MassMatrix(mMass[0], CurrentProcessInfo);
                (*itElem)->CalculateLocalVelocityContribution(mDamp[0], RHS_Contribution, CurrentProcessInfo);

                (*itElem)->EquationIdVector(EquationId, CurrentProcessInfo);

                //adding the dynamic contributions (statics is already included)
                AddDynamicsToLHS(LHS_Contribution, mDamp[0], mMass[0], CurrentProcessInfo);
                AddDynamicsToRHS((*itElem), RHS_Contribution, mDamp[0], mMass[0], CurrentProcessInfo);

                GeometryType& rGeom = (*itElem)->GetGeometry();
                unsigned int NumNodes = rGeom.PointsNumber();
                unsigned int Dimension = rGeom.WorkingSpaceDimension();
                unsigned int index = 0;

                for (unsigned int i = 0; i < NumNodes; i++)
                {
                    rGeom[i].FastGetSolutionStepValue(REACTION_X,0) -= RHS_Contribution[index++];
                    rGeom[i].FastGetSolutionStepValue(REACTION_Y,0) -= RHS_Contribution[index++];
                    if (Dimension == 3) rGeom[i].FastGetSolutionStepValue(REACTION_Z,0) -= RHS_Contribution[index++];
                    index++; // skip pressure dof
                }
            }

            rModelPart.GetCommunicator().AssembleCurrentData(REACTION);
        }

        //************************************************************************************************
        //************************************************************************************************

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
        double mBetaNewmark;
        double mGammaNewmark;
        double mMeshVelocity;

        double ma0;
        double ma1;
        double ma2;
        double ma3;
        double ma4;
        double ma5;
        double mam;

        std::vector< Matrix > mMass;
        std::vector< Matrix > mDamp;
        std::vector< Vector > mvel;
        std::vector< Vector > macc;
        std::vector< Vector > maccold;

        /*@} */
        /**@name Protected Operators*/
        /*@{ */

        //*********************************************************************************
        //Updating first time Derivative
        //*********************************************************************************
        void UpdateDisplacement(array_1d<double, 3 > & CurrentDisplacement,
                                const array_1d<double, 3 > & OldDisplacement,
                                const array_1d<double, 3 > & OldVelocity,
                                const array_1d<double, 3 > & OldAcceleration,
                                const array_1d<double, 3 > & CurrentAcceleration)
        {
            noalias(CurrentDisplacement) = OldDisplacement + ma3 * OldVelocity + ma4 * OldAcceleration + ma5*CurrentAcceleration;
        }



        //**************************************************************************

        void UpdateAcceleration(array_1d<double, 3 > & CurrentAcceleration,
                                const array_1d<double, 3 > & DeltaVel,
                                const array_1d<double, 3 > & OldAcceleration)
        {
            noalias(CurrentAcceleration) = ma0 * DeltaVel + ma2*OldAcceleration;
        }




        //****************************************************************************

        /**
        Kdyn = a0*M + D + a1*K
         */
        void AddDynamicsToLHS(LocalSystemMatrixType& LHS_Contribution,
                              LocalSystemMatrixType& D,
                              LocalSystemMatrixType& M,
                              ProcessInfo& CurrentProcessInfo)
        {

            //multipling time scheme factor
            LHS_Contribution *= ma1;

            // adding mass contribution to the dynamic stiffness
            if (M.size1() != 0) // if M matrix declared
            {
                noalias(LHS_Contribution) += mam*M;
            }

            //adding  damping contribution
            if (D.size1() != 0) // if M matrix declared
            {
                noalias(LHS_Contribution) += D;
            }
        }





        //****************************************************************************

        /**
        bdyn = b - M*acc - D*vel
					
         */
        void AddDynamicsToRHS(Element::Pointer rCurrentElement,
                              LocalSystemVectorType& RHS_Contribution,
                              LocalSystemMatrixType& D,
                              LocalSystemMatrixType& M,
                              ProcessInfo& CurrentProcessInfo)
        {
            //adding inertia contributionDISPLACEMENT

            if (M.size1() != 0) {
                int k = OpenMPUtils::ThisThread();
                rCurrentElement->GetSecondDerivativesVector(macc[k], 0);
                (macc[k]) *= (1.00 - mAlphaBossak);
                rCurrentElement->GetSecondDerivativesVector(maccold[k], 1);
                noalias(macc[k]) += mAlphaBossak * maccold[k];
                noalias(RHS_Contribution) -= prod(M, macc[k]);
            }
        }

        void AddDynamicsToRHS(
                              Condition::Pointer rCurrentElement,
                              LocalSystemVectorType& RHS_Contribution,
                              LocalSystemMatrixType& D,
                              LocalSystemMatrixType& M,
                              ProcessInfo& CurrentProcessInfo)
        {
            //adding inertia contributionDISPLACEMENT
            if (M.size1() != 0)
            {
                int k = OpenMPUtils::ThisThread();
                rCurrentElement->GetSecondDerivativesVector(macc[k], 0);
                (macc[k]) *= (1.00 - mAlphaBossak);
                rCurrentElement->GetSecondDerivativesVector(maccold[k], 1);
                noalias(macc[k]) += mAlphaBossak * maccold[k];

                noalias(RHS_Contribution) -= prod(M, macc[k]);
            }
        }

        /*@} */
        /**@name Protected Operations*/
        /*@{ */

        /// Normalize a vector.
        /**
          @param rThis the vector
          @return Original norm of the input vector
          */
        template< class TVectorType >
        double Normalize(TVectorType& rThis)
        {
            double Norm = 0;
            for(typename TVectorType::iterator iComponent = rThis.begin(); iComponent < rThis.end(); ++iComponent)
                Norm += (*iComponent)*(*iComponent);
            Norm = sqrt(Norm);
            for(typename TVectorType::iterator iComponent = rThis.begin(); iComponent < rThis.end(); ++iComponent)
                *iComponent /= Norm;
            return Norm;
        }


        /// Compute a rotation matrix to transform values from the cartesian base to one oriented with the node's normal
        /**
          The normal is read from solution step data NORMAL. Use NormalCalculationUtils::CalculateOnSimplex to
          obtain and store the nodal normal from the normals of the model's conditons.
          @param rRot The rotation matrix (output)
          @param rThisPoint The point used to orient the new coordinate system.
          @see NormalCalculationUtils
          */
        template<class TMatrixType>
        void RotationOperator(TMatrixType& rRot,
                              GeometryType::PointType& rThisPoint)
        {
            typedef boost::numeric::ublas::matrix_row<TMatrixType> ThisRowType;
            // Get the normal evaluated at the node
            const array_1d<double,3>& rNormal = rThisPoint.FastGetSolutionStepValue(NORMAL);

            if(mDomainSize == 3)
            {
                // Define the new coordinate system, where the first vector is aligned with the normal
                ThisRowType rN(rRot,0);
                for( size_t i = 0; i < 3; ++i)
                    rN[i] = rNormal[i];
                this->Normalize(rN);

                // To choose the remaining two vectors, we project the first component of the cartesian base to the tangent plane
                ThisRowType rT1(rRot,1);
                rT1(0) = 1.0;
                rT1(1) = 0.0;
                rT1(2) = 0.0;

                double dot = TDenseSpace::Dot(rN,rT1);

                // It is possible that the normal is aligned with (1,0,0), resulting in norm(rT1) = 0
                // If this is the case, repeat the procedure using (0,1,0)
                if ( fabs(dot) > 0.99 )
                {
                    rT1(0) = 0.0;
                    rT1(1) = 1.0;
                    rT1(2) = 0.0;

                    dot = TDenseSpace::Dot(rN,rT1);
                }

                // calculate projection and normalize
                rT1 -= dot * rN;
                this->Normalize(rT1);

                // The third base component is choosen as N x T1, which is normalized by construction
                ThisRowType rT2(rRot,2);
                rT2(0) = rN(1)*rT1(2) - rN(2)*rT1(1);
                rT2(1) = rN(2)*rT1(0) - rN(0)*rT1(2);
                rT2(2) = rN(0)*rT1(1) - rN(1)*rT1(0);
            }
            else //if(mDomainSize == 2)
            {
                /* The basis for the new coordinate system is (normal,tangent)
                   Tangent vector is chosen (-normal_y, normal_x) so that the resulting base
                   is right-handed.
                 */
                ThisRowType rN(rRot,0);
                ThisRowType rT(rRot,1);

                rN[0] = rNormal[0];
                rN[1] = rNormal[1];
                this->Normalize(rN);
                rT[0] = -rN[1];
                rT[1] = rN[0];
            }

        }


        /// Rotate the local system contributions so that they are oriented with each node's normal
        /**
          @param rLocalMatrix Local system matrix
          @param rLocalVector Local RHS vector
          @param rGeometry A reference to the element's (or condition's) geometry
          */
        void Rotate(LocalSystemMatrixType& rLocalMatrix,
                    LocalSystemVectorType& rLocalVector,
                    GeometryType& rGeometry)
        {
            const size_t BlockSize = mDomainSize + 1; // Number of rows associated to each node. Assuming vx,vy,vz,p (this is what ASGS and VMS elements do)
            const size_t LocalSize = rLocalVector.size(); // We expect this to work both with elements (4 nodes) and conditions (3 nodes)

            LocalSystemMatrixType Rotation = IdentityMatrix(LocalSize,LocalSize);
            bool NeedRotation = false;

            size_t Index = 0;

            for(size_t j = 0; j < rGeometry.PointsNumber(); ++j)
            {
                if( rGeometry[j].GetValue(IS_STRUCTURE) != 0.0 )
                {
/*					if(rGeometry[j].Id() == 1309)
						KRATOS_ERROR(std::logic_error,"found node 1309","")*/
					
                    NeedRotation = true;
                    MatrixBlockType Block(Rotation,range(Index,Index+mDomainSize),range(Index,Index+mDomainSize));
                    this->RotationOperator<MatrixBlockType>(Block,rGeometry[j]);
					
/*					if(rGeometry[j].Id() == 1309)
						KRATOS_WATCH(Block)*/
                }

                Index += BlockSize;
            }
            if(NeedRotation)
            {
//                LocalSystemMatrixType tmp = boost::numeric::ublas::prod(rLocalMatrix,boost::numeric::ublas::trans(Rotation));
//                rLocalMatrix = boost::numeric::ublas::prod(Rotation,tmp);
                this->ApplyRotation(rLocalMatrix,Rotation);

                LocalSystemVectorType aaa = boost::numeric::ublas::prod(Rotation,rLocalVector);
                noalias(rLocalVector) = aaa;
            }
        }

        /// RHS only version of Rotate
        void Rotate(LocalSystemVectorType& rLocalVector,
                    GeometryType& rGeometry)
        {
            const size_t BlockSize = mDomainSize + 1; // Number of rows associated to each node. Assuming vx,vy,vz,p (this is what ASGS and VMS elements do)
            const size_t LocalSize = rLocalVector.size(); // We expect this to work both with elements (4 nodes) and conditions (3 nodes)

            LocalSystemMatrixType Rotation = IdentityMatrix(LocalSize,LocalSize);
            bool NeedRotation = false;

            size_t Index = 0;

            for(size_t j = 0; j < rGeometry.PointsNumber(); ++j)
            {
                if( rGeometry[j].GetValue(IS_STRUCTURE) != 0.0 )
                {
                    NeedRotation = true;
                    MatrixBlockType Block(Rotation,range(Index,Index+mDomainSize),range(Index,Index+mDomainSize));
                    this->RotationOperator<MatrixBlockType>(Block,rGeometry[j]);
                }
                Index += BlockSize;
            }

            if(NeedRotation)
            {
                LocalSystemVectorType Tmp = boost::numeric::ublas::prod(Rotation,rLocalVector);
                noalias(rLocalVector) = Tmp;
            }
        }

        /// Apply slip boundary conditions to the rotated local contributions.
        /** This function takes the local system contributions rotated so each
          node's velocities are expressed using a base oriented with its normal
          and imposes that the normal velocity is equal to the mesh velocity in
          the normal direction.
          */
        void ApplySlipCondition(LocalSystemMatrixType& rLocalMatrix,
                                LocalSystemVectorType& rLocalVector,
                                GeometryType& rGeometry)
        {
            const size_t BlockSize = mDomainSize + 1; // Number of rows associated to each node. Assuming vx,vy,vz,p (this is what ASGS and VMS elements do)
            const size_t LocalSize = rLocalVector.size(); // We expect this to work both with elements (4 nodes) and conditions (3 nodes)

            for(size_t itNode = 0; itNode < rGeometry.PointsNumber(); ++itNode)
            {
                if( rGeometry[itNode].GetValue(IS_STRUCTURE) != 0.0 )
                {
                    // We fix the first dof (normal velocity) for each rotated block
                    size_t j = itNode * BlockSize;
                    //const double k = rLocalMatrix(j,j)+rLocalMatrix(j,j+1)+rLocalMatrix(j,j+2);

                    // If the mesh is moving, we must impose v_normal = vmesh_normal
                    array_1d<double,3> VMesh = rGeometry[itNode].FastGetSolutionStepValue(MESH_VELOCITY);
                    VMesh -= rGeometry[itNode].FastGetSolutionStepValue(VELOCITY);
                    array_1d<double,3> rN = rGeometry[itNode].FastGetSolutionStepValue(NORMAL);
                    this->Normalize(rN);

                    for( size_t i = 0; i < j; ++i) // Skip term (i,j)
                    {
                        rLocalMatrix(i,j) = 0.0;
                        rLocalMatrix(j,i) = 0.0;
                    }
                    for( size_t i = j+1; i < LocalSize; ++i)
                    {
                        rLocalMatrix(i,j) = 0.0;
                        rLocalMatrix(j,i) = 0.0;
                    }

                    rLocalVector(j) = inner_prod(rN,VMesh);
                    rLocalMatrix(j,j) = 1.0;
                }
            }
        }

        /// RHS only version of ApplySlipCondition
        void ApplySlipCondition(LocalSystemVectorType& rLocalVector,
                                GeometryType& rGeometry)
        {
            const size_t BlockSize = mDomainSize + 1; // Number of rows associated to each node. Assuming vx,vy,vz,p (this is what ASGS and VMS elements do)

            for(size_t itNode = 0; itNode < rGeometry.PointsNumber(); ++itNode)
            {
                if( rGeometry[itNode].GetValue(IS_STRUCTURE) != 0.0 )
                {
                    // We fix the first dof (normal velocity) for each rotated block
                    size_t j = itNode * BlockSize;

                    // If the mesh is moving, we must impose v_normal = vmesh_normal
                    array_1d<double,3> VMesh = rGeometry[itNode].FastGetSolutionStepValue(MESH_VELOCITY);
                    VMesh -= rGeometry[itNode].FastGetSolutionStepValue(VELOCITY);
                    array_1d<double,3> rN = rGeometry[itNode].FastGetSolutionStepValue(NORMAL);
					this->Normalize(rN);
					
                    rLocalVector[j] = inner_prod(rN,VMesh);
                }
            }
        }

        /// Transform nodal velocities to the rotated coordinates (aligned with each node's normal)
        void RotateVelocities(ModelPart& rModelPart)
        {
            LocalSystemMatrixType Rotation = IdentityMatrix(mDomainSize,mDomainSize);
            LocalSystemVectorType Vel(mDomainSize);
            LocalSystemVectorType Tmp(mDomainSize);

            for (ModelPart::NodeIterator itNode = rModelPart.NodesBegin(); itNode != rModelPart.NodesEnd(); ++itNode)
            {
                if(itNode->GetValue(IS_STRUCTURE) != 0.0)
                {
                    this->RotationOperator<LocalSystemMatrixType>(Rotation,*itNode);
                    array_1d<double,3>& rVelocity = itNode->FastGetSolutionStepValue(VELOCITY);
                    for(size_t i = 0; i < mDomainSize; i++) Vel[i] = rVelocity[i];
                    noalias(Tmp) = boost::numeric::ublas::prod(Rotation,Vel);
                    for(size_t i = 0; i < mDomainSize; i++) rVelocity[i] = Tmp[i];
                }
            }
        }

        /// Transform nodal velocities from the rotated system to the original one
        void RecoverVelocities(ModelPart& rModelPart)
        {
            LocalSystemMatrixType Rotation = IdentityMatrix(mDomainSize,mDomainSize);
            LocalSystemVectorType Vel(mDomainSize);
            LocalSystemVectorType Tmp(mDomainSize);

            for (ModelPart::NodeIterator itNode = rModelPart.NodesBegin(); itNode != rModelPart.NodesEnd(); ++itNode)
            {
                if(itNode->GetValue(IS_STRUCTURE) != 0.0)
                {
                    this->RotationOperator<LocalSystemMatrixType>(Rotation,*itNode);
                    array_1d<double,3>& rVelocity = itNode->FastGetSolutionStepValue(VELOCITY);
                    for(size_t i = 0; i < mDomainSize; i++) Vel[i] = rVelocity[i];
                    noalias(Tmp) = boost::numeric::ublas::prod(boost::numeric::ublas::trans(Rotation),Vel);
                    for(size_t i = 0; i < mDomainSize; i++) rVelocity[i] = Tmp[i];
                }
            }
        }

        /// Transform a local contribution from cartesian coordinates to rotated ones
        void ApplyRotation(LocalSystemMatrixType& rMatrix,
                           const LocalSystemMatrixType& rRotation)
        {
            // compute B = R*A*transpose(R)
            const size_t BlockSize = mDomainSize + 1;
            const size_t LocalSize = rMatrix.size1();
            const size_t NumBlocks = LocalSize / BlockSize;
            LocalSystemMatrixType Tmp = ZeroMatrix(LocalSize,LocalSize);

            for (size_t iBlock = 0; iBlock < NumBlocks; iBlock++)
            {
                for (size_t jBlock = 0; jBlock < NumBlocks; jBlock++)
                {
                    for (size_t i = iBlock*BlockSize; i < (iBlock+1)*BlockSize; i++)
                    {
                        for(size_t j = jBlock*BlockSize; j < (jBlock+1)*BlockSize; j++)
                        {
                            double& tij = Tmp(i,j);
                            for(size_t k = iBlock*BlockSize; k < (iBlock+1)*BlockSize; k++)
                            {
                                for(size_t l = jBlock*BlockSize; l < (jBlock+1)*BlockSize; l++)
                                {
                                    tij += rRotation(i,k)*rMatrix(k,l)*rRotation(j,l);
                                }
                            }
                        }
                    }
                }
            }
            noalias(rMatrix) = Tmp;
        }

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

        const size_t mDomainSize;

        Process::Pointer mpTurbulenceModel;

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

