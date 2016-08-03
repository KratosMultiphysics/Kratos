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

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/deprecated_variables.h"
#include "solving_strategies/schemes/scheme.h"
#include "includes/variables.h"
#include "includes/cfd_variables.h"
#include "containers/array_1d.h"
#include "utilities/openmp_utils.h"
#include "utilities/coordinate_transformation_utilities.h"
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
          mRotationTool(DomainSize,DomainSize+1,IS_STRUCTURE,0.0), // Second argument is number of matrix rows per node: monolithic elements have velocity and pressure dofs.
          mrPeriodicIdVar(Kratos::Variable<int>::StaticObject())
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
        }


        /** Constructor without a turbulence model with periodic conditions
         */
        ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent(
            double NewAlphaBossak,
            double MoveMeshStrategy,
            unsigned int DomainSize,
            const Variable<int>& rPeriodicIdVar)
        :
          Scheme<TSparseSpace, TDenseSpace>(),
          mRotationTool(DomainSize,DomainSize+1,IS_STRUCTURE,0.0), // Second argument is number of matrix rows per node: monolithic elements have velocity and pressure dofs.
          mrPeriodicIdVar(rPeriodicIdVar)
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
        }


        /** Constructor without a turbulence model
         */
        ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent(
            double NewAlphaBossak,
            double MoveMeshStrategy,
            unsigned int DomainSize,
            Variable<double>& rSlipVar)
        :
          Scheme<TSparseSpace, TDenseSpace>(),
          mRotationTool(DomainSize,DomainSize+1,rSlipVar,0.0), // Second argument is number of matrix rows per node: monolithic elements have velocity and pressure dofs.
          mrPeriodicIdVar(Kratos::Variable<int>::StaticObject())
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
          mRotationTool(DomainSize,DomainSize+1,IS_STRUCTURE,0.0), // Second argument is number of matrix rows per node: monolithic elements have velocity and pressure dofs
          mrPeriodicIdVar(Kratos::Variable<int>::StaticObject()),
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

            mRotationTool.RotateVelocities(r_model_part);

            BasicUpdateOperations(r_model_part, rDofSet, A, Dv, b);

            mRotationTool.RecoverVelocities(r_model_part);

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
            if (rModelPart.GetCommunicator().MyPID() == 0)
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

            if (rModelPart.GetCommunicator().MyPID() == 0)
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
            (rCurrentElement)->CalculateMassMatrix(mMass[k], CurrentProcessInfo);
            (rCurrentElement)->CalculateLocalVelocityContribution(mDamp[k], RHS_Contribution, CurrentProcessInfo);

            (rCurrentElement)->EquationIdVector(EquationId, CurrentProcessInfo);

            //adding the dynamic contributions (statics is already included)

            AddDynamicsToLHS(LHS_Contribution, mDamp[k], mMass[k], CurrentProcessInfo);
            AddDynamicsToRHS(rCurrentElement, RHS_Contribution, mDamp[k], mMass[k], CurrentProcessInfo);

            // If there is a slip condition, apply it on a rotated system of coordinates
            mRotationTool.Rotate(LHS_Contribution,RHS_Contribution,rCurrentElement->GetGeometry());
            mRotationTool.ApplySlipCondition(LHS_Contribution,RHS_Contribution,rCurrentElement->GetGeometry());

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
            (rCurrentElement)->CalculateMassMatrix(mMass[k], CurrentProcessInfo);

            (rCurrentElement)->CalculateLocalVelocityContribution(mDamp[k], RHS_Contribution, CurrentProcessInfo);

            (rCurrentElement)->EquationIdVector(EquationId, CurrentProcessInfo);

            //adding the dynamic contributions (static is already included)

            AddDynamicsToRHS(rCurrentElement, RHS_Contribution, mDamp[k], mMass[k], CurrentProcessInfo);

            // If there is a slip condition, apply it on a rotated system of coordinates
            mRotationTool.Rotate(RHS_Contribution,rCurrentElement->GetGeometry());
            mRotationTool.ApplySlipCondition(RHS_Contribution,rCurrentElement->GetGeometry());
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
            (rCurrentCondition)->CalculateMassMatrix(mMass[k], CurrentProcessInfo);
            //(rCurrentCondition)->CalculateDampingMatrix(VelocityBossakAuxiliaries::mDamp,CurrentProcessInfo);
            (rCurrentCondition)->CalculateLocalVelocityContribution(mDamp[k], RHS_Contribution, CurrentProcessInfo);
            (rCurrentCondition)->EquationIdVector(EquationId, CurrentProcessInfo);


            AddDynamicsToLHS(LHS_Contribution, mDamp[k], mMass[k], CurrentProcessInfo);

            AddDynamicsToRHS(rCurrentCondition, RHS_Contribution, mDamp[k], mMass[k], CurrentProcessInfo);

            // Rotate contributions (to match coordinates for slip conditions)
            mRotationTool.Rotate(LHS_Contribution,RHS_Contribution,rCurrentCondition->GetGeometry());
            mRotationTool.ApplySlipCondition(LHS_Contribution,RHS_Contribution,rCurrentCondition->GetGeometry());

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
            (rCurrentCondition)->CalculateMassMatrix(mMass[k],rCurrentProcessInfo);
            //(rCurrentCondition)->CalculateDampingMatrix(VelocityBossakAuxiliaries::mDamp,CurrentProcessInfo);
            (rCurrentCondition)->CalculateLocalVelocityContribution(mDamp[k], RHS_Contribution,rCurrentProcessInfo);
            (rCurrentCondition)->EquationIdVector(EquationId,rCurrentProcessInfo);

            //adding the dynamic contributions (static is already included)
            AddDynamicsToRHS(rCurrentCondition, RHS_Contribution, mDamp[k], mMass[k],rCurrentProcessInfo);

            // Rotate contributions (to match coordinates for slip conditions)
            mRotationTool.Rotate(RHS_Contribution,rCurrentCondition->GetGeometry());
            mRotationTool.ApplySlipCondition(RHS_Contribution,rCurrentCondition->GetGeometry());



            KRATOS_CATCH("");
        }
        //*************************************************************************************
        //*************************************************************************************

        virtual void InitializeSolutionStep(ModelPart& r_model_part,
                                            TSystemMatrixType& A,
                                            TSystemVectorType& Dx,
                                            TSystemVectorType& b)
        {
            ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

            Scheme<TSparseSpace, TDenseSpace>::InitializeSolutionStep(r_model_part, A, Dx, b);

            double DeltaTime = CurrentProcessInfo[DELTA_TIME];

            if (DeltaTime == 0)
                KRATOS_THROW_ERROR(std::logic_error, "detected delta_time = 0 in the Bossak Scheme ... check if the time step is created correctly for the current model part", "");

            //initializing constants
            ma0 = 1.0 / (mGammaNewmark * DeltaTime);
            ma1 = DeltaTime * mBetaNewmark / mGammaNewmark;
            ma2 = (-1 + mGammaNewmark) / mGammaNewmark;
            ma3 = DeltaTime;
            ma4 = pow(DeltaTime, 2)*(-2.0 * mBetaNewmark + 1.0) / 2.0;
            ma5 = pow(DeltaTime, 2) * mBetaNewmark;
            mam = (1.0 - mAlphaBossak) / (mGammaNewmark * DeltaTime);
        }
        //*************************************************************************************
        //*************************************************************************************

        virtual void InitializeNonLinIteration(ModelPart& r_model_part,
                                               TSystemMatrixType& A,
                                               TSystemVectorType& Dx,
                                               TSystemVectorType& b)
        {
            KRATOS_TRY

            if (mpTurbulenceModel != 0) // If not null
                mpTurbulenceModel->Execute();

            KRATOS_CATCH("")
        }

        virtual void FinalizeNonLinIteration(ModelPart &rModelPart, TSystemMatrixType &A, TSystemVectorType &Dx, TSystemVectorType &b)
        {
            ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

            //if orthogonal subscales are computed
            if (CurrentProcessInfo[OSS_SWITCH] == 1.0) {
                if (rModelPart.GetCommunicator().MyPID() == 0)
                    std::cout << "Computing OSS projections" << std::endl;
                for (typename ModelPart::NodesContainerType::iterator ind = rModelPart.NodesBegin(); ind != rModelPart.NodesEnd(); ind++) {

                    noalias(ind->FastGetSolutionStepValue(ADVPROJ)) = ZeroVector(3);

                    ind->FastGetSolutionStepValue(DIVPROJ) = 0.0;

                    ind->FastGetSolutionStepValue(NODAL_AREA) = 0.0;


                }//end of loop over nodes

                //loop on nodes to compute ADVPROJ   CONVPROJ NODALAREA
                array_1d<double, 3 > output;


                for (typename ModelPart::ElementsContainerType::iterator elem = rModelPart.ElementsBegin(); elem != rModelPart.ElementsEnd(); elem++)
                {
                    elem->Calculate(ADVPROJ, output, CurrentProcessInfo);
                }

                rModelPart.GetCommunicator().AssembleCurrentData(NODAL_AREA);
                rModelPart.GetCommunicator().AssembleCurrentData(DIVPROJ);
                rModelPart.GetCommunicator().AssembleCurrentData(ADVPROJ);

                // Correction for periodic conditions
                this->PeriodicConditionProjectionCorrection(rModelPart);

                for (typename ModelPart::NodesContainerType::iterator ind = rModelPart.NodesBegin(); ind != rModelPart.NodesEnd(); ind++)
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
                (*itElem)->CalculateMassMatrix(mMass[0], CurrentProcessInfo);
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

            // Base scheme calls FinalizeSolutionStep method of elements and conditions
            Scheme<TSparseSpace, TDenseSpace>::FinalizeSolutionStep(rModelPart, A, Dx, b);
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


        /** On periodic boundaries, the nodal area and the values to project need to take into account contributions from elements on
         * both sides of the boundary. This is done using the conditions and the non-historical nodal data containers as follows:\n
         * 1- The partition that owns the PeriodicCondition adds the values on both nodes to their non-historical containers.\n
         * 2- The non-historical containers are added across processes, transmiting the right value from the condition owner to all partitions.\n
         * 3- The value on all periodic nodes is replaced by the one received in step 2.
         */
        void PeriodicConditionProjectionCorrection(ModelPart& rModelPart)
        {
            if (mrPeriodicIdVar.Key() != 0)
            {
                for (typename ModelPart::ConditionIterator itCond = rModelPart.ConditionsBegin(); itCond != rModelPart.ConditionsEnd(); itCond++ )
                {
                    ModelPart::ConditionType::GeometryType& rGeom = itCond->GetGeometry();
                    if (rGeom.PointsNumber() == 2)
                    {
                        Node<3>& rNode0 = rGeom[0];
                        int Node0Pair = rNode0.FastGetSolutionStepValue(mrPeriodicIdVar);

                        Node<3>& rNode1 = rGeom[1];
                        int Node1Pair = rNode1.FastGetSolutionStepValue(mrPeriodicIdVar);

                        // If the nodes are marked as a periodic pair (this is to avoid acting on two-noded conditions that are not PeriodicCondition)
                        if ( ( static_cast<int>(rNode0.Id()) == Node1Pair ) && (static_cast<int>(rNode1.Id()) == Node0Pair ) )
                        {
                            double NodalArea = rNode0.FastGetSolutionStepValue(NODAL_AREA) + rNode1.FastGetSolutionStepValue(NODAL_AREA);
                            array_1d<double,3> AdvProj = rNode0.FastGetSolutionStepValue(ADVPROJ) + rNode1.FastGetSolutionStepValue(ADVPROJ);
                            double DivProj = rNode0.FastGetSolutionStepValue(DIVPROJ) + rNode1.FastGetSolutionStepValue(DIVPROJ);

                            rNode0.GetValue(NODAL_AREA) = NodalArea;
                            rNode0.GetValue(ADVPROJ) = AdvProj;
                            rNode0.GetValue(DIVPROJ) = DivProj;

                            rNode1.GetValue(NODAL_AREA) = NodalArea;
                            rNode1.GetValue(ADVPROJ) = AdvProj;
                            rNode1.GetValue(DIVPROJ) = DivProj;
                        }
                    }
                    else if (rGeom.PointsNumber() == 4)
                    {
                        double NodalArea = 0.0;
                        array_1d<double,3> AdvProj(3,0.0);
                        double DivProj = 0.0;
                        for ( unsigned int i = 0; i < 4; i++ )
                        {
                            NodalArea += rGeom[i].FastGetSolutionStepValue(NODAL_AREA);
                            AdvProj += rGeom[i].FastGetSolutionStepValue(ADVPROJ);
                            DivProj += rGeom[i].FastGetSolutionStepValue(DIVPROJ);
                        }

                        for ( unsigned int i = 0; i < 4; i++ )
                        {
                            rGeom[i].GetValue(NODAL_AREA) = NodalArea;
                            rGeom[i].GetValue(ADVPROJ) = AdvProj;
                            rGeom[i].GetValue(DIVPROJ) = DivProj;
                        }

                    }
                }

                rModelPart.GetCommunicator().AssembleNonHistoricalData(NODAL_AREA);
                rModelPart.GetCommunicator().AssembleNonHistoricalData(ADVPROJ);
                rModelPart.GetCommunicator().AssembleNonHistoricalData(DIVPROJ);

                for (typename ModelPart::NodeIterator itNode = rModelPart.NodesBegin(); itNode != rModelPart.NodesEnd(); itNode++)
                {
                    if (itNode->GetValue(NODAL_AREA) != 0.0)
                    {
                        itNode->FastGetSolutionStepValue(NODAL_AREA) = itNode->GetValue(NODAL_AREA);
                        itNode->FastGetSolutionStepValue(ADVPROJ) = itNode->GetValue(ADVPROJ);
                        itNode->FastGetSolutionStepValue(DIVPROJ) = itNode->GetValue(DIVPROJ);

                        // reset for next iteration
                        itNode->GetValue(NODAL_AREA) = 0.0;
                        itNode->GetValue(ADVPROJ) = array_1d<double,3>(3,0.0);
                        itNode->GetValue(DIVPROJ) = 0.0;
                    }
                }
            }
        }






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
        Kdyn = am*M + D + a1*K
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

