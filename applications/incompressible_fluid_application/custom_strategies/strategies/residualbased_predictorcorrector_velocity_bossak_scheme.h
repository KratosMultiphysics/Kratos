/*
==============================================================================
KratosStructuralApplication
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
/* *********************************************************
 *
 *   Last Modified by:    $Author: Kazem $
 *   Date:                $Date: 2008-07-25 14:48:17 $
 *   Revision:            $Revision: 1.1 $
 *
 * ***********************************************************/


#if !defined(KRATOS_RESIDUALBASED_PREDICTOR_CORRECTOR_VELOCITY_BOSSAK_SCHEME )
#define  KRATOS_RESIDUALBASED_PREDICTOR_CORRECTOR_VELOCITY_BOSSAK_SCHEME


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "includes/variables.h"
#include "containers/array_1d.h"
#include "utilities/openmp_utils.h"
#include "includes/cfd_variables.h"

namespace Kratos
{

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

/** Short class definition.

  This class provides the implementation of the basic tasks that are needed by the solution strategy.
  It is intended to be the place for tailoring the solution strategies to problem specific tasks.

        Detail class definition.

          \URL[Example of use html]{ extended_documentation/no_ex_of_use.html}

                \URL[Example of use pdf]{ extended_documentation/no_ex_of_use.pdf}

                  \URL[Example of use doc]{ extended_documentation/no_ex_of_use.doc}

                        \URL[Example of use ps]{ extended_documentation/no_ex_of_use.ps}


                                \URL[Extended documentation html]{ extended_documentation/no_ext_doc.html}

                                  \URL[Extended documentation pdf]{ extended_documentation/no_ext_doc.pdf}

                                        \URL[Extended documentation doc]{ extended_documentation/no_ext_doc.doc}

                                          \URL[Extended documentation ps]{ extended_documentation/no_ext_doc.ps}


 */
template<class TSparseSpace,
         class TDenseSpace //= DenseSpace<double>
         >
class ResidualBasedPredictorCorrectorVelocityBossakScheme : public Scheme<TSparseSpace, TDenseSpace>
{
public:
    /**@name Type Definitions */
    /*@{ */

    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedPredictorCorrectorVelocityBossakScheme);

    typedef Scheme<TSparseSpace, TDenseSpace> BaseType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename Element::DofsVectorType DofsVectorType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;


    /*@} */
    /**@name Life Cycle
     */
    /*@{ */

    /** Constructor.
     */
    ResidualBasedPredictorCorrectorVelocityBossakScheme(double NewAlphaBossak, double MoveMeshStrategy)
        : Scheme<TSparseSpace, TDenseSpace>()
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



// mAlphaBossak= 0.0;
// mGammaNewmark= 1.0;
        //mGammaNewmark = 1.0;
        //mBetaNewmark = 0.5;
        //mAlphaBossak = 0.0;
        //sizing work matrices
        //mMass.resize(10,10);
        //mDamp.resize(10,10);

        //mvel.resize(10,false);
        //macc.resize(10,false);
        //maccold.resize(10,false);


        std::cout << "using the velocity Bossak Time Integration Scheme" << std::endl;
    }

    /** Destructor.
     */
    virtual ~ResidualBasedPredictorCorrectorVelocityBossakScheme()
    {
    }


    /*@} */
    /**@name Operators
     */
    /*@{ */

    /**
            Performing the update of the solution.
     */

    //***************************************************************************

    void Update(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dv,
        TSystemVectorType& b
    ) override
    {
        KRATOS_TRY

        BasicUpdateOperations(r_model_part, rDofSet, A, Dv, b);

        AdditionalUpdateOperations(r_model_part, rDofSet, A, Dv, b);

        KRATOS_CATCH("")
    }

    void BasicUpdateOperations(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b)
    {
        mpDofUpdater->UpdateDofs(rDofSet,Dx);
    }


    void AdditionalUpdateOperations(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dv,
        TSystemVectorType& b
    )
    {
        KRATOS_TRY

        int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector NodePartition;
        OpenMPUtils::DivideInPartitions(rModelPart.Nodes().size(),NumThreads,NodePartition);

        //updating time derivatives (nodally for efficiency)
        #pragma omp parallel
        {
            array_1d<double, 3 > DeltaVel;

            int k = OpenMPUtils::ThisThread();

            ModelPart::NodeIterator NodesBegin = rModelPart.NodesBegin() + NodePartition[k];
            ModelPart::NodeIterator NodesEnd = rModelPart.NodesBegin() + NodePartition[k+1];

            for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; itNode++)
            {
                noalias(DeltaVel) = (itNode)->FastGetSolutionStepValue(VELOCITY) - (itNode)->FastGetSolutionStepValue(VELOCITY, 1);

                array_1d<double, 3 > & CurrentAcceleration = (itNode)->FastGetSolutionStepValue(ACCELERATION, 0);
                array_1d<double, 3 > & OldAcceleration = (itNode)->FastGetSolutionStepValue(ACCELERATION, 1);

                UpdateAcceleration(CurrentAcceleration, DeltaVel, OldAcceleration);

                if (mMeshVelocity == 2)//Lagrangian
                {
                    array_1d<double, 3 > & CurrentDisplacement = (itNode)->FastGetSolutionStepValue(DISPLACEMENT, 0);
                    array_1d<double, 3 > & OldDisplacement = (itNode)->FastGetSolutionStepValue(DISPLACEMENT, 1);

                    array_1d<double, 3 > & OldVelocity = (itNode)->FastGetSolutionStepValue(VELOCITY, 1);

                    noalias(itNode->FastGetSolutionStepValue(MESH_VELOCITY) ) = itNode->FastGetSolutionStepValue(VELOCITY);
                    UpdateDisplacement(CurrentDisplacement, OldDisplacement, OldVelocity, OldAcceleration, CurrentAcceleration);
                }
            }
        }

        KRATOS_CATCH("")

    }


    //***************************************************************************
    //predicts the solution at the current step as
    // v = vold

    void Predict(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dv,
        TSystemVectorType& b
    ) override
    {
//             std::cout << "prediction" << std::endl;

        int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector NodePartition;
        OpenMPUtils::DivideInPartitions(rModelPart.Nodes().size(),NumThreads,NodePartition);

        #pragma omp parallel
        {
            //array_1d<double, 3 > DeltaDisp;

            int k = OpenMPUtils::ThisThread();

            ModelPart::NodeIterator NodesBegin = rModelPart.NodesBegin() + NodePartition[k];
            ModelPart::NodeIterator NodesEnd = rModelPart.NodesBegin() + NodePartition[k+1];

            for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; itNode++)
            {
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

                    noalias(itNode->FastGetSolutionStepValue(MESH_VELOCITY) ) = itNode->FastGetSolutionStepValue(VELOCITY);
                    UpdateDisplacement(CurrentDisplacement, OldDisplacement, OldVelocity, OldAcceleration, CurrentAcceleration);
                }
            }
        }

//             std::cout << "end of prediction" << std::endl;

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
    void CalculateSystemContributions(
        Element::Pointer rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo
    ) override
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



        KRATOS_CATCH("")

    }

    void Calculate_RHS_Contribution(
        Element::Pointer rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo) override
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

    }

    /** functions totally analogous to the precedent but applied to
    the "condition" objects
     */
    virtual void Condition_CalculateSystemContributions(
        Condition::Pointer rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo) override
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

        KRATOS_CATCH("")
    }

    void Condition_Calculate_RHS_Contribution(
        Condition::Pointer rCurrentCondition,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        int k = OpenMPUtils::ThisThread();

        //KRATOS_WATCH("CONDITION LOCALVELOCITYCONTRIBUTION IS NOT DEFINED");
        //Initializing the non linear iteration for the current condition
        (rCurrentCondition) -> InitializeNonLinearIteration(CurrentProcessInfo);

        //basic operations for the element considered
        (rCurrentCondition)->CalculateRightHandSide(RHS_Contribution, CurrentProcessInfo);
        (rCurrentCondition)->CalculateMassMatrix(mMass[k], CurrentProcessInfo);
        //(rCurrentCondition)->CalculateDampingMatrix(VelocityBossakAuxiliaries::mDamp,CurrentProcessInfo);
        (rCurrentCondition)->CalculateLocalVelocityContribution(mDamp[k], RHS_Contribution, CurrentProcessInfo);
        (rCurrentCondition)->EquationIdVector(EquationId, CurrentProcessInfo);

        //adding the dynamic contributions (static is already included)

        AddDynamicsToRHS(rCurrentCondition, RHS_Contribution, mDamp[k], mMass[k], CurrentProcessInfo);

        KRATOS_CATCH("")
    }
    //*************************************************************************************
    //*************************************************************************************

    void InitializeSolutionStep(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    ) override
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

    void InitializeNonLinIteration(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        //if orthogonal subscales are computed
        if (CurrentProcessInfo[OSS_SWITCH] == 1.0)
        {
            std::cout << ">>>>>>>>>>>>>>>Using OSS<<<<<<<<<<<<<<<<<<<" << std::endl;
            for (typename ModelPart::NodesContainerType::iterator ind = r_model_part.NodesBegin(); ind != r_model_part.NodesEnd(); ind++)
            {

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
    //************************************************************************************************
    //************************************************************************************************

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

    std::vector< Matrix >mMass;
    std::vector< Matrix >mDamp;
    std::vector< Vector >mvel;
    std::vector< Vector >macc;
    std::vector< Vector >maccold;





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


    void AddDynamicsToLHS(
        LocalSystemMatrixType& LHS_Contribution,
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
    void AddDynamicsToRHS(
        Element::Pointer rCurrentElement,
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

        //adding damping contribution
        //damping contribution

//             if (D.size1() != 0) {
//                 rCurrentElement->GetFirstDerivativesVector(VelocityBossakAuxiliaries::mvel, 0);
//                 noalias(RHS_Contribution) -= prod(D, VelocityBossakAuxiliaries::mvel);
//             }


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

        //adding damping contribution
        //damping contribution

//             if (D.size1() != 0) {
//                 rCurrentElement->GetFirstDerivativesVector(VelocityBossakAuxiliaries::mvel, 0);
//                 noalias(RHS_Contribution) -= prod(D, VelocityBossakAuxiliaries::mvel);
//             }


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

    typename TSparseSpace::DofUpdaterPointerType mpDofUpdater = TSparseSpace::CreateDofUpdater();

    /*		Matrix mMass;
                    Matrix mDamp;

                    Vector mvel;
                    Vector macc;
                    Vector maccold;

                    DofsVectorType mElementalDofList;
     */

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

