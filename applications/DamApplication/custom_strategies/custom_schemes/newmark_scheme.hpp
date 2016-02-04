//
//   Project Name:
//   Last modified by:    $Author:
//   Date:                $Date:
//   Revision:            $Revision:
//

#if !defined(KRATOS_NEWMARK_SCHEME )
#define  KRATOS_NEWMARK_SCHEME

/* System includes */

/* External includes */
#include "boost/smart_ptr.hpp"

/* Project includes */
#include "containers/array_1d.h"
#include "includes/variables.h"
#include "includes/element.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"

#include "dam_application.h"

namespace Kratos
{

template<class TSparseSpace, class TDenseSpace>

class NewmarkScheme : public Scheme<TSparseSpace,TDenseSpace>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( NewmarkScheme );

    typedef Scheme<TSparseSpace,TDenseSpace>                      BaseType;
    typedef typename BaseType::DofsArrayType                 DofsArrayType;
    typedef typename Element::DofsVectorType                DofsVectorType;
    typedef typename BaseType::TSystemMatrixType         TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType         TSystemVectorType;
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
    typedef ModelPart::ElementsContainerType             ElementsArrayType;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    //Constructor
    NewmarkScheme(bool is_dynamic) : Scheme<TSparseSpace,TDenseSpace>()
    {
        mIsDynamic = is_dynamic;

        //Allocate auxiliary memory
        int NumThreads = OpenMPUtils::GetNumThreads();
        mMassMatrix.resize(NumThreads);
        mAccelerationVector.resize(NumThreads);
        mDampingMatrix.resize(NumThreads);
        mVelocityVector.resize(NumThreads);
    }

    //------------------------------------------------------------------------------------

    //Destructor
    virtual ~NewmarkScheme() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    int Check(ModelPart& r_model_part)
    {
        KRATOS_TRY

        int err = Scheme<TSparseSpace, TDenseSpace>::Check(r_model_part);
        if(err!=0) return err;

        //check for variables keys (verify that the variables are correctly initialized)
        if(DISPLACEMENT.Key() == 0)
            KRATOS_THROW_ERROR( std::invalid_argument,"DISPLACEMENT has Key zero! (check if the application is correctly registered", "" )
        if(VELOCITY.Key() == 0)
            KRATOS_THROW_ERROR( std::invalid_argument,"VELOCITY has Key zero! (check if the application is correctly registered", "" )
        if(ACCELERATION.Key() == 0)
            KRATOS_THROW_ERROR( std::invalid_argument,"ACCELERATION has Key zero! (check if the application is correctly registered", "" )

        //check that variables are correctly allocated
        for(ModelPart::NodesContainerType::iterator it=r_model_part.NodesBegin(); it!=r_model_part.NodesEnd(); it++)
        {
            if(it->SolutionStepsDataHas(DISPLACEMENT) == false)
                KRATOS_THROW_ERROR( std::logic_error, "DISPLACEMENT variable is not allocated for node ", it->Id() )
            if(it->SolutionStepsDataHas(VELOCITY) == false)
                KRATOS_THROW_ERROR( std::logic_error, "DISPLACEMENT variable is not allocated for node ", it->Id() )
            if(it->SolutionStepsDataHas(ACCELERATION) == false)
                KRATOS_THROW_ERROR( std::logic_error, "DISPLACEMENT variable is not allocated for node ", it->Id() )
        }

        //check that dofs exist
        for(ModelPart::NodesContainerType::iterator it=r_model_part.NodesBegin(); it!=r_model_part.NodesEnd(); it++)
        {
            if(it->HasDofFor(DISPLACEMENT_X) == false)
                KRATOS_THROW_ERROR( std::invalid_argument,"missing DISPLACEMENT_X dof on node ",it->Id() )
            if(it->HasDofFor(DISPLACEMENT_Y) == false)
                KRATOS_THROW_ERROR( std::invalid_argument,"missing DISPLACEMENT_Y dof on node ",it->Id() )
            if(it->HasDofFor(DISPLACEMENT_Z) == false)
                KRATOS_THROW_ERROR( std::invalid_argument,"missing DISPLACEMENT_Z dof on node ",it->Id() )
        }

        //check for minimum value of the buffer index.
        if (r_model_part.GetBufferSize() < 2)
            KRATOS_THROW_ERROR( std::logic_error, "insufficient buffer size. Buffer size should be greater than 2. Current size is", r_model_part.GetBufferSize() )

        if (r_model_part.GetProcessInfo()[DELTA_TIME] < 1e-25)
            KRATOS_THROW_ERROR( std::logic_error, "Detected DELTA_TIME < 1e-25 in the Solution Scheme. DELTA_TIME should be greater than 0", "" )

        return 0;

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Initialize(ModelPart& r_model_part)
    {
        KRATOS_TRY

        mBeta = r_model_part.GetProcessInfo()[NEWMARK_BETA];
        mGamma = r_model_part.GetProcessInfo()[NEWMARK_GAMMA];

        this->mSchemeIsInitialized = true;

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void InitializeNonLinIteration(ModelPart& r_model_part,TSystemMatrixType& A,TSystemVectorType& Dx,TSystemVectorType& b)
    {
        KRATOS_TRY

        //initialize non linear iteration for all of the elements
        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
        ElementsArrayType& pElements = r_model_part.Elements();

        for (ElementsArrayType::iterator it = pElements.begin(); it != pElements.end(); ++it)
        {
            (it) -> InitializeNonLinearIteration(CurrentProcessInfo);
        }

        /*
        ConditionsArrayType& pConditions = r_model_part.Conditions();
        for (ConditionsArrayType::iterator it = pConditions.begin(); it != pConditions.end(); ++it)
        {
            (it) -> InitializeNonLinearIteration(CurrentProcessInfo);
        }
        */

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Predict(ModelPart& r_model_part,DofsArrayType& rDofSet,TSystemMatrixType& A,TSystemVectorType& Dx,TSystemVectorType& b)
    {
        array_1d<double,3> DeltaDisplacement;
        double DeltaTime = r_model_part.GetProcessInfo()[DELTA_TIME];

        for (ModelPart::NodeIterator i = r_model_part.NodesBegin();i != r_model_part.NodesEnd(); ++i)
        {
            //Update imposed conditions (linearly incremented) 
            array_1d<double,3>& CurrentDisplacement = (i)->FastGetSolutionStepValue(DISPLACEMENT);
            noalias(CurrentDisplacement) += (i)->FastGetSolutionStepValue(IMPOSED_DISPLACEMENT);

            array_1d<double,3>& PointLoad = (i)->FastGetSolutionStepValue(POINT_LOAD);
            noalias(PointLoad) += (i)->FastGetSolutionStepValue(IMPOSED_POINT_LOAD);

            array_1d<double,3>& LineLoad = (i)->FastGetSolutionStepValue(LINE_LOAD);
            noalias(LineLoad) += (i)->FastGetSolutionStepValue(IMPOSED_LINE_LOAD);

            array_1d<double,3>& SurfaceLoad = (i)->FastGetSolutionStepValue(SURFACE_LOAD);
            noalias(SurfaceLoad) += (i)->FastGetSolutionStepValue(IMPOSED_SURFACE_LOAD);

            double& NormalContactStress = (i)->FastGetSolutionStepValue(NORMAL_CONTACT_STRESS);
            NormalContactStress += (i)->FastGetSolutionStepValue(IMPOSED_NORMAL_STRESS);

            double& TangentialContactStress = (i)->FastGetSolutionStepValue(TANGENTIAL_CONTACT_STRESS);
            TangentialContactStress += (i)->FastGetSolutionStepValue(IMPOSED_TANGENTIAL_STRESS);
            
            double& Temperature = (i)->FastGetSolutionStepValue(TEMPERATURE);
            Temperature += (i)->FastGetSolutionStepValue(IMPOSED_TEMPERATURE);


            //Predict Acceleration and Velocity
            noalias(DeltaDisplacement)               = CurrentDisplacement - (i)->FastGetSolutionStepValue(DISPLACEMENT, 1);
            array_1d<double,3>& CurrentAcceleration  = (i)->FastGetSolutionStepValue(ACCELERATION);
            array_1d<double,3>& PreviousAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION, 1);
            array_1d<double,3>& CurrentVelocity      = (i)->FastGetSolutionStepValue(VELOCITY);
            array_1d<double,3>& PreviousVelocity     = (i)->FastGetSolutionStepValue(VELOCITY, 1);

            UpdateAcceleration(CurrentAcceleration, DeltaDisplacement, PreviousVelocity, PreviousAcceleration, DeltaTime);
            UpdateVelocity(CurrentVelocity, PreviousVelocity, PreviousAcceleration, CurrentAcceleration, DeltaTime);
        }
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateSystemContributions(Element::Pointer rCurrentElement,LocalSystemMatrixType& LHS_Contribution,LocalSystemVectorType& RHS_Contribution,
                                        Element::EquationIdVectorType& EquationId,ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        (rCurrentElement) -> CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);

        if(mIsDynamic)
        {
            (rCurrentElement) -> CalculateMassMatrix(mMassMatrix[thread],CurrentProcessInfo);

            (rCurrentElement) -> CalculateDampingMatrix(mDampingMatrix[thread],CurrentProcessInfo);

            AddDynamicsToLHS (LHS_Contribution, mMassMatrix[thread], mDampingMatrix[thread], CurrentProcessInfo);

            AddDynamicsToRHS (rCurrentElement, RHS_Contribution, mMassMatrix[thread], mDampingMatrix[thread], CurrentProcessInfo);
        }

        (rCurrentElement) -> EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Condition_CalculateSystemContributions(Condition::Pointer rCurrentCondition,LocalSystemMatrixType& LHS_Contribution,LocalSystemVectorType& RHS_Contribution,
                                                Element::EquationIdVectorType& EquationId,ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY

        (rCurrentCondition) -> CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);

        (rCurrentCondition) -> EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Update(ModelPart& r_model_part,DofsArrayType& rDofSet,TSystemMatrixType& A,TSystemVectorType& Dx,TSystemVectorType& b )
    {
        KRATOS_TRY

        //Update Displacement and Pressure (DOFs)
        for (typename DofsArrayType::iterator i_dof = rDofSet.begin(); i_dof != rDofSet.end(); ++i_dof)
            if (i_dof->IsFree())
                i_dof->GetSolutionStepValue() += Dx[i_dof->EquationId()];

        //Update Acceleration and Velocity
        array_1d<double,3> DeltaDisplacement;
        double DeltaTime = r_model_part.GetProcessInfo()[DELTA_TIME];

        for (ModelPart::NodeIterator i = r_model_part.NodesBegin(); i != r_model_part.NodesEnd(); ++i)
        {
            noalias(DeltaDisplacement)               = (i)->FastGetSolutionStepValue(DISPLACEMENT) - (i)->FastGetSolutionStepValue(DISPLACEMENT, 1);
            array_1d<double,3>& CurrentAcceleration  = (i)->FastGetSolutionStepValue(ACCELERATION);
            array_1d<double,3>& PreviousAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION, 1);
            array_1d<double,3>& CurrentVelocity      = (i)->FastGetSolutionStepValue(VELOCITY);
            array_1d<double,3>& PreviousVelocity     = (i)->FastGetSolutionStepValue(VELOCITY, 1);

            UpdateAcceleration(CurrentAcceleration, DeltaDisplacement, PreviousVelocity, PreviousAcceleration, DeltaTime);
            UpdateVelocity(CurrentVelocity, PreviousVelocity, PreviousAcceleration, CurrentAcceleration, DeltaTime);
        }

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Calculate_RHS_Contribution(Element::Pointer rCurrentElement,LocalSystemVectorType& RHS_Contribution,Element::EquationIdVectorType& EquationId,
                                    ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        (rCurrentElement) -> CalculateRightHandSide(RHS_Contribution,CurrentProcessInfo);

        if(mIsDynamic)
        {
            (rCurrentElement) -> CalculateMassMatrix(mMassMatrix[thread], CurrentProcessInfo);

            (rCurrentElement) -> CalculateDampingMatrix(mDampingMatrix[thread],CurrentProcessInfo);

            AddDynamicsToRHS (rCurrentElement, RHS_Contribution, mMassMatrix[thread], mDampingMatrix[thread], CurrentProcessInfo);
        }

        (rCurrentElement) -> EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Condition_Calculate_RHS_Contribution(Condition::Pointer rCurrentCondition,LocalSystemVectorType& RHS_Contribution,
                                            Element::EquationIdVectorType& EquationId,ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY

        (rCurrentCondition) -> CalculateRightHandSide(RHS_Contribution, CurrentProcessInfo);

        (rCurrentCondition) -> EquationIdVector(EquationId, CurrentProcessInfo);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

    bool mIsDynamic;

    double mBeta, mGamma;

    std::vector< Matrix > mMassMatrix;

    std::vector< Vector > mAccelerationVector;

    std::vector< Matrix > mDampingMatrix;

    std::vector< Vector > mVelocityVector;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void AddDynamicsToLHS(LocalSystemMatrixType& LHS_Contribution,LocalSystemMatrixType& M,LocalSystemMatrixType& C,ProcessInfo& CurrentProcessInfo)
    {
        double DeltaTime = CurrentProcessInfo[DELTA_TIME];

        // adding mass contribution
        if (M.size1() != 0)
            noalias(LHS_Contribution) += 1/(mBeta*DeltaTime*DeltaTime)*M;

        // adding damping contribution
        if (C.size1() != 0)
            noalias(LHS_Contribution) += mGamma/(mBeta*DeltaTime)*C;
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    inline void UpdateVelocity(array_1d<double,3>& CurrentVelocity,const array_1d<double,3>& PreviousVelocity,const array_1d<double,3>& PreviousAcceleration,
                                const array_1d<double,3>& CurrentAcceleration, double DeltaTime)
    {
        noalias(CurrentVelocity) = PreviousVelocity + (1-mGamma)*DeltaTime*PreviousAcceleration + mGamma*DeltaTime*CurrentAcceleration;
    }

    //------------------------------------------------------------------------------------

    inline void UpdateAcceleration(array_1d<double,3>& CurrentAcceleration,const array_1d<double,3>& DeltaDisplacement,const array_1d<double,3>& PreviousVelocity,
                                    const array_1d<double,3>& PreviousAcceleration, double DeltaTime)
    {
        noalias(CurrentAcceleration) = 1/(mBeta*DeltaTime*DeltaTime)*(DeltaDisplacement - DeltaTime*PreviousVelocity - (0.5-mBeta)*DeltaTime*DeltaTime*PreviousAcceleration);
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void AddDynamicsToRHS(Element::Pointer rCurrentElement,LocalSystemVectorType& RHS_Contribution,
                            LocalSystemMatrixType& M,LocalSystemMatrixType& C,ProcessInfo& CurrentProcessInfo)
    {
        int thread = OpenMPUtils::ThisThread();

        //adding inertia contribution
        if (M.size1() != 0)
        {
            rCurrentElement->GetSecondDerivativesVector(mAccelerationVector[thread], 0);

            noalias(RHS_Contribution) -= prod(M, mAccelerationVector[thread]);
        }

        //adding damping contribution
        if (C.size1() != 0)
        {
            rCurrentElement->GetFirstDerivativesVector(mVelocityVector[thread], 0);

            noalias(RHS_Contribution) -= prod(C, mVelocityVector[thread]);
        }
    }

}; // Class NewmarkScheme
}  // namespace Kratos

#endif // KRATOS_NEWMARK_SCHEME defined
