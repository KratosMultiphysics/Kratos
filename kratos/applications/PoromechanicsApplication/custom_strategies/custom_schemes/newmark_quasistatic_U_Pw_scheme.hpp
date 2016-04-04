//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_NEWMARK_QUASISTATIC_U_PW_SCHEME )
#define  KRATOS_NEWMARK_QUASISTATIC_U_PW_SCHEME

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"

// Application includes
#include "poromechanics_application_variables.h"

namespace Kratos
{

template<class TSparseSpace, class TDenseSpace>

class NewmarkQuasistaticUPwScheme : public Scheme<TSparseSpace,TDenseSpace>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( NewmarkQuasistaticUPwScheme );

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
    NewmarkQuasistaticUPwScheme(ModelPart& r_model_part, double beta, double gamma, double theta) : Scheme<TSparseSpace,TDenseSpace>()
    {   
        mDeltaTime = r_model_part.GetProcessInfo()[DELTA_TIME];
        mBeta = beta;
        mGamma = gamma;
        mTheta = theta;
        r_model_part.GetProcessInfo()[NEWMARK_COEFFICIENT_U] = gamma/(beta*mDeltaTime);
        r_model_part.GetProcessInfo()[NEWMARK_COEFFICIENT_P] = 1.0/(theta*mDeltaTime);
    }

    //------------------------------------------------------------------------------------
    
    //Destructor
    virtual ~NewmarkQuasistaticUPwScheme() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    int Check(ModelPart& r_model_part)
    {
        KRATOS_TRY
        
        //check for variables keys (verify that the variables are correctly initialized)
        if(DISPLACEMENT.Key() == 0)
            KRATOS_THROW_ERROR( std::invalid_argument,"DISPLACEMENT has Key zero! (check if the application is correctly registered", "" )
        if(VELOCITY.Key() == 0)
            KRATOS_THROW_ERROR( std::invalid_argument,"VELOCITY has Key zero! (check if the application is correctly registered", "" )
        if(ACCELERATION.Key() == 0)
            KRATOS_THROW_ERROR( std::invalid_argument,"ACCELERATION has Key zero! (check if the application is correctly registered", "" )
            
        if(WATER_PRESSURE.Key() == 0)
            KRATOS_THROW_ERROR( std::invalid_argument, "WATER_PRESSURE has Key zero! (check if the application is correctly registered", "" )
        if(DT_WATER_PRESSURE.Key() == 0)
            KRATOS_THROW_ERROR( std::invalid_argument, "DT_WATER_PRESSURE has Key zero! (check if the application is correctly registered", "" )

        //check that variables are correctly allocated
        for(ModelPart::NodesContainerType::iterator it=r_model_part.NodesBegin(); it!=r_model_part.NodesEnd(); it++)
        {
            if(it->SolutionStepsDataHas(DISPLACEMENT) == false)
                KRATOS_THROW_ERROR( std::logic_error, "DISPLACEMENT variable is not allocated for node ", it->Id() )
            if(it->SolutionStepsDataHas(VELOCITY) == false)
                KRATOS_THROW_ERROR( std::logic_error, "DISPLACEMENT variable is not allocated for node ", it->Id() )
            if(it->SolutionStepsDataHas(ACCELERATION) == false)
                KRATOS_THROW_ERROR( std::logic_error, "DISPLACEMENT variable is not allocated for node ", it->Id() )
            if(it->SolutionStepsDataHas(WATER_PRESSURE) == false)
                KRATOS_THROW_ERROR( std::logic_error, "WATER_PRESSURE variable is not allocated for node ", it->Id() )
            if(it->SolutionStepsDataHas(DT_WATER_PRESSURE) == false)
                KRATOS_THROW_ERROR( std::logic_error, "DT_WATER_PRESSURE variable is not allocated for node ", it->Id() )

            if(it->HasDofFor(DISPLACEMENT_X) == false)
                KRATOS_THROW_ERROR( std::invalid_argument,"missing DISPLACEMENT_X dof on node ",it->Id() )
            if(it->HasDofFor(DISPLACEMENT_Y) == false)
                KRATOS_THROW_ERROR( std::invalid_argument,"missing DISPLACEMENT_Y dof on node ",it->Id() )
            if(it->HasDofFor(DISPLACEMENT_Z) == false)
                KRATOS_THROW_ERROR( std::invalid_argument,"missing DISPLACEMENT_Z dof on node ",it->Id() )
            if(it->HasDofFor(WATER_PRESSURE) == false)
                KRATOS_THROW_ERROR( std::invalid_argument,"missing WATER_PRESSURE dof on node ",it->Id() )
        }

        //check for minimum value of the buffer index.
        if (r_model_part.GetBufferSize() < 2)
            KRATOS_THROW_ERROR( std::logic_error, "insufficient buffer size. Buffer size should be greater than 2. Current size is", r_model_part.GetBufferSize() )
        
        // Check DeltaTime
        if (mDeltaTime < 1.0e-25)
            KRATOS_THROW_ERROR( std::logic_error, "Detected DELTA_TIME < 1e-25 in the Solution Scheme. DELTA_TIME should be greater than 0.0", "" )
        
        // Check beta, gamma and theta
        if(mBeta <= 0.0 || mGamma<= 0.0 || mTheta <= 0.0)
            KRATOS_THROW_ERROR( std::invalid_argument,"Some of the scheme variables: beta, gamma or theta has an invalid value ", "" )
            
        return 0;
        
        KRATOS_CATCH( "" )
    }
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Predict(ModelPart& r_model_part,DofsArrayType& rDofSet,TSystemMatrixType& A,TSystemVectorType& Dx,TSystemVectorType& b)
    {
        this->UpdateVariablesDerivatives(r_model_part);
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

    void CalculateSystemContributions(Element::Pointer rCurrentElement,LocalSystemMatrixType& LHS_Contribution,LocalSystemVectorType& RHS_Contribution,
                                        Element::EquationIdVectorType& EquationId,ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY

        (rCurrentElement) -> CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);

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
        {
            if (i_dof->IsFree())
                i_dof->GetSolutionStepValue() += Dx[i_dof->EquationId()];
        }

        this->UpdateVariablesDerivatives(r_model_part);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Calculate_RHS_Contribution(Element::Pointer rCurrentElement,LocalSystemVectorType& RHS_Contribution,Element::EquationIdVectorType& EquationId,
                                    ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY

        (rCurrentElement) -> CalculateRightHandSide(RHS_Contribution,CurrentProcessInfo);

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
        
    double mBeta;
    double mGamma;
    double mTheta;
    double mDeltaTime;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    inline void UpdateVariablesDerivatives(ModelPart& r_model_part)
    {
        KRATOS_TRY

        //Update Acceleration, Velocity and DtPressure
        
        array_1d<double,3> DeltaDisplacement;
        array_1d<double,3> PreviousAcceleration;
        array_1d<double,3> PreviousVelocity;
        double DeltaPressure;
        double PreviousDtPressure;
        
        for (ModelPart::NodeIterator i = r_model_part.NodesBegin();i != r_model_part.NodesEnd(); ++i)
        {
            array_1d<double,3>& CurrentAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION);
            array_1d<double,3>& CurrentVelocity = (i)->FastGetSolutionStepValue(VELOCITY);
            noalias(DeltaDisplacement) = (i)->FastGetSolutionStepValue(DISPLACEMENT) - (i)->FastGetSolutionStepValue(DISPLACEMENT, 1);
            noalias(PreviousAcceleration) = (i)->FastGetSolutionStepValue(ACCELERATION, 1);
            noalias(PreviousVelocity) = (i)->FastGetSolutionStepValue(VELOCITY, 1);
            
            noalias(CurrentAcceleration) = 1.0/(mBeta*mDeltaTime*mDeltaTime)*(DeltaDisplacement - mDeltaTime*PreviousVelocity - (0.5-mBeta)*mDeltaTime*mDeltaTime*PreviousAcceleration);
            noalias(CurrentVelocity) = PreviousVelocity + (1.0-mGamma)*mDeltaTime*PreviousAcceleration + mGamma*mDeltaTime*CurrentAcceleration;
            
            double& CurrentDtPressure = (i)->FastGetSolutionStepValue(DT_WATER_PRESSURE);
            DeltaPressure = (i)->FastGetSolutionStepValue(WATER_PRESSURE) - (i)->FastGetSolutionStepValue(WATER_PRESSURE, 1);
            PreviousDtPressure = (i)->FastGetSolutionStepValue(DT_WATER_PRESSURE, 1);

            CurrentDtPressure = 1.0/(mTheta*mDeltaTime)*(DeltaPressure - (1.0-mTheta)*mDeltaTime*PreviousDtPressure);
        }

        KRATOS_CATCH( "" )
    }

}; // Class NewmarkQuasistaticUPwScheme
}  // namespace Kratos

#endif // KRATOS_NEWMARK_QUASISTATIC_U_PW_SCHEME defined
