//
//   Project Name:        			KratosDamApplication $
//   Last Modified by:    $Author:    	  Lorenzo Gracia $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_DAM_P_SCHEME )
#define  KRATOS_DAM_P_SCHEME

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "includes/variables.h"
#include "containers/array_1d.h"
#include "includes/element.h"

// Application includes
#include "dam_application_variables.h"

namespace Kratos
{

template<class TSparseSpace, class TDenseSpace>

class DamPScheme : public Scheme<TSparseSpace,TDenseSpace>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( DamPScheme );

    typedef Scheme<TSparseSpace,TDenseSpace>                      BaseType;
    typedef typename BaseType::DofsArrayType                 DofsArrayType;
    typedef typename BaseType::TSystemMatrixType         TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType         TSystemVectorType;
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ///Constructor
    DamPScheme(double beta, double gamma): Scheme<TSparseSpace,TDenseSpace>()
    {
        mBeta = beta;
        mGamma = gamma;

    }

    //------------------------------------------------------------------------------------

    ///Destructor
    virtual ~DamPScheme() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    int Check(ModelPart& r_model_part) override
    {
        KRATOS_TRY

        //check for variables keys (verify that the variables are correctly initialized)
        if(PRESSURE.Key() == 0)
            KRATOS_THROW_ERROR( std::invalid_argument, "PRESSURE has Key zero! (check if the application is correctly registered", "" )
        if(Dt_PRESSURE.Key() == 0)
            KRATOS_THROW_ERROR( std::invalid_argument, "Dt_PRESSURE has Key zero! (check if the application is correctly registered", "" )
        if(Dt2_PRESSURE.Key() == 0)
            KRATOS_THROW_ERROR( std::invalid_argument, "Dt2_PRESSURE has Key zero! (check if the application is correctly registered", "" )
        if ( VELOCITY_PRESSURE_COEFFICIENT.Key() == 0 )
            KRATOS_THROW_ERROR( std::invalid_argument, "VELOCITY_PRESSURE_COEFFICIENT has Key zero! (check if the application is correctly registered", "" )
        if ( ACCELERATION_PRESSURE_COEFFICIENT.Key() == 0 )
            KRATOS_THROW_ERROR( std::invalid_argument, "ACCELERATION_PRESSURE_COEFFICIENT has Key zero! (check if the application is correctly registered", "" )

        //check that variables are correctly allocated
        for(ModelPart::NodesContainerType::iterator it=r_model_part.NodesBegin(); it!=r_model_part.NodesEnd(); it++)
        {
            if(it->SolutionStepsDataHas(PRESSURE) == false)
                KRATOS_THROW_ERROR( std::logic_error, "PRESSURE variable is not allocated for node ", it->Id() )
            if(it->SolutionStepsDataHas(Dt_PRESSURE) == false)
                KRATOS_THROW_ERROR( std::logic_error, "Dt_PRESSURE variable is not allocated for node ", it->Id() )
            if(it->SolutionStepsDataHas(Dt2_PRESSURE) == false)
                KRATOS_THROW_ERROR( std::logic_error, "Dt2_PRESSURE variable is not allocated for node ", it->Id() )

            if(it->HasDofFor(PRESSURE) == false)
                KRATOS_THROW_ERROR( std::invalid_argument,"missing PRESSURE dof on node ",it->Id() )
        }

        //check for minimum value of the buffer index.
        if (r_model_part.GetBufferSize() < 2)
            KRATOS_THROW_ERROR( std::logic_error, "insufficient buffer size. Buffer size should be greater than 2. Current size is", r_model_part.GetBufferSize() )

        // Check beta, gamma and theta
        if(mBeta <= 0.0 || mGamma<= 0.0)
            KRATOS_THROW_ERROR( std::invalid_argument,"Some of the scheme variables: beta or  gamma has an invalid value ", "" )

        return 0;

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Initialize(ModelPart& r_model_part) override
    {
        KRATOS_TRY

        mDeltaTime = r_model_part.GetProcessInfo()[DELTA_TIME];
        r_model_part.GetProcessInfo()[VELOCITY_PRESSURE_COEFFICIENT] = mGamma/(mBeta*mDeltaTime);
        r_model_part.GetProcessInfo()[ACCELERATION_PRESSURE_COEFFICIENT] = 1.0/(mBeta*mDeltaTime*mDeltaTime);

        BaseType::mSchemeIsInitialized = true;

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void InitializeSolutionStep(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        mDeltaTime = r_model_part.GetProcessInfo()[DELTA_TIME];
        r_model_part.GetProcessInfo()[VELOCITY_PRESSURE_COEFFICIENT] = mGamma/(mBeta*mDeltaTime);
        r_model_part.GetProcessInfo()[ACCELERATION_PRESSURE_COEFFICIENT] = 1.0/(mBeta*mDeltaTime*mDeltaTime);

        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        int NElems = static_cast<int>(r_model_part.Elements().size());
        ModelPart::ElementsContainerType::iterator el_begin = r_model_part.ElementsBegin();

        #pragma omp parallel for
        for(int i = 0; i < NElems; i++)
        {
            ModelPart::ElementsContainerType::iterator itElem = el_begin + i;
            itElem -> InitializeSolutionStep(CurrentProcessInfo);
        }

        int NCons = static_cast<int>(r_model_part.Conditions().size());
        ModelPart::ConditionsContainerType::iterator con_begin = r_model_part.ConditionsBegin();

        #pragma omp parallel for
        for(int i = 0; i < NCons; i++)
        {
            ModelPart::ConditionsContainerType::iterator itCond = con_begin + i;
            itCond -> InitializeSolutionStep(CurrentProcessInfo);
        }

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Predict(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {

      KRATOS_TRY

        // Updating  DtPressure and Dt2Pressure
        double DeltaPressure;

        const int NNodes = static_cast<int>(r_model_part.Nodes().size());
        ModelPart::NodesContainerType::iterator node_begin = r_model_part.NodesBegin();

        #pragma omp parallel for private(DeltaPressure)
        for(int i = 0; i < NNodes; i++)
        {
            ModelPart::NodesContainerType::iterator itNode = node_begin + i;

            // Terms related to Pressure field
            double& CurrentDt2Pressure = itNode->FastGetSolutionStepValue(Dt2_PRESSURE);
            double& CurrentDtPressure = itNode->FastGetSolutionStepValue(Dt_PRESSURE);
            DeltaPressure = itNode->FastGetSolutionStepValue(PRESSURE) - itNode->FastGetSolutionStepValue(PRESSURE, 1);
            const double& PreviousDt2Pressure = itNode->FastGetSolutionStepValue(Dt2_PRESSURE, 1);
            const double& PreviousDtPressure = itNode->FastGetSolutionStepValue(Dt_PRESSURE, 1);

            CurrentDt2Pressure = 1.0/(mBeta*mDeltaTime*mDeltaTime)*(DeltaPressure - mDeltaTime*PreviousDtPressure - (0.5-mBeta)*mDeltaTime*mDeltaTime*PreviousDt2Pressure);
            CurrentDtPressure = PreviousDtPressure + (1.0-mGamma)*mDeltaTime*PreviousDt2Pressure + mGamma*mDeltaTime*CurrentDt2Pressure;

        }

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void InitializeNonLinIteration(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        int NElems = static_cast<int>(r_model_part.Elements().size());
        ModelPart::ElementsContainerType::iterator el_begin = r_model_part.ElementsBegin();

        #pragma omp parallel for
        for(int i = 0; i < NElems; i++)
        {
            ModelPart::ElementsContainerType::iterator itElem = el_begin + i;
            itElem -> InitializeNonLinearIteration(CurrentProcessInfo);
        }

        int NCons = static_cast<int>(r_model_part.Conditions().size());
        ModelPart::ConditionsContainerType::iterator con_begin = r_model_part.ConditionsBegin();

        #pragma omp parallel for
        for(int i = 0; i < NCons; i++)
        {
            ModelPart::ConditionsContainerType::iterator itCond = con_begin + i;
            itCond -> InitializeNonLinearIteration(CurrentProcessInfo);
        }

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void FinalizeNonLinIteration(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        int NElems = static_cast<int>(r_model_part.Elements().size());
        ModelPart::ElementsContainerType::iterator el_begin = r_model_part.ElementsBegin();

        #pragma omp parallel for
        for(int i = 0; i < NElems; i++)
        {
            ModelPart::ElementsContainerType::iterator itElem = el_begin + i;
            itElem -> FinalizeNonLinearIteration(CurrentProcessInfo);
        }

        int NCons = static_cast<int>(r_model_part.Conditions().size());
        ModelPart::ConditionsContainerType::iterator con_begin = r_model_part.ConditionsBegin();

        #pragma omp parallel for
        for(int i = 0; i < NCons; i++)
        {
            ModelPart::ConditionsContainerType::iterator itCond = con_begin + i;
            itCond -> FinalizeNonLinearIteration(CurrentProcessInfo);
        }

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Note: this is in a parallel loop

    void CalculateSystemContributions(
        Element::Pointer rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        (rCurrentElement) -> CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);

        (rCurrentElement) -> EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Note: this is in a parallel loop

    void Condition_CalculateSystemContributions(
        Condition::Pointer rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        (rCurrentCondition) -> CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);

        (rCurrentCondition) -> EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Note: this is in a parallel loop

    void Calculate_RHS_Contribution(
        Element::Pointer rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        (rCurrentElement) -> CalculateRightHandSide(RHS_Contribution,CurrentProcessInfo);

        (rCurrentElement) -> EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Note: this is in a parallel loop

    void Condition_Calculate_RHS_Contribution(
        Condition::Pointer rCurrentCondition,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        (rCurrentCondition) -> CalculateRightHandSide(RHS_Contribution, CurrentProcessInfo);

        (rCurrentCondition) -> EquationIdVector(EquationId, CurrentProcessInfo);


        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Note: this is in a parallel loop

    void Calculate_LHS_Contribution(
        Element::Pointer rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        (rCurrentElement) -> CalculateLeftHandSide(LHS_Contribution, CurrentProcessInfo);

        (rCurrentElement) -> EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Note: this is in a parallel loop

    void Condition_Calculate_LHS_Contribution(
        Condition::Pointer rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        (rCurrentCondition) -> CalculateLeftHandSide(LHS_Contribution, CurrentProcessInfo);

        (rCurrentCondition) -> EquationIdVector(EquationId, CurrentProcessInfo);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Update(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector DofSetPartition;
        OpenMPUtils::DivideInPartitions(rDofSet.size(), NumThreads, DofSetPartition);

        #pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();

            typename DofsArrayType::iterator DofsBegin = rDofSet.begin() + DofSetPartition[k];
            typename DofsArrayType::iterator DofsEnd = rDofSet.begin() + DofSetPartition[k+1];

            //Update Displacement and Pressure (DOFs)
            for (typename DofsArrayType::iterator itDof = DofsBegin; itDof != DofsEnd; ++itDof)
            {
                if (itDof->IsFree())
                    itDof->GetSolutionStepValue() += TSparseSpace::GetValue(Dx, itDof->EquationId());
            }
        }

        this->UpdateVariablesDerivatives(r_model_part);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

    double mBeta;
    double mGamma;
    double mDeltaTime;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    inline void UpdateVariablesDerivatives(ModelPart& r_model_part)
    {
        KRATOS_TRY

        // Updating DtPressure and Dt2Pressure
        double DeltaPressure;

        const int NNodes = static_cast<int>(r_model_part.Nodes().size());
        ModelPart::NodesContainerType::iterator node_begin = r_model_part.NodesBegin();

        #pragma omp parallel for private(DeltaPressure)
        for(int i = 0; i < NNodes; i++)
        {
            ModelPart::NodesContainerType::iterator itNode = node_begin + i;

            // Terms related to Pressure field

            double& CurrentDt2Pressure = itNode->FastGetSolutionStepValue(Dt2_PRESSURE);
            double& CurrentDtPressure = itNode->FastGetSolutionStepValue(Dt_PRESSURE);
            DeltaPressure = itNode->FastGetSolutionStepValue(PRESSURE) - itNode->FastGetSolutionStepValue(PRESSURE, 1);
            const double& PreviousDt2Pressure = itNode->FastGetSolutionStepValue(Dt2_PRESSURE, 1);
            const double& PreviousDtPressure = itNode->FastGetSolutionStepValue(Dt_PRESSURE, 1);

            CurrentDt2Pressure = 1.0/(mBeta*mDeltaTime*mDeltaTime)*(DeltaPressure - mDeltaTime*PreviousDtPressure - (0.5-mBeta)*mDeltaTime*mDeltaTime*PreviousDt2Pressure);
            CurrentDtPressure = PreviousDtPressure + (1.0-mGamma)*mDeltaTime*PreviousDt2Pressure + mGamma*mDeltaTime*CurrentDt2Pressure;

        }

        KRATOS_CATCH( "" )
    }


}; // Class DamPScheme
}  // namespace Kratos

#endif // KRATOS_DAM_P_SCHEME defined

