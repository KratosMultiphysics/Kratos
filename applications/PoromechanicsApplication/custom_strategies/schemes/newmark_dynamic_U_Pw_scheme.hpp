//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_NEWMARK_DYNAMIC_U_PW_SCHEME )
#define  KRATOS_NEWMARK_DYNAMIC_U_PW_SCHEME

// Application includes
#include "custom_strategies/schemes/newmark_quasistatic_U_Pw_scheme.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

template<class TSparseSpace, class TDenseSpace>

class NewmarkDynamicUPwScheme : public NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( NewmarkDynamicUPwScheme );

    typedef Scheme<TSparseSpace,TDenseSpace>                      BaseType;
    typedef typename BaseType::DofsArrayType                 DofsArrayType;
    typedef typename BaseType::TSystemMatrixType         TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType         TSystemVectorType;
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
    using NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>::mDeltaTime;
    using NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>::mBeta;
    using NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>::mGamma;
    using NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>::mTheta;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ///Constructor
    NewmarkDynamicUPwScheme(double beta, double gamma, double theta,
                            double rayleigh_m, double rayleigh_k) : NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>(beta, gamma, theta)
    {
        mRayleighAlpha = rayleigh_m;
        mRayleighBeta = rayleigh_k;
                
        //Allocate auxiliary memory
        int NumThreads = OpenMPUtils::GetNumThreads();
        mMassMatrix.resize(NumThreads);
        mAccelerationVector.resize(NumThreads);
        mDampingMatrix.resize(NumThreads);
        mVelocityVector.resize(NumThreads);
    }

    //------------------------------------------------------------------------------------
    
    ///Destructor
    ~NewmarkDynamicUPwScheme() override {}
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    int Check(ModelPart& r_model_part) override
    {
        KRATOS_TRY
        
        int ierr = NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>::Check(r_model_part);
        if(ierr != 0) return ierr;

        if ( RAYLEIGH_ALPHA.Key() == 0 )
            KRATOS_THROW_ERROR( std::invalid_argument, "RAYLEIGH_ALPHA Key is 0. Check if all applications were correctly registered.", "" )
        if ( RAYLEIGH_BETA.Key() == 0 )
            KRATOS_THROW_ERROR( std::invalid_argument, "RAYLEIGH_BETA Key is 0. Check if all applications were correctly registered.", "" )
            
        // Check rayleigh coefficients
        if( mRayleighAlpha < 0.0 || mRayleighBeta < 0.0 )
            KRATOS_THROW_ERROR( std::invalid_argument,"Some of the rayleigh coefficients has an invalid value ", "" )
            
        return ierr;
        
        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Initialize(ModelPart& r_model_part) override
    {
        KRATOS_TRY
        
        NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>::Initialize(r_model_part);
        
        r_model_part.GetProcessInfo()[RAYLEIGH_ALPHA] = mRayleighAlpha;
        r_model_part.GetProcessInfo()[RAYLEIGH_BETA] = mRayleighBeta;
        
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
        // Predict Displacements on free nodes and update Acceleration, Velocity and DtPressure
        
        array_1d<double,3> DeltaDisplacement;
        double DeltaPressure;
        
        const int NNodes = static_cast<int>(r_model_part.Nodes().size());
        ModelPart::NodesContainerType::iterator node_begin = r_model_part.NodesBegin();
        
        #pragma omp parallel for private(DeltaDisplacement,DeltaPressure)
        for(int i = 0; i < NNodes; i++)
        {
            ModelPart::NodesContainerType::iterator itNode = node_begin + i;
            
            array_1d<double,3>& CurrentDisplacement = itNode->FastGetSolutionStepValue(DISPLACEMENT);
            array_1d<double,3>& CurrentAcceleration = itNode->FastGetSolutionStepValue(ACCELERATION);
            array_1d<double,3>& CurrentVelocity = itNode->FastGetSolutionStepValue(VELOCITY);
            
            const array_1d<double,3>& PreviousDisplacement = itNode->FastGetSolutionStepValue(DISPLACEMENT, 1);
            const array_1d<double,3>& PreviousAcceleration = itNode->FastGetSolutionStepValue(ACCELERATION, 1);
            const array_1d<double,3>& PreviousVelocity = itNode->FastGetSolutionStepValue(VELOCITY, 1);

            if (itNode -> IsFixed(ACCELERATION_X))
            {
                CurrentDisplacement[0] = PreviousDisplacement[0] + mDeltaTime * PreviousVelocity[0] + std::pow(mDeltaTime, 2) * ( ( 0.5 - mBeta) * PreviousAcceleration[0] + mBeta * CurrentAcceleration[0] );
            }
            else if (itNode -> IsFixed(VELOCITY_X))
            {
                CurrentDisplacement[0] = PreviousDisplacement[0] + mDeltaTime*(mBeta/mGamma*(CurrentVelocity[0]-PreviousVelocity[0])+PreviousVelocity[0]);
            }
            else if (itNode -> IsFixed(DISPLACEMENT_X) == false)
            {
                CurrentDisplacement[0] = PreviousDisplacement[0] + mDeltaTime * PreviousVelocity[0] + 0.5 * std::pow(mDeltaTime, 2) * PreviousAcceleration[0];
            }

            if (itNode -> IsFixed(ACCELERATION_Y))
            {
                CurrentDisplacement[1] = PreviousDisplacement[1] + mDeltaTime * PreviousVelocity[1] + std::pow(mDeltaTime, 2) * ( ( 0.5 - mBeta) * PreviousAcceleration[1] + mBeta * CurrentAcceleration[1] );
            }
            else if (itNode -> IsFixed(VELOCITY_Y))
            {
                CurrentDisplacement[1] = PreviousDisplacement[1] + mDeltaTime*(mBeta/mGamma*(CurrentVelocity[1]-PreviousVelocity[1])+PreviousVelocity[1]);
            }
            else if (itNode -> IsFixed(DISPLACEMENT_Y) == false)
            {
                CurrentDisplacement[1] = PreviousDisplacement[1] + mDeltaTime * PreviousVelocity[1] + 0.5 * std::pow(mDeltaTime, 2) * PreviousAcceleration[1];
            }

            // For 3D cases
            if (itNode -> HasDofFor(DISPLACEMENT_Z))
            {
                if (itNode -> IsFixed(ACCELERATION_Z))
                {
                    CurrentDisplacement[2] = PreviousDisplacement[2] + mDeltaTime * PreviousVelocity[2] + std::pow(mDeltaTime, 2) * ( ( 0.5 - mBeta) * PreviousAcceleration[2] + mBeta * CurrentAcceleration[2] );
                }
                else if (itNode -> IsFixed(VELOCITY_Z))
                {
                    CurrentDisplacement[2] = PreviousDisplacement[2] + mDeltaTime*(mBeta/mGamma*(CurrentVelocity[2]-PreviousVelocity[2])+PreviousVelocity[2]);
                }
                else if (itNode -> IsFixed(DISPLACEMENT_Z) == false)
                {
                    CurrentDisplacement[2] = PreviousDisplacement[2] + mDeltaTime * PreviousVelocity[2] + 0.5 * std::pow(mDeltaTime, 2) * PreviousAcceleration[2];
                }
            }
            
            noalias(DeltaDisplacement) = CurrentDisplacement - PreviousDisplacement;
            
            noalias(CurrentAcceleration) = 1.0/(mBeta*mDeltaTime*mDeltaTime)*(DeltaDisplacement - mDeltaTime*PreviousVelocity - (0.5-mBeta)*mDeltaTime*mDeltaTime*PreviousAcceleration);
            noalias(CurrentVelocity) = PreviousVelocity + (1.0-mGamma)*mDeltaTime*PreviousAcceleration + mGamma*mDeltaTime*CurrentAcceleration;
            
            double& CurrentDtPressure = itNode->FastGetSolutionStepValue(DT_WATER_PRESSURE);
            DeltaPressure = itNode->FastGetSolutionStepValue(WATER_PRESSURE) - itNode->FastGetSolutionStepValue(WATER_PRESSURE, 1);
            const double& PreviousDtPressure = itNode->FastGetSolutionStepValue(DT_WATER_PRESSURE, 1);

            CurrentDtPressure = 1.0/(mTheta*mDeltaTime)*(DeltaPressure - (1.0-mTheta)*mDeltaTime*PreviousDtPressure);
        }
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

        int thread = OpenMPUtils::ThisThread();

        (rCurrentElement) -> CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);

        (rCurrentElement) -> CalculateMassMatrix(mMassMatrix[thread],CurrentProcessInfo);

        (rCurrentElement) -> CalculateDampingMatrix(mDampingMatrix[thread],CurrentProcessInfo);

        this->AddDynamicsToLHS (LHS_Contribution, mMassMatrix[thread], mDampingMatrix[thread], CurrentProcessInfo);

        this->AddDynamicsToRHS (rCurrentElement, RHS_Contribution, mMassMatrix[thread], mDampingMatrix[thread], CurrentProcessInfo);

        (rCurrentElement) -> EquationIdVector(EquationId,CurrentProcessInfo);

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

        int thread = OpenMPUtils::ThisThread();

        (rCurrentElement) -> CalculateRightHandSide(RHS_Contribution,CurrentProcessInfo);

        (rCurrentElement) -> CalculateMassMatrix(mMassMatrix[thread], CurrentProcessInfo);

        (rCurrentElement) -> CalculateDampingMatrix(mDampingMatrix[thread],CurrentProcessInfo);

        this->AddDynamicsToRHS (rCurrentElement, RHS_Contribution, mMassMatrix[thread], mDampingMatrix[thread], CurrentProcessInfo);

        (rCurrentElement) -> EquationIdVector(EquationId,CurrentProcessInfo);

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

        int thread = OpenMPUtils::ThisThread();

        (rCurrentElement) -> CalculateLeftHandSide(LHS_Contribution,CurrentProcessInfo);

        (rCurrentElement) -> CalculateMassMatrix(mMassMatrix[thread],CurrentProcessInfo);

        (rCurrentElement) -> CalculateDampingMatrix(mDampingMatrix[thread],CurrentProcessInfo);

        this->AddDynamicsToLHS (LHS_Contribution, mMassMatrix[thread], mDampingMatrix[thread], CurrentProcessInfo);

        (rCurrentElement) -> EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:
    
    /// Member Variables
    
    double mRayleighAlpha;
    double mRayleighBeta;
    
    std::vector< Matrix > mMassMatrix;
    std::vector< Vector > mAccelerationVector;
    std::vector< Matrix > mDampingMatrix;
    std::vector< Vector > mVelocityVector;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void AddDynamicsToLHS(LocalSystemMatrixType& LHS_Contribution,LocalSystemMatrixType& M,LocalSystemMatrixType& C,ProcessInfo& CurrentProcessInfo)
    {        
        // adding mass contribution
        if (M.size1() != 0)
            noalias(LHS_Contribution) += 1.0/(mBeta*mDeltaTime*mDeltaTime)*M;

        // adding damping contribution
        if (C.size1() != 0)
            noalias(LHS_Contribution) += mGamma/(mBeta*mDeltaTime)*C;
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

}; // Class NewmarkDynamicUPwScheme
}  // namespace Kratos

#endif // KRATOS_NEWMARK_DYNAMIC_U_PW_SCHEME defined
