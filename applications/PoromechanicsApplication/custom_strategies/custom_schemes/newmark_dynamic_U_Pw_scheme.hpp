//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_NEWMARK_DYNAMIC_U_PW_SCHEME )
#define  KRATOS_NEWMARK_DYNAMIC_U_PW_SCHEME

// Application includes
#include "custom_strategies/custom_schemes/newmark_quasistatic_U_Pw_scheme.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

template<class TSparseSpace, class TDenseSpace>

class NewmarkDynamicUPwScheme : public NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( NewmarkDynamicUPwScheme );

    typedef Scheme<TSparseSpace,TDenseSpace>                      BaseType;
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
    using NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>::mDeltaTime;
    using NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>::mBeta;
    using NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>::mGamma;
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    //Constructor
    NewmarkDynamicUPwScheme(ModelPart& r_model_part, double beta, double gamma, double theta,
                            double rayleigh_m, double rayleigh_k) : NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>(r_model_part, beta, gamma, theta)
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
    
    //Destructor
    virtual ~NewmarkDynamicUPwScheme() {}
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    int Check(ModelPart& r_model_part)
    {
        KRATOS_TRY
        
        int err = NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>::Check(r_model_part);
        if(err != 0) return err;

        if ( RAYLEIGH_ALPHA.Key() == 0 )
            KRATOS_THROW_ERROR( std::invalid_argument, "RAYLEIGH_ALPHA has Key zero! (check if the application is correctly registered", "" )
        if ( RAYLEIGH_BETA.Key() == 0 )
            KRATOS_THROW_ERROR( std::invalid_argument, "RAYLEIGH_BETA has Key zero! (check if the application is correctly registered", "" )
            
        // Check rayleigh coefficients
        if( mRayleighAlpha < 0.0 || mRayleighBeta < 0.0 )
            KRATOS_THROW_ERROR( std::invalid_argument,"Some of the rayleigh coefficients has an invalid value ", "" )
            
        return err;
        
        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Initialize(ModelPart& r_model_part)
    {
        KRATOS_TRY
        
        NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>::Initialize(r_model_part);
        
        r_model_part.GetProcessInfo()[RAYLEIGH_ALPHA] = mRayleighAlpha;
        r_model_part.GetProcessInfo()[RAYLEIGH_BETA] = mRayleighBeta;
        
        KRATOS_CATCH("")
    }
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateSystemContributions(Element::Pointer rCurrentElement,LocalSystemMatrixType& LHS_Contribution,LocalSystemVectorType& RHS_Contribution,
                                        Element::EquationIdVectorType& EquationId,ProcessInfo& CurrentProcessInfo)
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

    void Calculate_RHS_Contribution(Element::Pointer rCurrentElement,LocalSystemVectorType& RHS_Contribution,Element::EquationIdVectorType& EquationId,
                                    ProcessInfo& CurrentProcessInfo)
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

protected:
    
    // Member Variables
    
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
