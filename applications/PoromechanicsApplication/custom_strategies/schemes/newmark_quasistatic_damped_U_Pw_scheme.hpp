//
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_NEWMARK_QUASISTATIC_DAMPED_U_PW_SCHEME )
#define  KRATOS_NEWMARK_QUASISTATIC_DAMPED_U_PW_SCHEME

// Application includes
#include "custom_strategies/schemes/newmark_quasistatic_U_Pw_scheme.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

template<class TSparseSpace, class TDenseSpace>

class NewmarkQuasistaticDampedUPwScheme : public NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( NewmarkQuasistaticDampedUPwScheme );

    typedef Scheme<TSparseSpace,TDenseSpace>                      BaseType;
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
    using NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>::mDeltaTime;
    using NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>::mBeta;
    using NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>::mGamma;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ///Constructor
    NewmarkQuasistaticDampedUPwScheme(double beta, double gamma, double theta)
        : NewmarkQuasistaticUPwScheme<TSparseSpace,TDenseSpace>(beta, gamma, theta)
    {
        //Allocate auxiliary memory
        int NumThreads = ParallelUtilities::GetNumThreads();
        mDampingMatrix.resize(NumThreads);
        mVelocityVector.resize(NumThreads);
    }

    //------------------------------------------------------------------------------------

    ///Destructor
    ~NewmarkQuasistaticDampedUPwScheme() override {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Note: this is in a parallel loop

    void CalculateSystemContributions(
        Element& rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        rCurrentElement.CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);

        rCurrentElement.CalculateDampingMatrix(mDampingMatrix[thread],CurrentProcessInfo);

        this->AddDampingToLHS (LHS_Contribution, mDampingMatrix[thread], CurrentProcessInfo);

        this->AddDampingToRHS (rCurrentElement, RHS_Contribution, mDampingMatrix[thread], CurrentProcessInfo);

        rCurrentElement.EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Note: this is in a parallel loop

    void CalculateRHSContribution(
        Element& rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        rCurrentElement.CalculateRightHandSide(RHS_Contribution,CurrentProcessInfo);

        rCurrentElement.CalculateDampingMatrix(mDampingMatrix[thread],CurrentProcessInfo);

        this->AddDampingToRHS (rCurrentElement, RHS_Contribution, mDampingMatrix[thread], CurrentProcessInfo);

        rCurrentElement.EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Note: this is in a parallel loop

    void CalculateLHSContribution(
        Element& rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        rCurrentElement.CalculateLeftHandSide(LHS_Contribution,CurrentProcessInfo);

        rCurrentElement.CalculateDampingMatrix(mDampingMatrix[thread],CurrentProcessInfo);

        this->AddDampingToLHS (LHS_Contribution, mDampingMatrix[thread], CurrentProcessInfo);

        rCurrentElement.EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

    std::vector< Matrix > mDampingMatrix;
    std::vector< Vector > mVelocityVector;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void AddDampingToLHS(LocalSystemMatrixType& LHS_Contribution,LocalSystemMatrixType& C,const ProcessInfo& CurrentProcessInfo)
    {
        // adding damping contribution
        if (C.size1() != 0)
            noalias(LHS_Contribution) += mGamma/(mBeta*mDeltaTime)*C;
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void AddDampingToRHS(Element& rCurrentElement,LocalSystemVectorType& RHS_Contribution,
                            LocalSystemMatrixType& C,const ProcessInfo& CurrentProcessInfo)
    {
        int thread = OpenMPUtils::ThisThread();

        //adding damping contribution
        if (C.size1() != 0)
        {
            rCurrentElement.GetFirstDerivativesVector(mVelocityVector[thread], 0);

            noalias(RHS_Contribution) -= prod(C, mVelocityVector[thread]);
        }
    }

}; // Class NewmarkQuasistaticDampedUPwScheme
}  // namespace Kratos

#endif // KRATOS_NEWMARK_QUASISTATIC_DAMPED_U_PW_SCHEME defined
