//
//   Project Name:        KratosPoromechanicsApplication $
//   Author:              $Author:    Ignasi de Pouplana $
//   Last Modified by:    $Author:    Danilo Cavalcanti  $
//   Date:                $Date:               June 2023 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_PORO_NEWMARK_QUASISTATIC_DAMPED_U_PL_PG_SCHEME )
#define  KRATOS_PORO_NEWMARK_QUASISTATIC_DAMPED_U_PL_PG_SCHEME

// Application includes
#include "custom_strategies/schemes/two-phase_flow/poro_newmark_quasistatic_U_Pl_Pg_scheme.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

template<class TSparseSpace, class TDenseSpace>

class PoroNewmarkQuasistaticDampedUPlPgScheme : public PoroNewmarkQuasistaticUPlPgScheme<TSparseSpace,TDenseSpace>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( PoroNewmarkQuasistaticDampedUPlPgScheme );

    typedef Scheme<TSparseSpace,TDenseSpace>                      BaseType;
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
    using PoroNewmarkQuasistaticUPlPgScheme<TSparseSpace,TDenseSpace>::mDeltaTime;
    using PoroNewmarkQuasistaticUPlPgScheme<TSparseSpace,TDenseSpace>::mTheta_u;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ///Constructor
    PoroNewmarkQuasistaticDampedUPlPgScheme(double theta_u,
                                            double theta_pl,
                                            double theta_pg,
                                            double beta,
                                            double gamma) :
                                            PoroNewmarkQuasistaticUPlPgScheme<TSparseSpace,TDenseSpace>(
                                                theta_u,
                                                theta_pl,
                                                theta_pg,
                                                beta,
                                                gamma)
    {
        //Allocate auxiliary memory
        int NumThreads = ParallelUtilities::GetNumThreads();
        mDampingMatrix.resize(NumThreads);
        mVelocityVector.resize(NumThreads);
    }

    //------------------------------------------------------------------------------------

    ///Destructor
    ~PoroNewmarkQuasistaticDampedUPlPgScheme() override {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Note: this is in a parallel loop

    void CalculateSystemContributions(
        Element& rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        rCurrentElement.CalculateLocalSystem(LHS_Contribution,RHS_Contribution,rCurrentProcessInfo);

        rCurrentElement.CalculateDampingMatrix(mDampingMatrix[thread],rCurrentProcessInfo);

        this->AddDampingToLHS (LHS_Contribution, mDampingMatrix[thread], rCurrentProcessInfo);

        this->AddDampingToRHS (rCurrentElement, RHS_Contribution, mDampingMatrix[thread], rCurrentProcessInfo);

        rCurrentElement.EquationIdVector(EquationId,rCurrentProcessInfo);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Note: this is in a parallel loop

    void CalculateRHSContribution(
        Element& rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        rCurrentElement.CalculateRightHandSide(RHS_Contribution,rCurrentProcessInfo);

        rCurrentElement.CalculateDampingMatrix(mDampingMatrix[thread],rCurrentProcessInfo);

        this->AddDampingToRHS (rCurrentElement, RHS_Contribution, mDampingMatrix[thread], rCurrentProcessInfo);

        rCurrentElement.EquationIdVector(EquationId,rCurrentProcessInfo);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Note: this is in a parallel loop

    void CalculateLHSContribution(
        Element& rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        rCurrentElement.CalculateLeftHandSide(LHS_Contribution,rCurrentProcessInfo);

        rCurrentElement.CalculateDampingMatrix(mDampingMatrix[thread],rCurrentProcessInfo);

        this->AddDampingToLHS (LHS_Contribution, mDampingMatrix[thread], rCurrentProcessInfo);

        rCurrentElement.EquationIdVector(EquationId,rCurrentProcessInfo);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

    std::vector< Matrix > mDampingMatrix;
    std::vector< Vector > mVelocityVector;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void AddDampingToLHS(LocalSystemMatrixType& LHS_Contribution,LocalSystemMatrixType& C,const ProcessInfo& rCurrentProcessInfo)
    {
        // adding damping contribution
        if (C.size1() != 0)
            noalias(LHS_Contribution) += 1.0/(mTheta_u*mDeltaTime)*C;
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void AddDampingToRHS(Element& rCurrentElement,LocalSystemVectorType& RHS_Contribution,
                            LocalSystemMatrixType& C,const ProcessInfo& rCurrentProcessInfo)
    {
        int thread = OpenMPUtils::ThisThread();

        //adding damping contribution
        if (C.size1() != 0)
        {
            rCurrentElement.GetFirstDerivativesVector(mVelocityVector[thread], 0);

            noalias(RHS_Contribution) -= prod(C, mVelocityVector[thread]);
        }
    }

}; // Class PoroNewmarkQuasistaticDampedUPlPgScheme
}  // namespace Kratos

#endif // KRATOS_PORO_NEWMARK_QUASISTATIC_DAMPED_U_PL_PG_SCHEME defined
