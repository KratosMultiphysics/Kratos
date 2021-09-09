//
//   Project Name:        KratosDamApplication   $
//   Last Modified by:    $Author:Ignasi de Pouplana $
//   Date:                $Date:    February 2017$
//   Revision:            $Revision:         1.0 $
//

#if !defined(KRATOS_TRILINOS_INCREMENTAL_UPDATE_STATIC_DAMPED_SCHEME )
#define  KRATOS_TRILINOS_INCREMENTAL_UPDATE_STATIC_DAMPED_SCHEME

// Application includes
#include "custom_strategies/schemes/trilinos_residual_based_bossak_displacement_scheme.h"
#include "dam_application_variables.h"

namespace Kratos
{

template<class TSparseSpace, class TDenseSpace>

class TrilinosIncrementalUpdateStaticDampedScheme : public TrilinosResidualBasedBossakDisplacementScheme<TSparseSpace,TDenseSpace>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( TrilinosIncrementalUpdateStaticDampedScheme );

    typedef Scheme<TSparseSpace,TDenseSpace>     BaseType;
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
    typedef ResidualBasedBossakDisplacementScheme<TSparseSpace,TDenseSpace> GrandMotherType;
    using GrandMotherType::mMatrix;
    using GrandMotherType::mAlpha;
    using GrandMotherType::mNewmark;
    using GrandMotherType::mVector;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ///Constructor
    TrilinosIncrementalUpdateStaticDampedScheme(double rAlpham = 0.0)
        : TrilinosResidualBasedBossakDisplacementScheme<TSparseSpace,TDenseSpace>(rAlpham) {}

    //------------------------------------------------------------------------------------

    ///Destructor
    virtual ~TrilinosIncrementalUpdateStaticDampedScheme() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateSystemContributions(
        Element& rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY;

        int thread = OpenMPUtils::ThisThread();

        rCurrentElement.CalculateLocalSystem(LHS_Contribution, RHS_Contribution, CurrentProcessInfo);

        rCurrentElement.CalculateDampingMatrix(mMatrix.D[thread], CurrentProcessInfo);

        this->AddDampingToLHS (LHS_Contribution, mMatrix.D[thread], CurrentProcessInfo);

        this->AddDampingToRHS (rCurrentElement, RHS_Contribution, mMatrix.D[thread], CurrentProcessInfo);

        rCurrentElement.EquationIdVector(EquationId, CurrentProcessInfo);

        KRATOS_CATCH( "" );
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Calculate_RHS_Contribution(
        Element& rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo)
    {

        KRATOS_TRY;

        int thread = OpenMPUtils::ThisThread();

        rCurrentElement.CalculateRightHandSide(RHS_Contribution, CurrentProcessInfo);

        rCurrentElement.CalculateDampingMatrix(mMatrix.D[thread], CurrentProcessInfo);

        this->AddDampingToRHS (rCurrentElement, RHS_Contribution, mMatrix.D[thread], CurrentProcessInfo);

        rCurrentElement.EquationIdVector(EquationId, CurrentProcessInfo);

        KRATOS_CATCH( "" );
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateSystemContributions(
        Condition& rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        rCurrentCondition.CalculateLocalSystem(LHS_Contribution,RHS_Contribution, CurrentProcessInfo);

        rCurrentCondition.EquationIdVector(EquationId, CurrentProcessInfo);

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Calculate_RHS_Contribution(
        Condition& rCurrentCondition,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY

        rCurrentCondition.CalculateRightHandSide(RHS_Contribution, CurrentProcessInfo);

        rCurrentCondition.EquationIdVector(EquationId, CurrentProcessInfo);

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void AddDampingToLHS(
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemMatrixType& D,
        const ProcessInfo& CurrentProcessInfo)
    {
        // Adding  damping contribution
        if (D.size1() != 0) // if D matrix declared
        {
            noalias(LHS_Contribution) += D * (1.0 - mAlpha.f) * mNewmark.c1;
        }
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void AddDampingToRHS(
        Element& rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        LocalSystemMatrixType& D,
        const ProcessInfo& CurrentProcessInfo)
    {
        int thread = OpenMPUtils::ThisThread();

        // Adding damping contribution
        if (D.size1() != 0)
        {
            rCurrentElement.GetFirstDerivativesVector(mVector.v[thread], 0);

            noalias(RHS_Contribution) -= prod(D, mVector.v[thread]);
        }
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

}; // Class TrilinosIncrementalUpdateStaticDampedScheme
}  // namespace Kratos

#endif // KRATOS_TRILINOS_INCREMENTAL_UPDATE_STATIC_DAMPED_SCHEME defined