// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Armin Geiser, https://github.com/armingeiser
//

// System includes
#if !defined(ADJOINT_SEMI_ANALYTIC_BASE_CONDITION )
#define  ADJOINT_SEMI_ANALYTIC_BASE_CONDITION

// System includes

// External includes

// Project includes
#include "includes/condition.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/** \brief AdjointSemiAnalyticBaseCondition
 *
 * This condition is a wrapper for a primal condition to calculate condition derivatives using
 * finite differencing (adjoint semi analytic approach). It is designed to be used in adjoint
 * sensitivity analysis
 */

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION)  AdjointSemiAnalyticBaseCondition
    : public Condition
{
public:
    ///@name Type Definitions
    ///@{
    /// Counted pointer of AdjointSemiAnalyticBaseCondition
    KRATOS_CLASS_POINTER_DEFINITION( AdjointSemiAnalyticBaseCondition );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    AdjointSemiAnalyticBaseCondition(): Condition()
    {
    }

    AdjointSemiAnalyticBaseCondition(
        Condition::Pointer pPrimalCondition
        ): Condition( pPrimalCondition->Id(), pPrimalCondition->pGetGeometry(), pPrimalCondition->pGetProperties() ),
           mpPrimalCondition(pPrimalCondition)
        {
        }

    /// Destructor.
    ~AdjointSemiAnalyticBaseCondition() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Informations
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo ) override
    {
        KRATOS_ERROR << "EquationIdVector of the base class called!" << std::endl;
    }

    void GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& rCurrentProcessInfo ) override
    {
        KRATOS_ERROR << "GetDofList of the base class called!" << std::endl;
    }

    IntegrationMethod GetIntegrationMethod() override
    {
        return mpPrimalCondition->GetIntegrationMethod();
    }

    void GetValuesVector(Vector& rValues, int Step = 0 ) override
    {
        KRATOS_ERROR << "GetValuesVector of the base class called!" << std::endl;
    }

    void Initialize() override
    {
        mpPrimalCondition->Initialize();
    }

    void ResetConstitutiveLaw() override
    {
        mpPrimalCondition->ResetConstitutiveLaw();
    }

    void CleanMemory() override
    {
        mpPrimalCondition->CleanMemory();
    }

    void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->InitializeSolutionStep(rCurrentProcessInfo);
    }

    void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->InitializeNonLinearIteration(rCurrentProcessInfo);
    }

    void FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->FinalizeNonLinearIteration(rCurrentProcessInfo);
    }

    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->FinalizeSolutionStep(rCurrentProcessInfo);
    }

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
				      VectorType& rRightHandSideVector,
				      ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->CalculateLocalSystem(rLeftHandSideMatrix,
				    rRightHandSideVector,
				    rCurrentProcessInfo);
    }

    void CalculateLocalSystem(std::vector< MatrixType >& rLeftHandSideMatrices,
                                      const std::vector< Variable< MatrixType > >& rLHSVariables,
                                      std::vector< VectorType >& rRightHandSideVectors,
                                      const std::vector< Variable< VectorType > >& rRHSVariables,
                                      ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->CalculateLocalSystem(rLeftHandSideMatrices,
                                      rLHSVariables,
                                      rRightHandSideVectors,
                                      rRHSVariables,
                                      rCurrentProcessInfo);
    }

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
				       ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->CalculateLeftHandSide(rLeftHandSideMatrix,
					rCurrentProcessInfo);
    }

    void CalculateLeftHandSide(std::vector< MatrixType >& rLeftHandSideMatrices,
					const std::vector< Variable< MatrixType > >& rLHSVariables,
					ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->CalculateLeftHandSide(rLeftHandSideMatrices,
					rLHSVariables,
					rCurrentProcessInfo);
    }

    void CalculateRightHandSide(VectorType& rRightHandSideVector,
					ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->CalculateRightHandSide(rRightHandSideVector,
					rCurrentProcessInfo);
    }

    void CalculateRightHandSide(std::vector< VectorType >& rRightHandSideVectors,
					const std::vector< Variable< VectorType > >& rRHSVariables,
					ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->CalculateRightHandSide(rRightHandSideVectors,
					rRHSVariables,
					rCurrentProcessInfo);
    }

    void CalculateFirstDerivativesContributions(MatrixType& rLeftHandSideMatrix,
							VectorType& rRightHandSideVector,
							ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->CalculateFirstDerivativesContributions(rLeftHandSideMatrix,
							rRightHandSideVector,
							rCurrentProcessInfo);
    }

    void CalculateFirstDerivativesLHS(MatrixType& rLeftHandSideMatrix,
					      ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->CalculateFirstDerivativesLHS(rLeftHandSideMatrix,
					    rCurrentProcessInfo);
    }

    void CalculateFirstDerivativesRHS(VectorType& rRightHandSideVector,
					      ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->CalculateFirstDerivativesRHS(rRightHandSideVector,
					      rCurrentProcessInfo);
    }

    void CalculateSecondDerivativesContributions(MatrixType& rLeftHandSideMatrix,
							 VectorType& rRightHandSideVector,
							 ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->CalculateSecondDerivativesContributions(rLeftHandSideMatrix,
							 rRightHandSideVector,
							 rCurrentProcessInfo);
    }

    void CalculateSecondDerivativesLHS(MatrixType& rLeftHandSideMatrix,
					       ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->CalculateSecondDerivativesLHS(rLeftHandSideMatrix,
					       rCurrentProcessInfo);
    }

    void CalculateSecondDerivativesRHS(VectorType& rRightHandSideVector,
					       ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->CalculateSecondDerivativesRHS(rRightHandSideVector,
					       rCurrentProcessInfo);
    }

    void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->CalculateMassMatrix(rMassMatrix, rCurrentProcessInfo);
    }

    void CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->CalculateDampingMatrix(rDampingMatrix, rCurrentProcessInfo);
    }

    void Calculate(const Variable<double >& rVariable,
			   double& Output,
			   const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_ERROR << "Calculate of the adjoint base condition is called!" << std::endl;
    }

    void Calculate(const Variable< array_1d<double,3> >& rVariable,
			   array_1d<double,3>& Output,
			   const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_ERROR << "Calculate of the adjoint base condition is called!" << std::endl;
    }

    void Calculate(const Variable<Vector >& rVariable,
			   Vector& Output,
			   const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_ERROR << "Calculate of the adjoint base condition is called!" << std::endl;
    }

    void Calculate(const Variable<Matrix >& rVariable,
			   Matrix& Output,
			   const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_ERROR << "Calculate of the adjoint base condition is called!" << std::endl;
    }


    void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
					      std::vector<double>& rOutput,
					      const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_ERROR << "CalculateOnIntegrationPoints of the adjoint base condition is called!" << std::endl;
    }

    void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
					      std::vector< array_1d<double, 3 > >& Output,
					      const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_ERROR << "CalculateOnIntegrationPoints of the adjoint base condition is called!" << std::endl;
    }

    void CalculateOnIntegrationPoints(const Variable<Vector >& rVariable,
					      std::vector< Vector >& Output,
					      const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_ERROR << "CalculateOnIntegrationPoints of the adjoint base condition is called!" << std::endl;
    }

    void CalculateOnIntegrationPoints(const Variable<Matrix >& rVariable,
					      std::vector< Matrix >& Output,
					      const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_ERROR << "CalculateOnIntegrationPoints of the adjoint base condition is called!" << std::endl;
    }

    int Check( const ProcessInfo& rCurrentProcessInfo ) override
    {
        KRATOS_ERROR << "Check of the base class called!" << std::endl;
    }

    /**
     * Calculates the pseudo-load contribution of the condition w.r.t.  a scalar design variable.
     */
    void CalculateSensitivityMatrix(const Variable<double>& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_ERROR << "CalculateSensitivityMatrix of the base class called!" << std::endl;
    }

    /**
     * Calculates the pseudo-load contribution of the condition w.r.t.  a vector design variable.
     */
    void CalculateSensitivityMatrix(const Variable<array_1d<double,3> >& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_ERROR << "CalculateSensitivityMatrix of the base class called!" << std::endl;
    }

    Condition::Pointer pGetPrimalCondition()
    {
        return mpPrimalCondition;
    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{

    /**
     * pointer to the primal condition
     */
    Condition::Pointer mpPrimalCondition;

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{



    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition );
        rSerializer.save("mpPrimalCondition", mpPrimalCondition);
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition );
        rSerializer.load("mpPrimalCondition", mpPrimalCondition);
    }

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class AdjointSemiAnalyticBaseCondition

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.

#endif // ADJOINT_SEMI_ANALYTIC_BASE_CONDITION  defined


