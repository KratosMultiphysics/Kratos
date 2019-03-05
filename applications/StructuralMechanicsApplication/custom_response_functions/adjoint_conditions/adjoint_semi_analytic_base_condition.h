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

template <typename TPrimalCondition>
class AdjointSemiAnalyticBaseCondition
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

    AdjointSemiAnalyticBaseCondition(IndexType NewId = 0)
    : Condition(NewId), mPrimalCondition(NewId, pGetGeometry())
    {
    }

    AdjointSemiAnalyticBaseCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : Condition(NewId, pGeometry), mPrimalCondition(NewId, pGeometry)
    {
    }

    AdjointSemiAnalyticBaseCondition(IndexType NewId,
                        GeometryType::Pointer pGeometry,
                        PropertiesType::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties), mPrimalCondition(NewId, pGeometry, pProperties)
    {
    }

    virtual ~AdjointSemiAnalyticBaseCondition() {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Informations
    ///@{

    ///@}
    ///@name Operations
    ///@{

    virtual Condition::Pointer Create(IndexType NewId,
                              NodesArrayType const& ThisNodes,
                              PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_shared<AdjointSemiAnalyticBaseCondition<TPrimalCondition>>(
            NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    virtual Condition::Pointer Create(IndexType NewId,
                              GeometryType::Pointer pGeometry,
                              PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_shared<AdjointSemiAnalyticBaseCondition<TPrimalCondition>>(
            NewId, pGeometry, pProperties);
    }

    virtual void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo ) override
    {
        KRATOS_ERROR << "EquationIdVector of the base class called!" << std::endl;
    }

    virtual void GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& rCurrentProcessInfo ) override
    {
        KRATOS_ERROR << "GetDofList of the base class called!" << std::endl;
    }

    IntegrationMethod GetIntegrationMethod() override
    {
        return mPrimalCondition.GetIntegrationMethod();
    }

    virtual void GetValuesVector(Vector& rValues, int Step = 0 ) override
    {
        KRATOS_ERROR << "GetValuesVector of the base class called!" << std::endl;
    }

    void Initialize() override
    {
        mPrimalCondition.Initialize();
    }

    void ResetConstitutiveLaw() override
    {
        mPrimalCondition.ResetConstitutiveLaw();
    }

    void CleanMemory() override
    {
        mPrimalCondition.CleanMemory();
    }

    void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) override
    {
        mPrimalCondition.InitializeSolutionStep(rCurrentProcessInfo);
    }

    void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override
    {
        mPrimalCondition.InitializeNonLinearIteration(rCurrentProcessInfo);
    }

    void FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override
    {
        mPrimalCondition.FinalizeNonLinearIteration(rCurrentProcessInfo);
    }

    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override
    {
        mPrimalCondition.FinalizeSolutionStep(rCurrentProcessInfo);
    }

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
				      VectorType& rRightHandSideVector,
				      ProcessInfo& rCurrentProcessInfo) override
    {
        mPrimalCondition.CalculateLocalSystem(rLeftHandSideMatrix,
				    rRightHandSideVector,
				    rCurrentProcessInfo);
    }

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
				       ProcessInfo& rCurrentProcessInfo) override
    {
        mPrimalCondition.CalculateLeftHandSide(rLeftHandSideMatrix,
					rCurrentProcessInfo);
    }

    void CalculateLeftHandSide(std::vector< MatrixType >& rLeftHandSideMatrices,
					const std::vector< Variable< MatrixType > >& rLHSVariables,
					ProcessInfo& rCurrentProcessInfo) override
    {
        mPrimalCondition.CalculateLeftHandSide(rLeftHandSideMatrices,
					rLHSVariables,
					rCurrentProcessInfo);
    }

    void CalculateRightHandSide(VectorType& rRightHandSideVector,
					ProcessInfo& rCurrentProcessInfo) override
    {
        mPrimalCondition.CalculateRightHandSide(rRightHandSideVector,
					rCurrentProcessInfo);
    }

    void CalculateFirstDerivativesContributions(MatrixType& rLeftHandSideMatrix,
							VectorType& rRightHandSideVector,
							ProcessInfo& rCurrentProcessInfo) override
    {
        mPrimalCondition.CalculateFirstDerivativesContributions(rLeftHandSideMatrix,
							rRightHandSideVector,
							rCurrentProcessInfo);
    }

    void CalculateFirstDerivativesLHS(MatrixType& rLeftHandSideMatrix,
					      ProcessInfo& rCurrentProcessInfo) override
    {
        mPrimalCondition.CalculateFirstDerivativesLHS(rLeftHandSideMatrix,
					    rCurrentProcessInfo);
    }

    void CalculateFirstDerivativesRHS(VectorType& rRightHandSideVector,
					      ProcessInfo& rCurrentProcessInfo) override
    {
        mPrimalCondition.CalculateFirstDerivativesRHS(rRightHandSideVector,
					      rCurrentProcessInfo);
    }

    void CalculateSecondDerivativesContributions(MatrixType& rLeftHandSideMatrix,
							 VectorType& rRightHandSideVector,
							 ProcessInfo& rCurrentProcessInfo) override
    {
        mPrimalCondition.CalculateSecondDerivativesContributions(rLeftHandSideMatrix,
							 rRightHandSideVector,
							 rCurrentProcessInfo);
    }

    void CalculateSecondDerivativesLHS(MatrixType& rLeftHandSideMatrix,
					       ProcessInfo& rCurrentProcessInfo) override
    {
        mPrimalCondition.CalculateSecondDerivativesLHS(rLeftHandSideMatrix,
					       rCurrentProcessInfo);
    }

    void CalculateSecondDerivativesRHS(VectorType& rRightHandSideVector,
					       ProcessInfo& rCurrentProcessInfo) override
    {
        mPrimalCondition.CalculateSecondDerivativesRHS(rRightHandSideVector,
					       rCurrentProcessInfo);
    }

    void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) override
    {
        mPrimalCondition.CalculateMassMatrix(rMassMatrix, rCurrentProcessInfo);
    }

    void CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo) override
    {
        mPrimalCondition.CalculateDampingMatrix(rDampingMatrix, rCurrentProcessInfo);
    }

    virtual void Calculate(const Variable<double >& rVariable,
			   double& Output,
			   const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_ERROR << "Calculate of the adjoint base condition is called!" << std::endl;
    }

    virtual void Calculate(const Variable< array_1d<double,3> >& rVariable,
			   array_1d<double,3>& Output,
			   const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_ERROR << "Calculate of the adjoint base condition is called!" << std::endl;
    }

    virtual void Calculate(const Variable<Vector >& rVariable,
			   Vector& Output,
			   const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_ERROR << "Calculate of the adjoint base condition is called!" << std::endl;
    }

    virtual void Calculate(const Variable<Matrix >& rVariable,
			   Matrix& Output,
			   const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_ERROR << "Calculate of the adjoint base condition is called!" << std::endl;
    }


    virtual void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
					      std::vector<double>& rOutput,
					      const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_ERROR << "CalculateOnIntegrationPoints of the adjoint base condition is called!" << std::endl;
    }

    virtual void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
					      std::vector< array_1d<double, 3 > >& Output,
					      const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_ERROR << "CalculateOnIntegrationPoints of the adjoint base condition is called!" << std::endl;
    }

    virtual void CalculateOnIntegrationPoints(const Variable<Vector >& rVariable,
					      std::vector< Vector >& Output,
					      const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_ERROR << "CalculateOnIntegrationPoints of the adjoint base condition is called!" << std::endl;
    }

    virtual void CalculateOnIntegrationPoints(const Variable<Matrix >& rVariable,
					      std::vector< Matrix >& Output,
					      const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_ERROR << "CalculateOnIntegrationPoints of the adjoint base condition is called!" << std::endl;
    }

    virtual int Check( const ProcessInfo& rCurrentProcessInfo ) override
    {
        KRATOS_ERROR << "Check of the base class called!" << std::endl;
    }

    /**
     * Calculates the pseudo-load contribution of the condition w.r.t.  a scalar design variable.
     */
    virtual void CalculateSensitivityMatrix(const Variable<double>& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_ERROR << "CalculateSensitivityMatrix of the base class called!" << std::endl;
    }

    /**
     * Calculates the pseudo-load contribution of the condition w.r.t.  a vector design variable.
     */
    virtual void CalculateSensitivityMatrix(const Variable<array_1d<double,3> >& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_ERROR << "CalculateSensitivityMatrix of the base class called!" << std::endl;
    }

    Condition::Pointer pGetPrimalCondition()
    {
        return mPrimalCondition;
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

    TPrimalCondition mPrimalCondition;

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
        rSerializer.save("mPrimalCondition", mPrimalCondition);
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition );
        rSerializer.load("mPrimalCondition", mPrimalCondition);
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


