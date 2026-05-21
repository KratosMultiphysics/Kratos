// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Armin Geiser, https://github.com/armingeiser
//

// System includes
#pragma once

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
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( AdjointSemiAnalyticBaseCondition );

    ///@}
    ///@name Life Cycle
    ///@{

    AdjointSemiAnalyticBaseCondition(IndexType NewId = 0)
    : Condition(NewId),
      mpPrimalCondition(Kratos::make_intrusive<TPrimalCondition>(NewId, pGetGeometry()))
    {
    }

    AdjointSemiAnalyticBaseCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : Condition(NewId, pGeometry),
      mpPrimalCondition(Kratos::make_intrusive<TPrimalCondition>(NewId, pGeometry))
    {
    }

    AdjointSemiAnalyticBaseCondition(IndexType NewId,
                        GeometryType::Pointer pGeometry,
                        PropertiesType::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties),
      mpPrimalCondition(Kratos::make_intrusive<TPrimalCondition>(NewId, pGeometry, pProperties))
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

    Condition::Pointer Create(IndexType NewId,
                              NodesArrayType const& ThisNodes,
                              PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<AdjointSemiAnalyticBaseCondition<TPrimalCondition>>(
            NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    Condition::Pointer Create(IndexType NewId,
                              GeometryType::Pointer pGeometry,
                              PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<AdjointSemiAnalyticBaseCondition<TPrimalCondition>>(
            NewId, pGeometry, pProperties);
    }

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo ) const override;

    void GetDofList(DofsVectorType& ElementalDofList, const ProcessInfo& rCurrentProcessInfo ) const override;

    IntegrationMethod GetIntegrationMethod() const override
    {
        return mpPrimalCondition->GetIntegrationMethod();
    }

    void GetValuesVector(Vector& rValues, int Step = 0 ) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->Initialize(rCurrentProcessInfo);
    }

    void ResetConstitutiveLaw() override
    {
        mpPrimalCondition->ResetConstitutiveLaw();
    }

    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->InitializeSolutionStep(rCurrentProcessInfo);
    }

    void InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->InitializeNonLinearIteration(rCurrentProcessInfo);
    }

    void FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->FinalizeNonLinearIteration(rCurrentProcessInfo);
    }

    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->FinalizeSolutionStep(rCurrentProcessInfo);
    }

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
				      VectorType& rRightHandSideVector,
				      const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->CalculateLocalSystem(rLeftHandSideMatrix,
				    rRightHandSideVector,
				    rCurrentProcessInfo);
    }

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
				       const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->CalculateLeftHandSide(rLeftHandSideMatrix,
					rCurrentProcessInfo);
    }

    void CalculateRightHandSide(VectorType& rRightHandSideVector,
					const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->CalculateRightHandSide(rRightHandSideVector,
					rCurrentProcessInfo);
    }

    void CalculateFirstDerivativesContributions(MatrixType& rLeftHandSideMatrix,
							VectorType& rRightHandSideVector,
							const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->CalculateFirstDerivativesContributions(rLeftHandSideMatrix,
							rRightHandSideVector,
							rCurrentProcessInfo);
    }

    void CalculateFirstDerivativesLHS(MatrixType& rLeftHandSideMatrix,
					      const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->CalculateFirstDerivativesLHS(rLeftHandSideMatrix,
					    rCurrentProcessInfo);
    }

    void CalculateFirstDerivativesRHS(VectorType& rRightHandSideVector,
					      const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->CalculateFirstDerivativesRHS(rRightHandSideVector,
					      rCurrentProcessInfo);
    }

    void CalculateSecondDerivativesContributions(MatrixType& rLeftHandSideMatrix,
							 VectorType& rRightHandSideVector,
							 const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->CalculateSecondDerivativesContributions(rLeftHandSideMatrix,
							 rRightHandSideVector,
							 rCurrentProcessInfo);
    }

    void CalculateSecondDerivativesLHS(MatrixType& rLeftHandSideMatrix,
					       const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->CalculateSecondDerivativesLHS(rLeftHandSideMatrix,
					       rCurrentProcessInfo);
    }

    void CalculateSecondDerivativesRHS(VectorType& rRightHandSideVector,
					       const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->CalculateSecondDerivativesRHS(rRightHandSideVector,
					       rCurrentProcessInfo);
    }

    void CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->CalculateMassMatrix(rMassMatrix, rCurrentProcessInfo);
    }

    void CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo) override
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

    // Results calculation on integration points
    void CalculateOnIntegrationPoints(
        const Variable<bool>& rVariable,
        std::vector<bool>& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(
        const Variable<array_1d<double, 3>>& rVariable,
        std::vector<array_1d<double, 3>>& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(
        const Variable<array_1d<double, 4>>& rVariable,
        std::vector<array_1d<double, 4>>& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(
        const Variable<array_1d<double, 6>>& rVariable,
        std::vector<array_1d<double, 6>>& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(
        const Variable<array_1d<double, 9>>& rVariable,
        std::vector<array_1d<double, 9>>& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(
        const Variable<Matrix>& rVariable,
        std::vector<Matrix>& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;

    int Check( const ProcessInfo& rCurrentProcessInfo ) const override;

    /**
     * Calculates the pseudo-load contribution of the condition w.r.t.  a scalar design variable.
     */
    void CalculateSensitivityMatrix(const Variable<double>& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Calculates the pseudo-load contribution of the condition w.r.t.  a vector design variable.
     */
    void CalculateSensitivityMatrix(const Variable<array_1d<double,3> >& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) override;

    Condition::Pointer pGetPrimalCondition()
    {
        return mpPrimalCondition;
    }

    const Condition::Pointer pGetPrimalCondition() const
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

    /**
     * Get the perturbation size for a scalar variable
     */
    double GetPerturbationSize(const Variable<double>& rDesignVariable, const ProcessInfo& rCurrentProcessInfo) const;

    /**
     * Get the perturbation size for a vector variable
     */
    double GetPerturbationSize(const Variable<array_1d<double,3>>& rDesignVariable, const ProcessInfo& rCurrentProcessInfo) const;

    /**
     * Get the perturbation size modification factor for a scalar variable.
     * The computed factor reflects the current value of the property (design variable).
     * Note: This approach is only based on experience.
     * This can be overwritten by derived classes.
     */
    virtual double GetPerturbationSizeModificationFactor(const Variable<double>& rVariable) const;

    /**
     * Get the perturbation size modification factor for a vector variable.
     * The computed factor reflects the size of the geometry.
     * Note: This approach is only based on experience.
     * This can be overwritten by derived classes.
     */
    virtual double GetPerturbationSizeModificationFactor(const Variable<array_1d<double,3>>& rDesignVariable) const;


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


