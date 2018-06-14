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
 * This consition a wrapper for a primal condition to calculate condition derivatives using
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
    AdjointSemiAnalyticBaseCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry
        ): Condition( NewId, pGeometry )
    {
    }

    AdjointSemiAnalyticBaseCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties,
        Condition::Pointer pPrimalCondition
        ): Condition( NewId, pGeometry )
        {
            mpPrimalCondition = pPrimalCondition;
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

    // not a virtual function in condition.h
    // SizeType WorkingSpaceDimension() const override
    // {
    //      return mpPrimalCondition->WorkingSpaceDimension();
    // }

    ///@}
    ///@name Operations
    ///@{

    virtual Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties, Condition::Pointer pPrimalCondition ) const
    {
        return Kratos::make_shared<AdjointSemiAnalyticBaseCondition>( NewId, GetGeometry().Create( ThisNodes ), pProperties, pPrimalCondition );
    }

    // TODO add missing create and clone methods

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

    void GetFirstDerivativesVector(Vector& values, int Step = 0) override
    {
        mpPrimalCondition->GetFirstDerivativesVector(values, Step);
    }

    void GetSecondDerivativesVector(Vector& values, int Step = 0) override
    {
        mpPrimalCondition->GetSecondDerivativesVector(values, Step);
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

    // TODO explicit functions lines 668-698

    void Calculate(const Variable<double >& rVariable,
			   double& Output,
			   const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->Calculate(rVariable,
			   Output,
			   rCurrentProcessInfo);
    }

    void Calculate(const Variable< array_1d<double,3> >& rVariable,
			   array_1d<double,3>& Output,
			   const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->Calculate(rVariable,
			   Output,
			   rCurrentProcessInfo);
    }

    void Calculate(const Variable<Vector >& rVariable,
			   Vector& Output,
			   const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->Calculate(rVariable,
			   Output,
			   rCurrentProcessInfo);
    }

    void Calculate(const Variable<Matrix >& rVariable,
			   Matrix& Output,
			   const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->Calculate(rVariable,
			   Output,
			   rCurrentProcessInfo);
    }


    void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
					      std::vector<double>& rOutput,
					      const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->CalculateOnIntegrationPoints(rVariable,
					      rOutput,
					      rCurrentProcessInfo);
    }

    void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
					      std::vector< array_1d<double, 3 > >& Output,
					      const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->CalculateOnIntegrationPoints(rVariable,
					      Output,
					      rCurrentProcessInfo);
    }

    void CalculateOnIntegrationPoints(const Variable<Vector >& rVariable,
					      std::vector< Vector >& Output,
					      const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->CalculateOnIntegrationPoints(rVariable,
					      Output,
					      rCurrentProcessInfo);
    }

    void CalculateOnIntegrationPoints(const Variable<Matrix >& rVariable,
					      std::vector< Matrix >& Output,
					      const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->CalculateOnIntegrationPoints(rVariable,
					      Output,
					      rCurrentProcessInfo);
    }

    void SetValueOnIntegrationPoints(const Variable<double>& rVariable,
					     std::vector<double>& rValues,
					     const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->SetValueOnIntegrationPoints(rVariable,
					     rValues,
					     rCurrentProcessInfo);
    }

    void SetValueOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
					     std::vector<array_1d<double, 3 > > rValues,
					     const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->SetValueOnIntegrationPoints(rVariable,
					     rValues,
					     rCurrentProcessInfo);
    }

    void SetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
					     std::vector<Vector>& rValues,
					     const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->SetValueOnIntegrationPoints(rVariable,
					     rValues,
					     rCurrentProcessInfo);
    }

    void SetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
					     std::vector<Matrix>& rValues,
					     const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->SetValueOnIntegrationPoints(rVariable,
					     rValues,
					     rCurrentProcessInfo);
    }


    void GetValueOnIntegrationPoints(const Variable<double>& rVariable,
					     std::vector<double>& rValues,
					     const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->GetValueOnIntegrationPoints(rVariable,
					     rValues,
					     rCurrentProcessInfo);
    }

    void GetValueOnIntegrationPoints(const Variable<int>& rVariable,
					     std::vector<int>& rValues,
					     const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->GetValueOnIntegrationPoints(rVariable,
					     rValues,
					     rCurrentProcessInfo);
    }

    void GetValueOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
					     std::vector<array_1d<double, 3 > >& rValues,
					     const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->GetValueOnIntegrationPoints(rVariable,
					     rValues,
					     rCurrentProcessInfo);
    }

    void GetValueOnIntegrationPoints(const Variable<array_1d<double, 6 > >& rVariable,
					     std::vector<array_1d<double, 6 > >& rValues,
					     const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->GetValueOnIntegrationPoints(rVariable,
					     rValues,
					     rCurrentProcessInfo);
    }

    void GetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
					     std::vector<Vector>& rValues,
					     const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->GetValueOnIntegrationPoints(rVariable,
					     rValues,
					     rCurrentProcessInfo);
    }

    void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
					     std::vector<Matrix>& rValues,
					     const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalCondition->GetValueOnIntegrationPoints(rVariable,
					     rValues,
					     rCurrentProcessInfo);
    }

    int Check( const ProcessInfo& rCurrentProcessInfo ) override
    {
        KRATOS_ERROR << "Check of the base class called!" << std::endl;
    }

    void CalculateSensitivityMatrix(const Variable<double>& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_ERROR << "CalculateSensitivityMatrix of the base class called!" << std::endl;
    }

    void CalculateSensitivityMatrix(const Variable<array_1d<double,3> >& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_ERROR << "CalculateSensitivityMatrix of the base class called!" << std::endl;
    }


    // PropertiesType::Pointer pGetProperties() override
    // {
    //     return mpPrimalCondition->pGetProperties();
    // }

    // const PropertiesType::Pointer pGetProperties() const override
    // {
    //     return mpPrimalCondition->pGetProperties();
    // }

    // PropertiesType& GetProperties() override
    // {
    //     return mpPrimalCondition->GetProperties();
    // }

    // PropertiesType const& GetProperties() const override
    // {
    //     return mpPrimalCondition->GetProperties();
    // }

    // void SetProperties(PropertiesType::Pointer pProperties) override
    // {
    //     mpPrimalCondition->SetProperties(pProperties);
    // }


    // DataValueContainer& Data() override
    // {
    //     return mpPrimalCondition->Data();
    // }

    // DataValueContainer const& GetData() const override
    // {
    //     return mpPrimalCondition->GetData();
    // }

    // void SetData(DataValueContainer const& rThisData) override
    // {
    //     mpPrimalCondition->SetData(rThisData);
    // }

    // not a virtual function in condition.h
    // template<class TDataType> bool Has(const Variable<TDataType>& rThisVariable) const override
    // {
    //     return mpPrimalCondition->Has(rThisVariable);
    // }

    // template<class TAdaptorType> bool Has(
    //     const VariableComponent<TAdaptorType>& rThisVariable) const override
    // {
    //     return mpPrimalCondition->Has(rThisVariable);
    // }

    // template<class TVariableType> void SetValue(
    //     const TVariableType& rThisVariable,
    //     typename TVariableType::Type const& rValue) override
    // {
    //     return mpPrimalCondition->SetValue(
    //         rThisVariable,
    //         rValue);
    // }

    // template<class TVariableType> typename TVariableType::Type& GetValue(
    //     const TVariableType& rThisVariable) override
    // {
    //     return mpPrimalCondition->GetValue(rThisVariable);
    // }

    // template<class TVariableType> typename TVariableType::Type const& GetValue(
    //     const TVariableType& rThisVariable) const override
    // {
    //     return mpPrimalCondition->GetValue(rThisVariable);
    // }

    // Flags& GetFlags() override
    // {
    //     return mpPrimalCondition->GetFlags();
    // }

    // Flags const& GetFlags() const override
    // {
    //     return mpPrimalCondition->GetFlags();
    // }

    // void SetFlags(Flags const& rThisFlags) override
    // {
    //     mpPrimalCondition->SetFlags(rThisFlags);
    // }

    // bool HasProperties() const override
    // {
    //     return mpPrimalCondition->HasProperties();
    // }

    // std::string Info() const override
    // {
    //     return mpPrimalCondition->Info();
    // }

    // void PrintInfo(std::ostream& rOStream) const override
    // {
    //     mpPrimalCondition->PrintInfo(rOStream);
    // }

    // void PrintData(std::ostream& rOStream) const override
    // {
    //     mpPrimalCondition->PrintData(rOStream);
    // }

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

    // A protected default constructor necessary for serialization
    AdjointSemiAnalyticBaseCondition(): Condition(){};

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
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition );
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


