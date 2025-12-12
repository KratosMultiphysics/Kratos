// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

#pragma once

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_elements/U_Pw_small_strain_element.h"
#include "custom_utilities/element_utilities.hpp"
#include "custom_utilities/stress_strain_utilities.h"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) TransientPwElement : public UPwSmallStrainElement<TDim, TNumNodes>
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(TransientPwElement);

    using BaseType = UPwSmallStrainElement<TDim, TNumNodes>;

    using IndexType            = std::size_t;
    using PropertiesType       = Properties;
    using NodeType             = Node;
    using GeometryType         = Geometry<NodeType>;
    using NodesArrayType       = GeometryType::PointsArrayType;
    using VectorType           = Vector;
    using MatrixType           = Matrix;
    using DofsVectorType       = Element::DofsVectorType;
    using EquationIdVectorType = Element::EquationIdVectorType;

    /// The definition of the sizetype
    using SizeType = std::size_t;
    using BaseType::mConstitutiveLawVector;
    using BaseType::mRetentionLawVector;

    using ElementVariables = typename BaseType::ElementVariables;

    explicit TransientPwElement(IndexType NewId = 0) : BaseType(NewId) {}

    /// Constructor using an array of nodes
    TransientPwElement(IndexType                          NewId,
                       const NodesArrayType&              ThisNodes,
                       std::unique_ptr<StressStatePolicy> pStressStatePolicy,
                       std::unique_ptr<IntegrationCoefficientModifier> pCoefficientModifier = nullptr)
        : BaseType(NewId, ThisNodes, std::move(pStressStatePolicy), std::move(pCoefficientModifier))
    {
    }

    /// Constructor using Geometry
    TransientPwElement(IndexType                          NewId,
                       GeometryType::Pointer              pGeometry,
                       std::unique_ptr<StressStatePolicy> pStressStatePolicy,
                       std::unique_ptr<IntegrationCoefficientModifier> pCoefficientModifier = nullptr)
        : BaseType(NewId, pGeometry, std::move(pStressStatePolicy), std::move(pCoefficientModifier))
    {
    }

    /// Constructor using Properties
    TransientPwElement(IndexType                          NewId,
                       GeometryType::Pointer              pGeometry,
                       PropertiesType::Pointer            pProperties,
                       std::unique_ptr<StressStatePolicy> pStressStatePolicy,
                       std::unique_ptr<IntegrationCoefficientModifier> pCoefficientModifier = nullptr)
        : BaseType(NewId, pGeometry, pProperties, std::move(pStressStatePolicy), std::move(pCoefficientModifier))
    {
    }

    ~TransientPwElement() = default;

    Element::Pointer Create(IndexType               NewId,
                            NodesArrayType const&   ThisNodes,
                            PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    void InitializeNonLinearIteration(const ProcessInfo&) override;

    void FinalizeNonLinearIteration(const ProcessInfo&) override;

    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    void GetValuesVector(Vector& rValues, int Step = 0) const override;

    void GetFirstDerivativesVector(Vector& rValues, int Step = 0) const override;

    void GetSecondDerivativesVector(Vector& rValues, int Step = 0) const override;

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                                      std::vector<double>&    rOutput,
                                      const ProcessInfo&      rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3>>& rVariable,
                                      std::vector<array_1d<double, 3>>&    rOutput,
                                      const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                      std::vector<Matrix>&    rOutput,
                                      const ProcessInfo&      rCurrentProcessInfo) override;

    // Turn back information as a string.
    std::string Info() const override;

    // Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

protected:
    void CalculateAll(MatrixType&        rLeftHandSideMatrix,
                      VectorType&        rRightHandSideVector,
                      const ProcessInfo& CurrentProcessInfo,
                      bool               CalculateStiffnessMatrixFlag,
                      bool               CalculateResidualVectorFlag) override;

    void InitializeElementVariables(ElementVariables& rVariables, const ProcessInfo& CurrentProcessInfo) override;

    void CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables) override;

    void CalculateAndAddRHS(VectorType& rRightHandSideVector, ElementVariables& rVariables, unsigned int) override;

    void CalculateKinematics(ElementVariables& rVariables, unsigned int IntegrationPointIndex) override;

    void CalculateAndAddCompressibilityMatrix(MatrixType&             rLeftHandSideMatrix,
                                              const ElementVariables& rVariables) override;
    void CalculateAndAddPermeabilityFlow(VectorType&             rRightHandSideVector,
                                         const ElementVariables& rVariables) override;
    void CalculateAndAddFluidBodyFlow(VectorType& rRightHandSideVector, const ElementVariables& rVariables) override;
    void CalculateAndAddCompressibilityFlow(VectorType&             rRightHandSideVector,
                                            const ElementVariables& rVariables) override;

    std::size_t GetNumberOfDOF() const override;

private:
    [[nodiscard]] DofsVectorType GetDofs() const override;

    /// Serialization

    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;

    // Assignment operator.
    TransientPwElement& operator=(TransientPwElement const& rOther);

    // Copy constructor.
    TransientPwElement(TransientPwElement const& rOther);

}; // Class TransientPwElement

} // namespace Kratos
