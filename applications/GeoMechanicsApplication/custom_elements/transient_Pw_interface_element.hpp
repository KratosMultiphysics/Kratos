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

#if !defined(KRATOS_GEO_TRANSIENT_PW_INTERFACE_ELEMENT_H_INCLUDED )
#define  KRATOS_GEO_TRANSIENT_PW_INTERFACE_ELEMENT_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_elements/U_Pw_small_strain_interface_element.hpp"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(GEO_MECHANICS_APPLICATION) TransientPwInterfaceElement :
    public UPwSmallStrainInterfaceElement<TDim,TNumNodes>
{

public:

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( TransientPwInterfaceElement );

    using BaseType = UPwSmallStrainInterfaceElement<TDim, TNumNodes>;

    using IndexType = std::size_t;
    using PropertiesType = Properties;
    using NodeType = Node;
    using GeometryType = Geometry<NodeType>;
    using NodesArrayType = GeometryType::PointsArrayType;
    using VectorType = Vector;
    using MatrixType = Matrix;

    using DofsVectorType = Element::DofsVectorType;
    using EquationIdVectorType = Element::EquationIdVectorType;

    /// The definition of the sizetype
    using SizeType = std::size_t;

    using BaseType::mRetentionLawVector;
    using BaseType::mThisIntegrationMethod;
    using BaseType::CalculateRetentionResponse;

    using InterfaceElementVariables = typename BaseType::InterfaceElementVariables;
    using SFGradAuxVariables = typename BaseType::SFGradAuxVariables;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    TransientPwInterfaceElement(IndexType NewId = 0) : UPwSmallStrainInterfaceElement<TDim,TNumNodes>( NewId ) {}

    /// Constructor using an array of nodes
    TransientPwInterfaceElement(IndexType NewId,
                                const NodesArrayType& ThisNodes) : UPwSmallStrainInterfaceElement<TDim,TNumNodes>(NewId, ThisNodes) {}

    /// Constructor using Geometry
    TransientPwInterfaceElement(IndexType NewId,
                                GeometryType::Pointer pGeometry) : UPwSmallStrainInterfaceElement<TDim,TNumNodes>(NewId, pGeometry) {}

    /// Constructor using Properties
    TransientPwInterfaceElement(IndexType NewId,
                                GeometryType::Pointer pGeometry,
                                PropertiesType::Pointer pProperties)
                                : UPwSmallStrainInterfaceElement<TDim,TNumNodes>( NewId, pGeometry, pProperties )
    {}

    /// Destructor
    ~TransientPwInterfaceElement() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Element::Pointer Create(IndexType NewId,
                            NodesArrayType const& ThisNodes,
                            PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId,
                            GeometryType::Pointer pGeom,
                            PropertiesType::Pointer pProperties) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    void GetDofList(DofsVectorType& rElementalDofList, 
                    const ProcessInfo& rCurrentProcessInfo) const override;

    void EquationIdVector( EquationIdVectorType& rResult,
                           const ProcessInfo& rCurrentProcessInfo ) const override;

    void GetValuesVector(Vector& rValues, int Step = 0) const override;

    void GetFirstDerivativesVector(Vector& rValues, int Step = 0) const override;

    void GetSecondDerivativesVector(Vector& rValues, int Step = 0) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateMassMatrix(MatrixType& rMassMatrix,
                             const ProcessInfo& rCurrentProcessInfo) override;

    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;
    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    void CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                      std::vector<Matrix>& rValues,
                                      const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                                     std::vector<double>& rValues,
                                     const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<array_1d<double,3>>& rVariable,
                                      std::vector<array_1d<double,3>>& rValues,
                                      const ProcessInfo& rCurrentProcessInfo) override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateOnLobattoIntegrationPoints(const Variable<array_1d<double,3>>& rVariable,
                                            std::vector<array_1d<double,3>>& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnLobattoIntegrationPoints(const Variable<Matrix>& rVariable,
                                             std::vector<Matrix>& rOutput,
                                             const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateAll( MatrixType& rLeftHandSideMatrix,
                       VectorType& rRightHandSideVector,
                       const ProcessInfo& CurrentProcessInfo,
                       const bool CalculateStiffnessMatrixFlag,
                       const bool CalculateResidualVectorFlag) override;

    void InitializeElementVariables(InterfaceElementVariables& rVariables,
                                    const GeometryType& Geom,
                                    const PropertiesType& Prop,
                                    const ProcessInfo& CurrentProcessInfo) override;

    void CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix,
                            InterfaceElementVariables& rVariables) override;

    void CalculateAndAddCompressibilityMatrix(MatrixType& rLeftHandSideMatrix,
                                              InterfaceElementVariables& rVariables) override;

    void CalculateAndAddPermeabilityMatrix(MatrixType& rLeftHandSideMatrix,
                                           InterfaceElementVariables& rVariables) override;

    void CalculateAndAddRHS(VectorType& rRightHandSideVector,
                            InterfaceElementVariables& rVariables,
                            unsigned int GPoint) override;

    void CalculateAndAddCompressibilityFlow(VectorType& rRightHandSideVector,
                                            InterfaceElementVariables& rVariables) override;

    void CalculateAndAddPermeabilityFlow(VectorType& rRightHandSideVector,
                                         InterfaceElementVariables& rVariables) override;

    void CalculateAndAddFluidBodyFlow(VectorType& rRightHandSideVector,
                                      InterfaceElementVariables& rVariables) override;

    unsigned int GetNumberOfDOF() const override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Member Variables

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
    }

    /// Assignment operator.
    TransientPwInterfaceElement & operator=(TransientPwInterfaceElement const& rOther);

    /// Copy constructor.
    TransientPwInterfaceElement(TransientPwInterfaceElement const& rOther);

}; // Class TransientPwInterfaceElement

} // namespace Kratos

#endif // KRATOS_GEO_TRANSIENT_PW_INTERFACE_ELEMENT_H_INCLUDED  defined
