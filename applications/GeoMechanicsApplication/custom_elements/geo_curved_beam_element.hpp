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

#if !defined(KRATOS_GEO_CURVED_BEAM_ELEMENT_H_INCLUDED)
#define KRATOS_GEO_CURVED_BEAM_ELEMENT_H_INCLUDED

// Project includes

// Application includes
#include "custom_elements/geo_structural_base_element.hpp"

namespace Kratos
{
/**
 * @class GeoCurvedBeamElement
 *
 * @brief This is a geometrically non-linear (curved) beam element.
 *        The formulation can be found in papers written by Karan S. Surana, e.g:
 *        "1. Geometrically non-linear formulation for the axisymmetric shell elements"
 *        "2. Geometrically non-linear formulation for two dimensional curved beam elements"
 *        Discriptions of beam elements can be found in the following book, chapter 9.
 *        For 2D curved beams, see section 9.4:
 *        "Non-linear Finite element analysis of solids and structures" by De Borst et al.
 *
 * @author Vahid Galavi
 */

template <unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoCurvedBeamElement
    : public GeoStructuralBaseElement<TDim, TNumNodes>
{
public:
    using IndexType                   = std::size_t;
    using PropertiesType              = Properties;
    using NodeType                    = Node;
    using GeometryType                = Geometry<NodeType>;
    using NodesArrayType              = GeometryType::PointsArrayType;
    using VectorType                  = Vector;
    using MatrixType                  = Matrix;
    using ShapeFunctionsGradientsType = GeometryData::ShapeFunctionsGradientsType;

    /// The definition of the sizetype
    using SizeType = std::size_t;

    using GeoStructuralBaseElement<TDim, TNumNodes>::mThisIntegrationMethod;
    using GeoStructuralBaseElement<TDim, TNumNodes>::mConstitutiveLawVector;
    using GeoStructuralBaseElement<TDim, TNumNodes>::mStressVector;
    using GeoStructuralBaseElement<TDim, TNumNodes>::N_DOF_ELEMENT;
    using GeoStructuralBaseElement<TDim, TNumNodes>::VoigtSize;
    using GeoStructuralBaseElement<TDim, TNumNodes>::GetNodalDofValuesVector;
    using ElementVariables = typename GeoStructuralBaseElement<TDim, TNumNodes>::ElementVariables;

    KRATOS_CLASS_POINTER_DEFINITION(GeoCurvedBeamElement);

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    GeoCurvedBeamElement(IndexType NewId = 0) : GeoStructuralBaseElement<TDim, TNumNodes>(NewId) {}

    /// Constructor using an array of nodes
    GeoCurvedBeamElement(IndexType NewId, const NodesArrayType& ThisNodes)
        : GeoStructuralBaseElement<TDim, TNumNodes>(NewId, ThisNodes)
    {
    }

    /// Constructor using Geometry
    GeoCurvedBeamElement(IndexType NewId, GeometryType::Pointer pGeometry)
        : GeoStructuralBaseElement<TDim, TNumNodes>(NewId, pGeometry)
    {
    }

    /// Constructor using Properties
    GeoCurvedBeamElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : GeoStructuralBaseElement<TDim, TNumNodes>(NewId, pGeometry, pProperties)
    {
        mThisIntegrationMethod = this->GetIntegrationMethod();
    }

    /// Destructor
    virtual ~GeoCurvedBeamElement() {}

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Element::Pointer Create(IndexType               NewId,
                            NodesArrayType const&   ThisNodes,
                            PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                      std::vector<Matrix>&    rOutput,
                                      const ProcessInfo&      rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3>>& rVariable,
                                      std::vector<array_1d<double, 3>>&    rOutput,
                                      const ProcessInfo& rCurrentProcessInfo) override;

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:
    /// Member Variables

    SizeType GetCrossNumberIntegrationPoints() const override;
    SizeType GetAlongNumberIntegrationPoints() const override;

    double CalculateIntegrationCoefficient(unsigned int GPointCross, double detJ, double weight) const;

    virtual void CalculateBMatrix(Matrix&                                  B,
                                  unsigned int                             GPointCross,
                                  const BoundedMatrix<double, TDim, TDim>& InvertDetJacobian,
                                  ElementVariables&                        rVariables) const;

    virtual void CalculateLocalBMatrix(Matrix&                                  B,
                                       unsigned int                             GPointCross,
                                       const BoundedMatrix<double, TDim, TDim>& InvertDetJacobian,
                                       ElementVariables&                        rVariables) const;

    void CalculateStrainVector(ElementVariables& rVariables) const;

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    virtual void CalculateAll(MatrixType&        rLeftHandSideMatrix,
                              VectorType&        rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo,
                              const bool         CalculateStiffnessMatrixFlag,
                              const bool         CalculateResidualVectorFlag) override;

    virtual void CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables) const;

    virtual void CalculateAndAddRHS(VectorType&       rRightHandSideVector,
                                    ElementVariables& rVariables,
                                    unsigned int      GPoint) const;

    virtual void CalculateLocalInternalForce(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo);

    virtual void CalculateTransformationMatrix(Matrix& TransformationMatrix, const Matrix& GradNe) const;

    virtual void CalculateNodalCrossDirection(Matrix& NodalCrossDirection) const override;

    virtual double CalculateAngleAtGaussPoint(const Matrix& GradNe) const;

    virtual double CalculateAngleAtNode(unsigned int GPoint,
                                        const BoundedMatrix<double, TNumNodes, TNumNodes>& DN_DXContainer) const;

    virtual void CalculateJacobianMatrix(unsigned int            GPointCross,
                                         const ElementVariables& rVariables,
                                         BoundedMatrix<double, TDim, TDim>& DeterminantJacobian) const;

    void CalculateAndAddBodyForce(VectorType& rRightHandSideVector, ElementVariables& rVariables) const;

    void CalculateAndAddStiffnessForce(VectorType&       rRightHandSideVector,
                                       ElementVariables& rVariables,
                                       unsigned int      GPoint) const;
    void SetRotationalInertiaVector(const PropertiesType& Prop, Vector& RotationalInertia) const;

    void InitializeElementVariables(ElementVariables&            rVariables,
                                    ConstitutiveLaw::Parameters& rConstitutiveParameters,
                                    const GeometryType&          Geom,
                                    const PropertiesType&        Prop,
                                    const ProcessInfo& rCurrentProcessInfo) const override;

    void InterpolateOnOutputPoints(Matrix& Values) const;
    void InterpolateOnOutputPoints(Vector& Values) const;

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    /// Assignment operator.
    GeoCurvedBeamElement& operator=(GeoCurvedBeamElement const& rOther);

    /// Copy constructor.
    GeoCurvedBeamElement(GeoCurvedBeamElement const& rOther);

    /// Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element)
    }

    // const values
    static constexpr SizeType N_DOF_NODE_DISP = TDim;
    static constexpr SizeType N_DOF_NODE_ROT  = (TDim == 2 ? 1 : 3);
    static constexpr SizeType N_DOF_NODE      = N_DOF_NODE_DISP + N_DOF_NODE_ROT;
    static constexpr SizeType N_POINT_CROSS   = 2;

}; // Class GeoCurvedBeamElement

} // namespace Kratos

#endif // KRATOS_GEO_CURVED_BEAM_ELEMENT_H_INCLUDED  defined
