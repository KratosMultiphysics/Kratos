// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Klaus B. Sautter,
//                   Vahid Galavi
//

#if !defined(KRATOS_GEO_CABLE_ELEMENT_H_INCLUDED)
#define KRATOS_GEO_CABLE_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes

#include "custom_elements/geo_truss_element.hpp"

namespace Kratos
{
/**
 * @class GeoCableElement
 *
 * @brief This is a cable element inheriting from the GeoTrussElement
 *
 * @author Vahid Galavi
 */
template <unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoCableElement : public GeoTrussElement<TDim, TNumNodes>
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GeoCableElement);

    using BaseType = GeoTrussElement<TDim, TNumNodes>;

    using GeometryType   = Element::GeometryType;
    using NodesArrayType = Element::NodesArrayType;
    using PropertiesType = Element::PropertiesType;
    using IndexType      = Element::IndexType;
    using SizeType       = Element::SizeType;
    using MatrixType     = Element::MatrixType;
    using VectorType     = Element::VectorType;

    using FullDofMatrixType = typename GeoTrussElementBase<TDim, TNumNodes>::FullDofMatrixType;
    using FullDofVectorType = typename GeoTrussElementBase<TDim, TNumNodes>::FullDofVectorType;

    using GeoTrussElementBase<TDim, TNumNodes>::mpConstitutiveLaw;
    using GeoTrussElement<TDim, TNumNodes>::mInternalStressesFinalizedPrevious;
    using GeoTrussElement<TDim, TNumNodes>::mInternalStresses;

    GeoCableElement(){};
    GeoCableElement(IndexType NewId, GeometryType::Pointer pGeometry);
    GeoCableElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ~GeoCableElement() override;

    /**
     * @brief Creates a new element
     * @param NewId The Id of the new created element
     * @param pGeom The pointer to the geometry of the element
     * @param pProperties The pointer to property
     * @return The pointer to the created element
     */
    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

    /**
     * @brief Creates a new element
     * @param NewId The Id of the new created element
     * @param ThisNodes The array containing nodes
     * @param pProperties The pointer to property
     * @return The pointer to the created element
     */
    Element::Pointer Create(IndexType               NewId,
                            NodesArrayType const&   ThisNodes,
                            PropertiesType::Pointer pProperties) const override;

    void CreateElementStiffnessMatrix(MatrixType&        rLocalStiffnessMatrix,
                                      const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief This function updates the internal normal force w.r.t. the current deformations
     * @param rinternalForces The current updated internal forces
     */
    void UpdateInternalForces(BoundedVector<double, TDim * TNumNodes>& rinternalForces,
                              const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3>>& rVariable,
                                      std::vector<array_1d<double, 3>>&    rOutput,
                                      const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
                                      std::vector<Vector>&    rOutput,
                                      const ProcessInfo&      rCurrentProcessInfo) override;

private:
    // boolean for the cable --> does not resist to compression
    bool mIsCompressed;

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
        rSerializer.save("mIscompressed", mIsCompressed);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
        rSerializer.load("mIscompressed", mIsCompressed);
    }
};
} // namespace Kratos

#endif
