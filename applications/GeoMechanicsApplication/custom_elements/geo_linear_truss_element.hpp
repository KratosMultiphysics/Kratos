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

#if !defined(KRATOS_GEO_LINEAR_TRUSS_ELEMENT_H_INCLUDED)
#define KRATOS_GEO_LINEAR_TRUSS_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/geo_linear_truss_element_base.hpp"

namespace Kratos
{
/**
 * @class GeoLinearTrussElement
 *
 * @brief This is a linear truss element
 *
 * @author Vahid Galavi
 */

template <unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoLinearTrussElement
    : public GeoTrussElementLinearBase<TDim, TNumNodes>
{
protected:
    // const values
    static constexpr int mStressVectorSize                  = 1;
    Vector               mInternalStresses                  = ZeroVector(mStressVectorSize);
    Vector               mInternalStressesFinalized         = ZeroVector(mStressVectorSize);
    Vector               mInternalStressesFinalizedPrevious = ZeroVector(mStressVectorSize);

public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GeoLinearTrussElement);

    using BaseType          = GeoTrussElementLinearBase<TDim, TNumNodes>;
    using GeometryType      = Element::GeometryType;
    using NodesArrayType    = Element::NodesArrayType;
    using PropertiesType    = Element::PropertiesType;
    using IndexType         = Element::IndexType;
    using SizeType          = Element::SizeType;
    using MatrixType        = Element::MatrixType;
    using VectorType        = Element::VectorType;
    using FullDofMatrixType = typename GeoTrussElementBase<TDim, TNumNodes>::FullDofMatrixType;
    using FullDofVectorType = typename GeoTrussElementBase<TDim, TNumNodes>::FullDofVectorType;

    using GeoTrussElementBase<TDim, TNumNodes>::mpConstitutiveLaw;

    GeoLinearTrussElement(){};
    GeoLinearTrussElement(IndexType NewId, GeometryType::Pointer pGeometry);
    GeoLinearTrussElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ~GeoLinearTrussElement() override;

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

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief This function updates the internal normal force w.r.t. the current deformations
     * @param rinternalForces The current updated internal forces
     */
    void UpdateInternalForces(FullDofVectorType& rInternalForces, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3>>& rVariable,
                                      std::vector<array_1d<double, 3>>&    rOutput,
                                      const ProcessInfo& rCurrentProcessInfo) override;

    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    void ResetConstitutiveLaw() override;

private:
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
        rSerializer.save("InternalStresses", mInternalStresses);
        rSerializer.save("InternalStressesFinalized", mInternalStressesFinalized);
        rSerializer.save("InternalStressesFinalizedPrevious", mInternalStressesFinalizedPrevious);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
        rSerializer.load("InternalStresses", mInternalStresses);
        rSerializer.load("InternalStressesFinalized", mInternalStressesFinalized);
        rSerializer.load("InternalStressesFinalizedPrevious", mInternalStressesFinalizedPrevious);
    }
};
} // namespace Kratos

#endif
