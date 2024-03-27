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

#if !defined(KRATOS_GEO_TRUSS_ELEMENT_LINEAR_BASE_H_INCLUDED)
#define KRATOS_GEO_TRUSS_ELEMENT_LINEAR_BASE_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/geo_truss_element_base.hpp"
#include "includes/define.h"
#include "includes/variables.h"

namespace Kratos
{
/**
 * @class GeoTrussElementLinearBase
 *
 * @brief This is a linear truss element inheriting from GeoTrussElementBase
 *
 * @author Klaus B Sautter, Vahid Galavi
 */

template <unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoTrussElementLinearBase
    : public GeoTrussElementBase<TDim, TNumNodes>
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GeoTrussElementLinearBase);

    using BaseType          = GeoTrussElementBase<TDim, TNumNodes>;
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

    GeoTrussElementLinearBase(){};
    GeoTrussElementLinearBase(IndexType NewId, GeometryType::Pointer pGeometry);
    GeoTrussElementLinearBase(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ~GeoTrussElementLinearBase() override;

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

    void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief This function adds forces from prestressing to the force vector
     * @param rRightHandSideVector The right hand side of the problem
     */
    void AddPrestressLinear(VectorType& rRightHandSideVector);

    void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3>>& rVariable,
                                      std::vector<array_1d<double, 3>>&    rOutput,
                                      const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
                                      std::vector<Vector>&    rOutput,
                                      const ProcessInfo&      rCurrentProcessInfo) override;

    /**
     * @brief This function calculates the total stiffness matrix for the element
     */
    void CreateElementStiffnessMatrix(MatrixType&        rLocalStiffnessMatrix,
                                      const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief This function calculates the original nodal postion for the transformation matrix
     * @param rReferenceCoordinates The original coordinates
     */
    void WriteTransformationCoordinates(FullDofVectorType& rReferenceCoordinates) override;

    /**
     * @brief This function calculates the current linear-Lagrange strain
     */
    double CalculateLinearStrain();

    /**
     * @brief This function updates the internal forces
     * @param rinternalForces The internal forces
     */

    void UpdateInternalForces(FullDofVectorType& rInternalForces, const ProcessInfo& rCurrentProcessInfo) override;

    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

private:
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
    }
};

} // namespace Kratos

#endif
