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

#if !defined(KRATOS_GEO_CR_BEAM_ELEMENT_LINEAR_3D2N_H_INCLUDED)
#define KRATOS_GEO_CR_BEAM_ELEMENT_LINEAR_3D2N_H_INCLUDED

// System includes

// External includes

// Project includes
#include "../StructuralMechanicsApplication/custom_elements/cr_beam_element_linear_3D2N.hpp"
#include "custom_elements/geo_cr_beam_element_3D2N.hpp"
#include "includes/define.h"
#include "includes/variables.h"

namespace Kratos
{
/**
 * @class GeoCrBeamElementLinear3D2N
 *
 * @brief This is a linear 3D-2node beam element with 3 translational dofs and 3 rotational dof per node inheriting from CrBeamElement3D2N
 *
 * @author Vahid Galavi
 */

class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoCrBeamElementLinear3D2N : public CrBeamElement3D2N
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GeoCrBeamElementLinear3D2N);

    using BaseType             = CrBeamElement3D2N;
    using GeometryType         = BaseType::GeometryType;
    using NodesArrayType       = BaseType::NodesArrayType;
    using PropertiesType       = BaseType::PropertiesType;
    using IndexType            = BaseType::IndexType;
    using SizeType             = BaseType::SizeType;
    using MatrixType           = BaseType::MatrixType;
    using VectorType           = BaseType::VectorType;
    using EquationIdVectorType = BaseType::EquationIdVectorType;
    using DofsVectorType       = BaseType::DofsVectorType;

    GeoCrBeamElementLinear3D2N(){};
    GeoCrBeamElementLinear3D2N(IndexType NewId, GeometryType::Pointer pGeometry);
    GeoCrBeamElementLinear3D2N(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ~GeoCrBeamElementLinear3D2N() override;

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

    void CalculateLocalSystem(MatrixType&        rLeftHandSideMatrix,
                              VectorType&        rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3>>& rVariable,
                                      std::vector<array_1d<double, 3>>&    rOutput,
                                      const ProcessInfo& rCurrentProcessInfo) override;

    void ResetConstitutiveLaw() override;

protected:
    Vector mInternalGlobalForcesFinalized         = ZeroVector(msElementSize);
    Vector mInternalGlobalForcesFinalizedPrevious = ZeroVector(msElementSize);

private:
    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
};

} // namespace Kratos

#endif
