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

#if !defined(KRATOS_GEO_CR_BEAM_ELEMENT_3D2N_H_INCLUDED)
#define KRATOS_GEO_CR_BEAM_ELEMENT_3D2N_H_INCLUDED

// System includes

// External includes

// Project includes
#include "../StructuralMechanicsApplication/custom_elements/cr_beam_element_3D2N.hpp"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "includes/variables.h"

namespace Kratos
{
/**
 * @class GeoCrBeamElement3D2N
 *
 * @brief This is a 3D-2node beam element with 3 translational dofs and 3 rotational dof per node,
 *        based on CrBeamElement3D2N element in Structural Mechanics application, modified to account for reset displacements
 *
 * @author Vahid Galavi
 */

class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoCrBeamElement3D2N : public CrBeamElement3D2N
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GeoCrBeamElement3D2N);

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

    GeoCrBeamElement3D2N(){};
    GeoCrBeamElement3D2N(IndexType NewId, GeometryType::Pointer pGeometry);
    GeoCrBeamElement3D2N(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ~GeoCrBeamElement3D2N() override;

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

    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;
    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    void ConstCalculateRightHandSide(VectorType&        rRightHandSideVector,
                                     const ProcessInfo& rCurrentProcessInfo) const override;

    void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3>>& rVariable,
                                      std::vector<array_1d<double, 3>>&    rOutput,
                                      const ProcessInfo& rCurrentProcessInfo) override;

    void ResetConstitutiveLaw() override;

protected:
    Vector mLocalForcesFinalized         = ZeroVector(msElementSize);
    Vector mLocalForcesFinalizedPrevious = ZeroVector(msElementSize);

private:
    bool mIsInitialization = false;

    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
};

} // namespace Kratos

#endif
