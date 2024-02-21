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

// System includes

// External includes

// Project includes
#include "custom_elements/geo_cr_beam_element_2D2N.hpp"
#include "geo_mechanics_application_variables.h"
#include "includes/define.h"

namespace Kratos
{
GeoCrBeamElement2D2N::GeoCrBeamElement2D2N(IndexType NewId, GeometryType::Pointer pGeometry)
    : CrBeamElement2D2N(NewId, pGeometry)
{
}

//----------------------------------------------------------------------------------------------------
GeoCrBeamElement2D2N::GeoCrBeamElement2D2N(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : CrBeamElement2D2N(NewId, pGeometry, pProperties)
{
}

//----------------------------------------------------------------------------------------------------
Element::Pointer GeoCrBeamElement2D2N::Create(IndexType               NewId,
                                              NodesArrayType const&   rThisNodes,
                                              PropertiesType::Pointer pProperties) const
{
    const GeometryType& rGeom = GetGeometry();
    return Kratos::make_intrusive<GeoCrBeamElement2D2N>(NewId, rGeom.Create(rThisNodes), pProperties);
}

//----------------------------------------------------------------------------------------------------
Element::Pointer GeoCrBeamElement2D2N::Create(IndexType               NewId,
                                              GeometryType::Pointer   pGeom,
                                              PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<GeoCrBeamElement2D2N>(NewId, pGeom, pProperties);
}

//----------------------------------------------------------------------------------------------------
GeoCrBeamElement2D2N::~GeoCrBeamElement2D2N() {}

//----------------------------------------------------------------------------------------------------
void GeoCrBeamElement2D2N::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (mIsInitialization) {
        if (rCurrentProcessInfo.Has(RESET_DISPLACEMENTS)) {
            if (rCurrentProcessInfo[RESET_DISPLACEMENTS])
                noalias(mInternalGlobalForcesFinalizedPrevious) = mInternalGlobalForcesFinalized;
            else noalias(mInternalGlobalForcesFinalized) = mInternalGlobalForcesFinalizedPrevious;
        } else {
            noalias(mInternalGlobalForcesFinalized)         = ZeroVector(msElementSize);
            noalias(mInternalGlobalForcesFinalizedPrevious) = ZeroVector(msElementSize);
        }
    }

    mIsInitialization = false;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
void GeoCrBeamElement2D2N::ResetConstitutiveLaw()
{
    KRATOS_TRY

    mInternalGlobalForcesFinalized         = ZeroVector(msElementSize);
    mInternalGlobalForcesFinalizedPrevious = ZeroVector(msElementSize);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------
void GeoCrBeamElement2D2N::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    mIsInitialization = true;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------
void GeoCrBeamElement2D2N::CalculateLocalSystem(MatrixType&        rLeftHandSideMatrix,
                                                VectorType&        rRightHandSideVector,
                                                const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    // t
    mDeformationForces = CalculateInternalStresses_DeformationModes();

    // qe
    Vector nodal_forces = ZeroVector(msElementSize);
    nodal_forces        = ReturnElementForces_Local();

    // q
    GlobalizeVector(nodal_forces);
    mInternalGlobalForces = nodal_forces;

    // Kt
    CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);

    // residual >>> r = f_ext - f_int
    rRightHandSideVector = ZeroVector(msElementSize);
    noalias(rRightHandSideVector) -= (mInternalGlobalForces + mInternalGlobalForcesFinalizedPrevious);

    noalias(rRightHandSideVector) += CalculateBodyForces();

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------
void GeoCrBeamElement2D2N::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    // t
    mDeformationForces = CalculateInternalStresses_DeformationModes();

    // qe
    Vector nodal_forces = ZeroVector(msElementSize);
    nodal_forces        = ReturnElementForces_Local();

    // q
    GlobalizeVector(nodal_forces);
    mInternalGlobalForces = nodal_forces;

    // residual >>> r = f_ext - f_int
    rRightHandSideVector = ZeroVector(msElementSize);
    noalias(rRightHandSideVector) -= (mInternalGlobalForces + mInternalGlobalForcesFinalizedPrevious);

    noalias(rRightHandSideVector) += CalculateBodyForces();
    KRATOS_CATCH("")
}

/////////////////////////////////////////////////
///////////// CUSTOM FUNCTIONS --->>
/////////////////////////////////////////////////

void GeoCrBeamElement2D2N::CalculateOnIntegrationPoints(const Variable<array_1d<double, 3>>& rVariable,
                                                        std::vector<array_1d<double, 3>>& rOutput,
                                                        const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    // Element with two nodes can only represent results at one node
    const auto&                                     r_geometry = GetGeometry();
    const GeometryType::IntegrationPointsArrayType& r_integration_points =
        r_geometry.IntegrationPoints(Kratos::GeometryData::IntegrationMethod::GI_GAUSS_3);
    const SizeType write_points_number = r_integration_points.size();
    if (rOutput.size() != write_points_number) {
        rOutput.resize(write_points_number);
    }

    BoundedMatrix<double, msElementSize, msElementSize> transformation_matrix = CreateRotationMatrix();
    Vector stress = mInternalGlobalForces + mInternalGlobalForcesFinalizedPrevious;
    stress        = prod(trans(transformation_matrix), stress);

    // rOutput[GP 1,2,3][x,y,z]

    if (rVariable == MOMENT) {
        rOutput[0][0] = 0.00;
        rOutput[1][0] = 0.00;
        rOutput[2][0] = 0.00;

        rOutput[0][1] = 0.00;
        rOutput[1][1] = 0.00;
        rOutput[2][1] = 0.00;

        rOutput[0][2] = 1.0 * stress[2] * 0.75 - stress[5] * 0.25;
        rOutput[1][2] = 1.0 * stress[2] * 0.50 - stress[5] * 0.50;
        rOutput[2][2] = 1.0 * stress[2] * 0.25 - stress[5] * 0.75;
    } else if (rVariable == FORCE) {
        rOutput[0][0] = -1.0 * stress[0] * 0.75 + stress[3] * 0.25;
        rOutput[1][0] = -1.0 * stress[0] * 0.50 + stress[3] * 0.50;
        rOutput[2][0] = -1.0 * stress[0] * 0.25 + stress[3] * 0.75;

        rOutput[0][1] = -1.0 * stress[1] * 0.75 + stress[4] * 0.25;
        rOutput[1][1] = -1.0 * stress[1] * 0.50 + stress[4] * 0.50;
        rOutput[2][1] = -1.0 * stress[1] * 0.25 + stress[4] * 0.75;

        rOutput[0][2] = 0.00;
        rOutput[1][2] = 0.00;
        rOutput[2][2] = 0.00;
    } else if (rVariable == INTEGRATION_COORDINATES) {
        Point global_point;
        for (IndexType point_number = 0; point_number < write_points_number; ++point_number) {
            r_geometry.GlobalCoordinates(global_point, r_integration_points[point_number]);
            rOutput[point_number] = global_point.Coordinates();
        }
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------
void GeoCrBeamElement2D2N::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    noalias(mInternalGlobalForcesFinalized) = mInternalGlobalForces + mInternalGlobalForcesFinalizedPrevious;
}

//----------------------------------------------------------------------------------------------------
void GeoCrBeamElement2D2N::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, CrBeamElement2D2N)
    rSerializer.save("InternalGlobalForcesFinalized", mInternalGlobalForcesFinalized);
    rSerializer.save("InternalGlobalForcesFinalizedPrevious", mInternalGlobalForcesFinalizedPrevious);
}

//----------------------------------------------------------------------------------------------------
void GeoCrBeamElement2D2N::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, CrBeamElement2D2N)
    rSerializer.load("InternalGlobalForcesFinalized", mInternalGlobalForcesFinalized);
    rSerializer.load("InternalGlobalForcesFinalizedPrevious", mInternalGlobalForcesFinalizedPrevious);
}
} // namespace Kratos.
