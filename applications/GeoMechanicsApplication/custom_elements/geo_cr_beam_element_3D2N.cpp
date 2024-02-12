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
#include "custom_elements/geo_cr_beam_element_3D2N.hpp"
#include "geo_mechanics_application_variables.h"
#include "includes/define.h"

namespace Kratos
{
GeoCrBeamElement3D2N::GeoCrBeamElement3D2N(IndexType NewId, GeometryType::Pointer pGeometry)
    : CrBeamElement3D2N(NewId, pGeometry)
{
}

//----------------------------------------------------------------------------------------------------
GeoCrBeamElement3D2N::GeoCrBeamElement3D2N(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : CrBeamElement3D2N(NewId, pGeometry, pProperties)
{
}

//----------------------------------------------------------------------------------------------------
Element::Pointer GeoCrBeamElement3D2N::Create(IndexType               NewId,
                                              NodesArrayType const&   rThisNodes,
                                              PropertiesType::Pointer pProperties) const
{
    const GeometryType& rGeom = GetGeometry();
    return Kratos::make_intrusive<GeoCrBeamElement3D2N>(NewId, rGeom.Create(rThisNodes), pProperties);
}

//----------------------------------------------------------------------------------------------------
Element::Pointer GeoCrBeamElement3D2N::Create(IndexType               NewId,
                                              GeometryType::Pointer   pGeom,
                                              PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<GeoCrBeamElement3D2N>(NewId, pGeom, pProperties);
}

//----------------------------------------------------------------------------------------------------
GeoCrBeamElement3D2N::~GeoCrBeamElement3D2N() {}

//-------------------------------------------------------------------------------------------------
void GeoCrBeamElement3D2N::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    mIsInitialization = true;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
void GeoCrBeamElement3D2N::ResetConstitutiveLaw()
{
    KRATOS_TRY

    mLocalForcesFinalized         = ZeroVector(msElementSize);
    mLocalForcesFinalizedPrevious = ZeroVector(msElementSize);

    KRATOS_CATCH("")
}

//-------------------------------------------------------------------------------------------------
void GeoCrBeamElement3D2N::ConstCalculateRightHandSide(VectorType&        rRightHandSideVector,
                                                       const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY
    Vector nodal_forces_local_qe = CalculateLocalNodalForces();

    BoundedMatrix<double, msElementSize, msElementSize> total_rotation_matrix =
        GetTransformationMatrixGlobal();

    // Nodal element forces global
    Vector internal_forces         = prod(total_rotation_matrix, nodal_forces_local_qe);
    Vector internal_previous_force = prod(total_rotation_matrix, mLocalForcesFinalizedPrevious);

    rRightHandSideVector = ZeroVector(msElementSize);
    noalias(rRightHandSideVector) -= internal_forces;
    noalias(rRightHandSideVector) -= internal_previous_force;

    // add bodyforces
    noalias(rRightHandSideVector) += CalculateBodyForces();
    KRATOS_CATCH("")
}

//-------------------------------------------------------------------------------------------------
void GeoCrBeamElement3D2N::CalculateOnIntegrationPoints(const Variable<array_1d<double, 3>>& rVariable,
                                                        std::vector<array_1d<double, 3>>& rOutput,
                                                        const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    // element with two nodes can only represent results at one node
    const unsigned int& write_points_number =
        GetGeometry().IntegrationPointsNumber(Kratos::GeometryData::IntegrationMethod::GI_GAUSS_3);
    if (rOutput.size() != write_points_number) {
        rOutput.resize(write_points_number);
    }

    // rOutput[GP 1,2,3][x,y,z]

    if (rVariable == MOMENT) {
        Vector nodal_forces_local_qe = CalculateLocalNodalForces();
        nodal_forces_local_qe += mLocalForcesFinalizedPrevious;
        rOutput[0][0] = -1.0 * nodal_forces_local_qe[3] * 0.75 + nodal_forces_local_qe[9] * 0.25;
        rOutput[1][0] = -1.0 * nodal_forces_local_qe[3] * 0.50 + nodal_forces_local_qe[9] * 0.50;
        rOutput[2][0] = -1.0 * nodal_forces_local_qe[3] * 0.25 + nodal_forces_local_qe[9] * 0.75;

        rOutput[0][1] = -1.0 * nodal_forces_local_qe[4] * 0.75 + nodal_forces_local_qe[10] * 0.25;
        rOutput[1][1] = -1.0 * nodal_forces_local_qe[4] * 0.50 + nodal_forces_local_qe[10] * 0.50;
        rOutput[2][1] = -1.0 * nodal_forces_local_qe[4] * 0.25 + nodal_forces_local_qe[10] * 0.75;

        rOutput[0][2] = 1.0 * nodal_forces_local_qe[5] * 0.75 - nodal_forces_local_qe[11] * 0.25;
        rOutput[1][2] = 1.0 * nodal_forces_local_qe[5] * 0.50 - nodal_forces_local_qe[11] * 0.50;
        rOutput[2][2] = 1.0 * nodal_forces_local_qe[5] * 0.25 - nodal_forces_local_qe[11] * 0.75;
    } else if (rVariable == FORCE) {
        Vector nodal_forces_local_qe = CalculateLocalNodalForces();
        nodal_forces_local_qe += mLocalForcesFinalizedPrevious;
        rOutput[0][0] = -1.0 * nodal_forces_local_qe[0] * 0.75 + nodal_forces_local_qe[6] * 0.25;
        rOutput[1][0] = -1.0 * nodal_forces_local_qe[0] * 0.50 + nodal_forces_local_qe[6] * 0.50;
        rOutput[2][0] = -1.0 * nodal_forces_local_qe[0] * 0.25 + nodal_forces_local_qe[6] * 0.75;

        rOutput[0][1] = -1.0 * nodal_forces_local_qe[1] * 0.75 + nodal_forces_local_qe[7] * 0.25;
        rOutput[1][1] = -1.0 * nodal_forces_local_qe[1] * 0.50 + nodal_forces_local_qe[7] * 0.50;
        rOutput[2][1] = -1.0 * nodal_forces_local_qe[1] * 0.25 + nodal_forces_local_qe[7] * 0.75;

        rOutput[0][2] = -1.0 * nodal_forces_local_qe[2] * 0.75 + nodal_forces_local_qe[8] * 0.25;
        rOutput[1][2] = -1.0 * nodal_forces_local_qe[2] * 0.50 + nodal_forces_local_qe[8] * 0.50;
        rOutput[2][2] = -1.0 * nodal_forces_local_qe[2] * 0.25 + nodal_forces_local_qe[8] * 0.75;
    } else if (rVariable == LOCAL_AXIS_1) {
        BoundedMatrix<double, msElementSize, msElementSize> rotation_matrix = GetTransformationMatrixGlobal();
        for (SizeType i = 0; i < msDimension; ++i) {
            rOutput[1][i] = column(rotation_matrix, 0)[i];
        }
    } else if (rVariable == LOCAL_AXIS_2) {
        BoundedMatrix<double, msElementSize, msElementSize> rotation_matrix = GetTransformationMatrixGlobal();
        for (SizeType i = 0; i < msDimension; ++i) {
            rOutput[1][i] = column(rotation_matrix, 1)[i];
        }
    } else if (rVariable == LOCAL_AXIS_3) {
        BoundedMatrix<double, msElementSize, msElementSize> rotation_matrix = GetTransformationMatrixGlobal();
        for (SizeType i = 0; i < msDimension; ++i) {
            rOutput[1][i] = column(rotation_matrix, 2)[i];
        }
    }

    KRATOS_CATCH("")
}

//-------------------------------------------------------------------------------------------------
void GeoCrBeamElement3D2N::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (mIsInitialization) {
        if (rCurrentProcessInfo.Has(RESET_DISPLACEMENTS)) {
            if (rCurrentProcessInfo[RESET_DISPLACEMENTS])
                noalias(mLocalForcesFinalizedPrevious) = mLocalForcesFinalized;
            else noalias(mLocalForcesFinalized) = mLocalForcesFinalizedPrevious;
        } else {
            noalias(mLocalForcesFinalized)         = ZeroVector(msElementSize);
            noalias(mLocalForcesFinalizedPrevious) = ZeroVector(msElementSize);
        }
    }
    mIsInitialization = false;

    KRATOS_CATCH("")
}

//-------------------------------------------------------------------------------------------------
void GeoCrBeamElement3D2N::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    noalias(mLocalForcesFinalized) = CalculateLocalNodalForces() + mLocalForcesFinalizedPrevious;

    KRATOS_CATCH("")
}

//-------------------------------------------------------------------------------------------------
void GeoCrBeamElement3D2N::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, CrBeamElement3D2N)
    rSerializer.save("LocalForcesFinalized", mLocalForcesFinalized);
    rSerializer.save("LocalForcesFinalizedPrevious", mLocalForcesFinalizedPrevious);
}

//-------------------------------------------------------------------------------------------------
void GeoCrBeamElement3D2N::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, CrBeamElement3D2N)
    rSerializer.load("LocalForcesFinalized", mLocalForcesFinalized);
    rSerializer.load("LocalForcesFinalizedPrevious", mLocalForcesFinalizedPrevious);
}

} // namespace Kratos.
