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
#include "custom_elements/geo_cr_beam_element_linear_2D2N.hpp"
#include "custom_utilities/static_condensation_utility.h"
#include "geo_mechanics_application_variables.h"
#include "includes/define.h"

namespace Kratos
{
GeoCrBeamElementLinear2D2N::GeoCrBeamElementLinear2D2N(IndexType NewId, GeometryType::Pointer pGeometry)
    : CrBeamElementLinear2D2N(NewId, pGeometry)
{
}

//-------------------------------------------------------------------------------------------------
GeoCrBeamElementLinear2D2N::GeoCrBeamElementLinear2D2N(IndexType               NewId,
                                                       GeometryType::Pointer   pGeometry,
                                                       PropertiesType::Pointer pProperties)
    : CrBeamElementLinear2D2N(NewId, pGeometry, pProperties)
{
}

//-------------------------------------------------------------------------------------------------
Element::Pointer GeoCrBeamElementLinear2D2N::Create(IndexType               NewId,
                                                    NodesArrayType const&   rThisNodes,
                                                    PropertiesType::Pointer pProperties) const
{
    const GeometryType& rGeom = GetGeometry();
    return Kratos::make_intrusive<GeoCrBeamElementLinear2D2N>(NewId, rGeom.Create(rThisNodes), pProperties);
}

//-------------------------------------------------------------------------------------------------
Element::Pointer GeoCrBeamElementLinear2D2N::Create(IndexType               NewId,
                                                    GeometryType::Pointer   pGeom,
                                                    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<GeoCrBeamElementLinear2D2N>(NewId, pGeom, pProperties);
}

//-------------------------------------------------------------------------------------------------
GeoCrBeamElementLinear2D2N::~GeoCrBeamElementLinear2D2N() {}

//----------------------------------------------------------------------------------------
void GeoCrBeamElementLinear2D2N::ResetConstitutiveLaw()
{
    KRATOS_TRY

    mInternalGlobalForcesFinalized         = ZeroVector(msElementSize);
    mInternalGlobalForcesFinalizedPrevious = ZeroVector(msElementSize);

    KRATOS_CATCH("")
}

//-------------------------------------------------------------------------------------------------
void GeoCrBeamElementLinear2D2N::CalculateLocalSystem(MatrixType&        rLeftHandSideMatrix,
                                                      VectorType&        rRightHandSideVector,
                                                      const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    // Kt
    CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
    Vector nodal_deformation = ZeroVector(msElementSize);
    GetValuesVector(nodal_deformation);

    rRightHandSideVector           = ZeroVector(msElementSize);
    noalias(mInternalGlobalForces) = prod(rLeftHandSideMatrix, nodal_deformation);

    noalias(rRightHandSideVector) -= (mInternalGlobalForces + mInternalGlobalForcesFinalizedPrevious);

    noalias(rRightHandSideVector) += CalculateBodyForces();

    KRATOS_CATCH("")
}

//-------------------------------------------------------------------------------------------------
void GeoCrBeamElementLinear2D2N::CalculateRightHandSide(VectorType&        rRightHandSideVector,
                                                        const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    Vector nodal_deformation = ZeroVector(msElementSize);
    GetValuesVector(nodal_deformation);
    rRightHandSideVector           = ZeroVector(msElementSize);
    noalias(mInternalGlobalForces) = prod(mK_Master, nodal_deformation);

    noalias(rRightHandSideVector) -= (mInternalGlobalForces + mInternalGlobalForcesFinalizedPrevious);
    noalias(rRightHandSideVector) += CalculateBodyForces();
    KRATOS_CATCH("")
}

/////////////////////////////////////////////////
///////////// CUSTOM FUNCTIONS --->>
/////////////////////////////////////////////////

void GeoCrBeamElementLinear2D2N::CalculateOnIntegrationPoints(const Variable<array_1d<double, 3>>& rVariable,
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

    Matrix left_hand_side_matrix = CreateElementStiffnessMatrix_Total();

    Vector nodal_deformation = ZeroVector(msElementSize);
    GetValuesVector(nodal_deformation);

    BoundedMatrix<double, msElementSize, msElementSize> transformation_matrix = CreateRotationMatrix();
    // calculate local displacements
    nodal_deformation = prod((trans(transformation_matrix)), nodal_deformation);

    //// start static back condensation
    if (Has(CONDENSED_DOF_LIST)) {
        Vector           dof_list_input = GetValue(CONDENSED_DOF_LIST);
        std::vector<int> dofList(dof_list_input.size());
        for (SizeType i = 0; i < dof_list_input.size(); ++i) {
            dofList[i] = dof_list_input[i];
        }
        Vector nodal_deformation_temp = nodal_deformation;
        GeoStaticCondensationUtility::ConvertingCondensation(
            *this, nodal_deformation_temp, nodal_deformation, dofList, left_hand_side_matrix);
    }
    //// end static back condensation

    Vector stress = prod(left_hand_side_matrix, nodal_deformation);
    stress += mInternalGlobalForcesFinalizedPrevious;

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
    }
    if (rVariable == FORCE) {
        rOutput[0][0] = -1.0 * stress[0] * 0.75 + stress[3] * 0.25;
        rOutput[1][0] = -1.0 * stress[0] * 0.50 + stress[3] * 0.50;
        rOutput[2][0] = -1.0 * stress[0] * 0.25 + stress[3] * 0.75;

        rOutput[0][1] = -1.0 * stress[1] * 0.75 + stress[4] * 0.25;
        rOutput[1][1] = -1.0 * stress[1] * 0.50 + stress[4] * 0.50;
        rOutput[2][1] = -1.0 * stress[1] * 0.25 + stress[4] * 0.75;

        rOutput[0][2] = 0.00;
        rOutput[1][2] = 0.00;
        rOutput[2][2] = 0.00;
    }

    KRATOS_CATCH("")
}

//-------------------------------------------------------------------------------------------------
void GeoCrBeamElementLinear2D2N::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
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

//----------------------------------------------------------------------------------------------------
void GeoCrBeamElementLinear2D2N::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    mIsInitialization = true;

    KRATOS_CATCH("")
}

//-------------------------------------------------------------------------------------------------
void GeoCrBeamElementLinear2D2N::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    noalias(mInternalGlobalForcesFinalized) = mInternalGlobalForces + mInternalGlobalForcesFinalizedPrevious;

    KRATOS_CATCH("")
}

//-------------------------------------------------------------------------------------------------
void GeoCrBeamElementLinear2D2N::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, CrBeamElementLinear2D2N)
    rSerializer.save("InternalGlobalForcesFinalized", mInternalGlobalForcesFinalized);
    rSerializer.save("InternalGlobalForcesFinalizedPrevious", mInternalGlobalForcesFinalizedPrevious);
}

void GeoCrBeamElementLinear2D2N::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, CrBeamElementLinear2D2N)
    rSerializer.load("InternalGlobalForcesFinalized", mInternalGlobalForcesFinalized);
    rSerializer.load("InternalGlobalForcesFinalizedPrevious", mInternalGlobalForcesFinalizedPrevious);
}

} // namespace Kratos.
