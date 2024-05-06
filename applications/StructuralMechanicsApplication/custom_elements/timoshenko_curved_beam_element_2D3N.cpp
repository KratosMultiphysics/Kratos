// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:  Alejandro Cornejo
//
//

// System includes

// External includes

// Project includes

// Application includes
#include "custom_elements/timoshenko_curved_beam_element_2D3N.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

void LinearTimoshenkoCurvedBeamElement2D3N::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Initialization should not be done again in a restart!
    if (!rCurrentProcessInfo[IS_RESTARTED]) {
        if (this->UseGeometryIntegrationMethod()) {
            if (GetProperties().Has(INTEGRATION_ORDER) ) {
                mThisIntegrationMethod = static_cast<GeometryData::IntegrationMethod>(GetProperties()[INTEGRATION_ORDER] - 1);
            } else {
                mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_5;
            }
        }

        const auto& r_integration_points = this->IntegrationPoints(mThisIntegrationMethod);

        // Constitutive Law initialisation
        if (mConstitutiveLawVector.size() != r_integration_points.size())
            mConstitutiveLawVector.resize(r_integration_points.size());
        InitializeMaterial();
    }
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::InitializeMaterial()
{
    KRATOS_TRY

    if (GetProperties()[CONSTITUTIVE_LAW] != nullptr) {
        const auto& r_geometry   = GetGeometry();
        const auto& r_properties = GetProperties();
        auto N_values            = Vector();
        for (IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number) {
            mConstitutiveLawVector[point_number] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
            mConstitutiveLawVector[point_number]->InitializeMaterial(r_properties, r_geometry, N_values);
        }
    } else
        KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer LinearTimoshenkoCurvedBeamElement2D3N::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    LinearTimoshenkoCurvedBeamElement2D3N::Pointer p_new_elem = Kratos::make_intrusive<LinearTimoshenkoCurvedBeamElement2D3N>
        (NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    // Currently selected integration methods
    p_new_elem->SetIntegrationMethod(mThisIntegrationMethod);

    // The vector containing the constitutive laws
    p_new_elem->SetConstitutiveLawVector(mConstitutiveLawVector);

    return p_new_elem;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    const auto& r_geometry = this->GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dofs_per_node = DoFperNode; // u, v, theta

    IndexType local_index = 0;

    if (rResult.size() != dofs_per_node * number_of_nodes)
        rResult.resize(dofs_per_node * number_of_nodes, false);

    const IndexType xpos    = this->GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);
    const IndexType rot_pos = this->GetGeometry()[0].GetDofPosition(ROTATION_Z);

    for (IndexType i = 0; i < number_of_nodes; ++i) {
        rResult[local_index++] = r_geometry[i].GetDof(DISPLACEMENT_X, xpos    ).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(DISPLACEMENT_Y, xpos + 1).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(ROTATION_Z    , rot_pos ).EquationId();
    }
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    const auto& r_geom = GetGeometry();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType dofs_per_node = DoFperNode; // u, v, theta
    rElementalDofList.resize(dofs_per_node * number_of_nodes);

    for (IndexType i = 0; i < number_of_nodes; ++i) {
        const SizeType index = i * dofs_per_node;
        rElementalDofList[index]     = r_geom[i].pGetDof(DISPLACEMENT_X);
        rElementalDofList[index + 1] = r_geom[i].pGetDof(DISPLACEMENT_Y);
        rElementalDofList[index + 2] = r_geom[i].pGetDof(ROTATION_Z    );
    }
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

const double LinearTimoshenkoCurvedBeamElement2D3N::GetJacobian()
{
    // GlobalSizeVector r_N;
    // GetFirstDerivativesNu0ShapeFunctionsValues(r_N, )
    // const double dx_dxi = 
    // return std::sqrt()
    return 0.0;
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::GetShapeFunctionsValues(
    GlobalSizeVector& rN,
    const double J,
    const double ShearFactor,
    const double xi
    )
{
    // todo
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::GetFirstDerivativesShapeFunctionsValues(
    GlobalSizeVector& rN,
    const double J,
    const double ShearFactor,
    const double xi
    )
{
    // todo
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::GetSecondDerivativesShapeFunctionsValues(
    GlobalSizeVector& rN,
    const double J,
    const double ShearFactor,
    const double xi
    )
{
    // todo
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::GetThirdDerivativesShapeFunctionsValues(
    GlobalSizeVector& rN,
    const double J,
    const double ShearFactor,
    const double xi
    )
{
    // todo
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::GetNThetaShapeFunctionsValues(
    GlobalSizeVector& rN,
    const double J,
    const double ShearFactor,
    const double k0,
    const double xi
    )
{
    // todo
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::GetFirstDerivativesNThetaShapeFunctionsValues(
    GlobalSizeVector& rN,
    const double J,
    const double ShearFactor,
    const double k0,
    const double xi
    )
{
    // todo
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::GetNu0ShapeFunctionsValues(
    GlobalSizeVector& rNu,
    const double xi
    )
{
    rNu.clear();
    rNu[0] = 0.5 * xi * (xi - 1.0);
    rNu[3] = 0.5 * xi * (xi + 1.0);
    rNu[6] = 1.0 - std::pow(xi, 2);
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::GetFirstDerivativesNu0ShapeFunctionsValues(
    GlobalSizeVector& rNu,
    const double J,
    const double xi
    )
{
    GetLocalFirstDerivativesNu0ShapeFunctionsValues(rNu, xi);
    rNu /= J;
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::GetLocalFirstDerivativesNu0ShapeFunctionsValues(
    GlobalSizeVector& rNu,
    const double xi
    )
{
    rNu.clear();
    rNu[0] = xi - 0.5;
    rNu[3] = xi + 0.5;
    rNu[6] = -2.0 * xi;
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::GetLocalSecondDerivativesNu0ShapeFunctionsValues(
    GlobalSizeVector& rNu,
    const double xi
    )
{
    rNu.clear();
    rNu[0] = 1.0;
    rNu[3] = 1.0;
    rNu[6] = -2.0;
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::GetSecondDerivativesNu0ShapeFunctionsValues(
    GlobalSizeVector& rNu,
    const double J,
    const double xi
    )
{
    GetLocalSecondDerivativesNu0ShapeFunctionsValues(rNu, xi);
    rNu /= std::pow(J, 2);
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::GetNodalValuesVector(GlobalSizeVector& rNodalValues)
{
    const auto &r_geom = GetGeometry();
    const auto& r_displ_0 = r_geom[0].FastGetSolutionStepValue(DISPLACEMENT);
    const auto& r_displ_1 = r_geom[1].FastGetSolutionStepValue(DISPLACEMENT);
    const auto& r_displ_2 = r_geom[2].FastGetSolutionStepValue(DISPLACEMENT);

    rNodalValues[0] = r_displ_0[0];
    rNodalValues[1] = r_displ_0[1];
    rNodalValues[2] = r_geom[0].FastGetSolutionStepValue(ROTATION_Z);
    rNodalValues[3] = r_displ_1[0];
    rNodalValues[4] = r_displ_1[1];
    rNodalValues[5] = r_geom[1].FastGetSolutionStepValue(ROTATION_Z);
    rNodalValues[6] = r_displ_2[0];
    rNodalValues[7] = r_displ_2[1];
    rNodalValues[8] = r_geom[2].FastGetSolutionStepValue(ROTATION_Z);
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::CalculateGeneralizedStrainsVector(
    VectorType& rStrain,
    const double J,
    const double ShearFactor,
    const double xi,
    const GlobalSizeVector &rNodalValues
    )
{
    // todo
}

/***********************************************************************************/
/***********************************************************************************/

double LinearTimoshenkoCurvedBeamElement2D3N::CalculateAxialStrain(
    const double J,
    const double ShearFactor,
    const double xi,
    const GlobalSizeVector& rNodalValues
    )
{
    return 0.0;
}

/***********************************************************************************/
/***********************************************************************************/

double LinearTimoshenkoCurvedBeamElement2D3N::CalculateShearStrain(
    const double J,
    const double ShearFactor,
    const double xi,
    const GlobalSizeVector& rNodalValues
    )
{
    return 0.0;
}

/***********************************************************************************/
/***********************************************************************************/

double LinearTimoshenkoCurvedBeamElement2D3N::CalculateBendingCurvature(
    const double J,
    const double ShearFactor,
    const double xi,
    const GlobalSizeVector& rNodalValues
    )
{
    return 0.0;
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 3> LinearTimoshenkoCurvedBeamElement2D3N::GetLocalAxesBodyForce(
    const Element &rElement,
    const GeometryType::IntegrationPointsArrayType &rIntegrationPoints,
    const IndexType PointNumber
    )
{
    return array_1d<double, 3>();
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::CalculateLocalSystem(
    MatrixType& rLHS,
    VectorType& rRHS,
    const ProcessInfo& rProcessInfo
    )
{

}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::CalculateLeftHandSide(
    MatrixType& rLHS,
    const ProcessInfo& rProcessInfo
    )
{

}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::CalculateRightHandSide(
    VectorType& rRHS,
    const ProcessInfo& rProcessInfo
    )
{

}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rProcessInfo
    )
{

}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::CalculateOnIntegrationPoints(
    const Variable<ConstitutiveLaw::Pointer>& rVariable,
    std::vector<ConstitutiveLaw::Pointer>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (rVariable == CONSTITUTIVE_LAW) {
        const SizeType integration_points_number = mConstitutiveLawVector.size();
        if (rValues.size() != integration_points_number) {
            rValues.resize(integration_points_number);
        }
        for (IndexType point_number = 0; point_number < integration_points_number; ++point_number) {
            rValues[point_number] = mConstitutiveLawVector[point_number];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

int LinearTimoshenkoCurvedBeamElement2D3N::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    return mConstitutiveLawVector[0]->Check(GetProperties(), GetGeometry(), rCurrentProcessInfo);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    int IntMethod = int(this->GetIntegrationMethod());
    rSerializer.save("IntegrationMethod",IntMethod);
    rSerializer.save("ConstitutiveLawVector", mConstitutiveLawVector);
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    int IntMethod;
    rSerializer.load("IntegrationMethod",IntMethod);
    mThisIntegrationMethod = IntegrationMethod(IntMethod);
    rSerializer.load("ConstitutiveLawVector", mConstitutiveLawVector);
}

/***********************************************************************************/
/***********************************************************************************/

} // Namespace Kratos
