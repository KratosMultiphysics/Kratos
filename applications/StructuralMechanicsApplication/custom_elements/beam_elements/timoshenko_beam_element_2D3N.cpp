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
#include "timoshenko_beam_element_3D2N.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer LinearTimoshenkoBeamElement3D2N::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    LinearTimoshenkoBeamElement3D2N::Pointer p_new_elem = Kratos::make_intrusive<LinearTimoshenkoBeamElement3D2N>
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

void LinearTimoshenkoBeamElement3D2N::GetNodalValuesVector(
    VectorType& rNodalValues
    ) const
{




}

/***********************************************************************************/
/***********************************************************************************/

double LinearTimoshenkoBeamElement3D2N::CalculateAxialStrain(
    const double Length,
    const double Phi,
    const double xi,
    const VectorType& rNodalValues
    ) const
{
    return 0.0;


}

/***********************************************************************************/
/***********************************************************************************/

double LinearTimoshenkoBeamElement3D2N::CalculateShearStrainXY(
    const double Length,
    const double Phi,
    const double xi,
    const VectorType& rNodalValues
    ) const
{
    // VectorType N_derivatives(6), N_theta(6);
    // GetFirstDerivativesShapeFunctionsValues(N_derivatives, Length, Phi, xi);
    // GetNThetaShapeFunctionsValues(N_theta, Length, Phi, xi);
    // const VectorType N_s = N_derivatives - N_theta;
    // return N_s[0] * rNodalValues[1] + N_s[1] * rNodalValues[2] + N_s[4] * rNodalValues[4] +
    //        N_s[5] * rNodalValues[5] + N_s[2] * rNodalValues[7] + N_s[3] * rNodalValues[8];
    return 0.0;
}

/***********************************************************************************/
/***********************************************************************************/

double LinearTimoshenkoBeamElement3D2N::CalculateShearStrainXZ(
    const double Length,
    const double Phi,
    const double xi,
    const VectorType& rNodalValues
    ) const
{
    // VectorType N_derivatives(6), N_theta(6);
    // GetFirstDerivativesShapeFunctionsValues(N_derivatives, Length, Phi, xi);
    // GetNThetaShapeFunctionsValues(N_theta, Length, Phi, xi);
    // const VectorType N_s = N_derivatives - N_theta;
    // return N_s[0] * rNodalValues[1] + N_s[1] * rNodalValues[2] + N_s[4] * rNodalValues[4] +
    //        N_s[5] * rNodalValues[5] + N_s[2] * rNodalValues[7] + N_s[3] * rNodalValues[8];
    return 0.0;
}

/***********************************************************************************/
/***********************************************************************************/

double LinearTimoshenkoBeamElement3D2N::CalculateBendingCurvatureX(
    const double Length,
    const double Phi,
    const double xi,
    const VectorType& rNodalValues
    ) const
{
    // VectorType N_theta_derivatives(6);
    // GetFirstDerivativesNThetaShapeFunctionsValues(N_theta_derivatives, Length, Phi, xi);
    // return N_theta_derivatives[0] * rNodalValues[1] + N_theta_derivatives[1] * rNodalValues[2] +
    //        N_theta_derivatives[4] * rNodalValues[4] + N_theta_derivatives[5] * rNodalValues[5] +
    //        N_theta_derivatives[2] * rNodalValues[7] + N_theta_derivatives[3] * rNodalValues[8];
    return 0.0;
}

/***********************************************************************************/
/***********************************************************************************/

double LinearTimoshenkoBeamElement3D2N::CalculateBendingCurvatureY(
    const double Length,
    const double Phi,
    const double xi,
    const VectorType& rNodalValues
    ) const
{
    // VectorType N_theta_derivatives(6);
    // GetFirstDerivativesNThetaShapeFunctionsValues(N_theta_derivatives, Length, Phi, xi);
    // return N_theta_derivatives[0] * rNodalValues[1] + N_theta_derivatives[1] * rNodalValues[2] +
    //        N_theta_derivatives[4] * rNodalValues[4] + N_theta_derivatives[5] * rNodalValues[5] +
    //        N_theta_derivatives[2] * rNodalValues[7] + N_theta_derivatives[3] * rNodalValues[8];
    return 0.0;
}

/***********************************************************************************/
/***********************************************************************************/

double LinearTimoshenkoBeamElement3D2N::CalculateBendingCurvatureZ(
    const double Length,
    const double Phi,
    const double xi,
    const VectorType& rNodalValues
    ) const
{
    // VectorType N_theta_derivatives(6);
    // GetFirstDerivativesNThetaShapeFunctionsValues(N_theta_derivatives, Length, Phi, xi);
    // return N_theta_derivatives[0] * rNodalValues[1] + N_theta_derivatives[1] * rNodalValues[2] +
    //        N_theta_derivatives[4] * rNodalValues[4] + N_theta_derivatives[5] * rNodalValues[5] +
    //        N_theta_derivatives[2] * rNodalValues[7] + N_theta_derivatives[3] * rNodalValues[8];
    return 0.0;
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement3D2N::RotateLHS(
    MatrixType& rLHS,
    const GeometryType& rGeometry
    )
{
    // const double angle = GetAngle();

    // if (std::abs(angle) > std::numeric_limits<double>::epsilon()) {
    //     BoundedMatrix<double, 3, 3> T, Tt;
    //     BoundedMatrix<double, 9, 9> global_size_T, aux_product;
    //     StructuralMechanicsElementUtilities::BuildRotationMatrixForBeam(T, angle);
    //     StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor2D3NBeam(T, global_size_T);
    //     noalias(aux_product) = prod(rLHS, trans(global_size_T));
    //     noalias(rLHS) = prod(global_size_T, aux_product);
    // }
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement3D2N::RotateRHS(
    VectorType& rRHS,
    const GeometryType& rGeometry
    )
{
    // const double angle = GetAngle();
    // if (std::abs(angle) > std::numeric_limits<double>::epsilon()) {
    //     BoundedMatrix<double, 3, 3> T;
    //     BoundedMatrix<double, 9, 9> global_size_T;
    //     BoundedVector<double, 9> local_rhs;
    //     noalias(local_rhs) = rRHS;
    //     StructuralMechanicsElementUtilities::BuildRotationMatrixForBeam(T, angle);
    //     StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor2D3NBeam(T, global_size_T);

    //     noalias(rRHS) = prod(global_size_T, local_rhs);
    // }
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement3D2N::RotateAll(
    MatrixType& rLHS,
    VectorType& rRHS,
    const GeometryType& rGeometry
    )
{
    // const double angle = GetAngle();
    // if (std::abs(angle) > std::numeric_limits<double>::epsilon()) {
    //     BoundedMatrix<double, 3, 3> T;
    //     BoundedMatrix<double, 9, 9> global_size_T, aux_product;
    //     BoundedVector<double, 9> local_rhs;
    //     StructuralMechanicsElementUtilities::BuildRotationMatrixForBeam(T, angle);
    //     StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor2D3NBeam(T, global_size_T);

    //     noalias(local_rhs) = rRHS;
    //     noalias(rRHS) = prod(global_size_T, local_rhs);

    //     noalias(aux_product) = prod(rLHS, trans(global_size_T));
    //     noalias(rLHS) = prod(global_size_T, aux_product);
    // }
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement3D2N::save(
    Serializer& rSerializer
    ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement3D2N::load(
    Serializer& rSerializer
    )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
}

/***********************************************************************************/
/***********************************************************************************/

} // Namespace Kratos
