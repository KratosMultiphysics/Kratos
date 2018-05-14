// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Vicente Mataix Ferrándiz
//

// System includes

// External includes


// Project includes
#include "custom_elements/total_lagrangian.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

TotalLagrangian::TotalLagrangian(IndexType NewId, GeometryType::Pointer pGeometry)
    : BaseSolidElement(NewId, pGeometry)
{
    // DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

TotalLagrangian::TotalLagrangian( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
        : BaseSolidElement( NewId, pGeometry, pProperties )
{
    //DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer TotalLagrangian::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_shared<TotalLagrangian>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

void TotalLagrangian::CalculateAll( 
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag 
    )
{
    KRATOS_TRY;
    const auto& r_geom = GetGeometry();
    const unsigned int mat_size = r_geom.PointsNumber() * r_geom.WorkingSpaceDimension();
    if (CalculateStiffnessMatrixFlag == true)
    {
        if (rLeftHandSideMatrix.size1() != mat_size || rLeftHandSideMatrix.size2() != mat_size)
            rLeftHandSideMatrix.resize(mat_size, mat_size, false);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);
    }
    if (CalculateResidualVectorFlag == true)
    {
        if (rRightHandSideVector.size() != mat_size)
            rRightHandSideVector.resize(mat_size, false);
        noalias(rRightHandSideVector) = ZeroVector(mat_size);
    }

    LargeDisplacementDeformationVariables deformation_vars(r_geom, IsAxisymmetric());
    ConstitutiveVariables cl_vars(GetStrainSize());
    for (unsigned int g = 0; g < r_geom.IntegrationPointsNumber(); ++g)
    {
        double weight = GetIntegrationWeight(r_geom.IntegrationPoints(), g, deformation_vars.DetJ0(g)); 
        CalculateStressAndConstitutiveMatrix(deformation_vars, g, cl_vars, rCurrentProcessInfo);
        if (CalculateStiffnessMatrixFlag == true)
            CalculateStiffnessMatrix(deformation_vars, cl_vars, g, weight, rLeftHandSideMatrix);
        if (CalculateResidualVectorFlag == true)
            CalculateAndAddResidualVector(
                rRightHandSideVector, row(r_geom.ShapeFunctionsValues(), g),
                deformation_vars.B(g), rCurrentProcessInfo,
                GetBodyForce(r_geom.IntegrationPoints(), g), cl_vars.StressVector, weight);
    }
    KRATOS_CATCH("")
}

// This function is implemented for compatibility with base solid element.
void TotalLagrangian::CalculateKinematicVariables(KinematicVariables& rThisKinematicVariables,
                                                  const unsigned int PointNumber,
                                                  const GeometryType::IntegrationMethod& rIntegrationMethod)
{
    KRATOS_TRY;
    const auto& r_geom = GetGeometry();
    #ifdef KRATOS_DEBUG
    KRATOS_ERROR_IF(rIntegrationMethod != r_geom.GetDefaultIntegrationMethod())
        << "Non-default integration method.\n";
    #endif
    noalias(rThisKinematicVariables.N) = row(r_geom.ShapeFunctionsValues(), PointNumber);
    LargeDisplacementDeformationVariables deformation_vars(r_geom, IsAxisymmetric());
    rThisKinematicVariables.detJ0 = deformation_vars.DetJ0(PointNumber);
    noalias(rThisKinematicVariables.F) = deformation_vars.F(PointNumber);
    noalias(rThisKinematicVariables.B) = deformation_vars.B(PointNumber);
    rThisKinematicVariables.detF = deformation_vars.DetF(PointNumber);
    KRATOS_CATCH("");
}

double TotalLagrangian::GetIntegrationWeight(const GeometryType::IntegrationPointsArrayType& rThisIntegrationPoints,
                                             const unsigned int PointNumber,
                                             const double detJ)
{
    KRATOS_TRY;
    double weight = BaseSolidElement::GetIntegrationWeight(
        rThisIntegrationPoints, PointNumber, detJ);
    if (GetGeometry().WorkingSpaceDimension() == 2 && GetProperties().Has(THICKNESS))
        weight *= GetProperties()[THICKNESS];
    return weight;
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void TotalLagrangian::CalculateStressAndConstitutiveMatrix(
    LargeDisplacementDeformationVariables& rDeformationVars,
    std::size_t IntegrationIndex,
    ConstitutiveVariables& rOutput,
    ProcessInfo const& rCurrentProcessInfo)
{
    KRATOS_TRY;
    ConstitutiveLaw::Parameters cl_params(GetGeometry(), GetProperties(), rCurrentProcessInfo);
    cl_params.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS |
                               ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR |
                               ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    cl_params.SetDeformationGradientF(rDeformationVars.F(IntegrationIndex));
    cl_params.SetDeterminantF(rDeformationVars.DetF(IntegrationIndex));
    cl_params.SetStrainVector(rDeformationVars.StrainVector(IntegrationIndex));
    cl_params.SetStressVector(rOutput.StressVector);
    cl_params.SetConstitutiveMatrix(rOutput.D);
    mConstitutiveLawVector[IntegrationIndex]->CalculateMaterialResponse(
        cl_params, GetStressMeasure());
    KRATOS_CATCH("");
}

void TotalLagrangian::CalculateStress(Vector& rStrain,
                                      std::size_t IntegrationIndex,
                                      Vector& rStress,
                                      ProcessInfo const& rCurrentProcessInfo)
{
    KRATOS_TRY;
    if (rStress.size() != rStrain.size())
        rStress.resize(rStrain.size(), false);
    ConstitutiveLaw::Parameters cl_params(GetGeometry(), GetProperties(), rCurrentProcessInfo);
    cl_params.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS | ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    cl_params.SetStrainVector(rStrain);
    cl_params.SetStressVector(rStress);
    mConstitutiveLawVector[IntegrationIndex]->CalculateMaterialResponse(cl_params, GetStressMeasure());
    KRATOS_CATCH("");
}

void TotalLagrangian::CalculateStiffnessMatrix(LargeDisplacementDeformationVariables& rDeformationVars,
                                               ConstitutiveVariables& rConstitutiveVars,
                                               std::size_t IntegrationIndex,
                                               double Weight,
                                               MatrixType& rStiffnessMatrix)
{
    KRATOS_TRY;
    /* Material stiffness matrix */
    this->CalculateAndAddKm(
        rStiffnessMatrix, rDeformationVars.B(IntegrationIndex), rConstitutiveVars.D, Weight);
    /* Geometric stiffness matrix */
    this->CalculateAndAddKg(rStiffnessMatrix, rDeformationVars.DN_DX0(IntegrationIndex),
                            rConstitutiveVars.StressVector, Weight);
    KRATOS_CATCH("")
}

std::size_t TotalLagrangian::GetStrainSize() const
{
    return GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();
}

bool TotalLagrangian::IsAxisymmetric() const
{
    return (GetStrainSize() == 4);
}

void TotalLagrangian::CalculateInternalForceSensitivityContribution(
    Vector& rResidualSensitivity,
    std::size_t IntegrationIndex,
    ShapeParameter Deriv,
    LargeDisplacementDeformationVariables& rDeformationVars,
    LargeDisplacementSensitivityVariables& rSensitivityVars,
    Vector const& rStressVector,
    Vector const& rStressSensitivityVector)
{
    KRATOS_TRY;
    const auto& r_geom = GetGeometry();
    const double weight = r_geom.IntegrationPoints()[IntegrationIndex].Weight() *
                          rDeformationVars.DetJ0(IntegrationIndex);
    const double weight_deriv = r_geom.IntegrationPoints()[IntegrationIndex].Weight() *
                                rSensitivityVars.DetJ0(IntegrationIndex, Deriv);
    const Matrix& rB = rDeformationVars.B(IntegrationIndex);
    const Matrix& rB_deriv = rSensitivityVars.B(IntegrationIndex, Deriv);
    if (rResidualSensitivity.size() != rB.size2())
        rResidualSensitivity.resize(rB.size2(), false);
    noalias(rResidualSensitivity) = -weight_deriv * prod(trans(rB), rStressVector);
    noalias(rResidualSensitivity) -= weight * prod(trans(rB_deriv), rStressVector);
    noalias(rResidualSensitivity) -= weight * prod(trans(rB), rStressSensitivityVector);
    KRATOS_CATCH("");
}

void TotalLagrangian::CalculateAndAddExternalForceSensitivityContribution(
    Vector& rResidualSensitivity,
    std::size_t IntegrationIndex,
    ShapeParameter Deriv,
    LargeDisplacementSensitivityVariables& rSensitivityVars,
    Vector const& rN,
    Vector const& rBodyForce,
    ProcessInfo const& rCurrentProcessInfo)
{
    KRATOS_TRY;
    const auto& r_geom = GetGeometry();
    const double weight_deriv = r_geom.IntegrationPoints()[IntegrationIndex].Weight() *
                                rSensitivityVars.DetJ0(IntegrationIndex, Deriv);
    CalculateAndAddExtForceContribution(rN, rCurrentProcessInfo, rBodyForce,
                                        rResidualSensitivity, weight_deriv);
    KRATOS_CATCH("");
}
/***********************************************************************************/
/***********************************************************************************/

int TotalLagrangian::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    int ier = BaseSolidElement::Check(rCurrentProcessInfo);

    return ier;

    KRATOS_CATCH( "" );
}

void TotalLagrangian::CalculateSensitivityMatrix(const Variable<array_1d<double, 3>>& rDesignVariable,
                                                 Matrix& rOutput,
                                                 const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const auto& r_geom = GetGeometry();
    if (rDesignVariable == SHAPE_SENSITIVITY)
    {
        const std::size_t ws_dim = r_geom.WorkingSpaceDimension();
        const std::size_t nnodes = r_geom.PointsNumber();
        Vector residual_deriv(nnodes * ws_dim), N, body_force, stress_vector, stress_vector_deriv;
        if (rOutput.size1() != residual_deriv.size() || rOutput.size2() != residual_deriv.size())
            rOutput.resize(residual_deriv.size(), residual_deriv.size(), false);
        LargeDisplacementDeformationVariables deformation_vars(r_geom, IsAxisymmetric());
        LargeDisplacementSensitivityVariables sensitivity_vars(r_geom);
        for (std::size_t g = 0; g < r_geom.IntegrationPointsNumber(); ++g)
        {
            CalculateStress(deformation_vars.StrainVector(g), g, stress_vector,
                            rCurrentProcessInfo);
            N = row(r_geom.ShapeFunctionsValues(), g);
            body_force = GetBodyForce(r_geom.IntegrationPoints(), g);
            for (auto s = ShapeParameter::Sequence(nnodes, ws_dim); s; ++s)
            {
                const auto& deriv = s.CurrentValue();
                // Assumes constant constitutive matrix wrt design parameter.
                CalculateStress(sensitivity_vars.StrainVector(g, deriv), g,
                                stress_vector_deriv, rCurrentProcessInfo);
                CalculateInternalForceSensitivityContribution(
                    residual_deriv, g, deriv, deformation_vars,
                    sensitivity_vars, stress_vector, stress_vector_deriv);
                CalculateAndAddExternalForceSensitivityContribution(
                    residual_deriv, g, deriv, sensitivity_vars, N, body_force,
                    rCurrentProcessInfo);
                for (std::size_t k = 0; k < residual_deriv.size(); ++k)
                    rOutput(k, deriv.NodeIndex * ws_dim + deriv.Direction) =
                        residual_deriv(k);
            }
        }
    }
    else
        KRATOS_ERROR << "Unsupported variable: " << rDesignVariable << std::endl;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void TotalLagrangian::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseSolidElement );
}

/***********************************************************************************/
/***********************************************************************************/

void TotalLagrangian::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseSolidElement );
}

} // Namespace Kratos
