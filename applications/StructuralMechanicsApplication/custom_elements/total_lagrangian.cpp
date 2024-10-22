// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "custom_elements/total_lagrangian.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/structural_mechanics_math_utilities.hpp"
#include "utilities/geometrical_sensitivity_utility.h"

namespace Kratos
{
namespace
{
class LargeDisplacementKinematics
{
public:
    LargeDisplacementKinematics(Element::GeometryType const& rGeom) : mrGeom(rGeom) { mThisIntegrationMethod = rGeom.GetDefaultIntegrationMethod();}
    LargeDisplacementKinematics(Element::GeometryType const& rGeom, Element::GeometryType::IntegrationMethod ThisIntegrationMethod) : mrGeom(rGeom), mThisIntegrationMethod(ThisIntegrationMethod) {}

    void DeformationGradient(std::size_t IntegrationPoint, Matrix& rF)
    {
        if (IntegrationPoint != mCurrentIntegrationPoint)
            Recalculate(IntegrationPoint);
        GeometryUtils::DeformationGradient(mJ, mInvJ0, rF);
    }

    void ShapeFunctionsGradients(std::size_t IntegrationPoint, Matrix& rDN_DX0)
    {
        if (IntegrationPoint != mCurrentIntegrationPoint)
            Recalculate(IntegrationPoint);
        Matrix const& rDN_De = mrGeom.ShapeFunctionsLocalGradients(mThisIntegrationMethod)[IntegrationPoint];
        GeometryUtils::ShapeFunctionsGradients(rDN_De, mInvJ0, rDN_DX0);
    }

    double DetJ0(std::size_t IntegrationPoint)
    {
        if (IntegrationPoint != mCurrentIntegrationPoint)
            Recalculate(IntegrationPoint);
        return mDetJ0;
    }

private:
    void Recalculate(std::size_t IntegrationPoint)
    {
        // Calculate mInvJ0 and mDetJ0.
        GeometryUtils::JacobianOnInitialConfiguration(
            mrGeom, mrGeom.IntegrationPoints(mThisIntegrationMethod)[IntegrationPoint], mJ);
        MathUtils<double>::InvertMatrix(mJ, mInvJ0, mDetJ0);
        // Calculate mJ.
        mrGeom.Jacobian(mJ, IntegrationPoint, mThisIntegrationMethod);
        // Update current integration point.
        mCurrentIntegrationPoint = IntegrationPoint;
    }

    Element::GeometryType const& mrGeom; /// Currently geometry
    Element::GeometryType::IntegrationMethod mThisIntegrationMethod; /// Currently selected integration methods
    Matrix mJ;
    Matrix mInvJ0;
    double mDetJ0 = 0.0;
    std::size_t mCurrentIntegrationPoint = -1;
};

void CalculateGreenLagrangeStrainSensitivity(Matrix const& rF,
                                             Matrix const& rF_Deriv,
                                             Matrix& rE_Deriv)
{
    rE_Deriv = 0.5 * (prod(trans(rF_Deriv), rF) + prod(trans(rF), rF_Deriv));
}
}

/***********************************************************************************/
/***********************************************************************************/

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
    return Kratos::make_intrusive<TotalLagrangian>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

//************************************************************************************
//************************************************************************************

Element::Pointer TotalLagrangian::Create( IndexType NewId,  GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<TotalLagrangian>( NewId, pGeom, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

TotalLagrangian::~TotalLagrangian()
{
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer TotalLagrangian::Clone (
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    TotalLagrangian::Pointer p_new_elem = Kratos::make_intrusive<TotalLagrangian>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    // Currently selected integration methods
    p_new_elem->SetIntegrationMethod(BaseType::mThisIntegrationMethod);

    // The vector containing the constitutive laws
    p_new_elem->SetConstitutiveLawVector(BaseType::mConstitutiveLawVector);

    return p_new_elem;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void TotalLagrangian::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag
    )
{
    KRATOS_TRY;

    const SizeType number_of_nodes = this->GetGeometry().size();
    const SizeType dimension = this->GetGeometry().WorkingSpaceDimension();
    const auto strain_size = GetStrainSize();
    const bool is_rotated = IsElementRotated();

    KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
    ConstitutiveVariables this_constitutive_variables(strain_size);

    // Resizing as needed the LHS
    const SizeType mat_size = number_of_nodes * dimension;

    if ( CalculateStiffnessMatrixFlag == true ) { // Calculation of the matrix is required
        if ( rLeftHandSideMatrix.size1() != mat_size )
            rLeftHandSideMatrix.resize( mat_size, mat_size, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( mat_size, mat_size ); //resetting LHS
    }

    // Resizing as needed the RHS
    if ( CalculateResidualVectorFlag == true ) { // Calculation of the matrix is required
        if ( rRightHandSideVector.size() != mat_size )
            rRightHandSideVector.resize( mat_size, false );

        noalias(rRightHandSideVector) = ZeroVector( mat_size ); //resetting RHS
    }

    // Reading integration points
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    // Set constitutive law flags:
    Flags& ConstitutiveLawOptions=Values.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    if ( CalculateStiffnessMatrixFlag ) {
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    } else {
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
    }

    // If strain has to be computed inside of the constitutive law with PK2
    Values.SetStrainVector(this_constitutive_variables.StrainVector); //this is the input  parameter

    // Some declarations
    array_1d<double, 3> body_force;
    double int_to_reference_weight;

    // Computing in all integrations points
    for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
        // Contribution to external forces
        noalias(body_force) = this->GetBodyForce(integration_points, point_number);

        // Compute element kinematics B, F, DN_DX ...
        this->CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());

        // Compute material response
        this->CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points, this->GetStressMeasure(), is_rotated);

        // Calculating weights for integration on the reference configuration
        int_to_reference_weight = GetIntegrationWeight(integration_points, point_number, this_kinematic_variables.detJ0);

        if ( dimension == 2 && this->GetProperties().Has( THICKNESS ))
            int_to_reference_weight *= this->GetProperties()[THICKNESS];

        if ( CalculateStiffnessMatrixFlag == true ) { // Calculation of the matrix is required
            // Contributions to stiffness matrix calculated on the reference config
            /* Material stiffness matrix */
            this->CalculateAndAddKm( rLeftHandSideMatrix, this_kinematic_variables.B, this_constitutive_variables.D, int_to_reference_weight );

            /* Geometric stiffness matrix */
            this->CalculateAndAddKg( rLeftHandSideMatrix, this_kinematic_variables.DN_DX, this_constitutive_variables.StressVector, int_to_reference_weight );
        }

        if ( CalculateResidualVectorFlag == true ) { // Calculation of the matrix is required
            this->CalculateAndAddResidualVector(rRightHandSideVector, this_kinematic_variables, rCurrentProcessInfo, body_force, this_constitutive_variables.StressVector, int_to_reference_weight);
        }
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void TotalLagrangian::CalculateKinematicVariables(
    KinematicVariables& rThisKinematicVariables,
    const IndexType PointNumber,
    const GeometryType::IntegrationMethod& rIntegrationMethod
    )
{
    // Shape functions
    rThisKinematicVariables.N = row(GetGeometry().ShapeFunctionsValues(rIntegrationMethod), PointNumber);

    rThisKinematicVariables.detJ0 = this->CalculateDerivativesOnReferenceConfiguration(rThisKinematicVariables.J0, rThisKinematicVariables.InvJ0, rThisKinematicVariables.DN_DX, PointNumber, rIntegrationMethod);
    KRATOS_ERROR_IF(rThisKinematicVariables.detJ0 < 0.0) << "WARNING:: ELEMENT ID: " << this->Id() << " INVERTED. DETJ0: " << rThisKinematicVariables.detJ0 << std::endl;

    Matrix J;
    J = this->GetGeometry().Jacobian(J, PointNumber, rIntegrationMethod);
    if (IsAxissymmetric()) {
        CalculateAxisymmetricF(J, rThisKinematicVariables.InvJ0,
                               rThisKinematicVariables.N,
                               rThisKinematicVariables.F);
        CalculateAxisymmetricB(
            rThisKinematicVariables.B, rThisKinematicVariables.F,
            rThisKinematicVariables.DN_DX, rThisKinematicVariables.N);
    } else {
        GeometryUtils::DeformationGradient(J, rThisKinematicVariables.InvJ0,
                                           rThisKinematicVariables.F);
        CalculateB(rThisKinematicVariables.B, rThisKinematicVariables.F,
                   rThisKinematicVariables.DN_DX);
    }

    rThisKinematicVariables.detF = MathUtils<double>::Det(rThisKinematicVariables.F);
}

/***********************************************************************************/
/***********************************************************************************/

void TotalLagrangian::CalculateB(Matrix& rB, Matrix const& rF, const Matrix& rDN_DX)
{
    KRATOS_TRY;

    if (GetGeometry().WorkingSpaceDimension() == 2)
        Calculate2DB(rB, rF, rDN_DX);
    else
        Calculate3DB(rB, rF, rDN_DX);

    KRATOS_CATCH("");
}

void TotalLagrangian::Calculate2DB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    KRATOS_TRY

    const SizeType number_of_nodes = this->GetGeometry().PointsNumber();
    const SizeType dimension = this->GetGeometry().WorkingSpaceDimension();

    for (IndexType i = 0; i < number_of_nodes; ++i)
    {
        const auto index = dimension * i;
        rB(0, index + 0) = rF(0, 0) * rDN_DX(i, 0);
        rB(0, index + 1) = rF(1, 0) * rDN_DX(i, 0);
        rB(1, index + 0) = rF(0, 1) * rDN_DX(i, 1);
        rB(1, index + 1) = rF(1, 1) * rDN_DX(i, 1);
        rB(2, index + 0) = rF(0, 0) * rDN_DX(i, 1) + rF(0, 1) * rDN_DX(i, 0);
        rB(2, index + 1) = rF(1, 0) * rDN_DX(i, 1) + rF(1, 1) * rDN_DX(i, 0);
    }

    KRATOS_CATCH( "" )
}

void TotalLagrangian::Calculate3DB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    KRATOS_TRY

    const SizeType number_of_nodes = this->GetGeometry().PointsNumber();
    const SizeType dimension = this->GetGeometry().WorkingSpaceDimension();

    for (IndexType i = 0; i < number_of_nodes; ++i)
    {
        const auto index = dimension * i;
        rB(0, index + 0) = rF(0, 0) * rDN_DX(i, 0);
        rB(0, index + 1) = rF(1, 0) * rDN_DX(i, 0);
        rB(0, index + 2) = rF(2, 0) * rDN_DX(i, 0);
        rB(1, index + 0) = rF(0, 1) * rDN_DX(i, 1);
        rB(1, index + 1) = rF(1, 1) * rDN_DX(i, 1);
        rB(1, index + 2) = rF(2, 1) * rDN_DX(i, 1);
        rB(2, index + 0) = rF(0, 2) * rDN_DX(i, 2);
        rB(2, index + 1) = rF(1, 2) * rDN_DX(i, 2);
        rB(2, index + 2) = rF(2, 2) * rDN_DX(i, 2);
        rB(3, index + 0) = rF(0, 0) * rDN_DX(i, 1) + rF(0, 1) * rDN_DX(i, 0);
        rB(3, index + 1) = rF(1, 0) * rDN_DX(i, 1) + rF(1, 1) * rDN_DX(i, 0);
        rB(3, index + 2) = rF(2, 0) * rDN_DX(i, 1) + rF(2, 1) * rDN_DX(i, 0);
        rB(4, index + 0) = rF(0, 1) * rDN_DX(i, 2) + rF(0, 2) * rDN_DX(i, 1);
        rB(4, index + 1) = rF(1, 1) * rDN_DX(i, 2) + rF(1, 2) * rDN_DX(i, 1);
        rB(4, index + 2) = rF(2, 1) * rDN_DX(i, 2) + rF(2, 2) * rDN_DX(i, 1);
        rB(5, index + 0) = rF(0, 2) * rDN_DX(i, 0) + rF(0, 0) * rDN_DX(i, 2);
        rB(5, index + 1) = rF(1, 2) * rDN_DX(i, 0) + rF(1, 0) * rDN_DX(i, 2);
        rB(5, index + 2) = rF(2, 2) * rDN_DX(i, 0) + rF(2, 0) * rDN_DX(i, 2);
    }

    KRATOS_CATCH("")
}

void TotalLagrangian::CalculateAxisymmetricB(Matrix& rB,
                                             const Matrix& rF,
                                             const Matrix& rDN_DX,
                                             const Vector& rN)
{
    KRATOS_TRY

    const SizeType number_of_nodes = this->GetGeometry().PointsNumber();
    const SizeType dimension = this->GetGeometry().WorkingSpaceDimension();
    double radius = 0.0;
    radius = StructuralMechanicsMathUtilities::CalculateRadius(rN, this->GetGeometry());
    for (IndexType i = 0; i < number_of_nodes; ++i)
    {
        const SizeType index = dimension * i;

        rB(0, index + 0) = rF(0, 0) * rDN_DX(i, 0);
        rB(0, index + 1) = rF(1, 0) * rDN_DX(i, 0);
        rB(1, index + 1) = rF(0, 1) * rDN_DX(i, 1);
        rB(1, index + 1) = rF(1, 1) * rDN_DX(i, 1);
        rB(2, index + 0) = rN[i] / radius;
        rB(3, index + 0) = rF(0, 0) * rDN_DX(i, 1) + rF(0, 1) * rDN_DX(i, 0);
        rB(3, index + 1) = rF(1, 0) * rDN_DX(i, 1) + rF(1, 1) * rDN_DX(i, 0);
    }

    KRATOS_CATCH("")
}

void TotalLagrangian::CalculateAxisymmetricF(Matrix const& rJ,
                                             Matrix const& rInvJ0,
                                             Vector const& rN,
                                             Matrix& rF)
{
    KRATOS_TRY;

    GeometryUtils::DeformationGradient(rJ, rInvJ0, rF);
    BoundedMatrix<double, 2, 2> F2x2 = rF;
    rF.resize(3, 3, false);
    for (unsigned i = 0; i < 2; ++i)
    {
        for (unsigned j = 0; j < 2; ++j)
            rF(i, j) = F2x2(i, j);
        rF(i, 2) = rF(2, i) = 0.0;
    }
    const double current_radius = StructuralMechanicsMathUtilities::CalculateRadius(
        rN, this->GetGeometry(), Current);
    const double initial_radius = StructuralMechanicsMathUtilities::CalculateRadius(
        rN, this->GetGeometry(), Initial);
    rF(2, 2) = current_radius / initial_radius;

    KRATOS_CATCH("");
}

void TotalLagrangian::CalculateStress(Vector& rStrain,
                                      std::size_t IntegrationPoint,
                                      Vector& rStress,
                                      ProcessInfo const& rCurrentProcessInfo)
{
    KRATOS_TRY;
    ConstitutiveLaw::Parameters cl_params(GetGeometry(), GetProperties(), rCurrentProcessInfo);
    cl_params.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS | ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    cl_params.SetStrainVector(rStrain);
    cl_params.SetStressVector(rStress);
    mConstitutiveLawVector[IntegrationPoint]->CalculateMaterialResponse(cl_params, GetStressMeasure());
    KRATOS_CATCH("");
}

void TotalLagrangian::CalculateStress(Matrix const& rF,
                                      std::size_t IntegrationPoint,
                                      Vector& rStress,
                                      ProcessInfo const& rCurrentProcessInfo)
{
    KRATOS_TRY;
    Vector strain(mConstitutiveLawVector[IntegrationPoint]->GetStrainSize());
    CalculateStrain(rF, *mConstitutiveLawVector[IntegrationPoint], strain, rCurrentProcessInfo);
    ConstitutiveLaw::Parameters cl_params(GetGeometry(), GetProperties(), rCurrentProcessInfo);
    cl_params.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS | ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    cl_params.SetStrainVector(strain);
    cl_params.SetStressVector(rStress);
    mConstitutiveLawVector[IntegrationPoint]->CalculateMaterialResponse(cl_params, GetStressMeasure());
    KRATOS_CATCH("");
}

void TotalLagrangian::CalculateStrain(Matrix const& rF,
                                      std::size_t IntegrationPoint,
                                      Vector& rStrain,
                                      ProcessInfo const& rCurrentProcessInfo)
{
    KRATOS_TRY;
    ConstitutiveLaw::Parameters cl_params(GetGeometry(), GetProperties(), rCurrentProcessInfo);
    cl_params.SetDeformationGradientF(rF);
    cl_params.SetStrainVector(rStrain);
    mConstitutiveLawVector[IntegrationPoint]->CalculateMaterialResponse(cl_params, GetStressMeasure());
    KRATOS_CATCH("");
}

void TotalLagrangian::CalculateShapeSensitivity(ShapeParameter Deriv,
                                                Matrix& rDN_DX0,
                                                Matrix& rDN_DX0_Deriv,
                                                Matrix& rF_Deriv,
                                                double& rDetJ0_Deriv,
                                                std::size_t IntegrationPointIndex)
{
    KRATOS_TRY;
    const unsigned ws_dim = GetGeometry().WorkingSpaceDimension();
    const unsigned ls_dim = GetGeometry().LocalSpaceDimension();
    Matrix J0(ws_dim, ls_dim);
    GeometryUtils::JacobianOnInitialConfiguration(
        GetGeometry(), GetGeometry().IntegrationPoints(this->GetIntegrationMethod())[IntegrationPointIndex], J0);
    const Matrix& rDN_De = GetGeometry().ShapeFunctionsLocalGradients()[IntegrationPointIndex];
    auto sensitivity_utility =
        GeometricalSensitivityUtility(J0, rDN_De);
    sensitivity_utility.CalculateSensitivity(Deriv, rDetJ0_Deriv, rDN_DX0_Deriv);
    rF_Deriv.resize(ws_dim, ws_dim, false);
    rF_Deriv.clear();
    for (unsigned i = 0; i < ws_dim; ++i)
        for (unsigned j = 0; j < ws_dim; ++j)
        {
            for (unsigned k = 0; k < GetGeometry().PointsNumber(); ++k)
                rF_Deriv(i, j) +=
                    GetGeometry()[k].Coordinates()[i] * rDN_DX0_Deriv(k, j);
        }
    for (unsigned j = 0; j < ws_dim; ++j)
      rF_Deriv(Deriv.Direction, j) += rDN_DX0(Deriv.NodeIndex, j);
    KRATOS_CATCH("");
}

void TotalLagrangian::CalculateBSensitivity(Matrix const& rDN_DX,
                                            Matrix const& rF,
                                            Matrix const& rDN_DX_Deriv,
                                            Matrix const& rF_Deriv,
                                            Matrix& rB_Deriv)
{
    KRATOS_TRY;
    const unsigned dimension = GetGeometry().WorkingSpaceDimension();
    const auto strain_size = GetStrainSize();
    rB_Deriv.resize(strain_size, dimension * GetGeometry().PointsNumber(), false);
    Matrix B_deriv1(rB_Deriv.size1(), rB_Deriv.size2());
    Matrix B_deriv2(rB_Deriv.size1(), rB_Deriv.size2());

    CalculateB(B_deriv1, rF_Deriv, rDN_DX);
    CalculateB(B_deriv2, rF, rDN_DX_Deriv);
    rB_Deriv = B_deriv1 + B_deriv2;
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

std::size_t TotalLagrangian::GetStrainSize() const
{
    return GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();
}

/***********************************************************************************/
/***********************************************************************************/


bool TotalLagrangian::IsAxissymmetric() const
{
    return (GetStrainSize() == 4);
}
/***********************************************************************************/
/***********************************************************************************/

void TotalLagrangian::CalculateSensitivityMatrix(
    const Variable<array_1d<double, 3>>& rDesignVariable,
    Matrix& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    const auto& r_geom = GetGeometry();
    if (rDesignVariable == SHAPE_SENSITIVITY) {
        const std::size_t ws_dim = r_geom.WorkingSpaceDimension();
        const std::size_t nnodes = r_geom.PointsNumber();
        const std::size_t mat_dim = nnodes * ws_dim;
        rOutput.resize(mat_dim, mat_dim, false);
        rOutput.clear();
        Matrix F, F_deriv, DN_DX0_deriv, strain_tensor_deriv, DN_DX0, B, B_deriv;
        Matrix M_deriv;
        const auto strain_size = GetStrainSize();
        B.resize(strain_size, ws_dim * nnodes, false);
        Vector strain_vector_deriv(strain_size);
        Vector stress_vector(strain_size), stress_vector_deriv(strain_size);
        Vector residual_deriv(ws_dim * nnodes);
        Vector body_force, acceleration;
        double detJ0_deriv;
        LargeDisplacementKinematics large_disp_kinematics(r_geom, this->GetIntegrationMethod());
        for (std::size_t g = 0; g < r_geom.IntegrationPointsNumber(); ++g) {
            large_disp_kinematics.DeformationGradient(g, F);
            large_disp_kinematics.ShapeFunctionsGradients(g, DN_DX0);
            CalculateB(B, F, DN_DX0);
            CalculateStress(F, g, stress_vector, rCurrentProcessInfo);
            double weight = GetIntegrationWeight(
                r_geom.IntegrationPoints(this->GetIntegrationMethod()), g, large_disp_kinematics.DetJ0(g));
            const Vector& rN = row(GetGeometry().ShapeFunctionsValues(), g);
            body_force = GetBodyForce(r_geom.IntegrationPoints(this->GetIntegrationMethod()), g);

            for (auto s = ShapeParameter::Sequence(nnodes, ws_dim); s; ++s) {
                const auto& deriv = s.CurrentValue();
                CalculateShapeSensitivity(deriv, DN_DX0, DN_DX0_deriv, F_deriv, detJ0_deriv, g);
                CalculateGreenLagrangeStrainSensitivity(F, F_deriv, strain_tensor_deriv);
                noalias(strain_vector_deriv) =
                    MathUtils<double>::StrainTensorToVector(strain_tensor_deriv);
                CalculateStress(strain_vector_deriv, g, stress_vector_deriv, rCurrentProcessInfo);
                CalculateBSensitivity(DN_DX0, F, DN_DX0_deriv, F_deriv, B_deriv);
                const double weight_deriv =
                    GetIntegrationWeight(r_geom.IntegrationPoints(this->GetIntegrationMethod()), g, detJ0_deriv);
                noalias(residual_deriv) = -weight_deriv * prod(trans(B), stress_vector);
                noalias(residual_deriv) -= weight * prod(trans(B_deriv), stress_vector);
                noalias(residual_deriv) -= weight * prod(trans(B), stress_vector_deriv);
                CalculateAndAddExtForceContribution(
                    rN, rCurrentProcessInfo, body_force, residual_deriv, weight_deriv);
                for (std::size_t k = 0; k < residual_deriv.size(); ++k)
                    rOutput(deriv.NodeIndex * ws_dim + deriv.Direction, k) +=
                        residual_deriv(k);
            }
        }

        for (auto s = ShapeParameter::Sequence(nnodes, ws_dim); s; ++s) {
            const auto& deriv = s.CurrentValue();
            CalculateShapeGradientOfMassMatrix(M_deriv, deriv);
            GetSecondDerivativesVector(acceleration);
            noalias(residual_deriv) = -prod(M_deriv, acceleration);
            for (std::size_t k = 0; k < residual_deriv.size(); ++k)
                rOutput(deriv.NodeIndex * ws_dim + deriv.Direction, k) +=
                    residual_deriv(k);
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


