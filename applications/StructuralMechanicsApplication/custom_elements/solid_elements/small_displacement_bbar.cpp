// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Manuel Caicedo
//                   Marcelo Raschi
//                   Javier Mroginski


// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_elements/small_displacement_bbar.h"

namespace Kratos
{
SmallDisplacementBbar::SmallDisplacementBbar(IndexType NewId,
                                             GeometryType::Pointer pGeometry)
        : BaseSolidElement(NewId, pGeometry)
{
    // DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

SmallDisplacementBbar::SmallDisplacementBbar(IndexType NewId,
                                             GeometryType::Pointer pGeometry,
                                             PropertiesType::Pointer pProperties)
        : BaseSolidElement(NewId, pGeometry, pProperties)
{
}

//************************************************************************************
//************************************************************************************

Element::Pointer SmallDisplacementBbar::Create(IndexType NewId,
                                               NodesArrayType const& ThisNodes,
                                               PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<SmallDisplacementBbar>(
            NewId, GetGeometry().Create(ThisNodes), pProperties);
}

//************************************************************************************
//************************************************************************************

Element::Pointer SmallDisplacementBbar::Create(IndexType NewId,
                                               GeometryType::Pointer pGeom,
                                               PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<SmallDisplacementBbar>(
            NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

//************************************************************************************
//************************************************************************************

SmallDisplacementBbar::~SmallDisplacementBbar()
{
}

//************************************************************************************
//************************************************************************************

bool SmallDisplacementBbar::UseElementProvidedStrain() const
{
    return true;
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementBbar::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag
    )
{
    KRATOS_TRY;

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType strain_size = (dimension == 3) ? 6 : 4; // necessary include component zz in the computation of kinematic variables
    const bool is_rotated = IsElementRotated();

    KinematicVariablesBbar this_kinematic_variables(strain_size, dimension, number_of_nodes);
    ConstitutiveVariables this_constitutive_variables(strain_size);

    // Resizing as needed the LHS
    const SizeType mat_size = number_of_nodes * dimension;

    if ( CalculateStiffnessMatrixFlag ) { // Calculation of the matrix is required
        if ( rLeftHandSideMatrix.size1() != mat_size )
            rLeftHandSideMatrix.resize( mat_size, mat_size, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( mat_size, mat_size ); //resetting LHS
    }

    // Resizing as needed the RHS
    if ( CalculateResidualVectorFlag ) { // Calculation of the matrix is required
        if ( rRightHandSideVector.size() != mat_size )
            rRightHandSideVector.resize( mat_size, false );

        noalias(rRightHandSideVector) = ZeroVector( mat_size ); //resetting RHS
    }

    // Reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points =
        GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

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
    Values.SetStrainVector(this_constitutive_variables.StrainVector); //this is the input parameter

    // compute Hydrostatic B-Matrix
    this->SmallDisplacementBbar::CalculateHydrostaticDeformationMatrix(this_kinematic_variables);

    // Some declarations
    array_1d<double, 3> body_force;
    double int_to_reference_weight;

    // Computing in all integrations points
    for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
        // Contribution to external forces
        noalias(body_force) = this->GetBodyForce(integration_points, point_number);

        // Compute element kinematics B, F, DN_DX ...
        CalculateKinematicVariablesBbar(this_kinematic_variables, point_number, integration_points);

        // Compute material response
        CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables,
                                       Values, point_number, integration_points,
                                       GetStressMeasure(), is_rotated);

        // Calculating weights for integration on the reference configuration
        int_to_reference_weight = GetIntegrationWeight(integration_points, point_number,
                                                              this_kinematic_variables.detJ0);

        if ( dimension == 2 && GetProperties().Has( THICKNESS ))
            int_to_reference_weight *= GetProperties()[THICKNESS];

        if ( CalculateStiffnessMatrixFlag ) { //calculation of the matrix is required
            // Contributions to stiffness matrix calculated on the reference config
            this->CalculateAndAddKm( rLeftHandSideMatrix, this_kinematic_variables.B,
                                     this_constitutive_variables.D, int_to_reference_weight );
        }

        if ( CalculateResidualVectorFlag ) { //calculation of the matrix is required
            this->CalculateAndAddResidualVector(rRightHandSideVector, this_kinematic_variables,
                                                rCurrentProcessInfo, body_force,
                                                this_constitutive_variables.StressVector, int_to_reference_weight);
        }
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementBbar::CalculateKinematicVariables(
        KinematicVariables& rThisKinematicVariables,
        const IndexType PointNumber,
        const GeometryType::IntegrationMethod& rIntegrationMethod
)
{
    const GeometryType::IntegrationPointsArrayType& r_integration_points =
            GetGeometry().IntegrationPoints(rIntegrationMethod);
    // Shape functions
    rThisKinematicVariables.N = GetGeometry().ShapeFunctionsValues(rThisKinematicVariables.N, r_integration_points[PointNumber].Coordinates());

    rThisKinematicVariables.detJ0 = CalculateDerivativesOnReferenceConfiguration(
            rThisKinematicVariables.J0,
            rThisKinematicVariables.InvJ0,
            rThisKinematicVariables.DN_DX,
            PointNumber,
            rIntegrationMethod);

    KRATOS_ERROR_IF(rThisKinematicVariables.detJ0 < 0.0)
            << "WARNING:: ELEMENT ID: " << this->Id()
            << " INVERTED. DETJ0: " << rThisKinematicVariables.detJ0
            << std::endl;

    // Compute B
    CalculateB( rThisKinematicVariables.B,
                rThisKinematicVariables.DN_DX);

    // Compute equivalent F
    Vector displacements;
    GetValuesVector(displacements);
    Vector strain_vector = prod(rThisKinematicVariables.B, displacements);
    rThisKinematicVariables.F = ComputeEquivalentF(strain_vector);
    rThisKinematicVariables.detF = MathUtils<double>::Det(rThisKinematicVariables.F);
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementBbar::CalculateKinematicVariablesBbar(
    KinematicVariablesBbar& rThisKinematicVariables,
    const IndexType PointNumber,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints
    )
{
    // Shape functions
    rThisKinematicVariables.N = GetGeometry().ShapeFunctionsValues(
            rThisKinematicVariables.N,
            IntegrationPoints[PointNumber].Coordinates());

    rThisKinematicVariables.detJ0 = CalculateDerivativesOnReferenceConfiguration(
            rThisKinematicVariables.J0,
            rThisKinematicVariables.InvJ0,
            rThisKinematicVariables.DN_DX,
            PointNumber,
            this->GetIntegrationMethod());

    KRATOS_ERROR_IF(rThisKinematicVariables.detJ0 < 0.0)
                << "WARNING:: ELEMENT ID: " << this->Id()
                << " INVERTED. DETJ0: " << rThisKinematicVariables.detJ0
                << std::endl;

    // Compute Bbar
    CalculateBbar(rThisKinematicVariables.B,
                  rThisKinematicVariables.Bh,
                  rThisKinematicVariables.DN_DX,
                  IntegrationPoints,
                  PointNumber);

    // Compute equivalent F
    Vector displacements;
    GetValuesVector(displacements);
    Vector strain_vector = prod(rThisKinematicVariables.B, displacements);
    rThisKinematicVariables.F = ComputeEquivalentF(strain_vector);
    rThisKinematicVariables.detF = MathUtils<double>::Det(rThisKinematicVariables.F);
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementBbar::SetConstitutiveVariables(
    KinematicVariables& rThisKinematicVariables,
    ConstitutiveVariables& rThisConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType PointNumber,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints
    )
{
    // Displacements vector
    Vector displacements;
    GetValuesVector(displacements);

    // Compute strain
    noalias(rThisConstitutiveVariables.StrainVector) = prod(rThisKinematicVariables.B, displacements);

    // Compute equivalent F
    rThisKinematicVariables.F = ComputeEquivalentF(rThisConstitutiveVariables.StrainVector);

    // Here we essentially set the input parameters
    rThisKinematicVariables.detF = MathUtils<double>::Det(rThisKinematicVariables.F);
    rValues.SetShapeFunctionsValues(rThisKinematicVariables.N); // shape functions
    rValues.SetDeterminantF(rThisKinematicVariables.detF); //assuming the determinant is computed somewhere else
    rValues.SetDeformationGradientF(rThisKinematicVariables.F); //F computed somewhere else

    // Here we set the space on which the results shall be written
    rValues.SetConstitutiveMatrix(rThisConstitutiveVariables.D); //assuming the determinant is computed somewhere else
    rValues.SetStressVector(rThisConstitutiveVariables.StressVector); //F computed somewhere else
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementBbar::CalculateB(
        Matrix& rB,
        const Matrix& DN_DX
        )
{
    KRATOS_TRY;

    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    rB.clear();

    if (dimension == 2) {
        for ( IndexType i = 0; i < number_of_nodes; ++i ) {
            unsigned int index = 2 * i;
            rB(0, index + 0) = DN_DX(i, 0);
            rB(0, index + 1) = 0.0;
            rB(1, index + 0) = 0.0;
            rB(1, index + 1) = DN_DX(i, 1);
            rB(2, index + 0) = 0.0;
            rB(2, index + 1) = 0.0;
            rB(3, index + 0) = DN_DX(i, 1);
            rB(3, index + 1) = DN_DX(i, 0);
        }
    }
    else {
        for ( IndexType i = 0; i < number_of_nodes; ++i ) {
            unsigned int index = 3 * i;
            rB(0, index + 0) = DN_DX(i, 0);
            rB(1, index + 1) = DN_DX(i, 1);
            rB(2, index + 2) = DN_DX(i, 2);
            rB(3, index + 0) = DN_DX(i, 1);
            rB(3, index + 1) = DN_DX(i, 0);
            rB(4, index + 1) = DN_DX(i, 2);
            rB(4, index + 2) = DN_DX(i, 1);
            rB(5, index + 0) = DN_DX(i, 2);
            rB(5, index + 2) = DN_DX(i, 0);
        }
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************
void SmallDisplacementBbar::CalculateBbar(
            Matrix &rB,
            Vector &rBh,
            const Matrix &rDN_DX,
            const GeometryType::IntegrationPointsArrayType &IntegrationPoints,
            const IndexType PointNumber
    )
{
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int strain_size = (dimension == 3) ? 6 : 4; // added component zz, necessary for plasticity.

    if (rB.size1() != strain_size || rB.size2() != dimension * number_of_nodes)
        rB.resize(strain_size, dimension * number_of_nodes, false);

    Matrix rBn;
    rBn.resize(strain_size, dimension * number_of_nodes, false);
    noalias(rBn) = ZeroMatrix(strain_size, dimension * number_of_nodes);

    this->CalculateB(rB, rDN_DX);

    if (dimension == 2) {
        rBn(0, 0) = 2. / 3. * rB(0, 0);
        rBn(0, 1) = -1. / 3. * rB(1, 1);
        rBn(0, 2) = 2. / 3. * rB(0, 2);
        rBn(0, 3) = -1. / 3. * rB(1, 3);
        rBn(0, 4) = 2. / 3. * rB(0, 4);
        rBn(0, 5) = -1. / 3. * rB(1, 5);
        rBn(0, 6) = 2. / 3. * rB(0, 6);
        rBn(0, 7) = -1. / 3. * rB(1, 7);

        rBn(1, 0) = -1. / 3. * rB(0, 0);
        rBn(1, 1) = 2. / 3. * rB(1, 1);
        rBn(1, 2) = -1. / 3. * rB(0, 2);
        rBn(1, 3) = 2. / 3. * rB(1, 3);
        rBn(1, 4) = -1. / 3. * rB(0, 4);
        rBn(1, 5) = 2. / 3. * rB(1, 5);
        rBn(1, 6) = -1. / 3. * rB(0, 6);
        rBn(1, 7) = 2. / 3. * rB(1, 7);

        rBn(2, 0) = -1. / 3. * rB(0, 0);
        rBn(2, 1) = -1. / 3. * rB(1, 1);
        rBn(2, 2) = -1. / 3. * rB(0, 2);
        rBn(2, 3) = -1. / 3. * rB(1, 3);
        rBn(2, 4) = -1. / 3. * rB(0, 4);
        rBn(2, 5) = -1. / 3. * rB(1, 5);
        rBn(2, 6) = -1. / 3. * rB(0, 6);
        rBn(2, 7) = -1. / 3. * rB(1, 7);

        for(IndexType i = 0; i < number_of_nodes * dimension; i++) {
            rBn(0, i) += 1. / 3. * rBh(i);
            rBn(1, i) += 1. / 3. * rBh(i);
            rBn(2, i) += 1. / 3. * rBh(i);
            rBn(3, i) = rB(3, i);
        }
    }
    else {
        for(IndexType i = 0; i < number_of_nodes; i++) {
            unsigned int index = 3 * i;
            rBn(0, index + 0) = 2. / 3. * rB(0, index + 0);
            rBn(1, index + 0) = -1. / 3. * rB(0, index + 0);
            rBn(2, index + 0) = -1. / 3. * rB(0, index + 0);
            rBn(0, index + 1) = -1. / 3. * rB(1, index + 1);
            rBn(1, index + 1) = 2. / 3. * rB(1, index + 1);
            rBn(2, index + 1) = -1. / 3. * rB(1, index + 1);
            rBn(0, index + 2) = -1. / 3. * rB(2, index + 2);
            rBn(1, index + 2) = -1. / 3. * rB(2, index + 2);
            rBn(2, index + 2) = 2. / 3. * rB(2, index + 2);
        }
        for(IndexType i = 0; i < number_of_nodes * dimension; i++) {
            rBn(0, i) += 1. / 3. * rBh(i);
            rBn(1, i) += 1. / 3. * rBh(i);
            rBn(2, i) += 1. / 3. * rBh(i);
            rBn(3, i) = rB(3, i);
            rBn(4, i) = rB(4, i);
            rBn(5, i) = rB(5, i);
        }
    }

    rB = rBn;

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementBbar::CalculateHydrostaticDeformationMatrix(KinematicVariablesBbar& rVariables)
{
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());
    rVariables.Bh.resize(dimension * number_of_nodes, false);
    noalias(rVariables.Bh) = ZeroVector(dimension * number_of_nodes);

    if (dimension == 2) {
        double TotalArea = 0.0;
        for (IndexType PointNumber = 0; PointNumber < integration_points.size(); PointNumber++) {
            // Compute element kinematics B, F, DN_DX ...
            CalculateKinematicVariables(rVariables, PointNumber, this->GetIntegrationMethod());
            double IntegrationWeight = rVariables.detJ0 * integration_points[PointNumber].Weight();
            TotalArea += IntegrationWeight;

            for(IndexType i = 0; i < number_of_nodes * 2; i++) {
                // Bh = Bh + sum(Bs(1:3, : )) * wg(iPG) * detJ;
                rVariables.Bh(i) += (rVariables.B(0, i) + rVariables.B(1, i)) * IntegrationWeight;
            }
        }
        for(IndexType i = 0; i < number_of_nodes * 2; i++) {
            rVariables.Bh(i) /= TotalArea;
        }
    }
    else {
        double TotalVolume = 0.0;
        for (IndexType PointNumber = 0; PointNumber < integration_points.size(); PointNumber++) {
            // Compute element kinematics B, F, DN_DX ...
            CalculateKinematicVariables(rVariables, PointNumber, this->GetIntegrationMethod());
            double IntegrationWeight = rVariables.detJ0 * integration_points[PointNumber].Weight();
            TotalVolume += IntegrationWeight;

            for(IndexType i = 0; i < number_of_nodes; i++) {
                // Bh = Bh + sum(Bs(1:3, : )) * wg(iPG) * detJ;
                rVariables.Bh(i * 3 + 0) += rVariables.B(0, i * 3 + 0) * IntegrationWeight;
                rVariables.Bh(i * 3 + 1) += rVariables.B(1, i * 3 + 1) * IntegrationWeight;
                rVariables.Bh(i * 3 + 2) += rVariables.B(2, i * 3 + 2) * IntegrationWeight;
            }
        }
        for(IndexType i = 0; i < number_of_nodes * 3; i++) {
            rVariables.Bh(i) /= TotalVolume;
        }
    }
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

Matrix SmallDisplacementBbar::ComputeEquivalentF(const Vector& rStrainTensor)
{
    const SizeType dim = GetGeometry().WorkingSpaceDimension();
    Matrix F(dim,dim);

    if(dim == 2) {
        F(0,0) = 1.0 + rStrainTensor(0);
        F(0,1) = 0.5 * rStrainTensor(2);
        F(1,0) = 0.5 * rStrainTensor(2);
        F(1,1) = 1.0 + rStrainTensor(1);
    } else {
        F(0,0) = 1.0 + rStrainTensor(0);
        F(0,1) = 0.5 * rStrainTensor(3);
        F(0,2) = 0.5 * rStrainTensor(5);
        F(1,0) = 0.5 * rStrainTensor(3);
        F(1,1) = 1.0 + rStrainTensor(1);
        F(1,2) = 0.5 * rStrainTensor(4);
        F(2,0) = 0.5 * rStrainTensor(5);
        F(2,1) = 0.5 * rStrainTensor(4);
        F(2,2) = 1.0 + rStrainTensor(2);
    }

    return F;
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementBbar::CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rOutput,
        const ProcessInfo& rCurrentProcessInfo)
{
    const SizeType number_of_integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod()).size();
    if (rOutput.size() != number_of_integration_points)
        rOutput.resize(number_of_integration_points);

    if (rVariable == STRAIN_ENERGY) {
        const SizeType number_of_nodes = GetGeometry().size();
        const SizeType dimension = GetGeometry().WorkingSpaceDimension();
        const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

        KinematicVariablesBbar this_kinematic_variables(strain_size, dimension, number_of_nodes);
        ConstitutiveVariables this_constitutive_variables(strain_size);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions=Values.GetOptions();
        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, false);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        // Reading integration points
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

        // If strain has to be computed inside of the constitutive law with PK2
        Values.SetStrainVector(this_constitutive_variables.StrainVector); //this is the input  parameter

        // compute Hydrostatic B-Matrix
        this->SmallDisplacementBbar::CalculateHydrostaticDeformationMatrix(this_kinematic_variables);

        for(IndexType point_number = 0; point_number < integration_points.size(); point_number++) {
            // Compute element kinematics B, F, DN_DX ...
            CalculateKinematicVariablesBbar(this_kinematic_variables, point_number, integration_points);

            // Compute material response
            SetConstitutiveVariables(this_kinematic_variables, this_constitutive_variables,
                                           Values, point_number, integration_points);

            double StrainEnergy = 0.0;
            mConstitutiveLawVector[point_number]->CalculateValue(Values, STRAIN_ENERGY, StrainEnergy);
            rOutput[point_number] = StrainEnergy;
        }
    }
    else if (rVariable == VON_MISES_STRESS) {
        const SizeType number_of_nodes = GetGeometry().size();
        const SizeType dimension = GetGeometry().WorkingSpaceDimension();
        const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

        KinematicVariablesBbar this_kinematic_variables(strain_size, dimension, number_of_nodes);
        ConstitutiveVariables this_constitutive_variables(strain_size);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions=Values.GetOptions();
        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        Values.SetStrainVector(this_constitutive_variables.StrainVector);

        // Reading integration points
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( this->GetIntegrationMethod());

        // compute Hydrostatic B-Matrix
        this->SmallDisplacementBbar::CalculateHydrostaticDeformationMatrix(this_kinematic_variables);

        for(IndexType point_number = 0; point_number < integration_points.size(); point_number++) {
            // Compute element kinematics B, F, DN_DX ...
            CalculateKinematicVariablesBbar(this_kinematic_variables, point_number, integration_points);
            // Compute material response
            CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables,
                                           Values, point_number, integration_points,
                                           GetStressMeasure(), false);

            const Matrix stress_tensor = MathUtils<double>::StressVectorToTensor( this_constitutive_variables.StressVector );

            double sigma_equivalent = 0.0;

            if (dimension == 2) {
                sigma_equivalent = std::pow((stress_tensor(0,0) - stress_tensor(1,1)), 2.0) +
                                   3*(stress_tensor(0,1) * stress_tensor(1,0));
            }
            else {
                sigma_equivalent = 0.5*(std::pow((stress_tensor(0,0) - stress_tensor(1,1)), 2.0) +
                                        std::pow((stress_tensor(1,1) - stress_tensor(2,2)), 2.0) +
                                        std::pow((stress_tensor(2,2) - stress_tensor(0,0)), 2.0) +
                                        6*(stress_tensor(0,1) * stress_tensor(1,0) +
                                           stress_tensor(1,2) * stress_tensor(2,1) +
                                           stress_tensor(2,0) * stress_tensor(0,2)));
            }

            if( sigma_equivalent < 0.0 )
                rOutput[point_number] = 0.0;
            else
                rOutput[point_number] = std::sqrt(sigma_equivalent);
        }
    }
    else
    {
        BaseSolidElement::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
    }
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementBbar::CalculateOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
    )
{
    const SizeType number_of_integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod()).size();
    if ( rOutput.size() != number_of_integration_points )
        rOutput.resize( number_of_integration_points );

    if ( rVariable == CAUCHY_STRESS_VECTOR || rVariable == PK2_STRESS_VECTOR ) {
        // Create and initialize element variables:
        const bool is_rotated = IsElementRotated();
        const SizeType number_of_nodes = GetGeometry().size();
        const SizeType dimension = GetGeometry().WorkingSpaceDimension();
        const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

        KinematicVariablesBbar this_kinematic_variables(strain_size, dimension, number_of_nodes);
        ConstitutiveVariables this_constitutive_variables(strain_size);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions=Values.GetOptions();
        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        Values.SetStrainVector(this_constitutive_variables.StrainVector);

        // Reading integration points
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

        // compute Hydrostatic B-Matrix
        this->CalculateHydrostaticDeformationMatrix(this_kinematic_variables);

        // Reading integration points
        for ( IndexType point_number = 0; point_number < integration_points.size(); point_number++ )
        {
            // Compute element kinematics B, F, DN_DX ...
            CalculateKinematicVariablesBbar(this_kinematic_variables, point_number, integration_points);
            //call the constitutive law to update material variables
            if( rVariable == CAUCHY_STRESS_VECTOR) {
                // Compute material response
                CalculateConstitutiveVariables(this_kinematic_variables,
                                               this_constitutive_variables,
                                               Values,
                                               point_number,
                                               integration_points,
                                               ConstitutiveLaw::StressMeasure_Cauchy, is_rotated);
            }
            else {
                // Compute material response
                CalculateConstitutiveVariables(this_kinematic_variables,
                                               this_constitutive_variables,
                                               Values,
                                               point_number,
                                               integration_points,
                                               ConstitutiveLaw::StressMeasure_PK2, is_rotated);
            }

            if ( rOutput[point_number].size() != strain_size )
                rOutput[point_number].resize( strain_size, false );

            noalias(rOutput[point_number]) = this_constitutive_variables.StressVector;
        }
    }
    else if(rVariable == GREEN_LAGRANGE_STRAIN_VECTOR  || rVariable == ALMANSI_STRAIN_VECTOR) {
        // Create and initialize element variables:
        const SizeType number_of_nodes = GetGeometry().size();
        const SizeType dimension = GetGeometry().WorkingSpaceDimension();
        const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

        KinematicVariablesBbar this_kinematic_variables(strain_size, dimension, number_of_nodes);
        ConstitutiveVariables this_constitutive_variables(strain_size);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags &ConstitutiveLawOptions=Values.GetOptions();
        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, false);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        Values.SetStrainVector(this_constitutive_variables.StrainVector);

        const ConstitutiveLaw::StressMeasure this_stress_measure = rVariable == GREEN_LAGRANGE_STRAIN_VECTOR ? ConstitutiveLaw::StressMeasure_PK2 : ConstitutiveLaw::StressMeasure_Kirchhoff;

        // Reading integration points
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

        // compute Hydrostatic B-Matrix
        this->CalculateHydrostaticDeformationMatrix(this_kinematic_variables);

        //reading integration points
        for ( IndexType point_number = 0; point_number < integration_points.size(); point_number++ ) {
            // Compute element kinematics B, F, DN_DX ...
            CalculateKinematicVariablesBbar(this_kinematic_variables, point_number, integration_points);
            // Compute material response
            CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables,
                                           Values, point_number, integration_points, this_stress_measure, false);

            if ( rOutput[point_number].size() != strain_size)
                rOutput[point_number].resize( strain_size, false );

            noalias(rOutput[point_number]) = this_constitutive_variables.StrainVector;
        }
    }
    else {
        BaseSolidElement::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
    }
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementBbar::CalculateOnIntegrationPoints(
        const Variable<Matrix>& rVariable,
        std::vector<Matrix>& rOutput,
        const ProcessInfo& rCurrentProcessInfo)
{
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType number_of_integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod()).size();

    if (rOutput.size() != number_of_integration_points)
        rOutput.resize(number_of_integration_points);

    if ( rVariable == CAUCHY_STRESS_TENSOR || rVariable == PK2_STRESS_TENSOR ) {
        std::vector<Vector> stress_vector;

        if( rVariable == CAUCHY_STRESS_TENSOR )
            this->CalculateOnIntegrationPoints( CAUCHY_STRESS_VECTOR, stress_vector, rCurrentProcessInfo );
        else
            this->CalculateOnIntegrationPoints( PK2_STRESS_VECTOR, stress_vector, rCurrentProcessInfo );

        // Loop integration points
        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); point_number++ ) {
            if ( rOutput[point_number].size2() != dimension )
                rOutput[point_number].resize( dimension, dimension, false );

            noalias(rOutput[point_number]) = MathUtils<double>::StressVectorToTensor(stress_vector[point_number]);
        }
    }
    else if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR  || rVariable == ALMANSI_STRAIN_TENSOR) {
        std::vector<Vector> strain_vector;
        if( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR )
            CalculateOnIntegrationPoints( GREEN_LAGRANGE_STRAIN_VECTOR, strain_vector, rCurrentProcessInfo );
        else
            CalculateOnIntegrationPoints( ALMANSI_STRAIN_VECTOR, strain_vector, rCurrentProcessInfo );

        // Loop integration points
        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); point_number++ ) {
            if ( rOutput[point_number].size2() != dimension )
                rOutput[point_number].resize( dimension, dimension, false );

            noalias(rOutput[point_number]) = MathUtils<double>::StrainVectorToTensor(strain_vector[point_number]);
        }
    }
    else if (rVariable == CONSTITUTIVE_MATRIX) {
        // Create and initialize element variables:
        const SizeType number_of_nodes = GetGeometry().size();
        const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();
        const bool is_rotated = IsElementRotated();

        KinematicVariablesBbar this_kinematic_variables(strain_size, dimension, number_of_nodes);
        ConstitutiveVariables this_constitutive_variables(strain_size);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions=Values.GetOptions();
        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, false);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

        Values.SetStrainVector(this_constitutive_variables.StrainVector);
        Values.SetConstitutiveMatrix(this_constitutive_variables.D); //this is the output parameter

        // Read integration points
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

        // Compute Hydrostatic B-Matrix
        this->SmallDisplacementBbar::CalculateHydrostaticDeformationMatrix(this_kinematic_variables);

        // Read integration points
        for(IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); point_number++) {
            // Compute element kinematics B, F, DN_DX ...
            CalculateKinematicVariablesBbar(this_kinematic_variables, point_number, integration_points);
            // Compute material response
            CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables,
                                           Values, point_number, integration_points,
                                           GetStressMeasure(), is_rotated);

            if( rOutput[point_number].size2() != this_constitutive_variables.D.size2() )
                rOutput[point_number].resize( this_constitutive_variables.D.size1() , this_constitutive_variables.D.size2() , false );

            noalias(rOutput[point_number]) = this_constitutive_variables.D;
        }
    }
    else if ( rVariable == DEFORMATION_GRADIENT ) { // VARIABLE SET FOR TRANSFER PURPOUSES
        // Create and initialize element variables:
        const SizeType number_of_nodes = GetGeometry().size();
        const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

        KinematicVariablesBbar this_kinematic_variables(strain_size, dimension, number_of_nodes);
        ConstitutiveVariables this_constitutive_variables(strain_size);

        // Read integration points
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);

        // Compute Hydrostatic B-Matrix
        this->SmallDisplacementBbar::CalculateHydrostaticDeformationMatrix(this_kinematic_variables);

        // Read integration points
        for(IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); point_number++ ) {
            // Compute element kinematics B, F, DN_DX ...
            CalculateKinematicVariablesBbar(this_kinematic_variables, point_number, integration_points);
            // Compute material response
            CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables,
                                           Values, point_number, integration_points,
                                           GetStressMeasure(), false);

            if( rOutput[point_number].size2() != this_kinematic_variables.F.size2() )
                rOutput[point_number].resize( this_kinematic_variables.F.size1() , this_kinematic_variables.F.size2() , false );

            noalias(rOutput[point_number]) = this_kinematic_variables.F;
        }
    }
    else {
        BaseSolidElement::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
    }
}



//************************************************************************************
//************************************************************************************

void SmallDisplacementBbar::CalculateAndAddResidualVector(
    VectorType& rRightHandSideVector,
    const KinematicVariables& rThisKinematicVariables,
    const ProcessInfo& rCurrentProcessInfo,
    const array_1d<double, 3>& rBodyForce,
    const Vector& rStressVector,
    const double IntegrationWeight
    ) const
{
    KRATOS_TRY

    // Operation performed: rRightHandSideVector += ExtForce * IntegrationWeight
    this->CalculateAndAddExtForceContribution(
            rThisKinematicVariables.N,
            rCurrentProcessInfo,
            rBodyForce,
            rRightHandSideVector,
            IntegrationWeight);

    // Operation performed: rRightHandSideVector -= IntForce * IntegrationWeight
    noalias(rRightHandSideVector) -= IntegrationWeight * prod(trans(rThisKinematicVariables.B), rStressVector);

   KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementBbar::FinalizeSolutionStep( const ProcessInfo& rCurrentProcessInfo )
{
    // Create and initialize element variables:
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();
    const bool is_rotated = IsElementRotated();

    KinematicVariablesBbar this_kinematic_variables(strain_size, dimension, number_of_nodes);
    ConstitutiveVariables this_constitutive_variables(strain_size);

    // Create constitutive law parameters:
    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    // Set constitutive law flags:
    Flags& ConstitutiveLawOptions=Values.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

    Values.SetStrainVector(this_constitutive_variables.StrainVector);
    Values.SetStressVector(this_constitutive_variables.StressVector);
    Values.SetConstitutiveMatrix(this_constitutive_variables.D);

    // Read integration points
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

    // Compute Hydrostatic B-Matrix
    this->SmallDisplacementBbar::CalculateHydrostaticDeformationMatrix(this_kinematic_variables);

    // Read integration points
    for(IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); point_number++) {
        // Compute element kinematics B, F, DN_DX ...
        CalculateKinematicVariablesBbar(this_kinematic_variables, point_number, integration_points);

        // Compute constitutive law variables
        SetConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points);

        if (is_rotated)
            RotateToLocalAxes(Values, this_kinematic_variables);

        // Call the constitutive law to update material variables
        mConstitutiveLawVector[point_number]->FinalizeMaterialResponse(Values, GetStressMeasure());

        // TODO: To be deprecated
        mConstitutiveLawVector[point_number]->FinalizeSolutionStep(
                GetProperties(),
                GetGeometry(),
                row(GetGeometry().ShapeFunctionsValues(), point_number),
                rCurrentProcessInfo);
    }
}

//************************************************************************************
//************************************************************************************

    void SmallDisplacementBbar::save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseSolidElement);
    }

//************************************************************************************
//************************************************************************************

    void SmallDisplacementBbar::load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseSolidElement);
    }

} // Namespace Kratos
