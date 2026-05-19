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
#include "utilities/math_utils.h"

// Application includes
#include "small_displacement.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"
#include "structural_mechanics_application_variables.h"
#include "custom_response_functions/response_utilities/finite_difference_utility.h"


namespace Kratos
{
SmallDisplacement::SmallDisplacement( IndexType NewId, GeometryType::Pointer pGeometry )
    : BaseSolidElement( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

SmallDisplacement::SmallDisplacement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
        : BaseSolidElement( NewId, pGeometry, pProperties )
{
    //DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer SmallDisplacement::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<SmallDisplacement>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer SmallDisplacement::Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<SmallDisplacement>( NewId, pGeom, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

SmallDisplacement::~SmallDisplacement()
{
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer SmallDisplacement::Clone (
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    SmallDisplacement::Pointer p_new_elem = Kratos::make_intrusive<SmallDisplacement>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
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


//void SmallDisplacement::ComputeStiffnessDerivative(
//        std::span<const IAdjoint::DynamicVariable> Variables,
//        const ProcessInfo& rProcessInfo,
//        Matrix& rOutput) const {
//            KRATOS_TRY
//
//            KRATOS_CATCH("")
//}


/***********************************************************************************/
/***********************************************************************************/

bool SmallDisplacement::UseElementProvidedStrain() const
{
    return true;
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacement::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag
    )
{
    KRATOS_TRY;

    auto& r_geometry = this->GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();
    const SizeType strain_size = GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize();
    const bool is_rotated = IsElementRotated();

    KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
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
    const GeometryType::IntegrationPointsArrayType& integration_points = this->IntegrationPoints(this->GetIntegrationMethod());

    ConstitutiveLaw::Parameters Values(r_geometry,GetProperties(),rCurrentProcessInfo);

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
        CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());

        // Compute material response
        CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points, GetStressMeasure(), is_rotated);

        // Calculating weights for integration on the reference configuration
        int_to_reference_weight = GetIntegrationWeight(integration_points, point_number, this_kinematic_variables.detJ0);

        if ( dimension == 2 && GetProperties().Has( THICKNESS ))
            int_to_reference_weight *= GetProperties()[THICKNESS];

        if ( CalculateStiffnessMatrixFlag ) { // Calculation of the matrix is required
            // Contributions to stiffness matrix calculated on the reference config
            this->CalculateAndAddKm( rLeftHandSideMatrix, this_kinematic_variables.B, this_constitutive_variables.D, int_to_reference_weight );
        }

        if ( CalculateResidualVectorFlag ) { // Calculation of the matrix is required
            this->CalculateAndAddResidualVector(rRightHandSideVector, this_kinematic_variables, rCurrentProcessInfo, body_force, this_constitutive_variables.StressVector, int_to_reference_weight);
        }
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacement::CalculateKinematicVariables(
    KinematicVariables& rThisKinematicVariables,
    const IndexType PointNumber,
    const GeometryType::IntegrationMethod& rIntegrationMethod
    )
{
    const auto& r_geometry = GetGeometry();

    const GeometryType::IntegrationPointsArrayType& r_integration_points = this->IntegrationPoints(rIntegrationMethod);
    // Shape functions
    rThisKinematicVariables.N = r_geometry.ShapeFunctionsValues(rThisKinematicVariables.N, r_integration_points[PointNumber].Coordinates());

    rThisKinematicVariables.detJ0 = CalculateDerivativesOnReferenceConfiguration(rThisKinematicVariables.J0, rThisKinematicVariables.InvJ0, rThisKinematicVariables.DN_DX, PointNumber, rIntegrationMethod);

    KRATOS_ERROR_IF(rThisKinematicVariables.detJ0 < 0.0) << "WARNING:: ELEMENT ID: " << this->Id() << " INVERTED. DETJ0: " << rThisKinematicVariables.detJ0 << std::endl;

    // Compute B
    CalculateB( rThisKinematicVariables.B, rThisKinematicVariables.DN_DX, r_integration_points, PointNumber );

    // Compute equivalent F
    GetValuesVector(rThisKinematicVariables.Displacements);
    Vector strain_vector(mConstitutiveLawVector[0]->GetStrainSize());
    noalias(strain_vector) = prod(rThisKinematicVariables.B, rThisKinematicVariables.Displacements);
    ComputeEquivalentF(rThisKinematicVariables.F, strain_vector);
    rThisKinematicVariables.detF = MathUtils<double>::Det(rThisKinematicVariables.F);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacement::GetStiffnessInfluencingVariables(
    std::vector<IAdjoint::DynamicVariable>& rOutput,
    const ProcessInfo& rProcessInfo) const {
        // Collect variables from the material.
        this->GetProperties().GetValue(CONSTITUTIVE_LAW)->GetInfluencingVariables<IAdjoint::ResidualTerm::Stiffness>(
            rOutput,
            rProcessInfo);

        // Conditional thickness.
        const std::size_t dimension_count = this->GetGeometry().WorkingSpaceDimension();
        if (dimension_count == 2 && this->GetProperties().Has(THICKNESS)) {
            const bool has_thickness = std::find_if(
                rOutput.begin(),
                rOutput.end(),
                [] (const IAdjoint::DynamicVariable& r_variable) {return r_variable.Key() == THICKNESS.Key();}
            ) != rOutput.end();
            if (!has_thickness)
                rOutput.push_back(THICKNESS);
        }

        // Collect displacement and shape variables.
        const std::size_t node_count = this->GetGeometry().size();
        const std::array<std::array<const VariableData*,3>,2> variable_sets {
            std::array<const VariableData*,3> {
                &DISPLACEMENT_X,
                &DISPLACEMENT_Y,
                &DISPLACEMENT_Z},
            std::array<const VariableData*,3> {
                &SHAPE_X,
                &SHAPE_Y,
                &SHAPE_Z}};
        for (const auto& r_variables : variable_sets)
            for (std::size_t i_node=0ul; i_node<node_count; ++i_node)
                for (std::size_t i_dimension=0ul; i_dimension<dimension_count; ++i_dimension)
                    rOutput.push_back(IAdjoint::DynamicVariable(
                        *r_variables[i_dimension],
                        i_node));
}

void SmallDisplacement::GetLoadInfluencingVariables(
    std::vector<IAdjoint::DynamicVariable>& rOutput,
    const ProcessInfo& rProcessInfo) const {
        // Collect variables from the stiffness term.
        this->GetStiffnessInfluencingVariables(
            rOutput,
            rProcessInfo);

        // Collect variables from the material.
        {
            std::vector<IAdjoint::DynamicVariable> buffer;
            this->GetProperties().GetValue(CONSTITUTIVE_LAW)->GetInfluencingVariables<IAdjoint::ResidualTerm::Load>(
                buffer,
                rProcessInfo);
            buffer.erase(
                std::remove_if(
                    buffer.begin(),
                    buffer.end(),
                    [&rOutput] (const IAdjoint::DynamicVariable& r_variable) -> bool {
                        return std::find(
                            rOutput.begin(),
                            rOutput.end(),
                            r_variable
                        ) != rOutput.end();}),
                buffer.end());
            rOutput.insert(
                rOutput.end(),
                buffer.begin(),
                buffer.end());
        }

        // Collect inertial variables (from GetBodyForce).
        const Properties& r_properties = this->GetProperties();
        if (r_properties.Has(DENSITY)) rOutput.emplace_back(DENSITY);
        if (r_properties.Has(VOLUME_ACCELERATION)) rOutput.emplace_back(VOLUME_ACCELERATION);
        if (this->GetGeometry()[0].SolutionStepsDataHas(VOLUME_ACCELERATION)) {
            const std::size_t dimension_count = this->GetGeometry().WorkingSpaceDimension();
            const std::array<const VariableData*,3> components {
                &VOLUME_ACCELERATION_X,
                &VOLUME_ACCELERATION_Y,
                &VOLUME_ACCELERATION_Z};
            for (std::size_t i_node=0ul; i_node<this->GetGeometry().size(); ++i_node)
                for (std::size_t i_dimension=0ul; i_dimension<dimension_count; ++i_dimension)
                    rOutput.emplace_back(*components[i_dimension], i_node);
        }
}


void SmallDisplacement::ComputeStiffnessDerivative(
    Matrix& rOutput,
    std::span<const IAdjoint::DynamicVariable> Variables,
    const Vector& rValues,
    const ProcessInfo& rProcessInfo,
    int iBuffer) const {
        KRATOS_TRY
            // Define the finite differencing utility
            // that will be used for approximating derivatives.
            using Utility = AdjointFiniteDifferenceUtility<
                IAdjoint::ResidualTerm::Stiffness,
                Element>;

            // Set a default for the perturbation's magnitude.
            // This will be used for all variables that don't
            // have special treatment implemented.
            // Note that finding a suitable perturbation size
            // is key to computing an accurate approximation
            // robustly, so do try to implement special treatment
            // for every variable you plan to compute derivatives
            // for.
            double default_perturbation_magnitude = 1e-6;
            double shape_perturbation_magnitude = default_perturbation_magnitude;

            if (this->Has(PERTURBATION_SIZE)) {
                default_perturbation_magnitude = this->GetValue(PERTURBATION_SIZE);
                shape_perturbation_magnitude = default_perturbation_magnitude;
            } else if (this->GetProperties().Has(PERTURBATION_SIZE)) {
                default_perturbation_magnitude = this->GetProperties()[PERTURBATION_SIZE];
                shape_perturbation_magnitude = default_perturbation_magnitude;
            } else if (rProcessInfo.Has(PERTURBATION_SIZE)) {
                default_perturbation_magnitude = rProcessInfo[PERTURBATION_SIZE];
                shape_perturbation_magnitude = default_perturbation_magnitude;
            } else {
                // Assign a default perturbation magnitude for nodal coordinates if necesssary.
                const bool require_shape_derivatives = std::find_if(
                    Variables.begin(),
                    Variables.end(),
                    [] (const IAdjoint::DynamicVariable& rVariable) -> bool {
                        return rVariable.SourceKey() == SHAPE.SourceKey()
                            || rVariable.SourceKey() == DISPLACEMENT.SourceKey();}
                ) != Variables.end();
                if (require_shape_derivatives) {
                    shape_perturbation_magnitude = std::sqrt(default_perturbation_magnitude) * this->GetGeometry().DomainSize();
                } // if require_shape_derivatives
            }

            // Assemble a list of perturbations. These include the
            // variable to be perturbed, where it is stored (context)
            // and how big the perturbation should be.
            std::vector<Utility::Perturbation> perturbations;
            perturbations.reserve(Variables.size());
            for (const IAdjoint::DynamicVariable& r_variable : Variables) {
                bool found_variable = true;

                switch (r_variable.Key()) {
                    // Buffered variables in nodes.
                    case DISPLACEMENT_X.Key():
                    case DISPLACEMENT_Y.Key():
                    case DISPLACEMENT_Z.Key():
                        perturbations.push_back({
                            .mVariable = r_variable,
                            .mContext = Globals::DataLocation::NodeHistorical,
                            .mMagnitude = shape_perturbation_magnitude});
                        break;

                    // Unbuffered variables in nodes.
                    case SHAPE_X.Key():
                    case SHAPE_Y.Key():
                    case SHAPE_Z.Key():
                        perturbations.push_back({
                            .mVariable = r_variable,
                            .mContext = Globals::DataLocation::NodeNonHistorical,
                            .mMagnitude = shape_perturbation_magnitude});
                        break;

                    // Variables in Properties (that map to elements):
                    case THICKNESS.Key(): {
                        const double perturbation_magnitude = default_perturbation_magnitude * this->GetProperties().Data().template GetValue<double>(r_variable);
                        perturbations.push_back({
                            .mVariable = r_variable,
                            .mContext = Globals::DataLocation::Element,
                            .mMagnitude = perturbation_magnitude});
                        break;
                    }

                    // Unknown variable.
                    default:
                        found_variable = false;
                } // switch r_variable.Key()

                // Try subcomponents if the variable is unknown by the element.
                if (!found_variable) {
                    if (this->GetProperties().Has(CONSTITUTIVE_LAW)) {
                        std::vector<IAdjoint::DynamicVariable> constitutive_law_variables;
                        this->GetProperties()[CONSTITUTIVE_LAW]->GetInfluencingVariables<IAdjoint::ResidualTerm::Stiffness>(
                            constitutive_law_variables,
                            rProcessInfo);
                        found_variable = std::find(
                            constitutive_law_variables.begin(),
                            constitutive_law_variables.end(),
                            r_variable) != constitutive_law_variables.end();
                        if (found_variable) {
                            const double reference_value = this->GetProperties().Data().template GetValue<double>(r_variable);
                            const double perturbation_magnitude = reference_value
                                ? reference_value * default_perturbation_magnitude
                                : default_perturbation_magnitude;
                            perturbations.push_back({
                                .mVariable = r_variable,
                                .mContext = Globals::DataLocation::Element,
                                .mMagnitude = perturbation_magnitude});
                        }
                    }
                } // if not found_variable

                KRATOS_ERROR_IF_NOT(found_variable)
                    << "the stiffness term "
                    << "of element " << this->Id() << " "
                    << "does not depend on " << r_variable.Name();
            } // for r_variable in Variables

            // Instantiate the utility and compute finite differences.
            Utility utility;
            utility.FiniteDifferenceDerivative(
                *this,
                rValues,
                perturbations,
                rOutput,
                iBuffer,
                rProcessInfo);
        KRATOS_CATCH("")
}


void SmallDisplacement::ComputeLoadDerivative(
    Matrix& rOutput,
    std::span<const IAdjoint::DynamicVariable> Variables,
    const ProcessInfo& rProcessInfo,
    int iBuffer) const {
        KRATOS_TRY
            // Define the finite differencing utility
            // that will be used for approximating derivatives.
            using Utility = AdjointFiniteDifferenceUtility<
                IAdjoint::ResidualTerm::Load,
                Element>;

            // Set a default for the perturbation's magnitude.
            // This will be used for all variables that don't
            // have special treatment implemented.
            // Note that finding a suitable perturbation size
            // is key to computing an accurate approximation
            // robustly, so do try to implement special treatment
            // for every variable you plan to compute derivatives
            // for.
            double default_perturbation_magnitude = 1e-6;
            double shape_perturbation_magnitude = default_perturbation_magnitude;

            if (this->Has(PERTURBATION_SIZE)) {
                default_perturbation_magnitude = this->GetValue(PERTURBATION_SIZE);
                shape_perturbation_magnitude = default_perturbation_magnitude;
            } else if (this->GetProperties().Has(PERTURBATION_SIZE)) {
                default_perturbation_magnitude = this->GetProperties()[PERTURBATION_SIZE];
                shape_perturbation_magnitude = default_perturbation_magnitude;
            } else if (rProcessInfo.Has(PERTURBATION_SIZE)) {
                default_perturbation_magnitude = rProcessInfo[PERTURBATION_SIZE];
                shape_perturbation_magnitude = default_perturbation_magnitude;
            } else {
                // Assign a default perturbation magnitude for nodal coordinates if necesssary.
                const bool require_shape_derivatives = std::find_if(
                    Variables.begin(),
                    Variables.end(),
                    [] (const IAdjoint::DynamicVariable& rVariable) -> bool {
                        return rVariable.SourceKey() == SHAPE.SourceKey()
                            || rVariable.SourceKey() == DISPLACEMENT.SourceKey();}
                ) != Variables.end();
                if (require_shape_derivatives) {
                    shape_perturbation_magnitude = std::sqrt(default_perturbation_magnitude) * this->GetGeometry().DomainSize();
                } // if require_shape_derivatives
            }

            // Assemble a list of perturbations. These include the
            // variable to be perturbed, where it is stored (context)
            // and how big the perturbation should be.
            std::vector<Utility::Perturbation> perturbations;
            perturbations.reserve(Variables.size());
            for (const IAdjoint::DynamicVariable& r_variable : Variables) {
                bool found_variable = true;

                switch (r_variable.Key()) {
                    // Buffered variables in nodes.
                    case DISPLACEMENT_X.Key():
                    case DISPLACEMENT_Y.Key():
                    case DISPLACEMENT_Z.Key():
                        perturbations.push_back({
                            .mVariable = r_variable,
                            .mContext = Globals::DataLocation::NodeHistorical,
                            .mMagnitude = shape_perturbation_magnitude});
                        break;

                    // Unbuffered variables in nodes.
                    case SHAPE_X.Key():
                    case SHAPE_Y.Key():
                    case SHAPE_Z.Key():
                        perturbations.push_back({
                            .mVariable = r_variable,
                            .mContext = Globals::DataLocation::NodeNonHistorical,
                            .mMagnitude = shape_perturbation_magnitude});
                        break;

                    // Variables in Properties (that map to elements):
                    case DENSITY.Key():
                    case THICKNESS.Key(): {
                        const double perturbation_magnitude = default_perturbation_magnitude * this->GetProperties().Data().template GetValue<double>(r_variable);
                        perturbations.push_back({
                            .mVariable = r_variable,
                            .mContext = Globals::DataLocation::Element,
                            .mMagnitude = perturbation_magnitude});
                        break;
                    }

                    // Ambiguous variables.
                    // These may be defined on the element
                    // itself or on its nodes.
                    case VOLUME_ACCELERATION.Key():
                        if (0 <= r_variable.GetDynamicIndex()) {
                            perturbations.push_back({
                                .mVariable = r_variable,
                                .mContext = Globals::DataLocation::NodeHistorical,
                                .mMagnitude = default_perturbation_magnitude});
                        } else {
                            perturbations.push_back({
                                .mVariable = r_variable,
                                .mContext = Globals::DataLocation::Element,
                                .mMagnitude = default_perturbation_magnitude});
                        }
                        break;

                    // Unknown variable.
                    default:
                        found_variable = false;
                } // switch r_variable.Key()

                // Try subcomponents if the variable is unknown by the element.
                if (!found_variable) {
                    if (this->GetProperties().Has(CONSTITUTIVE_LAW)) {
                        std::vector<IAdjoint::DynamicVariable> constitutive_law_variables;
                        this->GetProperties()[CONSTITUTIVE_LAW]->GetInfluencingVariables<IAdjoint::ResidualTerm::Stiffness>(
                            constitutive_law_variables,
                            rProcessInfo);
                        found_variable = std::find(
                            constitutive_law_variables.begin(),
                            constitutive_law_variables.end(),
                            r_variable) != constitutive_law_variables.end();
                        if (found_variable) {
                            const double reference_value = this->GetProperties().Data().template GetValue<double>(r_variable);
                            const double perturbation_magnitude = reference_value
                                ? reference_value * default_perturbation_magnitude
                                : default_perturbation_magnitude;
                            perturbations.push_back({
                                .mVariable = r_variable,
                                .mContext = Globals::DataLocation::Element,
                                .mMagnitude = perturbation_magnitude});
                        }
                    }
                } // if not found_variable

                KRATOS_ERROR_IF_NOT(found_variable)
                    << "the load term "
                    << "of element " << this->Id() << " "
                    << "does not depend on " << r_variable.Name();
            } // for r_variable in Variables

            // Instantiate the utility and compute finite differences.
            Utility utility;
            utility.FiniteDifferenceDerivative(
                *this,
                perturbations,
                rOutput,
                iBuffer,
                rProcessInfo);
        KRATOS_CATCH("")
}


/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacement::SetConstitutiveVariables(
    KinematicVariables& rThisKinematicVariables,
    ConstitutiveVariables& rThisConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType PointNumber,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints
    )
{
    const auto& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();
    const SizeType mat_size = number_of_nodes * dimension;

    // Displacements vector
    Vector displacements(mat_size);
    GetValuesVector(displacements);

    // Compute strain
    noalias(rThisConstitutiveVariables.StrainVector) = prod(rThisKinematicVariables.B, displacements);

    // Here we essentially set the input parameters
    rValues.SetShapeFunctionsValues(rThisKinematicVariables.N); // shape functions
    rValues.SetDeterminantF(rThisKinematicVariables.detF); //assuming the determinant is computed somewhere else
    rValues.SetDeformationGradientF(rThisKinematicVariables.F); //F computed somewhere else

    // Here we set the space on which the results shall be written
    rValues.SetConstitutiveMatrix(rThisConstitutiveVariables.D); //assuming the determinant is computed somewhere else
    rValues.SetStressVector(rThisConstitutiveVariables.StressVector); //F computed somewhere else
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacement::CalculateB(
    Matrix& rB,
    const Matrix& rDN_DX,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
    const IndexType PointNumber
    ) const
{
    KRATOS_TRY;

    StructuralMechanicsElementUtilities::CalculateB(*this, rDN_DX, rB);

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacement::ComputeEquivalentF(
    Matrix& rF,
    const Vector& rStrainTensor
    ) const
{
    StructuralMechanicsElementUtilities::ComputeEquivalentF(*this, rStrainTensor, rF);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseSolidElement );
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseSolidElement );
}

} // Namespace Kratos


