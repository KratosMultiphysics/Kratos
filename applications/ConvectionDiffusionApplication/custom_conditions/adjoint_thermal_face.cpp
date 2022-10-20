// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Jordi Cotela
//

#include "adjoint_thermal_face.h"

#include "convection_diffusion_application_variables.h"

#include "includes/checks.h"
#include "includes/convection_diffusion_settings.h"
#include "utilities/line_sensitivity_utility.h"
#include "utilities/math_utils.h"

namespace Kratos
{

AdjointThermalFace::AdjointThermalFace(IndexType NewId, typename GeometryType::Pointer pGeometry):
    ThermalFace(NewId, pGeometry)
{}

AdjointThermalFace::AdjointThermalFace(
    IndexType NewId, typename GeometryType::Pointer pGeometry, Properties::Pointer pProperties):
    ThermalFace(NewId, pGeometry, pProperties)
{}

AdjointThermalFace::~AdjointThermalFace() {}

Condition::Pointer AdjointThermalFace::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<AdjointThermalFace>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}

Condition::Pointer AdjointThermalFace::Create(
    IndexType NewId,
    typename GeometryType::Pointer pGeometry,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<AdjointThermalFace>(NewId, pGeometry, pProperties);
}

void AdjointThermalFace::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Delegating LHS matrix to base class
    ThermalFace::CalculateLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);

    // Setting the RHS vector to zero
    noalias(rRightHandSideVector) = ZeroVector(rLeftHandSideMatrix.size2());
}

void AdjointThermalFace::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    const Geometry<Node<3>>& r_geom = this->GetGeometry();
    const unsigned int num_nodes = r_geom.PointsNumber();

    if (rRightHandSideVector.size() != num_nodes)
    {
        rRightHandSideVector.resize(num_nodes,false);
    }

    noalias(rRightHandSideVector) = ZeroVector(num_nodes);
}

void AdjointThermalFace::GetValuesVector(Vector& rValues, int Step) const
{
    const GeometryType& r_geom = this->GetGeometry();
    const unsigned int num_nodes = r_geom.PointsNumber();

    if (rValues.size() != num_nodes)
    {
        rValues.resize(num_nodes,false);
    }

    for (unsigned int i = 0; i < num_nodes; i++)
    {
        rValues[i] = r_geom[i].FastGetSolutionStepValue(ADJOINT_HEAT_TRANSFER, Step);
    }
}

void AdjointThermalFace::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const GeometryType& r_geom = this->GetGeometry();
    const unsigned int num_nodes = r_geom.PointsNumber();

    if (rResult.size() != num_nodes)
    {
        rResult.resize(num_nodes,false);
    }

    for (unsigned int i = 0; i < num_nodes; i++)
    {
        rResult[i] = r_geom[i].GetDof(ADJOINT_HEAT_TRANSFER).EquationId();
    }
}

void AdjointThermalFace::GetDofList(
    DofsVectorType& rConditionDofList, const ProcessInfo& rCurrentProcessInfo) const
{
    const GeometryType& r_geom = this->GetGeometry();
    const unsigned int num_nodes = r_geom.PointsNumber();

    if (rConditionDofList.size() != num_nodes)
    {
        rConditionDofList.resize(num_nodes);
    }

    for (unsigned int i = 0; i < num_nodes; i++)
    {
        rConditionDofList[i] = r_geom[i].pGetDof(ADJOINT_HEAT_TRANSFER);
    }
}

int AdjointThermalFace::Check(const ProcessInfo& rProcessInfo) const
{
    KRATOS_TRY
    KRATOS_ERROR_IF_NOT(rProcessInfo.Has(CONVECTION_DIFFUSION_SETTINGS)) << "No CONVECTION_DIFFUSION_SETTINGS defined in ProcessInfo." << std::endl;
    ConvectionDiffusionSettings::Pointer p_settings = rProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
    auto& r_settings = *p_settings;

    KRATOS_ERROR_IF_NOT(r_settings.IsDefinedSurfaceSourceVariable()) << "No Surface Source Variable defined in provided CONVECTION_DIFFUSION_SETTINGS." << std::endl;

    const Variable<double>& r_surface_source_variable = r_settings.GetSurfaceSourceVariable();

    const GeometryType& r_geom = this->GetGeometry();
    const unsigned int num_nodes = r_geom.PointsNumber();
    for (unsigned int i = 0; i < num_nodes; i++)
    {
        const Node<3>& r_node = r_geom[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_HEAT_TRANSFER, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_surface_source_variable, r_node);
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_HEAT_TRANSFER, r_node);
    }

    KRATOS_CATCH("")
    return ThermalFace::Check(rProcessInfo);
}

std::string AdjointThermalFace::Info() const
{
    std::stringstream buffer;
    buffer << "AdjointThermalFace #" << this->Id();
    return buffer.str();
}

void AdjointThermalFace::PrintInfo(std::ostream& rOStream) const
{
    const GeometryType& r_geom = this->GetGeometry();
    const unsigned int dimension = r_geom.WorkingSpaceDimension();
    const unsigned int num_nodes = r_geom.PointsNumber();
    rOStream << "AdjointThermalFace" << dimension << "D" << num_nodes << "N";
}

void AdjointThermalFace::CalculateSensitivityMatrix(
    const Variable<array_1d<double, 3>>& rDesignVariable,
    Matrix& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    const GeometryType& r_geom = this->GetGeometry();
    const unsigned int dimension = r_geom.WorkingSpaceDimension();
    const unsigned int local_dimension = r_geom.LocalSpaceDimension();
    const unsigned int num_nodes = r_geom.PointsNumber();
    const unsigned int sensitivity_size = dimension * num_nodes;

    if (rOutput.size1() != sensitivity_size || rOutput.size2() != num_nodes)
    {
        rOutput.resize(sensitivity_size,num_nodes,false);
    }
    noalias(rOutput) = ZeroMatrix(sensitivity_size,num_nodes);

    const auto integration_method = this->GetIntegrationMethod();
    const auto integration_points = r_geom.IntegrationPoints(integration_method);
    const unsigned int num_integration_points = integration_points.size();


    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
    auto& r_settings = *p_settings;

    const Variable<double>& r_flux_variable = r_settings.GetSurfaceSourceVariable();
    const Variable<double>& r_unknown_variable = r_settings.GetUnknownVariable();

    Vector nodal_flux = ZeroVector(num_nodes);
    Vector nodal_unknown = ZeroVector(num_nodes);
    for (unsigned int i = 0; i < num_nodes; i++)
    {
        nodal_flux[i] = r_geom[i].FastGetSolutionStepValue(r_flux_variable);
        nodal_unknown[i] = r_geom[i].FastGetSolutionStepValue(r_unknown_variable);
    }

    const Properties& r_properties = this->GetProperties();
    const double ambient_temperature = r_properties.GetValue(AMBIENT_TEMPERATURE);
    const double convection_coefficient = r_properties.GetValue(CONVECTION_COEFFICIENT);
    const double emissivity = r_properties.GetValue(EMISSIVITY);

    if (rDesignVariable == SHAPE_SENSITIVITY)
    {
        Matrix shape_function_local_gradients(num_nodes,local_dimension);
        Matrix jacobian(dimension,local_dimension);

        Matrix shape_functions = r_geom.ShapeFunctionsValues(integration_method);

        for (unsigned int g = 0; g < num_integration_points; g++)
        {
            noalias(shape_function_local_gradients) = r_geom.ShapeFunctionLocalGradient(g, integration_method);
            noalias(jacobian) = this->GetJacobian(integration_method, g);
            LineSensitivityUtility sensitivity_utility(jacobian,shape_function_local_gradients);

            Vector N = row(shape_functions, g);
            double q_gauss = inner_prod(N,nodal_flux);
            double value_gauss = inner_prod(N, nodal_unknown);
            q_gauss -= convection_coefficient*(value_gauss - ambient_temperature); // add flux contribution from convection condition
            q_gauss -= emissivity * StefanBoltzmann * (std::pow(value_gauss,4) - std::pow(ambient_temperature,4)); // flux contribution from radiation

            const double weight = integration_points[g].Weight();

            for (auto s = ShapeParameter::Sequence(num_nodes, dimension); s; ++s)
            {
                const auto& deriv = s.CurrentValue();
                double jacobian_sensitivity;
                sensitivity_utility.CalculateSensitivity(deriv, jacobian_sensitivity);

                // d/dX_l (w * J * N_i * N_j * q_j) = w * N_i * N_j * q_j * dJ/dX_l
                // Note that N_j * q_j = q_gauss
                for (unsigned int i = 0; i < num_nodes; i++)
                {
                    rOutput(deriv.NodeIndex * dimension + deriv.Direction, i) -= weight * N[i] * q_gauss * jacobian_sensitivity;
                }
            }
        }
    }
    else
    {
        KRATOS_ERROR << "Unsupported design variable " << rDesignVariable << std::endl;
    }

    KRATOS_CATCH("")
}

void AdjointThermalFace::CalculateSensitivityMatrix(
    const Variable<double>& rDesignVariable,
    Matrix& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType& r_geom = this->GetGeometry();
    const unsigned int num_nodes = r_geom.PointsNumber();

    if (rDesignVariable == CONDUCTIVITY_SENSITIVITY) {

        const unsigned int sensitivity_size =  num_nodes;

        if (rOutput.size1() != sensitivity_size || rOutput.size2() != num_nodes) {
            rOutput.resize(sensitivity_size,num_nodes,false);
        }
        noalias(rOutput) = ZeroMatrix(sensitivity_size,num_nodes);
    }
    else if ( ThermalFace::GetProperties().Has(rDesignVariable) ) {
        // define working variables
        Vector RHS;
        Vector RHS_perturbed;

        ThermalFace::CalculateRightHandSide(RHS, rCurrentProcessInfo);
        if ( (rOutput.size1() != 1) || (rOutput.size2() != RHS.size() ) ) {
            rOutput.resize(1, RHS.size(), false);
        }

        // Save property pointer
        Properties::Pointer p_global_properties = ThermalFace::pGetProperties();

        // Create new property and assign it to the condition
        Properties::Pointer p_local_property(Kratos::make_shared<Properties>(Properties(*p_global_properties)));
        ThermalFace::SetProperties(p_local_property);

        // perturb the design variable
        const double current_property_value = ThermalFace::GetProperties()[rDesignVariable];
        double alpha = 1e-5; //MFusseder TODO
        double perturbation_size = alpha;
        if (current_property_value > 0 ){ //MFusseder TODO rework
            perturbation_size = current_property_value * alpha;
        }
        p_local_property->SetValue(rDesignVariable, (current_property_value + perturbation_size));

        // Compute RHS after perturbation
        ThermalFace::CalculateRightHandSide(RHS_perturbed, rCurrentProcessInfo);

        // Compute derivative of RHS w.r.t. design variable with finite differences
        double element_area = r_geom.Area();
        for(IndexType i = 0; i < RHS_perturbed.size(); ++i) {
            rOutput(0, i) = (RHS_perturbed[i] - RHS[i]) / perturbation_size; //TODO MFusseder check for sign
            rOutput(0, i) /= element_area;
        }

        // Give thermal face condition original properties back
        ThermalFace::SetProperties(p_global_properties);
    }
    else {
        rOutput = ZeroMatrix(1, num_nodes);
    }

    KRATOS_CATCH("")
}

void AdjointThermalFace::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    const GeometryType& r_geometry = this->GetGeometry();
    const unsigned int num_nodes = r_geometry.PointsNumber();
    const auto &r_gauss_pts = r_geometry.IntegrationPoints(this->GetIntegrationMethod());
    const unsigned int n_gauss = r_gauss_pts.size();
    const MatrixType N_container = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

    rValues.resize(n_gauss);

    if( rVariable == PSEUDO_THERMAL_FACE_CONVECTION_COEFFICIENT || rVariable == PSEUDO_THERMAL_FACE_EMISSIVITY || rVariable == PSEUDO_THERMAL_FACE_AMBIENT_TEMPERATURE) {

        Vector node_multipliers = ZeroVector(num_nodes);

        ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
        auto& r_settings = *p_settings;
        const Variable<double>& r_unknown_variable = r_settings.GetUnknownVariable();
        Vector nodal_unknown = ZeroVector(num_nodes);
        for (unsigned int i = 0; i < num_nodes; i++) {
            nodal_unknown[i] = r_geometry[i].FastGetSolutionStepValue(r_unknown_variable);
        }

        const Properties& r_properties = this->GetProperties();
        const double ambient_temperature = r_properties.GetValue(AMBIENT_TEMPERATURE);
        const double convection_coefficient = r_properties.GetValue(CONVECTION_COEFFICIENT);
        const double emissivity = r_properties.GetValue(EMISSIVITY);

        if (rVariable == PSEUDO_THERMAL_FACE_CONVECTION_COEFFICIENT && r_properties.Has(CONVECTION_COEFFICIENT)) {
            for (unsigned int i = 0; i < num_nodes; i++) {
                node_multipliers[i] = -1 * (nodal_unknown[i] - ambient_temperature); //TODO MFusseder check for sign
            }
        }
        else if (rVariable == PSEUDO_THERMAL_FACE_EMISSIVITY && r_properties.Has(EMISSIVITY)) {
            for (unsigned int i = 0; i < num_nodes; i++) {
                 node_multipliers[i] = -1 * (StefanBoltzmann * (std::pow(nodal_unknown[i], 4) - pow(ambient_temperature, 4))); //TODO MFusseder check for sign
            }
        }
        else {
            for (unsigned int i = 0; i < num_nodes; i++) {
                if (r_properties.Has(AMBIENT_TEMPERATURE)) {
                    node_multipliers[i] = -1 * (- convection_coefficient - emissivity * StefanBoltzmann * 4 * pow(ambient_temperature, 3)); //TODO MFusseder check for sign
                }
            }
        }

        // Gauss pts. loop
        for (unsigned int g = 0; g < n_gauss; g++) {
            Vector N = row(N_container, g);
            rValues[g] = inner_prod(N, node_multipliers); //TODO MFusseder check for sign
        }
    }
    else if(rVariable == CONVECTION_COEFFICIENT_GP_SENSITIVITY || rVariable == EMISSIVITY_GP_SENSITIVITY || rVariable == AMBIENT_TEMPERATURE_GP_SENSITIVITY) {
        Vector gauss_pts_J_det = ZeroVector(n_gauss);
        r_geometry.DeterminantOfJacobian(gauss_pts_J_det, this->GetIntegrationMethod());

        std::vector<double> pseudo_quantity;
        if( rVariable == CONVECTION_COEFFICIENT_GP_SENSITIVITY) {
            this->CalculateOnIntegrationPoints(PSEUDO_THERMAL_FACE_CONVECTION_COEFFICIENT, pseudo_quantity, rCurrentProcessInfo);
        }
        else if( rVariable == EMISSIVITY_GP_SENSITIVITY) {
            this->CalculateOnIntegrationPoints(PSEUDO_THERMAL_FACE_EMISSIVITY, pseudo_quantity, rCurrentProcessInfo);
        }
        else {
            this->CalculateOnIntegrationPoints(PSEUDO_THERMAL_FACE_AMBIENT_TEMPERATURE, pseudo_quantity, rCurrentProcessInfo);
        }

        Vector adjoint_variables = ZeroVector(num_nodes);
        this->GetValuesVector(adjoint_variables, 0);

        // Gauss pts. loop
        for (unsigned int g = 0; g < n_gauss; g++) {
            Vector N = row(N_container, g);
            double adj_variable_gp = inner_prod(N, adjoint_variables);
            rValues[g] = adj_variable_gp * pseudo_quantity[g];
        }
    }
    else if (this->Has(rVariable)) {
        rValues[0] = this->GetValue(rVariable);
        for (unsigned int g = 1; g < n_gauss; ++g) {
            rValues[g] = rValues[0];
        }
    }
    else {
        //KRATOS_ERROR << "Unsupported output variable." << std::endl;
    }

}


typename AdjointThermalFace::MatrixType AdjointThermalFace::GetJacobian(
    GeometryData::IntegrationMethod QuadratureOrder,
    unsigned int IntegrationPointIndex) const
{
    const auto& r_geometry = this->GetGeometry();
    const auto& rDN_De = r_geometry.ShapeFunctionLocalGradient(IntegrationPointIndex, QuadratureOrder);
    MatrixType jacobian(r_geometry.WorkingSpaceDimension(), r_geometry.LocalSpaceDimension());
    MatrixType coordinates(r_geometry.WorkingSpaceDimension(), r_geometry.PointsNumber());

    for (unsigned int i = 0; i < r_geometry.PointsNumber(); i++)
    {
        const auto& r_coordinates = r_geometry[i].Coordinates();
        for (unsigned int d = 0; d < r_geometry.WorkingSpaceDimension(); d++)
        {
            coordinates(d,i) = r_coordinates[d];
        }
    }

    noalias(jacobian) = prod(coordinates,rDN_De);
    return jacobian;
}

}