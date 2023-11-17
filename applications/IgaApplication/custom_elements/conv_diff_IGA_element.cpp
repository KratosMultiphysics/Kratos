// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolò Antonelli
//

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/convection_diffusion_settings.h"
#include "utilities/geometry_utilities.h"

// Application includes
#include "custom_elements/conv_diff_IGA_element.h"
#include "utilities/math_utils.h"

// NICO
// #include "convection_diffusion_application.h"


namespace Kratos
{

template<std::size_t TDim>
ConvDiffIGAElement<TDim>::ConvDiffIGAElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : Element(
        NewId,
        pGeometry)
{
}

template<std::size_t TDim>
ConvDiffIGAElement<TDim>::ConvDiffIGAElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : Element(
        NewId,
        pGeometry,
        pProperties)
{
}

template<std::size_t TDim>
Element::Pointer ConvDiffIGAElement<TDim>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<ConvDiffIGAElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

template<std::size_t TDim>
Element::Pointer ConvDiffIGAElement<TDim>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<ConvDiffIGAElement>(NewId, pGeom, pProperties);
}

// Deconstructor
template<std::size_t TDim>
ConvDiffIGAElement<TDim>::~ConvDiffIGAElement()
{
}




template<std::size_t TDim>
void ConvDiffIGAElement<TDim>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                            VectorType& rRightHandSideVector,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const ProcessInfo& r_process_info = rCurrentProcessInfo;
    ConvectionDiffusionSettings::Pointer p_settings = r_process_info[CONVECTION_DIFFUSION_SETTINGS];
    auto& r_settings = *p_settings;

    const Variable<double>& r_unknown_var = r_settings.GetUnknownVariable(); // Temperature
    const Variable<double>& r_diffusivity_var = r_settings.GetDiffusionVariable(); // Conductivity
    const Variable<double>& r_volume_source_var = r_settings.GetVolumeSourceVariable(); // HeatFlux

    const auto& r_geometry = GetGeometry();
    const unsigned int number_of_points = r_geometry.size();  // = 9
    const unsigned int dim = r_geometry.WorkingSpaceDimension();

    //Element variables
    KRATOS_WATCH(rCurrentProcessInfo)
    double theta = rCurrentProcessInfo[TIME_INTEGRATION_THETA];
    double dyn_st_beta = 0.0 ; // rCurrentProcessInfo[DYNAMIC_TAU]; // Why can't have access to DYNAMIC_TAU?
    const double delta_t = rCurrentProcessInfo[DELTA_TIME];
    double dt_inv = 1.0 / delta_t;
    double lumping_factor = 1.00 / double(number_of_points);

    double conductivity = 0.0;
    double specific_heat = 0.0;
    double density = 0.0;
    double beta = 0.0;
    double div_v = 0.0;
    ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
    const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();
    // exit(0);

    //resizing as needed the LHS
    if(rLeftHandSideMatrix.size1() != number_of_points)
        rLeftHandSideMatrix.resize(number_of_points,number_of_points,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_points,number_of_points); //resetting LHS
    
    // resizing as needed the RHS
    if(rRightHandSideVector.size() != number_of_points)
        rRightHandSideVector.resize(number_of_points,false);
    noalias(rRightHandSideVector) = ZeroVector(number_of_points); //resetting RHS

    // reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());
    const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    const Matrix& N_gausspoint = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

    
    // Initialize Jacobian
    // Matrix J0;
    // Initialize DN_DX
    Matrix DN_DX(number_of_points,dim);
    Matrix DN_DPSI(number_of_points,dim);
    Matrix InvJ0(dim,dim);
    Vector temp(number_of_points);

    Vector heat_flux_local(number_of_points);
    Vector nodal_conductivity(number_of_points);
    for(unsigned int node_element = 0; node_element<number_of_points; node_element++)
    {
        heat_flux_local[node_element] = r_geometry[node_element].FastGetSolutionStepValue(r_volume_source_var);
        nodal_conductivity[node_element] = r_geometry[node_element].FastGetSolutionStepValue(r_diffusivity_var);
        // heat_flux_local[node_element] = -2.0*r_geometry.Center().X()*(r_geometry.Center().X()-2.0) - 2.0*r_geometry.Center().Y()*(r_geometry.Center().Y()-2.0) ;
        // heat_flux_local[node_element] = -0.0 ;

    }
    
    std::vector<array_1d<double,3>> v(number_of_points);
    std::vector<array_1d<double,3>> vold(number_of_points);
    for(std::size_t i_point = 0; i_point < integration_points.size(); ++i_point) {
        Matrix J0;
        r_geometry.Jacobian(J0, i_point, this->GetIntegrationMethod());
        double DetJ0;
        // Calculating inverse jacobian and jacobian determinant
        MathUtils<double>::InvertMatrix(J0,InvJ0,DetJ0);
        noalias(DN_DX) = prod(DN_De[i_point],InvJ0);
        // Velocity
        v[i_point] = ZeroVector(3) ;
        vold[i_point]=  ZeroVector(3);
        const Variable<array_1d<double, 3 > >& rVelocityVar = my_settings->GetVelocityVariable();
        v[i_point] = GetGeometry()[i_point].FastGetSolutionStepValue(rVelocityVar);
        vold[i_point] = GetGeometry()[i_point].FastGetSolutionStepValue(rVelocityVar,1);
    }
    KRATOS_WATCH(v)
    KRATOS_WATCH(v[0])
    // CONVECTION
    // Computing the divergence
    for (unsigned int i = 0; i < number_of_points; i++) {
        for(unsigned int k=0; k<TDim; k++)
        {
            div_v += DN_DX(i,k)*(v[i][k] * theta + vold[i][k] * (1.0-theta));
        }
    }
    //Some auxilary definitions
    Matrix aux1 = ZeroMatrix(number_of_points, number_of_points);
    Matrix aux2 = ZeroMatrix(number_of_points, number_of_points);
    Matrix tmp(number_of_points, TDim);

    // Calculate h
    double h =0.0;
    for(unsigned int i=0; i<number_of_points; i++){
        double h_inv = 0.0;
        for(unsigned int k=0; k<TDim; k++)  { h_inv += DN_DX(i,k)*DN_DX(i,k); }
        h += 1.0/h_inv;
    }
    h = sqrt(h)/static_cast<double>(number_of_points);

    for(unsigned int igauss=0; igauss<number_of_points; igauss++){
        auto N = row(N_gausspoint,igauss); // these are the N which correspond to the gauss point "i_point"

        //obtain the velocity in the middle of the tiem step
        array_1d<double, TDim > vel_gauss=ZeroVector(TDim);
        for (unsigned int i = 0; i < number_of_points; i++)
        {
                for(unsigned int k=0; k<TDim; k++)
                vel_gauss[k] += N[i] * (v[i][k] * theta + vold[i][k] * (1.0-theta) );
        }
        const double norm_vel = norm_2(vel_gauss);
        Vector a_dot_grad = prod(DN_DX, vel_gauss);

        // CalculateTau 
        double inv_tau = dyn_st_beta * dt_inv;
        inv_tau += 2.0 * norm_vel / h + beta*  div_v;
        // Dynamic and convection terms are multiplyied by density*specific_heat to have consistent dimensions
        inv_tau *= density * specific_heat;
        inv_tau += 4.0 * conductivity / (h*h);
        inv_tau = std::max(inv_tau, 1e-2);
        const double tau = (density * specific_heat) / inv_tau;
 
        //terms multiplying dphi/dt (aux1)
        noalias(aux1) += (1.0 + tau * beta * div_v) * outer_prod(N, N);
        noalias(aux1) +=  tau * outer_prod(a_dot_grad, N);

        //terms which multiply the gradient of phi
        noalias(aux2) += (1.0 + tau * beta* div_v) * outer_prod(N, a_dot_grad);
        noalias(aux2) += tau * outer_prod(a_dot_grad, a_dot_grad);
    }
    KRATOS_WATCH(aux1)
    KRATOS_WATCH(aux2)

    // //adding the second and third term in the formulation
    // noalias(rLeftHandSideMatrix)  = (dt_inv*density*specific_heat + theta*beta*div_v)*aux1;
    // noalias(rRightHandSideVector) = (dt_inv*density*specific_heat - (1.0-theta)*beta*div_v)*prod(aux1,phi_old);

    // //adding the diffusion
    // noalias(rLeftHandSideMatrix)  += (conductivity * theta * prod(DN_DX, trans(DN_DX)))*static_cast<double>(number_of_points);
    // noalias(rRightHandSideVector) -= prod((conductivity * (1.0-theta) * prod(DN_DX, trans(DN_DX))),phi_old)*static_cast<double>(number_of_points) ;

    // //terms in aux2
    // noalias(rLeftHandSideMatrix) += density*specific_heat*theta*aux2;
    // noalias(rRightHandSideVector) -= density*specific_heat*(1.0-theta)*prod(aux2,phi_old);

    // // volume source terms (affecting the RHS only)
    // noalias(rRightHandSideVector) += prod(aux1, volumetric_source);

    // //take out the dirichlet part to finish computing the residual
    // noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, phi);
    exit(0);

    

    for(std::size_t i_point = 0; i_point < integration_points.size(); ++i_point)
    {
        const IndexType IntegrationPointIndex = i_point ;
        Matrix J0;
        
        r_geometry.Jacobian(J0, IntegrationPointIndex, this->GetIntegrationMethod());
        double DetJ0;

        // Calculating inverse jacobian and jacobian determinant
        MathUtils<double>::InvertMatrix(J0,InvJ0,DetJ0);

        // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
        noalias(DN_DX) = prod(DN_De[i_point],InvJ0);
        // Calculating the PARAMETER SPACE derivatives
        Matrix Identity_Matrix = ZeroMatrix(2,2);
        Identity_Matrix(0,0) = 1.0;
        Identity_Matrix(1,1) = 1.0;
        noalias(DN_DPSI) = prod(DN_De[i_point],Identity_Matrix);


        // NEW!!  ->  WATCH OUT WHEN YOU DEAL WITH 2D problems, without this it invents the third component of DN_DX
        for (size_t i = 0; i < DN_DX.size1(); ++i) {
            DN_DX(i, 2) = 0.0;
            DN_DPSI(i, 2) = 0.0;
        }

        auto N = row(N_gausspoint,i_point); // these are the N which correspond to the gauss point "i_point"
        const double IntToReferenceWeight = integration_points[i_point].Weight() * abs(DetJ0); // I have added a minus here because the determinant is negative
        // Watch out DetJ0 = -1

        // Also the conductivity is multiplied by the value of the shape function at the GaussPoint
        const double conductivity_gauss = inner_prod(N, nodal_conductivity);
        noalias(rLeftHandSideMatrix) += IntToReferenceWeight * conductivity_gauss * prod(DN_DX, trans(DN_DX)); //


        // Calculating the local RHS
        const double qgauss = inner_prod(N, heat_flux_local);
        noalias(rRightHandSideVector) += IntToReferenceWeight * qgauss * N;
    }


    // RHS = ExtForces - K*temp;
    for (unsigned int i = 0; i < number_of_points; i++) {
        temp[i] = r_geometry[i].GetSolutionStepValue(r_unknown_var);

        // std::ofstream outputFile("Id_active_control_points.txt", std::ios::app);
        // outputFile << r_geometry[i].GetId() <<"\n";
        // outputFile.close();
    }

    // RHS -= K*temp
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp);

    KRATOS_CATCH("")
}


// From classical Laplacian
template<std::size_t TDim>
void ConvDiffIGAElement<TDim>::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{   
    KRATOS_WATCH('LHS___111')
    exit(0) ;
    VectorType temp(0);
    CalculateLocalSystem(rLeftHandSideMatrix, temp, rCurrentProcessInfo);
}


// From classical Laplacian
template<std::size_t TDim>
void ConvDiffIGAElement<TDim>::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType temp(0,0);
    CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
}





template<std::size_t TDim>
void ConvDiffIGAElement<TDim>::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const
{
    const ProcessInfo& r_process_info = rCurrentProcessInfo;
    ConvectionDiffusionSettings::Pointer p_settings = r_process_info[CONVECTION_DIFFUSION_SETTINGS];
    auto& r_settings = *p_settings;

    const Variable<double>& r_unknown_var = r_settings.GetUnknownVariable();

    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    if(rResult.size() != number_of_nodes)
        rResult.resize(number_of_nodes);
    for (unsigned int i=0; i<number_of_nodes; i++)
    {
        rResult[i] = GetGeometry()[i].GetDof(r_unknown_var).EquationId();
    }
}

//************************************************************************************
//************************************************************************************
template<std::size_t TDim>
void ConvDiffIGAElement<TDim>::GetDofList(DofsVectorType& ElementalDofList,const ProcessInfo& rCurrentProcessInfo) const
{
    const ProcessInfo& r_process_info = rCurrentProcessInfo;
    ConvectionDiffusionSettings::Pointer p_settings = r_process_info[CONVECTION_DIFFUSION_SETTINGS];
    auto& r_settings = *p_settings;

    const Variable<double>& r_unknown_var = r_settings.GetUnknownVariable();

    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    if(ElementalDofList.size() != number_of_nodes)
        ElementalDofList.resize(number_of_nodes);
    for (unsigned int i=0; i<number_of_nodes; i++)
    {
        ElementalDofList[i] = GetGeometry()[i].pGetDof(r_unknown_var);
    }
}


template<std::size_t TDim>
int ConvDiffIGAElement<TDim>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(CONVECTION_DIFFUSION_SETTINGS)) << "No CONVECTION_DIFFUSION_SETTINGS defined in ProcessInfo." << std::endl;
    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
    auto& r_settings = *p_settings;

    KRATOS_ERROR_IF_NOT(r_settings.IsDefinedUnknownVariable()) << "No Unknown Variable defined in provided CONVECTION_DIFFUSION_SETTINGS." << std::endl;
    KRATOS_ERROR_IF_NOT(r_settings.IsDefinedDiffusionVariable()) << "No Diffusion Variable defined in provided CONVECTION_DIFFUSION_SETTINGS." << std::endl;
    KRATOS_ERROR_IF_NOT(r_settings.IsDefinedVolumeSourceVariable()) << "No Volume Source Variable defined in provided CONVECTION_DIFFUSION_SETTINGS." << std::endl;

    const Variable<double>& r_unknown_var = r_settings.GetUnknownVariable();
    const Variable<double>& r_diffusivity_var = r_settings.GetDiffusionVariable();
    const Variable<double>& r_volume_source_var = r_settings.GetVolumeSourceVariable();

    const auto& r_geom = GetGeometry();

    for (unsigned int i = 0; i < r_geom.PointsNumber(); i++)
    {
        const auto& r_node = r_geom[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_unknown_var, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_diffusivity_var, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_volume_source_var, r_node);

        KRATOS_CHECK_DOF_IN_NODE(r_unknown_var, r_node);
    }

    return Element::Check(rCurrentProcessInfo);
}


template<std::size_t TDim>
Element::IntegrationMethod ConvDiffIGAElement<TDim>::GetIntegrationMethod() const
{
    return GeometryData::IntegrationMethod::GI_GAUSS_1;
}




template<std::size_t TDim>
void ConvDiffIGAElement<TDim>::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
    {
    const auto& r_geometry = GetGeometry();
    const SizeType nb_nodes = r_geometry.size();

    // Integration Points
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints();
    // Shape function values
    const Matrix& r_N = r_geometry.ShapeFunctionsValues();

    const ProcessInfo& r_process_info = rCurrentProcessInfo;
    ConvectionDiffusionSettings::Pointer p_settings = r_process_info[CONVECTION_DIFFUSION_SETTINGS];
    auto& r_settings = *p_settings;
    const Variable<double>& r_unknown_var = r_settings.GetUnknownVariable(); // Temperature

    double x_coord_gauss_point = 0;
    double y_coord_gauss_point = 0;
    double rOutput = 0;

    for (IndexType i = 0; i < nb_nodes; ++i)
    {
        // KRATOS_WATCH(r_geometry[i])
        double output_solution_step_value = r_geometry[i].GetSolutionStepValue(r_unknown_var);
        rOutput += r_N(0, i) * output_solution_step_value;
        x_coord_gauss_point += r_N(0, i) * r_geometry[i].X();
        y_coord_gauss_point += r_N(0, i) * r_geometry[i].Y();
    }        

    std::ofstream output_file("output.txt", std::ios::app);
    if (output_file.is_open()) {
        output_file << std::scientific << std::setprecision(14); // Set precision to 10^-14
        output_file << rOutput << " " << r_geometry.Center().X() << " " << r_geometry.Center().Y() << std::endl;
        output_file.close();
    }
    }

template<std::size_t TDim>
void ConvDiffIGAElement<TDim>::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo){}


template class ConvDiffIGAElement<2>;
template class ConvDiffIGAElement<3>;

} // Namespace Kratos