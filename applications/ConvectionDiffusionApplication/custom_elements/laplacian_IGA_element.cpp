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
#include "custom_elements/laplacian_IGA_element.h"

#include "utilities/math_utils.h"

namespace Kratos
{

template<std::size_t TDim>
LaplacianIGAElement<TDim>::LaplacianIGAElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : Element(
        NewId,
        pGeometry)
{
}

template<std::size_t TDim>
LaplacianIGAElement<TDim>::LaplacianIGAElement(
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
Element::Pointer LaplacianIGAElement<TDim>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<LaplacianIGAElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

template<std::size_t TDim>
Element::Pointer LaplacianIGAElement<TDim>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<LaplacianIGAElement>(NewId, pGeom, pProperties);
}

// Deconstructor
template<std::size_t TDim>
LaplacianIGAElement<TDim>::~LaplacianIGAElement()
{
}


// From classical Laplacian
template<std::size_t TDim>
void LaplacianIGAElement<TDim>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
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
    const unsigned int number_of_points = r_geometry.size();
    const unsigned int dim = r_geometry.WorkingSpaceDimension();
    // const unsigned int dim = 2;

    // KRATOS_WATCH(number_of_points) // number_of_points = 9
    
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

    // const Matrix& r_DN_De = r_geometry.ShapeFunctionDerivatives(1, 0);
    // KRATOS_WATCH(r_DN_De)
    
    // KRATOS_WATCH(integration_points)
    // KRATOS_WATCH(DN_De)
    
    Element::GeometryType::JacobiansType J0;
    Matrix DN_DX(number_of_points,dim);
    Matrix InvJ0(dim,dim);
    Vector temp(number_of_points);

    Vector heat_flux_local(number_of_points);
    Vector nodal_conductivity(number_of_points);
    for(unsigned int node_element = 0; node_element<number_of_points; node_element++)
    {
        heat_flux_local[node_element] = r_geometry[node_element].FastGetSolutionStepValue(r_volume_source_var);
        nodal_conductivity[node_element] = r_geometry[node_element].FastGetSolutionStepValue(r_diffusivity_var);
    }

    r_geometry.Jacobian(J0,this->GetIntegrationMethod());
    double DetJ0;
    // KRATOS_WATCH(integration_points.size())
    // exit(0) ;

    for(std::size_t i_point = 0; i_point < integration_points.size(); ++i_point)
    {
        // Calculating inverse jacobian and jacobian determinant
        MathUtils<double>::InvertMatrix(J0[i_point],InvJ0,DetJ0);
        // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
        noalias(DN_DX) = prod(DN_De[i_point],InvJ0);

        // NEW ! - 2D 
        for (size_t i = 0; i < DN_DX.size1(); ++i) {
            DN_DX(i, 2) = 0.0;
        }

        auto N = row(N_gausspoint,i_point); //these are the N which correspond to the gauss point "i_point"
        const double IntToReferenceWeight = integration_points[i_point].Weight() * DetJ0;
        const double conductivity_gauss = inner_prod(N, nodal_conductivity);
        noalias(rLeftHandSideMatrix) += IntToReferenceWeight * conductivity_gauss * prod(DN_DX, trans(DN_DX)); //

        // Calculating the local RHS
        const double qgauss = inner_prod(N, heat_flux_local);

        noalias(rRightHandSideVector) += IntToReferenceWeight*qgauss*N;
    }

    // exit(0); 

    // RHS = ExtForces - K*temp;
    for (unsigned int i = 0; i < number_of_points; i++)
        temp[i] = r_geometry[i].GetSolutionStepValue(r_unknown_var);
    
    // RHS -= K*temp
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp);

    // KRATOS_WATCH(rRightHandSideVector)
    // KRATOS_WATCH(rLeftHandSideMatrix)
    // exit(0);
    KRATOS_CATCH("")
}


// From classical Laplacian
template<std::size_t TDim>
void LaplacianIGAElement<TDim>::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{   
    KRATOS_WATCH('LHS111')
    VectorType temp(0);
    CalculateLocalSystem(rLeftHandSideMatrix, temp, rCurrentProcessInfo);
}




// From classical Laplacian
template<std::size_t TDim>
void LaplacianIGAElement<TDim>::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    // KRATOS_WATCH('RHS111')
    MatrixType temp(0,0);
    CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
}





template<std::size_t TDim>
void LaplacianIGAElement<TDim>::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const
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
void LaplacianIGAElement<TDim>::GetDofList(DofsVectorType& ElementalDofList,const ProcessInfo& rCurrentProcessInfo) const
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
int LaplacianIGAElement<TDim>::Check(const ProcessInfo& rCurrentProcessInfo) const
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
Element::IntegrationMethod LaplacianIGAElement<TDim>::GetIntegrationMethod() const
{
    return GeometryData::IntegrationMethod::GI_GAUSS_1;
}



template<std::size_t TDim>
void LaplacianIGAElement<TDim>::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
    {
    const auto& r_geometry = GetGeometry();
    const SizeType nb_nodes = r_geometry.size();

    // Integration Points
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints();
    // Shape function values
    const Matrix& r_N = r_geometry.ShapeFunctionsValues();

    // KRATOS_WATCH(r_N)  
    // const auto& integration_method = r_geometry.GetDefaultIntegrationMethod();
    // KRATOS_WATCH(r_geometry)
    // KRATOS_WATCH(this->GetIntegrationMethod())
    // const Matrix& r_N2 = r_geometry.ShapeFunctionsValues(r_geometry[0]) ;
    // KRATOS_WATCH(r_N2)

    // exit(0) ;

    std::vector<double> rOutput;

    const ProcessInfo& r_process_info = rCurrentProcessInfo;
    ConvectionDiffusionSettings::Pointer p_settings = r_process_info[CONVECTION_DIFFUSION_SETTINGS];
    auto& r_settings = *p_settings;
    const Variable<double>& r_unknown_var = r_settings.GetUnknownVariable(); // Temperature

    if (rOutput.size() != integration_points.size())
        rOutput.resize(integration_points.size());

    // KRATOS_WATCH(integration_points.size())

    for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number)
    {
        rOutput[point_number] = 0.0;
        for (IndexType i = 0; i < nb_nodes; ++i)
        {
            // KRATOS_WATCH(r_geometry[i])
            double output_solution_step_value = r_geometry[i].FastGetSolutionStepValue(r_unknown_var);
            rOutput[point_number] += r_N(point_number, i) * output_solution_step_value;
            // KRATOS_WATCH(r_N(point_number, i))
            // KRATOS_WATCH(r_geometry.GetPoint(i)) // works
           
        }
        // exit(0) ;
        // KRATOS_WATCH(rOutput[point_number])
        
    }
    // double output = rOutput[0] ;
    // KRATOS_WATCH(rOutput)
    // KRATOS_WATCH(r_geometry.Center())
    // std::ofstream output_file("output.txt", std::ios::app);
    // if (output_file.is_open()) {
    //     output_file << output << " " << r_geometry.Center().X() << " " << r_geometry.Center().Y() << std::endl;
    //     output_file.close();
    // }
    }


template class LaplacianIGAElement<2>;
template class LaplacianIGAElement<3>;

} // Namespace Kratos