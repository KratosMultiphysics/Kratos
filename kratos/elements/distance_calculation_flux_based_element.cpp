//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Daniel Diez
//                   Pablo Becker
//

#include "utilities/geometry_utilities.h"
#include "includes/kratos_flags.h"
#include "utilities/element_size_calculator.h"
#include "distance_calculation_flux_based_element.h"
#include "includes/checks.h"


namespace Kratos
{

/////specifications must be on top

//gauss points for the 3D4N Tetra
template <>
void DistanceCalculationFluxBasedElement<3,4>::GetSimplexShapeFunctionsOnGauss(BoundedMatrix<double,4, 4>& Ncontainer)
{
    Ncontainer(0,0) = 0.58541020; Ncontainer(0,1) = 0.13819660; Ncontainer(0,2) = 0.13819660; Ncontainer(0,3) = 0.13819660;
    Ncontainer(1,0) = 0.13819660; Ncontainer(1,1) = 0.58541020; Ncontainer(1,2) = 0.13819660; Ncontainer(1,3) = 0.13819660;
    Ncontainer(2,0) = 0.13819660; Ncontainer(2,1) = 0.13819660; Ncontainer(2,2) = 0.58541020; Ncontainer(2,3) = 0.13819660;
    Ncontainer(3,0) = 0.13819660; Ncontainer(3,1) = 0.13819660; Ncontainer(3,2) = 0.13819660; Ncontainer(3,3) = 0.58541020;
}

//gauss points for the 2D3N triangle
template <>
void DistanceCalculationFluxBasedElement<2,3>::GetSimplexShapeFunctionsOnGauss(BoundedMatrix<double,3,3>& Ncontainer)
{
    const double one_sixt = 1.0/6.0;
    const double two_third = 2.0/3.0;
    Ncontainer(0,0) = one_sixt; Ncontainer(0,1) = one_sixt; Ncontainer(0,2) = two_third;
    Ncontainer(1,0) = one_sixt; Ncontainer(1,1) = two_third; Ncontainer(1,2) = one_sixt;
    Ncontainer(2,0) = two_third; Ncontainer(2,1) = one_sixt; Ncontainer(2,2) = one_sixt;
}

//unused, but needed for compilation
template <unsigned int TDim , unsigned int TNumNodes >
void DistanceCalculationFluxBasedElement<TDim,TNumNodes>::GetSimplexShapeFunctionsOnGauss(BoundedMatrix<double,TNumNodes,TNumNodes>& Ncontainer)
{
    KRATOS_ERROR << "CALLING TEMPLATED GetSimplexShapeFunctionsOnGauss" << std::endl;
}


//2d3n triangles
template<>
void DistanceCalculationFluxBasedElement<2,3>::CalculateGaussPointsData(
    const GeometryType& rGeometry,
    BoundedVector<double,3> &rGaussWeights,
    BoundedMatrix<double,3,3> &rNContainer,
    array_1d<BoundedMatrix<double, 3, 2>, 3>& rDN_DXContainer 
    )
{
    BoundedMatrix<double,3,2> DN_DX; // Gradients matrix
    array_1d<double,3> N;            // Position of the gauss point
    double area;
    GeometryUtils::CalculateGeometryData(rGeometry, DN_DX, N, area);
    GetSimplexShapeFunctionsOnGauss(rNContainer);
    for (unsigned int i = 0; i < 3; ++i) {
        rGaussWeights[i] = area/3.0;
        rDN_DXContainer[i] = DN_DX;
    }
}

//3d4n tetras
template<>
void DistanceCalculationFluxBasedElement<3,4>::CalculateGaussPointsData(
    const GeometryType& rGeometry,
    BoundedVector<double,4> &rGaussWeights,
    BoundedMatrix<double,4,4> &rNContainer,
    array_1d<BoundedMatrix<double, 4, 3>, 4>& rDN_DXContainer 
    )
{
    BoundedMatrix<double,4,3> DN_DX; // Gradients matrix
    array_1d<double,4> N;            // Position of the gauss point
    double volume;
    GeometryUtils::CalculateGeometryData(rGeometry, DN_DX, N, volume);
    GetSimplexShapeFunctionsOnGauss(rNContainer);
    for (unsigned int i = 0; i < 4; ++i) {
        rGaussWeights[i] = volume/4.0;
        rDN_DXContainer[i] = DN_DX;
    }

    
}

//generic. lower performance but valid for all geoms.
template <unsigned int TDim , unsigned int TNumNodes >
void DistanceCalculationFluxBasedElement<TDim,TNumNodes>::CalculateGaussPointsData(
    const GeometryType& rGeometry,
    BoundedVector<double,TNumNodes> &rGaussWeights,
    BoundedMatrix<double,TNumNodes,TNumNodes> &rNContainer,
    array_1d<BoundedMatrix<double, TNumNodes, TDim>, TNumNodes>& rDN_DXContainer 
    )
{
    Vector det_j_vector(TNumNodes);
    const auto integration_method = GeometryData::IntegrationMethod::GI_EXTENDED_GAUSS_2;
    rNContainer = rGeometry.ShapeFunctionsValues(integration_method);
    ShapeFunctionsGradientsType DN_DXContainer_unbounded;
    //getting information from geometry in DenseVector<Matrix> format
    rGeometry.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer_unbounded, det_j_vector, integration_method);


    const unsigned int number_of_gauss_points = rGeometry.IntegrationPointsNumber(integration_method);
    const auto& r_integration_points = rGeometry.IntegrationPoints(integration_method);

    for (unsigned int g = 0; g < number_of_gauss_points; ++g){
        rGaussWeights[g] = det_j_vector[g] * r_integration_points[g].Weight();
        rDN_DXContainer[g] = DN_DXContainer_unbounded[g];
    }
}


//////////////////////////////////////////////////////////////////////////////
// Public Operations

template <unsigned int TDim , unsigned int TNumNodes >
void DistanceCalculationFluxBasedElement<TDim, TNumNodes >::CalculateLocalSystem(
    MatrixType &rLeftHandSideMatrix,
    VectorType &rRightHandSideVector,
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    // Resize and intialize output
    if (rLeftHandSideMatrix.size1() != LocalSize || rLeftHandSideMatrix.size2() != LocalSize)
        rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);

    if (rRightHandSideVector.size() != LocalSize)
        rRightHandSideVector.resize(LocalSize, false);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);
    noalias(rRightHandSideVector) = ZeroVector(LocalSize);

    const unsigned int step = rCurrentProcessInfo[FRACTIONAL_STEP];
    
    if(step == 1){ //solve a transcient diffusion problem
        CalculatePotentialFlowSystem(
            rLeftHandSideMatrix, 
            rRightHandSideVector,
            rCurrentProcessInfo);
    }
    else if (step == 2){ // solve convection + source
        CalculateDistanceSystem(
            rLeftHandSideMatrix, 
            rRightHandSideVector,
            rCurrentProcessInfo);
    }
    else{
        KRATOS_ERROR << "FRACTIONAL_STEP must be 1 or 2" << std::endl;
    }

    KRATOS_CATCH("");
}


template <unsigned int TDim , unsigned int TNumNodes >
void DistanceCalculationFluxBasedElement<TDim, TNumNodes >::CalculatePotentialFlowSystem(
        MatrixType &rLeftHandSideMatrix,
        VectorType &rRightHandSideVector,
        const ProcessInfo &rCurrentProcessInfo)
{
    BoundedVector<double,TNumNodes> gauss_weights;
    BoundedMatrix<double,TNumNodes,TNumNodes> N_container;
    array_1d<BoundedMatrix<double, TNumNodes, TDim>, TNumNodes> DN_DX_container;

    CalculateGaussPointsData( GetGeometry(), gauss_weights, N_container, DN_DX_container);

    array_1d<double, TNumNodes > nodal_values;
    for(unsigned int i=0; i<TNumNodes; i++){
        nodal_values[i] = GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
    }

    //computing density so that Fourier number of the complete domain = 1
    //this ensures that all regions notice the effect of the fixed nodes.
    const double domain_length = rCurrentProcessInfo[CHARACTERISTIC_LENGTH];
    const double density  = 1.0/ (domain_length*domain_length);

    
    //Looping Gauss Points
    const std::size_t n_gauss = gauss_weights.size();
    for (unsigned int igauss = 0; igauss<n_gauss; igauss++) {
        array_1d<double, TNumNodes > N = row(N_container, igauss);
        noalias(rLeftHandSideMatrix) += gauss_weights[igauss]*prod(DN_DX_container[igauss],trans(DN_DX_container[igauss]));
        const double mass_factor = density*gauss_weights[igauss];
        const double d_gauss = inner_prod(N,nodal_values);
        const double initial_value = d_gauss > 0.0 ? domain_length : -domain_length ;
        for(unsigned int j=0; j<TNumNodes; j++){
            rLeftHandSideMatrix(j,j)+= mass_factor*N[j];
            rRightHandSideVector[j] += mass_factor*N[j]*initial_value;
        }
    }
    
    //substracting previous solution
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,nodal_values);
}

    
template <unsigned int TDim, unsigned int TNumNodes >
void DistanceCalculationFluxBasedElement<TDim, TNumNodes >::CalculateDistanceSystem(
        MatrixType &rLeftHandSideMatrix,
        VectorType &rRightHandSideVector,
        const ProcessInfo &rCurrentProcessInfo)
{
    BoundedVector<double,TNumNodes> gauss_weights;
    BoundedMatrix<double,TNumNodes,TNumNodes> N_container;
    array_1d<BoundedMatrix<double, TNumNodes, TDim>, TNumNodes> DN_DX_container;

    CalculateGaussPointsData( GetGeometry(), gauss_weights, N_container, DN_DX_container);

    array_1d<double, TNumNodes > nodal_values;
    array_1d< array_1d<double, TDim >, TNumNodes> v; //convection velocity
    array_1d< array_1d<double, TDim >, TNumNodes> v_unit; //to decide if convection must be turned off
    array_1d<double, TDim > avg_v_unit = ZeroVector(TDim); //to decide if convection must be turned off
    bool has_fixed_node = false;

    for(unsigned int i=0; i<TNumNodes; i++){
        nodal_values[i] = GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
        const auto& rVel_i = GetGeometry()[i].GetValue(POTENTIAL_GRADIENT);
        for(unsigned int j=0; j<TDim; j++){
            v[i][j] = rVel_i[j];
        }
        v_unit[i] =  v[i]/norm_2(v[i]); // |v_unit| = 1
        avg_v_unit += v_unit[i];
        has_fixed_node = has_fixed_node || GetGeometry()[i].IsFixed(DISTANCE);
    }
    avg_v_unit/=TNumNodes;

    //computing element size(for Tau)
    const double h = ElementSizeCalculator<TDim,TNumNodes>::AverageElementSize(this->GetGeometry());

    //DECIDE IF CONVECTION IS ACTIVE/INACTIVE
   //we decide by checking if flow converges to /diverges from this element
   //if all velocity vectors are oriented in the same direction, then convection is valid
   //but if the element acts as a "sink" or "source", convection is turned off to avoid instabiliities.
    bool add_convection = true;
    double average_unit_vel_misaligment = 0.0;
    for(unsigned int i=0; i<TNumNodes; i++) {
        average_unit_vel_misaligment += norm_2( v_unit[i]  - avg_v_unit );
    }
    average_unit_vel_misaligment/=TNumNodes;
    if(average_unit_vel_misaligment>0.33) add_convection=false;

    //checking if the flow comes from face that has neighbour element (only if it is not inlet)
    if( !has_fixed_node )
    {   
        double value_lowest_prod = 0.0;
        unsigned int face_lowest_prod = 0;
        const auto elem_boundaries = GetGeometry().LocalSpaceDimension()==3 ? GetGeometry().GenerateFaces() : GetGeometry().GenerateEdges() ;

        for(unsigned int i=0; i<elem_boundaries.size(); i++) {
            const auto cond_center = elem_boundaries[i].Center();
            const array_1d<double,3> unit_normal = elem_boundaries[i].UnitNormal(cond_center);
            const double face_result = inner_prod(avg_v_unit,unit_normal);
            if(face_result<value_lowest_prod){
                value_lowest_prod = face_result;
                face_lowest_prod=i;
            }
        }
        if( this->GetValue(NEIGHBOUR_ELEMENTS)(face_lowest_prod).get()==nullptr){
            add_convection=false;
        }
    }

    //TERMS TO BE ASSEMBLED
    //terms which multiply the gradient of phi
    BoundedMatrix<double,TNumNodes, TNumNodes> aux2 = ZeroMatrix(TNumNodes, TNumNodes); //terms multiplying phi
    
    //source term
    array_1d<double, TNumNodes> rhs_volumetric_heat = ZeroVector(TNumNodes);

    //Gauss point container (using nodal integration)
    const bounded_matrix<double, TNumNodes, TNumNodes> Ncontainer=IdentityMatrix(TNumNodes); 

    //Looping Gauss Points
    const std::size_t n_gauss = gauss_weights.size();
    for (unsigned int igauss = 0; igauss<n_gauss; igauss++) {
        //Getting the correct GP
        array_1d<double, TNumNodes > N = row(Ncontainer, igauss);

        //Velocity in GP
        array_1d<double, TDim > vel_gauss = ZeroVector(TDim);
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            for (unsigned int k = 0; k<TDim; k++)
                vel_gauss[k] += N[i] * v[i][k];
        }

        const double d_gauss = inner_prod(N,nodal_values);

        if (d_gauss <= 0.0) {
            vel_gauss = -vel_gauss;
        }

        const double norm_vel = norm_2(vel_gauss);

        array_1d<double, TNumNodes > a_dot_grad = prod(DN_DX_container[igauss], vel_gauss);

        const double tau_denom = std::max(2.0 * norm_vel / h ,  1e-3); 
        const double tau = 1.0 / (tau_denom);

        const double source_term = d_gauss > 0 ? norm_vel : -norm_vel;

        if(add_convection){
            //convection + convection stabilization
            noalias(aux2) += outer_prod(N, a_dot_grad)*gauss_weights[igauss];
            noalias(aux2) += tau*outer_prod(a_dot_grad, a_dot_grad)*gauss_weights[igauss];

            //source
            rhs_volumetric_heat += source_term*N*gauss_weights[igauss]; 
            rhs_volumetric_heat += (tau*source_term)*a_dot_grad*gauss_weights[igauss];
        } else {
            noalias(aux2) += tau*prod(DN_DX_container[igauss],trans(DN_DX_container[igauss]))*gauss_weights[igauss];
        }

        //spherical diffusion to enhance stabilization
        double spherical_diffusion = 0.1*tau*(norm_vel*norm_vel);
        noalias(aux2) += spherical_diffusion*prod(DN_DX_container[igauss],trans(DN_DX_container[igauss]))*gauss_weights[igauss];
    }

    //adding to system LHS
    noalias(rLeftHandSideMatrix) += aux2;

    //adding to system RHS
    noalias(rRightHandSideVector) += rhs_volumetric_heat; //external forces


    //take out the dirichlet part to finish computing the residual
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, nodal_values);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Inquiry

template <unsigned int TDim , unsigned int TNumNodes >
int DistanceCalculationFluxBasedElement<TDim, TNumNodes >::Check(const ProcessInfo &rCurrentProcessInfo) const
{
    KRATOS_TRY;
    // Generic geometry check
    int out = Element::Check(rCurrentProcessInfo);
    if (out != 0) {
        return out;
    }
    
    KRATOS_ERROR_IF( !rCurrentProcessInfo.Has(CHARACTERISTIC_LENGTH) || rCurrentProcessInfo[CHARACTERISTIC_LENGTH]<std::numeric_limits<double>::epsilon() ) << "CHARACTERISTIC_LENGTH is zero or undefined" << std::endl;

    const auto& r_geometry = this->GetGeometry();

    const unsigned int n_expected_neighbours = r_geometry.LocalSpaceDimension()==3 ? r_geometry.FacesNumber() : r_geometry.EdgesNumber() ;
    KRATOS_ERROR_IF( !this->Has(NEIGHBOUR_ELEMENTS) || this->GetValue(NEIGHBOUR_ELEMENTS).size()!=n_expected_neighbours ) << "Neighbour elements not defined" << std::endl;

    for(unsigned int i=0; i<TNumNodes; ++i){
        const Node<3>& rNode = r_geometry[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISTANCE,rNode);
        KRATOS_CHECK_DOF_IN_NODE(DISTANCE,rNode);
    }

    // If this is a 2D problem, check that nodes are in XY plane
    if ( TDim == 2){
        for (unsigned int i=0; i<TNumNodes; ++i) {
            if (std::abs(r_geometry[i].Z())>1e-9)
                KRATOS_ERROR << "Node " << r_geometry[i].Id() << "has non-zero Z coordinate." << std::endl;
        }
    }

    return out;

    KRATOS_CATCH("");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public I/O

template <unsigned int TDim , unsigned int TNumNodes >
std::string DistanceCalculationFluxBasedElement<TDim, TNumNodes >::Info() const
{
    std::stringstream buffer;
    buffer << "DistanceCalculationFluxBasedElement" << TDim << "D" << TNumNodes << "N #" << this->Id();
    return buffer.str();
}

template <unsigned int TDim, unsigned int TNumNodes >
void DistanceCalculationFluxBasedElement<TDim,TNumNodes>::PrintInfo(
    std::ostream &rOStream) const
{
    rOStream << this->Info() << std::endl;
}

template <unsigned int TDim , unsigned int TNumNodes >
void DistanceCalculationFluxBasedElement< TDim, TNumNodes  >::EquationIdVector(EquationIdVectorType &rResult, const ProcessInfo &rCurrentProcessInfo) const
{
    const auto& r_geometry = this->GetGeometry();

    unsigned int LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    for (unsigned int i = 0; i < TNumNodes; ++i){
        rResult[LocalIndex++] = r_geometry[i].GetDof(DISTANCE).EquationId();
    }
}


template <unsigned int TDim , unsigned int TNumNodes >
void DistanceCalculationFluxBasedElement< TDim , TNumNodes >::GetDofList(DofsVectorType &rElementalDofList, const ProcessInfo &rCurrentProcessInfo) const
{
    const auto& r_geometry = this->GetGeometry();

     if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

     unsigned int LocalIndex = 0;

     for (unsigned int i = 0; i < TNumNodes; ++i){
        rElementalDofList[LocalIndex++] = r_geometry[i].pGetDof(DISTANCE);
     }
}


template <unsigned int TDim , unsigned int TNumNodes >
void DistanceCalculationFluxBasedElement< TDim, TNumNodes >::AddExplicitContribution(const ProcessInfo& rCurrentProcessInfo)
{
    auto& rGeometry = this->GetGeometry();

    BoundedVector<double,TNumNodes> gauss_weights;
    BoundedMatrix<double,TNumNodes,TNumNodes> N_container;
    array_1d<BoundedMatrix<double, TNumNodes, TDim>, TNumNodes> DN_DX_container;

    CalculateGaussPointsData( rGeometry , gauss_weights, N_container, DN_DX_container);


    //get the nodal values
    array_1d<double, TNumNodes > nodal_values;
    for(unsigned int i=0; i<TNumNodes; i++){
        nodal_values[i] = GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
    }

    BoundedMatrix<double, TNumNodes, TDim> avg_DN_DX = ZeroMatrix(TNumNodes, TDim);
    for(unsigned int i=0; i<TNumNodes; i++){
        avg_DN_DX += DN_DX_container[i] /TNumNodes;
    }
    
    const array_1d<double, TDim> avg_grad = prod(trans(avg_DN_DX), nodal_values);

    //saving data
    const double vol_factor =  rGeometry.DomainSize()/TNumNodes;
    for (unsigned int j = 0; j < TNumNodes; j++){ //looping 4 nodes of the elem:
        rGeometry[j].SetLock();
        rGeometry[j].GetValue(NODAL_VOLUME)+=vol_factor;
        rGeometry[j].GetValue(POTENTIAL_GRADIENT)+= avg_grad*vol_factor;      
        rGeometry[j].UnSetLock();
    }
}


// Class template instantiation
template class  DistanceCalculationFluxBasedElement<2,3>;
template class  DistanceCalculationFluxBasedElement<3,4>;
template class  DistanceCalculationFluxBasedElement<3,8>;


} // namespace Kratos