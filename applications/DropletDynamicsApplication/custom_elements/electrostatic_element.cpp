//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Christian Narváez-Muñoz
//

// System includes


// External includes
#include "utilities/element_size_calculator.h"

// Include Base h
#include "includes/checks.h"
#include "includes/variables.h"
#include "utilities/geometry_utilities.h"
#include "custom_elements/electrostatic_element.h"
#include "utilities/math_utils.h"
#include "utilities/enrichment_utilities.h"
#include "utilities/enrich_2d_2dofs.h"


namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * Constructor.
 */
template< unsigned int TDim >
ElectrostaticElement<TDim>::ElectrostaticElement(IndexType NewId)
    : Element(NewId) 
{
}

/**
 * Constructor using an array of nodes
 */
template< unsigned int TDim >
ElectrostaticElement<TDim>::ElectrostaticElement(IndexType NewId, const NodesArrayType& ThisNodes)
    : Element(NewId, ThisNodes) 
{
}

/**
 * Constructor using Geometry
 */
template< unsigned int TDim >
ElectrostaticElement<TDim>::ElectrostaticElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry) 
{
}

/**
 * Constructor using Properties
 */
template< unsigned int TDim >
ElectrostaticElement<TDim>::ElectrostaticElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties) 
{
}

/**
 * Copy Constructor
 */
template< unsigned int TDim >
ElectrostaticElement<TDim>::ElectrostaticElement(ElectrostaticElement const& rOther)
    : Element(rOther) 
{
}

/**
 * Destructor
 */
template< unsigned int TDim >
ElectrostaticElement<TDim>::~ElectrostaticElement()
{
}

///@}
///@name Operators
///@{

/// Assignment operator.
template< unsigned int TDim >
ElectrostaticElement<TDim> & ElectrostaticElement<TDim>::operator=(ElectrostaticElement const& rOther)
{
    BaseType::operator=(rOther);
    Flags::operator =(rOther);
    // mpProperties = rOther.mpProperties;
    return *this;
}

///@}
///@name Operations
///@{

/**
 * ELEMENTS inherited from this class have to implement next
 * Create and Clone methods: MANDATORY
 */

/**
 * creates a new element pointer
 * @param NewId: the ID of the new element
 * @param ThisNodes: the nodes of the new element
 * @param pProperties: the properties assigned to the new element
 * @return a Pointer to the new element
 */
template< unsigned int TDim >
Element::Pointer ElectrostaticElement<TDim>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<ElectrostaticElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

/**
 * creates a new element pointer
 * @param NewId: the ID of the new element
 * @param pGeom: the geometry to be employed
 * @param pProperties: the properties assigned to the new element
 * @return a Pointer to the new element
 */
template< unsigned int TDim >
Element::Pointer ElectrostaticElement<TDim>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<ElectrostaticElement>(NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

/**
 * creates a new element pointer and clones the previous element data
 * @param NewId: the ID of the new element
 * @param ThisNodes: the nodes of the new element
 * @param pProperties: the properties assigned to the new element
 * @return a Pointer to the new element
 */
template< unsigned int TDim >
Element::Pointer ElectrostaticElement<TDim>::Clone(IndexType NewId, NodesArrayType const& ThisNodes) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<ElectrostaticElement>(NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    KRATOS_CATCH("");
}

/**
 * this determines the elemental equation ID vector for all elemental
 * DOFs
 * @param rResult: the elemental equation ID vector
 * @param rCurrentProcessInfo: the current process info instance
 */
template< unsigned int TDim >
void ElectrostaticElement<TDim>::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& CurrentProcessInfo) const
{
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    if (rResult.size() != number_of_nodes)
        rResult.resize(number_of_nodes, false);

    for (unsigned int i = 0; i < number_of_nodes; i++)
        rResult[i] = GetGeometry()[i].GetDof(EPOTENTIAL).EquationId();


}

/**
 * determines the elemental list of DOFs
 * @param ElementalDofList: the list of DOFs
 * @param rCurrentProcessInfo: the current process info instance
 */
template< unsigned int TDim >
void ElectrostaticElement<TDim>::GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& CurrentProcessInfo) const
{
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    if (rElementalDofList.size() != number_of_nodes)
        rElementalDofList.resize(number_of_nodes);

      for (unsigned int i = 0; i < number_of_nodes; i++)
        rElementalDofList[i] = GetGeometry()[i].pGetDof(EPOTENTIAL);


}

/**
 * ELEMENTS inherited from this class have to implement next
 * CalculateLocalSystem, CalculateLeftHandSide and CalculateRightHandSide methods
 * they can be managed internally with a private method to do the same calculations
 * only once: MANDATORY
 */

/**
 * this is called during the assembling process in order
 * to calculate all elemental contributions to the global system
 * matrix and the right hand side
 * @param rLeftHandSideMatrix: the elemental left hand side matrix
 * @param rRightHandSideVector: the elemental right hand side
 * @param rCurrentProcessInfo: the current process info instance
 */
template< >
void ElectrostaticElement<2>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    /// general sets  /////
    const int num_dim  = 2;
    const int num_nodes  = num_dim + 1;
    const int num_dof = num_nodes; 
    const int num_enrch = num_nodes+1; 
    const int num_edges = 3; 
    GeometryData::ShapeFunctionsGradientsType DN_DX;
    Matrix N;
    Vector DetJ;
    Vector weights;
    
    //////////////////////////////////////////////////////////////////////////
    ///// matrix /////////////////////////////
    Vector distance(num_dof);
    Vector distances(num_dof);   
    Vector values(num_dof);
    BoundedMatrix<double,num_dof,num_dof> lhs = ZeroMatrix(num_dof,num_dof);
    Vector rhs = ZeroVector(num_dof);
    //array_1d<double,3> surface_sources;
    /////////////////////////////////////////////////////////////////////////////////
    const double permittivity_neg = this->GetProperties().GetValue(PERMITTIVITYNEG);
    const double permittivity_pos = this->GetProperties().GetValue(PERMITTIVITYPOS);
    const double conductivity_neg = this->GetProperties().GetValue(CONDUCTIVITYNEG);
    const double conductivity_pos = this->GetProperties().GetValue(CONDUCTIVITYPOS);
    ///////////////////////////////////////////////////////////////////////////////////
    const GeometryData::IntegrationMethod integration_method = GeometryData::IntegrationMethod::GI_GAUSS_1;
    GeometryType::Pointer p_geometry = this->pGetGeometry();
    const unsigned int number_of_gauss_points = p_geometry->IntegrationPointsNumber(integration_method);
    p_geometry->ShapeFunctionsIntegrationPointsGradients(DN_DX, DetJ, integration_method);
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    if(rLeftHandSideMatrix.size1() != num_dof)
        rLeftHandSideMatrix.resize(num_dof,num_dof,false); //resizing the system in case it does not have the right size 

    if(rRightHandSideVector.size() != num_dof)
        rRightHandSideVector.resize(num_dof,false);

    for (unsigned int i_node=0; i_node < num_nodes; ++i_node){
        const double dist = (*p_geometry)[i_node].FastGetSolutionStepValue(DISTANCE);
        if (dist == 0.0){                   
            (*p_geometry)[i_node].FastGetSolutionStepValue(DISTANCE) = 1.0e-10;
        }
    }
    
    unsigned int nneg=0, npos=0;

    for (unsigned int i_node=0; i_node < num_nodes; ++i_node){
        const double dist = (*p_geometry)[i_node].FastGetSolutionStepValue(DISTANCE);
        distance[i_node] = dist;
        values[i_node] = (*p_geometry)[i_node].FastGetSolutionStepValue(EPOTENTIAL);
        
        if (dist > 0.0){
            npos += 1;
        }
        else{
            nneg += 1;
        }
    }
    
    if (nneg == num_nodes || npos == num_nodes){
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = p_geometry->IntegrationPoints(integration_method);
        if (weights.size() != number_of_gauss_points) {
            weights.resize(number_of_gauss_points,false);
        }
        //calculating the area of the elementes
        for (unsigned int gp = 0; gp < number_of_gauss_points; gp++){
            weights[gp] = DetJ[gp] * IntegrationPoints[gp].Weight();    
           
            if (npos == 3){ //same results as enriched_n_element code
                lhs = -conductivity_pos*weights[gp]*prod(DN_DX[gp],trans(DN_DX[gp])); 
                //lhs = -permittivity_pos*weights[gp]*prod(DN_DX[gp],trans(DN_DX[gp])); 
            }else {
                lhs = -conductivity_neg*weights[gp]*prod(DN_DX[gp],trans(DN_DX[gp]));
                //lhs = -permittivity_neg*weights[gp]*prod(DN_DX[gp],trans(DN_DX[gp]));
            }
        }

        noalias(rLeftHandSideMatrix) = lhs;
        noalias(rRightHandSideVector) = rhs - prod(rLeftHandSideMatrix,values);  
    } else {  
        //////////////////////////////////////////////////////////////////////////////////////////

        Vector cut_edges = ZeroVector(3);
        // Implementing the boundary condition on distance gradient
        
        const auto& neighbour_elems = this->GetValue(NEIGHBOUR_ELEMENTS);

        for (unsigned int i_ne = 0; i_ne < 3; i_ne++){
                cut_edges(i_ne)=i_ne;
        }

        ///////////////////////// SPLITTED SHAPE FUNCTIONS  //////////////////////////////////////  
        ModifiedShapeFunctions::Pointer p_modified_sh_func =
            Kratos::make_shared<Triangle2D3ModifiedShapeFunctions>(p_geometry, distance);

        Matrix neg_N, pos_N;
        GeometryType::ShapeFunctionsGradientsType neg_DN_DX, pos_DN_DX;
        Vector neg_weights, pos_weights;

        Matrix int_N;
        GeometryType::ShapeFunctionsGradientsType int_DN_DX;
        Vector int_weights; 

        ModifiedShapeFunctions::AreaNormalsContainerType int_normal, pos_extr_normal, neg_extr_normal;

        //ModifiedShapeFunctions::ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues(
        Matrix pos_extr_N, neg_extr_N ; 
        GeometryType::ShapeFunctionsGradientsType pos_extr_DN_DX, neg_extr_DN_DX;
        Vector pos_extr_weights, neg_extr_weights; 
        


        ///////////////////////// SPLITTED SHAPE FUNCTIONS  //////////////////////////////////////  
        p_modified_sh_func->ComputePositiveSideShapeFunctionsAndGradientsValues(
            pos_N,        // N
            pos_DN_DX,    // DN_DX
            pos_weights,  // weight * detJ
            integration_method);

        p_modified_sh_func->ComputeNegativeSideShapeFunctionsAndGradientsValues(
            neg_N,        // N
            neg_DN_DX,    // DN_DX
            neg_weights,  // weight * detJ
            integration_method);
            //KRATOS_WATCH(neg_weights)       
        ////////////////////////////////////////////////////////////////////////////////////////////////
        p_modified_sh_func->ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
            int_N,
            int_DN_DX,
            int_weights,
            integration_method);  

        p_modified_sh_func->ComputeNegativeSideInterfaceAreaNormals(
            int_normal,  // it is not normalized 
            integration_method);
        /////////////////////////////////////////////////////////////////////////////////////////////////
        //NORMAL CONTRIBUTION 
        BoundedMatrix<double, 1, 1 > K_enrich;
        noalias(K_enrich) = ZeroMatrix(1,1);
        BoundedMatrix<double, 1, 3 > B_enrich;
        noalias(B_enrich) = ZeroMatrix(1,3);
        //////////////////////////////
        // FLUX CONTRIBUTION 
        BoundedMatrix<double, 1, 1 > K_flux_enrich;
        noalias(K_flux_enrich) = ZeroMatrix(1,1);
        BoundedMatrix<double, 1, 3 > Bij_flux_enrich;
        noalias(Bij_flux_enrich) = ZeroMatrix(1,3);
        BoundedMatrix<double, 1, 3 > Bji_flux_enrich;
        noalias(Bji_flux_enrich) = ZeroMatrix(1,3);
        /////////////////////////////////////////////////////////////////////////////////////////////////
        for (unsigned int n = 0 ; n < cut_edges.size() ; n++ ){
            if (cut_edges[n] < 3)
            {
                unsigned int cut_edge = cut_edges[n];
                unsigned int pos_FaceId = cut_edge;
                unsigned int neg_FaceId = cut_edge; 
                unsigned int pos_normal_FaceId = cut_edge;
                unsigned int neg_normal_FaceId = cut_edge;
                ////////////////////////////////////////////////////////////////////////////////////////////////////////////
                p_modified_sh_func->ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues(
                    pos_extr_N,
                    pos_extr_DN_DX,
                    pos_extr_weights,
                    pos_FaceId,
                    integration_method);  

                p_modified_sh_func->ComputePositiveExteriorFaceAreaNormals(
                    pos_extr_normal,  // it is not normalized 
                    pos_normal_FaceId,
                    integration_method);
                
                for (unsigned int gp_epos = 0; gp_epos < pos_extr_weights.size() ; gp_epos++){
                    double term2 = 0;
                    Vector grad_term1 = ZeroVector(num_dim);
                    Vector grad_term2 = ZeroVector(num_dim);  
                    Vector n_dn_dx = ZeroVector (3);
                    Vector n_dn_dx_enrch = ZeroVector(1);
                    double Nenr_1 = 0;
                    double Nenr_term2 = 0;

                    //// normalized 
                    Vector normal_pos_side = ZeroVector(num_dim);
                    normal_pos_side = normal_pos_side + pos_extr_normal[gp_epos];
                    normal_pos_side = (1.0 / norm_2(normal_pos_side)) * normal_pos_side;
                                
                    for (unsigned int i = 0; i < num_dof; i++){   
                        term2 += pos_extr_N(gp_epos,i)*distance(i);          
                        for (unsigned int dim = 0; dim < num_dim; dim++){
                            grad_term1[dim] += pos_extr_DN_DX[0](i,dim)*std::abs(distance(i));
                            grad_term2[dim] += pos_extr_DN_DX[0](i,dim)*distance(i);
                        }    
                        Nenr_1 += pos_extr_N(gp_epos,i)*std::abs(distance(i)); // first term of the enrichment function 
                        Nenr_term2 += pos_extr_N(gp_epos,i)*(distance(i)); // second term of the enrichment function 
                    }  
        
                    double sign = 1.0;
                    if (term2 < 0.0) sign = -1.0;
                    for (unsigned int dim = 0; dim < num_dim; ++dim){
                        n_dn_dx_enrch[0] += normal_pos_side(dim)*(grad_term1[dim] - sign*grad_term2[dim]); // added n * grad dn_dx enriched
                    }

                    n_dn_dx = prod(normal_pos_side, trans(pos_extr_DN_DX[0])); // vector (n).(DN_DX)
                    const double pos_Nenr = Nenr_1 - std::abs(Nenr_term2);
                    
                    // components for codensation
                    for (unsigned int i = 0; i < num_dof; ++i){
                        Bij_flux_enrich(0,i) -= conductivity_pos*pos_extr_weights[gp_epos]* pos_Nenr * n_dn_dx(i); // adde
                        Bji_flux_enrich(0,i) -= conductivity_pos*pos_extr_weights[gp_epos]* pos_extr_N(gp_epos,i) * n_dn_dx_enrch[0]; // added  
                        //Bij_flux_enrich(0,i) -= permittivity_pos*pos_extr_weights[gp_epos]* pos_Nenr * n_dn_dx(i); // adde
                        //Bji_flux_enrich(0,i) -= permittivity_pos*pos_extr_weights[gp_epos]* pos_extr_N(gp_epos,i) * n_dn_dx_enrch[0]; // added  
                    }
                
                    K_flux_enrich(0,0) -= conductivity_pos*pos_extr_weights[gp_epos]*pos_Nenr*n_dn_dx_enrch[0];
                    //K_flux_enrich(0,0) -= permittivity_pos*pos_extr_weights[gp_epos]*pos_Nenr*n_dn_dx_enrch[0];
                
                }
                
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
                p_modified_sh_func->ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues(
                    neg_extr_N,
                    neg_extr_DN_DX,
                    neg_extr_weights,
                    neg_FaceId,
                    integration_method);   

                p_modified_sh_func->ComputeNegativeExteriorFaceAreaNormals(
                    neg_extr_normal,  // it is not normalized 
                    neg_normal_FaceId,
                    integration_method);
                ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //// NEGATIVE SIDE
                for (unsigned int gp_eneg = 0; gp_eneg < neg_extr_weights.size() ; gp_eneg++){
                    double term2 = 0;
                    Vector grad_term1 = ZeroVector(num_dim);
                    Vector grad_term2 = ZeroVector(num_dim);  
                    Vector n_dn_dx = ZeroVector (3);
                    Vector n_dn_dx_enrch = ZeroVector(1);   
                    double Nenr_1 = 0;
                    double Nenr_term2 = 0;

                    //// normalized 
                    Vector normal_neg_side = ZeroVector(num_dim);
                    normal_neg_side = normal_neg_side + neg_extr_normal[gp_eneg];
                    normal_neg_side = (1.0 / norm_2(normal_neg_side)) * normal_neg_side;
                
                    for (unsigned int i = 0; i < num_dof; i++){   
                        term2 += neg_extr_N(gp_eneg,i)*distance(i);          
                        for (unsigned int dim = 0; dim < num_dim; dim++){
                            grad_term1[dim] += neg_extr_DN_DX[0](i,dim)*std::abs(distance(i));
                            grad_term2[dim] += neg_extr_DN_DX[0](i,dim)*distance(i);
                        }    
                        Nenr_1 += neg_extr_N(gp_eneg,i)*std::abs(distance(i)); // first term of the enrichment function 
                        Nenr_term2 += neg_extr_N(gp_eneg,i)*(distance(i)); // second term of the enrichment function 
                    }  
                
                    double sign = 1.0;
                    if (term2 < 0.0) sign = -1.0;

                    for (unsigned int dim = 0; dim < num_dim; ++dim){
                        n_dn_dx_enrch[0] += normal_neg_side(dim)*(grad_term1[dim] - sign*grad_term2[dim]); // added n * grad dn_dx enriched
                    }

                    n_dn_dx = prod(normal_neg_side, trans(neg_extr_DN_DX[0])); // vector (n).(DN_DX)
                    const double neg_Nenr = Nenr_1 - std::abs(Nenr_term2);

                    // components for codensation
                    for (unsigned int i = 0; i < num_dof; ++i){
                        Bij_flux_enrich(0,i) -= conductivity_neg*neg_extr_weights[gp_eneg]* neg_Nenr * n_dn_dx(i); // adde
                        Bji_flux_enrich(0,i) -= conductivity_neg*neg_extr_weights[gp_eneg]* neg_extr_N(gp_eneg,i) * n_dn_dx_enrch[0]; // added  
                        //Bij_flux_enrich(0,i) -= permittivity_neg*neg_extr_weights[gp_eneg]* neg_Nenr * n_dn_dx(i); // adde
                        //Bji_flux_enrich(0,i) -= permittivity_neg*neg_extr_weights[gp_eneg]* neg_extr_N(gp_eneg,i) * n_dn_dx_enrch[0]; // added  
                    }
            
                    K_flux_enrich(0,0) -= conductivity_neg*neg_extr_weights[gp_eneg]*neg_Nenr*n_dn_dx_enrch[0];
                    //K_flux_enrich(0,0) -= permittivity_neg*neg_extr_weights[gp_eneg]*neg_Nenr*n_dn_dx_enrch[0];
                }
                //KRATOS_WATCH("end flux_neg contribution")
            }
        }

        /////////////////////////////////////////////////////////////////////////////////////////////////
        const int num_pos_gp = pos_weights.size();
        const int num_neg_gp = neg_weights.size();
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //// POSITIVE SIDE 
        Vector pos_Nenr = ZeroVector(pos_weights.size());
        BoundedMatrix<double, 2, 2 > pos_grad_enrich;
        noalias(pos_grad_enrich) = ZeroMatrix(1,1);

        for (unsigned int pos_e = 0; pos_e < num_pos_gp; pos_e++){   
            lhs -= conductivity_pos*pos_weights(pos_e)*prod(pos_DN_DX[pos_e],trans(pos_DN_DX[pos_e]));
            //lhs -= permittivity_pos*pos_weights(pos_e)*prod(pos_DN_DX[pos_e],trans(pos_DN_DX[pos_e]));

            double term2 = 0;
            Vector grad_term1 = ZeroVector(num_dim);
            Vector grad_term2 = ZeroVector(num_dim);  
            double Nenr_1 = 0;
            double Nenr_term2 = 0;

            for (unsigned int i = 0; i < num_dof; i++){   
                term2 += pos_N(pos_e,i)*distance(i);          
                for (unsigned int dim = 0; dim < num_dim; dim++){
                    grad_term1[dim] += pos_DN_DX[pos_e](i,dim)*std::abs(distance(i));
                    grad_term2[dim] += pos_DN_DX[pos_e](i,dim)*distance(i);
                }    
            }  

          
            // components for codensation
            for (unsigned int dim = 0; dim < num_dim; ++dim){
                double sign = 1.0;
                if (term2 < 0.0) sign = -1.0;
                const double enrich_DN_DX_gp_dim = sign*grad_term2[dim];
                for (unsigned int i = 0; i < num_dof; ++i){
                    B_enrich(0, i) -= conductivity_pos*pos_weights(pos_e)
                        *pos_DN_DX[pos_e](i,dim)*( grad_term1[dim] - enrich_DN_DX_gp_dim );
/*                     B_enrich(0, i) -= permittivity_pos*pos_weights(pos_e)
                        *pos_DN_DX[pos_e](i,dim)*( grad_term1[dim] - enrich_DN_DX_gp_dim ); */
                    
                }
                K_enrich(0,0) -= conductivity_pos*pos_weights(pos_e)
                    *( grad_term1[dim] - enrich_DN_DX_gp_dim )*( grad_term1[dim] - enrich_DN_DX_gp_dim );
                pos_grad_enrich (pos_e,dim) = grad_term1[dim] - enrich_DN_DX_gp_dim ;
/*                                 K_enrich(0,0) -= permittivity_pos*pos_weights(pos_e)
                    *( grad_term1[dim] - enrich_DN_DX_gp_dim )*( grad_term1[dim] - enrich_DN_DX_gp_dim );
                pos_grad_enrich (pos_e,dim) = grad_term1[dim] - enrich_DN_DX_gp_dim ; */
            }    
        }
        
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //// NEGATIVE SIDE
        Vector neg_Nenr = ZeroVector(neg_weights.size());
        BoundedMatrix<double, 2, 2 > neg_grad_enrich;
        noalias(neg_grad_enrich) = ZeroMatrix(1,1);
        for (unsigned int neg_e = 0; neg_e < num_neg_gp; neg_e++){           
            lhs -= conductivity_neg*neg_weights(neg_e)*prod(neg_DN_DX[neg_e],trans(neg_DN_DX[neg_e]));
            //lhs -= permittivity_neg*neg_weights(neg_e)*prod(neg_DN_DX[neg_e],trans(neg_DN_DX[neg_e]));

            double term2 = 0;
            Vector grad_term1 = ZeroVector(num_dim);
            Vector grad_term2 = ZeroVector(num_dim);  
            double Nenr_1 = 0;
            double Nenr_term2 = 0;
            for (unsigned int i = 0; i < num_dof; i++){ 
                term2 += neg_N(neg_e,i)*distance(i);         
                for (unsigned int dim = 0; dim < num_dim; dim++){
                    grad_term1[dim] += neg_DN_DX[neg_e](i,dim)*std::abs(distance(i));
                    grad_term2[dim] += neg_DN_DX[neg_e](i,dim)* distance(i);
                }
            }  
           
        
            // components for codensation
            for (unsigned int dim = 0; dim < num_dim; ++dim){
                double sign = 1.0;
                if (term2 < 0.0) sign = -1.0;
                const double enrich_DN_DX_gp_dim = sign*grad_term2[dim];
                for (unsigned int i = 0; i < num_dof; ++i){
                    B_enrich(0, i) -= conductivity_neg*neg_weights(neg_e)
                        *neg_DN_DX[neg_e](i,dim)*( grad_term1[dim] - enrich_DN_DX_gp_dim );            
/*                                             B_enrich(0, i) -= permittivity_neg*neg_weights(neg_e)
                        *neg_DN_DX[neg_e](i,dim)*( grad_term1[dim] - enrich_DN_DX_gp_dim );          */       
                }
                K_enrich(0,0) -= conductivity_neg*neg_weights(neg_e)
                    *( grad_term1[dim] - enrich_DN_DX_gp_dim )*( grad_term1[dim] - enrich_DN_DX_gp_dim );
                neg_grad_enrich (neg_e,dim) = grad_term1[dim] - enrich_DN_DX_gp_dim ;

/*                                 K_enrich(0,0) -= permittivity_neg*neg_weights(neg_e)
                    *( grad_term1[dim] - enrich_DN_DX_gp_dim )*( grad_term1[dim] - enrich_DN_DX_gp_dim );
                neg_grad_enrich (neg_e,dim) = grad_term1[dim] - enrich_DN_DX_gp_dim ; */
            }
        }
    
        ///adding flux contribution

        Bij_flux_enrich = B_enrich - Bij_flux_enrich;
        Bji_flux_enrich = B_enrich ;//- Bji_flux_enrich;
        K_flux_enrich = K_enrich - K_flux_enrich;

        // saving enrich data 
        const double inv_K_enrich = K_flux_enrich(0,0);
        this-> SetValue(INV_K_ENRICH,(inv_K_enrich));

        Vector bij = ZeroVector(3);
        Vector bji = ZeroVector(3);

        for (unsigned int i = 0; i < num_nodes; ++i){   
            
            bij(i) = Bij_flux_enrich(0,i);//B_enrich(0,i);      
            bji(i) = Bji_flux_enrich(0,i);//B_enrich(0,i);
        }  

        this-> SetValue(BIJ_ENRICH_ROW,(bij));
        this-> SetValue(BJI_ENRICH_ROW,(bji));

         for (unsigned int i = 0; i < num_dim; ++i){
            (this->GetValue(POS_GRAD_ENRICH))(i) = pos_grad_enrich(0,i);
            (this->GetValue(NEG_GRAD_ENRICH))(i) = neg_grad_enrich(0,i);
        } 
        /// CONDENSATION 
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //noalias(rLeftHandSideMatrix) = lhs - 1/K_enrich(0,0)*prod(trans(B_enrich),B_enrich);  //Please check the condensation
        noalias(rLeftHandSideMatrix) = lhs - 1/K_flux_enrich(0,0)*prod(trans(Bij_flux_enrich),Bji_flux_enrich);
        noalias(rRightHandSideVector) = rhs - prod(rLeftHandSideMatrix,values);
        
        (this->GetValue(RHS_ENRICH))= rRightHandSideVector;
    } 

    KRATOS_CATCH("");
    
}

/**
 * this is called during the assembling process in order
 * to calculate all elemental contributions to the global system
 * matrix and the right hand side
 * @param rLeftHandSideMatrix: the elemental left hand side matrix
 * @param rRightHandSideVector: the elemental right hand side
 * @param rCurrentProcessInfo: the current process info instance
 */
template< >
void ElectrostaticElement<3>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    /// general sets  /////
    const int num_dim  = 3;
    const int num_nodes  = num_dim + 1;
    const int num_dof = num_nodes; 
    const int num_enrch = num_nodes+1; 
    GeometryData::ShapeFunctionsGradientsType DN_DX;
    Matrix N;
    Vector DetJ;
    Vector weights;
    
    //////////////////////////////////////////////////////////////////////////
    ///// matrix /////////////////////////////
    Vector distance(num_dof);
    Vector distances(num_dof);   
    Vector values(num_dof);
    BoundedMatrix<double,num_dof,num_dof> lhs = ZeroMatrix(num_dof,num_dof);
    Vector rhs = ZeroVector(num_dof);
    //array_1d<double,3> surface_sources;
    /////////////////////////////////////////////////////////////////////////////////
    
    const double conductivity_pos = GetProperties()[CONDUCTIVITYPOS];
    const double conductivity_neg = GetProperties()[CONDUCTIVITYNEG];
    const double permittivity_pos = GetProperties()[PERMITTIVITYPOS];
    const double permittivity_neg = GetProperties()[PERMITTIVITYNEG];
  
    ///////////////////////////////////////////////////////////////////////////////////
    const GeometryData::IntegrationMethod integration_method = GeometryData::IntegrationMethod::GI_GAUSS_1;
    GeometryType::Pointer p_geometry = this->pGetGeometry();
    const unsigned int number_of_gauss_points = p_geometry->IntegrationPointsNumber(integration_method);
    p_geometry->ShapeFunctionsIntegrationPointsGradients(DN_DX, DetJ, integration_method);
    ////////////////////////////////////////////////////////////////////////////////////////////////////
 
    if(rLeftHandSideMatrix.size1() != num_dof)
        rLeftHandSideMatrix.resize(num_dof,num_dof,false); //resizing the system in case it does not have the right size 

    if(rRightHandSideVector.size() != num_dof)
        rRightHandSideVector.resize(num_dof,false);

    for (unsigned int i_node=0; i_node < num_nodes; ++i_node){
        const double dist = (*p_geometry)[i_node].FastGetSolutionStepValue(DISTANCE);

        if (dist == 0.0){                   
            (*p_geometry)[i_node].FastGetSolutionStepValue(DISTANCE) = 1.0e-10;
        }
    }
    
    unsigned int nneg=0, npos=0;
    
    for (unsigned int i_node=0; i_node < num_nodes; ++i_node){
        const double dist = (*p_geometry)[i_node].FastGetSolutionStepValue(DISTANCE);
        distance[i_node] = dist;
        values[i_node] = (*p_geometry)[i_node].FastGetSolutionStepValue(EPOTENTIAL);
        
        if (dist > 0.0){
            npos += 1;
        }
        else{
            nneg += 1;
        }
    }
    
    if (nneg == num_nodes || npos == num_nodes){
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = p_geometry->IntegrationPoints(integration_method);
        if (weights.size() != number_of_gauss_points) {
            weights.resize(number_of_gauss_points,false);
        }

        //calculating the area of the elementes
        for (unsigned int gp = 0; gp < number_of_gauss_points; gp++){
            weights[gp] = DetJ[gp] * IntegrationPoints[gp].Weight();    
           
            if (npos == num_nodes){ //same results as enriched_n_element code
                lhs = -conductivity_pos*weights[gp]*prod(DN_DX[gp],trans(DN_DX[gp])); 
                //lhs = -permittivity_pos*weights[gp]*prod(DN_DX[gp],trans(DN_DX[gp])); 
            }else {
                lhs = -conductivity_neg*weights[gp]*prod(DN_DX[gp],trans(DN_DX[gp]));
                //lhs = -permittivity_neg*weights[gp]*prod(DN_DX[gp],trans(DN_DX[gp]));
            }
        }

        noalias(rLeftHandSideMatrix) = lhs;
        noalias(rRightHandSideVector) = rhs - prod(rLeftHandSideMatrix,values);  
    } else {  
        //////////////////////////////////////////////////////////////////////////////////////////
        
        Vector cut_edges = ZeroVector(4);
        // Implementing the boundary condition on distance gradient
        
        //const auto& neighbour_elems = this->GetValue(NEIGHBOUR_ELEMENTS);

         for (unsigned int i_ne = 0; i_ne < 4; i_ne++){
                cut_edges(i_ne)=i_ne;
        } 
        ///////////////////////// SPLITTED SHAPE FUNCTIONS  //////////////////////////////////////  
        ModifiedShapeFunctions::Pointer p_modified_sh_func =
            Kratos::make_shared<Tetrahedra3D4ModifiedShapeFunctions>(p_geometry, distance);

        Matrix neg_N, pos_N;
        GeometryType::ShapeFunctionsGradientsType neg_DN_DX, pos_DN_DX;
        Vector neg_weights, pos_weights;

        Matrix int_N;
        GeometryType::ShapeFunctionsGradientsType int_DN_DX;
        Vector int_weights; 

        ModifiedShapeFunctions::AreaNormalsContainerType int_normal, pos_extr_normal, neg_extr_normal;

        //ModifiedShapeFunctions::ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues(
        Matrix pos_extr_N, neg_extr_N ; 
        GeometryType::ShapeFunctionsGradientsType pos_extr_DN_DX, neg_extr_DN_DX;
        Vector pos_extr_weights, neg_extr_weights; 
        


        ///////////////////////// SPLITTED SHAPE FUNCTIONS  //////////////////////////////////////  
        p_modified_sh_func->ComputePositiveSideShapeFunctionsAndGradientsValues(
            pos_N,        // N
            pos_DN_DX,    // DN_DX
            pos_weights,  // weight * detJ
            integration_method);

        p_modified_sh_func->ComputeNegativeSideShapeFunctionsAndGradientsValues(
            neg_N,        // N
            neg_DN_DX,    // DN_DX
            neg_weights,  // weight * detJ
            integration_method);
            //KRATOS_WATCH(neg_weights)       
        /////////////////////////////////////////////////////////////////////////////////////////////////
        p_modified_sh_func->ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
            int_N,
            int_DN_DX,
            int_weights,
            integration_method);  

        p_modified_sh_func->ComputeNegativeSideInterfaceAreaNormals(
            int_normal,  // it is not normalized 
            integration_method);
        /////////////////////////////////////////////////////////////////////////////////////////////////
        //NORMAL CONTRIBUTION 
        BoundedMatrix<double, 1, 1 > K_enrich;
        noalias(K_enrich) = ZeroMatrix(1,1);
        BoundedMatrix<double, 1, 4 > B_enrich;
        noalias(B_enrich) = ZeroMatrix(1,4);
        //////////////////////////////
        // FLUX CONTRIBUTION 
        BoundedMatrix<double, 1, 1 > K_flux_enrich;
        noalias(K_flux_enrich) = ZeroMatrix(1,1);
        BoundedMatrix<double, 1, 4 > Bij_flux_enrich;
        noalias(Bij_flux_enrich) = ZeroMatrix(1,4);
        BoundedMatrix<double, 1, 4 > Bji_flux_enrich;
        noalias(Bji_flux_enrich) = ZeroMatrix(1,4);
        /////////////////////////////////////////////////////////////////////////////////////////////////
        for (unsigned int n = 0 ; n < cut_edges.size() ; n++ ){
            if (cut_edges[n] < 4)
            {
                unsigned int cut_edge = cut_edges[n];
                unsigned int pos_FaceId = cut_edge;
                unsigned int neg_FaceId = cut_edge; 
                unsigned int pos_normal_FaceId = cut_edge;
                unsigned int neg_normal_FaceId = cut_edge;
                ////////////////////////////////////////////////////////////////////////////////////////////////////////////
                p_modified_sh_func->ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues(
                    pos_extr_N,
                    pos_extr_DN_DX,
                    pos_extr_weights,
                    pos_FaceId,
                    integration_method);  

                p_modified_sh_func->ComputePositiveExteriorFaceAreaNormals(
                    pos_extr_normal,  // it is not normalized 
                    pos_normal_FaceId,
                    integration_method);
                
                for (unsigned int gp_epos = 0; gp_epos < pos_extr_weights.size() ; gp_epos++){
                    double term2 = 0;
                    Vector grad_term1 = ZeroVector(num_dim);
                    Vector grad_term2 = ZeroVector(num_dim);  
                    Vector n_dn_dx = ZeroVector (3);
                    Vector n_dn_dx_enrch = ZeroVector(1);
                    double Nenr_1 = 0;
                    double Nenr_term2 = 0;

                    //// normalized 
                    Vector normal_pos_side = ZeroVector(num_dim);
                    normal_pos_side = normal_pos_side + pos_extr_normal[gp_epos];
                    normal_pos_side = (1.0 / norm_2(normal_pos_side)) * normal_pos_side;
                                
                    for (unsigned int i = 0; i < num_dof; i++){   
                        term2 += pos_extr_N(gp_epos,i)*distance(i);          
                        for (unsigned int dim = 0; dim < num_dim; dim++){
                            grad_term1[dim] += pos_extr_DN_DX[gp_epos](i,dim)*std::abs(distance(i));
                            grad_term2[dim] += pos_extr_DN_DX[gp_epos](i,dim)*distance(i);
                        }    
                        Nenr_1 += pos_extr_N(gp_epos,i)*std::abs(distance(i)); // first term of the enrichment function 
                        Nenr_term2 += pos_extr_N(gp_epos,i)*(distance(i)); // second term of the enrichment function 
                    }  
        
                    double sign = 1.0;
                    if (term2 < 0.0) sign = -1.0;
                    for (unsigned int dim = 0; dim < num_dim; ++dim){
                        n_dn_dx_enrch[0] += normal_pos_side(dim)*(grad_term1[dim] - sign*grad_term2[dim]); // added n * grad dn_dx enriched
                    }

                    n_dn_dx = prod(normal_pos_side, trans(pos_extr_DN_DX[gp_epos])); // vector (n).(DN_DX)
                    const double pos_Nenr = Nenr_1 - std::abs(Nenr_term2);
                    
                    // components for codensation
                    for (unsigned int i = 0; i < num_dof; ++i){
                        Bij_flux_enrich(0,i) -= conductivity_pos*pos_extr_weights[gp_epos]* pos_Nenr * n_dn_dx(i); // adde
                        Bji_flux_enrich(0,i) -= conductivity_pos*pos_extr_weights[gp_epos]* pos_extr_N(gp_epos,i) * n_dn_dx_enrch[0]; // added  
                        //Bij_flux_enrich(0,i) -= permittivity_pos*pos_extr_weights[gp_epos]* pos_Nenr * n_dn_dx(i); // adde
                        //Bji_flux_enrich(0,i) -= permittivity_pos*pos_extr_weights[gp_epos]* pos_extr_N(gp_epos,i) * n_dn_dx_enrch[0]; // added  
                    }
                
                    K_flux_enrich(0,0) -= conductivity_pos*pos_extr_weights[gp_epos]*pos_Nenr*n_dn_dx_enrch[0];
                    //K_flux_enrich(0,0) -= permittivity_pos*pos_extr_weights[gp_epos]*pos_Nenr*n_dn_dx_enrch[0];
                
                }
                
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
                p_modified_sh_func->ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues(
                    neg_extr_N,
                    neg_extr_DN_DX,
                    neg_extr_weights,
                    neg_FaceId,
                    integration_method);   

                p_modified_sh_func->ComputeNegativeExteriorFaceAreaNormals(
                    neg_extr_normal,  // it is not normalized 
                    neg_normal_FaceId,
                    integration_method);
                ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //// NEGATIVE SIDE
                for (unsigned int gp_eneg = 0; gp_eneg < neg_extr_weights.size() ; gp_eneg++){
                    double term2 = 0;
                    Vector grad_term1 = ZeroVector(num_dim);
                    Vector grad_term2 = ZeroVector(num_dim);  
                    Vector n_dn_dx = ZeroVector (3);
                    Vector n_dn_dx_enrch = ZeroVector(1);   
                    double Nenr_1 = 0;
                    double Nenr_term2 = 0;

                    //// normalized 
                    Vector normal_neg_side = ZeroVector(num_dim);
                    normal_neg_side = normal_neg_side + neg_extr_normal[gp_eneg];
                    normal_neg_side = (1.0 / norm_2(normal_neg_side)) * normal_neg_side;
                
                    for (unsigned int i = 0; i < num_dof; i++){   
                        term2 += neg_extr_N(gp_eneg,i)*distance(i);          
                        for (unsigned int dim = 0; dim < num_dim; dim++){
                            grad_term1[dim] += neg_extr_DN_DX[gp_eneg](i,dim)*std::abs(distance(i));
                            grad_term2[dim] += neg_extr_DN_DX[gp_eneg](i,dim)*distance(i);
                        }    
                        Nenr_1 += neg_extr_N(gp_eneg,i)*std::abs(distance(i)); // first term of the enrichment function 
                        Nenr_term2 += neg_extr_N(gp_eneg,i)*(distance(i)); // second term of the enrichment function 
                    }  
                
                    double sign = 1.0;
                    if (term2 < 0.0) sign = -1.0;

                    for (unsigned int dim = 0; dim < num_dim; ++dim){
                        n_dn_dx_enrch[0] += normal_neg_side(dim)*(grad_term1[dim] - sign*grad_term2[dim]); // added n * grad dn_dx enriched
                    }

                    n_dn_dx = prod(normal_neg_side, trans(neg_extr_DN_DX[gp_eneg])); // vector (n).(DN_DX)
                    const double neg_Nenr = Nenr_1 - std::abs(Nenr_term2);

                    // components for codensation
                    for (unsigned int i = 0; i < num_dof; ++i){
                        Bij_flux_enrich(0,i) -= conductivity_neg*neg_extr_weights[gp_eneg]* neg_Nenr * n_dn_dx(i); // adde
                        Bji_flux_enrich(0,i) -= conductivity_neg*neg_extr_weights[gp_eneg]* neg_extr_N(gp_eneg,i) * n_dn_dx_enrch[0]; // added  
                        //Bij_flux_enrich(0,i) -= permittivity_neg*neg_extr_weights[gp_eneg]* neg_Nenr * n_dn_dx(i); // adde
                        //Bji_flux_enrich(0,i) -= permittivity_neg*neg_extr_weights[gp_eneg]* neg_extr_N(gp_eneg,i) * n_dn_dx_enrch[0]; // added  
                    }
            
                    K_flux_enrich(0,0) -= conductivity_neg*neg_extr_weights[gp_eneg]*neg_Nenr*n_dn_dx_enrch[0];
                    //K_flux_enrich(0,0) -= permittivity_neg*neg_extr_weights[gp_eneg]*neg_Nenr*n_dn_dx_enrch[0];
                }
                //KRATOS_WATCH("end flux_neg contribution")
            }
        }
        /////////////////////////////////////////////////////////////////////////////////////////////////       
        const int num_pos_gp = pos_weights.size();
        const int num_neg_gp = neg_weights.size();
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //// POSITIVE SIDE 
        Vector pos_Nenr = ZeroVector(pos_weights.size());
        Vector pos_grad_enrich = ZeroVector(3);

        for (unsigned int pos_e = 0; pos_e < num_pos_gp; pos_e++){   
            lhs -= conductivity_pos*pos_weights(pos_e)*prod(pos_DN_DX[pos_e],trans(pos_DN_DX[pos_e]));
        }
 
        for (unsigned int pos_e = 0; pos_e < num_pos_gp; pos_e++){
            double term2 = 0;
            Vector grad_term1 = ZeroVector(num_dim);
            Vector grad_term2 = ZeroVector(num_dim);  
            double Nenr_1 = 0;
            double Nenr_term2 = 0;
            for (unsigned int i = 0; i < num_dof; i++){   
                term2 += pos_N(pos_e,i)*distance(i);          
                for (unsigned int dim = 0; dim < num_dim; dim++){
                    grad_term1[dim] += pos_DN_DX[pos_e](i,dim)*std::abs(distance(i));
                    grad_term2[dim] += pos_DN_DX[pos_e](i,dim)*distance(i);
                }    
            }  
            // components for codensation
            for (unsigned int dim = 0; dim < num_dim; ++dim){
                double sign = 1.0;
                if (term2 < 0.0) sign = -1.0;
                const double enrich_DN_DX_gp_dim = sign*grad_term2[dim];
                for (unsigned int i = 0; i < num_dof; ++i){
                    B_enrich(0, i) -= conductivity_pos*pos_weights(pos_e)
                        *pos_DN_DX[pos_e](i,dim)*( grad_term1[dim] - enrich_DN_DX_gp_dim );
                }
                K_enrich(0,0) -= conductivity_pos*pos_weights(pos_e)
                    *( grad_term1[dim] - enrich_DN_DX_gp_dim )*( grad_term1[dim] - enrich_DN_DX_gp_dim );
                
                pos_grad_enrich (dim) = grad_term1[dim] - enrich_DN_DX_gp_dim ; 
            }   
        } 

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //// NEGATIVE SIDE
        Vector neg_Nenr = ZeroVector(neg_weights.size());
        Vector neg_grad_enrich = ZeroVector(3);
        
        for (unsigned int neg_e = 0; neg_e < num_neg_gp; neg_e++){           
            lhs -= conductivity_neg*neg_weights(neg_e)*prod(neg_DN_DX[neg_e],trans(neg_DN_DX[neg_e]));
        }

        for (unsigned int neg_e = 0; neg_e < num_neg_gp; neg_e++){    
            double term2 = 0;
            Vector grad_term1 = ZeroVector(num_dim);
            Vector grad_term2 = ZeroVector(num_dim);  
            double Nenr_1 = 0;
            double Nenr_term2 = 0;
            for (unsigned int i = 0; i < num_dof; i++){ 
                term2 += neg_N(neg_e,i)*distance(i);         
                for (unsigned int dim = 0; dim < num_dim; dim++){
                    grad_term1[dim] += neg_DN_DX[neg_e](i,dim)*std::abs(distance(i));
                    grad_term2[dim] += neg_DN_DX[neg_e](i,dim)* distance(i);
                }
            }  
        
            // components for codensation
            for (unsigned int dim = 0; dim < num_dim; ++dim){
                double sign = 1.0;
                if (term2 < 0.0) sign = -1.0;
                const double enrich_DN_DX_gp_dim = sign*grad_term2[dim];
                for (unsigned int i = 0; i < num_dof; ++i){
                    B_enrich(0, i) -= conductivity_neg*neg_weights(neg_e)
                        *neg_DN_DX[neg_e](i,dim)*( grad_term1[dim] - enrich_DN_DX_gp_dim );               
                }
                K_enrich(0,0) -= conductivity_neg*neg_weights(neg_e)
                    *( grad_term1[dim] - enrich_DN_DX_gp_dim )*( grad_term1[dim] - enrich_DN_DX_gp_dim ); 
                
                neg_grad_enrich (dim) = grad_term1[dim] - enrich_DN_DX_gp_dim ; 
            } 
        } 

        Bij_flux_enrich = B_enrich - Bij_flux_enrich;
        Bji_flux_enrich = B_enrich ;//- Bji_flux_enrich;
        K_flux_enrich = K_enrich - K_flux_enrich;

        const double inv_K_enrich = K_flux_enrich(0,0);
        this-> SetValue(INV_K_ENRICH,(inv_K_enrich));
        
        Vector bij = ZeroVector(4);
        Vector bji = ZeroVector(4);

        for (unsigned int i = 0; i < num_nodes; ++i){   
            
            bij(i) = Bij_flux_enrich(0,i);//B_enrich(0,i);      
            bji(i) = Bji_flux_enrich(0,i);//B_enrich(0,i);
        }  

        this-> SetValue(BIJ_ENRICH_ROW,(bij));
        this-> SetValue(BJI_ENRICH_ROW,(bji));

        for (unsigned int i = 0; i < num_dim; ++i){
            (this->GetValue(POS_GRAD_ENRICH))(i) = pos_grad_enrich(i);
            (this->GetValue(NEG_GRAD_ENRICH))(i) = neg_grad_enrich(i);
        }
                      
        /// CONDENSATION 
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //noalias(rLeftHandSideMatrix) = lhs - 1/K_enrich(0,0)*prod(trans(B_enrich),B_enrich);  //Please check the condensation
        noalias(rLeftHandSideMatrix) = lhs - 1/K_flux_enrich(0,0) * prod(trans(Bij_flux_enrich),Bji_flux_enrich);
        noalias(rRightHandSideVector) = rhs - prod(rLeftHandSideMatrix,values);
        
        (this->GetValue(RHS_ENRICH))= rRightHandSideVector;
    } 

    KRATOS_CATCH("");
}


template< unsigned int TDim >
void ElectrostaticElement<TDim>::GetValueOnIntegrationPoints(const Variable<array_1d<double,3>>& rVariable, 
std::vector<array_1d<double,3> >& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    if(rVariable==EFIELD)
    {
    CalculateOnIntegrationPoints (rVariable, rValues, rCurrentProcessInfo);    
    } 
    else if(rVariable==EFORCE)
    {
    CalculateOnIntegrationPoints (rVariable, rValues, rCurrentProcessInfo);    
    }  
    else if(rVariable==SCHARGE)
    {
    CalculateOnIntegrationPoints (rVariable, rValues, rCurrentProcessInfo);    
    }   
    else if(rVariable==EFIELDPOS)
    {
    CalculateOnIntegrationPoints (rVariable, rValues, rCurrentProcessInfo);    
    }  
    else if(rVariable==EFIELDNEG)
    {
    CalculateOnIntegrationPoints (rVariable, rValues, rCurrentProcessInfo);    
    } 
}


///////************************* ELECTRIC FIELD and FORCE CALCULATION ************//////////////////////////////
template< >
void ElectrostaticElement<2>::CalculateOnIntegrationPoints(const Variable<array_1d<double,3>>& rVariable, 
std::vector<array_1d<double,3> >& Output, const ProcessInfo& rCurrentProcessInfo)
{ 
    const int num_dim  = 2;
    const int num_nodes  = num_dim + 1;
    const int num_dof = num_nodes; 
    //KRATOS_WATCH("h1")
    GeometryData::ShapeFunctionsGradientsType DN_DX;
    Matrix N;
    Vector DetJ;
    Vector weights;
    Vector distance(num_dof);
    Vector values(num_dof);
    Vector efield = ZeroVector(num_dof);
    Vector efield_pos_side = ZeroVector(num_dof);
    Vector efield_neg_side = ZeroVector(num_dof);
    Vector n_Fe = ZeroVector(2);
    Vector n_Te_pos = ZeroVector(3);
    Vector n_Te_neg = ZeroVector(3);
    Vector n_Te = ZeroVector(3);
    Vector charge_density = ZeroVector(2);
    /////////////////////////////////////////////////////////////////////////////////
    const double permittivity_neg = this->GetProperties().GetValue(PERMITTIVITYNEG);
    const double permittivity_pos = this->GetProperties().GetValue(PERMITTIVITYPOS);
    const double conductivity_neg = this->GetProperties().GetValue(CONDUCTIVITYNEG);
    const double conductivity_pos = this->GetProperties().GetValue(CONDUCTIVITYPOS);
    ///////////////////////////////////////////////////////////////////////////////////
    const GeometryData::IntegrationMethod integration_method = GeometryData::IntegrationMethod::GI_GAUSS_1;
    GeometryType::Pointer p_geometry = this->pGetGeometry();
    const unsigned int number_of_gauss_points = p_geometry->IntegrationPointsNumber(integration_method);
    p_geometry->ShapeFunctionsIntegrationPointsGradients(DN_DX, DetJ, integration_method);

    for (unsigned int i_node=0; i_node < num_nodes; ++i_node){
        const double dist = (*p_geometry)[i_node].FastGetSolutionStepValue(DISTANCE);
        if (dist == 0.0){
            (*p_geometry)[i_node].FastGetSolutionStepValue(DISTANCE) = -1.0e-10;
        }
    }
    
    unsigned int nneg=0, npos=0;

    for (unsigned int i_node=0; i_node < num_nodes; ++i_node){
        const double dist = (*p_geometry)[i_node].FastGetSolutionStepValue(DISTANCE);
        distance[i_node] = dist;
        values[i_node] = (*p_geometry)[i_node].FastGetSolutionStepValue(EPOTENTIAL);
        
        if (dist > 0.0){
            npos += 1;
        }
        else{
            nneg += 1;
        }
    }

    if (nneg == num_nodes || npos == num_nodes){
        
        efield = -prod(trans(DN_DX[0]),values); 
    } else {   
        ///////////////////////// SPLITTED SHAPE FUNCTIONS  //////////////////////////////////////  
        ModifiedShapeFunctions::Pointer p_modified_sh_func =
            Kratos::make_shared<Triangle2D3ModifiedShapeFunctions>(p_geometry, distance);

        Matrix neg_N, pos_N;
        GeometryType::ShapeFunctionsGradientsType neg_DN_DX, pos_DN_DX;
        Vector neg_weights, pos_weights;

        Matrix int_N;
        GeometryType::ShapeFunctionsGradientsType int_DN_DX;
        Vector int_weights;

        ModifiedShapeFunctions::AreaNormalsContainerType int_normal;

        ///////////////////////// SPLITTED SHAPE FUNCTIONS  //////////////////////////////////////  
        p_modified_sh_func->ComputePositiveSideShapeFunctionsAndGradientsValues(
            pos_N,        // N
            pos_DN_DX,    // DN_DX
            pos_weights,  // weight * detJ
            integration_method);

        p_modified_sh_func->ComputeNegativeSideShapeFunctionsAndGradientsValues(
            neg_N,        // N
            neg_DN_DX,    // DN_DX
            neg_weights,  // weight * detJ
            integration_method);
            //KRATOS_WATCH(neg_weights)       
        /////////////////////////////////////////////////////////////////////////////////////////////////
        p_modified_sh_func->ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
                int_N,
                int_DN_DX,
                int_weights,
                integration_method);  

        p_modified_sh_func->ComputeNegativeSideInterfaceAreaNormals(
                int_normal,  // it is not normalized 
                integration_method);
                            
        //////////////////////////////////////////////////////////////////////////////////////////////// 
        //DECONDENSATION 
        double k_enrich = this->GetValue(INV_K_ENRICH); //reading information from the main loop 
        Vector bij_enrich_row = ZeroVector(num_nodes);
        bij_enrich_row = this->GetValue(BIJ_ENRICH_ROW); //reading information from the main loop 

        double phi_discontinuous = 0;
        for (unsigned int i_int = 0; i_int<3;i_int++){
            phi_discontinuous -= bij_enrich_row(i_int)*values(i_int); //B_enr*phi
        }
        phi_discontinuous /= k_enrich; //(B_enr*phi/K_enr) -> phi(decondensed)

        Vector N_enr_1 = ZeroVector(2); // enric shape fuctions at the interface

        for (unsigned int i = 0 ; i < 3 ; i++){

            N_enr_1 (0) += int_N(i)*std::abs(distance(i)); // term 1 = N*abs(distance)
            N_enr_1 (1) += int_N(i)*distance(i); //term 2 = N*(distance)
           
        }
        double phi_interface = int_N(0)*values(0) + int_N(1)*values(1) + int_N(2)*values(2); //phi at the interface
        double N_enr = N_enr_1(0) - std::abs(N_enr_1(1)); // N*abs(distance) + abs(N*(distance))
        double phi_int_enrich = phi_interface + N_enr*phi_discontinuous; // enriched phi at the interface 
       
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // ELECTRIC FIELD     
            
        efield_pos_side = - prod(trans(pos_DN_DX[0]),values); // change the 
        efield_neg_side = - prod(trans(neg_DN_DX[0]),values); 

        

          for (unsigned int ef = 0 ; ef < num_nodes ; ef++){
            efield_pos_side(ef) -= this->GetValue(POS_GRAD_ENRICH)(ef) * phi_discontinuous; 
            efield_neg_side(ef) -= this->GetValue(NEG_GRAD_ENRICH)(ef) * phi_discontinuous;
        } 


        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //////// ELECTRIC FORCE (Fe=)


        Vector normal_efield = ZeroVector(2);
        normal_efield = normal_efield + int_normal[0];
        normal_efield = (1.0 / norm_2(normal_efield)) * normal_efield;


        n_Te_pos = permittivity_pos * (inner_prod (normal_efield,efield_pos_side) * efield_pos_side - 0.5 * inner_prod(efield_pos_side,efield_pos_side)*normal_efield);
        n_Te_neg = permittivity_neg * (inner_prod (normal_efield,efield_neg_side) * efield_neg_side - 0.5 * inner_prod(efield_neg_side,efield_neg_side)*normal_efield);
        
        n_Te = - n_Te_neg + n_Te_pos;

        charge_density[0] = - permittivity_neg * inner_prod (normal_efield,efield_neg_side) + permittivity_pos * inner_prod (normal_efield,efield_pos_side);      

    }

    this->SetValue(EFORCE,n_Te);

    const int total_gp = p_geometry->IntegrationPointsNumber(integration_method);
    
    for (unsigned int i = 0; i < total_gp; i++ ) {
    
        if(rVariable == SCHARGE)//electric force form the negative side 
        {
                noalias(Output[i]) = charge_density;
        }
        else if(rVariable == EFORCE)
        {
                noalias(Output[i]) = n_Te;
                Output[i][2]=0.0;
        }
        else if(rVariable == EFIELDPOS)
        {
                noalias(Output[i]) = efield_pos_side;
                Output[i][2]=0.0;
        }
        else if(rVariable == EFIELDNEG)
        {
                noalias(Output[i]) = efield_neg_side;
                Output[i][2]=0.0;
        }
        else
        if(rVariable == EFIELD)//electric force form the negative side 
        {
                noalias(Output[i]) = efield + efield_pos_side;
                Output[i][2]=0.0;
        }                                                                                                                                                                                                                                                   
        
    }   
}


///////************************* ELECTRIC FIELD and FORCE CALCULATION ************//////////////////////////////
template< >
void ElectrostaticElement<3>::CalculateOnIntegrationPoints(const Variable<array_1d<double,3>>& rVariable, 
std::vector<array_1d<double,3> >& Output, const ProcessInfo& rCurrentProcessInfo)
{ 
const int num_dim  = 3;
    const int num_nodes  = num_dim + 1;
    const int num_dof = num_nodes; 


    GeometryData::ShapeFunctionsGradientsType DN_DX;
    Matrix N;
    Vector DetJ;
    Vector weights;

    Vector distance(num_dof);
    Vector values(num_dof);
    Vector efield = ZeroVector(3);
    Vector efield_pos_side = ZeroVector(3);
    Vector efield_neg_side = ZeroVector(3);
    Vector n_Fe = ZeroVector(3);
    Vector n_Te_pos = ZeroVector(3);
    Vector n_Te_neg = ZeroVector(3);
    Vector n_Te = ZeroVector(3);
    Vector charge_density = ZeroVector(3);
    /////////////////////////////////////////////////////////////////////////////////
    const double conductivity_pos = GetProperties()[CONDUCTIVITYPOS];
    const double conductivity_neg = GetProperties()[CONDUCTIVITYNEG];
    const double permittivity_pos = GetProperties()[PERMITTIVITYPOS];
    const double permittivity_neg = GetProperties()[PERMITTIVITYNEG];
    ///////////////////////////////////////////////////////////////////////////////////
    const GeometryData::IntegrationMethod integration_method = GeometryData::IntegrationMethod::GI_GAUSS_1;
    GeometryType::Pointer p_geometry = this->pGetGeometry();
    const unsigned int number_of_gauss_points = p_geometry->IntegrationPointsNumber(integration_method);
    p_geometry->ShapeFunctionsIntegrationPointsGradients(DN_DX, DetJ, integration_method);

    for (unsigned int i_node=0; i_node < num_nodes; ++i_node){
        const double dist = (*p_geometry)[i_node].FastGetSolutionStepValue(DISTANCE);
        if (dist == 0.0){
            (*p_geometry)[i_node].FastGetSolutionStepValue(DISTANCE) = -1.0e-10;
        }
    }
    
    unsigned int nneg=0, npos=0;

    for (unsigned int i_node=0; i_node < num_nodes; ++i_node){
        const double dist = (*p_geometry)[i_node].FastGetSolutionStepValue(DISTANCE);
        distance[i_node] = dist;
        values[i_node] = (*p_geometry)[i_node].FastGetSolutionStepValue(EPOTENTIAL);
        
        if (dist > 0.0){
            npos += 1;
        }
        else{
            nneg += 1;
        }
    }

    if (nneg == num_nodes || npos == num_nodes){
        
        efield = -prod(trans(DN_DX[0]),values); 

    } else {   
        ///////////////////////// SPLITTED SHAPE FUNCTIONS  //////////////////////////////////////  
        ModifiedShapeFunctions::Pointer p_modified_sh_func =
            Kratos::make_shared<Tetrahedra3D4ModifiedShapeFunctions>(p_geometry, distance);

        Matrix neg_N, pos_N;
        GeometryType::ShapeFunctionsGradientsType neg_DN_DX, pos_DN_DX;
        Vector neg_weights, pos_weights;

        Matrix int_N;
        GeometryType::ShapeFunctionsGradientsType int_DN_DX;
        Vector int_weights;

        ModifiedShapeFunctions::AreaNormalsContainerType int_normal;

        ///////////////////////// SPLITTED SHAPE FUNCTIONS  //////////////////////////////////////  
        p_modified_sh_func->ComputePositiveSideShapeFunctionsAndGradientsValues(
            pos_N,        // N
            pos_DN_DX,    // DN_DX
            pos_weights,  // weight * detJ
            integration_method);

        p_modified_sh_func->ComputeNegativeSideShapeFunctionsAndGradientsValues(
            neg_N,        // N
            neg_DN_DX,    // DN_DX
            neg_weights,  // weight * detJ
            integration_method);
            //KRATOS_WATCH(neg_weights)       
        /////////////////////////////////////////////////////////////////////////////////////////////////
        p_modified_sh_func->ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
                int_N,
                int_DN_DX,
                int_weights,
                integration_method);  

        p_modified_sh_func->ComputeNegativeSideInterfaceAreaNormals(
                int_normal,  // it is not normalized 
                integration_method);
                            
        //////////////////////////////////////////////////////////////////////////////////////////////// 
        //DECONDENSATION 
        double k_enrich = this->GetValue(INV_K_ENRICH); //reading information from the main loop 
        Vector bij_enrich_row = ZeroVector(num_nodes);
        bij_enrich_row = this->GetValue(BIJ_ENRICH_ROW); //reading information from the main loop 

        double phi_discontinuous = 0;
        for (unsigned int i_int = 0; i_int < num_nodes ;i_int++){
            phi_discontinuous -= bij_enrich_row(i_int)*values(i_int); //B_enr*phi
        }
        phi_discontinuous /= k_enrich; //(B_enr*phi/K_enr) -> phi(decondensed)

   
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // ELECTRIC FIELD     

        efield_pos_side = - prod(trans(pos_DN_DX[0]),values);
        efield_neg_side = - prod(trans(neg_DN_DX[0]),values); //considering linear elements (one GP)
        
        for (unsigned int ef = 0 ; ef < 3 ; ef++){
            efield_pos_side(ef) -= this->GetValue(POS_GRAD_ENRICH)(ef) * phi_discontinuous; 
            efield_neg_side(ef) -= this->GetValue(NEG_GRAD_ENRICH)(ef) * phi_discontinuous;
        }   

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //////// ELECTRIC FORCE (Fe=)


        Vector normal_efield = ZeroVector(3);
        
        for (unsigned int n_size = 0; n_size < int_normal.size(); n_size++){
            normal_efield = normal_efield + int_normal[n_size];
        }
        //normal_efield = normal_efield + int_normal[0]; //check how many gauss points 
        normal_efield = (1.0 / norm_2(normal_efield)) * normal_efield;


        n_Te_pos = permittivity_pos * (inner_prod (normal_efield,efield_pos_side) * efield_pos_side - 0.5 * inner_prod(efield_pos_side,efield_pos_side)*normal_efield);

        n_Te_neg = permittivity_neg * (inner_prod (normal_efield,efield_neg_side) * efield_neg_side - 0.5 * inner_prod(efield_neg_side,efield_neg_side)*normal_efield);
        
        n_Te = - n_Te_neg + n_Te_pos;

        //charge_density[0] = - permittivity_neg * inner_prod (normal_efield,efield_neg_side) + permittivity_pos * inner_prod (normal_efield,efield_pos_side);
       

    }

    this->SetValue(EFORCE,n_Te);

    const int total_gp = p_geometry->IntegrationPointsNumber(integration_method);
    
    for (unsigned int i = 0; i < total_gp; i++ ) {
    
        if(rVariable == SCHARGE)//electric force form the negative side 
        {
                noalias(Output[i]) = charge_density;
        }
        else if(rVariable == EFORCE)
        {
                noalias(Output[i]) = n_Te;
                Output[i][2]=0.0;
        }
        else if(rVariable == EFIELDPOS)
        {
                noalias(Output[i]) = efield_pos_side;
                Output[i][2]=0.0;
        }
        else if(rVariable == EFIELDNEG)
        {
                noalias(Output[i]) = efield_neg_side;
                Output[i][2]=0.0;
        }
        else
        if(rVariable == EFIELD)//electric force form the negative side 
        {
                noalias(Output[i]) = efield + efield_pos_side;
                Output[i][2]=0.0;
        }                                                                                                                                                                                                                                                   
        
    }  
}


template< unsigned int TDim >
void ElectrostaticElement<TDim>::Calculate(const Variable<array_1d<double, 3 > >& rVariable,
array_1d<double, 3 > & Output, const ProcessInfo& rCurrentProcessInfo)
{
    
}


/**
 * this is called during the assembling process in order
 * to calculate the elemental left hand side matrix only
 * @param rLeftHandSideMatrix: the elemental left hand side matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
template< unsigned int TDim >
void ElectrostaticElement<TDim>::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental right hand side vector only
 * @param rRightHandSideVector: the elemental right hand side vector
 * @param rCurrentProcessInfo: the current process info instance
 */
template< unsigned int TDim >
void ElectrostaticElement<TDim>::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
}

/**
 * this is called during the assembling process in order
 * to calculate the first derivatives contributions for the LHS and RHS
 * @param rLeftHandSideMatrix: the elemental left hand side matrix
 * @param rRightHandSideVector: the elemental right hand side
 * @param rCurrentProcessInfo: the current process info instance
 */
template< unsigned int TDim >
void ElectrostaticElement<TDim>::CalculateFirstDerivativesContributions(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != 0)
        rLeftHandSideMatrix.resize(0, 0, false);
    if (rRightHandSideVector.size() != 0)
        rRightHandSideVector.resize(0, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental left hand side matrix for the first derivatives constributions
 * @param rLeftHandSideMatrix: the elemental left hand side matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
template< unsigned int TDim >
void ElectrostaticElement<TDim>::CalculateFirstDerivativesLHS(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != 0)
        rLeftHandSideMatrix.resize(0, 0, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental right hand side vector for the first derivatives constributions
 * @param rRightHandSideVector: the elemental right hand side vector
 * @param rCurrentProcessInfo: the current process info instance
 */
template< unsigned int TDim >
void ElectrostaticElement<TDim>::CalculateFirstDerivativesRHS(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    if (rRightHandSideVector.size() != 0)
        rRightHandSideVector.resize(0, false);
}

/**
 * ELEMENTS inherited from this class must implement this methods
 * if they need to add dynamic element contributions
 * note: second derivatives means the accelerations if the displacements are the dof of the analysis
 * note: time integration parameters must be set in the rCurrentProcessInfo before calling these methods
 * CalculateSecondDerivativesContributions,
 * CalculateSecondDerivativesLHS, CalculateSecondDerivativesRHS methods are : OPTIONAL
 */


/**
 * this is called during the assembling process in order
 * to calculate the second derivative contributions for the LHS and RHS
 * @param rLeftHandSideMatrix: the elemental left hand side matrix
 * @param rRightHandSideVector: the elemental right hand side
 * @param rCurrentProcessInfo: the current process info instance
 */
template< unsigned int TDim >
void ElectrostaticElement<TDim>::CalculateSecondDerivativesContributions(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != 0)
        rLeftHandSideMatrix.resize(0, 0, false);
    if (rRightHandSideVector.size() != 0)
        rRightHandSideVector.resize(0, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental left hand side matrix for the second derivatives constributions
 * @param rLeftHandSideMatrix: the elemental left hand side matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
template< unsigned int TDim >
void ElectrostaticElement<TDim>::CalculateSecondDerivativesLHS(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != 0)
        rLeftHandSideMatrix.resize(0, 0, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental right hand side vector for the second derivatives constributions
 * @param rRightHandSideVector: the elemental right hand side vector
 * @param rCurrentProcessInfo: the current process info instance
 */
template< unsigned int TDim >
void ElectrostaticElement<TDim>::CalculateSecondDerivativesRHS(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rRightHandSideVector.size() != 0)
        rRightHandSideVector.resize(0, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental mass matrix
 * @param rMassMatrix: the elemental mass matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
template< unsigned int TDim >
void ElectrostaticElement<TDim>::CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    if (rMassMatrix.size1() != 0)
        rMassMatrix.resize(0, 0, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental damping matrix
 * @param rDampingMatrix: the elemental damping matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
template< unsigned int TDim >
void ElectrostaticElement<TDim>::CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    if (rDampingMatrix.size1() != 0)
        rDampingMatrix.resize(0, 0, false);
}

/**
 * This method provides the place to perform checks on the completeness of the input
 * and the compatibility with the problem options as well as the contitutive laws selected
 * It is designed to be called only once (or anyway, not often) typically at the beginning
 * of the calculations, so to verify that nothing is missing from the input
 * or that no common error is found.
 * @param rCurrentProcessInfo
 * this method is: MANDATORY
 */
template< unsigned int TDim >
int ElectrostaticElement<TDim>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

       // Base class checks for positive Jacobian and Id > 0
      int ierr = Element::Check(rCurrentProcessInfo);
      if(ierr != 0) return ierr;
  
      // Check that all required variables have been registered
      KRATOS_CHECK_VARIABLE_KEY(EPOTENTIAL)
 
      unsigned const int number_of_points = GetGeometry().size();  //added cornejo
      // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
      for ( unsigned int i = 0; i < number_of_points; i++ )
      {
          auto &rnode = this->GetGeometry()[i];
          KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(EPOTENTIAL,rnode)
          KRATOS_CHECK_DOF_IN_NODE(EPOTENTIAL,rnode)
      }
  
      return ierr;
  
      KRATOS_CATCH("");
}

///@}
///@name Access
///@{


///@}
///@name Inquiry
///@{


///@}
///@name Input and output
///@{

/// Turn back information as a string.
template< unsigned int TDim >
std::string ElectrostaticElement<TDim>::Info() const {
    std::stringstream buffer;
    buffer << "ElectricEnrichedElement #" << Id();
    return buffer.str();
}

/// Print information about this object.
template< unsigned int TDim >
void ElectrostaticElement<TDim>::PrintInfo(std::ostream& rOStream) const {
    rOStream << "ElectricEnrichedElement #" << Id();
}

/// Print object's data.
template< unsigned int TDim >
void ElectrostaticElement<TDim>::PrintData(std::ostream& rOStream) const {
    pGetGeometry()->PrintData(rOStream);
}

///@}
///@name Friends
///@{

///@}

///@name Protected static Member Variables
///@{

///@}
///@name Protected member Variables
///@{

///@}
///@name Protected Operators
///@{

///@}
///@name Protected Operations
///@{

///@}
///@name Protected  Access
///@{

///@}
///@name Protected Inquiry
///@{

///@}
///@name Protected LifeCycle
///@{

///@}

///@name Static Member Variables
///@{

///@}
///@name Member Variables
///@{

///@}
///@name Private Operators
///@{

///@}
///@name Private Operations
///@{

///@}
///@name Serialization
///@{
template< unsigned int TDim >
void ElectrostaticElement<TDim>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );

    // List
    // To be completed with the class member list
}
template< unsigned int TDim >
void ElectrostaticElement<TDim>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element );

    // List
    // To be completed with the class member list
}

///@}
///@name Private  Access
///@{

///@}
///@name Private Inquiry
///@{

///@}
///@name Un accessible methods
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template< unsigned int TDim >
inline std::istream & operator >> (std::istream& rIStream, ElectrostaticElement<TDim>& rThis);

/// output stream function
template< unsigned int TDim >
inline std::ostream & operator << (std::ostream& rOStream, const ElectrostaticElement<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

template class ElectrostaticElement<2>;
template class ElectrostaticElement<3>;
} // namespace Kratos.
