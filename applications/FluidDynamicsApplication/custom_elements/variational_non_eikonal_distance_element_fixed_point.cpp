//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    KratosAppGenerator
//

// System includes
// see the header file

// External includes
// see the header file

// Include Base headers
// see the header file

// Project includes
#include "custom_elements/variational_non_eikonal_distance_element.h"
#include "custom_utilities/element_size_calculator.h"

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
VariationalNonEikonalDistanceElement::VariationalNonEikonalDistanceElement(IndexType NewId)
    : Element(NewId) 
{
}

/**
 * Constructor using an array of nodes
 */
VariationalNonEikonalDistanceElement::VariationalNonEikonalDistanceElement(IndexType NewId, const NodesArrayType& ThisNodes)
    : Element(NewId, ThisNodes) 
{
}

/**
 * Constructor using Geometry
 */
VariationalNonEikonalDistanceElement::VariationalNonEikonalDistanceElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry) 
{
}

/**
 * Constructor using Properties
 */
VariationalNonEikonalDistanceElement::VariationalNonEikonalDistanceElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties) 
{
}

/**
 * Copy Constructor
 */
VariationalNonEikonalDistanceElement::VariationalNonEikonalDistanceElement(VariationalNonEikonalDistanceElement const& rOther)
    : Element(rOther) 
{
}

/**
 * Destructor
 */
VariationalNonEikonalDistanceElement::~VariationalNonEikonalDistanceElement()
{
}

///@}
///@name Operators
///@{

/// Assignment operator.
VariationalNonEikonalDistanceElement & VariationalNonEikonalDistanceElement::operator=(VariationalNonEikonalDistanceElement const& rOther)
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
Element::Pointer VariationalNonEikonalDistanceElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<VariationalNonEikonalDistanceElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

/**
 * creates a new element pointer
 * @param NewId: the ID of the new element
 * @param pGeom: the geometry to be employed
 * @param pProperties: the properties assigned to the new element
 * @return a Pointer to the new element
 */
Element::Pointer VariationalNonEikonalDistanceElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<VariationalNonEikonalDistanceElement>(NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

/**
 * creates a new element pointer and clones the previous element data
 * @param NewId: the ID of the new element
 * @param ThisNodes: the nodes of the new element
 * @param pProperties: the properties assigned to the new element
 * @return a Pointer to the new element
 */
Element::Pointer VariationalNonEikonalDistanceElement::Clone(IndexType NewId, NodesArrayType const& ThisNodes) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<VariationalNonEikonalDistanceElement>(NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    KRATOS_CATCH("");
}

/**
 * this determines the elemental equation ID vector for all elemental
 * DOFs
 * @param rResult: the elemental equation ID vector
 * @param rCurrentProcessInfo: the current process info instance
 */
void VariationalNonEikonalDistanceElement::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    //KRATOS_INFO("VariationalNonEikonalDistanceElement") << "EquationIdVector" << std::endl;

    const int num_dim  = 3;
    const int num_nodes  = num_dim + 1;

    const int num_dof = num_nodes; //(num_dim + 1)*num_nodes;
    if (rResult.size() != num_dof){
        rResult.resize(num_dof, false);
    }

    for(unsigned int i=0; i<num_nodes; i++){
        rResult[i]  =  this->GetGeometry()[i].GetDof(DISTANCE_AUX2).EquationId();
        /* rResult[i*num_dim + num_nodes + 0]  =  this->GetGeometry()[i].GetDof(DISTANCE_GRADIENT_X).EquationId();
        rResult[i*num_dim + num_nodes + 1]  =  this->GetGeometry()[i].GetDof(DISTANCE_GRADIENT_Y).EquationId();
        rResult[i*num_dim + num_nodes + 2]  =  this->GetGeometry()[i].GetDof(DISTANCE_GRADIENT_Z).EquationId(); */

        /* rResult[i*(num_dim+1)]  =  this->GetGeometry()[i].GetDof(DISTANCE_AUX2).EquationId();
        rResult[i*(num_dim+1) + 1]  =  this->GetGeometry()[i].GetDof(DISTANCE_GRADIENT_X).EquationId();
        rResult[i*(num_dim+1) + 2]  =  this->GetGeometry()[i].GetDof(DISTANCE_GRADIENT_Y).EquationId();
        rResult[i*(num_dim+1) + 3]  =  this->GetGeometry()[i].GetDof(DISTANCE_GRADIENT_Z).EquationId(); */
    }
}

/**
 * determines the elemental list of DOFs
 * @param ElementalDofList: the list of DOFs
 * @param rCurrentProcessInfo: the current process info instance
 */
void VariationalNonEikonalDistanceElement::GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& CurrentProcessInfo)
{
    //KRATOS_INFO("VariationalNonEikonalDistanceElement") << "GetDofList" << std::endl;

    const int num_dim  = 3;
    const int num_nodes  = num_dim + 1;

    const int num_dof = num_nodes; //(num_dim + 1)*num_nodes;
    if (rElementalDofList.size() != num_dof){
        rElementalDofList.resize(num_dof);
    }

    for(unsigned int i=0; i<num_nodes; i++){
        rElementalDofList[i] = this->GetGeometry()[i].pGetDof(DISTANCE_AUX2);
        /* rElementalDofList[i*num_dim + num_nodes + 0] = this->GetGeometry()[i].pGetDof(DISTANCE_GRADIENT_X);
        rElementalDofList[i*num_dim + num_nodes + 1] = this->GetGeometry()[i].pGetDof(DISTANCE_GRADIENT_Y);
        rElementalDofList[i*num_dim + num_nodes + 2] = this->GetGeometry()[i].pGetDof(DISTANCE_GRADIENT_Z); */

        /* rElementalDofList[i*(num_dim+1)] = this->GetGeometry()[i].pGetDof(DISTANCE_AUX2);
        rElementalDofList[i*(num_dim+1) + 1] = this->GetGeometry()[i].pGetDof(DISTANCE_GRADIENT_X);
        rElementalDofList[i*(num_dim+1) + 2] = this->GetGeometry()[i].pGetDof(DISTANCE_GRADIENT_Y);
        rElementalDofList[i*(num_dim+1) + 3] = this->GetGeometry()[i].pGetDof(DISTANCE_GRADIENT_Z); */
    }
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
void VariationalNonEikonalDistanceElement::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //KRATOS_INFO("VariationalNonEikonalDistanceElement") << "Here 1" << std::endl;

    const int num_dim  = 3;
    const int num_nodes  = num_dim + 1;
    const int num_dof = num_nodes; //(num_dim + 1)*num_nodes;

    const double element_size = ElementSizeCalculator<3,4>::AverageElementSize(this->GetGeometry());

    //const double tau = 0.5;
    //const double penalty_curvature = 0.0;//1.0e-3; // Not possible for curvature itself since normalized DISTANCE_GRADIENT is needed.
    const double penalty_phi0 = 1.0e9;//0.0;//
    const double source_coeff = 3.0/(8.0*element_size);//1.0e0; //Usually we have ~10 elements across a Radius
    //const double dissipative_coefficient = 1.0e-8;

    GeometryData::ShapeFunctionsGradientsType DN_DX;
    Matrix N;
    Vector DetJ;
    Vector weights;

    //Vector curvatures(num_nodes);
    Vector distances0(num_nodes);
    //Matrix distance_gradient0(num_nodes, num_dim);
    Vector values(num_dof);
    VectorType grad_phi_avg;    //Recovered (lumped-mass) gradient
    VectorType grad_phi_old;    //Based on the unmodified DISTANCE
    VectorType grad_phi;        //Standard gradient

    //BoundedMatrix<double,num_dof,num_dof> tempPhi; // Only the 4 rows associated with Phi will be filled
    //tempPhi = ZeroMatrix(num_dof,num_dof);
    //BoundedMatrix<double,num_dof,num_dof> tempGradPhi; // Only 3*4 rows associated with GradPhi will be filled
    //tempGradPhi = ZeroMatrix(num_dof,num_dof);

    BoundedMatrix<double,num_dof,num_dof> lhs = ZeroMatrix(num_dof,num_dof);
    Vector rhs = ZeroVector(num_dof);

    const GeometryData::IntegrationMethod integration_method = GeometryData::GI_GAUSS_2;
    //const GeometryType& r_geometry = this->GetGeometry();
    GeometryType::Pointer p_geometry = this->pGetGeometry();
    const unsigned int number_of_gauss_points = p_geometry->IntegrationPointsNumber(integration_method);

    // Getting data for the given geometry
    p_geometry->ShapeFunctionsIntegrationPointsGradients(DN_DX, DetJ, integration_method);

    if (N.size1() != number_of_gauss_points || N.size2() != num_nodes) {
        N.resize(number_of_gauss_points,num_nodes,false);
    }
    N = p_geometry->ShapeFunctionsValues(integration_method);

    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = p_geometry->IntegrationPoints(integration_method);
    if (weights.size() != number_of_gauss_points) {
        weights.resize(number_of_gauss_points,false);
    }
    for (unsigned int gp = 0; gp < number_of_gauss_points; gp++){
        weights[gp] = DetJ[gp] * IntegrationPoints[gp].Weight();
    }

    unsigned int nneg=0, npos=0;

    for (unsigned int i_node=0; i_node < num_nodes; ++i_node){
        //curvatures(i_node) = (*p_geometry)[i_node].FastGetSolutionStepValue(CURVATURE);
        //KRATOS_INFO("VariationalNonEikonalDistanceElement, storing curvature") << curvatures(i_node) << std::endl;

        const double dist0 = (*p_geometry)[i_node].FastGetSolutionStepValue(DISTANCE);
        distances0(i_node) = dist0;

        if (dist0 > 0.0) npos += 1;
        else nneg += 1;

        const double dist_dof = (*p_geometry)[i_node].FastGetSolutionStepValue(DISTANCE_AUX2);
        values(i_node) = dist_dof;
        /* values(i_node*num_dim + num_nodes + 0) = (*p_geometry)[i_node].FastGetSolutionStepValue(DISTANCE_GRADIENT_X);
        values(i_node*num_dim + num_nodes + 1) = (*p_geometry)[i_node].FastGetSolutionStepValue(DISTANCE_GRADIENT_Y);
        values(i_node*num_dim + num_nodes + 2) = (*p_geometry)[i_node].FastGetSolutionStepValue(DISTANCE_GRADIENT_Z); */

        // values(i_node*(num_dim + 1) + 0) = (*p_geometry)[i_node].FastGetSolutionStepValue(DISTANCE_AUX2);
        // values(i_node*(num_dim + 1) + 1) = (*p_geometry)[i_node].FastGetSolutionStepValue(DISTANCE_GRADIENT_X);
        // values(i_node*(num_dim + 1) + 2) = (*p_geometry)[i_node].FastGetSolutionStepValue(DISTANCE_GRADIENT_Y);
        // values(i_node*(num_dim + 1) + 3) = (*p_geometry)[i_node].FastGetSolutionStepValue(DISTANCE_GRADIENT_Z);

        // distance_gradient0(i_node, 0) = (*p_geometry)[i_node].GetValue(DISTANCE_GRADIENT_X);
        // distance_gradient0(i_node, 1) = (*p_geometry)[i_node].GetValue(DISTANCE_GRADIENT_Y);
        // distance_gradient0(i_node, 2) = (*p_geometry)[i_node].GetValue(DISTANCE_GRADIENT_Z);
    }

    // num_dof = 4*num_nodes
    if(rLeftHandSideMatrix.size1() != num_dof)
        rLeftHandSideMatrix.resize(num_dof,num_dof,false); //resizing the system in case it does not have the right size 

    if(rRightHandSideVector.size() != num_dof)
        rRightHandSideVector.resize(num_dof,false);

    //KRATOS_INFO("VariationalNonEikonalDistanceElement") << "Here 2" << std::endl;

    const unsigned int step = rCurrentProcessInfo[FRACTIONAL_STEP];

    double diffusion = 0.0;

    for (unsigned int gp = 0; gp < number_of_gauss_points; gp++){
        grad_phi_avg = ZeroVector(num_dim);
        for (unsigned int i_node = 0; i_node < num_nodes; i_node++){
            grad_phi_avg += GetGeometry()[i_node].GetValue(DISTANCE_GRADIENT)*N(gp, i_node);
        }
        grad_phi = prod(trans(DN_DX[gp]),values);
        grad_phi_old = prod(trans(DN_DX[gp]),distances0);

        const double norm_grad_phi_avg = norm_2(grad_phi_avg);

        //KRATOS_WATCH(norm_grad_phi_avg)
        if (norm_grad_phi_avg > 1.0){
            diffusion = 1.0/norm_grad_phi_avg;
        } else{
            diffusion = (3.0 - 2.0*norm_grad_phi_avg)*norm_grad_phi_avg;
        }

        for (unsigned int i_node = 0; i_node < num_nodes; i_node++){
            for (unsigned int j_node = 0; j_node < num_nodes; j_node++){
                //tempPhi: LHS associated with the dissipative smoother
                //tempPhi(i_node*(num_dim + 1) + 0, j_node*(num_dim + 1) + 0) +=
                //    weights(gp) /* * (1/dissipative_coefficient) */ * N(gp, i_node) * N(gp, j_node);

                //rhs: RHS associated with the dissipative smoother
                //rhs(i_node*(num_dim + 1) + 0) +=
                //    weights(gp) /* * (1/dissipative_coefficient) */ * distances0(j_node) * N(gp, j_node) * N(gp, i_node);

                for (unsigned int k_dim = 0; k_dim < num_dim; k_dim++){

                    //tempPhi: LHS for the first 4 rows
                    /* tempPhi(i_node, j_node) += weights(gp) * tau * (DN_DX[gp])(i_node, k_dim) * (DN_DX[gp])(j_node, k_dim);
                    tempPhi(i_node, j_node*num_dim + num_nodes + k_dim) +=
                        weights(gp) * (1-tau) * N(gp, i_node) * (DN_DX[gp])(j_node, k_dim);

                    //tempGradPhi: LHS for the next 3*4 rows
                    tempGradPhi(i_node*num_dim + num_nodes + k_dim, j_node*num_dim + num_nodes + k_dim) +=
                        weights(gp) * N(gp, i_node) * N(gp, j_node);
                    tempGradPhi(i_node*num_dim + num_nodes + k_dim, j_node) -=
                        weights(gp) * N(gp, i_node) * (DN_DX[gp])(j_node, k_dim); */

                    //tempPhi: LHS for the first 4 rows
                    // tempPhi(i_node*(num_dim + 1) + 0, j_node*(num_dim + 1) + 0) +=
                    //     weights(gp) * (/* dissipative_coefficient + */ tau) * (DN_DX[gp])(i_node, k_dim) * (DN_DX[gp])(j_node, k_dim);
                    // tempPhi(i_node*(num_dim + 1) + 0, j_node*(num_dim + 1) + k_dim + 1) +=
                    //     weights(gp) * (1 - tau) * N(gp, i_node) * (DN_DX[gp])(j_node, k_dim);

                    //tempGradPhi: LHS for the next 3*4 rows
                    // tempGradPhi(i_node*(num_dim + 1) + k_dim + 1, j_node*(num_dim + 1) + k_dim + 1) +=
                    //     weights(gp) * N(gp, i_node) * N(gp, j_node);
                    // tempGradPhi(i_node*(num_dim + 1) + k_dim + 1, j_node*(num_dim + 1) + 0) -=
                    //     weights(gp) * N(gp, i_node) * (DN_DX[gp])(j_node, k_dim);

                    lhs(i_node, j_node) += weights(gp) * (DN_DX[gp])(i_node, k_dim) * (DN_DX[gp])(j_node, k_dim);

                    if (step > 1){
                        rhs[i_node] += diffusion * weights(gp) * (DN_DX[gp])(i_node, k_dim) * N(gp, j_node) * grad_phi_avg[k_dim];//grad_phi[k_dim];//
                    }
                }
            }
        }
    }

    //KRATOS_INFO("VariationalNonEikonalDistanceElement") << "Here 3" << std::endl;

    if (npos != 0 && nneg != 0)
    {
        //KRATOS_INFO("VariationalNonEikonalDistanceElement, nneg, npos") << nneg << ", " << npos << std::endl;
        //KRATOS_INFO("VariationalNonEikonalDistanceElement, distances0")
        //    << "0: " << distances0(0)
        //    << ", 1: " << distances0(1)
        //    << ", 2: " << distances0(2)
        //    << ", 3: " << distances0(3)
        //    << std::endl;

        ModifiedShapeFunctions::Pointer p_modified_sh_func =
            Kratos::make_shared<Tetrahedra3D4ModifiedShapeFunctions>(p_geometry, distances0);

        Matrix int_N;
        GeometryType::ShapeFunctionsGradientsType int_DN_DX;
        Vector int_weights;

        if (step == 1){
            //KRATOS_WATCH(nneg)
            Matrix neg_N, pos_N;
            GeometryType::ShapeFunctionsGradientsType neg_DN_DX, pos_DN_DX;
            Vector neg_weights, pos_weights;

            p_modified_sh_func->ComputePositiveSideShapeFunctionsAndGradientsValues(
                pos_N,        // N
                pos_DN_DX,    // DN_DX
                pos_weights,  // weight * detJ
                integration_method);

            const std::size_t number_of_pos_gauss_points = pos_weights.size();

            for (unsigned int pos_gp = 0; pos_gp < number_of_pos_gauss_points; pos_gp++){
                for (unsigned int i_node = 0; i_node < num_nodes; i_node++){
                    //rhs(i_node*(num_dim + 1) + 0) += pos_weights(pos_gp) * pos_N(pos_gp, i_node);
                    rhs(i_node) += 0.0;//source_coeff * pos_weights(pos_gp) * pos_N(pos_gp, i_node);
                    // There is no need for the positive source term
                }
            }

            p_modified_sh_func->ComputeNegativeSideShapeFunctionsAndGradientsValues(
                neg_N,        // N
                neg_DN_DX,    // DN_DX
                neg_weights,  // weight * detJ
                integration_method);

            const std::size_t number_of_neg_gauss_points = neg_weights.size();

            for (unsigned int neg_gp = 0; neg_gp < number_of_neg_gauss_points; neg_gp++){
                for (unsigned int i_node = 0; i_node < num_nodes; i_node++){
                    //rhs(i_node*(num_dim + 1) + 0) -= neg_weights(neg_gp) * neg_N(neg_gp, i_node);
                    rhs(i_node) -= 1.0e0 * source_coeff * neg_weights(neg_gp) * neg_N(neg_gp, i_node);
                }
            }
        } //else{

            //KRATOS_INFO("VariationalNonEikonalDistanceElement") << "Here 4" << std::endl;

            p_modified_sh_func->ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
                int_N,
                int_DN_DX,
                int_weights,
                integration_method);

            const std::size_t number_of_int_gauss_points = int_weights.size();

            for (unsigned int int_gp = 0; int_gp < number_of_int_gauss_points; int_gp++){
                for (unsigned int i_node = 0; i_node < num_nodes; i_node++){
                    for (unsigned int j_node = 0; j_node < num_nodes; j_node++){
                        //tempPhi: LHS for 4 rows associated with Phi
                        // tempPhi(i_node*(num_dim + 1) + 0, j_node*(num_dim + 1) + 0) +=
                        //     penalty_phi0 * int_weights(int_gp) * int_N(int_gp, i_node) * int_N(int_gp, j_node);

                        /* for (unsigned int k_dim = 0; k_dim < num_dim; k_dim++){
                            //tempGradPhi: LHS for 3*4 rows associated with GradPhi
                            // tempGradPhi(i_node*(num_dim + 1) + k_dim + 1, j_node*(num_dim + 1) + k_dim + 1) +=
                            //     penalty_curvature * int_weights(int_gp) * int_N(int_gp, i_node) * (int_DN_DX[int_gp])(j_node,k_dim);

                            // rhs(i_node*(num_dim + 1) + k_dim + 1) +=
                            //     penalty_curvature*int_weights(int_gp)*(int_DN_DX[int_gp])(j_node,k_dim)*distance_gradient0(j_node, k_dim)*int_N(int_gp, i_node);

                        } */

                        lhs(i_node, j_node) += penalty_phi0*int_weights(int_gp)*int_N(int_gp, i_node)*int_N(int_gp, j_node);
                    }
                }
            }

            std::vector<unsigned int> contact_line_faces;
            std::vector<unsigned int> contact_line_indices;

            std::vector<MatrixType> contact_N;                                   //std::vector for multiple contact lines
            std::vector<GeometryType::ShapeFunctionsGradientsType> contact_DN_DX;//std::vector for multiple contact lines
            std::vector<Kratos::Vector> contact_weights;                         //std::vector for multiple contact lines

            auto& neighbour_elems = this->GetValue(NEIGHBOUR_ELEMENTS);

            for (unsigned int i_cl = 0; i_cl < contact_line_faces.size(); i_cl++){
                if (neighbour_elems[ contact_line_faces[i_cl] ].Id() == this->Id() ){
                    contact_line_indices.push_back(i_cl);
                }
            }

            // Call the Contact Line negative side shape functions calculator
            p_modified_sh_func->ComputeContactLineNegativeSideShapeFunctionsAndGradientsValues(
                contact_line_indices, //ADDED
                contact_N,
                contact_DN_DX,
                contact_weights,
                GeometryData::GI_GAUSS_2);

            for (unsigned int i_cl = 0; i_cl < contact_weights.size(); i_cl++){
                const std::size_t number_of_contact_gauss_points = (contact_weights[i_cl]).size();
                for (unsigned int contact_gp = 0; contact_gp < number_of_contact_gauss_points; contact_gp++){
                    for (unsigned int i_node = 0; i_node < num_nodes; i_node++){
                        for (unsigned int j_node = 0; j_node < num_nodes; j_node++){

                            lhs(i_node, j_node) += 0.0;
                                1.0e0 * penalty_phi0*(contact_weights[i_cl])(contact_gp)*(contact_N[i_cl])(contact_gp, i_node)*(contact_N[i_cl])(contact_gp, j_node);

                        }
                    }
                }
            }
        //}

    } else if (step == 1){
        //KRATOS_WATCH(nneg)
        double source;
        if (npos != 0)
            source = 0.0;//1.0; // There is no need to add positive source term
        else
            source = -1.0e0;

        for (unsigned int gp = 0; gp < number_of_gauss_points; gp++){
            for (unsigned int i_node = 0; i_node < num_nodes; i_node++){
                //rhs(i_node*(num_dim + 1) + 0) += weights(gp) * source * N(gp, i_node);
                rhs(i_node) += source_coeff * source * weights(gp) * N(gp, i_node);
            }
        }
    }

    if (step == 1){
        GeometryType::GeometriesArrayType faces = GetGeometry().GenerateFaces();

        unsigned int i_face = 0;
        const unsigned int num_faces = num_nodes; //for simplex elements
        const unsigned int num_face_nodes = num_nodes - 1;

        while (i_face < num_faces) {

            GeometryType& r_face = faces[i_face];
            unsigned int boundary_face_node = 0;

            double contact_angle = 0.0;
            double contact_angle_weight = 0.0;
            Vector solid_normal = ZeroVector(num_dim);

            for (unsigned int i=0; i < num_face_nodes; ++i){
                if ( r_face[i].Is(BOUNDARY) ){
                    boundary_face_node++;
                    const double contact_angle_i = r_face[i].FastGetSolutionStepValue(CONTACT_ANGLE);
                    if (contact_angle_i > 1.0e-12)
                    {
                        contact_angle += contact_angle_i;
                        contact_angle_weight += 1.0;
                    }
                    solid_normal += r_face[i].FastGetSolutionStepValue(NORMAL);
                }
            }

            if (boundary_face_node == num_face_nodes){
                double minus_cos_contact_angle = 0.0;
                const double norm_solid_normal = Kratos::norm_2(solid_normal);
                solid_normal = (1.0/norm_solid_normal)*solid_normal;

                const unsigned int num_int_pts = (faces[i_face]).IntegrationPointsNumber(integration_method);

                //std::vector< IntegrationPoint<3> > face_gauss_pts;
                auto face_gauss_pts = (faces[i_face]).IntegrationPoints(integration_method);

                VectorType face_jacobian;
                (faces[i_face]).DeterminantOfJacobian(face_jacobian, integration_method);

                // Get the original geometry shape function and gradients values over the intersection
                for (unsigned int i_gauss = 0; i_gauss < num_int_pts; ++i_gauss) {
                    // Store the Gauss points weights
                    const double face_weight = face_jacobian(i_gauss) * face_gauss_pts[i_gauss].Weight();

                    // Compute the global coordinates of the face Gauss pt.
                    GeometryType::CoordinatesArrayType global_coords = ZeroVector(num_dim);
                    global_coords = (faces[i_face]).GlobalCoordinates(global_coords, face_gauss_pts[i_gauss].Coordinates());

                    // Compute the parent geometry local coordinates of the Gauss pt.
                    GeometryType::CoordinatesArrayType loc_coords = ZeroVector(num_dim);
                    loc_coords = GetGeometry().PointLocalCoordinates(loc_coords, global_coords);

                    // // Compute shape function values
                    // // Obtain the parent subgeometry shape function values

                    Kratos::Vector face_shape_func = ZeroVector(num_nodes);
                    face_shape_func = GetGeometry().ShapeFunctionsValues(face_shape_func, loc_coords);

            //         /* double det_jac;
            //         VectorType face_jacobian_inv;
            //         MathUtils<double>::GeneralizedInvertMatrix( face_jacobian(i_gauss), face_jacobian_inv, det_jac );
            //         Matrix& face_shape_func_derivative = =  prod( GetGeometry().ShapeFunctionsLocalGradients(face_shape_func_derivative, loc_coords),
            //             face_jacobian_inv); */ // This code is not working, it is only as a hint for further development

                    for (unsigned int i_node = 0; i_node < num_nodes; i_node++){
                        if (contact_angle_weight > 0.0){
                            minus_cos_contact_angle = -std::cos(contact_angle/contact_angle_weight);
                        } else{
                            const double norm_grad_phi_old = norm_2(grad_phi_old); // For linear elements, grad_phi_old is constant, however, this is not the correct way to calculate it
                            minus_cos_contact_angle = Kratos::inner_prod(solid_normal,1.0/norm_grad_phi_old*grad_phi_old);
                                //GetGeometry()[i_node].FastGetSolutionStepValue(DISTANCE_GRADIENT));
                        }

                        rhs(i_node) += minus_cos_contact_angle * face_weight * face_shape_func(i_node);
                        //In case we use this process for the sake of redistancing, this RHS contribution should be turned off.
                    }
                }
            }
            i_face++;
        }
    }

    noalias(rLeftHandSideMatrix) = lhs;// tempPhi + tempGradPhi;
    noalias(rRightHandSideVector) = rhs - prod(rLeftHandSideMatrix,values); // Reducing the contribution of the last known values

    KRATOS_CATCH("");
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental left hand side matrix only
 * @param rLeftHandSideMatrix: the elemental left hand side matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
void VariationalNonEikonalDistanceElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental right hand side vector only
 * @param rRightHandSideVector: the elemental right hand side vector
 * @param rCurrentProcessInfo: the current process info instance
 */
void VariationalNonEikonalDistanceElement::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
}

/**
 * this is called during the assembling process in order
 * to calculate the first derivatives contributions for the LHS and RHS
 * @param rLeftHandSideMatrix: the elemental left hand side matrix
 * @param rRightHandSideVector: the elemental right hand side
 * @param rCurrentProcessInfo: the current process info instance
 */
void VariationalNonEikonalDistanceElement::CalculateFirstDerivativesContributions(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
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
void VariationalNonEikonalDistanceElement::CalculateFirstDerivativesLHS(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
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
void VariationalNonEikonalDistanceElement::CalculateFirstDerivativesRHS(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
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
void VariationalNonEikonalDistanceElement::CalculateSecondDerivativesContributions(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
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
void VariationalNonEikonalDistanceElement::CalculateSecondDerivativesLHS(
    MatrixType& rLeftHandSideMatrix,
    ProcessInfo& rCurrentProcessInfo)
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
void VariationalNonEikonalDistanceElement::CalculateSecondDerivativesRHS(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
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
void VariationalNonEikonalDistanceElement::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
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
void VariationalNonEikonalDistanceElement::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo)
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
int VariationalNonEikonalDistanceElement::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //KRATOS_INFO("VariationalNonEikonalDistanceElement") << "Check" << std::endl;

    KRATOS_ERROR_IF(this->Id() < 1) <<"VariationalNonEikonalDistanceElement found with Id 0 or negative" << std::endl;

    KRATOS_ERROR_IF(this->GetGeometry().Area() <= 0) << "On VariationalNonEikonalDistanceElement -> "
        << this->Id() <<  "; Area cannot be less than or equal to 0" << std::endl;

      // Base class checks for positive Jacobian and Id > 0
      int ierr = Element::Check(rCurrentProcessInfo);
      if(ierr != 0) return ierr;

      // Check that all required variables have been registered
      KRATOS_CHECK_VARIABLE_KEY(DISTANCE)
      KRATOS_CHECK_VARIABLE_KEY(DISTANCE_AUX2)
      /* KRATOS_CHECK_VARIABLE_KEY(DISTANCE_GRADIENT_X)
      KRATOS_CHECK_VARIABLE_KEY(DISTANCE_GRADIENT_Y)
      KRATOS_CHECK_VARIABLE_KEY(DISTANCE_GRADIENT_Z) */

      unsigned const int number_of_points = GetGeometry().size();
      // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
      for ( unsigned int i = 0; i < number_of_points; i++ )
      {
          Node<3> &rnode = this->GetGeometry()[i];
          /* KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISTANCE,rnode)
          KRATOS_CHECK_DOF_IN_NODE(DISTANCE,rnode) */
          KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISTANCE_AUX2,rnode)
          KRATOS_CHECK_DOF_IN_NODE(DISTANCE_AUX2,rnode)
          /* KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISTANCE_GRADIENT_X,rnode)
          KRATOS_CHECK_DOF_IN_NODE(DISTANCE_GRADIENT_X,rnode)
          KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISTANCE_GRADIENT_Y,rnode)
          KRATOS_CHECK_DOF_IN_NODE(DISTANCE_GRADIENT_Y,rnode)
          KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISTANCE_GRADIENT_Z,rnode)
          KRATOS_CHECK_DOF_IN_NODE(DISTANCE_GRADIENT_Z,rnode) */
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

std::string VariationalNonEikonalDistanceElement::Info() const {
    std::stringstream buffer;
    buffer << "VariationalNonEikonalDistanceElement #" << Id();
    return buffer.str();
}

/// Print information about this object.

void VariationalNonEikonalDistanceElement::PrintInfo(std::ostream& rOStream) const {
    rOStream << "VariationalNonEikonalDistanceElement #" << Id();
}

/// Print object's data.

void VariationalNonEikonalDistanceElement::PrintData(std::ostream& rOStream) const {
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

void VariationalNonEikonalDistanceElement::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );

    // List
    // To be completed with the class member list
}

void VariationalNonEikonalDistanceElement::load(Serializer& rSerializer)
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
inline std::istream & operator >> (std::istream& rIStream, VariationalNonEikonalDistanceElement& rThis);

/// output stream function
inline std::ostream & operator << (std::ostream& rOStream, const VariationalNonEikonalDistanceElement& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.