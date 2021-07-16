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
#include "custom_elements/surface_smoothing_element.h"

#define PI 3.14159265358979

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
SurfaceSmoothingElement::SurfaceSmoothingElement(IndexType NewId)
    : Element(NewId) 
{
}

/**
 * Constructor using an array of nodes
 */
SurfaceSmoothingElement::SurfaceSmoothingElement(IndexType NewId, const NodesArrayType& ThisNodes)
    : Element(NewId, ThisNodes) 
{
}

/**
 * Constructor using Geometry
 */
SurfaceSmoothingElement::SurfaceSmoothingElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry) 
{
}

/**
 * Constructor using Properties
 */
SurfaceSmoothingElement::SurfaceSmoothingElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties) 
{
}

/**
 * Copy Constructor
 */
SurfaceSmoothingElement::SurfaceSmoothingElement(SurfaceSmoothingElement const& rOther)
    : Element(rOther) 
{
}

/**
 * Destructor
 */
SurfaceSmoothingElement::~SurfaceSmoothingElement()
{
}

///@}
///@name Operators
///@{

/// Assignment operator.
SurfaceSmoothingElement & SurfaceSmoothingElement::operator=(SurfaceSmoothingElement const& rOther)
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
Element::Pointer SurfaceSmoothingElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<SurfaceSmoothingElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

/**
 * creates a new element pointer
 * @param NewId: the ID of the new element
 * @param pGeom: the geometry to be employed
 * @param pProperties: the properties assigned to the new element
 * @return a Pointer to the new element
 */
Element::Pointer SurfaceSmoothingElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<SurfaceSmoothingElement>(NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

/**
 * creates a new element pointer and clones the previous element data
 * @param NewId: the ID of the new element
 * @param ThisNodes: the nodes of the new element
 * @param pProperties: the properties assigned to the new element
 * @return a Pointer to the new element
 */
Element::Pointer SurfaceSmoothingElement::Clone(IndexType NewId, NodesArrayType const& ThisNodes) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<SurfaceSmoothingElement>(NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    KRATOS_CATCH("");
}

/**
 * this determines the elemental equation ID vector for all elemental
 * DOFs
 * @param rResult: the elemental equation ID vector
 * @param rCurrentProcessInfo: the current process info instance
 */
void SurfaceSmoothingElement::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    const int num_dim  = 3;
    const int num_nodes  = num_dim + 1;

    // num_dof = num_nodes
    if (rResult.size() != num_nodes){
        rResult.resize(num_nodes, false);
    }

    for(unsigned int i=0; i<num_nodes; i++){
        rResult[i]  =  this->GetGeometry()[i].GetDof(DISTANCE_AUX).EquationId();
    }
}

/**
 * determines the elemental list of DOFs
 * @param ElementalDofList: the list of DOFs
 * @param rCurrentProcessInfo: the current process info instance
 */
void SurfaceSmoothingElement::GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& CurrentProcessInfo)
{
    const int num_dim  = 3;
    const int num_nodes  = num_dim + 1;

    // num_dof = num_nodes
    if (rElementalDofList.size() != num_nodes){
        rElementalDofList.resize(num_nodes);
    }

    for(unsigned int i=0; i<num_nodes; i++){
        rElementalDofList[i] = this->GetGeometry()[i].pGetDof(DISTANCE_AUX);
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
void SurfaceSmoothingElement::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const int num_dim  = 3;
    const int num_nodes  = num_dim + 1;
    const unsigned int num_faces = num_nodes; //Simplex element

    const double dt = rCurrentProcessInfo.GetValue(DELTA_TIME);

    //GeometryType::Pointer p_geom = this->pGetGeometry();
    //const double he = ElementSizeCalculator<3,4>::AverageElementSize(*p_geom);
    //const double he = ElementSizeCalculator<3,4>::AverageElementSize(GetGeometry()); //this->GetValue(ELEMENT_H);
    //const double epsilon = 2.0e3*dt*he*he;//1.0e0*dt*he;//1.0e4*dt*he*he;
    //KRATOS_INFO("smoothing coefficient:") << epsilon << std::endl;

    BoundedMatrix<double,num_nodes,num_dim> DN_DX;  // Gradients matrix 
    array_1d<double,num_nodes> N; //dimension = number of nodes . Position of the gauss point 

    array_1d<double,num_nodes> PHIdof; //dimension = number of DOFs . . since we are using a residualbased approach
    array_1d<double,num_nodes> PHIold; //dimension = number of DOFs . . since we are using a residualbased approach
    array_1d<VectorType,num_nodes> GradPHIold; 

    BoundedMatrix<double,num_nodes,num_nodes> tempM;
    tempM = ZeroMatrix(num_nodes,num_nodes);
    BoundedMatrix<double,num_nodes,num_nodes> tempMlumped;
    tempMlumped = ZeroMatrix(num_nodes,num_nodes);
    BoundedMatrix<double,num_nodes,num_nodes> tempA;
    tempA = ZeroMatrix(num_nodes,num_nodes);

    array_1d<double,num_nodes> tempLaplacianRHS = ZeroVector(num_nodes); 
    array_1d<double,num_nodes> tempBCRHS = ZeroVector(num_nodes); 

    // num_dof = num_nodes
    if(rLeftHandSideMatrix.size1() != num_nodes)
        rLeftHandSideMatrix.resize(num_nodes,num_nodes,false); //resizing the system in case it does not have the right size 

    if(rRightHandSideVector.size() != num_nodes)
        rRightHandSideVector.resize(num_nodes,false);

    // Getting data for the given geometry
    double area;
    GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, area); //asking for gradients and other info
    const double he = ElementSizeCalculator<3,4>::GradientsElementSize(DN_DX);
    const double epsilon = 5.0e2*dt*he*he;//1.0e0*dt*he;//1.0e4*dt*he*he;

    const double zeta = 5.0e-1;//1.0;//0.7;//
    const double gamma = 0.072;//0.0426;//0.0311;//
    const double micro_length_scale = 1.0e-9;

    const double theta_advancing = 149.0*PI/180.0;
    const double theta_receding = 115.0*PI/180.0;
    // const double cos_theta_s = -0.4539905;///* 0.5299192642332 */-0.25881904510252076;//0.779337965;//
    // const double theta_s = std::acos(cos_theta_s);

    // Main loop (one Gauss point)
    //const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();
    //noalias(rLeftHandSideMatrix) = prod(DN_DX, Matrix(prod(D, trans(DN_DX))));  // Bt D B

    // Subtracting the dirichlet term
    // RHS -= LHS*DUMMY_UNKNOWNs

    /* bool has_solid_contact = false;
    VectorType solid_normal = ZeroVector(num_dim);
    for(unsigned int i = 0; i<num_nodes; i++){
        if (GetGeometry()[i].GetValue(IS_STRUCTURE) == 1.0)
        {
            has_solid_contact = true;
            solid_normal = GetGeometry()[i].FastGetSolutionStepValue(NORMAL);
            const double norm = Kratos::norm_2(solid_normal);
            solid_normal = (1.0/norm)*solid_normal;
            //KRATOS_INFO("Smoothing, found contact") << solid_normal << std::endl;
            break;
        }
    } 

    VectorType n_dot_grad = ZeroVector(num_nodes);
    if (has_solid_contact){
        for(unsigned int i = 0; i<num_nodes; i++){
            for (unsigned int k = 0; k<num_dim; k++){
                n_dot_grad(i) += solid_normal(k)*DN_DX(i,k);
            }
        }
    } */

    VectorType grad_phi_old = ZeroVector(num_dim);
    for(unsigned int i = 0; i<num_nodes; i++){
            for (unsigned int k = 0; k<num_dim; k++){
                grad_phi_old(k) += GetGeometry()[i].FastGetSolutionStepValue(DISTANCE)*DN_DX(i,k);
            }
    }

    for(unsigned int i = 0; i<num_nodes; i++){
        PHIold[i] = GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
        PHIdof[i] = GetGeometry()[i].FastGetSolutionStepValue(DISTANCE_AUX);
        GradPHIold[i] = grad_phi_old;//GetGeometry()[i].FastGetSolutionStepValue(DISTANCE_GRADIENT);
    }

    for(unsigned int i = 0; i<num_nodes; i++){
        tempMlumped(i,i) = area*N[i];

        for(unsigned int j = 0; j<num_nodes; j++){
            tempM(i,j) = area*N[i]*N[j];

            for (unsigned int k = 0; k<num_dim; k++){
                tempA(i,j) += area*epsilon*DN_DX(i,k)*DN_DX(j,k);

                tempLaplacianRHS[i] += area*epsilon*DN_DX(i,k)*N[j]*(GradPHIold[j])[k];
            }

            //tempA(i,j) -= 0.5*area*epsilon*n_dot_grad(i)*n_dot_grad(j);

        }
    }

    //KRATOS_INFO("SurfaceSmoothingElement") << "Start BC" << std::endl;
    ///////////////////////////////////////////////////////////////////////////////
    GeometryType::GeometriesArrayType faces = GetGeometry().GenerateFaces();

    bool not_found_surface = true;
    unsigned int i_face = 0;

    while (i_face < num_faces && not_found_surface) {
        GeometryType& r_face = faces[i_face];
        unsigned int contact_node = 0;
        const unsigned int num_face_nodes = num_nodes - 1;
        unsigned int n_pos = 0;

        double positive_viscosity = 0.0;
        double negative_viscosity = 0.0;

        for (unsigned int i=0; i < num_face_nodes; ++i){
            if ( r_face[i].GetValue(IS_STRUCTURE) == 1.0 ){
                contact_node++;
            }
            if ( r_face[i].FastGetSolutionStepValue(DISTANCE) > 0.0 ){
                n_pos++;
                positive_viscosity = r_face[i].FastGetSolutionStepValue(DYNAMIC_VISCOSITY);
            } else{
                negative_viscosity = r_face[i].FastGetSolutionStepValue(DYNAMIC_VISCOSITY);
            }
        }

        if (contact_node == num_face_nodes){
            not_found_surface = false;

            //MatrixType FaceShapeFunctions;
            //GeometryType::ShapeFunctionsGradientsType FaceShapeFunctionsGradients;
            //Kratos::Vector FaceWeights;

            //VectorType solid_normal = ZeroVector(num_dim);
            //VectorType distance_grad = ZeroVector(num_dim);

            auto IntegrationMethod = GeometryData::GI_GAUSS_1;
            const unsigned int num_int_pts = (faces[i_face]).IntegrationPointsNumber(IntegrationMethod);

            //FaceShapeFunctions.resize(n_int_pts, n_nodes, false);
            //FaceShapeFunctionsGradients.resize(n_int_pts, false);
            //FaceWeights.resize(n_int_pts, false);

            std::vector < GeometryType::CoordinatesArrayType > face_gauss_pts_gl_coords, face_gauss_pts_loc_coords;
            face_gauss_pts_gl_coords.clear();
            face_gauss_pts_loc_coords.clear();
            face_gauss_pts_gl_coords.reserve(num_int_pts);
            face_gauss_pts_loc_coords.reserve(num_int_pts);

            //std::vector< IntegrationPoint<3> > face_gauss_pts;
            auto face_gauss_pts = (faces[i_face]).IntegrationPoints(IntegrationMethod);

            VectorType face_jacobians;
            (faces[i_face]).DeterminantOfJacobian(face_jacobians, IntegrationMethod);

            // Get the original geometry shape function and gradients values over the intersection
            for (unsigned int i_gauss = 0; i_gauss < num_int_pts; ++i_gauss) {
                // Store the Gauss points weights
               const double face_weight = face_jacobians(i_gauss) * face_gauss_pts[i_gauss].Weight();

                // Compute the global coordinates of the face Gauss pt.
                GeometryType::CoordinatesArrayType global_coords = ZeroVector(num_dim);
                global_coords = (faces[i_face]).GlobalCoordinates(global_coords, face_gauss_pts[i_gauss].Coordinates());

                // Compute the parent geometry local coordinates of the Gauss pt.
                GeometryType::CoordinatesArrayType loc_coords = ZeroVector(num_dim);
                loc_coords = GetGeometry().PointLocalCoordinates(loc_coords, global_coords);

                // Compute shape function values
                // Obtain the parent subgeometry shape function values
                double det_jac;
                Kratos::Vector face_shape_func;
                face_shape_func = GetGeometry().ShapeFunctionsValues(face_shape_func, loc_coords);

                double temp_value = 0.0;

                for (unsigned int i = 0; i < num_nodes; i++){

                    Vector solid_normal = ZeroVector(num_dim);
                    Vector corrected_gradient = ZeroVector(num_dim);

                    if (GetGeometry()[i].GetValue(IS_STRUCTURE) == 1.0){
                        if (n_pos > 0 && n_pos < num_face_nodes){ //cut solid surface element
                            solid_normal = GetGeometry()[i].FastGetSolutionStepValue(NORMAL);
                            const double norm = Kratos::norm_2(solid_normal);
                            solid_normal = (1.0/norm)*solid_normal;

                            const double norm_grad_phi = norm_2(GradPHIold[i]);
                            const Vector normal = GradPHIold[i]/norm_grad_phi;

                            Vector contact_tangential = ZeroVector(num_dim);
                            MathUtils<double>::UnitCrossProduct(contact_tangential, normal, solid_normal);
                            Vector slip_vector = ZeroVector(num_dim);
                            MathUtils<double>::UnitCrossProduct(slip_vector, solid_normal, contact_tangential);

                            const double slip_velocity = inner_prod(slip_vector,
                                GetGeometry()[i].FastGetSolutionStepValue(VELOCITY));

                            Vector contact_vector_macro = ZeroVector(num_dim);
                            MathUtils<double>::UnitCrossProduct(contact_vector_macro, contact_tangential, normal);
                            const double cos_theta_macro = inner_prod(slip_vector,contact_vector_macro);
                            const double theta_macro = GetGeometry()[i].FastGetSolutionStepValue(CONTACT_ANGLE);//std::acos(cos_theta_macro);

                            double theta_equilibrium = theta_receding;
                            if (theta_macro > theta_equilibrium){
                                if (theta_macro >= theta_advancing){
                                    theta_equilibrium = theta_advancing;
                                } else {
                                    theta_equilibrium = theta_macro;
                                }
                            }
                            double cos_theta_equilibrium = std::cos(theta_equilibrium);

                            const double cos_theta_d = cos_theta_equilibrium - zeta/gamma * slip_velocity;//Check the sign of slip velocity

                            KRATOS_WARNING_IF("SurfaceSmooting", std::abs(cos_theta_d) > 1.0)
                                << "cos_theta_d is larger than one." << std::endl;

                            double theta_d = 0.0;
                            if (std::abs(cos_theta_d) <= 1.0){
                                theta_d = std::acos(cos_theta_d);
                            } else if (cos_theta_d > 1.0){
                                theta_d = 0.0;
                            } else { //if (cos_theta_d < -1.0){
                                theta_d = PI;
                            }

                            const double effective_viscosity = 0.5*(positive_viscosity + negative_viscosity);
                            const double capilary_number = effective_viscosity*slip_velocity/gamma;

                            if ( std::abs(theta_d - theta_equilibrium) < 6.0e-1 &&
                                capilary_number < 3.0e-1){

                                double contact_angle_macro = 0.0;

                                const double cubic_contact_angle_macro = std::pow(theta_d, 3.0)
                                    + 9*capilary_number*std::log(he/micro_length_scale);

                                KRATOS_WARNING_IF("SurfaceSmooting", cubic_contact_angle_macro < 0.0 ||
                                    cubic_contact_angle_macro > 31.0)
                                    << "Hydrodynamics theory failed to estimate micro contact-angle (large slip velocity)." 
                                    << std::endl;

                                if (cubic_contact_angle_macro >= 0.0 &&
                                        cubic_contact_angle_macro <= 31.0) //std::pow(PI, 3.0))
                                    contact_angle_macro = std::pow(cubic_contact_angle_macro, 1.0/3.0);
                                else if (cubic_contact_angle_macro < 0.0)
                                    contact_angle_macro = 0.0; //contact_angle_equilibrium;
                                else //if (cubic_contact_angle_micro_gp > 31.0){
                                    contact_angle_macro = PI; //contact_angle_equilibrium;
                                //}

                                const double cos_contact_angle_macro = std::cos(contact_angle_macro);
                                const double sin_contact_angle_macro = std::sqrt( 1.0 -
                                    cos_contact_angle_macro*cos_contact_angle_macro );

                                corrected_gradient = norm_grad_phi*( -cos_contact_angle_macro*solid_normal
                                    + sin_contact_angle_macro*slip_vector );

                            } /* else */ {
                                if (std::abs(cos_theta_d) <= 1.0){
                                    const double sin_theta_d = std::sqrt( 1.0 - cos_theta_d*cos_theta_d );
                                    corrected_gradient = norm_grad_phi*( -cos_theta_d*solid_normal + sin_theta_d*slip_vector );
                                } else if (cos_theta_d > 1.0){
                                    corrected_gradient = -norm_grad_phi*solid_normal;
                                } else //if (cos_theta_d < -1.0){
                                    corrected_gradient = norm_grad_phi*solid_normal;
                                //}
                            }

                            //corrected_gradient = GradPHIold[i];

                        } else { //not a cut solid surface element
                            corrected_gradient = GradPHIold[i];
                        }
                    }

                    temp_value += face_shape_func(i)*Kratos::inner_prod(solid_normal, corrected_gradient);
                }

                for (unsigned int i = 0; i < num_nodes; i++){
                    tempBCRHS[i] += epsilon * temp_value * face_weight * face_shape_func(i);
                }
            }
        }

        i_face++;
    }
    ///////////////////////////////////////////////////////////////////////////////
    //KRATOS_INFO("SurfaceSmoothingElement") << "End BC" << std::endl;

    //const unsigned int step = rCurrentProcessInfo[FRACTIONAL_STEP];

    //if (step == 1){
        noalias(rLeftHandSideMatrix) = tempM + tempA; // + tempMlumped;
        noalias(rRightHandSideVector) = /* -tempLaplacianRHS + */ tempBCRHS + prod(tempM,PHIold) - prod(rLeftHandSideMatrix,PHIdof);
    //} else {
    //    noalias(rLeftHandSideMatrix) = tempM; //+ tempA; // + tempMlumped;
    //    noalias(rRightHandSideVector) = tempLaplacianRHS - tempBCRHS + prod(tempM,PHIold) - prod(rLeftHandSideMatrix,PHIdof);
    //}  

    KRATOS_CATCH("");
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental left hand side matrix only
 * @param rLeftHandSideMatrix: the elemental left hand side matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
void SurfaceSmoothingElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental right hand side vector only
 * @param rRightHandSideVector: the elemental right hand side vector
 * @param rCurrentProcessInfo: the current process info instance
 */
void SurfaceSmoothingElement::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
}

/**
 * this is called during the assembling process in order
 * to calculate the first derivatives contributions for the LHS and RHS
 * @param rLeftHandSideMatrix: the elemental left hand side matrix
 * @param rRightHandSideVector: the elemental right hand side
 * @param rCurrentProcessInfo: the current process info instance
 */
void SurfaceSmoothingElement::CalculateFirstDerivativesContributions(
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
void SurfaceSmoothingElement::CalculateFirstDerivativesLHS(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
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
void SurfaceSmoothingElement::CalculateFirstDerivativesRHS(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
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
void SurfaceSmoothingElement::CalculateSecondDerivativesContributions(
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
void SurfaceSmoothingElement::CalculateSecondDerivativesLHS(
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
void SurfaceSmoothingElement::CalculateSecondDerivativesRHS(
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
void SurfaceSmoothingElement::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
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
void SurfaceSmoothingElement::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo)
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
int SurfaceSmoothingElement::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(this->Id() < 1) <<"SurfaceSmoothingElement found with Id 0 or negative" << std::endl;

    KRATOS_ERROR_IF(this->GetGeometry().Area() <= 0) << "On SurfaceSmoothingElement -> "
        << this->Id() <<  "; Area cannot be less than or equal to 0" << std::endl;
  
      // Base class checks for positive Jacobian and Id > 0
      int ierr = Element::Check(rCurrentProcessInfo);
      if(ierr != 0) return ierr;
  
      // Check that all required variables have been registered
      KRATOS_CHECK_VARIABLE_KEY(DISTANCE)
  
      unsigned const int number_of_points = GetGeometry().size(); 
      // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
      for ( unsigned int i = 0; i < number_of_points; i++ )
      {
          Node<3> &rnode = this->GetGeometry()[i];
          KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISTANCE,rnode)
          KRATOS_CHECK_DOF_IN_NODE(DISTANCE,rnode)
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

std::string SurfaceSmoothingElement::Info() const {
    std::stringstream buffer;
    buffer << "SurfaceSmoothingElement #" << Id();
    return buffer.str();
}

/// Print information about this object.

void SurfaceSmoothingElement::PrintInfo(std::ostream& rOStream) const {
    rOStream << "SurfaceSmoothingElement #" << Id();
}

/// Print object's data.

void SurfaceSmoothingElement::PrintData(std::ostream& rOStream) const {
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

void SurfaceSmoothingElement::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );

    // List
    // To be completed with the class member list
}

void SurfaceSmoothingElement::load(Serializer& rSerializer)
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
inline std::istream & operator >> (std::istream& rIStream, SurfaceSmoothingElement& rThis);

/// output stream function
inline std::ostream & operator << (std::ostream& rOStream, const SurfaceSmoothingElement& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.