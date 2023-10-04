//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Mohammad R. Hashemi
//

// Project includes
#include "utilities/element_size_calculator.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "droplet_dynamics_application_variables.h"
#include "custom_elements/droplet_dynamics_distance_smoothing_element.h"

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

    constexpr static std::array<unsigned int, 4> mNode0ID3D = {2, 0, 0, 0};
    constexpr static std::array<unsigned int, 4> mNode1ID3D = {3, 3, 1, 2};
    constexpr static std::array<unsigned int, 4> mNode2ID3D = {1, 2, 3, 1};

    constexpr static std::array<unsigned int, 3> mNode0ID2D = {1, 2, 0};
    constexpr static std::array<unsigned int, 3> mNode1ID2D = {2, 0, 1};

/**
 * Constructor.
 */
template< unsigned int TDim >
DistanceSmoothingElement<TDim>::DistanceSmoothingElement(IndexType NewId)
    : Element(NewId)
{
}

/**
 * Constructor using an array of nodes
 */
template< unsigned int TDim >
DistanceSmoothingElement<TDim>::DistanceSmoothingElement(IndexType NewId, const NodesArrayType& ThisNodes)
    : Element(NewId, ThisNodes)
{
}

/**
 * Constructor using Geometry
 */
template< unsigned int TDim >
DistanceSmoothingElement<TDim>::DistanceSmoothingElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
}

/**
 * Constructor using Properties
 */
template< unsigned int TDim >
DistanceSmoothingElement<TDim>::DistanceSmoothingElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
}

/**
 * Copy Constructor
 */
template< unsigned int TDim >
DistanceSmoothingElement<TDim>::DistanceSmoothingElement(DistanceSmoothingElement const& rOther)
    : Element(rOther)
{
}

/**
 * Destructor
 */
template< unsigned int TDim >
DistanceSmoothingElement<TDim>::~DistanceSmoothingElement()
{
}

///@}
///@name Operators
///@{

/// Assignment operator.
template< unsigned int TDim >
DistanceSmoothingElement<TDim> & DistanceSmoothingElement<TDim>::operator=(DistanceSmoothingElement const& rOther)
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
Element::Pointer DistanceSmoothingElement<TDim>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<DistanceSmoothingElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
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
Element::Pointer DistanceSmoothingElement<TDim>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<DistanceSmoothingElement>(NewId, pGeom, pProperties);
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
Element::Pointer DistanceSmoothingElement<TDim>::Clone(IndexType NewId, NodesArrayType const& ThisNodes) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<DistanceSmoothingElement>(NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    KRATOS_CATCH("");
}

/**
 * this determines the elemental equation ID vector for all elemental
 * DOFs
 * @param rResult: the elemental equation ID vector
 * @param rCurrentProcessInfo: the current process info instance
 */
template< unsigned int TDim >
void DistanceSmoothingElement<TDim>::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& CurrentProcessInfo) const
{
    const unsigned int num_nodes = TDim + 1;

    // num_dof = num_nodes
    if (rResult.size() != num_nodes){
        rResult.resize(num_nodes, false);
    }

    for(unsigned int i=0; i<num_nodes; i++){
        rResult[i]  =  this->GetGeometry()[i].GetDof(DISTANCE).EquationId();
    }
}

/**
 * determines the elemental list of DOFs
 * @param ElementalDofList: the list of DOFs
 * @param rCurrentProcessInfo: the current process info instance
 */
template< unsigned int TDim >
void DistanceSmoothingElement<TDim>::GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& CurrentProcessInfo) const
{
    const unsigned int num_nodes = TDim + 1;

    // Here, num_dof = num_nodes
    if (rElementalDofList.size() != num_nodes){
        rElementalDofList.resize(num_nodes);
    }

    for(unsigned int i=0; i<num_nodes; i++){
        rElementalDofList[i] = this->GetGeometry()[i].pGetDof(DISTANCE);
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
template< >
void DistanceSmoothingElement<2>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr static unsigned int num_dim  = 2;
    constexpr static unsigned int num_nodes  = num_dim + 1;
    constexpr static unsigned int num_faces = num_nodes;
    constexpr static unsigned int num_face_nodes = num_nodes - 1;

    const double dt = rCurrentProcessInfo.GetValue(DELTA_TIME);
    const auto& geometry = this->GetGeometry();
    const double scale = rCurrentProcessInfo.GetValue(NODAL_AREA);

    const double element_size = ElementSizeCalculator<3,4>::AverageElementSize(this->GetGeometry()); //MinimumElementSize!!?
    const double tolerance = 1.0e-3*element_size;

    // const double he = ElementSizeCalculator<num_dim,num_nodes>::AverageElementSize(geometry);
    // const double epsilon = (rCurrentProcessInfo.GetValue(SMOOTHING_COEFFICIENT))*dt*he*he;

    BoundedMatrix<double,num_nodes,num_dim> DN_DX;  // Gradients matrix
    array_1d<double,num_nodes> N; // dimension = number of nodes . Position of the gauss point

    array_1d<double,num_nodes> PHI_dof; //dimension = number of DOFs, needed since we are using a residualbased approach
    array_1d<double,num_nodes> PHI_old; //dimension = number of DOFs
    // array_1d<VectorType,num_nodes> grad_phi_old;

    BoundedMatrix<double,num_nodes,num_nodes> tempM = ZeroMatrix(num_nodes,num_nodes);
    BoundedMatrix<double,num_nodes,num_nodes> tempA = ZeroMatrix(num_nodes,num_nodes);

    array_1d<double,num_nodes> tempBCRHS = ZeroVector(num_nodes);

    if(rLeftHandSideMatrix.size1() != num_nodes) // be carefull that num_dof = num_nodes
        rLeftHandSideMatrix.resize(num_nodes,num_nodes,false); //resizing the system in case it does not have the right size

    if(rRightHandSideVector.size() != num_nodes)
        rRightHandSideVector.resize(num_nodes,false);

    // Getting data for the given geometry
    double volume;
    GeometryUtils::CalculateGeometryData(geometry, DN_DX, N, volume); //asking for gradients and other info
     double area;
    GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, area); //asking for gradients and other info
    const double he = ElementSizeCalculator<num_dim,num_nodes>::GradientsElementSize(DN_DX);//AverageElementSize(geometry);
    const double epsilon = 1.0e-10;//(rCurrentProcessInfo.GetValue(SMOOTHING_COEFFICIENT))*dt*he*he;
    VectorType grad_phi_old = ZeroVector(num_dim);
    //VectorType grad_phi_avg = ZeroVector(num_dim);
    for(unsigned int i = 0; i<num_nodes; i++){
        for (unsigned int k = 0; k<num_dim; k++){
            grad_phi_old(k) += GetGeometry()[i].FastGetSolutionStepValue(DISTANCE)*DN_DX(i,k);
        }
    }

    for(unsigned int i = 0; i<num_nodes; i++){
        PHI_old[i] = scale*geometry[i].FastGetSolutionStepValue(DISTANCE, 1);
        PHI_dof[i] = geometry[i].FastGetSolutionStepValue(DISTANCE);
        // grad_phi_old[i] = geometry[i].FastGetSolutionStepValue(DISTANCE_GRADIENT);
    }

    for(unsigned int i = 0; i<num_nodes; i++){
        for(unsigned int j = 0; j<num_nodes; j++){
            tempM(i,j) = area*N[i]*N[j];

            for (unsigned int k = 0; k<num_dim; k++){
                tempA(i,j) += area*epsilon*DN_DX(i,k)*DN_DX(j,k);
            }
        }
    }

    // Implementing the boundary condition on distance gradient
    const auto& neighbour_elems = this->GetValue(NEIGHBOUR_ELEMENTS);

    for (unsigned int i_ne = 0; i_ne < num_faces; i_ne++){
        if(nullptr == neighbour_elems(i_ne).get()){
            auto outer_face = Line3D2< GeometryType::PointType >(
                                    geometry.pGetPoint(mNode0ID2D[i_ne]),
                                    geometry.pGetPoint(mNode1ID2D[i_ne]));

            unsigned int contact_node = 0;
            unsigned int n_pos = 0;
            double positive_viscosity = 0.0;
            double negative_viscosity = 0.0;

            double contact_angle = 0.0;
            double contact_angle_weight = 0.0;
            Vector solid_normal = ZeroVector(num_dim);
            for (unsigned int i = 0; i < num_face_nodes; ++i)
            {
                if (outer_face[i].Is(CONTACT))
                {
                    contact_node++;
                    const double contact_angle_i = outer_face[i].FastGetSolutionStepValue(CONTACT_ANGLE);
                    if (contact_angle_i > 1.0e-12)
                    {
                        contact_angle += contact_angle_i;
                        contact_angle_weight += 1.0;
                    }
                    solid_normal += outer_face[i].FastGetSolutionStepValue(NORMAL);
                }
                if (outer_face[i].FastGetSolutionStepValue(DISTANCE) > 0.0)
                {
                    n_pos++;
                    positive_viscosity = outer_face[i].FastGetSolutionStepValue(DYNAMIC_VISCOSITY);
                }
                else
                {
                    negative_viscosity = outer_face[i].FastGetSolutionStepValue(DYNAMIC_VISCOSITY);
                }
            }

            if (contact_node == num_face_nodes){

                double minus_cos_contact_angle = 0.0;
                const double norm_solid_normal = Kratos::norm_2(solid_normal);
                solid_normal = (1.0 / norm_solid_normal) * solid_normal;

                auto IntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_1;
                auto face_gauss_pts = outer_face.IntegrationPoints(IntegrationMethod);
                const unsigned int num_int_pts = face_gauss_pts.size();

                VectorType face_jacobians;
                outer_face.DeterminantOfJacobian(face_jacobians, IntegrationMethod);

                Kratos::Vector face_shape_func;
                GeometryType::CoordinatesArrayType global_coords;
                GeometryType::CoordinatesArrayType loc_coords;
                VectorType solid_normal = ZeroVector(num_dim);

                // Get the original geometry shape function and gradients values over the intersection
                for (unsigned int i_gauss = 0; i_gauss < num_int_pts; ++i_gauss) {
                    // Store the Gauss points weights
                    const double face_weight = face_jacobians(i_gauss) * face_gauss_pts[i_gauss].Weight();

                    // Compute the global coordinates of the face Gauss pt.
                    global_coords = outer_face.GlobalCoordinates(global_coords, face_gauss_pts[i_gauss].Coordinates());

                    // Compute the parent geometry local coordinates of the Gauss pt.
                    loc_coords = geometry.PointLocalCoordinates(loc_coords, global_coords);

                    // Compute shape function values
                    // Obtain the parent subgeometry shape function values
                    face_shape_func = geometry.ShapeFunctionsValues(face_shape_func, loc_coords);

                    // double temp_value = 0.0;

//                     for (unsigned int i = 0; i < num_nodes; i++){
//                         if (geometry[i].Is(CONTACT) )
//                         {
//                             solid_normal = geometry[i].FastGetSolutionStepValue(NORMAL);
//                             const double norm = Kratos::norm_2(solid_normal);
// #ifdef KRATOS_DEBUG
//                             KRATOS_WARNING_IF("DistanceSmoothingElement", norm < 1.0e-12) << "WARNING: Normal close to zero" <<std::endl;
// #endif
//                             solid_normal = (1.0/norm)*solid_normal;
//                         }

//                         temp_value += face_shape_func(i)*Kratos::inner_prod(solid_normal, grad_phi_old[i]);
//                     }

                    for (unsigned int i = 0; i < num_nodes; i++){
                        // tempBCRHS[i] += epsilon * temp_value * face_weight * face_shape_func(i);
                        VectorType grad_phi_avg_i = GetGeometry()[i].FastGetSolutionStepValue(DISTANCE_GRADIENT);
                    const double norm_grad_phi_avg_i = norm_2( grad_phi_avg_i );
                    const double distance_i = GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);

                    /* const double normalDist = 1.0/(5.0*element_size)* std::abs(distance_i);
                    double theta = 1.0;
                    if (normalDist < 1.0){
                        theta = 1.0 - std::max(0.0, std::min( normalDist*normalDist*(3.0 - 2.0*normalDist), 1.0 ));//0.0; for small unpinned hydrophobic droplet //1.0; for pinned droplet
                    } */
                    const double theta = 0.0; //std::min( std::exp( -( (std::abs(distance_i) + tolerance)/(20.0*element_size) - 1.0 ) ), 1.0); //1.0;

                    if (contact_angle_weight > 0.0){
                        minus_cos_contact_angle = -theta*norm_grad_phi_avg_i*std::cos(contact_angle/contact_angle_weight) +
                        (1.0-theta)*Kratos::inner_prod(solid_normal,grad_phi_avg_i)/* /norm_grad_phi_avg_i */;
                        /* minus_cos_contact_angle = -std::cos(contact_angle/contact_angle_weight);
                        minus_cos_contact_angle = minus_cos_contact_angle*norm_grad_phi_avg_i; */
                    } else{
                        minus_cos_contact_angle = Kratos::inner_prod(solid_normal,grad_phi_avg_i)/* /norm_grad_phi_avg_i */;
                    }

                    // if (contact_angle_weight > 0.0){
                    //     minus_cos_contact_angle = -std::cos(contact_angle/contact_angle_weight);
                    // } else{
                    //     minus_cos_contact_angle = Kratos::inner_prod(solid_normal, //grad_phi_old);
                    //         GetGeometry()[i].FastGetSolutionStepValue(DISTANCE_GRADIENT));
                    // }
                    //KRATOS_WATCH(rCurrentProcessInfo[NODAL_AREA])
                    tempBCRHS[i] += epsilon * (scale*minus_cos_contact_angle) * face_weight * face_shape_func(i);
                    }
                }
            }

        }
    }

    noalias(rLeftHandSideMatrix) = tempM + tempA;
    noalias(rRightHandSideVector) = tempBCRHS + prod(tempM,PHI_old) - prod(rLeftHandSideMatrix,PHI_dof);

    KRATOS_CATCH("");
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
void DistanceSmoothingElement<3>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr static unsigned int num_dim  = 3;
    constexpr static unsigned int num_nodes  = num_dim + 1;
    constexpr static unsigned int num_faces = num_nodes;
    constexpr static unsigned int num_face_nodes = num_nodes - 1;

    const double dt = rCurrentProcessInfo.GetValue(DELTA_TIME);
    const auto& geometry = this->GetGeometry();
    const double he = ElementSizeCalculator<num_dim,num_nodes>::AverageElementSize(geometry);
    const double epsilon = (rCurrentProcessInfo.GetValue(SMOOTHING_COEFFICIENT))*dt*he*he;

    BoundedMatrix<double,num_nodes,num_dim> DN_DX;  // Gradients matrix
    array_1d<double,num_nodes> N; // dimension = number of nodes . Position of the gauss point

    array_1d<double,num_nodes> phi_dof; //dimension = number of DOFs, needed since we are using a residualbased approach
    array_1d<double,num_nodes> phi_old; //dimension = number of DOFs
    array_1d<VectorType,num_nodes> grad_phi_old;

    BoundedMatrix<double,num_nodes,num_nodes> tempM = ZeroMatrix(num_nodes,num_nodes);
    BoundedMatrix<double,num_nodes,num_nodes> tempA = ZeroMatrix(num_nodes,num_nodes);

    array_1d<double,num_nodes> tempBCRHS = ZeroVector(num_nodes);

    if(rLeftHandSideMatrix.size1() != num_nodes) // be carefull that num_dof = num_nodes
        rLeftHandSideMatrix.resize(num_nodes,num_nodes,false); //resizing the system in case it does not have the right size

    if(rRightHandSideVector.size() != num_nodes)
        rRightHandSideVector.resize(num_nodes,false);

    // Getting data for the given geometry
    double volume;
    GeometryUtils::CalculateGeometryData(geometry, DN_DX, N, volume); //asking for gradients and other info

    for(unsigned int i = 0; i<num_nodes; i++){
        phi_old[i] = geometry[i].FastGetSolutionStepValue(DISTANCE, 1);
        phi_dof[i] = geometry[i].FastGetSolutionStepValue(DISTANCE);
        grad_phi_old[i] = geometry[i].FastGetSolutionStepValue(DISTANCE_GRADIENT);
    }

    for(unsigned int i = 0; i<num_nodes; i++){
        for(unsigned int j = 0; j<num_nodes; j++){
            tempM(i,j) = volume*N[i]*N[j];

            for (unsigned int k = 0; k<num_dim; k++){
                tempA(i,j) += volume*epsilon*DN_DX(i,k)*DN_DX(j,k);
            }
        }
    }

    // Implementing the boundary condition on distance gradient: can be implemented
    // by a custom condition for a more general case
    const auto& neighbour_elems = this->GetValue(NEIGHBOUR_ELEMENTS);

    for (unsigned int i_ne = 0; i_ne < num_faces; i_ne++){
        if(nullptr == neighbour_elems(i_ne).get()){
            auto outer_face = Triangle3D3< GeometryType::PointType >(
                                    geometry.pGetPoint(mNode0ID3D[i_ne]),
                                    geometry.pGetPoint(mNode1ID3D[i_ne]),
                                    geometry.pGetPoint(mNode2ID3D[i_ne]));

            unsigned int contact_node = 0;
            for (unsigned int i=0; i < num_face_nodes; ++i){
                if ( outer_face[i].Is(CONTACT) ){
                    contact_node++;
                }
            }

            if (contact_node == num_face_nodes){
                auto IntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_1;
                auto face_gauss_pts = outer_face.IntegrationPoints(IntegrationMethod);
                const unsigned int num_int_pts = face_gauss_pts.size();

                VectorType face_jacobians;
                outer_face.DeterminantOfJacobian(face_jacobians, IntegrationMethod);

                Kratos::Vector face_shape_func;
                GeometryType::CoordinatesArrayType global_coords;
                GeometryType::CoordinatesArrayType loc_coords;
                VectorType solid_normal = ZeroVector(num_dim);

                // Get the original geometry shape function and gradients values over the intersection
                for (unsigned int i_gauss = 0; i_gauss < num_int_pts; ++i_gauss) {
                    // Store the Gauss points weights
                    const double face_weight = face_jacobians(i_gauss) * face_gauss_pts[i_gauss].Weight();

                    // Compute the global coordinates of the face Gauss pt.
                    global_coords = outer_face.GlobalCoordinates(global_coords, face_gauss_pts[i_gauss].Coordinates());

                    // Compute the parent geometry local coordinates of the Gauss pt.
                    loc_coords = geometry.PointLocalCoordinates(loc_coords, global_coords);

                    // Compute shape function values
                    // Obtain the parent subgeometry shape function values
                    face_shape_func = geometry.ShapeFunctionsValues(face_shape_func, loc_coords);

                    double temp_value = 0.0;

                    for (unsigned int i = 0; i < num_nodes; i++){
                        if (geometry[i].Is(CONTACT) )
                        {
                            solid_normal = geometry[i].FastGetSolutionStepValue(NORMAL);
                            const double norm = Kratos::norm_2(solid_normal);
#ifdef KRATOS_DEBUG
                            KRATOS_WARNING_IF("DistanceSmoothingElement", norm < 1.0e-12) << "WARNING: Normal close to zero" <<std::endl;
#endif
                            solid_normal = (1.0/norm)*solid_normal;
                        }

                        temp_value += face_shape_func(i)*Kratos::inner_prod(solid_normal, grad_phi_old[i]);
                    }

                    for (unsigned int i = 0; i < num_nodes; i++){
                        tempBCRHS[i] += epsilon * temp_value * face_weight * face_shape_func(i);
                    }
                }
            }

        }
    }

    noalias(rLeftHandSideMatrix) = tempM + tempA;
    noalias(rRightHandSideVector) = tempBCRHS + prod(tempM,phi_old) - prod(rLeftHandSideMatrix,phi_dof);

    KRATOS_CATCH("");
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental left hand side matrix only
 * @param rLeftHandSideMatrix: the elemental left hand side matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
template< unsigned int TDim >
void DistanceSmoothingElement<TDim>::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental right hand side vector only
 * @param rRightHandSideVector: the elemental right hand side vector
 * @param rCurrentProcessInfo: the current process info instance
 */
template< unsigned int TDim >
void DistanceSmoothingElement<TDim>::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
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
void DistanceSmoothingElement<TDim>::CalculateFirstDerivativesContributions(
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
void DistanceSmoothingElement<TDim>::CalculateFirstDerivativesLHS(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
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
void DistanceSmoothingElement<TDim>::CalculateFirstDerivativesRHS(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
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
void DistanceSmoothingElement<TDim>::CalculateSecondDerivativesContributions(
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
void DistanceSmoothingElement<TDim>::CalculateSecondDerivativesLHS(
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
void DistanceSmoothingElement<TDim>::CalculateSecondDerivativesRHS(
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
void DistanceSmoothingElement<TDim>::CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo)
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
void DistanceSmoothingElement<TDim>::CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo)
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
int DistanceSmoothingElement<TDim>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    KRATOS_ERROR_IF(this->GetGeometry().Area() <= 0) << "On DistanceSmoothingElement -> "
        << this->Id() <<  "; Area cannot be less than or equal to 0" << std::endl;

    // Base class checks for domain size and Id > 0
    int ierr = Element::Check(rCurrentProcessInfo);
    if(ierr != 0) return ierr;

    unsigned const int number_of_points = GetGeometry().size();
    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for ( unsigned int i = 0; i < number_of_points; i++ )
    {
        auto& rnode = pGetGeometry()->GetPoint(i);
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

template< unsigned int TDim >
std::string DistanceSmoothingElement<TDim>::Info() const {
    std::stringstream buffer;
    buffer << "DistanceSmoothingElement #" << Id();
    return buffer.str();
}

/// Print information about this object.

template< unsigned int TDim >
void DistanceSmoothingElement<TDim>::PrintInfo(std::ostream& rOStream) const {
    rOStream << "DistanceSmoothingElement #" << Id();
}

/// Print object's data.

template< unsigned int TDim >
void DistanceSmoothingElement<TDim>::PrintData(std::ostream& rOStream) const {
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
void DistanceSmoothingElement<TDim>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );

    // List
    // To be completed with the class member list
}

template< unsigned int TDim >
void DistanceSmoothingElement<TDim>::load(Serializer& rSerializer)
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
inline std::istream & operator >> (std::istream& rIStream, DistanceSmoothingElement<TDim>& rThis);

/// output stream function
template< unsigned int TDim >
inline std::ostream & operator << (std::ostream& rOStream, const DistanceSmoothingElement<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

template class DistanceSmoothingElement<2>;
template class DistanceSmoothingElement<3>;

} // namespace Kratos.
