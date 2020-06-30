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
#include "custom_elements/surface_smoothing_element.h"

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

    constexpr static std::array<unsigned int, 4> mNode0ID = {2, 0, 0, 0};
    constexpr static std::array<unsigned int, 4> mNode1ID = {3, 3, 1, 2};
    constexpr static std::array<unsigned int, 4> mNode2ID = {1, 2, 3, 1};

    constexpr static unsigned int num_dim  = 3;
    constexpr static unsigned int num_nodes  = num_dim + 1;
    constexpr static unsigned int num_faces = num_nodes;
    constexpr static unsigned int num_face_nodes = num_nodes - 1;

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
void SurfaceSmoothingElement::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& CurrentProcessInfo) const
{
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
void SurfaceSmoothingElement::GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& CurrentProcessInfo) const
{
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
void SurfaceSmoothingElement::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const double dt = rCurrentProcessInfo.GetValue(DELTA_TIME);
    const auto& geometry = this->GetGeometry();
    const double he = ElementSizeCalculator<3,4>::AverageElementSize(geometry);
    const double epsilon = 1.0e2*dt*he*he;

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
        if (neighbour_elems[ i_ne ].Id() == this->Id() ){
            auto outer_face = Triangle3D3< GeometryType::PointType >(
                                    geometry.pGetPoint(mNode0ID[i_ne]),
                                    geometry.pGetPoint(mNode1ID[i_ne]),
                                    geometry.pGetPoint(mNode2ID[i_ne]));

            unsigned int contact_node = 0;
            for (unsigned int i=0; i < num_face_nodes; ++i){
                if ( outer_face[i].Is(CONTACT) ){
                    contact_node++;
                }
            }

            if (contact_node == num_face_nodes){
                auto IntegrationMethod = GeometryData::GI_GAUSS_1;
                auto face_gauss_pts = outer_face.IntegrationPoints(IntegrationMethod);
                const unsigned int num_int_pts = face_gauss_pts.size();

                std::vector < GeometryType::CoordinatesArrayType > face_gauss_pts_gl_coords, face_gauss_pts_loc_coords;
                face_gauss_pts_gl_coords.clear();
                face_gauss_pts_loc_coords.clear();
                face_gauss_pts_gl_coords.reserve(num_int_pts);
                face_gauss_pts_loc_coords.reserve(num_int_pts);

                VectorType face_jacobians;
                outer_face.DeterminantOfJacobian(face_jacobians, IntegrationMethod);

                Kratos::Vector face_shape_func;
                GeometryType::CoordinatesArrayType global_coords = ZeroVector(num_dim);
                GeometryType::CoordinatesArrayType loc_coords = ZeroVector(num_dim);
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
void SurfaceSmoothingElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental right hand side vector only
 * @param rRightHandSideVector: the elemental right hand side vector
 * @param rCurrentProcessInfo: the current process info instance
 */
void SurfaceSmoothingElement::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
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
int SurfaceSmoothingElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    KRATOS_ERROR_IF(this->GetGeometry().Area() <= 0) << "On SurfaceSmoothingElement -> "
        << this->Id() <<  "; Area cannot be less than or equal to 0" << std::endl;

    // Base class checks for domain size and Id > 0
    int ierr = Element::Check(rCurrentProcessInfo);
    if(ierr != 0) return ierr;

    // Check that all required variables have been registered
    KRATOS_CHECK_VARIABLE_KEY(DISTANCE)
    KRATOS_CHECK_VARIABLE_KEY(DISTANCE_GRADIENT)

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
