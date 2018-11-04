//
//   Project Name:        Kratos
//   Last Modified by:    $Author: gcasas $
//   Date:                $Date: 2016-03-12
//

#if !defined(KRATOS_COMPUTE_GRADIENT_POULIOT_2012_H_INCLUDED )
#define  KRATOS_COMPUTE_GRADIENT_POULIOT_2012_H_INCLUDED

// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "containers/array_1d.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "geometries/geometry.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"

// Application includes
#include "includes/variables.h"
#include "fluid_dynamics_application_variables.h"
#include "calculate_component_gradient_simplex_element.h"

namespace Kratos
{

///@addtogroup SwimmingDEMApplication
///@{

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

/// A post-processing element to recover the Laplacian from the velocity solution.
/**
 */
template< unsigned int TDim,
          unsigned int TNumNodes = TDim + 1 >
class KRATOS_API(SWIMMING_DEM_APPLICATION) ComputeGradientPouliot2012 : public ComputeComponentGradientSimplex<TDim, TNumNodes>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ComputeGradientPouliot2012
    KRATOS_CLASS_POINTER_DEFINITION(ComputeGradientPouliot2012);

    typedef ComputeComponentGradientSimplex<TDim, TNumNodes> BaseType;
    /// Node type (default is: Node<3>)
    typedef Node <3> NodeType;

    /// Geometry type (using with given NodeType)
    typedef Geometry<NodeType> GeometryType;

    /// Definition of nodes container type, redefined from GeometryType
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

    /// Vector type for local contributions to the linear system
    typedef Vector VectorType;

    /// Matrix type for local contributions to the linear system
    typedef Matrix MatrixType;

    typedef Properties PropertiesType;

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    typedef std::vector<std::size_t> EquationIdVectorType;

    typedef std::vector< Dof<double>::Pointer > DofsVectorType;

    typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;

    typedef VectorMap<IndexType, DataValueContainer> SolutionStepsElementalDataContainerType;

    /// Type for shape function values container
    typedef Kratos::Vector ShapeFunctionsType;

    /// Type for a matrix containing the shape function gradients
    typedef Kratos::Matrix ShapeFunctionDerivativesType;

    /// Type for an array of shape function gradient matrices
    typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionDerivativesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    //Constructors.

    /// Default constuctor.
    /**
     * @param NewId Index number of the new element (optional)
     */
    ComputeGradientPouliot2012(IndexType NewId = 0) :
        BaseType(NewId)
    {}

    /// Constructor using an array of nodes.
    /**
     * @param NewId Index of the new element
     * @param ThisNodes An array containing the nodes of the new element
     */
    ComputeGradientPouliot2012(IndexType NewId, const NodesArrayType& ThisNodes) :
        BaseType(NewId, ThisNodes)
    {}	

    /// Constructor using a geometry object.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     */
    ComputeGradientPouliot2012(IndexType NewId, GeometryType::Pointer pGeometry) :
        BaseType(NewId, pGeometry)
    {}

    /// Constuctor using geometry and properties.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     * @param pProperties Pointer to the element's properties
     */
    ComputeGradientPouliot2012(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) :
       BaseType(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    virtual ~ComputeGradientPouliot2012()
    {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /// Create a new element of this type
    /**
     * Returns a pointer to a new ComputeGradientPouliot2012 element, created using given input
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
                            PropertiesType::Pointer pProperties) const override
    {
        return Element::Pointer(new ComputeGradientPouliot2012(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
    }

    /// Calculate the element's local contribution to the system for the current step. It is a combination of the
    /// RHS of the Base class (weighed by a small parameter to stabilize the system) and the contribution described in
    /// the Pouliot 2012 paper
    virtual void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

    /// Provides the global indices for each one of this element's local rows
    /**
     * this determines the elemental equation ID vector for all elemental
     * DOFs
     * @param rResult A vector containing the global Id of each row
     * @param rCurrentProcessInfo the current process info object (unused)
     */
    virtual void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

    /// Returns a list of the element's Dofs
    /**
     * @param ElementalDofList the list of DOFs
     * @param rCurrentProcessInfo the current process info instance
     */
    virtual void GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo) override;

    /// Checks the input and that all required Kratos variables have been registered.
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo The ProcessInfo of the ModelPart that contains this element.
     * @return 0 if no errors were found.
     */
    virtual int Check(const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "ComputeGradientPouliot2012 #" << this->Id();
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ComputeGradientPouliot2012" << TDim << "D";
    }

//        /// Print object's data.
//        virtual void PrintData(std::ostream& rOStream) const;

    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
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

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    ///@}
    ///@name Private Operators
    ///@{

    /// Computes local contributions to the mass matrix
    /**
     * Provides the local contributions to the mass matrix, which is defined here
     * as the matrix associated to velocity derivatives. Note that the mass
     * matrix implemented here is lumped.
     * @param rMassMatrix Will be filled with the elemental mass matrix
     * @param rCurrentProcessInfo the current process info instance
     */

    virtual void AddPouliot2012LHS(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);

    virtual void AddPouliot2012StabilizationLHS(const double epsilon, MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);

    virtual void AddFEMLaplacianStabilizationLHS(const double epsilon, MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);

    void AssembleEdgeLHSContribution(const unsigned int edge[2], const array_1d<double, 3>& edge_normalized_vector, MatrixType& rLeftHandSideMatrix);

    virtual void AddPouliot2012RHS(VectorType& F, ProcessInfo& rCurrentProcessInfo);

    virtual void AddStabilizationRHSContribution(VectorType& F,
                                                 const array_1d<double, TNumNodes>& rShapeFunc,
                                                 const double Weight);

    void CalculateStabilizationRHS(const double epsilon, VectorType& F, ProcessInfo& rCurrentProcessInfo);

    void AssembleEdgeRHSContributionX(const unsigned int edge[2], const double h_edge_inv, const array_1d<double, 3>& edge_normalized_vector, VectorType& F);

    void AssembleEdgeRHSContributionY(const unsigned int edge[2], const double h_edge_inv, const array_1d<double, 3>& edge_normalized_vector, VectorType& F);

    void AssembleEdgeRHSContributionZ(const unsigned int edge[2], const double h_edge_inv, const array_1d<double, 3>& edge_normalized_vector, VectorType& F);


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    ComputeGradientPouliot2012 & operator=(ComputeGradientPouliot2012 const& rOther);

    /// Copy constructor.
    ComputeGradientPouliot2012(ComputeGradientPouliot2012 const& rOther);

    ///@}

}; // Class ComputeGradientPouliot2012

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< unsigned int TDim >
inline std::istream& operator >>(std::istream& rIStream,
                                 ComputeGradientPouliot2012<TDim>& rThis)
{
    return rIStream;
}

/// output stream function
template< unsigned int TDim >
inline std::ostream& operator <<(std::ostream& rOStream,
                                 const ComputeGradientPouliot2012<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} // Fluid Dynamics Application group

} // namespace Kratos.

#endif // KRATOS_COMPUTE_GRADIENT_POULIOT_2012_H_INCLUDED  defined
