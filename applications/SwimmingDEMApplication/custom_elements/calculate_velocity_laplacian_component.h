//
//   Project Name:        Kratos
//   Last Modified by:    $Author: gcasas $
//   Date:                $Date: 2016-03-12
//

#if !defined(KRATOS_COMPUTE_VELOCITY_LAPLACIAN_COMPONENT_SIMPLEX_ELEMENT_H_INCLUDED )
#define  KRATOS_COMPUTE_VELOCITY_LAPLACIAN_COMPONENT_SIMPLEX_ELEMENT_H_INCLUDED

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
#include "calculate_mat_deriv_simplex_element.h"

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
class KRATOS_API(SWIMMING_DEM_APPLICATION) ComputeVelocityLaplacianComponentSimplex : public ComputeMaterialDerivativeSimplex<TDim, TNumNodes>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ComputeVelocityLaplacianComponentSimplex
    KRATOS_CLASS_POINTER_DEFINITION(ComputeVelocityLaplacianComponentSimplex);

    typedef ComputeMaterialDerivativeSimplex<TDim, TNumNodes> BaseType;

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
    ComputeVelocityLaplacianComponentSimplex(IndexType NewId = 0) :
        BaseType(NewId), mCurrentComponent('X')
    {}

    /// Constructor using an array of nodes.
    /**
     * @param NewId Index of the new element
     * @param ThisNodes An array containing the nodes of the new element
     */
    ComputeVelocityLaplacianComponentSimplex(IndexType NewId, const NodesArrayType& ThisNodes) :
        BaseType(NewId, ThisNodes), mCurrentComponent('X')
    {}

    /// Constructor using a geometry object.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     */
    ComputeVelocityLaplacianComponentSimplex(IndexType NewId, GeometryType::Pointer pGeometry) :
        BaseType(NewId, pGeometry), mCurrentComponent('X')
    {}

    /// Constuctor using geometry and properties.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     * @param pProperties Pointer to the element's properties
     */
    ComputeVelocityLaplacianComponentSimplex(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) :
       BaseType(NewId, pGeometry, pProperties), mCurrentComponent('X')
    {}

    /// Destructor.
    virtual ~ComputeVelocityLaplacianComponentSimplex()
    {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /// Create a new element of this type
    /**
     * Returns a pointer to a new ComputeVelocityLaplacianComponentSimplex element, created using given input
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
                            PropertiesType::Pointer pProperties) const override
    {
        return Element::Pointer(new ComputeVelocityLaplacianComponentSimplex(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
    }

    /// Calculate the element's local contribution to the system for the current step.
    virtual void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      const ProcessInfo& rCurrentProcessInfo) override;


    /// Provides the global indices for each one of this element's local rows
    /**
     * this determines the elemental equation ID vector for all elemental
     * DOFs
     * @param rResult A vector containing the global Id of each row
     * @param rCurrentProcessInfo the current process info object (unused)
     */
    virtual void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const override;

    /// Returns a list of the element's Dofs
    /**
     * @param ElementalDofList the list of DOFs
     * @param rCurrentProcessInfo the current process info instance
     */
    virtual void GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& rCurrentProcessInfo) const override;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Elemental Data
    ///@{

    /// Checks the input and that all required Kratos variables have been registered.
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo The ProcessInfo of the ModelPart that contains this element.
     * @return 0 if no errors were found.
     */
    virtual int Check(const ProcessInfo& rCurrentProcessInfo) const override;

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
        buffer << "ComputeVelocityLaplacianComponentSimplex #" << this->Id();
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ComputeVelocityLaplacianComponentSimplex" << TDim << "D";
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
    char mCurrentComponent;

    void CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLumpedMassMatrix(MatrixType& rLHSMatrix, const double Mass) override;

    void AddConsistentMassMatrixContribution(MatrixType& rLHSMatrix, const array_1d<double,TNumNodes>& rShapeFunc, const double Weight) override;

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

    void AddIntegrationPointRHSContribution(VectorType& F,
                         const array_1d<double,TNumNodes>& rShapeFunc,
                         const BoundedMatrix<double, TNumNodes, TDim>& rShapeDeriv,
                         const double Weight) override;

//    Variable<array_1d<double, 3> > GetCurrentVariable();

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
    ComputeVelocityLaplacianComponentSimplex & operator=(ComputeVelocityLaplacianComponentSimplex const& rOther);

    /// Copy constructor.
    ComputeVelocityLaplacianComponentSimplex(ComputeVelocityLaplacianComponentSimplex const& rOther);

    ///@}

}; // Class ComputeVelocityLaplacianComponentSimplex

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< unsigned int TDim >
inline std::istream& operator >>(std::istream& rIStream,
                                 ComputeVelocityLaplacianComponentSimplex<TDim>& rThis)
{
    return rIStream;
}

/// output stream function
template< unsigned int TDim >
inline std::ostream& operator <<(std::ostream& rOStream,
                                 const ComputeVelocityLaplacianComponentSimplex<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} // Fluid Dynamics Application group

} // namespace Kratos.

#endif // KRATOS_COMPUTE_VELOCITY_LAPLACIAN_COMPONENT_SIMPLEX_ELEMENT_H_INCLUDED  defined
