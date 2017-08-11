//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#ifndef KRATOS_FLUID_ELEMENT_H
#define KRATOS_FLUID_ELEMENT_H

#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "geometries/geometry.h"

#include "includes/cfd_variables.h"
#include "fluid_dynamics_application_variables.h"

#include "custom_elements/integration_point_data.h"

namespace Kratos
{

///@addtogroup FluidDynamicsApplication
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

template< class TElementData >
class FluidElement : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FluidElement
    KRATOS_CLASS_POINTER_DEFINITION(FluidElement);

    /// Node type (default is: Node<3>)
    typedef Node<3> NodeType;

    /// Geometry type (using with given NodeType)
    typedef Geometry<NodeType> GeometryType;

    /// Definition of nodes container type, redefined from GeometryType
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

    /// Vector type for local contributions to the linear system
    typedef Vector VectorType;

    /// Matrix type for local contributions to the linear system
    typedef Matrix MatrixType;

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
    FluidElement(IndexType NewId = 0);

    /// Constructor using an array of nodes.
    /**
     * @param NewId Index of the new element
     * @param ThisNodes An array containing the nodes of the new element
     */
    FluidElement(IndexType NewId, const NodesArrayType& ThisNodes);

    /// Constructor using a geometry object.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     */
    FluidElement(IndexType NewId, GeometryType::Pointer pGeometry);

    /// Constuctor using geometry and properties.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     * @param pProperties Pointer to the element's properties
     */
    FluidElement(IndexType NewId, GeometryType::Pointer pGeometry, Properties::Pointer pProperties);

    /// Destructor.
    virtual ~FluidElement();

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    /// Create a new element of this type
    /**
     * Returns a pointer to a new FluidElement element, created using given input
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId,
                            NodesArrayType const& ThisNodes,
                            Properties::Pointer pProperties) const override;

    /// Create a new element of this type using given geometry
    /**
     * Returns a pointer to a new FluidElement element, created using given input
     * @param NewId: the ID of the new element
     * @param pGeom: a pointer to the geomerty to be used to create the element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId,
                            GeometryType::Pointer pGeom,
                            Properties::Pointer pProperties) const override;

    /**
     * @brief CalculateLocalSystem Return empty matrices and vectors of appropriate size.
     * This element does not have a local contribution in terms of displacements, but the scheme may
     * require a proper-sized matrix, even if it is empty.
     * @param rLeftHandSideMatrix Local finite element system matrix (output)
     * @param rRightHandSideVector Local finite element residual vector (output)
     * @param rCurrentProcessInfo Current ProcessInfo values (input)
     */
    virtual void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      ProcessInfo& rCurrentProcessInfo);

    /**
     * @brief CalculateLeftHandSide Return an empty matrix of appropriate size.
     * This element does not have a local contribution in terms of displacements, but the scheme may
     * require a proper-sized matrix, even if it is empty.
     * @param rLeftHandSideMatrix Local finite element system matrix (output)
     * @param rCurrentProcessInfo Current ProcessInfo values (input)
     */
    virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                       ProcessInfo& rCurrentProcessInfo);

    /**
     * @brief CalculateRightHandSide Return an empty matrix of appropriate size.
     * This element does not have a local contribution in terms of displacements, but the scheme may
     * require a proper-sized matrix, even if it is empty.
     * @param rRightHandSideVector Local finite element residual vector (output)
     * @param rCurrentProcessInfo Current ProcessInfo values (input)
     */
    virtual void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                        ProcessInfo& rCurrentProcessInfo);

    /**
     * @brief CalculateLocalVelocityContribution Calculate the local contribution in terms of velocity and pressure.
     * @param rDampMatrix Local finite element system matrix (output)
     * @param rRightHandSideVector Local finite element residual vector (output)
     * @param rCurrentProcessInfo Current ProcessInfo values (input)
     */
    virtual void CalculateLocalVelocityContribution(MatrixType &rDampMatrix,
                                                    VectorType &rRightHandSideVector,
                                                    ProcessInfo &rCurrentProcessInfo);

    /**
     * @brief MassMatrix Calculate the local mass matrix.
     * @param rMassMatrix Local mass matrix (output)
     * @param rCurrentProcessInfo Current ProcessInfo values (input)
     */
    virtual void CalculateMassMatrix(MatrixType &rMassMatrix,
                                     ProcessInfo &rCurrentProcessInfo);


    virtual void Calculate(const Variable<double>& rVariable,
                           double& rOutput,
                           const ProcessInfo& rCurrentProcessInfo);


    virtual void Calculate(const Variable<array_1d<double, 3 > >& rVariable,
                           array_1d<double, 3 > & rOutput,
                           const ProcessInfo& rCurrentProcessInfo);

    /**
     * @brief EquationIdVector Returns the global system rows corresponding to each local row.
     * @param rResult rResult[i] is the global index of local row i (output)
     * @param rCurrentProcessInfo Current ProcessInfo values (input)
     */
    virtual void EquationIdVector(EquationIdVectorType& rResult,
                                  ProcessInfo& rCurrentProcessInfo);

    /**
     * @brief GetDofList Returns a list of the element's Dofs.
     * @param rElementalDofList List of DOFs. (output)
     * @param rCurrentProcessInfo Current ProcessInfo instance. (input)
     */
    virtual void GetDofList(DofsVectorType& rElementalDofList,
                            ProcessInfo& rCurrentProcessInfo);


    /**
     * @brief GetFirstDerivativesVector Returns VELOCITY_X, VELOCITY_Y, (VELOCITY_Z,) PRESSURE for each node.
     * @param Values Vector of nodal unknowns
     * @param Step Get result from 'Step' steps back, 0 is current step. (Must be smaller than buffer size)
     */
    virtual void GetFirstDerivativesVector(Vector& Values, int Step = 0);



    /**
     * @brief Returns ACCELERATION_X, ACCELERATION_Y, (ACCELERATION_Z,) 0 for each node
     * @param Values Vector of nodal second derivatives
     * @param Step Get result from 'Step' steps back, 0 is current step. (Must be smaller than buffer size)
     */
    virtual void GetSecondDerivativesVector(Vector& Values, int Step = 0);


    /**
     * @brief GetIntegrationMethod Return the integration order to be used.
     * @return Gauss Order
     */
    virtual GeometryData::IntegrationMethod GetIntegrationMethod() const;


    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{

    virtual int Check(const ProcessInfo &rCurrentProcessInfo);

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const;


    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const;


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


    /// Determine integration point weights and shape funcition derivatives from the element's geometry.
    virtual void CalculateGeometryData(Vector& rGaussWeights,
                                       Matrix& rNContainer,
                                       ShapeFunctionDerivativesArrayType& rDN_DX);

    void EvaluateInPoint(
        double& rResult,
        const typename TElementData::ScalarDataType& rNodalValues,
        const ShapeFunctionsType& rShapeFunc);

    void EvaluateInPoint(
        array_1d<double,3>& rResult,
        const typename TElementData::VectorDataType& rNodalValues,
        const ShapeFunctionsType& rShapeFunc);


    /// Characteristic element size h to be used in stabilization parameters.
    virtual double ElementSize();

    virtual void CalculateStaticTau(double Density,
                                    double KinematicVisc,
                                    const array_1d<double,3> &Velocity,
                                    double ElemSize,
                                    const ProcessInfo& rProcessInfo,
                                    double &TauOne,
                                    double &TauTwo);

    /**
     * @brief EffectiveViscosity Evaluate the total kinematic viscosity at a given integration point.
     * This function is used to implement Smagorinsky type LES or non-Newtonian dynamics in derived classes.
     * @param rData TElementData instance with information about nodal values
     * @param rN Shape function values at integration point
     * @param rDN_DX Shape function derivatives at integration point
     * @param ElemSize Characteristic length representing the element (for Smagorinsky, this is the filter width)
     * @param rCurrentProcessInfo
     * @return Kinematic viscosity at the integration point.
     */
    virtual double EffectiveViscosity(
        const TElementData& rData,
        const IntegrationPointData<TElementData>& rIPData,
        double ElemSize,
        const ProcessInfo &rCurrentProcessInfo);


    void ResolvedConvectiveVelocity(
        const TElementData& rData,
        const ShapeFunctionsType &rN,
        array_1d<double,3> &rConvVel);


    void FullConvectiveVelocity(array_1d<double,3>& rConvVel,
                                const ShapeFunctionsType& rN,
                                const array_1d<double,3>& rSubscaleVel);


    /**
     * @brief Write the convective operator evaluated at this point (for each nodal funciton) to an array
     * Evaluate the convective operator for each node's shape function at an arbitrary point
     * @param rResult: Output vector
     * @param rVelocity: Velocity evaluated at the integration point
     * @param rShapeDeriv: Derivatives of shape functions evaluated at the integration point
     */
    void ConvectionOperator(Vector& rResult,
                            const array_1d<double,3>& rConvVel,
                            const ShapeFunctionDerivativesType& DN_DX);


    virtual void AddSystemTerms(
        const TElementData& rData,
        const IntegrationPointData<TElementData>& rIPData,
        const ProcessInfo& rProcessInfo,
        MatrixType& rLHS,
        VectorType& rRHS) = 0;


    virtual void AddMassTerms(
        const TElementData& rData,
        const IntegrationPointData<TElementData>& rIPData,
        MatrixType& rMassMatrix) = 0;


    virtual void AddMassStabilization(
        const TElementData& rData,
        const IntegrationPointData<TElementData>& rIPData,
        const ProcessInfo& rProcessInfo,
        MatrixType& rMassMatrix) = 0;

    virtual void CalculateProjections() = 0;

    void IntegrationPointVorticity(const ShapeFunctionDerivativesType& rDN_DX,
                                   array_1d<double,3> &rVorticity) const;


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

    virtual void save(Serializer& rSerializer) const;

    virtual void load(Serializer& rSerializer);

    ///@}
    ///@name Private Operators
    ///@{


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
    FluidElement& operator=(FluidElement const& rOther);

    /// Copy constructor.
    FluidElement(FluidElement const& rOther);

    ///@}


}; // Class FluidElement

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< class TElementData >
inline std::istream& operator >>(std::istream& rIStream,
                                 FluidElement<TElementData>& rThis)
{
    return rIStream;
}

/// output stream function
template< class TElementData >
inline std::ostream& operator <<(std::ostream& rOStream,
                                 const FluidElement<TElementData>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} // Fluid Dynamics Application group

} // namespace Kratos.

#endif // KRATOS_FLUID_ELEMENT_H
