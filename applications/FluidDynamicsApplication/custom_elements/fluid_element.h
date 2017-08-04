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
    FluidElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

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
                            PropertiesType::Pointer pProperties) const;


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


    virtual void GetValueOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
                                             std::vector<array_1d<double, 3 > >& rValues,
                                             const ProcessInfo& rCurrentProcessInfo);


    virtual void GetValueOnIntegrationPoints(const Variable<double>& rVariable,
                                             std::vector<double>& rValues,
                                             const ProcessInfo& rCurrentProcessInfo);


    virtual void GetValueOnIntegrationPoints(const Variable<array_1d<double, 6 > >& rVariable,
                                             std::vector<array_1d<double, 6 > >& rValues,
                                             const ProcessInfo& rCurrentProcessInfo);


    virtual void GetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
                                             std::vector<Vector>& rValues,
                                             const ProcessInfo& rCurrentProcessInfo);


    virtual void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                             std::vector<Matrix>& rValues,
                                             const ProcessInfo& rCurrentProcessInfo);


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
    

    /**
     * @brief EvaluateInPoint Interpolate nodal data inside the element.
     * Evaluate a nodal variable in the point where the form functions take the
     * values given by rShapeFunc and write the result to rResult.
     * This is an auxiliary function used to compute values in integration points.
     * @param rResult The variable where the value will be added to
     * @param rVar The nodal variable to be read
     * @param rShapeFunc The values of the form functions in the point
     */
    template< class TVariableType >
    void EvaluateInPoint(TVariableType& rResult,
                         const Kratos::Variable<TVariableType>& Var,
                         const ShapeFunctionsType& rShapeFunc)
    {
        GeometryType& rGeom = this->GetGeometry();
        const unsigned int NumNodes = rGeom.PointsNumber();

        rResult = rShapeFunc[0] * rGeom[0].FastGetSolutionStepValue(Var);

        for(unsigned int i = 1; i < NumNodes; i++)
        {
            rResult += rShapeFunc[i] * rGeom[i].FastGetSolutionStepValue(Var);
        }
    }


    /**
     * @brief EvaluateInPoint Interpolate nodal data inside the element.
     * Evaluate a nodal variable in the point where the form functions take the
     * values given by rShapeFunc and write the result to rResult.
     * This is an auxiliary function used to compute values in integration points.
     * @param rResult The variable where the value will be added to
     * @param rVar The nodal variable to be read
     * @param rShapeFunc The values of the form functions in the point
     * @param Step Number of time steps back
     */
    template< class TVariableType >
    void EvaluateInPoint(TVariableType& rResult,
                         const Kratos::Variable<TVariableType>& Var,
                         const ShapeFunctionsType& rShapeFunc,
                         const IndexType Step)
    {
        GeometryType& rGeom = this->GetGeometry();
        const unsigned int NumNodes = rGeom.PointsNumber();

        rResult = rShapeFunc[0] * rGeom[0].FastGetSolutionStepValue(Var,Step);

        for(unsigned int i = 1; i < NumNodes; i++)
        {
            rResult += rShapeFunc[i] * rGeom[i].FastGetSolutionStepValue(Var,Step);
        }
    }


    /// Characteristic element size h to be used in stabilization parameters.
    virtual double ElementSize();

    virtual void CalculateStaticTau(double Density,
                                    double KinematicVisc,
                                    const array_1d<double,3> &Velocity,
                                    double ElemSize,
                                    const ProcessInfo& rProcessInfo,
                                    double &TauOne,
                                    double &TauTwo);


    virtual void ASGSMomentumResidual(double GaussIndex,
                                      const ShapeFunctionsType &rN,
                                      const ShapeFunctionDerivativesType &rDN_DX,
                                      array_1d<double,3>& rMomentumRes);


    virtual void ASGSMassResidual(double GaussIndex,
                                  const ShapeFunctionsType &rN,
                                  const ShapeFunctionDerivativesType &rDN_DX,
                                  double& rMomentumRes);


    virtual void OSSMomentumResidual(double GaussIndex,
                                     const ShapeFunctionsType &rN,
                                     const ShapeFunctionDerivativesType &rDN_DX,
                                     array_1d<double,3>& rMomentumRes);

    virtual void OSSMassResidual(double GaussIndex,
                                 const ShapeFunctionsType& rN,
                                 const ShapeFunctionDerivativesType& rDN_DX,
                                 double& rMassRes);


    virtual void MomentumProjTerm(double GaussIndex,
                                  const ShapeFunctionsType &rN,
                                  const ShapeFunctionDerivativesType &rDN_DX,
                                  array_1d<double,3>& rMomentumRHS);


    virtual void MassProjTerm(double GaussIndex,
                              const ShapeFunctionsType &rN,
                              const ShapeFunctionDerivativesType &rDN_DX,
                              double& rMassRHS);


    /**
     * @brief EffectiveViscosity Evaluate the total kinematic viscosity at a given integration point.
     * This function is used to implement Smagorinsky type LES or non-Newtonian dynamics in derived classes.
     * @param rN Shape function values at integration point
     * @param rDN_DX Shape function derivatives at integration point
     * @param ElemSize Characteristic length representing the element (for Smagorinsky, this is the filter width)
     * @param rCurrentProcessInfo
     * @return Kinematic viscosity at the integration point.
     */
    virtual double EffectiveViscosity(const ShapeFunctionsType &rN,
                                      const ShapeFunctionDerivativesType &rDN_DX,
                                      double ElemSize,
                                      const ProcessInfo &rCurrentProcessInfo);


    void ResolvedConvectiveVelocity(array_1d<double,3>& rConvVel,
                                    const ShapeFunctionsType& rN);


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
        unsigned int GaussIndex,
        double GaussWeight,
        const ShapeFunctionsType& rN,
        const ShapeFunctionDerivativesType& rDN_DX,
        const ProcessInfo& rProcessInfo,
        MatrixType& rLHS,
        VectorType& rRHS) = 0;


    virtual void AddMassTerms(
        const TElementData& rData,
        double GaussWeight,
        const ShapeFunctionsType& rN,
        MatrixType& rMassMatrix) = 0;


    virtual void AddMassStabilization(
        const TElementData& rData,
        unsigned int GaussIndex,
        double GaussWeight,
        const ShapeFunctionsType& rN,
        const ShapeFunctionDerivativesType& rDN_DX,
        const ProcessInfo& rProcessInfo,
        MatrixType& rMassMatrix) = 0;


    virtual void SubscaleVelocity(unsigned int GaussIndex,
                                  const ShapeFunctionsType& rN,
                                  const ShapeFunctionDerivativesType& rDN_DX,
                                  const ProcessInfo& rProcessInfo,
                                  array_1d<double,3>& rVelocitySubscale);

    virtual void SubscalePressure(unsigned int GaussIndex,
                                  const ShapeFunctionsType& rN,
                                  const ShapeFunctionDerivativesType& rDN_DX,
                                  const ProcessInfo& rProcessInfo,
                                  double &rPressureSubscale);

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
