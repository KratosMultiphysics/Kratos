//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#ifndef KRATOS_DYNSS_H
#define KRATOS_DYNSS_H

#include "includes/define.h"
#include "custom_elements/dss.h"
#include "includes/serializer.h"
#include "geometries/geometry.h"

#include "stabilized_cfd_application_variables.h"
#include "../FluidDynamicsApplication/fluid_dynamics_application_variables.h"

namespace Kratos
{

///@addtogroup StabilizedCFDApplication
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

template< unsigned int TDim >
class DynSS : public DSS<TDim>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DSS
    KRATOS_CLASS_POINTER_DEFINITION(DynSS);

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

    typedef Element::PropertiesType PropertiesType;

    ///@}
    ///@name Life Cycle
    ///@{

    //Constructors.

    /// Default constuctor.
    /**
     * @param NewId Index number of the new element (optional)
     */
    DynSS(IndexType NewId = 0);

    /// Constructor using an array of nodes.
    /**
     * @param NewId Index of the new element
     * @param ThisNodes An array containing the nodes of the new element
     */
    DynSS(IndexType NewId, const NodesArrayType& ThisNodes);

    /// Constructor using a geometry object.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     */
    DynSS(IndexType NewId, GeometryType::Pointer pGeometry);

    /// Constuctor using geometry and properties.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     * @param pProperties Pointer to the element's properties
     */
    DynSS(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~DynSS();

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    /// Create a new element of this type
    /**
     * Returns a pointer to a new DSS element, created using given input
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId,
                            NodesArrayType const& ThisNodes,
                            PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId,
                           GeometryType::Pointer pGeom,
                           PropertiesType::Pointer pProperties) const override;

    /// Allocate memory for small scale values to be tracked.
    void Initialize() override;

    /// Update the values of tracked small scale quantities.
    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

    /// Predict the value of the small scale velocity for the current iteration.
    void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;


    ///@}
    ///@name Access
    ///@{


    void GetValueOnIntegrationPoints(const Variable<double>& rVariable,
                                     std::vector<double>& rValues,
                                     const ProcessInfo& rCurrentProcessInfo) override;

    void GetValueOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
                                     std::vector<array_1d<double, 3 > >& rValues,
                                     const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * TO COPY SUBSCALE VALUES AFTER SERIALIZATION
     */
    void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                  std::vector<double>& rOutput,
                  const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
                  std::vector< array_1d<double, 3 > >& rOutput,
                  const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<Vector >& rVariable,
                  std::vector< Vector >& rOutput,
                  const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<Matrix >& rVariable,
                  std::vector< Matrix >& rOutput,
                  const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * TO COPY SUBSCALE VALUES AFTER SERIALIZATION
     */
    void SetValueOnIntegrationPoints(const Variable<double>& rVariable,
                 std::vector<double>& rValues,
                 const ProcessInfo& rCurrentProcessInfo) override;

    void SetValueOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
                 std::vector<array_1d<double, 3 > > rValues,
                 const ProcessInfo& rCurrentProcessInfo) override;

    void SetValueOnIntegrationPoints(const Variable<array_1d<double, 6 > >& rVariable,
                 std::vector<array_1d<double, 6 > > rValues,
                 const ProcessInfo& rCurrentProcessInfo) override;

    void SetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
                 std::vector<Vector>& rValues,
                 const ProcessInfo& rCurrentProcessInfo) override;

    void SetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
                 std::vector<Matrix>& rValues,
                 const ProcessInfo& rCurrentProcessInfo) override;

    void SetValueOnIntegrationPoints(const Variable<ConstitutiveLaw::Pointer>& rVariable,
                 std::vector<ConstitutiveLaw::Pointer>& rValues,
                 const ProcessInfo& rCurrentProcessInfo) override;


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;


    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;


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

    ///@TODO: should not be in protected
    std::vector< array_1d<double,3> > mPredSsVel;
    std::vector< array_1d<double,3> > mOldSsVel;

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    virtual void CalculateTau(double Density,
                              double KinematicVisc,
                              const array_1d<double,3> &Velocity,
                              const ProcessInfo& rProcessInfo,
                              double ElementSize,
                              double &TauOne,
                              double &TauTwo,
                              double &TauP);


    virtual void DASGSMomentumResidual(const ShapeFunctionsType &rN,
                                      const ShapeFunctionDerivativesType &rDN_DX,
                                      const array_1d<double,3> &rConvVel,
                                      array_1d<double,3>& rMomentumRes);


    void ASGSMassResidual(double GaussIndex,
                          const ShapeFunctionsType &rN,
                          const ShapeFunctionDerivativesType &rDN_DX,
                          double& rMomentumRes) override;


    virtual void DOSSMomentumResidual(const ShapeFunctionsType &rN,
                                     const ShapeFunctionDerivativesType &rDN_DX,
                                     const array_1d<double,3> &rConvVel,
                                     array_1d<double,3>& rMomentumRes);

    void OSSMassResidual(double GaussIndex,
                         const ShapeFunctionsType& rN,
                         const ShapeFunctionDerivativesType& rDN_DX,
                         double& rMassRes) override;


    virtual void DynamicMomentumProjTerm(const ShapeFunctionsType &rN,
                                  const ShapeFunctionDerivativesType &rDN_DX,
                                  const array_1d<double,3> &rConvVel,
                                  array_1d<double,3>& rMomentumRHS);


    void MassProjTerm(double GaussIndex,
                     const ShapeFunctionsType &rN,
                     const ShapeFunctionDerivativesType &rDN_DX,
                     double& rMassRHS) override;


    virtual void FullConvectiveVelocity(array_1d<double,3>& rConvVel,
                                        const ShapeFunctionsType& rN,
                                        const array_1d<double,3>& rSubscaleVel);



    void AddSystemTerms(unsigned int GaussIndex,
                        double GaussWeight,
                        const ShapeFunctionsType& rN,
                        const ShapeFunctionDerivativesType& rDN_DX,
                        const ProcessInfo& rProcessInfo,
                        MatrixType& rLHS,
                        VectorType& rRHS) override;


    void AddMassStabilization(unsigned int GaussIndex,
                              double GaussWeight,
                              const ShapeFunctionsType& rN,
                              const ShapeFunctionDerivativesType& rDN_DX,
                              const ProcessInfo& rProcessInfo,
                              MatrixType& rMassMatrix) override;


    virtual void DenseSystemSolve(const Matrix &rA,
                                  const Vector &rB,
                                  Vector &rX) const;

    void SubscaleVelocity(unsigned int GaussIndex,
                          const ShapeFunctionsType& rN,
                          const ShapeFunctionDerivativesType& rDN_DX,
                          const ProcessInfo& rProcessInfo,
                          array_1d<double,3>& rVelocitySubscale) override;

    void SubscalePressure(unsigned int GaussIndex,
                          const ShapeFunctionsType& rN,
                          const ShapeFunctionDerivativesType& rDN_DX,
                          const ProcessInfo& rProcessInfo,
                          double &rPressureSubscale) override;



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

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

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
    DynSS& operator=(DynSS const& rOther);

    /// Copy constructor.
    DynSS(DynSS const& rOther);

    ///@}


}; // Class DynSS

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< unsigned int TDim >
inline std::istream& operator >>(std::istream& rIStream,
                                 DynSS<TDim>& rThis)
{
    return rIStream;
}

/// output stream function
template< unsigned int TDim >
inline std::ostream& operator <<(std::ostream& rOStream,
                                 const DynSS<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} // Fluid Dynamics Application group

} // namespace Kratos.

#endif // KRATOS_DYNSS_H
