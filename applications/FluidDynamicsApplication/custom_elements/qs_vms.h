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

#ifndef KRATOS_QS_VMS_H
#define KRATOS_QS_VMS_H

#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "geometries/geometry.h"

#include "includes/cfd_variables.h"
#include "custom_elements/fluid_element.h"
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
class QSVMS : public FluidElement<TElementData>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of QSVMS
    KRATOS_CLASS_POINTER_DEFINITION(QSVMS);

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

    constexpr static unsigned int Dim = FluidElement<TElementData>::Dim;
    constexpr static unsigned int NumNodes = FluidElement<TElementData>::NumNodes;
    constexpr static unsigned int BlockSize = FluidElement<TElementData>::BlockSize;
    constexpr static unsigned int LocalSize = FluidElement<TElementData>::LocalSize;
    constexpr static unsigned int StrainSize = FluidElement<TElementData>::StrainSize;

    ///@}
    ///@name Life Cycle
    ///@{

    //Constructors.

    /// Default constuctor.
    /**
     * @param NewId Index number of the new element (optional)
     */
    QSVMS(IndexType NewId = 0);

    /// Constructor using an array of nodes.
    /**
     * @param NewId Index of the new element
     * @param ThisNodes An array containing the nodes of the new element
     */
    QSVMS(IndexType NewId, const NodesArrayType& ThisNodes);

    /// Constructor using a geometry object.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     */
    QSVMS(IndexType NewId, GeometryType::Pointer pGeometry);

    /// Constuctor using geometry and properties.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     * @param pProperties Pointer to the element's properties
     */
    QSVMS(IndexType NewId, GeometryType::Pointer pGeometry, Properties::Pointer pProperties);

    /// Destructor.
    ~QSVMS() override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    /// Create a new element of this type
    /**
     * Returns a pointer to a new QSVMS element, created using given input
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId,
                            NodesArrayType const& ThisNodes,
                            Properties::Pointer pProperties) const override;

    /// Create a new element of this type using given geometry
    /**
     * Returns a pointer to a new FluidElement element, created using given input
     * @param NewId the ID of the new element
     * @param pGeom a pointer to the geomerty to be used to create the element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId,
                            GeometryType::Pointer pGeom,
                            Properties::Pointer pProperties) const override;


    void Calculate(
        const Variable<double>& rVariable,
        double& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;


    void Calculate(
        const Variable<array_1d<double, 3 > >& rVariable,
        array_1d<double, 3 > & rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;


    void Calculate(
        const Variable<Vector >& rVariable,
        Vector& Output,
        const ProcessInfo& rCurrentProcessInfo) override;

    void Calculate(
        const Variable<Matrix >& rVariable,
        Matrix& Output,
        const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Access
    ///@{

    void GetValueOnIntegrationPoints(Variable<array_1d<double, 3>> const& rVariable,
                                     std::vector<array_1d<double, 3>>& rValues,
                                     ProcessInfo const& rCurrentProcessInfo) override;

    void GetValueOnIntegrationPoints(Variable<double> const& rVariable,
                                     std::vector<double>& rValues,
                                     ProcessInfo const& rCurrentProcessInfo) override;

    void GetValueOnIntegrationPoints(Variable<array_1d<double, 6>> const& rVariable,
                                     std::vector<array_1d<double, 6>>& rValues,
                                     ProcessInfo const& rCurrentProcessInfo) override;

    void GetValueOnIntegrationPoints(Variable<Vector> const& rVariable,
                                     std::vector<Vector>& rValues,
                                     ProcessInfo const& rCurrentProcessInfo) override;

    void GetValueOnIntegrationPoints(Variable<Matrix> const& rVariable,
                                     std::vector<Matrix>& rValues,
                                     ProcessInfo const& rCurrentProcessInfo) override;

    ///@}
    ///@name Inquiry
    ///@{

    int Check(const ProcessInfo &rCurrentProcessInfo) override;

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


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    // Protected interface of FluidElement ////////////////////////////////////

    void AddTimeIntegratedSystem(
        TElementData& rData,
        MatrixType& rLHS,
        VectorType& rRHS) override;

    void AddTimeIntegratedLHS(
        TElementData& rData,
        MatrixType& rLHS) override;

    void AddTimeIntegratedRHS(
        TElementData& rData,
        VectorType& rRHS) override;

    void AddVelocitySystem(
        TElementData& rData,
        MatrixType& rLocalLHS,
        VectorType& rLocalRHS) override;

    void AddMassLHS(
        TElementData& rData,
        MatrixType& rMassMatrix) override;

    // This function integrates the traction over a cut. It is only required to implement embedded formulations
    void AddBoundaryTraction(
        TElementData& rData,
        const Vector& rUnitNormal,
        MatrixType& rLHS,
        VectorType& rRHS) override;

    // Implementation details of QSVMS ////////////////////////////////////////

    void AddMassStabilization(
        TElementData& rData,
        MatrixType& rMassMatrix);

    void AddViscousTerm(
        const TElementData& rData,
        BoundedMatrix<double,LocalSize,LocalSize>& rLHS,
        VectorType& rRHS);

    /**
     * @brief EffectiveViscosity Evaluate the total kinematic viscosity at a given integration point.
     * This function is used to implement Smagorinsky type LES or non-Newtonian dynamics in derived classes.
     * @param rData TElementData instance with information about nodal values
     * @param ElemSize Characteristic length representing the element (for Smagorinsky, this is the filter width)
     * @return Kinematic viscosity at the integration point.
     */
    KRATOS_DEPRECATED virtual double EffectiveViscosity(
        TElementData& rData,
        double ElementSize);

    virtual void CalculateTau(
        const TElementData& rData,
        const array_1d<double,3> &Velocity,
        double &TauOne,
        double &TauTwo) const;

    virtual void CalculateProjections(const ProcessInfo &rCurrentProcessInfo);

    virtual void MomentumProjTerm(
        const TElementData& rData,
        const array_1d<double,3>& rConvectionVelocity,
        array_1d<double,3>& rMomentumRHS) const;

    virtual void MassProjTerm(
        const TElementData& rData,
        double& rMassRHS) const;

    virtual void SubscaleVelocity(
        const TElementData& rData,
        array_1d<double,3>& rVelocitySubscale) const;

    virtual void SubscalePressure(
        const TElementData& rData,
        double &rPressureSubscale) const;

    virtual void AlgebraicMomentumResidual(
        const TElementData& rData,
        const array_1d<double,3> &rConvectionVelocity,
        array_1d<double,3>& rResidual) const;

    virtual void AlgebraicMassResidual(
        const TElementData& rData,
        double& rMomentumRes) const;

    virtual void OrthogonalMomentumResidual(
        const TElementData& rData,
        const array_1d<double,3> &rConvectionVelocity,
        array_1d<double,3>& rResidual) const;

    virtual void OrthogonalMassResidual(
        const TElementData& rData,
        double& rMassRes) const;

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
    QSVMS& operator=(QSVMS const& rOther);

    /// Copy constructor.
    QSVMS(QSVMS const& rOther);

    ///@}


}; // Class QSVMS

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< class TElementData >
inline std::istream& operator >>(std::istream& rIStream,
                                 QSVMS<TElementData>& rThis)
{
    return rIStream;
}

/// output stream function
template< class TElementData >
inline std::ostream& operator <<(std::ostream& rOStream,
                                 const QSVMS<TElementData>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} // Fluid Dynamics Application group

} // namespace Kratos.

#endif // KRATOS_QS_VMS_H
