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
//                   Ignasi de Pouplana
//

#ifndef KRATOS_FIC_FLUID_ELEMENT_H
#define KRATOS_FIC_FLUID_ELEMENT_H

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

/// A detailed description of the formulation can be found in:
/// A FIC-based stabilized finite element formulation for turbulent flows, 2017
/// https://www.sciencedirect.com/science/article/pii/S004578251630799X

template< class TElementData >
class FIC : public FluidElement<TElementData>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FIC
    KRATOS_CLASS_POINTER_DEFINITION(FIC);

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
    FIC(IndexType NewId = 0);

    /// Constructor using an array of nodes.
    /**
     * @param NewId Index of the new element
     * @param ThisNodes An array containing the nodes of the new element
     */
    FIC(IndexType NewId, const NodesArrayType& ThisNodes);

    /// Constructor using a geometry object.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     */
    FIC(IndexType NewId, GeometryType::Pointer pGeometry);

    /// Constuctor using geometry and properties.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     * @param pProperties Pointer to the element's properties
     */
    FIC(IndexType NewId, GeometryType::Pointer pGeometry, Properties::Pointer pProperties);

    /// Destructor.
    ~FIC() override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    /// Create a new element of this type
    /**
     * Returns a pointer to a new FIC element, created using given input
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

    void CalculateOnIntegrationPoints(
        Variable<array_1d<double, 3>> const& rVariable,
        std::vector<array_1d<double, 3>>& rValues,
        ProcessInfo const& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(
        Variable<double> const& rVariable,
        std::vector<double>& rValues,
        ProcessInfo const& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(
        Variable<array_1d<double, 6>> const& rVariable,
        std::vector<array_1d<double, 6>>& rValues,
        ProcessInfo const& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(
        Variable<Vector> const& rVariable,
        std::vector<Vector>& rValues,
        ProcessInfo const& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(
        Variable<Matrix> const& rVariable,
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
    void AddBoundaryIntegral(
        TElementData& rData,
        const Vector& rUnitNormal,
        MatrixType& rLHS,
        VectorType& rRHS) override;

    // Implementation details of FIC ////////////////////////////////////////

    void AddMassStabilization(
        TElementData& rData,
        MatrixType& rMassMatrix);

    void AddViscousTerm(
        const TElementData& rData,
        BoundedMatrix<double,LocalSize,LocalSize>& rLHS,
        VectorType& rRHS);

    virtual void CalculateTau(
        const TElementData& rData,
        const array_1d<double,3> &Velocity,
        double &TauIncompr,
        double &TauMomentum,
        array_1d<double,3> &TauGrad) const;

    virtual void CalculateTauGrad(
        const TElementData& rData,
        array_1d<double,3> &TauGrad) const;

    virtual void AlgebraicMomentumResidual(
        const TElementData& rData,
        const Vector& rConvection,
        array_1d<double,3>& rResidual) const;

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
    FIC& operator=(FIC const& rOther);

    /// Copy constructor.
    FIC(FIC const& rOther);

    ///@}


}; // Class FIC

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< class TElementData >
inline std::istream& operator >>(std::istream& rIStream,
                                 FIC<TElementData>& rThis)
{
    return rIStream;
}

/// output stream function
template< class TElementData >
inline std::ostream& operator <<(std::ostream& rOStream,
                                 const FIC<TElementData>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


namespace Internals {

template <class TElementData, bool TDataKnowsAboutTimeIntegration>
class FICSpecializedAddTimeIntegratedSystem {
   public:
    static void AddSystem(FIC<TElementData>* pElement,
        TElementData& rData, Kratos::Matrix& rLHS, Vector& rRHS);
};

template <class TElementData>
class FICSpecializedAddTimeIntegratedSystem<TElementData, true> {
   public:
    static void AddSystem(FIC<TElementData>* pElement,
        TElementData& rData, Kratos::Matrix& rLHS, Vector& rRHS);
};

template <class TElementData>
class FICSpecializedAddTimeIntegratedSystem<TElementData, false> {
   public:
    static void AddSystem(FIC<TElementData>* pElement,
        TElementData& rData, Kratos::Matrix& rLHS, Vector& rRHS);
};

///@} // Fluid Dynamics Application group

} // namespace Internals

} // namespace Kratos.

#endif // KRATOS_FIC_FLUID_ELEMENT_H
