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

#ifndef KRATOS_D_VMS_H
#define KRATOS_D_VMS_H

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

// Forward decalration of auxiliary class
namespace Internals {
template <class TElementData, bool TDataKnowsAboutTimeIntegration>
class SpecializedAddTimeIntegratedSystem;
}

template< class TElementData >
class DVMS : public FluidElement<TElementData>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DVMS
    KRATOS_CLASS_POINTER_DEFINITION(DVMS);

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
    DVMS(IndexType NewId = 0);

    /// Constructor using an array of nodes.
    /**
     * @param NewId Index of the new element
     * @param ThisNodes An array containing the nodes of the new element
     */
    DVMS(IndexType NewId, const NodesArrayType& ThisNodes);

    /// Constructor using a geometry object.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     */
    DVMS(IndexType NewId, GeometryType::Pointer pGeometry);

    /// Constuctor using geometry and properties.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     * @param pProperties Pointer to the element's properties
     */
    DVMS(IndexType NewId, GeometryType::Pointer pGeometry, Properties::Pointer pProperties);

    /// Destructor.
    ~DVMS() override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    /// Create a new element of this type
    /**
     * Returns a pointer to a new DVMS element, created using given input
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

    /// Set up the element.
    /** Allocate the subscale velocity containers and let base class initialize the constitutive law */
    void Initialize() override;

    /// Update the values of tracked small scale quantities.
    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

    /// Predict the value of the small scale velocity for the current iteration.
    void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;

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
    void AddBoundaryIntegral(
        TElementData& rData,
        const Vector& rUnitNormal,
        MatrixType& rLHS,
        VectorType& rRHS) override;

    // Implementation details of DVMS /////////////////////////////////////////

    void AddMassStabilization(
        TElementData& rData,
        MatrixType& rMassMatrix);

    void AddViscousTerm(
        const TElementData& rData,
        boost::numeric::ublas::bounded_matrix<double,LocalSize,LocalSize>& rLHS,
        VectorType& rRHS);

    void CalculateProjections(const ProcessInfo &rCurrentProcessInfo);

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
    
    constexpr static double mTauC1 = 8.0;
    constexpr static double mTauC2 = 2.0;
    constexpr static double mSubscalePredictionVelocityTolerance = 1e-14;
    constexpr static double mSubscalePredictionResidualTolerance = 1e-14;
    constexpr static unsigned int mSubscalePredictionMaxIterations = 10;

    ///@}
    ///@name Member Variables
    ///@{

    // Velocity subscale history, stored at integration points
    std::vector< array_1d<double,Dim> > mPredictedSubscaleVelocity;
    std::vector< array_1d<double,Dim> > mOldSubscaleVelocity;

    ///@}
    ///@name Friends
    ///@{

    friend class Internals::SpecializedAddTimeIntegratedSystem<TElementData, TElementData::ElementManagesTimeIntegration>;

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

    void CalculateTau(
        const TElementData& rData,
        const array_1d<double,3> &Velocity,
        double &TauOne,
        double &TauTwo,
        double &TauP) const;

    void MomentumProjTerm(
        const TElementData& rData,
        const array_1d<double,3>& rConvectionVelocity,
        array_1d<double,3>& rMomentumRHS) const;

    void MassProjTerm(
        const TElementData& rData,
        double& rMassRHS) const;
    
    void SubscaleVelocity(
        const TElementData& rData,
        const unsigned int GaussPointIndex,
        array_1d<double,3>& rVelocitySubscale);

    void SubscalePressure(
        const TElementData& rData,
        const unsigned int GaussPointIndex,
        double &rPressureSubscale);

    void UpdateSubscaleVelocityPrediction(
        const TElementData& rData,
        const unsigned int GaussPointIndex);

    void ASGSMomentumResidual(
        const TElementData& rData,
        const array_1d<double,3> &rConvectionVelocity,
        array_1d<double,3>& rResidual) const;

    void ASGSMassResidual(
        const TElementData& rData,
        double& rResidual) const;

    void OSSMomentumResidual(
        const TElementData& rData,
        const array_1d<double,3> &rConvectionVelocity,
        array_1d<double,3>& rResidual) const;

    void OSSMassResidual(
        const TElementData& rData,
        double& rResidual) const;

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
    DVMS& operator=(DVMS const& rOther);

    /// Copy constructor.
    DVMS(DVMS const& rOther);

    ///@}


}; // Class DVMS

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< class TElementData >
inline std::istream& operator >>(std::istream& rIStream,
                                 DVMS<TElementData>& rThis)
{
    return rIStream;
}

/// output stream function
template< class TElementData >
inline std::ostream& operator <<(std::ostream& rOStream,
                                 const DVMS<TElementData>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} // Fluid Dynamics Application group

namespace Internals {

template <class TElementData, bool TDataKnowsAboutTimeIntegration>
class SpecializedAddTimeIntegratedSystem {
   public:
    static void AddSystem(DVMS<TElementData>* pElement,
        TElementData& rData, Matrix& rLHS, Vector& rRHS);
};

template <class TElementData>
class SpecializedAddTimeIntegratedSystem<TElementData, true> {
   public:
    static void AddSystem(DVMS<TElementData>* pElement,
        TElementData& rData, Matrix& rLHS, Vector& rRHS);
};

template <class TElementData>
class SpecializedAddTimeIntegratedSystem<TElementData, false> {
   public:
    static void AddSystem(DVMS<TElementData>* pElement,
        TElementData& rData, Matrix& rLHS, Vector& rRHS);
};

} // namespace Internals

} // namespace Kratos.

#endif // KRATOS_D_VMS_H
