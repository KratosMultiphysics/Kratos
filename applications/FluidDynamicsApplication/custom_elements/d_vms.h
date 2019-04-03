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

//#define KRATOS_D_VMS_SUBSCALE_ERROR_INSTRUMENTATION

#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "geometries/geometry.h"

#include "includes/cfd_variables.h"
#include "custom_elements/qs_vms.h"
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
class DVMS : public QSVMS<TElementData>
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

    constexpr static unsigned int Dim = QSVMS<TElementData>::Dim;
    constexpr static unsigned int NumNodes = QSVMS<TElementData>::NumNodes;
    constexpr static unsigned int BlockSize = QSVMS<TElementData>::BlockSize;
    constexpr static unsigned int LocalSize = QSVMS<TElementData>::LocalSize;
    constexpr static unsigned int StrainSize = QSVMS<TElementData>::StrainSize;

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
     * Returns a pointer to a new DVMS element, created using given input
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

    void AddVelocitySystem(
        TElementData& rData,
        MatrixType& rLocalLHS,
        VectorType& rLocalRHS) override;

    void AddMassLHS(
        TElementData& rData,
        MatrixType& rMassMatrix) override;

    // Implementation details of DVMS /////////////////////////////////////////

    void AddMassStabilization(
        TElementData& rData,
        MatrixType& rMassMatrix);

    void CalculateProjections(const ProcessInfo &rCurrentProcessInfo) override;

    void CalculateStabilizationParameters(
        const TElementData& rData,
        const array_1d<double,3> &Velocity,
        double &TauOne,
        double &TauTwo,
        double &TauP) const;

    void SubscaleVelocity(
        const TElementData& rData,
        array_1d<double,3>& rVelocitySubscale) const override;

    void SubscalePressure(
        const TElementData& rData,
        double &rPressureSubscale) const override;

    array_1d<double,3> FullConvectiveVelocity(
        const TElementData& rData) const;


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
    DenseVector< array_1d<double,Dim> > mPredictedSubscaleVelocity;
    DenseVector< array_1d<double,Dim> > mOldSubscaleVelocity;

    #ifdef KRATOS_D_VMS_SUBSCALE_ERROR_INSTRUMENTATION
    std::vector< double > mSubscaleIterationError;
    std::vector< unsigned int > mSubscaleIterationCount;
    #endif

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

    void UpdateSubscaleVelocityPrediction(
        const TElementData& rData);

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

} // namespace Kratos.

#endif // KRATOS_D_VMS_H
