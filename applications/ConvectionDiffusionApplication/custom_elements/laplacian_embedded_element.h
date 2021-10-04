// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Franziska Wahl
//

#ifndef KRATOS_LAPLACIAN_EMBEDDED_ELEMENT_H_INCLUDED
#define  KRATOS_LAPLACIAN_EMBEDDED_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "modified_shape_functions/modified_shape_functions.h"

// Application includes
#include "laplacian_element.h"

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

// Forward declaration of data container class
namespace LaplacianEmbeddedInternals {
    template<std::size_t TDim> class EmbeddedElementData;
} 

template<std::size_t TDim>
class LaplacianEmbeddedElement : public LaplacianElement
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of LaplacianEmbeddedElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(LaplacianEmbeddedElement);

    typedef LaplacianElement BaseType;
    typedef GeometryData::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;

    static constexpr std::size_t NumNodes = TDim + 1;

    using NodalScalarData = array_1d<double,NumNodes>;
    using EmbeddedElementData = LaplacianEmbeddedInternals::EmbeddedElementData < TDim >;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    LaplacianEmbeddedElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry);

    LaplacianEmbeddedElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~LaplacianEmbeddedElement();

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const override;

    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{


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

    /**
     * @brief Intersected element data structure initialization
     * This method sets the data structure geometry fields (shape functions, gradients) for an
     * intersected element. To do that, the modified shape functions utility is firstly created and then called
     * to perform all operations on the positive side of the element.
     */
    void InitializeGeometryData(
        EmbeddedElementData& rData);

    /**
     * @brief For an intersected element, normalize the interface normals
     * This method normalizes the interface normals for an intersected element.
     * @param rNormals interface normals container
     * @param Tolerance tolerance to avoid division by 0 when normalizing
     */
    void NormalizeInterfaceNormals(
        typename EmbeddedElementData::InterfaceNormalsType& rNormals,
        double Tolerance) const;

    /**
     * @brief Calculation of positive side for intersected elements
     * This method calculates a volume integral of the positive side
     * to the local system of an intersected element
     */
    void AddPositiveElementSide(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const EmbeddedElementData& rData);

    /**
     * @brief Calculation of the interface terms for intersected elements
     * This method calculates the interface terms on the positive side of an intersected element
     * by performing an interface integral.
     */
    void AddPositiveInterfaceTerms(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const EmbeddedElementData& rData);

    /**
     * @brief Calculation of the Nitsche boundary terms for intersected elements
     * This method calculates the Nitsche boundary terms on the positive side of an intersected element
     * by performing interface integrals for the weak imposition of a Dirichlet boundary condition.
     */
    void AddNitscheBoundaryTerms(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const EmbeddedElementData& rData);

    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{

    // Protected default constructor necessary for serialization
    LaplacianEmbeddedElement() : LaplacianElement()
    {
    }

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

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, LaplacianElement);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, LaplacianElement);
    }

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
    //LaplacianEmbeddedElement& operator=(const LaplacianEmbeddedElement& rOther);

    /// Copy constructor.
    //LaplacianEmbeddedElement(const LaplacianEmbeddedElement& rOther);

    ///@}

}; // Class LaplacianEmbeddedElement

namespace LaplacianEmbeddedInternals {

template <size_t TDim, size_t TNumNodes>
ModifiedShapeFunctions::Pointer GetContinuousShapeFunctionCalculator(
    const Element &rElement,
    const Vector &rNodalDistances);

template<std::size_t TDim>
class EmbeddedElementData
{
public:
    ///@name Type Definitions
    ///@{

    static constexpr std::size_t NumNodes = TDim + 1;

    typedef GeometryData::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;
    typedef std::vector<array_1d<double,3>> InterfaceNormalsType;
    typedef array_1d<double,NumNodes> NodalScalarData;

    ///@}
    ///@name Public Members
    ///@{

    double PenaltyCoefficient;

    NodalScalarData NodalDistances;

    Matrix PositiveSideN;
    ShapeFunctionsGradientsType PositiveSideDNDX;
    Vector PositiveSideWeights;

    Matrix PositiveInterfaceN;
    ShapeFunctionsGradientsType PositiveInterfaceDNDX;
    Vector PositiveInterfaceWeights;
    InterfaceNormalsType PositiveInterfaceUnitNormals;

    std::vector< size_t > PositiveIndices;

    size_t NumPositiveNodes;
    size_t NumNegativeNodes;

    ///@}
    ///@name Public Operations
    ///@{

    /**
     * @brief Split element data container initialization
     * This method initializes the embedded formulation data container. This implies to get the nodal distances.
     */
    void Initialize(
        const Element& rElement
    );

    /**
     * @brief Checks if the current element is intersected
     * Checks if the current element is intersected by checking the number of positive and negative distance nodes.
     * @return true if the element is intersected
     * @return false if the element is not intersected
     */
    bool IsSplit();

    ///@}
};

} //namespace LaplacianEmbeddedInternals

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
				    LaplacianEmbeddedElement& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
				    const LaplacianEmbeddedElement& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

}  // namespace Kratos.

#endif // KRATOS_LAPLACIAN_EMBEDDED_ELEMENT_H_INCLUDED  defined
