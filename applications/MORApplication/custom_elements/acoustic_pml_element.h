// KRATOS
//
//  License:		 BSD License
//					 license: ../../license.txt
//
//  Main authors:
//
//


#if !defined(KRATOS_ACOUSTIC_PML_ELEMENT_H_INCLUDED )
#define  KRATOS_ACOUSTIC_PML_ELEMENT_H_INCLUDED

// System includes


// External includes

// Project includes
#include "includes/define.h"
#include "mor_application_variables.h"
#include "acoustic_element.h"
#include "includes/ublas_complex_interface.h"

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

/**
 * @class AcousticElement
 * @ingroup MORApplication
 * @brief Acoustic element
 * @details
 * @author
 * @author
 */

class AcousticPMLElement
    : public AcousticElement
{
public:
    ///@name Type Definitions
    ///@{
    ///Reference type definition for constitutive laws
    typedef ConstitutiveLaw ConstitutiveLawType;
    ///Pointer type for constitutive laws
    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;
    ///Type definition for integration methods
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    /// Type for shape function values container
    typedef Kratos::Matrix ShapeFunctionsType;

    /// Type for a matrix containing the shape function gradients
    typedef Kratos::Matrix ShapeFunctionDerivativesType;

    /// Type for an array of shape function gradient matrices
    typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionDerivativesArrayType;

    /// The base element type
    typedef AcousticElement BaseType;


    /// The definition of the index type
    typedef std::size_t IndexType;

    /// The definition of the sizetype
    typedef std::size_t SizeType;

    /// THed definition of the coordinates array
    typedef BaseType CoordinatesArrayType;


    /// Counted pointer of AcousticElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(AcousticPMLElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    AcousticPMLElement(IndexType NewId, GeometryType::Pointer pGeometry);
    AcousticPMLElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    // Copy constructor
    AcousticPMLElement(AcousticPMLElement const& rOther)
        :BaseType(rOther)
    {};

    /// Destructor.
    ~AcousticPMLElement() override;

    ///@}
    ///@name Operators
    ///@{
    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Creates a new element
     * @param NewId The Id of the new created element
     * @param pGeom The pointer to the geometry of the element
     * @param pProperties The pointer to property
     * @return The pointer to the created element
     */
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const override;

    /**
     * @brief Creates a new element
     * @param NewId The Id of the new created element
     * @param ThisNodes The array containing nodes
     * @param pProperties The pointer to property
     * @return The pointer to the created element
     */
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
        ) const override;

    /**
     * @brief It creates a new element pointer and clones the previous element data
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone (
        IndexType NewId,
        NodesArrayType const& rThisNodes
        ) const override;

    /**
     * @brief This function provides the place to perform checks on the completeness of the input.
     * @details It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo The current process info instance
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) override;

    // ///@}
    // ///@name Access
    // ///@{

    // ///@}
    // ///@name Inquiry
    // ///@{
    // ///@}
    // ///@name Input and output
    // ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Acoustic PML Element #" << Id() << "\n";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Acoustic PML Element #" << Id() << "\n";
    }

    void EquationIdVector(
            EquationIdVectorType& rResult,
            ProcessInfo& rCurrentProcessInfo) override;

    void GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo) override;

    void Initialize() override;


protected:
    
    // ///@}
    // ///@name Protected Operations
    // ///@{

    void CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo) override;

    void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

    void ComplexJacobian(ComplexMatrix& rResult, IndexType IntegrationPointIndex, IntegrationMethod ThisMethod, ProcessInfo& rCurrentProcessInfo);

    void CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, bool CalculateStiffnessMatrixFlag, bool CalculateResidualVectorFlag);

    void CalculateLocalSystem(
      MatrixType& rLeftHandSideMatrix,
      VectorType& rRightHandSideVector,
      ProcessInfo& rCurrentProcessInfo) override;

    ///@name Serialization
    ///@{
    friend class Serializer;


    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    AcousticPMLElement(): AcousticElement()
    {

    }

}; // Class AcousticElement


} // namespace Kratos.
#endif
