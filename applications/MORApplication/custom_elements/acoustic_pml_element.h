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
    typedef Kratos::Vector ShapeFunctionsType;

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
    // Element::Pointer Create(
    //     IndexType NewId,
    //     GeometryType::Pointer pGeom,
    //     PropertiesType::Pointer pProperties
    //     ) const override;

    // /**
    //  * @brief Creates a new element
    //  * @param NewId The Id of the new created element
    //  * @param ThisNodes The array containing nodes
    //  * @param pProperties The pointer to property
    //  * @return The pointer to the created element
    //  */
    // Element::Pointer Create(
    //     IndexType NewId,
    //     NodesArrayType const& ThisNodes,
    //     PropertiesType::Pointer pProperties
    //     ) const override;

    // /**
    //  * @brief It creates a new element pointer and clones the previous element data
    //  * @param NewId the ID of the new element
    //  * @param ThisNodes the nodes of the new element
    //  * @param pProperties the properties assigned to the new element
    //  * @return a Pointer to the new element
    //  */
    // Element::Pointer Clone (
    //     IndexType NewId,
    //     NodesArrayType const& rThisNodes
    //     ) const override;

    // /**
    //  * @brief This function provides the place to perform checks on the completeness of the input.
    //  * @details It is designed to be called only once (or anyway, not often) typically at the beginning
    //  * of the calculations, so to verify that nothing is missing from the input
    //  * or that no common error is found.
    //  * @param rCurrentProcessInfo The current process info instance
    //  */
    // int Check(const ProcessInfo& rCurrentProcessInfo) override;

    // ///@}
    // ///@name Access
    // ///@{

    // ///@}
    // ///@name Inquiry
    // ///@{
    // ///@}
    // ///@name Input and output
    // ///@{

    // /// Turn back information as a string.
    // std::string Info() const override
    // {
    //     std::stringstream buffer;
    //     buffer << "Acoustic PML Element #" << Id() << "\n";
    //     return buffer.str();
    // }

    // /// Print information about this object.
    // void PrintInfo(std::ostream& rOStream) const override
    // {
    //     rOStream << "Acoustic PML Element #" << Id() << "\n";
    // }

    // // Print object's data.
    // void PrintData(std::ostream& rOStream) const override
    // {
    //     pGetGeometry()->PrintData(rOStream);
    // }

    // void EquationIdVector(
    //         EquationIdVectorType& rResult,
    //         ProcessInfo& rCurrentProcessInfo) override;

    // void GetDofList(
    //     DofsVectorType& rElementalDofList,
    //     ProcessInfo& rCurrentProcessInfo) override;

    // void Initialize() override;


protected:
    // //* //@name Protected static Member Variables
    // ///@{
    // struct KinematicVariables
    // {
    //     Vector  N;
    //     Matrix  B;
    //     double  detF;
    //     Matrix  F;
    //     double  detJ0;
    //     Matrix  J0;
    //     Matrix  InvJ0;
    //     Matrix  DN_DX;
    //     Vector Displacements;

    //     /**
    //      * The default constructor
    //      * @param NumberOfNodes The number of nodes in the element
    //      */
    //     KinematicVariables(const SizeType NumberOfNodes, const SizeType Dimension)
    //     {
    //         detF = 1.0;
    //         detJ0 = 1.0;
    //         N = ZeroVector(NumberOfNodes);//@name Protected static Member Variables
    // ///@{
    // struct KinematicVariables
    // {
    //     Vector  N;
    //     Matrix  B;
    //     double  detF;
    //     Matrix  F;
    //     double  detJ0;
    //     Matrix  J0;
    //     Matrix  InvJ0;
    //     Matrix  DN_DX;
    //     Vector Displacements;

    //     /**
    //      * The default constructor
    //      * @param NumberOfNodes The number of nodes in the element
    //      */
    //     KinematicVariables(const SizeType NumberOfNodes, const SizeType Dimension)
    //     {
    //         detF = 1.0;
    //         detJ0 = 1.0;
    //         N = ZeroVector(NumberOfNodes);
    //         B = ZeroMatrix(Dimension * NumberOfNodes);
    //         F = IdentityMatrix(Dimension);
    //         DN_DX = ZeroMatrix(NumberOfNodes, Dimension);
    //         J0 = ZeroMatrix(Dimension, Dimension);
    //         InvJ0 = ZeroMatrix(Dimension, Dimension);
    //         Displacements = ZeroVector(Dimension * NumberOfNodes);
    //     }
    // };
    // ///@}
    // ///@name Protected member Variables
    // ///@{

    // IntegrationMethod mThisIntegrationMethod;

    // ///@}
    // ///@name Protected Operators
    // ///@{

    // AcousticElement() : Element()
    // {
    // }

    // ///@}
    // ///@name Protected Operations
    // ///@{

    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        ProcessInfo& rCurrentProcessInfo) override;

    // ShapeFunctionDerivativesArrayType& ShapeFunctionsIntegrationPointsGradientsImag(
    //     ShapeFunctionDerivativesArrayType& rResult, CoordinatesArrayType imag_coordinates, Vector& determinants_of_jacobian, IntegrationMethod ThisMethod ) const
    // {
    //     const GeometryType& geom = GetGeometry();
    //     const unsigned int integration_points_number = geom.IntegrationPointsNumber( ThisMethod );

    //     if ( integration_points_number == 0 )
    //         KRATOS_ERROR << "This integration method is not supported " << geom << std::endl;

    //     if ( rResult.size() != integration_points_number )
    //         rResult.resize(  geom.IntegrationPointsNumber( ThisMethod ), false  );
    //     if ( determinants_of_jacobian.size() != integration_points_number )
    //         determinants_of_jacobian.resize(  geom.IntegrationPointsNumber( ThisMethod ), false  );

    //     //calculating the local gradients
    //     const AcousticPMLElement::ShapeFunctionDerivativesArrayType& DN_De = geom.ShapeFunctionsLocalGradients( ThisMethod );

    //     //loop over all integration points
    //     Matrix J(geom.WorkingSpaceDimension(),geom.LocalSpaceDimension());
    //     Matrix Jinv(geom.WorkingSpaceDimension(),geom.LocalSpaceDimension());
    //     double DetJ;
    //     for ( unsigned int pnt = 0; pnt < integration_points_number; pnt++ )
    //     {
    //         if(rResult[pnt].size1() != geom.WorkingSpaceDimension() ||  rResult[pnt].size2() != geom.LocalSpaceDimension())
    //             rResult[pnt].resize( (geom).size(), geom.LocalSpaceDimension(), false );
    //         geom.Jacobian(J, imag_coordinates);
    //         MathUtils<double>::GeneralizedInvertMatrix( J, Jinv, DetJ );
    //         noalias(rResult[pnt]) =  prod( DN_De[pnt], Jinv );
    //         determinants_of_jacobian[pnt] = DetJ;
    //     }

    //     return rResult;
    // }

    void ComplexJacobian( ComplexMatrix& rResult, IndexType IntegrationPointIndex, IntegrationMethod ThisMethod );


    // void CalculateRightHandSide(
    //   VectorType& rRightHandSideVector,
    //   ProcessInfo& rCurrentProcessInfo) override;

    // void CalculateLocalSystem(
    //   MatrixType& rLeftHandSideMatrix,
    //   VectorType& rRightHandSideVector,
    //   ProcessInfo& rCurrentProcessInfo) override;

    // void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) override;

    // void CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo) override;


    // void SetIntegrationMethod(const IntegrationMethod& ThisIntegrationMethod)
    // {
    //      mThisIntegrationMethod = ThisIntegrationMethod;
    // }
    //         J0 = ZeroMatrix(Dimension, Dimension);
    //         InvJ0 = ZeroMatrix(Dimension, Dimension);
    //         Displacements = ZeroVector(Dimension * NumberOfNodes);
    //     }
    // };
    // ///@}
    // ///@name Protected member Variables
    // ///@{

    // IntegrationMethod mThisIntegrationMethod;

    // ///@}
    // ///@name Protected Operators
    // ///@{

    // AcousticElement() : Element()
    // {
    // }

    // ///@}
    // ///@name Protected Operations
    // ///@{

    // void CalculateLeftHandSide(
    // MatrixType& rLeftHandSideMatrix,
    // ProcessInfo& rCurrentProcessInfo) override;

    // void CalculateRightHandSide(
    //   VectorType& rRightHandSideVector,
    //   ProcessInfo& rCurrentProcessInfo) override;

    // void CalculateLocalSystem(
    //   MatrixType& rLeftHandSideMatrix,
    //   VectorType& rRightHandSideVector,
    //   ProcessInfo& rCurrentProcessInfo) override;

    // void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) override;

    // void CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo) override;


    // void SetIntegrationMethod(const IntegrationMethod& ThisIntegrationMethod)
    // {
    //      mThisIntegrationMethod = ThisIntegrationMethod;
    // } */

    // double CalculateDerivativesOnReferenceConfiguration(
    //     Matrix& rJ0,
    //     Matrix& rInvJ0,
    //     Matrix& rDN_DX,
    //     const IndexType PointNumber,
    //     IntegrationMethod ThisIntegrationMethod) const;

    // void CalculateKinematicVariables(
    //     KinematicVariables& rThisKinematicVariables,
    //     const IndexType PointNumber,
    //     const GeometryType::IntegrationMethod& rIntegrationMethod
    //     );
    // AcousticPMLElement& operator=(const AcousticPMLElement& rOther);
    // /// Copy constructor.
    // AcousticPMLElement(const AcousticPMLElement& rOther); Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  AcSizeType
    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    // A private default constructor necessary for serialization

    //void save(Serializer& rSerializer) const override;
    //void load(Serializer& rSerializer) override;

    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    AcousticPMLElement(): AcousticElement()
    {

    }

    /// Assignment operator.
    //AcousticPMLElement& operator=(const AcousticPMLElement& rOther);
    /// Copy constructor.
    //AcousticPMLElement(const AcousticPMLElement& rOther);
    ///@}

}; // Class AcousticElement

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif  // KRATOS_ACOUSTIC_ELEMENT_H_INCLUDED  defined
