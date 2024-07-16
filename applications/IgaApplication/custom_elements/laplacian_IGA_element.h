// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#if !defined(KRATOS_LAPLACIAN_IGA_ELEMENT_H_INCLUDED )
#define  KRATOS_LAPLACIAN_IGA_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"


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

class LaplacianIGAElement : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of LaplacianIGAElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(LaplacianIGAElement);

    typedef Element BaseType;

    // static constexpr std::size_t NumNodes = TDim + 1;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    LaplacianIGAElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry);

    LaplacianIGAElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~LaplacianIGAElement();

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



    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const override;

    void GetDofList(DofsVectorType& ElementalDofList, const ProcessInfo& CurrentProcessInfo) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    LaplacianIGAElement() : Element()
    {
    }

    // New functions NICO



    // void CalculateRightHandSide(
    //     VectorType& rRightHandSideVector,
    //     const ProcessInfo& rCurrentProcessInfo) override
    // {
    //     const SizeType number_of_nodes = GetGeometry().size();
    //     const SizeType mat_size = number_of_nodes * 3;

    //     if (rRightHandSideVector.size() != mat_size)
    //         rRightHandSideVector.resize(mat_size);
    //     noalias(rRightHandSideVector) = ZeroVector(mat_size);

    //     MatrixType left_hand_side_matrix;

    //     KRATOS_WATCH("RHS1")
        
    //     // CalculateAll(left_hand_side_matrix, rRightHandSideVector,
    //     //     rCurrentProcessInfo, false, true);
    //     MatrixType temp(0,0);
    //     CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
    // }


    // void CalculateLeftHandSide(
    //     MatrixType& rLeftHandSideMatrix,
    //     const ProcessInfo& rCurrentProcessInfo) override
    // {
    //     KRATOS_WATCH("LHS0")
    //     const SizeType number_of_nodes = GetGeometry().size();
    //     const SizeType mat_size = number_of_nodes * 3;

    //     VectorType right_hand_side_vector;

    //     if (rLeftHandSideMatrix.size1() != mat_size)
    //         rLeftHandSideMatrix.resize(mat_size, mat_size);
    //     noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);

    //     // CalculateAll(rLeftHandSideMatrix, right_hand_side_vector,
    //     //     rCurrentProcessInfo, true, false);
        
    //     KRATOS_WATCH("LHS1')

    //     VectorType temp(0);
    //     CalculateLocalSystem(rLeftHandSideMatrix, temp, rCurrentProcessInfo);
    // }


    // void CalculateLocalSystem(
    //     MatrixType& rLeftHandSideMatrix,
    //     VectorType& rRightHandSideVector,
    //     const ProcessInfo& rCurrentProcessInfo) override
    // {
    //     const SizeType number_of_nodes = GetGeometry().size();
    //     const SizeType mat_size = number_of_nodes * 3;

    //     if (rRightHandSideVector.size() != mat_size)
    //         rRightHandSideVector.resize(mat_size);
    //     noalias(rRightHandSideVector) = ZeroVector(mat_size);

    //     if (rLeftHandSideMatrix.size1() != mat_size)
    //         rLeftHandSideMatrix.resize(mat_size, mat_size);
    //     noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);

    //     CalculateAll(rLeftHandSideMatrix, rRightHandSideVector,
    //         rCurrentProcessInfo, true, true);
    // }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{
    
    IntegrationMethod GetIntegrationMethod() const override;

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

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    
    // /// Calculates LHS and RHS dependent on flags
    // void CalculateAll(
    //     MatrixType& rLeftHandSideMatrix,
    //     VectorType& rRightHandSideVector,
    //     const ProcessInfo& rCurrentProcessInfo,
    //     const bool CalculateStiffnessMatrixFlag,
    //     const bool CalculateResidualVectorFlag
    // ) const;


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    }

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    // std::vector<std::size_t> GetSurrogateFacesIds();

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
    //LaplacianIGAElement& operator=(const LaplacianIGAElement& rOther);

    /// Copy constructor.
    //LaplacianIGAElement(const LaplacianIGAElement& rOther);

    ///@}

}; // Class LaplacianIGAElement

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
				    LaplacianIGAElement& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
				    const LaplacianIGAElement& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);
      return rOStream;
    }*/
///@}

}  // namespace Kratos.

#endif // KRATOS_LAPLACIAN_SHIFTED_BOUNDARY_ELEMENT_H_INCLUDED  defined