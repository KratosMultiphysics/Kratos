//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

#pragma once

// System includes
#include "includes/define.h"
#include "includes/condition.h"

// External includes

// Project includes
#include "iga_application_variables.h"

namespace Kratos
{
/// Condition for penalty support condition
class KRATOS_API(IGA_APPLICATION) SBMLaplacianCondition
    : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer definition of SBMLaplacianCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SBMLaplacianCondition);

    /// Size types
    using SizeType = std::size_t;
    using IndexType = std::size_t;

    // enum
    enum class BoundaryConditionType {
        Dirichlet,
        Neumann,
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with Id and geometry
    SBMLaplacianCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : Condition(NewId, pGeometry)
    {};

    /// Constructor with Id, geometry and property
    SBMLaplacianCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : Condition(NewId, pGeometry, pProperties)
    {};

    /// Default constructor
    SBMLaplacianCondition() : Condition()
    {};

    /// Destructor
    virtual ~SBMLaplacianCondition() override
    {};

    ///@}
    ///@name Life Cycle
    ///@{

    /// Create with Id, pointer to geometry and pointer to property
    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
    ) const override
    {
        return Kratos::make_intrusive<SBMLaplacianCondition>(
            NewId, pGeom, pProperties);
    };

    /// Create with Id, pointer to geometry and pointer to property
    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
    ) const override
    {
        return Kratos::make_intrusive<SBMLaplacianCondition>(
            NewId, GetGeometry().Create(ThisNodes), pProperties);
    };

    ///@}
    ///@name Operations
    ///@{

    /**
    * @brief This is called during the assembling process in order
    *        to calculate the condition right hand side matrix
    * @param rLeftHandSideMatrix the condition right hand side matrix
    * @param rCurrentProcessInfo the current process info
    */
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
    * @brief This is called during the assembling process in order
    *        to calculate the condition left hand side matrix
    * @param rLeftHandSideMatrix the condition left hand side matrix
    * @param rCurrentProcessInfo the current process info
    */
    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief This function provides a more general interface to the element.
     * @details It is designed so that rLHSvariables and rRHSvariables are
     *          passed to the element thus telling what is the desired output
     * @param rLeftHandSideMatrix container with the output Left Hand Side matrix
     * @param rRightHandSideVector container for the desired RHS output
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;
    /**
    * @brief Sets on rResult the ID's of the element degrees of freedom
    * @param rResult The vector containing the equation id
    * @param rCurrentProcessInfo The current process info instance
    */
    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const override;

    /**
    * @brief Sets on rConditionDofList the degrees of freedom of the considered element geometry
    * @param rElementalDofList The vector containing the dof of the element
    * @param rCurrentProcessInfo The current process info instance
    */
    void GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const override;

    ///@}
    ///@name Check
    ///@{

    /// Performs check if Penalty factor is provided.
    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    BoundaryConditionType GetBoundaryConditionType(const std::string& rType);

    /**
     * @brief Compute the factorial of a positive integer n
     * 
     * @param n 
     * @return unsigned long long 
     */
    unsigned long long factorial(IndexType n); 

    /**
     * @brief compute the Taylor expansion for apply the Shifted Boundary Method in 2D
     * @param derivative 
     * @param dx 
     * @param k 
     * @param dy 
     * @param n_k 
     * @return double 
     */
    double computeTaylorTerm(
        double derivative, 
        double dx, IndexType k, 
        double dy, IndexType n_k);

    /**
     * @brief compute the Taylor expansion for apply the Shifted Boundary Method in 3D
     * @param derivative 
     * @param dx 
     * @param k 
     * @param dy 
     * @param n_k 
     * @return double 
     */
    double computeTaylorTerm3D(
        double derivative, 
        double dx, IndexType k_x, 
        double dy, IndexType k_y, 
        double dz, IndexType k_z);

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "\"SBMLaplacianCondition\" #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "\"SBMLaplacianCondition\" #" << Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        pGetGeometry()->PrintData(rOStream);
    }

    ///@}

private:
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
    }

    virtual void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
    }
    
    /**
     * @brief Calculate All for Dirichlet boundary conditions
     * 
     * @param rLeftHandSideMatrix 
     * @param rRightHandSideVector 
     * @param rCurrentProcessInfo 
     */
    void CalculateAllDirichlet(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
    );

    /**
     * @brief Calculate All for Neumann boundary conditions
     * 
     * @param rLeftHandSideMatrix 
     * @param rRightHandSideVector 
     * @param rCurrentProcessInfo 
     */
    void CalculateAllNeumann(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
    );

    // sbm variables
    array_1d<double, 3> mNormalParameterSpace;
    Matrix mHsum = ZeroMatrix(1, this->GetGeometry().size());
    Vector mDistanceVector;
    std::vector<Matrix> mShapeFunctionDerivatives;
    IndexType mBasisFunctionsOrder;

    ///@}

}; // Class SBMLaplacianCondition

}  // namespace Kratos.

