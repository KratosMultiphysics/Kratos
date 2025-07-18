//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Minas Apostolakis
//                   
//                   
//

#if !defined(KRATOS_COUPLING_SOLIDSHELL3P_PENALTY_CONDITION_H_INCLUDED )
#define  KRATOS_COUPLING_SOLIDSHELL3P_PENALTY_CONDITION_H_INCLUDED

// System includes
#include "includes/define.h"
#include "includes/condition.h"

// External includes

// Project includes
#include "iga_application_variables.h"
#include "custom_utilities/iga_flags.h"
#include "geometries/coupling_geometry.h"


namespace Kratos
{
    /// Penalty factor based coupling condition.
/** This condition can be used to apply continuity between different
*   discretizations with the penalty approach.
*
*   The aproach is described in https://doi.org/10.1186/s40323-018-0109-4
*   Eq 15 ff
*
*   The condition needs a PENALTY as parameter in the Properties.
*   The Geometry needs to be of type CouplingMasterSlave and must have
*   at least one slave geometry.
*   The continuities can be enabled or disabled with the
*   FIX_DISPLACEMENT_{dir} flags.
*/
class CouplingSolidShell3pPenaltyCondition
    : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of CouplingSolidShell3pPenaltyCondition
    KRATOS_CLASS_POINTER_DEFINITION(CouplingSolidShell3pPenaltyCondition);

    /// Size types
    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with Id and geometry : exaclty as CouplingPenaltyCondition
    CouplingSolidShell3pPenaltyCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry) :
        Condition(NewId, pGeometry)
    {
    };
    
    /// Constructor with Id, geometry and property : exaclty as CouplingPenaltyCondition
    CouplingSolidShell3pPenaltyCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties) :
        Condition(NewId, pGeometry, pProperties)
    {   

    };

    /// Default constructor
    CouplingSolidShell3pPenaltyCondition()
        : Condition()
    {
    };

    /// Destructor.
    virtual ~CouplingSolidShell3pPenaltyCondition() = default;

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
       return Kratos::make_intrusive<CouplingSolidShell3pPenaltyCondition>(
           NewId, pGeom, pProperties);
   };
   
   /// Create with Id, pointer to geometry and pointer to property
   Condition::Pointer Create(
       IndexType NewId,
       NodesArrayType const& ThisNodes,
       PropertiesType::Pointer pProperties
   ) const override
   {
       return Kratos::make_intrusive< CouplingSolidShell3pPenaltyCondition >(
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
        const ProcessInfo& rCurrentProcessInfo) override
    {
        MatrixType left_hand_side_matrix = Matrix(0, 0);
    
        CalculateAll(left_hand_side_matrix, rRightHandSideVector,
            rCurrentProcessInfo, false, true);
    }

    /**
    * @brief This is called during the assembling process in order
    *        to calculate the condition left hand side matrix
    * @param rLeftHandSideMatrix the condition left hand side matrix
    * @param rCurrentProcessInfo the current process info
    */
    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        VectorType right_hand_side_vector = Vector(0);
    
        CalculateAll(rLeftHandSideMatrix, right_hand_side_vector,
            rCurrentProcessInfo, true, false);
    }

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
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateAll(rLeftHandSideMatrix, rRightHandSideVector,
            rCurrentProcessInfo, true, true);
    }

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
    * @brief Sets on rElementalDofList the degrees of freedom of the considered element geometry
    * @param rElementalDofList The vector containing the dof of the element
    * @param rCurrentProcessInfo The current process info instance
    */
    void GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const override;

    /// Calculates left (K) and right (u) hand sides, according to the flags
    void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    );

    void DeterminantOfJacobianInitial(
        const GeometryType& rGeometry,
        Vector& rDeterminantOfJacobian);

    void OutOfPlaneDeformationFirstVariation(
        Matrix& OutOfPlaneDeformationWholeMatrix,
        const size_t& mat_size,
        const double theta,
        const array_1d<double, 3>& A1,
        const array_1d<double, 3>& A2,
        const Matrix& shape_functions_gradients_slave);

    array_1d<double, 3> Calculate_Phi_r_cross_A3(
        const array_1d<double, 3>& N_theta1, 
        const array_1d<double, 3>& N_theta2,
        const array_1d<double, 3>& A1, 
        const array_1d<double, 3>& A2);


    ///@}
    ///@name Check
    ///@{

    /// Performs check if Penalty factor is provided.
    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "\"CouplingSolidShell3pPenaltyCondition\" #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "\"CouplingSolidShell3pPenaltyCondition\" #" << Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {
        pGetGeometry()->PrintData(rOStream);
    }

    ///@}

private:

    std::vector<array_1d<double, 3>> _g1, _g2, _g3;
    std::vector<double> _theta3;
    

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

    ///@}

}; // Class CouplingSolidShell3pPenaltyCondition

}  // namespace Kratos.

#endif // KRATOS_COUPLING_PENALTY_CONDITION_H_INCLUDED  defined