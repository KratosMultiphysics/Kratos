//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Tobias Teschemacher
//


#if !defined(KRATOS_POST_PROCESS_CONDITION_H_INCLUDED )
#define  KRATOS_POST_PROCESS_CONDITION_H_INCLUDED


// System includes

// External includes

// Project includes
#include "iga_application_variables.h"
#include "custom_utilities/geometry_utilities/iga_geometry_utilities.h"

#include "iga_base_condition.h"

#include "includes/checks.h"
#include "includes/define.h"
#include "includes/variables.h"

namespace Kratos
{
/**
* @class PostProcessCondition
* @ingroup IgaApplication
* @brief This class is used to measure on specific points on the
*        NURBS-described surface.
*/
class PostProcessCondition
    : public Condition
{
public:
    ///@name Type Definitions
    ///@{
    /// Counted pointer of PostProcessCondition
    KRATOS_CLASS_POINTER_DEFINITION( PostProcessCondition );
    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor void
    PostProcessCondition()
    {};

    // Constructor using an array of nodes
    PostProcessCondition(IndexType NewId, GeometryType::Pointer pGeometry)
        :Condition(NewId, pGeometry)
    {};

    // Constructor using an array of nodes with properties
    PostProcessCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        :Condition(NewId, pGeometry, pProperties)
    {};

    // Destructor
    ~PostProcessCondition() override
    {};

    ///@}
    ///@name Operations
    ///@{

    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<PostProcessCondition >(
            NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
    ) const override
    {
        return Kratos::make_intrusive<PostProcessCondition>(
            NewId, pGeom, pProperties);
    };

    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo)
    {
        MatrixType left_hand_side_matrix = Matrix(0, 0);

        CalculateAll(left_hand_side_matrix, rRightHandSideVector,
            rCurrentProcessInfo, false, true);
    };

    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        ProcessInfo& rCurrentProcessInfo)
    {
        VectorType right_hand_side_vector = Vector(0);

        CalculateAll(rLeftHandSideMatrix, right_hand_side_vector,
            rCurrentProcessInfo, true, false);
    };

    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo)
    {
        CalculateAll(rLeftHandSideMatrix, rRightHandSideVector,
            rCurrentProcessInfo, true, true);
    };

    /**
    * This functions calculates both the RHS and the LHS
    * @param rLeftHandSideMatrix: The LHS
    * @param rRightHandSideVector: The RHS
    * @param rCurrentProcessInfo: The current process info instance
    * @param CalculateStiffnessMatrixFlag: The flag to set if compute the LHS
    * @param CalculateResidualVectorFlag: The flag to set if compute the RHS
    */
    void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag);

    //void CalculateOnIntegrationPoints(
    //    const Variable<bool>& rVariable,
    //    std::vector<bool>& rOutput,
    //    const ProcessInfo& rCurrentProcessInfo);

    //void CalculateOnIntegrationPoints(
    //    const Variable<int>& rVariable,
    //    std::vector<int>& rOutput,
    //    const ProcessInfo& rCurrentProcessInfo);

    //void CalculateOnIntegrationPoints(
    //    const Variable<double>& rVariable,
    //    std::vector<double>& rOutput,
    //    const ProcessInfo& rCurrentProcessInfo);

    //void CalculateOnIntegrationPoints(
    //    const Variable<array_1d<double, 3>>& rVariable,
    //    std::vector<array_1d<double, 3>>& rOutput,
    //    const ProcessInfo& rCurrentProcessInfo);

    //void CalculateOnIntegrationPoints(
    //    const Variable<array_1d<double, 6>>& rVariable,
    //    std::vector<array_1d<double, 6>>& rOutput,
    //    const ProcessInfo& rCurrentProcessInfo);

    //void CalculateOnIntegrationPoints(
    //    const Variable<Vector>& rVariable,
    //    std::vector<Vector>& rOutput,
    //    const ProcessInfo& rCurrentProcessInfo);

    void CalculateOnIntegrationPoints(
        const Variable<Matrix>& rVariable,
        std::vector<Matrix>& rOutput,
        const ProcessInfo& rCurrentProcessInfo);

    /**
    * @brief Sets on rResult the ID's of the element degrees of freedom
    * @param rResult The vector containing the equation id
    * @param rCurrentProcessInfo The current process info instance
    */
    void EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo
    ) override;

    /**
    * @brief Sets on rElementalDofList the degrees of freedom of the considered element geometry
    * @param rElementalDofList The vector containing the dof of the element
    * @param rCurrentProcessInfo The current process info instance
    */
    void GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo
    ) override;

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "\"PostProcessCondition\" #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "\"PostProcessCondition\" #" << Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const {
        pGetGeometry()->PrintData(rOStream);
    }
    ///@}

protected:

    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

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
    ///@}
}; // Class PostProcessCondition

}  // namespace Kratos.

#endif // KRATOS_POST_PROCESS_CONDITION_H_INCLUDED  defined


