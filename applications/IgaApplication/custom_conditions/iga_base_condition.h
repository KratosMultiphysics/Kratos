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


#if !defined(KRATOS_IGA_BASE_CONDITION_H_INCLUDED )
#define  KRATOS_IGA_BASE_CONDITION_H_INCLUDED


// System includes

// External includes

// Project includes
#include "iga_application_variables.h"
#include "custom_utilities/geometry_utilities/iga_geometry_utilities.h"

#include "includes/condition.h"

#include "includes/checks.h"
#include "includes/define.h"
#include "includes/variables.h"

namespace Kratos
{
/**
* @class IgaBaseCondition
* @ingroup IgaApplication
* @brief This is base clase used to define discrete elements
*/
class IgaBaseCondition
    : public Condition
{
public:
    ///@name Type Definitions
    ///@{
    /// Counted pointer of IgaBaseCondition
    KRATOS_CLASS_POINTER_DEFINITION( IgaBaseCondition );
    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor void
    IgaBaseCondition()
    {};

    // Constructor using an array of nodes
    IgaBaseCondition(IndexType NewId, GeometryType::Pointer pGeometry) :Condition(NewId, pGeometry)
    {};

    // Constructor using an array of nodes with properties
    IgaBaseCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        :Condition(NewId, pGeometry, pProperties)
    {};

    // Destructor
    ~IgaBaseCondition() override
    {};

    ///@}
    ///@name Operations
    ///@{

    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override
    {
        KRATOS_ERROR << "Trying to create a \"IgaBaseCondition\"" << std::endl;
    }

    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
    ) const override
    {
        KRATOS_ERROR << "Trying to create a \"IgaBaseCondition\"" << std::endl;
    };

    /**
    * @brief This function provides a more general interface to the element.
    * @details It is designed so that rLHSvariables and rRHSvariables are passed to the element thus telling what is the desired output
    * @param rLeftHandSideMatrix container with the output Left Hand Side matrix
    * @param rRightHandSideVector container for the desired RHS output
    * @param rCurrentProcessInfo the current process info instance
    */
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo
    ) override;

    /**
    * @brief This is called during the assembling process in order to calculate the elemental right hand side vector only
    * @param rRightHandSideVector the elemental right hand side vector
    * @param rCurrentProcessInfo the current process info instance
    */
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo
    ) override;

    /**
    * @brief This is called during the assembling process in order to calculate the elemental left hand side matrix only
    * @param rLeftHandSideMatrix the elemental left hand side matrix
    * @param rCurrentProcessInfo the current process info instance
    */
    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        ProcessInfo& rCurrentProcessInfo
    ) override;

    /**
    * @brief Sets on rValues the nodal displacements
    * @param rValues The values of displacements
    * @param Step The step to be computed
    */
    void GetValuesVector(
        Vector& rValues,
        int Step = 0
    ) override;

    /**
    * @brief Sets on rValues the nodal velocities
    * @param rValues The values of velocities
    * @param Step The step to be computed
    */
    void GetFirstDerivativesVector(
        Vector& rValues,
        int Step = 0
    ) override;

    /**
    * @brief Sets on rValues the nodal accelerations
    * @param rValues The values of accelerations
    * @param Step The step to be computed
    */
    void GetSecondDerivativesVector(
        Vector& rValues,
        int Step = 0
    ) override;


    /********************************************************************/
    /*    Calculate                                                     */
    /********************************************************************/
    /**
    * @brief Calculate a double array_1d on the Condition
    * @param rVariable The variable we want to get
    * @param rOutput The values obtained int the integration points
    * @param rCurrentProcessInfo the current process info instance
    */
    void Calculate(
        const Variable<array_1d<double, 3>>& rVariable,
        array_1d<double, 3>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
    ) override;

    /**
    * @brief Calculate a Vector Variable on the Condition
    * @param rVariable The variable we want to get
    * @param rOutput The values obtained int the integration points
    * @param rCurrentProcessInfo the current process info instance
    */
    void Calculate(
        const Variable<Vector>& rVariable,
        Vector& rOutput,
        const ProcessInfo& rCurrentProcessInfo
    ) override;


    /********************************************************************/
    /*    SetValuesOnIntegrationPoints                                   */
    /********************************************************************/
    /**
    * @brief Set a Constitutive Law Value on the Condition
    * @param rVariable The variable we want to set
    * @param rValues The values to set in the integration points
    * @param rCurrentProcessInfo the current process info instance
    */
    void SetValueOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues,
        const ProcessInfo& rCurrentProcessInfo
    ) override;



    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "\"IgaBaseCondition\" #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "\"IgaBaseCondition\" #" << Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const {
        pGetGeometry()->PrintData(rOStream);
    }

    /**
    * This method provides the place to perform checks on the completeness of the input
    * and the compatibility with the problem options as well as the contitutive laws selected
    * It is designed to be called only once (or anyway, not often) typically at the beginning
    * of the calculations, so to verify that nothing is missing from the input
    * or that no common error is found.
    * @param rCurrentProcessInfo
    * this method is: MANDATORY
    */

    virtual int Check(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        const SizeType number_of_nodes = this->GetGeometry().size();
        const SizeType dimension = this->GetGeometry().WorkingSpaceDimension();

        // Verify that the variables are correctly initialized
        KRATOS_CHECK_VARIABLE_KEY(SHAPE_FUNCTION_VALUES)

        // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
        for (IndexType i = 0; i < number_of_nodes; i++) {
            NodeType &rnode = this->GetGeometry()[i];
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT, rnode)

            KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, rnode)
            KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, rnode)
            KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z, rnode)
        }

        return 0;

        KRATOS_CATCH("");
    }
protected:

    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
    * This functions calculates both the RHS and the LHS
    * @param rLeftHandSideMatrix: The LHS
    * @param rRightHandSideVector: The RHS
    * @param rCurrentProcessInfo: The current process info instance
    * @param CalculateStiffnessMatrixFlag: The flag to set if compute the LHS
    * @param CalculateResidualVectorFlag: The flag to set if compute the RHS
    */
    virtual void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag);

    /** Gets the number of dofs per node.
    * If a derived condition has more degrees of freedom this function has to be overridden
    *
    * @return Number of dofs per node.
    */
    virtual inline std::size_t DofsPerNode() const
    {
        KRATOS_ERROR << "Trying to call DofsPerNode from \"IgaBaseCondition\". This function has to be derived in each class" << std::endl;
    }

    /** Gets the number of nodes.
    *
    * @return Number of nodes.
    */
    std::size_t inline NumberOfNodes() const
    {
        return GetGeometry().size();
    }

    /** Gets the number of degrees of freedom.
    *
    * @return Number of degrees of freedom.
    */
    std::size_t inline NumberOfDofs() const
    {
        return NumberOfNodes() * DofsPerNode();
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
    ///@}
}; // Class IgaBaseCondition

}  // namespace Kratos.

#endif // KRATOS_IGA_BASE_CONDITION_H_INCLUDED  defined


