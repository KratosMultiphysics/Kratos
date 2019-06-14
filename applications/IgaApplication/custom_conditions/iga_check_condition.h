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


#if !defined(KRATOS_IGA_CHECK_CONDITION_H_INCLUDED )
#define  KRATOS_IGA_CHECK_CONDITION_H_INCLUDED


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
* @class IgaCheckCondition
* @ingroup IgaApplication
* @brief This class is used to measure on specific points on the
*        NURBS-described surface.
*/
class IgaCheckCondition
    : public IgaBaseCondition
{
public:
    ///@name Type Definitions
    ///@{
    /// Counted pointer of IgaCheckCondition
    KRATOS_CLASS_POINTER_DEFINITION( IgaCheckCondition );
    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor void
    IgaCheckCondition()
    {};

    // Constructor using an array of nodes
    IgaCheckCondition(IndexType NewId, GeometryType::Pointer pGeometry)
        :IgaBaseCondition(NewId, pGeometry)
    {};

    // Constructor using an array of nodes with properties
    IgaCheckCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        :IgaBaseCondition(NewId, pGeometry, pProperties)
    {};

    // Destructor
    ~IgaCheckCondition() override
    {};

    ///@}
    ///@name Operations
    ///@{

    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<IgaCheckCondition >(
            NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
    ) const override
    {
        return Kratos::make_intrusive<IgaCheckCondition>(
            NewId, pGeom, pProperties);
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
        const bool CalculateResidualVectorFlag) override;


    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "\"IgaCheckCondition\" #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "\"IgaCheckCondition\" #" << Id();
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

    /** Gets the number of dofs per node.
    * If a derived condition has more degrees of freedom this function has to be overridden
    *
    * @return Number of dofs per node.
    */
    inline std::size_t DofsPerNode() const override
    {
        return 3;
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
}; // Class IgaCheckCondition

}  // namespace Kratos.

#endif // KRATOS_IGA_CHECK_CONDITION_H_INCLUDED  defined


