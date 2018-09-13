// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes
#if !defined(KRATOS_POINT_LOAD_CONDITION_H_INCLUDED )
#define  KRATOS_POINT_LOAD_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_conditions/base_load_condition.h"
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

/// Short class definition.
/** Detail class definition.
*/

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION)  PointLoadCondition
    : public BaseLoadCondition
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of PointLoadCondition
    KRATOS_CLASS_POINTER_DEFINITION( PointLoadCondition );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PointLoadCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry
        );

    PointLoadCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties
        );

    /// Destructor.
    ~PointLoadCondition() override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const override;

    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
        ) const override;

    /**
     * Check if Rotational Dof existant
     */
    bool HasRotDof() override {return false;};

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
//      virtual String Info() const;

    /// Print information about this object.
//      virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
//      virtual void PrintData(std::ostream& rOStream) const;


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
        bool CalculateStiffnessMatrixFlag,
        bool CalculateResidualVectorFlag
        ) override;

    /**
     * It calcules the integration load for the point load
     */
    virtual double GetPointLoadIntegrationWeight();

    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{

    // A protected default constructor necessary for serialization
    PointLoadCondition() {};

    ///@}

private:
    ///@name Static Member Variables
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
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseLoadCondition );
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseLoadCondition );
    }

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    //PointLoadCondition& operator=(const PointLoadCondition& rOther);

    /// Copy constructor.
    //PointLoadCondition(const PointLoadCondition& rOther);


    ///@}

}; // Class PointLoadCondition

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
        PointLoadCondition& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
        const PointLoadCondition& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

}  // namespace Kratos.

#endif // KRATOS_POINT_LOAD_CONDITION_H_INCLUDED  defined


