//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_POINT_LOAD_CONDITION_H_INCLUDED)
#define  KRATOS_POINT_LOAD_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_conditions/load_conditions/load_condition.hpp"


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

/// Point Load Condition for 3D and 2D geometries. (base class)

class KRATOS_API(SOLID_MECHANICS_APPLICATION) PointLoadCondition
    : public LoadCondition
{
public:

    ///@name Type Definitions
    ///@{
    // Counted pointer of PointLoadCondition
    KRATOS_CLASS_POINTER_DEFINITION( PointLoadCondition );
    ///@}

    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PointLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry );

    PointLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    /// Copy constructor
    PointLoadCondition( PointLoadCondition const& rOther);

    /// Destructor
    ~PointLoadCondition() override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * creates a new condition pointer
     * @param NewId: the ID of the new condition
     * @param ThisNodes: the nodes of the new condition
     * @param pProperties: the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Create(IndexType NewId,
			      NodesArrayType const& ThisNodes,
			      PropertiesType::Pointer pProperties ) const override;


    /**
     * clones the selected condition variables, creating a new one
     * @param NewId: the ID of the new condition
     * @param ThisNodes: the nodes of the new condition
     * @param pProperties: the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Clone(IndexType NewId,
			     NodesArrayType const& ThisNodes) const override;


    //************************************************************************************
    //************************************************************************************
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    int Check( const ProcessInfo& rCurrentProcessInfo ) override;

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

    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Point Load Condition #" << Id();
        return buffer.str();
    }

    /// Print information about this object.

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Point Load Condition #" << Id();
    }

    /// Print object's data.

    void PrintData(std::ostream& rOStream) const override
    {
        pGetGeometry()->PrintData(rOStream);
    }

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
    PointLoadCondition() {};
    ///@}
    ///@name Protected Operators
    ///@{
    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * Initialize System Matrices
     */
    void InitializeConditionVariables(ConditionVariables& rVariables,
					    const ProcessInfo& rCurrentProcessInfo) override;


    /**
     * Calculate Condition Kinematics
     */
    void CalculateKinematics(ConditionVariables& rVariables,
				     const double& rPointNumber) override;

    /**
     * Calculate the External Load of the Condition
     */
    void CalculateExternalLoad(ConditionVariables& rVariables) override;


    /**
     * Calculates the condition contributions
     */
    void CalculateConditionSystem(LocalSystemComponents& rLocalSystem,
					  const ProcessInfo& rCurrentProcessInfo) override;


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


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

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;


}; // class PointLoadCondition.

} // namespace Kratos.

#endif // KRATOS_POINT_LOAD_CONDITION_H_INCLUDED defined
