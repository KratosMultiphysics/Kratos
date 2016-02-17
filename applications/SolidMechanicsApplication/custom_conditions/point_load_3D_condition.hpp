//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_POINT_LOAD_3D_CONDITION_H_INCLUDED )
#define  KRATOS_POINT_LOAD_3D_CONDITION_H_INCLUDED



// System includes

// External includes

// Project includes
#include "custom_conditions/force_load_condition.hpp"


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

/// Force Load Condition for 3D and 2D geometries. (base class)

/**
 * Implements a Force Load definition for structural analysis.
 * This works for arbitrary geometries in 3D and 2D (base class)
 */
class KRATOS_API(SOLID_MECHANICS_APPLICATION) PointLoad3DCondition
    : public ForceLoadCondition
{
public:

    ///@name Type Definitions
    ///@{
    // Counted pointer of PointLoad3DCondition
    KRATOS_CLASS_POINTER_DEFINITION( PointLoad3DCondition );
    ///@}

    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PointLoad3DCondition( IndexType NewId, GeometryType::Pointer pGeometry );

    PointLoad3DCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    /// Copy constructor
    PointLoad3DCondition( PointLoad3DCondition const& rOther);

    /// Destructor
    virtual ~PointLoad3DCondition();

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
			      PropertiesType::Pointer pProperties ) const;


    /**
     * clones the selected condition variables, creating a new one
     * @param NewId: the ID of the new condition
     * @param ThisNodes: the nodes of the new condition
     * @param pProperties: the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Clone(IndexType NewId, 
			     NodesArrayType const& ThisNodes) const;


    //************************************************************************************
    //************************************************************************************
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    virtual int Check( const ProcessInfo& rCurrentProcessInfo );

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

    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Point Load 3D Condition #" << Id();
        return buffer.str();
    }

    /// Print information about this object.

    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Point Load 3D Condition #" << Id();
    }

    /// Print object's data.

    virtual void PrintData(std::ostream& rOStream) const
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
    PointLoad3DCondition() {};
    ///@}
    ///@name Protected Operators
    ///@{
    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * Initialize System Matrices
     */
    virtual void InitializeGeneralVariables(GeneralVariables& rVariables, 
					    const ProcessInfo& rCurrentProcessInfo);


    /**
     * Calculate Condition Kinematics
     */
    virtual void CalculateKinematics(GeneralVariables& rVariables, 
				     const double& rPointNumber);

    /**
     * Calculation of the Vector Force of the Condition
     */
    virtual Vector& CalculateVectorForce(Vector& rVectorForce, GeneralVariables& rVariables);


    /**
     * Calculation of the Integration Weight
     */
    virtual double& CalculateIntegrationWeight(double& rIntegrationWeight);


    /**
     * Calculates the condition contributions
     */
    virtual void CalculateConditionSystem(LocalSystemComponents& rLocalSystem,
					  const ProcessInfo& rCurrentProcessInfo);


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

    virtual void save(Serializer& rSerializer) const;

    virtual void load(Serializer& rSerializer);


}; // class PointLoad3DCondition.

} // namespace Kratos.

#endif // KRATOS_POINT_LOAD_3D_CONDITION_H_INCLUDED defined 
