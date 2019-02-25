//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:                LMonforte $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

#if !defined(KRATOS_EP_AXISYM_POINT_RIGID_CONTACT_PENALTY_2D_CONDITION_H_INCLUDED )
#define  KRATOS_EP_AXISYM_POINT_RIGID_CONTACT_PENALTY_2D_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_conditions/rigid_contact/EP_point_rigid_contact_penalty_2D_condition.hpp"

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
class KRATOS_API(CONTACT_MECHANICS_APPLICATION) EPAxisymPointRigidContactPenalty2DCondition
    : public EPPointRigidContactPenalty2DCondition
{
public:

   ///@name Type Definitions

    ///Tensor order 1 definition
    typedef Vector VectorType;

    ///@{
    // Counted pointer of PointRigidContactPenalty2DCondition
    KRATOS_CLASS_POINTER_DEFINITION( EPAxisymPointRigidContactPenalty2DCondition );
    ///@}

    ///@}
    ///@name Life Cycle
    ///@{

    /// Serialization constructor
    EPAxisymPointRigidContactPenalty2DCondition(){};

    /// Default constructor
    EPAxisymPointRigidContactPenalty2DCondition(IndexType NewId, GeometryType::Pointer pGeometry);

    EPAxisymPointRigidContactPenalty2DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    EPAxisymPointRigidContactPenalty2DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties, SpatialBoundingBox::Pointer pRigidWall);


    /// Copy constructor
    EPAxisymPointRigidContactPenalty2DCondition( EPAxisymPointRigidContactPenalty2DCondition const& rOther);


    /// Destructor.
    virtual ~EPAxisymPointRigidContactPenalty2DCondition();


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
    Condition::Pointer Create(IndexType NewId, NodesArrayType const&
                              ThisNodes,  PropertiesType::Pointer pProperties) const override;


    /**
     * clones the selected condition variables, creating a new one
     * @param NewId: the ID of the new condition
     * @param ThisNodes: the nodes of the new condition
     * @param pProperties: the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Clone(IndexType NewId,
			     NodesArrayType const& ThisNodes) const override;


    ///@}
    ///@name Access
    ///@{
    ///@}
    ///@name Inquiry
    ///@{
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


     /**
     * Calculate Condition Kinematics
     */
    virtual void CalculateKinematics(ConditionVariables& rVariables,
				     const ProcessInfo& rCurrentProcessInfo,
				     const double& rPointNumber) override;


    /**
     * Calculation and addition of the matrices of the LHS
     */
    virtual void CalculateAndAddLHS(LocalSystemComponents& rLocalSystem,
                                    ConditionVariables& rVariables,
                                    double& rIntegrationWeight) override;

    /**
     * Calculation and addition of the vectors of the RHS
     */
    virtual void CalculateAndAddRHS(LocalSystemComponents& rLocalSystem,
                                    ConditionVariables& rVariables,
                                    double& rIntegrationWeight) override;

    /**
     * Calculation of the contidion radius (axisymmetry)
     */
    void CalculateRadius(double & rCurrentRadius,
			 double & rReferenceRadius);


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

    virtual void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, EPPointRigidContactPenalty2DCondition )
    }

    virtual void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, EPPointRigidContactPenalty2DCondition )
    }


}; // Class EPAxisymPointRigidContactPenalty2DCondition

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
				    EPAxisymPointRigidContactPenalty2DCondition& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
				    const EPAxisymPointRigidContactPenalty2DCondition& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

}  // namespace Kratos.

#endif // KRATOS_AXISYM_POINT_RIGID_CONTACT_PENALTY_2D_CONDITION_H_INCLUDED  defined
