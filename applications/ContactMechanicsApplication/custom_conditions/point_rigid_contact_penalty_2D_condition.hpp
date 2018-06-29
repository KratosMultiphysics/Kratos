//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

#if !defined(KRATOS_POINT_RIGID_CONTACT_PENALTY_2D_CONDITION_H_INCLUDED )
#define  KRATOS_POINT_RIGID_CONTACT_PENALTY_2D_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_conditions/point_rigid_contact_penalty_3D_condition.hpp"

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
class KRATOS_API(CONTACT_MECHANICS_APPLICATION) PointRigidContactPenalty2DCondition
    : public PointRigidContactPenalty3DCondition
{
protected:

    typedef struct
     {
        bool   Slip;
        double Sign;

        double DeltaTime;
        double PreviousTangentForceModulus;

        double FrictionCoefficient;
        double DynamicFrictionCoefficient;
        double StaticFrictionCoefficient;

     } TangentialContactVariables;



public:

   ///@name Type Definitions

    ///Tensor order 1 definition
    typedef BoundedVector<double, 3>    PointType;

    ///@{
    // Counted pointer of PointRigidContactPenalty3DCondition
    KRATOS_CLASS_POINTER_DEFINITION( PointRigidContactPenalty2DCondition );
    ///@}

    ///@}
    ///@name Life Cycle
    ///@{

    /// Serialization constructor
    PointRigidContactPenalty2DCondition(){};

    /// Default constructor
    PointRigidContactPenalty2DCondition(IndexType NewId, GeometryType::Pointer pGeometry);

    PointRigidContactPenalty2DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    PointRigidContactPenalty2DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties, SpatialBoundingBox::Pointer pRigidWall);


    /// Copy constructor
    PointRigidContactPenalty2DCondition( PointRigidContactPenalty2DCondition const& rOther);


    /// Destructor.
    virtual ~PointRigidContactPenalty2DCondition();


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
     * Calculation of the Integration Weight
     */
    double& CalculateIntegrationWeight(double& rIntegrationWeight) override;



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

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, PointRigidContactPenalty3DCondition )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, PointRigidContactPenalty3DCondition )
    }


}; // Class PointRigidContactPenalty2DCondition

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
				    PointRigidContactPenalty2DCondition& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
				    const PointRigidContactPenalty2DCondition& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

}  // namespace Kratos.

#endif // KRATOS_POINT_RIGID_CONTACT_PENALTY_2D_CONDITION_H_INCLUDED  defined
