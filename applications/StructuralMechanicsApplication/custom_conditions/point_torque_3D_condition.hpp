//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_POINT_TORQUE_3D_CONDITION_H_INCLUDED )
#define  KRATOS_POINT_TORQUE_3D_CONDITION_H_INCLUDED


// System includes

// External includes

// Project includes
#include "includes/condition.h"
#include "utilities/math_utils.h"

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
class PointTorque3DCondition
    : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of PointTorque3DCondition
    KRATOS_CLASS_POINTER_DEFINITION( PointTorque3DCondition );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PointTorque3DCondition(IndexType NewId, GeometryType::Pointer pGeometry);
    PointTorque3DCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

    ///Copy constructor
    PointTorque3DCondition( PointTorque3DCondition const& rOther);

    /// Destructor.
    virtual ~PointTorque3DCondition();


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    //************* STARTING - ENDING  METHODS


    /**
     * Called at the beginning of each solution step
     */
    void Initialize() override;


    /**
     * Called at the beginning of each iteration
     */
    void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;


    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
    //virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

    void GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo) override;


    /**
     * Get on rVariable a double Value
     */
    void GetValueOnIntegrationPoints( const Variable<double>& rVariable,
                                      std::vector<double>& rValues,
                                      const ProcessInfo& rCurrentProcessInfo ) override;

    /**
     * Calculate a double Variable
     */
    void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                                      std::vector<double>& rOutput,
                                      const ProcessInfo& rCurrentProcessInfo) override;


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

    friend class Serializer;

    // A private default constructor necessary for serialization
    PointTorque3DCondition() {};

    virtual void save(Serializer& rSerializer) const override;

    virtual void load(Serializer& rSerializer)override;

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
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    //PointTorque3DCondition& operator=(const PointTorque3DCondition& rOther);

    /// Copy constructor.
    //PointTorque3DCondition(const PointTorque3DCondition& rOther);


    ///@}

}; // Class PointTorque3DCondition

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
				    PointTorque3DCondition& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
				    const PointTorque3DCondition& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

}  // namespace Kratos.

#endif // KRATOS_POINT_TORQUE_3D_CONDITION_H_INCLUDED  defined 


