//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_POINT_RIGID_CONTACT_PENALTY_2D_CONDITION_H_INCLUDED )
#define  KRATOS_POINT_RIGID_CONTACT_PENALTY_2D_CONDITION_H_INCLUDED



// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "utilities/math_utils.h"

#include "custom_modelers/rigid_tool_bounding_box.hpp"

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
class PointRigidContactPenalty2DCondition
    : public Condition
{
public:
    ///@name Type Definitions
    ///@{
  
    ///Tensor order 1 definition
    typedef Vector VectorType;
  

protected:


   typedef struct
    {
        VectorType Normal;        //normal direction
        VectorType Tangent;       //tangent direction
	   
    } SurfaceVector;


    typedef struct
    {
        double Normal;        //normal component 
        double Tangent;       //tangent component
        	   
    } SurfaceScalar;


    typedef struct
    {
      Flags           Options;               //calculation options
      
      //Geometrical gaps:
      SurfaceScalar   Gap;                   //normal and tangential gap

      //Friction:
      double          FrictionCoefficient;   //total friction coeffitient mu

      SurfaceScalar   Penalty;               //Penalty Parameter normal and tangent

      //Geometric variables
      SurfaceVector   Surface;               //normal and tangent vector to the surface
      
    } ContactVariables;


public:

    /// Counted pointer of PointRigidContactPenalty2DCondition
    KRATOS_CLASS_POINTER_DEFINITION(PointRigidContactPenalty2DCondition);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
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

    Condition::Pointer Create(IndexType NewId, NodesArrayType const&
                              ThisNodes,  PropertiesType::Pointer pProperties) const;

  
    void InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo);
  
    void FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo);
  
  
    void SetRigidWall(SpatialBoundingBox::Pointer pRigidWall);


    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType&
                              rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo&
                                rCurrentProcessInfo);
    
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo&
                          rCurrentProcessInfo);

    void GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo&
                    CurrentProcessInfo);

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

    // A protected default constructor necessary for serialization
    PointRigidContactPenalty2DCondition() {};

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

    SpatialBoundingBox::Pointer mpRigidWall;
  

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition );
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition );
    }

    ///@}
    ///@name Private Operators
    ///@{

    void CalculateContactFactors(ContactVariables &rContact);

    inline MatrixType outer_prod_2(const array_1d<double, 3>& a, const array_1d<double, 3>& b);


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
    //PointRigidContactPenalty2DCondition& operator=(const PointRigidContactPenalty2DCondition& rOther);

    /// Copy constructor.
    //PointRigidContactPenalty2DCondition(const PointRigidContactPenalty2DCondition& rOther);


    ///@}

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



