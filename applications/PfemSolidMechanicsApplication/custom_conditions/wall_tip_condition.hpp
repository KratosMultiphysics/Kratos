//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Last modified by:    $Author:                JMCarbonell $
//   Date:                $Date:                    July 2013 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_WALL_TIP_CONDITION_H_INCLUDED )
#define  KRATOS_WALL_TIP_CONDITION_H_INCLUDED



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

class WallTipCondition
    : public Condition
{
public:

    ///@name Type Definitions
    ///@{

    /// Counted pointer of WallTipCondition
    KRATOS_CLASS_POINTER_DEFINITION( WallTipCondition );


    struct WallTipVariables
    {
      bool    wall_tip_set;
      double  Radius;
      array_1d<double,3> Center;
      array_1d<double,3> Velocity;
      array_1d<double,3> Position;
    };


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    WallTipCondition( IndexType NewId, GeometryType::Pointer pGeometry );
    WallTipCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties );

    /// Destructor.
    virtual ~WallTipCondition();


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Condition::Pointer Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const;

    void Initialize ();

    void InitializeSolutionStep ( ProcessInfo& CurrentProcessInfo );
    void FinalizeSolutionStep   ( ProcessInfo& CurrentProcessInfo );    


    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    virtual int Check( const ProcessInfo& rCurrentProcessInfo ){ return 1; };
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
    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Tool Tip 2D Condition #" << Id();
    }


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

    WallTipVariables mWallTip;


    ///@}
    ///@name Private Operators
    ///@{
    void CalculateConditionalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
                                   ProcessInfo& rCurrentProcessInfo,
                                   bool CalculateStiffnessMatrixFlag,
                                   bool CalculateResidualVectorFlag );

    void CalculateAndSubKp(Matrix& rK,
                           const Matrix& rDN_De,
                           const Vector& rN,
                           double rPressure,
                           double rIntegrationWeight
                          );


    void CalculateAndAddFacePressure (Vector& rF,
                                      const Vector& rN,
                                      Vector& rNormal,
                                      double rPressure,
                                      double rIntegrationWeight );

    void CalculateAndAddFaceForce(Vector& rF,
                                  const Vector& rN,
                                  Vector& rForce,
                                  double rIntegrationWeight );

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

    // A private default constructor necessary for serialization
    WallTipCondition() {};

    virtual void save( Serializer& rSerializer ) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition );
    }

    virtual void load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition );
	mWallTip.wall_tip_set =false;
    }


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    //WallTipCondition& operator=(const WallTipCondition& rOther);

    /// Copy constructor.
    //WallTipCondition(const WallTipCondition& rOther);


    ///@}

}; // Class WallTipCondition

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
        WallTipCondition& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
        const WallTipCondition& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

}  // namespace Kratos.

#endif // KRATOS_WALL_TIP_CONDITION_H_INCLUDED  defined 


