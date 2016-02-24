//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2013 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_WALL_CONDITION_H_INCLUDED )
#define  KRATOS_WALL_CONDITION_H_INCLUDED



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

class WallCondition
    : public Condition
{
public:

    ///@name Type Definitions
    ///@{

    /// Counted pointer of WallCondition
    KRATOS_CLASS_POINTER_DEFINITION( WallCondition );


    struct WallVariables
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
    WallCondition( IndexType NewId, GeometryType::Pointer pGeometry );
    WallCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties );

    /// Destructor.
    virtual ~WallCondition();


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Condition::Pointer Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const;

    Condition::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const;


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
    // virtual String Info() const;


    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Wall Condition #" << Id();
    }


    /// Print object's data.
    //  virtual void PrintData(std::ostream& rOStream) const;


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

    WallVariables mWallVariables;


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

    // A private default constructor necessary for serialization
    WallCondition() {};

    virtual void save( Serializer& rSerializer ) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition )
    }

    virtual void load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition )
	mWallVariables.wall_tip_set =false;
    }


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    //WallCondition& operator=(const WallCondition& rOther);

    /// Copy constructor.
    //WallCondition(const WallCondition& rOther);


    ///@}

}; // Class WallCondition

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
        WallCondition& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
        const WallCondition& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

}  // namespace Kratos.

#endif // KRATOS_WALL_CONDITION_H_INCLUDED  defined 


