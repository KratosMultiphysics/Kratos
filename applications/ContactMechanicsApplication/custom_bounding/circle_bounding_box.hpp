//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

#if !defined(KRATOS_CIRCLE_BOUNDING_BOX_H_INCLUDED )
#define  KRATOS_CIRCLE_BOUNDING_BOX_H_INCLUDED

// External includes

// System includes

// Project includes
#include "custom_bounding/sphere_bounding_box.hpp"


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

    This Box represents a 2D wall composed by a circle
    
    A convexity parameter is given to determine which side of each nose is considered 
    the internal or external boundary

    This bounding box is essentially used for rigid wall contact purposes
*/

class KRATOS_API(CONTACT_MECHANICS_APPLICATION) CircleBoundingBox
  : public SphereBoundingBox
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CircleBoundingBox
    KRATOS_CLASS_POINTER_DEFINITION( CircleBoundingBox );

    //typedef bounded_vector<double, 3>                     PointType;
    typedef array_1d<double, 3>                             PointType;
    typedef ModelPart::NodeType                              NodeType;
    typedef ModelPart::NodesContainerType          NodesContainerType;
    typedef NodesContainerType::Pointer     NodesContainerTypePointer;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CircleBoundingBox() : SphereBoundingBox()
    {
      KRATOS_TRY

      std::cout<< "Calling Rigid Circle Wall BBX empty constructor" <<std::endl;

      KRATOS_CATCH("")
    }


    //**************************************************************************
    //**************************************************************************

    CircleBoundingBox(Parameters CustomParameters) 
      : SphereBoundingBox(CustomParameters)
    {
      KRATOS_TRY

      KRATOS_CATCH("") 
    }

    //**************************************************************************
    //**************************************************************************


    // General Wall constructor
    CircleBoundingBox(PointType Center,
		      double Radius,
		      PointType Velocity,
		      int Convexity) 
      : SphereBoundingBox(Center, Radius, Velocity, Convexity)
    {
      KRATOS_TRY

      KRATOS_CATCH("") 
    }

    

    //**************************************************************************
    //**************************************************************************

    /// Assignment operator.
    CircleBoundingBox& operator=(CircleBoundingBox const& rOther)
    {
      KRATOS_TRY

      SphereBoundingBox::operator=(rOther);
      return *this;

      KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************

    /// Copy constructor.
    CircleBoundingBox(CircleBoundingBox const& rOther) 
      :SphereBoundingBox(rOther)
    {
    }


    //**************************************************************************
    //**************************************************************************

    /// Destructor.
    virtual ~CircleBoundingBox() {};


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
  

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
        return "CircleBoundingBox";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << this->mBox.UpperPoint << " , " << this->mBox.LowerPoint;
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


    ///@}


}; // Class CircleBoundingBox

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_CIRCLE_BOUNDING_BOX_H_INCLUDED  defined 


