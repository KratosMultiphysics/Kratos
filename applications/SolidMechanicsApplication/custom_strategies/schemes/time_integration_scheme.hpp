//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            November 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_TIME_INTEGRATION_SCHEME )
#define  KRATOS_TIME_INTEGRATION_SCHEME

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/variables.h"

namespace Kratos
{
  ///@addtogroup ConstitutiveModelsApplication
  ///@{

  
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
   * This class performs predict and update of dofs variables, their time derivatives and time integrals      
   */
  class KRATOS_API(SOLID_MECHANICS_APPLICATION) TimeIntegrationScheme
  {
  public:
 
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( TimeIntegrationScheme );

    ///@}
    ///@name Life Cycle
    ///@{

    
    /// Default Constructor.
    TimeIntegrationScheme(double rAlpham = 0.0) {}

    /// Copy Constructor.
    TimeIntegrationScheme(TimeIntegrationScheme& rOther) {}

    /// Clone
    TimeIntegrationScheme::Pointer Clone()
    {
      return TimeIntegrationScheme::Pointer( new TimeIntegrationScheme(*this) );
    }

    /// Destructor.
    ~TimeIntegrationScheme(){}

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
        std::stringstream buffer;
        buffer << "TimeIntegrationScheme";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "TimeIntegrationScheme";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
      rOStream << "TimeIntegrationScheme Data";     
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
    ///@name Serialization
    ///@{
  
    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{
  
    ///@}
  
  }; // Class TimeIntegrationScheme
  
  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
                                    TimeIntegrationScheme& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
                                    const TimeIntegrationScheme& rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream <<" : " << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }
  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_TIME_INTEGRATION_SCHEME defined
