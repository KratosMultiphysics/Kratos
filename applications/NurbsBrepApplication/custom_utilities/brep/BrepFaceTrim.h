#if !defined(KRATOS_BREP_FACE_TRIM_H_INCLUDED )
#define  KRATOS_BREP_FACE_TRIM_H_INCLUDED


// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>
#include <cmath>
#include <math.h>


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
  class BrepFaceTrim
  {
  public:
    ///@name Type Definitions
    ///@{

    //typedef std::vector<TrimmingCurve> TrimmingCurveVector;
    
    /// Pointer definition of KratosNurbsBrepApplication
    //KRATOS_CLASS_POINTER_DEFINITION(BrepFaceTrim);

    ///@}
    ///@name Life Cycle 
    ///@{ 

    /// Constructor.
    BrepFaceTrim(unsigned int face_id, 
      unsigned int trim_index, 
      bool relative_direction);

    /// Destructor.
    virtual ~BrepFaceTrim();

    /// Copy constructor.
    //BrepFaceTrim(BrepFaceTrim const& rOther);

    /// Assignment operator.
    //BrepFaceTrim& operator=(BrepFaceTrim const& rOther);
    ///@} 
  protected:

  private:
    ///@name Member Variables
    ///@{ 
    unsigned int m_face_id;
    unsigned int m_trim_index;
    bool m_relative_direction;
    ///@}    
     
    ///@name Un accessible methods 
    ///@{ 



    ///@}    

  }; // Class BrepFaceTrim 

}  // namespace Kratos.

#endif // KRATOS_BREP_FACE_TRIM_H_INCLUDED  defined