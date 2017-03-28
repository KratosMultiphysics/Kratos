#if !defined(KRATOS_BREP_BOUNDARY_LOOP_H_INCLUDED )
#define  KRATOS_BREP_BOUNDARY_LOOP_H_INCLUDED



// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>
#include <cmath>
#include <math.h>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "BrepTrimmingCurve.h"


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
  class BrepBoundaryLoop
  {
  public:
    ///@name Type Definitions
    ///@{

    typedef std::vector<BrepTrimmingCurve> BrepTrimmingCurveVector;
    
    /// Pointer definition of KratosNurbsBrepApplication
    KRATOS_CLASS_POINTER_DEFINITION(BrepBoundaryLoop);

    ///@}
    ///@name Life Cycle 
    ///@{ 

    //TrimmingCurveVector& GetTrimmingCurves();
    bool& IsOuterLoop();
    std::vector<array_1d<double, 2>> GetBoundaryPolygon();

    /// Constructor.
    BrepBoundaryLoop(BrepTrimmingCurveVector& brep_trimming_curves, bool is_outer_loop);

    /// Destructor.
    virtual ~BrepBoundaryLoop();

    /// Copy constructor.
    //BrepBoundaryLoop(BrepBoundaryLoop const& rOther);

    /// Assignment operator.
    //BrepBoundaryLoop& operator=(BrepBoundaryLoop const& rOther);
    ///@} 
  protected:

  private:
    ///@name Member Variables
    ///@{ 

    BrepTrimmingCurveVector m_brep_trimming_curves;
    bool m_is_outer_loop;

  
    ///@name Un accessible methods 
    ///@{ 



    ///@}    

  }; // Class BrepBoundaryLoop 

}  // namespace Kratos.

#endif // KRATOS_BREP_BOUNDARY_LOOP_H_INCLUDED  defined