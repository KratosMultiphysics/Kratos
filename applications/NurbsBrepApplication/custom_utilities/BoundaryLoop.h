#if !defined(KRATOS_BREP_MODEL_GEOMETRY_READER_H_INCLUDED )
#define  KRATOS_BREP_MODEL_GEOMETRY_READER_H_INCLUDED



// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>
#include <cmath>
#include <math.h>

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------
//#include <boost/python.hpp>
//#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/vector.hpp>
//#include <boost/numeric/ublas/io.hpp>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "TrimmingCurve.h"


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
  class BoundaryLoop
  {
  public:
    ///@name Type Definitions
    ///@{

    typedef std::vector<TrimmingCurve> TrimmingCurveVector;
    
    /// Pointer definition of KratosNurbsTestcaseApplication
    KRATOS_CLASS_POINTER_DEFINITION(BoundaryLoop);

    ///@}
    ///@name Life Cycle 
    ///@{ 

    //TrimmingCurveVector& GetTrimmingCurves();
    bool& IsOuterLoop();
    std::vector<array_1d<double, 2>> BoundaryLoop::GetBoundaryPolygon();

    /// Constructor.
    BoundaryLoop(TrimmingCurveVector& trimming_curves, bool is_outer_loop);

    /// Destructor.
    virtual ~BoundaryLoop();

    /// Copy constructor.
    //BoundaryLoop(BoundaryLoop const& rOther);

    /// Assignment operator.
    //BoundaryLoop& operator=(BoundaryLoop const& rOther);
    ///@} 
  protected:

  private:
    ///@name Member Variables
    ///@{ 

    TrimmingCurveVector m_trimming_curves;
    bool m_is_outer_loop;

  
    ///@name Un accessible methods 
    ///@{ 



    ///@}    

  }; // Class BoundaryLoop 

}  // namespace Kratos.

#endif // KRATOS_BOUNDARY_LOOP_APPLICATION_H_INCLUDED  defined