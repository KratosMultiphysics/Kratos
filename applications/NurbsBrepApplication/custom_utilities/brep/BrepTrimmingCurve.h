#if !defined(KRATOS_BREP_TRIMMING_CURVE_H_INCLUDED )
#define  KRATOS_BREP_TRIMMING_CURVE_H_INCLUDED


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
#include "nurbs_brep_application.h"
#include "nurbs_brep_application_variables.h"

#include "../nurbs_utilities.h"


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
  class BrepTrimmingCurve
  {
  public:
    ///@name Type Definitions
    ///@{

    typedef std::vector<array_1d<double, 4>> ControlPointVector;
    
    /// Pointer definition of KratosNurbsBrepApplication
    KRATOS_CLASS_POINTER_DEFINITION(BrepTrimmingCurve);

    ///@}
    ///@name Life Cycle 
    ///@{ 
    std::vector<array_1d<double, 2>> CreatePolygon(unsigned int& number_polygon_points);
    void EvaluateCurvePoint(Point<3>& rCurvePoint, double parameter_u);
    unsigned int& GetIndex();

    //TODO: you need to give reading access to your internals through the Calculate function
    /// Constructor.
    //TODO: pass by reference not by value
    //TODO: why control points have size 4??? pass them as kratos nodes
    BrepTrimmingCurve(unsigned int trim_index, bool curve_direction, Vector& knot_vector_u,
      unsigned int p, ControlPointVector& control_points,
      Vector& active_range);

    /// Destructor.
    virtual ~BrepTrimmingCurve();

    /// Copy constructor.
    //BrepTrimmingCurve(BrepTrimmingCurve const& rOther);

    /// Assignment operator.
    //BrepTrimmingCurve& operator=(BrepTrimmingCurve const& rOther);
    ///@} 
  protected:

  private:
    ///@name Member Variables
    ///@{ 
    unsigned int m_trim_index;
    bool m_curve_direction;
    Vector m_knot_vector_u;
    unsigned int m_p;
    ControlPointVector m_control_points;
    Vector m_active_range;
    ///@}    
     
    ///@name Un accessible methods 
    ///@{ 



    ///@}    

  }; // Class BrepTrimmingCurve 

}  // namespace Kratos.

#endif // KRATOS_BREP_TRIMMING_CURVE_H_INCLUDED  defined