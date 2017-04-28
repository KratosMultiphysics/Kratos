#if !defined(KRATOS_KNOT_SPAN_2D_H_INCLUDED )
#define  KRATOS_KNOT_SPAN_2D_H_INCLUDED


// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>
#include <cmath>
#include <math.h>
#include <vector>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "../integration_utilities.h"
//#include "../../kratos/includes/node.h"
//
#include "nurbs_brep_application.h"
#include "nurbs_brep_application_variables.h"

// ==============================================================================

namespace Kratos
{

  ///@name Kratos Globals
  ///@{ 
  ///@} 
  ///@name Type Definitions
  ///@{ 
  typedef Node<3> NodeType;
  typedef std::vector<NodeType::Pointer> NodeVector;
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
  class KnotSpan2d : public IndexedObject, public Flags
  {
  public:
    ///@name Type Definitions
    ///@{

    
    /// Pointer definition of KratosNurbsTestcaseApplication
    //KRATOS_CLASS_POINTER_DEFINITION(KnotSpan2d);

    ///@}
    ///@name Life Cycle 
    ///@{ 
    std::vector<array_1d<double, 3>> KnotSpan2d::getIntegrationPointsInFullGaussianDomain();
    std::vector<array_1d<double, 3>> KnotSpan2d::getIntegrationPointsInParameterDomain();


    //TODO: you need to give reading access to your internals through the Calculate function
    /// Constructor.
    KnotSpan2d(unsigned int knot_span_2d_id,
      bool is_untrimmed,
      int p, int q,
      Vector parameter_span_u,
      Vector parameter_span_v);

    /// Destructor.
    virtual ~KnotSpan2d();

    /// Copy constructor.
    //KnotSpan2d(KnotSpan2d const& rOther);

    /// Assignment operator.
    //KnotSpan2d& operator=(KnotSpan2d const& rOther);
    ///@} 
  protected:

  private:

    ///@name Private methods
    ///@{ 

    ///@} 
    ///@name Member Variables
    ///@{ 
    bool m_is_untrimmed;
    int m_p;
    int m_q;
    Vector m_parameter_span_u;
    Vector m_parameter_span_v;

    ///@}    

  }; // Class KnotSpan2d 

}  // namespace Kratos.

#endif // KRATOS_KNOT_SPAN_2D_H_INCLUDED  defined