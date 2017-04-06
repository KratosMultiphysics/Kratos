#if !defined(KRATOS_KNOT_SPAN_1D_H_INCLUDED )
#define  KRATOS_KNOT_SPAN_1D_H_INCLUDED


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
  class KnotSpan1d : public IndexedObject, public Flags
  {
  public:
    ///@name Type Definitions
    ///@{

    
    /// Pointer definition of KratosNurbsTestcaseApplication
    //KRATOS_CLASS_POINTER_DEFINITION(KnotSpan1d);

    ///@}
    ///@name Life Cycle 
    ///@{ 
    std::vector<array_1d<double, 3>> KnotSpan1d::getIntegrationPointsInFullGaussianDomain();
    std::vector<array_1d<double, 3>> KnotSpan1d::getIntegrationPointsInParameterDomain();


    //TODO: you need to give reading access to your internals through the Calculate function
    /// Constructor.
    KnotSpan1d(unsigned int knot_span_1d_id,
      int p,
      Vector parameter_span_u);

    /// Destructor.
    virtual ~KnotSpan1d();

    /// Copy constructor.
    //KnotSpan1d(KnotSpan1d const& rOther);

    /// Assignment operator.
    //KnotSpan1d& operator=(KnotSpan1d const& rOther);
    ///@} 
  protected:

  private:

    ///@name Private methods
    ///@{ 

    ///@} 
    ///@name Member Variables
    ///@{ 
    int m_p;
    Vector m_parameter_span_u;

    ///@}    

  }; // Class KnotSpan1d 

}  // namespace Kratos.

#endif // KRATOS_KNOT_SPAN_1D_H_INCLUDED  defined