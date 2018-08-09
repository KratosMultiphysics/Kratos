#if !defined(KRATOS_KNOT_SPAN_2D_N_INTEGRATE_H_INCLUDED )
#define  KRATOS_KNOT_SPAN_2D_N_INTEGRATE_H_INCLUDED


// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>
#include <cmath>
#include <math.h>
#include <vector>

#include <ios>
#include <fstream>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "../integration_utilities.h"
#include "../Polygon.h"

//#include "../../kratos/includes/node.h"
//
#include "nurbs_brep_application.h"
#include "nurbs_brep_application_variables.h"

#include "KnotSpan2d.h"

// ==============================================================================

namespace Kratos
{

  ///@name Kratos Globals
  ///@{ 
  ///@} 
  ///@name Type Definitions
  ///@{ 
  //typedef Node<3> NodeType;
  //typedef std::vector<NodeType::Pointer> NodeVector;
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
  class KnotSpan2dNIntegrate : KnotSpan2d
  {
  public:
    ///@name Type Definitions
    ///@{

    
    /// Pointer definition of KratosNurbsTestcaseApplication
    //KRATOS_CLASS_POINTER_DEFINITION(KnotSpan2dNIntegrate);

    ///@}
    ///@name Life Cycle 
    ///@{ 
    std::vector<array_1d<double, 3>> getIntegrationPointsInFullGaussianDomain() override;
    std::vector<array_1d<double, 3>> getIntegrationPointsInParameterDomain() override;

    std::vector<array_1d<double, 3>> pointsByQuad(const Matrix &quad);
    std::vector<array_1d<double, 3>> pointsByTriangle(const Matrix &triangle);

    //bool triangulate();
    //void CALLBACK onTessBeginData(GLenum type, void *polygonData);
    //void CALLBACK onTessVertexData(void *vertexData, void *polygonData);
    //void CALLBACK onTessEndData(void *polygonData);

    //std::vector<Matrix> tessellate(const std::vector<Vector> &paths);

    //TODO: you need to give reading access to your internals through the Calculate function
    /// Constructor.
    KnotSpan2dNIntegrate(
      unsigned int knot_span_2d_id,
      bool is_untrimmed,
      int p, int q,
      Vector parameter_span_u,
      Vector parameter_span_v,
      Polygon polygon);

    /// Destructor.
    virtual ~KnotSpan2dNIntegrate();

    /// Copy constructor.
    //KnotSpan2dNIntegrate(KnotSpan2dNIntegrate const& rOther);

    /// Assignment operator.
    //KnotSpan2dNIntegrate& operator=(KnotSpan2dNIntegrate const& rOther);
    ///@} 
  protected:

  private:

    ///@name Private methods
    ///@{ 

    ///@} 
    ///@name Member Variables
    ///@{ 
    Polygon m_polygon;
    //bool m_is_untrimmed;
    //int m_p;
    //int m_q;
    //Vector m_parameter_span_u;
    //Vector m_parameter_span_v;

    ///@}    

  }; // Class KnotSpan2dNIntegrate 

}  // namespace Kratos.

#endif // KRATOS_KNOT_SPAN_2D_N_INTEGRATE_H_INCLUDED  defined