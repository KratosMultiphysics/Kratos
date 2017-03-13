#if !defined(KRATOS_EDGE_H_INCLUDED )
#define  KRATOS_EDGE_H_INCLUDED


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
#include "FaceTrim.h"
#include "../../kratos/includes/node.h"

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
  class Edge : public IndexedObject, public Flags
  {
  public:
    ///@name Type Definitions
    ///@{

    typedef std::vector<FaceTrim> FaceTrimVector;
    typedef std::vector<Vector> ParameterVector;
    
    /// Pointer definition of KratosNurbsTestcaseApplication
    //KRATOS_CLASS_POINTER_DEFINITION(Edge);

    ///@}
    ///@name Life Cycle 
    ///@{ 



    //TODO: you need to give reading access to your internals through the Calculate function
    /// Constructor.
    Edge(unsigned int edge_id,
      ParameterVector& boundary_vertices,
      FaceTrimVector& face_trims_vector);

    /// Destructor.
    virtual ~Edge();

    /// Copy constructor.
    //Edge(Edge const& rOther);

    /// Assignment operator.
    //Edge& operator=(Edge const& rOther);
    ///@} 
  protected:

  private:

    ///@name Private methods
    ///@{ 

    ///@} 
    ///@name Member Variables
    ///@{ 

    ParameterVector m_boundary_vertices;
    FaceTrimVector m_face_trims_vector;

    ///@}    

  }; // Class Edge 

}  // namespace Kratos.

#endif // KRATOS_EDGE_H_INCLUDED  defined