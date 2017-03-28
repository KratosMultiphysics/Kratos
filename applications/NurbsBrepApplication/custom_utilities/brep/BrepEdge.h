#if !defined(KRATOS_BREP_EDGE_H_INCLUDED )
#define  KRATOS_BREP_EDGE_H_INCLUDED


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
#include "BrepFaceTrim.h"
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
  class BrepEdge : public IndexedObject, public Flags
  {
  public:
    ///@name Type Definitions
    ///@{

    typedef std::vector<BrepFaceTrim> BrepFaceTrimVector;
    typedef std::vector<Vector> ParameterVector;
    
    /// Pointer definition of KratosNurbsTestcaseApplication
    //KRATOS_CLASS_POINTER_DEFINITION(BrepEdge);

    ///@}
    ///@name Life Cycle 
    ///@{ 



    //TODO: you need to give reading access to your internals through the Calculate function
    /// Constructor.
    BrepEdge(unsigned int edge_id,
      ParameterVector& boundary_vertices,
      BrepFaceTrimVector& brep_face_trims_vector);

    /// Destructor.
    virtual ~BrepEdge();

    /// Copy constructor.
    //BrepEdge(BrepEdge const& rOther);

    /// Assignment operator.
    //BrepEdge& operator=(BrepEdge const& rOther);
    ///@} 
  protected:

  private:

    ///@name Private methods
    ///@{ 

    ///@} 
    ///@name Member Variables
    ///@{ 

    ParameterVector m_boundary_vertices;
    BrepFaceTrimVector m_brep_face_trims_vector;

    ///@}    

  }; // Class BrepEdge 

}  // namespace Kratos.

#endif // KRATOS_BREP_EDGE_H_INCLUDED  defined