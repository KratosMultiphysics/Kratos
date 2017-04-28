#if !defined(KRATOS_BREP_MODEL_H_INCLUDED )
#define  KRATOS_BREP_MODEL_H_INCLUDED


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
#include "BrepFace.h"
#include "BrepEdge.h"

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
  class BrepModel : public IndexedObject, public Flags
  {
  public:
    ///@name Type Definitions
    ///@{

    typedef std::vector<BrepFace> FacesVector;
    typedef std::vector<BrepEdge> EdgesVector;
    
    /// Pointer definition of KratosNurbsTestcaseApplication
    KRATOS_CLASS_POINTER_DEFINITION(BrepModel);

    ///@}
    ///@name Life Cycle 
    ///@{ 
    FacesVector& GetFaceVector();
    EdgesVector& GetEdgeVector();


    //TODO: you need to give reading access to your internals through the Calculate function
    /// Constructor.
    BrepModel(unsigned int& brep_id, 
      FacesVector& faces, 
      EdgesVector& edges);

    /// Destructor.
    virtual ~BrepModel();

    /// Copy constructor.
    //BrepModel(BrepModel const& rOther);

    /// Assignment operator.
    //BrepModel& operator=(BrepModel const& rOther);
    ///@} 
  protected:

  private:

    ///@name Private methods
    ///@{ 

    ///@} 
    ///@name Member Variables
    ///@{ 
    FacesVector m_brep_faces;
    EdgesVector m_brep_edges;
    ///@}    

  }; // Class BrepModel 

}  // namespace Kratos.

#endif // KRATOS_BREP_MODEL_H_INCLUDED  defined