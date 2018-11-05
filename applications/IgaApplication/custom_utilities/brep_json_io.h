#if !defined(KRATOS_BREP_JSON_IO_H_INCLUDED )
#define  KRATOS_BREP_JSON_IO_H_INCLUDED


// System includes

// External includes

// Project includes
#include "brep_topology/brep_model.h"
#include "brep_topology/brep_trimming_curve.h"
#include "brep_topology/brep_boundary_loop.h"
#include "brep_topology/brep_face.h"
#include "brep_topology/brep_edge.h"

#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "includes/node.h"

namespace Kratos
{

  ///@name Kratos Classes
  ///@{
  /// Short class definition.
  /** Gives IO capabilities for Nurbs based Brep models in the JSON format defined in 
  https://amses-journal.springeropen.com/articles/10.1186/s40323-018-0109-4.
  */
    class BrepJsonIO
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosNurbsTestcaseApplication
    KRATOS_CLASS_POINTER_DEFINITION(BrepJsonIO);

    ///@}
    ///@name Life Cycle
    ///@{
    void WriteIntegrationDomainJson(
        ModelPart& rModelPart, 
        const std::string& rOutputFileName);

    std::vector<BrepModel> ImportNurbsBrepGeometry(
        ModelPart& model_part, 
        Parameters& rNurbsBrepGeometryJson);

    /// Constructor.
    BrepJsonIO() {};

    /// Destructor.
    virtual ~BrepJsonIO() {};

    ///@}
  protected:

  private:

  }; // Class BrepJsonIO
  ///@}
  ///@name Type Definitions
  ///@{
  ///@}
  ///@name Input and output
  ///@{
  ///@}

}  // namespace Kratos.

#endif // KRATOS_BREP_JSON_IO_H_INCLUDED  defined