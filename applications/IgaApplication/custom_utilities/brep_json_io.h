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
    /// Pointer definition of KratosNurbsTestcaseApplication
    KRATOS_CLASS_POINTER_DEFINITION(BrepJsonIO);

    void WriteIntegrationDomainJson(
        ModelPart& rModelPart, 
        const std::string& rOutputFileName);

    void ExportNurbsGeometry(std::vector<BrepModel>  m_brep_model_vector);

    std::vector<BrepModel> ImportNurbsBrepGeometry(
        ModelPart& rModelPart,
        Parameters rNurbsBrepGeometryJson);

    /// Constructor.
    BrepJsonIO(const int EchoLevel = 0)
        :mEchoLevel(EchoLevel)
    {};

    /// Destructor.
    virtual ~BrepJsonIO() {};

  private:

      void ImportBrepEdges(
          const Parameters& rEdges,
          std::vector<BrepEdge>& rEdgesVector,
          ModelPart& rModelPart);

      void ImportBrepVertices(
          const Parameters& rVertices,
          std::vector<BrepVertex>& rVerticesVector,
          ModelPart& rModelPart);

      void ImportTrimmingCurve(
          const Parameters& rTrimmingCurve,
          std::vector<BrepTrimmingCurve>& rrTrimmingCurves);

      int mEchoLevel;
  }; // Class BrepJsonIO
}  // namespace Kratos.

#endif // KRATOS_BREP_JSON_IO_H_INCLUDED  defined