#if !defined(KRATOS_EMBEDDED_IGA_TESSELLATION_H_INCLUDED )
#define  KRATOS_EMBEDDED_IGA_TESSELLATION_H_INCLUDED

extern "C" 
{
    #ifdef SINGLE
        #define REAL float
    #else /* not SINGLE */
        #define REAL double
    #endif /* not SINGLE */
    void triangulate(char *, struct triangulateio *, struct triangulateio *,struct triangulateio *);    
}

// System includes

// External includes
#include "custom_utilities/anurbs.h"

// Project includes
#include "iga_application_variables.h"
#include "custom_utilities/nurbs_brep_modeler.h"

namespace Kratos
{
class EmbeddedIgaTessellation
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosNurbsTestcaseApplication
    KRATOS_CLASS_POINTER_DEFINITION(EmbeddedIgaTessellation);

    ///@}
    ///@name functions
    ///@{
    
    static void CreateTessellation1D(
        const double mTessellationTolerance,
        const BrepEdge& rCurveGeometry,
        std::vector<array_1d<double, 3>>& rPolygon);

    static void CreateTessellation2D(
        const double mTessellationTolerance,
        const BrepFace& rFaceGeometry,
        std::vector<std::vector<array_1d<double, 2>>>& rOuterPolygon,
        std::vector<std::vector<array_1d<double, 2>>>& rInnerPolygon);
    
    ///@}
    ///@name Life Cycle
    ///@{
    /// Constructor.
    EmbeddedIgaTessellation();

    /// Destructor.
    virtual ~EmbeddedIgaTessellation()
    {};

    ///@}
protected:

private:
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

};

}  // namespace Kratos.
#endif // KRATOS_EMBEDDED_IGA_TESSELLATION_H_INCLUDED defined


