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
#include "triangle.h"  

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
    
    std::vector<std::vector<double>> CreateTessellation(
        const BrepFace& rFaceGeometry); 
    
    void InitTriangulationDataStructure(triangulateio& tr)
    {
        tr.pointlist                  = (REAL*) NULL;
        tr.pointattributelist         = (REAL*) NULL;
        tr.pointmarkerlist            = (int*) NULL;
        tr.numberofpoints             = 0;
        tr.numberofpointattributes    = 0;
        tr.trianglelist               = (int*) NULL;
        tr.triangleattributelist      = (REAL*) NULL;
        tr.trianglearealist           = (REAL*) NULL;
        tr.neighborlist               = (int*) NULL;
        tr.numberoftriangles          = 0;
        tr.numberofcorners            = 3;
        tr.numberoftriangleattributes = 0;
        tr.segmentlist                = (int*) NULL;
        tr.segmentmarkerlist          = (int*) NULL;
        tr.numberofsegments           = 0;
        tr.holelist                   = (REAL*) NULL;
        tr.numberofholes              = 0;
        tr.regionlist                 = (REAL*) NULL;
        tr.numberofregions            = 0;
        tr.edgelist                   = (int*) NULL;
        tr.edgemarkerlist             = (int*) NULL;
        tr.normlist                   = (REAL*) NULL;
        tr.numberofedges              = 0;
    };  

    
    
    
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


