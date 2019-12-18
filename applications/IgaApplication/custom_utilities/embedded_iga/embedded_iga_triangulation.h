#if !defined(KRATOS_EMBEDDED_IGA_TRIANGULATION_H_INCLUDED )
#define  KRATOS_EMBEDDED_IGA_TRIANGULATION_H_INCLUDED


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
#include "triangle.h"  

// Project includes
#include "iga_application_variables.h"
#include "custom_utilities/nurbs_brep_modeler.h"
#include "custom_utilities/embedded_iga/embedded_iga_error_estimation.h"
#include "custom_external_libraries/polylabel/include/mapbox/polylabel.hpp"

namespace Kratos
{
class EmbeddedIgaTriangulation
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosNurbsTestcaseApplication
    KRATOS_CLASS_POINTER_DEFINITION(EmbeddedIgaTriangulation);

    ///@}
    ///@name functions
    ///@{
    
    void CreateTriangulation(
        const double mTriangulationTolerance,
        const double mInitialTriangleArea,
        const int mMaxTriangulationIterations,
        const int mEchoLevel,
        const BrepFace& rFaceGeometry,
        const std::vector<std::vector<array_1d<double,2>>>& rOuterPolygon,
        const std::vector<std::vector<array_1d<double,2>>>& rInnerPolygon,
        std::vector<Matrix>& rTriangulation_xyz);
    
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

    void CleanTriangulationDataStructure( triangulateio& tr )
    {
        if(tr.pointlist != NULL) 
        {
            free(tr.pointlist);
            tr.pointlist = nullptr; 
        }
        if(tr.pointattributelist != NULL) 
        {
            free(tr.pointattributelist);
            tr.pointattributelist = nullptr;
        }
        if(tr.pointmarkerlist != NULL) 
        {
            free(tr.pointmarkerlist);
            tr.pointmarkerlist = nullptr; 
        }
        if(tr.trianglelist != NULL)
        {
            free(tr.trianglelist);
            tr.trianglelist = nullptr; 
        } 
        if(tr.triangleattributelist != NULL) 
        {
            free(tr.triangleattributelist);
            tr.triangleattributelist = nullptr; 
        }
        if(tr.trianglearealist != NULL) 
        {
            free(tr.trianglearealist);
            tr.trianglearealist = nullptr; 
        }
        if(tr.neighborlist != NULL) 
        {
            free(tr.neighborlist);
            tr.neighborlist = nullptr; 
        }
        if(tr.segmentlist != NULL) 
        {
            free(tr.segmentlist);
            tr.segmentlist = nullptr; 
        }
        if(tr.segmentmarkerlist != NULL) 
        {
            free(tr.segmentmarkerlist);
            tr.segmentmarkerlist = nullptr; 
        }
        if(tr.holelist != NULL) 
        {
            free(tr.holelist);
            tr.holelist = nullptr; 
        }
        if(tr.regionlist != NULL) 
        {
            free(tr.regionlist);
            tr.regionlist = nullptr; 
        }
        if(tr.edgelist != NULL)    
        {
            free(tr.edgelist);
            tr.edgelist = nullptr; 
        }
        if(tr.edgemarkerlist != NULL) 
        {
            free(tr.edgemarkerlist);
            tr.edgemarkerlist = nullptr; 
        }
        if(tr.normlist != NULL) 
        {
            free(tr.normlist);
            tr.normlist = nullptr; 
        }
    };

    ///@}
    ///@name Life Cycle
    ///@{
    /// Constructor.
    EmbeddedIgaTriangulation();

    /// Destructor.
    virtual ~EmbeddedIgaTriangulation()
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
#endif // KRATOS_EMBEDDED_IGA_TRIANGULATION_H_INCLUDED defined