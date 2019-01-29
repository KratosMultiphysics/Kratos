#if !defined(KRATOS_EMBEDDED_IGA_MODELER_H_INCLUDED )
#define  KRATOS_EMBEDDED_IGA_MODELER_H_INCLUDED

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
#include "anurbs.h"
#include "triangle.h"   

// Project includes
#include "containers/model.h"
#include "includes/model_part.h"
#include "iga_application_variables.h"
#include "nurbs_brep_modeler.h"
#include "embedded_iga_tessellation.h"
#include "embedded_iga_triangulation.h"
#include "embedded_iga_error_estimation.h"


namespace Kratos
{
    class EmbeddedIgaModeler : public NurbsBrepModeler
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of KratosNurbsTestcaseApplication
        KRATOS_CLASS_POINTER_DEFINITION(EmbeddedIgaModeler);

        ///@}
        ///@name functions
        ///@{
        
        void CreateElements2D(
            ModelPart& rSkinModelPart);
        
        std::vector<Matrix> TriangulateEmpire(); 
        
        
        std::vector<std::vector<double>> PrintCurveTessellationPoints(); 
        std::vector<std::vector<double>> PrintTriangulationPoints(); 
        std::vector<std::vector<double>> PrintParameterCurveTessellationPoints(); 
        std::vector<std::vector<double>> PrintGaussPoints();
        std::vector<std::vector<double>> PrintMappedGaussPoints(); 

        std::vector<std::vector<double>> Triangulate();  
        std::vector<std::vector<double>> MapTriangulationVertices(); 



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
        EmbeddedIgaModeler(ModelPart& rModelPart);

        /// Destructor.
        virtual ~EmbeddedIgaModeler() override
        {};

        ///@}
    protected:

    private:
        ///@name Member Variables
        ///@{
        ModelPart&                 m_model_part;

        ///@}
        ///@name Private Operations
        ///@{

        ///@}
        ///@name Un accessible methods
        ///@{

        ///@}

    };

}  // namespace Kratos.
#endif // KRATOS_EMBEDDED_IGA_MODELER_H_INCLUDED defined


