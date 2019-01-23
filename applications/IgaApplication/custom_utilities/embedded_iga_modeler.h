#if !defined(KRATOS_EMBEDDED_IGA_MODELER_H_INCLUDED )
#define  KRATOS_EMBEDDED_IGA_MODELER_H_INCLUDED

// System includes

// External includes
#include "anurbs.h"

// Project includes
#include "containers/model.h"
#include "includes/model_part.h"
#include "iga_application_variables.h"
#include "nurbs_brep_modeler.h"
#include "embedded_iga_triangulation.h"
#include "embedded_iga_error_estimation.h"
#include "meshing_application.h"

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
        
        void CreateTessellationCurve(
            ANurbs::Pointer<ANurbs::CurveTessellation3D>& rTessellation); 
        
        void CreateTessellationParameterCurve(
            std::vector<array_1d<double, 3> >& rPolygon);

        void CreateElements2D(
            ModelPart& rSkinModelPart);
        
        std::vector<Matrix> Triangulate(); 
        
        
        std::vector<std::vector<double>> PrintCurveTessellationPoints(); 
        std::vector<std::vector<double>> PrintTriangulationPoints(); 
        std::vector<std::vector<double>> PrintParameterCurveTessellationPoints(); 
        std::vector<std::vector<double>> PrintGaussPoints();
        std::vector<std::vector<double>> PrintMappedGaussPoints(); 

        void TestTriangle(); 
        

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
#endif // KRATOS_EMBEDDED_IGA_TRIANGULATION_H_INCLUDED defined


