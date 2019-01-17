#if !defined(KRATOS_EMBEDDED_IGA_MODELER_H_INCLUDED )
#define  KRATOS_EMBEDDED_IGA_MODELER_H_INCLUDED

// System includes

// External includes
#include "anurbs.h"
#include "containers/model.h"
// #include "triangle.h"

// Project includes
#include "iga_application_variables.h"
#include "nurbs_brep_modeler.h"
#include "includes/model_part.h"

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
        
        void CreateTessellationCurve(ANurbs::Pointer<ANurbs::CurveTessellation3D>& rTessellation); 
        void CreateTessellationParameterCurve(); 
        void CreateElements2D(ModelPart& rSkinModelPart); 
        std::vector<double> PrintNodesX();
        std::vector<double> PrintNodesY();

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
#endif // KRATOS_EMBEDDED_IGA_H_INCLUDED defined


