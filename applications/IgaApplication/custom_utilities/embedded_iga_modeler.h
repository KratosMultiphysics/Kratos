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
#include "embedded_iga/embedded_iga_tessellation.h"
#include "embedded_iga/embedded_iga_triangulation.h"
#include "embedded_iga/embedded_iga_error_estimation.h"
#include "embedded_iga/embedded_iga_mapper.h"


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

        void CreateElements3D(ModelPart& rSkinModelPart);

        std::vector<std::vector<double>> PrintMappedPoints();
        std::vector<std::vector<double>> PrintParametricTessellation(); 
        std::vector<std::vector<double>> PrintParametricTriangulation();
        std::vector<std::vector<double>> TestCreateElements3D();




            
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


