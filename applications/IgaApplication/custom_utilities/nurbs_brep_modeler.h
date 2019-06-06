#if !defined(KRATOS_NURBS_BREP_MODELER_H_INCLUDED )
#define  KRATOS_NURBS_BREP_MODELER_H_INCLUDED

// System includes
#include <vector>
#include <fstream>

// Project includes
#include "iga_application_variables.h"

#include "brep_json_io.h"
#include "brep_topology/brep_model.h"

#include "integration_utilities/iga_integration_utilities.h"

#include "spatial_containers/bins_dynamic_objects.h"
#include "custom_utilities/search_utilities/bins_iga_configure.h"
#include "custom_utilities/search_utilities/bins_iga_object.h"

#include "includes/model_part.h"

namespace Kratos
{
    struct ElementConditionParameter {
        std::string type;
        std::string name;
        int properties_id;
        int shape_function_derivatives_order;
        std::vector<std::string> variables;

        ElementConditionParameter(
            const std::string& rType,
            const std::string& rName,
            const int& rPropertiesId,
            const int& rShapeFunctionDerivativesOrder,
            std::vector<std::string> rVariables)
        {
            type = rType;
            name = rName;
            properties_id = rPropertiesId;
            shape_function_derivatives_order = rShapeFunctionDerivativesOrder;
            variables = rVariables;
        }
    };

    class NurbsBrepModeler
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of KratosNurbsTestcaseApplication
        KRATOS_CLASS_POINTER_DEFINITION(NurbsBrepModeler);

        ///@}
        ///@name functions
        ///@{

        /**
        * Imports and adds a brep model with the use of the BrepJSON_IO
        * @param rBrepJSON_IO the IO reader
        */
        void ImportGeometry(
            BrepJsonIO& rBrepJsonIO,
            Parameters& rNurbsBrepGeometryJson);

        // here shall be added the functionality for more import option
        // especially directly creating the model in python/ Rhino/ GiD...
        //void ImportGeometry(std::vector<BrepModel>& rBrepModel);

        void ImportModelPart(ModelPart& model_part, Parameters& rModelPartParameters);

        void GetInterfaceConditions(
            ModelPart& rModelPart,
            ModelPart& rIgaModelPart,
            ModelPart& rInterfaceConditionsModelPart,
            const std::string& rConditionName);

        void GetInterfaceConditionsDEM(
            ModelPart& rExternalModelPart,
            ModelPart& rIgaModelPart,
            ModelPart& rInterfaceConditionsModelPart,
            const std::string& rConditionName,
            const double ShapeFunctionDerivativesOrder,
            const double SearchRadius,
            const double Accuracy,
            const double Tolerance,
            const double NumberOfIterations);

        void ExportGeometry();

        void getTolerance();

        ///@}
        ///@name Life Cycle
        ///@{
        /// Constructor.
        NurbsBrepModeler(ModelPart& rModelPart);

        /// Destructor.
        virtual ~NurbsBrepModeler()
        { };

        ///@}
    protected:
        ModelPart&                 m_model_part;
        std::vector<BrepModel>     m_brep_model_vector;

    private:
        ///@name Member Variables
        ///@{

        ///@}
        ///@name Private Operations
        ///@{

        const BrepVertex& GetBrepVertex(int& rBrepId) const;
        const BrepEdge& GetBrepEdge(int& rBrepId) const;
        const BrepFace& GetBrepFace(int& rBrepId) const;
        ///@}
        ///@name Un accessible methods
        ///@{

        ///@}
    }; // Class NurbsBrepModeler

}  // namespace Kratos.
#endif // KRATOS_NURBS_BREP_MODELER_APPLICATION_H_INCLUDED defined


