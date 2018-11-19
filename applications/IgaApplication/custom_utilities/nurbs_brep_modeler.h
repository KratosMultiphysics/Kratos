#if !defined(KRATOS_NURBS_BREP_MODELER_H_INCLUDED )
#define  KRATOS_NURBS_BREP_MODELER_H_INCLUDED

// System includes
#include <vector>

// Project includes
#include "iga_application_variables.h"

#include "brep_json_io.h"
#include "brep_topology/brep_model.h"
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

    private:
        ///@name Member Variables
        ///@{
        ModelPart&                 m_model_part;
        std::vector<BrepModel>     m_brep_model_vector;

        ///@}
        ///@name Private Operations
        ///@{

        ///@}
        ///@name Un accessible methods
        ///@{

        ///@}
    }; // Class NurbsBrepModeler

}  // namespace Kratos.
#endif // KRATOS_NURBS_BREP_MODELER_APPLICATION_H_INCLUDED defined


