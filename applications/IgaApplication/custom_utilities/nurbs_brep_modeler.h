#if !defined(KRATOS_NURBS_BREP_MODELER_H_INCLUDED )
#define  KRATOS_NURBS_BREP_MODELER_H_INCLUDED

// System includes

// Project includes
#include "brep_json_io.h"
#include "brep_topology/brep_model.h"

namespace Kratos
{
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

        ///@} 
        ///@name Life Cycle 
        ///@{ 
        /// Constructor.
        NurbsBrepModeler::NurbsBrepModeler(ModelPart& rModelPart)
            : m_model_part(rModelPart)
        { };

        /// Destructor.
        NurbsBrepModeler::~NurbsBrepModeler()
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
#endif // KRATOS_NURBS_BREP_MODELER_APPLICATION_H_INCLUDED  defined 


