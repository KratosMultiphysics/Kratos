

#if !defined(KRATOS_ADD_VISUALIZATION_UTILITIES_INCLUDED )
#define KRATOS_ADD_VISUALIZATION_UTILITIES_INCLUDED

// System includes
#include <string>
#include <iostream> 
#include <algorithm>

// Project includes 
#include "includes/define.h"
#include "pfem_2_application_variables.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h" 
#include "includes/ublas_interface.h"
#include "includes/variables.h" 
#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "geometries/point_2d.h"


namespace Kratos
{
	//this class is to be modified by the user to customize the interpolation process
	//template< unsigned int TDim>
	class VisualizationUtilities
	{
	public:
	
		KRATOS_CLASS_POINTER_DEFINITION(VisualizationUtilities);

		VisualizationUtilities()
		{
			KRATOS_TRY	
			std::cout << "Hello, I am the constructor of the Visualization 2d Utility" << std::endl;
			KRATOS_CATCH("")	
		}
		

		~VisualizationUtilities()
		{}

		
    //**********************************************************************************************
    //**********************************************************************************************


    void VisualizationModelPart(ModelPart& rCompleteModelPart, ModelPart& rEulerianModelPart, ModelPart & rLagrangianModelPart)
    {
        KRATOS_TRY;

        rCompleteModelPart.Elements() = rEulerianModelPart.Elements();
        rCompleteModelPart.Nodes() = rEulerianModelPart.Nodes();

		KRATOS_WATCH("hola11")

        unsigned int id;
        if(rEulerianModelPart.Nodes().size()!= 0)
            id = (rEulerianModelPart.Nodes().end() - 1)->Id() + 1;
        else
            id = 1;

        //preallocate the memory needed
        int tot_nodes = rEulerianModelPart.Nodes().size() + rLagrangianModelPart.Nodes().size();
        rCompleteModelPart.Nodes().reserve( tot_nodes );

		KRATOS_WATCH("hola12")

        //note that here we renumber the nodes
        KRATOS_WATCH(rLagrangianModelPart.Nodes().size())
        for (ModelPart::NodesContainerType::iterator node_it = rLagrangianModelPart.NodesBegin();
                node_it != rLagrangianModelPart.NodesEnd(); node_it++)
        {
            node_it->SetId(id++);
            rCompleteModelPart.AddNode(*(node_it.base()));
        }
        
        KRATOS_WATCH("hola13")

        KRATOS_CATCH("");
    }

    //**********************************************************************************************
    //**********************************************************************************************
		
		
	protected:


	private:


	};

}  // namespace Kratos.

#endif // KRATOS_ADD_VISUALIZATION_UTILITIES_INCLUDED  defined 


