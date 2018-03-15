//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Anna Rehr
//

#if !defined(KRATOS_FORMFINDING_UTILIY_H_INCLUDED )
#define  KRATOS_FORMFINDING_UTILIY_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/properties.h"


namespace Kratos
{
	namespace FormfindingUtilities
	{
		void PrintModelPart(const ModelPart& rModelPart){
			// erase nodal data
        	mModelPart.GetNodalSolutionStepVariablesList().clear();

        	// erase elemental data
        	for( auto& ele: mModelPart.Elements()){
        	    const Variable<Matrix> variable = KratosComponents<Variable<Matrix>>::Get("MEMBRANE_PRESTRESS");
        	    if(ele.Has(variable)){
        	        const Matrix membrane_prestress(ele.GetValue(MEMBRANE_PRESTRESS));
        	        ele.Data().Clear();
        	        ele.SetValue(MEMBRANE_PRESTRESS, membrane_prestress);
        	    }
        	    else
        	        ele.Data().Clear();

        	}

        	// erase conditional data
        	for( auto& cond: mModelPart.Conditions()){
        	        cond.Data().Clear();
        	}

			// erase properties

        	ModelPartIO model_part_io("formfinding_out", IO::WRITE);
        	model_part_io.WriteModelPart(mModelPart);
			
		}

		void PrintElementalData(ModelPart& rModelPart){
			const TVariableType variable = KratosComponents<Matrix>::Get(MEMBRANE_PRESTRESS);
        	(*mpStream) << "Begin "<<rObjectName<<"alData "<<variable.Name()<<std::endl;
        	for(auto& object : rThisObjectContainer){
        	    if(object.Has(variable)){
        	        (*mpStream)<<object.Id()<<"\t"<<object.GetValue(variable)<<std::endl;
        	    }
        	}
        	(*mpStream)<<"End "<<rObjectName<<"alData\n"<<std::endl;
		}


	}  // namespace Formfinding Utilities.
  
}  // namespace Kratos.

#endif // KRATOS_FORMFINDING_UTILIY_H_INCLUDED  defined 


