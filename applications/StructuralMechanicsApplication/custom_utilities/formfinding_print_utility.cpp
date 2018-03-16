//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher 
//
	        

// System includes


// External includes 


// Project includes
#include "includes/checks.h"
#include "formfinding_print_utility.h"
#include "structural_mechanics_application_variables.h"
#include "includes/model_part_io.h"


namespace Kratos
{

	FormfindingPrintUtility::FormfindingPrintUtility(const ModelPart& rModelPart, const Parameters rParameter): 
		mModelPart(rModelPart)
		{}
    void FormfindingPrintUtility::PrintModelPart(){
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
		for( auto& prop: mModelPart.rProperties()){
			prop.Data().Clear();
		}

    	ModelPartIO model_part_io("formfinding_out", IO::WRITE);
    	model_part_io.WriteModelPart(mModelPart);

		ModelPartIO model_part_io_prestress("prestress_data", IO::WRITE);
    	model_part_io_prestress.WriteDataBlock(mModelPart.Elements(), "Element");
		
	}
}  // namespace Kratos.


