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
#include "formfinding_IO_utility.h"
#include "structural_mechanics_application_variables.h"
#include "includes/model_part_io.h"
#include "includes/model_part.h"


namespace Kratos
{

	FormfindingIOUtility::FormfindingIOUtility(ModelPart& rModelPart, const Parameters rParameter): 
		mModelPart(rModelPart)
		{}
    void FormfindingIOUtility::PrintModelPart(){
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

		// Write ModelPart
    	ModelPartIO model_part_io("formfinding_out", IO::WRITE);
    	model_part_io.WriteModelPart(mModelPart);
		
	}

	void FormfindingIOUtility::PrintPrestressData(){

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

		// write prestress data
		ModelPartIO model_part_io_prestress("prestress_data", IO::WRITE);
    	model_part_io_prestress.WriteDataBlock(mModelPart.Elements(), "Element");

	}

	void FormfindingIOUtility::ReadPrestressData(){
		ModelPartIO model_part_io("prestress_data", IO::READ);
		model_part_io.ReadInitialValues(mModelPart);
		std::cout<<"Prestress read"<<std::endl;
	}
}  // namespace Kratos.


