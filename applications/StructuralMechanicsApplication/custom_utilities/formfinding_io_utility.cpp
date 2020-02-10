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


// System includes


// External includes


// Project includes
#include "formfinding_io_utility.h"
#include "structural_mechanics_application_variables.h"
#include "includes/model_part_io.h"
#include "includes/model_part.h"


namespace Kratos
{

	FormfindingIOUtility::FormfindingIOUtility(ModelPart& rModelPart, const Parameters rParameter):
		mModelPart(rModelPart)
		{}
    void FormfindingIOUtility::PrintModelPart()
    {
		KRATOS_INFO("FormfindingIOUtility") << "Attention: Removing internal modelpart data. "
			<< "The modelpart will not work in the same way as before." << std::endl;
		// erase nodal data
    	mModelPart.GetNodalSolutionStepVariablesList().clear();

    	// erase elemental data
    	for( auto& ele: mModelPart.Elements()){
    	    if(ele.Has(MEMBRANE_PRESTRESS)){
    	        const Matrix membrane_prestress(ele.GetValue(MEMBRANE_PRESTRESS));
    	        ele.Data().Clear();
    	        ele.SetValue(MEMBRANE_PRESTRESS, membrane_prestress);
    	    }
    	    else
    	        ele.Data().Clear();
    	}

    	// erase conditional data
    	for (auto& cond: mModelPart.Conditions())
			cond.Data().Clear();

		// erase properties
		for (auto& prop: mModelPart.rProperties())
			prop.Data().Clear();

		// Write ModelPart
    	ModelPartIO model_part_io("formfinding_out", IO::WRITE);
    	model_part_io.WriteModelPart(mModelPart);
	}

	void FormfindingIOUtility::PrintPrestressData()
    {
        KRATOS_ERROR << "This function is currently not working, use \"PrintModelPart\" instead"  << std::endl;
		// erase elemental data
    	for( auto& ele: mModelPart.Elements()){
    	    if(ele.Has(MEMBRANE_PRESTRESS)){
    	        const Matrix membrane_prestress(ele.GetValue(MEMBRANE_PRESTRESS));
    	        ele.Data().Clear();
    	        ele.SetValue(MEMBRANE_PRESTRESS, membrane_prestress);
    	    }
    	    else
    	        ele.Data().Clear();
    	}

		// write prestress data
		ModelPartIO model_part_io_prestress("prestress_data", IO::WRITE);
    	// model_part_io_prestress.WriteDataBlock(mModelPart.Elements(), "Element");

	}

	void FormfindingIOUtility::ReadPrestressData()
    {
		ModelPartIO model_part_io("prestress_data", IO::READ);
		model_part_io.ReadInitialValues(mModelPart);
	}
}  // namespace Kratos.


