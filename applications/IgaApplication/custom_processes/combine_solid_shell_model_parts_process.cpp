//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

// System includes
#include <iostream>
#include <fstream>
#include "containers/model.h"
// External includes

// Project includes
#include "combine_solid_shell_model_parts_process.h"
#include "utilities/parallel_utilities.h"
#include "processes/fast_transfer_between_model_parts_process.h"

namespace Kratos
{
	void CombineSolidShellModelPartsProcess::ExecuteInitialize(){
		
		ModelPart& SolidModelPart = _Model.GetModelPart("NurbsMesh");
		ModelPart& ShellModelPart = _Model.GetModelPart("IgaModelPart");
		ModelPart& CoupledModelPart = _Model.GetModelPart("CoupledSolidShellModelPart");

		size_t  CoupleIdCounter = 1;
		ReorderGeometryIds(CoupleIdCounter, SolidModelPart);
		ReorderGeometryIds(CoupleIdCounter, ShellModelPart);

		CoupleIdCounter = 1;
		ReorderNodeIds(CoupleIdCounter, SolidModelPart);
		ReorderNodeIds(CoupleIdCounter, ShellModelPart);

		CoupleIdCounter = 1;
		ReorderElementIds(CoupleIdCounter, SolidModelPart);
		ReorderElementIds(CoupleIdCounter, ShellModelPart);

		CoupleIdCounter = 1;
		ReorderConditionIds(CoupleIdCounter, SolidModelPart);
		ReorderConditionIds(CoupleIdCounter, ShellModelPart);

		CoupleIdCounter = 1;
		ReorderMasterSlaveConstraintIds(CoupleIdCounter, SolidModelPart);
		ReorderMasterSlaveConstraintIds(CoupleIdCounter, ShellModelPart);

		CoupleIdCounter = 1;
		ReorderPropertyIds(CoupleIdCounter, SolidModelPart);
		ReorderPropertyIds(CoupleIdCounter, ShellModelPart);

		RecursiveAddEntities(CoupledModelPart.CreateSubModelPart("NurbsMesh"), SolidModelPart);
		auto& NodalSolutionStepVariablesList = SolidModelPart.GetNodalSolutionStepVariablesList();
		for (auto& r_var : NodalSolutionStepVariablesList) {
			CoupledModelPart.GetSubModelPart("NurbsMesh").AddNodalSolutionStepVariable(r_var);
		}

		RecursiveAddEntities(CoupledModelPart.CreateSubModelPart("IgaModelPart"), ShellModelPart);
		NodalSolutionStepVariablesList = ShellModelPart.GetNodalSolutionStepVariablesList();
		for (auto& r_var : NodalSolutionStepVariablesList) {
			CoupledModelPart.GetSubModelPart("IgaModelPart").AddNodalSolutionStepVariable(r_var);
		}

		CoupledModelPart.SetNodalSolutionStepVariablesList();

		auto pProperties = *(CoupledModelPart.PropertiesBegin().base());

		// Iterate over conditions
		block_for_each(CoupledModelPart.Conditions(), [&pProperties](Condition& rCondition) {
			rCondition.SetProperties(pProperties);
			});

		//// Iterate over elements
		//&pProperties;
		//block_for_each(CoupledModelPart.Elements(), [&pProperties](Element& rElement) {
		//	rElement.SetProperties(pProperties);
		//	});
	}

	void CombineSolidShellModelPartsProcess::RecursiveAddEntities(ModelPart& rTwinModelPart, ModelPart& rOriginModelPart)
	{
		// Lambda to transfer the model parts
		auto transfer_lambda = [](ModelPart& rTwinModelPart, ModelPart& rOriginModelPart) {
			FastTransferBetweenModelPartsProcess(rTwinModelPart, rOriginModelPart).Execute();

			// Copy properties
			for (auto it_prop = rOriginModelPart.PropertiesBegin(); it_prop < rOriginModelPart.PropertiesEnd(); it_prop++) {
				if (!rTwinModelPart.HasProperties(it_prop->Id())) {
					rTwinModelPart.AddProperties(*(it_prop.base()));
					//std::cout << "Minas: Adding properties for " << rTwinModelPart.GetRootModelPart().Name() << "." << rTwinModelPart.Name() << std::endl;
				}
			}
			};

		// Recursively add of ModelParts to the list
		if (rOriginModelPart.NumberOfSubModelParts() > 0) {
			for (auto& r_sub_model_part : rOriginModelPart.SubModelParts()) {
				auto& r_sub_destination_model_part = rTwinModelPart.CreateSubModelPart(r_sub_model_part.Name());
				RecursiveAddEntities(r_sub_destination_model_part, r_sub_model_part);
			}
		}
		transfer_lambda(rTwinModelPart, rOriginModelPart);
	}

	void CombineSolidShellModelPartsProcess::ReorderGeometryIds(size_t& CoupleIdCounter, ModelPart& mModelPart) {

		IndexType counter = 0;
		for (auto it = mModelPart.Geometries().begin(); it != mModelPart.Geometries().end(); it++) {
			it->SetId(CoupleIdCounter + counter);
			counter++;
		}
		CoupleIdCounter = mModelPart.NumberOfGeometries() + 1;
	}

	void CombineSolidShellModelPartsProcess::ReorderNodeIds(size_t& CoupleIdCounter, ModelPart& mModelPart){
		
		IndexType counter = 0;
		for (auto it = mModelPart.Nodes().begin(); it != mModelPart.Nodes().end(); it++) {
			it->SetId(CoupleIdCounter + counter);
			counter++;
		}
		CoupleIdCounter = mModelPart.NumberOfNodes() + 1;
	}

	void CombineSolidShellModelPartsProcess::ReorderElementIds(size_t& CoupleIdCounter, ModelPart& mModelPart) {

		IndexType counter = 0;
		for (auto it = mModelPart.Elements().begin(); it != mModelPart.Elements().end(); it++) {
			it->SetId(CoupleIdCounter + counter);
			counter++;
		}
		CoupleIdCounter = mModelPart.NumberOfElements() + 1;
	}

	void CombineSolidShellModelPartsProcess::ReorderConditionIds(size_t& CoupleIdCounter, ModelPart& mModelPart) {

		IndexType counter = 0;
		for (auto it = mModelPart.Conditions().begin(); it != mModelPart.Conditions().end(); it++) {
			it->SetId(CoupleIdCounter + counter);
			counter++;
		}
		CoupleIdCounter = mModelPart.NumberOfConditions() + 1;
	}

	void CombineSolidShellModelPartsProcess::ReorderMasterSlaveConstraintIds(size_t& CoupleIdCounter, ModelPart& mModelPart) {

		IndexType counter = 0;
		for (auto it = mModelPart.MasterSlaveConstraints().begin(); it != mModelPart.MasterSlaveConstraints().end(); it++) {
			it->SetId(CoupleIdCounter + counter);
			counter++;
		}
		CoupleIdCounter = mModelPart.NumberOfMasterSlaveConstraints() + 1;
	}

	void CombineSolidShellModelPartsProcess::ReorderPropertyIds(size_t& CoupleIdCounter, ModelPart& mModelPart) {

		IndexType counter = 0;
		for (auto it = mModelPart.PropertiesBegin(); it != mModelPart.PropertiesEnd(); it++) {
			//std::cout << "Property of " << mModelPart.Name() << " with Id: " << it->Id();
			it->SetId(CoupleIdCounter + counter);
			//std::cout << " will have Id: " << it->Id() << std::endl;
			counter++;
		}
		CoupleIdCounter = mModelPart.NumberOfProperties() + 1;
	}
	} // namespace Kratos
