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
#include <string>
#include "containers/model.h"

// External includes

// Project includes
#include "combine_solid_shell_model_parts_process.h"
#include "processes/fast_transfer_between_model_parts_process.h"


namespace Kratos
{
	void CombineSolidShellModelPartsProcess::ExecuteInitialize() {

		// 1. Get (trimming) curve from shell model part that lies on the coupling inerface of the shell
		// Should be done before CombineModelParts(), because CombineModelParts() reorders geometry IDs
		ModelPart& IgaModelPart = _Model.GetModelPart("IgaModelPart");
		Geometry<Node>& ShellCurveOnCouplingInterface = IgaModelPart.GetGeometry(6);

		// 2. Get integration points of the solid part model lying on the coupling inerface
		ModelPart& SolidModelPart_NeumannBC = _Model.GetModelPart("NurbsMesh").GetSubModelPart("Neumann_BC");

		// 3. Project each integration point from Solid to Shell midsurface 
		const SizeType num_integration_points = SolidModelPart_NeumannBC.Conditions().size();
		PointerVector<Geometry<Node>> master_quadrature_points_geometry; 
		CoordinatesArrayType integration_points_global_coords_vector(3);
		IntegrationPointsArrayType integration_points_slave_local_coords_vector(num_integration_points);
		CoordinatesArrayType local_slave_coords = ZeroVector(3);
		size_t ii = 0;

			// Check  Read normals Vectors for Solid Triangles on the Interface
			//vector<array_1d<double, 3>> vectorsNomralVectors(SolidModelPart_NeumannBC.NumberOfConditions());
			//ReadNormalVectors("Queso_NormalVectors_Interface.txt", vectorsNomralVectors, SolidModelPart_NeumannBC);

		for (auto it_condition = SolidModelPart_NeumannBC.ConditionsBegin(); it_condition != SolidModelPart_NeumannBC.ConditionsEnd(); ++it_condition) {
			master_quadrature_points_geometry.push_back((it_condition->pGetGeometry()));

			integration_points_global_coords_vector[0] = master_quadrature_points_geometry(ii)->IntegrationPoints()[0].X();
			integration_points_global_coords_vector[1] = master_quadrature_points_geometry(ii)->IntegrationPoints()[0].Y();
			integration_points_global_coords_vector[2] = master_quadrature_points_geometry(ii)->IntegrationPoints()[0].Z();

			ShellCurveOnCouplingInterface.ProjectionPointGlobalToLocalSpace(integration_points_global_coords_vector, local_slave_coords);
			integration_points_slave_local_coords_vector[ii][0] = local_slave_coords[0];
			integration_points_slave_local_coords_vector[ii][1] = local_slave_coords[1];
			integration_points_slave_local_coords_vector[ii][2] = local_slave_coords[2];

			std::vector<array_1d<double, 3>> derivatives(3);
			ShellCurveOnCouplingInterface.GlobalSpaceDerivatives(derivatives, local_slave_coords, 0);
			CoordinatesArrayType& rProjectedPointGlobalCoordinates = derivatives[0];

			ii += 1;
		}

		// 4. Create CreateQuadraturePointGeometries for the shell
		IndexType NumberOfShapeFunctionDerivatives = 3;
		IntegrationInfo rIntegrationInfo = ShellCurveOnCouplingInterface.GetDefaultIntegrationInfo();
		PointerVector<Geometry<Node>> slave_quadrature_points_geometry(num_integration_points);

		ShellCurveOnCouplingInterface.CreateQuadraturePointGeometries(
			slave_quadrature_points_geometry, // GeometriesArrayType
			NumberOfShapeFunctionDerivatives, // IndexType
			integration_points_slave_local_coords_vector,
			rIntegrationInfo); //IntegrationInfo is not used, yet has to be passed to function

		// 5. Create CoupingGeometries
		PointerVector<CouplingGeometry<Node>> Coupled_Quadrature_Point_Geometries(num_integration_points);
		for (SizeType i = 0; i < num_integration_points; ++i) {
			Coupled_Quadrature_Point_Geometries(i) = Kratos::make_shared<CouplingGeometry<Node>>(master_quadrature_points_geometry(i), slave_quadrature_points_geometry(i));
		}

		// 6. Remove conditions which belong to Neumann_BC (NurbsMesh)
		std::vector<size_t> ConditionsIds;
		for (auto& cond_it = _Model.GetModelPart("NurbsMesh").ConditionsBegin(); cond_it != _Model.GetModelPart("NurbsMesh").ConditionsEnd(); ++cond_it) {
			if (_Model.GetModelPart("NurbsMesh").GetSubModelPart("Neumann_BC").HasCondition(cond_it->Id())) {
				ConditionsIds.push_back(cond_it->Id());
			}
		}
		for (auto& cond_it = ConditionsIds.begin(); cond_it != ConditionsIds.end(); ++cond_it) {
			_Model.GetModelPart("NurbsMesh").RemoveConditionFromAllLevels(*cond_it);
		}

		// 7. Combine Model Parts merged the two model parts to one without deleted the original ones
		CombineModelParts();

		// 8. Remove submodel parts which are not longer needed.
		_Model.GetModelPart("CoupledSolidShellModelPart").GetSubModelPart("NurbsMesh").RemoveSubModelPart("Neumann_BC");
		//_Model.GetModelPart("CoupledSolidShellModelPart").GetSubModelPart("IgaModelPart").RemoveSubModelPart("Load_3"); 

		// 9. Assign each coupling geometry to a coupling (penalty/nitsche) condition
		ModelPart& interface_model_part = _Model.GetModelPart("CoupledSolidShellModelPart").GetSubModelPart("CouplingInterface");

		std::string name = _Parameters["coupling_method"].GetString();;
		const Condition& rReferenceCondition = KratosComponents<Condition>::Get(name);

		ModelPart::ConditionsContainerType new_condition_list;
		IndexType condition_id = interface_model_part.GetParentModelPart().NumberOfConditions() + 1; 

		for (SizeType i = 0; i < Coupled_Quadrature_Point_Geometries.size(); ++i) {
			PropertiesPointerType pProperty;
			new_condition_list.push_back(
				rReferenceCondition.Create(condition_id,
					Coupled_Quadrature_Point_Geometries(i), // 
					pProperty));

			for (SizeType ii = 0; i < Coupled_Quadrature_Point_Geometries(i)->size(); ++ii) {
				interface_model_part.AddNode(Coupled_Quadrature_Point_Geometries(i)->pGetPoint(ii));
			}

			condition_id++;
		}

		// 10. Add conditions to interface model part
		interface_model_part.AddConditions(new_condition_list.begin(), new_condition_list.end());

		// 11. Assign the p_properties to the model part's conditions.
		auto& r_conditions_array = interface_model_part.Conditions();
		auto& p_prop = interface_model_part.pGetProperties(3); // TO DO: ID is preset to three (StructuralMaterials.json). Generalize it later
		block_for_each(
			r_conditions_array,
			[&p_prop](Condition& rCondition)
			{ rCondition.SetProperties(p_prop); }
		);

		// 11.5) Assign properties from master and slave geometries to coupling condition (as in nitsche_stabilization_model_part_process.cpp)
		SizeType prop_id = (interface_model_part.ConditionsBegin()->pGetProperties())->Id();
		if (name == "CouplingSolidShellNitscheCondition") {
			Properties::Pointer master_properties = _Model.GetModelPart("CoupledSolidShellModelPart").GetSubModelPart("NurbsMesh").ElementsBegin()->pGetProperties();
			Properties::Pointer slave_properties = _Model.GetModelPart("CoupledSolidShellModelPart").GetSubModelPart("IgaModelPart").ElementsBegin()->pGetProperties();
			interface_model_part.pGetProperties(prop_id)->AddSubProperties(master_properties);
			interface_model_part.pGetProperties(prop_id)->AddSubProperties(slave_properties);
		}

		// 12. Set Flags
		const auto flagX = IgaFlags::FIX_DISPLACEMENT_X;
		const auto flagY = IgaFlags::FIX_DISPLACEMENT_Y;
		const auto flagZ = IgaFlags::FIX_DISPLACEMENT_Z;
		const bool thisbool = true;
		VariableUtils varUtility;
		varUtility.SetFlag(flagX, thisbool, interface_model_part.Conditions());
		varUtility.SetFlag(flagY, thisbool, interface_model_part.Conditions());
		varUtility.SetFlag(flagZ, thisbool, interface_model_part.Conditions());

		// 13. Write points to VTK
		std::string filename = "data/Intergration_Points_on_Solid1_Coupling_Surface.vtk";
		size_t selectpart = 0;
		writeVTKFile(filename, Coupled_Quadrature_Point_Geometries, selectpart);

		filename = "data/Intergration_Points_on_Shell1_Coupling_Surface.vtk";
		selectpart = 1;
		writeVTKFile(filename, Coupled_Quadrature_Point_Geometries, selectpart);

		
	}

	void CombineSolidShellModelPartsProcess::CombineModelParts() {

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
			if (!it->IsIdGeneratedFromString()) {
				it->SetId(CoupleIdCounter + counter);
			}
			counter++;
		}
		CoupleIdCounter = mModelPart.NumberOfGeometries() + 1;
	}

	void CombineSolidShellModelPartsProcess::ReorderNodeIds(size_t& CoupleIdCounter, ModelPart& mModelPart) {

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

	void CombineSolidShellModelPartsProcess::writeVTKFile(const std::string& filename, const PointerVector<CouplingGeometry<Node>> Coupled_Quadrature_Point_Geometries, size_t selectpart) {


		// Open the file in write mode, which creates or overwrites the file
		std::ofstream file(filename, std::ios::trunc);
		if (!file.is_open()) {
			std::cerr << "Error opening file: " << filename << std::endl;
			return;
		}
		// No need to check if the file is open; proceed directly to writing
		auto this_size = Coupled_Quadrature_Point_Geometries.size();
		file << "# vtk DataFile Version 3.0\n";
		file << "Vertices example\n";
		file << "ASCII\n";
		file << "DATASET POLYDATA\n";
		file << "POINTS " << this_size << " double\n";

		for (size_t i = 0; i < Coupled_Quadrature_Point_Geometries.size(); i++) {
			auto& rCoupled_Quadrature_Point_Geometry = *(Coupled_Quadrature_Point_Geometries(i));

			auto& Geometry = *(rCoupled_Quadrature_Point_Geometry.pGetGeometryPart(selectpart));
			auto& integrationpoint = Geometry.IntegrationPoints();
			CoordinatesArrayType XYZ_coos;
			auto point = Geometry.Center();
			file << point.X() << " " << point.Y() << " " << point.Z() << "\n";

			//}
		}

		file.close();
		std::cout << "VTK file written successfully: " << filename << std::endl;
	}

	void CombineSolidShellModelPartsProcess::ReadNormalVectors(const std::string& filename, vector<array_1d<double, 3>>& vectors, ModelPart& SolidModelPart) {

		// Open the text file named 
		std::ifstream inputFile("Queso_NormalVectors_Interface.txt");
		KRATOS_ERROR_IF_NOT(inputFile.is_open()) << ":: ReadNormalVectors :: Error opening the file!" << std::endl;
		std::string line;

		// Count how many nomral vectors are associated with the corresponding integration point 
		size_t NumberOfAssociatedNormalVector = 0;

		while (std::getline(inputFile, line)) {
			// Read each line of Queso_NormalVectors_Interface
			// First 3 numbers are components of normal vector
			//Last 3 numbers are the coordinates of the integration point
			std::stringstream iss(line);
			double num;
			array_1d<double, 3> normal;
			array_1d<double, 3> Coordinates;
			for (int i = 0; i < 6 && iss >> num; ++i) {
				if (i < 3) normal[i] = num;
				if (i >= 3 && i < 6) Coordinates[i - 3] = num;
			}

			// Match normal vector with the correct Integration Point
			double XX, YY, ZZ, diff;
			size_t cond_it = -1;

			for (auto it_condition = SolidModelPart.ConditionsBegin(); it_condition != SolidModelPart.ConditionsEnd(); ++it_condition) {
				cond_it += 1;
				XX = (SolidModelPart.ConditionsBegin() + cond_it)->pGetGeometry()->IntegrationPoints()[0].X();
				YY = (SolidModelPart.ConditionsBegin() + cond_it)->pGetGeometry()->IntegrationPoints()[0].Y();
				ZZ = (SolidModelPart.ConditionsBegin() + cond_it)->pGetGeometry()->IntegrationPoints()[0].Z();

				diff = abs(Coordinates[0] - XX) + abs(Coordinates[1] - YY) + abs(Coordinates[2] - ZZ);
				if (diff < 1E-14) {
					vectors[cond_it] = normal;
					NumberOfAssociatedNormalVector += 1;
					break;
				}
			}
		}
		KRATOS_ERROR_IF(NumberOfAssociatedNormalVector != SolidModelPart.NumberOfConditions()) << "Not all vectors where associated with the corresponding normals";
		inputFile.close();

	}


} // namespace Kratos
