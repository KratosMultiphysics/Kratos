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
#include "processes/fast_transfer_between_model_parts_process.h"


namespace Kratos
{
	void CombineSolidShellModelPartsProcess::ExecuteInitialize() {

		// Get (trimming) curve from shell model part that lies on the coupling inerface 
		// Should be done before CombineModelParts(), because CombineModelParts() 
		// reorders geometry IDs
		ModelPart& IgaModelPart = _Model.GetModelPart("IgaModelPart");
		Geometry<Node>& ShellCurveOnCouplingInterface = IgaModelPart.GetGeometry(6);

		// Get integration points of the solid part model lying on the coupling inerface 
		ModelPart& SolidModelPart_NeumannBC = _Model.GetModelPart("NurbsMesh").GetSubModelPart("Neumann_BC");
				
		// For the solid part, there is 1 quadrature point per condition 
		const SizeType num_integration_points = SolidModelPart_NeumannBC.Conditions().size();
		PointerVector<Geometry<Node>> master_quadrature_points_geometry; //TO DO: geometry type may not be correct

		PointerVector<Geometry<Node>> slave_quadrature_points_geometry(num_integration_points); //TO DO: idem
		CoordinatesArrayType integration_points_global_coords_vector(3);
		IntegrationPointsArrayType integration_points_slave_local_coords_vector(num_integration_points);
		CoordinatesArrayType local_slave_coords = ZeroVector(3); 
		size_t i = 0;
			
		for (auto it_condition = SolidModelPart_NeumannBC.ConditionsBegin(); it_condition != SolidModelPart_NeumannBC.ConditionsEnd(); ++it_condition) {
			master_quadrature_points_geometry.push_back( (it_condition->pGetGeometry() )) ;

			integration_points_global_coords_vector[0] = master_quadrature_points_geometry(i)->IntegrationPoints()[0].X();
			integration_points_global_coords_vector[1] = master_quadrature_points_geometry(i)->IntegrationPoints()[0].Y();
			integration_points_global_coords_vector[2] = master_quadrature_points_geometry(i)->IntegrationPoints()[0].Z();

			ShellCurveOnCouplingInterface.ProjectionPointGlobalToLocalSpace(integration_points_global_coords_vector, local_slave_coords);
			integration_points_slave_local_coords_vector[i][0] = local_slave_coords[0];
			integration_points_slave_local_coords_vector[i][1] = local_slave_coords[1];
			integration_points_slave_local_coords_vector[i][2] = local_slave_coords[2];

			std::vector<array_1d<double, 3>> derivatives(3);
			ShellCurveOnCouplingInterface.GlobalSpaceDerivatives(derivatives, local_slave_coords, 0);
			CoordinatesArrayType& rProjectedPointGlobalCoordinates = derivatives[0];
			
			i += 1;
		}
		IndexType NumberOfShapeFunctionDerivatives = 2;

		IntegrationInfo rIntegrationInfo = ShellCurveOnCouplingInterface.GetDefaultIntegrationInfo();
		
		ShellCurveOnCouplingInterface.CreateQuadraturePointGeometries(
				slave_quadrature_points_geometry, // GeometriesArrayType
				NumberOfShapeFunctionDerivatives, // IndexType
				integration_points_slave_local_coords_vector,
				rIntegrationInfo); //IntegrationInfo is not used, yet has to be passed to function
		
		PointerVector<CouplingGeometry<Node>> Coupled_Quadrature_Point_Geometries(num_integration_points);
		
		for (SizeType i = 0; i < num_integration_points; ++i) { 
			Coupled_Quadrature_Point_Geometries(i) = Kratos::make_shared<CouplingGeometry<Node>>(master_quadrature_points_geometry(i), slave_quadrature_points_geometry(i)); 
			if (i < 1) {
				std::cout << " Print the geometry Info for master of coupled geometry" << master_quadrature_points_geometry(i)->Info() << std::endl;
				std::cout << " Print the geometry Info for slave of coupled geometry" << slave_quadrature_points_geometry(i)->Info() << std::endl;

			}
		}		

		// Remove conditions which belong to Neumann_BC (NurbsMesh)
		std::vector<size_t> ConditionsIds;
		for (auto& cond_it = _Model.GetModelPart("NurbsMesh").ConditionsBegin(); cond_it != _Model.GetModelPart("NurbsMesh").ConditionsEnd(); ++cond_it) {
			if (_Model.GetModelPart("NurbsMesh").GetSubModelPart("Neumann_BC").HasCondition(cond_it->Id())) {
				ConditionsIds.push_back(cond_it->Id());
			}
		}
		for (auto& cond_it = ConditionsIds.begin(); cond_it != ConditionsIds.end(); ++cond_it) {
			_Model.GetModelPart("NurbsMesh").RemoveConditionFromAllLevels(*cond_it);
		}

		//Combine Model Parts merged the two model parts to one without deleted the original ones
		CombineModelParts();

		_Model.GetModelPart("CoupledSolidShellModelPart").GetSubModelPart("NurbsMesh").RemoveSubModelPart("Neumann_BC");
		//_Model.GetModelPart("CoupledSolidShellModelPart").GetSubModelPart("IgaModelPart").RemoveSubModelPart("Load_3");
		
		// assign each coupling geometry to the coupling condition
		ModelPart& interface_model_part = _Model.GetModelPart("CoupledSolidShellModelPart").GetSubModelPart("CouplingInterface");
		std::string name = "CouplingSolidShellPenaltyCondition";
		//std::string name = "CouplingPenaltyCondition";
		const Condition& rReferenceCondition = KratosComponents<Condition>::Get(name);

		ModelPart::ConditionsContainerType new_condition_list;
		IndexType condition_id = interface_model_part.GetParentModelPart().NumberOfConditions() + 1; // To be deleted

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

		interface_model_part.AddConditions(new_condition_list.begin(), new_condition_list.end());

		// Assign the p_properties to the model part's elements and conditions.
		auto& r_conditions_array = interface_model_part.Conditions();

		auto& p_prop = interface_model_part.pGetProperties(3); // TO DO: ID is preset to three (StructuralMaterials.json). Generalize it later

		block_for_each(
			r_conditions_array,
			[&p_prop](Condition& rCondition)
			{ rCondition.SetProperties(p_prop); }
		);

		// Set Flags
		const auto flagX = IgaFlags::FIX_DISPLACEMENT_X;
		const auto flagY = IgaFlags::FIX_DISPLACEMENT_Y;
		const auto flagZ = IgaFlags::FIX_DISPLACEMENT_Z;
		const bool thisbool = true;
		VariableUtils varUtility;
		varUtility.SetFlag(flagX, thisbool, interface_model_part.Conditions());
		varUtility.SetFlag(flagY, thisbool, interface_model_part.Conditions());
		varUtility.SetFlag(flagZ, thisbool, interface_model_part.Conditions());
				
		// Write points to VTK
		std::string filename = "data/Intergration_Points_on_Solid1_Coupling_Surface.vtk";
		size_t selectpart = 0;
		writeVTKFile(filename, Coupled_Quadrature_Point_Geometries, selectpart);

		filename = "data/Intergration_Points_on_Shell1_Coupling_Surface.vtk";
		selectpart = 1;
		writeVTKFile(filename, Coupled_Quadrature_Point_Geometries, selectpart);

	}

	void CombineSolidShellModelPartsProcess::CombineModelParts(){
		
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

	void CombineSolidShellModelPartsProcess::writeVTKFile(const std::string& filename,const PointerVector<CouplingGeometry<Node>> Coupled_Quadrature_Point_Geometries, size_t selectpart) {
		
		
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
			//Geometry.GlobalCoordinates(XYZ_coos,integrationpoint[0]);
			auto point = Geometry.Center();

			std::vector<array_1d<double, 3>> derivatives(3);
			//auto ParentGeometry = Geometry.GetGeometryParent(0);
			//ParentGeometry.GlobalSpaceDerivatives(derivatives, integrationpoint[0], 0);
			CoordinatesArrayType& rProjectedPointGlobalCoordinates = derivatives[0];
			//for (auto& it_intpoints = integrationpoint.begin(); it_intpoints != integrationpoint.end(); it_intpoints++) {
				//file << XYZ_coos[0] << " " << XYZ_coos[1] << " " << XYZ_coos[2] << "\n";
			file << point.X() << " " << point.Y() << " " << point.Z() << "\n";

			//}
		}

		file.close();
		std::cout << "VTK file written successfully: " << filename << std::endl;
	}

	//void CombineSolidShellModelPartsProcess::projectIntergrationPointsToMidsurface(GeometryType& ShellCurveOnCouplingInterface,
	//	GeometryType::IntegrationPointsArrayType& IntegrationPointsArray, GeometryType::IntegrationPointsArrayType& GlobalCoordinates_IntegrationPoints_Shell)
	//{
	//	double X, Y, Z;
	//	for (auto& it_IntegrationPoints = IntegrationPointsArray.begin(); it_IntegrationPoints != IntegrationPointsArray.end(); it_IntegrationPoints++) {
	//		// The coordinates from IntegrationPointType must be converted GeometryType::CoordinatesArrayType to be passed to ProjectedPointGlobalToLocalSpace
	//		X = it_IntegrationPoints->X();
	//		Y = it_IntegrationPoints->Y();
	//		Z = it_IntegrationPoints->Z();
	//		std::array<double, 3> array1 = { X,Y,Z };
	//		CoordinatesArrayType IntegrationPointGlobalCoos(3, array1);
	//		array_1d<double, 3> ProjectedPointLocalCoordinates(1, { 0 });

	//		ShellCurveOnCouplingInterface.ProjectionPointGlobalToLocalSpace(IntegrationPointGlobalCoos, ProjectedPointLocalCoordinates);
	//		std::vector<array_1d<double, 3>> derivatives(3);
	//		ShellCurveOnCouplingInterface.GlobalSpaceDerivatives(derivatives, ProjectedPointLocalCoordinates, 0);
	//		CoordinatesArrayType& rProjectedPointGlobalCoordinates = derivatives[0];

	//		GeometryType::IntegrationPointType GlobalCoordinate_IntegrationPoint(rProjectedPointGlobalCoordinates[0], rProjectedPointGlobalCoordinates[1],
	//			rProjectedPointGlobalCoordinates[2], 0);
	//		GlobalCoordinates_IntegrationPoints_Shell.push_back(GlobalCoordinate_IntegrationPoint);

	//	}
	//}

	//void CombineSolidShellModelPartsProcess::GetIntergrationPointsFromInterfacePart
	//(ModelPart& mModelPart, GeometryType::IntegrationPointsArrayType& IntegrationPointsArray) {
	//	for (auto& it_condition = mModelPart.ConditionsBegin(); it_condition != mModelPart.ConditionsEnd(); it_condition++) {
	//		auto cGeometry = it_condition->GetGeometry();
	//		const GeometryType::IntegrationPointsArrayType& IntegrationPointsArrayPerCondition = it_condition->GetGeometry().IntegrationPoints();
	//		for (auto& it_intpoints = IntegrationPointsArrayPerCondition.begin(); it_intpoints != IntegrationPointsArrayPerCondition.end(); it_intpoints++) {
	//			IntegrationPointsArray.push_back(*it_intpoints);
	//		}
	//		for (auto& it_geometry = cGeometry.begin(); it_geometry != cGeometry.end(); it_geometry++) {
	//			auto cGeometry = it_condition->GetGeometry();
	//		}
	//	}

	//}

	} // namespace Kratos
