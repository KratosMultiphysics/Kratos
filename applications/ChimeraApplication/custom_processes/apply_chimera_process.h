// ==============================================================================
//  ChimeraApplication
//
//  License:         BSD License
//                   license: ChimeraApplication/license.txt
//
//  Main authors:    Aditya Ghantasala, https://github.com/adityaghantasala
//
// ==============================================================================
//

#if !defined(KRATOS_CUSTOM_APPLY_CHIMERA_USING_CHIMERA_H_INCLUDED)
#define KRATOS_CUSTOM_APPLY_CHIMERA_USING_CHIMERA_H_INCLUDED

// System includes

#include <string>
#include <iostream>
#include <algorithm>

// External includes
#include "includes/kratos_flags.h"

#include "utilities/binbased_fast_point_locator.h"

// Project includes

#include "includes/define.h"
#include "processes/process.h"
#include "includes/kratos_flags.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "geometries/geometry_data.h"
#include "includes/variables.h"
#include "utilities/math_utils.h"
#include "includes/kratos_parameters.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "elements/distance_calculation_element_simplex.h"

// Application includes
#include "chimera_application_variables.h"
#include "custom_utilities/multipoint_constraint_data.hpp"
#include "custom_processes/custom_calculate_signed_distance_process.h"
#include "custom_hole_cutting_process.h"
#include "custom_utilities/vtk_output.hpp"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.

template <unsigned int TDim>

class ApplyChimeraProcess : public Process
{
  public:
	///@name Type Definitions
	///@{

	///@}
	///@name Pointer Definitions
	/// Pointer definition of ApplyChimeraProcess
	KRATOS_CLASS_POINTER_DEFINITION(ApplyChimeraProcess);
	typedef ProcessInfo::Pointer ProcessInfoPointerType;
	typedef typename BinBasedFastPointLocator<TDim>::Pointer BinBasedPointLocatorPointerType;
	typedef ModelPart::ConditionsContainerType ConditionsArrayType;
	typedef std::pair<unsigned int, unsigned int> SlavePairType;
	typedef Kratos::MpcData::MasterDofWeightMapType MasterDofWeightMapType;
	typedef ProcessInfo ProcessInfoType;
	typedef MpcData::Pointer MpcDataPointerType;
	typedef std::vector<MpcDataPointerType> *MpcDataPointerVectorType;
	typedef Dof<double> DofType;
	typedef std::vector<DofType> DofVectorType;
	typedef MpcData::VariableComponentType VariableComponentType;
	typedef unsigned int IndexType;
	typedef MpcData::VariableType VariableType;
	typedef Node<3> NodeType;

	///@}
	///@name Life Cycle
	///@{

	ApplyChimeraProcess(ModelPart &MainModelPart, Parameters rParameters) : Process(Flags()), mrMainModelPart(MainModelPart), m_parameters(rParameters)
	{

		Parameters default_parameters(R"(
            {
                "process_name":"apply_chimera_process",
                    "background": {
									"model_part_name":"GENERIC_background",
									"pressure_coupling":"",
									"type" : "nearest_element",
									"IsWeak" : true
                                  },
                	"patch" :     {
									"model_part_name":"GENERIC_patch",
									"pressure_coupling" : "one",
									"type" : "nearest_element",
									"IsWeak" : true
                                  },
					"pressure_coupling_node" : 0.0,			  
                    "patch_boundary_model_part_name":"GENERIC_patchBoundary",
                    "overlap_distance":0.045
            })");

		m_background_model_part_name = m_parameters["background"]["model_part_name"].GetString();
		m_patch_model_part_name = m_parameters["patch"]["model_part_name"].GetString();
		m_patch_boundary_model_part_name = m_parameters["patch_boundary_model_part_name"].GetString();
		m_type_patch = m_parameters["patch"]["type"].GetString();
		m_type_background = m_parameters["background"]["type"].GetString();
		m_overlap_distance = m_parameters["overlap_distance"].GetDouble();

		ModelPart &rBackgroundModelPart = mrMainModelPart.GetSubModelPart(m_background_model_part_name);
		ModelPart &rPatchModelPart = mrMainModelPart.GetSubModelPart(m_patch_model_part_name);

		ProcessInfoPointerType info = mrMainModelPart.pGetProcessInfo();
		if (info->GetValue(MPC_DATA_CONTAINER) == NULL)
			info->SetValue(MPC_DATA_CONTAINER, new std::vector<MpcDataPointerType>());

		this->pBinLocatorForBackground = BinBasedPointLocatorPointerType(new BinBasedFastPointLocator<TDim>(rBackgroundModelPart));
		this->pBinLocatorForPatch = BinBasedPointLocatorPointerType(new BinBasedFastPointLocator<TDim>(rPatchModelPart));
		this->pMpcPatch = MpcDataPointerType(new MpcData(m_type_patch));
		this->pMpcBackground = MpcDataPointerType(new MpcData(m_type_background));
		this->pHoleCuttingProcess = CustomHoleCuttingProcess::Pointer(new CustomHoleCuttingProcess());
		this->pCalculateDistanceProcess = typename CustomCalculateSignedDistanceProcess<TDim>::Pointer(new CustomCalculateSignedDistanceProcess<TDim>());

		this->pMpcPatch->SetName(m_patch_model_part_name);
		this->pMpcBackground->SetName(m_background_model_part_name);
		this->pMpcPatch->SetActive(true);
		this->pMpcBackground->SetActive(true);

		MpcDataPointerVectorType mpcDataVector = info->GetValue(MPC_DATA_CONTAINER);
		(*mpcDataVector).push_back(pMpcPatch);
		(*mpcDataVector).push_back(pMpcBackground);
	}

	/// Destructor.
	virtual ~ApplyChimeraProcess()
	{
		Clear();
	}

	///@}
	///@name Operators
	///@{

	void operator()()
	{
		Execute();
	}

	///@}
	///@name Operations
	///@{

	virtual void Execute()
	{
	}

	virtual void Clear()
	{
		pMpcBackground->Clear();
		pMpcPatch->Clear();
		std::cout << "Chimera process is cleared" << std::endl;
	}

	void ExecuteBeforeSolutionLoop() override
	{
	}

	void ExecuteInitializeSolutionStep() override
	{
		KRATOS_TRY;
		// Actual execution of the functionality of this class
		FormulateChimera();

		KRATOS_CATCH("");
	}

	void ExecuteFinalizeSolutionStep() override
	{

		Clear();
		//for multipatch
		for (ModelPart::ElementsContainerType::iterator it = mrMainModelPart.ElementsBegin(); it != mrMainModelPart.ElementsEnd(); ++it)
			it->Set(VISITED, false);
	}

	void ExecuteBeforeOutputStep() override
	{
	}

	void ExecuteAfterOutputStep() override
	{
	}

	void ExecuteFinalize() override
	{
	}

	void ApplyMpcConstraint(ModelPart &rBoundaryModelPart, BinBasedPointLocatorPointerType &pBinLocator, MpcDataPointerType pMpc, std::string pressure_coupling)
	{

		//loop over nodes and find the triangle in which it falls, than do interpolation
		array_1d<double, TDim + 1> N;
		const int max_results = 10000;
		typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);
		const int n_boundary_nodes = rBoundaryModelPart.Nodes().size();

		unsigned int counter = 0;

#pragma omp parallel for firstprivate(results, N)
		//MY NEW LOOP: reset the visited flag
		for (int i = 0; i < n_boundary_nodes; i++)

		{
			ModelPart::NodesContainerType::iterator iparticle = rBoundaryModelPart.NodesBegin() + i;
			Node<3>::Pointer p_boundary_node = *(iparticle.base());
			p_boundary_node->Set(VISITED, false);
		}

		for (int i = 0; i < n_boundary_nodes; i++)

		{
			ModelPart::NodesContainerType::iterator iparticle = rBoundaryModelPart.NodesBegin() + i;
			Node<3>::Pointer p_boundary_node = *(iparticle.base());

			typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();

			Element::Pointer pElement;

			bool is_found = false;
			is_found = pBinLocator->FindPointOnMesh(p_boundary_node->Coordinates(), N, pElement, result_begin, max_results);

			// Initialise the boundary nodes dofs to 0 at ever time steps

			p_boundary_node->GetDof(VELOCITY_X).GetSolutionStepValue(0) = 0.0;
			p_boundary_node->GetDof(VELOCITY_Y).GetSolutionStepValue(0) = 0.0;

			if (TDim == 3)
				p_boundary_node->GetDof(VELOCITY_Z).GetSolutionStepValue(0) = 0.0;

			if (pressure_coupling == "all")
				p_boundary_node->GetDof(PRESSURE).GetSolutionStepValue(0) = 0.0;

			if (is_found == true)
			{
				Geometry<Node<3>> &geom = pElement->GetGeometry();

				for (unsigned int i = 0; i < geom.size(); i++)
				{
					//Interpolation of velocity
					p_boundary_node->GetDof(VELOCITY_X).GetSolutionStepValue(0) += geom[i].GetDof(VELOCITY_X).GetSolutionStepValue(0) * N[i];
					p_boundary_node->GetDof(VELOCITY_Y).GetSolutionStepValue(0) += geom[i].GetDof(VELOCITY_Y).GetSolutionStepValue(0) * N[i];

					//Define master slave relation for velocity
					AddMasterSlaveRelationWithNodesAndVariableComponents(pMpc, geom[i], VELOCITY_X, *p_boundary_node, VELOCITY_X, N[i]);
					AddMasterSlaveRelationWithNodesAndVariableComponents(pMpc, geom[i], VELOCITY_Y, *p_boundary_node, VELOCITY_Y, N[i]);
					if (TDim == 3)
					{
						//Interpolation of velocity
						p_boundary_node->GetDof(VELOCITY_Z).GetSolutionStepValue(0) += geom[i].GetDof(VELOCITY_Z).GetSolutionStepValue(0) * N[i];

						//Define master slave relation for velocity
						AddMasterSlaveRelationWithNodesAndVariableComponents(pMpc, geom[i], VELOCITY_Z, *p_boundary_node, VELOCITY_Z, N[i]);
					}

					if (pressure_coupling == "all")
					{
						//Interpolation of pressure
						p_boundary_node->GetDof(PRESSURE).GetSolutionStepValue(0) += geom[i].GetDof(PRESSURE).GetSolutionStepValue(0) * N[i];

						//Defining master slave relation for pressure
						AddMasterSlaveRelationWithNodesAndVariable(pMpc, geom[i], PRESSURE, *p_boundary_node, PRESSURE, N[i]);

						counter++;
					}
				} // end of loop over host element nodes

			} // if (is_found = true)

			// Setting the buffer 1 same buffer 0
			p_boundary_node->GetDof(VELOCITY_X).GetSolutionStepValue(1) = p_boundary_node->GetDof(VELOCITY_X).GetSolutionStepValue(0);
			p_boundary_node->GetDof(VELOCITY_Y).GetSolutionStepValue(1) = p_boundary_node->GetDof(VELOCITY_Y).GetSolutionStepValue(0);
			if (TDim == 3)
				p_boundary_node->GetDof(VELOCITY_Z).GetSolutionStepValue(1) = p_boundary_node->GetDof(VELOCITY_Z).GetSolutionStepValue(0);

			if (pressure_coupling == "all")
				p_boundary_node->GetDof(PRESSURE).GetSolutionStepValue(1) = p_boundary_node->GetDof(PRESSURE).GetSolutionStepValue(0);

		} // end of loop over boundary nodes

		if (pressure_coupling == "one")
		{

			Node<3>::Pointer p_boundary_node;
			if (m_parameters["pressure_coupling_node"].GetDouble() == 0.0)
			{

				ModelPart::NodesContainerType::iterator iparticle = rBoundaryModelPart.NodesBegin();
				p_boundary_node = *(iparticle.base());
			}

			else
			{
				unsigned int node_num = m_parameters["pressure_coupling_node"].GetDouble();

				p_boundary_node = rBoundaryModelPart.pGetNode(node_num);
			}

			typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();

			Element::Pointer pElement;

			bool is_found = false;
			is_found = pBinLocator->FindPointOnMesh(p_boundary_node->Coordinates(), N, pElement, result_begin, max_results);

			//Initialsing pressure to 0 at every time steps
			p_boundary_node->GetDof(PRESSURE).GetSolutionStepValue(0) = 0.0;

			if (is_found == true)
			{
				Geometry<Node<3>> &geom = pElement->GetGeometry();
				for (unsigned int i = 0; i < geom.size(); i++)
				{
					// Interpolation of pressure
					p_boundary_node->GetDof(PRESSURE).GetSolutionStepValue(0) += geom[i].GetDof(PRESSURE).GetSolutionStepValue(0) * N[i];

					// Defining master slave relation for one pressure node
					AddMasterSlaveRelationWithNodesAndVariable(pMpc, geom[i], PRESSURE, *p_boundary_node, PRESSURE, N[i]);
					counter++;
				}
			}
			// Setting the buffer 1 same buffer 0
			p_boundary_node->GetDof(PRESSURE).GetSolutionStepValue(1) = p_boundary_node->GetDof(PRESSURE).GetSolutionStepValue(0);

			std::cout << "Coordinates of node that are pressure coupled" << std::endl;
			std::cout << p_boundary_node->X() << "," << p_boundary_node->Y() << "," << p_boundary_node->Z() << std::endl;

		} // end of if (pressure_coupling == "one")

		counter /= TDim + 1;

		std::cout << counter << " pressure nodes from " << rBoundaryModelPart.Name() << " is coupled" << std::endl;
	}

	void ApplyMpcConstraintConservative(ModelPart &rBoundaryModelPart, BinBasedPointLocatorPointerType &pBinLocator, MpcDataPointerType pMpc, std::string pressure_coupling)
	{

		double rtMinvR = 0;
		DofVectorType slaveDofVector;
		double R = 0;

		ApplyMpcConstraint(rBoundaryModelPart, pBinLocator, pMpc, pressure_coupling);
		std::vector<VariableComponentType> dofComponentVector = {VELOCITY_X, VELOCITY_Y, VELOCITY_Z};

		// Calculation of Rt*Minv*R and assignment of nodalnormals to the slave dofs
		for (ModelPart::NodesContainerType::iterator inode = rBoundaryModelPart.NodesBegin(); inode != rBoundaryModelPart.NodesEnd(); ++inode)
		{
			double Minode = inode->FastGetSolutionStepValue(NODAL_MASS);

			for (unsigned int i = 0; i < TDim; i++)
			{

				double rIdof = inode->FastGetSolutionStepValue(NORMAL)[i];
				DofType &slaveDOF = inode->GetDof(dofComponentVector[i]);
				AddNodalNormalSlaveRelationWithDofs(pMpc, inode->GetDof(dofComponentVector[i]), rIdof);
				slaveDofVector.push_back(slaveDOF);
				rtMinvR += (rIdof * rIdof) / Minode;
				R += rIdof;
			}

			AddNodalNormalSlaveRelationWithDofs(pMpc, inode->GetDof(PRESSURE), 0);
		}

		SetRtMinvR(pMpc, rtMinvR);
		std::cout << "RtMRinv of " << rBoundaryModelPart.Name() << rtMinvR << std::endl;

		CalculateConservativeCorrections(mrMainModelPart, pMpc);
		ApplyConservativeCorrections(mrMainModelPart, pMpc);
	}

	//Apply Chimera with or without overlap
	void FormulateChimera()
	{

		ModelPart &rBackgroundModelPart = mrMainModelPart.GetSubModelPart(m_background_model_part_name);
		ModelPart &rPatchBoundaryModelPart = mrMainModelPart.GetSubModelPart(m_patch_boundary_model_part_name);

		this->pBinLocatorForBackground->UpdateSearchDatabase();
		this->pBinLocatorForPatch->UpdateSearchDatabase();

		for (ModelPart::ElementsContainerType::iterator it = mrMainModelPart.ElementsBegin(); it != mrMainModelPart.ElementsEnd(); ++it)
		{
			if (!it->Is(VISITED)) //for multipatch
				it->Set(ACTIVE, true);
		}

		const double epsilon = 1e-12;
		if (m_overlap_distance < epsilon)
		{
			KRATOS_THROW_ERROR("", "Overlap distance should be a positive number \n", "");
		}

		if (m_overlap_distance > epsilon)

		{

			ModelPart::Pointer pHoleModelPart = ModelPart::Pointer(new ModelPart("HoleModelpart"));
			ModelPart::Pointer pHoleBoundaryModelPart = ModelPart::Pointer(new ModelPart("HoleBoundaryModelPart"));

			this->pCalculateDistanceProcess->CalculateSignedDistance(rBackgroundModelPart, rPatchBoundaryModelPart);
			this->pHoleCuttingProcess->CreateHoleAfterDistance(rBackgroundModelPart, *pHoleModelPart, *pHoleBoundaryModelPart, m_overlap_distance);

			//for multipatch
			for (ModelPart::ElementsContainerType::iterator it = pHoleModelPart->ElementsBegin(); it != pHoleModelPart->ElementsEnd(); ++it)
				it->Set(VISITED, true);

			CalculateNodalAreaAndNodalMass(rPatchBoundaryModelPart, 1);
			CalculateNodalAreaAndNodalMass(*pHoleBoundaryModelPart, -1);

			pMpcPatch->SetIsWeak(m_parameters["patch"]["IsWeak"].GetBool());
			pMpcBackground->SetIsWeak(m_parameters["background"]["IsWeak"].GetBool());

			std::string pr_coupling_patch = m_parameters["patch"]["pressure_coupling"].GetString();
			std::string pr_coupling_background = m_parameters["background"]["pressure_coupling"].GetString();

			if (m_type_patch == "nearest_element")
			{
				ApplyMpcConstraint(rPatchBoundaryModelPart, pBinLocatorForBackground, pMpcPatch, pr_coupling_patch);
				std::cout << "Patch boundary coupled with background" << std::endl;
			}

			else if (m_type_patch == "conservative")
			{

				ApplyMpcConstraintConservative(rPatchBoundaryModelPart, pBinLocatorForBackground, pMpcPatch, pr_coupling_patch);
				std::cout << "Patch boundary coupled with background using conservative approach" << std::endl;
			}

			if (m_type_background == "nearest_element")
			{
				ApplyMpcConstraint(*pHoleBoundaryModelPart, pBinLocatorForPatch, pMpcBackground, pr_coupling_background);
				std::cout << "HoleBoundary  coupled with patch" << std::endl;
			}

			else if (m_type_background == "conservative")
			{

				ApplyMpcConstraintConservative(*pHoleBoundaryModelPart, pBinLocatorForPatch, pMpcBackground, pr_coupling_background);
				std::cout << "HoleBoundary  coupled with patch using conservative approach" << std::endl;
			}
		}
	}

	void SetOverlapDistance(double distance)
	{

		this->m_overlap_distance = distance;
	}

	void SetType(std::string model_part_string, std::string type)
	{
		KRATOS_TRY
		{
			if ((type != "nearest_element") && (type != "conservative"))
				KRATOS_THROW_ERROR("", "Second argument should be either nearest_element or conservative \n", "");

			if (model_part_string == "patch")
			{
				this->m_type_patch = type;
				pMpcPatch->SetType(this->m_type_patch);
			}

			else if (model_part_string == "background")
			{
				this->m_type_background = type;
				pMpcBackground->SetType(this->m_type_background);
			}

			else

				KRATOS_THROW_ERROR("", "First argument should be either patch or background \n", "");
		}

		KRATOS_CATCH("")
	}

	void CalculateNodalAreaAndNodalMass(ModelPart &rBoundaryModelPart, int sign)
	{
		KRATOS_TRY
		ConditionsArrayType &rConditions = rBoundaryModelPart.Conditions();
		//resetting the normals and calculating centre point

		array_1d<double, 3> zero;
		array_1d<double, 3> centre;
		unsigned int n_nodes = rBoundaryModelPart.Nodes().size();

		zero[0] = 0.0;
		zero[1] = 0.0;
		zero[2] = 0.0;

		centre[0] = 0.0;
		centre[1] = 0.0;
		centre[2] = 0.0;

		for (ConditionsArrayType::iterator it = rConditions.begin();
			 it != rConditions.end(); it++)
		{
			Element::GeometryType &rNodes = it->GetGeometry();
			for (unsigned int in = 0; in < rNodes.size(); in++)
			{
				noalias((rNodes[in]).GetSolutionStepValue(NORMAL)) = zero;
			}
		}

		for (ModelPart::NodesContainerType::iterator inode = rBoundaryModelPart.NodesBegin(); inode != rBoundaryModelPart.NodesEnd(); ++inode)
		{

			centre += inode->Coordinates();
		}

		centre = centre / n_nodes;

		//calculating the normals and storing on the conditions
		array_1d<double, 3> An;
		if (TDim == 2)
		{
			for (ConditionsArrayType::iterator it = rConditions.begin();
				 it != rConditions.end(); it++)
			{
				if (it->GetGeometry().PointsNumber() == 2)
					CalculateNormal2D(it, An, centre, sign);
			}
		}
		else if (TDim == 3)
		{
			array_1d<double, 3> v1;
			array_1d<double, 3> v2;
			for (ConditionsArrayType::iterator it = rConditions.begin();
				 it != rConditions.end(); it++)
			{
				//calculate the normal on the given condition
				if (it->GetGeometry().PointsNumber() == 3)
					CalculateNormal3D(it, An, v1, v2, centre, sign);
			}
		}

		//adding the normals to the nodes
		for (ConditionsArrayType::iterator it = rConditions.begin();
			 it != rConditions.end(); it++)
		{
			Geometry<Node<3>> &pGeometry = (it)->GetGeometry();
			double coeff = 1.00 / pGeometry.size();
			const array_1d<double, 3> &normal = it->GetValue(NORMAL);
			double nodal_mass = MathUtils<double>::Norm3(normal);

			for (unsigned int i = 0; i < pGeometry.size(); i++)
			{
				noalias(pGeometry[i].FastGetSolutionStepValue(NORMAL)) += coeff * normal;
				pGeometry[i].FastGetSolutionStepValue(NODAL_MASS) += coeff * nodal_mass;
			}
		}

		KRATOS_CATCH("")
	}

	void CalculateNormal2D(ConditionsArrayType::iterator it, array_1d<double, 3> &An, array_1d<double, 3> &centre, int sign)
	{
		Geometry<Node<3>> &pGeometry = (it)->GetGeometry();
		array_1d<double, 3> rVector;

		An[0] = pGeometry[1].Y() - pGeometry[0].Y();
		An[1] = -(pGeometry[1].X() - pGeometry[0].X());
		An[2] = 0.00;

		rVector[0] = centre[0] - pGeometry[0].X();
		rVector[1] = centre[1] - pGeometry[0].Y();
		rVector[2] = 0.00;

		array_1d<double, 3> &normal = (it)->GetValue(NORMAL);
		noalias(normal) = An;

		if ((MathUtils<double>::Dot(An, rVector) > 0))
			normal = -1 * normal;

		normal = normal * sign;

		// 				(it)->SetValue(NORMAL,An);
	}

	void CalculateNormal3D(ConditionsArrayType::iterator it, array_1d<double, 3> &An,
						   array_1d<double, 3> &v1, array_1d<double, 3> &v2, array_1d<double, 3> &centre, int sign)
	{
		Geometry<Node<3>> &pGeometry = (it)->GetGeometry();
		array_1d<double, 3> rVector;

		v1[0] = pGeometry[1].X() - pGeometry[0].X();
		v1[1] = pGeometry[1].Y() - pGeometry[0].Y();
		v1[2] = pGeometry[1].Z() - pGeometry[0].Z();

		v2[0] = pGeometry[2].X() - pGeometry[0].X();
		v2[1] = pGeometry[2].Y() - pGeometry[0].Y();
		v2[2] = pGeometry[2].Z() - pGeometry[0].Z();

		rVector[0] = centre[0] - pGeometry[0].X();
		rVector[1] = centre[1] - pGeometry[0].Y();
		rVector[2] = centre[2] - pGeometry[0].Z();

		MathUtils<double>::CrossProduct(An, v1, v2);
		An *= 0.5;

		array_1d<double, 3> &normal = (it)->GetValue(NORMAL);
		noalias(normal) = An;

		if ((MathUtils<double>::Dot(An, rVector) > 0))
			normal = -1 * normal;

		normal = normal * sign;
		// 				noalias((it)->GetValue(NORMAL)) = An;
	}

	void CalculateConservativeCorrections(ModelPart &r_model_part, MpcDataPointerType pMpc)
	{

		double nodalMass;
		unsigned int slaveNodeId;
		unsigned int slaveNodeIdOther;
		unsigned int slaveDofKey;
		unsigned int slaveDofKeyOther;
		double slaveDofValueOther;
		SlavePairType slaveDofMap;
		SlavePairType slaveDofMapOther;
		double RtMinvR = pMpc->RtMinvR;
		double NodalNormalComponent;
		double NodalNormalComponentOther;
		//std::cout << " RtMinvR " << RtMinvR << std::endl;
		std::vector<double> VectorOfconstants;
		unsigned int slaveIndex = 0;

		for (auto slaveMasterDofMap : pMpc->mDofConstraints)
		{
			slaveDofMap = slaveMasterDofMap.first;
			slaveNodeId = slaveDofMap.first;
			slaveDofKey = slaveDofMap.second;
			Node<3> &slaveNode = r_model_part.Nodes()[slaveNodeId];
			Node<3>::DofsContainerType::iterator idof = slaveNode.GetDofs().find(slaveDofKey);
			nodalMass = slaveNode.FastGetSolutionStepValue(NODAL_MASS);
			NodalNormalComponent = pMpc->mSlaveDofToNodalNormalMap[slaveDofMap];
			VectorOfconstants.push_back(0.0);
			for (auto slaveMasterDofMapOther : pMpc->mDofConstraints)
			{

				slaveDofMapOther = slaveMasterDofMapOther.first;
				slaveNodeIdOther = slaveDofMapOther.first;
				slaveDofKeyOther = slaveDofMapOther.second;
				Node<3> &slaveNodeOther = r_model_part.Nodes()[slaveNodeIdOther];
				Node<3>::DofsContainerType::iterator idofOther = slaveNodeOther.GetDofs().find(slaveDofKeyOther);
				slaveDofValueOther = idofOther->GetSolutionStepValue();
				NodalNormalComponentOther = pMpc->mSlaveDofToNodalNormalMap[slaveDofMapOther];
				VectorOfconstants[slaveIndex] -= ((NodalNormalComponent * NodalNormalComponentOther) / (nodalMass * RtMinvR)) * slaveDofValueOther; // correction for zero flux

			} // slaveMasterDofMapOher loop

			slaveIndex++;

		} // slaveMasterDofMap loop

		slaveIndex = 0;

		//Applying correction in the slaveDofValue
		for (auto slaveMasterDofMap : pMpc->mDofConstraints)
		{

			slaveDofMap = slaveMasterDofMap.first;
			slaveNodeId = slaveDofMap.first;
			slaveDofKey = slaveDofMap.second;
			Node<3> &slaveNode = r_model_part.Nodes()[slaveNodeId];
			Node<3>::DofsContainerType::iterator idof = slaveNode.GetDofs().find(slaveDofKey);
			unsigned int slaveEquationId = idof->EquationId();

			pMpc->mSlaveEquationIdConstantsMap[slaveEquationId] = VectorOfconstants[slaveIndex];

			slaveIndex++;

		} // slaveMasterDofMap loop

		//std::cout << "Conservative correction to the velocity field applied" << std::endl;
		std::cout << "Conservative Correction of " << pMpc->mName << " is calculated" << std::endl;
	}

	void ApplyConservativeCorrections(ModelPart &r_model_part, MpcDataPointerType pMpc)
	{

		for (auto slaveMasterDofMap : pMpc->mDofConstraints)
		{
			SlavePairType slaveDofMap = slaveMasterDofMap.first;
			unsigned int slaveNodeId = slaveDofMap.first;
			unsigned int slaveDofKey = slaveDofMap.second;
			NodeType &node = r_model_part.Nodes()[slaveNodeId];
			Node<3>::DofsContainerType::iterator it = node.GetDofs().find(slaveDofKey);
			unsigned int slaveEquationId = it->EquationId();

			it->GetSolutionStepValue(0) += pMpc->mSlaveEquationIdConstantsMap[slaveEquationId];
			it->GetSolutionStepValue(1) += pMpc->mSlaveEquationIdConstantsMap[slaveEquationId];
		}

		std::cout << "Conservative Correction of " << pMpc->mName << " is applied" << std::endl;
	}

	// Functions which use two variable components

	/**
		Applies the MPC condition using two nodes, one as master and other as slave, and with the given weight
		@arg MasterNode 
        @arg MasterVariable 
        @arg SlaveNode 
        @arg SlaveVariable
        @arg weight
		*/
	void AddMasterSlaveRelationWithNodesAndVariableComponents(MpcDataPointerType pMpc, Node<3> &MasterNode, VariableComponentType &MasterVariable, Node<3> &SlaveNode, VariableComponentType &SlaveVariable, double weight, double constant = 0.0)
	{
		SlaveNode.Set(SLAVE);
		DofType &pointerSlaveDOF = SlaveNode.GetDof(SlaveVariable);
		DofType &pointerMasterDOF = MasterNode.GetDof(MasterVariable);
		AddMasterSlaveRelationWithDofs(pMpc, pointerSlaveDOF, pointerMasterDOF, weight, constant);
	}

	void AddMasterSlaveRelationWithNodeIdsAndVariableComponents(MpcDataPointerType pMpc, IndexType MasterNodeId, VariableComponentType &MasterVariable, IndexType SlaveNodeId, VariableComponentType &SlaveVariable, double weight, double constant = 0.0)
	{
		Node<3> &SlaveNode = mrMainModelPart.Nodes()[SlaveNodeId];
		Node<3> &MasterNode = mrMainModelPart.Nodes()[MasterNodeId];
		SlaveNode.Set(SLAVE);
		DofType &pointerSlaveDOF = SlaveNode.GetDof(SlaveVariable);
		DofType &pointerMasterDOF = MasterNode.GetDof(MasterVariable);
		AddMasterSlaveRelationWithDofs(pMpc, pointerSlaveDOF, pointerMasterDOF, weight, constant);
	}

	// Functions with use two variables
	void AddMasterSlaveRelationWithNodesAndVariable(MpcDataPointerType pMpc, Node<3> &MasterNode, VariableType &MasterVariable, Node<3> &SlaveNode, VariableType &SlaveVariable, double weight, double constant = 0.0)
	{
		SlaveNode.Set(SLAVE);
		DofType &pointerSlaveDOF = SlaveNode.GetDof(SlaveVariable);
		DofType &pointerMasterDOF = MasterNode.GetDof(MasterVariable);
		AddMasterSlaveRelationWithDofs(pMpc, pointerSlaveDOF, pointerMasterDOF, weight, constant);
	}

	void AddMasterSlaveRelationWithNodeIdsAndVariable(MpcDataPointerType pMpc, IndexType MasterNodeId, VariableType &MasterVariable, IndexType SlaveNodeId, VariableType &SlaveVariable, double weight, double constant = 0.0)
	{
		Node<3> &SlaveNode = mrMainModelPart.Nodes()[SlaveNodeId];
		Node<3> &MasterNode = mrMainModelPart.Nodes()[MasterNodeId];
		SlaveNode.Set(SLAVE);
		DofType &pointerSlaveDOF = SlaveNode.GetDof(SlaveVariable);
		DofType &pointerMasterDOF = MasterNode.GetDof(MasterVariable);
		AddMasterSlaveRelationWithDofs(pMpc, pointerSlaveDOF, pointerMasterDOF, weight, constant);
	}

	// Default functions
	/**
		Applies the MPC condition using DOFs, one as master and other as slave, and with the given weight
		@arg slaveDOF 
        @arg masterDOF 
        @arg weight
		*/
	void AddMasterSlaveRelationWithDofs(MpcDataPointerType pMpc, DofType slaveDOF, DofType masterDOF, double masterWeight, double constant = 0.0)
	{
		pMpc->AddConstraint(slaveDOF, masterDOF, masterWeight, constant);
	}

	void AddNodalNormalSlaveRelationWithDofs(MpcDataPointerType pMpc, DofType slaveDOF, double nodalNormalComponent = 0.0)
	{
		pMpc->AddNodalNormalToSlaveDof(slaveDOF, nodalNormalComponent);
	}

	/**
		Activates the constraint set or deactivates
		@arg isActive true/false
		*/
	void SetActive(bool isActive = true)
	{
		pMpcPatch->SetActive(isActive);
		pMpcBackground->SetActive(isActive);
	}

	void SetRtMinvR(MpcDataPointerType pMpc, double value)
	{

		pMpc->RtMinvR = value;
	}

	void PrintGIDMesh(ModelPart &rmodel_part)
	{
		std::ofstream myfile;
		myfile.open(rmodel_part.Name() + ".post.msh");
		myfile << "MESH \"leaves\" dimension 2 ElemType Line Nnode 2" << std::endl;
		myfile << "# color 96 96 96" << std::endl;
		myfile << "Coordinates" << std::endl;
		myfile << "# node number coordinate_x coordinate_y coordinate_z  " << std::endl;

		for (unsigned int i = 0; i < rmodel_part.Nodes().size(); i++)
		{
			ModelPart::NodesContainerType::iterator iparticle = rmodel_part.NodesBegin() + i;
			Node<3>::Pointer p_node = *(iparticle.base());
			myfile << p_node->Id() << "  " << p_node->Coordinates()[0] << "  " << p_node->Coordinates()[1] << "  " << p_node->Coordinates()[2] << std::endl;
		}

		myfile << "end coordinates" << std::endl;
		myfile << "elements" << std::endl;
		myfile << "# element node_1 node_2 material_number" << std::endl;

		for (ConditionsArrayType::iterator it = rmodel_part.Conditions().begin();
			 it != rmodel_part.Conditions().end(); it++)
		{

			myfile << it->Id() << "  ";
			for (unsigned int i = 0; i < it->GetGeometry().PointsNumber(); i++)
				myfile << (it->GetGeometry()[i]).Id() << "  ";

			myfile << std::endl;
		}

		myfile << "end elements" << std::endl;
	}

	virtual std::string Info() const
	{
		return "ApplyChimeraProcess";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const
	{
		rOStream << "ApplyChimeraProcess";
	}

	/// Print object's data.
	virtual void PrintData(std::ostream &rOStream) const
	{

		std::cout << "\nNumber of  patch slave nodes :: " << std::endl;
		pMpcPatch->GetInfo();

		std::cout << "\nNumber of  background slave nodes :: " << std::endl;
		pMpcBackground->GetInfo();
	}

	///@}
	///@name Friends
	///@{

	///@}

  protected:
	///@name Protected static Member Variables
	///@{

	///@}
	///@name Protected member Variables
	///@{

	///@}
	///@name Protected Operators
	///@{

	///@}
	///@name Protected Operations
	///@{

	///@}
	///@name Protected  Access
	///@{

	///@}
	///@name Protected Inquiry
	///@{

	///@}
	///@name Protected LifeCycle
	///@{

	///@}

  private:
	///@name Static Member Variables
	///@{

	///@}
	///@name Member Variables
	///@{

	//ModelPart &mrBackGroundModelPart;
	//ModelPart &mrPatchSurfaceModelPart;
	BinBasedPointLocatorPointerType pBinLocatorForBackground; // Template argument 3 stands for 3D case
	BinBasedPointLocatorPointerType pBinLocatorForPatch;
	MpcDataPointerType pMpcPatch;
	MpcDataPointerType pMpcBackground;
	CustomHoleCuttingProcess::Pointer pHoleCuttingProcess;
	typename CustomCalculateSignedDistanceProcess<TDim>::Pointer pCalculateDistanceProcess;
	ModelPart &mrMainModelPart;
	double m_overlap_distance;

	Parameters m_parameters;
	std::string m_background_model_part_name;
	std::string m_patch_boundary_model_part_name;
	std::string m_patch_model_part_name;
	std::string m_type_patch;
	std::string m_type_background;

	// epsilon
	//static const double epsilon;

	///@}
	///@name Private Operators
	///@{

	///@}
	///@name Private Operations
	///@{

	///@}
	///@name Private  Access
	///@{

	///@}
	///@name Private Inquiry
	///@{

	///@}
	///@name Un accessible methods
	///@{

	/// Assignment operator.
	ApplyChimeraProcess &operator=(ApplyChimeraProcess const &rOther);

	/// Copy constructor.
	//CustomExtractVariablesProcess(CustomExtractVariablesProcess const& rOther);

	///@}
}; // Class CustomExtractVariablesProcess

} // namespace Kratos.

#endif // KRATOS_CUSTOM_EXTRACT_VARIABLES_PROCESS_H_INCLUDED  defined
