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

#if !defined(KRATOS_CUSTOM_CALCULATE_SIGNED_DISTANCE_PROCESS_H_INCLUDED)
#define KRATOS_CUSTOM_CALCULATE_SIGNED_DISTANCE_PROCESS_H_INCLUDED

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

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "elements/distance_calculation_element_simplex.h"
#include "utilities/parallel_levelset_distance_calculator.h"
#include "custom_processes/calculate_signed_distance_to_2d_condition_skin_process.h"
#include "processes/calculate_signed_distance_to_3d_condition_skin_process.h"

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
template <std::size_t TDim>
class CustomCalculateSignedDistanceProcess
{
  public:
	///@name Type Definitions
	///@{

	///@}
	///@name Pointer Definitions
	/// Pointer definition of CustomExtractVariablesProcess
	KRATOS_CLASS_POINTER_DEFINITION(CustomCalculateSignedDistanceProcess);

	///@}
	///@name Life Cycle
	///@{

	CustomCalculateSignedDistanceProcess()
	{
		this->pBinLocator = NULL;
		this->pDistanceCalculator = NULL;

		p2DSignedDistanceCalculator = NULL;
		p3DSignedDistanceCalculator = NULL;
	}

	/// Destructor.
	virtual ~CustomCalculateSignedDistanceProcess()
	{
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
	}

	void ExtractDistance(ModelPart &fromPatchModelPart, ModelPart &toBackgroundModelPart, ModelPart &patchBoundaryModelPart)
	{
		this->pDistanceCalculator = typename ParallelDistanceCalculator<TDim>::Pointer(new ParallelDistanceCalculator<TDim>());
		this->pBinLocator = typename BinBasedFastPointLocator<TDim>::Pointer(new BinBasedFastPointLocator<TDim>(fromPatchModelPart));

		const int n_patch_boundary_nodes = patchBoundaryModelPart.Nodes().size();
		//Set the boundary node distance to 0 and is_visited to 1.0
		for (int i = 0; i < n_patch_boundary_nodes; i++)
		{

			ModelPart::NodesContainerType::iterator it = patchBoundaryModelPart.NodesBegin() + i;
			double &is_visited = it->GetValue(IS_VISITED);
			double &distance = it->FastGetSolutionStepValue(DISTANCE);

			//Set the IS_VISITED to 1.0 and DISTANCE to 0.0
			is_visited = 1.0;
			distance = 0.0;
		}

		std::size_t max_level = 100;
		double max_distance = 1000;

		this->pDistanceCalculator->CalculateDistancesLagrangianSurface(fromPatchModelPart, DISTANCE, NODAL_AREA, max_level, max_distance);

		//For signed distance
		/*std::size_t n_patch_nodes = fromPatchModelPart.Nodes().size();
		for(int i = 0; i < n_patch_nodes; i++ ) {

			ModelPart::NodesContainerType::iterator it=patchBoundaryModelPart.NodesBegin()+i;

            double& distance = it->FastGetSolutionStepValue(DISTANCE);
			distance = -distance;


		}*/

		/*
		  This part of the code below is adapted from "MappingPressureToStructure" function of class CalculateSignedDistanceTo3DSkinProcess
		 */

		{
			//loop over nodes and find the tetra/triangle in which it falls, than do interpolation
			array_1d<double, TDim + 1> N;
			const int max_results = 10000;
			typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);
			const int n_background_nodes = toBackgroundModelPart.Nodes().size();

#pragma omp parallel for firstprivate(results, N)
			//MY NEW LOOP: reset the visited flag
			for (int i = 0; i < n_background_nodes; i++)
			{
				ModelPart::NodesContainerType::iterator iparticle = toBackgroundModelPart.NodesBegin() + i;
				typename Node<3>::Pointer p_background_node = *(iparticle.base());
				p_background_node->Set(VISITED, false);
			}
			for (int i = 0; i < n_background_nodes; i++)
			{
				ModelPart::NodesContainerType::iterator iparticle = toBackgroundModelPart.NodesBegin() + i;
				typename Node<3>::Pointer p_background_node = *(iparticle.base());
				typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();
				Element::Pointer pElement;

				bool is_found = this->pBinLocator->FindPointOnMesh(p_background_node->Coordinates(), N, pElement, result_begin, max_results);

				if (is_found == true)
				{

					Geometry<Node<3>> &geom = pElement->GetGeometry();

					double &nodalDistance = p_background_node->FastGetSolutionStepValue(DISTANCE);

					// Do mapping
					MapDistancetoNode((*p_background_node), geom, nodalDistance, N);
					//nodalVariableValues is passed by reference

					p_background_node->Set(VISITED);
				}
			}
		}
	}

	void CalculateSignedDistanceOnModelPart(ModelPart &fromPatchModelPart,ModelPart &patchBoundaryModelPart)
	{

		this->pDistanceCalculator = typename ParallelDistanceCalculator<TDim>::Pointer(new ParallelDistanceCalculator<TDim>());

		const int n_patch_boundary_nodes = patchBoundaryModelPart.Nodes().size();
		//Set the boundary node distance to 0 and is_visited to 1.0
		for (int i = 0; i < n_patch_boundary_nodes; i++)
		{

			ModelPart::NodesContainerType::iterator it = patchBoundaryModelPart.NodesBegin() + i;
			double &is_visited = it->GetValue(IS_VISITED);
			double &distance = it->FastGetSolutionStepValue(DISTANCE);

			//Set the IS_VISITED to 1.0 and DISTANCE to 0.0
			is_visited = 1.0;
			distance = 0.0;
		}

		std::size_t max_level = 100;
		double max_distance = 1000;

		this->pDistanceCalculator->CalculateDistancesLagrangianSurface(fromPatchModelPart, DISTANCE, NODAL_AREA, max_level, max_distance);

		//For signed distance
		std::size_t n_patch_nodes = fromPatchModelPart.Nodes().size();
		for (int i = 0; i < n_patch_nodes; i++)
		{

			ModelPart::NodesContainerType::iterator it = fromPatchModelPart.NodesBegin() + i;

			double &distance = it->FastGetSolutionStepValue(DISTANCE);
			distance = -distance;
		}

	}

	void CalculateSignedDistance(ModelPart &toBackgroundModelPart, ModelPart &patchBoundaryModelPart)
	{

		if (TDim == 2)
		{
			// Implemented in the custom_processes
			p2DSignedDistanceCalculator = CalculateSignedDistanceTo2DConditionSkinProcess::Pointer(new CalculateSignedDistanceTo2DConditionSkinProcess(patchBoundaryModelPart, toBackgroundModelPart));
			p2DSignedDistanceCalculator->Execute();
		}

		if (TDim == 3)
		{
			// From Core ?? Improve performance and algorithm based on CalculateSignedDistanceToSkinProcess
			p3DSignedDistanceCalculator = CalculateSignedDistanceTo3DConditionSkinProcess::Pointer(new CalculateSignedDistanceTo3DConditionSkinProcess(patchBoundaryModelPart, toBackgroundModelPart));


			p3DSignedDistanceCalculator->Execute();

		}

		std::size_t max_level = 100;
		double max_distance = 200;
		//CorrectSign(toBackgroundModelPart);

		pDistanceCalculator->CalculateDistances(toBackgroundModelPart, DISTANCE, NODAL_AREA, max_level, max_distance);
	}

	void MapDistancetoNode(const Node<3> &pNode,
						   Geometry<Node<3>> &geom, double &rVariable, const array_1d<double, TDim + 1> &N)
	{
		rVariable = 0;
		for (int i = 0; i < geom.size(); i++)
		{

			rVariable += N[i] * geom[i].FastGetSolutionStepValue(DISTANCE);
		}
	}

 	void CorrectSign(ModelPart &rModelPart)
	 {
		for (ModelPart::ElementsContainerType::iterator it = rModelPart.ElementsBegin(); it != rModelPart.ElementsEnd(); ++it)
		{
			bool is_split = it->GetValue(SPLIT_ELEMENT);
            if (is_split == false)
			{
				double elementDistance = 0.0;
				Geometry<Node<3>> &geom = it->GetGeometry();
				int positive =0;
				int negative =0;
				for (int j = 0; j < geom.size(); j++)
				{
					bool newsign = it->GetGeometry()[j].FastGetSolutionStepValue(DISTANCE)<0 ? 0:1;
					if(newsign == true )
						positive ++;
					else
						negative++;
				}
				for (int j = 0; j < geom.size(); j++)
				{
					elementDistance = fabs (it->GetGeometry()[j].FastGetSolutionStepValue(DISTANCE));
					if(positive > negative)
						it->GetGeometry()[j].FastGetSolutionStepValue(DISTANCE)= elementDistance;
					else
						it->GetGeometry()[j].FastGetSolutionStepValue(DISTANCE)= -elementDistance;
				}
			}
		}
	 }


	///@}
	///@name Access
	///@{

	///@}
	///@name Inquiry
	///@{

	///@}
	///@name Input and output
	///@{

	/// Turn back information as a string.
	virtual std::string Info() const
	{
		return "CustomCalculateSignedDistanceProcess";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const
	{
		rOStream << "CustomCalculateSignedDistanceProcess";
	}

	/// Print object's data.
	virtual void PrintData(std::ostream &rOStream) const
	{
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
	typename ParallelDistanceCalculator<TDim>::Pointer pDistanceCalculator;
	//ModelPart &mrBackGroundModelPart;
	//ModelPart &mrPatchSurfaceModelPart;
	typename BinBasedFastPointLocator<TDim>::Pointer pBinLocator; // Template argument 3 stands for 3D case
	CalculateSignedDistanceTo2DConditionSkinProcess::Pointer p2DSignedDistanceCalculator;
	CalculateSignedDistanceTo3DConditionSkinProcess::Pointer p3DSignedDistanceCalculator;

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
	CustomCalculateSignedDistanceProcess<TDim> &operator=(CustomCalculateSignedDistanceProcess<TDim> const &rOther);

	/// Copy constructor.
	//CustomCalculateSignedDistanceProcess(CustomCalculateSignedDistanceProcess<TDim> const& rOther);

	///@}

}; // Class CustomExtractVariablesProcess

} // namespace Kratos.

#endif // KRATOS_CUSTOM_EXTRACT_VARIABLES_PROCESS_H_INCLUDED  defined
