// ==============================================================================
//  KratosTopologyOptimizationApplication
//
//  License:         BSD License
//                   license: TopologyOptimizationApplication/license.txt
//
//  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
//                   Octaviano Malfavón Farías
//                   Eric Gonzales
//
// ==============================================================================

#if !defined(KRATOS_STRUCTURE_ADJOINT_SENSITIVITY_STRATEGY_H_INCLUDED)
#define  KRATOS_STRUCTURE_ADJOINT_SENSITIVITY_STRATEGY_H_INCLUDED

// System includes
#include <iostream>
#include <string>
#include <algorithm>

// External includes
#include <boost/python.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "../custom_elements/small_displacement_simp_element.hpp"
// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

// Application includes
#include "topology_optimization_application.h"


namespace Kratos {

///@addtogroup TopologyOptimizationApplication
///@{

///@name Kratos Classes
///@{

/// Solution strategy to calculate the sensitivities.
/// Derives from the previously defined Solving Strategy

template<class TSparseSpace, class TDenseSpace, class TLinearSolver>
class StructureAdjointSensitivityStrategy: public SolvingStrategy<TSparseSpace, TDenseSpace,TLinearSolver>
{
public:

	///@name Type Definitions
	///@{

	KRATOS_CLASS_POINTER_DEFINITION(StructureAdjointSensitivityStrategy);

	typedef SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver> BaseType;
	typedef typename Scheme<TSparseSpace,TDenseSpace>::Pointer SchemePointerType;
	typedef typename BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>::Pointer BuilderAndSolverPointerType;

	///@}
	///@name Life Cycle
	///@{

	StructureAdjointSensitivityStrategy( ModelPart& rStructureModelPart,
			typename TLinearSolver::Pointer pNewLinearSolver,
			const int dimension = 3)
	: BaseType(rStructureModelPart),
	  mr_structure_model_part(rStructureModelPart),
	  m_dimension(dimension)
	{}


	virtual ~StructureAdjointSensitivityStrategy()
	{}

	///@}
	///@name Operations
	///@{

	// ---------------------------------------------------------------------------------------------------------------------------------------------
	// --------------------------------- COMPUTE SENSITIVITIES  ------------------------------------------------------------------------------------
	// ---------------------------------------------------------------------------------------------------------------------------------------------

	/// Computes DCDX sensitivities from the adjoint solution
	void ComputeStrainEnergySensitivities()
	{
		KRATOS_TRY;

		double Out = 0.0;

		clock_t begin = clock();

		for ( ModelPart::ElementIterator element_i = mr_structure_model_part.ElementsBegin(); element_i!= mr_structure_model_part.ElementsEnd();
				element_i++ )
		{
			const ProcessInfo& ConstProcessInfo= mr_structure_model_part.GetProcessInfo();
			element_i->Calculate(DCDX, Out, ConstProcessInfo);
		}

		clock_t end = clock();
		std::cout << "  Objective Function sensitivities computed  [ spent time =  " << double(end - begin) / CLOCKS_PER_SEC << " ] " << std::endl;

		KRATOS_CATCH("");
	}


	/// Computes DVDX sensitivities from the adjoint solution
	void ComputeVolumeFractionSensitivities()
	{
		KRATOS_TRY;

		double Out = 0.0;

		clock_t begin = clock();

		for ( ModelPart::ElementIterator element_i = mr_structure_model_part.ElementsBegin(); element_i!= mr_structure_model_part.ElementsEnd();
				element_i++ )
		{
			const ProcessInfo& ConstProcessInfo= mr_structure_model_part.GetProcessInfo();
			element_i->Calculate(DVDX, Out, ConstProcessInfo);
		}

		clock_t end = clock();
		std::cout << "  Volume fraction sensitivities computed     [ spent time =  " << double(end - begin) / CLOCKS_PER_SEC << " ] " << std::endl;

		KRATOS_CATCH("");
	}

	///@}

private:

	///@name Member Variables
	///@{

	ModelPart& mr_structure_model_part;
	ModelPart* mpAdjointModelPart;
	typename BaseType::Pointer mpStrategy;
	int m_dimension;

	///@}
	///@name Private Operations
	///@{

	///@}
}; // class StructureAdjointSensitivityStrategy

///@} // Kratos classes
///@} // AdjointStructureApplication group
}

#endif	/* KRATOS_STRUCTURE_ADJOINT_SENSITIVITY_STRATEGY_H_INCLUDED */
