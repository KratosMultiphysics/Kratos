// ==============================================================================
/*
 KratosTopologyOptimizationApplication
 A library based on:
 Kratos
 A General Purpose Software for Multi-Physics Finite Element Analysis
 (Released on march 05, 2007).

 Copyright (c) 2016: Daniel Baumgaertner
                     daniel.baumgaertner@tum.de
                     Chair of Structural Analysis
                     Technische Universitaet Muenchen
                     Arcisstrasse 21 80333 Munich, Germany

 Permission is hereby granted, free  of charge, to any person obtaining
 a  copy  of this  software  and  associated  documentation files  (the
 "Software"), to  deal in  the Software without  restriction, including
 without limitation  the rights to  use, copy, modify,  merge, publish,
 distribute,  sublicense and/or  sell copies  of the  Software,  and to
 permit persons to whom the Software  is furnished to do so, subject to
 the following condition:

 Distribution of this code for  any  commercial purpose  is permissible
 ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

 The  above  copyright  notice  and  this permission  notice  shall  be
 included in all copies or substantial portions of the Software.

 THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
 EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
 CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
 TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
 SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
//==============================================================================
//
//   Project Name:        KratosTopology                        $
//   Last modified by:	  $Author:   daniel.baumgaertner@tum.de $
// 						  $Co-Author: Octaviano Malfavón Farías $
//   Date:                $Date:                    August 2016 $
//   Revision:            $Revision:                        0.0 $
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
	ModelPart::Pointer mpAdjointModelPart;
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
