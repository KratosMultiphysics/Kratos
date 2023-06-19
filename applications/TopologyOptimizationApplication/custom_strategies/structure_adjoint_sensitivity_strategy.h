//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
//                   Octaviano Malfavón Farías
//                   Eric Gonzales
//					 Philipp Hofer
//					 Erich Wehrle
//
// ==============================================================================

#if !defined(KRATOS_STRUCTURE_ADJOINT_SENSITIVITY_STRATEGY_H_INCLUDED)
#define  KRATOS_STRUCTURE_ADJOINT_SENSITIVITY_STRATEGY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "custom_elements/small_displacement_simp_element.h"
#include "utilities/builtin_timer.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "topology_optimization_application.h"


namespace Kratos {

///@addtogroup TopologyOptimizationApplication
///@{

///@name Kratos Classes
///@{

/// Solution strategy to calculate the sensitivities.
/// Derives from the previously defined Solving Strategy

template    <class TSparseSpace, 
            class TDenseSpace, 
            class TLinearSolver
            >
class StructureAdjointSensitivityStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace>
{
public:

    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(StructureAdjointSensitivityStrategy);

    typedef SolvingStrategy<TSparseSpace,TDenseSpace> BaseType;

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


    ///virtual ~StructureAdjointSensitivityStrategy()
    ~StructureAdjointSensitivityStrategy()	override
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
        int i= 0;

        BuiltinTimer timer;
        const ProcessInfo& ConstProcessInfo= mr_structure_model_part.GetProcessInfo();

        block_for_each(mr_structure_model_part.Elements(), [&](Element& element_i)
        {
            element_i.Calculate(DCDX, Out, ConstProcessInfo);
            i++;
        });

        KRATOS_INFO("[TopOpt]") << "  Objective Function sensitivities computed  [ spent time =  " << timer.ElapsedSeconds() << " ] " << std::endl;

        KRATOS_CATCH("");
    }


    /// Computes DVDX sensitivities from the adjoint solution
    void ComputeVolumeFractionSensitivities()
    {
        KRATOS_TRY;

        double Out = 0.0;

        BuiltinTimer timer;
        const ProcessInfo& ConstProcessInfo= mr_structure_model_part.GetProcessInfo();

        block_for_each(mr_structure_model_part.Elements(), [&](Element& element_i)
        {
            element_i.Calculate(DVDX, Out, ConstProcessInfo);
        });

        KRATOS_INFO("[TopOpt]") << "  Volume fraction sensitivities computed     [ spent time =  " << timer.ElapsedSeconds() << " ] " << std::endl;

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
