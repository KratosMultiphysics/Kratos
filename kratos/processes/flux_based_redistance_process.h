//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Daniel Diez
//                   Pablo Becker
//

#pragma once


// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "elements/distance_calculation_flux_based_element.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
#include "solving_strategies/convergencecriterias/displacement_criteria.h"
#include "processes/generic_find_elements_neighbours_process.h"
#include "utilities/variable_utils.h"
#include "spatial_containers/spatial_containers.h" 
#include "modeler/connectivity_preserve_modeler.h"
#include "utilities/merge_variable_lists_utility.h"
#include "variational_distance_calculation_process.h"
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
/**takes a model part full of SIMPLICIAL ELEMENTS (triangles and tetras) and computes a field which resembles the filltime
     * This field is to be used to compute the optimization problem of gate location
     * 
    */
template <unsigned int TDim, unsigned int TNumNodes, class TSparseSpace, class TDenseSpace, class TLinearSolver>
class FluxBasedRedistanceProcess : public VariationalDistanceCalculationProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{

    typedef Scheme<TSparseSpace, TDenseSpace> SchemeType;
    typedef ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace,TLinearSolver> SolvingStrategyType;
    typedef typename ConvergenceCriteria< TSparseSpace, TDenseSpace >::Pointer ConvergenceCriteriaPointerType;


    ///@}
    ///@name Pointer Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(FluxBasedRedistanceProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    FluxBasedRedistanceProcess(
        ModelPart &rBaseModelPart,
        typename TLinearSolver::Pointer pLinearSolver,
        Parameters Settings) :
        VariationalDistanceCalculationProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>(rBaseModelPart)
    {
        KRATOS_TRY

        //by default we compute the pseudofilltime
        Parameters default_parameters(R"(
            {
                "echo_level"     : 0,
				"redistance_model_part_name" : "RedistanceCalculationPart"
            }  )");
        Settings.ValidateAndAssignDefaults(default_parameters);

		mAuxModelPartName = Settings["redistance_model_part_name"].GetString();

        //checking model part is correct;
        ValidateInput();
		//TODO:override this function so that non-simplex geometries can be used
		KRATOS_ERROR_IF(rBaseModelPart.ElementsBegin()->GetGeometry().size()!=TDim+1 )<< "Only simplex geometries supported" << std::endl;

        // Generate an auxilary model part and populate it by elements of type DistanceCalculationElementSimplex
        ReGenerateDistanceModelPart(rBaseModelPart);
        ModelPart& r_distance_model_part = mrModel.CreateModelPart( mAuxModelPartName );

        //compute density to get a Fo=1
        mDomainLength = CalculateDomainLength();
        r_distance_model_part.GetProcessInfo().SetValue(CHARACTERISTIC_LENGTH, mDomainLength);

		//finding neigh elems to see if elems have other elems in the upwind direction
        auto neigh_proc = GenericFindElementalNeighboursProcess(r_distance_model_part);
        neigh_proc.ExecuteInitialize();

        auto p_builder_solver = Kratos::make_shared<ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> >(pLinearSolver);
        InitializeSolutionStrategy(p_builder_solver);


        KRATOS_CATCH("")
    }

    /// Destructor.
    virtual ~FluxBasedRedistanceProcess()
    {
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    //TODO Charlie: This has to gone if you derive from process.
    /// Turn back information as a string.
    virtual std::string Info() const override
    {
        return "FluxBasedRedistanceProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "FluxBasedRedistanceProcess";
    }

    /**
    * @brief Computes the distance by solving two implicit systems of equations
    * the first system returns a potential field used to compute the pseudovelocity
    * and the second system uses the previous field to compute the pseudo filltime
    */
    void Execute() override
    {
        KRATOS_TRY;

		if(mDistancePartIsInitialized == false){
            ReGenerateDistanceModelPart(mrBaseModelPart);
        }

		this->ResetDistanceValuesAndFixities();

        ModelPart& r_distance_model_part = mrModel.GetModelPart( mAuxModelPartName );

        KRATOS_INFO("FluxBasedRedistanceProcess") << "Solving first redistance step\n";
		r_distance_model_part.GetProcessInfo().SetValue(REDISTANCE_STEP, 1);
        mpSolvingStrategy->Solve();

        //step 1.2: compute velocity using the gradient of the potential field:
        ComputeVelocities();

        //step 2: compute distance
        r_distance_model_part.GetProcessInfo().SetValue(REDISTANCE_STEP, 2);
        KRATOS_INFO("FluxBasedRedistanceProcess") << "Solving second redistance step\n";
        mpSolvingStrategy->Solve();

        VariableUtils().ApplyFixity(DISTANCE, false, r_distance_model_part.Nodes());

        KRATOS_CATCH("")
    }


    ///@}

protected:
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

    double mDomainLength;

    ///@}
    ///@name Protected Operators
    ///@{

    /// Constructor without linear solver for derived classes
    FluxBasedRedistanceProcess(
        ModelPart &rBaseModelPart) : mrBaseModelPart(rBaseModelPart), mrModel(rBaseModelPart.GetModel())
    {
        mDistancePartIsInitialized = false;
    }

    ///@}
    ///@name Protected Operations
    ///@{

    virtual void ReGenerateDistanceModelPart(ModelPart &rBaseModelPart) override
    {
        KRATOS_TRY

        if(mrModel.HasModelPart( mAuxModelPartName ))
            mrModel.DeleteModelPart( mAuxModelPartName );

        // Ensure that the nodes have distance as a DOF
        VariableUtils().AddDof<Variable<double> >(DISTANCE, rBaseModelPart);

        // Generate
        ModelPart& r_distance_model_part = mrModel.CreateModelPart( mAuxModelPartName );

        //varibles
        MergeVariableListsUtility().Merge(r_distance_model_part, rBaseModelPart);  

        // Generating the model part
        Element::Pointer p_element = Kratos::make_intrusive<DistanceCalculationFluxBasedElement<TDim,TNumNodes>>();
        
		ConnectivityPreserveModeler modeler;
        modeler.GenerateModelPart(rBaseModelPart, r_distance_model_part, *p_element);

        mDistancePartIsInitialized = true;

        KRATOS_CATCH("")
    }

    /**
         * @brief Computes the nodal velocities 
         * The direction of the velocity is defined as the gradient of the DISTANCE
         * While the modulus of the velocity is proportional to its thickness.
         * NOTE: It is not a real velocity since it does not take into account flowrate or flowing section
         * It is only meant to provide a faster moving flow front in thicker regions.
         */
    void ComputeVelocities()
    {
        KRATOS_TRY

        ModelPart& r_distance_model_part = mrModel.GetModelPart( mAuxModelPartName );

        const ProcessInfo &rCurrentProcessInfo = r_distance_model_part.GetProcessInfo();

        //not using variable utils to do the two tasks in the same loop
        block_for_each(r_distance_model_part.Nodes(), [&](Node<3> &rNode) {
            rNode.SetValue(NODAL_VOLUME, 0.0);
            rNode.SetValue(POTENTIAL_GRADIENT, ZeroVector(3));
        });

        block_for_each(r_distance_model_part.Elements(), [&](ModelPart::ElementType &rElement) {
            rElement.AddExplicitContribution(rCurrentProcessInfo);
        });

        //not using variable utils to do the two tasks in the same loop
        block_for_each(r_distance_model_part.Nodes(), [&](Node<3> &rNode) {
            array_1d<double, 3> &vel = rNode.GetValue(POTENTIAL_GRADIENT);
            vel /= rNode.GetValue(NODAL_VOLUME);
        });

        KRATOS_CATCH("")
    }

    double CalculateDomainLength()
    {
        typedef CombinedReduction<MinReduction<double>,
                                  MinReduction<double>,
                                  MinReduction<double>,
                                  MaxReduction<double>,
                                  MaxReduction<double>,
                                  MaxReduction<double>>
            MultipleReduction;
        double x_min, y_min, z_min, x_max, y_max, z_max;

        ModelPart& r_distance_model_part = mrModel.CreateModelPart( mAuxModelPartName );

        std::tie(
            x_min,
            y_min,
            z_min,
            x_max,
            y_max,
            z_max) = block_for_each<MultipleReduction>(r_distance_model_part.Nodes(), [&](Node<3> &rNode) {
            const array_1d<double, 3> coord = rNode.Coordinates();
            return std::make_tuple(
                coord[0],
                coord[1],
                coord[2],
                coord[0],
                coord[1],
                coord[2]);
        });

        return std::sqrt((x_max - x_min) * (x_max - x_min) + (y_max - y_min) * (y_max - y_min) + (z_max - z_min) * (z_max - z_min));
    }

    ///@}



private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

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
    FluxBasedRedistanceProcess &operator=(FluxBasedRedistanceProcess const &rOther);

    ///@}

}; // Class FluxBasedRedistanceProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// Input stream function

///@}

} // namespace Kratos.
