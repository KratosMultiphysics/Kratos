//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Ruben Zorrilla
//
//


#pragma once

// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes

// Project includes
#include "includes/define.h"
#include "containers/model.h"
#include "includes/kratos_flags.h"
#include "elements/distance_calculation_element_simplex.h"
#include "linear_solvers/linear_solver.h"
#include "processes/process.h"
#include "modeler/connectivity_preserve_modeler.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "utilities/variable_utils.h"

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
/**takes a model part full of SIMPLICIAL ELEMENTS (triangles and tetras) and recomputes a signed distance function
mantaining as much as possible the position of the zero of the function prior to the call.

This is achieved by minimizing the function  ( 1 - norm( gradient( distance ) )**2
with the restriction that "distance" is a finite elment function
*/
template< unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver >
class VariationalDistanceCalculationProcess : public Process
{
public:

    KRATOS_DEFINE_LOCAL_FLAG(PERFORM_STEP1);
    KRATOS_DEFINE_LOCAL_FLAG(DO_EXPENSIVE_CHECKS);
    KRATOS_DEFINE_LOCAL_FLAG(CALCULATE_EXACT_DISTANCES_TO_PLANE);

    ///@name Type Definitions
    ///@{

    typedef Scheme< TSparseSpace,  TDenseSpace > SchemeType;
    typedef typename SchemeType::Pointer SchemePointerType;
    typedef typename BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>::Pointer BuilderSolverPointerType;
    typedef ImplicitSolvingStrategy< TSparseSpace, TDenseSpace, TLinearSolver > SolvingStrategyType;

    ///@}
    ///@name Pointer Definitions

    /// Pointer definition of VariationalDistanceCalculationProcess
    KRATOS_CLASS_POINTER_DEFINITION(VariationalDistanceCalculationProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /**This process recomputed the distance function mantaining the zero of the existing distance distribution
     * for this reason the DISTANCE should be initialized to values distinct from zero in at least some portions of the domain
     * alternatively, the DISTANCE shall be fixed to zero at least on some nodes, and the process will compute a positive distance
     * respecting that zero
     * @param base_model_parr - is the model part on the top of which the calculation will be performed
     * @param plinear_solver  - linear solver to be used internally
     * @max_iterations        - maximum number of iteration to be employed in the nonlinear optimization process.
     *                        - can also be set to 0 if a (very) rough approximation is enough
     *
     * EXAMPLE OF USAGE FROM PYTHON:
     *
     class distance_linear_solver_settings:
         solver_type = "AMGCL"
         tolerance = 1E-3
         max_iteration = 200
         scaling = False
         krylov_type = "CG"
         smoother_type = "SPAI0"
         verbosity = 0

     import linear_solver_factory
     distance_linear_solver = linear_solver_factory.ConstructSolver(distance_linear_solver_settings)

     max_iterations=1
     distance_calculator = VariationalDistanceCalculationProcess2D(fluid_model_part, distance_linear_solver, max_iterations)
     distance_calculator.Execute()
     */

    /**
     * @brief Construct a new Variational Distance Calculation Process object
     * This process recomputes the distance function mantaining the zero of the existing distance distribution, stored in DISTANCE
     * For this reason the DISTANCE should be initialized to values distinct from zero in at least some portions of the domain
     * Alternatively, the DISTANCE shall be fixed to zero at least on some nodes, and the process will compute a positive distance
     * respecting that zero
     * @param rModel The model container
     * @param pLinearSolver Pointer to the linear solver to be used internally
     * @param ThisParameters Process settings to be validated
     */
    VariationalDistanceCalculationProcess(
        Model& rModel,
        typename TLinearSolver::Pointer pLinearSolver,
        Parameters ThisParameters)
        : VariationalDistanceCalculationProcess(
            rModel,
            pLinearSolver,
            Kratos::make_shared<ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> >(pLinearSolver),
            ThisParameters)
    {
    }

    /**
     * @brief Construct a new Variational Distance Calculation Process object
     * This process recomputes the distance function mantaining the zero of the existing distance distribution, stored in DISTANCE
     * For this reason the DISTANCE should be initialized to values distinct from zero in at least some portions of the domain
     * Alternatively, the DISTANCE shall be fixed to zero at least on some nodes, and the process will compute a positive distance
     * respecting that zero
     * Note that this constructor with custom builder and solver is to be used in the trilinos version, since the trilinos builder and
     * solver needs additional data (the EpetraComm).
     * @param rModel The model container
     * @param pLinearSolver Pointer to the linear solver to be used internally
     * @param pBuilderAndSolver Pointer to the custom linear and solver to be used internally
     * @param ThisParameters Process settings to be validated
     */
    VariationalDistanceCalculationProcess(
        Model& rModel,
        typename TLinearSolver::Pointer pLinearSolver,
        BuilderSolverPointerType pBuilderAndSolver,
        Parameters ThisParameters)
        : mDistancePartIsInitialized(false)
        , mrModel(rModel)
        , mrBaseModelPart(rModel.GetModelPart(ThisParameters["model_part_name"].GetString()))
    {
        // Check and assign settings
        ThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());
        if (ThisParameters["calculate_exact_distances_to_plane"].GetBool()) {
            mOptions = CALCULATE_EXACT_DISTANCES_TO_PLANE;
        } else {
            mOptions = CALCULATE_EXACT_DISTANCES_TO_PLANE.AsFalse();
        }
        mMaxIterations = ThisParameters["max_iterations"].GetInt();
        mAuxModelPartName = ThisParameters["auxiliary_model_part_name"].GetString();
        mCoefficient1 = ThisParameters["variational_redistance_coefficient_1"].GetDouble();
        mCoefficient2 = ThisParameters["variational_redistance_coefficient_2"].GetDouble();

        // Check that the process input is valid
        ValidateInput();

        // Generate an auxilary model part and populate it by elements of type DistanceCalculationElementSimplex
        ReGenerateDistanceModelPart(mrBaseModelPart);

        // Set provided builder and solver and initialize the solution strategy
        InitializeSolutionStrategy(pBuilderAndSolver);
        mpSolvingStrategy->SetEchoLevel(ThisParameters["echo_level"].GetInt());
    }

    VariationalDistanceCalculationProcess(
        ModelPart& rBaseModelPart,
        typename TLinearSolver::Pointer pLinearSolver,
        unsigned int MaxIterations = 10,
        Flags Options = CALCULATE_EXACT_DISTANCES_TO_PLANE.AsFalse(),
        std::string AuxPartName = "RedistanceCalculationPart",
        double Coefficient1 = 0.01,
        double Coefficient2 = 0.1)
    :
        mDistancePartIsInitialized(false),
        mMaxIterations(MaxIterations),
        mrModel( rBaseModelPart.GetModel() ),
        mrBaseModelPart (rBaseModelPart),
        mOptions( Options ),
        mAuxModelPartName( AuxPartName ),
        mCoefficient1(Coefficient1),
        mCoefficient2(Coefficient2)
    {
        KRATOS_TRY

        KRATOS_WARNING("VariationalDistanceCalculationProcess") << "This constructor is deprecated, please use the Parameters-based one." << std::endl;

        ValidateInput();

        // Generate an auxilary model part and populate it by elements of type DistanceCalculationElementSimplex
        ReGenerateDistanceModelPart(rBaseModelPart);

        auto p_builder_solver = Kratos::make_shared<ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> >(pLinearSolver);

        InitializeSolutionStrategy(p_builder_solver);

        KRATOS_CATCH("")
    }

    /// Constructor with custom Builder And Solver
    /** To be used in the trilinos version, since the trilinos builder and
     *  solver needs additional data (the EpetraComm).
     *  @param rBaseModelPart Reference ModelPart for distance calculation.
     *  @param pLinearSolver Linear solver for the distance system.
     *  @param MaxIterations Maximum number of non-linear optimization iterations.
     *  @param Options Configuration flags for the procedure.
     *  @param AuxPartName Name to be used for the internal distance calculation ModelPart.
     */
    VariationalDistanceCalculationProcess(
        ModelPart& rBaseModelPart,
        typename TLinearSolver::Pointer pLinearSolver,
        BuilderSolverPointerType pBuilderAndSolver,
        unsigned int MaxIterations = 10,
        Flags Options = CALCULATE_EXACT_DISTANCES_TO_PLANE.AsFalse(),
        std::string AuxPartName = "RedistanceCalculationPart",
        double Coefficient1 = 0.01,
        double Coefficient2 = 0.1)
    :
        mDistancePartIsInitialized(false),
        mMaxIterations(MaxIterations),
        mrModel( rBaseModelPart.GetModel() ),
        mrBaseModelPart (rBaseModelPart),
        mOptions( Options ),
        mAuxModelPartName( AuxPartName ),
        mCoefficient1(Coefficient1),
        mCoefficient2(Coefficient2)
    {
        KRATOS_TRY

        KRATOS_WARNING("VariationalDistanceCalculationProcess") << "This constructor is deprecated, please use the Parameters-based one." << std::endl;

        ValidateInput();

        // Generate an auxilary model part and populate it by elements of type DistanceCalculationElementSimplex
        ReGenerateDistanceModelPart(rBaseModelPart);

        InitializeSolutionStrategy(pBuilderAndSolver);

        KRATOS_CATCH("")
    }

    /// Destructor.
    ~VariationalDistanceCalculationProcess() override
    {
        Clear();
    };

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

    void Execute() override
    {
        KRATOS_TRY;

        if(mDistancePartIsInitialized == false){
            ReGenerateDistanceModelPart(mrBaseModelPart);
        }

        ModelPart& r_distance_model_part = mrModel.GetModelPart( mAuxModelPartName );

        // TODO: check flag    PERFORM_STEP1
        // Step1 - solve a poisson problem with a source term which depends on the sign of the existing distance function
        r_distance_model_part.pGetProcessInfo()->SetValue(FRACTIONAL_STEP,1);

        // Unfix the distances
        const int nnodes = static_cast<int>(r_distance_model_part.NumberOfNodes());

        block_for_each(r_distance_model_part.Nodes(), [](Node& rNode){
            double& d = rNode.FastGetSolutionStepValue(DISTANCE);

            // Free the DISTANCE values
            rNode.Free(DISTANCE);
            // Set the fix flag to 0
            rNode.Set(BLOCKED, false);

            // Save the distances
            rNode.SetValue(DISTANCE, d);

            if(d == 0){
                d = 1.0e-15;
                rNode.Set(BLOCKED, true);
                rNode.Fix(DISTANCE);
            } else {
                if(d > 0.0){
                    d = 1.0e15; // Set to a large number, to make sure that that the minimal distance is computed according to CaculateTetrahedraDistances
                } else {
                    d = -1.0e15;
                }
            }
        });

        block_for_each(r_distance_model_part.Elements(), [this](Element& rElem){
            array_1d<double,TDim+1> distances;
            auto& geom = rElem.GetGeometry();

            for(unsigned int i=0; i<TDim+1; i++){
                distances[i] = geom[i].GetValue(DISTANCE);
            }

            const array_1d<double,TDim+1> original_distances = distances;

            // The element is cut by the interface
            if(this->IsSplit(distances)){
                // Compute the unsigned distance using GeometryUtils
                if (mOptions.Is(CALCULATE_EXACT_DISTANCES_TO_PLANE)) {
                    GeometryUtils::CalculateExactDistancesToPlane(geom, distances);
                }
                else {
                    if constexpr (TDim==3){
                        GeometryUtils::CalculateTetrahedraDistances(geom, distances);
                    }
                    else {
                        GeometryUtils::CalculateTriangleDistances(geom, distances);
                    }
                }

                // Assign the sign using the original distance values
                for(unsigned int i = 0; i < TDim+1; ++i){
                    if(original_distances[i] < 0){
                        distances[i] = -distances[i];
                    }
                }

                for(unsigned int i = 0; i < TDim+1; ++i){
                    double &d = geom[i].FastGetSolutionStepValue(DISTANCE);
                    geom[i].SetLock();
                    if(std::abs(d) > std::abs(distances[i])){
                        d = distances[i];
                    }
                    geom[i].Set(BLOCKED, true);
                    geom[i].Fix(DISTANCE);
                    geom[i].UnSetLock();
                }
            }
        });

        // SHALL WE SYNCHRONIZE SOMETHING IN HERE?¿?¿??¿ WE'VE CHANGED THE NODAL DISTANCE VALUES FROM THE ELEMENTS...
        this->SynchronizeFixity();
        this->SynchronizeDistance();

        // Compute the maximum and minimum distance for the fixed nodes
        double max_dist = 0.0;
        double min_dist = 0.0;
        for(int i_node = 0; i_node < nnodes; ++i_node){
            auto it_node = r_distance_model_part.NodesBegin() + i_node;
            if(it_node->IsFixed(DISTANCE)){
                const double& d = it_node->FastGetSolutionStepValue(DISTANCE);
                if(d > max_dist){
                    max_dist = d;
                }
                if(d < min_dist){
                    min_dist = d;
                }
            }
        }

        // Synchronize the maximum and minimum distance values
        const auto &r_communicator = r_distance_model_part.GetCommunicator().GetDataCommunicator();
        max_dist = r_communicator.MaxAll(max_dist);
        min_dist = r_communicator.MinAll(min_dist);

        // Assign the max dist to all of the non-fixed positive nodes
        // and the minimum one to the non-fixed negatives
        block_for_each(r_distance_model_part.Nodes(), [&min_dist, &max_dist](Node& rNode){
            if(!rNode.IsFixed(DISTANCE)){
                double& d = rNode.FastGetSolutionStepValue(DISTANCE);
                if(d>0){
                    d = max_dist;
                } else {
                    d = min_dist;
                }
            }
        });
        mpSolvingStrategy->Solve();

        // Step2 - minimize the target residual
        r_distance_model_part.pGetProcessInfo()->SetValue(FRACTIONAL_STEP,2);
        for(unsigned int it = 0; it<mMaxIterations; it++){
             mpSolvingStrategy->Solve();
        }

        // Unfix the distances
        VariableUtils().ApplyFixity(DISTANCE, false, r_distance_model_part.Nodes());
        VariableUtils().SetFlag(BOUNDARY, false, r_distance_model_part.Nodes());
        VariableUtils().SetFlag(BLOCKED, false, r_distance_model_part.Nodes());

        KRATOS_CATCH("")
    }

    void Clear() override
    {
        if(mrModel.HasModelPart( mAuxModelPartName ))
            mrModel.DeleteModelPart( mAuxModelPartName );
        mDistancePartIsInitialized = false;

        mpSolvingStrategy->Clear();

    }

    const Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"({
            "model_part_name" : "",
            "auxiliary_model_part_name" : "RedistanceCalculationPart",
            "echo_level" : 0,
            "max_iterations" : 10,
            "calculate_exact_distances_to_plane" : false,
            "variational_redistance_coefficient_1" : 0.01,
            "variational_redistance_coefficient_2" : 0.1
        })");

        return default_parameters;
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
    std::string Info() const override
    {
        return "VariationalDistanceCalculationProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "VariationalDistanceCalculationProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
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

    bool mDistancePartIsInitialized;
    unsigned int mMaxIterations;

    Model& mrModel;
    ModelPart& mrBaseModelPart;
    Flags mOptions;
    std::string mAuxModelPartName;

    double mCoefficient1;
    double mCoefficient2;

    typename SolvingStrategyType::UniquePointer mpSolvingStrategy;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    void ValidateInput()
    {
        const DataCommunicator& r_comm = mrBaseModelPart.GetCommunicator().GetDataCommunicator();
        int num_elements = mrBaseModelPart.NumberOfElements();
        int num_nodes = mrBaseModelPart.NumberOfNodes();

        if (num_elements > 0)
        {
            const auto geometry_family = mrBaseModelPart.ElementsBegin()->GetGeometry().GetGeometryFamily();
            KRATOS_ERROR_IF( (TDim == 2) && (geometry_family != GeometryData::KratosGeometryFamily::Kratos_Triangle) )
            << "In 2D the element type is expected to be a triangle." << std::endl;
            KRATOS_ERROR_IF( (TDim == 3) && (geometry_family != GeometryData::KratosGeometryFamily::Kratos_Tetrahedra) )
            << "In 3D the element type is expected to be a tetrahedron" << std::endl;
        }

        KRATOS_ERROR_IF(r_comm.SumAll(num_nodes) == 0) << "The model part has no nodes." << std::endl;
        KRATOS_ERROR_IF(r_comm.SumAll(num_elements) == 0) << "The model Part has no elements." << std::endl;

        // Check that required nodal variables are present
        VariableUtils().CheckVariableExists<Variable<double > >(DISTANCE, mrBaseModelPart.Nodes());
    }

    void InitializeSolutionStrategy(BuilderSolverPointerType pBuilderAndSolver)
    {
        // Generate a linear strategy
        auto p_scheme = Kratos::make_shared< ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace,TDenseSpace > >();

        ModelPart& r_distance_model_part = mrModel.GetModelPart( mAuxModelPartName );

        bool CalculateReactions = false;
        bool ReformDofAtEachIteration = false;
        bool CalculateNormDxFlag = false;

        mpSolvingStrategy = Kratos::make_unique<ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver> >(
            r_distance_model_part,
            p_scheme,
            pBuilderAndSolver,
            CalculateReactions,
            ReformDofAtEachIteration,
            CalculateNormDxFlag);

        // TODO: check flag DO_EXPENSIVE_CHECKS
        mpSolvingStrategy->Check();
    }

    virtual void ReGenerateDistanceModelPart(ModelPart& rBaseModelPart)
    {
        KRATOS_TRY

        if(mrModel.HasModelPart( mAuxModelPartName ))
            mrModel.DeleteModelPart( mAuxModelPartName );

        // Ensure that the nodes have distance as a DOF
        VariableUtils().AddDof<Variable<double> >(DISTANCE, rBaseModelPart);

        // Generate
        ModelPart& r_distance_model_part = mrModel.CreateModelPart( mAuxModelPartName );

        Element::Pointer p_distance_element = Kratos::make_intrusive<DistanceCalculationElementSimplex<TDim> >();

        r_distance_model_part.GetNodalSolutionStepVariablesList() = rBaseModelPart.GetNodalSolutionStepVariablesList();

        ConnectivityPreserveModeler modeler;
        modeler.GenerateModelPart(rBaseModelPart, r_distance_model_part, *p_distance_element);

        // Using the conditions to mark the boundary with the flag boundary
        // Note that we DO NOT add the conditions to the model part
        VariableUtils().SetFlag<ModelPart::NodesContainerType>(BOUNDARY, false, r_distance_model_part.Nodes());
        // Note that above we have assigned the same geometry. Thus the flag is
        // set in the distance model part despite we are iterating the base one
        for (auto it_cond = rBaseModelPart.ConditionsBegin(); it_cond != rBaseModelPart.ConditionsEnd(); ++it_cond){
            Geometry< Node >& geom = it_cond->GetGeometry();
            for(unsigned int i=0; i<geom.size(); i++){
                geom[i].Set(BOUNDARY,true);
            }
        }

        rBaseModelPart.GetCommunicator().SynchronizeOrNodalFlags(BOUNDARY);

        r_distance_model_part.GetProcessInfo().SetValue(VARIATIONAL_REDISTANCE_COEFFICIENT_FIRST, mCoefficient1);
        r_distance_model_part.GetProcessInfo().SetValue(VARIATIONAL_REDISTANCE_COEFFICIENT_SECOND, mCoefficient2);

        mDistancePartIsInitialized = true;

        KRATOS_CATCH("")
    }

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

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    bool IsSplit(const array_1d<double,TDim+1> &rDistances){
        unsigned int positives = 0, negatives = 0;

        for(unsigned int i = 0; i < TDim+1; ++i){
            if(rDistances[i] >= 0){
                ++positives;
            } else {
                ++negatives;
            }
        }

        if (positives > 0 && negatives > 0){
            return true;
        }

        return false;
    }

    void SynchronizeDistance(){
        ModelPart& r_distance_model_part = mrModel.GetModelPart( mAuxModelPartName );
        auto &r_communicator = r_distance_model_part.GetCommunicator();

        // Only required in the MPI case
        if(r_communicator.TotalProcesses() != 1){
            int nnodes = static_cast<int>(r_distance_model_part.NumberOfNodes());

            // Set the distance absolute value
            #pragma omp parallel for
            for(int i_node = 0; i_node < nnodes; ++i_node){
                auto it_node = r_distance_model_part.NodesBegin() + i_node;
                it_node->FastGetSolutionStepValue(DISTANCE) = std::abs(it_node->FastGetSolutionStepValue(DISTANCE));
            }

            // Synchronize the unsigned value to minimum
            r_communicator.SynchronizeCurrentDataToMin(DISTANCE);

            // Set the distance sign again by retrieving it from the non-historical database
            #pragma omp parallel for
            for(int i_node = 0; i_node < nnodes; ++i_node){
                auto it_node = r_distance_model_part.NodesBegin() + i_node;
                if(it_node->GetValue(DISTANCE) < 0.0){
                    it_node->FastGetSolutionStepValue(DISTANCE) = -it_node->FastGetSolutionStepValue(DISTANCE);
                }
            }
        }
    }

    void SynchronizeFixity(){
        ModelPart& r_distance_model_part = mrModel.GetModelPart( mAuxModelPartName );
        auto &r_communicator = r_distance_model_part.GetCommunicator();

        // Only required in the MPI case
        if(r_communicator.TotalProcesses() != 1){
            int nnodes = static_cast<int>(r_distance_model_part.NumberOfNodes());

            // Synchronize the fixity flag variable to minium
            // (true means fixed and false means free)
            r_communicator.SynchronizeOrNodalFlags(BLOCKED);

            // Set the fixity according to the synchronized flag
            #pragma omp parallel for
            for(int i_node = 0; i_node < nnodes; ++i_node){
                auto it_node = r_distance_model_part.NodesBegin() + i_node;
                if (it_node->Is(BLOCKED)){
                    it_node->Fix(DISTANCE);
                }
            }
        }
    }

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
    VariationalDistanceCalculationProcess& operator=(VariationalDistanceCalculationProcess const& rOther);

    /// Copy constructor.
    //VariationalDistanceCalculationProcess(VariationalDistanceCalculationProcess const& rOther);

    ///@}
}; // Class VariationalDistanceCalculationProcess

//avoiding using the macro since this has a template parameter. If there was no template plase use the KRATOS_CREATE_LOCAL_FLAG macro
template< unsigned int TDim,class TSparseSpace, class TDenseSpace, class TLinearSolver >
const Kratos::Flags VariationalDistanceCalculationProcess<TDim,TSparseSpace,TDenseSpace,TLinearSolver>::PERFORM_STEP1(Kratos::Flags::Create(0));

template< unsigned int TDim,class TSparseSpace, class TDenseSpace, class TLinearSolver >
const Kratos::Flags VariationalDistanceCalculationProcess<TDim,TSparseSpace,TDenseSpace,TLinearSolver>::DO_EXPENSIVE_CHECKS(Kratos::Flags::Create(1));

template< unsigned int TDim,class TSparseSpace, class TDenseSpace, class TLinearSolver >
const Kratos::Flags VariationalDistanceCalculationProcess<TDim,TSparseSpace,TDenseSpace,TLinearSolver>::CALCULATE_EXACT_DISTANCES_TO_PLANE(Kratos::Flags::Create(2));

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
inline std::istream& operator >> (std::istream& rIStream,
                                  VariationalDistanceCalculationProcess<TDim,TSparseSpace,TDenseSpace,TLinearSolver>& rThis);

/// output stream function
template< unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const VariationalDistanceCalculationProcess<TDim,TSparseSpace,TDenseSpace,TLinearSolver>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

