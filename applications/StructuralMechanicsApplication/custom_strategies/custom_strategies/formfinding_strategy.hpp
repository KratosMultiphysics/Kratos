// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Iñigo López Canalejo
//                   Klaus B. Sautter
//

#if !defined(KRATOS_FORMFINDING_STRATEGY )
#define  KRATOS_FORMFINDING_STRATEGY
// System includes

// External includes

// Project includes
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
#include "custom_utilities/project_vector_on_surface_utility.h"
#include "includes/model_part_io.h"
#include "includes/kratos_filesystem.h"
#include "input_output/vtk_output.h"

namespace Kratos
{
    template<class TSparseSpace,
    class TDenseSpace, // = DenseSpace<double>,
    class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
    >
    class FormfindingStrategy
        : public ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
    {
        public:
            ///@name Type Definitions
            ///@{
            typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> TConvergenceCriteriaType;

            // Counted pointer of ClassName
            KRATOS_CLASS_POINTER_DEFINITION(FormfindingStrategy);

            typedef ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
            typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;
            typedef typename BaseType::TSchemeType TSchemeType;
            typedef GidIO<> IterationIOType;
            typedef IterationIOType::Pointer IterationIOPointerType;
            typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
            typedef typename BaseType::TSystemVectorType TSystemVectorType;

            ///@}
            ///@name Life Cycle

            ///@{

            /**
            * Constructor.
            */

            FormfindingStrategy(
                ModelPart& model_part,
                typename TSchemeType::Pointer pScheme,
                typename TLinearSolver::Pointer pNewLinearSolver,
                typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
                ModelPart& rFormFindingModelPart,
                Parameters ProjectionSetting,
                int MaxIterations = 30,
                bool CalculateReactions = false,
                bool ReformDofSetAtEachStep = false,
                bool MoveMeshFlag = false
                )
                : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, pScheme,
                    pNewLinearSolver,
                    pNewConvergenceCriteria,
                    MaxIterations,
                    CalculateReactions,
                    ReformDofSetAtEachStep,
                    MoveMeshFlag),
                    mProjectionSettings(ProjectionSetting),
                    mrFormFindingModelPart(rFormFindingModelPart)
                 {
                    InitializeIterationIO();
                 }

            // constructor with Builder and Solver
            FormfindingStrategy(
                ModelPart& model_part,
                typename TSchemeType::Pointer pScheme,
                typename TLinearSolver::Pointer pNewLinearSolver,
                typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
                typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
                ModelPart& rFormFindingModelPart,
                Parameters ProjectionSetting,
                int MaxIterations = 30,
                bool CalculateReactions = false,
                bool ReformDofSetAtEachStep = false,
                bool MoveMeshFlag = false
                )
                : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, pScheme,
                    pNewLinearSolver,pNewConvergenceCriteria,pNewBuilderAndSolver,MaxIterations,CalculateReactions,ReformDofSetAtEachStep,
                    MoveMeshFlag),
                    mProjectionSettings(ProjectionSetting),
                    mrFormFindingModelPart(rFormFindingModelPart)
                 {
                    InitializeIterationIO();
                 }

            // Destructor
            ~FormfindingStrategy() = default;

        private:
            bool SolveSolutionStep() override
            {

                /* KRATOS_ERROR_IF_NOT(mpIterationIO) << " IterationIO is uninitialized!" << std::endl;
                mpIterationIO->InitializeResults(0.0, BaseType::GetModelPart().GetMesh()); */

                BaseType::SolveSolutionStep();

                /* mpIterationIO->FinalizeResults(); */

                return true;
            }

            void UpdateDatabase(
                    TSystemMatrixType& A,
                    TSystemVectorType& Dx,
                    TSystemVectorType& b,
                    const bool MoveMesh
                ) override
                {
                    BaseType::UpdateDatabase(A,Dx, b, MoveMesh);
                    for(auto& r_node : mrFormFindingModelPart.Nodes()){
                        // Updating reference
                        const array_1d<double, 3>& disp = r_node.FastGetSolutionStepValue(DISPLACEMENT);
                        array_1d<double, 3>& disp_non_historical = r_node.GetValue(DISPLACEMENT);

                        disp_non_historical = disp_non_historical + disp;
                        r_node.GetInitialPosition() += disp;

                        r_node.FastGetSolutionStepValue(DISPLACEMENT) = ZeroVector(3);
                    }
                    ProjectVectorOnSurfaceUtility::Execute(mrFormFindingModelPart, mProjectionSettings);
                }

                FormfindingStrategy(const FormfindingStrategy& Other) {};

                void EchoInfo(const unsigned int IterationNumber) override
                {
                    BaseType::EchoInfo(IterationNumber);

                    /* KRATOS_ERROR_IF_NOT(mpIterationIO) << " IterationIO is uninitialized!" << std::endl;
                    mpIterationIO->WriteNodalResults(DISPLACEMENT, BaseType::GetModelPart().Nodes(), IterationNumber, 0); */

                    Parameters vtk_params( R"({
                        "file_format"                        : "binary",
                        "output_precision"                   : 7,
                        "output_control_type"                : "step",
                        "save_output_files_in_folder"        : true,
                        "folder_name"                        : "formfinding_results_vtk",
                        "nodal_data_value_variables"         : ["DISPLACEMENT"]
                    })");

                    const int max_prefix = int(std::floor(std::log10(BaseType::mMaxIterationNumber)))+1;
                    std::stringstream postfix;
                    postfix << std::setw(max_prefix) << std::setfill('0') << IterationNumber;

                    VtkOutput(BaseType::GetModelPart(), vtk_params).PrintOutput("formfinding_"+postfix.str());
                }
                

                void InitializeIterationIO()
                {
                    if (Kratos::filesystem::exists("formfinding_results_vtk")){
                        Kratos::filesystem::remove_all("formfinding_results_vtk"); 
                    }       
                    Kratos::filesystem::create_directory("formfinding_results_vtk");
                }


                void FinalizeSolutionStep() override
                {
                    BaseType::FinalizeSolutionStep();

                    // write new mdpa
                    // erase properties
                    for (auto& prop: BaseType::GetModelPart().rProperties())
                        prop.Data().Clear();

                    ModelPartIO model_part_io("formfinding_result_model", IO::WRITE);
                    model_part_io.WriteModelPart(BaseType::GetModelPart());


                    
                }


                /* IterationIOPointerType mpIterationIO; */
                Parameters mProjectionSettings;
                ModelPart& mrFormFindingModelPart;
        }; /* Class FormfindingStrategy */
} /* namespace Kratos. */

#endif /* KRATOS_FORMFINDING_UPDATED_REFERENCE_STRATEGY defined */
