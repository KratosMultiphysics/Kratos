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
#include "includes/gid_io.h"
#include "structural_mechanics_application_variables.h"

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
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    typedef GidIO<> IterationIOType;
    //typedef IterationIOType::Pointer IterationIOPointerType;
    typedef Kratos::unique_ptr<IterationIOType> IterationIOPointerType;

    ///@}
    ///@name Life Cycle

    ///@{

    /**
    * Constructor.
    */

    // constructor with Builder and Solver
    FormfindingStrategy(
        ModelPart& model_part,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        ModelPart& rFormFindingModelPart,
        bool WriteFormFoundGeometryFile,
        const std::string& rPrintingFormat,
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
        mrFormFindingModelPart(rFormFindingModelPart),
        mPrintingFormat(rPrintingFormat),
        mWriteFormFoundGeometryFile(WriteFormFoundGeometryFile)
    {
        InitializeIterationIO();
    }

    // Destructor
    ~FormfindingStrategy() = default;

    static void WriteFormFoundMdpa(ModelPart& rModelPart)
    {
        Matrix output_matrix;
        const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
        for(auto& r_element : rModelPart.Elements()){
            r_element.Calculate(MEMBRANE_PRESTRESS,output_matrix,r_process_info);
            r_element.Data().Clear();
            r_element.SetValue(MEMBRANE_PRESTRESS,output_matrix);
        }

        rModelPart.Conditions().clear();
        rModelPart.rProperties().clear();
        rModelPart.GetNodalSolutionStepVariablesList().clear();

        ModelPartIO model_part_io("formfinding_result_model", IO::WRITE);
        model_part_io.WriteModelPart(rModelPart);
    }

private:
    void UpdateDatabase(
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b,
        const bool MoveMesh) override
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

        PrintResults();
    }

    void EchoInfo(const unsigned int IterationNumber) override
    {
        BaseType::EchoInfo(IterationNumber);
        mIterationNumber = IterationNumber;
    }

    void PrintResults()
    {
        if (mPrintingFormat=="all"){
            PrintVtkFiles(mIterationNumber);
            PrintGiDFiles(mIterationNumber);
        }
        else if (mPrintingFormat=="vtk") PrintVtkFiles(mIterationNumber);
        else if (mPrintingFormat=="gid") PrintGiDFiles(mIterationNumber);
        else;
    }

    void FinalizeSolutionStep() override
    {
        BaseType::FinalizeSolutionStep();
        if (mPrintingFormat=="all" || mPrintingFormat=="gid") mpIterationIO->FinalizeResults();
    }



    void PrintVtkFiles(const int rIterationNumber)
    {
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
        postfix << std::setw(max_prefix) << std::setfill('0') << rIterationNumber;
        VtkOutput(BaseType::GetModelPart(), vtk_params).PrintOutput("formfinding_"+postfix.str());
    }

    void PrintGiDFiles(const int rIterationNumber)
    {
        double solution_tag = rIterationNumber;
        mpIterationIO->WriteNodalResultsNonHistorical(DISPLACEMENT,BaseType::GetModelPart().Nodes(),solution_tag);
    }

    void InitializeIterationIO()
    {
        // check user input for visualization of results
        std::vector<std::string> printing_possibilities {"all","vtk","gid","none"};
        if (std::find(printing_possibilities.begin(), printing_possibilities.end(), mPrintingFormat)==printing_possibilities.end()){
            KRATOS_ERROR << "Chosen printing format :" << mPrintingFormat << " is not available. Please use: " << printing_possibilities << std::endl;
        }

        // initialize i/o
        if (mPrintingFormat=="all" || mPrintingFormat=="vtk"){
            if (Kratos::filesystem::exists("formfinding_results_vtk")){
                Kratos::filesystem::remove_all("formfinding_results_vtk");
            }
            Kratos::filesystem::create_directory("formfinding_results_vtk");
            PrintVtkFiles(mIterationNumber);
        }

        if (mPrintingFormat=="all" || mPrintingFormat=="gid"){
            mpIterationIO = Kratos::make_unique<IterationIOType>(
                "formfinding_iterations",
                GiD_PostBinary,
                MultiFileFlag::SingleFile,
                WriteDeformedMeshFlag::WriteUndeformed,
                WriteConditionsFlag::WriteConditions);

            mpIterationIO->InitializeMesh(0.0);
            mpIterationIO->WriteMesh(BaseType::GetModelPart().GetMesh());
            mpIterationIO->WriteNodeMesh(BaseType::GetModelPart().GetMesh());
            mpIterationIO->FinalizeMesh();
            mpIterationIO->InitializeResults(0.0, BaseType::GetModelPart().GetMesh());
        }
    }


    IterationIOPointerType mpIterationIO;
    Parameters mProjectionSettings;
    ModelPart& mrFormFindingModelPart;
    std::string mPrintingFormat;
    int mIterationNumber = 0;
    bool mWriteFormFoundGeometryFile = true;

}; /* Class FormfindingStrategy */
} /* namespace Kratos. */

#endif /* KRATOS_FORMFINDING_UPDATED_REFERENCE_STRATEGY defined */
