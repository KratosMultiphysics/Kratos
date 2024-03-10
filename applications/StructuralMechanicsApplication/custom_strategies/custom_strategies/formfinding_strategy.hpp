// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Iñigo López Canalejo
//                   Klaus B. Sautter
//

#pragma once

// System includes
#include <filesystem>

// External includes

// Project includes
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
#include "custom_utilities/project_vector_on_surface_utility.h"
#include "includes/model_part_io.h"
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

    // Base class definition
    typedef ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    /// Definition of the current scheme
    typedef FormfindingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> ClassType;

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
     * @brief Default constructor
     */
    explicit FormfindingStrategy() 
    {
    }

    /**
    * Constructor.
    */

    // constructor with Builder and Solver
    FormfindingStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
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
    : BaseType(rModelPart, pScheme,
        pNewConvergenceCriteria,pNewBuilderAndSolver,MaxIterations,CalculateReactions,ReformDofSetAtEachStep,
        MoveMeshFlag),
        mProjectionSettings(ProjectionSetting),
        mpFormFindingModelPart(&rFormFindingModelPart),
        mPrintingFormat(rPrintingFormat),
        mWriteFormFoundGeometryFile(WriteFormFoundGeometryFile)
    {
        InitializeIterationIO();
    }

    /**
     * @brief Default constructor. (with parameters)
     * @param rModelPart The model part of the problem
     * @param ThisParameters The configuration parameters
     */
    explicit FormfindingStrategy(ModelPart& rModelPart, Parameters ThisParameters)
        : BaseType(rModelPart)
    {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);

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
            r_element.GetData().Clear();
            r_element.SetValue(MEMBRANE_PRESTRESS,output_matrix);
        }

        rModelPart.Conditions().clear();
        rModelPart.rProperties().clear();
        rModelPart.GetNodalSolutionStepVariablesList().clear();

        ModelPartIO model_part_io("formfinding_result_model", IO::WRITE);
        model_part_io.WriteModelPart(rModelPart);
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name"                             : "formfinding_strategy",
            "printing_format"                  : "all",
            "write_formfound_geometry_file"    : true,
            "formfinding_model_part_name"      : "",
            "projection_settings": {
                "model_part_name"             : "Structure",
                "echo_level"                  : 0,
                "projection_type"             : "planar",
                "global_direction"            : [1,0,0],
                "variable_name"               : "PLEASE_SPECIFY",
                "visualize_in_vtk"            : false,
                "method_specific_settings"    : { },
                "check_local_space_dimension" : false
            }
        })");

        // Getting base class default parameters
        const Parameters base_default_parameters = BaseType::GetDefaultParameters();
        default_parameters.RecursivelyAddMissingParameters(base_default_parameters);
        return default_parameters;
    }

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        return "formfinding_strategy";
    }
    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

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

    /**
     * @brief This method assigns settings to member variables
     * @param ThisParameters Parameters that are assigned to the member variables
     */
    void AssignSettings(const Parameters ThisParameters) override
    {
        BaseType::AssignSettings(ThisParameters);
        mpFormFindingModelPart = &BaseType::GetModelPart().GetSubModelPart(ThisParameters["formfinding_model_part_name"].GetString());
        mPrintingFormat = ThisParameters["printing_format"].GetString();
        mWriteFormFoundGeometryFile = ThisParameters["write_formfound_geometry_file"].GetBool();
        mProjectionSettings = ThisParameters["projection_settings"];
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
    void UpdateDatabase(
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b,
        const bool MoveMesh) override
    {
        BaseType::UpdateDatabase(A,Dx, b, MoveMesh);
        for(auto& r_node : mpFormFindingModelPart->Nodes()){
            // Updating reference
            const array_1d<double, 3>& disp = r_node.FastGetSolutionStepValue(DISPLACEMENT);
            array_1d<double, 3>& disp_non_historical = r_node.GetValue(DISPLACEMENT);

            disp_non_historical = disp_non_historical + disp;
            r_node.GetInitialPosition() += disp;

            r_node.FastGetSolutionStepValue(DISPLACEMENT) = ZeroVector(3);
        }
        ProjectVectorOnSurfaceUtility::Execute(*mpFormFindingModelPart, mProjectionSettings);

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
            "output_path"                        : "formfinding_results_vtk",
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
            if (std::filesystem::exists("formfinding_results_vtk")){
                std::filesystem::remove_all("formfinding_results_vtk");
            }
            std::filesystem::create_directory("formfinding_results_vtk");
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
    ModelPart* mpFormFindingModelPart;
    std::string mPrintingFormat;
    int mIterationNumber = 0;
    bool mWriteFormFoundGeometryFile = true;

}; /* Class FormfindingStrategy */
} /* namespace Kratos. */
