//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Altug Emiroglu, http://github.com/emiroglu
//

// System includes

// External includes

// Project includes
#include "utilities/constraint_utilities.h"
#include "custom_processes/postprocess_rom_basis_process.h"
#include "rom_application_variables.h"
#include "custom_io/gid_rom_basis_io.h"
#include "custom_io/vtk_rom_basis_output.h"

namespace Kratos
{

namespace { // helpers namespace

// helper struct to wrap the different outputs
struct BaseRomBasisOutputWrapper
{
    virtual ~BaseRomBasisOutputWrapper() = default;

    virtual void PrintOutput(
        const std::string& rLabel,
        const int AnimationStep,
        const std::vector<Variable<double>>& rRequestedDoubleResults,
        const std::vector<Variable<array_1d<double,3>>>& rRequestedVectorResults) = 0;
};

struct GidRomBasisOutputWrapper : public BaseRomBasisOutputWrapper
{
public:
    GidRomBasisOutputWrapper(ModelPart& rModelPart, Parameters OutputParameters)
        : mrModelPart(rModelPart)
    {
        KRATOS_TRY
        std::string result_file_name = OutputParameters["result_file_name"].GetString();

        if (result_file_name == "") { // use the name of the ModelPart in case nothing was assigned
            result_file_name = mrModelPart.Name();
        }

        result_file_name += "_RomBasis_";

        const std::string file_label = OutputParameters["file_label"].GetString();
        if (file_label == "step") {
            result_file_name += std::to_string(mrModelPart.GetProcessInfo()[STEP]);
        } else if (file_label == "time") {
            result_file_name += std::to_string(mrModelPart.GetProcessInfo()[TIME]);
        } else {
            KRATOS_ERROR << "\"file_label\" can only be \"step\" or \"time\"" << std::endl;
        }

        if (OutputParameters["save_output_files_in_folder"].GetBool()) {
            result_file_name = OutputParameters["folder_name"].GetString() + "/" + result_file_name;
        }

        const auto post_mode = OutputParameters["result_file_format_use_ascii"].GetBool() ? GiD_PostAscii : GiD_PostBinary;

        mpGidRomBasisIO = Kratos::make_unique<GidRomBasisIO>(
            result_file_name,
            post_mode,
            MultiFileFlag::SingleFile,
            WriteDeformedMeshFlag::WriteUndeformed,
            WriteConditionsFlag::WriteConditions);

        // deliberately rewritting the mesh in case the geometry is updated
        mpGidRomBasisIO->InitializeMesh(0.0);
        mpGidRomBasisIO->WriteMesh(mrModelPart.GetMesh());
        mpGidRomBasisIO->WriteNodeMesh(mrModelPart.GetMesh());
        mpGidRomBasisIO->FinalizeMesh();
        mpGidRomBasisIO->InitializeResults(0.0, mrModelPart.GetMesh());
        KRATOS_CATCH("")
    }

    ~GidRomBasisOutputWrapper()
    {
        mpGidRomBasisIO->FinalizeResults();
    }

    void PrintOutput(
        const std::string& rLabel,
        const int AnimationStep,
        const std::vector<Variable<double>>& rRequestedDoubleResults,
        const std::vector<Variable<array_1d<double,3>>>& rRequestedVectorResults) override
    {
        KRATOS_TRY
        for (const auto& r_variable : rRequestedDoubleResults) {
            mpGidRomBasisIO->WriteRomBasis(mrModelPart, r_variable, rLabel, AnimationStep);
        }

        for (const auto& r_variable : rRequestedVectorResults) {
            mpGidRomBasisIO->WriteRomBasis(mrModelPart, r_variable, rLabel, AnimationStep);
        }
        KRATOS_CATCH("")
    }

private:
    Kratos::unique_ptr<GidRomBasisIO> mpGidRomBasisIO;
    ModelPart& mrModelPart;

};

struct VtkRomBasisOutputWrapper : public BaseRomBasisOutputWrapper
{
public:
    VtkRomBasisOutputWrapper(ModelPart& rModelPart, Parameters OutputParameters)
    {
        KRATOS_TRY
        Parameters vtk_parameters(Parameters(R"({
            "file_format" : "binary"
        })" ));

        if (OutputParameters["result_file_format_use_ascii"].GetBool()) {
            vtk_parameters["file_format"].SetString("ascii");
        }

        mpVtkRomBasisOutput = Kratos::make_unique<VtkRomBasisOutput>(rModelPart, OutputParameters, vtk_parameters);
        KRATOS_CATCH("")
    }

    void PrintOutput(
        const std::string& rLabel,
        const int AnimationStep,
        const std::vector<Variable<double>>& rRequestedDoubleResults,
        const std::vector<Variable<array_1d<double,3>>>& rRequestedVectorResults) override
    {
        KRATOS_TRY
        mpVtkRomBasisOutput->PrintRomBasisOutput(rLabel, AnimationStep, rRequestedDoubleResults, rRequestedVectorResults);
        KRATOS_CATCH("")
    }


private:
    Kratos::unique_ptr<VtkRomBasisOutput> mpVtkRomBasisOutput;

};

} // helpers namespace

PostprocessRomBasisProcess::PostprocessRomBasisProcess(ModelPart& rModelPart,
                                                                Parameters OutputParameters)
                                                                : mrModelPart(rModelPart),
                                                                mOutputParameters(OutputParameters)
{
    KRATOS_TRY
    Parameters default_parameters(R"(
        {
            "result_file_name"              : "Structure",
            "file_format"                   : "gid",
            "file_label"                    : "step",
            "result_file_format_use_ascii"  : false,
            "folder_name"                   : "RomBasis",
            "save_output_files_in_folder"   : true,
            "animation_steps"               : 1,
            "list_of_result_variables"      : ["DISPLACEMENT"]
        }  )"
    );

    mOutputParameters.RecursivelyValidateAndAssignDefaults(default_parameters);
    KRATOS_CATCH("")
}

void PostprocessRomBasisProcess::ExecuteFinalizeSolutionStep()
{
    KRATOS_TRY
    Kratos::unique_ptr<BaseRomBasisOutputWrapper> p_rom_basis_io_wrapper;
    const std::string file_format(mOutputParameters["file_format"].GetString());
    if (file_format == "vtk") {
        p_rom_basis_io_wrapper = Kratos::make_unique<VtkRomBasisOutputWrapper>(mrModelPart, mOutputParameters);
    } else if (file_format == "gid") {
        p_rom_basis_io_wrapper = Kratos::make_unique<GidRomBasisOutputWrapper>(mrModelPart, mOutputParameters);
    } else {
        KRATOS_ERROR << "\"file_format\" can only be \"vtk\" or \"gid\"" << std::endl;
    }

    const auto& r_eigenvalue_vector = mrModelPart.GetProcessInfo()[EIGENVALUE_VECTOR];
    const auto nodes_begin = mrModelPart.NodesBegin();
    const SizeType num_rom_basis = r_eigenvalue_vector.size();
    const SizeType num_animation_steps = mOutputParameters["animation_steps"].GetInt();

    std::vector<Variable<double>> requested_double_results;
    std::vector<Variable<array_1d<double,3>>> requested_vector_results;
    GetVariables(requested_double_results, requested_vector_results);

    for (SizeType i_step=0; i_step<num_animation_steps; ++i_step) {
        const double cos_angle = std::cos(2 * Globals::Pi * i_step / num_animation_steps);

        for (SizeType i_rom_basis=0; i_rom_basis<num_rom_basis; ++i_rom_basis) {
            const std::string label = GetLabel(i_rom_basis);

            #pragma omp parallel for
            for (int i_node=0; i_node<static_cast<int>(mrModelPart.NumberOfNodes()); ++i_node) {
                // Copy the basis vector to the Solutionstepvariable. Credit to Michael Andre
                DofsContainerType& r_node_dofs = (nodes_begin+i_node)->GetDofs();
                const Matrix& r_node_rom_basis = (nodes_begin+i_node)->GetValue(ROM_BASIS);

                KRATOS_ERROR_IF_NOT(r_node_dofs.size() == r_node_rom_basis.size1())
                    << "Number of results on node #" << (nodes_begin+i_node)->Id() << " is wrong" << std::endl;

                SizeType l = 0;
                for (auto& r_dof : r_node_dofs) {
                    r_dof->GetSolutionStepValue(0) = cos_angle * r_node_rom_basis(l++,i_rom_basis);
                }

            }

            p_rom_basis_io_wrapper->PrintOutput(label, i_step, requested_double_results, requested_vector_results);
        }
    }
    KRATOS_CATCH("")
}

void PostprocessRomBasisProcess::GetVariables(std::vector<Variable<double>>& rRequestedDoubleResults,
                                                std::vector<Variable<array_1d<double,3>>>& rRequestedVectorResults) const
{
    KRATOS_TRY
    for (SizeType i=0; i<mOutputParameters["list_of_result_variables"].size(); ++i) {
        const std::string variable_name = mOutputParameters["list_of_result_variables"].GetArrayItem(i).GetString();

        if( KratosComponents< Variable<double> >::Has( variable_name ) ) {
            const Variable<double >& variable = KratosComponents< Variable<double > >::Get(variable_name);

            KRATOS_ERROR_IF_NOT(mrModelPart.GetNodalSolutionStepVariablesList().Has( variable ))
                << "Requesting RomBasis for a Variable that is not in the ModelPart: "
                << variable << std::endl;

            rRequestedDoubleResults.push_back(variable);
        } else if (KratosComponents< Variable< array_1d<double, 3> > >::Has(variable_name) ) {
            const Variable<array_1d<double,3> >& variable = KratosComponents< Variable<array_1d<double,3> > >::Get(variable_name);

            KRATOS_ERROR_IF_NOT(mrModelPart.GetNodalSolutionStepVariablesList().Has( variable ))
                << "Requesting RomBasis for a Variable that is not in the ModelPart: "
                << variable << std::endl;

            rRequestedVectorResults.push_back(variable);
        } else {
            KRATOS_ERROR << "Invalid Type of Variable, name: " << variable_name << std::endl;
        }
    }
    KRATOS_CATCH("")
}

std::string PostprocessRomBasisProcess::GetLabel(const std::size_t RomBasisIndex) const
{
    KRATOS_TRY
    std::stringstream parser;
    parser << (RomBasisIndex + 1);
    std::string label = parser.str();

    return label;
    KRATOS_CATCH("")
}

}  // namespace Kratos.
