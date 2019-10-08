//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//

// System includes

// External includes

// Project includes
#include "utilities/constraint_utilities.h"
#include "custom_processes/postprocess_eigenvalues_process.h"
#include "structural_mechanics_application_variables.h"
#include "custom_io/gid_eigen_io.h"
#include "custom_io/vtk_eigen_output.h"

namespace Kratos
{

namespace { // helpers namespace

// helper struct to wrap the different outputs
struct BaseEigenOutputWrapper
{
    virtual ~BaseEigenOutputWrapper() = default;

    virtual void PrintOutput(
        const std::string& rLabel,
        const int AnimationStep,
        const std::vector<Variable<double>>& rRequestedDoubleResults,
        const std::vector<Variable<array_1d<double,3>>>& rRequestedVectorResults) = 0;
};

struct GidEigenOutputWrapper : public BaseEigenOutputWrapper
{
public:
    GidEigenOutputWrapper(ModelPart& rModelPart, Parameters OutputParameters)
        : mrModelPart(rModelPart)
    {
        std::string result_file_name = OutputParameters["result_file_name"].GetString();

        if (result_file_name == "") { // use the name of the ModelPart in case nothing was assigned
            result_file_name = mrModelPart.Name();
        }

        result_file_name += "_EigenResults_";

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

        mpGidEigenIO = Kratos::make_unique<GidEigenIO>(
            result_file_name,
            post_mode,
            MultiFileFlag::SingleFile,
            WriteDeformedMeshFlag::WriteUndeformed,
            WriteConditionsFlag::WriteConditions);

        // deliberately rewritting the mesh in case the geometry is updated
        mpGidEigenIO->InitializeMesh(0.0);
        mpGidEigenIO->WriteMesh(mrModelPart.GetMesh());
        mpGidEigenIO->WriteNodeMesh(mrModelPart.GetMesh());
        mpGidEigenIO->FinalizeMesh();
        mpGidEigenIO->InitializeResults(0.0, mrModelPart.GetMesh());
    }

    ~GidEigenOutputWrapper()
    {
        mpGidEigenIO->FinalizeResults();
    }

    void PrintOutput(
        const std::string& rLabel,
        const int AnimationStep,
        const std::vector<Variable<double>>& rRequestedDoubleResults,
        const std::vector<Variable<array_1d<double,3>>>& rRequestedVectorResults) override
    {
        for (const auto& r_variable : rRequestedDoubleResults) {
            mpGidEigenIO->WriteEigenResults(mrModelPart, r_variable, rLabel, AnimationStep);
        }

        for (const auto& r_variable : rRequestedVectorResults) {
            mpGidEigenIO->WriteEigenResults(mrModelPart, r_variable, rLabel, AnimationStep);
        }
    }

private:
    Kratos::unique_ptr<GidEigenIO> mpGidEigenIO;
    ModelPart& mrModelPart;

};

struct VtkEigenOutputWrapper : public BaseEigenOutputWrapper
{
public:
    VtkEigenOutputWrapper(ModelPart& rModelPart, Parameters OutputParameters)
    {
        Parameters vtk_parameters(Parameters(R"({
            "file_format" : "binary"
        })" ));

        if (OutputParameters["result_file_format_use_ascii"].GetBool()) {
            vtk_parameters["file_format"].SetString("ascii");
        }

        mpVtkEigenOutput = Kratos::make_unique<VtkEigenOutput>(rModelPart, OutputParameters, vtk_parameters);
    }

    void PrintOutput(
        const std::string& rLabel,
        const int AnimationStep,
        const std::vector<Variable<double>>& rRequestedDoubleResults,
        const std::vector<Variable<array_1d<double,3>>>& rRequestedVectorResults) override
    {
        mpVtkEigenOutput->PrintEigenOutput(rLabel, AnimationStep, rRequestedDoubleResults, rRequestedVectorResults);
    }


private:
    Kratos::unique_ptr<VtkEigenOutput> mpVtkEigenOutput;

};

} // helpers namespace

PostprocessEigenvaluesProcess::PostprocessEigenvaluesProcess(ModelPart& rModelPart,
                                                                Parameters OutputParameters)
                                                                : mrModelPart(rModelPart),
                                                                mOutputParameters(OutputParameters)
{
    Parameters default_parameters(R"(
        {
            "result_file_name"              : "Structure",
            "file_format"                   : "vtk",
            "file_label"                    : "step",
            "result_file_format_use_ascii"  : false,
            "folder_name"                   : "EigenResults",
            "save_output_files_in_folder"   : true,
            "animation_steps"               : 20,
            "label_type"                    : "frequency",
            "list_of_result_variables"      : ["DISPLACEMENT"]
        }  )"
    );

    mOutputParameters.RecursivelyValidateAndAssignDefaults(default_parameters);
}

void PostprocessEigenvaluesProcess::ExecuteFinalizeSolutionStep()
{
    Kratos::unique_ptr<BaseEigenOutputWrapper> p_eigen_io_wrapper;
    const std::string file_format(mOutputParameters["file_format"].GetString());
    if (file_format == "vtk") {
        p_eigen_io_wrapper = Kratos::make_unique<VtkEigenOutputWrapper>(mrModelPart, mOutputParameters);
    } else if (file_format == "gid") {
        p_eigen_io_wrapper = Kratos::make_unique<GidEigenOutputWrapper>(mrModelPart, mOutputParameters);
    } else {
        KRATOS_ERROR << "\"file_format\" can only be \"vtk\" or \"gid\"" << std::endl;
    }

    const auto& eigenvalue_vector = mrModelPart.GetProcessInfo()[EIGENVALUE_VECTOR];
    // Note: this is omega^2
    const SizeType num_eigenvalues = eigenvalue_vector.size();
    const auto nodes_begin = mrModelPart.NodesBegin();
    const SizeType num_animation_steps = mOutputParameters["animation_steps"].GetInt();

    std::vector<Variable<double>> requested_double_results;
    std::vector<Variable<array_1d<double,3>>> requested_vector_results;
    GetVariables(requested_double_results, requested_vector_results);


    for (SizeType i=0; i<num_animation_steps; ++i) {
        const double cos_angle = std::cos(2 * Globals::Pi * i / num_animation_steps);

        for (SizeType j=0; j<num_eigenvalues; ++j) {
            const std::string label = GetLabel(j, eigenvalue_vector[j]);

            #pragma omp parallel for
            for (int k=0; k<static_cast<int>(mrModelPart.NumberOfNodes()); ++k) {
                // Copy the eigenvector to the Solutionstepvariable. Credit to Michael Andre
                DofsContainerType& r_node_dofs = (nodes_begin+k)->GetDofs();
                const Matrix& r_node_eigenvectors = (nodes_begin+k)->GetValue(EIGENVECTOR_MATRIX);

                KRATOS_ERROR_IF_NOT(r_node_dofs.size() == r_node_eigenvectors.size2())
                    << "Number of results on node #" << (nodes_begin+i)->Id() << " is wrong" << std::endl;

                SizeType l = 0;
                for (auto& r_dof : r_node_dofs) {
                    r_dof->GetSolutionStepValue(0) = cos_angle * r_node_eigenvectors(j,l++);
                }
            }

            // Reconstruct the animation on slave-dofs
            if (mrModelPart.NumberOfMasterSlaveConstraints() > 0) {
                ConstraintUtilities::ResetSlaveDofs(mrModelPart);
                ConstraintUtilities::ApplyConstraints(mrModelPart);
            }

            p_eigen_io_wrapper->PrintOutput(label, i, requested_double_results, requested_vector_results);
        }
    }
}

void PostprocessEigenvaluesProcess::GetVariables(std::vector<Variable<double>>& rRequestedDoubleResults,
                                                std::vector<Variable<array_1d<double,3>>>& rRequestedVectorResults) const
{
    for (SizeType i=0; i<mOutputParameters["list_of_result_variables"].size(); ++i) {
        const std::string variable_name = mOutputParameters["list_of_result_variables"].GetArrayItem(i).GetString();

        if( KratosComponents< Variable<double> >::Has( variable_name ) ) {
            const Variable<double > variable = KratosComponents< Variable<double > >::Get(variable_name);

            KRATOS_ERROR_IF_NOT(mrModelPart.GetNodalSolutionStepVariablesList().Has( variable ))
                << "Requesting EigenResults for a Variable that is not in the ModelPart: "
                << variable << std::endl;

            rRequestedDoubleResults.push_back(variable);
        } else if (KratosComponents< Variable< array_1d<double, 3> > >::Has(variable_name) ) {
            const Variable<array_1d<double,3> > variable = KratosComponents< Variable<array_1d<double,3> > >::Get(variable_name);

            KRATOS_ERROR_IF_NOT(mrModelPart.GetNodalSolutionStepVariablesList().Has( variable ))
                << "Requesting EigenResults for a Variable that is not in the ModelPart: "
                << variable << std::endl;

            rRequestedVectorResults.push_back(variable);
        } else if (KratosComponents< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > >::Has(variable_name) ) {
            KRATOS_ERROR << "Vector Components cannot be querried, name: " << variable_name << std::endl;
        } else {
            KRATOS_ERROR << "Invalid Type of Variable, name: " << variable_name << std::endl;
        }
    }
}

std::string PostprocessEigenvaluesProcess::GetLabel(const int NumberOfEigenvalue,
                                                    const double EigenvalueSolution) const
{
    double label_number;

    std::stringstream parser;
    parser << (NumberOfEigenvalue + 1);
    std::string label = parser.str();

    const std::string lable_type = mOutputParameters["label_type"].GetString();

    if (lable_type == "angular_frequency") {
        label += "_EigenValue_[rad/s]_";
        label_number = std::sqrt(EigenvalueSolution);
    } else if (lable_type == "frequency") {
        label += "_EigenFrequency_[Hz]_";
        label_number = std::sqrt(EigenvalueSolution) / (2 * Globals::Pi);
    } else {
        KRATOS_ERROR << "The requested label_type \"" << lable_type << "\" is not available!\n"
                        << "Available options are: \"angular_frequency\", \"frequency\"" << std::endl;
    }

    // reset the stringstream
    parser.str( std::string() );
    parser.clear();

    parser << label_number;
    return label + parser.str();
}

}  // namespace Kratos.
