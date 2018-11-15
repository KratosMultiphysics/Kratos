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
#include "custom_processes/postprocess_eigenvalues_process.h"
#include "structural_mechanics_application_variables.h"


namespace Kratos
{

    PostprocessEigenvaluesProcess::PostprocessEigenvaluesProcess(ModelPart& rModelPart,
                                                                 Parameters OutputParameters)
                                                                 : mrModelPart(rModelPart),
                                                                   mOutputParameters(OutputParameters)
    {
        Parameters default_parameters(R"(
            {
                "result_file_name"              : "Structure",
                "file_label"                    : "step",
                "result_file_format_use_ascii"  : false,
                "animation_steps"               : 20,
                "label_type"                    : "frequency",
                "list_of_result_variables"      : ["DISPLACEMENT"]
            }  )"
        );

        mOutputParameters.RecursivelyValidateAndAssignDefaults(default_parameters);
    }

    void PostprocessEigenvaluesProcess::ExecuteFinalizeSolutionStep()
    {
        std::string result_file_name = mOutputParameters["result_file_name"].GetString();

        if (result_file_name == "") { // use the name of the ModelPart in case nothing was assigned
            result_file_name = mrModelPart.Name();
        }

        result_file_name += "_EigenResults_";

        const std::string file_label = mOutputParameters["file_label"].GetString();

        if (file_label == "step") {
            result_file_name += std::to_string(mrModelPart.GetProcessInfo()[STEP]);
        } else if (file_label == "time") {
            result_file_name += std::to_string(mrModelPart.GetProcessInfo()[TIME]);
        } else {
            KRATOS_ERROR << "\"file_label\" can only be \"step\" or \"time\"" << std::endl;
        }

        auto post_mode = GiD_PostBinary;
        if (mOutputParameters["result_file_format_use_ascii"].GetBool()) { // this format is only needed for testing
            post_mode = GiD_PostAscii;
        }

       const auto p_gid_eigen_io = Kratos::make_unique<GidEigenIO>(
            result_file_name,
            post_mode,
            MultiFileFlag::SingleFile,
            WriteDeformedMeshFlag::WriteUndeformed,
            WriteConditionsFlag::WriteConditions);

        KRATOS_ERROR_IF_NOT(p_gid_eigen_io) << "EigenIO could not be initialized!" << std::endl;

        p_gid_eigen_io->InitializeMesh(0.0);
        p_gid_eigen_io->WriteMesh(mrModelPart.GetMesh());
        p_gid_eigen_io->WriteNodeMesh(mrModelPart.GetMesh());
        p_gid_eigen_io->FinalizeMesh();

        const auto& eigenvalue_vector = mrModelPart.GetProcessInfo()[EIGENVALUE_VECTOR];
        // Note: this is omega^2
        const SizeType num_eigenvalues = eigenvalue_vector.size();

        const SizeType num_animation_steps = mOutputParameters["animation_steps"].GetInt();
        double angle = 0.0;
        std::string label = "";

        std::vector<Variable<double>> requested_double_results;
        std::vector<Variable<array_1d<double,3>>> requested_vector_results;
        GetVariables(requested_double_results, requested_vector_results);

        p_gid_eigen_io->InitializeResults(0.0, mrModelPart.GetMesh());

        for (SizeType i=0; i < num_animation_steps; ++i)
        {
            angle = 2 * Globals::Pi * i / num_animation_steps;

            for(SizeType j=0; j<num_eigenvalues; ++j)
            {
                label = GetLabel(j, eigenvalue_vector[j]);

                for(auto& r_node : mrModelPart.Nodes())
                {
                    // Copy the eigenvector to the Solutionstepvariable. Credit to Michael Andre
                    DofsContainerType& r_node_dofs = r_node.GetDofs();
                    Matrix& r_node_eigenvectors = r_node.GetValue(EIGENVECTOR_MATRIX);

                    KRATOS_ERROR_IF_NOT(r_node_dofs.size() == r_node_eigenvectors.size2())
                        << "Number of results on node " << r_node.Id() << " is wrong" << std::endl;

                    SizeType k = 0;
                    for (auto& r_dof : r_node_dofs)
                        r_dof.GetSolutionStepValue(0) = std::cos(angle) * r_node_eigenvectors(j,k++);
                }

                for (const auto& variable : requested_double_results) {
                    p_gid_eigen_io->WriteEigenResults(mrModelPart, variable, label, i);
                }

                for (const auto& variable : requested_vector_results) {
                    p_gid_eigen_io->WriteEigenResults(mrModelPart, variable, label, i);
                }
            }
        }
        p_gid_eigen_io->FinalizeResults();
    }

    void PostprocessEigenvaluesProcess::GetVariables(std::vector<Variable<double>>& rRequestedDoubleResults,
                                                     std::vector<Variable<array_1d<double,3>>>& rRequestedVectorResults)
    {
        std::string variable_name;

        for (SizeType i=0; i<mOutputParameters["list_of_result_variables"].size(); ++i)
        {
            variable_name = mOutputParameters["list_of_result_variables"].GetArrayItem(i).GetString();

            if( KratosComponents< Variable<double> >::Has( variable_name ) ) //case of double variable
            {
                const Variable<double > variable = KratosComponents< Variable<double > >::Get(variable_name);

                KRATOS_ERROR_IF_NOT(mrModelPart.GetNodalSolutionStepVariablesList().Has( variable ))
                    << "Requesting EigenResults for a Variable that is not in the ModelPart: "
                    << variable << std::endl;

                rRequestedDoubleResults.push_back(variable);
            }
            else if (KratosComponents< Variable< array_1d<double, 3> > >::Has(variable_name) ) //case of component variable
            {
                const Variable<array_1d<double,3> > variable = KratosComponents< Variable<array_1d<double,3> > >::Get(variable_name);

                KRATOS_ERROR_IF_NOT(mrModelPart.GetNodalSolutionStepVariablesList().Has( variable ))
                    << "Requesting EigenResults for a Variable that is not in the ModelPart: "
                    << variable << std::endl;

                rRequestedVectorResults.push_back(variable);
            }
            else if (KratosComponents< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > >::Has(variable_name) ) //case of component variable
            {
                KRATOS_ERROR << "Vector Components cannot be querried, name: " << variable_name << std::endl;
            }
            else
            {
                KRATOS_ERROR << "Invalid Type of Variable, name: " << variable_name << std::endl;
            }
        }
    }

    std::string PostprocessEigenvaluesProcess::GetLabel(const int NumberOfEigenvalue,
                                                        const double EigenvalueSolution)
    {
        double label_number;

        std::stringstream parser;
        parser << (NumberOfEigenvalue + 1);
        std::string label = parser.str();

        const std::string lable_type = mOutputParameters["label_type"].GetString();

        if (lable_type == "angular_frequency")
        {
            label += "_EigenValue_[rad/s]_";
            label_number = std::sqrt(EigenvalueSolution);
        }
        else if (lable_type == "frequency")
        {
            label += "_EigenFrequency_[Hz]_";
            label_number = std::sqrt(EigenvalueSolution) / (2 * Globals::Pi);
        }
        else
        {
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
