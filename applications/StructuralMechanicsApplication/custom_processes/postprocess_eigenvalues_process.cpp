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

    PostprocessEigenvaluesProcess::PostprocessEigenvaluesProcess(ModelPart &rModelPart,
                                                                 Parameters OutputParameters) 
                                                                 : mrModelPart(rModelPart), 
                                                                   mOutputParameters(OutputParameters)
    {
        Parameters default_parameters(R"(
            {
                "result_file_name" : "",
                "animation_steps"   :  1,
                "label_type" : "angular_frequency",
                "list_of_result_variables" : []
            }  )"
        );

        mOutputParameters.RecursivelyValidateAndAssignDefaults(default_parameters);
    }


    void PostprocessEigenvaluesProcess::ExecuteInitialize()
    {
        std::string result_file_name = mOutputParameters["result_file_name"].GetString();

        if (result_file_name == "") // use the name of the ModelPart in case nothing was assigned
            result_file_name = mrModelPart.Name();
        
        result_file_name += "_EigenResults";
        
        mpGidEigenIO = GidEigenIO::Pointer (new GidEigenIO( 
                            result_file_name,
                            GiD_PostBinary,
                            MultiFileFlag::SingleFile,
                            WriteDeformedMeshFlag::WriteUndeformed,
                            WriteConditionsFlag::WriteConditions) );

        KRATOS_ERROR_IF_NOT(mpGidEigenIO) << "EigenIO could not be initialized!" << std::endl;
    }

    void PostprocessEigenvaluesProcess::ExecuteBeforeSolutionLoop()
    {
        KRATOS_ERROR_IF_NOT(mpGidEigenIO) << " EigenIO is uninitialized!" << std::endl;

        mpGidEigenIO->InitializeMesh(0.0);
        mpGidEigenIO->WriteMesh(mrModelPart.GetMesh());
        mpGidEigenIO->WriteNodeMesh(mrModelPart.GetMesh());
        mpGidEigenIO->FinalizeMesh();
    }

    void PostprocessEigenvaluesProcess::ExecuteFinalize()
    {
        KRATOS_ERROR_IF_NOT(mpGidEigenIO) << " EigenIO is uninitialized!" << std::endl;

        const auto& eigenvalue_vector = mrModelPart.GetProcessInfo()[EIGENVALUE_VECTOR];
        // Note: this is omega^2
        const SizeType num_eigenvalues = eigenvalue_vector.size();

        const SizeType num_animation_steps = mOutputParameters["animation_steps"].GetInt();
        double angle = 0.0;
        std::string label = "";

        std::vector<Variable<double>> requested_double_results;
        std::vector<Variable<array_1d<double,3>>> requested_vector_results;
        GetVariables(requested_double_results, requested_vector_results);

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
                
                for (const auto& variable : requested_double_results)
                    mpGidEigenIO->WriteEigenResults(mrModelPart, variable, label, i);                
            
                for (const auto& variable : requested_vector_results)
                    mpGidEigenIO->WriteEigenResults(mrModelPart, variable, label, i);               
            }
        }
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
            else if( KratosComponents< Variable< array_1d<double, 3> > >::Has(variable_name) ) //case of component variable
            {
                const Variable<array_1d<double,3> > variable = KratosComponents< Variable<array_1d<double,3> > >::Get(variable_name);

                KRATOS_ERROR_IF_NOT(mrModelPart.GetNodalSolutionStepVariablesList().Has( variable )) 
                    << "Requesting EigenResults for a Variable that is not in the ModelPart: " 
                    << variable << std::endl;

                rRequestedVectorResults.push_back(variable);       
            }
            else if( KratosComponents< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > >::Has(variable_name) ) //case of component variable
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
        std::string label;
        double label_number;

        const std::string lable_type = mOutputParameters["label_type"].GetString();

        if (lable_type == "angular_frequency")
        {
            label = "EigenValue_[rad/s]_";
            label_number = std::sqrt(EigenvalueSolution);
        }
        else if (lable_type == "frequency")
        {
            label = "EigenFrequency_[Hz]_";
            label_number = std::sqrt(EigenvalueSolution) / (2 * Globals::Pi);
        }
        else if (lable_type == "step")
        {
            label = "NoOfEigenValue_";
            label_number = NumberOfEigenvalue;
        }
        else
        {
            KRATOS_ERROR << "Wrong \"lable_type\"! Available options are: \"angular_frequency\","
                         << "\"frequency\", \"step\"" << std::endl;
        }

        std::stringstream strstr;
        strstr << label_number;
        return label + strstr.str();   
    }

}  // namespace Kratos.
