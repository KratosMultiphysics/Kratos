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
                "use_eigenfrequency_in_label" : false,
                "eigen_results" : []
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
        const SizeType num_eigenvalues = eigenvalue_vector.size();

        const SizeType num_animation_steps = mOutputParameters["animation_steps"].GetInt();
        double angle = 0.0;
        std::string label = "";

        const SizeType num_requested_results = mOutputParameters["eigen_results"].size();
        std::vector<std::string> requested_results(num_requested_results);

        for (SizeType i=0; i<num_requested_results; ++i)
            requested_results[i] = mOutputParameters["eigen_results"].GetArrayItem(i).GetString();

        for (SizeType i=0; i < num_animation_steps; ++i)
        {
            angle = 2 * Globals::Pi * i / num_animation_steps;

            for(SizeType j=0; j<num_eigenvalues; ++j)
            {
                label = GetLabel(eigenvalue_vector[j]);

                for(auto& node : mrModelPart.Nodes())
                {
                    // Copy the eigenvector to the Solutionstepvariable. Credit to Michael Andre
                    DofsContainerType& r_node_dofs = node.GetDofs();
                    Matrix& r_node_eigenvectors = node.GetValue(EIGENVECTOR_MATRIX);

                    KRATOS_ERROR_IF_NOT(r_node_dofs.size() == r_node_eigenvectors.size2()) 
                        << "Number of results on node " << node.Id() << " is wrong" << std::endl;

                    SizeType k = 0;
                    for (auto& r_dof : r_node_dofs)
                        r_dof.GetSolutionStepValue(0) = std::cos(angle) * r_node_eigenvectors(j,k++);
                }
                
                for (const auto& var_name : requested_results)
                    WrapperForIOCall(var_name, label, i);                
            }
        }
    }

    std::string PostprocessEigenvaluesProcess::GetLabel(double LabelNumber)
    {
        std::string label = "EigenValue_";

        // Compute the Eigenfrequency from the Eigenvalue (if the User specified it)
        if (mOutputParameters["use_eigenfrequency_in_label"].GetBool())
        {
            label = "EigenFrequency_";
            LabelNumber = std::sqrt(LabelNumber) / (2 * Globals::Pi);
        }

        std::stringstream strstr;
        strstr << std::fixed << std::setprecision(3) << LabelNumber;
        return label + strstr.str();   
    }

    void PostprocessEigenvaluesProcess::WrapperForIOCall(const std::string VariableName, 
                                                         const std::string Label,
                                                         const SizeType AnimationStepNumber)
    {
        if( KratosComponents< Variable<double> >::Has( VariableName ) ) //case of double variable
        {
            const Variable<double > r_variable = KratosComponents< Variable<double > >::Get(VariableName);
            mpGidEigenIO->WriteEigenResults(mrModelPart, 
                                            r_variable, 
                                            Label, 
                                            AnimationStepNumber); 
        }
        else if( KratosComponents< Variable< array_1d<double, 3> > >::Has(VariableName) ) //case of component variable
        {
            const Variable<array_1d<double,3> > r_variable = KratosComponents< Variable<array_1d<double,3> > >::Get(VariableName);
            mpGidEigenIO->WriteEigenResults(mrModelPart, 
                                            r_variable, 
                                            Label, 
                                            AnimationStepNumber);            
        }
        else if( KratosComponents< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > >::Has(VariableName) ) //case of component variable
        {
            KRATOS_ERROR << "Vector Components cannot be querried!" << std::endl;
        }
        else
        {
            KRATOS_ERROR << "Invalid Type of Variable" << std::endl;
        }
        

    }

}  // namespace Kratos.
