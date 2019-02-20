// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder
//


// System includes

// External includes

// Project includes
#include "utilities/openmp_utils.h"
#include "utilities/variable_utils.h"
#include "direct_sensitivity_postprocess.h"
#include "includes/kratos_parameters.h"
#include "utilities/math_utils.h"
//#include <iomanip>

namespace Kratos
{

    /// Constructor.
    DirectSensitivityPostprocess::DirectSensitivityPostprocess(ModelPart& rModelPart, DirectSensitivityVariable& rDirectSensitivityVariable,
                                                Parameters SensitivitySettings)
      : mrModelPart(rModelPart) , mrDirectSensitivityVariable(rDirectSensitivityVariable)
    {
        KRATOS_TRY;

        Parameters default_settings(R"(
        {
            "sensitivity_model_part_name": "PLEASE_SPECIFY_SENSITIVITY_MODEL_PART",
            "build_mode": "static"
        })");

        SensitivitySettings.ValidateAndAssignDefaults(default_settings);

        auto sensitivity_model_part_name =
            SensitivitySettings["sensitivity_model_part_name"].GetString();
        if (sensitivity_model_part_name !=
            "PLEASE_SPECIFY_SENSITIVITY_MODEL_PART")
        {
            mpSensitivityModelPart = &mrModelPart.GetSubModelPart(sensitivity_model_part_name);
        }
        else
        {
            mpSensitivityModelPart = &mrModelPart;
        }

        mBuildMode = SensitivitySettings["build_mode"].GetString(); 

        //mVariableType = std::string("element_data_type");

        mVariableType = mrDirectSensitivityVariable.GetDesignVariableType();  

        //design_variable_name = std::string("I22");

        mDesignVariableName = mrDirectSensitivityVariable.GetDesignVariableName();       
        
        KRATOS_CATCH("");
    }



    /// Destructor.
    DirectSensitivityPostprocess::~DirectSensitivityPostprocess()
    {
    }



    void DirectSensitivityPostprocess::Initialize()
    {
        /*KRATOS_TRY;

        this->Clear();

        // Initialize flags.
        VariableUtils().SetNonHistoricalVariable(UPDATE_SENSITIVITIES, true, mpSensitivityModelPart->Nodes());
        VariableUtils().SetNonHistoricalVariable(UPDATE_SENSITIVITIES, true, mpSensitivityModelPart->Elements());
        VariableUtils().SetNonHistoricalVariable(UPDATE_SENSITIVITIES, true, mpSensitivityModelPart->Conditions());

        KRATOS_CATCH("");*/
    }



    void DirectSensitivityPostprocess::Clear()
    {
        /*KRATOS_TRY;    

        // Reset flags.
        VariableUtils().SetNonHistoricalVariable(UPDATE_SENSITIVITIES, false, mrModelPart.Nodes());
        VariableUtils().SetNonHistoricalVariable(UPDATE_SENSITIVITIES, false, mrModelPart.Elements());
        VariableUtils().SetNonHistoricalVariable(UPDATE_SENSITIVITIES, false, mrModelPart.Conditions());

        this->SetAllSensitivityVariablesToZero();         

        KRATOS_CATCH("");*/
    }



    void DirectSensitivityPostprocess::SetAllSensitivityVariablesToZero()
    {
        KRATOS_TRY;
                
        if (mVariableType == "nodal_data_type")
        {
            if (KratosComponents<Variable<double>>::Has(mDesignVariableName))                                          
                VariableUtils().SetToZero_ScalarVar( ReadScalarSensitivityVariables(mDesignVariableName), mrModelPart.Nodes() );
            else if (KratosComponents<Variable<array_1d<double,3>>>::Has(mDesignVariableName))
                VariableUtils().SetToZero_VectorVar( ReadVectorSensitivityVariables(mDesignVariableName), mrModelPart.Nodes() ); 
            else
                KRATOS_ERROR << "Unsupported variable: " <<  mDesignVariableName << "." << std::endl;
        }              
          
        // Set element sensitivity result variables to zero.
        if (mVariableType == "element_data_type")
        {
            #pragma omp parallel for
            for (int i = 0; i< static_cast<int> (mrModelPart.NumberOfElements()); ++i)
            {
                auto it = mrModelPart.ElementsBegin() + i;
                
                if ( KratosComponents<Variable<double>>::Has(mDesignVariableName) )
                    it->SetValue(ReadScalarSensitivityVariables(mDesignVariableName), ReadScalarSensitivityVariables(mDesignVariableName).Zero());
                else if ( KratosComponents<Variable<array_1d<double,3>>>::Has(mDesignVariableName) )
                    it->SetValue(ReadVectorSensitivityVariables(mDesignVariableName), ReadVectorSensitivityVariables(mDesignVariableName).Zero());  
                else
                    KRATOS_ERROR << "Unsupported variable: " <<  mDesignVariableName << "." << std::endl;              
            }
        }
        
        // Set conditional sensitivity result variables to zero.
        if (mVariableType == "condition_data_type")
        {
            #pragma omp parallel for
            for (int i = 0; i< static_cast<int> (mrModelPart.NumberOfConditions()); ++i)
            {
                auto it = mrModelPart.ConditionsBegin() + i;
                
                if ( KratosComponents<Variable<double>>::Has(mDesignVariableName) )
                    it->SetValue(ReadScalarSensitivityVariables(mDesignVariableName), ReadScalarSensitivityVariables(mDesignVariableName).Zero());
                else if ( KratosComponents<Variable<array_1d<double,3>>>::Has(mDesignVariableName) )
                    it->SetValue(ReadVectorSensitivityVariables(mDesignVariableName), ReadVectorSensitivityVariables(mDesignVariableName).Zero());
                else
                    KRATOS_ERROR << "Unsupported variable: " <<  mDesignVariableName << "." << std::endl;
            }
        }

        KRATOS_CATCH("");
    }    



    void DirectSensitivityPostprocess::UpdateSensitivities(DirectSensitivityResponseFunction& Resp)
    {
        KRATOS_TRY;

        std::cout << "in UpdateSensitivities()" << std::endl;

        if (mBuildMode == "static")
        {
            // overwrite existing.
            this->SetAllSensitivityVariablesToZero();
        }
        else
        {
            KRATOS_ERROR << "Unsupported \"build_mode\": " << mBuildMode << std::endl;
        }
        
        // Define response variables ("MOMENT", "DISPLACEMENT" etc.) for which the sensitivity should get computed
        std::string response_variable_name = Resp.GetResponseVariableName();
        std::cout << response_variable_name << std::endl; 

        // Define response output variables ("MOMENT_SENSITIVITY", "DISPLACEMENT_SENSITIVITY" etc.) to save the results 
        std::string response_output_variable_name = response_variable_name + std::string("_SENSITIVITY");
                
        //  Update the sensitivies        
        if ( KratosComponents<Variable<array_1d<double,3>>>::Has(response_variable_name) )
        {
            if ( KratosComponents<Variable<array_1d<double,3>>>::Has(response_output_variable_name) )
            { 
                if (Resp.GetEvaluationFlag() == "on_gauss_point")
                    this->UpdateSensitivityOnGaussPoint( Resp, ReadVectorSensitivityVariables(response_variable_name), 
                                                ReadVectorSensitivityVariables(response_output_variable_name) );
                if (Resp.GetEvaluationFlag() == "on_node")
                    this->UpdateSensitivityOnNode( Resp, ReadVectorSensitivityVariables(std::string("ADJOINT_") + response_variable_name), 
                                                ReadVectorSensitivityVariables(response_output_variable_name) );
            }
            else
                KRATOS_ERROR << "No matching output variable for " << response_variable_name << " exist." << std::endl;
        }
        else
            KRATOS_ERROR << "Unsupported variable: " <<  response_variable_name << "." << std::endl;
    
        KRATOS_CATCH("");
    }



    template <typename TDataType>
    void DirectSensitivityPostprocess::UpdateSensitivityOnGaussPoint(DirectSensitivityResponseFunction& Resp, Variable<TDataType> const& rResponseVariable, 
                        Variable<TDataType> const& rOutputVariable)
    {
        KRATOS_TRY;
        
        ProcessInfo& r_process_info = mrModelPart.GetProcessInfo();
        const int num_threads = 1;
        std::vector<std::vector<TDataType>> sensitivity_vector(num_threads);
        std::vector<std::vector<std::vector<TDataType>>> response_displacement_gradient(num_threads);
        std::vector<std::vector<TDataType>> response_sensitivity_gradient(num_threads);
        std::vector<Vector> displacement_gradient(num_threads);
        std::vector<std::vector<TDataType>> scalar_product(num_threads);
        
        int k = 0;
        
        for (auto& elem_i : mrModelPart.Elements())
        {   
            std::cout << "UpdateElement: ELEMENT: "<< elem_i.Id() << std::endl;
            //Element::GeometryType& r_geom = elem_i.GetGeometry();
            
            // No element contribution if the design variable is an condition data type
            if(mVariableType == "condition_data_type")
                break;  

            // Calculate derivative of the response function wrt. the displacements             
            Resp.CalculateGradient(elem_i, rResponseVariable, response_displacement_gradient[k], r_process_info);

            //OutputUtility::OutputOnTerminal("dg/du", response_displacement_gradient[k]);

            // Calculate derivative of the response function wrt. the design variable
            Resp.CalculatePartialSensitivity(elem_i, mrDirectSensitivityVariable, rResponseVariable,
                                                response_sensitivity_gradient[k], r_process_info);            

            // Size sensitivity vector
            const SizeType num_traced_integr_pts = response_displacement_gradient[k][0].size();
            sensitivity_vector[k].resize(num_traced_integr_pts);
            scalar_product[k].resize(num_traced_integr_pts);
                                           
            // Get the displacement vector derived wrt. the design parameter
            elem_i.GetValuesVector(displacement_gradient[k]);
                        
            //OutputUtility::OutputOnTerminal("du/ds", displacement_gradient[k]);
                        
            if( (response_displacement_gradient[k].size() > 0) && (displacement_gradient[k].size() > 0) )
            {
                KRATOS_ERROR_IF_NOT( response_displacement_gradient[k].size() == displacement_gradient[k].size() ) << 
                    "Sizes of the response vector derived wrt. the displacement" <<
                    " and of the displacement vector derived wrt. the design parameter do not match!" << std::endl;
                 
                VectorMath::ScalarProduct(response_displacement_gradient[k], displacement_gradient[k], scalar_product[k]); 

                //OutputUtility::OutputOnTerminal("ScalarProduct", scalar_product[k]);

                VectorMath::SetToZero(sensitivity_vector[k]);

                VectorMath::Addition(sensitivity_vector[k], scalar_product[k]);
            }   
            if( response_sensitivity_gradient[k].size() > 0 )
            {
                KRATOS_ERROR_IF_NOT( response_sensitivity_gradient[k].size() == sensitivity_vector[k].size() ) << 
                    "Sizes of the sensitivity_vector and the response sensitivity gradient do not match!" << std::endl;
                
                //OutputUtility::OutputOnTerminal("dg/ds", response_sensitivity_gradient[k]);
                
                VectorMath::MultiplyByFactor(response_sensitivity_gradient[k], -1);

                VectorMath::Addition(sensitivity_vector[k], response_sensitivity_gradient[k]);

                OutputUtility::OutputOnTerminal("sensitivity_vector", sensitivity_vector[k]);
            }
            else
                OutputUtility::OutputOnTerminal("sensitivity_vector", sensitivity_vector[k]); 

            this->AssembleElementSensitivityContribution(rOutputVariable, sensitivity_vector[k], elem_i);  
        }

        KRATOS_CATCH("");
    }
    


    template <typename TDataType>
    void DirectSensitivityPostprocess::UpdateSensitivityOnNode(DirectSensitivityResponseFunction& Resp, Variable<TDataType> const& rResponseVariable, 
                        Variable<TDataType> const& rOutputVariable)
    {
        KRATOS_TRY;
        
        ProcessInfo& r_process_info = mrModelPart.GetProcessInfo();
        const int num_threads = 1;
        std::vector<TDataType> sensitivity_vector(num_threads);
        std::vector<std::vector<TDataType>> response_displacement_gradient(num_threads);
        std::vector<TDataType> response_sensitivity_gradient(num_threads);
        std::vector<TDataType> displacement_gradient(num_threads);
        std::vector<TDataType> scalar_product(num_threads);
        
        int k = 0;
        
        for (auto& node_i : mrModelPart.Nodes())
        {   
            std::cout << "UpdateNode: Node: "<< node_i.Id() << std::endl; 

            // Calculate derivative of the response function wrt. the displacements             
            Resp.CalculateGradient(node_i, rResponseVariable, response_displacement_gradient[k], r_process_info);

            //OutputUtility::OutputOnTerminal("dg/du", response_displacement_gradient[k]);

            // Calculate derivative of the response function wrt. the design variable
            Resp.CalculatePartialSensitivity(node_i, mrDirectSensitivityVariable, rResponseVariable,
                                                response_sensitivity_gradient[k], r_process_info);            

            // Size sensitivity vector
            
            sensitivity_vector[k].resize(0);
            scalar_product[k].resize(0);
                                           
            // Get the displacement vector derived wrt. the design parameter
            displacement_gradient[k] = node_i.FastGetSolutionStepValue(rResponseVariable);
            
            //OutputUtility::OutputOnTerminal("du/ds", displacement_gradient[k]);
                        
            if( (response_displacement_gradient[k].size() > 0) && (displacement_gradient[k].size() > 0) )
            {
                KRATOS_ERROR_IF_NOT( response_displacement_gradient[k].size() == displacement_gradient[k].size() ) << 
                    "Sizes of the response vector derived wrt. the displacement" <<
                    " and of the displacement vector derived wrt. the design parameter do not match!" << std::endl;
                 
                VectorMath::ScalarProduct(response_displacement_gradient[k], displacement_gradient[k], scalar_product[k]); 

                //OutputUtility::OutputOnTerminal("ScalarProduct", scalar_product[k]);

                VectorMath::SetToZero(sensitivity_vector[k]);

                VectorMath::Addition(sensitivity_vector[k], scalar_product[k]);
            }   
            
            if( response_sensitivity_gradient[k].size() > 0 )
            {
                KRATOS_ERROR_IF_NOT( response_sensitivity_gradient[k].size() == sensitivity_vector[k].size() ) << 
                    "Sizes of the sensitivity_vector and the response sensitivity gradient do not match!" << std::endl;
                
                //OutputUtility::OutputOnTerminal("dg/ds", response_sensitivity_gradient[k]);
                
                VectorMath::MultiplyByFactor(response_sensitivity_gradient[k], -1);

                VectorMath::Addition(sensitivity_vector[k], response_sensitivity_gradient[k]);

                OutputUtility::OutputOnTerminal("sensitivity_vector", sensitivity_vector[k]);
            }
            else
                OutputUtility::OutputOnTerminal("sensitivity_vector", sensitivity_vector[k]); 

            this->AssembleNodalSensitivityContribution(rOutputVariable, sensitivity_vector[k], node_i);  
        }



        KRATOS_CATCH("");
    }    
  

    void DirectSensitivityPostprocess::AssembleNodalSensitivityContribution(Variable<array_1d<double,3>> const& rSensitivityVariable,
                                              array_1d<double,3> const& rSensitivityVector, Node<3>& rNode)
    {
        array_1d<double,3>& r_sensitivity = rNode.FastGetSolutionStepValue(rSensitivityVariable);
        rNode.SetLock();
        for (IndexType i = 0; i < 3; ++i)
            r_sensitivity[i] += rSensitivityVector[i];
        rNode.UnSetLock();               
    }

    void DirectSensitivityPostprocess::AssembleElementSensitivityContribution(Variable<array_1d<double, 3>> const& rVariable,
                                                std::vector<array_1d<double,3>>& rSensitivityVector, Element& rElement)
    {
        ProcessInfo& r_process_info = mrModelPart.GetProcessInfo();

        rElement.SetValueOnIntegrationPoints(rVariable, rSensitivityVector, r_process_info);
        
        //rElement.GetValue(rVariable) += rSensitivityVector[0];
    }
    

    Variable<double> DirectSensitivityPostprocess::ReadScalarSensitivityVariables(std::string const& rVariableName)
    {
        const Variable<double>& r_variable = KratosComponents<Variable<double>>::Get(rVariableName);
            return r_variable;                      
    }

    Variable<array_1d<double,3>> DirectSensitivityPostprocess::ReadVectorSensitivityVariables(std::string const& rVariableName)
    {
        const Variable<array_1d<double,3>>& r_variable =
                KratosComponents<Variable<array_1d<double,3>>>::Get(rVariableName);
            return r_variable;                  
    }
    

};


