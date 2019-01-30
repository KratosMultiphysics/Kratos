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


namespace Kratos
{

    /// Constructor.
    DirectSensitivityPostprocess::DirectSensitivityPostprocess(ModelPart& rModelPart, AdjointResponseFunction& rResponseFunction, 
                                DirectSensitivityVariable& rDirectSensitivityVariable, Parameters SensitivitySettings)
      : mrModelPart(rModelPart) , mrResponseFunction(rResponseFunction) , mrDirectSensitivityVariable(rDirectSensitivityVariable)
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

        std::cout << "In direct sensitivity postprocess!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl; 

        //MBraun TODO: evalutate if gradient mode should also a memeber of this class?
        
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



    void DirectSensitivityPostprocess::UpdateSensitivities()
    {
        KRATOS_TRY;

        std::cout << "In UpdateSensitivities!!!!!!!!" << std::endl;

        if (mBuildMode == "static")
        {
            // overwrite existing.
            this->SetAllSensitivityVariablesToZero();
        }
        else
        {
            KRATOS_ERROR << "Unsupported \"build_mode\": " << mBuildMode << std::endl;
        }

        std::string output_variable_name = mDesignVariableName + std::string("_SENSITIVITY");
                
        //  Update the sensitivies

        if ( KratosComponents<Variable<double>>::Has(mDesignVariableName) )
        {   
            if ( KratosComponents<Variable<double>>::Has(output_variable_name) )
            {
                this->UpdateElementContributionToSensitivity( ReadScalarSensitivityVariables(mDesignVariableName), 
                                                        ReadScalarSensitivityVariables(output_variable_name) );
                this->UpdateConditionContributionToSensitivity( ReadScalarSensitivityVariables(mDesignVariableName), 
                                                        ReadScalarSensitivityVariables(output_variable_name) );
            }
            else
                KRATOS_ERROR << "No matching output variable for " << mDesignVariableName << "exist." << std::endl;
        }
        else if ( KratosComponents<Variable<array_1d<double,3>>>::Has(mDesignVariableName) )
        {
            if ( KratosComponents<Variable<double>>::Has(output_variable_name) )
            {
            this->UpdateElementContributionToSensitivity( ReadVectorSensitivityVariables(mDesignVariableName), 
                                                        ReadVectorSensitivityVariables(output_variable_name) );
            this->UpdateConditionContributionToSensitivity( ReadVectorSensitivityVariables(mDesignVariableName), 
                                                        ReadVectorSensitivityVariables(output_variable_name) );
            }
            else
                KRATOS_ERROR << "No matching output variable for " << mDesignVariableName << "exist." << std::endl;
        }
        else
                KRATOS_ERROR << "Unsupported variable: " <<  mDesignVariableName << "." << std::endl;

        KRATOS_CATCH("");
    }



    template <typename TDataType>
    void DirectSensitivityPostprocess::UpdateElementContributionToSensitivity(Variable<TDataType> const& rSensitivityVariable, 
                        Variable<TDataType> const& rOutputVariable)
    {
        KRATOS_TRY;
        
        ProcessInfo& r_process_info = mrModelPart.GetProcessInfo();
        const int num_threads = 1;
        std::vector<Vector> sensitivity_vector(num_threads);
        std::vector<Vector> response_displacement_gradient(num_threads);
        std::vector<Vector> response_sensitivity_gradient(num_threads);
        std::vector<Vector> displacement_gradient(num_threads);
        std::vector<Matrix> LHS(num_threads);
        std::vector<Matrix> helper_matrix(num_threads);
        int k = 0;

        
        for (auto& elem_i : mrModelPart.Elements())
        {   
            std::cout << "UpdateElementContributionToSensitivity: ELEMENT: "<< elem_i.Id() << std::endl;
            Element::GeometryType& r_geom = elem_i.GetGeometry();
            
            // No element contribution if the design variable is an condition data type
            if(mVariableType == "condition_data_type")
                break;  

            elem_i.CalculateLeftHandSide(LHS[k], r_process_info);            

            if(helper_matrix[k].size1() != 1 || helper_matrix[k].size2() != LHS[k].size2())
                helper_matrix[k].resize(1, LHS[k].size2() ,false);  

            mrResponseFunction.CalculateGradient(elem_i, LHS[k], response_displacement_gradient[k], r_process_info);            

            mrResponseFunction.CalculatePartialSensitivity(elem_i, rSensitivityVariable, helper_matrix[k], response_sensitivity_gradient[k], r_process_info);
                     
            // Get the displacement vector derived wrt. the design parameter
            elem_i.GetValuesVector(displacement_gradient[k]);
            
            // Output of dg/du und du/ds
            for(IndexType i = 0; i < displacement_gradient[k].size(); ++i)
                std::cout << "du/ds: " << displacement_gradient[k][i] << "  ::  " << "dg/du: " << response_displacement_gradient[k][i] << std::endl; 

            // Sizing of the sensitivity vector 
            if(sensitivity_vector[k].size() != 1)
                sensitivity_vector[k].resize(1, false);
            
            if( (response_displacement_gradient[k].size() > 0) && (displacement_gradient[k].size() > 0) )
            {
                KRATOS_ERROR_IF_NOT( response_displacement_gradient[k].size() == displacement_gradient[k].size() ) << 
                    "Sizes of the response vector derived wrt. the displacement" <<
                    " and of the displacement vector derived wrt. the design parameter do not match!" << std::endl;

                                
                
                double product = MathUtils<double>::Dot(response_displacement_gradient[k], displacement_gradient[k]);
                sensitivity_vector[k][0] = product;
                for(IndexType i = 0; i < sensitivity_vector[k].size(); ++i)
                    std::cout << "product: "<< sensitivity_vector[k][i] << std::endl;
            }                  

            if( response_sensitivity_gradient[k].size() > 0)
            {
                // Output of dg/ds
                for(IndexType i = 0; i < response_sensitivity_gradient.size(); ++i)  
                    std::cout << "dg/ds: "<< response_sensitivity_gradient[k][i] << std::endl;

                KRATOS_ERROR_IF_NOT( response_sensitivity_gradient[k].size() == sensitivity_vector[k].size() ) << 
                    "Sizes of the sensitivity_vector and the response sensitivity gradient do not match!" << std::endl;
                
                sensitivity_vector[k][0] += response_sensitivity_gradient[k][0];                
            } 

            // Output of sensitivity_vector
            for(IndexType i = 0; i < sensitivity_vector[k].size(); ++i)
                std::cout << "Sensitivity_vector:  "<< sensitivity_vector[k][i] << std::endl;

            // Assembling 
            if( (response_displacement_gradient[k].size() > 0) || (response_sensitivity_gradient[k].size() > 0) )
            {
                if(mVariableType == "element_data_type")
                this->AssembleElementSensitivityContribution(rOutputVariable, sensitivity_vector[k], elem_i);

                if(mVariableType == "nodal_data_type")
                this->AssembleNodalSensitivityContribution(rOutputVariable, sensitivity_vector[k], r_geom);
            }
        }

        KRATOS_CATCH("");
    }



    template <typename TDataType>
    void DirectSensitivityPostprocess::UpdateConditionContributionToSensitivity(Variable<TDataType> const& rSensitivityVariable, 
                        Variable<TDataType> const& rOutputVariable)
    {
        KRATOS_TRY;

        std::cout << "In UpdateConditionContributionToSensitivity!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
        
        ProcessInfo& r_process_info = mrModelPart.GetProcessInfo();
        const int num_threads = 1;
        std::vector<Vector> sensitivity_vector(num_threads);
        std::vector<Vector> response_displacement_gradient(num_threads);
        std::vector<Vector> response_sensitivity_gradient(num_threads);
        std::vector<Vector> displacement_gradient(num_threads);
        std::vector<Matrix> LHS(num_threads);
        std::vector<Matrix> helper_matrix(num_threads);
        int k = 0;

        for (auto& cond_i : mrModelPart.Conditions())
        {
            std::cout << "UpdateConditionContributionToSensitivity: Condition: "<< cond_i.Id() << std::endl;
            Condition::GeometryType& r_geom = cond_i.GetGeometry();
            
            // No condition contribution if the design variable is an element data type
            if(mVariableType == "element_data_type")
                break;  

            cond_i.CalculateLeftHandSide(LHS[k], r_process_info);

            if(helper_matrix[k].size1() != 1 || helper_matrix[k].size2() != LHS[k].size2())
                helper_matrix[k].resize(1, LHS[k].size2() ,false);    

            mrResponseFunction.CalculateGradient(cond_i, LHS[k], response_displacement_gradient[k], r_process_info);            

            mrResponseFunction.CalculatePartialSensitivity(cond_i, rSensitivityVariable, helper_matrix[k], response_sensitivity_gradient[k], r_process_info);
                     
            // Get the displacement vector derived wrt. the design parameter
            cond_i.GetValuesVector(displacement_gradient[k]);
            
            // Output of dg/du und du/ds                 
            for(IndexType i = 0; i < displacement_gradient[k].size(); ++i)
                std::cout << "du/ds: " << displacement_gradient[k][i] << "  ::  " << "dg/du: " << response_displacement_gradient[k][i] << std::endl; 

            // Sizing of the sensitivity vector
            if(sensitivity_vector[k].size() != 1)
                sensitivity_vector[k].resize(1, false);
            
            if( (response_displacement_gradient[k].size() > 0) && (displacement_gradient[k].size() > 0) )
            {
                KRATOS_ERROR_IF_NOT( response_displacement_gradient[k].size() == displacement_gradient[k].size() ) << 
                    "Sizes of the response vector derived wrt. the displacement" <<
                    " and of the displacement vector derived wrt. the design parameter do not match!" << std::endl;
                
                double product = MathUtils<double>::Dot(response_displacement_gradient[k], displacement_gradient[k]);
                sensitivity_vector[k][0] = product;
                for(IndexType i = 0; i < sensitivity_vector[k].size(); ++i)
                    std::cout << "product: "<< sensitivity_vector[k][i] << std::endl;
            }                  

            if( response_sensitivity_gradient[k].size() > 0)
            {
                // Output of dg/ds
                for(IndexType i = 0; i < response_sensitivity_gradient.size(); ++i)  
                    std::cout << "dg/ds: "<< response_sensitivity_gradient[k][i] << std::endl;

                KRATOS_ERROR_IF_NOT( response_sensitivity_gradient[k].size() == sensitivity_vector[k].size() ) << 
                    "Sizes of the sensitivity_vector and the response sensitivity gradient do not match!" << std::endl;
                
                sensitivity_vector[k][0] += response_sensitivity_gradient[k][0];               
            } 

            // Output of sensitivity_vector
            for(IndexType i = 0; i < sensitivity_vector[k].size(); ++i)
                std::cout << "Sensitivity_vector:  "<< sensitivity_vector[k][i] << std::endl;

            // Assembling
            if( (response_displacement_gradient[k].size() > 0) || (response_sensitivity_gradient[k].size() > 0) )
            {
                if(mVariableType == "condition_data_type")
                this->AssembleConditionSensitivityContribution(rOutputVariable, sensitivity_vector[k], cond_i);

                if(mVariableType == "nodal_data_type")
                this->AssembleNodalSensitivityContribution(rOutputVariable, sensitivity_vector[k], r_geom);
            }
        }
        
        KRATOS_CATCH("");
    }


    void DirectSensitivityPostprocess::AssembleNodalSensitivityContribution(Variable<double> const& rSensitivityVariable,
                                              Vector const& rSensitivityVector, Element::GeometryType& rGeom)
    {
        IndexType index = 0;
        for (IndexType i_node = 0; i_node < rGeom.PointsNumber(); ++i_node)
        {
            if (rGeom[i_node].GetValue(UPDATE_SENSITIVITIES) == true)
            {
                double& r_sensitivity = rGeom[i_node].FastGetSolutionStepValue(rSensitivityVariable);
                rGeom[i_node].SetLock();
                r_sensitivity += rSensitivityVector[index++];
                rGeom[i_node].UnSetLock();
            }
            else
                ++index;
        }
    }


    void DirectSensitivityPostprocess::AssembleNodalSensitivityContribution(Variable<array_1d<double, 3>> const& rSensitivityVariable,
                                              Vector const& rSensitivityVector, Element::GeometryType& rGeom)
    {
        IndexType index = 0;
        for (IndexType i_node = 0; i_node < rGeom.PointsNumber(); ++i_node)
        {
            if (rGeom[i_node].GetValue(UPDATE_SENSITIVITIES) == true)
            {
                array_1d<double, 3>& r_sensitivity =
                    rGeom[i_node].FastGetSolutionStepValue(rSensitivityVariable);
                rGeom[i_node].SetLock();
                for (IndexType d = 0; d < rGeom.WorkingSpaceDimension(); ++d)
                    r_sensitivity[d] += rSensitivityVector[index++];
                rGeom[i_node].UnSetLock();
            }
            else
                index += rGeom.WorkingSpaceDimension();
        }
    }


    void DirectSensitivityPostprocess::AssembleElementSensitivityContribution(Variable<double> const& rVariable,
                                                Vector const& rSensitivityVector, Element& rElement)
    {
        std::cout << "In AssembleElementSensitivityContribution!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;

        KRATOS_DEBUG_ERROR_IF(rSensitivityVector.size() != 1) << "rSensitivityVector.size() = " << rSensitivityVector.size() << std::endl;
        rElement.GetValue(rVariable) += rSensitivityVector[0];
    }


    void DirectSensitivityPostprocess::AssembleElementSensitivityContribution(Variable<array_1d<double, 3>> const& rVariable,
                                                Vector const& rSensitivityVector, Element& rElement)
    {
        array_1d<double, 3>& r_sensitivity = rElement.GetValue(rVariable);
        const auto ws_dim = rElement.GetGeometry().WorkingSpaceDimension();
        KRATOS_DEBUG_ERROR_IF(rSensitivityVector.size() != ws_dim) << "rSensitivityVector.size() = " << rSensitivityVector.size() << std::endl;
        for (unsigned d = 0; d < ws_dim; ++d)
            r_sensitivity[d] += rSensitivityVector[d];
    }


    void DirectSensitivityPostprocess::AssembleConditionSensitivityContribution(Variable<double> const& rSensitivityVariable,
                                              Vector const& rSensitivityVector,
                                              Condition& rCondition)
    {
        KRATOS_DEBUG_ERROR_IF(rSensitivityVector.size() != 1) << "rSensitivityVector.size() = " << rSensitivityVector.size() << std::endl;
        rCondition.GetValue(rSensitivityVariable) += rSensitivityVector[0];
    }


    void DirectSensitivityPostprocess::AssembleConditionSensitivityContribution(Variable<array_1d<double, 3>> const& rSensitivityVariable,
                                              Vector const& rSensitivityVector,
                                              Condition& rCondition)
    {
        array_1d<double, 3>& r_sensitivity = rCondition.GetValue(rSensitivityVariable);
        const auto ws_dim = rCondition.GetGeometry().WorkingSpaceDimension();
        KRATOS_DEBUG_ERROR_IF(rSensitivityVector.size() != ws_dim) << "rSensitivityVector.size() = " << rSensitivityVector.size() << std::endl;
        for (unsigned d = 0; d < ws_dim; ++d)
            r_sensitivity[d] += rSensitivityVector[d];
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


