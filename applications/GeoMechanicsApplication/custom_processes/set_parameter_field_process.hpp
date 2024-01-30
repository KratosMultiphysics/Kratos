//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Aron Noordam
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"

#include "utilities/function_parser_utility.h"
#include "utilities/mortar_utilities.h"
namespace Kratos {

    ///@name Kratos Globals
    ///@{

    ///@name Kratos Classes
    ///@{


    /**
     * @class SetParameterFieldProcess
     * @ingroup GeoMechanicsApplication
     * @brief Process to set a parameter field
     * @details This process sets values from a custom parameter field per individual element within the model part. the possibilities are: setting a parameter field from an
     * input function; setting a parameter field with a user defined python script; setting a parameter field from a json file.
     * @author Aron Noordam
    */
    class KRATOS_API(GEO_MECHANICS_APPLICATION) SetParameterFieldProcess : public Process
    {
    public:
        ///@name Type Definitions
        ///@{

        using SizeType = std::size_t;

        ///@}
        ///@name Pointer Definitions
        /// Pointer definition of SetMovingLoadProcess
        KRATOS_CLASS_POINTER_DEFINITION(SetParameterFieldProcess);

        ///@}
        ///@name Life Cycle
        ///@{

        SetParameterFieldProcess(ModelPart & rModelPart,
            const Parameters & rParameters);

        ///@}
        ///@name Operations
        ///@{

        /**
         * \brief  Initializes the set parameter field process. Checks if the value that needs to be changed is a UMAT_PARAMETERS(so Vector) or double.
         */
        void ExecuteInitialize() override;


        ///@}
        ///@name Input and output
        ///@{

        /// Turn back information as a string.
        std::string Info() const override {
            return "SetParameterFieldProcess";
        }
        ///@}

    private:
        ///@name Member Variables
        ///@{

        ModelPart& mrModelPart;
        Parameters mParameters;


        ///@}
        ///@name Private Operations
        ///@{

        /**
         * \brief Sets the value of a certain variable for only 1 element
         * \param rElement element for which a variable should be altered
         * \param rVar variable type which is to be altered
         * \param Value new value for the to be altered variable
         */
        template <typename T>
        static void SetValueAtElement(Element & rElement, const Variable<T>&rVar, const T & Value)
        {
            Properties& r_prop = rElement.GetProperties();

            // Copies properties
            Properties::Pointer p_new_prop = Kratos::make_shared<Properties>(r_prop);

            // Adds new properties to the element
            p_new_prop->SetValue(rVar, Value);
            rElement.SetProperties(p_new_prop);
        }

        /**
        * \brief Sets the parameter field, using an input function
        * \param rVar variable type which is to be used to generate the parameter field
        */
        void SetParameterFieldUsingInputFunction(const Variable<double>&rVar);
        void SetParameterFieldUsingInputFunction(const Variable<Vector>&rVar);

        /**
        * \brief Sets the parameter field, using a custom Parameter structure
        * \param rVar variable type which is to be used to generate the parameter field
        */
        void SetParameterFieldUsingParametersClass(const Variable<double>&rVar, const Parameters & rParameters);
        void SetParameterFieldUsingParametersClass(const Variable<Vector>&rVar, const Parameters & rParameters);

        /**
        * \brief Sets the parameter field, using an input json file
        * \param rVar variable type which is to be used to generate the parameter field
        */
        template <typename T>
        void SetParameterFieldUsingJsonFile(const Variable<T>&rVar)
        {
            // Read json string in field parameters file, create Parameters
            const std::string& field_file_name = mParameters["dataset_file_name"].GetString();
            KRATOS_ERROR_IF_NOT(std::filesystem::exists(field_file_name))
                << "The parameter field file specified with name \"" << field_file_name << "\" does not exist!" << std::endl;

            std::ifstream ifs(field_file_name);
            this->SetParameterFieldUsingParametersClass(rVar, Parameters{ ifs });
        }

        /**
        * \brief Sets the parameter field, using a json string
        * \param rVar variable type which is to be used to generate the parameter field
        */
        template <typename T>
        void SetParameterFieldUsingJsonString(const Variable<T>&rVar)
        {
            // get new data from the data set
            const std::string& r_dataset = mParameters["dataset"].GetString();
            Parameters new_data{ r_dataset };
            this->SetParameterFieldUsingParametersClass(rVar, new_data);
        }

        /**
        * \brief  Checks what type of input field is given and generates the parameter field.
        * \param rVar variable type which is to be used to generate the parameter field
        */
        template <typename T>
        void SetParameterFieldForVariableType(const Variable<T>&r_var)
        {
            if (mParameters["func_type"].GetString() == "input") {
                this->SetParameterFieldUsingInputFunction(r_var);
            }
            else if (mParameters["func_type"].GetString() == "json_string") {
                this->SetParameterFieldUsingJsonString(r_var);
            }
            else if (mParameters["func_type"].GetString() == "json_file") {
                this->SetParameterFieldUsingJsonFile(r_var);
            }
        }

        std::vector<IndexType> GetVectorIndices() const;

        ///@}
        ///@name Serialization
        ///@{

        ///@}

    }; // Class SetParameterFieldProcess

    ///@}

}  // namespace Kratos.