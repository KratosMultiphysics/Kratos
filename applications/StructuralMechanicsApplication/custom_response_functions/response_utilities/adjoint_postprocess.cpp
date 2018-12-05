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
#include "adjoint_postprocess.h"
#include "utilities/variable_utils.h"
#include "adjoint_postprocess.h"


namespace Kratos
{

    /// Constructor.
    AdjointPostprocess::AdjointPostprocess(ModelPart& rModelPart, AdjointStructuralResponseFunction& rResponseFunction, Parameters SensitivitySettings)
      : mrModelPart(rModelPart) , mrResponseFunction(rResponseFunction)
    {
        KRATOS_TRY;

        Parameters default_settings(R"(
        {
            "sensitivity_model_part_name": "PLEASE_SPECIFY_SENSITIVITY_MODEL_PART",
            "nodal_sensitivity_variables": [],
            "element_sensitivity_variables": [],
            "condition_sensitivity_variables": [],
            "build_mode": "static"
        })");

        SensitivitySettings.ValidateAndAssignDefaults(default_settings);

        auto sensitivity_model_part_name =
            SensitivitySettings["sensitivity_model_part_name"].GetString();
        if (sensitivity_model_part_name !=
            "PLEASE_SPECIFY_SENSITIVITY_MODEL_PART")
        {
            mpSensitivityModelPart =
                &mrModelPart.GetSubModelPart(sensitivity_model_part_name);
        }
        else
        {
            mpSensitivityModelPart = &mrModelPart;
        }

        mBuildMode = SensitivitySettings["build_mode"].GetString();

        this->ReadDesignVariables(mNodalSensitivityScalarVariables, mNodalSensitivityVectorVariables, SensitivitySettings["nodal_sensitivity_variables"]);
        this->ReadDesignVariables(mElementSensitivityScalarVariables, mElementSensitivityVectorVariables, SensitivitySettings["element_sensitivity_variables"]);
        this->ReadDesignVariables(mConditionSensitivityScalarVariables, mConditionSensitivityVectorVariables, SensitivitySettings["condition_sensitivity_variables"]);

        //MFusseder TODO: evalutate if gradient mode should also a memeber of this class?
        KRATOS_CATCH("");
    }

    /// Destructor.
    AdjointPostprocess::~AdjointPostprocess()
    {
    }

    void AdjointPostprocess::Initialize()
    {
        KRATOS_TRY;

        this->Clear();

        // Initialize flags.
        VariableUtils().SetNonHistoricalVariable(UPDATE_SENSITIVITIES, true, mpSensitivityModelPart->Nodes());
        VariableUtils().SetNonHistoricalVariable(UPDATE_SENSITIVITIES, true, mpSensitivityModelPart->Elements());
        VariableUtils().SetNonHistoricalVariable(UPDATE_SENSITIVITIES, true, mpSensitivityModelPart->Conditions());

        KRATOS_CATCH("");
    }

    void AdjointPostprocess::Clear()
    {
        KRATOS_TRY;

        // Reset flags.
        VariableUtils().SetNonHistoricalVariable(UPDATE_SENSITIVITIES, false, mrModelPart.Nodes());
        VariableUtils().SetNonHistoricalVariable(UPDATE_SENSITIVITIES, false, mrModelPart.Elements());
        VariableUtils().SetNonHistoricalVariable(UPDATE_SENSITIVITIES, false, mrModelPart.Conditions());

        this->SetAllSensitivityVariablesToZero();

        KRATOS_CATCH("");
    }

    void AdjointPostprocess::SetAllSensitivityVariablesToZero()
    {
        KRATOS_TRY;

        // Set nodal sensitivity result variables to zero.
        for (const auto& variable_pair : mNodalSensitivityScalarVariables)
            VariableUtils().SetToZero_ScalarVar(variable_pair[1], mrModelPart.Nodes());
        for (const auto& variable_pair : mNodalSensitivityVectorVariables)
            VariableUtils().SetToZero_VectorVar(variable_pair[1], mrModelPart.Nodes());
        // Set elemental sensitivity result variables to zero.
        for (const auto& variable_pair : mElementSensitivityScalarVariables)
        {
            #pragma omp parallel for
            for (int i = 0; i< static_cast<int> (mrModelPart.NumberOfElements()); ++i)
            {
                auto it = mrModelPart.ElementsBegin() + i;
                it->SetValue(variable_pair[1], variable_pair[1].Zero());
            }
        }
        for (const auto& variable_pair : mElementSensitivityVectorVariables)
        {
             #pragma omp parallel for
            for (int i = 0; i< static_cast<int> (mrModelPart.NumberOfElements()); ++i)
            {
                auto it = mrModelPart.ElementsBegin() + i;
                it->SetValue(variable_pair[1], variable_pair[1].Zero());
            }
        }
        // Set conditional sensitivity result variables to zero.
        for (const auto& variable_pair : mConditionSensitivityScalarVariables)
        {
            #pragma omp parallel for
            for (int i = 0; i< static_cast<int> (mrModelPart.NumberOfConditions()); ++i)
            {
                auto it = mrModelPart.ConditionsBegin() + i;
                const SizeType number_of_nodes = it->GetGeometry().size();
                for(IndexType j = 0; j < number_of_nodes; ++j)
                    it->GetGeometry()[j].FastGetSolutionStepValue(variable_pair[1]) = variable_pair[1].Zero();
            }
        }
        for (const auto& variable_pair : mConditionSensitivityVectorVariables)
        {
            #pragma omp parallel for
            for (int i = 0; i< static_cast<int> (mrModelPart.NumberOfConditions()); ++i)
            {
                auto it = mrModelPart.ConditionsBegin() + i;
                const SizeType number_of_nodes = it->GetGeometry().size();
                for(IndexType j = 0; j < number_of_nodes; ++j)
                    it->GetGeometry()[j].FastGetSolutionStepValue(variable_pair[1]) = variable_pair[1].Zero();
            }
        }

        KRATOS_CATCH("");
    }

    void AdjointPostprocess::UpdateSensitivities()
    {
        KRATOS_TRY;

        if (mBuildMode == "static")
        {
            // overwrite existing.
            this->SetAllSensitivityVariablesToZero();
        }
        else
        {
            KRATOS_ERROR << "Unsupported \"build_mode\": " << mBuildMode << std::endl;
        }

        for (const auto& variable_pair : mNodalSensitivityScalarVariables)
            this->UpdateNodalSensitivities(variable_pair[0], variable_pair[1]);
        for (const auto& variable_pair : mNodalSensitivityVectorVariables)
            this->UpdateNodalSensitivities(variable_pair[0], variable_pair[1]);
        for (const auto& variable_pair : mElementSensitivityScalarVariables)
            this->UpdateElementSensitivities(variable_pair[0], variable_pair[1]);
        for (const auto& variable_pair : mElementSensitivityVectorVariables)
            this->UpdateElementSensitivities(variable_pair[0], variable_pair[1]);
        for (const auto& variable_pair : mConditionSensitivityScalarVariables)
            this->UpdateConditionSensitivities(variable_pair[0], variable_pair[1]);
        for (const auto& variable_pair : mConditionSensitivityVectorVariables)
            this->UpdateConditionSensitivities(variable_pair[0], variable_pair[1]);

        KRATOS_CATCH("");
    }

    template <typename TDataType>
    void AdjointPostprocess::UpdateNodalSensitivities(Variable<TDataType> const& rSensitivityVariable, Variable<TDataType> const& rOutputVariable)
    {
        KRATOS_TRY;

        this->BuildNodalSolutionStepElementContributions(rSensitivityVariable, rOutputVariable);

        this->BuildNodalSolutionStepConditionContributions(rSensitivityVariable, rOutputVariable);

        mrModelPart.GetCommunicator().AssembleCurrentData(rSensitivityVariable);

        KRATOS_CATCH("");
    }

    template <typename TDataType>
    void AdjointPostprocess::BuildNodalSolutionStepElementContributions(Variable<TDataType> const& rSensitivityVariable, Variable<TDataType> const& rOutputVariable)
    {
        KRATOS_TRY;

        ProcessInfo& r_process_info = mrModelPart.GetProcessInfo();
        const int num_threads = 1;
        std::vector<Vector> sensitivity_vector(num_threads);
        std::vector<Vector> response_gradient(num_threads);
        std::vector<Vector> adjoint_vector(num_threads);
        std::vector<Matrix> sensitivity_matrix(num_threads);

        int k = 0;

        for (auto& elem_i : mrModelPart.Elements())
        {
            Element::GeometryType& r_geom = elem_i.GetGeometry();
            bool update_sensitivities = false;
            for (IndexType i_node = 0; i_node < r_geom.PointsNumber(); ++i_node)
                if (r_geom[i_node].GetValue(UPDATE_SENSITIVITIES))
                {
                    update_sensitivities = true;
                    break;
                }

            if (!update_sensitivities) // true for most elements
                continue;

            // Compute the pseudo load
            elem_i.CalculateSensitivityMatrix(
                rSensitivityVariable, sensitivity_matrix[k], r_process_info);

            // This part of the sensitivity is computed from the objective
            // with primal variables treated as constant.
            mrResponseFunction.CalculateSensitivityGradient(
                elem_i, rSensitivityVariable, sensitivity_matrix[k],
                response_gradient[k], r_process_info);

            if( (response_gradient[k].size() > 0) && (sensitivity_matrix[k].size1() > 0) )
            {
                KRATOS_ERROR_IF_NOT( response_gradient[k].size() ==
                    sensitivity_matrix[k].size1() ) << "Sizes of sensitivity" <<
                        "matrix and response gradient do not match!" << std::endl;
            }

            if(sensitivity_matrix[k].size1() > 0)
            {
                if (sensitivity_vector[k].size() != sensitivity_matrix[k].size1())
                    sensitivity_vector[k].resize(sensitivity_matrix[k].size1(), false);

                // Get the adjoint displacement field
                elem_i.GetValuesVector(adjoint_vector[k]);

                // Compute the adjoint variable times the sensitivity_matrix (pseudo load)
                noalias(sensitivity_vector[k]) = prod(sensitivity_matrix[k], adjoint_vector[k]);
            }

            if(response_gradient[k].size() > 0)
            {
                if (sensitivity_vector[k].size() != response_gradient[k].size())
                    sensitivity_vector[k].resize(response_gradient[k].size(), false);

                // Add the partial response gradient
                noalias(sensitivity_vector[k]) += response_gradient[k];
            }

            if( (response_gradient[k].size() > 0) || (sensitivity_matrix[k].size1() > 0) )
            {
                this->AssembleNodalSensitivityContribution(
                    rOutputVariable, sensitivity_vector[k], r_geom);
            }
        }

        KRATOS_CATCH("");
    }
    template <typename TDataType>
    void AdjointPostprocess::BuildNodalSolutionStepConditionContributions(Variable<TDataType> const& rSensitivityVariable, Variable<TDataType> const& rOutputVariable)
    {
        KRATOS_TRY;

        ProcessInfo& r_process_info = mrModelPart.GetProcessInfo();
        const int num_threads = 1;
        std::vector<Vector> sensitivity_vector(num_threads);
        std::vector<Vector> response_gradient(num_threads);
        std::vector<Vector> adjoint_vector(num_threads);
        std::vector<Matrix> sensitivity_matrix(num_threads);

        int k = 0;

        // Assemble condition contributions.
        for (auto& cond_i : mrModelPart.Conditions())
        {
            Condition::GeometryType& r_geom = cond_i.GetGeometry();
            bool update_sensitivities = false;
            for (IndexType i_node = 0; i_node < r_geom.PointsNumber(); ++i_node)
                if (r_geom[i_node].GetValue(UPDATE_SENSITIVITIES) == true)
                {
                    update_sensitivities = true;
                    break;
                }

            if (update_sensitivities == false)
                continue;

            // This is multiplied with the adjoint to compute sensitivity
            // contributions from the condition.
            cond_i.CalculateSensitivityMatrix(
                rSensitivityVariable, sensitivity_matrix[k], r_process_info);

            // This part of the sensitivity is computed from the objective
            // with primal variables treated as constant.
            mrResponseFunction.CalculateSensitivityGradient(
                cond_i, rSensitivityVariable, sensitivity_matrix[k],
                response_gradient[k], r_process_info);

            if( (response_gradient[k].size() > 0) && (sensitivity_matrix[k].size1() > 0) )
            {
                KRATOS_ERROR_IF_NOT( response_gradient[k].size() ==
                    sensitivity_matrix[k].size1() ) << "Sizes of sensitivity" <<
                        "matrix and response gradient do not match!" << std::endl;
            }

            if(sensitivity_matrix[k].size1() > 0)
            {
                if (sensitivity_vector[k].size() != sensitivity_matrix[k].size1())
                    sensitivity_vector[k].resize(sensitivity_matrix[k].size1(), false);

                // Get the adjoint displacement field
                cond_i.GetValuesVector(adjoint_vector[k]);

                // Compute the adjoint variable times the sensitivity_matrix (pseudo load)
                noalias(sensitivity_vector[k]) = prod(sensitivity_matrix[k], adjoint_vector[k]);
            }

            if(response_gradient[k].size() > 0)
            {
                if (sensitivity_vector[k].size() != response_gradient[k].size())
                    sensitivity_vector[k].resize(response_gradient[k].size(), false);

                // Add the partial response gradient
                noalias(sensitivity_vector[k]) += response_gradient[k];
            }

            if( (response_gradient[k].size() > 0) || (sensitivity_matrix[k].size1() > 0) )
                this->AssembleNodalSensitivityContribution(rOutputVariable, sensitivity_vector[k], r_geom);
        }

        KRATOS_CATCH("");
    }

    template <typename TDataType>
    void AdjointPostprocess::UpdateElementSensitivities(Variable<TDataType> const& rSensitivityVariable, Variable<TDataType> const& rOutputVariable)
    {
        KRATOS_TRY;

        ProcessInfo& r_process_info = mrModelPart.GetProcessInfo();
        const int num_threads = OpenMPUtils::GetNumThreads();
        std::vector<Vector> sensitivity_vector(num_threads);
        std::vector<Vector> response_gradient(num_threads);
        std::vector<Vector> adjoint_vector(num_threads);
        std::vector<Matrix> sensitivity_matrix(num_threads);

        #pragma omp parallel for
        for (int i = 0; i< static_cast<int> (mrModelPart.NumberOfElements()); ++i)
        {
            const unsigned int k = OpenMPUtils::ThisThread();
            auto it = mrModelPart.ElementsBegin() + i;

            if (!(it->GetValue(UPDATE_SENSITIVITIES)))
                continue;

            // Compute the pseudo load
            it->CalculateSensitivityMatrix(
                rSensitivityVariable, sensitivity_matrix[k], r_process_info);

            // This part of the sensitivity is computed from the objective
            // with primal variables treated as constant.
            mrResponseFunction.CalculateSensitivityGradient(
                *it, rSensitivityVariable, sensitivity_matrix[k],
                    response_gradient[k], r_process_info);

            if( (response_gradient[k].size() > 0) && (sensitivity_matrix[k].size1() > 0) )
            {
                KRATOS_ERROR_IF_NOT( response_gradient[k].size() ==
                    sensitivity_matrix[k].size1() ) << "Sizes of sensitivity" <<
                        "matrix and response gradient do not match!" << std::endl;
            }

            if(sensitivity_matrix[k].size1() > 0)
            {
                if (sensitivity_vector[k].size() != sensitivity_matrix[k].size1())
                    sensitivity_vector[k].resize(sensitivity_matrix[k].size1(), false);

                // Get the adjoint displacement field
                it->GetValuesVector(adjoint_vector[k]);

                // Compute the adjoint variable times the sensitivity_matrix (pseudo load)
                noalias(sensitivity_vector[k]) = prod(sensitivity_matrix[k], adjoint_vector[k]) ;
            }

            if(response_gradient[k].size() > 0)
            {
                if (sensitivity_vector[k].size() != response_gradient[k].size())
                    sensitivity_vector[k].resize(response_gradient[k].size(), false);

                // Add the partial response gradient
                noalias(sensitivity_vector[k]) += response_gradient[k];
            }

            if( (response_gradient[k].size() > 0) || (sensitivity_matrix[k].size1() > 0) )
            {
                this->AssembleElementSensitivityContribution(
                            rOutputVariable, sensitivity_vector[k], *it);
            }
        }

        mrModelPart.GetCommunicator().AssembleCurrentData(rSensitivityVariable);

        KRATOS_CATCH("");
    }

    template <typename TDataType>
    void AdjointPostprocess::UpdateConditionSensitivities(Variable<TDataType> const& rSensitivityVariable, Variable<TDataType> const& rOutputVariable)
    {
        KRATOS_TRY;

        ProcessInfo& r_process_info = mrModelPart.GetProcessInfo();
        const int num_threads = OpenMPUtils::GetNumThreads();
        std::vector<Vector> sensitivity_vector(num_threads);
        std::vector<Vector> response_gradient(num_threads);
        std::vector<Vector> adjoint_vector(num_threads);
        std::vector<Matrix> sensitivity_matrix(num_threads);

        //  Assemble condition contributions.
        #pragma omp parallel for
        for (int i = 0; i< static_cast<int> (mrModelPart.NumberOfConditions()); ++i)
        {
            const unsigned int k = OpenMPUtils::ThisThread();
            auto it = mrModelPart.ConditionsBegin() + i;

            if (!(it->GetValue(UPDATE_SENSITIVITIES)))
                continue;

            // Compute the pseudo load
            it->CalculateSensitivityMatrix(
                rSensitivityVariable, sensitivity_matrix[k], r_process_info);

            // This part of the sensitivity is computed from the objective
            // with primal variables treated as constant.
            mrResponseFunction.CalculateSensitivityGradient(
                *it, rSensitivityVariable, sensitivity_matrix[k],
                response_gradient[k], r_process_info);

            if( (response_gradient[k].size() > 0) && (sensitivity_matrix[k].size1() > 0) )
            {
                KRATOS_ERROR_IF_NOT( response_gradient[k].size() ==
                    sensitivity_matrix[k].size1() ) << "Sizes of sensitivity" <<
                        "matrix and response gradient do not match!" << std::endl;
            }

            if(sensitivity_matrix[k].size1() > 0)
            {
                if (sensitivity_vector[k].size() != sensitivity_matrix[k].size1())
                    sensitivity_vector[k].resize(sensitivity_matrix[k].size1(), false);

                // Get the adjoint displacement field
                it->GetValuesVector(adjoint_vector[k]);

                // Compute the adjoint variable times the sensitivity_matrix (pseudo load)
                noalias(sensitivity_vector[k]) = prod(sensitivity_matrix[k], adjoint_vector[k]);
            }

            if(response_gradient[k].size() > 0)
            {
                if (sensitivity_vector[k].size() != response_gradient[k].size())
                    sensitivity_vector[k].resize(response_gradient[k].size(), false);

                // Add the partial response gradient
                noalias(sensitivity_vector[k]) += response_gradient[k];
            }

            if( (response_gradient[k].size() > 0) || (sensitivity_matrix[k].size1() > 0) )
            {
                Condition::GeometryType& r_geom = it->GetGeometry();
                this->AssembleConditionSensitivityContribution(
                            rOutputVariable, sensitivity_vector[k], r_geom);
            }
        }

        mrModelPart.GetCommunicator().AssembleCurrentData(rSensitivityVariable);

        KRATOS_CATCH("");
    }

    void AdjointPostprocess::AssembleNodalSensitivityContribution(Variable<double> const& rSensitivityVariable,
                                              Vector const& rSensitivityVector,
                                              Element::GeometryType& rGeom)
    {
        IndexType index = 0;
        for (IndexType i_node = 0; i_node < rGeom.PointsNumber(); ++i_node)
        {
            if (rGeom[i_node].GetValue(UPDATE_SENSITIVITIES) == true)
            {
                double& r_sensitivity =
                    rGeom[i_node].FastGetSolutionStepValue(rSensitivityVariable);
                rGeom[i_node].SetLock();
                r_sensitivity += rSensitivityVector[index++];
                rGeom[i_node].UnSetLock();
            }
            else
                ++index;
        }
    }

    void AdjointPostprocess::AssembleNodalSensitivityContribution(Variable<array_1d<double, 3>> const& rSensitivityVariable,
                                              Vector const& rSensitivityVector,
                                              Element::GeometryType& rGeom)
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

    void AdjointPostprocess::AssembleElementSensitivityContribution(Variable<double> const& rVariable,
                                                Vector const& rSensitivityVector,
                                                Element& rElement)
    {
        KRATOS_DEBUG_ERROR_IF(rSensitivityVector.size() != 1) << "rSensitivityVector.size() = " << rSensitivityVector.size() << std::endl;
        rElement.GetValue(rVariable) += rSensitivityVector[0];
    }

    void AdjointPostprocess::AssembleElementSensitivityContribution(Variable<array_1d<double, 3>> const& rVariable,
                                                Vector const& rSensitivityVector,
                                                Element& rElement)
    {
        array_1d<double, 3>& r_sensitivity = rElement.GetValue(rVariable);
        const auto ws_dim = rElement.GetGeometry().WorkingSpaceDimension();
        KRATOS_DEBUG_ERROR_IF(rSensitivityVector.size() != ws_dim) << "rSensitivityVector.size() = " << rSensitivityVector.size() << std::endl;
        for (unsigned d = 0; d < ws_dim; ++d)
            r_sensitivity[d] += rSensitivityVector[d];
    }

    void AdjointPostprocess::AssembleConditionSensitivityContribution(Variable<double> const& rSensitivityVariable,
                                              Vector const& rSensitivityVector,
                                              Element::GeometryType& rGeom)
    {
        IndexType index = 0;
        for (IndexType i_node = 0; i_node < rGeom.PointsNumber(); ++i_node)
        {
            double& r_sensitivity =
                rGeom[i_node].FastGetSolutionStepValue(rSensitivityVariable);
            rGeom[i_node].SetLock();
            r_sensitivity += rSensitivityVector[index++];
            rGeom[i_node].UnSetLock();
        }
    }

    void AdjointPostprocess::AssembleConditionSensitivityContribution(Variable<array_1d<double, 3>> const& rSensitivityVariable,
                                              Vector const& rSensitivityVector,
                                              Element::GeometryType& rGeom)
    {
        IndexType index = 0;
        for (IndexType i_node = 0; i_node < rGeom.PointsNumber(); ++i_node)
        {
            array_1d<double, 3>& r_sensitivity =
                rGeom[i_node].FastGetSolutionStepValue(rSensitivityVariable);
            rGeom[i_node].SetLock();
            for (IndexType d = 0; d < rGeom.WorkingSpaceDimension(); ++d)
                r_sensitivity[d] += rSensitivityVector[index++];
            rGeom[i_node].UnSetLock();
        }
    }

    void AdjointPostprocess::ReadDesignVariables(std::vector<std::vector<Variable<double>>>& rScalarDesignVariables,
        std::vector<std::vector<Variable<array_1d<double,3>>>>& rVectorDesignVariables, Parameters DesignVariableSettings)
    {
        for (IndexType i = 0; i < DesignVariableSettings.size(); ++i)
        {
            const std::string variable_label = DesignVariableSettings[i].GetString();
            const std::string output_variable_label = variable_label + "_SENSITIVITY";
            std::vector<Variable<double>> helper_scalar_variables;
            std::vector<Variable<array_1d<double,3>>> helper_vector_variables;

            if (KratosComponents<Variable<double>>::Has(variable_label))
            {
                const Variable<double>& r_variable =
                    KratosComponents<Variable<double>>::Get(variable_label);

                helper_scalar_variables.push_back(r_variable);

                if (KratosComponents<Variable<double>>::Has(output_variable_label))
                {
                    const Variable<double>& r_output_variable =
                        KratosComponents<Variable<double>>::Get(output_variable_label);
                    helper_scalar_variables.push_back(r_output_variable);
                }
                else
                    KRATOS_ERROR << "Unsupported output variable: " << output_variable_label << "." << std::endl;

                rScalarDesignVariables.push_back(helper_scalar_variables);
            }
            else if (KratosComponents<Variable<array_1d<double,3>>>::Has(variable_label))
            {
                const Variable<array_1d<double, 3>>& r_variable =
                    KratosComponents<Variable<array_1d<double, 3>>>::Get(variable_label);

                helper_vector_variables.push_back(r_variable);

                if (KratosComponents<Variable<array_1d<double,3>>>::Has(output_variable_label))
                {
                    const Variable<array_1d<double, 3>>& r_output_variable =
                        KratosComponents<Variable<array_1d<double, 3>>>::Get(output_variable_label);
                    helper_vector_variables.push_back(r_output_variable);
                }
                else
                    KRATOS_ERROR << "Unsupported output variable: " << output_variable_label << "." << std::endl;

                rVectorDesignVariables.push_back(helper_vector_variables);
            }
            else
                KRATOS_ERROR << "Unsupported variable: " << variable_label << "." << std::endl;
        }
    }

};


