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

        mSensitivityModelPartName = SensitivitySettings["sensitivity_model_part_name"].GetString();

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

    ModelPart& AdjointPostprocess::GetModelPart()
    {
      return mrModelPart;
    }

    ModelPart& AdjointPostprocess::GetModelPart() const
    {
      return mrModelPart;
    }

    void AdjointPostprocess::Initialize()
    {
        KRATOS_TRY;

        this->Clear();
        this->Check();

        ModelPart& r_model_part = this->GetModelPart();
        ModelPart& r_sensitivity_model_part = r_model_part.GetSubModelPart(mSensitivityModelPartName);

        // Initialize flags.
        VariableUtils().SetNonHistoricalVariable(UPDATE_SENSITIVITIES, true, r_sensitivity_model_part.Nodes());
        VariableUtils().SetNonHistoricalVariable(UPDATE_SENSITIVITIES, true, r_sensitivity_model_part.Elements());
        VariableUtils().SetNonHistoricalVariable(UPDATE_SENSITIVITIES, true, r_sensitivity_model_part.Conditions());

        KRATOS_CATCH("");
    }

    void AdjointPostprocess::Check()
    {
        KRATOS_TRY;

        ModelPart& r_model_part = this->GetModelPart();

        KRATOS_ERROR_IF_NOT(r_model_part.HasSubModelPart(mSensitivityModelPartName))
            << "No sub model part \"" << mSensitivityModelPartName << "\"" << std::endl;

        KRATOS_CATCH("");
    }

    void AdjointPostprocess::Clear()
    {
        KRATOS_TRY;

        ModelPart& r_model_part = this->GetModelPart();

        // Reset flags.
        VariableUtils().SetNonHistoricalVariable(UPDATE_SENSITIVITIES, false, r_model_part.Nodes());
        VariableUtils().SetNonHistoricalVariable(UPDATE_SENSITIVITIES, false, r_model_part.Elements());
        VariableUtils().SetNonHistoricalVariable(UPDATE_SENSITIVITIES, false, r_model_part.Conditions());

        this->SetAllSensitivityVariablesToZero();

        KRATOS_CATCH("");
    }

    void AdjointPostprocess::SetAllSensitivityVariablesToZero()
    {
        KRATOS_TRY;

        ModelPart& r_model_part = this->GetModelPart();

        // Set nodal sensitivity result variables to zero.
        for (const auto& variable_pair : mNodalSensitivityScalarVariables)
            VariableUtils().SetToZero_ScalarVar(variable_pair[1], r_model_part.Nodes());
        for (const auto& variable_pair : mNodalSensitivityVectorVariables)
            VariableUtils().SetToZero_VectorVar(variable_pair[1], r_model_part.Nodes());
        // Set elemental sensitivity result variables to zero.
        for (const auto& variable_pair : mElementSensitivityScalarVariables)
        {
            #pragma omp parallel for
            for (int i = 0; i< static_cast<int> (r_model_part.NumberOfElements()); ++i)
            {
                auto it = r_model_part.ElementsBegin() + i;
                it->SetValue(variable_pair[1], variable_pair[1].Zero());
            }
        }
        for (const auto& variable_pair : mElementSensitivityVectorVariables)
        {
             #pragma omp parallel for
            for (int i = 0; i< static_cast<int> (r_model_part.NumberOfElements()); ++i)
            {
                auto it = r_model_part.ElementsBegin() + i;
                it->SetValue(variable_pair[1], variable_pair[1].Zero());
            }
        }
        // Set conditional sensitivity result variables to zero.
        for (const auto& variable_pair : mConditionSensitivityScalarVariables)
        {
            #pragma omp parallel for
            for (int i = 0; i< static_cast<int> (r_model_part.NumberOfConditions()); ++i)
            {
                auto it = r_model_part.ConditionsBegin() + i;
                const SizeType number_of_nodes = it->GetGeometry().size();
                for(IndexType j = 0; j < number_of_nodes; ++j)
                    it->GetGeometry()[j].FastGetSolutionStepValue(variable_pair[1]) = variable_pair[1].Zero();
            }
        }
        for (const auto& variable_pair : mConditionSensitivityVectorVariables)
        {
            #pragma omp parallel for
            for (int i = 0; i< static_cast<int> (r_model_part.NumberOfConditions()); ++i)
            {
                auto it = r_model_part.ConditionsBegin() + i;
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

        ModelPart& r_model_part = this->GetModelPart();
        ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
        const int num_threads = 1;
        std::vector<Vector> sensitivity_vector(num_threads);
        std::vector<Vector> response_gradient(num_threads);
        std::vector<Vector> adjoint_vector(num_threads);
        std::vector<Matrix> sensitivity_matrix(num_threads);

        int k = 0;

        for (auto& elem_i : r_model_part.Elements())
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

        // Assemble condition contributions.
        for (auto& cond_i : r_model_part.Conditions())
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

        r_model_part.GetCommunicator().AssembleCurrentData(rSensitivityVariable);

        KRATOS_CATCH("");
    }

    template <typename TDataType>
    void AdjointPostprocess::UpdateElementSensitivities(Variable<TDataType> const& rSensitivityVariable, Variable<TDataType> const& rOutputVariable)
    {
        KRATOS_TRY;

        ModelPart& r_model_part = this->GetModelPart();
        ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
        const int num_threads = OpenMPUtils::GetNumThreads();
        std::vector<Vector> sensitivity_vector(num_threads);
        std::vector<Vector> response_gradient(num_threads);
        std::vector<Vector> adjoint_vector(num_threads);
        std::vector<Matrix> sensitivity_matrix(num_threads);

        #pragma omp parallel for
        for (int i = 0; i< static_cast<int> (r_model_part.NumberOfElements()); ++i)
        {
            const unsigned int k = OpenMPUtils::ThisThread();
            auto it = r_model_part.ElementsBegin() + i;

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

        r_model_part.GetCommunicator().AssembleCurrentData(rSensitivityVariable);

        KRATOS_CATCH("");
    }

    template <typename TDataType>
    void AdjointPostprocess::UpdateConditionSensitivities(Variable<TDataType> const& rSensitivityVariable, Variable<TDataType> const& rOutputVariable)
    {
        KRATOS_TRY;

        ModelPart& r_model_part = this->GetModelPart();
        ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
        const int num_threads = OpenMPUtils::GetNumThreads();
        std::vector<Vector> sensitivity_vector(num_threads);
        std::vector<Vector> response_gradient(num_threads);
        std::vector<Vector> adjoint_vector(num_threads);
        std::vector<Matrix> sensitivity_matrix(num_threads);

        //  Assemble condition contributions.
        #pragma omp parallel for
        for (int i = 0; i< static_cast<int> (r_model_part.NumberOfConditions()); ++i)
        {
            const unsigned int k = OpenMPUtils::ThisThread();
            auto it = r_model_part.ConditionsBegin() + i;

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

        r_model_part.GetCommunicator().AssembleCurrentData(rSensitivityVariable);

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

    void AdjointPostprocess::AssembleElementSensitivityContribution(Variable<double> const& rSensitivityVariable,
                                                Vector const& rSensitivityVector,
                                                Element& rElem)
    {
        rElem.SetValue(rSensitivityVariable , rSensitivityVector[0]);
        // attention: one has to ensure that element is able to print the variable type later on his Gauss-Points
    }

    void AdjointPostprocess::AssembleElementSensitivityContribution(Variable<array_1d<double, 3>> const& rSensitivityVariable,
                                                Vector const& rSensitivityVector,
                                                Element& rElem)
    {
        rElem.SetValue(rSensitivityVariable , rSensitivityVector);
        // attention: one has to ensure that element is able to print the variable type later on his Gauss-Points
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


