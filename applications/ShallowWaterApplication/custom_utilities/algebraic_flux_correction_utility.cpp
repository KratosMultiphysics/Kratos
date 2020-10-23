//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

// System includes


// External includes


// Project includes
#include "algebraic_flux_correction_utility.h"
#include "shallow_water_application_variables.h"
#include "processes/calculate_nodal_area_process.h"
#include "processes/find_global_nodal_neighbours_process.h"
#include "utilities/variable_utils.h"
#include "includes/model_part.h"
#include "custom_utilities/shallow_water_utilities.h"

namespace Kratos
{
    AlgebraicFluxCorrectionUtility::AlgebraicFluxCorrectionUtility(
        ModelPart& rModelPart,
        Parameters ThisParameters) :
        mrModelPart(rModelPart)
    {
        ThisParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());
        mMaximumIterations = ThisParameters["maximum_iterations"].GetInt();
        mRebuildLevel = ThisParameters["rebuild_level"].GetInt();

        Check();
        GenerateVariablesList(ThisParameters["limiting_variables"].GetStringArray());
        InitializeNonhistoricalVariables();
        ResizeNodalAndElementalVectors();
        const auto& r_data_communicator = mrModelPart.GetCommunicator().GetDataCommunicator();
        FindGlobalNodalNeighboursProcess(r_data_communicator, mrModelPart).Execute();
        GetElementalDofList();
        AssembleElementalMassMatrices();
    }

    void AlgebraicFluxCorrectionUtility::InitializeCorrection()
    {
        if (mRebuildLevel > 0) {
            ResizeNodalAndElementalVectors();
            const auto& r_data_communicator = mrModelPart.GetCommunicator().GetDataCommunicator();
            FindGlobalNodalNeighboursProcess(r_data_communicator, mrModelPart).Execute();
            GetElementalDofList();
            AssembleElementalMassMatrices();
        }
    }

    void AlgebraicFluxCorrectionUtility::ApplyCorrection()
    {
        ComputeElementalAlgebraicFluxCorrections();
        ComputeLimiters();
        AssembleLimitedCorrections();

        size_t iteration = 1;
        while (iteration < mMaximumIterations)
        {
            ComputeRejectedAlgebraicFluxCorrections();
            GetLowOrderValues();
            ComputeLimiters();
            AssembleLimitedCorrections();
            iteration++;
        }
    }

    void AlgebraicFluxCorrectionUtility::GetHighOrderValues()
    {
        // Update the variables if need
        UpdateLimitingVariables();

        // Elemental values
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfElements()); ++i)
        {
            Vector& r_element_values_vector = *(mHighOrderValues.begin() + i);
            (mrModelPart.ElementsBegin() + i)->GetValuesVector(r_element_values_vector, 0);
        }

        // Nodal limiting values
        size_t variable_counter = 0;
        for (auto p_var : mLimitingVariables)
        {
            Vector& nodal_values_vector = mLimitingNodalHighOrderValues[variable_counter++];
            #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfNodes()); ++i)
            {
                double& r_nodal_value = *(nodal_values_vector.begin() + i);
                r_nodal_value = (mrModelPart.NodesBegin() + i)->FastGetSolutionStepValue(*p_var);
            }
        }
    }

    void AlgebraicFluxCorrectionUtility::GetLowOrderValues()
    {
        // Update the variables if need
        UpdateLimitingVariables();

        // Elemental values
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfElements()); ++i)
        {
            Vector& r_element_values_vector = *(mLowOrderValues.begin() + i);
            (mrModelPart.ElementsBegin() + i)->GetValuesVector(r_element_values_vector, 0);
        }

        // Nodal values
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfNodes()); ++i)
        {
            Vector& lo_values = *(mDofsLowOrderValues.begin() + i);
            const auto& dofs = (mrModelPart.NodesBegin() + i)->GetDofs();
            lo_values.resize(dofs.size(), false);
            for (size_t d = 0; d < dofs.size(); ++d)
            {
                lo_values[d] = (*dofs[d])(0);
            }
        }

        // Nodal limiting values
        size_t variable_counter = 0;
        for (auto p_var : mLimitingVariables)
        {
            Vector& nodal_values_vector = mLimitingNodalLowOrderValues[variable_counter++];
            #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfNodes()); ++i)
            {
                double& r_nodal_value = *(nodal_values_vector.begin() + i);
                r_nodal_value = (mrModelPart.NodesBegin() + i)->FastGetSolutionStepValue(*p_var);
            }
        }
    }

    void AlgebraicFluxCorrectionUtility::UpdateLimitingVariables()
    {
        for (auto p_var : mLimitingVariables)
        {
            if (*p_var == INTERNAL_ENERGY)
            {
                ShallowWaterUtilities().ComputeEnergy(mrModelPart);
            }
            else if (*p_var == FREE_SURFACE_ELEVATION)
            {
                ShallowWaterUtilities().ComputeFreeSurfaceElevation(mrModelPart);
            }
        }
    }

    void AlgebraicFluxCorrectionUtility::ComputeElementalAlgebraicFluxCorrections()
    {
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfElements()); ++i)
        {
            auto& corrections = *(mAlgebraicFluxCorrections.begin() + i);
            auto ho_values = *(mHighOrderValues.begin() + i);
            auto lo_values = *(mLowOrderValues.begin() + i);
            auto mass_matrix = *(mElementalMassMatrices.begin() + i);
            corrections.resize(ho_values.size(), false);
            for (size_t j = 0; j < ho_values.size(); ++j)
            {
                corrections[j] = (ho_values[j] - lo_values[j]) * mass_matrix[j];
            }
        }
    }

    void AlgebraicFluxCorrectionUtility::AssembleElementalMassMatrices()
    {
        const auto& r_process_info = mrModelPart.GetProcessInfo();
        const size_t domain_size = r_process_info.GetValue(DOMAIN_SIZE);
        CalculateNodalAreaProcess<CalculateNodalAreaSettings::SaveAsNonHistoricalVariable>(mrModelPart, domain_size).Execute();

        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfElements()); ++i)
        {
            auto it_elem = mrModelPart.ElementsBegin() + i;
            Matrix mass_matrix;
            it_elem->CalculateMassMatrix(mass_matrix, r_process_info);
            const size_t local_size = mass_matrix.size1();
            const size_t number_of_nodes = it_elem->GetGeometry().PointsNumber();
            const size_t block_size = local_size / number_of_nodes;
            const auto& r_geom = it_elem->GetGeometry();
            Vector lumped_mass_matrix = ZeroVector(local_size);
            for (size_t i_node = 0; i_node < number_of_nodes; ++i_node)
            {
                for (size_t i_dof = 0; i_dof < block_size; ++ i_dof)
                {
                    const size_t row = i_node*block_size + i_dof;
                    for (size_t col = 0; col < local_size; ++col)
                    {
                        lumped_mass_matrix[row] += mass_matrix(row, col);
                    }
                    lumped_mass_matrix[row] /= r_geom[i_node].GetValue(NODAL_AREA);
                }
            }
            *(mElementalMassMatrices.begin() + i) = lumped_mass_matrix;
        }
    }

    void AlgebraicFluxCorrectionUtility::ComputeLimiters()
    {
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfElements()); ++i)
        {
            *(mElementalLimiters.begin() + i) = 1.0;
        }

        size_t variable_counter = 0;
        for (auto p_var : mLimitingVariables)
        {
            const Vector& high_order_values = mLimitingNodalHighOrderValues[variable_counter];
            const Vector& low_order_values = mLimitingNodalLowOrderValues[variable_counter];
            ComputeLimiterSingleVariable(*p_var, high_order_values, low_order_values);
            ++variable_counter;
        }
    }

    void AlgebraicFluxCorrectionUtility::ComputeLimiterSingleVariable(
        const Variable<double>& rVariable,
        const Vector& rHighOrderValues,
        const Vector& rLowOrderValues)
    {
        // first step: Get the nodal positive and negative contributions
        const Vector nodal_contributions = rHighOrderValues - rLowOrderValues;
        const size_t number_of_nodes = mrModelPart.NumberOfNodes();
        Vector nodal_positive_contributions(number_of_nodes);
        Vector nodal_negative_contributions(number_of_nodes);
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(number_of_nodes); ++i)
        {
            auto it_node = mrModelPart.NodesBegin() + i;
            const double aec = *(nodal_contributions.begin() + i);
            it_node->GetValue(ALGEBRAIC_CONTRIBUTION) = aec;
            *(nodal_positive_contributions.begin() + i) = std::max(0.0, aec);
            *(nodal_negative_contributions.begin() + i) = std::min(0.0, aec);

            const double u_n = it_node->FastGetSolutionStepValue(rVariable,1);
            const double u_l = *(rLowOrderValues.begin() + i);
            const double u_max = std::max(u_l, u_n);
            const double u_min = std::min(u_l, u_n);
            it_node->GetValue(MAXIMUM_VALUE) = u_max;
            it_node->GetValue(MINIMUM_VALUE) = u_min;
        }

        // second step: get the maximum and minimum increments.
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(number_of_nodes); ++i)
        {
            auto it_node = mrModelPart.NodesBegin() + i;
            auto neighbor_nodes = it_node->GetValue(NEIGHBOUR_NODES);
            double u_max = it_node->GetValue(MAXIMUM_VALUE);
            double u_min = it_node->GetValue(MINIMUM_VALUE);
            for (size_t j = 0; j < neighbor_nodes.size(); ++j)
            {
                u_max = std::max(u_max, neighbor_nodes[j].GetValue(MAXIMUM_VALUE));
                u_min = std::min(u_min, neighbor_nodes[j].GetValue(MINIMUM_VALUE));
            }
            const double u_l = *(rLowOrderValues.begin() + i);
            const double nodal_max_increment = u_max - u_l;
            const double nodal_min_increment = u_min - u_l;

            // using the maximum/minimum increments to compute the positive/negative ratios
            const double positive_contribution = *(nodal_positive_contributions.begin() + i);
            if (positive_contribution > mEpsilon) {
                it_node->GetValue(POSITIVE_RATIO) = std::min(1.0, nodal_max_increment / positive_contribution);
            } else {
                it_node->GetValue(POSITIVE_RATIO) = 0.0;
            }
            const double negative_contribution = *(nodal_negative_contributions.begin() + i);
            if (negative_contribution < -mEpsilon) {
                it_node->GetValue(NEGATIVE_RATIO) = std::min(1.0, nodal_min_increment / negative_contribution);
            } else {
                it_node->GetValue(NEGATIVE_RATIO) = 0.0;
            }
        }

        // third step: get the elemental limiter
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfElements()); ++i)
        {
            auto it_elem = mrModelPart.ElementsBegin() + i;
            auto& geom = it_elem->GetGeometry();
            double& c = *(mElementalLimiters.begin() + i);
            for (auto& r_node : geom)
            {
                if (r_node.GetValue(ALGEBRAIC_CONTRIBUTION) > 0.0){
                    c = std::min(c, r_node.GetValue(POSITIVE_RATIO));
                } else {
                    c = std::min(c, r_node.GetValue(NEGATIVE_RATIO));
                }
            }
        }
    }

    void AlgebraicFluxCorrectionUtility::AssembleLimitedCorrections()
    {
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfNodes()); ++i)
        {
            const Vector& lo_values = *(mDofsLowOrderValues.begin() + i);
            auto& dofs = (mrModelPart.NodesBegin() + i)->GetDofs();
            for (size_t d = 0; d < dofs.size(); ++d)
            {
                (*dofs[d])(0) = lo_values[d];
            }
        }

        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfElements()); ++i)
        {
            const Vector& corrections = *(mAlgebraicFluxCorrections.begin() + i);
            double limiter = *(mElementalLimiters.begin() + i);
            DofsVectorType dofs = *(mElementalDofs.begin() + i);
            for (size_t d = 0; d < dofs.size(); ++d)
            {
                #pragma omp atomic
                (*dofs[d])(0) += limiter * corrections[d];
            }
        }
    }

    void AlgebraicFluxCorrectionUtility::ComputeRejectedAlgebraicFluxCorrections()
    {
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfElements()); ++i)
        {
            Vector& corrections = *(mAlgebraicFluxCorrections.begin() + i);
            double limiter = *(mElementalLimiters.begin() + i);
            corrections *= (1 - limiter);
        }
    }

    void AlgebraicFluxCorrectionUtility::InitializeNonhistoricalVariables()
    {
        VariableUtils().SetNonHistoricalVariableToZero(MAXIMUM_VALUE, mrModelPart.Nodes());
        VariableUtils().SetNonHistoricalVariableToZero(MINIMUM_VALUE, mrModelPart.Nodes());
        VariableUtils().SetNonHistoricalVariableToZero(POSITIVE_RATIO, mrModelPart.Nodes());
        VariableUtils().SetNonHistoricalVariableToZero(NEGATIVE_RATIO, mrModelPart.Nodes());
        VariableUtils().SetNonHistoricalVariableToZero(ALGEBRAIC_CONTRIBUTION, mrModelPart.Nodes());
        VariableUtils().SetNonHistoricalVariableToZero(MAXIMUM_VALUE, mrModelPart.Elements());
        VariableUtils().SetNonHistoricalVariableToZero(MINIMUM_VALUE, mrModelPart.Elements());
    }

    void AlgebraicFluxCorrectionUtility::ResizeNodalAndElementalVectors()
    {
        size_t number_of_elements = mrModelPart.NumberOfElements();
        size_t number_of_nodes = mrModelPart.NumberOfNodes();
        size_t number_of_variables = mLimitingVariables.size();
        mHighOrderValues.resize(number_of_elements);
        mLowOrderValues.resize(number_of_elements);
        mAlgebraicFluxCorrections.resize(number_of_elements);
        mElementalMassMatrices.resize(number_of_elements);
        mElementalLimiters.resize(number_of_elements);
        mElementalDofs.resize(number_of_elements);
        mDofsLowOrderValues.resize(number_of_nodes);
        mLimitingNodalHighOrderValues.resize(number_of_variables);
        mLimitingNodalLowOrderValues.resize(number_of_variables);
        for (size_t i = 0; i < number_of_variables; ++i)
        {
            mLimitingNodalHighOrderValues[i].resize(number_of_nodes);
            mLimitingNodalLowOrderValues[i].resize(number_of_nodes);
        }
    }

    void AlgebraicFluxCorrectionUtility::GenerateVariablesList(const std::vector<std::string>& rVariablesNames)
    {
        size_t n_variables = rVariablesNames.size();
        mLimitingVariables.reserve(n_variables);

        for (size_t v = 0; v < n_variables; ++v)
        {
            const std::string variable_name = rVariablesNames[v];

            if(KratosComponents<Variable<double>>::Has(variable_name)){
                const auto& r_var = KratosComponents<Variable<double>>::Get(variable_name);
                mLimitingVariables.push_back(&r_var);
            } else {
                KRATOS_ERROR << "AlgebraicFluxCorrectionUtility. Only double variables are allowed in the limiting variables list." << std::endl;
            }
        }
    }

    const Parameters AlgebraicFluxCorrectionUtility::GetDefaultParameters() const
    {
        const Parameters default_parameters = Parameters(R"(
        {
            "name"               : "algebraic_flux_correction_utility",
            "limiting_variables" : ["VARIABLE_NAME"],
            "maximum_iterations" : 1,
            "rebuild_level"      : 0
        })");
        return default_parameters;
    }

    int AlgebraicFluxCorrectionUtility::Check()
    {
        const auto domain_size = mrModelPart.GetProcessInfo().GetValue(DOMAIN_SIZE);
        KRATOS_ERROR_IF(((domain_size != 2) && (domain_size != 3))) << "Domain size must be 2 or 3" << std::endl;
        return 1;
    }

    void AlgebraicFluxCorrectionUtility::GetElementalDofList()
    {
        const ProcessInfo& r_process_info = mrModelPart.GetProcessInfo();
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfElements()); ++i)
        {
            (mrModelPart.ElementsBegin() + i)->GetDofList(*(mElementalDofs.begin() + i), r_process_info);
        }
    }
}
