//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolò Antonelli
//

#pragma once

// Base class include
#include "custom_strategies/schemes/bdf2_turbulent_scheme.h"
#include "../applications/LinearSolversApplication/custom_utilities/compute_div_sigma_utility.h"

namespace Kratos {

template<class TSparseSpace, class TDenseSpace>
class BDF2HigherOrderVMSScheme : public BDF2TurbulentScheme<TSparseSpace, TDenseSpace>
{
public:
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(BDF2HigherOrderVMSScheme);

    /// Base type
    using BaseType = BDF2TurbulentScheme<TSparseSpace, TDenseSpace>;
    using TSystemMatrixType = typename TSparseSpace::MatrixType;
    using TSystemVectorType = typename TSparseSpace::VectorType;
    using LocalSystemMatrixType = typename BaseType::LocalSystemMatrixType;
    using LocalSystemVectorType = typename BaseType::LocalSystemVectorType;

    /// Constructor: same as parent
    BDF2HigherOrderVMSScheme()
    : BaseType() {}

    BDF2HigherOrderVMSScheme(Process::Pointer pTurbulenceModel)
    : BaseType(pTurbulenceModel) {}

    BDF2HigherOrderVMSScheme(const Kratos::Variable<int>& rPeriodicVar)
    : BaseType(rPeriodicVar) {}

    /// Destructor
    ~BDF2HigherOrderVMSScheme() override = default;

    /// Info
    std::string Info() const override {
        return "BDF2HigherOrderVMSScheme";
    }

    // Override condition assembly
    void CalculateSystemContributions(
        Condition& rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY
        rCurrentCondition.EquationIdVector(rEquationId, rCurrentProcessInfo);
        rCurrentCondition.CalculateLocalSystem(LHS_Contribution, RHS_Contribution, rCurrentProcessInfo);
        KRATOS_CATCH("")
    }

    void CalculateRHSContribution(
        Condition& rCurrentCondition,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY
        rCurrentCondition.EquationIdVector(rEquationId, rCurrentProcessInfo);
        rCurrentCondition.CalculateRightHandSide(RHS_Contribution, rCurrentProcessInfo);
        KRATOS_CATCH("")
    }

    /// Set the time iteration coefficients
    void InitializeNonLinIteration(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        // Call the base class implementation first
        BDF2TurbulentScheme<TSparseSpace, TDenseSpace>::InitializeNonLinIteration(rModelPart, A, Dx, b);

        const int nelements = static_cast<int>(rModelPart.Elements().size());
        const ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
        ModelPart::ElementsContainerType::iterator el_begin = rModelPart.ElementsBegin();
        ModelPart::ConditionsContainerType::iterator cond_begin = rModelPart.ConditionsBegin();

        // Find the first element with LocalSpaceDimension() != 1 to use as reference
        int first_valid_index = 0;
        for (; first_valid_index < nelements; ++first_valid_index) {
            auto it_elem = el_begin + first_valid_index;
            if (it_elem->GetGeometry().LocalSpaceDimension() != 1) break;
        }
        KRATOS_ERROR_IF(first_valid_index >= nelements)
            << "BDF2HigherOrderVMSScheme: no 2D/3D elements found to use as reference." << std::endl;

        // KRATOS_INFO("BDF2HigherOrderVMSScheme") << "Number of Gauss points per knot span: " << gauss_point_per_knot_span << std::endl;
        // KRATOS_INFO("BDF2HigherOrderVMSScheme") << "Number of control points: " << number_of_control_points << std::endl;

        std::vector<std::vector<double>> collected_stress;
        std::vector<std::vector<double>> collected_shape_functions;
        std::vector<std::vector<double>> collected_shape_functions_dx;
        std::vector<std::vector<double>> collected_shape_functions_dy;
        std::vector<array_1d<double, 3>> collected_coordinates;
        unsigned int collected_count = 0;

        // Iterate over all elements, skipping any with LocalSpaceDimension()==1
        for (int k = first_valid_index; k < nelements; ++k) {
            auto it_elem = el_begin + k;
            if (it_elem->GetGeometry().LocalSpaceDimension() == 1) {
                continue; // skip 1D (inner-loop) elements wherever they appear
            }

            const unsigned int gauss_point_per_knot_span = (it_elem)->GetGeometry().size();
            const unsigned int number_of_control_points  = (it_elem)->GetGeometry().size();

            if (it_elem->IsActive()) {
                std::vector<Vector> stress_vector;
                std::vector<Matrix> D_constitutive_matrix;
                it_elem->CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, stress_vector, CurrentProcessInfo);

                // Collect
                std::vector<double> stress_vector_temp(3);
                stress_vector_temp[0] = stress_vector[0][0];
                stress_vector_temp[1] = stress_vector[0][1];
                stress_vector_temp[2] = stress_vector[0][2];
                
                collected_stress.push_back(stress_vector_temp);
                array_1d<double, 3> center_coord = it_elem->GetGeometry().Center();
                collected_coordinates.push_back(center_coord);

                // Collect the shape function values and their derivatives at the Gauss points
                auto integration_method = it_elem->GetGeometry().GetDefaultIntegrationMethod();
                const auto& N_gausspoint = it_elem->GetGeometry().ShapeFunctionsValues(integration_method);
                auto shape_function_value = row(N_gausspoint,0);
                // Retrieve the shape function gradients (derivatives)
                const GeometryData::ShapeFunctionsGradientsType& DN_De = it_elem->GetGeometry().ShapeFunctionsLocalGradients(integration_method);
                // Assuming gauss_point_per_knot_span = 9, collect 9 basis function values
                std::vector<double> element_shape_functions;
                std::vector<double> element_shape_functions_dx;
                std::vector<double> element_shape_functions_dy;

                for (std::size_t cp = 0; cp < number_of_control_points; ++cp) {
                    // Collect the shape function value for the current Gauss point (assuming 1 value per Gauss point)
                    double shape_function_value_cp = shape_function_value(cp);
                    element_shape_functions.push_back(shape_function_value_cp);
                    auto DN_DX = DN_De[0];
                    element_shape_functions_dx.push_back(DN_DX(cp, 0)); // Derivative w.r.t. local x (ξ)
                    element_shape_functions_dy.push_back(DN_DX(cp, 1)); // Derivative w.r.t. local y (η)
                }
                // Store the collected shape function values and their derivatives
                collected_shape_functions.push_back(element_shape_functions);
                collected_shape_functions_dx.push_back(element_shape_functions_dx);
                collected_shape_functions_dy.push_back(element_shape_functions_dy);        

                // Increase the count of collected Gauss points
                collected_count++;

                // When the Gauss points are collected, call the utility and reset the containers
                if (collected_count == gauss_point_per_knot_span) {
                    // Call the utility passing the data of the GPs of the current knot span
                    ComputeDivSigmaUtility div_sigma_utility;
                    
                    div_sigma_utility.SetInputData(collected_stress, collected_coordinates, collected_shape_functions, 
                        collected_shape_functions_dx, collected_shape_functions_dy);

                    std::vector<std::vector<double>> collected_divergence = div_sigma_utility.ComputeDivergence() ;

                    // setValue of DIVERGENCE(SIGMA) to the Gauss Point
                    for (size_t i = 0; i < gauss_point_per_knot_span; ++i) {
                        Vector divergence_value(2);
                        divergence_value[0] = collected_divergence[i][0]; // div_sigma_1
                        divergence_value[1] = collected_divergence[i][1]; // div_sigma_2
                        
                        auto it_element_in_knot_span = it_elem - (gauss_point_per_knot_span-1) + i;
                        // Apply the divergence values back to the Gauss point
                        it_element_in_knot_span->SetValue(RECOVERED_STRESS, divergence_value);
                    }                    

                    // Reset the collections for the next set of gauss points
                    collected_stress.clear();
                    collected_coordinates.clear();
                    collected_shape_functions.clear();
                    collected_shape_functions_dx.clear();
                    collected_shape_functions_dy.clear();
                    collected_count = 0;
                }
            }
        }
    }
};

/// Input/output operators
template<class TSparseSpace, class TDenseSpace>
inline std::istream& operator>>(std::istream& rIStream, BDF2HigherOrderVMSScheme<TSparseSpace, TDenseSpace>& rThis) {
    return rIStream;
}

template<class TSparseSpace, class TDenseSpace>
inline std::ostream& operator<<(std::ostream& rOStream, const BDF2HigherOrderVMSScheme<TSparseSpace, TDenseSpace>& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos