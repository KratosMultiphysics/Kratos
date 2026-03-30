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

// System includes

// Base class include
#include "custom_strategies/schemes/bdf2_turbulent_scheme.h"
#include "custom_utilities/compute_div_sigma_utility.h"
#include "includes/kratos_components.h"

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
        const auto& r_divergence_stress_variable = KratosComponents<Variable<Vector>>::Get("DIVERGENCE_STRESS");

        const int nelements = static_cast<int>(rModelPart.Elements().size());
        const ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
        ModelPart::ElementsContainerType::iterator el_begin = rModelPart.ElementsBegin();

        // Find the first element with LocalSpaceDimension() != 1 to use as reference
        int first_valid_index = 0;
        for (; first_valid_index < nelements; ++first_valid_index) {
            auto it_elem = el_begin + first_valid_index;
            if (it_elem->GetGeometry().LocalSpaceDimension() != 1) break;
        }
        KRATOS_ERROR_IF(first_valid_index >= nelements)
            << "BDF2HigherOrderVMSScheme: no 2D/3D elements found to use as reference." << std::endl;

        const unsigned int gauss_point_per_knot_span =
            (el_begin + first_valid_index)->GetGeometry().size();
        const unsigned int number_of_control_points =
            (el_begin + first_valid_index)->GetGeometry().size();
        const unsigned int reference_dim =
            (el_begin + first_valid_index)->GetGeometry().WorkingSpaceDimension();
        const unsigned int reference_basis_order =
            (reference_dim == 3)
                ? static_cast<unsigned int>(std::cbrt(number_of_control_points) - 1.0)
                : static_cast<unsigned int>(std::sqrt(number_of_control_points) - 1.0);

        if (reference_basis_order <= 1) {
            return;
        }

        const unsigned int stress_size = (reference_dim == 3) ? 6 : 3;
        Matrix collected_stress(gauss_point_per_knot_span, stress_size);
        Matrix collected_shape_functions(gauss_point_per_knot_span, number_of_control_points);
        Matrix collected_shape_functions_dx(gauss_point_per_knot_span, number_of_control_points);
        Matrix collected_shape_functions_dy(gauss_point_per_knot_span, number_of_control_points);
        Matrix collected_shape_functions_dz;
        if (reference_dim == 3) {
            collected_shape_functions_dz.resize(gauss_point_per_knot_span, number_of_control_points, false);
        }
        std::vector<Element*> collected_elements;
        collected_elements.reserve(gauss_point_per_knot_span);
        unsigned int collected_count = 0;

        // Iterate over all elements, skipping any with LocalSpaceDimension()==1
        for (int k = first_valid_index; k < nelements; ++k) {
            auto it_elem = el_begin + k;
            if (it_elem->GetGeometry().LocalSpaceDimension() == 1) {
                continue; // skip 1D (inner-loop) elements wherever they appear
            }

            const unsigned int gauss_point_per_knot_span = (it_elem)->GetGeometry().size();
            const unsigned int number_of_control_points  = (it_elem)->GetGeometry().size();
            const unsigned int elem_dim = it_elem->GetGeometry().WorkingSpaceDimension();

            if (it_elem->IsActive()) {
                std::vector<Vector> stress_vector;
                it_elem->CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, stress_vector, CurrentProcessInfo);
                KRATOS_ERROR_IF(stress_vector.empty())
                    << "BDF2HigherOrderVMSScheme: element " << it_elem->Id()
                    << " returned no CAUCHY_STRESS_VECTOR values." << std::endl;
                KRATOS_ERROR_IF(stress_vector[0].size() < ((elem_dim == 3) ? 6 : 3))
                    << "BDF2HigherOrderVMSScheme: element " << it_elem->Id()
                    << " returned CAUCHY_STRESS_VECTOR of size " << stress_vector[0].size()
                    << " for dimension " << elem_dim << "." << std::endl;
                KRATOS_ERROR_IF(elem_dim != reference_dim)
                    << "BDF2HigherOrderVMSScheme: mixed element dimensions are not supported." << std::endl;
                KRATOS_ERROR_IF(gauss_point_per_knot_span != collected_stress.size1())
                    << "BDF2HigherOrderVMSScheme: inconsistent knot-span size." << std::endl;
                KRATOS_ERROR_IF(number_of_control_points != collected_shape_functions.size2())
                    << "BDF2HigherOrderVMSScheme: inconsistent number of control points." << std::endl;

                const unsigned int row_index = collected_count;
                for (std::size_t s = 0; s < stress_size; ++s) {
                    collected_stress(row_index, s) = stress_vector[0][s];
                }
                collected_elements.push_back(&(*it_elem));

                // Collect the shape function values and their derivatives at the Gauss points
                auto integration_method = it_elem->GetGeometry().GetDefaultIntegrationMethod();
                const auto& N_gausspoint = it_elem->GetGeometry().ShapeFunctionsValues(integration_method);
                auto shape_function_value = row(N_gausspoint,0);
                // Retrieve the shape function gradients (derivatives)
                const GeometryData::ShapeFunctionsGradientsType& DN_De = it_elem->GetGeometry().ShapeFunctionsLocalGradients(integration_method);
                const auto& DN_DX = DN_De[0];

                for (std::size_t cp = 0; cp < number_of_control_points; ++cp) {
                    collected_shape_functions(row_index, cp) = shape_function_value(cp);
                    collected_shape_functions_dx(row_index, cp) = DN_DX(cp, 0); // Derivative w.r.t. local x (ξ)
                    collected_shape_functions_dy(row_index, cp) = DN_DX(cp, 1); // Derivative w.r.t. local y (η)
                    if (elem_dim == 3) {
                        collected_shape_functions_dz(row_index, cp) = DN_DX(cp, 2); // Derivative w.r.t. local z (ζ)
                    }
                }

                // Increase the count of collected Gauss points
                collected_count++;

                // When the Gauss points are collected, call the utility and reset the containers
                if (collected_count == gauss_point_per_knot_span) {
                    // Call the utility passing the data of the GPs of the current knot span
                    ComputeDivSigmaUtility div_sigma_utility;
                    Matrix collected_divergence;

                    if (elem_dim == 3) {
                        collected_divergence = div_sigma_utility.ComputeDivergence(
                            collected_stress,
                            collected_shape_functions,
                            collected_shape_functions_dx,
                            collected_shape_functions_dy,
                            collected_shape_functions_dz);
                    } else {
                        collected_divergence = div_sigma_utility.ComputeDivergence(
                            collected_stress,
                            collected_shape_functions,
                            collected_shape_functions_dx,
                            collected_shape_functions_dy);
                    }

                    KRATOS_ERROR_IF(collected_divergence.size1() != collected_elements.size())
                        << "BDF2HigherOrderVMSScheme: divergence recovery size mismatch. Got "
                        << collected_divergence.size1() << " values for "
                        << collected_elements.size() << " collected elements." << std::endl;
                    KRATOS_ERROR_IF(collected_divergence.size2() != elem_dim)
                        << "BDF2HigherOrderVMSScheme: divergence recovery component mismatch. Got "
                        << collected_divergence.size2() << " values for dimension " << elem_dim << "."
                        << std::endl;

                    // setValue of DIVERGENCE(SIGMA) to the Gauss Point
                    for (size_t i = 0; i < collected_elements.size(); ++i) {
                        Vector divergence_value(elem_dim);
                        for (unsigned int d = 0; d < elem_dim; ++d) {
                            divergence_value[d] = collected_divergence(i, d);
                        }

                        collected_elements[i]->SetValue(r_divergence_stress_variable, divergence_value);
                    }

                    // Reset the collections for the next set of gauss points
                    collected_elements.clear();
                    collected_count = 0;
                }
            }
        }

        KRATOS_ERROR_IF(collected_count != 0)
            << "BDF2HigherOrderVMSScheme: collected " << collected_count
            << " active non-1D elements that do not complete a knot-span block of size "
            << gauss_point_per_knot_span << '.' << std::endl;
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
