//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Reza Najian Asl
//                   Suneth Warnakulasuriya
//

// System includes
#include <sstream>
#include <type_traits>
#include <functional>

// Project includes
#include "expression/variable_expression_io.h"
#include "includes/define.h"
#include "includes/variables.h"
#include "utilities/element_size_calculator.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/variable_utils.h"

// Application includes
#include "custom_utilities/optimization_utils.h"
#include "custom_utilities/properties_variable_expression_io.h"
#include "optimization_application_variables.h"

// Include base h
#include "mass_response_utils.h"

namespace Kratos
{

///@name Kratos Classes
///@{

bool MassResponseUtils::HasVariableInProperties(
    const ModelPart& rModelPart,
    const Variable<double>& rVariable)
{
    KRATOS_TRY

    bool local_has_variable = false;
    if (rModelPart.NumberOfElements() > 0) {
        const auto& r_properties = rModelPart.Elements().front().GetProperties();
        local_has_variable = r_properties.Has(rVariable);
    }
    return rModelPart.GetCommunicator().GetDataCommunicator().OrReduceAll(local_has_variable);

    KRATOS_CATCH("");
}

void MassResponseUtils::Check(const ModelPart& rModelPart)
{
    const auto& r_data_communicator = rModelPart.GetCommunicator().GetDataCommunicator();

    if (!OptimizationUtils::IsVariableExistsInAllContainerProperties(rModelPart.Elements(), DENSITY, r_data_communicator)) {
        KRATOS_ERROR << "Some elements' properties in " << rModelPart.FullName()
                     << " does not have DENSITY variable.";
    }

    if (OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties(rModelPart.Elements(), THICKNESS, r_data_communicator) &&
        OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties(rModelPart.Elements(), CROSS_AREA, r_data_communicator)) {
        KRATOS_ERROR << rModelPart.FullName() << " has elements consisting THICKNESS and CROSS_AREA. "
                     << "Please break down this response to SumResponseFunction where each sub "
                     << "response function only has elements with either THICKNESS or CROSS_AREA.";
    }

    if (OptimizationUtils::GetContainerEntityGeometryType(rModelPart.Elements(), r_data_communicator) == GeometryData::KratosGeometryType::Kratos_generic_type) {
        KRATOS_ERROR << rModelPart.FullName() << " has elements with different geometry types. "
                     << "Please break down this response to SumResponseFunction where each sub "
                     << "response function only has elements with one geometry type.";
    }
}

double MassResponseUtils::CalculateValue(const ModelPart& rModelPart)
{
    KRATOS_TRY

    if (rModelPart.GetCommunicator().GlobalNumberOfElements() == 0) {
        return 0.0;
    }

    KRATOS_ERROR_IF_NOT(HasVariableInProperties(rModelPart, DENSITY))
        << "DENSITY is not found in element properties of " << rModelPart.FullName() << ".\n";

    KRATOS_ERROR_IF(HasVariableInProperties(rModelPart, THICKNESS) && HasVariableInProperties(rModelPart, CROSS_AREA))
        << rModelPart.FullName()
        << " has elements with properties having both THICKNESS and CROSS_AREA. Please separate the model part such that either one of them is present in elemental properties.\n";

    std::function<double(const ModelPart::ElementType&)> get_thickness;
    if (HasVariableInProperties(rModelPart, THICKNESS)) {
        get_thickness = [](const ModelPart::ElementType& rElement) -> double { return rElement.GetProperties()[THICKNESS]; };
    } else {
        get_thickness = [](const ModelPart::ElementType& rElement) -> double { return 1.0; };
    }

    std::function<double(const ModelPart::ElementType&)> get_cross_area;
    if (HasVariableInProperties(rModelPart, CROSS_AREA)) {
        get_cross_area = [](const ModelPart::ElementType& rElement) -> double { return rElement.GetProperties()[CROSS_AREA]; };
    } else {
        get_cross_area = [](const ModelPart::ElementType& rElement) -> double { return 1.0; };
    }

    const double local_mass = block_for_each<SumReduction<double>>(rModelPart.Elements(), [&](const auto& rElement) {
        return rElement.GetGeometry().DomainSize() * rElement.GetProperties()[DENSITY] * get_thickness(rElement) * get_cross_area(rElement);
    });

    return rModelPart.GetCommunicator().GetDataCommunicator().SumAll(local_mass);

    KRATOS_CATCH("")
}

void MassResponseUtils::CalculateGradient(
    const PhysicalFieldVariableTypes& rPhysicalVariable,
    ModelPart& rGradientRequiredModelPart,
    ModelPart& rGradientComputedModelPart,
    std::vector<ContainerExpressionType>& rListOfContainerExpressions,
    const double PerturbationSize)
{
    KRATOS_TRY

    std::visit([&](auto pVariable) {
        if (*pVariable == DENSITY) {
            // clears the existing values
            block_for_each(rGradientRequiredModelPart.Elements(), [](auto& rElement) { rElement.GetProperties().SetValue(DENSITY_SENSITIVITY, 0.0); });

            // computes density sensitivty and store it within each elements' properties
            CalculateMassDensityGradient(rGradientComputedModelPart, DENSITY_SENSITIVITY);
        } else if (*pVariable == THICKNESS) {
            // clears the existing values
            block_for_each(rGradientRequiredModelPart.Elements(), [](auto& rElement) { rElement.GetProperties().SetValue(THICKNESS_SENSITIVITY, 0.0); });

            // computes density sensitivty and store it within each elements' properties
            CalculateMassThicknessGradient(rGradientComputedModelPart, THICKNESS_SENSITIVITY);
        } else if (*pVariable == CROSS_AREA) {
            // clears the existing values
            block_for_each(rGradientRequiredModelPart.Elements(), [](auto& rElement) { rElement.GetProperties().SetValue(CROSS_AREA_SENSITIVITY, 0.0); });

            // computes density sensitivty and store it within each elements' properties
            CalculateMassCrossAreaGradient(rGradientComputedModelPart, CROSS_AREA_SENSITIVITY);
        } else if (*pVariable == SHAPE) {
            // clears the existing values
            VariableUtils().SetNonHistoricalVariableToZero(SHAPE_SENSITIVITY, rGradientRequiredModelPart.Nodes());

            // computes density sensitivty and store it within each elements' properties
            CalculateMassShapeGradient(rGradientComputedModelPart, SHAPE_SENSITIVITY, PerturbationSize);
        } else {
            KRATOS_ERROR
                << "Unsupported sensitivity w.r.t. " << pVariable->Name()
                << " requested. Followings are supported sensitivity variables:"
                << "\n\t" << DENSITY.Name()
                << "\n\t" << THICKNESS.Name()
                << "\n\t" << CROSS_AREA.Name()
                << "\n\t" << SHAPE.Name();
        }

        // now fill the container expressions
        for (auto& p_container_expression : rListOfContainerExpressions) {
            std::visit([pVariable](auto& pContainerExpression){
                using container_type = std::decay_t<decltype(*pContainerExpression)>;

                if (*pVariable == SHAPE) {
                    if constexpr(std::is_same_v<container_type, ContainerExpression<ModelPart::NodesContainerType>>) {
                        VariableExpressionIO::Read(*pContainerExpression, &SHAPE_SENSITIVITY, false);
                    } else {
                        KRATOS_ERROR << "Requesting sensitivity w.r.t. "
                                        "SHAPE for a Container expression "
                                        "which is not a NodalExpression. [ "
                                        "Requested container expression = "
                                        << *pContainerExpression << " ].\n";
                    }
                } else {
                    if constexpr(std::is_same_v<container_type, ContainerExpression<ModelPart::ElementsContainerType>>) {
                        const auto& sensitivity_variable = KratosComponents<Variable<double>>::Get(pVariable->Name() + "_SENSITIVITY");
                        PropertiesVariableExpressionIO::Read(*pContainerExpression, &sensitivity_variable);
                    } else {
                        KRATOS_ERROR << "Requesting sensitivity w.r.t. "
                                     << pVariable->Name()
                                     << " for a Container expression "
                                        "which is not an ElementExpression. [ "
                                        "Requested container expression = "
                                     << *pContainerExpression << " ].\n";
                    }
                }


            }, p_container_expression);
        }
    }, rPhysicalVariable);

    KRATOS_CATCH("");
}

void MassResponseUtils::CalculateMassShapeGradient(
    ModelPart& rModelPart,
    const Variable<array_1d<double, 3>>& rOutputGradientVariable,
    const double PerturbationSize)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(rModelPart.NumberOfElements() == 0) << rModelPart.FullName() << " does not contain any elements.\n";

    KRATOS_ERROR_IF_NOT(HasVariableInProperties(rModelPart, DENSITY))
        << "DENSITY is not found in element properties of " << rModelPart.FullName() << ".\n";

    KRATOS_ERROR_IF(HasVariableInProperties(rModelPart, THICKNESS) && HasVariableInProperties(rModelPart, CROSS_AREA))
        << rModelPart.FullName()
        << " has elements with properties having both THICKNESS and CROSS_AREA. Please separate the model part such that either one of them is present in elemental properties.\n";

    std::function<double(const ModelPart::ElementType&)> get_thickness;
    if (HasVariableInProperties(rModelPart, THICKNESS)) {
        get_thickness = [](const ModelPart::ElementType& rElement) -> double { return rElement.GetProperties()[THICKNESS]; };
    } else {
        get_thickness = [](const ModelPart::ElementType& rElement) -> double { return 1.0; };
    }

    std::function<double(const ModelPart::ElementType&)> get_cross_area;
    if (HasVariableInProperties(rModelPart, CROSS_AREA)) {
        get_cross_area = [](const ModelPart::ElementType& rElement) -> double { return rElement.GetProperties()[CROSS_AREA]; };
    } else {
        get_cross_area = [](const ModelPart::ElementType& rElement) -> double { return 1.0; };
    }

    using VolumeDerivativeMethodType = std::function<double(IndexType, IndexType, GeometryType&)>;

    VolumeDerivativeMethodType volume_derivative_method;
    bool is_analytical_derivatives_used = true;
    switch (rModelPart.Elements().begin()->GetGeometry().GetGeometryType()) {
        case GeometryData::KratosGeometryType::Kratos_Line2D2:
            volume_derivative_method = [](IndexType NodeIndex, IndexType DirectionIndex,  GeometryType& rGeometry) {
                const double lx = rGeometry[0].X() - rGeometry[1].X();
                const double lx_derivative = ((NodeIndex == 0) - (NodeIndex == 1)) * (DirectionIndex == 0);

                const double ly = rGeometry[0].Y() - rGeometry[1].Y();
                const double ly_derivative = ((NodeIndex == 0) - (NodeIndex == 1)) * (DirectionIndex == 1);

                const double length = lx * lx + ly * ly;
                const double length_derivative = 2 * lx * lx_derivative + 2 * ly * ly_derivative;

                return 0.5 * length_derivative / std::sqrt( length );
            };
            break;
        case GeometryData::KratosGeometryType::Kratos_Line3D2:
            volume_derivative_method = [](IndexType NodeIndex, IndexType DirectionIndex,  GeometryType& rGeometry) {
                const double lx = rGeometry[0].X() - rGeometry[1].X();
                const double lx_derivative = ((NodeIndex == 0) - (NodeIndex == 1)) * (DirectionIndex == 0);

                const double ly = rGeometry[0].Y() - rGeometry[1].Y();
                const double ly_derivative = ((NodeIndex == 0) - (NodeIndex == 1)) * (DirectionIndex == 1);

                const double lz = rGeometry[0].Z() - rGeometry[1].Z();
                const double lz_derivative = ((NodeIndex == 0) - (NodeIndex == 1)) * (DirectionIndex == 2);

                const double length = lx * lx + ly * ly + lz * lz;
                const double length_derivative = 2 * lx * lx_derivative + 2 * ly * ly_derivative + 2 * lz * lz_derivative;

                return 0.5 * length_derivative / std::sqrt( length );
            };
            break;
        case GeometryData::KratosGeometryType::Kratos_Triangle2D3:
            volume_derivative_method = [](IndexType NodeIndex, IndexType DirectionIndex,  GeometryType& rGeometry) {
                return 2.0 * ElementSizeCalculator<2, 3>::AverageElementSize(rGeometry) * ElementSizeCalculator<2, 3>::AverageElementSizeDerivative(NodeIndex, DirectionIndex, rGeometry);
            };
            break;
        case GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4:
            volume_derivative_method = [](IndexType NodeIndex, IndexType DirectionIndex,  GeometryType& rGeometry) {
                return 2.0 * ElementSizeCalculator<2, 4>::AverageElementSize(rGeometry) * ElementSizeCalculator<2, 4>::AverageElementSizeDerivative(NodeIndex, DirectionIndex, rGeometry);
            };
            break;
        case GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4:
            volume_derivative_method = [](IndexType NodeIndex, IndexType DirectionIndex,  GeometryType& rGeometry) {
                return 3.0 * std::pow(ElementSizeCalculator<3, 4>::AverageElementSize(rGeometry), 2) * ElementSizeCalculator<3, 4>::AverageElementSizeDerivative(NodeIndex, DirectionIndex, rGeometry);
            };
            break;
        case GeometryData::KratosGeometryType::Kratos_Prism3D6:
            volume_derivative_method = [](IndexType NodeIndex, IndexType DirectionIndex,  GeometryType& rGeometry) {
                return 3.0 * std::pow(ElementSizeCalculator<3, 6>::AverageElementSize(rGeometry), 2) * ElementSizeCalculator<3, 6>::AverageElementSizeDerivative(NodeIndex, DirectionIndex, rGeometry);
            };
            break;
        default:
            is_analytical_derivatives_used = false;
            volume_derivative_method = [PerturbationSize](IndexType NodeIndex, IndexType DirectionIndex, GeometryType& rGeometry) {
                auto& coordinates = rGeometry[NodeIndex].Coordinates();
                coordinates[DirectionIndex] += PerturbationSize;
                const double perturbed_domain_size = rGeometry.DomainSize();
                coordinates[DirectionIndex] -= PerturbationSize;
                return perturbed_domain_size;
            };
            break;
    }

    if (is_analytical_derivatives_used) {
        block_for_each(rModelPart.Elements(), [&](auto& rElement){
            auto& r_geometry = rElement.GetGeometry();
            const IndexType dimension = r_geometry.WorkingSpaceDimension();

            const double density = rElement.GetProperties()[DENSITY];
            const double thickness = get_thickness(rElement);
            const double cross_area = get_cross_area(rElement);

            for (IndexType c = 0; c < r_geometry.PointsNumber(); ++c) {
                auto& r_derivative_value = r_geometry[c].GetValue(rOutputGradientVariable);

                for (IndexType k = 0; k < dimension; ++k) {
                    const double derivative_value = volume_derivative_method(c, k, r_geometry) * thickness * density * cross_area;
                    AtomicAdd(r_derivative_value[k], derivative_value);
                }
            }
        });
    } else {
        block_for_each(rModelPart.Elements(), Node::Pointer(), [&](auto& rElement, auto& pThreadLocalNode){
            if (!pThreadLocalNode) {
                pThreadLocalNode = Kratos::make_intrusive<Node>(1, 0, 0, 0);
            }

            auto& r_geometry = rElement.GetGeometry();
            const IndexType dimension = r_geometry.WorkingSpaceDimension();

            const double density = rElement.GetProperties()[DENSITY];
            const double thickness = get_thickness(rElement);
            const double cross_area = get_cross_area(rElement);
            const double initial_domain_size = r_geometry.DomainSize();

            for (IndexType c = 0; c < r_geometry.PointsNumber(); ++c) {
                auto& r_derivative_value = r_geometry[c].GetValue(rOutputGradientVariable);

                // get the geometry node pointer
                auto& p_node = r_geometry(c);

                // now copy the node data to thread local node using the operator= in Node
                (*pThreadLocalNode) = (*p_node);

                // now swap entity node with the thread local node
                std::swap(p_node, pThreadLocalNode);

                for (IndexType k = 0; k < dimension; ++k) {
                    const double perturbed_domain_size = volume_derivative_method(c, k, r_geometry);
                    const double domain_size_derivative = (perturbed_domain_size - initial_domain_size) * thickness * density * cross_area / PerturbationSize;
                    r_derivative_value[k] += domain_size_derivative;
                }

                // revert back the node change.
                std::swap(p_node, pThreadLocalNode);
            }
        });
    }

    rModelPart.GetCommunicator().AssembleNonHistoricalData(rOutputGradientVariable);

    KRATOS_CATCH("")
}

void MassResponseUtils::CalculateMassDensityGradient(
    ModelPart& rModelPart,
    const Variable<double>& rOutputGradientVariable)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(HasVariableInProperties(rModelPart, DENSITY))
        << "DENSITY is not found in element properties of " << rModelPart.FullName() << ".\n";

    KRATOS_ERROR_IF(HasVariableInProperties(rModelPart, THICKNESS) && HasVariableInProperties(rModelPart, CROSS_AREA))
        << rModelPart.FullName()
        << " has elements with properties having both THICKNESS and CROSS_AREA. Please separate the model part such that either one of them is present in elemental properties.\n";

    std::function<double(const ModelPart::ElementType&)> get_thickness;
    if (HasVariableInProperties(rModelPart, THICKNESS)) {
        get_thickness = [](const ModelPart::ElementType& rElement) -> double { return rElement.GetProperties()[THICKNESS]; };
    } else {
        get_thickness = [](const ModelPart::ElementType& rElement) -> double { return 1.0; };
    }

    std::function<double(const ModelPart::ElementType&)> get_cross_area;
    if (HasVariableInProperties(rModelPart, CROSS_AREA)) {
        get_cross_area = [](const ModelPart::ElementType& rElement) -> double { return rElement.GetProperties()[CROSS_AREA]; };
    } else {
        get_cross_area = [](const ModelPart::ElementType& rElement) -> double { return 1.0; };
    }

    block_for_each(rModelPart.Elements(), [&](auto& rElement) {
        rElement.GetProperties().SetValue(rOutputGradientVariable, rElement.GetGeometry().DomainSize() * get_thickness(rElement) * get_cross_area(rElement));
    });

    KRATOS_CATCH("")
}

void MassResponseUtils::CalculateMassThicknessGradient(
    ModelPart& rModelPart,
    const Variable<double>& rOutputGradientVariable)
{
    CalculateMassGeometricalPropertyGradient(rModelPart, THICKNESS, CROSS_AREA, rOutputGradientVariable);
}

void MassResponseUtils::CalculateMassCrossAreaGradient(
    ModelPart& rModelPart,
    const Variable<double>& rOutputGradientVariable)
{
    CalculateMassGeometricalPropertyGradient(rModelPart, CROSS_AREA, THICKNESS, rOutputGradientVariable);
}

void MassResponseUtils::CalculateMassGeometricalPropertyGradient(
    ModelPart& rModelPart,
    const Variable<double>& rGeometricalPropertyGradientVariable,
    const Variable<double>& rGeometricalConflictingPropertyGradientVariable,
    const Variable<double>& rOutputGradientVariable)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(HasVariableInProperties(rModelPart, DENSITY))
        << "DENSITY is not found in element properties of " << rModelPart.FullName() << ".\n";

    KRATOS_ERROR_IF_NOT(HasVariableInProperties(rModelPart, rGeometricalPropertyGradientVariable))
        << rGeometricalPropertyGradientVariable.Name() << " is not found in element properties of "
        << rModelPart.FullName() << " which is required to compute sensitivities w.r.t. "
        << rGeometricalPropertyGradientVariable.Name() << ".\n";

    KRATOS_ERROR_IF(HasVariableInProperties(rModelPart, rGeometricalConflictingPropertyGradientVariable))
        << rModelPart.FullName() << " has elements with properties having both "
        << rGeometricalPropertyGradientVariable.Name() << " and " << rGeometricalConflictingPropertyGradientVariable.Name()
        << ". Please separate the model part such that either one of them is present in elemental properties.\n";

    block_for_each(rModelPart.Elements(), [&](auto& rElement) {
        rElement.GetProperties().SetValue(rOutputGradientVariable, rElement.GetGeometry().DomainSize() * rElement.GetProperties()[DENSITY]);
    });

    KRATOS_CATCH("")
}

///@}
}