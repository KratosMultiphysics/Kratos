//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "utilities/mortar_utilities.h"
#include "utilities/math_utils.h"
#include "utilities/variable_utils.h"

namespace Kratos
{
namespace MortarUtilities
{
bool LengthCheck(
    const GeometryPointType& rGeometryLine,
    const double Tolerance
    ) {
    const double lx = rGeometryLine[0].X() - rGeometryLine[1].X();
    const double ly = rGeometryLine[0].Y() - rGeometryLine[1].Y();

    const double length = std::sqrt(lx * lx + ly * ly);

    return (length < Tolerance) ? true : false;
}

/***********************************************************************************/
/***********************************************************************************/

bool HeronCheck(const GeometryPointType& rGeometryTriangle)
{
    return HeronCheck(rGeometryTriangle[0], rGeometryTriangle[1], rGeometryTriangle[2]);
}

/***********************************************************************************/
/***********************************************************************************/

bool HeronCheck(
    const PointType& rPointOrig1,
    const PointType& rPointOrig2,
    const PointType& rPointOrig3
    )
{
    const double a = MathUtils<double>::Norm3(rPointOrig1.Coordinates()-rPointOrig2.Coordinates());
    const double b = MathUtils<double>::Norm3(rPointOrig2.Coordinates()-rPointOrig3.Coordinates());
    const double c = MathUtils<double>::Norm3(rPointOrig3.Coordinates()-rPointOrig1.Coordinates());

    const double s = 0.5 * (a + b + c);
    const double A2 = s * (s - a) * (s - b) * (s - c);

    const bool Check = A2 <= 0.0 ? true : false;  // We consider as bad shaped the ones with no area or negative A2 (semiperimeter smaller than any side)

//         // Debug
//         KRATOS_INFO("Check") << Check << " A2: " << A2 << std::endl;
//         if (Check == true) {
//             KRATOS_WARNING("Bad shape") << "Warning:: The triangle is in bad shape" << std::endl;
//             KRATOS_INFO("Mathematica triangle") << "Graphics3D[{EdgeForm[Thick],Triangle[{{" << rPointOrig1.X() << "," << rPointOrig1.Y() << "," << rPointOrig1.Z()  << "},{" << rPointOrig2.X() << "," << rPointOrig2.Y() << "," << rPointOrig2.Z()  << "},{" << rPointOrig3.X() << "," << rPointOrig3.Y() << "," << rPointOrig3.Z()  << "}}]}]" << std::endl;
//         }

    return Check;
}

/***********************************************************************************/
/***********************************************************************************/

void RotatePoint(
    PointType& rPointToRotate,
    const PointType& rPointReferenceRotation,
    const array_1d<double, 3>& rSlaveTangentXi,
    const array_1d<double, 3>& rSlaveTangentEta,
    const bool Inversed
    )
{
    // We move to the (0,0,0)
    PointType aux_point_to_rotate;
    aux_point_to_rotate.Coordinates() = rPointToRotate.Coordinates() - rPointReferenceRotation.Coordinates();

    BoundedMatrix<double, 3, 3> rotation_matrix = ZeroMatrix(3, 3);

    if (Inversed == false) {
        for (IndexType i = 0; i < 3; ++i) {
            rotation_matrix(0, i) = rSlaveTangentXi[i];
            rotation_matrix(1, i) = rSlaveTangentEta[i];
        }
    } else {
        for (IndexType i = 0; i < 3; ++i) {
            rotation_matrix(i, 0) = rSlaveTangentXi[i];
            rotation_matrix(i, 1) = rSlaveTangentEta[i];
        }
    }

    rPointToRotate.Coordinates() = prod(rotation_matrix, aux_point_to_rotate) + rPointReferenceRotation.Coordinates();
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double,3> GaussPointUnitNormal(
    const Vector& rN,
    const GeometryType& rGeometry
    ) {
    array_1d<double,3> r_normal = ZeroVector(3);
    for( IndexType i_node = 0; i_node < rGeometry.PointsNumber(); ++i_node )
        r_normal += rN[i_node] * rGeometry[i_node].FastGetSolutionStepValue(NORMAL);

    const double this_norm = norm_2(r_normal);

    KRATOS_DEBUG_ERROR_IF(this_norm < std::numeric_limits<double>::epsilon()) << "Zero norm r_normal vector. Norm:" << this_norm << std::endl;

    r_normal /= this_norm;

    return r_normal;
}

/***********************************************************************************/
/***********************************************************************************/

void ComputeNodesMeanNormalModelPart(
    ModelPart& rModelPart,
    const bool ComputeConditions
    )
{
    // Check NORMAL is available
    KRATOS_ERROR_IF_NOT(rModelPart.HasNodalSolutionStepVariable(NORMAL)) << "NORMAL is not available on the solution step data variable database" << std::endl;

    // We iterate over nodes
    auto& r_nodes_array = rModelPart.Nodes();
    const auto it_node_begin = r_nodes_array.begin();
    const int num_nodes = static_cast<int>(r_nodes_array.size());

    // Auxiliar zero array
    const array_1d<double, 3> zero_array = ZeroVector(3);

    // Reset NORMAL
    VariableUtils().SetVariable(NORMAL, zero_array, r_nodes_array);

    // Declare auxiliar coordinates
    CoordinatesArrayType aux_coords;

    if (ComputeConditions) {
        // Sum all the nodes normals
        auto& r_conditions_array = rModelPart.Conditions();
        const auto it_cond_begin = r_conditions_array.begin();

        #pragma omp parallel for firstprivate(aux_coords)
        for (int i = 0; i < static_cast<int>(r_conditions_array.size()); ++i) {
            auto it_cond = it_cond_begin + i;
            const GeometryType& r_geometry = it_cond->GetGeometry();

            // Avoid not "flat" conditions
            if (r_geometry.WorkingSpaceDimension() != r_geometry.LocalSpaceDimension() + 1) {
                continue;
            }

            // Set condition normal
            r_geometry.PointLocalCoordinates(aux_coords, r_geometry.Center());
            it_cond->SetValue(NORMAL, r_geometry.UnitNormal(aux_coords));
        }

        // Adding the normal contribution of each node
        for (Condition& r_cond : r_conditions_array) {
            GeometryType& r_geometry = r_cond.GetGeometry();

            // Avoid not "flat" conditions
            if (r_geometry.WorkingSpaceDimension() != r_geometry.LocalSpaceDimension() + 1) {
                continue;
            }

            // Iterate over nodes
            for (NodeType& r_node : r_geometry) {
                r_geometry.PointLocalCoordinates(aux_coords, r_node.Coordinates());
                noalias(r_node.FastGetSolutionStepValue(NORMAL)) += r_geometry.UnitNormal(aux_coords);
            }
        }
    } else {
        auto& r_elements_array = rModelPart.Elements();
        const auto it_elem_begin = r_elements_array.begin();

        #pragma omp parallel for firstprivate(aux_coords)
        for (int i = 0; i < static_cast<int>(r_elements_array.size()); ++i) {
            auto it_elem = it_elem_begin + i;
            const GeometryType& r_geometry = it_elem->GetGeometry();

            // Avoid not "flat" elements
            if (r_geometry.WorkingSpaceDimension() != r_geometry.LocalSpaceDimension() + 1) {
                continue;
            }

            // Set elemition normal
            r_geometry.PointLocalCoordinates(aux_coords, r_geometry.Center());
            it_elem->SetValue(NORMAL, r_geometry.UnitNormal(aux_coords));
        }

        // Adding the normal contribution of each node
        for (Element& r_elem : r_elements_array) {
            GeometryType& r_geometry = r_elem.GetGeometry();

            // Avoid not "flat" elements
            if (r_geometry.WorkingSpaceDimension() != r_geometry.LocalSpaceDimension() + 1) {
                continue;
            }

            // Iterate over nodes
            for (NodeType& r_node : r_geometry) {
                r_geometry.PointLocalCoordinates(aux_coords, r_node.Coordinates());
                noalias(r_node.FastGetSolutionStepValue(NORMAL)) += r_geometry.UnitNormal(aux_coords);
            }
        }
    }

    #pragma omp parallel for
    for (int i = 0; i < num_nodes; ++i) {
        auto it_node = it_node_begin + i;

        array_1d<double, 3>& r_normal = it_node->FastGetSolutionStepValue(NORMAL);
        const double norm_normal = norm_2(r_normal);

        if (norm_normal > std::numeric_limits<double>::epsilon()) r_normal /= norm_normal;
        else KRATOS_ERROR_IF(it_node->Is(INTERFACE)) << "ERROR:: ZERO NORM NORMAL IN NODE: " << it_node->Id() << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void ComputeNodesTangentModelPart(
    ModelPart& rModelPart,
    const Variable<array_1d<double, 3>>* pSlipVariable,
    const double SlipCoefficient,
    const bool SlipAlways
    )
{
    // Check NORMAL is available
    KRATOS_ERROR_IF_NOT(rModelPart.HasNodalSolutionStepVariable(NORMAL)) << "NORMAL is not available on the solution step data variable database" << std::endl;

    // The current dimension
    const std::size_t domain_size = rModelPart.GetProcessInfo()[DOMAIN_SIZE];

    // Checking if we have LM
    const bool has_lm = rModelPart.HasNodalSolutionStepVariable(VECTOR_LAGRANGE_MULTIPLIER);
    KRATOS_ERROR_IF(!has_lm && pSlipVariable == NULL) << "Slip variable not defined" <<  std::endl;

    // We iterate over nodes
    auto& r_nodes_array = rModelPart.Nodes();
    const auto it_node_begin = r_nodes_array.begin();
    const int num_nodes = static_cast<int>(r_nodes_array.size());

    #pragma omp parallel for
    for(int i = 0; i < num_nodes; ++i) {
        auto it_node = it_node_begin + i;

        // Computing only slave nodes
        if (it_node->Is(SLAVE)) {
            if (has_lm && !SlipAlways) {
                ComputeTangentNodeWithLMAndSlip(*it_node, 0, pSlipVariable, SlipCoefficient, domain_size);
            } else {
                ComputeTangentNodeWithSlip(*it_node, 0, pSlipVariable, SlipCoefficient, domain_size);
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void ComputeNodesTangentFromNormalModelPart(ModelPart& rModelPart)
{
    // Check NORMAL is available
    KRATOS_ERROR_IF_NOT(rModelPart.HasNodalSolutionStepVariable(NORMAL)) << "NORMAL is not available on the solution step data variable database" << std::endl;

    // The current dimension
    const std::size_t domain_size = rModelPart.GetProcessInfo()[DOMAIN_SIZE];

    // We iterate over nodes
    auto& r_nodes_array = rModelPart.Nodes();
    const auto it_node_begin = r_nodes_array.begin();
    const int num_nodes = static_cast<int>(r_nodes_array.size());

    #pragma omp parallel for
    for(int i = 0; i < num_nodes; ++i) {
        auto it_node = it_node_begin + i;

        // Computing only slave nodes
        if (it_node->Is(SLAVE)) {
            const array_1d<double, 3>& r_normal = it_node->FastGetSolutionStepValue(NORMAL);
            ComputeTangentsFromNormal(*it_node, r_normal, domain_size);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void ComputeTangentsFromNormal(
    NodeType& rNode,
    const array_1d<double, 3>& rNormal,
    const std::size_t Dimension
    )
{
    array_1d<double, 3> tangent_xi, tangent_eta;
    MathUtils<double>::OrthonormalBasis(rNormal, tangent_xi, tangent_eta);
    if (Dimension == 3) {
        rNode.SetValue(TANGENT_XI, tangent_xi);
        rNode.SetValue(TANGENT_ETA, tangent_eta);
    } else {
        if (std::abs(tangent_xi[2]) > std::numeric_limits<double>::epsilon()) {
            rNode.SetValue(TANGENT_XI, tangent_eta);
        } else {
            rNode.SetValue(TANGENT_XI, tangent_xi);
        }
        tangent_eta[0] = 0.0;
        tangent_eta[1] = 0.0;
        tangent_eta[2] = 1.0;
        rNode.SetValue(TANGENT_ETA, tangent_eta);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void ComputeTangentNodeWithLMAndSlip(
    NodeType& rNode,
    const std::size_t StepLM,
    const Variable<array_1d<double, 3>>* pSlipVariable,
    const double SlipCoefficient,
    const std::size_t Dimension
    )
{
    // Zero tolerance
    const double zero_tolerance = std::numeric_limits<double>::epsilon();

    array_1d<double, 3> tangent;
    const array_1d<double, 3>& r_lm = rNode.FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER, StepLM);
    const array_1d<double, 3>& r_normal = rNode.FastGetSolutionStepValue(NORMAL, StepLM);
    if (rNode.Is(ACTIVE) && (rNode.IsNot(SLIP) || pSlipVariable == NULL)) { // STICK
        if (norm_2(r_lm) > zero_tolerance) { // Non zero LM vector
            noalias(tangent) = r_lm - inner_prod(r_lm, r_normal) * r_normal;
            if (norm_2(tangent) > zero_tolerance) {
                noalias(tangent) = tangent/norm_2(tangent);
                rNode.SetValue(TANGENT_XI, tangent);
            } else {
                ComputeTangentsFromNormal(rNode, r_normal, Dimension);
            }
        } else { // In case of zero LM
            ComputeTangentsFromNormal(rNode, r_normal, Dimension);
        }
    } else if (rNode.Is(ACTIVE) && pSlipVariable != NULL) { // SLIP
        // Nodal length
        const double length = rNode.FastGetSolutionStepValue(NODAL_H);
        // Slip value
        const array_1d<double, 3>& r_slip = rNode.FastGetSolutionStepValue(*pSlipVariable, StepLM)/rNode.GetValue(NODAL_AREA);
        if (norm_2(r_slip) > 1.0e-5 * length) { // Non zero slip vector
            noalias(tangent) = r_slip - inner_prod(r_slip, r_normal) * r_normal;
            if (norm_2(tangent) > zero_tolerance) {
                noalias(tangent) = SlipCoefficient * tangent/norm_2(tangent);
                rNode.SetValue(TANGENT_XI, tangent);
            } else {
                ComputeTangentsFromNormal(rNode, r_normal, Dimension);
            }
        } else if (norm_2(r_lm) > zero_tolerance) { // Non zero LM vector
            const array_1d<double, 3> tangent_component = r_lm - inner_prod(r_lm, r_normal) * r_normal;
            if (norm_2(tangent_component) > zero_tolerance) {
                noalias(tangent) = tangent_component/norm_2(tangent_component);
                rNode.SetValue(TANGENT_XI, tangent);
            } else {
                ComputeTangentsFromNormal(rNode, r_normal, Dimension);
            }
        } else { // In case of zero LM or zero weighted slip
            ComputeTangentsFromNormal(rNode, r_normal, Dimension);
        }
    } else { // Default
        ComputeTangentsFromNormal(rNode, r_normal, Dimension);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void ComputeTangentNodeWithSlip(
    NodeType& rNode,
    const std::size_t StepLM,
    const Variable<array_1d<double, 3>>* pSlipVariable,
    const double SlipCoefficient,
    const std::size_t Dimension
    )
{
    // Zero tolerance
    const double zero_tolerance = std::numeric_limits<double>::epsilon();

    if (pSlipVariable != NULL) { // SLIP
        // Slip value
        const array_1d<double, 3>& r_slip = rNode.FastGetSolutionStepValue(*pSlipVariable, StepLM)/rNode.GetValue(NODAL_AREA);
        if (norm_2(r_slip) > zero_tolerance) { // Non zero slip vector
            const array_1d<double, 3>& r_normal = rNode.FastGetSolutionStepValue(NORMAL, StepLM);
            const array_1d<double, 3> tangent_component = r_slip - inner_prod(r_slip, r_normal) * r_normal;
            if (norm_2(tangent_component) > zero_tolerance) {
                const array_1d<double, 3> tangent = SlipCoefficient * tangent_component/norm_2(tangent_component);
                rNode.SetValue(TANGENT_XI, tangent);
            } else {
                const array_1d<double, 3>& r_normal = rNode.FastGetSolutionStepValue(NORMAL, StepLM);
                array_1d<double, 3> tangent_xi, tangent_eta;
                MathUtils<double>::OrthonormalBasis(r_normal, tangent_xi, tangent_eta);
                if (Dimension == 3) {
                    rNode.SetValue(TANGENT_XI, tangent_xi);
                } else {
                    if (std::abs(tangent_xi[2]) > zero_tolerance) {
                        rNode.SetValue(TANGENT_XI, tangent_eta);
                    } else {
                        rNode.SetValue(TANGENT_XI, tangent_xi);
                    }
                }
            }
        } else { // In case of zero LM or zero weighted slip
            const array_1d<double, 3>& r_normal = rNode.FastGetSolutionStepValue(NORMAL, StepLM);
            array_1d<double, 3> tangent_xi, tangent_eta;
            MathUtils<double>::OrthonormalBasis(r_normal, tangent_xi, tangent_eta);
            if (Dimension == 3) {
                rNode.SetValue(TANGENT_XI, tangent_xi);
            } else {
                if (std::abs(tangent_xi[2]) > zero_tolerance) {
                    rNode.SetValue(TANGENT_XI, tangent_eta);
                } else {
                    rNode.SetValue(TANGENT_XI, tangent_xi);
                }
            }
        }
    } else { // Default
        const array_1d<double, 3>& r_normal = rNode.FastGetSolutionStepValue(NORMAL, StepLM);
        array_1d<double, 3> tangent_xi, tangent_eta;
        MathUtils<double>::OrthonormalBasis(r_normal, tangent_xi, tangent_eta);
        if (Dimension == 3) {
            rNode.SetValue(TANGENT_XI, tangent_xi);
        } else {
            if (std::abs(tangent_xi[2]) > zero_tolerance) {
                rNode.SetValue(TANGENT_XI, tangent_eta);
            } else {
                rNode.SetValue(TANGENT_XI, tangent_xi);
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ResetValue<Variable<double>, MortarUtilitiesSettings::SaveAsHistoricalVariable>(
    ModelPart& rThisModelPart,
    const Variable<double>& rThisVariable
    )
{
    auto& r_nodes_array = rThisModelPart.Nodes();
    VariableUtils().SetVariable(rThisVariable, 0.0, r_nodes_array);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ResetValue<Variable<array_1d<double, 3>>, MortarUtilitiesSettings::SaveAsHistoricalVariable>(
    ModelPart& rThisModelPart,
    const Variable<array_1d<double, 3>>& rThisVariable
    )
{
    auto& r_nodes_array = rThisModelPart.Nodes();
    VariableUtils().SetVariable(rThisVariable, ZeroVector(3), r_nodes_array);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ResetValue<Variable<double>, MortarUtilitiesSettings::SaveAsNonHistoricalVariable>(
    ModelPart& rThisModelPart,
    const Variable<double>& rThisVariable
    )
{
    auto& r_nodes_array = rThisModelPart.Nodes();
    VariableUtils().SetNonHistoricalVariable(rThisVariable, 0.0, r_nodes_array);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ResetValue<Variable<array_1d<double, 3>>, MortarUtilitiesSettings::SaveAsNonHistoricalVariable>(
    ModelPart& rThisModelPart,
    const Variable<array_1d<double, 3>>& rThisVariable
    )
{
    const array_1d<double, 3> zero_array = ZeroVector(3);
    auto& r_nodes_array = rThisModelPart.Nodes();
    VariableUtils().SetNonHistoricalVariable(rThisVariable, zero_array, r_nodes_array);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ResetAuxiliarValue<Variable<double>>(ModelPart& rThisModelPart)
{
    auto& r_nodes_array = rThisModelPart.Nodes();
    VariableUtils().SetNonHistoricalVariable(NODAL_MAUX, 0.0, r_nodes_array);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ResetAuxiliarValue<Variable<array_1d<double, 3>>>(ModelPart& rThisModelPart)
{
    const array_1d<double, 3> zero_array = ZeroVector(3);
    auto& r_nodes_array = rThisModelPart.Nodes();
    VariableUtils().SetNonHistoricalVariable(NODAL_VAUX, zero_array, r_nodes_array);
}

/***********************************************************************************/
/***********************************************************************************/

template< >
const std::string GetAuxiliarVariable<Variable<double>>()
{
    return "NODAL_MAUX";
}

/***********************************************************************************/
/***********************************************************************************/

template< >
const std::string GetAuxiliarVariable<Variable<array_1d<double, 3>>>()
{
    return "NODAL_VAUX";
}

/***********************************************************************************/
/***********************************************************************************/

template< >
double GetAuxiliarValue<Variable<double>>(
    NodeType::Pointer pThisNode,
    const std::size_t iSize
    )
{
    return pThisNode->GetValue(NODAL_MAUX);
}

/***********************************************************************************/
/***********************************************************************************/

template< >
double GetAuxiliarValue<Variable<array_1d<double, 3>>>(
    NodeType::Pointer pThisNode,
    const std::size_t iSize
    )
{
    switch ( iSize ) {
        case 0:
            return pThisNode->GetValue(NODAL_VAUX_X);
        case 1:
            return pThisNode->GetValue(NODAL_VAUX_Y);
        case 2:
            return pThisNode->GetValue(NODAL_VAUX_Z);
        default:
            return 0.0;
    }

    return 0.0;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MatrixValue<Variable<double>, MortarUtilitiesSettings::SaveAsHistoricalVariable>(
    const GeometryType& rThisGeometry,
    const Variable<double>& rThisVariable,
    Matrix& rThisValue
    )
{
    if (rThisValue.size1() != rThisGeometry.size() || rThisValue.size2() != 1)
        rThisValue.resize(rThisGeometry.size(), 1, false);

    for (IndexType i_node = 0; i_node < rThisGeometry.size(); ++i_node)
        rThisValue(i_node, 0) = rThisGeometry[i_node].FastGetSolutionStepValue(rThisVariable);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MatrixValue<Variable<array_1d<double, 3>>, MortarUtilitiesSettings::SaveAsHistoricalVariable>(
    const GeometryType& rThisGeometry,
    const Variable<array_1d<double, 3>>& rThisVariable,
    Matrix& rThisValue
    )
{
    const std::size_t num_nodes = rThisGeometry.size();
    const std::size_t dimension = rThisGeometry.WorkingSpaceDimension();
    if (rThisValue.size1() != num_nodes || rThisValue.size2() != dimension)
        rThisValue.resize(num_nodes, dimension, false);

    for (IndexType i_node = 0; i_node < num_nodes; ++i_node) {
        const auto& r_value = rThisGeometry[i_node].FastGetSolutionStepValue(rThisVariable);
        for (IndexType i_dim = 0; i_dim < dimension; ++i_dim)
            rThisValue(i_node, i_dim) = r_value[i_dim];
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MatrixValue<Variable<double>, MortarUtilitiesSettings::SaveAsNonHistoricalVariable>(
    const GeometryType& rThisGeometry,
    const Variable<double>& rThisVariable,
    Matrix& rThisValue
    )
{
    if (rThisValue.size1() != rThisGeometry.size() || rThisValue.size2() != 1)
        rThisValue.resize(rThisGeometry.size(), 1, false);

    for (IndexType i_node = 0; i_node < rThisGeometry.size(); ++i_node)
        rThisValue(i_node, 0) = rThisGeometry[i_node].GetValue(rThisVariable);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void MatrixValue<Variable<array_1d<double, 3>>, MortarUtilitiesSettings::SaveAsNonHistoricalVariable>(
    const GeometryType& rThisGeometry,
    const Variable<array_1d<double, 3>>& rThisVariable,
    Matrix& rThisValue
    )
{
    const std::size_t num_nodes = rThisGeometry.size();
    const std::size_t dimension = rThisGeometry.WorkingSpaceDimension();
    if (rThisValue.size1() != num_nodes || rThisValue.size2() != dimension)
        rThisValue.resize(num_nodes, dimension, false);

    for (IndexType i_node = 0; i_node < num_nodes; ++i_node) {
        const auto& r_value = rThisGeometry[i_node].GetValue(rThisVariable);
        for (IndexType i_dim = 0; i_dim < dimension; ++i_dim)
            rThisValue(i_node, i_dim) = r_value[i_dim];
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AddValue<Variable<double>, MortarUtilitiesSettings::SaveAsHistoricalVariable>(
    GeometryType& rThisGeometry,
    const Variable<double>& rThisVariable,
    const Matrix& rThisValue
    )
{
    for (IndexType i_node = 0; i_node < rThisGeometry.size(); ++i_node) {
        double& r_aux_value = rThisGeometry[i_node].FastGetSolutionStepValue(rThisVariable);
        #pragma omp atomic
        r_aux_value += rThisValue(i_node, 0);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AddValue<Variable<array_1d<double, 3>>, MortarUtilitiesSettings::SaveAsHistoricalVariable>(
    GeometryType& rThisGeometry,
    const Variable<array_1d<double, 3>>& rThisVariable,
    const Matrix& rThisValue
    )
{
    for (IndexType i_node = 0; i_node < rThisGeometry.size(); ++i_node) {
        auto& r_aux_vector = rThisGeometry[i_node].FastGetSolutionStepValue(rThisVariable);
        for (IndexType i_dim = 0; i_dim < rThisGeometry.WorkingSpaceDimension(); ++i_dim) {
            double& r_aux_value = r_aux_vector[i_dim];
            #pragma omp atomic
            r_aux_value += rThisValue(i_node, i_dim);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AddValue<Variable<double>, MortarUtilitiesSettings::SaveAsNonHistoricalVariable>(
    GeometryType& rThisGeometry,
    const Variable<double>& rThisVariable,
    const Matrix& rThisValue
    )
{
    for (IndexType i_node = 0; i_node < rThisGeometry.size(); ++i_node) {
        double& r_aux_value = rThisGeometry[i_node].GetValue(rThisVariable);
        #pragma omp atomic
        r_aux_value += rThisValue(i_node, 0);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AddValue<Variable<array_1d<double, 3>>, MortarUtilitiesSettings::SaveAsNonHistoricalVariable>(
    GeometryType& rThisGeometry,
    const Variable<array_1d<double, 3>>& rThisVariable,
    const Matrix& rThisValue
    )
{
    for (IndexType i_node = 0; i_node < rThisGeometry.size(); ++i_node) {
        auto& aux_vector = rThisGeometry[i_node].GetValue(rThisVariable);
        for (IndexType i_dim = 0; i_dim < rThisGeometry.WorkingSpaceDimension(); ++i_dim) {
            double& r_aux_value = aux_vector[i_dim];
            #pragma omp atomic
            r_aux_value += rThisValue(i_node, i_dim);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AddAreaWeightedNodalValue<Variable<double>, MortarUtilitiesSettings::SaveAsHistoricalVariable>(
    NodeType::Pointer pThisNode,
    const Variable<double>& rThisVariable,
    const double RefArea,
    const double Tolerance
    )
{
    double area_coeff = pThisNode->GetValue(NODAL_AREA);
    const bool null_area = (std::abs(area_coeff) < RefArea * Tolerance);
    area_coeff = null_area ? 0.0 : 1.0/area_coeff;
    pThisNode->FastGetSolutionStepValue(rThisVariable) += area_coeff * pThisNode->GetValue(NODAL_MAUX);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AddAreaWeightedNodalValue<Variable<array_1d<double, 3>>, MortarUtilitiesSettings::SaveAsHistoricalVariable>(
    NodeType::Pointer pThisNode,
    const Variable<array_1d<double, 3>>& rThisVariable,
    const double RefArea,
    const double Tolerance
    )
{
    double area_coeff = pThisNode->GetValue(NODAL_AREA);
    const bool null_area = (std::abs(area_coeff) < RefArea * Tolerance);
    area_coeff = null_area ? 0.0 : 1.0/area_coeff;
    pThisNode->FastGetSolutionStepValue(rThisVariable) += area_coeff * pThisNode->GetValue(NODAL_VAUX);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AddAreaWeightedNodalValue<Variable<double>, MortarUtilitiesSettings::SaveAsNonHistoricalVariable>(
    NodeType::Pointer pThisNode,
    const Variable<double>& rThisVariable,
    const double RefArea,
    const double Tolerance
    )
{
    double area_coeff = pThisNode->GetValue(NODAL_AREA);
    const bool null_area = (std::abs(area_coeff) < RefArea * Tolerance);
    area_coeff = null_area ? 0.0 : 1.0/area_coeff;
    pThisNode->GetValue(rThisVariable) += area_coeff * pThisNode->GetValue(NODAL_MAUX);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AddAreaWeightedNodalValue<Variable<array_1d<double, 3>>, MortarUtilitiesSettings::SaveAsNonHistoricalVariable>(
    NodeType::Pointer pThisNode,
    const Variable<array_1d<double, 3>>& rThisVariable,
    const double RefArea,
    const double Tolerance
    )
{
    double area_coeff = pThisNode->GetValue(NODAL_AREA);
    const bool null_area = (std::abs(area_coeff) < RefArea * Tolerance);
    area_coeff = null_area ? 0.0 : 1.0/area_coeff;
    pThisNode->GetValue(rThisVariable) += area_coeff * pThisNode->GetValue(NODAL_VAUX);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void UpdateDatabase<Variable<double>, MortarUtilitiesSettings::SaveAsHistoricalVariable>(
    ModelPart& rThisModelPart,
    const Variable<double>& rThisVariable,
    Vector& rDx,
    const std::size_t Index,
    IntMap& rConectivityDatabase
    )
{
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rDx.size()); ++i) {
        auto p_node = rThisModelPart.pGetNode(rConectivityDatabase[i]);
        p_node->FastGetSolutionStepValue(rThisVariable) += rDx[i];
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void UpdateDatabase<Variable<array_1d<double, 3>>, MortarUtilitiesSettings::SaveAsHistoricalVariable>(
    ModelPart& rThisModelPart,
    const Variable<array_1d<double, 3>>& rThisVariable,
    Vector& rDx,
    const std::size_t Index,
    IntMap& rConectivityDatabase
    )
{
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rDx.size()); ++i) {
        auto p_node = rThisModelPart.pGetNode(rConectivityDatabase[i]);
        auto& r_value = p_node->FastGetSolutionStepValue(rThisVariable);
        r_value[Index] += rDx[i];
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void UpdateDatabase<Variable<double>, MortarUtilitiesSettings::SaveAsNonHistoricalVariable>(
    ModelPart& rThisModelPart,
    const Variable<double>& rThisVariable,
    Vector& rDx,
    const std::size_t Index,
    IntMap& rConectivityDatabase
    )
{
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rDx.size()); ++i) {
        auto p_node = rThisModelPart.pGetNode(rConectivityDatabase[i]);
        p_node->GetValue(rThisVariable) += rDx[i];
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void UpdateDatabase<Variable<array_1d<double, 3>>, MortarUtilitiesSettings::SaveAsNonHistoricalVariable>(
    ModelPart& rThisModelPart,
    const Variable<array_1d<double, 3>>& rThisVariable,
    Vector& rDx,
    const std::size_t Index,
    IntMap& rConectivityDatabase
    )
{
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(rDx.size()); ++i) {
        auto p_node = rThisModelPart.pGetNode(rConectivityDatabase[i]);
        auto& value = p_node->GetValue(rThisVariable);
        value[Index] += rDx[i];
    }
}
} // namespace MortarUtilities
} // namespace Kratos
