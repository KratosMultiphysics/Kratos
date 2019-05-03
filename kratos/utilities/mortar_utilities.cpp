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

void ComputeNodesMeanNormalModelPart(ModelPart& rModelPart)
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
    VariableUtils().SetVectorVar(NORMAL, zero_array, r_nodes_array);

    // Sum all the nodes normals
    auto& r_conditions_array = rModelPart.Conditions();
    const auto it_cond_begin = r_conditions_array.begin();

    // Declare auxiliar coordinates
    CoordinatesArrayType aux_coords;

    #pragma omp parallel for firstprivate(aux_coords)
    for(int i = 0; i < static_cast<int>(r_conditions_array.size()); ++i) {
        auto it_cond = it_cond_begin + i;
        const GeometryType& r_geometry = it_cond->GetGeometry();

        // Set condition normal
        r_geometry.PointLocalCoordinates(aux_coords, r_geometry.Center());
        it_cond->SetValue(NORMAL, r_geometry.UnitNormal(aux_coords));
    }

    // Adding the normal contribution of each node
    for(Condition& r_cond : r_conditions_array) {
        GeometryType& r_geometry = r_cond.GetGeometry();

        // Iterate over nodes
        for (NodeType& r_node : r_geometry) {
            r_geometry.PointLocalCoordinates(aux_coords, r_node.Coordinates());
            noalias(r_node.FastGetSolutionStepValue(NORMAL)) += r_geometry.UnitNormal(aux_coords);
        }
    }

    #pragma omp parallel for
    for(int i = 0; i < num_nodes; ++i) {
        auto it_node = it_node_begin + i;

        array_1d<double, 3>& r_normal = it_node->FastGetSolutionStepValue(NORMAL);
        const double norm_normal = norm_2(r_normal);

        if (norm_normal > std::numeric_limits<double>::epsilon()) r_normal /= norm_normal;
        else KRATOS_ERROR_IF(it_node->Is(INTERFACE)) << "ERROR:: ZERO NORM NORMAL IN NODE: " << it_node->Id() << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ResetValue<Variable<double>, Historical>(
    ModelPart& rThisModelPart,
    Variable<double>& rThisVariable
    )
{
    auto& r_nodes_array = rThisModelPart.Nodes();
    VariableUtils().SetScalarVar(rThisVariable, 0.0, r_nodes_array);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ResetValue<Variable<array_1d<double, 3>>, Historical>(
    ModelPart& rThisModelPart,
    Variable<array_1d<double, 3>>& rThisVariable
    )
{
    auto& r_nodes_array = rThisModelPart.Nodes();
    VariableUtils().SetVectorVar(rThisVariable, ZeroVector(3), r_nodes_array);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ResetValue<Variable<double>, NonHistorical>(
    ModelPart& rThisModelPart,
    Variable<double>& rThisVariable
    )
{
    auto& r_nodes_array = rThisModelPart.Nodes();
    VariableUtils().SetNonHistoricalVariable(rThisVariable, 0.0, r_nodes_array);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ResetValue<Variable<array_1d<double, 3>>, NonHistorical>(
    ModelPart& rThisModelPart,
    Variable<array_1d<double, 3>>& rThisVariable
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
Variable<double> GetAuxiliarVariable<Variable<double>>()
{
    return NODAL_MAUX;
}

/***********************************************************************************/
/***********************************************************************************/

template< >
Variable<array_1d<double, 3>> GetAuxiliarVariable<Variable<array_1d<double, 3>>>()
{
    return NODAL_VAUX;
}

/***********************************************************************************/
/***********************************************************************************/

template< >
double GetAuxiliarValue<Variable<double>>(
    NodeType::Pointer pThisNode,
    unsigned int iSize
    )
{
    return pThisNode->GetValue(NODAL_MAUX);
}

/***********************************************************************************/
/***********************************************************************************/

template< >
double GetAuxiliarValue<Variable<array_1d<double, 3>>>(
    NodeType::Pointer pThisNode,
    unsigned int iSize
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
void MatrixValue<Variable<double>, Historical>(
    GeometryType& rThisGeometry,
    Variable<double>& rThisVariable,
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
void MatrixValue<Variable<array_1d<double, 3>>, Historical>(
    GeometryType& rThisGeometry,
    Variable<array_1d<double, 3>>& rThisVariable,
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
void MatrixValue<Variable<double>, NonHistorical>(
    GeometryType& rThisGeometry,
    Variable<double>& rThisVariable,
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
void MatrixValue<Variable<array_1d<double, 3>>, NonHistorical>(
    GeometryType& rThisGeometry,
    Variable<array_1d<double, 3>>& rThisVariable,
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
void AddValue<Variable<double>, Historical>(
    GeometryType& rThisGeometry,
    Variable<double>& rThisVariable,
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
void AddValue<Variable<array_1d<double, 3>>, Historical>(
    GeometryType& rThisGeometry,
    Variable<array_1d<double, 3>>& rThisVariable,
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
void AddValue<Variable<double>, NonHistorical>(
    GeometryType& rThisGeometry,
    Variable<double>& rThisVariable,
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
void AddValue<Variable<array_1d<double, 3>>, NonHistorical>(
    GeometryType& rThisGeometry,
    Variable<array_1d<double, 3>>& rThisVariable,
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
void AddAreaWeightedNodalValue<Variable<double>, Historical>(
    NodeType::Pointer pThisNode,
    Variable<double>& rThisVariable,
    const double RefArea,
    const double Tolerance
    )
{
    double area_coeff = pThisNode->GetValue(NODAL_AREA);
    const bool null_area = (std::abs(area_coeff) < RefArea * Tolerance);
#ifdef KRATOS_DEBUG
    KRATOS_WARNING_IF("WARNING:: NODE OF NULL AREA.", null_area && pThisNode->Is(SLAVE)) << " ID: " << pThisNode->Id() << std::endl;
#endif
    area_coeff = null_area ? 0.0 : 1.0/area_coeff;
    double& r_aux_value = pThisNode->FastGetSolutionStepValue(rThisVariable);
    #pragma omp atomic
    r_aux_value += area_coeff * pThisNode->GetValue(NODAL_MAUX);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AddAreaWeightedNodalValue<Variable<array_1d<double, 3>>, Historical>(
    NodeType::Pointer pThisNode,
    Variable<array_1d<double, 3>>& rThisVariable,
    const double RefArea,
    const double Tolerance
    )
{
    double area_coeff = pThisNode->GetValue(NODAL_AREA);
    const bool null_area = (std::abs(area_coeff) < RefArea * Tolerance);
#ifdef KRATOS_DEBUG
    KRATOS_WARNING_IF("WARNING:: NODE OF NULL AREA.", null_area && pThisNode->Is(SLAVE)) << " ID: " << pThisNode->Id() << std::endl;
#endif
    area_coeff = null_area ? 0.0 : 1.0/area_coeff;
    auto& aux_vector = pThisNode->FastGetSolutionStepValue(rThisVariable);
    const auto& nodal_vaux = pThisNode->GetValue(NODAL_VAUX);
    for (IndexType i = 0; i < 3; ++i) {
        double& r_aux_value = aux_vector[i];
        #pragma omp atomic
        r_aux_value += area_coeff * nodal_vaux[i];
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AddAreaWeightedNodalValue<Variable<double>, NonHistorical>(
    NodeType::Pointer pThisNode,
    Variable<double>& rThisVariable,
    const double RefArea,
    const double Tolerance
    )
{
    double area_coeff = pThisNode->GetValue(NODAL_AREA);
    const bool null_area = (std::abs(area_coeff) < RefArea * Tolerance);
#ifdef KRATOS_DEBUG
    KRATOS_WARNING_IF("WARNING:: NODE OF NULL AREA.", null_area && pThisNode->Is(SLAVE)) << " ID: " << pThisNode->Id() << std::endl;
#endif
    area_coeff = null_area ? 0.0 : 1.0/area_coeff;
    double& r_aux_value = pThisNode->GetValue(rThisVariable);
    #pragma omp atomic
    r_aux_value += area_coeff * pThisNode->GetValue(NODAL_MAUX);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AddAreaWeightedNodalValue<Variable<array_1d<double, 3>>, NonHistorical>(
    NodeType::Pointer pThisNode,
    Variable<array_1d<double, 3>>& rThisVariable,
    const double RefArea,
    const double Tolerance
    )
{
    double area_coeff = pThisNode->GetValue(NODAL_AREA);
    const bool null_area = (std::abs(area_coeff) < RefArea * Tolerance);
#ifdef KRATOS_DEBUG
    KRATOS_WARNING_IF("WARNING:: NODE OF NULL AREA.", null_area && pThisNode->Is(SLAVE)) << " ID: " << pThisNode->Id() << std::endl;
#endif
    area_coeff = null_area ? 0.0 : 1.0/area_coeff;
    auto& aux_vector = pThisNode->GetValue(rThisVariable);
    const auto& nodal_vaux = pThisNode->GetValue(NODAL_VAUX);
    for (IndexType i = 0; i < 3; ++i) {
        double& r_aux_value = aux_vector[i];
        #pragma omp atomic
        r_aux_value += area_coeff * nodal_vaux[i];
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void UpdateDatabase<Variable<double>, Historical>(
    ModelPart& rThisModelPart,
    Variable<double>& rThisVariable,
    Vector& Dx,
    unsigned int Index,
    IntMap& ConectivityDatabase
    )
{
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(Dx.size()); ++i) {
        auto p_node = rThisModelPart.pGetNode(ConectivityDatabase[i]);
        p_node->FastGetSolutionStepValue(rThisVariable) += Dx[i];
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void UpdateDatabase<Variable<array_1d<double, 3>>, Historical>(
    ModelPart& rThisModelPart,
    Variable<array_1d<double, 3>>& rThisVariable,
    Vector& Dx,
    unsigned int Index,
    IntMap& ConectivityDatabase
    )
{
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(Dx.size()); ++i) {
        auto p_node = rThisModelPart.pGetNode(ConectivityDatabase[i]);
        auto& r_value = p_node->FastGetSolutionStepValue(rThisVariable);
        r_value[Index] += Dx[i];
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void UpdateDatabase<Variable<double>, NonHistorical>(
    ModelPart& rThisModelPart,
    Variable<double>& rThisVariable,
    Vector& Dx,
    unsigned int Index,
    IntMap& ConectivityDatabase
    )
{
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(Dx.size()); ++i) {
        auto p_node = rThisModelPart.pGetNode(ConectivityDatabase[i]);
        p_node->GetValue(rThisVariable) += Dx[i];
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void UpdateDatabase<Variable<array_1d<double, 3>>, NonHistorical>(
    ModelPart& rThisModelPart,
    Variable<array_1d<double, 3>>& rThisVariable,
    Vector& Dx,
    unsigned int Index,
    IntMap& ConectivityDatabase
    )
{
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(Dx.size()); ++i) {
        auto p_node = rThisModelPart.pGetNode(ConectivityDatabase[i]);
        auto& value = p_node->GetValue(rThisVariable);
        value[Index] += Dx[i];
    }
}
} // namespace MortarUtilities
} // namespace Kratos
