// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Navaneeth K Narayanan
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_conditions/hydrostatic_load_condition.h"
#include "utilities/math_utils.h"
#include "utilities/integration_utilities.h"
#include "includes/deprecated_variables.h"
#include "custom_utilities/volume_calculation_under_plane_utility.h"

namespace Kratos
{

//******************************* CONSTRUCTOR ****************************************
//************************************************************************************

HydrostaticLoadCondition::HydrostaticLoadCondition()
{
}

//***********************************************************************************
//***********************************************************************************

HydrostaticLoadCondition::HydrostaticLoadCondition(
    IndexType NewId,
    GeometryType::Pointer pGeometry

    )
    : SurfaceLoadCondition3D(NewId, pGeometry)
{
}

//***********************************************************************************
//***********************************************************************************

HydrostaticLoadCondition::HydrostaticLoadCondition(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : SurfaceLoadCondition3D(NewId, pGeometry, pProperties)
{
}

//********************************* CREATE *******************************************
//************************************************************************************

Condition::Pointer HydrostaticLoadCondition::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_shared<HydrostaticLoadCondition>(NewId, pGeom, pProperties);
}

//***********************************************************************************
//***********************************************************************************

Condition::Pointer HydrostaticLoadCondition::Create(
    IndexType NewId,
    NodesArrayType const &ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_shared<HydrostaticLoadCondition>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

//******************************* DESTRUCTOR *****************************************
//************************************************************************************

HydrostaticLoadCondition::~HydrostaticLoadCondition()
{
}

//***********************************************************************************
//***********************************************************************************

void HydrostaticLoadCondition::CalculateAndSubKpSym(
    Matrix &K,
    const array_1d<double, 3> &rGe,
    const array_1d<double, 3> &rGn,
    const Matrix &rDN_De,
    const Vector &rN,
    const double Pressure,
    const double Weight)
{
    KRATOS_TRY;

    Matrix Kij(3, 3);

    double coeff;
    BoundedMatrix<double, 3, 3> cross_tangent_xi, cross_tangent_eta;
    const SizeType number_of_nodes = GetGeometry().size();
    MakeCrossMatrix(cross_tangent_xi, rGe);
    MakeCrossMatrix(cross_tangent_eta, rGn);

    for (IndexType m = 0; m < number_of_nodes; ++m)
    {
        const IndexType RowIndex = m * 3;

        for (IndexType n = 0; n < number_of_nodes; ++n)
        {
            const IndexType ColIndex = n * 3;

            coeff = Pressure * rN[m] * rDN_De(n, 1) * Weight;
            noalias(Kij) = coeff * cross_tangent_xi;

            coeff = Pressure * rN[n] * rDN_De(m, 1) * Weight; // sym part
            Kij -= coeff * cross_tangent_xi;

            coeff = Pressure * rN[m] * rDN_De(n, 0) * Weight;
            Kij -= coeff * cross_tangent_eta;

            coeff = Pressure * rN[n] * rDN_De(m, 0) * Weight; // sym part
            Kij += coeff * cross_tangent_eta;

            Kij *= 0.5;

            //TAKE CARE: the load correction matrix should be SUBTRACTED not added
            MathUtils<double>::SubtractMatrix(K, Kij, RowIndex, ColIndex);
        }
    }

    KRATOS_CATCH("");
}

//***********************************************************************************
//***********************************************************************************

void HydrostaticLoadCondition::CalculateAndSubKpHydrostatic(
    Matrix &rK,
    const Vector &rN,
    const array_1d<double, 3> &rNormal,
    const double SpecificWeight,
    const array_1d<double, 3> &rW,
    const double Weight)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    Matrix Kij(3, 3);
    double coeff;

    /// Tangent stiffness accounting for change in pressure due to deformation.

    for (unsigned int m = 0; m < number_of_nodes; ++m)
    {
        const unsigned int RowIndex = m * 3;

        for (unsigned int n = 0; n < number_of_nodes; ++n)
        {
            const unsigned int ColIndex = n * 3;

            coeff = -SpecificWeight * rN[m] * rN[n] * Weight;
            DyadicProduct(Kij, rNormal, rW);
            Kij *= -coeff; // (-ve)because normal  is flipped in the implementation following surface_load_condition

            MathUtils<double>::SubtractMatrix(rK, Kij, RowIndex, ColIndex);
        }
    }

    KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void HydrostaticLoadCondition::CalculateAndSubKpHydrostaticSym(
    Matrix &rK,
    const Vector &rN,
    const array_1d<double, 3> &rNormal,
    const double SpecificWeight,
    const array_1d<double, 3> &rW,
    const double Weight)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    Matrix Kij(3, 3);
    Matrix Kij_(3, 3);
    double coeff;

    /// Tangent stiffness accounting for change in pressure due to deformation.

    for (unsigned int m = 0; m < number_of_nodes; ++m)
    {
        const unsigned int RowIndex = m * 3;

        for (unsigned int n = 0; n < number_of_nodes; ++n)
        {
            const unsigned int ColIndex = n * 3;

            coeff = -0.5 * SpecificWeight * rN[m] * rN[n] * Weight;
            DyadicProduct(Kij, rNormal, rW);
            DyadicProduct(Kij_, rW, rNormal);
            Kij += Kij_;
            Kij *= -coeff; // because normal is flipped in the implementation following surface_load_condition

            MathUtils<double>::SubtractMatrix(rK, Kij, RowIndex, ColIndex);
        }
    }

    KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void HydrostaticLoadCondition::DyadicProduct(Matrix &rM, const array_1d<double, 3> &rU, const array_1d<double, 3> &rV)
{
    KRATOS_TRY;

    if (rU.size() == rV.size())
    {

        for (IndexType i = 0; i < rU.size(); ++i)
        {

            for (IndexType j = 0; j < rU.size(); ++j)
            {

                rM(i, j) = rU[i] * rV[j];
            }
        }
    }
    else
        KRATOS_ERROR << "The size of the two input vectors don't match" << std::endl;

    KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

unsigned int HydrostaticLoadCondition::NumberOfCommonElements(NodeType &rNodeM, NodeType &rNodeN)

{

    WeakPointerVector<Element> &node_m_elements = rNodeM.GetValue(NEIGHBOUR_ELEMENTS);
    WeakPointerVector<Element> &node_n_elements = rNodeN.GetValue(NEIGHBOUR_ELEMENTS);
    const SizeType size_m = node_m_elements.size();
    const SizeType size_n = node_n_elements.size();
    std::vector<unsigned int> element_1_ids;
    std::vector<unsigned int> element_2_ids;
    std::vector<unsigned int> element_ids;

    for (IndexType i = 0; i < size_m; ++i)
    {
        element_1_ids.push_back(node_m_elements[i].Id());
    }

    std::sort(element_1_ids.begin(), element_1_ids.end());

    for (IndexType i = 0; i < size_n; ++i)
    {
        element_2_ids.push_back(node_n_elements[i].Id());
    }

    std::sort(element_2_ids.begin(), element_2_ids.end());

    std::set_intersection(element_1_ids.begin(), element_1_ids.end(), element_2_ids.begin(), element_2_ids.end(), std::back_inserter(element_ids));

    return element_ids.size();
}

//***********************************************************************************
//***********************************************************************************

//***********************************************************************************
//***********************************************************************************

void HydrostaticLoadCondition::CalculateAll(
    MatrixType &rLeftHandSideMatrix,
    VectorType &rRightHandSideVector,
    const ProcessInfo &rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag)
{
    KRATOS_TRY;

    GeometryType &geom = GetGeometry();
    const unsigned int number_of_nodes = geom.size();
    const unsigned int mat_size = number_of_nodes * 3;

    //Resizing as needed the LHS
    if (CalculateStiffnessMatrixFlag == true) //calculation of the matrix is required
    {
        if (rLeftHandSideMatrix.size1() != mat_size)
        {
            rLeftHandSideMatrix.resize(mat_size, mat_size, false);
        }

        noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size); //resetting LHS
    }

    // Resizing as needed the RHS
    if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
    {
        if (rRightHandSideVector.size() != mat_size)
        {
            rRightHandSideVector.resize(mat_size, false);
        }

        rRightHandSideVector = ZeroVector(mat_size); //resetting RHS
    }

    array_1d<double, 3> distances_vector;
    bool is_split = false;
    bool is_negative = false;

    for (unsigned int i = 0; i < geom.size(); ++i)
        distances_vector(i) = geom[i].FastGetSolutionStepValue(DISTANCE);

    SetValue(ELEMENTAL_DISTANCES, distances_vector);

    Vector &r_elemental_distances = GetValue(ELEMENTAL_DISTANCES);

    DivideTriangle2D3 triangle_splitter(geom, r_elemental_distances);

    VolumeCalculationUnderPlaneUtility::IsNegativeOrSplit(geom, is_split, is_negative);

    IntegrationMethod integration_method = IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation(geom);

    GeometryType::JacobiansType J;
    J = geom.Jacobian(J, integration_method);

    array_1d<double, 3> ge, gn;
    ge[0] = J[0](0, 0);
    gn[0] = J[0](0, 1);
    ge[1] = J[0](1, 0);
    gn[1] = J[0](1, 1);
    ge[2] = J[0](2, 0);
    gn[2] = J[0](2, 1); // For triangle

    array_1d<double, 3> normal;
    MathUtils<double>::UnitCrossProduct(normal, gn, ge); // Should it be UnitCrossProduct(ge,gn)??

    //const double intersected_area = pGetProperties()->GetValue(FREE_SURFACE_AREA);
    array_1d<double, 3> w = pGetProperties()->GetValue(FREE_SURFACE_NORMAL);
    double norm_w = norm_2(w);

    if (norm_w > std::numeric_limits<double>::epsilon())
        w = w / norm_w;
    else
        w *= 0;

    Vector height_on_nodes(number_of_nodes, 0.0);
    Vector pressure_on_nodes(number_of_nodes, 0.0);
    double specific_wt = pGetProperties()->GetValue(SPECIFIC_WEIGHT);

    for (IndexType i = 0; i < height_on_nodes.size(); ++i)
    {
        if (geom[i].SolutionStepsDataHas(DISTANCE))
        {
            height_on_nodes[i] = geom[i].FastGetSolutionStepValue(DISTANCE);

            pressure_on_nodes[i] = -specific_wt * height_on_nodes[i];
        }
    }

    if (is_split)
    {
        // Call the divide geometry method

        triangle_splitter.GenerateDivision();
        // Calculating actual jacobian

        for (unsigned int i = 0; i < triangle_splitter.mNegativeSubdivisions.size(); ++i)
        {

            IndexedPointGeometryType indexed_subgeom = *(triangle_splitter.mNegativeSubdivisions[i]);

            Node<3>::Pointer node_i = Node<3>::Pointer(new Node<3>(indexed_subgeom[0].Id(), indexed_subgeom[0].X(), indexed_subgeom[0].Y(), indexed_subgeom[0].Z()));

            Node<3>::Pointer node_j = Node<3>::Pointer(new Node<3>(indexed_subgeom[1].Id(), indexed_subgeom[1].X(), indexed_subgeom[1].Y(), indexed_subgeom[1].Z()));
            Node<3>::Pointer node_k = Node<3>::Pointer(new Node<3>(indexed_subgeom[2].Id(), indexed_subgeom[2].X(), indexed_subgeom[2].Y(), indexed_subgeom[2].Z()));

            Triangle3D3<Node<3>>::Pointer subgeom = Triangle3D3<Node<3>>::Pointer(new Triangle3D3<Node<3>>(node_i, node_j, node_k));

            CalculateAllInSplitAndNegativeDistanceConditions(*subgeom, integration_method, ge, gn, normal, pressure_on_nodes, w,
                                                             rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
                                                             CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
        }
    }
    else
    {

        if (is_negative)
        {

            CalculateAllInSplitAndNegativeDistanceConditions(geom, integration_method, ge, gn, normal, pressure_on_nodes, w,
                                                             rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
                                                             CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
        }
    }

    KRATOS_CATCH("")
}

void HydrostaticLoadCondition::CalculateAllInSplitAndNegativeDistanceConditions(
    const GeometryType &rSubGeom,
    IntegrationMethod &rIntegrationMethod,
    const array_1d<double, 3> &rGe,
    const array_1d<double, 3> &rGn,
    const array_1d<double, 3> &rNormal,
    Vector &rPressureOnNodes,
    const array_1d<double, 3> &rW,
    MatrixType &rLeftHandSideMatrix,
    VectorType &rRightHandSideVector,
    const ProcessInfo &rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag)
{

    KRATOS_TRY;
    GeometryType &geom = GetGeometry();
    const unsigned int number_of_nodes = geom.size();
    array_1d<double, 3> global_coords = ZeroVector(3);
    array_1d<double, 3> local_coords = ZeroVector(3);
    Vector N(geom.size());

    Vector jacobians_values;
    rSubGeom.DeterminantOfJacobian(jacobians_values, rIntegrationMethod);
    Matrix DN_De;
    const GeometryType::IntegrationPointsArrayType &integration_points = rSubGeom.IntegrationPoints(rIntegrationMethod);

    const double specific_wt = pGetProperties()->GetValue(SPECIFIC_WEIGHT);

    for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number)
    {
        rSubGeom.GlobalCoordinates(global_coords, integration_points[point_number].Coordinates());
        geom.PointLocalCoordinates(local_coords, global_coords);
        geom.ShapeFunctionsValues(N, local_coords);
        geom.ShapeFunctionsLocalGradients(DN_De, local_coords);

        const double det_j = jacobians_values[point_number];
        const double integration_weight = GetIntegrationWeight(integration_points, point_number, det_j);
        double integration_weight_without_detj = integration_weight / det_j;

        // Calculating the pressure on the gauss point
        double pressure = 0.0;
        for (IndexType ii = 0; ii < number_of_nodes; ++ii)
        {
            pressure += N[ii] * rPressureOnNodes[ii];
        }

        // LEFT HAND SIDE MATRIX
        if (CalculateStiffnessMatrixFlag == true)
        {
            if (fabs(pressure) > std::numeric_limits<double>::epsilon())
            {

                if (pGetProperties()->GetValue(USE_HYDROSTATIC_MATRIX))
                    CalculateAndSubKpSym(rLeftHandSideMatrix, rGe, rGn, DN_De, N, pressure, integration_weight_without_detj);
                else
                    CalculateAndSubKp(rLeftHandSideMatrix, rGe, rGn, DN_De, N, pressure, integration_weight_without_detj);
            }

            if ((specific_wt > std::numeric_limits<double>::epsilon()) && pGetProperties()->GetValue(USE_HYDROSTATIC_MATRIX))
            //if (specific_wt > std::numeric_limits<double>::epsilon())
            {
                CalculateAndSubKpHydrostaticSym(rLeftHandSideMatrix, N, rNormal, specific_wt, rW, integration_weight);
            }
        }

        // RIGHT HAND SIDE VECTOR
        if (CalculateResidualVectorFlag == true)
        {
            if (pressure != 0.0)
            {
                if (pGetProperties()->GetValue(DO_UPDATE_FROM_DEL_VOL))
                {
                    double del_p;
                    double del_vol = pGetProperties()->GetValue(FLUID_VOLUME) - pGetProperties()->GetValue(CURRENT_FLUID_VOLUME);
                    double free_surface_area = pGetProperties()->GetValue(FREE_SURFACE_AREA);
                    double specific_weight = pGetProperties()->GetValue(SPECIFIC_WEIGHT);
                    if (pGetProperties()->GetValue(FLUID_VOLUME) < std::numeric_limits<double>::epsilon() || free_surface_area < std::numeric_limits<double>::epsilon())
                        del_p = 0.0;
                    else
                        del_p = del_vol * specific_weight / free_surface_area;

                    pressure = pressure + del_p;
                }

                CalculateAndAddPressureForce(rRightHandSideVector, N, rNormal, pressure, integration_weight, rCurrentProcessInfo);

                if (pGetProperties()->GetValue(ADD_RHS_FOR_RANK_ONE_UPDATE))
                {

                    CalculateAndAddPressureForce(rRightHandSideVector, N, rNormal, 1.0, integration_weight, rCurrentProcessInfo);
                }
            }
        }
    }
    KRATOS_CATCH("")
}

} // namespace Kratos
