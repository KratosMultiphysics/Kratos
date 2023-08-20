// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:
//
// ==============================================================================

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <functional>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "geometry_utilities.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/key_hash.h"
#include "shape_optimization_application.h"
#include "utilities/variable_utils.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "includes/global_pointer_variables.h"
#include "utilities/element_size_calculator.h"
#include "utilities/atomic_utilities.h"
#include "geometries/geometry_data.h"
#include "processes/calculate_nodal_area_process.h"
#include "spatial_containers/spatial_containers.h"
#include "custom_utilities/heat_method_utilities.h"
#include "linear_solvers/linear_solver.h"




// ==============================================================================

namespace Kratos
{

Vector HeatMethodUtilities::row(Matrix M, unsigned int i)
{
    Vector A = ZeroVector(M.size2());
    for (SizeType k = 0; k < M.size2(); ++k)
    {
        A(k) = M(i, k);
    }
    return A;
}

Vector HeatMethodUtilities::col(Matrix M, unsigned int i)
{
    Vector A = ZeroVector(M.size1());
    for (SizeType k = 0; k < M.size1(); ++k)
    {
        A(k) = M(k, i);
    }
    return A;
}

double HeatMethodUtilities::ComputeAngle(Vector v1, Vector v2)
{   
    KRATOS_TRY
    double length1 = norm_2(v1);
    double length2 = norm_2(v2);
    double dot_product = MathUtils<double>::Dot3(v1, v2);
    return acos((dot_product / (length1 * length2)));
    KRATOS_CATCH("");
}

double HeatMethodUtilities::TriangleArea(array_1d<array_3d, 3> vertices)
{
    Vector edge1, edge2;
    edge1 = vertices(1) - vertices(0);
    edge2 = vertices(2) - vertices(0);
    Vector a = MathUtils<double>::CrossProduct(edge1, edge2);
    return 0.5 * norm_2(a);
}

void HeatMethodUtilities::ElementsArea(Vector& elements_area, Matrix V, Matrix F)
{
    elements_area = ZeroVector(mrModelPart.Conditions().size());
    array_1d<array_3d, 3> vertices;
    for (SizeType i = 0; i < elements_area.size(); ++i)
    {
        vertices(0) = HeatMethodUtilities::row(V, F(i, 0));
        vertices(1) = HeatMethodUtilities::row(V, F(i, 1));
        vertices(2) = HeatMethodUtilities::row(V, F(i, 2));
        elements_area(i) = HeatMethodUtilities::TriangleArea(vertices);
    }
}

void HeatMethodUtilities::VerticesMatrix(Matrix& V)
{
    V = ZeroMatrix(mrModelPart.Nodes().size(), 3);
    int m = 0;
    for (auto& node_i : mrModelPart.Nodes())
    {
        //int i = node_i.Id();
        V(m, 0) = node_i.X();
        V(m, 1) = node_i.Y();
        V(m, 2) = node_i.Z();
        ++m;
    }
}

void HeatMethodUtilities::FacesMatrix(Matrix_int& F)
{
    F = ZeroMatrix_int(mrModelPart.Conditions().size(), 3);
    int m = 0;
    for(auto& cond_i : mrModelPart.Conditions())
    {
        for (SizeType i = 0; i < cond_i.GetGeometry().PointsNumber(); ++i)
        {
            F(m, i) = cond_i.GetGeometry()[i].Id() - 1;
        }
        ++m;
    }
}

void HeatMethodUtilities::SourceNodes(Vector_int& source_nodes)
{
    ModelPart& r_root_model_part = mrModelPart.GetRootModelPart();
    ModelPart& r_start_edge = r_root_model_part.GetSubModelPart("start_edge");
    source_nodes = ZeroVector_int(r_start_edge.Nodes().size());
    int m = 0;
    for (auto& node_i : r_start_edge.Nodes())
    {
        int i = node_i.Id();
        source_nodes(m) = i;
        ++m;
    }
}

void HeatMethodUtilities::BoundaryNodes(Vector_int& boundary_nodes)
{
    ModelPart& r_root_model_part = mrModelPart.GetRootModelPart();
    ModelPart& r_end_edge = r_root_model_part.GetSubModelPart("end_edge");
    boundary_nodes = ZeroVector_int(r_end_edge.Nodes().size());
    int m = 0;
    for (auto& node_i : r_end_edge.Nodes())
    {
        int i = node_i.Id();
        boundary_nodes(m) = i;
        ++m;
    }
}

void HeatMethodUtilities::NodesLabel(Vector_int& nodes_label, Vector_int source_nodes, Vector_int boundary_nodes)
{
    nodes_label = ZeroVector_int(mrModelPart.Nodes().size());
    for (SizeType i = 0; i < source_nodes.size(); ++i)
    {
        nodes_label(source_nodes(i) - 1) = 1;
    }

    for (SizeType i = 0; i < boundary_nodes.size(); ++i)
    {
        nodes_label(boundary_nodes(i) - 1) = 2;
    }
}

void HeatMethodUtilities::HeatEquationMapping(Vector_int& heat_equation_mapping, Vector_int nodes_label, Vector_int boundary_nodes)
{
    const int heat_equation_size = mrModelPart.Nodes().size() - boundary_nodes.size();
    heat_equation_mapping = ZeroVector_int(heat_equation_size);
    int m = 0;
	for (SizeType i = 0; i < mrModelPart.Nodes().size(); ++i)
	{
		if ((nodes_label(i) == 0) || (nodes_label(i) == 1))
		{
			heat_equation_mapping(m) = i;
			++m;
		}
	}
}

void HeatMethodUtilities::HeatEquationInverseMapping(Vector_int& heat_equation_inverse_mapping, Vector_int nodes_label)
{
    heat_equation_inverse_mapping = ZeroVector_int(mrModelPart.Nodes().size());
    int m = 0;
	for (SizeType i = 0; i < mrModelPart.Nodes().size(); ++i)
	{
		int j = i - m;
		if (nodes_label(i) == 2)
		{
			++m;
		}
		heat_equation_inverse_mapping(i) = j;
	}
}

void HeatMethodUtilities::VerticesToFaces(std::vector<std::vector<unsigned int>>& VF, Matrix_int F)
{
    VF.resize(mrModelPart.Nodes().size());
    for (SizeType i = 0; i < F.size1(); ++i)
	{
		for (SizeType j = 0; j < F.size2(); ++j)
		{
			VF[F(i, j)].push_back(i);
		}
	}
}

void HeatMethodUtilities::VerticesToVertices(std::vector<std::vector<unsigned int>>& VV, std::vector<std::vector<unsigned int>> VF, Matrix_int F)
{
    VV.resize(mrModelPart.Nodes().size());
    for (SizeType i = 0; i < VV.size(); ++i)
	{
		for (SizeType j = 0; j < VF[i].size(); ++j)
		{
			for (int k = 0; k < 3; ++k)
			{
				unsigned int value = F(VF[i][j], k);
				if (value != i)
					VV[i].push_back(F(VF[i][j], k));
			}
		}
		std::vector<unsigned int>::iterator ip;
		std::sort(VV[i].begin(), VV[i].end());
		ip = std::unique(VV[i].begin(), VV[i].begin() + VV[i].size());
		VV[i].resize(std::distance(VV[i].begin(), ip));
	}
}

void HeatMethodUtilities::ConstructLaplacian(Matrix& laplacian_matrix, Matrix V, Matrix_int F, std::vector<std::vector<unsigned int>> VV, std::vector<std::vector<unsigned int>> VF)
{
    const SizeType number_of_nodes = mrModelPart.Nodes().size();
    laplacian_matrix = ZeroMatrix(number_of_nodes, number_of_nodes);
    Vector v1, v2, v3;
    for (SizeType i = 0; i < number_of_nodes; ++i) {
		for (SizeType j = 0; j < VV[i].size(); ++j) {
			for (SizeType k = 0; k < VF[i].size(); ++k) {
				std::vector<unsigned int> vec{ F(VF[i][k],0), F(VF[i][k],1), F(VF[i][k],2) };
				std::vector<unsigned int>::iterator it;
				it = std::find(vec.begin(), vec.end(), VV[i][j]);
				if (it != vec.end()) {
					for (int l = 0; l < 3; ++l) {
						if (F(VF[i][k], l) == i) { v1 = HeatMethodUtilities::row(V, F(VF[i][k], l)); }
						else if (F(VF[i][k], l) == VV[i][j]) { v2 = HeatMethodUtilities::row(V, F(VF[i][k], l)); }
						else { v3 = HeatMethodUtilities::row(V, F(VF[i][k], l)); }
					}
					double angle0 = HeatMethodUtilities::ComputeAngle(v1 - v3, v2 - v3);
					laplacian_matrix(i,VV[i][j]) += 0.5 * cos(angle0) / sin(angle0);
                    laplacian_matrix(i,i) -= 0.5 * cos(angle0) / sin(angle0);
				}
			}
		}
	}
}

void HeatMethodUtilities::ConstructK (Matrix& K, Matrix laplacian_matrix, Vector_int heat_equation_mapping, Vector_int heat_equation_inverse_mapping, Vector_int nodes_label) {
    K = ZeroMatrix(heat_equation_mapping.size());
    for (SizeType i = 0; i < K.size1(); ++i) {
		for (SizeType j = 0; j < laplacian_matrix.size1(); ++j) {
			unsigned int m = heat_equation_inverse_mapping(j);
			if ((nodes_label(j) == 0) || (nodes_label(j) == 1)) {
				K(i, m) = laplacian_matrix(heat_equation_mapping(i), j);
			}
		}
    }
}

void HeatMethodUtilities::ConstructM (Matrix& M, Vector elements_area, Vector_int heat_equation_mapping, std::vector<std::vector<unsigned int>> VF) {
    M = ZeroMatrix(heat_equation_mapping.size());
    for (unsigned int i = 0; i < heat_equation_mapping.size(); ++i) {
		unsigned int m = heat_equation_mapping(i);
		for (unsigned int j = 0; j < VF[m].size(); ++j) {
			M(i, i) += elements_area(VF[m][j]) / 3;
		}
	}
}

void HeatMethodUtilities::ConstructU0 (Vector& U0, Vector_int heat_equation_mapping, Vector_int nodes_label) {
    U0 = ZeroVector(heat_equation_mapping.size());
    for (SizeType i = 0; i < heat_equation_mapping.size(); ++i) {
		if (nodes_label(heat_equation_mapping(i)) == 1) {
			U0(i) = 1.0;
		}
    }
}

void HeatMethodUtilities::ComputeLaplacian(LinearSolver<DenseSpace, DenseSpace>& rSolver) {
    KRATOS_TRY

    // checks if object is constructed
    KRATOS_INFO("ShapeOpt") << "Heat method created for centerline mapper..." << std::endl;
    
    
    Matrix V;
    Matrix_int F;
    HeatMethodUtilities::VerticesMatrix(V);
    HeatMethodUtilities::FacesMatrix(F);
    
    Vector_int source_nodes, boundary_nodes;
    HeatMethodUtilities::SourceNodes(source_nodes);
    HeatMethodUtilities::BoundaryNodes(boundary_nodes);

    Vector_int nodes_label;
    HeatMethodUtilities::NodesLabel(nodes_label, source_nodes, boundary_nodes);
    
    Vector_int heat_equation_mapping, heat_equation_inverse_mapping;
    HeatMethodUtilities::HeatEquationMapping(heat_equation_mapping, nodes_label, boundary_nodes);
    HeatMethodUtilities::HeatEquationInverseMapping(heat_equation_inverse_mapping, nodes_label);

    Vector elements_area;
    HeatMethodUtilities::ElementsArea(elements_area, V, F);

    std::vector<std::vector<unsigned int>> VF, VV;
    HeatMethodUtilities::VerticesToFaces(VF, F);
    HeatMethodUtilities::VerticesToVertices(VV, VF, F);

    Matrix laplacian_matrix;
    HeatMethodUtilities::ConstructLaplacian(laplacian_matrix, V, F, VV, VF);

    Matrix K0;
    HeatMethodUtilities::ConstructK(K0, laplacian_matrix, heat_equation_mapping, heat_equation_inverse_mapping, nodes_label);

    Matrix M;
    HeatMethodUtilities::ConstructM(M, elements_area, heat_equation_mapping, VF);
    
    Vector U0;
    HeatMethodUtilities::ConstructU0(U0, heat_equation_mapping, nodes_label);

    double t = 1000;
    Matrix K = M - t * K0;

    Vector U(U0.size());

    std::cout << "K: " << K.size1() << ", " << K.size2() << std::endl;
    std::cout << "U: " << U.size() << std::endl;
    std::cout << "U0: " << U0.size() << std::endl; 

    Matrix A = ZeroMatrix(3);

    std::cout << "SOLVER INFO:  " <<rSolver.Info() << std::endl;
    rSolver.Solve(K, U,U0);
    std::cout << "U:\n";
    for (SizeType i = 0; i < 100; ++i)
    {
        std::cout << ", " << U(i);
    }
    
    Vector U_true = ZeroVector(V.size1());

    for (SizeType i = 0; i < U.size(); ++i) {
        U_true(heat_equation_mapping(i)) = U(i);
    }

    block_for_each(mrModelPart.Nodes(), [&](NodeType &rNode) {
        double& r_heat_distance = rNode.FastGetSolutionStepValue(HEAT_DISTANCE);
        SizeType i = rNode.Id();
        r_heat_distance = U_true(i-1);
        array_3d& r_heat_gradient = rNode.FastGetSolutionStepValue(HEAT_GRADIENT);
        r_heat_gradient(0) = 1;
        r_heat_gradient(1) = 2;
        r_heat_gradient(2) = 3;
    });

    KRATOS_CATCH("");
}

}  // namespace Kratos.