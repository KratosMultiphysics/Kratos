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
void test_function();
void HeatMethodUtilities::hmt() {
    std::cout << "hmt worked!" << std::endl;
}
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
        source_nodes(m) = i - 1;
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
        boundary_nodes(m) = i - 1;
        ++m;
    }
}

void HeatMethodUtilities::NodesLabel(Vector_int& nodes_label, Vector_int source_nodes, Vector_int boundary_nodes)
{
    nodes_label = ZeroVector_int(mrModelPart.Nodes().size());
    for (SizeType i = 0; i < source_nodes.size(); ++i)
    {
        nodes_label(source_nodes(i)) = 1;
    }

    for (SizeType i = 0; i < boundary_nodes.size(); ++i)
    {
        nodes_label(boundary_nodes(i)) = 2;
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

void HeatMethodUtilities::ConstructLaplacian(CompressedMatrix& laplacian_matrix, Matrix V, Matrix_int F, std::vector<std::vector<unsigned int>> VV, std::vector<std::vector<unsigned int>> VF)
{
    const SizeType number_of_nodes = mrModelPart.Nodes().size();
    CompressedMatrix LC(number_of_nodes, number_of_nodes);
    laplacian_matrix = LC;
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

void HeatMethodUtilities::ConstructK (CompressedMatrix& K, CompressedMatrix laplacian_matrix, Vector_int heat_equation_mapping, Vector_int heat_equation_inverse_mapping, Vector_int nodes_label) {
    CompressedMatrix K_hehe(heat_equation_mapping.size(), heat_equation_mapping.size());
    K = K_hehe;
    for (SizeType i = 0; i < K.size1(); ++i) {
		for (SizeType j = 0; j < laplacian_matrix.size1(); ++j) {
			unsigned int m = heat_equation_inverse_mapping(j);
			if ((nodes_label(j) == 0) || (nodes_label(j) == 1)) {
				K(i, m) = laplacian_matrix(heat_equation_mapping(i), j);
			}
		}
    }
}

void HeatMethodUtilities::ConstructM (CompressedMatrix& M, Vector elements_area, Vector_int heat_equation_mapping, std::vector<std::vector<unsigned int>> VF) {
    CompressedMatrix M_hehe(heat_equation_mapping.size(), heat_equation_mapping.size());
    M = M_hehe;
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

void HeatMethodUtilities::ComputeHeatGradient (Matrix& heat_gradient, Matrix V, Matrix F, Vector HeatField) {
    heat_gradient = ZeroMatrix(F.size1(), 3);
    for (SizeType i = 0; i < heat_gradient.size1(); ++i) {
        int i1 = F(i, 0), i2 = F(i, 1), i3 = F(i, 2);
        array_1d<Vector, 3> vertices, edges;
        vertices(0) = HeatMethodUtilities::row(V, i1);
		vertices(1) = HeatMethodUtilities::row(V, i2);
		vertices(2) = HeatMethodUtilities::row(V, i3);
		Vector elementCenter = (vertices(0) + vertices(1) + vertices(2)) / 3;
		edges(0) = vertices(2) - vertices(1);
		edges(1) = vertices(0) - vertices(2);
		edges(2) = vertices(1) - vertices(0);
        Vector normal_vector = MathUtils<double>::CrossProduct(edges(0), edges(1));
        normal_vector /= norm_2(normal_vector);
        Vector normal_vector_at_vertice_1 = mrModelPart.Nodes()[i1 + 1].FastGetSolutionStepValue(NORMALIZED_SURFACE_NORMAL);
        if (MathUtils<double>::Dot3(normal_vector, normal_vector_at_vertice_1) < 0)
        {
            normal_vector *= -1;
        }
        for (int j = 0; j < 3; ++j) {
            Vector vec = vertices(j) - elementCenter;
            Vector test_edge = MathUtils<double>::CrossProduct(edges(j), vec);
			if (MathUtils<double>::Dot3(normal_vector, edges(j)) < 0)
			{
                edges(j) *= -1;
            }
            array_3d grad = HeatField(F(i, j)) * MathUtils<double>::CrossProduct(normal_vector, edges(j));
            for (int k = 0; k <3; ++k) {
                heat_gradient(i, k) += grad(k);
            }
        }
        for (unsigned int k = 0; k <3 ; ++k) {
            heat_gradient(i, k) /= -norm_2(row(heat_gradient, i));
        }
	}
}

void HeatMethodUtilities::ComputeHeatDivergence (Vector& heat_divergence, Matrix heat_gradient, Matrix V, Matrix F, std::vector<std::vector<unsigned int>> VF) {
    heat_divergence = ZeroVector(mrModelPart.Nodes().size());
	for (SizeType i = 0; i < heat_divergence.size(); ++i)
	{
		array_1d<Vector, 3> vertices;
		vertices(0) = HeatMethodUtilities::row(V, i);
		for (SizeType j = 0; j < VF[i].size(); ++j) {
			int index = VF[i][j];
			int m = 1;
			for (SizeType k = 0; k < 3; ++k) {
				if (F(index, k) != i) {
					vertices(m) = HeatMethodUtilities::row(V, (F(index, k)));
					m ++;
				}
			}
			array_1d<Vector, 2> edges;
			edges(0) = vertices(1) - vertices(0);
			edges(1) = vertices(2) - vertices(0);
			array_1d<double, 2> inte_angles;
			inte_angles(0) = HeatMethodUtilities::ComputeAngle(vertices(0) - vertices(2), vertices(1) - vertices(2));
			inte_angles(1) = HeatMethodUtilities::ComputeAngle(vertices(0) - vertices(1), vertices(2) - vertices(1));
            for (int l = 0; l < 2; l++) {
                heat_divergence(i) += 0.5 * cotan(inte_angles(l)) * dot(edges(l), row(heat_gradient, index));
            }
		}
	}
}

void HeatMethodUtilities::ComputeHeatField(LinearSolver<SparseSpaceType, LocalSpaceType>& rSolver) {
    KRATOS_TRY

    // checks if object is constructed
    KRATOS_INFO("ShapeOpt") << "Heat method created for centerline mapper..." << std::endl;
    
    HeatMethodUtilities hmt1(mrModelPart);
    Matrix V;
    Matrix_int F;
    hmt1.VerticesMatrix(V);
    hmt1.FacesMatrix(F);
    
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

    CompressedMatrix laplacian_matrix;
    HeatMethodUtilities::ConstructLaplacian(laplacian_matrix, V, F, VV, VF);

    CompressedMatrix K0;
    HeatMethodUtilities::ConstructK(K0, laplacian_matrix, heat_equation_mapping, heat_equation_inverse_mapping, nodes_label);

    CompressedMatrix M;
    HeatMethodUtilities::ConstructM(M, elements_area, heat_equation_mapping, VF);
    M *= 1/3;
    
    Vector U0;
    HeatMethodUtilities::ConstructU0(U0, heat_equation_mapping, nodes_label);

    double t = 100000;
    CompressedMatrix K = M - t * K0;

    Vector U(U0.size());

    std::cout << "SOLVER INFO:  " <<rSolver.Info() << std::endl;
    rSolver.Solve(K, U, U0);
    
    Vector heat_field = ZeroVector(V.size1());

    for (SizeType i = 0; i < U.size(); ++i) {
        heat_field(heat_equation_mapping(i)) = U(i);
    }

    block_for_each(mrModelPart.Nodes(), [&](NodeType &rNode) {
        double& r_heat_gradient = rNode.FastGetSolutionStepValue(HEAT_FIELD);
        SizeType i = rNode.Id();
        r_heat_gradient = heat_field(i-1);
    });
    
    Matrix heat_gradient_hehe;
    HeatMethodUtilities::ComputeHeatGradient(heat_gradient_hehe, V, F, heat_field);

    Vector heat_divergence;
    HeatMethodUtilities::ComputeHeatDivergence(heat_divergence, heat_gradient_hehe, V, F, VF);

    Vector heat_distance(heat_divergence.size());

    std::cout << "SOLVER INFO:  " <<rSolver.Info() << std::endl;
    rSolver.Solve(laplacian_matrix, heat_distance, heat_divergence);

    double s = 0;
	for (SizeType i = 0; i < source_nodes.size(); ++i)
	{
		s += heat_distance(source_nodes(i));
	}
	s = s / source_nodes.size();

	for (SizeType i = 0; i < heat_distance.size(); ++i)
	{
		heat_distance(i) -= s;
	}

    for (SizeType i = 0; i < nodes_label.size(); ++i)
    {
        std::cout << "heat_distance " << i << ":  " << heat_distance(i) << std::endl;
    }

    for (SizeType i = 0; i < source_nodes.size(); ++i)
    {
        std::cout << "source nodes:  " << source_nodes(i) << std::endl;
    }

    for (SizeType i = 0; i < source_nodes.size(); ++i)
    {
        std::cout << "heat at source:  " << heat_field(source_nodes(i)) << std::endl;
    }

    for (SizeType i = 0; i < source_nodes.size(); ++i)
    {
        std::cout << "distance at source:  " << heat_distance(source_nodes(i)) << std::endl;
    }

    block_for_each(mrModelPart.Nodes(), [&](NodeType &rNode) {
        double& r_heat_distance = rNode.FastGetSolutionStepValue(HEAT_DISTANCE);
        SizeType i = rNode.Id();
        r_heat_distance = heat_distance(i-1);
        // array_3d& r_heat_gradient = rNode.FastGetSolutionStepValue(HEAT_GRADIENT);
        // r_heat_gradient(0) = 1;
        // r_heat_gradient(1) = 2;
        // r_heat_gradient(2) = 3;
    });

    test_function();
    HeatMethodUtilities aa(mrModelPart);
    aa.hmt();
    KRATOS_CATCH("");
}

void test_function() {
    std::cout << "test function worked!" << std::endl;
}

}  // namespace Kratos.