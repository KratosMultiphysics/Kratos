// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl
//
// ==============================================================================

// System includes

// External includes

// Project includes
#include "custom_utilities/response_functions/thickness_response_function_utility.h"
#include "custom_utilities/search_based_functions.h"
#include "shape_optimization_application.h"
#include "utilities/variable_utils.h"
#include "utilities/integration_utilities.h"
#include "utilities/geometry_utilities.h"

// ==============================================================================

namespace Kratos
{

ThicknessResponseFunctionUtility::ThicknessResponseFunctionUtility(ModelPart& rModelPart, Parameters ResponseSettings)
	: mrModelPart(rModelPart)
{
	const std::size_t domain_size = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];
	KRATOS_ERROR_IF(domain_size != 3) << "ThicknessResponseFunctionUtility can only be used on 3D geometries!" << std::endl;

	mMinThickness = ResponseSettings["min_thickness"].GetDouble();

	const std::string gradient_mode = ResponseSettings["gradient_mode"].GetString();
	if (gradient_mode == "finite_differencing")
	{
		mDelta = ResponseSettings["step_size"].GetDouble();
	}
	else
		KRATOS_ERROR << "Specified gradient_mode '" << gradient_mode << "' not recognized. The only option is: finite_differencing" << std::endl;


	// Create search tree
	mListOfNodesInModelPart.resize(mrModelPart.Nodes().size());
	int counter = 0;
	for (ModelPart::NodesContainerType::iterator node_it = mrModelPart.NodesBegin(); node_it != mrModelPart.NodesEnd(); ++node_it)
	{
		NodeTypePointer pnode = *(node_it.base());
		mListOfNodesInModelPart[counter++] = pnode;
	}

	mpSearchTree = Kratos::shared_ptr<KDTree>(new KDTree(mListOfNodesInModelPart.begin(), mListOfNodesInModelPart.end(), mBucketSize));	
}

double ThicknessResponseFunctionUtility::CalculateValue()
{
	KRATOS_TRY;

	const double value = block_for_each<SumReduction<double>>(mrModelPart.Conditions(), [&](Condition& rCond) {
		array_3d  rQualPointCoords(3,0.0);
		if(FindConditionQaulPoint(rCond,rQualPointCoords))
			return CalculateConditionValue(rCond,rQualPointCoords);
		else
			return 0.0;
		
	});

	mValue = value;

	return mValue;

	KRATOS_CATCH("");
}

void ThicknessResponseFunctionUtility::CalculateGradient()
{
	KRATOS_TRY;
	// First gradients are initialized
	VariableUtils().SetHistoricalVariableToZero(SHAPE_SENSITIVITY, mrModelPart.Nodes());
	VariableUtils().SetHistoricalVariableToZero(NODAL_THICKNESS, mrModelPart.Nodes());


	for (auto& cond_i : mrModelPart.Conditions()){

		array_3d  rQualPointCoords(3,0.0);
		if(FindConditionQaulPoint(cond_i,rQualPointCoords)){
			const double g_i = CalculateConditionValue(cond_i,rQualPointCoords);
			// Compute sensitivities using finite differencing in the three spatial direction
			array_3d gradient(3, 0.0);

			for (auto& node_i : cond_i.GetGeometry()){

				// Apply pertubation in X-direction
				double g_i_after_fd = 0.0;
				node_i.X() += mDelta;
				node_i.X0() += mDelta;
				g_i_after_fd = CalculateConditionValue(cond_i,rQualPointCoords);
				gradient[0] = (g_i_after_fd - g_i) / mDelta;
				node_i.X() -= mDelta;
				node_i.X0() -= mDelta;

				// Apply pertubation in Y-direction
				g_i_after_fd = 0.0;
				node_i.Y() += mDelta;
				node_i.Y0() += mDelta;
				g_i_after_fd = CalculateConditionValue(cond_i,rQualPointCoords);
				gradient[1] = (g_i_after_fd - g_i) / mDelta;
				node_i.Y() -= mDelta;
				node_i.Y0() -= mDelta;

				// Apply pertubation in Z-direction
				g_i_after_fd = 0.0;
				node_i.Z() += mDelta;
				node_i.Z0() += mDelta;
				g_i_after_fd = CalculateConditionValue(cond_i,rQualPointCoords);
				gradient[2] = (g_i_after_fd - g_i) / mDelta;
				node_i.Z() -= mDelta;
				node_i.Z0() -= mDelta;

				// Add to aggregated sensitivities
				noalias(node_i.FastGetSolutionStepValue(SHAPE_SENSITIVITY)) += gradient;
				node_i.FastGetSolutionStepValue(NODAL_THICKNESS) = g_i;
			}			
		}
	}

	KRATOS_CATCH("");
}

double ThicknessResponseFunctionUtility::CalculateConditionValue(const Condition& rFace,array_3d & rQualPointCoords)
{
	// face normal

    auto integration_method = rFace.GetGeometry().GetDefaultIntegrationMethod();
    const auto& integration_points = rFace.GetGeometry().IntegrationPoints( integration_method );
	Vector GaussPtsJDet = ZeroVector(integration_points.size());
    rFace.GetGeometry().DeterminantOfJacobian(GaussPtsJDet, integration_method); 



	double val=0;
	for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
		const array_3d gp_local_coords = integration_points[point_number];
		array_3d gp_global_coords;
		NodeType gp_node(0,gp_global_coords);
		rFace.GetGeometry().GlobalCoordinates(gp_global_coords,gp_local_coords);

		array_3d dist_vec = gp_global_coords-rQualPointCoords;
		double dist = MathUtils<double>::Norm3(dist_vec);
		double delta_dist = (dist - mMinThickness);
		val += (1.0/(1.0+std::exp(20.0*delta_dist))) * integration_points[point_number].Weight() * GaussPtsJDet[point_number];

	}
	
	return val;
}

bool ThicknessResponseFunctionUtility::FindConditionQaulPoint(const Condition& rFace,array_3d & rQualPointCoords)
{

	bool is_cond_active = false;
	Vector face_normal = ZeroVector(3);
	for (auto& node_i : rFace.GetGeometry()){
		face_normal += node_i.FastGetSolutionStepValue(NORMALIZED_SURFACE_NORMAL);
	}
	face_normal/=rFace.GetGeometry().size();


	NodeType center_node(0,rFace.GetGeometry().Center());

	NodeVector neighbor_nodes( mMaxNeighborNodes );
	DoubleVector resulting_squared_distances( mMaxNeighborNodes, 0.0 );
	unsigned int number_of_neighbors = mpSearchTree->SearchInRadius( center_node,
																	mMinThickness,
																	neighbor_nodes.begin(),
																	resulting_squared_distances.begin(),
																	mMaxNeighborNodes );
	double min_neighbour_dist = 2 * mMinThickness;
	for(IndexType npoint_number = 0; npoint_number < number_of_neighbors; ++npoint_number ){
		ModelPart::NodeType& neighbor_node = *neighbor_nodes[npoint_number];
		const array_3d delta = center_node.Coordinates() - neighbor_node.Coordinates();
		double dist = MathUtils<double>::Norm3(delta);
		const array_3d& neighbor_normal = neighbor_node.FastGetSolutionStepValue(NORMALIZED_SURFACE_NORMAL);
		const double projected_normals = inner_prod(face_normal, neighbor_normal);

		if(dist<min_neighbour_dist && projected_normals<0){
			min_neighbour_dist = dist;
			rQualPointCoords = neighbor_node.Coordinates();
			is_cond_active = true;
		}
		

	}
	return is_cond_active;
}

} // namespace Kratos.
