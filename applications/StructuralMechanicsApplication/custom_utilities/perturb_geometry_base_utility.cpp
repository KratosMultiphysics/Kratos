// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Manuel Messmer
//

// System includes

// External includes

// Project includes
#include "custom_utilities/perturb_geometry_base_utility.h"
#include "utilities/builtin_timer.h"
#include "utilities/parallel_utilities.h"
#include "utilities/mortar_utilities.h"
#include "utilities/normal_calculation_utils.h"

namespace Kratos
{

PerturbGeometryBaseUtility::PerturbGeometryBaseUtility( ModelPart& rInitialModelPart, Parameters Settings) :
        mrInitialModelPart(rInitialModelPart),
        mCorrelationLength(Settings["correlation_length"].GetDouble()),
        mTruncationError(Settings["truncation_error"].GetDouble()),
        mEchoLevel(Settings["echo_level"].GetInt()),
        mMaximalDisplacement(Settings["max_displacement"].GetDouble())
    {
        KRATOS_TRY
        NormalCalculationUtils().CalculateUnitNormals<Element>(mrInitialModelPart, true);
        mpPerturbationMatrix = TDenseSpaceType::CreateEmptyMatrixPointer();
        KRATOS_CATCH("")
    }

void PerturbGeometryBaseUtility::ApplyRandomFieldVectorsToGeometry( ModelPart& rThisModelPart, const std::vector<double>& variables )
{
    KRATOS_TRY;
    BuiltinTimer assemble_eigenvectors_timer;
    DenseMatrixType& rPerturbationMatrix = *mpPerturbationMatrix;
    const int num_of_random_variables = variables.size();
    const int num_of_eigenvectors = rPerturbationMatrix.size2();
    const int num_of_nodes = rThisModelPart.NumberOfNodes();

    KRATOS_WARNING_IF("PerturbGeometryBaseUtility",
                    num_of_random_variables != num_of_eigenvectors)
        << "Number of random variables does not match number of eigenvectors: "
        << "Number of random variables: " << num_of_random_variables << ", "
        << "Number of eigenvectors: " << num_of_eigenvectors
        << std::endl;

    std::vector<double> random_field(num_of_nodes,0.0);

    // Fill random field vector
    IndexPartition<unsigned int>(num_of_nodes).for_each(
        [num_of_random_variables, &random_field, &variables, &rPerturbationMatrix](unsigned int i){
            for( int j = 0; j < num_of_random_variables; j++){
                random_field[i] += variables[j]* rPerturbationMatrix(i,j);
            }
        }
    );

    // Shift random field to zero mean
    double shifted_mean = 1.0 / num_of_nodes * std::accumulate(random_field.begin(), random_field.end(), 0.0);
    std::for_each( random_field.begin(), random_field.end(),
        [shifted_mean](double& element) {element -= shifted_mean; } );

    // Find maximum deviation
    double max_disp = *std::max_element(random_field.begin(), random_field.end());
    double min_disp = *std::min_element(random_field.begin(), random_field.end());
    double max_abs_disp = std::max( std::abs(max_disp), std::abs(min_disp) );

    // Scale random field to maximum displacement
    double multiplier = mMaximalDisplacement/max_abs_disp;
    std::for_each( random_field.begin(), random_field.end(),
        [multiplier](double& element) {element *= multiplier; } );

    // Apply random field to geometry
    array_1d<double, 3> normal;
    const auto it_node_original_begin = mrInitialModelPart.NodesBegin();
    const auto it_node_perturb_begin = rThisModelPart.NodesBegin();

    #pragma omp parallel for private(normal)
    for( int i = 0; i < num_of_nodes; i++){
        auto it_node_original = it_node_original_begin + i;
        auto it_node_perturb = it_node_perturb_begin + i;

        // Get normal direction of original model
        normal =  it_node_original->FastGetSolutionStepValue(NORMAL);
        noalias(it_node_perturb->GetInitialPosition().Coordinates()) += normal*random_field[i];
        noalias(it_node_perturb->Coordinates()) += normal*random_field[i];
        //For visualization
        // it_node_perturb->FastGetSolutionStepValue(DISPLACEMENT_X) += normal(0)*random_field[i];
        // it_node_perturb->FastGetSolutionStepValue(DISPLACEMENT_Y) += normal(1)*random_field[i];
        // it_node_perturb->FastGetSolutionStepValue(DISPLACEMENT_Z) += normal(2)*random_field[i];
    }

    KRATOS_INFO_IF("PerturbGeometryBaseUtility: Apply Random Field to Geometry Time", mEchoLevel > 0)
        << assemble_eigenvectors_timer.ElapsedSeconds() << std::endl;

    KRATOS_CATCH("")
}

double PerturbGeometryBaseUtility::CorrelationFunction( ModelPart::NodeType& itNode1, ModelPart::NodeType& itNode2, double CorrelationLength)
{
    array_1d<double, 3> coorrdinate;
    coorrdinate = itNode1.GetInitialPosition().Coordinates() - itNode2.GetInitialPosition().Coordinates();

    double norm = std::sqrt( coorrdinate(0)*coorrdinate(0) + coorrdinate(1)*coorrdinate(1) + coorrdinate(2)*coorrdinate(2) );

    return( std::exp( - norm*norm / (CorrelationLength*CorrelationLength) ) );
}

} // namespace Kratos