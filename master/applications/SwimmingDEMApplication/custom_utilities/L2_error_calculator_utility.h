//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Joaquin Gonzalez-Usua
//

#if !defined(KRATOS_L2_ERROR_CALCULATOR_UTILITY)
#define KRATOS_L2_ERROR_CALCULATOR_UTILITY

// External includes
#include <omp.h>
#include <vector>

// Project includes
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"
#include "includes/kratos_parameters.h"
#include "includes/process_info.h"

namespace Kratos
{
class L2ErrorNormCalculator
{

public:

    typedef Node NodeType;
    typedef Properties PropertiesType;
    typedef Geometry<NodeType> GeometryType;
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;
    typedef ModelPart::ElementsContainerType::iterator  ElementIterator;
    typedef ModelPart::NodesContainerType::iterator     NodeIterator;
    typedef Kratos::Vector ShapeFunctionsType;

    KRATOS_CLASS_POINTER_DEFINITION(L2ErrorNormCalculator);

L2ErrorNormCalculator()
{}

virtual ~L2ErrorNormCalculator(){}

void ComputeDofsErrors(ModelPart& r_model_part)
{
    block_for_each(r_model_part.Nodes(), [&](Node& r_node)
    {
        auto& r_vectorial_error = r_node.FastGetSolutionStepValue(VECTORIAL_ERROR);
        r_vectorial_error = r_node.FastGetSolutionStepValue(VELOCITY) - r_node.FastGetSolutionStepValue(EXACT_VELOCITY);

        auto& r_scalar_error = r_node.FastGetSolutionStepValue(SCALAR_ERROR);
        r_scalar_error = r_node.FastGetSolutionStepValue(PRESSURE) - r_node.FastGetSolutionStepValue(EXACT_PRESSURE);
    });
}

double GetL2VectorErrorNorm(ModelPart& r_model_part)
{
    double total_area = 0.0, result = 0.0;
    const unsigned int dim = r_model_part.GetProcessInfo()[DOMAIN_SIZE];

    for (const auto& r_node : r_model_part.Nodes())
    {
        double local_area = r_node.FastGetSolutionStepValue(NODAL_AREA);
        total_area += local_area;
        double error_x = r_node.FastGetSolutionStepValue(VECTORIAL_ERROR_X);
        double error_y = r_node.FastGetSolutionStepValue(VECTORIAL_ERROR_Y);
        double error_z = 0.0;
        if (dim == 3) error_z += r_node.FastGetSolutionStepValue(VECTORIAL_ERROR_Z);
        double squared_modulus = std::pow(error_x,2)+std::pow(error_y,2)+std::pow(error_z,2);
        result += squared_modulus * local_area;
    }

    return std::sqrt(result/total_area);
}


double GetL2ScalarErrorNorm(ModelPart& r_model_part)
{
    double total_area = 0.0, result = 0.0;

    for (const auto& r_node : r_model_part.Nodes())
    {
        double local_area = r_node.FastGetSolutionStepValue(NODAL_AREA);
        total_area += local_area;
        double error_p = r_node.FastGetSolutionStepValue(SCALAR_ERROR);
        result += std::pow(error_p,2) * local_area;
    }

    return std::sqrt(result/total_area);
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

}; // Class L2ErrorNormCalculator

} // namespace Kratos.

#endif // KRATOS_L2_ERROR_CALCULATOR_UTILITY  defined