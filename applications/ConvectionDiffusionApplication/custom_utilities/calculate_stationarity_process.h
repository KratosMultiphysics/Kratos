//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Aitor Bazán Escoda
//

#if !defined(KRATOS_CALCULATE_STATIONARITY_PROCESS)
#define KRATOS_CALCULATE_STATIONARITY_PROCESS

#include <string>
#include <iostream>
#include <algorithm>

// External includes
#include <omp.h>
#include <vector>
#include <tuple>


// Project includes
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"
#include "includes/kratos_parameters.h"
#include "includes/process_info.h"
#include "convection_diffusion_application.h"

namespace Kratos
{
class CalculateStationarityProcess
{

public:

    typedef Node < 3 > NodeType;
    typedef Properties PropertiesType;
    typedef Geometry<NodeType> GeometryType;
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;
    typedef ModelPart::ElementsContainerType::iterator  ElementIterator;
    typedef ModelPart::NodesContainerType::iterator     NodeIterator;
    typedef Kratos::Vector ShapeFunctionsType;
    
    KRATOS_CLASS_POINTER_DEFINITION(CalculateStationarityProcess);

CalculateStationarityProcess()
{}

virtual ~CalculateStationarityProcess()
{}


void ComputeDofsErrorsCD(ModelPart& r_model_part)
{
    block_for_each(r_model_part.Nodes(), [&](Node<3>& r_node)
    {
        auto& r_vectorial_error = r_node.FastGetSolutionStepValue(VECTORIAL_STAT_ERROR);
        r_vectorial_error = r_node.FastGetSolutionStepValue(VELOCITY) - r_node.FastGetSolutionStepValue(VELOCITY,1);
        // KRATOS_WATCH(r_vectorial_error)
        auto& r_scalar_error = r_node.FastGetSolutionStepValue(SCALAR_STAT_ERROR);
        r_scalar_error = r_node.FastGetSolutionStepValue(TEMPERATURE) - r_node.FastGetSolutionStepValue(TEMPERATURE,1);
        // KRATOS_WATCH(r_scalar_error)
    });
}

double GetL2VectorErrorNormCD(ModelPart& r_model_part)
{
    double total_area = 0.0, result = 0.0;
    const unsigned int dim = r_model_part.GetProcessInfo()[DOMAIN_SIZE];

    for (const auto& r_node : r_model_part.Nodes())
    {
        double local_area = r_node.FastGetSolutionStepValue(NODAL_AREA);
        total_area += local_area;
        // KRATOS_WATCH(local_area)
        double error_x = r_node.FastGetSolutionStepValue(VECTORIAL_STAT_ERROR_X);
        // KRATOS_WATCH(error_x)
        double error_y = r_node.FastGetSolutionStepValue(VECTORIAL_STAT_ERROR_Y);
        // KRATOS_WATCH(error_y)
        double error_z = 0.0;
        if (dim == 3) error_z += r_node.FastGetSolutionStepValue(VECTORIAL_STAT_ERROR_Z);
        double squared_modulus = std::pow(error_x,2)+std::pow(error_y,2)+std::pow(error_z,2);
        // KRATOS_WATCH(squared_modulus)
        result += squared_modulus * local_area;

    }
    return std::sqrt(result/total_area);
}

double GetL2ScalarErrorNormCD(ModelPart& r_model_part)
{
    double total_area = 0.0, result = 0.0;

    for (const auto& r_node : r_model_part.Nodes())
    {
        double local_area = r_node.FastGetSolutionStepValue(NODAL_AREA);
        total_area += local_area;
        double error_p = r_node.FastGetSolutionStepValue(SCALAR_STAT_ERROR);
        result += std::pow(error_p,2) * local_area;
        // KRATOS_WATCH(result)
        // KRATOS_WATCH(total_area)
    }

    return std::sqrt(result/total_area);
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

}; // Class CalculateStationarityProcess

} // namespace Kratos.

#endif // KRATOS_STATIONARITY_ERROR_CALCULATOR_UTILITY  defined