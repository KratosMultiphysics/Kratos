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

class ErrorNormCalculator
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
    typedef Element::Pointer  ElementPointerType;
    typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionDerivativesArrayType;

    KRATOS_CLASS_POINTER_DEFINITION(ErrorNormCalculator);

ErrorNormCalculator()
{}

virtual ~ErrorNormCalculator(){}

void ComputeDofsErrors(ModelPart& r_model_part);

double GetL2VectorErrorNorm(ModelPart& r_model_part, const Variable<array_1d<double,3>>& rVariable, const Variable<array_1d<double,3>>& rExact_Variable);

double GetL2ScalarErrorNorm(ModelPart& r_model_part, const Variable<double>& rVariable, const Variable<double>& rExact_Variable);

double GetH1ScalarErrorSemiNorm(ModelPart& r_model_part);

double GetH1VectorErrorSemiNorm(ModelPart& r_model_part);

double GetL2MaterialAccelNormUsingNodeDerivatives(ModelPart& r_model_part, const Variable<array_1d<double,3>>& rExact_Variable);

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

}; // Class L2ErrorNormCalculator

} // namespace Kratos.

#endif // KRATOS_L2_ERROR_CALCULATOR_UTILITY  defined