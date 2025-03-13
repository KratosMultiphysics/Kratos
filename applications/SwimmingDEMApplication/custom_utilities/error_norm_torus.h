//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Aniol Sala
//

#if !defined(KRATOS_ERROR_TORUS)
#define KRATOS_ERROR_TORUS

// External includes
#include <omp.h>
#include <vector>

// Project includes
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"
#include "includes/kratos_parameters.h"
#include "includes/process_info.h"

#include "swimming_DEM_application.h"

namespace Kratos
{

class ErrorNormTorusCalculator
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

    KRATOS_CLASS_POINTER_DEFINITION(ErrorNormTorusCalculator);

ErrorNormTorusCalculator(const double major_radius, const double minor_radius, const double center_velocity): mMajorRadius(major_radius), mMinorRadius(minor_radius), mCenterVelocity(center_velocity)
{}

virtual ~ErrorNormTorusCalculator(){}

double getDistanceToCenter(const array_1d<double, 3>& coor);
double getVelocityModule(const array_1d<double, 3>& coor);
double getThetaAngle(const array_1d<double, 3>& coor);
void CalculateVelocity(const array_1d<double, 3>& coor, array_1d<double, 3>& vel);
void CalculateMaterialAcceleration(const array_1d<double, 3>& coor, array_1d<double, 3>& accel);

double getMaterialAccelerationL2NormExactGaussPointsValues(ModelPart& r_model_part, const Variable<array_1d<double,3>>& rExact_Variable);
double getMaterialAccelerationL2NormUsingShapeFunctionsExactNodesValues(ModelPart& r_model_part, const Variable<array_1d<double,3>>& rExact_Variable);
double getMaterialAccelerationL2NormUsingShapeFunctionsExactGaussPointsValues(ModelPart& r_model_part);

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

private:
const double mMajorRadius;
const double mMinorRadius;
const double mCenterVelocity;

}; // Class ErrorNormTorusCalculator

} // namespace Kratos.

#endif // KRATOS_ERROR_TORUS defined