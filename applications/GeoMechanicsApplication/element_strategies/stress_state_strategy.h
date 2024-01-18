// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#pragma once

class Matrix;
class Vector;

namespace Kratos
{

class StressStateStrategy
{
public:
    virtual void CalculateBMatrix(Matrix& rB, const Matrix& GradNpT, const Vector& Np) = 0;
    virtual ~StressStateStrategy()                                                     = default;
};

} // namespace Kratos
