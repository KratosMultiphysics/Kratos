//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
//

#ifndef NLOPT_OPTIMIZER_H
#define NLOPT_OPTIMIZER_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>
#include <functional>


// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "algorithm_base.h"

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------
#include "custom_external_libraries/header/nlopt.hpp"

// ==============================================================================

namespace Kratos {

typedef struct {
    double a, b;
} my_constraint_data;

double myvconstraint(const std::vector<double> &x, std::vector<double> &grad, void *data)
{
  my_constraint_data *d = reinterpret_cast<my_constraint_data*>(data);
  double a = d->a, b = d->b;
  if (!grad.empty()) {
    grad[0] = 3 * a * (a*x[0] + b) * (a*x[0] + b);
    grad[1] = -1.0;
  }
  return ((a*x[0] + b) * (a*x[0] + b) * (a*x[0] + b) - x[1]);
}




class KRATOS_API(OPTIMIZATION_APPLICATION) NLOptOptimizer 
{
    public:
        KRATOS_CLASS_POINTER_DEFINITION(NLOptOptimizer);

      void Set_Function(std::function <double(std::vector<double> x, std::vector<double> grad, void *data)> myfunc );
      int Optimize();

    private:

    std ::function <const double(std::vector<double> &x, std::vector<double> &grad, void *data)> callback_;

       
    }; // Class NLOptOptimizer

} // namespace Kratos

#endif // NLOPT_OPTIMIZER_H
