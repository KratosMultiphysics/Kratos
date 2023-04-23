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

#include "nlopt_optimizer.h"
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include "custom_external_libraries/header/nlopt.hpp"

namespace Kratos {
 

void NLOptOptimizer::Set_Function(std:: function <double( std::vector<double> x, std::vector<double> grad, void *data)> myfunc ){
   callback_ = std::move(myfunc);
}


int NLOptOptimizer::Optimize(){ 
    std :: vector <double> v1; 
    v1.push_back(1);
    v1.push_back(1);                     
    callback_(v1,v1,NULL);
    nlopt::opt opt("LD_MMA", 2);
    std::vector<double> lb(2);
    lb[0] = -HUGE_VAL; lb[1] = 0;
    opt.set_lower_bounds(lb);
    opt.set_min_objective(callback_, NULL);
    my_constraint_data data[2] = { {2,0}, {-1,1} };
    opt.add_inequality_constraint(myvconstraint, &data[0], 1e-8);
    opt.add_inequality_constraint(myvconstraint, &data[1], 1e-8);
    opt.set_xtol_rel(1e-4);

  // try setting an algorithm parameter: */
  opt.set_param("inner_maxeval", 123);
  if (opt.get_param("inner_maxeval", 1234) != 123 || opt.get_param("not a param", 1234) != 1234 ||
      opt.num_params() != 1 || std::string(opt.nth_param(0)) != "inner_maxeval") {
    std::cerr << "failed to retrieve nlopt parameter" << std::endl;
    return EXIT_FAILURE;
  }

  std::vector<double> x(2);
  x[0] = 1.234; x[1] = 5.678;
  double minf;

  try{
    opt.optimize(x, minf);
    std::cerr << "found minimum at f(" << x[0] << "," << x[1] << ") = "
              << std::setprecision(10) << minf <<std::endl;
    return std::fabs(minf - 0.5443310474) < 1e-3 ? EXIT_SUCCESS : EXIT_FAILURE;
  }
  catch(std::exception &e) {
    std::cerr << "nlopt failed: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  }
}

