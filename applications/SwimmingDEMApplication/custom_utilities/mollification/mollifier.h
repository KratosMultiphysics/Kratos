#ifndef KRATOS_MOLLIFIER_H
#define KRATOS_MOLLIFIER_H
// System includes
#include <string>
#include <iostream>
#include <stdlib.h>
#include "density_function.h"

namespace Kratos
{

template< unsigned int TDim>
class Mollifier: public DensityFunction<TDim>
{
public:
KRATOS_CLASS_POINTER_DEFINITION(Mollifier);

Mollifier():{}

virtual ~Mollifier(){}

void ComputeWeights(std::vector<double> & distances, std::vector<double> & weights){}

private:

}; // class DensityFunctionPolynomial

} // namespace Kratos.
#endif // KRATOS_MOLLIFIER_H
