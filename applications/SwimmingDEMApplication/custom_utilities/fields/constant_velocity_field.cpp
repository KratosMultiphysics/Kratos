#include "constant_velocity_field.h"

namespace Kratos
{

void ConstantVelocityField::Evaluate(const double time, const array_1d<double, 3>& coor, array_1d<double, 3>& vector, const int i_thread)
{
    vector[0] = mVx;
    vector[1] = mVy;
    vector[2] = mVz;
}



} // namespace Kratos.




