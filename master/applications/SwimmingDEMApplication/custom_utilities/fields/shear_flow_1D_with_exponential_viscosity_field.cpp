#include "shear_flow_1D_with_exponential_viscosity_field.h"

namespace Kratos
{

void ShearFlow1DWithExponentialViscosityField::Evaluate(const double time, const array_1d<double, 3>& coor, array_1d<double, 3>& vector, const int i_thread)
{
//    if (coor[2] <= 0.0){
//        vector[0] = 0.0;
//        vector[1] = 0.0;
//        vector[2] = 0.0;
//    }
//    else {
//        if (fabs(mZmax) < 1e-7){
//            vector[0] = mUFarField;
//        }
//        else {
//            vector[0] = mUFarField * std::min(1.0, fabs(coor[2] / mZmax));
//            vector[1] = 0;
//            vector[2] = 0;
//        }
//    }

    vector[0] = mUFarField;
    vector[1] = 0.0;
    vector[2] = 0.0;
}

void ShearFlow1DWithExponentialViscosityField::SetRimZoneThickness(const double z_max)
{
    mZmax = z_max;
}

} // namespace Kratos.




