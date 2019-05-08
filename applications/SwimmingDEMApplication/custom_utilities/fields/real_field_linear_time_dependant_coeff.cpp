#include "real_field_linear_time_dependant_coeff.h"

namespace Kratos
{

// Values

double LinearRealField::Evaluate(const double time,
                                 const array_1d<double, 3>& coor,
                                 const int i_thread)
{
    return mX0 + mpFx->Evaluate(time) * coor[0] + mY0 + mpFy->Evaluate(time) * coor[1] + mZ0 + mpFz->Evaluate(time) * coor[2];
}

double LinearRealField::Evaluate(const double time,
                const DenseVector<double>& coor,
                                 const int i_thread)
{
    return mX0 + mpFx->Evaluate(time) * coor[0] + mY0 + mpFy->Evaluate(time) * coor[1] + mZ0 + mpFz->Evaluate(time) * coor[2];
}


double LinearRealField::CalculateTimeDerivative(const double time,
                                                const array_1d<double, 3>& coor,
                                                const int i_thread)
{
    return mpFx->CalculateDerivative(time) * coor[0] + mY0 + mpFy->CalculateDerivative(time) * coor[1] + mZ0 + mpFz->CalculateDerivative(time) * coor[2];;
}

double LinearRealField::CalculateTimeDerivative(const double time,
                               const DenseVector<double>& coor,
                                                const int i_thread)
{
    return mpFx->CalculateDerivative(time) * coor[0] + mY0 + mpFy->CalculateDerivative(time) * coor[1] + mZ0 + mpFz->CalculateDerivative(time) * coor[2];;
}

void LinearRealField::CalculateGradient(const double time,
                       const array_1d<double, 3>& coor,
                       array_1d<double, 3>& gradient,
                       const int i_thread)
{
    gradient[0] = mpFx->Evaluate(time);
    gradient[1] = mpFy->Evaluate(time);
    gradient[2] = mpFz->Evaluate(time);
}

void LinearRealField::CalculateGradient(const double time,
                       const DenseVector<double>& coor,
                       DenseVector<double>& gradient,
                       const int i_thread)
{
    gradient[0] = mpFx->Evaluate(time);
    gradient[1] = mpFy->Evaluate(time);
    gradient[2] = mpFz->Evaluate(time);
}


// template<>
// double Evaluate(const double time,
//                 const array_1d<double, 3>& coor,
//                 const int i_thread)
// {

// }

} // namespace Kratos.




