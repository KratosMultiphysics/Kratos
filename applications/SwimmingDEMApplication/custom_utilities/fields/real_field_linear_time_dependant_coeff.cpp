#include "real_field_linear_time_dependant_coeff.h"

namespace Kratos
{

// Values

double LinearRealField::Evaluate(const double time,
                                 const array_1d<double, 3>& coor,
                                 const int i_thread)
{
    return mpB->Evaluate(time) + mA[0]->Evaluate(time) * coor[0] + mA[1]->Evaluate(time) * coor[1] + mA[2]->Evaluate(time) * coor[2];
}

double LinearRealField::Evaluate(const double time,
                const DenseVector<double>& coor,
                                 const int i_thread)
{
    return mpB->Evaluate(time) + mA[0]->Evaluate(time) * coor[0] + mA[1]->Evaluate(time) * coor[1] + mA[2]->Evaluate(time) * coor[2];
}


double LinearRealField::CalculateTimeDerivative(const double time,
                                                const array_1d<double, 3>& coor,
                                                const int i_thread)
{
    return mA[0]->CalculateDerivative(time) * coor[0] + mA[1]->CalculateDerivative(time) * coor[1] + mA[2]->CalculateDerivative(time) * coor[2];;
}

double LinearRealField::CalculateTimeDerivative(const double time,
                               const DenseVector<double>& coor,
                                                const int i_thread)
{
    return mA[0]->CalculateDerivative(time) * coor[0] + mA[1]->CalculateDerivative(time) * coor[1] + mA[2]->CalculateDerivative(time) * coor[2];;
}

void LinearRealField::CalculateGradient(const double time,
                       const array_1d<double, 3>& coor,
                       array_1d<double, 3>& gradient,
                       const int i_thread)
{
    gradient[0] = mA[0]->Evaluate(time);
    gradient[1] = mA[1]->Evaluate(time);
    gradient[2] = mA[2]->Evaluate(time);
}

void LinearRealField::CalculateGradient(const double time,
                       const DenseVector<double>& coor,
                       DenseVector<double>& gradient,
                       const int i_thread)
{
    gradient[0] = mA[0]->Evaluate(time);
    gradient[1] = mA[1]->Evaluate(time);
    gradient[2] = mA[2]->Evaluate(time);
}

// template<>
// double Evaluate(const double time,
//                 const array_1d<double, 3>& coor,
//                 const int i_thread)
// {

// }

} // namespace Kratos.




