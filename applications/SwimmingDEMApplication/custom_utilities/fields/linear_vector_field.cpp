#include "linear_vector_field.h"

namespace Kratos
{

void LinearVectorField::Evaluate(const double time,
                                 const array_1d<double, 3>& coor,
                                 array_1d<double, 3>& vector,
                                 const int i_thread)
{
    noalias(vector) = prod(mA, coor) + mb;
}

void  LinearVectorField::Evaluate(const double time,
              const DenseVector<double>& coor,
              DenseVector<double>& vector,
              const int i_thread)
{
    vector = prod(mA, coor) + mb;
}

void LinearVectorField::CalculateGradient(const double time,
                                          const array_1d<double, 3>& coor,
                                          array_1d< array_1d<double, 3>, 3>& gradient,
                                          const int i_thread)
{

}

void LinearVectorField::CalculateGradient(const double time,
                                          const DenseVector<double>& coor,
                                          DenseMatrix<double>& gradient,
                                          const int i_thread)
{
    for (auto i = 0; i < 3; ++i){
        for (auto j = 0; j < 3; ++j){
            gradient(i, j) = mA(i, j);
        }
    }
}

double LinearVectorField::CalculateDivergence(const double time,
                                              const array_1d<double, 3>& coor,
                                              const int i_thread)
{
    return mA(0, 0) + mA(1, 1) + mA(2, 2);
}

double LinearVectorField::CalculateDivergence(const double time,
                                              const DenseVector<double>& coor,
                                              const int i_thread)
{
    return mA(0, 0) + mA(1, 1) + mA(2, 2);
}

void LinearVectorField::CalculateRotational(const double time,
                                            const array_1d<double, 3>& coor,
                                            array_1d<double, 3>& rot,
                                            const int i_thread)
{
    rot[0] = mA(2, 1) - mA(1, 2);
    rot[1] = mA(0, 2) - mA(2, 0);
    rot[2] = mA(1, 0) - mA(0, 1);
}

void LinearVectorField::CalculateRotational(const double time,
                                            const DenseVector<double>& coor,
                                            DenseVector<double>& rot,
                                            const int i_thread)
{
    rot[0] = mA(2, 1) - mA(1, 2);
    rot[1] = mA(0, 2) - mA(2, 0);
    rot[2] = mA(1, 0) - mA(0, 1);
}

} // namespace Kratos.




