#include "velocity_field.h"

namespace Kratos
{
void VelocityField::Evaluate(const double time, const array_1d<double, 3>& coor, array_1d<double, 3>& vector)
{
    UpdateCoordinates(time, coor);
    vector[0] = U0();
    vector[1] = U1();
    vector[2] = U2();
} // namespace Kratos.

void VelocityField::CalculateTimeDerivative(const double time, const array_1d<double, 3>& coor, array_1d<double, 3>& deriv)
{
    UpdateCoordinates(time, coor);
    deriv[0] = U0DT();
    deriv[1] = U1DT();
    deriv[2] = U2DT();
}

void VelocityField::CalculateGradient(const double time, const array_1d<double, 3>& coor, array_1d< array_1d<double, 3>, 3>& gradient)
{
    UpdateCoordinates(time, coor);
    gradient[0][0] = U0D0();
    gradient[0][1] = U0D1();
    gradient[0][2] = U0D2();
    gradient[1][0] = U1D0();
    gradient[1][1] = U1D1();
    gradient[1][2] = U1D2();
    gradient[2][0] = U2D0();
    gradient[2][1] = U2D1();
    gradient[2][2] = U2D2();
}

void VelocityField::CalculateDivergence(const double time, const array_1d<double, 3>& coor, double& div)
{
    UpdateCoordinates(time, coor);
    div = U0D0() + U1D1() + U2D2();
}

void VelocityField::CalculateRotational(const double time, const array_1d<double, 3>& coor, array_1d<double, 3>& rot)
{
    UpdateCoordinates(time, coor);
    rot[0] = U2D1() - U1D2();
    rot[1] = U0D2() - U2D0();
    rot[2] = U1D0() - U0D1();
}

void VelocityField::CalculateLaplacian(const double time, const array_1d<double, 3>& coor, array_1d<double, 3>& lapl)
{
    UpdateCoordinates(time, coor);
    lapl[0] = U0D0D0() + U0D1D1() + U0D2D2();
    lapl[1] = U1D0D0() + U1D1D1() + U1D2D2();
    lapl[2] = U2D0D0() + U2D1D1() + U2D2D2();
}

void VelocityField::CalculateMaterialAcceleration(const double time, const array_1d<double, 3>& coor, array_1d<double, 3>& accel)
{
    UpdateCoordinates(time, coor);
    double u0 = U0();
    double u1 = U1();
    double u2 = U2();
    array_1d<double, 3> vel_rate;
    array_1d< array_1d<double, 3>, 3> grad;
    CalculateTimeDerivative(time, coor, vel_rate);
    CalculateGradient(time, coor, grad);

    accel[0] = vel_rate[0] + u0 * grad[0][0] + u1 * grad[0][1] + u2 * grad[0][2];
    accel[1] = vel_rate[1] + u0 * grad[1][0] + u1 * grad[1][1] + u2 * grad[1][2];
    accel[2] = vel_rate[2] + u0 * grad[2][0] + u1 * grad[2][1] + u2 * grad[2][2];
}
} // namespace Kratos.




