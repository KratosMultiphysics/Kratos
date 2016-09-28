#include "velocity_field.h"

namespace Kratos
{
void VelocityField::Evaluate(const double time, const array_1d<double, 3>& coor, array_1d<double, 3>& vector, const unsigned int i)
{
    UpdateCoordinates(time, coor, i);
    vector[0] = U0(i);
    vector[1] = U1(i);
    vector[2] = U2(i);
}

void VelocityField::Evaluate(const double time, const vector<double>& coor, vector<double>& result, const unsigned int i)
{
    UpdateCoordinates(time, coor, i);
    result[0] = U0(i);
    result[1] = U1(i);
    result[2] = U2(i);
}

void VelocityField::CalculateTimeDerivative(const double time, const array_1d<double, 3>& coor, array_1d<double, 3>& deriv, const unsigned int i)
{
    UpdateCoordinates(time, coor, i);
    deriv[0] = U0DT(i);
    deriv[1] = U1DT(i);
    deriv[2] = U2DT(i);
}

void VelocityField::CalculateTimeDerivative(const double time, const vector<double>& coor, vector<double>& deriv, const unsigned int i)
{
    UpdateCoordinates(time, coor, i);
    deriv[0] = U0DT(i);
    deriv[1] = U1DT(i);
    deriv[2] = U2DT(i);
}


void VelocityField::CalculateGradient(const double time, const array_1d<double, 3>& coor, array_1d< array_1d<double, 3>, 3>& gradient, const unsigned int i)
{
    UpdateCoordinates(time, coor, i);
    gradient[0][0] = U0D0(i);
    gradient[0][1] = U0D1(i);
    gradient[0][2] = U0D2(i);
    gradient[1][0] = U1D0(i);
    gradient[1][1] = U1D1(i);
    gradient[1][2] = U1D2(i);
    gradient[2][0] = U2D0(i);
    gradient[2][1] = U2D1(i);
    gradient[2][2] = U2D2(i);
}

void VelocityField::CalculateDivergence(const double time, const array_1d<double, 3>& coor, double& div, const unsigned int i)
{
    UpdateCoordinates(time, coor, i);
    div = U0D0(i) + U1D1(i) + U2D2(i);
}

double VelocityField::CalculateDivergence(const double time, const vector<double>& coor, const unsigned int i)
{
    UpdateCoordinates(time, coor, i);
    double div = U0D0(i) + U1D1(i) + U2D2(i);
    return div;
}

void VelocityField::CalculateRotational(const double time, const array_1d<double, 3>& coor, array_1d<double, 3>& rot, const unsigned int i)
{
    UpdateCoordinates(time, coor, i);
    rot[0] = U2D1(i) - U1D2(i);
    rot[1] = U0D2(i) - U2D0(i);
    rot[2] = U1D0(i) - U0D1(i);
}

void VelocityField::CalculateRotational(const double time, const vector<double>& coor, vector<double>& rot, const unsigned int i)
{
    UpdateCoordinates(time, coor, i);
    rot[0] = U2D1(i) - U1D2(i);
    rot[1] = U0D2(i) - U2D0(i);
    rot[2] = U1D0(i) - U0D1(i);
}

void VelocityField::CalculateLaplacian(const double time, const array_1d<double, 3>& coor, array_1d<double, 3>& lapl, const unsigned int i)
{
    UpdateCoordinates(time, coor, i);
    lapl[0] = U0D0D0(i) + U0D1D1(i) + U0D2D2(i);
    lapl[1] = U1D0D0(i) + U1D1D1(i) + U1D2D2(i);
    lapl[2] = U2D0D0(i) + U2D1D1(i) + U2D2D2(i);
}

void VelocityField::CalculateLaplacian(const double time, const vector<double>& coor, vector<double>& lapl, const unsigned int i)
{
    UpdateCoordinates(time, coor, i);
    lapl[0] = U0D0D0(i) + U0D1D1(i) + U0D2D2(i);
    lapl[1] = U1D0D0(i) + U1D1D1(i) + U1D2D2(i);
    lapl[2] = U2D0D0(i) + U2D1D1(i) + U2D2D2(i);
}

void VelocityField::CalculateMaterialAcceleration(const double time, const array_1d<double, 3>& coor, array_1d<double, 3>& accel, const unsigned int i)
{
    UpdateCoordinates(time, coor, i);
    double u0 = U0(i);
    double u1 = U1(i);
    double u2 = U2(i);
    array_1d<double, 3> vel_rate;
    array_1d< array_1d<double, 3>, 3> grad;
    CalculateTimeDerivative(time, coor, vel_rate);
    CalculateGradient(time, coor, grad);

    accel[0] = vel_rate[0] + u0 * grad[0][0] + u1 * grad[0][1] + u2 * grad[0][2];
    accel[1] = vel_rate[1] + u0 * grad[1][0] + u1 * grad[1][1] + u2 * grad[1][2];
    accel[2] = vel_rate[2] + u0 * grad[2][0] + u1 * grad[2][1] + u2 * grad[2][2];
}

void VelocityField::CalculateMaterialAcceleration(const double time, const vector<double>& coor, vector<double>& accel, const unsigned int i)
{
    UpdateCoordinates(time, coor, i);
    double u0 = U0(i);
    double u1 = U1(i);
    double u2 = U2(i);
    array_1d<double, 3> vel_rate;
    array_1d< array_1d<double, 3>, 3> grad;
    CalculateTimeDerivative(time, coor, vel_rate);
    CalculateGradient(time, coor, grad);

    accel[0] = vel_rate[0] + u0 * grad[0][0] + u1 * grad[0][1] + u2 * grad[0][2];
    accel[1] = vel_rate[1] + u0 * grad[1][0] + u1 * grad[1][1] + u2 * grad[1][2];
    accel[2] = vel_rate[2] + u0 * grad[2][0] + u1 * grad[2][1] + u2 * grad[2][2];
}

void VelocityField::ImposeFieldOnNodes(ModelPart& r_model_part, const VariablesList& variables_to_be_imposed)
{
    bool must_impose_fluid_velocity           = VariableIsInList(variables_to_be_imposed, FLUID_VEL_PROJECTED);
    bool must_impose_fluid_acceleration       = VariableIsInList(variables_to_be_imposed, FLUID_ACCEL_PROJECTED);
    bool must_impose_fluid_velocity_laplacian = VariableIsInList(variables_to_be_imposed, FLUID_VEL_LAPL_PROJECTED);
    const double time = r_model_part.GetProcessInfo()[TIME];

    #pragma omp parallel firstprivate(must_impose_fluid_velocity, must_impose_fluid_acceleration, must_impose_fluid_velocity_laplacian, time)
    {

        #pragma omp barrier

        unsigned int thread_number = omp_get_thread_num();

        #pragma omp for
        for (int i = 0; i < (int)r_model_part.Nodes().size(); i++){
            ModelPart::NodesContainerType::iterator i_particle = r_model_part.NodesBegin() + i;
            Node<3>::Pointer p_node = *(i_particle.base());
            const array_1d<double, 3>& coor = p_node->Coordinates();

            if (must_impose_fluid_velocity){
                array_1d<double, 3> fluid_vel;
                Evaluate(time, coor, fluid_vel, thread_number);
                array_1d<double, 3>& fluid_vel_projected = p_node->FastGetSolutionStepValue(FLUID_VEL_PROJECTED);
                noalias(fluid_vel_projected) = fluid_vel;
            }

            if (must_impose_fluid_acceleration){
                array_1d<double, 3> fluid_accel;
                CalculateMaterialAcceleration(time, coor, fluid_accel, thread_number);
                array_1d<double, 3>& fluid_accel_projected = p_node->FastGetSolutionStepValue(FLUID_ACCEL_PROJECTED);
                noalias(fluid_accel_projected) = fluid_accel;
            }

            if (must_impose_fluid_velocity_laplacian){
                array_1d<double, 3> fluid_laplacian;
                CalculateLaplacian(time, coor, fluid_laplacian, thread_number);
                array_1d<double, 3>& fluid_laplacian_projected = p_node->FastGetSolutionStepValue(FLUID_VEL_LAPL_PROJECTED);
                noalias(fluid_laplacian_projected) = fluid_laplacian;
            }
        }
    }
}

bool VelocityField::VariableIsInList(const VariablesList var_list, const VariableData& var)
{
    for (unsigned int i = 0; i != var_list.size(); ++i){

        if (*var_list[i] == var){
            return true;
        }
    }

    return false;
}

} // namespace Kratos.

