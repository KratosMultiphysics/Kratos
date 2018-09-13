#include "velocity_field.h"

namespace Kratos
{
void VelocityField::Evaluate(const double time,
                             const array_1d<double, 3>& coor,
                             array_1d<double, 3>& vector,
                             const int i_thread)
{
    UpdateCoordinates(time, coor, i_thread);
    vector[0] = U0(i_thread);
    vector[1] = U1(i_thread);
    vector[2] = U2(i_thread);
}

void VelocityField::Evaluate(const double time,
                             const DenseVector<double>& coor,
                             DenseVector<double>& result,
                             const int i_thread)
{
    UpdateCoordinates(time, coor, i_thread);
    result[0] = U0(i_thread);
    result[1] = U1(i_thread);
    result[2] = U2(i_thread);
}

void VelocityField::CalculateTimeDerivative(const double time,
                                            const array_1d<double, 3>& coor,
                                            array_1d<double, 3>& deriv,
                                            const int i_thread)
{
    UpdateCoordinates(time, coor, i_thread);
    deriv[0] = U0DT(i_thread);
    deriv[1] = U1DT(i_thread);
    deriv[2] = U2DT(i_thread);
}

void VelocityField::CalculateTimeDerivative(const double time,
                                            const DenseVector<double>& coor,
                                            DenseVector<double>& deriv,
                                            const int i_thread)
{
    UpdateCoordinates(time, coor, i_thread);
    deriv[0] = U0DT(i_thread);
    deriv[1] = U1DT(i_thread);
    deriv[2] = U2DT(i_thread);
}

void VelocityField::CalculateGradient(const double time,
                                      const array_1d<double, 3>& coor,
                                      array_1d< array_1d<double, 3>, 3>& gradient,
                                      const int i_thread)
{
    UpdateCoordinates(time, coor, i_thread);
    gradient[0][0] = U0D0(i_thread);
    gradient[0][1] = U0D1(i_thread);
    gradient[0][2] = U0D2(i_thread);
    gradient[1][0] = U1D0(i_thread);
    gradient[1][1] = U1D1(i_thread);
    gradient[1][2] = U1D2(i_thread);
    gradient[2][0] = U2D0(i_thread);
    gradient[2][1] = U2D1(i_thread);
    gradient[2][2] = U2D2(i_thread);
}

void VelocityField::CalculateGradient(const double time,
                                      const array_1d<double, 3>& coor,
                                      DenseVector< double>& gradient_x,
                                      DenseVector< double>& gradient_y,
                                      DenseVector< double>& gradient_z,
                                      const int i_thread)
{
    UpdateCoordinates(time, coor, i_thread);
    gradient_x[0] = U0D0(i_thread);
    gradient_x[1] = U0D1(i_thread);
    gradient_x[2] = U0D2(i_thread);
    gradient_y[0] = U1D0(i_thread);
    gradient_y[1] = U1D1(i_thread);
    gradient_y[2] = U1D2(i_thread);
    gradient_z[0] = U2D0(i_thread);
    gradient_z[1] = U2D1(i_thread);
    gradient_z[2] = U2D2(i_thread);
}

double VelocityField::CalculateDivergence(const double time, const array_1d<double, 3>& coor, const int i_thread)
{
    UpdateCoordinates(time, coor, i_thread);
    double div = U0D0(i_thread) + U1D1(i_thread) + U2D2(i_thread);
    return div;
}

double VelocityField::CalculateDivergence(const double time, const DenseVector<double>& coor, const int i_thread)
{
    UpdateCoordinates(time, coor, i_thread);
    double div = U0D0(i_thread) + U1D1(i_thread) + U2D2(i_thread);
    return div;
}

void VelocityField::CalculateRotational(const double time,
                                        const array_1d<double, 3>& coor,
                                        array_1d<double, 3>& rot,
                                        const int i_thread)
{
    UpdateCoordinates(time, coor, i_thread);
    rot[0] = U2D1(i_thread) - U1D2(i_thread);
    rot[1] = U0D2(i_thread) - U2D0(i_thread);
    rot[2] = U1D0(i_thread) - U0D1(i_thread);
}

void VelocityField::CalculateRotational(const double time,
                                        const DenseVector<double>& coor,
                                        DenseVector<double>& rot,
                                        const int i_thread)
{
    UpdateCoordinates(time, coor, i_thread);
    rot[0] = U2D1(i_thread) - U1D2(i_thread);
    rot[1] = U0D2(i_thread) - U2D0(i_thread);
    rot[2] = U1D0(i_thread) - U0D1(i_thread);
}

void VelocityField::CalculateLaplacian(const double time,
                                       const array_1d<double, 3>& coor,
                                       array_1d<double, 3>& lapl,
                                       const int i_thread)

{
    UpdateCoordinates(time, coor, i_thread);
    lapl[0] = U0D0D0(i_thread) + U0D1D1(i_thread) + U0D2D2(i_thread);
    lapl[1] = U1D0D0(i_thread) + U1D1D1(i_thread) + U1D2D2(i_thread);
    lapl[2] = U2D0D0(i_thread) + U2D1D1(i_thread) + U2D2D2(i_thread);
}

void VelocityField::CalculateLaplacian(const double time,
                                       const DenseVector<double>& coor,
                                       DenseVector<double>& lapl,
                                       const int i_thread)
{
    UpdateCoordinates(time, coor, i_thread);
    lapl[0] = U0D0D0(i_thread) + U0D1D1(i_thread) + U0D2D2(i_thread);
    lapl[1] = U1D0D0(i_thread) + U1D1D1(i_thread) + U1D2D2(i_thread);
    lapl[2] = U2D0D0(i_thread) + U2D1D1(i_thread) + U2D2D2(i_thread);
}

void VelocityField::CalculateMaterialAcceleration(const double time,
                                                  const array_1d<double, 3>& coor,
                                                  array_1d<double, 3>& accel,
                                                  const int i_thread)
{
    UpdateCoordinates(time, coor, i_thread);
    double u0 = U0(i_thread);
    double u1 = U1(i_thread);
    double u2 = U2(i_thread);
    array_1d<double, 3> vel_rate;
    array_1d< array_1d<double, 3>, 3> grad;
    CalculateTimeDerivative(time, coor, vel_rate, i_thread);
    CalculateGradient(time, coor, grad, i_thread);

    accel[0] = vel_rate[0] + u0 * grad[0][0] + u1 * grad[0][1] + u2 * grad[0][2];
    accel[1] = vel_rate[1] + u0 * grad[1][0] + u1 * grad[1][1] + u2 * grad[1][2];
    accel[2] = vel_rate[2] + u0 * grad[2][0] + u1 * grad[2][1] + u2 * grad[2][2];
}

void VelocityField::CalculateMaterialAcceleration(const double time,
                                                  const DenseVector<double>& coor,
                                                  DenseVector<double>& accel,
                                                  const int i_thread)
{
    UpdateCoordinates(time, coor, i_thread);
    double u0 = U0(i_thread);
    double u1 = U1(i_thread);
    double u2 = U2(i_thread);
    array_1d<double, 3> vel_rate;
    array_1d< array_1d<double, 3>, 3> grad;
    CalculateTimeDerivative(time, coor, vel_rate, i_thread);
    CalculateGradient(time, coor, grad, i_thread);

    accel[0] = vel_rate[0] + u0 * grad[0][0] + u1 * grad[0][1] + u2 * grad[0][2];
    accel[1] = vel_rate[1] + u0 * grad[1][0] + u1 * grad[1][1] + u2 * grad[1][2];
    accel[2] = vel_rate[2] + u0 * grad[2][0] + u1 * grad[2][1] + u2 * grad[2][2];
}

void VelocityField::CalculateConvectiveDerivative(const double time,
                                                  const array_1d<double, 3>& coor,
                                                  array_1d<double, 3>& accel,
                                                  const int i_thread)
{
    UpdateCoordinates(time, coor, i_thread);
    double u0 = U0(i_thread);
    double u1 = U1(i_thread);
    double u2 = U2(i_thread);
    array_1d< array_1d<double, 3>, 3> grad;
    CalculateGradient(time, coor, grad, i_thread);

    accel[0] = u0 * grad[0][0] + u1 * grad[0][1] + u2 * grad[0][2];
    accel[1] = u0 * grad[1][0] + u1 * grad[1][1] + u2 * grad[1][2];
    accel[2] = u0 * grad[2][0] + u1 * grad[2][1] + u2 * grad[2][2];
}

void VelocityField::CalculateConvectiveDerivative(const double time,
                                                  const DenseVector<double>& coor,
                                                  DenseVector<double>& accel,
                                                  const int i_thread)
{
    UpdateCoordinates(time, coor, i_thread);
    double u0 = U0(i_thread);
    double u1 = U1(i_thread);
    double u2 = U2(i_thread);
    array_1d< array_1d<double, 3>, 3> grad;
    CalculateGradient(time, coor, grad, i_thread);

    accel[0] = u0 * grad[0][0] + u1 * grad[0][1] + u2 * grad[0][2];
    accel[1] = u0 * grad[1][0] + u1 * grad[1][1] + u2 * grad[1][2];
    accel[2] = u0 * grad[2][0] + u1 * grad[2][1] + u2 * grad[2][2];
}

void VelocityField::CalculateAccelerationFollowingTheParticle(const double time,
                                                              const array_1d<double, 3>& coor,
                                                              array_1d<double, 3>& accel,
                                                              const array_1d<double, 3>& particle_vel,
                                                              const int i_thread)
{
    UpdateCoordinates(time, coor, i_thread);
    double u0 = particle_vel[0];
    double u1 = particle_vel[1];
    double u2 = particle_vel[2];
    array_1d<double, 3> vel_rate;
    array_1d< array_1d<double, 3>, 3> grad;
    CalculateTimeDerivative(time, coor, vel_rate, i_thread);
    CalculateGradient(time, coor, grad, i_thread);

    accel[0] = vel_rate[0] + u0 * grad[0][0] + u1 * grad[0][1] + u2 * grad[0][2];
    accel[1] = vel_rate[1] + u0 * grad[1][0] + u1 * grad[1][1] + u2 * grad[1][2];
    accel[2] = vel_rate[2] + u0 * grad[2][0] + u1 * grad[2][1] + u2 * grad[2][2];
}

void VelocityField::ImposeFieldOnNodes(ModelPart& r_model_part, const VariablesList& variables_to_be_imposed)
{
    bool must_impose_fluid_velocity                            = VariableIsInList(variables_to_be_imposed, FLUID_VEL_PROJECTED);
    bool must_impose_fluid_acceleration                        = VariableIsInList(variables_to_be_imposed, FLUID_ACCEL_PROJECTED);
    bool must_impose_fluid_acceleration_following_the_particle = VariableIsInList(variables_to_be_imposed, FLUID_ACCEL_FOLLOWING_PARTICLE_PROJECTED);
    bool must_impose_fluid_velocity_laplacian                  = VariableIsInList(variables_to_be_imposed, FLUID_VEL_LAPL_PROJECTED);
    const double time = r_model_part.GetProcessInfo()[TIME];

    #pragma omp parallel for firstprivate(must_impose_fluid_velocity, must_impose_fluid_acceleration, must_impose_fluid_velocity_laplacian, time)
    for (int i = 0; i < (int)r_model_part.Nodes().size(); i++){
        int thread_number = OpenMPUtils::ThisThread();
        ModelPart::NodesContainerType::iterator i_particle = r_model_part.NodesBegin() + i;
        Node<3>::Pointer p_node = *(i_particle.base());
        const array_1d<double, 3>& coor = p_node->Coordinates();
        UpdateCoordinates(time, coor, thread_number);
        LockCoordinates(thread_number);

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

        if (must_impose_fluid_acceleration_following_the_particle){
            const array_1d<double, 3> particle_vel = p_node->FastGetSolutionStepValue(VELOCITY);
            array_1d<double, 3> fluid_accel;
            CalculateAccelerationFollowingTheParticle(time, coor, fluid_accel, particle_vel, thread_number);
            array_1d<double, 3>& fluid_accel_projected = p_node->FastGetSolutionStepValue(FLUID_ACCEL_FOLLOWING_PARTICLE_PROJECTED);
            noalias(fluid_accel_projected) = fluid_accel;
        }

        if (must_impose_fluid_velocity_laplacian){
            array_1d<double, 3> fluid_laplacian;
            CalculateLaplacian(time, coor, fluid_laplacian, thread_number);
            array_1d<double, 3>& fluid_laplacian_projected = p_node->FastGetSolutionStepValue(FLUID_VEL_LAPL_PROJECTED);
            noalias(fluid_laplacian_projected) = fluid_laplacian;
        }

        UnlockCoordinates(thread_number);
    }
}

void VelocityField::ImposeVelocityOnNodes(ModelPart& r_model_part, const Variable<array_1d<double, 3> >& container_variable)
{
    const double time = r_model_part.GetProcessInfo()[TIME];
    int thread_number = OpenMPUtils::ThisThread();

    for (int i = 0; i < (int)r_model_part.Nodes().size(); ++i){
        ModelPart::NodesContainerType::iterator i_particle = r_model_part.NodesBegin() + i;
        Node<3>::Pointer p_node = *(i_particle.base());
        const array_1d<double, 3>& coor = p_node->Coordinates();
        array_1d<double, 3> fluid_vel;
        Evaluate(time, coor, fluid_vel, thread_number);
        array_1d<double, 3>& slip_vel = p_node->FastGetSolutionStepValue(container_variable);
        noalias(slip_vel) = fluid_vel;
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
