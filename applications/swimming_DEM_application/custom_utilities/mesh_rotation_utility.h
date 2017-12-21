#if !defined(KRATOS_MESH_ROTATION_UTILITY)
#define KRATOS_MESH_ROTATION_UTILITY

// /* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

// Project includes
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"
#include "includes/kratos_parameters.h"

namespace Kratos
{
class MeshRotationUtility
{

public:

KRATOS_CLASS_POINTER_DEFINITION(MeshRotationUtility);

MeshRotationUtility(Parameters & r_parameters)
{
    mOmega = r_parameters["angular_velocity_magnitude"].GetDouble();
    mAInit = r_parameters["frame_rotation_axis_initial_point"].GetVector();
    mAFinal = r_parameters["frame_rotation_axis_final_point"].GetVector();
    mAxisVersor = CalculateNormalized(mAFinal - mAInit);

    // Constructing auxiliary Rodrigues matrices
    mI = identity_matrix<double>(3);

    for (int i = 0; i < 3; ++i){
        for (int j = 0; j < 3; ++j){
            mUU(i, j) = mAxisVersor[i] * mAxisVersor[j];
        }
    }

    mUx = ZeroMatrix(3, 3);
    mUx(0, 1) = - mAxisVersor[2];
    mUx(0, 2) =   mAxisVersor[1];
    mUx(1, 0) =   mAxisVersor[2];
    mUx(1, 2) = - mAxisVersor[0];
    mUx(2, 0) = - mAxisVersor[1];
    mUx(2, 1) =   mAxisVersor[0];
}

virtual ~MeshRotationUtility(){}

void RotateMesh(ModelPart& r_model_part, double time)
{
    CalculateRodriguesMatrices(time);

    const int nnodes = r_model_part.Nodes().size();

    if (nnodes > 0){
        auto it_begin = r_model_part.NodesBegin();
        array_1d<double, 3> P0;
        array_1d<double, 3> P;       

        #pragma omp parallel for private(P0, P)
        for (int i = 0; i < nnodes; ++i){
            auto it = it_begin + i;
            RotateNode(*(it.base()), P0, P);
            array_1d<double, 3>& displacement = it->FastGetSolutionStepValue(DISPLACEMENT);
            noalias(displacement) = P - P0;
            array_1d<double, 3>& velocity = it->FastGetSolutionStepValue(MESH_VELOCITY);
            noalias(velocity) = prod(mRp, P0 - mAInit);
        }
    }
}

void RotateDEMMesh(ModelPart& r_model_part, double time)
{
    CalculateRodriguesMatrices(time);
    const int nnodes = r_model_part.Nodes().size();
    if (nnodes > 0){
        auto it_begin = r_model_part.NodesBegin();
        array_1d<double, 3> P0;
        array_1d<double, 3> P;       

        #pragma omp parallel for private(P0, P)
        for (int i = 0; i < nnodes; ++i){
            auto it = it_begin + i;
            RotateNode(*(it.base()), P0, P);
            array_1d<double, 3>& displacement = it->FastGetSolutionStepValue(DISPLACEMENT);
            noalias(displacement) = P - P0;
            it->Fix(DISPLACEMENT_X);       
            it->Fix(DISPLACEMENT_Y); 
            it->Fix(DISPLACEMENT_Z);      
            array_1d<double, 3>& velocity = it->FastGetSolutionStepValue(VELOCITY);
            noalias(velocity) = prod(mRp, P0 - mAInit);
            it->Fix(VELOCITY_X);
            it->Fix(VELOCITY_Y);
            it->Fix(VELOCITY_Z);
        }
    }
}

private:

void CalculateRodriguesMatrices(const double time)
{
    double sin_theta = std::sin(time * mOmega);
    double cos_theta = std::cos(time * mOmega);

    // Rotation matrix
    noalias(mR) = cos_theta * mI + sin_theta * mUx + (1.0 - cos_theta) * mUU;

    // Rotation matrix derivative (derivative of R with respect to time)
    noalias(mRp) = - mOmega * sin_theta * mI + mOmega * cos_theta * mUx + mOmega * sin_theta * mUU;
}

void RotateNode(Kratos::Node<3>::Pointer p_node, array_1d<double, 3>& P0, array_1d<double, 3>& P)
{
    P0[0] = p_node->X0();
    P0[1] = p_node->Y0();
    P0[2] = p_node->Z0();
    noalias(P) = mAInit + prod(mR, P0 - mAInit);
    p_node->X() = P[0];
    p_node->Y() = P[1];
    p_node->Z() = P[2];
}

array_1d<double, 3> CalculateNormalized(array_1d<double, 3>&& vector)
{
    const double norm = norm_2(vector);

    if (norm > 0.0){
        vector *= 1.0 / norm;
    }

    return vector;
}

double mOmega;
array_1d<double, 3> mAInit;
array_1d<double, 3> mAFinal;
array_1d<double, 3> mAxisVersor;
boost::numeric::ublas::bounded_matrix<double, 3, 3> mI;
boost::numeric::ublas::bounded_matrix<double, 3, 3> mUU;
boost::numeric::ublas::bounded_matrix<double, 3, 3> mUx;
boost::numeric::ublas::bounded_matrix<double, 3, 3> mR;
boost::numeric::ublas::bounded_matrix<double, 3, 3> mRp;

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

}; // Class MeshRotationUtility

} // namespace Kratos.

#endif // KRATOS_MESH_ROTATION_UTILITY  defined

