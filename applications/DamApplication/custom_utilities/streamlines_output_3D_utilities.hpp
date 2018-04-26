//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Lorenzo Gracia
//
//

#if !defined(KRATOS_STREAMLINES_OUTPUT_3D_UTILITIES)
#define KRATOS_STREAMLINES_OUTPUT_3D_UTILITIES

// System includes
#include <cmath>

// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "utilities/openmp_utils.h"
#include "utilities/math_utils.h"
#include "custom_utilities/solid_mechanics_math_utilities.hpp"

// Application includes
#include "dam_application_variables.h"

namespace Kratos
{

class StreamlinesOutput3DUtilities
{

  public:
    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    StreamlinesOutput3DUtilities() {}

    ///------------------------------------------------------------------------------------

    /// Destructor
    ~StreamlinesOutput3DUtilities() {}

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ComputeOutputStep(ModelPart &r_model_part, const int dimension)
    {

        KRATOS_TRY;

        const int NNodes = static_cast<int>(r_model_part.Nodes().size());
        ModelPart::NodesContainerType::iterator node_begin = r_model_part.NodesBegin();

        array_1d<double, 3> S_Si;
        array_1d<double, 3> S_Sii;
        array_1d<double, 3> S_Siii;
        array_1d<double, 3> Tangential_components;
        array_1d<double, 3> Vsi;
        array_1d<double, 3> Vsiii;
        array_1d<double, 3> Vsi_pos;
        array_1d<double, 3> Vsiii_pos;
        array_1d<double, 3> PrincipalStresses;

#pragma omp parallel for
        for (int i = 0; i < NNodes; i++)
        {

            ModelPart::NodesContainerType::iterator itNode = node_begin + i;

            const Matrix &NodalStress = itNode->FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR);

            Tangential_components[0] = NodalStress(1, 0); // Sxy
            Tangential_components[1] = NodalStress(2, 0); // Sxz
            Tangential_components[2] = NodalStress(2, 1); // Syz

            //We compute the principal Stresses
            Vector PrincipalStresses(dimension);
            noalias(PrincipalStresses) = SolidMechanicsMathUtilities<double>::EigenValuesDirectMethod(NodalStress);

            for (unsigned int i = 0; i < 3; i++)
            {
                S_Si[i] = NodalStress(i, i) - PrincipalStresses[0];
                S_Sii[i] = NodalStress(i, i) - PrincipalStresses[1];
                S_Siii[i] = NodalStress(i, i) - PrincipalStresses[2];
            }

            noalias(Vsi) = this->ComputeVaps(Tangential_components, S_Siii, S_Sii, PrincipalStresses);
            noalias(Vsiii) = this->ComputeVaps(Tangential_components, S_Sii, S_Si, PrincipalStresses);

            for (unsigned int a = 0; a < 3; a++)
            {
                Vsi_pos[a] = Vsi(a) * PrincipalStresses[0];
                Vsiii_pos[a] = Vsiii(a) * PrincipalStresses[2];
            }
            array_1d<double, 3> &StreamlineVi = itNode->FastGetSolutionStepValue(Vi_POSITIVE);
            noalias(StreamlineVi) = Vsi_pos;
            array_1d<double, 3> &StreamlineViii = itNode->FastGetSolutionStepValue(Viii_POSITIVE);
            noalias(StreamlineViii) = Vsiii_pos;
        }

        KRATOS_CATCH("");
    }

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  protected:
    /// Member Variables

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    array_1d<double, 3> ComputeVaps(const array_1d<double, 3> &Tangential_components, const array_1d<double, 3> &S_Siii, const array_1d<double, 3> &S_Sii, const Vector &PrincipalStresses)
    {

        array_1d<double, 3> vap1;
        array_1d<double, 3> vap2;
        array_1d<double, 3> vap3;
        array_1d<double, 3> Result;
        BoundedMatrix<double, 3, 3> AutovectorMatrix;

        // Tangential Components
        const double &Sxy = Tangential_components[0];
        const double &Sxz = Tangential_components[1];
        const double &Syz = Tangential_components[2];

        AutovectorMatrix(0, 0) = S_Sii[0] * S_Siii[0] + Sxy * Sxy + Sxz * Sxz;
        AutovectorMatrix(1, 1) = S_Sii[1] * S_Siii[1] + Sxy * Sxy + Syz * Syz;
        AutovectorMatrix(2, 2) = S_Sii[2] * S_Siii[2] + Sxz * Sxz + Syz * Syz;
        AutovectorMatrix(0, 1) = S_Sii[0] * Sxy + Sxy * S_Siii[1] + Sxz * Syz;
        AutovectorMatrix(1, 0) = AutovectorMatrix(0, 1);
        AutovectorMatrix(0, 2) = S_Sii[0] * Sxz + Sxy * Syz + S_Siii[2] * Sxz;
        AutovectorMatrix(2, 0) = AutovectorMatrix(0, 2);
        AutovectorMatrix(2, 1) = Sxy * Sxz + S_Sii[1] * Syz + Syz * S_Siii[2];
        AutovectorMatrix(1, 2) = AutovectorMatrix(2, 1);

        for (unsigned int j = 0; j < 3; j++)
        {
            vap1[j] = AutovectorMatrix(0, j);
            vap2[j] = AutovectorMatrix(1, j);
            vap3[j] = AutovectorMatrix(2, j);
        }

        // Normalization of the vector
        const double norm_vap1 = sqrt(vap1[0] * vap1[0] + vap1[1] * vap1[1] + vap1[2] * vap1[2]);
        const double norm_vap2 = sqrt(vap2[0] * vap2[0] + vap2[1] * vap2[1] + vap2(2) * vap2[2]);
        const double norm_vap3 = sqrt(vap3[0] * vap3[0] + vap3[1] * vap3[1] + vap3[2] * vap3[2]);

        if (norm_vap1 > 1.0e-12)
        {
            vap1[0] *= 1.0 / norm_vap1;
            vap1[1] *= 1.0 / norm_vap1;
            vap1[2] *= 1.0 / norm_vap1;
        }
        else
        {
            noalias(vap1) = ZeroVector(3);
        }

        if (norm_vap2 > 1.0e-12)
        {
            vap2[0] *= 1.0 / norm_vap2;
            vap2[1] *= 1.0 / norm_vap2;
            vap2[2] *= 1.0 / norm_vap2;
        }
        else
        {
            noalias(vap2) = ZeroVector(3);
        }

        if (norm_vap3 > 1.0e-12)
        {
            vap3[0] *= 1.0 / norm_vap3;
            vap3[1] *= 1.0 / norm_vap3;
            vap3[2] *= 1.0 / norm_vap3;
        }
        else
        {
            noalias(vap3) = ZeroVector(3);
        }

        // Checking if z-component is positive
        for (unsigned int l = 0; l < 3; l++)
        {
            if (vap1[2] > 0.0)
            {
                Result[l] = vap1[l];
            }
            else if (vap2[2] > 0.0)
            {
                Result[l] = vap2[l];
            }
            else if (vap3[2] > 0.0)
            {
                Result[l] = vap3[l];
            }
            else
            {
                Result[l] = vap1[l];
            }
        }

        return Result;
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

}; //Class

} /* namespace Kratos.*/

#endif /* KRATOS_STREAMLINES_OUTPUT_3D_UTILITIES defined */
