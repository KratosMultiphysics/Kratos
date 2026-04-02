//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Max Friedrichs Dachale, Juan Ignacio Camarotti, Ricky Aristio
//

// System includes

// External includes

// Project includes

// Application includes
#include "custom_elements/beam_bernoulli_element_3D.h"

namespace Kratos
{
    void BeamBernoulliElement3D::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        KRATOS_TRY;
        const auto& r_geometry = this->GetGeometry();
        const SizeType number_of_control_points = r_geometry.size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(4 * number_of_control_points);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_X));
        }
        KRATOS_CATCH("")
    }

    void BeamBernoulliElement3D::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        KRATOS_TRY;

        const auto& r_geometry = this->GetGeometry();
        const SizeType number_of_control_points = r_geometry.size();

        if (rResult.size() != 4 * number_of_control_points)
            rResult.resize(4 * number_of_control_points, false);

        const IndexType pos = r_geometry[0].GetDofPosition(DISPLACEMENT_X);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const IndexType index = i * 4;
            rResult[index] = r_geometry[i].GetDof(DISPLACEMENT_X, pos).EquationId();
            rResult[index + 1] = r_geometry[i].GetDof(DISPLACEMENT_Y, pos + 1).EquationId();
            rResult[index + 2] = r_geometry[i].GetDof(DISPLACEMENT_Z, pos + 2).EquationId();
            rResult[index + 3] = r_geometry[i].GetDof(ROTATION_X, pos + 3).EquationId();
        }

        KRATOS_CATCH("")
    };

    void BeamBernoulliElement3D::ComputeGAxial(
        const IndexType IntegrationPointIndex,
        Matrix& rGAxial,
        KinematicVariables& rKinematicVariables) const
    {
        KRATOS_TRY;
        KRATOS_CATCH("");
    }

} // Namespace Kratos