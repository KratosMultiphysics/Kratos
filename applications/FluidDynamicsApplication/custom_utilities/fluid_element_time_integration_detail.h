//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#if !defined(KRATOS_FLUID_ELEMENT_TIME_INTEGRATION_DETAIL_H)
#define KRATOS_FLUID_ELEMENT_TIME_INTEGRATION_DETAIL_H

// External includes

// Project includes
#include "includes/define.h"
#include "custom_elements/fluid_element.h"

namespace Kratos
{

namespace Internals {

template <class TElementData, bool TDataKnowsAboutTimeIntegration>
class FluidElementTimeIntegrationDetail {
public:
    static void AddTimeIntegratedSystem(
        FluidElement<TElementData>* pElement,
        TElementData& rData,
        Kratos::Matrix& rLHS,
        Kratos::Vector& rRHS);
};

// For Standard data: Time integration is not available ///////////////////////////////////////////

template <class TElementData>
class FluidElementTimeIntegrationDetail<TElementData, false> {
public:
    static void AddTimeIntegratedSystem(
        FluidElement<TElementData>* pElement,
        TElementData& rData,
        Kratos::Matrix& rLHS,
        Kratos::Vector& rRHS)
    {
        KRATOS_TRY;
        KRATOS_ERROR << "Trying to use time-integrated element functions with a "
                        "data type that does not know previous time step data"
                     << std::endl;
        KRATOS_CATCH("");
    }
};

// Specialized time integration (BDF2) ////////////////////////////////////////////////////////////

template <class TElementData>
class FluidElementTimeIntegrationDetail<TElementData, true> {
public:
    static void AddTimeIntegratedSystem(
        FluidElement<TElementData>* pElement,
        TElementData& rData,
        Kratos::Matrix& rLHS,
        Kratos::Vector& rRHS)
    {
        Kratos::Matrix mass_matrix = ZeroMatrix(rLHS.size1(),rLHS.size2());
        Kratos::Matrix velocity_lhs = ZeroMatrix(rLHS.size1(),rLHS.size2());

        pElement->AddVelocitySystem(rData,velocity_lhs,rRHS);
        pElement->AddMassLHS(rData,mass_matrix);

        noalias(rLHS) += rData.bdf0*mass_matrix + velocity_lhs;

        Vector acceleration = ZeroVector(rRHS.size());

        int row = 0;
        const auto& r_velocities = rData.Velocity;
        const auto& r_velocities_step1 = rData.Velocity_OldStep1;
        const auto& r_velocities_step2 = rData.Velocity_OldStep2;

        for (unsigned int i = 0; i < TElementData::NumNodes; ++i) {
            for (unsigned int d = 0; d < TElementData::Dim; ++d)  {
                // Velocity dofs
                acceleration[row] = rData.bdf0*r_velocities(i,d);
                acceleration[row] += rData.bdf1*r_velocities_step1(i,d);
                acceleration[row] += rData.bdf2*r_velocities_step2(i,d);
                ++row;
            }
            ++row; // Pressure dof
        }

        noalias(rRHS) -= prod(mass_matrix,acceleration);
    }
};

} // namespace Internals

} // namespace Kratos

#endif // KRATOS_FLUID_ELEMENT_TIME_INTEGRATION_DETAIL_H