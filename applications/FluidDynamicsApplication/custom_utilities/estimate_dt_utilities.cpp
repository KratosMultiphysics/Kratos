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
//                   Ruben Zorrilla
//
//


// System includes


// External includes


// Project includes
#include "includes/cfd_variables.h"
#include "utilities/geometry_utilities.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "estimate_dt_utilities.h"


namespace Kratos
{

    template<unsigned int TDim>
    void EstimateDtUtility<TDim>::SetCFL(const double CFL)
    {
        mCFL = CFL;
    }

    template<unsigned int TDim>
    void EstimateDtUtility<TDim>::SetDtMin(const double DtMin)
    {
        mDtMin = DtMin;
    }

    template<unsigned int TDim>
    void EstimateDtUtility<TDim>::SetDtMax(const double DtMax)
    {
        mDtMax = DtMax;
    }

    template<unsigned int TDim>
    double EstimateDtUtility<TDim>::EstimateDt() const
    {
        KRATOS_TRY;

        // Obtain the maximum CFL
        GeometryDataContainer geometry_info;
        const double current_dt = mrModelPart.GetProcessInfo().GetValue(DELTA_TIME);
        const double current_cfl = block_for_each<MaxReduction<double>>(mrModelPart.Elements(), geometry_info, [&](Element& rElement, GeometryDataContainer& rGeometryInfo) -> double {
            return CalculateElementCFL(rElement, rGeometryInfo, current_dt);
        });

        double NewDt = 0.0;

        // Avoid division by 0 when the maximum CFL number is close to 0 (e.g. problem initialization)
        if (current_cfl < 1e-10)
        {
            KRATOS_INFO("EstimateDtUtility") << "Setting minimum delta time " << mDtMin << " as current time step." << std::endl;
            NewDt = mDtMin;
        }
        else
        {
            // Compute new Dt
            NewDt = mCFL * current_dt / current_cfl;
            // Limit max and min Dt
            if (NewDt > mDtMax)
            {
                NewDt = mDtMax;
            }
            else if (NewDt < mDtMin)
            {
                NewDt = mDtMin;
            }
        }

        // Perform MPI sync if needed
        NewDt = mrModelPart.GetCommunicator().GetDataCommunicator().MinAll(NewDt);

        return NewDt;

        KRATOS_CATCH("")
    }

    template<unsigned int TDim>
    void EstimateDtUtility<TDim>::CalculateLocalCFL(ModelPart& rModelPart)
    {
        KRATOS_TRY;

        GeometryDataContainer geometry_info;
        const double current_dt = rModelPart.GetProcessInfo().GetValue(DELTA_TIME);
        block_for_each(rModelPart.Elements(), geometry_info, [&](Element& rElement, GeometryDataContainer& rGeometryInfo){
            const double element_cfl = EstimateDtUtility<TDim>::CalculateElementCFL(rElement, rGeometryInfo, current_dt);
            rElement.SetValue(CFL_NUMBER, element_cfl);
        });

        KRATOS_CATCH("")
    }

    template<unsigned int TDim>
    void EstimateDtUtility<TDim>::CalculateLocalCFL()
    {
        EstimateDtUtility<TDim>::CalculateLocalCFL(mrModelPart);
    }

    template<unsigned int TDim>
    double EstimateDtUtility<TDim>::CalculateElementCFL(
        Element &rElement,
        GeometryDataContainer& rGeometryInfo,
        double Dt)
    {
        double Proj = 0.0;

        // Get the element's geometric parameters
        const auto& r_geometry = rElement.GetGeometry();
        GeometryUtils::CalculateGeometryData(r_geometry, rGeometryInfo.DN_DX, rGeometryInfo.N, rGeometryInfo.Area);

        // Elemental Velocity
        array_1d<double,3> ElementVel = rGeometryInfo.N[0]*r_geometry[0].FastGetSolutionStepValue(VELOCITY);
        for (unsigned int i = 1; i < TDim+1; ++i)
            ElementVel += rGeometryInfo.N[i]*r_geometry[i].FastGetSolutionStepValue(VELOCITY);

        // Calculate u/h as the maximum projection of the velocity along element heights
        for (unsigned int i = 0; i < TDim+1; ++i)
        {
            for (unsigned int d = 0; d < TDim; ++d)
                Proj += ElementVel[d]*rGeometryInfo.DN_DX(i,d);
            Proj = std::abs(Proj);
        }

        return Proj*Dt;
    }

template class EstimateDtUtility<2>;
template class EstimateDtUtility<3>;

} // namespace Kratos.