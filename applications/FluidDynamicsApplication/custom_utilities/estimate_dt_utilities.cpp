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
#include "utilities/element_size_calculator.h"
#include "utilities/geometry_utilities.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "estimate_dt_utilities.h"


namespace Kratos
{

    void EstimateDtUtility::SetCFL(const double CFL)
    {
        mCFL = CFL;
    }

    void EstimateDtUtility::SetDtMin(const double DtMin)
    {
        mDtMin = DtMin;
    }

    void EstimateDtUtility::SetDtMax(const double DtMax)
    {
        mDtMax = DtMax;
    }

    double EstimateDtUtility::EstimateDt() const
    {
        KRATOS_TRY;

        // Get the projected element size function according to the corresponding geometry
        // Note that in here it is assumed that all the elements in the model part feature the same geometry
        const auto& r_geom = mrModelPart.ElementsBegin()->GetGeometry();
        ElementSizeFunctionType projected_h_func = EstimateDtUtility::GetMinimumElementSizeFunction(r_geom);

        // Obtain the maximum CFL
        const double current_dt = mrModelPart.GetProcessInfo().GetValue(DELTA_TIME);
        const double current_cfl = block_for_each<MaxReduction<double>>(mrModelPart.Elements(), [&](Element& rElement) -> double {
            return CalculateElementCFL(rElement, projected_h_func, current_dt);
        });

        // Calculate the new time increment from the maximum local CFL in the mesh
        double new_dt = 0.0;
        if (current_cfl < 1e-10) {
            // Avoid division by 0 when the maximum CFL number is close to 0 (e.g. problem initialization)
            KRATOS_INFO("EstimateDtUtility") << "Setting minimum delta time " << mDtMin << " as current time step." << std::endl;
            new_dt = mDtMin;
        } else {
            // Compute new Dt
            new_dt = mCFL * current_dt / current_cfl;
            // Limit max and min Dt
            if (new_dt > mDtMax) {
                new_dt = mDtMax;
            } else if (new_dt < mDtMin) {
                new_dt = mDtMin;
            }
        }

        // Perform MPI sync if needed
        new_dt = mrModelPart.GetCommunicator().GetDataCommunicator().MinAll(new_dt);

        return new_dt;

        KRATOS_CATCH("")
    }

    void EstimateDtUtility::CalculateLocalCFL(ModelPart& rModelPart)
    {
        KRATOS_TRY;

        // Get the projected element size function according to the corresponding geometry
        // Note that in here it is assumed that all the elements in the model part feature the same geometry
        const auto& r_geom = rModelPart.ElementsBegin()->GetGeometry();
        ElementSizeFunctionType projected_h_func = EstimateDtUtility::GetMinimumElementSizeFunction(r_geom);

        // Calculate the CFL number in each element
        const double current_dt = rModelPart.GetProcessInfo().GetValue(DELTA_TIME);
        block_for_each(rModelPart.Elements(), [&](Element& rElement){
            const double element_cfl = EstimateDtUtility::CalculateElementCFL(rElement, projected_h_func, current_dt);
            rElement.SetValue(CFL_NUMBER, element_cfl);
        });

        KRATOS_CATCH("")
    }

    double EstimateDtUtility::CalculateElementCFL(
        const Element &rElement,
        const ElementSizeFunctionType& rElementSizeCalculator,
        const double Dt)
    {
        // Calculate the midpoint velocity
        const auto& r_geometry = rElement.GetGeometry();
        const unsigned int n_nodes = r_geometry.PointsNumber();
        array_1d<double,3> element_vel = r_geometry[0].FastGetSolutionStepValue(VELOCITY);
        for (unsigned int i = 1; i < n_nodes; ++i) {
            element_vel += r_geometry[i].FastGetSolutionStepValue(VELOCITY);
        }
        element_vel /= static_cast<double>(n_nodes);

        // Calculate element CFL
        const double h_min = rElementSizeCalculator(r_geometry);
        const double elem_cfl = norm_2(element_vel) * Dt / h_min;

        return elem_cfl;
    }

    typename EstimateDtUtility::ElementSizeFunctionType EstimateDtUtility::GetMinimumElementSizeFunction(const Geometry<Node<3>>& rGeometry)
    {
        ElementSizeFunctionType projected_h_func;
        const auto geometry_type = rGeometry.GetGeometryType();
        switch (geometry_type) {
            case GeometryData::Kratos_Triangle2D3:
                projected_h_func = [&](const Geometry<Node<3>>& rGeometry){return ElementSizeCalculator<2,3>::MinimumElementSize(rGeometry);};
                break;
            case GeometryData::Kratos_Quadrilateral2D4:
                projected_h_func = [&](const Geometry<Node<3>>& rGeometry){return ElementSizeCalculator<2,4>::MinimumElementSize(rGeometry);};
                break;
            case GeometryData::Kratos_Tetrahedra3D4:
                projected_h_func = [&](const Geometry<Node<3>>& rGeometry){return ElementSizeCalculator<3,4>::MinimumElementSize(rGeometry);};
                break;
            case GeometryData::Kratos_Quadrilateral3D8:
                projected_h_func = [&](const Geometry<Node<3>>& rGeometry){return ElementSizeCalculator<3,8>::MinimumElementSize(rGeometry);};
                break;
            default:
                KRATOS_ERROR << "Non supported geometry type." << std::endl;
        }

        return projected_h_func;
    }

    void EstimateDtUtility::save(Serializer& rSerializer) const
    {
        rSerializer.save("mCFL", mCFL);
        rSerializer.save("mDtMax", mDtMax);
        rSerializer.save("mDtMin", mDtMin);
        rSerializer.save("mrModelPart", mrModelPart);
    }

    void EstimateDtUtility::load(Serializer& rSerializer)
    {
        rSerializer.load("mCFL", mCFL);
        rSerializer.load("mDtMax", mDtMax);
        rSerializer.load("mDtMin", mDtMin);
        rSerializer.load("mrModelPart", mrModelPart);
    }

} // namespace Kratos.