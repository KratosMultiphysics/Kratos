//
// Author:
//
// Project includes
#include "glued_to_wall_scheme.h"

namespace Kratos {

    GluedToWallScheme::GluedToWallScheme(Condition* p_wall, SphericParticle* p_sphere) {

        mCondition = p_wall;

        array_1d<double, 3> wall_to_sphere_node;
        noalias(wall_to_sphere_node) = p_sphere->GetGeometry()[0].Coordinates() - p_wall->GetGeometry()[0].Coordinates();

        array_1d<double, 3> unitary_normal_to_wall;
        mCondition->GetGeometry().UnitNormal(unitary_normal_to_wall);

        mDistanceSignedWithNormal = DEM_INNER_PRODUCT_3(wall_to_sphere_node, unitary_normal_to_wall);

        DEM_MULTIPLY_BY_SCALAR_3(unitary_normal_to_wall, mDistanceSignedWithNormal);
        array_1d<double, 3> inner_point;
        noalias(inner_point) = wall_to_sphere_node - unitary_normal_to_wall;

        array_1d<double, 3> point_local_coordinates;
        mCondition->GetGeometry().PointLocalCoordinates(point_local_coordinates, inner_point);
        mShapeFunctionsValues.resize(3);
        mCondition->GetGeometry().ShapeFunctionsValues(mShapeFunctionsValues, point_local_coordinates);

    }

    void GluedToWallScheme::SetTranslationalIntegrationSchemeInProperties(Properties::Pointer pProp, bool verbose) const {
//         if(verbose) KRATOS_INFO("DEM") << "Assigning SymplecticEulerScheme to properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_TRANSLATIONAL_INTEGRATION_SCHEME_POINTER, this->CloneShared());
    }

    void GluedToWallScheme::SetRotationalIntegrationSchemeInProperties(Properties::Pointer pProp, bool verbose) const {
//         if(verbose) KRATOS_INFO("DEM") << "Assigning SymplecticEulerScheme to properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_ROTATIONAL_INTEGRATION_SCHEME_POINTER, this->CloneShared());
    }

    void GluedToWallScheme::Move(Node<3> & i, const double delta_t, const double force_reduction_factor, const int StepFlag) {
        if (i.Is(DEMFlags::BELONGS_TO_A_CLUSTER)) return;

        array_1d<double, 3> inner_point;
        noalias(inner_point) = mShapeFunctionsValues[0] * mCondition->GetGeometry()[0] + mShapeFunctionsValues[1] * mCondition->GetGeometry()[1] + mShapeFunctionsValues[2] * mCondition->GetGeometry()[2];

        array_1d<double, 3> unitary_normal_to_wall;
        mCondition->GetGeometry().UnitNormal(unitary_normal_to_wall);
        DEM_MULTIPLY_BY_SCALAR_3(unitary_normal_to_wall, mDistanceSignedWithNormal);

        array_1d<double, 3> & coordinates = i.Coordinates();
        array_1d<double, 3> old_coordinates;
        noalias(old_coordinates) = coordinates;
        noalias(coordinates) = inner_point + unitary_normal_to_wall;
        array_1d<double, 3> & delta_displacement = i.FastGetSolutionStepValue(DELTA_DISPLACEMENT);
        noalias(delta_displacement) = coordinates - old_coordinates;
        noalias(i.FastGetSolutionStepValue(DISPLACEMENT)) += delta_displacement;
    }

    void GluedToWallScheme::Rotate(Node<3> & i, const double delta_t, const double moment_reduction_factor, const int StepFlag) {
        if (i.Is(DEMFlags::BELONGS_TO_A_CLUSTER)) return;

    }

} //namespace Kratos
