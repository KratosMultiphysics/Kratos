//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta, Bodhinanda Chandra
//
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes

// Project includes
#include "includes/define.h"
#include "includes/define_python.h"

#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_constitutive_laws_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"

#include "custom_elements/updated_lagrangian.hpp"
#include "custom_elements/updated_lagrangian_UP.hpp"
#include "custom_elements/updated_lagrangian_quadrilateral.hpp"

#include "geometries/triangle_3d_3.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/hexahedra_3d_8.h"

#include "particle_mechanics_application.h"

namespace Kratos
{

namespace Python
{

using namespace pybind11;

Element::Pointer CreateUpdatedLagragian2D3N()
{
    return Kratos::make_shared<UpdatedLagrangian>( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) );
}

Element::Pointer CreateUpdatedLagragianUP2D3N()
{
    return Kratos::make_shared<UpdatedLagrangianUP>( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) );
}

Element::Pointer CreateUpdatedLagragian3D4N()
{
    return Kratos::make_shared<UpdatedLagrangian>( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) );
}
Element::Pointer CreateUpdatedLagragian2D4N()
{
    return Kratos::make_shared<UpdatedLagrangianQuadrilateral>( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) );

}
Element::Pointer CreateUpdatedLagragian3D8N()
{
    return Kratos::make_shared<UpdatedLagrangianQuadrilateral>( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8 ) ) ) );
}


PYBIND11_MODULE(KratosParticleMechanicsApplication, m)
{
    class_<KratosParticleMechanicsApplication,
        KratosParticleMechanicsApplication::Pointer,
        KratosApplication >(m, "KratosParticleMechanicsApplication")
        .def(init<>())
        ;

    AddCustomStrategiesToPython(m);
    AddCustomUtilitiesToPython(m);
    AddCustomConstitutiveLawsToPython(m);
    AddCustomProcessesToPython(m);

    m.def("CreateUpdatedLagragian2D3N", &CreateUpdatedLagragian2D3N);
    m.def("CreateUpdatedLagragianUP2D3N", &CreateUpdatedLagragianUP2D3N);
    m.def("CreateUpdatedLagragian3D4N", &CreateUpdatedLagragian3D4N);
    m.def("CreateUpdatedLagragian2D4N", &CreateUpdatedLagragian2D4N);
    m.def("CreateUpdatedLagragian3D8N", &CreateUpdatedLagragian3D8N);

    // Registering variables in python
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,  GAUSS_COORD )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, COUNTER);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MP_NUMBER);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MP_BOOL);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, WEIGHT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MP_MASS);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MP_DENSITY);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MP_VOLUME);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MP_KINETIC_ENERGY);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MP_STRAIN_ENERGY);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MP_TOTAL_ENERGY);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MP_PRESSURE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MP_JACOBIAN);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MP_DELTA_PLASTIC_STRAIN );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MP_DELTA_PLASTIC_VOLUMETRIC_STRAIN );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MP_DELTA_PLASTIC_DEVIATORIC_STRAIN );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MP_EQUIVALENT_PLASTIC_STRAIN );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MP_ACCUMULATED_PLASTIC_VOLUMETRIC_STRAIN );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MP_ACCUMULATED_PLASTIC_DEVIATORIC_STRAIN );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MP_CONSTITUTIVE_PRESSURE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MP_MATERIAL_ID);

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, AUX_VELOCITY);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, AUX_ACCELERATION);

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, MP_DISPLACEMENT);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, MP_VELOCITY);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, MP_ACCELERATION);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, AUX_MP_VELOCITY);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, AUX_MP_ACCELERATION);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, MP_VOLUME_ACCELERATION);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, NODAL_INTERNAL_FORCE);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, DISPLACEMENT_AUX);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MP_CAUCHY_STRESS_VECTOR);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MP_ALMANSI_STRAIN_VECTOR);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PREVIOUS_MP_CAUCHY_STRESS_VECTOR);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PREVIOUS_MP_ALMANSI_STRAIN_VECTOR);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MP_CONSTITUTIVE_MATRIX);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, NODAL_MOMENTUM);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, NODAL_INERTIA);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,  DILATANCY_COEFFICIENT );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NODAL_MPRESSURE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, AUX_PRESSURE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, AUX_MP_PRESSURE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,  AUX_R )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,  AUX_T )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,  AUX_R_VEL )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,  AUX_T_VEL )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,  AUX_R_ACC )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,  AUX_T_ACC )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,  NODAL_LUMPED_MASS )

    // Nodal load variables
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,  POINT_LOAD )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,  LINE_LOAD )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,  SURFACE_LOAD )

    // Condition load variables
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,POINT_LOADS_VECTOR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,LINE_LOADS_VECTOR )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,SURFACE_LOADS_VECTOR )
}


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
