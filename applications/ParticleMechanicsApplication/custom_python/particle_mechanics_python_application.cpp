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
#include "custom_elements/updated_lagrangian_axisymmetry.hpp"

#include "geometries/triangle_3d_3.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/hexahedra_3d_8.h"

#include "particle_mechanics_application.h"

namespace Kratos{
namespace Python{

    namespace py = pybind11;

    // Triangular and Tetrahedral 2D and 3D
    Element::Pointer CreateUpdatedLagragian2D3N()
    {
        return Kratos::make_intrusive<UpdatedLagrangian>( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) );
    }

    Element::Pointer CreateUpdatedLagragianUP2D3N()
    {
        return Kratos::make_intrusive<UpdatedLagrangianUP>( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) );
    }

    Element::Pointer CreateUpdatedLagragian3D4N()
    {
        return Kratos::make_intrusive<UpdatedLagrangian>( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) );
    }

    // Quadrilateral and Hexahedral 2D and 3D
    Element::Pointer CreateUpdatedLagragian2D4N()
    {
        return Kratos::make_intrusive<UpdatedLagrangianQuadrilateral>( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) );

    }
    Element::Pointer CreateUpdatedLagragian3D8N()
    {
        return Kratos::make_intrusive<UpdatedLagrangianQuadrilateral>( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8 ) ) ) );
    }

    // Axis Symmetry Element 2D (Triangular and Quadrilateral)
    Element::Pointer CreateUpdatedLagragianAxis2D3N()
    {
        return Kratos::make_intrusive<UpdatedLagrangianAxisymmetry>( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) );
    }

    Element::Pointer CreateUpdatedLagragianAxis2D4N()
    {
        return Kratos::make_intrusive<UpdatedLagrangianAxisymmetry>( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) );
    }

    PYBIND11_MODULE(KratosParticleMechanicsApplication, m)
    {
        py::class_<KratosParticleMechanicsApplication,
            KratosParticleMechanicsApplication::Pointer,
            KratosApplication >(m, "KratosParticleMechanicsApplication")
            .def(py::init<>())
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
        m.def("CreateUpdatedLagragianAxis2D3N", &CreateUpdatedLagragianAxis2D3N);
	    m.def("CreateUpdatedLagragianAxis2D4N", &CreateUpdatedLagragianAxis2D4N);

        // Registering variables in python
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, MP_COORD);
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MP_MASS);
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MP_DENSITY);
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MP_VOLUME);
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MP_POTENTIAL_ENERGY);
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MP_KINETIC_ENERGY);
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MP_STRAIN_ENERGY);
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MP_TOTAL_ENERGY);
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MP_PRESSURE);
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MP_DELTA_PLASTIC_STRAIN );
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MP_DELTA_PLASTIC_VOLUMETRIC_STRAIN );
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MP_DELTA_PLASTIC_DEVIATORIC_STRAIN );
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MP_EQUIVALENT_PLASTIC_STRAIN );
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MP_ACCUMULATED_PLASTIC_VOLUMETRIC_STRAIN );
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MP_ACCUMULATED_PLASTIC_DEVIATORIC_STRAIN );
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MP_MATERIAL_ID);
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PARTICLES_PER_ELEMENT);
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, IGNORE_GEOMETRIC_STIFFNESS);

        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, MPC_COORD);
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MPC_CONDITION_ID);
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MPC_IS_NEUMANN);
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MPC_AREA);
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, MPC_NORMAL);
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, MPC_DISPLACEMENT);
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, MPC_IMPOSED_DISPLACEMENT);
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, MPC_VELOCITY);
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, MPC_IMPOSED_VELOCITY);
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, MPC_ACCELERATION);
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, MPC_IMPOSED_ACCELERATION);
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, MPC_CONTACT_FORCE);
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PARTICLES_PER_CONDITION);

        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, MP_DISPLACEMENT);
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, MP_VELOCITY);
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, MP_ACCELERATION);
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, MP_VOLUME_ACCELERATION);
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, NODAL_INTERNAL_FORCE);
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MP_CAUCHY_STRESS_VECTOR);
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MP_ALMANSI_STRAIN_VECTOR);
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, NODAL_MOMENTUM);
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, NODAL_INERTIA);
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NODAL_MPRESSURE);
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PRESSURE_REACTION);

        // Essential Boundary variables
        KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PENALTY_FACTOR);

        // Nodal load variables
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, POINT_LOAD )
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, LINE_LOAD )
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, SURFACE_LOAD )
    }

}  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
