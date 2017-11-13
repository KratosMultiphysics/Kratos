/*
==============================================================================
KratosParticleMechanicsApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

//
//   Project Name:        Kratos
//   Last modified by:    $Author:  ilaria $
//   Date:                $Date: July 2015$
//   Revision:            $Revision: 1.3 $
//
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"

#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_constitutive_laws_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"

#include "custom_elements/updated_lagrangian.hpp"
#include "custom_elements/updated_lagrangian_UP.hpp"
#include "custom_elements/updated_lagrangian_quadrilateral.hpp"
//#include "custom_elements/updated_lagrangian_UP_quadrilateral.hpp"

//#include "custom_elements/total_lagrangian.hpp"
#include "geometries/triangle_3d_3.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/tetrahedra_3d_4.h"


#include "particle_mechanics_application.h"

namespace Kratos
{

namespace Python
{

using namespace boost::python;

Element::Pointer CreateUpdatedLagragian2D3N()
{
    UpdatedLagrangian::Pointer NewElement(
        new UpdatedLagrangian( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ));
    return NewElement;

}

Element::Pointer CreateUpdatedLagragianUP2D3N()
{
    UpdatedLagrangianUP::Pointer NewElement(
        new UpdatedLagrangianUP( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ));
    return NewElement;

}

Element::Pointer CreateUpdatedLagragian3D4N()
{
    UpdatedLagrangian::Pointer NewElement(
        new UpdatedLagrangian( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ));
    return NewElement;

}
Element::Pointer CreateUpdatedLagragian2D4N()
{
    UpdatedLagrangianQuadrilateral::Pointer NewElement(
        new UpdatedLagrangianQuadrilateral( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ));
    return NewElement;

}


//Element::Pointer CreateUpdatedLagragianUP2D4N()
//{
//UpdatedLagrangianUPQuadrilateral::Pointer NewElement(
//new UpdatedLagrangianUPQuadrilateral( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ));
//return NewElement;

//}



BOOST_PYTHON_MODULE(KratosParticleMechanicsApplication)
{

    class_<KratosParticleMechanicsApplication,
           KratosParticleMechanicsApplication::Pointer,
           bases<KratosApplication>, boost::noncopyable >("KratosParticleMechanicsApplication")
           ;

    AddCustomStrategiesToPython();
    AddCustomUtilitiesToPython();
    AddCustomConstitutiveLawsToPython();
    AddCustomProcessesToPython();

    def("CreateUpdatedLagragian2D3N", &CreateUpdatedLagragian2D3N);
    def("CreateUpdatedLagragianUP2D3N", &CreateUpdatedLagragianUP2D3N);
    def("CreateUpdatedLagragian3D4N", &CreateUpdatedLagragian3D4N);
    def("CreateUpdatedLagragian2D4N", &CreateUpdatedLagragian2D4N);
//	def("CreateUpdatedLagragianUP2D4N", &CreateUpdatedLagragianUP2D4N);

    //def("CreateTotalLagragian2D3N", &CreateTotalLagragian2D3N);
    //def("CreateTotalLagragian3D4N", &CreateTotalLagragian3D4N);


    //registering variables in python

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( GAUSS_COORD )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(COUNTER);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(MP_NUMBER);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(MP_BOOL);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(WEIGHT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(MP_MASS);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(MP_DENSITY);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(MP_VOLUME);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(MP_KINETIC_ENERGY);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(MP_STRAIN_ENERGY);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(MP_TOTAL_ENERGY);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(MP_PRESSURE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(MP_JACOBIAN);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(MP_EQUIVALENT_PLASTIC_STRAIN);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(MP_CONSTITUTIVE_PRESSURE);
    //KRATOS_REGISTER_IN_PYTHON_VARIABLE(NODAL_MASS);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(AUX_VELOCITY);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(AUX_ACCELERATION);

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(MP_DISPLACEMENT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(MP_VELOCITY);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(MP_ACCELERATION);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(AUX_MP_VELOCITY);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(AUX_MP_ACCELERATION);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(MP_VOLUME_ACCELERATION);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(NODAL_INTERNAL_FORCE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DISPLACEMENT_AUX);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(MP_CAUCHY_STRESS_VECTOR);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(MP_ALMANSI_STRAIN_VECTOR);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(PREVIOUS_MP_CAUCHY_STRESS_VECTOR);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(PREVIOUS_MP_ALMANSI_STRAIN_VECTOR);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(MP_CONSTITUTIVE_MATRIX);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(NODAL_MOMENTUM);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(NODAL_INERTIA);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( DILATANCY_COEFFICIENT );
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(NODAL_MPRESSURE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(AUX_PRESSURE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(AUX_MP_PRESSURE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( AUX_R )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( AUX_T )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( AUX_R_VEL )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( AUX_T_VEL )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( AUX_R_ACC )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( AUX_T_ACC )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( NODAL_LUMPED_MASS )

}


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
