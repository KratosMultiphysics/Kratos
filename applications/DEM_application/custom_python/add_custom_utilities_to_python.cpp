/*
==============================================================================
KratosTestApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


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
//   Last modified by:    $Author: G.Casas$
//   Date:                $Date: 2011-6-27 16:07:33 $
//   Revision:            $Revision: 1.2 $
//
//

// System includes 

// External includes 
#include <boost/python.hpp>

// Project includes

#include "includes/model_part.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_strategies/explicit_solver.h"
//#include "custom_strategies/explicit_solver_with_rotation.h"
#include "custom_utilities/create_and_destroy.h"

//namespace Kratos{
//
//namespace Python{
//
//void  AddCustomUtilitiesToPython(){
//    using namespace boost::python;
//    class_<Neighbours_Calculator<2, ParticleType, ParticlePointerType, ParticleVectorType, ParticleWeakVectorType, ParticlePointerVectorType,
//        ParticleWeakIteratorType, ParticleIteratorType, ParticlePointerIteratorType, DistanceVectorType, DistanceIteratorType> >("Neighbours_Calculator_2D", init<>())
//    .def("Search_Neighbours",&Neighbours_Calculator<2, ParticleType, ParticlePointerType, ParticleVectorType, ParticleWeakVectorType, ParticlePointerVectorType,
//        ParticleWeakIteratorType, ParticleIteratorType, ParticlePointerIteratorType, DistanceVectorType, DistanceIteratorType>::Search_Neighbours);
//    class_<Explicit_Solver<2 > >("Explicit_Solver", init<int, double, double, ModelPart& >())
//    .def("Search_Neighbours",&Explicit_Solver<2 >::Search_Neighbours)
//    .def("Calculate_Forces",&Explicit_Solver<2 >::Calculate_Forces)
//    .def("Evolve_Motion",&Explicit_Solver<2 >::Evolve_Motion)
//    .def("Estimate_Time_Step",&Explicit_Solver<2 >::Estimate_Time_Step);
//    }
//}  // namespace Python.
//
//} // Namespace Kratos

namespace Kratos{

namespace Python{

typedef ModelPart::NodesContainerType::iterator PointIterator;
typedef std::vector<array_1d<double, 3 > > ComponentVectorType;
typedef std::vector<array_1d<double, 3 > >::iterator ComponentIteratorType;
/*
typedef CircularParticle CircleType;
typedef CircularParticle::Pointer CirclePointerType;
typedef CircularParticle::ParticleWeakVectorType CircleWeakVectorType;
typedef CircularParticle::ParticleWeakIteratorType CircleWeakIteratorType;
typedef CircularParticle::DistanceVectorType CircleDistanceVectorType;
typedef CircularParticle::DistanceIteratorType CircleDistanceIteratorType;
typedef std::vector<CircleType> CircleVectorType;
typedef std::vector<CircleType>::iterator CircleIteratorType;
typedef std::vector<CircleType::Pointer> CirclePointerVectorType;
typedef std::vector<CircleType::Pointer>::iterator CirclePointerIteratorType;

typedef SphericParticle SphereType;
typedef SphericParticle::Pointer SpherePointerType;
typedef SphericParticle::ParticleWeakVectorType SphereWeakVectorType;
typedef SphericParticle::ParticleWeakIteratorType SphereWeakIteratorType;
typedef SphericParticle::DistanceVectorType SphereDistanceVectorType;
typedef SphericParticle::DistanceIteratorType SphereDistanceIteratorType;
typedef std::vector<SphereType> SphereVectorType;
typedef std::vector<SphereType>::iterator SphereIteratorType;
typedef std::vector<SphereType::Pointer> SpherePointerVectorType;
typedef std::vector<SphereType::Pointer>::iterator SpherePointerIteratorType;

typedef CircularHertzianParticle HertzianCircleType;
typedef CircularHertzianParticle::Pointer HertzianCirclePointerType;
typedef CircularHertzianParticle::ParticleWeakVectorType HertzianCircleWeakVectorType;
typedef CircularHertzianParticle::ParticleWeakIteratorType HertzianCircleWeakIteratorType;
typedef CircularHertzianParticle::DistanceVectorType HertzianCircleDistanceVectorType;
typedef CircularHertzianParticle::DistanceIteratorType HertzianCircleDistanceIteratorType;
typedef std::vector<HertzianCircleType> HertzianCircleVectorType;
typedef std::vector<HertzianCircleType>::iterator HertzianCircleIteratorType;
typedef std::vector<HertzianCircleType::Pointer> HertzianCirclePointerVectorType;
typedef std::vector<HertzianCircleType::Pointer>::iterator HertzianCirclePointerIteratorType;
*/
typedef SphericHertzianParticle HertzianSphereType;
typedef SphericHertzianParticle::Pointer HertzianSpherePointerType;
typedef SphericHertzianParticle::ParticleWeakVectorType HertzianSphereWeakVectorType;
typedef SphericHertzianParticle::ParticleWeakIteratorType HertzianSphereWeakIteratorType;
typedef SphericHertzianParticle::DistanceVectorType HertzianSphereDistanceVectorType;
typedef SphericHertzianParticle::DistanceIteratorType HertzianSphereDistanceIteratorType;
typedef std::vector<HertzianSphereType> HertzianSphereVectorType;
typedef std::vector<HertzianSphereType>::iterator HertzianSphereIteratorType;
typedef std::vector<HertzianSphereType::Pointer> HertzianSpherePointerVectorType;
typedef std::vector<HertzianSphereType::Pointer>::iterator HertzianSpherePointerIteratorType;
/*
typedef SphericRotatingParticle RotatingSphereType;
typedef SphericRotatingParticle::Pointer RotatingSpherePointerType;
typedef SphericRotatingParticle::ParticleWeakVectorType RotatingSphereWeakVectorType;
typedef SphericRotatingParticle::ParticleWeakIteratorType RotatingSphereWeakIteratorType;
typedef SphericRotatingParticle::DistanceVectorType RotatingSphereDistanceVectorType;
typedef SphericRotatingParticle::DistanceIteratorType RotatingSphereDistanceIteratorType;
typedef std::vector<RotatingSphereType> RotatingSphereVectorType;
typedef std::vector<RotatingSphereType>::iterator RotatingSphereIteratorType;
typedef std::vector<RotatingSphereType::Pointer> RotatingSpherePointerVectorType;
typedef std::vector<RotatingSphereType::Pointer>::iterator RotatingSpherePointerIteratorType;

typedef SphericRotatingHertzianParticle RotatingHertzianSphereType;
typedef SphericRotatingHertzianParticle::Pointer RotatingHertzianSpherePointerType;
typedef SphericRotatingHertzianParticle::ParticleWeakVectorType RotatingHertzianSphereWeakVectorType;
typedef SphericRotatingHertzianParticle::ParticleWeakIteratorType RotatingHertzianSphereWeakIteratorType;
typedef SphericRotatingHertzianParticle::DistanceVectorType RotatingHertzianSphereDistanceVectorType;
typedef SphericRotatingHertzianParticle::DistanceIteratorType RotatingHertzianSphereDistanceIteratorType;
typedef std::vector<RotatingHertzianSphereType> RotatingHertzianSphereVectorType;
typedef std::vector<RotatingHertzianSphereType>::iterator RotatingHertzianSphereIteratorType;
typedef std::vector<RotatingHertzianSphereType::Pointer> RotatingHertzianSpherePointerVectorType;
typedef std::vector<RotatingHertzianSphereType::Pointer>::iterator RotatingHertzianSpherePointerIteratorType;

*/
void  AddCustomUtilitiesToPython(){
    using namespace boost::python;

    /*
    class_<Neighbours_Calculator<2, CircleType, CirclePointerType, CircleVectorType, CircleWeakVectorType, CirclePointerVectorType,
        CircleWeakIteratorType, CircleIteratorType, CirclePointerIteratorType, CircleDistanceVectorType, CircleDistanceIteratorType> >("Circles_Neighbours_Calculator", init<>())
    .def("Search_Neighbours",&Neighbours_Calculator<2, CircleType, CirclePointerType, CircleVectorType, CircleWeakVectorType, CirclePointerVectorType,
        CircleWeakIteratorType, CircleIteratorType, CirclePointerIteratorType, CircleDistanceVectorType, CircleDistanceIteratorType>::Search_Neighbours);

     class_<Neighbours_Calculator<2, HertzianCircleType, HertzianCirclePointerType, HertzianCircleVectorType, HertzianCircleWeakVectorType, HertzianCirclePointerVectorType,
        HertzianCircleWeakIteratorType, HertzianCircleIteratorType, HertzianCirclePointerIteratorType, HertzianCircleDistanceVectorType, HertzianCircleDistanceIteratorType> >("Circles_Hertzian_Neighbours_Calculator", init<>())
    .def("Search_Neighbours",&Neighbours_Calculator<2, HertzianCircleType, HertzianCirclePointerType, HertzianCircleVectorType, HertzianCircleWeakVectorType, HertzianCirclePointerVectorType,
        HertzianCircleWeakIteratorType, HertzianCircleIteratorType, HertzianCirclePointerIteratorType, HertzianCircleDistanceVectorType, HertzianCircleDistanceIteratorType>::Search_Neighbours);

    class_<Neighbours_Calculator<3, SphereType, SpherePointerType, SphereVectorType, SphereWeakVectorType, SpherePointerVectorType,
        SphereWeakIteratorType, SphereIteratorType, SpherePointerIteratorType, SphereDistanceVectorType, SphereDistanceIteratorType> >("Spheres_Neighbours_Calculator", init<>())
    .def("Search_Neighbours",&Neighbours_Calculator<3, SphereType, SpherePointerType, SphereVectorType, SphereWeakVectorType, SpherePointerVectorType,
        SphereWeakIteratorType, SphereIteratorType, SpherePointerIteratorType, SphereDistanceVectorType, SphereDistanceIteratorType>::Search_Neighbours);
*/
    class_<Neighbours_Calculator<3, HertzianSphereType, HertzianSpherePointerType, HertzianSphereVectorType, HertzianSphereWeakVectorType, HertzianSpherePointerVectorType,
        HertzianSphereWeakIteratorType, HertzianSphereIteratorType, HertzianSpherePointerIteratorType, HertzianSphereDistanceVectorType, HertzianSphereDistanceIteratorType> >("Spheres_Hertzian_Neighbours_Calculator", init<>())
    .def("Search_Neighbours",&Neighbours_Calculator<3, HertzianSphereType, HertzianSpherePointerType, HertzianSphereVectorType, HertzianSphereWeakVectorType, HertzianSpherePointerVectorType,
        HertzianSphereWeakIteratorType, HertzianSphereIteratorType, HertzianSpherePointerIteratorType, HertzianSphereDistanceVectorType, HertzianSphereDistanceIteratorType>::Search_Neighbours);
/*
     class_<Neighbours_Calculator<3, RotatingSphereType, RotatingSpherePointerType, RotatingSphereVectorType, RotatingSphereWeakVectorType, RotatingSpherePointerVectorType,
        RotatingSphereWeakIteratorType, RotatingSphereIteratorType, RotatingSpherePointerIteratorType, RotatingSphereDistanceVectorType, RotatingSphereDistanceIteratorType> >("Rotating_Spheres_Neighbours_Calculator", init<>())
    .def("Search_Neighbours",&Neighbours_Calculator<3, RotatingSphereType, RotatingSpherePointerType, RotatingSphereVectorType, RotatingSphereWeakVectorType, RotatingSpherePointerVectorType,
        RotatingSphereWeakIteratorType, RotatingSphereIteratorType, RotatingSpherePointerIteratorType, RotatingSphereDistanceVectorType, RotatingSphereDistanceIteratorType> ::Search_Neighbours);

    class_<Neighbours_Calculator<3, RotatingHertzianSphereType, RotatingHertzianSpherePointerType, RotatingHertzianSphereVectorType, RotatingHertzianSphereWeakVectorType, RotatingHertzianSpherePointerVectorType,
        RotatingHertzianSphereWeakIteratorType, RotatingHertzianSphereIteratorType, RotatingHertzianSpherePointerIteratorType, RotatingHertzianSphereDistanceVectorType, RotatingHertzianSphereDistanceIteratorType> >("Rotating_Hertzian_Spheres_Neighbours_Calculator", init<>())
    .def("Search_Neighbours",&Neighbours_Calculator<3, RotatingHertzianSphereType, RotatingHertzianSpherePointerType, RotatingHertzianSphereVectorType, RotatingHertzianSphereWeakVectorType, RotatingHertzianSpherePointerVectorType,
        RotatingHertzianSphereWeakIteratorType, RotatingHertzianSphereIteratorType, RotatingHertzianSpherePointerIteratorType, RotatingHertzianSphereDistanceVectorType, RotatingHertzianSphereDistanceIteratorType> ::Search_Neighbours);
*/
 /*
    class_<Explicit_Solver<2, CircleType > >("Circles_Explicit_Solver", init<int, double, double, ModelPart& >())
    .def("Search_Neighbours",&Explicit_Solver<2, CircleType >::Search_Neighbours)
    .def("Calculate_Forces",&Explicit_Solver<2, CircleType >::Calculate_Forces)
    .def("Evolve_Motion",&Explicit_Solver<2, CircleType >::Evolve_Motion)
    .def("Estimate_Time_Step",&Explicit_Solver<2, CircleType >::Estimate_Time_Step_Circles)
    .def("Get_List_Of_Particle_Pointers",&Explicit_Solver<2, CircleType >::GetListOfParticlePointers, return_internal_reference<>());

    class_<Explicit_Solver<3, SphereType > >("Spheres_Explicit_Solver", init<int, double, double, ModelPart& >())
    .def("Search_Neighbours",&Explicit_Solver<3, SphereType >::Search_Neighbours)
    .def("Get_List_Of_Particle_Pointers",&Explicit_Solver<3, SphereType >::GetListOfParticlePointers, return_internal_reference<>())
    .def("Calculate_Forces",&Explicit_Solver<3, SphereType >::Calculate_Forces)
    .def("Evolve_Motion",&Explicit_Solver<3, SphereType >::Evolve_Motion)
    .def("Estimate_Time_Step",&Explicit_Solver<3, SphereType >::Estimate_Time_Step_Spheres)
    .def("Get_List_Of_Particle_Pointers",&Explicit_Solver<3, SphereType >::GetListOfParticlePointers, return_internal_reference<>());

    class_<Explicit_Solver<2, HertzianCircleType > >("Circles_Hertzian_Explicit_Solver", init<int, double, double, ModelPart& >())
    .def("Search_Neighbours",&Explicit_Solver<2, HertzianCircleType >::Search_Neighbours)
    .def("Calculate_Forces",&Explicit_Solver<2, HertzianCircleType >::Calculate_Forces)
    .def("Evolve_Motion",&Explicit_Solver<2, HertzianCircleType >::Evolve_Motion)
    .def("Estimate_Time_Step",&Explicit_Solver<2, HertzianCircleType >::Estimate_Time_Step_Hertzian_Circles)
    .def("Get_List_Of_Particle_Pointers", &Explicit_Solver < 2, HertzianCircleType >::GetListOfParticlePointers, return_internal_reference<>());
*/

    class_<Explicit_Solver<3, HertzianSphereType > >("Spheres_Hertzian_Explicit_Solver", init<int, double, double, ModelPart& >())
    .def("Search_Neighbours",&Explicit_Solver<3, HertzianSphereType >::Search_Neighbours)
    .def("Set_Initial_Contacts",&Explicit_Solver<3, HertzianSphereType >::Set_Initial_Contacts)
    .def("Calculate_Forces",&Explicit_Solver<3, HertzianSphereType >::Calculate_Forces)
    .def("Evolve_Motion",&Explicit_Solver<3, HertzianSphereType >::Evolve_Motion)
    .def("Estimate_Time_Step",&Explicit_Solver<3, HertzianSphereType >::Estimate_Time_Step_Hertzian_Spheres)
    .def("Get_List_Of_Particle_Pointers", &Explicit_Solver < 3, HertzianSphereType >::GetListOfParticlePointers, return_internal_reference<>());
/*
    class_<Explicit_Solver_With_Rotation<3, RotatingSphereType > >("Rotating_Spheres_Explicit_Solver", init<int, double, double, ModelPart& >())
    .def("Search_Neighbours",&Explicit_Solver_With_Rotation<3, RotatingSphereType >::Search_Neighbours)
    .def("Update_Contacting_Neighbours",&Explicit_Solver_With_Rotation<3, RotatingSphereType >::Update_Contacting_Neighbours)
//    .def("Evolve_Displacements",&Explicit_Solver_With_Rotation<3, RotatingSphereType >::Evolve_Displacements)
    .def("Calculate_Forces",&Explicit_Solver_With_Rotation<3, RotatingSphereType >::Calculate_Forces)
    .def("Evolve_Motion",&Explicit_Solver_With_Rotation<3, RotatingSphereType >::Evolve_Motion)
    .def("Estimate_Time_Step",&Explicit_Solver_With_Rotation<3, RotatingSphereType >::Estimate_Time_Step_Spheres)
    .def("Get_List_Of_Particle_Pointers", &Explicit_Solver_With_Rotation< 3, RotatingSphereType >::GetListOfParticlePointers, return_internal_reference<>());

    class_<Explicit_Solver_With_Rotation<3, RotatingHertzianSphereType > >("Rotating_Hertzian_Spheres_Explicit_Solver", init<int, double, double, ModelPart& >())
    .def("Search_Neighbours",&Explicit_Solver_With_Rotation<3, RotatingHertzianSphereType >::Search_Neighbours)
    .def("Update_Contacting_Neighbours",&Explicit_Solver_With_Rotation<3, RotatingHertzianSphereType >::Update_Contacting_Neighbours)
//    .def("Evolve_Displacements",&Explicit_Solver_With_Rotation<3, RotatingHertzianSphereType >::Evolve_Displacements)
    .def("Calculate_Forces",&Explicit_Solver_With_Rotation<3, RotatingHertzianSphereType >::Calculate_Forces)
    .def("Evolve_Motion",&Explicit_Solver_With_Rotation<3, RotatingHertzianSphereType >::Evolve_Motion)
    .def("Estimate_Time_Step",&Explicit_Solver_With_Rotation<3, RotatingHertzianSphereType >::Estimate_Time_Step_Spheres)
    .def("Get_List_Of_Particle_Pointers", &Explicit_Solver_With_Rotation< 3, RotatingHertzianSphereType >::GetListOfParticlePointers, return_internal_reference<>());

    class_<Particle_Creator_Destructor<3, CircleType, CirclePointerType, CircleVectorType, CircleWeakVectorType, CirclePointerVectorType,
        CircleWeakIteratorType, CircleIteratorType, CirclePointerIteratorType, CircleDistanceVectorType, CircleDistanceIteratorType> >("Circles_Creator_Destructor", init<>())
    .def("Calculate_Surrounding_Bounding_Box", &Particle_Creator_Destructor<3, CircleType, CirclePointerType, CircleVectorType, CircleWeakVectorType, CirclePointerVectorType,
        CircleWeakIteratorType, CircleIteratorType, CirclePointerIteratorType, CircleDistanceVectorType, CircleDistanceIteratorType>::CalculateSurroundingBoundingBox)
    .def("Destroy_Distant_Particles", &Particle_Creator_Destructor<3, CircleType, CirclePointerType, CircleVectorType, CircleWeakVectorType, CirclePointerVectorType,
        CircleWeakIteratorType, CircleIteratorType, CirclePointerIteratorType, CircleDistanceVectorType, CircleDistanceIteratorType>::DestroyDistantParticles);
    
     class_<Particle_Creator_Destructor<3, HertzianCircleType, HertzianCirclePointerType, HertzianCircleVectorType, HertzianCircleWeakVectorType, HertzianCirclePointerVectorType,
        HertzianCircleWeakIteratorType, HertzianCircleIteratorType, HertzianCirclePointerIteratorType, HertzianCircleDistanceVectorType, HertzianCircleDistanceIteratorType> >("Circles_Hertzian_Creator_Destructor", init<>())
    .def("Calculate_Surrounding_Bounding_Box", &Particle_Creator_Destructor<3, HertzianCircleType, HertzianCirclePointerType, HertzianCircleVectorType, HertzianCircleWeakVectorType, HertzianCirclePointerVectorType,
        HertzianCircleWeakIteratorType, HertzianCircleIteratorType, HertzianCirclePointerIteratorType, HertzianCircleDistanceVectorType, HertzianCircleDistanceIteratorType>::CalculateSurroundingBoundingBox)
    .def("Destroy_Distant_Particles", &Particle_Creator_Destructor<3, HertzianCircleType, HertzianCirclePointerType, HertzianCircleVectorType, HertzianCircleWeakVectorType, HertzianCirclePointerVectorType,
        HertzianCircleWeakIteratorType, HertzianCircleIteratorType, HertzianCirclePointerIteratorType, HertzianCircleDistanceVectorType, HertzianCircleDistanceIteratorType>::DestroyDistantParticles);

    class_<Particle_Creator_Destructor<3, SphereType, SpherePointerType, SphereVectorType, SphereWeakVectorType, SpherePointerVectorType,
        SphereWeakIteratorType, SphereIteratorType, SpherePointerIteratorType, SphereDistanceVectorType, SphereDistanceIteratorType> >("Spheres_Creator_Destructor", init<>())
    .def("Calculate_Surrounding_Bounding_Box", &Particle_Creator_Destructor<3, SphereType, SpherePointerType, SphereVectorType, SphereWeakVectorType, SpherePointerVectorType,
        SphereWeakIteratorType, SphereIteratorType, SpherePointerIteratorType, SphereDistanceVectorType, SphereDistanceIteratorType>::CalculateSurroundingBoundingBox)
    .def("Destroy_Distant_Particles", &Particle_Creator_Destructor<3, SphereType, SpherePointerType, SphereVectorType, SphereWeakVectorType, SpherePointerVectorType,
        SphereWeakIteratorType, SphereIteratorType, SpherePointerIteratorType, SphereDistanceVectorType, SphereDistanceIteratorType>::DestroyDistantParticles);
*/
    class_<Particle_Creator_Destructor<3, HertzianSphereType, HertzianSpherePointerType, HertzianSphereVectorType, HertzianSphereWeakVectorType, HertzianSpherePointerVectorType,
        HertzianSphereWeakIteratorType, HertzianSphereIteratorType, HertzianSpherePointerIteratorType, HertzianSphereDistanceVectorType, HertzianSphereDistanceIteratorType> >("Spheres_Hertzian_Creator_Destructor", init<>())
    .def("Calculate_Surrounding_Bounding_Box", &Particle_Creator_Destructor<3, HertzianSphereType, HertzianSpherePointerType, HertzianSphereVectorType, HertzianSphereWeakVectorType, HertzianSpherePointerVectorType,
        HertzianSphereWeakIteratorType, HertzianSphereIteratorType, HertzianSpherePointerIteratorType, HertzianSphereDistanceVectorType, HertzianSphereDistanceIteratorType>::CalculateSurroundingBoundingBox)
    .def("Destroy_Distant_Particles", &Particle_Creator_Destructor<3, HertzianSphereType, HertzianSpherePointerType, HertzianSphereVectorType, HertzianSphereWeakVectorType, HertzianSpherePointerVectorType,
        HertzianSphereWeakIteratorType, HertzianSphereIteratorType, HertzianSpherePointerIteratorType, HertzianSphereDistanceVectorType, HertzianSphereDistanceIteratorType>::DestroyDistantParticles);
/*
    class_<Particle_Creator_Destructor<3, RotatingSphereType, RotatingSpherePointerType, RotatingSphereVectorType, RotatingSphereWeakVectorType, RotatingSpherePointerVectorType,
        RotatingSphereWeakIteratorType, RotatingSphereIteratorType, RotatingSpherePointerIteratorType, RotatingSphereDistanceVectorType, RotatingSphereDistanceIteratorType> >("Rotating_Spheres_Creator_Destructor", init<>())
    .def("Calculate_Surrounding_Bounding_Box", &Particle_Creator_Destructor<3, RotatingSphereType, RotatingSpherePointerType, RotatingSphereVectorType, RotatingSphereWeakVectorType, RotatingSpherePointerVectorType,
        RotatingSphereWeakIteratorType, RotatingSphereIteratorType, RotatingSpherePointerIteratorType, RotatingSphereDistanceVectorType, RotatingSphereDistanceIteratorType>::CalculateSurroundingBoundingBox)
    .def("Destroy_Distant_Particles", &Particle_Creator_Destructor<3, RotatingSphereType, RotatingSpherePointerType, RotatingSphereVectorType, RotatingSphereWeakVectorType, RotatingSpherePointerVectorType,
        RotatingSphereWeakIteratorType, RotatingSphereIteratorType, RotatingSpherePointerIteratorType, RotatingSphereDistanceVectorType, RotatingSphereDistanceIteratorType>::DestroyDistantParticles);

    class_<Particle_Creator_Destructor<3, RotatingHertzianSphereType, RotatingHertzianSpherePointerType, RotatingHertzianSphereVectorType, RotatingHertzianSphereWeakVectorType, RotatingHertzianSpherePointerVectorType,
        RotatingHertzianSphereWeakIteratorType, RotatingHertzianSphereIteratorType, RotatingHertzianSpherePointerIteratorType, RotatingHertzianSphereDistanceVectorType, RotatingHertzianSphereDistanceIteratorType> >("Rotating_Hertzian_Spheres_Creator_Destructor", init<>())
    .def("Calculate_Surrounding_Bounding_Box", &Particle_Creator_Destructor<3, RotatingHertzianSphereType, RotatingHertzianSpherePointerType, RotatingHertzianSphereVectorType, RotatingHertzianSphereWeakVectorType, RotatingHertzianSpherePointerVectorType,
        RotatingHertzianSphereWeakIteratorType, RotatingHertzianSphereIteratorType, RotatingHertzianSpherePointerIteratorType, RotatingHertzianSphereDistanceVectorType, RotatingHertzianSphereDistanceIteratorType>::CalculateSurroundingBoundingBox)
    .def("Destroy_Distant_Particles", &Particle_Creator_Destructor<3, RotatingHertzianSphereType, RotatingHertzianSpherePointerType, RotatingHertzianSphereVectorType, RotatingHertzianSphereWeakVectorType, RotatingHertzianSpherePointerVectorType,
        RotatingHertzianSphereWeakIteratorType, RotatingHertzianSphereIteratorType, RotatingHertzianSpherePointerIteratorType, RotatingHertzianSphereDistanceVectorType, RotatingHertzianSphereDistanceIteratorType>::DestroyDistantParticles);
*/
    /*
    class_<std::vector<CircleType::Pointer > >("CirclesPointersVector", init<>());
    class_<std::vector<HertzianCircleType::Pointer > >("CirclesHertzianPointersVector", init<>());
    class_<std::vector<SphereType::Pointer > >("SpheresPointersVectors", init<>());
   */
       class_<std::vector<HertzianSphereType::Pointer > >("SpheresHertzianPointersVector", init<>());
    /*
     * class_<std::vector<RotatingSphereType::Pointer > >("RotatingSpheresPointersVector", init<>());
    class_<std::vector<RotatingHertzianSphereType::Pointer > >("RotatingHertzianSpheresPointersVector", init<>());
     *
     */

    }
}  // namespace Python.

} // Namespace Kratos
