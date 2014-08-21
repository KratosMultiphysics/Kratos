//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Last modified by:    $Author:                JMCarbonell $
//   Date:                $Date:                    July 2013 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes

// External includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_modelers_to_python.h"

// Meshers
#include "custom_modelers/triangular_mesh_2D_modeler.hpp"
#include "custom_modelers/contact_domain_2D_modeler.hpp"

// Bounding Boxes
#include "custom_modelers/spatial_bounding_box.hpp"
#include "custom_modelers/rigid_wall_bounding_box.hpp"


namespace Kratos
{

namespace Python
{
  //////////////////////////////////////////////////////////////////////////////////////////
  //											//
  //				ADAPTIVE 3D MESHER					//
  //											//
  //////////////////////////////////////////////////////////////////////////////////////////

  //tetrahedron mesher

  //////////////////////////////////////////////////////////////////////////////////////////
  //											//
  //				ADAPTIVE 2D MESHER					//
  //											//
  //////////////////////////////////////////////////////////////////////////////////////////

  //triangle mesher
  void GenerateTriangleMesh(TriangularMesh2DModeler& Mesher, ModelPart& model_part)
  {
    Mesher.GenerateMesh(model_part);
  }

  void GenerateTriangleContactMesh(ContactDomain2DModeler& Mesher, ModelPart& model_part, char* ElementName, char* ConditionName, bool constrained, double alpha_shape, double h_factor, double offset_factor, double penalty_parameter, double stability_parameter, bool friction_active, double nu_static, double nu_dynamic)
  {
    Mesher.GenerateContactMesh(model_part,
			       KratosComponents<Element>::Get(ElementName),
			       KratosComponents<Condition>::Get(ConditionName),
			       constrained, alpha_shape, h_factor, offset_factor,
			       penalty_parameter, stability_parameter,
			       friction_active, nu_static, nu_dynamic);
  }

  void SetRemeshDataOnMesher(TriangularMesh2DModeler& Mesher, char* ElementName, char* ConditionName,bool remesh, bool constrained, bool mesh_smoothing, bool jacobi_smoothing, bool avoid_tip_elements, double alpha_shape, double offset_factor, int domain)
  {
    //this is the maximum number of parameters to compile successfully
    Mesher.SetRemeshData (KratosComponents<Element>::Get(ElementName),
			  KratosComponents<Condition>::Get(ConditionName),
			  remesh, constrained, mesh_smoothing, jacobi_smoothing, 
			  avoid_tip_elements, alpha_shape, offset_factor, domain );
  }

  void SetRefineDataOnMesher(TriangularMesh2DModeler& Mesher, bool refine, double h_factor, double critical_dissipation, double critical_radius, double reference_error, int domain)
  {
    //this is the maximum number of parameters to compile successfully
    Mesher.SetRefineData ( refine, h_factor, critical_dissipation,
			   critical_radius, reference_error, domain);
  }

  void SetRigidWall(TriangularMesh2DModeler& Mesher,RigidWallBoundingBox::Pointer pRigidWall)
  {
    Mesher.SetRigidWall(pRigidWall);
  }

  void SetRefiningBox(TriangularMesh2DModeler& Mesher,double radius, Vector center, Vector velocity)
  {
    Mesher.SetRefiningBox(radius,center,velocity);
  }

  void SetInitialDataOnMesher(TriangularMesh2DModeler& Mesher,int number_of_domains)
  {
    Mesher.Initialize(number_of_domains);
  }

  void TransferContactBoundaryData(ContactDomain2DModeler& Mesher, ModelPart& model_part, bool initial)
  {
    Mesher.TransferContactBoundaryData(model_part,initial);
  }


  //////////////////////////////////////////////////////////////////////////////////////////
  //											//
  //				CONTACT MESHER			                        //
  //											//
  //////////////////////////////////////////////////////////////////////////////////////////



  void  AddCustomModelersToPython()
  {

    using namespace boost::python;
    //class that allows 3D adaptive remeshing (inserting and erasing nodes)

    
    //class that allows 2D adaptive remeshing (inserting and erasing nodes)
    class_<TriangularMesh2DModeler >("TriangularMesh2DModeler",
				   init< >())
      .def("SetInitialMeshData",SetInitialDataOnMesher)
      .def("SetRemeshData",SetRemeshDataOnMesher)
      .def("SetRefineData",SetRefineDataOnMesher)
      .def("SetRigidWall",SetRigidWall)
      .def("SetRefiningBox",SetRefiningBox)
      .def("GenerateMesh",GenerateTriangleMesh)
      .def("SetEchoLevel",&TriangularMesh2DModeler::SetEchoLevel)
      
      ;

    //class that allows 2D contact domain spatial search
    class_<ContactDomain2DModeler >("ContactDomain2DModeler",
				    init< >())
      .def("GenerateContactMesh",GenerateTriangleContactMesh)
      .def("TransferContactData",TransferContactBoundaryData)
      .def("SetEchoLevel",&ContactDomain2DModeler::SetEchoLevel)
      ;


    //bounding box set
    class_<SpatialBoundingBox, boost::noncopyable > 
      ( "SpatialBoundingBox", 
	init<Vector, double, Vector>() )
      .def("SetDimension",&SpatialBoundingBox::SetDimension)
      ;

    class_<RigidWallBoundingBox, boost::noncopyable > 
      ( "RigidWallBoundingBox", 
	init<int, Vector, Vector, Vector, Vector, Matrix, Vector, Vector, Vector>() )
      .def("SetAxisymmetric",&RigidWallBoundingBox::SetAxisymmetric)
      .def("SetDimension",&RigidWallBoundingBox::SetDimension)
      ;
	
     
  }

}  // namespace Python.

} // Namespace Kratos

