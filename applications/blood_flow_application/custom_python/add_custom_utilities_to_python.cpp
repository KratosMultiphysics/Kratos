/*
==============================================================================
KratosBloodFlowApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2013
Pooyan Dadvand, Riccardo Rossi, Eduardo Soudah
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
esoudah@cimne.upc.edu
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
//   Last modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.2 $
//
//


// System includes 

// External includes 
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "spaces/ublas_space.h"
#include "custom_utilities/artery_time_integrator.h"



namespace Kratos
{
	
namespace Python
{
  
  void CreateNewCondition2nodes(char* ConditionName, ModelPart& model_part, int cond_id, int prop_id,  Node<3>::Pointer pnode1, Node<3>::Pointer pnode2)
  {
      Condition::GeometryType::PointsArrayType plist_of_nodes;
      plist_of_nodes.push_back(pnode1);
      plist_of_nodes.push_back(pnode2);      
      //Condition::GeometryType::Pointer pgeom = new Line3D2<Node<3> >( pnode1, pnode2 );
      Properties::Pointer properties = model_part.GetMesh().pGetProperties(prop_id);
      Condition::Pointer pnew_cond = (KratosComponents<Condition>::Get(ConditionName)).Create(cond_id, plist_of_nodes, properties );      
      model_part.AddCondition( pnew_cond );
  }
  
    void CreateNewCondition1node(char* ConditionName, ModelPart& model_part, int cond_id, int prop_id,  Node<3>::Pointer pnode1)
  {
      Condition::GeometryType::PointsArrayType plist_of_nodes;
      plist_of_nodes.push_back(pnode1);
      
      //Condition::GeometryType::Pointer pgeom = new Line3D2<Node<3> >( pnode1, pnode2 );
      Properties::Pointer properties = model_part.GetMesh().pGetProperties(prop_id);
      Condition::Pointer pnew_cond = (KratosComponents<Condition>::Get(ConditionName)).Create(cond_id, plist_of_nodes, properties );
      
      model_part.AddCondition( pnew_cond );
  }
	
  void  AddCustomUtilitiesToPython()
  {
	using namespace boost::python;


    class_<ArteryTimeIntegrator > ("ArteryTimeIntegrator", init<>())
    .def("Initialize", &ArteryTimeIntegrator::Initialize)
    .def("SolveStep", &ArteryTimeIntegrator::SolveStep)
    .def("EstimateDeltaTime", &ArteryTimeIntegrator::EstimateDeltaTime)
    .def("Element_minLength", &ArteryTimeIntegrator::Element_minLength)
    .def("CheckCardiacCovergence", &ArteryTimeIntegrator::CheckCardiacCovergence)
    .def("ComputePressure", &ArteryTimeIntegrator::ComputePressure)    
    ;
    
    def("CreateNewCondition",CreateNewCondition1node);
    def("CreateNewCondition",CreateNewCondition2nodes);
    
    
  }
	




}  // namespace Python.

} // Namespace Kratos

