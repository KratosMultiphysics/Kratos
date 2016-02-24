//
//   Project Name:        KratosPfemBaseApplication $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:      February 2016 $
//   Revision:            $Revision:            0.0 $
//
//

// System includes

// External includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_bounding_to_python.h"

// Meshers

// Bounding Boxes
#include "custom_bounding/spatial_bounding_box.hpp"

namespace Kratos
{

namespace Python
{

  typedef SpatialBoundingBox                        BoundingBoxBaseType;
  typedef SpatialBoundingBox::Pointer                BoundingBoxPointer;
  typedef std::vector<SpatialBoundingBox::Pointer> BoundingBoxContainer;

  void Push_Back_Bounding_Box( BoundingBoxContainer& ThisBoundingBoxContainer,
			       BoundingBoxPointer ThisBoundingBox )
  {
    ThisBoundingBoxContainer.push_back( ThisBoundingBox );
  }


  void  AddCustomBoundingToPython()
  {

    using namespace boost::python;
    //class that allows 3D adaptive remeshing (inserting and erasing nodes)

    
    //class that allows 2D adaptive remeshing (inserting and erasing nodes)


    //class that allows 2D contact domain spatial search


    //bounding box container
    class_< BoundingBoxContainer >( "BoundingBoxContainer", init<>() )
    .def( "PushBack", Push_Back_Bounding_Box )
    ;

    //bounding box set
    class_<SpatialBoundingBox, SpatialBoundingBox::Pointer, boost::noncopyable > 
      ( "SpatialBoundingBox", 
	init<Vector, double, Vector>() )
      .def("SetDimension",&SpatialBoundingBox::SetDimension)
      ;
     
  }

}  // namespace Python.

} // Namespace Kratos

