//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Last modified by:    $Author:                JMCarbonell $
//   Date:                $Date:                    July 2013 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_MODELER_H_INCLUDED )
#define  KRATOS_MODELER_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <stdlib.h>

#include <boost/timer.hpp>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/model_part.h"
#include "spatial_containers/spatial_containers.h"
#include "containers/variables_list_data_value_container.h"

#include "processes/find_nodal_h_process.h"

#include "custom_processes/elemental_neighbours_search_process.hpp"
#include "custom_processes/nodal_neighbours_search_process.hpp"

#include "custom_utilities/boundary_normals_calculation_utilities.hpp"
#include "custom_utilities/mesh_data_transfer_utilities.hpp"
#include "custom_utilities/modeler_utilities.hpp"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class Modeler
{
public:

    /**
     * Flags related to the meshing parameters
     */

    //meshing options
    KRATOS_DEFINE_LOCAL_FLAG (REMESH);
    KRATOS_DEFINE_LOCAL_FLAG (RECONNECT);
    KRATOS_DEFINE_LOCAL_FLAG (REFINE_MESH);

    KRATOS_DEFINE_LOCAL_FLAG (CONSTRAINED_MESH);
    KRATOS_DEFINE_LOCAL_FLAG (BOUNDARIES_SEARCH);
    KRATOS_DEFINE_LOCAL_FLAG (NEIGHBOURS_SEARCH);

    KRATOS_DEFINE_LOCAL_FLAG (SET_DOF);

    KRATOS_DEFINE_LOCAL_FLAG (CONTACT_SEARCH);

    //refining options
    KRATOS_DEFINE_LOCAL_FLAG (SELECT_ELEMENTS);
    KRATOS_DEFINE_LOCAL_FLAG (PASS_ALPHA_SHAPE);

    KRATOS_DEFINE_LOCAL_FLAG (REFINE_INSERT_NODES);
    KRATOS_DEFINE_LOCAL_FLAG (REFINE_ADD_NODES);

    KRATOS_DEFINE_LOCAL_FLAG (REFINE_ELEMENTS);
    KRATOS_DEFINE_LOCAL_FLAG (REFINE_BOUNDARY);

    KRATOS_DEFINE_LOCAL_FLAG (REMOVE_NODES);
    KRATOS_DEFINE_LOCAL_FLAG (REMOVE_ON_BOUNDARY);

    KRATOS_DEFINE_LOCAL_FLAG (CRITERION_ERROR);
    KRATOS_DEFINE_LOCAL_FLAG (CRITERION_ENERGY);
    KRATOS_DEFINE_LOCAL_FLAG (CRITERION_DISTANCE);

    KRATOS_DEFINE_LOCAL_FLAG (ENGAGED_NODES);
    KRATOS_DEFINE_LOCAL_FLAG (WALL_TIP);

    ///@name Type Definitions
    ///@{

    /// Pointer definition of Modeler
    KRATOS_CLASS_POINTER_DEFINITION(Modeler);

    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    //typedef ModelPart::GeometricalDataContainerType GeometricalDataContainerType;

    //typedef GeometricalDataContainerType::GeometricalDataType GeometricalDataType;

    //typedef ModelPart::GeometryType GeometryType;

    //typedef ModelPart::GeometriesContainerType GeometriesContainerType;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Modeler() {}

    /// Destructor.
    virtual ~Modeler() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    /*	  void CollapsePoints(ModelPart& rThisModelPart, double Tolerance)
    {
      double distance;
      Point3D<Node<3> >::Pointer p_founded_point;

      GeometriesContainerType& geometries = rThisModelPart.Geometries();
      GeometriesContainerType::PointsArrayType& r_points_array = geometries.Points();

      if(geometries.NumberOfGeometries() == 0)
    	  return;

      for(GeometriesContainerType::GeometryIterator i_geometry = geometries.GeometriesBegin() ;
    	  i_geometry != geometries.GeometriesEnd() ; i_geometry++)
      {
      // At this moment a brute force search is used
    	  for(GeometryType::iterator i_point = i_geometry->begin() ; i_point != i_geometry->end() ; i_point++)
    	  {
    		  bool founded = false;
    		  for(GeometriesContainerType::PointIterator j_point = r_points_array.begin() ;
    			  j_point != r_points_array.end() ; j_point++)
    		  {
    			distance = j_point->Distance(*i_point);
    			if(distance < Tolerance)
    			{
    				founded = true;
    				p_founded_point = *(j_point.base());
    				break;
    			}
    		  }
    		  if(founded)
    		  {
    			  *(i_point.base()) = p_founded_point->pGetPoint(0);
    		  }
    		  else
    		  {
    			  r_points_array.push_back(Point3D<Node<3> >(*(i_point.base())));
    		  }
    	  }
      }

    }
    */


    virtual void GenerateMesh(ModelPart& ThisModelPart, Element const& rReferenceElement, Condition const& rReferenceBoundaryCondition)
    {
        KRATOS_ERROR(std::logic_error, "This modeler CAN NOT be used for mesh generation.", "");
    }

    virtual void GenerateNodes(ModelPart& ThisModelPart)
    {
        KRATOS_ERROR(std::logic_error, "This modeler CAN NOT be used for node generation.", "");
    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "Modeler";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    Modeler& operator=(Modeler const& rOther);

    /// Copy constructor.
    Modeler(Modeler const& rOther);


    ///@}

}; // Class Modeler

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  Modeler& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const Modeler& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_MODELER_H_INCLUDED  defined 


