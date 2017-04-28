//
//   Project Name:        Kratos       
//   Last Modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_NURBS_BREP_MODELER_H_INCLUDED )
#define  KRATOS_NURBS_BREP_MODELER_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <ctime>

// External includes 
#include <boost/numeric/ublas/io.hpp>

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"
#include "includes/kratos_parameters.h"

#include "../../kratos/spatial_containers/spatial_containers.h"
//#include "../../kratos/utilities/binbased_fast_point_locator.h"
//#include "../../kratos/includes/kratos_flags.h"

//#include "nurbs_utilities.h"

#include "BrepModelGeometryReader.h"
#include "brep/BrepModel.h"

namespace Kratos
{

  ///@name Kratos Globals
  ///@{ 






  ///@} 
  ///@name Type Definitions
  ///@{ 
/**
 * Typedefs for search
 */
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
  class NurbsBrepModeler
  {
  public:
    ///@name Type Definitions
    ///@{
    // Variables definition 
    typedef Node<3> NodeType;
    typedef std::vector<NodeType::Pointer> NodeVector;
    typedef std::vector<NodeType::Pointer>::iterator NodeIterator;
    typedef std::vector<double>                          DistanceVector;
    typedef std::vector<double>::iterator                DistanceIterator;

    typedef std::vector<BrepFace> BrepFacesVector;
    typedef std::vector<BrepEdge> BrepEdgesVector;

    typedef std::vector<BrepModel> BrepModelVector;

    //Search Tree
    // KDtree definitions
    typedef Bucket< 3, NodeType, NodeVector, NodeType::Pointer, NodeIterator, DistanceIterator > BucketType;
    typedef Tree< KDTreePartition<BucketType> > tree;

    /// Pointer definition of KratosNurbsTestcaseApplication
    KRATOS_CLASS_POINTER_DEFINITION(NurbsBrepModeler);

    ///@}
    ///@name functions 
    ///@{ 
    void CreateIntegrationDomain(const int& shapefunction_order, ModelPart& model_part);

    ///@} 
    ///@name Life Cycle 
    ///@{ 
        
   
    /// Constructor.
    NurbsBrepModeler(BrepModelGeometryReader& brep_model_geometry_reader, ModelPart& model_part);

        //NurbsBrepModeler();
    /// Destructor.
    virtual ~NurbsBrepModeler();
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
    //       static const ApplicationCondition  msApplicationCondition; 
    ///@} 
    ///@name Member Variables
    ///@{ 
    ModelPart m_model_part;
    BrepModelVector m_brep_model_vector;
    //TreeType m_search_tree;

    ///@} 
    ///@name Private Operators
    ///@{ 

    ///@} 
    ///@name Private Operations
    ///@{ 
    void CreateMeshedPoints(ModelPart& model_part);
    BrepFace& GetFace(const unsigned int face_id);
    //Tree< KDTreePartition<BucketType> > CreateSearchTree(ModelPart model_part);
    void MapNode(const Node<3>::Pointer& node, Node<3>::Pointer& node_on_geometry);
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
    NurbsBrepModeler& operator=(NurbsBrepModeler const& rOther);

    /// Copy constructor.
    NurbsBrepModeler(NurbsBrepModeler const& rOther);
    ///@}    

  }; // Class NurbsBrepModeler 

}  // namespace Kratos.

#endif // KRATOS_NURBS_BREP_MODELER_APPLICATION_H_INCLUDED  defined 


