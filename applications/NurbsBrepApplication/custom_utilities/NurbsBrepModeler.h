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

#include "BrepModelGeometryReader.h"
#include "brep/BrepModel.h"

namespace Kratos
{
  class NurbsBrepModeler
  {
  public:
    ///@name Type Definitions
    ///@{
    // Variables definition 
    typedef Node<3> NodeType;
    typedef std::vector<NodeType::Pointer>           NodeVector;
    typedef std::vector<NodeType::Pointer>::iterator NodeIterator;
    typedef std::vector<double>                      DistanceVector;
    typedef std::vector<double>::iterator            DistanceIterator;

    typedef std::vector<BrepFace>  BrepFacesVector;
    typedef std::vector<BrepEdge>  BrepEdgesVector;

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

    void LoadGeometry(BrepModelGeometryReader& rBrepModelGeometryReader);

    //void CreateIntegrationDomain(const int& shapefunction_order, ModelPart& model_part);
    void CreateIntegrationDomain(
        Parameters& rIntegrationDomainParameter,
        ModelPart& rIntegrationDomainModelPart);

    void CreateIntegrationDomainElements(
	    Parameters& rIntegrationDomainParameter, ModelPart& rModelPart);

    /* Applies all geometry refinement operations on the patches which are defined in the
    *  RefinementParameters.
    *  @param[in] rRefinementParameters */
    void ApplyGeometryRefinement(Parameters& rRefinementParameters);

    /* Computes the spatial area of each model part and prints it respectiveley */
    void ComputeArea(ModelPart& rModelPart);

    void MapNode(const Node<3>::Pointer& node, Node<3>::Pointer& node_on_geometry, ModelPart& rSearchModelPart);
    void GetInterfaceConditions(ModelPart& rParticleModelPart, ModelPart& rConditionModelPart, ModelPart& rSearchModelPart);
    void GetInterfaceConditionsAdvanced(ModelPart& rParticleModelPart, ModelPart& rIGAModelPart, ModelPart& rInterfaceConditionsModelPart);
    
    //void SelectNewInterfaces(
    //    const BinsIgaConfigure::ResultContainerType& rNeighborResults, 
    //    const int& rNumberOfResults, 
    //    const double radius, 
    //    std::vector<Node<3>::Pointer>& rInterfaceList, 
    //    ModelPart& rModelPart);

    void CalculateSurfaceNormal(
        Vector& rNeighborResults, 
        const Node<3>::Pointer rNode,
        ModelPart& rModelPart);

    void GetUpdatedLocation(ModelPart& rIGAModelPart);
    void GetUpdatedLocationNewModelPart(ModelPart& rIGAModelPart, ModelPart& rIGAModelPartIntegrationDomain);
    ///@} 
    ///@name Life Cycle 
    ///@{ 
    /// Constructor.
    NurbsBrepModeler(Kratos::shared_ptr<ModelPart> rModelPart);

    /// Destructor.
    virtual ~NurbsBrepModeler();
    ///@} 
  protected:
  private:
    ///@name Member Variables
    ///@{ 
	// should be ModelPart::Pointer ??
    Kratos::shared_ptr<ModelPart>       mp_model_part;
    BrepModelVector m_brep_model_vector;

	double m_model_tolerance;
    ///@} 
    ///@name Private Operations
    ///@{ 
    void CreateMeshedPoints(ModelPart& model_part);
    BrepFace& GetFace(const unsigned int face_id);

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


