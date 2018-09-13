//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

#if !defined(KRATOS_RIGID_BODY_ELEMENT_CREATION_UTILITY_H_INCLUDED )
#define  KRATOS_RIGID_BODY_ELEMENT_CREATION_UTILITY_H_INCLUDED


// External includes

// System includes

// Project includes
#include "includes/model_part.h"
#include "geometries/point_2d.h"
#include "geometries/point_3d.h"

#include "custom_utilities/rigid_body_utilities.hpp"

#include "custom_bounding/spatial_bounding_box.hpp"

#include "custom_elements/translatory_rigid_body_element.hpp"

#include "contact_mechanics_application_variables.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/// Rigid body element build processes in Kratos.
/**
 * Builds and element defined by its center of mass, and the properties of weight and inertia tensor
 * The RigidBodyElement is defined from a mesh of RigidBodyGeometricalElements
 * The RigidBodyProperties given by the RigidBodyBoundingBox, calculated or passed as given data to it.
 */
class KRATOS_API(CONTACT_MECHANICS_APPLICATION) RigidBodyElementCreationUtility
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION( RigidBodyElementCreationUtility );

    typedef ModelPart::NodeType                   NodeType;
    typedef ModelPart::ElementType             ElementType;
    typedef ModelPart::PropertiesType       PropertiesType;
    typedef ElementType::GeometryType         GeometryType;
    typedef Point2D<ModelPart::NodeType>       Point2DType;
    typedef Point3D<ModelPart::NodeType>       Point3DType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    RigidBodyElementCreationUtility() {}


    /// Destructor.
    virtual ~RigidBodyElementCreationUtility() {}


    ///@}
    ///@name Operators
    ///@{


private:

    //************************************************************************************
    //************************************************************************************

    void CalculateRigidBodyParameters(ModelPart& rModelPart, Vector& rCenterOfGravity, Matrix& rInertiaTensor, Matrix& rLocalAxesMatrix, double& rMass)
    {

      KRATOS_TRY

      RigidBodyUtilities RigidBodyUtils;

      rMass             =  RigidBodyUtils.MassCalculation(rModelPart);
      rCenterOfGravity  =  RigidBodyUtils.CalculateCenterOfMass(rModelPart);
      rInertiaTensor    =  RigidBodyUtils.CalculateInertiaTensor(rModelPart);


      //set inertia tensor in main axes (local inertia tensor means main axes)
      Matrix MainAxes    = ZeroMatrix(3,3);
      Matrix MainInertia = rInertiaTensor;
      RigidBodyUtils.InertiaTensorToMainAxes(MainInertia, MainAxes);

      // std::cout<<" Main Axes "<<MainAxes<<std::endl;
      // std::cout<<" Main Inertia "<<MainInertia<<std::endl;
      // std::cout<<" Inertia Tensor "<<InertiaTensor<<std::endl;

      rLocalAxesMatrix = IdentityMatrix(3);

      // main axes given in rows
      for(unsigned int i=0; i<3; i++)
    	{
    	  Vector Axis = ZeroVector(3);
    	  for(unsigned int j=0; j<3; j++)
    	    {
    	      Axis[j] = MainAxes(i,j);
    	    }

    	  double norm = norm_2(Axis);
    	  if( norm != 0)
    	    Axis/=norm;

    	  for(unsigned int j=0; j<3; j++)
    	    {
    	      rLocalAxesMatrix(j,i) = Axis[j]; //column disposition
    	    }
    	}

      rInertiaTensor = MainInertia;

      // rVolumeAcceleration = RigidBodyUtils.GetVolumeAcceleration(rModelPart);
      // rElasticModulus     = RigidBodyUtils.GetElasticModulus(rModelPart);

      std::cout<<"  [ Mass "<<rMass<<" ]"<<std::endl;
      std::cout<<"  [ CenterOfGravity "<<rCenterOfGravity<<" ]"<<std::endl;
      std::cout<<"  [ InertiaTensor "<<rInertiaTensor<<" ]"<<std::endl;

      KRATOS_CATCH("")

    }


    // //************************************************************************************
    // //************************************************************************************

    void CreateNode (NodeType::Pointer& Node, ModelPart& rModelPart, const Vector& rPoint, unsigned int& nodeId, bool& rBodyIsFixed)
    {
      KRATOS_TRY

      Node = rModelPart.CreateNewNode( nodeId, rPoint[0], rPoint[1], rPoint[2]);

      rModelPart.AddNode( Node );

      //generating the dofs
      NodeType::DofsContainerType& reference_dofs = (rModelPart.NodesBegin())->GetDofs();


      for(NodeType::DofsContainerType::iterator iii = reference_dofs.begin(); iii != reference_dofs.end(); iii++)
      	{
      	  NodeType::DofType& rDof = *iii;
      	  Node->pAddDof( rDof );
      	}


      if( rBodyIsFixed ){

    	//fix dofs:
    	NodeType::DofsContainerType& new_dofs = Node->GetDofs();

    	for(NodeType::DofsContainerType::iterator iii = new_dofs.begin(); iii != new_dofs.end(); iii++)
      	{
      	  NodeType::DofType& rDof = *iii;
    	  rDof.FixDof(); // dofs fixed
      	}

      }
      else{

    	//free dofs:
    	NodeType::DofsContainerType& new_dofs = Node->GetDofs();

    	for(NodeType::DofsContainerType::iterator iii = new_dofs.begin(); iii != new_dofs.end(); iii++)
      	{
      	  NodeType::DofType& rDof = *iii;
    	  rDof.FreeDof(); // dofs free
      	}

      }

      //generating step data
      // unsigned int buffer_size = (rModelPart.NodesBegin())->GetBufferSize();
      // unsigned int step_data_size = rModelPart.GetNodalSolutionStepDataSize();
      // for(unsigned int step = 0; step<buffer_size; step++)
      // 	{
      // 	  double* NodeData = Node->SolutionStepData().Data(step);
      // 	  double* ReferenceData = (rModelPart.NodesBegin())->SolutionStepData().Data(step);

      // 	  //copying this data in the position of the vector we are interested in
      // 	  for(unsigned int j= 0; j<step_data_size; j++)
      // 	    {
      // 	      NodeData[j] = ReferenceData[j];
      // 	    }
      // 	}

      KRATOS_CATCH("")

    }


    // ///@}
    // ///@name Operations
    // ///@{

 public:


    void CreateRigidBodyElement(ModelPart& rMainModelPart,
    			        SpatialBoundingBox::Pointer pRigidBodyBox,
    				Parameters CustomParameters)
    {

      KRATOS_TRY

      Parameters DefaultParameters( R"(
            {
                "rigid_body_element_type": "TranslatoryRigidElement3D1N",
                "fixed_body": true,
                "compute_body_parameters": false,
                "rigid_body_model_part_name": "RigidBodyDomain",
                "rigid_body_parameters":{
                    "center_of_gravity": [0.0 ,0.0, 0.0],
                    "mass": 0.0,
                    "main_inertias": [0.0, 0.0, 0.0],
                    "main_axes": [ [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0] ]
                }
            }  )" );


      //validate against defaults -- this also ensures no type mismatch
      CustomParameters.ValidateAndAssignDefaults(DefaultParameters);

      bool BodyIsFixed = CustomParameters["fixed_body"].GetBool();

      //create properties for the rigid body
      unsigned int NumberOfProperties = rMainModelPart.NumberOfProperties();

      PropertiesType::Pointer pProperties = Kratos::make_shared<PropertiesType>(NumberOfProperties);

      double Mass = 0;
      Vector CenterOfGravity   = ZeroVector(3);
      Matrix InertiaTensor     = ZeroMatrix(3);
      Matrix LocalAxesMatrix   = IdentityMatrix(3);

      bool ComputeBodyParameters = CustomParameters["compute_body_parameters"].GetBool();
      ModelPart& rRigidBodyModelPart = rMainModelPart.GetSubModelPart(CustomParameters["rigid_body_model_part_name"].GetString());

      if( ComputeBodyParameters ){

    	this->CalculateRigidBodyParameters( rRigidBodyModelPart, CenterOfGravity, InertiaTensor, LocalAxesMatrix, Mass );
      }
      else{

    	Parameters RigidBodyProperties = CustomParameters["rigid_body_parameters"];

    	Mass = RigidBodyProperties["mass"].GetDouble();

    	unsigned int size = RigidBodyProperties["main_inertias"].size();

    	for( unsigned int i=0; i<size; i++ )
    	  {
    	    Parameters LocalAxesRow = RigidBodyProperties["main_axes"][i];

    	    CenterOfGravity[i]     = RigidBodyProperties["center_of_gravity"][i].GetDouble();
    	    InertiaTensor(i,i)     = RigidBodyProperties["main_inertias"][i].GetDouble();

    	    LocalAxesMatrix(0,i)   = LocalAxesRow[0].GetDouble(); //column disposition
    	    LocalAxesMatrix(1,i)   = LocalAxesRow[1].GetDouble();
    	    LocalAxesMatrix(2,i)   = LocalAxesRow[2].GetDouble();
    	  }

	std::cout<<"  [ Mass "<<Mass<<" ]"<<std::endl;
	std::cout<<"  [ CenterOfGravity "<<CenterOfGravity<<" ]"<<std::endl;
	std::cout<<"  [ InertiaTensor "<<InertiaTensor<<" ]"<<std::endl;


      }

      pProperties->SetValue(NODAL_MASS, Mass);
      pProperties->SetValue(LOCAL_INERTIA_TENSOR, InertiaTensor);
      pProperties->SetValue(LOCAL_AXES_MATRIX, LocalAxesMatrix);

      //add properties to model part
      rMainModelPart.AddProperties(pProperties);

      // create node for the rigid body center of gravity:
      unsigned int LastNodeId  = rMainModelPart.Nodes().back().Id() + 1;

      NodeType::Pointer NodeCenterOfGravity;
      this->CreateNode( NodeCenterOfGravity, rMainModelPart, CenterOfGravity, LastNodeId, BodyIsFixed);

      //Set this node to the boundary model_part where it belongs to
      unsigned int RigidBodyNodeId = rRigidBodyModelPart.Nodes().back().Id();

      for(ModelPart::SubModelPartIterator i_mp= rMainModelPart.SubModelPartsBegin(); i_mp!=rMainModelPart.SubModelPartsEnd(); i_mp++)
	{

	  if(i_mp->Is(BOUNDARY)){
	    std::cout<<" boundary model part "<<i_mp->Name()<<std::endl;
	    for(ModelPart::NodesContainerType::iterator i_node = i_mp->NodesBegin(); i_node!= i_mp->NodesEnd(); i_node++)
	      {
		if( i_node->Id() == RigidBodyNodeId ){
		  i_mp->AddNode(NodeCenterOfGravity);
		  std::cout<<" node set "<<std::endl;
		  break;
		}
	      }
	  }

	}

      // set node variables
      NodeCenterOfGravity->GetSolutionStepValue(VOLUME_ACCELERATION) = rRigidBodyModelPart.Nodes().back().GetSolutionStepValue(VOLUME_ACCELERATION);

      // set node flags
      NodeCenterOfGravity->Set(MASTER);
      NodeCenterOfGravity->Set(RIGID);

      // set node to the spatial bounding box
      pRigidBodyBox->SetRigidBodyCenter(NodeCenterOfGravity);


      // create rigid body element:
      unsigned int LastElementId = rMainModelPart.Elements().back().Id() + 1;

      std::string ElementName = CustomParameters["rigid_body_element_type"].GetString();

      // always a 3D point
      GeometryType::Pointer pGeometry = Kratos::make_shared<Point3DType>(NodeCenterOfGravity);

      ModelPart::NodesContainerType::Pointer pNodes =  rRigidBodyModelPart.pNodes();

      ElementType::Pointer pRigidBodyElement;

      //Rigid Body Element:
      if( ElementName == "RigidBodyElement3D1N" || ElementName == "RigidBodyElement2D1N" ){
	//std::cout<<" RigidBodyElement "<<std::endl;
    	pRigidBodyElement = Kratos::make_shared<RigidBodyElement>(LastElementId, pGeometry, pProperties, pNodes);
      }
      else if( ElementName == "TranslatoryRigidBodyElement3D1N" || ElementName == "TranslatoryRigidBodyElement2D1N"){
	//std::cout<<" TranslatoryRigidBodyElement "<<std::endl;
	//ElementType::Pointer pElem = KratosComponents<Element>::Get("TranslatoryRigidBodyElement")
    	pRigidBodyElement = Kratos::make_shared<TranslatoryRigidBodyElement>(LastElementId, pGeometry, pProperties, pNodes);
      }
      else if( ElementName == "RigidBodyEMCElement3D1N" || ElementName == "RigidBodyEMCElement2D1N" ){
	//std::cout<<" RigidBodyEMCElement "<<std::endl;
    	//pRigidBodyElement = Kratos::make_shared<RigidBodyEMCElement>(LastElementId, pGeometry, pProperties, pNodes);
      }

      // once conventional constructor and registered
      // ElementType::NodesArrayType ElementNodes;
      // ElementNodes.push_back(NodeCenterOfGravity);
      // ElementType const& rCloneElement = KratosComponents<ElementType>::Get(ElementName);
      // ElementType::Pointer pRigidBodyElement = rCloneElement.Create(LastElementId, ElementNodes, pProperties);
      // rRigidBodyModelPart.AddElement(pRigidBodyElement);

      // other posibility
      // std::vector<int> NodeIds;
      // NodeIds.push_back(LastNodeId);
      // rRigidBodyModelPart.CreateNewElement(ElementName,LastElementId, NodeIds, pProperties);

      rRigidBodyModelPart.AddElement(pRigidBodyElement);
      rRigidBodyModelPart.AddNode(NodeCenterOfGravity);

      //add rigid body element to computing model part:
      for(ModelPart::SubModelPartIterator i_mp= rMainModelPart.SubModelPartsBegin() ; i_mp!=rMainModelPart.SubModelPartsEnd(); i_mp++)
	{
	  if( (i_mp->Is(ACTIVE)) ){ //computing_domain
	    pRigidBodyElement->Set(ACTIVE,true);
	    rMainModelPart.GetSubModelPart(i_mp->Name()).AddElement(pRigidBodyElement);
	    rMainModelPart.GetSubModelPart(i_mp->Name()).AddNode(NodeCenterOfGravity);
	  }
	}

      std::cout<<"  [ "<<ElementName<<" Created : [NodeId:"<<LastNodeId<<"] [ElementId:"<<LastElementId<<"] CG("<<NodeCenterOfGravity->X()<<","<<NodeCenterOfGravity->Y()<<","<<NodeCenterOfGravity->Z()<<") ]"<<std::endl;


      KRATOS_CATCH( "" )

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
        return "RigidBodyElementCreationUtility";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "RigidBodyElementCreationUtility";
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

    /// Assignment operator.
    RigidBodyElementCreationUtility& operator=(RigidBodyElementCreationUtility const& rOther);

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
    ///@name Unaccessible methods
    ///@{

    ///@}

    ///@}

}; // Class Process

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  RigidBodyElementCreationUtility& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const RigidBodyElementCreationUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_RIGID_BODY_ELEMENT_CREATION_UTILITY_H_INCLUDED  defined
