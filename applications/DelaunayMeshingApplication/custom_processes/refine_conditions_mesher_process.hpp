//
//   Project Name:        KratosDelaunayMeshingApplication $
//   Created by:          $Author:             JMCarbonell $
//   Last modified by:    $Co-Author:                      $
//   Date:                $Date:                April 2018 $
//   Revision:            $Revision:                   0.0 $
//
//

#if !defined(KRATOS_REFINE_CONDITIONS_MESHER_PROCESS_H_INCLUDED )
#define  KRATOS_REFINE_CONDITIONS_MESHER_PROCESS_H_INCLUDED


// External includes

// System includes

// Project includes
#include "includes/model_part.h"
#include "custom_utilities/mesher_utilities.hpp"
#include "custom_processes/mesher_process.hpp"

///VARIABLES used:
//StepData: NODAL_H, NORMAL, CONTACT_FORCE, DISPLACEMENT
//Flags:    (checked) BOUNDARY,
//          (set)     BOUNDARY(nodes), TO_ERASE(conditions), NEW_ENTITY(conditions,nodes)(set),
//          (modified)
//          (reset)
//(set):=(set in this process)

namespace Kratos
{

///@name Kratos Classes
///@{

/// Refine Mesh Boundary Process
/** The process labels the boundary conditions
    Dependencies: RemoveMeshNodesProcess.Execute()  is needed as a previous step

    Determines if new conditions must be inserted in boundary.
    If boundary must to be kept (CONSTRAINED),
    New conditions will be rebuild (splitting the old ones and inserting new nodes)
    Old conditions will be erased at the end.

*/

class RefineConditionsMesherProcess
  : public MesherProcess
 {
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION( RefineConditionsMesherProcess );

    typedef ModelPart::NodeType                   NodeType;
    typedef ModelPart::ConditionType         ConditionType;
    typedef ModelPart::PropertiesType       PropertiesType;
    typedef ConditionType::GeometryType       GeometryType;
    typedef PointerVector<NodeType>        PointsArrayType;

    typedef PointerVectorSet<ConditionType, IndexedObject> ConditionsContainerType;
    typedef ConditionsContainerType::iterator                    ConditionIterator;
    typedef ConditionsContainerType::const_iterator      ConditionConstantIterator;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    RefineConditionsMesherProcess(ModelPart& rModelPart,
			      MesherUtilities::MeshingParameters& rRemeshingParameters,
			      int EchoLevel)
      : mrModelPart(rModelPart),
	mrRemesh(rRemeshingParameters)
    {
      mEchoLevel = EchoLevel;
    }

    /// Destructor.
    virtual ~RefineConditionsMesherProcess() {}


    ///@}
    ///@name Operators
    ///@{

    /// This operator is provided to call the process as a function and simply calls the Execute method.
    void operator()()
    {
        Execute();
    }


    ///@}
    ///@name Operations
    ///@{


    /// Execute method is used to execute the Process algorithms.
    void Execute()  override
    {

      KRATOS_TRY

      if( this->mEchoLevel > 0 ){
        std::cout<<" [ REFINE BOUNDARY : "<<std::endl;
	//std::cout<<"   Nodes and Conditions : "<<mrModelPart.Nodes().size()<<", "<<mrModelPart.Conditions().size()<<std::endl;
      }


      if( mrModelPart.Name() != mrRemesh.SubModelPartName )
	std::cout<<" ModelPart Supplied do not corresponds to the Meshing Domain: ("<<mrModelPart.Name()<<" != "<<mrRemesh.SubModelPartName<<")"<<std::endl;


      mrRemesh.Info->InsertedBoundaryConditions = mrModelPart.NumberOfConditions();
      mrRemesh.Info->InsertedBoundaryNodes = mrModelPart.NumberOfNodes();


     //if the insert switches are activated, we check if the boundaries got too coarse
     if( (mrRemesh.Refine->RefiningOptions.Is(MesherUtilities::REFINE_INSERT_NODES) || mrRemesh.Refine->RefiningOptions.Is(MesherUtilities::REFINE_ADD_NODES)) && mrRemesh.Refine->RefiningOptions.Is(MesherUtilities::REFINE_BOUNDARY) )
     {

        std::vector<NodeType::Pointer>  list_of_nodes;
        std::vector<ConditionType::Pointer> list_of_conditions;

	unsigned int conditions_size = mrModelPart.NumberOfConditions();

	list_of_nodes.reserve(conditions_size);
	list_of_conditions.reserve(conditions_size);

	this->SelectBoundaryToRefine(mrModelPart); //conditions (TO_REFINE)

	this->GenerateNewNodes(mrModelPart, list_of_nodes, list_of_conditions); //points (NEW_ENTITY)

	this->GenerateNewConditions(mrModelPart, list_of_nodes, list_of_conditions);// new conditions(NEW_ENTITY)  //old conditions (TO_ERASE)

	this->SetNodesToModelPart(mrModelPart, list_of_nodes);

	this->SetConditionsToModelPart(mrModelPart, list_of_conditions);

	//new conditions are added to model part and later added to condition model parts via (NEW_ENTITY) and (MODEL_PART_NAME)


     } // REFINE END;


     mrRemesh.Info->InsertedBoundaryNodes = mrModelPart.NumberOfNodes()-mrRemesh.Info->InsertedBoundaryNodes;

     if( this->mEchoLevel > 0 ){
        std::cout<<"   [ CONDITIONS ( total : "<<mrModelPart.NumberOfConditions()<<" ) ]"<<std::endl;
        std::cout<<"   [ NODES      ( inserted : "<<mrRemesh.Info->InsertedBoundaryNodes<<" total: "<<mrModelPart.NumberOfNodes()<<" ) ]"<<std::endl;

	if( this->mEchoLevel >=1 ){
	  mrRemesh.Refine->Info.BoundaryConditionsRefined.EchoStats();
	}

        std::cout<<"   REFINE BOUNDARY ]; "<<std::endl;
     }

     KRATOS_CATCH(" ")

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
    std::string Info() const override
    {
        return "RefineConditionsMesherProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "RefineConditionsMesherProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
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

    ModelPart& mrModelPart;

    MesherUtilities::MeshingParameters& mrRemesh;

    MesherUtilities mMesherUtilities;

    int mEchoLevel;

    ///@}
    ///@name Protected Operators
    ///@{

    bool RefineOnThreshold(ConditionType::Pointer& pCondition, ProcessInfo& rCurrentProcessInfo, double& critical_size)
    {
      KRATOS_TRY

      if( pCondition->GetValue(MASTER_ELEMENTS).size() > 0 ){

	Element::ElementType& MasterElement = pCondition->GetValue(MASTER_ELEMENTS).back();

	std::vector<double> Value;

	MasterElement.GetValueOnIntegrationPoints(mrRemesh.Refine->GetThresholdVariable(),Value,rCurrentProcessInfo);

	//calculate threshold value (plastic power)
	double threshold_value = 0;

	for(std::vector<double>::iterator v = Value.begin(); v!=Value.end(); ++v)
	  threshold_value += *v;

	threshold_value /= double(Value.size());
	threshold_value *= MasterElement.GetGeometry().DomainSize();

	//calculate condition length
	double face_size = mMesherUtilities.CalculateBoundarySize(pCondition->GetGeometry());

	if( threshold_value > mrRemesh.Refine->ReferenceThreshold && face_size > critical_size )
	  return true;
      }

      return false;

      KRATOS_CATCH( "" )
    }

    //*******************************************************************************************
    //*******************************************************************************************


    bool RefineOnDistance(ConditionType::Pointer& pCondition, double& critical_size)
    {
      KRATOS_TRY

      //calculate condition length
      double face_size = mMesherUtilities.CalculateBoundarySize(pCondition->GetGeometry());

      if( face_size > critical_size )
	return true;

      return false;

      KRATOS_CATCH( "" )
    }

    //*******************************************************************************************
    //*******************************************************************************************

    bool RefineBoundaryCondition(ConditionType::Pointer& pCondition, ProcessInfo& rCurrentProcessInfo)
    {
      KRATOS_TRY

      bool refine_condition = false;

      //THRESHOLD VALUE INSERT
      double size_for_threshold_face  = 2.50 * mrRemesh.Refine->CriticalSide;

      if ( mrRemesh.Refine->RefiningOptions.Is(MesherUtilities::REFINE_BOUNDARY_ON_THRESHOLD) )
	refine_condition = this->RefineOnThreshold(pCondition, rCurrentProcessInfo, size_for_threshold_face);

      if( refine_condition ){

	//check a the critical distance to not refine without limits
	double size_for_boundary_threshold = mrRemesh.Refine->CriticalSide;
	refine_condition = this->RefineOnDistance(pCondition, size_for_boundary_threshold);

	if( refine_condition ){
	  mrRemesh.Refine->Info.BoundaryConditionsRefined.on_threshold++;
	  return true;
	}

      }

      refine_condition = false;

      //DISTANCE VALUE INSERT
      double size_for_boundary_face   = 3.50 * mrRemesh.Refine->CriticalSide;

      if ( mrRemesh.Refine->RefiningOptions.Is(MesherUtilities::REFINE_BOUNDARY_ON_DISTANCE) )
	refine_condition = this->RefineOnDistance(pCondition, size_for_boundary_face);

      if( refine_condition ){
	mrRemesh.Refine->Info.BoundaryConditionsRefined.on_distance++;
	return true;
      }

      return false;

      KRATOS_CATCH( "" )
    }


    //*******************************************************************************************
    //*******************************************************************************************

    bool RefineContactCondition(ConditionType::Pointer& pCondition, ProcessInfo& rCurrentProcessInfo)
    {
      KRATOS_TRY

      if ( mrRemesh.Refine->RefiningOptions.Is(MesherUtilities::REFINE_BOUNDARY_ON_DISTANCE)
	   || ( mrRemesh.Refine->RefiningOptions.Is(MesherUtilities::REFINE_BOUNDARY_ON_THRESHOLD) ) ){

	bool refine_condition    = false;
	bool curved_contact      = false;
	bool contact_semi_active = false;
	std::vector<bool> semi_active_nodes;

	bool contact_active = mMesherUtilities.CheckContactActive(pCondition->GetGeometry(), contact_semi_active, semi_active_nodes);

	double factor = 3.50;

	if( contact_semi_active ){

	  std::vector<array_1d<double,3> > contact_normals;

	  curved_contact = mMesherUtilities.CheckContactCurvature(pCondition->GetGeometry(), contact_normals);

	  //FACTOR VALUE INSERT plane contact transition
	  factor = 2.0;

	  //FACTOR VALUE INSERT curved contact transition
	  if( curved_contact )
	    factor = 0.75;


	  if( contact_active ){

	    //FACTOR VALUE INSERT plane contact
	    factor = 1.5;

	    //FACTOR VALUE INSERT curved contact
	    if( curved_contact )
	      factor = 0.5;

	  }

	}

	if( contact_active || contact_semi_active ){

	  double size_for_boundary_contact_face  = factor * mrRemesh.Refine->CriticalSide;
	  refine_condition = this->RefineOnDistance(pCondition, size_for_boundary_contact_face);

	  if( refine_condition ){

	    mrRemesh.Refine->Info.BoundaryConditionsRefined.on_distance++;
	    if(contact_active || contact_semi_active){
	      mrRemesh.Refine->Info.BoundaryConditionsRefined.in_contact++;
	      if(curved_contact)
		mrRemesh.Refine->Info.BoundaryConditionsRefined.in_concave_boundary++;
	    }
	    return true;
	  }

	}

      }

      return false;

      KRATOS_CATCH( "" )
    }


    //*******************************************************************************************
    //*******************************************************************************************

    void SetNodalPosition(ConditionType::Pointer& pCondition, ProcessInfo& rCurrentProcessInfo, double& xc, double& yc, double& zc)
    {
      KRATOS_TRY

      bool contact_semi_active = false;
      std::vector<bool> semi_active_nodes;

      bool contact_active = mMesherUtilities.CheckContactActive(pCondition->GetGeometry(), contact_semi_active, semi_active_nodes);

      if( contact_semi_active || contact_active ){ // if contact is semi_ative or active

	array_1d<double,3> position_correction;
	position_correction.clear();

	if( pCondition->GetGeometry().size() == 2 && pCondition->GetGeometry().WorkingSpaceDimension() == 2 ){ //line

	  //circle interpolation
	  //this->CircleInterpolation(pCondition, position_correction);

	  //hermite interpolation
	  this->HermiteInterpolation(pCondition, position_correction);
	}


	if( pCondition->GetGeometry().size() == 3 && pCondition->GetGeometry().WorkingSpaceDimension() == 3 ){ //triangle

	  //hermite interpolation
	  this->HermiteTriangleInterpolation(pCondition, position_correction);

	}

	//correct only if not self contact
	MesherUtilities MesherUtils;
	if( !MesherUtils.CheckSubdomain(pCondition->GetGeometry()) )
	  {
	    xc += position_correction[0];
	    yc += position_correction[1];
	    zc += position_correction[2];
	  }

      }

      KRATOS_CATCH( "" )
    }


    //*******************************************************************************************
    //*******************************************************************************************

    void CircleInterpolation(ConditionType::Pointer& pCondition, array_1d<double,3 >& rVector)
    {
      KRATOS_TRY

      bool is_curved = false;
      std::vector<array_1d<double, 3> > normals;

      is_curved = mMesherUtilities.CheckContactCurvature(pCondition->GetGeometry(), normals);

      array_1d<double,3> normal_direction;
      normal_direction.clear();

      array_1d<double,3> tangent_direction;
      tangent_direction.clear();

      if( is_curved ){

	//compute the saggita
	double projection = 0.0;
	for( unsigned int i = 0; i<3; i++ )
	  projection += normals[0][i] * normals[1][i];

	projection = std::sqrt(projection);

	double angle = std::acos(projection);

	double face_size = mMesherUtilities.CalculateBoundarySize(pCondition->GetGeometry());

	double sagitta = 0.5 * face_size * std::tan(0.25*angle);

	//correction vector according to contact curvature
	normal_direction = normals[0] +  normals[1];

	double modulus = norm_2(normal_direction);
	if( modulus )
	  normal_direction /= modulus;

	normal_direction *= sagitta;

	//check correct curvature convexity
	tangent_direction = pCondition->GetGeometry()[0]-pCondition->GetGeometry()[1];
	modulus = norm_2(tangent_direction);
	if( modulus )
	  tangent_direction /= modulus;

	//note:  two possible directions depending on the NORMAL_CONTACT definition
	if( inner_prod( normals[0], tangent_direction) > 0 )
	  normal_direction *= (-1.0);
      }

      rVector = normal_direction;

      KRATOS_CATCH( "" )
    }


    //*******************************************************************************************
    //*******************************************************************************************

    void HermiteInterpolation(ConditionType::Pointer& pCondition, array_1d<double,3 >& rVector)
    {
      KRATOS_TRY

      bool is_curved = false;
      std::vector<array_1d<double, 3> > normals;

      is_curved = mMesherUtilities.CheckContactCurvature(pCondition->GetGeometry(), normals);

      if( is_curved ){

	this->HermiteInterpolation( pCondition->GetGeometry()[0], pCondition->GetGeometry()[1], normals[0], normals[1], rVector, 0.5);
      }



      KRATOS_CATCH( "" )
    }



    //*******************************************************************************************
    //*******************************************************************************************

    void HermiteTriangleInterpolation(ConditionType::Pointer& pCondition, array_1d<double,3 >& rVector)
    {
      KRATOS_TRY

      bool is_curved = false;
      std::vector<array_1d<double, 3> > normals;

      is_curved = mMesherUtilities.CheckContactCurvature(pCondition->GetGeometry(), normals);

      double baricenter = 2.0/3.0;
      array_1d<double, 3> MidPoint;
      MidPoint.clear();

      array_1d<double, 3> Normal;
      Normal.clear();

      array_1d<double, 3> Curve;
      Curve.clear();

      if( is_curved ){

	//first curve (node to face midpoint)
	MidPoint = 0.5 * (pCondition->GetGeometry()[2] - pCondition->GetGeometry()[1]);
	Normal   = 0.5 * (normals[2] + normals[1]);

	this->HermiteInterpolation( pCondition->GetGeometry()[0], MidPoint, normals[0], Normal, Curve, baricenter);

	rVector = Curve;

	//second curve (node to face midpoint)
	MidPoint = 0.5 * (pCondition->GetGeometry()[2] - pCondition->GetGeometry()[0]);
	Normal   = 0.5 * (normals[2] + normals[0]);

	Curve.clear();
	this->HermiteInterpolation( pCondition->GetGeometry()[1], MidPoint, normals[1], Normal, Curve, baricenter);

	rVector += Curve;

	//first curve (node to face midpoint)
	MidPoint = 0.5 * (pCondition->GetGeometry()[1] - pCondition->GetGeometry()[0]);
	Normal   = 0.5 * (normals[1] + normals[0]);

	this->HermiteInterpolation( pCondition->GetGeometry()[2], MidPoint, normals[2], Normal, Curve, baricenter);

	rVector += Curve;

	rVector *= 1.0/3.0;
      }


      KRATOS_CATCH( "" )
    }


    //*******************************************************************************************
    //*******************************************************************************************

    void HermiteInterpolation(const array_1d<double,3 >& rP1, const array_1d<double,3 >& rP2, const array_1d<double,3 >& rN1, const array_1d<double,3 >& rN2, array_1d<double,3 >& rD, double s)
    {

       KRATOS_TRY

       //compute points distance 1-2
       array_1d<double, 3> T1 = rP2 - rP1;

       //compute tangents
       double projection = 0.0;
       for( unsigned int i = 0; i<3; i++ )
	 projection += T1[i] * rN1[i];

       T1  -= ( projection * rN1 );

       double modulus = norm_2(T1);
       if( modulus )
	 T1 /= modulus;

       //compute points distance 2-1
       array_1d<double, 3> T2 = rP1 - rP2;

       //compute tangents
       projection = 0.0;
       for( unsigned int i = 0; i<3; i++ )
	 projection += T2[i] * rN2[i];

       T2 -= ( projection * rN2 );

       modulus = norm_2(T2);
       if( modulus )
	 T2 /= modulus;

       T2 *= (-1);

       //compute normalized s-point position s in [0,1]
       array_1d<double, 3> M = s * rP2 - s * rP1;

       modulus  = norm_2(rP2-rP1);
       if( modulus )
	 projection = norm_2(M)/modulus;

       //hermite basis functions
       double h00 = 2.0 * projection * projection * projection - 3.0 * projection * projection + 1.0;
       double h10 = projection * projection * projection - 2.0 * projection * projection + projection;
       double h01 = -2.0 * projection * projection * projection + 3.0 * projection * projection;
       double h11 = projection * projection * projection - projection * projection;

       //hermite interpolation polinomial
       rD  = h00 * rP1;
       rD += h10 * modulus * T1;
       rD += h01 * rP2;
       rD += h11 * modulus * T2;

       //increment of position
       M = (s * rP2 + (1-s) * rP1);

       rD -= M;

       //check correct curvature convexity
       T1 = rP2 - rP1;
       modulus = norm_2(T1);
       if( modulus )
	 T1 /= modulus;

       //note:  two possible directions depending on the NORMAL_CONTACT definition
       double d1 = inner_prod( rN1, T1);
       double d2 = inner_prod( rN2, -T1);
       if( d1 * d2 > 0 ){
	 if( d1 > 0 )
	   rD *= (-1.0);
       }
       else{
	 rD *= 0.0;
       }

       KRATOS_CATCH( "" )

    }

    //*******************************************************************************************
    //*******************************************************************************************

    void SelectBoundaryToRefine(ModelPart& rModelPart)
    {

      KRATOS_TRY

      ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();
      bool refine_condition = false;


      mrRemesh.Refine->Info.BoundaryConditionsRefined.Initialize();

      //LOOP TO CONSIDER ALL SUBDOMAIN CONDITIONS

      bool refine_candidate = false;
      for(ModelPart::ConditionsContainerType::iterator i_cond = rModelPart.ConditionsBegin(); i_cond!= rModelPart.ConditionsEnd(); ++i_cond)
	{

	  refine_candidate = false;
	  i_cond->Set(TO_REFINE, false);

	  if( mrRemesh.Options.Is(MesherUtilities::CONSTRAINED) ){
	    if( i_cond->Is(BOUNDARY) ) //ONLY SET TO THE BOUNDARY SKIN CONDITIONS (CompositeCondition)
	      refine_candidate = true;
	    else
	      refine_candidate = false;
	  }
	  else{
	    refine_candidate = true;
	  }


	  if( refine_candidate ){
	    if (mrRemesh.Refine->RefiningBoxSetFlag == true ){
	      refine_candidate = mMesherUtilities.CheckConditionInBox(*(i_cond.base()), *(mrRemesh.Refine->RefiningBox), rCurrentProcessInfo);
	    }
	  }


	  if( refine_candidate ){

	    //double condition_radius = 0;
	    if( i_cond->IsNot(TO_ERASE) ){

	      refine_condition = this->RefineBoundaryCondition(*(i_cond.base()), rCurrentProcessInfo);

	      if( refine_condition ){
		i_cond->Set(TO_REFINE);
	      }
	      else{
		refine_condition = this->RefineContactCondition(*(i_cond.base()), rCurrentProcessInfo);

		if( refine_condition )
		  i_cond->Set(TO_REFINE);
	      }

	    }
	    else{
	      if( this->mEchoLevel > 0 )
		std::cout<<" Condition "<<i_cond->Id()<<" TO_ERASE "<<std::endl;
	    }

	  }
        }


      KRATOS_CATCH( "" )
    }

    //*******************************************************************************************
    //*******************************************************************************************

    void GenerateNewNodes(ModelPart& rModelPart, std::vector<NodeType::Pointer>& list_of_nodes, std::vector<ConditionType::Pointer>& list_of_conditions)
    {
      KRATOS_TRY

      ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

      MeshDataTransferUtilities DataTransferUtilities;

      NodeType::Pointer pNode;

      //center
      double xc = 0;
      double yc = 0;
      double zc = 0;

      //radius
      double radius = 0;

      //assign data to dofs
      NodeType::DofsContainerType& ReferenceDofs = rModelPart.Nodes().front().GetDofs();

      VariablesList& VariablesList = rModelPart.GetNodalSolutionStepVariablesList();


      std::vector<double> ShapeFunctionsN;

      unsigned int id = MesherUtilities::GetMaxNodeId(rModelPart) + 1;

      unsigned int size  = 0;
      unsigned int count = 0;

      for(ModelPart::ConditionsContainerType::iterator i_cond = rModelPart.ConditionsBegin(); i_cond!= rModelPart.ConditionsEnd(); ++i_cond)
	{
	  if( i_cond->Is(TO_REFINE) )
	    {

	      Geometry< Node<3> >& rGeometry = i_cond->GetGeometry();

	      size = rGeometry.size();

	      ShapeFunctionsN.resize(size);


	      if( size == 2 )
		DataTransferUtilities.CalculateCenterAndSearchRadius( rGeometry[0].X(), rGeometry[0].Y(),
								      rGeometry[1].X(), rGeometry[1].Y(),
								      xc,yc,radius);


	      if( size == 3 )
		DataTransferUtilities.CalculateCenterAndSearchRadius( rGeometry[0].X(), rGeometry[0].Y(), rGeometry[0].Z(),
								      rGeometry[1].X(), rGeometry[1].Y(), rGeometry[1].Z(),
								      rGeometry[2].X(), rGeometry[2].Y(), rGeometry[2].Z(),
								      xc,yc,zc,radius);


	      this->SetNodalPosition(*(i_cond.base()), rCurrentProcessInfo, xc, yc, zc);

	      //create a new node
	      pNode = Kratos::make_shared< NodeType >( id, xc, yc, zc );

	      //giving model part variables list to the node
	      pNode->SetSolutionStepVariablesList(&VariablesList);

	      //set buffer size
	      pNode->SetBufferSize(rModelPart.GetBufferSize());

	      //generating the dofs
	      for(Node<3>::DofsContainerType::iterator i_dof = ReferenceDofs.begin(); i_dof != ReferenceDofs.end(); ++i_dof)
		{
		  NodeType::DofType& rDof = *i_dof;
		  NodeType::DofType::Pointer pNewDof = pNode->pAddDof( rDof );

		  count = 0;
		  for( unsigned int i = 0; i<size; i++ )
		    {
		      if(rGeometry[i].IsFixed(rDof.GetVariable()))
			count++;
		    }

		  if( count == size )
		    (pNewDof)->FixDof();
		  else
		    (pNewDof)->FreeDof();
		}

	      std::fill(ShapeFunctionsN.begin(), ShapeFunctionsN.end(), 1.0/double(size));

	      double alpha = 1;
	      DataTransferUtilities.Interpolate( rGeometry, ShapeFunctionsN, VariablesList, pNode, alpha );

	      //set flags
	      pNode->Set(NEW_ENTITY);
	      pNode->Set(BOUNDARY);

	      //set variables
	      this->SetNewNodeVariables(rModelPart, *(i_cond.base()), pNode);

	      list_of_nodes.push_back(pNode);
	      list_of_conditions.push_back(*(i_cond.base()));

	      id++;

	    }

	}

      KRATOS_CATCH( "" )
    }


    //*******************************************************************************************
    //*******************************************************************************************

    virtual void SetNewNodeVariables(ModelPart& rModelPart, ConditionType::Pointer& pCondition, NodeType::Pointer& pNode)
    {
      KRATOS_TRY

      //set variables:
      Geometry< Node<3> >& rGeometry = pCondition->GetGeometry();


      //set model part
      pNode->SetValue(MODEL_PART_NAME,rModelPart.Name());

      //set nodal_h
      //pNode->FastGetSolutionStepValue(NODAL_H) = mrRemesh.Refine->CriticalSide; //too small problems
      pNode->FastGetSolutionStepValue(NODAL_H) = rGeometry.DomainSize();

      //set normal
      noalias(pNode->FastGetSolutionStepValue(NORMAL)) = pCondition->GetValue(NORMAL);

      //set original position
      const array_1d<double,3>& Displacement = pNode->FastGetSolutionStepValue(DISPLACEMENT);
      pNode->X0() = pNode->X() - Displacement[0];
      pNode->Y0() = pNode->Y() - Displacement[1];
      pNode->Z0() = pNode->Z() - Displacement[2];

      //set contact force
      unsigned int count = 0;
      for( unsigned int i = 0; i<rGeometry.size(); i++ )
	{
	  if( norm_2(rGeometry[i].FastGetSolutionStepValue(CONTACT_FORCE)) == 0 )
	    count++;
	}

      if( count )
	pNode->FastGetSolutionStepValue(CONTACT_FORCE).clear();


      KRATOS_CATCH( "" )
    }


    //*******************************************************************************************
    //*******************************************************************************************

    void GenerateNewConditions(ModelPart& rModelPart, std::vector<NodeType::Pointer >& list_of_nodes, std::vector<ConditionType::Pointer>& list_of_conditions)
    {
      KRATOS_TRY

      std::vector<ConditionType::Pointer> list_of_new_conditions;

      unsigned int id = MesherUtilities::GetMaxConditionId(rModelPart) + 1;

      ConditionType::Pointer pCondition;

      int size = 0;

      unsigned int counter = 0;

      for(std::vector<ConditionType::Pointer>::iterator i_cond = list_of_conditions.begin(); i_cond!= list_of_conditions.end(); ++i_cond)
	{
	  Geometry< Node<3> >& rGeometry = (*i_cond)->GetGeometry();

	  size = rGeometry.size();

	  PointsArrayType Nodes(size);

	  if( size == 2 ){

	    //new condition 1
	    Nodes(0) = rGeometry(0);
	    Nodes(1) = list_of_nodes[counter];

	    pCondition = (*i_cond)->Clone(id, Nodes);

	    //set flags
	    pCondition->Set(NEW_ENTITY);

	    SetNewConditionVariables((*i_cond), pCondition);

	    id++;

	    list_of_new_conditions.push_back(pCondition);

	    //new condition 2
	    Nodes(0) = list_of_nodes[counter];
	    Nodes(1) = rGeometry(1);

	    pCondition = (*i_cond)->Clone(id, Nodes);

	    //set flags
	    pCondition->Set(NEW_ENTITY);

	    //set variables
	    this->SetNewConditionVariables((*i_cond), pCondition);

	    id++;

	    list_of_new_conditions.push_back(pCondition);
	  }

	  if( size == 3 ){

	    //new condition 1
	    Nodes(0) = rGeometry(0);
	    Nodes(1) = rGeometry(1);
	    Nodes(2) = list_of_nodes[counter];

	    pCondition = (*i_cond)->Clone(id, Nodes);

	    //set flags
	    pCondition->Set(NEW_ENTITY);

	    SetNewConditionVariables((*i_cond), pCondition);

	    id++;

	    list_of_new_conditions.push_back(pCondition);

	    //new condition 2
	    Nodes(0) = rGeometry(1);
	    Nodes(1) = rGeometry(2);
	    Nodes(2) = list_of_nodes[counter];


	    pCondition = (*i_cond)->Clone(id, Nodes);

	    //set flags
	    pCondition->Set(NEW_ENTITY);

	    //set variables
	    SetNewConditionVariables((*i_cond), pCondition);

	    id++;

	    list_of_new_conditions.push_back(pCondition);


	    //new condition 3
	    Nodes(0) = rGeometry(2);
	    Nodes(1) = rGeometry(0);
	    Nodes(2) = list_of_nodes[counter];


	    pCondition = (*i_cond)->Clone(id, Nodes);

	    //set flags
	    pCondition->Set(NEW_ENTITY);

	    //set variables
	    SetNewConditionVariables((*i_cond), pCondition);

	    id++;

	    list_of_new_conditions.push_back(pCondition);


	  }

	  // once the condition is refined set to erase
	  (*i_cond)->Set(TO_ERASE);
	  (*i_cond)->Set(TO_REFINE, false);

	  ConditionsContainerType& ChildrenConditions = (*i_cond)->GetValue(CHILDREN_CONDITIONS);

	  for (ConditionConstantIterator cn = ChildrenConditions.begin() ; cn != ChildrenConditions.end(); ++cn)
	    {
	      cn->Set(TO_ERASE);
	    }


	  counter++;

	}

      //update the list of old conditions with the list of new conditions
      list_of_conditions = list_of_new_conditions;

      KRATOS_CATCH( "" )
    }



    //*******************************************************************************************
    //*******************************************************************************************

    virtual void SetNewConditionVariables(ConditionType::Pointer& pOldCondition, ConditionType::Pointer& pNewCondition)
    {
      KRATOS_TRY

      //set variables
      pNewCondition->SetValue( MASTER_NODES          , pOldCondition->GetValue(MASTER_NODES)         );
      pNewCondition->SetValue( NORMAL                , pOldCondition->GetValue(NORMAL)               );
      pNewCondition->SetValue( CAUCHY_STRESS_VECTOR  , pOldCondition->GetValue(CAUCHY_STRESS_VECTOR) );
      pNewCondition->SetValue( DEFORMATION_GRADIENT  , pOldCondition->GetValue(DEFORMATION_GRADIENT) );

      KRATOS_CATCH( "" )
    }


    //*******************************************************************************************
    //*******************************************************************************************

    void SetNodesToModelPart(ModelPart& rModelPart, std::vector<NodeType::Pointer>& list_of_nodes)
    {
      KRATOS_TRY

      if(list_of_nodes.size()){

	//add new conditions: ( SOLID body model part )
	for(std::vector<NodeType::Pointer>::iterator i_node = list_of_nodes.begin(); i_node!= list_of_nodes.end(); ++i_node)
	  {
	    rModelPart.Nodes().push_back(*(i_node));
	  }

      }

      KRATOS_CATCH( "" )
    }
    //*******************************************************************************************
    //*******************************************************************************************

    void SetConditionsToModelPart(ModelPart& rModelPart, std::vector<ConditionType::Pointer>& list_of_conditions)
    {
      KRATOS_TRY

      if(list_of_conditions.size()){

	//clear erased conditions: ( SOLID body model part )
	this->CleanModelPartConditions(rModelPart);

	//add new conditions: ( SOLID body model part )
	for(std::vector<ConditionType::Pointer>::iterator i_cond = list_of_conditions.begin(); i_cond!= list_of_conditions.end(); ++i_cond)
	  {
	    rModelPart.Conditions().push_back(*(i_cond));
	  }

      }

       if( this->mEchoLevel > 0 )
	std::cout<<"   [ CONDITIONS ( inserted : "<<list_of_conditions.size()<<" ) ]"<<std::endl;

       mrRemesh.Info->InsertedBoundaryConditions = list_of_conditions.size();

      //renumerate conditions
      // unsigned int id=1;
      // for(ModelPart::ConditionsContainerType::iterator i_cond = rModelPart.ConditionsBegin(); i_cond!= rModelPart.ConditionsEnd(); ++i_cond)
      // 	{
      // 	  i_cond->SetId(id);
      // 	  id++;
      // 	}


      KRATOS_CATCH( "" )
    }


    //*******************************************************************************************
    //*******************************************************************************************

    void CleanModelPartConditions(ModelPart& rModelPart)
    {
      KRATOS_TRY


      //clean old conditions (TO_ERASE) and add new conditions (NEW_ENTITY)
      ModelPart::ConditionsContainerType PreservedConditions;
      PreservedConditions.reserve(rModelPart.Conditions().size());
      PreservedConditions.swap(rModelPart.Conditions());

      for(ModelPart::ConditionsContainerType::iterator i_cond = PreservedConditions.begin(); i_cond!= PreservedConditions.end(); ++i_cond)
    	{
    	  if(i_cond->IsNot(TO_ERASE))
    	    rModelPart.Conditions().push_back(*(i_cond.base()));
    	}

      if( this->mEchoLevel > 0 )
	std::cout<<"   [ CONDITIONS ( erased : "<<PreservedConditions.size()-rModelPart.Conditions().size()<<" ) ]"<<std::endl;


      KRATOS_CATCH( "" )
    }

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

    ///@name Private Static Member Variables
    ///@{

    ///@}
    ///@name Private Member Variables
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
    RefineConditionsMesherProcess& operator=(RefineConditionsMesherProcess const& rOther);


    /// this function is a private function


    /// Copy constructor.
    //Process(Process const& rOther);


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
                                  RefineConditionsMesherProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const RefineConditionsMesherProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_REFINE_CONDITIONS_MESHER_PROCESS_H_INCLUDED  defined
