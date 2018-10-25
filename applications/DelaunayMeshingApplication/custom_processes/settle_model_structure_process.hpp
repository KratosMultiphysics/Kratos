//
//   Project Name:        KratosDelaunayMeshingApplication $
//   Created by:          $Author:             JMCarbonell $
//   Last modified by:    $Co-Author:                      $
//   Date:                $Date:                April 2018 $
//   Revision:            $Revision:                   0.0 $
//
//

#if !defined(KRATOS_SETTLE_MODEL_STRUCTURE_PROCESS_H_INCLUDED )
#define  KRATOS_SETTLE_MODEL_STRUCTURE_PROCESS_H_INCLUDED


// System includes

// External includes

// Project includes
#include "custom_utilities/mesher_utilities.hpp"
#include "custom_processes/nodal_neighbours_search_process.hpp"
#include "custom_processes/build_model_part_boundary_process.hpp"

#include "delaunay_meshing_application_variables.h"

///VARIABLES used:
//Data:
//StepData:
//Flags:    (checked)
//          (set)
//          (modified)
//          (reset)


namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{
typedef  ModelPart::ConditionType                                ConditionType;
typedef  ModelPart::NodesContainerType                      NodesContainerType;
typedef  ModelPart::ElementsContainerType                ElementsContainerType;
typedef  ModelPart::ConditionsContainerType            ConditionsContainerType;
typedef  ModelPart::MeshType::GeometryType::PointsArrayType    PointsArrayType;

typedef PointerVectorSet<ConditionType, IndexedObject> ConditionsContainerType;
typedef ConditionsContainerType::iterator                    ConditionIterator;
typedef ConditionsContainerType::const_iterator      ConditionConstantIterator;

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
class SettleModelStructureProcess
    : public Process
{
 public:
  ///@name Type Definitions
  ///@{

  /// Pointer definition of SettleModelStructureProcess
  KRATOS_CLASS_POINTER_DEFINITION( SettleModelStructureProcess );

  ///@}
  ///@name Life Cycle
  ///@{

  /// Default constructor.
  SettleModelStructureProcess(ModelPart& rMainModelPart, Flags Options, int EchoLevel = 0)
      : mrMainModelPart(rMainModelPart)
  {
    mOptions   = Options;
    mEchoLevel = EchoLevel;
  }

  /// Destructor.
  virtual ~SettleModelStructureProcess()
  {
  }

  ///@}
  ///@name Operators
  ///@{

  void operator()()
  {
    Execute();
  }

  ///@}
  ///@name Operations
  ///@{

  /// this function is designed for being called at the beginning of the computations
  /// right after reading the model and the groups
  void ExecuteInitialize() override
  {
    KRATOS_TRY

    //Clean Nodal and Conditional Flags
    this->CleanMeshFlags(mrMainModelPart);

    //Sort Conditions
    this->SortModelPartConditions();

    KRATOS_CATCH(" ")
  }


  /// this function is designed for being called at the end of the computations
  /// right after reading the model and the groups
  void ExecuteFinalize() override
  {
    KRATOS_TRY

    // Restore meshes coherency
    this->BuildModelPartStructure();

    // Perform searches for next processes (i.e. contact search)
    this->PerformModelSearches();

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
    return "SettleModelStructureProcess";
  }

  /// Print information about this object.
  void PrintInfo(std::ostream& rOStream) const override
  {
    rOStream << "SettleModelStructureProcess";
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

  ModelPart& mrMainModelPart;

  Flags mOptions;

  int mEchoLevel;

  ///@}
  ///@name Protected Operators
  ///@{

  ///@}
  ///@name Protected Operations
  ///@{


  //*******************************************************************************************
  //*******************************************************************************************

  void SortModelPartConditions()
  {

    KRATOS_TRY

    //Sort Conditions
    unsigned int consecutive_index = 1;
    for(ModelPart::ConditionsContainerType::iterator i_cond = mrMainModelPart.ConditionsBegin(); i_cond!=mrMainModelPart.ConditionsEnd(); ++i_cond)
      (i_cond)->SetId(++consecutive_index);

    mrMainModelPart.Conditions().Sort();
    // mrMainModelPart.Conditions().Unique();

    KRATOS_CATCH(" ")
  }


  //*******************************************************************************************
  //*******************************************************************************************

  void SortModelPartElements()
  {

    KRATOS_TRY

    // Sort Elements
    unsigned int consecutive_index = 1;
    for(ModelPart::ElementsContainerType::iterator i_elem = mrMainModelPart.ElementsBegin(); i_elem!=mrMainModelPart.ElementsEnd(); ++i_elem)
      i_elem->SetId(++consecutive_index);

    mrMainModelPart.Elements().Sort();
    // mrMainModelPart.Elements().Unique();

    KRATOS_CATCH(" ")

  }


  //*******************************************************************************************
  //*******************************************************************************************

  void SortModelPartNodes()
  {

    KRATOS_TRY

    //Sort Nodes, set STRUCTURE nodes at end
    unsigned int consecutive_index = 1;
    unsigned int reverse_index = mrMainModelPart.Nodes().size();
    for(ModelPart::NodesContainerType::iterator i_node = mrMainModelPart.NodesBegin(); i_node!=mrMainModelPart.NodesEnd(); ++i_node)
    {
      if(i_node->IsNot(STRUCTURE) ){
        i_node->SetId(++consecutive_index);
      }
      else{
        i_node->SetId(reverse_index--);
      }
    }

    mrMainModelPart.Nodes().Sort();
    // mrMainModelPart.Nodes().Unique();

    KRATOS_CATCH(" ")

  }


  //*******************************************************************************************
  //*******************************************************************************************

  void BuildModelPartStructure()
  {

    KRATOS_TRY

    //Once all model parts are build, the main model part must be reconstructed coherently
    unsigned int NumberOfSubModelParts=mrMainModelPart.NumberOfSubModelParts();
    if(NumberOfSubModelParts>0){
      this->BuildTotalModelPart(mrMainModelPart, mEchoLevel);
    }
    else{
      this->CleanMeshFlags(mrMainModelPart);
    }


    KRATOS_CATCH(" ")

  }

  //*******************************************************************************************
  //*******************************************************************************************

  virtual void PerformModelSearches()
  {
    KRATOS_TRY

    //NODAL NEIGHBOURS SEARCH
    NodalNeighboursSearchProcess FindNeighbours(mrMainModelPart);
    FindNeighbours.Execute();

    //NODAL_H SEARCH
    //FindNodalHProcess FindNodalH(mrMainModelPart);
    //FindNodalH.Execute();

    //CONDITIONS MASTER_ELEMENTS and MASTER_NODES SEARCH
    BuildModelPartBoundaryProcess BuildBoundaryProcess(mrMainModelPart, mrMainModelPart.Name(), mEchoLevel);
    BuildBoundaryProcess.SearchConditionMasters();

    //BOUNDARY NORMALS SEARCH and SHRINKAGE FACTOR
    BoundaryNormalsCalculationUtilities BoundaryComputation;
    BoundaryComputation.CalculateWeightedBoundaryNormals(mrMainModelPart, mEchoLevel);

    KRATOS_CATCH(" ")
  }


  //*******************************************************************************************
  //*******************************************************************************************

  virtual void BuildTotalModelPart(ModelPart& rModelPart, int EchoLevel)
  {

    KRATOS_TRY

    //Mesh Id=0

    if( EchoLevel > 0 )
      std::cout<<"   [ START MODEL PART ["<<rModelPart.Name()<<"] [Elems=:"<<rModelPart.NumberOfElements()<<"|Nodes="<<rModelPart.NumberOfNodes()<<"|Conds="<<rModelPart.NumberOfConditions()<<"] ] "<<std::endl;

    rModelPart.Nodes().clear();
    rModelPart.Elements().clear();

    //contact conditions are located on Mesh_0
    ModelPart::ConditionsContainerType PreservedConditions;

    unsigned int nodeId=1;
    unsigned int elemId=1;
    unsigned int condId=1;

    this->BuildBodyModelParts(rModelPart, PreservedConditions, nodeId, elemId, condId);

    this->BuildBoundaryModelParts(rModelPart,PreservedConditions, nodeId, elemId, condId);

    this->BuildContactModelParts(rModelPart, PreservedConditions, nodeId, elemId, condId);

    //now set new conditions
    rModelPart.Conditions().swap(PreservedConditions);

    //Sort
    rModelPart.Nodes().Sort();
    rModelPart.Elements().Sort();
    rModelPart.Conditions().Sort();

    //Unique
    rModelPart.Nodes().Unique();
    rModelPart.Elements().Unique();
    rModelPart.Conditions().Unique();

    //Sort Again to have coherent numeration for nodes (mesh with shared nodes)
    unsigned int consecutive_index = 1;
    for(ModelPart::NodesContainerType::iterator in = rModelPart.NodesBegin() ; in != rModelPart.NodesEnd() ; ++in)
      in->SetId(++consecutive_index);

    this->BuildComputingDomain(rModelPart, EchoLevel);

    if( EchoLevel > 0 )
      std::cout<<"   [ END MODEL PART ["<<rModelPart.Name()<<"] [Elems=:"<<rModelPart.NumberOfElements()<<"|Nodes="<<rModelPart.NumberOfNodes()<<"|Conds="<<rModelPart.NumberOfConditions()<<"] ] "<<std::endl;


    KRATOS_CATCH(" ")

  }


  //*******************************************************************************************
  //*******************************************************************************************

  void BuildBodyModelParts(ModelPart& rModelPart, ModelPart::ConditionsContainerType& rPreservedConditions, unsigned int& rNodeId, unsigned int& rElemId, unsigned int& rCondId)
  {
    KRATOS_TRY

    //Add Fluid Bodies modelparts to main modelpart  flags: ( FLUID / NOT_ACTIVE / NOT_BOUNDARY )
    //Add Solid Bodies modelparts to main modelpart  flags: ( SOLID / NOT_ACTIVE / NOT_BOUNDARY )
    //Add Rigid Bodies modelparts to main modelpart  flags: ( RIGID / NOT_ACTIVE / NOT_BOUNDARY )

    const array_1d<double,3> ZeroNormal(3,0.0);

    for(ModelPart::SubModelPartIterator i_mp= rModelPart.SubModelPartsBegin() ; i_mp!=rModelPart.SubModelPartsEnd(); ++i_mp)
    {

      bool add_to_main_model_part = false;

      if( i_mp->IsNot(ACTIVE) && i_mp->IsNot(BOUNDARY) ){ //only the domains (no computing, no boundary)

        if( i_mp->Is(SOLID) || i_mp->Is(FLUID) || i_mp->Is(RIGID) || i_mp->IsNot(CONTACT) )
          add_to_main_model_part = true;
      }


      if( add_to_main_model_part  ){


        if( mEchoLevel > 0 )
          std::cout<<"    [ SUBMODEL PART ["<<i_mp->Name()<<"] [Elems="<<i_mp->NumberOfElements()<<"|Nodes="<<i_mp->NumberOfNodes()<<"|Conds="<<i_mp->NumberOfConditions()<<"] ";//Clean Nodes when redefining the main model part:
        ModelPart::NodesContainerType temporal_nodes;
        temporal_nodes.reserve(i_mp->Nodes().size());
        temporal_nodes.swap(i_mp->Nodes());


        if( !i_mp->NumberOfElements() ){ // usually rigid domains (however rigid domain have rigid surface elements)

          for(ModelPart::NodesContainerType::iterator i_node = temporal_nodes.begin() ; i_node != temporal_nodes.end() ; ++i_node)
          {
            if( i_node->IsNot(TO_ERASE) ){
              (i_mp->Nodes()).push_back(*(i_node.base()));
              (rModelPart.Nodes()).push_back(*(i_node.base()));
              rModelPart.Nodes().back().SetId(rNodeId);
              rNodeId+=1;
            }
          }
        }
        else{

          //reset domain flags in nodes before new assignment
          if( i_mp->Is(FLUID) )
            for(ModelPart::NodesContainerType::iterator i_node = temporal_nodes.begin() ; i_node != temporal_nodes.end() ; ++i_node)
              i_node->Set(FLUID,false);

          if( i_mp->Is(SOLID) )
            for(ModelPart::NodesContainerType::iterator i_node = temporal_nodes.begin() ; i_node != temporal_nodes.end() ; ++i_node)
              i_node->Set(SOLID,false);

          //set new element and domain flags to nodes
          for(ModelPart::ElementsContainerType::iterator i_elem = i_mp->ElementsBegin() ; i_elem != i_mp->ElementsEnd() ; ++i_elem)
          {
            if( i_elem->IsNot(TO_ERASE) ){ //at this point any element must be TO_ERASE

              PointsArrayType& vertices=i_elem->GetGeometry().Points();

              if(i_mp->Is(SOLID)){

                for(unsigned int i=0; i<vertices.size(); ++i)
                {
                  //vertices[i].Set(BLOCKED,true);
                  vertices[i].Set(SOLID,true);
                }
              }
              else if(i_mp->Is(FLUID)){

                for(unsigned int i=0; i<vertices.size(); ++i)
                {
                  //vertices[i].Set(BLOCKED,true);
                  vertices[i].Set(FLUID,true);
                }

              }

              (rModelPart.Elements()).push_back(*(i_elem.base()));
              rModelPart.Elements().back().SetId(rElemId);
              rElemId+=1;

            }
            else{
              std::cout<<" ELEMENT TO_ERASE must be RELEASED in a previous stage "<<std::endl;
            }
          }

          for(ModelPart::NodesContainerType::iterator i_node = temporal_nodes.begin() ; i_node != temporal_nodes.end() ; ++i_node)
          {

            if( i_node->IsNot(TO_ERASE) ){
              (i_mp->Nodes()).push_back(*(i_node.base()));
              (rModelPart.Nodes()).push_back(*(i_node.base()));
              rModelPart.Nodes().back().SetId(rNodeId);
              ++rNodeId;
            }

            i_node->Set(TO_REFINE,false);

            if ( i_node->IsNot(BOUNDARY)) {
                  i_node->Set(NEW_ENTITY,false);
            }

            // if( (i_node->Is(BLOCKED) || i_node->Is(ISOLATED) ) && i_node->IsNot(TO_ERASE) ){
            //   (i_mp->Nodes()).push_back(*(i_node.base()));
            //   (rModelPart.Nodes()).push_back(*(i_node.base()));
            //   rModelPart.Nodes().back().SetId(rNodeId);
            //   ++rNodeId;
            // }

            // if( i_node->Is(BLOCKED) || i_node->Is(RIGID) ){ //all nodes belonging to an element had been blocked previously

            //   i_node->Set(ISOLATED,false);   //reset isolated

            //   if( i_node->IsNot(TO_ERASE) ){

            //     (i_mp->Nodes()).push_back(*(i_node.base()));
            //     (rModelPart.Nodes()).push_back(*(i_node.base()));
            //     rModelPart.Nodes().back().SetId(rNodeId);
            //     rNodeId+=1;

            //   }

            // }
            // else{

            //   if( i_node->IsNot(SOLID) )
            //     i_node->Set(ISOLATED,true);

            //   if( mOptions.Is(MesherUtilities::KEEP_ISOLATED_NODES) && i_node->IsNot(TO_ERASE) ){

            //     (i_mp->Nodes()).push_back(*(i_node.base()));
            //     (rModelPart.Nodes()).push_back(*(i_node.base()));
            //     rModelPart.Nodes().back().SetId(rNodeId);
            //     rNodeId+=1;

            //   }

            // }

            // i_node->Set(TO_REFINE,false);  //reset if was labeled to refine (to not duplicate boundary conditions)
            // i_node->Set(BLOCKED,false);
            // i_node->Set(NEW_ENTITY,false);

            if(i_node->IsNot(BOUNDARY)){

              if( i_node->SolutionStepsDataHas(CONTACT_FORCE) )
                noalias(i_node->GetSolutionStepValue(CONTACT_FORCE)) = ZeroNormal;
            }

          }

        }

        for(ModelPart::ConditionsContainerType::iterator i_cond = i_mp->ConditionsBegin() ; i_cond != i_mp->ConditionsEnd() ; ++i_cond)
        {

          if( i_cond->IsNot(TO_ERASE) ){
            i_cond->Set(TO_REFINE,false);  //reset if was labeled to refine (to not duplicate boundary conditions)
            rPreservedConditions.push_back(*(i_cond.base()));
            rPreservedConditions.back().SetId(rCondId);
            rCondId+=1;
          }
        }

        if( mEchoLevel > 0 )
          std::cout<<" / [Elems="<<i_mp->NumberOfElements()<<"|Nodes="<<i_mp->NumberOfNodes()<<"|Conds="<<i_mp->NumberOfConditions()<<"] ] "<<std::endl;

      }

    }

    KRATOS_CATCH(" ")
  }


  //*******************************************************************************************
  //*******************************************************************************************

  void BuildBoundaryModelParts(ModelPart& rModelPart, ModelPart::ConditionsContainerType& rPreservedConditions, unsigned int& rNodeId, unsigned int& rElemId, unsigned int& rCondId)
  {

    KRATOS_TRY

    unsigned int body_model_part_conditions = rPreservedConditions.size();

    if(body_model_part_conditions > 0){

      //add new conditions: ( BOUNDARY model parts )
      for(ModelPart::SubModelPartIterator i_mp= mrMainModelPart.SubModelPartsBegin() ; i_mp!=mrMainModelPart.SubModelPartsEnd(); ++i_mp)
      {
        if( i_mp->Is(BOUNDARY) ){ //boundary model part

          if( mEchoLevel > 0 )
            std::cout<<"    [ SUBMODEL PART ["<<i_mp->Name()<<"]  initial [Elems="<<i_mp->NumberOfElements()<<"|Nodes="<<i_mp->NumberOfNodes()<<"|Conds="<<i_mp->NumberOfConditions()<<"] ]"<<std::endl;

          this->CleanModelPartConditions(*i_mp);

          for(ModelPart::ConditionsContainerType::iterator i_cond = rPreservedConditions.begin(); i_cond!= rPreservedConditions.end(); ++i_cond)
          {
            ConditionsContainerType& ChildrenConditions = i_cond->GetValue(CHILDREN_CONDITIONS);

            //this conditions are cloned, then the id has no coherence, must be renumbered at the end of the assignation
            for (ConditionConstantIterator cn = ChildrenConditions.begin() ; cn != ChildrenConditions.end(); ++cn)
            {
              if( cn->GetValue(MODEL_PART_NAME) == i_mp->Name() ){
                i_mp->Conditions().push_back(*(cn.base()));

                if( i_cond->Is(NEW_ENTITY) ){
                  for(unsigned int i=0; i<i_cond->GetGeometry().size(); ++i)
                  {
                    if( i_cond->GetGeometry()[i].Is(NEW_ENTITY) ){
                      (i_mp->Nodes()).push_back(i_cond->GetGeometry()(i));
                      //i_cond->GetGeometry()[i].Set(NEW_ENTITY,false); //reset if was new
                    }
                  }
                  //i_cond->Set(NEW_ENTITY,false); //reset here if the condition is inserted
                }

              }

            }

          }
        }

      }

    }

    // Set new nodes to the dirichlet sub model parts (works in 2D. not shure in 3D).  Must be reviewed
    for(ModelPart::SubModelPartIterator i_model_part = mrMainModelPart.SubModelPartsBegin() ; i_model_part != mrMainModelPart.SubModelPartsEnd(); ++i_model_part)
    {
      if( i_model_part->IsNot(BOUNDARY) &&  i_model_part->IsNot(ACTIVE) && i_model_part->IsNot(RIGID) ){

        for(ModelPart::NodesContainerType::iterator i_node = i_model_part->NodesBegin() ; i_node != i_model_part->NodesEnd(); ++i_node)
        {

          if( i_node->Is(BOUNDARY) && i_node->Is(NEW_ENTITY) ){

            // Generate a list of neighbour nodes
            unsigned int NodeId = i_node->Id();

            std::vector<int> list_of_neighbour_nodes;

            for( ModelPart::ConditionsContainerType::iterator j_cond = rPreservedConditions.begin(); j_cond != rPreservedConditions.end(); ++j_cond)
            {

              bool node_belongs_to_condition = false;
              Geometry< Node<3 > > & rjGeom = j_cond->GetGeometry();

              if ( j_cond->Is(NEW_ENTITY) ){
                for ( unsigned int j = 0; j < rjGeom.size() ; ++j) {
                  if ( rjGeom[j].Id() == NodeId) {
                    node_belongs_to_condition = true;
                    break;
                  }
                }

                if (node_belongs_to_condition){
                  for (unsigned int j = 0; j < rjGeom.size() ; ++j)  {
                    list_of_neighbour_nodes.push_back( rjGeom[j].Id() );
                  }
                }

              }
            }

            if(list_of_neighbour_nodes.size() == 0){
              std::cout << " something wierd, this new node does not have any new neighbour: " << NodeId << std::endl;
              // aqui falta un continue o algu ( no un break)
              continue;
            }


            // unique and sort
            std::sort(list_of_neighbour_nodes.begin(), list_of_neighbour_nodes.end() );
            std::vector<int>::iterator new_end = std::unique( list_of_neighbour_nodes.begin(), list_of_neighbour_nodes.end() );
            list_of_neighbour_nodes.resize( std::distance( list_of_neighbour_nodes.begin(), new_end) );

            for (ModelPart::SubModelPartIterator i_mp = mrMainModelPart.SubModelPartsBegin(); i_mp != mrMainModelPart.SubModelPartsEnd(); ++i_mp)
            {
              if ( i_mp->Is(BOUNDARY) && i_mp->IsNot(CONTACT) && (i_mp->NumberOfConditions() == 0) )
              {
                unsigned int counter = 0;

                for (unsigned int ii = 0; ii < list_of_neighbour_nodes.size(); ++ii)
                {
                  unsigned int target = list_of_neighbour_nodes[ii];

                  for (ModelPart::NodesContainerType::iterator i_node = i_mp->NodesBegin(); i_node != i_mp->NodesEnd(); ++i_node)
                  {
                    if ( i_node->Id() == target) {
                      counter++;
                    }
                  }
                }

                if ( counter == list_of_neighbour_nodes.size()-1)
                  i_mp->Nodes().push_back( *(i_node.base() ) );
              }
            }

          }
        }


      }

    }

    //reset NEW_ENTITIES in conditions
    for( ModelPart::ConditionsContainerType::iterator j_cond = rPreservedConditions.begin(); j_cond != rPreservedConditions.end(); ++j_cond)
    {
      j_cond->Set(NEW_ENTITY,false);
      for(unsigned int j=0; j<j_cond->GetGeometry().size(); ++j)
      {
        if( j_cond->GetGeometry()[j].Is(NEW_ENTITY) ){
          j_cond->GetGeometry()[j].Set(NEW_ENTITY,false); //reset if was new
        }
      }
    }


    //add new nodes: ( BOUNDARY model parts ) and remove erased nodes
    for(ModelPart::SubModelPartIterator i_mp= mrMainModelPart.SubModelPartsBegin() ; i_mp!=mrMainModelPart.SubModelPartsEnd(); ++i_mp)
    {
      if( i_mp->Is(BOUNDARY) ){ //boundary model part

        ModelPart::NodesContainerType temporal_nodes;
        temporal_nodes.reserve(i_mp->Nodes().size());
        temporal_nodes.swap(i_mp->Nodes());

        for(ModelPart::NodesContainerType::iterator i_node = temporal_nodes.begin() ; i_node != temporal_nodes.end() ; ++i_node)
        {
          if( i_node->IsNot(TO_ERASE) )
            (i_mp->Nodes()).push_back(*(i_node.base()));
        }

        if( mEchoLevel > 0 )
          std::cout<<"    [ SUBMODEL PART ["<<i_mp->Name()<<"]  final [Elems="<<i_mp->NumberOfElements()<<"|Nodes="<<i_mp->NumberOfNodes()<<"|Conds="<<i_mp->NumberOfConditions()<<"] ] "<<std::endl;
      }
    }


    //add boundary domain conditions to preserved conditions
    for(ModelPart::SubModelPartIterator i_mp= rModelPart.SubModelPartsBegin() ; i_mp!=rModelPart.SubModelPartsEnd(); ++i_mp)
    {
      if( i_mp->Is(BOUNDARY) ){ //boundary model part

        for(ModelPart::ConditionsContainerType::iterator i_cond = i_mp->ConditionsBegin() ; i_cond != i_mp->ConditionsEnd() ; ++i_cond)
        {
          if( i_cond->IsNot(TO_ERASE) ){
            i_cond->Set(NEW_ENTITY,false); //reset here if the condition is inserted
            rPreservedConditions.push_back(*(i_cond.base()));
            rPreservedConditions.back().SetId(rCondId);
            rCondId+=1;
          }
        }
      }

    }

    KRATOS_CATCH( "" )
  }



  //*******************************************************************************************
  //*******************************************************************************************

  void BuildContactModelParts(ModelPart& rModelPart, ModelPart::ConditionsContainerType& rPreservedConditions, unsigned int& rNodeId, unsigned int& rElemId, unsigned int& rCondId)
  {
    KRATOS_TRY

    //Add Contact modelparts to main modelpart  flags: ( CONTACT ) in contact model parts keep only nodes and contact conditions // after that a contact search will be needed

    //if contact condition has the same geometry size as an elements printing ids will coincide,
    //renumber conditions with rElemId instead of rCondId :: in order to ensure it check maximun and apply it
    unsigned int rContactId = rCondId;
    if( rElemId > rCondId )
      rContactId = rElemId;


    for(ModelPart::SubModelPartIterator i_mp= rModelPart.SubModelPartsBegin() ; i_mp!=rModelPart.SubModelPartsEnd(); ++i_mp)
    {

      if( i_mp->Is(CONTACT) ){ //keep only contact conditions

        if( mEchoLevel > 0 )
          std::cout<<"    [ SUBMODEL PART ["<<i_mp->Name()<<"] [Elems="<<i_mp->NumberOfElements()<<"|Nodes="<<i_mp->NumberOfNodes()<<"|Conds="<<i_mp->NumberOfConditions()<<"] ";

        i_mp->Elements().clear();

        //Clean Nodes when redefining the main model part:
        ModelPart::NodesContainerType temporal_nodes;
        temporal_nodes.reserve(i_mp->Nodes().size());
        temporal_nodes.swap(i_mp->Nodes());


        for(ModelPart::NodesContainerType::iterator i_node = temporal_nodes.begin() ; i_node != temporal_nodes.end() ; ++i_node)
        {
          if( i_node->IsNot(TO_ERASE) ){
            (i_mp->Nodes()).push_back(*(i_node.base()));
          }
        }

        ModelPart::ConditionsContainerType temporal_conditions;
        temporal_conditions.reserve(i_mp->Conditions().size());
        temporal_conditions.swap(i_mp->Conditions());

        for(ModelPart::ConditionsContainerType::iterator i_cond = temporal_conditions.begin() ; i_cond != temporal_conditions.end() ; ++i_cond)
        {

          if( i_cond->Is(CONTACT) ){ //keep only contact conditions

            if( i_cond->IsNot(TO_ERASE) ){ //it can not be to erase

              (i_mp->Conditions()).push_back(*(i_cond.base()));
              rPreservedConditions.push_back(*(i_cond.base()));
              rPreservedConditions.back().SetId(rContactId);
              rContactId+=1;

            }

          }

        }

        if( mEchoLevel > 0 )
          std::cout<<" / [Elems="<<i_mp->NumberOfElements()<<"|Nodes="<<i_mp->NumberOfNodes()<<"|Conds="<<i_mp->NumberOfConditions()<<"] ] "<<std::endl;


      }

    }


    if( rElemId > rCondId )
      rElemId = rContactId;
    else
      rCondId = rContactId;


    KRATOS_CATCH(" ")
  }


  //*******************************************************************************************
  //*******************************************************************************************

  void CleanModelPartConditions(ModelPart& rModelPart)
  {

    KRATOS_TRY

    if( rModelPart.Is(BOUNDARY) )
      rModelPart.Conditions().clear();

    //clean old conditions (TO_ERASE) and add new conditions (NEW_ENTITY)
    // ModelPart::ConditionsContainerType PreservedConditions;
    // PreservedConditions.reserve(rModelPart.Conditions().size());
    // PreservedConditions.swap(rModelPart.Conditions());

    // for(ModelPart::ConditionsContainerType::iterator i_cond = PreservedConditions.begin(); i_cond!= PreservedConditions.end(); ++i_cond)
    // 	{
    // 	  if(i_cond->IsNot(TO_ERASE))
    // 	    rModelPart.Conditions().push_back(*(i_cond.base()));
    // 	}



    KRATOS_CATCH( "" )
  }


  //*******************************************************************************************
  //*******************************************************************************************

  virtual void BuildComputingDomain (ModelPart& rModelPart, int EchoLevel)
  {
    KRATOS_TRY

    std::string ComputingModelPartName;
    for(ModelPart::SubModelPartIterator i_mp= rModelPart.SubModelPartsBegin() ; i_mp!=rModelPart.SubModelPartsEnd(); ++i_mp)
    {
      if( i_mp->Is(ACTIVE) ){ //computing_domain
        ComputingModelPartName = i_mp->Name();
      }
    }


    ModelPart& rComputingModelPart = rModelPart.GetSubModelPart(ComputingModelPartName);

    rComputingModelPart.Nodes().clear();
    rComputingModelPart.Elements().clear();
    rComputingModelPart.Conditions().clear();

    //add all needed computing entities (elements, nodes, conditions)

    for(ModelPart::SubModelPartIterator i_mp= rModelPart.SubModelPartsBegin() ; i_mp!=rModelPart.SubModelPartsEnd(); ++i_mp)
    {
      if( (i_mp->IsNot(BOUNDARY) && i_mp->IsNot(ACTIVE)) || (i_mp->Is(RIGID)) ){

        if( i_mp->IsNot(CONTACT)  ){

          for(ModelPart::NodesContainerType::iterator i_node = i_mp->NodesBegin() ; i_node != i_mp->NodesEnd() ; ++i_node)
          {
            (rComputingModelPart.Nodes()).push_back(*(i_node.base()));
          }

          for(ModelPart::ConditionsContainerType::iterator i_cond = i_mp->ConditionsBegin() ; i_cond != i_mp->ConditionsEnd() ; ++i_cond)
          {
            (rComputingModelPart.Conditions()).push_back(*(i_cond.base()));
          }

          for(ModelPart::ElementsContainerType::iterator i_elem = i_mp->ElementsBegin() ; i_elem != i_mp->ElementsEnd() ; ++i_elem)
          {
            (rComputingModelPart.Elements()).push_back(*(i_elem.base()));
          }

        }
      }
    }

    //add all contact conditions

    for(ModelPart::SubModelPartIterator i_mp= rModelPart.SubModelPartsBegin() ; i_mp!=rModelPart.SubModelPartsEnd(); ++i_mp)
    {
      if( i_mp->Is(CONTACT)  ){

        for(ModelPart::ConditionsContainerType::iterator i_cond = i_mp->ConditionsBegin() ; i_cond != i_mp->ConditionsEnd() ; ++i_cond)
        {
          if( i_cond->Is(CONTACT) )
            (rComputingModelPart.Conditions()).push_back(*(i_cond.base()));
        }

      }

    }


    //Sort
    rComputingModelPart.Nodes().Sort();
    // rComputingModelPart.Elements().Sort();
    // rComputingModelPart.Conditions().Sort();

    //Unique
    rComputingModelPart.Nodes().Unique();
    // rComputingModelPart.Elements().Unique();
    // rComputingModelPart.Conditions().Unique();

    if( EchoLevel > 0 )
      std::cout<<"    [ SUBMODEL PART ["<<rComputingModelPart.Name()<<"] [Elems="<<rComputingModelPart.NumberOfElements()<<"|Nodes="<<rComputingModelPart.NumberOfNodes()<<"|Conds="<<rComputingModelPart.NumberOfConditions()<<"] ] "<<std::endl;


    KRATOS_CATCH(" ")
  }


  //*******************************************************************************************
  //*******************************************************************************************

  void CleanMeshFlags(ModelPart& rModelPart)
  {

    KRATOS_TRY

    for(ModelPart::NodesContainerType::const_iterator i_node = rModelPart.NodesBegin(); i_node != rModelPart.NodesEnd(); ++i_node)
    {

      i_node->Set(NEW_ENTITY,false); //reset here if the node is labeled as insert
      i_node->Set(TO_REFINE,false);  //reset here if the node is labeled as refine (to not duplicate

    }

    for(ModelPart::ConditionsContainerType::iterator i_cond = rModelPart.ConditionsBegin() ; i_cond != rModelPart.ConditionsEnd() ; ++i_cond)
    {
      i_cond->Set(NEW_ENTITY,false); //reset here if the node is inserted
    }

    KRATOS_CATCH(" ")
  }


  //*******************************************************************************************
  //*******************************************************************************************


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
  SettleModelStructureProcess& operator=(SettleModelStructureProcess const& rOther);

  /// Copy constructor.
  //SettleModelStructureProcess(SettleModelStructureProcess const& rOther);


  ///@}

}; // Class SettleModelStructureProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  SettleModelStructureProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const SettleModelStructureProcess& rThis)
{
  rThis.PrintInfo(rOStream);
  rOStream << std::endl;
  rThis.PrintData(rOStream);

  return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_SETTLE_MODEL_STRUCTURE_PROCESS_H_INCLUDED  defined
