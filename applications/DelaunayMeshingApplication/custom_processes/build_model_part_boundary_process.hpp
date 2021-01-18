//
//   Project Name:        KratosDelaunayMeshingApplication $
//   Created by:          $Author:             JMCarbonell $
//   Last modified by:    $Co-Author:                      $
//   Date:                $Date:                April 2018 $
//   Revision:            $Revision:                   0.0 $
//
//

#if !defined(KRATOS_BUILD_MODEL_PART_BOUNDARY_PROCESS_H_INCLUDED )
#define  KRATOS_BUILD_MODEL_PART_BOUNDARY_PROCESS_H_INCLUDED


// System includes

// External includes

// Project includes
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"

#include "geometries/line_2d_2.h"
#include "geometries/line_2d_3.h"
#include "geometries/line_3d_2.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/triangle_3d_6.h"

#include "custom_conditions/composite_condition.hpp"
#include "custom_utilities/boundary_normals_calculation_utilities.hpp"
#include "custom_utilities/mesher_utilities.hpp"
#include "custom_processes/mesher_process.hpp"

#include "delaunay_meshing_application_variables.h"

///VARIABLES used:
//Data:     MASTER_ELEMENTS(set), MASTER_NODES(set), NEIGHBOUR_ELEMENTS(checked)
//StepData:
//Flags:    (checked) CONTACT and RIGID
//          (set)     BOUNDARY and FREE_SURFACE (set)
//          (modified)
//          (reset)
//(set):=(set in this process)

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{
typedef ModelPart::NodesContainerType NodesContainerType;
typedef ModelPart::ElementsContainerType ElementsContainerType;
typedef ModelPart::ConditionsContainerType ConditionsContainerType;

typedef GlobalPointersVector<Node<3> > NodeWeakPtrVectorType;
typedef GlobalPointersVector<Element> ElementWeakPtrVectorType;
typedef GlobalPointersVector<Condition> ConditionWeakPtrVectorType;
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
class BuildModelPartBoundaryProcess
    : public MesherProcess
{
 public:
  ///@name Type Definitions
  ///@{

  typedef ModelPart::ConditionType         ConditionType;
  typedef ConditionType::GeometryType       GeometryType;

  /// Pointer definition of BuildModelPartBoundaryProcess
  KRATOS_CLASS_POINTER_DEFINITION( BuildModelPartBoundaryProcess );

  ///@}
  ///@name Life Cycle
  ///@{

  /// Default constructor.
  BuildModelPartBoundaryProcess(ModelPart& rModelPart,
                                std::string const rModelPartName,
                                int EchoLevel = 0)
      : mrModelPart(rModelPart)
  {
    mModelPartName = rModelPartName;
    mEchoLevel     = EchoLevel;
  }

  /// Destructor.
  virtual ~BuildModelPartBoundaryProcess()
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

  void Execute() override
  {

    KRATOS_TRY

    bool success=false;

    double begin_time = OpenMPUtils::GetCurrentTime();

    unsigned int NumberOfSubModelParts=mrModelPart.NumberOfSubModelParts();

    this->ResetFreeSurfaceFlag(mrModelPart);

    if( mModelPartName == mrModelPart.Name() ){

      for(auto& i_mp : mrModelPart.SubModelParts())
      {

        if( mEchoLevel >= 1 )
          std::cout<<" [ Construct Boundary on ModelPart ["<<i_mp.Name()<<"] ]"<<std::endl;

        success=UniqueSkinSearch(i_mp);

        if(!success)
        {
          std::cout<<"  ERROR: BOUNDARY CONSTRUCTION FAILED ModelPart : ["<<i_mp.Name()<<"] "<<std::endl;
        }
        else
        {
          if( mEchoLevel >= 1 ){
            double end_time = OpenMPUtils::GetCurrentTime();
            std::cout<<" [ Performed in Time = "<<end_time-begin_time<<" ]"<<std::endl;
          }
          //PrintSkin(i_mp);
        }
      }
    }
    else{

      if( mEchoLevel >= 1 )
        std::cout<<" [ Construct Boundary on ModelPart["<<mModelPartName<<"] ]"<<std::endl;

      ModelPart& rModelPart = mrModelPart.GetSubModelPart(mModelPartName);
      success=UniqueSkinSearch(rModelPart);

      if(!success)
      {
        std::cout<<"  ERROR: BOUNDARY CONSTRUCTION FAILED on ModelPart : ["<<rModelPart.Name()<<"] "<<std::endl;
      }
      else
      {
        if( mEchoLevel >= 1 ){
          double end_time = OpenMPUtils::GetCurrentTime();
          std::cout<<" [ Performed in Time = "<<end_time-begin_time<<" ]"<<std::endl;
        }
        //PrintSkin(rModelPart);
      }
    }

    if( NumberOfSubModelParts > 1 ){
      SetMainModelPartConditions();
      SetComputingModelPart();
    }

    //ComputeBoundaryNormals BoundUtils;
    BoundaryNormalsCalculationUtilities BoundaryComputation;
    if( mModelPartName == mrModelPart.Name() ){
      BoundaryComputation.CalculateWeightedBoundaryNormals(mrModelPart, mEchoLevel);
    }
    else{
      ModelPart& rModelPart = mrModelPart.GetSubModelPart(mModelPartName);
      BoundaryComputation.CalculateWeightedBoundaryNormals(rModelPart, mEchoLevel);
    }


    if( mEchoLevel >= 1 )
      std::cout<<"  Boundary Normals Computed and Assigned ] "<<std::endl;

    KRATOS_CATCH( "" )

  }

  void CheckMasterElement(int Id, const Element::WeakPointer& old_nelem, const Element::WeakPointer& new_nelem)
  {
    if( mEchoLevel >= 1 ){
      KRATOS_ERROR << "sorry weak pointers are not any longer supported " << std::endl;
      // if(!old_nelem.expired()){
      //   if(old_nelem->Id() != new_nelem->Id())
      //     std::cout<<"Condition "<<Id<<" WARNING: master elements ("<<old_nelem->Id()<<" != "<<new_nelem->Id()<<")"<<std::endl;
      // }
      // else{
      //   std::cout<<"Condition "<<Id<<" WARNING: master elements (expired != "<<new_nelem->Id()<<")"<<std::endl;
      // }
    }

  }


  void CheckMasterNode(int Id, const Node<3>::WeakPointer& old_nnode, const Node<3>::Pointer& new_nnode)
  {

    if( mEchoLevel >= 1 ){
      KRATOS_ERROR << "sorry weak pointers are not any longer supported " << std::endl;

      // if(!old_nnode.expired()){
      //   if(old_nnode->Id() != new_nnode->Id())
      //     std::cout<<"Condition "<<Id<<" WARNING: master nodes ("<<old_nnode->Id()<<" != "<<new_nnode->Id()<<")"<<std::endl;
      // }
      // else{
      //   std::cout<<"Condition "<<Id<<" WARNING: master nodes (expired != "<<new_nnode->Id()<<")"<<std::endl;
      // }
    }

  }


  bool SearchConditionMasters()
  {

    KRATOS_TRY

    int composite_conditions = 0;
    int total_conditions = 0;
    int counter = 0;

    bool found=false;

    if(mEchoLevel > 1)
      std::cout<<" [ START SEARCH CONDITIONS MASTERS : "<<std::endl;

    for(auto& i_cond : mrModelPart.Conditions())
    {
      if(i_cond.Is(BOUNDARY)) //composite condition
        ++composite_conditions;

      if(mEchoLevel > 1){
        KRATOS_ERROR << "sorry weak pointers are not any longer supported " << std::endl;

        // std::cout<<" Masters::Condition ("<<i_cond.Id()<<")";
        // auto& MasterElements = i_cond.GetValue(MASTER_ELEMENTS);
        // if(MasterElements.size()!=0){
        //   if(!MasterElements(0).expired())
        //     std::cout<<" ME="<<MasterElements.front().Id()<<" (size:"<<MasterElements.size()<<")";
        // }
        // auto& MasterNodes = i_cond.GetValue(MASTER_NODES);
        // if(MasterNodes.size()!=0){
        //   if(!MasterNodes(0).expired())
        //     std::cout<<" MN= "<<MasterNodes.front().Id()<<" (size:"<<MasterNodes.size()<<")";
        // }
        // std::cout<<std::endl;
      }

      //********************************************************************

      DenseMatrix<unsigned int> lpofa; //connectivities of points defining faces
      DenseVector<unsigned int> lnofa; //number of points defining faces

      GeometryType& cGeometry = i_cond.GetGeometry();
      unsigned int size = cGeometry.size();

      bool perform_search = true;
      for(unsigned int i=0; i<size; ++i)
      {
        if( cGeometry[i].Is(RIGID) ) //if is a rigid wall do not search else do search
          perform_search = false;
      }


      if( i_cond.Is(CONTACT) )
        perform_search = false;

      //********************************************************************
      found=false;

      if( perform_search )
      {

        if(size == 2){

          auto& nElements1 = cGeometry[0].GetValue(NEIGHBOUR_ELEMENTS);
          auto& nElements2 = cGeometry[1].GetValue(NEIGHBOUR_ELEMENTS);

          if( nElements1.size() == 0 || nElements2.size() == 0 )
            std::cout<<" NO SIZE in NEIGHBOUR_ELEMENTS "<<std::endl;

          KRATOS_ERROR << "sorry weak pointers are not any longer supported " << std::endl;


          for(auto i_nelem(nElements1.begin()); i_nelem != nElements1.end(); ++i_nelem)
          {
            KRATOS_ERROR << "sorry weak pointers are not any longer supported " << std::endl;
            //if(!i_nelem.base()->expired())
            {

              for(auto j_nelem(nElements2.begin()); j_nelem != nElements2.end(); ++j_nelem)
              {
                //if(!j_nelem.base()->expired())
                KRATOS_ERROR << "sorry weak pointers are not any longer supported " << std::endl;
                {

                  if(j_nelem->Id() == i_nelem->Id() && !found){

                    auto& rMasterElements = i_cond.GetValue(MASTER_ELEMENTS);

                    if(rMasterElements.size()){
                      this->CheckMasterElement(i_cond.Id(),rMasterElements(0),*i_nelem.base());
                      rMasterElements.clear();
                    }
                    else{
                      if(mEchoLevel >= 1)
                        std::cout<<" First Assignment of Master Elements "<<i_cond.Id()<<std::endl;
                    }

                    rMasterElements.push_back(*i_nelem.base());

                    GeometryType& eGeometry = i_nelem->GetGeometry();

                    //get matrix nodes in faces
                    eGeometry.NodesInFaces(lpofa);
                    eGeometry.NumberNodesInFaces(lnofa);

                    int node = 0;
                    for (unsigned int iface=0; iface<eGeometry.size(); ++iface)
                    {
                      MesherUtilities MesherUtils;
                      found = MesherUtils.FindCondition(cGeometry,eGeometry,lpofa,lnofa,iface);

                      if(found){
                        node=iface;
                        break;
                      }
                    }

                    if(found){
                      NodeWeakPtrVectorType& rMasterNodes = i_cond.GetValue(MASTER_NODES);

                      if(rMasterNodes.size()){
                        const auto& pMaster = rMasterNodes(0);
                        const auto& pother = eGeometry(lpofa(0,node));
                        this->CheckMasterNode(i_cond.Id(),pMaster,pother);
                        rMasterNodes.clear();
                      }
                      else{
                        if(mEchoLevel >= 1)
                          std::cout<<" First Assignment of Master Nodes "<<i_cond.Id()<<std::endl;
                      }

                      rMasterNodes.push_back(eGeometry(lpofa(0,node)));
                    }
                    else{
                      std::cout<<" 2N Geometry MASTER_NODE not FOUND : something is wrong "<<std::endl;
                    }
                  }
                }
              }
            }
          }
        }
        if(size == 3){

          auto& nElements1 = cGeometry[0].GetValue(NEIGHBOUR_ELEMENTS);
          auto& nElements2 = cGeometry[1].GetValue(NEIGHBOUR_ELEMENTS);
          auto& nElements3 = cGeometry[2].GetValue(NEIGHBOUR_ELEMENTS);

          if( nElements1.size() == 0 || nElements2.size() == 0 || nElements3.size() == 0 )
            KRATOS_WARNING("")<<" NO SIZE in NEIGHBOUR_ELEMENTS "<<std::endl;

          for(auto i_nelem(nElements1.begin()); i_nelem != nElements1.end(); ++i_nelem)
          {
            KRATOS_ERROR << "sorry weak pointers are not any longer supported " << std::endl;
            //if(!i_nelem.base()->expired())
            {
              for(auto j_nelem(nElements2.begin()); j_nelem != nElements2.end(); ++j_nelem)
              {
                KRATOS_ERROR << "sorry weak pointers are not any longer supported " << std::endl;
                //if(!j_nelem.base()->expired())
                {
                  if(j_nelem->Id() == i_nelem->Id() && !found)
                  {
                    for(auto k_nelem(nElements3.begin()); k_nelem != nElements3.end(); ++k_nelem)
                    {
                      KRATOS_ERROR << "sorry weak pointers are not any longer supported " << std::endl;
                      //if(!k_nelem.base()->expired())
                      {
                        if(k_nelem->Id() == i_nelem->Id() && !found)
                        {
                          auto& rMasterElements = i_cond.GetValue(MASTER_ELEMENTS);

                          if(rMasterElements.size()){
                            this->CheckMasterElement(i_cond.Id(),rMasterElements(0),*i_nelem.base());
                            rMasterElements.clear();
                          }
                          else{
                            if(mEchoLevel >= 1)
                              std::cout<<" First Assignment of Master Elements "<<i_cond.Id()<<std::endl;
                          }

                          rMasterElements.push_back(*i_nelem.base());

                          GeometryType& eGeometry = i_nelem->GetGeometry();

                          //get matrix nodes in faces
                          eGeometry.NodesInFaces(lpofa);
                          eGeometry.NumberNodesInFaces(lnofa);

                          int node = 0;
                          for(unsigned int iface=0; iface<eGeometry.size(); ++iface)
                          {
                            MesherUtilities MesherUtils;
                            found = MesherUtils.FindCondition(cGeometry,eGeometry,lpofa,lnofa,iface);

                            if(found){
                              node=iface;
                              break;
                            }
                          }

                          if(found){

                            NodeWeakPtrVectorType& rMasterNodes = i_cond.GetValue(MASTER_NODES);

                            if(rMasterNodes.size()){
                              this->CheckMasterNode(i_cond.Id(),rMasterNodes(0),eGeometry(lpofa(0,node)));
                              rMasterNodes.clear();
                            }
                            else{
                              if(mEchoLevel >= 1)
                                std::cout<<" First Assignment of Master Nodes "<<i_cond.Id()<<std::endl;
                            }

                            rMasterNodes.push_back(eGeometry(lpofa(0,node)));

                          }
                          else{
                            std::cout<<" 3N Geometry MASTER_NODE not FOUND : something is wrong "<<std::endl;
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }

        total_conditions++;
      }

      //********************************************************************

      if(mEchoLevel > 1){
        KRATOS_ERROR << "sorry weak pointers are not any longer supported " << std::endl;
        // std::cout<<" SearchResult::Condition ("<<i_cond.Id()<<")";
        // auto MasterElements = i_cond.GetValue(MASTER_ELEMENTS);
        // if(MasterElements.size()!=0){
        //   if(!MasterElements(0).expired())
        //     std::cout<<" ME="<<MasterElements.front().Id()<<" (size:"<<MasterElements.size()<<")";
        // }
        // auto MasterNodes = i_cond.GetValue(MASTER_NODES);
        // if(MasterNodes.size()!=0){
        //   if(!MasterNodes(0).expired())
        //     std::cout<<" MN= "<<MasterNodes.front().Id()<<" (size:"<<MasterNodes.size()<<")";
        // }
        // std::cout<<std::endl;
      }

      if(found)
        ++counter;

    }

    if(counter == total_conditions){
      if(mEchoLevel >= 1)
        std::cout<<"   Condition Masters (ModelPart "<<mrModelPart.Name()<<"): LOCATED ["<<counter<<"]"<<std::endl;
      found=true;
    }
    else{
      if(mEchoLevel >= 1)
        std::cout<<"   Condition Masters (ModelPart "<<mrModelPart.Name()<<"): not LOCATED ["<<counter-total_conditions<<"]"<<std::endl;
      found=false;
    }

    if(counter != composite_conditions)
      if(mEchoLevel >= 1)
        std::cout<<"   Condition Masters (ModelPart "<<mrModelPart.Name()<<"): LOCATED ["<<counter<<"] COMPOSITE ["<<composite_conditions<<"] NO MATCH"<<std::endl;

    if(mEchoLevel > 1)
      std::cout<<"   END SEARCH CONDITIONS MASTERS ] "<<found<<std::endl;

    return found;

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
  std::string Info() const override
  {
    return "BuildModelPartBoundaryProcess";
  }

  /// Print information about this object.
  void PrintInfo(std::ostream& rOStream) const override
  {
    rOStream << "BuildModelPartBoundaryProcess";
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

  std::string mModelPartName;

  int mEchoLevel;


  ///@}
  ///@name Protected Operators
  ///@{


  ///@}
  ///@name Protected Operations
  ///@{

  bool ClearMasterEntities( ModelPart& rModelPart, ModelPart::ConditionsContainerType& rTemporaryConditions)
  {
    KRATOS_TRY

    for(auto& i_cond : rTemporaryConditions)
    {
      i_cond.GetValue(MASTER_ELEMENTS).clear();
      i_cond.GetValue(MASTER_NODES).clear();

      KRATOS_ERROR_IF(i_cond.GetValue(MASTER_ELEMENTS).size() || i_cond.GetValue(MASTER_NODES).size())<<
          " Master Entities not cleared "<<std::endl;
    }

    return true;

    KRATOS_CATCH( "" )
  }

  bool UniqueSkinSearch( ModelPart& rModelPart )
  {

    KRATOS_TRY

    if( mEchoLevel > 0 ){
      std::cout<<" [ Initial Conditions : "<<rModelPart.Conditions().size()<<std::endl;
    }

    if( !rModelPart.Elements().size() || (rModelPart.Is(ACTIVE)) ){
      if( mEchoLevel > 0 ){
        std::cout<<" [ Final Conditions   : "<<rModelPart.Conditions().size()<<std::endl;
      }
      return true;
    }

    //check if a remesh process has been performed and there is any node to erase
    bool any_node_to_erase = false;
    for(auto& i_node : rModelPart.Nodes())
    {
      if( any_node_to_erase == false )
        if( i_node.Is(TO_ERASE) )
          any_node_to_erase = true;
    }

    //swap conditions for a temporary use
    unsigned int ConditionId=1;
    ModelPart::ConditionsContainerType TemporaryConditions;

    //if there are no conditions check main modelpart mesh conditions
    if( !rModelPart.Conditions().size() ){

      ModelPart::ConditionsContainerType& rConditions = rModelPart.GetParentModelPart().Conditions();

      for(auto i_cond(rConditions.begin()); i_cond != rConditions.end(); ++i_cond)
      {
        TemporaryConditions.push_back(*i_cond.base());
        i_cond->SetId(ConditionId);
        ConditionId++;
      }

    }
    else{

      TemporaryConditions.reserve(rModelPart.Conditions().size());
      TemporaryConditions.swap(rModelPart.Conditions());

      //set consecutive ids in the mesh conditions
      if( any_node_to_erase ){
        for(auto& i_cond : TemporaryConditions)
        {
          GeometryType& cGeometry = i_cond.GetGeometry();
          for(unsigned int i=0; i<cGeometry.size(); ++i)
          {
            if(cGeometry[i].Is(TO_ERASE)){
              i_cond.Set(TO_ERASE);
              break;
            }
          }

          i_cond.SetId(ConditionId);
          ++ConditionId;
        }
      }
      else{
        for(auto& i_cond : TemporaryConditions)
        {
          i_cond.SetId(ConditionId);
          ++ConditionId;
        }
      }

    }

    //control the previous mesh conditions
    std::vector<int> PreservedConditions( TemporaryConditions.size() + 1 );
    std::fill( PreservedConditions.begin(), PreservedConditions.end(), 0 );

    //build new skin for the Modelpart
    this->BuildCompositeConditions(rModelPart, TemporaryConditions, PreservedConditions, ConditionId);

    //add other conditions out of the skin space dimension
    this->AddOtherConditions(rModelPart, TemporaryConditions, PreservedConditions, ConditionId);

    return true;

    KRATOS_CATCH( "" )
  }

  virtual void SetBoundaryFlags( ModelPart& rModelPart )
  {

    KRATOS_TRY

    //clear nodal boundary flag
    for(auto& i_elem : rModelPart.Elements())
    {
      GeometryType& eGeometry = i_elem.GetGeometry();

      for(unsigned int j=0; j<eGeometry.size(); ++j)
      {
        eGeometry[j].Set(BOUNDARY,false);
        if(rModelPart.Is(FLUID)){
          eGeometry[j].Set(FREE_SURFACE,false);
        }
      }
    }

    for(auto& i_elem : rModelPart.Elements())
    {
      GeometryType& eGeometry = i_elem.GetGeometry();

      const unsigned int dimension = eGeometry.WorkingSpaceDimension();

      if( eGeometry.FacesNumber() >= (dimension+1) ){ //3 or 4

        DenseMatrix<unsigned int> lpofa; //connectivities of points defining faces
        DenseVector<unsigned int> lnofa; //number of points defining faces

        //get matrix nodes in faces
        eGeometry.NodesInFaces(lpofa);
        eGeometry.NumberNodesInFaces(lnofa);

        auto& nElements = i_elem.GetValue(NEIGHBOUR_ELEMENTS);

        //loop on neighbour elements of an element
        unsigned int iface=0;
        for(auto& i_nelem : nElements)
        {
          unsigned int NumberNodesInFace = lnofa[iface];
          if(i_nelem.Id() == i_elem.Id())
          {
            //if no neighbour is present => the face is free surface
            unsigned int rigid_nodes = 0;
            unsigned int inlet_nodes = 0;
            unsigned int free_surface_nodes = 0;
            for(unsigned int j=1; j<=NumberNodesInFace; ++j)
            {
              eGeometry[lpofa(j,iface)].Set(BOUNDARY,true);
              if(rModelPart.Is(FLUID)){
                if(eGeometry[lpofa(j,iface)].Is(RIGID) || eGeometry[lpofa(j,iface)].Is(SOLID)){
                  ++rigid_nodes;
                }
                else if(eGeometry[lpofa(j,iface)].Is(INLET)){
                  ++inlet_nodes;
                }
                else{
                  ++free_surface_nodes;
                }
              }
              //std::cout<<" node ["<<j<<"]"<<eGeometry[lpofa(j,iface)].Id()<<std::endl;
            }

            if(rModelPart.Is(FLUID)){
              if( (free_surface_nodes>0 && (rigid_nodes>0 || inlet_nodes>0)) || (rigid_nodes==0 && inlet_nodes==0) ){
                for(unsigned int j=1; j<=NumberNodesInFace; ++j)
                {
                  eGeometry[lpofa(j,iface)].Set(FREE_SURFACE,true);
                }
              }
            }

          } //end face condition

          ++iface;
        } //end loop neighbours

      }
      else{
        //set nodes to BOUNDARY for elements outside of the working space dimension
        for(unsigned int j=0; j<eGeometry.size(); ++j)
        {
          eGeometry[j].Set(BOUNDARY,true);
        }
      }
    }


    KRATOS_CATCH( "" )
  }

  virtual bool BuildCompositeConditions( ModelPart& rModelPart, ModelPart::ConditionsContainerType& rTemporaryConditions, std::vector<int>& rPreservedConditions, unsigned int& rConditionId )
  {

    KRATOS_TRY

    //master conditions must be deleted and set them again in the build
    this->ClearMasterEntities(rModelPart, rTemporaryConditions);

    //properties to be used in the generation
    int number_properties = rModelPart.GetParentModelPart().NumberOfProperties();

    if(number_properties<0)
      KRATOS_ERROR<<" number of properties is "<<number_properties<<std::endl;

    Properties::Pointer properties = rModelPart.GetParentModelPart().GetMesh().pGetProperties(number_properties-1);

    //clear nodal boundary flag
    for(auto& i_elem : rModelPart.Elements())
    {
      GeometryType& eGeometry = i_elem.GetGeometry();

      for(unsigned int j=0; j<eGeometry.size(); ++j)
      {
        eGeometry[j].Set(BOUNDARY,false);
        if(rModelPart.Is(FLUID)){
          eGeometry[j].Set(FREE_SURFACE,false);
        }
      }
    }

    ElementsContainerType& rElements = rModelPart.Elements();
    rConditionId=0;
    for(auto i_elem(rElements.begin()); i_elem != rElements.end(); ++i_elem)
    {
      GeometryType& eGeometry = i_elem->GetGeometry();

      const unsigned int dimension = eGeometry.WorkingSpaceDimension();

      if( eGeometry.FacesNumber() >= (dimension+1) ){ //3 or 4

        //********************************************************************
        /*each face is opposite to the corresponding node number so in 2D triangle
          0 ----- 1 2
          1 ----- 2 0
          2 ----- 0 1
        */

        /*each face is opposite to the corresponding node number so in 3D tetrahedron
          0 ----- 1 2 3
          1 ----- 2 0 3
          2 ----- 0 1 3
          3 ----- 0 2 1
        */
        //********************************************************************

        //finding boundaries and creating the "skin"
        DenseMatrix<unsigned int> lpofa; //connectivities of points defining faces
        DenseVector<unsigned int> lnofa; //number of points defining faces

        auto& nElements = i_elem->GetValue(NEIGHBOUR_ELEMENTS);

        //get matrix nodes in faces
        eGeometry.NodesInFaces(lpofa);
        eGeometry.NumberNodesInFaces(lnofa);

        //loop on neighbour elements of an element
        unsigned int iface=0;
        for(auto& i_nelem : nElements)
        {
          unsigned int NumberNodesInFace = lnofa[iface];
          if(i_nelem.Id() == i_elem->Id())
          {
            //if no neighbour is present => the face is free surface
            unsigned int rigid_nodes = 0;
            unsigned int inlet_nodes = 0;
            unsigned int free_surface_nodes = 0;
            for(unsigned int j=1; j<=NumberNodesInFace; ++j)
            {
              eGeometry[lpofa(j,iface)].Set(BOUNDARY,true);
              if(rModelPart.Is(FLUID)){
                if(eGeometry[lpofa(j,iface)].Is(RIGID) || eGeometry[lpofa(j,iface)].Is(SOLID)){
                  ++rigid_nodes;
                }
                else if(eGeometry[lpofa(j,iface)].Is(INLET)){
                  ++inlet_nodes;
                }
                else{
                  ++free_surface_nodes;
                }
              }
              //std::cout<<" node ["<<j<<"]"<<eGeometry[lpofa(j,iface)].Id()<<std::endl;
            }

            if(rModelPart.Is(FLUID)){
              if( (free_surface_nodes>0 && (rigid_nodes>0 || inlet_nodes>0)) || (rigid_nodes==0 && inlet_nodes==0) ){
                for(unsigned int j=1; j<=NumberNodesInFace; ++j)
                {
                  eGeometry[lpofa(j,iface)].Set(FREE_SURFACE,true);
                }
              }
            }

            //1.- create geometry: points array and geometry type
            Condition::NodesArrayType        FaceNodes;
            Condition::GeometryType::Pointer ConditionVertices;

            FaceNodes.reserve(NumberNodesInFace);

            for(unsigned int j=1; j<=NumberNodesInFace; ++j)
            {
              FaceNodes.push_back(eGeometry(lpofa(j,iface)));
            }


            if( NumberNodesInFace == 2 ){
              if( dimension == 2 )
                ConditionVertices = Kratos::make_shared< Line2D2< Node<3> > >(FaceNodes);
              else
                ConditionVertices = Kratos::make_shared< Line3D2< Node<3> > >(FaceNodes);

            }
            else if ( NumberNodesInFace == 3 ){
              if( dimension == 2 )
                ConditionVertices = Kratos::make_shared< Line2D3< Node<3> > >(FaceNodes);
              else
                ConditionVertices = Kratos::make_shared< Triangle3D3< Node<3> > >(FaceNodes);
            }

            rConditionId +=1;

            //Create a composite condition
            CompositeCondition::Pointer p_cond = Kratos::make_intrusive< CompositeCondition >(rConditionId,ConditionVertices,properties);

            bool condition_found = false;
            bool point_condition = false;

            // Search for existing conditions: start
            for(auto i_cond(rTemporaryConditions.begin()); i_cond != rTemporaryConditions.end(); ++i_cond)
            {
              GeometryType& cGeometry = i_cond->GetGeometry();

              MesherUtilities MesherUtils;
              condition_found = MesherUtils.FindCondition(cGeometry,eGeometry,lpofa,lnofa,iface);

              if( condition_found ){

                p_cond->AddChild(*i_cond.base());

                rPreservedConditions[i_cond->Id()] += 1;

                if( cGeometry.PointsNumber() == 1 )
                  point_condition = true;
              }

            }
            // Search for existing conditions: end

            if( !point_condition ){
              // usually one MasterElement and one MasterNode for 2D and 3D simplex
              // can be more than one in other geometries -> it has to be extended to that cases

              // std::cout<<" ID "<<p_cond->Id()<<" MASTER ELEMENT "<<i_elem->Id()<<std::endl;
              // std::cout<<" MASTER NODE "<<eGeometry[lpofa(0,iface)].Id()<<" or "<<eGeometry[lpofa(NumberNodesInFace,iface)].Id()<<std::endl;

              auto& MasterElements = p_cond->GetValue(MASTER_ELEMENTS);
              MasterElements.push_back(*i_elem.base());

              NodeWeakPtrVectorType& MasterNodes = p_cond->GetValue(MASTER_NODES);
              MasterNodes.push_back(eGeometry(lpofa(0,iface)));

            }

            rModelPart.Conditions().push_back(p_cond);
            // Set new conditions: end

          } //end face condition

          ++iface;
        } //end loop neighbours

      }
      else{
        //set nodes to BOUNDARY for elements outside of the working space dimension
        for(unsigned int j=0; j<eGeometry.size(); ++j)
        {
          eGeometry[j].Set(BOUNDARY,true);
        }
      }
    }


    return true;

    KRATOS_CATCH( "" )
  }

  bool FindConditionID(ModelPart& rModelPart, GeometryType& rGeometry)
  {
    KRATOS_TRY

    //check if the condition exists and belongs to the modelpart checking node Ids
    for(unsigned int i=0; i<rGeometry.size(); ++i)
    {
      for(auto& i_node : rModelPart.Nodes())
      {
        if( rGeometry[i].Id() == i_node.Id() )
          return true;
      }
    }

    return false;

    KRATOS_CATCH( "" )
  }

  void PrintSkin(ModelPart& rModelPart)
  {

    KRATOS_TRY

    //PRINT SKIN:
    std::cout<<" CONDITIONS: geometry nodes ("<<rModelPart.Conditions().size()<<")"<<std::endl;

    for(auto& i_cond : rModelPart.Conditions())
    {
      GeometryType& cGeometry = i_cond.GetGeometry();
      std::cout<<"["<<i_cond.Id()<<"]:"<<std::endl;
      //i_cond.PrintInfo(std::cout);
      std::cout<<"( ";
      for(unsigned int i = 0; i < cGeometry.size(); ++i)
      {
        std::cout<< cGeometry[i].Id()<<", ";
      }
      std::cout<<" ): ";

      i_cond.GetValue(MASTER_ELEMENTS).front().PrintInfo(std::cout);

      std::cout<<std::endl;
    }
    std::cout<<std::endl;

    KRATOS_CATCH( "" )

  }

  bool AddOtherConditions(ModelPart& rModelPart, ModelPart::ConditionsContainerType& rTemporaryConditions, std::vector<int>& rPreservedConditions, unsigned int& rConditionId)
  {

    KRATOS_TRY

    unsigned int counter = 0;

    //add all previous conditions not found in the skin search:
    for(auto& i_cond : rTemporaryConditions)
    {
      if(rPreservedConditions[i_cond.Id()] == 0){ //I have not used the condition and any node of the condition

        if(this->CheckAcceptedCondition(rModelPart, i_cond)){

          GeometryType& cGeometry = i_cond.GetGeometry();

          Condition::NodesArrayType FaceNodes;

          FaceNodes.reserve(cGeometry.size() );

          for(unsigned int j=0; j<cGeometry.size(); ++j)
          {
            FaceNodes.push_back(cGeometry(j));
          }

          rPreservedConditions[i_cond.Id()] += 1;

          rConditionId +=1;

          Condition::Pointer p_cond = i_cond.Clone(rConditionId, FaceNodes);

          this->AddConditionToModelPart(rModelPart,p_cond);

          counter++;
        }
      }
    }

    //control if all previous conditions have been added:
    bool all_assigned = true;
    unsigned int lost_conditions = 0;
    for(unsigned int i=1; i<rPreservedConditions.size(); ++i)
    {
      if( rPreservedConditions[i] == 0 ){
        all_assigned = false;
        lost_conditions++;
      }
    }

    if( mEchoLevel >= 1 ){
      std::cout<<"   Final Conditions   : "<<rModelPart.NumberOfConditions()<<std::endl;
      if(all_assigned == true)
        std::cout<<"   ALL_PREVIOUS_CONDITIONS_RELOCATED "<<std::endl;
      else
        std::cout<<"   SOME_PREVIOUS_CONDITIONS_ARE_LOST [lost_conditions:"<<lost_conditions<<"]"<<std::endl;
    }

    return all_assigned;

    KRATOS_CATCH( "" )

  }

  virtual bool CheckAcceptedCondition(ModelPart& rModelPart, Condition& rCondition)
  {
    KRATOS_TRY

    GeometryType& cGeometry = rCondition.GetGeometry();

    bool accepted = FindConditionID(rModelPart, cGeometry);

    return accepted;

    KRATOS_CATCH( "" )
  }

  virtual void AddConditionToModelPart(ModelPart& rModelPart, Condition::Pointer pCondition)
  {
    KRATOS_TRY

    rModelPart.AddCondition(pCondition);
    //rModelPart.Conditions().push_back(pCondition);

    KRATOS_CATCH( "" )
  }

  void ResetFreeSurfaceFlag(ModelPart& rModelPart)
  {
    //reset the boundary flag in all nodes (set old_entity in edge nodes to be recognized)
    if(rModelPart.Is(FLUID)){

      for(auto& i_node : rModelPart.Nodes())
      {
        i_node.Set(FREE_SURFACE,false);
      }

    }
  }

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

  void SetComputingModelPart()
  {
    KRATOS_TRY

    std::string ComputingModelPartName;
    for(auto& i_mp : mrModelPart.SubModelParts())
    {
      if(i_mp.Is(ACTIVE) && i_mp.IsNot(THERMAL))
        ComputingModelPartName = i_mp.Name();
    }

    ModelPart& rModelPart = mrModelPart.GetSubModelPart(ComputingModelPartName);

    if( mEchoLevel >= 1 ){
      std::cout<<" ["<<ComputingModelPartName<<" :: CONDITIONS [OLD:"<<rModelPart.NumberOfConditions();
    }

    ModelPart::ConditionsContainerType KeepConditions;

    for(auto& i_mp : mrModelPart.SubModelParts())
    {
      if( i_mp.NumberOfElements() && ComputingModelPartName != i_mp.Name() ){ //conditions of model_parts with elements only

        ModelPart::ConditionsContainerType& rConditions = i_mp.Conditions();

        for(auto i_cond(rConditions.begin()); i_cond != rConditions.end(); ++i_cond)
        {
          // i_cond->PrintInfo(std::cout);
          // std::cout<<" -- "<<std::endl;
          KeepConditions.push_back(*i_cond.base());
          // KeepConditions.back().PrintInfo(std::cout);
          // std::cout<<std::endl;
        }
      }
    }

    rModelPart.Conditions().swap(KeepConditions);

    if( mEchoLevel >= 1 ){
      std::cout<<" / NEW:"<<rModelPart.NumberOfConditions()<<"] "<<std::endl;
    }

    KRATOS_CATCH( "" )

  }

  void SetMainModelPartConditions()
  {

    KRATOS_TRY

    if( mEchoLevel >= 1 ){
      std::cout<<" ["<<mrModelPart.Name()<<" :: CONDITIONS [OLD:"<<mrModelPart.NumberOfConditions();
    }

    //contact conditions are located on Mesh_0
    ModelPart::ConditionsContainerType KeepConditions;

    unsigned int condId=1;
    if( mModelPartName == mrModelPart.Name() ){

      for(auto& i_mp : mrModelPart.SubModelParts())
      {
        if( !(i_mp.Is(ACTIVE)) && !(i_mp.Is(CONTACT)) ){
          //std::cout<<" ModelPartName "<<i_mp.Name()<<" conditions "<<i_mp.NumberOfConditions()<<std::endl;
          ModelPart::ConditionsContainerType& rConditions = i_mp.Conditions();

          for(auto i_cond(rConditions.begin()); i_cond != rConditions.end(); ++i_cond)
          {
            // i_cond.PrintInfo(std::cout);
            // std::cout<<" -- "<<std::endl;
            KeepConditions.push_back(*(i_cond.base()));
            KeepConditions.back().SetId(condId);
            condId+=1;
            // KeepConditions.back().PrintInfo(std::cout);
            // std::cout<<std::endl;
          }
        }
      }
    }

    ModelPart::ConditionsContainerType& rConditions = mrModelPart.Conditions();

    for(auto i_cond(rConditions.begin()); i_cond != rConditions.end(); ++i_cond)
    {
      if(i_cond->Is(CONTACT)){
        KeepConditions.push_back(*i_cond.base());
        KeepConditions.back().SetId(condId);
        condId+=1;
        //std::cout<<" -- "<<std::endl;
        //KeepConditions.back().PrintInfo(std::cout);
        //std::cout<<std::endl;
      }

    }


    mrModelPart.Conditions().swap(KeepConditions);

    if( mEchoLevel >= 1 ){
      std::cout<<" / NEW:"<<mrModelPart.NumberOfConditions()<<"] "<<std::endl;
    }

    KRATOS_CATCH( "" )
  }


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
  BuildModelPartBoundaryProcess& operator=(BuildModelPartBoundaryProcess const& rOther);

  /// Copy constructor.
  //BuildModelPartBoundaryProcess(BuildModelPartBoundaryProcess const& rOther);


  ///@}

}; // Class BuildModelPartBoundaryProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  BuildModelPartBoundaryProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const BuildModelPartBoundaryProcess& rThis)
{
  rThis.PrintInfo(rOStream);
  rOStream << std::endl;
  rThis.PrintData(rOStream);

  return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_BUILD_MODEL_PART_BOUNDARY_PROCESS_H_INCLUDED  defined
