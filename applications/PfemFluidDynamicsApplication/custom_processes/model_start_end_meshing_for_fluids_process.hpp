//
//   Project Name:        KratosPfemFluidApplication $
//   Created by:          $Author:      JMCarbonell  $
//   Last modified by:    $Co-Author:                $
//   Date:                $Date:      February 2016  $
//   Revision:            $Revision:            0.0  $
//
//

#if !defined(KRATOS_MODEL_START_END_MESHING_FOR_FLUIDS_PROCESS_H_INCLUDED)
#define KRATOS_MODEL_START_END_MESHING_FOR_FLUIDS_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_processes/settle_model_structure_process.hpp"

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
  typedef ModelPart::ConditionType ConditionType;
  typedef ModelPart::NodesContainerType NodesContainerType;
  typedef ModelPart::ElementsContainerType ElementsContainerType;
  typedef ModelPart::ConditionsContainerType ConditionsContainerType;
  typedef ModelPart::MeshType::GeometryType::PointsArrayType PointsArrayType;

  typedef PointerVectorSet<ConditionType, IndexedObject> ConditionsContainerType;
  typedef ConditionsContainerType::iterator ConditionIterator;
  typedef ConditionsContainerType::const_iterator ConditionConstantIterator;

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
  class ModelStartEndMeshingForFluidsProcess
      : public SettleModelStructureProcess
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ModelStartEndMeshingForFluidsProcess
    KRATOS_CLASS_POINTER_DEFINITION(ModelStartEndMeshingForFluidsProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ModelStartEndMeshingForFluidsProcess(ModelPart &rMainModelPart,
                                         Flags Options,
                                         int EchoLevel = 0)
        : SettleModelStructureProcess(rMainModelPart, Options, EchoLevel)
    {
    }

    /// Destructor.
    virtual ~ModelStartEndMeshingForFluidsProcess()
    {
    }

    void operator()()
    {
      Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    void Execute() override{

    };

    ///@}
    ///@name Operators
    ///@{

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
      return "ModelStartEndMeshingForFluidsProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override
    {
      rOStream << "ModelStartEndMeshingForFluidsProcess";
    }

  protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    //*******************************************************************************************
    //*******************************************************************************************

    void BuildTotalModelPart(ModelPart &rModelPart, int EchoLevel) override
    {

      KRATOS_TRY

      //Mesh Id=0

      if (EchoLevel > 1)
        std::cout << "   [ START MODEL PART [" << rModelPart.Name() << "] [Elems=:" << rModelPart.NumberOfElements() << "|Nodes=" << rModelPart.NumberOfNodes() << "|Conds=" << rModelPart.NumberOfConditions() << "] ] " << std::endl;

      rModelPart.Nodes().clear();
      rModelPart.Elements().clear();

      //contact conditions are located on Mesh_0
      // ModelPart::ConditionsContainerType PreservedConditions;

      unsigned int nodeId = 1;
      unsigned int elemId = 1;
      // unsigned int condId=1;

      for (ModelPart::SubModelPartIterator i_mp = rModelPart.SubModelPartsBegin(); i_mp != rModelPart.SubModelPartsEnd(); i_mp++)
      {
        // if( (i_mp->Is(SOLID) && i_mp->IsNot(ACTIVE)) || (i_mp->Is(BOUNDARY) && i_mp->Is(RIGID)) ){ //only the solid domains (no computing) and the rigid body domains (rigid)
        if ((i_mp->Is(SOLID) && i_mp->IsNot(ACTIVE)) || (i_mp->Is(FLUID) && i_mp->IsNot(ACTIVE)) || (i_mp->Is(BOUNDARY) && i_mp->Is(RIGID)))
        { //only the solid domains (no computing) and the rigid body domains (rigid)

          if (EchoLevel > 1)
            std::cout << "    [ SUBMODEL PART [" << i_mp->Name() << "] [Elems=" << i_mp->NumberOfElements() << "|Nodes=" << i_mp->NumberOfNodes() << "|Conds=" << i_mp->NumberOfConditions() << "] ] " << std::endl;

          //Clean Nodes when redefining the main model part:
          const array_1d<double, 3> ZeroNormal(3, 0.0);
          ModelPart::NodesContainerType temporal_nodes;
          temporal_nodes.reserve(i_mp->Nodes().size());
          temporal_nodes.swap(i_mp->Nodes());

          if (!i_mp->NumberOfElements())
          {
            for (ModelPart::NodesContainerType::iterator i_node = temporal_nodes.begin(); i_node != temporal_nodes.end(); i_node++)
            {
              if (i_node->IsNot(TO_ERASE))
              {
                (i_mp->Nodes()).push_back(*(i_node.base()));
                (rModelPart.Nodes()).push_back(*(i_node.base()));
                rModelPart.Nodes().back().SetId(nodeId);
                nodeId += 1;
              }
            }
          }
          else
          {

            for (ModelPart::ElementsContainerType::iterator i_elem = i_mp->ElementsBegin(); i_elem != i_mp->ElementsEnd(); i_elem++)
            {
              if (i_elem->IsNot(TO_ERASE))
              {

                PointsArrayType &vertices = i_elem->GetGeometry().Points();
                for (unsigned int i = 0; i < vertices.size(); i++)
                {
                  vertices[i].Set(BLOCKED);
                  if (i_mp->Is(FLUID))
                  {
                    vertices[i].Set(FLUID);
                  }
                }
                i_elem->SetId(elemId);
                (rModelPart.Elements()).push_back(*(i_elem.base()));
                rModelPart.Elements().back().SetId(elemId);
                elemId += 1;
              }
            }

            for (ModelPart::NodesContainerType::iterator i_node = temporal_nodes.begin(); i_node != temporal_nodes.end(); i_node++)
            {
              //i_node->PrintInfo(std::cout);
              //std::cout<<std::endl;

              if (i_node->Is(BLOCKED) || i_node->Is(RIGID))
              {
                if (i_node->Is(RIGID) && i_node->IsNot(BLOCKED))
                {
                  // double pressureRigid=i_node->FastGetSolutionStepValue(PRESSURE);
                  i_node->FastGetSolutionStepValue(PRESSURE) = 0;
                  if (i_mp->Is(FLUID))
                  {
                    i_node->Reset(FLUID); //reset isolated
                  }
                  // std::cout<<" fluid 1. SET PRESSURE 0 TO ISOLATED NODE ("<<nodeId<<") its pressure was "<<pressureRigid<<std::endl;
                }
                i_node->Reset(ISOLATED);   //reset isolated
                i_node->Reset(NEW_ENTITY); //reset if was new
                i_node->Reset(TO_REFINE);  //reset if was labeled to refine (to not duplicate boundary conditions)
                i_node->Reset(BLOCKED);

                if (i_node->IsNot(TO_ERASE))
                {

                  (i_mp->Nodes()).push_back(*(i_node.base()));
                  (rModelPart.Nodes()).push_back(*(i_node.base()));
                  rModelPart.Nodes().back().SetId(nodeId);
                  nodeId += 1;
                }
              }
              else
              {

                i_node->Set(ISOLATED);
                i_node->Reset(BOUNDARY);
                i_node->Reset(NEW_ENTITY); //reset if was new
                i_node->Reset(TO_REFINE);  //reset if was labeled to refine (to not duplicate boundary conditions)
                i_node->Reset(BLOCKED);

                if (mOptions.Is(MesherUtilities::KEEP_ISOLATED_NODES) && i_node->IsNot(TO_ERASE))
                {
                  i_node->FastGetSolutionStepValue(PRESSURE) = 0;
                  (i_mp->Nodes()).push_back(*(i_node.base()));
                  (rModelPart.Nodes()).push_back(*(i_node.base()));
                  rModelPart.Nodes().back().SetId(nodeId);
                  nodeId += 1;
                }
              }

            }
          }

          // for(ModelPart::ConditionsContainerType::iterator i_cond = i_mp->ConditionsBegin() ; i_cond != i_mp->ConditionsEnd() ; i_cond++)
          //   {
          // 	if( i_cond->IsNot(TO_ERASE) ){
          // 	  i_cond->Reset(NEW_ENTITY); //reset here if the condition is inserted
          // 	  PreservedConditions.push_back(*(i_cond.base()));
          // 	  PreservedConditions.back().SetId(condId);
          // 	  condId+=1;

          // 	  Geometry< Node<3> >& rGeometry = i_cond->GetGeometry();
          // 	  unsigned int NumNodes=rGeometry.size();
          // 	  unsigned int freeSurfaceNodes=0;
          // 	  unsigned int rigidNodes=0;
          // 	  for (unsigned int n = 0; n < NumNodes; ++n)
          // 	    {

          // 	      if(rGeometry[n].Is(RIGID) || rGeometry[n].Is(SOLID)){
          // 	  	// std::cout<<"rigid node! "<<rGeometry[n].X()<<" "<<rGeometry[n].Y()<<std::endl;
          // 	  	rigidNodes++;
          // 	      }else {
          // 	  	freeSurfaceNodes++;
          // 	      }
          // 	    }
          // 	  if((freeSurfaceNodes>0 && rigidNodes>0) || rigidNodes==0){
          // 	    for (unsigned int n = 0; n < NumNodes; ++n)
          // 	      {
          // 		rGeometry[n].Set(FREE_SURFACE);
          // 	      }
          // 	  }

          // 	}
          //   }
        }

        else
        {
          // Skipping the Computing Model Part
          if (!((i_mp->Is(ACTIVE) && i_mp->Is(SOLID)) || (i_mp->Is(ACTIVE) && i_mp->Is(FLUID))))
          {
            if (i_mp->NumberOfElements())
            {
              // Adding the remaining elements to the Main Mode Part
              for (ModelPart::ElementsContainerType::iterator i_elem = i_mp->ElementsBegin(); i_elem != i_mp->ElementsEnd(); i_elem++)
              {
                i_elem->SetId(elemId);
                (rModelPart.Elements()).push_back(*(i_elem.base()));
                rModelPart.Elements().back().SetId(elemId);
                elemId += 1;
              }
            }
          }
        }
      }

      // this->BuildBoundaryModelParts(rModelPart,PreservedConditions, nodeId, elemId, condId);

      // for(ModelPart::SubModelPartIterator i_mp= rModelPart.SubModelPartsBegin() ; i_mp!=rModelPart.SubModelPartsEnd(); i_mp++)
      // 	{
      // 	  if( i_mp->Is(BOUNDARY) ){ //boundary model part

      // 	    for(ModelPart::ConditionsContainerType::iterator i_cond = i_mp->ConditionsBegin() ; i_cond != i_mp->ConditionsEnd() ; i_cond++)
      // 	      {
      // 		if( i_cond->IsNot(TO_ERASE) ){
      // 		  i_cond->Reset(NEW_ENTITY); //reset here if the condition is inserted
      // 		  PreservedConditions.push_back(*(i_cond.base()));
      // 		  PreservedConditions.back().SetId(condId);
      // 		  condId+=1;
      // 		}
      // 	      }
      // 	  }

      // 	}

      // for(ModelPart::ConditionsContainerType::iterator i_cond = rModelPart.ConditionsBegin(); i_cond!= rModelPart.ConditionsEnd(); i_cond++)
      // 	{
      // 	  if(i_cond->Is(CONTACT)){
      // 	    PreservedConditions.push_back(*(i_cond.base()));
      // 	    PreservedConditions.back().SetId(condId);
      // 	    condId+=1;
      // 	  }
      // 	}

      // rModelPart.Conditions().swap(PreservedConditions);

      //Sort
      rModelPart.Nodes().Sort();
      rModelPart.Elements().Sort();
      // rModelPart.Conditions().Sort();

      //Unique
      rModelPart.Nodes().Unique();
      rModelPart.Elements().Unique();
      // rModelPart.Conditions().Unique();

      //Sort Again to have coherent numeration for nodes (mesh with shared nodes)
      unsigned int consecutive_index = 1;
      for (ModelPart::NodesContainerType::iterator in = rModelPart.NodesBegin(); in != rModelPart.NodesEnd(); in++)
        in->SetId(consecutive_index++);

      this->BuildComputingDomain(rModelPart, EchoLevel);

      if (EchoLevel > 1)
        std::cout << "   [ END MODEL PART [" << rModelPart.Name() << "] [Elems=:" << rModelPart.NumberOfElements() << "|Nodes=" << rModelPart.NumberOfNodes() << "|Conds=" << rModelPart.NumberOfConditions() << "] ] " << std::endl;

      KRATOS_CATCH(" ")
    }

    //*******************************************************************************************
    //*******************************************************************************************

    void BuildComputingDomain(ModelPart &rModelPart, int EchoLevel) override
    {
      KRATOS_TRY

      std::string ComputingModelPartName;
      for (ModelPart::SubModelPartIterator i_mp = rModelPart.SubModelPartsBegin(); i_mp != rModelPart.SubModelPartsEnd(); i_mp++)
      {
        // if( (i_mp->Is(ACTIVE) && i_mp->Is(SOLID)) ){ //solid_computing_domain
        if ((i_mp->Is(ACTIVE) && i_mp->Is(SOLID)) || (i_mp->Is(ACTIVE) && i_mp->Is(FLUID)))
        { // solid_computing_domain and fluid_computing_domain
          ComputingModelPartName = i_mp->Name();
        }
      }

      ModelPart &rComputingModelPart = rModelPart.GetSubModelPart(ComputingModelPartName);

      rComputingModelPart.Nodes().clear();
      rComputingModelPart.Elements().clear();
      // rComputingModelPart.Conditions().clear();

      for (ModelPart::SubModelPartIterator i_mp = rModelPart.SubModelPartsBegin(); i_mp != rModelPart.SubModelPartsEnd(); i_mp++)
      {
        // if( (i_mp->Is(SOLID) && i_mp->IsNot(ACTIVE)) || (i_mp->Is(BOUNDARY) && i_mp->Is(RIGID)) ){
        if (((i_mp->Is(ACTIVE) && i_mp->Is(SOLID)) || (i_mp->Is(BOUNDARY) && i_mp->Is(RIGID))) || (i_mp->Is(FLUID) && i_mp->IsNot(ACTIVE)))
        {

          for (ModelPart::NodesContainerType::iterator i_node = i_mp->NodesBegin(); i_node != i_mp->NodesEnd(); i_node++)
          {
            rComputingModelPart.Nodes().push_back(*(i_node.base()));
            // rComputingModelPart.AddNode(*(i_node.base())); // very slow!
          }

          // for(ModelPart::ConditionsContainerType::iterator i_cond = i_mp->ConditionsBegin() ; i_cond != i_mp->ConditionsEnd() ; i_cond++)
          //   {
          // 	rComputingModelPart.AddCondition(*(i_cond.base()));
          //   }

          for (ModelPart::ElementsContainerType::iterator i_elem = i_mp->ElementsBegin(); i_elem != i_mp->ElementsEnd(); i_elem++)
          {
            rComputingModelPart.AddElement(*(i_elem.base()));
          }
        }
      }

      //Sort
      rComputingModelPart.Nodes().Sort();
      rComputingModelPart.Elements().Sort();
      // rComputingModelPart.Conditions().Sort();

      //Unique
      rComputingModelPart.Nodes().Unique();
      rComputingModelPart.Elements().Unique();
      // rComputingModelPart.Conditions().Unique();

      if (EchoLevel > 1)
        std::cout << "    [ SUBMODEL PART [" << rComputingModelPart.Name() << "] [Elems=" << rComputingModelPart.NumberOfElements() << "|Nodes=" << rComputingModelPart.NumberOfNodes() << "|Conds=" << rComputingModelPart.NumberOfConditions() << "] ] " << std::endl;

      KRATOS_CATCH(" ")
    }

    void PerformModelSearches() override
    {

      KRATOS_TRY

      //NODAL NEIGHBOURS SEARCH
      NodalNeighboursSearchProcess FindNeighbours(mrMainModelPart);
      FindNeighbours.Execute();

      //BOUNDARY NORMALS SEARCH and SHRINKAGE FACTOR
      BoundaryNormalsCalculationUtilities BoundaryComputation;
      BoundaryComputation.CalculateWeightedBoundaryNormals(mrMainModelPart, mEchoLevel);

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
    ModelStartEndMeshingForFluidsProcess &operator=(ModelStartEndMeshingForFluidsProcess const &rOther);

    /// Copy constructor.
    //ModelStartEndMeshingForFluidsProcess(ModelStartEndMeshingForFluidsProcess const& rOther);

    ///@}

  }; // Class ModelStartEndMeshingForFluidsProcess

  ///@}

  ///@name Type Definitions
  ///@{

  ///@}
  ///@name Input and output
  ///@{

  /// input stream function
  inline std::istream &operator>>(std::istream &rIStream,
                                  ModelStartEndMeshingForFluidsProcess &rThis);

  /// output stream function
  inline std::ostream &operator<<(std::ostream &rOStream,
                                  const ModelStartEndMeshingForFluidsProcess &rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }
  ///@}

} // namespace Kratos.

#endif // KRATOS_MODEL_START_END_MESHING_FOR_FLUIDS_PROCESS_H_INCLUDED  defined
