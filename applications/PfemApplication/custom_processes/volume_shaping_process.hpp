//
//   Project Name:        KratosPfemApplication     $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:        August 2018 $
//   Revision:            $Revision:            0.0 $
//
//


#if !defined(KRATOS_VOLUME_SHAPING_PROCESS_H_INCLUDED )
#define  KRATOS_VOLUME_SHAPING_PROCESS_H_INCLUDED

// External includes

// System includes

// Project includes
#include "includes/model_part.h"
#include "custom_utilities/mesher_utilities.hpp"
#include "processes/process.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/// Move free surface to restore volume losses
/** The process moves free surface nodes artificially in order to recover the loss of volume
 */

class VolumeShapingProcess : public Process
{
 public:
  ///@name Type Definitions
  ///@{

  /// Pointer definition of Process
  KRATOS_CLASS_POINTER_DEFINITION( VolumeShapingProcess );

  typedef ModelPart::NodeType                   NodeType;
  typedef ModelPart::ConditionType         ConditionType;
  typedef ModelPart::PropertiesType       PropertiesType;
  typedef ConditionType::GeometryType       GeometryType;

  ///@}
  ///@name Life Cycle
  ///@{

  /// Default constructor.
  VolumeShapingProcess(ModelPart& rModelPart, Parameters rParameters)
      : Process(Flags()), mrModelPart(rModelPart)
  {

    Parameters default_parameters( R"(
        {
             "variable_name": "VOLUME_WEAR",
             "flags_list": [],
             "properties": {}
        }  )" );


    // Validate against defaults -- this ensures no type mismatch
    rParameters.ValidateAndAssignDefaults(default_parameters);

    mVariableName = rParameters["variable_name"].GetString();

    // check variable to store volume shaping magnitude
    if( rModelPart.GetNodalSolutionStepVariablesList().Has( KratosComponents< Variable<double> >::Get( mVariableName ) ) == false )
    {
      KRATOS_ERROR << "trying to set a variable that is not in the model_part - variable name is " << mVariableName << std::endl;
    }

    // set properties
    Parameters variables = rParameters["properties"];
    for(auto iter = variables.begin(); iter != variables.end(); ++iter) {

      std::string variable_name = iter.name();
      const Parameters value = variables.GetValue(variable_name);

      if(KratosComponents<Variable<double> >::Has(variable_name)) {
        const Variable<double>& variable = KratosComponents<Variable<double>>().Get(variable_name);
        mProperties.SetValue(variable, value.GetDouble());
      } else if(KratosComponents<Variable<bool> >::Has(variable_name)) {
        const Variable<bool>& variable = KratosComponents<Variable<bool>>().Get(variable_name);
        mProperties.SetValue(variable, value.GetBool());
      } else if(KratosComponents<Variable<int> >::Has(variable_name)) {
        const Variable<int>& variable = KratosComponents<Variable<int>>().Get(variable_name);
        mProperties.SetValue(variable, value.GetInt());
      } else if(KratosComponents<Variable<std::string> >::Has(variable_name)) {
        const Variable<std::string>& variable = KratosComponents<Variable<std::string>>().Get(variable_name);
        mProperties.SetValue(variable, value.GetString());
      } else {
        KRATOS_ERROR << "Value type for Variable: "<<variable_name<<" not defined";
      }
    }

    for(unsigned int i=0; i<rParameters["flags_list"].size(); ++i)
    {
      mControlFlags.push_back(KratosComponents< Flags >::Get( rParameters["flags_list"][i].GetString() ));
    }

  }


  /// Destructor.
  virtual ~VolumeShapingProcess() {}


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

  /// Execute method is used to execute the AssignPropertiesToNodesProcess algorithms.
  void Execute()  override
  {
  }

  /// this function is designed for being called at the beginning of the computations
  /// right after reading the model and the groups
  void ExecuteInitialize() override
  {
  }

  /// this function is designed for being execute once before the solution loop but after all of the
  /// solvers where built
  void ExecuteBeforeSolutionLoop() override
  {
  }

  /// this function will be executed at every time step BEFORE performing the solve phase
  void ExecuteInitializeSolutionStep() override
  {
  }

  /// this function will be executed at every time step AFTER performing the solve phase
  void ExecuteFinalizeSolutionStep() override
  {
    KRATOS_TRY

    //shape volume
    this->ShapeVolume(mrModelPart);

    KRATOS_CATCH("")
  }


  /// this function will be executed at every time step BEFORE  writing the output
  void ExecuteBeforeOutputStep() override
  {
  }


  /// this function will be executed at every time step AFTER writing the output
  void ExecuteAfterOutputStep() override
  {
  }


  /// this function is designed for being called at the end of the computations
  /// right after reading the model and the groups
  void ExecuteFinalize() override
  {
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
    return "VolumeShapingProcess";
  }

  /// Print information about this object.
  void PrintInfo(std::ostream& rOStream) const override
  {
    rOStream << "VolumeShapingProcess";
  }

  /// Print object's data.s
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

  ModelPart& mrModelPart;

  std::vector<Flags> mControlFlags;

  Properties mProperties;

  std::string mVariableName;

  struct ShapingVariables{

    double TotalVolume;
    double TotalVolumeLoss;
    double TotalSurface;
    double OffsetFactor;

    std::vector<unsigned int>  VisitedNodesIds;
    std::vector<double> NodalSurface;
    std::vector<double> NodalVolume;
    std::vector<double> NodalVolumeLoss;
  };

  ///@}
  ///@name Private Operators
  ///@{
  ///@}
  ///@name Private Operations
  ///@{

  //*******************************************************************************************
  //*******************************************************************************************

  void ShapeVolume(ModelPart& rModelPart)
  {

    KRATOS_TRY

    double Tolerance = 1e-14;
    unsigned int NumberOfIterations = 5;
    unsigned int iteration = -1;

    this->SetVisitedEntities(rModelPart);

    ShapingVariables Variables;

    Variables.TotalVolume  = this->ComputeTotalVolume(rModelPart);
    Variables.TotalSurface = this->ComputeTotalSurface(rModelPart);

    this->ComputeNodalVolumeAndSurface(rModelPart,Variables);

    Variables.TotalVolumeLoss = this->ComputeVolumeLosses(rModelPart,Variables);

    Variables.OffsetFactor = 1;

    double VolumeIncrement = 0;

    double Ratio = fabs(Variables.TotalVolume-VolumeIncrement/(Variables.TotalVolume-Variables.TotalVolumeLoss));
    double Error = 1;

    while(++iteration < NumberOfIterations && ( Ratio < 0.9 || Ratio > 1 ) && Error > Tolerance)
    {

      this->MoveSurface(rModelPart,Variables);

      VolumeIncrement = Variables.TotalVolume - this->ComputeTotalVolume(rModelPart);

      if( VolumeIncrement > 0 )
        Variables.OffsetFactor*=(1.0-(Variables.TotalVolumeLoss/VolumeIncrement));

      if( (VolumeIncrement-Variables.TotalVolumeLoss) * Variables.OffsetFactor > 0 )
        Variables.OffsetFactor*=(-1);

      Ratio = fabs(Variables.TotalVolume-VolumeIncrement/(Variables.TotalVolume-Variables.TotalVolumeLoss));

      Error = fabs(Variables.TotalVolumeLoss-VolumeIncrement);
      //std::cout<<" Iteration: "<<iteration<<" Ratio in Volume "<<Ratio<<" Volume Increment "<<VolumeIncrement<<std::endl;
    }


    this->ResetVisitedEntities(rModelPart);

    std::cout<<" Surface shaping perfomed in "<<iteration<<" iterations : Error "<<Error<<std::endl;
    std::cout<<" TotalVolume "<<Variables.TotalVolume<<" TotalVolumeLoss "<<Variables.TotalVolumeLoss<<" VolumeIncrement "<<VolumeIncrement<<std::endl;

    KRATOS_CATCH("")
  }


  //*******************************************************************************************
  //*******************************************************************************************

  void MoveSurface(ModelPart& rModelPart, const ShapingVariables& rVariables)
  {
    KRATOS_TRY

    ModelPart::NodesContainerType::iterator it_begin = rModelPart.NodesBegin();
    int NumberOfNodes = rModelPart.NumberOfNodes();

    #pragma omp parallel for
    for (int i=0; i<NumberOfNodes; ++i)
    {
      ModelPart::NodesContainerType::iterator i_node = it_begin + i;
      if(i_node->Is(VISITED)){
        unsigned int id = rVariables.VisitedNodesIds[i_node->Id()];
        const array_1d<double,3>& rNormal = i_node->FastGetSolutionStepValue(NORMAL);
        double& rShrinkFactor = i_node->FastGetSolutionStepValue(SHRINK_FACTOR);
        double rOffset = -rVariables.OffsetFactor * rShrinkFactor * rVariables.NodalVolumeLoss[id] * rVariables.NodalSurface[id] / (100 * rVariables.TotalVolume * rVariables.TotalSurface);
        //std::cout<<" rOffset "<<rOffset<<" "<<rShrinkFactor<<" "<<rVariables.NodalVolumeLoss[id]<<" "<<rVariables.NodalSurface[id]<<" "<<rVariables.TotalSurface<<std::endl;

        i_node->Coordinates() += rOffset * rNormal;
        i_node->FastGetSolutionStepValue(DISPLACEMENT) += rOffset * rNormal;
        i_node->FastGetSolutionStepValue(DISPLACEMENT,1) += rOffset * rNormal;
      }
    }

    KRATOS_CATCH("")
  }

  //*******************************************************************************************
  //*******************************************************************************************

  void SetVisitedEntities(ModelPart& rModelPart)
  {
    KRATOS_TRY

    this->SetVisitedNodes(rModelPart);
    this->SetVisitedElements(rModelPart);

    KRATOS_CATCH("")
  }


  //*******************************************************************************************
  //*******************************************************************************************

  void ResetVisitedEntities(ModelPart& rModelPart)
  {
    KRATOS_TRY

    this->ResetVisitedNodes(rModelPart);
    this->ResetVisitedElements(rModelPart);

    KRATOS_CATCH("")
  }

  //*******************************************************************************************
  //*******************************************************************************************

  void SetVisitedNodes(ModelPart& rModelPart)
  {
    KRATOS_TRY

    for(ModelPart::NodesContainerType::iterator i_node = rModelPart.NodesBegin() ; i_node != rModelPart.NodesEnd() ; ++i_node)
    {
      if(this->MatchControlFlags(*(i_node.base()))){
        i_node->Set(VISITED,true);
      }
    }

    KRATOS_CATCH("")
  }


  //*******************************************************************************************
  //*******************************************************************************************

  void ResetVisitedNodes(ModelPart& rModelPart)
  {
    KRATOS_TRY

    for(ModelPart::NodesContainerType::iterator i_node = rModelPart.NodesBegin() ; i_node != rModelPart.NodesEnd() ; ++i_node)
    {
      i_node->Set(VISITED,false);
    }

    KRATOS_CATCH("")
  }


  //*******************************************************************************************
  //*******************************************************************************************

  void SetVisitedElements(ModelPart& rModelPart)
  {
    KRATOS_TRY

    for(ModelPart::ElementsContainerType::iterator i_elem = rModelPart.ElementsBegin() ; i_elem != rModelPart.ElementsEnd() ; ++i_elem)
    {
      for(unsigned int i=0; i<i_elem->GetGeometry().size(); ++i)
      {
        if(this->MatchControlFlags(i_elem->GetGeometry()(i))){
          i_elem->Set(VISITED,true);
          break;
        }
      }
    }

    KRATOS_CATCH("")
  }


  //*******************************************************************************************
  //*******************************************************************************************

  void ResetVisitedElements(ModelPart& rModelPart)
  {
    KRATOS_TRY

    for(ModelPart::ElementsContainerType::iterator i_elem = rModelPart.ElementsBegin() ; i_elem != rModelPart.ElementsEnd() ; ++i_elem)
    {
      i_elem->Set(VISITED,false);
    }

    KRATOS_CATCH("")
  }

  //*******************************************************************************************
  //*******************************************************************************************

  double ComputeTotalVolume(ModelPart& rModelPart)
  {
    KRATOS_TRY

    const unsigned int dimension = rModelPart.GetProcessInfo()[SPACE_DIMENSION];
    double TotalVolume = 0;

    if( dimension == 2 ){

      ModelPart::ElementsContainerType::iterator it_begin = rModelPart.ElementsBegin();
      int NumberOfElements = rModelPart.NumberOfElements();

      #pragma omp parallel for reduction(+:TotalVolume)
      for (int i=0; i<NumberOfElements; ++i)
      {
        ModelPart::ElementsContainerType::iterator i_elem = it_begin + i;
        if( i_elem->GetGeometry().Dimension() == 2 && i_elem->Is(VISITED) )
          TotalVolume += i_elem->GetGeometry().Area();
      }
    }
    else{ //dimension == 3

      ModelPart::ElementsContainerType::iterator it_begin = rModelPart.ElementsBegin();
      int NumberOfElements = rModelPart.NumberOfElements();

      #pragma omp parallel for reduction(+:TotalVolume)
      for (int i=0; i<NumberOfElements; ++i)
      {
        ModelPart::ElementsContainerType::iterator i_elem = it_begin + i;
        if( i_elem->GetGeometry().Dimension() == 3 && i_elem->Is(VISITED) )
          TotalVolume += i_elem->GetGeometry().Volume();
      }
    }

    return TotalVolume;

    KRATOS_CATCH("")

  }


  //*******************************************************************************************
  //*******************************************************************************************

  double ComputeTotalSurface(ModelPart& rModelPart)
  {
    KRATOS_TRY

    const unsigned int dimension = rModelPart.GetProcessInfo()[SPACE_DIMENSION];
    double TotalSurface = 0;

    if( dimension == 2 ){

      ModelPart::ElementsContainerType::iterator it_begin = rModelPart.ElementsBegin();
      int NumberOfElements = rModelPart.NumberOfElements();

      #pragma omp parallel for reduction(+:TotalSurface)
      for (int i=0; i<NumberOfElements; ++i)
      {
        ModelPart::ElementsContainerType::iterator i_elem = it_begin + i;
        GeometryType& rGeometry = i_elem->GetGeometry();
        if( rGeometry.Dimension() == 2 && i_elem->Is(VISITED) ){
          for(unsigned int j=0; j<rGeometry.size()-1; ++j)
          {
            if(rGeometry[j].Is(VISITED)){
              for(unsigned int k=j+1; k<rGeometry.size(); ++k)
              {
                if(rGeometry[k].Is(VISITED)){
                  TotalSurface += norm_2( rGeometry[k].Coordinates() - rGeometry[j].Coordinates() );
                }
              }

            }
          }
        }
      }
    }
    else{ //dimension == 3

      DenseMatrix<unsigned int> lpofa; //connectivities of points defining faces
      DenseVector<unsigned int> lnofa; //number of points defining faces

      ModelPart::ElementsContainerType::iterator it_begin = rModelPart.ElementsBegin();
      int NumberOfElements = rModelPart.NumberOfElements();

      #pragma omp parallel for private(lpofa,lnofa) reduction(+:TotalSurface)
      for (int i=0; i<NumberOfElements; ++i)
      {
        ModelPart::ElementsContainerType::iterator i_elem = it_begin + i;

        GeometryType& rGeometry = i_elem->GetGeometry();

        if( rGeometry.Dimension() == 3 && i_elem->Is(VISITED) ){

          rGeometry.NodesInFaces(lpofa);
          rGeometry.NumberNodesInFaces(lnofa);

          for(unsigned int iface=0; iface<rGeometry.FacesNumber(); ++iface)
          {
            unsigned int free_surface = 0;
            for(unsigned int j=1; j<=lnofa[iface]; ++j)
              if(rGeometry[j].Is(VISITED))
                ++free_surface;

            if(free_surface==lnofa[iface])
              TotalSurface+=Compute3DArea(rGeometry[lpofa(1,iface)].Coordinates(),
                                          rGeometry[lpofa(2,iface)].Coordinates(),
                                          rGeometry[lpofa(3,iface)].Coordinates());
          }
        }
      }
    }

    return TotalSurface;

    KRATOS_CATCH("")

  }


  //*******************************************************************************************
  //*******************************************************************************************

  void ComputeNodalVolumeAndSurface(ModelPart& rModelPart, ShapingVariables& rVariables)
  {
    KRATOS_TRY

    unsigned int MaxNodeId = MesherUtilities::GetMaxNodeId(rModelPart);
    rVariables.VisitedNodesIds.resize(MaxNodeId+1);
    std::fill( rVariables.VisitedNodesIds.begin(), rVariables.VisitedNodesIds.end(), 0 );

    std::vector<WeakPointerVector<Element> > Neighbours(rModelPart.NumberOfNodes()+1);
    unsigned int id = 1;
    for(ModelPart::ElementsContainerType::iterator i_elem = rModelPart.ElementsBegin(); i_elem!=rModelPart.ElementsEnd(); ++i_elem)
    {
      if(i_elem->Is(VISITED)){

        Element::GeometryType& pGeometry = i_elem->GetGeometry();

        for(unsigned int i = 0; i < pGeometry.size(); ++i)
        {
          if( pGeometry[i].Is(VISITED) ){

            if(rVariables.VisitedNodesIds[pGeometry[i].Id()]==0){
              rVariables.VisitedNodesIds[pGeometry[i].Id()]=id;
              Neighbours[id].push_back( Element::WeakPointer( *(i_elem.base()) ) );
              ++id;
            }
            else{
              Neighbours[rVariables.VisitedNodesIds[pGeometry[i].Id()]].push_back( Element::WeakPointer( *(i_elem.base()) ) );
            }

          }
        }
      }
    }

    rVariables.NodalSurface.resize(id);
    std::fill( rVariables.NodalSurface.begin(), rVariables.NodalSurface.end(), 0 );

    rVariables.NodalVolume.resize(id);
    std::fill( rVariables.NodalVolume.begin(), rVariables.NodalVolume.end(), 0 );

    const unsigned int dimension = rModelPart.GetProcessInfo()[SPACE_DIMENSION];
    double TotalSurface = 0;

    if( dimension == 2 ){

      ModelPart::NodesContainerType::iterator it_begin = rModelPart.NodesBegin();
      int NumberOfNodes = rModelPart.NumberOfNodes();

      #pragma omp parallel for reduction(+:TotalSurface)
      for (int i=0; i<NumberOfNodes; ++i)
      {
        ModelPart::NodesContainerType::iterator i_node = it_begin + i;

        if( i_node->Is(VISITED) ){
          unsigned int id = rVariables.VisitedNodesIds[i_node->Id()];
          WeakPointerVector<Element>& rE = Neighbours[id];
          for(WeakPointerVector<Element >::iterator ie= rE.begin(); ie!=rE.end(); ++ie)
          {
            GeometryType& rGeometry = ie->GetGeometry();

            rVariables.NodalVolume[id] += ie->GetGeometry().Area() / double(ie->GetGeometry().size());

            if( rGeometry.Dimension() == 2 ){
              for(unsigned int j=0; j<rGeometry.size()-1; ++j)
              {
                if(rGeometry[j].Is(VISITED)){
                  for(unsigned int k=j+1; k<rGeometry.size(); ++k)
                  {
                    if(rGeometry[k].Is(VISITED) && (rGeometry[j].Id() == i_node->Id() || rGeometry[k].Id() == i_node->Id()) ){
                      rVariables.NodalSurface[id] += 0.5 * norm_2( rGeometry[k].Coordinates() - rGeometry[j].Coordinates() ); //linear elements 2 noded sides
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    else{ //dimension == 3

      DenseMatrix<unsigned int> lpofa; //connectivities of points defining faces
      DenseVector<unsigned int> lnofa; //number of points defining faces

      ModelPart::NodesContainerType::iterator it_begin = rModelPart.NodesBegin();
      int NumberOfNodes = rModelPart.NumberOfNodes();

      double athird = 1.0/3.0;

      #pragma omp parallel for private(lpofa,lnofa)
      for (int i=0; i<NumberOfNodes; ++i)
      {
        ModelPart::NodesContainerType::iterator i_node = it_begin + i;

        if( i_node->Is(VISITED) ){
          unsigned int id = rVariables.VisitedNodesIds[i_node->Id()];
          WeakPointerVector<Element>& rE = Neighbours[id];
          for(WeakPointerVector<Element >::iterator ie= rE.begin(); ie!=rE.end(); ++ie)
          {
            GeometryType& rGeometry = ie->GetGeometry();

            rVariables.NodalVolume[id] += ie->GetGeometry().Volume() / double(ie->GetGeometry().size());

            if( rGeometry.Dimension() == 3 ){

              rGeometry.NodesInFaces(lpofa);
              rGeometry.NumberNodesInFaces(lnofa);

              for(unsigned int iface=0; iface<rGeometry.FacesNumber(); ++iface)
              {
                unsigned int visited = 0;
                bool selected = false;
                for(unsigned int j=1; j<=lnofa[iface]; ++j){
                  if(rGeometry[j].Is(VISITED))
                    ++visited;
                  if(rGeometry[i].Id() == i_node->Id())
                    selected = true;
                }

                if(visited==lnofa[iface] && selected)
                  rVariables.NodalSurface[id]+= athird * Compute3DArea(rGeometry[lpofa(1,iface)].Coordinates(),
                                                                       rGeometry[lpofa(2,iface)].Coordinates(),
                                                                       rGeometry[lpofa(3,iface)].Coordinates());
              }

            }
          }
        }
      }

    }

    KRATOS_CATCH("")

  }

  //*******************************************************************************************
  //*******************************************************************************************

  double ComputeVolumeLosses(ModelPart& rModelPart, ShapingVariables& rVariables)
  {
    KRATOS_TRY

    double TotalVolumeLoss = 0;

    Variable<double> VolumeLossVariable = KratosComponents< Variable<double> >::Get(mVariableName);

    double ArchardCoefficient = 0;
    if( mProperties.Has(WEAR_COEFFICIENT) && mProperties.Has(INDENTATION_HARDNESS) ){
      ArchardCoefficient = mProperties[WEAR_COEFFICIENT] / mProperties[INDENTATION_HARDNESS];
    }
    else{
      KRATOS_WARNING("Wear calculation not possible") << " WEAR_COEFFICIENT and INDENTATION_HARDNESS variables not DEFINED " << std::endl;
    }

    ModelPart::NodesContainerType::iterator it_begin = rModelPart.NodesBegin();
    int NumberOfNodes = rModelPart.NumberOfNodes();

    rVariables.NodalVolumeLoss.resize(rVariables.NodalVolume.size());
    std::fill( rVariables.NodalVolumeLoss.begin(), rVariables.NodalVolumeLoss.end(), 0 );

    #pragma omp parallel for reduction(+:TotalVolumeLoss)
    for (int i=0; i<NumberOfNodes; ++i)
    {
      ModelPart::NodesContainerType::iterator i_node = it_begin + i;

      unsigned int id = rVariables.VisitedNodesIds[i_node->Id()];
      if( i_node->Is(VISITED) && id != 0 ){

        const array_1d<double,3>& Normal = i_node->FastGetSolutionStepValue(NORMAL);
        const array_1d<double,3>& ContactForce = i_node->FastGetSolutionStepValue(CONTACT_FORCE);

        array_1d<double,3> DeltaDisplacement = i_node->FastGetSolutionStepValue(DISPLACEMENT) - i_node->FastGetSolutionStepValue(DISPLACEMENT,1);

        rVariables.NodalVolumeLoss[id] = ArchardCoefficient * fabs(inner_prod(ContactForce,Normal)) * norm_2( DeltaDisplacement - inner_prod(DeltaDisplacement,Normal) * Normal);

        i_node->FastGetSolutionStepValue(VolumeLossVariable) += rVariables.NodalVolumeLoss[id];

        if( rVariables.NodalVolumeLoss[id] > rVariables.NodalVolume[id] )
          rVariables.NodalVolumeLoss[id] = rVariables.NodalVolume[id];

        TotalVolumeLoss += rVariables.NodalVolumeLoss[id];

      }
    }

    return TotalVolumeLoss;

    KRATOS_CATCH("")
  }


  //*******************************************************************************************
  //*******************************************************************************************


  bool MatchControlFlags(const Node<3>::Pointer& pNode)
  {

    for(unsigned int i = 0; i<mControlFlags.size(); i++)
    {
      if( pNode->IsNot(mControlFlags[i]) )
        return false;
    }

    return true;

  }


  double Compute3DArea(array_1d<double,3> PointA, array_1d<double,3> PointB, array_1d<double,3> PointC){
    double a = MathUtils<double>::Norm3(PointA - PointB);
    double b = MathUtils<double>::Norm3(PointB - PointC);
    double c = MathUtils<double>::Norm3(PointC - PointA);
    double s = (a+b+c) / 2.0;
    double Area=std::sqrt(s*(s-a)*(s-b)*(s-c));
    return Area;
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
  VolumeShapingProcess& operator=(VolumeShapingProcess const& rOther);


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
                                  VolumeShapingProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const VolumeShapingProcess& rThis)
{
  rThis.PrintInfo(rOStream);
  rOStream << std::endl;
  rThis.PrintData(rOStream);

  return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_VOLUME_SHAPING_PROCESS_H_INCLUDED  defined
