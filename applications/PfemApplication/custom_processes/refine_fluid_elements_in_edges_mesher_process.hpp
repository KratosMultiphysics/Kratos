//
//   Project Name:        KratosPfemFluidApplication $
//   Created by:          $Author:       JMCarbonell $
//   Last modified by:    $Co-Author:                $
//   Date:                $Date:           July 2018 $
//   Revision:            $Revision:             0.0 $
//
//

#if !defined(KRATOS_REFINE_FLUID_ELEMENTS_IN_EDGES_MESHER_PROCESS_H_INCLUDED )
#define  KRATOS_REFINE_FLUID_ELEMENTS_IN_EDGES_MESHER_PROCESS_H_INCLUDED


// External includes

// System includes

// Project includes
#include "custom_processes/refine_elements_in_edges_mesher_process.hpp"

namespace Kratos
{

///@name Kratos Classes
///@{

/// Refine Mesh Elements Process 2D and 3D
/** The process inserts nodes in the rigid edge elements to be splitted in the remeshing process
*/

class RefineFluidElementsInEdgesMesherProcess
    : public RefineElementsInEdgesMesherProcess
{
 public:
  ///@name Type Definitions
  ///@{

  /// Pointer definition of Process
  KRATOS_CLASS_POINTER_DEFINITION( RefineFluidElementsInEdgesMesherProcess );

  typedef ModelPart::ConditionType         ConditionType;
  typedef ModelPart::PropertiesType       PropertiesType;
  typedef ConditionType::GeometryType       GeometryType;

  typedef RefineElementsInEdgesMesherProcess    BaseType;

  ///@}
  ///@name Life Cycle
  ///@{

  /// Default constructor.
  RefineFluidElementsInEdgesMesherProcess(ModelPart& rModelPart,
                                          MesherUtilities::MeshingParameters& rRemeshingParameters,
                                          int EchoLevel) : BaseType(rModelPart,rRemeshingParameters,EchoLevel)
  {
  }


  /// Destructor.
  virtual ~RefineFluidElementsInEdgesMesherProcess() {}


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
    return "RefineFluidElementsInEdgesMesherProcess";
  }

  /// Print information about this object.
  void PrintInfo(std::ostream& rOStream) const override
  {
    rOStream << "RefineFluidElementsInEdgesMesherProcess";
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
  ///@}
  ///@name Protected Operators
  ///@{
  ///@}
  ///@name Protected Operations
  ///@{
  ///@}


  //**************************************************************************
  //**************************************************************************

  void SelectFullBoundaryEdgedElements(ModelPart& rModelPart,
                                       ModelPart::ElementsContainerType& rBoundaryEdgedElements) override
  {
    KRATOS_TRY

    bool is_full_rigid_boundary = false;
    bool is_full_fluid_boundary = false;
    for(ModelPart::ElementsContainerType::iterator i_elem = rModelPart.ElementsBegin();
        i_elem != rModelPart.ElementsEnd(); ++i_elem)
    {
      Geometry< Node<3> >& rGeometry = i_elem->GetGeometry();

      is_full_rigid_boundary = true;
      for(unsigned int i=0; i<rGeometry.size(); ++i)
      {
        if( rGeometry[i].IsNot(BOUNDARY) || rGeometry[i].IsNot(RIGID) ){
          is_full_rigid_boundary = false;
          break;
        }

      }

      is_full_fluid_boundary = true;
      for(unsigned int i=0; i<rGeometry.size(); ++i)
      {
        if( rGeometry[i].IsNot(FREE_SURFACE) || rGeometry[i].Is(RIGID) || rGeometry[i].Is(SOLID) ){
          is_full_fluid_boundary = false;
          break;
        }

        // this condition fails
        // if ( rGeometry[i].GetValue(NEIGHBOUR_NODES).size() > rGeometry.size()-1 ){
        //   is_full_fluid_boundary = false;
        //   break;
        // }

        WeakPointerVector<Node<3> >& rN = rGeometry[i].GetValue(NEIGHBOUR_NODES);
        for(unsigned int j = 0; j < rN.size(); ++j)
        {
          if( rN[j].Is(SOLID) || rN[j].Is(RIGID) ){
            is_full_fluid_boundary = false;
            break;
          }
        }
      }
      // if( is_full_rigid_boundary )
      //   std::cout<<" is full rigid boundary "<<std::endl;

      // if( is_full_fluid_boundary )
      //   std::cout<<" is full fluid boundary "<<std::endl;

      if( is_full_rigid_boundary || is_full_fluid_boundary ){
        rBoundaryEdgedElements.push_back(*(i_elem.base()));
      }

    }

    KRATOS_CATCH( "" )
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
  ///@name Private Static Member Variables
  ///@{
  ///@}
  ///@name Private Static Member Variables
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
  RefineFluidElementsInEdgesMesherProcess& operator=(RefineFluidElementsInEdgesMesherProcess const& rOther);


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
                                  RefineFluidElementsInEdgesMesherProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const RefineFluidElementsInEdgesMesherProcess& rThis)
{
  rThis.PrintInfo(rOStream);
  rOStream << std::endl;
  rThis.PrintData(rOStream);

  return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_REFINE_FLUID_ELEMENTS_IN_EDGES_MESHER_PROCESS_H_INCLUDED  defined
