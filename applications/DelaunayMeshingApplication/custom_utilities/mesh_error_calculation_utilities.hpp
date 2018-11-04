//
//   Project Name:        KratosDelaunayMeshingApplication $
//   Created by:          $Author:             JMCarbonell $
//   Last modified by:    $Co-Author:                      $
//   Date:                $Date:                April 2018 $
//   Revision:            $Revision:                   0.0 $
//
//

#if !defined(KRATOS_MESH_ERROR_CALCULATION_UTILITIES_H_INCLUDED )
#define  KRATOS_MESH_ERROR_CALCULATION_UTILITIES_H_INCLUDED


// System includes

// External includes
#include "utilities/timer.h"

// Project includes
#include "includes/variables.h"
#include "includes/model_part.h"
#include "custom_utilities/mesher_utilities.hpp"

#include "delaunay_meshing_application_variables.h"

///VARIABLES used:
//Data:     NEIGHBOUR_ELEMENTS
//StepData:
//Flags:    (checked) NEW_ENTITY
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
  /** Computes mesh error for a given variable. Setting the error to elements or nodes for 2D and 3D
   */
  class MeshErrorCalculationUtilities
  {
  public:

    ///@name Type Definitions
    ///@{

    typedef ModelPart::ElementsContainerType                  ElementsContainerType;
    typedef ModelPart::NodesContainerType                        NodesContainerType;
    typedef ModelPart::MeshType::GeometryType::PointsArrayType      PointsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MeshErrorCalculationUtilities(){};

    /// Destructor.
    ~MeshErrorCalculationUtilities(){};


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    //**************************************************************************
    //**************************************************************************


    void NodalErrorCalculation(ModelPart& rModelPart,std::vector<double>& rNodalError,std::vector<int>& rIds,const Variable<double> &rVariable)
    {
      KRATOS_TRY

      std::vector<double> ElementalError;
      std::vector<int>    elems_ids;

      ElementalErrorCalculation(rModelPart,ElementalError,elems_ids,rVariable);

      if( !rNodalError.size() )
	rNodalError.resize(rModelPart.NumberOfNodes()+1);

      std::fill( rNodalError.begin(), rNodalError.end(), 100 );

      if( !rIds.size() )
	rIds.resize(MesherUtilities::GetMaxNodeId(rModelPart)+1);

      std::fill( rIds.begin(), rIds.end(), 0 );

      double NodalMeanError  = 0;

      int id=1;
      for(ModelPart::NodesContainerType::iterator in = rModelPart.NodesBegin() ; in != rModelPart.NodesEnd() ; ++in)
	{
	  if( in->IsNot(NEW_ENTITY) ){// && in->IsNot(STRUCTURE)){

	    WeakPointerVector<Element >& neighb_elems = in->GetValue(NEIGHBOUR_ELEMENTS);

	    NodalMeanError  = 0;


	    for(WeakPointerVector< Element >::iterator ne = neighb_elems.begin(); ne!=neighb_elems.end(); ++ne)
	      {
		NodalMeanError += ElementalError[elems_ids[ne->Id()]];
	      }

	    rIds[in->Id()]   = id;
	    rNodalError[id]  = NodalMeanError / double(neighb_elems.size());

	  }
	  else{

	    rIds[in->Id()]=id;
	    rNodalError[id] = 0;
	  }


	  id++;
	}

      KRATOS_CATCH( "" )
    }


    //**************************************************************************
    //**************************************************************************


    void ElementalErrorCalculation(ModelPart& rModelPart,std::vector<double>& rElementalError,std::vector<int>& rIds,const Variable<double> rVariable)
    {
      KRATOS_TRY

      ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();


      std::vector<double>  ElementVariable(MesherUtilities::GetMaxElementId(rModelPart)+1);
      std::fill( ElementVariable.begin(), ElementVariable.end(), 0 );

      std::vector<int> elems_ids;
      elems_ids.resize(MesherUtilities::GetMaxElementId(rModelPart)+1); //mesh 0
      std::fill( elems_ids.begin(), elems_ids.end(), 0 );

      double VariableMax = std::numeric_limits<double>::min();
      double VariableMin = std::numeric_limits<double>::max();

      std::vector<double> Value(1);

      unsigned int id=1;
      for(ModelPart::ElementsContainerType::const_iterator ie = rModelPart.ElementsBegin(); ie != rModelPart.ElementsEnd(); ++ie)
	{
	  (ie)->GetValueOnIntegrationPoints(rVariable,Value,CurrentProcessInfo);

	  elems_ids[ie->Id()] = id;
	  ElementVariable[id] = Value[0];

	  if(ElementVariable[id]>VariableMax)
	    VariableMax = ElementVariable[id];

	  if(ElementVariable[id]<VariableMin)
	    VariableMin = ElementVariable[id];

	  // if( ElementVariable[id] == 0 )
	  //   std::cout<<" ELEMENT ["<<ie->Id()<<"] : "<<ElementVariable[id]<<std::endl;

	  id++;
	}


      std::vector<double> NodalError(rModelPart.NumberOfNodes()+1);
      std::fill( NodalError.begin(), NodalError.end(), 0 );

      std::vector<int> nodes_ids (MesherUtilities::GetMaxNodeId(rModelPart)+1); //mesh 0
      std::fill( nodes_ids.begin(), nodes_ids.end(), 0 );

      double PatchSize  = 0;
      double PatchError = 0;

      double Size  = 0;
      double Error = 0;

      id=1;
      for(ModelPart::NodesContainerType::const_iterator in = rModelPart.NodesBegin(); in!=rModelPart.NodesEnd(); ++in)
	{

	  if(in->IsNot(NEW_ENTITY) ){// && in->IsNot(STRUCTURE)){

	    PatchSize  = 0;
	    PatchError = 0;

	    WeakPointerVector<Element >& neighb_elems = in->GetValue(NEIGHBOUR_ELEMENTS);


	    for(WeakPointerVector< Element >::iterator ne = neighb_elems.begin(); ne!=neighb_elems.end(); ++ne)
	      {

		Geometry<Node<3> >& pGeom = (ne)->GetGeometry();

		Size   = pGeom.DomainSize();  //Area(); or Volume();
		Error  = ElementVariable[elems_ids[ne->Id()]] * Size;

		PatchSize  += Size;
		PatchError += Error;

	      }

	    if(PatchSize!=0){
	      nodes_ids[in->Id()]=id;
	      NodalError[id] = PatchError/PatchSize;
	      //std::cout<<" Node ["<<in->Id()<<"] : ( NodalError: "<<NodalError[id]<<", PatchError: "<<" PatchSize:"<<PatchSize<<std::endl;
	    }
	    else{
	      std::cout<<" WARNING : Size surrounding node: "<<in->Id()<<" is null "<<std::endl;
	    }

	  }
	  else{

	    nodes_ids[in->Id()]=id;
	    NodalError[id]=0;

	  }

	  id++;


	}


      rElementalError.resize(MesherUtilities::GetMaxElementId(rModelPart)+1);
      std::fill( rElementalError.begin(), rElementalError.end(), 100 );

      rIds.resize(MesherUtilities::GetMaxElementId(rModelPart)+1); //mesh 0
      std::fill( rIds.begin(), rIds.end(), 0 );

      double VariableVariation = VariableMax - VariableMin;

      if(VariableVariation == 0){
	VariableVariation = 1;
	std::cout<<"   WARNING: "<<rVariable<<" min-max errors are the same ( MinVar= "<<VariableMin<<", MaxVar= "<<VariableMax<<" )"<<std::endl;
      }
      else{

	if( GetEchoLevel() > 0 )
	  std::cout<<"   Variable errors ( MinVar= "<<VariableMin<<", MaxVar= "<<VariableMax<<" )"<<std::endl;

	id=1;
	for(ModelPart::ElementsContainerType::const_iterator ie = rModelPart.ElementsBegin(); ie != rModelPart.ElementsEnd(); ++ie)
	  {

	    PointsArrayType& vertices=ie->GetGeometry().Points();


	    PatchError = 0;

	    unsigned int NumberOfVertices =vertices.size();
	    for(unsigned int i=0; i<NumberOfVertices; ++i)
	      {
		PatchError += NodalError[nodes_ids[vertices[i].Id()]];
	      }

	    if(NumberOfVertices!=0){
	      PatchError /= double(NumberOfVertices);
	    }
	    else{
	      std::cout<<" ERROR ME: Number of Vertices of the Element: "<<ie->Id()<<" is null "<<std::endl;
	    }

	    rIds[ie->Id()]=id;
	    rElementalError[id] = fabs( ( PatchError - ElementVariable[id] ) / VariableVariation ) * 100;

	    id++;
	  }
      }
      KRATOS_CATCH( "" )

	};


    /**
     * level of echo for the error calculation
     */
    virtual void SetEchoLevel(int Level)
    {
      mEchoLevel = Level;
    }

    int GetEchoLevel()
    {
      return mEchoLevel;
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
    int mEchoLevel;

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
    ///@name Unaccessible methods
    ///@{

    ///@}


  }; // Class MeshErrorCalculationUtilities

  ///@}
  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  ///@}

} // namespace Kratos.

#endif // KRATOS_MESH_ERROR_CALCULATION_UTILITIES_H_INCLUDED defined
