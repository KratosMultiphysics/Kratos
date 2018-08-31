//
//   Project Name:        KratosDelaunayMeshingApplication $
//   Created by:          $Author:             JMCarbonell $
//   Last modified by:    $Co-Author:                      $
//   Date:                $Date:                April 2018 $
//   Revision:            $Revision:                   0.0 $
//
//

#if !defined( KRATOS_REFINE_ELEMENTS_ON_SIZE_MESHER_PROCESS_H_INCLUDED )
#define KRATOS_REFINE_ELEMENTS_ON_SIZE_MESHER_PROCESS_H_INCLUDED


// External includes

// System includes

// Project includes
#include "containers/variables_list_data_value_container.h"
#include "spatial_containers/spatial_containers.h"

#include "includes/model_part.h"
#include "custom_utilities/mesher_utilities.hpp"
#include "custom_processes/mesher_process.hpp"

///VARIABLES used:
//Data:
//StepData: NODAL_H, CONTACT_FORCE
//Flags:    (checked) TO_REFOME, BOUNDARY, NEW_ENTITY
//          (set)
//          (modified)
//          (reset)
//(set):=(set in this process)

namespace Kratos
{

///@name Kratos Classes
///@{

/// Refine Mesh Elements Process 2D and 3D
/** The process labels the elements to be refined in the mesher
    it applies a size constraint to elements that must be refined.

*/
class RefineElementsOnSizeMesherProcess
  : public MesherProcess
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION( RefineElementsOnSizeMesherProcess );

    typedef ModelPart::ConditionType         ConditionType;
    typedef ModelPart::PropertiesType       PropertiesType;
    typedef ConditionType::GeometryType       GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    RefineElementsOnSizeMesherProcess(ModelPart& rModelPart,
				    MesherUtilities::MeshingParameters& rRemeshingParameters,
				    int EchoLevel)
      : mrModelPart(rModelPart),
	mrRemesh(rRemeshingParameters)
    {
      mEchoLevel = EchoLevel;
    }


    /// Destructor.
    virtual ~RefineElementsOnSizeMesherProcess() {}


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
    void Execute() override
    {
      KRATOS_TRY

      if( mEchoLevel > 0 ){
	std::cout<<" [ SELECT ELEMENTS TO REFINE : "<<std::endl;
	//std::cout<<"   refine selection "<<std::endl;
      }


      //***SIZES :::: parameters do define the tolerance in mesh size:
      double size_for_inside_elements   = 0.75 * mrRemesh.Refine->CriticalRadius;
      double size_for_boundary_elements = 1.50 * mrRemesh.Refine->CriticalRadius;

      double nodal_h_refining_factor     = 0.75;
      double nodal_h_non_refining_factor = 2.00;

      ProcessInfo& CurrentProcessInfo = mrModelPart.GetProcessInfo();

      unsigned int refine_on_size      = 0;
      unsigned int refine_on_threshold = 0;

      int id = 0;
      if(mrRemesh.Refine->RefiningOptions.Is(MesherUtilities::REFINE_ELEMENTS)
	 && mrRemesh.Refine->RefiningOptions.Is(MesherUtilities::REFINE_ADD_NODES) )
	{

	  ModelPart::ElementsContainerType::iterator element_begin = mrModelPart.ElementsBegin();

	  unsigned int nds = (*element_begin).GetGeometry().size();

	  MesherUtilities::MeshContainer& InMesh = mrRemesh.InMesh;

	  InMesh.CreateElementList(mrRemesh.Info->NumberOfElements, nds); //number of preserved elements
	  InMesh.CreateElementSizeList(mrRemesh.Info->NumberOfElements);

	  int& OutNumberOfElements = mrRemesh.OutMesh.GetNumberOfElements();

	  int* InElementList        = mrRemesh.InMesh.GetElementList();
	  double* InElementSizeList = mrRemesh.InMesh.GetElementSizeList();

	  int* OutElementList       = mrRemesh.OutMesh.GetElementList();

	  ModelPart::NodesContainerType::iterator nodes_begin = mrModelPart.NodesBegin();

	  //PREPARE THE NODAL_H as a variable to control the automatic point insertion
	  //**************************************************************************

	  if(mrRemesh.Refine->RefiningOptions.IsNot(MesherUtilities::REFINE_BOUNDARY)){

	    for(unsigned int i = 0; i<mrModelPart.Nodes().size(); i++)
	      {
		////Assign a huge NODAL_H to the free surface nodes, so that there no nodes will be added
		// if ( (nodes_begin + i)->Is(FREE_SURFACE))
		// {
		// 	double & nodal_h=(nodes_begin + i)->FastGetSolutionStepValue(NODAL_H);
		// 	nodal_h*=nodal_h_non_refining_factor;
		// }

		//Assign a huge NODAL_H to the Boundary nodes, so that there no nodes will be added
		if ( (nodes_begin + i)->Is(BOUNDARY))
		  {
		    double & nodal_h=(nodes_begin + i)->FastGetSolutionStepValue(NODAL_H);
		    nodal_h*=nodal_h_non_refining_factor;
		  }


	      }

	  }

	  //SET THE REFINED ELEMENTS AND THE AREA (NODAL_H)
	  //*********************************************************************
	  mrRemesh.Info->CriticalElements = 0;

	  for(int el = 0; el< OutNumberOfElements; el++)
	    {
	      if(mrRemesh.PreservedElements[el])
		{

		  double prescribed_h      = 0;
		  bool   dissipative       = false;  //dissipative means reference threshold is overwhelmed
		  bool   refine_size       = false;

		  unsigned int    count_dissipative = 0;
		  unsigned int    count_boundary_inserted = 0;
		  unsigned int    count_boundary = 0;
		  unsigned int    count_contact_boundary = 0;

		  Geometry<Node<3> > vertices;

		  for(unsigned int pn=0; pn<nds; pn++)
		    {

		      InElementList[id*nds+pn]= OutElementList[el*nds+pn];

		      vertices.push_back(*(nodes_begin + OutElementList[el*nds+pn]-1).base());

		      prescribed_h += (nodes_begin + OutElementList[el*nds+pn]-1)->FastGetSolutionStepValue(NODAL_H);

		      if((nodes_begin + OutElementList[el*nds+pn]-1)->Is(TO_REFINE))
			count_dissipative+=1;

		      if((nodes_begin + OutElementList[el*nds+pn]-1)->Is(BOUNDARY)){

			count_boundary+=1;

			if((nodes_begin + OutElementList[el*nds+pn]-1)->Is(NEW_ENTITY))
			  count_boundary_inserted+=1;


			if( (nodes_begin + OutElementList[el*nds+pn]-1)->SolutionStepsDataHas(CONTACT_FORCE) ){
			  array_1d<double, 3 > & ContactForceNormal = (nodes_begin + OutElementList[el*nds+pn]-1)->FastGetSolutionStepValue(CONTACT_FORCE);
			  if( norm_2(ContactForceNormal) )
			    count_contact_boundary+=1;
			}

		      }

		    }

		  // if(count_dissipative>0)
		  //   std::cout<<" Count REFINE nodes "<<count_dissipative<<std::endl;
		  // std::cout<<"   prescribed_h (el:"<<el<<") = "<<prescribed_h<<std::endl;

		  bool refine_candidate = true;
		  if (mrRemesh.Refine->RefiningBoxSetFlag == true ){
		    refine_candidate = mMesherUtilities.CheckVerticesInBox(vertices,*(mrRemesh.Refine->RefiningBox),CurrentProcessInfo);
		  }



		  double element_size = 0;
		  double element_radius = mMesherUtilities.CalculateElementRadius(vertices, element_size);


		  //calculate the prescribed h
		  prescribed_h *= 0.3333;

		  double h = mrRemesh.AlphaParameter * prescribed_h;
		  double element_ideal_radius = 0;
		  if( nds == 3 ){ //if h is the height of a equilateral triangle, the area is sqrt(3)*h*h/4
		    element_ideal_radius = sqrt(3.0) * 0.25 * ( h * h );
		  }

		  if( nds == 4 ){//if h is the height of a regular tetrahedron, the volume is h*h*h/(6*sqrt(2))
		    //element_ideal_radius = (27.0/16.0)* sqrt(1.0/108.0) * ( h * h* h );
		    element_ideal_radius = ( h * h * h )/( 6.0 * sqrt(2.0) );
		  }


		  //std::cout<<"   prescribed_h (el:"<<el<<") = "<<prescribed_h<<std::endl;

		  if( refine_candidate ){


		    //********* PLASTIC POWER ENERGY REFINEMENT CRITERION (A)
		    if(count_dissipative>=nds-1){
		      dissipative = true;
		      //Set Critical Elements
		      mrRemesh.Info->CriticalElements += 1;
		      //std::cout<<" Dissipative True "<<std::endl;
		    }


		    //********* SIZE REFINEMENT CRITERION (B)

		    // if(dissipative)
		    //   std::cout<<" element_radius "<<element_radius<<" CriticalRadius "<<size_for_inside_elements<<std::endl;

		    double critical_size = size_for_inside_elements;
		    if(count_boundary>=nds-1)
		      critical_size = size_for_boundary_elements;


		    if( element_radius > critical_size ){
		      refine_size = true;
		    }


		    //Also a criteria for the CriticalDissipation (set in nodes)
		    if(mrRemesh.Refine->RefiningOptions.Is(MesherUtilities::REFINE_ELEMENTS_ON_THRESHOLD)
		       && mrRemesh.Refine->RefiningOptions.Is(MesherUtilities::REFINE_ELEMENTS_ON_DISTANCE))
		      {
			// if( dissipative )
			// 	std::cout<<" [ Refine Criteria ["<<id<<"] : (dissipative: "<<dissipative<<"; size: "<<refine_size<<"; prescribed h: "<<prescribed_h<<") ]"<<std::endl;


			//********* THRESHOLD REFINEMENT CRITERION (A)
			if( (dissipative == true && refine_size == true) )
			  {
			    InElementSizeList[id] = nodal_h_refining_factor * element_size;

			    refine_on_threshold += 1;
			    //std::cout<<" Area Factor Refine DISSIPATIVE :"<<InElementSizeList[id]<<std::endl;
			  }
			else if (refine_size == true)
			  {


			    InElementSizeList[id] = element_ideal_radius;

			    if( count_boundary_inserted && count_contact_boundary){

			      InElementSizeList[id] = nodal_h_refining_factor * element_ideal_radius;

			      std::cout<<" count boundary inserted-contact on "<<std::endl;
			    }

			    refine_on_size +=1;
			    //std::cout<<" Area Factor Refine SIZE :"<<InElementSizeList[id]<<std::endl;

			  }
			else{

			  InElementSizeList[id] = nodal_h_non_refining_factor * element_size;
			  //std::cout<<" Area Factor Refine NO :"<<InElementSizeList[id]<<std::endl;
			}


			//std::cout<<" [ AREA: "<<Area<<"; NEW AREA: "<<InElementSizeList[id]<<" ]"<<std::endl;
		      }
		    else if(mrRemesh.Refine->RefiningOptions.Is(MesherUtilities::REFINE_ELEMENTS_ON_DISTANCE))
		      {

			//********* SIZE REFINEMENT CRITERION (B)
			if( refine_size == true ){

			  InElementSizeList[id] = nodal_h_refining_factor * element_ideal_radius;
			}
			else{

			  //InElementSizeList[id] = element_size;
			  InElementSizeList[id] = element_ideal_radius;
			}

		      }
		    else{

		      //InElementSizeList[id] = element_ideal_radius;
		      InElementSizeList[id] = nodal_h_non_refining_factor * element_size;

		    }


		    //std::cout<<"   mod_prescribed_h (el:"<<el<<") = "<<prescribed_h<<" [ Triangle Area: "<<InElementSizeList[id]<<" ]"<<std::endl;

		  }
		  else{

		    //InElementSizeList[id] = element_ideal_radius;
		    InElementSizeList[id] = nodal_h_non_refining_factor * element_size;

		  }


		  id += 1;

		}



	    }



	  //RESTORE THE NODAL_H ON BOUNDARY
	  //*********************************************************************
	  if(mrRemesh.Refine->RefiningOptions.IsNot(MesherUtilities::REFINE_BOUNDARY)){

	    for(unsigned int i = 0; i<mrModelPart.Nodes().size(); i++)
	      {
		// //Unassign the NODAL_H of the free surface nodes
		// if ( (nodes_begin + i)->Is(FREE_SURFACE))
		// {
		// 	double & nodal_h=(nodes_begin + i)->FastGetSolutionStepValue(NODAL_H);
		// 	nodal_h/=nodal_h_non_refining_factor;
		// }

		//Unassign the NODAL_H of the Boundary nodes
		if ( (nodes_begin + i)->Is(BOUNDARY))
		  {
		    double & nodal_h=(nodes_begin + i)->FastGetSolutionStepValue(NODAL_H);
		    nodal_h/=nodal_h_non_refining_factor;
		  }


	      }
	  }
	}
      else{

	  ModelPart::ElementsContainerType::iterator element_begin = mrModelPart.ElementsBegin();

	  unsigned int nds = (*element_begin).GetGeometry().size();

	  MesherUtilities::MeshContainer& InMesh = mrRemesh.InMesh;

	  InMesh.CreateElementList(mrRemesh.Info->NumberOfElements, nds); //number of preserved elements
	  InMesh.CreateElementSizeList(mrRemesh.Info->NumberOfElements);

	  int& OutNumberOfElements = mrRemesh.OutMesh.GetNumberOfElements();

	  int* InElementList        = mrRemesh.InMesh.GetElementList();
	  double* InElementSizeList = mrRemesh.InMesh.GetElementSizeList();

	  int* OutElementList       = mrRemesh.OutMesh.GetElementList();


	  ModelPart::NodesContainerType::iterator nodes_begin = mrModelPart.NodesBegin();

	  for(int el = 0; el< OutNumberOfElements; el++)
	    {
	      if(mrRemesh.PreservedElements[el])
		{
		  Geometry<Node<3> > vertices;
		  for(unsigned int pn=0; pn<nds; pn++)
		    {

		      InElementList[id*nds+pn]= OutElementList[el*nds+pn];
		      vertices.push_back(*(nodes_begin + OutElementList[el*nds+pn]-1).base());
		    }

		  double element_size = 0;
		  mMesherUtilities.CalculateElementRadius(vertices, element_size);

		  InElementSizeList[id] = nodal_h_non_refining_factor * element_size;

		  id++;
		}

	    }

      }


      if( mEchoLevel > 0 ){
	std::cout<<"   Visited Elements: "<<id<<" [threshold:"<<refine_on_threshold<<"/size:"<<refine_on_size<<"]"<<std::endl;
	std::cout<<"   SELECT ELEMENTS TO REFINE ]; "<<std::endl;
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
        return "RefineElementsOnSizeMesherProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "RefineElementsOnSizeMesherProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }


    ///@}
    ///@name Friends
    ///@{

    ///@}


private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Static Member Variables
    ///@{
    ModelPart& mrModelPart;

    MesherUtilities::MeshingParameters& mrRemesh;

    MesherUtilities mMesherUtilities;

    int mEchoLevel;

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}



    /// Assignment operator.
    RefineElementsOnSizeMesherProcess& operator=(RefineElementsOnSizeMesherProcess const& rOther);


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
                                  RefineElementsOnSizeMesherProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const RefineElementsOnSizeMesherProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_REFINE_ELEMENTS_ON_SIZE_MEHSER_PROCESS_H_INCLUDED defined

