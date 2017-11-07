//
//   Project Name:        KratosPfemApplication     $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:      February 2016 $
//   Revision:            $Revision:            0.0 $
//
//

#if !defined(KRATOS_PRINT_OUTPUT_MESH_PROCESS_H_INCLUDED )
#define  KRATOS_PRINT_OUTPUT_MESH_PROCESS_H_INCLUDED


// External includes
#include <iostream>
#include <fstream>

// System includes

// Project includes
#include "containers/variables_list_data_value_container.h"
#include "spatial_containers/spatial_containers.h"

#include "includes/model_part.h"
#include "custom_utilities/mesh_error_calculation_utilities.hpp"
#include "custom_utilities/modeler_utilities.hpp"


namespace Kratos
{

  ///@name Kratos Classes
  ///@{


  class PrintOutputMeshProcess
    : public Process
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION( PrintOutputMeshProcess );


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PrintOutputMeshProcess(ModelPart& rModelPart,
			   ModelerUtilities::MeshingParameters& rRemeshingParameters,
			   int EchoLevel) 
      : mrModelPart(rModelPart),
	mrRemesh(rRemeshingParameters)
    {
      mEchoLevel = EchoLevel;
    }


    /// Destructor.
    virtual ~PrintOutputMeshProcess() {}


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
    virtual void Execute()
    {
      KRATOS_TRY

      if( mEchoLevel > 0 ){
	std::cout<<" [ PRINT OUTPUT FROM MESHER: "<<std::endl;
	//std::cout<<"   Nodes before erasing : "<<mrModelPart.Nodes().size()<<std::endl;
      }

      PrintPointsXYZ();
      
      PrintMesh();
      
      std::cout<<"   PRINT OUTPUT FROM MESHER ]; "<<std::endl;	
      

      KRATOS_CATCH(" ")
    }


    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    virtual void ExecuteInitialize()
    {
    }

    /// this function is designed for being execute once before the solution loop but after all of the
    /// solvers where built
    virtual void ExecuteBeforeSolutionLoop()
    {
    }

    /// this function will be executed at every time step BEFORE performing the solve phase
    virtual void ExecuteInitializeSolutionStep()
    {	
    }

    /// this function will be executed at every time step AFTER performing the solve phase
    virtual void ExecuteFinalizeSolutionStep()
    {
    }

    /// this function will be executed at every time step BEFORE  writing the output
    virtual void ExecuteBeforeOutputStep()
    {
    }

    /// this function will be executed at every time step AFTER writing the output
    virtual void ExecuteAfterOutputStep()
    {
    }

    /// this function is designed for being called at the end of the computations
    /// right after reading the model and the groups
    virtual void ExecuteFinalize()
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
    virtual std::string Info() const
    {
      return "PrintOutputMeshProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
      rOStream << "PrintOutputMeshProcess";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
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
 
    ModelerUtilities::MeshingParameters& mrRemesh;

    int mEchoLevel;

    ///@}
    ///@name Un accessible methods
    ///@{


    void PrintPointsXYZ()
    {
      KRATOS_TRY

      const int& step = mrModelPart.GetProcessInfo()[STEP];
      
      const int& dimension = mrModelPart.GetProcessInfo()[SPACE_DIMENSION];
      
      std::string FileName;
      FileName += mrModelPart.Name();
      FileName += "_points_output_";
      FileName += std::to_string(step);
      FileName += ".txt";
      
      std::ofstream File;

      File.open(FileName);
	
      double* OutPointList   = mrRemesh.OutMesh.GetPointList();
      int& OutNumberOfPoints = mrRemesh.OutMesh.GetNumberOfPoints();

      
      unsigned int base = 0;
      for(unsigned int pn=0; pn<OutNumberOfPoints; pn++)
	{
	  std::string Point;
	  for(unsigned int i=0; i<dimension; i++)
	    {	      
	      Point += " ";
	      Point += std::to_string(OutPointList[base]);
	      base++;
	    }

	  Point += " \n";

	  File << Point;
	}

      File.close();
      
      KRATOS_CATCH(" ")
    }


    void PrintMesh()
    {
      KRATOS_TRY

      const int& step = mrModelPart.GetProcessInfo()[STEP];

      std::string FileName;
      FileName += mrModelPart.Name();
      FileName += "_mesh_output_";
      FileName += std::to_string(step);
      FileName += ".msh";
      
      std::ofstream File;

      File.open(FileName);
	
      const int& dimension = mrModelPart.GetProcessInfo()[SPACE_DIMENSION];
      
      if( dimension == 3 ) //number of nodes of a tetrahedron
	File <<  "mesh dimension 3 elemtype tetrahedra nnode 4 \n";
      else if( dimension == 2 ) //number of nodes of a triangle
	File <<  "mesh dimension 2 elemtype triangle nnode 3 \n";

      
      // write node coordinates
      PrintNodes(File);

      // write element connectivities
      PrintElements(File);

      File.close();
	
      KRATOS_CATCH(" ")
    }

    
    void PrintNodes(std::ofstream& File)
    {
      KRATOS_TRY
	
      File << "\n";
      File << "coordinates  \n";
      File << "\n";
      
      double* OutPointList   = mrRemesh.OutMesh.GetPointList();
      int& OutNumberOfPoints = mrRemesh.OutMesh.GetNumberOfPoints();

      const int& dimension = mrModelPart.GetProcessInfo()[SPACE_DIMENSION];
      
      unsigned int base = 0;
      for(unsigned int pn=0; pn<OutNumberOfPoints; pn++)
	{
	  std::string Point(std::to_string( pn+1 ));
	  
	  for(unsigned int i=0; i<dimension; i++)
	    {	      
	      Point += " ";
	      Point += std::to_string(OutPointList[base]);
	      base++;
	    }

	  Point += " \n";

	  File << Point;
	}


      File << "\n";
      File << "end coordinates  \n";
      File << "\n";
         
      KRATOS_CATCH(" ")
    }


    void PrintElements(std::ofstream& File)
    {
      KRATOS_TRY

      File << "\n";
      File << "elements  \n";
      File << "\n";
      
      int* OutElementList = mrRemesh.OutMesh.GetElementList();
      int& OutNumberOfElements = mrRemesh.OutMesh.GetNumberOfElements();

      const int& dimension = mrModelPart.GetProcessInfo()[SPACE_DIMENSION];

      int nds = 3; //number of nodes of a triangle
      if( dimension == 3 ) //number of nodes of a tetrahedron
	nds = 4;

      
      for(unsigned int el=0; el<OutNumberOfElements; el++)
	{
	  std::string Element(std::to_string( el+1 ));
	  
	  for(unsigned int pn=0; pn<nds; pn++)
	    {
	      Element += " ";
	      Element += std::to_string(OutElementList[el*nds+pn]);
	    }

	  Element += " \n";

	  File << Element;
	}


      File << "\n";
      File << "end elements  \n";
      File << "\n";
	
      KRATOS_CATCH(" ")
    }
  
    /// Assignment operator.
    PrintOutputMeshProcess& operator=(PrintOutputMeshProcess const& rOther);


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
				    PrintOutputMeshProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const PrintOutputMeshProcess& rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }
  ///@}


}  // namespace Kratos.

#endif // KRATOS_PRINT_OUTPUT_MESH_PROCESS_H_INCLUDED  defined 


