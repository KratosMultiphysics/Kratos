//
//   Project Name:        KratosDelaunayMeshingApplication $
//   Created by:          $Author:             JMCarbonell $
//   Last modified by:    $Co-Author:                      $
//   Date:                $Date:                April 2018 $
//   Revision:            $Revision:                   0.0 $
//
//

#if !defined(KRATOS_PRINT_MESH_OUTPUT_MESHER_PROCESS_H_INCLUDED )
#define  KRATOS_PRINT_MESH_OUTPUT_MESHER_PROCESS_H_INCLUDED


// External includes
#include <fstream>

// System includes

// Project includes
#include "containers/variables_list_data_value_container.h"
#include "spatial_containers/spatial_containers.h"

#include "includes/model_part.h"
#include "custom_utilities/mesh_error_calculation_utilities.hpp"
#include "custom_utilities/mesher_utilities.hpp"
#include "custom_processes/mesher_process.hpp"

namespace Kratos
{

  ///@name Kratos Classes
  ///@{


  class PrintMeshOutputMesherProcess
    : public MesherProcess
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION( PrintMeshOutputMesherProcess );


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PrintMeshOutputMesherProcess(ModelPart& rModelPart,
			   MesherUtilities::MeshingParameters& rRemeshingParameters,
			   std::string FileName,
			   int EchoLevel)
      : mrModelPart(rModelPart),
	mrRemesh(rRemeshingParameters)
    {
      mFileName  = FileName;
      mEchoLevel = EchoLevel;
    }


    /// Destructor.
    virtual ~PrintMeshOutputMesherProcess() {}


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
	 std::cout<<" [ PRINT IN/OUT MESHER: ("<<mFileName<<") "<<std::endl;
	//std::cout<<"   Nodes before erasing : "<<mrModelPart.Nodes().size()<<std::endl;
      }

      PrintPointsXYZ();

      PrintMesh();

      std::cout<<"   PRINT IN/OUT MESHER ]; "<<std::endl;


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
      return "PrintMeshOutputMesherProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << "PrintMeshOutputMesherProcess";
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

    std::string mFileName;

    int mEchoLevel;

    ///@}
    ///@name Un accessible methods
    ///@{


    void PrintPointsXYZ()
    {
      KRATOS_TRY

      const int& step = mrModelPart.GetProcessInfo()[STEP];

      const unsigned int& dimension = mrModelPart.GetProcessInfo()[SPACE_DIMENSION];

      std::string FileName;
      FileName += mrModelPart.Name();
      FileName += "_points_";
      FileName += mFileName;
      FileName += "_";
      FileName += std::to_string(step);
      FileName += ".txt";

      std::ofstream File;

      File.open(FileName);

      double* PointList;
      unsigned int NumberOfPoints;
      if( mFileName == "input" ){
	  PointList = mrRemesh.InMesh.GetPointList();
	  NumberOfPoints = mrRemesh.InMesh.GetNumberOfPoints();
      }
      else{
	  PointList = mrRemesh.OutMesh.GetPointList();
	  NumberOfPoints = mrRemesh.OutMesh.GetNumberOfPoints();
      }

      unsigned int base = 0;
      for(unsigned int pn=0; pn<NumberOfPoints; pn++)
	{
	  std::string Point;
	  for(unsigned int i=0; i<dimension; i++)
	    {
	      Point += " ";
	      Point += std::to_string(PointList[base]);
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
      FileName += "_mesh_";
      FileName += mFileName;
      FileName += "_";
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

      double* PointList;
      unsigned int NumberOfPoints;
      if( mFileName == "input" ){
	  PointList = mrRemesh.InMesh.GetPointList();
	  NumberOfPoints = mrRemesh.InMesh.GetNumberOfPoints();
      }
      else{
	  PointList = mrRemesh.OutMesh.GetPointList();
	  NumberOfPoints = mrRemesh.OutMesh.GetNumberOfPoints();
      }

      const unsigned int& dimension = mrModelPart.GetProcessInfo()[SPACE_DIMENSION];

      unsigned int base = 0;
      for(unsigned int pn=0; pn<NumberOfPoints; pn++)
	{
	  std::string Point(std::to_string( pn+1 ));

	  for(unsigned int i=0; i<dimension; i++)
	    {
	      Point += " ";
	      Point += std::to_string(PointList[base]);
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

      int* ElementList;
      unsigned int NumberOfElements;
      if( mFileName == "input" ){
	  ElementList = mrRemesh.InMesh.GetElementList();
	  NumberOfElements = mrRemesh.InMesh.GetNumberOfElements();
      }
      else{
	  ElementList = mrRemesh.OutMesh.GetElementList();
	  NumberOfElements = mrRemesh.OutMesh.GetNumberOfElements();
      }


      const unsigned int& dimension = mrModelPart.GetProcessInfo()[SPACE_DIMENSION];

      unsigned int nds = 3; //number of nodes of a triangle
      if( dimension == 3 ) //number of nodes of a tetrahedron
	nds = 4;


      for(unsigned int el=0; el<NumberOfElements; el++)
	{
	  std::string Element(std::to_string( el+1 ));

	  for(unsigned int pn=0; pn<nds; pn++)
	    {
	      Element += " ";
	      Element += std::to_string(ElementList[el*nds+pn]);
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
    PrintMeshOutputMesherProcess& operator=(PrintMeshOutputMesherProcess const& rOther);


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
				    PrintMeshOutputMesherProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const PrintMeshOutputMesherProcess& rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }
  ///@}


}  // namespace Kratos.

#endif // KRATOS_PRINT_OUTPUT_MESH_PROCESS_H_INCLUDED  defined


