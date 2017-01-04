//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   
//

#if !defined(KRATOS_STRUCTURED_MESH_GENERATOR_PROCESS_H_INCLUDED )
#define  KRATOS_STRUCTURED_MESH_GENERATOR_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "processes/process.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"


namespace Kratos
{
  ///@addtogroup KratosCore
  ///@{

  ///@name Kratos Classes
  ///@{

  /// Short class definition.
  /** Detail class definition.
  */
  class StructuredMeshGeneratorProcess : public Process
    {
    public:
		using GeometryType = Geometry<Point<3> >;
      ///@name Type Definitions
      ///@{

      /// Pointer definition of StructuredMeshGeneratorProcess
      KRATOS_CLASS_POINTER_DEFINITION(StructuredMeshGeneratorProcess);

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor is deleted.
      StructuredMeshGeneratorProcess() = delete;

	  /// Constructor to be used. Takes the geometry to be meshed and ModelPart to be filled
	  StructuredMeshGeneratorProcess(GeometryType& rGeometry, ModelPart& rOutputModelPart, Parameters TheParameters);

	  /// The object is not copyable.
	  StructuredMeshGeneratorProcess(StructuredMeshGeneratorProcess const& rOther) = delete;

      /// Destructor.
      virtual ~StructuredMeshGeneratorProcess() ;

      ///@}
      ///@name Operators
      ///@{

	  /// It is not assignable.
	  StructuredMeshGeneratorProcess& operator=(StructuredMeshGeneratorProcess const& rOther) = delete;

      ///@}
      ///@name Operations
      ///@{

	  void Execute() override;


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
      virtual std::string Info() const;

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const;

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const;


      ///@}
      ///@name Friends
      ///@{


      ///@}

      private:
      ///@name Static Member Variables
      ///@{


      ///@}
      ///@name Member Variables
      ///@{
		  GeometryType& mrGeometry;

		  std::size_t mNumberOfDivisions;
		  std::size_t mStartNodeId;
		  std::size_t mStartElementId;
		  std::size_t mStartConditionId;
		  std::size_t mElementPropertiesId;
		  std::size_t mConditiongPropertiesId;
		  std::string mElementName;
		  std::string mConditionName;
		  bool mCrateSkinSubModelPart;
		  ModelPart& mrOutputModelPart;


      ///@}
      ///@name Private Operations
      ///@{

		  void Generate2DMesh();

		  void Generate3DMesh();

		  void GenerateNodes2D(Point<3> const& rMinPoint, Point<3> const& rMaxPoint);

		  void GenerateNodes3D(Point<3> const& rMinPoint, Point<3> const& rMaxPoint);

		  void GenerateTriangularElements();

		  void GenerateTetrahedraElements();

		  void CreateCellTetrahedra(std::size_t I, std::size_t J, std::size_t K, Properties::Pointer pProperties);

		  std::size_t GetNodeId(std::size_t I, std::size_t J, std::size_t K);

		  void GetLocalCoordinatesRange(Point<3>& rMinPoint, Point<3>& rMaxPoint);

      ///@}
      ///@name Private  Access
      ///@{


      ///@}
      ///@name Private Inquiry
      ///@{


      ///@}
      ///@name Un accessible methods
      ///@{



      ///@}

    }; // Class StructuredMeshGeneratorProcess

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    StructuredMeshGeneratorProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const StructuredMeshGeneratorProcess& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_STRUCTURED_MESH_GENERATOR_PROCESS_H_INCLUDED  defined
