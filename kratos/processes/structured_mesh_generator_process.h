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
  class KRATOS_API(KRATOS_CORE) StructuredMeshGeneratorProcess : public Process
    {
    public:
        using GeometryType = Geometry<Node<3> >;
      ///@name Type Definitions
      ///@{

      /// Pointer definition of StructuredMeshGeneratorProcess
      KRATOS_CLASS_POINTER_DEFINITION(StructuredMeshGeneratorProcess);

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor is deleted.
      StructuredMeshGeneratorProcess() = delete;

      /// Constructors to be used. They take the geometry to be meshed and ModelPart to be filled. The second constructor is
      /// provided for the Python interface.
      StructuredMeshGeneratorProcess(const GeometryType& rGeometry, ModelPart& rOutputModelPart, Parameters& TheParameters);

      StructuredMeshGeneratorProcess(GeometryType::Pointer pGeometry, ModelPart& rOutputModelPart, Parameters& TheParameters):
          StructuredMeshGeneratorProcess(*pGeometry, rOutputModelPart, TheParameters){KRATOS_WATCH(mNumberOfDivisions)};

      /// The object is not copyable.
	  StructuredMeshGeneratorProcess(StructuredMeshGeneratorProcess const& rOther) = delete;

      /// Destructor.
      ~StructuredMeshGeneratorProcess() override ;

      ///@}
      ///@name Operators
      ///@{

	  /// It is not assignable.
	  StructuredMeshGeneratorProcess& operator=(StructuredMeshGeneratorProcess const& rOther) = delete;

      ///@}
      ///@name Operations
      ///@{

	  void Execute() override;

      int Check() override;

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
      std::string Info() const override;

      /// Print information about this object.
      void PrintInfo(std::ostream& rOStream) const override;

      /// Print object's data.
      void PrintData(std::ostream& rOStream) const override;


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
          const GeometryType& mrGeometry;
		  std::size_t mNumberOfDivisions;
		  std::size_t mStartNodeId;
		  std::size_t mStartElementId;
		  std::size_t mStartConditionId;
		  std::size_t mElementPropertiesId;
		  std::size_t mConditiongPropertiesId;
		  std::string mElementName;
		  std::string mConditionName;
          bool mCreateSkinSubModelPart;
		  ModelPart& mrOutputModelPart;


      ///@}
      ///@name Private Operations
      ///@{

		  void Generate2DMesh();

          void Generate3DMesh();

		  void GenerateNodes2D(Point const& rMinPoint, Point const& rMaxPoint);

		  void GenerateNodes3D(Point const& rMinPoint, Point const& rMaxPoint);

		  void GenerateTriangularElements();

		  void GenerateTetrahedraElements();

		  void CreateCellTetrahedra(std::size_t I, std::size_t J, std::size_t K, Properties::Pointer pProperties);

		  std::size_t GetNodeId(std::size_t I, std::size_t J, std::size_t K);

		  void GetLocalCoordinatesRange(Point& rMinPoint, Point& rMaxPoint);

          bool CheckDomainGeometry();

          bool CheckDomainGeometryConnectivityForQuadrilateral2D4();

          bool CheckDomainGeometryConnectivityForHexahedra3D8();

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
