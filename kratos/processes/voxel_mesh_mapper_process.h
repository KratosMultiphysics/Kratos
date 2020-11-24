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

#if !defined(KRATOS_VOXEL_MESH_MAPPER_PROCESS_H_INCLUDED )
#define  KRATOS_VOXEL_MESH_MAPPER_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <algorithm>


// External includes


// Project includes
#include "processes/process.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "geometries/hexahedra_3d_8.h"
#include "processes/voxel_mesh_generator_process.h"


namespace Kratos
{

  ///@addtogroup KratosCore
  ///@{

  ///@name Kratos Classes
  ///@{

  /// Short class definition.
  /** Detail class definition.
  */
  class KRATOS_API(KRATOS_CORE) VoxelMeshMapperProcess : public VoxelMeshGeneratorProcess
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of VoxelMeshMapperProcess
      KRATOS_CLASS_POINTER_DEFINITION(VoxelMeshMapperProcess);


      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor is deleted.
      VoxelMeshMapperProcess() = delete;

      /// Constructors to be used. They take the geometry to be meshed and Input file containing the results to be Mapped. The second constructor is
      /// provided for the Python interface.
      VoxelMeshMapperProcess(std::string InputFileName, ModelPart& rVolumePart,
        ModelPart& rSkinPart, Parameters& TheParameters);

      /// The object is not copyable.
	  VoxelMeshMapperProcess(VoxelMeshMapperProcess const& rOther) = delete;

      /// Destructor.
      ~VoxelMeshMapperProcess() override ;

      ///@}
      ///@name Operators
      ///@{

	  /// It is not assignable.
	  VoxelMeshMapperProcess& operator=(VoxelMeshMapperProcess const& rOther) = delete;

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

      protected:
      ///@name Static Member Variables
      ///@{


      ///@}
      ///@name Member Variables
      ///@{

        Internals::CartesianMeshColors mInputMesh;

      ///@}
      ///@name Private Operations
      ///@{
        
        void ReadInputFile(std::string InputFileName);
        void CenterToNodalCoordinates(std::vector<double>& rCoordinates);

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

    }; // Class VoxelMeshMapperProcess

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    VoxelMeshMapperProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const VoxelMeshMapperProcess& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_VOXEL_MESH_MAPPER_PROCESS_H_INCLUDED  defined
