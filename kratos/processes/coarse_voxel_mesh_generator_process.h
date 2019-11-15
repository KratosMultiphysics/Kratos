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

#if !defined(KRATOS_COARSE_VOXEL_MESH_GENERATOR_PROCESS_H_INCLUDED )
#define  KRATOS_COARSE_VOXEL_MESH_GENERATOR_PROCESS_H_INCLUDED



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
#include "processes/find_intersected_geometrical_objects_process.h"
#include "processes/internals/cartesian_mesh_colors.h"
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
  class KRATOS_API(KRATOS_CORE) CoarseVoxelMeshGeneratorProcess : public VoxelMeshGeneratorProcess
    {
    public:
        using GeometryType = Geometry<Node<3> >;
      ///@name Type Definitions
      ///@{

      /// Pointer definition of CoarseVoxelMeshGeneratorProcess
      KRATOS_CLASS_POINTER_DEFINITION(CoarseVoxelMeshGeneratorProcess);

      typedef Element::GeometryType IntersectionGeometryType;
      typedef std::vector<std::pair<double, IntersectionGeometryType*> > IntersectionsContainerType;

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor is deleted.
      CoarseVoxelMeshGeneratorProcess() = delete;

      CoarseVoxelMeshGeneratorProcess(Point const& MinPoint, Point const& MaxPoint,
        ModelPart& rVolumePart,
        ModelPart& rSkinPart, Parameters& TheParameters);

      CoarseVoxelMeshGeneratorProcess(std::vector<double> const& XCoordinates, std::vector<double> const& YCoordinates, std::vector<double> const& ZCoordinates,
        ModelPart& rVolumePart,
        ModelPart& rSkinPart, Parameters& TheParameters);

      /// The object is not copyable.
	  CoarseVoxelMeshGeneratorProcess(CoarseVoxelMeshGeneratorProcess const& rOther) = delete;

      /// Destructor.
      ~CoarseVoxelMeshGeneratorProcess() override ;

      ///@}
      ///@name Operators
      ///@{

	  /// It is not assignable.
	  CoarseVoxelMeshGeneratorProcess& operator=(CoarseVoxelMeshGeneratorProcess const& rOther) = delete;

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
      ///@name Private Operations
      ///@{

        void Generate3DCoarseMesh();

      ///@}

    }; // Class CoarseVoxelMeshGeneratorProcess

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    CoarseVoxelMeshGeneratorProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const CoarseVoxelMeshGeneratorProcess& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_COARSE_VOXEL_MESH_GENERATOR_PROCESS_H_INCLUDED  defined
