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

#if !defined(KRATOS_VOXEL_MESH_GENERATOR_PROCESS_H_INCLUDED )
#define  KRATOS_VOXEL_MESH_GENERATOR_PROCESS_H_INCLUDED



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


namespace Kratos
{

  ///@addtogroup KratosCore
  ///@{

  ///@name Kratos Classes
  ///@{

  /// Short class definition.
  /** Detail class definition.
  */
  class KRATOS_API(KRATOS_CORE) VoxelMeshGeneratorProcess : public Process
    {
    public:
        using GeometryType = Geometry<Node<3> >;
      ///@name Type Definitions
      ///@{

      /// Pointer definition of VoxelMeshGeneratorProcess
      KRATOS_CLASS_POINTER_DEFINITION(VoxelMeshGeneratorProcess);
      using ConfigurationType = Internals::DistanceSpatialContainersConfigure;
      using CellType = OctreeBinaryCell<ConfigurationType>;
      using OctreeType = OctreeBinary<CellType>;
      using CellNodeDataType = ConfigurationType::cell_node_data_type;

      typedef Element::GeometryType IntersectionGeometryType;
      typedef std::vector<std::pair<double, IntersectionGeometryType*> > IntersectionsContainerType;

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor is deleted.
      VoxelMeshGeneratorProcess() = delete;

      /// Constructors to be used. They take the geometry to be meshed and ModelPart to be filled. The second constructor is
      /// provided for the Python interface.
      VoxelMeshGeneratorProcess(Point const& MinPoint, Point const& MaxPoint,
        ModelPart& rVolumePart,
        ModelPart& rSkinPart, Parameters& TheParameters);

      /// The object is not copyable.
	  VoxelMeshGeneratorProcess(VoxelMeshGeneratorProcess const& rOther) = delete;

      /// Destructor.
      ~VoxelMeshGeneratorProcess() override ;

      ///@}
      ///@name Operators
      ///@{

	  /// It is not assignable.
	  VoxelMeshGeneratorProcess& operator=(VoxelMeshGeneratorProcess const& rOther) = delete;

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
        Internals::CartesianMeshColors mColors;
          const Point mMinPoint;
          const Point mMaxPoint;
		  array_1d<std::size_t,3> mNumberOfDivisions;
		  std::size_t mStartNodeId;
		  std::size_t mStartElementId;
		  std::size_t mStartConditionId;
		  std::size_t mElementPropertiesId;
		  std::size_t mConditiongPropertiesId;
		  std::string mElementName;
		  std::string mConditionName;
          bool mCreateSkinSubModelPart;
		  ModelPart& mrVolumePart;
		  ModelPart& mrSkinPart;
          std::vector<bool> mCellIsEmpty;
          array_1d<double,3> mCellSizes;

          const double mExtraRaysEpsilon = 1.0e-8;
          Parameters mColoringParameters;
          std::string mEntitiesToGenerate;
          bool mCoarseMeshType;


      ///@}
      ///@name Private Operations
      ///@{

          void Generate3DMesh();

          void GenerateNodes3D(Point const& rMinPoint, Point const& rMaxPoint);

          void GenerateCenterOfElements(Point const& rMinPoint, Point const& rMaxPoint);

          Node<3>::Pointer pGetNode(std::size_t I, std::size_t J, std::size_t K);
          Node<3>& GetCenterNode(std::size_t I, std::size_t J, std::size_t K);
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

    }; // Class VoxelMeshGeneratorProcess

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    VoxelMeshGeneratorProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const VoxelMeshGeneratorProcess& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_VOXEL_MESH_GENERATOR_PROCESS_H_INCLUDED  defined
