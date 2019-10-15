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

#if !defined(KRATOS_VOXEL_MESH_COLORING_PROCESS_H_INCLUDED )
#define  KRATOS_VOXEL_MESH_COLORING_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "processes/process.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "geometries/hexahedra_3d_8.h"
#include "processes/find_intersected_geometrical_objects_process.h"
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
  class KRATOS_API(KRATOS_CORE) VoxelMeshColoringProcess : public Process
    {
    public:
        using GeometryType = Geometry<Node<3> >;
      ///@name Type Definitions
      ///@{

      /// Pointer definition of VoxelMeshColoringProcess
      KRATOS_CLASS_POINTER_DEFINITION(VoxelMeshColoringProcess);
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
      VoxelMeshColoringProcess() = delete;

      /// Constructors to be used. They take the geometry to be meshed and ModelPart to be filled. The second constructor is
      /// provided for the Python interface.
      VoxelMeshColoringProcess(Internals::CartesianMeshColors& MeshColors,  array_1d<std::size_t,3> const& NumberOfDivisions, 
        ModelPart& rVolumePart,
        ModelPart& rSkinPart, Parameters& TheParameters);

      /// The object is not copyable.
	  VoxelMeshColoringProcess(VoxelMeshColoringProcess const& rOther) = delete;

      /// Destructor.
      ~VoxelMeshColoringProcess() override ;

      ///@}
      ///@name Operators
      ///@{

	  /// It is not assignable.
	  VoxelMeshColoringProcess& operator=(VoxelMeshColoringProcess const& rOther) = delete;

      ///@}
      ///@name Operations
      ///@{

	  void Initialize();

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
        Internals::CartesianMeshColors& mColors;
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
          FindIntersectedGeometricalObjectsProcess mFindIntersectedObjectsProcess;
          std::vector<bool> mCellIsEmpty;
          array_1d<double,3> mCellSizes;

          const double mExtraRaysEpsilon = 1.0e-8;

          double mInsideColor;
          double mOutsideColor;
          bool mApplyOutsideColor;
          std::string mColoringEntities;


      ///@}
      ///@name Private Operations
      ///@{

          void CalculateVoxelsColor(ModelPart& TheSubModelPart, int TheColor);

          void InitializeRays(array_1d< std::size_t, 3 > const& MinCellPosition, array_1d< std::size_t, 3 > const& MaxCellPosition);

          void CalculateRayDistances(ModelPart& TheSubModelPart, int TheColor);

          void GenerateTriangularElements();

          void GenerateTetrahedraElements();

          void CreateCellTetrahedra(std::size_t I, std::size_t J, std::size_t K, Properties::Pointer pProperties);

          std::size_t GetNodeId(std::size_t I, std::size_t J, std::size_t K);

          double DistancePositionInSpace(const Node<3> &rNode);
	
          void GetRayIntersections(const double* ray, const unsigned int direction,
		        std::vector<std::pair<double,Element::GeometryType*> >& rIntersections);

          int GetCellIntersections(
                OctreeType::cell_type* cell,
                const double* ray,
                OctreeType::key_type* ray_key,
                const unsigned int direction,
                std::vector<std::pair<double, Element::GeometryType*> > &rIntersections);

          int ComputeRayIntersection(
                Element::GeometryType& rGeometry,
                const double* pRayPoint1,
                const double* pRayPoint2,
                double* pIntersectionPoint);

          void ComputeExtraRayColors(
                  const double Epsilon,
                  const double RayPerturbation,
                  const array_1d<double,3> &rCoords,
                  array_1d<double,3> &rDistances);

          void GetExtraRayOrigins(
                  const double RayEpsilon,
                  const array_1d<double,3> &rCoords,
                  std::vector<array_1d<double,3>> &rExtraRayOrigs);

          void CorrectExtraRayOrigin(double* ExtraRayCoords);

          void MarkIntersectedCells(Element::GeometryType&);

          template<typename TPointsContainerType>
          void CalculateMinMaxCellsPositions(TPointsContainerType const& Points, array_1d< std::size_t, 3 >& MinCellPosition, array_1d< std::size_t, 3 >& MaxCellPosition);

          std::size_t CalculateCellIndex(Point const &ThePoint ) const;

          std::size_t CalculatePosition( double Coordinate, int ThisDimension ) const;

          bool CellIntersectGeometry(std::size_t Ix, std::size_t Iy, std::size_t Iz, Element::GeometryType& TheGeometry);

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

    }; // Class VoxelMeshColoringProcess

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    VoxelMeshColoringProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const VoxelMeshColoringProcess& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_VOXEL_MESH_COLORING_PROCESS_H_INCLUDED  defined
