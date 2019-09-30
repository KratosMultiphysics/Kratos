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


// External includes


// Project includes
#include "processes/process.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "geometries/hexahedra_3d_8.h"
#include "processes/find_intersected_geometrical_objects_process.h"


namespace Kratos
{
  namespace Internals
  {
    class CartesianRay{
      int mDirection;
      Point mPoint1;
      Point mPoint2;
      std::vector<double> mIntersectionPoints;
    public:
      CartesianRay() = delete;

      CartesianRay(int Direction, Point const& Point1, Point const& Point2): mDirection(Direction), mPoint1(Point1), mPoint2(Point2) {}

      template<typename TGeometryType>
      void AddIntersection(TGeometryType const& rGeometry, double Tolerance){
        // Call the line - triangle intersection util
        array_1d<double,3> intersection_point = ZeroVector(3);
        const int is_intersected = IntersectionUtilities::ComputeTriangleLineIntersection(
          rGeometry,
          mPoint1,
          mPoint2,
          intersection_point,
          Tolerance);

          if(is_intersected == 1){ // There is an intersection but not coplanar
            mIntersectionPoints.push_back(intersection_point[mDirection]);
          }
      }

      void CollapseIntersectionPoints(double Tolerance){
        if (!mIntersectionPoints.empty()) {
            // Sort
            std::sort(mIntersectionPoints.begin(), mIntersectionPoints.end());
            // Unique
            auto i_begin = mIntersectionPoints.begin();
            auto i_intersection = mIntersectionPoints.begin();
            while (++i_begin != mIntersectionPoints.end()) {
                // considering the very near points as the same points
                if (std::abs(*i_begin - *i_intersection) > Tolerance) // if the hit points are far enough they are not the same
                    *(++i_intersection) = *i_begin;
            }
          mIntersectionPoints.resize((++i_intersection) - mIntersectionPoints.begin());
        }
      }

      std::vector<double> const& GetIntersectionPoints() const {return mIntersectionPoints;}
    };

    class CartesianMeshColors{
      array_1d<DenseVector<double>, 3> mCoordinates;
      DenseVector<array_1d<double,3>> mColors;
     public:
      CartesianMeshColors(){}
      DenseVector<double> const& GetCoordinates(int Index) const {return mCoordinates[Index];}
      void SetCoordinates(array_1d<DenseVector<double>, 3>&& TheCoordinates){
        mCoordinates = TheCoordinates;
        mColors.resize(mCoordinates[0].size()*mCoordinates[1].size()*mCoordinates[2].size());
      }

      array_1d<double, 3>& GetColor(std::size_t I, std::size_t J, std::size_t K){
        const std::size_t index = I + J * mCoordinates[0].size() + K * mCoordinates[1].size() * mCoordinates[0].size();
        return mColors[index];
      }

      Point GetPoint(std::size_t I, std::size_t J, std::size_t K){
        return Point(mCoordinates[0][I], mCoordinates[1][J], mCoordinates[2][K]);
      }

    };
  }

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
