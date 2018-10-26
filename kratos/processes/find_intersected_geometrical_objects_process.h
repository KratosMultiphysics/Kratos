//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan dadvand
//

#if !defined(KRATOS_FIND_INTERSECTED_GEOMETRICAL_OBJECTS_PROCESS_H_INCLUDED )
#define  KRATOS_FIND_INTERSECTED_GEOMETRICAL_OBJECTS_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "processes/process.h"
#include "includes/model_part.h"
#include "spatial_containers/octree_binary.h"


namespace Kratos
{
	namespace Internals {
		class DistanceSpatialContainersConfigure
		{
		public:
			class CellNodeData
			{
				double mDistance;
				double mCoordinates[3];
				std::size_t mId;
			public:
				double& Distance() { return mDistance; }
				double& X() { return mCoordinates[0]; }
				double& Y() { return mCoordinates[1]; }
				double& Z() { return mCoordinates[2]; }
				double& Coordinate(int i) { return mCoordinates[i]; }
				std::size_t& Id() { return mId; }
			};




			///@name Type Definitions
			///@{

			enum {
				Dimension = 3,
				DIMENSION = 3,
				MAX_LEVEL = 12,
				MIN_LEVEL = 2    // this cannot be less than 2!!!
			};

			typedef Point			                                PointType;  /// always the point 3D
			typedef std::vector<double>::iterator                   DistanceIteratorType;
			typedef ModelPart::ElementsContainerType::ContainerType ContainerType;
			typedef ContainerType::value_type                       PointerType;
			typedef ContainerType::iterator                         IteratorType;
			typedef ModelPart::ElementsContainerType::ContainerType ResultContainerType;
			typedef ResultContainerType::value_type                 ResultPointerType;
			typedef ResultContainerType::iterator                   ResultIteratorType;

			typedef Element::Pointer                                        pointer_type;
			typedef CellNodeData                cell_node_data_type;
			typedef std::vector<CellNodeData*> data_type;

			typedef std::vector<PointerType>::iterator             PointerTypeIterator;




			/// Pointer definition of DistanceSpatialContainersConfigure
			KRATOS_CLASS_POINTER_DEFINITION(DistanceSpatialContainersConfigure);

			///@}
			///@name Life Cycle
			///@{

			/// Default constructor.
			DistanceSpatialContainersConfigure() {}

			/// Destructor.
			virtual ~DistanceSpatialContainersConfigure() {}


			///@}
			///@name Operators
			///@{


			///@}
			///@name Operations
			///@{

			static data_type* AllocateData() {
				return new data_type(27, (CellNodeData*)NULL);
			}

			static void CopyData(data_type* source, data_type* destination) {
				*destination = *source;
			}

			static void DeleteData(data_type* data) {
				delete data;
			}

			static inline void CalculateBoundingBox(const PointerType& rObject, PointType& rLowPoint, PointType& rHighPoint)
			{
				rHighPoint = rObject->GetGeometry().GetPoint(0);
				rLowPoint = rObject->GetGeometry().GetPoint(0);

				for (unsigned int point = 0; point<rObject->GetGeometry().PointsNumber(); point++)
				{
					for (std::size_t i = 0; i<3; i++)
					{
						rLowPoint[i] = (rLowPoint[i]  >  rObject->GetGeometry().GetPoint(point)[i]) ? rObject->GetGeometry().GetPoint(point)[i] : rLowPoint[i];
						rHighPoint[i] = (rHighPoint[i] <  rObject->GetGeometry().GetPoint(point)[i]) ? rObject->GetGeometry().GetPoint(point)[i] : rHighPoint[i];
					}
				}
			}

			static inline void GetBoundingBox(const PointerType rObject, double* rLowPoint, double* rHighPoint)
			{

				for (std::size_t i = 0; i<3; i++)
				{
					rLowPoint[i] = rObject->GetGeometry().GetPoint(0)[i];
					rHighPoint[i] = rObject->GetGeometry().GetPoint(0)[i];
				}

				for (unsigned int point = 0; point<rObject->GetGeometry().PointsNumber(); point++)
				{
					for (std::size_t i = 0; i<3; i++)
					{
						rLowPoint[i] = (rLowPoint[i]  >  rObject->GetGeometry().GetPoint(point)[i]) ? rObject->GetGeometry().GetPoint(point)[i] : rLowPoint[i];
						rHighPoint[i] = (rHighPoint[i] <  rObject->GetGeometry().GetPoint(point)[i]) ? rObject->GetGeometry().GetPoint(point)[i] : rHighPoint[i];
					}
				}
			}

			static inline bool Intersection(const PointerType& rObj_1, const PointerType& rObj_2)
			{
				Element::GeometryType& geom_1 = rObj_1->GetGeometry();
				Element::GeometryType& geom_2 = rObj_2->GetGeometry();
				return  geom_1.HasIntersection(geom_2);

			}


			static inline bool  IntersectionBox(const PointerType& rObject, const PointType& rLowPoint, const PointType& rHighPoint)
			{
				return rObject->GetGeometry().HasIntersection(rLowPoint, rHighPoint);
			}


			static  inline bool  IsIntersected(const Element::Pointer rObject, double Tolerance, const double* rLowPoint, const double* rHighPoint)
			{
				Point low_point(rLowPoint[0] - Tolerance, rLowPoint[1] - Tolerance, rLowPoint[2] - Tolerance);
				Point high_point(rHighPoint[0] + Tolerance, rHighPoint[1] + Tolerance, rHighPoint[2] + Tolerance);

				KRATOS_THROW_ERROR(std::logic_error, "Not Implemented method", "")
					//return HasIntersection(rObject->GetGeometry(), low_point, high_point);
			}


			///@}
			///@name Input and output
			///@{

			/// Turn back information as a string.
			virtual std::string Info() const
			{
				return " Spatial Containers Configure";
			}

			/// Print information about this object.
			virtual void PrintInfo(std::ostream& rOStream) const {}

			/// Print object's data.
			virtual void PrintData(std::ostream& rOStream) const {}


			///@}

		protected:

		private:

			/// Assignment operator.
			DistanceSpatialContainersConfigure& operator=(DistanceSpatialContainersConfigure const& rOther);

			/// Copy constructor.
			DistanceSpatialContainersConfigure(DistanceSpatialContainersConfigure const& rOther);


		}; // Class DistanceSpatialContainersConfigure

	} // manespace Internals

  ///@addtogroup KratosCore
  ///@{

  ///@name Kratos Classes
  ///@{

  /// This class takes two modelparts and marks the intersected ones with SELECTED flag.
  /** It creates a spatial datastructure and search for interaction.
  It also provides some helper methods for derived classes to check
  individual element or condition interesections.
  */
  class KRATOS_API(KRATOS_CORE) FindIntersectedGeometricalObjectsProcess : public Process
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of FindIntersectedGeometricalObjectsProcess
      KRATOS_CLASS_POINTER_DEFINITION(FindIntersectedGeometricalObjectsProcess);

	  using ConfigurationType = Internals::DistanceSpatialContainersConfigure;
	  using CellType = OctreeBinaryCell<ConfigurationType>;
	  using OctreeType = OctreeBinary<CellType>;
	  using CellNodeDataType = ConfigurationType::cell_node_data_type;

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor.
	  FindIntersectedGeometricalObjectsProcess() = delete;

	  /// Copy constructor.
	  FindIntersectedGeometricalObjectsProcess(FindIntersectedGeometricalObjectsProcess const& rOther) = delete;

	  /// Constructor to be used.
	  FindIntersectedGeometricalObjectsProcess(ModelPart& rPart1, ModelPart& rPart2);


	  /// Destructor.
	  ~FindIntersectedGeometricalObjectsProcess() override {}

	  ///@name Member Variables
	  ///@{

	  std::vector<PointerVector<GeometricalObject>> mIntersectedObjects;

      ///@}
      ///@name Operations
      ///@{

	  virtual void Initialize();

	  virtual void FindIntersectedSkinObjects(std::vector<PointerVector<GeometricalObject>>& rResults);

	  virtual void FindIntersections();

	  virtual std::vector<PointerVector<GeometricalObject>>& GetIntersections();

	  virtual ModelPart& GetModelPart1();

	  virtual OctreeBinary<OctreeBinaryCell<Internals::DistanceSpatialContainersConfigure>>* GetOctreePointer();

	  virtual void Clear();

	  void Execute() override;

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

    private:
      ///@name Static Member Variables
      ///@{


      ///@}
      ///@name Member Variables
      ///@{

		ModelPart& mrModelPart1;
		ModelPart& mrModelPart2;
		OctreeType mOctree;

      ///@}
      ///@name Private Operations
      ///@{

		void GenerateOctree();
		void SetOctreeBoundingBox();
		void MarkIfIntersected(Element& rElement1, std::vector<OctreeType::cell_type*>& leaves);
		bool HasIntersection(Element::GeometryType& rFirstGeometry, Element::GeometryType& rSecondGeometry);
		bool HasIntersection2D(Element::GeometryType& rFirstGeometry, Element::GeometryType& rSecondGeometry);
		bool HasIntersection3D(Element::GeometryType& rFirstGeometry, Element::GeometryType& rSecondGeometry);
		void FindIntersectedSkinObjects(Element& rElement1, std::vector<OctreeType::cell_type*>& leaves, PointerVector<GeometricalObject>& rResults);


      ///@}
      ///@name Un accessible methods
      ///@{

      /// Assignment operator.
      FindIntersectedGeometricalObjectsProcess& operator=(FindIntersectedGeometricalObjectsProcess const& rOther);


      ///@}

    }; // Class FindIntersectedGeometricalObjectsProcess

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    FindIntersectedGeometricalObjectsProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const FindIntersectedGeometricalObjectsProcess& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_FIND_INTERSECTED_GEOMETRICAL_OBJECTS_PROCESS_H_INCLUDED  defined
