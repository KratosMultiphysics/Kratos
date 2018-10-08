//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		     BSD License
//					         Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//


#if !defined(KRATOS_CALCULATE_DISTANCE_TO_SKIN_PROCESS_H_INCLUDED )
#define  KRATOS_CALCULATE_DISTANCE_TO_SKIN_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "processes/find_intersected_geometrical_objects_process.h"
#include "processes/calculate_discontinuous_distance_to_skin_process.h"

namespace Kratos
{
  ///@addtogroup Kratos Core
  ///@{

  ///@name Kratos Classes
  ///@{

  /// Calculates the nodal distances using elemental discontinuous distances.
  /** This class calculates the nodal distances as a minimum elemental distances connected to it.
  */
  template<std::size_t TDim = 3>
  class KRATOS_API(KRATOS_CORE) CalculateDistanceToSkinProcess : public CalculateDiscontinuousDistanceToSkinProcess<TDim>
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of CalculateDistanceToSkinProcess
      KRATOS_CLASS_POINTER_DEFINITION(CalculateDistanceToSkinProcess);

      //TODO: These using statements have been included to make the old functions able to compile. It is still pending to update them.
      using ConfigurationType = Internals::DistanceSpatialContainersConfigure;
	  using CellType = OctreeBinaryCell<ConfigurationType>;
	  using OctreeType = OctreeBinary<CellType>;
	  using CellNodeDataType = ConfigurationType::cell_node_data_type;

      ///@}
      ///@name Life Cycle
      ///@{

	  /// Constructor to be used.
	  CalculateDistanceToSkinProcess(ModelPart& rVolumePart, ModelPart& rSkinPart);

	  /// Destructor.
      ~CalculateDistanceToSkinProcess() override;

	  ///@}
	  ///@name Deleted
	  ///@{

      /// Default constructor.
      CalculateDistanceToSkinProcess() = delete;;

	  /// Copy constructor.
	  CalculateDistanceToSkinProcess(CalculateDistanceToSkinProcess const& rOther) = delete;

	  /// Assignment operator.
	  CalculateDistanceToSkinProcess& operator=(CalculateDistanceToSkinProcess const& rOther) = delete;

	  ///@}
      ///@name Operations
      ///@{
      void Initialize() override;

      void CalculateDistances(std::vector<PointerVector<GeometricalObject>>& rIntersectedObjects) override;

      double CalculateDistanceToNode(Element &rElement1, const int NodeIndex, PointerVector<GeometricalObject> &rIntersectedObjects, const double Epsilon);

      void CalculateElementalDistances(std::vector<PointerVector<GeometricalObject>> &rIntersectedObjects);

      virtual void InitializeNodalDistances();

      virtual void CalculateNodalDistances();

      virtual void CalculateNodesDistances(); //TODO: This method has been adapted from the previous implementation. It is still pending to update it.

      virtual void CalculateNodeDistance(Node<3>& rNode); //TODO: This method has been adapted from the previous implementation. It is still pending to update it.

      virtual double DistancePositionInSpace(double* pCoords); //TODO: This method has been adapted from the previous implementation. It is still pending to update it.

      virtual void GetRayIntersections(double* ray, int direction, std::vector<std::pair<double,Element::GeometryType*> >& intersections); //TODO: This method has been adapted from the previous implementation. It is still pending to update it.

      virtual int GetCellIntersections(OctreeType::cell_type* cell, double* ray,
                                       OctreeType::key_type* ray_key, int direction,
                                       std::vector<std::pair<double, Element::GeometryType*> >& intersections); //TODO: This method has been adapted from the previous implementation. It is still pending to update it.


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


      ///@}
      ///@name Private Operators
      ///@{


      ///@}
      ///@name Private Operations
      ///@{

      double inline CalculatePointDistance(
        const Element::GeometryType &rIntObjGeom,
        const Point &rDistancePoint);

      int ComputeRayIntersection(
        Element::GeometryType& rGeometry,
        const double* pRayPoint1,
        const double* pRayPoint2,
        double* pIntersectionPoint);

      //TODO: This method has been adapted from the previous implementation. It is still pending to update it.
      int IntersectionTriangleSegment(
        Element::GeometryType& rGeometry,
        const double* RayPoint1,
        const double* RayPoint2,
        double* IntersectionPoint);

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

    }; // Class CalculateDistanceToSkinProcess

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    CalculateDistanceToSkinProcess<>& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const CalculateDistanceToSkinProcess<>& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CALCULATE_DISTANCE_TO_SKIN_PROCESS_H_INCLUDED  defined
