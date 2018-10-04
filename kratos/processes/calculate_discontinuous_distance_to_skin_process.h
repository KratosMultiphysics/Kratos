//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		     BSD License
//					         Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand, Ruben Zorrilla
//

#if !defined(KRATOS_CALCULATE_DISCONTINUOUS_DISTANCE_TO_SKIN_PROCESS_H_INCLUDED )
#define  KRATOS_CALCULATE_DISCONTINUOUS_DISTANCE_TO_SKIN_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/checks.h"
#include "processes/process.h"
#include "processes/find_intersected_geometrical_objects_process.h"


namespace Kratos
{
  ///@addtogroup Kratos Core
  ///@{

  ///@name Kratos Classes
  ///@{

  /// This only calculates the distance. Calculating the inside outside should be done by a derived class of this.
  /** This process takes a volume model part (with tetrahedra mesh) and a skin model part (with triangle mesh) and
      and calcualtes the distance to the skin for all the elements and nodes of the volume model part.
  */
  template<std::size_t TDim = 3>
  class KRATOS_API(KRATOS_CORE) CalculateDiscontinuousDistanceToSkinProcess : public Process
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of CalculateDiscontinuousDistanceToSkinProcess
      KRATOS_CLASS_POINTER_DEFINITION(CalculateDiscontinuousDistanceToSkinProcess);

      ///@}
      ///@name Life Cycle
      ///@{

	  /// Constructor to be used.
	  CalculateDiscontinuousDistanceToSkinProcess(ModelPart& rVolumePart, ModelPart& rSkinPart);

	  /// Destructor.
      ~CalculateDiscontinuousDistanceToSkinProcess() override;


      ///@}
      ///@name Deleted
      ///@{

	  /// Default constructor.
	  CalculateDiscontinuousDistanceToSkinProcess() = delete;

	  /// Copy constructor.
	  CalculateDiscontinuousDistanceToSkinProcess(Process const& rOther) = delete;

	  /// Assignment operator.
	  CalculateDiscontinuousDistanceToSkinProcess& operator=(CalculateDiscontinuousDistanceToSkinProcess const& rOther) = delete;

	  /// Copy constructor.
	  CalculateDiscontinuousDistanceToSkinProcess(CalculateDiscontinuousDistanceToSkinProcess const& rOther);

      FindIntersectedGeometricalObjectsProcess mFindIntersectedObjectsProcess;

      ///@}
      ///@name Operations
      ///@{

      virtual void Initialize();

      virtual void FindIntersections();

      virtual std::vector<PointerVector<GeometricalObject>>& GetIntersections();

      virtual void CalculateDistances(std::vector<PointerVector<GeometricalObject>>& rIntersectedObjects);

      virtual void Clear();

	  void Execute() override;

      ///@}
      ///@name Access
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

    private:

      ///@name Member Variables
      ///@{

        ModelPart& mrSkinPart;
        ModelPart& mrVolumePart;

      ///@}
      ///@name Private Operations
      ///@{

		void CalculateElementalDistances(Element& rElement1, PointerVector<GeometricalObject>& rIntersectedObjects);

		double CalculateDistanceToNode(Element& rElement1, int NodeIndex, PointerVector<GeometricalObject>& rIntersectedObjects, const double Epsilon);

    unsigned int ComputeEdgesIntersections(
      Element& rElement1, 
      const PointerVector<GeometricalObject>& rIntersectedObjects,
      std::vector<unsigned int> &rCutEdgesVector,
      std::vector<array_1d <double,3> > &rIntersectionPointsArray);

    int ComputeEdgeIntersection(
      const Element::GeometryType& rIntObjGeometry,
      const Element::NodeType& rEdgePoint1,
      const Element::NodeType& rEdgePoint2, 
      Point& rIntersectionPoint);

    void ComputeIntersectionNormal(
      Element::GeometryType& rGeometry,
      const Vector& rElementalDistances,
      array_1d<double,3> &rNormal);

    void ComputePlaneApproximation(
      const Element& rElement1,
      const std::vector< array_1d<double,3> >& rPointsCoord,
      array_1d<double,3>& rPlaneBasePointCoords,
      array_1d<double,3>& rPlaneNormal);

    void CorrectDistanceOrientation(
      Element::GeometryType& rGeometry,
      const PointerVector<GeometricalObject>& rIntersectedObjects,
      Vector& rElementalDistances
    );

    void inline ComputeIntersectionNormalFromGeometry(
      const Element::GeometryType &rGeometry,
      array_1d<double,3> &rIntObjNormal);

    Plane3D inline SetIntersectionPlane(
      const std::vector<array_1d<double,3>> &rIntPtsVector);

      ///@}

    }; // Class CalculateDiscontinuousDistanceToSkinProcess

  ///@}

  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    CalculateDiscontinuousDistanceToSkinProcess<>& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const CalculateDiscontinuousDistanceToSkinProcess<>& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CALCULATE_DISCONTINUOUS_DISTANCE_TO_SKIN_PROCESS_H_INCLUDED  defined
