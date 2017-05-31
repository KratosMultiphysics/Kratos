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

#if !defined(KRATOS_CALCULATE_DISCONTINUOUS_DISTANCE_TO_SKIN_PROCESS_H_INCLUDED )
#define  KRATOS_CALCULATE_DISCONTINUOUS_DISTANCE_TO_SKIN_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "processes/find_intersected_geometrical_objects_process.h"
#include "includes/checks.h"


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
  class KRATOS_API(KRATOS_CORE) CalculateDiscontinuousDistanceToSkinProcess : public FindIntersectedGeometricalObjectsProcess
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
      virtual ~CalculateDiscontinuousDistanceToSkinProcess();


      ///@}
      ///@name Deleted 
      ///@{

	  /// Default constructor.
	  CalculateDiscontinuousDistanceToSkinProcess() = delete;

	  /// Copy constructor.
	  CalculateDiscontinuousDistanceToSkinProcess(FindIntersectedGeometricalObjectsProcess const& rOther) = delete;

	  /// Assignment operator.
	  CalculateDiscontinuousDistanceToSkinProcess& operator=(CalculateDiscontinuousDistanceToSkinProcess const& rOther) = delete;

	  /// Copy constructor.
	  CalculateDiscontinuousDistanceToSkinProcess(CalculateDiscontinuousDistanceToSkinProcess const& rOther);

      ///@}
      ///@name Operations
      ///@{

	  virtual void Execute() override;

      ///@}
      ///@name Access
      ///@{

	  ModelPart& GetSkinRepresentation() { return mSkinRepresentation; }

      ///@}
      ///@name Input and output
      ///@{

      /// Turn back information as a string.
      virtual std::string Info() const override;

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const override;

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const override;

      ///@}

    private:

		// TODO: I should move this class to a separate file but is out of scope of this branch
		class Plane3D {
		public:
			using VectorType = array_1d<double, 3>;
			using PointType = Point<3>;

			Plane3D(VectorType const& TheNormal, double DistanceToOrigin) :mNormal(TheNormal), mD(DistanceToOrigin) {}
			Plane3D() = delete;
			Plane3D(PointType const& Point1, PointType const& Point2, PointType const& Point3) {
				VectorType v1 = Point2 - Point1;
				VectorType v2 = Point3 - Point1;
				MathUtils<double>::CrossProduct(mNormal, v1, v2);
				auto normal_length = norm_2(mNormal);
				KRATOS_DEBUG_CHECK_GREATER(normal_length, std::numeric_limits<double>::epsilon());
				mNormal /= normal_length;
				mD = -inner_prod(mNormal, Point1);
			}
			VectorType const& GetNormal() { return mNormal; }
			double GetDistance() { return mD; }
			double CalculateSignedDistance(PointType const& ThePoint) {
				return inner_prod(mNormal, ThePoint) + mD;
			}

		private:
			VectorType mNormal;
			double mD;
		};

      ///@name Member Variables
      ///@{

		ModelPart mSkinRepresentation;

      ///@}
      ///@name Private Operations
      ///@{

		void CalculateElementalDistances(Element& rElement1, PointerVector<GeometricalObject>& rIntersectedObjects);
		
		double CalculateDistanceToNode(Element& rElement1, int NodeIndex, PointerVector<GeometricalObject>& rIntersectedObjects, const double Epsilon);

      ///@}

    }; // Class CalculateDiscontinuousDistanceToSkinProcess

  ///@}

  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    CalculateDiscontinuousDistanceToSkinProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const CalculateDiscontinuousDistanceToSkinProcess& rThis)
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
