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

#if !defined(KRATOS_CALCULATE_DISTANCE_TO_SKIN_PROCESS_H_INCLUDED )
#define  KRATOS_CALCULATE_DISTANCE_TO_SKIN_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "processes/find_intersected_geometrical_objects_process.h"


namespace Kratos
{
  ///@addtogroup Kratos Core
  ///@{

  ///@name Kratos Globals
  ///@{

  ///@}
  ///@name Type Definitions
  ///@{

  ///@}
  ///@name  Enum's
  ///@{

  ///@}
  ///@name  Functions
  ///@{

  ///@}
  ///@name Kratos Classes
  ///@{

  /// This only calculates the distance. Calculating the inside outside should be done by a derived class of this.
  /** This process takes a volume model part (with tetrahedra mesh) and a skin model part (with triangle mesh) and
      and calcualtes the distance to the skin for all the elements and nodes of the volume model part.
  */
  class KRATOS_API(KRATOS_CORE) CalculateDistanceToSkinProcess : public FindIntersectedGeometricalObjectsProcess
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of CalculateDistanceToSkinProcess
      KRATOS_CLASS_POINTER_DEFINITION(CalculateDistanceToSkinProcess);

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor.
      CalculateDistanceToSkinProcess() = delete;

	  /// Copy constructor.
	  CalculateDistanceToSkinProcess(FindIntersectedGeometricalObjectsProcess const& rOther) = delete;

	  /// Constructor to be used.
	  CalculateDistanceToSkinProcess(ModelPart& rVolumePart, ModelPart& rSkinPart);

	  /// Destructor.
      virtual ~CalculateDistanceToSkinProcess();


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{

	  virtual void Execute() override;

      ///@}
      ///@name Access
      ///@{

	  ModelPart& GetSkinRepresentation() { return mSkinRepresentation; }


      ///@}
      ///@name Inquiry
      ///@{


      ///@}
      ///@name Input and output
      ///@{

      /// Turn back information as a string.
      virtual std::string Info() const;

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const;

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const;


      ///@}
      ///@name Friends
      ///@{


      ///@}

    protected:
      ///@name Protected static Member Variables
      ///@{


      ///@}
      ///@name Protected member Variables
      ///@{


      ///@}
      ///@name Protected Operators
      ///@{


      ///@}
      ///@name Protected Operations
      ///@{


      ///@}
      ///@name Protected  Access
      ///@{


      ///@}
      ///@name Protected Inquiry
      ///@{


      ///@}
      ///@name Protected LifeCycle
      ///@{


      ///@}

    private:

		class LineSegment : public std::array<Point<3>::Pointer, 2>{
		public:
			LineSegment(Point<3>::Pointer pPoint1, Point<3>::Pointer pPoint2) {
				this->operator[](0) = pPoint1;
				this->operator[](1) = pPoint2;
			}

			Point<3> const& GetPoint(int PointId) {
				return *(this->operator[](PointId).get());
			}

			int TriangleIntersectionPoint(Element::GeometryType& rGeometry, Point<3>& IntersectionPoint)
			{
				// This is the adaption of the implemnetation provided in:
				// http://www.softsurfer.com/Archive/algorithm_0105/algorithm_0105.htm#intersect_RayTriangle()
				// with additional segment range check

				const double epsilon = 1.00e-12;

				Point<3> const& r_point_1 = GetPoint(0);
				Point<3> const& r_point_2 = GetPoint(1);


				array_1d<double, 3>    u, v, n;             // triangle vectors
				array_1d<double, 3>    dir, w0, w;          // ray vectors
				double     r, a, b;             // params to calc ray-plane intersect


												// get triangle edge vectors and plane normal
				u = rGeometry[1] - rGeometry[0];
				v = rGeometry[2] - rGeometry[0];

				MathUtils<double>::CrossProduct(n, u, v);             // cross product

				if (norm_2(n) == 0)            // triangle is degenerate
					return -1;                 // do not deal with this case

				for (int i = 0; i < 3; i++)
				{
					dir[i] = r_point_2[i] - r_point_1[i];             // ray direction vector
					w0[i] = r_point_1[i] - rGeometry[0][i];
				}

				a = -inner_prod(n, w0);
				b = inner_prod(n, dir);

				if (fabs(b) < epsilon) {     // ray is parallel to triangle plane
					if (a == 0)                // ray lies in triangle plane
						return 2;
					else return 0;             // ray disjoint from plane 
				}

				// get intersect point of ray with triangle plane
				r = a / b;
				if (r < 0.0)                   // ray goes away from triangle
					return 0;                  // => no intersect
				if (r > 1.0) // for a segment, also test if (r > 1.0) => no intersect
					return 0; 

				for (int i = 0; i < 3; i++)
					IntersectionPoint[i] = r_point_1[i] + r * dir[i];           // intersect point of ray and plane

																				// is I inside T?
				double    uu, uv, vv, wu, wv, D;
				uu = inner_prod(u, u);
				uv = inner_prod(u, v);
				vv = inner_prod(v, v);


				for (int i = 0; i < 3; i++)
					w[i] = IntersectionPoint[i] - rGeometry[0][i];


				wu = inner_prod(w, u);
				wv = inner_prod(w, v);
				D = uv * uv - uu * vv;

				// get and test parametric coords
				double s, t;
				s = (uv * wv - vv * wu) / D;
				if (s < 0.0 - epsilon || s > 1.0 + epsilon)        // I is outside T
					return 0;
				t = (uv * wu - uu * wv) / D;
				if (t < 0.0 - epsilon || (s + t) > 1.0 + epsilon)  // I is outside T
					return 0;

				return 1;                      // I is in T

			}

		};

      ///@name Static Member Variables
      ///@{


      ///@}
      ///@name Member Variables
      ///@{

		ModelPart mSkinRepresentation;

      ///@}
      ///@name Private Operators
      ///@{


      ///@}
      ///@name Private Operations
      ///@{

		void CalculateElementDistance(Element& rElement1, PointerVector<GeometricalObject>& rIntersectedObjects);
		int CalculateIntersectionPoints(LineSegment& rSegment, PointerVector<GeometricalObject>& rIntersectedObjects, std::vector<Point<3> >& IntersectionPoints);

      ///@}
      ///@name Private  Access
      ///@{


      ///@}
      ///@name Private Inquiry
      ///@{


      ///@}
      ///@name Un accessible methods
      ///@{

      /// Assignment operator.
      CalculateDistanceToSkinProcess& operator=(CalculateDistanceToSkinProcess const& rOther);

      /// Copy constructor.
      CalculateDistanceToSkinProcess(CalculateDistanceToSkinProcess const& rOther);


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
				    CalculateDistanceToSkinProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const CalculateDistanceToSkinProcess& rThis)
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
