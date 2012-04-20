/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
mossaiby@yahoo.com
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/


// System includes 

// External includes 
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "spaces/ublas_space.h"

#include "custom_python/add_custom_utilities_to_python.h"



// #include "spatial_containers/spatial_containers.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "includes/model_part.h"
#include "custom_utilities/bins_dynamic_mpi.h"
#include "custom_utilities/bins_dynamic_objects_mpi.h"

  template< std::size_t dim_type>
  class Point {

	public:

		double       coord[dim_type];
		std::size_t  id;
		std::size_t  tag;
		//int id;

		double& operator[](std::size_t i) {return coord[i];}

		double const & operator[](std::size_t i) const {return coord[i];}

		void RandomCoord(){
		  for(std::size_t i = 0 ; i < dim_type ; i++)
			  coord[i] = double(rand())/RAND_MAX;
		}

		void operator=(Point<dim_type> const& Other){
		  for(std::size_t i = 0; i < dim_type; i++)
			  coord[i] = Other.coord[i];
		}
  };

  template< std::size_t dim_type >
  std::ostream & operator<<( std::ostream& rOut, Point<dim_type> & rPoint){
	  rOut << "(" << rPoint.id << ") ";
	for(std::size_t i = 0 ; i < dim_type ; i++)
		rOut << rPoint[i] << " "; 
	return rOut; 
  };

  template< std::size_t dim_type >
  std::istream & operator>>( std::istream& rIn, Point<dim_type> & rPoint){
	for(std::size_t i = 0 ; i < dim_type ; i++)
		rIn >> rPoint[i]; 
	return rIn; 
  };

  template< class T, std::size_t dim >
  class PointDistance{
	public:
		double operator()( T const& p1, T const& p2 ){
		  double dist = 0.0;
		  for( std::size_t i = 0 ; i < dim ; i++){
			  double tmp = p1[i] - p2[i];
			  dist += tmp*tmp;
		  }
		  return sqrt(dist);
		}
  };

  template< class T, std::size_t dim >
  class PointDistance2{
	public:
		double operator()( T const& p1, T const& p2 ){
		  double dist = 0.0;
		  for( std::size_t i = 0 ; i < dim ; i++){
			  double tmp = p1[i] - p2[i];
			  dist += tmp*tmp;
		  }
		  return dist;
		}
  };

  template< std::size_t dim >
  bool LowerPoint( Point<dim> const& reference, Point<dim> const& new_ ){
	for(std::size_t i = 0 ; i < dim ; i++)
		if( reference[i] < new_[i] )
		  return false;
	return true;
  };

  template< std::size_t dim >
  bool UpperPoint( Point<dim> const& reference, Point<dim> const& new_ ){
	for(std::size_t i = 0 ; i < dim ; i++)
		if( reference[i] > new_[i] )
		  return false;
	return true;
  };

  // Spheres

  template< std::size_t dim_type>
  class Sphere {

	public:

	  Point<dim_type> center;
	  double          radius;

		double& operator[](std::size_t i) {return center[i];}

		double const & operator[](std::size_t i) const {return center[i];}

		void RandomCoord(){
		  center.RandomCoord();
		}

		void operator=(Sphere<dim_type> const& Other){
		  center = Other.center;
		  radius = Other.radius;
		}
  };


  template< std::size_t dim_type >
  std::ostream & operator<<( std::ostream& rOut, Sphere<dim_type> & rSphere){
	rOut << rSphere.center << " " << rSphere.radius;
	return rOut; 
  };

  template< std::size_t dim_type >
  std::istream & operator>>( std::istream& rIn, Sphere<dim_type> & rSphere){
	rIn >> rSphere.center.id >> rSphere.center >> rSphere.radius;
	return rIn; 
  };

  template< std::size_t dim >
  double DistanceSphereSphere( const Sphere<dim>* sphere1, const Sphere<dim>* sphere2 )
  {
	double dist_center = PointDistance2<Point<dim>,dim>()(sphere1->center,sphere2->center);
	double radius = (sphere1->radius + sphere2->radius);
	return abs(dist_center - radius*radius);
  }

  template< std::size_t dim >
  bool IntersectionSphereSphere( const Sphere<dim>* sphere1, const Sphere<dim>* sphere2 )
  {
	double dist_center = PointDistance2<Point<dim>,dim>()(sphere1->center,sphere2->center);
	double radius = (sphere1->radius + sphere2->radius) * 1.01;
	return (dist_center < radius*radius);
  }

  template< std::size_t dim >
  void SphereBoundingBox( const Sphere<dim>* sphere, Point<dim>& Low, Point<dim>& High )
  {
	for(std::size_t i = 0 ; i < dim ; i++)
	{
	  High[i] = sphere->center[i] + sphere->radius;
	  Low[i]  = sphere->center[i] - sphere->radius;
	}
  }
  
  template< std::size_t dim >
  void SphereBoundingBox( const Sphere<dim>* sphere, Point<dim>& Low, Point<dim>& High, std::size_t Radius )
  {
	for(std::size_t i = 0 ; i < dim ; i++)
	{
	  High[i] = sphere->center[i] + Radius;
	  Low[i]  = sphere->center[i] - Radius;
	}
  }

  template< std::size_t dim >
  bool IntersectionSphereBox( Sphere<dim> const* sphere, Point<dim> const& LowerPoint, Point<dim> const& HighPoint )
  {
	double dmin = 0.00;
	double tmp, radius2;
	for(std::size_t i = 0 ; i < dim ; i++)
	{
	  tmp = std::max(LowerPoint[i] - sphere->center[i],0.00) + std::max(sphere->center[i] - HighPoint[i],0.00);
	  if( tmp > sphere->radius )
		return false;
	  dmin += tmp*tmp;
	}
	radius2 = sphere->radius * sphere->radius;
	return (dmin <= radius2);
  };

  template< class T>
  class Pair
  {
	public:
	  Pair(T const& d1, T const& d2){
		data[0] = d1;
		data[1] = d2;
	  }

	  Pair( Pair<T> const& p ) {
		data[0] = p.data[0];
		data[1] = p.data[1];
	  }

	  Pair(){}

	  ~Pair(){}

	  T const& operator[](std::size_t index) const {
		return data[index];
	  }

	  T& operator[](std::size_t index) {
		return data[index];
	  }

	  Pair<T>& operator = ( Pair<T> const& p ){
		data[0] = p.data[0];
		data[1] = p.data[1];
		return *this;
	  }

	  bool operator==( Pair<T> const& p ){
		return ( (data[0]==p.data[0]) && (data[1]==p.data[1]) );
	  }

	private:
	  T data[2];

  };
  
  class SpheresSpatialConfigure
  {
	public:
	  

	static const std::size_t Dim = 3;
	static const std::size_t Dimension = 3;
	
	typedef std::size_t  SizeType;
	typedef std::size_t  IndexType;

	typedef Point<Dim> PointType;    //typedef typename TConfigure::PointType  PointType;
	typedef Sphere<Dim> SphereType;

	typedef SphereType*          			PointerType;   //    PtrSphereType;
	typedef PointerType*            		ContainerType; //    SphereVector;
	typedef PointerType*            		IteratorType;  //    SphereIterator;

	typedef double*                 		DistanceContainerType; // DistanceVector;
	typedef double*                 		DistanceIteratorType;  // DistanceIterator;

	typedef ContainerType           		ResultContainerType;
	typedef IteratorType            		ResultIteratorType;

	typedef Pair<PointerType> 				ContactPairType;

	typedef std::vector<ContactPairType>  	ContainerContactType;
	typedef ContainerContactType::iterator  IteratorContactType;

	static inline bool Intersection(const PointerType& p1, const PointerType& p2){
	  return IntersectionSphereSphere<Dim>(p1,p2);
	}

	static inline void CalculateBoundingBox( PointerType& p, PointType& pmin, PointType& pmax){
	  SphereBoundingBox<Dim>(p,pmin,pmax);
	}
	
	static inline void CalculateBoundingBox( PointerType& p, PointType& pmin, PointType& pmax, SizeType Radius){
	  SphereBoundingBox<Dim>(p,pmin,pmax,Radius);
	}

	static inline bool IntersectionBox(const PointerType& p, const PointType& pmin, const PointType& pmax){
	  return IntersectionSphereBox<Dim>(p,pmin,pmax);
	}
	
	static inline double Distance(const PointerType& p1, const PointerType& p2){
	  return DistanceSphereSphere<Dim>(p1,p2);
	}


  };

namespace Kratos
{
	namespace Python
	{	
		void AddCustomUtilitiesToPython()
		{
	  
			static const std::size_t Dim = 3;

			typedef Kratos::Point<Dim> 					PointType;
			typedef Kratos::Tetrahedra3D4<PointType>	ObjectType;

			typedef PointType*							PtrPointType;
			typedef ObjectType*							PtrObjectType;

			typedef PtrPointType*						PointVector;
			typedef PtrPointType*						PointIterator;

			typedef PtrObjectType*						ObjectVector;
			typedef PtrObjectType*						ObjectIterator;

			typedef double*								DistanceVector;
			typedef double*								DistanceIterator;
	  
			typedef Kratos::SearchUtils::SquaredDistanceFunction<Dim,PointType> DistanceFunction;
			typedef Kratos::BinsDynamicMpi< Dim, PointType, PointVector, PtrPointType, PointIterator, DistanceIterator, DistanceFunction > BinsDynamicMpi;
			typedef Kratos::BinsObjectDynamicMpi<SpheresSpatialConfigure> BinsObjectDynamicMpi;
			using namespace boost::python;
		  
			class_< BinsDynamicMpi, boost::noncopyable > ("BinsDynamicMpi", init<Kratos::ModelPart * , Kratos::ModelPart *, double>() )
				.def("MultiSearchInRadiusTest"	, &BinsDynamicMpi::MultiSearchInRadiusTest)
				.def("MPISingleSearchInRadiusTest", &BinsDynamicMpi::MPISingleSearchInRadiusTest)
				.def("MPIMultiSearchInRadiusTest", &BinsDynamicMpi::MPIMultiSearchInRadiusTest)
			;
		  
			class_< BinsObjectDynamicMpi, boost::noncopyable > ("BinsObjectDynamicMpi", init<Kratos::ModelPart * , Kratos::ModelPart *, double>())
				.def("SingleSearchObjectsInRadiusTest", &BinsObjectDynamicMpi::SingleSearchObjectsInRadiusTest)
			;
	}
	  
  }  // namespace Python.

} // Namespace Kratos

