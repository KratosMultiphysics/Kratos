/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
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
 
//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: janosch $
//   Date:                $Date: 2008-10-23 11:30:49 $
//   Revision:            $Revision: 1.12 $
//
//


#if !defined(KRATOS_GEOMETRY_H_INCLUDED )
#define  KRATOS_GEOMETRY_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <sstream>
#include <cstddef>


// External includes
#include <boost/array.hpp> 


// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "containers/pointer_vector.h"
#include "integration/integration_point.h"
#include "geometries/geometry_data.h"
#include "spatial_containers/bounding_box.h"
#include "utilities/math_utils.h"


namespace Kratos
{
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
    
    ///Geometry base class.
    /** As a base class Geometry has all the common
     * interface of Kratos' geometries. Also it contains array of
     * pointers to its points, reference to shape functions values in
     * all integrations points and also local gradients of shape
     * functions evaluated in all integrations points.
     * 
     * Geometry is a template class with just one template parameter: 
     * - TPointType which reperesent the type of the point this geometry
     * type contain and build on.
     * 
     * @see Point
     * @see Node
     * @see Formulation
     * @see GeometryAndFormulationElement
     */
    template<class TPointType> 
    class Geometry : public PointerVector<TPointType>
    {
    public:
      ///@}
      ///@name Type Definitions
      ///@{
      
      /// This Geometry type.
      typedef Geometry<TPointType> GeometryType;

      /// Pointer definition of Geometry
      KRATOS_CLASS_POINTER_DEFINITION(Geometry);
  
      /** Base type for geometry.
      */
      typedef PointerVector<TPointType> BaseType;

  
      /** The bounding box */
      typedef BoundingBox<TPointType, GeometryType>  BoundingBoxType;      

      /** Array of counted pointers to point. This type used to hold
	  geometry's points.
      */
      typedef PointerVector<TPointType> PointsArrayType;

      /** Integration methods implemented in geometry. 
      */
      typedef GeometryData::IntegrationMethod IntegrationMethod;

     /** A Vector of counted pointers to Geometries. Used for
	  returning edges of the geometry.
      */
      typedef PointerVector<GeometryType> GeometriesArrayType;
      
      /** Redefinition of geometry template parameter TPointType as this geometry point type.
       */
      typedef TPointType PointType;
      
      /** Type used for indexing in geometry class.std::size_t used for indexing
	  point or integration point access methods and also all other
	  methods which need point or integration point index.
      */
      typedef std::size_t IndexType;
      
      
      /** This typed used to return size or dimension in
	  geometry. Dimension, WorkingDimension, PointsNumber and
	  ... return this type as their results.
      */
      typedef std::size_t SizeType;

	  typedef typename PointType::CoordinatesArrayType CoordinatesArrayType;
      
     
      /** This type used for representing an integration point in
	  geometry. This integration point is a point with an
	  additional weight component.
      */
      typedef IntegrationPoint<3> IntegrationPointType;
      
      /** A Vector of IntegrationPointType which used to hold
	  integration points related to an integration
	  method. IntegrationPoints functions used this type to return
	  their results.
      */
      typedef std::vector<IntegrationPointType> IntegrationPointsArrayType;
      
      /** A Vector of IntegrationPointsArrayType which used to hold
	  integration points related to different integration method
	  implemented in geometry.
      */
      typedef boost::array<IntegrationPointsArrayType, GeometryData::NumberOfIntegrationMethods> IntegrationPointsContainerType;
      
      /** A third order tensor used as shape functions' values
	  continer.
      */
      typedef boost::array<Matrix, GeometryData::NumberOfIntegrationMethods> ShapeFunctionsValuesContainerType;
      
      /** A fourth order tensor used as shape functions' local
	  gradients container in geometry.
      */
	  typedef GeometryData::ShapeFunctionsLocalGradientsContainerType ShapeFunctionsLocalGradientsContainerType;
      
      /** A third order tensor to hold jacobian matrices evaluated at
	  integration points. Jacobian and InverseOfJacobian functions
	  return this type as their result.
      */
      typedef boost::numeric::ublas::vector<Matrix > JacobiansType;
      
      /** A third order tensor to hold shape functions'  gradients. 
	  ShapefunctionsGradients function return this
	  type as its result.
      */
      typedef GeometryData::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;
      
      /** A third order tensor to hold shape functions' local second derivatives.
	  ShapefunctionsLocalGradients function return this
	  type as its result.
      */
      typedef GeometryData::ShapeFunctionsSecondDerivativesType ShapeFunctionsSecondDerivativesType;
      
      /** A fourth order tensor to hold shape functions' local third order derivatives
       */
      typedef GeometryData::ShapeFunctionsThirdDerivativesType ShapeFunctionsThirdDerivativesType;
      
      /** Type of the normal vector used for normal to edges in geomety.
       */
      typedef boost::numeric::ublas::vector<double> NormalType;

      
      typedef typename BaseType::iterator              iterator;
    typedef typename BaseType::const_iterator          const_iterator;                
    typedef typename BaseType::reverse_iterator        reverse_iterator;
    typedef typename BaseType::const_reverse_iterator  const_reverse_iterator;              

    typedef typename BaseType::ptr_iterator ptr_iterator;
    typedef typename BaseType::ptr_const_iterator ptr_const_iterator;
    typedef typename BaseType::ptr_reverse_iterator ptr_reverse_iterator;
    typedef typename BaseType::ptr_const_reverse_iterator ptr_const_reverse_iterator;
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      Geometry() : mpGeometryData(0)
	{
	  
	}
      /** Complete argument constructor. This constructor gives a
	  complete set of arguments to pass all the initial value of
	  all the member variables of geometry class. Also it has
	  default value for integration variables to make it usefull
	  in the case of constructing new geometry without mapping and
	  integrating properties.

	  @param ThisPoints Vector of pointers to points which this
	  geometry constructing on them. Points must have dimension
	  equal or greater than working space dimension though there
	  is no control on it.  

	  @param ThisDefaultMethod Default integration method. Its
	  default value is gaussian integration with orden one which
	  make no deference while in this condition there is no shape
	  function database exist and integrating is not possible
	  including by default method.

	  @param ThisIntegrationPoints All the integration points in
	  all methods. This is a Vector of IntegrationPointsArrayType
	  and It must have at least four component correspounding to
	  four integration method defined now. If there is some
	  geometry which don't have all this method implemented
	  related points Vector must exist but with zero size. For
	  example if a geometry don't have gaussian orden one
	  ThisIntegrationPoints[GI_GAUSS_1] must be an empty
	  IntegrationPointsArrayType.

	  @param ThisShapeFunctionsValues Values of all the shape
	  functions evaluated in all integrations points of all
	  integration methods. It's a three dimensional array \f$
	  F_{ijk} \f$ where i = GI_GAUSS_1,..., GI_GAUSS_4 and j is
	  the integration point index and k is the shape function
	  index. In the other word component \f$ f_{ijk} \f$ is the
	  value of the shape function related to node k evaluated in
	  integration point j of i integration method point set. Again
	  if there is some integration method unsupported an empty
	  Matrix must assigned to related place. For example if a
	  geometry don't have gaussian orden four
	  ThisShapeFunctionsValues[GI_GAUSS_4] must be an empty
	  Matrix.

	  @param ThisShapeFunctionsLocalGradients Values of local
	  gradients respected to all local coordinates of all the
	  shape functions evaluated in all integrations points of all
	  integration methods. It's a four dimensional array \f$
	  F_{ijkh} \f$ where i = GI_GAUSS_1,..., GI_GAUSS_4 and j is
	  the integration point index and k is the shape function
	  index and h is local coordinate index. In the other word
	  component \f$ f_{ijkh} \f$ is the value of h'th component of
	  local gradient of the shape function related to node k
	  evaluated in integration point j of i integration method
	  point set. Again if there is some integration method
	  unsupported an empty ShapeFunctionsGradientsType must
	  assigned to related place. For example if a geometry don't
	  have gaussian orden two ThisShapeFunctionsValues[GI_GAUSS_2]
	  must be an empty ShapeFunctionsGradientsType.
      */
      Geometry(const PointsArrayType& ThisPoints,
	       GeometryData const* pThisGeometryData= 0) 
	: BaseType(ThisPoints)
	, mpGeometryData(pThisGeometryData)
	{
	}

//       Geometry(const PointsArrayType& ThisPoints,
// 	       GeometryData const& ThisGeometryData= msEmptyGeometryData) 
// 	: BaseType(ThisPoints)
// 	, mGeometryData(ThisGeometryData)
// 	{
// 	}

      /** Copy constructor.
	  Construct this geometry as a copy of given geometry.

	  @note This copy constructor don't copy the points and new
	  geometry shares points with given source geometry. It's
	  obvious that any change to this new geometry's point affect
	  source geometry's points too.
      */
      Geometry(const Geometry& rOther)
	: BaseType(rOther)
	, mpGeometryData(rOther.mpGeometryData)
	{
	}


      /** Copy constructor from a geometry with other point type.
	  Construct this geometry as a copy of given geometry which
	  has different type of points. The given goemetry's
	  TOtherPointType* must be implicity convertible to this
	  geometry PointType.

	  @note This copy constructor don't copy the points and new
	  geometry shares points with given source geometry. It's
	  obvious that any change to this new geometry's point affect
	  source geometry's points too.
      */
      template<class TOtherPointType> Geometry(Geometry<TOtherPointType> const & rOther)
	: BaseType(rOther.begin(), rOther.end())
	, mpGeometryData(rOther.mpGeometryData)
	{
	}

      /// Destructor. Do nothing!!!
      virtual ~Geometry(){}
      
	virtual GeometryData::KratosGeometryFamily GetGeometryFamily(){return GeometryData::Kratos_generic_family; }
	virtual GeometryData::KratosGeometryType GetGeometryType(){return GeometryData::Kratos_generic_type; }

      ///@}
      ///@name Operators 
      ///@{
      
      /** Assignment operator.

      @note This operator don't copy the points and this
      geometry shares points with given source geometry. It's
      obvious that any change to this geometry's point affect
      source geometry's points too.
      
      @see Clone
      @see ClonePoints
      */
      Geometry& operator=(const Geometry& rOther)
	{
	  BaseType::operator=(rOther);
	  mpGeometryData = rOther.mpGeometryData;

	  return *this;
	}

      /** Assignment operator for geometries with different point type.

      @note This operator don't copy the points and this
      geometry shares points with given source geometry. It's
      obvious that any change to this geometry's point affect
      source geometry's points too.
      
      @see Clone
      @see ClonePoints
      */
      template<class TOtherPointType> 
	Geometry& operator=(Geometry<TOtherPointType> const & rOther)
	{
	  this->clear();
	  for(typename Geometry<TOtherPointType>::ptr_const_iterator i = rOther.ptr_begin() ; i != rOther.ptr_end() ; ++i)
	    push_back(typename PointType::Pointer(new PointType(**i)));

	  mpGeometryData = rOther.mpGeometryData;

	  return *this;
	}
      
      ///@}
      ///@name Operations
      ///@{
      
      virtual Pointer Create(PointsArrayType const& ThisPoints) const
	{
	  return Pointer(new Geometry(ThisPoints, mpGeometryData));
	}

      /** This methods will create a duplicate of all its points and
	  substitute them with its points. */
      void ClonePoints()
	{
	  for(ptr_iterator i = this->ptr_begin() ; i != this->ptr_end() ; i++)
	    *i = typename PointType::Pointer(new PointType(**i));
	}

      virtual boost::shared_ptr< Geometry< Point<3> > > Clone() const
	{
	  Geometry< Point<3> >::PointsArrayType NewPoints;

	  //making a copy of the nodes TO POINTS (not Nodes!!!)
	  for(IndexType i = 0 ; i < this->size() ; i++)
	    NewPoints.push_back((*this)[i]);
	    //NewPoints[i] = typename Point<3>::Pointer(new Point<3>(*mPoints[i]));

	  //creating a geometry with the new points
	  boost::shared_ptr< Geometry< Point<3> > > p_clone(new Geometry< Point<3> >(NewPoints));
	  p_clone->ClonePoints();

	  return p_clone;
	}
      
      //lumping factors for the calculation of the lumped mass matrix
      virtual Vector& LumpingFactors(Vector& rResult)  const
	{ 
	  KRATOS_ERROR(std::logic_error, "Called the virtual function for LumpingFactors" , *this);
	  return rResult; 
	}

      ///@}
      ///@name Informations
      ///@{

      /** Dimension of the geometry for example a triangle2d is a 2
	  dimensional shape

	  @return SizeType, dimension of this geometry. 
	  @see WorkingSpaceDimension()
	  @see LocalSpaceDimension()
      */
      inline SizeType Dimension() const
	{
	  return mpGeometryData->Dimension();
	}

      /** Working space dimension. for example a triangle is a 2
	  dimensional shape but can be used in 3 dimensional space.

	  @return SizeType, working space dimension of this geometry.
	  @see Dimension()
	  @see LocalSpaceDimension()
      */
      inline SizeType WorkingSpaceDimension() const
	{
	  return mpGeometryData->WorkingSpaceDimension();
	}

      /** Local space dimension. for example a triangle is a 2
	  dimensional shape but can have 3 dimensional area
	  coordinates l1, l2, l3.

	  @return SizeType, local space dimension of this geometry.
	  @see Dimension()
	  @see WorkingSpaceDimension()
      */
      inline SizeType LocalSpaceDimension() const
	{
	  return mpGeometryData->LocalSpaceDimension();
	}

      /** Returns number of the points which this geometry has.

      @return SizeType, number of the points in this geometry.
      */
      SizeType PointsNumber() const
	{
	  return this->size();
	}

      /** This method calculate and return Length or charactereistic
	  length of this geometry depending to it's dimension. For one
	  dimensional geometry for example Line it returns length of it
	  and for the other geometries it gives Characteristic length
	  otherwise.

	  @return double value contains length or Characteristic
	  length
	  @see Area()
	  @see Volume()
	  @see DomainSize()
      */
      virtual double Length() const
	{
	  return 0;
	}

      /** This method calculate and return area or surface area of
	  this geometry depending to it's dimension. For one dimensional
	  geometry it returns zero, for two dimensional it gives area
	  and for three dimensional geometries it gives surface area.

	  @return double value contains area or surface
	  area.
	  @see Length()
	  @see Volume()
	  @see DomainSize()
      */
      virtual double Area() const
	{
	  return 0;
	}

      /** This method calculate and return volume of this
	  geometry. For one and two dimensional geometry it returns
	  zero and for three dimensional it gives volume of geometry.

	  @return double value contains volume.
	  @see Length()
	  @see Area()
	  @see DomainSize()
      */
      virtual double Volume() const
	{
	  return 0;
	}


      /** This method calculate and return length, area or volume of
	  this geometry depending to it's dimension. For one dimensional
	  geometry it returns its length, for two dimensional it gives area
	  and for three dimensional geometries it gives its volume.

	  @return double value contains length, area or volume.
	  @see Length()
	  @see Area()
	  @see Volume()
      */
      virtual double DomainSize() const
	{
	  return 0;
	}



/*     
     virtual typename  BoundingBoxType::Pointer Create
               (const TPointType& LowPoint, const TPointType& HighPoint) //  GeometryType* ThisGeom )
	{
            
                return typename BoundingBoxType::Pointer(new BoundingBoxType(LowPoint, HighPoint, this)) ; 	
                //mpBondingBoxType = new BoundingBoxType(LowPoint, HighPoint, this ) ;
                //return mpBondingBoxType;    
	}
*/     
        

     virtual void Bounding_Box(BoundingBoxType& rResult) const
             {   
                
                Bounding_Box(rResult.LowPoint(), rResult.HighPoint());  
             }

 
     virtual void Bounding_Box(TPointType& rLowPoint, TPointType& rHighPoint) const 
            { 
               
                               
		rHighPoint = this->GetPoint(0);
		rLowPoint  = this->GetPoint(0);           
                 
               for (unsigned int point = 1; point<PointsNumber(); point++)      
                 {
                    for(unsigned int i = 0; i<3; i++)
                      {
                         this->GetPoint(point)[i];  
                         rHighPoint[i] =  (rHighPoint[i] < this->GetPoint(point)[i] ) ? this->GetPoint(point)[i] : rHighPoint[i];
                         rLowPoint[i]  =  (rLowPoint[i]  > this->GetPoint(point)[i] ) ? this->GetPoint(point)[i] : rLowPoint[i]; 
                      }
                  }
            } 
       
     

      /** Calculates center of this geometry by a simple averaging algorithm.
	  Each center point component calculated using:
	  \f[
	  c_i = \sum_j^n(x_i^j) / n
	  \f]

	  where \f$ c_i \f$ is component i of center point and \f$
	  X_i^j \f$ is component i of j'th point of geometry and n is
	  number of the points in this geometry.

	  @return PointType which is the calculated center of this geometry.
      */
      virtual PointType Center() const
	{
	  const SizeType points_number = this->size();
	  if(points_number == 0)
	    return PointType();

	  PointType result = (*this)[0];
	  for(IndexType i = 1 ; i < points_number ; i++)
	    result += (*this)[i];

	  double temp = 1.0/ double(points_number);
	  result *= temp;

	  return result;
	}

      ///@}
      ///@name Access
      ///@{ 
      
      /** A constant access method to the Vector of the points stored in
	  this geometry.

	  @return A constant reference to PointsArrayType contains
	  pointers to the points.
      */
      const PointsArrayType& Points() const
	{
	  return *this;
	}

      /** An access method to the Vector of the points stored in
	  this geometry.

	  @return A reference to PointsArrayType contains pointers to
	  the points.
      */
      PointsArrayType& Points()
	{
	  return *this;
	}

      /** A constant access method to the i'th points stored in
	  this geometry.

	  @return A constant counted pointer to i'th point of 
	  geometry.
      */
      const typename TPointType::Pointer pGetPoint(const int Index) const
	{
	  KRATOS_TRY_LEVEL_3
	    return (*this)(Index);
	  KRATOS_CATCH_LEVEL_3(*this)
	    }

      /** An access method to the i'th points stored in
	  this geometry.

	  @return A counted pointer to i'th point of 
	  geometry.
      */
      typename TPointType::Pointer pGetPoint(const int Index)
	{
	  KRATOS_TRY_LEVEL_3
	    return (*this)(Index);
	  KRATOS_CATCH_LEVEL_3(*this);
	}
      
      /** A constant access method to the i'th points stored in
	  this geometry.

	  @return A constant counted pointer to i'th point of 
	  geometry.
      */
      TPointType const& GetPoint(const int Index) const
	{
	  KRATOS_TRY_LEVEL_3
	    return (*this)[Index];
	  KRATOS_CATCH_LEVEL_3(*this)
	    }


      /** An access method to the i'th points stored in
	  this geometry.

	  @return A counted pointer to i'th point of 
	  geometry.
      */
      TPointType& GetPoint(const int Index)
	{
	  KRATOS_TRY_LEVEL_3
	    return (*this)[Index];
	  KRATOS_CATCH_LEVEL_3(*this);
	}
    
       /**
        * Returns a matrix of the local coordinates of all points 
        * @param rResult a Matrix that will be overwritten by the results
        * @return the coordinates of all points of the current geometry 
        */
       virtual Matrix& PointsLocalCoordinates( Matrix& rResult ) const 
       {
           KRATOS_ERROR(std::logic_error,
                        "Calling base class PointsLocalCoordinates method instead of drived class one. Please check the definition of derived class." , *this);
           return rResult;
       }
	   
	   /**
		* Returns the local coordinates of a given arbitrary point
		*/
	   virtual CoordinatesArrayType& PointLocalCoordinates( CoordinatesArrayType& rResult, 
			   const CoordinatesArrayType& rPoint )
	   {
		   Matrix J = ZeroMatrix( LocalSpaceDimension(), LocalSpaceDimension() );
		   if( rResult.size() != LocalSpaceDimension() )
			   rResult.resize( LocalSpaceDimension() );
		   //starting with xi = 0
		   rResult = ZeroVector( LocalSpaceDimension() );
		   Vector DeltaXi = ZeroVector( LocalSpaceDimension() );
		   CoordinatesArrayType CurrentGlobalCoords( ZeroVector(3) );
		   //Newton iteration:
		   double tol = 1.0e-8;
		   int maxiter = 1000;
		   for( int k=0; k<maxiter; k++ )
		   {
			   CurrentGlobalCoords = ZeroVector(3);
			   GlobalCoordinates( CurrentGlobalCoords, rResult );
			   noalias( CurrentGlobalCoords ) = rPoint - CurrentGlobalCoords;
			   InverseOfJacobian( J, rResult );
			   noalias(DeltaXi) = prod(J,CurrentGlobalCoords);
			   noalias( rResult ) += DeltaXi;
			   if( MathUtils<double>::Norm3(DeltaXi) > 30 )
			   {
				   break;
			   }
			   if( MathUtils<double>::Norm3(DeltaXi) < tol )
			   {
				   break;
			   }
		   }
		   return( rResult );
	   }
	   
	   /**
		* Returns whether given arbitrary point is inside the Geometry and the respective
        * local point for the given global point
		*/
	   virtual bool IsInside( const CoordinatesArrayType& rPoint, CoordinatesArrayType& rResult )
	   {
           KRATOS_ERROR(std::logic_error,
                        "Calling base class IsInside method instead of derived class one. Please check the definition of derived class." , *this);
           return false;
       }
       
       ///@}
      ///@name Inquiry
      ///@{
      
      /** This method confirm you if this geometry has a specific
	  integration method or not. This method will be usefull to
	  control the geometry before intagrating using a specific
	  method. In Geometry class this method controls if the
	  integration points vector respecting to this method is empty
	  or not.

	  @return bool true if this integration method exist and false if this
	  method is not imeplemented for this geometry.
      */
	bool HasIntegrationMethod(IntegrationMethod ThisMethod) const
	{
	  return (mpGeometryData->HasIntegrationMethod(ThisMethod));
	}

	/**
	* @return default integration method 
	*/

	IntegrationMethod GetDefaultIntegrationMethod()
	{
		return mpGeometryData->DefaultIntegrationMethod();
	}
      /** This method is to know if this geometry is symmetric or
	  not.

	  @todo Making some method related to symmetry axis and more...

	  @return bool true if this geometry is symmetric and false if
	  it's not.
      */
      virtual bool IsSymmetric() const
	{
	  return false;
	}

      ///@}
      ///@name Edge
      ///@{

      /** This method gives you number of all edges of this
       * geometry. For example, for a hexahedron, this would be
       * 12
       * 
       * @return SizeType containes number of this geometry edges.
       * @see EdgesNumber()
       * @see Edges()
       * @see FacesNumber()
       * @see Faces()
       */
      // will be used by refinement algorithm, thus uncommented. janosch.
      virtual SizeType EdgesNumber() const
	{
	  KRATOS_ERROR(std::logic_error,
		       "Calling base class Edge method instead of drived class one. Please check the definition of derived class." , *this);

	  return SizeType();
	}

      /** This method gives you all edges of this geometry. This
	  method will gives you all the edges with one dimension less
	  than this geometry. for example a triangle would return
	  three lines as its edges or a tetrahedral would return four
	  triangle as its edges but won't return its six edge
	  lines by this method.

	  @return GeometriesArrayType containes this geometry edges.
	  @see EdgesNumber()
	  @see Edge()
      */
      // will be used by refinement algorithm, thus uncommented. janosch.
      virtual GeometriesArrayType Edges(void)
	{
	  KRATOS_ERROR(std::logic_error,
		       "Calling base class Edge method instead of drived class one. Please check the definition of derived class." , *this);

	  return GeometriesArrayType();
	}

      /** This method gives you an edge of this geometry which holds
	  given points. This method will gives you an edge with
	  dimension related to given points number. for example a
	  tetrahedral would return a triangle for given three points or
	  return an edge line for given two nodes by this method.

	  @return Geometry::Pointer to this geometry specific edge.
	  @see EdgesNumber()
	  @see Edges()
      */
      // Commented for possible change in Edge interface of geometry. Pooyan.
//       virtual Pointer Edge(const PointsArrayType& EdgePoints)
// 	{
// 	  KRATOS_ERROR(std::logic_error,
// 		       "Calling base class Edge method instead of drived class one. Please check the definition of derived class." , *this);

// 	}
    
    /**
     * Returns the number of faces of the current geometry.
     * This is only implemented for 3D geometries, since 2D geometries
     * only have edges but no faces
     * @see EdgesNumber
     * @see Edges
     * @see Faces
     */
    virtual SizeType FacesNumber() const
    {
        KRATOS_ERROR(std::logic_error,
                     "Calling base class Edge method instead of drived class one. Please check the definition of derived class." , *this);

        return SizeType();
    }
    
    /**
     * Returns all faces of the current geometry.
     * This is only implemented for 3D geometries, since 2D geometries
     * only have edges but no faces
     * @see EdgesNumber
     * @see Edges
     * @see FacesNumber
     */
    virtual GeometriesArrayType Faces(void)
    {
        KRATOS_ERROR(std::logic_error,
                     "Calling base class Face method instead of drived class one. Please check the definition of derived class." , *this);

        return GeometriesArrayType();
    }

      /** This method gives you an edge of this geometry related to
	  given index. The numbering order of each geometries edges is
	  depended to type of that geometry.

	  @return Geometry::Pointer to this geometry specific edge.
	  @see EdgesNumber()
	  @see Edges()
      */
      // Commented for possible change in Edge interface of geometry. Pooyan.
//       virtual Pointer Edge(IndexType EdgeIndex)
// 	{
// 	  KRATOS_ERROR(std::logic_error,
// 		       "Calling base class Edge method instead of drived class one. Please check the definition of derived class." , *this);

// 	}

      /** This method gives you normal edge of this geometry which holds
	  given points. 

	  @return NormalType which is normal to this geometry specific edge.
	  @see Edge()
      */
      // Commented for possible change in Edge interface of geometry. Pooyan.
//       virtual NormalType NormalEdge(const PointsArrayType& EdgePoints)
// 	{
// 	  KRATOS_ERROR(std::logic_error,
// 		       "Calling base class NormalEdge method instead of drived class one. Please check the definition of derived class." , *this);

// 	  return NormalType();
// 	}

      /** This method gives you normal to edge of this geometry related to
	  given index. The numbering order of each geometries edges is
	  depended to type of that geometry.

	  @return NormalType which is normal to this geometry specific edge.
	  @see Edge()
      */
      // Commented for possible change in Edge interface of geometry. Pooyan.
//       virtual NormalType NormalEdge(IndexType EdgeIndex)
// 	{
// 	  KRATOS_ERROR(std::logic_error,
// 		       "Calling base class NormalEdge method instead of drived class one. Please check the definition of derived class." , *this);

// 	  return NormalType();
// 	}

      ///@}      
      ///@name Integration Points
      ///@{

      /** Number of integtation points for default integration
	  method. This method just call IntegrationPointsNumber(enum
	  IntegrationMethod ThisMethod) with default integration
	  method.

	  @return SizeType which is the number of integration points
	  for default integrating method.
      */
      SizeType IntegrationPointsNumber() const
	{
	  return mpGeometryData->IntegrationPoints().size();
	}

      /** Number of integtation points for given integration
	  method. This method use integration points data base to
	  obtain size of the integration points Vector respected to
	  given method.

	  @return SizeType which is the number of integration points
	  for given integrating method.
      */
      SizeType IntegrationPointsNumber(IntegrationMethod ThisMethod) const
	{
	  return mpGeometryData->IntegrationPointsNumber(ThisMethod);
	}


      /** Integtation points for default integration
	  method. This method just call IntegrationPoints(enum
	  IntegrationMethod ThisMethod) with default integration
	  method.

	  @return const IntegrationPointsArrayType which is Vector of integration points
	  for default integrating method.
      */
      const IntegrationPointsArrayType& IntegrationPoints() const
	{
	  return mpGeometryData->IntegrationPoints();
	}

      /** Integtation points for given integration
	  method. This method use integration points data base to
	  obtain integration points Vector respected to
	  given method.

	  @return const IntegrationPointsArrayType which is Vector of integration points
	  for default integrating method.
      */
      const IntegrationPointsArrayType& IntegrationPoints(IntegrationMethod ThisMethod) const
	{
	  return mpGeometryData->IntegrationPoints(ThisMethod);
	}

      ///@}      
      ///@name Jacobian 
      ///@{

	virtual CoordinatesArrayType& GlobalCoordinates(CoordinatesArrayType& rResult, CoordinatesArrayType const& LocalCoordinates)
	{
		noalias(rResult)= ZeroVector(3);
		for(IndexType i = 0 ; i < this->size() ; i++)
			noalias(rResult) += ShapeFunctionValue(i,LocalCoordinates) * GetPoint(i);
	  	return rResult;
	}

      /** Jacobians for default integration method. This method just
	  call Jacobian(enum IntegrationMethod ThisMethod) with
	  default integration method.

	  @return JacobiansType a Vector of jacobian
	  matrices \f$ J_i \f$ where \f$ i=1,2,...,n \f$ is the integration
	  point index of default integration method.

	  @see DeterminantOfJacobian
	  @see InverseOfJacobian
      */
      JacobiansType& Jacobian(JacobiansType& rResult) const
	{
	  Jacobian(rResult, mpGeometryData->DefaultIntegrationMethod());
	  return rResult;
	}

      /** Jacobians for given  method. This method
	  calculate jacobians matrices in all integrations points of
	  given integration method.

	  @param ThisMethod integration method which jacobians has to
	  be calculated in its integration points.

	  @return JacobiansType a Vector of jacobian
	  matrices \f$ J_i \f$ where \f$ i=1,2,...,n \f$ is the integration
	  point index of given integration method.

	  @see DeterminantOfJacobian
	  @see InverseOfJacobian
      */
      virtual JacobiansType& Jacobian(JacobiansType& rResult, IntegrationMethod ThisMethod) const
	{
	  KRATOS_ERROR(std::logic_error,
		       "Calling base class Jacobian method instead of drived class one. Please check the definition of derived class." , *this);
	  return rResult;
	}

      /** Jacobian in specific integration point of default integration method. This method just
	  call Jacobian(IndexType IntegrationPointIndex, enum IntegrationMethod ThisMethod) with
	  default integration method.

	  @param IntegrationPointIndex index of integration point which jacobians has to
	  be calculated in it.

	  @return Matrix<double> Jacobian matrix \f$ J_i \f$ where \f$
	  i \f$ is the given integration point index of default
	  integration method.

	  @see DeterminantOfJacobian
	  @see InverseOfJacobian
      */
      Matrix& Jacobian(Matrix& rResult, IndexType IntegrationPointIndex) const
	{
	  Jacobian(rResult, IntegrationPointIndex, mpGeometryData->DefaultIntegrationMethod());
	  return rResult;
	}

      /** Jacobian in specific integration point of given integration
	  method. This method calculate jacobian matrix in given
	  integration point of given integration method.

	  @param IntegrationPointIndex index of integration point which jacobians has to
	  be calculated in it.

	  @param ThisMethod integration method which jacobians has to
	  be calculated in its integration points.

	  @return Matrix<double> Jacobian matrix \f$ J_i \f$ where \f$
	  i \f$ is the given integration point index of given
	  integration method.

	  @see DeterminantOfJacobian
	  @see InverseOfJacobian
      */
      virtual Matrix& Jacobian(Matrix& rResult, IndexType IntegrationPointIndex, IntegrationMethod ThisMethod) const
	{
	  KRATOS_ERROR(std::logic_error,
		       "Calling base class Jacobian method instead of drived class one. Please check the definition of derived class." , *this);
	  return rResult;
	}

      /** Jacobian in given point. This method calculate jacobian
	  matrix in given point.

	  @param rPoint point which jacobians has to
	  be calculated in it.

	  @return Matrix of double which is jacobian matrix \f$ J \f$ in given point.

	  @see DeterminantOfJacobian
	  @see InverseOfJacobian
      */
      virtual Matrix& Jacobian(Matrix& rResult, const CoordinatesArrayType& rCoordinates) const
	{
	  KRATOS_ERROR(std::logic_error,
		       "Calling base class Jacobian method instead of drived class one. Please check the definition of derived class." , *this);
	  return rResult;
	}

      /** Determinant of jacobians for default integration method. This method just
	  call DeterminantOfJacobian(enum IntegrationMethod ThisMethod) with
	  default integration method.

	  @return Vector of double which is vector of determinants of
	  jacobians \f$ |J|_i \f$ where \f$ i=1,2,...,n \f$ is the
	  integration point index of default integration method.

	  @see Jacobian
	  @see InverseOfJacobian
      */
      Vector& DeterminantOfJacobian(Vector& rResult) const
	{
	  DeterminantOfJacobian(rResult, mpGeometryData->DefaultIntegrationMethod());
	  return rResult;
	}

      /** Determinant of jacobians for given integration method. This
	  method calculate determinant of jacobian in all
	  integrations points of given integration method.

	  @return Vector of double which is vector of determinants of
	  jacobians \f$ |J|_i \f$ where \f$ i=1,2,...,n \f$ is the
	  integration point index of given integration method.

	  @see Jacobian
	  @see InverseOfJacobian
      */
      virtual Vector& DeterminantOfJacobian(Vector& rResult, IntegrationMethod ThisMethod) const
	{
	  KRATOS_ERROR(std::logic_error,
		       "Calling base class DeterminantOfJacobian method instead of drived class one. Please check the definition of derived class." , *this);
	  return rResult;
	}

      /** Determinant of jacobian in specific integration point of
	  default integration method. This method just call
	  DeterminantOfJacobian(IndexType IntegrationPointIndex, enum
	  IntegrationMethod ThisMethod) with default integration
	  method.

	  @param IntegrationPointIndex index of integration point
	  which determinant jacobians has to be calculated in it.

	  @return Determinamt of jacobian matrix \f$ |J|_i \f$ where \f$
	  i \f$ is the given integration point index of default
	  integration method.

	  @see Jacobian
	  @see InverseOfJacobian
      */
      double DeterminantOfJacobian(IndexType IntegrationPointIndex) const
	{
	  return DeterminantOfJacobian(IntegrationPointIndex, mpGeometryData->DefaultIntegrationMethod());
	}

      /** Determinant of jacobian in specific integration point of
	  given integration method. This method calculate determinant
	  of jacobian in given integration point of given integration
	  method.

	  @param IntegrationPointIndex index of integration point which jacobians has to
	  be calculated in it.

	  @param IntegrationPointIndex index of integration point
	  which determinant of jacobians has to be calculated in it.

	  @return Determinamt of jacobian matrix \f$ |J|_i \f$ where \f$
	  i \f$ is the given integration point index of given
	  integration method.

	  @see Jacobian
	  @see InverseOfJacobian
      */
      virtual double DeterminantOfJacobian(IndexType IntegrationPointIndex, IntegrationMethod ThisMethod) const
	{
	  KRATOS_ERROR(std::logic_error,
		       "Calling base class DeterminantOfJacobian method instead of drived class one. Please check the definition of derived class." , *this);

	  return 0;
	}


      /** Determinant of jacobian in given point. This method calculate determinant of jacobian
	  matrix in given point.

	  @param rPoint point which determinant of jacobians has to
	  be calculated in it.

	  @return Determinamt of jacobian matrix \f$ |J| \f$ in given
	  point.

	  @see DeterminantOfJacobian
	  @see InverseOfJacobian
      */
      virtual double DeterminantOfJacobian(const Point<3>& rPoint) const
	{
	  KRATOS_ERROR(std::logic_error,
		       "Calling base class DeterminantOfJacobian method instead of drived class one. Please check the definition of derived class." , *this);

	  return 0;
	}


      /** Inverse of jacobians for default integration method. This method just
	  call InverseOfJacobian(enum IntegrationMethod ThisMethod) with
	  default integration method.

	  @return Inverse of jacobian
	  matrices \f$ J_i^{-1} \f$ where \f$ i=1,2,...,n \f$ is the integration
	  point index of default integration method.

	  @see Jacobian
	  @see DeterminantOfJacobian
      */
      JacobiansType& InverseOfJacobian(JacobiansType& rResult) const
	{
	  InverseOfJacobian(rResult, mpGeometryData->DefaultIntegrationMethod());
	  return rResult;
	}

      /** Inverse of jacobians for given integration method. This method
	  calculate inverse of jacobians matrices in all integrations points of
	  given integration method.

	  @param ThisMethod integration method which inverse of jacobians has to
	  be calculated in its integration points.

	  @return Inverse of jacobian
	  matrices \f$ J^{-1}_i \f$ where \f$ i=1,2,...,n \f$ is the integration
	  point index of given integration method.

	  @see Jacobian
	  @see DeterminantOfJacobian
      */
      virtual JacobiansType& InverseOfJacobian(JacobiansType& rResult, IntegrationMethod ThisMethod) const
	{
	  KRATOS_ERROR(std::logic_error,
		       "Calling base class InverseOfJacobian method instead of drived class one. Please check the definition of derived class." , *this);
	  return rResult;
	}

      /** Inverse of jacobian in specific integration point of default integration method. This method just
	  call InverseOfJacobian(IndexType IntegrationPointIndex, enum IntegrationMethod ThisMethod) with
	  default integration method.

	  @param IntegrationPointIndex index of integration point which inverse of jacobians has to
	  be calculated in it.

	  @return Inverse of jacobian matrix \f$ J^{-1}_i \f$ where \f$
	  i \f$ is the given integration point index of default
	  integration method.

	  @see Jacobian
	  @see DeterminantOfJacobian
      */
      Matrix& InverseOfJacobian(Matrix& rResult, IndexType IntegrationPointIndex) const
	{
	  InverseOfJacobian(rResult, IntegrationPointIndex, mpGeometryData->DefaultIntegrationMethod());
	  return rResult;
	}

      /** Inverse of jacobian in specific integration point of given integration
	  method. This method calculate Inverse of jacobian matrix in given
	  integration point of given integration method.

	  @param IntegrationPointIndex index of integration point which inverse of jacobians has to
	  be calculated in it.

	  @param ThisMethod integration method which inverse of jacobians has to
	  be calculated in its integration points.

	  @return Inverse of jacobian matrix \f$ J^{-1}_i \f$ where \f$
	  i \f$ is the given integration point index of given
	  integration method.

	  @see Jacobian
	  @see DeterminantOfJacobian
      */
      virtual Matrix& InverseOfJacobian(Matrix& rResult, IndexType IntegrationPointIndex, IntegrationMethod ThisMethod) const
	{
	  KRATOS_ERROR(std::logic_error,
		       "Calling base class InverseOfJacobian method instead of drived class one. Please check the definition of derived class." , *this);
	  return rResult;
	}

      /** Inverse of jacobian in given point. This method calculate inverse of jacobian
	  matrix in given point.

	  @param rPoint point which inverse of jacobians has to
	  be calculated in it.

	  @return Inverse of jacobian matrix \f$ J^{-1} \f$ in given point.

	  @see DeterminantOfJacobian
	  @see InverseOfJacobian
      */
      virtual Matrix& InverseOfJacobian(Matrix& rResult, const CoordinatesArrayType& rCoordinates) const
	{
	  KRATOS_ERROR(std::logic_error,
		       "Calling base class InverseOfJacobian method instead of drived class one. Please check the definition of derived class." , *this);
	  return rResult;
	}



      ///@}
      ///@name Shape Function
      ///@{ 

      /** This method gives all shape functions values evaluated in all
	  integration points of default integration method. It just
	  call ShapeFunctionsValues(enum IntegrationMethod ThisMethod)
	  with default integration method.There is no calculation and
	  it just give it from shape functions values container.

	  \note There is no control if the return matrix is empty or not!

	  @return Matrix of values of shape functions \f$ F_{ij} \f$
	  where i is the integration point index and j is the shape
	  function index. In other word component \f$ f_{ij} \f$ is value
	  of the shape function corresponding to node j evaluated in
	  integration point i of default integration method.

	  @see ShapeFunctionValue
	  @see ShapeFunctionsLocalGradients
	  @see ShapeFunctionLocalGradient
      */ 
      const Matrix& ShapeFunctionsValues() const
	{return mpGeometryData->ShapeFunctionsValues();} 

      /** This method gives all shape functions values evaluated in all
	  integration points of given integration method. There is no
	  calculation and it just give it from shape functions values
	  container. 

	  \note There is no control if the return matrix is empty or not!

	  @param ThisMethod integration method which shape functions
	  evaluated in its integration points.

	  @return Matrix of values of shape functions \f$ F_{ij} \f$
	  where i is the integration point index and j is the shape
	  function index. In other word component \f$ f_{ij} \f$ is value
	  of the shape function corresponding to node j evaluated in
	  integration point i of given integration method.

	  @see ShapeFunctionValue
	  @see ShapeFunctionsLocalGradients
	  @see ShapeFunctionLocalGradient
      */ 
      const Matrix& ShapeFunctionsValues(IntegrationMethod ThisMethod)  const
	{return mpGeometryData->ShapeFunctionsValues(ThisMethod);} 

      /** This method gives value of given shape function evaluated in
	  given integration point of default integration method. It just
	  call ShapeFunctionValue(IndexType IntegrationPointIndex,
	  IndexType ShapeFunctionIndex, enum IntegrationMethod
	  ThisMethod) with default integration method. There is no
	  calculation and it just give it from shape functions values
	  container if they are existing. Otherwise it gives you error
	  which this value is not exist.

	  @param IntegrationPointIndex index of integration point
	  which shape functions evaluated in it.

	  @param ShapeFunctionIndex index of node which correspounding
	  shape function evaluated in given integration point.

	  @return Value of given shape function in given integration
	  point of default integration method.

	  @see ShapeFunctionsValues
	  @see ShapeFunctionsLocalGradients
	  @see ShapeFunctionLocalGradient
      */ 
      double ShapeFunctionValue(IndexType IntegrationPointIndex, IndexType ShapeFunctionIndex) const
	{return mpGeometryData->ShapeFunctionValue(IntegrationPointIndex, ShapeFunctionIndex);} 

      /** This method gives value of given shape function evaluated in given
	  integration point of given integration method. There is no
	  calculation and it just give it from shape functions values
	  container if they are existing. Otherwise it gives you error
	  which this value is not exist.

	  @param IntegrationPointIndex index of integration point
	  which shape functions evaluated in it.

	  @param ShapeFunctionIndex index of node which correspounding
	  shape function evaluated in given integration point.

	  @param ThisMethod integration method which shape function
	  evaluated in its integration point.

	  @return Value of given shape function in given integration
	  point of given integration method.

	  @see ShapeFunctionsValues
	  @see ShapeFunctionsLocalGradients
	  @see ShapeFunctionLocalGradient
      */ 
      double ShapeFunctionValue(IndexType IntegrationPointIndex, IndexType ShapeFunctionIndex, IntegrationMethod ThisMethod) const
	{return mpGeometryData->ShapeFunctionValue(IntegrationPointIndex, ShapeFunctionIndex, ThisMethod);} 

      /** This method gives value of given shape function evaluated in given
	  point.

	  @param rPoint Point of evaluation of the shape
	  function. This point must be in local coordinate.

	  @param ShapeFunctionIndex index of node which correspounding
	  shape function evaluated in given integration point.

	  @return Value of given shape function in given point.

	  @see ShapeFunctionsValues
	  @see ShapeFunctionsLocalGradients
	  @see ShapeFunctionLocalGradient
      */ 
      virtual double ShapeFunctionValue(IndexType ShapeFunctionIndex, const CoordinatesArrayType& rCoordinates) const
	{
	  KRATOS_ERROR(std::logic_error,
		       "Calling base class ShapeFunctionValue method instead of drived class one. Please check the definition of derived class." , *this);

	  return 0;
	}

      /** This method gives all shape functions gradients evaluated in all
	  integration points of default integration method. It just
	  call ShapeFunctionsLocalGradients(enum IntegrationMethod ThisMethod)
	  with default integration method. There is no calculation and
	  it just give it from shape functions values container.

	  \note There is no control if there is any gradient calculated or not!

	  @return shape functions' gradients \f$ F_{ijk} \f$ where i
	  is the integration point index and j is the shape function
	  index and k is local coordinate index. In other word
	  component \f$ f_{ijk} \f$ is k'th component of gradient of
	  the shape function corresponding to node j evaluated in
	  integration point i of default integration method.

	  @see ShapeFunctionsValues
	  @see ShapeFunctionValue
	  @see ShapeFunctionLocalGradient
      */ 

      const ShapeFunctionsGradientsType& ShapeFunctionsLocalGradients() const
	{return mpGeometryData->ShapeFunctionsLocalGradients();} 

      /** This method gives all shape functions gradients evaluated in
	  all integration points of given integration method. There is
	  no calculation and it just give it from shape functions
	  values container.

	  \note There is no control if there is any gradient calculated or not!

	  @param ThisMethod integration method which shape functions
	  gradients evaluated in its integration points.

	  @return shape functions' gradients \f$ F_{ijk} \f$ where i
	  is the integration point index and j is the shape function
	  index and k is local coordinate index. In other word
	  component \f$ f_{ijk} \f$ is k'th component of gradient of
	  the shape function corresponding to node j evaluated in
	  integration point i of given integration method.

	  @see ShapeFunctionsValues
	  @see ShapeFunctionValue
	  @see ShapeFunctionLocalGradient
      */ 
      const ShapeFunctionsGradientsType& ShapeFunctionsLocalGradients(IntegrationMethod ThisMethod) const
	{return mpGeometryData->ShapeFunctionsLocalGradients(ThisMethod);} 

      /** This method gives gradient of given shape function evaluated in
	  given integration point of default integration method. It just
	  call ShapeFunctionLocalGradient(IndexType IntegrationPointIndex,
	  IndexType ShapeFunctionIndex, enum IntegrationMethod
	  ThisMethod) with default integration method. There is no
	  calculation and it just give it from shape functions values
	  container if they are existing. Otherwise it gives you error
	  which this value is not exist.

	  @param IntegrationPointIndex index of integration point
	  which shape function gradient evaluated in it.

	  @param ShapeFunctionIndex index of node which correspounding
	  shape function gradient evaluated in given integration point.

	  @return Gradient of given shape function in given integration
	  point of default integration method.

	  @see ShapeFunctionsValues
	  @see ShapeFunctionValue
	  @see ShapeFunctionsLocalGradients
      */ 
      Matrix& ShapeFunctionLocalGradient(IndexType IntegrationPointIndex)  const
	{return mpGeometryData->ShapeFunctionLocalGradient(IntegrationPointIndex);} 

      /** This method gives gradient of given shape function evaluated
	  in given integration point of given integration
	  method. There is no calculation and it just give it from
	  shape functions values container if they are
	  existing. Otherwise it gives you error which this value is
	  not exist.

	  @param IntegrationPointIndex index of integration point
	  which shape function gradient evaluated in it.

	  @param ShapeFunctionIndex index of node which correspounding
	  shape function gradient evaluated in given integration point.

	  @param ThisMethod integration method which shape function gradient
	  evaluated in its integration points.

	  @return Gradient of given shape function in given integration
	  point of given integration method.

	  @see ShapeFunctionsValues
	  @see ShapeFunctionValue
	  @see ShapeFunctionsLocalGradients
      */ 
 //     Matrix& ShapeFunctionLocalGradient(IndexType IntegrationPointIndex, IndexType ShapeFunctionIndex, IntegrationMethod ThisMethod)  const
//	{return mpGeometryData->ShapeFunctionLocalGradient(IntegrationPointIndex, ThisMethod);} 
      
      
      /** This method gives gradient of all shape functions evaluated
       * in given point.
       * There is no calculation and it just give it from
       * shape functions values container if they are
       * existing. Otherwise it gives you error which this value is
       * not exist.
       * 
       * @param rResult the given Container that will be overwritten by the solution
       * @param rPoint the given local coordinates the gradients will be evaluated for 
       * @return a matrix of gradients for each shape function 
       */
      virtual Matrix& ShapeFunctionsLocalGradients( Matrix& rResult, const CoordinatesArrayType& rPoint ) const 
      {
          KRATOS_ERROR(std::logic_error,
                       "Calling base class ShapeFunctionsLocalGradients method instead of derived class one. Please check the definition of derived class." , *this);
          return rResult; 
      }
      
      /** This method gives second order derivatives of all shape
       * functions evaluated in given point.
       * 
       * @param rResult the given container will be overwritten by the results
       * @param rPoint the given local coordinates the derivatives will be evaluated for.
       * @return a third order tensor containing the second order derivatives for each shape function 
       */
      virtual ShapeFunctionsSecondDerivativesType& ShapeFunctionsSecondDerivatives( ShapeFunctionsSecondDerivativesType& rResult, const CoordinatesArrayType& rPoint ) const 
      {
          KRATOS_ERROR(std::logic_error,
                       "Calling base class ShapeFunctionsSecondDerivatives method instead of derived class one. Please check the definition of derived class." , *this);
          return rResult; 
      }

      /** This method gives third order derivatives of all shape
       * functions evaluated in given point.
       * 
       * @param rResult the given container will be overwritten by the results
       * @param rPoint the given local coordinates the derivatives will be evaluated for.
       * @return a fourth order tensor containing the second order derivatives for each shape function 
       */
      virtual ShapeFunctionsThirdDerivativesType& ShapeFunctionsThirdDerivatives( ShapeFunctionsThirdDerivativesType& rResult, const CoordinatesArrayType& rPoint ) const 
      {
          KRATOS_ERROR(std::logic_error,
                       "Calling base class ShapeFunctionsThirdDerivatives method instead of derived class one. Please check the definition of derived class." , *this);
          return rResult; 
      }


     ShapeFunctionsGradientsType& ShapeFunctionsIntegrationPointsGradients(ShapeFunctionsGradientsType& rResult) const
	{
	  ShapeFunctionsIntegrationPointsGradients(rResult, mpGeometryData->DefaultIntegrationMethod());
	  return rResult;
	}

      virtual ShapeFunctionsGradientsType& ShapeFunctionsIntegrationPointsGradients(ShapeFunctionsGradientsType& rResult, IntegrationMethod ThisMethod) const
	{
	  KRATOS_ERROR(std::logic_error,
		       "Calling base class ShapeFunctionsGaussPointsGradients method instead of derived class one. Please check the definition of derived class." , *this);
	  return rResult;
	}


     boost::numeric::ublas::vector<Matrix> const& MassFactors() const
       {
	 return mpGeometryData->MassFactors();
       }
      
     boost::numeric::ublas::vector<Matrix> const& MassFactors(IntegrationMethod ThisMethod) const
       {
	 return  mpGeometryData->MassFactors(ThisMethod);
       }
      
     Matrix const& MassFactors(IndexType IntegrationPointIndex) const
       {
	 return mpGeometryData->MassFactors(IntegrationPointIndex);
       }
      
     Matrix const& MassFactors(IndexType IntegrationPointIndex, IntegrationMethod ThisMethod) const
       {
	 return mpGeometryData->MassFactors(IntegrationPointIndex, ThisMethod);
       }
      
      ///@}      
      ///@name Input and output
      ///@{

      /** Turn back information as a string.

      @return String contains information about this geometry.
      @see PrintData()
      @see PrintInfo()
      */
      virtual std::string Info() const
	{
	  std::stringstream buffer;
	  buffer << Dimension()
			 << " dimensional geometry in " 
			 << WorkingSpaceDimension()
			 << "D space";
	  return buffer.str();
	}

      /** Print information about this object.

      @param rOStream Stream to print into it.
      @see PrintData()
      @see Info()
      */
      virtual void PrintInfo(std::ostream& rOStream) const
	{
	  rOStream << Dimension()  << " dimensional geometry in " << WorkingSpaceDimension() << "D space";
	}

      /** Print geometry's data into given stream. Prints it's points
	  by the order they stored in the geometry and then center
	  point of geometry.

	  @param rOStream Stream to print into it.
	  @see PrintInfo()
	  @see Info()
      */
      virtual void PrintData(std::ostream& rOStream) const
	{
	  if(mpGeometryData)
	    mpGeometryData->PrintData(rOStream);
	  rOStream << std::endl;
	  rOStream << std::endl;
	  for(unsigned int i = 0 ; i < this->size() ; ++i)
	    {
	      rOStream << "    Point " << i+1 << "\t : ";
	      (*this)[i].PrintData(rOStream);
	      rOStream << std::endl;
	    }

	  rOStream << "    Center\t : ";
	  Center().PrintData(rOStream);
	  rOStream << std::endl;
	  rOStream << std::endl;
	  rOStream << "    Length\t : " << Length() << std::endl;
	  rOStream << "    Area\t : " << Area() << std::endl;
	  rOStream << "    Volume\t : " << Volume();
// 	  for(unsigned int i = 0 ; i < mPoints.size() ; ++i)
// 	    {
// 	      rOStream << "    Point " << i+1 << "\t            : ";
// 	      mPoints[i].PrintData(rOStream);
// 	      rOStream << std::endl;
// 	    }

// 	  rOStream << "    Center\t            : ";
// 	  Center().PrintData(rOStream);
// 	  rOStream << std::endl;
// 	  rOStream << std::endl;
// 	  rOStream << "    Length                  : " << Length() << std::endl;
// 	  rOStream << "    Area                    : " << Area() << std::endl;
// 	  rOStream << "    Volume                  : " << Volume();
	}      
      
            
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
      
      /** Protected Constructor.
	  Avoids object to be created Except for derived classes   
      */

            
      ///@}
  
    

     private:
      ///@name Static Member Variables 
      ///@{ 
        
// 	  static const GeometryData msEmptyGeometryData;
        
      ///@} 
      ///@name Member Variables 
      ///@{ 
        
	  GeometryData const* mpGeometryData;
         
           

      ///@} 
      ///@name Private Operators
      ///@{ 
        
        
      ///@} 
      ///@name Private Operations
      ///@{ 

	  static const GeometryData GenerateEmptyGeometryData()
	    {
	      IntegrationPointsContainerType integration_points = {};
	      ShapeFunctionsValuesContainerType shape_functions_values = {};	  
	      ShapeFunctionsLocalGradientsContainerType shape_functions_local_gradients = {};	  
	      return GeometryData(2,
							  2,
							  2,
							  GeometryData::GI_GAUSS_1,
							  integration_points,
							  shape_functions_values,
							  shape_functions_local_gradients);
	    }
        
        
      ///@} 
      ///@name Private  Access 
      ///@{ 
        
        
      ///@}    
      ///@name Private Inquiry 
      ///@{ 
        
        
      ///@}    
      ///@name Private Friends
      ///@{
      
      template<class TOtherPointType> friend class Geometry;

      ///@}    
      ///@name Un accessible methods 
      ///@{ 
      
        
      ///@}    
        
    }; // Class Geometry 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
      /// input stream function
      template<class TPointType>
	inline std::istream& operator >> (std::istream& rIStream, 
					  Geometry<TPointType>& rThis);

      /// output stream function
      template<class TPointType>
	inline std::ostream& operator << (std::ostream& rOStream, 
					  const Geometry<TPointType>& rThis)
	{
	  rThis.PrintInfo(rOStream);
	  rOStream << std::endl;
	  rThis.PrintData(rOStream);

	  return rOStream;
	}
      ///@} 

//        template<class TPointType>  
// /*  	const GeometryData Geometry<TPointType>::msEmptyGeometryData = GeometryData(TPointType::Dimension(), TPointType::Dimension(), TPointType::Dimension());  */
//  	const GeometryData Geometry<TPointType>::msEmptyGeometryData = Geometry<TPointType>::GenerateEmptyGeometryData(); 

  
}  // namespace Kratos.

#endif // KRATOS_GEOMETRY_H_INCLUDED  defined 


