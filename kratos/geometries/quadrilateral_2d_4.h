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
//   Last Modified by:    $Author: nelson $
//   Date:                $Date: 2009-01-21 09:56:09 $
//   Revision:            $Revision: 1.15 $
//
//
#if !defined(KRATOS_QUADRILATERAL_2D_4_H_INCLUDED )
#define  KRATOS_QUADRILATERAL_2D_4_H_INCLUDED



// System includes


// External includes
#include <boost/array.hpp> 


// Project includes
#include "includes/define.h"
#include "geometries/geometry.h"
#include "geometries/line_2d_2.h"
#include "integration/quadrature.h"
#include "integration/quadrilateral_gaussian_integration_points.h"


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
  
    /**
     * A four node quadrilateral geometry. While the shape functions are only defined in
     * 2D it is possible to define an arbitrary orientation in space. Thus it can be used for
     * defining surfaces on 3D elements.
     */
    template<class TPointType> class Quadrilateral2D4 
    : public Geometry<TPointType>
    {
        public:
        ///@}
        ///@name Type Definitions
        ///@{
      
        /**
         * Geometry as base class.
         */ 
        typedef Geometry<TPointType> BaseType; 
        typedef BoundingBox<TPointType, BaseType>  BoundingBoxType;  


        /**
         * Type of edge geometry
         */
        typedef Line2D2<TPointType> EdgeType;
                 
        /**
         * Pointer definition of Quadrilateral2D4
         */ 
        KRATOS_CLASS_POINTER_DEFINITION(Quadrilateral2D4);
        
        /** 
         * Integration methods implemented in geometry. 
         */
        typedef GeometryData::IntegrationMethod IntegrationMethod;
        
        /** 
         * A Vector of counted pointers to Geometries.  
         * Used for returning edges of the geometry.
         */
        typedef typename BaseType::GeometriesArrayType GeometriesArrayType;
      
        /**  
         * Redefinition of template parameter TPointType.
         */
        typedef TPointType PointType;
        
        /** 
         * Type used for indexing in geometry class.
         * std::size_t used for indexing
         * point or integration point access methods and also all other
         * methods which need point or integration point index.
         */
        typedef typename BaseType::IndexType IndexType;
        
        /** 
         * This type is used to return size or dimension in
         * geometry. Dimension, WorkingDimension, PointsNumber and
         * ... return this type as their results.
         */
        typedef typename BaseType::SizeType SizeType;
      
        /** 
         * Array of counted pointers to point. 
         * This type used to hold geometry's points.
         */
        typedef  typename BaseType::PointsArrayType PointsArrayType;
        
        /**
         * Array of coordinates. Can be Nodes, Points or IntegrationPoints
         */
        typedef typename BaseType::CoordinatesArrayType CoordinatesArrayType;
      
        /** 
         * This type used for representing an integration point in geometry. 
         * This integration point is a point with an additional weight component.
         */
        typedef typename BaseType::IntegrationPointType IntegrationPointType;
      
        /** 
         * A Vector of IntegrationPointType which used to hold
         * integration points related to an integration
         * method. 
         * IntegrationPoints functions used this type to return
         * their results.
         */
        typedef typename BaseType::IntegrationPointsArrayType IntegrationPointsArrayType;
      
        /** 
         * A Vector of IntegrationPointsArrayType which used to hold
         * integration points related to different integration method
         * implemented in geometry.
         */
        typedef typename BaseType::IntegrationPointsContainerType IntegrationPointsContainerType;
      
        /** 
         * A third order tensor used as shape functions' values
         * container.
         */
        typedef typename BaseType::ShapeFunctionsValuesContainerType ShapeFunctionsValuesContainerType;
      
        /** 
         * A fourth order tensor used as shape functions' local
         * gradients container in geometry.
         */
        typedef typename BaseType::ShapeFunctionsLocalGradientsContainerType ShapeFunctionsLocalGradientsContainerType;
      
        /** 
         * A third order tensor to hold jacobian matrices evaluated at
         * integration points. Jacobian and InverseOfJacobian functions
         * return this type as their result.
         */
        typedef typename BaseType::JacobiansType JacobiansType;
      
        /** 
         * A third order tensor to hold shape functions' local
         * gradients. ShapefunctionsLocalGradients function return this
         * type as its result.
         */
        typedef typename BaseType::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;
        
        /** 
         * A third order tensor to hold shape functions' local second derivatives.
         * ShapefunctionsLocalGradients function return this
         * type as its result.
         */
        typedef typename BaseType::ShapeFunctionsSecondDerivativesType
                ShapeFunctionsSecondDerivativesType;
      
        /** 
        * A third order tensor to hold shape functions' local second derivatives.
        * ShapefunctionsLocalGradients function return this
        * type as its result.
        */
        typedef typename BaseType::ShapeFunctionsThirdDerivativesType
                ShapeFunctionsThirdDerivativesType;
      
        /** 
         * Type of the normal vector used for normal to edges in geomety.
         */
        typedef typename BaseType::NormalType NormalType;
      	
        ///@}
        ///@name Life Cycle 
        ///@{ 

        Quadrilateral2D4( const PointType& FirstPoint, 
                        const PointType& SecondPoint, 
                        const PointType& ThirdPoint, 
                        const PointType& FourthPoint ) 
            : BaseType(PointsArrayType(), &msGeometryData)
        {
            this->Points().push_back(typename PointType::Pointer(new PointType(FirstPoint)));
            this->Points().push_back(typename PointType::Pointer(new PointType(SecondPoint)));
            this->Points().push_back(typename PointType::Pointer(new PointType(ThirdPoint)));
            this->Points().push_back(typename PointType::Pointer(new PointType(FourthPoint)));
        }
        
        Quadrilateral2D4( typename PointType::Pointer pFirstPoint, 
                          typename PointType::Pointer pSecondPoint, 
                          typename PointType::Pointer pThirdPoint, 
                          typename PointType::Pointer pFourthPoint) 
            : BaseType(PointsArrayType(), &msGeometryData)
        {
            this->Points().push_back(pFirstPoint);
            this->Points().push_back(pSecondPoint);
            this->Points().push_back(pThirdPoint);
            this->Points().push_back(pFourthPoint);
        }

	Quadrilateral2D4( const PointsArrayType& ThisPoints ) 
            : BaseType(ThisPoints, &msGeometryData)
        {
            if( this->PointsNumber() != 4)
                KRATOS_ERROR(std::invalid_argument,
                             "Invalid points number. Expected 4, given " , this->PointsNumber());
        }
        
        /** 
         * Copy constructor.
         * Construct this geometry as a copy of given geometry.
         * 
         * @note This copy constructor does not copy the points and new
         * geometry shares points with given source geometry. It is
         * obvious that any change to this new geometry's point affect
         * source geometry's points too.
         */
        Quadrilateral2D4(Quadrilateral2D4 const& rOther)
            : BaseType(rOther)
        {
        }
        
        /** 
         * Copy constructor from a geometry with other point type.
         * Construct this geometry as a copy of given geometry which
         * has different type of points. The given goemetry's
         * TOtherPointType* must be implicity convertible to this
         * geometry PointType.
         * 
         * @note This copy constructor does not copy the points and new
         * geometry shares points with given source geometry. It is
         * obvious that any change to this new geometry's point affect
         * source geometry's points too.
         */
        template<class TOtherPointType> Quadrilateral2D4( Quadrilateral2D4<TOtherPointType> const& rOther )
            : BaseType(rOther)
        {
	}
 
        /**
         * Destructor. Does nothing!!!
         */
        virtual ~Quadrilateral2D4(){}
        
	    GeometryData::KratosGeometryFamily GetGeometryFamily(){return GeometryData::Kratos_Quadrilateral; }
	    GeometryData::KratosGeometryType GetGeometryType(){return GeometryData::Kratos_Quadrilateral2D4; }

        ///@}
        ///@name Operators 
        ///@{
      
        /** 
         * Assignment operator.
         * 
         * @note This operator don't copy the points and this
         * geometry shares points with given source geometry. It's
         * obvious that any change to this geometry's point affect
         * source geometry's points too.
         * 
         * @see Clone
         * @see ClonePoints
         */
        Quadrilateral2D4& operator=(const Quadrilateral2D4& rOther)
        {
            BaseType::operator=(rOther);
            return *this;
        }
        
        /** 
         * Assignment operator for geometries with different point type.
         * 
         * @note This operator don't copy the points and this
         * geometry shares points with given source geometry. It's
         * obvious that any change to this geometry's point affect
         * source geometry's points too.
         * 
         * @see Clone
         * @see ClonePoints
         */
        template<class TOtherPointType> 
                Quadrilateral2D4& operator=(Quadrilateral2D4<TOtherPointType> const & rOther)
        {
            BaseType::operator=(rOther);
            return *this;
        }
        
        ///@}
        ///@name Operations
        ///@{
      
        typename BaseType::Pointer Create( PointsArrayType const& ThisPoints ) const
        {
            return typename BaseType::Pointer( new Quadrilateral2D4(ThisPoints) );
	}
 
        virtual boost::shared_ptr< Geometry< Point<3> > > Clone() const
        {
            Geometry< Point<3> >::PointsArrayType NewPoints;
            
            //making a copy of the nodes TO POINTS (not Nodes!!!)
            for(IndexType i = 0 ; i < this->Points().size() ; i++)
                NewPoints.push_back(this->Points()[i]);
            
            //creating a geometry with the new points
            boost::shared_ptr< Geometry< Point<3> > > p_clone(new Quadrilateral2D4< Point<3> >(NewPoints));
            p_clone->ClonePoints();
            
            return p_clone;
        }
        
        /**
         * returns the local coordinates of all nodes of the current geometry
         * @param rResult a Matrix object that will be overwritten by the result
         * @return the local coordinates of all nodes
         */
        virtual Matrix& PointsLocalCoordinates( Matrix& rResult ) const
        {
            noalias(rResult) = ZeroMatrix(4,2);
            rResult(0,0) = -1.0;
            rResult(0,1) = -1.0;
            rResult(1,0) =  1.0;
            rResult(1,1) = -1.0;
            rResult(2,0) =  1.0;
            rResult(2,1) =  1.0;
            rResult(3,0) = -1.0;
            rResult(3,1) =  1.0;
            return rResult;
        }
        
        //lumping factors for the calculation of the lumped mass matrix
        virtual Vector& LumpingFactors(Vector& rResult) const
        { 
            rResult.resize(4, false);
            std::fill(rResult.begin(), rResult.end(), 1.00 / 4.00);
            return rResult;
        }
        
        ///@}
        ///@name Information
        ///@{
         
        /** 
         * This method calculates and returns Length or charactereistic
         * length of this geometry depending on it's dimension. 
         * For one dimensional geometry for example Line it returns 
         * length of it and for the other geometries it gives Characteristic 
         * length otherwise.
         * In the current geometry this function returns the determinant of
         * jacobian
         * 
         * @return double value contains length or Characteristic
         * length
         * @see Area()
         * @see Volume()
         * @see DomainSize()
         */
        /**
         * :TODO: could be replaced by something more suitable
         * (comment by janosch)
         */
        virtual double Length() const
	{										
	//return sqrt(fabs( DeterminantOfJacobian(PointType())));
           double length = 0.000; 
           length = sqrt(fabs(Area()));
           return length;  

        }
        
        /** This method calculates and returns area or surface area of
         * this geometry depending to it's dimension. For one dimensional
         * geometry it returns zero, for two dimensional it gives area
         * and for three dimensional geometries it gives surface area.
         * 
         * @return double value contains area or surface
         * area.N
         * @see Length()
         * @see Volume()
         * @see DomainSize()
         */
        /**
         * :TODO: could be replaced by something more suitable
         * (comment by janosch)
         */
        virtual double Area() const
        {				

	  Vector temp;
	  DeterminantOfJacobian(temp, msGeometryData.DefaultIntegrationMethod());
	  const IntegrationPointsArrayType& integration_points = this->IntegrationPoints(msGeometryData.DefaultIntegrationMethod());
          double Area=0.00;
	  for(unsigned int i=0;i<integration_points.size();i++)
		{
			Area+= temp[i]*integration_points[i].Weight();
		}
          //KRATOS_WATCH(temp)
		return Area;
         
            //return fabs(DeterminantOfJacobian(PointType())) * 0.5;                                                       
        }
        
        
        /** This method calculates and returns length, area or volume of
         * this geometry depending to it's dimension. For one dimensional
         * geometry it returns its length, for two dimensional it gives area
         * and for three dimensional geometries it gives its volume.
         * 
         * @return double value contains length, area or volume.
         * @see Length()
         * @see Area()
         * @see Volume()
         */
        /**
         * :TODO: could be replaced by something more suitable
         * (comment by janosch)
         */
        virtual double DomainSize() const
        {
            return fabs(DeterminantOfJacobian(PointType())) * 0.5;
        }
		
		/**
		 * Returns whether given arbitrary point is inside the Geometry
		 */
		virtual bool IsInside( const CoordinatesArrayType& rPoint, CoordinatesArrayType& rResult )
		{
			this->PointLocalCoordinates( rResult, rPoint );
			if( rResult[0] >= -1.0 && rResult[0] <= 1.0 )
				if( rResult[1] >= -1.0 && rResult[1] <= 1.0 )
					return true;
			return false;
		}
        


        /** 
           Geometry as derivate class class.      

        */
/*
      virtual typename  BoundingBoxType::Pointer Create
               (const TPointType& LowPoint, const TPointType& HighPoint) 
	  {
                //mpBondingBoxType = new BoundingBoxType(LowPoint, HighPoint, this ) ;
                //return mpBondingBoxType;  
		return typename BoundingBoxType::Pointer(new BoundingBoxType(LowPoint, HighPoint, this)) ;    
	  }  
*/

       //typedef Quadrilateral2D4<TPointType> Quadrilateral2D4Type;    
       virtual void Bounding_Box(BoundingBoxType& rResult) const
             {              
                //KRATOS_WATCH("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")  
                //rResult.GetGeom() = this; 
                BaseType::Bounding_Box(rResult.LowPoint(), rResult.HighPoint());  
             }
     



        ///@}      
        ///@name Jacobian
        ///@{
        /**
         *	TODO: implemented but not yet tested
         */
        /** 
         * Jacobians for given method. 
         * This method calculates jacobians matrices in all
         * integrations points of given integration method.
         * 
         * @param ThisMethod integration method which jacobians has to
         * be calculated in its integration points.
         * 
         * @return JacobiansType a Vector of jacobian
         * matrices \f$ J_i \f$ where \f$ i=1,2,...,n \f$ is the integration
         * point index of given integration method.
         * 
         * @see DeterminantOfJacobian
         * @see InverseOfJacobian
         */
        virtual JacobiansType& Jacobian( JacobiansType& rResult, 
                                         IntegrationMethod ThisMethod) const
        {
            //getting derivatives of shape functions
            ShapeFunctionsGradientsType shape_functions_gradients =
                    CalculateShapeFunctionsIntegrationPointsLocalGradients(ThisMethod);
            //getting values of shape functions
            Matrix shape_functions_values =
                    CalculateShapeFunctionsIntegrationPointsValues(ThisMethod);
            //workaround by riccardo...
            if(rResult.size() != this->IntegrationPointsNumber(ThisMethod))
            {
                // KLUDGE: While there is a bug in ublas 
                // vector resize, I have to put this beside resizing!!
                JacobiansType temp(this->IntegrationPointsNumber(ThisMethod));
                rResult.swap(temp);
            }
            //loop over all integration points
            for( unsigned int pnt=0; pnt < this->IntegrationPointsNumber(ThisMethod); pnt++ )
            {
                //defining single jacobian matrix
                Matrix jacobian = ZeroMatrix(2,2);
                //loop over all nodes
                for( unsigned int i=0; i<this->PointsNumber(); i++ )
                {
                    jacobian(0,0) += 
                            (this->GetPoint(i).X())*(shape_functions_gradients[pnt](i,0));
                    jacobian(0,1) += 
                            (this->GetPoint(i).X())*(shape_functions_gradients[pnt](i,1));
                    jacobian(1,0) += 
                            (this->GetPoint(i).Y())*(shape_functions_gradients[pnt](i,0));
                    jacobian(1,1) += 
                            (this->GetPoint(i).Y())*(shape_functions_gradients[pnt](i,1));
                }
                rResult[pnt] = jacobian;
            }//end of loop over all integration points
            return rResult;
        }
        
        /**
         *	TODO: implemented but not yet tested
         */
        /** 
         * Jacobian in specific integration point of given integration
         * method. This method calculate jacobian matrix in given
         * integration point of given integration method.
         * 
         * @param IntegrationPointIndex index of integration point which jacobians has to
         * be calculated in it.
         * 
         * @param ThisMethod integration method which jacobians has to
         * be calculated in its integration points.
         * 
         * @return Matrix<double> Jacobian matrix \f$ J_i \f$ where \f$
         * i \f$ is the given integration point index of given
         * integration method.
         * 
         * @see DeterminantOfJacobian
         * @see InverseOfJacobian
         */
        virtual Matrix& Jacobian( Matrix& rResult, 
                                  IndexType IntegrationPointIndex, 
                                  IntegrationMethod ThisMethod) const
        {
            //setting up size of jacobian matrix
            rResult.resize(2,2);
            //derivatives of shape functions
            ShapeFunctionsGradientsType shape_functions_gradients =
                    CalculateShapeFunctionsIntegrationPointsLocalGradients(ThisMethod);
            Matrix ShapeFunctionsGradientInIntegrationPoint =
                    shape_functions_gradients(IntegrationPointIndex);
            //values of shape functions in integration points
            vector<double> ShapeFunctionValuesInIntegrationPoint = ZeroVector(4);
            /*vector<double>*/ ShapeFunctionValuesInIntegrationPoint = 
                    row(CalculateShapeFunctionsIntegrationPointsValues( 
					ThisMethod ),IntegrationPointIndex);
            
            //Elements of jacobian matrix (e.g. J(1,1) = dX1/dXi1)
            //loop over all nodes
            for( unsigned int i=0; i<this->PointsNumber(); i++ )
            {
                rResult(0,0) +=
                        (this->GetPoint(i).X())*(ShapeFunctionsGradientInIntegrationPoint(i,0));
                rResult(0,1) +=
                        (this->GetPoint(i).X())*(ShapeFunctionsGradientInIntegrationPoint(i,1));
                rResult(1,0) +=
                        (this->GetPoint(i).Y())*(ShapeFunctionsGradientInIntegrationPoint(i,0));
                rResult(1,1) +=
                        (this->GetPoint(i).Y())*(ShapeFunctionsGradientInIntegrationPoint(i,1));
            }
            return rResult;
        }
        
        /**
         *	TODO: implemented but not yet tested
         */
     	/**
         * Jacobian in given point. This method calculate jacobian
         * matrix in given point.
         * 
         * @param rPoint point which jacobians has to
	 *	be calculated in it.
	 *
	 *	@return Matrix of double which is jacobian matrix \f$ J \f$ in given point.
	 *
	 *	@see DeterminantOfJacobian
	 *	@see InverseOfJacobian
     	 */
     	virtual Matrix& Jacobian(Matrix& rResult, const CoordinatesArrayType& rPoint) const
	{
            //setting up size of jacobian matrix
            rResult.resize(2,2);
            //derivatives of shape functions
            Matrix shape_functions_gradients;
            shape_functions_gradients = ShapeFunctionsLocalGradients( 
                    shape_functions_gradients, rPoint );
            //Elements of jacobian matrix (e.g. J(1,1) = dX1/dXi1)
            //loop over all nodes
            for( unsigned int i=0; i<this->PointsNumber(); i++ )
            {
                rResult(0,0) += (this->GetPoint(i).X())*(shape_functions_gradients(i,0));
                rResult(0,1) += (this->GetPoint(i).X())*(shape_functions_gradients(i,1));
                rResult(1,0) += (this->GetPoint(i).Y())*(shape_functions_gradients(i,0));
                rResult(1,1) += (this->GetPoint(i).Y())*(shape_functions_gradients(i,1));
            }
            return rResult;
        }
        
        /**
         *	TODO: implemented but not yet tested
         */
        /** 
         * Determinant of jacobians for given integration method. 
         * This method calculates determinant of jacobian in all
         * integrations points of given integration method.
         * 
         * @return Vector of double which is vector of determinants of
         * jacobians \f$ |J|_i \f$ where \f$ i=1,2,...,n \f$ is the
         * integration point index of given integration method.
         * 
         * @see Jacobian
         * @see InverseOfJacobian
         */
        virtual Vector& DeterminantOfJacobian( Vector& rResult, 
                                               IntegrationMethod ThisMethod) const
        {
            //workaround by riccardo
            if(rResult.size() != this->IntegrationPointsNumber(ThisMethod))
            {
																rResult.resize(this->IntegrationPointsNumber(ThisMethod), false);
                // KLUDGE: While there is a bug in ublas 
                // vector resize, I have to put this beside resizing!!
//                Vector temp(this->IntegrationPointsNumber(ThisMethod));
//                rResult.swap(temp);
            }
												
            //for all integration points
            for( unsigned int pnt=0; pnt<this->IntegrationPointsNumber(ThisMethod); pnt++ )
            {
                rResult[pnt] = DeterminantOfJacobian( pnt, ThisMethod );
            }
            return rResult;
        }
        
        /**
         *	TODO: implemented but not yet tested
         */
        /** 
         * Determinant of jacobian in specific integration point of
         * given integration method. This method calculate determinant
         * of jacobian in given integration point of given integration
         * method.
         * 
         * @param IntegrationPointIndex index of integration point which jacobians has to
         * be calculated in it.
         * 
         * @param IntegrationPointIndex index of integration point
         * which determinant of jacobians has to be calculated in it.
         * 
         * @return Determinamt of jacobian matrix \f$ |J|_i \f$ where \f$
         * i \f$ is the given integration point index of given
         * integration method.
         * 
         * @see Jacobian
         * @see InverseOfJacobian
         */
        virtual double DeterminantOfJacobian( IndexType IntegrationPointIndex, 
                                              IntegrationMethod ThisMethod) const
        {
            Matrix jacobian = ZeroMatrix(2,2);
            jacobian = Jacobian( jacobian, IntegrationPointIndex, ThisMethod );
            return( (jacobian(0,0)*jacobian(1,1))-(jacobian(0,1)*jacobian(1,0)) );
        }
        
        /** 
         *	TODO: implemented but not yet tested
         */
        /** 
         * Determinant of jacobian in given point. 
         * This method calculate determinant of jacobian
         * matrix in given point.
         * @param rPoint point which determinant of jacobians has to
         * be calculated in it.
         * 
         * @return Determinamt of jacobian matrix \f$ |J| \f$ in given
         * point.
         * 
         * @see DeterminantOfJacobian
         * @see InverseOfJacobian
         * 
         * KLUDGE: PointType needed for proper functionality
         * KLUDGE: works only with explicitly generated Matrix object
         */
        /**
         * :TODO: needs to be changed to Point<3> again. As PointType can
         * be a Node with unique ID or an IntegrationPoint or any arbitrary
         * point in space this needs to be reviewed
         * (comment by janosch)
         */
        virtual double DeterminantOfJacobian( const CoordinatesArrayType& rPoint ) const
        {
            Matrix jacobian = ZeroMatrix(2,2);
            jacobian = Jacobian( jacobian, rPoint );
            return( (jacobian(0,0)*jacobian(1,1))-(jacobian(0,1)*jacobian(1,0)));
															
	}
 
        /**
         *	TODO: implemented but not yet tested
         */
        /** 
         * Inverse of jacobians for given integration method. 
         * This method calculates inverse of jacobians matrices 
         * in all integrations points of
         * given integration method.
         * 
         * @param ThisMethod integration method which inverse of jacobians has to
         * be calculated in its integration points.
         * 
         * @return Inverse of jacobian
         * matrices \f$ J^{-1}_i \f$ where \f$ i=1,2,...,n \f$ is the integration
         * point index of given integration method.
         * 
         * @see Jacobian
         * @see DeterminantOfJacobian
         * 
         * KLUDGE: works only with explicitly generated Matrix object
         */
        virtual JacobiansType& InverseOfJacobian( JacobiansType& rResult, 
                IntegrationMethod ThisMethod ) const
        {
            //workaround by riccardo
            if(rResult.size() != this->IntegrationPointsNumber(ThisMethod))
            {
                // KLUDGE: While there is a bug in ublas
                // vector resize, I have to put this beside resizing!!
                JacobiansType temp(this->IntegrationPointsNumber(ThisMethod));
                rResult.swap(temp);
            }
            //loop over all integration points
            for( unsigned int pnt=0; pnt<this->IntegrationPointsNumber(ThisMethod); pnt++ )
            {
                Matrix tempMatrix(2,2);
                rResult[pnt] = InverseOfJacobian( tempMatrix, pnt, ThisMethod );
            }
            return rResult;
        }
        
        /**
         *	TODO: implemented but not yet tested
         */
        /** 
         * Inverse of jacobian in specific integration point of given integration
         * method. This method calculate Inverse of jacobian matrix in given
         * integration point of given integration method.
         * 
         * @param IntegrationPointIndex index of integration point
         * which inverse of jacobians has to
         * be calculated in it.
         * @param ThisMethod integration method which inverse of jacobians has to
         * be calculated in its integration points.
         * 
         * @return Inverse of jacobian matrix \f$ J^{-1}_i \f$ where \f$
         * i \f$ is the given integration point index of given
         * integration method.
         * 
         * @see Jacobian
         * @see DeterminantOfJacobian
         * 
         * KLUDGE: works only with explicitly generated Matrix object
         */
        virtual Matrix& InverseOfJacobian( Matrix& rResult, 
                                           IndexType IntegrationPointIndex, 
                                           IntegrationMethod ThisMethod ) const
        {
            //current jacobian
            Matrix tempMatrix = ZeroMatrix(2,2);
            tempMatrix = Jacobian( tempMatrix, IntegrationPointIndex, ThisMethod );
            //determinant of jacobian
            double det_j = DeterminantOfJacobian( IntegrationPointIndex, ThisMethod );
            //checking for singularity
            if(det_j == 0.00)
                KRATOS_ERROR(std::runtime_error, 
                             "Zero determinant of jacobian during inversion of matrix!" , 
                             *this);
            //setting up result matrix
            rResult.resize(2,2);
            //filling matrix
            rResult(0,0) = (tempMatrix(1,1))/(det_j);
            rResult(1,0) = -(tempMatrix(1,0))/(det_j);
            rResult(0,1) = -(tempMatrix(0,1))/(det_j);
            rResult(1,1) = (tempMatrix(0,0))/(det_j);
            return rResult;
        }
        
        /**
         *	TODO: implemented but not yet tested
         */
        /** 
         * Inverse of jacobian in given point.
         * This method calculates inverse of jacobian
         * matrix in given point.
         * @param rPoint point which inverse of jacobians has to
         * be calculated in it.
         * @return Inverse of jacobian matrix \f$ J^{-1} \f$ in given point.
         * 
         * @see DeterminantOfJacobian
         * @see InverseOfJacobian
         * 
         * KLUDGE: works only with explicitly generated Matrix object
         */
        virtual Matrix& InverseOfJacobian( Matrix& rResult, 
                                           const CoordinatesArrayType& rPoint ) const
        {
            //current jacobian
            Matrix tempMatrix = ZeroMatrix(2,2);
            tempMatrix = Jacobian( tempMatrix, rPoint );
            //deteminant of Jacobian
            double det_j = DeterminantOfJacobian( rPoint );
            //checking for singularity
            if(det_j == 0.00)
                KRATOS_ERROR(std::runtime_error, 
                             "Zero determinant of jacobian during inversion of matrix!", 
                             *this);
            //setting up result matrix
            rResult.resize(2,2);
            //filling matrix
            rResult(0,0) = (tempMatrix(1,1))/(det_j);
            rResult(1,0) = -(tempMatrix(1,0))/(det_j);
            rResult(0,1) = -(tempMatrix(0,1))/(det_j);
            rResult(1,1) = (tempMatrix(0,0))/(det_j);
            return rResult;
        }
        
        /** This method gives you number of all edges of this
        geometry. This method will gives you number of all the edges
        with one dimension less than this geometry. for example a
        triangle would return three or a tetrahedral would return
        four but won't return nine related to its six edge lines.

        @return SizeType containes number of this geometry edges.
        @see Edges()
        @see Edge()
        */
        virtual SizeType EdgesNumber() const
        {
            return 4;
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
        virtual GeometriesArrayType Edges(void)
        {
            GeometriesArrayType edges = GeometriesArrayType();
            edges.push_back( EdgeType( this->pGetPoint(0), this->pGetPoint(1) ) );
            edges.push_back( EdgeType( this->pGetPoint(1), this->pGetPoint(2) ) );
            edges.push_back( EdgeType( this->pGetPoint(2), this->pGetPoint(3) ) );
            edges.push_back( EdgeType( this->pGetPoint(3), this->pGetPoint(0) ) );
            return edges;
        }
        
        ///@}
        ///@name Shape Function
        ///@{ 
        
        /**
         *	TODO: implemented but not yet tested
         */
        /**
         * Calculates the value of a given shape function at a given point.
         * 
         * @param ShapeFunctionIndex The number of the desired shape function
         * @param rPoint the given point in local coordinates at which the value of the shape
         * function is calculated
         * 
         * @return the value of the shape function at the given point 
         */
     	virtual double ShapeFunctionValue( IndexType ShapeFunctionIndex, 
                                         const CoordinatesArrayType& rPoint) const
        {
            switch( ShapeFunctionIndex )
            {
                case 0:
                    return( 0.25*(1.0-rPoint[0])*(1.0-rPoint[1]) );
                case 1:
                    return( 0.25*(1.0+rPoint[0])*(1.0-rPoint[1]) );
                case 2:
                    return( 0.25*(1.0+rPoint[0])*(1.0+rPoint[1]) );
                case 3:
                    return( 0.25*(1.0-rPoint[0])*(1.0+rPoint[1]) );
                default:
                    KRATOS_ERROR(std::logic_error,
                                 "Wrong index of shape function!" , 
                                 *this);
            }
            return 0;
        }
        
        /**
         *	TODO: implemented but not yet tested
         */
        /**
         * Calculates the Gradients of the shape functions.
         * Calculates the gradients of the shape functions with regard to 
         * the global coordinates in all
         * integration points (\f$ \frac{\partial N^i}{\partial X_j} \f$)
         * 
         * @param rResult a container which takes the calculated gradients
         * @param ThisMethod the given IntegrationMethod
         * 
         * @return the gradients of all shape functions with regard to the global coordinates
         * KLUDGE: method call only works with explicit JacobiansType rather than creating
         * JacobiansType within argument list
	 */
        virtual ShapeFunctionsGradientsType& ShapeFunctionsIntegrationPointsGradients(
                ShapeFunctionsGradientsType& rResult, 
                IntegrationMethod ThisMethod) const
        {
            const unsigned int integration_points_number =
                    msGeometryData.IntegrationPointsNumber(ThisMethod);
            if(integration_points_number == 0)
                KRATOS_ERROR(std::logic_error,
                             "This integration method is not supported" , *this);
            //workaround by riccardo
            if(rResult.size() != integration_points_number)
            {
                // KLUDGE: While there is a bug in ublas
                // vector resize, I have to put this beside resizing!!
                ShapeFunctionsGradientsType temp(integration_points_number);
                rResult.swap(temp);
            }
            //calculating the local gradients
            ShapeFunctionsGradientsType locG =
                    CalculateShapeFunctionsIntegrationPointsLocalGradients(ThisMethod);
            //getting the inverse jacobian matrices
            JacobiansType temp( integration_points_number );
            JacobiansType invJ = InverseOfJacobian( temp, ThisMethod );
            //loop over all integration points
            for( unsigned int pnt=0; pnt<integration_points_number; pnt++ )
            {
                rResult[pnt].resize(4,2);
                for( int i=0; i<4; i++ )
                {
                    for( int j=0; j<2; j++ )
                    {
                        rResult[pnt](i,j) = 
                                (locG[pnt](i,0)*invJ[pnt](j,0))
                                +(locG[pnt](i,1)*invJ[pnt](j,1));
                    }
                } 
            }//end of loop over integration points
            return rResult;
        }
        
        ///@}
        ///@name Input and output
        ///@{
        
        /** 
         * Turn back information as a string.
         * 
         * @return String contains information about this geometry.
         * @see PrintData()
         * @see PrintInfo()
         */
        virtual std::string Info() const
        {
            return "2 dimensional quadrilateral with four nodes in 2D space";
        }
        
        /** 
         * Print information about this object.
         * @param rOStream Stream to print into it.
         * @see PrintData()
         * @see Info()
         */
        virtual void PrintInfo(std::ostream& rOStream) const
        {
            rOStream << "2 dimensional quadrilateral with four nodes in 2D space";
        }
        
        /** 
         * Print geometry's data into given stream. 
         * Prints it's points
         * by the order they stored in the geometry and then center
         * point of geometry.
         * 
         * @param rOStream Stream to print into it.
         * @see PrintInfo()
         * @see Info()
         */
        /**
         * :TODO: needs to be reviewed because it is not properly implemented yet
         * (comment by janosch)
         */
        virtual void PrintData(std::ostream& rOStream) const
        {
            BaseType::PrintData(rOStream);
            std::cout << std::endl;
            Matrix jacobian;
            Jacobian(jacobian, PointType());
            rOStream << "    Jacobian in the origin\t : " << jacobian;
        }
        
        /**
         * Calculates the local gradients for all integration points for
         * given integration method
         */
        virtual ShapeFunctionsGradientsType ShapeFunctionsLocalGradients( 
                IntegrationMethod ThisMethod)
        {
            ShapeFunctionsGradientsType localGradients 
                    = CalculateShapeFunctionsIntegrationPointsLocalGradients( ThisMethod );
            const int integration_points_number 
                    = msGeometryData.IntegrationPointsNumber(ThisMethod);
            ShapeFunctionsGradientsType Result(integration_points_number);
            for( int pnt=0; pnt<integration_points_number; pnt++ )
            {
                Result[pnt] = localGradients[pnt];
            }
            return Result;
        }
        
        /**
         * Calculates the local gradients for all integration points for the
         * default integration method
         */
        virtual ShapeFunctionsGradientsType ShapeFunctionsLocalGradients()
        {
            IntegrationMethod ThisMethod = msGeometryData.DefaultIntegrationMethod();
            ShapeFunctionsGradientsType localGradients 
                    = CalculateShapeFunctionsIntegrationPointsLocalGradients( ThisMethod );
            const int integration_points_number 
                    = msGeometryData.IntegrationPointsNumber(ThisMethod);
            ShapeFunctionsGradientsType Result(integration_points_number);
            for( int pnt=0; pnt<integration_points_number; pnt++ )
            {
                Result[pnt] = localGradients[pnt];
            }
            return Result;
        }
        
        /**
         * Calculates the gradients in terms of local coordinates 
         * of all shape functions in a given point.
         * 
         * @param rPoint the current point at which the gradients are calculated in local
         * coordinates
         * @return the gradients of all shape functions 
         * \f$ \frac{\partial N^i}{\partial \xi_j} \f$
         */
        virtual Matrix& ShapeFunctionsLocalGradients( Matrix& rResult, 
                const CoordinatesArrayType& rPoint ) const 
        {
            rResult.resize(4,2);
            noalias(rResult) = ZeroMatrix(4,2);
            rResult(0,0) = -0.25*(1.0-rPoint[1]);
            rResult(0,1) = -0.25*(1.0-rPoint[0]);
            rResult(1,0) = 0.25*(1.0-rPoint[1]);
            rResult(1,1) = -0.25*(1.0+rPoint[0]);
            rResult(2,0) = 0.25*(1.0+rPoint[1]);
            rResult(2,1) = 0.25*(1.0+rPoint[0]);
            rResult(3,0) = -0.25*(1.0+rPoint[1]);
            rResult(3,1) = 0.25*(1.0-rPoint[0]);
            return rResult;
        }
        
        
        
        /**
         * returns the shape function gradients in an arbitrary point, 
         * given in local coordinates
         * 
         * @param rResult the matrix of gradients, 
         * will be overwritten with the gradients for all
         * shape functions in given point
         * @param rPoint the given point the gradients are calculated in
         */
        virtual Matrix& ShapeFunctionsGradients( Matrix& rResult, PointType& rPoint )
        {
            rResult.resize(4,2);
            rResult(0,0) = -0.25*(1.0-rPoint.Y());
            rResult(0,1) = -0.25*(1.0-rPoint.X());
            rResult(1,0) = 0.25*(1.0-rPoint.Y());
            rResult(1,1) = -0.25*(1.0+rPoint.X());
            rResult(2,0) = 0.25*(1.0+rPoint.Y());
            rResult(2,1) = 0.25*(1.0+rPoint.X());
            rResult(3,0) = -0.25*(1.0+rPoint.Y());
            rResult(3,1) = 0.25*(1.0-rPoint.X());
            return rResult;
        }
        
        /**
         * returns the second order derivatives of all shape functions 
         * in given arbitrary pointers
         * @param rResult a third order tensor which contains the second derivatives
         * @param rPoint the given point the second order derivatives are calculated in
         */
        virtual ShapeFunctionsSecondDerivativesType& ShapeFunctionsSecondDerivatives( ShapeFunctionsSecondDerivativesType& rResult, const CoordinatesArrayType& rPoint ) const
        {
            if(rResult.size() != this->PointsNumber())
            {
                // KLUDGE: While there is a bug in
                // ublas vector resize, I have to put this beside resizing!!
                ShapeFunctionsGradientsType temp(this->PointsNumber());
                rResult.swap(temp);
            }
            rResult[0].resize( 2, 2 );
            rResult[1].resize( 2, 2 );
            rResult[2].resize( 2, 2 );
            rResult[3].resize( 2, 2 );
            rResult[0](0,0) = 0.0;
            rResult[0](0,1) = 0.25;
            rResult[0](1,0) = 0.25;
            rResult[0](1,1) = 0.0;
            rResult[1](0,0) = 0.0;
            rResult[1](0,1) = -0.25;
            rResult[1](1,0) = -0.25;
            rResult[1](1,1) = 0.0;
            rResult[2](0,0) = 0.0;
            rResult[2](0,1) = 0.25;
            rResult[2](1,0) = 0.25;
            rResult[2](1,1) = 0.0;
            rResult[3](0,0) = 0.0;
            rResult[3](0,1) = -0.25;
            rResult[3](1,0) = -0.25;
            rResult[3](1,1) = 0.0;
            return rResult;
        }
        /**
         * returns the third order derivatives of all shape functions 
         * in given arbitrary pointers
         * @param rResult a fourth order tensor which contains the third derivatives
         * @param rPoint the given point the third order derivatives are calculated in
         */
        virtual ShapeFunctionsThirdDerivativesType& ShapeFunctionsThirdDerivatives( ShapeFunctionsThirdDerivativesType& rResult, const CoordinatesArrayType& rPoint ) const
        {
           if(rResult.size() != this->PointsNumber())
            {
                // KLUDGE: While there is a bug in
                // ublas vector resize, I have to put this beside resizing!!
//                 ShapeFunctionsGradientsType 
                ShapeFunctionsThirdDerivativesType temp(this->PointsNumber());
                rResult.swap(temp);
            }
            for( IndexType i=0; i<rResult.size(); i++ )
            {
                vector<Matrix> temp(this->PointsNumber());
                rResult[i].swap(temp);
            }
			for( unsigned int i=0; i<this->PointsNumber(); i++ )
			{
				for( unsigned int j=0; j<2; j++)
				{
					rResult[i][j].resize(2,2);
					noalias(rResult[i][j]) = ZeroMatrix(2,2);
				}
			}
//             rResult(0,0).resize( 2, 2);
//             rResult(0,1).resize( 2, 2);
//             rResult(1,0).resize( 2, 2);
//             rResult(1,1).resize( 2, 2);
//             rResult(2,0).resize( 2, 2);
//             rResult(2,1).resize( 2, 2);
//             rResult(3,0).resize( 2, 2);
//             rResult(3,0).resize( 2, 2);
            
            for(int i=0; i<4; i++)
            {
				rResult[i][0](0,0)=0.0;
				rResult[i][0](0,1)=0.0;
				rResult[i][0](1,0)=0.0;
				rResult[i][0](1,1)=0.0;
				rResult[i][1](0,0)=0.0;
				rResult[i][1](0,1)=0.0;
				rResult[i][1](1,0)=0.0;
				rResult[i][1](1,1)=0.0;
            }
            
            return rResult;
        }
        
        ///@}
        ///@name Friends
        ///@{
            
        ///@}

        protected:
        /**
         * There are no protected members in class Quadrilateral2D4
         */
        
        private:
        ///@name Static Member Variables 
        ///@{ 
        static const GeometryData msGeometryData;
        
        ///@} 
        ///@name Member Variables 
        ///@{ 
        

        ///@} 
        ///@name Private Operators
        ///@{ 
        
        
        ///@} 
        ///@name Private Operations
        ///@{
        
        /**
         * TODO: implemented but not yet tested
         */
        /**
         * Calculates the values of all shape function in all integration points.
         * Integration points are expected to be given in local coordinates
         * @param ThisMethod the current integration method
         * @return the matrix of values of every shape function in each integration point
         * :KLUDGE: number of points is hard-coded -> be careful if you want to copy and paste!
         */
        static Matrix CalculateShapeFunctionsIntegrationPointsValues( 
                typename BaseType::IntegrationMethod ThisMethod)
        {
            IntegrationPointsContainerType all_integration_points = 
                    AllIntegrationPoints();
            IntegrationPointsArrayType integration_points = 
                    all_integration_points[ThisMethod];
            //number of integration points
            const int integration_points_number = integration_points.size();
            //number of nodes in current geometry
            const int points_number = 4;
            //setting up return matrix
            Matrix shape_function_values(integration_points_number, points_number );
            //loop over all integration points
            for(int pnt = 0; pnt < integration_points_number; pnt++)
            {
                shape_function_values(pnt,0) =
                        0.25*(1.0-integration_points[pnt].X())
                        *(1.0-integration_points[pnt].Y());
                shape_function_values(pnt,1) =
                        0.25*(1.0+integration_points[pnt].X())
                        *(1.0-integration_points[pnt].Y());
                shape_function_values(pnt,2) =
                        0.25*(1.0+integration_points[pnt].X())
                        *(1.0+integration_points[pnt].Y());
                shape_function_values(pnt,3) =
                        0.25*(1.0-integration_points[pnt].X())
                        *(1.0+integration_points[pnt].Y());
            }
            return shape_function_values;
        }
        
        /**
         *	TODO: implemented but not yet tested
         */
        /**
         * Calculates the local gradients of all shape functions 
         * in all integration points.
         * Integration points are expected to be given in local coordinates
         *
         * @param ThisMethod the current integration method
         * @return the vector of the gradients of all shape functions 
         * in each integration point
         */
        static ShapeFunctionsGradientsType
                CalculateShapeFunctionsIntegrationPointsLocalGradients( 
                typename BaseType::IntegrationMethod ThisMethod)
        {
            IntegrationPointsContainerType all_integration_points = 
                    AllIntegrationPoints();
            IntegrationPointsArrayType integration_points = 
                    all_integration_points[ThisMethod];
            //number of integration points
            const int integration_points_number = integration_points.size();
            ShapeFunctionsGradientsType d_shape_f_values(integration_points_number);
            //initialising container
            //std::fill(d_shape_f_values.begin(), d_shape_f_values.end(), Matrix(4,2));
            //loop over all integration points
            for( int pnt=0; pnt<integration_points_number; pnt++ )
            {
                Matrix result(4,2);
                result(0,0) = -0.25*(1.0-integration_points[pnt].Y());
                result(0,1) = -0.25*(1.0-integration_points[pnt].X());
                result(1,0) = 0.25*(1.0-integration_points[pnt].Y());
                result(1,1) = -0.25*(1.0+integration_points[pnt].X());
                result(2,0) = 0.25*(1.0+integration_points[pnt].Y());
                result(2,1) = 0.25*(1.0+integration_points[pnt].X());
                result(3,0) = -0.25*(1.0+integration_points[pnt].Y());
                result(3,1) = 0.25*(1.0-integration_points[pnt].X());
                d_shape_f_values[pnt] = result;
            }
            return d_shape_f_values;
        }
        
        /**
         *	TODO: testing
         */
        static const IntegrationPointsContainerType AllIntegrationPoints()
        {
            IntegrationPointsContainerType integration_points = 
			{{
                			Quadrature< QuadrilateralGaussianIntegrationPoints1, 
                            2, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                            Quadrature<QuadrilateralGaussianIntegrationPoints2, 
                            2, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                            Quadrature<QuadrilateralGaussianIntegrationPoints3,
                            2, IntegrationPoint<3> >::GenerateIntegrationPoints()
			}};
            return integration_points;
        }
        
        /**
         *	TODO: testing
         */
        static const ShapeFunctionsValuesContainerType AllShapeFunctionsValues()
		{
            ShapeFunctionsValuesContainerType shape_functions_values = 
			{{
                Quadrilateral2D4<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                        GeometryData::GI_GAUSS_1),
                Quadrilateral2D4<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                        GeometryData::GI_GAUSS_2),
                Quadrilateral2D4<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                        GeometryData::GI_GAUSS_3)
			}};
            return shape_functions_values;
        }
        
        /**
         *	TODO: testing
         */
        static const ShapeFunctionsLocalGradientsContainerType 
                AllShapeFunctionsLocalGradients()
        {
            ShapeFunctionsLocalGradientsContainerType shape_functions_local_gradients = 
			{{
                Quadrilateral2D4<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(GeometryData::GI_GAUSS_1),
                Quadrilateral2D4<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(GeometryData::GI_GAUSS_2),
                Quadrilateral2D4<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(GeometryData::GI_GAUSS_3)
			}};
            return shape_functions_local_gradients;
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
      
        template<class TOtherPointType> friend class Quadrilateral2D4;
        
        ///@}
        ///@name Un accessible methods 
        ///@{ 
      
        Quadrilateral2D4();

            
        
        ///@}
    }; // Class Geometry 

    ///@}
    ///@name Type Definitions
    ///@{ 
  
  
    ///@} 
    ///@name Input and output 
    ///@{ 
    /**
     * input stream functions
     */
    template<class TPointType> inline std::istream& operator >> ( 
                     std::istream& rIStream, 
                     Quadrilateral2D4<TPointType>& rThis );
    /**
     * output stream functions
     */
    template<class TPointType> inline std::ostream& operator << ( 
            std::ostream& rOStream, 
            const Quadrilateral2D4<TPointType>& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);
        return rOStream;
    }
    
    ///@} 
    
    template<class TPointType> const 
            GeometryData Quadrilateral2D4<TPointType>::msGeometryData( 
                2, 2, 2,
                GeometryData::GI_GAUSS_2,
                Quadrilateral2D4<TPointType>::AllIntegrationPoints(),
                Quadrilateral2D4<TPointType>::AllShapeFunctionsValues(),
                AllShapeFunctionsLocalGradients()
                                                                     );
}// namespace Kratos.

#endif // KRATOS_QUADRILATERAL_2D_4_H_INCLUDED  defined 

