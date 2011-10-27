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
//   Date:                $Date: 2009-01-22 18:44:36 $
//   Revision:            $Revision: 1.4 $
//
//


#if !defined(KRATOS_LINE_2D_3_H_INCLUDED )
#define  KRATOS_LINE_2D_3_H_INCLUDED



// System includes


// External includes
#include <boost/array.hpp>


// Project includes
#include "includes/define.h"
#include "geometries/geometry.h"
#include "integration/quadrature.h"
#include "integration/gauss_legendre_integration_points.h"


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
    */
    template<class TPointType>

    class Line2D3 : public Geometry<TPointType>
    {
        public:
            ///@}
            ///@name Type Definitions
            ///@{

            /// Geometry as base class.
            typedef Geometry<TPointType> BaseType;

            /// Pointer definition of Line2D3
            KRATOS_CLASS_POINTER_DEFINITION( Line2D3 );

            /** Integration methods implemented in geometry.
            */
            typedef GeometryData::IntegrationMethod IntegrationMethod;

            /** A Vector of counted pointers to Geometries. Used for
            returning edges of the geometry.
             */
            typedef typename BaseType::GeometriesArrayType GeometriesArrayType;

            /** Redefinition of template parameter TPointType.
             */
            typedef TPointType PointType;

            /** Type used for indexing in geometry class.std::size_t used for indexing
            point or integration point access methods and also all other
            methods which need point or integration point index.
            */
            typedef typename BaseType::IndexType IndexType;


            /** This typed used to return size or dimension in
            geometry. Dimension, WorkingDimension, PointsNumber and
            ... return this type as their results.
            */
            typedef typename BaseType::SizeType SizeType;

            /** Array of counted pointers to point. This type used to hold
            geometry's points.
            */
            typedef  typename BaseType::PointsArrayType PointsArrayType;

            /** This type used for representing an integration point in
            geometry. This integration point is a point with an
            additional weight component.
            */
            typedef typename BaseType::IntegrationPointType IntegrationPointType;

            /** A Vector of IntegrationPointType which used to hold
            integration points related to an integration
            method. IntegrationPoints functions used this type to return
            their results.
            */
            typedef typename BaseType::IntegrationPointsArrayType IntegrationPointsArrayType;

            /** A Vector of IntegrationPointsArrayType which used to hold
            integration points related to different integration method
            implemented in geometry.
            */
            typedef typename BaseType::IntegrationPointsContainerType IntegrationPointsContainerType;

            /** A third order tensor used as shape functions' values
            continer.
            */
            typedef typename BaseType::ShapeFunctionsValuesContainerType ShapeFunctionsValuesContainerType;

            /** A fourth order tensor used as shape functions' local
            gradients container in geometry.
            */
            typedef typename BaseType::ShapeFunctionsLocalGradientsContainerType ShapeFunctionsLocalGradientsContainerType;

            /** A third order tensor to hold jacobian matrices evaluated at
            integration points. Jacobian and InverseOfJacobian functions
            return this type as their result.
            */
            typedef typename BaseType::JacobiansType JacobiansType;

            /** A third order tensor to hold shape functions' local
            gradients. ShapefunctionsLocalGradients function return this
            type as its result.
            */
            typedef typename BaseType::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;

            /** Type of the normal vector used for normal to edges in geomety.
             */
            typedef typename BaseType::NormalType NormalType;

            /**
             * Type of coordinates array
             */
            typedef typename BaseType::CoordinatesArrayType CoordinatesArrayType;

            ///@}
            ///@name Life Cycle
            ///@{

            Line2D3( const PointType& FirstPoint, const PointType& SecondPoint, const PointType& ThirdPoint )
                    : BaseType( PointsArrayType(), &msGeometryData )
            {
                BaseType::Points().push_back( typename PointType::Pointer( new PointType( FirstPoint ) ) );
                BaseType::Points().push_back( typename PointType::Pointer( new PointType( SecondPoint ) ) );
                BaseType::Points().push_back( typename PointType::Pointer( new PointType( ThirdPoint ) ) );
            }

            Line2D3( typename PointType::Pointer pFirstPoint, typename PointType::Pointer pSecondPoint,
                     typename PointType::Pointer pThirdPoint )
                    : BaseType( PointsArrayType(), &msGeometryData )
            {
                BaseType::Points().push_back( pFirstPoint );
                BaseType::Points().push_back( pSecondPoint );
                BaseType::Points().push_back( pThirdPoint );
            }

            Line2D3( const PointsArrayType& ThisPoints )
                    : BaseType( ThisPoints, &msGeometryData )
            {
                if ( BaseType::PointsNumber() != 3 )
                    KRATOS_ERROR( std::invalid_argument,
                                  "Invalid points number. Expected 3, given " , BaseType::PointsNumber() );
            }

            /** Copy constructor.
             * Construct this geometry as a copy of given geometry.
             * @note This copy constructor don't copy the points and new
             * geometry shares points with given source geometry. It's
             * obvious that any change to this new geometry's point affect
             * source geometry's points too.
             */
            Line2D3( Line2D3 const& rOther )
                    : BaseType( rOther )
            {
            }

            /** Copy constructor from a geometry with other point type.
             * Construct this geometry as a copy of given geometry which
             * has different type of points. The given goemetry's
             * TOtherPointType* must be implicity convertible to this
             * geometry PointType.
             * @note This copy constructor don't copy the points and new
             * geometry shares points with given source geometry. It's
             * obvious that any change to this new geometry's point affect
             * source geometry's points too.
             */
            template<class TOtherPointType> Line2D3( Line2D3<TOtherPointType> const& rOther )
                    : BaseType( rOther )
            {
            }

            /// Destructor. Do nothing!!!
            virtual ~Line2D3() {}

            GeometryData::KratosGeometryFamily GetGeometryFamily() {return GeometryData::Kratos_Linear; }

            GeometryData::KratosGeometryType GetGeometryType() {return GeometryData::Kratos_Line2D3; }

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
            Line2D3& operator=( const Line2D3& rOther )
            {
                BaseType::operator=( rOther );
                return *this;
            }

            /** Assignment operator for geometries with different point type.
             * @note This operator don't copy the points and this
             * geometry shares points with given source geometry. It's
             * obvious that any change to this geometry's point affect
             * source geometry's points too.
             * @see Clone
             * @see ClonePoints
             */
            template<class TOtherPointType>
            Line2D3& operator=( Line2D3<TOtherPointType> const & rOther )
            {
                BaseType::operator=( rOther );
                return *this;
            }

            ///@}
            ///@name Operations
            ///@{

            typename BaseType::Pointer Create( PointsArrayType const& ThisPoints ) const
            {
                return typename BaseType::Pointer( new Line2D3( ThisPoints ) );
            }

            virtual boost::shared_ptr< Geometry< Point<3> > > Clone() const
            {
                Geometry< Point<3> >::PointsArrayType NewPoints;
                //making a copy of the nodes TO POINTS (not Nodes!!!)

                for ( IndexType i = 0 ; i < BaseType::Points().size() ; i++ )
                    NewPoints.push_back( BaseType::Points()[i] );

                //creating a geometry with the new points
                boost::shared_ptr< Geometry< Point<3> > > p_clone( new Line2D3< Point<3> >( NewPoints ) );

                p_clone->ClonePoints();

                return p_clone;
            }

            //lumping factors for the calculation of the lumped mass matrix
            virtual Vector& LumpingFactors( Vector& rResult ) const
            {
                rResult.resize( 3, false );
                rResult[0] = 0.25;
                rResult[1] = 0.5;
                rResult[2] = 0.25;
                return rResult;
            }

            ///@}
            ///@name Informations
            ///@{

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
                double lenght = pow( BaseType::GetPoint( 0 ).X() - BaseType::GetPoint( 2 ).X(), 2 ) + pow( BaseType::GetPoint( 0 ).Y() - BaseType::GetPoint( 2 ).Y(), 2 );
                return sqrt( lenght );
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
                return 0.00;
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
                double lenght = pow( BaseType::GetPoint( 0 ).X() - BaseType::GetPoint( 2 ).X(), 2 ) + pow( BaseType::GetPoint( 0 ).Y() - BaseType::GetPoint( 2 ).Y(), 2 );
                return sqrt( lenght );
            }



            /**
             * Returns whether given arbitrary point is inside the Geometry
             */
            virtual bool IsInside( const CoordinatesArrayType& rPoint, CoordinatesArrayType& rResult )
            {
                this->PointLocalCoordinates( rResult, rPoint );

                if ( fabs( rResult[0] ) < 1 + 1.0e-8 )
                    return true;

                return false;
            }

            ///@}
            ///@name Jacobian
            ///@{

            /** Jacobians for given  method. This method
             * calculate jacobians matrices in all integrations points of
             * given integration method.
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
            virtual JacobiansType& Jacobian( JacobiansType& rResult, IntegrationMethod ThisMethod ) const
            {
                //getting derivatives of shape functions
                ShapeFunctionsGradientsType shape_functions_gradients =
                    CalculateShapeFunctionsIntegrationPointsLocalGradients( ThisMethod );
                //getting values of shape functions
                Matrix shape_functions_values =
                    CalculateShapeFunctionsIntegrationPointsValues( ThisMethod );

                if ( rResult.size() != this->IntegrationPointsNumber( ThisMethod ) )
                {
                    JacobiansType temp( this->IntegrationPointsNumber( ThisMethod ) );
                    rResult.swap( temp );
                }

                //loop over all integration points
                for ( unsigned  int pnt = 0; pnt < this->IntegrationPointsNumber( ThisMethod ); pnt++ )
                {
                    //defining single jacobian matrix
                    Matrix jacobian = ZeroMatrix( 2, 1 );
                    //loop over all nodes

                    for ( unsigned  int i = 0; i < this->PointsNumber(); i++ )
                    {
                        jacobian( 0, 0 ) += ( this->GetPoint( i ).X() ) * ( shape_functions_gradients[pnt]( i, 0 ) );
                        jacobian( 1, 0 ) += ( this->GetPoint( i ).Y() ) * ( shape_functions_gradients[pnt]( i, 0 ) );
                    }

                    rResult[pnt] = jacobian;
                }//end of loop over all integration points

                return rResult;
            }

            /** Jacobian in specific integration point of given integration
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
             * @see DeterminantOfJacobian
             * @see InverseOfJacobian
             */
            virtual Matrix& Jacobian( Matrix& rResult, IndexType IntegrationPointIndex, IntegrationMethod ThisMethod ) const
            {
                //setting up size of jacobian matrix
                rResult.resize( 2, 1 );
                //derivatives of shape functions
                ShapeFunctionsGradientsType shape_functions_gradients =
                    CalculateShapeFunctionsIntegrationPointsLocalGradients( ThisMethod );
                Matrix ShapeFunctionsGradientInIntegrationPoint =
                    shape_functions_gradients( IntegrationPointIndex );
                //values of shape functions in integration points
                vector<double> ShapeFunctionValuesInIntegrationPoint = ZeroVector( 3 );
                ShapeFunctionValuesInIntegrationPoint = row( CalculateShapeFunctionsIntegrationPointsValues( ThisMethod ),
                                                        IntegrationPointIndex );

                //Elements of jacobian matrix (e.g. J(1,1) = dX1/dXi1)
                //loop over all nodes

                for ( unsigned int i = 0; i < this->PointsNumber(); i++ )
                {
                    rResult( 0, 0 ) += ( this->GetPoint( i ).X() ) * ( ShapeFunctionsGradientInIntegrationPoint( i, 0 ) );
                    rResult( 1, 0 ) += ( this->GetPoint( i ).Y() ) * ( ShapeFunctionsGradientInIntegrationPoint( i, 0 ) );
                }

                return rResult;
            }

            /** Jacobian in given point. This method calculate jacobian
             * matrix in given point.
             *
             * @param rPoint point which jacobians has to
             * be calculated in it.
             *
             * @return Matrix of double which is jacobian matrix \f$ J \f$ in given point.
             *
             * @see DeterminantOfJacobian
             * @see InverseOfJacobian
             */
            virtual Matrix& Jacobian( Matrix& rResult, const CoordinatesArrayType& rPoint ) const
            {
                //setting up size of jacobian matrix
                rResult.resize( 2, 1 );
                //derivatives of shape functions
                Matrix shape_functions_gradients;
                shape_functions_gradients = ShapeFunctionsLocalGradients( shape_functions_gradients, rPoint );
                //Elements of jacobian matrix (e.g. J(1,1) = dX1/dXi1)
                //loop over all nodes

                for ( unsigned int i = 0; i < this->PointsNumber(); i++ )
                {
                    rResult( 0, 0 ) += ( this->GetPoint( i ).X() ) * ( shape_functions_gradients( i, 0 ) );
                    rResult( 1, 0 ) += ( this->GetPoint( i ).Y() ) * ( shape_functions_gradients( i, 0 ) );
                }

                return rResult;
            }

            /** Determinant of jacobians for given integration method. This
             * method calculate determinant of jacobian in all
             * integrations points of given integration method.
             *
             * @return Vector of double which is vector of determinants of
             * jacobians \f$ |J|_i \f$ where \f$ i=1,2,...,n \f$ is the
             * integration point index of given integration method.
             *
             * @see Jacobian
             * @see InverseOfJacobian
             */
            virtual Vector& DeterminantOfJacobian( Vector& rResult, IntegrationMethod ThisMethod ) const
            {
                KRATOS_ERROR( std::logic_error, "Jacobian is not square" , "" );
            }

            /** Determinant of jacobian in specific integration point of
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
            virtual double DeterminantOfJacobian( IndexType IntegrationPointIndex, IntegrationMethod ThisMethod ) const
            {
                KRATOS_ERROR( std::logic_error, "Jacobian is not square" , "" );
            }

            /** Determinant of jacobian in given point. This method calculate determinant of jacobian
             * matrix in given point.
             *
             * @param rPoint point which determinant of jacobians has to
             * be calculated in it.
             *
             * @return Determinamt of jacobian matrix \f$ |J| \f$ in given
             * point.
             *
             * @see DeterminantOfJacobian
             * @see InverseOfJacobian
             */
            virtual double DeterminantOfJacobian( const CoordinatesArrayType& rPoint ) const
            {
                KRATOS_ERROR( std::logic_error, "Jacobian is not square" , "" );
            }

            /** Inverse of jacobians for given integration method. This method
             * calculate inverse of jacobians matrices in all integrations points of
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
             */
            virtual JacobiansType& InverseOfJacobian( JacobiansType& rResult, IntegrationMethod ThisMethod ) const
            {
                KRATOS_ERROR( std::logic_error, "Jacobian is not square" , "" );
            }

            /** Inverse of jacobian in specific integration point of given integration
             * method. This method calculate Inverse of jacobian matrix in given
             * integration point of given integration method.
             *
             * @param IntegrationPointIndex index of integration point which inverse of jacobians has to
             * be calculated in it.
             *
             * @param ThisMethod integration method which inverse of jacobians has to
             * be calculated in its integration points.
             *
             * @return Inverse of jacobian matrix \f$ J^{-1}_i \f$ where \f$
             * i \f$ is the given integration point index of given
             * integration method.
             *
             * @see Jacobian
             * @see DeterminantOfJacobian
             */
            virtual Matrix& InverseOfJacobian( Matrix& rResult, IndexType IntegrationPointIndex,
                                               IntegrationMethod ThisMethod ) const
            {
                KRATOS_ERROR( std::logic_error, "Jacobian is not square" , "" );
            }

            /** Inverse of jacobian in given point. This method calculate inverse of jacobian
             * matrix in given point.
             *
             * @param rPoint point which inverse of jacobians has to
             * be calculated in it.
             *
             * @return Inverse of jacobian matrix \f$ J^{-1} \f$ in given point.
             *
             * @see DeterminantOfJacobian
             * @see InverseOfJacobian
             */
            virtual Matrix& InverseOfJacobian( Matrix& rResult, const CoordinatesArrayType& rPoint ) const
            {
                KRATOS_ERROR( std::logic_error, "Jacobian is not square" , "" );
            }

            ///@}
            ///@name Shape Function
            ///@{

            virtual double ShapeFunctionValue( IndexType ShapeFunctionIndex,
                                               const CoordinatesArrayType& rPoint ) const
            {
                switch ( ShapeFunctionIndex )
                {
                    case 0:
                        return( 0.5*( rPoint[0] - 1.0 )*rPoint[0] );
                    case 1:
                        return( 1.0 -rPoint[0]*rPoint[0] );
                    case 2:
                        return( 0.5*( rPoint[0] + 1.0 )*rPoint[0] );
                    default:
                        KRATOS_ERROR( std::logic_error,
                                      "Wrong index of shape function!" ,
                                      *this );
                }

                return 0;
            }



            virtual ShapeFunctionsGradientsType& ShapeFunctionsIntegrationPointsGradients( ShapeFunctionsGradientsType& rResult, IntegrationMethod ThisMethod ) const
            {
                KRATOS_ERROR( std::logic_error, "Jacobian is not square" , "" );
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
                return "1 dimensional line with 3 nodes in 2D space";
            }

            /** Print information about this object.

            @param rOStream Stream to print into it.
            @see PrintData()
            @see Info()
            */
            virtual void PrintInfo( std::ostream& rOStream ) const
            {
                rOStream << "1 dimensional line with 3 nodes in 2D space";
            }

            /** Print geometry's data into given stream. Prints it's points
            by the order they stored in the geometry and then center
            point of geometry.

            @param rOStream Stream to print into it.
            @see PrintInfo()
            @see Info()
            */
            virtual void PrintData( std::ostream& rOStream ) const
            {
                BaseType::PrintData( rOStream );
                std::cout << std::endl;
                Matrix jacobian;
                Jacobian( jacobian, PointType() );
                rOStream << "    Jacobian\t : " << jacobian;
            }

            /**
             * Calculates the local gradients for all integration points for
             * given integration method
             */
            virtual ShapeFunctionsGradientsType ShapeFunctionsLocalGradients(
                IntegrationMethod ThisMethod )
            {
                ShapeFunctionsGradientsType localGradients
                = CalculateShapeFunctionsIntegrationPointsLocalGradients( ThisMethod );
                const int integration_points_number
                = msGeometryData.IntegrationPointsNumber( ThisMethod );
                ShapeFunctionsGradientsType Result( integration_points_number );

                for ( int pnt = 0; pnt < integration_points_number; pnt++ )
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
                = msGeometryData.IntegrationPointsNumber( ThisMethod );
                ShapeFunctionsGradientsType Result( integration_points_number );

                for ( int pnt = 0; pnt < integration_points_number; pnt++ )
                {
                    Result[pnt] = localGradients[pnt];
                }

                return Result;
            }

            /**
             * Calculates the gradients in terms of local coordinates
             * of all shape functions in a given point.
             *
             * @param rPoint the current point at which the gradients are calculated
             * @return the gradients of all shape functions
             * \f$ \frac{\partial N^i}{\partial \xi_j} \f$
             */
            virtual Matrix& ShapeFunctionsLocalGradients( Matrix& rResult,
                    const CoordinatesArrayType& rPoint ) const
            {
                //setting up result matrix
                rResult.resize( 3, 1 );
                noalias( rResult ) = ZeroMatrix( 3, 1 );
                rResult( 0, 0 ) = rPoint[0] - 0.5;
                rResult( 1, 0 ) = -2.0 * rPoint[0];
                rResult( 2, 0 ) = rPoint[0] + 0.5;
                return( rResult );
            }

            /**
             * returns the local coordinates of all nodes of the current geometry
             * @param rResult a Matrix object that will be overwritten by the result
             * @return the local coordinates of all nodes
             */
            virtual Matrix& PointsLocalCoordinates( Matrix& rResult ) const
            {
                rResult.resize( 3, 1 );
                noalias( rResult ) = ZeroMatrix( 3, 1 );
                rResult( 0, 0 ) = -1.0;
                rResult( 1, 0 ) =  0.0;
                rResult( 2, 0 ) =  1.0;
                return rResult;
            }

            /**
             * returns the shape function gradients in an arbitrary point,
             * given in local coordinates
             * @param rResult the matrix of gradients, will be overwritten
             * with the gradients for all
             * shape functions in given point
             * @param rPoint the given point the gradients are calculated in
             */
            virtual Matrix& ShapeFunctionsGradients( Matrix& rResult, CoordinatesArrayType& rPoint )
            {
                rResult.resize( 3, 1 );
                noalias( rResult ) = ZeroMatrix( 3, 1 );

                rResult( 0, 0 ) = rPoint[0] - 0.5;
                rResult( 1, 0 ) = -2.0 * rPoint[0];
                rResult( 2, 0 ) = rPoint[0] + 0.5;
                return rResult;
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


            ///@}

        private:
            ///@name Static Member Variables
            ///@{

            static const GeometryData msGeometryData;

            ///@}
            ///@name Member Variables
            ///@{


            ///@}
            ///@name Serialization
            ///@{

            friend class Serializer;

            virtual void save( Serializer& rSerializer ) const
            {
                rSerializer.save( "Name", "Line2D3" );
                KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, PointsArrayType );
            }

            virtual void load( Serializer& rSerializer )
            {
                KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, PointsArrayType );
            }

            Line2D3(): BaseType( PointsArrayType(), &msGeometryData ) {}

            ///@}
            ///@name Private Operators
            ///@{


            ///@}
            ///@name Private Operations
            ///@{

            static Matrix CalculateShapeFunctionsIntegrationPointsValues( typename BaseType::IntegrationMethod ThisMethod )
            {
                const IntegrationPointsContainerType& all_integration_points = AllIntegrationPoints();
                const IntegrationPointsArrayType& IntegrationPoints = all_integration_points[ThisMethod];
                int integration_points_number = IntegrationPoints.size();
                Matrix N( integration_points_number, 3 );

                for ( int it_gp = 0; it_gp < integration_points_number; it_gp++ )
                {
                    double e = IntegrationPoints[it_gp].X();
                    N( it_gp, 0 ) = 0.5 * ( e - 1 ) * e;
                    N( it_gp, 1 ) = 1.0 - e * e;
                    N( it_gp, 2 ) = 0.5 * ( 1 + e ) * e;
                }

                return N;
            }

            static ShapeFunctionsGradientsType CalculateShapeFunctionsIntegrationPointsLocalGradients(
                typename BaseType::IntegrationMethod ThisMethod )
            {
                const IntegrationPointsContainerType& all_integration_points = AllIntegrationPoints();
                const IntegrationPointsArrayType& IntegrationPoints = all_integration_points[ThisMethod];
                ShapeFunctionsGradientsType DN_De( IntegrationPoints.size() );
                std::fill( DN_De.begin(), DN_De.end(), Matrix( 3, 1 ) );

                for ( unsigned int it_gp = 0; it_gp < IntegrationPoints.size(); it_gp++ )
                {
                    double e = IntegrationPoints[it_gp].X();
                    DN_De[it_gp]( 0, 0 ) = e - 0.5;
                    DN_De[it_gp]( 1, 0 ) = -2.0 * e;
                    DN_De[it_gp]( 2, 0 ) = e + 0.5;
                }

                return DN_De;
            }

            static const IntegrationPointsContainerType AllIntegrationPoints()
            {
                IntegrationPointsContainerType integration_points = {{
                        Quadrature<GaussLegendreIntegrationPoints1, 1, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                        Quadrature<GaussLegendreIntegrationPoints2, 1, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                        Quadrature<GaussLegendreIntegrationPoints3, 1, IntegrationPoint<3> >::GenerateIntegrationPoints()
                    }
                };
                return integration_points;
            }

            static const ShapeFunctionsValuesContainerType AllShapeFunctionsValues()
            {
                ShapeFunctionsValuesContainerType shape_functions_values = {{
                        Line2D3<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::GI_GAUSS_1 ),
                        Line2D3<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::GI_GAUSS_2 ),
                        Line2D3<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::GI_GAUSS_3 )
                    }
                };
                return shape_functions_values;
            }

            static const ShapeFunctionsLocalGradientsContainerType AllShapeFunctionsLocalGradients()
            {
                ShapeFunctionsLocalGradientsContainerType shape_functions_local_gradients = {{
                        Line2D3<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_GAUSS_1 ),
                        Line2D3<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_GAUSS_2 ),
                        Line2D3<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_GAUSS_3 )
                    }
                };
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

            template<class TOtherPointType> friend class Line2D3;

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
    inline std::istream& operator >> ( std::istream& rIStream,
                                       Line2D3<TPointType>& rThis );

    /// output stream function
    template<class TPointType>
    inline std::ostream& operator << ( std::ostream& rOStream,
                                       const Line2D3<TPointType>& rThis )
    {
        rThis.PrintInfo( rOStream );
        rOStream << std::endl;
        rThis.PrintData( rOStream );

        return rOStream;
    }

    ///@}

//   template<class TPointType>
//   const typename Line2D3<TPointType>::IntegrationPointsContainerType Line2D3<TPointType>::msIntegrationPoints = {
//    Quadrature<TriangleGaussianIntegrationPoints<1>, 2, IntegrationPoint<3> >::GenerateIntegrationPoints(),
//    Quadrature<TriangleGaussianIntegrationPoints<2>, 2, IntegrationPoint<3> >::GenerateIntegrationPoints(),
//    Quadrature<TriangleGaussianIntegrationPoints<3>, 2, IntegrationPoint<3> >::GenerateIntegrationPoints()
//   };


//   template<class TPointType>
//   const typename Line2D3<TPointType>::ShapeFunctionsValuesContainerType
//   Line2D3<TPointType>::msShapeFunctionsValues = {
//    Line2D3<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(GeometryData::GI_GAUSS_1),
//    Line2D3<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(GeometryData::GI_GAUSS_2),
//    Line2D3<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(GeometryData::GI_GAUSS_3)
//   };


    //template<class TPointType>
    //const typename GeometryData::ShapeFunctionsLocalGradientsContainerType
    //Line2D3<TPointType>::msShapeFunctionsLocalGradients = {
    // Line2D3<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(GeometryData::GI_GAUSS_1),
    // Line2D3<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(GeometryData::GI_GAUSS_2),
    // Line2D3<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(GeometryData::GI_GAUSS_3)
    //};

    template<class TPointType>
    const GeometryData Line2D3<TPointType>::msGeometryData( 2,
            2,
            1,
            GeometryData::GI_GAUSS_2,
            Line2D3<TPointType>::AllIntegrationPoints(),
            Line2D3<TPointType>::AllShapeFunctionsValues(),
            AllShapeFunctionsLocalGradients() );

}  // namespace Kratos.

#endif // KRATOS_LINE_2D_3_H_INCLUDED  defined 

