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
//   Date:                $Date: 2008-05-06 13:42:56 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_LINE_3D_2_H_INCLUDED )
#define  KRATOS_LINE_3D_2_H_INCLUDED



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

    class Line3D2 : public Geometry<TPointType>
    {

        public:
            ///@}
            ///@name Type Definitions
            ///@{

            /// Geometry as base class.
            typedef Geometry<TPointType> BaseType;

            /// Pointer definition of Line3D2
            KRATOS_CLASS_POINTER_DEFINITION( Line3D2 );

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

            Line3D2( const PointType& FirstPoint, const PointType& SecondPoint )
                    : BaseType( PointsArrayType(), &msGeometryData )
            {
                BaseType::Points().push_back( typename PointType::Pointer( new PointType( FirstPoint ) ) );
                BaseType::Points().push_back( typename PointType::Pointer( new PointType( SecondPoint ) ) );
            }

            Line3D2( typename PointType::Pointer pFirstPoint, typename PointType::Pointer pSecondPoint )
                    : BaseType( PointsArrayType(), &msGeometryData )
            {
                BaseType::Points().push_back( pFirstPoint );
                BaseType::Points().push_back( pSecondPoint );
            }

            Line3D2( const PointsArrayType& ThisPoints )
                    : BaseType( ThisPoints, &msGeometryData )
            {
                if ( BaseType::PointsNumber() != 2 )
                    KRATOS_ERROR( std::invalid_argument,
                                  "Invalid points number. Expected 2, given " , BaseType::PointsNumber() );
            }

            /** Copy constructor.
            Construct this geometry as a copy of given geometry.

            @note This copy constructor don't copy the points and new
            geometry shares points with given source geometry. It's
            obvious that any change to this new geometry's point affect
            source geometry's points too.
            */
            Line3D2( Line3D2 const& rOther )
                    : BaseType( rOther )
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
            template<class TOtherPointType> Line3D2( Line3D2<TOtherPointType> const& rOther )
                    : BaseType( rOther )
            {
            }

            /// Destructor. Do nothing!!!
            virtual ~Line3D2() {}

            GeometryData::KratosGeometryFamily GetGeometryFamily()
            {
                return GeometryData::Kratos_Linear;
            }

            GeometryData::KratosGeometryType GetGeometryType()
            {
                return GeometryData::Kratos_Line3D2;
            }

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
            Line3D2& operator=( const Line3D2& rOther )
            {
                BaseType::operator=( rOther );

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
            Line3D2& operator=( Line3D2<TOtherPointType> const & rOther )
            {
                BaseType::operator=( rOther );

                return *this;
            }

            ///@}
            ///@name Operations
            ///@{

            typename BaseType::Pointer Create( PointsArrayType const& ThisPoints ) const
            {
                return typename BaseType::Pointer( new Line3D2( ThisPoints ) );
            }

            virtual boost::shared_ptr< Geometry< Point<3> > > Clone() const
            {
                Geometry< Point<3> >::PointsArrayType NewPoints;

                //making a copy of the nodes TO POINTS (not Nodes!!!)

                for ( IndexType i = 0 ; i < BaseType::Points().size() ; i++ )
                    NewPoints.push_back( BaseType::Points()[i] );

                //creating a geometry with the new points
                boost::shared_ptr< Geometry< Point<3> > > p_clone( new Line3D2< Point<3> >( NewPoints ) );

                p_clone->ClonePoints();

                return p_clone;
            }

            //lumping factors for the calculation of the lumped mass matrix
            virtual Vector& LumpingFactors( Vector& rResult ) const
            {
                rResult.resize( 2, false );
                rResult[0] = 0.5;
                rResult[1] = 0.5;
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
                double lenght = pow( BaseType::GetPoint( 0 ).X() - BaseType::GetPoint( 1 ).X(), 2 ) + pow( BaseType::GetPoint( 0 ).Y() - BaseType::GetPoint( 1 ).Y(), 2 ) +
                                pow( BaseType::GetPoint( 0 ).Z() - BaseType::GetPoint( 1 ).Z(), 2 );
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
                double lenght = pow( BaseType::GetPoint( 0 ).X() - BaseType::GetPoint( 1 ).X(), 2 ) + pow( BaseType::GetPoint( 0 ).Y() - BaseType::GetPoint( 1 ).Y(), 2 );
                return sqrt( lenght );
            }

//      virtual void Bounding_Box(BoundingBox<TPointType, BaseType>& rResult) const
//              {
//                 //rResult.Geometry() = *(this);
//                 BaseType::Bounding_Box(rResult.LowPoint(), rResult.HighPoint());
//              }

            ///@}
            ///@name Jacobian
            ///@{


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
            virtual JacobiansType& Jacobian( JacobiansType& rResult, IntegrationMethod ThisMethod ) const
            {
                Matrix jacobian( 3, 1 );
                jacobian( 0, 0 ) = ( this->GetPoint( 1 ).X() - this->GetPoint( 0 ).X() ) * 0.5; //on the Gauss points (J is constant at each element)
                jacobian( 1, 0 ) = ( this->GetPoint( 1 ).Y() - this->GetPoint( 0 ).Y() ) * 0.5;
                jacobian( 2, 0 ) = ( this->GetPoint( 1 ).Z() - this->GetPoint( 0 ).Z() ) * 0.5;

                if ( rResult.size() != BaseType::IntegrationPointsNumber( ThisMethod ) )
                {
                    // KLUDGE: While there is a bug in ublas vector resize, I have to put this beside resizing!!
                    JacobiansType temp( BaseType::IntegrationPointsNumber( ThisMethod ) );
                    rResult.swap( temp );
                }

                std::fill( rResult.begin(), rResult.end(), jacobian );

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
            virtual Matrix& Jacobian( Matrix& rResult, IndexType IntegrationPointIndex, IntegrationMethod ThisMethod ) const
            {
                rResult.resize( 3, 1, false );
                //on the Gauss points (J is constant at each element)
                rResult( 0, 0 ) = ( this->GetPoint( 1 ).X() - this->GetPoint( 0 ).X() ) * 0.5;
                rResult( 1, 0 ) = ( this->GetPoint( 1 ).Y() - this->GetPoint( 0 ).Y() ) * 0.5;
                rResult( 2, 0 ) = ( this->GetPoint( 1 ).Z() - this->GetPoint( 0 ).Z() ) * 0.5;
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
            virtual Matrix& Jacobian( Matrix& rResult, const CoordinatesArrayType& rPoint ) const
            {
                rResult.resize( 3, 1, false );
                //on the Gauss points (J is constant at each element)
                rResult( 0, 0 ) = ( this->GetPoint( 1 ).X() - this->GetPoint( 0 ).X() ) * 0.5;
                rResult( 1, 0 ) = ( this->GetPoint( 1 ).Y() - this->GetPoint( 0 ).Y() ) * 0.5;
                rResult( 2, 0 ) = ( this->GetPoint( 1 ).Z() - this->GetPoint( 0 ).Z() ) * 0.5;
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
            virtual Vector& DeterminantOfJacobian( Vector& rResult, IntegrationMethod ThisMethod ) const
            {
                KRATOS_ERROR( std::logic_error, "Jacobian is not square" , "" );
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
            virtual double DeterminantOfJacobian( IndexType IntegrationPointIndex, IntegrationMethod ThisMethod ) const
            {
                KRATOS_ERROR( std::logic_error, "Jacobian is not square" , "" );
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
            virtual double DeterminantOfJacobian( const CoordinatesArrayType& rPoint ) const
            {
                KRATOS_ERROR( std::logic_error, "Jacobian is not square" , "" );
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
            virtual JacobiansType& InverseOfJacobian( JacobiansType& rResult, IntegrationMethod ThisMethod ) const
            {
                KRATOS_ERROR( std::logic_error, "Jacobian is not square" , "" );
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
            virtual Matrix& InverseOfJacobian( Matrix& rResult, IndexType IntegrationPointIndex, IntegrationMethod ThisMethod ) const
            {
                KRATOS_ERROR( std::logic_error, "Jacobian is not square" , "" );
            }

            /** Inverse of jacobian in given point. This method calculate inverse of jacobian
            matrix in given point.

            @param rPoint point which inverse of jacobians has to
            be calculated in it.

            @return Inverse of jacobian matrix \f$ J^{-1} \f$ in given point.

            @see DeterminantOfJacobian
            @see InverseOfJacobian
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
                        return( 0.5*( 1.0 - rPoint[0] ) );

                    case 1:
                        return( 0.5*( 1.0 + rPoint[0] ) );

                    default:
                        KRATOS_ERROR( std::logic_error,
                                      "Wrong index of shape function!" ,
                                      *this );
                }

                return 0;
            }

            virtual Matrix& ShapeFunctionsLocalGradients( Matrix& rResult,
                    const CoordinatesArrayType& rPoint ) const
            {
                //setting up result matrix
                rResult.resize( 2, 1 );
                noalias( rResult ) = ZeroMatrix( 2, 1 );
                rResult( 0, 0 ) = -0.5;
                rResult( 1, 0 ) = 0.5;
                return( rResult );
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
                return "1 dimensional line with 2 nodes in 3D space";
            }

            /** Print information about this object.

            @param rOStream Stream to print into it.
            @see PrintData()
            @see Info()
            */
            virtual void PrintInfo( std::ostream& rOStream ) const
            {
                rOStream << "1 dimensional line with 2 nodes in 3D space";
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
                Matrix N( integration_points_number, 2 );

                for ( int it_gp = 0; it_gp < integration_points_number; it_gp++ )
                {
                    double e = IntegrationPoints[it_gp].X();
                    N( it_gp, 0 ) = 0.5 * ( 1 - e );
                    N( it_gp, 1 ) = 0.5 * ( 1 + e );
                }

                return N;
            }

            static ShapeFunctionsGradientsType CalculateShapeFunctionsIntegrationPointsLocalGradients( typename BaseType::IntegrationMethod ThisMethod )
            {
                const IntegrationPointsContainerType& all_integration_points = AllIntegrationPoints();
                const IntegrationPointsArrayType& IntegrationPoints = all_integration_points[ThisMethod];
                ShapeFunctionsGradientsType DN_De( IntegrationPoints.size() );
                std::fill( DN_De.begin(), DN_De.end(), Matrix( 2, 1 ) );

                for ( unsigned int it_gp = 0; it_gp < IntegrationPoints.size(); it_gp++ )
                {
                    // double e = IntegrationPoints[it_gp].X();
                    DN_De[it_gp]( 0, 0 ) = -0.5;
                    DN_De[it_gp]( 1, 0 ) =  0.5;
                }

                return DN_De;

            }

            static const IntegrationPointsContainerType AllIntegrationPoints()
            {
                IntegrationPointsContainerType integration_points = {{
                        Quadrature<GaussLegendreIntegrationPoints1, 1, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                        Quadrature<GaussLegendreIntegrationPoints2, 1, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                        Quadrature<GaussLegendreIntegrationPoints3, 1, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                        Quadrature<GaussLegendreIntegrationPoints4, 1, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                        Quadrature<GaussLegendreIntegrationPoints5, 1, IntegrationPoint<3> >::GenerateIntegrationPoints()
                    }
                };
                return integration_points;
            }

            static const ShapeFunctionsValuesContainerType AllShapeFunctionsValues()
            {
                ShapeFunctionsValuesContainerType shape_functions_values = {{
                        Line3D2<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::GI_GAUSS_1 ),
                        Line3D2<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::GI_GAUSS_2 ),
                        Line3D2<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::GI_GAUSS_3 ),
                        Line3D2<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::GI_GAUSS_4 ),
                        Line3D2<TPointType>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::GI_GAUSS_5 )
                    }
                };
                return shape_functions_values;
            }

            static const ShapeFunctionsLocalGradientsContainerType AllShapeFunctionsLocalGradients()
            {
                ShapeFunctionsLocalGradientsContainerType shape_functions_local_gradients = {{
                        Line3D2<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_GAUSS_1 ),
                        Line3D2<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_GAUSS_2 ),
                        Line3D2<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_GAUSS_3 ),
                        Line3D2<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_GAUSS_4 ),
                        Line3D2<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_GAUSS_5 ),

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

            template<class TOtherPointType> friend class Line3D2;

            ///@}
            ///@name Un accessible methods
            ///@{

            Line3D2();



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
                                       Line3D2<TPointType>& rThis );

    /// output stream function
    template<class TPointType>
    inline std::ostream& operator << ( std::ostream& rOStream,
                                       const Line3D2<TPointType>& rThis )
    {
        rThis.PrintInfo( rOStream );
        rOStream << std::endl;
        rThis.PrintData( rOStream );

        return rOStream;
    }

    ///@}

//   template<class TPointType>
//   const typename Line3D2<TPointType>::IntegrationPointsContainerType Line3D2<TPointType>::msIntegrationPoints = {
//    Quadrature<TriangleGaussianIntegrationPoints<1>, 2, IntegrationPoint<3> >::GenerateIntegrationPoints(),
//    Quadrature<TriangleGaussianIntegrationPoints<2>, 2, IntegrationPoint<3> >::GenerateIntegrationPoints(),
//    Quadrature<TriangleGaussianIntegrationPoints<3>, 2, IntegrationPoint<3> >::GenerateIntegrationPoints()
//   };


//   template<class TPointType>
//   const typename Line3D2<TPointType>::ShapeFunctionsValuesContainerType
//   Line3D2<TPointType>::msShapeFunctionsValues = {
//    Line3D2<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(GeometryData::GI_GAUSS_1),
//    Line3D2<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(GeometryData::GI_GAUSS_2),
//    Line3D2<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(GeometryData::GI_GAUSS_3)
//   };


    //template<class TPointType>
    //const typename GeometryData::ShapeFunctionsLocalGradientsContainerType
    //Line3D2<TPointType>::msShapeFunctionsLocalGradients = {
    // Line3D2<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(GeometryData::GI_GAUSS_1),
    // Line3D2<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(GeometryData::GI_GAUSS_2),
    // Line3D2<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(GeometryData::GI_GAUSS_3)
    //};

    template<class TPointType>
    const GeometryData Line3D2<TPointType>::msGeometryData( 3,
            3,
            1,
            GeometryData::GI_GAUSS_1,
            Line3D2<TPointType>::AllIntegrationPoints(),
            Line3D2<TPointType>::AllShapeFunctionsValues(),
            AllShapeFunctionsLocalGradients() );

}  // namespace Kratos.

#endif // KRATOS_LINE_3D_2_H_INCLUDED  defined 

