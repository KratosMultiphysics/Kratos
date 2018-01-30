//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Janosch Stascheit
//                   Felix Nagel
//  contributors:    Hoang Giang Bui
//                   Josep Maria Carbonell
//

#if !defined(KRATOS_NURBS_2D_3_H_INCLUDED ) //modified by Matthias
#define  KRATOS_NURBS_2D_3_H_INCLUDED  //modified by Matthias

// System includes

// External includes

// Project includes
#include "includes/node.h"
#include "geometries/nurbs_base_geometry.h"
#include "integration/quadrilateral_gauss_legendre_integration_points.h"

//#include "opennurbs.h"

//Constants needed by OpenNURBS-Library
#define ON_SQRT_EPSILON 1.490116119385000000e-8

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

template<class TPointType> class NurbsPatchGeometry2D  //modified by Matthias
    : public NurbsPatchGeometry<TPointType>
{
public:
    ///@}
    ///@name Type Definitions
    ///@{

    /**
     * Geometry as base class.
     */
    typedef NurbsPatchGeometry<TPointType> BaseType;


    /**
     * Pointer definition of NurbsSurfaceGeometry3D
     */
    KRATOS_CLASS_POINTER_DEFINITION( NurbsPatchGeometry2D ); //modified by Matthias

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
    * A third order tensor to hold shape functions' local third derivatives.
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

    /**
     * A standard constructor which will create an empty NURBS-surface
     */

    NurbsPatchGeometry2D():
        BaseType( PointsArrayType(), &msGeometryData),
        mKnotsXi(),
        mKnotsEta(),
        mWeights(),
        mPolynomialDegreeP(),
        mPolynomialDegreeQ(),
        mNumberOfCPsU(),
        mNumberOfCPsV()
    {}

    /**
     * A constructor which reads in a NURBS-Element.
     * The needed information to identify a NURBS element are:
     *
     * @param ControlPoints: A PointerVector of pointers to node<3>, which contain x-y-z- coordinate of each Control Point
     * @param Weights: A Vector which stores each weight of the Control Point (same order as in PointerVector ControlPoints)
     * @param KnotsXi: A Vector of doubles which contains the Knot Vector in Xi- direction
     * @param KnotsEta: A Vector of doubles which contains the Knot Vector in Eta- direction
     * @param p: A Double which is the polynomial degree in Xi- direction
     * @param q: A Double which is the polynomial degree in Eta- direction
     * @param NumberOfCpXi: Number of Control Points in Xi- direction
     * @param NumberOfCpEta: Number of Control Points in Eta- direction
     * @param LowerElementBoundaryXi: Lowest possible value of Xi for the element
     * @param HigherElementBoundaryXi: Highest possible value of Xi for the element
     * @param LowerElementBoundaryEta: Lowest possible value of Eta for the element
     * @param HigherElementBoundarEta: Highest possible value of Eta for the element
     *
     * @note: Elements in the frame of IGA are defined on knot spans (meaning to consecutive knots which have different values)
     *
     */


    NurbsPatchGeometry2D( const PointerVector< Node<3> > ControlPoints,
                            const Vector Weights,
                            const Vector KnotsXi,
                            const Vector KnotsEta,
                            int p,
                            int q,
                            int NumberOfCpXi,
                            int NumberOfCpEta,
                            double LowerElementBoundaryXi,
                            double HigherElementBoundaryXi,
                            double LowerElementBoundaryEta,
                            double HigherElementBoundarEta)
           : BaseType( PointsArrayType(ControlPoints), &msGeometryData ),  mKnotsXi(KnotsXi), mKnotsEta(KnotsEta) ,mWeights(Weights), mPolynomialDegreeP(p), mPolynomialDegreeQ(q),
             mNumberOfCPsU(NumberOfCpXi), mNumberOfCPsV(NumberOfCpEta)
    {
        mLowerXi = LowerElementBoundaryXi;
        mUpperXi = HigherElementBoundaryXi;
        mLowerEta = LowerElementBoundaryEta;
        mUpperEta = HigherElementBoundarEta;
        mNumberOfGeometries = 1;
    }


    /**
     * A constructor which reads in a NURBS-patch.
     * The needed information to identify a NURBS patch are:
     *
     * @param ControlPoints: A PointerVector of pointers to node<3>, which contain x-y-z- coordinate of each Control Point
     * @param Weights: A Vector which stores each weight of the Control Point (same order as in PointerVector ControlPoints)
     * @param KnotsXi: A Vector of doubles which contains the Knot Vector in Xi- direction
     * @param KnotsEta: A Vector of doubles which contains the Knot Vector in Eta- direction
     * @param p: A Double which is the polynomial degree in Xi- direction
     * @param q: A Double which is the polynomial degree in Eta- direction
     * @param NumberOfCpXi: Number of Control Points in Xi- direction
     * @param NumberOfCpEta: Number of Control Points in Eta- direction
     *
     * @note:   As the same class for elements and patches is used, the number of elements will be initialized to 1
     *          but once calling NurbsPatchGeometry2D::DefineGeometries() the number of elements will be set to the
     *          correct value.
     *
     *          The values for Upper/Lower Xi/Eta are initialized to first and last entry of the corresponding Knot Vector
     */

    NurbsPatchGeometry2D( const PointerVector< Node<3> > ControlPoints,
                            const Vector Weights,
                            const Vector KnotsXi,
                            const Vector KnotsEta,
                            int p,
                            int q,
                            int NumberOfCpXi,
                            int NumberOfCpEta)
           : BaseType( PointsArrayType(ControlPoints), &msGeometryData ),  mKnotsXi(KnotsXi), mKnotsEta(KnotsEta) ,mWeights(Weights), mPolynomialDegreeP(p), mPolynomialDegreeQ(q),
             mNumberOfCPsU(NumberOfCpXi), mNumberOfCPsV(NumberOfCpEta)
    {
        if (mKnotsXi.size() != 0)
        {
        mLowerXi = mKnotsXi[0];
        mUpperXi = mKnotsXi[mKnotsXi.size()-1];
        mLowerEta = mKnotsEta[0];
        mUpperEta = mKnotsEta[mKnotsEta.size()-1];
        mNumberOfGeometries = 1;
        }
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
    NurbsPatchGeometry2D( NurbsPatchGeometry2D const& rOther )
        : BaseType( rOther ),
          mKnotsXi(rOther.mKnotsXi),
          mKnotsEta(rOther.mKnotsEta),
          mWeights(rOther.mWeights),
          mPolynomialDegreeP(rOther.mPolynomialDegreeP),
          mPolynomialDegreeQ(rOther.mPolynomialDegreeQ),
          mNumberOfCPsU(rOther.mNumberOfCPsU),
          mNumberOfCPsV(rOther.mNumberOfCPsV)
 {
        std::cout << "Copy" << std::endl;
     if (mKnotsXi.size() != 0)
     {
     mLowerXi = mKnotsXi[0];
     mUpperXi = mKnotsXi[mKnotsXi.size()-1];
     mLowerEta = mKnotsEta[0];
     mUpperEta = mKnotsEta[mKnotsEta.size()-1];
     }
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
    template<class TOtherPointType> NurbsPatchGeometry2D( NurbsPatchGeometry2D<TOtherPointType> const& rOther )
        : BaseType( rOther )
    {
    }

    /**
     * Destructor. Does nothing!!!
     */
    virtual ~NurbsPatchGeometry2D() {}

    GeometryData::KratosGeometryFamily GetGeometryFamily() const override
    {
        return GeometryData::Kratos_Quadrilateral;
    }

    GeometryData::KratosGeometryType GetGeometryType() const override
    {
        return GeometryData::Kratos_Triangle3D3;
    }

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



    NurbsPatchGeometry2D& operator=( const NurbsPatchGeometry2D& rOther )

    {
        //std::cout << "Equal" << std::endl;
        BaseType::operator=( rOther );
        this->mKnotsXi = (rOther.mKnotsXi);
        this->mKnotsEta = (rOther.mKnotsEta);
        this->mWeights = (rOther.mWeights);
        this->mPolynomialDegreeP = (rOther.mPolynomialDegreeP);
        this->mPolynomialDegreeQ = (rOther.mPolynomialDegreeQ);
        this->mNumberOfCPsU = (rOther.mNumberOfCPsU);
        this->mNumberOfCPsV = (rOther.mNumberOfCPsV);

                if (this->mKnotsXi.size() != 0)
                {
                this->mLowerXi = this->mKnotsXi[0];
                this->mUpperXi = this->mKnotsXi[mKnotsXi.size()-1];
                this->mLowerEta = this->mKnotsEta[0];
                this->mUpperEta = this->mKnotsEta[mKnotsEta.size()-1];
                }
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
    NurbsPatchGeometry2D& operator=( NurbsPatchGeometry2D<TOtherPointType> const & rOther )
    {
        BaseType::operator=( rOther );
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{






    typename BaseType::Pointer Create( PointerVector< Point > ,Vector const &weights, Vector const& KnotsXi, Vector const& KnotsEta ) const 
    {
        return typename BaseType::Pointer( new NurbsPatchGeometry2D(PointerVector< Point >(),KnotsXi, KnotsEta) );
    }


    //lumping factors for the calculation of the lumped mass matrix
    virtual Vector& LumpingFactors( Vector& rResult ) const override
    {
        rResult.resize( 3, false );
        std::fill( rResult.begin(), rResult.end(), 1.00 / 3.00 );
        return rResult;
    }

    ///@}
    ///@name Information
    ///@{





    ///@}
    ///@name Jacobian
    ///@{
    /**
     * TODO: implemented but not yet tested
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
                                     IntegrationMethod ThisMethod ) const override
    {
        const unsigned int integration_points_number =
            msGeometryData.IntegrationPointsNumber( ThisMethod );

        if ( integration_points_number == 0 )
            KRATOS_ERROR << "This integration method is not supported" << *this << std::endl;


        //calculating the local gradients and the shape functions derivatives
        ShapeFunctionsGradientsType d_Xi_shape_f_values( integration_points_number );
        ShapeFunctionsGradientsType d_Eta_shape_f_values( integration_points_number );

        //Now also gives Shape Function Values
        rResult = CalculateShapeFunctionsIntegrationPointsLocalGradients( ThisMethod,
                                                                          d_Xi_shape_f_values,
                                                                          d_Eta_shape_f_values);
        return rResult;
    }





    /**
     * TODO: implemented but not yet tested
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
                              IntegrationMethod ThisMethod ) const override
    {

        IntegrationPointsContainerType all_integration_points = AllIntegrationPoints();
        IntegrationPointsArrayType integration_points = all_integration_points[ThisMethod];
        double XiLocalCoordinates,EtaLocalCoordinates;
        rResult.resize(2,2);
        rResult.resize(2,2, false);
        XiLocalCoordinates = (mUpperXi - mLowerXi) / 2 * (1+integration_points[IntegrationPointIndex].X()) + mLowerXi;
        EtaLocalCoordinates = (mUpperEta - mLowerEta) / 2 * (1+integration_points[IntegrationPointIndex].Y())+mLowerEta;
        rResult = ElementGeometryDerivatives(XiLocalCoordinates,EtaLocalCoordinates);

        return rResult;
    }

    /**
     * TODO: implemented but not yet tested
     */
    /**
       * Jacobian in given point. This method calculate jacobian
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
    virtual Matrix& Jacobian( Matrix& rResult,
                              const CoordinatesArrayType& rPoint ) const
    {

        rResult.resize(2,2);
        rResult.resize(2,2, false);
        rResult = ElementGeometryDerivatives(rPoint[0],rPoint[1]);
        return rResult;
    }

    /**
     * TODO: implemented but not yet tested
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
                                           IntegrationMethod ThisMethod ) const override
    {
        JacobiansType Jacobian2D;
        Jacobian(Jacobian2D, ThisMethod );
        const unsigned int integration_points_number =  msGeometryData.IntegrationPointsNumber( ThisMethod );
        rResult.resize(integration_points_number);

        for (unsigned int i=0; i< integration_points_number; i++)
        {
            rResult[i] = Jacobian2D[i](0,0)*Jacobian2D[i](1,1)-Jacobian2D[i](1,0) * Jacobian2D[i](0,1);
        }


        return rResult;
    }

    /**
     * TODO: implemented but not yet tested
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
                                          IntegrationMethod ThisMethod ) const override
    {
        Matrix Jacobian2D(2,2);
        Jacobian2D = ZeroMatrix(2,2);
        Jacobian( Jacobian2D,IntegrationPointIndex,ThisMethod);
        double DetJacobian(0);
        DetJacobian = Jacobian2D(0,0)*Jacobian2D(1,1)*Jacobian2D(1,0)*Jacobian2D(0,1);
        return DetJacobian;
    }


    /**
     * TODO: implemented but not yet tested
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
    virtual double DeterminantOfJacobian( const CoordinatesArrayType& rPoint ) const override
    {
        Matrix Jacobian2D(2,2);
        Jacobian2D = ZeroMatrix(2,2);
        Jacobian(Jacobian2D, rPoint );
        double DetJacobian(0);
        DetJacobian = Jacobian2D(0,0)*Jacobian2D(1,1)*Jacobian2D(1,0)*Jacobian2D(0,1);
        return DetJacobian;
    }

    /**
     * TODO: implemented but not yet tested
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
            IntegrationMethod ThisMethod ) const override
    {
        JacobiansType Jacobian2D;
         Jacobian( Jacobian2D,ThisMethod );
         const unsigned int integration_points_number = msGeometryData.IntegrationPointsNumber( ThisMethod );
         double det(0);

         for (unsigned int i=0; i<integration_points_number ; i++)
         {
             MathUtils<double>::InvertMatrix2(Jacobian2D(i),rResult(i),det);

         }


         return rResult;
    }

    /**
     * TODO: implemented but not yet tested
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
                                       IntegrationMethod ThisMethod ) const override
    {
        Matrix Jacobian2D(2,2);
        Jacobian2D = ZeroMatrix(2,2);
        Jacobian(Jacobian2D,IntegrationPointIndex,ThisMethod );
        double det(0);
        MathUtils<double>::InvertMatrix2(Jacobian2D,rResult,det);

        return rResult;
    }

    /**
     * TODO: implemented but not yet tested
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
                                       const CoordinatesArrayType& rPoint ) const override
    {
        Matrix Jacobian2D(2,2);
        Jacobian2D = ZeroMatrix(2,2);
        Jacobian(Jacobian2D,rPoint );
        double det(0);
        MathUtils<double>::InvertMatrix2(Jacobian2D,rResult,det);
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
    virtual SizeType EdgesNumber() const override
    {
        KRATOS_ERROR << "Nurbs_2d::EdgesNumber" << "No Edges defined for NURBS-surfaces" << std::endl;
        return 0;
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
    virtual GeometriesArrayType Edges( void ) override
    {
        GeometriesArrayType edges = GeometriesArrayType();
        KRATOS_ERROR << "Nurbs_2d::Edges" << "No Edges defined for NURBS-surfaces"  << std::endl;
        return edges;

    }

    virtual SizeType FacesNumber() const
    {
        return 0;
    }

    /**
     * Returns all faces of the current geometry.
     * This is only implemented for 3D geometries, since 2D geometries
     * only have edges but no faces
     * @see EdgesNumber
     * @see Edges
     * @see FacesNumber
    */
    virtual GeometriesArrayType Faces( void ) override
    {
        return GeometriesArrayType();
    }


    /**
     * ConvertArrayToVector
     *
     * This is a function which will be used to transform a c-array to a Kratos::Vector.
     * Will be needed to use the OpenNURBS library (wrapper)
     *
     * @param array: specifies which array will be transformed
     * @param size: size of the array
     * @return: Vector which contains the same information as the array passed to this function
     */

    Vector ConvertArrayToVector(double *array, int size)
    {
        Vector vector(size);
        for (int i=0; i<size; i++)
        {
            vector(i) = array[i];
        }

        return vector;
    }




    /**
     * WriteKnots: Will print all Knots in Xi as well as in Eta direction on the console
     *
     */
    void WriteKnots()
    {
        std::cout<< "Knot Vector u-direction: [ ";
        for (int i=0 ; i<mKnotsXi.size();i++)
        {
          std::cout<< mKnotsXi[i]<<" ";
        }
          std::cout<< "]"<<std::endl;

          std::cout<< "Knot Vector v-direction: [ ";
          for (int i=0 ; i<mKnotsEta.size();i++)
          {
            std::cout<< mKnotsEta[i]<<" ";
          }
            std::cout<< "]"<<std::endl;
    }

    /**
     * WriteCPs: Will write all Control Points on the console
     */
    void WriteCPs()
    {
        double x(0),y(0),z(0);
        for (int i=0 ; i<this->Points().size();i++)
        {
          std::cout<< "CPs-ID: "<<i+1<<std::endl;

            x= this->Points()[i].X();
            y= this->Points()[i].Y();

          std::cout<< "Coordinates: x = "<< x <<" y = "<< y<< " z = "<< z <<std::endl;
        }

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
        return "2 dimensional NURBS-surface";
    }

    /**
     * Print information about this object.
     * @param rOStream Stream to print into it.
     * @see PrintData()
     * @see Info()
     */
    virtual void PrintInfo( std::ostream& rOStream ) const
    {
        rOStream << "2 dimensional NURBS-surface";
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
    virtual void PrintData( std::ostream& rOStream ) const
    {
        BaseType::PrintData( rOStream );
        std::cout << std::endl;
        Matrix jacobian;
        Jacobian( jacobian, PointType() );
        rOStream << "    Jacobian in the origin\t : " << jacobian;
    }



    /**
     * ShapeFunctionsValues: Will compute all shape functions values at all Gauss Points
     * @param ThisMethod: Specifies which Gauss Scheme is used
     * @return  The Shape Functions Values will be stored in a Matrix \f$ F_i_j \f$ where
     *          i is the given integration point index of the integration method and j defines
     *          the shape function index.
     */
    const Matrix& ShapeFunctionsValues(IntegrationMethod ThisMethod) const
    {
        IntegrationPointsContainerType all_integration_points =
            AllIntegrationPoints();
        IntegrationPointsArrayType integration_points =
            all_integration_points[ThisMethod];
        //number of integration points
        const unsigned int integration_points_number = integration_points.size();
        double XiLocalCoordinates,EtaLocalCoordinates;
        Matrix rResult[integration_points_number][(mPolynomialDegreeP+1)*(mPolynomialDegreeQ+1)];
        Matrix LocalShapeFunctionsValues[mPolynomialDegreeP+1][mPolynomialDegreeQ+1];
        for (unsigned int pnt = 0; pnt < integration_points_number; pnt++ )
        {
            //Map from Gaußdomain to Parameterdomain
            XiLocalCoordinates = (mUpperXi - mLowerXi) / 2 * (1+integration_points[pnt].X())+mLowerXi;
            EtaLocalCoordinates = (mUpperEta - mLowerEta) / 2 * (1+integration_points[pnt].Y())+mLowerEta;
            //Initialize Zero Matrix
            LocalShapeFunctionsValues = ZeroMatrix(mPolynomialDegreeP+1,mPolynomialDegreeQ+1);
            //Calculate Shapefunction values at Gaußpoint
            LocalShapeFunctionsValues = ElementGeometryNurbsFunctionsValues(XiLocalCoordinates,EtaLocalCoordinates);

            //Assemble shapefunctionvalues of Gaußpoints to ONE Matrix -> rResult
            for (unsigned int i=0; i<mPolynomialDegreeQ+1; i++)
            {
                for (unsigned int j=0; j<mPolynomialDegreeP+1; j++)
                {
                    rResult[pnt][i+j*(mPolynomialDegreeP+1)] = LocalShapeFunctionsValues[i][j];
                }
            }

        }

        return rResult;
    }


    /**
     * @brief ShapeFunctionsValues: Computes the potentially non-zero shape functions values
     *          at a given point in the local parameter space
     * @param rResult: Vector which will be filled with the shape functions values
     * @param Xi: Xi-Coordinate in the local parameter space where the shape functions will be evaluated
     * @param Eta: Eta-Coordinate in the local parameter space where the shape functions will be evaluated
     * @return Vector of doubles where F_i provides the value of the i-th shape function
     *
     * @note: Potentially at a B-Spline/NURBS Basis Function (p+1) shape functions are unequal to 0. As
     *          We are dealing with Surfaces (Vector Product of two NURBS-Basis Functions) (p+1) * (q+1)
     *          shape functions are potentially unequal to 0.
     */

    Vector& ShapeFunctionsValues(Vector &rResult,
                                const CoordinatesArrayType& rCoordinates) const
    {
        const double& Xi = rCoordinates[0];
        const double& Eta = rCoordinates[1];
        
        rResult.resize((mPolynomialDegreeP+1)*(mPolynomialDegreeQ+1));
        rResult.resize((mPolynomialDegreeP+1)*(mPolynomialDegreeQ+1), false);
        rResult = ZeroVector((mPolynomialDegreeP+1)*(mPolynomialDegreeQ+1));
        Matrix ShapeFunctionsValue = ElementGeometryNurbsFunctionsValues(Xi,Eta);
        for (unsigned int i=0;i<mPolynomialDegreeP+1;i++)
        {
            for (int j=0;j<mPolynomialDegreeQ+1;j++)
            {
                rResult[i+(mPolynomialDegreeP+1)*j] = ShapeFunctionsValue(i,j);
            }
        }
        return rResult;
    }




    /**
     * TODO: implemented but not yet tested
     */
    /**
     * Calculates the Gradients of the shape functions.
     * Calculates the gradients of the shape functions with regard to
     * the global coordinates in all
     * integration points (\f$ \frac{\partial N^i}{\partial X_j} \f$)
     *
     * @param rResult a container which takes the calculated gradients
     * @param ThisMethod the given IntegrationMethod
     * @param ShapeFunctionsValues a Matrix which will be filled with all shape functions values
     *        at all integration points
     * @param determinants_of_jacobian a Vector which will hold all determinants at the integration points
     * @return the gradients of all shape functions with regard to the global coordinates
     * KLUDGE: method call only works with explicit JacobiansType rather than creating
     * JacobiansType within argument list
     *
     * @note: The shape functions values will be calculated here as well, because the superclass
     *        geometry.h doesn't provide any virtual function to calculate all shape functions values
     *        at all integration points.
    */

    virtual ShapeFunctionsGradientsType& ShapeFunctionsIntegrationPointsGradients( ShapeFunctionsGradientsType& rResult,
                                                                                   Vector& determinants_of_jacobian,
                                                                                   IntegrationMethod ThisMethod,
                                                                                   Matrix& ShapeFunctionsValues) const
    {
        const unsigned int integration_points_number = msGeometryData.IntegrationPointsNumber( ThisMethod );

        if ( integration_points_number == 0 )
            KRATOS_ERROR << "This integration method is not supported" << *this << std::endl;

        //workaround by riccardo
        if ( rResult.size() != integration_points_number )
        {
            // KLUDGE: While there is a bug in ublas
            // vector resize, I have to put this beside resizing!!
            ShapeFunctionsGradientsType temp( integration_points_number );
            rResult.swap( temp );
        }

        if (determinants_of_jacobian.size() != integration_points_number )
        {
            determinants_of_jacobian.resize(integration_points_number);
            determinants_of_jacobian.resize(integration_points_number, false);
        }

        //Allocating the memory for the Values of the shape functions
        ShapeFunctionsValues.resize( integration_points_number,(mPolynomialDegreeP+1)*(mPolynomialDegreeQ+1));
        ShapeFunctionsValues.resize( integration_points_number,(mPolynomialDegreeP+1)*(mPolynomialDegreeQ+1), false);
        ShapeFunctionsValues = ZeroMatrix(integration_points_number,(mPolynomialDegreeP+1)*(mPolynomialDegreeQ+1));

        //calculating the local gradients and the shape functions derivatives
        ShapeFunctionsGradientsType d_Xi_shape_f_values( integration_points_number );
        ShapeFunctionsGradientsType d_Eta_shape_f_values( integration_points_number );

        //Now also gives Shape Function Values
        ShapeFunctionsGradientsType locG = CalculateShapeFunctionsIntegrationPointsLocalGradients( ThisMethod,
                                                                                                   d_Xi_shape_f_values,
                                                                                                   d_Eta_shape_f_values,
                                                                                                   ShapeFunctionsValues);
        JacobiansType invJ (integration_points_number);

       // JacobiansType as well as ShapeFunctionsGradientsType contain a VECTOR of Matrices, therefore just one element
       // of the element can be passed to the function InvertMatrix2

        //loop over all integration points
        for ( unsigned int pnt = 0; pnt < integration_points_number; pnt++ )
        {


            //getting the inverse jacobian matrices
            MathUtils<double>::InvertMatrix2(locG(pnt),invJ(pnt),determinants_of_jacobian[pnt]);
            //Applying the determinant of the coordinate transformation from gauß to parameter domain
            determinants_of_jacobian[pnt] = determinants_of_jacobian[pnt] * (mUpperXi-mLowerXi)*0.25*(mUpperEta-mLowerEta);


            //Allocating the memory for the Derivatives of the shape functions
            rResult[pnt].resize( this->Points().size(),2 );
            rResult[pnt].resize( this->Points().size(),2 , false);
            rResult[pnt] = ZeroMatrix(this->Points().size(),2);


            for ( int i = 0; i < mPolynomialDegreeQ+1; i++ )
            {
                for ( int j = 0; j < mPolynomialDegreeP+1; j++ )
                {
                   for (int k = 0; k<2; k++)
                   {
                        rResult[pnt]( i*(mPolynomialDegreeP+1)+j,k ) =
                            ( d_Xi_shape_f_values[pnt]( j,i ) * invJ[pnt](k,0))//( 0, k ) )
                            + ( d_Eta_shape_f_values[pnt]( j,i ) * invJ[pnt](k,1));//( 1, k ) );
                    }
                }
            }

        }//end of loop over integration points

        return rResult;
    }

    virtual ShapeFunctionsGradientsType& ShapeFunctionsIntegrationPointsGradients(
            ShapeFunctionsGradientsType& rResult,
            Vector& determinants_of_jacobian,
            IntegrationMethod ThisMethod ) const
    {
        const unsigned int integration_points_number = msGeometryData.IntegrationPointsNumber( ThisMethod );
        Matrix ShapeFunctionsValues;
        ShapeFunctionsValues.resize( integration_points_number,(mPolynomialDegreeP+1)*(mPolynomialDegreeQ+1));
        ShapeFunctionsValues.resize( integration_points_number,(mPolynomialDegreeP+1)*(mPolynomialDegreeQ+1),false);
        ShapeFunctionsValues = ZeroMatrix(integration_points_number,(mPolynomialDegreeP+1)*(mPolynomialDegreeQ+1));
        rResult = ShapeFunctionsIntegrationPointsGradients(rResult,
                                                           determinants_of_jacobian,
                                                           ThisMethod,
                                                           ShapeFunctionsValues);
                return rResult;
    }

    virtual ShapeFunctionsGradientsType& ShapeFunctionsIntegrationPointsGradients(
        ShapeFunctionsGradientsType& rResult,
        IntegrationMethod ThisMethod ) const
    {
        Vector DeterminantsOfJacobian;
        this->ShapeFunctionsIntegrationPointsGradients(rResult,DeterminantsOfJacobian,ThisMethod);
        return rResult;
    }


    /**
     * Calculates the gradients in terms of local coordinates
     * of all shape functions at given Xi and Eta.
     *
     * @param Xi defines the Xi-value where the shape functions derivatives are evaluated
     * @param Eta defines the Eta-value where the shape functions derivatives are evaluated
     * @return the gradients of all shape functions
     *
     * @note: For exact information what is returned, refer to ElementGeometryDerivatives()
     */

        Matrix& ShapeFunctionsLocalGradients(   double Xi,
                                                double Eta,
                                                Matrix &rResult) const
    {
        rResult.resize(mPolynomialDegreeP+1,mPolynomialDegreeQ+1);
        rResult.resize(mPolynomialDegreeP+1,mPolynomialDegreeQ+1,false);
        rResult = ElementGeometryDerivatives(Xi,Eta);
        return rResult;

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
     * returns the shape function gradients in an arbitrary point,
     * given in local coordinates
     *
     * @param rResult the matrix of gradients,
     * will be overwritten with the gradients for all
     * shape functions in given point
     * @param rPoint the given point the gradients are calculated in
     */
    virtual Matrix& ShapeFunctionsGradients( Matrix& rResult,
                                             PointType& rPoint )
    {
        rResult.resize(mPolynomialDegreeP+1,mPolynomialDegreeQ+1);
        rResult.resize(mPolynomialDegreeP+1,mPolynomialDegreeQ+1,false);
        rResult = ElementGeometryDerivatives(rPoint[0],rPoint[1]);
        return rResult;
    }

	/**
	* returns the shape function gradients in an arbitrary point,
	* given in local coordinates
	*
	* @param rResult the matrix of gradients,
	* will be overwritten with the gradients for all
	* shape functions in given point
	* @param rCoordinates the given point the gradients are calculated in
	*/
	virtual Matrix& ShapeFunctionsGradients(Matrix& rResult,
		const CoordinatesArrayType& rCoordinates)
	{
		rResult.resize(mPolynomialDegreeP + 1, mPolynomialDegreeQ + 1);
		rResult = ElementGeometryDerivatives(rCoordinates[0], rCoordinates[1]);
		return rResult;
	}

	/**
	* @brief ShapeFunctionsDerivativesValues: 
	* @param rResult: Vector which will be filled with the shape functions values
	* @param rCoordinates: the position of the evaluation. rCoordinates[0] is the Xi / rCoordinates[1] is the Eta
	* @return Vector of doubles where F_i provides the value of the i-th shape function
	*
	* @note: Potentially at a B-Spline/NURBS Basis Function (p+1) shape functions are unequal to 0. As
	*          We are dealing with Surfaces (Vector Product of two NURBS-Basis Functions) (p+1) * (q+1)
	*          shape functions are potentially unequal to 0.
	*/

	Matrix& ShapeFunctionsDerivativesValues(Matrix &rResult,
		const CoordinatesArrayType& rCoordinates) const
	{
		//std::cout<<"Vector& ShapeFunctionsValues (nurbs_2d.h)"<<std::endl;
		const double& Xi = rCoordinates[0];
		const double& Eta = rCoordinates[1];

		Matrix NurbsBasisFunctionDerivativesXi = ZeroMatrix((mPolynomialDegreeP + 1), (mPolynomialDegreeQ + 1));
		Matrix NurbsBasisFunctionDerivativesEta = ZeroMatrix((mPolynomialDegreeP + 1), (mPolynomialDegreeQ + 1));

		rResult.resize((mPolynomialDegreeP + 1)*(mPolynomialDegreeQ + 1),2);
		rResult = ZeroMatrix((mPolynomialDegreeP + 1)*(mPolynomialDegreeQ + 1),2);
		
		Matrix Jacobian2D = ElementGeometryDerivatives(Xi, Eta, NurbsBasisFunctionDerivativesXi, NurbsBasisFunctionDerivativesEta);

		for (unsigned int i = 0; i<mPolynomialDegreeP + 1; i++)
		{
			for (int j = 0; j<mPolynomialDegreeQ + 1; j++)
			{
				rResult(i + (mPolynomialDegreeP + 1)*j,0) = NurbsBasisFunctionDerivativesXi(i, j);
				rResult(i + (mPolynomialDegreeP + 1)*j,1) = NurbsBasisFunctionDerivativesEta(i, j);
			}
		}
		return rResult;
	}


	Matrix ShapeFunctionsSecondDerivativesValues(Matrix &rResult,
		const CoordinatesArrayType& rCoordinates)
	{
		//std::cout<<"Vector& ShapeFunctionsValues (nurbs_2d.h)"<<std::endl;
		const double& Xi = rCoordinates[0];
		const double& Eta = rCoordinates[1];
		Vector NurbsBasisFunction = ZeroVector((mPolynomialDegreeP + 1)*(mPolynomialDegreeQ + 1));

		Matrix NurbsBasisFunctionDerivativesXi = ZeroMatrix((mPolynomialDegreeP + 1), (mPolynomialDegreeQ + 1));
		Matrix NurbsBasisFunctionDerivativesEta = ZeroMatrix((mPolynomialDegreeP + 1), (mPolynomialDegreeQ + 1));

		Matrix NurbsBasisFunctionSecondDerivativesXi = ZeroMatrix((mPolynomialDegreeP + 1), (mPolynomialDegreeQ + 1));
		Matrix NurbsBasisFunctionSecondDerivativesEta = ZeroMatrix((mPolynomialDegreeP + 1), (mPolynomialDegreeQ + 1));
		Matrix NurbsBasisFunctionSecondDerivativesXiEta = ZeroMatrix((mPolynomialDegreeP + 1), (mPolynomialDegreeQ + 1));


		rResult.resize((mPolynomialDegreeP + 1)*(mPolynomialDegreeQ + 1), 5);
		rResult = ZeroMatrix((mPolynomialDegreeP + 1)*(mPolynomialDegreeQ + 1), 5);

		ElementGeometrySecondDerivatives(Xi, Eta,
			NurbsBasisFunctionSecondDerivativesXi, 
			NurbsBasisFunctionSecondDerivativesEta,
			NurbsBasisFunctionSecondDerivativesXiEta, 
			NurbsBasisFunctionDerivativesXi, 
			NurbsBasisFunctionDerivativesEta, 
			NurbsBasisFunction);


		for (unsigned int i = 0; i<mPolynomialDegreeP + 1; i++)
		{
			for (int j = 0; j<mPolynomialDegreeQ + 1; j++)
			{
				rResult(i + (mPolynomialDegreeP + 1)*j, 0) = NurbsBasisFunctionDerivativesXi(i, j);
				rResult(i + (mPolynomialDegreeP + 1)*j, 1) = NurbsBasisFunctionDerivativesEta(i, j);
				rResult(i + (mPolynomialDegreeP + 1)*j, 2) = NurbsBasisFunctionSecondDerivativesXi(i, j);
				rResult(i + (mPolynomialDegreeP + 1)*j, 3) = NurbsBasisFunctionSecondDerivativesEta(i, j);
				rResult(i + (mPolynomialDegreeP + 1)*j, 4) = NurbsBasisFunctionSecondDerivativesXiEta(i, j);
			}
		}
		//KRATOS_WATCH(rResult)

		return rResult;
	}

    /**
     * returns the second order derivatives of all shape functions
     * in given arbitrary pointers
     * @param rResult a third order tensor which contains the second derivatives
     * @param rPoint the given point the second order derivatives are calculated in
     */
    //virtual ShapeFunctionsSecondDerivativesType& ShapeFunctionsSecondDerivatives( int order,  Vector N, Matrix *rResult, double t ) const
    virtual Matrix& ShapeFunctionsSecondDerivatives( int order,
                                                     Vector N,
                                                     Matrix *rResult,
                                                     double t ) const 
    {
        KRATOS_ERROR << "Nurbs_2d::ShapeFunctionsSecondDerivatives" << "Second order derivatives not yet implemented" << std::endl;
        return *rResult;
    }

    /**
     * returns the third order derivatives of all shape functions
     * in given arbitrary pointers
     * @param rResult a fourth order tensor which contains the third derivatives
     * @param rPoint the given point the third order derivatives are calculated in
     */
    //virtual ShapeFunctionsThirdDerivativesType& ShapeFunctionsThirdDerivatives(int order,  Vector N, Matrix *rResult, double t ) const
    virtual Matrix& ShapeFunctionsThirdDerivatives(int order,
                                                   Vector N,
                                                   Matrix *rResult,
                                                   double t ) const 
    {
            KRATOS_ERROR << "Nurbs_2d::ShapeFunctionsThirdDerivatives" << "Third order derivatives not yet implemented" << std::endl;
            return *rResult;
    }







        /**
     * DefineGeometries: In this function the both knot vectors in Xi and Eta direction
     * will be analysed. For each knot span combination (Xi x Eta) there will be generated
     * one new geometry of the type NurbsPatchGeometry2D. Each geometry will contain
     * exactly the information which is needed to compute the shape functions and its
     * derivatives in a knot span. Furthermore the upper and lower limits for the Xi- and
     * Eta-values are stored.
     *
     * @param GeometryContainer: a container which holds all created geometries.
     */
    void DefineGeometries(std::vector< Geometry< Node<3> >::Pointer>& GeometryContainer)
    {
        /**
         * loop to count all knot spans (Xi- and Eta-direction)
         */

        int NumberOfSpansXi(0), NumberOfSpansEta(0), IndexU(0), IndexV(0);
		Vector KnotsXi(2 * mPolynomialDegreeP);
		Vector KnotsEta(2 * mPolynomialDegreeQ);
		Vector VectorKnotsXi(2*mPolynomialDegreeP), VectorKnotsEta(2*mPolynomialDegreeQ);
        double xi,eta;

        for (int i = 1; i <= mKnotsXi.size() - 2*mPolynomialDegreeP; i++)
        {
            if (mKnotsXi[mPolynomialDegreeP+i] != mKnotsXi[mPolynomialDegreeP+i-1])
            {
                NumberOfSpansXi++;
            }
        }

        for (int i = 1; i <= mKnotsEta.size() - 2*mPolynomialDegreeQ; i++)
        {
            if (mKnotsEta[mPolynomialDegreeQ+i] != mKnotsEta[mPolynomialDegreeQ+i-1])
            {
                NumberOfSpansEta++;
            }
        }


        mNumberOfGeometries = NumberOfSpansXi*NumberOfSpansEta;
        GeometryContainer.resize(mNumberOfGeometries);
        Vector WeightsElement((mPolynomialDegreeP+1)*(mPolynomialDegreeQ+1));
        int counter = 0;


        for (int i = 1; i <= mKnotsXi.size() - 2*mPolynomialDegreeP; i++)
        {
            if (mKnotsXi[mPolynomialDegreeP+i] != mKnotsXi[mPolynomialDegreeP+i-1])
            {
                for (int j = 1; j <= mKnotsEta.size() - 2*mPolynomialDegreeQ; j++)
                {
                    if (mKnotsEta[mPolynomialDegreeQ+j] != mKnotsEta[mPolynomialDegreeQ+j-1])
                    {
                        mLowerXi = mKnotsXi(mPolynomialDegreeP+i-1);
                        mUpperXi = mKnotsXi(mPolynomialDegreeP+i);
                        mLowerEta = mKnotsEta(mPolynomialDegreeQ+j-1);
                        mUpperEta = mKnotsEta(mPolynomialDegreeQ+j);
                        xi = (mKnotsXi(mPolynomialDegreeP+i)+mKnotsXi(mPolynomialDegreeP+i-1))/2;
                        eta = (mKnotsEta(mPolynomialDegreeQ+j)+mKnotsEta(mPolynomialDegreeQ+j-1))/2;
                        IndexU = FindKnotsToEvaluate(mPolynomialDegreeP,xi,1,&KnotsXi(0))-mPolynomialDegreeP;
                        IndexV = FindKnotsToEvaluate(mPolynomialDegreeQ,eta,2,&KnotsEta(0))-mPolynomialDegreeQ;
                        PointerVector< Node<3> > ControlPointsElement;
                        for(int k=0; k <= mPolynomialDegreeQ; k++)
                        {
                            for(int l=0; l <= mPolynomialDegreeP; l++)
                            {
                                ControlPointsElement.push_back(this->Points()((IndexV+k)*(mNumberOfCPsU)+(IndexU+l)));
                                WeightsElement[k *(mPolynomialDegreeP+1) + l] = mWeights[(IndexV+k)*(mNumberOfCPsU)+(IndexU+l)];
								
                            }
                        }
                        GeometryContainer[counter] = Geometry< Node<3> >::Pointer( new NurbsPatchGeometry2D(ControlPointsElement,
                                                                                                            WeightsElement,
																											KnotsXi,
                                                                                                            KnotsEta,
                                                                                                            mPolynomialDegreeP,
                                                                                                            mPolynomialDegreeQ,
                                                                                                            mPolynomialDegreeP+1,
                                                                                                            mPolynomialDegreeQ+1,
                                                                                                            mKnotsXi(mPolynomialDegreeP+i-1),
                                                                                                            mKnotsXi(mPolynomialDegreeP+i),
                                                                                                            mKnotsEta(mPolynomialDegreeQ+j-1),
                                                                                                            mKnotsEta(mPolynomialDegreeQ+j)) );
                        counter++;
                    }
                }
            }
        }
    }




    int GeometryNumber()
    {
    return mNumberOfGeometries;
    }


    /**
     * FindGeometryId: This function will provide the id of the geometry corresponding to an arbitrary
     * Xi and Eta combination.
     *
     * @param Xi: Value of the Xi- coordinate for which the right geometry is searched
     * @param Eta: Value of the Eta- coordinate for which the right geometry is searched
     * @return: Integer which gives the correct Id to the searched geometry.
     */





    int FindGeometryId(double Xi, double Eta)
    {
        int NumberOfSpansXi(0), ElementId(0), IdOfSpansU(0), IdOfSpansV(0), counter(0);

        for (int i = 0; i <= mKnotsXi.size() - (2*mPolynomialDegreeP+1); i++)
        {
            if (mKnotsXi[(mPolynomialDegreeP+1)+i] != mKnotsXi[(mPolynomialDegreeP+1)+i-1])
            {
                NumberOfSpansXi++;
                if (Xi <= mKnotsXi[(mPolynomialDegreeP+1)+i])
                {
                   if (counter == 0)
                   {
                      IdOfSpansU = NumberOfSpansXi;
                   }

                   counter = 1;
                }
            }
        }

        for (int i = 0; i <= mKnotsEta.size() - (2*mPolynomialDegreeQ+1); i++)
        {
            if (mKnotsEta[(mPolynomialDegreeQ+1)+i] != mKnotsEta[(mPolynomialDegreeQ+1)+i-1])
            {
                IdOfSpansV++;
                if (Eta <= mKnotsEta[(mPolynomialDegreeQ+1)+i])
                {
                    break;
                }
            }
        }


        ElementId = IdOfSpansU + (IdOfSpansV-1) * NumberOfSpansXi;
        return ElementId;

    }



    /**
     * @brief GeometryConvertToGlobalCoordinates will return the global coordinates
     * of a point in the parameter space
     * @param LocalCoordinates: Provides the local coordinates of the point
     * @param GlobalCoordinates: Provides the global coordinates of the point
     * @return: GlobalCoordinates
     */
    CoordinatesArrayType& GeometryConvertToGlobalCoordinates(CoordinatesArrayType &LocalCoordinates,
                                                 CoordinatesArrayType &GlobalCoordinates)
    {
        Matrix NurbsBasisFunction(mPolynomialDegreeP+1,mPolynomialDegreeQ+1);
        NurbsBasisFunction = ZeroMatrix(mPolynomialDegreeP+1,mPolynomialDegreeQ+1);
        NurbsBasisFunction = ElementGeometryNurbsFunctionsValues(LocalCoordinates[0],LocalCoordinates[1]);

        double x(0),y(0);

        for (int i=0; i<mPolynomialDegreeP+1; i++)
        {
            for (int j=0; j<mPolynomialDegreeQ+1; j++)
            {
                x += NurbsBasisFunction(i,j)* this->Points()[i+mNumberOfCPsU*j].X();
                y += NurbsBasisFunction(i,j)* this->Points()[i+mNumberOfCPsU*j].Y();

            }
        }
        GlobalCoordinates[0] = x;
        GlobalCoordinates[1] = y;


        return GlobalCoordinates;
    }


    /**
     * @brief BaseVectors: Will compute the base vectors of a NURBS-surface at a given
     * point in local coordinates.
     *
     * @param Xi: Local Xi-coordinate
     * @param Eta: Local Eta-coordinate
     * @param gXi: Vector for the base vector in Xi-direction
     * @param gEta: Vector for the base vector in Eta-direction
     */

    void BaseVectors(double Xi, double Eta, Vector &gXi, Vector &gEta) const
    {
        gXi.resize(3);
        gXi.resize(3,false);
        gXi = ZeroVector(3);
        gEta.resize(3);
        gEta.resize(3,false);
        gEta = ZeroVector(3);

        Matrix XiDerivatives, EtaDerivatives;
        ElementGeometryDerivatives(Xi,Eta,XiDerivatives,EtaDerivatives);

        for (int i=0;
             i<mPolynomialDegreeP+1;
             i++)
        {
            for(int j=0;
                j<mPolynomialDegreeQ+1;
                j++)
            {
                gXi[0] += XiDerivatives(i,j)*this->Points()[i+(mPolynomialDegreeP+1)*j].X();
                gXi[1] += XiDerivatives(i,j)*this->Points()[i+j*(mPolynomialDegreeP+1)].Y();
                gEta[0] += EtaDerivatives(i,j)*this->Points()[i+j*(mPolynomialDegreeP+1)].X();
                gEta[1] += EtaDerivatives(i,j)*this->Points()[i+j*(mPolynomialDegreeP+1)].Y();

            }

        }
    }


    /**
     * @brief NormalVectorToNurbsSurface: Will calculate the normal vector of a NURBS-surface
     * at a given point in local coordinates
     *
     * @param Xi: Local Xi-coordinate
     * @param Eta: Local Eta-coordinate
     * @return: 3D-Vector which contains the normal vector
     */

    Vector NormalVectorToNurbsSurface(double Xi, double Eta)const
    {
        Vector gXi,gEta;
        BaseVectors(Xi,Eta,gXi,gEta);
        Vector Normal = MathUtils<double>::UnitCrossProduct(gXi,gEta);
        return Normal;
    }




    /**
     * PointLocalCoordinates will return the closest Point (in local coordinates) from an provided
     * Point to the element geometry. If the Point does not lie on the geometry, it will check the
     * boundary values of the Xi-/Eta- span of the element and assign the largest, respectively smallest
     * possible Xi-/Eta- values to the return value. This assigning of values depends on the question if
     * the computed Xi-/Eta- values are smaller than the lower boundary of the Xi-/Eta- span or larger
     * than the upper boundary of the Xi-/Eta- span.
     *
     * @param rPoint contains the global coordinates for which the closest point on the surface is searched.
     *
     * @param rResult provides the closest point in local coordinates and also contains a check, if the
     * point lies inside or outside the geometry (rResult[2] = 0 means inside, rResult[2] = 1 means outside
     * of the geometry).
     *
     * @return rResult
     */

    virtual CoordinatesArrayType& PointLocalCoordinates( CoordinatesArrayType& rResult,
            const CoordinatesArrayType& rPoint )
    {

        boost::numeric::ublas::bounded_matrix<double,3,2> DN;

        double tol = 1.0e-8;
        int maxiter = 1000;

        Matrix J = ZeroMatrix( 2, 2 );
        Matrix invJ = ZeroMatrix( 2, 2 );

        //starting with xi in the middle of the element
        rResult = ZeroVector( 3 );
        rResult[0] = (mUpperXi-mLowerXi)/2.0;
        rResult[1] = (mUpperEta-mLowerEta)/2.0;
        Vector DeltaXi = ZeroVector( 2 );
        array_1d<double,3> CurrentGlobalCoords;


        //Newton iteration:
        for ( int k = 0; k < maxiter; k++ )
        {
            noalias(CurrentGlobalCoords) = ZeroVector( 3 );
            CurrentGlobalCoords = GeometryConvertToGlobalCoordinates(rResult,CurrentGlobalCoords);
            noalias( CurrentGlobalCoords ) = rPoint - CurrentGlobalCoords;

            //DN is the Transpose of the Jacobian and the Jacobian itself is provided by the function ElementGeometryDerivatives()
            Matrix Jacobian = ElementGeometryDerivatives(rResult[0],rResult[1]);
            DN = trans(Jacobian);
            noalias(J) = prod(trans(DN),DN);
            Vector res = prod(trans(DN),CurrentGlobalCoords);

            //deteminant of Jacobian
            const double det_j = J( 0, 0 ) * J( 1, 1 ) - J( 0, 1 ) * J( 1, 0 );

            //filling matrix
            invJ( 0, 0 ) = ( J( 1, 1 ) ) / ( det_j );
            invJ( 1, 0 ) = -( J( 1, 0 ) ) / ( det_j );
            invJ( 0, 1 ) = -( J( 0, 1 ) ) / ( det_j );
            invJ( 1, 1 ) = ( J( 0, 0 ) ) / ( det_j );

            //computing the change in local coordinates
            DeltaXi( 0 ) = invJ( 0, 0 ) * res[0] + invJ( 0, 1 ) * res[1];
            DeltaXi( 1 ) = invJ( 1, 0 ) * res[0] + invJ( 1, 1 ) * res[1];

            rResult[0] += DeltaXi[0];
            rResult[1] += DeltaXi[1];
            rResult[2] = 0.0;

            if ( norm_2( DeltaXi ) > 300 )
            {
                res[0] = 0.0;
                res[1] = 0.0;
                res[2] = 0.0;
                std::cout << "detJ =" << det_j << "DeltaX = " << DeltaXi << " stopping calculation and assigning the baricenter" << std::endl;
                break;
                //KRATOS_ERROR << "Computation of local coordinates failed at iteration" << k << std::endl;
            }

            if ( norm_2( DeltaXi ) < tol )
            {
                break;
            }
            KRATOS_WATCH(rResult);
        }

         CurrentGlobalCoords = GeometryConvertToGlobalCoordinates(rResult,CurrentGlobalCoords);
         KRATOS_WATCH(CurrentGlobalCoords);

         //Check if Xi/Eta lie inside of the element boundaries and assigning the closest possible local coordinates of the element
         //rResult[2] is indicator if or not the Point lies inside of the element boundaries (0 = yes, 1 = no)

         if(rResult[0] < mLowerXi)
         {
             rResult[2]=1;
             rResult[0] = mLowerXi;
         }
         else if(rResult[0] > mUpperXi)
         {
             rResult[2]=1;
             rResult[0] = mUpperXi;
         }
         else if(rResult[1] < mLowerEta )
         {
             rResult[2]=1;
             rResult[1] = mLowerEta;
         }
         else if(rResult[1] > mUpperEta)
         {
             rResult[2]=1;
             rResult[1] = mUpperEta;
         }


        return( rResult );
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    /**
     * There are no protected members in class NurbsPatchGeometry2D
     */

private:
    ///@name Static Member Variables
    ///@{
    static const GeometryData msGeometryData;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save( Serializer& rSerializer ) const
    {

        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, PointsArrayType );
    }

    virtual void load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, PointsArrayType );
    }

    //NurbsPatchGeometry2D(): BaseType( PointsArrayType(), &msGeometryData ) {}

    ///@}
    ///@name Member Variables
    ///@{

    Vector mKnotsXi;
    Vector mKnotsEta;
    Vector mWeights;
    int mPolynomialDegreeP,mPolynomialDegreeQ,mNumberOfGeometries;
    int mNumberOfCPsU, mNumberOfCPsV;
    double mUpperXi,mLowerXi,mUpperEta,mLowerEta;
    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    /**
     * @brief FindKnotSpan returns the id of the knot which is the lower knot value of the knot span in which the searched Xi/Eta
     * value lies in.
     * @param n: n = m - p - 1 with m = size of the knot vector; p = polynomial degree
     * @param p: polynomial degree
     * @param u: Xi/Eta Value for which the id of the knot span shall be deduced
     * @param Knots: Regarded Knot Vector (either Xi- or Eta-direction)
     * @return int which deduces the id of the knot span
     */
    int FindKnotSpan(int n,
                     int p,
                     double u,
                     const Vector &Knots)const
    {

        {
            if (u == Knots[n+1]) return(n);
            int low(p),high(n+1);
            int mid = low+high/2;
            while  (u<Knots(mid) || u >= Knots(mid+1))
            {
                if(u<Knots(mid)) high = mid;
                else low = mid;
                mid = (low + high)/2;
            }
            return(mid);
        }

    }


    /**
     * @brief FindKnotsToEvaluate: Function to get all Knots which are needed for calculation of the shape functions
     * at a specific point in the parameter space
     * @param d: polynomial degree
     * @param t: Xi/Eta- Position for which the reduced Knot Vector is needed
     * @param dimension: either 1 (Xi-direction) or 2 (Eta-direction)
     * @param knots: c-array which after calling the function contains the desired knot values
     * @return the index of the knot which is the lower knot value of the knot span in which the searched Xi/Eta
     * value lies in.
     */
    int FindKnotsToEvaluate(int d,
                            double t,
                            int dimension,
                            double *knots)const
    {
        int position(0);
        switch(dimension)
        {
        case 1:
            for (int i=d;i != (mKnotsXi.size()-d);i++)
            {
                if(t < mKnotsXi[i] )
                {
                    position = i-1;

                    std::cout << "The Element with the following knot vector has been evaluated: (u-Direction)" << std::endl <<"[ ";
                    for (int j=0 ;j != (2*d);j++)
                    {

                        knots[j] = mKnotsXi(i-d+j);
                        std::cout << knots[j] << " ";
                    }
                    std::cout << "] at Xi = "<< t << std::endl;
                    i=mKnotsXi.size()-d-1; //exit loop
                    return position;
                }

            }
            break;

        case 2:
            for (int i=d;i != (mKnotsEta.size()-d);i++)
            {
                if(t < mKnotsEta[i] )
                {
                    position = i-1;
                    std::cout << "The Element with the following knot vector has been evaluated: (v-Direction)" << std::endl <<"[ ";
                    for (int j=0 ;j != (2*d);j++)
                    {

                        knots[j] = mKnotsEta(i-d+j);
                        std::cout << knots[j] << " ";
                    }
                    std::cout << "] at Eta = "<< t << std::endl;
                    i=mKnotsEta.size()-d-1; //exit loop
                    return position;
                }
            }
            break;

        default:
            std::cout<<"No or wrong dimension submitted to the function \"FindKnotsToEvaluate\" "<<std::endl;
        }
        return 0;

    }



    /**
     * ElementGeometryDerivatives calculates the local derivatives of x and y
     * in the ElementGeometry. Therefore it needs the already reduced ElementGeometry
     * obtained by calling the function DefineGeometries() and a point in the parameter space
     * described by its local coordinates Xi and Eta.
     * The returned matrix F_i_j will have the dimension (p+1) x (q+1), where
     * p:= Polynomial Degree in Xi-direction and q:= Polynomial Degree in Eta-direction.
     * Therefore the calculated Matrix holds the following data:
     * F_0_0 = dx/dXi
     * F_0_1 = dx/dEta
     * F_1_0 = dy/dXi
     * F_1_1 = dy/dEta
     *
     * @param Xi: Local Xi- Coordinate where the derivatives shall be evaluated
     * @param Eta: Local Eta- Coordinate where the derivatives shall be evaluated
     * @param NurbsBasisFunctionDerivativesXi: will contain the NURBS-Basis Functions Xi-Derivatives
     * @param NurbsBasisFunctionDerivativesEta: will contain the NURBS-Basis Functions Eta-Derivatives
     * @param NurbsFunctionsValues: will contain the NURBS-Basis Functions Values
     * @return: Local derivatives of x and y.
     *
     * @note: The NurbsBasisFunctionDerivativesXi/Eta R_i_j are the NURBS local derivatives.
     *          i = Control Point index in Xi-direction
     *          j = Control Point index in Eta-direction
     *
     */
	Matrix ElementGeometryDerivatives(double Xi,
		double Eta,
		Matrix &NurbsBasisFunctionDerivativesXi,
		Matrix &NurbsBasisFunctionDerivativesEta,
		Vector &NurbsFunctionsValues)const
	{
		//std::cout<<"ElementGeometryDerivatives "<<std::endl;
		NurbsBasisFunctionDerivativesXi.resize(mPolynomialDegreeP + 1, mPolynomialDegreeQ + 1);
		NurbsBasisFunctionDerivativesEta.resize(mPolynomialDegreeP + 1, mPolynomialDegreeQ + 1);
		NurbsFunctionsValues.resize((mPolynomialDegreeP + 1)*(mPolynomialDegreeQ + 1));
		NurbsBasisFunctionDerivativesXi = ZeroMatrix(mPolynomialDegreeP + 1, mPolynomialDegreeQ + 1);
		NurbsBasisFunctionDerivativesEta = ZeroMatrix(mPolynomialDegreeP + 1, mPolynomialDegreeQ + 1);
		NurbsFunctionsValues = ZeroVector((mPolynomialDegreeP + 1)*(mPolynomialDegreeQ + 1));
		int der_count(1);
		double SumDenominatorXiEta(0), SumNumeratorXi(0), SumNumeratorEta(0);
		Matrix dN_Xi(mPolynomialDegreeP + 1,mPolynomialDegreeP + 1);
		Matrix dN_Eta(mPolynomialDegreeQ + 1,mPolynomialDegreeQ + 1);

		// Call the function of open NURBS to calculate Bspline shape functions;
		ON_ShapeFunctionValue(mPolynomialDegreeP + 1, &mKnotsXi[0], Xi, &dN_Xi(0,0));
		ON_ShapeFunctionValue(mPolynomialDegreeQ + 1, &mKnotsEta[0], Eta, &dN_Eta(0,0));
		Vector N_Xi(mPolynomialDegreeP + 1), N_Eta(mPolynomialDegreeQ + 1);

		for (int i = 0; i<mPolynomialDegreeP + 1; i++)
		{
			N_Xi[i] = dN_Xi(0,i);
		}

		for (int j = 0; j<mPolynomialDegreeQ + 1; j++)
		{
			N_Eta[j] = dN_Eta(0,j);
		}

		Matrix ShapeFunctionsValuesMatrixForm(mPolynomialDegreeP + 1, mPolynomialDegreeQ + 1);

		// Modify Bsplines to NURBS taking into account the weights;
		ShapeFunctionsValuesMatrixForm = BSplinesToNurbs(N_Xi, N_Eta);


		for (int j = 0; j<mPolynomialDegreeQ + 1; j++)
		{
			for (int i = 0; i< mPolynomialDegreeP + 1; i++)
			{
				NurbsFunctionsValues[j*(mPolynomialDegreeP + 1) + i] = ShapeFunctionsValuesMatrixForm(i, j);
			}
		}

		ON_EvaluateNurbsBasisDerivatives(mPolynomialDegreeQ + 1, &mKnotsEta[0], der_count, &dN_Eta(0,0));
		ON_EvaluateNurbsBasisDerivatives(mPolynomialDegreeP + 1, &mKnotsXi[0], der_count, &dN_Xi(0,0));
		Matrix Jacobian2D(2, 2);
		Jacobian2D = ZeroMatrix(2, 2);

		for (int i = 0; i<mPolynomialDegreeP + 1; i++)
		{
			for (int j = 0; j<mPolynomialDegreeQ + 1; j++)
			{
				SumDenominatorXiEta += dN_Xi(0,i) * dN_Eta(0,j) * mWeights[i + mNumberOfCPsU*j];
				SumNumeratorXi += dN_Xi(1,i) * dN_Eta(0,j) * mWeights[i + mNumberOfCPsU*j];
				SumNumeratorEta += dN_Xi(0,i) * dN_Eta(1,j) * mWeights[i + mNumberOfCPsU*j];
			}
		}

		double invSumDenominatorXiEta(1 / (SumDenominatorXiEta*SumDenominatorXiEta));

		for (int j = 0; j<mPolynomialDegreeQ + 1; j++)
		{
			for (int i = 0; i<mPolynomialDegreeP + 1; i++)
			{
				NurbsBasisFunctionDerivativesXi(i, j) = (SumDenominatorXiEta * dN_Xi(1,i) * dN_Eta(0,j) * mWeights[i + mNumberOfCPsU*j] - dN_Xi(0,i) * dN_Eta(0,j) * mWeights[i + mNumberOfCPsU*j] * SumNumeratorXi) * invSumDenominatorXiEta;
				NurbsBasisFunctionDerivativesEta(i, j) = (SumDenominatorXiEta * dN_Xi(0,i) * dN_Eta(1,j) * mWeights[i + mNumberOfCPsU*j] - dN_Xi(0,i) * dN_Eta(0,j) * mWeights[i + mNumberOfCPsU*j] * SumNumeratorEta) * invSumDenominatorXiEta;
				Jacobian2D(0, 0) += NurbsBasisFunctionDerivativesXi(i, j)* this->Points()[i + (mPolynomialDegreeP + 1)*j].X();
				Jacobian2D(0, 1) += NurbsBasisFunctionDerivativesXi(i, j)* this->Points()[i + (mPolynomialDegreeP + 1)*j].Y();
				Jacobian2D(1, 0) += NurbsBasisFunctionDerivativesEta(i, j)* this->Points()[i + (mPolynomialDegreeP + 1)*j].X();
				Jacobian2D(1, 1) += NurbsBasisFunctionDerivativesEta(i, j)* this->Points()[i + (mPolynomialDegreeP + 1)*j].Y();
			}
		}

		return Jacobian2D;
	}


	void ElementGeometrySecondDerivatives(
		double Xi,
		double Eta,
		Matrix &NurbsBasisFunctionSecondDerivativesXi,
		Matrix &NurbsBasisFunctionSecondDerivativesEta,
		Matrix &NurbsBasisFunctionSecondDerivativesXiEta,
		Matrix &NurbsBasisFunctionDerivativesXi,
		Matrix &NurbsBasisFunctionDerivativesEta,
		Vector &NurbsFunctionsValues)const
	{
		NurbsBasisFunctionSecondDerivativesXi.resize(mPolynomialDegreeP + 1, mPolynomialDegreeQ + 1);
		NurbsBasisFunctionSecondDerivativesEta.resize(mPolynomialDegreeP + 1, mPolynomialDegreeQ + 1);
		NurbsBasisFunctionSecondDerivativesXiEta.resize((mPolynomialDegreeP + 1),(mPolynomialDegreeQ + 1));
		NurbsBasisFunctionDerivativesXi.resize(mPolynomialDegreeP + 1, mPolynomialDegreeQ + 1);
		NurbsBasisFunctionDerivativesEta.resize(mPolynomialDegreeP + 1, mPolynomialDegreeQ + 1);
		NurbsFunctionsValues.resize((mPolynomialDegreeP + 1)*(mPolynomialDegreeQ + 1));

		NurbsBasisFunctionSecondDerivativesXi = ZeroMatrix(mPolynomialDegreeP + 1, mPolynomialDegreeQ + 1);
		NurbsBasisFunctionSecondDerivativesEta = ZeroMatrix(mPolynomialDegreeP + 1, mPolynomialDegreeQ + 1);
		NurbsBasisFunctionSecondDerivativesXiEta = ZeroMatrix((mPolynomialDegreeP + 1),(mPolynomialDegreeQ + 1));
		NurbsBasisFunctionDerivativesXi = ZeroMatrix(mPolynomialDegreeP + 1, mPolynomialDegreeQ + 1);
		NurbsBasisFunctionDerivativesEta = ZeroMatrix(mPolynomialDegreeP + 1, mPolynomialDegreeQ + 1);
		NurbsFunctionsValues = ZeroVector((mPolynomialDegreeP + 1)*(mPolynomialDegreeQ + 1));

		int der_count = 2;
		double SumDenominatorXiEta(0), SumNumeratorXi(0), SumNumeratorEta(0);
		double SumDenominatorSecondXiEta(0), SumNumeratorSecondXi(0), SumNumeratorSecondEta(0); //, SumNumeratorSecondXiEta(0);

		Matrix dN_Xi = ZeroMatrix(mPolynomialDegreeP + 2, mPolynomialDegreeP + 1);
		Matrix dN_Eta = ZeroMatrix(mPolynomialDegreeQ + 2, mPolynomialDegreeQ + 1);

		// Call the function of open NURBS to calculate Bspline shape functions;
		ON_ShapeFunctionValue(mPolynomialDegreeP + 1, &mKnotsXi[0], Xi, &dN_Xi(0, 0));
		ON_ShapeFunctionValue(mPolynomialDegreeQ + 1, &mKnotsEta[0], Eta, &dN_Eta(0, 0));
		Vector N_Xi = ZeroVector(mPolynomialDegreeP + 1);
		Vector N_Eta = ZeroVector(mPolynomialDegreeQ + 1);

		for (int i = 0; i<mPolynomialDegreeP + 1; i++)
		{
			N_Xi[i] = dN_Xi(0, i);
		}

		for (int j = 0; j<mPolynomialDegreeQ + 1; j++)
		{
			N_Eta[j] = dN_Eta(0, j);
		}

		Matrix ShapeFunctionsValuesMatrixForm = ZeroMatrix(mPolynomialDegreeP + 1, mPolynomialDegreeQ + 1);

		// Modify Bsplines to NURBS taking into account the weights;
		ShapeFunctionsValuesMatrixForm = BSplinesToNurbs(N_Xi, N_Eta);

		for (int j = 0; j<mPolynomialDegreeQ + 1; j++)
		{
			for (int i = 0; i< mPolynomialDegreeP + 1; i++)
			{
				NurbsFunctionsValues[j*(mPolynomialDegreeP + 1) + i] = ShapeFunctionsValuesMatrixForm(i, j);
			}
		}

		ON_EvaluateNurbsBasisDerivatives(mPolynomialDegreeQ + 1, &mKnotsEta[0], der_count, &dN_Eta(0, 0));
		ON_EvaluateNurbsBasisDerivatives(mPolynomialDegreeP + 1, &mKnotsXi[0], der_count, &dN_Xi(0, 0));
		//Matrix Jacobian2D(2, 2);
		//Jacobian2D = ZeroMatrix(2, 2);

		//KRATOS_WATCH(dN_Eta)
		//KRATOS_WATCH(dN_Xi)

		//KRATOS_WATCH(mWeights)
		for (int j = 0; j<mPolynomialDegreeQ + 1; j++)
		{
			for (int i = 0; i<mPolynomialDegreeP + 1; i++)
			{
				//std::cout <<  mWeights[i + (mPolynomialDegreeP + 1)*j] << std::endl;

				SumDenominatorXiEta += dN_Xi(0, i) * dN_Eta(0, j) * mWeights[i + (mPolynomialDegreeP + 1)*j];
				//std::cout << dN_Xi(0, i) * dN_Eta(0, j) * mWeights[i + (mPolynomialDegreeP + 1)*j] << std::endl;

				SumNumeratorXi += dN_Xi(1, i) * dN_Eta(0, j) * mWeights[i + (mPolynomialDegreeP + 1)*j];
				SumNumeratorEta += dN_Xi(0, i) * dN_Eta(1, j) * mWeights[i + (mPolynomialDegreeP + 1)*j];

				//std::cout << dN_Xi(1, i) * dN_Eta(0, j) * mWeights[i + (mPolynomialDegreeP + 1)*j] << std::endl;
				//std::cout << dN_Xi(0, i) * dN_Eta(1, j) * mWeights[i + (mPolynomialDegreeP + 1)*j] << std::endl;


				SumNumeratorSecondXi += dN_Xi(2, i)*dN_Eta(0, j)*mWeights[i + (mPolynomialDegreeP + 1)*j];
				SumNumeratorSecondEta += dN_Xi(0, i)*dN_Eta(2, j)*mWeights[i + (mPolynomialDegreeP + 1)*j];
				SumDenominatorSecondXiEta += dN_Xi(1, i)*dN_Eta(1, j)*mWeights[i + (mPolynomialDegreeP + 1)*j];
			}
		}
		
		double invSumDenominator = 1/ SumDenominatorXiEta;
		double invSumDenominatorXiEta(1 / (SumDenominatorXiEta*SumDenominatorXiEta));
		double invSumDenominatorXiEta2 = 1 / (SumDenominatorXiEta*SumDenominatorXiEta*SumDenominatorXiEta);

		for (int j = 0; j<mPolynomialDegreeQ + 1; j++)
		{
			for (int i = 0; i<mPolynomialDegreeP + 1; i++)
			{
				NurbsBasisFunctionDerivativesXi(i, j) = (SumDenominatorXiEta * dN_Xi(1, i) * dN_Eta(0, j) * mWeights[i + (mPolynomialDegreeP + 1)*j] 
					- dN_Xi(0, i) * dN_Eta(0, j) * mWeights[i + (mPolynomialDegreeP + 1)*j] * SumNumeratorXi) * invSumDenominatorXiEta;
				NurbsBasisFunctionDerivativesEta(i, j) = (SumDenominatorXiEta * dN_Xi(0, i) * dN_Eta(1, j) * mWeights[i + (mPolynomialDegreeP + 1)*j]
					- dN_Xi(0, i) * dN_Eta(0, j) * mWeights[i + (mPolynomialDegreeP + 1)*j] * SumNumeratorEta) * invSumDenominatorXiEta;

				NurbsBasisFunctionSecondDerivativesXi(i, j) = dN_Xi(2, i)*dN_Eta(0, j)*mWeights[i + (mPolynomialDegreeP + 1)*j] * invSumDenominator
					- 2.0*dN_Xi(1, i) * dN_Eta(0, j) * mWeights[i + (mPolynomialDegreeP + 1)*j] * SumNumeratorXi * invSumDenominatorXiEta
					- dN_Xi(0, i) * dN_Eta(0, j) * mWeights[i + (mPolynomialDegreeP + 1)*j] * SumDenominatorSecondXiEta * invSumDenominatorXiEta
					+ 2.0 * dN_Xi(0, i) * dN_Eta(0, j) * mWeights[i + (mPolynomialDegreeP + 1)*j] * SumNumeratorXi * SumNumeratorXi * invSumDenominatorXiEta2;


				NurbsBasisFunctionSecondDerivativesEta(i, j) = dN_Xi(0, i)*dN_Eta(2, j)*mWeights[i + (mPolynomialDegreeP + 1)*j] * invSumDenominator
					- 2.0*dN_Xi(0, i) * dN_Eta(1, j) * mWeights[i + (mPolynomialDegreeP + 1)*j] * SumNumeratorEta * invSumDenominatorXiEta
					- dN_Xi(0, i) * dN_Eta(0, j) * mWeights[i + (mPolynomialDegreeP + 1)*j] * SumNumeratorSecondEta * invSumDenominatorXiEta 
					+ 2.0*dN_Xi(0, i) * dN_Eta(0, j) * mWeights[i + (mPolynomialDegreeP + 1)*j] * SumNumeratorEta * SumNumeratorEta * invSumDenominatorXiEta2;

				NurbsBasisFunctionSecondDerivativesXiEta(i, j) = dN_Xi(1, i)*dN_Eta(1, j)*mWeights[i + (mPolynomialDegreeP + 1)*j] * invSumDenominator
					- dN_Xi(1, i) * dN_Eta(0, j) * mWeights[i + (mPolynomialDegreeP + 1)*j] *SumNumeratorEta * invSumDenominatorXiEta
					- dN_Xi(0, i) * dN_Eta(1, j) * mWeights[i + (mPolynomialDegreeP + 1)*j] *SumNumeratorXi * invSumDenominatorXiEta
					- dN_Xi(0, i) * dN_Eta(0, j) * mWeights[i + (mPolynomialDegreeP + 1)*j] * SumDenominatorSecondXiEta * invSumDenominatorXiEta
					+ 2.0* dN_Xi(0, i) * dN_Eta(0, j) * mWeights[i + (mPolynomialDegreeP + 1)*j] * SumNumeratorXi * SumNumeratorEta * invSumDenominatorXiEta2;

				//Jacobian2D(0, 0) += NurbsBasisFunctionDerivativesXi(i, j)* this->Points()[i + (mPolynomialDegreeP + 1)*j].X();
				//Jacobian2D(0, 1) += NurbsBasisFunctionDerivativesXi(i, j)* this->Points()[i + (mPolynomialDegreeP + 1)*j].Y();
				//Jacobian2D(1, 0) += NurbsBasisFunctionDerivativesEta(i, j)* this->Points()[i + (mPolynomialDegreeP + 1)*j].X();
				//Jacobian2D(1, 1) += NurbsBasisFunctionDerivativesEta(i, j)* this->Points()[i + (mPolynomialDegreeP + 1)*j].Y();
			}
		}
		//KRATOS_WATCH(dN_Eta)
		//KRATOS_WATCH(dN_Xi)
		//KRATOS_WATCH(NurbsBasisFunctionDerivativesXi)
		//KRATOS_WATCH(NurbsBasisFunctionDerivativesEta)
		//KRATOS_WATCH(NurbsBasisFunctionSecondDerivativesXi)
		//KRATOS_WATCH(NurbsBasisFunctionSecondDerivativesEta)
		//KRATOS_WATCH(NurbsBasisFunctionSecondDerivativesXiEta)
		//return Jacobian2D;
	}



    /**
     * ElementGeometryDerivatives calculates the local derivatives of x and y
     * in the ElementGeometry. Therefore it needs the already reduced ElementGeometry
     * obtained by calling the function DefineGeometries() and a point in the parameter space
     * described by its local coordinates Xi and Eta.
     * The returned matrix F_i_j will have the dimension (p+1) x (q+1), where
     * p:= Polynomial Degree in Xi-direction and q:= Polynomial Degree in Eta-direction.
     * Therefore the calculated Matrix holds the following data:
     * F_0_0 = dx/dXi
     * F_0_1 = dx/dEta
     * F_1_0 = dy/dXi
     * F_1_1 = dy/dEta
     *
     * @param Xi: Local Xi- Coordinate where the derivatives shall be evaluated
     * @param Eta: Local Eta- Coordinate where the derivatives shall be evaluated
     * @param NurbsBasisFunctionDerivativesXi: will contain the NURBS-Basis Functions Xi-Derivatives
     * @param NurbsBasisFunctionDerivativesEta: will contain the NURBS-Basis Functions Eta-Derivatives
     * @return: Local derivatives of x and y.
     *
     * @note: The NurbsBasisFunctionDerivativesXi/Eta R_i_j are the NURBS local derivatives.
     *          i = Control Point index in Xi-direction
     *          j = Control Point index in Eta-direction
     *
     */

    Matrix ElementGeometryDerivatives(double Xi,
                                      double Eta,
                                      Matrix &NurbsBasisFunctionDerivativesXi,
                                      Matrix &NurbsBasisFunctionDerivativesEta)const
    {
        Matrix Jacobian2D;
        Vector NurbsFunctionsValues;
        Jacobian2D = ElementGeometryDerivatives(Xi,Eta,NurbsBasisFunctionDerivativesXi,NurbsBasisFunctionDerivativesEta,NurbsFunctionsValues);

        return Jacobian2D;
    }





    /**
     * ElementGeometryDerivatives calculates the local derivatives of x and y
     * in the ElementGeometry. Therefore it needs the already reduced ElementGeometry
     * obtained by calling the function DefineGeometries() and a point in the parameter space
     * described by its local coordinates Xi and Eta.
     * The returned matrix F_i_j will have the dimension (p+1) x (q+1), where
     * p:= Polynomial Degree in Xi-direction and q:= Polynomial Degree in Eta-direction.
     * Therefore the calculated Matrix holds the following data:
     * F_0_0 = dx/dXi
     * F_0_1 = dx/dEta
     * F_1_0 = dy/dXi
     * F_1_1 = dy/dEta
     *
     * @param Xi: Local Xi- Coordinate where the derivatives shall be evaluated
     * @param Eta: Local Eta- Coordinate where the derivatives shall be evaluated
     * @return: Local derivatives of x and y.
     */

    Matrix ElementGeometryDerivatives(double Xi,
                                      double Eta)const
    {
        Matrix NurbsBasisFunctionDerivativesXi, NurbsBasisFunctionDerivativesEta, Jacobian2D;
        Vector NurbsFunctionsValues;
        Jacobian2D = ElementGeometryDerivatives(Xi,Eta,NurbsBasisFunctionDerivativesXi,NurbsBasisFunctionDerivativesEta,NurbsFunctionsValues);
        return Jacobian2D;
    }




	     /**
     * BSplinesToNurbs function takes as input the
     * evaluated B-Spline Shape Functions and modifies them to NURBS-
     * Basis Functions taking into account the weights
     * Thereby the NURBS-Functions are stored in a Matrix NurbsBasisFunction_i_j where
     *          i = Control Point index in Xi-direction
     *          j = Control Point index in Eta-direction
     *
     * @param N_Xi are the B-Spline Shape Functions Values in Xi-direction
     * @param N_Eta are the B-Spline Shape Functions Values in Eta-direction
     * @return NurbsBasisFunction
     */
    Matrix BSplinesToNurbs(Vector N_Xi,Vector N_Eta)const
    {

        //std::cout<<"BSplinesToNurbs (nurbs_2d.h)"<<std::endl;
        Matrix NurbsBasisFunction(mPolynomialDegreeP+1,mPolynomialDegreeQ+1);

        double sum(0);

        for (int i=0; i<mPolynomialDegreeP+1; i++)
        {
            for (int j=0; j<mPolynomialDegreeQ+1; j++)
            {
                sum += N_Xi[i]*N_Eta[j]*mWeights[i+mNumberOfCPsU*j];
                NurbsBasisFunction(i,j) = N_Xi[i]*N_Eta[j]*mWeights[i+mNumberOfCPsU*j];
            }
        }
        const double inv_sum = 1.0/sum;
        for (int i=0; i<mPolynomialDegreeP+1; i++)
        {
            for (int j=0; j<mPolynomialDegreeQ+1; j++)
            {
                NurbsBasisFunction(i,j) = NurbsBasisFunction(i,j)*inv_sum;
            }
        }

        return NurbsBasisFunction;
    }






    /**
     * ElementGeometryNurbsFunctionsValues here takes as input the OpenNURBS results from
     * evaluating the B-Spline Shape Functions and combines the B-Splines to NURBS-
     * Basis Functions.
     * Thereby the NURBS-Functions are stored in a Matrix F_i_j where
     *          i = Control Point index in Xi-direction
     *          j = Control Point index in Eta-direction
     *
     * @param N_Xi are the B-Spline Shape Functions Values in Xi-direction
     * @param N_Eta are the B-Spline Shape Functions Values in Eta-direction
     * @return F_i_j
     */
    Matrix ElementGeometryNurbsFunctionsValues(Vector N_Xi,
                                               Vector N_Eta)const
    {

        Matrix NurbsBasisFunction(mPolynomialDegreeP+1,mPolynomialDegreeQ+1);

        double sum(0);

        for (int i=0; i<mPolynomialDegreeP+1; i++)
        {
            for (int j=0; j<mPolynomialDegreeQ+1; j++)
            {
                sum += N_Xi[i]*N_Eta[j]*mWeights[i+mNumberOfCPsU*j];
                NurbsBasisFunction(i,j) = N_Xi[i]*N_Eta[j]*mWeights[i+mNumberOfCPsU*j];
            }
        }
        const double inv_sum = 1.0/sum;
        for (int i=0; i<mPolynomialDegreeP+1; i++)
        {
            for (int j=0; j<mPolynomialDegreeQ+1; j++)
            {
                    NurbsBasisFunction(i,j) = NurbsBasisFunction(i,j)*inv_sum;
            }
        }

        return NurbsBasisFunction;
    }




    /**
     * ElementGeometryNurbsFunctionsValues here takes as input a local point
     * in the parameter space and returns all Shape Functions Values (NURBS).
     * Thereby the NURBS-Functions are stored in a Matrix F_i_j where
     *          i = Control Point index in Xi-direction
     *          j = Control Point index in Eta-direction
     *
     * @param Xi is the Xi-Coordinate of the point to be evaluated in the parameter space
     * @param Eta is the Eta-Coordinate of the point to be evaluated in the parameter space
     * @return F_i_j
     */



    Matrix ElementGeometryNurbsFunctionsValues(const double Xi,const double Eta)const
    {
        Matrix N_Xi(mPolynomialDegreeP+1,mPolynomialDegreeP+1);
        Matrix N_Eta(mPolynomialDegreeQ+1,mPolynomialDegreeQ+1);
        ON_ShapeFunctionValue(mPolynomialDegreeP+1, &mKnotsXi[0], Xi, &N_Xi(0,0));
        ON_ShapeFunctionValue(mPolynomialDegreeQ+1, &mKnotsEta[0], Eta, &N_Eta(0,0));

        Matrix NurbsBasisFunction(mPolynomialDegreeP+1,mPolynomialDegreeQ+1);
        double sum(0);

        for (int i=0; i<mPolynomialDegreeP+1; i++)
        {
            for (int j=0; j<mPolynomialDegreeQ+1; j++)
            {
                sum += N_Xi(0,i)*N_Eta(0,j)*mWeights[i+mNumberOfCPsU*j];
                NurbsBasisFunction(i,j) = N_Xi(0,i)*N_Eta(0,j)*mWeights[i+mNumberOfCPsU*j];
            }
        }
        const double inv_sum = 1.0/sum;
        for (int i=0; i<mPolynomialDegreeP+1; i++)
        {
            for (int j=0; j<mPolynomialDegreeQ+1; j++)
            {
                    NurbsBasisFunction(i,j) = NurbsBasisFunction(i,j)*inv_sum;
            }
        }

        return NurbsBasisFunction;
    }







    /**
     * CalculateShapeFunctionsIntegrationPointsLocalGradients is a function to compute
     * the local derivatives dx/dXi, dx/dEta, dy/dXi and dy/dEta for each integration
     * point.
     *
     * @param ThisMethod is the integration method used
     * @param d_Xi_shape_f_values container of the Shape Functions derivatives dR/dXi
     *        where R are the NURBS-shape functions.
     * @param d_Eta_shape_f_values container of the Shape Functions derivatives dR/dEta
     *        where R are the NURBS-shape functions.
     * @param ShapeFunctionsValues container of the NURBS Shape functions values
     *
     * @return a Vector of Matrices. Each Matrix represents the local derivatives at one gauss point
     *  return[gausspoint](0,0) = dx/dXi
     *  return[gausspoint](0,1) = dx/dEta
     *  return[gausspoint](1,0) = dy/dXi
     *  return[gausspoint](1,1) = dy/dEta

     * @note ShapeFunctionsValues, d_Xi_shape_f_values and d_Eta_shape_f_values also contain information
     * at each gauss point.
     * Thereby the storage scheme is for all three of them like F_i_j where
     * i = Control Point index in Xi-direction
     * j = Control Point index in Eta-direction
     */


    ShapeFunctionsGradientsType
    CalculateShapeFunctionsIntegrationPointsLocalGradients(
            typename BaseType::IntegrationMethod ThisMethod,
            ShapeFunctionsGradientsType &d_Xi_shape_f_values,
            ShapeFunctionsGradientsType &d_Eta_shape_f_values,
            Matrix &ShapeFunctionsValues) const
    {
        IntegrationPointsContainerType all_integration_points =
            AllIntegrationPoints();
        IntegrationPointsArrayType integration_points =
            all_integration_points[ThisMethod];
        //number of integration points
        const int integration_points_number = integration_points.size();

        ShapeFunctionsGradientsType d_shape_f_values( integration_points_number );
        Vector ShapeFunctionsValuesAtGaussPoint(this->PointsNumber());
        double XiLocalCoordinates,EtaLocalCoordinates;

        for ( int pnt = 0; pnt < integration_points_number; pnt++ )
        {
            d_Xi_shape_f_values[pnt].resize( mPolynomialDegreeP+1, mPolynomialDegreeQ+1,false);
            d_Eta_shape_f_values[pnt].resize( mPolynomialDegreeP+1, mPolynomialDegreeQ+1,false);
            XiLocalCoordinates = (mUpperXi - mLowerXi) / 2 * (1+integration_points[pnt].X())+mLowerXi;
            EtaLocalCoordinates = (mUpperEta - mLowerEta) / 2 * (1+integration_points[pnt].Y())+mLowerEta;
            d_shape_f_values[pnt] = ElementGeometryDerivatives(XiLocalCoordinates,
                                                                        EtaLocalCoordinates,
                                                                        d_Xi_shape_f_values[pnt],
                                                                        d_Eta_shape_f_values[pnt],
                                                                        ShapeFunctionsValuesAtGaussPoint);
            for(int j=0;j<this->PointsNumber(); j++)
            {
                    ShapeFunctionsValues(pnt,j) =  ShapeFunctionsValuesAtGaussPoint[j];
            }

        }

        return d_shape_f_values;
    }




    /**
     * CalculateShapeFunctionsIntegrationPointsLocalGradients is a function to compute
     * the local derivatives dx/dXi, dx/dEta, dy/dXi and dy/dEta for each integration
     * point.
     *
     * @param ThisMethod is the integration method used
     * @param d_Xi_shape_f_values container of the Shape Functions derivatives dR/dXi
     *        where R are the NURBS-shape functions.
     * @param d_Eta_shape_f_values container of the Shape Functions derivatives dR/dEta
     *        where R are the NURBS-shape functions.
     *
     * @return a Vector of Matrices. Each Matrix represents the local derivatives at one gauss point
     *  return[gausspoint](0,0) = dx/dXi
     *  return[gausspoint](0,1) = dx/dEta
     *  return[gausspoint](1,0) = dy/dXi
     *  return[gausspoint](1,1) = dy/dEta

     * @note d_Xi_shape_f_values and d_Eta_shape_f_values also contain information
     * at each gauss point.
     * Thereby the storage scheme is for both of them like F_i_j where
     * i = Control Point index in Xi-direction
     * j = Control Point index in Eta-direction
     */

    ShapeFunctionsGradientsType
    CalculateShapeFunctionsIntegrationPointsLocalGradients(
            typename BaseType::IntegrationMethod ThisMethod,
            ShapeFunctionsGradientsType &d_Xi_shape_f_values,
            ShapeFunctionsGradientsType &d_Eta_shape_f_values) const
    {

        Matrix ShapeFunctionsValues;
        IntegrationPointsContainerType all_integration_points =AllIntegrationPoints();
        IntegrationPointsArrayType integration_points = all_integration_points[ThisMethod];
        const int integration_points_number = integration_points.size();
        ShapeFunctionsGradientsType d_shape_f_values( integration_points_number );
        d_shape_f_values = CalculateShapeFunctionsIntegrationPointsLocalGradients(
                ThisMethod,
                d_Xi_shape_f_values,
                d_Eta_shape_f_values,
                ShapeFunctionsValues);

        return d_shape_f_values;
    }





    /**
     * CalculateShapeFunctionsIntegrationPointsLocalGradients is a function to compute
     * the local derivatives dx/dXi, dx/dEta, dy/dXi and dy/dEta for each integration
     * point.
     *
     * @param ThisMethod is the integration method used
     *
     * @return a Vector of Matrices. Each Matrix represents the local derivatives at one gauss point
     *  return[gausspoint](0,0) = dx/dXi
     *  return[gausspoint](0,1) = dx/dEta
     *  return[gausspoint](1,0) = dy/dXi
     *  return[gausspoint](1,1) = dy/dEta
     */

    ShapeFunctionsGradientsType
    CalculateShapeFunctionsIntegrationPointsLocalGradients(
        typename BaseType::IntegrationMethod ThisMethod ) const
    {

        Matrix ShapeFunctionsValues;

        IntegrationPointsContainerType all_integration_points =AllIntegrationPoints();
        IntegrationPointsArrayType integration_points = all_integration_points[ThisMethod];
        const int integration_points_number = integration_points.size();
        ShapeFunctionsGradientsType d_shape_f_values( integration_points_number );
        ShapeFunctionsGradientsType d_Xi_shape_f_values( integration_points_number );
        ShapeFunctionsGradientsType d_Eta_shape_f_values( integration_points_number );
        d_shape_f_values = CalculateShapeFunctionsIntegrationPointsLocalGradients(
                ThisMethod,
                d_Xi_shape_f_values,
                d_Eta_shape_f_values,
                ShapeFunctionsValues);

        return d_shape_f_values;
    }


    /**
     * TODO: implemented but not yet tested
     */
    /**
     * Calculates the values of all shape function in all integration points.
     * Integration points are expected to be given in local coordinates
     * @param ThisMethod the current integration method
     * @return the matrix of values of every shape function in each integration point
     */
    Matrix CalculateShapeFunctionsIntegrationPointsValues(
        typename BaseType::IntegrationMethod ThisMethod )
    {
        IntegrationPointsContainerType all_integration_points =
            AllIntegrationPoints();
        IntegrationPointsArrayType integration_points =
            all_integration_points[ThisMethod];
        //number of integration points
        const int integration_points_number = integration_points.size();
        //number of nodes in current geometry
        const int points_number = this->Points().size();
        //setting up return matrix
        Matrix shape_function_values( integration_points_number, points_number );
        //loop over all integration points

    Matrix NurbsBasisFunction(mPolynomialDegreeP+1,mPolynomialDegreeQ+1);
    double XiLocalCoordinates(0), EtaLocalCoordinates(0);

        for ( int pnt = 0; pnt < integration_points_number; pnt++ )
        {
            XiLocalCoordinates = (mUpperXi - mLowerXi) / 2 * (1+integration_points[pnt].X()) + mLowerXi;
            EtaLocalCoordinates = (mUpperEta - mLowerEta) / 2 * (1+integration_points[pnt].Y()) + mLowerEta;
            NurbsBasisFunction = ElementGeometryNurbsFunctionsValues(XiLocalCoordinates,EtaLocalCoordinates);

            for (int i = 0; i < mPolynomialDegreeP+1; i++)
            {
               for (int j=0; j < mPolynomialDegreeQ+1; j++)
                {
                  shape_function_values( pnt, j*(mPolynomialDegreeQ+1) + i ) = NurbsBasisFunction(i,j);
                }
            }

        }

        return shape_function_values;
    }





    /**
     * TODO: testing
     */
    static const IntegrationPointsContainerType AllIntegrationPoints()
    {
        IntegrationPointsContainerType integration_points =
        {
            {
                Quadrature < QuadrilateralGaussLegendreIntegrationPoints1,
                2, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < QuadrilateralGaussLegendreIntegrationPoints2,
                2, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < QuadrilateralGaussLegendreIntegrationPoints3,
                2, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < QuadrilateralGaussLegendreIntegrationPoints4,
                2, IntegrationPoint<3> >::GenerateIntegrationPoints()
            }
        };

        return integration_points;
    }

    /**
     * TODO: testing
     */


    static const ShapeFunctionsValuesContainerType AllShapeFunctionsValues()
    {
        ShapeFunctionsValuesContainerType shape_functions_values =
        {
            {

            }
        };
        return shape_functions_values;
    }

    /**
     * TODO: testing
     */


    static const ShapeFunctionsLocalGradientsContainerType
    AllShapeFunctionsLocalGradients()
    {
        ShapeFunctionsLocalGradientsContainerType shape_functions_local_gradients =
        {
            {
            }
        };
        return shape_functions_local_gradients;
    }



    /**
     * ON_ShapeFunctionValue is one of the two used OpenNURBS functions.
     * It is used to compute all potentially non-zero B-Spline(!) shape
     * functions values at a point (Xi or Eta) in the Xi- /Eta - parameter
     * domain respectively.
     *
     * @param order is the polynomial degree +1
     * @param knot is the already reduced(! -> @note) knot vector but has to be provided as c-array
     * @param t the Xi/Eta-Value where the shape functions are evaluated
     * @param N container of the shape functions values
     *
     * @return true if the computation run without problems
     *
     * @note The returned Shape Function Values will be contained in
     * N[0][j] with j=0 ... order-1 = 0 ... polynomial degree
     * N[1][j] will give the evaluated B-Spline progenitor (polynomial degree -1)
     * N[i][j] will give the j-th shape function value with polynomial degree = "highest" polynomial degree - i
     *
     * Reduced knot vector means, that only the knots needed for computing the
     * shape functions are passed (length = 2*polynomial degree). These reduced knot
     * vector can be obtained via the function FindKnotsToEvaluate(...)
     *
     */
    bool ON_ShapeFunctionValue( int order,
                                const double *knot,
                                double t,
                                double *N )const
    {
    /*****************************************************************************
    Evaluate B-spline basis functions

    INPUT:
      order >= 1
        d = degree = order - 1
      knot[]
        array of length 2*d.
        Generally, knot[0] <= ... <= knot[d-1] < knot[d] <= ... <= knot[2*d-1].
      N[]
        array of length order*order

    OUTPUT:
      If "N" were declared as double N[order][order], then

                     k
        N[d-k][i] = N (t) = value of i-th degree k basis function.
                     i
      where 0 <= k <= d and k <= i <= d.

        In particular, N[0], ..., N[d] - values of degree d basis functions.
      The "lower left" triangle is not initialized.

      Actually, the above is true when knot[d-1] <= t < knot[d].  Otherwise, the
      value returned is the value of the polynomial that agrees with N_i^k on the
      half open domain [ knot[d-1], knot[d] )

    COMMENTS:
      If a degree d NURBS has n control points, then the TL knot vector has
      length d+n-1. ( Most literature, including DeBoor and The NURBS Book,
      duplicate the TL start and end knot and have knot vectors of length
      d+n+1. )

      Assume C is a B-spline of degree d (order=d+1) with n control vertices
      (n>=d+1) and knot[] is its knot vector.  Then

        C(t) = Sum( 0 <= i < n, N_{i}(t) * C_{i} )

      where N_{i} are the degree d b-spline basis functions and C_{i} are the control
      vertices.  The knot[] array length d+n-1 and satisfies

        knot[0] <= ... <= knot[d-1] < knot[d]
        knot[n-2] < knot[n-1] <= ... <= knot[n+d-2]
        knot[i] < knot[d+i] for 0 <= i < n-1
        knot[i] <= knot[i+1] for 0 <= i < n+d-2

      The domain of C is [ knot[d-1], knot[n-1] ].

      The support of N_{i} is [ knot[i-1], knot[i+d] ).

      If d-1 <= k < n-1 and knot[k] <= t < knot[k+1], then
      N_{i}(t) = 0 if i <= k-d
               = 0 if i >= k+2
               = B[i-k+d-1] if k-d+1 <= i <= k+1, where B[] is computed by the call
                 TL_EvNurbBasis( d+1, knot+k-d+1, t, B );

      If 0 <= j < n-d, 0 <= m <= d, knot[j+d-1] <= t < knot[j+d], and B[] is
      computed by the call

        TL_EvNurbBasis( d+1, knot+j, t, B ),

      then

        N_{j+m}(t) = B[m].

    EXAMPLE:
    REFERENCE:
      The NURBS book
    RELATED FUNCTIONS:
      TL_EvNurbBasis
      TL_EvNurbBasisDer
    *****************************************************************************/
      register double a0, a1, x, y;
      const double *k0;
      double *t_k, *k_t, *N0;
      const int d = order-1;
      register int j, r;

      t_k = (double*)alloca( d<<4 );
          k_t = t_k + d;

      if (knot[d-1] == knot[d]) {
            /* value is defined to be zero on empty spans */
        memset( N, 0, order*order*sizeof(*N) );
        return true;
      }

      N  += order*order-1;
        N[0] = 1.0;
      knot += d;
      k0 = knot - 1;

      for (j = 0; j < d; j++ ) {
            N0 = N;
        N -= order+1;
        t_k[j] = t - *k0--;
        k_t[j] = *knot++ - t;

        x = 0.0;
        for (r = 0; r <= j; r++) {
          a0 = t_k[j-r];
          a1 = k_t[r];
          y = N0[r]/(a0 + a1);
          N[r] = x + a1*y;
          x = a0*y;
        }

        N[r] = x;
      }

      //   16 September 2003 Dale Lear (at Chuck's request)
      //   When t is at an end knot, do a check to
      //   get exact values of basis functions.
      //   The problem being that a0*y above can
      //   fail to be one by a bit or two when knot
      //   values are large.
      x = 1.0-ON_SQRT_EPSILON;
      if ( N[0] > x )
      {
        if ( N[0] != 1.0 && N[0] < 1.0 + ON_SQRT_EPSILON )
        {
          r = 1;
          for ( j = 1; j <= d && r; j++ )
          {
            if ( N[j] != 0.0 )
              r = 0;
          }
          if (r)
            N[0] = 1.0;
        }
      }
      else if ( N[d] > x )
      {
        if ( N[d] != 1.0 && N[d] < 1.0 + ON_SQRT_EPSILON )
        {
          r = 1;
          for ( j = 0; j < d && r; j++ )
          {
            if ( N[j] != 0.0 )
              r = 0;
          }
          if (r)
            N[d] = 1.0;
        }
      }

      return true; }




    /**
     * ON_EvaluateNurbsBasisDerivatives is one of the two used OpenNURBS functions.
     * It is used to compute all potentially non-zero B-Spline(!) shape
     * functions derivatives at a point (Xi or Eta) in the Xi- /Eta - parameter
     * domain respectively
     *
     * @param order is the polynomial degree +1
     * @param knot is the already reduced(!) knot vector but has to be provided as c-array
     * @param der_count is the derivative depth until which the B-Splines(!) will be derivated
     * @param N has to be the output of the function ON_ShapeFunctionValue()!!! An will contain
     *        the derivatives after the execution of the function.
     *
     * @note N[order][order] is the maximum needed memory for this function (if der_count = polynomial degree)
     * where
     * N[0][j] for j=0...order-1 = 0...polynomial degree will contain the j-th shape function value
     * N[1][j] contains the 1st derivative of the shape functions
     * N[i][j] contains the i-th derivative of the j-th shape function
     *
     * All derivatives ara evaluated at the same point in the parameter space as the input parameter N
     * was evaluated in the function ON_ShapeFunctionValue()
     *
     * Reduced knot vector means, that only the knots needed for computing the
     * shape functions are passed (length = 2*polynomial degree). These reduced knot
     * vectors can be obtained via the function FindKnotsToEvaluate(...)
     */
    bool ON_EvaluateNurbsBasisDerivatives( int order, const double* knot,
                           int der_count, double* N )const


    {



        /* INPUT:
         *   Results of the call
         *      TL_EvNurbBasis( order, knot, t, N );  (initializes N[] )
         *   are sent to
         *      TL_EvNurbBasisDer( order, knot, der_count, N ),
         *   where 1 <= der_count < order
         *
         * OUTPUT:
       *  If "N" were declared as double N[order][order], then
         *
       *                                    d
       *    N[d-k][i] = k-th derivative of N (t)
       *                                    i
       *
         *  where 0 <= k <= d and 0 <= i <= d.
         *
         * In particular,
         *   N[0], ..., N[d] - values of degree d basis functions.
         *   N[order], ..., N[order_d] - values of first derivative.
         *
       * Actually, the above is true when knot[d-1] <= t < knot[d].  Otherwise, the
       * values returned are the values of the polynomials that agree with N_i^k on the
       * half open domain [ knot[d-1], knot[d] )
         *
         * Ref: The NURBS Book
         */
        double dN, c;
        const double *k0, *k1;
        //const Vector &k0, &k1;
        double *a0, *a1, *ptr, **dk;
        //Vector &a0, &a1, &ptr;
        //Matrix dk;
        int i, j, k, jmax;

        const int d = order-1;
        const int Nstride = -der_count*order;

        /* workspaces for knot differences and coefficients
         *
         * a0[] and a1[] have order doubles
         *
         * dk[0] = array of d knot differences
         * dk[1] = array of (d-1) knot differences
         *
         * dk[der_count-1] = 1.0/(knot[d] - knot[d-1])
         * dk[der_count] = dummy pointer to make loop efficient
         */
        dk = (double**)alloca( (der_count+1) << 3 ); /* << 3 in case pointers are 8 bytes long */
        a0 = (double*)alloca( (order*(2 + ((d+1)>>1))) << 3 ); /* d for a0, d for a1, d*order/2 for dk[]'s and slop to avoid /2 */
        a1 = a0 + order;

        /* initialize reciprocal of knot differences */
        dk[0] = a1 + order;
        for (k = 0; k < der_count; k++) {
            j = d-k;
            k0 = knot++;
            k1 = k0 + j;
            for (i = 0; i < j; i++)
                dk[k][i] = 1.0/(*k1++ - *k0++);
            dk[k+1] = dk[k] + j;
        }
        dk--;
        /* dk[1] = 1/{t[d]-t[0], t[d+1]-t[1], ..., t[2d-2] - t[d-2], t[2d-1] - t[d-1]}
         *       = diffs needed for 1rst derivative
         * dk[2] = 1/{t[d]-t[1], t[d+1]-t[2], ..., t[2d-2] - t[d-1]}
         *       = diffs needed for 2nd derivative
         * ...
         * dk[d] = 1/{t[d] - t[d-1]}
         *       = diff needed for d-th derivative
         *
         * d[k][n] = 1.0/( t[d+n] - t[k-1+n] )
         */

        N += order;
        /* set N[0] ,..., N[d] = 1rst derivatives,
         * N[order], ..., N[order+d] = 2nd, etc.
         */
        for ( i=0; i<order; i++) {
            a0[0] = 1.0;
            for (k = 1; k <= der_count; k++) {
                /* compute k-th derivative of N_i^d up to d!/(d-k)! scaling factor */
                dN = 0.0;
                j = k-i;
                if (j <= 0) {
                    dN = (a1[0] = a0[0]*dk[k][i-k])*N[i];
                    j = 1;
                }
                jmax = d-i;
                if (jmax < k) {
                    while (j <= jmax) {
                        dN += (a1[j] = (a0[j] - a0[j-1])*dk[k][i+j-k])*N[i+j];
                        j++;
                    }
                }
                else {
                    /* sum j all the way to j = k */
                    while (j < k) {
                        dN += (a1[j] = (a0[j] - a0[j-1])*dk[k][i+j-k])*N[i+j];
                        j++;
                    }
                    dN += (a1[k] = -a0[k-1]*dk[k][i])*N[i+k];
                }

                /* d!/(d-k)!*dN = value of k-th derivative */
                N[i] = dN;
                N += order;
                /* a1[] s for next derivative = linear combination
                 * of a[]s used to compute this derivative.
                 */
                ptr = a0; a0 = a1; a1 = ptr;
            }
            N += Nstride;
        }

        /* apply d!/(d-k)! scaling factor */
        dN = c = (double)d;
        k = der_count;
        while (k--) {
            i = order;
            while (i--)
                *N++ *= c;
            dN -= 1.0;
            c *= dN;
        }
      return true;
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

    template<class TOtherPointType> friend class NurbsPatchGeometry2D;

    ///@}
    ///@name Un accessible methods
    ///@{



    ///@}
}; // Class Geometry



template<class TPointType> const
GeometryData NurbsPatchGeometry2D<TPointType>::msGeometryData(
    2, 2, 2,
    GeometryData::GI_GAUSS_1,
    NurbsPatchGeometry2D<TPointType>::AllIntegrationPoints(),
    NurbsPatchGeometry2D<TPointType>::AllShapeFunctionsValues(),
    AllShapeFunctionsLocalGradients()
);

}// namespace Kratos.

#endif
