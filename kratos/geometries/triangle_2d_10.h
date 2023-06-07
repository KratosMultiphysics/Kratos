//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Mohamed Nabi
//
//
//  Contributors:
//
//

#pragma once

// System includes

// External includes

// Project includes
#include "geometries/line_2d_4.h"
#include "integration/triangle_gauss_legendre_integration_points.h"

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
 * @class Triangle2D10
 * @ingroup KratosCore
 * @brief A ten node 2D triangular geometry with cubic shape functions
 * @details While the shape functions are only defined in 2D it is possible to define an arbitrary orientation in space. Thus it can be used for defining surfaces on 3D elements.
 * The node ordering corresponds with:
 *          2
 *          |`\
 *          7   6
 *          |    `\
 *          8   9  `5
 *          |        `\
 *          0---3--4---1
 * @author Mohamed Nabi
 * @author
 * @author
 */
template<class TPointType>
class Triangle2D10 : public Geometry<TPointType>
{
public:
    ///@}
    ///@name Type Definitions
    ///@{

    /**
     * Geometry as base class.
     */
    typedef Geometry<TPointType> BaseType;

    /**
     * Type of edge geometry
     */
    typedef Line2D4<TPointType> EdgeType;

    /**
     * Pointer definition of Triangle2D10
     */
    KRATOS_CLASS_POINTER_DEFINITION(Triangle2D10);

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
     * Type of the normal vector used for normal to edges in geometry.
     */
    typedef typename BaseType::NormalType NormalType;

    ///@}
    ///@name Life Cycle
    ///@{

    Triangle2D10(const PointType& Point01, const PointType& Point02, const PointType& Point03,
                 const PointType& Point04, const PointType& Point05, const PointType& Point06,
                 const PointType& Point07, const PointType& Point08, const PointType& Point09,
                 const PointType& Point10) : BaseType(PointsArrayType(), &msGeometryData)
    {
        auto& r_points = this->Points();
        r_points.push_back(typename PointType::Pointer(new PointType(Point01)));
        r_points.push_back(typename PointType::Pointer(new PointType(Point02)));
        r_points.push_back(typename PointType::Pointer(new PointType(Point03)));
        r_points.push_back(typename PointType::Pointer(new PointType(Point04)));
        r_points.push_back(typename PointType::Pointer(new PointType(Point05)));
        r_points.push_back(typename PointType::Pointer(new PointType(Point06)));
        r_points.push_back(typename PointType::Pointer(new PointType(Point07)));
        r_points.push_back(typename PointType::Pointer(new PointType(Point08)));
        r_points.push_back(typename PointType::Pointer(new PointType(Point09)));
        r_points.push_back(typename PointType::Pointer(new PointType(Point10)));
    }

    Triangle2D10(typename PointType::Pointer pPoint01, typename PointType::Pointer pPoint02,
                 typename PointType::Pointer pPoint03, typename PointType::Pointer pPoint04,
                 typename PointType::Pointer pPoint05, typename PointType::Pointer pPoint06,
                 typename PointType::Pointer pPoint07, typename PointType::Pointer pPoint08,
                 typename PointType::Pointer pPoint09, typename PointType::Pointer pPoint10)
                 : BaseType(PointsArrayType(), &msGeometryData)
    {
        auto& r_points = this->Points();
        r_points.push_back(pPoint01);
        r_points.push_back(pPoint02);
        r_points.push_back(pPoint03);
        r_points.push_back(pPoint04);
        r_points.push_back(pPoint05);
        r_points.push_back(pPoint06);
        r_points.push_back(pPoint07);
        r_points.push_back(pPoint08);
        r_points.push_back(pPoint09);
        r_points.push_back(pPoint10);
    }

    Triangle2D10(const PointsArrayType& ThisPoints) : BaseType(ThisPoints, &msGeometryData)
    {
        if (this->PointsNumber() != 10)
        {
            KRATOS_ERROR << "Invalid points number. Expected 10, given " << this->PointsNumber()
                << std::endl;
        }
    }

    /// Constructor with Geometry Id
    explicit Triangle2D10(const IndexType GeometryId, const PointsArrayType& rThisPoints)
        : BaseType(GeometryId, rThisPoints, &msGeometryData)
    {
        KRATOS_ERROR_IF(this->PointsNumber() != 10) << "Invalid points number. Expected 10, given "
                                                    << this->PointsNumber() << std::endl;
    }

    /// Constructor with Geometry Name
    explicit Triangle2D10(const std::string& rGeometryName, const PointsArrayType& rThisPoints)
        : BaseType(rGeometryName, rThisPoints, &msGeometryData)
    {
        KRATOS_ERROR_IF(this->PointsNumber() != 10) << "Invalid points number. Expected 10, given "
                                                    << this->PointsNumber() << std::endl;
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
    Triangle2D10(Triangle2D10 const& rOther) : BaseType(rOther)
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
    template<class TOtherPointType> Triangle2D10(Triangle2D10<TOtherPointType> const& rOther)
        : BaseType(rOther)
    {
    }

    /**
     * Destructor. Does nothing!!!
     */
    ~Triangle2D10() override {}

    GeometryData::KratosGeometryFamily GetGeometryFamily() const override
    {
        return GeometryData::KratosGeometryFamily::Kratos_Triangle;
    }

    GeometryData::KratosGeometryType GetGeometryType() const override
    {
        return GeometryData::KratosGeometryType::Kratos_Triangle2D10;
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
    Triangle2D10& operator=(const Triangle2D10& rOther)
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
    Triangle2D10& operator=(Triangle2D10<TOtherPointType> const& rOther)
    {
        BaseType::operator=(rOther);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Creates a new geometry pointer
     * @param rThisPoints the nodes of the new geometry
     * @return Pointer to the new geometry
     */
    typename BaseType::Pointer Create(PointsArrayType const& rThisPoints) const override
    {
        return typename BaseType::Pointer(new Triangle2D10(rThisPoints));
    }

    /**
     * @brief Creates a new geometry pointer
     * @param NewGeometryId the ID of the new geometry
     * @param rThisPoints the nodes of the new geometry
     * @return Pointer to the new geometry
     */
    typename BaseType::Pointer Create(const IndexType NewGeometryId,
        PointsArrayType const& rThisPoints) const override
    {
        return typename BaseType::Pointer(new Triangle2D10(NewGeometryId, rThisPoints));
    }

    /**
     * @brief Creates a new geometry pointer
     * @param rGeometry reference to an existing geometry
     * @return Pointer to the new geometry
     */
    typename BaseType::Pointer Create(const BaseType& rGeometry) const override
    {
        auto p_geometry = typename BaseType::Pointer(new Triangle2D10(rGeometry.Points()));
        p_geometry->SetData(rGeometry.GetData());
        return p_geometry;
    }

    /**
     * @brief Creates a new geometry pointer
     * @param NewGeometryId the ID of the new geometry
     * @param rGeometry reference to an existing geometry
     * @return Pointer to the new geometry
     */
    typename BaseType::Pointer Create(const IndexType NewGeometryId, const BaseType& rGeometry) const override
    {
        auto p_geometry = typename BaseType::Pointer(new Triangle2D10(NewGeometryId, rGeometry.Points()));
        p_geometry->SetData(rGeometry.GetData());
        return p_geometry;
    }

    /**
     * returns the local coordinates of all nodes of the current geometry
     * @param rResult a Matrix object that will be overwritten by the result
     * @return the local coordinates of all nodes
     */
    Matrix& PointsLocalCoordinates(Matrix& rResult) const override
    {
        rResult.resize(10, 2, false);
        //
        const double oneThird = 1.0 / 3.0;
        const double twoThird = 2.0 / 3.0;
        //
        rResult(0, 0) = 0.0;
        rResult(0, 1) = 0.0;
        //
        rResult(1, 0) = 1.0;
        rResult(1, 1) = 0.0;
        //
        rResult(2, 0) = 0.0;
        rResult(2, 1) = 1.0;
        //
        rResult(3, 0) = oneThird;
        rResult(3, 1) = 0.0;
        //
        rResult(4, 0) = twoThird;
        rResult(4, 1) = 0.0;
        //
        rResult(5, 0) = twoThird;
        rResult(5, 1) = oneThird;
        //
        rResult(6, 0) = oneThird;
        rResult(6, 1) = twoThird;
        //
        rResult(7, 0) = 0.0;
        rResult(7, 1) = twoThird;
        //
        rResult(8, 0) = 0.0;
        rResult(8, 1) = oneThird;
        //
        rResult(9, 0) = oneThird;
        rResult(9, 1) = oneThird;
        //
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
    double Length() const override
    {
        return std::sqrt(std::abs(Area()));
    }

    /** This method calculates and returns area or surface area of
     * this geometry depending to it's dimension. For one dimensional
     * geometry it returns zero, for two dimensional it gives area
     * and for three dimensional geometries it gives surface area.
     *
     * @return double value contains area or surfacede
     * area.
     * @see Length()
     * @see Volume()
     * @see DomainSize()
     */
    double Area() const override
    {
        Vector temp;
        this->DeterminantOfJacobian(temp, msGeometryData.DefaultIntegrationMethod());
        const IntegrationPointsArrayType& integration_points = this->IntegrationPoints(msGeometryData.DefaultIntegrationMethod());
        double area = 0.00;
        for (unsigned int i = 0; i < integration_points.size(); ++i)
        {
            area += temp[i] * integration_points[i].Weight();
        }
        return area;
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
    double DomainSize() const override
    {
        return Area();
    }

    /**
     * @brief Returns whether given arbitrary point is inside the Geometry and the respective
     * local point for the given global point
     * @param rPoint The point to be checked if is inside o note in global coordinates
     * @param rResult The local coordinates of the point
     * @param Tolerance The  tolerance that will be considered to check if the point is inside or not
     * @return True if the point is inside, false otherwise
     */
    bool IsInside(const CoordinatesArrayType& rPoint, CoordinatesArrayType& rResult,
        const double Tolerance = std::numeric_limits<double>::epsilon()) const override
    {
        this->PointLocalCoordinates(rResult, rPoint);
        if ((rResult[0] >= (0.0 - Tolerance)) && (rResult[0] <= (1.0 + Tolerance))) {
            if ((rResult[1] >= (0.0 - Tolerance)) && (rResult[1] <= (1.0 + Tolerance))) {
                if ((rResult[0] + rResult[1]) <= (1.0 + Tolerance)) {
                    return true;
                }
            }
        }
        return false;
    }

    ///@}
    ///@name Shape Function
    ///@{

    /**
     * @brief Returns vector of shape function values at local coordinate.
     *
     * For a definition of the shape functions see, e.g.,
     * https://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/IFEM.Ch18.d/IFEM.Ch18.pdf.
     */
    Vector& ShapeFunctionsValues(Vector& rResult, const CoordinatesArrayType& rCoordinates) const override
    {
        if (rResult.size() != 10) rResult.resize(10, false);
        const double xi = rCoordinates[0];
        const double et = rCoordinates[1];
        const double zt = 1.0 - xi - et;
        //
        rResult[0] = zt * (3.0 * zt - 1.0) * (3.0 * zt - 2.0) * 0.5;
        rResult[1] = xi * (3.0 * xi - 1.0) * (3.0 * xi - 2.0) * 0.5;
        rResult[2] = et * (3.0 * et - 1.0) * (3.0 * et - 2.0) * 0.5;
        rResult[3] = xi * zt * (3.0 * zt - 1.0) * 4.5;
        rResult[4] = xi * zt * (3.0 * xi - 1.0) * 4.5;
        rResult[5] = xi * et * (3.0 * xi - 1.0) * 4.5;
        rResult[6] = xi * et * (3.0 * et - 1.0) * 4.5;
        rResult[7] = et * zt * (3.0 * et - 1.0) * 4.5;
        rResult[8] = et * zt * (3.0 * zt - 1.0) * 4.5;
        rResult[9] = xi * et * zt * 27.0;
        //
        return rResult;
    }

    /**
     * Calculates the value of a given shape function at a given point.
     *
     * @param ShapeFunctionIndex The number of the desired shape function
     * @param rPoint the given point in local coordinates at which the value of the shape
     * function is calculated
     *
     * @return the value of the shape function at the given point
     */
    double ShapeFunctionValue(IndexType ShapeFunctionIndex, const CoordinatesArrayType& rPoint) const override
    {
        const double xi = rPoint[0];
        const double et = rPoint[1];
        const double zt = 1.0 - xi - et;
        double shape = 0.0;
        //
        switch (ShapeFunctionIndex)
        {
        case 0:
            shape = zt * (3.0 * zt - 1.0) * (3.0 * zt - 2.0) * 0.5;
            break;
        case 1:
            shape = xi * (3.0 * xi - 1.0) * (3.0 * xi - 2.0) * 0.5;
            break;
        case 2:
            shape = et * (3.0 * et - 1.0) * (3.0 * et - 2.0) * 0.5;
            break;
        case 3:
            shape = xi * zt * (3.0 * zt - 1.0) * 4.5;
            break;
        case 4:
            shape = xi * zt * (3.0 * xi - 1.0) * 4.5;
            break;
        case 5:
            shape = xi * et * (3.0 * xi - 1.0) * 4.5;
            break;
        case 6:
            shape = xi * et * (3.0 * et - 1.0) * 4.5;
            break;
        case 7:
            shape = et * zt * (3.0 * et - 1.0) * 4.5;
            break;
        case 8:
            shape = et * zt * (3.0 * zt - 1.0) * 4.5;
            break;
        case 9:
            shape = xi * et * zt * 27.0;
            break;
        default:
            KRATOS_ERROR << "Wrong index of shape function!" << *this << std::endl;
            break;
        }
        //
        return shape;
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
    std::string Info() const override
    {
        return "2 dimensional triangle with ten nodes in 2D space";
    }

    /**
     * Print information about this object.
     * @param rOStream Stream to print into it.
     * @see PrintData()
     * @see Info()
     */
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "2 dimensional triangle with ten nodes in 2D space";
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
    void PrintData(std::ostream& rOStream) const override
    {
        PrintInfo(rOStream);
        BaseType::PrintData(rOStream);
        std::cout << std::endl;
        Matrix jacobian;
        this->Jacobian(jacobian, PointType());
        rOStream << "    Jacobian in the origin\t : " << jacobian;
    }

    ///@}
    ///@name Edge
    ///@{

    /**
     * @brief This method gives you number of all edges of this geometry.
     * @details For example, for a hexahedron, this would be 12
     * @return SizeType containes number of this geometry edges.
     * @see EdgesNumber()
     * @see Edges()
     * @see GenerateEdges()
     * @see FacesNumber()
     * @see Faces()
     * @see GenerateFaces()
     */
    SizeType EdgesNumber() const override
    {
        return 3;
    }

    /**
     * @brief This method gives you all edges of this geometry.
     * @details This method will gives you all the edges with one dimension less than this geometry.
     * For example a triangle would return three lines as its edges or a tetrahedral would return four triangle as its edges but won't return its six edge lines by this method.
     * @return GeometriesArrayType containes this geometry edges.
     * @see EdgesNumber()
     * @see Edge()
     */
    GeometriesArrayType GenerateEdges() const override
    {
        GeometriesArrayType edges = GeometriesArrayType();
        edges.push_back(Kratos::make_shared<EdgeType>(this->pGetPoint(0), this->pGetPoint(1), this->pGetPoint(3), this->pGetPoint(4)));
        edges.push_back(Kratos::make_shared<EdgeType>(this->pGetPoint(1), this->pGetPoint(2), this->pGetPoint(5), this->pGetPoint(6)));
        edges.push_back(Kratos::make_shared<EdgeType>(this->pGetPoint(2), this->pGetPoint(0), this->pGetPoint(7), this->pGetPoint(8)));
        return edges;
    }

    SizeType FacesNumber() const override
    {
        return 3;
    }

    //Connectivities of faces required
    void NumberNodesInFaces(DenseVector<unsigned int>& NumberNodesInFaces) const override
    {
        if (NumberNodesInFaces.size() != 3) NumberNodesInFaces.resize(3, false);
        NumberNodesInFaces[0] = 4;
        NumberNodesInFaces[1] = 4;
        NumberNodesInFaces[2] = 4;
    }

    void NodesInFaces(DenseMatrix<unsigned int>& NodesInFaces) const override
    {
        // faces in columns
        if (NodesInFaces.size1() != 5 || NodesInFaces.size2() != 3)
            NodesInFaces.resize(5, 3, false);
        //
        //compatible to Edges() function. Triangle library considers a different order
        NodesInFaces(0, 0) = 0;//face or master node
        NodesInFaces(1, 0) = 1;
        NodesInFaces(2, 0) = 5;
        NodesInFaces(3, 0) = 6;
        NodesInFaces(4, 0) = 2;
        //
        NodesInFaces(0, 1) = 1;//face or master node
        NodesInFaces(1, 1) = 2;
        NodesInFaces(2, 1) = 7;
        NodesInFaces(3, 1) = 8;
        NodesInFaces(4, 1) = 0;
        //
        NodesInFaces(0, 2) = 2;//face or master node
        NodesInFaces(1, 2) = 0;
        NodesInFaces(2, 2) = 3;
        NodesInFaces(3, 2) = 4;
        NodesInFaces(4, 2) = 1;
    }

    /**
     * Calculates the local gradients for all integration points for
     * given integration method
     */
    virtual ShapeFunctionsGradientsType ShapeFunctionsLocalGradients(IntegrationMethod ThisMethod)
    {
        ShapeFunctionsGradientsType localGradients
            = CalculateShapeFunctionsIntegrationPointsLocalGradients(ThisMethod);
        const int integration_points_number = msGeometryData.IntegrationPointsNumber(ThisMethod);
        ShapeFunctionsGradientsType Result(integration_points_number);
        for (int pnt = 0; pnt < integration_points_number; ++pnt)
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
            = CalculateShapeFunctionsIntegrationPointsLocalGradients(ThisMethod);
        const int integration_points_number = msGeometryData.IntegrationPointsNumber(ThisMethod);
        ShapeFunctionsGradientsType Result(integration_points_number);
        for (int pnt = 0; pnt < integration_points_number; ++pnt)
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
    Matrix& ShapeFunctionsLocalGradients(Matrix& rResult, const CoordinatesArrayType& rPoint) const override
    {
        rResult.resize(10, 2, false);
        const double xi = rPoint[0];
        const double et = rPoint[1];
        const double zt = 1.0 - xi - et;
        //
        noalias(rResult) = ZeroMatrix(10, 2);
        //
        rResult(0, 0) = -4.5 * zt * (3.0 * zt - 2.0) - 1.0;
        rResult(0, 1) = -4.5 * zt * (3.0 * zt - 2.0) - 1.0;
        //
        rResult(1, 0) = 4.5 * xi * (3.0 * xi - 2.0) + 1.0;
        rResult(1, 1) = 0.0;
        //
        rResult(2, 0) = 0.0;
        rResult(2, 1) = 4.5 * et * (3.0 * et - 2.0) + 1.0;
        //
        rResult(3, 0) = 4.5 * (zt * (3.0 * zt - 1.0) - xi * (6.0 * zt - 1.0));
        rResult(3, 1) = -4.5 * xi * (6.0 * zt - 1.0);
        //
        rResult(4, 0) = 4.5 * (zt * (6.0 * xi - 1.0) - xi * (3.0 * xi - 1.0));
        rResult(4, 1) = -4.5 * xi * (3.0 * xi - 1.0);
        //
        rResult(5, 0) = 4.5 * et * (6.0 * xi - 1.0);
        rResult(5, 1) = 4.5 * xi * (3.0 * xi - 1.0);
        //
        rResult(6, 0) = 4.5 * et * (3.0 * et - 1.0);
        rResult(6, 1) = 4.5 * xi * (6.0 * et - 1.0);
        //
        rResult(7, 0) = -4.5 * et * (3.0 * et - 1.0);
        rResult(7, 1) = 4.5 * (zt * (6.0 * et - 1.0) - et * (3.0 * et - 1.0));
        //
        rResult(8, 0) = -4.5 * et * (6.0 * zt - 1.0);
        rResult(8, 1) = 4.5 * (zt * (3.0 * zt - 1.0) - et * (6.0 * zt - 1.0));
        //
        rResult(9, 0) = 27.0 * et * (zt - xi);
        rResult(9, 1) = 27.0 * xi * (zt - et);
        //
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
    virtual Matrix& ShapeFunctionsGradients(Matrix& rResult, const CoordinatesArrayType& rPoint)
    {
        rResult.resize(10, 2, false);
        //
        const double xi = rPoint[0];
        const double et = rPoint[1];
        const double zt = 1.0 - xi - et;
        //
        noalias(rResult) = ZeroMatrix(10, 2);
        //
        rResult(0, 0) = -4.5 * zt * (3.0 * zt - 2.0) - 1.0;
        rResult(0, 1) = -4.5 * zt * (3.0 * zt - 2.0) - 1.0;
        //
        rResult(1, 0) = 4.5 * xi * (3.0 * xi - 2.0) + 1.0;
        rResult(1, 1) = 0.0;
        //
        rResult(2, 0) = 0.0;
        rResult(2, 1) = 4.5 * et * (3.0 * et - 2.0) + 1.0;
        //
        rResult(3, 0) = 4.5 * (zt * (3.0 * zt - 1.0) - xi * (6.0 * zt - 1.0));
        rResult(3, 1) = -4.5 * xi * (6.0 * zt - 1.0);
        //
        rResult(4, 0) = 4.5 * (zt * (6.0 * xi - 1.0) - xi * (3.0 * xi - 1.0));
        rResult(4, 1) = -4.5 * xi * (3.0 * xi - 1.0);
        //
        rResult(5, 0) = 4.5 * et * (6.0 * xi - 1.0);
        rResult(5, 1) = 4.5 * xi * (3.0 * xi - 1.0);
        //
        rResult(6, 0) = 4.5 * et * (3.0 * et - 1.0);
        rResult(6, 1) = 4.5 * xi * (6.0 * et - 1.0);
        //
        rResult(7, 0) = -4.5 * et * (3.0 * et - 1.0);
        rResult(7, 1) = 4.5 * (zt * (6.0 * et - 1.0) - et * (3.0 * et - 1.0));
        //
        rResult(8, 0) = -4.5 * et * (6.0 * zt - 1.0);
        rResult(8, 1) = 4.5 * (zt * (3.0 * zt - 1.0) - et * (6.0 * zt - 1.0));
        //
        rResult(9, 0) = 27.0 * et * (zt - xi);
        rResult(9, 1) = 27.0 * xi * (zt - et);
        //
        return rResult;
    }

    /**
     * returns the second order derivatives of all shape functions
     * in given arbitrary pointers
     * @param rResult a third order tensor which contains the second derivatives
     * @param rPoint the given point the second order derivatives are calculated in
     */
    ShapeFunctionsSecondDerivativesType& ShapeFunctionsSecondDerivatives(
        ShapeFunctionsSecondDerivativesType& rResult, const CoordinatesArrayType& rPoint) const override
    {
        if (rResult.size() != this->PointsNumber())
        {
            ShapeFunctionsGradientsType temp(this->PointsNumber());
            rResult.swap(temp);
        }
        //
        rResult[0].resize(2, 2, false);
        rResult[1].resize(2, 2, false);
        rResult[2].resize(2, 2, false);
        rResult[3].resize(2, 2, false);
        rResult[4].resize(2, 2, false);
        rResult[5].resize(2, 2, false);
        rResult[6].resize(2, 2, false);
        rResult[7].resize(2, 2, false);
        rResult[8].resize(2, 2, false);
        rResult[9].resize(2, 2, false);
        //
        const double xi = rPoint[0];
        const double et = rPoint[1];
        const double zt = 1.0 - xi - et;
        //
        rResult[0](0, 0) = 9.0 * (3.0 * zt - 1.0);
        rResult[0](0, 1) = 9.0 * (3.0 * zt - 1.0);
        rResult[0](1, 0) = 9.0 * (3.0 * zt - 1.0);
        rResult[0](1, 1) = 9.0 * (3.0 * zt - 1.0);
        //
        rResult[1](0, 0) = 9.0 * (3.0 * xi - 1.0);
        rResult[1](0, 1) = 0.0;
        rResult[1](1, 0) = 0.0;
        rResult[1](1, 1) = 0.0;
        //
        rResult[2](0, 0) = 0.0;
        rResult[2](0, 1) = 0.0;
        rResult[2](1, 0) = 0.0;
        rResult[2](1, 1) = 9.0 * (3.0 * et - 1.0);
        //
        rResult[3](0, 0) = 9.0 * (3.0 * xi - 6.0 * zt + 1.0);
        rResult[3](0, 1) = 4.5 * (6.0 * xi - 6.0 * zt + 1.0);
        rResult[3](1, 0) = 4.5 * (6.0 * xi - 6.0 * zt + 1.0);
        rResult[3](1, 1) = 27.0 * xi;
        //
        rResult[4](0, 0) = 9.0 * (3.0 * zt - 6.0 * xi + 1.0);
        rResult[4](0, 1) = -4.5 * (6.0 * xi - 1.0);
        rResult[4](1, 0) = -4.5 * (6.0 * xi - 1.0);
        rResult[4](1, 1) = 0.0;
        //
        rResult[5](0, 0) = 27.0 * et;
        rResult[5](0, 1) = 4.5 * (6.0 * xi - 1.0);
        rResult[5](1, 0) = 4.5 * (6.0 * xi - 1.0);
        rResult[5](1, 1) = 0.0;
        //
        rResult[6](0, 0) = 0.0;
        rResult[6](0, 1) = 4.5 * (6.0 * et - 1.0);
        rResult[6](1, 0) = 4.5 * (6.0 * et - 1.0);
        rResult[6](1, 1) = 27.0 * xi;
        //
        rResult[7](0, 0) = 0.0;
        rResult[7](0, 1) = 4.5 * (6.0 * et - 1.0);
        rResult[7](1, 0) = 4.5 * (6.0 * et - 1.0);
        rResult[7](1, 1) = 9.0 * (3.0 * zt - 6.0 * et + 1.0);
        //
        rResult[8](0, 0) = 27.0 * et;
        rResult[8](0, 1) = 4.5 * (6.0 * et - 6.0 * zt + 1.0);
        rResult[8](1, 0) = 4.5 * (6.0 * et - 6.0 * zt + 1.0);
        rResult[8](1, 1) = 9.0 * (3.0 * et - 6.0 * zt + 1.0);
        //
        rResult[9](0, 0) = -54.0 * et;
        rResult[9](0, 1) = -27.0 * (xi - zt + et);
        rResult[9](1, 0) = -27.0 * (xi - zt + et);
        rResult[9](1, 1) = -54.0 * xi;
        //
        return rResult;
    }

    /**
     * returns the third order derivatives of all shape functions
     * in given arbitrary pointers
     * @param rResult a fourth order tensor which contains the third derivatives
     * @param rPoint the given point the third order derivatives are calculated in
     */
    ShapeFunctionsThirdDerivativesType& ShapeFunctionsThirdDerivatives(
        ShapeFunctionsThirdDerivativesType& rResult, const CoordinatesArrayType& rPoint) const override
    {
        if (rResult.size() != this->PointsNumber())
        {
            ShapeFunctionsThirdDerivativesType temp(this->PointsNumber());
            rResult.swap(temp);
        }
        //
        for (IndexType i = 0; i < rResult.size(); ++i)
        {
            DenseVector<Matrix> temp(this->PointsNumber());
            rResult[i].swap(temp);
        }
        //
        rResult[0][0].resize(2, 2, false);
        rResult[0][1].resize(2, 2, false);
        rResult[1][0].resize(2, 2, false);
        rResult[1][1].resize(2, 2, false);
        rResult[2][0].resize(2, 2, false);
        rResult[2][1].resize(2, 2, false);
        rResult[3][0].resize(2, 2, false);
        rResult[3][1].resize(2, 2, false);
        rResult[4][0].resize(2, 2, false);
        rResult[4][1].resize(2, 2, false);
        rResult[5][0].resize(2, 2, false);
        rResult[5][1].resize(2, 2, false);
        rResult[6][0].resize(2, 2, false);
        rResult[6][1].resize(2, 2, false);
        rResult[7][0].resize(2, 2, false);
        rResult[7][1].resize(2, 2, false);
        rResult[8][0].resize(2, 2, false);
        rResult[8][1].resize(2, 2, false);
        rResult[9][0].resize(2, 2, false);
        rResult[9][1].resize(2, 2, false);
        //
        rResult[0][0](0, 0) = -27.0;
        rResult[0][0](0, 1) = -27.0;
        rResult[0][0](1, 0) = -27.0;
        rResult[0][0](1, 1) = -27.0;
        rResult[0][1](0, 0) = -27.0;
        rResult[0][1](0, 1) = -27.0;
        rResult[0][1](1, 0) = -27.0;
        rResult[0][1](1, 1) = -27.0;
        //
        rResult[1][0](0, 0) = 27.0;
        rResult[1][0](0, 1) = 0.0;
        rResult[1][0](1, 0) = 0.0;
        rResult[1][0](1, 1) = 0.0;
        rResult[1][1](0, 0) = 0.0;
        rResult[1][1](0, 1) = 0.0;
        rResult[1][1](1, 0) = 0.0;
        rResult[1][1](1, 1) = 0.0;
        //
        rResult[2][0](0, 0) = 0.0;
        rResult[2][0](0, 1) = 0.0;
        rResult[2][0](1, 0) = 0.0;
        rResult[2][0](1, 1) = 0.0;
        rResult[2][1](0, 0) = 0.0;
        rResult[2][1](0, 1) = 0.0;
        rResult[2][1](1, 0) = 0.0;
        rResult[2][1](1, 1) = 27.0;
        //
        rResult[3][0](0, 0) = 81.0;
        rResult[3][0](0, 1) = 54.0;
        rResult[3][0](1, 0) = 54.0;
        rResult[3][0](1, 1) = 27.0;
        rResult[3][1](0, 0) = 54.0;
        rResult[3][1](0, 1) = 27.0;
        rResult[3][1](1, 0) = 27.0;
        rResult[3][1](1, 1) = 0.0;
        //
        rResult[4][0](0, 0) = -81.0;
        rResult[4][0](0, 1) = -27.0;
        rResult[4][0](1, 0) = -27.0;
        rResult[4][0](1, 1) = 0.0;
        rResult[4][1](0, 0) = -27.0;
        rResult[4][1](0, 1) = 0.0;
        rResult[4][1](1, 0) = 0.0;
        rResult[4][1](1, 1) = 0.0;
        //
        rResult[5][0](0, 0) = 0.0;
        rResult[5][0](0, 1) = 27.0;
        rResult[5][0](1, 0) = 27.0;
        rResult[5][0](1, 1) = 0.0;
        rResult[5][1](0, 0) = 27.0;
        rResult[5][1](0, 1) = 0.0;
        rResult[5][1](1, 0) = 0.0;
        rResult[5][1](1, 1) = 0.0;
        //
        rResult[6][0](0, 0) = 0.0;
        rResult[6][0](0, 1) = 0.0;
        rResult[6][0](1, 0) = 0.0;
        rResult[6][0](1, 1) = 27.0;
        rResult[6][1](0, 0) = 0.0;
        rResult[6][1](0, 1) = 27.0;
        rResult[6][1](1, 0) = 27.0;
        rResult[6][1](1, 1) = 0.0;
        //
        rResult[7][0](0, 0) = 0.0;
        rResult[7][0](0, 1) = 0.0;
        rResult[7][0](1, 0) = 0.0;
        rResult[7][0](1, 1) = -27.0;
        rResult[7][1](0, 0) = 0.0;
        rResult[7][1](0, 1) = -27.0;
        rResult[7][1](1, 0) = -27.0;
        rResult[7][1](1, 1) = -81.0;
        //
        rResult[8][0](0, 0) = 0.0;
        rResult[8][0](0, 1) = 27.0;
        rResult[8][0](1, 0) = 27.0;
        rResult[8][0](1, 1) = 54.0;
        rResult[8][1](0, 0) = 27.0;
        rResult[8][1](0, 1) = 54.0;
        rResult[8][1](1, 0) = 54.0;
        rResult[8][1](1, 1) = 81.0;
        //
        rResult[9][0](0, 0) = 0.0;
        rResult[9][0](0, 1) = -54.0;
        rResult[9][0](1, 0) = -54.0;
        rResult[9][0](1, 1) = -54.0;
        rResult[9][1](0, 0) = -54.0;
        rResult[9][1](0, 1) = -54.0;
        rResult[9][1](1, 0) = -54.0;
        rResult[9][1](1, 1) = 0.0;
        //
        return rResult;
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    /**
     * There are no protected members in class Triangle2D10
     */

private:
    ///@name Static Member Variables
    ///@{

    static const GeometryData msGeometryData;

    static const GeometryDimension msGeometryDimension;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
    }

    Triangle2D10() : BaseType(PointsArrayType(), &msGeometryData) {}

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
     * Calculates the values of all shape function in all integration points.
     * Integration points are expected to be given in local coordinates
     * @param ThisMethod the current integration method
     * @return the matrix of values of every shape function in each integration point
     * :KLUDGE: number of points is hard-coded -> be careful if you want to copy and paste!
     */
    static Matrix CalculateShapeFunctionsIntegrationPointsValues(typename BaseType::IntegrationMethod ThisMethod)
    {
        IntegrationPointsContainerType all_integration_points = AllIntegrationPoints();
        IntegrationPointsArrayType integration_points = all_integration_points[static_cast<int>(ThisMethod)];
        //number of integration points
        const int integration_points_number = integration_points.size();
        //number of nodes in current geometry
        const int points_number = 10;
        //setting up return matrix
        Matrix shape_function_values(integration_points_number, points_number);
        //loop over all integration points
        for (int pnt = 0; pnt < integration_points_number; ++pnt)
        {
            const double xi = integration_points[pnt].X();
            const double et = integration_points[pnt].Y();
            const double zt = 1.0 - xi - et;
            //
            shape_function_values(pnt, 0) = zt * (3.0 * zt - 1.0) * (3.0 * zt - 2.0) * 0.5;
            shape_function_values(pnt, 1) = xi * (3.0 * xi - 1.0) * (3.0 * xi - 2.0) * 0.5;
            shape_function_values(pnt, 2) = et * (3.0 * et - 1.0) * (3.0 * et - 2.0) * 0.5;
            shape_function_values(pnt, 3) = xi * zt * (3.0 * zt - 1.0) * 4.5;
            shape_function_values(pnt, 4) = xi * zt * (3.0 * xi - 1.0) * 4.5;
            shape_function_values(pnt, 5) = xi * et * (3.0 * xi - 1.0) * 4.5;
            shape_function_values(pnt, 6) = xi * et * (3.0 * et - 1.0) * 4.5;
            shape_function_values(pnt, 7) = et * zt * (3.0 * et - 1.0) * 4.5;
            shape_function_values(pnt, 8) = et * zt * (3.0 * zt - 1.0) * 4.5;
            shape_function_values(pnt, 9) = xi * et * zt * 27.0;
        }
        return shape_function_values;
    }

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
        CalculateShapeFunctionsIntegrationPointsLocalGradients(typename BaseType::IntegrationMethod ThisMethod)
    {
        IntegrationPointsContainerType all_integration_points = AllIntegrationPoints();
        IntegrationPointsArrayType integration_points = all_integration_points[static_cast<int>(ThisMethod)];
        //number of integration points
        const int integration_points_number = integration_points.size();
        ShapeFunctionsGradientsType d_shape_f_values(integration_points_number);
        //loop over all integration points
        for (int pnt = 0; pnt < integration_points_number; ++pnt)
        {
            Matrix result(10, 2);
            const double xi = integration_points[pnt].X();
            const double et = integration_points[pnt].Y();
            const double zt = 1.0 - xi - et;
            //
            noalias(result) = ZeroMatrix(10, 2);
            //
            result(0, 0) = -4.5 * zt * (3.0 * zt - 2.0) - 1.0;
            result(0, 1) = -4.5 * zt * (3.0 * zt - 2.0) - 1.0;
            //
            result(1, 0) = 4.5 * xi * (3.0 * xi - 2.0) + 1.0;
            result(1, 1) = 0.0;
            //
            result(2, 0) = 0.0;
            result(2, 1) = 4.5 * et * (3.0 * et - 2.0) + 1.0;
            //
            result(3, 0) = 4.5 * (zt * (3.0 * zt - 1.0) - xi * (6.0 * zt - 1.0));
            result(3, 1) = -4.5 * xi * (6.0 * zt - 1.0);
            //
            result(4, 0) = 4.5 * (zt * (6.0 * xi - 1.0) - xi * (3.0 * xi - 1.0));
            result(4, 1) = -4.5 * xi * (3.0 * xi - 1.0);
            //
            result(5, 0) = 4.5 * et * (6.0 * xi - 1.0);
            result(5, 1) = 4.5 * xi * (3.0 * xi - 1.0);
            //
            result(6, 0) = 4.5 * et * (3.0 * et - 1.0);
            result(6, 1) = 4.5 * xi * (6.0 * et - 1.0);
            //
            result(7, 0) = -4.5 * et * (3.0 * et - 1.0);
            result(7, 1) = 4.5 * (zt * (6.0 * et - 1.0) - et * (3.0 * et - 1.0));
            //
            result(8, 0) = -4.5 * et * (6.0 * zt - 1.0);
            result(8, 1) = 4.5 * (zt * (3.0 * zt - 1.0) - et * (6.0 * zt - 1.0));
            //
            result(9, 0) = 27.0 * et * (zt - xi);
            result(9, 1) = 27.0 * xi * (zt - et);
            //
            d_shape_f_values[pnt] = result;
        }
        return d_shape_f_values;
    }

    static const IntegrationPointsContainerType AllIntegrationPoints()
    {
        IntegrationPointsContainerType integration_points =
        {
            {
                Quadrature<TriangleGaussLegendreIntegrationPoints1, 2, IntegrationPoint<3>>::GenerateIntegrationPoints(),
                Quadrature<TriangleGaussLegendreIntegrationPoints2, 2, IntegrationPoint<3>>::GenerateIntegrationPoints(),
                Quadrature<TriangleGaussLegendreIntegrationPoints3, 2, IntegrationPoint<3>>::GenerateIntegrationPoints(),
                Quadrature<TriangleGaussLegendreIntegrationPoints4, 2, IntegrationPoint<3>>::GenerateIntegrationPoints(),
                Quadrature<TriangleGaussLegendreIntegrationPoints5, 2, IntegrationPoint<3>>::GenerateIntegrationPoints()
            }
        };
        return integration_points;
    }

    static const ShapeFunctionsValuesContainerType AllShapeFunctionsValues()
    {
        ShapeFunctionsValuesContainerType shape_functions_values =
        {
            {
                Triangle2D10<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(GeometryData::IntegrationMethod::GI_GAUSS_1),
                Triangle2D10<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(GeometryData::IntegrationMethod::GI_GAUSS_2),
                Triangle2D10<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(GeometryData::IntegrationMethod::GI_GAUSS_3),
                Triangle2D10<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(GeometryData::IntegrationMethod::GI_GAUSS_4),
                Triangle2D10<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(GeometryData::IntegrationMethod::GI_GAUSS_5),
            }
        };
        return shape_functions_values;
    }

    static const ShapeFunctionsLocalGradientsContainerType AllShapeFunctionsLocalGradients()
    {
        ShapeFunctionsLocalGradientsContainerType shape_functions_local_gradients =
        {
            {
                Triangle2D10<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(GeometryData::IntegrationMethod::GI_GAUSS_1),
                Triangle2D10<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(GeometryData::IntegrationMethod::GI_GAUSS_2),
                Triangle2D10<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(GeometryData::IntegrationMethod::GI_GAUSS_3),
                Triangle2D10<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(GeometryData::IntegrationMethod::GI_GAUSS_4),
                Triangle2D10<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(GeometryData::IntegrationMethod::GI_GAUSS_5),
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

    template<class TOtherPointType> friend class Triangle2D10;

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

/**
 * input stream functions
 */
template<class TPointType>
inline std::istream& operator >> (std::istream& rIStream, Triangle2D10<TPointType>& rThis);

/**
 * output stream functions
 */
template<class TPointType>
inline std::ostream& operator << (std::ostream& rOStream, const Triangle2D10<TPointType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

///@}

template<class TPointType> const
GeometryData Triangle2D10<TPointType>::msGeometryData(
    &msGeometryDimension,
    GeometryData::IntegrationMethod::GI_GAUSS_4,
    Triangle2D10<TPointType>::AllIntegrationPoints(),
    Triangle2D10<TPointType>::AllShapeFunctionsValues(),
    AllShapeFunctionsLocalGradients());

template<class TPointType>
const GeometryDimension Triangle2D10<TPointType>::msGeometryDimension(2, 2);

}// namespace Kratos.