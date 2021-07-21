//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Manuel Messmer
//

#if !defined(KRATOS_SURFACE_IN_NURBS_VOLUME_INCLUDE_H )
#define  KRATOS_SURFACE_IN_NURBS_VOLUME_INCLUDE_H

// Project includes
#include "geometries/geometry.h"

#include "geometries/nurbs_volume_geometry.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_interval.h"

#include "utilities/nurbs_utilities/projection_nurbs_geometry_utilities.h"

namespace Kratos {

template <int TWorkingSpaceDimension, class TVolumeContainerPointType>
class SurfaceInNurbsVolumeGeometry : public Geometry<typename TVolumeContainerPointType::value_type>
{
public:
    ///@name Type Definitions
    ///@{

    typedef typename TVolumeContainerPointType::value_type NodeType;

    typedef Geometry<NodeType> GeometryType;
    typedef typename GeometryType::Pointer GeometryPointerType;

    typedef Geometry<NodeType> BaseType;

    typedef typename BaseType::IndexType IndexType;
    typedef typename BaseType::SizeType SizeType;

    typedef NurbsVolumeGeometry<TVolumeContainerPointType> NurbsVolumeType;
    typedef typename NurbsVolumeType::Pointer NurbsVolumePointerType;

    typedef typename BaseType::CoordinatesArrayType CoordinatesArrayType;
    typedef typename BaseType::PointsArrayType PointsArrayType;
    typedef typename BaseType::GeometriesArrayType GeometriesArrayType;
    typedef typename BaseType::IntegrationPointsArrayType IntegrationPointsArrayType;

    // using base class functionalities.
    using BaseType::CreateQuadraturePointGeometries;
    using BaseType::pGetPoint;
    using BaseType::GetPoint;

    static constexpr IndexType VOLUME_INDEX = -1;

    /// Counted pointer of SurfaceInNurbsVolumeGeometry
    KRATOS_CLASS_POINTER_DEFINITION(SurfaceInNurbsVolumeGeometry);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    SurfaceInNurbsVolumeGeometry(
        NurbsVolumePointerType pVolume,
        GeometryPointerType pSurface)
        : BaseType(PointsArrayType(), &msGeometryData)
        , mpNurbsVolume(pVolume)
        , mpSurface(pSurface)
    {
    }

    /// Default constructor
    SurfaceInNurbsVolumeGeometry()
        : BaseType(PointsArrayType(), &msGeometryData)
    {};

    /// Copy constructor
    SurfaceInNurbsVolumeGeometry(SurfaceInNurbsVolumeGeometry const& rOther)
        : BaseType(rOther)
        , mpNurbsVolume(rOther.mpNurbsVolume)
        , mpSurface(rOther.mpSurface)
    {
    }

    /// Copy constructor, with different point type.
    template<class TOtherVolumeContainerPointType>
    SurfaceInNurbsVolumeGeometry(SurfaceInNurbsVolumeGeometry<TWorkingSpaceDimension, TOtherVolumeContainerPointType> const& rOther)
        : BaseType(rOther, &msGeometryData)
        , mpNurbsVolume(rOther.mpNurbsVolume)
        , mpSurface(rOther.mpSurface)
    {
    }

    /// Destructor
    ~SurfaceInNurbsVolumeGeometry() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator
    SurfaceInNurbsVolumeGeometry& operator=(const SurfaceInNurbsVolumeGeometry& rOther)
    {
        BaseType::operator=(rOther);
        mpNurbsVolume = rOther.mpNurbsVolume;
        mpSurface = rOther.mpSurface;
        return *this;
    }

    /// @brief Assignment operator for geometries with different point type.
    template<class TOtherVolumeContainerPointType>
    SurfaceInNurbsVolumeGeometry& operator=(
        SurfaceInNurbsVolumeGeometry<TWorkingSpaceDimension, TOtherVolumeContainerPointType> const & rOther)
    {
        BaseType::operator=(rOther);
        mpNurbsVolume = rOther.mpNurbsVolume;
        mpSurface = rOther.mpSurface;
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    typename BaseType::Pointer Create(
        TVolumeContainerPointType const& ThisPoints) const override
    {
        KRATOS_ERROR << "SurfaceInNurbsVolumeGeometry cannot be created with 'PointsArrayType const& ThisPoints'. "
            << "'Create' is not allowed as it would not contain the required pointers to the surface and the curve."
            << std::endl;
    }

    ///@}
    ///@name Access to Geometry Parts
    ///@{

    /**
    * @brief This function returns the pointer of the geometry
    *        which is corresponding to the BACKGROUND_GEOMETRY_INDEX.
    * @param Index: BACKGROUND_GEOMETRY_INDEX.
    * @return pointer of geometry, corresponding to the index.
    */
    GeometryPointerType pGetGeometryPart(const IndexType Index) override
    {
        const auto& const_this = *this;
        return std::const_pointer_cast<GeometryType>(
            const_this.pGetGeometryPart(Index));
    }

    /**
    * @brief This function returns the pointer of the geometry
    *        which is corresponding to the BACKGROUND_GEOMETRY_INDEX.
    * @param Index: BACKGROUND_GEOMETRY_INDEX.
    * @return pointer of geometry, corresponding to the index.
    */
    const GeometryPointerType pGetGeometryPart(const IndexType Index) const override
    {
        if (Index == GeometryType::BACKGROUND_GEOMETRY_INDEX)
            return mpNurbsVolume;

        KRATOS_ERROR << "Index " << Index << " not existing in SurfaceInNurbsVolume: "
            << this->Id() << std::endl;
    }

    /**
    * @brief This function is used to check if the index is
    *        GeometryType::BACKGROUND_GEOMETRY_INDEX.
    * @param Index of the geometry part.
    * @return true if GeometryType::BACKGROUND_GEOMETRY_INDEX.
    */
    bool HasGeometryPart(const IndexType Index) const override
    {
        return Index == GeometryType::BACKGROUND_GEOMETRY_INDEX;
    }

    ///@}
    ///@name Mathematical Informations
    ///@{

    /// Return polynomial degree of the embedded surface
    SizeType PolynomialDegree(IndexType LocalDirectionIndex) const override
    {
        return mpSurface->PolynomialDegree(LocalDirectionIndex);
    }

    ///@}
    ///@name Set/ Get functions
    ///@{

    /// Returns the const NurbsSurface::Pointer of this brep.
    GeometryPointerType pGetSurface()
    {
        return mpSurface;
    }

    /// Returns the const NurbsSurface::Pointer of this brep.
    const GeometryPointerType pGetSurface() const
    {
        return mpSurface;
    }

    ///@}
    ///@name Geometrical Informations
    ///@{

    /// Provides the center of the underlying surface
    Point Center() const override
    {
        return mpSurface->Center();
    }

    /// Computes the area of the underlying surface
    double Area() const override
    {
        return mpSurface->Area();
    }

    ///@}
    ///@name Jacobian
    ///@{

    double DeterminantOfJacobian(
        const CoordinatesArrayType& rPoint) const override
    {
        std::vector<CoordinatesArrayType> global_space_derivatives(3);
        this->GlobalSpaceDerivatives(
            global_space_derivatives, rPoint, 1);

        BoundedMatrix<double,3,2> global_tangents;

        column(global_tangents,0) = global_space_derivatives[1];
        column(global_tangents,1) = global_space_derivatives[2];

        return MathUtils<double>::GeneralizedDet( global_tangents );
    }

    ///@}
    ///@name Quadrature Point Geometries
    ///@{

    /**
     * @brief This method creates a list of quadrature point geometries
     *        from a list of integration points.
     *
     * @param rResultGeometries list of quadrature point geometries.
     * @param rDummyIntegrationPoints list of integration points.
     * @param NumberOfShapeFunctionDerivatives the number provided
     *        derivatives of shape functions in the system.
     *
     * @see quadrature_point_geometry.h
     */
    void CreateQuadraturePointGeometries(
        GeometriesArrayType& rResultGeometries,
        IndexType NumberOfShapeFunctionDerivatives,
        const IntegrationPointsArrayType& rDummyIntegrationPoints) override
    {
        // shape function container.
        NurbsVolumeShapeFunction shape_function_container(
            mpNurbsVolume->PolynomialDegreeU(), mpNurbsVolume->PolynomialDegreeV(), mpNurbsVolume->PolynomialDegreeW(),
            NumberOfShapeFunctionDerivatives);

        auto default_method = this->GetDefaultIntegrationMethod();
        SizeType num_nonzero_cps = shape_function_container.NumberOfNonzeroControlPoints();

        Matrix N(1, num_nonzero_cps);
        DenseVector<Matrix> shape_function_derivatives(NumberOfShapeFunctionDerivatives - 1);

        for (IndexType i = 0; i < NumberOfShapeFunctionDerivatives - 1; ++i) {
            const IndexType num_derivatives = (2 + i) * (3 + i) / 2;
            shape_function_derivatives[i].resize(num_nonzero_cps, num_derivatives);
        }

        // Get integration points from surface
        const SizeType number_of_points = mpSurface->IntegrationPointsNumber();
        auto surface_method = mpSurface->GetDefaultIntegrationMethod();
        const IntegrationPointsArrayType& integration_points = mpSurface->IntegrationPoints(surface_method);

        // Resize containers.
        if (rResultGeometries.size() != integration_points.size() ){
            rResultGeometries.resize(integration_points.size());
        }

        for (IndexType point_index = 0; point_index < number_of_points; ++point_index)
        {

            array_1d<double,3> global_coordinates;
            mpSurface->GlobalCoordinates(global_coordinates, integration_points[point_index]);

            shape_function_container.ComputeBSplineShapeFunctionValues(
                mpNurbsVolume->KnotsU(), mpNurbsVolume->KnotsV(),  mpNurbsVolume->KnotsW(),
                global_coordinates[0], global_coordinates[1], global_coordinates[2]);

            /// Get List of Control Points
            PointsArrayType nonzero_control_points(num_nonzero_cps);
            auto cp_indices = shape_function_container.ControlPointIndices(
                mpNurbsVolume->NumberOfControlPointsU(), mpNurbsVolume->NumberOfControlPointsV(), mpNurbsVolume->NumberOfControlPointsW());
            for (IndexType j = 0; j < num_nonzero_cps; j++) {
                nonzero_control_points(j) = mpNurbsVolume->pGetPoint(cp_indices[j]);
            }
            /// Get Shape Functions N
            if (NumberOfShapeFunctionDerivatives >= 0) {
                for (IndexType j = 0; j < num_nonzero_cps; j++) {
                    N(0, j) = shape_function_container(j, 0);
                }
            }

            /// Get Shape Function Derivatives DN_De, ...
            if (NumberOfShapeFunctionDerivatives > 0) {
                IndexType shape_derivative_index = 1;
                for (IndexType n = 0; n < NumberOfShapeFunctionDerivatives - 1; ++n) {
                    const IndexType num_derivatives = (2 + n) * (3 + n) / 2;
                    for (IndexType k = 0; k < num_derivatives; ++k) {
                        for (IndexType j = 0; j < num_nonzero_cps; ++j) {
                            shape_function_derivatives[n](j, k) = shape_function_container(j, shape_derivative_index + k);
                        }
                    }
                    shape_derivative_index += num_derivatives;
                }
            }

            const double determinant_volume = mpNurbsVolume->DeterminantOfJacobian(global_coordinates);
            // Get area of surface in global space
            const double weight = mpSurface->Area() * determinant_volume;
            IntegrationPoint<3> tmp_integration_point(global_coordinates[0], global_coordinates[1], global_coordinates[2], weight);

            GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> data_container(
                default_method, tmp_integration_point,
                N, shape_function_derivatives);

            Matrix jacobian_surface;
            mpSurface->Jacobian(jacobian_surface, integration_points[point_index]);
            BoundedMatrix<double,3,2> local_tangents;
            local_tangents = jacobian_surface;

            rResultGeometries(point_index) = CreateQuadraturePointsUtility<NodeType>::CreateQuadraturePointSurfaceInVolume(
                data_container, nonzero_control_points,
                local_tangents, this);
        }
    }


    ///@}
    ///@name Operation within Global Space
    ///@{

    /*
    * @brief This method maps from dimension space to working space.
    * From Piegl and Tiller, The NURBS Book, Algorithm A3.1/ A4.1
    * @param rResult array_1d<double, 3> with the coordinates in working space
    * @param LocalCoordinates The local coordinates in dimension space
    * @return array_1d<double, 3> with the coordinates in working space
    * @see PointLocalCoordinates
    */
    CoordinatesArrayType& GlobalCoordinates(
        CoordinatesArrayType& rResult,
        const CoordinatesArrayType& rLocalCoordinates
    ) const override
    {
        // Compute the coordinates of the embedded curve in the parametric space of the surface
        CoordinatesArrayType result_local = mpSurface->GlobalCoordinates(rResult, rLocalCoordinates);

        // Compute and return the coordinates of the surface in the geometric space
        return mpNurbsVolume->GlobalCoordinates(rResult, result_local);
    }

/**
    * @brief This method maps from dimension space to working space and computes the
    *        number of derivatives at the dimension parameter.
    * @param rCoordinates The local coordinates in dimension space
    * @param DerivativeOrder Number of computed derivatives
    * @return std::vector<array_1d<double, 3>> with the coordinates in working space
    * @see PointLocalCoordinates
    */
    void GlobalSpaceDerivatives(
        std::vector<CoordinatesArrayType>& rGlobalSpaceDerivatives,
        const CoordinatesArrayType& rCoordinates,
        const SizeType DerivativeOrder) const override
    {

        KRATOS_ERROR_IF( DerivativeOrder > 1 )  << "Surface in Nurbs Volume: Global Space Derivatives are not yet "
            << "implemented for derivative order > 1." << std::endl;

        // Check size of output
        if (rGlobalSpaceDerivatives.size() != 3) {
            rGlobalSpaceDerivatives.resize(3);
        }

        // Compute the gradients of the embedded curve in the parametric space of the surface
        array_1d<double,3> global_coordinates;
        mpSurface->GlobalCoordinates(global_coordinates, rCoordinates);

        // Compute the gradients of the volume in the geometric space
        std::vector<CoordinatesArrayType> volume_derivatives(1 + DerivativeOrder*3);
        mpNurbsVolume->GlobalSpaceDerivatives(volume_derivatives, global_coordinates, DerivativeOrder);

        rGlobalSpaceDerivatives[0] = global_coordinates;

        BoundedMatrix<double,3,3> volume_jacobian;
        noalias(column(volume_jacobian,0)) =  volume_derivatives[1];
        noalias(column(volume_jacobian,1)) =  volume_derivatives[2];
        noalias(column(volume_jacobian,2)) =  volume_derivatives[3];

        Matrix jacobian_surface;
        mpSurface->Jacobian(jacobian_surface, rCoordinates);
        BoundedMatrix<double,3,2> surface_jacobian;
        surface_jacobian = jacobian_surface;

        BoundedMatrix<double,3,2> global_tangents = prod(volume_jacobian, surface_jacobian);
        rGlobalSpaceDerivatives[1] = column(global_tangents,0);
        rGlobalSpaceDerivatives[2] = column(global_tangents,1);


    }

    ///@}
    ///@name Information
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "2 dimensional nurbs curve on 3D surface.";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "2 dimensional nurbs curve on 3D surface.";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
private:
    ///@name Private Static Member Variables
    ///@{

    static const GeometryData msGeometryData;

    static const GeometryDimension msGeometryDimension;

    ///@}
    ///@name Private Member Variables
    ///@{

    NurbsVolumePointerType mpNurbsVolume;
    GeometryPointerType mpSurface;

    ///@}
    ///@name Private Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
        rSerializer.save("pNurbsVolume", mpNurbsVolume);
        rSerializer.save("pSurface", mpSurface);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
        rSerializer.load("pNurbsSurface", mpNurbsVolume);
        rSerializer.load("pSurface", mpSurface);
    }

    ///@}

}; // class SurfaceInNurbsVolumeGeometry

template<int TWorkingSpaceDimension, class TVolumeContainerPointType>
const GeometryData SurfaceInNurbsVolumeGeometry<TWorkingSpaceDimension, TVolumeContainerPointType>::msGeometryData(
    &msGeometryDimension,
    GeometryData::GI_GAUSS_1,
    {}, {}, {});

template<int TWorkingSpaceDimension, class TVolumeContainerPointType>
const GeometryDimension SurfaceInNurbsVolumeGeometry<TWorkingSpaceDimension, TVolumeContainerPointType>::msGeometryDimension(
    2, TWorkingSpaceDimension, 3);

} // namespace Kratos

#endif // KRATOS_SURFACE_IN_NURBS_VOLUME_INCLUDE_H defined
