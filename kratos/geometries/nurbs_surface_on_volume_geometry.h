//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Andrea Gorgi Nicoló Antonelli
//

#if !defined(KRATOS_NURBS_SURFACE_ON_VOLUME_INCLUDED_H )
#define  KRATOS_NURBS_SURFACE_ON_VOLUME_INCLUDED_H

// Project includes
#include "geometries/geometry.h"

#include "geometries/nurbs_volume_geometry.h"
#include "geometries/nurbs_surface_geometry.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_interval.h"

namespace Kratos {

template <int TWorkingSpaceDimension, class TVolumeContainerPointType>
class NurbsSurfaceOnVolumeGeometry : public Geometry<typename TVolumeContainerPointType::value_type>
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

    /// Counted pointer of NurbsSurfaceOnVolumeGeometry
    KRATOS_CLASS_POINTER_DEFINITION(NurbsSurfaceOnVolumeGeometry);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    NurbsSurfaceOnVolumeGeometry(
        NurbsVolumePointerType pVolume,
        GeometryPointerType pSurface)
        : BaseType(PointsArrayType(), &msGeometryData)
        , mpNurbsVolume(pVolume)
        , mpSurface(pSurface)
    {
    }

    /// Default constructor
    NurbsSurfaceOnVolumeGeometry()
        : BaseType(PointsArrayType(), &msGeometryData)
    {};

    /// Copy constructor
    NurbsSurfaceOnVolumeGeometry(NurbsSurfaceOnVolumeGeometry const& rOther)
        : BaseType(rOther)
        , mpNurbsVolume(rOther.mpNurbsVolume)
        , mpSurface(rOther.mpSurface)
    {
    }

    /// Copy constructor, with different point type.
    template<class TOtherVolumeContainerPointType>
    NurbsSurfaceOnVolumeGeometry(NurbsSurfaceOnVolumeGeometry<TWorkingSpaceDimension, TOtherVolumeContainerPointType> const& rOther)
        : BaseType(rOther, &msGeometryData)
        , mpNurbsVolume(rOther.mpNurbsVolume)
        , mpSurface(rOther.mpSurface)
    {
    }

    /// Destructor
    ~NurbsSurfaceOnVolumeGeometry() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator
    NurbsSurfaceOnVolumeGeometry& operator=(const NurbsSurfaceOnVolumeGeometry& rOther)
    {
        BaseType::operator=(rOther);
        mpNurbsVolume = rOther.mpNurbsVolume;
        mpSurface = rOther.mpSurface;
        return *this;
    }

    /// @brief Assignment operator for geometries with different point type.
    template<class TOtherVolumeContainerPointType>
    NurbsSurfaceOnVolumeGeometry& operator=(
        NurbsSurfaceOnVolumeGeometry<TWorkingSpaceDimension, TOtherVolumeContainerPointType> const & rOther)
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
        KRATOS_ERROR << "NurbsSurfaceOnVolumeGeometry cannot be created with 'PointsArrayType const& ThisPoints'. "
            << "'Create' is not allowed as it would not contain the required pointers to the surface and the nurbs volume."
            << std::endl;
    }

    ///@}
    ///@name Access to Geometry Parts
    ///@{

    /**
    * @brief This function returns the pointer of the geometry
    *        which is corresponding to the BACKGROUND_GEOMETRY_INDEX.
    * @param Index: BACKGROUND_GEOMETRY_INDEX.
    * @return PPointer of geometry, corresponding to the index.
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
    * @return Pointer of geometry, corresponding to the index.
    */
    const GeometryPointerType pGetGeometryPart(const IndexType Index) const override
    {
        if (Index == GeometryType::BACKGROUND_GEOMETRY_INDEX)
            return mpNurbsVolume;

        KRATOS_ERROR << "Index " << Index << " not existing in NurbsSurfaceOnVolume with ID: "
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
    ///@name Set / Calculate access
    ///@{

    void Calculate(
        const Variable<array_1d<double, 3>>& rVariable,
        array_1d<double, 3>& rOutput
        ) const override
    {
        if (rVariable == PARAMETER_2D_COORDINATES) {
            array_1d<double, 3>& local_parameter = rOutput;
            mpSurface->GlobalCoordinates(rOutput, local_parameter);
        }
    }

    ///@}
    ///@name Mathematical Informations
    ///@{

    /// Return polynomial degree of the curve
    SizeType PolynomialDegree(IndexType LocalDirectionIndex) const override
    {
        return mpSurface->PolynomialDegree(0) + mpSurface->PolynomialDegree(1);
    }

    ///@}
    ///@name Set/ Get functions
    ///@{

    /// Returns the NurbsSurface::Pointer of this brep.
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

    /// Provides the center of the underlying surface.
    Point Center() const override
    {
        Point local_center = mpSurface->Center();
        Point global_center;
        mpNurbsVolume->GlobalCoordinates(global_center, local_center);

        return global_center;
    }

    /// Computes the area of the underlying surface.
    double Area() const override
    {
        IntegrationPointsArrayType integration_points;
        IntegrationInfo integration_info = this->GetDefaultIntegrationInfo();
        CreateIntegrationPoints(integration_points, integration_info);

        double area = 0.0;
        for (IndexType i = 0; i < integration_points.size(); ++i) {
            const double determinant_jacobian = DeterminantOfJacobian(integration_points[i]);
            area += integration_points[i].Weight() * determinant_jacobian;
        }

        return area;
    }

    ///@}
    ///@name Jacobian
    ///@{

    double DeterminantOfJacobian(
        const CoordinatesArrayType& rPoint) const override
    {
        std::vector<CoordinatesArrayType> global_space_derivatives(3);
        this->GlobalSpaceDerivatives( global_space_derivatives, rPoint, 1);

        BoundedMatrix<double,3,2> global_tangents;

        column(global_tangents,0) = global_space_derivatives[1];
        column(global_tangents,1) = global_space_derivatives[2];

        return MathUtils<double>::GeneralizedDet( global_tangents );
    }

    ///@}
    ///@name Integration Info
    ///@{

    /// Provides the default integration dependent on the polynomial degree.
    IntegrationInfo GetDefaultIntegrationInfo() const override
    {
        
        // IndexType p = mpNurbsVolume->PolynomialDegreeU() + mpNurbsVolume->PolynomialDegreeV() +1;

        return IntegrationInfo(
            { mpNurbsVolume->PolynomialDegreeU() + 1, mpNurbsVolume->PolynomialDegreeV() + 1 },
            { IntegrationInfo::QuadratureMethod::GAUSS, IntegrationInfo::QuadratureMethod::GAUSS });
    }

    // IntegrationInfo GetDefaultIntegrationInfo() const override
    // {
    //     IndexType p = mpNurbsVolume->PolynomialDegreeU() + mpNurbsVolume->PolynomialDegreeV() + mpNurbsVolume->PolynomialDegreeW() + 1;
    //     return IntegrationInfo(2, p, IntegrationInfo::QuadratureMethod::GAUSS);
    // }


    ///@}
    ///@name Quadrature Point Geometries
    ///@{

    /**
     * @brief Returns the integration points of the surface.
     * @param result integration points.
     **/
    void CreateIntegrationPoints(
        IntegrationPointsArrayType& rIntegrationPoints,
        IntegrationInfo& rIntegrationInfo) const override
    {
        const auto surface_method = mpSurface->GetDefaultIntegrationMethod();
        rIntegrationPoints = mpSurface->IntegrationPoints(surface_method);
    }

    /**
     * @brief This method creates a list of quadrature point geometries
     *        from a list of integration points.
     *
     * @param rResultGeometries List of quadrature point geometries.
     * @param rIntegrationPoints List of integration points.
     * @param NumberOfShapeFunctionDerivatives The number provided
     *        derivatives of shape functions in the system.
     *
     * @see quadrature_point_geometry.h
     */
    void CreateQuadraturePointGeometries(
        GeometriesArrayType& rResultGeometries,
        IndexType NumberOfShapeFunctionDerivatives,
        const IntegrationPointsArrayType& rIntegrationPoints,
        IntegrationInfo& rIntegrationInfo) override
    {
        // Shape function container.
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
        const SizeType number_of_points = rIntegrationPoints.size();

        // Resize containers.
        if (rResultGeometries.size() != rIntegrationPoints.size() ){
            rResultGeometries.resize(rIntegrationPoints.size());
        }

        for (IndexType point_index = 0; point_index < number_of_points; ++point_index)
        {

            array_1d<double,3> global_coordinates;
            mpSurface->GlobalCoordinates(global_coordinates, rIntegrationPoints[point_index]);

            if (isExitingDirectionSBM > -1) {
                // SURROGATE BOUNDARY Gauss Point which is EXITING)
                IndexType SpanU = NurbsUtilities::GetLowerSpan(mpNurbsVolume->PolynomialDegreeU(), mpNurbsVolume->KnotsU(), global_coordinates[0]);
                IndexType SpanV = NurbsUtilities::GetLowerSpan(mpNurbsVolume->PolynomialDegreeV(), mpNurbsVolume->KnotsV(), global_coordinates[1]);
                IndexType SpanW = NurbsUtilities::GetLowerSpan(mpNurbsVolume->PolynomialDegreeW(), mpNurbsVolume->KnotsW(), global_coordinates[2]);

                // Understand if the knot-boundary is along U or V
                switch (isExitingDirectionSBM) {
                    case 0: {
                        SpanU++;
                        break;
                    }
                    case 1: {
                        SpanV++;
                        break;
                    }
                    case 2: {
                        SpanW++;
                        break;
                    }
                    default: {
                        KRATOS_ERROR << "WRONG EXITING DIRECTION: not external face and not parallel to any direction" << std::endl;
                    }
                }

                shape_function_container.ComputeBSplineShapeFunctionValuesAtSpan(
                            mpNurbsVolume->KnotsU(), mpNurbsVolume->KnotsV(),  mpNurbsVolume->KnotsW(), SpanU, SpanV, SpanW,
                            global_coordinates[0], global_coordinates[1], global_coordinates[2]);
            } else {
                // The Gauss point is ENTERING or it is an external one
                shape_function_container.ComputeBSplineShapeFunctionValues(
                                        mpNurbsVolume->KnotsU(), mpNurbsVolume->KnotsV(),  mpNurbsVolume->KnotsW(),
                                        global_coordinates[0], global_coordinates[1], global_coordinates[2]);
            }
            

            /// Get List of Control Pointsx
            PointsArrayType nonzero_control_points(num_nonzero_cps);
            auto cp_indices = shape_function_container.ControlPointIndices(
                mpNurbsVolume->NumberOfControlPointsU(), mpNurbsVolume->NumberOfControlPointsV(), mpNurbsVolume->NumberOfControlPointsW());
            for (IndexType j = 0; j < num_nonzero_cps; j++) {
                nonzero_control_points(j) = mpNurbsVolume->pGetPoint(cp_indices[j]);
            }

            /// Get Shape Functions N
            for (IndexType j = 0; j < num_nonzero_cps; j++) {
                N(0, j) = shape_function_container(j, 0);
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

            GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> data_container(
                default_method, rIntegrationPoints[point_index],
                N, shape_function_derivatives);

            Matrix jacobian_surface;
            mpSurface->Jacobian(jacobian_surface, rIntegrationPoints[point_index]);

            // Feature to pass the normal at each integration point
            if (msNormal.size()==0) KRATOS_ERROR << "[NURBS_SURFACE_ON_VOLUME]:: Normal not defined" << std::endl;
            rResultGeometries(point_index) = CreateQuadraturePointsUtility<NodeType>::CreateQuadraturePointSurfaceInVolume(
                data_container, nonzero_control_points,
                jacobian_surface, msNormal, this);
        }
    }

    ///@}
    ///@name Operation within Global Space
    ///@{

    /**
    * @brief This method maps from dimension space to working space.
    * @param rResult array_1d<double, 3> with the coordinates in working space
    * @param LocalCoordinates The local coordinates in dimension space
    * @return array_1d<double, 3> with the coordinates in working space
    * @see PointLocalCoordinates
    **/
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

        KRATOS_ERROR_IF( DerivativeOrder > 1 )  << "NurbsSurfaceOnVolume :: Global Space Derivatives are not yet "
            << "implemented for derivative orders > 1." << std::endl;

        // Check size of output
        if (rGlobalSpaceDerivatives.size() != 3) {
            rGlobalSpaceDerivatives.resize(3);
        }

        // Map rCoordinates into local space.
        array_1d<double,3> global_coordinates;
        mpSurface->GlobalCoordinates(global_coordinates, rCoordinates);

        // Compute the gradients of the volume in the geometric space
        std::vector<CoordinatesArrayType> volume_derivatives(1 + DerivativeOrder*3);
        mpNurbsVolume->GlobalSpaceDerivatives(volume_derivatives, global_coordinates, DerivativeOrder);

        rGlobalSpaceDerivatives[0] = global_coordinates;

        BoundedMatrix<double,3,3> volume_jacobian;
        column(volume_jacobian,0) =  volume_derivatives[1];
        column(volume_jacobian,1) =  volume_derivatives[2];
        column(volume_jacobian,2) =  volume_derivatives[3];

        Matrix jacobian_surface;
        mpSurface->Jacobian(jacobian_surface, rCoordinates);

        BoundedMatrix<double,3,2> global_tangents = prod(volume_jacobian, jacobian_surface);
        rGlobalSpaceDerivatives[1] = column(global_tangents,0);
        rGlobalSpaceDerivatives[2] = column(global_tangents,1);
    }


    /* @brief  Provides intersections of the nurbs curve with the knots of the surface.
     * @return vector of interval limitations.
     */
    void SpansLocalSpace(std::vector<double>& rSpansU, std::vector<double>& rSpansV,
        double Start_U, double End_U, double Start_V, double End_V) const
    {
        std::vector<double> volume_spans_u;
        std::vector<double> volume_spans_v;
        std::vector<double> volume_spans_w;
        mpNurbsVolume->SpansLocalSpace(volume_spans_u, 0);
        mpNurbsVolume->SpansLocalSpace(volume_spans_v, 1);
        mpNurbsVolume->SpansLocalSpace(volume_spans_w, 2);

        CoordinatesArrayType physical_coord_lower_left = ZeroVector(3);
        CoordinatesArrayType local_coord_lower_left = ZeroVector(3);
        local_coord_lower_left[0] = Start_U;
        local_coord_lower_left[1] = Start_V;
        mpSurface->GlobalCoordinates(physical_coord_lower_left, local_coord_lower_left);

        CoordinatesArrayType physical_coord_lower_right = ZeroVector(3);
        CoordinatesArrayType local_coord_lower_right = ZeroVector(3);
        local_coord_lower_right[0] = End_U;
        local_coord_lower_right[1] = Start_V;
        // KRATOS_WATCH(physical_coord_lower_right)
        mpSurface->GlobalCoordinates(physical_coord_lower_right, local_coord_lower_right);
        // KRATOS_WATCH(physical_coord_lower_right)

        CoordinatesArrayType physical_coord_upper_left = ZeroVector(3);
        CoordinatesArrayType local_coord_upper_left = ZeroVector(3);
        local_coord_upper_left[0] = Start_U;
        local_coord_upper_left[1] = End_V;
        // KRATOS_WATCH(local_coord_lower_left)
        mpSurface->GlobalCoordinates(physical_coord_upper_left, local_coord_upper_left);
        // KRATOS_WATCH(physical_coord_upper_left)


        // Compute the intersection rSpansU and rSpansV

        // LOWER SEGMENT
        const double toll = 1e-10;
        const double toll2 = 1e-15;
        
        // Compute constant coordinate -> understand where the face is in volume parameter space

        // Scale factor between volume coordinate and surface one
        double physical_length_segment = norm_2(physical_coord_lower_left-physical_coord_lower_right);
        double parameter_length_segment = norm_2(local_coord_lower_left-local_coord_lower_right);
        double scale_factor = parameter_length_segment/physical_length_segment;

        std::vector<double> knot_interval(2);
        // std::vector<double> intersection_knot_values_in_volume_lower ; // RSPANSU
        if (std::abs(physical_coord_lower_left[0]-physical_coord_lower_right[0]) > toll) {// left or right face
            // direction_lower_segment = 0; // compare to knot_spans_vector_u
            knot_interval[0] = physical_coord_lower_left[0]-toll2; 
            knot_interval[1] = physical_coord_lower_right[0]+toll2;
            std::sort(knot_interval.begin(), knot_interval.end());
            // Compare with volume_spans_u
            for (int i = 0; i < volume_spans_u.size(); i++) {
                double curr_knot_value = volume_spans_u[i];
                if (curr_knot_value < knot_interval[0]) {continue;}
                if (std::abs(curr_knot_value - knot_interval[0]) < toll2*10) knot_interval[0] = curr_knot_value;
                if (curr_knot_value > knot_interval[1]) {break;}
                if (std::abs(curr_knot_value - knot_interval[1]) < toll2*10) knot_interval[1] = curr_knot_value;
                double knot_value_in_surface_parameter = Start_U + (curr_knot_value-knot_interval[0]) * scale_factor;
                rSpansU.push_back(knot_value_in_surface_parameter);
                
            }
            
        } else if (std::abs(physical_coord_lower_left[1]-physical_coord_lower_right[1]) > toll) {// front or back face
            // direction_lower_segment = 1;
            knot_interval[0] = physical_coord_lower_left[1]-toll2; 
            knot_interval[1] = physical_coord_lower_right[1]+toll2;
            std::sort(knot_interval.begin(), knot_interval.end());
            // Compare with volume_spans_v
            for (int i = 0; i < volume_spans_v.size(); i++) {
                double curr_knot_value = volume_spans_v[i];
                if (curr_knot_value < knot_interval[0]) {continue;}
                if (std::abs(curr_knot_value - knot_interval[0]) < toll2*10) knot_interval[0] = curr_knot_value;
                if (curr_knot_value > knot_interval[1]) {break;}
                if (std::abs(curr_knot_value - knot_interval[1]) < toll2*10) knot_interval[1] = curr_knot_value;
                double knot_value_in_surface_parameter = Start_U + (curr_knot_value-knot_interval[0]) * scale_factor;
                rSpansU.push_back(knot_value_in_surface_parameter);
                
            }
        } else if (std::abs(physical_coord_lower_left[2]-physical_coord_lower_right[2]) > toll) {// bottom or top face
            // direction_lower_segment = 2;
            knot_interval[0] = physical_coord_lower_left[2]-toll2; 
            knot_interval[1] = physical_coord_lower_right[2]+toll2;
            std::sort(knot_interval.begin(), knot_interval.end());
            // Compare with volume_spans_w
            for (int i = 0; i < volume_spans_w.size(); i++) {
                double curr_knot_value = volume_spans_w[i];
                if (curr_knot_value < knot_interval[0]) {continue;}
                if (std::abs(curr_knot_value - knot_interval[0]) < toll2*10) knot_interval[0] = curr_knot_value;
                if (curr_knot_value > knot_interval[1]) {break;}
                if (std::abs(curr_knot_value - knot_interval[1]) < toll2*10) knot_interval[1] = curr_knot_value;
                double knot_value_in_surface_parameter = Start_U + (curr_knot_value-knot_interval[0]) * scale_factor;
                rSpansU.push_back(knot_value_in_surface_parameter);
                
            }
        } else {
            KRATOS_ERROR << "FACE NOT PARALLEL TO ANY VOLUME KNOT VECTOR";
        }
        // //------------------------------------------------------------------
        

       
        
        // LEFT SEGMENT
        // Compute constant coordinate -> understand where the face is in volume parameter space

        physical_length_segment = norm_2(physical_coord_lower_left-physical_coord_upper_left);
        parameter_length_segment = norm_2(local_coord_lower_left-local_coord_upper_left);
        scale_factor = parameter_length_segment/physical_length_segment;
        // int direction_left_segment;
        std::vector<double> intersection_knot_values_in_volume_left ;
        if (std::abs(physical_coord_lower_left[0]-physical_coord_upper_left[0]) > toll) {// left or right face
            // direction_left_segment = 0; // compare to knot_spans_vector_u
            knot_interval[0]  = physical_coord_lower_left[0]-toll2; 
            knot_interval[1] = physical_coord_upper_left[0]+toll2;
            std::sort(knot_interval.begin(), knot_interval.end());
            // Compare with volume_spans_u
            for (int i = 0; i < volume_spans_u.size(); i++) {
                double curr_knot_value = volume_spans_u[i];
                if (curr_knot_value < knot_interval[0]) {continue;}
                if (std::abs(curr_knot_value - knot_interval[0]) < toll2*10) knot_interval[0] = curr_knot_value;
                if (curr_knot_value > knot_interval[1]) {break;}
                if (std::abs(curr_knot_value - knot_interval[1]) < toll2*10) knot_interval[1] = curr_knot_value;
                double knot_value_in_surface_parameter = Start_V + (curr_knot_value-knot_interval[0]) * scale_factor;
                rSpansV.push_back(knot_value_in_surface_parameter);
                
            }
            
        } else if (std::abs(physical_coord_lower_left[1]-physical_coord_upper_left[1]) > toll) {// front or back face
            // direction_left_segment = 1;
            knot_interval[0] = physical_coord_lower_left[1]-toll2; 
            knot_interval[1] = physical_coord_upper_left[1]+toll2;
            std::sort(knot_interval.begin(), knot_interval.end());
            // Compare with volume_spans_v
            for (int i = 0; i < volume_spans_v.size(); i++) {
                double curr_knot_value = volume_spans_v[i];
                if (curr_knot_value < knot_interval[0]) {continue;}
                if (std::abs(curr_knot_value - knot_interval[0]) < toll2*10) knot_interval[0] = curr_knot_value;
                if (curr_knot_value > knot_interval[1]) {break;}
                if (std::abs(curr_knot_value - knot_interval[1]) < toll2*10) knot_interval[1] = curr_knot_value;
                double knot_value_in_surface_parameter = Start_V + (curr_knot_value-knot_interval[0]) * scale_factor;
                rSpansV.push_back(knot_value_in_surface_parameter);
                
            }
        } else if (std::abs(physical_coord_lower_left[2]-physical_coord_upper_left[2]) > toll) {// bottom or top face
            // direction_left_segment = 2;
            knot_interval[0] = physical_coord_lower_left[2]-toll2; 
            knot_interval[1] = physical_coord_upper_left[2]+toll2;
            std::sort(knot_interval.begin(), knot_interval.end());
            // Compare with volume_spans_w
            for (int i = 0; i < volume_spans_w.size(); i++) {
                double curr_knot_value = volume_spans_w[i];
                if (curr_knot_value < knot_interval[0]) {continue;}
                if (std::abs(curr_knot_value - knot_interval[0]) < toll2*10) knot_interval[0] = curr_knot_value;
                if (curr_knot_value > knot_interval[1]) {break;}
                if (std::abs(curr_knot_value - knot_interval[1]) < toll2*10) knot_interval[1] = curr_knot_value;
                double knot_value_in_surface_parameter = Start_V + (curr_knot_value-knot_interval[0]) * scale_factor;
                rSpansV.push_back(knot_value_in_surface_parameter);
                
            }
        } else {
            KRATOS_ERROR << "FACE NOT PARALLEL TO ANY VOLUME KNOT VECTOR";
        }

        
        if (rSpansU.size()>2) KRATOS_WATCH(rSpansU)
        if (rSpansV.size()>2) KRATOS_WATCH(rSpansV)
        KRATOS_ERROR_IF(rSpansU.size()<2) << "WRONG NUMBER OF INTERSECTION (<2) FOUND" << std::endl;
        KRATOS_ERROR_IF(rSpansV.size()<2) << "WRONG NUMBER OF INTERSECTION (<2) FOUND" << std::endl;

        // CurveAxisIntersection<CurveNodeType>::ComputeAxisIntersection(
        //     rSpans,
        //     *(mpNurbsCurve.get()), Start, End,
        //     surface_spans_u, surface_spans_v,
        //     1e-6);
    }

    ///@}
    ///@name Kratos Geometry Families
    ///@{

    GeometryData::KratosGeometryFamily GetGeometryFamily() const override
    {
        return GeometryData::KratosGeometryFamily::Kratos_Nurbs;
    }

    GeometryData::KratosGeometryType GetGeometryType() const override
    {
        return GeometryData::KratosGeometryType::Kratos_Nurbs_Surface_On_Volume;
    }

    ///@}
    ///@name Information
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "2 dimensional surface in 3D nurbs volume.";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "2 dimensional surface in 3D nurbs volume.";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    void SetIsExitingDirectionSBM(int isExitingDirection) 
    {
        isExitingDirectionSBM = isExitingDirection;
    }

    void SetNormalSBM(Vector Normal) 
    {
        msNormal = Normal;
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

    int isExitingDirectionSBM = -1;
    Vector msNormal;

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
        rSerializer.load("pNurbsVolume", mpNurbsVolume);
        rSerializer.load("pSurface", mpSurface);
    }

    ///@}

}; // class NurbsSurfaceOnVolumeGeometry

template<int TWorkingSpaceDimension, class TVolumeContainerPointType>
const GeometryData NurbsSurfaceOnVolumeGeometry<TWorkingSpaceDimension, TVolumeContainerPointType>::msGeometryData(
    &msGeometryDimension,
    GeometryData::IntegrationMethod::GI_GAUSS_1,
    {}, {}, {});

template<int TWorkingSpaceDimension, class TVolumeContainerPointType>
const GeometryDimension NurbsSurfaceOnVolumeGeometry<TWorkingSpaceDimension, TVolumeContainerPointType>::msGeometryDimension(TWorkingSpaceDimension, 3);

} // namespace Kratos

#endif // KRATOS_NURBS_SURFACE_ON_VOLUME_INCLUDED_H defined
