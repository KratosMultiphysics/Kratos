// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Navaneeth K Narayanan
//

// System includes

#if !defined(KRATOS_VOLUME_UNDER_PLANE_UTILITY_H_INCLUDED)
#define KRATOS_VOLUME_UNDER_PLANE_UTILITY_H_INCLUDED

/* System includes */

/* External includes */
#include "math.h"
#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/math_utils.h"
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "utilities/divide_geometry.h"
#include "utilities/divide_triangle_2d_3.h"
#include "geometries/point.h"
#include "geometries/triangle_3d_3.h"
#include "includes/condition.h"
#include "includes/node.h"
#include "includes/deprecated_variables.h"
namespace Kratos
{
/**
 * auxiliary functions for calculating volume of a 3D structure
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) VolumeCalculationUnderPlaneUtility
{

  public:
    /**
     * Type Definitions
     */

    KRATOS_CLASS_POINTER_DEFINITION(VolumeCalculationUnderPlaneUtility);

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    typedef Geometry<Node<3>> GeometryType;
    typedef GeometryType::Pointer GeometryPointerType;
    typedef GeometryType::CoordinatesArrayType CoordinatesArrayType;
    typedef GeometryData::IntegrationMethod IntegrationMethodType;
    typedef IntegrationPoint<3> IntegrationPointType;
    typedef std::vector<IntegrationPointType> IntegrationPointsArrayType;
    typedef DivideGeometry::IndexedPointGeometryType IndexedPointGeometryType;

    /**
     * Life Cycle
     */

    /**
     * Constructor.
     */
    VolumeCalculationUnderPlaneUtility(Vector Centre, double Radius, Vector Normal) : mCentre(Centre), mRadius(Radius)
    {

        double norm = norm_2(Normal);
        if (norm < std::numeric_limits<double>::epsilon())
        {
            Normal *= 0;
            KRATOS_WARNING("The plane normal is set to zero");
        }
        else
            mNormal = Normal / norm;

        mRefCentre = mCentre;
        mVolume = 0.0;
        mIntersectedArea = 0.0;
    }

    //Default constructor

    VolumeCalculationUnderPlaneUtility()

    {
        mCentre = ZeroVector(3);
        mRefCentre = ZeroVector(3);
        mRadius = 0.0;

        mNormal = ZeroVector(3);
        mNormal[2] = 1.0;
        mVolume = 0.0;
        mIntersectedArea = 0.0;
    }

    /**
     * Destructor.
     */
    virtual ~VolumeCalculationUnderPlaneUtility() {}

    /**
     * Operations
     */

    /**
       * This function calculates the volume of a  surface under a plane
       * @param rModelPart the model part
       */
    // VM !!!!!!!!!

    double CalculateVolume(ModelPart &rModelPart)

    {
        KRATOS_TRY;

        std::cout << "Calculating volume enclosed by plane" << std::endl;
        const double dist_limit = pow(10, std::numeric_limits<double>::digits10);
        double volume = 0.0;
        double int_distance_dot_n = 0.0;
        double int_area_dot_n_plane = 0.0;
        double intersected_area = 0.0;
        double node_h_dist;
        double node_v_distance;
        array_1d<double, 3> node_vector;

        IntegrationMethodType IntegrationMethod = GeometryData::GI_GAUSS_1;

        // Check if the first node has DISPLACEMENT variable
        bool HasDisplacement = rModelPart.NodesBegin()->SolutionStepsDataHas(DISPLACEMENT);

        if (!(HasDisplacement))
            KRATOS_WARNING("The model doesn't have DISPLACEMENT as variable. The predicted displacement will be set to zero");

        for (ModelPart::NodesContainerType::iterator i_node = rModelPart.NodesBegin(); i_node != rModelPart.NodesEnd(); i_node++)
        {

            i_node->FastGetSolutionStepValue(DISTANCE) = 1.0 + dist_limit; // Initialise with max_distance
            i_node->FastGetSolutionStepValue(NORMAL) = ZeroVector(3);      // initialize nodal normal to zero
        }

        // Calculation of vertical distance of nodes of the model part from the cutting plane
        if (mRadius < std::numeric_limits<double>::epsilon())
        {

            for (ModelPart::NodesContainerType::iterator i_node = rModelPart.NodesBegin(); i_node != rModelPart.NodesEnd(); i_node++)
            {

                node_v_distance = CalculateDistanceFromPlane(i_node->Coordinates());
                i_node->FastGetSolutionStepValue(DISTANCE) = node_v_distance;
            }
        }

        else
        {

            for (ModelPart::NodesContainerType::iterator i_node = rModelPart.NodesBegin(); i_node != rModelPart.NodesEnd(); i_node++)
            {

                node_v_distance = CalculateDistanceFromPlane(i_node->Coordinates());
                node_vector = i_node->Coordinates() - mCentre;
                node_h_dist = norm_2(node_vector - node_v_distance * mNormal);

                if (node_h_dist < mRadius)
                    i_node->FastGetSolutionStepValue(DISTANCE) = node_v_distance;
            }
        }

        for (ModelPart::ConditionsContainerType::iterator cond = rModelPart.ConditionsBegin(); cond != rModelPart.ConditionsEnd(); cond++)
        {

            array_1d<double, 3> distances_vector;

            bool is_split = false;
            bool is_negative = false;

            GeometryType &geom = cond->GetGeometry();

            if (geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle3D3)
            {

                for (IndexType i = 0; i < geom.size(); ++i)
                    distances_vector(i) = geom[i].FastGetSolutionStepValue(DISTANCE);

                cond->SetValue(ELEMENTAL_DISTANCES, distances_vector);

                Vector &r_elemental_distances = cond->GetValue(ELEMENTAL_DISTANCES);

                DivideTriangle2D3 triangle_splitter(geom, r_elemental_distances);

                VolumeCalculationUnderPlaneUtility::IsNegativeOrSplit(geom, is_split, is_negative);

                if (is_split)
                {
                    // Call the divide geometry method

                    triangle_splitter.GenerateDivision();

                    for (IndexType i = 0; i < triangle_splitter.mNegativeSubdivisions.size(); i++)
                    {

                        IndexedPointGeometryType indexed_subgeom = *(triangle_splitter.mNegativeSubdivisions[i]);

                        Node<3>::Pointer node_i = Node<3>::Pointer(new Node<3>(indexed_subgeom[0].Id(), indexed_subgeom[0].X(), indexed_subgeom[0].Y(), indexed_subgeom[0].Z()));

                        Node<3>::Pointer node_j = Node<3>::Pointer(new Node<3>(indexed_subgeom[1].Id(), indexed_subgeom[1].X(), indexed_subgeom[1].Y(), indexed_subgeom[1].Z()));
                        Node<3>::Pointer node_k = Node<3>::Pointer(new Node<3>(indexed_subgeom[2].Id(), indexed_subgeom[2].X(), indexed_subgeom[2].Y(), indexed_subgeom[2].Z()));

                        Triangle3D3<Node<3>>::Pointer subgeom = Triangle3D3<Node<3>>::Pointer(new Triangle3D3<Node<3>>(node_i, node_j, node_k));

                        VolumeCalculationUnderPlaneUtility::CalculateIntDistanceDotN(*subgeom, IntegrationMethod, int_distance_dot_n);
                        VolumeCalculationUnderPlaneUtility::CalculateIntAreaDotNplane(*subgeom, IntegrationMethod, int_area_dot_n_plane);
                        VolumeCalculationUnderPlaneUtility::CalculateAndAssignNodalNormal(*subgeom, geom, IntegrationMethod);
                        /* if (HasDisplacement)
                            VolumeCalculationUnderPlaneUtility::CalculateDisplacementDotN(*subgeom, geom, IntegrationMethod, dV); //TODO: Use the nodal normals to predict */
                    }
                }
                else
                {

                    if (is_negative)
                    {

                        VolumeCalculationUnderPlaneUtility::CalculateIntDistanceDotN(geom, IntegrationMethod, int_distance_dot_n);
                        VolumeCalculationUnderPlaneUtility::CalculateIntAreaDotNplane(geom, IntegrationMethod, int_area_dot_n_plane);
                        VolumeCalculationUnderPlaneUtility::CalculateAndAssignNodalNormal(geom, geom, IntegrationMethod);
                        /* if (HasDisplacement)
                            VolumeCalculationUnderPlaneUtility::CalculateDisplacementDotN(geom, geom, IntegrationMethod, dV); //TODO: Use the nodal normals to predict */
                    }
                }
            }
        }

        volume = int_distance_dot_n;
        intersected_area = -int_area_dot_n_plane;

        mVolume = volume;

        mIntersectedArea = intersected_area;

        return volume;

        KRATOS_CATCH(" ");
    }

    void CalculateVolumeForConditions(WeakPointerVector<Condition> &rConditonWeakPointersVector, WeakPointerVector<Node<3>> &rNodeWeakPointersVector)

    {
        KRATOS_TRY;

        std::cout << "Calculating volume enclosed by plane" << std::endl;
        const double dist_limit = pow(10, std::numeric_limits<double>::digits10);
        double volume = 0.0;
        double int_distance_dot_n = 0.0;
        double int_area_dot_n_plane = 0.0;
        double intersected_area = 0.0;
        double node_h_dist;
        double node_v_distance;
        array_1d<double, 3> node_vector;

        IntegrationMethodType IntegrationMethod = GeometryData::GI_GAUSS_1;

        for (WeakPointerVector<Node<3>>::iterator i_node = rNodeWeakPointersVector.begin(); i_node != rNodeWeakPointersVector.end(); i_node++)
        {

            i_node->FastGetSolutionStepValue(DISTANCE) = 1.0 + dist_limit; // Initialise with max_distance
            i_node->FastGetSolutionStepValue(NORMAL) = ZeroVector(3);
        }

        // Calculation of vertical distance of nodes of the model part from the cutting plane
        if (mRadius < std::numeric_limits<double>::epsilon())
        {

            for (WeakPointerVector<Node<3>>::iterator i_node = rNodeWeakPointersVector.begin(); i_node != rNodeWeakPointersVector.end(); i_node++)
            {

                node_v_distance = CalculateDistanceFromPlane(i_node->Coordinates());
                i_node->FastGetSolutionStepValue(DISTANCE) = node_v_distance;
            }
        }

        else
        {

            for (WeakPointerVector<Node<3>>::iterator i_node = rNodeWeakPointersVector.begin(); i_node != rNodeWeakPointersVector.end(); i_node++)
            {

                node_v_distance = VolumeCalculationUnderPlaneUtility::CalculateDistanceFromPlane(i_node->Coordinates());
                node_vector = i_node->Coordinates() - mCentre;
                node_h_dist = norm_2(node_vector - node_v_distance * mNormal);

                if (node_h_dist < mRadius)
                    i_node->FastGetSolutionStepValue(DISTANCE) = node_v_distance;
            }
        }

        for (WeakPointerVector<Condition>::iterator cond = rConditonWeakPointersVector.begin(); cond != rConditonWeakPointersVector.end(); cond++)
        {

            array_1d<double, 3> distances_vector;

            bool is_split = false;
            bool is_negative = false;

            GeometryType &geom = cond->GetGeometry();

            if (geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle3D3)
            {

                for (IndexType i = 0; i < geom.size(); ++i)
                    distances_vector(i) = geom[i].FastGetSolutionStepValue(DISTANCE);

                cond->SetValue(ELEMENTAL_DISTANCES, distances_vector);

                Vector &r_elemental_distances = cond->GetValue(ELEMENTAL_DISTANCES);

                DivideTriangle2D3 triangle_splitter(geom, r_elemental_distances);

                VolumeCalculationUnderPlaneUtility::IsNegativeOrSplit(geom, is_split, is_negative);

                if (is_split)
                {
                    // Call the divide geometry method

                    triangle_splitter.GenerateDivision();

                    for (IndexType i = 0; i < triangle_splitter.mNegativeSubdivisions.size(); i++)
                    {

                        IndexedPointGeometryType indexed_subgeom = *(triangle_splitter.mNegativeSubdivisions[i]);

                        Node<3>::Pointer node_i = Node<3>::Pointer(new Node<3>(indexed_subgeom[0].Id(), indexed_subgeom[0].X(), indexed_subgeom[0].Y(), indexed_subgeom[0].Z()));

                        Node<3>::Pointer node_j = Node<3>::Pointer(new Node<3>(indexed_subgeom[1].Id(), indexed_subgeom[1].X(), indexed_subgeom[1].Y(), indexed_subgeom[1].Z()));
                        Node<3>::Pointer node_k = Node<3>::Pointer(new Node<3>(indexed_subgeom[2].Id(), indexed_subgeom[2].X(), indexed_subgeom[2].Y(), indexed_subgeom[2].Z()));

                        Triangle3D3<Node<3>>::Pointer subgeom = Triangle3D3<Node<3>>::Pointer(new Triangle3D3<Node<3>>(node_i, node_j, node_k));

                        VolumeCalculationUnderPlaneUtility::CalculateIntDistanceDotN(*subgeom, IntegrationMethod, int_distance_dot_n);
                        VolumeCalculationUnderPlaneUtility::CalculateIntAreaDotNplane(*subgeom, IntegrationMethod, int_area_dot_n_plane);
                        VolumeCalculationUnderPlaneUtility::CalculateAndAssignNodalNormal(*subgeom, geom, IntegrationMethod);
                    }
                }
                else
                {

                    if (is_negative)
                    {

                        VolumeCalculationUnderPlaneUtility::CalculateIntDistanceDotN(geom, IntegrationMethod, int_distance_dot_n);
                        VolumeCalculationUnderPlaneUtility::CalculateIntAreaDotNplane(geom, IntegrationMethod, int_area_dot_n_plane);
                        VolumeCalculationUnderPlaneUtility::CalculateAndAssignNodalNormal(geom, geom, IntegrationMethod);
                    }
                }
            }
        }

        volume = int_distance_dot_n;
        intersected_area = -int_area_dot_n_plane;

        mVolume = volume;

        mIntersectedArea = intersected_area;

        KRATOS_CATCH(" ");
    }

    void CalculateIntDistanceDotN(GeometryType &rGeom, const IntegrationMethodType &IntegrationMethod, double &rIntDistanceDotN)

    {
        KRATOS_TRY;
        const unsigned int n_int_pts = rGeom.IntegrationPointsNumber(IntegrationMethod);

        double v_distance = 0.0;
        array_1d<double, 3> distance_vector;
        CoordinatesArrayType proj_global_coords;
        double distance_dot_n;

        // Get the  Gauss points coordinates
        IntegrationPointsArrayType gauss_pts;
        gauss_pts = rGeom.IntegrationPoints(IntegrationMethod);

        Vector jacobians_values;
        rGeom.DeterminantOfJacobian(jacobians_values, IntegrationMethod);

        // Degubbing navaneeth

        array_1d<double, 3> area_normal;
        double area;
        // Compute the x dot n at every gauss points
        for (IndexType i_gauss = 0; i_gauss < n_int_pts; ++i_gauss)
        {

            area_normal = rGeom.AreaNormal(gauss_pts[i_gauss].Coordinates());
            area = norm_2(area_normal);
            if (area > std::numeric_limits<double>::epsilon())
                area_normal /= area;
            else
            {

                area_normal *= 0;
            }

            // Compute the global coordinates of the Gauss pt.
            CoordinatesArrayType global_coords = ZeroVector(3);

            rGeom.GlobalCoordinates(global_coords, gauss_pts[i_gauss].Coordinates());

            v_distance = CalculateDistanceFromPlane(global_coords);

            noalias(distance_vector) = v_distance * mNormal;

            distance_dot_n = MathUtils<double>::Dot(distance_vector, area_normal);

            rIntDistanceDotN += distance_dot_n * gauss_pts[i_gauss].Weight() * jacobians_values[i_gauss];
        }
        KRATOS_CATCH("")
    }

    void CalculateIntAreaDotNplane(GeometryType &rGeom, const IntegrationMethodType &IntegrationMethod, double &rIntAreaDotNplane)

    {
        KRATOS_TRY;
        const unsigned int n_int_pts = rGeom.IntegrationPointsNumber(IntegrationMethod);

        double area_dot_n_normal;

        // Get the  Gauss points coordinates
        IntegrationPointsArrayType gauss_pts;
        gauss_pts = rGeom.IntegrationPoints(IntegrationMethod);

        Vector jacobians_values;
        rGeom.DeterminantOfJacobian(jacobians_values, IntegrationMethod);

        // Compute the area dot plane normal at every gauss points
        for (IndexType i_gauss = 0; i_gauss < n_int_pts; ++i_gauss)
        {

            array_1d<double, 3> area_normal = rGeom.AreaNormal(gauss_pts[i_gauss].Coordinates());

            area_dot_n_normal = MathUtils<double>::Dot(area_normal, mNormal);

            rIntAreaDotNplane += area_dot_n_normal;
        }
        KRATOS_CATCH("")
    }

    /* void CalculateDisplacementDotN(GeometryType &rSubGeom, GeometryType &rGeom, const IntegrationMethodType &IntegrationMethod, double &rdV)
    {
        KRATOS_TRY;

        const unsigned int n_int_pts = rSubGeom.IntegrationPointsNumber(IntegrationMethod);
        double displacement_dot_n = 0.0;
        IntegrationPointsArrayType gauss_pts;

        // Gauss points from subgeom
        gauss_pts = rSubGeom.IntegrationPoints(IntegrationMethod);
        Vector jacobians_values;
        rSubGeom.DeterminantOfJacobian(jacobians_values, IntegrationMethod);
        array_1d<double, 3> area_normal;
        double area;

        for (IndexType i_gauss = 0; i_gauss < n_int_pts; ++i_gauss)
        {
            CoordinatesArrayType global_coords = ZeroVector(3);
            CoordinatesArrayType local_coords = ZeroVector(3);

            // Degubbing navaneeth
            area_normal = rSubGeom.AreaNormal(gauss_pts[i_gauss].Coordinates());
            area = norm_2(area_normal);

            if (area > std::numeric_limits<double>::epsilon())
                area_normal /= area;
            else
            {

                area_normal *= 0;
            }

            // Interpolate the displacement on the gauss points of subgeometry
            rSubGeom.GlobalCoordinates(global_coords, gauss_pts[i_gauss].Coordinates());
            rGeom.PointLocalCoordinates(local_coords, global_coords);
            Vector N(rGeom.size());
            array_1d<double, 3> displacement = ZeroVector(3);
            rGeom.ShapeFunctionsValues(N, local_coords);
            for (IndexType i_node = 0; i_node < rGeom.size(); i_node++)
            {

                displacement += N[i_node] * (rGeom[i_node].GetSolutionStepValue(DISPLACEMENT, 0) - rGeom[i_node].GetSolutionStepValue(DISPLACEMENT, 1));
            }

            // Integration displacement dot normal
            displacement_dot_n = MathUtils<double>::Dot(displacement, area_normal);
            rdV += displacement_dot_n * gauss_pts[i_gauss].Weight() * jacobians_values[i_gauss];
        }
        KRATOS_CATCH("")
    } */

    void CalculateAndAssignNodalNormal(GeometryType &rSubGeom, GeometryType &rGeom, const IntegrationMethodType &IntegrationMethod)
    {
        KRATOS_TRY;

        const unsigned int n_int_pts = rSubGeom.IntegrationPointsNumber(IntegrationMethod);
        SizeType number_of_nodes = rGeom.size();
        IntegrationPointsArrayType gauss_pts;

        // Gauss points from subgeom
        gauss_pts = rSubGeom.IntegrationPoints(IntegrationMethod);
        Vector jacobians_values;
        rSubGeom.DeterminantOfJacobian(jacobians_values, IntegrationMethod);
        array_1d<double, 3> area_normal;
        double area;
        Vector N(number_of_nodes);
        CoordinatesArrayType global_coords = ZeroVector(3);
        CoordinatesArrayType local_coords = ZeroVector(3);

        for (IndexType i_node = 0; i_node < number_of_nodes; ++i_node)
        {

            for (IndexType i_gauss = 0; i_gauss < n_int_pts; ++i_gauss)
            {

                area_normal = rSubGeom.AreaNormal(gauss_pts[i_gauss].Coordinates());
                area = norm_2(area_normal);

                if (area > std::numeric_limits<double>::epsilon())
                    area_normal /= area;
                else
                {

                    area_normal *= 0;
                }
                rSubGeom.GlobalCoordinates(global_coords, gauss_pts[i_gauss].Coordinates());
                rGeom.PointLocalCoordinates(local_coords, global_coords);
                rGeom.ShapeFunctionsValues(N, local_coords);

                //#pragma omp atomic
                noalias(rGeom[i_node].FastGetSolutionStepValue(NORMAL)) += area_normal * N[i_node] * gauss_pts[i_gauss].Weight() * jacobians_values[i_gauss];
            }
        }

        KRATOS_CATCH("")
    }

    static void IsNegativeOrSplit(GeometryType &rGeom, bool &rIsSplit, bool &rIsNegative)
    {

        unsigned int n_pos = 0, n_neg = 0, n_zero = 0;
        double nodal_distance = 0.0;
        const double dist_limit = pow(10, std::numeric_limits<double>::digits10);

        for (IndexType i_node = 0; i_node < rGeom.size(); ++i_node)
        {
            nodal_distance = rGeom[i_node].FastGetSolutionStepValue(DISTANCE);

            if (nodal_distance < 0.0 - std::numeric_limits<double>::epsilon())
            {
                n_neg++;
            }
            else if ((nodal_distance > 0.0 + std::numeric_limits<double>::epsilon()) && (nodal_distance < dist_limit))
            {
                n_pos++;
            }

            else if (nodal_distance > dist_limit)
            {
                n_zero++;
            }
            else
            {
                n_zero++;
                n_neg++;
            }
        }

        if ((n_pos > 0) && (n_neg > 0))
            rIsSplit = true;

        else
            rIsSplit = false;

        if (n_zero > 0)
            rIsSplit = false;

        if (n_neg == rGeom.size())
            rIsNegative = true;

        else
            rIsNegative = false;
    }

    double CalculateDistanceFromPlane(const CoordinatesArrayType &rPoint)
    {

        return inner_prod(mNormal, rPoint) - inner_prod(mNormal, mCentre);
    }
    Vector GetCentre() const
    {
        return mCentre;
    }

    double GetVolume() const
    {
        return mVolume;
    }

    double GetIntersectedArea() const
    {
        return mIntersectedArea;
    }

    double &GetPlaneRadius()
    {
        return this->mRadius;
    }

    Vector &GetPlaneNormal()
    {
        return this->mNormal;
    }

    Vector &GetPlaneCentre()
    {
        return this->mCentre;
    }

    void SetPlaneParameters(Vector rCentre, double rRadius, Vector rNormal)
    {

        mCentre = rCentre;
        mRadius = rRadius;
        mNormal = rNormal;
        mRefCentre = mCentre;
    }

    void UpdatePositionOfPlaneBasedOnTargetVolume(ModelPart &rModelPart, const double VolTarget, const double MaxRelRes = 1E-6, const unsigned int MaxIterations = 20)
    {
        KRATOS_TRY;

        double vol = CalculateVolume(rModelPart);
        double area = mIntersectedArea;
        double vol_inter = 0.0;
        array_1d<double, 3> displacement_vector;
        double movement;
        bool HasDisplacement = rModelPart.NodesBegin()->SolutionStepsDataHas(DISPLACEMENT);

        ///For checking the nodal normal implementation

        double displacement_value = 0.0;

        if (HasDisplacement)
        {
            for (ModelPart::NodesContainerType::iterator i_node = rModelPart.NodesBegin(); i_node != rModelPart.NodesEnd(); i_node++)
            {
                displacement_value += MathUtils<double>::Dot((i_node->GetSolutionStepValue(DISPLACEMENT, 0) - i_node->GetSolutionStepValue(DISPLACEMENT, 1)), i_node->FastGetSolutionStepValue(NORMAL));
            }
        }

        //Predictor for the plane's position

        if (area <= std::numeric_limits<double>::epsilon())
        {
            KRATOS_WARNING("Predicted displacement is Nan") << "The plane position is reinitialzed to :: " << mRefCentre[0] << ", " << mRefCentre[1] << ", " << mRefCentre[2] << std::endl;
            displacement_value = MathUtils<double>::Dot((mRefCentre - mCentre), mNormal);
        }

        else
        {
            //displacement_value /= area;
            displacement_value = 0.0;
        }

        displacement_vector = displacement_value * mNormal;
        UpdatePlaneCentre(displacement_vector);
        vol = CalculateVolume(rModelPart);
        area = mIntersectedArea;
        double vol_residual = VolTarget - vol;
        double vol_inter_residual = 0.0;
        double vol_rel_residual = vol_residual / VolTarget;
        unsigned int iteration_nr = 0;
        movement = displacement_value;

        std::cout << "Volume after predictor step :: " << vol << std::endl;
        std::cout << "Movement after predictor step  :: " << movement << std::endl;

        // Corrector

        /*   if (type == "NewtonRaphson")
        {
            while (fabs(vol_rel_residual) >= MaxRelRes && iteration_nr < MaxIterations)
            {
                if (area <= std::numeric_limits<double>::epsilon())
                    KRATOS_ERROR << "Intersected area is zero, Newton Raphson method failed. Try with the different initial position of the plane" << std::endl;

                displacement_value = vol_residual / area;
                displacement_vector = displacement_value * mNormal;
                UpdatePlaneCentre(displacement_vector);
                movement += displacement_value;

                vol = CalculateVolume(rModelPart);
                area = mIntersectedArea;
                vol_residual = (VolTarget - vol);
                vol_rel_residual = vol_residual / VolTarget;
                iteration_nr += 1;

                std::cout << "Iteration Nr. :: " << iteration_nr << std::endl;
                std::cout << "Relative Volume residual :: " << vol_rel_residual << std::endl;
                std::cout << "Volume :: " << vol << std::endl;
                std::cout << "Movement :: " << movement << std::endl;

                if (iteration_nr == MaxIterations)
                    KRATOS_WARNING("Max iterations reached. Newton Raphson method didn't converge");
            }
        }

        else if (type == "LeapFroggingNewton") */
        // {
        while (std::fabs(vol_rel_residual) >= MaxRelRes && std::fabs(vol_residual - vol_inter_residual) >= std::numeric_limits<double>::epsilon() && iteration_nr < MaxIterations)
        {
            if (area <= std::numeric_limits<double>::epsilon())
                KRATOS_ERROR << "Intersected area is zero, Leapfrogging Newton method failed. Try with the different initial position of the plane" << std::endl;

            displacement_value = vol_residual / area;
            displacement_vector = displacement_value * mNormal;
            UpdatePlaneCentre(displacement_vector);
            vol_inter = CalculateVolume(rModelPart);
            displacement_vector *= -1;
            UpdatePlaneCentre(displacement_vector);
            vol_inter_residual = (VolTarget - vol_inter);

            displacement_value = (vol_residual * vol_residual) / ((vol_residual - vol_inter_residual) * area);
            displacement_vector = displacement_value * mNormal;
            UpdatePlaneCentre(displacement_vector);
            movement += displacement_value;

            vol = CalculateVolume(rModelPart);
            area = mIntersectedArea;
            vol_residual = (VolTarget - vol);
            vol_rel_residual = vol_residual / VolTarget;

            iteration_nr += 1;

            std::cout << "Iteration Nr. :: " << iteration_nr << std::endl;
            std::cout << "Relative Volume residual :: " << vol_rel_residual << std::endl;
            std::cout << "Movement :: " << movement << std::endl;
            std::cout << "Centre :: " << mCentre[0] << ", " << mCentre[1] << ", " << mCentre[2] << std::endl;
            std::cout << "Volume :: " << vol << std::endl;

            if (iteration_nr == MaxIterations)
                KRATOS_WARNING("Max iterations reached. Leap Frogging Newton method didn't converge");
        }
        // }

        /*   else
            KRATOS_ERROR << "String specifying NewtonRaphson or LeapFroggingNewton is required as argument" << std::endl; */

        KRATOS_CATCH(" ");
    }

    void UpdatePositionOfPlaneBasedOnTargetVolumeForConditions(WeakPointerVector<Condition> &rConditionWeakPointersVector, WeakPointerVector<Node<3>> &rNodeWeakPointersVector, const double VolTarget, const double MaxRelRes = 1E-6, const unsigned int MaxIterations = 20)
    {
        KRATOS_TRY;

        CalculateVolumeForConditions(rConditionWeakPointersVector, rNodeWeakPointersVector);
        double vol = GetVolume();
        double area = mIntersectedArea;

        double vol_inter = 0.0;
        array_1d<double, 3> displacement_vector;
        double movement;

        ///For checking the nodal normal implementation

        double displacement_value = 0.0;
        bool HasDisplacement = rNodeWeakPointersVector.begin()->SolutionStepsDataHas(DISPLACEMENT);
        if (HasDisplacement)
        {
            for (WeakPointerVector<Node<3>>::iterator i_node = rNodeWeakPointersVector.begin(); i_node != rNodeWeakPointersVector.end(); i_node++)
            {

                displacement_value += MathUtils<double>::Dot((i_node->GetSolutionStepValue(DISPLACEMENT, 0) - i_node->GetSolutionStepValue(DISPLACEMENT, 1)), i_node->FastGetSolutionStepValue(NORMAL));
            }
        }

        //Predictor for the plane's position
        if (area < std::numeric_limits<double>::epsilon())
        {

            KRATOS_WARNING("Predicted displacement is Nan") << "The plane position is reinitialzed to :: " << mRefCentre[0] << ", " << mRefCentre[1] << ", " << mRefCentre[2] << std::endl;
            displacement_value = MathUtils<double>::Dot((mRefCentre - mCentre), mNormal);
        }

        else
        {

            displacement_value /= area;
        }

        displacement_vector = displacement_value * mNormal;
        UpdatePlaneCentre(displacement_vector);
        CalculateVolumeForConditions(rConditionWeakPointersVector, rNodeWeakPointersVector);
        vol = GetVolume();
        area = mIntersectedArea;
        double vol_residual = VolTarget - vol;
        double vol_inter_residual = 0.0;
        double vol_rel_residual = vol_residual / VolTarget;
        unsigned int iteration_nr = 0;
        movement = displacement_value;

        //std::cout << "Volume after predictor step :: " << vol << std::endl;
        std::cout << "Movement after predictor step (ref) :: " << movement << std::endl;

        // Corrector

        while (std::fabs(vol_rel_residual) >= MaxRelRes && std::fabs(vol_residual - vol_inter_residual) >= std::numeric_limits<double>::epsilon() && iteration_nr < MaxIterations)
        {

            if (area <= std::numeric_limits<double>::epsilon())
                KRATOS_ERROR << "Intersected area is zero, Leapfrogging Newton method failed. Try with the different initial position of the plane" << std::endl;

            displacement_value = vol_residual / area;
            displacement_vector = displacement_value * mNormal;
            UpdatePlaneCentre(displacement_vector);
            CalculateVolumeForConditions(rConditionWeakPointersVector, rNodeWeakPointersVector);
            vol_inter = GetVolume();
            displacement_vector *= -1;
            UpdatePlaneCentre(displacement_vector);
            vol_inter_residual = (VolTarget - vol_inter);

            displacement_value = (vol_residual * vol_residual) / ((vol_residual - vol_inter_residual) * area);
            displacement_vector = displacement_value * mNormal;
            UpdatePlaneCentre(displacement_vector);
            movement += displacement_value;

            CalculateVolumeForConditions(rConditionWeakPointersVector, rNodeWeakPointersVector);
            vol = GetVolume();
            area = mIntersectedArea;
            vol_residual = (VolTarget - vol);
            vol_rel_residual = vol_residual / VolTarget;

            iteration_nr += 1;

            std::cout << "Iteration Nr. :: " << iteration_nr << std::endl;
            std::cout << "Relative Volume residual :: " << vol_rel_residual << std::endl;
            std::cout << "Movement :: " << movement << std::endl;
            std::cout << "Centre :: " << mCentre[0] << ", " << mCentre[1] << ", " << mCentre[2] << std::endl;
            std::cout << "Volume :: " << vol << std::endl;

            if (iteration_nr == MaxIterations)
                KRATOS_WARNING("Max iterations reached. Leap Frogging Newton method didn't converge");
        }
        // }

        /*   else
            KRATOS_ERROR << "String specifying NewtonRaphson or LeapFroggingNewton is required as argument" << std::endl; */

        KRATOS_CATCH(" ");
    }

    void UpdatePlaneCentre(array_1d<double, 3> &rDisplacement_vector)
    {

        mCentre += rDisplacement_vector;
    }

  private:
    Vector mCentre;
    double mRadius;
    Vector mNormal;

    Vector mRefCentre;

    double mVolume;
    double mIntersectedArea;

}; // namespace Kratos
} /* namespace Kratos.*/

#endif /* KRATOS_VOLUME_UTILITY_H_INCLUDED  defined */
