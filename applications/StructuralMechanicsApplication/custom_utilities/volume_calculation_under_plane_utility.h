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
#include "geometries/plane.h"
#include "utilities/divide_triangle_2d_3.h"
#include "geometries/point.h"
#include "geometries/triangle_3d_3.h"
#include "includes/element.h"
#include "includes/node.h"
namespace Kratos
{
/**
 * auxiliary functions for calculating volume of a 3D structure
 */

const double tol = 1E-12;

class VolumeCalculationUnderPlaneUtility
{
  public:
    /**
     * Type Definitions
     */

    KRATOS_CLASS_POINTER_DEFINITION(VolumeCalculationUnderPlaneUtility);

    typedef Geometry<Node<3>> GeometryType;
    typedef GeometryType::Pointer GeometryPointerType;
    typedef GeometryType::CoordinatesArrayType CoordinatesArrayType;
    typedef GeometryData::IntegrationMethod IntegrationMethodType;
    typedef IntegrationPoint<3> IntegrationPointType;
    typedef std::vector<IntegrationPointType> IntegrationPointsArrayType;
    typedef boost::array<IntegrationPointsArrayType, GeometryData::NumberOfIntegrationMethods> IntegrationPointsContainerType;
    typedef DivideGeometry::IndexedPointGeometryType IndexedPointGeometryType;

    /**
     * Life Cycle
     */

    /**
     * Constructor.
     */
    VolumeCalculationUnderPlaneUtility(ModelPart &r_plane_model_part) : m_plane_model_part(r_plane_model_part)
    {

        m_volume = 0.0;
        m_intersected_area = 0.0;
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
       * @param mr_model_part the model part
       */
    // VM !!!!!!!!!

    double CalculateVolume(ModelPart &r_model_part)

    {
        KRATOS_TRY;

        std::cout << "Calculating volume enclosed by plane" << std::endl;
        double volume = 0.0;
        double int_distance_dot_n = 0.0;
        double int_area_dot_n_plane = 0.0;
        double intersected_area = 0.0;

        IntegrationMethodType IntegrationMethod = GeometryData::GI_GAUSS_1;

        // Calculation of cutting plane parameters
        ModelPart::ConditionsContainerType::iterator icond = m_plane_model_part.ConditionsBegin();

        CoordinatesArrayType p0 = icond->GetGeometry()[0].Coordinates();
        CoordinatesArrayType p1 = icond->GetGeometry()[1].Coordinates();
        CoordinatesArrayType p2 = icond->GetGeometry()[2].Coordinates();

        Plane cutting_plane(p0, p1, p2);

        // Calculation of distance of nodes of the model part from the cutting plane

        for (ModelPart::NodesContainerType::iterator inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); inode++)
        {

            double distance = cutting_plane.DistanceTo(inode->Coordinates());

            inode->FastGetSolutionStepValue(DISTANCE) = distance;
        }

        for (ModelPart::ElementsContainerType::iterator ielem = r_model_part.ElementsBegin(); ielem != r_model_part.ElementsEnd(); ielem++)
        {

            array_1d<double, 3> distances_vector;

            bool is_split = false;
            bool is_negative = false;

            GeometryType &geom = ielem->GetGeometry();

            for (unsigned int i = 0; i < geom.size(); ++i)
                distances_vector(i) = geom[i].FastGetSolutionStepValue(DISTANCE);

            ielem->SetValue(ELEMENTAL_DISTANCES, distances_vector);

            Vector &r_elemental_distances = ielem->GetValue(ELEMENTAL_DISTANCES);

            DivideTriangle2D3 triangle_splitter(geom, r_elemental_distances);

            IsNegativeOrSplit(geom, is_split, is_negative);

            if (is_split)
            {
                IsNegativeOrSplit(geom, is_split, is_negative);
                // Call the divide geometry method
                triangle_splitter.GenerateDivision();

                for (unsigned int i = 0; i < triangle_splitter.mNegativeSubdivisions.size(); i++)
                {

                    IndexedPointGeometryType indexed_subgeom = *(triangle_splitter.mNegativeSubdivisions[i]);

                    Node<3>::Pointer node_i = Node<3>::Pointer(new Node<3>(indexed_subgeom[0].Id(), indexed_subgeom[0].X(), indexed_subgeom[0].Y(), indexed_subgeom[0].Z()));

                    Node<3>::Pointer node_j = Node<3>::Pointer(new Node<3>(indexed_subgeom[1].Id(), indexed_subgeom[1].X(), indexed_subgeom[1].Y(), indexed_subgeom[1].Z()));
                    Node<3>::Pointer node_k = Node<3>::Pointer(new Node<3>(indexed_subgeom[2].Id(), indexed_subgeom[2].X(), indexed_subgeom[2].Y(), indexed_subgeom[2].Z()));

                    Triangle3D3<Node<3>>::Pointer subgeom = Triangle3D3<Node<3>>::Pointer(new Triangle3D3<Node<3>>(node_i, node_j, node_k));

                    CalculateIntDistanceDotN(*subgeom, cutting_plane, IntegrationMethod, int_distance_dot_n);
                    CalculateIntAreaDotNplane(*subgeom, cutting_plane, IntegrationMethod, int_area_dot_n_plane);
                }
            }
            else
            {

                if (is_negative)
                {

                    CalculateIntDistanceDotN(geom, cutting_plane, IntegrationMethod, int_distance_dot_n);
                    CalculateIntAreaDotNplane(geom, cutting_plane, IntegrationMethod, int_area_dot_n_plane);
                }
            }
        }

        volume = int_distance_dot_n;
        intersected_area = -int_area_dot_n_plane;

        m_volume = volume;
        m_intersected_area = intersected_area;
        return volume;

        KRATOS_CATCH(" ");
    }

    void CalculateIntDistanceDotN(GeometryType &r_geom, Plane &r_cutting_plane, const IntegrationMethodType IntegrationMethod, double &r_int_distance_dot_n)

    {

        const unsigned int n_int_pts = r_geom.IntegrationPointsNumber(IntegrationMethod);

        double distance = 0.0;
        array_1d<double, 3> distance_vector;
        CoordinatesArrayType proj_global_coords;
        double distance_dot_n;

        // Get the  Gauss points coordinates
        IntegrationPointsArrayType gauss_pts;
        gauss_pts = r_geom.IntegrationPoints(IntegrationMethod);

        Vector jacobians_values;
        r_geom.DeterminantOfJacobian(jacobians_values, IntegrationMethod);

        // Compute the x dot n at every gauss points
        for (unsigned int i_gauss = 0; i_gauss < n_int_pts; ++i_gauss)
        {

            array_1d<double, 3> area_normal = r_geom.UnitNormal(gauss_pts[i_gauss].Coordinates());

            // Compute the global coordinates of the Gauss pt.
            CoordinatesArrayType global_coords = ZeroVector(3);

            global_coords = r_geom.GlobalCoordinates(global_coords, gauss_pts[i_gauss].Coordinates());

            distance = r_cutting_plane.DistanceTo(global_coords);

            noalias(distance_vector) = distance * r_cutting_plane.mNormal;

            distance_dot_n = MathUtils<double>::Dot(distance_vector, area_normal);

            r_int_distance_dot_n += distance_dot_n * gauss_pts[i_gauss].Weight() * jacobians_values[i_gauss];
        }
    }

    void CalculateIntAreaDotNplane(GeometryType &r_geom, Plane &r_cutting_plane, const IntegrationMethodType IntegrationMethod, double &r_int_area_dot_n_plane)

    {

        const unsigned int n_int_pts = r_geom.IntegrationPointsNumber(IntegrationMethod);

        double area_dot_n_normal;

        // Get the  Gauss points coordinates
        IntegrationPointsArrayType gauss_pts;
        gauss_pts = r_geom.IntegrationPoints(IntegrationMethod);

        Vector jacobians_values;
        r_geom.DeterminantOfJacobian(jacobians_values, IntegrationMethod);

        // Compute the area dot plane normal at every gauss points
        for (unsigned int i_gauss = 0; i_gauss < n_int_pts; ++i_gauss)
        {

            array_1d<double, 3> area_normal = r_geom.AreaNormal(gauss_pts[i_gauss].Coordinates());

            area_dot_n_normal = MathUtils<double>::Dot(area_normal, r_cutting_plane.mNormal);

            r_int_area_dot_n_plane += area_dot_n_normal;
        }
    }

    void IsNegativeOrSplit(GeometryType &r_geom, bool &r_is_split, bool &r_is_negative)
    {

        unsigned int n_pos = 0, n_neg = 0, n_zero = 0;
        double nodal_distance = 0.0;

        for (unsigned int i = 0; i < r_geom.size(); ++i)
        {
            nodal_distance = r_geom[i].FastGetSolutionStepValue(DISTANCE);

            if (nodal_distance < 0.0 - tol)
            {
                n_neg++;
            }
            else if (nodal_distance > 0.0 + tol)
            {
                n_pos++;
            }

            else
            {
                n_zero++;
                n_neg++;
            }
        }

        if ((n_pos > 0) && (n_neg > 0))
            r_is_split = true;

        else
            r_is_split = false;

        if (n_zero > 0)
            r_is_split = false;

        if (n_neg == r_geom.size())
            r_is_negative = true;

        else
            r_is_negative = false;
    }

    double GetVolume()
    {
        return m_volume;
    }

    double GetIntersectedArea()
    {
        return m_intersected_area;
    }

    ModelPart &GetPlaneModelPart()
    {
        return this->m_plane_model_part;
    }

    void UpdatePositionOfPlaneBasedOnTargetVolume(ModelPart &r_model_part, double v_target, double max_rel_res, std::string type, unsigned int max_iterations)
    {

        double vol_residual = v_target - CalculateVolume(r_model_part);
        double vol_inter_residual = 0.0;
        double vol_rel_residual = vol_residual / v_target;
        double displacement_value = 0.0;

        ModelPart::ConditionsContainerType::iterator icond = m_plane_model_part.ConditionsBegin();

        CoordinatesArrayType p0 = icond->GetGeometry()[0].Coordinates();
        CoordinatesArrayType p1 = icond->GetGeometry()[1].Coordinates();
        CoordinatesArrayType p2 = icond->GetGeometry()[2].Coordinates();
        Plane cutting_plane(p0, p1, p2);
        double vol = 0.0;
        double vol_inter = 0.0;
        double area = 0.0;
        array_1d<double, 3> displacement_vector;
        unsigned int iteration_nr = 0;
        double movement = 0.0;

        if (type == "NewtonRaphson")
        {
            while (fabs(vol_rel_residual) >= max_rel_res && iteration_nr <= max_iterations)
            {
                vol = CalculateVolume(r_model_part);
                area = m_intersected_area;
                if (area <= tol)
                    KRATOS_ERROR << "Intersected area is zero, Newton Raphson method failed. Try with the different initial position of the plane" << std::endl;
                vol_residual = (v_target - vol);
                vol_rel_residual = vol_residual / v_target;
                
                displacement_value = vol_residual / area;
                displacement_vector = displacement_value * cutting_plane.mNormal;
                UpdatePlaneModelPart(displacement_vector);
                movement += displacement_value;
                std::cout << "Iteration Nr. :: " << iteration_nr << std::endl;
                std::cout << "Relative Volume residual :: " << vol_rel_residual << std::endl;
                std::cout << "Volume :: " << vol << std::endl;
                std::cout << "Movement :: " << movement << std::endl;
                if (iteration_nr == max_iterations)
                    std::cout<<"Max iterations reached. Newton Raphson method didn't converge"<<::std::endl;
                iteration_nr += 1;
            }
        }

        else if (type == "LeapFroggingNewton")
        {
            while (fabs(vol_rel_residual) >= max_rel_res && fabs(vol_residual - vol_inter_residual) >= tol && iteration_nr <= max_iterations)
            {
                vol = CalculateVolume(r_model_part);
                area = m_intersected_area;
                if (area <= tol)
                    KRATOS_ERROR << "Intersected area is zero, Leapfrogging Newton method failed. Try with the different initial position of the plane" << std::endl;
                vol_residual = (v_target - vol);
                vol_rel_residual = vol_residual / v_target;
                displacement_value = vol_residual / area;
                displacement_vector = displacement_value * cutting_plane.mNormal;
                UpdatePlaneModelPart(displacement_vector);
                vol_inter = CalculateVolume(r_model_part);
                displacement_vector *= -1;
                UpdatePlaneModelPart(displacement_vector);
                vol_inter_residual = (v_target - vol_inter);
                
                displacement_value = (vol_residual * vol_residual) / ((vol_residual - vol_inter_residual) * area);
                displacement_vector = displacement_value * cutting_plane.mNormal;
                UpdatePlaneModelPart(displacement_vector);
                movement += displacement_value;
                std::cout << "Iteration Nr. :: " << iteration_nr << std::endl;
                std::cout << "Relative Volume residual :: " << vol_rel_residual << std::endl;
                std::cout << "Movement :: " << movement << std::endl;
                std::cout << "Volume :: " << vol << std::endl;

                if (iteration_nr == max_iterations)
                    std::cout<<"Max iterations reached. Leap Frogging Newton method didn't converge"<<::std::endl;


                iteration_nr += 1;
                
            }
        }

        else
            KRATOS_ERROR << "String specifying NewtonRaphson or LeapFroggingNewton is required as argument" << std::endl;
    }

    void UpdatePlaneModelPart(array_1d<double, 3> &r_displacement_vector)
    {

        for (ModelPart::NodesContainerType::iterator inode = m_plane_model_part.NodesBegin(); inode != m_plane_model_part.NodesEnd(); inode++)
        {
            inode->Coordinates() += r_displacement_vector;
        }
    }



    /* double CalculateIntersectedAreaOnPlane(ModelPart &r_model_part, ModelPart &r_plane_model_part)

    {
        KRATOS_TRY;

        std::cout << "Calculating surface area intersected by closed surface on plane" << std::endl;
        double area = 0.0;
        double int_area_dot_n_plane = 0.0;

        IntegrationMethodType IntegrationMethod = GeometryData::GI_GAUSS_1;

        // Calculation of cutting plane parameters
        ModelPart::ConditionsContainerType::iterator icond = r_plane_model_part.ConditionsBegin();

        CoordinatesArrayType p0 = icond->GetGeometry()[0].Coordinates();
        CoordinatesArrayType p1 = icond->GetGeometry()[1].Coordinates();
        CoordinatesArrayType p2 = icond->GetGeometry()[2].Coordinates();

        Plane cutting_plane(p0, p1, p2);

        // Calculation of distance of nodes of the model part from the cutting plane

        for (ModelPart::NodesContainerType::iterator inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); inode++)
        {

            double distance = cutting_plane.DistanceTo(inode->Coordinates());

            inode->FastGetSolutionStepValue(DISTANCE) = distance;
        }

        for (ModelPart::ElementsContainerType::iterator ielem = r_model_part.ElementsBegin(); ielem != r_model_part.ElementsEnd(); ielem++)
        {

            array_1d<double, 3> distances_vector;

            bool is_split = false;
            bool is_negative = false;

            GeometryType &geom = ielem->GetGeometry();

            for (unsigned int i = 0; i < geom.size(); ++i)
                distances_vector(i) = geom[i].FastGetSolutionStepValue(DISTANCE);

            ielem->SetValue(ELEMENTAL_DISTANCES, distances_vector);

            Vector &r_elemental_distances = ielem->GetValue(ELEMENTAL_DISTANCES);

            DivideTriangle2D3 triangle_splitter(geom, r_elemental_distances);

            IsNegativeOrSplit(geom, is_split, is_negative);

            if (is_split)
            {
                IsNegativeOrSplit(geom, is_split, is_negative);
                // Call the divide geometry method
                triangle_splitter.GenerateDivision();

                for (unsigned int i = 0; i < triangle_splitter.mNegativeSubdivisions.size(); i++)
                {

                    IndexedPointGeometryType indexed_subgeom = *(triangle_splitter.mNegativeSubdivisions[i]);

                    Node<3>::Pointer node_i = Node<3>::Pointer(new Node<3>(indexed_subgeom[0].Id(), indexed_subgeom[0].X(), indexed_subgeom[0].Y(), indexed_subgeom[0].Z()));

                    Node<3>::Pointer node_j = Node<3>::Pointer(new Node<3>(indexed_subgeom[1].Id(), indexed_subgeom[1].X(), indexed_subgeom[1].Y(), indexed_subgeom[1].Z()));
                    Node<3>::Pointer node_k = Node<3>::Pointer(new Node<3>(indexed_subgeom[2].Id(), indexed_subgeom[2].X(), indexed_subgeom[2].Y(), indexed_subgeom[2].Z()));

                    Triangle3D3<Node<3>>::Pointer subgeom = Triangle3D3<Node<3>>::Pointer(new Triangle3D3<Node<3>>(node_i, node_j, node_k));

                    CalculateIntAreaDotNplane(*subgeom, cutting_plane, IntegrationMethod, int_area_dot_n_plane);
                }
            }
            else
            {

                if (is_negative)
                    CalculateIntAreaDotNplane(geom, cutting_plane, IntegrationMethod, int_area_dot_n_plane);
            }
        }

        area = -int_area_dot_n_plane;

        return area;

        KRATOS_CATCH(" ");
    } */

    /* 
    double CalculateVolumeEnclosedByClosedSurface(ModelPart &r_model_part)

    {

        IntegrationMethodType IntegrationMethod = GeometryData::GI_GAUSS_1;

        double volume = 0.0;
        double int_coord_dot_n = 0.0;
        CoordinatesArrayType centre = ZeroVector(3);
        double n_nodes = r_model_part.Nodes().size();

        for (ModelPart::NodesContainerType::iterator inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); inode++)
        {

            noalias(centre) += inode->Coordinates();
        }

        centre = centre / n_nodes;

        for (ModelPart::ElementsContainerType::iterator ielem = r_model_part.ElementsBegin(); ielem != r_model_part.ElementsEnd(); ielem++)
        {

            GeometryType &geom = ielem->GetGeometry();

            CalculateIntCoordDotN(geom, centre, IntegrationMethod, int_coord_dot_n);
        }

        volume = int_coord_dot_n / 3;

        return volume;
    } */

    /*     void CalculateIntCoordDotN(GeometryType &r_geom, const CoordinatesArrayType centre, const IntegrationMethodType IntegrationMethod, double &r_int_coord_dot_n)

    {

        const unsigned int n_int_pts = r_geom.IntegrationPointsNumber(IntegrationMethod);

        double coord_dot_n;

        // Get the  Gauss points coordinates
        IntegrationPointsArrayType gauss_pts;
        gauss_pts = r_geom.IntegrationPoints(IntegrationMethod);

        Vector jacobians_values;
        r_geom.DeterminantOfJacobian(jacobians_values, IntegrationMethod);

        // Compute the x dot n at every gauss points
        for (unsigned int i_gauss = 0; i_gauss < n_int_pts; ++i_gauss)
        {

            array_1d<double, 3> area_normal = r_geom.UnitNormal(gauss_pts[i_gauss].Coordinates());

            // Compute the global coordinates of the Gauss pt.
            CoordinatesArrayType global_coords = ZeroVector(3);

            global_coords = r_geom.GlobalCoordinates(global_coords, gauss_pts[i_gauss].Coordinates());

            CorrectNormalDirection(area_normal, global_coords, centre);

            coord_dot_n = MathUtils<double>::Dot(global_coords, area_normal);

            r_int_coord_dot_n += coord_dot_n * gauss_pts[i_gauss].Weight() * jacobians_values[i_gauss];
        }
    } */

    /*         void CorrectNormalDirection(array_1d<double, 3> &r_normal, const CoordinatesArrayType r_cordinate, const CoordinatesArrayType r_centre)
    {

        CoordinatesArrayType r = r_cordinate - r_centre;

        if (MathUtils<double>::Dot(r, r_normal) > 0)
            r_normal *= 1;
        else
            r_normal *= -1;
    } */

    

    // VM !!!!!!!!!!

  private:
    double m_volume;
    double m_intersected_area;
    ModelPart &m_plane_model_part;
}; //class VolumeCalculationUnderPlaneUtility
} /* namespace Kratos.*/

#endif /* KRATOS_VOLUME_UTILITY_H_INCLUDED  defined */
