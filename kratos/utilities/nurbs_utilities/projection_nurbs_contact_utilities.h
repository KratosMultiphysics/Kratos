//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Andrea Gorgi
//

#if !defined(KRATOS_PROJECTION_NURBS_CONTACT_UTILITIES_H_INCLUDED)
#define KRATOS_PROJECTION_NURBS_CONTACT_UTILITIES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "containers/array_1d.h"

namespace Kratos
{
template<class TPointType, class TSurfaceContainerPointType> 
class ProjectionNurbsContactUtilities
{
public:

    /// Geometry as base class.
    typedef Geometry<TPointType> BaseType;
    typedef Geometry<TPointType> GeometryType;

    ///definition of the geometry type with given NodeType
    typedef Geometry<Node> GeometryNodeType;

    /// Pointer definition of CouplingGeometry
    KRATOS_CLASS_POINTER_DEFINITION( ProjectionNurbsContactUtilities );

    typedef TPointType PointType;

    typedef typename BaseType::IndexType IndexType;
    typedef typename BaseType::SizeType SizeType;

    typedef typename BaseType::IntegrationPointsArrayType IntegrationPointsArrayType;
    typedef typename BaseType::PointsArrayType PointsArrayType;
    typedef typename BaseType::CoordinatesArrayType CoordinatesArrayType;
    typedef std::vector<CoordinatesArrayType> CoordinatesArrayVectorType;
    typedef PointerVector<GeometryType> GeometriesArrayType;

    typedef NurbsSurfaceGeometry<3, TSurfaceContainerPointType> NurbsSurfaceType;
    typedef typename TSurfaceContainerPointType::value_type NodeType;

    typedef typename NurbsSurfaceType::Pointer NurbsSurfaceTypePointer;


    // MODIFIED --------------------------------------
    typedef PointerVector<Node> ContainerNodeType;
    typedef PointerVector<Point> ContainerEmbeddedNodeType;

    typedef BrepCurveOnSurface<ContainerNodeType, ContainerEmbeddedNodeType> BrepCurveOnSurfaceType;

    typedef DenseVector<typename BrepCurveOnSurfaceType::Pointer> BrepCurveOnSurfaceArrayType;

    /*
    * @brief Returns the projection of a point onto a Nurbs curve
    *        geometry using the Newton-Rapshon iterative method
    * @param rProjectedPointLoaclCoordinates Intial guess for the Newton-Rapshon algorithm
    *        overwritten by the local coordinates of the projected point onto
    *        the Nurbs curve geometry
    * @param rPoint The point to be projected onto the Nurbs curve geometry
    *        This is overwritten by the Cartesian coordinates of the projected
    *        point in case the projection is successful
    * @param rProjectedPointGlobalCoordinates The projection onto the Nurbs curve geometry
    * @param rNurbsCurve The Nurbs curve geometry onto which the point is
    *        to be projected
    * @param MaxIterations Maximum number of iterations for the Newton-Rapshon
    *        algorithm
    * @param Accuracy Accuracy for the the Newton-Rapshon algorithm
    */
    static bool NewtonRaphsonCurve(
        CoordinatesArrayType& rProjectedPointLocalCoordinates,
        const CoordinatesArrayType& rPointGlobalCoordinatesCoordinates,
        CoordinatesArrayType& rProjectedPointGlobalCoordinates,
        const Geometry<TPointType>& rGeometry,
        double& distance,
        const int MaxIterations = 20,
        const double Accuracy = 1e-6)
    {
        // Intialize variables
        double residual, delta_t;

        std::vector<array_1d<double, 3>> derivatives(3);
        array_1d<double, 3> distance_vector;

        bool projection_reset_to_boundary = false;
        const double Acc = 1e-6;
        // Loop over all Newton-Raphson iterations
        for (int i = 0; i < MaxIterations; ++i)
        {
            // Compute the position, the base and the acceleration vector
            rGeometry.GlobalSpaceDerivatives(
                derivatives,
                rProjectedPointLocalCoordinates,
                2);
            rProjectedPointGlobalCoordinates = derivatives[0];

            // Compute the distance vector between the point and its
            // projection on the curve
            distance_vector = rProjectedPointGlobalCoordinates - rPointGlobalCoordinatesCoordinates;

            distance = norm_2(distance_vector);
            if (norm_2(distance_vector) < Acc) // Acc
                return true;

            // Compute the residual
            residual = inner_prod(distance_vector, derivatives[1]);
            if (std::abs(residual) < Acc) // Acc
                return true;

            // Compute the increment
            delta_t = residual / (inner_prod(derivatives[2], distance_vector) + pow(norm_2(derivatives[1]), 2));

            // Increment the parametric coordinate
            rProjectedPointLocalCoordinates[0] -= delta_t;

            // Check if the increment is too small and if yes return true
            if (norm_2(delta_t * derivatives[1]) < Acc) // Acc
                return true;

            // Check if the parameter gets out of its interval of definition and if so clamp it
            // back to the boundaries
            int check = rGeometry.ClosestPointLocalToLocalSpace(
                rProjectedPointLocalCoordinates, rProjectedPointLocalCoordinates);
            if (check == 0) {
                if (projection_reset_to_boundary) { return false; }
                else { projection_reset_to_boundary = true; }
            }
        }

        // Return false if the Newton-Raphson iterations did not converge
        return false;
    }

    /*
    * @brief Returns the projection of a point onto a Nurbs surface
    *        geometry using the Newton-Rapshon iterative method
    * @param rProjectedPointLocalCoordinates Intial guess for the Newton-Rapshon algorithm
    *        overwritten by the local coordinates of the projected point onto
    *        the Nurbs surface geometry
    * @param rPoint The point to be projected onto the Nurbs surface geometry
    *        This is overwritten by the Cartesian coordinates of the projected
    *        point in case the projection is successful
    * @param rResult The projection onto the Nurbs surface geometry
    * @param rNurbsCurve The Nurbs curve geometry onto which the point is
    *        to be projected
    * @param MaxIterations Maximum number of iterations for the Newton-Rapshon
    *        algorithm
    * @param Accuracy Accuracy for the the Newton-Rapshon algorithm
    */
    template <int TDimension>
    static bool NewtonRaphsonSurface(
        CoordinatesArrayType& rProjectedPointLocalCoordinates,
        const CoordinatesArrayType& rPointGlobalCoordinates,
        CoordinatesArrayType& rProjectedPointGlobalCoordinates,
        const NurbsSurfaceGeometry<TDimension, TPointType>& rNurbsSurface,
        const int MaxIterations = 20,
        const double Accuracy = 1e-6)
    {
        // Initialize variables
        bool is_first_row_zero, is_second_row_zero, is_first_column_zero, is_second_column_zero, is_system_invertible;
        double d_u = 0.0;
        double d_v = 0.0;
        double xi_cos, eta_cos, residual_u, residual_v, j_00, j_01, j_11, det_j;

        // Loop over all the Newton-Raphson iterations
        for (int i = 0; i < MaxIterations; i++) {

            // Compute the position, the base and the acceleration vectors
            std::vector<array_1d<double, 3>> s;
            rNurbsSurface.GlobalSpaceDerivatives(s, rProjectedPointLocalCoordinates, 2);
            rProjectedPointGlobalCoordinates = s[0];

            // Compute the distance vector
            const array_1d<double, 3> distance_vector = s[0] - rPointGlobalCoordinates;

            // Compute the distance
            const double distance = norm_2(distance_vector);
            if (distance < Accuracy)
                return true;

            // Compute the residuals along both parametric directions
            residual_u = -inner_prod(s[1], distance_vector);
            residual_v = -inner_prod(s[2], distance_vector);

            // Compute the cosine with respect to the u-parametric coordinate
            xi_cos = std::abs(residual_u) / norm_2(s[1]) / norm_2(distance_vector);

            // Compute the cosine with respect to the v-parametric coordinate
            eta_cos = std::abs(residual_v) / norm_2(s[2]) / norm_2(distance_vector);

            // Check the orthogonality condition
            if (xi_cos < Accuracy && eta_cos < Accuracy)
                return true;

            // Compute the Jacobian of the nonlinear system
            j_00 = inner_prod(s[1], s[1]) + inner_prod(s[3], distance_vector);
            j_01 = inner_prod(s[1], s[2]) + inner_prod(s[4], distance_vector);
            j_11 = inner_prod(s[2], s[2]) + inner_prod(s[5], distance_vector);

            // Check for singularities otherwise update the parametric coordinates as usual
            is_first_row_zero = false;
            if ((std::abs(j_00) < Accuracy && std::abs(j_01) < Accuracy)) {
                is_first_row_zero = true;
            }
            is_second_row_zero = false;
            if (std::abs(j_01) < Accuracy && fabs(j_11) < Accuracy) {
                is_second_row_zero = true;
            }
            is_first_column_zero = false;
            if ((std::abs(j_00) < Accuracy && std::abs(j_01) < Accuracy)) {
                is_first_column_zero = true;
            }
            is_second_column_zero = false;
            if ((std::abs(j_01) < Accuracy && std::abs(j_11) < Accuracy)) {
                is_second_column_zero = true;
            }

            // Check if the system is solvable by checking the condition of the diagonal entries
            is_system_invertible = true;
            if (is_first_row_zero || is_second_row_zero || is_first_column_zero || is_second_column_zero) {
                is_system_invertible = false;
            }

            // Solve the 2x2 linear equation system and take into account special cases where singularities occur
            if (is_system_invertible) {
                det_j = j_00 * j_11 - j_01 * j_01;
                d_u = -(residual_v * j_01 - residual_u * j_11) / det_j;
                d_v = -(residual_u * j_01 - residual_v * j_00) / det_j;
            }
            else {
                if (is_first_row_zero) {
                    d_u = residual_v / j_11;
                    d_v = 0.0;
                }
                else if (is_second_row_zero) {
                    d_u = residual_u / j_00;
                    d_v = 0.0;
                }
                else if (is_first_column_zero) {
                    d_v = (residual_u + residual_v) / (j_01 + j_11);
                    d_u = 0.0;
                }
                else if (is_second_column_zero) {
                    d_u = (residual_u + residual_v) / (j_00 + j_01);
                    d_v = 0.0;
                }
            }

            // Check if the step size is too small
            if (norm_2(d_u * s[1] + d_v * s[2]) < Accuracy)
                return true;

            // Update the parametric coordinates
            rProjectedPointLocalCoordinates[0] += d_u;
            rProjectedPointLocalCoordinates[1] += d_v;

            // Check if the parametric coordinates get out of their interval of definition
            // and if so clamp them back to their boundaries
            rNurbsSurface.DomainIntervalU().IsInside(rProjectedPointLocalCoordinates[0]);
            rNurbsSurface.DomainIntervalV().IsInside(rProjectedPointLocalCoordinates[1]);
        }

        return false;
    }


    //--------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------
    static void GetDisplacement(NurbsSurfaceTypePointer rpNurbsSurfacePaired, CoordinatesArrayType& rParameterPointCoord, Vector& displacement) {
        // COMPUTE DISPLACEMENT AND DISPLACMENT DERIVATIVE IN THE PARAMETER SPACE (OF THE SURFACE)

        NurbsSurfaceShapeFunction shape_function_container(
            rpNurbsSurfacePaired->PolynomialDegreeU(), rpNurbsSurfacePaired->PolynomialDegreeV(),
            1);


        // GET DERIVATIVES
        shape_function_container.ComputeBSplineShapeFunctionValues(
                            rpNurbsSurfacePaired->KnotsU(), rpNurbsSurfacePaired->KnotsV(),
                            rParameterPointCoord[0], rParameterPointCoord[1]);

        SizeType num_nonzero_cps = shape_function_container.NumberOfNonzeroControlPoints();
    
        /// Get List of Control Points
        PointsArrayType nonzero_control_points(num_nonzero_cps);
        auto cp_indices = shape_function_container.ControlPointIndices(
            rpNurbsSurfacePaired->NumberOfControlPointsU(), rpNurbsSurfacePaired->NumberOfControlPointsV());

        for (IndexType j = 0; j < num_nonzero_cps; j++) {
            nonzero_control_points(j) = rpNurbsSurfacePaired->pGetPoint(cp_indices[j]);
        }

        displacement = ZeroVector(3); //dim

        // COMPUTE DISPLACEMENT
        for (IndexType j = 0; j < num_nonzero_cps; j++) {
                
            displacement += nonzero_control_points[j].GetSolutionStepValue(DISPLACEMENT) * shape_function_container(j,0);
        }
    }

    /**
     * @brief Get the Displacement Derivatives On Curve object
     * 
     * @param rpNurbsSurfacePaired 
     * @param rParameterPointCoord 
     * @param n_derivatives 
     * @param displacement 
     * @param displacement_derivatives 
     * @param local_space_derivatives_on_curve 
     */
    static void GetDisplacementDerivativesOnCurve(NurbsSurfaceTypePointer rpNurbsSurfacePaired, 
                                                    CoordinatesArrayType& rParameterPointCoord,
                                                    const SizeType& n_derivatives, 
                                                    std::vector<Vector>& displacement_derivatives,
                                                    const std::vector<CoordinatesArrayType>& local_space_derivatives_on_curve) 
    {

        DenseVector<Matrix> displacement_derivatives_global; // [[du/dx, du/dy]
                                                              //  [dv/dx, dv/dy]]

        displacement_derivatives.resize(n_derivatives+1); 
        displacement_derivatives[0].resize(3);

        GetDisplacementDerivatives(rpNurbsSurfacePaired, rParameterPointCoord, n_derivatives, displacement_derivatives[0], displacement_derivatives_global);

        

        if (n_derivatives > 2) KRATOS_ERROR << "ERROR: du/dt not implemented for third order derivatives or higher \n";

        for (IndexType i_der = 1; i_der < n_derivatives+1; i_der++)
        {
            displacement_derivatives[i_der].resize(3);
            for (IndexType i_dim = 0; i_dim < 2; i_dim++)
            {
                // du/dt^n = sum(0,k) (n)  
                //                    (k)

                if (i_der == 1) {
                    displacement_derivatives[1][i_dim] = displacement_derivatives_global[0](i_dim, 0)*local_space_derivatives_on_curve[1][0] + displacement_derivatives_global[0](i_dim, 1)*local_space_derivatives_on_curve[1][1];
                }
                else if (i_der == 2)
                {
                    displacement_derivatives[2][i_dim] = displacement_derivatives_global[1](i_dim, 0) * pow(local_space_derivatives_on_curve[1][0],2) 
                                                         + 2*displacement_derivatives_global[1](i_dim, 1) * local_space_derivatives_on_curve[1][0] * local_space_derivatives_on_curve[1][1]
                                                         + displacement_derivatives_global[1](i_dim, 2) * pow(local_space_derivatives_on_curve[1][1], 2)
                                                         + displacement_derivatives_global[0](i_dim, 0) * local_space_derivatives_on_curve[2][0] 
                                                         + displacement_derivatives_global[0](i_dim, 1) * local_space_derivatives_on_curve[2][1];
                }

                
                
            }
        }


    }


    static void GetDisplacementDerivatives(NurbsSurfaceTypePointer rpNurbsSurfacePaired, 
                                    CoordinatesArrayType& rParameterPointCoord,
                                    const SizeType& n_derivatives, 
                                    Vector& displacement,
                                    DenseVector<Matrix>& displacement_derivatives) {
        // COMPUTE DISPLACEMENT AND DISPLACMENT DERIVATIVE IN THE PARAMETER SPACE (OF THE SURFACE)

        NurbsSurfaceShapeFunction shape_function_container(
            rpNurbsSurfacePaired->PolynomialDegreeU(), rpNurbsSurfacePaired->PolynomialDegreeV(),
            n_derivatives);

        // GET DERIVATIVES
        shape_function_container.ComputeBSplineShapeFunctionValues(
                            rpNurbsSurfacePaired->KnotsU(), rpNurbsSurfacePaired->KnotsV(),
                            rParameterPointCoord[0], rParameterPointCoord[1]);

        SizeType num_nonzero_cps = shape_function_container.NumberOfNonzeroControlPoints();
    
        /// Get List of Control Points
        PointsArrayType nonzero_control_points(num_nonzero_cps);
        auto cp_indices = shape_function_container.ControlPointIndices(
            rpNurbsSurfacePaired->NumberOfControlPointsU(), rpNurbsSurfacePaired->NumberOfControlPointsV());

        for (IndexType j = 0; j < num_nonzero_cps; j++) {
            nonzero_control_points(j) = rpNurbsSurfacePaired->pGetPoint(cp_indices[j]);
        }

        displacement = ZeroVector(3); //dim
        displacement_derivatives.resize(n_derivatives);

        // COMPUTE DISPLACEMENT
        for (IndexType j = 0; j < num_nonzero_cps; j++) {
                
            displacement += nonzero_control_points[j].GetSolutionStepValue(DISPLACEMENT) * shape_function_container(j,0);

        }

        // COMPUTE DERIVATIVES
        for (IndexType i = 0; i < n_derivatives; i++) {
            displacement_derivatives[i] = ZeroMatrix(3, i+2);
        }

        IndexType shape_derivative_index = 1;
        for (IndexType i_der = 0; i_der < n_derivatives; i_der++) 
        {
            for (IndexType k = 0; k < i_der + 2; k++) {
                for (IndexType j = 0; j < num_nonzero_cps; j++) {
                    const Vector curr_displacement_coefficient = nonzero_control_points[j].GetSolutionStepValue(DISPLACEMENT);
                    for (IndexType i_dim = 0; i_dim < 3; i_dim++) 
                    {
                        displacement_derivatives[i_der](i_dim, k) += curr_displacement_coefficient[i_dim] * shape_function_container(j,shape_derivative_index+k);

                    }
                }
            }
            shape_derivative_index += i_der + 2;
        }
    }


    ///@name Projection functionalaties
    ///@{

    static bool GetProjection(const NurbsSurfaceTypePointer rpNurbsSurfaceParent, 
                              const NurbsSurfaceTypePointer rpNurbsSurfacePaired, 
                              const CoordinatesArrayType& parentPointLocalCoord, 
                                BrepCurveOnSurfaceType &parent_geometry, 
                                BrepCurveOnSurfaceType &paired_geometry, 
                                double& localProjection, double& distance) {
        
        CoordinatesArrayType parentPointGlobalCoord;
        CoordinatesArrayType parentPointParamCoord;
        CoordinatesArrayType parentPointGlobalCoord_updated;

        array_1d<double, 3> old_normal = parent_geometry.Normal(parentPointLocalCoord);

        parent_geometry.GlobalCoordinates(parentPointGlobalCoord, parentPointLocalCoord);

        parent_geometry.LocalCoordinates(parentPointParamCoord, parentPointLocalCoord);

        // GET DERIVATIVES

        IntegrationPoint<1> integrationPointParent(parentPointLocalCoord[0]);
        std::vector<CoordinatesArrayType> local_space_derivatives(3);
        // CoordinatesArrayType param_coordinates(2);

        parent_geometry.LocalSpaceDerivatives( local_space_derivatives, integrationPointParent, 1);

        std::vector<CoordinatesArrayType> surface_global_space_derivatives(2);
        rpNurbsSurfaceParent->GlobalSpaceDerivatives(surface_global_space_derivatives, parentPointParamCoord,2);

        Matrix Jacobian = ZeroMatrix(2,2);
        Jacobian(0,0) = surface_global_space_derivatives[1][0];
        Jacobian(0,1) = surface_global_space_derivatives[2][0];
        Jacobian(1,0) = surface_global_space_derivatives[1][1];
        Jacobian(1,1) = surface_global_space_derivatives[2][1];

        Vector tangent_undeformed_physical(2); 
        tangent_undeformed_physical[0] = local_space_derivatives[1][0]; tangent_undeformed_physical[1] = local_space_derivatives[1][1];


        tangent_undeformed_physical = prod(Jacobian,tangent_undeformed_physical);

        //***** */

        std::vector<Vector> displacement_derivatives;

        GetDisplacementDerivativesOnCurve(rpNurbsSurfaceParent, parentPointParamCoord, 1, displacement_derivatives,
                                          local_space_derivatives); //get [du/dt, dv/dt]

        Vector normal_physical_deformed = ZeroVector(3); 
        normal_physical_deformed[0] = displacement_derivatives[1][1] + tangent_undeformed_physical[1];
        normal_physical_deformed[1] = -(displacement_derivatives[1][0] + tangent_undeformed_physical[0]);
        normal_physical_deformed[2] = 0.0;

        CoordinatesArrayType new_normal = (normal_physical_deformed)/ norm_2(normal_physical_deformed);


        parentPointGlobalCoord_updated = parentPointGlobalCoord + displacement_derivatives[0];

        // find slave normal in that point
        // CoordinatesArrayType normal = parent_geometry.UnitNormal(slavePointLocalCoord);
        
        const double toll = 1e-8;
        double res = toll + 1;
        int it = 0;
        const int itMax = 100;
        const double number_of_initial_guesses = 50;

        // starting point NEWTON RAPHSON
        CoordinatesArrayType globCoordPaired;
        CoordinatesArrayType paramCoordPaired;
        CoordinatesArrayType globCoordPaired_updated;

        Vector displacement_paired;
        std::vector<Vector> displacement_derivatives_paired;

        Vector displacement_parent;
        CoordinatesArrayType t(3); 

        Vector b_rep_interval;
        paired_geometry.DomainInterval(b_rep_interval);

        bool has_converged_at_least_once = false;
        distance = 1e12;
        
        //first trial
        // paired_geometry.ProjectionPointGlobalToLocalSpace(slavePoint_updated, t, 1e-12);

        // KRATOS_WATCH(t)
        int larger_direction = 1;
        int smaller_direction = 0;

        //

        if (std::abs(new_normal[0]) > std::abs(new_normal[1])) { larger_direction = 0; smaller_direction = 1;}

        for (IndexType i_guess = 0; i_guess < number_of_initial_guesses; i_guess++)
        {   
            res = toll + 1;
            it = 0;
            t[0] = b_rep_interval[0] +  (b_rep_interval[1]- b_rep_interval[0])/(number_of_initial_guesses-1) * float(i_guess);

            paired_geometry.GlobalCoordinates(globCoordPaired, t);

            paired_geometry.LocalCoordinates(paramCoordPaired, t);

            //
            
            //
            GetDisplacement(rpNurbsSurfacePaired, paramCoordPaired, displacement_paired); //get slave displacement
            globCoordPaired_updated = globCoordPaired + displacement_paired;
            //

            double s_line = (globCoordPaired_updated[larger_direction] - parentPointGlobalCoord_updated[larger_direction])/new_normal[larger_direction];
            double smaller_direction_line_current = s_line*new_normal[smaller_direction] + parentPointGlobalCoord_updated[smaller_direction];
            res = smaller_direction_line_current - globCoordPaired_updated[smaller_direction];

            while (std::abs(res) > toll && it < itMax) {

                //----------------
                std::vector<CoordinatesArrayType> rGlobalSpaceDerivatives;

                paired_geometry.GlobalSpaceDerivatives(rGlobalSpaceDerivatives, t, 1);

                paired_geometry.LocalCoordinates(paramCoordPaired, t);

                std::vector<CoordinatesArrayType> rLocalSpaceDerivatives;

                paired_geometry.LocalSpaceDerivatives(rLocalSpaceDerivatives, t, 1);
                
                // NEW

                GetDisplacementDerivativesOnCurve(rpNurbsSurfacePaired, paramCoordPaired, 1, displacement_derivatives_paired,
                                                rLocalSpaceDerivatives);
                //

                Vector paired_derivative_updated = rGlobalSpaceDerivatives[1] + displacement_derivatives_paired[1];

                const double f_der = new_normal[smaller_direction]/new_normal[larger_direction] * paired_derivative_updated[larger_direction] - paired_derivative_updated[smaller_direction];

                if (std::abs(f_der) < 1e-10) break;
                t[0] -= res/f_der;

                if (t[0] < b_rep_interval[0] || t[0] > b_rep_interval[1]) break;

                paired_geometry.GlobalCoordinates(globCoordPaired, t);

                paired_geometry.LocalCoordinates(paramCoordPaired, t);

                GetDisplacement(rpNurbsSurfacePaired, paramCoordPaired, displacement_paired);

                globCoordPaired_updated = globCoordPaired + displacement_paired;
                //
                

                s_line = (globCoordPaired_updated[larger_direction] - parentPointGlobalCoord_updated[larger_direction])/new_normal[larger_direction];

                smaller_direction_line_current = s_line*new_normal[smaller_direction] + parentPointGlobalCoord_updated[smaller_direction];

                res = smaller_direction_line_current - globCoordPaired_updated[smaller_direction];

                it += 1;
            }

            if (std::abs(res) < toll) {
                has_converged_at_least_once = true; 

                double curr_distance = norm_2(parentPointGlobalCoord_updated-globCoordPaired_updated);
                if (curr_distance < distance) {
                    distance = curr_distance;
                    localProjection = t[0];
                }
            }
        }
        
        


        if (!has_converged_at_least_once) 
        {
            // std::string name_output_file; 
            // if (parentPointGlobalCoord[0] > 5.0+1e-12) {
            //     name_output_file = "txt_files/error_projection_contact_slave.txt";
            // } else {
            //     name_output_file = "txt_files/error_projection_contact_master.txt";
            // }

            // std::ofstream outputFile(name_output_file, std::ios::app);
            // if (!outputFile.is_open())
            // {
            //     std::cerr << "Failed to open the file for writing." << std::endl;
            //     return 3;
            // }

            // t[0] = localProjection;
            // paired_geometry.GlobalCoordinates(globCoordPaired, t);
            // paired_geometry.LocalCoordinates(paramCoordPaired, t);

            // GetDisplacement(rpNurbsSurfacePaired, paramCoordPaired, displacement_paired);

            // globCoordPaired_updated = globCoordPaired + displacement_paired;

            // outputFile << std::setprecision(14); // Set precision to 10^-14
            // outputFile << parentPointGlobalCoord_updated[0] << "  " << parentPointGlobalCoord_updated[1] << "  "
            //             << globCoordPaired_updated[0] << "  " << globCoordPaired_updated[1] << "  " 
            //             << new_normal[0] << "  " << new_normal[1] << "\n";
            return 0;
        }
        else
        {
            std::string name_output_file; 
            if (parentPointGlobalCoord[0] > 5.0+1e-12) {
                name_output_file = "txt_files/projection_contact_slave.txt";
            } else {
                name_output_file = "txt_files/projection_contact_master.txt";
            }

            std::ofstream outputFile(name_output_file, std::ios::app);
            if (!outputFile.is_open())
            {
                std::cerr << "Failed to open the file for writing." << std::endl;
                return 3;
            }

            t[0] = localProjection;
            paired_geometry.GlobalCoordinates(globCoordPaired, t);
            paired_geometry.LocalCoordinates(paramCoordPaired, t);

            GetDisplacement(rpNurbsSurfacePaired, paramCoordPaired, displacement_paired);

            globCoordPaired_updated = globCoordPaired + displacement_paired;

            outputFile << std::setprecision(14); // Set precision to 10^-14
            outputFile << parentPointGlobalCoord_updated[0] << "  " << parentPointGlobalCoord_updated[1] << "  "
                        << globCoordPaired_updated[0] << "  " << globCoordPaired_updated[1] << "  " 
                        << new_normal[0] << "  " << new_normal[1] << "\n";

            // outputFile << parentPointGlobalCoord[0] << "  " << parentPointGlobalCoord[1] << "  "
            //             << globCoordPaired[0] << "  " << globCoordPaired[1] << "  " 
            //             << new_normal[0] << "  " << new_normal[1] << "\n";
            outputFile.close();
            
        }

        return 1;
    }    
    /// }




    /*
    * @brief Returns the projection of a point onto a Nurbs curve
    *        geometry using the Newton-Rapshon iterative method
    * @param rProjectedPointLoaclCoordinates Intial guess for the Newton-Rapshon algorithm
    *        overwritten by the local coordinates of the projected point onto
    *        the Nurbs curve geometry
    * @param rPoint The point to be projected onto the Nurbs curve geometry
    *        This is overwritten by the Cartesian coordinates of the projected
    *        point in case the projection is successful
    * @param rProjectedPointGlobalCoordinates The projection onto the Nurbs curve geometry
    * @param rNurbsCurve The Nurbs curve geometry onto which the point is
    *        to be projected
    * @param MaxIterations Maximum number of iterations for the Newton-Rapshon
    *        algorithm
    * @param Accuracy Accuracy for the the Newton-Rapshon algorithm
    */
    static bool NewtonRaphsonCurveOnDeformed(
        const NurbsSurfaceTypePointer rpNurbsSurfaceParent, 
        const NurbsSurfaceTypePointer rpNurbsSurfacePaired, 
        CoordinatesArrayType& rProjectedPointLocalCoordinates,
        const CoordinatesArrayType& rPointLocalCoordinates, // curve on surface local coord (curve parameter space)
        CoordinatesArrayType& rProjectedPointGlobalCoordinates,
        BrepCurveOnSurfaceType &rParentGeometry, BrepCurveOnSurfaceType &rPairedGeometry, 
        double& distance,
        const int rNumberOfInitialGuesses,
        const int MaxIterations = 20,
        const double Accuracy = 1e-6)
    {

        // Intialize variables
        double residual, delta_t;

        CoordinatesArrayType projected_point_local_coordinates(3);

        std::vector<array_1d<double, 3>> global_derivatives(3);

        std::vector<array_1d<double, 3>> local_derivatives(3);

        std::vector<Vector> derivatives_on_deformed;
        array_1d<double, 3> distance_vector;

        const double Acc = Accuracy;

        CoordinatesArrayType point_local_coordinates(3); //parameter
        CoordinatesArrayType point_global_coordinates(3);
        Vector point_global_displacement;
        rParentGeometry.GlobalCoordinates(point_global_coordinates, rPointLocalCoordinates);

        rParentGeometry.LocalCoordinates(point_local_coordinates, rPointLocalCoordinates);
        GetDisplacement(rpNurbsSurfaceParent, point_local_coordinates, point_global_displacement); 


        CoordinatesArrayType point_deformed_global_coordinates = point_global_coordinates + point_global_displacement;

        point_deformed_global_coordinates[2] = 0.0;

        Vector projected_point_global_displacement;
        CoordinatesArrayType projected_point_deformed_global_coordinates(3);

        Vector b_rep_interval;
        rPairedGeometry.DomainInterval(b_rep_interval);
        bool is_converged_at_least_once = false;
        double best_distance = 1e12;

        // KRATOS_WATCH(point_deformed_global_coordinates)
        // KRATOS_WATCH(point_deformed_global_coordinates)
        for (IndexType i_guess = 0; i_guess < rNumberOfInitialGuesses; i_guess++)
        {   
            residual = Acc + 1;
            CoordinatesArrayType t = ZeroVector(3);
            t[0] = b_rep_interval[0] +  (b_rep_interval[1]- b_rep_interval[0])/(rNumberOfInitialGuesses-1) * float(i_guess);

            bool projection_reset_to_boundary = false;

            // Loop over all Newton-Raphson iterations
            for (int i = 0; i < MaxIterations; ++i)
            {
                // Compute the position, the base and the acceleration vector
                rPairedGeometry.GlobalSpaceDerivatives(
                    global_derivatives,
                    t,
                    2);
                    
                rProjectedPointGlobalCoordinates = global_derivatives[0];

                rPairedGeometry.LocalCoordinates(projected_point_local_coordinates, t);

                projected_point_local_coordinates[2] = 0.0;
                
                std::vector<array_1d<double, 3>> derivatives_updated(3);

                rPairedGeometry.LocalSpaceDerivatives(
                    local_derivatives,
                    t,
                    2);

                // NEW
                GetDisplacementDerivativesOnCurve(rpNurbsSurfacePaired, projected_point_local_coordinates, 2, derivatives_on_deformed,
                                                  local_derivatives);
                
                for (int i_der = 0; i_der < 3; i_der++){
                    for (int i_dim = 0; i_dim < 2; i_dim++) {
                        derivatives_updated[i_der][i_dim] =  global_derivatives[i_der][i_dim] + derivatives_on_deformed[i_der][i_dim];
                    }
                }

                derivatives_updated[0][2] = 0.0; derivatives_updated[1][2] = 0.0; derivatives_updated[2][2] = 0.0;

                // GetDisplacement(rpNurbsSurfacePaired, projected_point_local_coordinates, projected_point_global_displacement);  //???

                projected_point_deformed_global_coordinates = rProjectedPointGlobalCoordinates + derivatives_on_deformed[0];

                projected_point_deformed_global_coordinates[2] = 0.0;

                // Compute the distance vector between the point and its
                // projection on the curve
                distance_vector = projected_point_deformed_global_coordinates - point_deformed_global_coordinates;
                distance = norm_2(distance_vector);

                if (distance < Acc) // Acc
                {
                    rProjectedPointLocalCoordinates = t;
                    return true;
                }

                // Compute the residual
                residual = inner_prod(distance_vector, derivatives_updated[1]);
                if (std::abs(residual) < Acc) // Acc
                {
                    if (std::isnan(residual)) break;
                    else 
                    {
                        is_converged_at_least_once = true;

                        if (distance < best_distance) {
                            best_distance = distance;
                            rProjectedPointLocalCoordinates = t;
                        }
                        break;
                    }
                }

                // Compute the increment
                delta_t = residual / (inner_prod(derivatives_updated[2], distance_vector) + pow(norm_2(derivatives_updated[1]), 2));

                // Increment the parametric coordinate
                t[0] -= delta_t;

                // Check if the increment is too small and if yes return true
                if (norm_2(delta_t * derivatives_updated[1]) < Acc) // Acc
                {
                    double to_check = norm_2(delta_t * derivatives_updated[1]);
                    if (std::isnan(to_check)) break;
                    else 
                    {
                        rProjectedPointLocalCoordinates = t;
                        return true;
                    }
                }
                    

                // Check if the parameter gets out of its interval of definition and if so clamp it
                // back to the boundaries
                // CoordinatesArrayType closest_t = ZeroVector(3);
                // int check = rPairedGeometry.ClosestPointLocalToLocalSpace(
                //     t, t,2);

                int check = 0;
                double b_rep_interval_min; double b_rep_interval_max; 
                // KRATOS_WATCH(t)
                if (b_rep_interval[0] < b_rep_interval[1]) {
                    b_rep_interval_min = b_rep_interval[0];
                    b_rep_interval_max = b_rep_interval[1];
                }
                else {
                    b_rep_interval_min = b_rep_interval[1];
                    b_rep_interval_max = b_rep_interval[0];
                }
                //----
                if (t[0] > b_rep_interval_min && t[0] < b_rep_interval_max) check = 1;
                else if (std::abs(t[0] - b_rep_interval_min) < 1e-6) {
                    // t[0] = b_rep_interval_min; 
                    check = 2;
                }
                else if (std::abs(t[0] - b_rep_interval_max) < 1e-6) {
                    // t[0] = b_rep_interval_max; 
                    check = 2;
                }

                if (check == 0) {
                    // if (projection_reset_to_boundary) { return false; }
                    // else { projection_reset_to_boundary = true; }

                    break;
                }
            }
        }

        if (is_converged_at_least_once) 
        {
            distance = best_distance;
            return true;
        }
        // Return false if the Newton-Raphson iterations did not converge
        return false;
    }



    static bool GetProjectionOnPairedGeometry(NurbsSurfaceTypePointer rpNurbsSurfaceParent,
                                              NurbsSurfaceTypePointer rpNurbsSurfacePaired,
                                              BrepCurveOnSurfaceType& rParentGeometry,
                                              const BrepCurveOnSurfaceArrayType& rPairedGeometryList,
                                              const CoordinatesArrayType& rParentPointLocal,
                                              CoordinatesArrayType& ProjectionOnPairedGeometry,
                                              int &best_brep_id_slave,
                                              int rNumberInitialGuesses = 25,
                                              int rMaxIt = 50,
                                              double toll = 1e-9
                                              ) 
    {
        // MASTER
        std::vector<CoordinatesArrayType> parameter_space_derivatives_master;
        std::vector<CoordinatesArrayType> physical_space_derivatives_master;

        rParentGeometry.LocalSpaceDerivatives(
                parameter_space_derivatives_master,
                rParentPointLocal,
                1);

        rParentGeometry.GlobalSpaceDerivatives(
                physical_space_derivatives_master,
                rParentPointLocal,
                1);

        // SLAVE
        CoordinatesArrayType master_quadrature_point_parameter(3);
        master_quadrature_point_parameter[0] = parameter_space_derivatives_master[0][0]; master_quadrature_point_parameter[1] = parameter_space_derivatives_master[0][1];
        master_quadrature_point_parameter[2] = 0.0;

        CoordinatesArrayType master_quadrature_point_physical(3);
        master_quadrature_point_physical[0] = physical_space_derivatives_master[0][0]; master_quadrature_point_physical[1] = physical_space_derivatives_master[0][1];
        master_quadrature_point_physical[2] = 0.0;

        Vector displacement_on_projection; Vector displacement_on_master_quadrature_point;
        GetDisplacement(rpNurbsSurfaceParent, master_quadrature_point_parameter, displacement_on_master_quadrature_point);
        CoordinatesArrayType master_quadrature_point_deformed = master_quadrature_point_physical + displacement_on_master_quadrature_point;

        //---------------------------------------------------------
        double best_distance = 1e16;
        CoordinatesArrayType best_projected_point_on_slave_local;
        best_brep_id_slave = -1;
        bool isConvergedAtLeastOnce = false;
        CoordinatesArrayType best_projected_on_slave;
        for (IndexType i_brep_s = 0; i_brep_s < rPairedGeometryList.size(); i_brep_s++) {
            //---------------------------------------------------------------------------------------------------

            CoordinatesArrayType local_coord_projected_on_slave; //first trial
            local_coord_projected_on_slave[0] = 0;
            CoordinatesArrayType rProjectedPointGlobalCoordinates;
            double current_distance;
            

            Vector b_rep_interval;
            rPairedGeometryList[i_brep_s]->DomainInterval(b_rep_interval);

            bool isConverged = false;
            isConverged = NewtonRaphsonCurveOnDeformed(rpNurbsSurfaceParent, rpNurbsSurfacePaired,
                                        local_coord_projected_on_slave,
                                        rParentPointLocal,
                                        rProjectedPointGlobalCoordinates,
                                        rParentGeometry, *rPairedGeometryList[i_brep_s], 
                                        current_distance,
                                        rNumberInitialGuesses,
                                        rMaxIt, toll);

            // double temp_loc_coord = 0.0;
            // isConverged = ProjectionNurbsContactUtilities<TPointType, TSurfaceContainerPointType>::GetProjection(rpNurbsSurfaceParent, rpNurbsSurfacePaired, 
            //                                         rParentPointLocal, 
            //                                         rParentGeometry, *rPairedGeometryList[i_brep_s],
            //                                         temp_loc_coord, current_distance);
            // local_coord_projected_on_slave[0] = temp_loc_coord;

            if (isConverged) {
                isConvergedAtLeastOnce = true;

                if (current_distance < best_distance) {
                    best_distance = current_distance;
                    best_projected_point_on_slave_local = local_coord_projected_on_slave;

                    
                    best_brep_id_slave = i_brep_s;

                    // if (std::abs(local_coord_projected_on_slave[0] - b_rep_interval[0]) < 1e-12
                    //     || std::abs(local_coord_projected_on_slave[0] - b_rep_interval[1]) < 1e-12) {
                    //         best_projected_point_on_slave_local[0] += ((b_rep_interval[0]-local_coord_projected_on_slave[0]) + (b_rep_interval[1]-local_coord_projected_on_slave[0]))
                    //                                                     / (std::abs(b_rep_interval[1]-b_rep_interval[0])) *1e-9;
                    //     }
                }
            } else {
                Vector interval; 
                rPairedGeometryList[i_brep_s]->DomainInterval(interval);
                local_coord_projected_on_slave[0] = interval[0];
                CoordinatesArrayType point_projected_on_slave_parameter(3);        
                rPairedGeometryList[i_brep_s]->LocalCoordinates(point_projected_on_slave_parameter, local_coord_projected_on_slave);

                CoordinatesArrayType point_projected_on_slave_physical(3);        
                rPairedGeometryList[i_brep_s]->GlobalCoordinates(point_projected_on_slave_physical, local_coord_projected_on_slave);

                // NEW
                ProjectionNurbsContactUtilities<TPointType, TSurfaceContainerPointType>::GetDisplacement(rpNurbsSurfacePaired, point_projected_on_slave_parameter, displacement_on_projection);

                CoordinatesArrayType point_projected_on_slave_deformed = point_projected_on_slave_physical + displacement_on_projection;
                
                current_distance = norm_2(master_quadrature_point_deformed-point_projected_on_slave_deformed);

                if (current_distance < best_distance) {
                    best_distance = current_distance;
                    best_projected_point_on_slave_local = local_coord_projected_on_slave;

                    // best_projected_point_on_slave_local[0] -= 1e-9;
                    best_brep_id_slave = i_brep_s;
                }
                //*******************************************************
                // SECOND VERTEX
                //***************************************************** */ */
                local_coord_projected_on_slave[0] = interval[1];
                rPairedGeometryList[i_brep_s]->LocalCoordinates(point_projected_on_slave_parameter, local_coord_projected_on_slave);

                rPairedGeometryList[i_brep_s]->GlobalCoordinates(point_projected_on_slave_physical, local_coord_projected_on_slave);

                    // NEW
                ProjectionNurbsContactUtilities<TPointType, TSurfaceContainerPointType>::GetDisplacement(rpNurbsSurfacePaired, point_projected_on_slave_parameter, displacement_on_projection);
                point_projected_on_slave_deformed = point_projected_on_slave_physical + displacement_on_projection;

                current_distance = norm_2(master_quadrature_point_deformed-point_projected_on_slave_deformed);


                if (current_distance < best_distance) {
                    best_distance = current_distance;
                    best_projected_point_on_slave_local = local_coord_projected_on_slave;

                    // best_projected_point_on_slave_local[0] -= 1e-9;
                    best_brep_id_slave = i_brep_s;
                }
            }

        }

        ProjectionOnPairedGeometry = best_projected_point_on_slave_local;
        
        if (best_distance < 1e-1) isConvergedAtLeastOnce = true;
        return isConvergedAtLeastOnce;
    }

};
} // namespace Kratos

#endif // KRATOS_PROJECTION_NURBS_CONTACT_UTILITIES_H_INCLUDED
