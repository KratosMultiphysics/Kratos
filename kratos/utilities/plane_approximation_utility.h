//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#if !defined(KRATOS_PLANE_APPROXIMATION_UTILITY )
#define  KRATOS_PLANE_APPROXIMATION_UTILITY

/* System includes */

/* External includes */

/* Project includes */
#include "utilities/math_utils.h"

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
 * @class PlaneApproximationUtility
 * @ingroup KratosCore
 * @brief Utility to compute an approximation plane for a set of points
 * @details Provided a set of points, this utility computes the plane that best
 * approximate them in a least square sense. The base point of the plane is the
 * set of points average coordinates while the plane normal is taken as the
 * eigenvector associated to the minimum eigenvalue of the matrix A. Matrix A is
 * trans(P)*P, being P the matrix that stores in each row the difference between 
 * the base point and each one of the points to approximate.
 * @author Ruben Zorrilla
 */
class PlaneApproximationUtility
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of PlaneApproximationUtility
    KRATOS_CLASS_POINTER_DEFINITION( PlaneApproximationUtility );
    
    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief Default constructor
     */
    PlaneApproximationUtility(){}

    /// Destructor
    virtual ~PlaneApproximationUtility(){}

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Provided a list of point coordinates, computes 
     * the best plane approximation in a least squares sense.
     * @param rPointsCoords Vector containing the set of point coordinates
     * @return rPlaneBasePointCoords Plane base point
     * @return rPlaneNormal Plane normal
     */
    static void ComputePlaneApproximation(
        const std::vector< array_1d<double,3> > &rPointsCoords,
        array_1d<double,3> &rPlaneBasePointCoords,
        array_1d<double,3> &rPlaneNormal)
    {
        rPlaneBasePointCoords = GetPlaneBasePoint(rPointsCoords);
        rPlaneNormal = GetPlaneNormal(rPointsCoords, rPlaneBasePointCoords);
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

private:

    ///@name Private static Member Variables
    ///@{

    ///@}
    ///@name Private member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This function computes the plane base point
     * @param rPointsCoords Vector containing the set of point coordinates
     * @return base_point Plane base point
     */
    static array_1d<double, 3> GetPlaneBasePoint(
        const std::vector< array_1d<double,3> > &rPointsCoords) 
    {
        const unsigned int n_points = rPointsCoords.size();
        KRATOS_ERROR_IF(n_points == 0) << "No base point can be computed. Points container is empty." << std::endl;

        array_1d<double, 3> base_point = ZeroVector(3);
        for (unsigned int j = 0; j < n_points; ++j){
            base_point += rPointsCoords[j];
        }
        base_point /= n_points;

        return base_point;
    }

    /**
     * @brief This function sets the matrix A as trans(P)*P
     * @param rPointsCoords Vector containing the set of point coordinates
     * @return rA A matrix
     */
    static bounded_matrix<double,3,3> SetMatrixA(
        const std::vector< array_1d< double,3 > > &rPointsCoords,
        const array_1d<double,3> &rPlaneBasePointCoords) 
    {
        bounded_matrix<double,3,3> A = ZeroMatrix(3);
        const unsigned int n_points = rPointsCoords.size();

        for (unsigned int i = 0; i < 3; ++i){
            const double base_i = rPlaneBasePointCoords(i);
            for (unsigned int j = 0; j < 3; ++j){
                const double base_j = rPlaneBasePointCoords(j);
                for (unsigned int m = 0; m < n_points; ++m){
                    const double pt_m_i = rPointsCoords[m](i);
                    const double pt_m_j = rPointsCoords[m](j);
                    A(i,j) += (base_i - pt_m_i)*(base_j - pt_m_j);
                }
            }
        }

        return A;
    }

    /**
     * @brief This function computes using the inverse power method the minimal eigenvalue
     * @param InputMatrix The matrix to compute the eigenvalue
     * @return condition_number The condition number
     */
    static array_1d<double, 3> GetPlaneNormal(
        const std::vector< array_1d<double,3> > &rPointsCoords,
        const array_1d<double,3> &rPlaneBasePointCoords) 
    {
        // Solve the A matrix eigenvalue problem 
        bounded_matrix<double, 3, 3> eigenval_mat, eigenvector_mat;
        bool converged = MathUtils<double>::EigenSystem<3>(SetMatrixA(rPointsCoords, rPlaneBasePointCoords), eigenvector_mat, eigenval_mat);
        KRATOS_ERROR_IF(!converged) << "Plane normal can't be computed. Eigenvalue problem did not converge." << std::endl;

        // Find the minimum eigenvalue
        double min_eigval = 0.0;
        unsigned int min_eigval_id = 0;
        for (unsigned int i = 0; i < 3; ++i){
            if (eigenval_mat(i,i) < min_eigval){
                min_eigval_id = i;
                min_eigval = std::min(min_eigval, eigenval_mat(i,i));
            }
        }

        // Set as plane normal the eigenvector associated to the minimum eigenvalue
        array_1d<double,3> plane_normal;
        for (unsigned int i = 0; i < 3; ++i){
            plane_normal(i) = eigenvector_mat(min_eigval_id,i);
        }

        return plane_normal;
    }

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Private LifeCycle
    ///@{

    ///@}
    ///@name Unaccessible methods
    ///@{

}; /* Class PlaneApproximationUtility */

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

}  /* namespace Kratos.*/

#endif /* KRATOS_PLANE_APPROXIMATION_UTILITY  defined */

