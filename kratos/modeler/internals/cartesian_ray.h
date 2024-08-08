//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

#pragma once

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include  "utilities/intersection_utilities.h"


namespace Kratos
{
namespace Internals
{

///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/// A cartesian ray to be used in ray casting operations
/** This class represents a cartesian ray in 3D space defined by a direction and Point1 and Point2
*/
template<typename TGeometryType>
class CartesianRay
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CartesianRay
    KRATOS_CLASS_POINTER_DEFINITION(CartesianRay);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CartesianRay(): mIsValid(true), mDirection(0), mPoint1(), mPoint2() {}

    /// Constructor with all needed parameters
    CartesianRay(int Direction, Point const& Point1, Point const& Point2): mIsValid(true), mDirection(Direction), mPoint1(Point1), mPoint2(Point2) {}

    // Copy constructor
    CartesianRay(CartesianRay const& Other): mIsValid(Other.mIsValid), mDirection(Other.mDirection), mPoint1(Other.mPoint1), mPoint2(Other.mPoint2), mIntersections(Other.mIntersections) {}

    /// Destructor.
    virtual ~CartesianRay(){}

    ///@}
    ///@name Operators
    ///@{

    CartesianRay& operator=(CartesianRay const& Other){
        mIsValid = Other.mIsValid;
        mDirection = Other.mDirection;
        mPoint1 = Other.mPoint1;
        mPoint2 = Other.mPoint2;
        mIntersections = Other.mIntersections;

        return *this;
    }


    ///@}
    ///@name Operations
    ///@{

    void AddIntersection(TGeometryType const& rGeometry, double Tolerance){

        array_1d<double,3> intersection_point = ZeroVector(3);
        const int is_intersected = ComputeTriangleRayIntersection(
          rGeometry,
          mPoint1,
          mPoint2,
          intersection_point,
          Tolerance);

        if(is_intersected == 1){ // There is an intersection but not coplanar
            mIntersections.push_back(std::make_pair(intersection_point[mDirection], &rGeometry));
        }
    }

    void CollapseIntersectionPoints(double Tolerance){

        if (mIntersections.size() < 2) {
            return;
        }
        // Sort
        std::sort(mIntersections.begin(), mIntersections.end());
        auto i_unique_end = mIntersections.begin();
        auto i_patch_begin = mIntersections.begin();
        auto i_patch_end = mIntersections.begin();
        for(auto i_intersection = mIntersections.begin() + 1 ; i_intersection != mIntersections.end() ; i_intersection++) {
            if (std::abs(i_intersection->first - i_patch_begin->first) < Tolerance){ // we are in a patch so let's continue until the patch finishes
                i_patch_end = i_intersection;
                if(i_intersection + 1 != mIntersections.end()) { // if i_intersection is not the last one we keep continuing
                    continue;
                }
            }
            if((i_patch_begin == i_patch_end) || // The previous patch was only one intersection and we add it.
               (CheckPassingThroughByExtraRays(i_patch_begin, i_patch_end+1, Tolerance,  2.00*Tolerance))) { // more than one intersection to be checked with extra rays
                if(i_unique_end != i_patch_begin) {
                    *i_unique_end = std::move(*i_patch_begin);
                }
                i_unique_end++;
            }
            if((i_intersection != i_patch_end)&&(i_intersection + 1 == mIntersections.end())) { // Adding the last intersection
                if(i_unique_end != i_intersection) {
                    *i_unique_end = std::move(*i_intersection);
                }
                i_unique_end++;
            }
            i_patch_begin = i_intersection;
            i_patch_end = i_intersection;
        }
        auto new_size = std::distance(mIntersections.begin(), i_unique_end);
        mIntersections.resize(new_size);
    }

    void CalculateColor(std::vector<double> const& Coordinates, int InsideColor, int OutsideColor, std::vector<double>& ResultingColors, double NearEnough){

        bool is_inside=false;

        if(ResultingColors.size() != Coordinates.size()){
            ResultingColors.resize(Coordinates.size());
        }

        std::size_t current_index=0;
        const std::size_t size = Coordinates.size();

        for(auto& i_intersection : mIntersections){
            while(current_index < size){
                if((i_intersection.first - Coordinates[current_index]) > NearEnough){
                    ResultingColors[current_index++] = (is_inside) ? InsideColor : OutsideColor;
                } else if((i_intersection.first - Coordinates[current_index]) > -NearEnough){ // Considering very near to wall as inside
                    ResultingColors[current_index++] =  InsideColor;
                }
                else{
                    break;
                }
            }
            is_inside = !is_inside;
        }

        mIsValid = !is_inside; // Ray is not valid if is_inside is true after the last intersection

        // now coloring the points after the last intersection as outside
        while(current_index < size){
            ResultingColors[current_index++] = OutsideColor;
        }
    }

    void MarkIntersectedIntervals(std::vector<double> const& Coordinates, int InsideColor, int OutsideColor, std::vector<double>& ResultingColors, double NearEnough){

        const std::size_t size = Coordinates.size() - 1;

        if(ResultingColors.size() != Coordinates.size()){
            ResultingColors.resize(size);
        }

        std::size_t current_index=0;

        for(auto& i_intersection : mIntersections){
            while(current_index < size){
                if((i_intersection.first - Coordinates[current_index+1]) > NearEnough){
                    ResultingColors[current_index++] = OutsideColor;
                }
                else{
                    ResultingColors[current_index++] = InsideColor;
                    break;
                }
            }
        }

        // now coloring the points after the last intersection as outside
        while(current_index < size){
            ResultingColors[current_index++] = OutsideColor;
        }
    }


    ///@}
    ///@name Access
    ///@{

    std::vector<std::pair<double, const TGeometryType*>> const& GetIntersections() const {return mIntersections;}

    Point& GetPoint1(){return mPoint1;}

    Point& GetPoint2(){return mPoint2;}


    ///@}
    ///@name Inquiry
    ///@{

    bool IsValid(){return mIsValid;}

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "CartesianRay" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "CartesianRay";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

    ///@}
    ///@name Friends
    ///@{


    ///@}
private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    bool mIsValid;
    int mDirection;
    Point mPoint1;
    Point mPoint2;
    std::vector<std::pair<double, const TGeometryType*>> mIntersections;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    bool CheckPassingThroughByExtraRays(typename std::vector<std::pair<double, const TGeometryType*>>::iterator Begin, typename std::vector<std::pair<double, const TGeometryType*>>::iterator End, double Tolerance, double Delta){
        const std::array<double, 8> delta_u{Delta, Delta, 0.00, -Delta, -Delta, -Delta, 0.00, Delta};
        const std::array<double, 8> delta_v{0.00, Delta, Delta, Delta, 0.00, -Delta, -Delta, -Delta};
        const std::array<std::size_t, 4> axes{1,2,0,1};

        const std::size_t i_u = axes[mDirection];
        const std::size_t i_v = axes[mDirection+1];

        std::size_t no_hit_cases = 0;

        for(std::size_t i_ray = 0 ; i_ray < 8 ; i_ray++){
            CartesianRay extra_ray(mDirection, mPoint1, mPoint2);
            extra_ray.mPoint1[i_u] += delta_u[i_ray];
            extra_ray.mPoint1[i_v] += delta_v[i_ray];
            extra_ray.mPoint2[i_u] += delta_u[i_ray];
            extra_ray.mPoint2[i_v] += delta_v[i_ray];

           for(auto i_intersection = Begin ; i_intersection != End; i_intersection++){
                extra_ray.AddIntersection(*(i_intersection->second), Tolerance);
            }
            if(extra_ray.mIntersections.size() == 0){
                no_hit_cases++;
            }
            if(no_hit_cases > 4) // more than half
                return false;
        }
        return true;
    }

        /**
     * Find the 3D intersection of a ray with a triangle ignoring the coplanar and degenerated situations
     * Todo: To be passed to the intersection_utilities.h in core
     * @param rTriangleGeometry Is the triangle to intersect
     * @param rLinePoint1 Coordinates of the first point of the intersecting line
     * @param rLinePoint2 Coordinates of the second point of the intersecting line
     * @return rIntersectionPoint The intersection point coordinates
     * @return The intersection type index:
     * -1 (the triangle is degenerate)
     * 0 (disjoint - no intersection)
     * 1 (intersect in a unique point)
     * 2 (are in the same plane)
     */

    static int ComputeTriangleRayIntersection(
        const TGeometryType& rTriangleGeometry,
        const array_1d<double,3>& rLinePoint1,
        const array_1d<double,3>& rLinePoint2,
        array_1d<double,3>& rIntersectionPoint,
        const double epsilon = 1e-12) {

            // This is the adaption of the implemnetation of Moller-Trumbore algorithm provided in:
            // https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm

            const array_1d<double,3> edge1 = rTriangleGeometry[1] - rTriangleGeometry[0];
            const array_1d<double,3> edge2 = rTriangleGeometry[2] - rTriangleGeometry[0];
            const array_1d<double,3> line_vector = (rLinePoint2 - rLinePoint1)/norm_2(rLinePoint2 - rLinePoint1);
            array_1d<double,3> h, h_norm;
            const array_1d<double,3> edge1_normalized = edge1/norm_2(edge1);
            const array_1d<double,3> edge2_normalized = edge2/norm_2(edge2);
            MathUtils<double>::CrossProduct(h_norm, line_vector, edge2_normalized);
            MathUtils<double>::CrossProduct(h, line_vector, edge2);
            const double a_norm = inner_prod(edge1_normalized,h_norm);
            const double a = inner_prod(edge1,h);
            if (a_norm > -epsilon && a_norm < epsilon) {
                return 0;    // This ray is parallel to this triangle.
            }
            const double f = 1.0/a;
            const array_1d<double,3> s = rLinePoint1 - rTriangleGeometry[0];
            const double u = f * inner_prod(s,h);
            if (u < -epsilon || u > 1.0 + epsilon) {
                return 0;
            }
            array_1d<double,3> q;
            MathUtils<double>::CrossProduct(q, s, edge1);
            const double v = f * inner_prod(line_vector,q);
            if (v < -epsilon || u + v > 1.0 + epsilon) {
                return 0;
            }
            // At this stage we can compute t to find out where the intersection point is on the line.
            float t = f * inner_prod(edge2,q);
            if (t > epsilon && t < 1.00/epsilon) {
                // ray intersection
                rIntersectionPoint = rLinePoint1 + line_vector * t;
                return 1;
            }
            else {
                // This means that there is a line intersection but not a ray intersection.
                return 0;
            }
        }



    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class CartesianRay

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<typename TGeometryType>
inline std::istream& operator >> (std::istream& rIStream,
                CartesianRay<TGeometryType>& rThis){
                    return rIStream;
                }

/// output stream function
template<typename TGeometryType>
inline std::ostream& operator << (std::ostream& rOStream,
                const CartesianRay<TGeometryType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

} // namespace Internals

}  // namespace Kratos.

