//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

#if !defined(KRATOS_INTERNALS_CARTESIAN_RAY_H_INCLUDED )
#define  KRATOS_INTERNALS_CARTESIAN_RAY_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"


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
		const double relative_tolerance = 1.0e-12*std::sqrt(rGeometry.Length());
        const int is_intersected = IntersectionUtilities::ComputeTriangleRayIntersection(
          rGeometry,
          mPoint1,
          mPoint2,
          intersection_point,
          relative_tolerance);

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
        // Unique
        // auto i_unique = mIntersections.begin();
        // auto i_patch_begin = mIntersections.begin();
        // auto i_patch_end = mIntersections.begin();
        // bool patch_started = false;
        // for(auto i_intersection = mIntersections.begin() + 1 ; i_intersection != mIntersections.end() ; i_intersection++){
        //     if (std::abs(i_intersection->first - i_patch_begin->first) < Tolerance) { // There is a near hit so we start the patch
        //        if(i_intersection + 1 == mIntersections.end()){ // we are at the end so we should finish the patch and check if it is really pass through a duplicated surface or just passing tangent to model
        //             i_patch_end = i_intersection;
        //             if(CheckPassingThroughByExtraRays(i_patch_begin, i_patch_end, Tolerance,  2.00*Tolerance)) {
        //                 if(i_unique != i_patch_begin) {
        //                     *i_unique=std::move(*i_patch_begin); 
        //                 }
        //                 i_unique++;
        //             }
        //         }
        //         else {
        //             patch_started = true;
        //         }
        //     }
        //     else { // the hit points are far enough 
        //         if(patch_started){ // then we should finish the patch and check if it is really pass through a duplicated surface or just passing tangent to model
        //             i_patch_end = i_intersection;
        //             if(CheckPassingThroughByExtraRays(i_patch_begin, i_patch_end, Tolerance,  2.00*Tolerance)) {
        //                 if(i_unique != i_patch_begin) {
        //                     *i_unique=std::move(*i_patch_begin); 
        //                 }
        //                 i_unique++;
        //             }
        //             patch_started = false;
        //         }
        //         else{
        //             i_patch_begin = i_intersection; // This would be the next patch begin or not depending on the next intersection.
        //             if(i_unique != i_intersection) {
        //                 *i_unique=std::move(*i_intersection);
        //             } 
        //             i_unique++;
        //         }
        //     }
        // }
        // auto new_size = std::distance(mIntersections.begin(), i_unique);
        std::size_t new_size = 0;
        auto i_begin = mIntersections.begin();
        auto i_intersection = mIntersections.begin();
        while (++i_begin != mIntersections.end()) {
            // considering the very near points as the same points
            if (std::abs(i_begin->first - i_intersection->first) > Tolerance) { // if the hit points are far enough they are not the same 
                if(new_size == 0)
                    new_size++; // we should consider the first intersection then.    
                if(++i_intersection != i_begin)           
                    *i_intersection = std::move(*i_begin);
                new_size++;
            }
            else{ // Now there are near hits, so we check if it is really pass through a duplicated surface or just passing tangent to model
                // Getting the patch of geometries with near hit
                auto i_patch_begin = i_intersection;
                auto i_patch_end = i_begin;
                while(++i_begin != mIntersections.end()){
                    if (std::abs(i_begin->first - i_intersection->first) > Tolerance) { // This is the end of the patch.               
                        break;
                    }
                }
                
                i_patch_end = i_begin;

                if(CheckPassingThroughByExtraRays(i_patch_begin, i_patch_end, Tolerance,  2.00*Tolerance)) {
                    if(++i_intersection != i_begin)           
                        *i_intersection = std::move(*i_begin);
                    new_size++;
                }

                if(i_begin == mIntersections.end())
                    break;
            }
        }
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
            is_inside = (is_inside) ? false:true;
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
                    ResultingColors[current_index] = OutsideColor;
                    current_index++;
                }
                else{
                    ResultingColors[current_index] = InsideColor;
                    current_index++;
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
        const std::array<int, 4> axes{1,2,0,1};

        const int i_u = axes[mDirection];
        const int i_v = axes[mDirection+1];

        std::size_t no_hit_cases = 0;

        for(int i_ray = 0 ; i_ray < 8 ; i_ray++){
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

#endif // KRATOS_INTERNALS_CARTESIAN_RAY_H_INCLUDED  defined


