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
#include <algorithm>

// External includes

// Project includes
#include "processes/process.h"
#include "includes/model_part.h"
#include "modeler/internals/cartesian_ray.h"
#include "utilities/parallel_utilities.h"

namespace Kratos
{
  namespace Internals
  {
    class CartesianMeshColors{
      double mTolerance;
      Point mMinPoint;
      Point mMaxPoint;
      array_1d<std::vector<double>, 3> mNodalCoordinates;
      array_1d<std::vector<double>, 3> mElementCenterCoordinates;
      std::vector<array_1d<double,3>> mNodalRayColors;
      std::vector<array_1d<double,3>> mElementalRayColors;
      std::vector<double> mNodalColors;
      std::vector<double> mElementalColors;
      std::vector<array_1d<double,6>> mElementalFaceColors;
      DenseMatrix<Internals::CartesianRay<Element::GeometryType>> mXYRays;
      DenseMatrix<Internals::CartesianRay<Element::GeometryType>> mXZRays;
      DenseMatrix<Internals::CartesianRay<Element::GeometryType>> mYZRays;

      inline std::size_t ElementIndex(std::size_t I, std::size_t J, std::size_t K) const {
        return I + J * mElementCenterCoordinates[0].size() + K * mElementCenterCoordinates[1].size() * mElementCenterCoordinates[0].size();
      }
     public:
      CartesianMeshColors(): mTolerance(1e-12), mMinPoint(0.00, 0.00, 0.00), mMaxPoint(0.00, 0.00, 0.00){}

      std::vector<double> const& GetNodalCoordinates(int Index) const {return mNodalCoordinates[Index];}

      std::vector<double> const& GetElementCenterCoordinates(int Index) const {return mElementCenterCoordinates[Index];}

    template<typename TPointsContainerType>
      void ExtendBoundingBox(TPointsContainerType const& Points, double Margine){
        for(auto const& point : Points){
            for(std::size_t i = 0; i<3; i++)
            {
                mMinPoint[i] =  (mMinPoint[i] >  point[i] - Margine) ?  point[i] - Margine : mMinPoint[i];
                mMaxPoint[i] =  (mMaxPoint[i] <  point[i] + Margine) ?  point[i] + Margine : mMaxPoint[i];
            }
        }
      }

    void SetCoordinates(std::vector<double> const& rNodalXCoordinates, std::vector<double> const& rNodalYCoordinates, std::vector<double> const& rNodalZCoordinates){

        KRATOS_ERROR_IF((rNodalXCoordinates.size() < 2) || (rNodalYCoordinates.size() < 2) || (rNodalZCoordinates.size() < 2))
            << "The coordinates should have at least two entries defining the bounding box." << std::endl;

        mNodalCoordinates[0] = rNodalXCoordinates;
        mNodalCoordinates[1] = rNodalYCoordinates;
        mNodalCoordinates[2] = rNodalZCoordinates;

        for(std::size_t i = 0 ; i < 3 ; i++){
            for(std::size_t j = 0 ; j < mNodalCoordinates[i].size() - 1 ; j++){
                double center_coordinate = (mNodalCoordinates[i][j] + mNodalCoordinates[i][j+1]) * .5;
                mElementCenterCoordinates[i].push_back(center_coordinate);
            }
        }

        mXYRays.resize(mNodalCoordinates[0].size(),mNodalCoordinates[1].size(),false);
        mXZRays.resize(mNodalCoordinates[0].size(),mNodalCoordinates[2].size(),false);
        mYZRays.resize(mNodalCoordinates[1].size(),mNodalCoordinates[2].size(),false);

        mNodalRayColors.resize(mNodalCoordinates[0].size()*mNodalCoordinates[1].size()*mNodalCoordinates[2].size());
        mNodalColors.resize(mNodalCoordinates[0].size()*mNodalCoordinates[1].size()*mNodalCoordinates[2].size());
        mElementalRayColors.resize((mNodalCoordinates[0].size() - 1) * (mNodalCoordinates[1].size() - 1)*(mNodalCoordinates[2].size() - 1));
        mElementalColors.resize((mNodalCoordinates[0].size() - 1) * (mNodalCoordinates[1].size() - 1)*(mNodalCoordinates[2].size() - 1));
        mElementalFaceColors.resize((mNodalCoordinates[0].size() - 1) * (mNodalCoordinates[1].size() - 1)*(mNodalCoordinates[2].size() - 1));

        Point min_point( mNodalCoordinates[0][0] - mTolerance, mNodalCoordinates[1][0] - mTolerance, mNodalCoordinates[2][0] - mTolerance);
        Point max_point( mNodalCoordinates[0].back() + mTolerance, mNodalCoordinates[1].back() + mTolerance, mNodalCoordinates[2].back() + mTolerance);

        mMinPoint = min_point;
        mMaxPoint = max_point;

        SetAllColors(0.00);
     }

    void SetAllColors(double TheColor){

        const std::size_t number_of_nodes = mNodalColors.size();
        for(std::size_t i = 0 ; i < number_of_nodes ; i++ ){
            mNodalRayColors[i] = ScalarVector(3,TheColor);
            mNodalColors[i] = TheColor;
        }

        const std::size_t number_of_elements = mElementalColors.size();
        for(std::size_t i = 0 ; i < number_of_elements ; i++ ){
            mElementalRayColors[i] = ScalarVector(3,TheColor);
            mElementalColors[i] = TheColor;
            mElementalFaceColors[i] = ScalarVector(6,TheColor);
        }
    }

    double& GetNodalColor(std::size_t I, std::size_t J, std::size_t K){
        const std::size_t index = I + J * mNodalCoordinates[0].size() + K * mNodalCoordinates[1].size() * mNodalCoordinates[0].size();
        return mNodalColors[index];
    }

    double& GetElementalColor(std::size_t I, std::size_t J, std::size_t K){
        return mElementalColors[ElementIndex(I, J, K)];
    }

    double GetElementalColor(std::size_t I, std::size_t J, std::size_t K) const {
        return mElementalColors[ElementIndex(I, J, K)];
    }

    array_1d<double,6>& GetElementalFaceColor(std::size_t I, std::size_t J, std::size_t K){
        return mElementalFaceColors[ElementIndex(I, J, K)];
    }

    const array_1d<double,6>& GetElementalFaceColor(std::size_t I, std::size_t J, std::size_t K) const {
        return mElementalFaceColors[ElementIndex(I, J, K)];
    }

    std::vector<double>& GetNodalColors(){ return mNodalColors;}

    std::vector<double>& GetElementalColors(){return mElementalColors;}

    std::vector<array_1d<double,6>>& GetElementalFaceColors(){return mElementalFaceColors;}

    Point GetPoint(std::size_t I, std::size_t J, std::size_t K) const {
        return Point(mNodalCoordinates[0][I], mNodalCoordinates[1][J], mNodalCoordinates[2][K]);
    }

    Point GetPoint(double I, double J, double K) const {
        // The coordinates of the point are interpolated between those of floor(I, J, K) and ceil(I, J, K)
        std::size_t i(I), j(J), k(K);
        std::vector<double> start{mNodalCoordinates[0][i], mNodalCoordinates[1][j], mNodalCoordinates[2][k]};
        return Point(
            start[0] + (I-i)*(mNodalCoordinates[0][i+1]-start[0]),
            start[1] + (J-j)*(mNodalCoordinates[1][j+1]-start[1]),
            start[2] + (K-k)*(mNodalCoordinates[2][k+1]-start[2]));
    }

    Point GetCenterOfElement(std::size_t I, std::size_t J, std::size_t K) const {
        return Point(mElementCenterCoordinates[0][I], mElementCenterCoordinates[1][J], mElementCenterCoordinates[2][K]);
    }

    CartesianRay<Element::GeometryType>& GetXYRay(std::size_t I, std::size_t J){
        return mXYRays(I,J);
    }

    array_1d<std::size_t,3> GetGridSize() const {
        array_1d<std::size_t,3> tmp{mNodalCoordinates[0].size(), mNodalCoordinates[1].size(), mNodalCoordinates[2].size()};
        return tmp;
    }

    void InitializeRays(array_1d< std::size_t, 3 > const& MinRayPosition, array_1d< std::size_t, 3 > const& MaxRayPosition, const std::string& EntititesToColor){

        std::vector<double>* p_x_coordinates = &(mNodalCoordinates[0]);
        std::vector<double>* p_y_coordinates = &(mNodalCoordinates[1]);
        std::vector<double>* p_z_coordinates = &(mNodalCoordinates[2]);

        if((EntititesToColor == "center_of_elements") || (EntititesToColor == "face_of_elements")){
            p_x_coordinates = &(mElementCenterCoordinates[0]);
            p_y_coordinates = &(mElementCenterCoordinates[1]);
            p_z_coordinates = &(mElementCenterCoordinates[2]);
        }


        IndexPartition<std::size_t>(MaxRayPosition[0]-MinRayPosition[0]).for_each([&](const std::size_t p) {
            std::size_t i = MinRayPosition[0] + p;
            for(std::size_t j = MinRayPosition[1] ; j < MaxRayPosition[1] ; j++){
                mXYRays(i,j) = Internals::CartesianRay<Element::GeometryType>(2,
                        Point((*p_x_coordinates)[i], (*p_y_coordinates)[j], mMinPoint[2]),
                        Point((*p_x_coordinates)[i], (*p_y_coordinates)[j], mMaxPoint[2]));
            }
        });
        IndexPartition<std::size_t>(MaxRayPosition[0]-MinRayPosition[0]).for_each([&](const std::size_t p) {
            std::size_t i = MinRayPosition[0] + p;
            for(std::size_t k = MinRayPosition[2] ; k < MaxRayPosition[2] ; k++){
                mXZRays(i,k) = Internals::CartesianRay<Element::GeometryType>(1,
                        Point((*p_x_coordinates)[i], mMinPoint[1], (*p_z_coordinates)[k]),
                        Point((*p_x_coordinates)[i], mMaxPoint[1], (*p_z_coordinates)[k]));
            }
        });
        IndexPartition<std::size_t>(MaxRayPosition[1]-MinRayPosition[1]).for_each([&](const std::size_t p) {
            std::size_t j = MinRayPosition[1] + p;
            for(std::size_t k = MinRayPosition[2] ; k < MaxRayPosition[2] ; k++){
                mYZRays(j,k) = Internals::CartesianRay<Element::GeometryType>(0,
                        Point(mMinPoint[0], (*p_y_coordinates)[j], (*p_z_coordinates)[k]),
                        Point(mMaxPoint[0], (*p_y_coordinates)[j], (*p_z_coordinates)[k]));
            }
        });
    }

        void AddGeometry(Element::GeometryType const& rGeometry, bool IsNodal){
                array_1d< std::size_t, 3 > min_position;
                array_1d< std::size_t, 3 > max_position;
                if(IsNodal) {
                    CalculateMinMaxNodePositions(rGeometry, min_position, max_position);
                }
                else
                {
                    CalculateMinMaxCenterOfElementPositions(rGeometry, min_position, max_position);
                }

                // #pragma omp parallel for
                for(std::size_t i = min_position[0] ; i < max_position[0] ; i++){
                    for(std::size_t j = min_position[1] ; j < max_position[1] ; j++){
                        mXYRays(i,j).AddIntersection(rGeometry, mTolerance);
                    }
                }

                // #pragma omp parallel for
                for(std::size_t i = min_position[0] ; i < max_position[0] ; i++){
                    for(std::size_t k = min_position[2] ; k < max_position[2] ; k++){
                        mXZRays(i,k).AddIntersection(rGeometry, mTolerance);
                    }
                }

                // #pragma omp parallel for
                for(std::size_t j = min_position[1] ; j < max_position[1] ; j++){
                    for(std::size_t k = min_position[2] ; k < max_position[2] ; k++){
                        mYZRays(j,k).AddIntersection(rGeometry, mTolerance);
                    }
                }

        }

        void CalculateNodalRayColors(array_1d< std::size_t, 3 > const& MinRayPosition, array_1d< std::size_t, 3 > const& MaxRayPosition, int InsideColor, int OutsideColor){
            std::vector<double> colors;
            // #pragma omp parallel for
            for(std::size_t i = MinRayPosition[0] ; i < MaxRayPosition[0] ; i++){
                for(std::size_t j = MinRayPosition[1] ; j < MaxRayPosition[1] ; j++){
                    auto& ray = mXYRays(i,j);
                    ray.CollapseIntersectionPoints(mTolerance);
                    ray.CalculateColor(mNodalCoordinates[2], InsideColor, OutsideColor, colors, mTolerance);
                    for(std::size_t k = MinRayPosition[2] ; k < MaxRayPosition[2] ; k++){
                        if(colors[k] == InsideColor){
                            GetNodalRayColor(i,j,k)[0] = colors[k];
                        }
                    }
                }
            }
            // #pragma omp parallel for
            for(std::size_t i = MinRayPosition[0] ; i < MaxRayPosition[0] ; i++){
                for(std::size_t k = MinRayPosition[2] ; k < MaxRayPosition[2] ; k++){
                    auto& ray = mXZRays(i,k);
                    ray.CollapseIntersectionPoints(mTolerance);
                    ray.CalculateColor(mNodalCoordinates[1], InsideColor, OutsideColor, colors, mTolerance);
                    for(std::size_t j = MinRayPosition[1] ; j < MaxRayPosition[1] ; j++){
                        if(colors[j] == InsideColor){
                            GetNodalRayColor(i,j,k)[1] = colors[j];
                        }
                    }
                }
            }
            // #pragma omp parallel for
            for(std::size_t j = MinRayPosition[1] ; j < MaxRayPosition[1] ; j++){
                for(std::size_t k = MinRayPosition[2] ; k < MaxRayPosition[2] ; k++){
                    mYZRays(j,k).CollapseIntersectionPoints(mTolerance);
                    mYZRays(j,k).CalculateColor(mNodalCoordinates[0], InsideColor, OutsideColor, colors, mTolerance);
                    for(std::size_t i = MinRayPosition[0] ; i < MaxRayPosition[0] ; i++){
                        if(colors[i] == InsideColor){
                            GetNodalRayColor(i,j,k)[2] = colors[i];
                        }
                    }
                }
            }

            // #pragma omp parallel for
            for(std::size_t k = MinRayPosition[2] ; k < MaxRayPosition[2] ; k++){
                for(std::size_t j = MinRayPosition[1] ; j < MaxRayPosition[1] ; j++){
                    for(std::size_t i = MinRayPosition[0] ; i < MaxRayPosition[0] ; i++){
                        auto& ray_colors = GetNodalRayColor(i,j,k);
                        std::size_t n_inside = 0;
                        std::size_t n_outside = 0;
                        for(std::size_t dim = 0 ; dim < 3 ; dim++){
                            if(ray_colors[dim] == InsideColor){
                                n_inside++;
                            }
                            else {
                                n_outside++;
                            }
                        }
                        GetNodalColor(i,j,k) = (n_inside > n_outside) ? InsideColor : OutsideColor;
                    }
                }
            }
       }

        void CalculateElementalRayColors(array_1d< std::size_t, 3 > const& MinRayPosition, array_1d< std::size_t, 3 > const& MaxRayPosition, int InsideColor, int OutsideColor){
            std::vector<double> colors;
            std::vector<double> x_coordinates(mNodalCoordinates[0].size() - 1);
            std::vector<double> y_coordinates(mNodalCoordinates[1].size() - 1);
            std::vector<double> z_coordinates(mNodalCoordinates[2].size() - 1);

            for(std::size_t i = 0 ; i < x_coordinates.size() ; i++){
                x_coordinates[i] = (mNodalCoordinates[0][i] +  mNodalCoordinates[0][i+1]) * 0.5;
            }

            for(std::size_t i = 0 ; i < y_coordinates.size() ; i++){
                y_coordinates[i] = (mNodalCoordinates[1][i] +  mNodalCoordinates[1][i+1]) * 0.5;
            }

            for(std::size_t i = 0 ; i < z_coordinates.size() ; i++){
                z_coordinates[i] = (mNodalCoordinates[2][i] +  mNodalCoordinates[2][i+1]) * 0.5;
            }

            // #pragma omp parallel for
            for(std::size_t i = MinRayPosition[0] ; i < MaxRayPosition[0] ; i++){
                for(std::size_t j = MinRayPosition[1] ; j < MaxRayPosition[1] ; j++){
                    auto& ray = mXYRays(i,j);
                    ray.CollapseIntersectionPoints(mTolerance);
                    ray.CalculateColor(z_coordinates, InsideColor, OutsideColor, colors, mTolerance);
                    for(std::size_t k = MinRayPosition[2] ; k < MaxRayPosition[2] ; k++){
                        if(colors[k] == InsideColor){
                            GetElementalRayColor(i,j,k)[0] = colors[k];
                        }
                    }
                }
            }
            // #pragma omp parallel for
            for(std::size_t i = MinRayPosition[0] ; i < MaxRayPosition[0] ; i++){
                for(std::size_t k = MinRayPosition[2] ; k < MaxRayPosition[2] ; k++){
                    auto& ray = mXZRays(i,k);
                    ray.CollapseIntersectionPoints(mTolerance);
                    ray.CalculateColor(y_coordinates, InsideColor, OutsideColor, colors, mTolerance);
                    for(std::size_t j = MinRayPosition[1] ; j < MaxRayPosition[1] ; j++){
                        if(colors[j] == InsideColor){
                            GetElementalRayColor(i,j,k)[1] = colors[j];
                        }
                    }
                }
            }
            // #pragma omp parallel for
            for(std::size_t j = MinRayPosition[1] ; j < MaxRayPosition[1] ; j++){
                for(std::size_t k = MinRayPosition[2] ; k < MaxRayPosition[2] ; k++){
                    mYZRays(j,k).CollapseIntersectionPoints(mTolerance);
                    mYZRays(j,k).CalculateColor(x_coordinates, InsideColor, OutsideColor, colors, mTolerance);
                    for(std::size_t i = MinRayPosition[0] ; i < MaxRayPosition[0] ; i++){
                        if(colors[i] == InsideColor){
                            GetElementalRayColor(i,j,k)[2] = colors[i];
                        }
                    }
                }
            }
            // #pragma omp parallel for
            for(std::size_t k = MinRayPosition[2] ; k < MaxRayPosition[2] ; k++){
                for(std::size_t j = MinRayPosition[1] ; j < MaxRayPosition[1] ; j++){
                    for(std::size_t i = MinRayPosition[0] ; i < MaxRayPosition[0] ; i++){
                        auto& ray_colors = GetElementalRayColor(i,j,k);
                        std::size_t n_inside = 0;
                        std::size_t n_outside = 0;
                        for(std::size_t dim = 0 ; dim < 3 ; dim++){
                            if(ray_colors[dim] == InsideColor){
                                n_inside++;
                            }
                            else {
                                n_outside++;
                            }
                        }
                        auto& color = GetElementalColor(i,j,k);
                        color = (n_inside > n_outside) ? InsideColor : color;
                    }
                }
            }

        }

        void CalculateElementalFaceColors(array_1d< std::size_t, 3 > const& MinRayPosition, array_1d< std::size_t, 3 > const& MaxRayPosition, int FaceColor, int OutsideColor, int VolumeColor){
            if(FaceColor == OutsideColor)
                return; // Nothing to do!


            const std::size_t size_x = mElementCenterCoordinates[0].size();
            const std::size_t size_y = mElementCenterCoordinates[1].size();
            const std::size_t size_z = mElementCenterCoordinates[2].size();

            std::vector<double> x_coordinates(size_x + 2);
            std::vector<double> y_coordinates(size_y + 2);
            std::vector<double> z_coordinates(size_z + 2);

            x_coordinates.front() = mMinPoint[0];
            y_coordinates.front() = mMinPoint[1];
            z_coordinates.front() = mMinPoint[2];

            for(std::size_t i = 0 ; i < size_x ; i++){
                x_coordinates[i+1] = mElementCenterCoordinates[0][i];
            }

            for(std::size_t i = 0 ; i < size_y ; i++){
                y_coordinates[i+1] = mElementCenterCoordinates[1][i];
            }

            for(std::size_t i = 0 ; i < size_z ; i++){
                z_coordinates[i+1] = mElementCenterCoordinates[2][i];
            }

            x_coordinates.back() = mMaxPoint[0];
            y_coordinates.back() = mMaxPoint[1];
            z_coordinates.back() = mMaxPoint[2];

            IndexPartition<std::size_t>(MaxRayPosition[0]-MinRayPosition[0]).for_each(
                std::vector<double>(), [&](const std::size_t p, std::vector<double>& rColors) {
                    std::size_t i = MinRayPosition[0] + p;
                    for(std::size_t j = MinRayPosition[1] ; j < MaxRayPosition[1] ; j++){
                        auto& ray = mXYRays(i,j);
                        ray.CollapseIntersectionPoints(mTolerance);
                        ray.MarkIntersectedIntervals(z_coordinates, FaceColor, OutsideColor, rColors, mTolerance);

                        double previous_center_color = OutsideColor;
                        for(std::size_t k = 0 ; k < size_z; k++){
                            auto next_center_color = GetElementalColor(i,j,k);
                            auto interval_color = rColors[k];
                            if(interval_color == FaceColor){
                                if((previous_center_color == VolumeColor) && (next_center_color != VolumeColor)){
                                    GetElementalFaceColor(i,j,k-1)[5] = FaceColor;   // [-x,-y,-z,x,y,z]
                                }
                                else if((previous_center_color != VolumeColor) && (next_center_color == VolumeColor)){
                                    GetElementalFaceColor(i,j,k)[2] = FaceColor;   // [-x,-y,-z,x,y,z]
                                }
                                else if((k > 1) && (rColors[k-2] == VolumeColor)){
                                    GetElementalFaceColor(i,j,k-2)[5] = FaceColor;   // [-x,-y,-z,x,y,z]
                                }
                                else if((k < size_z - 1) && (rColors[k+1] == VolumeColor)){
                                    GetElementalFaceColor(i,j,k+1)[2] = FaceColor;   // [-x,-y,-z,x,y,z]
                                }
                                else{
                                    GetElementalFaceColor(i,j,k)[2] = -FaceColor;  // Error color is the -FaceColor
                                    // KRATOS_WARNING("CartesianMesh") << "The given interface is not in the interface of the volume for cell [" << i << "," << j << "," << k << "]" << std::endl;
                                }
                            }
                            previous_center_color = next_center_color;
                        }
                        if(rColors.back() == FaceColor){
                            if(previous_center_color == VolumeColor){
                                GetElementalFaceColor(i,j,size_z-1)[5] = FaceColor;   // [-x,-y,-z,x,y,z]
                            }
                            else{
                                GetElementalFaceColor(i,j,size_z-1)[5] = -FaceColor;   // Error color is the -FaceColor
                                // KRATOS_WARNING("CartesianMesh") << "The given interface is not in the interface of the volume for cell [" << i << "," << j << "," << size_z-1 << "]" << std::endl;
                            }
                        }
                    }
                }
            );

            IndexPartition<std::size_t>(MaxRayPosition[0]-MinRayPosition[0]).for_each(
                std::vector<double>(), [&](const std::size_t p, std::vector<double>& rColors) {
                    std::size_t i = MinRayPosition[0] + p;
                    for(std::size_t k = MinRayPosition[2] ; k < MaxRayPosition[2] ; k++){
                        auto& ray = mXZRays(i,k);
                        ray.CollapseIntersectionPoints(mTolerance);
                        ray.MarkIntersectedIntervals(y_coordinates, FaceColor, OutsideColor, rColors, mTolerance);
                        double previous_center_color = OutsideColor;
                        for(std::size_t j = 0 ; j < size_y ; j++){
                            auto next_center_color = GetElementalColor(i,j,k);
                            auto interval_color = rColors[j];
                            if(interval_color == FaceColor){
                                if((previous_center_color == VolumeColor) && (next_center_color != VolumeColor)){
                                    GetElementalFaceColor(i,j-1,k)[4] = FaceColor;   // [-x,-y,-z,x,y,z]
                                }
                                else if((previous_center_color != VolumeColor) && (next_center_color == VolumeColor)){
                                    GetElementalFaceColor(i,j,k)[1] = FaceColor;   // [-x,-y,-z,x,y,z]
                                }
                                else if((j > 1) && (rColors[j-2] == VolumeColor)){
                                    GetElementalFaceColor(i,j-2,k)[4] = FaceColor;   // [-x,-y,-z,x,y,z]
                                }
                                else if((j < size_y - 1) && (rColors[j+1] == VolumeColor)){
                                    GetElementalFaceColor(i,j+1,k)[1] = FaceColor;   // [-x,-y,-z,x,y,z]
                                }
                                else{
                                    GetElementalFaceColor(i,j,k)[1] = -FaceColor;   // Error color is the -FaceColor
                                    // KRATOS_WARNING("CartesianMesh") << "The given interface is not in the interface of the volume for cell [" << i << "," << j << "," << k << "]" << std::endl;
                                }
                            }
                            previous_center_color = next_center_color;
                        }
                        if(rColors.back() == FaceColor){
                            if(previous_center_color == VolumeColor){
                                GetElementalFaceColor(i, size_y-1,k)[4] = FaceColor;   // [-x,-y,-z,x,y,z]
                            }
                            else{
                                GetElementalFaceColor(i, size_y-1,k)[4] = -FaceColor;   // Error color is the -FaceColor
                                // KRATOS_WARNING("CartesianMesh") << "The given interface is not in the interface of the volume for cell [" << i << "," << size_y-1 << "," << k << "]" << std::endl;
                            }
                        }
                    }
                }
            );

            IndexPartition<std::size_t>(MaxRayPosition[1]-MinRayPosition[1]).for_each(
                std::vector<double>(), [&](const std::size_t p, std::vector<double>& rColors) {
                    std::size_t j = MinRayPosition[1] + p;
                    for(std::size_t k = MinRayPosition[2] ; k < MaxRayPosition[2] ; k++){
                        auto& ray= mYZRays(j,k);
                        ray.CollapseIntersectionPoints(mTolerance);
                        ray.MarkIntersectedIntervals(x_coordinates, FaceColor, OutsideColor, rColors, mTolerance);
                        double previous_center_color = OutsideColor;
                        for(std::size_t i = 0 ; i < size_x ; i++){
                            auto next_center_color = GetElementalColor(i,j,k);
                            auto interval_color = rColors[i];
                            if(interval_color == FaceColor){
                                if((previous_center_color == VolumeColor) && (next_center_color != VolumeColor)){
                                    GetElementalFaceColor(i-1,j,k)[3] = FaceColor;   // [-x,-y,-z,x,y,z]
                                }
                                else if((previous_center_color != VolumeColor) && (next_center_color == VolumeColor)){
                                    GetElementalFaceColor(i,j,k)[0] = FaceColor;   // [-x,-y,-z,x,y,z]
                                }
                                else if((i > 1) && (rColors[i-2] == VolumeColor)){
                                    GetElementalFaceColor(i-2,j,k)[3] = FaceColor;   // [-x,-y,-z,x,y,z]
                                }
                                else if((i < size_x - 1) && (rColors[i+1] == VolumeColor)){
                                    GetElementalFaceColor(i+1,j,k)[0] = FaceColor;   // [-x,-y,-z,x,y,z]
                                }
                                else{
                                    GetElementalFaceColor(i,j,k)[0] = -FaceColor;   // Error color is the -FaceColor
                                    // KRATOS_WARNING("VoxelMesher") << "The given interface is not in the interface of the volume for cell [" << i << "," << j << "," << k << "]" << std::endl;
                                }
                            }
                            previous_center_color = next_center_color;
                        }
                        if(rColors.back() == FaceColor){
                            if(previous_center_color == VolumeColor){
                                GetElementalFaceColor(size_x-1,j,k)[3] = FaceColor;   // [-x,-y,-z,x,y,z]
                            }
                            else{
                                GetElementalFaceColor(size_x-1,j,k)[3] = -FaceColor;   // Error color is the -FaceColor
                                // KRATOS_WARNING << "The given interface is not in the interface of the volume for cell [" << size_x-1 << "," << j << "," << k << "]" << std::endl;
                            }
                        }
                    }
                }
            );
        }

        void CalculateElementalFaceColorsBetweenColors(int FaceColor, int OutsideColor, int VolumeColor){
            if(FaceColor == OutsideColor)
                return; // Nothing to do!

            const std::size_t size_x = mElementCenterCoordinates[0].size();
            const std::size_t size_y = mElementCenterCoordinates[1].size();
            const std::size_t size_z = mElementCenterCoordinates[2].size();

            // First we look in the Z direction
            IndexPartition<std::size_t>(size_x).for_each([&](const std::size_t i) {
                for(std::size_t j = 0 ; j < size_y ; j++){
                    for(std::size_t k = 1 ; k < size_z; k++){
                        auto next_center_color = GetElementalColor(i,j,k);
                        auto current_color = GetElementalColor(i,j,k-1);
                        if((current_color == VolumeColor) && (next_center_color == OutsideColor)){
                            GetElementalFaceColor(i,j,k-1)[5] = FaceColor;   // [-x,-y,-z,x,y,z]
                        }
                        else if((current_color == OutsideColor) && (next_center_color == VolumeColor)){
                            GetElementalFaceColor(i,j,k)[2] = FaceColor;   // [-x,-y,-z,x,y,z]
                        }
                    }
                }
            });

            // Now we look in the Y Direction
            IndexPartition<std::size_t>(size_x).for_each([&](const std::size_t i) {
                for(std::size_t j = 1 ; j < size_y ; j++){
                    for(std::size_t k = 0 ; k < size_z; k++){
                        auto next_center_color = GetElementalColor(i,j,k);
                        auto current_color = GetElementalColor(i,j-1,k);
                        if((current_color == VolumeColor) && (next_center_color == OutsideColor)){
                            GetElementalFaceColor(i,j-1,k)[4] = FaceColor;   // [-x,-y,-z,x,y,z]
                        }
                        else if((current_color == OutsideColor) && (next_center_color == VolumeColor)){
                            GetElementalFaceColor(i,j,k)[1] = FaceColor;   // [-x,-y,-z,x,y,z]
                        }

                    }
                }
            });
            // Now we look in the X Direction
            IndexPartition<std::size_t>(size_x-1).for_each([&](const std::size_t i) {
                for(std::size_t j = 0 ; j < size_y ; j++){
                    for(std::size_t k = 0 ; k < size_z; k++){
                        auto next_center_color = GetElementalColor(i+1,j,k);
                        auto current_color = GetElementalColor(i,j,k);
                        if((current_color == VolumeColor) && (next_center_color == OutsideColor)){
                            GetElementalFaceColor(i,j,k)[3] = FaceColor;   // [-x,-y,-z,x,y,z]
                        }
                        else if((current_color == OutsideColor) && (next_center_color == VolumeColor)){
                            GetElementalFaceColor(i+1,j,k)[0] = FaceColor;   // [-x,-y,-z,x,y,z]
                        }
                    }
                }
            });

        }

        template<typename TPointsContainerType>
        void CalculateOuterMinMaxNodePositions(TPointsContainerType const& Points, array_1d< std::size_t, 3 >& MinNodePosition, array_1d< std::size_t, 3 >& MaxNodePosition) const{
            if(Points.empty())
                return;

            Point min_point;
            Point max_point;

            CalculateMinMaxPoints(Points, min_point, max_point);

            for(std::size_t i = 0; i < 3; i++ ) {
                auto min_position = CalculateNodePosition( min_point[i] - mTolerance, i );
                MinNodePosition[ i ] = (min_position == 0) ? 0 : min_position - 1;
                MaxNodePosition[ i ] = CalculateNodePosition( max_point[i] + mTolerance, i );
            }
        }

        template<typename TPointsContainerType>
        void CalculateMinMaxNodePositions(TPointsContainerType const& Points, array_1d< std::size_t, 3 >& MinNodePosition, array_1d< std::size_t, 3 >& MaxNodePosition){
            if(Points.empty())
                return;

            Point min_point;
            Point max_point;

            CalculateMinMaxPoints(Points, min_point, max_point);

            for(std::size_t i = 0; i < 3; i++ ) {
                MinNodePosition[ i ] = CalculateNodePosition( min_point[i], i );
                MaxNodePosition[ i ] = CalculateNodePosition( max_point[i], i ) + 1;
            }
        }

        template<typename TPointsContainerType>
        void CalculateMinMaxPoints(TPointsContainerType const& Points, Point& MinPoint, Point& MaxPoint) const {
            if(Points.empty())
                return;

            MinPoint = *(Points.begin());
            MaxPoint = *(Points.begin());
            for(auto const& point : Points){
                for(std::size_t i = 0; i<3; i++)
                {
                    MinPoint[i] =  (MinPoint[i] >  point[i] ) ?  point[i] : MinPoint[i];
                    MaxPoint[i] =  (MaxPoint[i] <  point[i] ) ?  point[i] : MaxPoint[i];
                }
            }
        }

        template<typename TPointsContainerType>
        void CalculateMinMaxCenterOfElementPositions(TPointsContainerType const& Points, array_1d< std::size_t, 3 >& MinNodePosition, array_1d< std::size_t, 3 >& MaxNodePosition){
            if(Points.empty())
                return;

            Point min_point = *(Points.begin());
            Point max_point = *(Points.begin());
            for(auto const& point : Points){
                for(std::size_t i = 0; i<3; i++)
                {
                    min_point[i] =  (min_point[i] >  point[i] ) ?  point[i] : min_point[i];
                    max_point[i] =  (max_point[i] <  point[i] ) ?  point[i] : max_point[i];
                }
            }

            for(std::size_t i = 0; i < 3; i++ ) {
                MinNodePosition[ i ] = CalculateCenterOfElementPosition( min_point[i], i );
                MaxNodePosition[ i ] = CalculateCenterOfElementPosition( max_point[i], i ) + 1;
            }
        }

        std::size_t CalculateNodePosition( double Coordinate, int ThisDimension ) const {
            auto const& coordinates = mNodalCoordinates[ThisDimension];
            auto i_min = std::lower_bound(coordinates.begin(), coordinates.end(), Coordinate);
            if(i_min == coordinates.end())
                return coordinates.size() - 1;

            return std::distance(coordinates.begin(), i_min);
        }

        std::size_t CalculateCenterOfElementPosition( double Coordinate, int ThisDimension ) const {
            auto const& coordinates = mElementCenterCoordinates[ThisDimension];
            auto i_min = std::lower_bound(coordinates.begin(), coordinates.end(), Coordinate);
            if(i_min == coordinates.end())
                return coordinates.size() - 1;

            return std::distance(coordinates.begin(), i_min);
         }



        void WriteParaViewVTR(std::string const& Filename) {
            std::ofstream output_file(Filename);

            output_file << "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
            output_file << "<RectilinearGrid WholeExtent=\"0 " << mNodalCoordinates[0].size() - 1  << " 0 " << mNodalCoordinates[1].size() - 1
                        << " 0 " << mNodalCoordinates[2].size() - 1 << "\">" << std::endl;
            output_file << "<Piece Extent=\"0 " << mNodalCoordinates[0].size()  - 1 << " 0 " << mNodalCoordinates[1].size() - 1
                        << " 0 " <<  mNodalCoordinates[2].size() - 1 << "\">" << std::endl;
            output_file << "<Coordinates> " << std::endl;
            output_file << "<DataArray type=\"Float64\" Name=\"coordinates\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
            for (auto i_x : mNodalCoordinates[0]) {
                output_file << i_x << " ";
            }
            output_file << "</DataArray> " << std::endl;
            output_file << "<DataArray type=\"Float64\" Name=\"coordinates\" NumberOfComponents=\"1\" format=\"ascii\">"
                        << std::endl;
            for (auto i_y : mNodalCoordinates[1]) {
                output_file << i_y << " ";
            }
            output_file << "</DataArray> " << std::endl;
            output_file << "<DataArray type=\"Float64\" Name=\"coordinates\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
            for (auto i_z : mNodalCoordinates[2]) {
                output_file << i_z << " ";
            }
            output_file << "</DataArray> " << std::endl;
            output_file << "</Coordinates> " << std::endl;

            output_file << "<PointData Scalars=\"" << "Colors" << "\">" << std::endl;
            for(std::size_t i = 0 ; i < 3 ; i++){
                output_file << "<DataArray type=\"Float64\" Name=\"" << "RayColor" << i << "\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;

                for (auto& color: mNodalRayColors) {
                    output_file << color[i] << " ";
                }
                output_file << std::endl;
                output_file << "</DataArray> " << std::endl;
            }
            output_file << "</PointData> " << std::endl;
            output_file << "<CellData Scalars=\"" << "Color" << "\">" << std::endl;
            for(std::size_t i = 0 ; i < 3 ; i++){
                output_file << "<DataArray type=\"Float64\" Name=\"" << "RayColor" << i << "\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;

                for (auto& color: mElementalRayColors) {
                    output_file << color[i] << " ";
                }
                output_file << std::endl;
                output_file << "</DataArray> " << std::endl;
            }
            for(std::size_t i = 0 ; i < 6 ; i++){
                output_file << "<DataArray type=\"Float64\" Name=\"" << "FaceColor" << i << "\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;

                for (auto& color: mElementalFaceColors) {
                    output_file << color[i] << " ";
                }
                output_file << std::endl;
                output_file << "</DataArray> " << std::endl;
            }
            output_file << "<DataArray type=\"Float64\" Name=\"" << "ElementColor" << "\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;

            for (auto& color: mElementalColors) {
                    output_file << color << " ";
                }
            output_file << std::endl;
            output_file << "</DataArray> " << std::endl;
            output_file << "</CellData> " << std::endl;

            output_file << "</Piece>" << std::endl;
            output_file << "</RectilinearGrid>" << std::endl;
            output_file << "</VTKFile>" << std::endl;
        }

        void PrintAllVoxelColorsToAsciiFile(const std::string& filename) {
            std::ofstream file_out(filename);

            const std::size_t ncols = mElementCenterCoordinates[0].size();
            file_out << "NCOLS "   << ncols <<"\n";
            const std::size_t nrows = mElementCenterCoordinates[1].size();
            file_out << "NROWS "   << nrows <<"\n";
            const std::size_t nlayers = mElementCenterCoordinates[2].size();
            file_out << "NLAYERS " << nlayers <<"\n";
            file_out << "XLLCORNER " << mNodalCoordinates[0][0] <<"\n";
            file_out << "YLLCORNER " << mNodalCoordinates[1][0] <<"\n";
            file_out << "ZLLCORNER " << mNodalCoordinates[2][0] <<"\n";

            file_out << "XCELLSIZE " << (mNodalCoordinates[0].back() - mNodalCoordinates[0][0]) / ncols <<"\n";
            file_out << "YCELLSIZE " << (mNodalCoordinates[1].back() - mNodalCoordinates[1][0]) / nrows <<"\n";
            file_out << "ZCELLSIZE " << (mNodalCoordinates[2].back() - mNodalCoordinates[2][0]) / nlayers <<"\n";
            file_out << "NODATA_VALUE -9999 \n";
            for(std::size_t k=0; k<nlayers; ++k) {
                for(std::size_t j=nrows; j>0; --j) { // ESRI grids consider Y axis downwards (https://en.wikipedia.org/wiki/Esri_grid)
                    for(std::size_t i=0; i<ncols; ++i) {
                        file_out << GetElementalColor(i, j, k) << "  ";
                    }
                    file_out << "\n";
                }
            }
        }
    private:
        array_1d<double, 3>& GetNodalRayColor(std::size_t I, std::size_t J, std::size_t K){
            const std::size_t index = I + J * mNodalCoordinates[0].size() + K * mNodalCoordinates[1].size() * mNodalCoordinates[0].size();
            return mNodalRayColors[index];
        }

        array_1d<double, 3>& GetElementalRayColor(std::size_t I, std::size_t J, std::size_t K){
            const std::size_t index = I + J * (mNodalCoordinates[0].size() - 1) + K * (mNodalCoordinates[1].size() - 1) * (mNodalCoordinates[0].size() - 1);
            return mElementalRayColors[index];
        }

    };
  }
}  // namespace Kratos.
