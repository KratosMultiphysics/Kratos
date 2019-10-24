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
//

#if !defined(KRATOS_INTERNALS_CARTESIAN_MESH_COLORS_H_INCLUDED )
#define  KRATOS_INTERNALS_CARTESIAN_MESH_COLORS_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <algorithm>


// External includes


// Project includes
#include "processes/process.h"
#include "includes/model_part.h"
#include "processes/internals/cartesian_ray.h"



namespace Kratos
{
  namespace Internals
  {
    class CartesianMeshColors{
      double mTolerance;
      array_1d<std::vector<double>, 3> mNodalCoordinates;
      array_1d<std::vector<double>, 3> mElementCenterCoordinates;
      std::vector<array_1d<double,3>> mNodalRayColors;
      std::vector<array_1d<double,3>> mElementalRayColors;
      std::vector<double> mNodalColors;
      std::vector<double> mElementalColors;
      DenseMatrix<Internals::CartesianRay<Element::GeometryType>> mXYRays;
      DenseMatrix<Internals::CartesianRay<Element::GeometryType>> mXZRays;
      DenseMatrix<Internals::CartesianRay<Element::GeometryType>> mYZRays;
     public:
      CartesianMeshColors(): mTolerance(1e-12){}

      std::vector<double> const& GetNodalCoordinates(int Index) const {return mNodalCoordinates[Index];}

      void SetCoordinates(array_1d<std::vector<double>, 3>& TheNodalCoordinates){

        KRATOS_ERROR_IF((TheNodalCoordinates[0].size() < 2) || (TheNodalCoordinates[1].size() < 2) || (TheNodalCoordinates[2].size() < 2)) 
            << "The coordinates should have at least two entries defining the bounding box." << std::endl;

        mNodalCoordinates = TheNodalCoordinates;

        for(int i = 0 ; i < 3 ; i++){
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

        SetAllColors(0.00);
     }

    void SetAllColors(double TheColor){

        const int number_of_nodes = static_cast<int>(mNodalColors.size());
        for(int i = 0 ; i < number_of_nodes ; i++ ){
            mNodalRayColors[i] = ScalarVector(3,TheColor);
            mNodalColors[i] = TheColor;
        }

        const int number_of_elements = static_cast<int>(mElementalColors.size());
        for(int i = 0 ; i < number_of_elements ; i++ ){
            mElementalRayColors[i] = ScalarVector(3,TheColor);
            mElementalColors[i] = TheColor;
        }
    }

    double& GetNodalColor(std::size_t I, std::size_t J, std::size_t K){
        const std::size_t index = I + J * mNodalCoordinates[0].size() + K * mNodalCoordinates[1].size() * mNodalCoordinates[0].size();
        return mNodalColors[index];
    }

    double& GetElementalColor(std::size_t I, std::size_t J, std::size_t K){
        const std::size_t index = I + J * mElementCenterCoordinates[0].size() + K * mElementCenterCoordinates[1].size() * mElementCenterCoordinates[0].size();
        return mElementalColors[index];
    }

    Point GetPoint(std::size_t I, std::size_t J, std::size_t K){
        return Point(mNodalCoordinates[0][I], mNodalCoordinates[1][J], mNodalCoordinates[2][K]);
    }

    Point GetCenterOfElement(std::size_t I, std::size_t J, std::size_t K){
        return Point(mElementCenterCoordinates[0][I], mElementCenterCoordinates[1][J], mElementCenterCoordinates[2][K]);
    }
      
    CartesianRay<Element::GeometryType>& GetXYRay(std::size_t I, std::size_t J){
        return mXYRays(I,J);
    }

    void InitializeRays(array_1d< std::size_t, 3 > const& MinRayPosition, array_1d< std::size_t, 3 > const& MaxRayPosition, std::string EntititesToColor){

        if(EntititesToColor == "center_of_elements"){
            for(std::size_t i = MinRayPosition[0] ; i < MaxRayPosition[0] ; i++){
                for(std::size_t j = MinRayPosition[1] ; j < MaxRayPosition[1] ; j++){
                mXYRays(i,j) = Internals::CartesianRay<Element::GeometryType>(2, GetCenterOfElement(i,j,0), GetCenterOfElement(i,j,mNodalCoordinates[2].size() - 2));
                }
            }
            for(std::size_t i = MinRayPosition[0] ; i < MaxRayPosition[0] ; i++){
                for(std::size_t k = MinRayPosition[2] ; k < MaxRayPosition[2] ; k++){
                mXZRays(i,k) = Internals::CartesianRay<Element::GeometryType>(1, GetCenterOfElement(i,0, k), GetCenterOfElement(i,mNodalCoordinates[1].size() - 2,k));
                }
            }
            for(std::size_t j = MinRayPosition[1] ; j < MaxRayPosition[1] ; j++){
                for(std::size_t k = MinRayPosition[2] ; k < MaxRayPosition[2] ; k++){
                mYZRays(j,k) = Internals::CartesianRay<Element::GeometryType>(0, GetCenterOfElement(0,j,k), GetCenterOfElement(mNodalCoordinates[0].size() - 2,j,k));
                }
            }
        }
        else if(EntititesToColor == "nodes"){
            for(std::size_t i = MinRayPosition[0] ; i < MaxRayPosition[0] ; i++){
                for(std::size_t j = MinRayPosition[1] ; j < MaxRayPosition[1] ; j++){
                mXYRays(i,j) = Internals::CartesianRay<Element::GeometryType>(2, GetPoint(i,j,0), GetPoint(i,j,mNodalCoordinates[2].size() - 1));
                }
            }
            for(std::size_t i = MinRayPosition[0] ; i < MaxRayPosition[0] ; i++){
                for(std::size_t k = MinRayPosition[2] ; k < MaxRayPosition[2] ; k++){
                mXZRays(i,k) = Internals::CartesianRay<Element::GeometryType>(1, GetPoint(i,0, k), GetPoint(i,mNodalCoordinates[1].size() - 1,k));
                }
            }
            for(std::size_t j = MinRayPosition[1] ; j < MaxRayPosition[1] ; j++){
                for(std::size_t k = MinRayPosition[2] ; k < MaxRayPosition[2] ; k++){
                mYZRays(j,k) = Internals::CartesianRay<Element::GeometryType>(0, GetPoint(0,j,k), GetPoint(mNodalCoordinates[0].size() - 1,j,k));
                }
            }
        }
        else{
            KRATOS_ERROR << "Undefined entities type: \"" << EntititesToColor << "\" The possible options are:  \"ndoe\" and  \"center_of_elements\"" << std::endl;
        }
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
                
                for(std::size_t i = min_position[0] ; i < max_position[0] ; i++){
                    for(std::size_t j = min_position[1] ; j < max_position[1] ; j++){
                        mXYRays(i,j).AddIntersection(rGeometry, mTolerance);
                    }
                }
                for(std::size_t i = min_position[0] ; i < max_position[0] ; i++){
                    for(std::size_t k = min_position[2] ; k < max_position[2] ; k++){
                        mXZRays(i,k).AddIntersection(rGeometry, mTolerance);
                    }
                }
                for(std::size_t j = min_position[1] ; j < max_position[1] ; j++){
                    for(std::size_t k = min_position[2] ; k < max_position[2] ; k++){
                        mYZRays(j,k).AddIntersection(rGeometry, mTolerance);
                    }
                }
        
        }

        void CalculateNodalRayColors(array_1d< std::size_t, 3 > const& MinRayPosition, array_1d< std::size_t, 3 > const& MaxRayPosition, int InsideColor, int OutsideColor){
            std::vector<double> colors;
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
            for(std::size_t j = MinRayPosition[1] ; j < MaxRayPosition[1] ; j++){
                for(std::size_t k = MinRayPosition[2] ; k < MaxRayPosition[2] ; k++){
                    mYZRays(j,k).CollapseIntersectionPoints(mTolerance);
                    mYZRays(j,k).CalculateColor(mNodalCoordinates[0], InsideColor, OutsideColor, colors, mTolerance);
                    for(std::size_t i = MinRayPosition[0] ; i < MaxRayPosition[0] ; i++){
                        if(colors[i] == InsideColor){
                            GetNodalRayColor(i,j,k)[2] = colors[i];
                            // KRATOS_WATCH(GetColor(i,j,k)[2]);
                        }
                    }
                }
            }

            for(std::size_t k = MinRayPosition[2] ; k < MaxRayPosition[2] ; k++){
                for(std::size_t j = MinRayPosition[1] ; j < MaxRayPosition[1] ; j++){
                    for(std::size_t i = MinRayPosition[0] ; i < MaxRayPosition[0] ; i++){
                        auto& ray_colors = GetNodalRayColor(i,j,k);
                        std::size_t n_inside = 0;
                        std::size_t n_outside = 0;
                        for(int dim = 0 ; dim < 3 ; dim++){
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

            for(std::size_t i = MinRayPosition[0] ; i < MaxRayPosition[0] ; i++){
                for(std::size_t j = MinRayPosition[1] ; j < MaxRayPosition[1] ; j++){
                    auto& ray = mXYRays(i,j);
                    ray.CollapseIntersectionPoints(mTolerance);
                    ray.CalculateColor(z_coordinates, InsideColor, OutsideColor, colors, mTolerance);
                    for(std::size_t k = MinRayPosition[2] ; k < MaxRayPosition[2] ; k++){
                        if(colors[k] == InsideColor){
                            GetElementalRayColor(i,j,k)[0] = colors[k];
                            // KRATOS_WATCH(color);
                        }
                    }
                }
            }
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
            for(std::size_t j = MinRayPosition[1] ; j < MaxRayPosition[1] ; j++){
                for(std::size_t k = MinRayPosition[2] ; k < MaxRayPosition[2] ; k++){
                    mYZRays(j,k).CollapseIntersectionPoints(mTolerance);
                    mYZRays(j,k).CalculateColor(x_coordinates, InsideColor, OutsideColor, colors, mTolerance);
                    for(std::size_t i = MinRayPosition[0] ; i < MaxRayPosition[0] ; i++){
                        if(colors[i] == InsideColor){
                            GetElementalRayColor(i,j,k)[2] = colors[i];
                            // KRATOS_WATCH(GetColor(i,j,k)[2]);
                        }
                    }
                }
            }
            for(std::size_t k = MinRayPosition[2] ; k < MaxRayPosition[2] ; k++){
                for(std::size_t j = MinRayPosition[1] ; j < MaxRayPosition[1] ; j++){
                    for(std::size_t i = MinRayPosition[0] ; i < MaxRayPosition[0] ; i++){
                        auto& ray_colors = GetElementalRayColor(i,j,k);
                        std::size_t n_inside = 0;
                        std::size_t n_outside = 0;
                        for(int dim = 0 ; dim < 3 ; dim++){
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

        template<typename TPointsContainerType>
        void CalculateMinMaxNodePositions(TPointsContainerType const& Points, array_1d< std::size_t, 3 >& MinNodePosition, array_1d< std::size_t, 3 >& MaxNodePosition){
            if(Points.empty())
                return;

            Point min_point;
            Point max_point;
            max_point = *(Points.begin());
            min_point = *(Points.begin());
            for(auto const& point : Points){
                for(std::size_t i = 0; i<3; i++)
                {
                    min_point[i] =  (min_point[i] >  point[i] ) ?  point[i] : min_point[i];
                    max_point[i] =  (max_point[i] <  point[i] ) ?  point[i] : max_point[i];
                }
            }
                
            for(int i = 0; i < 3; i++ ) {
                MinNodePosition[ i ] = CalculateNodePosition( min_point[i], i );
                MaxNodePosition[ i ] = CalculateNodePosition( max_point[i], i ) + 1;
            }
        }

        template<typename TPointsContainerType>
        void CalculateMinMaxCenterOfElementPositions(TPointsContainerType const& Points, array_1d< std::size_t, 3 >& MinNodePosition, array_1d< std::size_t, 3 >& MaxNodePosition){
            if(Points.empty())
                return;

            Point min_point;
            Point max_point;
            max_point = *(Points.begin());
            min_point = *(Points.begin());
            for(auto const& point : Points){
                for(std::size_t i = 0; i<3; i++)
                {
                    min_point[i] =  (min_point[i] >  point[i] ) ?  point[i] : min_point[i];
                    max_point[i] =  (max_point[i] <  point[i] ) ?  point[i] : max_point[i];
                }
            }
                
            for(int i = 0; i < 3; i++ ) {
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
            for(int i = 0 ; i < 3 ; i++){
                output_file << "<DataArray type=\"Float64\" Name=\"" << "RayColor" << i << "\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
                
                for (auto& color: mNodalRayColors) {
                    output_file << color[i] << std::endl;
                }
                    
                output_file << "</DataArray> " << std::endl;
            }
            output_file << "</PointData> " << std::endl;
            output_file << "<CellData Scalars=\"" << "Color" << "\">" << std::endl;
            for(int i = 0 ; i < 3 ; i++){
                output_file << "<DataArray type=\"Float64\" Name=\"" << "RayColor" << i << "\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
                
                for (auto& color: mElementalRayColors) {
                    output_file << color[i] << std::endl;
                }
                    
                output_file << "</DataArray> " << std::endl;
            }
            output_file << "<DataArray type=\"Float64\" Name=\"" << "RayColor" << "\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
                
            for (auto& color: mElementalColors) {
                output_file << color << std::endl;
            }
                    
            output_file << "</DataArray> " << std::endl;
            output_file << "</CellData> " << std::endl;

            output_file << "</Piece>" << std::endl;
            output_file << "</RectilinearGrid>" << std::endl;
            output_file << "</VTKFile>" << std::endl;
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

#endif // KRATOS_INTERNALS_CARTESIAN_MESH_COLORS_H_INCLUDED  defined
