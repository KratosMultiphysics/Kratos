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



namespace Kratos
{
  namespace Internals
  {
    template<typename TGeometryType> 
    class CartesianRay{
      int mDirection;
      Point mPoint1;
      Point mPoint2;
      std::vector<std::pair<double, const TGeometryType*>> mIntersections;
    public:
      CartesianRay(): mDirection(0), mPoint1(), mPoint2() {}

      CartesianRay(int Direction, Point const& Point1, Point const& Point2): mDirection(Direction), mPoint1(Point1), mPoint2(Point2) {}

      CartesianRay(CartesianRay const& Other): mDirection(Other.mDirection), mPoint1(Other.mPoint1), mPoint2(Other.mPoint2), mIntersections(Other.mIntersections) {}

      CartesianRay& operator=(CartesianRay const& Other){
          mDirection=Other.mDirection;
          mPoint1 = Other.mPoint1;
          mPoint2 = Other.mPoint2;
          mIntersections = Other.mIntersections;

          return *this;
      }

      void AddIntersection(TGeometryType const& rGeometry, double Tolerance){
        // Call the line - triangle intersection util
        array_1d<double,3> intersection_point = ZeroVector(3);
        const int is_intersected = IntersectionUtilities::ComputeTriangleLineIntersection(
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
        if (mIntersections.empty()) {
            return;
        }

        std::size_t new_size = 0;
        // Sort
        std::sort(mIntersections.begin(), mIntersections.end());
        // Unique
        auto i_begin = mIntersections.begin();
        auto i_intersection = mIntersections.begin();
        while (++i_begin != mIntersections.end()) {
            // considering the very near points as the same points
            if (std::abs(i_begin->first - i_intersection->first) > Tolerance) { // if the hit points are far enough they are not the same                
                *(++i_intersection) = *i_begin;
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

                if(CheckPassingThroughByExtraRays(i_patch_begin, i_patch_end, Tolerance,  10.00*Tolerance)) {
                    *(++i_intersection) = *i_begin;
                    new_size++;
                }
            }
            if(i_begin == mIntersections.end())
                break;
        }
        mIntersections.resize(new_size);
      }

      void CalculateColor(std::vector<double> const& Coordinates, int InsideColor, int OutsideColor, std::vector<double>& ResultingColors){
        bool is_inside=false;

        if(ResultingColors.size() != Coordinates.size())
        ResultingColors.resize(Coordinates.size());
    
        std::size_t current_index=0;
        const std::size_t size = Coordinates.size();

        for(auto& i_intersection : mIntersections){
            while(current_index < size){
                if(Coordinates[current_index] < i_intersection.first){
                    ResultingColors[current_index++] = (is_inside) ? InsideColor : OutsideColor;
                }
                else{
                    current_index++;
                    break;
                }
            }
            is_inside = (is_inside) ? false:true;
        }
        // now coloring the points after the last intersection as outside
        while(current_index < size){
            ResultingColors[current_index++] = OutsideColor;
        }
      }

      std::vector<std::pair<double, const TGeometryType*>> const& GetIntersections() const {return mIntersections;}

      Point& GetPoint1(){return mPoint1;}

      Point& GetPoint2(){return mPoint2;}

    private:

      bool CheckPassingThroughByExtraRays(typename std::vector<std::pair<double, const TGeometryType*>>::iterator Begin, typename std::vector<std::pair<double, const TGeometryType*>>::iterator End, double Tolerance, double Delta){
        std::array<double, 8> delta_u{Delta, Delta, 0.00, -Delta, -Delta, -Delta, 0.00, Delta};
        std::array<double, 8> delta_v{0.00, Delta, Delta, Delta, 0.00, -Delta, -Delta, -Delta};

        std::size_t no_hit_cases = 0;

        for(int i_ray = 0 ; i_ray < 8 ; i_ray++){
            CartesianRay extra_ray(mDirection, mPoint1, mPoint2);
            if(mDirection == 0) {
                extra_ray.mPoint1[1] += delta_u[i_ray];
                extra_ray.mPoint1[2] += delta_v[i_ray];
                extra_ray.mPoint2[1] += delta_u[i_ray];
                extra_ray.mPoint2[2] += delta_v[i_ray];
            }
            else if(mDirection == 1) {
                extra_ray.mPoint1[0] += delta_u[i_ray];
                extra_ray.mPoint1[2] += delta_v[i_ray];
                extra_ray.mPoint2[0] += delta_u[i_ray];
                extra_ray.mPoint2[2] += delta_v[i_ray];
            }
            else if(mDirection == 2) {
                extra_ray.mPoint1[0] += delta_u[i_ray];
                extra_ray.mPoint1[1] += delta_v[i_ray];
                extra_ray.mPoint2[0] += delta_u[i_ray];
                extra_ray.mPoint2[1] += delta_v[i_ray];
            }
            for(auto i_intersection = Begin ; i_intersection != End; i_intersection++){
                extra_ray.AddIntersection(*(i_intersection->second), Tolerance);
            }
            if(extra_ray.mIntersections.size() == 0){
                no_hit_cases++;
            }
            if(no_hit_cases > 1) // I admit only 1 ray passing through a gap
                return false;
        }
        return true;
      }
    };

    class CartesianMeshColors{
      double mTolerance;
      array_1d<std::vector<double>, 3> mCoordinates;
      std::vector<array_1d<double,3>> mColors;
      DenseMatrix<Internals::CartesianRay<Element::GeometryType>> mXYRays;
      DenseMatrix<Internals::CartesianRay<Element::GeometryType>> mXZRays;
      DenseMatrix<Internals::CartesianRay<Element::GeometryType>> mYZRays;
     public:
      CartesianMeshColors(): mTolerance(1e-6){}

      std::vector<double> const& GetCoordinates(int Index) const {return mCoordinates[Index];}

      void SetCoordinates(array_1d<std::vector<double>, 3>& TheCoordinates){
        mCoordinates = TheCoordinates;

        mXYRays.resize(mCoordinates[0].size(),mCoordinates[1].size(),false);
        mXZRays.resize(mCoordinates[0].size(),mCoordinates[2].size(),false);
        mYZRays.resize(mCoordinates[1].size(),mCoordinates[2].size(),false); 

        mColors.resize(mCoordinates[0].size()*mCoordinates[1].size()*mCoordinates[2].size());

        for(auto& i_color : mColors){
            i_color = ScalarVector(3,0.00);
        }
      }

      array_1d<double, 3>& GetColor(std::size_t I, std::size_t J, std::size_t K){
        const std::size_t index = I + J * mCoordinates[0].size() + K * mCoordinates[1].size() * mCoordinates[0].size();
        return mColors[index];
      }

      Point GetPoint(std::size_t I, std::size_t J, std::size_t K){
        return Point(mCoordinates[0][I], mCoordinates[1][J], mCoordinates[2][K]);
      }
      
      CartesianRay<Element::GeometryType>& GetXYRay(std::size_t I, std::size_t J){
          return mXYRays(I,J);
      }

        void InitializeRays(array_1d< std::size_t, 3 > const& MinRayPosition, array_1d< std::size_t, 3 > const& MaxRayPosition){

        array_1d<std::size_t, 3> number_of_rays = MaxRayPosition - MinRayPosition;
        number_of_rays += scalar_vector<std::size_t>(3,1);

        for(std::size_t i = MinRayPosition[0] ; i < MaxRayPosition[0] ; i++){
            for(std::size_t j = MinRayPosition[1] ; j < MaxRayPosition[1] ; j++){
            mXYRays(i,j) = Internals::CartesianRay<Element::GeometryType>(2, GetPoint(i,j,MinRayPosition[2]), GetPoint(i,j,MaxRayPosition[2] - 1));
            }
        }
        for(std::size_t i = MinRayPosition[0] ; i < MaxRayPosition[0] ; i++){
            for(std::size_t k = MinRayPosition[2] ; k < MaxRayPosition[2] ; k++){
            mXZRays(i,k) = Internals::CartesianRay<Element::GeometryType>(1, GetPoint(i,MinRayPosition[1], k), GetPoint(i,MaxRayPosition[1] - 1,k));
            }
        }
        for(std::size_t j = MinRayPosition[1] ; j < MaxRayPosition[1] ; j++){
            for(std::size_t k = MinRayPosition[2] ; k < MaxRayPosition[2] ; k++){
            mYZRays(j,k) = Internals::CartesianRay<Element::GeometryType>(0, GetPoint(MinRayPosition[0],j,k), GetPoint(MaxRayPosition[0] - 1,j,k));
            }
        }
        }

        void AddGeometry(Element::GeometryType const& rGeometry){
                array_1d< std::size_t, 3 > min_position;
                array_1d< std::size_t, 3 > max_position;
                CalculateMinMaxCellsPositions(rGeometry, min_position, max_position);
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

        void CalculateColors(array_1d< std::size_t, 3 > const& MinRayPosition, array_1d< std::size_t, 3 > const& MaxRayPosition, int InsideColor, int OutsideColor){
            std::vector<double> colors;
            for(std::size_t i = MinRayPosition[0] ; i < MaxRayPosition[0] ; i++){
                for(std::size_t j = MinRayPosition[1] ; j < MaxRayPosition[1] ; j++){
                    mXYRays(i,j).CollapseIntersectionPoints(mTolerance);
                    mXYRays(i,j).CalculateColor(mCoordinates[2], InsideColor, OutsideColor, colors);
                    for(std::size_t k = MinRayPosition[2] ; k < MaxRayPosition[2] ; k++){
                        if(colors[k] == InsideColor){
                            auto& color=GetColor(i,j,k)[0];
                            GetColor(i,j,k)[0] = colors[k];
                            KRATOS_WATCH(color);
                        }
                    }
                }
            }
            for(std::size_t i = MinRayPosition[0] ; i < MaxRayPosition[0] ; i++){
                for(std::size_t k = MinRayPosition[2] ; k < MaxRayPosition[2] ; k++){
                    mXZRays(i,k).CollapseIntersectionPoints(mTolerance);
                    mXZRays(i,k).CalculateColor(mCoordinates[1], InsideColor, OutsideColor, colors);
                    for(std::size_t j = MinRayPosition[1] ; j < MaxRayPosition[1] ; j++){
                        if(colors[k] == InsideColor){
                            GetColor(i,j,k)[1] = colors[j];
                        }
                    }
                }
            }
            for(std::size_t j = MinRayPosition[1] ; j < MaxRayPosition[1] ; j++){
                for(std::size_t k = MinRayPosition[2] ; k < MaxRayPosition[2] ; k++){
                    mYZRays(j,k).CollapseIntersectionPoints(mTolerance);
                    mYZRays(j,k).CalculateColor(mCoordinates[0], InsideColor, OutsideColor, colors);
                    for(std::size_t i = MinRayPosition[0] ; i < MaxRayPosition[0] ; i++){
                        if(colors[k] == InsideColor){
                            GetColor(i,j,k)[2] = colors[i];
                            KRATOS_WATCH(GetColor(i,j,k)[2]);
                        }
                    }
                }
            }
        }

        template<typename TPointsContainerType>
        void CalculateMinMaxCellsPositions(TPointsContainerType const& Points, array_1d< std::size_t, 3 >& MinCellPosition, array_1d< std::size_t, 3 >& MaxCellPosition){
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
                
            for ( int i = 0; i < 3; i++ ) {
                MinCellPosition[ i ] = CalculatePosition( min_point[i], i );
                MaxCellPosition[ i ] = CalculatePosition( max_point[i], i ) + 1;
            }
        }

        std::size_t CalculatePosition( double Coordinate, int ThisDimension ) const {
            auto const& coordinates = mCoordinates[ThisDimension];
            auto i_min = std::lower_bound(coordinates.begin(), coordinates.end(), Coordinate);
            if(i_min == coordinates.end())
            return coordinates.size() - 1;

            return std::distance(coordinates.begin(), i_min);
        }

        

        void WriteParaViewVTS(std::string const& Filename) {
            std::ofstream output_file(Filename + ".vtr");

            output_file << "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
            output_file << "<RectilinearGrid WholeExtent=\"0 " << mCoordinates[0].size() - 1  << " 0 " << mCoordinates[1].size() - 1 
                        << " 0 " << mCoordinates[2].size() - 1 << "\">" << std::endl;
            output_file << "<Piece Extent=\"0 " << mCoordinates[0].size()  - 1 << " 0 " << mCoordinates[1].size() - 1 
                        << " 0 " <<  mCoordinates[2].size() - 1 << "\">" << std::endl;
            output_file << "<Coordinates> " << std::endl;
            output_file << "<DataArray type=\"Float64\" Name=\"coordinates\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
            for (auto i_x : mCoordinates[0]) {
                output_file << i_x << " ";
            }
            output_file << "</DataArray> " << std::endl;
            output_file << "<DataArray type=\"Float64\" Name=\"coordinates\" NumberOfComponents=\"1\" format=\"ascii\">"
                        << std::endl;
            for (auto i_y : mCoordinates[1]) {
                output_file << i_y << " ";
            }
            output_file << "</DataArray> " << std::endl;
            output_file << "<DataArray type=\"Float64\" Name=\"coordinates\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
            for (auto i_z : mCoordinates[2]) {
                output_file << i_z << " ";
            }
            output_file << "</DataArray> " << std::endl;
            output_file << "</Coordinates> " << std::endl;

            for(int i = 0 ; i < 3 ; i++){
                output_file << "<PointData Scalars=\"" << "Color" << i << "\">" << std::endl;
                output_file << "<DataArray type=\"Float64\" Name=\"" << "Color" << i << "\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
                
                for (auto& color: mColors) {
                    output_file << color[i] << std::endl;
                }
                    
                output_file << "</DataArray> " << std::endl;
                output_file << "</PointData> " << std::endl;
            }

            output_file << "</Piece>" << std::endl;
            output_file << "</RectilinearGrid>" << std::endl;
            output_file << "</VTKFile>" << std::endl;
        }

    };
  }
}  // namespace Kratos.

#endif // KRATOS_INTERNALS_CARTESIAN_MESH_COLORS_H_INCLUDED  defined
