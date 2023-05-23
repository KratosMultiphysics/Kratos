//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Fabian Meister
//

// Include base h
#include "generic_vtk_output.h"


namespace Kratos
{

        GenericVtkOutput::GenericVtkOutput(){};

        void GenericVtkOutput::outputStructuredPoints(const std::string &path,
                                                      int numPointsX, int numPointsY, int numPointsZ,
                                                      float spaceBetweenPointsX = 1., float spaceBetweenPointsY = 1., float spaceBetweenPointsZ = 1.,
                                                      const std::string scalarDataName = "scalars",
                                                      const std::string vectorDataName = "vectors",
                                                      const std::string tensorDataName = "tensors",
                                                      const std::vector<float> &scalarData = {},
                                                      const std::vector<std::vector<float>> &vectorData = {{}},
                                                      const std::vector<std::vector<std::vector<float>>> &tensorData = {{{}}})
        {
                auto output_file = createFileStreamFromPath(path);

                output_file << "# vtk DataFile Version 4.0" << std::endl;
                output_file << "vtk output" << std::endl;
                output_file << "ASCII" << std::endl;
                output_file << "DATASET STRUCTURED_POINTS" << std::endl;
                output_file << "DIMENSIONS " << numPointsX << " " << numPointsY << " " << numPointsZ << std::endl;
                output_file << "ORIGIN 0.0 0.0 0.0" << std::endl;
                output_file << "SPACING " << spaceBetweenPointsX << " " << spaceBetweenPointsY << " " << spaceBetweenPointsZ << std::endl;

                output_file << std::endl;

                int numOfPoints = numPointsX*numPointsY*numPointsZ;
                output_file << "POINT_DATA " << numOfPoints << std::endl;


                if(scalarData.size() > 0){
                        if(scalarData.size() != numOfPoints){
                                std::printf("The number of scalar data values (%d) does not match the number of points (%d)",scalarData.size(),numOfPoints);
                        }
                        output_file << "SCALARS "<< scalarDataName <<" float 1" << std::endl;
                        output_file << "LOOKUP_TABLE default" << std::endl;
                        for (auto value:scalarData ){
                                output_file << value << " ";
                        }
                        output_file << std::endl;
                        output_file << std::endl;
                }

                if(vectorData.size() > 0){
                        if(vectorData.size() != numOfPoints){
                                std::printf("The number of vector data values (%d) does not match the number of points (%d)",vectorData.size(),numOfPoints);
                        }
                        output_file << "VECTORS "<< vectorDataName <<" float" << std::endl;
                        int counter = 0;
                        for (auto valuesPerPoint:vectorData ){
                                counter = 0;
                                for(auto value:valuesPerPoint){
                                        output_file << value << " ";
                                        counter ++;
                                }
                                for (counter; counter<3; counter++){
                                        output_file << "0.0 ";
                                }
                                output_file << std::endl;
                        }
                        output_file << std::endl;
                }

                if(tensorData.size() > 0){
                        if(tensorData.size() != numOfPoints){
                                std::printf("The number of tensor data values (%d) does not match the number of points (%d)",tensorData.size(),numOfPoints);
                        }
                        output_file << "TENSORS "<< tensorDataName <<" float" << std::endl;
                        int counter_first_dim = 0;
                        int counter_sec_dim = 0;
                        for (auto valuesPerPoint:tensorData ){
                                counter_first_dim = 0;
                                for(auto values_first_dim:valuesPerPoint){
                                        counter_sec_dim = 0;
                                        counter_first_dim++;
                                        for(auto value:values_first_dim){
                                                output_file << value << " ";
                                                counter_sec_dim ++;
                                        }
                                        for (counter_sec_dim; counter_sec_dim<3; counter_sec_dim++){
                                                output_file << "0.0 ";
                                        }
                                        // output_file << std::endl;
                                }
                                for (counter_first_dim; counter_first_dim<3; counter_first_dim++){
                                        output_file << "0.0 0.0 0.0";
                                }

                                output_file << std::endl;
                        }
                        output_file << std::endl;
                }


                output_file.close();
        };

        std::ofstream GenericVtkOutput::createFileStreamFromPath(const std::string &path){
                auto posOfFileNameStart = path.find_last_of("/");
                if(posOfFileNameStart < path.length()){
                        std::string directories = path.substr( 0,posOfFileNameStart);
                        std::filesystem::create_directories( directories);
                }

                std::remove(path.c_str());
                std::ofstream output_file;
                output_file.open (path, std::fstream::app);
                return output_file;
        }
}