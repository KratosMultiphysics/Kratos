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

// System includes
#include <iostream>
#include <fstream>
#include <filesystem>
#include <ios>
#include <vector>


// External includes

// Project includes

namespace Kratos
{
        /** \brief VtkOutput
         * A simple class that has functionality to write vtk outputs for container variable data - extention of the kratos core vtk output
         * @see : https://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
         */

        /**
         * @brief Print the given Data as a VTK file
         */
        class GenericVtkOutput
        {
        public:
                explicit GenericVtkOutput();

                /// Destructor.
                virtual ~GenericVtkOutput() = default;

                void outputStructuredPoints(const std::string &path,
                int numPointsX, int numPointsY, int numPointsZ,
                float spaceBetweenPointsX, float spaceBetweenPointsY, float spaceBetweenPointsZ,
                const std::string scalarDataName,
                const std::string vectorDataName,
                const std::string tensorDataName,
                const std::vector<float> &scalarData,
                const std::vector<std::vector<float>> &vectorData,
                const std::vector<std::vector<std::vector<float>>> &tensorData);
        private:
                std::ofstream createFileStreamFromPath(const std::string &path);
        };
}