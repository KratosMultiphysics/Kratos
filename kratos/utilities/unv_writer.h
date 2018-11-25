#ifndef KRATOSMULTIPHYSICS_UNV_WRITER_H_INCLUDED
#define KRATOSMULTIPHYSICS_UNV_WRITER_H_INCLUDED

#include <iostream>
#include <fstream>
#include <includes/model_part.h>
#include "includes/exception.h"

namespace Kratos {


class UnvWriter {
public:

    UnvWriter(Kratos::ModelPart &modelPart, const std::string &outFileWithoutExtension)
            : mrOutputModelPart(modelPart),
              mOutputFileName(outFileWithoutExtension + ".unv") {
    }


    void writeMesh() {
        initializeOutputFile();
        writeNodes();
        writeElements();
    }

    void initializeOutputFile() {
        std::ofstream outputFile;
        outputFile.open(mOutputFileName, std::ios::out | std::ios::trunc);
        outputFile.close();
    }

    void writeNodes() {
        std::ofstream outputFile;
        outputFile.open(mOutputFileName, std::ios::out | std::ios::app);

        outputFile << std::scientific;
        outputFile << std::setprecision(15);

        const int dataSetNumberForNodes = 2411;
        const int exportCoordinateSystemNumber = 0;
        const int displacementCoordinateSystemNumber = 0;
        const int color = 0;


        outputFile << std::setw(6) << "-1" << "\n";
        outputFile << std::setw(6) << dataSetNumberForNodes << "\n";


        for (auto &node_i : mrOutputModelPart.Nodes()) {
            int node_label = node_i.Id();
            double x_coordinate = node_i.X();
            double y_coordinate = node_i.Y();
            double z_coordinate = node_i.Z();
            outputFile << std::setw(10) << node_label << std::setw(10) << exportCoordinateSystemNumber
                       << std::setw(10)
                       << displacementCoordinateSystemNumber << std::setw(10) << color << "\n";
            outputFile << std::setw(25) << x_coordinate << std::setw(25) << y_coordinate << std::setw(25)
                       << z_coordinate << "\n";
        }
        outputFile << std::setw(6) << "-1" << "\n";

        outputFile.close();
    }


    void writeElements() {
        std::ofstream outputFile;
        outputFile.open(mOutputFileName, std::ios::out | std::ios::app);

        const int dataSetNumberForElements = 2412;
        const int physicalPropertyTableNumber = 1;
        const int materialPropertyTableNumber = 1;
        const int color = 0;

        outputFile << std::setw(6) << "-1" << "\n";
        outputFile << std::setw(6) << dataSetNumberForElements << "\n";

        for (auto &element : mrOutputModelPart.Elements()) {
            const int elementLabel = element.Id();
            Kratos::ModelPart::ConditionType::GeometryType elementGeometry = element.GetGeometry();
            // Write triangles
            if (elementGeometry.size() == 3 && elementGeometry.Dimension() == 2) {
                const int feDescriptorId = 41; // Plane Stress Linear Triangle
                const int numberOfNodes = 3;
                outputFile << std::setw(10) << elementLabel;
                outputFile << std::setw(10) << feDescriptorId;
                outputFile << std::setw(10) << physicalPropertyTableNumber;
                outputFile << std::setw(10) << materialPropertyTableNumber;
                outputFile << std::setw(10) << color;
                outputFile << std::setw(10) << numberOfNodes << "\n";
                outputFile << std::setw(10) << elementGeometry[0].Id();
                outputFile << std::setw(10) << elementGeometry[1].Id();
                outputFile << std::setw(10) << elementGeometry[2].Id() << "\n";
            }
                // Write tetrahedras
            else if (elementGeometry.size() == 4 && elementGeometry.Dimension() == 3) {
                const int feDescriptorId = 111; // Solid linear tetrahedron
                const int numberOfNodes = 4;
                outputFile << std::setw(10) << elementLabel;
                outputFile << std::setw(10) << feDescriptorId;
                outputFile << std::setw(10) << physicalPropertyTableNumber;
                outputFile << std::setw(10) << materialPropertyTableNumber;
                outputFile << std::setw(10) << color;
                outputFile << std::setw(10) << numberOfNodes << "\n";
                outputFile << std::setw(10) << elementGeometry[0].Id();
                outputFile << std::setw(10) << elementGeometry[1].Id();
                outputFile << std::setw(10) << elementGeometry[2].Id();
                outputFile << std::setw(10) << elementGeometry[3].Id() << "\n";
            }
        }
        outputFile << std::setw(6) << "-1" << "\n";
        outputFile.close();
    }

private:
    Kratos::ModelPart &mrOutputModelPart;
    std::string mOutputFileName;
};
}

#endif //KRATOSMULTIPHYSICS_UNV_WRITER_H_INCLUDED
