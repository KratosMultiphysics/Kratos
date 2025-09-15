//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos Roig
//
//

/* System includes */
#include <iostream>
#include <fstream>

/* External includes */

/* Project includes */
#include "includes/model_part.h"
#include "includes/exception.h"

#include "input_output/unv_output.h"

namespace Kratos {

UnvOutput::UnvOutput(Kratos::ModelPart &modelPart, const std::string &outFileWithoutExtension)
  : mrOutputModelPart(modelPart),
    mOutputFileName(outFileWithoutExtension + ".unv") {
}


void UnvOutput::WriteMesh() {
    WriteNodes();
    WriteElements();
}

void UnvOutput::InitializeOutputFile() {
    std::ofstream outputFile;
    outputFile.open(mOutputFileName, std::ios::out | std::ios::trunc);
    outputFile.close();
}

void UnvOutput::WriteNodes() {
    std::ofstream outputFile;
    outputFile.open(mOutputFileName, std::ios::out | std::ios::app);

    outputFile << std::scientific;
    outputFile << std::setprecision(15);

    const int exportCoordinateSystemNumber = 0;
    const int displacementCoordinateSystemNumber = 0;
    const int color = 0;

    outputFile << std::setw(6) << "-1" << "\n";
    outputFile << std::setw(6) << as_integer(DatasetID::NODES_DATASET) << "\n";

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

void UnvOutput::WriteElements() {
    std::ofstream outputFile;
    outputFile.open(mOutputFileName, std::ios::out | std::ios::app);

    const int physicalPropertyTableNumber = 1;
    const int materialPropertyTableNumber = 1;
    const int color = 0;

    outputFile << std::setw(6) << "-1" << "\n";
    outputFile << std::setw(6) << as_integer(DatasetID::ELEMENTS_DATASET) << "\n";

    for (auto &element : mrOutputModelPart.Elements()) {
        const int elementLabel = element.Id();
        Kratos::ModelPart::ConditionType::GeometryType elementGeometry = element.GetGeometry();
        // Write triangles
        if (elementGeometry.size() == 3 && elementGeometry.WorkingSpaceDimension() == 2) {
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
        else if (elementGeometry.size() == 4 && elementGeometry.WorkingSpaceDimension() == 3) {
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

void UnvOutput::WriteNodalResults(const Variable<bool>& rVariable, const double timeStep) {
    WriteNodalResultRecords(rVariable, 1, timeStep);
}

void UnvOutput::WriteNodalResults(const Variable<int>& rVariable, const double timeStep) {
    WriteNodalResultRecords(rVariable, 1, timeStep);
}   

void UnvOutput::WriteNodalResults(const Variable<double>& rVariable, const double timeStep) {
    WriteNodalResultRecords(rVariable, 1, timeStep);
} 

void UnvOutput::WriteNodalResults(const Variable<array_1d<double,3>>& rVariable, const double timeStep) {
    WriteNodalResultRecords(rVariable, 3, timeStep);
}

void UnvOutput::WriteNodalResults(const Variable<Vector>& rVariable, const double timeStep) {
    KRATOS_ERROR << "Dynamic Vector results are not yet supported in UNV" << std::endl;
}

void UnvOutput::WriteNodalResults(const Variable<Matrix>& rVariable, const double timeStep) {
    KRATOS_ERROR << "Matrix results are not yet supported in UNV" << std::endl;
}

UnvOutput::DataCharacteristics UnvOutput::GetDataType(const Variable<bool>& rVariable) {
    return UnvOutput::DataCharacteristics::SCALAR;
}

UnvOutput::DataCharacteristics UnvOutput::GetDataType(const Variable<int>& rVariable) {
    return UnvOutput::DataCharacteristics::SCALAR;
}   

UnvOutput::DataCharacteristics UnvOutput::GetDataType(const Variable<double>& rVariable) {
    return UnvOutput::DataCharacteristics::SCALAR;
} 

UnvOutput::DataCharacteristics UnvOutput::GetDataType(const Variable<array_1d<double,3>>& rVariable) {
    return UnvOutput::DataCharacteristics::D3_DOF_GLOBAL_TRANSLATION_VECTOR;
}

UnvOutput::DataCharacteristics UnvOutput::GetDataType(const Variable<Vector>& rVariable) {
    KRATOS_ERROR << "Dynamic Vector results are not yet supported in UNV" << std::endl;
} 

UnvOutput::DataCharacteristics UnvOutput::GetDataType(const Variable<Matrix>& rVariable) {
    KRATOS_ERROR << "Matrix results are not yet supported in UNV" << std::endl;
}

void UnvOutput::WriteNodalResultValues(std::ofstream &outputFile, const Node& node, const Variable<bool>& rVariable) {
    outputFile << std::setw(13) << node.FastGetSolutionStepValue(rVariable) << "\n";
}

void UnvOutput::WriteNodalResultValues(std::ofstream &outputFile, const Node& node, const Variable<int>& rVariable) {
    outputFile << std::setw(13) << node.FastGetSolutionStepValue(rVariable) << "\n";
}

void UnvOutput::WriteNodalResultValues(std::ofstream &outputFile, const Node& node, const Variable<double>& rVariable) {
    outputFile << std::setw(13) << node.FastGetSolutionStepValue(rVariable) << "\n";
}

void UnvOutput::WriteNodalResultValues(std::ofstream &outputFile, const Node& node, const Variable<array_1d<double,3>>& rVariable) {
    auto & temp = node.FastGetSolutionStepValue(rVariable);

    outputFile << std::setw(13) << temp[0];
    outputFile << std::setw(13) << temp[1];
    outputFile << std::setw(13) << temp[2];
    outputFile << "\n";
}

void UnvOutput::WriteNodalResultValues(std::ofstream &outputFile, const Node& node, const Variable<Vector>& rVariable) {
    KRATOS_ERROR << "Dynamic Vector results are not yet supported by in UNV" << std::endl;
}

void UnvOutput::WriteNodalResultValues(std::ofstream &outputFile, const Node& node, const Variable<Matrix>& rVariable) {
    KRATOS_ERROR << "Matrix results are not yet supported by in UNV" << std::endl;
}

} // namespace Kratos
