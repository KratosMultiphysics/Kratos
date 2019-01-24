#ifndef KRATOSMULTIPHYSICS_UNV_WRITER_H_INCLUDED
#define KRATOSMULTIPHYSICS_UNV_WRITER_H_INCLUDED

#include <iostream>
#include <fstream>
#include <includes/model_part.h>
#include "includes/exception.h"

namespace Kratos {


class UnvWriter {
public:
    
    enum class DatasetLocation {
        DATA_AT_NODES               = 1,
        DATA_AT_ELEMENTS            = 2,
        DATA_AT_NODES_ON_ELEMENTS   = 3,
        DATA_AT_POINTS              = 5,
        DATA_ON_ELEMENTS_AT_NODES   = 6
    };

    enum class ModelType {
        UNKNOWN         = 0,
        STRUCTURAL      = 1,
        HEAT_TRANSFER   = 2,
        FLUID_FLOW      = 3
    };

    enum class AnalysisType {
        UNKNOWN                             = 0,
        STATIC                              = 1,
        NORMAL_MODE                         = 2,
        COMPLEX_EIGENVALUE_FIRST_ORDER      = 3,
        TRANSIENT                           = 4,
        FREQUENCY_RESPONSE                  = 5,
        BUCKLING                            = 6,
        COMPLEX_EIGENVALUE_SECOND_ORDER     = 7,
        STATIC_NON_LINEAR                   = 9,
        CRAIG_BAMPTON_CONSTRAINT_MODES      = 10,
        EQUIVALENT_ATTACHMENT_MODES         = 11,
        EFFECTIVE_MASS_MODES                = 12,
        EFFECTIVE_MASS_MATRIX               = 13,
        EFFECTIVE_MASS_MATRIX_COPY          = 14,   // This record is duplicared intentionally
        DISTRIBUTED_LOAD_LOAD_DISTRIBUTION  = 15,
        DISTRIBUTED_LOAD_ATTACHMENT_MODES   = 16
    };

    // 3_DOF_* Records have their name changed to D3_DOF to avoid conflicts.
    enum class DataCharacteristics {
        UNKNOWN = 0,
        SCALAR = 1,
        D3_DOF_GLOBAL_TRANSLATION_VECTOR = 2,
        D3_DOF_GLOBAL_TRANSLATION_ROTATION_VECTOR = 3,
        SYMMETRIC_GLOBAL_TENSOR = 4,
        STRESS_RESULTANTS = 6
    };

    // This record is to big, I will consider moving these enums to another file.
    // enum class ResultType {

    // };

    enum class DataType {
        INTEGER = 1,
        SINGLE_PRECISION_FLOATING_POINT = 2,
        DOUBLE_PRECISION_FLOATING_POINT = 4,
        SINGLE_PRECISION_COMPLEX = 5,
        DOUBLE_PRECISION_COMPLEX = 6
    };

    UnvWriter(Kratos::ModelPart &modelPart, const std::string &outFileWithoutExtension)
            : mrOutputModelPart(modelPart),
              mOutputFileName(outFileWithoutExtension + ".unv") {
    }


    void writeMesh() {
        initializeOutputFile();
        writeNodes();
        writeElements();
        writeNodalResults();
    }

    template <typename Enumeration>
    auto as_integer(Enumeration const value)
        -> typename std::underlying_type<Enumeration>::type
    {
        return static_cast<typename std::underlying_type<Enumeration>::type>(value);
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

    // Partially extracted from: http://users.ices.utexas.edu
    // # beginning of dataset
    // # type of dataset: data at mesh entities
    // # R.  1: unique number of dataset (dataset_label)
    // # R.  2: text describing content (dataset_name)
    // # R.  3: data belongs to: nodes, elements,...
    // #        (dataset_location)
    // # R.  4: user-specified text (id_lines_1_to_5[0])
    // # R.  5: user-specified text (id_lines_1_to_5[1])
    // # R.  6: user-specified text (id_lines_1_to_5[2])
    // # R.  7: user-specified text (id_lines_1_to_5[3])
    // # R.  8: user-specified text (id_lines_1_to_5[4])
    // # R.  9: (model_type) (analysis_type)
    // #        (data_characteristic) (result_type)
    // #        (data_type) (nvaldc)
    // # R. 10: analysis-specific data (record_10)
    // # R. 11: analysis-specific data (record_11)
    // # R. 12: analysis-specific data (record_12)
    // # R. 13: analysis-specific data (record_13)

    // Fordes

    void writeNodalResults() {
        std::ofstream outputFile;
        outputFile.open(mOutputFileName, std::ios::out | std::ios::app);

        const int dataSetNumberForResults = 2414;
        std::string dataSetLabel = 1;
        std::string dataSetName = "NodalResults";   
        const int physicalPropertyTableNumber = 1;
        const int materialPropertyTableNumber = 1;
        const int color = 0;

        outputFile << std::setw(6) << "-1" << "\n";                                                 // Begin block
        outputFile << std::setw(6) << dataSetNumberForResults << "\n";                              // DatasetID

        outputFile << std::setw(10) << dataSetLabel << "\n";                                        // Record 1
        outputFile << std::setw(6) << dataSetName << "\n";                                          // Record 2
        outputFile << std::setw(10) << as_integer(DatasetLocation::DATA_AT_NODES) << "\n";          // Record 3

        outputFile << std::setw(6) << dataSetName << "\n";                                          // Record 4
        outputFile << std::setw(6) << 'None' << "\n";                                               // Record 5
        outputFile << std::setw(6) << 'NONE' << "\n";                                               // Record 6
        outputFile << std::setw(6) << 'NONE' << "\n";                                               // Record 7
        outputFile << std::setw(6) << 'NONE' << "\n";                                               // Record 8
        
        // ModelType, AnalysisType, DataCharacteristic, ResultType, DataType, NumberOfDataValues    // Record 9
        outputFile << std::setw(6);
        outputFile << as_integer(ModelType::STRUCTURAL); 
        outputFile << as_integer(AnalysisType::STATIC);
        outputFile << as_integer(DataCharacteristics::SCALAR);
        outputFile << 5;
        outputFile << as_integer(DataType::SINGLE_PRECISION_FLOATING_POINT);
        outputFile << 1; 
        outputFile << "\n";

        // ????
        outputFile << std::setw(6) << 0 << 0 << 0 << 0 << 0 << 1 << 0 << 0 << "\n";                 // Record 10
        outputFile << std::setw(6) << 0 << 0 << "\n";                                               // Record 11
        outputFile << std::setw(6) << 0 << "\n";                                                    // Record 12
        outputFile << std::setw(6) << 0 << "\n";                                                    // Record 13

        // Data at nodes:
        for (auto &node_i : mrOutputModelPart.Nodes()) {
            int node_label = node_i.Id();
            outputFile << std::setw(6) << node_label << "\n";                                       // Record 14 - Node Number
            outputFile << std::setw(6) << node_i.FastGetSolutionStepValue(TEMPERATURE) << "\n";     // Record 15 - NumberOfDataValues' data of the node
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
