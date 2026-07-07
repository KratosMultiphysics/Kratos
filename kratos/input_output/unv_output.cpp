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

// System includes

// External includes

// Project includes
#include "input_output/unv_output.h"

namespace Kratos {

UnvOutput::UnvOutput(
    Kratos::ModelPart& rModelPart, 
    const std::string& rOutFileWithoutExtension
    )
  : mrOutputModelPart(rModelPart),
    mOutputFileName(rOutFileWithoutExtension + ".unv") {
}

void UnvOutput::WriteMesh() {
    // Write the nodes
    WriteNodes();

    // Write the geometry (elements or conditions) based on availability
    if (mrOutputModelPart.Elements().size() > 0) {
        // Write the elements if they exist in the model part
        WriteElements();
    } else if (mrOutputModelPart.Conditions().size() > 0) {
        KRATOS_WARNING("UnvOutput") << "No elements found in the model part. Writing conditions instead." << std::endl;
        
        // Write the conditions if they exist in the model part
        WriteConditions();
    } else {
        KRATOS_WARNING("UnvOutput") << "No elements or conditions found in the model part. No mesh will be written." << std::endl;
    }
}

void UnvOutput::InitializeOutputFile() {
    std::ofstream rOutputFile;
    rOutputFile.open(mOutputFileName, std::ios::out | std::ios::trunc);
    rOutputFile.close();
}

void UnvOutput::WriteNodes() {
    std::ofstream rOutputFile;
    rOutputFile.open(mOutputFileName, std::ios::out | std::ios::app);

    rOutputFile << std::scientific;
    rOutputFile << std::setprecision(15);

    const int export_coordinate_system = 0;
    const int displacement_coordinate_system_number = 0;
    const int color = 0;

    rOutputFile << std::setw(6) << "-1" << "\n";
    rOutputFile << std::setw(6) << as_integer(DatasetID::NODES_DATASET) << "\n";

    int node_label;
    double x_coordinate, y_coordinate, z_coordinate;
    for (auto& r_node : mrOutputModelPart.Nodes()) {
        node_label = r_node.Id();
        x_coordinate = r_node.X();
        y_coordinate = r_node.Y();
        z_coordinate = r_node.Z();
        rOutputFile << std::setw(10) << node_label << std::setw(10) << export_coordinate_system
                    << std::setw(10)
                    << displacement_coordinate_system_number << std::setw(10) << color << "\n";
        rOutputFile << std::setw(25) << x_coordinate << std::setw(25) << y_coordinate << std::setw(25)
                    << z_coordinate << "\n";
    }
    rOutputFile << std::setw(6) << "-1" << "\n";

    rOutputFile.close();
}

void UnvOutput::WriteElements() {
    std::ofstream rOutputFile;
    rOutputFile.open(mOutputFileName, std::ios::out | std::ios::app);

    const int physical_property_table_number = 1;
    const int material_property_table_number = 1;
    const int color = 0;

    rOutputFile << std::setw(6) << "-1" << "\n";
    rOutputFile << std::setw(6) << as_integer(DatasetID::ELEMENTS_DATASET) << "\n";

    for (auto& r_element : mrOutputModelPart.Elements()) {
        const int element_label = r_element.Id();
        auto& r_geometry = r_element.GetGeometry();
        // Write triangles
        if (r_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle2D3) {
            const int fe_descriptor_id = 41; // Plane Stress Linear Triangle
            const int number_of_nodes = 3;
            rOutputFile << std::setw(10) << element_label;
            rOutputFile << std::setw(10) << fe_descriptor_id;
            rOutputFile << std::setw(10) << physical_property_table_number;
            rOutputFile << std::setw(10) << material_property_table_number;
            rOutputFile << std::setw(10) << color;
            rOutputFile << std::setw(10) << number_of_nodes << "\n";
            rOutputFile << std::setw(10) << r_geometry[0].Id();
            rOutputFile << std::setw(10) << r_geometry[1].Id();
            rOutputFile << std::setw(10) << r_geometry[2].Id() << "\n";
        } else if (r_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4) {
            const int fe_descriptor_id = 111; // Solid linear tetrahedron
            const int number_of_nodes = 4;
            rOutputFile << std::setw(10) << element_label;
            rOutputFile << std::setw(10) << fe_descriptor_id;
            rOutputFile << std::setw(10) << physical_property_table_number;
            rOutputFile << std::setw(10) << material_property_table_number;
            rOutputFile << std::setw(10) << color;
            rOutputFile << std::setw(10) << number_of_nodes << "\n";
            rOutputFile << std::setw(10) << r_geometry[0].Id();
            rOutputFile << std::setw(10) << r_geometry[1].Id();
            rOutputFile << std::setw(10) << r_geometry[2].Id();
            rOutputFile << std::setw(10) << r_geometry[3].Id() << "\n";
        } else {
            KRATOS_ERROR << "Element with ID " << element_label << " has unsupported geometry for UNV output." << std::endl;
        }
    }
    rOutputFile << std::setw(6) << "-1" << "\n";
    rOutputFile.close();
}

void UnvOutput::WriteConditions() {
    std::ofstream rOutputFile;
    rOutputFile.open(mOutputFileName, std::ios::out | std::ios::app);

    const int physical_property_table_number = 1;
    const int material_property_table_number = 1;
    const int color = 0;

    rOutputFile << std::setw(6) << "-1" << "\n";
    rOutputFile << std::setw(6) << as_integer(DatasetID::ELEMENTS_DATASET) << "\n";

    for (auto& r_condition : mrOutputModelPart.Conditions()) {
        const int element_label = r_condition.Id();
        auto& r_geometry = r_condition.GetGeometry();
        // Write lines and triangles
        if (r_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle3D3) {
            const int fe_descriptor_id = 21; // Linear beam
            const int number_of_nodes = 2;
            rOutputFile << std::setw(10) << element_label;
            rOutputFile << std::setw(10) << fe_descriptor_id;
            rOutputFile << std::setw(10) << physical_property_table_number;
            rOutputFile << std::setw(10) << material_property_table_number;
            rOutputFile << std::setw(10) << color;
            rOutputFile << std::setw(10) << number_of_nodes << "\n";
            rOutputFile << std::setw(10) << r_geometry[0].Id();
            rOutputFile << std::setw(10) << r_geometry[1].Id();
        } else if (r_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle3D3) {
            const int fe_descriptor_id = 91; // Thin Shell Linear Triangle
            const int number_of_nodes = 3;
            rOutputFile << std::setw(10) << element_label;
            rOutputFile << std::setw(10) << fe_descriptor_id;
            rOutputFile << std::setw(10) << physical_property_table_number;
            rOutputFile << std::setw(10) << material_property_table_number;
            rOutputFile << std::setw(10) << color;
            rOutputFile << std::setw(10) << number_of_nodes << "\n";
            rOutputFile << std::setw(10) << r_geometry[0].Id();
            rOutputFile << std::setw(10) << r_geometry[1].Id();
            rOutputFile << std::setw(10) << r_geometry[2].Id() << "\n";
        } else {
            KRATOS_ERROR << "Element with ID " << element_label << " has unsupported geometry for UNV output." << std::endl;
        }
    }
    rOutputFile << std::setw(6) << "-1" << "\n";
    rOutputFile.close();
}

void UnvOutput::WriteNodalResults(const Variable<bool>& rVariable, const double TimeStep) {
    WriteNodalResultRecords<Variable<bool>, WriteType::HISTORICAL>(rVariable, 1, TimeStep);
}

void UnvOutput::WriteNodalResults(const Variable<int>& rVariable, const double TimeStep) {
    WriteNodalResultRecords<Variable<int>, WriteType::HISTORICAL>(rVariable, 1, TimeStep);
}   

void UnvOutput::WriteNodalResults(const Variable<double>& rVariable, const double TimeStep) {
    WriteNodalResultRecords<Variable<double>, WriteType::HISTORICAL>(rVariable, 1, TimeStep);
} 

void UnvOutput::WriteNodalResults(const Variable<array_1d<double,3>>& rVariable, const double TimeStep) {
    WriteNodalResultRecords<Variable<array_1d<double,3>>, WriteType::HISTORICAL>(rVariable, 3, TimeStep);
}

void UnvOutput::WriteNodalResults(const Variable<Vector>& rVariable, const double TimeStep) {
    KRATOS_ERROR << "Dynamic Vector results are not yet supported in UNV" << std::endl;
}

void UnvOutput::WriteNodalResults(const Variable<Matrix>& rVariable, const double TimeStep) {
    KRATOS_ERROR << "Matrix results are not yet supported in UNV" << std::endl;
}

void UnvOutput::WriteNodalNonHistoricalResults(const Variable<bool>& rVariable, const double TimeStep) {
    WriteNodalResultRecords<Variable<bool>, WriteType::NON_HISTORICAL>(rVariable, 1, TimeStep);
}

void UnvOutput::WriteNodalNonHistoricalResults(const Variable<int>& rVariable, const double TimeStep) {
    WriteNodalResultRecords<Variable<int>, WriteType::NON_HISTORICAL>(rVariable, 1, TimeStep);
}   

void UnvOutput::WriteNodalNonHistoricalResults(const Variable<double>& rVariable, const double TimeStep) {
    WriteNodalResultRecords<Variable<double>, WriteType::NON_HISTORICAL>(rVariable, 1, TimeStep);
} 

void UnvOutput::WriteNodalNonHistoricalResults(const Variable<array_1d<double,3>>& rVariable, const double TimeStep) {
    WriteNodalResultRecords<Variable<array_1d<double,3>>, WriteType::NON_HISTORICAL>(rVariable, 3, TimeStep);
}

void UnvOutput::WriteNodalNonHistoricalResults(const Variable<Vector>& rVariable, const double TimeStep) {
    KRATOS_ERROR << "Dynamic Vector results are not yet supported in UNV" << std::endl;
}

void UnvOutput::WriteNodalNonHistoricalResults(const Variable<Matrix>& rVariable, const double TimeStep) {
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

void UnvOutput::WriteNodalResultValues(std::ofstream& rOutputFile, const Node& rNode, const Variable<bool>& rVariable) {
    rOutputFile << std::setw(13) << rNode.FastGetSolutionStepValue(rVariable) << "\n";
}

void UnvOutput::WriteNodalResultValues(std::ofstream& rOutputFile, const Node& rNode, const Variable<int>& rVariable) {
    rOutputFile << std::setw(13) << rNode.FastGetSolutionStepValue(rVariable) << "\n";
}

void UnvOutput::WriteNodalResultValues(std::ofstream& rOutputFile, const Node& rNode, const Variable<double>& rVariable) {
    rOutputFile << std::setw(13) << rNode.FastGetSolutionStepValue(rVariable) << "\n";
}

void UnvOutput::WriteNodalResultValues(std::ofstream& rOutputFile, const Node& rNode, const Variable<array_1d<double,3>>& rVariable) {
    const auto& r_temp = rNode.FastGetSolutionStepValue(rVariable);

    rOutputFile << std::setw(13) << r_temp[0];
    rOutputFile << std::setw(13) << r_temp[1];
    rOutputFile << std::setw(13) << r_temp[2];
    rOutputFile << "\n";
}

void UnvOutput::WriteNodalResultValues(std::ofstream& rOutputFile, const Node& rNode, const Variable<Vector>& rVariable) {
    KRATOS_ERROR << "Dynamic Vector results are not yet supported by in UNV" << std::endl;
}

void UnvOutput::WriteNodalResultValues(std::ofstream& rOutputFile, const Node& rNode, const Variable<Matrix>& rVariable) {
    KRATOS_ERROR << "Matrix results are not yet supported by in UNV" << std::endl;
}

void UnvOutput::WriteNodalNonHistoricalResultValues(std::ofstream& rOutputFile, const Node& rNode, const Variable<bool>& rVariable) {
    rOutputFile << std::setw(13) << rNode.GetValue(rVariable) << "\n";
}

void UnvOutput::WriteNodalNonHistoricalResultValues(std::ofstream& rOutputFile, const Node& rNode, const Variable<int>& rVariable) {
    rOutputFile << std::setw(13) << rNode.GetValue(rVariable) << "\n";
}

void UnvOutput::WriteNodalNonHistoricalResultValues(std::ofstream& rOutputFile, const Node& rNode, const Variable<double>& rVariable) {
    rOutputFile << std::setw(13) << rNode.GetValue(rVariable) << "\n";
}

void UnvOutput::WriteNodalNonHistoricalResultValues(std::ofstream& rOutputFile, const Node& rNode, const Variable<array_1d<double,3>>& rVariable) {
    const auto& r_temp = rNode.GetValue(rVariable);

    rOutputFile << std::setw(13) << r_temp[0];
    rOutputFile << std::setw(13) << r_temp[1];
    rOutputFile << std::setw(13) << r_temp[2];
    rOutputFile << "\n";
}

void UnvOutput::WriteNodalNonHistoricalResultValues(std::ofstream& rOutputFile, const Node& rNode, const Variable<Vector>& rVariable) {
    KRATOS_ERROR << "Dynamic Vector results are not yet supported by in UNV" << std::endl;
}

void UnvOutput::WriteNodalNonHistoricalResultValues(std::ofstream& rOutputFile, const Node& rNode, const Variable<Matrix>& rVariable) {
    KRATOS_ERROR << "Matrix results are not yet supported by in UNV" << std::endl;
}

} // namespace Kratos
