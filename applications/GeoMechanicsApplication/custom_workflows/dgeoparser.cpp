// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Jonathan Nuttall
//

#pragma once

namespace Kratos
{
    void KratosGeoParser::parseMaterial(Model& model, std::string filepath)
    {
        std::string parameters = "{ \"Parameters\" : { \"materials_filename\" :\"" + filepath + "\"}}";
        Parameters material_file{ parameters };
        ReadMaterialsUtility(material_file, model);
    }

    void KratosGeoParser::parseMesh(ModelPart& model_part, std::string filepath)
    {
        // Parses MDPA file into model_part
        std::ifstream input(filepath);
        bool read_properties = false;
        bool read_nodes = false;
        bool read_elements = false;

        bool read_subparts = false;
        bool read_subparts_table = false;
        bool read_subparts_nodes = false;
        bool read_subparts_elements = false;
        bool read_subparts_conditions = false;

        std::string element_type;
        std::string part_name;
        std::string nodeStr;

        for (std::string line; getline(input, line);)
        {

            //===================== Properties =========================
            if (line.substr(0, 16) == "Begin Properties")
            {
                read_properties = true;
                std::size_t found = line.find_last_of(" ");
                int property_id = stoi(line.substr(found + 1));
                model_part.CreateNewProperties(property_id);
                continue;
            }
            if (line == "End Properties")
            {
                read_properties = false;
                continue;
            }
            if (read_properties)
            {
                KRATOS_ERROR << "Reading Properties - Not Implemented " << std::endl;
            }
            //=====================   Nodes   =========================
            if (line == "Begin Nodes")
            {
                read_nodes = true;
                continue;
            }
            if (line == "End Nodes")
            {
                read_nodes = false;
                continue;
            }
            if (read_nodes)
            {
                std::istringstream iss(line);
                int nodeId;
                double x, y, z;
                iss >> nodeId >> x >> y >> z;
                model_part.CreateNewNode(nodeId, x, y, z);
            }
            //====================   Element   ==========================
            if (line.substr(0, 14) == "Begin Elements")
            {
                read_elements = true;
                std::size_t found = line.find_last_of(" ");
                element_type = line.substr(found + 1);

                std::size_t posD = element_type.find_last_of("D");
                nodeStr = element_type.substr(posD + 1);

                continue;
            }
            if (line == "End Elements")
            {
                read_elements = false;
                continue;
            }

            if (read_elements)
            {
                unsigned long elementId, propertyId, node1, node2, node3;
                std::vector<ModelPart::IndexType> element_nodes;
                std::istringstream iss(line);

                if (nodeStr == "4N")
                {
                    unsigned long node4;
                    iss >> elementId >> propertyId >> node1 >> node2 >> node3 >> node4;
                    element_nodes = { node1, node2, node3, node4 };
                }
                else if (nodeStr == "3N")
                {
                    iss >> elementId >> propertyId >> node1 >> node2 >> node3;
                    element_nodes = { node1, node2, node3 };
                }
                else
                {
                    KRATOS_ERROR << "Element Type Unknown / Not Implemented " << std::endl;
                }
                auto p_elem_prop = model_part.pGetProperties(propertyId);
                model_part.CreateNewElement(element_type, elementId, element_nodes, p_elem_prop);
            }

            //===================== Properties =========================

            if (line.substr(0, 18) == "Begin SubModelPart")
            {
                read_subparts = true;
                std::size_t found = line.find_last_of(" ");
                part_name = line.substr(found + 1);
                model_part.CreateSubModelPart(part_name);
                continue;
            }
            if (line == "End SubModelPart")
            {
                read_subparts = false;
                continue;
            }
            if (read_subparts)
            {
                auto subpart = model_part.pGetSubModelPart(part_name);
                //===========  Sub-Tables  ===============
                if (line == "  Begin SubModelPartTables")
                {
                    read_subparts_table = true;
                    continue;
                }
                if (line == "  End SubModelPartTables")
                {
                    read_subparts_table = false;
                    continue;
                }
                if (read_subparts_table)
                {
                    KRATOS_ERROR << "Subpart Tables - Not Implemented " << std::endl;
                }

                //===========  Sub-Nodes  ===============

                if (line == "  Begin SubModelPartNodes")
                {
                    read_subparts_nodes = true;
                    continue;
                }
                if (line == "  End SubModelPartNodes")
                {
                    read_subparts_nodes = false;
                    continue;
                }
                if (read_subparts_nodes)
                {
                    auto node = model_part.pGetNode(stoi(line));
                    subpart->AddNode(node);
                }

                //===========  Sub-Elements  ===============

                if (line == "  Begin SubModelPartElements")
                {
                    read_subparts_elements = true;
                    continue;
                }
                if (line == "  End SubModelPartElements")
                {
                    read_subparts_elements = false;
                    continue;
                }
                if (read_subparts_elements)
                {
                    auto element = model_part.pGetElement(stoi(line));
                    subpart->AddElement(element);
                }

                //===========  Sub-Elements  ===============

                if (line == "  Begin SubModelPartConditions")
                {
                    read_subparts_conditions = true;
                    continue;
                }
                if (line == "  End SubModelPartConditions")
                {
                    read_subparts_conditions = false;
                    continue;
                }
                if (read_subparts_conditions)
                {
                    KRATOS_ERROR << "Subpart Conditions - Not Implemented " << std::endl;
                }
            }
        }
        input.close();
    }

    std::vector<std::shared_ptr<Process>> KratosGeoParser::parseProcess(ModelPart& model_part, Parameters projFile)
    {
        // Currently: In DGeoflow only fixed hydrostatic head has been , also need load of gravity.

        std::vector<std::shared_ptr<Process>> processes;

        auto constraints_processes = projFile["processes"]["constraints_process_list"];
        for (Parameters process : constraints_processes)
        {
            // we only support fixed hydrostatic head
            auto name = process["Parameters"]["model_part_name"].GetString();
            auto pressure_type = process["Parameters"]["fluid_pressure_type"].GetString();

            std::size_t found = name.find_last_of(".");
            std::string subname = name.substr(found + 1);

            ModelPart& part = model_part.GetSubModelPart(subname);

            if (pressure_type == "Uniform")
            {
                auto value = process["Parameters"]["value"].GetDouble();
                processes.push_back(make_shared<GeoFlowApplyConstantScalarValueProcess>(GeoFlowApplyConstantScalarValueProcess(part, WATER_PRESSURE,
                    value, 0, GeoFlowApplyConstantScalarValueProcess::VARIABLE_IS_FIXED)));
            }
            else if (pressure_type == "Hydrostatic")
            {
                auto cProcesses = process.Clone();
                cProcesses["Parameters"].RemoveValue("fluid_pressure_type");
                processes.push_back(make_shared<GeoFlowApplyConstantHydrostaticPressureProcess>(GeoFlowApplyConstantHydrostaticPressureProcess(part, cProcesses["Parameters"])));
            }
            else
            {
                KRATOS_ERROR << "Reading Processing - Not Implemented - Pressure_type" << std::endl;
            }
        }

        auto loads_processes = projFile["processes"]["loads_process_list"];
        // Should only have one.
        auto name = loads_processes.GetArrayItem(0)["Parameters"]["model_part_name"].GetString();
        std::size_t found = name.find_last_of(".");
        std::string subname = name.substr(found + 1);
        ModelPart& part = model_part.GetSubModelPart(subname);
        processes.push_back(make_shared<ApplyConstantScalarValueProcess>(ApplyConstantScalarValueProcess(part, VOLUME_ACCELERATION_X,
            0.0, 0, ApplyConstantScalarValueProcess::VARIABLE_IS_FIXED)));

        processes.push_back(make_shared<ApplyConstantScalarValueProcess>(ApplyConstantScalarValueProcess(part, VOLUME_ACCELERATION_Y, -9.81,
            0, ApplyConstantScalarValueProcess::VARIABLE_IS_FIXED)));

        processes.push_back(make_shared<Process>(ApplyConstantScalarValueProcess(part, VOLUME_ACCELERATION_Z, 0.0,
            0, ApplyConstantScalarValueProcess::VARIABLE_IS_FIXED)));

        return processes;
    }
}
