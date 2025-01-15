//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:   Andrea Gorgi

// Project includes
#include "iga_contact_process.h"
#include "geometries/nurbs_volume_geometry.h"

#include "utilities/variable_utils.h"

namespace Kratos
{

    IgaContactProcess::IgaContactProcess(
        Model& rModel, Parameters ThisParameters) : mpModel(&rModel), mParameters(ThisParameters)
    {
        
        mEchoLevel = mParameters["echo_level"].GetInt();

        // KRATOS_ERROR_IF_NOT(mParameters.Has("contact_model_part_name"))
        //     << "Missing \"contact_model_part_name\" section" << std::endl;

        // ModelPart& contact_model_part = mpModel->HasModelPart(mParameters["contact_model_part_name"].GetString())
        //                             ? mpModel->GetModelPart(mParameters["contact_model_part_name"].GetString())
        //                             : mpModel->CreateModelPart(mParameters["contact_model_part_name"].GetString());


        // std::string contact_model_part_name = "IgaModelPart.ContactInterface";

        // // Model& mainModel = r_model_part->GetModel();

        // ModelPart& contact_model_part = mpModel->HasModelPart(contact_model_part_name)
        //                               ? mpModel->GetModelPart(contact_model_part_name)
        //                               : mpModel->CreateModelPart(contact_model_part_name);


        // PHYSICS HERE OR IN IGA_MODELER???????????
        // const std::string& rDataFileName = mParameters.Has("physics_file_name")
        //     ? mParameters["physics_file_name"].GetString()
        //     : "physics.iga.json";

        //-----------------------------------------------------------------------------------------------
        // Obtain SLAVE interface b_reps
        const std::string slave_model_part_name = mParameters["contact_parameters"]["slave_model_part"]["sub_model_part_name"].GetString();

        KRATOS_ERROR_IF_NOT(mpModel->HasModelPart(slave_model_part_name)) << "ERROR: SLAVE MODEL PART " 
                                                << slave_model_part_name << "NOT CREATED" << std::endl; 

        ModelPart& slave_model_part = mpModel->GetModelPart(slave_model_part_name);
        GeometriesArrayType geometry_list_slave;
        GetCadGeometryList(geometry_list_slave, slave_model_part, mParameters["contact_parameters"]["slave_model_part"]);

        const IndexType slave_property_id = mParameters["contact_parameters"]["slave_model_part"]["property_id"].GetInt();
        mpPropSlave = slave_model_part.pGetProperties(slave_property_id);
    
        // Obtain MASTER interface b_reps
        const std::string master_model_part_name = mParameters["contact_parameters"]["master_model_part"]["sub_model_part_name"].GetString();

        KRATOS_ERROR_IF_NOT(mpModel->HasModelPart(master_model_part_name)) << "ERROR: MASTER MODEL PART " 
                                                << master_model_part_name << "NOT CREATED" << std::endl; 

        ModelPart& master_model_part = mpModel->GetModelPart(master_model_part_name);
        GeometriesArrayType geometry_list_master;
        GetCadGeometryList(geometry_list_master, master_model_part, mParameters["contact_parameters"]["master_model_part"]);

        const IndexType master_property_id = mParameters["contact_parameters"]["master_model_part"]["property_id"].GetInt();
        mpPropMaster = master_model_part.pGetProperties(master_property_id);

        mpCouplingGeometry = Kratos::make_shared<NurbsCouplingGeometry2D<PointType, PointerVector<NodeType>>>(geometry_list_slave, geometry_list_master);
   }

    void IgaContactProcess::Execute(){

        std::string main_model_part_name = "IgaModelPart";

        std::string contact_sub_model_part_name = "ContactInterface";

        // Model& mainModel = r_model_part->GetModel();

        ModelPart& main_model_part = mpModel->HasModelPart(main_model_part_name)
                                      ? mpModel->GetModelPart(main_model_part_name)
                                      : mpModel->CreateModelPart(main_model_part_name);

        // if (main_model_part.HasSubModelPart(contact_sub_model_part_name)) main_model_part.RemoveSubModelPart(contact_sub_model_part_name);

        
        ModelPart& contact_sub_model_part = main_model_part.HasSubModelPart(contact_sub_model_part_name)
                                            ? main_model_part.GetSubModelPart(contact_sub_model_part_name)
                                            : main_model_part.CreateSubModelPart(contact_sub_model_part_name);


        // KRATOS_WATCH(main_model_part.Conditions().size())
        // KRATOS_WATCH(contact_sub_model_part.Conditions().size())
        
        ConditionsArrayType& r_conditions_array = contact_sub_model_part.Conditions();
        KRATOS_TRACE_IF("Empty model part", r_conditions_array.size() == 0) << "YOUR CONTACT MODEL PART IS EMPTY" << std::endl;
        block_for_each(r_conditions_array, [&](Condition& rCond) {
                rCond.Set(TO_ERASE, true);
        });

        contact_sub_model_part.RemoveConditionsFromAllLevels(TO_ERASE);

        // KRATOS_WATCH(contact_sub_model_part.Conditions().size())
        // KRATOS_WATCH(main_model_part.Conditions().size())

        //---------------------------------------------------------------
        KRATOS_WATCH("ÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇ")
        GeometriesArrayType geometries;
        SizeType shape_function_derivatives_order = 1;
        if (mParameters.Has("shape_function_derivatives_order")) {
            shape_function_derivatives_order = mParameters["shape_function_derivatives_order"].GetInt();
        }
        else {
            KRATOS_INFO_IF("CreateQuadraturePointGeometries", mEchoLevel > 4)
                << "shape_function_derivatives_order is not provided and thus being considered as 1. " << std::endl;
        }

        std::string quadrature_method = mParameters.Has("quadrature_method")
            ? mParameters["integration_rule"].GetString()
            : "GAUSS";
        IntegrationInfo integration_info = mpCouplingGeometry->GetDefaultIntegrationInfo();
        for (IndexType i = 0; i < integration_info.LocalSpaceDimension(); ++i) {
            if (quadrature_method == "GAUSS") {
                integration_info.SetQuadratureMethod(0, IntegrationInfo::QuadratureMethod::GAUSS);
            }
            else if (quadrature_method == "EXTENDED_GAUSS") {
                integration_info.SetQuadratureMethod(0, IntegrationInfo::QuadratureMethod::EXTENDED_GAUSS);
            }
            else if (quadrature_method == "GRID") {
                integration_info.SetQuadratureMethod(0, IntegrationInfo::QuadratureMethod::GRID);
            }
            else {
                KRATOS_INFO("CreateQuadraturePointGeometries") << "Quadrature method: " << quadrature_method
                    << " is not available. Available options are \"GAUSS\" and \"GRID\". Default quadrature method is being considered." << std::endl;
            }
        }

        mpCouplingGeometry->CreateQuadraturePointGeometries(geometries, shape_function_derivatives_order, integration_info);

        SizeType id = 1;
        // if (contact_model_part.GetRootModelPart().Conditions().size() > 0)
        //     id = contact_model_part.GetRootModelPart().Conditions().back().Id() + 1;
        if (contact_sub_model_part.GetRootModelPart().Conditions().size() > 0)
            id = contact_sub_model_part.GetRootModelPart().Conditions().back().Id() + 1;
        
        KRATOS_ERROR_IF_NOT(mParameters.Has("name"))
            << "\"name\" need to be specified." << std::endl;
        std::string name = mParameters["name"].GetString();


        this->CreateConditions(
                        geometries.ptr_begin(), geometries.ptr_end(),
                        contact_sub_model_part, name, id, mpPropMaster, mpPropSlave);

    }



     void IgaContactProcess::CreateConditions(
        typename GeometriesArrayType::ptr_iterator rGeometriesBegin,
        typename GeometriesArrayType::ptr_iterator rGeometriesEnd,
        ModelPart& rModelPart,
        std::string& rConditionName,
        SizeType& rIdCounter,
        PropertiesPointerType pPropMaster,
        PropertiesPointerType pPropSlave) const
    {
        // const SupportContact2DCondition rReferenceCondition = SupportContact2DCondition();
        // const Condition& rReferenceCondition = KratosComponents<Condition>::Get(rConditionName);

        ModelPart::ConditionsContainerType new_condition_list;

        KRATOS_INFO_IF("CreateConditions", mEchoLevel > 2)
            << "Creating conditions of type " << rConditionName
            << " in " << rModelPart.Name() << "-SubModelPart." << std::endl;

        for (auto it = rGeometriesBegin; it != rGeometriesEnd; ++it)
        {
            new_condition_list.push_back(
                Kratos::make_intrusive<SupportContact2DCondition>(
                rIdCounter, (*it), pPropMaster, pPropSlave));
            for (SizeType i = 0; i < (*it)->size(); ++i) {
                // rModelPart.AddNode((*it)->pGetPoint(i));
                rModelPart.Nodes().push_back((*it)->pGetPoint(i));
            }
            rIdCounter++;
        }

        // INITIALIZE THE CONDITIONS
        // const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // for (int i = 0; i < new_condition_list.size(); i++) 
        // {
        //     new_condition_list[i].Initialize(r_current_process_info);
        // }
        
        rModelPart.AddConditions(new_condition_list.begin(), new_condition_list.end());

        EntitiesUtilities::InitializeEntities<Condition>(rModelPart);

    }

    ///@name CAD functionalities
    ///@{

    void IgaContactProcess::GetCadGeometryList(
        GeometriesArrayType& rGeometryList,
        ModelPart& rModelPart,
        const Parameters rParameters) const
    {

        static int starting_brep_ids;

        // get surrogate model part name
        std::string surrogate_model_part_name;
        if (!mParameters.Has("surrogate_model_part_name")) surrogate_model_part_name = "surrogate_model_part";
        else {
            surrogate_model_part_name = mParameters["surrogate_model_part_name"].GetString();
        }

        if (rParameters.Has("brep_id")) {
            rGeometryList.push_back(rModelPart.pGetGeometry(rParameters["brep_id"].GetInt()));
        }
        if (rParameters.Has("brep_ids")) {
            for (SizeType i = 0; i < rParameters["brep_ids"].size(); ++i) {
                rGeometryList.push_back(rModelPart.pGetGeometry(rParameters["brep_ids"][i].GetInt()));
            }
            // int lastIndex
            starting_brep_ids = rParameters["brep_ids"][rParameters["brep_ids"].size()-1].GetInt() + 1;

            // SBM THINGS
            // // OUTER
            // std::string conditionName = rParameters["iga_model_part"].GetString();
            // if (conditionName.rfind("SBM", 0) == 0) { 
            //     ModelPart& surrogateModelPart_outer = mpModel->GetModelPart(surrogate_model_part_name +"_outer");
            //     if (surrogateModelPart_outer.Conditions().size() > 0) {
            //         // 2D
            //         if (surrogateModelPart_outer.GetCondition(1).GetGeometry().size() == 2) {
            //             int sizeSurrogateLoop_outer = surrogateModelPart_outer.Nodes().size();
            //             for (int j = 0; j < (sizeSurrogateLoop_outer-1); ++j) {
            //                 // Add the brep_ids of the internal boundary for SBMLaplacianCondition
            //                 rGeometryList.push_back(rModelPart.pGetGeometry(starting_brep_ids));
            //                 starting_brep_ids++;
            //             }
                        
            //         } else { // 3D
            //             int sizeSurrogateLoop = surrogateModelPart_outer.Conditions().size(); //lastSurrogateNode.Id() - firstSurrogateNode.Id() + 1 ;
            //             for (SizeType j = 0; j < sizeSurrogateLoop-1; ++j) {
            //                 // Add the brep_ids of the internal boundary for SBMLaplacianCondition
            //                 rGeometryList.push_back(rModelPart.pGetGeometry(starting_brep_ids));
            //                 starting_brep_ids++;
            //             }
            //         }     
            //     }
            // }
        }
        // else { // SBM THINGS
        //     // INNER   
        //     ModelPart& surrogateModelPart_inner = mpModel->GetModelPart(surrogate_model_part_name + "_inner");
        //     if (surrogateModelPart_inner.Elements().size() > 0) {

        //         int sizeSurrogateLoop = surrogateModelPart_inner.Nodes().size();

        //         for (int iel = 1; iel < surrogateModelPart_inner.Elements().size()+1; iel++) {
        //             // Each element in the surrogate_model_part represents a surrogate boundary loop. First "node" is the initial ID of the first surrogate node and
        //             // the second "node" is the last surrogate node of that loop. (We have done this in the case we have multiple surrogate boundaries and 1 model part)
        //             Node& firstSurrogateNode = surrogateModelPart_inner.pGetElement(iel)->GetGeometry()[0]; // Element 1 because is the only surrogate loop
        //             Node& lastSurrogateNode = surrogateModelPart_inner.pGetElement(iel)->GetGeometry()[1];  // Element 1 because is the only surrogate loop
        //             int sizeSurrogateLoop = lastSurrogateNode.Id() - firstSurrogateNode.Id() + 1 ;

        //             for (SizeType j = 0; j < sizeSurrogateLoop; ++j) {
        //                 // Add the brep_ids of the internal boundary for SBMLaplacianCondition
        //                 rGeometryList.push_back(rModelPart.pGetGeometry(starting_brep_ids));
        //                 starting_brep_ids++;
        //             }
        //         }
        //     } else { // 3D
        //         int sizeSurrogateLoop = surrogateModelPart_inner.Conditions().size(); //lastSurrogateNode.Id() - firstSurrogateNode.Id() + 1 ;
        //         for (SizeType j = 0; j < sizeSurrogateLoop; ++j) {
        //             // Add the brep_ids of the internal boundary for SBMLaplacianCondition
        //             rGeometryList.push_back(rModelPart.pGetGeometry(starting_brep_ids));
        //             starting_brep_ids++;
        //         }
        //     }
        // }
        if (rParameters.Has("brep_name")) {
            rGeometryList.push_back(rModelPart.pGetGeometry(rParameters["brep_name"].GetString()));
        }
        if (rParameters.Has("brep_names")) {
            for (SizeType i = 0; i < rParameters["brep_names"].size(); ++i) {
                rGeometryList.push_back(rModelPart.pGetGeometry(rParameters["brep_names"][i].GetString()));
            }
        }
        KRATOS_ERROR_IF(rGeometryList.size() == 0)
            << "Empty geometry list. Either \"brep_id\", \"brep_ids\", \"brep_name\" or \"brep_names\" are the possible options." << std::endl;
    }

    ///@}

    /// Reads in a json formatted file and returns its KratosParameters instance.
    Parameters IgaContactProcess::ReadParamatersFile(
        const std::string& rDataFileName) const
    {
        // Check if rDataFileName ends with ".cad.json" and add it if needed.
        const std::string data_file_name = (rDataFileName.compare(rDataFileName.size() - 9, 9, ".iga.json") != 0)
            ? rDataFileName + ".iga.json"
            : rDataFileName;

        std::ifstream infile(data_file_name);
        KRATOS_ERROR_IF_NOT(infile.good()) << "Physics fil: "
            << data_file_name << " cannot be found." << std::endl;
        KRATOS_INFO_IF("ReadParamatersFile", mEchoLevel > 3)
            << "Reading file: \"" << data_file_name << "\"" << std::endl;

        std::stringstream buffer;
        buffer << infile.rdbuf();

        return Parameters(buffer.str());
    }
} // End namespace Kratos
