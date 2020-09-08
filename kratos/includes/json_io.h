//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//











#if !defined(KRATOS_JSON_IO_H_INCLUDED )
#define  KRATOS_JSON_IO_H_INCLUDED



// System includes
#include <string>
// #include <fstream>
#include <cstdio>
#include <iostream>
#include <set>
#include <map>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/io.h"
#include "utilities/timer.h"
#include "containers/flags.h"
#include "includes/mesh.h"

#include "rapidjson/filereadstream.h"
#include "rapidjson/document.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// An IO class for reading and writing a modelpart
/** This class reads and writes all modelpart data including the meshes.
*/
class KratosJsonIO : public IO
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosJsonIO
    KRATOS_CLASS_POINTER_DEFINITION(KratosJsonIO);

    typedef IO BaseType;

    typedef BaseType::NodeType NodeType;

    typedef BaseType::MeshType MeshType;

    typedef BaseType::NodesContainerType NodesContainerType;

    typedef BaseType::PropertiesContainerType PropertiesContainerType;

    typedef BaseType::ElementsContainerType ElementsContainerType;

    typedef BaseType::ConditionsContainerType ConditionsContainerType;

    typedef BaseType::ConnectivitiesContainerType ConnectivitiesContainerType;

//     typedef std::vector<std::ofstream*> OutputFilesContainerType;

    typedef std::size_t SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with  filenames.
    KratosJsonIO(
        std::string const& Filename,
        const Flags Options = IO::READ | IO::IGNORE_VARIABLES_ERROR.AsFalse())
        :IO()
    {
        mFilename = Filename;
        mOptions = Options;
    }

    virtual void ReadModelPart(ModelPart & rThisModelPart)
    {
        //read mFilename
        FILE* fp = fopen(mFilename.c_str(), "rb"); // non-Windows use "r"
        char readBuffer[65536];
        rapidjson::FileReadStream is(fp, readBuffer, sizeof(readBuffer));
        rapidjson::Document d;
        d.ParseStream<0, rapidjson::UTF8<>, rapidjson::FileReadStream>(is);
        fclose(fp);

        bool consecutive_reordering = true; //TODO: pass this as a flag in the options
        if(consecutive_reordering)
            EnforceConsecutiveOrdering(d);



        //construct the model part
        FillModelPart(d, rThisModelPart);

    }




    /// Constructor with filenames.

    /// Destructor.
    virtual ~KratosJsonIO() {};


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{



    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "KratosJsonIO";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "KratosJsonIO";
    }


    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

protected:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    std::string mFilename;
    Flags mOptions;


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{
    void FillModelPart(rapidjson::Document& d, ModelPart & rThisModelPart)
    {
        KRATOS_TRY

        //get the section of the model part related to nodes
        rapidjson::Value& json_model_part = d["model_part"]; // Using a reference for consecutive access is handy and faster.

        //////////////////////////////////////////////read Nodes
        if(json_model_part.HasMember("Nodes") == false) KRATOS_THROW_ERROR(std::invalid_argument, "missing the subsection Nodes of the model_part","")
        {
            rapidjson::Value& json_nodes = json_model_part["Nodes"]; // Using a reference for consecutive access is handy and faster.
            if(json_nodes.IsArray() == false) KRATOS_THROW_ERROR(std::invalid_argument, "the Nodes section is not an array as expected!!","")

                std::map<std::size_t, std::size_t> nodes_id_map; //create a list of nodes -- Create New node uses insert hence it is slow if input is not ordered
            for (rapidjson::SizeType i = 0; i < json_nodes.Size(); i++) // rapidjson uses SizeType instead of size_t.
            {
                std::size_t node_id =  json_nodes[i][0].GetInt() ;
                nodes_id_map[ node_id ] = i; //store in the map the position in the array
            }

            for(std::map<std::size_t, std::size_t>::iterator it = nodes_id_map.begin(); it!= nodes_id_map.end(); it++)
            {
                std::size_t pos = it->second;
//                 KRATOS_WATCH(it->first);
//                 KRATOS_WATCH(it->second);
                const rapidjson::Value& node_input = json_nodes[pos];             //obtain the row with the node data
                rThisModelPart.CreateNewNode( node_input[0].GetInt(), node_input[1].GetDouble(), node_input[2].GetDouble(), node_input[3].GetDouble());
            }
            //here nodes are not any more needed in the json so remove them to free the memory
            d.RemoveMember(json_nodes);
        }

        //////////////////////////////////////////////read Properties
        if(json_model_part.HasMember("Properties") == false) KRATOS_THROW_ERROR(std::invalid_argument, "missing the subsection Properties of the model_part","")
        {
            const rapidjson::Value& json_properties = json_model_part["Properties"]; // Using a reference for consecutive access is handy and faster.
            if(json_properties.IsObject() == false) KRATOS_THROW_ERROR(std::invalid_argument, "the Properties section is not a json object as expected!!","")
                for (rapidjson::Value::ConstMemberIterator itr = json_properties.MemberBegin(); itr != json_properties.MemberEnd(); ++itr)
                {
                    unsigned int prop_id = atoi( itr->name.GetString() );
                    if(itr->value.HasMember("Variables") == false) KRATOS_THROW_ERROR(std::invalid_argument, "missing the subsection Variables of the Properties subsection of model_part","");

                    const rapidjson::Value& prop_variables = itr->value["Variables"];
                    for (rapidjson::Value::ConstMemberIterator iii = prop_variables.MemberBegin(); iii != prop_variables.MemberEnd(); ++iii)
                    {
                        if( CheckAndAssignPropertyValue< Variable<double> >( rThisModelPart.pGetProperties(prop_id), iii->name.GetString(), iii->value) ) {}
                        else if( CheckAndAssignPropertyValue< Variable<bool> >( rThisModelPart.pGetProperties(prop_id), iii->name.GetString(), iii->value) ) {}
                        else if( CheckAndAssignPropertyValue< Variable<int> >( rThisModelPart.pGetProperties(prop_id), iii->name.GetString(), iii->value) ) {}
                        else if( CheckAndAssignPropertyValue< Variable<unsigned int> >( rThisModelPart.pGetProperties(prop_id), iii->name.GetString(), iii->value) ) {}
                        else
                        {
                            typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > component_type;
                            if( CheckAndAssignPropertyValue< component_type >( rThisModelPart.pGetProperties(prop_id), iii->name.GetString(), iii->value) ) {}
                            else
                            {
                                KRATOS_THROW_ERROR(std::invalid_argument, "can not read a variable type from the property block -- variable is ",iii->name.GetString());
                            }
                        }
                    }
                }
        }

        //read Elements
        if(json_model_part.HasMember("Elements"))
        {
            const rapidjson::Value& json_elements = json_model_part["Elements"]; // Using a reference for consecutive access is handy and faster.
            if(json_elements.IsObject() == false) KRATOS_THROW_ERROR(std::invalid_argument, "the Elements section is not a json object as expected!!","")
                for (rapidjson::Value::ConstMemberIterator itr = json_elements.MemberBegin(); itr != json_elements.MemberEnd(); ++itr)
                {
                    std::string element_name = itr->name.GetString();
                    KRATOS_WATCH(element_name);
                    if(itr->value["connectivity"].IsArray() == false) KRATOS_THROW_ERROR(std::invalid_argument, "the connectivity section of the Elements is not an array as expected!!","")

                        const rapidjson::Value& connectivity = itr->value["connectivity"];
                    for (SizeType i = 0; i < connectivity.Size(); i++)
                    {
                        const rapidjson::Value& local_data = connectivity[i];
                        unsigned int elem_id = local_data[0].GetInt();
                        unsigned int prop_id = local_data[1].GetInt();

                        unsigned int nnodes = local_data.Size()-2;
                        std::vector< ModelPart::IndexType > node_ids(nnodes);
                        for(unsigned int j=0; j<nnodes; j++) node_ids[j]= local_data[j+2].GetInt() ;

                        rThisModelPart.CreateNewElement(element_name, elem_id, node_ids, rThisModelPart.pGetProperties(prop_id) );
                    }
                }
        }



        //read Conditions
        if(json_model_part.HasMember("Conditions"))
        {
            const rapidjson::Value& json_conditions = json_model_part["Conditions"]; // Using a reference for consecutive access is handy and faster.
            if(json_conditions.IsObject() == false) KRATOS_THROW_ERROR(std::invalid_argument, "the Conditions section is not a json object as expected!!","")
                for (rapidjson::Value::ConstMemberIterator itr = json_conditions.MemberBegin(); itr != json_conditions.MemberEnd(); ++itr)
                {
                    std::string condition_name = itr->name.GetString();
                    if(itr->value["connectivity"].IsArray() == false) KRATOS_THROW_ERROR(std::invalid_argument, "the connectivity section of the Conditions is not an array as expected!!","")

                        const rapidjson::Value& connectivity = itr->value["connectivity"];
                    for (SizeType i = 0; i < connectivity.Size(); i++)
                    {
                        const rapidjson::Value& local_data = connectivity[i];
                        unsigned int elem_id = local_data[0].GetInt();
                        unsigned int prop_id = local_data[1].GetInt();

                        unsigned int nnodes = local_data.Size()-2;
                        std::vector< ModelPart::IndexType > node_ids(nnodes);
                        for(unsigned int j=0; j<nnodes; j++) node_ids[j]= local_data[j+2].GetInt() ;

                        rThisModelPart.CreateNewCondition(condition_name, elem_id, node_ids, rThisModelPart.pGetProperties(prop_id) );
                    }
                }
        }

        //read Meshes
        if(json_model_part.HasMember("Meshes"))
        {
            const rapidjson::Value& json_meshes = json_model_part["Meshes"]; // Using a reference for consecutive access is handy and faster.
            if(json_meshes.IsObject() == false) KRATOS_THROW_ERROR(std::invalid_argument, "the Meshes section is not a json object as expected!!","")
                for (rapidjson::Value::ConstMemberIterator itr = json_meshes.MemberBegin(); itr != json_meshes.MemberEnd(); ++itr)
                {
                    //get the mesh
                    unsigned int mesh_id = atoi( itr->name.GetString() );

                    //TODO: this is UGLY - there should be a better way to create a mesh
                    unsigned int number_of_meshes = rThisModelPart.GetMeshes().size();
                    ModelPart::MeshType empty_mesh;
                    for(ModelPart::IndexType i = number_of_meshes ; i < mesh_id + 1 ; i++)
                        rThisModelPart.GetMeshes().push_back(empty_mesh.Clone());

                    ModelPart::MeshType& mesh = rThisModelPart.GetMesh(mesh_id);

                    if((itr->value).HasMember("NodePointers"))
                    {
                        const rapidjson::Value& json_nodepointers = (itr->value)["NodePointers"];
                        if(json_nodepointers.IsArray() == false) KRATOS_THROW_ERROR(std::invalid_argument, "the NodePointers section is not a json array as expected!!","")
                            mesh.Nodes().reserve( json_nodepointers.Size() );
                        KRATOS_WATCH( json_nodepointers.Size() )
                        for(rapidjson::SizeType i=0; i<json_nodepointers.Size(); i++) mesh.Nodes().push_back( rThisModelPart.pGetNode(  json_nodepointers[i].GetInt()  ));
                    }

                    if((itr->value).HasMember("ElementPointers"))
                    {
                        const rapidjson::Value& json_elementpointers = (itr->value)["ElementPointers"];
                        if(json_elementpointers.IsArray() == false) KRATOS_THROW_ERROR(std::invalid_argument, "the ElementPointers section is not a json array as expected!!","")
                            mesh.Elements().reserve( json_elementpointers.Size() );
                        for(rapidjson::SizeType i=0; i<json_elementpointers.Size(); i++) mesh.Elements().push_back( rThisModelPart.pGetElement( json_elementpointers[i].GetInt() ));
                    }

                    if((itr->value).HasMember("ConditionPointers"))
                    {
                        const rapidjson::Value& json_conditionpointers = (itr->value)["ConditionPointers"];
                        if(json_conditionpointers.IsArray() == false) KRATOS_THROW_ERROR(std::invalid_argument, "the ConditionPointers section is not a json array as expected!!","")
                            mesh.Conditions().reserve( json_conditionpointers.Size() );
                        for(rapidjson::SizeType i=0; i<json_conditionpointers.Size(); i++) mesh.Conditions().push_back( rThisModelPart.pGetCondition( json_conditionpointers[i].GetInt() ));
                    }

                    mesh.Nodes().Unique();
                    mesh.Elements().Unique();
                    mesh.Conditions().Unique();
                }


        }


        //read NodalValues
        if(json_model_part.HasMember("NodalData"))
        {
            {
                const rapidjson::Value& json_nodaldata = json_model_part["NodalData"];
                if(json_nodaldata.IsObject() == false) KRATOS_THROW_ERROR(std::invalid_argument, "the NodalData section is not a json object as expected!!","")
                    for (rapidjson::Value::ConstMemberIterator itr = json_nodaldata.MemberBegin(); itr != json_nodaldata.MemberEnd(); ++itr)
                    {
                        std::string variable_name = itr->name.GetString();

                        if(itr->value.HasMember("stepindex") == false) KRATOS_THROW_ERROR(std::invalid_argument, "missing the subsection stepindex of the NodalData subsection of model_part","");
                        unsigned int stepindex =  itr->value["stepindex"].GetInt();

                        if(itr->value.HasMember("data") == false) KRATOS_THROW_ERROR(std::invalid_argument, "missing the subsection data of the NodalData subsection of model_part","");
                        const rapidjson::Value& json_nodaldata = itr->value["data"];

//                     for (rapidjson::Value::ConstMemberIterator iii = json_nodaldata.MemberBegin(); iii != json_nodaldata.MemberEnd(); ++iii)
//                     {
//                     for(SizeType i = 0; i<json_nodaldata.Size(); i++)
                        {
                            if( CheckAndAssignNodalValue< Variable<double> >( rThisModelPart, variable_name, json_nodaldata, stepindex) ) {}
                            else if( CheckAndAssignNodalValue_ErrorIfFixed< Variable<bool> >( rThisModelPart, variable_name, json_nodaldata,stepindex) ) {}
                            else if( CheckAndAssignNodalValue_ErrorIfFixed< Variable<int> >( rThisModelPart, variable_name, json_nodaldata,stepindex) ) {}
                            else if( CheckAndAssignNodalValue_ErrorIfFixed< Variable<unsigned int> >( rThisModelPart, variable_name, json_nodaldata,stepindex) ) {}
                            else
                            {
                                typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > component_type;
                                if( CheckAndAssignNodalComponentValue< component_type >( rThisModelPart, variable_name, json_nodaldata,stepindex) ) {}
                                else
                                {
                                    KRATOS_THROW_ERROR(std::invalid_argument, "can not read a variable type from the property block -- variable is ",variable_name);
                                }
                            }
                        }
                    }
            }
        }


        //read ElementalValues

        //read ElementalConditions

        KRATOS_CATCH("")
    }




    //this function changes the json data structure -- completely independent from the kratos data structure
    void EnforceConsecutiveOrdering(rapidjson::Document& d)
    {

        if(d["model_part"].HasMember("Nodes") == false) KRATOS_THROW_ERROR(std::invalid_argument, "missing the subsection Nodes of the model_part","");

            rapidjson::Value& json_nodes = d["model_part"]["Nodes"]; // Using a reference for consecutive access is handy and faster.
        if(json_nodes.IsArray() == false) KRATOS_THROW_ERROR(std::invalid_argument, "the Nodes section is not an array as expected!!","")


        ///////////////////////////////////////////////////
        //first ensure that nodes are ordered consecutively
        ///////////////////////////////////////////////////        //obtain the list of consecutive ids
        std::set<std::size_t> original_ids_set; //create a list of nodes -- Create New node uses insert hence it is slow if input is not ordered
        for (rapidjson::SizeType i = 0; i < json_nodes.Size(); i++) // rapidjson uses SizeType instead of size_t.
        {
            std::size_t node_id = json_nodes[i][0].GetInt();
            original_ids_set.insert( node_id ); //store in the map the position in the array
        }

        //generate a map assigning to every original id a new consecutive_id
        std::map < std::size_t, std::size_t > consecutive_nodeids_map;
        std::size_t consecutive_id = 1;
        for(std::set<std::size_t>::iterator it = original_ids_set.begin(); it!= original_ids_set.end(); it++)
        {
            consecutive_nodeids_map[ *it ] = consecutive_id;
            consecutive_id++;
        }

        ///////////////////////////////////////////////////
        //now order consecutively elements - reuse original_ids_set since it is not needed anymore - note that it requires two loops!
        ///////////////////////////////////////////////////
        original_ids_set.clear();
        if(d["model_part"].HasMember("Elements"))
        {
            rapidjson::Value& json_elements = d["model_part"]["Elements"]; // Using a reference for consecutive access is handy and faster.
            if(json_elements.IsObject() == false) KRATOS_THROW_ERROR(std::invalid_argument, "the Elements section is not a json object as expected!!","");
                for (rapidjson::Value::MemberIterator itr = json_elements.MemberBegin(); itr != json_elements.MemberEnd(); ++itr)
                {

                    rapidjson::Value& connectivity = itr->value["connectivity"];

                    for (SizeType i = 0; i < connectivity.Size(); i++)
                    {
                        rapidjson::Value& local_data = connectivity[i];
                        std::size_t original_id = local_data[0].GetInt();
                        original_ids_set.insert(original_id);
                    }
                }
        }
        std::map < std::size_t, std::size_t > consecutive_elementids_map;
        consecutive_id = 1;
        for(std::set<std::size_t>::iterator it = original_ids_set.begin(); it!= original_ids_set.end(); it++)
        {
            consecutive_elementids_map[ *it ] = consecutive_id;
            consecutive_id++;
        }

        ///////////////////////////////////////////////////
        //finally order consecutively conditions - reuse original_ids_set since it is not needed anymore - note that it requires two loops!
        ///////////////////////////////////////////////////
        original_ids_set.clear();
        if(d["model_part"].HasMember("Conditions"))
        {
            rapidjson::Value& json_conditions = d["model_part"]["Conditions"]; // Using a reference for consecutive access is handy and faster.
            if(json_conditions.IsObject() == false) KRATOS_THROW_ERROR(std::invalid_argument, "the Conditions section is not a json object as expected!!","");
                for (rapidjson::Value::MemberIterator itr = json_conditions.MemberBegin(); itr != json_conditions.MemberEnd(); ++itr)
                {

                    rapidjson::Value& connectivity = itr->value["connectivity"];

                    for (SizeType i = 0; i < connectivity.Size(); i++)
                    {
                        rapidjson::Value& local_data = connectivity[i];
                        std::size_t original_id = local_data[0].GetInt();
                        original_ids_set.insert(original_id);
                    }
                }
        }
        std::map < std::size_t, std::size_t > consecutive_conditionids_map;
        consecutive_id = 1;
        for(std::set<std::size_t>::iterator it = original_ids_set.begin(); it!= original_ids_set.end(); it++)
        {
            consecutive_conditionids_map[ *it ] = consecutive_id;
            consecutive_id++;
        }


        //now modify the json data inplace for the new consecutive_nodeids_map
        for(rapidjson::SizeType i = 0; i < json_nodes.Size(); i++) // rapidjson uses SizeType instead of size_t.
        {
            json_nodes[i][0].SetInt( consecutive_nodeids_map[ json_nodes[i][0].GetInt() ] );
        }

        if(d["model_part"].HasMember("NodalData"))
        {
            {
                rapidjson::Value& json_nodaldata = d["model_part"]["NodalData"];
                if(json_nodaldata.IsObject() == false) KRATOS_THROW_ERROR(std::invalid_argument, "the NodalData section is not a json object as expected!!","")
                    for (rapidjson::Value::MemberIterator itr = json_nodaldata.MemberBegin(); itr != json_nodaldata.MemberEnd(); ++itr)
                    {
                        if(itr->value.HasMember("data") == false) KRATOS_THROW_ERROR(std::invalid_argument, "missing the subsection data of the NodalData subsection of model_part","");
                        rapidjson::Value& json_nodaldata = itr->value["data"];

                        for(rapidjson::SizeType i = 0; i<json_nodaldata.Size(); i++)
                        {
                            rapidjson::Value& value = json_nodaldata[i];

                            value[0].SetInt( consecutive_nodeids_map[ value[0].GetInt() ] );
                        }
                    }
            }
        }


        //loop over node pointers
        if(d["model_part"].HasMember("Meshes"))
        {
            rapidjson::Value& json_meshes = d["model_part"]["Meshes"]; // Using a reference for consecutive access is handy and faster.
            if(json_meshes.IsObject() == false) KRATOS_THROW_ERROR(std::invalid_argument, "the Meshes section is not a json object as expected!!","")
                for (rapidjson::Value::MemberIterator itr = json_meshes.MemberBegin(); itr != json_meshes.MemberEnd(); ++itr)
                {
                    //get the mesh
                    if((itr->value).HasMember("NodePointers"))
                    {
                        rapidjson::Value& json_nodepointers = (itr->value)["NodePointers"];
                        if(json_nodepointers.IsArray() == false) KRATOS_THROW_ERROR(std::invalid_argument, "the NodePointers section is not a json array as expected!!","")
                            for(rapidjson::SizeType i=0; i<json_nodepointers.Size(); i++) json_nodepointers[i].SetInt(  consecutive_nodeids_map[ json_nodepointers[i].GetInt()]  );
                    }

                    if((itr->value).HasMember("ElementPointers"))
                    {
                        rapidjson::Value& json_elementpointers = (itr->value)["ElementPointers"];
                        if(json_elementpointers.IsArray() == false) KRATOS_THROW_ERROR(std::invalid_argument, "the ElementPointers section is not a json array as expected!!","")
                        for(rapidjson::SizeType i=0; i<json_elementpointers.Size(); i++) json_elementpointers[i].SetInt( consecutive_elementids_map[ json_elementpointers[i].GetInt()] );
                    }

                    if((itr->value).HasMember("ConditionPointers"))
                    {
                        rapidjson::Value& json_conditionpointers = (itr->value)["ConditionPointers"];
                        if(json_conditionpointers.IsArray() == false) KRATOS_THROW_ERROR(std::invalid_argument, "the ConditionPointers section is not a json array as expected!!","")
                        for(rapidjson::SizeType i=0; i<json_conditionpointers.Size(); i++) json_conditionpointers[i].SetInt( consecutive_elementids_map[ json_conditionpointers[i].GetInt()] );
                    }
                }
        }

        //loop over nodal data


        //now loop over all Nodes Elements and Conditions and change the ids to the one contained in the values of nodes_id_map
        if(d["model_part"].HasMember("Elements"))
        {
            rapidjson::Value& json_elements = d["model_part"]["Elements"]; // Using a reference for consecutive access is handy and faster.
            if(json_elements.IsObject() == false) KRATOS_THROW_ERROR(std::invalid_argument, "the Elements section is not a json object as expected!!","")
                for (rapidjson::Value::MemberIterator itr = json_elements.MemberBegin(); itr != json_elements.MemberEnd(); ++itr)
                {

                    rapidjson::Value& connectivity = itr->value["connectivity"];

                    for (SizeType i = 0; i < connectivity.Size(); i++)
                    {
                        rapidjson::Value& local_data = connectivity[i];
                        local_data[0].SetInt( consecutive_elementids_map[ local_data[0].GetInt() ]) ; //change the id of the element
                        unsigned int nnodes = local_data.Size()-2;
                        for(unsigned int j=0; j<nnodes; j++) local_data[j+2].SetInt( consecutive_nodeids_map[ local_data[j+2].GetInt() ]) ;

                    }
                }
        }

        //now loop over all Nodes Elements and Conditions and change the ids to the one contained in the values of nodes_id_map
        if(d["model_part"].HasMember("Conditions"))
        {
            rapidjson::Value& json_elements = d["model_part"]["Conditions"]; // Using a reference for consecutive access is handy and faster.
            if(json_elements.IsObject() == false) KRATOS_THROW_ERROR(std::invalid_argument, "the Conditions section is not a json object as expected!!","")
                for (rapidjson::Value::MemberIterator itr = json_elements.MemberBegin(); itr != json_elements.MemberEnd(); ++itr)
                {
                    std::string element_name = itr->name.GetString();
                    if(itr->value["connectivity"].IsArray() == false) KRATOS_THROW_ERROR(std::invalid_argument, "the connectivity section of the Conditions is not an array as expected!!","");

                    rapidjson::Value& connectivity = itr->value["connectivity"];
                    for (SizeType i = 0; i < connectivity.Size(); i++)
                    {
                        rapidjson::Value& local_data = connectivity[i];
                        local_data[0].SetInt( consecutive_conditionids_map[ local_data[0].GetInt() ]) ; //change the id of the condition
                        unsigned int nnodes = local_data.Size()-2;
                        for(unsigned int j=0; j<nnodes; j++) local_data[j+2].SetInt( consecutive_nodeids_map[ local_data[j+2].GetInt() ]) ;

                    }
                }
        }



    }

    template< class TVarType >
    bool CheckAndAssignNodalValue( ModelPart& r_model_part,  std::string variable_name, const rapidjson::Value& json_value, unsigned int stepindex)
    {
        if( KratosComponents< TVarType >::Has( variable_name ) )
        {
            const TVarType& rVar = KratosComponents< TVarType >::Get( variable_name );
            if( r_model_part.GetNodalSolutionStepVariablesList().Has( rVar ) == false )
                KRATOS_THROW_ERROR(std::runtime_error,"trying to assign a variable that is not in the model_part - variable name is ",variable_name);

            for(rapidjson::SizeType i = 0; i<json_value.Size(); i++) //loop over all the nodes with prescribed velocity
            {
                const rapidjson::Value& value = json_value[i];

                unsigned int node_id = value[0].GetInt();
                bool is_fixed = static_cast<bool>(value[1].GetInt());
                double double_value = value[2].GetDouble();

                ModelPart::NodesContainerType::iterator itnode = r_model_part.Nodes().find( node_id );
                if(itnode == r_model_part.NodesEnd() ) KRATOS_THROW_ERROR(std::runtime_error, "node not found when assigning NodalData --> Node Id is ", node_id);

                itnode->GetSolutionStepValue(rVar, stepindex) = double_value;
                if(is_fixed) itnode->Fix(rVar);
            }

            return true;
        }
        return false; //variable was not of the given type and nothing was done
    }

    template< class TVarType >
    bool CheckAndAssignNodalComponentValue( ModelPart& r_model_part,  std::string variable_name, const rapidjson::Value& json_value, unsigned int stepindex)
    {
        if( KratosComponents< TVarType >::Has( variable_name ) )
        {
            const TVarType& rVar = KratosComponents< TVarType >::Get( variable_name );

            const std::string base_variable_name = rVar.GetSourceVariable().Name();
            const  Variable<array_1d<double,3> >& rVectorBaseVar = KratosComponents< Variable<array_1d<double,3> > >::Get( base_variable_name );

            if( r_model_part.GetNodalSolutionStepVariablesList().Has( rVectorBaseVar ) == false )
                KRATOS_THROW_ERROR(std::runtime_error,"trying to assign a variable that is not in the model_part - variable name is ",variable_name);

            for(rapidjson::SizeType i = 0; i<json_value.Size(); i++) //loop over all the nodes with prescribed velocity
            {
                const rapidjson::Value& value = json_value[i];

                unsigned int node_id = value[0].GetInt();
                bool is_fixed = static_cast<bool>(value[1].GetInt());
                double double_value = value[2].GetDouble();

                ModelPart::NodesContainerType::iterator itnode = r_model_part.Nodes().find( node_id );
                if(itnode == r_model_part.NodesEnd() ) KRATOS_THROW_ERROR(std::runtime_error, "node not found when assigning NodalData --> Node Id is ", node_id);

                itnode->GetSolutionStepValue(rVar, stepindex) = double_value;
                if(is_fixed) itnode->Fix(rVar);
            }

            return true;
        }
        return false; //variable was not of the given type and nothing was done
    }

    template< class TVarType >
    bool CheckAndAssignNodalValue_ErrorIfFixed( ModelPart& r_model_part,  std::string variable_name, const rapidjson::Value& json_value, unsigned int stepindex)
    {
        if( KratosComponents< TVarType >::Has( variable_name ) )
        {
            const TVarType& rVar = KratosComponents< TVarType >::Get( variable_name );
            if( r_model_part.GetNodalSolutionStepVariablesList().Has( rVar ) == false )
                KRATOS_THROW_ERROR(std::runtime_error,"trying to assign a variable that is not in the model_part - variable name is ",variable_name);

            for(rapidjson::SizeType i = 0; i<json_value.Size(); i++) //loop over all the nodes with prescribed velocity
            {
                const rapidjson::Value& value = json_value[i];

                unsigned int node_id = value[0].GetInt();
                bool is_fixed = static_cast<bool>(value[1].GetInt());
                double double_value = value[2].GetDouble();

                ModelPart::NodesContainerType::iterator itnode = r_model_part.Nodes().find( node_id );
                if(itnode == r_model_part.NodesEnd() ) KRATOS_THROW_ERROR(std::runtime_error, "node not found when assigning NodalData --> Node Id is ", node_id);

                itnode->GetSolutionStepValue(rVar, stepindex) = double_value;
                if(is_fixed)
                {
                    std::string err_msg = std::string("variable name is ")+variable_name; //+std::string(" node Id is ") + std::string(itoa(node_id);
                    KRATOS_THROW_ERROR(std::runtime_error, "sorry this variable is not of double or component type and can not be fixed -- ", err_msg);
                }
            }

            return true;
        }
        return false; //variable was not of the given type and nothing was done
    }


    template< class TVarType >
    bool CheckAndAssignPropertyValue( Properties::Pointer pProperty,  std::string variable_name, const rapidjson::Value& json_value)
    {
        if( KratosComponents< TVarType >::Has( variable_name ) )
        {
            const TVarType& rVar = KratosComponents< TVarType >::Get( variable_name );
            pProperty->GetValue( rVar ) = json_value.GetDouble();
            return true;
        }
        return false; //variable was not of the given type and nothing was done
    }


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    KratosJsonIO& operator=(KratosJsonIO const& rOther);

    /// Copy constructor.
    KratosJsonIO(KratosJsonIO const& rOther);


    ///@}

}; // Class KratosJsonIO

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


//   /// input stream function
//   inline std::istream& operator >> (std::istream& rIStream,
//                  KratosJsonIO& rThis);

//   /// output stream function
//   inline std::ostream& operator << (std::ostream& rOStream,
//                  const KratosJsonIO& rThis)
//     {
//       rThis.PrintInfo(rOStream);
//       rOStream << std::endl;
//       rThis.PrintData(rOStream);

//       return rOStream;
//     }
//   ///@}


}  // namespace Kratos.

#endif // KRATOS_JSON_IO_H_INCLUDED  defined
