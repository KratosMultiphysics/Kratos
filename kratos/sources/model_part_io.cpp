//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//

// System includes
#include<unordered_set>

// External includes

// Project includes
#include "includes/model_part_io.h"
#include "input_output/logger.h"
#include "utilities/quaternion.h"
#include "utilities/openmp_utils.h"
#include "utilities/compare_elements_and_conditions_utility.h"

namespace Kratos
{
/// Constructor with  filenames.
ModelPartIO::ModelPartIO(std::string const& Filename, const Flags Options)
    : mNumberOfLines(1)
    , mBaseFilename(Filename)
    , mFilename(Filename + ".mdpa")
    , mOptions(Options)
{
    Kratos::shared_ptr<std::fstream> pFile = Kratos::make_shared<std::fstream>();
    std::fstream::openmode OpenMode;

    // Set the mode
    if (mOptions.Is(IO::READ))
    {
        OpenMode = std::fstream::in;
    }
    else if (mOptions.Is(IO::APPEND))
    {
        OpenMode = std::fstream::in | std::fstream::app;
    }
    else if (mOptions.Is(IO::WRITE))
    {
        OpenMode = std::fstream::out;
    }
    else
    {
        // If none of the READ, WRITE or APPEND are defined we will take READ as
        // default.
        OpenMode = std::fstream::in;
    }

    pFile->open(mFilename.c_str(), OpenMode);

    KRATOS_ERROR_IF_NOT(pFile->is_open()) << "Error opening mdpa file : " << mFilename.c_str() << std::endl;

    // Store the pointer as a regular std::iostream
    mpStream = pFile;

    if (mOptions.IsNot(IO::SKIP_TIMER)) Timer::SetOuputFile(Filename + ".time");
}

/// Constructor with stream
ModelPartIO::ModelPartIO(Kratos::shared_ptr<std::iostream> Stream, const Flags Options)
    : mNumberOfLines(1)
    , mOptions(Options)
{
    // Check if the pointer is valid
    KRATOS_ERROR_IF(Stream == nullptr) << "Error: ModelPartIO Stream is invalid " << std::endl;

    // Check if the pointer was .reset() or never initialized and if its a NULL pointer)
    KRATOS_ERROR_IF(Stream == nullptr || Stream == Kratos::shared_ptr<std::iostream>(NULL)) << "Error: ModelPartIO Stream is invalid " << std::endl;

    mpStream = Stream;
}

/// Destructor.
ModelPartIO::~ModelPartIO() {
    if (mOptions.IsNot(IO::SKIP_TIMER)) Timer::CloseOuputFile();
}

///@}
///@name Operators
///@{


///@}
///@name Operations
///@{

bool ModelPartIO::ReadNode(NodeType& rThisNode)
{
    KRATOS_ERROR << "Calling base class member. Please check the definition of derived class." << std::endl;
}

bool ModelPartIO::ReadNodes(NodesContainerType& rThisNodes)
{
    KRATOS_TRY
    ResetInput();
    std::string word;
    while(true)
    {
        ReadWord(word);
        if(mpStream->eof())
            break;
        ReadBlockName(word);
        if(word == "Nodes")
            ReadNodesBlock(rThisNodes);
        else
            SkipBlock(word);
    }

    return true;
    KRATOS_CATCH("")
}

std::size_t ModelPartIO::ReadNodesNumber()
{
    KRATOS_TRY;
    ResetInput();
    std::string word;
    std::size_t num_nodes = 0;
    while(true)
    {
        ReadWord(word);
        if (mpStream->eof())
            break;
        ReadBlockName(word);
        if(word == "Nodes")
            num_nodes += CountNodesInBlock();
        else
            SkipBlock(word);
    }

    return num_nodes;
    KRATOS_CATCH("");
}

void ModelPartIO::WriteNodes(NodesContainerType const& rThisNodes)
{
    (*mpStream) << "Begin Nodes" << std::endl;
    for(NodesContainerType::const_iterator it_node = rThisNodes.begin() ; it_node != rThisNodes.end() ; it_node++)
        (*mpStream) << "\t" << it_node->Id() << "\t" << it_node->X()  << "\t" << it_node->Y() << "\t" << it_node->Z() << std::endl;
    (*mpStream) << "End Nodes" << std::endl << std::endl;
}

void ModelPartIO::ReadProperties(Properties& rThisProperties)
{
    KRATOS_ERROR << "Calling base class member. Please check the definition of derived class" << std::endl;
}

void ModelPartIO::ReadProperties(PropertiesContainerType& rThisProperties)
{
    KRATOS_TRY
    ResetInput();
    std::string word;
    while(true)
    {
        ReadWord(word);
        if(mpStream->eof())
            break;
        ReadBlockName(word);
        if(word == "Properties")
            ReadPropertiesBlock(rThisProperties);
        else
            SkipBlock(word);
    }
    KRATOS_CATCH("")
}

void ModelPartIO::WriteProperties(PropertiesContainerType& rThisProperties)
{
    std::string aux_string;
    const std::string string_to_remove = "This properties contains 0 tables";

    // We write at least one empty property
    if (rThisProperties.size() == 0) {
        (*mpStream) << "Begin Properties 0" << std::endl;
        (*mpStream) << "End Properties" << std::endl << std::endl;
    }

    // General case where the at leat one property is defined
    for (auto i_properties = rThisProperties.begin() ; i_properties != rThisProperties.end() ; ++i_properties) {
        std::ostringstream aux_ostream;
        (*mpStream) << "Begin Properties " << i_properties->Id() << std::endl;
        i_properties->PrintData(aux_ostream);

        aux_string = aux_ostream.str();

        // We remove the line of Constitutive Laws and we add it manually. We do this because when calling the Info() method to the data_value_container it returns the address of the pointer, and we are interested in the name
        if (i_properties->Has(CONSTITUTIVE_LAW)) {
            // First we remove the entire line (we will add it by ourselves later)
            std::string::size_type it_constitutive_law_begin = aux_string.find("CONSTITUTIVE_LAW");

            if (it_constitutive_law_begin != std::string::npos) {
            std::string::size_type it_constitutive_law_end = aux_string.find("\n", it_constitutive_law_begin);
                aux_string.erase(it_constitutive_law_begin, it_constitutive_law_end);
            }

            // Now we look for the constitutive law name, we do something similar to elements and conditions. We iterate over the database of the Kratos Components until the type is the same. The IsSameType should be implemented properly in the CL class
            const ConstitutiveLaw::Pointer p_law = i_properties->GetValue(CONSTITUTIVE_LAW);
            auto components_cl = KratosComponents<ConstitutiveLaw>::GetComponents();
            std::string cl_name = "";
            for (const auto& comp_cl : components_cl) {
                if (p_law->HasSameType(p_law.get(), comp_cl.second)) {
                    cl_name = comp_cl.first;
                    break;
                }
            }
            if (cl_name != "") aux_string += "CONSTITUTIVE_LAW " + cl_name + "\n";
        }


        std::string::size_type it_to_remove = aux_string.find(string_to_remove);

        if (it_to_remove != std::string::npos) {
            aux_string.erase(it_to_remove, string_to_remove.length());
        }

        aux_string.erase(std::remove(aux_string.begin(), aux_string.end(), ':'), aux_string.end());

        (*mpStream) << aux_string << std::endl;
        (*mpStream) << "End Properties" << std::endl << std::endl;
    }
}

void ModelPartIO::ReadElement(NodesContainerType& rThisNodes, PropertiesContainerType& rThisProperties, Element::Pointer& pThisElements)
{
    KRATOS_ERROR << "Calling base class member. Please check the definition of derived class" << std::endl;
}

void ModelPartIO::ReadElements(NodesContainerType& rThisNodes, PropertiesContainerType& rThisProperties, ElementsContainerType& rThisElements)
{
    KRATOS_TRY
    ResetInput();
    std::string word;
    while(true)
    {
        ReadWord(word);
        if(mpStream->eof())
            break;
        ReadBlockName(word);
        if(word == "Elements")
            ReadElementsBlock(rThisNodes,rThisProperties,rThisElements);
        else
            SkipBlock(word);
    }
    KRATOS_CATCH("")
}

std::size_t  ModelPartIO::ReadElementsConnectivities(ConnectivitiesContainerType& rElementsConnectivities)
{
    KRATOS_TRY
    std::size_t number_of_elements = 0;
    ResetInput();
    std::string word;
    while(true)
    {
        ReadWord(word);
        if(mpStream->eof())
            break;
        ReadBlockName(word);
        if(word == "Elements")
            number_of_elements += ReadElementsConnectivitiesBlock(rElementsConnectivities);
        else
            SkipBlock(word);
    }
    return number_of_elements;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void ModelPartIO::WriteElements(ElementsContainerType const& rThisElements)
{
    // We are going to procede like the following, we are going to iterate over all the elements and compare with the components, we will save the type and we will compare until we get that the type of element has changed
    if (rThisElements.size() > 0) {
        std::string element_name;

        auto it_element = rThisElements.begin();
        auto elements_components = KratosComponents<Element>::GetComponents();

        // Fisrt we do the first element
        CompareElementsAndConditionsUtility::GetRegisteredName(*it_element, element_name);

        (*mpStream) << "Begin Elements\t" << element_name << std::endl;
        (*mpStream) << "\t" << rThisElements.begin()->Id() << "\t" << (rThisElements.begin()->pGetProperties())->Id() << "\t";
        for (std::size_t i_node = 0; i_node < rThisElements.begin()->GetGeometry().size(); i_node++)
            (*mpStream) << rThisElements.begin()->GetGeometry()[i_node].Id() << "\t";
        (*mpStream) << std::endl;

        // Now we iterate over all the elements
        for(std::size_t i = 1; i < rThisElements.size(); i++) {
            auto it_elem_previous = rThisElements.begin() + i - 1;
            auto it_elem_current = rThisElements.begin() + i;

            if(GeometricalObject::IsSame(*it_elem_previous, *it_elem_current)) {
                (*mpStream) << "\t" << it_elem_current->Id() << "\t" << (it_elem_current->pGetProperties())->Id() << "\t";
                for (std::size_t i_node = 0; i_node < it_elem_current->GetGeometry().size(); i_node++)
                    (*mpStream) << it_elem_current->GetGeometry()[i_node].Id() << "\t";
                (*mpStream) << std::endl;
            } else {
                (*mpStream) << "End Elements" << std::endl << std::endl;;

                CompareElementsAndConditionsUtility::GetRegisteredName(*it_elem_current, element_name);

                (*mpStream) << "Begin Elements\t" << element_name << std::endl;
                (*mpStream) << "\t" << it_elem_current->Id() << "\t" << (it_elem_current->pGetProperties())->Id() << "\t";
                for (std::size_t i_node = 0; i_node < it_elem_current->GetGeometry().size(); i_node++)
                    (*mpStream) << it_elem_current->GetGeometry()[i_node].Id() << "\t";
                (*mpStream) << std::endl;
            }
        }

        (*mpStream) << "End Elements" << std::endl << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void ModelPartIO::ReadConditions(NodesContainerType& rThisNodes, PropertiesContainerType& rThisProperties, ConditionsContainerType& rThisConditions)
{
    KRATOS_TRY
    ResetInput();
    std::string word;
    while(true)
    {
        ReadWord(word);
        if(mpStream->eof())
            break;
        ReadBlockName(word);
        if(word == "Conditions")
            ReadConditionsBlock(rThisNodes,rThisProperties,rThisConditions);
        else
            SkipBlock(word);
    }
    KRATOS_CATCH("")
}

std::size_t  ModelPartIO::ReadConditionsConnectivities(ConnectivitiesContainerType& rConditionsConnectivities)
{
    KRATOS_TRY
    std::size_t number_of_elements = 0;
    ResetInput();
    std::string word;
    while(true)
    {
        ReadWord(word);
        if(mpStream->eof())
            break;
        ReadBlockName(word);
        if(word == "Conditions")
            number_of_elements += ReadConditionsConnectivitiesBlock(rConditionsConnectivities);
        else
            SkipBlock(word);
    }
    return number_of_elements;
    KRATOS_CATCH("")
}

void ModelPartIO::WriteConditions(ConditionsContainerType const& rThisConditions)
{
    // We are going to procede like the following, we are going to iterate over all the conditions and compare with the components, we will save the type and we will compare until we get that the type of condition has changed

    if (rThisConditions.size() > 0) {
        std::string condition_name;

        auto it_condition = rThisConditions.begin();
        auto conditions_components = KratosComponents<Condition>::GetComponents();

        // Fisrt we do the first condition
        CompareElementsAndConditionsUtility::GetRegisteredName(*it_condition, condition_name);

        (*mpStream) << "Begin Conditions\t" << condition_name << std::endl;
        (*mpStream) << "\t" << rThisConditions.begin()->Id() << "\t" << (rThisConditions.begin()->pGetProperties())->Id() << "\t";
        for (std::size_t i_node = 0; i_node < rThisConditions.begin()->GetGeometry().size(); i_node++)
            (*mpStream) << rThisConditions.begin()->GetGeometry()[i_node].Id() << "\t";
        (*mpStream) << std::endl;

        // Now we iterate over all the conditions
        for(std::size_t i = 1; i < rThisConditions.size(); i++) {
            auto it_cond_previous = rThisConditions.begin() + i - 1;
            auto it_cond_current = rThisConditions.begin() + i;

            if(GeometricalObject::IsSame(*it_cond_previous, *it_cond_current)) {
                (*mpStream) << "\t" << it_cond_current->Id() << "\t" << (it_cond_current->pGetProperties())->Id() << "\t";
                for (std::size_t i_node = 0; i_node < it_cond_current->GetGeometry().size(); i_node++)
                    (*mpStream) << it_cond_current->GetGeometry()[i_node].Id() << "\t";
                (*mpStream) << std::endl;
            } else {
                (*mpStream) << "End Conditions" << std::endl << std::endl;;

                CompareElementsAndConditionsUtility::GetRegisteredName(*it_cond_current, condition_name);

                (*mpStream) << "Begin Conditions\t" << condition_name << std::endl;
                (*mpStream) << "\t" << it_cond_current->Id() << "\t" << (it_cond_current->pGetProperties())->Id() << "\t";
                for (std::size_t i_node = 0; i_node < it_cond_current->GetGeometry().size(); i_node++)
                    (*mpStream) << it_cond_current->GetGeometry()[i_node].Id() << "\t";
                (*mpStream) << std::endl;
            }
        }

        (*mpStream) << "End Conditions" << std::endl << std::endl;
    }
}

void ModelPartIO::ReadInitialValues(ModelPart& rThisModelPart)
{
    KRATOS_TRY

    ElementsContainerType& rThisElements = rThisModelPart.Elements();
    ConditionsContainerType& rThisConditions = rThisModelPart.Conditions();


    ResetInput();
    std::string word;
    while(true)
    {
        ReadWord(word);
        if(mpStream->eof())
            break;
        ReadBlockName(word);
        if(word == "NodalData")
            ReadNodalDataBlock(rThisModelPart);
        else if(word == "ElementalData")
            ReadElementalDataBlock(rThisElements);
        else if(word == "ConditionalData")
            ReadConditionalDataBlock(rThisConditions);
        else
            SkipBlock(word);
    }
    KRATOS_CATCH("")
}

//       void ReadGeometries(NodesContainerType& rThisNodes, GeometriesContainerType& rResults);

void ModelPartIO::ReadMesh(MeshType & rThisMesh)
{
    KRATOS_ERROR << "ModelPartIO does not implement this method." << std::endl;
}

void ModelPartIO::WriteMesh(MeshType & rThisMesh)
{
    WriteProperties(rThisMesh.Properties());
    WriteNodes(rThisMesh.Nodes());
    WriteElements(rThisMesh.Elements());
    WriteConditions(rThisMesh.Conditions());
}

void ModelPartIO::ReadModelPart(ModelPart & rThisModelPart)
{
    KRATOS_TRY

    if (mOptions.IsNot(IO::SKIP_TIMER)) Timer::Start("Reading Input");

    ResetInput();
    std::string word;
    while(true)
    {
        ReadWord(word);
        if(mpStream->eof())
            break;
        ReadBlockName(word);
        if(word == "ModelPartData") {
            if (mOptions.IsNot(IO::MESH_ONLY))
                ReadModelPartDataBlock(rThisModelPart);
        } else if(word == "Table") {
            if (mOptions.IsNot(IO::MESH_ONLY))
                ReadTableBlock(rThisModelPart.Tables());
        } else if(word == "Properties") {
            ReadPropertiesBlock(rThisModelPart.rProperties());
        } else if(word == "Nodes") {
            ReadNodesBlock(rThisModelPart);
        } else if(word == "Elements") {
            ReadElementsBlock(rThisModelPart);
        } else if(word == "Conditions") {
            ReadConditionsBlock(rThisModelPart);
        } else if(word == "NodalData") {
            if (mOptions.IsNot(IO::MESH_ONLY))
                ReadNodalDataBlock(rThisModelPart);
        } else if(word == "ElementalData") {
            if (mOptions.IsNot(IO::MESH_ONLY))
                ReadElementalDataBlock(rThisModelPart.Elements());
        } else if (word == "ConditionalData") {
            if (mOptions.IsNot(IO::MESH_ONLY))
                ReadConditionalDataBlock(rThisModelPart.Conditions());
        } else if(word == "CommunicatorData") {
            if (mOptions.IsNot(IO::MESH_ONLY)) {
                ReadCommunicatorDataBlock(rThisModelPart.GetCommunicator(), rThisModelPart.Nodes());
                //Adding the elements and conditions to the communicator
                rThisModelPart.GetCommunicator().LocalMesh().Elements() = rThisModelPart.Elements();
                rThisModelPart.GetCommunicator().LocalMesh().Conditions() = rThisModelPart.Conditions();
            }
        } else if (word == "Mesh") {
            ReadMeshBlock(rThisModelPart);
        } else if (word == "SubModelPart") {
            ReadSubModelPartBlock(rThisModelPart, rThisModelPart);
        }
    }
    KRATOS_INFO("ModelPartIO") << "  [Total Lines Read : " << mNumberOfLines<<"]" << std::endl;
    if (mOptions.IsNot(IO::SKIP_TIMER)) Timer::Stop("Reading Input");
    KRATOS_CATCH("")
}

void ModelPartIO::WriteModelPart(ModelPart & rThisModelPart)
{
    if (mOptions.IsNot(IO::SKIP_TIMER)) Timer::Start("Writing Output");

    // Setting the buffer size
//     size_t size_buffer = 4096; // Look to modify this
//     char Buffer[size_buffer];
//     mpStream->rdbuf()->pubsetbuf(Buffer, size_buffer);
//
//     WriteModelPartDataBlock(rThisModelPart); // TODO: FINISH ME

    if (mOptions.IsNot(IO::MESH_ONLY))
        WriteTableBlock(rThisModelPart.Tables());
    WriteMesh(rThisModelPart.GetMesh());
    if (mOptions.IsNot(IO::MESH_ONLY)) {
        WriteNodalDataBlock(rThisModelPart); // TODO: FINISH ME
        WriteDataBlock(rThisModelPart.Elements(), "Element");
        WriteDataBlock(rThisModelPart.Conditions(),"Condition");
    }
//     WriteCommunicatorDataBlock(); // TODO: FINISH ME
//     WriteMeshBlock(rThisModelPart); // TODO: FINISH ME
    WriteSubModelPartBlock(rThisModelPart, "");

    KRATOS_INFO("ModelPartIO") << "  [Total Lines Wrote : " << mNumberOfLines<<"]" << std::endl;

    if (mOptions.IsNot(IO::SKIP_TIMER)) Timer::Stop("Writing Output");
}


/// Read the input file and create the nodal connectivities graph, stored in CSR format.
/**
    * This function produces input for Metis' nodal graph partitioning algorithms.
    * The nodal graph is stored as a (compressed) matrix where index (i,j) is non-zero if
    * there is an edge in the mesh joining nodes i and j (note that nodes are numbered from zero here,
    * to make integration with Metis simpler).
    * @param NodeIndices After call, will point to C array of size NumNodes+1 containing the
    * first CSR array: entries related to node k are stored between positions (*NodeIndices)[k]
    * and (*NodeIndices)[k+1] of *NodeConnectivities.
    * @param NodeConnectivities After call, will point to a C array of size (*NodeIndices)[NumNodes].
    * entries between (*NodeIndices)[k] and (*NodeIndices)[k+1] are a list of all nodes connected
    * to node k (counting from 0).
    * @return Number of nodes.
    */
std::size_t ModelPartIO::ReadNodalGraph(ConnectivitiesContainerType& aux_connectivities)
{
    // 1. Define an auxiliary vector of vectors
    //ConnectivitiesContainerType aux_connectivities(0);

    // 2. Fill the auxiliary vector by reading elemental and conditional connectivities
    ResetInput();
    std::string word;
    while(true)
    {
        ReadWord(word);
        if(mpStream->eof())
            break;
        ReadBlockName(word);
        if (word == "Elements")
            FillNodalConnectivitiesFromElementBlock(aux_connectivities);
        else if (word == "Conditions")
            FillNodalConnectivitiesFromConditionBlock(aux_connectivities);
        else
            SkipBlock(word);
    }

    // Check that node indices are consecutive
    unsigned int n = 0;
    for (ConnectivitiesContainerType::iterator inode = aux_connectivities.begin(); inode != aux_connectivities.end(); inode++)
    {
        n++;
        if (inode->size() == 0)
        {
            std::stringstream msg;
            msg << "Nodes are not consecutively numbered. Node " << n << " was not found in mdpa file." << std::endl;
            KRATOS_ERROR << msg.str() << std::endl;
        }
    }

    // 3. Sort each entry in the auxiliary connectivities vector, remove duplicates
    SizeType num_entries = 0;
    for (ConnectivitiesContainerType::iterator it = aux_connectivities.begin(); it != aux_connectivities.end(); it++)
    {
        std::sort(it->begin(),it->end());
        std::vector<SizeType>::iterator unique_end = std::unique(it->begin(),it->end());
        it->resize(unique_end - it->begin());
        num_entries += it->size();
    }
    SizeType num_nodes = aux_connectivities.size();

    /*// 4. Write connectivity data in CSR format
    SizeType num_nodes = aux_connectivities.size();
    *NodeIndices = new int[num_nodes+1];
    (*NodeIndices)[0] = 0;
    *NodeConnectivities = new int[num_entries];

    SizeType i = 0;
    SizeType aux_index = 0;

    for (ConnectivitiesContainerType::iterator it = aux_connectivities.begin(); it != aux_connectivities.end(); it++)
    {
        for (std::vector<SizeType>::iterator entry_it = it->begin(); entry_it != it->end(); entry_it++)
            (*NodeConnectivities)[aux_index++] = (*entry_it - 1); // substract 1 to make Ids start from 0
        (*NodeIndices)[i++] = aux_index;
    }*/

    return num_nodes;
}

void ModelPartIO::DivideInputToPartitions(SizeType NumberOfPartitions, GraphType const& DomainsColoredGraph,
                                        PartitionIndicesType const& NodesPartitions,
                                        PartitionIndicesType const& ElementsPartitions,
                                        PartitionIndicesType const& ConditionsPartitions,
                                        PartitionIndicesContainerType const& NodesAllPartitions,
                                        PartitionIndicesContainerType const& ElementsAllPartitions,
                                        PartitionIndicesContainerType const& ConditionsAllPartitions)
{
    KRATOS_TRY
    ResetInput();
    std::string word;
    OutputFilesContainerType output_files;

    for(SizeType i = 0 ; i < NumberOfPartitions ; i++)
    {
        std::stringstream buffer;
        buffer << mBaseFilename << "_" << i << ".mdpa";
        std::ofstream* p_ofstream = new std::ofstream(buffer.str().c_str());
        KRATOS_ERROR_IF_NOT(!(*p_ofstream)) << "Error opening mdpa file : " << buffer.str() << std::endl;

        output_files.push_back(p_ofstream);
    }

    while(true)
    {
        ReadWord(word);
        if(mpStream->eof())
            break;
        ReadBlockName(word);
        if(word == "ModelPartData")
            DivideModelPartDataBlock(output_files);
        else if(word == "Table")
            DivideTableBlock(output_files);
        else if(word == "Properties")
            DividePropertiesBlock(output_files);
        else if(word == "Nodes")
            DivideNodesBlock(output_files, NodesAllPartitions);
        else if(word == "Elements")
            DivideElementsBlock(output_files, ElementsAllPartitions);
        else if(word == "Conditions")
            DivideConditionsBlock(output_files, ConditionsAllPartitions);
        else if(word == "NodalData")
            DivideNodalDataBlock(output_files, NodesAllPartitions);
        else if(word == "ElementalData")
            DivideElementalDataBlock(output_files, ElementsAllPartitions);
        else if(word == "ConditionalData")
            DivideConditionalDataBlock(output_files, ConditionsAllPartitions);
        else if (word == "Mesh")
            DivideMeshBlock(output_files, NodesAllPartitions, ElementsAllPartitions, ConditionsAllPartitions);
        else if (word == "SubModelPart")
            DivideSubModelPartBlock(output_files, NodesAllPartitions, ElementsAllPartitions, ConditionsAllPartitions);

    }

    WritePartitionIndices(output_files, NodesPartitions, NodesAllPartitions);

    WriteCommunicatorData(output_files, NumberOfPartitions, DomainsColoredGraph, NodesPartitions, ElementsPartitions, ConditionsPartitions, NodesAllPartitions, ElementsAllPartitions, ConditionsAllPartitions);
    KRATOS_INFO("ModelPartIO") << "  [Total Lines Read : " << mNumberOfLines<<"]" << std::endl;

    for(SizeType i = 0 ; i < NumberOfPartitions ; i++)
        delete output_files[i];
    KRATOS_CATCH("")
}

void ModelPartIO::DivideInputToPartitions(
    Kratos::shared_ptr<std::iostream> * Streams,
    SizeType NumberOfPartitions, GraphType const& DomainsColoredGraph,
    PartitionIndicesType const& NodesPartitions,
    PartitionIndicesType const& ElementsPartitions,
    PartitionIndicesType const& ConditionsPartitions,
    PartitionIndicesContainerType const& NodesAllPartitions,
    PartitionIndicesContainerType const& ElementsAllPartitions,
    PartitionIndicesContainerType const& ConditionsAllPartitions) {

    KRATOS_TRY
    ResetInput();
    std::string word;
    OutputFilesContainerType output_files;

    for(SizeType i = 0 ; i < NumberOfPartitions ; i++)
    {
        output_files.push_back(static_cast<std::ostream *>(&*Streams[i]));
    }

    while(true)
    {
        ReadWord(word);
        if(mpStream->eof())
            break;
        ReadBlockName(word);
        if(word == "ModelPartData")
            DivideModelPartDataBlock(output_files);
        else if(word == "Table")
            DivideTableBlock(output_files);
        else if(word == "Properties")
            DividePropertiesBlock(output_files);
        else if(word == "Nodes")
            DivideNodesBlock(output_files, NodesAllPartitions);
        else if(word == "Elements")
            DivideElementsBlock(output_files, ElementsAllPartitions);
        else if(word == "Conditions")
            DivideConditionsBlock(output_files, ConditionsAllPartitions);
        else if(word == "NodalData")
            DivideNodalDataBlock(output_files, NodesAllPartitions);
        else if(word == "ElementalData")
            DivideElementalDataBlock(output_files, ElementsAllPartitions);
        else if(word == "ConditionalData")
            DivideConditionalDataBlock(output_files, ConditionsAllPartitions);
        else if (word == "Mesh")
            DivideMeshBlock(output_files, NodesAllPartitions, ElementsAllPartitions, ConditionsAllPartitions);
        else if (word == "SubModelPart")
            DivideSubModelPartBlock(output_files, NodesAllPartitions, ElementsAllPartitions, ConditionsAllPartitions);

    }

    WritePartitionIndices(output_files, NodesPartitions, NodesAllPartitions);

    WriteCommunicatorData(output_files, NumberOfPartitions, DomainsColoredGraph, NodesPartitions, ElementsPartitions, ConditionsPartitions, NodesAllPartitions, ElementsAllPartitions, ConditionsAllPartitions);
    KRATOS_INFO("ModelPartIO") << "  [Total Lines Read : " << mNumberOfLines<<"]" << std::endl;

    // for(SizeType i = 0 ; i < NumberOfPartitions ; i++)
    //     delete output_files[i];
    KRATOS_CATCH("")
}

std::string& ModelPartIO::ReadBlockName(std::string& rBlockName)
{
    KRATOS_TRY

    CheckStatement("Begin", rBlockName);
    ReadWord(rBlockName);

    return rBlockName;

    KRATOS_CATCH("")
}


void ModelPartIO::SkipBlock(std::string const& BlockName)
{
    KRATOS_TRY

    std::string word;
    int number_of_nested_blocks = 0;

    while(!mpStream->eof())
    {
        ReadWord(word);
        if(word == "End")
        {
            ReadWord(word);
            if(number_of_nested_blocks == 0){
                    if(CheckStatement(word , BlockName))
                        break;
                }
        }
        else if(word == "Begin")
        {
            number_of_nested_blocks++;
        }
    }

    KRATOS_CATCH("")
}

bool ModelPartIO::CheckEndBlock(std::string const& BlockName, std::string& rWord)
{
    if(rWord == "End")
    {
        ReadWord(rWord);
        CheckStatement(BlockName, rWord);
        return true;
    }

    return false;
}

void ModelPartIO::ReadSubModelPartDataBlock(ModelPart& rModelPart)
{
    KRATOS_TRY

    ReadModelPartDataBlock(rModelPart, true);

    KRATOS_CATCH("")
}

void ModelPartIO::ReadModelPartDataBlock(ModelPart& rModelPart, const bool is_submodelpart)
{
    KRATOS_TRY

    std::string variable_name;

    while(!mpStream->eof())
    {
        ReadWord(variable_name);
        if(!is_submodelpart){
            if(CheckEndBlock("ModelPartData", variable_name))
                break;
        }
        else {
            if(CheckEndBlock("SubModelPartData", variable_name))
                break;
        }
        if(KratosComponents<Variable<double> >::Has(variable_name))
        {
            std::string value;
            double temp;

            ReadWord(value); // reading value
            ExtractValue(value,temp);
            rModelPart[KratosComponents<Variable<double> >::Get(variable_name)] = temp;
        }
        else if(KratosComponents<Variable<bool> >::Has(variable_name))
        {
            std::string value;
            bool temp;

            ReadWord(value); // reading value
            ExtractValue(value,temp);
            rModelPart[KratosComponents<Variable<bool> >::Get(variable_name)] = temp;
        }
        else if(KratosComponents<Variable<int> >::Has(variable_name))
        {
            std::string value;
            int temp;

            ReadWord(value); // reading value
            ExtractValue(value,temp);
            rModelPart[KratosComponents<Variable<int> >::Get(variable_name)] = temp;
        }
        else if(KratosComponents<Variable<array_1d<double, 3> > >::Has(variable_name))
        {
            Vector temp_vector; // defining a Vector because for array_1d the operator >> is not defined yet!
            ReadVectorialValue(temp_vector);
            rModelPart[KratosComponents<Variable<array_1d<double,3> > >::Get(variable_name)] = temp_vector;
        }
        else if(KratosComponents<Variable<Quaternion<double> > >::Has(variable_name))
        {
            Vector temp_vector; // defining a Vector because for Quaternion the operator >> is not defined yet!
            ReadVectorialValue(temp_vector);
            rModelPart[KratosComponents<Variable<Quaternion<double> > >::Get(variable_name)] = temp_vector;
        }
        else if(KratosComponents<Variable<Matrix> >::Has(variable_name))
        {
            ReadVectorialValue(rModelPart[KratosComponents<Variable<Matrix> >::Get(variable_name)]);
        }
        else if(KratosComponents<Variable<std::string> >::Has(variable_name))
        {
            std::string value;
    std::string  temp;

            ReadWord(value); // reading value
            ExtractValue(value,temp);
            rModelPart[KratosComponents<Variable<std::string> >::Get(variable_name)] = temp;
        }
        else
        {
            std::stringstream buffer;
            buffer << variable_name << " is not a valid variable!!!" << std::endl;
            buffer << " [Line " << mNumberOfLines << " ]";
            KRATOS_ERROR << buffer.str() << std::endl;;
        }


    }

    KRATOS_CATCH("")
}

void ModelPartIO::WriteModelPartDataBlock(ModelPart& rModelPart, const bool is_submodelpart)
{
    KRATOS_TRY;

    (*mpStream) << "Begin ModelPartData" << std::endl;
    // TODO: Finish me!!!!!
    (*mpStream) << "End ModelPartData" << std::endl;

    KRATOS_CATCH("");
}


template<class TablesContainerType>
void ModelPartIO::ReadTableBlock(TablesContainerType& rTables)
{
    KRATOS_TRY

    ModelPart::TableType temp_table;

    //SizeType table_id;
    std::string word;

    std::string variable_name;
    ReadWord(variable_name);

    if(!KratosComponents<VariableData>::Has(variable_name))

    {
        std::stringstream buffer;
        buffer << variable_name << " is not a valid argument variable!!! Table only accepts double arguments." << std::endl;
        buffer << " [Line " << mNumberOfLines << " ]";
        KRATOS_ERROR << buffer.str() << std::endl;;

    }

    VariableData const& r_x_variable = KratosComponents<VariableData>::Get(variable_name);

    ReadWord(variable_name);

    if(!KratosComponents<VariableData>::Has(variable_name))
    {
        std::stringstream buffer;
        buffer << variable_name << " is not a valid value variable!!! Table only accepts double values." << std::endl;
        buffer << " [Line " << mNumberOfLines << " ]";
        KRATOS_ERROR << buffer.str() << std::endl;;

    }
    VariableData const& r_y_variable = KratosComponents<VariableData>::Get(variable_name);

    while(!mpStream->eof())
    {
        double x;
        double y;
        ReadWord(word);
        if(CheckEndBlock("Table", word))
            break;

        ExtractValue(word, x);
        ReadWord(word);
        ExtractValue(word, y);

        temp_table.insert(x,y);
    }

    rTables.SetTable(r_x_variable, r_y_variable, temp_table);

    KRATOS_CATCH("")
}

void ModelPartIO::ReadTableBlock(ModelPart::TablesContainerType& rTables)
{
    KRATOS_TRY

    ModelPart::TableType temp_table;

    SizeType table_id;
    std::string word;

    ReadWord(word);
    ExtractValue(word, table_id);

    std::string variable_name;

    ReadWord(variable_name);

    //if(!KratosComponents<Variable<double> >::Has(variable_name))
    //{
    //    std::stringstream buffer;
    //    buffer << variable_name << " is not a valid argument variable!!! Table only accepts double arguments." << std::endl;
    //    buffer << " [Line " << mNumberOfLines << " ]";
    //    KRATOS_ERROR << buffer.str() << std::endl;;

    //}

    ReadWord(variable_name);

    //if(!KratosComponents<Variable<double> >::Has(variable_name))
    //{
    //    std::stringstream buffer;
    //    buffer << variable_name << " is not a valid value variable!!! Table only accepts double values." << std::endl;
    //    buffer << " [Line " << mNumberOfLines << " ]";
    //    KRATOS_ERROR << buffer.str() << std::endl;;

    //}

    while(!mpStream->eof())
    {
        double x;
        double y;
        ReadWord(word);
        if(CheckEndBlock("Table", word))
            break;

        ExtractValue(word, x);
        ReadWord(word);
        ExtractValue(word, y);

        temp_table.insert(x,y);
    }

    rTables.insert(table_id, temp_table);

    KRATOS_CATCH("")
}

template<class TablesContainerType>
void ModelPartIO::WriteTableBlock(TablesContainerType& rTables)
{
    std::string variable1, variable2; // NOTE: NOT POSSIBLE TO KNOW

    // "SOLUTION" // FIXME: We need to think about the how to make possible to define the variables described in the table
    variable1 = "DISTANCE";
    variable2 = "TIME";

    auto numTables = rTables.end() - rTables.begin();

    for(unsigned int i = 0; i < numTables; i++)
    {
        auto itTable = rTables.begin() + i;

        const auto& Data = itTable->Data();
        std::size_t size = Data.size();

        (*mpStream) << "Begin Table" << i << "\t" << variable1 << "\t" << variable2 << " // NOTE: The variables does not correspond with the real ones. Right now the KRATOS Table's does not store the variables"<< std::endl;

        for(std::size_t j = 1 ; j < size ; j++)
        {
            const auto X = Data[j].first;
            const auto Y = (Data[j].second)[0];

            (*mpStream) << X << "\t" << Y << std::endl;
        }

        (*mpStream) << "End Table" << std::endl;
    }
}

void ModelPartIO::WriteTableBlock(ModelPart::TablesContainerType& rTables)
{
    std::string variable1, variable2; // NOTE: NOT POSSIBLE TO KNOW

    // "SOLUTION" // FIXME: We need to think about the how to make possible to define the variables described in the table
    variable1 = "DISTANCE";
    variable2 = "TIME";

    auto numTables = rTables.end() - rTables.begin();

    for(unsigned int i = 0; i < numTables; i++)
    {
        auto itTable = rTables.begin() + i;

        const auto& Data = itTable->Data();
        std::size_t size = Data.size();

        (*mpStream) << "Begin Table" << i << "\t" << variable1 << "\t" << variable2 << " // NOTE: The variables does not correspond with the real ones. Right now the KRATOS Table's does not store the variables"<< std::endl;

        for(std::size_t j = 1 ; j < size ; j++)
        {
            const auto X = Data[j].first;
            const auto Y = (Data[j].second)[0];

            (*mpStream) << X << "\t" << Y << std::endl;
        }

        (*mpStream) << "End Table" << std::endl;
    }
}

void ModelPartIO::ReadNodesBlock(NodesContainerType& rThisNodes)
{
KRATOS_TRY

    SizeType temp_id;
    double x,y,z;
    std::string word;

    SizeType number_of_nodes_read = 0;

    KRATOS_INFO("ModelPartIO") << "  [Reading Nodes    : ";

    while(!mpStream->eof())
    {
        ReadWord(word);
        if(CheckEndBlock("Nodes", word))
            break;

        ExtractValue(word, temp_id);
        ReadWord(word);
        ExtractValue(word, x);
        ReadWord(word);
        ExtractValue(word, y);
        ReadWord(word);
        ExtractValue(word, z);
        NodeType::Pointer temp_node = Kratos::make_shared< NodeType >( ReorderedNodeId(temp_id), x, y, z);
        temp_node->X0() = temp_node->X();
        temp_node->Y0() = temp_node->Y();
        temp_node->Z0() = temp_node->Z();

        rThisNodes.push_back(temp_node);
        number_of_nodes_read++;
    }
    KRATOS_INFO("") << number_of_nodes_read << " nodes read]" << std::endl;

    unsigned int numer_of_nodes_read = rThisNodes.size();
    rThisNodes.Unique();
    KRATOS_WARNING_IF("ModelPartIO", rThisNodes.size() != numer_of_nodes_read) << "attention! we read " << numer_of_nodes_read << " but there are only " << rThisNodes.size() << " non repeated nodes" << std::endl;

    KRATOS_CATCH("")
}

void ModelPartIO::ReadNodesBlock(ModelPart& rModelPart)
{
    KRATOS_TRY
/*
NodeType temp_node;
    SizeType temp_id;

    // Giving model part's variables list to the node
    temp_node.SetSolutionStepVariablesList(&rModelPart.GetNodalSolutionStepVariablesList());

    //set buffer size
    temp_node.SetBufferSize(rModelPart.GetBufferSize());


    std::string word;

    SizeType number_of_nodes_read = 0;

    KRATOS_INFO("ModelPartIO") << "  [Reading Nodes    : ";

    while(!mpStream->eof())
    {
        ReadWord(word);
        if(CheckEndBlock("Nodes", word))
            break;

        ExtractValue(word, temp_id);
        temp_node.SetId(ReorderedNodeId(temp_id));
        ReadWord(word);
        ExtractValue(word, temp_node.X());
        ReadWord(word);
        ExtractValue(word, temp_node.Y());
        ReadWord(word);
        ExtractValue(word, temp_node.Z());

        temp_node.X0() = temp_node.X();
        temp_node.Y0() = temp_node.Y();
        temp_node.Z0() = temp_node.Z();


        rModelPart.Nodes().push_back(temp_node);
        number_of_nodes_read++;
    }
    KRATOS_INFO("") << number_of_nodes_read << " nodes read]" << std::endl;

    unsigned int numer_of_nodes_read = rModelPart.Nodes().size();
    rModelPart.Nodes().Unique();
    KRATOS_WARNING_IF("ModelPartIO", rModelPart.Nodes().size() != numer_of_nodes_read) << "attention! we read " << numer_of_nodes_read << " but there are only " << rModelPart.Nodes().size() << " non repeated nodes" << std::endl;
*/
SizeType id;
double x;
double y;
double z;

std::string word;

SizeType number_of_nodes_read = 0;
    const unsigned int old_size = rModelPart.Nodes().size();

    typedef std::map< unsigned int, array_1d<double,3> > map_type;
    map_type read_coordinates;

    KRATOS_INFO("ModelPartIO") << "  [Reading Nodes    : ";

while(!mpStream->eof())
{
    ReadWord(word);
    if(CheckEndBlock("Nodes", word))
    break;

    ExtractValue(word, id);
    ReadWord(word);
    ExtractValue(word, x);
    ReadWord(word);
    ExtractValue(word, y);
    ReadWord(word);
    ExtractValue(word, z);

        array_1d<double,3> coords;
        coords[0]=x;
        coords[1]=y;
        coords[2]=z;
        read_coordinates[ReorderedNodeId(id)] = coords;
    number_of_nodes_read++;
}


//make this to construct the nodes "in parallel" - the idea is that first touch is being done in parallel but the reading is actually sequential
const int nnodes = read_coordinates.size();
const int nthreads = OpenMPUtils::GetNumThreads();
    std::vector<int> partition;
    OpenMPUtils::DivideInPartitions(nnodes, nthreads, partition);

    map_type::const_iterator it = read_coordinates.begin();
    for(int ithread=0; ithread<nthreads; ithread++)
    {
        #pragma omp parallel
        {
            //note that the reading is only done by one of the threads
            if(OpenMPUtils::ThisThread() == ithread)
            {
                for(int i=partition[ithread]; i<partition[ithread+1]; i++)
                {
                    const unsigned int node_id = it->first;
                    const array_1d<double,3>& coords = it->second;
                    rModelPart.CreateNewNode(node_id,coords[0],coords[1],coords[2]);
                    it++;
                }
            }
        }
    }

KRATOS_INFO("") << number_of_nodes_read << " nodes read]" << std::endl;
KRATOS_WARNING_IF("ModelPartIO", rModelPart.Nodes().size() - old_size != number_of_nodes_read) << "attention! we read " << number_of_nodes_read << " but there are only " << rModelPart.Nodes().size() - old_size<< " non repeated nodes" << std::endl;

    KRATOS_CATCH("")
}

std::size_t ModelPartIO::CountNodesInBlock()
{
    KRATOS_TRY;

    std::vector<SizeType> found_ids;

    SizeType temp_id;

    std::string word;

    SizeType number_of_nodes_read = 0;

//KRATOS_INFO("ModelPartIO") << "  [Reading Nodes    : ";

    while(!mpStream->eof())
    {
        ReadWord(word);
        if(CheckEndBlock("Nodes", word))
            break;

        ExtractValue(word, temp_id);
        found_ids.push_back(temp_id);

        ReadWord(word); // skip X coordinate
        ReadWord(word); // skip Y
        ReadWord(word); // skip Z

        number_of_nodes_read++;
    }
    //KRATOS_INFO("") << number_of_nodes_read << " nodes read]" << std::endl;

    // Error check: look for duplicate nodes
    std::sort(found_ids.begin(),found_ids.end());
    std::vector<std::size_t>::iterator unique_end = std::unique(found_ids.begin(),found_ids.end());
    std::size_t number_of_unique_nodes = std::distance(found_ids.begin(),unique_end);

    KRATOS_WARNING_IF("ModelPartIO", number_of_unique_nodes != number_of_nodes_read) << "attention! we read " << number_of_nodes_read << " but there are only " << number_of_unique_nodes << " non repeated nodes" << std::endl;

    return number_of_nodes_read;

    KRATOS_CATCH("");
}

void ModelPartIO::ReadPropertiesBlock(PropertiesContainerType& rThisProperties)
{
    KRATOS_TRY

    Properties::Pointer props = Kratos::make_shared<Properties>();
    Properties& temp_properties = *props;
    //Properties temp_properties;

    std::string word;
    std::string variable_name;

    SizeType temp_properties_id;

    ReadWord(word);
    ExtractValue(word, temp_properties_id);
    temp_properties.SetId(temp_properties_id);

    while(!mpStream->eof())
    {
        ReadWord(variable_name);
        if(CheckEndBlock("Properties", variable_name))
            break;

        if(variable_name == "Begin") // here we have some nested block.
        {
            ReadBlockName(variable_name);
            if(variable_name == "Table") // At this moment the only supported nested block is a table
                ReadTableBlock(temp_properties);
        }
        else if(KratosComponents<Variable<std::string> >::Has(variable_name))
        {
            std::string value;
            std::string  temp;

            ReadWord(value); // reading value
            ExtractValue(value,temp);
            temp_properties[KratosComponents<Variable<std::string> >::Get(variable_name)] = temp;
        }
        else if(KratosComponents<Variable<double> >::Has(variable_name))
        {
            std::string value;
            double temp;

            ReadWord(value); // reading value
            ExtractValue(value,temp);
            temp_properties[KratosComponents<Variable<double> >::Get(variable_name)] = temp;
        }
        else if(KratosComponents<Variable<int> >::Has(variable_name))
        {
            std::string value;
            int temp;

            ReadWord(value); // reading value
            ExtractValue(value,temp);
            temp_properties[KratosComponents<Variable<int> >::Get(variable_name)] = temp;
        }
        else if(KratosComponents<Variable<bool> >::Has(variable_name))
        {
            std::string value;
            bool temp;

            ReadWord(value); // reading value
            ExtractValue(value,temp);
            temp_properties[KratosComponents<Variable<bool> >::Get(variable_name)] = temp;
        }
        else if(KratosComponents<Variable<array_1d<double, 3> > >::Has(variable_name))
        {
            Vector temp_vector; // defining a Vector because for array_1d the operator >> is not defined yet!
            ReadVectorialValue(temp_vector);
            temp_properties[KratosComponents<Variable<array_1d<double,3> > >::Get(variable_name)] = temp_vector;
        }
        else if(KratosComponents<Variable<Quaternion<double> > >::Has(variable_name))
        {
            Vector temp_vector; // defining a Vector because for Quaternion the operator >> is not defined yet!
            ReadVectorialValue(temp_vector);
            temp_properties[KratosComponents<Variable<Quaternion<double> > >::Get(variable_name)] = temp_vector;
        }
        else if(KratosComponents<Variable<Vector> >::Has(variable_name))
        {
            ReadVectorialValue(temp_properties[KratosComponents<Variable<Vector> >::Get(variable_name)]);
        }
        else if(KratosComponents<Variable<Matrix> >::Has(variable_name))
        {
            ReadVectorialValue(temp_properties[KratosComponents<Variable<Matrix> >::Get(variable_name)]);
        }
        else if(KratosComponents<Variable<ConstitutiveLaw::Pointer> >::Has(variable_name))
        {
            ReadConstitutiveLawValue(temp_properties[KratosComponents<Variable<ConstitutiveLaw::Pointer> >::Get(variable_name)]);
        }
        else
        {
            std::stringstream buffer;
            buffer << variable_name << " is not a valid variable!!!" << std::endl;
            buffer << " [Line " << mNumberOfLines << " ]";
            KRATOS_ERROR << buffer.str() << std::endl;;
        }

    }

    rThisProperties.push_back(props);
//         rThisProperties.push_back(temp_properties);

    KRATOS_CATCH("")
}

void ModelPartIO::ReadElementsBlock(ModelPart& rModelPart)
{
    KRATOS_TRY

    SizeType id;
    SizeType properties_id;
    SizeType node_id;
    SizeType number_of_read_elements = 0;


    std::string word;
    std::string element_name;

    ReadWord(element_name);
    KRATOS_INFO("ModelPartIO") << "  [Reading Elements : ";

    if(!KratosComponents<Element>::Has(element_name))
    {
        std::stringstream buffer;
        buffer << "Element " << element_name << " is not registered in Kratos.";
        buffer << " Please check the spelling of the element name and see if the application which containing it, is registered corectly.";
        buffer << " [Line " << mNumberOfLines << " ]";
        KRATOS_ERROR << buffer.str() << std::endl;;
        return;
    }

    Element const& r_clone_element = KratosComponents<Element>::Get(element_name);
    SizeType number_of_nodes = r_clone_element.GetGeometry().size();
    Element::NodesArrayType temp_element_nodes;
    ModelPart::ElementsContainerType aux_elements;

    while(!mpStream->eof())
    {
        ReadWord(word); // Reading the element id or End
        if(CheckEndBlock("Elements", word))
            break;

        ExtractValue(word,id);
        ReadWord(word); // Reading the properties id;
        ExtractValue(word, properties_id);
        Properties::Pointer p_temp_properties = *(FindKey(rModelPart.rProperties(), properties_id, "Properties").base());
        temp_element_nodes.clear();
        for(SizeType i = 0 ; i < number_of_nodes ; i++)
        {
            ReadWord(word); // Reading the node id;
            ExtractValue(word, node_id);
            temp_element_nodes.push_back( *(FindKey(rModelPart.Nodes(), ReorderedNodeId(node_id), "Node").base()));
        }

        aux_elements.push_back(r_clone_element.Create(ReorderedElementId(id), temp_element_nodes, p_temp_properties));
        number_of_read_elements++;

    }
    KRATOS_INFO("") << number_of_read_elements << " elements read] [Type: " <<element_name << "]" << std::endl;
    aux_elements.Unique();

    rModelPart.AddElements(aux_elements.begin(), aux_elements.end());

    KRATOS_CATCH("")
}

void ModelPartIO::ReadElementsBlock(NodesContainerType& rThisNodes, PropertiesContainerType& rThisProperties, ElementsContainerType& rThisElements)
{
    KRATOS_TRY

    SizeType id;
    SizeType properties_id;
    SizeType node_id;
    SizeType number_of_read_elements = 0;


    std::string word;
    std::string element_name;

    ReadWord(element_name);
    KRATOS_INFO("ModelPartIO") << "  [Reading Elements : ";

    if(!KratosComponents<Element>::Has(element_name))
    {
        std::stringstream buffer;
        buffer << "Element " << element_name << " is not registered in Kratos.";
        buffer << " Please check the spelling of the element name and see if the application which containing it, is registered corectly.";
        buffer << " [Line " << mNumberOfLines << " ]";
        KRATOS_ERROR << buffer.str() << std::endl;;
        return;
    }

    Element const& r_clone_element = KratosComponents<Element>::Get(element_name);
    SizeType number_of_nodes = r_clone_element.GetGeometry().size();
    Element::NodesArrayType temp_element_nodes;


    while(!mpStream->eof())
    {
        ReadWord(word); // Reading the element id or End
        if(CheckEndBlock("Elements", word))
            break;

        ExtractValue(word,id);
        ReadWord(word); // Reading the properties id;
        ExtractValue(word, properties_id);
        Properties::Pointer p_temp_properties = *(FindKey(rThisProperties, properties_id, "Properties").base());
        temp_element_nodes.clear();
        for(SizeType i = 0 ; i < number_of_nodes ; i++)
        {
            ReadWord(word); // Reading the node id;
            ExtractValue(word, node_id);
            temp_element_nodes.push_back( *(FindKey(rThisNodes, ReorderedNodeId(node_id), "Node").base()));
        }

        rThisElements.push_back(r_clone_element.Create(ReorderedElementId(id), temp_element_nodes, p_temp_properties));
        number_of_read_elements++;

    }
    KRATOS_INFO("") << number_of_read_elements << " elements read] [Type: " <<element_name << "]" << std::endl;
    rThisElements.Unique();

    KRATOS_CATCH("")
}


void ModelPartIO::ReadConditionsBlock(ModelPart& rModelPart)
{
    ReadConditionsBlock(rModelPart.Nodes(), rModelPart.rProperties(), rModelPart.Conditions());
}

void ModelPartIO::ReadConditionsBlock(NodesContainerType& rThisNodes, PropertiesContainerType& rThisProperties, ConditionsContainerType& rThisConditions)
{
    KRATOS_TRY

    SizeType id;
    SizeType properties_id;
    SizeType node_id;
    SizeType number_of_read_conditions = 0;


    std::string word;
    std::string condition_name;

    ReadWord(condition_name);
    KRATOS_INFO("ModelPartIO") << "  [Reading Conditions : ";

    if(!KratosComponents<Condition>::Has(condition_name))
    {
        std::stringstream buffer;
        buffer << "Condition " << condition_name << " is not registered in Kratos.";
        buffer << " Please check the spelling of the condition name and see if the application containing it is registered corectly.";
        buffer << " [Line " << mNumberOfLines << " ]";
        KRATOS_ERROR << buffer.str() << std::endl;;
        return;
    }

    Condition const& r_clone_condition = KratosComponents<Condition>::Get(condition_name);
    SizeType number_of_nodes = r_clone_condition.GetGeometry().size();
    Condition::NodesArrayType temp_condition_nodes;

    while(!mpStream->eof())
    {
        ReadWord(word); // Reading the condition id or End
        if(CheckEndBlock("Conditions", word))
            break;

        ExtractValue(word,id);
        ReadWord(word); // Reading the properties id;
        ExtractValue(word, properties_id);
        Properties::Pointer p_temp_properties = *(FindKey(rThisProperties, properties_id, "Properties").base());
        temp_condition_nodes.clear();
        for(SizeType i = 0 ; i < number_of_nodes ; i++)
        {
            ReadWord(word); // Reading the node id;
            ExtractValue(word, node_id);
            temp_condition_nodes.push_back( *(FindKey(rThisNodes, ReorderedNodeId(node_id), "Node").base()));
        }

        rThisConditions.push_back(r_clone_condition.Create(ReorderedConditionId(id), temp_condition_nodes, p_temp_properties));
        number_of_read_conditions++;
    }
    KRATOS_INFO("") << number_of_read_conditions << " conditions read] [Type: " << condition_name << "]" << std::endl;
    rThisConditions.Unique();

    KRATOS_CATCH("")
}


void ModelPartIO::ReadNodalDataBlock(ModelPart& rThisModelPart)
{
    KRATOS_TRY

    NodesContainerType& rThisNodes = rThisModelPart.Nodes();

    typedef VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > array_1d_component_type;

    std::string variable_name;

    ReadWord(variable_name);

    VariablesList rThisVariables = rThisModelPart.GetNodalSolutionStepVariablesList();


    if(KratosComponents<Flags >::Has(variable_name))
    {
        ReadNodalFlags(rThisNodes, static_cast<Flags const& >(KratosComponents<Flags >::Get(variable_name)));
    }
    else if(KratosComponents<Variable<int> >::Has(variable_name))
    {
        const bool has_been_added = rThisVariables.Has(KratosComponents<Variable<int> >::Get(variable_name)) ;
        if( !has_been_added && mOptions.Is(IGNORE_VARIABLES_ERROR) ) {
            KRATOS_WARNING("ModelPartIO") <<"WARNING: Skipping NodalData block. Variable "<<variable_name<<" has not been added to ModelPart '"<<rThisModelPart.Name()<<"'"<<std::endl<<std::endl;
            SkipBlock("NodalData");
        }
        else {
            KRATOS_ERROR_IF_NOT(has_been_added) << "The nodal solution step container does not have this variable: " << variable_name << "." << std::endl;
            ReadNodalScalarVariableData(rThisNodes, static_cast<Variable<int> const& >(KratosComponents<Variable<int> >::Get(variable_name)));
        }
    }
    else if(KratosComponents<Variable<double> >::Has(variable_name))
    {
        const bool has_been_added = rThisVariables.Has(KratosComponents<Variable<double> >::Get(variable_name)) ;
        if( !has_been_added && mOptions.Is(IGNORE_VARIABLES_ERROR) ) {
            KRATOS_WARNING("ModelPartIO")<<"WARNING: Skipping NodalData block. Variable "<<variable_name<<" has not been added to ModelPart '"<<rThisModelPart.Name()<<"'"<<std::endl<<std::endl;
            SkipBlock("NodalData");
        }
        else {
            KRATOS_ERROR_IF_NOT(has_been_added) << "The nodal solution step container does not have this variable: " << variable_name << "." << std::endl;
            ReadNodalDofVariableData(rThisNodes, static_cast<Variable<double> const& >(KratosComponents<Variable<double> >::Get(variable_name)));
        }
    }
    else if(KratosComponents<array_1d_component_type>::Has(variable_name))
    {
        const bool has_been_added = rThisVariables.Has(KratosComponents<array_1d_component_type>::Get(variable_name).GetSourceVariable()) ;
        if( !has_been_added && mOptions.Is(IGNORE_VARIABLES_ERROR) ) {
            KRATOS_WARNING("ModelPartIO") <<"WARNING: Skipping NodalData block. Variable "<<variable_name<<" has not been added to ModelPart '"<<rThisModelPart.Name()<<"'"<<std::endl<<std::endl;
            SkipBlock("NodalData");
        }
        else {
            KRATOS_ERROR_IF_NOT(has_been_added) << "The nodal solution step container does not have this variable: " << variable_name << "." << std::endl;
            ReadNodalDofVariableData(rThisNodes, static_cast<array_1d_component_type const& >(KratosComponents<array_1d_component_type>::Get(variable_name)));
        }
    }
    else if(KratosComponents<Variable<array_1d<double, 3> > >::Has(variable_name))
    {
        const bool has_been_added = rThisVariables.Has(KratosComponents<Variable<array_1d<double, 3> > >::Get(variable_name)) ;
        if( !has_been_added && mOptions.Is(IGNORE_VARIABLES_ERROR) ) {
            KRATOS_WARNING("ModelPartIO")<<"WARNING: Skipping NodalData block. Variable "<<variable_name<<" has not been added to ModelPart '"<<rThisModelPart.Name()<<"'"<<std::endl<<std::endl;
        }
        else {
            KRATOS_ERROR_IF_NOT(has_been_added) << "The nodal solution step container does not have this variable: " << variable_name << "." << std::endl;
            ReadNodalVectorialVariableData(rThisNodes, static_cast<Variable<array_1d<double, 3> > const& >(KratosComponents<Variable<array_1d<double, 3> > >::Get(variable_name)), Vector(3));
        }
    }
    else if(KratosComponents<Variable<Quaternion<double> > >::Has(variable_name))
    {
        const bool has_been_added = rThisVariables.Has(KratosComponents<Variable<Quaternion<double> > >::Get(variable_name)) ;
        if( !has_been_added && mOptions.Is(IGNORE_VARIABLES_ERROR) ) {
            KRATOS_WARNING("ModelPartIO")<<"WARNING: Skipping NodalData block. Variable "<<variable_name<<" has not been added to ModelPart '"<<rThisModelPart.Name()<<"'"<<std::endl<<std::endl;
        }
        else {
            KRATOS_ERROR_IF_NOT(has_been_added) << "The nodal solution step container does not have this variable: " << variable_name << "." << std::endl;
            ReadNodalVectorialVariableData(rThisNodes, static_cast<Variable<Quaternion<double> > const& >(KratosComponents<Variable<Quaternion<double> > >::Get(variable_name)), Vector(4));
        }
    }
    else if(KratosComponents<Variable<Matrix> >::Has(variable_name))
    {
        ReadNodalVectorialVariableData(rThisNodes, static_cast<Variable<Matrix > const& >(KratosComponents<Variable<Matrix> >::Get(variable_name)), Matrix(3,3));
    }
    else if(KratosComponents<Variable<Vector> >::Has(variable_name))
    {
        ReadNodalVectorialVariableData(rThisNodes, static_cast<Variable<Vector > const& >(KratosComponents<Variable<Vector> >::Get(variable_name)), Vector(3));
    }
    else if(KratosComponents<VariableData>::Has(variable_name))
    {
        std::stringstream buffer;
        buffer << variable_name << " is not supported to be read by this IO or the type of variable is not registered correctly" << std::endl;
        buffer << " [Line " << mNumberOfLines << " ]";
        KRATOS_ERROR << buffer.str() << std::endl;;
    }
    else
    {
        std::stringstream buffer;
        buffer << variable_name << " is not a valid variable!!!" << std::endl;
        buffer << " [Line " << mNumberOfLines << " ]";
        KRATOS_ERROR << buffer.str() << std::endl;;
    }

    KRATOS_CATCH("")
}

void ModelPartIO::WriteNodalDataBlock(ModelPart& rThisModelPart)
{
    KRATOS_TRY

    VariablesList r_this_variables = rThisModelPart.GetNodalSolutionStepVariablesList();
    NodesContainerType& r_this_nodes = rThisModelPart.Nodes();

    typedef VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > array_1d_component_type;

    std::string variable_name;

    // FIXME: Maybe there is a better way (I get confused with to much KratosComponents)
    for(std::size_t i = 0; i < r_this_variables.size(); i++)
    {
        auto it_var = r_this_variables.begin() + i;

        variable_name = it_var->Name();

        if(KratosComponents<Flags >::Has(variable_name))
        {
            (*mpStream) << "Begin NodalData\t" << variable_name << std::endl;
            auto Variable = static_cast<Flags const& >(KratosComponents<Flags>::Get(variable_name));
            for(std::size_t j = 0; j < r_this_nodes.size(); j++)
            {
                auto it_node = r_this_nodes.begin() + j;

                const bool is_fixed = it_node->Is(Variable); // FIXME: I don't know if they are set this way
                (*mpStream) << it_node->Id() <<"\t" << is_fixed << "\t" << it_node->Is(Variable) << std::endl;
            }
            (*mpStream) << "End NodalData" << std::endl << std::endl;
        }
        else if(KratosComponents<Variable<int>>::Has(variable_name))
        {
            (*mpStream) << "Begin NodalData\t" << variable_name << std::endl;
            auto Variable = static_cast<Kratos::Variable<int> const& >(KratosComponents<Kratos::Variable<int> >::Get(variable_name));
            for(std::size_t j = 0; j < r_this_nodes.size(); j++)
            {
                auto it_node = r_this_nodes.begin() + j;

                const bool is_fixed = it_node->IsFixed(Variable);
                (*mpStream) << it_node->Id() <<"\t" << is_fixed << "\t" << it_node->FastGetSolutionStepValue(Variable, 0) << std::endl;
            }
            (*mpStream) << "End NodalData" << std::endl << std::endl;
        }
        else if(KratosComponents<Variable<double>>::Has(variable_name))
        {
            (*mpStream) << "Begin NodalData\t" << variable_name << std::endl;
            auto Variable = static_cast<Kratos::Variable<double> const& >(KratosComponents<Kratos::Variable<double> >::Get(variable_name));
            for(std::size_t j = 0; j < r_this_nodes.size(); j++)
            {
                auto it_node = r_this_nodes.begin() + j;

                const bool is_fixed = it_node->IsFixed(Variable);
                (*mpStream) << it_node->Id() <<"\t" << is_fixed << "\t" << it_node->FastGetSolutionStepValue(Variable, 0) << std::endl;
            }
            (*mpStream) << "End NodalData" << std::endl << std::endl;
        }
        else if(KratosComponents<array_1d_component_type>::Has(variable_name))
        {
            (*mpStream) << "Begin NodalData\t" << variable_name << std::endl;
            auto Variable = static_cast<array_1d_component_type const& >(KratosComponents<array_1d_component_type >::Get(variable_name));
            for(std::size_t j = 0; j < r_this_nodes.size(); j++)
            {
                auto it_node = r_this_nodes.begin() + j;

                const bool is_fixed = it_node->IsFixed(Variable);
                (*mpStream) << it_node->Id() <<"\t" << is_fixed << "\t" << it_node->FastGetSolutionStepValue(Variable, 0) << std::endl;
            }
            (*mpStream) << "End NodalData" << std::endl << std::endl;
        }
        else if(KratosComponents<Variable<array_1d<double, 3> > >::Has(variable_name))
        {
            if(KratosComponents<array_1d_component_type>::Has(variable_name + "_X")) // To check if it defined by components or as a vector
            {
                (*mpStream) << "Begin NodalData\t" << variable_name << "_X" << std::endl;
                auto VariableX = static_cast<array_1d_component_type const& >(KratosComponents<array_1d_component_type >::Get(variable_name+"_X"));
                for(std::size_t j = 0; j < r_this_nodes.size(); j++)
                {
                    auto it_node = r_this_nodes.begin() + j;

                    const bool is_fixed = it_node->IsFixed(VariableX);
                    (*mpStream) << it_node->Id() <<"\t" << is_fixed << "\t" << it_node->FastGetSolutionStepValue(VariableX, 0) << std::endl;
                }
                (*mpStream) << "End NodalData" << std::endl << std::endl;

                (*mpStream) << "Begin NodalData\t" << variable_name << "_Y" << std::endl;
                auto VariableY = static_cast<array_1d_component_type const& >(KratosComponents<array_1d_component_type >::Get(variable_name+"_Y"));
                for(std::size_t j = 0; j < r_this_nodes.size(); j++)
                {
                    auto it_node = r_this_nodes.begin() + j;

                    const bool is_fixed = it_node->IsFixed(VariableY);
                    (*mpStream) << it_node->Id() <<"\t" << is_fixed << "\t" << it_node->FastGetSolutionStepValue(VariableY, 0) << std::endl;
                }
                (*mpStream) << "End NodalData" << std::endl << std::endl;

                (*mpStream) << "Begin NodalData\t" << variable_name << "_Z" << std::endl;
                auto VariableZ = static_cast<array_1d_component_type const& >(KratosComponents<array_1d_component_type >::Get(variable_name+"_Z"));
                for(std::size_t j = 0; j < r_this_nodes.size(); j++)
                {
                    auto it_node = r_this_nodes.begin() + j;

                    const bool is_fixed = it_node->IsFixed(VariableZ);
                    (*mpStream) << it_node->Id() <<"\t" << is_fixed << "\t" << it_node->FastGetSolutionStepValue(VariableZ, 0) << std::endl;
                }
                (*mpStream) << "End NodalData" << std::endl << std::endl;
            }
            else
            {
                KRATOS_WARNING("ModelPartIO") << variable_name << " is not a valid variable for output!!!" << std::endl;
//                 (*mpStream) << "Begin NodalData\t" << variable_name << std::endl;
//                 auto Variable = KratosComponents<array_1d<double, 3>>::Get(variable_name);
//                 // TODO: Finish me
//                 (*mpStream) << "End NodalData" << std::endl << std::endl;
            }
        }
//             else if(KratosComponents<Variable<Quaternion<double> > >::Has(variable_name))
//             {
//                 (*mpStream) << "Begin NodalData\t" << variable_name << std::endl;
//                 auto Variable = KratosComponents<Quaternion<double>>::Get(variable_name);
//                 // TODO: Finish me
//                 (*mpStream) << "End NodalData" << std::endl << std::endl;
//             }
//             else if(KratosComponents<Variable<Matrix> >::Has(variable_name))
//             {
//                 (*mpStream) << "Begin NodalData\t" << variable_name << std::endl;
//                 auto Variable = KratosComponents<Matrix>::Get(variable_name);
//                 // TODO: Finish me
//                 (*mpStream) << "End NodalData" << std::endl << std::endl;
//             }
//             else if(KratosComponents<Variable<Vector> >::Has(variable_name))
//             {
//                 (*mpStream) << "Begin NodalData\t" << variable_name << std::endl;
//                 auto Variable = KratosComponents<Matrix>::Get(variable_name);
//                 // TODO: Finish me
//                 (*mpStream) << "End NodalData" << std::endl << std::endl;
//             }
        else
        {
            KRATOS_WARNING("ModelPartIO") << variable_name << " is not a valid variable for output!!!" << std::endl;
        }

    }

    KRATOS_CATCH("")
}
template<class TObjectsContainerType>
void ModelPartIO::WriteDataBlock(const TObjectsContainerType& rThisObjectContainer, const std::string& rObjectName){
    typedef VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > array_1d_component_type;
    std::unordered_set<std::string> variables;

    for(auto& object :rThisObjectContainer){
        for(auto& var:object.GetData()){
            auto const& is_included = variables.find(var.first->Name());
            if(is_included == variables.end()){
                variables.insert(var.first->Name());
                // determine variable type
                if(KratosComponents<Variable<bool>>::Has(var.first->Name())){
                    WriteDataBlock<Variable<bool>, TObjectsContainerType>(rThisObjectContainer, var.first, rObjectName);
                }
                else if(KratosComponents<Variable<int>>::Has(var.first->Name())){
                    WriteDataBlock<Variable<int>, TObjectsContainerType>(rThisObjectContainer, var.first, rObjectName);
                }
                else if(KratosComponents<Variable<double>>::Has(var.first->Name())){
                    WriteDataBlock<Variable<double>, TObjectsContainerType>(rThisObjectContainer, var.first, rObjectName);
                }
                else if(KratosComponents<array_1d_component_type>::Has(var.first->Name())){
                    WriteDataBlock<array_1d_component_type>(rThisObjectContainer, var.first, rObjectName);
                }
                else if(KratosComponents<Variable<array_1d<double,3>>>::Has(var.first->Name())){
                    WriteDataBlock<Variable<array_1d<double,3>>, TObjectsContainerType>(rThisObjectContainer, var.first, rObjectName);
                }
                else if(KratosComponents<Variable<Quaternion<double>>>::Has(var.first->Name())){
                    WriteDataBlock<Variable<Quaternion<double>>, TObjectsContainerType>(rThisObjectContainer, var.first, rObjectName);
                }
                else if(KratosComponents<Variable<Vector>>::Has(var.first->Name())){
                    WriteDataBlock<Variable<Vector>, TObjectsContainerType>(rThisObjectContainer, var.first, rObjectName);
                }
                else if(KratosComponents<Variable<Matrix>>::Has(var.first->Name())){
                    WriteDataBlock<Variable<Matrix>, TObjectsContainerType>(rThisObjectContainer, var.first, rObjectName);
                }
                else
                    KRATOS_WARNING("ModelPartIO") << var.first->Name() << " is not a valid variable for output!!!" << std::endl;

            }
        }
    }
}

template<class TVariableType, class TObjectsContainerType>
void ModelPartIO::WriteDataBlock(const TObjectsContainerType& rThisObjectContainer,const VariableData* rVariable, const std::string& rObjectName){
    const TVariableType variable = KratosComponents<TVariableType>::Get(rVariable->Name());
    (*mpStream) << "Begin "<<rObjectName<<"alData "<<variable.Name()<<std::endl;
    for(auto& object : rThisObjectContainer){
        if(object.Has(variable)){
            (*mpStream)<<object.Id()<<"\t"<<object.GetValue(variable)<<std::endl;
        }
    }
    (*mpStream)<<"End "<<rObjectName<<"alData\n"<<std::endl;
}

template<class TVariableType>
void ModelPartIO::ReadNodalDofVariableData(NodesContainerType& rThisNodes, TVariableType& rVariable)
{
    KRATOS_TRY

    SizeType id;
    bool is_fixed;
    double nodal_value;

    std::string value;

    while(!mpStream->eof())
    {
        ReadWord(value); // reading id
        if(CheckEndBlock("NodalData", value))
            break;

        ExtractValue(value, id);
        typename NodesContainerType::iterator it_node = FindKey(rThisNodes, ReorderedNodeId(id), "Node");

        // reading is_fixed
        ReadWord(value);
        ExtractValue(value, is_fixed);
        if(is_fixed)
            it_node->Fix(rVariable);

        // reading nodal_value
        ReadWord(value);
        ExtractValue(value, nodal_value);

        it_node->GetSolutionStepValue(rVariable, 0) =  nodal_value;
    }

    KRATOS_CATCH("")
}

void ModelPartIO::ReadNodalFlags(NodesContainerType& rThisNodes, Flags const& rFlags)
{

    KRATOS_TRY

    SizeType id;

    std::string value;

    while(!mpStream->eof())
    {
        ReadWord(value); // reading id
        if(CheckEndBlock("NodalData", value))
            break;

        ExtractValue(value, id);

        FindKey(rThisNodes, ReorderedNodeId(id), "Node")->Set(rFlags);
    }

    KRATOS_CATCH("")
}

template<class TVariableType>
void ModelPartIO::ReadNodalScalarVariableData(NodesContainerType& rThisNodes, TVariableType& rVariable)
{
    KRATOS_TRY

    SizeType id;
    bool is_fixed;
    typename TVariableType::Type nodal_value;

    std::string value;

    while(!mpStream->eof())
    {
        ReadWord(value); // reading id
        if(CheckEndBlock("NodalData", value))
            break;

        ExtractValue(value, id);

        // reading is_fixed
        ReadWord(value);
        ExtractValue(value, is_fixed);
        if(is_fixed)
        {
            std::stringstream buffer;
            buffer << "Only double variables or components can be fixed.";
            buffer <<  " [Line " << mNumberOfLines << " ]";
            KRATOS_ERROR << buffer.str() << std::endl;;
        }

        // reading nodal_value
        ReadWord(value);
        ExtractValue(value, nodal_value);

        FindKey(rThisNodes, ReorderedNodeId(id), "Node")->GetSolutionStepValue(rVariable, 0) =  nodal_value;
    }

    KRATOS_CATCH("")
}

template<class TVariableType, class TDataType>
void ModelPartIO::ReadNodalVectorialVariableData(NodesContainerType& rThisNodes, TVariableType& rVariable, TDataType Dummy)
{
    KRATOS_TRY

    SizeType id;
    bool is_fixed;
    TDataType nodal_value;

    std::string value;

    while(!mpStream->eof())
    {
        ReadWord(value); // reading id
        if(CheckEndBlock("NodalData", value))
            break;

        ExtractValue(value, id);

        // reading is_fixed
        ReadWord(value);
        ExtractValue(value, is_fixed);
        if(is_fixed)
        {
            std::stringstream buffer;
            buffer << "Only double variables or components can be fixed.";
            buffer <<  " [Line " << mNumberOfLines << " ]";
            KRATOS_ERROR << buffer.str() << std::endl;;
        }

        // reading nodal_value
        ReadVectorialValue(nodal_value);

        FindKey(rThisNodes, ReorderedNodeId(id), "Node")->GetSolutionStepValue(rVariable, 0) =  nodal_value;
    }

    KRATOS_CATCH("")
}

void ModelPartIO::ReadElementalDataBlock(ElementsContainerType& rThisElements)
{
    KRATOS_TRY

    typedef VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > array_1d_component_type;

    std::string variable_name;

    ReadWord(variable_name);

    if(KratosComponents<Variable<bool> >::Has(variable_name))
    {
        ReadElementalScalarVariableData(rThisElements, static_cast<Variable<bool> const& >(KratosComponents<Variable<bool> >::Get(variable_name)));
    }
    else if(KratosComponents<Variable<int> >::Has(variable_name))
    {
        ReadElementalScalarVariableData(rThisElements, static_cast<Variable<int> const& >(KratosComponents<Variable<int> >::Get(variable_name)));
    }
    else if(KratosComponents<Variable<double> >::Has(variable_name))
    {
        ReadElementalScalarVariableData(rThisElements, static_cast<Variable<double> const& >(KratosComponents<Variable<double> >::Get(variable_name)));
    }
    else if(KratosComponents<array_1d_component_type>::Has(variable_name))
    {
        ReadElementalScalarVariableData(rThisElements, static_cast<array_1d_component_type const& >(KratosComponents<array_1d_component_type>::Get(variable_name)));
    }
    else if(KratosComponents<Variable<array_1d<double, 3> > >::Has(variable_name))
    {
        ReadElementalVectorialVariableData(rThisElements, static_cast<Variable<array_1d<double, 3> > const& >(KratosComponents<Variable<array_1d<double, 3> > >::Get(variable_name)), Vector(3));
    }
    else if(KratosComponents<Variable<Quaternion<double> > >::Has(variable_name))
    {
        ReadElementalVectorialVariableData(rThisElements, static_cast<Variable<Quaternion<double> > const& >(KratosComponents<Variable<Quaternion<double> > >::Get(variable_name)), Vector(4));
    }
    else if(KratosComponents<Variable<Matrix> >::Has(variable_name))
    {
        ReadElementalVectorialVariableData(rThisElements, static_cast<Variable<Matrix > const& >(KratosComponents<Variable<Matrix> >::Get(variable_name)), Matrix(3,3));
    }
    else if(KratosComponents<Variable<Vector> >::Has(variable_name))
    {
        ReadElementalVectorialVariableData(rThisElements, static_cast<Variable<Vector > const& >(KratosComponents<Variable<Vector> >::Get(variable_name)), Vector(3));
    }
    else
    {
        std::stringstream buffer;
        buffer << variable_name << " is not a valid variable!!!" << std::endl;
        buffer << " [Line " << mNumberOfLines << " ]";
        KRATOS_ERROR << buffer.str() << std::endl;
    }

    KRATOS_CATCH("")
}

template<class TVariableType>
void ModelPartIO::ReadElementalScalarVariableData(ElementsContainerType& rThisElements, TVariableType& rVariable)
{
    KRATOS_TRY

    SizeType id;
    double elemental_value;

    std::string value;

    while(!mpStream->eof())
    {
        ReadWord(value); // reading id
        if(CheckEndBlock("ElementalData", value))
            break;

        ExtractValue(value, id);

        // reading nodal_value
        ReadWord(value);
        ExtractValue(value, elemental_value);

        ModelPart::ElementIterator i_result = rThisElements.find(ReorderedElementId(id));
        if(i_result != rThisElements.end())
            i_result->GetValue(rVariable) =  elemental_value;
        else
            KRATOS_WARNING("ModelPartIO")  << "WARNING! Assigning " << rVariable.Name() << " to not existing element #" << id << " [Line " << mNumberOfLines << " ]" << std::endl;
    }

    KRATOS_CATCH("")
}


template<class TVariableType, class TDataType>
void ModelPartIO::ReadElementalVectorialVariableData(ElementsContainerType& rThisElements, TVariableType& rVariable, TDataType Dummy)
{
    KRATOS_TRY

    SizeType id;
    TDataType elemental_value;

    std::string value;

    while(!mpStream->eof())
    {
        ReadWord(value); // reading id
        if(CheckEndBlock("ElementalData", value))
            break;

        ExtractValue(value, id);

        // reading nodal_value
        ReadVectorialValue(elemental_value);

        ModelPart::ElementIterator i_result = rThisElements.find(ReorderedElementId(id));
        if(i_result != rThisElements.end())
            i_result->GetValue(rVariable) =  elemental_value;
        else
            KRATOS_WARNING("ModelPartIO")  << "WARNING! Assigning " << rVariable.Name() << " to not existing element #" << id << " [Line " << mNumberOfLines << " ]" << std::endl;
    }

    KRATOS_CATCH("")
}

void ModelPartIO::ReadConditionalDataBlock(ConditionsContainerType& rThisConditions)
{
    KRATOS_TRY

    typedef VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > array_1d_component_type;

    std::string variable_name;

    ReadWord(variable_name);

    if(KratosComponents<Variable<double> >::Has(variable_name))
    {
        ReadConditionalScalarVariableData(rThisConditions, static_cast<Variable<double> const& >(KratosComponents<Variable<double> >::Get(variable_name)));
    }
    else if(KratosComponents<Variable<bool> >::Has(variable_name))
    {
        ReadConditionalScalarVariableData(rThisConditions, static_cast<Variable<bool> const& >(KratosComponents<Variable<bool> >::Get(variable_name)));
    }
    else if(KratosComponents<Variable<int> >::Has(variable_name))
    {
        ReadConditionalScalarVariableData(rThisConditions, static_cast<Variable<int> const& >(KratosComponents<Variable<int> >::Get(variable_name)));
    }
    else if(KratosComponents<array_1d_component_type>::Has(variable_name))
    {
        ReadConditionalScalarVariableData(rThisConditions, static_cast<array_1d_component_type const& >(KratosComponents<array_1d_component_type>::Get(variable_name)));
    }
    else if(KratosComponents<Variable<array_1d<double, 3> > >::Has(variable_name))
    {
        ReadConditionalVectorialVariableData(rThisConditions, static_cast<Variable<array_1d<double, 3> > const& >(KratosComponents<Variable<array_1d<double, 3> > >::Get(variable_name)), Vector(3));
    }
    else if(KratosComponents<Variable<Quaternion<double> > >::Has(variable_name))
    {
        ReadConditionalVectorialVariableData(rThisConditions, static_cast<Variable<Quaternion<double> > const& >(KratosComponents<Variable<Quaternion<double> > >::Get(variable_name)), Vector(4));
    }
    else if(KratosComponents<Variable<Matrix> >::Has(variable_name))
    {
        ReadConditionalVectorialVariableData(rThisConditions, static_cast<Variable<Matrix > const& >(KratosComponents<Variable<Matrix> >::Get(variable_name)), Matrix(3,3));
    }
    else if(KratosComponents<Variable<Vector> >::Has(variable_name))
    {
        ReadConditionalVectorialVariableData(rThisConditions, static_cast<Variable<Vector > const& >(KratosComponents<Variable<Vector> >::Get(variable_name)), Vector(3));
    }
    else
    {
        std::stringstream buffer;
        buffer << variable_name << " is not a valid variable!!!" << std::endl;
        buffer << " [Line " << mNumberOfLines << " ]";
        KRATOS_ERROR << buffer.str() << std::endl;
    }

    KRATOS_CATCH("")
}

template<class TVariableType>
void ModelPartIO::ReadConditionalScalarVariableData(ConditionsContainerType& rThisConditions, TVariableType& rVariable)
{
    KRATOS_TRY

    SizeType id;
    double conditional_value;

    std::string value;

    while(!mpStream->eof())
    {
        ReadWord(value); // reading id
        if(CheckEndBlock("ConditionalData", value))
            break;

        ExtractValue(value, id);

        // reading nodal_value
        ReadWord(value);
        ExtractValue(value, conditional_value);

        ModelPart::ConditionIterator i_result = rThisConditions.find(ReorderedConditionId(id));
        if(i_result != rThisConditions.end())
            i_result->GetValue(rVariable) =  conditional_value;
        else
            KRATOS_WARNING("ModelPartIO")  << "WARNING! Assigning " << rVariable.Name() << " to not existing condition #" << id << " [Line " << mNumberOfLines << " ]" << std::endl;
    }

    KRATOS_CATCH("")
}


template<class TVariableType, class TDataType>
void ModelPartIO::ReadConditionalVectorialVariableData(ConditionsContainerType& rThisConditions, TVariableType& rVariable, TDataType Dummy)
{
    KRATOS_TRY

    SizeType id;
    TDataType conditional_value;

    std::string value;

    while(!mpStream->eof())
    {
        ReadWord(value); // reading id
        if(CheckEndBlock("ConditionalData", value))
            break;

        ExtractValue(value, id);

        // reading nodal_value
        ReadVectorialValue(conditional_value);

        ModelPart::ConditionIterator i_result = rThisConditions.find(ReorderedConditionId(id));
        if(i_result != rThisConditions.end())
            i_result->GetValue(rVariable) =  conditional_value;
        else
            KRATOS_WARNING("ModelPartIO")  << "WARNING! Assigning " << rVariable.Name() << " to not existing condition #" << id << " [Line " << mNumberOfLines << " ]" << std::endl;
    }

    KRATOS_CATCH("")
}


ModelPartIO::SizeType ModelPartIO::ReadElementsConnectivitiesBlock(ConnectivitiesContainerType& rThisConnectivities)
{
    KRATOS_TRY

    SizeType id;
    SizeType node_id;
    SizeType number_of_connectivities = 0;

    std::string word;
    std::string element_name;

    ReadWord(element_name);
    if(!KratosComponents<Element>::Has(element_name))
    {
        std::stringstream buffer;
        buffer << "Element " << element_name << " is not registered in Kratos.";
        buffer << " Please check the spelling of the element name and see if the application containing it is registered corectly.";
        buffer << " [Line " << mNumberOfLines << " ]";
        KRATOS_ERROR << buffer.str() << std::endl;
        return number_of_connectivities;
    }

    Element const& r_clone_element = KratosComponents<Element>::Get(element_name);
    SizeType number_of_nodes = r_clone_element.GetGeometry().size();
    ConnectivitiesContainerType::value_type temp_element_nodes;

    while(!mpStream->eof())
    {
        ReadWord(word); // Reading the element id or End
        if(CheckEndBlock("Elements", word))
            break;

        ExtractValue(word,id);
        ReadWord(word); // Reading the properties id;
        temp_element_nodes.clear();
        for(SizeType i = 0 ; i < number_of_nodes ; i++)
        {
            ReadWord(word); // Reading the node id;
            ExtractValue(word, node_id);
            temp_element_nodes.push_back(ReorderedNodeId(node_id));
        }
        const int index = ReorderedElementId(id) - 1;
        const int size = rThisConnectivities.size();
        if(index == size)  // I do push back instead of resizing to size+1
            rThisConnectivities.push_back(temp_element_nodes);
        else if(index < size)
            rThisConnectivities[index]= temp_element_nodes;
        else
        {
            rThisConnectivities.resize(index+1);
            rThisConnectivities[index]= temp_element_nodes;

        }
        number_of_connectivities++;
    }

    return number_of_connectivities;

    KRATOS_CATCH("")
}


ModelPartIO::SizeType ModelPartIO::ReadConditionsConnectivitiesBlock(ConnectivitiesContainerType& rThisConnectivities)
{
    KRATOS_TRY

    SizeType id;
    SizeType node_id;
    SizeType number_of_connectivities = 0;

    std::string word;
    std::string condition_name;

    ReadWord(condition_name);
    if(!KratosComponents<Condition>::Has(condition_name))
    {
        std::stringstream buffer;
        buffer << "Condition " << condition_name << " is not registered in Kratos.";
        buffer << " Please check the spelling of the condition name and see if the application containing it is registered corectly.";
        buffer << " [Line " << mNumberOfLines << " ]";
        KRATOS_ERROR << buffer.str() << std::endl;
        return number_of_connectivities;
    }

    Condition const& r_clone_condition = KratosComponents<Condition>::Get(condition_name);
    SizeType number_of_nodes = r_clone_condition.GetGeometry().size();
    ConnectivitiesContainerType::value_type temp_condition_nodes;

    while(!mpStream->eof())
    {
        ReadWord(word); // Reading the condition id or End
        if(CheckEndBlock("Conditions", word))
            break;

        ExtractValue(word,id);
        ReadWord(word); // Reading the properties id;
        temp_condition_nodes.clear();
        for(SizeType i = 0 ; i < number_of_nodes ; i++)
        {
            ReadWord(word); // Reading the node id;
            ExtractValue(word, node_id);
            temp_condition_nodes.push_back(ReorderedNodeId(node_id));
        }

        const int index = ReorderedConditionId(id) - 1;
        const int size = rThisConnectivities.size();
        if(index == size)  // I do push back instead of resizing to size+1
            rThisConnectivities.push_back(temp_condition_nodes);
        else if(index < size)
            rThisConnectivities[index]= temp_condition_nodes;
        else
        {
            rThisConnectivities.resize(index+1);
            rThisConnectivities[index]= temp_condition_nodes;

        }
        number_of_connectivities++;
    }

    return number_of_connectivities;

    KRATOS_CATCH("")
}

void ModelPartIO::FillNodalConnectivitiesFromElementBlock(ConnectivitiesContainerType& rNodalConnectivities)
{
    KRATOS_TRY;

    SizeType id;
    SizeType node_id;
    SizeType position;
    SizeType used_size = rNodalConnectivities.size();
    SizeType reserved_size = (rNodalConnectivities.capacity() > 0) ? rNodalConnectivities.capacity() : 1;

    std::string word;
    std::string element_name;

    ReadWord(element_name);
    if(!KratosComponents<Element>::Has(element_name))
    {
        std::stringstream buffer;
        buffer << "Element " << element_name << " is not registered in Kratos.";
        buffer << " Please check the spelling of the element name and see if the application containing it is registered corectly.";
        buffer << " [Line " << mNumberOfLines << " ]";
        KRATOS_ERROR << buffer.str() << std::endl;
    }

    Element const& r_clone_element = KratosComponents<Element>::Get(element_name);
    SizeType n_nodes_in_elem = r_clone_element.GetGeometry().size();
    ConnectivitiesContainerType::value_type temp_element_nodes;

    while(!mpStream->eof())
    {
        ReadWord(word); // Reading the element id or End
        if(CheckEndBlock("Elements", word))
            break;

        ExtractValue(word,id);
        ReadWord(word); // Reading the properties id;
        temp_element_nodes.clear();
        for(SizeType i = 0 ; i < n_nodes_in_elem ; i++)
        {
            ReadWord(word); // Reading the node id;
            ExtractValue(word, node_id);
            temp_element_nodes.push_back(ReorderedNodeId(node_id));
        }

        for (SizeType i = 0; i < n_nodes_in_elem; i++)
        {
            position = temp_element_nodes[i]-1; // Ids start from 1, position in rNodalConnectivities starts from 0
            if (position >= used_size)
            {
                used_size = position+1;
                if (position >= reserved_size)
                {
                    reserved_size = (used_size > reserved_size) ? 2*used_size : 2*reserved_size;
                    rNodalConnectivities.reserve(reserved_size);
                }
                rNodalConnectivities.resize(used_size);
            }

            for (SizeType j = 0; j < i; j++)
                rNodalConnectivities[position].push_back(temp_element_nodes[j]);
            for (SizeType j = i+1; j < n_nodes_in_elem; j++)
                rNodalConnectivities[position].push_back(temp_element_nodes[j]);
        }
    }

    KRATOS_CATCH("");
}

void ModelPartIO::FillNodalConnectivitiesFromConditionBlock(ConnectivitiesContainerType& rNodalConnectivities)
{
    KRATOS_TRY;

    SizeType id;
    SizeType node_id;
    SizeType position;
    SizeType used_size = rNodalConnectivities.size();
    SizeType reserved_size = (rNodalConnectivities.capacity() > 0) ? rNodalConnectivities.capacity() : 1;

    std::string word;
    std::string condition_name;

    ReadWord(condition_name);
    if(!KratosComponents<Condition>::Has(condition_name))
    {
        std::stringstream buffer;
        buffer << "Condition " << condition_name << " is not registered in Kratos.";
        buffer << " Please check the spelling of the condition name and see if the application containing it is registered corectly.";
        buffer << " [Line " << mNumberOfLines << " ]";
        KRATOS_ERROR << buffer.str() << std::endl;;
    }

    Condition const& r_clone_condition = KratosComponents<Condition>::Get(condition_name);
    SizeType n_nodes_in_cond = r_clone_condition.GetGeometry().size();
    ConnectivitiesContainerType::value_type temp_condition_nodes;

    while(!mpStream->eof())
    {
        ReadWord(word); // Reading the condition id or End
        if(CheckEndBlock("Conditions", word))
            break;

        ExtractValue(word,id);
        ReadWord(word); // Reading the properties id;
        temp_condition_nodes.clear();
        for(SizeType i = 0 ; i < n_nodes_in_cond ; i++)
        {
            ReadWord(word); // Reading the node id;
            ExtractValue(word, node_id);
            temp_condition_nodes.push_back(ReorderedNodeId(node_id));
        }

        for (SizeType i = 0; i < n_nodes_in_cond; i++)
        {
            position = temp_condition_nodes[i]-1; // Ids start from 1, position in rNodalConnectivities starts from 0
            if (position >= used_size)
            {
                used_size = position+1;
                if (position >= reserved_size)
                {
                    reserved_size = (used_size > reserved_size) ? 2*used_size : 2*reserved_size;
                    rNodalConnectivities.reserve(reserved_size);
                }
                rNodalConnectivities.resize(used_size);
            }

            for (SizeType j = 0; j < i; j++)
                rNodalConnectivities[position].push_back(temp_condition_nodes[j]);
            for (SizeType j = i+1; j < n_nodes_in_cond; j++)
                rNodalConnectivities[position].push_back(temp_condition_nodes[j]);
        }
    }

    KRATOS_CATCH("");
}


void ModelPartIO::ReadCommunicatorDataBlock(Communicator& rThisCommunicator, NodesContainerType& rThisNodes)
{
    KRATOS_TRY

    std::string word;
    while(true)
    {
        ReadWord(word);
        if(mpStream->eof())
            break;
        if(CheckEndBlock("CommunicatorData", word))
            break;
        if(word == "NEIGHBOURS_INDICES")
        {
            ReadVectorialValue(rThisCommunicator.NeighbourIndices());
        }
        else if(word == "NUMBER_OF_COLORS")
        {
            ReadWord(word);
            SizeType number_of_colors;
            ExtractValue(word, number_of_colors);
            rThisCommunicator.SetNumberOfColors(number_of_colors);
        }
        else
        {
            ReadBlockName(word);
            if(word == "LocalNodes")
            {
                ReadCommunicatorLocalNodesBlock(rThisCommunicator, rThisNodes);
            }
            else if(word == "GhostNodes")
            {
                ReadCommunicatorGhostNodesBlock(rThisCommunicator, rThisNodes);
            }
            else
            {
                SkipBlock(word);
            }
        }
    }

    return ;

    KRATOS_CATCH("")
}

void ModelPartIO::ReadCommunicatorLocalNodesBlock(Communicator& rThisCommunicator, NodesContainerType& rThisNodes)
{
    KRATOS_TRY

    SizeType interface_id;
    SizeType node_id;

    std::string word;
    std::string condition_name;

    ReadWord(word); // reading the interface id
    ExtractValue(word,interface_id);

    if(interface_id > rThisCommunicator.GetNumberOfColors())
    {
        std::stringstream buffer;
        buffer << "Interface " << interface_id << " is not valid.";
        buffer << " The number of colors is " << rThisCommunicator.GetNumberOfColors() << " and the interface id must be les than or equal to number of colors" ;
        buffer << " [Line " << mNumberOfLines << " ]";
        KRATOS_ERROR << buffer.str() << std::endl;;
    }

    Communicator::MeshType* p_local_mesh;
    Communicator::MeshType* p_interface_mesh;

    if(interface_id == 0)
    {
        p_local_mesh = &(rThisCommunicator.LocalMesh());
        p_interface_mesh = &(rThisCommunicator.InterfaceMesh());
    }
    else
    {
        p_local_mesh = &(rThisCommunicator.LocalMesh(interface_id-1));
        p_interface_mesh = &(rThisCommunicator.InterfaceMesh(interface_id-1));
    }

    NodesContainerType aux_local;
    NodesContainerType aux_interface;

    while(!mpStream->eof())
    {
        ReadWord(word); // Reading the node id or End
        if(CheckEndBlock("LocalNodes", word))
            break;

        ExtractValue(word,node_id);
        NodesContainerType::iterator it_node = FindKey(rThisNodes, ReorderedNodeId(node_id), "Node");
        auto p_node = *(it_node.base());
        aux_local.push_back(p_node);
        aux_interface.push_back(p_node);
    }

    for(auto it = aux_local.begin(); it!= aux_local.end(); it++)
        p_local_mesh->Nodes().push_back(*(it.base()));

    for(auto it = aux_interface.begin(); it!= aux_interface.end(); it++)
        p_interface_mesh->Nodes().push_back(*(it.base()));

    p_local_mesh->Nodes().Unique();
    p_interface_mesh->Nodes().Unique();

    KRATOS_CATCH("")
}


void ModelPartIO::ReadCommunicatorGhostNodesBlock(Communicator& rThisCommunicator, NodesContainerType& rThisNodes)
{
    KRATOS_TRY

    SizeType interface_id;
    SizeType node_id;


    std::string word;
    std::string condition_name;

    ReadWord(word); // reading the interface id
    ExtractValue(word,interface_id);


    if(interface_id > rThisCommunicator.GetNumberOfColors())
    {
        std::stringstream buffer;
        buffer << "Interface " << interface_id << " is not valid.";
        buffer << " The number of colors is " << rThisCommunicator.GetNumberOfColors() << " and the interface id must be les than or equal to number of colors" ;
        buffer << " [Line " << mNumberOfLines << " ]";
        KRATOS_ERROR << buffer.str() << std::endl;;
    }

    Communicator::MeshType* p_ghost_mesh;
    Communicator::MeshType* p_interface_mesh;

    if(interface_id == 0)
    {
        p_ghost_mesh = &(rThisCommunicator.GhostMesh());
        p_interface_mesh = &(rThisCommunicator.InterfaceMesh());
    }
    else
    {
        p_ghost_mesh = &(rThisCommunicator.GhostMesh(interface_id-1));
        p_interface_mesh = &(rThisCommunicator.InterfaceMesh(interface_id-1));
    }

    while(!mpStream->eof())
    {
        ReadWord(word); // Reading the node id or End
        if(CheckEndBlock("GhostNodes", word))
            break;

        ExtractValue(word,node_id);
        NodesContainerType::iterator it_node = FindKey(rThisNodes, ReorderedNodeId(node_id), "Node");
        p_ghost_mesh->Nodes().push_back(*(it_node.base()));
        p_interface_mesh->Nodes().push_back(*(it_node.base()));
    }

    p_ghost_mesh->Nodes().Unique();
    p_interface_mesh->Nodes().Unique();

    KRATOS_CATCH("")
}

void ModelPartIO::ReadMeshBlock(ModelPart& rModelPart)
{
    KRATOS_TRY

    std::string word;
    SizeType mesh_id;

    ReadWord(word);
    ExtractValue(word, mesh_id);

    SizeType number_of_meshes = rModelPart.NumberOfMeshes();

    // This would be a case of error in reading.
    KRATOS_ERROR_IF(mesh_id > 1000000) << "Too large mesh id : " << mesh_id << std::endl;

    // This would be a case of error in reading.
    KRATOS_ERROR_IF(mesh_id == 0) << "The mesh zero is the reference mesh and already created. You cannot create a mesh 0 with mesh block." << std::endl;

    // adding necessary meshes to the model part.
    MeshType empty_mesh;
    for(SizeType i = number_of_meshes ; i < mesh_id + 1 ; i++)
        rModelPart.GetMeshes().push_back(Kratos::make_shared<MeshType>(empty_mesh.Clone()));

    MeshType& mesh = rModelPart.GetMesh(mesh_id);

    while(true)
    {
        ReadWord(word);

        if(mpStream->eof())
        {
            break;
        }

        if(CheckEndBlock("Mesh", word))
        {
                break;
        }

        ReadBlockName(word);
        if(word == "MeshData")
        {
            ReadMeshDataBlock(mesh);
        }
        else if(word == "MeshNodes")
        {
            ReadMeshNodesBlock(rModelPart, mesh);
        }

        else if(word == "MeshElements")
        {
            ReadMeshElementsBlock(rModelPart, mesh);
        }

        else if(word == "MeshConditions")
        {
            ReadMeshConditionsBlock(rModelPart, mesh);
        }

//             else if(word == "MeshProperties")
//                 ReadMeshPropertiesBlock(rModelPart, mesh);

        else
            {
                SkipBlock(word);
            }
    }

    KRATOS_CATCH("")

}


void ModelPartIO::ReadMeshDataBlock(MeshType& rMesh)
{
    KRATOS_TRY

    std::string variable_name;

    while(!mpStream->eof())
    {
        ReadWord(variable_name);
        if(CheckEndBlock("MeshData", variable_name))
            break;
        if(KratosComponents<Variable<double> >::Has(variable_name))
        {
            std::string value;
            double temp;

            ReadWord(value); // reading value
            ExtractValue(value,temp);
            rMesh[KratosComponents<Variable<double> >::Get(variable_name)] = temp;
        }
        else if(KratosComponents<Variable<bool> >::Has(variable_name))
        {
            std::string value;
            bool temp;

            ReadWord(value); // reading value
            ExtractValue(value,temp);
            rMesh[KratosComponents<Variable<bool> >::Get(variable_name)] = temp;
        }
        else if(KratosComponents<Variable<int> >::Has(variable_name))
        {
            std::string value;
            int temp;

            ReadWord(value); // reading value
            ExtractValue(value,temp);
            rMesh[KratosComponents<Variable<int> >::Get(variable_name)] = temp;
        }
        else if(KratosComponents<Variable<array_1d<double, 3> > >::Has(variable_name))
        {
            Vector temp_vector; // defining a Vector because for array_1d the operator >> is not defined yet!
            ReadVectorialValue(temp_vector);
            rMesh[KratosComponents<Variable<array_1d<double,3> > >::Get(variable_name)] = temp_vector;
        }
        else if(KratosComponents<Variable<Quaternion<double> > >::Has(variable_name))
        {
            Vector temp_vector; // defining a Vector because for Quaternion the operator >> is not defined yet!
            ReadVectorialValue(temp_vector);
            rMesh[KratosComponents<Variable<Quaternion<double> > >::Get(variable_name)] = temp_vector;
        }
        else if(KratosComponents<Variable<Matrix> >::Has(variable_name))
        {
            ReadVectorialValue(rMesh[KratosComponents<Variable<Matrix> >::Get(variable_name)]);
        }
        else if(KratosComponents<Variable<std::string> >::Has(variable_name))
        {
            std::string value;
    std::string  temp;

            ReadWord(value); // reading value
            ExtractValue(value,temp);
            rMesh[KratosComponents<Variable<std::string> >::Get(variable_name)] = temp;
        }
        else
        {
            std::stringstream buffer;
            buffer << variable_name << " is not a valid variable!!!" << std::endl;
            buffer << " [Line " << mNumberOfLines << " ]";
            KRATOS_ERROR << buffer.str() << std::endl;;
        }


    }

    KRATOS_CATCH("")
}


void ModelPartIO::ReadMeshNodesBlock(ModelPart& rModelPart, MeshType& rMesh)
{
    KRATOS_TRY

    SizeType node_id;

    std::string word;


    while(!mpStream->eof())
    {
        ReadWord(word); // Reading the node id or End
        if(CheckEndBlock("MeshNodes", word))
            break;

        ExtractValue(word,node_id);
        NodesContainerType::iterator it_node = FindKey(rModelPart.Nodes(), ReorderedNodeId(node_id), "Node");
        rMesh.Nodes().push_back(*(it_node.base()));
    }

    rMesh.Nodes().Sort();
    KRATOS_CATCH("")
}

void ModelPartIO::ReadMeshElementsBlock(ModelPart& rModelPart, MeshType& rMesh)
{
    KRATOS_TRY

    SizeType element_id;

    std::string word;


    while(!mpStream->eof())
    {
        ReadWord(word); // Reading the element id or End
        if(CheckEndBlock("MeshElements", word))
            break;

        ExtractValue(word,element_id);
        ElementsContainerType::iterator i_element = FindKey(rModelPart.Elements(), ReorderedElementId(element_id), "Element");
        rMesh.Elements().push_back(*(i_element.base()));
    }

    rMesh.Elements().Sort();
    KRATOS_CATCH("")
}

void ModelPartIO::ReadMeshConditionsBlock(ModelPart& rModelPart, MeshType& rMesh)
{
    KRATOS_TRY

    SizeType condition_id;

    std::string word;


    while(!mpStream->eof())
    {
        ReadWord(word); // Reading the element id or End
        if(CheckEndBlock("MeshConditions", word))
            break;

        ExtractValue(word,condition_id);
        ConditionsContainerType::iterator i_condition = FindKey(rModelPart.Conditions(), ReorderedConditionId(condition_id), "Condition");
        rMesh.Conditions().push_back(*(i_condition.base()));
    }

    rMesh.Conditions().Sort();
    KRATOS_CATCH("")
}

void ModelPartIO::ReadMeshPropertiesBlock(ModelPart& rModelPart, MeshType& rMesh)
{
    KRATOS_TRY

    Properties::Pointer props = Kratos::make_shared<Properties>();
    Properties& temp_properties = *props;
//         Properties temp_properties;

    std::string word;
    std::string variable_name;

    SizeType temp_properties_id;

    ReadWord(word);
    ExtractValue(word, temp_properties_id);
    temp_properties.SetId(temp_properties_id);


    while(!mpStream->eof())
    {
        ReadWord(variable_name);
        if(CheckEndBlock("MeshProperties", variable_name))
            break;

    if(KratosComponents<Variable<std::string> >::Has(variable_name))
        {
            std::string value;
            std::string  temp;

            ReadWord(value); // reading value
            ExtractValue(value,temp);
            temp_properties[KratosComponents<Variable<std::string> >::Get(variable_name)] = temp;
        }
    else if(KratosComponents<Variable<double> >::Has(variable_name))
        {
            std::string value;
            double temp;

            ReadWord(value); // reading value
            ExtractValue(value,temp);
            temp_properties[KratosComponents<Variable<double> >::Get(variable_name)] = temp;
        }
        else if(KratosComponents<Variable<int> >::Has(variable_name))
        {
            std::string value;
            int temp;

            ReadWord(value); // reading value
            ExtractValue(value,temp);
            temp_properties[KratosComponents<Variable<int> >::Get(variable_name)] = temp;
        }
        else if(KratosComponents<Variable<bool> >::Has(variable_name))
        {
            std::string value;
            bool temp;

            ReadWord(value); // reading value
            ExtractValue(value,temp);
            temp_properties[KratosComponents<Variable<bool> >::Get(variable_name)] = temp;
        }
        else if(KratosComponents<Variable<array_1d<double, 3> > >::Has(variable_name))
        {
            Vector temp_vector; // defining a Vector because for array_1d the operator >> is not defined yet!
            ReadVectorialValue(temp_vector);
            temp_properties[KratosComponents<Variable<array_1d<double,3> > >::Get(variable_name)] = temp_vector;
        }
        else if(KratosComponents<Variable<Quaternion<double> > >::Has(variable_name))
        {
            Vector temp_vector; // defining a Vector because for Quaternion the operator >> is not defined yet!
            ReadVectorialValue(temp_vector);
            temp_properties[KratosComponents<Variable<Quaternion<double> > >::Get(variable_name)] = temp_vector;
        }
        else if(KratosComponents<Variable<Vector> >::Has(variable_name))
        {
            ReadVectorialValue(temp_properties[KratosComponents<Variable<Vector> >::Get(variable_name)]);
        }
        else if(KratosComponents<Variable<Matrix> >::Has(variable_name))
        {
            ReadVectorialValue(temp_properties[KratosComponents<Variable<Matrix> >::Get(variable_name)]);
        }
        else
        {
            std::stringstream buffer;
            buffer << variable_name << " is not a valid variable!!!" << std::endl;
            buffer << " [Line " << mNumberOfLines << " ]";
            KRATOS_ERROR << buffer.str() << std::endl;;
        }

    }

    rMesh.Properties().push_back(props);
//         rMesh.Properties().push_back(temp_properties);

    KRATOS_CATCH("")
}

void ModelPartIO::ReadSubModelPartBlock(ModelPart& rMainModelPart, ModelPart& rParentModelPart)
{
    KRATOS_TRY

    std::string word;

    ReadWord(word); // Reading the name of the sub model part

    ModelPart& r_sub_model_part = rParentModelPart.CreateSubModelPart(word);

    while (true)
    {
        ReadWord(word);
        //if (mpStream->eof())
        //    break;
        if (CheckEndBlock("SubModelPart", word))
            break;

        ReadBlockName(word);
        if (word == "SubModelPartData")
            ReadSubModelPartDataBlock(r_sub_model_part);
        else if (word == "SubModelPartTables")
            ReadSubModelPartTablesBlock(rMainModelPart, r_sub_model_part);
        else if (word == "SubModelPartProperties")
            ReadSubModelPartPropertiesBlock(rMainModelPart, r_sub_model_part);
        else if (word == "SubModelPartNodes")
            ReadSubModelPartNodesBlock(rMainModelPart, r_sub_model_part);
        else if (word == "SubModelPartElements")
            ReadSubModelPartElementsBlock(rMainModelPart, r_sub_model_part);
        else if (word == "SubModelPartConditions")
            ReadSubModelPartConditionsBlock(rMainModelPart, r_sub_model_part);
        // TODO: Add the following blocks. Pooyan.
        //else if (word == "CommunicatorData")
        //{
        //    ReadCommunicatorDataBlock(rThisModelPart.GetCommunicator(), rThisModelPart.Nodes());
        //    //Adding the elements and conditions to the communicator
        //    rThisModelPart.GetCommunicator().LocalMesh().Elements() = rThisModelPart.Elements();
        //    rThisModelPart.GetCommunicator().LocalMesh().Conditions() = rThisModelPart.Conditions();
        //}
        //else if (word == "Mesh")
        //    ReadMeshBlock(rThisModelPart);
        else if (word == "SubModelPart")
            ReadSubModelPartBlock(rMainModelPart, r_sub_model_part);
    }

    KRATOS_CATCH("")
}

void ModelPartIO::WriteSubModelPartBlock(
    ModelPart& rMainModelPart,
    const std::string& InitialTabulation) {

    KRATOS_TRY;

    const std::vector<std::string> sub_model_part_names = rMainModelPart.GetSubModelPartNames();

    for (unsigned int i_sub = 0; i_sub < sub_model_part_names.size(); i_sub++) {

        const std::string sub_model_part_name = sub_model_part_names[i_sub];
        ModelPart& r_sub_model_part = rMainModelPart.GetSubModelPart(sub_model_part_name);

        (*mpStream) << InitialTabulation << "Begin SubModelPart\t" << sub_model_part_name << std::endl;

        // Submodelpart data section
        (*mpStream) << InitialTabulation << "\tBegin SubModelPartData" << std::endl;
        // VARIABLE_NAME value // TODO: Finish me
        (*mpStream) << InitialTabulation  << "\tEnd SubModelPartData" << std::endl;

        // Submodelpart tables section
        (*mpStream) << InitialTabulation  << "\tBegin SubModelPartTables" << std::endl;
//                     ModelPart::TablesContainerType& rThisTables = rMainModelPart.Tables();
//                     auto numTables = rThisTables.end() - rThisTables.begin();
//                     for(unsigned int i = 0; i < numTables; i++)
//                     {
//                         auto itTable = rThisTables.begin() + i;
//                         (*mpStream) << InitialTabulation << "\t" << itTable->Id() << std::endl; //FIXME: Tables does not have Id() Whyyyyy?
//                     }
        (*mpStream) << InitialTabulation << "\tEnd SubModelPartTables" << std::endl;

        // Submodelpart nodes section
        (*mpStream) << InitialTabulation << "\tBegin SubModelPartNodes" << std::endl;
        NodesContainerType& rThisNodes = r_sub_model_part.Nodes();
        auto numNodes = rThisNodes.end() - rThisNodes.begin();
        for(unsigned int i = 0; i < numNodes; i++) {
            auto itNode = rThisNodes.begin() + i;
            (*mpStream) << InitialTabulation << "\t\t" << itNode->Id() << std::endl;
        }
        (*mpStream) << InitialTabulation << "\tEnd SubModelPartNodes" << std::endl;

        // Submodelpart elements section
        (*mpStream) << InitialTabulation << "\tBegin SubModelPartElements" << std::endl;
        ElementsContainerType& rThisElements = r_sub_model_part.Elements();
        auto num_elements = rThisElements.end() - rThisElements.begin();
        for(unsigned int i = 0; i < num_elements; i++) {
            auto itElem = rThisElements.begin() + i;
            (*mpStream) << InitialTabulation << "\t\t" << itElem->Id() << std::endl;
        }
        (*mpStream) << InitialTabulation << "\tEnd SubModelPartElements" << std::endl;

        // Submodelpart conditions section
        (*mpStream) << InitialTabulation << "\tBegin SubModelPartConditions" << std::endl;
        ConditionsContainerType& rThisConditions= r_sub_model_part.Conditions();
        auto numConditions = rThisConditions.end() - rThisConditions.begin();
        for(unsigned int i = 0; i < numConditions; i++) {
            auto itCond = rThisConditions.begin() + i;
            (*mpStream) << InitialTabulation << "\t\t" << itCond->Id() << std::endl;
        }
        (*mpStream) << InitialTabulation << "\tEnd SubModelPartConditions" << std::endl;

        // Write the subsubmodelparts
        WriteSubModelPartBlock(r_sub_model_part, InitialTabulation+"\t");

        (*mpStream) << InitialTabulation << "End SubModelPart\t" << std::endl << std::endl;
    }

    KRATOS_CATCH("");
}

void  ModelPartIO::ReadSubModelPartTablesBlock(ModelPart& rMainModelPart, ModelPart& rSubModelPart)
{
    KRATOS_TRY

    SizeType table_id;
    std::string word;

    while (!mpStream->eof())
    {
        ReadWord(word); // Reading the node id or End
        if (CheckEndBlock("SubModelPartTables", word))
            break;

        ExtractValue(word, table_id);
        ModelPart::TablesContainerType::iterator i_table = FindKey(rMainModelPart.Tables(), table_id, "Table");
        rSubModelPart.AddTable((i_table.base())->first, (i_table.base())->second);
    }
    KRATOS_CATCH("")
}

void  ModelPartIO::ReadSubModelPartPropertiesBlock(ModelPart& rMainModelPart, ModelPart& rSubModelPart)
{
    KRATOS_TRY

    SizeType properties_id;
    std::string word;

    while (!mpStream->eof())
    {
        ReadWord(word); // Reading the node id or End
        if (CheckEndBlock("SubModelPartProperties", word))
            break;

        ExtractValue(word, properties_id);
        PropertiesContainerType::iterator i_properties = FindKey(rMainModelPart.rProperties(), properties_id, "Properties");
        rSubModelPart.AddProperties(*(i_properties.base()));
    }
    KRATOS_CATCH("")
}

void  ModelPartIO::ReadSubModelPartNodesBlock(ModelPart& rMainModelPart, ModelPart& rSubModelPart)
{
    KRATOS_TRY

    SizeType node_id;
    std::string word;

            std::vector<SizeType> ordered_ids;

    while (!mpStream->eof())
    {
        ReadWord(word); // Reading the node id or End
        if (CheckEndBlock("SubModelPartNodes", word))
            break;

        ExtractValue(word, node_id);
                    ordered_ids.push_back(ReorderedNodeId(node_id));
    }

    std::sort(ordered_ids.begin(), ordered_ids.end());
            rSubModelPart.AddNodes(ordered_ids);

    KRATOS_CATCH("")
}

void  ModelPartIO::ReadSubModelPartElementsBlock(ModelPart& rMainModelPart, ModelPart& rSubModelPart)
{
    KRATOS_TRY

    SizeType element_id;
    std::string word;
            std::vector<SizeType> ordered_ids;

    while (!mpStream->eof())
    {
        ReadWord(word); // Reading the node id or End
        if (CheckEndBlock("SubModelPartElements", word))
            break;

        ExtractValue(word, element_id);
                    ordered_ids.push_back(ReorderedElementId(element_id));
    }
    std::sort(ordered_ids.begin(), ordered_ids.end());
    rSubModelPart.AddElements(ordered_ids);

    KRATOS_CATCH("")
}

void  ModelPartIO::ReadSubModelPartConditionsBlock(ModelPart& rMainModelPart, ModelPart& rSubModelPart)
{
    KRATOS_TRY

    SizeType condition_id;
    std::string word;
            std::vector<SizeType> ordered_ids;

    while (!mpStream->eof())
    {
        ReadWord(word); // Reading the node id or End
        if (CheckEndBlock("SubModelPartConditions", word))
            break;

        ExtractValue(word, condition_id);
                    ordered_ids.push_back(ReorderedConditionId(condition_id));
    }
    std::sort(ordered_ids.begin(), ordered_ids.end());
    rSubModelPart.AddConditions(ordered_ids);

    KRATOS_CATCH("")
}


void ModelPartIO::DivideModelPartDataBlock(OutputFilesContainerType& OutputFiles)
{
    KRATOS_TRY
    std::string block;

    WriteInAllFiles(OutputFiles, "Begin ModelPartData\n");

    ReadBlock(block, "ModelPartData");
    WriteInAllFiles(OutputFiles, block);

    WriteInAllFiles(OutputFiles, "End ModelPartData\n");

    KRATOS_CATCH("")
}

void ModelPartIO::DivideTableBlock(OutputFilesContainerType& OutputFiles)
{
    KRATOS_TRY

    std::string block;

    WriteInAllFiles(OutputFiles, "Begin Table ");

    ReadBlock(block, "Table");
    WriteInAllFiles(OutputFiles, block);

    WriteInAllFiles(OutputFiles, "End Table\n");

    KRATOS_CATCH("")
}

void ModelPartIO::DividePropertiesBlock(OutputFilesContainerType& OutputFiles)
{
    KRATOS_TRY

    std::string block;

    WriteInAllFiles(OutputFiles, "Begin Properties ");

    ReadBlock(block, "Properties");
    WriteInAllFiles(OutputFiles, block);

    WriteInAllFiles(OutputFiles, "End Properties\n");

    KRATOS_CATCH("")
}

void ModelPartIO::DivideNodesBlock(OutputFilesContainerType& OutputFiles,
                        PartitionIndicesContainerType const& NodesAllPartitions)
{
    KRATOS_TRY

    std::string word;

    WriteInAllFiles(OutputFiles, "Begin Nodes \n");

    SizeType id;

    while(!mpStream->eof())
    {
        ReadWord(word);
        if(CheckEndBlock("Nodes", word))
            break;

        ExtractValue(word, id);

        if(ReorderedNodeId(id) > NodesAllPartitions.size())
        {
            std::stringstream buffer;
            buffer << "Invalid node id : " << id;
            buffer << " [Line " << mNumberOfLines << " ]";
            KRATOS_ERROR << buffer.str() << std::endl;;
        }

        std::stringstream node_data;
        node_data << ReorderedNodeId(id) << '\t'; // id
        ReadWord(word);
        node_data << word << '\t'; // x
        ReadWord(word);
        node_data << word << '\t'; // y
        ReadWord(word);
        node_data << word << '\n'; // z

        for(SizeType i = 0 ; i < NodesAllPartitions[ReorderedNodeId(id)-1].size() ; i++)
        {
            SizeType partition_id = NodesAllPartitions[ReorderedNodeId(id)-1][i];
            if(partition_id > OutputFiles.size())
            {
                std::stringstream buffer;
                buffer << "Invalid prtition id : " << partition_id;
                buffer << " for node " << id << " [Line " << mNumberOfLines << " ]";
                KRATOS_ERROR << buffer.str() << std::endl;;
            }

            *(OutputFiles[partition_id]) << node_data.str();
        }

    }

    WriteInAllFiles(OutputFiles, "End Nodes\n");

    KRATOS_CATCH("")
}

void ModelPartIO::DivideElementsBlock(OutputFilesContainerType& OutputFiles,
                            PartitionIndicesContainerType const& ElementsAllPartitions)
{
    KRATOS_TRY


    std::string word;
    std::string element_name;

    ReadWord(element_name);
    if(!KratosComponents<Element>::Has(element_name))
    {
        std::stringstream buffer;
        buffer << "Element " << element_name << " is not registered in Kratos.";
        buffer << " Please check the spelling of the element name and see if the application containing it is registered corectly.";
        buffer << " [Line " << mNumberOfLines << " ]";
        KRATOS_ERROR << buffer.str() << std::endl;;
        return;
    }

    Element const& r_clone_element = KratosComponents<Element>::Get(element_name);
    SizeType number_of_nodes = r_clone_element.GetGeometry().size();

    WriteInAllFiles(OutputFiles, "Begin Elements " +  element_name);

    SizeType id;

    while(!mpStream->eof())
    {
        ReadWord(word); // Reading the element id or End
        if(CheckEndBlock("Elements", word))
            break;

        ExtractValue(word,id);
        if(ReorderedElementId(id) > ElementsAllPartitions.size())
        {
            std::stringstream buffer;
            buffer << "Invalid element id : " << id;
            buffer << " [Line " << mNumberOfLines << " ]";
            KRATOS_ERROR << buffer.str() << std::endl;;
        }

        std::stringstream element_data;
        element_data << '\n' << ReorderedElementId(id) << '\t'; // id
        ReadWord(word); // Reading the properties id;
        element_data << word << '\t'; // properties id

        for(SizeType i = 0 ; i < number_of_nodes ; i++)
        {
            ReadWord(word); // Reading the node id;
            SizeType node_id;
            ExtractValue(word, node_id);
            element_data << ReorderedNodeId(node_id) << '\t'; // node id
        }


        for(SizeType i = 0 ; i < ElementsAllPartitions[ReorderedElementId(id)-1].size() ; i++)
        {
            SizeType partition_id = ElementsAllPartitions[ReorderedElementId(id)-1][i];
            if(partition_id > OutputFiles.size())
            {
                std::stringstream buffer;
                buffer << "Invalid prtition id : " << partition_id;
                buffer << " for node " << id << " [Line " << mNumberOfLines << " ]";
                KRATOS_ERROR << buffer.str() << std::endl;;
            }

            *(OutputFiles[partition_id]) << element_data.str();
        }

    }

    WriteInAllFiles(OutputFiles, "\nEnd Elements\n");

    KRATOS_CATCH("")
}



void ModelPartIO::DivideConditionsBlock(OutputFilesContainerType& OutputFiles,
                            PartitionIndicesContainerType const& ConditionsAllPartitions)
{
    KRATOS_TRY

    std::string word;
    std::string condition_name;

    ReadWord(condition_name);
    if(!KratosComponents<Condition>::Has(condition_name))
    {
        std::stringstream buffer;
        buffer << "Condition " << condition_name << " is not registered in Kratos.";
        buffer << " Please check the spelling of the condition name and see if the application containing it is registered corectly.";
        buffer << " [Line " << mNumberOfLines << " ]";
        KRATOS_ERROR << buffer.str() << std::endl;;
        return;
    }

    Condition const& r_clone_condition = KratosComponents<Condition>::Get(condition_name);
    SizeType number_of_nodes = r_clone_condition.GetGeometry().size();

    WriteInAllFiles(OutputFiles, "Begin Conditions " +  condition_name);

    SizeType id;

    while(!mpStream->eof())
    {
        ReadWord(word); // Reading the condition id or End
        if(CheckEndBlock("Conditions", word))
            break;

        ExtractValue(word,id);
        if(ReorderedConditionId(id) > ConditionsAllPartitions.size())
        {
            std::stringstream buffer;
            buffer << "Invalid condition id : " << id;
            buffer << " [Line " << mNumberOfLines << " ]";
            KRATOS_ERROR << buffer.str() << std::endl;;
        }

        std::stringstream condition_data;
        condition_data << '\n' << ReorderedConditionId(id) << '\t'; // id
        ReadWord(word); // Reading the properties id;
        condition_data << word << '\t'; // properties id

        for(SizeType i = 0 ; i < number_of_nodes ; i++)
        {
            ReadWord(word); // Reading the node id;
            SizeType node_id;
            ExtractValue(word, node_id);
            condition_data << ReorderedNodeId(node_id) << '\t'; // node id
        }


        for(SizeType i = 0 ; i < ConditionsAllPartitions[ReorderedConditionId(id)-1].size() ; i++)
        {
            SizeType partition_id = ConditionsAllPartitions[ReorderedConditionId(id)-1][i];
            if(partition_id > OutputFiles.size())
            {
                std::stringstream buffer;
                buffer << "Invalid prtition id : " << partition_id;
                buffer << " for node " << id << " [Line " << mNumberOfLines << " ]";
                KRATOS_ERROR << buffer.str() << std::endl;;
            }

            *(OutputFiles[partition_id]) << condition_data.str();
        }

    }

    WriteInAllFiles(OutputFiles, "\nEnd Conditions\n");

    KRATOS_CATCH("")
}


void ModelPartIO::DivideNodalDataBlock(OutputFilesContainerType& OutputFiles,
                            PartitionIndicesContainerType const& NodesAllPartitions)
{
    KRATOS_TRY

    typedef VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > array_1d_component_type;

    std::string word;

    WriteInAllFiles(OutputFiles, "Begin NodalData ");

    std::string variable_name;

    ReadWord(variable_name);

    WriteInAllFiles(OutputFiles, variable_name);
    WriteInAllFiles(OutputFiles, "\n");

    if(KratosComponents<Variable<double> >::Has(variable_name))
    {
        DivideDofVariableData(OutputFiles, NodesAllPartitions);
    }
    else if(KratosComponents<Variable<int> >::Has(variable_name))
    {
        DivideDofVariableData(OutputFiles, NodesAllPartitions);
    }
    else if(KratosComponents<Variable<bool> >::Has(variable_name))
    {
        DivideDofVariableData(OutputFiles, NodesAllPartitions);
    }
    else if(KratosComponents<array_1d_component_type>::Has(variable_name))
    {
        DivideDofVariableData(OutputFiles, NodesAllPartitions);
    }
    else if(KratosComponents<Variable<array_1d<double, 3> > >::Has(variable_name))
    {
        DivideVectorialVariableData(OutputFiles, NodesAllPartitions, "NodalData");
    }
    else if(KratosComponents<Variable<Quaternion<double> > >::Has(variable_name))
    {
        DivideVectorialVariableData(OutputFiles, NodesAllPartitions, "NodalData");
    }
    else if(KratosComponents<Variable<Matrix> >::Has(variable_name))
    {
        DivideVectorialVariableData(OutputFiles, NodesAllPartitions, "NodalData" );
    }
    else if(KratosComponents<VariableData>::Has(variable_name))
    {
        std::stringstream buffer;
        buffer << variable_name << " is not supported to be read by this IO or the type of variable is not registered correctly" << std::endl;
        buffer << " [Line " << mNumberOfLines << " ]";
        KRATOS_ERROR << buffer.str() << std::endl;;
    }
    else
    {
        std::stringstream buffer;
        buffer << variable_name << " is not a valid variable!!!" << std::endl;
        buffer << " [Line " << mNumberOfLines << " ]";
        KRATOS_ERROR << buffer.str() << std::endl;;
    }

    WriteInAllFiles(OutputFiles, "End NodalData\n");

    KRATOS_CATCH("")
}

void ModelPartIO::DivideDofVariableData(OutputFilesContainerType& OutputFiles,
                            PartitionIndicesContainerType const& NodesAllPartitions)
{
    KRATOS_TRY

    SizeType id;

    std::string word;


    while(!mpStream->eof())
    {
        ReadWord(word); // reading id
        if(CheckEndBlock("NodalData", word))
            break;

        ExtractValue(word, id);

        if(ReorderedNodeId(id) > NodesAllPartitions.size())
        {
            std::stringstream buffer;
            buffer << "Invalid node id : " << id;
            buffer << " [Line " << mNumberOfLines << " ]";
            KRATOS_ERROR << buffer.str() << std::endl;;
        }

        std::stringstream node_data;
        node_data << ReorderedNodeId(id) << '\t'; // id
        ReadWord(word);
        node_data << word << '\t'; // is fixed
        ReadWord(word);
        node_data << word << '\n'; // value

        for(SizeType i = 0 ; i < NodesAllPartitions[ReorderedNodeId(id)-1].size() ; i++)
        {
            SizeType partition_id = NodesAllPartitions[ReorderedNodeId(id)-1][i];
            if(partition_id > OutputFiles.size())
            {
                std::stringstream buffer;
                buffer << "Invalid prtition id : " << partition_id;
                buffer << " for node " << id << " [Line " << mNumberOfLines << " ]";
                KRATOS_ERROR << buffer.str() << std::endl;;
            }

            *(OutputFiles[partition_id]) << node_data.str();
        }
    }


    KRATOS_CATCH("")
}

void ModelPartIO::DivideVectorialVariableData(OutputFilesContainerType& OutputFiles,
                                    PartitionIndicesContainerType const& EntitiesPartitions,
                                    std::string BlockName)
{
    KRATOS_TRY

    SizeType id;

    std::string word;

    bool is_fixed;

    std::string value;

    while(!mpStream->eof())
    {
        ReadWord(word); // reading id
        if(CheckEndBlock(BlockName, word))
            break;

        ExtractValue(word, id);

        if(ReorderedNodeId(id) > EntitiesPartitions.size())
        {
            std::stringstream buffer;
            buffer << "Invalid node id : " << id;
            buffer << " [Line " << mNumberOfLines << " ]";
            KRATOS_ERROR << buffer.str() << std::endl;;
        }

        std::stringstream entity_data;
        entity_data << ReorderedNodeId(id) << '\t'; // id

        if(BlockName == "NodalData")
        {
            // reading is_fixed
            ReadWord(value);
            ExtractValue(value, is_fixed);
            if(is_fixed)
            {
                std::stringstream buffer;
                buffer << "Only double variables or components can be fixed.";
                buffer <<  " [Line " << mNumberOfLines << " ]";
                KRATOS_ERROR << buffer.str() << std::endl;;
            }
            entity_data << is_fixed << "\t"; // is_fixed
        }

        Vector temp_vector;
        ReadVectorialValue(temp_vector);

        for(SizeType i = 0 ; i < EntitiesPartitions[ReorderedNodeId(id)-1].size() ; i++)
        {
            SizeType partition_id = EntitiesPartitions[ReorderedNodeId(id)-1][i];
            if(partition_id > OutputFiles.size())
            {
                std::stringstream buffer;
                buffer << "Invalid prtition id : " << partition_id;
                buffer << " for node " << id << " [Line " << mNumberOfLines << " ]";
                KRATOS_ERROR << buffer.str() << std::endl;;
            }

            *(OutputFiles[partition_id]) << entity_data.str() << temp_vector << std::endl;
        }
    }


    KRATOS_CATCH("")
}


void ModelPartIO::DivideElementalDataBlock(OutputFilesContainerType& OutputFiles,
                                PartitionIndicesContainerType const& ElementsAllPartitions)
{
    KRATOS_TRY

    typedef VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > array_1d_component_type;

    std::string word;

    WriteInAllFiles(OutputFiles, "Begin ElementalData ");

    std::string variable_name;

    ReadWord(variable_name);

    WriteInAllFiles(OutputFiles, variable_name);
    WriteInAllFiles(OutputFiles, "\n");

    if(KratosComponents<Variable<double> >::Has(variable_name))
    {
        DivideScalarVariableData(OutputFiles, ElementsAllPartitions, "ElementalData");
    }
    else if(KratosComponents<Variable<bool> >::Has(variable_name))
    {
        DivideScalarVariableData(OutputFiles, ElementsAllPartitions, "ElementalData");
    }
    else if(KratosComponents<Variable<int> >::Has(variable_name))
    {
        DivideScalarVariableData(OutputFiles, ElementsAllPartitions, "ElementalData");
    }
    else if(KratosComponents<array_1d_component_type>::Has(variable_name))
    {
        DivideScalarVariableData(OutputFiles, ElementsAllPartitions, "ElementalData");
    }
    else if(KratosComponents<Variable<array_1d<double, 3> > >::Has(variable_name))
    {
        DivideVectorialVariableData(OutputFiles, ElementsAllPartitions, "ElementalData");
    }
    else if(KratosComponents<Variable<Quaternion<double> > >::Has(variable_name))
    {
        DivideVectorialVariableData(OutputFiles, ElementsAllPartitions, "ElementalData");
    }
    else if(KratosComponents<Variable<Matrix> >::Has(variable_name))
    {
        DivideVectorialVariableData(OutputFiles, ElementsAllPartitions, "ElementalData");
    }
    else if(KratosComponents<VariableData>::Has(variable_name))
    {
        std::stringstream buffer;
        buffer << variable_name << " is not supported to be read by this IO or the type of variable is not registered correctly" << std::endl;
        buffer << " [Line " << mNumberOfLines << " ]";
        KRATOS_ERROR << buffer.str() << std::endl;;
    }
    else
    {
        std::stringstream buffer;
        buffer << variable_name << " is not a valid variable!!!" << std::endl;
        buffer << " [Line " << mNumberOfLines << " ]";
        KRATOS_ERROR << buffer.str() << std::endl;;
    }

    WriteInAllFiles(OutputFiles, "End ElementalData\n");

    KRATOS_CATCH("")
}

void ModelPartIO::DivideScalarVariableData(OutputFilesContainerType& OutputFiles,
                                PartitionIndicesContainerType const& EntitiesPartitions,
                                std::string BlockName)
{
    KRATOS_TRY

    SizeType id;

    std::string word;


    while(!mpStream->eof())
    {
        ReadWord(word); // reading id
        if(CheckEndBlock(BlockName, word))
            break;

        ExtractValue(word, id);
        SizeType index = 0;
        if(BlockName == "ElementalData")
            index = ReorderedElementId(id);
        else if(BlockName == "ConditionalData")
            index = ReorderedConditionId(id);
        else
            KRATOS_ERROR << "Invalid block name :" << BlockName << std::endl;

        if(index > EntitiesPartitions.size())
        {
            std::stringstream buffer;
            buffer << "Invalid id : " << id;
            buffer << " [Line " << mNumberOfLines << " ]";
            KRATOS_ERROR << buffer.str() << std::endl;;
        }

        std::stringstream entity_data;
        entity_data << index << '\t'; // id
        ReadWord(word);
        entity_data << word <<'\n'; // value

        for(SizeType i = 0 ; i < EntitiesPartitions[index-1].size() ; i++)
        {
            SizeType partition_id = EntitiesPartitions[index-1][i];
            if(partition_id > OutputFiles.size())
            {
                std::stringstream buffer;
                buffer << "Invalid prtition id : " << partition_id;
                buffer << " for entity " << id << " [Line " << mNumberOfLines << " ]";
                KRATOS_ERROR << buffer.str() << std::endl;;
            }

            *(OutputFiles[partition_id]) << entity_data.str();
        }
    }


    KRATOS_CATCH("")
}


void ModelPartIO::DivideConditionalDataBlock(OutputFilesContainerType& OutputFiles,
                                PartitionIndicesContainerType const& ConditionsAllPartitions)
{
    KRATOS_TRY


    typedef VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > array_1d_component_type;

    std::string word;

    WriteInAllFiles(OutputFiles, "Begin ConditionalData ");

    std::string variable_name;

    ReadWord(variable_name);

    WriteInAllFiles(OutputFiles, variable_name);
    WriteInAllFiles(OutputFiles, "\n");

    if(KratosComponents<Variable<double> >::Has(variable_name))
    {
        DivideScalarVariableData(OutputFiles, ConditionsAllPartitions, "ConditionalData");
    }
    else if(KratosComponents<Variable<bool> >::Has(variable_name))
    {
        DivideScalarVariableData(OutputFiles, ConditionsAllPartitions, "ConditionalData");
    }
    else if(KratosComponents<Variable<int> >::Has(variable_name))
    {
        DivideScalarVariableData(OutputFiles, ConditionsAllPartitions, "ConditionalData");
    }
    else if(KratosComponents<array_1d_component_type>::Has(variable_name))
    {
        DivideScalarVariableData(OutputFiles, ConditionsAllPartitions, "ConditionalData");
    }
    else if(KratosComponents<Variable<array_1d<double, 3> > >::Has(variable_name))
    {
        DivideVectorialVariableData(OutputFiles, ConditionsAllPartitions, "ConditionalData");
    }
    else if(KratosComponents<Variable<Quaternion<double> > >::Has(variable_name))
    {
        DivideVectorialVariableData(OutputFiles, ConditionsAllPartitions, "ConditionalData");
    }
    else if(KratosComponents<Variable<Matrix> >::Has(variable_name))
    {
        DivideVectorialVariableData(OutputFiles, ConditionsAllPartitions, "ConditionalData");
    }
    else if(KratosComponents<VariableData>::Has(variable_name))
    {
        std::stringstream buffer;
        buffer << variable_name << " is not supported to be read by this IO or the type of variable is not registered correctly" << std::endl;
        buffer << " [Line " << mNumberOfLines << " ]";
        KRATOS_ERROR << buffer.str() << std::endl;;
    }
    else
    {
        std::stringstream buffer;
        buffer << variable_name << " is not a valid variable!!!" << std::endl;
        buffer << " [Line " << mNumberOfLines << " ]";
        KRATOS_ERROR << buffer.str() << std::endl;;
    }

    WriteInAllFiles(OutputFiles, "End ConditionalData\n");


    KRATOS_CATCH("")
}


void ModelPartIO::DivideMeshBlock(OutputFilesContainerType& OutputFiles,
                                        PartitionIndicesContainerType const& NodesAllPartitions,
                                        PartitionIndicesContainerType const& ElementsAllPartitions,
                                        PartitionIndicesContainerType const& ConditionsAllPartitions)
{
    KRATOS_TRY

    std::string word;
    ReadWord(word);

    word += "\n";


    WriteInAllFiles(OutputFiles, "Begin Mesh " + word);

    while(!mpStream->eof())
    {
        ReadWord(word);

        if(CheckEndBlock("Mesh", word))
            break;

        ReadBlockName(word);
        if(word == "MeshData")
            DivideMeshDataBlock(OutputFiles);
        else if(word == "MeshNodes")
            DivideMeshNodesBlock(OutputFiles, NodesAllPartitions);
        else if(word == "MeshElements")
            DivideMeshElementsBlock(OutputFiles, ElementsAllPartitions);
        else if(word == "MeshConditions")
            DivideMeshConditionsBlock(OutputFiles, ConditionsAllPartitions);
        else
            SkipBlock(word);
    }

    WriteInAllFiles(OutputFiles, "End Mesh\n");

    KRATOS_CATCH("")

}

void ModelPartIO::DivideSubModelPartBlock(OutputFilesContainerType& OutputFiles,
    PartitionIndicesContainerType const& NodesAllPartitions,
    PartitionIndicesContainerType const& ElementsAllPartitions,
    PartitionIndicesContainerType const& ConditionsAllPartitions)
{
    KRATOS_TRY

        std::string word;
    ReadWord(word);

    word += "\n";


    WriteInAllFiles(OutputFiles, "Begin SubModelPart " + word);

    while (!mpStream->eof())
    {
        ReadWord(word);

        if (CheckEndBlock("SubModelPart", word))
            break;

        ReadBlockName(word);
        if (word == "SubModelPartData")
            DivideSubModelPartDataBlock(OutputFiles);
        else if (word == "SubModelPartTables")
            DivideSubModelPartTableBlock(OutputFiles);
        else if (word == "SubModelPartNodes")
            DivideSubModelPartNodesBlock(OutputFiles, NodesAllPartitions);
        else if (word == "SubModelPartElements")
            DivideSubModelPartElementsBlock(OutputFiles, ElementsAllPartitions);
        else if (word == "SubModelPartConditions")
            DivideSubModelPartConditionsBlock(OutputFiles, ConditionsAllPartitions);
        else if (word == "SubModelPart")
            DivideSubModelPartBlock(OutputFiles, NodesAllPartitions, ElementsAllPartitions, ConditionsAllPartitions);
        else
            SkipBlock(word);
    }

    WriteInAllFiles(OutputFiles, "End SubModelPart\n");

    KRATOS_CATCH("")
}



void ModelPartIO::DivideMeshDataBlock(OutputFilesContainerType& OutputFiles)
{
    KRATOS_TRY
    std::string block;

    WriteInAllFiles(OutputFiles, "Begin MeshData");

    ReadBlock(block, "MeshData");
    WriteInAllFiles(OutputFiles, block);

    WriteInAllFiles(OutputFiles, "End MeshData\n");
    KRATOS_CATCH("")
}


void ModelPartIO::DivideMeshNodesBlock(OutputFilesContainerType& OutputFiles,
                                        PartitionIndicesContainerType const& NodesAllPartitions)
{
    KRATOS_TRY

    std::string word;

    WriteInAllFiles(OutputFiles, "Begin MeshNodes \n");

    SizeType id;

    while(!mpStream->eof())
    {
        ReadWord(word);

        if(CheckEndBlock("MeshNodes", word))
            break;

        ExtractValue(word, id);

        if(ReorderedNodeId(id) > NodesAllPartitions.size())
        {
            std::stringstream buffer;
            buffer << "Invalid node id : " << id;
            buffer << " [Line " << mNumberOfLines << " ]";
            KRATOS_ERROR << buffer.str() << std::endl;;
        }

        for(SizeType i = 0 ; i < NodesAllPartitions[ReorderedNodeId(id)-1].size() ; i++)
        {
            SizeType partition_id = NodesAllPartitions[ReorderedNodeId(id)-1][i];
            if(partition_id > OutputFiles.size())
            {
                std::stringstream buffer;
                buffer << "Invalid prtition id : " << partition_id;
                buffer << " for node " << id << " [Line " << mNumberOfLines << " ]";
                KRATOS_ERROR << buffer.str() << std::endl;;
            }

            *(OutputFiles[partition_id]) << ReorderedNodeId(id) << std::endl;
        }

    }

    WriteInAllFiles(OutputFiles, "End MeshNodes\n");

    KRATOS_CATCH("")
}


void ModelPartIO::DivideMeshElementsBlock(OutputFilesContainerType& OutputFiles,
                                        PartitionIndicesContainerType const& ElementsAllPartitions)
{
    KRATOS_TRY

    std::string word;

    WriteInAllFiles(OutputFiles, "Begin MeshElements \n");

    SizeType id;

    while(!mpStream->eof())
    {
        ReadWord(word);

        if(CheckEndBlock("MeshElements", word))
            break;

        ExtractValue(word, id);

        if(ReorderedElementId(id) > ElementsAllPartitions.size())
        {
            std::stringstream buffer;
            buffer << "Invalid element id : " << id;
            buffer << " [Line " << mNumberOfLines << " ]";
            KRATOS_ERROR << buffer.str() << std::endl;;
        }

        for(SizeType i = 0 ; i < ElementsAllPartitions[ReorderedElementId(id)-1].size() ; i++)
        {
            SizeType partition_id = ElementsAllPartitions[ReorderedElementId(id)-1][i];
            if(partition_id > OutputFiles.size())
            {
                std::stringstream buffer;
                buffer << "Invalid prtition id : " << partition_id;
                buffer << " for element " << id << " [Line " << mNumberOfLines << " ]";
                KRATOS_ERROR << buffer.str() << std::endl;;
            }

            *(OutputFiles[partition_id]) << ReorderedElementId(id) << std::endl;
        }

    }

    WriteInAllFiles(OutputFiles, "End MeshElements\n");

    KRATOS_CATCH("")
}

void ModelPartIO::DivideMeshConditionsBlock(OutputFilesContainerType& OutputFiles,
                                        PartitionIndicesContainerType const& ConditionsAllPartitions)
{
    KRATOS_TRY

    std::string word;

    WriteInAllFiles(OutputFiles, "Begin MeshConditions \n");

    SizeType id;

    while(!mpStream->eof())
    {
        ReadWord(word);

        if(CheckEndBlock("MeshConditions", word))
            break;

        ExtractValue(word, id);

        if(ReorderedConditionId(id) > ConditionsAllPartitions.size())
        {
            std::stringstream buffer;
            buffer << "Invalid condition id : " << id;
            buffer << " [Line " << mNumberOfLines << " ]";
            KRATOS_ERROR << buffer.str() << std::endl;;
        }

        for(SizeType i = 0 ; i < ConditionsAllPartitions[ReorderedConditionId(id)-1].size() ; i++)
        {
            SizeType partition_id = ConditionsAllPartitions[ReorderedConditionId(id)-1][i];
            if(partition_id > OutputFiles.size())
            {
                std::stringstream buffer;
                buffer << "Invalid prtition id : " << partition_id;
                buffer << " for condition " << id << " [Line " << mNumberOfLines << " ]";
                KRATOS_ERROR << buffer.str() << std::endl;;
            }

            *(OutputFiles[partition_id]) << ReorderedConditionId(id) << std::endl;
        }

    }

    WriteInAllFiles(OutputFiles, "End MeshConditions\n");

    KRATOS_CATCH("")
}




void ModelPartIO::DivideSubModelPartDataBlock(OutputFilesContainerType& OutputFiles)
{
    KRATOS_TRY
        std::string block;

    WriteInAllFiles(OutputFiles, "Begin SubModelPartData");

    ReadBlock(block, "SubModelPartData");
    WriteInAllFiles(OutputFiles, block);

    WriteInAllFiles(OutputFiles, "End SubModelPartData\n");
    KRATOS_CATCH("")
}

void ModelPartIO::DivideSubModelPartTableBlock(OutputFilesContainerType& OutputFiles)
{
    KRATOS_TRY
        std::string block;

    WriteInAllFiles(OutputFiles, "Begin SubModelPartTables");

    ReadBlock(block, "SubModelPartTables");
    WriteInAllFiles(OutputFiles, block);

    WriteInAllFiles(OutputFiles, "End SubModelPartTables\n");
    KRATOS_CATCH("")
}


void ModelPartIO::DivideSubModelPartNodesBlock(OutputFilesContainerType& OutputFiles,
    PartitionIndicesContainerType const& NodesAllPartitions)
{
    KRATOS_TRY

        std::string word;

    WriteInAllFiles(OutputFiles, "Begin SubModelPartNodes \n");

    SizeType id;

    while (!mpStream->eof())
    {
        ReadWord(word);

        if (CheckEndBlock("SubModelPartNodes", word))
            break;

        ExtractValue(word, id);

        if (ReorderedNodeId(id) > NodesAllPartitions.size())
        {
            std::stringstream buffer;
            buffer << "Invalid node id : " << id;
            buffer << " [Line " << mNumberOfLines << " ]";
            KRATOS_ERROR << buffer.str() << std::endl;;
        }

        for (SizeType i = 0; i < NodesAllPartitions[ReorderedNodeId(id) - 1].size(); i++)
        {
            SizeType partition_id = NodesAllPartitions[ReorderedNodeId(id) - 1][i];
            if (partition_id > OutputFiles.size())
            {
                std::stringstream buffer;
                buffer << "Invalid prtition id : " << partition_id;
                buffer << " for node " << id << " [Line " << mNumberOfLines << " ]";
                KRATOS_ERROR << buffer.str() << std::endl;;
            }

            *(OutputFiles[partition_id]) << ReorderedNodeId(id) << std::endl;
        }

    }

    WriteInAllFiles(OutputFiles, "End SubModelPartNodes\n");

    KRATOS_CATCH("")
}


void ModelPartIO::DivideSubModelPartElementsBlock(OutputFilesContainerType& OutputFiles,
    PartitionIndicesContainerType const& ElementsAllPartitions)
{
    KRATOS_TRY

        std::string word;

    WriteInAllFiles(OutputFiles, "Begin SubModelPartElements \n");

    SizeType id;

    while (!mpStream->eof())
    {
        ReadWord(word);

        if (CheckEndBlock("SubModelPartElements", word))
            break;

        ExtractValue(word, id);

        if (ReorderedElementId(id) > ElementsAllPartitions.size())
        {
            std::stringstream buffer;
            buffer << "Invalid element id : " << id;
            buffer << " [Line " << mNumberOfLines << " ]";
            KRATOS_ERROR << buffer.str() << std::endl;;
        }

        for (SizeType i = 0; i < ElementsAllPartitions[ReorderedElementId(id) - 1].size(); i++)
        {
            SizeType partition_id = ElementsAllPartitions[ReorderedElementId(id) - 1][i];
            if (partition_id > OutputFiles.size())
            {
                std::stringstream buffer;
                buffer << "Invalid prtition id : " << partition_id;
                buffer << " for element " << id << " [Line " << mNumberOfLines << " ]";
                KRATOS_ERROR << buffer.str() << std::endl;;
            }

            *(OutputFiles[partition_id]) << ReorderedElementId(id) << std::endl;
        }

    }

    WriteInAllFiles(OutputFiles, "End SubModelPartElements\n");

    KRATOS_CATCH("")
}

void ModelPartIO::DivideSubModelPartConditionsBlock(OutputFilesContainerType& OutputFiles,
    PartitionIndicesContainerType const& ConditionsAllPartitions)
{
    KRATOS_TRY

        std::string word;

    WriteInAllFiles(OutputFiles, "Begin SubModelPartConditions \n");

    SizeType id;

    while (!mpStream->eof())
    {
        ReadWord(word);

        if (CheckEndBlock("SubModelPartConditions", word))
            break;

        ExtractValue(word, id);

        if (ReorderedConditionId(id) > ConditionsAllPartitions.size())
        {
            std::stringstream buffer;
            buffer << "Invalid condition id : " << id;
            buffer << " [Line " << mNumberOfLines << " ]";
            KRATOS_ERROR << buffer.str() << std::endl;;
        }

        for (SizeType i = 0; i < ConditionsAllPartitions[ReorderedConditionId(id) - 1].size(); i++)
        {
            SizeType partition_id = ConditionsAllPartitions[ReorderedConditionId(id) - 1][i];
            if (partition_id > OutputFiles.size())
            {
                std::stringstream buffer;
                buffer << "Invalid prtition id : " << partition_id;
                buffer << " for condition " << id << " [Line " << mNumberOfLines << " ]";
                KRATOS_ERROR << buffer.str() << std::endl;;
            }

            *(OutputFiles[partition_id]) << ReorderedConditionId(id) << std::endl;
        }

    }

    WriteInAllFiles(OutputFiles, "End SubModelPartConditions\n");

    KRATOS_CATCH("")
}

void ModelPartIO::WritePartitionIndices(OutputFilesContainerType& OutputFiles, PartitionIndicesType const&  NodesPartitions, PartitionIndicesContainerType const& NodesAllPartitions)
{
    WriteInAllFiles(OutputFiles, "Begin NodalData PARTITION_INDEX\n");

    for(SizeType i_node = 0 ; i_node != NodesAllPartitions.size() ; i_node++)
    {
        for(SizeType i = 0 ; i < NodesAllPartitions[i_node].size() ; i++)
        {
            SizeType partition_id = NodesAllPartitions[i_node][i];
            if(partition_id > OutputFiles.size())
            {
                std::stringstream buffer;
                buffer << "Invalid prtition id : " << partition_id;
                buffer << " for node " << i_node+1 << " [Line " << mNumberOfLines << " ]";
                KRATOS_ERROR << buffer.str() << std::endl;;
            }

            const SizeType node_partition = NodesPartitions[i_node];
            *(OutputFiles[partition_id]) << i_node + 1 << "  0  " << node_partition << std::endl;
        }
    }


    WriteInAllFiles(OutputFiles, "End NodalData \n");

}


void ModelPartIO::WriteCommunicatorData(OutputFilesContainerType& OutputFiles, SizeType NumberOfPartitions, GraphType const& DomainsColoredGraph,
                            PartitionIndicesType const& NodesPartitions,
                            PartitionIndicesType const& ElementsPartitions,
                            PartitionIndicesType const& ConditionsPartitions,
                            PartitionIndicesContainerType const& NodesAllPartitions,
                            PartitionIndicesContainerType const& ElementsAllPartitions,
                            PartitionIndicesContainerType const& ConditionsAllPartitions)
{
    WriteInAllFiles(OutputFiles, "Begin CommunicatorData \n");


    // Writing the domains neighbours
    WriteInAllFiles(OutputFiles, "NEIGHBOURS_INDICES    ");
    for(SizeType i_partition = 0 ; i_partition < NumberOfPartitions ; i_partition++)
    {
        DenseVector<int> indices = row(DomainsColoredGraph, i_partition);
        *(OutputFiles[i_partition]) << indices << std::endl;
    }

    SizeType number_of_colors = 0;

    for(SizeType i_partition = 0 ; i_partition < DomainsColoredGraph.size1() ; i_partition++)
        for(SizeType i_interface = 0 ; i_interface < DomainsColoredGraph.size2() ; i_interface++)
            if(DomainsColoredGraph(i_partition, i_interface) >= 0)
                if(number_of_colors < i_interface)
                    number_of_colors = i_interface;

    number_of_colors++; // I have to add one to it to get the correct number

    // Writing the max colors
    for(SizeType i_partition = 0 ; i_partition < NumberOfPartitions ; i_partition++)
    {
        DenseVector<int> indices = row(DomainsColoredGraph, i_partition);
        *(OutputFiles[i_partition]) << "NUMBER_OF_COLORS    " << number_of_colors << std::endl;
    }


    // Writing the all local nodes
    WriteInAllFiles(OutputFiles, "    Begin LocalNodes 0\n");

    for(SizeType i = 0 ; i < NodesPartitions.size() ; i++)
        *(OutputFiles[NodesPartitions[i]]) << "    " << i+1 << std::endl;

    WriteInAllFiles(OutputFiles, "    End LocalNodes \n");


    std::vector<PartitionIndicesContainerType> local_nodes_indices(NumberOfPartitions, PartitionIndicesContainerType(number_of_colors));
    std::vector<PartitionIndicesContainerType> ghost_nodes_indices(NumberOfPartitions, PartitionIndicesContainerType(number_of_colors));

    DenseMatrix<int> interface_indices = scalar_matrix<int>(NumberOfPartitions, NumberOfPartitions, -1);

    for(SizeType i_partition = 0 ; i_partition < NumberOfPartitions ; i_partition++)
    {
        DenseVector<int> neighbours_indices = row(DomainsColoredGraph, i_partition);

        for(SizeType i = 0 ; i <  neighbours_indices.size() ; i++)
            if(SizeType(neighbours_indices[i]) < NumberOfPartitions)
                interface_indices(i_partition,neighbours_indices[i]) = i;
    }


    for(SizeType i = 0 ; i < NodesPartitions.size() ; i++)
    {
        const SizeType node_partition = NodesPartitions[i];
        const SizeType node_id = i + 1;

        PartitionIndicesType const& node_all_partitions = NodesAllPartitions[i];

        for(SizeType j = 0 ; j < node_all_partitions.size() ; j++)
        {
            SizeType i_node_partition = node_all_partitions[j];
            if(node_partition != i_node_partition)
            {
                SizeType local_interface_index = interface_indices(node_partition, i_node_partition);
                SizeType ghost_interface_index = interface_indices(i_node_partition, node_partition);
                local_nodes_indices[node_partition][local_interface_index].push_back(node_id);
                ghost_nodes_indices[i_node_partition][ghost_interface_index].push_back(node_id);
            }
        }

    }

    for(SizeType i_partition = 0 ; i_partition < NumberOfPartitions ; i_partition++)
    {
        PartitionIndicesContainerType& partition_local_nodes_indices = local_nodes_indices[i_partition];

        for(SizeType i_interface = 0 ; i_interface < partition_local_nodes_indices.size() ; i_interface++)
        {
            if(partition_local_nodes_indices[i_interface].size() > 0)
            {
                *(OutputFiles[i_partition]) << "    Begin LocalNodes " << i_interface + 1 << std::endl;
                for(SizeType i_interface_node = 0 ; i_interface_node < partition_local_nodes_indices[i_interface].size() ; i_interface_node++)
                    *(OutputFiles[i_partition]) << "    " << partition_local_nodes_indices[i_interface][i_interface_node] << std::endl;
                *(OutputFiles[i_partition]) << "    End LocalNodes " << std::endl;
            }
        }

        PartitionIndicesContainerType& partition_ghost_nodes_indices = ghost_nodes_indices[i_partition];

        std::set<unsigned int> all_ghost_nodes_indices;

        for(SizeType i_interface = 0 ; i_interface < partition_ghost_nodes_indices.size() ; i_interface++)
        {
            if(partition_ghost_nodes_indices[i_interface].size() > 0)
            {
                *(OutputFiles[i_partition]) << "    Begin GhostNodes " << i_interface + 1 << std::endl;
                for(SizeType i_interface_node = 0 ; i_interface_node < partition_ghost_nodes_indices[i_interface].size() ; i_interface_node++)
                {
                    *(OutputFiles[i_partition]) << "    " << partition_ghost_nodes_indices[i_interface][i_interface_node] << std::endl;
                    all_ghost_nodes_indices.insert(partition_ghost_nodes_indices[i_interface][i_interface_node]);
                }
                *(OutputFiles[i_partition]) << "    End GhostNodes "  << std::endl;
            }
        }

        *(OutputFiles[i_partition]) << "    Begin GhostNodes " << 0 << std::endl;
        for(std::set<unsigned int>::iterator id = all_ghost_nodes_indices.begin() ; id != all_ghost_nodes_indices.end() ; id++)
        {
            *(OutputFiles[i_partition]) << "    " << *id << std::endl;
        }
        *(OutputFiles[i_partition]) << "    End GhostNodes "  << std::endl;

    }

    WriteInAllFiles(OutputFiles, "End CommunicatorData \n");
}

void ModelPartIO::WriteCommunicatorLocalNodes(OutputFilesContainerType& OutputFiles, SizeType NumberOfPartitions, PartitionIndicesType const& NodesPartitions, PartitionIndicesContainerType const& NodesAllPartitions)
{
    WriteInAllFiles(OutputFiles, "    Begin LocalNodes 0\n");

    for(SizeType i = 0 ; i < NodesPartitions.size() ; i++)
        *(OutputFiles[NodesPartitions[i]]) << "    " << i+1 << std::endl;

    WriteInAllFiles(OutputFiles, "    End LocalNodes \n");

    PartitionIndicesContainerType local_nodes_indices(NumberOfPartitions);


}

void ModelPartIO::WriteInAllFiles(OutputFilesContainerType& OutputFiles, std::string const& ThisWord)
{
    for(SizeType i = 0 ; i < OutputFiles.size() ; i++)
        *(OutputFiles[i]) << ThisWord;
}


template<class TContainerType, class TKeyType>
typename TContainerType::iterator ModelPartIO::FindKey(TContainerType& ThisContainer , TKeyType ThisKey, std::string ComponentName)
{
    typename TContainerType::iterator i_result;
    if((i_result = ThisContainer.find(ThisKey)) == ThisContainer.end())
    {
        std::stringstream buffer;
        buffer << ComponentName << " #" << ThisKey << " is not found.";
        buffer << " [Line " << mNumberOfLines << " ]";
        KRATOS_ERROR << buffer.str() << std::endl;;
    }

    return i_result;
}



/**
* @note Basically it starts to read the character sequence until reaching a
*       "(" and then goes until corresponding ")" which means the vector or
*       matrix value is completely read. It can be used to read any kind of
*       vector or matrix with operator >> defined and writtern in following
*       format for a vector: [size] ( value1, value2,...., valueN )
*       format for a matrix: [size1,size2] (( )( )...( )) //look props read
*/
template<class TValueType>
TValueType& ModelPartIO::ReadVectorialValue(TValueType& rValue)
{
    std::stringstream value;

    char c = SkipWhiteSpaces();
    while((c != '(') && !mpStream->eof())
    {
        value << c;
        c = GetCharacter();
    }
    int open_parantesis = 1;
    while((open_parantesis != 0) && !mpStream->eof())
    {
        value << c;
        c = GetCharacter();
        if(c == '(')
            open_parantesis++;
        if(c == ')')
            open_parantesis--;
    }
    value << c; // adding the final parantesis

    value >>  rValue;

    return rValue;
}

template<class TValueType>
TValueType& ModelPartIO::ExtractValue(std::string rWord, TValueType & rValue)
{
    std::stringstream value(rWord);

    value >> rValue;

    return rValue;
}

void ModelPartIO::ReadConstitutiveLawValue(ConstitutiveLaw::Pointer& rValue) {
    std::string value;
    ReadWord(value);
    rValue = KratosComponents<ConstitutiveLaw>::Get(value).Clone();
}


ModelPartIO& ModelPartIO::ReadWord(std::string& Word)
{
    Word.clear();

    char c = SkipWhiteSpaces();
    while(!mpStream->eof() && !IsWhiteSpace(c))
    {
        Word += c;
        c = GetCharacter();
    }

    return *this;
}

ModelPartIO& ModelPartIO::ReadBlock(std::string& Block, std::string const& BlockName)
{
    Block.clear();
    std::vector<std::string> nested_block_names;
    nested_block_names.push_back(BlockName);

    char c = GetCharacter();
    std::string word;

    while(!mpStream->eof())
    {
        if(c == 'E')
        {
            word.clear();
            while(!mpStream->eof() && !IsWhiteSpace(c))
            {
                word += c;
                c = GetCharacter();
            }
            if (CheckEndBlock(nested_block_names.back(), word))
            {
                nested_block_names.pop_back();
                if(nested_block_names.empty())
                {
                    break;
                }
                else
                {
                    Block += "End ";
                }
            }

            Block += word;
        }
        else if (c == 'B')
        {
            word.clear();
            while (!mpStream->eof() && !IsWhiteSpace(c))
            {
                word += c;
                c = GetCharacter();
            }
            if (word == "Begin")
            {
                Block += word;
                Block += c;
                ReadWord(word);
                nested_block_names.push_back(word);
            }

            Block += word;
        }

        Block += c;

        c = GetCharacter();
    }

    return *this;
}


char ModelPartIO::SkipWhiteSpaces()
{
    char c = GetCharacter();
    while(IsWhiteSpace(c))
        c = GetCharacter();
    return c;
}

bool ModelPartIO::IsWhiteSpace(char C)
{
    return ((C == ' ') || (C == '\t') || (C == '\r') || (C == '\n'));
}

char ModelPartIO::GetCharacter() //Getting the next character skipping comments
{
    char c;
    if(mpStream->get(c))
    {
        if(c == '\n')
            mNumberOfLines++;
        else if(c == '/') // it may be a comment!
        {
            char next_c = mpStream->peek();
            if(next_c == '/') // it's a line comment
            {
                while((mpStream->get(c)) && (c != '\n')); // so going to the end of line
                if(!mpStream->eof())
                    mNumberOfLines++;
            }
            else if(next_c == '*') // it's a block comment
            {
                while((mpStream->get(c)) && !((c == '*') && (mpStream->peek() == '/'))) // so going to the end of block
                    if(c == '\n')
                        mNumberOfLines++;
                mpStream->get(c);
                c = GetCharacter(); // read a new character after comment
            }
        }
    }
    else
        c = 0;

    return c;

}

bool ModelPartIO::CheckStatement(std::string const& rStatement, std::string const& rGivenWord)
{
    bool result = false;
    if(rGivenWord != rStatement)
    {
        std::stringstream buffer;
        buffer << "A \"" << rStatement << "\" statement was expected but the given statement was \"";
        buffer <<  rGivenWord << "\"" << " [Line " << mNumberOfLines << " ]";
        KRATOS_ERROR << buffer.str() << std::endl;;
    }
    else
        result = true;

    return result;
}

void ModelPartIO::ResetInput()
{
    mpStream->clear();
    mpStream->seekg(0, std::ios_base::beg);
    mNumberOfLines = 1;
}

void ModelPartIO::SwapStreamSource(Kratos::shared_ptr<std::iostream> newStream)
{
    mpStream.swap(newStream);
}

inline void ModelPartIO::CreatePartition(unsigned int number_of_threads,const int number_of_rows, DenseVector<unsigned int>& partitions)
{
    partitions.resize(number_of_threads+1);
    int partition_size = number_of_rows / number_of_threads;
    partitions[0] = 0;
    partitions[number_of_threads] = number_of_rows;
    for(unsigned int i = 1; i<number_of_threads; i++)
        partitions[i] = partitions[i-1] + partition_size ;
}

/// Unaccessible assignment operator.
ModelPartIO& ModelPartIO::operator=(ModelPartIO const& rOther){return *this;}

/// Unaccessible copy constructor.
ModelPartIO::ModelPartIO(ModelPartIO const& rOther){}

ModelPartIO::SizeType ModelPartIO::ReorderedNodeId(ModelPartIO::SizeType NodeId)
{
    // The ModelPartIO does not reorder the nodes
    // This method is the one to be overriden by some reordering IO class
    return NodeId;
}
ModelPartIO::SizeType ModelPartIO::ReorderedElementId(ModelPartIO::SizeType ElementId)
{
    // The ModelPartIO does not reorder the elements
    // This method is the one to be overriden by some reordering IO class
    return ElementId;
}

ModelPartIO::SizeType ModelPartIO::ReorderedConditionId(ModelPartIO::SizeType ConditionId)
{
    // The ModelPartIO does not reorder the conditions
    // This method is the one to be overriden by some reordering IO class
    return ConditionId;
}

}  // namespace Kratos.
