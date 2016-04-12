//    |  /           | 
//    ' /   __| _` | __|  _ \   __| 
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/ 
//                   Multi-Physics  
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//












// Project includes
#include "includes/model_part_io.h"
#include "utilities/logger.h"



namespace Kratos
{
	

    /// Constructor with  filenames.
    ModelPartIO::ModelPartIO(std::string const& Filename, const Flags Options )
        : mNumberOfLines(1)
        , mBaseFilename(Filename)
        , mFilename(Filename + ".mdpa")
        , mOptions(Options)
    {
		if (mOptions.Is(IO::READ))
		{
			mFile.open(mFilename.c_str(), std::fstream::in);
			if(!(mFile.is_open()))
				KRATOS_THROW_ERROR(std::invalid_argument, "Error opening input file : ", mFilename.c_str());
		}
		else if (mOptions.Is(IO::APPEND))
		{
			mFile.open(mFilename.c_str(), std::fstream::in|std::fstream::app);
			if(!(mFile.is_open()))
				KRATOS_THROW_ERROR(std::invalid_argument, "Error opening input file : ", mFilename.c_str());
		}
		else if(mOptions.Is(IO::WRITE))
		{
			mFile.open(mFilename.c_str(), std::fstream::out);
			if(!(mFile.is_open()))
				KRATOS_THROW_ERROR(std::invalid_argument, "Error opening output file : ", mFilename.c_str());
		}
		else // if none of the READ, WRITE or APPEND are defined we will take READ as defualt.
		{
			mFile.open(mFilename.c_str(), std::fstream::in);
			if(!(mFile.is_open()))
				KRATOS_THROW_ERROR(std::invalid_argument, "Error opening input file : ", mFilename.c_str());
		}


        Timer::SetOuputFile(Filename + ".time");

    }


    /// Destructor.
    ModelPartIO::~ModelPartIO() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    bool ModelPartIO::ReadNode(NodeType& rThisNode)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class member. Please check the definition of derived class.", "")
    }

    bool ModelPartIO::ReadNodes(NodesContainerType& rThisNodes)
    {
        KRATOS_TRY
        ResetInput();
        std::string word;
        while(true)
        {
            ReadWord(word);
            if(mFile.eof())
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
            if (mFile.eof())
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
        mFile << "Begin Nodes" << std::endl;
        for(NodesContainerType::const_iterator i_node = rThisNodes.begin() ; i_node != rThisNodes.end() ; i_node++)
            mFile << i_node->Id() << "\t" << i_node->X()  << "\t" << i_node->Y() << "\t" << i_node->Z() << std::endl;
        mFile << "End Nodes" << std::endl;
    }

    void ModelPartIO::ReadProperties(Properties& rThisProperties)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class member. Please check the definition of derived class", "")
    }

    void ModelPartIO::ReadProperties(PropertiesContainerType& rThisProperties)
    {
        KRATOS_TRY
        ResetInput();
        std::string word;
        while(true)
        {
            ReadWord(word);
            if(mFile.eof())
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
        for(PropertiesContainerType::const_iterator i_properties = rThisProperties.begin() ; i_properties != rThisProperties.end() ; i_properties++)
        {
            mFile << "Begin Properties " << i_properties->Id() << std::endl;
            i_properties->PrintData(mFile);
            mFile << std::endl;
            mFile << "End Properties" << std::endl;
        }
    }

    void ModelPartIO::ReadElement(NodesContainerType& rThisNodes, PropertiesContainerType& rThisProperties, Element::Pointer& pThisElements)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Calling base class member. Please check the definition of derived class", "")
    }

    void ModelPartIO::ReadElements(NodesContainerType& rThisNodes, PropertiesContainerType& rThisProperties, ElementsContainerType& rThisElements)
    {
        KRATOS_TRY
        ResetInput();
        std::string word;
        while(true)
        {
            ReadWord(word);
            if(mFile.eof())
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
            if(mFile.eof())
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

    void ModelPartIO::WriteElements(ElementsContainerType const& rThisElements)
    {
		KRATOS_THROW_ERROR(std::logic_error, "This method has not been implemented yet!","");
        mFile << "Begin Elements" << std::endl;
        mFile << "End Elements" << std::endl;
    }

    void ModelPartIO::ReadConditions(NodesContainerType& rThisNodes, PropertiesContainerType& rThisProperties, ConditionsContainerType& rThisConditions)
    {
        KRATOS_TRY
        ResetInput();
        std::string word;
        while(true)
        {
            ReadWord(word);
            if(mFile.eof())
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
            if(mFile.eof())
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
		KRATOS_THROW_ERROR(std::logic_error, "This method has not been implemented yet!","");
        mFile << "Begin Conditions" << std::endl;
        mFile << "End Conditions" << std::endl;
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
            if(mFile.eof())
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
        KRATOS_THROW_ERROR(std::logic_error, "ModelPartIO does not implement this method.", "")
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

        Timer::Start("Reading Input");

        ResetInput();
        std::string word;
        while(true)
        {
            ReadWord(word);
            if(mFile.eof())
                break;
            ReadBlockName(word);
            if(word == "ModelPartData")
                ReadModelPartDataBlock(rThisModelPart);
            else if(word == "Table")
                ReadTableBlock(rThisModelPart.Tables());
            else if(word == "Properties")
                ReadPropertiesBlock(rThisModelPart.rProperties());
            else if(word == "Nodes")
                ReadNodesBlock(rThisModelPart);
            else if(word == "Elements")
                ReadElementsBlock(rThisModelPart);
            else if(word == "Conditions")
                ReadConditionsBlock(rThisModelPart);
            else if(word == "NodalData")
                ReadNodalDataBlock(rThisModelPart);
            else if(word == "ElementalData")
                ReadElementalDataBlock(rThisModelPart.Elements());
			else if (word == "ConditionalData")
				ReadConditionalDataBlock(rThisModelPart.Conditions());
			else if(word == "CommunicatorData")
            {
                ReadCommunicatorDataBlock(rThisModelPart.GetCommunicator(), rThisModelPart.Nodes());
                //Adding the elements and conditions to the communicator
                rThisModelPart.GetCommunicator().LocalMesh().Elements() = rThisModelPart.Elements();
                rThisModelPart.GetCommunicator().LocalMesh().Conditions() = rThisModelPart.Conditions();
            }
			else if (word == "Mesh")
				ReadMeshBlock(rThisModelPart);
			else if (word == "SubModelPart")
				ReadSubModelPartBlock(rThisModelPart, rThisModelPart);
		}
        std::cout << "  [Total Lines Read : " << mNumberOfLines<<"]";
        std::cout << std::endl;
	Timer::Stop("Reading Input");
        KRATOS_CATCH("")
    }



    void ModelPartIO::WriteModelPart(ModelPart & rThisModelPart)
    {
		KRATOS_THROW_ERROR(std::logic_error, "This method has not been implemented yet!","");
        mFile << "Begin ModelPartData" << std::endl;
        mFile << "End ModelPartData" << std::endl;
        WriteMesh(rThisModelPart.GetMesh());
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
    std::size_t ModelPartIO::ReadNodalGraph(int **NodeIndices, int **NodeConnectivities)
    {
        // 1. Define an auxiliary vector of vectors
        ConnectivitiesContainerType aux_connectivities(0);

        // 2. Fill the auxiliary vector by reading elemental and conditional connectivities
        ResetInput();
        std::string word;
        while(true)
        {
            ReadWord(word);
            if(mFile.eof())
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
                KRATOS_THROW_ERROR(std::runtime_error, msg.str(),"")
            }
        }

        // 3. Sort each entry in the auxiliary connectivities vector, remove duplicates
        SizeType num_entries = 0;
        for (ConnectivitiesContainerType::iterator it = aux_connectivities.begin(); it != aux_connectivities.end(); ++it)
        {
            std::sort(it->begin(),it->end());
            std::vector<SizeType>::iterator unique_end = std::unique(it->begin(),it->end());
            it->resize(unique_end - it->begin());
            num_entries += it->size();
        }

        // 4. Write connectivity data in CSR format
        SizeType num_nodes = aux_connectivities.size();
        *NodeIndices = new int[num_nodes+1];
        (*NodeIndices)[0] = 0;
        *NodeConnectivities = new int[num_entries];

        SizeType i = 0;
        SizeType aux_index = 0;

        for (ConnectivitiesContainerType::iterator it = aux_connectivities.begin(); it != aux_connectivities.end(); ++it)
        {
            for (std::vector<SizeType>::iterator entry_it = it->begin(); entry_it != it->end(); entry_it++)
                (*NodeConnectivities)[aux_index++] = (*entry_it - 1); // substract 1 to make Ids start from 0
            (*NodeIndices)[++i] = aux_index;
        }

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
            if(!(*p_ofstream))
                KRATOS_THROW_ERROR(std::invalid_argument, "Error opening output file : ", buffer.str());

            output_files.push_back(p_ofstream);
        }

        while(true)
        {
            ReadWord(word);
            if(mFile.eof())
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
        std::cout << "  [Total Lines Read : " << mNumberOfLines<<"]";
        std::cout << std::endl;

        for(SizeType i = 0 ; i < NumberOfPartitions ; i++)
            delete output_files[i];
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

        while(!mFile.eof())
        {
            ReadWord(word);
            if(word == "End")
            {
                ReadWord(word);
                if(word == BlockName)
                    break;
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

    void ModelPartIO::ReadModelPartDataBlock(ModelPart& rModelPart)
    {
        KRATOS_TRY

        std::string variable_name;

        while(!mFile.eof())
        {
            ReadWord(variable_name);
            if(CheckEndBlock("ModelPartData", variable_name))
                break;
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
            else if(KratosComponents<Variable<Matrix> >::Has(variable_name))
            {
                ReadVectorialValue(rModelPart[KratosComponents<Variable<Matrix> >::Get(variable_name)]);
            }
            else
            {
                std::stringstream buffer;
                buffer << variable_name << " is not a valid variable!!!" << std::endl;
                buffer << " [Line " << mNumberOfLines << " ]";
                KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
            }


        }

        KRATOS_CATCH("")
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
            KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");

        }

		VariableData const& r_x_variable = KratosComponents<VariableData>::Get(variable_name);

        ReadWord(variable_name);

        if(!KratosComponents<VariableData>::Has(variable_name))
        {
            std::stringstream buffer;
            buffer << variable_name << " is not a valid value variable!!! Table only accepts double values." << std::endl;
            buffer << " [Line " << mNumberOfLines << " ]";
            KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");

        }
		VariableData const& r_y_variable = KratosComponents<VariableData>::Get(variable_name);

        while(!mFile.eof())
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
        //    KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");

        //}

        ReadWord(variable_name);

        //if(!KratosComponents<Variable<double> >::Has(variable_name))
        //{
        //    std::stringstream buffer;
        //    buffer << variable_name << " is not a valid value variable!!! Table only accepts double values." << std::endl;
        //    buffer << " [Line " << mNumberOfLines << " ]";
        //    KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");

        //}

        while(!mFile.eof())
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

    void ModelPartIO::ReadNodesBlock(NodesContainerType& rThisNodes)
    {
        KRATOS_THROW_ERROR(std::logic_error,"it should not call this block","")
    }
/*    {
        KRATOS_TRY
//         KRATOS_THROW_ERROR(std::runtime_error, "this function shall be removed ReadNodesBlock(NodesContainerType& rThisNodes)","")

        NodeType temp_node;
        SizeType temp_id;

        std::string word;

        SizeType number_of_nodes_read = 0;

        std::cout << "  [Reading Nodes    : ";

        while(!mFile.eof())
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


            rThisNodes.push_back(temp_node);
            number_of_nodes_read++;
        }
        std::cout << number_of_nodes_read << " nodes read]" << std::endl;

        unsigned int numer_of_nodes_read = rThisNodes.size();
        rThisNodes.Unique();
        if(rThisNodes.size() != numer_of_nodes_read)
            std::cout << "attention! we read " << numer_of_nodes_read << " but there are only " << rThisNodes.size() << " non repeated nodes" << std::endl;

        KRATOS_CATCH("")
    }
*/
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

        std::cout << "  [Reading Nodes    : ";

        while(!mFile.eof())
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
        std::cout << number_of_nodes_read << " nodes read]" << std::endl;

        unsigned int numer_of_nodes_read = rModelPart.Nodes().size();
        rModelPart.Nodes().Unique();
        if(rModelPart.Nodes().size() != numer_of_nodes_read)
            std::cout << "attention! we read " << numer_of_nodes_read << " but there are only " << rModelPart.Nodes().size() << " non repeated nodes" << std::endl;
*/
	SizeType id;
	double x;
	double y;
	double z;

	std::string word;

	SizeType number_of_nodes_read = 0;
        
        typedef std::map< unsigned int, array_1d<double,3> > map_type;
        map_type read_coordinates;

        std::cout << "  [Reading Nodes    : ";

	while(!mFile.eof())
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
         read_coordinates[id] = coords;
	  number_of_nodes_read++;
	}
	
	for(map_type::const_iterator it = read_coordinates.begin(); it!=read_coordinates.end(); ++it)
        {
            const unsigned int node_id = it->first;
            const array_1d<double,3>& coords = it->second;
            rModelPart.CreateNewNode(node_id,coords[0],coords[1],coords[2]);
        }
        if(rModelPart.Nodes().size() != number_of_nodes_read)
            std::cout << "attention! we read " << number_of_nodes_read << " but there are only " << rModelPart.Nodes().size() << " non repeated nodes" << std::endl;

        KRATOS_CATCH("")
    }

    std::size_t ModelPartIO::CountNodesInBlock()
    {
        KRATOS_TRY;

        std::vector<SizeType> found_ids;

        SizeType temp_id;

        std::string word;

        SizeType number_of_nodes_read = 0;

	//std::cout << "  [Reading Nodes    : ";

        while(!mFile.eof())
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
        //std::cout << number_of_nodes_read << " nodes read]" << std::endl;

        // Error check: look for duplicate nodes
        std::sort(found_ids.begin(),found_ids.end());
        std::vector<std::size_t>::iterator unique_end = std::unique(found_ids.begin(),found_ids.end());
        std::size_t number_of_unique_nodes = std::distance(found_ids.begin(),unique_end);

        if(number_of_unique_nodes != number_of_nodes_read)
            std::cout << "attention! we read " << number_of_nodes_read << " but there are only " << number_of_unique_nodes << " non repeated nodes" << std::endl;

        return number_of_nodes_read;

        KRATOS_CATCH("");
    }

    void ModelPartIO::ReadPropertiesBlock(PropertiesContainerType& rThisProperties)
    {
        KRATOS_TRY

        Properties::Pointer props = boost::make_shared<Properties>();
        Properties& temp_properties = *props;
        //Properties temp_properties;
        
        std::string word;
        std::string variable_name;

        SizeType temp_properties_id;

        ReadWord(word);
        ExtractValue(word, temp_properties_id);
        temp_properties.SetId(temp_properties_id);

        while(!mFile.eof())
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
                KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
            }





        }

        rThisProperties.push_back(props);
//         rThisProperties.push_back(temp_properties);

        KRATOS_CATCH("")
    }

    void ModelPartIO::ReadElementsBlock(ModelPart& rModelPart)
    {
        ReadElementsBlock(rModelPart.Nodes(), rModelPart.rProperties(), rModelPart.Elements());
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
        std::cout << "  [Reading Elements : ";

        if(!KratosComponents<Element>::Has(element_name))
        {
            std::stringstream buffer;
            buffer << "Element " << element_name << " is not registered in Kratos.";
            buffer << " Please check the spelling of the element name and see if the application which containing it, is registered corectly.";
            buffer << " [Line " << mNumberOfLines << " ]";
            KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
            return;
        }

        Element const& r_clone_element = KratosComponents<Element>::Get(element_name);
        SizeType number_of_nodes = r_clone_element.GetGeometry().size();
        Element::NodesArrayType temp_element_nodes;

        while(!mFile.eof())
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
        std::cout << number_of_read_elements << " elements read] [Type: " <<element_name << "]" << std::endl;
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
        std::cout << "  [Reading Conditions : ";

        if(!KratosComponents<Condition>::Has(condition_name))
        {
            std::stringstream buffer;
            buffer << "Condition " << condition_name << " is not registered in Kratos.";
            buffer << " Please check the spelling of the condition name and see if the application containing it is registered corectly.";
            buffer << " [Line " << mNumberOfLines << " ]";
            KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
            return;
        }

        Condition const& r_clone_condition = KratosComponents<Condition>::Get(condition_name);
        SizeType number_of_nodes = r_clone_condition.GetGeometry().size();
        Condition::NodesArrayType temp_condition_nodes;

        while(!mFile.eof())
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
        std::cout << number_of_read_conditions << " conditions read] [Type: " << condition_name << "]" << std::endl;
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
        
        VariablesList r_modelpart_nodal_variables_list = rThisModelPart.GetNodalSolutionStepVariablesList();
        

        if(KratosComponents<Flags >::Has(variable_name))
        {
            ReadNodalFlags(rThisNodes, static_cast<Flags const& >(KratosComponents<Flags >::Get(variable_name)));
        }
        else if(KratosComponents<Variable<int> >::Has(variable_name))
        {
            bool has_been_added = r_modelpart_nodal_variables_list.Has(KratosComponents<Variable<int> >::Get(variable_name)) ;
			if( !has_been_added && mOptions.Is(IGNORE_VARIABLES_ERROR) ) {
                std::cout<<std::endl<<"WARNING: Skipping NodalData block. Variable "<<variable_name<<" has not been added to ModelPart '"<<rThisModelPart.Name()<<"'"<<std::endl<<std::endl;
                SkipBlock("NodalData");
            }
			else if (!has_been_added)
				KRATOS_THROW_ERROR(std::invalid_argument,"The nodal solution step container deos not have this variable: ", variable_name)
            else {
                ReadNodalScalarVariableData(rThisNodes, static_cast<Variable<int> const& >(KratosComponents<Variable<int> >::Get(variable_name)));
            }
        }
        else if(KratosComponents<Variable<double> >::Has(variable_name))
        {                    
            bool has_been_added = r_modelpart_nodal_variables_list.Has(KratosComponents<Variable<double> >::Get(variable_name)) ;
            if( !has_been_added && mOptions.Is(IGNORE_VARIABLES_ERROR) ) {
                std::cout<<std::endl<<"WARNING: Skipping NodalData block. Variable "<<variable_name<<" has not been added to ModelPart '"<<rThisModelPart.Name()<<"'"<<std::endl<<std::endl;
                SkipBlock("NodalData");
            }
			else if (!has_been_added)
				KRATOS_THROW_ERROR(std::invalid_argument,"The nodal solution step container deos not have this variable: ", variable_name)
            else {
                ReadNodalDofVariableData(rThisNodes, static_cast<Variable<double> const& >(KratosComponents<Variable<double> >::Get(variable_name)));                
            }
        }
        else if(KratosComponents<array_1d_component_type>::Has(variable_name))
        {
            ReadNodalDofVariableData(rThisNodes, static_cast<array_1d_component_type const& >(KratosComponents<array_1d_component_type>::Get(variable_name)));
        }
        else if(KratosComponents<Variable<array_1d<double, 3> > >::Has(variable_name))
        {         
            bool has_been_added = r_modelpart_nodal_variables_list.Has(KratosComponents<Variable<array_1d<double, 3> > >::Get(variable_name)) ;
            if( !has_been_added && mOptions.Is(IGNORE_VARIABLES_ERROR) ) {
                std::cout<<std::endl<<"WARNING: Skipping NodalData block. Variable "<<variable_name<<" has not been added to ModelPart '"<<rThisModelPart.Name()<<"'"<<std::endl<<std::endl;
            }
			else if (!has_been_added)
				KRATOS_THROW_ERROR(std::invalid_argument,"The nodal solution step container deos not have this variable: ", variable_name)
            else {
                ReadNodalVectorialVariableData(rThisNodes, static_cast<Variable<array_1d<double, 3> > const& >(KratosComponents<Variable<array_1d<double, 3> > >::Get(variable_name)), Vector(3));
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
            KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
        }
        else
        {
            std::stringstream buffer;
            buffer << variable_name << " is not a valid variable!!!" << std::endl;
            buffer << " [Line " << mNumberOfLines << " ]";
            KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
        }

        KRATOS_CATCH("")
    }

    template<class TVariableType>
    void ModelPartIO::ReadNodalDofVariableData(NodesContainerType& rThisNodes, TVariableType& rVariable)
    {
        KRATOS_TRY

        SizeType id;
        bool is_fixed;
        double nodal_value;

        std::string value;


        while(!mFile.eof())
        {
            ReadWord(value); // reading id
            if(CheckEndBlock("NodalData", value))
                break;

            ExtractValue(value, id);
            typename NodesContainerType::iterator i_node = FindKey(rThisNodes, ReorderedNodeId(id), "Node");

            // reading is_fixed
            ReadWord(value);
            ExtractValue(value, is_fixed);
            if(is_fixed)
                i_node->Fix(rVariable);

            // reading nodal_value
            ReadWord(value);
            ExtractValue(value, nodal_value);

            i_node->GetSolutionStepValue(rVariable, 0) =  nodal_value;
        }

        KRATOS_CATCH("")
    }


    void ModelPartIO::ReadNodalFlags(NodesContainerType& rThisNodes, Flags const& rFlags)
    {

        KRATOS_TRY

        SizeType id;

        std::string value;

        while(!mFile.eof())
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

        while(!mFile.eof())
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
                KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
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

        while(!mFile.eof())
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
                KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
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
        else if(KratosComponents<Variable<Matrix> >::Has(variable_name))
        {
            ReadElementalVectorialVariableData(rThisElements, static_cast<Variable<Matrix > const& >(KratosComponents<Variable<Matrix> >::Get(variable_name)), Matrix(3,3));
        }
        else
        {
            std::stringstream buffer;
            buffer << variable_name << " is not a valid variable!!!" << std::endl;
            buffer << " [Line " << mNumberOfLines << " ]";
            KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
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

        while(!mFile.eof())
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
                std::cout  << "WARNING! Assigning " << rVariable.Name() << " to not existing element #" << id << " [Line " << mNumberOfLines << " ]" << std::endl;
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

        while(!mFile.eof())
        {
            ReadWord(value); // reading id
            if(CheckEndBlock("ElementalData", value))
                break;

            ExtractValue(value, id);



            // reading nodal_value
            ReadVectorialValue(elemental_value);
            ExtractValue(value, elemental_value);

            ModelPart::ElementIterator i_result = rThisElements.find(ReorderedElementId(id));
            if(i_result != rThisElements.end())
                i_result->GetValue(rVariable) =  elemental_value;
            else
                std::cout  << "WARNING! Assigning " << rVariable.Name() << " to not existing element #" << id << " [Line " << mNumberOfLines << " ]" << std::endl;
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
        else if(KratosComponents<Variable<Matrix> >::Has(variable_name))
        {
            ReadConditionalVectorialVariableData(rThisConditions, static_cast<Variable<Matrix > const& >(KratosComponents<Variable<Matrix> >::Get(variable_name)), Matrix(3,3));
        }
        else
        {
            std::stringstream buffer;
            buffer << variable_name << " is not a valid variable!!!" << std::endl;
            buffer << " [Line " << mNumberOfLines << " ]";
            KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
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

        while(!mFile.eof())
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
                std::cout  << "WARNING! Assigning " << rVariable.Name() << " to not existing condition #" << id << " [Line " << mNumberOfLines << " ]" << std::endl;
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

        while(!mFile.eof())
        {
            ReadWord(value); // reading id
            if(CheckEndBlock("ConditionalData", value))
                break;

            ExtractValue(value, id);

            // reading nodal_value
            ReadVectorialValue(conditional_value);
            ExtractValue(value, conditional_value);

            ModelPart::ConditionIterator i_result = rThisConditions.find(ReorderedConditionId(id));
            if(i_result != rThisConditions.end())
                i_result->GetValue(rVariable) =  conditional_value;
            else
                std::cout  << "WARNING! Assigning " << rVariable.Name() << " to not existing condition #" << id << " [Line " << mNumberOfLines << " ]" << std::endl;
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
            KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
            return number_of_connectivities;
        }

        Element const& r_clone_element = KratosComponents<Element>::Get(element_name);
        SizeType number_of_nodes = r_clone_element.GetGeometry().size();
        ConnectivitiesContainerType::value_type temp_element_nodes;

        while(!mFile.eof())
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
            KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
            return number_of_connectivities;
        }

        Condition const& r_clone_condition = KratosComponents<Condition>::Get(condition_name);
        SizeType number_of_nodes = r_clone_condition.GetGeometry().size();
        ConnectivitiesContainerType::value_type temp_condition_nodes;

        while(!mFile.eof())
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
            KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
        }

        Element const& r_clone_element = KratosComponents<Element>::Get(element_name);
        SizeType n_nodes_in_elem = r_clone_element.GetGeometry().size();
        ConnectivitiesContainerType::value_type temp_element_nodes;

        while(!mFile.eof())
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
            KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
        }

        Condition const& r_clone_condition = KratosComponents<Condition>::Get(condition_name);
        SizeType n_nodes_in_cond = r_clone_condition.GetGeometry().size();
        ConnectivitiesContainerType::value_type temp_condition_nodes;

        while(!mFile.eof())
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

// 	KRATOS_WATCH("begin reading CommunicatorDataBlock")

        std::string word;
        while(true)
        {
            ReadWord(word);
            if(mFile.eof())
                break;
            if(CheckEndBlock("CommunicatorData", word))
                break;
            if(word == "NEIGHBOURS_INDICES")
                ReadVectorialValue(rThisCommunicator.NeighbourIndices());
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
                    ReadCommunicatorLocalNodesBlock(rThisCommunicator, rThisNodes);
                else if(word == "GhostNodes")
                    ReadCommunicatorGhostNodesBlock(rThisCommunicator, rThisNodes);
                else
                    SkipBlock(word);
            }
        }

// 	KRATOS_WATCH("finished reading CommunicatorDataBlock")

        return ;

        KRATOS_CATCH("")
    }

    void ModelPartIO::ReadCommunicatorLocalNodesBlock(Communicator& rThisCommunicator, NodesContainerType& rThisNodes)
    {
        KRATOS_TRY

// 	KRATOS_WATCH("begin reading CommunicatorLocalNodesBlock")

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
            KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
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

        while(!mFile.eof())
        {
            ReadWord(word); // Reading the node id or End
            if(CheckEndBlock("LocalNodes", word))
                break;

            ExtractValue(word,node_id);
            NodesContainerType::iterator i_node = FindKey(rThisNodes, ReorderedNodeId(node_id), "Node");
            p_local_mesh->Nodes().push_back(*(i_node.base()));
            p_interface_mesh->Nodes().push_back(*(i_node.base()));
        }

        p_local_mesh->Nodes().Sort();
        p_interface_mesh->Nodes().Sort();
// 	KRATOS_WATCH("finished reading CommunicatorLocalNodesBlock")

// 	KRATOS_WATCH(rThisCommunicator)
        KRATOS_CATCH("")
    }


    void ModelPartIO::ReadCommunicatorGhostNodesBlock(Communicator& rThisCommunicator, NodesContainerType& rThisNodes)
    {
        KRATOS_TRY

// 	KRATOS_WATCH("begin reading CommunicatorGhostNodesBlock")


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
            KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
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

        while(!mFile.eof())
        {
            ReadWord(word); // Reading the node id or End
            if(CheckEndBlock("GhostNodes", word))
                break;

            ExtractValue(word,node_id);
            NodesContainerType::iterator i_node = FindKey(rThisNodes, ReorderedNodeId(node_id), "Node");
            p_ghost_mesh->Nodes().push_back(*(i_node.base()));
            p_interface_mesh->Nodes().push_back(*(i_node.base()));
        }

        p_ghost_mesh->Nodes().Sort();
        p_interface_mesh->Nodes().Sort();


// 	KRATOS_WATCH(rThisCommunicator)
        KRATOS_CATCH("")

// 	KRATOS_WATCH("finished reading CommunicatorGhostNodesBlock")

    }

    void ModelPartIO::ReadMeshBlock(ModelPart& rModelPart)
    {
         KRATOS_TRY

        std::string word;
        SizeType mesh_id;

        ReadWord(word);
        ExtractValue(word, mesh_id);

        

        SizeType number_of_meshes = rModelPart.NumberOfMeshes();

        if(mesh_id > 1000000) // this would be a case of error in reading.
            KRATOS_THROW_ERROR(std::invalid_argument, "Too large mesh id :", mesh_id);

        if(mesh_id == 0) // this would be a case of error in reading.
            KRATOS_THROW_ERROR(std::invalid_argument, "The mesh zero is the reference mesh and already created. You cannot create a mesh 0 with mesh block.", "");

        // adding necessary meshes to the model part.
        MeshType empty_mesh;
        for(SizeType i = number_of_meshes ; i < mesh_id + 1 ; i++)
            rModelPart.GetMeshes().push_back(boost::make_shared<MeshType>(empty_mesh.Clone()));

        MeshType& mesh = rModelPart.GetMesh(mesh_id);

        while(true)
        {
            ReadWord(word);
            
            if(mFile.eof())
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

        while(!mFile.eof())
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
                KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
            }


        }

        KRATOS_CATCH("")
    }


    void ModelPartIO::ReadMeshNodesBlock(ModelPart& rModelPart, MeshType& rMesh)
    {
        KRATOS_TRY

        SizeType node_id;

        std::string word;


        while(!mFile.eof())
        {
            ReadWord(word); // Reading the node id or End
            if(CheckEndBlock("MeshNodes", word))
                break;

            ExtractValue(word,node_id);
            NodesContainerType::iterator i_node = FindKey(rModelPart.Nodes(), ReorderedNodeId(node_id), "Node");
            rMesh.Nodes().push_back(*(i_node.base()));
        }

        rMesh.Nodes().Sort();
        KRATOS_CATCH("")
    }

    void ModelPartIO::ReadMeshElementsBlock(ModelPart& rModelPart, MeshType& rMesh)
    {
        KRATOS_TRY

        SizeType element_id;

        std::string word;


        while(!mFile.eof())
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


        while(!mFile.eof())
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

        Properties::Pointer props = boost::make_shared<Properties>();
        Properties& temp_properties = *props;
//         Properties temp_properties;

        std::string word;
        std::string variable_name;

        SizeType temp_properties_id;

        ReadWord(word);
        ExtractValue(word, temp_properties_id);
        temp_properties.SetId(temp_properties_id);


        while(!mFile.eof())
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
                KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
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
			//if (mFile.eof())
			//	break;
			if (CheckEndBlock("SubModelPart", word))
				break;

			ReadBlockName(word);
			if (word == "SubModelPartData")
				ReadModelPartDataBlock(r_sub_model_part);
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
			//	ReadCommunicatorDataBlock(rThisModelPart.GetCommunicator(), rThisModelPart.Nodes());
			//	//Adding the elements and conditions to the communicator
			//	rThisModelPart.GetCommunicator().LocalMesh().Elements() = rThisModelPart.Elements();
			//	rThisModelPart.GetCommunicator().LocalMesh().Conditions() = rThisModelPart.Conditions();
			//}
			//else if (word == "Mesh")
			//	ReadMeshBlock(rThisModelPart);
			else if (word == "SubModelPart")
				ReadSubModelPartBlock(rMainModelPart, r_sub_model_part);
		}

		KRATOS_CATCH("")
	}
	void  ModelPartIO::ReadSubModelPartTablesBlock(ModelPart& rMainModelPart, ModelPart& rSubModelPart)
	{
		KRATOS_TRY

		SizeType table_id;
		std::string word;

		while (!mFile.eof())
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

		while (!mFile.eof())
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

		while (!mFile.eof())
		{
			ReadWord(word); // Reading the node id or End
			if (CheckEndBlock("SubModelPartNodes", word))
				break;

			ExtractValue(word, node_id);
			NodesContainerType::iterator i_node = FindKey(rMainModelPart.Nodes(), ReorderedNodeId(node_id), "Node");
			rSubModelPart.AddNode(*(i_node.base()));
		}
		KRATOS_CATCH("")
	}

	void  ModelPartIO::ReadSubModelPartElementsBlock(ModelPart& rMainModelPart, ModelPart& rSubModelPart)
	{
		KRATOS_TRY

		SizeType element_id;
		std::string word;

		while (!mFile.eof())
		{
			ReadWord(word); // Reading the node id or End
			if (CheckEndBlock("SubModelPartElements", word))
				break;

			ExtractValue(word, element_id);
			ElementsContainerType::iterator i_element = FindKey(rMainModelPart.Elements(), ReorderedElementId(element_id), "Element");
			rSubModelPart.AddElement(*(i_element.base()));
		}
		KRATOS_CATCH("")
	}

	void  ModelPartIO::ReadSubModelPartConditionsBlock(ModelPart& rMainModelPart, ModelPart& rSubModelPart)
	{
		KRATOS_TRY

		SizeType condition_id;
		std::string word;

		while (!mFile.eof())
		{
			ReadWord(word); // Reading the node id or End
			if (CheckEndBlock("SubModelPartConditions", word))
				break;

			ExtractValue(word, condition_id);
			ConditionsContainerType::iterator i_condition = FindKey(rMainModelPart.Conditions(), ReorderedConditionId(condition_id), "Condition");
			rSubModelPart.AddCondition(*(i_condition.base()));
		}
		KRATOS_CATCH("")
	}


    void ModelPartIO::DivideModelPartDataBlock(OutputFilesContainerType& OutputFiles)
    {
        KRATOS_TRY
        std::string block;

        WriteInAllFiles(OutputFiles, "Begin ModelPartData");

        ReadBlock(block, "ModelPartData");
        WriteInAllFiles(OutputFiles, block);

        WriteInAllFiles(OutputFiles, "End ModelPartData\n");
        KRATOS_WATCH("DivideModelPartDataBlock completed");
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
        KRATOS_WATCH("DivideTableBlock completed");
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
        KRATOS_WATCH("DividePropertiesBlock completed");
        KRATOS_CATCH("")
    }

    void ModelPartIO::DivideNodesBlock(OutputFilesContainerType& OutputFiles,
                          PartitionIndicesContainerType const& NodesAllPartitions)
    {
        KRATOS_TRY

        std::string word;

        WriteInAllFiles(OutputFiles, "Begin Nodes \n");

        SizeType id;

        while(!mFile.eof())
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
                KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
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
                    KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
                }

                *(OutputFiles[partition_id]) << node_data.str();
            }

        }

        WriteInAllFiles(OutputFiles, "End Nodes\n");
        KRATOS_WATCH("DivideNodesBlock completed");

        KRATOS_CATCH("")
    }

    void ModelPartIO::DivideElementsBlock(OutputFilesContainerType& OutputFiles,
                             PartitionIndicesContainerType const& ElementsAllPartitions)
    {
        KRATOS_TRY

        KRATOS_WATCH("DivideElementsBlock started");

        std::string word;
        std::string element_name;

        ReadWord(element_name);
        if(!KratosComponents<Element>::Has(element_name))
        {
            std::stringstream buffer;
            buffer << "Element " << element_name << " is not registered in Kratos.";
            buffer << " Please check the spelling of the element name and see if the application containing it is registered corectly.";
            buffer << " [Line " << mNumberOfLines << " ]";
            KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
            return;
        }

        Element const& r_clone_element = KratosComponents<Element>::Get(element_name);
        SizeType number_of_nodes = r_clone_element.GetGeometry().size();

        WriteInAllFiles(OutputFiles, "Begin Elements " +  element_name);

        SizeType id;

        while(!mFile.eof())
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
                KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
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
                    KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
                }

                *(OutputFiles[partition_id]) << element_data.str();
            }

        }

        WriteInAllFiles(OutputFiles, "\nEnd Elements\n");

        KRATOS_WATCH("DivideElementsBlock completed");

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
            KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
            return;
        }

        Condition const& r_clone_condition = KratosComponents<Condition>::Get(condition_name);
        SizeType number_of_nodes = r_clone_condition.GetGeometry().size();

        WriteInAllFiles(OutputFiles, "Begin Conditions " +  condition_name);

        SizeType id;

        while(!mFile.eof())
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
                KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
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
                    KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
                }

                *(OutputFiles[partition_id]) << condition_data.str();
            }

        }

        WriteInAllFiles(OutputFiles, "\nEnd Conditions\n");

        KRATOS_WATCH("DivideConditionsBlock completed");

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
        else if(KratosComponents<Variable<Matrix> >::Has(variable_name))
        {
            DivideVectorialVariableData(OutputFiles, NodesAllPartitions, "NodalData" );
        }
        else if(KratosComponents<VariableData>::Has(variable_name))
        {
            std::stringstream buffer;
            buffer << variable_name << " is not supported to be read by this IO or the type of variable is not registered correctly" << std::endl;
            buffer << " [Line " << mNumberOfLines << " ]";
            KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
        }
        else
        {
            std::stringstream buffer;
            buffer << variable_name << " is not a valid variable!!!" << std::endl;
            buffer << " [Line " << mNumberOfLines << " ]";
            KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
        }

        WriteInAllFiles(OutputFiles, "End NodalData\n");

        KRATOS_WATCH("DivideNodalDataBlock completed");

        KRATOS_CATCH("")
    }

    void ModelPartIO::DivideDofVariableData(OutputFilesContainerType& OutputFiles,
                               PartitionIndicesContainerType const& NodesAllPartitions)
    {
        KRATOS_TRY

        SizeType id;

        std::string word;


        while(!mFile.eof())
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
                KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
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
                    KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
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

        while(!mFile.eof())
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
                KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
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
                    KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
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
                    KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
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
        else if(KratosComponents<Variable<Matrix> >::Has(variable_name))
        {
            DivideVectorialVariableData(OutputFiles, ElementsAllPartitions, "ElementalData");
        }
        else if(KratosComponents<VariableData>::Has(variable_name))
        {
            std::stringstream buffer;
            buffer << variable_name << " is not supported to be read by this IO or the type of variable is not registered correctly" << std::endl;
            buffer << " [Line " << mNumberOfLines << " ]";
            KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
        }
        else
        {
            std::stringstream buffer;
            buffer << variable_name << " is not a valid variable!!!" << std::endl;
            buffer << " [Line " << mNumberOfLines << " ]";
            KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
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


        while(!mFile.eof())
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
				KRATOS_THROW_ERROR(std::logic_error, "Invalid block name :", BlockName);

            if(index > EntitiesPartitions.size())
            {
                std::stringstream buffer;
                buffer << "Invalid id : " << id;
                buffer << " [Line " << mNumberOfLines << " ]";
                KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
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
                    KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
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
        else if(KratosComponents<Variable<Matrix> >::Has(variable_name))
        {
            DivideVectorialVariableData(OutputFiles, ConditionsAllPartitions, "ConditionalData");
        }
        else if(KratosComponents<VariableData>::Has(variable_name))
        {
            std::stringstream buffer;
            buffer << variable_name << " is not supported to be read by this IO or the type of variable is not registered correctly" << std::endl;
            buffer << " [Line " << mNumberOfLines << " ]";
            KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
        }
        else
        {
            std::stringstream buffer;
            buffer << variable_name << " is not a valid variable!!!" << std::endl;
            buffer << " [Line " << mNumberOfLines << " ]";
            KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
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

        while(!mFile.eof())
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
		KRATOS_ERROR << "This Method is not implemented yet!" << std::endl;
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

        while(!mFile.eof())
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
                KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
            }

            for(SizeType i = 0 ; i < NodesAllPartitions[ReorderedNodeId(id)-1].size() ; i++)
            {
                SizeType partition_id = NodesAllPartitions[ReorderedNodeId(id)-1][i];
                if(partition_id > OutputFiles.size())
                {
                    std::stringstream buffer;
                    buffer << "Invalid prtition id : " << partition_id;
                    buffer << " for node " << id << " [Line " << mNumberOfLines << " ]";
                    KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
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

        while(!mFile.eof())
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
                KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
            }

            for(SizeType i = 0 ; i < ElementsAllPartitions[ReorderedElementId(id)-1].size() ; i++)
            {
                SizeType partition_id = ElementsAllPartitions[ReorderedElementId(id)-1][i];
                if(partition_id > OutputFiles.size())
                {
                    std::stringstream buffer;
                    buffer << "Invalid prtition id : " << partition_id;
                    buffer << " for element " << id << " [Line " << mNumberOfLines << " ]";
                    KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
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

        while(!mFile.eof())
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
                KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
            }

            for(SizeType i = 0 ; i < ConditionsAllPartitions[ReorderedConditionId(id)-1].size() ; i++)
            {
                SizeType partition_id = ConditionsAllPartitions[ReorderedConditionId(id)-1][i];
                if(partition_id > OutputFiles.size())
                {
                    std::stringstream buffer;
                    buffer << "Invalid prtition id : " << partition_id;
                    buffer << " for condition " << id << " [Line " << mNumberOfLines << " ]";
                    KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
                }

                *(OutputFiles[partition_id]) << ReorderedConditionId(id) << std::endl;
            }

        }

        WriteInAllFiles(OutputFiles, "End MeshConditions\n");

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
                    KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
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
            boost::numeric::ublas::vector<int> indices = row(DomainsColoredGraph, i_partition);
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
            boost::numeric::ublas::vector<int> indices = row(DomainsColoredGraph, i_partition);
            *(OutputFiles[i_partition]) << "NUMBER_OF_COLORS    " << number_of_colors << std::endl;
        }


        // Writing the all local nodes
        WriteInAllFiles(OutputFiles, "    Begin LocalNodes 0\n");

        for(SizeType i = 0 ; i < NodesPartitions.size() ; i++)
            *(OutputFiles[NodesPartitions[i]]) << "    " << i+1 << std::endl;

        WriteInAllFiles(OutputFiles, "    End LocalNodes \n");


        std::vector<PartitionIndicesContainerType> local_nodes_indices(NumberOfPartitions, PartitionIndicesContainerType(number_of_colors));
        std::vector<PartitionIndicesContainerType> ghost_nodes_indices(NumberOfPartitions, PartitionIndicesContainerType(number_of_colors));

        matrix<int> interface_indices = scalar_matrix<int>(NumberOfPartitions, NumberOfPartitions, -1);

        for(SizeType i_partition = 0 ; i_partition < NumberOfPartitions ; i_partition++)
        {
            vector<int> neighbours_indices = row(DomainsColoredGraph, i_partition);

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
            KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
        }

        return i_result;
    }




    // Basically it starts to read the character sequence until reaching a
    // "(" and then goes until corresponding ")" which means the vector or
    // matrix value is completely read. It can be used to read any kind of
    // vector or matrix with operator >> defined and writtern in following
    // format for a vector: [size] ( value1, value2,...., valueN )
    // format for a matrix: [size1,size2] ( )( )...( ) //look props read
    template<class TValueType>
    TValueType& ModelPartIO::ReadVectorialValue(TValueType& rValue)
    {
        std::stringstream value;

        char c = SkipWhiteSpaces();
        while((c != '(') && !mFile.eof())
        {
            value << c;
            c = GetCharacter();
        }
        int open_parantesis = 1;
        while((open_parantesis != 0) && !mFile.eof())
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


    ModelPartIO& ModelPartIO::ReadWord(std::string& Word)
    {
        Word.clear();

        char c = SkipWhiteSpaces();
        while(!mFile.eof() && !IsWhiteSpace(c))
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

        while(!mFile.eof())
        {
            if(c == 'E')
            {
                word.clear();
                while(!mFile.eof() && !IsWhiteSpace(c))
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
				while (!mFile.eof() && !IsWhiteSpace(c))
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
        if(mFile.get(c))
        {
            if(c == '\n')
                mNumberOfLines++;
            else if(c == '/') // it may be a comment!
            {
                char next_c = mFile.peek();
                if(next_c == '/') // it's a line comment
                {
                    while((mFile.get(c)) && (c != '\n')); // so going to the end of line
                    if(!mFile.eof())
                        mNumberOfLines++;
                }
                else if(next_c == '*') // it's a block comment
                {
                    while((mFile.get(c)) && !((c == '*') && (mFile.peek() == '/'))) // so going to the end of block
                        if(c == '\n')
                            mNumberOfLines++;
                    mFile.get(c);
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
            KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
        }
        else
            result = true;

        return result;
    }

    void ModelPartIO::ResetInput()
    {
        mFile.clear();
        mFile.seekg(0, std::ios_base::beg);
        mNumberOfLines = 1;
    }

    inline void ModelPartIO::CreatePartition(unsigned int number_of_threads,const int number_of_rows, vector<unsigned int>& partitions)
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

