// Kratos Multi-Physics
//
// Copyright (c) 2016 Pooyan Dadvand, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
// 	-	Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// 	-	Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
// 		in the documentation and/or other materials provided with the distribution.
// 	-	All advertising materials mentioning features or use of this software must display the following acknowledgement:
// 			This product includes Kratos Multi-Physics technology.
// 	-	Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED ANDON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
// THE USE OF THISSOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.




#if !defined(KRATOS_MODEL_PART_IO_H_INCLUDED )
#define  KRATOS_MODEL_PART_IO_H_INCLUDED



// System includes
#include <string>
#include <fstream>
#include <set>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/io.h"
#include "utilities/timer.h"
#include "containers/flags.h"


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
class ModelPartIO : public IO
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ModelPartIO
    KRATOS_CLASS_POINTER_DEFINITION(ModelPartIO);

   typedef IO BaseType;

    typedef BaseType::NodeType NodeType;

    typedef BaseType::MeshType MeshType;

    typedef BaseType::NodesContainerType NodesContainerType;

    typedef BaseType::PropertiesContainerType PropertiesContainerType;

    typedef BaseType::ElementsContainerType ElementsContainerType;

    typedef BaseType::ConditionsContainerType ConditionsContainerType;

    typedef BaseType::ConnectivitiesContainerType ConnectivitiesContainerType;

    typedef std::vector<std::ofstream*> OutputFilesContainerType;

    typedef std::size_t SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with  filenames.
	ModelPartIO(std::string const& Filename, const Flags Options = IO::READ|IO::NOT_IGNORE_VARIABLES_ERROR);


    /// Constructor with filenames.
//       ModelPartIO(std::string const& InputFilename, std::string const& OutputFilename)
// 	: mNumberOfLines(0), mInput(std::ifstream(InputFilename.c_str())), mOutput(std::ofstream(OutputFilename.c_str()))
//       {
//       }


    /// Destructor.
    virtual ~ModelPartIO();


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    virtual bool ReadNode(NodeType& rThisNode);

    virtual bool ReadNodes(NodesContainerType& rThisNodes);

    virtual std::size_t ReadNodesNumber();

    virtual void WriteNodes(NodesContainerType const& rThisNodes);

    virtual void ReadProperties(Properties& rThisProperties);

    virtual void ReadProperties(PropertiesContainerType& rThisProperties);

    virtual void WriteProperties(PropertiesContainerType& rThisProperties);

    virtual void ReadElement(NodesContainerType& rThisNodes, PropertiesContainerType& rThisProperties, Element::Pointer& pThisElements);

    virtual void ReadElements(NodesContainerType& rThisNodes, PropertiesContainerType& rThisProperties, ElementsContainerType& rThisElements);

    virtual std::size_t  ReadElementsConnectivities(ConnectivitiesContainerType& rElementsConnectivities);

    virtual void WriteElements(ElementsContainerType const& rThisElements);

    virtual void ReadConditions(NodesContainerType& rThisNodes, PropertiesContainerType& rThisProperties, ConditionsContainerType& rThisConditions);

    virtual std::size_t  ReadConditionsConnectivities(ConnectivitiesContainerType& rConditionsConnectivities);

    virtual void WriteConditions(ConditionsContainerType const& rThisConditions);

    virtual void ReadInitialValues(ModelPart& rThisModelPart);

//       void ReadGeometries(NodesContainerType& rThisNodes, GeometriesContainerType& rResults);

    virtual void ReadMesh(MeshType & rThisMesh);

    virtual void WriteMesh(MeshType & rThisMesh);

    virtual void ReadModelPart(ModelPart & rThisModelPart);



    virtual void WriteModelPart(ModelPart & rThisModelPart);


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
    virtual std::size_t ReadNodalGraph(int **NodeIndices, int **NodeConnectivities);

    virtual void DivideInputToPartitions(SizeType NumberOfPartitions, GraphType const& DomainsColoredGraph,
                                         PartitionIndicesType const& NodesPartitions,
                                         PartitionIndicesType const& ElementsPartitions,
                                         PartitionIndicesType const& ConditionsPartitions,
                                         PartitionIndicesContainerType const& NodesAllPartitions,
                                         PartitionIndicesContainerType const& ElementsAllPartitions,
                                         PartitionIndicesContainerType const& ConditionsAllPartitions);


    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

//       /// Turn back information as a string.
//       virtual std::string Info() const
//  {
//    return "ModelPartIO";
//  }

//       /// Print information about this object.
//       virtual void PrintInfo(std::ostream& rOStream) const
//  {
//    rOStream << "ModelPartIO";
//  }


//       /// Print object's data.
//       virtual void PrintData(std::ostream& rOStream) const
//  {
//  }


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

	virtual ModelPartIO::SizeType ReorderedNodeId(ModelPartIO::SizeType NodeId);
	virtual ModelPartIO::SizeType ReorderedElementId(ModelPartIO::SizeType ElementId);
	virtual ModelPartIO::SizeType ReorderedConditionId(ModelPartIO::SizeType ConditionId);

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

    SizeType mNumberOfLines;
    std::string mBaseFilename;
    std::string mFilename;
    std::fstream mFile;
    Flags mOptions;


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    std::string& ReadBlockName(std::string& rBlockName);


    void SkipBlock(std::string const& BlockName);

    bool CheckEndBlock(std::string const& BlockName, std::string& rWord);

    void ReadModelPartDataBlock(ModelPart& rModelPart);

	template<class TablesContainerType>
    void ReadTableBlock(TablesContainerType& rTables);

    void ReadTableBlock(ModelPart::TablesContainerType& rTables);

    void ReadNodesBlock(NodesContainerType& rThisNodes);

    void ReadNodesBlock(ModelPart& rModelPart);

    std::size_t CountNodesInBlock();

    void ReadPropertiesBlock(PropertiesContainerType& rThisProperties);

    void ReadElementsBlock(ModelPart& rModelPart);

    void ReadElementsBlock(NodesContainerType& rThisNodes, PropertiesContainerType& rThisProperties, ElementsContainerType& rThisElements);


    void ReadConditionsBlock(ModelPart& rModelPart);

    void ReadConditionsBlock(NodesContainerType& rThisNodes, PropertiesContainerType& rThisProperties, ConditionsContainerType& rThisConditions);


    void ReadNodalDataBlock(ModelPart& rThisModelPart);

    template<class TVariableType>
    void ReadNodalDofVariableData(NodesContainerType& rThisNodes, TVariableType& rVariable);


    void ReadNodalFlags(NodesContainerType& rThisNodes, Flags const& rFlags);

    template<class TVariableType>
    void ReadNodalScalarVariableData(NodesContainerType& rThisNodes, TVariableType& rVariable);



    template<class TVariableType, class TDataType>
    void ReadNodalVectorialVariableData(NodesContainerType& rThisNodes, TVariableType& rVariable, TDataType Dummy);

    void ReadElementalDataBlock(ElementsContainerType& rThisElements);

    template<class TVariableType>
    void ReadElementalScalarVariableData(ElementsContainerType& rThisElements, TVariableType& rVariable);


    template<class TVariableType, class TDataType>
    void ReadElementalVectorialVariableData(ElementsContainerType& rThisElements, TVariableType& rVariable, TDataType Dummy);
    void ReadConditionalDataBlock(ConditionsContainerType& rThisConditions);

    template<class TVariableType>
    void ReadConditionalScalarVariableData(ConditionsContainerType& rThisConditions, TVariableType& rVariable);


    template<class TVariableType, class TDataType>
    void ReadConditionalVectorialVariableData(ConditionsContainerType& rThisConditions, TVariableType& rVariable, TDataType Dummy);


    SizeType ReadElementsConnectivitiesBlock(ConnectivitiesContainerType& rThisConnectivities);


    SizeType ReadConditionsConnectivitiesBlock(ConnectivitiesContainerType& rThisConnectivities);

    void FillNodalConnectivitiesFromElementBlock(ConnectivitiesContainerType& rNodalConnectivities);

    void FillNodalConnectivitiesFromConditionBlock(ConnectivitiesContainerType& rNodalConnectivities);


    void ReadCommunicatorDataBlock(Communicator& rThisCommunicator, NodesContainerType& rThisNodes);

    void ReadCommunicatorLocalNodesBlock(Communicator& rThisCommunicator, NodesContainerType& rThisNodes);


    void ReadCommunicatorGhostNodesBlock(Communicator& rThisCommunicator, NodesContainerType& rThisNodes);

    void ReadMeshBlock(ModelPart& rModelPart);

	
    void ReadMeshDataBlock(MeshType& rMesh);


    void ReadMeshNodesBlock(ModelPart& rModelPart, MeshType& rMesh);

    void ReadMeshElementsBlock(ModelPart& rModelPart, MeshType& rMesh);

    void ReadMeshConditionsBlock(ModelPart& rModelPart, MeshType& rMesh);
    
    void ReadMeshPropertiesBlock(ModelPart& rModelPart, MeshType& rMesh);

    void DivideModelPartDataBlock(OutputFilesContainerType& OutputFiles);

    void DivideTableBlock(OutputFilesContainerType& OutputFiles);

    void DividePropertiesBlock(OutputFilesContainerType& OutputFiles);

    void DivideNodesBlock(OutputFilesContainerType& OutputFiles,
                          PartitionIndicesContainerType const& NodesAllPartitions);

    void DivideElementsBlock(OutputFilesContainerType& OutputFiles,
                             PartitionIndicesContainerType const& ElementsAllPartitions);



    void DivideConditionsBlock(OutputFilesContainerType& OutputFiles,
                               PartitionIndicesContainerType const& ConditionsAllPartitions);


    void DivideNodalDataBlock(OutputFilesContainerType& OutputFiles,
                              PartitionIndicesContainerType const& NodesAllPartitions);

    void DivideDofVariableData(OutputFilesContainerType& OutputFiles,
                               PartitionIndicesContainerType const& NodesAllPartitions);

    void DivideVectorialVariableData(OutputFilesContainerType& OutputFiles,
                                     PartitionIndicesContainerType const& EntitiesPartitions,
                                     std::string BlockName);


    void DivideElementalDataBlock(OutputFilesContainerType& OutputFiles,
                                  PartitionIndicesContainerType const& ElementsAllPartitions);

    void DivideScalarVariableData(OutputFilesContainerType& OutputFiles,
                                  PartitionIndicesContainerType const& EntitiesPartitions,
                                  std::string BlockName);


    void DivideConditionalDataBlock(OutputFilesContainerType& OutputFiles,
                                    PartitionIndicesContainerType const& ConditionsAllPartitions);


    void DivideMeshBlock(OutputFilesContainerType& OutputFiles,
                                         PartitionIndicesContainerType const& NodesAllPartitions,
                                         PartitionIndicesContainerType const& ElementsAllPartitions,
                                         PartitionIndicesContainerType const& ConditionsAllPartitions);

    void DivideMeshDataBlock(OutputFilesContainerType& OutputFiles);


    void DivideMeshNodesBlock(OutputFilesContainerType& OutputFiles,
                                         PartitionIndicesContainerType const& NodesAllPartitions);


    void DivideMeshElementsBlock(OutputFilesContainerType& OutputFiles,
                                         PartitionIndicesContainerType const& ElementsAllPartitions);

    void DivideMeshConditionsBlock(OutputFilesContainerType& OutputFiles,
                                         PartitionIndicesContainerType const& ConditionsAllPartitions);


    void WritePartitionIndices(OutputFilesContainerType& OutputFiles, PartitionIndicesType const&  NodesPartitions, PartitionIndicesContainerType const& NodesAllPartitions);


    void WriteCommunicatorData(OutputFilesContainerType& OutputFiles, SizeType NumberOfPartitions, GraphType const& DomainsColoredGraph,
                               PartitionIndicesType const& NodesPartitions,
                               PartitionIndicesType const& ElementsPartitions,
                               PartitionIndicesType const& ConditionsPartitions,
                               PartitionIndicesContainerType const& NodesAllPartitions,
                               PartitionIndicesContainerType const& ElementsAllPartitions,
                               PartitionIndicesContainerType const& ConditionsAllPartitions);

    void WriteCommunicatorLocalNodes(OutputFilesContainerType& OutputFiles, SizeType NumberOfPartitions, PartitionIndicesType const& NodesPartitions, PartitionIndicesContainerType const& NodesAllPartitions);

    void WriteInAllFiles(OutputFilesContainerType& OutputFiles, std::string const& ThisWord);


    template<class TContainerType, class TKeyType>
    typename TContainerType::iterator FindKey(TContainerType& ThisContainer , TKeyType ThisKey, std::string ComponentName);




    // Basically it starts to read the character sequence until reaching a
    // "(" and then goes until corresponding ")" which means the vector or
    // matrix value is completely read. It can be used to read any kind of
    // vector or matrix with operator >> defined and writtern in following
    // format for a vector: [size] ( value1, value2,...., valueN )
    // format for a matrix: [size1,size2] ( )( )...( ) //look props read
    template<class TValueType>
    TValueType& ReadVectorialValue(TValueType& rValue);

    template<class TValueType>
    TValueType& ExtractValue(std::string rWord, TValueType & rValue);


    ModelPartIO& ReadWord(std::string& Word);

    ModelPartIO& ReadBlock(std::string& Block, std::string const& BlockName);


    char SkipWhiteSpaces();

    bool IsWhiteSpace(char C);

    char GetCharacter();

    bool CheckStatement(std::string const& rStatement, std::string const& rGivenWord);

    void ResetInput();

    inline void CreatePartition(unsigned int number_of_threads,const int number_of_rows, vector<unsigned int>& partitions);



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
    ModelPartIO& operator=(ModelPartIO const& rOther);

    /// Copy constructor.
    ModelPartIO(ModelPartIO const& rOther);


    ///@}

}; // Class ModelPartIO

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


//   /// input stream function
//   inline std::istream& operator >> (std::istream& rIStream,
//                  ModelPartIO& rThis);

//   /// output stream function
//   inline std::ostream& operator << (std::ostream& rOStream,
//                  const ModelPartIO& rThis)
//     {
//       rThis.PrintInfo(rOStream);
//       rOStream << std::endl;
//       rThis.PrintData(rOStream);

//       return rOStream;
//     }
//   ///@}


}  // namespace Kratos.

#endif // KRATOS_MODEL_PART_IO_H_INCLUDED  defined 
