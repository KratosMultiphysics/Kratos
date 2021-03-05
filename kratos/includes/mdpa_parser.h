#include <set>
#include <vector>
#include <string>
#include <unordered_map>

#include <stdio.h>
#include <unistd.h>

#include "containers/array_1d.h"
#include "containers/model.h"
#include "containers/sparse_graph.h"

#include "includes/serializer.h"
#include "includes/element.h"
#include "includes/condition.h"

#include "includes/model_part.h"

#include "utilities/parallel_utilities.h"

namespace Kratos {

class SubMdpa {
    public:
        // Pointer definition of SubMdpa
        KRATOS_CLASS_POINTER_DEFINITION(SubMdpa);

        // Life cycle
        SubMdpa(std::string name) : mName(name) {}

        std::string                     mName;

        std::vector<std::size_t>        mNodeIds;
        std::vector<std::size_t>        mElemIds;
        std::vector<std::size_t>        mCondIds;
        std::vector<SubMdpa::Pointer>   mSubMdpa;
};

class MemoryMdpa {
    public:
        // Tuple Enums
        enum class NODE_TYPE {
            ID = 0,
            COORDS = 1,
        };

        enum class ENTITY_CONTAINER_TYPE {
            NAME = 0,
            DATA = 1,
        };

        enum class DATA_TYPE {
            ID = 0,
            VALUE = 1
        };

        enum class STEP_DATA_TYPE {
            ID = 0,
            FIXITY = 1,
            VALUE = 2
        };

        // Nodes
        typedef array_1d<double, 3>                                     CoordType;              // ( X, Y, Z )
        typedef std::tuple<std::size_t, CoordType>                      NodeType;               // ( Id, Coords )
        typedef std::vector<MemoryMdpa::NodeType>                       NodeContainerType;      // [ NodeType ]

        // Elements / conditions
        typedef std::tuple<int, int>                                    EntityMetadataType;     // ( Id, Prop ) 
        typedef std::vector<ModelPart::IndexType>                       EnityCoordTyp;          // [ Coords ]
        typedef std::tuple<
            std::vector<EntityMetadataType>,
            std::vector<EnityCoordTyp>
        >                                                               ConnectType;            // ( [ EntityMetadataType ], [ EnityCoordTyp ] )
        typedef std::unordered_map<std::string, ConnectType>            EntityContainerType;    // ( Entity Name, ConnectType )

        // typedef std::tuple<int, int, std::vector<ModelPart::IndexType>> ConnectType;            // ( Id, Porp_Id, [ Node1_Id, ... , NodeN_Id ] )
        // typedef std::vector<MemoryMdpa::ConnectType>                    ConnectVectorType;      // [ Connectivity ]
        // typedef std::unordered_map<std::string, ConnectVectorType>      EntityContainerType;    // ( Entity Name, ConnectVectorType )

        // SubModelparts
        typedef std::vector<SubMdpa::Pointer>                           SubMdpaContainerType;   // [ SubMdpaPtr ]

        // Data
        typedef std::tuple<int, double>                                 DataType;               // ( Id, value )
        typedef std::tuple<int, bool, double>                           StepDataType;           // ( Id, fixity, value )
        typedef std::unordered_map<std::string, DataType>               DataValuesType;         // [ DataType ]
        typedef std::unordered_map<std::string, StepDataType>           StepDataValuesType;     // [ StepDataType ]
        
        // Life cycle
        MemoryMdpa() {}

        friend class Serializer;

        void save(Serializer& rSerializer) const;

        void load(Serializer& rSerializer);

        std::string                     mProperties;            // String with the properties

        NodeContainerType               mNodes;                 // Id, x, y, z std::set<std::pair<std::size_t, Kratos::array_1d<3>>>
        EntityContainerType             mElems;                 // std::vector<std::pair<std::string, std::set<std::vector<int>>>>
                                                                // [ ( Element Name , ( id , property_id, node1 ... nodeN )) ]
        EntityContainerType             mConds;                 // std::vector<std::pair<std::string, std::set<std::vector<int>>>>
                                                                // [ ( Condition Name , ( id , property_id, node1 ... nodeN )) ]

        SubMdpaContainerType            mSubMdpa;               // std::vector<SubMdpa>

        DataValuesType                  mNodalSolutionStepData; // std::unordered_map<std::string, std::tuple<int, bool, double>> --> Variablename,  list(Id,fixiy,value)
        StepDataValuesType              mNodalData;             // std::unordered_map<std::string, std::tuple<int, double>> --> Variablename,  list(Id,value)

        int                             mNumNodes = 0;
        int                             mNumElems = 0;
        int                             mNumConds = 0;
};

class MdpaReader {
public:
    typedef Element     KratosElemType;
    typedef Condition   KratosCondType;

    // Pointer definition of MdpaReader
    KRATOS_CLASS_POINTER_DEFINITION(MdpaReader);

    MdpaReader() {};
    ~MdpaReader() {};

    enum LINE {
        SECTION = 1,    // Name of the section
        NAME = 2        // Name of the entity in case there is any
    };

    // trim from start (in place)
    static inline void ltrim(std::string &s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
            return !std::isspace(ch);
        }));
    }

    // trim from end (in place)
    static inline void rtrim(std::string &s) {
        s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
            return !std::isspace(ch);
        }).base(), s.end());
    }

    // trim to end (in place)
    static inline void trimComment(std::string &s, std::string delimiter) {
        s.erase(s.begin()+std::min(s.find(delimiter), s.size()), s.end()); // This is MIN and not FMIN !!!!!!DO NOT CHANGE IT!!!!!!
    }

    static inline void getTrimLine(std::istream &rFile, std::string &line) {
        std::getline(rFile, line);
        trimComment(line, "//");    // Strip comments
        ltrim(line);                // Need this for the tabs

        // I Don't think we need this but just in case.
        // rtrim(line);
    }

    void Read(const std::string mdpaFilename) {
        std::ifstream mdpa_file;

        mdpa_file.open(mdpaFilename);
        ReadFromStream(mdpa_file);
        mdpa_file.close();
    }

    template<class TDataType>
    void ReadValue(std::istream &rLine, TDataType &rData) {
        rLine >> rData;
    }

    template<class TDataType>
    void ExtValue(std::string &rWord, TDataType &rData) {
        std::stringstream word(rWord);
        word >> rData;
    }

    void ReadSection(std::istream &rFile, std::vector<std::string> &rSection) {
        std::string line, word;
        getTrimLine(rFile, line);
        std::cout << "[SECTION] " << line << std::endl;
        std::stringstream sline(line);

        for(int i = 0; i < 2; i++) {
            std::getline(sline, word, ' ');
            if(word[i] == '\n')
                break;
            rSection[i] = word;
        }
        
        if(rSection[1] == "Elements" || rSection[1] == "Conditions" || rSection[1] == "SubModelPart") {
            std::getline(sline, word, ' ');
            rSection[2] = word;
        }
    }

    void ReadDataBlock(MemoryMdpa &rMdpa, std::istream &rFile, std::string &rName) {
        std::string line;

        std::cout << "[DEBUG] Reading Data Block" << std::endl;

        while(!rFile.eof()) {
            getTrimLine(rFile, line);
            std::cout << "[DEBUG] \t" << line << std::endl;

            if(line == "End ModelPartData")
                break; // Block end.

            std::stringstream stream_line(line);
        }

        std::cout << "[DEBUG] Finished Data Block" << std::endl;
    }

    void ReadPropBlock(MemoryMdpa &rMdpa, std::istream &rFile) {
        std::string line;

        std::cout << "[DEBUG] Reading Properties Block" << std::endl;

        while(!rFile.eof()) {
            getTrimLine(rFile, line);

            if(line == "End Properties")
                break; // Block end.
        }

        std::cout << "[DEBUG] Finished Properties Block" << std::endl;
    }

    void ReadNodeBlock(MemoryMdpa &rMdpa, std::istream &rFile) {
        std::string line;

        std::cout << "[DEBUG] Reading Node Block" << std::endl;

        int id;
        array_1d<double, 3> coords;

        while(!rFile.eof()) {
            getTrimLine(rFile, line);

            if(line == "End Nodes")
                break; // Block end.

            sscanf(line.c_str(), "%i %lf %lf %lf", &id, &coords[0], &coords[1], &coords[2]);

            rMdpa.mNodes.push_back(MemoryMdpa::NodeType(id, coords));
        }

        rMdpa.mNumNodes += rMdpa.mNodes.size();

        std::cout << "[DEBUG] Finished Node Block" << std::endl;
    }

    // 2 Points (2D2N, ...)
    static inline void FastScan4(std::string& rLine, int& rId, int& rProp, std::vector<ModelPart::IndexType>& rData) {
        sscanf(rLine.c_str(), "%i %i %lu %lu", &rId, &rProp, &rData[0], &rData[1]);
    }

    // 3 Points (2D3N, 3D3N, ...)
    static inline void FastScan5(std::string& rLine, int& rId, int& rProp, std::vector<ModelPart::IndexType>& rData) {
        sscanf(rLine.c_str(), "%i %i %lu %lu %lu", &rId, &rProp, &rData[0], &rData[1], &rData[2]);
    }

    // 4 Points (2D4N, 3D4N, ...)
    static inline void FastScan6(std::string& rLine, int& rId, int& rProp, std::vector<ModelPart::IndexType>& rData) {
        sscanf(rLine.c_str(), "%i %i %lu %lu %lu %lu", &rId, &rProp, &rData[0], &rData[1], &rData[2], &rData[3]);
    }

    void ReadElemBlock(MemoryMdpa &rMdpa, std::istream &rFile, std::string &rName) {
        std::string line;

        std::cout << "[DEBUG] Reading Element Block: " << rName << std::endl;

        auto & elems = rMdpa.mElems[rName];
        KratosElemType const& base_elem = KratosComponents<KratosElemType>::Get(rName);
        std::size_t n_nodes = base_elem.GetGeometry().size();

        int id, prop;
        std::vector<ModelPart::IndexType> mElemData(n_nodes);

        while(!rFile.eof()) {
            getTrimLine(rFile, line);

            if(line == "End Elements")
                break; // Block end.

                 if(2 + n_nodes == 4) MdpaReader::FastScan4(line, id, prop, mElemData);
            else if(2 + n_nodes == 5) MdpaReader::FastScan5(line, id, prop, mElemData);
            else if(2 + n_nodes == 6) MdpaReader::FastScan6(line, id, prop, mElemData);
            else {
                // Failback
                sscanf(line.c_str(), "%i %i", &id, &prop);
                for(std::size_t i = 0; i < n_nodes; i++) {
                    sscanf(line.c_str(), "%lu", &mElemData[i]);
                }
            }

            std::get<0>(elems).push_back(MemoryMdpa::EntityMetadataType(id, prop));
            std::get<1>(elems).push_back(mElemData);
        }

        rMdpa.mNumElems += std::get<0>(elems).size();

        std::cout << "[DEBUG] Finished Element Block" << std::endl;
    }

    void ReadCondBlock(MemoryMdpa &rMdpa, std::istream &rFile, std::string &rName) {
        std::string line;

        std::cout << "[DEBUG] Reading Condition Block: " << rName << std::endl;

        auto & conds = rMdpa.mConds[rName];
        KratosCondType const& base_cond = KratosComponents<KratosCondType>::Get(rName);
        std::size_t n_nodes = base_cond.GetGeometry().size();

        int id, prop;
        std::vector<ModelPart::IndexType> mCondData(n_nodes);

        while(!rFile.eof()) {
            getTrimLine(rFile, line);

            if(line == "End Conditions")
                break; // Block end.

                 if(2 + n_nodes == 4) FastScan4(line, id, prop, mCondData);
            else if(2 + n_nodes == 5) FastScan5(line, id, prop, mCondData);
            else if(2 + n_nodes == 6) FastScan6(line, id, prop, mCondData);
            else {
                // Failback
                sscanf(line.c_str(), "%i %i", &id, &prop);
                for(std::size_t i = 0; i < n_nodes; i++) {
                    sscanf(line.c_str(), "%lu", &mCondData[i]);
                }
            }

            std::get<0>(conds).push_back(MemoryMdpa::EntityMetadataType(id, prop));
            std::get<1>(conds).push_back(mCondData);
        }

        rMdpa.mNumConds += std::get<0>(conds).size();

        std::cout << "[DEBUG] Finished Condition Block" << std::endl;
    }

    void ReadSubMNodeBlock(SubMdpa::Pointer pSubMdpa, std::istream &rFile) {
        std::string line;

        std::cout << "[DEBUG] Reading SubModelpart Node Block" << std::endl;

        int id;

        while(!rFile.eof()) {
            getTrimLine(rFile, line);

            if(line == "End SubModelPartNodes")
                break; // Block end.

            sscanf(line.c_str(), "%i", &id);

            pSubMdpa->mNodeIds.push_back(id);
        }

        std::cout << "[DEBUG] Finished SubModelpart Node Block" << std::endl;
    }

    void ReadSubMElemBlock(SubMdpa::Pointer pSubMdpa, std::istream &rFile) {
        std::string line;

        std::cout << "[DEBUG] Reading SubModelpart Elem Block" << std::endl;

        int id;

        while(!rFile.eof()) {
            getTrimLine(rFile, line);

            if(line == "End SubModelPartElements")
                break; // Block end.

            sscanf(line.c_str(), "%i", &id);

            pSubMdpa->mElemIds.push_back(id);
        }

        std::cout << "[DEBUG] Finished SubModelpart Elem Block" << std::endl;
    }

    void ReadSubMCondBlock(SubMdpa::Pointer pSubMdpa, std::istream &rFile) {
        std::string line;

        std::cout << "[DEBUG] Reading SubModelpart Cond Block" << std::endl;

        int id;

        while(!rFile.eof()) {
            getTrimLine(rFile, line);

            if(line == "End SubModelPartConditions")
                break; // Block end.

            sscanf(line.c_str(), "%i", &id);

            pSubMdpa->mCondIds.push_back(id);
        }

        std::cout << "[DEBUG] Finished SubModelpart Cond Block" << std::endl;
    }


    void ReadSubMBlock(SubMdpa::Pointer pSubMdpa, std::istream &rFile) {
        std::string line;
        std::vector<std::string> line_v(3);

        std::cout << "[DEBUG] Reading Submodelpart Block: " << pSubMdpa->mName << std::endl;

        while(!rFile.eof()) {
            getTrimLine(rFile, line);

            std::cout << "[DEBUG] Reading line subsection: " << line << std::endl;
            
            if(line == "End SubModelPart")
                break; // Block end.

            std::stringstream streamline(line);
            
            // Obtain the submodelpart section
            ReadSection(streamline, line_v);

                   if(line_v[MdpaReader::LINE::SECTION] == "SubModelPart") {
                pSubMdpa->mSubMdpa.push_back(std::make_shared<SubMdpa>(line_v[MdpaReader::LINE::NAME]));
                ReadSubMBlock(pSubMdpa->mSubMdpa.back(), rFile);
            } else if(line_v[MdpaReader::LINE::SECTION] == "SubModelPartNodes") {
                ReadSubMNodeBlock(pSubMdpa, rFile);
            } else if(line_v[MdpaReader::LINE::SECTION] == "SubModelPartElements") {
                ReadSubMElemBlock(pSubMdpa, rFile);
            } else if(line_v[MdpaReader::LINE::SECTION] == "SubModelPartConditions") {
                ReadSubMCondBlock(pSubMdpa, rFile);
            } 
        }

        std::cout << "[DEBUG] Finished Submodelpart Block" << std::endl;
    }

    void ReadFromStream(std::istream &rFile) {
        MemoryMdpa mdpa;
        std::vector<std::string> line_v(3);

        while(!rFile.eof()) {            
            // Obtain the section
            ReadSection(rFile, line_v);

            std::cout << "[DEBUG] Reading line " << line_v[0] << " " << line_v[1] << std::endl;

                   if(line_v[MdpaReader::LINE::SECTION] == "Properties") {
                ReadPropBlock(mdpa, rFile);
            } else if(line_v[MdpaReader::LINE::SECTION] == "Nodes") {
                ReadNodeBlock(mdpa, rFile);
            } else if(line_v[MdpaReader::LINE::SECTION] == "Elements") {
                ReadElemBlock(mdpa, rFile, line_v[MdpaReader::LINE::NAME]);
            } else if(line_v[MdpaReader::LINE::SECTION] == "Conditions") {
                ReadCondBlock(mdpa, rFile, line_v[MdpaReader::LINE::NAME]);
            } else if(line_v[MdpaReader::LINE::SECTION] == "SubModelPart") {
                mdpa.mSubMdpa.push_back(std::make_shared<SubMdpa>(line_v[MdpaReader::LINE::NAME]));
                ReadSubMBlock(mdpa.mSubMdpa.back(), rFile);
            } else if(line_v[MdpaReader::LINE::SECTION] == "ModelPartData") {
                ReadDataBlock(mdpa, rFile, line_v[MdpaReader::LINE::NAME]);
            }
        }

        PrintMdpaStats(mdpa);

        std::vector<ModelPart::IndexType> row_index, col_index;

        // usleep(2 * 1e6);

        GenerateCSR(mdpa, row_index, col_index);

        // usleep(2 * 1e6);

        // GenerateModelPart(mdpa);

        DynamicallyGenerateModelPart(mdpa);

        // usleep(2 * 1e6);
    }

    void GenerateModelPart(const MemoryMdpa& rMdpa) {

        Model model;
        ModelPart& model_part = model.CreateModelPart("ConsistentModelPart");

        model_part.AddNodalSolutionStepVariable(PRESSURE);
        model_part.AddNodalSolutionStepVariable(TEMPERATURE);
        model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

        std::map<int, ModelPart::NodeType::Pointer> nodeMap;

        for(auto & node: rMdpa.mNodes) {
            auto new_node = model_part.CreateNewNode(
                std::get<0>(node),
                std::get<1>(node)[0],
                std::get<1>(node)[1],
                std::get<1>(node)[2]
            );
            nodeMap[std::get<0>(node)] = new_node;
        }

        // TODO: OpenMP
        ModelPart::ElementsContainerType aux_elems;

        for(auto & entity: rMdpa.mElems) {
            auto elem_name = entity.first;
            auto elem_list = entity.second;

            ModelPart::ElementType const& r_base_elem = KratosComponents<ModelPart::ElementType>::Get(elem_name);
            Properties::Pointer p_prop = Kratos::make_shared<Properties>(0);

            std::size_t geom_size = r_base_elem.GetGeometry().size();

            Geometry< Node < 3 > >::PointsArrayType pElemNodes(geom_size);
            auto& pElemNodesContainer = pElemNodes.GetContainer();

            if(std::get<0>(elem_list).size() != std::get<1>(elem_list).size()) {
                KRATOS_ERROR << "Inconsistent element list size: " << std::get<0>(elem_list).size() << " indices found but " << std::get<0>(elem_list).size() << " coords found." << std::endl; 
            }

            auto meta_itr = std::get<0>(elem_list).begin();
            auto coords_itr = std::get<1>(elem_list).begin();
        
            for(unsigned int i = 0; i < std::get<0>(elem_list).size(); i++)  {
                for (unsigned int j = 0; j < geom_size; j++) {
                    pElemNodesContainer[j] = nodeMap[(*coords_itr)[j]];
                }

                aux_elems.push_back(r_base_elem.Create(std::get<0>(*meta_itr), pElemNodes, p_prop));

                meta_itr++;
                coords_itr++;
            }

            // for(auto & elem: elem_list) {
            //     auto& r_elem = std::get<2>(elem);
            //     for (unsigned int i = 0; i < geom_size; i++) {
            //         pElemNodesContainer[i] = nodeMap[r_elem[i]];
            //     }

            //     aux_elems.push_back(r_base_elem.Create(std::get<0>(elem), pElemNodes, p_prop));
            // }
        }


        /////////////////////////////

        aux_elems.Unique();
        model_part.AddElements(aux_elems.begin(), aux_elems.end());

        ModelPart::ConditionsContainerType aux_conds;
        for(auto & entity: rMdpa.mConds) {  
            auto cond_name = entity.first;
            auto cond_list = entity.second;

            ModelPart::ConditionType const& r_base_cond = KratosComponents<ModelPart::ConditionType>::Get(cond_name);
            Properties::Pointer p_prop = Kratos::make_shared<Properties>(0);

            std::size_t geom_size = r_base_cond.GetGeometry().size();
            
            Geometry< Node < 3 > >::PointsArrayType pCondNodes(r_base_cond.GetGeometry().size());
            auto& pCondNodesContainer = pCondNodes.GetContainer();

            if(std::get<0>(cond_list).size() != std::get<1>(cond_list).size()) {
                KRATOS_ERROR << "Inconsistent condition list size: " << std::get<0>(cond_list).size() << " indices found but " << std::get<0>(cond_list).size() << " coords found." << std::endl; 
            }

            auto meta_itr = std::get<0>(cond_list).begin();
            auto coords_itr = std::get<1>(cond_list).begin();
        
            for(unsigned int i = 0; i < std::get<0>(cond_list).size(); i++)  {
                for (unsigned int j = 0; j < geom_size; j++) {
                    pCondNodesContainer[j] = nodeMap[(*coords_itr)[j]];
                }

                aux_conds.push_back(r_base_cond.Create(std::get<0>(*meta_itr), pCondNodes, p_prop));

                meta_itr++;
                coords_itr++;
            }

            // for(auto & cond: cond_list) {
            //     auto& r_cond = std::get<2>(cond);
            //     for (unsigned int i = 0; i < r_base_cond.GetGeometry().size(); i++) {
            //         pCondNodesContainer[i] = nodeMap[r_cond[i]];
            //     }

            //     aux_conds.push_back(r_base_cond.Create(std::get<0>(cond), pCondNodes, p_prop));
            // }

            std::cout << "[DEBUG]: readed " << std::get<0>(cond_list).size() << " conditions" << std::endl;
        }

        aux_conds.Unique();
        model_part.AddConditions(aux_conds.begin(), aux_conds.end());
    }

    void DynamicallyGenerateModelPart(MemoryMdpa& rMdpa) {

        Model model;
        ModelPart& model_part = model.CreateModelPart("ConsistentModelPart");

        // model_part.AddNodalSolutionStepVariable(PRESSURE);
        // model_part.AddNodalSolutionStepVariable(TEMPERATURE);
        // model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

        std::map<int, ModelPart::NodeType::Pointer> nodeMap;

        for(auto & node: rMdpa.mNodes) {
            auto new_node = model_part.CreateNewNode(
                std::get<0>(node),
                std::get<1>(node)[0],
                std::get<1>(node)[1],
                std::get<1>(node)[2]
            );
            nodeMap[std::get<0>(node)] = new_node;
        }

        MemoryMdpa::NodeContainerType().swap(rMdpa.mNodes);

        // TODO: OpenMP
        ModelPart::ElementsContainerType aux_elems;

        for(auto & entity: rMdpa.mElems) {
            auto elem_name = entity.first;
            auto elem_list = entity.second;

            ModelPart::ElementType const& r_base_elem = KratosComponents<ModelPart::ElementType>::Get(elem_name);
            Properties::Pointer p_prop = Kratos::make_shared<Properties>(0);

            std::size_t geom_size = r_base_elem.GetGeometry().size();

            Geometry< Node < 3 > >::PointsArrayType pElemNodes(geom_size);
            auto& pElemNodesContainer = pElemNodes.GetContainer();

            if(std::get<0>(elem_list).size() != std::get<1>(elem_list).size()) {
                KRATOS_ERROR << "Inconsistent element list size: " << std::get<0>(elem_list).size() << " indices found but " << std::get<0>(elem_list).size() << " coords found." << std::endl; 
            }

            int icount = 0;
            int free_chunk = std::get<0>(elem_list).size();

            while(!std::get<0>(elem_list).empty()) {
                auto elem_m = std::get<0>(elem_list).back();
                auto elem_c = std::get<1>(elem_list).back();

                for (unsigned int i = 0; i < geom_size; i++) {
                    pElemNodesContainer[i] = nodeMap[elem_c[i]];
                }

                aux_elems.push_back(r_base_elem.Create(std::get<0>(elem_m), pElemNodes, p_prop));

                std::get<0>(elem_list).pop_back();
                std::get<1>(elem_list).pop_back();
            }

            std::get<0>(elem_list).resize(0);
            std::get<1>(elem_list).resize(0);
        }

        MemoryMdpa::EntityContainerType().swap(rMdpa.mElems);

        /////////////////////////////

        aux_elems.Unique();
        model_part.AddElements(aux_elems.begin(), aux_elems.end());
        ModelPart::ElementsContainerType().swap(aux_elems);

        ModelPart::ConditionsContainerType aux_conds;

        for(auto & entity: rMdpa.mConds) {  
            auto cond_name = entity.first;
            auto cond_list = entity.second;

            ModelPart::ConditionType const& r_base_cond = KratosComponents<ModelPart::ConditionType>::Get(cond_name);
            Properties::Pointer p_prop = Kratos::make_shared<Properties>(0);

            std::size_t geom_size = r_base_cond.GetGeometry().size();
            
            Geometry< Node < 3 > >::PointsArrayType pCondNodes(r_base_cond.GetGeometry().size());
            auto& pCondNodesContainer = pCondNodes.GetContainer();

            if(std::get<0>(cond_list).size() != std::get<1>(cond_list).size()) {
                KRATOS_ERROR << "Inconsistent condition list size: " << std::get<0>(cond_list).size() << " indices found but " << std::get<0>(cond_list).size() << " coords found." << std::endl; 
            }
        
            int icount = 0;
            int free_chunk = std::get<0>(cond_list).size();

            while(!std::get<0>(cond_list).empty()) {
                auto cond_m = std::get<0>(cond_list).back();
                auto cond_c = std::get<1>(cond_list).back();

                for (unsigned int i = 0; i < geom_size; i++) {
                    pCondNodesContainer[i] = nodeMap[cond_c[i]];
                }

                aux_conds.push_back(r_base_cond.Create(std::get<0>(cond_m), pCondNodes, p_prop));

                std::get<0>(cond_list).pop_back();
                std::get<1>(cond_list).pop_back();
            }

            std::cout << "[DEBUG]: readed " << std::get<0>(cond_list).size() << " conditions" << std::endl;

            std::get<0>(cond_list).resize(0);
            std::get<1>(cond_list).resize(0);
        }

        MemoryMdpa::EntityContainerType().swap(rMdpa.mConds);

        aux_conds.Unique();
        model_part.AddConditions(aux_conds.begin(), aux_conds.end());
        ModelPart::ConditionsContainerType().swap(aux_conds);
    }

    void GenerateCSR(const MemoryMdpa &rMdpa, std::vector<ModelPart::IndexType> &rRowIndex, std::vector<ModelPart::IndexType> &rColIndex) {
        SparseGraph<> Agraph;
        for(auto & entity: rMdpa.mElems) {
            for(const auto& c : std::get<1>(entity.second)) {
                Agraph.AddEntries(c);
            }
        }
        Agraph.Finalize();

        auto nrows = Agraph.ExportCSRArrays(rRowIndex, rColIndex);

        std::cout << "CSR Rows: " << rRowIndex.size() << std::endl;
        std::cout << "CSR Cols: " << rColIndex.size() << std::endl;
    }

    void PrintSubMdpa(SubMdpa::Pointer pSubMdpa, std::string printPadding, int to_print) {
        std::cout << printPadding << pSubMdpa->mName << std::endl;
        std::cout << printPadding << "-- Nodes: " << pSubMdpa->mNodeIds.size() << std::endl;
        for(std::size_t i = 0; i < std::fmin(to_print, pSubMdpa->mNodeIds.size()); i++) {
            std::cout << printPadding << "\t" << pSubMdpa->mNodeIds[i] << std::endl;
        }
        std::cout << printPadding << "\t..." << std::endl;
        for(std::size_t i = std::fmax(0, pSubMdpa->mNodeIds.size() - to_print); i < pSubMdpa->mNodeIds.size(); i++) {
            std::cout << printPadding << "\t" << pSubMdpa->mNodeIds[i] << std::endl;
        }

        std::cout << printPadding << "-- Elements: " << pSubMdpa->mElemIds.size() << std::endl;
        for(std::size_t i = 0; i < std::fmin(to_print, pSubMdpa->mElemIds.size()); i++) {
            std::cout << printPadding << "\t" << pSubMdpa->mElemIds[i] << std::endl;
        }
        std::cout << printPadding << "\t..." << std::endl;
        for(std::size_t i = std::fmax(0, pSubMdpa->mElemIds.size() - to_print); i < pSubMdpa->mElemIds.size(); i++) {
            std::cout << printPadding << "\t" << pSubMdpa->mElemIds[i] << std::endl;
        }

        std::cout << printPadding << "-- Conditions: " << pSubMdpa->mCondIds.size() << std::endl;
        for(std::size_t i = 0; i < std::fmin(to_print, pSubMdpa->mCondIds.size()); i++) {
            std::cout << printPadding << "\t" << pSubMdpa->mCondIds[i] << std::endl;
        }
        std::cout << printPadding << "\t..." << std::endl;
        for(std::size_t i = std::fmax(0, pSubMdpa->mCondIds.size() - to_print); i < pSubMdpa->mCondIds.size(); i++) {
            std::cout << printPadding << "\t" << pSubMdpa->mCondIds[i] << std::endl;
        }

        for(auto iter = pSubMdpa->mSubMdpa.begin(); iter != pSubMdpa->mSubMdpa.end(); ++iter) {
            PrintSubMdpa(*iter, printPadding+"\t", to_print);
        }
    }

    void PrintMdpaStats(const MemoryMdpa &rMdpa) {
        int to_print = 2;
        std::cout << "MdpaStats:" << std::endl;

        std::cout << "-- Nodes: " << rMdpa.mNodes.size() << std::endl;
        for(std::size_t i = 0; i < std::fmin(to_print, rMdpa.mNodes.size()); i++) {
            std::cout << "\t" << std::get<0>(rMdpa.mNodes[i]) << " " << std::get<1>(rMdpa.mNodes[i])[0] << " " << std::get<1>(rMdpa.mNodes[i])[1] << " " << std::get<1>(rMdpa.mNodes[i])[2] << std::endl;
        }
        std::cout << "\t..." << std::endl;
        for(std::size_t i = std::fmax(0, rMdpa.mNodes.size() - to_print); i < rMdpa.mNodes.size(); i++) {
            std::cout << "\t" << std::get<0>(rMdpa.mNodes[i]) << " " << std::get<1>(rMdpa.mNodes[i])[0] << " " << std::get<1>(rMdpa.mNodes[i])[1] << " " << std::get<1>(rMdpa.mNodes[i])[2] << std::endl;
        }

        std::cout << "-- Elements: " << std::endl;
        for(auto iter = rMdpa.mElems.begin(); iter != rMdpa.mElems.end(); ++iter) {
            std::cout << "--- "  << iter->first << ": " << std::get<0>(iter->second).size() << std::endl;
            auto meta_itr = std::get<0>(iter->second).begin();
            auto coords_itr = std::get<1>(iter->second).begin();
            for(std::size_t i = 0; i < std::fmin(to_print,  std::get<0>(iter->second).size()); i++) {
                std::cout << "\t" << std::get<0>(*(meta_itr+i)) << " " << *(coords_itr+i) << std::endl;
            }
            std::cout << "\t..." << std::endl;
            for(std::size_t i = std::fmax(0, std::get<0>(iter->second).size() - to_print); i < std::get<0>(iter->second).size(); i++) {
                std::cout << "\t" << std::get<0>(*(meta_itr+i)) << " " << *(coords_itr+i) << std::endl;
            }
        }

        // std::cout << "-- Conditions: " << std::endl;
        // for(auto iter = rMdpa.mConds.begin(); iter != rMdpa.mConds.end(); ++iter) {
        //     std::cout << "--- "  << iter->first << ": " << iter->second.size() << std::endl;
        //     for(std::size_t i = 0; i < std::fmin(to_print,  iter->second.size()); i++) {
        //         std::cout << "\t" << std::get<0>(iter->second[i]) << " " << std::get<2>(iter->second[i]) << std::endl;
        //     }
        //     std::cout << "\t..." << std::endl;
        //     for(std::size_t i = std::fmax(0, iter->second.size() - to_print); i < iter->second.size(); i++) {
        //         std::cout << "\t" << std::get<0>(iter->second[i]) << " " << std::get<2>(iter->second[i]) << std::endl;
        //     }
        // }

        std::cout << "-- SubModelParts" << std::endl;
        for(auto iter = rMdpa.mSubMdpa.begin(); iter != rMdpa.mSubMdpa.end(); ++iter) {
            PrintSubMdpa(*iter,"\t", to_print);
        }
    }
};

} // namespace Kratos

///////

// class MdpaReader {

//     void Read(std::string &MdpaFilename) {
//         MemoryMdpa my_dats;

//         Parse(my_data, mMMpda)
//         Reorder(my_data)

//         if(mpi)
//         {
//             local_data = DivdeAndCommunicate(my_data)  
//             CreateModelPart(local_data)
//         }
//         else
//             CreateModelPart(my_data)



//     }

//     void Parse

//     void Reorder

//     void DivideAndCommunicate

//     void CreateModelPart
 
// }

// MemoryMdpa a;
// MdpaReader m;
