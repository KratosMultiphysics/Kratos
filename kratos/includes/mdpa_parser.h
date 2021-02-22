#include <set>
#include <vector>
#include <string>
#include <unordered_map>
#include <stdio.h>

#include "containers/array_1d.h"

#include "includes/serializer.h"
#include "includes/element.h"
#include "includes/condition.h"

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
        typedef array_1d<double, 3>                                 CoordType;              // ( X, Y, Z )
        typedef std::tuple<std::size_t, CoordType>                  NodeType;               // ( Id, Coords )
        typedef std::vector<MemoryMdpa::NodeType>                   NodeContainerType;      // [ NodeType ]

        // Elements / conditions
        typedef std::vector<int>                                    ConnectType;            // [ Id, Porp_Id, Node1_Id, ... , NodeN_Id ]
        typedef std::vector<ConnectType>                            ConnectVectorType;      // [ Connectivity ]
        typedef std::unordered_map<std::string, ConnectVectorType>  EntityContainerType;    // ( Entity Name, Connectivity Vector )

        // SubModelparts
        typedef std::vector<SubMdpa::Pointer>                       SubMdpaContainerType;   // [ SubMdpaPtr ]

        // Data
        typedef std::tuple<int, double>                             DataType;               // ( Id, value )
        typedef std::tuple<int, bool, double>                       StepDataType;           // ( Id, fixity, value )
        typedef std::unordered_map<std::string, DataType>           DataValuesType;         // [ DataType ]
        typedef std::unordered_map<std::string, StepDataType>       StepDataValuesType;     // [ StepDataType ]
        
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

            // std::stringstream stream_line(line);

            // stream_line >> id;
            // stream_line >> coords[0];
            // stream_line >> coords[1];
            // stream_line >> coords[2];

            sscanf(line.c_str(), "%i %lf %lf %lf", &id, &coords[0], &coords[1], &coords[2]);

            rMdpa.mNodes.push_back(MemoryMdpa::NodeType(id, coords));
        }

        std::cout << "[DEBUG] Finished Node Block" << std::endl;
    }

    static inline void FastScan4(std::string& rLine, std::vector<int>& rData) {
        sscanf(rLine.c_str(), "%i %i %i %i", &rData[0], &rData[1], &rData[2], &rData[3]);
    }

    static inline void FastScan5(std::string& rLine, std::vector<int>& rData) {
        sscanf(rLine.c_str(), "%i %i %i %i %i", &rData[0], &rData[1], &rData[2], &rData[3], &rData[4]);
    }

    static inline void FastScan6(std::string& rLine, std::vector<int>& rData) {
        sscanf(rLine.c_str(), "%i %i %i %i %i %i", &rData[0], &rData[1], &rData[2], &rData[3], &rData[4], &rData[5]);
    }

    void ReadElemBlock(MemoryMdpa &rMdpa, std::istream &rFile, std::string &rName) {
        std::string line;

        std::cout << "[DEBUG] Reading Element Block: " << rName << std::endl;

        auto & elems = rMdpa.mElems[rName];
        KratosElemType const& base_elem = KratosComponents<KratosElemType>::Get(rName);
        std::size_t n_nodes = base_elem.GetGeometry().size();

        std::vector<int> mElemData(2 + n_nodes);

        while(!rFile.eof()) {
            getTrimLine(rFile, line);

            if(line == "End Elements")
                break; // Block end.

            // std::stringstream stream_line(line);

            // stream_line >> mElemData[0];        // Id
            // stream_line >> mElemData[1];        // Property Id

            // for(std::size_t i = 0; i < n_nodes; i++) {
            //     stream_line >> mElemData[2+i];  // Node Id
            // }

                 if(2 + n_nodes == 4) MdpaReader::FastScan4(line, mElemData);
            else if(2 + n_nodes == 5) MdpaReader::FastScan5(line, mElemData);
            else if(2 + n_nodes == 6) MdpaReader::FastScan6(line, mElemData);
            else {
                // Failback
                for(int i = 0; i < 2 + n_nodes; i++) {
                    sscanf(line.c_str(), "%i", &mElemData[i]);
                }
            }

            elems.push_back(mElemData);
        }

        std::cout << "[DEBUG] Finished Element Block" << std::endl;
    }

    void ReadCondBlock(MemoryMdpa &rMdpa, std::istream &rFile, std::string &rName) {
        std::string line;

        std::cout << "[DEBUG] Reading Condition Block: " << rName << std::endl;

        KratosCondType const& base_cond = KratosComponents<KratosCondType>::Get(rName);
        std::size_t n_nodes = base_cond.GetGeometry().size();

        std::vector<int> mCondData(2+n_nodes);

        while(!rFile.eof()) {
            getTrimLine(rFile, line);

            if(line == "End Conditions")
                break; // Block end.

            // std::stringstream stream_line(line);

            // stream_line >> mCondData[0];        // Id
            // stream_line >> mCondData[1];        // Property Id

            // for(std::size_t i = 0; i < n_nodes; i++) {
            //     stream_line >> mCondData[2+i];  // Node Id
            // }

                 if(2 + n_nodes == 4) FastScan4(line, mCondData);
            else if(2 + n_nodes == 5) FastScan5(line, mCondData);
            else if(2 + n_nodes == 6) FastScan6(line, mCondData);
            else {
                // Failback
                for(int i = 0; i < 2 + n_nodes; i++) {
                    sscanf(line.c_str(), "%i", &mCondData[i]);
                }
            }

            rMdpa.mConds[rName].push_back(mCondData);
        }

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

            // std::stringstream stream_line(line);
            // stream_line >> id;

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

            // std::stringstream stream_line(line);
            // stream_line >> id;

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

            // std::stringstream stream_line(line);
            // stream_line >> id;

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
    }

    void PrintSubMdpa(SubMdpa::Pointer pSubMpda, std::string printPadding, int to_print) {
        std::cout << printPadding << pSubMpda->mName << std::endl;
        std::cout << printPadding << "-- Nodes: " << pSubMpda->mNodeIds.size() << std::endl;
        for(std::size_t i = 0; i < std::fmin(to_print, pSubMpda->mNodeIds.size()); i++) {
            std::cout << printPadding << "\t" << pSubMpda->mNodeIds[i] << std::endl;
        }
        std::cout << printPadding << "\t..." << std::endl;
        for(std::size_t i = std::fmax(0, pSubMpda->mNodeIds.size() - to_print); i < pSubMpda->mNodeIds.size(); i++) {
            std::cout << printPadding << "\t" << pSubMpda->mNodeIds[i] << std::endl;
        }

        std::cout << printPadding << "-- Elements: " << pSubMpda->mElemIds.size() << std::endl;
        for(std::size_t i = 0; i < std::fmin(to_print, pSubMpda->mElemIds.size()); i++) {
            std::cout << printPadding << "\t" << pSubMpda->mElemIds[i] << std::endl;
        }
        std::cout << printPadding << "\t..." << std::endl;
        for(std::size_t i = std::fmax(0, pSubMpda->mElemIds.size() - to_print); i < pSubMpda->mElemIds.size(); i++) {
            std::cout << printPadding << "\t" << pSubMpda->mElemIds[i] << std::endl;
        }

        std::cout << printPadding << "-- Conditions: " << pSubMpda->mCondIds.size() << std::endl;
        for(std::size_t i = 0; i < std::fmin(to_print, pSubMpda->mCondIds.size()); i++) {
            std::cout << printPadding << "\t" << pSubMpda->mCondIds[i] << std::endl;
        }
        std::cout << printPadding << "\t..." << std::endl;
        for(std::size_t i = std::fmax(0, pSubMpda->mCondIds.size() - to_print); i < pSubMpda->mCondIds.size(); i++) {
            std::cout << printPadding << "\t" << pSubMpda->mCondIds[i] << std::endl;
        }

        for(auto iter = pSubMpda->mSubMdpa.begin(); iter != pSubMpda->mSubMdpa.end(); ++iter) {
            PrintSubMdpa(*iter, printPadding+"\t", to_print);
        }
    }

    void PrintMdpaStats(const MemoryMdpa &rMdpa) {
        int to_print = 2;
        std::cout << "MdpaStats:" << std::endl;

        std::cout << "-- Nodes: " << rMdpa.mNodes.size() << std::endl;
        for(std::size_t i = 0; i < std::fmin(to_print, rMdpa.mNodes.size()); i++) {
            std::cout << "\t" << i << " " << std::get<0>(rMdpa.mNodes[i]) << " " << std::get<1>(rMdpa.mNodes[i])[0] << " " << std::get<1>(rMdpa.mNodes[i])[1] << " " << std::get<1>(rMdpa.mNodes[i])[2] << std::endl;
        }
        std::cout << "\t..." << std::endl;
        for(std::size_t i = std::fmax(0, rMdpa.mNodes.size() - to_print); i < rMdpa.mNodes.size(); i++) {
            std::cout << "\t" << i << " " << std::get<0>(rMdpa.mNodes[i]) << " " << std::get<1>(rMdpa.mNodes[i])[0] << " " << std::get<1>(rMdpa.mNodes[i])[1] << " " << std::get<1>(rMdpa.mNodes[i])[2] << std::endl;
        }

        std::cout << "-- Elements: " << std::endl;
        for(auto iter = rMdpa.mElems.begin(); iter != rMdpa.mElems.end(); ++iter) {
            std::cout << "--- "  << iter->first << ": " << iter->second.size() << std::endl;
            for(std::size_t i = 0; i < std::fmin(to_print, rMdpa.mNodes.size()); i++) {
                std::cout << "\t" << i << " " << iter->second[i] << std::endl;
            }
            std::cout << "\t..." << std::endl;
            for(std::size_t i = std::fmax(0, iter->second.size() - to_print); i < iter->second.size(); i++) {
                std::cout << "\t" << i << " " << iter->second[i] << std::endl;
            }
        }

        std::cout << "-- Conditions: " << std::endl;
        for(auto iter = rMdpa.mConds.begin(); iter != rMdpa.mConds.end(); ++iter) {
            std::cout << "--- "  << iter->first << ": " << iter->second.size() << std::endl;
            for(std::size_t i = 0; i < std::fmin(to_print, rMdpa.mNodes.size()); i++) {
                std::cout << "\t" << i << " " << iter->second[i] << std::endl;
            }
            std::cout << "\t..." << std::endl;
            for(std::size_t i = std::fmax(0, iter->second.size() - to_print); i < iter->second.size(); i++) {
                std::cout << "\t" << i << " " << iter->second[i] << std::endl;
            }
        }

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
