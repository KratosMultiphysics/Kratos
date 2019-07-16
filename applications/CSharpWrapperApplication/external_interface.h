//    _____  _____ _                  __          __                                                 _ _           _   _
//   / ____|/ ____| |                 \ \        / /                               /\               | (_)         | | (_)
//  | |    | (___ | |__   __ _ _ __ _ _\ \  /\  / / __ __ _ _ __  _ __   ___ _ __ /  \   _ __  _ __ | |_  ___ __ _| |_ _  ___  _ __
//  | |     \___ \| '_ \ / _` | '__| '_ \ \/  \/ / '__/ _` | '_ \| '_ \ / _ \ '__/ /\ \ | '_ \| '_ \| | |/ __/ _` | __| |/ _ \| '_  |
//  | |____ ____) | | | | (_| | |  | |_) \  /\  /| | | (_| | |_) | |_) |  __/ | / ____ \| |_) | |_) | | | (_| (_| | |_| | (_) | | | |
//   \_____|_____/|_| |_|\__,_|_|  | .__/ \/  \/ |_|  \__,_| .__/| .__/ \___|_|/_/    \_\ .__/| .__/|_|_|\___\__,_|\__|_|\___/|_| |_|
//                                 | |                     | |   | |                    | |   | |
//                                 |_|                     |_|   |_|                    |_|   |_|
//
//
//  License: BSD License
//   license: CSharpWrapperApplication/license.txt
//
//  Main authors:    Hubert Balcerzak
//                   Riccardo Rossi
//

// System includes

// External includes

// Project includes
#include "kratos_internals.h"
#include "includes/node.h"
#include "mesh_converter.h"
#include "id_translator.h"

namespace CSharpKratosWrapper {

    using NodeType = Kratos::Node<3>;
    using ModelPart = Kratos::ModelPart;

    class CSharpInterface {

    public:
        static void init(const char* MDPAFilePath, const char* JSONFilePath = NULL);
        static void initWithSettings(const char* JSONFilePath = NULL);
        static void updateNodePos(const int nodeId, const float x, const float y, const float z);
        static void calculate();
        static float* getXCoordinates();
        static float* getYCoordinates();
        static float* getZCoordinates();
        static int getNodesCount();
        static int* getTriangles();
        static int getTrianglesCount();

    private:
        static int* pmKratosNodeIds;
        static int* pmUnityNodeIds;
        static int* pmTriangles;
        static int mTrianglesCount;
        static int mNodesCount;
        static float* pmXCoordinates;
        static float* pmYCoordinates;
        static float* pmZCoordinates;
        static KratosInternals mKratosInternals;
        static IdTranslator idTranslator;

        static std::vector<NodeType::Pointer> mFixedNodes;

        static void saveTriangles(MeshConverter& meshConverter);
        static void saveNodes(MeshConverter& meshConverter);
        static void retrieveNodesPos();
        static void freeNodes();
    };
}
