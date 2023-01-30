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

#if !defined(CSHARP_WRAPPER_APPLICATION_MESH_CONVERTER_H_INCLUDED )
#define  CSHARP_WRAPPER_APPLICATION_MESH_CONVERTER_H_INCLUDED

// System includes
#include <vector>

// External includes

// Project includes
#include "includes/model_part.h"

namespace CSharpKratosWrapper {

    struct face {
        int nodes[4];
    };

    class MeshConverter {
    public:
        /**
         * Generates surface mesh from elements. All elements must be tetrahedral
         * @param elements Vector of tetrahedral elements
         */
        void ProcessMesh(std::vector<Kratos::intrusive_ptr<Kratos::Element>> &elements);

        /**
         * Allows access to faces of generated mesh.
         * @see ProcessMesh
         */
        std::vector<face> &GetFaces();

        /**
        * Allows access to nodes of generated mesh.
        * @see ProcessMesh
        */
        std::vector<int> &GetNodes();

    private:
        std::vector<face> mFaces;
        std::vector<int> mNodes;

    };
}

#endif	/* CSHARP_WRAPPER_APPLICATION_MESH_CONVERTER_H_INCLUDED */
