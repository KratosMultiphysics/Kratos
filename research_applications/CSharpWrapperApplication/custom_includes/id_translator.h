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

#if !defined(CSHARP_WRAPPER_APPLICATION_ID_TRANSLATOR_H_INCLUDED )
#define  CSHARP_WRAPPER_APPLICATION_ID_TRANSLATOR_H_INCLUDED

// System includes

// External includes

// Project includes
#include <vector>

namespace CSharpKratosWrapper {

    class IdTranslator {
    public:
        /**
         * Initializes translator.
         * @param nodes Sorted vector of Kratos ids of surface nodes
         */
        void init(std::vector<int> &nodes);

        /**
         * Translates Kratos node id to surface id
         * @param kratosId Kratos node id. Must be on the surface.
         * @return surface id of this node
         * @see safeGetSurfaceId
         */
        int getSurfaceId(int kratosId);

        /**
         * Translates surface id to Kratos id.
         * @return Kratos node id
         * @see safeGetKratosId
         */
        int getKratosId(int surfaceId);

        /**
         *  Checks whether it is possible to translate given Kratos id.
         *  It is true if and only if corresponding node is on the surface.
         * @param kratosId Kratos id to check
         * @return true if corresponding surface Id exists
         */
        bool hasSurfaceId(int kratosId);

        /**
         *  Checks whether it is possible to translate given surface id.
         * @param surfaceId surface id to check
         * @return true if corresponding node exists
         */
        bool hasKratosId(int surfaceId);

        /**
         * Performs safety check and translates given Kratos id to surface id. Returns -1 if translation is not possible.
         * @param kratosId Kratos id to translate
         * @return surface node id or -1
         */
        int safeGetSurfaceId(int kratosId);

        /**
         * Performs safety check and translates given surface id to Kratos node id. Returns -1 if translation is not possible.
         * @param surfaceId surface id to translate
         * @return Kratos node id or -1
         */
        int safeGetKratosId(int surfaceId);

    private:
        int *pmKratosIds;
        int *pmSurfaceIds;
        int mSurfaceIdsSize;
        int mKratosIdsSize;
    };
}

#endif	/* CSHARP_WRAPPER_APPLICATION_ID_TRANSLATOR_H_INCLUDED */
