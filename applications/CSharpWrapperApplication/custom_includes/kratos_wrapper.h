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

#ifndef KRATOSMULTIPHYSICS_KRATOS_WRAPPER_H
#define KRATOSMULTIPHYSICS_KRATOS_WRAPPER_H

// System includes

// External includes

// Project includes
#include "kratos_internals.h"
#include "includes/node.h"
#include "model_part_wrapper.h"

namespace CSharpKratosWrapper {

    using NodeType = Kratos::Node<3>;
    using ModelPart = Kratos::ModelPart;

    class KratosWrapper {

    public:

        ~KratosWrapper();

        void init(const char *MDPAFilePath, const char *JSONFilePath = NULL);

        void initWithSettings(const char *JSONFilePath = NULL);

        void calculate();

        ModelPartWrapper* getRootModelPartWrapper();

    private:
        KratosInternals mKratosInternals;
        std::vector<NodeType::Pointer> mFixedNodes;
        ModelPartWrapper* pmMainModelPartWrapper;

        void freeNodes();

    };
}

#endif //KRATOSMULTIPHYSICS_KRATOS_WRAPPER_H
