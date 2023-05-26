//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes

// External includes

// Project includes
#include "input_output/base_64_encoded_output.h"

namespace Kratos 
{

Base64EncodedOutput::~Base64EncodedOutput()
{
    const IndexType padding = 3 - mByteTripletIndex;

    // Check if there are remaining bytes and padding needs to be added
    if (padding != 0 && mByteTripletIndex != 0) {
        // Fill the remaining bytes with null characters
        for (; mByteTripletIndex < 3; ++mByteTripletIndex) {
            mByteTriplet[mByteTripletIndex] = '\0';
        }

        // Encode the remaining byte triplet and write it to the output stream
        EncodeTriplet(mrOStream, mByteTriplet, padding);
    }

    mByteTripletIndex = 0; // Reset the byte triplet index for future encoding operations
}

} // namespace Kratos