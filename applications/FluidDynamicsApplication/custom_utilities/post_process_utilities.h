//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Angel Celigueta


#ifndef POST_PROCESS_UTILITIES_H
#define POST_PROCESS_UTILITIES_H

#include "includes/define.h"
#include "includes/variables.h"

namespace Kratos {

    class PostProcessUtilities {

    public:

        KRATOS_CLASS_POINTER_DEFINITION(PostProcessUtilities);

        /// Default constructor.

        PostProcessUtilities() {};

        /// Destructor.

        virtual ~PostProcessUtilities() {};

        void ComputeFlow(const ModelPart& rModelPart, double& flow);

    protected:

    private:


    }; // Class PostProcessUtilities

} // namespace Kratos.

#endif // POST_PROCESS_UTILITIES_H
