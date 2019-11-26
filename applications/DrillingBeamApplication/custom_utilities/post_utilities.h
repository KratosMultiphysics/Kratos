#ifndef POST_UTILITIES_H
#define POST_UTILITIES_H

#include "includes/define.h"
#include "includes/variables.h"
#include "includes/model_part.h"



namespace Kratos {


    class PostUtilities {

    public:

        typedef ModelPart::ElementsContainerType ElementsArrayType;
        typedef ModelPart::NodesContainerType NodesArrayType;

        KRATOS_CLASS_POINTER_DEFINITION(PostUtilities);

        /// Default constructor.

        PostUtilities() {};

        /// Destructor.

        virtual ~PostUtilities() {};

     
    void CreateSkinForBeam(ModelPart& r_model_part){

        // typedef ModelPart::NodesContainerType NodesArrayType;
        // NodesArrayType& pNodes = r_model_part.GetCommunicator().LocalMesh().Nodes();

        //#pragma omp parallel for
        // for (int k = 0; k < (int) pNodes.size(); k++) {

        //     ModelPart::NodeIterator i_iterator = pNodes.ptr_begin() + k;
        //     Node < 3 > & i = *i_iterator;

        //     //const double& displacement = i.FastGetSolutionStepValue(DISPLACEMENT);
        //     const array_1d<double, 3>& coords = i.Coordinates();

        //     KRATOS_WATCH(coords)
            
        // }//for Node
    }

    }; // Class PostUtilities

} // namespace Kratos.

#endif // POST_UTILITIES_H