#include <vector>
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
