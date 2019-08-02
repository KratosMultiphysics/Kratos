
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
