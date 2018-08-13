#include <vector>
#include "includes/model_part.h"

namespace CSharpKratosWrapper {

	struct face {
		int nodes[4];
	};

	class MeshConverter {
	public:
		void ProcessMesh(std::vector<Kratos::shared_ptr<Kratos::Element>>& elements);
		std::vector<face>& GetFaces();
		std::vector<int>& GetNodes();
	private:
		std::vector<face> mFaces;
		std::vector<int> mNodes;

	};
}