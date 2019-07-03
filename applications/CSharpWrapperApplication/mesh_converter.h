#include <vector>
#include "includes/model_part.h"

namespace CSharpKratosWrapper {

	struct face {
		int nodes[4];
	};

	class MeshConverter {
	public:
		void ProcessMesh(Kratos::ModelPart::ElementsContainerType::ContainerType &elements);

		std::vector<face> &GetFaces();

		std::vector<int> &GetNodes();

	private:
		std::vector<face> mFaces;
		std::vector<int> mNodes;

	};
}