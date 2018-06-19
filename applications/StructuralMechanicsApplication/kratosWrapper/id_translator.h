#include "mesh_structures.h"
#include <vector>

namespace KratosWrapper {

	class IdTranslator {
	public:
		void init(std::vector<int>& nodes);
		int getUnityId(int kratosId);
		int getKratosId(int unityId);
	private:
		int* pmKratosIds;
		int* pmUnityIds;
	};
}