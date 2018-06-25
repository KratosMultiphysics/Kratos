
#include <vector>

namespace CSharpKratosWrapper {

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