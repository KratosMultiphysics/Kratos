#include "includes/kernel.h"
#include "structural_mechanics_application.h"

namespace KratosWrapper {

	class Interface {
	public:
		static void init(char* path);
		static void updateNodesPos(float* xCoordinates, float* yCoorginates, float* zCoorfinates);
		static void calculate();
		static float* getXCoordinates();
		static float* getYCoordinates();
		static float* getZCoordinates();
		static int getNodesCount();
		static int* getTriangles();
		static int getTrianglesCount();
		
	private:
		static int* pmRequiredNodes;
		static int* pmTriangles;
		static int mTrianglesCount;
		static int mNodesCount;
		static float* pmXCoordinates;
		static float* pmYCoordinates;
		static float* pmZCoordinates;
		static Kratos::Kernel mKernel;
		static Kratos::KratosStructuralMechanicsApplication mApplication;
		static Kratos::ModelPart mMainModelPart;
	};
}