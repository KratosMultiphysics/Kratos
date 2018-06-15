#include "includes/kernel.h"
#include "structural_mechanics_application.h"
#include "mesh_converter.h"

namespace KratosWrapper {

	class Interface {
	public:
		static void init(char* mdpaPath);
		static void updateNodePos(int nodeId, float x, float y, float z);
		static void calculate();
		static float* getXCoordinates();
		static float* getYCoordinates();
		static float* getZCoordinates();
		static int getNodesCount();
		static int* getTriangles();
		static int getTrianglesCount();
		
	private:
		static int* pmKratosNodeIds;
		static int* pmTriangles;
		static int mTrianglesCount;
		static int mNodesCount;
		static float* pmXCoordinates;
		static float* pmYCoordinates;
		static float* pmZCoordinates;
		static Kratos::Kernel mKernel;
		static Kratos::KratosStructuralMechanicsApplication mApplication;
		static Kratos::ModelPart mMainModelPart;

		static void initInternals();
		static void loadMDPAFile(char* mdpaPath);
		static void saveTriangles(MeshConverter& meshConverter);
		static void saveNodes(MeshConverter& meshConverter);
	};
}