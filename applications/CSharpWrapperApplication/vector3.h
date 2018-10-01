
namespace CSharpKratosWrapper {

	class Vector3 {
	public:
		double x;
		double y;
		double z;
		Vector3(double x, double y, double z);
		Vector3();
		Vector3* cross(Vector3&);
		double dot(Vector3&);
		Vector3* add(Vector3&);
		Vector3* add(double otherX, double otherY, double otherZ);
		Vector3* sub(Vector3&);
		Vector3* sub(double otherX, double otherY, double otherZ);
	};

}