#include <Python.h>
#include <iostream>
extern int solve(std::string kratos_path,std::string model_1d_name,std::string model_3d_name);

int main(int argc, char *argv[])
{
	std::string kratos_path("\\home\\rrossi\\kratos\\");
	std::string model_1d_name("aaa");
	std::string model_3d_name("bbb");
	solve(kratos_path,model_1d_name,model_3d_name);
	return 0;
}