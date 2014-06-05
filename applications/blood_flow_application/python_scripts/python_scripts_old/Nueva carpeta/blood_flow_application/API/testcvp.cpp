//#include <iostream>
#include "cvp.h"

//extern int solve(std::string kratos_path,std::string model_1d_name,std::string model_3d_name);

int main(int argc, char *argv[])
{
	if(argc < 3) 
		return 0;

	//const char kratos_path[]="c:\\kratos_win32\\";
	const char model_1d_name[]="aaa";
	const char model_3d_name[]="bbb";
	solve(argv[1],model_1d_name,model_3d_name,argv[2]);
	return 0;
}