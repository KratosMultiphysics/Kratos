#if !defined(CVP_H_INCLUDED )
#define  CVP_H_INCLUDED


#ifdef _WINDLL
#define KRATOS_DLL_IMPORT __declspec(dllimport)
#else
#define KRATOS_DLL_IMPORT
#endif

KRATOS_DLL_IMPORT int solve(
	const char * kratos_path,
	const char * model_1d_name, 
	const char * model_3d_name,
	const char * script_name
	);
#undef KRATOS_DLL_IMPORT

#endif //CVP_H_INCLUDED  defined 

