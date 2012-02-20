/* gidpost 1.7 */
/* this should only be included in gidpostfor.c */

GIDPOST_API
FCALLSCFUN2(INT,GiD_OpenPostMeshFile,GID_OPENPOSTMESHFILE,gid_openpostmeshfile,STRING,INT);

GIDPOST_API
FCALLSCFUN0(INT,GiD_ClosePostMeshFile,GID_CLOSEPOSTMESHFILE,gid_closepostmeshfile);

GIDPOST_API
FCALLSCFUN7(INT,GiD_BeginMeshColor,GID_BEGINMESHCOLOR,gid_beginmeshcolor,STRING,INT,INT,INT,DOUBLE,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN4(INT,GiD_BeginMesh,GID_BEGINMESH,gid_beginmesh,STRING,INT,INT,INT);

GIDPOST_API
FCALLSCFUN0(INT,GiD_EndMesh,GID_ENDMESH,gid_endmesh);

GIDPOST_API
FCALLSCFUN0(INT,GiD_BeginCoordinates,GID_BEGINCOORDINATES,gid_begincoordinates);

GIDPOST_API
FCALLSCFUN0(INT,GiD_EndCoordinates,GID_ENDCOORDINATES,gid_endcoordinates);

GIDPOST_API
FCALLSCFUN1(INT,GiD_BeginMeshGroup,GID_BEGINMESHGROUP,gid_beginmeshgroup,STRING);

GIDPOST_API
FCALLSCFUN0(INT,GiD_EndMeshGroup,GID_ENDMESHGROUP,gid_endmeshgroup);

GIDPOST_API
FCALLSCFUN4(INT,GiD_WriteCoordinates,GID_WRITECOORDINATES,gid_writecoordinates,INT,DOUBLE,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN3(INT,GiD_WriteCoordinates2D,GID_WRITECOORDINATES2D,gid_writecoordinates2d,INT,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN0(INT,GiD_BeginElements,GID_BEGINELEMENTS,gid_beginelements);

GIDPOST_API
FCALLSCFUN0(INT,GiD_EndElements,GID_ENDELEMENTS,gid_endelements);

GIDPOST_API
FCALLSCFUN2(INT,GiD_WriteElement,GID_WRITEELEMENT,gid_writeelement,INT,INTV);

GIDPOST_API
FCALLSCFUN2(INT,GiD_WriteElementMat,GID_WRITEELEMENTMAT, gid_writeelementmat,INT,INTV);

GIDPOST_API
FCALLSCFUN3(INT,GiD_WriteSphere,GID_WRITESPHERE,gid_writesphere,INT,INT,DOUBLE);

GIDPOST_API
FCALLSCFUN4(INT,GiD_WriteSphereMat,GID_WRITESPHEREMAT,gid_writespheremat,INT,INT,DOUBLE,INT);

GIDPOST_API
FCALLSCFUN6(INT,GiD_WriteCircle,GID_WRITECIRCLE,gid_writecircle,INT,INT,DOUBLE,DOUBLE,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN7(INT,GiD_WriteCircleMat,GID_WRITECIRCLEMAT,gid_writecirclemat,INT,INT,DOUBLE,DOUBLE,DOUBLE,DOUBLE,INT);

GIDPOST_API
FCALLSCFUN2(INT,GiD_OpenPostResultFile,GID_OPENPOSTRESULTFILE,gid_openpostresultfile,STRING,INT);

GIDPOST_API
FCALLSCFUN0(INT,GiD_ClosePostResultFile,GID_CLOSEPOSTRESULTFILE,gid_closepostresultfile);

GIDPOST_API
FCALLSCFUN0(INT,GiD_FlushPostFile,GID_FLUSHPOSTFILE,gid_flushpostfile);

GIDPOST_API
FCALLSCFUN6(INT,GiD_BeginGaussPoint,GID_BEGINGAUSSPOINT,gid_begingausspoint,STRING,INT,STRING,INT,INT,INT);

GIDPOST_API
FCALLSCFUN0(INT,GiD_EndGaussPoint,GID_ENDGAUSSPOINT,gid_endgausspoint);

GIDPOST_API
FCALLSCFUN2(INT,GiD_WriteGaussPoint2D,GID_WRITEGAUSSPOINT2D,gid_writegausspoint2d,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN3(INT,GiD_WriteGaussPoint3D,GID_WRITEGAUSSPOINT3D,gid_writegausspoint3d,DOUBLE,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN1(INT,GiD_BeginRangeTable,GID_BEGINRANGETABLE,gid_beginrangetable,STRING);

GIDPOST_API
FCALLSCFUN0(INT,GiD_EndRangeTable,GID_ENDRANGETABLE,gid_endrangetable);

GIDPOST_API
FCALLSCFUN2(INT,GiD_WriteMinRange,GID_WRITEMINRANGE,gid_writeminrange,DOUBLE,STRING);

GIDPOST_API
FCALLSCFUN3(INT,GiD_WriteRange,GID_WRITERANGE,gid_writerange,DOUBLE,DOUBLE,STRING);

GIDPOST_API
FCALLSCFUN2(INT,GiD_WriteMaxRange,GID_WRITEMAXRANGE,gid_writemaxrange,DOUBLE,STRING);

GIDPOST_API
FCALLSCFUN1(INT,GiD_BeginOnMeshGroup,GID_BEGINONMESHGROUP,gid_beginonmeshgroup,STRING);

GIDPOST_API
FCALLSCFUN0(INT,GiD_EndOnMeshGroup,GID_ENDONMESHGROUP,gid_endonmeshgroup);

int GiD_BeginScalarResult(GP_CONST char * Result, GP_CONST char * Analysis, double step,
			  GiD_ResultLocation Where,
			  GP_CONST char * GaussPointsName, GP_CONST char * RangeTable,
			  GP_CONST char * Comp);

GIDPOST_API
FCALLSCFUN7(INT,GiD_BeginScalarResult,GID_BEGINSCALARRESULT,gid_beginscalarresult,STRING,STRING,DOUBLE,INT,STRING,STRING,STRING);

int GiD_BeginVectorResult(GP_CONST char * Result, GP_CONST char * Analysis, double step,
			  GiD_ResultLocation Where,
			  GP_CONST char * GaussPointsName, GP_CONST char * RangeTable,
			  GP_CONST char * Comp1, GP_CONST char * Comp2, GP_CONST char * Comp3,
			  GP_CONST char * Comp4);

GIDPOST_API
FCALLSCFUN10(INT,GiD_BeginVectorResult,GID_BEGINVECTORRESULT,gid_beginvectorresult,STRING,STRING,DOUBLE,INT,STRING,STRING,
	    STRING,STRING,STRING,STRING);

int GiD_Begin2DMatResult(GP_CONST char * Result, GP_CONST char * Analysis, double step,
			 GiD_ResultLocation Where,
			 GP_CONST char * GaussPointsName, GP_CONST char * RangeTable,
			 GP_CONST char * Comp1, GP_CONST char * Comp2,
			 GP_CONST char * Comp3);

GIDPOST_API
FCALLSCFUN9(INT,GiD_Begin2DMatResult,GID_BEGIN2DMATRESULT,gid_begin2dmatresult,STRING,STRING,DOUBLE,INT,STRING,STRING,
	    STRING,STRING,STRING);

int GiD_Begin3DMatResult(GP_CONST char * Result, GP_CONST char * Analysis, double step,
			 GiD_ResultLocation Where,
			 GP_CONST char * GaussPointsName, GP_CONST char * RangeTable,
			 GP_CONST char * Comp1, GP_CONST char * Comp2, GP_CONST char * Comp3,
			 GP_CONST char * Comp4, GP_CONST char * Comp5, GP_CONST char * Comp6);

GIDPOST_API
FCALLSCFUN12(INT,GiD_Begin3DMatResult,GID_BEGIN3DMATRESULT,gid_begin3dmatresult,STRING,STRING,DOUBLE,INT,STRING,STRING,
	    STRING,STRING,STRING,STRING,STRING,STRING);

int GiD_BeginPDMMatResult(GP_CONST char * Result, GP_CONST char * Analysis, double step,
			  GiD_ResultLocation Where,
			  GP_CONST char * GaussPointsName, GP_CONST char * RangeTable,
			  GP_CONST char * Comp1, GP_CONST char * Comp2, GP_CONST char * Comp3,
			  GP_CONST char * Comp4);

GIDPOST_API
FCALLSCFUN10(INT,GiD_BeginPDMMatResult,GID_BEGINPDMMATRESULT,gid_beginpdmmatresult,STRING,STRING,DOUBLE,INT,STRING,STRING,
	     STRING,STRING,STRING,STRING);

int GiD_BeginMainMatResult(GP_CONST char * Result, GP_CONST char * Analysis, double step,
			   GiD_ResultLocation Where,
			   GP_CONST char * GaussPointsName, GP_CONST char * RangeTable,
			   GP_CONST char * Comp1,  GP_CONST char * Comp2,  GP_CONST char * Comp3,
			   GP_CONST char * Comp4,  GP_CONST char * Comp5,  GP_CONST char * Comp6,
			   GP_CONST char * Comp7,  GP_CONST char * Comp8,  GP_CONST char * Comp9,
			   GP_CONST char * Comp10, GP_CONST char * Comp11, GP_CONST char * Comp12);

GIDPOST_API
FCALLSCFUN18(INT,GiD_BeginMainMatResult,GID_BEGINMAINMATRESULT,gid_beginmainmatresult,STRING,STRING,DOUBLE,INT,STRING,STRING,
	     STRING,STRING,STRING,STRING,STRING,STRING,STRING,STRING,STRING,STRING,STRING,STRING);

int GiD_BeginLAResult(GP_CONST char * Result, GP_CONST char * Analysis, double step,
		      GiD_ResultLocation Where,
		      GP_CONST char * GaussPointsName, GP_CONST char * RangeTable,
		      GP_CONST char * Comp1, GP_CONST char * Comp2, GP_CONST char * Comp3);

GIDPOST_API
FCALLSCFUN9(INT,GiD_BeginLAResult,GID_BEGINLARESULT,gid_beginlaresult,STRING,STRING,DOUBLE,INT,STRING,STRING,
	    STRING,STRING,STRING);

GIDPOST_API
FCALLSCFUN6(INT, GiD_BeginResultHeader, GID_BEGINRESULTHEADER, gid_beginresultheader,
	    STRING,  /* Result name */
	    STRING,  /* Analisys name */
	    DOUBLE,  /* step value */
	    INT,     /* result type */
	    INT,     /* result location */
	    STRING); /* location name */

GIDPOST_API
FCALLSCFUN4(INT, GiD_BeginResultGroup, GID_BEGINRESULTGROUP, gid_beginresultgroup,
	    STRING,  /* Analisys name */
	    DOUBLE,  /* step value */
	    INT,     /* result location */
	    STRING); /* location name */

GIDPOST_API
FCALLSCFUN2(INT, GiD_ResultDescription, GID_RESULTDESCRIPTION, gid_resultdescription,
	    STRING,  /* Result name */
	    INT);    /* result type */

GIDPOST_API
FCALLSCFUN3(INT, GiD_ResultDescriptionDim, GID_RESULTDESCRIPTIONDIM, gid_resultdescriptiondim,
	    STRING,  /* Result name */
	    INT,     /* result type */
            INT);    /* result dimension */

GIDPOST_API
FCALLSCFUN0(INT, GiD_ResultValues, GID_RESULTVALUES, gid_resultvalues);

GIDPOST_API
FCALLSCFUN1(INT, GiD_ResultRange, GID_RESULTRANGE, gid_resultrange, STRING);
/* This function is not supported yet inside GiD */
/* FCALLSCFUN1(INT, GiD_ResultUnit, GID_RESULTUNIT, gid_resultunit, STRING); */

int GiD_ScalarComp(GP_CONST char * Comp1);

GIDPOST_API
FCALLSCFUN1(INT, GiD_ScalarComp, GID_SCALARCOMP, gid_scalarcomp, STRING);

int GiD_VectorComp(GP_CONST char * Comp1, GP_CONST char * Comp2,
		   GP_CONST char * Comp3, GP_CONST char * Comp4);

GIDPOST_API
FCALLSCFUN4(INT, GiD_VectorComp, GID_VECTORCOMP, gid_vectorcomp,
	    STRING,  /* component 1 */
	    STRING,  /* component 2 */
	    STRING,  /* component 3 */
	    STRING); /* component 4 (modulus) */

int GiD_2DMatrixComp(GP_CONST char * Comp1, GP_CONST char * Comp2, GP_CONST char * Comp3);
FCALLSCFUN3(INT, GiD_2DMatrixComp, GID_2DMATRIXCOMP, gid_2dmatrixcomp,
	    STRING,  /* component Sxx */
	    STRING,  /* component Syy */
	    STRING); /* component Sxy */

int GiD_3DMatrixComp(GP_CONST char * Comp1, GP_CONST char * Comp2, GP_CONST char * Comp3,
		     GP_CONST char * Comp4, GP_CONST char * Comp5, GP_CONST char * Comp6);

GIDPOST_API
FCALLSCFUN6(INT, GiD_3DMatrixComp, GID_3DMATRIXCOMP, gid_3dmatrixcomp,
	    STRING,  /* component Sxx */
	    STRING,  /* component Syy */
	    STRING,  /* component Szz */
	    STRING,  /* component Sxy */
	    STRING,  /* component Syz */
	    STRING); /* component Sxz */

int GiD_PDMComp(GP_CONST char * Comp1, GP_CONST char * Comp2, GP_CONST char * Comp3, GP_CONST char * Comp4);

GIDPOST_API
FCALLSCFUN4(INT, GiD_PDMComp, GID_PDMCOMP, gid_pdmcomp,
	    STRING,  /* component Sxx */
	    STRING,  /* component Syy */
	    STRING,  /* component Sxy */
	    STRING); /* component Szz */

int GiD_MainMatrixComp(GP_CONST char * Comp1,  GP_CONST char * Comp2,  GP_CONST char * Comp3,
		       GP_CONST char * Comp4,  GP_CONST char * Comp5,  GP_CONST char * Comp6,
		       GP_CONST char * Comp7,  GP_CONST char * Comp8,  GP_CONST char * Comp9,
		       GP_CONST char * Comp10, GP_CONST char * Comp11, GP_CONST char * Comp12);

GIDPOST_API
FCALLSCFUN12(INT, GiD_MainMatrixComp, GID_MAINMATRIXCOMP, gid_mainmatrixcomp,
	     STRING,   /* component Si */
	     STRING,   /* component Sii */
	     STRING,   /* component Siii */
	     STRING,   /* component Vix */
	     STRING,   /* component Viy */
	     STRING,   /* component Viz */
	     STRING,   /* component Viix */
	     STRING,   /* component Viiy */
	     STRING,   /* component Viiz */
	     STRING,   /* component Viiix */
	     STRING,   /* component Viiiy */
	     STRING);  /* component Viiiz */

int GiD_LAComponents(GP_CONST char * Comp1, GP_CONST char * Comp2, GP_CONST char * Comp3);

GIDPOST_API
FCALLSCFUN3(INT, GiD_LAComponents, GID_LACOMPONENTS, gid_lacomponents,
	    STRING,  /* local axes name 1 */
	    STRING,  /* local axes name 2 */
	    STRING); /* local axes name 3 */

GIDPOST_API
FCALLSCFUN0(INT,GiD_EndResult,GID_ENDRESULT,gid_endresult);

GIDPOST_API
FCALLSCFUN2(INT, GiD_WriteScalar,GID_WRITESCALAR,gid_writescalar,INT,DOUBLE);

GIDPOST_API
FCALLSCFUN3(INT, GiD_Write2DVector,GID_WRITE2DVECTOR,gid_write2dvector,INT,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN4(INT, GiD_WriteVector,GID_WRITEVECTOR,gid_writevector,INT,DOUBLE,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN5(INT, GiD_WriteVectorModule,GID_WRITEVECTORMODULE,gid_writevectormodule,INT,DOUBLE,DOUBLE,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN4(INT, GiD_Write2DMatrix,GID_WRITE2DMATRIX,gid_write2dmatrix,INT,DOUBLE,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN7(INT, GiD_Write3DMatrix,GID_WRITE3DMATRIX,gid_write3dmatrix,INT,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN5(INT, GiD_WritePlainDefMatrix,GID_WRITEPLAINDEFMATRIX,gid_writeplaindefmatrix,INT,DOUBLE,DOUBLE,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN13(INT, GiD_WriteMainMatrix,GID_WRITEMAINMATRIX,gid_writemainmatrix,INT,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN4(INT, GiD_WriteLocalAxes,GID_WRITELOCALAXES,gid_writelocalaxes,INT,DOUBLE,DOUBLE,DOUBLE);
