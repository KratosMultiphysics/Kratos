/* gidpost 2.11 */
/* this should only be included in gidpostfor.c */

GIDPOST_API
FCALLSCFUN0(INT,GiD_PostInit,GID_POSTINIT,gid_postinit);

GIDPOST_API
FCALLSCFUN0(INT,GiD_PostDone,GID_POSTDONE,gid_postdone);

GIDPOST_API
FCALLSCFUN2(INT,GiD_OpenPostMeshFile,GID_OPENPOSTMESHFILE,gid_openpostmeshfile,STRING,INT);

GIDPOST_API
FCALLSCFUN2(INT,GiD_fOpenPostMeshFile,GID_FOPENPOSTMESHFILE,gid_fopenpostmeshfile,STRING,INT);

GIDPOST_API
FCALLSCFUN0(INT,GiD_ClosePostMeshFile,GID_CLOSEPOSTMESHFILE,gid_closepostmeshfile);

GIDPOST_API
FCALLSCFUN1(INT,GiD_fClosePostMeshFile,GID_FCLOSEPOSTMESHFILE,gid_fclosepostmeshfile,INT);

GIDPOST_API
FCALLSCFUN4(INT,GiD_BeginMesh,GID_BEGINMESH,gid_beginmesh,STRING,INT,INT,INT);

GIDPOST_API
FCALLSCFUN5(INT,GiD_fBeginMesh,GID_FBEGINMESH,gid_fbeginmesh,INT,STRING,INT,INT,INT);

GIDPOST_API
FCALLSCFUN7(INT,GiD_BeginMeshColor,GID_BEGINMESHCOLOR,gid_beginmeshcolor,STRING,INT,INT,INT,DOUBLE,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN8(INT,GiD_fBeginMeshColor,GID_FBEGINMESHCOLOR,gid_fbeginmeshcolor,INT,STRING,INT,INT,INT,DOUBLE,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN0(INT,GiD_EndMesh,GID_ENDMESH,gid_endmesh);

GIDPOST_API
FCALLSCFUN1(INT,GiD_fEndMesh,GID_FENDMESH,gid_fendmesh,INT);

GIDPOST_API
FCALLSCFUN1(INT,GiD_MeshUnit,GID_MESHUNIT,gid_meshunit,STRING);

GIDPOST_API
FCALLSCFUN2(INT,GiD_fMeshUnit,GID_FMESHUNIT,gid_fmeshunit,INT,STRING);

GIDPOST_API
FCALLSCFUN0(INT,GiD_BeginCoordinates,GID_BEGINCOORDINATES,gid_begincoordinates);

GIDPOST_API
FCALLSCFUN1(INT,GiD_fBeginCoordinates,GID_FBEGINCOORDINATES,gid_fbegincoordinates,INT);

GIDPOST_API
FCALLSCFUN0(INT,GiD_EndCoordinates,GID_ENDCOORDINATES,gid_endcoordinates);

GIDPOST_API
FCALLSCFUN1(INT,GiD_fEndCoordinates,GID_FENDCOORDINATES,gid_fendcoordinates,INT);

GIDPOST_API
FCALLSCFUN1(INT,GiD_BeginMeshGroup,GID_BEGINMESHGROUP,gid_beginmeshgroup,STRING);

GIDPOST_API
FCALLSCFUN2(INT,GiD_fBeginMeshGroup,GID_FBEGINMESHGROUP,gid_fbeginmeshgroup,INT,STRING);

GIDPOST_API
FCALLSCFUN0(INT,GiD_EndMeshGroup,GID_ENDMESHGROUP,gid_endmeshgroup);

GIDPOST_API
FCALLSCFUN1(INT,GiD_fEndMeshGroup,GID_FENDMESHGROUP,gid_fendmeshgroup,INT);

GIDPOST_API
FCALLSCFUN4(INT,GiD_WriteCoordinates,GID_WRITECOORDINATES,gid_writecoordinates,INT,DOUBLE,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN5(INT,GiD_fWriteCoordinates,GID_FWRITECOORDINATES,gid_fwritecoordinates,INT,INT,DOUBLE,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN3(INT,GiD_WriteCoordinates2D,GID_WRITECOORDINATES2D,gid_writecoordinates2d,INT,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN4(INT,GiD_fWriteCoordinates2D,GID_FWRITECOORDINATES2D,gid_fwritecoordinates2d,INT,INT,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN0(INT,GiD_BeginElements,GID_BEGINELEMENTS,gid_beginelements);

GIDPOST_API
FCALLSCFUN1(INT,GiD_fBeginElements,GID_FBEGINELEMENTS,gid_fbeginelements,INT);

GIDPOST_API
FCALLSCFUN0(INT,GiD_EndElements,GID_ENDELEMENTS,gid_endelements);

GIDPOST_API
FCALLSCFUN1(INT,GiD_fEndElements,GID_FENDELEMENTS,gid_fendelements,INT);

GIDPOST_API
FCALLSCFUN2(INT,GiD_WriteElement,GID_WRITEELEMENT,gid_writeelement,INT,INTV);

GIDPOST_API
FCALLSCFUN3(INT,GiD_fWriteElement,GID_FWRITEELEMENT,gid_fwriteelement,INT,INT,INTV);

GIDPOST_API
FCALLSCFUN2(INT,GiD_WriteElementMat,GID_WRITEELEMENTMAT, gid_writeelementmat,INT,INTV);

GIDPOST_API
FCALLSCFUN3(INT,GiD_fWriteElementMat,GID_FWRITEELEMENTMAT, gid_fwriteelementmat,INT,INT,INTV);

GIDPOST_API
FCALLSCFUN3(INT,GiD_WriteSphere,GID_WRITESPHERE,gid_writesphere,INT,INT,DOUBLE);

GIDPOST_API
FCALLSCFUN4(INT,GiD_fWriteSphere,GID_FWRITESPHERE,gid_fwritesphere,INT,INT,INT,DOUBLE);

GIDPOST_API
FCALLSCFUN4(INT,GiD_WriteSphereMat,GID_WRITESPHEREMAT,gid_writespheremat,INT,INT,DOUBLE,INT);

GIDPOST_API
FCALLSCFUN5(INT,GiD_fWriteSphereMat,GID_FWRITESPHEREMAT,gid_fwritespheremat,INT,INT,INT,DOUBLE,INT);

GIDPOST_API
FCALLSCFUN6(INT,GiD_WriteCircle,GID_WRITECIRCLE,gid_writecircle,INT,INT,DOUBLE,DOUBLE,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN7(INT,GiD_fWriteCircle,GID_FWRITECIRCLE,gid_fwritecircle,INT,INT,INT,DOUBLE,DOUBLE,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN7(INT,GiD_WriteCircleMat,GID_WRITECIRCLEMAT,gid_writecirclemat,INT,INT,DOUBLE,DOUBLE,DOUBLE,DOUBLE,INT);

GIDPOST_API
FCALLSCFUN8(INT,GiD_fWriteCircleMat,GID_FWRITECIRCLEMAT,gid_fwritecirclemat,INT,INT,INT,DOUBLE,DOUBLE,DOUBLE,DOUBLE,INT);

GIDPOST_API
FCALLSCFUN2(INT,GiD_OpenPostResultFile,GID_OPENPOSTRESULTFILE,gid_openpostresultfile,STRING,INT);

GIDPOST_API
FCALLSCFUN2(INT,GiD_fOpenPostResultFile,GID_FOPENPOSTRESULTFILE,gid_fopenpostresultfile,STRING,INT);

GIDPOST_API
FCALLSCFUN0(INT,GiD_ClosePostResultFile,GID_CLOSEPOSTRESULTFILE,gid_closepostresultfile);

GIDPOST_API
FCALLSCFUN1(INT,GiD_fClosePostResultFile,GID_FCLOSEPOSTRESULTFILE,gid_fclosepostresultfile,INT);

GIDPOST_API
FCALLSCFUN0(INT,GiD_FlushPostFile,GID_FLUSHPOSTFILE,gid_flushpostfile);

GIDPOST_API
FCALLSCFUN1(INT,GiD_fFlushPostFile,GID_FFLUSHPOSTFILE,gid_fflushpostfile,INT);

GIDPOST_API
FCALLSCFUN6(INT,GiD_BeginGaussPoint,GID_BEGINGAUSSPOINT,gid_begingausspoint,STRING,INT,STRING,INT,INT,INT);

GIDPOST_API
FCALLSCFUN7(INT,GiD_fBeginGaussPoint,GID_FBEGINGAUSSPOINT,gid_fbegingausspoint,INT,STRING,INT,STRING,INT,INT,INT);

GIDPOST_API
FCALLSCFUN0(INT,GiD_EndGaussPoint,GID_ENDGAUSSPOINT,gid_endgausspoint);

GIDPOST_API
FCALLSCFUN1(INT,GiD_fEndGaussPoint,GID_FENDGAUSSPOINT,gid_fendgausspoint,INT);

GIDPOST_API
FCALLSCFUN2(INT,GiD_WriteGaussPoint2D,GID_WRITEGAUSSPOINT2D,gid_writegausspoint2d,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN3(INT,GiD_fWriteGaussPoint2D,GID_FWRITEGAUSSPOINT2D,gid_fwritegausspoint2d,INT,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN3(INT,GiD_WriteGaussPoint3D,GID_WRITEGAUSSPOINT3D,gid_writegausspoint3d,DOUBLE,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN4(INT,GiD_fWriteGaussPoint3D,GID_FWRITEGAUSSPOINT3D,gid_fwritegausspoint3d,INT,DOUBLE,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN1(INT,GiD_BeginRangeTable,GID_BEGINRANGETABLE,gid_beginrangetable,STRING);

GIDPOST_API
FCALLSCFUN2(INT,GiD_fBeginRangeTable,GID_FBEGINRANGETABLE,gid_fbeginrangetable,INT,STRING);

GIDPOST_API
FCALLSCFUN0(INT,GiD_EndRangeTable,GID_ENDRANGETABLE,gid_endrangetable);

GIDPOST_API
FCALLSCFUN1(INT,GiD_fEndRangeTable,GID_FENDRANGETABLE,gid_fendrangetable,INT);

GIDPOST_API
FCALLSCFUN2(INT,GiD_WriteMinRange,GID_WRITEMINRANGE,gid_writeminrange,DOUBLE,STRING);

GIDPOST_API
FCALLSCFUN3(INT,GiD_fWriteMinRange,GID_FWRITEMINRANGE,gid_fwriteminrange,INT,DOUBLE,STRING);

GIDPOST_API
FCALLSCFUN3(INT,GiD_WriteRange,GID_WRITERANGE,gid_writerange,DOUBLE,DOUBLE,STRING);

GIDPOST_API
FCALLSCFUN4(INT,GiD_fWriteRange,GID_FWRITERANGE,gid_fwriterange,INT,DOUBLE,DOUBLE,STRING);

GIDPOST_API
FCALLSCFUN2(INT,GiD_WriteMaxRange,GID_WRITEMAXRANGE,gid_writemaxrange,DOUBLE,STRING);

GIDPOST_API
FCALLSCFUN3(INT,GiD_fWriteMaxRange,GID_FWRITEMAXRANGE,gid_fwritemaxrange,INT,DOUBLE,STRING);

GIDPOST_API
FCALLSCFUN1(INT,GiD_BeginOnMeshGroup,GID_BEGINONMESHGROUP,gid_beginonmeshgroup,STRING);

GIDPOST_API
FCALLSCFUN2(INT,GiD_fBeginOnMeshGroup,GID_FBEGINONMESHGROUP,gid_fbeginonmeshgroup,INT,STRING);

GIDPOST_API
FCALLSCFUN0(INT,GiD_EndOnMeshGroup,GID_ENDONMESHGROUP,gid_endonmeshgroup);

GIDPOST_API
FCALLSCFUN1(INT,GiD_fEndOnMeshGroup,GID_FENDONMESHGROUP,gid_fendonmeshgroup,INT);

int GiD_BeginScalarResult(GP_CONST char * Result, GP_CONST char * Analysis, double step,
			  GiD_ResultLocation Where,
			  GP_CONST char * GaussPointsName, GP_CONST char * RangeTable,
			  GP_CONST char * Comp);

int GiD_fBeginScalarResult(GiD_FILE fd,
                           GP_CONST char * Result, GP_CONST char * Analysis, double step,
                           GiD_ResultLocation Where,
                           GP_CONST char * GaussPointsName, GP_CONST char * RangeTable,
                           GP_CONST char * Comp);

GIDPOST_API
FCALLSCFUN7(INT,GiD_BeginScalarResult,GID_BEGINSCALARRESULT,gid_beginscalarresult,
            STRING,STRING,DOUBLE,INT,STRING,STRING,STRING);

GIDPOST_API
FCALLSCFUN8(INT,GiD_fBeginScalarResult,GID_FBEGINSCALARRESULT,gid_fbeginscalarresult,
            INT,STRING,STRING,DOUBLE,INT,STRING,STRING,STRING);

int GiD_BeginVectorResult(GP_CONST char * Result, GP_CONST char * Analysis, double step,
			  GiD_ResultLocation Where,
			  GP_CONST char * GaussPointsName, GP_CONST char * RangeTable,
			  GP_CONST char * Comp1, GP_CONST char * Comp2, GP_CONST char * Comp3,
			  GP_CONST char * Comp4);

int GiD_fBeginVectorResult(GiD_FILE fd,
                           GP_CONST char * Result, GP_CONST char * Analysis, double step,
                           GiD_ResultLocation Where,
                           GP_CONST char * GaussPointsName, GP_CONST char * RangeTable,
                           GP_CONST char * Comp1, GP_CONST char * Comp2, GP_CONST char * Comp3,
                           GP_CONST char * Comp4);

GIDPOST_API
FCALLSCFUN10(INT,GiD_BeginVectorResult,GID_BEGINVECTORRESULT,gid_beginvectorresult,
             STRING,STRING,DOUBLE,INT,STRING,STRING,STRING,STRING,STRING,STRING);

GIDPOST_API
FCALLSCFUN11(INT,GiD_fBeginVectorResult,GID_FBEGINVECTORRESULT,gid_fbeginvectorresult,
             INT,STRING,STRING,DOUBLE,INT,STRING,STRING,STRING,STRING,STRING,STRING);

int GiD_Begin2DMatResult(GP_CONST char * Result, GP_CONST char * Analysis, double step,
			 GiD_ResultLocation Where,
			 GP_CONST char * GaussPointsName, GP_CONST char * RangeTable,
			 GP_CONST char * Comp1, GP_CONST char * Comp2,
			 GP_CONST char * Comp3);

int GiD_fBegin2DMatResult(GiD_FILE fd,
                          GP_CONST char * Result, GP_CONST char * Analysis, double step,
                          GiD_ResultLocation Where,
                          GP_CONST char * GaussPointsName, GP_CONST char * RangeTable,
                          GP_CONST char * Comp1, GP_CONST char * Comp2,
                          GP_CONST char * Comp3);

GIDPOST_API
FCALLSCFUN9(INT,GiD_Begin2DMatResult,GID_BEGIN2DMATRESULT,gid_begin2dmatresult,STRING,STRING,DOUBLE,INT,STRING,STRING,
	    STRING,STRING,STRING);

GIDPOST_API
FCALLSCFUN10(INT,GiD_fBegin2DMatResult,GID_FBEGIN2DMATRESULT,gid_fbegin2dmatresult,
             INT,STRING,STRING,DOUBLE,INT,STRING,STRING,STRING,STRING,STRING);

int GiD_Begin3DMatResult(GP_CONST char * Result, GP_CONST char * Analysis, double step,
			 GiD_ResultLocation Where,
			 GP_CONST char * GaussPointsName, GP_CONST char * RangeTable,
			 GP_CONST char * Comp1, GP_CONST char * Comp2, GP_CONST char * Comp3,
			 GP_CONST char * Comp4, GP_CONST char * Comp5, GP_CONST char * Comp6);

int GiD_fBegin3DMatResult(GiD_FILE fd,
                          GP_CONST char * Result, GP_CONST char * Analysis, double step,
                          GiD_ResultLocation Where,
                          GP_CONST char * GaussPointsName, GP_CONST char * RangeTable,
                          GP_CONST char * Comp1, GP_CONST char * Comp2, GP_CONST char * Comp3,
                          GP_CONST char * Comp4, GP_CONST char * Comp5, GP_CONST char * Comp6);

GIDPOST_API
FCALLSCFUN12(INT,GiD_Begin3DMatResult,GID_BEGIN3DMATRESULT,gid_begin3dmatresult,STRING,STRING,DOUBLE,INT,STRING,STRING,
	    STRING,STRING,STRING,STRING,STRING,STRING);

GIDPOST_API
FCALLSCFUN13(INT,GiD_fBegin3DMatResult,GID_FBEGIN3DMATRESULT,gid_fbegin3dmatresult,
             INT,STRING,STRING,DOUBLE,INT,STRING,STRING,STRING,STRING,STRING,STRING,STRING,STRING);

int GiD_BeginPDMMatResult(GP_CONST char * Result, GP_CONST char * Analysis, double step,
			  GiD_ResultLocation Where,
			  GP_CONST char * GaussPointsName, GP_CONST char * RangeTable,
			  GP_CONST char * Comp1, GP_CONST char * Comp2, GP_CONST char * Comp3,
			  GP_CONST char * Comp4);

int GiD_fBeginPDMMatResult(GiD_FILE fd,
                           GP_CONST char * Result, GP_CONST char * Analysis, double step,
                           GiD_ResultLocation Where,
                           GP_CONST char * GaussPointsName, GP_CONST char * RangeTable,
                           GP_CONST char * Comp1, GP_CONST char * Comp2, GP_CONST char * Comp3,
                           GP_CONST char * Comp4);

GIDPOST_API
FCALLSCFUN10(INT,GiD_BeginPDMMatResult,GID_BEGINPDMMATRESULT,gid_beginpdmmatresult,
             STRING,STRING,DOUBLE,INT,STRING,STRING,STRING,STRING,STRING,STRING);

GIDPOST_API
FCALLSCFUN11(INT,GiD_fBeginPDMMatResult,FGID_BEGINPDMMATRESULT,gid_fbeginpdmmatresult,
             INT,STRING,STRING,DOUBLE,INT,STRING,STRING,STRING,STRING,STRING,STRING);

int GiD_BeginMainMatResult(GP_CONST char * Result, GP_CONST char * Analysis, double step,
			   GiD_ResultLocation Where,
			   GP_CONST char * GaussPointsName, GP_CONST char * RangeTable,
			   GP_CONST char * Comp1,  GP_CONST char * Comp2,  GP_CONST char * Comp3,
			   GP_CONST char * Comp4,  GP_CONST char * Comp5,  GP_CONST char * Comp6,
			   GP_CONST char * Comp7,  GP_CONST char * Comp8,  GP_CONST char * Comp9,
			   GP_CONST char * Comp10, GP_CONST char * Comp11, GP_CONST char * Comp12);

int GiD_fBeginMainMatResult(GiD_FILE fd,
                            GP_CONST char * Result, GP_CONST char * Analysis, double step,
                            GiD_ResultLocation Where,
                            GP_CONST char * GaussPointsName, GP_CONST char * RangeTable,
                            GP_CONST char * Comp1,  GP_CONST char * Comp2,  GP_CONST char * Comp3,
                            GP_CONST char * Comp4,  GP_CONST char * Comp5,  GP_CONST char * Comp6,
                            GP_CONST char * Comp7,  GP_CONST char * Comp8,  GP_CONST char * Comp9,
                            GP_CONST char * Comp10, GP_CONST char * Comp11, GP_CONST char * Comp12);

GIDPOST_API
FCALLSCFUN18(INT,GiD_BeginMainMatResult,GID_BEGINMAINMATRESULT,gid_beginmainmatresult,STRING,STRING,DOUBLE,INT,STRING,STRING,
	     STRING,STRING,STRING,STRING,STRING,STRING,STRING,STRING,STRING,STRING,STRING,STRING);

GIDPOST_API
FCALLSCFUN19(INT,GiD_fBeginMainMatResult,GID_FBEGINMAINMATRESULT,gid_fbeginmainmatresult,
             INT,STRING,STRING,DOUBLE,INT,STRING,STRING,STRING,STRING,STRING,STRING,STRING,STRING,STRING,STRING,STRING,STRING,STRING,STRING);

int GiD_BeginLAResult(GP_CONST char * Result, GP_CONST char * Analysis, double step,
		      GiD_ResultLocation Where,
		      GP_CONST char * GaussPointsName, GP_CONST char * RangeTable,
		      GP_CONST char * Comp1, GP_CONST char * Comp2, GP_CONST char * Comp3);

int GiD_fBeginLAResult(GiD_FILE fd,
                       GP_CONST char * Result, GP_CONST char * Analysis, double step,
                       GiD_ResultLocation Where,
                       GP_CONST char * GaussPointsName, GP_CONST char * RangeTable,
                       GP_CONST char * Comp1, GP_CONST char * Comp2, GP_CONST char * Comp3);

GIDPOST_API
FCALLSCFUN9(INT,GiD_BeginLAResult,GID_BEGINLARESULT,gid_beginlaresult,STRING,STRING,DOUBLE,INT,STRING,STRING,
	    STRING,STRING,STRING);

GIDPOST_API
FCALLSCFUN10(INT,GiD_fBeginLAResult,GID_FBEGINLARESULT,gid_fbeginlaresult,
             INT,STRING,STRING,DOUBLE,INT,STRING,STRING,STRING,STRING,STRING);

int GiD_BeginComplexScalarResult(GP_CONST char * Result, GP_CONST char * Analysis, double step,
		          GiD_ResultLocation Where,
		          GP_CONST char * GaussPointsName, GP_CONST char * RangeTable,
		          GP_CONST char * Re, GP_CONST char * Im);

int GiD_fBeginComplexScalarResult(GiD_FILE fd, GP_CONST char * Result, GP_CONST char * Analysis, double step,
		           GiD_ResultLocation Where,
		           GP_CONST char * GaussPointsName, GP_CONST char * RangeTable,
		           GP_CONST char * Re, GP_CONST char * Im);

GIDPOST_API
FCALLSCFUN8(INT,GiD_BeginComplexScalarResult,GID_BEGINCOMPLEXSCALARRESULT,gid_begincomplexscalarresult,STRING,STRING,DOUBLE,INT,STRING,STRING,
	    STRING,STRING);

GIDPOST_API
FCALLSCFUN9(INT,GiD_fBeginComplexScalarResult,GID_FBEGINCOMPLEXSCALARRESULT,gid_fbegincomplexscalarresult,
             INT,STRING,STRING,DOUBLE,INT,STRING,STRING,STRING,STRING);

int GiD_BeginComplexVectorResult(GP_CONST char * Result, GP_CONST char * Analysis, double step,
		          GiD_ResultLocation Where,
		          GP_CONST char * GaussPointsName, GP_CONST char * RangeTable,
		          GP_CONST char * Rex, GP_CONST char * Imx,
                          GP_CONST char * Rey, GP_CONST char * Imy,
                          GP_CONST char * Rez, GP_CONST char * Imz);

int GiD_fBeginComplexVectorResult(GiD_FILE fd, GP_CONST char * Result, GP_CONST char * Analysis, double step,
		          GiD_ResultLocation Where,
		          GP_CONST char * GaussPointsName, GP_CONST char * RangeTable,
		          GP_CONST char * Rex, GP_CONST char * Imx,
                          GP_CONST char * Rey, GP_CONST char * Imy,
                          GP_CONST char * Rez, GP_CONST char * Imz);

GIDPOST_API
FCALLSCFUN12(INT,GiD_BeginComplexVectorResult,GID_BEGINCOMPLEXVECTRORESULT,gid_begincomplexvectorresult,STRING,STRING,DOUBLE,INT,STRING,STRING,
	    STRING,STRING,STRING,STRING,STRING,STRING);

GIDPOST_API
FCALLSCFUN13(INT,GiD_fBeginComplexVectorResult,GID_FBEGINCOMPLEXVECTORRESULT,gid_fbegincomplexvectorresult,
             INT,STRING,STRING,DOUBLE,INT,STRING,STRING,STRING,STRING,STRING,STRING,STRING,STRING);


GIDPOST_API
FCALLSCFUN6(INT, GiD_BeginResultHeader, GID_BEGINRESULTHEADER, gid_beginresultheader,
	    STRING,  /* Result name */
	    STRING,  /* Analisys name */
	    DOUBLE,  /* step value */
	    INT,     /* result type */
	    INT,     /* result location */
	    STRING); /* location name */

GIDPOST_API
FCALLSCFUN7(INT, GiD_fBeginResultHeader, GID_FBEGINRESULTHEADER, gid_fbeginresultheader,
            INT,     /* FILE handler */
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
FCALLSCFUN5(INT, GiD_fBeginResultGroup, GID_FBEGINRESULTGROUP, gid_fbeginresultgroup,
            INT,     /* FILE handler */
	    STRING,  /* Analisys name */
	    DOUBLE,  /* step value */
	    INT,     /* result location */
	    STRING); /* location name */

GIDPOST_API
FCALLSCFUN2(INT, GiD_ResultDescription, GID_RESULTDESCRIPTION, gid_resultdescription,
	    STRING,  /* Result name */
	    INT);    /* result type */

GIDPOST_API
FCALLSCFUN3(INT, GiD_fResultDescription, GID_FRESULTDESCRIPTION, gid_fresultdescription,
            INT,     /* FILE handler */
	    STRING,  /* Result name */
	    INT);    /* result type */

GIDPOST_API
FCALLSCFUN3(INT, GiD_ResultDescriptionDim, GID_RESULTDESCRIPTIONDIM, gid_resultdescriptiondim,
	    STRING,  /* Result name */
	    INT,     /* result type */
            INT);    /* result dimension */

GIDPOST_API
FCALLSCFUN4(INT, GiD_fResultDescriptionDim, GID_FRESULTDESCRIPTIONDIM, gid_fresultdescriptiondim,
            INT,     /* FILE handler */
	    STRING,  /* Result name */
	    INT,     /* result type */
            INT);    /* result dimension */

GIDPOST_API
FCALLSCFUN0(INT, GiD_ResultValues, GID_RESULTVALUES, gid_resultvalues);

GIDPOST_API
FCALLSCFUN1(INT, GiD_fResultValues, GID_FRESULTVALUES, gid_fresultvalues, INT);

GIDPOST_API
FCALLSCFUN1(INT, GiD_ResultRange, GID_RESULTRANGE, gid_resultrange, STRING);

GIDPOST_API
FCALLSCFUN2(INT, GiD_fResultRange, GID_FRESULTRANGE, gid_fresultrange, INT,STRING);


int GiD_ScalarComp(GP_CONST char * Comp1);

int GiD_fScalarComp(GiD_FILE fd, GP_CONST char * Comp1);

GIDPOST_API
FCALLSCFUN1(INT, GiD_ScalarComp, GID_SCALARCOMP, gid_scalarcomp, STRING);

GIDPOST_API
FCALLSCFUN2(INT, GiD_fScalarComp, GID_FSCALARCOMP, gid_fscalarcomp, INT,STRING);

int GiD_VectorComp(GP_CONST char * Comp1, GP_CONST char * Comp2,
		   GP_CONST char * Comp3, GP_CONST char * Comp4);

int GiD_fVectorComp(GiD_FILE fd,
                    GP_CONST char * Comp1, GP_CONST char * Comp2,
                    GP_CONST char * Comp3, GP_CONST char * Comp4);

GIDPOST_API
FCALLSCFUN4(INT, GiD_VectorComp, GID_VECTORCOMP, gid_vectorcomp,
	    STRING,  /* component 1 */
	    STRING,  /* component 2 */
	    STRING,  /* component 3 */
	    STRING); /* component 4 (modulus) */

GIDPOST_API
FCALLSCFUN5(INT, GiD_fVectorComp, GID_FVECTORCOMP, gid_fvectorcomp,
            INT,     /* FILE handler */
	    STRING,  /* component 1 */
	    STRING,  /* component 2 */
	    STRING,  /* component 3 */
	    STRING); /* component 4 (modulus) */

int GiD_2DMatrixComp(GP_CONST char * Comp1, GP_CONST char * Comp2, GP_CONST char * Comp3);

int GiD_f2DMatrixComp(GiD_FILE fd,
                      GP_CONST char * Comp1, GP_CONST char * Comp2, GP_CONST char * Comp3);

FCALLSCFUN3(INT, GiD_2DMatrixComp, GID_2DMATRIXCOMP, gid_2dmatrixcomp,
	    STRING,  /* component Sxx */
	    STRING,  /* component Syy */
	    STRING); /* component Sxy */

int GiD_3DMatrixComp(GP_CONST char * Comp1, GP_CONST char * Comp2, GP_CONST char * Comp3,
		     GP_CONST char * Comp4, GP_CONST char * Comp5, GP_CONST char * Comp6);

int GiD_f3DMatrixComp(GiD_FILE fd,
                      GP_CONST char * Comp1, GP_CONST char * Comp2, GP_CONST char * Comp3,
                      GP_CONST char * Comp4, GP_CONST char * Comp5, GP_CONST char * Comp6);

GIDPOST_API
FCALLSCFUN6(INT, GiD_3DMatrixComp, GID_3DMATRIXCOMP, gid_3dmatrixcomp,
	    STRING,  /* component Sxx */
	    STRING,  /* component Syy */
	    STRING,  /* component Szz */
	    STRING,  /* component Sxy */
	    STRING,  /* component Syz */
	    STRING); /* component Sxz */

GIDPOST_API
FCALLSCFUN7(INT, GiD_f3DMatrixComp, GID_F3DMATRIXCOMP, gid_f3dmatrixcomp,
            INT,     /* FILE handler */
	    STRING,  /* component Sxx */
	    STRING,  /* component Syy */
	    STRING,  /* component Szz */
	    STRING,  /* component Sxy */
	    STRING,  /* component Syz */
	    STRING); /* component Sxz */

int GiD_PDMComp(GP_CONST char *Comp1,
                GP_CONST char *Comp2,
                GP_CONST char *Comp3,
                GP_CONST char *Comp4);

int GiD_fPDMComp(GiD_FILE fd,
                 GP_CONST char *Comp1,
                 GP_CONST char *Comp2,
                 GP_CONST char *Comp3,
                 GP_CONST char *Comp4);

GIDPOST_API
FCALLSCFUN4(INT, GiD_PDMComp, GID_PDMCOMP, gid_pdmcomp,
	    STRING,  /* component Sxx */
	    STRING,  /* component Syy */
	    STRING,  /* component Sxy */
	    STRING); /* component Szz */

GIDPOST_API
FCALLSCFUN5(INT, GiD_fPDMComp, GID_FPDMCOMP, gid_fpdmcomp,
            INT,     /* FILE handler */
	    STRING,  /* component Sxx */
	    STRING,  /* component Syy */
	    STRING,  /* component Sxy */
	    STRING); /* component Szz */

int GiD_MainMatrixComp(GP_CONST char * Comp1,  GP_CONST char * Comp2,  GP_CONST char * Comp3,
		       GP_CONST char * Comp4,  GP_CONST char * Comp5,  GP_CONST char * Comp6,
		       GP_CONST char * Comp7,  GP_CONST char * Comp8,  GP_CONST char * Comp9,
		       GP_CONST char * Comp10, GP_CONST char * Comp11, GP_CONST char * Comp12);

int GiD_fMainMatrixComp(GiD_FILE fd,
                       GP_CONST char * Comp1,  GP_CONST char * Comp2,  GP_CONST char * Comp3,
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

GIDPOST_API
FCALLSCFUN13(INT, GiD_fMainMatrixComp, GID_FMAINMATRIXCOMP, gid_fmainmatrixcomp,
             INT,     /* FILE handler */
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

int GiD_fLAComponents(GiD_FILE fd,
                      GP_CONST char * Comp1, GP_CONST char * Comp2, GP_CONST char * Comp3);

GIDPOST_API
FCALLSCFUN3(INT, GiD_LAComponents, GID_LACOMPONENTS, gid_lacomponents,
	    STRING,  /* local axes name 1 */
	    STRING,  /* local axes name 2 */
	    STRING); /* local axes name 3 */

GIDPOST_API
FCALLSCFUN4(INT, GiD_fLAComponents, GID_FLACOMPONENTS, gid_flacomponents,
            INT,     /* FILE handler */
	    STRING,  /* local axes name 1 */
	    STRING,  /* local axes name 2 */
	    STRING); /* local axes name 3 */

GIDPOST_API
FCALLSCFUN0(INT,GiD_EndResult,GID_ENDRESULT,gid_endresult);

GIDPOST_API
FCALLSCFUN1(INT,GiD_fEndResult,GID_FENDRESULT,gid_fendresult,INT);

GIDPOST_API
FCALLSCFUN1(INT,GiD_ResultUnit,GID_RESULTUNIT,gid_resultunit,STRING);

GIDPOST_API
FCALLSCFUN2(INT,GiD_fResultUnit,GID_FRESULTUNIT,gid_fresultunit,INT,STRING);

GIDPOST_API
FCALLSCFUN2(INT, GiD_WriteScalar,GID_WRITESCALAR,gid_writescalar,INT,DOUBLE);

GIDPOST_API
FCALLSCFUN3(INT, GiD_fWriteScalar,GID_FWRITESCALAR,gid_fwritescalar,
            INT,INT,DOUBLE);

GIDPOST_API
FCALLSCFUN3(INT, GiD_Write2DVector,GID_WRITE2DVECTOR,gid_write2dvector,
            INT,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN4(INT, GiD_fWrite2DVector,GID_FWRITE2DVECTOR,gid_fwrite2dvector,
            INT,INT,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN4(INT, GiD_WriteVector,GID_WRITEVECTOR,gid_writevector,
            INT,DOUBLE,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN5(INT, GiD_fWriteVector,GID_FWRITEVECTOR,gid_fwritevector,
            INT,INT,DOUBLE,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN5(INT, GiD_WriteVectorModule,GID_WRITEVECTORMODULE,gid_writevectormodule,
            INT,DOUBLE,DOUBLE,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN6(INT, GiD_fWriteVectorModule,GID_FWRITEVECTORMODULE,gid_fwritevectormodule,
            INT,INT,DOUBLE,DOUBLE,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN4(INT, GiD_Write2DMatrix,GID_WRITE2DMATRIX,gid_write2dmatrix,
            INT,DOUBLE,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN5(INT, GiD_fWrite2DMatrix,GID_FWRITE2DMATRIX,gid_fwrite2dmatrix,
            INT,INT,DOUBLE,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN7(INT, GiD_Write3DMatrix,GID_WRITE3DMATRIX,gid_write3dmatrix,
            INT,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN8(INT, GiD_fWrite3DMatrix,GID_FWRITE3DMATRIX,gid_fwrite3dmatrix,
            INT,INT,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN5(INT, GiD_WritePlainDefMatrix,GID_WRITEPLAINDEFMATRIX,gid_writeplaindefmatrix,
            INT,DOUBLE,DOUBLE,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN6(INT, GiD_fWritePlainDefMatrix,GID_FWRITEPLAINDEFMATRIX,gid_fwriteplaindefmatrix,
            INT,INT,DOUBLE,DOUBLE,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN13(INT, GiD_WriteMainMatrix,GID_WRITEMAINMATRIX,gid_writemainmatrix,
             INT,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN14(INT, GiD_fWriteMainMatrix,GID_FWRITEMAINMATRIX,gid_fwritemainmatrix,
             INT,INT,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN4(INT, GiD_WriteLocalAxes,GID_WRITELOCALAXES,gid_writelocalaxes,
            INT,DOUBLE,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN5(INT, GiD_fWriteLocalAxes,GID_FWRITELOCALAXES,gid_fwritelocalaxes,
            INT,INT,DOUBLE,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN3(INT, GiD_WriteComplexScalar,GID_WRITECOMPLEXSCALAR,gid_writecomplexscalar,
            INT,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN4(INT, GiD_fWriteComplexScalar,GID_FWRITECOMPLEXSCALAR,gid_fwritecomplexscalar,
            INT,INT,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN5(INT, GiD_Write2DComplexVector,GID_WRITE2DCOMPLEXVECTOR,gid_write2dcomplexvector,
            INT,DOUBLE,DOUBLE,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN6(INT, GiD_fWrite2DComplexVector,GID_FWRITE2DCOMPLEXVECTOR,gid_fwrite2dcomplexvector,
            INT,INT,DOUBLE,DOUBLE,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN7(INT, GiD_WriteComplexVector,GID_WRITECOMPLEXVECTOR,gid_writecomplexvector,
            INT,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE);

GIDPOST_API
FCALLSCFUN8(INT, GiD_fWriteComplexVector,GID_FWRITECOMPLEXVECTOR,gid_fwritecomplexvector,
            INT,INT,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE);


// GiD_fWriteCoordinatesBlock includes BeginCoordinates() and EndCoordinates()
// GIDPOST_API int GiD_fWriteCoordinatesBlock( GiD_FILE fd, int num_points, GP_CONST double *xyz_array );
GIDPOST_API
FCALLSCFUN3(INT,GiD_fWriteCoordinatesBlock,GID_FWRITECOORDINATESBLOCK,gid_fwritecoordinatesblock,
            INT,INT,DOUBLEV);
// GIDPOST_API int GiD_fWriteCoordinatesIdBlock( GiD_FILE fd, int num_points, GP_CONST int *list_ids, GP_CONST double *xyz_array );
GIDPOST_API
FCALLSCFUN4(INT,GiD_fWriteCoordinatesIdBlock,GID_FWRITECOORDINATESIDBLOCK,gid_fwritecoordinatesidblock,
            INT,INT,INTV,DOUBLEV);


// // GiD_fWriteElementsBlock includes BeginElements() and EndElements()
// GIDPOST_API int GiD_fWriteElementsBlock( GiD_FILE fd, int num_elements, GP_CONST int *connectivities );
GIDPOST_API
FCALLSCFUN3(INT,GiD_fWriteElementsBlock,GID_FWRITEELEMENTSBLOCK,gid_fwriteelementsblock, INT, INT, INTV);
// GIDPOST_API int GiD_fWriteElementsIdBlock( GiD_FILE fd, int num_elements, GP_CONST int *list_ids, GP_CONST int *connectivities );
GIDPOST_API
FCALLSCFUN4(INT,GiD_fWriteElementsIdBlock,GID_FWRITEELEMENTSIDBLOCK,gid_fwriteelementsidblock, INT, INT, INTV, INTV);
// GIDPOST_API int GiD_fWriteElementsMatBlock( GiD_FILE fd, int num_elements, GP_CONST int *connectivities, GP_CONST int *lst_material_id );
GIDPOST_API
FCALLSCFUN4(INT,GiD_fWriteElementsMatBlock,GID_FWRITEELEMENTSMATBLOCK,gid_fwriteelementsmatblock, INT, INT, INTV, INTV);
// GIDPOST_API int GiD_fWriteElementsIdMatBlock( GiD_FILE fd, int num_elements, GP_CONST int *list_ids, GP_CONST int *connectivities,
//                                              GP_CONST int *lst_material_id );
GIDPOST_API
FCALLSCFUN5(INT,GiD_fWriteElementsIdMatBlock,GID_FWRITEELEMENTSIDMATBLOCK,gid_fwriteelementsidmatblock, INT, INT, INTV, INTV, INTV);


// // GiD_fWriteResultBlock includes BeginResult()/BeginValues() and EndValues()/EndResult()
// GIDPOST_API
// int GiD_fWriteResultBlock( GiD_FILE fd, GP_CONST char *result_name, GP_CONST char *analysis_name, double step_value,
//                            GiD_ResultType result_type, GiD_ResultLocation result_location,
//                            GP_CONST char *gauss_point_name, GP_CONST char *range_table_name, int num_component_names,
//                            GP_CONST char *list_component_names[], 
// 			   GP_CONST char *unit_name,
// 			   int num_result_values, GP_CONST int *list_result_ids,
// 			   // list_component_values === 
// 			   //   num_component_values == 1 --> scalar_array
// 			   //   num_component_values == 2 --> VxVyVz_array
// 			   //   num_component_values == 6 --> SxxSyySxySzzSxzSyz_array
// 			   //   ...
//                            int num_component_values, GP_CONST double *list_component_values );
//                            
GIDPOST_API
FCALLSCFUN15(INT,GiD_fWriteResultBlock,GID_FWRITERESULTBLOCK,gid_fwriteresultblock, 
                INT, STRING, STRING, DOUBLE,
                INT, INT, STRING, STRING,
                INT, STRING, // should be STRINGV instead of STRING, but it seems to be unknown/undefined
                STRING,
                INT, INTV, INT, DOUBLEV);

