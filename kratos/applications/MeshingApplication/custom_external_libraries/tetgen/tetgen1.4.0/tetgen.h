///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TetGen                                                                    //
//                                                                           //
// A Quality Tetrahedral Mesh Generator and 3D Delaunay Triangulator         //
//                                                                           //
// Version 1.4                                                               //
// January 14, 2006                                                          //
//                                                                           //
// Copyright 2002, 2004, 2005, 2006                                          //
// Hang Si                                                                   //
// Rathausstr. 9, 10178 Berlin, Germany                                      //
// si@wias-berlin.de                                                         //
//                                                                           //
// You can obtain TetGen via internet: http://tetgen.berlios.de.  It may be  //
//   freely copied, modified, and redistributed under the copyright notices  //
//   given in the file LICENSE.                                              //
//                                                                           //
// TetGen computes Delaunay tetrahedralizations, constrained Delaunay tetra- //
//   hedralizations, and quality Delaunay tetrahedral meshes. The latter are //
//   nicely graded and whose tetrahedra have radius-edge ratio bounded. Such //
//   meshes are suitable for finite element and finite volume methods.       //
//                                                                           //
// TetGen incorporates a suit of geometrical and mesh generation algorithms. //
//   A brief description of algorithms used in TetGen is found in the first  //
//   section of the user's manual.  References are given for users who are   //
//   interesting in these approaches. The main references are given below:   //
//                                                                           //
//   The efficient Delaunay tetrahedralization algorithm is: H. Edelsbrunner //
//   and N. R. Shah, "Incremental Topological Flipping Works for Regular     //
//   Triangulations". Algorithmica 15: 223-241, 1996.                        //
//                                                                           //
//   The constrained Delaunay tetrahedralization algorithm is described in:  //
//   H. Si and K. Gaertner, "Meshing Piecewise Linear Complexs by Constrain- //
//   ed Delaunay Tetrahedralizations". Proceedings of the 14th International //
//   Mesh Roundtable. September 2005.                                        //
//                                                                           //
//   The Delaunay refinement algorithm is from: J. R. Shewchuk, "Tetrahedral //
//   Mesh Generation by Delaunay Refinement". Proceedings of the 14th Annual //
//   Symposium on Computational Geometry, pages 86-95, 1998.                 //
//                                                                           //
// The mesh data structure of TetGen is a combination of two types of mesh   //
//   data structures.  The tetrahedron-based mesh data structure introduced  //
//   by Shewchuk is eligible to implement tetrahedralization algorithms. The //
//   triangle-edge data structure provided by Muecke is used for represent-  //
//   ing boundary elements: subfaces and subsegments.  These data types are  //
//   handled through a set of fast mesh manipulation primitives.             //
//                                                                           //
//   J. R. Shewchuk, "Delaunay Refinement Mesh Generation". PhD thesis,      //
//   Carnegie Mellon University, 1997.                                       //
//                                                                           //
//   E. P. Muecke, "Shapes and Implementations in Three-Dimensional          //
//   Geometry". PhD thesis, Univ. of Illinois, Urbana, Illinois, 1993.       //
//                                                                           //
// The research of mesh generation is definitly on the move. Many State-of-  //
//   the-art algorithms need implementing and evaluating. I heartily welcome //
//   any new algorithm especially for generating quality conforming Delaunay //
//   meshes and anisotropic conforming Delaunay meshes.                      //
//                                                                           //
// TetGen is supported by the "pdelib" project of Weierstrass Institute for  //
//   Applied Analysis and Stochastics (WIAS) in Berlin.  It is a collection  //
//   of software components for solving non-linear partial differential      //
//   equations including 2D and 3D mesh generators, sparse matrix solvers,   //
//   and scientific visualization tools, etc.  For more information please   //
//   see: http://www.wias-berlin.de/software/pdelib.                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// tetgen.h                                                                  //
//                                                                           //
// Header file of the TetGen library. Also is the user-level header file.    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TetGen Library Overview                                                   //
//                                                                           //
// TetGen library is comprised by several data types and global functions.   //
//                                                                           //
// There are three main data types: tetgenio, tetgenbehavior, and tetgenmesh.//
// Tetgenio is used to pass data into and out of TetGen library; tetgenbeha- //
// vior keeps the runtime options and thus controls the behaviors of TetGen; //
// tetgenmesh, the biggest data type I've ever defined, contains mesh data   //
// structures and mesh traversing and transformation operators.  These data  //
// types are defined as C++ classes.                                         //
//                                                                           //
// There are few global functions. tetrahedralize() is provided for calling  //
// TetGen from another program. Two functions: orient3d() and insphere() are //
// incorporated from a public C code provided by Shewchuk.  They performing  //
// exact geometrical tests.                                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef tetgenH
#define tetgenH

// To compile TetGen as a library instead of an executable program, define
//   the TETLIBRARY symbol.

// #define TETLIBRARY

// Uncomment the following line to disable assert macros. These macros are
//   inserted in places where I hope to catch bugs.

// #define NDEBUG

// To insert lots of self-checks for internal errors, define the SELF_CHECK
//   symbol.  This will slow down the program significantly. 

// #define SELF_CHECK

// For single precision ( which will save some memory and reduce paging ),
//   define the symbol SINGLE by using the -DSINGLE compiler switch or by
//   writing "#define SINGLE" below.
//
// For double precision ( which will allow you to refine meshes to a smaller
//   edge length), leave SINGLE undefined.

// #define SINGLE

#ifdef SINGLE
  #define REAL float
#else
  #define REAL double
#endif 	// not defined SINGLE

// Here are the most general used head files for C/C++ programs.

#include <stdio.h>            // Standard IO: FILE, NULL, EOF, printf(), ...
#include <stdlib.h>        // Standard lib: abort(), system(), getenv(), ...
#include <string.h>         // String lib: strcpy(), strcat(), strcmp(), ...
#include <math.h>                     // Math lib: sin(), sqrt(), pow(), ...
#include <assert.h>
 
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The tetgenio data type                                                    //
//                                                                           //
// Used to pass data into and out of the library of TetGen.                  //
//                                                                           //
// If you want to program with the library of TetGen, it's necessary for you //
// to understand the tetgenio data type, while the other two data types can  //
// be hidden through calling the global function "tetrahedralize()".  As you //
// will see, that tetgenio is just a collection of arrays to storing points, //
// tetrahedra, (triangular) faces, boundary markers, and so forth.  They are //
// used to store data from input & output files of TetGen. If you understand //
// TetGen's file formats, then it is easy to understand these arrays.        //
//                                                                           //
// Once an object of tetgenio is declared (or created), all arrays of it are //
// automatically initialized to NULLs (by routine initialize()). Before they //
// can be used,  one has to first allocate enough memory for them, i.e., use //
// either 'malloc()' or 'new' operator. On deletion of the object, one needs //
// to free the memory occupied by these arrays.  Routine deinitialize() will //
// be automatically called. It will deallocate the memory for an array if it //
// is not a NULL. However, it assumes that the memory is allocated by 'new'  //
// (C++ operator).  If you use malloc(), you should free() them and set the  //
// pointers to NULLs before reaching deinitialize().                         //
//                                                                           //
// In all cases, the first item in any array is stored starting at index [0].//
// However, that item is item number `firstnumber' which may be '0' or '1'.  //
// Be sure to set the 'firstnumber' be '1' if your indices pointing into the //
// pointlist is starting from '1'. Default, it is initialized be '0'.        //
//                                                                           //
// Tetgenio also contains routines for reading and writing TetGen's files as //
// well.  Both the library of TetGen and TetView use these routines to parse //
// input files, i.e., .node, .poly, .smesh, .ele, .face, and .edge files.    //
// Other routines are provided mainly for debugging purpose.                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

class tetgenio {

  public:

    // Maximum number of characters in a file name (including the null).
    enum {FILENAMESIZE = 1024};

    // Maxi. numbers of chars in a line read from a file (incl. the null).
    enum {INPUTLINESIZE = 1024};

    // The polygon data structure.  A "polygon" is a planar polygon. It can
    //   be arbitrary shaped (convex or non-convex) and bounded by non-
    //   crossing segments, i.e., the number of vertices it has indictes the
    //   same number of edges.
    // 'vertexlist' is a list of vertex indices (integers), its length is
    //   indicated by 'numberofvertices'.  The vertex indices are odered in
    //   either counterclockwise or clockwise way.
    typedef struct {
      int *vertexlist;
      int numberofvertices;
    } polygon;

    static void init(polygon* p) {
      p->vertexlist = (int *) NULL;
      p->numberofvertices = 0;
    }

    // The facet data structure.  A "facet" is a planar facet.  It is used
    //   to represent a planar straight line graph (PSLG) in two dimension.
    //   A PSLG contains a list of polygons. It also may conatin holes in it,
    //   indicated by a list of hole points (their coordinates).
    typedef struct {
      polygon *polygonlist;
      int numberofpolygons;
      REAL *holelist;
      int numberofholes;
    } facet;

    static void init(facet* f) {
      f->polygonlist = (polygon *) NULL;
      f->numberofpolygons = 0;
      f->holelist = (REAL *) NULL;
      f->numberofholes = 0;
    }

    // The periodic boundary condition group data structure.  A "pbcgroup"
    //   contains the definition of a pbc and the list of pbc point pairs.
    //   'fmark1' and 'fmark2' are the facetmarkers of the two pbc facets f1
    //   and f2, respectively. 'transmat' is the transformation matrix which
    //   maps a point in f1 into f2.  An array of pbc point pairs are saved
    //   in 'pointpairlist'. The first point pair is at indices [0] and [1],
    //   followed by remaining pairs. Two integers per pair.
    typedef struct {
      int fmark1, fmark2;
      REAL transmat[4][4];
      int numberofpointpairs;
      int *pointpairlist;
    } pbcgroup;

  public:

    // Items are numbered starting from 'firstnumber' (0 or 1), default is 0.
    int firstnumber; 

    // Dimension of the mesh (2 or 3), default is 3.
    int mesh_dim;

    // `pointlist':  An array of point coordinates.  The first point's x
    //   coordinate is at index [0] and its y coordinate at index [1], its
    //   z coordinate is at index [2], followed by the coordinates of the
    //   remaining points.  Each point occupies three REALs. 
    // `pointattributelist':  An array of point attributes.  Each point's
    //   attributes occupy `numberofpointattributes' REALs. 
    // 'addpointlist':  An array of additional point coordinates.
    // 'addpointattributelist':  An array of attributes for addition points.
    // `pointmarkerlist':  An array of point markers; one int per point.
    REAL *pointlist;
    REAL *pointattributelist;
    REAL *addpointlist;
    REAL *addpointattributelist;
    int *pointmarkerlist;
    int numberofpoints;
    int numberofpointattributes;
    int numberofaddpoints;
 
    // `elementlist':  An array of element (triangle or tetrahedron) corners. 
    //   The first element's first corner is at index [0], followed by its
    //   other corners in counterclockwise order, followed by any other
    //   nodes if the element represents a nonlinear element.  Each element
    //   occupies `numberofcorners' ints.
    // `elementattributelist':  An array of element attributes.  Each
    //   element's attributes occupy `numberofelementattributes' REALs.
    // `elementconstraintlist':  An array of constraints, i.e. triangle's
    //   area or tetrahedron's volume; one REAL per element.  Input only.
    // `neighborlist':  An array of element neighbors; 3 or 4 ints per
    //   element.  Output only.
    int *tetrahedronlist;
    REAL *tetrahedronattributelist;
    REAL *tetrahedronvolumelist;
    int *neighborlist;
    int numberoftetrahedra;
    int numberofcorners;
    int numberoftetrahedronattributes;

    // `facetlist':  An array of facets.  Each entry is a structure of facet.
    // `facetmarkerlist':  An array of facet markers; one int per facet.
    facet *facetlist;
    int *facetmarkerlist;
    int numberoffacets;

    // `holelist':  An array of holes.  The first hole's x, y and z
    //   coordinates  are at indices [0], [1] and [2], followed by the
    //   remaining holes. Three REALs per hole. 
    REAL *holelist;
    int numberofholes;

    // `regionlist': An array of regional attributes and volume constraints.
    //   The first constraint's x, y and z coordinates are at indices [0],
    //   [1] and [2], followed by the regional attribute at index [3], foll-
    //   owed by the maximum volume at index [4]. Five REALs per constraint. 
    // Note that each regional attribute is used only if you select the `A'
    //   switch, and each volume constraint is used only if you select the
    //   `a' switch (with no number following).
    REAL *regionlist;
    int numberofregions;

    // `facetconstraintlist': An array of facet maximum area constraints.
    //   Two REALs per constraint. The first one is the facet marker (cast
    //   it to int), the second is its maximum area bound.
    // Note the 'facetconstraintlist' is used only for the 'q' switch. 
    REAL *facetconstraintlist;
    int numberoffacetconstraints;

    // `segmentconstraintlist': An array of seg maximum length constraints.
    //   Three REALs per constraint. The first two are the indices (pointing
    //   into 'pointlist') of the endpoints of the segment, the third is its
    //   maximum length bound.
    // Note the 'segmentconstraintlist' is used only for the 'q' switch. 
    REAL *segmentconstraintlist;
    int numberofsegmentconstraints;

    // `nodeconstraintlist':  An array of segment length constraints.  Two
    //   REALs per constraint. The first one is the index (pointing into
    //   'pointlist') of the node, the second is its edge length bound.
    // Note the 'nodeconstraintlist' is used only for the 'q' switch. 
    REAL *nodeconstraintlist;
    int numberofnodeconstraints;

    // 'pbcgrouplist':  An array of periodic boundary condition groups.
    pbcgroup *pbcgrouplist;
    int numberofpbcgroups;

    // `trifacelist':  An array of triangular face endpoints.  The first
    //   face's endpoints are at indices [0], [1] and [2], followed by the
    //   remaining faces.  Three ints per face.
    // `trifacemarkerlist':  An array of face markers; one int per face.
    int *trifacelist;
    int *trifacemarkerlist;
    int numberoftrifaces;

    // `edgelist':  An array of edge endpoints.  The first edge's endpoints
    //   are at indices [0] and [1], followed by the remaining edges.  Two
    //   ints per edge.
    // `edgemarkerlist':  An array of edge markers; one int per edge.
    int *edgelist;
    int *edgemarkerlist;
    int numberofedges;

  public:

    // Initialize routine.
    void initialize();
    void deinitialize();

    // Input & output routines.
    bool load_node_call(FILE* infile, int markers, char* nodefilename);
    bool load_node(char* filename);
    bool load_addnodes(char* filename);
    bool load_pbc(char* filename);
    bool load_var(char* filename);
    bool load_poly(char* filename);
    bool load_off(char* filename);
    bool load_ply(char* filename);
    bool load_stl(char* filename);
    bool load_medit(char* filename);
    bool load_plc(char* filename, int object);
    bool load_tetmesh(char* filename);
    void save_nodes(char* filename);
    void save_elements(char* filename);
    void save_faces(char* filename);
    void save_edges(char* filename);
    void save_neighbors(char* filename);
    void save_poly(char* filename);

    // Read line and parse string functions.
    char *readline(char* string, FILE* infile, int *linenumber);
    char *findnextfield(char* string);
    char *readnumberline(char* string, FILE* infile, char* infilename);
    char *findnextnumber(char* string);

    // Constructor and destructor.
    tetgenio() {initialize();}
    ~tetgenio() {deinitialize();}
};

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The tetgenbehavior data type                                              //
//                                                                           //
// Used to parse command line switches and file names.                       //
//                                                                           //
// It includes a list of variables corresponding to the commandline switches //
// for control the behavior of TetGen.  These varibales are all initialized  //
// to their default values.                                                  //
//                                                                           //
// parse_commandline() provides an simple interface to set the vaules of the //
// variables.  It accepts the standard parameters (e.g., 'argc' and 'argv')  //
// that pass to C/C++ main() function. Alternatively a string which contains //
// the command line options can be used as its parameter.                    //
//                                                                           //
// You don't need to understand this data type. It can be implicitly called  //
// by the global function "tetrahedralize()" defined below.  The necessary   //
// thing you need to know is the meaning of command line switches of TetGen. //
// They are described in the third section of the user's manual.             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

class tetgenbehavior {

  public:

    // Labels define the objects which are acceptable by TetGen. They are 
    //   recognized by the file extensions.
    //   - NODES, a list of nodes (.node); 
    //   - POLY, a piecewise linear complex (.poly or .smesh); 
    //   - OFF, a polyhedron (.off, Geomview's file format); 
    //   - PLY, a polyhedron (.ply, file format from gatech);
    //   - STL, a surface mesh (.stl, stereolithography format);
    //   - MEDIT, a surface mesh (.mesh, Medit's file format); 
    //   - MESH, a tetrahedral mesh (.ele).
    //   If no extension is available, the imposed commandline switch
    //   (-p or -r) implies the object. 

    enum objecttype {NONE, NODES, POLY, OFF, PLY, STL, MEDIT, MESH};

    // Variables of command line switches.  Each variable is corresponding
    //   to a specific switch and will be properly initialized.  Read the
    //   user's manul to find out the meaning of these switches.

    int plc;                                              // '-p' switch, 0.
    int refine;                                           // '-r' switch, 0.
    int quality;                                          // '-q' switch, 0.
    int varvolume;                         // '-a' switch without number, 0.
    int fixedvolume;                          // '-a' switch with number, 0.
    int insertaddpoints;                                  // '-i' switch, 0.
    int regionattrib;                                     // '-A' switch, 0.
    int detectinter;                                      // '-d' switch, 0.
    int offcenter;                                        // '-R' switch, 0.
    int zeroindex;                                        // '-z' switch, 0.
    int order;             // element order, specified after '-o' switch, 1.
    int facesout;                                         // '-f' switch, 0.
    int edgesout;                                         // '-e' switch, 0.
    int neighout;                                         // '-n' switch, 0.
    int meditview;                                        // '-g' switch, 0.
    int gidview;                                          // '-G' switch, 0.
    int geomview;                                         // '-O' switch, 0.
    int nobound;                                          // '-B' switch, 0.
    int nonodewritten;                                    // '-N' switch, 0.
    int noelewritten;                                     // '-E' switch, 0.
    int nofacewritten;                                    // '-F' switch, 0.
    int noiterationnum;                                   // '-I' switch, 0.
    int nomerge;           // count of how often '-M' switch is selected, 0.
    int nobisect;          // count of how often '-Y' switch is selected, 0.
    int noflip;                     // do not perform flips. '-Y' switch. 0.
    int nofreevert;                 // Remove free vertices. '-X' switch. 0.
    int nojettison;     // do not jettison redundants nodes. '-J' switch. 0.
    int steiner;                             // number after '-S' switch. 0.
    int dofullperturb;               // do full permutation. '-P' switch, 0.
    int docheck;                                          // '-C' switch, 0.
    int quiet;                                            // '-Q' switch, 0.
    int verbose;           // count of how often '-V' switch is selected, 0.
    int useshelles;            // '-p', '-r', '-q', '-d', or '-c' switch, 0.
    REAL minratio;                         // number after '-q' switch, 2.0.
    REAL goodratio;               // number calculated from 'minratio', 0.0. 
    REAL minangle;                             // minimum angle bound, 20.0.
    REAL goodangle;                      // cosine squared of minangle, 0.0.
    REAL maxvolume;                       // number after '-a' switch, -1.0.
    REAL maxdihedral;                    // number after '-s' switch, 175.0.
    REAL alpha3;                           // number after '-R' switch, 0.6.
    REAL epsilon;                       // number after '-T' switch, 1.0e-8.
    enum objecttype object;         // determined by -p, or -r switch. NONE.

    // Variables used to save command line switches and in/out file names.
    char commandline[1024];
    char infilename[1024];
    char outfilename[1024];

    tetgenbehavior();
    ~tetgenbehavior() {}

    void versioninfo();
    void syntax();
    void usage();

    // Command line parse routine.
    bool parse_commandline(int argc, char **argv);
    bool parse_commandline(char *switches) {
      return parse_commandline(0, &switches);
    }
};

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Geometric predicates                                                      //
//                                                                           //
// Return one of the values +1, 0, and -1 on basic geometric questions such  //
// as the orientation of point sets, in-circle, and in-sphere tests.  They   //
// are basic units for the composition of geometric algorithms.  TetGen uses //
// two 3D geometric predicates, which are the orientation test and the in-   //
// sphere test (e.g. the locally Deklaunay test).                            //
//                                                                           //
// Orientation test:  let a, b, c be a sequence of 3 non-collinear points in //
// R^3.  They defines a unique hypeplane H.  Let H+ and H- be the two spaces //
// separated by H, which are defined as follows (using the left-hand rule):  //
// make a fist using your left hand in such a way that your fingers follow   //
// the order of a, b and c, then your thumb is pointing to H+.  Given any    //
// point d in R^3, the orientation test returns +1 if d lies in H+, -1 if d  //
// lies in H-, or 0 if d lies on H.                                          //
//                                                                           //
// In-sphere test:  let a, b, c, d be 4 non-coplanar points in R^3.  They    //
// defines a unique circumsphere S.  Given any point e in R^3, the in-sphere //
// test returns +1 if e lies inside S, or -1 if e lies outside S, or 0 if e  //
// lies on S.                                                                //
//                                                                           //
// The correctness of geometric predicates is crucial for the control flow   //
// and hence for the correctness and robustness of an implementation of a    //
// geometric algorithm.  The following routines use arbitrary precision      //
// floating-point arithmetic. They are fast and robust. It is provided by J. //
// Schewchuk in public domain (http://www.cs.cmu.edu/~quake/robust.html).    //
// The source code are found in a separate file "predicates.cxx".            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

REAL exactinit();
REAL orient3d(REAL *pa, REAL *pb, REAL *pc, REAL *pd);
REAL insphere(REAL *pa, REAL *pb, REAL *pc, REAL *pd, REAL *pe);

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The tetgenmesh data type                                                  //
//                                                                           //
// Includes data types and mesh routines for creating tetrahedral meshes and //
// Delaunay tetrahedralizations, mesh input & output, and so on.             //
//                                                                           //
// An object of tetgenmesh can be used to store a triangular or tetrahedral  //
// mesh and its settings. TetGen's functions operates on one mesh each time. //
// This type allows reusing of the same function for different meshes.       //
//                                                                           //
// The mesh data structure (tetrahedron-based and triangle-edge data struct- //
// ures) are declared. There are other accessary data type defined as well,  //
// for efficient memory management and link list operations, etc.            //
//                                                                           //
// All algorithms TetGen used are implemented in this data type as member    //
// functions. References of these algorithms can be found in user's manual.  //
//                                                                           //
// It's not necessary to understand this type. There is a global function    //
// "tetrahedralize()" (defined at the end of this file) implicitly creates   //
// the object and calls its member functions according to the command line   //
// switches you specified.                                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

class tetgenmesh {

  public:

    // Maximum number of characters in a file name (including the null).
    enum {FILENAMESIZE = 1024};

    // For efficiency, a variety of data structures are allocated in bulk.
    //   The following constants determine how many of each structure is
    //   allocated at once.
    enum {VERPERBLOCK = 4092, SUBPERBLOCK = 4092, ELEPERBLOCK = 8188};

    // Used for the point location scheme of Mucke, Saias, and Zhu, to
    //   decide how large a random sample of tetrahedra to inspect.
    enum {SAMPLEFACTOR = 11};

    // Labels that signify two edge rings of a triangle defined in Muecke's
    //   triangle-edge data structure, one (CCW) traversing edges in count-
    //   erclockwise direction and one (CW) in clockwise direction.
    enum {CCW = 0, CW = 1};

    // Labels that signify whether a record consists primarily of pointers
    //   or of floating-point words.  Used to make decisions about data
    //   alignment.
    enum wordtype {POINTER, FLOATINGPOINT};

    // Labels that signify the type of a vertex. An UNUSEDVERTEX is a vertex
    //   read from input (.node file or tetgenio structure) or an isolated
    //   vertex (outside the mesh).  It is the default type for a newpoint.
    enum verttype {UNUSEDVERTEX, DUPLICATEDVERTEX, NACUTEVERTEX, ACUTEVERTEX,
                   FREESEGVERTEX, FACETVERTEX, FREESUBVERTEX, FREEVOLVERTEX,
                   DEADVERTEX = -32768};
 
    // Labels that signify the type of a subface.  A subface is SHARPSUB if
    //   it is in a facet which forms an acute dihedral angle (less than 90
    //   degree with other facets). It is SKINNYSUB if it has two segments
    //   form a small angle (less than 10 degree).
    enum shestype {NSHARPNSKINNY, SHARPSUB, SKINNYSUB, SHARPSKINNYSUB};

    // Labels that signify the type of flips can be applied on a face.
    //   A flipable face has the one of the types T23, T32, T22, and T44.
    //   Types UNFLIPABLE, NONCONVEX are unflipable.
    enum fliptype {T23, T32, T22, T44, UNFLIPABLE, FORBIDDENFACE,
                   FORBIDDENEDGE, NONCONVEX};

    // Labels that signify the result of triangle-triangle intersection test.
    //   Two triangles are DISJOINT, or adjoint at a vertex SHAREVERTEX, or
    //   adjoint at an edge SHAREEDGE, or coincident SHAREFACE or INTERSECT.
    enum interresult {DISJOINT, SHAREVERTEX, SHAREEDGE, SHAREFACE, INTERSECT};

    // Labels that signify the result of point location.  The result of a
    //   search indicates that the point falls inside a tetrahedron, inside
    //   a triangle, on an edge, on a vertex, or outside the mesh. 
    enum locateresult {INTETRAHEDRON, ONFACE, ONEDGE, ONVERTEX, OUTSIDE};

    // Labels that signify the result of vertex insertion.  The result
    //   indicates that the vertex was inserted with complete success, was
    //   inserted but encroaches upon a subsegment, was not inserted because
    //   it lies on a segment, or was not inserted because another vertex
    //   occupies the same location.
    enum insertsiteresult {SUCCESSINTET, SUCCESSONFACE, SUCCESSONEDGE,
                           DUPLICATEPOINT, OUTSIDEPOINT};

    // Labels that signify the result of direction finding.  The result
    //   indicates that a segment connecting the two query points accross
    //   an edge of the direction triangle/tetrahedron, across a face of
    //   the direction tetrahedron, along the left edge of the direction
    //   triangle/tetrahedron, along the right edge of the direction
    //   triangle/tetrahedron, or along the top edge of the tetrahedron.
    enum finddirectionresult {ACROSSEDGE, ACROSSFACE, LEFTCOLLINEAR,
                              RIGHTCOLLINEAR, TOPCOLLINEAR, BELOWHULL};

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The basic mesh element data structures                                    //
//                                                                           //
// There are four types of mesh elements: tetrahedra, subfaces, subsegments, //
// and points,  where subfaces and subsegments are triangles and edges which //
// appear on boundaries.  A tetrahedralization of a 3D point set comprises   //
// tetrahedra and points;  a surface mesh of a 3D domain comprises subfaces  //
// subsegments and points.  The elements of all the four types consist of a  //
// tetrahedral mesh of a 3D domain.  However, TetGen uses three data types:  //
// 'tetrahedron', 'shellface', and 'point'. A 'tetrahedron' is a tetrahedron;//
// while a 'shellface' can be either a subface or a subsegment; and a 'point'//
// is a point.  These three data types, linked by pointers comprise a mesh.  //
//                                                                           //
// A tetrahedron primarily consists of a list of 4 pointers to its corners,  //
// a list of 4 pointers to its adjoining tetrahedra, a list of 4 pointers to //
// its adjoining subfaces (when subfaces are needed). Optinoally, (depending //
// on the selected switches), it may contain an arbitrary number of user-    //
// defined floating-point attributes,  an optional maximum volume constraint //
// (for -a switch), and a pointer to a list of high-order nodes (-o2 switch).//
// Since the size of a tetrahedron is not determined until running time, it  //
// is not simply declared as a structure.                                    //
//                                                                           //
// The data structure of tetrahedron also stores the geometrical information.//
// Let t be a tetrahedron, v0, v1, v2, and v3 be the 4 nodes corresponding   //
// to the order of their storage in t.  v3 always has a negative orientation //
// with respect to v0, v1, v2 (ie,, v3 lies above the oriented plane passes  //
// through v0, v1, v2). Let the 4 faces of t be f0, f1, f2, and f3. Vertices //
// of each face are stipulated as follows: f0 (v0, v1, v2), f1 (v0, v3, v1), //
// f2 (v1, v3, v2), f3 (v2, v3, v0).                                         //
//                                                                           //
// A subface has 3 pointers to vertices, 3 pointers to adjoining subfaces, 3 //
// pointers to adjoining subsegments, 2 pointers to adjoining tetrahedra, a  //
// boundary marker(an integer). Like a tetrahedron, the pointers to vertices,//
// subfaces, and subsegments are ordered in a way that indicates their geom- //
// etric relation.  Let s be a subface, v0, v1 and v2 be the 3 nodes corres- //
// ponding to the order of their storage in s,  e0, e1 and e2 be the 3 edges,//
// then we have: e0 (v0, v1), e1 (v1, v2), e2 (v2, v0).                      //
//                                                                           //
// A subsegment has exactly the same data fields as a subface has, but only  //
// uses some of them. It has 2 pointers to its endpoints, 2 pointers to its  //
// adjoining (and collinear) subsegments, a pointer to a subface containing  //
// it (there may exist any number of subfaces having it, choose one of them  //
// arbitrarily). The geometric relation between its endpoints and adjoining  //
// subsegments is kept with respect to the storing order of its endpoints.   //
//                                                                           //
// The data structure of point is relatively simple.  A point is a list of   //
// floating-point numbers, starting with the x, y, and z coords, followed by //
// an arbitrary number of optional user-defined floating-point attributes,   //
// an integer boundary marker, an integer for the point type, and a pointer  //
// to a tetrahedron (used for speeding up point location).                   //
//                                                                           //
// For a tetrahedron on a boundary (or a hull) of the mesh, some or all of   //
// the adjoining tetrahedra may not be present. For an interior tetrahedron, //
// often no neighboring subfaces are present,  Such absent tetrahedra and    //
// subfaces are never represented by the NULL pointers; they are represented //
// by two special records: `dummytet', the tetrahedron fills "outer space",  //
// and `dummysh',  the vacuous subfaces which are omnipresent.               //
//                                                                           //
// Tetrahedra and adjoining subfaces are glued together through the pointers //
// saved in each data fields of them. Subfaces and adjoining subsegments are //
// connected in the same fashion.  However, there are no pointers directly   //
// gluing tetrahedra and adjoining subsegments.  For the purpose of saving   //
// space, the connections between tetrahedra and subsegments are entirely    //
// mediated through subfaces.  The following part explains how subfaces are  //
// connected in TetGen.                                                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The subface-subface and subface-subsegment connections                    //
//                                                                           //
// Adjoining subfaces sharing a common edge are connected in such a way that //
// they form a face ring around the edge. It is indeed a single linked list  //
// which is cyclic, e.g., one can start from any subface in it and traverse  //
// back. When the edge is not a subsegment, the ring only has two coplanar   //
// subfaces which are pointing to each other. Otherwise, the face ring may   //
// have any number of subfaces (and are not all coplanar).                   //
//                                                                           //
// How is the face ring formed?  Let s be a subsegment, f is one of subfaces //
// containing s as an edge.  The direction of s is stipulated from its first //
// endpoint to its second (according to their storage in s). Once the dir of //
// s is determined, the other two edges of f are oriented to follow this dir.//
// The "directional normal" N_f is a vector formed from any point in f and a //
// points orthogonally above f.                                              //
//                                                                           //
// The face ring of s is a cyclic ordered set of subfaces containing s, i.e.,//
// F(s) = {f1, f2, ..., fn}, n >= 1.  Where the order is defined as follows: //
// let fi, fj be two faces in F(s), the "normal-angle", NAngle(i,j) (range   //
// from 0 to 360 degree) is the angle between the N_fi and N_fj;  then fi is //
// in front of fj (or symbolically, fi < fj) if there exists another fk in   //
// F(s), and NAangle(k, i) < NAngle(k, j).  The face ring of s is: f1 < f2 < //
// ... < fn < f1.                                                            //
//                                                                           //
// The easiest way to imagine how a face ring is formed is to use the right- //
// hand rule.  Make a fist using your right hand with the thumb pointing to  //
// the direction of the subsegment. The face ring is connected following the //
// direction of your fingers.                                                //
//                                                                           //
// The subface and subsegment are also connected through pointers stored in  //
// their own data fields.  Every subface has a pointer to its adjoining sub- //
// segment. However, a subsegment only has one pointer to a subface which is //
// containing it. Such subface can be chosen arbitrarily, other subfaces are //
// found through the face ring.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

    // The tetrahedron data structure.  Fields of a tetrahedron contains:
    //   - a list of four adjoining tetrahedra;
    //   - a list of four vertices;
    //   - a list of four subfaces (optional, used for -p switch);
    //   - a list of user-defined floating-point attributes (optional);
    //   - a volume constraint (optional, used for -a switch);
    //   - a pointer to a list of high-ordered nodes (optional, -o2 switch);

    typedef REAL **tetrahedron;

    // The shellface data structure.  Fields of a shellface contains:
    //   - a list of three adjoining subfaces;
    //   - a list of three vertices;
    //   - a list of two adjoining tetrahedra;
    //   - a list of three adjoining subsegments;
    //   - a pointer to a badface containing it (used for -q);
    //   - an area constraint (optional, used for -q);
    //   - an integer for boundary marker;
    //   - an integer for type: SHARPSEGMENT, NONSHARPSEGMENT, ...;
    //   - an integer for pbc group (optional, if in->pbcgrouplist exists);

    typedef REAL **shellface;

    // The point data structure.  It is actually an array of REALs:
    //   - x, y and z coordinates;
    //   - a list of user-defined point attributes (optional);
    //   - an edge length constraint (optional, used for -q);
    //   - a pointer to a simplex (tet, tri, edge, or vertex);
    //   - a pointer to a parent (or duplicate) point;
    //   - a pointer to another pbc point (optional);
    //   - an integer for boundary marker;
    //   - an integer for verttype: INPUTVERTEX, FREEVERTEX, ...;

    typedef REAL *point;

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The mesh handle (triface, face) data types                                //
//                                                                           //
// Two special data types, 'triface' and 'face' are defined for maintaining  //
// and updating meshes. They are like pointers (or handles), which allow you //
// to hold one particular part of the mesh, i.e., a tetrahedron, a triangle, //
// an edge and a vertex.  However, these data types do not themselves store  //
// any part of the mesh. The mesh is made of the data types defined above.   //
//                                                                           //
// Muecke's "triangle-edge" data structure is the prototype for these data   //
// types.  It allows a universal representation for every tetrahedron,       //
// triangle, edge and vertex.  For understanding the following descriptions  //
// of these handle data structures,  readers are required to read both the   //
// introduction and implementation detail of "triangle-edge" data structure  //
// in Muecke's thesis.                                                       //
//                                                                           //
// A 'triface' represents a face of a tetrahedron and an oriented edge of    //
// the face simultaneously.  It has a pointer 'tet' to a tetrahedron, an     //
// integer 'loc' (range from 0 to 3) as the face index, and an integer 'ver' //
// (range from 0 to 5) as the edge version. A face of the tetrahedron can be //
// uniquly determined by the pair (tet, loc), and an oriented edge of this   //
// face can be uniquly determined by the triple (tet, loc, ver).  Therefore, //
// different usages of one triface are possible.  If we only use the pair    //
// (tet, loc), it refers to a face, and if we add the 'ver' additionally to  //
// the pair, it is an oriented edge of this face.                            //
//                                                                           //
// A 'face' represents a subface and an oriented edge of it simultaneously.  //
// It has a pointer 'sh' to a subface, an integer 'shver'(range from 0 to 5) //
// as the edge version.  The pair (sh, shver) determines a unique oriented   //
// edge of this subface.  A 'face' is also used to represent a subsegment,   //
// in this case, 'sh' points to the subsegment, and 'shver' indicates the    //
// one of two orientations of this subsegment, hence, it only can be 0 or 1. //
//                                                                           //
// Mesh navigation and updating are accomplished through a set of mesh       //
// manipulation primitives which operate on trifaces and faces.  They are    //
// introduced below.                                                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

    class triface {

      public:

        tetrahedron* tet;
        int loc, ver;

        // Constructors;
        triface() : tet(0), loc(0), ver(0) {}
        // Operators;
        triface& operator=(const triface& t) {
          tet = t.tet; loc = t.loc; ver = t.ver;
          return *this;
        }
        bool operator==(triface& t) {
          return tet == t.tet && loc == t.loc && ver == t.ver;
        }
        bool operator!=(triface& t) {
          return tet != t.tet || loc != t.loc || ver != t.ver;
        }
    };

    class face {

      public:

        shellface *sh;
        int shver;

        // Constructors;
        face() : sh(0), shver(0) {}
        // Operators;
        face& operator=(const face& s) {
          sh = s.sh; shver = s.shver;
          return *this;
        }
        bool operator==(face& s) {return (sh == s.sh) && (shver == s.shver);}
        bool operator!=(face& s) {return (sh != s.sh) || (shver != s.shver);}
    };

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The badface structure                                                     //
//                                                                           //
// A multiple usages structure. Despite of its name, a 'badface' can be used //
// to represent the following objects:                                       //
//   - a face of a tetrahedron which is (possibly) non-Delaunay;             //
//   - an encroached subsegment or subface;                                  //
//   - a bad-quality tetrahedron, i.e, has too large radius-edge ratio;      //
//   - a sliver, i.e., has good radius-edge ratio but nearly zero volume;    //
//   - a degenerate tetrahedron (see routine checkdegetet()).                //
//   - a recently flipped face (saved for undoing the flip later).           //
//                                                                           //
// It has the following fields:  'tt' holds a tetrahedron; 'ss' holds a sub- //
// segment or subface; 'cent' is the circumcent of 'tt' or 'ss', 'key' is a  //
// special value depending on the use, it can be either the square of the    //
// radius-edge ratio of 'tt' or the flipped type of 'tt';  'forg', 'fdest',  //
// 'fapex', and 'foppo' are vertices saved for checking the object in 'tt'   //
// or 'ss' is still the same when it was stored; 'noppo' is the fifth vertex //
// of a degenerate point set.  'previtem' and 'nextitem' implement a double  //
// link for managing many basfaces.                                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

    struct badface {
      triface tt; 
      face ss; 
      REAL key;
      REAL cent[3];
      point forg, fdest, fapex, foppo;
      point noppo;
      struct badface *previtem, *nextitem; 
    };

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The pbcdata structure                                                     //
//                                                                           //
// A pbcdata stores data of a periodic boundary condition defined on a pair  //
// of facets or segments. Let f1 and f2 define a pbcgroup. 'fmark' saves the //
// facet markers of f1 and f2;  'ss' contains two subfaces belong to f1 and  //
// f2, respectively.  Let s1 and s2 define a segment pbcgroup. 'segid' are   //
// the segment ids of s1 and s2; 'ss' contains two segments belong to s1 and //
// s2, respectively. 'transmat' are two transformation matrices. transmat[0] //
// transforms a point of f1 (or s1) into a point of f2 (or s2),  transmat[1] //
// does the inverse.                                                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

    struct pbcdata {
      int fmark[2];
      int segid[2];
      face ss[2];
      REAL transmat[2][4][4];
    };

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The list, link and queue data structures                                  //
//                                                                           //
// These data types are used to manipulate a set of (same-typed) data items. //
// For a given set S = {a, b, c, ...}, a list stores the elements of S in a  //
// piece of continuous memory. It allows quickly accessing each element of S,//
// thus is suitable for storing a fix-sized set.  While a link stores its    //
// elements incontinuously. It allows quickly inserting or deleting an item, //
// thus is suitable for storing a size-variable set.  A queue is basically a //
// special case of a link where one data element joins the link at the end   //
// and leaves in an ordered fashion at the other end.                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

    // The compfunc data type.  "compfunc" is a pointer to a linear-order
    //   function, which takes two 'void*' arguments and returning an 'int'. 
    //   
    // A function: int cmp(const T &, const T &),  is said to realize a
    //   linear order on the type T if there is a linear order <= on T such
    //   that for all x and y in T satisfy the following relation:
    //                 -1  if x < y.
    //   comp(x, y) =   0  if x is equivalent to y.
    //                 +1  if x > y.
    typedef int (*compfunc) (const void *, const void *);

    // The predefined compare functions for primitive data types.  They
    //   take two pointers of the corresponding date type, perform the
    //   comparation, and return -1, 0 or 1 indicating the default linear
    //   order of them.

    // Compare two 'integers'.
    static int compare_2_ints(const void* x, const void* y);
    // Compare two 'longs'. 
    static int compare_2_longs(const void* x, const void* y);
    // Compare two 'unsigned longs'. 
    static int compare_2_unsignedlongs(const void* x, const void* y);

    // The function used to determine the size of primitive data types and
    //   set the corresponding predefined linear order functions for them.
    static void set_compfunc(char* str, int* itembytes, compfunc* pcomp);

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// List data structure.                                                      //
//                                                                           //
// A 'list' is an array of items with automatically reallocation of memory.  //
// It behaves like an array.                                                 //
//                                                                           //
// 'base' is the starting address of the array;  The memory unit in list is  //
//   byte, i.e., sizeof(char). 'itembytes' is the size of each item in byte, //
//   so that the next item in list will be found at the next 'itembytes'     //
//   counted from the current position.                                      //
//                                                                           //
// 'items' is the number of items stored in list.  'maxitems' indicates how  //
//   many items can be stored in this list. 'expandsize' is the increasing   //
//   size (items) when the list is full.                                     //
//                                                                           //
// 'comp' is a pointer pointing to a linear order function for the list.     //
//   default it is set to 'NULL'.                                            //
//                                                                           //
// The index of list always starts from zero, i.e., for a list L contains    //
//   n elements, the first element is L[0], and the last element is L[n-1].  //
//   This feature lets lists likes C/C++ arrays.                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

    class list {

      public:

        char *base;
        int  itembytes;
        int  items, maxitems, expandsize;
        compfunc comp;

      public:

        list(int itbytes, compfunc pcomp, int mitems = 256, int exsize = 128) {
          listinit(itbytes, pcomp, mitems, exsize);
        }
        list(char* str, int mitems = 256, int exsize = 128) {
          set_compfunc(str, &itembytes, &comp);
          listinit(itembytes, comp, mitems, exsize);
        }
        ~list() { free(base); }

        void *operator[](int i) { return (void *) (base + i * itembytes); }

        void listinit(int itbytes, compfunc pcomp, int mitems, int exsize);
        void setcomp(compfunc compf) { comp = compf; }    
        void clear() { items = 0; }
        int  len() { return items; }
        void *append(void* appitem);
        void *insert(int pos, void* insitem);
        void del(int pos, int order);
        int  hasitem(void* checkitem);
        void sort();
    }; 

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Memorypool data structure.                                                //
//                                                                           //
// A type used to allocate memory.  (It is incorporated from Shewchuk's      //
// Triangle program)                                                         //
//                                                                           //
// firstblock is the first block of items. nowblock is the block from which  //
//   items are currently being allocated. nextitem points to the next slab   //
//   of free memory for an item. deaditemstack is the head of a linked list  //
//   (stack) of deallocated items that can be recycled.  unallocateditems is //
//   the number of items that remain to be allocated from nowblock.          //
//                                                                           //
// Traversal is the process of walking through the entire list of items, and //
//   is separate from allocation.  Note that a traversal will visit items on //
//   the "deaditemstack" stack as well as live items.  pathblock points to   //
//   the block currently being traversed.  pathitem points to the next item  //
//   to be traversed.  pathitemsleft is the number of items that remain to   //
//   be traversed in pathblock.                                              //
//                                                                           //
// itemwordtype is set to POINTER or FLOATINGPOINT, and is used to suggest   //
//   what sort of word the record is primarily made up of.  alignbytes       //
//   determines how new records should be aligned in memory.  itembytes and  //
//   itemwords are the length of a record in bytes (after rounding up) and   //
//   words.  itemsperblock is the number of items allocated at once in a     //
//   single block.  items is the number of currently allocated items.        //
//   maxitems is the maximum number of items that have been allocated at     //
//   once; it is the current number of items plus the number of records kept //
//   on deaditemstack.                                                       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

    class memorypool {

      public:

        void **firstblock, **nowblock;
        void *nextitem;
        void *deaditemstack;
        void **pathblock;
        void *pathitem;
        wordtype itemwordtype;
        int  alignbytes;
        int  itembytes, itemwords;
        int  itemsperblock;
        long items, maxitems;
        int  unallocateditems;
        int  pathitemsleft;

      public:

        memorypool();
        memorypool(int, int, enum wordtype, int);
        ~memorypool();
    
        void poolinit(int, int, enum wordtype, int);
        void restart();
        void *alloc();
        void dealloc(void*);
        void traversalinit();
        void *traverse();
    };  

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Link data structure.                                                      //
//                                                                           //
// A 'link' is a double linked nodes. It uses the memorypool data structure  //
// for memory management.  Following is an image of a link.                  //
//                                                                           //
//   head-> ____0____      ____1____      ____2____      _________<-tail     //
//         |__next___|--> |__next___|--> |__next___|--> |__NULL___|          //
//         |__NULL___|<-- |__prev___|<-- |__prev___|<-- |__prev___|          //
//         |         |    |_       _|    |_       _|    |         |          //
//         |         |    |_ Data1 _|    |_ Data2 _|    |         |          //
//         |_________|    |_________|    |_________|    |_________|          //
//                                                                           //
// The unit size for storage is size of pointer, which may be 4-byte (in 32- //
//   bit machine) or 8-byte (in 64-bit machine). The real size of an item is //
//   stored in 'linkitembytes'.                                              //
//                                                                           //
// 'head' and 'tail' are pointers pointing to the first and last nodes. They //
//   do not conatin data (See above).                                        //
//                                                                           //
// 'nextlinkitem' is a pointer pointing to a node which is the next one will //
//   be traversed. 'curpos' remembers the position (1-based) of the current  //
//   traversing node.                                                        //
//                                                                           //
// 'linkitems' indicates how many items in link. Note it is different with   //
//   'items' of memorypool.                                                  //
//                                                                           //
// The index of link starts from 1, i.e., for a link K contains n elements,  //
//   the first element of the link is K[1], and the last element is K[n].    //
//   See the above figure.                                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

    class link : public memorypool {

      public:

        void **head, **tail;
        void *nextlinkitem;
        int  linkitembytes;
        int  linkitems;
        int  curpos;
        compfunc comp;

      public:

        link(int _itembytes, compfunc _comp, int itemcount) {
          linkinit(_itembytes, _comp, itemcount);
        }
        link(char* str, int itemcount) {
          set_compfunc(str, &linkitembytes, &comp);
          linkinit(linkitembytes, comp, itemcount);
        }

        void linkinit(int _itembytes, compfunc _comp, int itemcount);
        void setcomp(compfunc compf) { comp = compf; }
        void rewind() { nextlinkitem = *head; curpos = 1; }
        void goend() { nextlinkitem = *(tail + 1); curpos = linkitems; }    
        long len() { return linkitems; }
        void clear();
        bool move(int numberofnodes);
        bool locate(int pos);
        void *add(void* newitem);
        void *insert(int pos, void* insitem);
        void *del(void* delitem);
        void *del(int pos);
        void *getitem();
        void *getnitem(int pos);
        int  hasitem(void* checkitem);
    };

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Queue data structure.                                                     //
//                                                                           //
// A 'queue' is a basically a link.  Following is an image of a queue.       //
//              ___________     ___________     ___________                  //
//   Pop() <-- |_         _|<--|_         _|<--|_         _| <-- Push()      //
//             |_  Data0  _|   |_  Data1  _|   |_  Data2  _|                 //
//             |___________|   |___________|   |___________|                 //
//              queue head                       queue tail                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

    class queue : public link {

      public:

        queue(int bytes, int count = 256) : link(bytes, NULL, count) {}
        queue(char* str, int count = 256) : link(str, count) {}

        int  empty() { return linkitems == 0; }
        void *push(void* newitem) { return link::add(newitem); }
        void *bot() { return link::getnitem(1); }
        void *pop() { return link::del(1); }
    };

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Following are variables used in 'tetgenmesh' for miscellaneous purposes.  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

    // Pointer to an object of 'tetgenio', which contains input data.
    tetgenio *in;

    // Pointer to an object of 'tetgenbehavor', which contains the user-
    //   defined command line swithes and filenames.
    tetgenbehavior *b;

    // Variables used to allocate and access memory for tetrahedra, subfaces
    //   subsegments, points, encroached subfaces, encroached subsegments,
    //   bad-quality tetrahedra, and so on.
    memorypool *tetrahedrons;
    memorypool *subfaces;
    memorypool *subsegs;
    memorypool *points;   
    memorypool *badsubsegs;
    memorypool *badsubfaces;
    memorypool *badtetrahedrons;
    memorypool *flipstackers;

    // Pointer to the 'tetrahedron' that occupies all of "outer space".
    tetrahedron *dummytet;
    tetrahedron *dummytetbase; // Keep base address so we can free it later.

    // Pointer to the omnipresent subface.  Referenced by any tetrahedron,
    //   or subface that isn't connected to a subface at that location.
    shellface *dummysh;
    shellface *dummyshbase;    // Keep base address so we can free it later.

    // A point above the plane in which the facet currently being used lies.
    //   It is used as a reference point for orient3d().
    point *facetabovepointarray, abovepoint;

    // Array (size = numberoftetrahedra * 6) for storing high-order nodes of
    //   tetrahedra (only used when -o2 switch is selected).
    point *highordertable;

    // Arrays for storing and searching pbc data. 'subpbcgrouptable', (size
    //   is numberofpbcgroups) for pbcgroup of subfaces. 'segpbcgrouptable',
    //   a list for pbcgroup of segments. Because a segment can have several
    //   pbcgroup incident on it, its size is unknown on input, it will be
    //   found in 'createsegpbcgrouptable()'.
    pbcdata *subpbcgrouptable;
    list *segpbcgrouptable;
    // A map for searching the pbcgroups of a given segment. 'idx2segpglist'
    //   (size = number of input segments + 1), and 'segpglist'.  
    int *idx2segpglist, *segpglist;

    // List used in Delaunay refinement for keeping tetrahedra which need to
    //   be classified by the quality measure.
    list *qualchecktetlist;

    // Array used in Delaunay refinement for keeping the protecting sphere
    //   radii of the input vertices.
    REAL *rpsarray;

    // Queues that maintain the bad (badly-shaped or too large) tetrahedra.
    //   The tails are pointers to the pointers that have to be filled in to
    //   enqueue an item.  The queues are ordered from 63 (highest priority)
    //   to 0 (lowest priority).
    badface *subquefront[2], **subquetail[2];
    badface *tetquefront[64], **tetquetail[64];

    // Pointer to a recently visited tetrahedron. Improves point location
    //   if proximate points are inserted sequentially.
    triface recenttet;

    REAL xmax, xmin, ymax, ymin, zmax, zmin;      // Bounding box of points.
    REAL longest;                       // The longest possible edge length.
    long hullsize;                        // Number of faces of convex hull.
    long insegment;                             // Number of input segments.
    int steinerleft;               // Number of Steiner points not yet used.
    int edgeboundindex;              // Index to find edge bound of a point.
    int pointmarkindex;         // Index to find boundary marker of a point.
    int point2simindex;      // Index to find a simplex adjacent to a point.
    int point2pbcptindex;           // Index to find a pbc point to a point.
    int highorderindex; // Index to find extra nodes for highorder elements.
    int elemattribindex;       // Index to find attributes of a tetrahedron.
    int volumeboundindex;    // Index to find volume bound of a tetrahedron.
    int shmarkindex;          // Index to find boundary marker of a subface.
    int areaboundindex;            // Index to find area bound of a subface.
    int checksubfaces;                // Are there subfaces in the mesh yet?
    int checkpbcs;                // Are there periodic boundary conditions?
    int varconstraint;  // Are there variant (node, seg, facet) constraints?
    int nonconvex;                            // Is current mesh non-convex?
    int dupverts;                          // Are there duplicated vertices?
    int unuverts;                              // Are there unused vertices?
    int relverts;                       // The number of relocated vertices.
    int jettisoninverts;         // The number of jettisoned input vertices.
    int symbolic;                             // Use symbolic insphere test.
    long samples;            // Number of random samples for point location.
    unsigned long randomseed;                 // Current random number seed.
    REAL macheps;                                    // The machine epsilon.
    REAL cosmaxdihed, cosmindihed; // The cosine values of max/min dihedral.
    int maxcavfaces, maxcavverts;         // The size of the largest cavity.
    int expcavcount;                      // The times of expanding cavitys.
    long abovecount;                   // Number of abovepoints calculation.
    long rejectsmalledgecount;       // Number of rejections of small edges.
    long bowatvolcount, bowatsubcount, bowatsegcount;     // Bowyer-Watsons.
    long r1count, r2count, r3count, r4count;   // Number of rules performed.
    long flip23s, flip32s, flip22s, flip44s;   // Number of flips performed.

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Fast lookup tables for mesh manipulation primitives.                      //
//                                                                           //
// Mesh manipulation primitives (given below) are basic operations on mesh   //
// data structures. They answer basic queries on mesh handles, such as "what //
// is the origin (or destination, or apex) of the face?", "what is the next  //
// (or previous) edge in the edge ring?", and "what is the next face in the  //
// face ring?", and so on.                                                   //
//                                                                           //
// The implementation of teste basic queries can take advangtage of the fact //
// that the mesh data structures additionally store geometric informations.  //
// For example, we have ordered the 4 vertices (from 0 to 3) and the 4 faces //
// (from 0 to 3) of a tetrahedron,  and for each face of the tetrahedron, a  //
// sequence of vertices has stipulated,  therefore the origin of any face of //
// the tetrahedron can be quickly determined by a table 'locver2org', which  //
// takes the index of the face and the edge version as inputs.  A list of    //
// fast lookup tables are defined below. They're just like global variables. //
// These tables are initialized at the runtime.                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

    // For enext() primitive, uses 'ver' as the index. 
    static int ve[6];

    // For org(), dest() and apex() primitives, uses 'ver' as the index.
    static int vo[6], vd[6], va[6];

    // For org(), dest() and apex() primitives, uses 'loc' as the first
    //   index and 'ver' as the second index.
    static int locver2org[4][6];
    static int locver2dest[4][6];
    static int locver2apex[4][6];

    // For oppo() primitives, uses 'loc' as the index.
    static int loc2oppo[4];

    // For fnext() primitives, uses 'loc' as the first index and 'ver' as
    //   the second index,  returns an array containing a new 'loc' and a
    //   new 'ver'. Note: Only valid for 'ver' equals one of {0, 2, 4}.
    static int locver2nextf[4][6][2];

    // For enumerating three edges of a triangle.
    static int plus1mod3[3];
    static int minus1mod3[3];

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Mesh manipulation primitives                                              //
//                                                                           //
// A serial of mesh operations such as topological maintenance,  navigation, //
// local modification, etc.,  is accomplished through a set of mesh manipul- //
// ation primitives. These primitives are indeed very simple functions which //
// take one or two handles ('triface's and 'face's) as parameters,  perform  //
// basic operations such as "glue two tetrahedra at a face",  "return the    //
// origin of a tetrahedron", "return the subface adjoining at the face of a  //
// tetrahedron", and so on.                                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

        // Primitives for tetrahedra.
    inline void decode(tetrahedron ptr, triface& t);
    inline tetrahedron encode(triface& t);
    inline void sym(triface& t1, triface& t2);
    inline void symself(triface& t);
    inline void bond(triface& t1, triface& t2);
    inline void dissolve(triface& t);
    inline point org(triface& t);
    inline point dest(triface& t);
    inline point apex(triface& t);
    inline point oppo(triface& t);
    inline void setorg(triface& t, point pointptr);
    inline void setdest(triface& t, point pointptr);
    inline void setapex(triface& t, point pointptr);
    inline void setoppo(triface& t, point pointptr);
    inline void esym(triface& t1, triface& t2);
    inline void esymself(triface& t);
    inline void enext(triface& t1, triface& t2);
    inline void enextself(triface& t);
    inline void enext2(triface& t1, triface& t2);
    inline void enext2self(triface& t);
    inline bool fnext(triface& t1, triface& t2);
    inline bool fnextself(triface& t);
    inline void enextfnext(triface& t1, triface& t2);
    inline void enextfnextself(triface& t);
    inline void enext2fnext(triface& t1, triface& t2);
    inline void enext2fnextself(triface& t);
    inline void infect(triface& t);
    inline void uninfect(triface& t);
    inline bool infected(triface& t);
    inline REAL elemattribute(tetrahedron* ptr, int attnum);
    inline void setelemattribute(tetrahedron* ptr, int attnum, REAL value);
    inline REAL volumebound(tetrahedron* ptr);
    inline void setvolumebound(tetrahedron* ptr, REAL value);
 
    // Primitives for subfaces and subsegments.
    inline void sdecode(shellface sptr, face& s);
    inline shellface sencode(face& s);
    inline void spivot(face& s1, face& s2);
    inline void spivotself(face& s);
    inline void sbond(face& s1, face& s2);
    inline void sbond1(face& s1, face& s2);
    inline void sdissolve(face& s);
    inline point sorg(face& s);
    inline point sdest(face& s);
    inline point sapex(face& s);
    inline void setsorg(face& s, point pointptr);
    inline void setsdest(face& s, point pointptr);
    inline void setsapex(face& s, point pointptr);
    inline void sesym(face& s1, face& s2);
    inline void sesymself(face& s);
    inline void senext(face& s1, face& s2);
    inline void senextself(face& s);
    inline void senext2(face& s1, face& s2);
    inline void senext2self(face& s);
    inline void sfnext(face&, face&);
    inline void sfnextself(face&);
    inline badface* shell2badface(face& s);
    inline void setshell2badface(face& s, badface* value);
    inline REAL areabound(face& s);
    inline void setareabound(face& s, REAL value);
    inline int shellmark(face& s);
    inline void setshellmark(face& s, int value);
    inline enum shestype shelltype(face& s);
    inline void setshelltype(face& s, enum shestype value); 
    inline int shellpbcgroup(face& s);
    inline void setshellpbcgroup(face& s, int value);
    inline void sinfect(face& s);
    inline void suninfect(face& s);
    inline bool sinfected(face& s);

    // Primitives for interacting tetrahedra and subfaces.
    inline void tspivot(triface& t, face& s);
    inline void stpivot(face& s, triface& t);
    inline void tsbond(triface& t, face& s);
    inline void tsdissolve(triface& t);    
    inline void stdissolve(face& s);

    // Primitives for interacting subfaces and subsegs.
    inline void sspivot(face& s, face& edge);
    inline void ssbond(face& s, face& edge);
    inline void ssdissolve(face& s);

    // Primitives for points.
    inline int  pointmark(point pt);
    inline void setpointmark(point pt, int value);
    inline enum verttype pointtype(point pt);
    inline void setpointtype(point pt, enum verttype value);
    inline tetrahedron point2tet(point pt);
    inline void setpoint2tet(point pt, tetrahedron value);
    inline shellface point2sh(point pt);
    inline void setpoint2sh(point pt, shellface value);
    inline point point2ppt(point pt);
    inline void setpoint2ppt(point pt, point value);
    inline point point2pbcpt(point pt);
    inline void setpoint2pbcpt(point pt, point value);
    inline REAL edgebound(point pt);
    inline void setedgebound(point pt, REAL value);

    // Advanced primitives.
    inline void adjustedgering(triface& t, int direction);
    inline void adjustedgering(face& s, int direction);
    inline bool isdead(triface* t);
    inline bool isdead(face* s);
    inline bool isfacehaspoint(triface* t, point testpoint);
    inline bool isfacehaspoint(face* t, point testpoint);
    inline bool isfacehasedge(face* s, point tend1, point tend2);
    inline bool issymexist(triface* t);
    bool getnextface(triface*, triface*);
    void getnextsface(face*, face*);
    void tsspivot(triface*, face*);
    void sstpivot(face*, triface*);   
    bool findorg(triface* t, point dorg);
    bool findorg(face* s, point dorg);
    void findedge(triface* t, point eorg, point edest);
    void findedge(face* s, point eorg, point edest);
    void findface(triface *fface, point forg, point fdest, point fapex);
    void getonextseg(face* s, face* lseg);
    void getseghasorg(face* sseg, point dorg);
    point getsubsegfarorg(face* sseg);
    point getsubsegfardest(face* sseg);
    void printtet(triface*);
    void printsh(face*);

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Triangle-triangle intersection test                                       //
//                                                                           //
// The triangle-triangle intersection test is implemented with exact arithm- //
// etic. It exactly tells whether or not two triangles in three dimensions   //
// intersect.  Before implementing this test myself,  I tried two C codes    //
// (implemented by Thomas Moeller and Philippe Guigue, respectively), which  //
// are all public available. However both of them failed frequently. Another //
// unconvenience is both codes only tell whether or not the two triangles    //
// intersect without distinguishing the cases whether they exactly intersect //
// in interior or they just share a vertex or share an edge. The two latter  //
// cases are acceptable and should return not intersection in TetGen.        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

    enum interresult edge_vert_col_inter(REAL*, REAL*, REAL*);
    enum interresult edge_edge_cop_inter(REAL*, REAL*, REAL*, REAL*, REAL*);
    enum interresult tri_vert_cop_inter(REAL*, REAL*, REAL*, REAL*, REAL*);
    enum interresult tri_edge_cop_inter(REAL*, REAL*, REAL*, REAL*, REAL*,
                                        REAL*);
    enum interresult tri_edge_inter_tail(REAL*, REAL*, REAL*, REAL*, REAL*,
                                         REAL, REAL);
    enum interresult tri_edge_inter(REAL*, REAL*, REAL*, REAL*, REAL*);
    enum interresult tri_tri_inter(REAL*, REAL*, REAL*, REAL*, REAL*, REAL*);

    // Geometric predicates
    REAL insphere_sos(REAL*, REAL*, REAL*, REAL*, REAL*, int, int, int, int,
                      int);
    bool iscollinear(REAL*, REAL*, REAL*, REAL eps);
    bool iscoplanar(REAL*, REAL*, REAL*, REAL*, REAL vol6, REAL eps);
    bool iscospheric(REAL*, REAL*, REAL*, REAL*, REAL*, REAL vol24, REAL eps);

    // Linear algebra functions
    inline REAL dot(REAL* v1, REAL* v2);
    inline void cross(REAL* v1, REAL* v2, REAL* n);
    bool lu_decmp(REAL lu[4][4], int n, int* ps, REAL* d, int N);
    void lu_solve(REAL lu[4][4], int n, int* ps, REAL* b, int N);

    // Geometric quantities calculators.
    inline REAL distance(REAL* p1, REAL* p2);
    REAL shortdistance(REAL* p, REAL* e1, REAL* e2);
    REAL shortdistance(REAL* p, REAL* e1, REAL* e2, REAL* e3);
    REAL interiorangle(REAL* o, REAL* p1, REAL* p2, REAL* n);
    void projpt2edge(REAL* p, REAL* e1, REAL* e2, REAL* prj);
    void projpt2face(REAL* p, REAL* f1, REAL* f2, REAL* f3, REAL* prj);
    void facenormal(REAL* pa, REAL* pb, REAL* pc, REAL* n, REAL* nlen);
    void edgeorthonormal(REAL* e1, REAL* e2, REAL* op, REAL* n);
    REAL facedihedral(REAL* pa, REAL* pb, REAL* pc1, REAL* pc2);
    void tetalldihedral(point, point, point, point, REAL dihed[6]);
    void tetallnormal(point, point, point, point, REAL N[4][3], REAL* volume);
    bool circumsphere(REAL*, REAL*, REAL*, REAL*, REAL* cent, REAL* radius);
    void inscribedsphere(REAL*, REAL*, REAL*, REAL*, REAL* cent, REAL* radius);
    void rotatepoint(REAL* p, REAL rotangle, REAL* p1, REAL* p2);
    void spherelineint(REAL* p1, REAL* p2, REAL* C, REAL R, REAL p[7]);
    void linelineint(REAL *p1,REAL *p2, REAL *p3, REAL *p4, REAL p[7]);
    void planelineint(REAL*, REAL*, REAL*, REAL*, REAL*, REAL*, REAL*);

    // Memory managment routines.
    void dummyinit(int, int);
    void initializepools();
    void tetrahedrondealloc(tetrahedron*);
    tetrahedron *tetrahedrontraverse();
    void shellfacedealloc(memorypool*, shellface*);
    shellface *shellfacetraverse(memorypool*);
    void badfacedealloc(memorypool*, badface*);
    badface *badfacetraverse(memorypool*);
    void pointdealloc(point);
    point pointtraverse();
    void maketetrahedron(triface*);
    void makeshellface(memorypool*, face*);
    void makepoint(point*);

    // Mesh items searching routines.
    void makepoint2tetmap();
    void makeindex2pointmap(point*& idx2verlist);
    void makesegmentmap(int*& idx2seglist, shellface**& segsperverlist);
    void makesubfacemap(int*& idx2facelist, shellface**& facesperverlist);
    void maketetrahedronmap(int*& idx2tetlist, tetrahedron**& tetsperverlist);

    // Point location routines.
    unsigned long randomnation(unsigned int choices);
    REAL distance2(tetrahedron* tetptr, point p);
    enum locateresult preciselocate(point searchpt, triface* searchtet, long);
    enum locateresult locate(point searchpt, triface* searchtet);
    enum locateresult adjustlocate(point, triface*, enum locateresult, REAL);
    enum locateresult locatesub(point searchpt, face* searchsh, int stopatseg);
    enum locateresult adjustlocatesub(point, face*, enum locateresult, REAL);
    enum locateresult locateseg(point searchpt, face* searchseg);
    enum locateresult adjustlocateseg(point, face*, enum locateresult, REAL);

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Mesh Local Transformation Operators                                       //
//                                                                           //
// These operators (including flips, insert & remove vertices and so on) are //
// used to transform (or replace) a set of mesh elements into another set of //
// mesh elements.                                                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

    // Mesh transformation routines.
    enum fliptype categorizeface(triface& horiz);
    void enqueueflipface(triface& checkface, queue* flipqueue);
    void enqueueflipedge(face& checkedge, queue* flipqueue);
    void flip23(triface* flipface, queue* flipqueue);
    void flip32(triface* flipface, queue* flipqueue);
    void flip22(triface* flipface, queue* flipqueue);
    void flip22sub(face* flipedge, queue* flipqueue);
    long flip(queue* flipqueue, badface **plastflip, bool, bool, bool);
    void undoflip(badface *lastflip);
    long flipsub(queue* flipqueue);

    void splittetrahedron(point newpoint, triface* splittet, queue* flipqueue);
    void unsplittetrahedron(triface* splittet);
    void splittetface(point newpoint, triface* splittet, queue* flipqueue);
    void unsplittetface(triface* splittet);
    void splitsubface(point newpoint, face* splitface, queue* flipqueue);
    void unsplitsubface(face* splitsh);
    void splittetedge(point newpoint, triface* splittet, queue* flipqueue);
    void unsplittetedge(triface* splittet);
    void splitsubedge(point newpoint, face* splitsh, queue* flipqueue);
    void unsplitsubedge(face* splitsh);
    enum insertsiteresult insertsite(point newpoint, triface* searchtet,
                                     bool approx, queue* flipqueue);
    void undosite(enum insertsiteresult insresult, triface* splittet, 
                  point torg, point tdest, point tapex, point toppo);

    // Delaunay tetrahedralization routines.
    void formstarpolyhedron(point pt, list* tetlist, list* verlist);
    void collectcavtets(point newpoint, list* cavtetlist);
    bool unifypoint(point testpt, triface*, enum locateresult, REAL);
    void closeopenface(triface* openface, queue* flipque);
    void inserthullsite(point inspoint, triface* horiz, queue* flipque);
    void incrflipdelaunay(triface*, point*, long, bool, bool, REAL, queue*);
    long delaunizevertices();

    // Surface triangulation routines.
    void formstarpolygon(point pt, list* trilist, list* verlist);
    void getfacetabovepoint(face* facetsh);
    void collectcavsubs(point newpoint, list* cavsublist);
    void collectvisiblesubs(int shmark, point inspoint, face* horiz, queue*);
    void incrflipdelaunaysub(int shmark, REAL eps, list*, int, REAL*, queue*);
    enum finddirectionresult finddirectionsub(face* searchsh, point tend);
    void insertsubseg(face* tri);
    bool scoutsegmentsub(face* searchsh, point tend);
    void flipedgerecursive(face* flipedge, queue* flipqueue);
    void constrainededge(face* startsh, point tend, queue* flipqueue);
    void recoversegment(point tstart, point tend, queue* flipqueue);
    void infecthullsub(memorypool* viri);
    void plaguesub(memorypool* viri);
    void carveholessub(int holes, REAL* holelist, memorypool* viri);
    void triangulate(int shmark, REAL eps, list* ptlist, list* conlist,
                     int holes, REAL* holelist, memorypool* viri, queue*);
    void retrievefacetsubs(list* newshlist, list* boundedgelist);
    void unifysegments();
    void mergefacets(queue* flipqueue);    
    void removefacetvertex(face* remsh, list* oldshlist, list* ptlist,
                           list* conlist, list* newshlist, list* boundedgelist,
                           memorypool* viri, queue* flipque);
    void removesegmentvertex(face* remseg, face* nremseg, face* newseg,
                             list* spinshlist, list* ptlist, list* conlist,
                             list* newshlist, list* boundedgelist,
                             list* newsegshlist, list* newshlistlist,
                             list* oldshlistlist, memorypool*, queue*);
    void removevolumevertex(point rempt, list* ptlist, list* newshlist,
                            list* oldtetlist, list* newtetlist,
                            list* frontlist, list* misfrontlist, queue*);
    void removefreepoints(list*, list*, memorypool*, queue*);
    long meshsurface();

    // Detect intersecting facets of PLC.
    void interecursive(shellface** subfacearray, int arraysize, int axis,
                       REAL bxmin, REAL bxmax, REAL bymin, REAL bymax,
                       REAL bzmin, REAL bzmax, int* internum);
    void detectinterfaces(); 

    // Periodic boundary condition supporting routines.
    void createsubpbcgrouptable();
    void getsubpbcgroup(face* pbcsub, pbcdata** pd, int *f1, int *f2);
    enum locateresult getsubpbcsympoint(point, face*, point, face*);
    void createsegpbcgrouptable();
    enum locateresult getsegpbcsympoint(point, face*, point, face*, int);

    // Vertex perturbation routines.
    REAL randgenerator(REAL range);
    bool checksub4cocir(face* testsub, REAL eps, bool once, bool enqflag);
    bool checktet4cosph(triface* testtet, REAL eps, bool once, bool enqflag);
    void tallcocirsubs(REAL eps, bool enqflag);
    void tallcosphtets(list* testtetlist, REAL eps, bool enqflag);
    bool tallencsegsfsubs(point testpt, list* cavsublist);
    bool tallencsubsfsubs(point testpt, list* cavtetlist);
    void collectflipedges(point inspoint, face* splitseg, queue* flipqueue);
    void perturbrepairencsegs(queue* flipqueue);
    void perturbrepairencsubs(REAL eps, list* cavsublist, queue* flipqueue);
    void perturbrepairbadtets(REAL eps, list* cavsublist, list*, queue*);
    void incrperturbvertices(REAL eps);

    // Segment recovery routines.
    void markacutevertices(REAL acuteangle);
    void initializerpsarray();
    enum finddirectionresult finddirection(triface* searchtet, point, long);
    void getsearchtet(point p1, point p2, triface* searchtet, point* tend);
    bool isedgeencroached(point p1, point p2, point testpt, bool degflag);
    point scoutrefpoint(triface* searchtet, point tend);
    point getsegmentorigin(face* splitseg);
    point getsplitpoint(face* splitseg, point refpoint);
    void delaunizesegments();

    // Facets recovery routines.
    bool insertsubface(face* insertsh, triface* searchtet);
    bool tritritest(triface* checktet, point p1, point p2, point p3);
    void initializecavity(list* floorlist, list* ceillist, list* frontlist);
    void delaunizecavvertices(triface*, list*, list*, list*, queue*);
    void retrievenewtets(list* newtetlist);
    void insertauxsubface(triface* front, triface* idfront);
    bool scoutfront(triface* front, triface* idfront, list* newtetlist);
    void gluefronts(triface* front, triface* front1);
    bool identifyfronts(list* frontlist, list* misfrontlist, list* newtetlist);
    void detachauxsubfaces(list* newtetlist);
    void expandcavity(list* frontlist, list* misfrontlist, list* newtetlist,
                      list* crosstetlist, queue* missingshqueue, queue*);
    void carvecavity(list* newtetlist, list* outtetlist, queue* flipque);
    void delaunizecavity(list* floorlist, list* ceillist, list* ceilptlist,
                         list* floorptlist, list* frontlist,list* misfrontlist,
                         list* newtetlist, list* crosstetlist, queue*, queue*);
    void formmissingregion(face* missingsh, list* missingshlist,
                           list* equatptlist, int* worklist);
    void formcavity(list* missingshlist, list* crossedgelist, 
                    list* equatptlist, list* crossshlist, list* crosstetlist,
                    list* belowfacelist, list* abovefacelist,
                    list* horizptlist, list* belowptlist, list* aboveptlist,
                    queue* missingshqueue, int* worklist);
    bool scoutcrossingedge(list* missingshlist, list* boundedgelist,
                           list* crossedgelist, int* worklist);
    void rearrangesubfaces(list* missingshlist, list* boundedgelist,
                           list* equatptlist, int* worklist);
    void insertallsubfaces(queue* missingshqueue);
    void constrainedfacets();

    // Carving out holes and concavities routines.
    void infecthull(memorypool *viri);
    void plague(memorypool *viri);
    void regionplague(memorypool *viri, REAL attribute, REAL volume);
    void removeholetets(memorypool *viri);
    void assignregionattribs();
    void carveholes();

    // Steiner points removing routines.
    void fillsteinerpolygon(list* oldshlist, list* boundedgelist);
    void orientnewsubfaces(list* newshlist, face* orientsh, REAL* norm);
    bool constrainedflip(triface* flipface, triface* front, queue* flipque);
    bool recoverfront(triface* front, list* newtetlist, queue* flipque);
    void repairflips(queue* flipque);
    bool constrainedcavity(triface* oldtet, list* floorlist, list* ceillist,
                           list* ptlist, list* frontlist, list* misfrontlist,
                           list* newtetlist, queue* flipque);
    void expandsteinercavity(point steinpt, REAL eps, list* frontlist, list*);
    bool relocatesteiner(point sp, point np, REAL* n, list* frontlist, list*);
    void bowatinsertsite(point steinpt, triface* oldtet, list*, list*, queue*);
    void removesubsteiner(face* steinsh, list* oldshlist, list* ptlist,
                          list* conlist, list* newshlist, list* boundedgelist,
                          list* oldtetlist, list* newtetlist, list* frontlist,
                          list*, memorypool*, queue*);
    void removesegsteiner(face* steinseg, list* spinshlist, list* ptlist,
                          list* conlist, list* newshlist, list* boundedgelist,
                          list* newsegshlist, list* newshlistlist,
                          list* oldshlistlist, list* oldtetlist,
                          list* newtetlist, list* frontlist,
                          list* misfrontlist, memorypool *viri, queue*);
    void removesteiners();

    // Mesh reconstruction rotuines.
    void assignvarconstraints(point *idx2verlist);
    long reconstructmesh();
    void insertaddpoints();

    // Mesh repair routines.
    bool checktet4dege(triface* testtet, REAL eps, list* degetetlist);
    bool finddiagonal(triface* akite);
    void removetetbypeeloff(triface *badtet, queue* flipqueue);
    void removetetbyflip32(triface *badtet, queue* flipqueue);
    bool removetetbyedgereduce(triface *badtet, queue* flipqueue);
    bool removekite(triface* akite, bool reduce, queue* flipqueue);
    void talldegetets(REAL eps, list* newtetlist, list* degetetlist);
    void repairdegetets(REAL eps, list* newtetlist, queue* flipqueue);

    // Delaunay refinement routines.
    void marksharpsubfaces(REAL dihedbound);
    void markskinnysubfaces(REAL anglebound);
    badface* dequeuebadtet();
    bool checkseg4encroach(face* testseg, point testpt, point*, bool enqflag);
    bool checksub4encroach(face* testsub, point testpt, bool enqflag);
    bool checkseg4badqual(face* testseg, bool enqflag);
    bool checksub4badqual(face* testsub, bool enqflag);
    bool checktet4badqual(triface* testtet, bool enqflag);
    bool checktet4sliver(triface* testtet, bool enqflag);
    bool checkseg4splitting(face* testseg, point* pencpt);
    bool checksub4splitting(face* testsub, bool bqual);
    bool tallencsegs(point testpt, list *cavtetlist, list* cavsublist);
    bool tallencsubs(point testpt, list *cavtetlist, list* cavsublist);
    void tallencsubsfseg(face* testseg);
    void tallbadsegs();
    void tallbadsubs();
    void tallbadtetrahedrons();
    void tallslivers();
    void doqualchecktetlist();
    void getsplitpoint1(face* splitseg, REAL* splitpt, point encpt);
    void collectvolcavtets(point, list*, list*);
    void collectsubcavtets(point, face*, list*, list*);
    void collectsegcavtets(point, face*, list*, list*, list*);
    void collectbowatcavsubs(face*, list*, list*);
    void formbowatcavity(point newpt, list* cavtetlist, list* ceillist,
                         list* cutceillist);
    void formbowatsubcavity(point newpt, list* spinshlist, list** cavsublists,
                            list** subceillists, bool cutchk);
    void bowatinsertsegsite(point newpt, face* splitseg, list* spinshlist);
    void bowatinsertsubsite(point newpt, face* splitseg, face* splitsegsh,
                            list* cavsublist, list* subceillist,
                            list* ceillist, bool chksubs);
    void bowatinsertvolsite(point newpt, face* splitseg, list* cavtetlist,
                            list* ceillist, list* spinshlist,
                            list** cavsublists, list** subceillists,
                            bool chksubs, bool chkbqual, bool chksliver);
    void repairencsegs(bool chksubs, bool chkbqual);
    void repairencsegs(list* cavtetlist, list* ceillist, list* cutceillist,
                       list* spinshlist, list* quadtetlist, bool chksubs,
                       bool chkbqual);
    void repairencsubs(list* cavtetlist, list* ceillist, list* cutceillist,
                       list* spinshlist, list* quadtetlist, bool chksubs);
    void repairbadtets(list* cavtetlist, list* ceillist, list* cutceillist,
                       list* spinshlist, list* quadtetlist);
    void repairslivers(list* cavtetlist, list* ceillist, list* cutceillist,
                       list* spinshlist);
    void enforcequality();

    // I/O routines
    void transfernodes();
    void jettisonnodes();
    void highorder();
    void outnodes(tetgenio* out);
    void outelements(tetgenio* out);
    void outfaces(tetgenio* out);
    void outhullfaces(tetgenio* out);
    void outsubfaces(tetgenio* out);
    void outsubsegments(tetgenio* out);
    void outneighbors(tetgenio* out);
    void outpbcnodes(tetgenio* out);
    void outsmesh(char* smfilename);
    void outmesh2medit(char* mfilename);
    void outmesh2gid(char* gfilename);
    void outmesh2off(char* ofilename);

    // User interaction routines.
    void internalerror();
    void checkmesh();
    void checkshells();
    void checkdelaunay(REAL eps, queue* flipqueue);
    void checkdegeneracy(REAL eps);
    void checkconforming();
    void algorithmicstatistics();
    void qualitystatistics();
    void statistics();

  public:

    // Constructor and destructor.
    tetgenmesh();
    ~tetgenmesh();

};                                               // End of class tetgenmesh.

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// tetrahedralize()    Interface for using TetGen's library to generate      //
//                     Delaunay tetrahedralizations, constrained Delaunay    //
//                     tetrahedralizations, quality tetrahedral meshes.      //
//                                                                           //
// Two functions (interfaces) are available.  The difference is only the way //
// of passing switches.  One directly accepts an object of 'tetgenbehavior', //
// while the other accepts a string which is the same as one can used in the //
// command line.  The latter may be more convenient for users who don't want //
// to kown the 'tetgenbehavir' structure.                                    //
//                                                                           //
// 'in' is an object of 'tetgenio' which contains a PLC you want to tetrahed-//
// ralize or a previously generated tetrahedral mesh you want to refine.  It //
// must not be a NULL. 'out' is another object of 'tetgenio' for storing the //
// generated tetrahedral mesh. It can be a NULL. If so, the output will be   //
// saved to file(s).                                                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void tetrahedralize(tetgenbehavior *b, tetgenio *in, tetgenio *out);
void tetrahedralize(char *switches, tetgenio *in, tetgenio *out);

#endif // #ifndef tetgenH
