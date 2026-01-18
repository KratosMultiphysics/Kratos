/* =============================================================================
**  This file is part of the mmg software package for the tetrahedral
**  mesh modification.
**  Copyright (c) Bx INP/CNRS/Inria/UBordeaux/UPMC, 2004-
**
**  mmg is free software: you can redistribute it and/or modify it
**  under the terms of the GNU Lesser General Public License as published
**  by the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  mmg is distributed in the hope that it will be useful, but WITHOUT
**  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
**  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
**  License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License and of the GNU General Public License along with mmg (in
**  files COPYING.LESSER and COPYING). If not, see
**  <http://www.gnu.org/licenses/>. Please read their terms carefully and
**  use this copy of the mmg distribution only if you accept them.
** =============================================================================
*/

/**
 * \file common/libmmgtypes.h
 * \ingroup API
 * \brief Types used throughout the Mmg libraries
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \date 01 2014
 * \copyright GNU Lesser General Public License.
 */
#include <stdint.h>
#include <stdarg.h>
#include <stddef.h>

#include "mmg/common/mmgcmakedefines.h"
#include "mmg/common/mmgversion.h"

#ifndef _LIBMMGTYPES_H
#define _LIBMMGTYPES_H

/**
 * \def MMG5_SUCCESS
 *
 * Return value for success.
 *
 */
#define MMG5_SUCCESS       0
/**
 * \def MMG5_LOWFAILURE
 *
 * Return value if the remesh process failed but we can save a conform
 * mesh.
 *
 */
#define MMG5_LOWFAILURE    1
/**
 * \def MMG5_STRONGFAILURE
 *
 * Return value if the remesh process failed and the mesh is
 * non-conform.
 *
 */
#define MMG5_STRONGFAILURE 2

/**
 * Implicit boundary in iso mode
 *
 */
#define MG_ISO    10

/**
 * Default reference to assign to positive domain in iso mode
 *
 */
#define MG_PLUS    2
/**
 * Default reference to assign to negative domain in iso mode
 *
 */
#define MG_MINUS   3

/**
 * \def MMG5_ARG_start
 *
 * To begin a list of variadic arguments (mandatory first arg for all our
 * variadic functions)
 *
 * \remark we cannot use an enum because used in
 * variadic functions).
 */
#define MMG5_ARG_start  1
/**
 * \def MMG5_ARG_ppMesh
 *
 * Pointer toward a MMG5_pMesh structure (for structure allocations purposes)
 *
 * \remark we cannot use an enum because used in
 * variadic functions).
 */
#define MMG5_ARG_ppMesh 2
/**
 * \def MMG5_ARG_ppLs
 *
 * Pointer toward a MMG5_pSol structure storing a level-set (for structure
 * allocations purposes)
 *
 * \remark we cannot use an enum because used in
 * variadic functions).
 */
#define MMG5_ARG_ppLs   3
/**
 * \def MMG5_ARG_ppMet
 *
 * Pointer toward a MMG5_pSol structure storing a metric (for structure
 * allocations purposes)
 *
 * \remark we cannot use an enum because used in
 * variadic functions).
 */
#define MMG5_ARG_ppMet  4
/**
 * \def MMG5_ARG_ppDisp
 *
 * Pointer toward a MMG5_pSol structure storing a displacement (for structure
 * allocations purposes)
 *
 * \remark we cannot use an enum because used in
 * variadic functions).
 */
#define MMG5_ARG_ppDisp 5
/**
 * \def MMG5_ARG_ppSols
 *
 * Pointer toward an array of MMG5_Sol structures storing a list of solutions
 * allocations purposes)
 *
 * \remark we cannot use an enum because used in
 * variadic functions).
 */
#define MMG5_ARG_ppSols  6
/**
 * \def MMG5_ARG_pMesh
 *
 * MMG5_pMesh structure
 *
 * \remark we cannot use an enum because used in
 * variadic functions).
 */
#define MMG5_ARG_pMesh  7
/**
 * \def MMG5_ARG_pMet
 *
 * MMG5_pSol structure storing a metric field
 *
 * \remark we cannot use an enum because used in
 * variadic functions).
 */
#define MMG5_ARG_pMet   8
/**
 * \def MMG5_ARG_pDisp
 *
 * MMG5_pSol structure storing a displacement field
 *
 * \remark we cannot use an enum because used in
 * variadic functions).
 */
#define MMG5_ARG_pDisp  9
/**
 * \def MMG5_ARG_end
 *
 * To end a list of variadic argument (mandatory last argument for all our
 * variadic functions)
 *
 * \remark we cannot use an enum because used in
 * variadic functions).
 */
#define MMG5_ARG_end    10

/**
 * \def MMG5_NSOLS_MAX
 *
 * Maximal number of solutions per entity
 *
 */
#define MMG5_NSOLS_MAX   100

/**
 * \def MMG5_FILENAME_LEN_MAX
 *
 * Maximal length of filenames
 *
 */
#define  MMG5_FILENAME_LEN_MAX 255

/**
 * \def MMG5_MMAT_NOSPLIT
 *
 * Entity that must not be splitted in multimat mode
 *
 */
#define MMG5_MMAT_NoSplit  0

/**
 * \def MMG5_MMAT_Split
 *
 * Entity that must be splitted in multimat mode
 *
 */
#define MMG5_MMAT_Split  1

/**
 * \enum MMG5_type
 * \brief Type of solutions.
 */
enum MMG5_type {
  MMG5_Notype, /*!< Undefined type (unusable) */
  MMG5_Scalar, /*!< Scalar solution */
  MMG5_Vector, /*!< Vectorial solution */
  MMG5_Tensor  /*!< Tensorial solution */
};

/**
 * \enum MMG5_entities
 * \brief Identifies the types of mesh entities.
 */
enum MMG5_entities {
  MMG5_Noentity, /*!< Undefined type (unusable) */
  MMG5_Vertex, /*!< Vertex entity */
  MMG5_Edg,  /*!< Edge entity */
  MMG5_Triangle, /*!< Triangle entity */
  MMG5_Tetrahedron, /*!< Tetra entity */
};

/**
 * \enum MMG5_Format
 * \brief Type of supported file format
 */
enum MMG5_Format {
  MMG5_FMT_MeditASCII, /*!< ASCII Medit (.mesh) */
  MMG5_FMT_MeditBinary, /*!< Binary Medit (.meshb) */
  MMG5_FMT_GmshASCII, /*!< ASCII Gmsh */
  MMG5_FMT_GmshBinary, /*!< Binary Gmsh */
  MMG5_FMT_VtkPvtp, /*!< VTK pvtp */
  MMG5_FMT_VtkPvtu, /*!< VTK pvtu */
  MMG5_FMT_VtkVtu, /*!< VTK vtu */
  MMG5_FMT_VtkVtp, /*!< VTK vtp */
  MMG5_FMT_VtkVtk, /*!< VTK vtk */
  MMG5_FMT_Tetgen, /*!< Tetgen or Triangle */
  MMG5_FMT_Unknown /*!< Unrecognized */
};

/**
 * \struct MMG5_Par
 * \brief Local parameters for a specific entity and reference
 *
 * This struct can store the local values (minimal and maximal edge lengths and
 * Hausdorff distance) associated to the given reference of an element of type
 * \a elt (point, edge... ).
 *
 */
typedef struct {
  double   hmin; /*!< minimal size for edges */
  double   hmax; /*!< maximal size for edges */
  double   hausd; /*!< Hausdorff value */
  MMG5_int ref; /*!< Reference value */
  int8_t   elt; /*!< Element type */
} MMG5_Par; typedef MMG5_Par * MMG5_pPar;

/**
 * \struct MMG5_Point
 * \brief Structure to store vertices of an MMG mesh.
 * \todo What to do with n[3], try to remove s.
 */
typedef struct {
  double   c[3]; /*!< Coordinates */
  double   n[3]; /*!< Unitary normal (regular points) or unitary tangent (ridge
                  * and ref points) for mmgs and unitary tangent (if needed) for
                  * mmg3d */
#ifdef USE_POINTMAP
  MMG5_int src; /*!< Source point in input mesh */
#endif
  MMG5_int ref; /*!< Reference of point */
  MMG5_int xp; /*!< Surface point number */
  MMG5_int tmp; /*!< Index of point in the saved mesh (we don't count
                  the unused points)*/
  MMG5_int flag; /*!< Flag to know if we have already treated the point */
  MMG5_int s;
  uint16_t  tag; /*!< Contains binary flags : if \f$tag=23=16+4+2+1\f$, then
                  the point is \a MG_REF, \a MG_GEO, \a MG_REQ and \a MG_BDY */
  int8_t   tagdel; /*!< Tag for delaunay */
} MMG5_Point;
typedef MMG5_Point * MMG5_pPoint;

/**
 * \struct MMG5_xPoint
 * \brief Structure to store surface vertices of an MMG mesh.
 */
typedef struct {
  double   n1[3],n2[3]; /*!< Normals at boundary vertex;
                          n1!=n2 if the vertex belong to a ridge */
  int8_t   nnor; /* By default 0; 1 if no normal available (internal NOM point) */
} MMG5_xPoint;
typedef MMG5_xPoint * MMG5_pxPoint;

/**
 * \struct MMG5_Edge
 * \brief Structure to store edges of am MMG mesh.
 */
typedef struct {
  MMG5_int a,b; /*!< Extremities of the edge */
  MMG5_int ref; /*!< Reference of the edge */
  MMG5_int base; /*!< 2Donly: used to store the tria+ tria edge indices
                   that allow to access to the edge */
  uint16_t  tag; /*!< Binary flags */
} MMG5_Edge;
typedef MMG5_Edge * MMG5_pEdge;

/**
 * \struct MMG5_Tria
 *
 * \brief Structure to store triangles of a MMG mesh.
 *
 * \remark Numbering convention
 * \verbatim
 *  Vertices            Edges                                  *
 *  2                    .                                     *
 *  |`\                  |`\                                   *
 *  |  `\                |  `\                                 *
 *  |    `\              1    `0                               *
 *  |      `\            |      `\                             *
 *  |        `\          |        `\                           *
 *  0----------1         .--- 2 ----.
 * \endverbatim
 *
 */
typedef struct {
  double   qual;   /*Quality of the triangle*/
  MMG5_int v[3]; /*!< Vertices of the triangle */
  MMG5_int ref; /*!< Reference of the triangle */
  MMG5_int base;
  MMG5_int cc; /*!< used to store the tetra + tetra face indices
                 that allow to access to the tria 4*k + i */
  MMG5_int edg[3]; /*!< edg[i] contains the ref of the \f$i^{th}\f$ edge
                     of triangle */
  MMG5_int flag;
  uint16_t  tag[3]; /*!< tag[i] contains the tag associated to the
                     \f$i^{th}\f$ edge of triangle */
  } MMG5_Tria;
typedef MMG5_Tria * MMG5_pTria;


/**
 * \struct MMG5_Quad
 *
 * \brief Structure to store quadrangles of an MMG mesh.
 *
 * \remark Numbering convention
 * \verbatim
 *  Vertices            Edges                                   *
 *                      .                                       *
 *  3----------2         +-----3----+                           *
 *  |          |         |          |                           *
 *  |          |         1          2                           *
 *  |          |         |          |                           *
 *  |          |         |          |                           *
 *  0----------1         +----0-----+                           *
 * \endverbatim
 *
 */
typedef struct {
  MMG5_int v[4]; /*!< Vertices of the quadrangle */
  MMG5_int ref; /*!< Reference of the quadrangle */
  MMG5_int base;
  MMG5_int edg[4]; /*!< edg[i] contains the ref of the \f$i^{th}\f$ edge
                     of quadrangle */
  uint16_t  tag[4]; /*!< tag[i] contains the tag associated to the
                     \f$i^{th}\f$ edge of quadrangle */
} MMG5_Quad;
typedef MMG5_Quad * MMG5_pQuad;


/**
 * \struct MMG5_Tetra
 *
 * \brief Structure to store tetrahedra of an MMG mesh.
 *
 * \remark The numbering conventions are illustrated below. Face i lies opposite to vertex i.
 * \verbatim
 *      Vertices                     Edges                       Faces           *
 *           3                          .                           .            *
 *         ,/|`\                      ,/|`\                       ,/|`\          *
 *       ,/  |  `\                  ,/  |  `\                   ,/  |  `\        *
 *     ,/    '.   `\              ,2    '.   `5               ,/    '.   `\      *
 *   ,/       |     `\          ,/       4     `\           ,/       1     `\    *
 * ,/         |       `\      ,/         |       `\       ,/         |   0   `\  *
 * 0-----------'.--------2    .--------1--'.--------.     .------2---'.--------. *
 * `\.         |      ,/      `\.         |      ,/       `\.         |      ,/  *
 *    `\.      |    ,/           `\.      |    ,3            `\.     3|    ,/    *
 *       `\.   '. ,/                `0.   '. ,/                 `\.   '. ,/      *
 *          `\. |/                     `\. |/                      `\. |/        *
 *             `1                         `.                          `.         *
 * \endverbatim
 *
 */
typedef struct {
  double   qual; /*!< Quality of the element */
  MMG5_int v[4]; /*!< Vertices of the tetrahedron */
  MMG5_int ref; /*!< Reference of the tetrahedron */
  MMG5_int base;
  MMG5_int mark; /*!< Used for delaunay */
  MMG5_int xt; /*!< Index of the surface \ref MMG5_xTetra associated to the
                 tetrahedron (only for tetrahedra that are adjacent to
                 surfaces) */
  MMG5_int flag;
  uint16_t  tag;
} MMG5_Tetra;
typedef MMG5_Tetra * MMG5_pTetra;

/**
 * \struct MMG5_xTetra
 * \brief Structure to store additional information for the surface tetrahedra of an MMG mesh.
 */
typedef struct {
  MMG5_int ref[4]; /*!< ref[i] is the reference of the opposite triangle to the
                     \f$i^{th}\f$ vertex of the tetrahedron;*/
  MMG5_int edg[6]; /*!< edg[i] contains the reference of the
                     \f$i^{th}\f$ edge of the tetrahedron */
  uint16_t  ftag[4]; /*!< ftag[i] contains the tag associated to the
                      \f$i^{th}\f$ face of the tetrahedron */
  uint16_t  tag[6]; /*!< tag[i] contains the tag associated to the
                     \f$i^{th}\f$ edge of the tetrahedron */
  int8_t   ori; /*!< Orientation of the triangles of the tetrahedron:
                  the $\f$i^{th}\f$ bit of ori is set to 0 when the
                  \f$i^{th}\f$ face is bad orientated */
} MMG5_xTetra;
typedef MMG5_xTetra * MMG5_pxTetra;

/**
 * \struct MMG5_Prism
 *
 * \brief Structure to store prsim of a MMG mesh.
 *
 * \warning prisms are not modified
 *
 * \remark Numbering convention
 * \verbatim
 *      Vertices                   Edges                  Faces          *
 *           3                       .                      .            *
 *         ,/|`\                   ,/|`\                  ,/|`\          *
 *       ,/  |  `\                6  |   7              ,/  |  `\        *
 *     ,/    |    `\           ,/    |    `\          ,/    1    `\      *
 *    4------+------5         .------8------.        .------+------.     *
 *    |      |      |         |      |      |        |      |      |     *
 *    |      |      |         |      2      |        |      |      |     *
 *    |      |      |         |      |      |        |      |      |     *
 *    |      |      |         |      |      |        |  4   |   3  |     *
 *    |      |      |         4      |      5        |      2      |     *
 *    |      0      |         |      .      |        |      .      |     *
 *    |    ,/ `\    |         |    ,/ `\    |        |    ,/ `\    |     *
 *    |  ,/     `\  |         |  ,0     `1  |        |  ,/     `\  |     *
 *    |,/         `\|         |,/         `\|        |,/     0   `\|     *
 *    1-------------2         .------3------.        .-------------.     *
 *
 * \endverbatim
 *
 */
typedef struct {
  MMG5_int v[6]; /*!< Vertices of the prism */
  MMG5_int ref; /*!< Reference of the prism */
  MMG5_int base;
  MMG5_int flag;
  MMG5_int xpr; /*!< Index of the surface \ref MMG5_xPrism associated to
                  the prism*/
  uint8_t   tag;
} MMG5_Prism;
typedef MMG5_Prism * MMG5_pPrism;

/**
 * \struct MMG5_xPrism
 * \brief Structure to store the surface prism of a MMG mesh.
 */
typedef struct {
  MMG5_int ref[5]; /*!< face references: ref[0]={0,1,2}, ref[1]={3,4,5},
                    * ref[2]={0,3,4,1}, ref[3]={0,2,5,1} */
  MMG5_int edg[9]; /*!< edges references:
                    * edg[0]={0,1},edg[1]={0,2},edg[2]={0,3},edg[3]={1,2},
                    * edg[4]={1,4},edg[5]={2,5},edg[6]={3,4},edg[7]={3,5},
                    * edg[8]={4,5}*/
  uint16_t  ftag[5]; /*!< ftag[i] contains the tag associated to the
                      \f$i^{th}\f$ face of the prism */
  uint16_t  tag[9]; /*!< tag[i] contains the tag associated to the
                     \f$i^{th}\f$ edge of the prism */
} MMG5_xPrism;
typedef MMG5_xPrism * MMG5_pxPrism;

/**
 * \struct MMG5_Mat
 * \brief To store user-defined references in the mesh (useful in LS mode)
 */
typedef struct {
  int8_t   dospl;
  MMG5_int ref,rin,rex;
} MMG5_Mat;
typedef MMG5_Mat * MMG5_pMat;

/**
 * \struct MMG5_InvMat
 * \brief To store lookup table for references in the mesh (useful in LS mode)
 */
typedef struct {
  MMG5_int offset;
  MMG5_int size;
  int      *lookup;
} MMG5_InvMat;
typedef MMG5_InvMat * MMG5_pInvMat;

/**
 * \struct MMG5_Info
 * \brief Structure to store input parameters of the job.
 */
typedef struct {
  MMG5_pPar     par;
  double        dhd,hmin,hmax,hsiz,hgrad,hgradreq,hausd;
  double        min[3],max[3],delta,ls,lxreg,rmc;
  MMG5_int      *br; /*!< list of based references to which an implicit surface can be attached */
  MMG5_int      isoref; /*!< isovalue reference in ls mode */
  MMG5_int      nsd; /*!< index of subdomain to save (0 by default == all subdomains are saved) */
  int           mem,npar,npari;
  int           nbr,nbri; /*!< number of based references for level-set (BC to which a material can be attached) */
  int           opnbdy; /*!< floating surfaces */
  int           renum; /*!< scotch renumbering */
  int           PROctree; /*!< octree to speedup delaunay insertion */
  int           nmati,nmat; /*!< number of materials in ls multimat mode */
  int           imprim; /*!< verbosity level */
  int8_t        nreg; /*!< normal regularization */
  int8_t        xreg; /*!< vertices regularization */
  int8_t        ddebug; /*!< debug mode if 1 */
  int8_t        badkal; /*!< 1 if the mesh contains a very bad element */
  int8_t        iso; /*!< level-set discretization mode */
  int8_t        isosurf; /*!< level-set discretization mode on the surface */
  int8_t        setfem; /*!< Enforce finite element mesh (try to avoid edges
                      * connecting 2 bdy points and tet with more than 1 bdy
                      * face) */
  int8_t        fem; /*!< internal value for fem / no fem mesh output */
  int8_t        lag; /*!< lagrangian mode */
  int8_t        parTyp; /*!< Contains binary flags to say which kind of local
                          param are setted: if \f$tag = 1+2+4\f$ then the point
                          is \a MG_Vert, MG_Tria and MG_Tetra */
  int8_t        sethmin; /*!< 1 if user set hmin, 0 otherwise (needed for multiple library calls) */
  int8_t        sethmax; /*!< 1 if user set hmin, 0 otherwise (needed for multiple library calls) */
  uint8_t       ani, optim, optimLES, noinsert, noswap, nomove, nosurf, nosizreq;
  uint8_t       metRidTyp;
  char          *fparam; /*!< name of the parameter file */
  /*!< metRidTyp
   * - in 3D: 0 for a classical storage of the aniso
   * metric at ridge, 1 for the Mmg storage (modified
   * by defsiz)
   * - in 2D: used to detect if we call assignEdge function for the first time inside the library */

  MMG5_pMat     mat;
  MMG5_InvMat   invmat;
} MMG5_Info;

/**
 * \struct MMG5_hgeom
 * \brief Cell of the hash table of geometric edges.
 */
typedef struct {
  MMG5_int     a; /*!< First extremity of edge */
  MMG5_int     b;  /*!< Second extremity of edge */
  MMG5_int     ref; /*!< Reference or idx (2D) of edge */
  MMG5_int     nxt; /*!< Next element of hash table */
  uint16_t     tag; /*!< tag of edge */
} MMG5_hgeom;

/**
 * \struct MMG5_HGeom
 * \brief Hash table to store geometric edges.
 */
typedef struct {
  MMG5_hgeom  *geom;
  MMG5_int    siz,max,nxt;
} MMG5_HGeom;


/**
 * \struct MMG5_hedge
 * \brief Used to hash edges (memory economy compared to \ref MMG5_hgeom).
 */
typedef struct {
  MMG5_int   a,b,nxt;
  MMG5_int   k; /*!< k = point along edge a b or triangle index */
  MMG5_int   s;
} MMG5_hedge;

/**
 * \struct MMG5_Hash
 * \brief Identic as \ref MMG5_HGeom but use \ref MMG5_hedge to store edges
 * instead of \ref MMG5_hgeom (memory economy).
 */
typedef struct {
  MMG5_int     siz,max,nxt;
  MMG5_hedge   *item;
} MMG5_Hash;

/**
 * \struct MMG5_Mesh
 * \brief MMG mesh structure.
 * \todo try to remove nc1;
 */
typedef struct {
  size_t    memMax; /*!< Maximum memory available */
  size_t    memCur; /*!< Current memory used */
  double    gap; /*!< Gap for table reallocation */
  int       ver; /*!< Version of the mesh file */
  int       dim; /*!< Dimension of the mesh */
  int       type; /*!< Type of the mesh */
  MMG5_int  npi,nti,nai,nei,np,na,nt,ne,npmax,namax,ntmax,nemax,xpmax,xtmax;
  MMG5_int  nquad,nprism; /*!< number of quadrangles and prisms */
  int       nsols; /*!< number of solutions (metric excluded) in the solution file (lower than \a NSOLS_MAX)*/
  MMG5_int  nc1;
  MMG5_int  base; /*!< Used with \a flag to know if an entity has been
                    treated */
  MMG5_int  mark; /*!< Flag for delaunay (to know if an entity has
                    been treated) */
  MMG5_int  xp,xt,xpr; /*!< Number of surfaces points, triangles/tetrahedra and prisms */
  MMG5_int  npnil; /*!< Index of first unused point */
  MMG5_int  nenil; /*!< Index of first unused element */
  MMG5_int  nanil; /*!< Index of first unused edge (2d only)*/
  MMG5_int  *adja; /*!< Table of tetrahedron adjacency: if
                    \f$adja[4*(i-1)+1+j]=4*k+l\f$ then the \f$i^{th}\f$ and
                    \f$k^th\f$ tetrahedra are adjacent and share their
                    faces \a j and \a l (resp.) */
  MMG5_int  *adjt; /*!< Table of triangles adjacency: if
                    \f$adjt[3*(i-1)+1+j]=3*k+l\f$ then the \f$i^{th}\f$ and
                    \f$k^th\f$ triangles are adjacent and share their
                    edges \a j and \a l (resp.) */
  MMG5_int  *adjapr; /*!< Table of prisms adjacency: if
                    \f$adjapr[5*(i-1)+1+j]=5*k+l\f$ then the \f$i^{th}\f$ and
                    \f$k^th\f$ prism are adjacent and share their
                    faces \a j and \a l (resp.) */
  MMG5_int  *adjq; /*!< Table of quadrangles adjacency: if
                    \f$adjq[4*(i-1)+1+j]=4*k+l\f$ then the \f$i^{th}\f$ and
                    \f$k^th\f$ quadrilaterals are adjacent and share their
                    edges \a j and \a l (resp.) */
  int       *ipar;  /*!< Store indices of the local parameters */
  MMG5_pPoint    point; /*!< Pointer toward the \ref MMG5_Point structure */
  MMG5_pxPoint   xpoint; /*!< Pointer toward the \ref MMG5_xPoint structure */
  MMG5_pTetra    tetra; /*!< Pointer toward the \ref MMG5_Tetra structure */
  MMG5_pxTetra   xtetra; /*!< Pointer toward the \ref MMG5_xTetra structure */
  MMG5_pPrism    prism; /*!< Pointer toward the \ref MMG5_Prism structure */
  MMG5_pxPrism   xprism; /*!< Pointer toward the \ref MMG5_pxPrism structure */
  MMG5_pTria     tria; /*!< Pointer toward the \ref MMG5_Tria structure */
  MMG5_pQuad     quadra; /*!< Pointer toward the \ref MMG5_Quad structure */
  MMG5_pEdge     edge; /*!< Pointer toward the \ref MMG5_Edge structure */
  MMG5_HGeom     htab; /*!< \ref MMG5_HGeom structure */
  MMG5_Info      info; /*!< \ref MMG5_Info structure */
  char           *namein; /*!< Input mesh name */
  char           *nameout; /*!< Output mesh name */

} MMG5_Mesh;
typedef MMG5_Mesh  * MMG5_pMesh;

/**
 * \struct MMG5_sol
 * \brief MMG Solution structure (for solution or metric).
 *
 */
typedef struct {
  int       ver; /* Version of the solution file */
  int       dim; /* Dimension of the solution file*/
  MMG5_int  np; /* Number of points of the solution */
  MMG5_int  npmax; /* Maximum number of points */
  MMG5_int  npi; /* Temporary number of points (internal use only) */
  int       size; /* Number of solutions per entity */
  int       type; /* Type of the solution (scalar, vectorial or tensorial) */
  int       entities; /* Type of the solution (scalar, vectorial of tensorial) */
  double    *m; /*!< Solution values */
  double    umin,umax; /*!<Min/max values for the solution */
  char      *namein; /*!< Input solution file name */
  char      *nameout; /*!< Output solution file name */
} MMG5_Sol;
typedef MMG5_Sol * MMG5_pSol;

#endif
