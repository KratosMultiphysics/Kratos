/* =============================================================================
**  This file is part of the mmg software package for the tetrahedral
**  mesh modification.
**  Copyright (c) Bx INP/Inria/UBordeaux/UPMC, 2004- .
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
 * Types
 */
#include <stdint.h>

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
 * Implicite domain ref in iso mode
 *
 */
#define MG_ISO    10

#include <stdarg.h>

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
 * \def MMG5_ARG_pMesh
 *
 * MMG5_pMesh structure
 *
 * \remark we cannot use an enum because used in
 * variadic functions).
 */
#define MMG5_ARG_pMesh  6
/**
 * \def MMG5_ARG_pMet
 *
 * MMG5_pSol structure storing a metric field
 *
 * \remark we cannot use an enum because used in
 * variadic functions).
 */
#define MMG5_ARG_pMet   7
/**
 * \def MMG5_ARG_pDisp
 *
 * MMG5_pSol structure storing a displacement field
 *
 * \remark we cannot use an enum because used in
 * variadic functions).
 */
#define MMG5_ARG_pDisp  8
/**
 * \def MMG5_ARG_end
 *
 * To end a list of variadic argument (mandatory last argument for all our
 * variadic functions)
 *
 * \remark we cannot use an enum because used in
 * variadic functions).
 */
#define MMG5_ARG_end    9

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
 * \brief Type of mesh entities.
 */
enum MMG5_entities {
  MMG5_Noentity, /*!< Undefined type (unusable) */
  MMG5_Vertex, /*!< Vertex entity */
  MMG5_Triangle, /*!< Triangle entity */
  MMG5_Tetrahedron, /*!< Tetra entity */
};

/**
 * \struct MMG5_Par
 * number) associated to a specific reference.
 *
 * Store the local values (minimal and maximal sizes and Hausdorff number)
 * associated to the given reference of an element of type \a elt (point,
 * edge... ).
 *
 */
typedef struct {
  double   hmin; /*!< minimal size for edges */
  double   hmax; /*!< maximal size for edges */
  double   hausd; /*!< Hausdorff value */
  int      ref; /*!< Reference value */
  char     elt; /*!< Element type */
} MMG5_Par; typedef MMG5_Par * MMG5_pPar;

/**
 * \struct MMG5_Point
 * \brief Structure to store points of a MMG mesh.
 * \todo What to do with n[3], try to remove s.
 */
typedef struct {
  double   c[3]; /*!< Coordinates of point */
  double   n[3]; /*!< Normal or Tangent for mmgs and Tangent (if needed) for mmg3d */
  int      ref; /*!< Reference of point */
  int      xp; /*!< Surface point number */
  int      tmp; /*!< Index of point in the saved mesh (we don't count
                  the unused points)*/
  int      flag; /*!< Flag to know if we have already treated the point */
  int      s;
  int16_t  tag; /*!< Contains binary flags : if \f$tag=23=16+4+2+1\f$, then
                  the point is \a MG_REF, \a MG_GEO, \a MG_REQ and \a MG_BDY */
  char     tagdel; /*!< Tag for delaunay */
} MMG5_Point;
typedef MMG5_Point * MMG5_pPoint;

/**
 * \struct MMG5_xPoint
 * \brief Structure to store surface points of a MMG mesh.
 */
typedef struct {
  double   n1[3],n2[3]; /*!< Normals at boundary vertex;
                          n1!=n2 if the vertex belong to a ridge */
} MMG5_xPoint;
typedef MMG5_xPoint * MMG5_pxPoint;

/**
 * \struct MMG5_Edge
 * \brief Structure to store edges of a MMG mesh.
 */
typedef struct {
  int      a,b; /*!< Extremities of the edge */
  int      ref; /*!< Reference of the edge */
  int      base; /*!< 2Donly: used to store the tria+ tria edge indices
                   that allow to access to the edge */
  int16_t  tag; /*!< Binary flags */
} MMG5_Edge;
typedef MMG5_Edge * MMG5_pEdge;

/**
 * \struct MMG5_Tria
 *
 * Structure to store triangles of a MMG mesh.
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
  int      v[3]; /*!< Vertices of the triangle */
  int      ref; /*!< Reference of the triangle */
  int      base;
  int      cc; /*!< used to store the tetra + tetra face indices
                 that allow to access to the tria */
  int      edg[3]; /*!< edg[i] contains the ref of the \f$i^{th}\f$ edge
                     of triangle */
  int      flag;
  int16_t  tag[3]; /*!< tag[i] contains the tag associated to the
                     \f$i^{th}\f$ edge of triangle */
} MMG5_Tria;
typedef MMG5_Tria * MMG5_pTria;


/**
 * \struct MMG5_Quad
 *
 * Structure to store quadrangles of a MMG mesh.
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
  int      v[4]; /*!< Vertices of the quadrangle */
  int      ref; /*!< Reference of the quadrangle */
  int      base;
  int      edg[4]; /*!< edg[i] contains the ref of the \f$i^{th}\f$ edge
                     of quadrangle */
  int16_t  tag[4]; /*!< tag[i] contains the tag associated to the
                     \f$i^{th}\f$ edge of quadrangle */
} MMG5_Quad;
typedef MMG5_Quad * MMG5_pQuad;


/**
 * \struct MMG5_Tetra
 *
 * Structure to store tetrahedra of a MMG mesh.
 *
 * \remark Numbering convention
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
  int      v[4]; /*!< Vertices of the tetrahedron */
  int      ref; /*!< Reference of the tetrahedron */
  int      base;
  int      mark; /*!< Used for delaunay */
  int      xt; /*!< Index of the surface \ref MMG5_xTetra associated to
                 the tetrahedron*/
  int      flag;
  int16_t  tag;
} MMG5_Tetra;
typedef MMG5_Tetra * MMG5_pTetra;

/**
 * \struct MMG5_xTetra
 * \brief Structure to store the surface tetrahedra of a MMG mesh.
 */
typedef struct {
  int      ref[4]; /*!< ref[i] is the reference of the opposite triangle to the
                     \f$i^{th}\f$ vertex of the tetrahedron;*/
  int      edg[6]; /*!< edg[i] contains the reference of the
                     \f$i^{th}\f$ edge of the tetrahedron */
  int16_t  ftag[4]; /*!< ftag[i] contains the tag associated to the
                      \f$i^{th}\f$ face of the tetrahedron */
  int16_t  tag[6]; /*!< tag[i] contains the tag associated to the
                     \f$i^{th}\f$ edge of the tetrahedron */
  char     ori; /*!< Orientation of the triangles of the tetrahedron:
                  the $\f$i^{th}\f$ bit of ori is set to 0 when the
                  \f$i^{th}\f$ face is bad orientated */
} MMG5_xTetra;
typedef MMG5_xTetra * MMG5_pxTetra;

/**
 * \struct MMG5_Prism
 *
 * Structure to store prsim of a MMG mesh.
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
  int      v[6]; /*!< Vertices of the prism */
  int      ref; /*!< Reference of the prism */
  int      base;
  int      flag;
  int      xpr; /*!< Index of the surface \ref MMG5_xPrism associated to
                  the prism*/
  char     tag;
} MMG5_Prism;
typedef MMG5_Prism * MMG5_pPrism;

/**
 * \struct MMG5_xPrism
 * \brief Structure to store the surface prism of a MMG mesh.
 */
typedef struct {
  int      ref[5]; /*!< face references: ref[0]={0,1,2}, ref[1]={3,4,5},
                    * ref[2]={0,3,4,1}, ref[3]={0,2,5,1} */
  int      edg[9]; /*!< edges references:
                    * edg[0]={0,1},edg[1]={0,2},edg[2]={0,3},edg[3]={1,2},
                    * edg[4]={1,4},edg[5]={2,5},edg[6]={3,4},edg[7]={3,5},
                    * edg[8]={4,5}*/
  int16_t ftag[5]; /*!< ftag[i] contains the tag associated to the
                      \f$i^{th}\f$ face of the prism */
  int16_t tag[9]; /*!< tag[i] contains the tag associated to the
                     \f$i^{th}\f$ edge of the prism */
} MMG5_xPrism;
typedef MMG5_xPrism * MMG5_pxPrism;
/**
 * \struct MMG5_Info
 * \brief Store input parameters of the run.
 */
typedef struct {
  MMG5_pPar     par;
  double        dhd,hmin,hmax,hgrad,hausd,min[3],max[3],delta,ls;
  int           mem,npar,npari;
  int           renum;
  int           octree;
  char          nreg;
  char          imprim,ddebug,badkal,iso,fem,lag;
  char          parTyp; /*!< Contains binary flags to say which kind of local
                          param are setted: if \f$tag = 1+2+4\f$ then the point
                          is \a MG_Vert, MG_Tria and MG_Tetra */
  unsigned char optim, optimLES, noinsert, noswap, nomove, nosurf;
} MMG5_Info;

/**
 * \struct MMG5_hgeom
 * \brief To store geometric edges.
 */
typedef struct {
  int     a; /*!< First extremity of edge */
  int     b;  /*!< Second extremity of edge */
  int     ref; /*!< Reference of edge */
  int     nxt; /*!< Next element of hash table */
  int16_t tag; /*!< tag of edge */
} MMG5_hgeom;

typedef struct {
  MMG5_hgeom  *geom;
  int         siz,max,nxt;
} MMG5_HGeom;

/**
 * \struct MMG5_Mesh
 * \brief MMG mesh structure.
 * \todo try to remove nc1;
 */
typedef struct {
  long long memMax; /*!< Maximum memory available */
  long long memCur; /*!< Current memory used */
  double    gap; /*!< Gap for table reallocation */
  int       ver; /*!< Version of the mesh file */
  int       dim; /*!< Dimension of the mesh */
  int       type; /*!< Type of the mesh */
  int       npi,nti,nai,nei,np,na,nt,ne,npmax,namax,ntmax,nemax,xpmax,xtmax;
  int       nquad,nprism; /* number of quadrangles and prisms */
  int       nc1;

  int       base; /*!< Used with \a flag to know if an entity has been
                    treated */
  int       mark; /*!< Flag for delaunay (to know if an entity has
                    been treated) */
  int       xp,xt,xpr; /*!< Number of surfaces points, triangles/tetrahedra and prisms */
  int       npnil; /*!< Index of first unused point */
  int       nenil; /*!< Index of first unused element */
  int       nanil; /*!< Index of first unused edge (2d only)*/
  int      *adja; /*!< Table of tetrahedron adjacency: if
                    \f$adja[4*i+1+j]=4*k+l\f$ then the \f$i^{th}\f$ and
                    \f$k^th\f$ tetrahedra are adjacent and share their
                    faces \a j and \a l (resp.) */
  int      *adjt; /*!< Table of triangles adjacency: if
                    \f$adjt[3*i+1+j]=3*k+l\f$ then the \f$i^{th}\f$ and
                    \f$k^th\f$ triangles are adjacent and share their
                    edges \a j and \a l (resp.) */
  int      *adjapr; /*!< Table of prisms adjacency: if
                    \f$adjapr[5*i+1+j]=5*k+l\f$ then the \f$i^{th}\f$ and
                    \f$k^th\f$ prism are adjacent and share their
                    faces \a j and \a l (resp.) */
  MMG5_pPoint    point; /*!< Pointer toward the \ref MMG5_Point structure */
  MMG5_pxPoint   xpoint; /*!< Pointer toward the \ref MMG5_xPoint structure */
  MMG5_pTetra    tetra; /*!< Pointer toward the \ref MMG5_Tetra structure */
  MMG5_pxTetra   xtetra; /*!< Pointer toward the \ref MMG5_xTetra structure */
  MMG5_pPrism    prism; /*!< Pointer toward the \ref MMG5_Prism structure */
  MMG5_pxPrism   xprism; /*!< Pointer toward the \ref MMG5_pxPrism structure */
  MMG5_pTria     tria; /*!< Pointer toward the \ref MMG5_Tria structure */
  MMG5_pQuad     quad; /*!< Pointer toward the \ref MMG5_Quad structure */
  MMG5_pEdge     edge; /*!< Pointer toward the \ref MMG5_Edge structure */
  MMG5_HGeom     htab; /*!< \ref MMG5_HGeom structure */
  MMG5_Info      info; /*!< \ref MMG5_Info structure */
  char     *namein; /*!< Input mesh name */
  char     *nameout; /*!< Output mesh name */

} MMG5_Mesh;
typedef MMG5_Mesh  * MMG5_pMesh;

/**
 * \struct MMG5_sol
 * \brief MMG Solution structure (for solution or metric).
 */
typedef struct {
  int       ver; /* Version of the solution file */
  int       dim; /* Dimension of the solution file*/
  int       np; /* Number of points of the solution */
  int       npmax; /* Maximum number of points */
  int       npi; /* Temporary number of points (internal use only) */
  int       size; /* Number of solutions per entity */
  int       type; /* Type of the solution (scalar, vectorial of tensorial) */
  double   *m; /*!< Solution values */
  char     *namein; /*!< Input solution file name */
  char     *nameout; /*!< Output solution file name */
} MMG5_Sol;
typedef MMG5_Sol * MMG5_pSol;

#endif
