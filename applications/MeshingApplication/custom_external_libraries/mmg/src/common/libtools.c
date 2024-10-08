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
 * \brief Functions needed by libraries API
 * \author Charles Dapogny (UPMC)
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \date 04 2015
 * \copyright GNU Lesser General Public License.
 **/

#include "mmgcommon_private.h"

/**
 * \param mesh pointer toward the mesh
 * \param dim string dontaining the dimension (3D,2D or S)
 *
 * Print MMG release and date
 */
void MMG5_version(MMG5_pMesh mesh,char *dim) {

  if ( mesh->info.imprim >= 0 ) {
#ifndef MMG_COMPARABLE_OUTPUT
    fprintf(stdout,"\n  %s\n   MODULE MMG%s: %s (%s)\n  %s\n",
            MG_STR,dim,MMG_VERSION_RELEASE,MMG_RELEASE_DATE,MG_STR);
#else
    fprintf(stdout,"\n  %s\n   MODULE MMG%s\n  %s\n",
            MG_STR,dim,MG_STR);
#endif

#if !defined _WIN32 && !defined MMG_COMPARABLE_OUTPUT
    fprintf(stdout,"     git branch: %s\n",MMG_GIT_BRANCH);
    fprintf(stdout,"     git commit: %s\n",MMG_GIT_COMMIT);
    fprintf(stdout,"     git date:   %s\n\n",MMG_GIT_DATE);
#endif
  }

}

/**
 * \param mesh pointer toward the mesh structure.
 * \return 0 if fail, 1 if success.
 *
 * Print the default parameters values.
 *
 */
void MMG5_mmgDefaultValues(MMG5_pMesh mesh) {

  fprintf(stdout,"\nDefault parameters values:\n");

  fprintf(stdout,"\n** Generic options :\n");
  fprintf(stdout,"verbosity                 (-v)      : %d\n",
          mesh->info.imprim);
  fprintf(stdout,"maximal memory size       (-m)      : %zu MB\n",
          mesh->memMax/MMG5_MILLION);


  fprintf(stdout,"\n**  Parameters\n");
  fprintf(stdout,"angle detection           (-ar)     : %lf\n",
          180/M_PI*acos(mesh->info.dhd) );
  fprintf(stdout,"minimal mesh size         (-hmin)   : %lf\n"
          "If not yet computed: 0.001 of "
          "the mesh bounding box if no metric is provided, 0.1 times the "
          "minimum of the metric sizes otherwise.\n",mesh->info.hmin);
  fprintf(stdout,"maximal mesh size         (-hmax)   : %lf\n"
          " If not yet computed: size of "
          "the mesh bounding box without metric, 10 times the maximum of the "
          "metric sizes otherwise.\n",mesh->info.hmax);
  fprintf(stdout,"Hausdorff distance        (-hausd)  : %lf\n",
          mesh->info.hausd);

  fprintf(stdout,"gradation control         (-hgrad)  : %lf\n",
          (mesh->info.hgrad < 0) ? mesh->info.hgrad : exp(mesh->info.hgrad) );

  fprintf(stdout,"gradation control for required entities (-hgradreq)  : %lf\n",
          (mesh->info.hgradreq < 0) ? mesh->info.hgradreq : exp(mesh->info.hgradreq) );
}

int MMG5_Set_multiMat(MMG5_pMesh mesh,MMG5_pSol sol,MMG5_int ref,
                      int split,MMG5_int rin,MMG5_int rex){
  MMG5_pMat mat;
  int k;

  (void)sol;

  if ( !mesh->info.nmat ) {
    fprintf(stderr,"\n  ## Error: %s: You must set the number of material",__func__);
    fprintf(stderr," with the MMG2D_Set_iparameters function before setting");
    fprintf(stderr," values in multi material structure. \n");
    return 0;
  }
  if ( mesh->info.nmati >= mesh->info.nmat ) {
    fprintf(stderr,"\n  ## Error: %s: unable to set a new material.\n",
            __func__);
    fprintf(stderr,"    max number of materials: %d\n",mesh->info.nmat);
    return 0;
  }
  if ( ref < 0 ) {
    fprintf(stderr,"\n  ## Error: %s: negative references are not allowed.\n",
            __func__);
    return 0;
  }

  for (k=0; k<mesh->info.nmati; k++) {
    mat = &mesh->info.mat[k];

    if ( mat->ref == ref ) {
      mat->dospl = split;
      if ( split ) {
        mat->rin   = rin;
        mat->rex   = rex;
      }
      else {
        mat->rin = mat->ref;
        mat->rex = mat->ref;
      }
      if ( (mesh->info.imprim > 5) || mesh->info.ddebug ) {
        fprintf(stderr,"\n  ## Warning: %s: new materials (interior, exterior)",
                __func__);
        fprintf(stderr," for material of ref %" MMG5_PRId "\n",ref);
      }
      return 1;
    }
  }

  if ( ( split != MMG5_MMAT_Split ) && ( split != MMG5_MMAT_NoSplit ) ) {
    fprintf(stderr,"\n ## Error: %s: unexpected value for the 'split' argument."
            " You must use the MMG5_MMAT_Split or MMG5_MMAT_NpSplit keywords \n",
            __func__);
    return 0;
  }

  mesh->info.mat[mesh->info.nmati].ref   = ref;
  mesh->info.mat[mesh->info.nmati].dospl = split;
  mesh->info.mat[mesh->info.nmati].rin   = rin;
  mesh->info.mat[mesh->info.nmati].rex   = rex;

  mesh->info.nmati++;

  /* Invert the table if all materials have been set */
  if( mesh->info.nmati == mesh->info.nmat )
    if( !MMG5_MultiMat_init(mesh) ) {
      fprintf(stderr,"\n ## Error: %s: unable to create lookup table for multiple materials.\n",
              __func__);
      return 0;
    }

  return 1;
}

int MMG5_Set_lsBaseReference(MMG5_pMesh mesh,MMG5_pSol sol,MMG5_int br) {

  (void)sol;

  if ( !mesh->info.nbr ) {
    fprintf(stderr,"\n  ## Error: %s: You must set the number of"
            " level-set based references",__func__);
    fprintf(stderr," with the MMG2D_Set_iparameters function before setting");
    fprintf(stderr," based references values. \n");
    return 0;
  }
  if ( mesh->info.nbri >= mesh->info.nbr ) {
    fprintf(stderr,"\n  ## Error: %s: unable to set a new level-set"
            " based reference.\n",__func__);
    fprintf(stderr,"    max number of level-set based references: %d\n",mesh->info.nbr);
    return 0;
  }
  if ( br < 0 ) {
    fprintf(stderr,"\n  ## Error: %s: negative references are not allowed.\n",
            __func__);
    return 0;
  }

  mesh->info.br[mesh->info.nbri] = br;
  mesh->info.nbri++;

  return 1;
}


/**
 * \param *prog pointer toward the program name.
 *
 * Print help for common options of the 3 codes (first section).
 *
 */
void MMG5_mmgUsage(char *prog) {
  fprintf(stdout,"\nUsage: %s [-v [n]] [opts..] filein [fileout]\n",prog);

  fprintf(stdout,"\n** Generic options\n");
  fprintf(stdout,"-h        Print this message\n");
  fprintf(stdout,"-v [n]    Tune level of verbosity, [-1..10]\n");
  fprintf(stdout,"-m [n]    Set maximal memory size to n Mbytes\n");
  fprintf(stdout,"-d        Turn on debug mode\n");
  fprintf(stdout,"-val      Print the default parameters values\n");
  fprintf(stdout,"-default  Save a local parameters file for default parameters"
          " values\n");

  fprintf(stdout,"\n**  File specifications\n");
  fprintf(stdout,"-in  file  input triangulation\n");
  fprintf(stdout,"-out file  output triangulation\n");
  fprintf(stdout,"-sol file  load solution or metric file\n");
  fprintf(stdout,"-met file  load metric file\n");

  fprintf(stdout,"\n**  Mode specifications (mesh adaptation by default)\n");
  fprintf(stdout,"-ls     val create mesh of isovalue val (0 if no argument provided)\n");
  fprintf(stdout,"-lssurf val split mesh boundaries on isovalue val (0 if no argument provided)\n");

}

/**
 *
 * Print help for common parameters options of the 3 codes (first section).
 *
 */
void MMG5_paramUsage1(void) {
  fprintf(stdout,"\n**  Parameters\n");
  fprintf(stdout,"-A           enable anisotropy (without metric file).\n");
  fprintf(stdout,"-ar     val  angle detection\n");
  fprintf(stdout,"-nr          no angle detection\n");
  fprintf(stdout,"-hausd  val  control Hausdorff distance\n");
  fprintf(stdout,"-hgrad  val  control gradation\n");
  fprintf(stdout,"-hmax   val  maximal mesh size\n");
  fprintf(stdout,"-hmin   val  minimal mesh size\n");
  fprintf(stdout,"-hsiz   val  constant mesh size\n");
  fprintf(stdout,"-rmc   [val] enable the removal of componants whose volume fraction is less than\n"
          "             val (1e-5 if not given) of the mesh volume (ls mode).\n");
}

/**
 *
 * Print help for common options of the 3 codes (second section).
 *
 */
void MMG5_paramUsage2(void) {

  fprintf(stdout,"-noinsert    no point insertion/deletion \n");
  fprintf(stdout,"-nomove      no point relocation\n");
  fprintf(stdout,"-noswap      no edge or face flipping\n");
  fprintf(stdout,"-nreg        normal regul.\n");
  fprintf(stdout,"-xreg        vertex regul.\n");
  fprintf(stdout,"-nsd    val  save the subdomain number val (0==all subdomain)\n");
  fprintf(stdout,"-optim       mesh optimization\n");

}

/**
 *
 * Print help for lagrangian motion option.
 *
 */
void MMG5_lagUsage(void) {

#ifdef USE_ELAS
  fprintf(stdout,"-lag [n]     lagrangian mesh displacement according to mode [0/1/2]\n");
  fprintf(stdout,"               0: displacement\n");
  fprintf(stdout,"               1: displacement + remeshing (swap and move)\n");
  fprintf(stdout,"               2: displacement + remeshing (split, collapse,"
          " swap and move)\n");
#endif
}

/**
 *
 * Print help for common options between 2D and 3D.
 *
 */
void MMG5_2d3dUsage(void) {

  fprintf(stdout,"-opnbdy      preserve input triangles at the interface of"
          " two domains of the same reference.\n");
}

/**
 *
 * Print help for advanced users of mmg.
 *
 */
void MMG5_advancedUsage(void) {

  fprintf(stdout,"\n**  Parameters for advanced users\n");
  fprintf(stdout,"-nosizreq       disable setting of required edge sizes over required vertices.\n");
  fprintf(stdout,"-hgradreq  val  control gradation from required entities toward others\n");

}

/**
 * \param mesh pointer toward mesh
 * \param pa pointer toward edge
 *
 * Clean tags linked to iso surface discretization (MG_CRN, MG_ISO) along edge.
 *
 */
static inline
void MMG5_Clean_isoTags(MMG5_pMesh mesh,MMG5_pEdge pa) {
  /* Remove MG_REQ and MG_CRN tags on ISO edges extremities */
  if ( MG_REQ & mesh->point[pa->a].tag ) {
    mesh->point[pa->a].tag &= ~MG_REQ;
  }
  if ( MG_REQ & mesh->point[pa->b].tag ) {
    mesh->point[pa->b].tag &= ~MG_REQ;
  }
  if ( MG_CRN & mesh->point[pa->a].tag ) {
    mesh->point[pa->a].tag &= ~MG_CRN;
  }
  if ( MG_CRN & mesh->point[pa->b].tag ) {
    mesh->point[pa->b].tag &= ~MG_CRN;
  }
}

/**
  * \param mesh pointer toward mesh
  * \param pa pointer toward mesh edge
  *
  * \return 1 if edge should be removed when cleaning old iso surface, 0 otherwise.
  *
  * Check if edge \a pa has to be removed from list of edges when cleaning old
  * isosurface.
  *
  */
static inline
int8_t MMG5_should_edge_be_removed(MMG5_pMesh mesh,MMG5_pEdge pa){
  int8_t to_remove;

  if ( !pa->a ) {
    to_remove = 1;
  }
  else {
    int8_t not_ridge = !(pa->tag & MG_GEO);
    to_remove = (MMG5_abs(pa->ref) == mesh->info.isoref) && not_ridge;
  }

  return to_remove;
}

/**
 * \param mesh pointer toward mesh
 * \param return 1 if successful, 0 if fail
 *
 * Clean edges belonging to isosurf, except for ridges.
 */
int MMG5_Clean_isoEdges(MMG5_pMesh mesh) {
  MMG5_int   k,nref;

  nref = 0;

  /** Deletion of edges that belong to isosurf */
  if ( mesh->edge ) {
    MMG5_int na = mesh->na;

    k  = 1;
    do {

      MMG5_pEdge pa = &mesh->edge[k];
      if ( !pa->a ) {
        continue;
      }

      if ( MMG5_abs(pa->ref) == mesh->info.isoref ) {
        /* Current tria will be suppressed */
        /* Remove MG_REQ and MG_CRN tags on ISO edges extremities */
        MMG5_Clean_isoTags(mesh,pa);

        /* Do not delete ridge */
        if ( !(pa->tag & MG_GEO) ) {
          /* search last non isosurf tria to fill empty position */
          MMG5_pEdge pa1 = &mesh->edge[mesh->na];
          assert( pa1 );

          int8_t to_remove = MMG5_should_edge_be_removed(mesh,pa1);

          while ( to_remove && (k < mesh->na) ) {
            if ( pa1->a ) {
              /* Remove MG_REQ and MG_CRN tags on ISO edges extremities */
              MMG5_Clean_isoTags(mesh,pa1);
            }
            --mesh->na;
            pa1 = &mesh->edge[mesh->na];
            to_remove = MMG5_should_edge_be_removed(mesh,pa1);
          }
          if ( pa != pa1 ) {
            /* We don't find any edge to keep after index k */
            memcpy(pa,pa1,sizeof(MMG5_Edge));
            --mesh->na;
          }
        }
      }

      /* Initially negative refs were used to mark isosurface: keep following
       * piece of code for retrocompatibility */
      if ( pa->ref < 0 ) {
        pa->ref = -pa->ref;
        ++nref;
      }
    }
    while ( ++k < mesh->na );

    /* At the end of the loop, either k==mesh->na, either k==mesh->na+1 (because
     * edg at idx mesh->na was iso or unused and element mesh->na+1 has been
     * copied into k) */
    assert ( (k==mesh->na) || (k==mesh->na+1) );

    /* Check if last edge is iso */
    MMG5_pEdge pa = &mesh->edge[mesh->na];
    if ( (!pa->a) || (MMG5_abs(pa->a) == mesh->info.isoref) ) {
      --mesh->na;
    }

    if ( mesh->info.imprim > 4 ) {
      fprintf(stdout,"     Deleted iso edges: %" MMG5_PRId "\n",na-mesh->na);
    }

    if( !mesh->na ) {
      MMG5_DEL_MEM(mesh,mesh->edge);
    }
    else if ( mesh->na < na ) {
      MMG5_ADD_MEM(mesh,(mesh->na-na)*sizeof(MMG5_Edge),"edges",
                   fprintf(stderr,"  Exit program.\n");
                   return 0);
      MMG5_SAFE_RECALLOC(mesh->edge,na+1,(mesh->na+1),MMG5_Edge,
                         "edges",return 0);
    }
  }

  return 1;
}
