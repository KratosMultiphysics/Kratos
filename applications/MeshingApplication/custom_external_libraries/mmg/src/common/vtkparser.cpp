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
 * \brief I/O at VTK format
 * \author Charles Dapogny (UPMC)
 * \author Pascal Frey (UPMC)
 * \author Algiane Froehly (Inria)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */
#include <vtkSmartPointer.h>
#include <vtkXMLReader.h>
#include <vtkXMLWriter.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkDataSetReader.h>
#include <vtkDataSet.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>
#include <vtkStructuredGrid.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkFieldData.h>
#include <vtkCellTypes.h>
#include <vtkDataArray.h>

#include "vtkparser.hpp"

/// @tparam TReader one of the VTK reader class.
/// @param filename name of the input file.
///
/// @return pointer toward a vtkDataSet object that contains the mesh and the
/// associated data.
///
/// Open and read a vtk xml file (.vtp, .vtu or .vtk). Store the data in a
/// vtkDataSet object.
///
///
template<class TReader>
vtkDataSet *MMG5_load_vtkXMLFile(const char*fileName)
{
  vtkSmartPointer<TReader> reader =
    vtkSmartPointer<TReader>::New();
  reader->SetFileName(fileName);

  reader->Update();
  if ( !reader->GetOutput() ) {
    throw "Unable to open file.";
  }

  reader->GetOutput()->Register(reader);

  return vtkDataSet::SafeDownCast(reader->GetOutput());
}

/// @param dataset pointer toward a vtkDataSet structure
/// @param mesh pointer toward a MMG5 mesh
/// @param ptMeditRef index of point data field that contains references
/// (field named medit:ref), -1 if no references
/// @param eltMeditRef index of a cell data field that contains references
/// (field named medit:ref), -1 if no references
/// @param nsols number of point data (except the medit:ref ones)
/// @param metricData 1 if file contains a metric data highlighted by the :metric name
///
/// @return 1 if success, -1 otherwise
///
/// Count the number of entities of eache type (points, triangles...) in the
/// mesh as well as the number of node data (solutions). If a data field name
/// contains the "medit:ref" string, it will be used as entity reference instead
/// of a solution.  If node references are detected, \a ptMeditRef is
/// incermented. If cell references are detected, \a eltMeditRef is
/// incremented.
///
/// @remark Note that Mmg supports only 1 reference by point and 1 reference per
/// cell.
///
///
static
int MMG5_count_vtkEntities ( vtkDataSet *dataset, MMG5_pMesh mesh,
                             int8_t *ptMeditRef, int8_t *eltMeditRef,
                             int *nsols, int8_t *metricData ) {

  static int8_t mmgWarn0 = 0;
  static int8_t mmgWarn1 = 0;

  // Count the number of entities in the mesh
  assert ( !mesh->np );
  assert ( !mesh->nt );
  assert ( !mesh->na );
  assert ( !mesh->nquad );
  assert ( !mesh->nprism );

  mesh->np           = dataset->GetNumberOfPoints();
  vtkIdType numCells = dataset->GetNumberOfCells();

  assert ( !mesh->npi );
  for ( vtkIdType i = 0; i < numCells; i++ ) {
    int typ = dataset->GetCellType(i);
    switch ( typ ) {
    case ( VTK_VERTEX ):
      ++mesh->npi;
      break;
    case ( VTK_POLY_LINE ):
      printf( "polylin %lld \n", dataset->GetCell(i)->GetNumberOfPoints() - 1 );

      printf( "%lld %lld %lld\n", dataset->GetCell(i)->GetPointId(0)+1,
              dataset->GetCell(i)->GetPointId(1)+1,dataset->GetCell(i)->GetPointId(2)+1);
      mesh->na += dataset->GetCell(i)->GetNumberOfPoints() - 1;
      break;
    case ( VTK_LINE ):
      ++mesh->na;
      break;
    case ( VTK_TRIANGLE ):
      ++mesh->nt;
      break;
    case ( VTK_QUAD):
      ++mesh->nquad;
      break;
    case ( VTK_TETRA):
      ++mesh->ne;
      break;
    case ( VTK_WEDGE):
      ++mesh->nprism;
      break;
    default:
      if ( !mmgWarn1 ) {
        printf( "  ## Warning:%s: unexpected element type (%d).",__func__,typ);
        mmgWarn1 = 1;
      }
    }
  }
  assert ( mesh->npi == mesh->np || !mesh->npi );
  mesh->npi = 0;

  // Count the number of point data and detect point references
  auto *pd = dataset->GetPointData();
  int npointData   = 0;
  int npointRef    = 0;
  int nmetricField = 0;

  *ptMeditRef = -1;
  if ( pd ) {
    npointData = pd->GetNumberOfArrays();
    for (int k = 0; k < npointData; k++) {
      if  ( strstr(pd->GetArrayName(k),":metric" ) ) {
        ++nmetricField;
      }
      else if ( strstr(pd->GetArrayName(k),"medit:ref" ) ) {
        (*ptMeditRef) = k;
        ++npointRef;
      }
    }

    if ( npointRef > 1 ) {
      printf( "  ## Warning:%s: %d reference fields detected (labelled 'medit:ref')."
              " Only the last is used, others are ignored.", __func__, npointRef);
    }
    if ( nmetricField > 1 ) {
      printf("  ## Error:%s: %d metric fields detected (labelled with a string"
             " containing the 'metric' keyword).\n"
             " Exit Program.\n",__func__,nmetricField);
      return -1;
    }
  }

  // Count the number of cell data and detect cell references
  auto *cd = dataset->GetCellData();
  int ncellData   = 0;
  int ncellRef    = 0;

  *eltMeditRef = -1;
  if ( cd ) {
    ncellData = cd->GetNumberOfArrays();
    for (int k = 0; k < ncellData; k++) {
      if (  strstr(cd->GetArrayName(k),"medit:ref" ) ) {
        *eltMeditRef = k;
        ++ncellRef;
      }
    }
  }

  // Count the number of field data
  auto *fd = dataset->GetFieldData();
  if ( fd->GetNumberOfArrays() ) {
    printf( "  ## Warning:%s: VTK field data not used in Mmg."
            " Ignored.\n",__func__ );
  }

  *nsols      = npointData + ncellData - npointRef - ncellRef;
  *metricData = ( nmetricField > 0 );

  return 1;
}

/// @param mesh pointer toward a MMG5 mesh
/// @param filename name of the input file.
/// @param dataset vtkdataset structure
/// @param ptMeditRef index of point data field that contains references
/// (field named medit:ref), -1 if no references
/// @param eltMeditRef index of a cell data field that contains references
/// (field named medit:ref), -1 if no references
/// @param nsols number of point data (except the medit:ref ones)
/// @param metricData 1 if file contains a metric data highlighted by the :metric name
///
/// @return 1 if success, 0 if fail to open/load the file, -1 otherwise;
///
/// I/O at Vtp VTK file format.
///
int MMG5_loadVtpMesh_part1(MMG5_pMesh mesh,const char *filename,vtkDataSet **dataset,
                           int8_t *ptMeditRef,int8_t *eltMeditRef,int *nsols,
                           int8_t *metricData) {

  (*nsols) = 0;
  (*ptMeditRef) = (*eltMeditRef) = -1;

  // Read all the data from the file
  try {
    (*dataset) = MMG5_load_vtkXMLFile<vtkXMLPolyDataReader> ( filename );
  }
  catch( ... ) {
    return 0;
  }

  // count the number of entities of each type
  int ier = MMG5_count_vtkEntities ( (*dataset),mesh,ptMeditRef,eltMeditRef,
                                     nsols,metricData );

  if ( ier != 1 ) {
    return -1;
  }
  return 1;
}

/// @param mesh pointer toward a MMG5 mesh
/// @param filename name of the input file.
/// @param dataset vtkdataset structure
/// @param ptMeditRef index of point data field that contains references
/// (field named medit:ref), -1 if no references
/// @param eltMeditRef index of a cell data field that contains references
/// (field named medit:ref), -1 if no references
/// @param nsols number of point data (except the medit:ref ones)
/// @param metricData 1 if file contains a metric data highlighted by the :metric name
///
/// @return 1 if success, 0 if fail to open/load the file, -1 otherwise;
///
/// I/O at Vtk VTK file format.
///
int MMG5_loadVtkMesh_part1(MMG5_pMesh mesh,const char *filename,vtkDataSet **dataset,
                           int8_t *ptMeditRef,int8_t *eltMeditRef,int *nsols,
                           int8_t *metricData ) {

  (*nsols) = 0;
  (*ptMeditRef) = (*eltMeditRef) = -1;

  // Read all the data from the file
  try {
    (*dataset) = MMG5_load_vtkXMLFile<vtkDataSetReader> ( filename );
  }
  catch( ... ) {
    return 0;
  }

  // count the number of entities of each type
  int ier = MMG5_count_vtkEntities ( (*dataset),mesh,ptMeditRef,eltMeditRef,
                                     nsols,metricData );

  if ( ier != 1 ) {
    return -1;
  }
  return 1;
}

/// @param mesh pointer toward a MMG5 mesh
/// @param filename pointer toward the filename
/// @param dataset vtkdataset structure
/// @param ptMeditRef index of point data field that contains references
/// (field named medit:ref), -1 if no references
/// @param eltMeditRef index of a cell data field that contains references
/// (field named medit:ref), -1 if no references
/// @param nsols number of point data (except the medit:ref ones)
/// @param metricData 1 if file contains a metric data highlighted by the :metric name
///
/// @return 1 if success, 0 if fail to open/load the file, -1 other errors;
///
/// I/O at Vtu VTK file format, part 1: file reading + count of the number of entities.
///
int MMG5_loadVtuMesh_part1(MMG5_pMesh mesh,const char *filename,vtkDataSet **dataset,
                           int8_t *ptMeditRef,int8_t *eltMeditRef,int *nsols,
                           int8_t *metricData) {

  (*nsols) = 0;
  (*ptMeditRef) = (*eltMeditRef) = -1;
  (*metricData) = 0;

  // Read all the data from the file
  try {
    (*dataset) = MMG5_load_vtkXMLFile<vtkXMLUnstructuredGridReader> ( filename );
  }
  catch ( ... ) {
    return 0;
  }

  // count the number of entities of each type
  int ier = MMG5_count_vtkEntities ( (*dataset),mesh,ptMeditRef,eltMeditRef,
                                     nsols,metricData );

  if ( ier != 1 ) {
    return -1;
  }

  return 1;
}



/// @param mesh pointer toward a MMG5 mesh
/// @param ptMeditRef 1 if a point data field contains references (field named medit:ref)
/// @param eltMeditRef 1 if a cell data field contains references (field named medit:ref)
/// @param nsols number of point data (except the medit:ref ones)
///
/// @return 1 if success, -1 if fail.
///
/// I/O at Vtu VTK file format, part 2: mesh and solution storing
///
int MMG5_loadVtkMesh_part2(MMG5_pMesh mesh,MMG5_pSol *sol,vtkDataSet **dataset,
                           int8_t ptMeditRef,int8_t eltMeditRef,int nsols) {
  vtkSmartPointer<vtkDataArray> ptar = NULL, car = NULL;
  int                           ier;
  MMG5_int                      nref = 0;
  static int8_t                 mmgWarn1 = 0;

  // Point transfers in Mmg data structure
  if ( ptMeditRef > -1 ) {
    auto *pd = (*dataset)->GetPointData();
    ptar = pd->GetArray(ptMeditRef);

    // Check the field name
    assert ( strstr( ptar->GetName(), "medit:ref" ) );

    // Check that we get 1 data only
    assert ( ptar->GetNumberOfComponents() == 1 );

    MMG5_int np = ptar->GetNumberOfTuples();
    if ( np != mesh->np ) {
      printf( "  ## Error: Point data size (%" MMG5_PRId ") differs from the number of"
              " vertices (%" MMG5_PRId ")\n",np,mesh->np);
      return -1;
    }
    // read vertices and vertices refs
    for ( vtkIdType k = 0; k < (*dataset)->GetNumberOfPoints(); k++ ) {
      MMG5_pPoint ppt = &mesh->point[k+1];
      (*dataset)->GetPoint(k,mesh->point[k+1].c);
      ppt->tag  = MG_NUL;
      ppt->tmp  = 0;
      ppt->ref  = ptar->GetTuple1(k);
      if ( ppt->ref < 0 ) {
        ppt->ref = -ppt->ref;
        ++nref;
      }
    }
  }
  else {
    // read vertices only
    for ( vtkIdType k = 0; k < (*dataset)->GetNumberOfPoints(); k++ ) {
      MMG5_pPoint ppt = &mesh->point[k+1];
      (*dataset)->GetPoint(k,mesh->point[k+1].c);
      ppt->tag  = MG_NUL;
      ppt->tmp  = 0;
      ppt->ref  = 0;
    }
  }

  // Cells transfer toward Mmg data structure
  assert ( !mesh->nei );
  assert ( !mesh->nti );
  assert ( !mesh->nai );
  assert ( (mesh->npi == mesh->np) || !mesh->npi );

  mesh->npi = 0;
  MMG5_int nqi   = 0;
  MMG5_int npri  = 0;
  MMG5_int na    = 0;
  MMG5_int nbl_a = 0;
  MMG5_int nt    = 0;
  MMG5_int nbl_t = 0;

  // Get pointer toward cells data containing element refs
  vtkIdType numCells = (*dataset)->GetNumberOfCells();

  if ( eltMeditRef > -1 ) {
    // If present, read elements references
    car = (*dataset)->GetCellData()->GetArray(eltMeditRef);

    // Check the field name
    assert ( strstr( car->GetName(), "medit:ref" ) );
    // Check that we get 1 data only
    assert ( car->GetNumberOfComponents() == 1 );

    MMG5_int ne = car->GetNumberOfTuples();
    if ( ne != numCells ) {
      printf( "  ## Error: Cell data size (%" MMG5_PRId ") differs from the number of"
              " cells (%lld)\n",ne,numCells);
      return -1;
    }
  }

  // Transfer cells to Mmg
  for ( vtkIdType k = 0; k < numCells; k++ ) {
    MMG5_pPoint ppt = NULL;
    MMG5_pEdge  pa  = NULL;
    MMG5_pTria  ptt = NULL;
    MMG5_pQuad  pq  = NULL;
    MMG5_pTetra pt  = NULL;
    MMG5_pPrism ppr = NULL;

    int typ = (*dataset)->GetCellType(k);
    MMG5_int ref = 0;

    switch ( typ ) {
    case ( VTK_VERTEX ):
      ppt = &mesh->point[++mesh->npi];
      ppt->ref = car ? car->GetTuple1(k) : 0;
      if ( ppt->ref < 0 ) {
        ppt->ref = -ppt->ref;
        ++nref;
      }
      break;
    case ( VTK_POLY_LINE ):
      int n;
      n = (*dataset)->GetCell(k)->GetNumberOfPoints() - 1;

      for ( int i=0; i<n; ++i ) {
        ++mesh->nai;
        ref = car ? car->GetTuple1(k) : 0;
      }
      /* Skip edges with iso ref */
      if ( mesh->info.iso &&  MMG5_abs(ref) == mesh->info.isoref ) {
        /* Skip this edge */
        ++nbl_a;
      }
      else {
        for ( int i=0; i < n; ++i ) {
          pa = &mesh->edge[++na];
          pa->a = (*dataset)->GetCell(k)->GetPointId(i)  +1;
          pa->b = (*dataset)->GetCell(k)->GetPointId(i+1)+1;
        }
        pa->ref = ref;
        pa->tag |= MG_REF;
        if ( pa->ref < 0 ) {
          pa->ref = -pa->ref;
          ++nref;
        }
      }

      break;
    case ( VTK_LINE ):

      ++mesh->nai;
      ref = car ? car->GetTuple1(k) : 0;

      // Skip edges with iso ref
      if ( mesh->info.iso &&  MMG5_abs(ref) == mesh->info.isoref ) {
        /* Skip this edge */
        ++nbl_a;
      }
      else {
        pa = &mesh->edge[++na];
        pa->a = (*dataset)->GetCell(k)->GetPointId(0)+1;
        pa->b = (*dataset)->GetCell(k)->GetPointId(1)+1;
        pa->ref = ref;
        pa->tag |= MG_REF;
        if ( pa->ref < 0 ) {
          pa->ref = -pa->ref;
          ++nref;
        }
      }
      assert( na+nbl_a <= mesh->na );
      break;

    case ( VTK_TRIANGLE ):
      ++mesh->nti;
      ref = car ? car->GetTuple1(k) : 0;

      // skip tria with iso ref in 3D
      if ( mesh->info.iso && MMG5_abs(ref) == mesh->info.isoref && mesh->dim == 3 ) {
        /* Skip this edge */
        ++nbl_t;
      }
      else {
        ptt = &mesh->tria[++nt];
        for ( int i=0; i<3; ++i ) {
          ptt->v[i] = (*dataset)->GetCell(k)->GetPointId(i)+1;
        }
        ptt->ref = ref;

        if ( ptt->ref < 0 ) {
          ptt->ref = -ptt->ref;
          ++nref;
        }
      }
      break;

    case ( VTK_QUAD ):
      pq = &mesh->quadra[++nqi];
      pq->ref = car ? car->GetTuple1(k) : 0;
      for ( int i=0; i<4; ++i ) {
        pq->v[i] = (*dataset)->GetCell(k)->GetPointId(i)+1;
      }
      if ( pq->ref < 0 ) {
        pq->ref = -pq->ref;
        ++nref;
      }
      break;

    case ( VTK_TETRA ):
      // Skip volume elts if called from mmgs
      if ( !mesh->tetra ) continue;

      pt = &mesh->tetra[++mesh->nei];
      pt->ref = car ? car->GetTuple1(k) : 0;
      for ( int i=0; i<4; ++i ) {
        pt->v[i] = (*dataset)->GetCell(k)->GetPointId(i)+1;
      }
      if ( pt->ref < 0 ) {
        pt->ref = -pt->ref;
        ++nref;
      }
      break;

    case ( VTK_WEDGE ):
      // Skip volume elts if called from mmgs
      if ( !mesh->prism ) continue;

      ppr = &mesh->prism[++npri];
      ppr->ref = car ? car->GetTuple1(k) : 0;
      for ( int i=0; i<6; ++i ) {
        ppr->v[i] = (*dataset)->GetCell(k)->GetPointId(i)+1;
      }
      if ( ppr->ref < 0 ) {
        ppr->ref = -ppr->ref;
        ++nref;
      }
      break;

    default:
      if ( !mmgWarn1 ) {
        printf("  ## Warning:%s: unexpected element type (%d).\n",
               __func__,typ);
        mmgWarn1 = 1;
      }
    }
  }

  // Checks
  assert ( mesh->nai  == mesh->na );
  assert ( mesh->nti  == mesh->nt );
  assert ( mesh->nei  == mesh->ne );
  assert ( nqi        == mesh->nquad );
  assert ( npri       == mesh->nprism );
  assert ( na + nbl_a == mesh->na );
  assert ( nt + nbl_t == mesh->nt );

  // Array reallocation if ISO refs has been skipped
  if (  mesh->info.iso ) {
    if ( mesh->nt ) {
      if( !nt )
        MMG5_DEL_MEM(mesh,mesh->tria);

      else if ( nt < mesh->nt ) {
        MMG5_ADD_MEM(mesh,(nt-mesh->nt)*sizeof(MMG5_Tria),"triangles",
                     fprintf(stderr,"  Exit program.\n");
                     return -1);
        MMG5_SAFE_RECALLOC(mesh->tria,mesh->nt+1,(nt+1),MMG5_Tria,"triangles",
                           return -1);
      }
      mesh->nt = nt;
    }
    if ( mesh->na ) {
      if( !na )
        MMG5_DEL_MEM(mesh,mesh->edge);
      else if ( na < mesh->na ) {
        MMG5_ADD_MEM(mesh,(na-mesh->na)*sizeof(MMG5_Edge),"edges",
                     fprintf(stderr,"  Exit program.\n");
                     return -1);
        MMG5_SAFE_RECALLOC(mesh->edge,mesh->na+1,(na+1),MMG5_Edge,"edges",
                           return -1);
      }
      mesh->na = na;
    }
  }

  ier = MMG5_check_readedMesh(mesh,nref);
  if ( ier < 1 ) return ier;

  if ( sol && *sol ) {
    // Read the solution at nodes
    // Init (*sol)[0] for the case where nsols=0
    MMG5_pSol psl = *sol;
    psl->ver = mesh->ver;
    psl->dim = mesh->dim;
    psl->type = 1;

    int isol = 0;
    if ( nsols ) {
      auto *pd = (*dataset)->GetPointData();

      auto *cd = (*dataset)->GetCellData();

      if ( pd ) {
        int npointData = pd->GetNumberOfArrays();

        for (int j = 0; j < npointData; j++) {
          char *ptr = NULL;
          bool metricData = 0;
          char chaine[MMG5_FILESTR_LGTH];
          strcpy(chaine,pd->GetArrayName(j));

          if  ( strstr(chaine,"medit:ref" ) ) {
            continue;
          }
          else if ( (ptr = strstr(chaine,":metric")) ) {
            *ptr = '\0';
            metricData = 1;
          }

          psl = *sol + isol;
          psl->ver = mesh->ver;
          psl->dim = mesh->dim;
          psl->type = 1;
          psl->entities = MMG5_Vertex;

          if ( !MMG5_Set_inputSolName(mesh,psl,chaine) ) {
            if ( !mmgWarn1 ) {
              mmgWarn1 = 1;
              fprintf(stderr,"\n  ## Warning: %s: unable to set solution name for"
                      " at least 1 solution.\n",__func__);
            }
          }

          auto ar = pd->GetArray(j);

          psl->np = ar->GetNumberOfTuples();
          if ( mesh->np != psl->np ) {
            fprintf(stderr,"  ** MISMATCHES DATA: THE NUMBER OF VERTICES IN "
                    "THE MESH (%" MMG5_PRId ") DIFFERS FROM THE NUMBER OF VERTICES IN "
                    "THE SOLUTION (%" MMG5_PRId ") \n",mesh->np,psl->np);
            return -1;
          }

          int ncp = ar->GetNumberOfComponents();

          if ( ncp == 1 ) {
            psl->size = 1;
            psl->type = 1;
          }
          else if ( ncp == 2 || ncp == 3 ) {
            assert ( ncp == mesh->dim );
            psl->size = ncp;
            psl->type = 2;
          }
          else if ( ncp == (mesh->dim * mesh->dim) ) {
            if ( metricData ) {
              psl->size = (psl->dim*(psl->dim+1))/2;
              psl->type = 3;
            }
            else {
              psl->size = ncp;
              psl->type = 4;
            }
          }
          else {
            fprintf(stderr,"  ** UNEXPECTED NUMBER OF COMPONENTS (%d). IGNORED \n",ncp);
            return -1;
          }

          // mem alloc
          if ( psl->m )  MMG5_DEL_MEM(mesh,psl->m);
          psl->npmax = mesh->npmax;

          MMG5_ADD_MEM(mesh,(psl->size*(psl->npmax+1))*sizeof(double),"initial solution",
                       fprintf(stderr,"  Exit program.\n");
                       return -1);
          MMG5_SAFE_CALLOC(psl->m,psl->size*(psl->npmax+1),double,return -1);

          switch ( psl->type ) {
          case ( 1 ): case ( 2 ):
            for (MMG5_int k=1; k<=psl->np; k++) {
              MMG5_int iadr = k*psl->size;
              ar->GetTuple ( k-1, &psl->m[iadr] );
            }
            break;

          case ( 3 ):
            // anisotropic sol
            double dbuf[9];

            for (MMG5_int k=1; k<=psl->np; k++) {
              ar->GetTuple ( k-1, dbuf );
              MMG5_int iadr = psl->size*k;

              if ( !metricData ) {
                // Non symmetric tensor
                if ( psl->dim ==2 ) {
                  psl->m[iadr] = dbuf[0];
                  psl->m[iadr+1] = dbuf[1];
                  psl->m[iadr+2] = dbuf[3];
                  psl->m[iadr+3] = dbuf[4];
                }
                else {
                  for ( int i=0 ; i<9 ; i++ ) {
                    psl->m[iadr+i] = dbuf[i];
                  }
                }
              }
              else {
                // Symmetric tensor
                if ( psl->dim ==2 ) {
                  assert ( dbuf[1] == dbuf[2] );

                  psl->m[iadr] = dbuf[0];
                  psl->m[iadr+1] = dbuf[1];
                  psl->m[iadr+2] = dbuf[3];
                }
                else {
                  assert ( dbuf[1]==dbuf[3] || dbuf[2]==dbuf[6] || dbuf[5]==dbuf[7] );

                  psl->m[iadr+0] = dbuf[0];
                  psl->m[iadr+1] = dbuf[1];
                  psl->m[iadr+2] = dbuf[2];
                  psl->m[iadr+3] = dbuf[4];
                  psl->m[iadr+4] = dbuf[5];
                  psl->m[iadr+5] = dbuf[8];
                }
              }
            }

            break;
          default:
            fprintf(stderr,"  ** UNEXPECTED METRIC TYPE (%d). EXIT PROGRAM \n",psl->type);
            return -1;
          }
          ++isol;
        }
      }

      if ( cd ) {
        int ncellData = cd->GetNumberOfArrays();

        for (int j = 0; j < ncellData; j++) {
          char *ptr = NULL;
          char chaine[MMG5_FILESTR_LGTH];
          strcpy(chaine,cd->GetArrayName(j));

          if  ( strstr(chaine,"medit:ref" ) ) {
            continue;
          }

          psl = *sol + isol;
          psl->ver = mesh->ver;
          psl->dim = mesh->dim;
          psl->type = 1;
          psl->entities = MMG5_Tetrahedron;

          if ( !MMG5_Set_inputSolName(mesh,psl,chaine) ) {
            if ( !mmgWarn1 ) {
              mmgWarn1 = 1;
              fprintf(stderr,"\n  ## Warning: %s: unable to set solution name for"
                      " at least 1 solution.\n",__func__);
            }
          }

          auto ar = cd->GetArray(j);

          psl->np = ar->GetNumberOfTuples();
          if ( numCells != psl->np ) {
            fprintf(stderr,"  ** MISMATCHES DATA: THE NUMBER OF ELEMENTS IN "
                    "THE MESH (%" MMG5_PRId ") DIFFERS FROM THE NUMBER OF CELLS DATA IN "
                    "THE SOLUTION (%" MMG5_PRId ") \n",mesh->ne,psl->np);
            return -1;
          }

          int ncp = ar->GetNumberOfComponents();

          if ( ncp == 1 ) {
            psl->size = 1;
            psl->type = 1;
          }
          else if ( ncp == 2 || ncp == 3 ) {
            assert ( ncp == mesh->dim );
            psl->size = ncp;
            psl->type = 2;
          }
          else if ( ncp == (mesh->dim * mesh->dim) ) {
            psl->size = ncp;
            psl->type = 3;
          }
          else {
            fprintf(stderr,"  ** UNEXPECTED NUMBER OF COMPONENTS (%d). IGNORED \n",ncp);
            return -1;
          }

          // mem alloc
          if ( psl->m )  MMG5_DEL_MEM(mesh,psl->m);
          psl->npmax = mesh->nemax;

          MMG5_ADD_MEM(mesh,(psl->size*(psl->npmax+1))*sizeof(double),"initial solution",
                       fprintf(stderr,"  Exit program.\n");
                       return -1);
          MMG5_SAFE_CALLOC(psl->m,psl->size*(psl->npmax+1),double,return -1);

          switch ( psl->type ) {
          case ( 1 ): case ( 2 ):
            for (MMG5_int k=1; k<=psl->np; k++) {
              MMG5_int iadr = k*psl->size;
              ar->GetTuple ( k-1, &psl->m[iadr] );
            }
            break;

          case ( 3 ):
            // anisotropic sol
            double dbuf[9];

            for (MMG5_int k=1; k<=psl->np; k++) {
              ar->GetTuple ( k-1, dbuf );
              MMG5_int iadr = psl->size*k;

              // Non symmetric tensor
              if ( psl->dim ==2 ) {
                psl->m[iadr] = dbuf[0];
                psl->m[iadr+1] = dbuf[1];
                psl->m[iadr+2] = dbuf[3];
                psl->m[iadr+3] = dbuf[4];
              }
              else {
                for ( int i=0 ; i<9 ; i++ ) {
                  psl->m[iadr+i] = dbuf[i];
                }
              }
            }

            break;
          default:
            fprintf(stderr,"  ** UNEXPECTED METRIC TYPE (%d). EXIT PROGRAM \n",psl->type);
            return -1;
          }

          ++isol;
        }
      }
    }
  }


  (*dataset)->Delete();

  return 1;
}
