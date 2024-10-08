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

#ifndef VTKPARSER_HPP
#define VTKPARSER_HPP

#ifdef USE_VTK

#include "mmgcommon_private.h"

#include <vtkSmartPointer.h>
#include <vtkXMLReader.h>
#include <vtkXMLWriter.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPUnstructuredGridWriter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLPPolyDataWriter.h>
#include <vtkDataSetReader.h>
#include <vtkDataSetWriter.h>
#include <vtkPDataSetWriter.h>
#include <vtkDataSet.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>
#include <vtkStructuredGrid.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkFieldData.h>
#include <vtkCellTypes.h>
#include <vtkDataArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>

#include <vtkLine.h>
#include <vtkTriangle.h>
#include <vtkQuad.h>
#include <vtkTetra.h>
#include <vtkWedge.h>
#include <vtkCellArray.h>
#include <typeinfo>

int MMG5_loadVtpMesh_part1(MMG5_pMesh,const char*,vtkDataSet**,int8_t*,int8_t*,int*,int8_t*);
int MMG5_loadVtuMesh_part1(MMG5_pMesh,const char*,vtkDataSet**,int8_t*,int8_t*,int*,int8_t*);
int MMG5_loadVtkMesh_part1(MMG5_pMesh,const char*,vtkDataSet**,int8_t*,int8_t*,int*,int8_t*);

int MMG5_loadVtkMesh_part2(MMG5_pMesh,MMG5_pSol*,vtkDataSet**,int8_t,int8_t,int);

/// @param d vtk data type in which we want to store the array \a ca
/// @param ca vtk cell array containing the lines connectivity
///
/// Store a list of vtk cells containing only lines into a vtkPolyData and reset
/// the list of cells. Otherwise, lines are casted into polygons by the SetPoly
/// method.
///
static void MMG5_internal_VTKSetLine(vtkSmartPointer<vtkPolyData> d,vtkSmartPointer<vtkCellArray> *ca) {

  d->SetLines(*ca);

  *ca = NULL;
  *ca = vtkSmartPointer< vtkCellArray >::New();

}

/// @param d vtk data type in which we want to store the array \a ca
/// @param ca vtk cell array containing the lines connectivity
///
/// Nothing to do here, the SetCell method allow to handle both lines and elements.
///
static void MMG5_internal_VTKSetLine(vtkSmartPointer<vtkUnstructuredGrid> d,vtkSmartPointer<vtkCellArray> *ca) {

}


/// @param t array of integer storing the cell types
/// @param d vtk data type in which we want to store the array \a ca
/// @param ca vtk cell array containing the cells connectivity
///
/// Store a list of vtk cells into a vtkPolyData
///
static void MMG5_internal_VTKSetCells(int * t,vtkSmartPointer<vtkPolyData> d,vtkSmartPointer<vtkCellArray> ca) {
  d->SetPolys(ca);
}

/// @param t array of integer storing the cell types
/// @param d vtk data type in which we want to store the array \a ca
/// @param ca vtk cell array containing the cells connectivity
///
/// Store a list of vtk cells into a vtkUnstructuredGrid
///
static void MMG5_internal_VTKSetCells(int * t,vtkSmartPointer<vtkUnstructuredGrid> d,vtkSmartPointer<vtkCellArray> ca) {
  d->SetCells(t,ca);
}

/// @param w vtk writer
/// @param binary 1 if we want to save in binary format
///
/// Try to set the suitable file format to the vtk writer
///
static void MMG5_internal_VTKbinary(vtkXMLUnstructuredGridWriter *w, int binary) {
  if ( binary ) {
    w->SetDataModeToBinary();
  }
  else {
    w->SetDataModeToAscii();
  }
}

/// @param w vtk writer
/// @param binary 1 if we want to save in binary format
///
/// Try to set the suitable file format to the vtk writer
///
static void MMG5_internal_VTKbinary(vtkXMLPolyDataWriter *w, int binary) {
  if ( binary ) {
    w->SetDataModeToBinary();
  }
  else {
    w->SetDataModeToAscii();
  }
}

/// @param w vtk writer
/// @param binary 1 if we want to save in binary format
///
/// Try to set the suitable file format to the vtk writer
///
static void MMG5_internal_VTKbinary(vtkDataSetWriter *w, int binary) {
  return;
}

/// @tparam T one of the VTK data class.
/// @tparam TWriter one of the VTK writer class.
/// @tparam PWriter one of the parallel VTK writer class.
///
/// @param mesh pointer toward a MMG5 mesh
/// @param sol pointer toward a MMG5 solution array
/// @param mfilename name of the master file to save (if call by master).
/// @param metricData 1 if sol contains a metric array
/// @param binary 1 to save file at binary format (if available in TWriter class)
/// @param npart number of parts of the saving
/// @param myid id of the process (save the file part number myid)
/// @param master id of the master process (save its part file + the master file)
///
/// @return 1 if success, 0 if fail.
///
/// Save a vtk file at .(p)vtp, .(p)vtu or .(p)vtk format.
///
template <class T, class TWriter, class PWriter>
int MMG5_saveVtkMesh_i(MMG5_pMesh mesh,MMG5_pSol *sol,
                       const char *mfilename,
                       int metricData,int binary,
                       int npart, int myid,int master) {
  int hasPointRef = 0, hasCellRef = 0;

  // Transfer points from Mmg to VTK
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

  MMG5_int np = 0;
  for ( MMG5_int k=1; k<=mesh->np; ++k ) {
    MMG5_pPoint ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) ) continue;

    ppt->tmp = np++;

    if ( mesh->dim == 2 ) ppt->c[2] = 0.;
    if ( ppt->ref )       hasPointRef = 1;

    points->InsertNextPoint(ppt->c[0], ppt->c[1], ppt->c[2]);
  }


  auto dataset = vtkSmartPointer<T> :: New();
  dataset->SetPoints(points);


  // VTK types
  vtkSmartPointer<vtkCellArray>  cellArray  = vtkSmartPointer<vtkCellArray>  ::New();
  vtkSmartPointer<vtkWedge>      wedge      = vtkSmartPointer<vtkWedge>      ::New();
  vtkSmartPointer<vtkTetra>      tetra      = vtkSmartPointer<vtkTetra>      ::New();
  vtkSmartPointer<vtkTriangle>   tria       = vtkSmartPointer<vtkTriangle>   ::New();
  vtkSmartPointer<vtkQuad>       quadra     = vtkSmartPointer<vtkQuad>       ::New();
  vtkSmartPointer<vtkLine>       edge       = vtkSmartPointer<vtkLine>       ::New();

  MMG5_int nc = mesh->na+mesh->nt+mesh->nquad+mesh->ne+mesh->nprism;
  MMG5_int ic = 0;

  int* types = NULL;
  MMG5_SAFE_MALLOC ( types, nc, int, return 0 );

  // transfer edges from Mmg to VTK
  for ( MMG5_int k=1; k<=mesh->na; ++k ) {
    MMG5_pEdge pa = &mesh->edge[k];
    if ( !pa || !pa->a ) continue;

    if ( pa->ref ) hasCellRef = 1;

    edge->GetPointIds()->SetId(0, mesh->point[pa->a].tmp);
    edge->GetPointIds()->SetId(1, mesh->point[pa->b].tmp);

    cellArray->InsertNextCell(edge);
    types[ic++] = VTK_LINE;
  }

  MMG5_internal_VTKSetLine(dataset, &cellArray);

 // Transfer trias from Mmg to VTK
 for ( MMG5_int k=1; k<=mesh->nt; ++k ) {
   MMG5_pTria  ptt = &mesh->tria[k];
   if ( !MG_EOK(ptt) ) continue;

   if ( ptt->ref ) hasCellRef = 1;

   for ( int i=0; i<3; ++i ) {
     tria->GetPointIds()->SetId(i, mesh->point[ptt->v[i]].tmp );
   }
   cellArray->InsertNextCell(tria);
   types[ic++] = VTK_TRIANGLE;
 }

 // Transfer quads from Mmg to VTK
 for ( MMG5_int k=1; k<=mesh->nquad; ++k ) {
   MMG5_pQuad  pq = &mesh->quadra[k];
   if ( !MG_EOK(pq) ) continue;

   if ( pq->ref ) hasCellRef = 1;

   for ( int i=0; i<4; ++i ) {
     quadra->GetPointIds()->SetId(i, mesh->point[pq->v[i]].tmp );
   }
   cellArray->InsertNextCell(quadra);
   types[ic++] = VTK_QUAD;
 }

  // Transfer tets from Mmg to VTK
  for ( MMG5_int k=1; k<=mesh->ne; ++k ) {
    MMG5_pTetra pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;

    if ( pt->ref ) hasCellRef = 1;

    for ( int i=0; i<4; ++i ) {
      tetra->GetPointIds()->SetId(i, mesh->point[pt->v[i]].tmp );
    }
    cellArray->InsertNextCell(tetra);
    types[ic++] = VTK_TETRA;
  }

  // Transfer prisms from Mmg to VTK
  for ( MMG5_int k=1; k<=mesh->nprism; ++k ) {
    MMG5_pPrism ppr = &mesh->prism[k];
    if ( !MG_EOK(ppr) ) continue;

    if ( ppr->ref ) hasCellRef = 1;

    for ( int i=0; i<6; ++i ) {
      wedge->GetPointIds()->SetId(i, mesh->point[ppr->v[i]].tmp );
    }
    cellArray->InsertNextCell(wedge);
    types[ic++] = VTK_WEDGE;
  }

  assert ( ic == nc );

  // Transfer cell array into data set
  MMG5_internal_VTKSetCells(types,dataset,cellArray);

  // Transfer references if needed (i.e. the mesh contains non 0 refs)
  if ( hasPointRef ) {
    auto *ar = vtkFloatArray::New();

    ar->SetNumberOfComponents(1);
    ar->SetNumberOfTuples(mesh->np);
    ar->SetName("medit:ref");

    for ( MMG5_int k = 0; k < mesh->np; k++ ) {
      MMG5_pPoint ppt = &mesh->point[k+1];
      if ( !MG_VOK(ppt) ) continue;

      ar->SetTuple1(ppt->tmp,ppt->ref);
    }

    dataset->GetPointData()->AddArray(ar);
  }
  if ( hasCellRef ) {
    auto *ar = vtkFloatArray::New();

    ar->SetNumberOfComponents(1);
    ar->SetNumberOfTuples(nc);
    ar->SetName("medit:ref");

    MMG5_int idx = 0;
    for ( MMG5_int k = 1; k <= mesh->na; k++ ) {
      if ( !mesh->edge[k].a ) continue;
      ar->SetTuple1(idx++,mesh->edge[k].ref);
    }
    for ( MMG5_int k = 1; k <= mesh->nt; k++ ) {
     if ( !MG_EOK(&mesh->tria[k]) ) continue;
      ar->SetTuple1(idx++,mesh->tria[k].ref);
    }
    for ( MMG5_int k = 1; k <= mesh->nquad; k++ ) {
      if ( !MG_EOK(&mesh->quadra[k]) ) continue;
      ar->SetTuple1(idx++,mesh->quadra[k].ref);
    }
    for ( MMG5_int k = 1; k <= mesh->ne; k++ ) {
      if ( !MG_EOK(&mesh->tetra[k]) ) continue;
      ar->SetTuple1(idx++,mesh->tetra[k].ref);
    }
    for ( MMG5_int k = 1; k <= mesh->nprism; k++ ) {
      if ( !MG_EOK(&mesh->prism[k]) ) continue;
      ar->SetTuple1(idx++,mesh->prism[k].ref);
    }

    dataset->GetCellData()->AddArray(ar);
  }

  // Transfer point solutions into data set
  MMG5_pSol   psl   = NULL;
  int         nsols;

  if ( metricData==1 ) {
    if ( sol && *sol && sol[0]->np ) {
      nsols = 1;
    }
    else {
      /* In analysis mode (-noinsert -noswap -nomove), metric is not allocated */
      nsols = 0;
    }
  }
  else {
    nsols = mesh->nsols;
  }

  static int mmgWarn = 0;
  for ( int isol=0; isol<nsols; ++isol) {
    psl = *sol + isol;

    if ( !psl->m ) {
      if ( !mmgWarn ) {
        mmgWarn = 1;
        fprintf(stderr, "  ## Warning: %s: missing data for at least 1 solution."
                " Skipped.\n",__func__);
      }
      continue;
    }

    int ncp;
    if ( psl->size == 1 ) {
      ncp = 1;
    }
    else if ( psl->size == psl->dim ) {
      ncp = psl->dim;
    }
    else {
      ncp = psl->dim*psl->dim;
    }

    auto *ar = vtkDoubleArray::New();

    ar->SetNumberOfComponents(ncp);
    ar->SetNumberOfTuples(mesh->np);

    if ( psl->namein ) {
      char *tmp = MMG5_Get_basename(psl->namein);
      char *data;

      MMG5_SAFE_CALLOC(data,strlen(tmp)+8,char,
                       MMG5_SAFE_FREE ( types ); return 0);

      strcpy(data,tmp);
      free(tmp); tmp = 0;

      if ( metricData ) {
        strcat ( data , ":metric");
      }
      ar->SetName(data);

      MMG5_DEL_MEM(mesh,data);
    }
    else {
      ar->SetName("no_name");
    }

    double* dfmt = NULL;
    MMG5_SAFE_MALLOC ( dfmt, ncp, double,  MMG5_SAFE_FREE ( types );return 0 );

    if ( psl->size!= (psl->dim*(psl->dim+1))/2 ) {
      /* scalar or vector field: data order and size isn't modified */
      for ( MMG5_int k=1; k<=mesh->np; k++) {
        MMG5_pPoint ppt = &mesh->point[k];
        if ( !MG_VOK(ppt) ) continue;

        MMG5_int    iadr  = k*psl->size;
        for ( int i=0; i<psl->size; ++i ) {
          dfmt[i] = psl->m[iadr+i];
        }

        if ( psl->dim==2 && ncp==3 ) {
          dfmt[2] = 0; // z-component for a vector field
        }
        ar->SetTuple(ppt->tmp,dfmt);
      }
    }
    else {
      /* Tensor field: we need to reconstruct the entire symetric matrix */
      for ( MMG5_int k=1; k<=mesh->np; k++) {
        MMG5_pPoint ppt = &mesh->point[k];
        if ( !MG_VOK(ppt) ) continue;

        MMG5_int iadr  = k*psl->size;
        double *d = &psl->m[iadr];
        if ( psl->dim==2 ) {
          dfmt[0] = d[0];
          dfmt[1] = d[1];
          dfmt[2] = d[1];
          dfmt[3] = d[2];
        }
        else {
          double dbuf[6] = {0,0,0,0,0,0};
          if ( metricData ) {
            assert(!mesh->nsols);
            MMG5_build3DMetric(mesh,psl,k,dbuf);
          }
          else {
            for ( int i=0; i<psl->size; i++)  dbuf[i] = psl->m[psl->size*k+i];
          }
          dfmt[0] = dbuf[0];
          dfmt[1] = dbuf[1];
          dfmt[2] = dbuf[2];
          dfmt[3] = dbuf[1];
          dfmt[4] = dbuf[3];
          dfmt[5] = dbuf[4];
          dfmt[6] = dbuf[2];
          dfmt[7] = dbuf[4];
          dfmt[8] = dbuf[5];
        }
        ar->SetTuple(ppt->tmp,dfmt);
      }
    }
    dataset->GetPointData()->AddArray(ar);

    MMG5_SAFE_FREE(dfmt);
  }

  if ( npart ) {
    // distributed IO
    vtkSmartPointer<PWriter> writer = vtkSmartPointer<PWriter>::New();

#if VTK_MAJOR_VERSION <= 5
    writer->SetInput(dataset);
#else
    writer->SetInputData(dataset);
#endif

    writer->SetFileName(mfilename);

    writer->SetNumberOfPieces(npart);
    writer->SetStartPiece(myid);
    writer->SetEndPiece(myid);
    writer->Write();
  }
  else {
    // centralized IO
    vtkSmartPointer<TWriter> writer = vtkSmartPointer<TWriter>::New();

    writer->SetFileName(mfilename);

#if VTK_MAJOR_VERSION <= 5
    writer->SetInput(dataset);
#else
    writer->SetInputData(dataset);
#endif

    //MMG5_internal_VTKbinary(writer,binary);
    writer->Write();
  }

  MMG5_SAFE_FREE(types);

  return 1;
}

/// @tparam T one of the VTK data class.
/// @tparam TWriter one of the VTK writer class.
/// @param mesh pointer toward a MMG5 mesh
/// @param sol pointer toward a MMG5 solution array
/// @param filename name of the input file.
/// @param metricData 1 if sol contains a metric array
/// @param binary 1 to save file at binary format (if available in TWriter class)
///
/// @return 1 if success, 0 if fail.
///
/// Save a vtk file at .vtp, .vtu or .vtk format.
///
template <class T, class TWriter,class PWriter>
int MMG5_saveVtkMesh(MMG5_pMesh mesh,MMG5_pSol *sol,const char *filename,
                     int metricData,int binary) {

  return MMG5_saveVtkMesh_i<T,TWriter,PWriter>( mesh,sol,filename,
                                                metricData,binary,0,0,0);

}

#endif

#endif
