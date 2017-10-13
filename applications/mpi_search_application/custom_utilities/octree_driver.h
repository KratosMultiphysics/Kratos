#pragma once

#include <math.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <stdint.h>
#include <limits.h>
#include <bitset>
#include <string.h>
#include "octree_binary.h"
#include "octree_binary_cell.h"
#include <assert.h>

using namespace std;
using namespace Kratos;

#define KRATOS_INDEPENDENT

//relation between lineal position of nodes and position of sons. Lineal positions are numbered following GiD criteria for hexas (right hand rule).
const int sons_positions[8] = {0,1,3,2,4,5,7,6};

const int OCTREE_MAX_LEVELS=30;//max number of levels for octree. This is to avoid infinite recursions.

struct DataInsideOctreeCell
{
	char is_outer_cell;
	int rank;
  //==0 means that the cell intersects triangles
  //==n>0 values indicate that the cell is outer, and it is n levels far from a 0 cell
  //==n<0 values indicate that the cell is inner, and it is abs(n) levels far from a 0 cell
  //GetRank and SetRank functions refer to the abs(is_outer_cell) value

public:
  DataInsideOctreeCell();
  ~DataInsideOctreeCell();

  char _GetIsOuterCell() const {return is_outer_cell;}
  void _SetIsOuterCell(char val) {is_outer_cell=val;}

};

class Node{
  double coords[3];

public:
	void SetCoord(const int ipos,const double val) {coords[ipos] = val;}
	double GetCoord(int ipos){
		return coords[ipos];
	};
	void MakeTest(double x,double y,double z){
		this->SetCoord(0,x);
		this->SetCoord(1,y);
		this->SetCoord(2,z);
	}
};

class Triangle
{
  Node* nodes[3];
public:

  Triangle(Node* nod1, Node* nod2, Node* nod3)  {
    nodes[0] = nod1; nodes[1] = nod2; nodes[2] = nod3;
  }
  ~Triangle(){
	}

  bool Intersects(const double* bbox_min, const double* bbox_max,const double tolerance);
  void CalcBoundingBox(double* min_corner, double* max_corner);
	double GetCoord(int inode,int ipos)
	{
		return nodes[inode]->GetCoord(ipos);
	}
};

class Boundary
{
	 Triangle *tri;
	 FILE *fp=NULL;
public:
	Boundary(const char name[100]){
		fp=fopen((char*)name,"r");
		if(!fp)
		{
			cout<<"File not opened-Program finished\n";
			assert(0);
		}
		ReadFileStl(fp,tri);
	}
	~Boundary(){
		delete tri;
	}
	void ReadFileStl(FILE *fp,Triangle *tri);
};

class OctreeConfigure {
public:
  enum {
    CHILDREN_NUMBER = 8,
    DIMENSION = 3,
    MAX_LEVEL = OCTREE_MAX_LEVELS,
    MIN_LEVEL = 2 // this cannot be less than 2!!! Pooyan
  };
  typedef DataInsideOctreeCell data_type;
  typedef double* point_coords;
  typedef Triangle* pointer_type;
  
  static data_type* AllocateData() {
    return new data_type;
  }

  static void CopyData(data_type* source, data_type* destination) {
    destination = source;
  }

  static void DeleteData(data_type* data) {
    delete data;
  }        

  static int IsIntersected(const pointer_type elem,const double tolerance,
    const point_coords min_coord, const point_coords max_coord){
      return elem->Intersects(min_coord,max_coord,tolerance);
  }

  static void GetBoundingBox(pointer_type elem,point_coords min_coord, point_coords max_coord){
    elem->CalcBoundingBox(min_coord,max_coord);
  }
};


typedef OctreeBinaryCell<OctreeConfigure> OctreeCell;
typedef OctreeBinary<OctreeCell> Octree;

typedef std::vector<OctreeCell*> OctreeCell_vector;

//class GiDOM_Octree. This class is like an interafce to use OctreeBinary function, taking into account the needs and data of the mesher
struct OctreeDriver
{
  Octree* octree_bin;
  int *rank_leaves;
public:

  OctreeDriver();
  ~OctreeDriver();

  Octree* GetOctree() const {return octree_bin;}

  bool CheckIsBalanced() const {return octree_bin->CheckConstrain2To1();}
  int BalanceOctree();

  OctreeCell* GetRoot() {return octree_bin->pGetRoot();};
  OctreeCell* GetOctreeCell(const double *coord) const {return octree_bin->pGetCellNormalized(coord);}
  OctreeCell* GetOctreeCell(Octree::key_type* keys) const {return octree_bin->pGetCell(keys);}

  int CalcAllLeavesVector(OctreeCell_vector*& all_leaves);
  void GetLeavesInBoundingBox(const double* coord1, const double* coord2,
    OctreeCell_vector& leaves) const {
      return octree_bin->GetLeavesInBoundingBoxNormalized(coord1,coord2,leaves);
  }

  int SubdivideCell(OctreeCell* cell) const; 

  int RefineOctreeWithSize(const double size) {     
    int fail=octree_bin->RefineWithUniformSizeNormalized(size);     
    return fail;
  }

  bool IsCellEmpty(const OctreeCell* cell) const
  {
    if (!cell->pGetObjects()) return true;
    return !cell->pGetObjects()->size();    
  }

  OctreeCell* GetNeighbour(const OctreeCell* cell, const int idir) const {    
    OctreeCell* ret= octree_bin->pGetNeighbourCell(cell,idir);     
    return ret;
  }

  OctreeCell* GetLeftNeighbour(const OctreeCell* cell) const {
    OctreeCell* ret= octree_bin->pGetLeftCell(cell);
    return ret;
  }

  OctreeCell* GetRightNeighbour(const OctreeCell* cell) const {
    OctreeCell* ret= octree_bin->pGetRightCell(cell);
    return ret;      
  }
  OctreeCell* GetTopNeighbour(const OctreeCell* cell) const {
    OctreeCell* ret= octree_bin->pGetTopCell(cell);
    return ret;      
  }
  OctreeCell* GetBottomNeighbour(const OctreeCell* cell) const {
    OctreeCell* ret= octree_bin->pGetBottomCell(cell);
    return ret;     
  }
  OctreeCell* GetFrontNeighbour(const OctreeCell* cell) const{
    OctreeCell* ret= octree_bin->pGetFrontCell(cell);
    return ret;
  }
  OctreeCell* GetBackNeighbour(const OctreeCell* cell) const{
    OctreeCell* ret= octree_bin->pGetBackCell(cell);
    return ret;
  }

	void GetCoordOctreePosition(const OctreeCell* cell,int ipos,double* coord_point) const;
	double GetCoordinate(Octree::key_type key) const {return GetOctree()->GetCoordinateNormalized(key);}
	double CalcSize(const OctreeCell* cell) const {return GetOctree()->CalcSizeNormalized(cell);}
	void PrintMeshGiD(int rank);
	void SubdivideCellsIntersected(const int n_levels,Triangle *tri, Octree *octree, int rank);
	void AsignRankLeaves( Octree *octree, 	OctreeCell_vector leaves);
	void BalanceOctreeMPI( Octree *octree,int rank);
};
size_t InvertKey(size_t key);
size_t GetMortonKey(size_t key_x,size_t key_y,size_t key_z);
void QuickSort(size_t *vec,int *position,int n);
void Sort_Single(size_t *vec,int *position,int L_izq,int L_der);
int GetRankLeaf(size_t *keys);


//this vector indicate the position (in father) corresponding to the iposition in the linear positions in childs
const int relative_pos_from_cell_to_son[8][8]=
{
  { 0,9,21,12,13,22,8,25  },//ison 0
  { 9,1,10,21,22,14,23,8  },//ison 1
  { 12,21,11,3,25,8,24,16 },//ison 2
  { 21,10,2,11,8,23,15,24 },//ison 3
  { 13,22,8,25,4,17,26,20 },//ison 4
  { 22,14,23,8,17,5,18,26 },//ison 5
  { 25,8,24,16,20,26,19,7 },//ison 6
  { 8,23,15,24,26,18,6,19 } //ison 7
};

const int aristes_cell[12][3]=
{
  { 0,9,1 },
  { 1,10,2},
  { 2,11,3},
  { 3,12,0},
  { 0,13,4},
  { 1,14,5},
  { 2,15,6},
  { 3,16,7},
  { 4,17,5},
  { 5,18,6},
  { 6,19,7},
  { 7,20,4} 
};

const int aristes_cell_faces[6][8]=
{//sorted following right hand rule (to match with the sorting of numbers of center of faces nodes
  { 0,9,1,10,2,11,3,12  },//bottom face
  { 0,13,4,17,5,14,1,9  },//back face
  { 1,14,5,18,6,15,2,10 },//right face
  { 2,15,6,19,7,16,3,11 },//front face
  { 3,16,7,20,4,13,0,12 },//left face
  { 4,20,7,19,6,18,5,17 } //upper face
};
