#include "octree_driver.h"
#include "geometric_operations.h"
OctreeDriver::OctreeDriver()
{
  octree_bin=new Octree();
}

OctreeDriver::~OctreeDriver()
{
  delete octree_bin;
}


int OctreeDriver::SubdivideCell(OctreeCell* cell) const
{
  assert(cell->GetLevel()>cell->MIN_LEVEL);
  if (!IsCellEmpty(cell)){
    cell->SubdivideCell();
    //TransferTrianglesToSons


  }
  else{
    cell->SubdivideCell();
  }
  return 0;
}
typedef std::size_t key_type;
//generate the coordinates of the node with position 'pos' that belongs to the cell.
void OctreeDriver::GetCoordOctreePosition(const OctreeCell* cell,int ipos,double* coord_point) const
{
  key_type keys[3];         
  cell->GetKey(ipos,keys);
  for (int i=0;i<3;i++){
    coord_point[i] = GetCoordinate(keys[i]);
  } 
}

int OctreeDriver::CalcAllLeavesVector(OctreeCell_vector*& all_leaves)
{//get a consecutive vector of all the leaves. The vector is cretaed inside the function.
  int num_leaves;

  if (!all_leaves) all_leaves=new OctreeCell_vector();
  else all_leaves->clear();
  octree_bin->GetAllLeavesVector(*all_leaves);
  num_leaves=(int)(all_leaves->size());

  return num_leaves;
}


int OctreeDriver::BalanceOctree()
{   
  octree_bin->Constrain2To1();
  assert(octree_bin->CheckConstrain2To1());

  return 0;
}


DataInsideOctreeCell::DataInsideOctreeCell()
{
  is_outer_cell=0;

}

DataInsideOctreeCell::~DataInsideOctreeCell()
{

}
void Boundary::ReadFileStl(FILE *fp,Triangle *tri)
{
	char trash[80];
	float rvalue;
	fscanf(fp,"%s",trash);
	fscanf(fp,"%s",trash);
	int cont=0;
	while(1)
	{
		fscanf(fp,"%s",trash);
		if(trash[0]=='e'&&
			 trash[1]=='n'&&
       trash[2]=='d')
		{
			cout<<"Ya se termino de leer el archivo"<<endl;
			break;
		}
		fscanf(fp,"%s",trash);
		for(int ipos=0;ipos<3;ipos++)
		{
			fscanf(fp,"%f",&rvalue);
		}
		fscanf(fp,"%s",trash);
		fscanf(fp,"%s",trash);
		Node nod[3];
		for(int ipos=0;ipos<3;ipos++)
		{
			fscanf(fp,"%s",trash);
			for(int jpos=0;jpos<3;jpos++)
			{
				fscanf(fp,"%f",&rvalue);
			  nod[ipos].SetCoord(jpos,rvalue);
			}
		}
		fscanf(fp,"%s",trash);
		fscanf(fp,"%s",trash);
		//ALLOCATE INFO FOR VECTOR tri



		cont++;
	}
	fclose(fp);
	assert(0);
}
bool Triangle::Intersects(const double* bbox_min, const double* bbox_max,const double tolerance){
	double boxcenter[3];
	double boxhalfsize[3];
	double triverts[3][3];
	boxcenter[0]   = 0.50 * ( bbox_min[0] + bbox_max[0] );
	boxcenter[1]   = 0.50 * ( bbox_min[1] + bbox_max[1] );
	boxcenter[2]   = 0.50 * ( bbox_min[2] + bbox_max[2] );

	boxhalfsize[0] = 0.50 * ( bbox_max[0] - bbox_min[0] );
	boxhalfsize[1] = 0.50 * ( bbox_max[1] - bbox_min[1] );
	boxhalfsize[2] = 0.50 * ( bbox_max[2] - bbox_min[2] );

	triverts[0][0]=GetCoord(0,0);	triverts[0][0]=GetCoord(0,0);	triverts[0][0]=GetCoord(0,0);
	for(int inode=0;inode<3;inode++){
		for(int ipos=0;ipos<3;ipos++){
			triverts[inode][ipos]=GetCoord(inode,ipos);
		}
	}


	return	DoesOverLapTriBox(boxcenter,boxhalfsize,triverts);
}
void Triangle::CalcBoundingBox(double* min_corner, double* max_corner){

	assert(0);
}
void OctreeDriver::PrintMeshGiD()
{
	ofstream rOStream;
	rOStream.open ("mesh.msh");
	int nleaves=0;
	OctreeCell_vector leaves;
	octree_bin->GetAllLeavesVector( leaves);
	nleaves=(int)(leaves.size());
	cout << "writing " << nleaves << " leaves" << endl;
	rOStream << "MESH \"leaves\" dimension 3 ElemType Hexahedra Nnode 8\n";
	rOStream << "# color 96 96 96" << endl;
	rOStream << "Coordinates" << endl;
	rOStream << "# node number coordinate_x coordinate_y coordinate_z  " << endl;
	size_t node_index = 1;
	for (size_t i = 0; i < nleaves; i++) {
		OctreeCell* cell = leaves[i];
		double min_point[3];
		double max_point[3];
		GetCoordOctreePosition(cell,0,min_point);
		GetCoordOctreePosition(cell,6,max_point);
		double cell_size = max_point[0]-min_point[0];

		for (size_t j = 0; j < 2; j++)
			for (size_t k = 0; k < 2; k++)
				for (size_t h = 0; h < 2; h++) {
					rOStream << node_index++ << "  " << min_point[0] + j * cell_size << "  " << min_point[1] + k * cell_size << "  " << min_point[2] + h * cell_size << endl;
				}
	}
	rOStream << "end coordinates" << endl;
	rOStream << "Elements" << endl;
	rOStream << "# element node_1 node_2 node_3 material_number" << endl;

	for (size_t i = 0; i < leaves.size(); i++) {
		if ((leaves[i]->pGetData()))
			rOStream << i + 1 << "  " << 8 * i + 1 << "  " << 8 * i + 2 << "  " << 8 * i + 4 << "  " << 8 * i + 3 << "  " << 8 * i + 5 << "  " << 8 * i + 6 << "  " << 8 * i + 8 << "  " << 8 * i + 7 << "  " << leaves[i]->GetLevel() + 100 << endl;
		else
			rOStream << i + 1 << "  " << 8 * i + 1 << "  " << 8 * i + 2 << "  " << 8 * i + 4 << "  " << 8 * i + 3 << "  " << 8 * i + 5 << "  " << 8 * i + 6 << "  " << 8 * i + 8 << "  " << 8 * i + 7 << "  " << int(leaves[i]->GetLevel()) << endl;

	}
	rOStream << "end elements" << endl;
	rOStream.close();
}
void OctreeDriver::SubdivideCellsIntersected(const int n_levels,Triangle *tri, Octree *octree)
{
	OctreeCell_vector leaves;
	int ilevel=0;
	while(ilevel<n_levels){
		ilevel++;

		octree->GetAllLeavesVector(leaves);//octree

		OctreeCell* leaf;
		const int n_leaves = leaves.size();
		const double tol = 1e-3;
		for (int ileaf=0;ileaf<n_leaves;ileaf++){
			leaf=leaves[ileaf];

			double min_coord[3];
			this->GetCoordOctreePosition(leaf,0,min_coord);//octree_driver
			double max_coord[3];
			this->GetCoordOctreePosition(leaf,6,max_coord);//octree_friver
			const int intersects =   tri->Intersects(min_coord,max_coord,tol);
			if (intersects){
				this->SubdivideCell(leaf);//octree_driver
			}
		}
		leaves.clear();
	}
}

