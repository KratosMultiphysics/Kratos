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
	int compare (const void * a, const void * b){
		return ( *(int*)a - *(int*)b );
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
void OctreeDriver::PrintMeshGiD(int rank)
{
	char rank_name[50];
	sprintf(rank_name,"%d",rank);
	ofstream rOStream;
	strcat(rank_name,"mesh.msh");
	rOStream.open (rank_name);
	int nleaves=0;
	OctreeCell_vector leaves;
	octree_bin->GetAllLeavesVector( leaves);
	nleaves=(int)(leaves.size());
	cout <<rank << " writing " << nleaves << " leaves" << endl;
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
void OctreeDriver::SubdivideCellsIntersected(const int n_levels,Triangle *tri, Octree *octree, int rank)
{
	/*
	 * Aqui se puede mejorar calculando en la primer division (los 
	 * primeros 8 hijos del root) los ranks de cada uno, y despues de 
	 * eso solo transferirlos a sus nuevos hijos (CHECAR CON ABEL COMO 
	 * AGREGAR RANK A LA CELDA Y TRANSMITIRLO A HIJOS)
	 * */
	
	OctreeCell_vector leaves;
	int ilevel=0;
	while(ilevel<n_levels){
		ilevel++;
		octree->GetAllLeavesVector(leaves);
		OctreeCell* leaf;
		const int n_leaves = leaves.size();
		AsignRankLeaves(octree,leaves);
		const double tol = 1e-3;
		for (int ileaf=0;ileaf<n_leaves;ileaf++){
			leaf=leaves[ileaf];
			double min_coord[3];
			this->GetCoordOctreePosition(leaf,0,min_coord);
			double max_coord[3];
			this->GetCoordOctreePosition(leaf,6,max_coord);
			const int intersects =   tri->Intersects(min_coord,max_coord,tol);
			if (intersects && (rank_leaves[ileaf+1]==rank)){
				this->SubdivideCell(leaf);
			}
		}
		leaves.clear();
		if(ilevel<n_levels){
			delete rank_leaves;
		}
	}
}
void OctreeDriver::AsignRankLeaves( Octree *octree, OctreeCell_vector leaves)
{
	const int n_leaves = leaves.size();
	rank_leaves = new int [n_leaves+1];
	rank_leaves[0]=n_leaves;
	for(int ileaf=0;ileaf<n_leaves;ileaf++){
		key_type keys[3];
		leaves[ileaf]->GetKey(0,keys);
		keys[0]>>=9;	keys[1]>>=9;	keys[2]>>=9;
		rank_leaves[ileaf+1]=GetRankLeaf(keys);
	}
}
size_t InvertKey(size_t key)
{
	for(int ibit=0;ibit<32;ibit++)
		key |= ( ( key & ( 1UL << ( 63-ibit ) ) ) >> (63 - ( 2*ibit ) ) );
	key<<=32;
	key>>=32;
	return key;
}
size_t GetMortonKey(size_t key_x,size_t key_y,size_t key_z)
{
	int limit=(sizeof(size_t)*CHAR_BIT)/3;
	size_t morton=0;
	for(int ibit=0;ibit<limit;ibit++)
		morton |= ( ( key_x & ( 1UL << ibit ) ) << ( ( 2*ibit )     )) | ( ( key_y & ( 1UL << ibit ) ) << ( ( 2*ibit ) + 1 )) | ( ( key_z & ( 1UL << ibit ) ) << ( ( 2*ibit ) + 2 ));
	return morton;
}
void Sort_Single(size_t *vec,int *position,int L_izq,int L_der)
{
    int izquierda,derecha,pivote,swap_int;
    size_t swap_dou;
    izquierda=L_izq;
    derecha=L_der;
    pivote=vec[(izquierda+derecha)/2];
    do{
        while(vec[izquierda]<pivote&&izquierda<L_der)
        {
            izquierda++;
        }
        while(pivote<vec[derecha]&&derecha>L_izq)
        {
            derecha--;
        }
        if(izquierda<=derecha)
        {
            swap_dou=vec[izquierda];
            vec[izquierda]=vec[derecha];
            vec[derecha]=swap_dou;
            swap_int=position[izquierda];
            position[izquierda]=position[derecha];
            position[derecha]=swap_int;
            izquierda++;
            derecha--;
        }
    }while(izquierda<=derecha);
    if(L_izq<derecha)
    {
        Sort_Single(vec,position,L_izq,derecha);
    }
    if(L_der>izquierda)
    {
        Sort_Single(vec,position,izquierda,L_der);
    }


}
void QuickSort(size_t *vec,int *position,int n)
{
    Sort_Single(vec,position,0,n-1);
}
int GetRankLeaf(size_t *keys)
{
	int rank=-1;
	size_t mask=(1UL<<19), v1=0,v2=0,v3=0;
	v1 = keys[ 0 ] & mask;
	v2 = keys[ 1 ] & mask;
	v3 = keys[ 2 ] & mask;
	if(v1==0){
		if(v2==0){
			if(v3==0){
				rank=0;
			}else{
				rank=4;
			}
		}else{
			if(v3==0){
				rank=3;
			}else{
				rank=7;
			}
		}
	}else
	{
		if(v2==0){
			if(v3==0){
				rank=1;
			}else{
				rank=5;
			}
		}else{
			if(v3==0){
				rank=2;
			}else{
				rank=6;
			}
		}
	}
	return rank;
}

void OctreeDriver::BalanceOctreeMPI(Octree *octree,int rank)
{
	OctreeCell_vector leaves;
	OctreeCell_vector next_leaves;
	//when the function will be at upper level (in mesher instead of octree) this vector (leaves) should be passed and copied instead of recomputed
	octree->GetAllLeavesVector(leaves);
	for (char i_level = octree_bin->MIN_LEVEL; i_level < octree_bin->ROOT_LEVEL - 1; i_level++) {
		for (size_t i_cell = 0; i_cell < leaves.size(); i_cell++) {
			key_type keys[3];
			leaves[i_cell]->GetKey(0,keys);
			keys[0]>>=9;	keys[1]>>=9;	keys[2]>>=9;
			if(GetRankLeaf(keys)==rank){
				OctreeCell* current_leaf = leaves[i_cell];
				if (current_leaf->GetLevel() == i_level) {
					key_type neighbour_key[3];
					//18 is the number of neighbours counting faces and edges of the cell
					for (int i_direction = 0; i_direction < 18; i_direction++) {
						if (current_leaf->GetNeighbourKey(i_direction, neighbour_key)) {
							OctreeCell* neighbour_cell = octree_bin->pGetCell(neighbour_key);
							key_type key_copy[3];
							key_copy[0]=neighbour_key[0]>>9;	key_copy[1]=neighbour_key[1]>>9;	key_copy[2]=neighbour_key[2]>>9;
							if(rank==GetRankLeaf(key_copy))
							{
								if (neighbour_cell->GetLevel() > i_level + 1) {
									OctreeCell* temp_neighbour_cell = neighbour_cell;
									for (char j_level = neighbour_cell->GetLevel(); j_level > i_level + 1; j_level--) {
										SubdivideCell(temp_neighbour_cell);
										temp_neighbour_cell->TransferObjectsToSonsNormalized();
										size_t child_index = temp_neighbour_cell->GetChildIndex(neighbour_key[0], neighbour_key[1], neighbour_key[2]);
										for (size_t j = 0; j < temp_neighbour_cell->CHILDREN_NUMBER; j++) {
											if (j != child_index) {
												next_leaves.push_back(temp_neighbour_cell->GetChildren() + j);
											}
										}
										temp_neighbour_cell = temp_neighbour_cell->GetChildren() + child_index;
										if (j_level == neighbour_cell->GetLevel() - 1) // the last loop we add all the child as leaf
											next_leaves.push_back(temp_neighbour_cell);
									}
								}
							}
						}
					}
				} else if (current_leaf->IsLeaf()) { // becuase it may be divided
					next_leaves.push_back(current_leaf);
				}
			}
		}
		leaves.swap(next_leaves);
		next_leaves.clear();
	}
	#ifdef KRATOS_INDEPENDENT
	#else
		KRATOS_WATCH(number_of_leaves_);
	#endif
}

		//bitset<64>x(morton_key[ileaf]);
		//cout<<x<<endl;















