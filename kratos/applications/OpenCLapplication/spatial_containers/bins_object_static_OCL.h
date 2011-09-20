//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: clabra $
//   Date:                $Date: 2007-03-27 17:02:19 $
//   Revision:            $Revision: 1.1.1.1 $
//
//

#if !defined(KRATOS_BINS_STATIC_OBJECTS_OCL_CONTAINER_H_INCLUDE)
#define KRATOS_BINS_STATIC_OBJECTS_OCL_CONTAINER_H_INCLUDE

#if defined( __GNUC__) && ! defined( __STRICT_ANSI__ )
#define KRATOS_OCL_4_ARRAY_X(Arr,p) 	(Arr[p].x)
#define KRATOS_OCL_4_ARRAY_Y(Arr,p) 	(Arr[p].y)
#define KRATOS_OCL_4_ARRAY_Z(Arr,p) 	(Arr[p].z)
#define KRATOS_OCL_4_ARRAY_W(Arr,p) 	(Arr[p].w)

#define KRATOS_OCL_4_X(Arr) 		(Arr.x)
#define KRATOS_OCL_4_Y(Arr)		(Arr.y)
#define KRATOS_OCL_4_Z(Arr)		(Arr.z)
#define KRATOS_OCL_4_W(Arr)		(Arr.w)

#define KRATOS_OCL_4_ITR_X(Arr) 	((*Arr).x)
#define KRATOS_OCL_4_ITR_Y(Arr)		((*Arr).y)
#define KRATOS_OCL_4_ITR_Z(Arr)		((*Arr).z)
#define KRATOS_OCL_4_ITR_W(Arr)		((*Arr).w)
#else
#define KRATOS_OCL_4_ARRAY_X(Arr,p) 	(Arr[p].s[0])
#define KRATOS_OCL_4_ARRAY_Y(Arr,p)	(Arr[p].s[1])
#define KRATOS_OCL_4_ARRAY_Z(Arr,p) 	(Arr[p].s[2])
#define KRATOS_OCL_4_ARRAY_W(Arr,p) 	(Arr[p].s[3])

#define KRATOS_OCL_4_X(Arr) 		(Arr.s[0])
#define KRATOS_OCL_4_Y(Arr) 		(Arr.s[1])
#define KRATOS_OCL_4_Z(Arr) 		(Arr.s[2])
#define KRATOS_OCL_4_W(Arr) 		(Arr.s[3])

#define KRATOS_OCL_4_ITR_X(Arr) 	((*Arr).s[0])
#define KRATOS_OCL_4_ITR_Y(Arr)		((*Arr).s[1])
#define KRATOS_OCL_4_ITR_Z(Arr)		((*Arr).s[2])
#define KRATOS_OCL_4_ITR_W(Arr)		((*Arr).s[3])
#endif

#define WORKGROUP_SIZE 	2
#define PARTICLE_BUFFER	0

#include "../../../kratos/spatial_containers/tree.h"
#include "../custom_utilities/opencl_interface.h"
#include "processes/node_erase_process.h"
#include "includes/model_part.h"
#include <malloc.h>
#include <stdio.h>
#include <string.h>

#include <omp.h>

namespace Kratos {



template<  std::size_t TDimension,
           class TPointType,
		   class TContainerType,
		   class TPointerType = typename TContainerType::value_type,
		   class TIteratorType = typename TContainerType::iterator,
		   class TDistanceIteratorType = typename std::vector<double>::iterator,
		   class TDistanceFunction = Kratos::SearchUtils::SquaredDistanceFunction<TDimension,TPointType> >
class BinsObjectStaticOCL : public TreeNode<TDimension,TPointType, TPointerType, TIteratorType, TDistanceIteratorType> 
{


   public:

      /// Pointer definition of Bins
      KRATOS_CLASS_POINTER_DEFINITION(BinsObjectStaticOCL);

      typedef TPointType                         PointType;
      typedef TContainerType                     ContainerType;
      typedef TIteratorType                      IteratorType;
      typedef TDistanceIteratorType              DistanceIteratorType;
      typedef TPointerType                       PointerType;
      typedef TDistanceFunction                  DistanceFunction;
      enum { Dimension = TDimension };

      typedef TreeNode<Dimension,PointType,PointerType,IteratorType,DistanceIteratorType> TreeNodeType;

      typedef typename TreeNodeType::SizeType        SizeType;
      typedef typename TreeNodeType::IndexType       IndexType;
      typedef typename TreeNodeType::CoordinateType  CoordinateType;

      typedef TreeNodeType LeafType;

      // not always PointVector == ContainerType ( if ContainerType = C array )
      typedef std::vector<PointerType>        PointVector;
      typedef typename PointVector::iterator  PointIterator;

      typedef std::vector<IteratorType>          IteratorVector;
      typedef typename IteratorVector::iterator  IteratorIterator;
      typedef typename IteratorVector::const_iterator IteratorConstIterator;

      typedef Tvector<IndexType,TDimension>   CellType;
    
      typedef typename TreeNodeType::IteratorIteratorType IteratorIteratorType;
      typedef typename TreeNodeType::SearchStructureType SearchStructureType;

      typedef Kratos::SearchUtils::SearchNearestInRange<PointType,PointerType,IteratorType,DistanceFunction,CoordinateType> SearchNearestInRange;
      typedef Kratos::SearchUtils::SearchRadiusInRange<PointType,IteratorType,DistanceIteratorType,DistanceFunction,SizeType,CoordinateType> SearchRadiusInRange;
      typedef Kratos::SearchUtils::SearchBoxInRange<PointType,IteratorType,SizeType,TDimension> SearchBoxInRange;
      
    public:

	//************************************************************************

	//Constructor
	BinsObjectStaticOCL() : mPointBegin(this->NullIterator()), mPointEnd(this->NullIterator()) {};

	//************************************************************************

	BinsObjectStaticOCL( IteratorType const& PointBegin, IteratorType const& PointEnd, OpenCL::DeviceGroup &DeviceGroup, SizeType BucketSize = 1 )
	    : OCLDeviceGroup(DeviceGroup), mPointBegin(PointBegin), mPointEnd(PointEnd)
	{
	    if(mPointBegin==mPointEnd)
		return;
	    
	    mNodeSize = SearchUtils::PointerDistance(mPointBegin, mPointEnd);
	    mTriangleSize = 0;
	    mProblemSize = mNodeSize;
	    mParticles = new cl_double4[mProblemSize];
	    
	    PointsTriangles = (cl_double4 *)malloc(sizeof(cl_double4) * mNodeSize);
	    
	    InitOCL();
	    
	    CalculateBoundingBox();
	    CalculateCellSize();
	    AllocateCellsContainer();
	}

	//************************************************************************
	 
	BinsObjectStaticOCL( IteratorType const& PointBegin, IteratorType const& PointEnd, cl_int4 * const& indexArray, int TSize, int sampleSize, OpenCL::DeviceGroup &DeviceGroup, SizeType BucketSize = 1 )
	    : OCLDeviceGroup(DeviceGroup), mPointBegin(PointBegin), mPointEnd(PointEnd)
	{
	    if(mPointBegin==mPointEnd)
		return;
	    
	    mNodeSize = SearchUtils::PointerDistance(mPointBegin, mPointEnd);
	    mTriangleSize = TSize;
	    
	    PointsTriangles = (cl_double4 *)malloc(sizeof(cl_double4) * mNodeSize); 
	    
	    InitOCL();
	    CalculateBoundingBox(indexArray);
	    CalculateCellSize();
	    AllocateCellsContainer();
	    
	}
	
	BinsObjectStaticOCL( ModelPart * StaticMesh, ModelPart * ParticMesh, OpenCL::DeviceGroup& DeviceGroup) 
	  : OCLDeviceGroup(DeviceGroup) /*, mStaticMesh(StaticMesh), mParticMesh(ParticMesh)*/
	{ 
	    std::cout << "Creadora por Modelpart" << std::endl;
	    
	    mStaticMesh = StaticMesh;
	    mParticMesh = ParticMesh;
	    
	    cl_int4      * mTriangleItr;
	    IteratorType   mNodeItr;
	    cl_double4   * mParticlesItr;
	    cl_double4   * mParticlesVelocityItr;
	    cl_double4   * mParticlesDisplaceItr;
	    cl_double4   * mParticlesForceItr;
	 
	    mNodeSize	  = StaticMesh->NumberOfNodes();
	    mTriangleSize = StaticMesh->NumberOfElements();
	    mProblemSize  = ParticMesh->NumberOfNodes();
	    
	    mPointBegin	    	= new PointType* [mNodeSize];
	    mPointEnd	    	= mPointBegin + mNodeSize;
	    mTriangles 	    	= (cl_int4 *   )malloc(sizeof(cl_int4)    * mTriangleSize);
	    mParticles 		= (cl_double4 *)malloc(sizeof(cl_double4) * (mProblemSize + PARTICLE_BUFFER) );
	    mParticlesVelocity 	= (cl_double4 *)malloc(sizeof(cl_double4) * (mProblemSize + PARTICLE_BUFFER) );
	    mParticlesVelocityOld = (cl_double4 *)malloc(sizeof(cl_double4) * (mProblemSize + PARTICLE_BUFFER) );
	    mParticlesDisplace 	= (cl_double4 *)malloc(sizeof(cl_double4) * (mProblemSize + PARTICLE_BUFFER) );
	    mParticlesForce 	= (cl_double4 *)malloc(sizeof(cl_double4) * (mProblemSize + PARTICLE_BUFFER) );
	    mParticlesPressure 	= (cl_double4 *)malloc(sizeof(cl_double4) * (mProblemSize + PARTICLE_BUFFER) );
	    mParticlesPress_Proj= (cl_double4 *)malloc(sizeof(cl_double4) * (mProblemSize + PARTICLE_BUFFER) );
	    
	    nodesV = (cl_double4 *)malloc(sizeof(cl_double4) * mNodeSize);
 	    nodesF = (cl_double4 *)malloc(sizeof(cl_double4) * mNodeSize);
	    nodesP = (cl_double4 *)malloc(sizeof(cl_double4) * mNodeSize);
 	    nodesR = (cl_double4 *)malloc(sizeof(cl_double4) * mNodeSize);
 	    nodesY = (    double *)malloc(sizeof(double)     * mNodeSize);
 	    nodesT = (    double *)malloc(sizeof(double)     * mNodeSize);
	    FixedV = (	     int *)malloc(sizeof(int) 	     * mNodeSize);
	    
	    std::cout << "\tNodes:     " << mNodeSize << std::endl;
	    std::cout << "\tElements:  " << mTriangleSize << std::endl;
	    std::cout << "\tParticles: " << mProblemSize << std::endl;
	    std::cout << "\tParticle Buffer: " << mProblemSize + PARTICLE_BUFFER << std::endl;
	    
	    //StaticMesh Nodes
	    std::cout << "Cargando Nodos" << std::endl;
	    mNodeItr = mPointBegin;
	    
	    for( ModelPart::NodesContainerType::iterator inode = mStaticMesh->NodesBegin(); inode != mStaticMesh->NodesEnd(); inode++, mNodeItr++) 
	    {
	      PointType auxPoint;

	      auxPoint[0] = inode->X();
	      auxPoint[1] = inode->Y();
	      auxPoint[2] = inode->Z();

	      (*mNodeItr) = new PointType(auxPoint);
	    }
	    
	    //StaticMesh Elements
	    std::cout << "Cargando Elementos" << std::endl;
	    
	    mTriangleItr = mTriangles;
	    for( ModelPart::ElementsContainerType::iterator iel = mStaticMesh->ElementsBegin(); iel != mStaticMesh->ElementsEnd(); iel++, mTriangleItr++) 
	    {
		KRATOS_OCL_4_ITR_X(mTriangleItr) = (int)iel->GetGeometry()[0].Id();
		KRATOS_OCL_4_ITR_Y(mTriangleItr) = (int)iel->GetGeometry()[1].Id();
		KRATOS_OCL_4_ITR_Z(mTriangleItr) = (int)iel->GetGeometry()[2].Id();
		KRATOS_OCL_4_ITR_W(mTriangleItr) = 0;//(int)iel->GetGeometry()[3].Id();
	    }	
	    
	    //Particles
// 	    std::cout << "Cargando Particulas" << std::endl;
// 	    mParticlesItr = mParticles + StaticMesh->NumberOfNodes();
// 	    for( ModelPart::ElementsContainerType::iterator iel = mStaticMesh->ElementsBegin(); iel != mStaticMesh->ElementsEnd(); iel++, mParticlesItr+=4) 
// 	    {
// 		Geometry<Node < 3 > >& geom = iel->GetGeometry();
// 		ComputeGaussPointPositions(geom,mParticlesItr);
// 	    }
	    
	    std::cout << "Filling Nodes" << std::endl;
	    
	    mParticlesItr = mParticles;
	    mParticlesVelocityItr = mParticlesVelocity;
	    mParticlesDisplaceItr = mParticlesDisplace;
	    mParticlesForceItr = mParticlesForce;
	    
	    for( ModelPart::NodesContainerType::iterator inode = mParticMesh->NodesBegin(); inode != mParticMesh->NodesEnd(); inode++, mParticlesItr++, mParticlesVelocityItr++, mParticlesDisplaceItr++, mParticlesForceItr++) 
	    {
	      	Kratos::array_1d<double, 3> & velocity = inode->FastGetSolutionStepValue(VELOCITY,1);
		Kratos::array_1d<double, 3> & displace = inode->FastGetSolutionStepValue(DISPLACEMENT);
	      
		noalias(inode->Coordinates()) = inode->GetInitialPosition();
		
	    	KRATOS_OCL_4_ITR_X(mParticlesItr) = inode->X();
		KRATOS_OCL_4_ITR_Y(mParticlesItr) = inode->Y();
		KRATOS_OCL_4_ITR_Z(mParticlesItr) = inode->Z();
		KRATOS_OCL_4_ITR_W(mParticlesItr) = inode->Id();
		
// 		std::cout << velocity << std::endl;
	
		KRATOS_OCL_4_ITR_X(mParticlesVelocityItr) = velocity[0];
		KRATOS_OCL_4_ITR_Y(mParticlesVelocityItr) = velocity[1];
		KRATOS_OCL_4_ITR_Z(mParticlesVelocityItr) = velocity[2];
		
		KRATOS_OCL_4_ITR_X(mParticlesDisplaceItr) = displace[0];
		KRATOS_OCL_4_ITR_Y(mParticlesDisplaceItr) = displace[1];
		KRATOS_OCL_4_ITR_Z(mParticlesDisplaceItr) = displace[2];
	    }
	    
// 	    mParticMesh->SetNodalSolutionStepVariablesList();
	    
	    std::cout << "Inicializando OCL" << std::endl;
	    InitOCL();
	    std::cout << "Calculando BoundingBox" << std::endl;
	    CalculateBoundingBox();
	    std::cout << "Calculando tamanyo de celdas" << std::endl;
	    CalculateCellSize();
	    std::cout << "Allocatando contenedores de celda" << std::endl;
	    AllocateCellsContainer();
	    std::cout << "Bins inicialiado correctamente" << std::endl;
	    
	    std::cout << std::endl;
	}

	//************************************************************************

	// destructor
	virtual ~BinsObjectStaticOCL(){
	    free(PointsTriangles);
	}

	//************************************************************************

	IteratorType Begin() { return mPointBegin; }

	//************************************************************************

	IteratorType End() { return mPointBegin; }

	//************************************************************************

	CoordinateType CellSize( SizeType const& iDim ) { return mCellSize[iDim]; }

	//************************************************************************

	SizeType NumCell( SizeType const& iDim ) { return mN[iDim]; }

	//************************************************************************
	 
    private:

	//************************************************************************
     
	void InitOCL() 
	{
// 	   OCL_program = OCLDeviceGroup.BuildProgramFromFile("binshashobjects.cl","-D WORKGROUP_SIZE=512");
	   OCL_program = OCLDeviceGroup.BuildProgramFromFile("binshashobjects_tetrahedron.cl","-D WORKGROUP_SIZE=512");

	   OCLGenerateBinsA     = OCLDeviceGroup.RegisterKernel(OCL_program,"GenerateBinsObjectsA");
	   OCLScan_Local1	= OCLDeviceGroup.RegisterKernel(OCL_program,"scanExclusiveLocal1");
	   OCLScan_Local2	= OCLDeviceGroup.RegisterKernel(OCL_program,"scanExclusiveLocal2");
// 	   OCLuniformUpdate	= OCLDeviceGroup.RegisterKernel(OCL_program,"uniformUpdate");
	   OCLGenerateBinsC     = OCLDeviceGroup.RegisterKernel(OCL_program,"GenerateBinsC");
	   
	   OCLUpdate		= OCLDeviceGroup.RegisterKernel(OCL_program,"Update");
	   OCLMove		= OCLDeviceGroup.RegisterKernel(OCL_program,"Move");
	   
	   OCLTransferA		= OCLDeviceGroup.RegisterKernel(OCL_program,"initializeNodes");
	   OCLTransferB 	= OCLDeviceGroup.RegisterKernel(OCL_program,"calculateField");
	   OCLTransferC		= OCLDeviceGroup.RegisterKernel(OCL_program,"updateField");
	   
	   std::cout << "OCL init ok" << std::endl;
	}
	 
	//************************************************************************
	
	void ComputeGaussPointPositions(Geometry< Node < 3 > >& geom, cl_double4 * searchPoint)
	{
	    double one_third = 1.0 / 3.0;
	    double one_sixt = 1.0 / 6.0;
	    double two_third = 2.0 * one_third;

	    //first
	    KRATOS_OCL_4_ARRAY_X(searchPoint,0) = one_sixt * geom[0].X() + one_sixt * geom[1].X() + two_third * geom[2].X();
	    KRATOS_OCL_4_ARRAY_Y(searchPoint,0) = one_sixt * geom[0].Y() + one_sixt * geom[1].Y() + two_third * geom[2].Y();
	    KRATOS_OCL_4_ARRAY_Z(searchPoint,0) = 0;//one_sixt * geom[0].Z() + one_sixt * geom[1].Z() + two_third * geom[2].Z();

	    //second
	    KRATOS_OCL_4_ARRAY_X(searchPoint,1) = two_third * geom[0].X() + one_sixt * geom[1].X() + one_sixt * geom[2].X();
	    KRATOS_OCL_4_ARRAY_Y(searchPoint,1) = two_third * geom[0].Y() + one_sixt * geom[1].Y() + one_sixt * geom[2].Y();
	    KRATOS_OCL_4_ARRAY_Z(searchPoint,1) = 0;//two_third * geom[0].Z() + one_sixt * geom[1].Z() + one_sixt * geom[2].Z();

	    //third
	    KRATOS_OCL_4_ARRAY_X(searchPoint,2) = one_sixt * geom[0].X() + two_third * geom[1].X() + one_sixt * geom[2].X();
	    KRATOS_OCL_4_ARRAY_Y(searchPoint,2) = one_sixt * geom[0].Y() + two_third * geom[1].Y() + one_sixt * geom[2].Y();
	    KRATOS_OCL_4_ARRAY_Z(searchPoint,2) = 0;//one_sixt * geom[0].Z() + two_third * geom[1].Z() + one_sixt * geom[2].Z();

	    //fourth
// 	    KRATOS_OCL_4_ARRAY_X(searchPoint,3) = one_third * geom[0].X() + one_third * geom[1].X() + one_third * geom[2].X();
// 	    KRATOS_OCL_4_ARRAY_Y(searchPoint,3) = one_third * geom[0].Y() + one_third * geom[1].Y() + one_third * geom[2].Y();
// 	    KRATOS_OCL_4_ARRAY_Z(searchPoint,3) = 0;//one_third * geom[0].Z() + one_third * geom[1].Z() + one_third * geom[2].Z();

	}
	
	//************************************************************************
	
	PointType BaryCenter(cl_int4 const& index4)
	{
	    IteratorType Point = mPointBegin;
	    PointType a,b,c,d,r;

 	    a = (**(Point = mPointBegin + KRATOS_OCL_4_X(index4)-1));
 	    b = (**(Point = mPointBegin + KRATOS_OCL_4_Y(index4)-1));
 	    c = (**(Point = mPointBegin + KRATOS_OCL_4_Z(index4)-1));
	    
	    //3D
	    if(KRATOS_OCL_4_W(index4) > 0)
	    {	

		d = (**(Point = mPointBegin + KRATOS_OCL_4_W(index4)-1));
		
		for(SizeType i = 0; i < TDimension; i++)
		{
		  r[i] = (a[i] + b[i] + c[i] + d[i]) / 4; // * 0.25f ? 
		}
	    }
	    //2D
	    else
	    {
		for(SizeType i = 0; i < TDimension; i++)
		{
		    r[i] = (a[i] + b[i] + c[i]) / 3; // * 0.33333333333333334f ? 
		}
	    }
	    return r;
	}

	//Calculate the boundingBox of all objects?
	void CalculateBoundingBox() 
	{
	    for(SizeType i = 0 ; i < TDimension ; i++)
	    {
		mMinPoint[i] = (**mPointBegin)[i];
		mMaxPoint[i] = (**mPointBegin)[i];
	    }
	   
	    for(IteratorType Point = mPointBegin ; Point != mPointEnd ; Point++)
	    {
		for(SizeType i = 0 ; i < TDimension ; i++)
		{
		    if( (**Point)[i] < mMinPoint[i] ) mMinPoint[i] = (**Point)[i];
		    if( (**Point)[i] > mMaxPoint[i] ) mMaxPoint[i] = (**Point)[i];
		}
	    }

	    // Check que se ha echo todo bien
	    std::cout << "Min: ";
	    for(SizeType i = 0 ; i < TDimension ; i++)
		std::cout << mMinPoint[i] << "  ";
	    std::cout << std::endl;
	    
	    std::cout << "Max: ";
	    for(SizeType i = 0 ; i < TDimension ; i++)
		std::cout << mMaxPoint[i] << "  ";
	    std::cout << std::endl;
	}
	 	
	//Calculates the bounding box for barycenters of tringle/tetraedre input
	void CalculateBoundingBox(cl_int4 * const& indexArray) 
	{
	    mTriangles = (cl_int4 *)malloc(sizeof(cl_int4) * mTriangleSize); 
	    
	    for(SizeType i = 0 ; i < TDimension ; i++)
	    {
		mMinPoint[i] = (**mPointBegin)[i];
		mMaxPoint[i] = (**mPointBegin)[i];
	    }
	   
	    for(IteratorType Point = mPointBegin ; Point != mPointEnd ; Point++)
	    {
		for(SizeType i = 0 ; i < TDimension ; i++)
		{
		    if( (**Point)[i] < mMinPoint[i] ) mMinPoint[i] = (**Point)[i];
		    if( (**Point)[i] > mMaxPoint[i] ) mMaxPoint[i] = (**Point)[i];
		}
	    }
	    
	    for(int i = 0; i < mTriangleSize; i++)
	    {
		mTriangles[i] = indexArray[i];
	    }
	   
	    // Check que se ha echo todo bien
	    std::cout << "Min: ";
	    for(SizeType i = 0 ; i < TDimension ; i++)
		std::cout << mMinPoint[i] << "  ";
	    std::cout << std::endl;
	    
	    std::cout << "Max: ";
	    for(SizeType i = 0 ; i < TDimension ; i++)
		std::cout << mMaxPoint[i] << "  ";
	    std::cout << std::endl;
	}

	void CalculateCellSize() 
	{ 
	    CoordinateType delta[TDimension];
	    CoordinateType alpha[TDimension];
	    CoordinateType mult_delta = 1.00;
	    SizeType index = 0;
	  
	    for(SizeType i = 0 ; i < TDimension ; i++) 
	    {
		delta[i] = (mMaxPoint[i] - mMinPoint[i]);
		if ( delta[i] > delta[index] )
		    index = i;
		delta[i] = (delta[i] == 0.00) ? 1.0f : delta[i];
	    }
	    
	    std::cout << "Delta: ";
	    for(SizeType i = 0 ; i < TDimension ; i++)
		std::cout << delta[i] << "  ";
	    std::cout << std::endl;

	    for(SizeType i = 0 ; i < TDimension ; i++){
		alpha[i] = delta[i] / delta[index];
		mult_delta *= alpha[i];
	    }
	    
	    std::cout << "MultDelta: " << mult_delta << std::endl;
	    
	    std::cout << "Alpha: ";
	    for(SizeType i = 0 ; i < TDimension ; i++)
		std::cout << alpha[i] << "  ";
	    std::cout << std::endl;

	    mN[index] = static_cast<SizeType>( pow(static_cast<CoordinateType>((mTriangleSize ? mTriangleSize : mNodeSize)/mult_delta), 1.00/TDimension)+1 );
	   
	    std::cout << "Index: " << index << std::endl;
	    std::cout << "mN[Index]: " << mN[0] << " " << mN[1] << " " << mN[2] << " " << std::endl;
 	    
	    for(SizeType i = 0 ; i < TDimension ; i++)
	    {
		if(i!=index) 
		{
		    mN[i] = static_cast<SizeType>(alpha[i] * mN[index]);
		    mN[i] = ( mN[i] == 0 ) ? 1 : mN[i];
		}
	    }
	    
	    std::cout << "mN[Index]: " << mN[0] << " " << mN[1] << " " << mN[2] << " " << std::endl;
	    
	    for(SizeType i = 0 ; i < TDimension ; i++)
	    {
		mCellSize[i] = delta[i] / mN[i];
		mInvCellSize[i] = 1.00 / mCellSize[i];
	    }
	}

	//************************************************************************
	
	void AllocateCellsContainer() 
	{
	    mCellSizeInt = 1;
	    for(SizeType i = 0 ; i < TDimension ; i++)
		  mCellSizeInt *= mN[i];
	    std::cout << "Size: \t" << mCellSizeInt << std::endl;
	    mIndexCell.resize(mCellSizeInt+1);
	    mIndexCellBegin = mIndexCell.begin();
	    mIndexCellEnd   = mIndexCell.end();
	}

	//************************************************************************
	
	IndexType CalculatePosition( CoordinateType const& ThisCoord, SizeType ThisDimension )
	{
	    CoordinateType d_index = (ThisCoord - mMinPoint[ThisDimension]) * mInvCellSize[ThisDimension];
	    IndexType index = static_cast<SizeType>( (d_index < 0.00) ? 0.00 : d_index );
	    
	    return  (index > mN[ThisDimension]-1) ? mN[ThisDimension]-1 : index;
	}

	//************************************************************************

	IndexType CalculateIndex( PointType const& ThisPoint )
	{
	    IndexType Index = 0;
	    for(SizeType iDim = TDimension-1 ; iDim > 0 ; iDim--)
	    {
		Index += CalculatePosition(ThisPoint[iDim],iDim);
		Index *= mN[iDim-1];
	    }
	    Index += CalculatePosition(ThisPoint[0],0);
	    
	    return Index;
	}

	//************************************************************************

	IndexType CalculateIndex( CellType const& ThisIndex )
	{
	    IndexType Index = 0;
	    for(SizeType iDim = TDimension-1 ; iDim > 0 ; iDim--)
	    {
		Index += ThisIndex[iDim];
		Index *= mN[iDim-1];
	    }
	    Index += ThisIndex[0];
	    
	    return Index;
	}

	//************************************************************************

	CellType CalculateCell( PointType const& ThisPoint )
	{
	    CellType Cell;
	    for(SizeType i = 0 ; i < TDimension ; i++)
		Cell[i] = CalculatePosition(ThisPoint[i],i);
	    
	    return Cell;
	}

	CellType CalculateCell( PointType const& ThisPoint, CoordinateType Radius )
	{
	    CellType Cell;
	    for(SizeType i = 0 ; i < TDimension ; i++)
		Cell[i] = CalculatePosition(ThisPoint[i]+Radius,i);
	    
	    return Cell;
	}

	//************************************************************************

    public:

	//************************************************************************
     
	void GenerateBins()
	{ 
      	    struct timespec begin;
	    struct timespec end;

	    int ArrayLenght = 1;
	    
	    while(ArrayLenght < mCellSizeInt+1)
		ArrayLenght <<= 1;
	    
	    std::cout << "ArrayLenght:\t" << ArrayLenght << std::endl; 
	    
	    int * BinsObjectContainer;
	    cl_double4 * MinPoint      = new cl_double4();
	    
	    int * IndexCellReference   = (int *)malloc(sizeof(int) * ArrayLenght);
	    int * IndexCellReferenceO  = (int *)malloc(sizeof(int) * ArrayLenght);
	    
	    double InvCellSize[TDimension];
	    double N[TDimension];
	        	    
	    clock_gettime( CLOCK_REALTIME, &begin );
	    OCL_IndexCellReference = OCLDeviceGroup.CreateBuffer(sizeof(cl_uint4)   * ArrayLenght / 4   , CL_MEM_READ_WRITE); 
	    OCL_IndexCellReferenceO= OCLDeviceGroup.CreateBuffer(sizeof(cl_uint4)   * ArrayLenght / 4   , CL_MEM_READ_WRITE);
	    OCL_InvCellSize        = OCLDeviceGroup.CreateBuffer(sizeof(double)     * TDimension        , CL_MEM_READ_ONLY );
	    OCL_N                  = OCLDeviceGroup.CreateBuffer(sizeof(double)     * TDimension        , CL_MEM_READ_ONLY );
	    OCL_MinPoint           = OCLDeviceGroup.CreateBuffer(sizeof(cl_double4)                     , CL_MEM_READ_ONLY );
	    OCL_TriangleList   	   = OCLDeviceGroup.CreateBuffer(sizeof(cl_int4)    * mTriangleSize     , CL_MEM_READ_ONLY );
	    OCL_PointsTriangle     = OCLDeviceGroup.CreateBuffer(sizeof(cl_double4) * mNodeSize         , CL_MEM_READ_ONLY );
	    
	    OCL_Scan_Buffer = OCLDeviceGroup.CreateBuffer(sizeof(uint) * ( 64 * ArrayLenght / (4 * OCLDeviceGroup.WorkGroupSizes[OCLScan_Local1][0])), CL_MEM_READ_ONLY );
	 
	    KRATOS_OCL_4_X(MinPoint[0]) = mMinPoint[0];
	    KRATOS_OCL_4_Y(MinPoint[0]) = mMinPoint[1];
	    KRATOS_OCL_4_Z(MinPoint[0]) = mMinPoint[2];
	    
	    for(int i = 0; i < ArrayLenght; i++)
	    {
		IndexCellReference[i] = 0;
		IndexCellReferenceO[i] = 0;
	    }

	    for(SizeType i = 0; i < TDimension; i++) 
	    {
		InvCellSize[i] = mInvCellSize[i];
		N[i] = mN[i];
	    }
	    clock_gettime( CLOCK_REALTIME, &end );
	    
	    std::cout << "Var filling:\t\t" << ((float)(end.tv_sec - begin.tv_sec) + (float)(end.tv_nsec-begin.tv_nsec)/1000000000) << std::endl;
	    
	    /////////////
	    PointsTriangles = (cl_double4 *)malloc(sizeof(cl_double4) * mNodeSize);
	    
	    int k = 0;
	    for( IteratorType Point = mPointBegin ; Point != mPointEnd ; Point++)
	    {
		KRATOS_OCL_4_ARRAY_X(PointsTriangles,k) = (**Point)[0];
		KRATOS_OCL_4_ARRAY_Y(PointsTriangles,k) = (**Point)[1];
		KRATOS_OCL_4_ARRAY_Z(PointsTriangles,k) = (**Point)[2];
		k++;
	    }
	    
	    OCLDeviceGroup.CopyBuffer(OCL_PointsTriangle, OpenCL::HostToDevice, OpenCL::VoidPList(1,PointsTriangles));
	    /////////////
	    
	    OCLDeviceGroup.CopyBuffer(OCL_TriangleList      , OpenCL::HostToDevice, OpenCL::VoidPList(1,mTriangles));
	    OCLDeviceGroup.CopyBuffer(OCL_IndexCellReference, OpenCL::HostToDevice, OpenCL::VoidPList(1,IndexCellReference));
	    OCLDeviceGroup.CopyBuffer(OCL_InvCellSize	    , OpenCL::HostToDevice, OpenCL::VoidPList(1,InvCellSize));
	    OCLDeviceGroup.CopyBuffer(OCL_N                 , OpenCL::HostToDevice, OpenCL::VoidPList(1,N));
	    OCLDeviceGroup.CopyBuffer(OCL_MinPoint          , OpenCL::HostToDevice, OpenCL::VoidPList(1,MinPoint));
	
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLGenerateBinsA, 0,  OCL_PointsTriangle);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLGenerateBinsA, 1,  OCL_TriangleList);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLGenerateBinsA, 2,  OCL_IndexCellReference);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLGenerateBinsA, 3,  OCL_InvCellSize);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLGenerateBinsA, 4,  OCL_N);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLGenerateBinsA, 5,  OCL_MinPoint);
	    OCLDeviceGroup.SetKernelArg(OCLGenerateBinsA, 6, mTriangleSize);
	    
	    std::cout << "Starting Kernel A..." << std::endl;
	    
	    clock_gettime( CLOCK_REALTIME, &begin );
	    OCLDeviceGroup.ExecuteKernel(OCLGenerateBinsA, mTriangleSize);
	    clock_gettime( CLOCK_REALTIME, &end );
	    
	    std::cout << "KERNEL A EXECUTED:\t\t" << ((float)(end.tv_sec - begin.tv_sec) + (float)(end.tv_nsec-begin.tv_nsec)/1000000000) << std::endl;
	    
	    // SCAN PHASE
	    
// 	    OCLDeviceGroup.SetBufferAsKernelArg(OCLScan_Local1, 0, OCL_IndexCellReferenceO);
// 	    OCLDeviceGroup.SetBufferAsKernelArg(OCLScan_Local1, 1, OCL_IndexCellReference);
// 	    OCLDeviceGroup.SetLocalMemAsKernelArg(OCLScan_Local1, 2, 2 * sizeof(uint) * OCLDeviceGroup.WorkGroupSizes[OCLScan_Local1][0]);
// 	    OCLDeviceGroup.SetKernelArg(OCLScan_Local1, 3, (uint)(OCLDeviceGroup.WorkGroupSizes[OCLScan_Local1][0] * 4));
// 	    
// 	    OCLDeviceGroup.SetBufferAsKernelArg(OCLScan_Local2, 0, OCL_Scan_Buffer);
// 	    OCLDeviceGroup.SetBufferAsKernelArg(OCLScan_Local2, 1, OCL_IndexCellReferenceO);
// 	    OCLDeviceGroup.SetBufferAsKernelArg(OCLScan_Local2, 2, OCL_IndexCellReference);
// 	    OCLDeviceGroup.SetLocalMemAsKernelArg(OCLScan_Local2, 3, 2 * sizeof(uint) * OCLDeviceGroup.WorkGroupSizes[OCLScan_Local1][0]);
// 	    OCLDeviceGroup.SetKernelArg(OCLScan_Local2, 4, (uint)(ArrayLenght / (OCLDeviceGroup.WorkGroupSizes[OCLScan_Local1][0] * 4)));
// 	    OCLDeviceGroup.SetKernelArg(OCLScan_Local2, 5, (uint)((OCLDeviceGroup.WorkGroupSizes[OCLScan_Local1][0] * 8 * 4) / (OCLDeviceGroup.WorkGroupSizes[OCLScan_Local1][0] * 4)));
// 	    
// 	    OCLDeviceGroup.SetBufferAsKernelArg(OCLuniformUpdate, 0, OCL_IndexCellReferenceO);
// 	    OCLDeviceGroup.SetBufferAsKernelArg(OCLuniformUpdate, 1, OCL_Scan_Buffer);
// 	    
// 	    std::cout << "Starting Kernel S1..." << std::endl;
// 	    OCLDeviceGroup.ExecuteKernel(OCLScan_Local1, ArrayLenght / 4);
// 	    
// 	    std::cout << "Starting Kernel S2..." << std::endl;
// 	    OCLDeviceGroup.ExecuteKernel(OCLScan_Local2, (ArrayLenght / (OCLDeviceGroup.WorkGroupSizes[OCLScan_Local1][0] * 4)) );
// 	    
// 	    std::cout << "Starting Kernel S3..." << std::endl;
// 	    OCLDeviceGroup.ExecuteKernel(OCLuniformUpdate, (ArrayLenght / 4));

	    OCLDeviceGroup.CopyBuffer(OCL_IndexCellReference, OpenCL::DeviceToHost, OpenCL::VoidPList(1,IndexCellReference));
	    
	    IndexCellReferenceO[0] = 0;
	    for(int i = 1; i < mCellSizeInt+1; i++) {
	      IndexCellReferenceO[i] += IndexCellReferenceO[i-1] + IndexCellReference[i-1];
	    }
	    
	    OCLDeviceGroup.CopyBuffer(OCL_IndexCellReferenceO, OpenCL::HostToDevice, OpenCL::VoidPList(1,IndexCellReferenceO));

	    mBinsObjectSize = IndexCellReferenceO[mCellSizeInt];

	    // END OF SCAN PHASE
	    
	    std::cout << "BinsObjectSize: " << mBinsObjectSize << std:: endl;
	    
	    OCL_BinsObjectContainer  = OCLDeviceGroup.CreateBuffer(sizeof(int) * mBinsObjectSize, CL_MEM_READ_WRITE);
	    BinsObjectContainer = (int *)malloc(sizeof(int) * mBinsObjectSize);
	    
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLGenerateBinsC, 0,  OCL_PointsTriangle);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLGenerateBinsC, 1,  OCL_TriangleList);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLGenerateBinsC, 2,  OCL_IndexCellReferenceO);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLGenerateBinsC, 3,  OCL_InvCellSize);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLGenerateBinsC, 4,  OCL_N);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLGenerateBinsC, 5,  OCL_MinPoint);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLGenerateBinsC, 6,  OCL_BinsObjectContainer);
	    OCLDeviceGroup.SetKernelArg(OCLGenerateBinsC, 7,  mTriangleSize);
	    
	    std::cout << "Starting Kernel C..." << std::endl;
	    
	    clock_gettime( CLOCK_REALTIME, &begin );
	    OCLDeviceGroup.ExecuteKernel(OCLGenerateBinsC, mTriangleSize);
	    clock_gettime( CLOCK_REALTIME, &end );
	    
	    std::cout << "KernelC Ex:\t\t" << ((float)(end.tv_sec - begin.tv_sec) + (float)(end.tv_nsec-begin.tv_nsec)/1000000000) << std::endl;
// 	    abort();
	    clock_gettime( CLOCK_REALTIME, &begin );
	    OCLDeviceGroup.CopyBuffer(OCL_BinsObjectContainer, OpenCL::DeviceToHost, OpenCL::VoidPList(1,BinsObjectContainer));
	    clock_gettime( CLOCK_REALTIME, &end );
	    
	    std::cout << "Kernel Copy:\t\t" << ((float)(end.tv_sec - begin.tv_sec) + (float)(end.tv_nsec-begin.tv_nsec)/1000000000) << std::endl;
	    
	    OCLDeviceGroup.CopyBuffer(OCL_IndexCellReferenceO, OpenCL::DeviceToHost, OpenCL::VoidPList(1,IndexCellReferenceO));
	    
	}

	//************************************************************************
	
	void LoadSample(/*cl_double4 * sample,
			cl_double4 * emitters,
			int sampleSize*/)
	{
// 	    mProblemSize = sampleSize;
// 
// 	    mParticles = sample;
// 	    mPointsEmitter  = emitters;

	    cl_double4   * mParticlesItr;
	    
	    //Velocity Field
	    std::cout << "Cargando campo de veclodiades" << std::endl;
	    mParticlesItr = mParticles;
	    
	    CopyStaticmeshData(1);
	    TransferStaticMeshToGPU();
	}
	
	void CopyStaticmeshData(const int direction)
	{
	    if (direction) 
	    {
		int i = 0;
		for( ModelPart::NodesContainerType::iterator inode = mStaticMesh->NodesBegin(); inode != mStaticMesh->NodesEnd(); inode++) 
		{	
		    array_1d<double, 3 > & velocity   = inode->FastGetSolutionStepValue(VELOCITY,1);
		    array_1d<double, 3 > & force      = inode->FastGetSolutionStepValue(FORCE);
		    array_1d<double, 3 > & press_proj = inode->FastGetSolutionStepValue(PRESS_PROJ);
		    
		    KRATOS_OCL_4_ARRAY_X(nodesV,i) = velocity[0];
		    KRATOS_OCL_4_ARRAY_Y(nodesV,i) = velocity[1];
		    KRATOS_OCL_4_ARRAY_Z(nodesV,i) = velocity[2];
		    KRATOS_OCL_4_ARRAY_W(nodesV,i) = inode->Id();
		    
		    KRATOS_OCL_4_ARRAY_X(nodesF,i) = force[0];
		    KRATOS_OCL_4_ARRAY_Y(nodesF,i) = force[1];
		    KRATOS_OCL_4_ARRAY_Z(nodesF,i) = force[2];
		    
		    KRATOS_OCL_4_ARRAY_X(nodesP,i) = press_proj[0];
		    KRATOS_OCL_4_ARRAY_Y(nodesP,i) = press_proj[1];
		    KRATOS_OCL_4_ARRAY_Z(nodesP,i) = press_proj[2];
		    
		    FixedV[i] = inode->IsFixed(VELOCITY_X);
		    
		    i++;
		}
	    } else {
	      
		for(int i = 0; i < mNodeSize; i++) 
		{
// 		  std::cout << KRATOS_OCL_4_ARRAY_W(nodesV,i) << std::endl;
		    Node<3>::Pointer inode = mStaticMesh->Nodes()(KRATOS_OCL_4_ARRAY_W(nodesV,i));

		    array_1d<double, 3 > & velocity   = inode->FastGetSolutionStepValue(VELOCITY);
// 		    array_1d<double, 3 > & force      = (inode)->FastGetSolutionStepValue(FORCE);
// 		    array_1d<double, 3 > & press_proj = (inode)->FastGetSolutionStepValue(PRESS_PROJ);
		    
		    velocity[0] = KRATOS_OCL_4_ARRAY_X(nodesV,i);
		    velocity[1] = KRATOS_OCL_4_ARRAY_Y(nodesV,i);
		    velocity[2] = KRATOS_OCL_4_ARRAY_Z(nodesV,i);
		    
// 		    force[0] = KRATOS_OCL_4_ARRAY_X(nodesF,i);
// 		    force[1] = KRATOS_OCL_4_ARRAY_Y(nodesF,i);
// 		    force[2] = KRATOS_OCL_4_ARRAY_Z(nodesF,i);
// 		    
// 		    press_proj[0] = KRATOS_OCL_4_ARRAY_X(nodesP,i);
// 		    press_proj[1] = KRATOS_OCL_4_ARRAY_Y(nodesP,i);
// 		    press_proj[2] = KRATOS_OCL_4_ARRAY_Z(nodesP,i);
		    
// 		    i++;
		}
	    }
	}
	
	void TransferStaticMeshToGPU()
	{
	    OCLDeviceGroup.CopyBuffer(OCL_NodesV, OpenCL::HostToDevice, OpenCL::VoidPList(1,&nodesV[0]));
	    OCLDeviceGroup.CopyBuffer(OCL_NodesF, OpenCL::HostToDevice, OpenCL::VoidPList(1,&nodesF[0]));
	    OCLDeviceGroup.CopyBuffer(OCL_NodesP, OpenCL::HostToDevice, OpenCL::VoidPList(1,&nodesP[0]));
	    OCLDeviceGroup.CopyBuffer(OCL_FixedV, OpenCL::HostToDevice, OpenCL::VoidPList(1,&FixedV[0]));
	}
	
	void TransferStaticMeshToCPU()
	{
	    OCLDeviceGroup.CopyBuffer(OCL_NodesV, OpenCL::DeviceToHost, OpenCL::VoidPList(1,&nodesV[0]));
// 	    OCLDeviceGroup.CopyBuffer(OCL_NodesF, OpenCL::DeviceToHost, OpenCL::VoidPList(1,&nodesF[0]));
// 	    OCLDeviceGroup.CopyBuffer(OCL_NodesP, OpenCL::DeviceToHost, OpenCL::VoidPList(1,&nodesP[0]));
	}
	
	void TransferParticMeshToGPU() {
	    cl_double4 * mParticlesItr 	       = mParticles;
	    cl_double4 * mParticlesVelocityItr = mParticlesVelocity;
	    cl_double4 * mParticlesDisplaceItr = mParticlesDisplace;
	    
	    mProblemSize = mParticMesh->NumberOfNodes();
	  
	    for( ModelPart::NodesContainerType::iterator inode = mParticMesh->NodesBegin(); inode != mParticMesh->NodesEnd(); inode++, mParticlesItr++, mParticlesVelocityItr++, mParticlesDisplaceItr++) 
	    {
	      	Kratos::array_1d<double, 3> & velocity = inode->FastGetSolutionStepValue(VELOCITY,1);
		Kratos::array_1d<double, 3> & displace = inode->FastGetSolutionStepValue(DISPLACEMENT,1);
		
// 		std::cout << velocity << std::endl;		
		noalias(inode->Coordinates()) = inode->GetInitialPosition();
		
	    	KRATOS_OCL_4_ITR_X(mParticlesItr) = inode->X();
		KRATOS_OCL_4_ITR_Y(mParticlesItr) = inode->Y();
		KRATOS_OCL_4_ITR_Z(mParticlesItr) = inode->Z();
		KRATOS_OCL_4_ITR_W(mParticlesItr) = inode->Id();
	
		KRATOS_OCL_4_ITR_X(mParticlesVelocityItr) = velocity[0];
		KRATOS_OCL_4_ITR_Y(mParticlesVelocityItr) = velocity[1];
		KRATOS_OCL_4_ITR_Z(mParticlesVelocityItr) = velocity[2];
		
		KRATOS_OCL_4_ITR_X(mParticlesDisplaceItr) = displace[0];
		KRATOS_OCL_4_ITR_Y(mParticlesDisplaceItr) = displace[1];
		KRATOS_OCL_4_ITR_Z(mParticlesDisplaceItr) = displace[2];
		
		inode->GetValue(ERASE_FLAG) = true;
	    }
	    
	    OCLDeviceGroup.CopyBuffer(OCL_Particles, OpenCL::HostToDevice, OpenCL::VoidPList(1,&mParticles[0]));
	    OCLDeviceGroup.CopyBuffer(OCL_ParticlesVelocity, OpenCL::HostToDevice, OpenCL::VoidPList(1,&mParticlesVelocity[0]));
	    OCLDeviceGroup.CopyBuffer(OCL_ParticlesDisplace, OpenCL::HostToDevice, OpenCL::VoidPList(1,&mParticlesDisplace[0]));
	}
	
	//************************************************************************
	  
	void AllocateOCLBuffers(int ConcurrentPoints/*, GLuint glBuffer*/)
	{
	    std::cout << "-Allocating Buffers-" << std::endl;
	  
	    int ConcurrentPointsReal = ConcurrentPoints;
	    
	    ConcurrentPointsReal = ((ConcurrentPointsReal + OCLDeviceGroup.WorkGroupSizes[OCLMove][0] - 1) / OCLDeviceGroup.WorkGroupSizes[OCLMove][0]) * OCLDeviceGroup.WorkGroupSizes[OCLMove][0];
	    
	    if(ConcurrentPoints > mProblemSize) 
	      ConcurrentPointsReal = mProblemSize;

	    std::cout << "\tAllocate Bufers: " << mProblemSize << std::endl;
	    std::cout << "\tAllocate Bufers Real: " << ConcurrentPointsReal << std::endl;
	    
	    ///////////////////////////////////////////////////////////////////////////
	    ConcurrentPointsReal = mProblemSize + PARTICLE_BUFFER;
	    ///////////////////////////////////////////////////////////////////////////
	    
	    OCL_Radius        = OCLDeviceGroup.CreateBuffer(sizeof(double) * 2 , CL_MEM_READ_ONLY);
	    OCL_Active	      = OCLDeviceGroup.CreateBuffer(sizeof(int) * 2    , CL_MEM_READ_WRITE);
	    OCL_Body_Force    = OCLDeviceGroup.CreateBuffer(sizeof(cl_double4) , CL_MEM_READ_WRITE);

	    OCL_Particles             = OCLDeviceGroup.CreateBuffer(sizeof(cl_double4) * ConcurrentPointsReal, CL_MEM_READ_WRITE);
	    OCL_ParticlesVelocity     = OCLDeviceGroup.CreateBuffer(sizeof(cl_double4) * ConcurrentPointsReal, CL_MEM_READ_WRITE);
	    OCL_ParticlesVelocityOld  = OCLDeviceGroup.CreateBuffer(sizeof(cl_double4) * ConcurrentPointsReal, CL_MEM_READ_WRITE);
	    OCL_ParticlesDisplace     = OCLDeviceGroup.CreateBuffer(sizeof(cl_double4) * ConcurrentPointsReal, CL_MEM_READ_WRITE);
	    OCL_ParticlesAcceleration = OCLDeviceGroup.CreateBuffer(sizeof(cl_double4) * ConcurrentPointsReal, CL_MEM_READ_WRITE);
	    OCL_ParticlesForce        = OCLDeviceGroup.CreateBuffer(sizeof(cl_double4) * ConcurrentPointsReal, CL_MEM_READ_WRITE);
	    OCL_ParticlesPress_Proj   = OCLDeviceGroup.CreateBuffer(sizeof(cl_double4) * ConcurrentPointsReal, CL_MEM_READ_WRITE);
	    OCL_ParticlesLock	      = OCLDeviceGroup.CreateBuffer(sizeof(int)        * ConcurrentPointsReal, CL_MEM_READ_WRITE);

            OCL_NodesV        = OCLDeviceGroup.CreateBuffer(sizeof(cl_double4) * mNodeSize, CL_MEM_READ_WRITE);
	    OCL_NodesF        = OCLDeviceGroup.CreateBuffer(sizeof(cl_double4) * mNodeSize, CL_MEM_READ_WRITE);
	    OCL_NodesP        = OCLDeviceGroup.CreateBuffer(sizeof(cl_double4) * mNodeSize, CL_MEM_READ_WRITE);
	    OCL_NodesR        = OCLDeviceGroup.CreateBuffer(sizeof(cl_double4) * mNodeSize, CL_MEM_READ_WRITE);
	    OCL_NodesY        = OCLDeviceGroup.CreateBuffer(sizeof(double)     * mNodeSize, CL_MEM_READ_WRITE);
	    OCL_NodesT        = OCLDeviceGroup.CreateBuffer(sizeof(double)     * mNodeSize, CL_MEM_READ_WRITE);
	    OCL_FixedV        = OCLDeviceGroup.CreateBuffer(sizeof(int)        * mNodeSize, CL_MEM_READ_WRITE);
	    OCL_FLAG          = OCLDeviceGroup.CreateBuffer(sizeof(int)        * mNodeSize, CL_MEM_READ_WRITE);
	    OCL_results       = OCLDeviceGroup.CreateBuffer(sizeof(int)        * ConcurrentPointsReal, CL_MEM_WRITE_ONLY);
	    OCL_outData       = OCLDeviceGroup.CreateBuffer(sizeof(int)        * ConcurrentPointsReal, CL_MEM_WRITE_ONLY);
	    OCL_distance      = OCLDeviceGroup.CreateBuffer(sizeof(double)     * ConcurrentPointsReal, CL_MEM_READ_WRITE);
	    OCL_NFunction     = OCLDeviceGroup.CreateBuffer(sizeof(cl_double4) * ConcurrentPointsReal, CL_MEM_READ_WRITE);
	    OCL_Dens	      = OCLDeviceGroup.CreateBuffer(sizeof(int)        * mTriangleSize, CL_MEM_READ_WRITE);
	}
	
	void InitializeBuffers(double Radius)
	{
	    std::cout << "-Initialize Buffers-" << std::endl;
	  
	    double HOST_memRadius[2];
	    
	    HOST_memRadius[0] = Radius;
	    HOST_memRadius[1] = Radius * Radius;
	    
	    int k = 0;
	    for( IteratorType Point = mPointBegin ; Point != mPointEnd ; Point++)
	    {
		KRATOS_OCL_4_ARRAY_X(PointsTriangles,k) = (**Point)[0];
		KRATOS_OCL_4_ARRAY_Y(PointsTriangles,k) = (**Point)[1];
		KRATOS_OCL_4_ARRAY_Z(PointsTriangles,k) = (**Point)[2];
		k++;
	    }
	    
	    OCLDeviceGroup.CopyBuffer(OCL_Radius        , OpenCL::HostToDevice, OpenCL::VoidPList(1,&HOST_memRadius));
	    OCLDeviceGroup.CopyBuffer(OCL_PointsTriangle, OpenCL::HostToDevice, OpenCL::VoidPList(1,PointsTriangles));
	    
	    int * dens = (int *)malloc(sizeof(int) * mTriangleSize);
	    int * particleLock = (int *)malloc(sizeof(int) * (mParticMesh->NumberOfNodes() + PARTICLE_BUFFER));
	    int * triangleFLAG = (int *)malloc(sizeof(int) * (mNodeSize));
	    
	    for(int i = 0 ; i < mTriangleSize; i++) {
	      dens[i] = 0 ;
	    }
	    
	    for(int i = 0 ; i < (mParticMesh->NumberOfNodes() + PARTICLE_BUFFER); i++) {
	      particleLock[i] = 1;
	    }
	    
	    for(int i = 0; i < mNodeSize; i++) {
	      triangleFLAG[i] = 0;
	    }
	    
	    OCLDeviceGroup.CopyBuffer(OCL_Dens, OpenCL::HostToDevice, OpenCL::VoidPList(1,dens));
	    
	    OCLDeviceGroup.CopyBuffer(OCL_Particles, OpenCL::HostToDevice, OpenCL::VoidPList(1,&mParticles[0]));
	    OCLDeviceGroup.CopyBuffer(OCL_ParticlesDisplace, OpenCL::HostToDevice, OpenCL::VoidPList(1,&mParticlesDisplace[0]));
	    OCLDeviceGroup.CopyBuffer(OCL_ParticlesVelocity, OpenCL::HostToDevice, OpenCL::VoidPList(1,&mParticlesVelocity[0]));
	    OCLDeviceGroup.CopyBuffer(OCL_ParticlesVelocityOld, OpenCL::HostToDevice, OpenCL::VoidPList(1,&mParticlesVelocityOld[0]));
	    OCLDeviceGroup.CopyBuffer(OCL_ParticlesForce, OpenCL::HostToDevice, OpenCL::VoidPList(1,&mParticlesForce[0]));
	    OCLDeviceGroup.CopyBuffer(OCL_ParticlesPress_Proj, OpenCL::HostToDevice, OpenCL::VoidPList(1,&mParticlesPress_Proj[0]));
	    OCLDeviceGroup.CopyBuffer(OCL_ParticlesLock, OpenCL::HostToDevice, OpenCL::VoidPList(1,&particleLock[0]));
	    OCLDeviceGroup.CopyBuffer(OCL_FLAG, OpenCL::HostToDevice, OpenCL::VoidPList(1,&triangleFLAG[0]));
	    
	}
         
	void SearchInRadiusOCL(double Radius, uint ConcurrentPoints, uint maxResults) 
	{
	    uint processed = 0;
	    
	    int result = 0;
	    
	    double HOST_memRadius  = Radius;
	    double HOST_memRadius2 = Radius * Radius;
	    
	    uint amount;
	
	    while (processed < mProblemSize)
	    {	  
		amount = (ConcurrentPoints > mProblemSize) ? mProblemSize : (processed + ConcurrentPoints) < mProblemSize ? ConcurrentPoints : mProblemSize - processed;
		
		int * results;
		int * resultPoints;
		
		results      = (int *)malloc(sizeof(int) * amount);
		resultPoints = (int *)malloc(sizeof(int) * amount * maxResults);

		OCLDeviceGroup.CopyBuffer(OCL_Radius   , OpenCL::HostToDevice, OpenCL::VoidPList(1,&HOST_memRadius));
		OCLDeviceGroup.CopyBuffer(OCL_Radius2  , OpenCL::HostToDevice, OpenCL::VoidPList(1,&HOST_memRadius2));
		OCLDeviceGroup.CopyBuffer(OCL_Particles, OpenCL::HostToDevice, OpenCL::VoidPList(1,&mParticles[processed]));
				
		OCLDeviceGroup.SetBufferAsKernelArg(OCLSearchInRadius, 0,  OCL_IndexCellReferenceO);
		OCLDeviceGroup.SetBufferAsKernelArg(OCLSearchInRadius, 1,  OCL_BinsObjectContainer);
		OCLDeviceGroup.SetBufferAsKernelArg(OCLSearchInRadius, 2 , OCL_InvCellSize);
		OCLDeviceGroup.SetBufferAsKernelArg(OCLSearchInRadius, 3 , OCL_N);
		OCLDeviceGroup.SetBufferAsKernelArg(OCLSearchInRadius, 4,  OCL_Radius);
		OCLDeviceGroup.SetBufferAsKernelArg(OCLSearchInRadius, 5,  OCL_Radius2);
		OCLDeviceGroup.SetBufferAsKernelArg(OCLSearchInRadius, 6,  OCL_Particles);
		OCLDeviceGroup.SetBufferAsKernelArg(OCLSearchInRadius, 7,  OCL_MinPoint);
		OCLDeviceGroup.SetBufferAsKernelArg(OCLSearchInRadius, 8,  OCL_outData);
		OCLDeviceGroup.SetBufferAsKernelArg(OCLSearchInRadius, 9,  OCL_results);
		OCLDeviceGroup.SetKernelArg(OCLSearchInRadius, 10, maxResults);

		OCLDeviceGroup.ExecuteKernel(OCLSearchInRadius, amount);
		
		OCLDeviceGroup.CopyBuffer(OCL_results, OpenCL::DeviceToHost, OpenCL::VoidPList(1,results));
		OCLDeviceGroup.CopyBuffer(OCL_outData, OpenCL::DeviceToHost, OpenCL::VoidPList(1,resultPoints));
		
		for(size_t j = 0; j < amount; j++)
		    result += results[j];
		
		free(results);
		free(resultPoints);
		
		processed += amount;
		//std::cout << "Total Results while " << processed << " points done: " << result << std::endl;
	    }
	    
	    std::cout << "Total Results: " << result << " (" << maxResults * SearchUtils::PointerDistance(mPointBegin, mPointEnd) << " stored)"<< std::endl;
	}
	
	void SearchTriangles(array_1d<double, 3 > & body_force, const double density, const double dt, const double substeps, const int ConcurrentPoints, const int use_eulerian,int copy_data) 
	{  
	    std::cout << "-Search Triangles-" << std::endl;
	  
	    int amount = 0;
	    int processed = 0;
	    int totalErrors = 0;
	    static int time = 0;
	    int activeParticles[2];
	    
	    cl_double4 * mParticlesItr;
	    cl_double4 * mParticlesVelocityItr;
	    cl_double4 * mParticlesDisplaceItr;
	    cl_double4 * mParticlesForceItr;
	    cl_double4 * mParticlesPress_ProjItr;
	    
	    struct timespec begin;
	    struct timespec end;
	    
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLUpdate, 0, OCL_Particles);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLUpdate, 1, OCL_PointsTriangle);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLUpdate, 2, OCL_ParticlesDisplace);
            OCLDeviceGroup.SetBufferAsKernelArg(OCLUpdate, 3, OCL_TriangleList);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLUpdate, 4, OCL_BinsObjectContainer);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLUpdate, 5, OCL_IndexCellReferenceO);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLUpdate, 6, OCL_Dens);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLUpdate, 7, OCL_ParticlesLock);

            OCLDeviceGroup.SetBufferAsKernelArg(OCLMove, 0,  OCL_IndexCellReferenceO);
            OCLDeviceGroup.SetBufferAsKernelArg(OCLMove, 1,  OCL_BinsObjectContainer);
            OCLDeviceGroup.SetBufferAsKernelArg(OCLMove, 2,  OCL_PointsTriangle);
            OCLDeviceGroup.SetBufferAsKernelArg(OCLMove, 3,  OCL_TriangleList);
            OCLDeviceGroup.SetBufferAsKernelArg(OCLMove, 4,  OCL_InvCellSize);
            OCLDeviceGroup.SetBufferAsKernelArg(OCLMove, 5,  OCL_N);
            OCLDeviceGroup.SetBufferAsKernelArg(OCLMove, 6,  OCL_Radius);
            OCLDeviceGroup.SetBufferAsKernelArg(OCLMove, 7,  OCL_MinPoint);
            OCLDeviceGroup.SetBufferAsKernelArg(OCLMove, 8,  OCL_NFunction);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLMove, 9,  OCL_Particles);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLMove, 10, OCL_ParticlesVelocity);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLMove, 11, OCL_ParticlesAcceleration);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLMove, 12, OCL_ParticlesDisplace);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLMove, 13, OCL_ParticlesForce);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLMove, 14, OCL_ParticlesPress_Proj);
 	    OCLDeviceGroup.SetBufferAsKernelArg(OCLMove, 15, OCL_NodesV);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLMove, 16, OCL_NodesF);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLMove, 17, OCL_NodesP);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLMove, 18, OCL_Dens);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLMove, 19, OCL_Body_Force);
	    OCLDeviceGroup.SetKernelArg(	OCLMove, 20, density);
	    OCLDeviceGroup.SetKernelArg(	OCLMove, 21, dt);
	    OCLDeviceGroup.SetKernelArg(	OCLMove, 22, substeps);
	    OCLDeviceGroup.SetKernelArg(	OCLMove, 23, use_eulerian);
	    
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLMove, 25, OCL_ParticlesVelocityOld);
	    OCLDeviceGroup.SetLocalMemAsKernelArg(OCLMove, 26, OCLDeviceGroup.WorkGroupSizes[OCLMove][0]*sizeof(cl_int4));
	    
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLTransferA, 0,  OCL_NodesV);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLTransferA, 1,  OCL_FixedV);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLTransferA, 2,  OCL_NodesY);
	    
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLTransferB, 0,  OCL_IndexCellReferenceO);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLTransferB, 1,  OCL_BinsObjectContainer);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLTransferB, 2,  OCL_PointsTriangle);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLTransferB, 3,  OCL_TriangleList);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLTransferB, 4 , OCL_InvCellSize);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLTransferB, 5 , OCL_N);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLTransferB, 6,  OCL_Radius);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLTransferB, 7,  OCL_MinPoint);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLTransferB, 8,  OCL_Particles);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLTransferB, 9,  OCL_ParticlesVelocity);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLTransferB, 10, OCL_NodesV);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLTransferB, 11, OCL_FixedV);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLTransferB, 12, OCL_NodesY);
	    OCLDeviceGroup.SetLocalMemAsKernelArg(OCLTransferB, 13, sizeof(int) * 256);
	    
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLTransferC, 0,  OCL_NodesV);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLTransferC, 1,  OCL_FixedV);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLTransferC, 2,  OCL_NodesY);
	    
	    clock_gettime( CLOCK_REALTIME, &begin );

	    //Parece que el modelpart de particulas que no pasan puede variar, asi que por ahora provamos a incrementar el tamanyo
	    
	    while (processed < mProblemSize)
	    {	 
		amount = (ConcurrentPoints > mProblemSize) ? mProblemSize : (processed + ConcurrentPoints) < mProblemSize ? ConcurrentPoints : mProblemSize - processed;
		
		amount = mProblemSize;
		activeParticles[0] = 0;
		activeParticles[1] = amount;
   
		OCLDeviceGroup.CopyBuffer(OCL_Active    , OpenCL::HostToDevice, OpenCL::VoidPList(1,&activeParticles));
		OCLDeviceGroup.CopyBuffer(OCL_Body_Force, OpenCL::HostToDevice, OpenCL::VoidPList(1,&body_force));
		
                OCLDeviceGroup.SetKernelArg(OCLMove, 24, amount);
		
		OCLDeviceGroup.SetKernelArg(OCLTransferB, 14, amount);
		
		OCLDeviceGroup.SetBufferAsKernelArg(OCLUpdate, 8, OCL_Active);
		OCLDeviceGroup.SetKernelArg(OCLUpdate, 9, mTriangleSize);
		
// 		OCLDeviceGroup.ExecuteKernel(OCLMove, (int)ceil(sqrt(amount)));
		std::cout << OCLDeviceGroup.WorkGroupSizes[OCLMove][0] << std::endl;
// 		OCLDeviceGroup.ExecuteKernel(OCLMove, OCLDeviceGroup.WorkGroupSizes[OCLMove][0]);
		OCLDeviceGroup.ExecuteKernel(OCLMove, amount);
		
// 		OCLDeviceGroup.ExecuteKernel(OCLUpdate, mTriangleSize);
		
		OCLDeviceGroup.CopyBuffer(OCL_Particles        , OpenCL::DeviceToHost, OpenCL::VoidPList(1,&mParticles[processed]));
		OCLDeviceGroup.CopyBuffer(OCL_ParticlesVelocity, OpenCL::DeviceToHost, OpenCL::VoidPList(1,&mParticlesVelocity[processed]));
		OCLDeviceGroup.CopyBuffer(OCL_ParticlesVelocityOld, OpenCL::DeviceToHost, OpenCL::VoidPList(1,&mParticlesVelocityOld[processed]));
		OCLDeviceGroup.CopyBuffer(OCL_ParticlesDisplace, OpenCL::DeviceToHost, OpenCL::VoidPList(1,&mParticlesDisplace[processed]));
		OCLDeviceGroup.CopyBuffer(OCL_ParticlesForce   , OpenCL::DeviceToHost, OpenCL::VoidPList(1,&mParticlesForce[processed]));

		OCLDeviceGroup.ExecuteKernel(OCLTransferA, mNodeSize);
		OCLDeviceGroup.ExecuteKernel(OCLTransferB, (int)ceil(sqrt(amount)));
		OCLDeviceGroup.ExecuteKernel(OCLTransferC, mNodeSize);
		
 		int errors = 0;
		double diff = 0;
		
		//Error Check End
		
		processed += amount;
		totalErrors += errors;
		
// 		ConcurrentPoints = activeParticles;
	    }
	    
	    clock_gettime( CLOCK_REALTIME, &end );
	    
	    std::cout << "Total Processed: " << processed << " Results: \033[41m" << processed - totalErrors << "\033[0;0m  ";// <<  std::endl;
	    std::cout << "\t\t\t" << ((float)(end.tv_sec - begin.tv_sec) + (float)(end.tv_nsec-begin.tv_nsec)/1000000000) << std::endl; /* MAIN */
    
    	    //Write back particles to ModelPart
	    std::cout << "Writing back Particles" << std::endl;
	    
	    mParticlesItr           = mParticles;
	    mParticlesVelocityItr   = mParticlesVelocity;
	    mParticlesDisplaceItr   = mParticlesDisplace;
	    mParticlesForceItr      = mParticlesForce;
	    mParticlesPress_ProjItr = mParticlesPress_Proj;
	    
	    //Always copy
// 	    for(int i = 0; i < mParticMesh->NumberOfNodes(); i++) 
// 	    {
// 		if (KRATOS_OCL_4_ARRAY_W(mParticles,i) >= 0) {
// 	      
// 		  Node<3>::Pointer inode = mParticMesh->Nodes()(KRATOS_OCL_4_ARRAY_W(mParticles,i));
// 		  
// 		  Kratos::array_1d<double, 3> & velocity = inode->FastGetSolutionStepValue(VELOCITY,1);
// 		  
// 		  velocity[0] = KRATOS_OCL_4_ARRAY_X(mParticlesVelocityOld,i);
// 		  velocity[1] = KRATOS_OCL_4_ARRAY_Y(mParticlesVelocityOld,i);
// 		  velocity[2] = KRATOS_OCL_4_ARRAY_Z(mParticlesVelocityOld,i);
// 		  
// 		  inode->GetValue(ERASE_FLAG) = false;
// 		
// 		} /*else {
// 		  
// 		  Node<3>::Pointer inode = mParticMesh->Nodes()(-KRATOS_OCL_4_ARRAY_W(mParticles,i));
// 		  inode->GetValue(ERASE_FLAG) = true;
// 		  
// 		}*/
// 	    }
	    
	    //Only copy in last transfer
	    if(copy_data || true)
	    for(int i = 0; i < mParticMesh->NumberOfNodes(); i++) 
	    {
// 		if (KRATOS_OCL_4_ARRAY_W(mParticles,i) >= 0) {
	      
		  Node<3>::Pointer inode = (KRATOS_OCL_4_ARRAY_W(mParticles,i) >= 0) ? mParticMesh->Nodes()(KRATOS_OCL_4_ARRAY_W(mParticles,i)) : mParticMesh->Nodes()(-KRATOS_OCL_4_ARRAY_W(mParticles,i));
		
		  Kratos::array_1d<double, 3> & velocity = inode->FastGetSolutionStepValue(VELOCITY);
		  Kratos::array_1d<double, 3> & displace = inode->FastGetSolutionStepValue(DISPLACEMENT);
		  Kratos::array_1d<double, 3> & force    = inode->FastGetSolutionStepValue(FORCE);
		  
		  inode->X() = KRATOS_OCL_4_ARRAY_X(mParticles,i);
		  inode->Y() = KRATOS_OCL_4_ARRAY_Y(mParticles,i);
		  inode->Z() = KRATOS_OCL_4_ARRAY_Z(mParticles,i);
		  
		  velocity[0] = KRATOS_OCL_4_ARRAY_X(mParticlesVelocity,i);
		  velocity[1] = KRATOS_OCL_4_ARRAY_Y(mParticlesVelocity,i);
		  velocity[2] = KRATOS_OCL_4_ARRAY_Z(mParticlesVelocity,i);
		  
		  displace[0] = KRATOS_OCL_4_ARRAY_X(mParticlesDisplace,i);
		  displace[1] = KRATOS_OCL_4_ARRAY_Y(mParticlesDisplace,i);
		  displace[2] = KRATOS_OCL_4_ARRAY_Z(mParticlesDisplace,i);
		  
		  force[0] = KRATOS_OCL_4_ARRAY_X(mParticlesForce,i);
		  force[1] = KRATOS_OCL_4_ARRAY_Y(mParticlesForce,i);
		  force[2] = KRATOS_OCL_4_ARRAY_Z(mParticlesForce,i);
		  
		  inode->GetValue(ERASE_FLAG) = false;
		
// 		} else {
// 		  
// 		  Node<3>::Pointer inode = mParticMesh->Nodes()(-KRATOS_OCL_4_ARRAY_W(mParticles,i));
// 		  inode->GetValue(ERASE_FLAG) = true;
// 		  
// 		}
	    }
	    
// 	    NodeEraseProcess(*mParticMesh).Execute();
	    
// 	    mProblemSize = activeParticles[1] > (PARTICLE_BUFFER + mStaticMesh->NumberOfNodes() + mStaticMesh->NumberOfElements() * 4) ? (PARTICLE_BUFFER + mStaticMesh->NumberOfNodes() + mStaticMesh->NumberOfElements() * 4) : activeParticles[1];
	}
         
	void SearchNearestOCL(double Radius, uint ConcurrentPoints) 
	{
	    uint processed = 0;
	    
	    int result = 0;
	    
	    double HOST_memRadius  = Radius;
	    double HOST_memRadius2 = Radius * Radius;
	    
	    uint amount;
	
	    while (processed < mProblemSize)
	    {	  
		amount = (ConcurrentPoints > mProblemSize) ? mProblemSize : ((processed + ConcurrentPoints) < mProblemSize) ? ConcurrentPoints : mProblemSize - processed;
		
		int * resultPoints;
		
		resultPoints = (int *)malloc(sizeof(int) * amount);

		OCLDeviceGroup.CopyBuffer(OCL_Radius   , OpenCL::HostToDevice, OpenCL::VoidPList(1,&HOST_memRadius));
		OCLDeviceGroup.CopyBuffer(OCL_Radius2  , OpenCL::HostToDevice, OpenCL::VoidPList(1,&HOST_memRadius2));
		OCLDeviceGroup.CopyBuffer(OCL_Particles, OpenCL::HostToDevice, OpenCL::VoidPList(1,&mParticles[processed]));

		OCLDeviceGroup.SetBufferAsKernelArg(OCLSearchNearest, 0,  OCL_IndexCellReferenceO);
		OCLDeviceGroup.SetBufferAsKernelArg(OCLSearchNearest, 1,  OCL_BinsObjectContainer);
		OCLDeviceGroup.SetBufferAsKernelArg(OCLSearchNearest, 2 , OCL_InvCellSize);
		OCLDeviceGroup.SetBufferAsKernelArg(OCLSearchNearest, 3 , OCL_N);
		OCLDeviceGroup.SetBufferAsKernelArg(OCLSearchNearest, 4,  OCL_Radius);
		OCLDeviceGroup.SetBufferAsKernelArg(OCLSearchNearest, 5,  OCL_Radius2);
		OCLDeviceGroup.SetBufferAsKernelArg(OCLSearchNearest, 6,  OCL_Particles);
		OCLDeviceGroup.SetBufferAsKernelArg(OCLSearchNearest, 7,  OCL_MinPoint);
		OCLDeviceGroup.SetBufferAsKernelArg(OCLSearchNearest, 8,  OCL_distance);
		OCLDeviceGroup.SetBufferAsKernelArg(OCLSearchNearest, 9,  OCL_results);

		OCLDeviceGroup.ExecuteKernel(OCLSearchNearest, amount);
	
		OCLDeviceGroup.CopyBuffer(OCL_results, OpenCL::DeviceToHost, OpenCL::VoidPList(1,resultPoints));

		result = resultPoints[amount-1];

		free(resultPoints);

		processed += amount;
	    }
	    
	    std::cout << "Nearest Point ID: " << result << std::endl;
	}
         
	 //************************************************************************
	 //************************************************************************

	 /// Turn back information as a string.
	 virtual std::string Info() const
	 {
	   return "BinsContainer";
	 }

	 /// Print information about this object.
	 virtual void PrintInfo(std::ostream& rOStream) const
	 {
	   rOStream << "BinsContainer";
	 }

	 /// Print object's data.
	 virtual void PrintData(std::ostream& rOStream, std::string const& Perfix = std::string()) const
	 {
	   rOStream << Perfix << "Bin[" << SearchUtils::PointerDistance(mPointBegin, mPointEnd) << "] : " << std::endl;
	   for(IteratorConstIterator i_cell = mIndexCell.begin() ; i_cell != mIndexCell.end()-1 ; i_cell++)
	   {
		 rOStream << Perfix << "[ " ;
		 for(IteratorType i_point = *i_cell ; i_point != *(i_cell+1) ; i_point++)
		   rOStream << **i_point << " ";
		 rOStream << "]" << std::endl;
	   }
	   rOStream << std::endl;
	 }

	 /// Print Size of Container
	 void PrintSize( std::ostream& rout ){
	   rout << " BinsSize: ";
	   for(SizeType i = 0 ; i < TDimension ; i++)
		 rout << "[" << mN[i] << "]";
	   rout << std::endl;
	 }

	 /// Print Limits Points of the Container
	 void PrintBox( std::ostream& rout ){
	   rout << " BinsBox: Min [";  mMinPoint.Print(rout);
	   rout <<       "];  Max [";  mMaxPoint.Print(rout);
	   rout <<       "];  Size ["; mCellSize.Print(rout);
	   rout << "]" << std::endl;
	 }

	 /// Assignment operator.
	 BinsObjectStaticOCL& operator=(BinsObjectStaticOCL const& rOther);

	 /// Copy constructor.
	 BinsObjectStaticOCL(BinsObjectStaticOCL const& rOther);

	 TPointType GetMinPoint()
	 {
	   TPointType point;
	   for(SizeType i = 0 ; i < TDimension ; i++)
		 point[i] = mMinPoint[i];
	   return point;
	 }

	 TPointType GetMaxPoint()
	 {
	   TPointType point;
	   for(SizeType i = 0 ; i < TDimension ; i++)
		 point[i] = mMaxPoint[i];
	   return point;
	 }

    private:
      
    // OCL mem and auxiliars
	OpenCL::DeviceGroup &OCLDeviceGroup;
	
	cl_int Err;
	
	cl_uint OCL_program;
	
	//Kernels
	cl_uint OCLGenerateBinsA;
	cl_uint OCLGenerateBinsB; 
	cl_uint OCLGenerateBinsB1; 
	cl_uint OCLGenerateBinsB2; 
	cl_uint OCLGenerateBinsB3;
	cl_uint OCLGenerateBinsC; 
	cl_uint OCLSearchInRadius;
	cl_uint OCLSearchTriangles;
	cl_uint OCLSearchTrianglesB;
	cl_uint OCLSearchNearest;
	cl_uint OCLUpdate;
	cl_uint OCLMove;
	cl_uint OCLScan_Local1;
	cl_uint OCLScan_Local2;
	cl_uint OCLuniformUpdate;
	
	cl_uint OCLTransferA;
	cl_uint OCLTransferB;
	cl_uint OCLTransferC;
	      
	//Kernel mem
	cl_uint OCL_PointsBins;
	cl_uint OCL_IndexCell;
	cl_uint OCL_IndexCellReference; 
	cl_uint OCL_IndexCellReferenceO;
	cl_uint OCL_InvCellSize;
	cl_uint OCL_N;
	cl_uint OCL_MinPoint; 
	cl_uint OCL_BinsObjectContainer; 
	cl_uint OCL_Radius; 
	cl_uint OCL_Radius2;
	cl_uint OCL_Active;
	cl_uint OCL_Scan_Buffer;
	
	cl_uint OCL_PointsTriangle;
	cl_uint OCL_TriangleList;
	cl_uint OCL_outData;
	cl_uint OCL_results; 
	cl_uint OCL_distance; 
	cl_uint OCL_NFunction;
	cl_uint OCL_Dens;
	
	cl_uint OCL_Particles;
	cl_uint OCL_ParticlesVelocity;
	cl_uint OCL_ParticlesVelocityOld;
	cl_uint OCL_ParticlesDisplace;
	cl_uint OCL_ParticlesAcceleration;
	cl_uint OCL_ParticlesForce;
	cl_uint OCL_ParticlesPress_Proj;
	cl_uint OCL_ParticlesLock;
        cl_uint OCL_NodesV;
	cl_uint OCL_FixedV;
	cl_uint OCL_NodesF;
	cl_uint OCL_NodesP;
	cl_uint OCL_NodesR;
	cl_uint OCL_NodesY;
	cl_uint OCL_NodesT;
	cl_uint OCL_Body_Force;
	cl_uint OCL_FLAG;
	
        cl_double4 * nodesV;
	int * FixedV;
	cl_double4 * nodesP;
	cl_double4 * nodesF;
	cl_double4 * nodesR;
	double     * nodesY;
	double     * nodesT;
	
	cl_double4 * PointsTriangles;
	
	int mCellSizeInt;
	int mBinsObjectSize;
	int mNodeSize;
	int mTriangleSize;
	int mProblemSize;
	
	//ModelPart
	ModelPart * mStaticMesh;
	ModelPart * mParticMesh;
    
        // Point and Barycenter Access Iterators (vector reordered!!)
	IteratorType mPointBegin;
	IteratorType mPointEnd;

        // Bin Parameters (Sizes,BoundingBox,...)
	PointType  mMinPoint;
	PointType  mMaxPoint;
	PointType  mCellSize;
	PointType  mInvCellSize;
	
	cl_double4 * mNodes;
	cl_double4 * mParticles;
	cl_double4 * mParticlesVelocity;
	cl_double4 * mParticlesVelocityOld;
	cl_double4 * mParticlesDisplace;
	cl_double4 * mParticlesForce;
	cl_double4 * mParticlesPressure;
	cl_double4 * mParticlesPress_Proj;
	cl_int4 * mTriangles;
	
	Tvector<SizeType,TDimension>  mN;

        // Bins Access Vector ( vector<Iterator> )
	IteratorVector           mIndexCell;
	IteratorVector           mIndexCellOCL;
	IteratorIterator         mIndexCellBegin;
	IteratorIterator         mIndexCellEnd;

    // Work Variables ( For non-copy of Search Variables )
    //SearchStructureType mBox;

   public:

	 static TreeNodeType* Construct(IteratorType PointsBegin, IteratorType PointsEnd, PointType MaxPoint, PointType MinPoint, SizeType BucketSize)
	 {
	   SizeType number_of_points = SearchUtils::PointerDistance(PointsBegin,PointsEnd);
	   if (number_of_points == 0)
		 return NULL;
	   else 
	   {
		 return new BinsObjectStaticOCL( PointsBegin, PointsEnd, MinPoint, MaxPoint, BucketSize );
	   }
	 }


};

template< std::size_t TDimension, class TPointType, class TContainerType, class TPointerType,
          class TIteratorType, class TDistanceIteratorType, class TDistanceFunction >
std::ostream & operator<<( std::ostream& rOStream, BinsObjectStaticOCL<TDimension,TPointType,TContainerType,TPointerType,TIteratorType,TDistanceIteratorType,TDistanceFunction>& rThis)
{
	rThis.PrintInfo(rOStream);
	rOStream << std::endl;
	rThis.PrintSize(rOStream);
	rThis.PrintData(rOStream);
	return rOStream;
}



}

#endif // KRATOS_BINS_CONTAINER_H_INCLUDE
