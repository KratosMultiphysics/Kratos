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

#define PARTICLE_BUFFER	76800

#include "../../../kratos/spatial_containers/tree.h"
#include "../custom_utilities/opencl_interface.h"
#include "processes/node_erase_process.h"
#include "includes/model_part.h"
#include "utilities/timer.h"
#include <malloc.h>
#include <stdio.h>
#include <string.h>
#include <queue>
#include "utilities/openmp_utils.h"


// Required for QueryPerformanceFrequency 
#if defined(_WIN64) || defined(_WIN32) || defined(WIN64) 
   #include <windows.h>
    #define CLOCK_REALTIME 1

	struct timespec
	{
		__int64 tv_sec,tv_nsec ;
	};
  
#endif

//#include <omp.h>
namespace Kratos {
#if defined(_WIN64) || defined(_WIN32) 
	void clock_gettime(int flag, timespec* t  )
	{
		LARGE_INTEGER frequency;
		LARGE_INTEGER start;
		LARGE_INTEGER end;

		//  Get the frequency
		QueryPerformanceFrequency(&frequency);

		//  Start timer
		QueryPerformanceCounter(&start);

		t->tv_sec = start.HighPart;
		t->tv_nsec = start.LowPart;
	}
#endif

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
	    mStaticMesh = StaticMesh;
	    mParticMesh = ParticMesh;
	    
	    cl_int4      * mTriangleItr;
	    IteratorType   mNodeItr;
	    cl_double4   * mParticlesItr;
	    cl_double4   * mParticlesVelocityItr;
	    cl_double4   * mParticlesDisplaceItr;
	    cl_double4   * mParticlesForceItr;
	    
	    std::cout << "Inicializando OCL" << std::endl;
	    InitOCL();
	 
	    mNodeSize	  = (StaticMesh->NodesEnd()-1)->Id();
	    mTriangleSize = (StaticMesh->ElementsEnd()-1)->Id();
	    mProblemSize  = ParticMesh->NumberOfNodes();
	    
	    mPointBegin = new PointType* [StaticMesh->NumberOfNodes()];
	    mPointEnd = mPointBegin + mNodeSize;
	    
	    mParticleBufferSize   = ceil(((float)mProblemSize / (float)PARTICLE_BUFFER)) * PARTICLE_BUFFER;
	    
	    recycleQueue = new std::queue<ModelPart::NodesContainerType::iterator>();
	    
	    mTriangleNodes 	  = (cl_int4 *   )malloc(sizeof(cl_int4)    * mTriangleSize	  );
	    mNodes		  = (cl_double4 *)malloc(sizeof(cl_double4) * mNodeSize           );
	    mParticles 		  = (cl_double4 *)malloc(sizeof(cl_double4) * mParticleBufferSize );
	    mParticlesVelocity 	  = (cl_double4 *)malloc(sizeof(cl_double4) * mParticleBufferSize );
	    mParticlesVelocityOld = (cl_double4 *)malloc(sizeof(cl_double4) * mParticleBufferSize );
	    mParticlesDisplace    = (cl_double4 *)malloc(sizeof(cl_double4) * mParticleBufferSize );
	    mParticlesForce 	  = (cl_double4 *)malloc(sizeof(cl_double4) * mParticleBufferSize );
	    mParticlesN 	  = (cl_double4 *)malloc(sizeof(cl_double4) * mParticleBufferSize );
	    mParticlesI		  = (int *	 )malloc(sizeof(int) 	    * mParticleBufferSize );
	    
	    mParticlesOnElement 	  = (int *	 )malloc(sizeof(int)	    * mStaticMesh->NumberOfElements()); 
	    
	    nodesV = (cl_double4 *)malloc(sizeof(cl_double4) * mNodeSize);
 	    nodesF = (cl_double4 *)malloc(sizeof(cl_double4) * mNodeSize);
	    nodesP = (cl_double4 *)malloc(sizeof(cl_double4) * mNodeSize);
 	    nodesR = (cl_double4 *)malloc(sizeof(cl_double4) * mNodeSize);
 	    nodesY = (    double *)malloc(sizeof(double)     * mNodeSize);
 	    nodesT = (    double *)malloc(sizeof(double)     * mNodeSize);
	    FixedV = (	     int *)malloc(sizeof(int) 	     * mNodeSize);
	    
	    std::cout << "\tNodes:     " << mNodeSize << std::endl;
	    std::cout << "\tElements:  " << mTriangleSize << std::endl;
	    std::cout << "\tParticles: " << mProblemSize << " + " <<  PARTICLE_BUFFER << std::endl;
	    
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
// 	      
	      int id = inode->Id()-1;
		
	      KRATOS_OCL_4_ARRAY_X(mNodes,id) = inode->X();
	      KRATOS_OCL_4_ARRAY_Y(mNodes,id) = inode->Y();
	      KRATOS_OCL_4_ARRAY_Z(mNodes,id) = inode->Z();
	      KRATOS_OCL_4_ARRAY_W(mNodes,id) = id+1;//(int)iel->GetGeometry()[3].Id();
	    }
	    
	    //StaticMesh Elements
	    std::cout << "Cargando Elementos" << std::endl;
	    
	    for( int i = 0; i < mTriangleSize; i++) {
		KRATOS_OCL_4_ARRAY_X(mTriangleNodes,i) = -1;
		KRATOS_OCL_4_ARRAY_Y(mTriangleNodes,i) = -1;
		KRATOS_OCL_4_ARRAY_Z(mTriangleNodes,i) = -1;
		KRATOS_OCL_4_ARRAY_W(mTriangleNodes,i) = -1;
	    }
	    
	    for( ModelPart::ElementsContainerType::iterator iel = mStaticMesh->ElementsBegin(); iel != mStaticMesh->ElementsEnd(); iel++, mTriangleItr++) 
	    {
		int id = iel->Id()-1;
		
		KRATOS_OCL_4_ARRAY_X(mTriangleNodes,id) = (int)iel->GetGeometry()[0].Id();
		KRATOS_OCL_4_ARRAY_Y(mTriangleNodes,id) = (int)iel->GetGeometry()[1].Id();
		KRATOS_OCL_4_ARRAY_Z(mTriangleNodes,id) = (int)iel->GetGeometry()[2].Id();
		KRATOS_OCL_4_ARRAY_W(mTriangleNodes,id) = 0;//(int)iel->GetGeometry()[3].Id();
	    }	
	    
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
	    Timer::SetOuputFile("cylinder_norberto_h.time");
	    Timer::PrintTimingInformation();
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
// 	  -cl-nv-maxrregcount=68
	   OCL_program = OCLDeviceGroup.BuildProgramFromFile("binshashobjects_tetrahedron.cl","-cl-nv-opt-level=3 -D WORKGROUP_SIZE=512");

	   OCLGenerateBinsA     = OCLDeviceGroup.RegisterKernel(OCL_program,"GenerateBinsObjectsA");
	   OCLScan_Local1	= OCLDeviceGroup.RegisterKernel(OCL_program,"scanExclusiveLocal1");
	   OCLScan_Local2	= OCLDeviceGroup.RegisterKernel(OCL_program,"scanExclusiveLocal2");
	   OCLGenerateBinsC     = OCLDeviceGroup.RegisterKernel(OCL_program,"GenerateBinsC");
	   OCLMove		= OCLDeviceGroup.RegisterKernel(OCL_program,"Move");
	   OCLMove2		= OCLDeviceGroup.RegisterKernel(OCL_program,"Move2");
	   OCLTransferB 	= OCLDeviceGroup.RegisterKernel(OCL_program,"calculateField");
	   OCLResetCounter	= OCLDeviceGroup.RegisterKernel(OCL_program,"resetCounter");
	   
	   std::cout << OCLDeviceGroup.WorkGroupSizes[OCLMove][0] << std::endl;
	   
	   std::cout << "OCL init ok" << std::endl;
	}
	 
	//************************************************************************
	
	void ComputeGaussPointPositions(Geometry< Node < 3 > >& geom, boost::numeric::ublas::bounded_matrix<double, 16, 3 > & pos, boost::numeric::ublas::bounded_matrix<double, 16, 3 > & N)
        {
            //lower diagonal terms
            double ypos = 1.0 / 12.0;
            int pos_counter = 0;
            for (unsigned int i = 0; i < 4; i++)
            {
                double xpos = 1.0 / 12.0;
                for (unsigned int j = 0; j < 4 - i; j++)
                {
                    double N1 = xpos;
                    double N2 = ypos;
                    double N3 = 1.0 - xpos - ypos;

                    pos(pos_counter, 0) = N1 * geom[0].X() + N2 * geom[1].X() + N3 * geom[2].X();
                    pos(pos_counter, 1) = N1 * geom[0].Y() + N2 * geom[1].Y() + N3 * geom[2].Y();
                    pos(pos_counter, 2) = N1 * geom[0].Z() + N2 * geom[1].Z() + N3 * geom[2].Z();

                    N(pos_counter, 0) = N1;
                    N(pos_counter, 1) = N2;
                    N(pos_counter, 2) = N3;

                    xpos += 1.0 / 4.0;
                    pos_counter += 1;

                }
                ypos += 1.0 / 4.0;
            }

            //lower diagonal terms
            ypos = 2.0 / 12.0;
            // pos_counter = 8;
            for (unsigned int i = 0; i < 3; i++)
            {
                double xpos = 2.0 / 12.0;
                for (unsigned int j = 0; j < 4 - i; j++)
                {
                    double N1 = xpos;
                    double N2 = ypos;
                    double N3 = 1.0 - xpos - ypos;

                    pos(pos_counter, 0) = N1 * geom[0].X() + N2 * geom[1].X() + N3 * geom[2].X();
                    pos(pos_counter, 1) = N1 * geom[0].Y() + N2 * geom[1].Y() + N3 * geom[2].Y();
                    pos(pos_counter, 2) = N1 * geom[0].Z() + N2 * geom[1].Z() + N3 * geom[2].Z();

                    N(pos_counter, 0) = N1;
                    N(pos_counter, 1) = N2;
                    N(pos_counter, 2) = N3;

                    xpos += 1.0 / 4.0;
                    pos_counter += 1;

                }
                ypos += 1.0 / 4.0;
            }
        }
        
        void ComputeGaussPointPositions(Geometry< Node < 3 > >& geom, boost::numeric::ublas::bounded_matrix<double, 4, 3 > & pos, boost::numeric::ublas::bounded_matrix<double, 4, 3 > & N)
        {
            double one_third = 1.0 / 3.0;
            double one_sixt = 1.0 / 6.0;
            double two_third = 2.0 * one_third;

            N(0, 0) = one_sixt;
            N(0, 1) = one_sixt;
            N(0, 2) = two_third;
            N(1, 0) = two_third;
            N(1, 1) = one_sixt;
            N(1, 2) = one_sixt;
            N(2, 0) = one_sixt;
            N(2, 1) = two_third;
            N(2, 2) = one_sixt;
            N(3, 0) = one_third;
            N(3, 1) = one_third;
            N(3, 2) = one_third;


            //first
            pos(0, 0) = one_sixt * geom[0].X() + one_sixt * geom[1].X() + two_third * geom[2].X();
            pos(0, 1) = one_sixt * geom[0].Y() + one_sixt * geom[1].Y() + two_third * geom[2].Y();
            pos(0, 2) = one_sixt * geom[0].Z() + one_sixt * geom[1].Z() + two_third * geom[2].Z();

            //second
            pos(1, 0) = two_third * geom[0].X() + one_sixt * geom[1].X() + one_sixt * geom[2].X();
            pos(1, 1) = two_third * geom[0].Y() + one_sixt * geom[1].Y() + one_sixt * geom[2].Y();
            pos(1, 2) = two_third * geom[0].Z() + one_sixt * geom[1].Z() + one_sixt * geom[2].Z();

            //third
            pos(2, 0) = one_sixt * geom[0].X() + two_third * geom[1].X() + one_sixt * geom[2].X();
            pos(2, 1) = one_sixt * geom[0].Y() + two_third * geom[1].Y() + one_sixt * geom[2].Y();
            pos(2, 2) = one_sixt * geom[0].Z() + two_third * geom[1].Z() + one_sixt * geom[2].Z();

            //fourth
            pos(3, 0) = one_third * geom[0].X() + one_third * geom[1].X() + one_third * geom[2].X();
            pos(3, 1) = one_third * geom[0].Y() + one_third * geom[1].Y() + one_third * geom[2].Y();
            pos(3, 2) = one_third * geom[0].Z() + one_third * geom[1].Z() + one_third * geom[2].Z();

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
	    mTriangleNodes = (cl_int4 *)malloc(sizeof(cl_int4) * mTriangleSize); 
	    
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
		mTriangleNodes[i] = indexArray[i];
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
	    
// 	    mN[0] *= 2;
// 	    mN[1] *= 2;
// 	    mN[2] *= 2;
	    
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
	    
	    OCL_Scan_Buffer = OCLDeviceGroup.CreateBuffer(sizeof(unsigned int) * ( 64 * ArrayLenght / (4 * OCLDeviceGroup.WorkGroupSizes[OCLScan_Local1][0])), CL_MEM_READ_ONLY );
	 
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
	    
	    PointsTriangles = (cl_double4 *)malloc(sizeof(cl_double4) * mNodeSize);
	    OCLDeviceGroup.CopyBuffer(OCL_PointsTriangle, OpenCL::HostToDevice, OpenCL::VoidPList(1,&mNodes[0]));
	    
	    OCLDeviceGroup.CopyBuffer(OCL_TriangleList      , OpenCL::HostToDevice, OpenCL::VoidPList(1,mTriangleNodes));
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
	
	void LoadSample()
	{
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
		for( ModelPart::NodesContainerType::iterator inode = mStaticMesh->NodesBegin(); inode != mStaticMesh->NodesEnd(); inode++) 
		{
		    int id = inode->Id()-1;
		  
		    array_1d<double, 3 > & velocity   = inode->FastGetSolutionStepValue(VELOCITY,1);
		    array_1d<double, 3 > & force      = inode->FastGetSolutionStepValue(FORCE);
		    array_1d<double, 3 > & press_proj = inode->FastGetSolutionStepValue(PRESS_PROJ);
		    
		    KRATOS_OCL_4_ARRAY_X(nodesV,id) = velocity[0];
		    KRATOS_OCL_4_ARRAY_Y(nodesV,id) = velocity[1];
		    KRATOS_OCL_4_ARRAY_Z(nodesV,id) = velocity[2];
		    KRATOS_OCL_4_ARRAY_W(nodesV,id) = id;
		    
		    KRATOS_OCL_4_ARRAY_X(nodesF,id) = force[0];
		    KRATOS_OCL_4_ARRAY_Y(nodesF,id) = force[1];
		    KRATOS_OCL_4_ARRAY_Z(nodesF,id) = force[2];
		    
		    KRATOS_OCL_4_ARRAY_X(nodesP,id) = press_proj[0];
		    KRATOS_OCL_4_ARRAY_Y(nodesP,id) = press_proj[1];
		    KRATOS_OCL_4_ARRAY_Z(nodesP,id) = press_proj[2];
		    
		    FixedV[id] = inode->IsFixed(VELOCITY_X);
		}
	    } else {
	      
		for( ModelPart::NodesContainerType::iterator inode = mStaticMesh->NodesBegin(); inode != mStaticMesh->NodesEnd(); inode++) 
		{
		    int id = inode->Id()-1;

		    array_1d<double, 3 > & velocity   = inode->FastGetSolutionStepValue(VELOCITY);
		    
		    velocity[0] = KRATOS_OCL_4_ARRAY_X(nodesV,id);
		    velocity[1] = KRATOS_OCL_4_ARRAY_Y(nodesV,id);
		    velocity[2] = KRATOS_OCL_4_ARRAY_Z(nodesV,id);
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
	}
	
	void TransferParticMeshToGPU() 
	{
	    if(mParticleBufferSize < static_cast<int>(mParticMesh->NumberOfNodes())) {
	      
		free(mParticles);
		free(mParticlesVelocity);
		free(mParticlesVelocityOld);
		free(mParticlesDisplace);
		free(mParticlesForce);
		free(mParticlesN);
		free(mParticlesI);
		
		while (mParticleBufferSize < static_cast<int>(mParticMesh->NumberOfNodes())) mParticleBufferSize += PARTICLE_BUFFER;
		
		mParticles 		  = (cl_double4 *)malloc(sizeof(cl_double4) * mParticleBufferSize );
		mParticlesVelocity 	  = (cl_double4 *)malloc(sizeof(cl_double4) * mParticleBufferSize );
		mParticlesVelocityOld 	  = (cl_double4 *)malloc(sizeof(cl_double4) * mParticleBufferSize );
		mParticlesDisplace    	  = (cl_double4 *)malloc(sizeof(cl_double4) * mParticleBufferSize );
		mParticlesForce 	  = (cl_double4 *)malloc(sizeof(cl_double4) * mParticleBufferSize );
		
		mParticlesN 	  	  = (cl_double4 *)malloc(sizeof(cl_double4) * mParticleBufferSize );
		mParticlesI		  = (int *	 )malloc(sizeof(int) 	    * mParticleBufferSize );
	    }
	  
	    cl_double4 * mParticlesItr 	       = mParticles;
	    cl_double4 * mParticlesVelocityItr = mParticlesVelocity;
	    cl_double4 * mParticlesDisplaceItr = mParticlesDisplace;
	  
	    for( ModelPart::NodesContainerType::iterator inode = mParticMesh->NodesBegin(); inode != mParticMesh->NodesEnd(); inode++) 
	    {
		if (!inode->GetValue(ERASE_FLAG)) {
		    Kratos::array_1d<double, 3> & velocity = inode->FastGetSolutionStepValue(VELOCITY,1);
		    Kratos::array_1d<double, 3> & displace = inode->FastGetSolutionStepValue(DISPLACEMENT,1);
			    
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
		    
		    mParticlesItr++, mParticlesVelocityItr++, mParticlesDisplaceItr++;
		}
		
// 		inode->GetValue(ERASE_FLAG) = false;
	    }
	}
	
	//************************************************************************
	  
	void AllocateOCLBuffers(int ConcurrentPoints/*, GLuint glBuffer*/)
	{
	    std::cout << "-Allocating Buffers-" << std::endl;
	  
	    int ConcurrentPointsReal = ConcurrentPoints;

	    std::cout << "\tAllocate Bufers: " << mProblemSize << std::endl;
	    std::cout << "\tAllocate Bufers Real: " << ConcurrentPointsReal << std::endl;
	    
	    ///////////////////////////////////////////////////////////////////////////
	    ConcurrentPointsReal = PARTICLE_BUFFER;
	    ///////////////////////////////////////////////////////////////////////////
	    
	    OCL_Radius        = OCLDeviceGroup.CreateBuffer(sizeof(double) * 2 , CL_MEM_READ_ONLY);
	    OCL_Body_Force    = OCLDeviceGroup.CreateBuffer(sizeof(cl_double4) , CL_MEM_READ_WRITE);

	    OCL_Particles             = OCLDeviceGroup.CreateBuffer(sizeof(cl_double4) 	* ConcurrentPointsReal, CL_MEM_READ_WRITE);
	    OCL_ParticlesVelocity     = OCLDeviceGroup.CreateBuffer(sizeof(cl_double4) 	* ConcurrentPointsReal, CL_MEM_READ_WRITE);
	    OCL_ParticlesVelocityOld  = OCLDeviceGroup.CreateBuffer(sizeof(cl_double4) 	* ConcurrentPointsReal, CL_MEM_READ_WRITE);
	    OCL_ParticlesDisplace     = OCLDeviceGroup.CreateBuffer(sizeof(cl_double4) 	* ConcurrentPointsReal, CL_MEM_READ_WRITE);
	    OCL_ParticlesDisplaceOld  = OCLDeviceGroup.CreateBuffer(sizeof(cl_double4) 	* ConcurrentPointsReal, CL_MEM_READ_WRITE);
	    OCL_ParticlesAcceleration = OCLDeviceGroup.CreateBuffer(sizeof(cl_double4) 	* ConcurrentPointsReal, CL_MEM_READ_WRITE);
	    OCL_ParticlesForce        = OCLDeviceGroup.CreateBuffer(sizeof(cl_double4) 	* ConcurrentPointsReal, CL_MEM_READ_WRITE);
	    OCL_ParticlesN            = OCLDeviceGroup.CreateBuffer(sizeof(cl_double4) 	* ConcurrentPointsReal, CL_MEM_READ_WRITE);
	    OCL_ParticlesLock	      = OCLDeviceGroup.CreateBuffer(sizeof(int)        	* ConcurrentPointsReal, CL_MEM_READ_WRITE);
	    OCL_ParticlesIndex        = OCLDeviceGroup.CreateBuffer(sizeof(int)        	* ConcurrentPointsReal, CL_MEM_READ_WRITE);
// 	    OCL_ElementContainer      = OCLDeviceGroup.CreateBuffer(sizeof(int)	       	* mTriangleSize, CL_MEM_READ_WRITE);

            OCL_NodesV        = OCLDeviceGroup.CreateBuffer(sizeof(cl_double4) * mNodeSize, CL_MEM_READ_WRITE);
	    OCL_NodesF        = OCLDeviceGroup.CreateBuffer(sizeof(cl_double4) * mNodeSize, CL_MEM_READ_WRITE);
	    OCL_NodesP        = OCLDeviceGroup.CreateBuffer(sizeof(cl_double4) * mNodeSize, CL_MEM_READ_WRITE);
	    OCL_NodesY        = OCLDeviceGroup.CreateBuffer(sizeof(double)     * mNodeSize, CL_MEM_READ_WRITE);
	    OCL_FixedV        = OCLDeviceGroup.CreateBuffer(sizeof(int)        * mNodeSize, CL_MEM_READ_WRITE);
	    OCL_results       = OCLDeviceGroup.CreateBuffer(sizeof(int)        * ConcurrentPointsReal, CL_MEM_WRITE_ONLY);
	    OCL_outData       = OCLDeviceGroup.CreateBuffer(sizeof(int)        * ConcurrentPointsReal, CL_MEM_WRITE_ONLY);
	    OCL_distance      = OCLDeviceGroup.CreateBuffer(sizeof(double)     * ConcurrentPointsReal, CL_MEM_READ_WRITE);
	}
	
	void InitializeBuffers(double Radius)
	{
	    std::cout << "-Initialize Buffers-" << std::endl;
	  
	    double HOST_memRadius[2];
	    
	    HOST_memRadius[0] = Radius;
	    HOST_memRadius[1] = Radius * Radius;
	    
	    OCLDeviceGroup.CopyBuffer(OCL_Radius        , OpenCL::HostToDevice, OpenCL::VoidPList(1,&HOST_memRadius));
	    OCLDeviceGroup.CopyBuffer(OCL_PointsTriangle, OpenCL::HostToDevice, OpenCL::VoidPList(1,&mNodes[0]));
	    
	    int * dens = (int *)malloc(sizeof(int) * mTriangleSize);
	    int * particleLock = (int *)malloc(sizeof(int) * (PARTICLE_BUFFER));
	    
	    for(int i = 0 ; i < mTriangleSize; i++) {
	      dens[i] = 0 ;
	    }
	    
	    for(int i = 0 ; i < PARTICLE_BUFFER; i++) {
	      particleLock[i] = 1;
	    }

	    OCLDeviceGroup.CopyBuffer(OCL_Particles, OpenCL::HostToDevice, OpenCL::VoidPList(1,&mParticles[0]));
	    OCLDeviceGroup.CopyBuffer(OCL_ParticlesDisplace, OpenCL::HostToDevice, OpenCL::VoidPList(1,&mParticlesDisplace[0]));
	    OCLDeviceGroup.CopyBuffer(OCL_ParticlesVelocity, OpenCL::HostToDevice, OpenCL::VoidPList(1,&mParticlesVelocity[0]));
	    OCLDeviceGroup.CopyBuffer(OCL_ParticlesVelocityOld, OpenCL::HostToDevice, OpenCL::VoidPList(1,&mParticlesVelocityOld[0]));
	    OCLDeviceGroup.CopyBuffer(OCL_ParticlesForce, OpenCL::HostToDevice, OpenCL::VoidPList(1,&mParticlesForce[0]));
	    OCLDeviceGroup.CopyBuffer(OCL_ParticlesLock, OpenCL::HostToDevice, OpenCL::VoidPList(1,&particleLock[0]));
	    
	}
         
	void SearchInRadiusOCL(double Radius, unsigned int ConcurrentPoints, unsigned int maxResults) 
	{
	    unsigned int processed = 0;
	    
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
	    }
	    
	    std::cout << "Total Results: " << result << " (" << maxResults * SearchUtils::PointerDistance(mPointBegin, mPointEnd) << " stored)"<< std::endl;
	}
	
	void SearchTriangles(array_1d<double, 3 > & body_force, const double density, const double dt, const double substeps, const int ConcurrentPoints, const int use_eulerian, const int copy_data, const int reseed) 
	{    
	    int amount = 0;
	    int processed = 0;
	    int oldParticleNum = 0;
	    
	    cl_double4 * mParticlesItr;
	    cl_double4 * mParticlesVelItr;
	    cl_double4 * mParticlesDisItr;
	    cl_double4 * mParticlesForItr;
	    
	    oldParticleNum = mParticMesh->NumberOfNodes();

            OCLDeviceGroup.SetBufferAsKernelArg(OCLMove, 0,  OCL_IndexCellReferenceO);
            OCLDeviceGroup.SetBufferAsKernelArg(OCLMove, 1,  OCL_BinsObjectContainer);
            OCLDeviceGroup.SetBufferAsKernelArg(OCLMove, 2,  OCL_PointsTriangle);
            OCLDeviceGroup.SetBufferAsKernelArg(OCLMove, 3,  OCL_TriangleList);
            OCLDeviceGroup.SetBufferAsKernelArg(OCLMove, 4,  OCL_InvCellSize);
            OCLDeviceGroup.SetBufferAsKernelArg(OCLMove, 5,  OCL_N);
            OCLDeviceGroup.SetBufferAsKernelArg(OCLMove, 6,  OCL_MinPoint);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLMove, 7,  OCL_Particles);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLMove, 8,  OCL_ParticlesVelocity);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLMove, 9, OCL_ParticlesDisplace);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLMove, 10, OCL_ParticlesForce);
 	    OCLDeviceGroup.SetBufferAsKernelArg(OCLMove, 11, OCL_NodesV);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLMove, 12, OCL_NodesF);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLMove, 13, OCL_NodesP);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLMove, 14, OCL_Body_Force);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLMove, 15, OCL_ParticlesVelocityOld);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLMove, 16, OCL_ParticlesDisplaceOld);
	    OCLDeviceGroup.SetKernelArg(	OCLMove, 17, 1/density);
	    OCLDeviceGroup.SetKernelArg(	OCLMove, 18, dt/substeps);
	    OCLDeviceGroup.SetKernelArg(	OCLMove, 19, substeps);
	    OCLDeviceGroup.SetKernelArg(	OCLMove, 20, use_eulerian);
	    
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLMove2, 0, OCL_Particles);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLMove2, 1, OCL_ParticlesVelocity);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLMove2, 2, OCL_ParticlesDisplace);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLMove2, 3, OCL_ParticlesVelocityOld);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLMove2, 4, OCL_ParticlesDisplaceOld);
	    OCLDeviceGroup.SetKernelArg(	OCLMove2, 5, dt); 
	    
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
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLTransferB, 11, OCL_ParticlesIndex);
	    OCLDeviceGroup.SetBufferAsKernelArg(OCLTransferB, 12, OCL_ParticlesN);
	    OCLDeviceGroup.SetKernelArg(	OCLTransferB, 13, PARTICLE_BUFFER);
	    OCLDeviceGroup.SetLocalMemAsKernelArg(OCLTransferB, 14, OCLDeviceGroup.WorkGroupSizes[OCLTransferB][0] * sizeof(int));
	    
	    while (processed < static_cast<int>(mParticMesh->NumberOfNodes()))
	    {	 
		amount = PARTICLE_BUFFER > mParticMesh->NumberOfNodes() ? mParticMesh->NumberOfNodes() : PARTICLE_BUFFER;
		
		Timer::Start("DeviceMemTransfer");
		OCLDeviceGroup.CopyBuffer(OCL_Particles, OpenCL::HostToDevice, OpenCL::VoidPList(1,&mParticles[processed]));
		OCLDeviceGroup.CopyBuffer(OCL_ParticlesDisplace, OpenCL::HostToDevice, OpenCL::VoidPList(1,&mParticlesDisplace[processed]));
		OCLDeviceGroup.CopyBuffer(OCL_ParticlesVelocity, OpenCL::HostToDevice, OpenCL::VoidPList(1,&mParticlesVelocity[processed]));
		OCLDeviceGroup.CopyBuffer(OCL_ParticlesVelocityOld, OpenCL::HostToDevice, OpenCL::VoidPList(1,&mParticlesVelocityOld[processed]));
		OCLDeviceGroup.CopyBuffer(OCL_ParticlesForce, OpenCL::HostToDevice, OpenCL::VoidPList(1,&mParticlesForce[processed]));
		OCLDeviceGroup.CopyBuffer(OCL_Body_Force, OpenCL::HostToDevice, OpenCL::VoidPList(1,&body_force));
		Timer::Stop("DeviceMemTransfer");
		
		Timer::Start("Gpu-Time");
		OCLDeviceGroup.ExecuteKernel(OCLMove, amount);
		OCLDeviceGroup.ExecuteKernel(OCLMove2, amount);
		Timer::Stop("Gpu-Time");
		
		Timer::Start("DeviceMemTransfer");
		OCLDeviceGroup.CopyBuffer(OCL_Particles           , OpenCL::DeviceToHost, OpenCL::VoidPList(1,&mParticles[processed]));
		OCLDeviceGroup.CopyBuffer(OCL_ParticlesVelocity   , OpenCL::DeviceToHost, OpenCL::VoidPList(1,&mParticlesVelocity[processed]));
		OCLDeviceGroup.CopyBuffer(OCL_ParticlesVelocityOld, OpenCL::DeviceToHost, OpenCL::VoidPList(1,&mParticlesVelocityOld[processed]));
		OCLDeviceGroup.CopyBuffer(OCL_ParticlesDisplace   , OpenCL::DeviceToHost, OpenCL::VoidPList(1,&mParticlesDisplace[processed]));
		OCLDeviceGroup.CopyBuffer(OCL_ParticlesForce      , OpenCL::DeviceToHost, OpenCL::VoidPList(1,&mParticlesForce[processed]));
		Timer::Stop("DeviceMemTransfer");
		
		Timer::Start("Gpu-Time");
		OCLDeviceGroup.ExecuteKernel(OCLTransferB, amount);
		Timer::Stop("Gpu-Time");
		
		Timer::Start("DeviceMemTransfer");
		OCLDeviceGroup.CopyBuffer(OCL_ParticlesN    , OpenCL::DeviceToHost, OpenCL::VoidPList(1,&mParticlesN[processed]));
		OCLDeviceGroup.CopyBuffer(OCL_ParticlesIndex, OpenCL::DeviceToHost, OpenCL::VoidPList(1,&mParticlesI[processed]));
		Timer::Stop("DeviceMemTransfer");

		processed += amount;
	    }
	    
	    Timer::Start("PythonToLibTransfer");
	    ModelPart::NodesContainerType::iterator inode = mParticMesh->NodesBegin();
	    for(int i = 0; inode != mParticMesh->NodesEnd(); inode++) 
	    {
		if (!inode->GetValue(ERASE_FLAG)) 
		{
		    if(use_eulerian) 
		    {
			Kratos::array_1d<double, 3> & velocity_old = inode->FastGetSolutionStepValue(VELOCITY,1);
			
			velocity_old[0] = KRATOS_OCL_4_ARRAY_X(mParticlesVelocityOld,i);
			velocity_old[1] = KRATOS_OCL_4_ARRAY_Y(mParticlesVelocityOld,i);
			velocity_old[2] = KRATOS_OCL_4_ARRAY_Z(mParticlesVelocityOld,i);
		    }
		    
		    i++;
		}
	    }
	    Timer::Stop("PythonToLibTransfer");

	    Timer::Start("PythonToLibTransfer");
	    if(copy_data) 
	    {
		mParticlesItr 	 = mParticles;
		mParticlesVelItr = mParticlesVelocity;
		mParticlesDisItr = mParticlesDisplace;
		mParticlesForItr = mParticlesForce;
		
		inode = mParticMesh->NodesBegin();

		for(int i = 0; inode != mParticMesh->NodesEnd(); inode++) 
		{
		    if (!inode->GetValue(ERASE_FLAG)) {
			Kratos::array_1d<double, 3> & velocity = inode->FastGetSolutionStepValue(VELOCITY);
			Kratos::array_1d<double, 3> & displace = inode->FastGetSolutionStepValue(DISPLACEMENT);
			Kratos::array_1d<double, 3> & force    = inode->FastGetSolutionStepValue(FORCE);
			
			velocity[0] = KRATOS_OCL_4_ARRAY_X(mParticlesVelocity,i);
			velocity[1] = KRATOS_OCL_4_ARRAY_Y(mParticlesVelocity,i);
			velocity[2] = KRATOS_OCL_4_ARRAY_Z(mParticlesVelocity,i);
			
			displace[0] = KRATOS_OCL_4_ARRAY_X(mParticlesDisplace,i);
			displace[1] = KRATOS_OCL_4_ARRAY_Y(mParticlesDisplace,i);
			displace[2] = KRATOS_OCL_4_ARRAY_Z(mParticlesDisplace,i);
			
			force[0] = KRATOS_OCL_4_ARRAY_X(mParticlesForce,i);
			force[1] = KRATOS_OCL_4_ARRAY_Y(mParticlesForce,i);
			force[2] = KRATOS_OCL_4_ARRAY_Z(mParticlesForce,i);
			
			inode->X() = KRATOS_OCL_4_ARRAY_X(mParticles,i);
			inode->Y() = KRATOS_OCL_4_ARRAY_Y(mParticles,i);
			inode->Z() = KRATOS_OCL_4_ARRAY_Z(mParticles,i);
			
			if(KRATOS_OCL_4_ARRAY_W(mParticles,i) < 0) 
			{
			    inode->SetValue(ERASE_FLAG,true);
			    inode->X() += 30;
			}
			
			i++;
		    }
		}
	    }
	    
	    Timer::Stop("PythonToLibTransfer");
	    
	    //Reseed
	    Timer::Start("Cpu-Time");
	    if(reseed) 
	    {
		mMinParticles = 4;
		
		int * mParticleIitr = mParticlesI;
		
		for (ModelPart::ElementsContainerType::iterator el_it = mStaticMesh->ElementsBegin(); el_it != mStaticMesh->ElementsEnd(); el_it++)
		{
		    el_it->SetValue(YOUNG_MODULUS,0.0);
		}
		
		for(ModelPart::NodesContainerType::iterator node = mParticMesh->NodesBegin(); node != mParticMesh->NodesEnd(); node++) 
		{
		    if(!node->GetValue(ERASE_FLAG)) 
		    {
			if((*mParticleIitr) > 0 ) {
			    ModelPart::ElementsContainerType::iterator el_it = mStaticMesh->ElementsBegin() + (*mParticleIitr);
			    double& counter = el_it->GetValue(YOUNG_MODULUS);
			    counter++;
			}
			mParticleIitr++;
		    }
		}

		mParticleIitr = mParticlesI;
		for(ModelPart::NodesContainerType::iterator node = mParticMesh->NodesBegin(); node != mParticMesh->NodesEnd(); node++) 
		{		  
		    if(!node->GetValue(ERASE_FLAG)) 
		    {
			
			if((*mParticleIitr) > 0) 
			{
			    ModelPart::ElementsContainerType::iterator el_it = mStaticMesh->ElementsBegin() + (*mParticleIitr);
			    if(!node->GetValue(ERASE_FLAG) && el_it->GetValue(YOUNG_MODULUS) < mMinParticles) {
				node->SetValue(ERASE_FLAG,true);
				node->X() += 30;
			    }
			}
			mParticleIitr++;
		    }
		}

		ModelPart::ElementsContainerType::iterator el_it = mStaticMesh->ElementsBegin();
		
// 		boost::numeric::ublas::bounded_matrix<double, 16, 3 > pos;
// 		boost::numeric::ublas::bounded_matrix<double, 16, 3 > Nnew;

		boost::numeric::ublas::bounded_matrix<double, 4, 3 > pos;
		boost::numeric::ublas::bounded_matrix<double, 4, 3 > Nnew;

		int id = (mParticMesh->NodesEnd() - 1)->Id();
		
		ModelPart::NodesContainerType::iterator no_it = mParticMesh->NodesBegin();
		
		for (ModelPart::ElementsContainerType::iterator el_it = mStaticMesh->ElementsBegin(); el_it != mStaticMesh->ElementsEnd(); el_it++)
		{
		    if (el_it->GetValue(YOUNG_MODULUS) < mMinParticles) 
		    {
			Geometry< Node<3> >& geom = el_it->GetGeometry();
			
			ComputeGaussPointPositions(geom, pos, Nnew);
			
// 			no_it = mParticMesh->NodesEnd();
			
			while(no_it != mParticMesh->NodesEnd() && !no_it->GetValue(ERASE_FLAG)) no_it++;
			
			for (unsigned int i = 0; i < pos.size1(); i++)
			{
			    if (no_it != mParticMesh->NodesEnd()) 
			    {
				no_it->SetValue(ERASE_FLAG,false);
				
				Kratos::array_1d<double, 3> newCoords;
				newCoords[0] = pos(i, 0);
				newCoords[1] = pos(i, 1);
				newCoords[2] = pos(i, 2);
				
				no_it->SetInitialPosition(newCoords);
				noalias(no_it->Coordinates()) = no_it->GetInitialPosition();
				
				array_1d<double, 3 > & dis = no_it->FastGetSolutionStepValue(DISPLACEMENT);
				noalias(dis) = ZeroVector(3);
				
				array_1d<double, 3 > & dis_old = no_it->FastGetSolutionStepValue(DISPLACEMENT,1);
				noalias(dis_old) = ZeroVector(3);
				
				array_1d<double, 3 > & force = no_it->FastGetSolutionStepValue(FORCE);
				noalias(force) = ZeroVector(3);

				array_1d<double, 3 > & vel = no_it->FastGetSolutionStepValue(VELOCITY);
				noalias(vel) = ZeroVector(3);
				for (unsigned int j = 0; j < 2 + 1; j++)
				  noalias(vel) += Nnew(i, j) * geom[j].FastGetSolutionStepValue(VELOCITY);
				  
				array_1d<double, 3 > & vel_old = no_it->FastGetSolutionStepValue(VELOCITY, 1);
				noalias(vel_old) = ZeroVector(3);
				for (unsigned int j = 0; j < 2 + 1; j++)
				  noalias(vel_old) += Nnew(i, j) * geom[j].FastGetSolutionStepValue(VELOCITY, 1);
			  
				while(no_it != mParticMesh->NodesEnd() && !no_it->GetValue(ERASE_FLAG)) no_it++;
			    } else {
				int node_id = id++;
				Node < 3 > ::Pointer pnode = mParticMesh->CreateNewNode(node_id, pos(i, 0), pos(i, 1), pos(i, 2));

				array_1d<double, 3 > & vel = pnode->FastGetSolutionStepValue(VELOCITY);
				noalias(vel) = ZeroVector(3);
				for (unsigned int j = 0; j < 2 + 1; j++)
				  noalias(vel) += Nnew(i, j) * geom[j].FastGetSolutionStepValue(VELOCITY);
				  
				array_1d<double, 3 > & vel_old = pnode->FastGetSolutionStepValue(VELOCITY, 1);
				noalias(vel_old) = ZeroVector(3);
				for (unsigned int j = 0; j < 2 + 1; j++)
				  noalias(vel_old) += Nnew(i, j) * geom[j].FastGetSolutionStepValue(VELOCITY, 1);
				
				no_it = mParticMesh->NodesEnd();
			    }
			}
		    }
		}
	    }
	    Timer::Stop("Cpu-Time");

	    //TransferToEulerianMeshShapeBased
	    Timer::Start("Cpu-Time");
	    for(int i = 0; i < mNodeSize; i++) 
	    {
		if (!FixedV[i]) 
		{
		    memset(&nodesV[i],0,sizeof(double)*3);
		    
		    nodesY[i] = 0;
		}
	    }	

	    for(int i = 0, j = 0; j < oldParticleNum; j++) 
	    {
		int tx, ty, tz;
		
		if ((mParticMesh->NodesBegin()+j)->GetValue(ERASE_FLAG) == false) 
		{
// 		    std::cout << "C-";
// 		    std::cout << mTriangleSize << " " << mParticlesI[i] << std::endl;
		  
		    int elementId = mParticlesI[i];
		      
		    if (elementId != -1) 
		    {
			tx = KRATOS_OCL_4_ARRAY_X(mTriangleNodes,(elementId))-1;
			ty = KRATOS_OCL_4_ARRAY_Y(mTriangleNodes,(elementId))-1;
			tz = KRATOS_OCL_4_ARRAY_Z(mTriangleNodes,(elementId))-1;
		      
			if(!FixedV[tx]) {
			    KRATOS_OCL_4_ARRAY_X(nodesV,tx) += KRATOS_OCL_4_ARRAY_X(mParticlesN,i) * KRATOS_OCL_4_ARRAY_X(mParticlesVelocity,i);
			    KRATOS_OCL_4_ARRAY_Y(nodesV,tx) += KRATOS_OCL_4_ARRAY_X(mParticlesN,i) * KRATOS_OCL_4_ARRAY_Y(mParticlesVelocity,i);
			    KRATOS_OCL_4_ARRAY_Z(nodesV,tx) += KRATOS_OCL_4_ARRAY_X(mParticlesN,i) * KRATOS_OCL_4_ARRAY_Z(mParticlesVelocity,i);
			    
			    nodesY[tx] += KRATOS_OCL_4_ARRAY_X(mParticlesN,i);
			}
			
			if(!FixedV[ty]) {
			    KRATOS_OCL_4_ARRAY_X(nodesV,ty) += KRATOS_OCL_4_ARRAY_Y(mParticlesN,i) * KRATOS_OCL_4_ARRAY_X(mParticlesVelocity,i);
			    KRATOS_OCL_4_ARRAY_Y(nodesV,ty) += KRATOS_OCL_4_ARRAY_Y(mParticlesN,i) * KRATOS_OCL_4_ARRAY_Y(mParticlesVelocity,i);
			    KRATOS_OCL_4_ARRAY_Z(nodesV,ty) += KRATOS_OCL_4_ARRAY_Y(mParticlesN,i) * KRATOS_OCL_4_ARRAY_Z(mParticlesVelocity,i);
			    
			    nodesY[ty] += KRATOS_OCL_4_ARRAY_Y(mParticlesN,i);
			}
			
			if(!FixedV[tz]) {
			    KRATOS_OCL_4_ARRAY_X(nodesV,tz) += KRATOS_OCL_4_ARRAY_Z(mParticlesN,i) * KRATOS_OCL_4_ARRAY_X(mParticlesVelocity,i);
			    KRATOS_OCL_4_ARRAY_Y(nodesV,tz) += KRATOS_OCL_4_ARRAY_Z(mParticlesN,i) * KRATOS_OCL_4_ARRAY_Y(mParticlesVelocity,i);
			    KRATOS_OCL_4_ARRAY_Z(nodesV,tz) += KRATOS_OCL_4_ARRAY_Z(mParticlesN,i) * KRATOS_OCL_4_ARRAY_Z(mParticlesVelocity,i);
			    
			    nodesY[tz] += KRATOS_OCL_4_ARRAY_Z(mParticlesN,i);
			} 
		     } 
		     
		     i++;
		}
	    }

	    for(int i = 0; i < mNodeSize; i++) 
	    { 
		if ( nodesY[i] != 0 && !FixedV[i]) 
		{
		    KRATOS_OCL_4_ARRAY_X(nodesV,i) /= nodesY[i];
		    KRATOS_OCL_4_ARRAY_Y(nodesV,i) /= nodesY[i];
		    KRATOS_OCL_4_ARRAY_Z(nodesV,i) /= nodesY[i];
		}
	    }
	    Timer::Stop("Cpu-Time");
	    
	    Timer::Start("Cpu-Time");
// 	    NodeEraseProcess(*mParticMesh).Execute();
	    Timer::Stop("Cpu-Time");
	    
// 	    ModelPart::NodesContainerType::iterator inode = mParticMesh->NodesBegin();
// 	    for(int i = 0; i < mParticMesh->NumberOfNodes(); i++, inode++) 
// 	    {
// 		Kratos::array_1d<double, 3> & velocity = inode->FastGetSolutionStepValue(VELOCITY,1);
// 		
// 		velocity[0] = KRATOS_OCL_4_ARRAY_X(mParticlesVelocityOld,i);
// 		velocity[1] = KRATOS_OCL_4_ARRAY_Y(mParticlesVelocityOld,i);
// 		velocity[2] = KRATOS_OCL_4_ARRAY_Z(mParticlesVelocityOld,i);
// 	    }
	    
	    std::cout << mParticMesh->NumberOfNodes() << std::endl;
	}
         
	void SearchNearestOCL(double Radius, unsigned int ConcurrentPoints) 
	{
	    unsigned int processed = 0;
	    
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
	cl_uint OCLMove;
	cl_uint OCLMove2;
	cl_uint OCLScan_Local1;
	cl_uint OCLScan_Local2;
	cl_uint OCLuniformUpdate;
	cl_uint OCLResetCounter;
	
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
	cl_uint OCL_Scan_Buffer;
	
	cl_uint OCL_PointsTriangle;
	cl_uint OCL_TriangleList;
// 	cl_uint OCL_ElementContainer;
	cl_uint OCL_outData;
	cl_uint OCL_results; 
	cl_uint OCL_distance; 
	
	cl_uint OCL_Particles;
	cl_uint OCL_ParticlesVelocity;
	cl_uint OCL_ParticlesVelocityOld;
	cl_uint OCL_ParticlesDisplace;
	cl_uint OCL_ParticlesDisplaceOld;
	cl_uint OCL_ParticlesAcceleration;
	cl_uint OCL_ParticlesIndex;
	cl_uint OCL_ParticlesN;
	cl_uint OCL_ParticlesForce;
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
	int mParticleBufferSize;
	int mMinParticles;
	
	int * BinsObjectContainer;
	int * mParticlesOnElement;
	
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
	cl_double4 * mParticlesN;
	
	int * mParticlesI;
	
	cl_int4 * mTriangleNodes;
	
	Tvector<SizeType,TDimension>  mN;

	IteratorVector    mIndexCell;
	IteratorVector    mIndexCellOCL;
	IteratorIterator  mIndexCellBegin;
	IteratorIterator  mIndexCellEnd;
	
	std::queue<ModelPart::NodesContainerType::iterator> * recycleQueue;

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
