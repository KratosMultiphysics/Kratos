//
//   Project Name:        Kratos
//   Last Modified by:    $Author: croig $
//   Date:                $Date: 2012-01-25 17:02:19 $
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

#define KRATOS_OCL_4_X(Arr)		(Arr.x)
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

#define PARTICLE_BUFFER	(116736)
#define MAX_BINS_ELEMENTS  64

#define KRATOS_OCL_MAP_PINNED_MEMORY_HTD(pointer_type,buffer,host_ptr,size)\
memcpy((void *)(host_ptr[0]),(const void *)(buffer),size*sizeof(pointer_type));

#define KRATOS_OCL_MAP_PINNED_MEMORY_DTH(pointer_type,buffer,host_ptr,size)\
memcpy((void *)(buffer),(const void *)(host_ptr[0]),size*sizeof(pointer_type));

#include "spatial_containers/tree.h"
#include "custom_utilities/opencl_interface.h"
#include "processes/node_erase_process.h"
#include "includes/model_part.h"
#include "includes/cfd_variables.h"
#include "utilities/timer.h"
#include "utilities/openmp_utils.h"

#include <malloc.h>
#include <stdio.h>
#include <string.h>
#include <queue>
#include <list>
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
namespace Kratos
{
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

//       typedef typename TreeNodeType::IteratorIteratorType IteratorIteratorType;
//       typedef typename TreeNodeType::SearchStructureType SearchStructureType;

//       typedef Kratos::SearchUtils::SearchNearestInRange<PointType,PointerType,IteratorType,DistanceFunction,CoordinateType> SearchNearestInRange;
//       typedef Kratos::SearchUtils::SearchRadiusInRange<PointType,IteratorType,DistanceIteratorType,DistanceFunction,SizeType,CoordinateType> SearchRadiusInRange;
//       typedef Kratos::SearchUtils::SearchBoxInRange<PointType,IteratorType,SizeType,TDimension> SearchBoxInRange;

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
        mParticlesDisplace    = (cl_double4 *)malloc(sizeof(cl_double4) * mParticleBufferSize );
        mParticlesIndex	  = (int *	 )malloc(sizeof(int) 	    * mParticleBufferSize );

        nodesV = (cl_double4 *)malloc(sizeof(cl_double4) * mNodeSize);
        nodesF = (cl_double4 *)malloc(sizeof(cl_double4) * mNodeSize);
        nodesP = (cl_double4 *)malloc(sizeof(cl_double4) * mNodeSize);
        nodesR = (cl_double4 *)malloc(sizeof(cl_double4) * mNodeSize);
        nodesY = (    double *)malloc(sizeof(double)     * mNodeSize);
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

        for(size_t i = 0; i < mTriangleSize; i++)
        {
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

        for( ModelPart::NodesContainerType::iterator inode = mParticMesh->NodesBegin(); inode != mParticMesh->NodesEnd(); inode++, mParticlesItr++, mParticlesVelocityItr++, mParticlesDisplaceItr++)
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
    virtual ~BinsObjectStaticOCL()
    {
        Timer::SetOuputFile("cylinder_norberto_h.time");
        Timer::PrintTimingInformation();
    }

    //************************************************************************

    IteratorType Begin()
    {
        return mPointBegin;
    }

    //************************************************************************

    IteratorType End()
    {
        return mPointBegin;
    }

    //************************************************************************

    CoordinateType CellSize( SizeType const& iDim )
    {
        return mCellSize[iDim];
    }

    //************************************************************************

    SizeType NumCell( SizeType const& iDim )
    {
        return mN[iDim];
    }

    //************************************************************************

private:

    //************************************************************************

    /// Initialize OpenCL interface
    void InitOCL()
    {
// 	  -cl-nv-maxrregcount=68
        OCL_program = OCLDeviceGroup.BuildProgramFromFile("binshashobjects_tetrahedron.cl","-D WORKGROUP_SIZE=512");

        OCLGenerateBinsA     	= OCLDeviceGroup.RegisterKernel(OCL_program,"GenerateBinsObjectsA");
        OCLGenerateBinsC     	= OCLDeviceGroup.RegisterKernel(OCL_program,"GenerateBinsObjectsC");
        OCLScan_Local1		= OCLDeviceGroup.RegisterKernel(OCL_program,"scanExclusiveLocal1");
        OCLScan_Local2		= OCLDeviceGroup.RegisterKernel(OCL_program,"scanExclusiveLocal2");
        OCLMove			= OCLDeviceGroup.RegisterKernel(OCL_program,"Move");
        OCLMove2			= OCLDeviceGroup.RegisterKernel(OCL_program,"Move2");
        OCLProject 			= OCLDeviceGroup.RegisterKernel(OCL_program,"calculateField");
        OCLResetCounter		= OCLDeviceGroup.RegisterKernel(OCL_program,"resetCounter");
        OCLCalculateParticleIndex 	= OCLDeviceGroup.RegisterKernel(OCL_program,"CalculateParticleIndex");

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

        for(SizeType i = 0 ; i < TDimension ; i++)
        {
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

// 	    for(int i = 1;;i*=2) if(mN[0] <= i) {mN[0] = i; break;}
// 	    for(int i = 1;;i*=2) if(mN[1] <= i) {mN[1] = i; break;}
// 	    for(int i = 1;;i*=2) if(mN[2] <= i) {mN[2] = i; break;}
        mN[0] *= 1;
        mN[1] *= 1;
        mN[2] = 1;

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

    /// Generates element bins from staticModelPart
    void GenerateBins()
    {
        Timer::Start("GenerateBins");

        unsigned int ArrayLenght = 1;

        while(ArrayLenght < mCellSizeInt+1)
            ArrayLenght <<= 1;

        cl_double4 MinPoint;

        int * IndexCellReference  = (int *)malloc(sizeof(int) * ArrayLenght);
        int * IndexCellReferenceO = (int *)malloc(sizeof(int) * ArrayLenght);

        double InvCellSize[TDimension];
        double N[TDimension];

        OCL_IndexCellReference = OCLDeviceGroup.CreateBuffer(sizeof(cl_uint4)   * ArrayLenght / 4   , CL_MEM_READ_WRITE);
        OCL_IndexCellReferenceO= OCLDeviceGroup.CreateBuffer(sizeof(cl_uint4)   * ArrayLenght / 4   , CL_MEM_READ_WRITE);
        OCL_InvCellSize        = OCLDeviceGroup.CreateBuffer(sizeof(double)     * TDimension        , CL_MEM_READ_ONLY );
        OCL_N                  = OCLDeviceGroup.CreateBuffer(sizeof(double)     * TDimension        , CL_MEM_READ_ONLY );
        OCL_MinPoint           = OCLDeviceGroup.CreateBuffer(sizeof(cl_double4)                     , CL_MEM_READ_ONLY );
        OCL_TriangleList   	   = OCLDeviceGroup.CreateBuffer(sizeof(cl_int4)    * mTriangleSize     , CL_MEM_READ_ONLY );
        OCL_PointsTriangle     = OCLDeviceGroup.CreateBuffer(sizeof(cl_double4) * mNodeSize         , CL_MEM_READ_ONLY );

        OCL_Scan_Buffer = OCLDeviceGroup.CreateBuffer(sizeof(unsigned int) * ( 64 * ArrayLenght / (4 * OCLDeviceGroup.WorkGroupSizes[OCLScan_Local1][0])), CL_MEM_READ_ONLY );

        KRATOS_OCL_4_X(MinPoint) = mMinPoint[0];
        KRATOS_OCL_4_Y(MinPoint) = mMinPoint[1];
        KRATOS_OCL_4_Z(MinPoint) = mMinPoint[2];

        for(size_t i = 0; i < ArrayLenght; i++)
        {
            IndexCellReference[i] = 0;
            IndexCellReferenceO[i] = 0;
        }

        for(size_t i = 0; i < TDimension; i++)
        {
            InvCellSize[i] = mInvCellSize[i];
            N[i] = mN[i];
        }

// 	    std::cout << "Var filling:\t\t" << ((float)(end.tv_sec - begin.tv_sec) + (float)(end.tv_nsec-begin.tv_nsec)/1000000000) << std::endl;

        PointsTriangles = (cl_double4 *)malloc(sizeof(cl_double4) * mNodeSize);
        OCLDeviceGroup.CopyBuffer(OCL_PointsTriangle, OpenCL::HostToDevice, OpenCL::VoidPList(1,&mNodes[0]));

        OCLDeviceGroup.CopyBuffer(OCL_TriangleList      , OpenCL::HostToDevice, OpenCL::VoidPList(1,mTriangleNodes));
        OCLDeviceGroup.CopyBuffer(OCL_IndexCellReference, OpenCL::HostToDevice, OpenCL::VoidPList(1,IndexCellReference));
        OCLDeviceGroup.CopyBuffer(OCL_InvCellSize	    , OpenCL::HostToDevice, OpenCL::VoidPList(1,InvCellSize));
        OCLDeviceGroup.CopyBuffer(OCL_N                 , OpenCL::HostToDevice, OpenCL::VoidPList(1,N));
        OCLDeviceGroup.CopyBuffer(OCL_MinPoint          , OpenCL::HostToDevice, OpenCL::VoidPList(1,&MinPoint));

        OCLDeviceGroup.SetBufferAsKernelArg(OCLGenerateBinsA, 0,  OCL_PointsTriangle);
        OCLDeviceGroup.SetBufferAsKernelArg(OCLGenerateBinsA, 1,  OCL_TriangleList);
        OCLDeviceGroup.SetBufferAsKernelArg(OCLGenerateBinsA, 2,  OCL_IndexCellReference);
        OCLDeviceGroup.SetBufferAsKernelArg(OCLGenerateBinsA, 3,  OCL_InvCellSize);
        OCLDeviceGroup.SetBufferAsKernelArg(OCLGenerateBinsA, 4,  OCL_N);
        OCLDeviceGroup.SetBufferAsKernelArg(OCLGenerateBinsA, 5,  OCL_MinPoint);
        OCLDeviceGroup.SetKernelArg(OCLGenerateBinsA, 6, mTriangleSize);

        OCLDeviceGroup.SetBufferAsKernelArg(OCLGenerateBinsC, 0,  OCL_PointsTriangle);
        OCLDeviceGroup.SetBufferAsKernelArg(OCLGenerateBinsC, 1,  OCL_TriangleList);
        OCLDeviceGroup.SetBufferAsKernelArg(OCLGenerateBinsC, 2,  OCL_IndexCellReferenceO);
        OCLDeviceGroup.SetBufferAsKernelArg(OCLGenerateBinsC, 3,  OCL_InvCellSize);
        OCLDeviceGroup.SetBufferAsKernelArg(OCLGenerateBinsC, 4,  OCL_N);
        OCLDeviceGroup.SetBufferAsKernelArg(OCLGenerateBinsC, 5,  OCL_MinPoint);
        OCLDeviceGroup.SetKernelArg(OCLGenerateBinsC, 7,  mTriangleSize);

        for(size_t i = 0; i < ArrayLenght; i++)
        {
            IndexCellReference[i] = 0;
            IndexCellReferenceO[i] = 0;
        }

        OCLDeviceGroup.CopyBuffer(OCL_IndexCellReference, OpenCL::HostToDevice, OpenCL::VoidPList(1,IndexCellReference));

        OCLDeviceGroup.ExecuteKernel(OCLGenerateBinsA, mTriangleSize);
        OCLDeviceGroup.Synchronize();
// 	    std::cout << "KERNEL A EXECUTED" << std::endl;

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

        int maxElementsInCell = IndexCellReference[0];
        int avgElementsInCell = IndexCellReference[0];

        IndexCellReferenceO[0] = 0;
        for(size_t i = 1; i < mCellSizeInt+1; i++)
        {
            avgElementsInCell += IndexCellReference[i];
            if (IndexCellReference[i] > maxElementsInCell) maxElementsInCell = IndexCellReference[i];
            IndexCellReferenceO[i] += IndexCellReferenceO[i-1] + IndexCellReference[i-1];
        }

        avgElementsInCell /= mCellSizeInt;

        OCLDeviceGroup.CopyBuffer(OCL_IndexCellReferenceO, OpenCL::HostToDevice, OpenCL::VoidPList(1,IndexCellReferenceO));

        mBinsObjectSize = IndexCellReferenceO[mCellSizeInt];

        // END OF SCAN PHASE

        OCL_BinsObjectContainer  = OCLDeviceGroup.CreateBuffer(sizeof(int) * mBinsObjectSize, CL_MEM_READ_WRITE);
        BinsObjectContainer  = (int *)malloc(sizeof(int) * mBinsObjectSize);
        BinsObjectContainerS = (int *)malloc(sizeof(int) * mBinsObjectSize);
        OCLDeviceGroup.SetBufferAsKernelArg(OCLGenerateBinsC, 6,  OCL_BinsObjectContainer);

        OCLDeviceGroup.ExecuteKernel(OCLGenerateBinsC, mTriangleSize);
        OCLDeviceGroup.Synchronize();

        OCLDeviceGroup.CopyBuffer(OCL_BinsObjectContainer, OpenCL::DeviceToHost, OpenCL::VoidPList(1,BinsObjectContainer));
        OCLDeviceGroup.CopyBuffer(OCL_IndexCellReferenceO, OpenCL::DeviceToHost, OpenCL::VoidPList(1,IndexCellReferenceO));

        int binsSize      = mN[0] * mN[1];
        int mBinsMaxLengh = mN[0] > mN[1] ? mN[0] : mN[1];
        int tmax = 0;

        int lastIndex = -1;

        level = 0;

        std::vector<std::list<cl_int4> > * mTriangleNodesSortedPointer = NULL;

        for(int i = 1; i < mBinsMaxLengh && tmax < MAX_BINS_ELEMENTS && i < 2; i*=2)
        {

            std::vector<std::list<cl_int4> > * tempTriangleNodesSorted = new std::vector<std::list<cl_int4> >(binsSize / i*i + 1);
            std::vector<uint> tempIndexCellReferenceElement = std::vector<uint>(binsSize / i*i + 1);

            for(size_t j = 0; j < mN[1] / i && tmax <= MAX_BINS_ELEMENTS; j++)
            {
                for(size_t k = 0; k < mN[0] / i && tmax <= MAX_BINS_ELEMENTS; k++)
                {
                    int tempMax = 0;
                    int binsIndex = j*mN[0]/i+k;

                    std::list<int> repeated;

                    for(int l = 0; l < i; l++)
                    {
                        int elementIndex = j*mN[0]*i+k*i+l*mN[0];

                        if((0 <= (elementIndex)) && ((elementIndex) < binsSize))
                        {
                            for(int m = 0; m < i; m++)
                            {
                                if((elementIndex+m)/mN[0] == elementIndex/mN[0])
                                {
// 					std::cout << "--" << elementIndex+m+1 << " " << elementIndex+m+2 << std::endl;
                                    for(int n = IndexCellReferenceO[elementIndex+m+1]; n < IndexCellReferenceO[elementIndex+m+2]; n++)
                                    {
                                        repeated.push_back(BinsObjectContainer[n]);
                                    }
                                }
                            }
                        }
                    }

                    repeated.sort();
                    repeated.unique();

                    (*tempTriangleNodesSorted)[binsIndex] = std::list<cl_int4>();

                    int storeIndex = lastIndex;
                    if (!(lastIndex != -1 && (*tempTriangleNodesSorted)[lastIndex].size()+repeated.size() < MAX_BINS_ELEMENTS ) )
                    {
                        storeIndex = binsIndex;
                    }
// 			storeIndex = binsIndex;s
// 				    std::cout << storeIndex << std::endl;

                    while(!repeated.empty())
                    {
                        (*tempTriangleNodesSorted)[storeIndex].push_back(mTriangleNodes[repeated.back()]);
                        repeated.pop_back();
                    }

                    tempMax = (*tempTriangleNodesSorted)[storeIndex].size();
                    tempIndexCellReferenceElement[storeIndex] = (*tempTriangleNodesSorted)[storeIndex].size();

                    tmax = tempMax > tmax ? tempMax : tmax;

                    lastIndex = storeIndex;
                }
            }

            std::cout << ((mN[1] / i)-1)*mN[0]/i+((mN[0] / i)-1) << "--" << tmax << std::endl;

            if(tmax <= MAX_BINS_ELEMENTS )
            {
                level++;

                mTriangleNodesSortedPointer = tempTriangleNodesSorted;
                IndexCellReferenceElement = tempIndexCellReferenceElement;
            }
// 		else
// 		{
            int size = (4 * binsSize) / (i*i);

            IndexCellReferenceElementO = std::vector<uint>(size + 1);
//
            for(int j = 0; j < size + 1; j++)
                IndexCellReferenceElementO[j] = 0;

            for(int j = 1; j < size + 1; j++)
            {
                IndexCellReferenceElementO[j] = IndexCellReferenceElementO[j-1] + IndexCellReferenceElement[j-1];
// 			std::cout << IndexCellReferenceElement[j-1] << " -- " << IndexCellReferenceElementO[j-1] << std::endl;
// 		    }
            }
        }

        macroCellSize = (int)(pow(2,level-1));
        mMacroBinsSize = binsSize / (macroCellSize * macroCellSize);
        mMacroBinsNumElements = IndexCellReferenceElementO[mMacroBinsSize];

        mMacroBins = std::vector<cl_int4>(mMacroBinsNumElements);

        for(int i = 0, j = 0; i < lastIndex+1; i++)
        {
            std::list<cl_int4> * elementList = &(*mTriangleNodesSortedPointer)[i];

            if (elementList != NULL && !elementList->empty())
            {
                if(elementList->size() != IndexCellReferenceElement[i])
                    std::cout << "Assertion fail: Size don't match" << elementList->size() << " " << IndexCellReferenceElement[i] << std::endl;

                while(!elementList->empty())
                {
                    mMacroBins[j++] = elementList->back();
                    elementList->pop_back();
                }
            }
        }

        OCL_IndexCellReferenceElementO = OCLDeviceGroup.CreateBuffer(sizeof(int) * mMacroBinsSize, CL_MEM_READ_WRITE);
        OCL_IndexCellReferenceSize     = OCLDeviceGroup.CreateBuffer(sizeof(int) * mMacroBinsSize, CL_MEM_READ_WRITE);

        OCL_MacroBins = OCLDeviceGroup.CreateBuffer(sizeof(cl_int4) * mMacroBinsNumElements , CL_MEM_READ_WRITE);

        OCLDeviceGroup.CopyBuffer(OCL_MacroBins                 , OpenCL::HostToDevice, OpenCL::VoidPList(1,&mMacroBins[0]));
        OCLDeviceGroup.CopyBuffer(OCL_IndexCellReferenceElementO, OpenCL::HostToDevice, OpenCL::VoidPList(1,&IndexCellReferenceElementO[0]));
        OCLDeviceGroup.CopyBuffer(OCL_IndexCellReferenceSize    , OpenCL::HostToDevice, OpenCL::VoidPList(1,&IndexCellReferenceElement[0]));

        mParticlesMatrix  = (cl_double4 **)malloc(sizeof(cl_double4 *) * (int)((mN[0] * mN[1]) / (macroCellSize * macroCellSize)));
        mParticlesMatrixRowCount = (uint *)malloc(sizeof(uint) * (mN[0] * mN[1]) / (macroCellSize * macroCellSize));

        for(int i = 0; i < lastIndex+1; i++)
            mParticlesMatrix[i] = (cl_double4 *)malloc(sizeof(cl_double4 ) * PARTICLE_BUFFER);

        OCLDeviceGroup.DeleteBuffer(OCL_IndexCellReference);

        Timer::Stop("GenerateBins");
    }

    //************************************************************************

    /// Load single particle into particle list
    void LoadSample()
    {
        cl_double4 * mParticlesItr;

        mParticlesItr = mParticles;
        CopyStaticmeshData(1);
        TransferStaticMeshToGPU();
    }

    /// Copy particles from staticModelPart to host memory
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
        }
        else
        {

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

    /// Transfer staticModelPart buffers back to GPU
    void TransferStaticMeshToGPU()
    {
        OCLDeviceGroup.CopyBuffer(OCL_NodesV, OpenCL::HostToDevice, OpenCL::VoidPList(1,&nodesV[0]));
        OCLDeviceGroup.CopyBuffer(OCL_NodesF, OpenCL::HostToDevice, OpenCL::VoidPList(1,&nodesF[0]));
        OCLDeviceGroup.CopyBuffer(OCL_NodesP, OpenCL::HostToDevice, OpenCL::VoidPList(1,&nodesP[0]));
        OCLDeviceGroup.CopyBuffer(OCL_FixedV, OpenCL::HostToDevice, OpenCL::VoidPList(1,&FixedV[0]));
    }

    /// Transfer staticModelPart buffers back to CPU
    void TransferStaticMeshToCPU()
    {
        OCLDeviceGroup.CopyBuffer(OCL_NodesV, OpenCL::DeviceToHost, OpenCL::VoidPList(1,&nodesV[0]));
    }

    /// Copy particles from particleModelPart to host memory
    void TransferParticMeshToGPU()
    {
        //Prevent malfunction while used with openMP
        int num_threads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::SetNumThreads(1);

        if(mParticleBufferSize < mParticMesh->NumberOfNodes())
        {
            free(mParticles);
            free(mParticlesVelocity);
            free(mParticlesDisplace);
            free(mParticlesIndex);

            while (mParticleBufferSize < mParticMesh->NumberOfNodes()) mParticleBufferSize += PARTICLE_BUFFER;

            mParticles 		  	= (cl_double4 *)malloc(sizeof(cl_double4) * mParticleBufferSize );
            mParticlesVelocity	= (cl_double4 *)malloc(sizeof(cl_double4) * mParticleBufferSize );
            mParticlesDisplace	= (cl_double4 *)malloc(sizeof(cl_double4) * mParticleBufferSize );
            mParticlesIndex		= (int *	   )malloc(sizeof(int) 	      * mParticleBufferSize );
        }

        cl_double4 * mParticlesItr 	       = mParticles;
        cl_double4 * mParticlesVelocityItr = mParticlesVelocity;
        cl_double4 * mParticlesDisplaceItr = mParticlesDisplace;

        int j = 0;
        for( ModelPart::NodesContainerType::iterator inode = mParticMesh->NodesBegin(); inode != mParticMesh->NodesEnd(); inode++)
        {
            if (!inode->Is(TO_ERASE))
            {
                Kratos::array_1d<double, 3> & velocity = inode->FastGetSolutionStepValue(VELOCITY,1);
                Kratos::array_1d<double, 3> & displace = inode->FastGetSolutionStepValue(DISPLACEMENT,1);

                noalias(inode->Coordinates()) = inode->GetInitialPosition();

                KRATOS_OCL_4_ITR_X(mParticlesItr) = inode->X();
                KRATOS_OCL_4_ITR_Y(mParticlesItr) = inode->Y();
                KRATOS_OCL_4_ITR_Z(mParticlesItr) = inode->Z();
                KRATOS_OCL_4_ITR_W(mParticlesItr) = (inode - mParticMesh->NodesBegin());

                KRATOS_OCL_4_ITR_X(mParticlesVelocityItr) = velocity[0];
                KRATOS_OCL_4_ITR_Y(mParticlesVelocityItr) = velocity[1];
                KRATOS_OCL_4_ITR_Z(mParticlesVelocityItr) = velocity[2];
                KRATOS_OCL_4_ITR_W(mParticlesVelocityItr) = inode->Id();

                KRATOS_OCL_4_ITR_X(mParticlesDisplaceItr) = displace[0];
                KRATOS_OCL_4_ITR_Y(mParticlesDisplaceItr) = displace[1];
                KRATOS_OCL_4_ITR_Z(mParticlesDisplaceItr) = displace[2];
                KRATOS_OCL_4_ITR_W(mParticlesDisplaceItr) = inode->Id();

                PointType auxPoint;

                auxPoint[0] = inode->X() + displace[0];
                auxPoint[1] = inode->Y() + displace[1];
                auxPoint[2] = inode->Z() + displace[2];

                mParticlesIndex[j] = CalculateIndex(auxPoint);

                mParticlesItr++;
                mParticlesVelocityItr++;
                mParticlesDisplaceItr++;
                j++;
            }
        }

        for(size_t i = 0; i < mMacroBinsSize; i++)
            mParticlesMatrixRowCount[i] = 0;

        for(size_t i = 0; i < mParticMesh->NumberOfNodes(); i++)
        {
            uint index = ((mParticlesIndex[i]/(mN[0]*macroCellSize)) * (mN[0]/macroCellSize)) + ((mParticlesIndex[i]%mN[0])/macroCellSize);

            while(IndexCellReferenceElement[index] == 0 && index > 0) index--;

            if (index > mMacroBinsSize)
            {
                std::cout << "Index out of bounds: " << index << std::endl;
                std::cout << "i: " << i << " " << mParticlesIndex[i] << std::endl;
                abort();
            }
            else
            {
                mParticlesMatrix[index][mParticlesMatrixRowCount[index]] = mParticles[i];
                mParticlesMatrixRowCount[index]++;
            }
        }

        mIndexSelectionSize = 0;

        for(size_t i = 0; i < mMacroBinsSize; i++)
            mIndexSelectionSize += (int)(mParticlesMatrixRowCount[i] / OCLDeviceGroup.WorkGroupSizes[OCLMove][0]) + 1;

        mIndexSelection = (int *)malloc(sizeof(int) * mIndexSelectionSize);

        OpenMPUtils::SetNumThreads(num_threads);
    }

    /// Allocate main gpu memory_objects
    void AllocateOCLBuffers(int ConcurrentPoints)
    {
        int ConcurrentPointsReal = ConcurrentPoints;

        ConcurrentPointsReal = PARTICLE_BUFFER;

        OCL_Radius        = OCLDeviceGroup.CreateBuffer(sizeof(double) * 2 , CL_MEM_READ_ONLY);

        OCL_Particles             		= OCLDeviceGroup.CreateBuffer(sizeof(cl_double4) * ConcurrentPointsReal, CL_MEM_READ_WRITE);
        OCL_ParticlesVelocity     		= OCLDeviceGroup.CreateBuffer(sizeof(cl_double4) * ConcurrentPointsReal, CL_MEM_READ_WRITE);
        OCL_ParticlesVelocityOld  		= OCLDeviceGroup.CreateBuffer(sizeof(cl_double4) * ConcurrentPointsReal, CL_MEM_READ_WRITE);
        OCL_ParticlesDisplace     		= OCLDeviceGroup.CreateBuffer(sizeof(cl_double4) * ConcurrentPointsReal, CL_MEM_READ_WRITE);
        OCL_ParticlesDisplaceOld  		= OCLDeviceGroup.CreateBuffer(sizeof(cl_double4) * ConcurrentPointsReal, CL_MEM_READ_WRITE);
        OCL_ParticlesAcceleration 		= OCLDeviceGroup.CreateBuffer(sizeof(cl_double4) * ConcurrentPointsReal, CL_MEM_READ_WRITE);
        OCL_ParticlesForce        		= OCLDeviceGroup.CreateBuffer(sizeof(cl_double4) * ConcurrentPointsReal, CL_MEM_READ_WRITE);
        OCL_ParticlesN            		= OCLDeviceGroup.CreateBuffer(sizeof(cl_double4) * ConcurrentPointsReal, CL_MEM_READ_WRITE);
        OCL_ParticlesIndex        		= OCLDeviceGroup.CreateBuffer(sizeof(int)        * ConcurrentPointsReal, CL_MEM_READ_WRITE);
        OCL_IndexTriangles			= OCLDeviceGroup.CreateBuffer(sizeof(cl_double4) * ConcurrentPointsReal, CL_MEM_READ_WRITE);

        OCL_Pinned_Particles      		= OCLDeviceGroup.CreateBuffer(sizeof(cl_double4) * ConcurrentPointsReal, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR);
        OCL_Pinned_ParticlesVelocity      	= OCLDeviceGroup.CreateBuffer(sizeof(cl_double4) * ConcurrentPointsReal, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR);
        OCL_Pinned_ParticlesVelocityOld     = OCLDeviceGroup.CreateBuffer(sizeof(cl_double4) * ConcurrentPointsReal, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR);
        OCL_Pinned_ParticlesDisplace     	= OCLDeviceGroup.CreateBuffer(sizeof(cl_double4) * ConcurrentPointsReal, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR);
        OCL_Pinned_ParticlesDisplaceOld     = OCLDeviceGroup.CreateBuffer(sizeof(cl_double4) * ConcurrentPointsReal, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR);

        OCL_Map_Particles 		 = OCLDeviceGroup.MapBuffer(OCL_Pinned_Particles,	    CL_MAP_READ | CL_MAP_WRITE);
        OCL_Map_ParticlesVelocity 	 = OCLDeviceGroup.MapBuffer(OCL_Pinned_ParticlesVelocity,   CL_MAP_READ | CL_MAP_WRITE);
        OCL_Map_ParticlesVelocityOld = OCLDeviceGroup.MapBuffer(OCL_Pinned_ParticlesVelocityOld,CL_MAP_READ | CL_MAP_WRITE);
        OCL_Map_ParticlesDisplace 	 = OCLDeviceGroup.MapBuffer(OCL_Pinned_ParticlesDisplace,   CL_MAP_READ | CL_MAP_WRITE);
        OCL_Map_ParticlesDisplaceOld = OCLDeviceGroup.MapBuffer(OCL_Pinned_ParticlesDisplaceOld,CL_MAP_READ | CL_MAP_WRITE);

        OCL_NodesV   = OCLDeviceGroup.CreateBuffer(sizeof(cl_double4) * mNodeSize, CL_MEM_READ_WRITE);
        OCL_NodesF   = OCLDeviceGroup.CreateBuffer(sizeof(cl_double4) * mNodeSize, CL_MEM_READ_WRITE);
        OCL_NodesP   = OCLDeviceGroup.CreateBuffer(sizeof(cl_double4) * mNodeSize, CL_MEM_READ_WRITE);
        OCL_NodesY   = OCLDeviceGroup.CreateBuffer(sizeof(double)     * mNodeSize, CL_MEM_READ_WRITE);
        OCL_FixedV   = OCLDeviceGroup.CreateBuffer(sizeof(int)        * mNodeSize, CL_MEM_READ_WRITE);
        OCL_results  = OCLDeviceGroup.CreateBuffer(sizeof(int)        * ConcurrentPointsReal, CL_MEM_WRITE_ONLY);
        OCL_outData  = OCLDeviceGroup.CreateBuffer(sizeof(int)        * ConcurrentPointsReal, CL_MEM_WRITE_ONLY);
        OCL_distance = OCLDeviceGroup.CreateBuffer(sizeof(double)     * ConcurrentPointsReal, CL_MEM_READ_WRITE);
    }

    /// Initialize buffers
    void InitializeBuffers(double Radius)
    {
        double HOST_memRadius[2];

        HOST_memRadius[0] = Radius;
        HOST_memRadius[1] = Radius * Radius;

        OCLDeviceGroup.CopyBuffer(OCL_Radius        , OpenCL::HostToDevice, OpenCL::VoidPList(1,&HOST_memRadius));
        OCLDeviceGroup.CopyBuffer(OCL_PointsTriangle, OpenCL::HostToDevice, OpenCL::VoidPList(1,&mNodes[0]));

        localBufferSize = PARTICLE_BUFFER * 2;

        mParticlesLocal 	        = (cl_double4 *)malloc(sizeof(cl_double4) * localBufferSize);
        mParticlesLocalVelocity     = (cl_double4 *)malloc(sizeof(cl_double4) * localBufferSize);
        mParticlesLocalVelocityOld  = (cl_double4 *)malloc(sizeof(cl_double4) * localBufferSize);
        mParticlesLocalDisplace     = (cl_double4 *)malloc(sizeof(cl_double4) * localBufferSize);
        mParticlesLocalForce 		= (cl_double4 *)malloc(sizeof(cl_double4) * localBufferSize);
        mParticlesLocalN            = (cl_double4 *)malloc(sizeof(cl_double4) * localBufferSize);
        mParticlesLocalI            = (int *       )malloc(sizeof(int)        * localBufferSize);
    }

    /// Search particles up to maxResults in Radius for each input
    void SearchInRadiusOCL(double Radius, unsigned int ConcurrentPoints, unsigned int maxResults)
    {
        unsigned int processed = 0;

        int result = 0;

        double HOST_memRadius  = Radius;
        double HOST_memRadius2 = Radius * Radius;

        unsigned int amount;

        while (processed < mProblemSize)
        {
            amount = (ConcurrentPoints > mProblemSize) ? mProblemSize : (processed + ConcurrentPoints) < mProblemSize ? ConcurrentPoints : mProblemSize - processed;

            int * results;
            int * resultPoints;

            results      = (int *)malloc(sizeof(int) * amount);
            resultPoints = (int *)malloc(sizeof(int) * amount * maxResults);

            OCLDeviceGroup.CopyBuffer(OCL_Radius   , OpenCL::HostToDevice, OpenCL::VoidPList(1,&HOST_memRadius));
            OCLDeviceGroup.CopyBuffer(OCL_Radius2  , OpenCL::HostToDevice, OpenCL::VoidPList(1,&HOST_memRadius2));
            OCLDeviceGroup.CopyBuffer(OCL_Particles, OpenCL::HostToDevice, OCL_Map_Particles);

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
            OCLDeviceGroup.Synchronize();

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

    /// Search particles
    void SearchParticles(array_1d<double, 3 > & body_force, const double density, const double dt, const double substeps, const int ConcurrentPoints, const int use_eulerian, const int copy_data, const int reseed)
    {
        //Prevent malfunction while used with openMP
        int num_threads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::SetNumThreads(1);

        unsigned int oldParticleNum = 0;

        KRATOS_OCL_4_X(mBody_Force) = body_force[0];
        KRATOS_OCL_4_Y(mBody_Force) = body_force[1];
        KRATOS_OCL_4_Z(mBody_Force) = body_force[2];
        KRATOS_OCL_4_X(mBody_Force) = 0;

        oldParticleNum = mParticMesh->NumberOfNodes();

        OCL_IndexSelection = OCLDeviceGroup.CreateBuffer(sizeof(int) * mIndexSelectionSize, CL_MEM_READ_WRITE);

        OCLDeviceGroup.SetBufferAsKernelArg(  OCLMove, 0,  OCL_PointsTriangle);
        OCLDeviceGroup.SetBufferAsKernelArg(  OCLMove, 1,  OCL_InvCellSize);
        OCLDeviceGroup.SetBufferAsKernelArg(  OCLMove, 2,  OCL_N);
        OCLDeviceGroup.SetBufferAsKernelArg(  OCLMove, 3,  OCL_MinPoint);
        OCLDeviceGroup.SetBufferAsKernelArg(  OCLMove, 4,  OCL_Particles);
        OCLDeviceGroup.SetBufferAsKernelArg(  OCLMove, 5,  OCL_ParticlesVelocity);
        OCLDeviceGroup.SetBufferAsKernelArg(  OCLMove, 6,  OCL_ParticlesDisplace);
        OCLDeviceGroup.SetBufferAsKernelArg(  OCLMove, 7,  OCL_ParticlesForce);
        OCLDeviceGroup.SetBufferAsKernelArg(  OCLMove, 8,  OCL_NodesV);
        OCLDeviceGroup.SetBufferAsKernelArg(  OCLMove, 9,  OCL_NodesF);
        OCLDeviceGroup.SetBufferAsKernelArg(  OCLMove, 10, OCL_NodesP);
        OCLDeviceGroup.SetKernelArg(          OCLMove, 11, mBody_Force);
        OCLDeviceGroup.SetBufferAsKernelArg(  OCLMove, 12, OCL_ParticlesVelocityOld);
        OCLDeviceGroup.SetBufferAsKernelArg(  OCLMove, 13, OCL_ParticlesDisplaceOld);
        OCLDeviceGroup.SetKernelArg(	  OCLMove, 14, 1/density);
        OCLDeviceGroup.SetKernelArg(	  OCLMove, 15, dt/substeps);
        OCLDeviceGroup.SetKernelArg(	  OCLMove, 16, substeps);
        OCLDeviceGroup.SetKernelArg(	  OCLMove, 17, use_eulerian);
        OCLDeviceGroup.SetLocalMemAsKernelArg(OCLMove, 18, sizeof(cl_double4) * MAX_BINS_ELEMENTS*4);
        OCLDeviceGroup.SetBufferAsKernelArg(  OCLMove, 19, OCL_IndexCellReferenceElementO);
        OCLDeviceGroup.SetBufferAsKernelArg(  OCLMove, 20, OCL_IndexCellReferenceSize);
        OCLDeviceGroup.SetBufferAsKernelArg(  OCLMove, 21, OCL_MacroBins);
        OCLDeviceGroup.SetBufferAsKernelArg(  OCLMove, 22, OCL_IndexSelection );
        OCLDeviceGroup.SetBufferAsKernelArg(  OCLMove, 23, OCL_IndexTriangles );

        OCLDeviceGroup.SetBufferAsKernelArg(OCLMove2, 0, OCL_Particles);
        OCLDeviceGroup.SetBufferAsKernelArg(OCLMove2, 1, OCL_ParticlesVelocity);
        OCLDeviceGroup.SetBufferAsKernelArg(OCLMove2, 2, OCL_ParticlesDisplace);
        OCLDeviceGroup.SetBufferAsKernelArg(OCLMove2, 3, OCL_ParticlesVelocityOld);
        OCLDeviceGroup.SetBufferAsKernelArg(OCLMove2, 4, OCL_ParticlesDisplaceOld);
        OCLDeviceGroup.SetKernelArg(	OCLMove2, 5, dt);

        OCLDeviceGroup.SetBufferAsKernelArg(OCLProject, 0,  OCL_IndexCellReferenceO);
        OCLDeviceGroup.SetBufferAsKernelArg(OCLProject, 1,  OCL_BinsObjectContainer);
        OCLDeviceGroup.SetBufferAsKernelArg(OCLProject, 2,  OCL_PointsTriangle);
        OCLDeviceGroup.SetBufferAsKernelArg(OCLProject, 3,  OCL_TriangleList);
        OCLDeviceGroup.SetBufferAsKernelArg(OCLProject, 4 , OCL_InvCellSize);
        OCLDeviceGroup.SetBufferAsKernelArg(OCLProject, 5 , OCL_N);
        OCLDeviceGroup.SetBufferAsKernelArg(OCLProject, 6,  OCL_Radius);
        OCLDeviceGroup.SetBufferAsKernelArg(OCLProject, 7,  OCL_MinPoint);
        OCLDeviceGroup.SetBufferAsKernelArg(OCLProject, 8,  OCL_Particles);
        OCLDeviceGroup.SetBufferAsKernelArg(OCLProject, 9,  OCL_ParticlesVelocity);
        OCLDeviceGroup.SetBufferAsKernelArg(OCLProject, 10, OCL_NodesV);
        OCLDeviceGroup.SetBufferAsKernelArg(OCLProject, 11, OCL_ParticlesIndex);
        OCLDeviceGroup.SetBufferAsKernelArg(OCLProject, 12, OCL_ParticlesN);
        OCLDeviceGroup.SetLocalMemAsKernelArg(OCLProject, 14, OCLDeviceGroup.WorkGroupSizes[OCLProject][0] * sizeof(int));

        unsigned int MoveKernelWorkGroupSzie = OCLDeviceGroup.WorkGroupSizes[OCLMove][0];
        unsigned int limitParticles 	= PARTICLE_BUFFER;
        unsigned int processed 		= 0;
        unsigned int macroCellItr    	= 0;

        unsigned int newlocalBufferSize = ceil(((float)mIndexSelectionSize * (float)MoveKernelWorkGroupSzie)/ (float)PARTICLE_BUFFER) * PARTICLE_BUFFER;

        if (newlocalBufferSize > localBufferSize)
        {
            localBufferSize = newlocalBufferSize;

            free(mParticlesLocal);
            free(mParticlesLocalVelocity);
            free(mParticlesLocalVelocityOld);
            free(mParticlesLocalDisplace);
            free(mParticlesLocalForce);
            free(mParticlesLocalN);
            free(mParticlesLocalI);

            mParticlesLocal 	    	= (cl_double4 *)malloc(sizeof(cl_double4) * localBufferSize);
            mParticlesLocalVelocity     = (cl_double4 *)malloc(sizeof(cl_double4) * localBufferSize);
            mParticlesLocalVelocityOld  = (cl_double4 *)malloc(sizeof(cl_double4) * localBufferSize);
            mParticlesLocalDisplace     = (cl_double4 *)malloc(sizeof(cl_double4) * localBufferSize);
            mParticlesLocalForce 	    = (cl_double4 *)malloc(sizeof(cl_double4) * localBufferSize);
            mParticlesLocalN	    	= (cl_double4 *)malloc(sizeof(cl_double4) * localBufferSize);
            mParticlesLocalI 	    	= (int *       )malloc(sizeof(int)        * localBufferSize);
        }

        int mIndexSelectionItr   = 0;
        int mIndexSelectionIndex = 0;

        unsigned int localIndex 	     = 0;

        while(processed < mParticMesh->NumberOfNodes())
        {
            unsigned int localParticles  = 0;

            while (macroCellItr < mMacroBinsSize && localParticles + mParticlesMatrixRowCount[macroCellItr] < limitParticles && macroCellItr < mMacroBinsSize)
            {
                unsigned int particleIndex = 0;

                for(particleIndex = 0; macroCellItr < mMacroBinsSize && particleIndex < mParticlesMatrixRowCount[macroCellItr] && localParticles < limitParticles; particleIndex++, localParticles++)
                {
                    mParticlesLocal[localIndex + localParticles] = mParticlesMatrix[macroCellItr][particleIndex];

                    ModelPart::NodesContainerType::iterator inode = mParticMesh->NodesBegin() + KRATOS_OCL_4_ARRAY_W(mParticlesLocal,localIndex + localParticles);

                    Kratos::array_1d<double, 3> & velocity = inode->FastGetSolutionStepValue(VELOCITY    ,1);
                    Kratos::array_1d<double, 3> & displace = inode->FastGetSolutionStepValue(DISPLACEMENT,1);

                    KRATOS_OCL_4_ARRAY_X(mParticlesLocalVelocity,localIndex + localParticles) = velocity[0];
                    KRATOS_OCL_4_ARRAY_Y(mParticlesLocalVelocity,localIndex + localParticles) = velocity[1];
                    KRATOS_OCL_4_ARRAY_Z(mParticlesLocalVelocity,localIndex + localParticles) = velocity[2];
                    KRATOS_OCL_4_ARRAY_W(mParticlesLocalVelocity,localIndex + localParticles) = inode->Id();

                    KRATOS_OCL_4_ARRAY_X(mParticlesLocalDisplace,localIndex + localParticles) = displace[0];
                    KRATOS_OCL_4_ARRAY_Y(mParticlesLocalDisplace,localIndex + localParticles) = displace[1];
                    KRATOS_OCL_4_ARRAY_Z(mParticlesLocalDisplace,localIndex + localParticles) = displace[2];
                    KRATOS_OCL_4_ARRAY_W(mParticlesLocalDisplace,localIndex + localParticles) = inode->Id();

                    if((localIndex + localParticles)%MoveKernelWorkGroupSzie == 0)
                    {
                        mIndexSelection[mIndexSelectionItr++] = macroCellItr;
                    }
                }

                /// We must fill the remaining buffer
                for(; (localIndex + localParticles)%MoveKernelWorkGroupSzie != 0 && localParticles < limitParticles; localParticles++)
                {
                    KRATOS_OCL_4_ARRAY_W(mParticlesLocal,localIndex + localParticles) = -2;
                }

                processed += mParticlesMatrixRowCount[macroCellItr];
                macroCellItr++;
            }

            if(localIndex + PARTICLE_BUFFER > localBufferSize)
            {
                std::cout << "Possible heap corruption detected, application will halt now." << localIndex + PARTICLE_BUFFER << " " << localBufferSize << std::endl;
                abort();
            }

            Timer::Start("ReadMem");
            KRATOS_OCL_MAP_PINNED_MEMORY_HTD(cl_double4,(&mParticlesLocal[localIndex])        ,OCL_Map_Particles        ,localParticles);
            KRATOS_OCL_MAP_PINNED_MEMORY_HTD(cl_double4,(&mParticlesLocalVelocity[localIndex]),OCL_Map_ParticlesVelocity,localParticles);
            KRATOS_OCL_MAP_PINNED_MEMORY_HTD(cl_double4,(&mParticlesLocalDisplace[localIndex]),OCL_Map_ParticlesDisplace,localParticles);
            Timer::Stop("ReadMem");

            Timer::Start("DeviceMemTransfer");
            OCLDeviceGroup.CopyBuffer(OCL_IndexSelection,           OpenCL::HostToDevice, OpenCL::VoidPList(1,&mIndexSelection[mIndexSelectionIndex]));
            OCLDeviceGroup.CopyBuffer(OCL_Particles,	 	OpenCL::HostToDevice, OCL_Map_Particles);
            OCLDeviceGroup.CopyBuffer(OCL_ParticlesDisplace, 	OpenCL::HostToDevice, OCL_Map_ParticlesDisplace);
            OCLDeviceGroup.CopyBuffer(OCL_ParticlesVelocity, 	OpenCL::HostToDevice, OCL_Map_ParticlesVelocity);
            Timer::Stop("DeviceMemTransfer");

            Timer::Start("Gpu-Time");
            OCLDeviceGroup.ExecuteKernel(OCLMove,  localParticles);
            OCLDeviceGroup.Synchronize();
            OCLDeviceGroup.ExecuteKernel(OCLMove2, localParticles);
            OCLDeviceGroup.Synchronize();
            Timer::Stop("Gpu-Time");

            Timer::Start("DeviceMemTransfer");
            OCLDeviceGroup.CopyBuffer(OCL_Particles           , OpenCL::DeviceToHost, OCL_Map_Particles);
            OCLDeviceGroup.CopyBuffer(OCL_ParticlesVelocity   , OpenCL::DeviceToHost, OCL_Map_ParticlesVelocity);
            OCLDeviceGroup.CopyBuffer(OCL_ParticlesVelocityOld, OpenCL::DeviceToHost, OCL_Map_ParticlesVelocityOld);
            OCLDeviceGroup.CopyBuffer(OCL_ParticlesDisplace   , OpenCL::DeviceToHost, OCL_Map_ParticlesDisplace);
            OCLDeviceGroup.CopyBuffer(OCL_ParticlesForce      , OpenCL::DeviceToHost, OpenCL::VoidPList(1,&mParticlesLocalForce[localIndex]));
            Timer::Stop("DeviceMemTransfer");

            Timer::Start("ReadMem");
            KRATOS_OCL_MAP_PINNED_MEMORY_DTH(cl_double4,(&mParticlesLocal[localIndex]),OCL_Map_Particles,localParticles);
            KRATOS_OCL_MAP_PINNED_MEMORY_DTH(cl_double4,(&mParticlesLocalVelocity[localIndex]),OCL_Map_ParticlesVelocity,localParticles);
            KRATOS_OCL_MAP_PINNED_MEMORY_DTH(cl_double4,(&mParticlesLocalVelocityOld[localIndex]),OCL_Map_ParticlesVelocityOld,localParticles);
            KRATOS_OCL_MAP_PINNED_MEMORY_DTH(cl_double4,(&mParticlesLocalDisplace[localIndex]),OCL_Map_ParticlesDisplace,localParticles);
            Timer::Stop("ReadMem");

            OCLDeviceGroup.SetKernelArg(OCLProject, 13, localParticles);

            Timer::Start("Gpu-Time");
            OCLDeviceGroup.ExecuteKernel(OCLProject, localParticles);
            OCLDeviceGroup.Synchronize();
            Timer::Stop("Gpu-Time");

            Timer::Start("DeviceMemTransfer");
            OCLDeviceGroup.CopyBuffer(OCL_ParticlesN    , OpenCL::DeviceToHost, OpenCL::VoidPList(1,&mParticlesLocalN[localIndex]));
            OCLDeviceGroup.CopyBuffer(OCL_ParticlesIndex, OpenCL::DeviceToHost, OpenCL::VoidPList(1,&mParticlesLocalI[localIndex]));
            Timer::Stop("DeviceMemTransfer");

            localIndex += localParticles;

            mIndexSelectionIndex = mIndexSelectionItr;
        }

        Timer::Start("PythonToLibTransfer");

        if(copy_data)
        {
            for(size_t i = 0;  i < localIndex; i++)
            {
                if (KRATOS_OCL_4_ARRAY_W(mParticlesLocal,i) != -2)
                {
                    ModelPart::NodesContainerType::iterator inode = mParticMesh->NodesBegin() + abs(KRATOS_OCL_4_ARRAY_W(mParticlesLocal,i));

                    Kratos::array_1d<double, 3> & oldvelocity = inode->FastGetSolutionStepValue(VELOCITY,1);
                    Kratos::array_1d<double, 3> & velocity = inode->FastGetSolutionStepValue(VELOCITY);
                    Kratos::array_1d<double, 3> & displace = inode->FastGetSolutionStepValue(DISPLACEMENT);
                    Kratos::array_1d<double, 3> & force    = inode->FastGetSolutionStepValue(FORCE);

                    oldvelocity[0] = KRATOS_OCL_4_ARRAY_X(mParticlesLocalVelocityOld,i);
                    oldvelocity[1] = KRATOS_OCL_4_ARRAY_Y(mParticlesLocalVelocityOld,i);
                    oldvelocity[2] = KRATOS_OCL_4_ARRAY_Z(mParticlesLocalVelocityOld,i);

                    velocity[0] = KRATOS_OCL_4_ARRAY_X(mParticlesLocalVelocity,i);
                    velocity[1] = KRATOS_OCL_4_ARRAY_Y(mParticlesLocalVelocity,i);
                    velocity[2] = KRATOS_OCL_4_ARRAY_Z(mParticlesLocalVelocity,i);

                    displace[0] = KRATOS_OCL_4_ARRAY_X(mParticlesLocalDisplace,i);
                    displace[1] = KRATOS_OCL_4_ARRAY_Y(mParticlesLocalDisplace,i);
                    displace[2] = KRATOS_OCL_4_ARRAY_Z(mParticlesLocalDisplace,i);

                    force[0] = KRATOS_OCL_4_ARRAY_X(mParticlesLocalForce,i);
                    force[1] = KRATOS_OCL_4_ARRAY_Y(mParticlesLocalForce,i);
                    force[2] = KRATOS_OCL_4_ARRAY_Z(mParticlesLocalForce,i);

                    inode->X() = KRATOS_OCL_4_ARRAY_X(mParticlesLocal,i);
                    inode->Y() = KRATOS_OCL_4_ARRAY_Y(mParticlesLocal,i);
                    inode->Z() = KRATOS_OCL_4_ARRAY_Z(mParticlesLocal,i);

                    if(KRATOS_OCL_4_ARRAY_W(mParticlesLocal,i) < 0)
                    {
                        inode->Set(TO_ERASE,true);
                        inode->X() += 30;
                    }
                }
            }
        }
        else
        {
            for(size_t i = 0;  i < localIndex; i++)
            {
                if (KRATOS_OCL_4_ARRAY_W(mParticlesLocal,i) != -2)
                {
                    ModelPart::NodesContainerType::iterator inode = mParticMesh->NodesBegin() + abs(KRATOS_OCL_4_ARRAY_W(mParticlesLocal,i));

                    Kratos::array_1d<double, 3> & oldvelocity = inode->FastGetSolutionStepValue(VELOCITY,1);

                    oldvelocity[0] = KRATOS_OCL_4_ARRAY_X(mParticlesLocalVelocityOld,i);
                    oldvelocity[1] = KRATOS_OCL_4_ARRAY_Y(mParticlesLocalVelocityOld,i);
                    oldvelocity[2] = KRATOS_OCL_4_ARRAY_Z(mParticlesLocalVelocityOld,i);
                }
            }
        }

        Timer::Stop("PythonToLibTransfer");

        Timer::Start("Cpu-Time");
        Timer::Start("Reseed");

        if(reseed)
        {
            OpenMPUtils::SetNumThreads(num_threads);
            mMaxParticles = 7;
            mMinParticles = 1;

            for (ModelPart::ElementsContainerType::iterator el_it = mStaticMesh->ElementsBegin(); el_it != mStaticMesh->ElementsEnd(); el_it++)
            {
                el_it->SetValue(YOUNG_MODULUS,0.0);
            }

            #pragma omp parallel for
            for(size_t i = 0; i < localIndex; i++)
            {
                if (KRATOS_OCL_4_ARRAY_W(mParticlesLocal,i) >= 0)
                {
                    ModelPart::NodesContainerType::iterator inode = mParticMesh->NodesBegin() + KRATOS_OCL_4_ARRAY_W(mParticlesLocal,i);

                    if(!inode->Is(TO_ERASE))
                    {
                        if((mParticlesLocalI[i]) > 0 )
                        {
                            ModelPart::ElementsContainerType::iterator el_it = mStaticMesh->ElementsBegin() + (mParticlesLocalI[i]);
                            double& counter = el_it->GetValue(YOUNG_MODULUS);

                            #pragma omp critical
                            counter+=1;
                        }
                    }
                }
            }

            #pragma omp parallel for
            for(size_t i = 0; i < localIndex; i++)
            {
                if (KRATOS_OCL_4_ARRAY_W(mParticlesLocal,i) >= 0)
                {
                    ModelPart::NodesContainerType::iterator inode = mParticMesh->NodesBegin() + KRATOS_OCL_4_ARRAY_W(mParticlesLocal,i);

                    if(!inode->Is(TO_ERASE))
                    {
                        if((mParticlesLocalI[i]) > 0 )
                        {
                            ModelPart::ElementsContainerType::iterator el_it = mStaticMesh->ElementsBegin() + (mParticlesLocalI[i]);

                            if(el_it->GetValue(YOUNG_MODULUS) < mMinParticles)
                            {
                                inode->Set(TO_ERASE,true);
                                inode->X() += 30;
                            }
                            else if(el_it->GetValue(YOUNG_MODULUS) > mMaxParticles)
                            {
                                #pragma omp critical
                                el_it->GetValue(YOUNG_MODULUS)--;

                                inode->Set(TO_ERASE,true);
                                inode->X() += 30;
                            }
                        }
                    }
                }
            }

            num_threads = OpenMPUtils::GetNumThreads();
            OpenMPUtils::SetNumThreads(1);

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

                    while(no_it != mParticMesh->NodesEnd() && !no_it->Is(TO_ERASE)) no_it++;

                    for (unsigned int i = 0; i < pos.size1(); i++)
                    {
                        if (no_it != mParticMesh->NodesEnd())
                        {
                            no_it->Set(TO_ERASE,false);

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

                            while(no_it != mParticMesh->NodesEnd() && !no_it->Is(TO_ERASE)) no_it++;
                        }
                        else
                        {
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
        Timer::Stop("Reseed");

        Timer::Start("TransferToEulerianMesh");

        for(size_t i = 0; i < mNodeSize; i++)
        {
            if (!FixedV[i])
            {
                memset(&nodesV[i],0,sizeof(double)*4);
            }

            nodesY[i] = 0;
        }

        for(size_t i = 0, j = 0; i < localIndex && j < oldParticleNum;)
        {
            int tx, ty, tz;

            if (i < localIndex && mParticlesLocalI[i] != -1 && KRATOS_OCL_4_ARRAY_W(mParticlesLocal,i) != -2 )
            {
                if ((mParticMesh->NodesBegin()+KRATOS_OCL_4_ARRAY_W(mParticlesLocal,i))->Is(TO_ERASE) == false)
                {
                    int elementId = mParticlesLocalI[i];

                    tx = KRATOS_OCL_4_ARRAY_X(mTriangleNodes,(elementId))-1;
                    ty = KRATOS_OCL_4_ARRAY_Y(mTriangleNodes,(elementId))-1;
                    tz = KRATOS_OCL_4_ARRAY_Z(mTriangleNodes,(elementId))-1;

                    if(!FixedV[tx])
                    {
                        KRATOS_OCL_4_ARRAY_X(nodesV,tx) += KRATOS_OCL_4_ARRAY_X(mParticlesLocalN,i) * KRATOS_OCL_4_ARRAY_X(mParticlesLocalVelocity,i);
                        KRATOS_OCL_4_ARRAY_Y(nodesV,tx) += KRATOS_OCL_4_ARRAY_X(mParticlesLocalN,i) * KRATOS_OCL_4_ARRAY_Y(mParticlesLocalVelocity,i);
                        KRATOS_OCL_4_ARRAY_Z(nodesV,tx) += KRATOS_OCL_4_ARRAY_X(mParticlesLocalN,i) * KRATOS_OCL_4_ARRAY_Z(mParticlesLocalVelocity,i);

                        nodesY[tx] += KRATOS_OCL_4_ARRAY_X(mParticlesLocalN,i);
                    }

                    if(!FixedV[ty])
                    {
                        KRATOS_OCL_4_ARRAY_X(nodesV,ty) += KRATOS_OCL_4_ARRAY_Y(mParticlesLocalN,i) * KRATOS_OCL_4_ARRAY_X(mParticlesLocalVelocity,i);
                        KRATOS_OCL_4_ARRAY_Y(nodesV,ty) += KRATOS_OCL_4_ARRAY_Y(mParticlesLocalN,i) * KRATOS_OCL_4_ARRAY_Y(mParticlesLocalVelocity,i);
                        KRATOS_OCL_4_ARRAY_Z(nodesV,ty) += KRATOS_OCL_4_ARRAY_Y(mParticlesLocalN,i) * KRATOS_OCL_4_ARRAY_Z(mParticlesLocalVelocity,i);

                        nodesY[ty] += KRATOS_OCL_4_ARRAY_Y(mParticlesLocalN,i);
                    }

                    if(!FixedV[tz])
                    {
                        KRATOS_OCL_4_ARRAY_X(nodesV,tz) += KRATOS_OCL_4_ARRAY_Z(mParticlesLocalN,i) * KRATOS_OCL_4_ARRAY_X(mParticlesLocalVelocity,i);
                        KRATOS_OCL_4_ARRAY_Y(nodesV,tz) += KRATOS_OCL_4_ARRAY_Z(mParticlesLocalN,i) * KRATOS_OCL_4_ARRAY_Y(mParticlesLocalVelocity,i);
                        KRATOS_OCL_4_ARRAY_Z(nodesV,tz) += KRATOS_OCL_4_ARRAY_Z(mParticlesLocalN,i) * KRATOS_OCL_4_ARRAY_Z(mParticlesLocalVelocity,i);

                        nodesY[tz] += KRATOS_OCL_4_ARRAY_Z(mParticlesLocalN,i);
                    }
                }
            }

            i++;
            while(i < localIndex && KRATOS_OCL_4_ARRAY_W(mParticlesLocal,i) == -2) i++;
        }

        for(size_t i = 0; i < mNodeSize; i++)
        {
            if (!FixedV[i] && nodesY[i] != 0)
            {
                KRATOS_OCL_4_ARRAY_X(nodesV,i) /= nodesY[i];
                KRATOS_OCL_4_ARRAY_Y(nodesV,i) /= nodesY[i];
                KRATOS_OCL_4_ARRAY_Z(nodesV,i) /= nodesY[i];
            }
        }

        Timer::Stop("TransferToEulerianMesh");
        Timer::Stop("Cpu-Time");

        std::cout << mParticMesh->NumberOfNodes() << std::endl;

        OCLDeviceGroup.DeleteBuffer(OCL_IndexSelection);

        free(mIndexSelection);

        //Restore default number of ompThreads
        OpenMPUtils::SetNumThreads(num_threads);
    }

    /// Search nearest point
    void SearchNearestOCL(double Radius, unsigned int ConcurrentPoints)
    {
        unsigned int processed = 0;

        int result = 0;

        double HOST_memRadius  = Radius;
        double HOST_memRadius2 = Radius * Radius;

        unsigned int amount;

        while (processed < mProblemSize)
        {
            amount = (ConcurrentPoints > mProblemSize) ? mProblemSize : ((processed + ConcurrentPoints) < mProblemSize) ? ConcurrentPoints : mProblemSize - processed;

            int * resultPoints;

            resultPoints = (int *)malloc(sizeof(int) * amount);

            OCLDeviceGroup.CopyBuffer(OCL_Radius   , OpenCL::HostToDevice, OpenCL::VoidPList(1,&HOST_memRadius));
            OCLDeviceGroup.CopyBuffer(OCL_Radius2  , OpenCL::HostToDevice, OpenCL::VoidPList(1,&HOST_memRadius2));
            OCLDeviceGroup.CopyBuffer(OCL_Particles, OpenCL::HostToDevice, OCL_Map_Particles);

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
            OCLDeviceGroup.Synchronize();

            OCLDeviceGroup.CopyBuffer(OCL_results, OpenCL::DeviceToHost, OpenCL::VoidPList(1,resultPoints));

            result = resultPoints[amount-1];

            free(resultPoints);

            processed += amount;
        }

// 	    std::cout << "Nearest Point ID: " << result << std::endl;
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
    void PrintSize( std::ostream& rout )
    {
        rout << " BinsSize: ";
        for(SizeType i = 0 ; i < TDimension ; i++)
            rout << "[" << mN[i] << "]";
        rout << std::endl;
    }

    /// Print Limits Points of the Container
    void PrintBox( std::ostream& rout )
    {
        rout << " BinsBox: Min [";
        mMinPoint.Print(rout);
        rout <<       "];  Max [";
        mMaxPoint.Print(rout);
        rout <<       "];  Size [";
        mCellSize.Print(rout);
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

    OpenCL::DeviceGroup &OCLDeviceGroup;

    ModelPart * mStaticMesh;
    ModelPart * mParticMesh;

    cl_int Err;

    cl_uint OCL_program;

    /// OpenCL Kernels
    cl_uint OCLGenerateBinsA;
    cl_uint OCLGenerateBinsC;
    cl_uint OCLGenerateBinsCentA;
    cl_uint OCLGenerateBinsCentC;
    cl_uint OCLSearchInRadius;
    cl_uint OCLSearchNearest;
    cl_uint OCLSearch;
    cl_uint OCLMove;
    cl_uint OCLMove2;
    cl_uint OCLScan_Local1;
    cl_uint OCLScan_Local2;
    cl_uint OCLuniformUpdate;
    cl_uint OCLResetCounter;
    cl_uint OCLCalculateParticleIndex;
    cl_uint OCLProject;
    cl_uint OCLProject2;

    /// OpenCL Kernel memoryObjects
    cl_uint OCL_PointsBins;
    cl_uint OCL_IndexCell;
    cl_uint OCL_IndexCellReference;
    cl_uint OCL_IndexCellReferenceO;
    cl_uint OCL_IndexCellReferenceElementO;
    cl_uint OCL_IndexCellReferenceSize;
    cl_uint OCL_InvCellSize;
    cl_uint OCL_N;
    cl_uint OCL_MinPoint;
    cl_uint OCL_BinsObjectContainer;
    cl_uint OCL_Radius;
    cl_uint OCL_Radius2;
    cl_uint OCL_Scan_Buffer;
    cl_uint OCL_PointsTriangle;
    cl_uint OCL_TriangleList;
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
    cl_uint OCL_NodesV;
    cl_uint OCL_NodesY;
    cl_uint OCL_FixedV;
    cl_uint OCL_NodesF;
    cl_uint OCL_NodesP;
    cl_uint OCL_NodesR;
    cl_uint OCL_MacroBins;
    cl_uint OCL_IndexSelection;
    cl_uint OCL_IndexTriangles;

    cl_double4 * nodesV;
    int * FixedV;
    cl_double4 * nodesP;
    cl_double4 * nodesF;
    cl_double4 * nodesR;
    double     * nodesY;

    cl_double4 * PointsTriangles;

    uint mCellSizeInt;
    uint mBinsObjectSize;
    uint mBinsObjectCenterSize;
    uint mNodeSize;
    uint mTriangleSize;
    uint mProblemSize;
    uint mParticleBufferSize;
    uint mMinParticles;
    uint mMaxParticles;
    uint megacellDim;
    uint mMacroBinsSize;
    uint mMacroBinsNumElements;

    cl_double4 mBody_Force;

    int * BinsObjectContainer;
    int * BinsObjectContainerS;
    int * BinsObjectCenterContainer;
    int * mBinsObjectIndirection;
    int * mTriangleIndirection;
    std::vector<cl_int4> mMacroBins;

    std::vector<uint> IndexCellReferenceElement;
    std::vector<uint> IndexCellReferenceElementO;

    // MacroCell
    uint macroCellSize;
    uint level;

    cl_int4 ** mTriangleNodesMatrix;

    cl_double4 ** mParticlesMatrix;
    uint * mParticlesMatrixRowCount;

    int mIndexSelectionSize;
    int * mIndexSelection;

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
    cl_double4 * mParticlesDisplace;

    cl_double4 * mParticlesLocal;
    cl_double4 * mParticlesLocalVelocity;
    cl_double4 * mParticlesLocalVelocityOld;
    cl_double4 * mParticlesLocalDisplace;
    cl_double4 * mParticlesLocalForce;
    cl_double4 * mParticlesLocalN;
    int        * mParticlesLocalI;

    uint localBufferSize;

    int * mParticlesI;
    int * mParticlesIndex;
    int * mParticleMap;

    cl_int4 * mTriangleNodes;
    std::vector<std::list<cl_int4> > mTriangleNodesSorted;

    Tvector<SizeType,TDimension>  mN;

    IteratorVector    mIndexCell;
    IteratorVector    mIndexCellOCL;
    IteratorIterator  mIndexCellBegin;
    IteratorIterator  mIndexCellEnd;

    /// OpenCL Pinned memory stuff
    cl_uint OCL_Pinned_Particles;
    cl_uint OCL_Pinned_ParticlesVelocity;
    cl_uint OCL_Pinned_ParticlesVelocityOld;
    cl_uint OCL_Pinned_ParticlesDisplace;
    cl_uint OCL_Pinned_ParticlesDisplaceOld;

    OpenCL::VoidPList OCL_Map_Particles;
    OpenCL::VoidPList OCL_Map_ParticlesVelocity;
    OpenCL::VoidPList OCL_Map_ParticlesVelocityOld;
    OpenCL::VoidPList OCL_Map_ParticlesDisplace;
    OpenCL::VoidPList OCL_Map_ParticlesDisplaceOld;

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
