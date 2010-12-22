
#include <fstream>
#include <iostream>
#include <list>
#include <time.h>
#include <string>
#include <cstdlib>

#define KRATOS_INDEPENDENT
#include "spatial_containers/spatial_containers.h"
#include "bins_static_OCL.h"
#include "timer.h"

#include <omp.h>

double rrandom(){
   return double(rand())/RAND_MAX;
};

template< std::size_t dim_type>
class Point {

   public:

      double       coord[dim_type];
      std::size_t  id;
      std::size_t  tag;
      //int id;

      double& operator[](std::size_t i) {return coord[i];}

      double const & operator[](std::size_t i) const {return coord[i];}

      void RandomCoord(){
         for(std::size_t i = 0 ; i < dim_type ; i++)
            coord[i] = rrandom();
      }

      void operator=(Point<dim_type> const& Other){
         for(std::size_t i = 0; i < dim_type; i++)
            coord[i] = Other.coord[i];
      }
};

template< std::size_t dim_type >
std::ostream & operator<<( std::ostream& rOut, Point<dim_type> & rPoint){
	 rOut << "(" << rPoint.id << ") ";
   for(std::size_t i = 0 ; i < dim_type ; i++)
      rOut << rPoint[i] << " "; 
   return rOut; 
};

template< std::size_t dim_type >
std::istream & operator>>( std::istream& rIn, Point<dim_type> & rPoint){
   for(std::size_t i = 0 ; i < dim_type ; i++)
      rIn >> rPoint[i]; 
   return rIn; 
};

template< class T, std::size_t dim >
class PointDistance{
   public:
      double operator()( T const& p1, T const& p2 ){
         double dist = 0.0;
         for( std::size_t i = 0 ; i < dim ; i++){
            double tmp = p1[i] - p2[i];
            dist += tmp*tmp;
         }
         return sqrt(dist);
      }
};

template< class T, std::size_t dim >
class PointDistance2{
   public:
      double operator()( T const& p1, T const& p2 ){
         double dist = 0.0;
         for( std::size_t i = 0 ; i < dim ; i++){
            double tmp = p1[i] - p2[i];
            dist += tmp*tmp;
         }
         return dist;
      }
};

template< std::size_t dim >
bool LowerPoint( Point<dim> const& reference, Point<dim> const& new_ ){
   for(std::size_t i = 0 ; i < dim ; i++)
      if( reference[i] < new_[i] )
         return false;
   return true;
};

template< std::size_t dim >
bool UpperPoint( Point<dim> const& reference, Point<dim> const& new_ ){
   for(std::size_t i = 0 ; i < dim ; i++)
      if( reference[i] > new_[i] )
         return false;
   return true;
};



template< class TreeType, class PointType, class IteratorType, class DistanceIterator >
void RunTest( char const Title[], IteratorType PBegin, IteratorType PEnd, IteratorType Results, DistanceIterator Distances, std::size_t MaxResults, PointType& search_point, double Radius, std::size_t numsearch )
{

   timer time1;
   std::size_t n;
   PointType* PNearest;

   std::size_t numsearch_nearest = 100 * numsearch;

   time1.restart();
   TreeType   nodes_tree( PBegin, PEnd );
   std::cout << Title << " \t" << time1 << "\t\t";
   
   time1.restart();
   for(std::size_t i = 0 ; i < numsearch ; i++)
      n = nodes_tree.SearchInRadius( search_point, Radius, Results, Distances, MaxResults );
    std::cout << time1  << "\t\t\t";
   

   time1.restart();

   for(std::size_t i = 0 ; i < numsearch_nearest ; i++)
      PNearest = nodes_tree.SearchNearestPoint( search_point, Distances[0] );
   std::cout << time1 << "\t\t\t";

  std::cout  << Radius << "\t" << n << "\t\t" << *PNearest << std::endl;
};



int main(int arg, char* argv[])
{

   static const std::size_t Dim = 3;
   
   typedef Point<Dim> PointType;

   typedef PointType*                          PtrPointType;
   typedef PtrPointType*                       PointVector;
   typedef PtrPointType*                       PointIterator;

   typedef double*                             DistanceVector;
   typedef double*                             DistanceIterator;

   typedef Kratos::Bins< Dim, PointType, PointVector, PtrPointType, PointIterator, DistanceIterator, PointDistance2<PointType,Dim> > StaticBins;
   typedef Kratos::BinsOCL< Dim, PointType, PointVector, PtrPointType, PointIterator, DistanceIterator, PointDistance2<PointType,Dim> > StaticBinsOCL;
   typedef Kratos::BinsDynamic< Dim, PointType, PointVector, PtrPointType, PointIterator, DistanceIterator, PointDistance2<PointType,Dim> > DynamicBins;

   typedef Kratos::Bucket<Dim,PointType,PointVector, PtrPointType, PointIterator, DistanceIterator, PointDistance2<PointType,Dim> > BucketType;
   typedef Kratos::Tree< Kratos::KDTreePartition<BucketType> > kd_tree;
   typedef Kratos::Tree< Kratos::KDTreePartitionAverageSplit<BucketType> > kd_tree_aversplit;
   typedef Kratos::Tree< Kratos::KDTreePartitionMidPointSplit<BucketType> > kd_tree_midpoint;
  
   typedef Kratos::Tree< Kratos::OCTreePartition<BucketType> > oct_tree;

   typedef Kratos::Tree< StaticBins > bins_tree;
   /*    typedef Kratos::Tree< DynamicBins> binsdyn_tree; */
  
   typedef Kratos::Tree< Kratos::KDTreePartitionAverageSplit<StaticBins>  > kdtree_bins_tree;
   /*    typedef Kratos::Tree< Kratos::KDTreePartition<DynamicBins> > kdtree_binsdyn_tree; */



  // set format
  //std::cout.setf(0,std::ios::floatfield);
  std::cout.precision(5);

  //std::cout << std::fixed;
  //std::cout << std::setprecision(5);
   
   PointVector points;
   std::string filename;
   
   std::list<cl_double4> TPoints;
   std::list<cl_int4> TTriangles;
   
   double radius  = 16.00;
   int rep = 10;
   int block = 1000;
   int maxresults = 10;
  
   if(arg == 1){
      std::cout << " argument not founded " << std::endl;
      return 0;
   }

   else if(arg == 3){
      radius = atof(argv[2]);
   }

   else if(arg == 4){
      radius = atof(argv[2]);
      rep = atof(argv[3]);
   }
   
   else if(arg == 5){
      radius = atof(argv[2]);
      rep = atof(argv[3]);
      block = atof(argv[4]);
   }
   
   else if(arg == 6){
      radius = atof(argv[2]);
      rep = atof(argv[3]);
      block = atof(argv[4]);
      maxresults = atof(argv[5]);
   }
   
   else {
      std::cout << "Invalid num of params: ./Filename (input file) [radius] [num of reps] [blockSize] [maxresults]" << std::endl;
      return 0;
   }

   filename = argv[1];

   std::ifstream input;
   input.open(filename.c_str());
   if (!input) {
      std::cout << "Cannot open data file" << std::endl;
      return 0;
   }
   
   PointType point;
   std::size_t npoints;
   std::size_t nTriangles;
   cl_double4 pointdata;
   cl_int4 indexdata;
   double radius3 = radius;
   std::string buffer;
   npoints = -1;
   
   //No me interesa todo eso ahora.
   for(int i = 0; i < 14; i++)
        input >> buffer;
   
  //Load point data into a buffer
  input >> buffer;
  for(int i = 1;buffer.find("End") == -1;i++)
  {
     input >> pointdata.x;
     input >> pointdata.y;
     input >> pointdata.z;
     pointdata.w = i;
     
     TPoints.push_back(pointdata);
     
     input >> buffer;
     //std::cout << "Data: ("<< pointdata.w << ") " << pointdata.x << " " << pointdata.y << " " << pointdata.z << " " << std::endl;
   }
   
   npoints = TPoints.size();

   std::cout << "ALL POINTS LOADED" << std::endl;
   
   //No me interesa todo eso ahora.
   for(int i = 0; i < 4; i++)
      input >> buffer;
   
   std::cout << buffer << std::endl;
   int mesh2d = buffer.find("2D");
   std::cout << "2D Mesh?: " << mesh2d << std::endl;
      
   //Load index data in to another vector
   input >> buffer;
   
   if(mesh2d != -1)
   {
      std::cout << "Loading 2D mesh: " << std::endl;
      for(int i = 1;buffer.find("End") == -1;i++)
      {
	  //Skip 0 value
	  input >> buffer;
	
	  //Load data
	  input >> indexdata.x;
	  input >> indexdata.y;
	  input >> indexdata.z;
	  indexdata.w = -1;
	  
	  //Tetrahedron 4th index
	  //input >> pointdata.w;
	  
	  TTriangles.push_back(indexdata);
	  
	  input >> buffer;
	  //std::cout << "TRIANGLE: ("<< pointdata.w << ") " << pointdata.x << " " << pointdata.y << " " << pointdata.z << " " << std::endl;
      }
   }
   else
   {
      std::cout << "Loading 3D mesh: " << std::endl;
      for(int i = 1;buffer.find("End") == -1;i++)
      {
	  //Skip 0 value
	  input >> buffer;
	
	  //Load data
	  input >> indexdata.x;
	  input >> indexdata.y;
	  input >> indexdata.z;
	  input >> indexdata.w;
	  
	  TTriangles.push_back(indexdata);
	  
	  input >> buffer;
	  //std::cout << "TRIANGLE: " << indexdata.x << " " << indexdata.y << " " << indexdata.z << " " << indexdata.w << std::endl;
      }
   }
   
   nTriangles = TTriangles.size();
   
   std::cout << "ALL TRIANGLES LOADED" << std::endl;
   
   PointVector PointsArray = new PointType* [TPoints.size()];
   
   //Provar a tratar un point como cl_int4 para pasar los 4 indicies de golpe. guarrrada , preguntar.
   cl_int4 * IndexsArray = (cl_int4 *)malloc(sizeof(cl_int4)*TTriangles.size());
   
   for(int i = 0; !TPoints.empty(); i++)
   {
      PointType auxPoint;
      cl_double4 pointdata = TPoints.front();
      auxPoint[0] = pointdata.x;
      auxPoint[1] = pointdata.y;
      auxPoint[2] = pointdata.z;
      auxPoint.id = pointdata.w;
      PointsArray[i] = new PointType(auxPoint);
      TPoints.pop_front();
   }
   
   std::cout << "ALL PONTS FILLED" << std::endl;
   
   //Asi se guardan de 4 en 4 i el ID coincide con el indice.
   for(int i = 0; !TTriangles.empty(); i++)
   {
      cl_int4 data = TTriangles.front();
      IndexsArray[i] = TTriangles.front();
      TTriangles.pop_front();
   }
   
   std::cout << "ALL INDEX FILLED" << std::endl;
   
   std::cout << "ALL ELEMENTS PROCESED AND LOADED INTO A CL_DOUBLE4 ARRAY" << std::endl;
   std::cout << "Check First and last elements..." << std::endl;
   
   struct timespec begin;
   struct timespec end;
   int sampleSize = 1000;
   
   clock_gettime( CLOCK_REALTIME, &begin );
   StaticBinsOCL binOCL(PointsArray,PointsArray+npoints,IndexsArray,nTriangles);
   clock_gettime( CLOCK_REALTIME, &end );
   
   std::cout << "Init bins:\t\t" << ((float)(end.tv_sec - begin.tv_sec) + (float)(end.tv_nsec-begin.tv_nsec)/1000000000) << std::endl;
   
   clock_gettime( CLOCK_REALTIME, &begin );
   binOCL.GenerateSampleInput(sampleSize);
   clock_gettime( CLOCK_REALTIME, &end );
   
   std::cout << "GenerateSample of: " << sampleSize * sampleSize << " elements:\t\t" << ((float)(end.tv_sec - begin.tv_sec) + (float)(end.tv_nsec-begin.tv_nsec)/1000000000) << std::endl;
   
   clock_gettime( CLOCK_REALTIME, &begin );
   binOCL.generateBins();
   clock_gettime( CLOCK_REALTIME, &end );
   
   std::cout << "Generate bins:\t\t" << ((float)(end.tv_sec - begin.tv_sec) + (float)(end.tv_nsec-begin.tv_nsec)/1000000000) << std::endl;
   
   clock_gettime( CLOCK_REALTIME, &begin );   
   binOCL.allocateOCLBuffers(block,maxresults);
   clock_gettime( CLOCK_REALTIME, &end );
   
   std::cout << "Allocate buffer:\t\t" << ((float)(end.tv_sec - begin.tv_sec) + (float)(end.tv_nsec-begin.tv_nsec)/1000000000) << std::endl;
   
   clock_gettime( CLOCK_REALTIME, &begin );
   for(int i = 0; i < rep; i++)
      binOCL.searchTriangles(radius3,block,maxresults);
   clock_gettime( CLOCK_REALTIME, &end );
   
   std::cout << "Perform Search:\t\t" << ((float)(end.tv_sec - begin.tv_sec) + (float)(end.tv_nsec-begin.tv_nsec)/1000000000) << std::endl;
   
   /*
   points = new PointType*[npoints];
   for(std::size_t i = 0; i < npoints; i++){
      input >> point;
      //point.RandomCoord();
      point.id = i+1;
      //points.push_back(new PointType(point));
      points[i] = new PointType(point);
   }

   PointType min_point(*points[0]);
   PointType max_point(*points[0]);
   PointType mid_point;
   for(std::size_t i = 0; i < npoints; i++){
      if( LowerPoint(min_point,*points[i]) ) min_point = *points[i];
      if( UpperPoint(max_point,*points[i]) ) max_point = *points[i];
   }

   for(std::size_t i = 0 ; i < Dim ; i++)
      mid_point.coord[i] = ( max_point[i] + min_point[i] ) / 2.00;
   //    PointType*     search_point = &mid_point;
   std::cout << std::endl;
   std::cout << "BoundingBox.min_point : " << min_point << std::endl;
   std::cout << "BoundingBox.max_point : " << max_point << std::endl;
   //    std::cout << " search_point : " << *search_point << std::endl;


   std::size_t maxloops = 1;
   std::size_t numsearch = std::min(100000,int(npoints));
   std::size_t numsearch_nearest = std::min(int(numsearch*100),int(npoints));
  
   std::cout << std::endl;
   std::cout << " Number of Points : " << npoints << std::endl;
   std::cout << " Number of Repetitions : " << numsearch << std::endl;
   std::cout << std::endl;

   timer time1;
   //PointIterator SearchPoint;
   
   PtrPointType  PNearest = NULL;
   PointVector results;
   results = new PointType*[npoints];
   DistanceVector distance;
   distance = new double[npoints];
   
   std::size_t max_results = npoints;
   std::size_t n = 0;
   double radius3 = radius;
   PointVector points1 = new PointType*[npoints];
   double AverDistance;
   std::size_t n_res;
   */
  
}

