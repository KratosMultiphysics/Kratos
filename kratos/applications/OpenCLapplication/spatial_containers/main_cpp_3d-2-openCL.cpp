
#include <fstream>
#include <iostream>
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
   
   double radius  = 16.00;
   int rep = 10;
   int block = 1000;
   int maxresults = 10;
  
   if(arg == 1){
      std::cout << " argument not founded " << std::endl;
      return 0;
   }

   if(arg == 3){
      radius = atof(argv[2]);
   }

   if(arg == 4){
      radius = atof(argv[2]);
      rep = atof(argv[3]);
   }
   
   if(arg == 5){
      radius = atof(argv[2]);
      rep = atof(argv[3]);
      block = atof(argv[4]);
   }
   
   if(arg == 6){
      radius = atof(argv[2]);
      rep = atof(argv[3]);
      block = atof(argv[4]);
      maxresults = atof(argv[5]);
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
   //int npoints = atoi(argv[1]);
   input >> npoints;

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
   /*    PointType*     search_point = &mid_point; */
   std::cout << std::endl;
   std::cout << "BoundingBox.min_point : " << min_point << std::endl;
   std::cout << "BoundingBox.max_point : " << max_point << std::endl;
   /*    std::cout << " search_point : " << *search_point << std::endl; */


   std::size_t maxloops = 1;
   std::size_t numsearch = std::min(100000,int(npoints));
   std::size_t numsearch_nearest = std::min(int(numsearch*100),int(npoints));
  
   std::cout << std::endl;
   std::cout << " Number of Points : " << npoints << std::endl;
   std::cout << " Number of Repetitions : " << numsearch << std::endl;
   std::cout << std::endl;

   timer time1;
   //PointIterator SearchPoint;

   std::cout << "\t\tGeneration\tsearch in radius\tsearch nearest point\tRadius\tNo. of Results\tNearestPoint" << std::endl;
   
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
   
  //**********************************************************************************
  /*CAPADO!!!
   time1.restart();
   for(std::size_t i = 0; i < npoints; i++)
     points1[i] = points[i];
   kd_tree_midpoint  nodes_tree_mp(points1,points1+npoints, 64 );
   //kd_tree  nodes_tree2(points1,points1+npoints, 2 );
   std::cout << "KdTreeMidPoint\t" << time1 << "\t\t";
    
   time1.restart();
   //CHARLIE -- for(int i = 0; i < 1200; i++) {
   n = 0;
   AverDistance = 0.00;
   for(std::size_t i = 0 ; i < npoints ; i++){
     n_res = nodes_tree_mp.SearchInRadius( *points[i], radius3, results, distance, max_results);
     n += n_res;
     for(std::size_t j = 0; j < n_res; j++)
       AverDistance += distance[j];
   }
   //CHARLIE }
   std::cout <<  time1  << "\t\t\t";

   time1.restart();
   for(std::size_t j = 0 ; j < maxloops ; j++){
     for(std::size_t i = 0 ; i < numsearch_nearest ; i++)
       PNearest = nodes_tree_mp.SearchNearestPoint( *points[i], distance[0] );
   }
   std::cout << time1 << "\t\t\t";

  std::cout  << radius << "\t" << n/npoints << "\t\t" << *PNearest << std::endl;
  std::cout  << "\t\t Total Distances = " << AverDistance << "(" << n << ")" << std::endl;
   */
   //**********************************************************************************

   time1.restart();
   //BinsContainer bin(points.begin(), points.end());
   for(std::size_t i = 0; i < npoints; i++)
     points1[i] = points[i];
   
   //Original
   StaticBins bin(points1,points1+npoints);
      
   //Ocl
   StaticBinsOCL binOCL(points1,points1+npoints);
      
   std::cout << "Bins\t\t" << time1 << "\t\t" << std::endl;
   time1.restart();
   n = 0;
   AverDistance = 0.00;

   
   
   for(std::size_t i = 0 ; i < npoints ; i++)
   {
     n_res = bin.SearchInRadius( *points[i], radius3, results, distance, max_results);
     n += n_res;
     for(std::size_t j = 0; j < n_res; j++)
       AverDistance += distance[j];
   }
   std::cout <<  time1  << "\t\t\t";
   
   time1.restart();/*
   for(std::size_t j = 0 ; j < maxloops ; j++){
     for(std::size_t i = 0 ; i < numsearch_nearest ; i++)
       PNearest = bin.SearchNearestPoint( *points[i], distance[0] );
   }*/
   std::cout << time1 << "\t\t\t";



   time1.restart();

   //RICCARDOs OPENMP

	











   //std::cout  << std::endl;
   std::cout  << radius << "\t" << n << "\t" << npoints << "\t\t" << 0/**PNearest*/ << std::endl;
   std::cout  << "\t\t Total Distances = " << AverDistance << "(" << n << ")" << std::endl;
  
   time1.restart();













   int omp_threads = omp_get_max_threads();
   
   // PointVector results = new PointType*[npoints];
   // DistanceVector distance = new double[npoints];
   
   DistanceVector *parallel_distance = new double *[omp_threads];
   PointVector *parallel_results = new PointType **[ omp_threads];
   #pragma omp paralle for
   for ( int i = 0; i < omp_threads; i++) {
     parallel_distance[i] = new double[ npoints]; 
     parallel_results[ i] = new PointType *[ npoints];
   }

    double t0 = omp_get_wtime();
   #pragma omp parallel for firstprivate(radius3,max_results) 
   for(std::size_t i = 0 ; i < npoints ; i++)
   {
     int my_id = omp_get_thread_num();
     n_res = bin.SearchInRadius( *points[i], radius3, parallel_results[my_id], parallel_distance[my_id], max_results);
/*     n += n_res;
     for(std::size_t j = 0; j < n_res; j++)
       AverDistance += distance[j];*/
   }
   std::cout << "parallel for time" <<  omp_get_wtime() - t0  << "\t\t\t";
   
   time1.restart();














   
    struct timespec begin;
    struct timespec end;


    std::cout << "\nTesting parallel point calculaton" << std::endl;
    std::cout << "Block\t\t" << "Reps\t\t" <<  std::endl;
    std::cout << block << "\t\t" << rep << std::endl;

    clock_gettime( CLOCK_REALTIME, &begin );
    for(std::size_t i = 0 ; i < npoints ; i++)
    {
	binOCL.prepareData(*points[i]);
    }  
   
    binOCL.allocateOCLBuffers(block,maxresults);
    
    for(int i = 0; i < rep; i++) 
    {
        binOCL.computeresultsN(radius3,block,maxresults);
    }
    
   clock_gettime( CLOCK_REALTIME, &end );
   
      std::cout << "BEGIN: " << begin.tv_sec << "." << begin.tv_nsec << "\n"
        << "END:    " << end.tv_sec << "." << end.tv_nsec << std::endl;
	std::cout << "TOTAL: " << ((float)(end.tv_sec - begin.tv_sec) + (float)(end.tv_nsec-begin.tv_nsec)/1000000000) << std::endl;

   //**********************************************************************************
   /*CAPADO!!!!
   //return 0;

   
   time1.restart();
   for(std::size_t i = 0; i < npoints; i++)
     points1[i] = points[i];
   DynamicBins nodes_treeDyn( points1, points1+npoints );
   std::cout << "Dynamic Bin\t" << time1 << "\t\t";
   
   time1.restart();
   n = 0;
   for(std::size_t i = 0 ; i < npoints ; i++){
     n += nodes_treeDyn.SearchInRadius( *points[i], radius3, results, distance, max_results);
   }
   std::cout <<  time1  << "\t\t\t";
   

   time1.restart();
   for(std::size_t j = 0 ; j < maxloops ; j++){
     for(std::size_t i = 0 ; i < numsearch_nearest ; i++)
       PNearest = nodes_treeDyn.SearchNearestPoint( *points[i], distance[0] );
   }
   std::cout << time1 << "\t\t\t";

  std::cout  << radius << "\t" << n/npoints << "\t\t" << *PNearest << std::endl;
 */
  /*    std::cout <<  "npoints = " << npoints <<  "  total results = " << n << std::endl; */
   //**********************************************************************************
   
   //return 0;
  
//CAPADO
/*
   time1.restart();
   for(std::size_t i = 0; i < npoints; i++)
     points1[i] = points[i];
   kd_tree  nodes_tree2(points1,points1+npoints, 64 );
   //kd_tree  nodes_tree2(points1,points1+npoints, 2 );
   std::cout << "KdTree\t\t" << time1 << "\t\t";
    
   time1.restart();
   n = 0;
   for(std::size_t i = 0 ; i < npoints ; i++)
     n += nodes_tree2.SearchInRadius( *points[i], radius3, results, distance, max_results);
   std::cout <<  time1  << "\t\t\t";

   time1.restart();
   for(std::size_t j = 0 ; j < maxloops ; j++){
     for(std::size_t i = 0 ; i < numsearch_nearest ; i++)
       PNearest = nodes_tree2.SearchNearestPoint( *points[i], distance[0] );
   }
   std::cout << time1 << "\t\t\t";

  std::cout  << radius << "\t" << n/npoints << "\t\t" << *PNearest << std::endl;
   
  //**********************************************************************************

   time1.restart();
   for(std::size_t i = 0; i < npoints; i++)
     points1[i] = points[i];
   kd_tree_aversplit  nodes_tree_av(points1,points1+npoints, 64 );
   //kd_tree  nodes_tree2(points1,points1+npoints, 2 );
   std::cout << "KdTreeAverSplit\t" << time1 << "\t\t";
    
   time1.restart();
   n = 0;
   for(std::size_t i = 0 ; i < npoints ; i++)
     n += nodes_tree_av.SearchInRadius( *points[i], radius3, results, distance, max_results);
   std::cout <<  time1  << "\t\t\t";

   time1.restart();
   for(std::size_t j = 0 ; j < maxloops ; j++){
     for(std::size_t i = 0 ; i < numsearch_nearest ; i++)
       PNearest = nodes_tree_av.SearchNearestPoint( *points[i], distance[0] );
   }
   std::cout << time1 << "\t\t\t";

  std::cout  << radius << "\t" << n/npoints << "\t\t" << *PNearest << std::endl;
   
  //**********************************************************************************

   time1.restart();
   for(std::size_t i = 0; i < npoints; i++)
     points1[i] = points[i];
   oct_tree  nodes_tree_oct(points1,points1+npoints, 64 );
   std::cout << "OctTree\t\t" << time1 << "\t\t";
    
   time1.restart();
   n = 0;
   for(std::size_t i = 0 ; i < npoints ; i++)
     n += nodes_tree_oct.SearchInRadius( *points[i], radius3, results, distance, max_results);
   std::cout <<  time1  << "\t\t\t";

   time1.restart();
   for(std::size_t j = 0 ; j < maxloops ; j++){
     for(std::size_t i = 0 ; i < numsearch_nearest ; i++)
       PNearest = nodes_tree_oct.SearchNearestPoint( *points[i], distance[0] );
   }
   std::cout << time1 << "\t\t\t";

  std::cout  << radius << "\t" << n/npoints << "\t\t" << *PNearest << std::endl;
   
  //**********************************************************************************
   
   //return 0;

   time1.restart();
   for(std::size_t i = 0; i < npoints; i++)
     points1[i] = points[i];
   //bins_tree  nodes_bintree2(points1,points1+npoints, std::size_t(npoints/2)+1 );
   kdtree_bins_tree  nodes_bintree2( points1, points1+npoints, kdtree_bins_tree::Partitions(2) );
   std::cout << "BinsKdTree2\t" << time1 << "\t\t";
    
   time1.restart();
   n = 0;
   for(std::size_t i = 0 ; i < npoints ; i++)
     n += nodes_bintree2.SearchInRadius( *points[i], radius3, results, distance, max_results);
   std::cout <<  time1  << "\t\t\t";

   time1.restart();
   for(std::size_t j = 0 ; j < maxloops ; j++){
     for(std::size_t i = 0 ; i < numsearch_nearest ; i++)
       PNearest = nodes_bintree2.SearchNearestPoint( *points[i], distance[0] );
   }
   std::cout << time1 << "\t\t\t";

  std::cout  << radius << "\t" << n/npoints << "\t\t" << *PNearest << std::endl;
   
  //**********************************************************************************
   
   //return 0;

   time1.restart();
   for(std::size_t i = 0; i < npoints; i++)
     points1[i] = points[i];
   kdtree_bins_tree  nodes_bintree4(points1,points1+npoints, kdtree_bins_tree::Partitions(4) );
   //kd_tree  nodes_tree2(points1,points1+npoints, 2 );
   std::cout << "BinsKdTree4\t" << time1 << "\t\t";
    
   time1.restart();
   n = 0;
   for(std::size_t i = 0 ; i < npoints ; i++)
     n += nodes_bintree4.SearchInRadius( *points[i], radius3, results, distance, max_results);
   std::cout <<  time1  << "\t\t\t";

   time1.restart();
   for(std::size_t j = 0 ; j < maxloops ; j++){
     for(std::size_t i = 0 ; i < numsearch_nearest ; i++)
       PNearest = nodes_bintree4.SearchNearestPoint( *points[i], distance[0] );
   }
   std::cout << time1 << "\t\t\t";

  std::cout  << radius << "\t" << n/npoints << "\t\t" << *PNearest << std::endl;
   
  //**********************************************************************************
   
   //return 0;

   time1.restart();
   for(std::size_t i = 0; i < npoints; i++)
     points1[i] = points[i];
   kdtree_bins_tree  nodes_bintree8(points1,points1+npoints, kdtree_bins_tree::Partitions(8) );
   //kd_tree  nodes_tree2(points1,points1+npoints, 2 );
   std::cout << "BinsKdTree8\t" << time1 << "\t\t";
    
   time1.restart();
   n = 0;
   for(std::size_t i = 0 ; i < npoints ; i++)
     n += nodes_bintree8.SearchInRadius( *points[i], radius3, results, distance, max_results);
   std::cout <<  time1  << "\t\t\t";

   time1.restart();
   for(std::size_t j = 0 ; j < maxloops ; j++){
     for(std::size_t i = 0 ; i < numsearch_nearest ; i++)
       PNearest = nodes_bintree8.SearchNearestPoint( *points[i], distance[0] );
   }
   std::cout << time1 << "\t\t\t";

  std::cout  << radius << "\t" << n/npoints << "\t\t" << *PNearest << std::endl;
   
  //**********************************************************************************
   
   //return 0;

   time1.restart();
   for(std::size_t i = 0; i < npoints; i++)
     points1[i] = points[i];
   kdtree_bins_tree  nodes_bintree16(points1,points1+npoints, kdtree_bins_tree::Partitions(16) );
   //kd_tree  nodes_tree2(points1,points1+npoints, 2 );
   std::cout << "BinsKdTree16\t" << time1 << "\t\t";
    
   time1.restart();
   n = 0;
   for(std::size_t i = 0 ; i < npoints ; i++)
     n += nodes_bintree16.SearchInRadius( *points[i], radius3, results, distance, max_results);
   std::cout <<  time1  << "\t\t\t";

   time1.restart();
   for(std::size_t j = 0 ; j < maxloops ; j++){
     for(std::size_t i = 0 ; i < numsearch_nearest ; i++)
       PNearest = nodes_bintree16.SearchNearestPoint( *points[i], distance[0] );
   }
   std::cout << time1 << "\t\t\t";

  std::cout  << radius << "\t" << n/npoints << "\t\t" << *PNearest << std::endl;
   
  //**********************************************************************************
   
   //return 0;

   time1.restart();
   for(std::size_t i = 0; i < npoints; i++)
     points1[i] = points[i];
   kdtree_bins_tree  nodes_bintree32(points1,points1+npoints, kdtree_bins_tree::Partitions(32) );
   //kd_tree  nodes_tree2(points1,points1+npoints, 2 );
   std::cout << "BinsKdTree32\t" << time1 << "\t\t";
    
   time1.restart();
   n = 0;
   for(std::size_t i = 0 ; i < npoints ; i++)
     n += nodes_bintree32.SearchInRadius( *points[i], radius3, results, distance, max_results);
   std::cout <<  time1  << "\t\t\t";

   time1.restart();
   for(std::size_t j = 0 ; j < maxloops ; j++){
     for(std::size_t i = 0 ; i < numsearch_nearest ; i++)
       PNearest = nodes_bintree32.SearchNearestPoint( *points[i], distance[0] );
   }
   std::cout << time1 << "\t\t\t";

  std::cout  << radius << "\t" << n/npoints << "\t\t" << *PNearest << std::endl;
   
  //**********************************************************************************
   
   //return 0;

   time1.restart();
   for(std::size_t i = 0; i < npoints; i++)
     points1[i] = points[i];
   kdtree_bins_tree  nodes_bintree64(points1,points1+npoints, kdtree_bins_tree::Partitions(64) );
   //kd_tree  nodes_tree2(points1,points1+npoints, 2 );
   std::cout << "BinsKdTree64\t" << time1 << "\t\t";
    
   time1.restart();
   n = 0;
   for(std::size_t i = 0 ; i < npoints ; i++)
     n += nodes_bintree64.SearchInRadius( *points[i], radius3, results, distance, max_results);
   std::cout <<  time1  << "\t\t\t";

   time1.restart();
   for(std::size_t j = 0 ; j < maxloops ; j++){
     for(std::size_t i = 0 ; i < numsearch_nearest ; i++)
       PNearest = nodes_bintree64.SearchNearestPoint( *points[i], distance[0] );
   }
   std::cout << time1 << "\t\t\t";

  std::cout  << radius << "\t" << n/npoints << "\t\t" << *PNearest << std::endl;
   
  //**********************************************************************************
  */
  return 0;
/*
   time1.restart();
   int NN = int(pow(double(npoints),1.00/Dim));
   double bucketsize = (max_point[0] - min_point[0] ) / double(NN);
   DynamicBins Tree33( min_point, max_point, bucketsize );
   PtrPointType res33;
   for(PointIterator pp = points; pp != points+npoints; pp++)
   {
     if(!(res33=Tree33.ExistPoint(*pp)))
       Tree33.AddPoint(*pp);
     else
       std::cout << "Duplicated node!!" << std::endl;
   }
   if(!Tree33.ExistPoint(*points))
     Tree33.AddPoint(*points);
   else
     std::cout << "Duplicated node!!" << std::endl;

   std::cout << "Dynamic Bin2\t" << time1 << "\t\t";
   
   time1.restart();
   for(std::size_t j = 0 ; j < maxloops ; j++){
     SearchPoint = points;
     for(std::size_t i = 0 ; i < numsearch ; i++){
       n = Tree33.SearchInRadius( *points[i], radius3, results, distance, max_results);
       SearchPoint++;
     }
   }
   std::cout <<  time1  << "\t\t\t";

   time1.restart();
   for(std::size_t j = 0 ; j < maxloops ; j++){
     SearchPoint = points;
     for(std::size_t i = 0 ; i < numsearch_nearest ; i++)
       PNearest = Tree33.SearchNearestPoint( *points[i], distance[0] );
   }
   std::cout << time1 << "\t\t\t";

  std::cout  << radius << "\t" << n/npoints << "\t\t" << *PNearest << std::endl;
*/
   //**********************************************************************************

   // delete points
   for(PointIterator i_point = points ; i_point != points+npoints ; i_point++)
      delete *i_point;
   
   return 0;
}

