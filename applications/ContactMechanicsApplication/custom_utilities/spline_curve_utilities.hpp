//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:              December 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

#if !defined(KRATOS_SPLINE_CURVE_UTILITIES_H_INCLUDED )
#define  KRATOS_SPLINE_CURVE_UTILITIES_H_INCLUDED


// System includes
#include <cmath>
#include <set>

#ifdef _OPENMP
#include <omp.h>
#endif

// External includes
#include "boost/smart_ptr.hpp"
#include "utilities/timer.h"

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"
#include "utilities/math_utils.h"
#include "spatial_containers/spatial_containers.h"
#include "contact_mechanics_application_variables.h"

namespace Kratos
{

  /// Short class definition.

  /** Computes the energy
   */

  class KRATOS_API(CONTACT_MECHANICS_APPLICATION) SplineCurveUtilities
  {
  public:

    ///@name Type Definitions
    ///@{

    typedef ModelPart::NodesContainerType          NodesContainerType;

    //definitions for spatial search
    typedef array_1d<double, 3>                             PointType;
    typedef Node<3>                                          NodeType;
    typedef NodeType::Pointer                         NodePointerType;
    typedef std::vector<NodePointerType>        NodePointerTypeVector;
    typedef std::vector<NodeType>                      NodeTypeVector;
    typedef NodePointerTypeVector::iterator       NodePointerIterator;
    typedef std::vector<double>                        DistanceVector;
    typedef std::vector<double>::iterator            DistanceIterator;
    typedef Bucket<3, NodeType, NodePointerTypeVector, NodePointerType, NodePointerIterator, DistanceIterator > BucketType;
    typedef Tree< KDTreePartition<BucketType> >            KdtreeType; //Kdtree

    //structure for a spline segment type
    typedef struct
    {
      int        id;       // Spline Knot number
      PointType  P0;       // First auxiliar point
      PointType  P1;       // First curve point
      PointType  P2;       // Second curve point
      PointType  P3;       // Second auxiliar point

    } SplineType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SplineCurveUtilities(){ mEchoLevel = 0;  mParallel = true; };

    SplineCurveUtilities(bool Parallel){ mEchoLevel = 0;  mParallel = Parallel; };

    /// Destructor.
    ~SplineCurveUtilities(){};


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    //************************************************************************************
    //************************************************************************************

    // Create Arch Length Parametrized Spline Curve with m cubic segments
    void CreateParametrizedCurve(NodePointerTypeVector& rGeneratrixPoints, NodePointerTypeVector& rKnotsList, int m)
    {

      //Definition of the Spline Curve Q(t): (Spline with rGeneratrixPoints.size() cubic segments with the different arch length)

      //the numeration and order of the GeneratrixPoints is very important for the interpolation
      unsigned int id = 0; //start with 0;

      //Set auxiliary control points (knots) using reflection at the begining ( 1 extra knot )
      NodesContainerType::iterator nodes_begin = rGeneratrixPoints.begin();

      double X = 2.0 * nodes_begin->X() - (nodes_begin + 1)->X();
      double Y = 2.0 * nodes_begin->Y() - (nodes_begin + 1)->Y();
      double Z = 2.0 * nodes_begin->Z() - (nodes_begin + 1)->Z();

      NodePointerType KnotPoint = NodePointerType( new NodeType(id,X,Y,Z) );
      rKnotsList.push_back( KnotPoint );
      id+=1;

      //Set generatrix control points (knots) from the initial control stations
      for(NodePointerTypeVector::iterator in = rGeneratrixPoints.begin(); in!=rGeneratrixPoints.end(); in++)
	{
	  KnotPoint = NodePointerType( new NodeType(id,(*in)->X(),(*in)->Y(),(*in)->Z()) );
	  rKnotsList.push_back( KnotPoint );
	  id++;
	}


      //Set auxiliary control points (knots) using reflection at the end ( 1 extra knots )
      NodesContainerType::iterator nodes_end = rGeneratrixPoints.end()-1;

      X = 2.0 * (nodes_end)->X() - (nodes_end - 1)->X();
      Y = 2.0 * (nodes_end)->Y() - (nodes_end - 1)->Y();
      Z = 2.0 * (nodes_end)->Z() - (nodes_end - 1)->Z();

      KnotPoint = NodePointerType( new NodeType(id,X,Y,Z) );
      rKnotsList.push_back( KnotPoint );


      //Definition of the Parametrized Spline Curve Q(s):

      //Compute an approximate Arch-Length Parametrized Curve

      double TotalLength = 0;
      std::vector<double> SegmentArchLengths;

      SplineType Spline;
      double S = 0;
      int size = rKnotsList.size()-2;

      //std::cout<<" knots "<<size<<std::endl;

      for(int i=1; i<size; i++)
	{

	  SetSpline(Spline, rKnotsList, i);

	  //Get Segment Arch-Length by the adaptative gaussian integration method
	  S = AdaptiveIntegration(Spline);

	  //std::cout<<" SegmentArchLength "<<S<<std::endl;

	  //Set Segment Length
	  SegmentArchLengths.push_back(S);

	  //Add Segment to the Total Arch Length
	  TotalLength += S;

	}

      int n_segments = SegmentArchLengths.size();
      //std::cout<<" TotalArchLength "<<TotalLength<<" Segments "<<n_segments<<std::endl;


      //Find m+1 equally spaced points along Q(s)
      double Length = TotalLength / double(m);

      NodePointerTypeVector NewKnotsList;

      nodes_begin = rKnotsList.begin();

      id = 0;
      X = nodes_begin->X();
      Y = nodes_begin->Y();
      Z = nodes_begin->Z();

      KnotPoint = NodePointerType( new NodeType(id,X,Y,Z) );
      NewKnotsList.push_back( KnotPoint );   //reserve space for the initial auxiliar node

      id+=1;
      X = (nodes_begin + 1)->X();
      Y = (nodes_begin + 1)->Y();
      Z = (nodes_begin + 1)->Z();

      KnotPoint = NodePointerType( new NodeType(id,X,Y,Z) );
      NewKnotsList.push_back( KnotPoint );  //first node

      id+=1; //start with 2;
      for(int i=1; i < m; i++)
	{

	  //Find the spline segment indexed by j [tj, tj+1] which satisfies the current segment length (i*Length)
	  S  = 0;
	  int j = 0;
	  for(int k = 0; k < n_segments; k++)
	    {
	      S += SegmentArchLengths[k];
	      //std::cout<<" k "<<k<<std::endl;
	      if( S >= i*Length ){
		j  = k+1;
		S -= SegmentArchLengths[k];
		break;
	      }

	    }

	  if( (i*Length - S) < 0 )
	    std::cout<<" Something is wrong in the index search KNOT["<<j<<"]: length "<<i*Length<<" S "<<S<<std::endl;


	  //std::cout<<" KNOT ["<<j<<"]"<<std::endl;


	  //Compute new knot using the Bisection method

	  int max_iters    = 500;
	  double tolerance = (i*Length - S) * 1e-3;
	  double error     = tolerance * 10;

	  SplineType FirstHalfSpline;
	  SplineType SecondHalfSpline;

	  // Start with the segment found previously
	  SetSpline(Spline, rKnotsList, j);

	  double DeltaS = 0;
	  int iters = 0;

	  double rj     = 0;
	  double left   = 0;
	  double right  = 1;
	  double middle = 0.5;

	  //std::cout<<" Tolerance "<<tolerance<<std::endl;

	  //std::cout<<" LENGTH "<<Length<<" S "<<S<<" i "<<i<<" diff "<<(i*Length - S)<<std::endl;

	  while ( fabs(error) > tolerance && iters < max_iters )
	    {
	      middle = 0.5 * (left + right);

	      //Get Segment Arch-Length by the numerical integration method
	      DeltaS  = ArchLengthGeometricIntegration(Spline, rj, middle);

	      //std::cout<<" ["<<left<<", "<<middle<<"]: "<<(i*Length - S)<<"("<<DeltaS<<") :: ";

	      if( (DeltaS) < (i*Length - S) ){ //solution in the second half
		left  = middle;
	      }
	      else{//solution in the first half
	       	right = middle;
     	      }

	      error  = (DeltaS) - (i*Length - S);

	      //std::cout<<DeltaS<<" :: "<<error<<std::endl;
	      //std::cout<<DeltaS<<std::endl;

	      if( left == middle )

	      iters++;
	    }

	  if( fabs(error) > tolerance || iters == max_iters )
	    std::cout<<" max iters reached in Bisection "<<iters<<" error: "<<error<<" ["<<tolerance<<"]"<<std::endl;

	  PointType Point;
	  Point = PointOnCurve(Point,Spline,middle);

	  //std::cout<<" Knot Point "<<Point<<" middle "<<middle<<std::endl;

	  X = Point[0];
	  Y = Point[1];
	  Z = Point[2];

	  KnotPoint = NodePointerType( new NodeType(id,X,Y,Z) );
	  NewKnotsList.push_back( KnotPoint );


	  //std::cout<<" X "<<X<<" Y "<<Y<<" Z "<<Z<<std::endl;

	  id++;

	}

      //last knot end
      nodes_end =  rKnotsList.end()-2;

      X = (nodes_end)->X();
      Y = (nodes_end)->Y();
      Z = (nodes_end)->Z();

      // std::cout<<"Last knot end:  "<<id<<" X "<<X<<" Y "<<Y<<" Z "<<Z<<std::endl;

      KnotPoint = NodePointerType( new NodeType(id,X,Y,Z) );
      NewKnotsList.push_back( KnotPoint );

      //reflection at the end ( 1 extra knots )
      id += 1;
      nodes_end =  NewKnotsList.end()-1;

      X = 2.0 * (nodes_end)->X() - (nodes_end - 1)->X();
      Y = 2.0 * (nodes_end)->Y() - (nodes_end - 1)->Y();
      Z = 2.0 * (nodes_end)->Z() - (nodes_end - 1)->Z();

      // std::cout<<"Last knot reflection 1: "<<id<<" X "<<X<<" Y "<<Y<<" Z "<<Z<<std::endl;

      KnotPoint = NodePointerType( new NodeType(id,X,Y,Z) );
      NewKnotsList.push_back( KnotPoint );

      //reflection at the end ( 2 extra knots ) :: needed for the tube surface mesh
      id += 1;
      nodes_end =  NewKnotsList.end()-1;

      X = 2.0 * (nodes_end)->X() - (nodes_end - 1)->X();
      Y = 2.0 * (nodes_end)->Y() - (nodes_end - 1)->Y();
      Z = 2.0 * (nodes_end)->Z() - (nodes_end - 1)->Z();

      // std::cout<<"Last knot reflection 2: "<<id<<" X "<<X<<" Y "<<Y<<" Z "<<Z<<std::endl;

      KnotPoint = NodePointerType( new NodeType(id,X,Y,Z) );
      NewKnotsList.push_back( KnotPoint );

      //first knot //reflection at the begining ( 1 extra knot )
      nodes_begin = NewKnotsList.begin() + 1;

      X = 2.0 * nodes_begin->X() - (nodes_begin + 1)->X();
      Y = 2.0 * nodes_begin->Y() - (nodes_begin + 1)->Y();
      Z = 2.0 * nodes_begin->Z() - (nodes_begin + 1)->Z();

      // std::cout<<"First  Node Reflection: X "<<X<<" Y "<<Y<<" Z "<<Z<<std::endl;
      // std::cout<<"Second Node Reflection: X "<<nodes_begin->X()<<" Y "<<nodes_begin->Y()<<" Z "<<nodes_begin->Z()<<std::endl;
      // std::cout<<"Third  Node Reflection: X "<<(nodes_begin + 1)->X()<<" Y "<<(nodes_begin + 1)->Y()<<" Z "<<(nodes_begin + 1)->Z()<<std::endl;

      id = 0;
      KnotPoint = NodePointerType( new NodeType(id,X,Y,Z) );
      NewKnotsList.front() =  KnotPoint;
      //NewKnotsList[0] = KnotPoint;

      // nodes_begin = NewKnotsList.begin() + 3;
      // std::cout<<"Forth Node Reflection: X "<<nodes_begin->X()<<" Y "<<nodes_begin->Y()<<" Z "<<nodes_begin->Z()<<std::endl;
      // std::cout<<"Fifth  Node Reflection: X "<<(nodes_begin + 1)->X()<<" Y "<<(nodes_begin + 1)->Y()<<" Z "<<(nodes_begin + 1)->Z()<<std::endl;

      //Define the knots of Q(s) using the m+1 equally spaced points  (Spline with m cubic segments with the same arch length)
      rKnotsList.swap(NewKnotsList);
      // rKnotsList.resize(0);
      // rKnotsList.clear();
      // for(int i=0; i<NewKnotsList.size(); i++)
      //   rKnotsList.push_back(NewKnotsList[i]);

      //Set correct Ids
      for(unsigned int i=0; i<rKnotsList.size(); i++){
	rKnotsList[i]->SetId(i);
	//std::cout<<" rKnotsList[i] "<<i<<": "<<rKnotsList[i]->Coordinates()<<std::endl;
      }

    }

    //************************************************************************************
    //************************************************************************************
    double AdaptiveIntegration(SplineType& rSpline)
    {

      //Compute the numerical integration of the arch length for the Spline segment or interval (rSpline)
      double S = 0;
      S = IntegrateSubInterval(rSpline);

      //std::cout<<" Initial Arch Length "<<S<<std::endl;

      //Start adaptative sub interval integration
      int max_iters    = 40;
      double tolerance = S * 1e-2;
      double error     = tolerance * 10;

      double Sk = 0;
      int iters = 1;

      std::vector<SplineType> Intervals;
      Intervals.push_back(rSpline);

      while( fabs(error) > tolerance && iters < max_iters )
	{
	  Sk = 0;

	  //Compute SubIntervals arch length Sk and store SubIntervals Splines
	  Sk = IntegrateSubInterval(rSpline, iters+1);

	  //Compute error
	  error = Sk - S;

	  //std::cout<<" SubInterval Arch Length "<<Sk<<" error "<<error<<" tol "<<tolerance<<std::endl;

	  //Update Total arch length S = Sk
	  S = Sk;

	  iters++;

	}

      if( fabs(error) > tolerance || iters == max_iters )
	std::cout<<" max iters reached in Adaptive Integration "<<iters<<" error: "<<error<<" ["<<tolerance<<"]"<<std::endl;


      return S;
    }



    //************************************************************************************
    //************************************************************************************

    //Divide the interval in n  sub-intervals (two halves) and compute the arch-lengths sum

    double IntegrateSubInterval(SplineType& rSpline, int n = 0)
    {

      double t = 1.0 / double(n + 1.0) ;

      double S  = 0;
      double a = 0;
      double b = t;

      for(int i=0; i<n+1; i++)
	{
	  S += ArchLengthGeometricIntegration(rSpline, a, b);
	  //std::cout<<" S "<<S<<" a "<<a<<" b "<<b<<std::endl;
	  a += t;
	  b += t;
	}

      return S;

    }


    //************************************************************************************
    //************************************************************************************
    double ArchLengthGeometricIntegration(SplineType& rSpline, double& a, double& b)
    {
      //apply numerical integration:: gaussian quadrature
      //return GaussianQuadratureIntegration(rSpline); //not implemented yet

      //apply numerical integration:: simpson rule
      return SimpsonRuleIntegration(rSpline, a, b, 7);

    }

    //************************************************************************************
    //************************************************************************************

    //(Simpson's 3/8 rule :: Cubic Interpolation)
    //Approximates the definite integral of f(t) (rSpline) from a to b by
    //the composite Simpson's rule, using n subintervals (n even number)
    double SimpsonRuleIntegration(SplineType& rSpline, double& a, double& b, double n = 1)
    {

      double IntegralValue = 0;

      double h = (b - a) / n;

      double s = SplineGeometricLength(rSpline, a) + SplineGeometricLength(rSpline, b);

      double t = 0;
      for(int i=1; i<n; i+=2)
	{
	  t = (a + i * h);
	  s += 4 * SplineGeometricLength(rSpline, t);
	}


      for(int i=2; i<n-1; i+=2)
	{
	  t = (a + i * h);
	  s += 2 * SplineGeometricLength(rSpline, t);
	}


      IntegralValue = (s * h / 3.0);

      return IntegralValue;

    }

    //************************************************************************************
    //************************************************************************************
    double SplineGeometricLength(SplineType& rSpline, double& t)
    {

      PointType PointA = ZeroVector(3);
      PointA = PointOnCurveFirstDerivative(rSpline, t);

      return (norm_2(PointA));
    }



    //************************************************************************************
    //************************************************************************************
    inline void SetSpline(SplineType& rOutputSpline,const SplineType& rInputSpline)
    {

      rOutputSpline.id = rInputSpline.id;

      rOutputSpline.P0 = rInputSpline.P0;
      rOutputSpline.P1 = rInputSpline.P1;
      rOutputSpline.P2 = rInputSpline.P2;
      rOutputSpline.P3 = rInputSpline.P3;

    }

    //************************************************************************************
    //************************************************************************************
    inline void SetSpline(SplineType& rSpline,const PointType&  P0,const PointType&  P1,const PointType&  P2,const PointType&  P3)
    {

      rSpline.id = 0;

      rSpline.P0 = P0;
      rSpline.P1 = P1;
      rSpline.P2 = P2;
      rSpline.P3 = P3;

    }

    //************************************************************************************
    //************************************************************************************
    inline void SetSpline(SplineType& rSpline, const NodePointerTypeVector& rKnotsList, int& id)
    {

      rSpline.id = id;

      rSpline.P0 = ZeroVector(3);
      rSpline.P0[0] = rKnotsList[id-1]->X();
      rSpline.P0[1] = rKnotsList[id-1]->Y();
      rSpline.P0[2] = rKnotsList[id-1]->Z();

      rSpline.P1 = ZeroVector(3);
      rSpline.P1[0] = rKnotsList[id]->X();
      rSpline.P1[1] = rKnotsList[id]->Y();
      rSpline.P1[2] = rKnotsList[id]->Z();

      rSpline.P2 = ZeroVector(3);
      rSpline.P2[0] = rKnotsList[id+1]->X();
      rSpline.P2[1] = rKnotsList[id+1]->Y();
      rSpline.P2[2] = rKnotsList[id+1]->Z();

      rSpline.P3 = ZeroVector(3);
      rSpline.P3[0] = rKnotsList[id+2]->X();
      rSpline.P3[1] = rKnotsList[id+2]->Y();
      rSpline.P3[2] = rKnotsList[id+2]->Z();

    }

   //************************************************************************************
    //************************************************************************************

    PointType& CalculatePointProjection(const PointType& rPoint, KdtreeType& rKnotsKdtree, const NodePointerTypeVector& rKnotsList, PointType& rPointProjection )
    {
      KRATOS_TRY

	//1.- Find the closest generatrix knot of the spline curve
	double PointDistance = 0;

        int id = rKnotsList.front()->Id(); // starting with the first Knot

	id = GetClosestKnotId( rPoint, rKnotsKdtree, PointDistance ); //starting with the closest Knot

        SplineType Spline;

	SetSpline(Spline, rKnotsList, id);

	//2.- Find the closest point on a spline curve:  (Sk := NormalizedArchLength)
	double Sk = CombinedMethod(rPoint, rKnotsList, Spline);

	//2.1-Projected Point:
	rPointProjection = PointOnCurve(rPointProjection,Spline,Sk);

	//std::cout<<"CD: KnotPoint "<<Spline.P1<<" ProjectedPoint "<<rPointProjection<<" BeamPoint "<<rPoint<<std::endl;

	return rPointProjection;

      KRATOS_CATCH( "" )
    }


    //************************************************************************************
    //************************************************************************************
    int GetClosestKnotId(const PointType& rPoint, KdtreeType& rKnotsKdtree, double& rPointDistance)
    {

      NodeType WorkPoint(0,rPoint[0],rPoint[1],rPoint[2]);
      rPointDistance = 0;

      NodePointerType NearestPoint = rKnotsKdtree.SearchNearestPoint(WorkPoint, rPointDistance);

      // get the id of the closest point i
      return NearestPoint->Id();

    }

    //************************************************************************************
    //************************************************************************************

    double CombinedMethod(const PointType& rPoint, const NodePointerTypeVector& rKnotsList, SplineType& rSpline, double s = 0.5)
    {
      //1.2.- Combined Method:
      int iters = 4;

      //1.2.1- Apply Quadratic Minimization Method
      double Sk = QuadraticMinimizationMethod(rPoint, rKnotsList, rSpline, iters);

      //std::cout<<" First prediction Sk "<<Sk<<std::endl;

      iters = 20;

      //1.2.2- Apply Newton's Method
      Sk = NewtonsMethod(rPoint, rKnotsList, rSpline, Sk, iters);



      return Sk;

    }


    //************************************************************************************
    //************************************************************************************

    double NewtonsMethod(const PointType& rPoint, const NodePointerTypeVector& rKnotsList, SplineType& rSpline, double Spredict = 0, double iters = 20, double s = 0.5)
    {

      //Initial normalized estimate on the selected spline
      double Sk = Spredict;
      int max_iters  = iters;
      int dist_iters = 25;

      int iter = 0;
      int segment_iter = 0;
      while( ((Sk < -0.1 || Sk >1) && segment_iter<max_iters) || segment_iter == 0 )
	{

	  if( Sk<0 ){ //previous adjacent segment

	    int new_id = rSpline.id - 1;

	    if(new_id <= 0)
	      new_id = 1;

	    //std::cout<<" NewId R "<<new_id<<std::endl;

	    Sk = 1 - Sk;

	    SetSpline(rSpline, rKnotsList, new_id);

	  }

	  if( Sk>1 ){ //posterior adjacent segment

	    int new_id = rSpline.id + 1;

	    if( new_id > (int)rKnotsList.size()-3 )
	      new_id = rKnotsList.size()-3;

	    //std::cout<<" NewId I "<<new_id<<std::endl;

	    Sk = Sk - 1;

	    SetSpline(rSpline, rKnotsList, new_id);

	  }

	  //std::cout<<" SPLINE ID "<<rSpline.id<<std::endl;

	  //Iterative Parameters
	  double tolerance = 1e-7; //1e-8

	  iter = 0;
	  double distance = tolerance*10;

	  while( iter<dist_iters && distance>=tolerance )
	    {

	      distance  = FirstDerivativeSquareDistancePointToSpline(rPoint,rSpline, Sk);
	      distance /= SecondDerivativeSquareDistancePointToSpline(rPoint,rSpline, Sk);

	      Sk -= distance;

	      //std::cout<<" Sk newton  "<<Sk<<" distance "<<distance<<std::endl;

	      distance = fabs(distance);

	      iter++;
	    }


	  // std::cout<<" iters "<<iter<<" distance "<<distance<<std::endl;

	  if( distance > tolerance ){
	    std::cout<<"BBX Tube Contact Search [Point:"<<rPoint<<",Distance:"<<distance<<",Tol:"<<tolerance<<"] ERROR"<<std::endl;
	    if( distance > 1e3*tolerance )
	      KRATOS_THROW_ERROR( std::logic_error," TUBE BBX::NewmarksMethod HAS NOT REACHED CONVERGENCE ", 0)
	    //std::cout<<" TUBE BBX::NewmarksMethod HAS NOT REACHED CONVERGENCE "<<std::endl;
	  }
	  else{
	    //std::cout<<" converged "<<Sk<<std::endl;
	  }

	  segment_iter++;
	}

      //std::cout<<" segment_iter "<<segment_iter<<" newtons iters "<<iter<<std::endl;

      if( Sk < 0 && Sk > -0.1)
	Sk = 0;

      return Sk;

    }

    //************************************************************************************
    //************************************************************************************

    double QuadraticMinimizationMethod(const PointType& rPoint, const NodePointerTypeVector& rKnotsList, SplineType& rSpline, double iters, double s = 0.5)
    {

      int max_iters = iters;
      double Skj = 0; //to start with the iteration

      int segment_iter = 0;
      while( ((Skj < 0 || Skj >1) && segment_iter<max_iters) || segment_iter == 0 )
	{

	  if( Skj<0 ){ //previous adjacent segment

	    int new_id = rSpline.id - 1;

	    if(new_id <= 0)
	      new_id = 1;

	    SetSpline(rSpline, rKnotsList, new_id);

	  }

	  if( Skj>1 ){ //posterior adjacent segment

	    int new_id = rSpline.id + 1;

	    if( new_id > (int)rKnotsList.size()-3 )
	      new_id = rKnotsList.size()-3;

	    SetSpline(rSpline, rKnotsList, new_id);

	  }

	  //std::cout<<" Q SPLINE ID "<<rSpline.id<<std::endl;

	  //Initial Segment estimate :: rSpline

	  //Initial three estimates for a normalized segment
	  Vector Estimates = ZeroVector(3);
	  Estimates[0] = 0.0; //S1
	  Estimates[1] = 0.5; //S2
	  Estimates[2] = 1.0; //S3

	  //1:
	  //Arch Length difference between the estimates splines
	  Vector Difference       = ZeroVector(3);

	  //Arch Length square difference between the estimates splines
	  Vector SquareDifference = ZeroVector(3);

	  //2:
	  //Square Distances to the estimates splines
	  Vector Function = ZeroVector(3);

	  //Polynomial for distance interpolation
	  Vector InterpolatedDistance = ZeroVector(4);


	  //Iterative Parameters
	  double tolerance = 1e-7; //1e-8;
	  double distance  = tolerance * 10;

	  int iter    = 0;
	  double Ski  = 0;

	  while( iter < max_iters && distance >= tolerance )
	    {

	      //Store previous Sk
	      Ski = Skj;

	      //Arch Length difference between the estimates splines
	      Difference       = CalculateArchLengthDifferences(Difference,Estimates);

	      //Arch Length square difference between the estimates splines
	      SquareDifference = CalculateSquareArchLengthDifferences(SquareDifference,Estimates);

	      //Square Distances to the estimates splines
	      Function[0] = SquareDistancePointToSpline(rPoint, rSpline, Estimates[0]);
	      Function[1] = SquareDistancePointToSpline(rPoint, rSpline, Estimates[1]);
	      Function[2] = SquareDistancePointToSpline(rPoint, rSpline, Estimates[2]);

	      Skj  = 0.5 * (SquareDifference[1]*Function[0] + SquareDifference[2]*Function[1] + SquareDifference[0]*Function[2]);
	      Skj /= (Difference[1]*Function[0] + Difference[2]*Function[1] + Difference[0]*Function[2]);

	      InterpolatedDistance[0] = EvaluateDistancePolynomial(Function, Estimates, Estimates[0]);
	      InterpolatedDistance[1] = EvaluateDistancePolynomial(Function, Estimates, Estimates[1]);
	      InterpolatedDistance[2] = EvaluateDistancePolynomial(Function, Estimates, Estimates[2]);

	      distance = EvaluateDistancePolynomial(Function, Estimates, Skj);

	      double larger_distance = fabs(InterpolatedDistance[0]);

	      int largest = 0;
	      for( unsigned int i=1; i<3; i++)
		{
		  if(larger_distance<fabs(InterpolatedDistance[i])){
		    larger_distance = fabs(InterpolatedDistance[i]);
		    largest = i;
		  }
		}

	      if(fabs(InterpolatedDistance[largest])>distance)
		{
		  Estimates[largest] = Skj;
		}


	      //tolerance
	      distance = fabs(Skj-Ski);

	      iter++;
	    }

	  //std::cout<<" iters "<<iter<<" Skj "<<Skj<<std::endl;

	  // if( distance > tolerance )
	  // 	std::cout<<" TUBE BBX::QuadraticMinimizationMethod HAS NOT REACHED CONVERGENCE "<<std::endl;

	  segment_iter++;
	}


      return Skj;

      //Projected Point:
      //PointType ProjectedPoint;
      //ProjectedPoint = PointOnCurve(ProjectedPoint,rSpline,Sk)


    }

    //************************************************************************************
    //************************************************************************************

    double EvaluateDistancePolynomial(const Vector& rFunction,const Vector& rEstimates, double& Sk)
    {
      //Polynomial to interpolate D(S) coefficients:
      Vector PolynomialBasis =  ZeroVector(3);

      //S1, S2, S3
      PolynomialBasis  = CalculateDistancePolynomialBasis(PolynomialBasis, rEstimates, Sk);
      double PolynomialValue = PolynomialBasis[0]*rFunction[0] + PolynomialBasis[1]*rFunction[1] + PolynomialBasis[2]*rFunction[2];

      return PolynomialValue;
    }

    //************************************************************************************
    //************************************************************************************

    Vector& CalculateArchLengthDifferences(Vector& rDifference, const Vector& rEstimates)
    {
      rDifference[0] = rEstimates[0]-rEstimates[1]; //12 (1-2)
      rDifference[1] = rEstimates[1]-rEstimates[2]; //23 (2-3)
      rDifference[2] = rEstimates[2]-rEstimates[0]; //31 (3-1)

      return rDifference;
    }

    //************************************************************************************
    //************************************************************************************

    Vector& CalculateSquareArchLengthDifferences(Vector& rSquareDifference, const Vector& rEstimates)
    {
      rSquareDifference[0] = rEstimates[0] * rEstimates[0] - rEstimates[1] * rEstimates[1]; //12 (1-2)
      rSquareDifference[1] = rEstimates[1] * rEstimates[1] - rEstimates[2] * rEstimates[2]; //23 (2-3)
      rSquareDifference[2] = rEstimates[2] * rEstimates[2] - rEstimates[0] * rEstimates[0]; //31 (3-1)

      return rSquareDifference;
    }

    //************************************************************************************
    //************************************************************************************

    // Values of the coefficients for the polynomial that interpolates the square distance function

    Vector& CalculateDistancePolynomialBasis(Vector& rPolynomialBasis, const Vector& rEstimates, double& rValue)
    {
      rPolynomialBasis[0] = (rValue-rEstimates[1]) * (rValue-rEstimates[2]) / ((rEstimates[0]-rEstimates[1]) * (rEstimates[0]-rEstimates[2]));
      rPolynomialBasis[1] = (rValue-rEstimates[0]) * (rValue-rEstimates[2]) / ((rEstimates[1]-rEstimates[0]) * (rEstimates[1]-rEstimates[2]));
      rPolynomialBasis[2] = (rValue-rEstimates[0]) * (rValue-rEstimates[1]) / ((rEstimates[2]-rEstimates[0]) * (rEstimates[2]-rEstimates[1]));

      return rPolynomialBasis;
    }


    //************************************************************************************
    //************************************************************************************


    double SquareDistancePointToSpline(const PointType& rPoint,const SplineType& rSpline, double& t)
    {
      //compute spline on t, normalized distance
      PointType SplinePoint = ZeroVector(3);
      SplinePoint = PointOnCurve(SplinePoint, rSpline, t);

      //compute square distance
      SplinePoint -= rPoint;
      double Distance = inner_prod(SplinePoint,SplinePoint);

      return Distance;
    }

    //************************************************************************************
    //************************************************************************************

    double FirstDerivativeSquareDistancePointToSpline(const PointType& rPoint,const SplineType& rSpline, double& t)
    {
      //compute spline on t, normalized distance
      PointType SplinePoint = ZeroVector(3);
      SplinePoint = PointOnCurve(SplinePoint, rSpline, t);

      SplinePoint -= rPoint;

      //compute spline first derivative on t, normalized distance
      PointType SplinePointFirstDerivative = PointOnCurveFirstDerivative(rSpline, t);

      double Distance = 2.0 * inner_prod(SplinePoint,SplinePointFirstDerivative);

      return Distance;
    }


    //************************************************************************************
    //************************************************************************************

    double SecondDerivativeSquareDistancePointToSpline(const PointType& rPoint,const SplineType& rSpline, double& t)
    {
      //compute spline on t, normalized distance
      PointType SplinePoint = ZeroVector(3);
      SplinePoint = PointOnCurve(SplinePoint, rSpline, t);

      SplinePoint -= rPoint;

      //compute spline first derivative on t, normalized distance
      PointType SplinePointFirstDerivative = PointOnCurveFirstDerivative(rSpline, t);

      //compute spline second derivative on t, normalized distance
      PointType SplinePointSecondDerivative = PointOnCurveSecondDerivative(rSpline, t);

      double Distance = inner_prod(SplinePointFirstDerivative,SplinePointFirstDerivative);
      Distance += inner_prod(SplinePoint,SplinePointSecondDerivative);
      Distance *= 2.0;

      return Distance;
    }



    //************************************************************************************
    //************************************************************************************

    /// Return a point on the curve between P1 and P2 with P0 and P3 describing curvature, at
    /// the normalized distance t, and the spline parameter s

    inline PointType& PointOnCurve(PointType& rPoint, const SplineType& rSpline, double& t, double s = 0.5)
    {
      Vector Basis = ZeroVector(4);

      Basis = SplineBasis( Basis, t, s );

      rPoint = Basis[0] * rSpline.P0 + Basis[1] * rSpline.P1 + Basis[2] * rSpline.P2 + Basis[3] * rSpline.P3;

      return rPoint;
    }



    //************************************************************************************
    //************************************************************************************

    /// Return a point on the curve between P1 and P2 with P0 and P3 describing curvature, at
    /// the normalized distance t, and the spline parameter s

    PointType PointOnCurveFirstDerivative(const SplineType& rSpline, double& t, double s = 0.5)
    {
      Vector Basis = ZeroVector(4);

      Basis = FirstDerivativeSplineBasis( Basis, t, s );

      PointType Result = Basis[0] * rSpline.P0 + Basis[1] * rSpline.P1 + Basis[2] * rSpline.P2 + Basis[3] * rSpline.P3;

      return Result;
    }


    //************************************************************************************
    //************************************************************************************

    /// Return a point on the curve between P1 and P2 with P0 and P3 describing curvature, at
    /// the normalized distance t, and the spline parameter s

    PointType PointOnCurveSecondDerivative(const SplineType& rSpline, double& t, double s = 0.5)
    {
      Vector Basis = ZeroVector(4);

      Basis = SecondDerivativeSplineBasis( Basis, t, s );

      PointType Result = Basis[0] * rSpline.P0 + Basis[1] * rSpline.P1 + Basis[2] * rSpline.P2 + Basis[3] * rSpline.P3;

      return Result;
    }


    //************************************************************************************
    //************************************************************************************

    /// Return the cubic basis for a Catmull-Rom spline
    /// set spline general coefficients for a given segment

    static inline std::vector<Vector>& SplineCoefficients(const SplineType& rSpline, std::vector<Vector>& rCoefficients, double s = 0.5)
    {
      if( rCoefficients.size() != 4 )
	rCoefficients.resize(4);

      rCoefficients[0] = (rSpline.P2 - rSpline.P0) * s + rSpline.P1 * (2.0 - s) + rSpline.P1 * (s - 2.0);
      rCoefficients[1] = (2.0 * rSpline.P0 - rSpline.P2) * s + rSpline.P1 * (s - 3.0) + rSpline.P1 * (3.0 - 2.0 * s);
      rCoefficients[2] = (rSpline.P1 - rSpline.P0) * s;
      rCoefficients[3] = rSpline.P0;

      return rCoefficients;
    }


    //************************************************************************************
    //************************************************************************************

    /// Return the cubic basis for a Catmull-Rom spline
    /// set a normalized distance t, and the spline parameter s

    static inline Vector& SplineBasis(Vector& Basis, double& t, double s = 0.5)
    {
      if( Basis.size() != 4 )
	Basis.resize(4,false);

      Basis[0] = ((-t + 2.0) * t - 1.0) * t * s;
      Basis[1] = ((((2.0/s - 1.0) * t + (1.0 - 3.0/s)) * t) * t + 1.0/s) * s;
      Basis[2] = (((1.0 - 2.0/s) * t + (3.0/s - 2.0)) * t + 1.0) * t * s;
      Basis[3] = ((t - 1.0) * t * t) * s;

      return Basis;
    }


    //************************************************************************************
    //************************************************************************************

    /// Return the cubic basis for a Catmull-Rom spline first derivative
    /// set a normalized distance t, and the spline parameter s

    static inline Vector& FirstDerivativeSplineBasis(Vector& Basis, double& t, double s = 0.5)
    {
      if( Basis.size() != 4 )
	Basis.resize(4,false);

      Basis[0] = ((-3.0 * t + 4.0) * t - 1.0) * s;
      Basis[1] = (((6.0/s - 3.0) * t + (2.0 - 6.0/s)) * t) * s;
      Basis[2] = (((3.0 - 6.0/s) * t + (6.0/s - 4.0)) * t + 1.0) * s;
      Basis[3] = (3.0 * t - 2.0) * t * s;

      return Basis;
    }

    //************************************************************************************
    //************************************************************************************

    /// Return the cubic basis for a Catmull-Rom spline second derivative
    /// set a normalized distance t, and the spline parameter s

    static inline Vector& SecondDerivativeSplineBasis(Vector& Basis, double& t, double s = 0.5)
    {
      if( Basis.size() != 4 )
	Basis.resize(4,false);

      Basis[0] = (-6.0 * t + 4.0) * s;
      Basis[1] = ((12.0/s - 6.0) * t + (2.0 - 6.0/s)) * s;
      Basis[2] = ((6.0 - 12.0/s) * t + (6.0/s - 4.0)) * s;
      Basis[3] = (6.0 * t - 2.0) * s;

      return Basis;
    }


    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{


    ///@}


  protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}
  private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    int mEchoLevel;

    bool mParallel;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Unaccessible methods
    ///@{

    ///@}


  }; // Class SplineCurveUtilities

  ///@}
  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  ///@}

} // namespace Kratos.

#endif // KRATOS_SPLINE_CURVE_UTILITIES_H_INCLUDED defined


