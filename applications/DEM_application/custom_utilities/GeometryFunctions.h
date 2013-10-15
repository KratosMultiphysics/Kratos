/* 
 * File:   GeometryFunctions.h
 * Author: msantasusana
 *
 * Created on 21 de mayo de 2012, 19:40
 */

#ifndef _GEOMETRYFUNCTIONS_H
#define	_GEOMETRYFUNCTIONS_H


#include <cmath>
/////////////////////****************************/////////////////////////////////
/////////////////////****************************/////////////////////////////////
/////////////////////May, 2012 , Created by Chun Feng/////////////////////////////////
/////////////////////****************************/////////////////////////////////
/////////////////////****************************/////////////////////////////////
namespace Kratos
{
    namespace GeometryFunctions
    {

    static inline int sign(double a)
    {
        double output;
        if(a < 0.0) output = -1.0;
        else if (a > 0.0) output = 1.0;
        else output = 0.0;
        return output;
    }
      
    static inline double min(double a, double b)
    {
        double output;
        if (a<=b) output = a;
        else output = b;
        return output;
    }

    static inline double max(double a, double b)
    {
        double output;
        if (a>=b) output = a;
        else output = b;
        return output;
    }

    static inline void norm(double Vector[3])
    {
            double distance = sqrt(Vector[0] * Vector[0] + Vector[1] * Vector[1] + Vector[2] * Vector[2]);
            Vector[0] = Vector[0] / distance;
            Vector[1] = Vector[1] / distance;
            Vector[2] = Vector[2] / distance;
    }
    
    static inline void norm( array_1d<double,3>& Vector, double& distance)
    {
            distance = sqrt(Vector[0] * Vector[0] + Vector[1] * Vector[1] + Vector[2] * Vector[2]);
            
            if(distance != 0.0)
            {
            Vector[0] = Vector[0] / distance;
            Vector[1] = Vector[1] / distance;
            Vector[2] = Vector[2] / distance;
            }

      }

    static inline void VectorGlobal2Local(double LocalCoordSystem[3][3], double GlobalVector[3], double LocalVector[3])
    {
    	int i,j;

	for (i=0; i<3; i++)
	{
            LocalVector[i] = 0.0;
            for (j=0; j<3; j++)
            {
                LocalVector[i]+=LocalCoordSystem[i][j]*GlobalVector[j];
            }
        }
    }

    static inline void VectorLocal2Global(double LocalCoordSystem[3][3], double LocalVector[3], double GlobalVector[3])
    {
    	int i,j;

	for (i=0; i<3; i++)
	{
            GlobalVector[i] = 0.0;
            for (j=0; j<3; j++)
            {
		GlobalVector[i]+=LocalCoordSystem[j][i]*LocalVector[j];
            }
	}
    }

    static inline double DotProduct(double Vector1[3], double Vector2[3])
    {
        return Vector1[0] * Vector2[0] + Vector1[1] * Vector2[1] + Vector1[2] * Vector2[2];
    }

    static inline double DotProduct(const array_1d<double,3> Vector1, const array_1d<double,3> Vector2)
    {
        return Vector1[0] * Vector2[0] + Vector1[1] * Vector2[1] + Vector1[2] * Vector2[2];
    }

    static inline void CrossProduct(double u[3], double v[3], double ReturnVector[3])
    {
    	ReturnVector[0] = u[1]*v[2] - u[2]*v[1];
	ReturnVector[1] = v[0]*u[2] - u[0]*v[2];
	ReturnVector[2] = u[0]*v[1] - u[1]*v[0];
    }

    static inline void CrossProduct( const array_1d<double,3>& u, const array_1d<double,3>& v, array_1d<double,3>& ReturnVector)
    {
    	ReturnVector[0] = u[1]*v[2] - u[2]*v[1];
	ReturnVector[1] = v[0]*u[2] - u[0]*v[2];
	ReturnVector[2] = u[0]*v[1] - u[1]*v[0];
    }

   //NOTE:: Modified by M. Santasusana Feb 2013 - simplification (may be the one proposed by F.Chun was for a more generalized case) 
    static inline void ComputeContactLocalCoordSystem(double NormalDirection[3], double LocalCoordSystem[3][3])  //inline: modifies the LocalCoordSystem as it were a reference
    {
        int ii;
        double Vector0[3] = {0.0},Vector1[3] = {0.0};

        norm(NormalDirection);

        if(fabs(NormalDirection[0])>=0.577)
        {
            Vector0[0]= - NormalDirection[1]/NormalDirection[0];
            Vector0[1]=1.0;
            Vector0[2]=0.0;
        }
        else if(fabs(NormalDirection[1])>=0.577)
        {

            Vector0[0]=0.0;
            Vector0[1]= - NormalDirection[2]/NormalDirection[1];
            Vector0[2]=1.0;
        }
        
        else
        {                   
            Vector0[0]= 1.0;
            Vector0[1]= 0.0;
            Vector0[2]= - NormalDirection[0] /NormalDirection[2];
        }

        norm(Vector0);
        CrossProduct(NormalDirection,Vector0,Vector1);
        norm(Vector1);
        
        for(ii=0;ii<3;ii++)
        {
            LocalCoordSystem[0][ii]=Vector0[ii];
            LocalCoordSystem[1][ii]=Vector1[ii];
            LocalCoordSystem[2][ii]=NormalDirection[ii];
        }

    }

     static inline double DistanceOfTwoPoint(double coord1[3], double coord2[3])
     {
         double dx = coord1[0] - coord2[0];
         double dy = coord1[1] - coord2[1];
         double dz = coord1[2] - coord2[2];

         return sqrt(dx * dx + dy * dy + dz * dz);
     }



     static inline void  TriAngleArea(double Coord1[3],double Coord2[3],double Coord3[3],double &area)
    {
	int k;
	double Vector1[3],Vector2[3],Vector0[3];
	for(k = 0;k < 3; k++)
	{
            Vector1[k] = Coord3[k] - Coord1[k];
	    Vector2[k] = Coord2[k] - Coord1[k];
	}

	CrossProduct(Vector1, Vector2, Vector0);
	area= sqrt(Vector0[0] * Vector0[0] + Vector0[1] * Vector0[1] + Vector0[2] * Vector0[2]) / 2.0;
     }

     //TriAngle Wight, coord1,coord2,coord3,testcoord,weight
     static inline void TriAngleWeight(double Coord1[3], double Coord2[3], double Coord3[3], double JudgeCoord[3], double Weight[3])
     {
        double area[3], s;
        TriAngleArea(Coord1, Coord2, JudgeCoord, area[0]);
        TriAngleArea(Coord2, Coord3, JudgeCoord, area[1]);
        TriAngleArea(Coord3, Coord1, JudgeCoord, area[2]);

        TriAngleArea(Coord1, Coord2, Coord3, s);
        /////s = area[0] + area[1] + area[2];

        Weight[0] = area[1] / s;
        Weight[1] = area[2] / s;
        Weight[2] = area[0] / s;
     }

     static inline  double DistancePointToPlane(double CoordInPlane[3], double PlaneUnitNormalVector[3], double TestCoord[3])
     {
         double Vector1[3] = {0.0};

         Vector1[0] = TestCoord[0] - CoordInPlane[0];
         Vector1[1] = TestCoord[1] - CoordInPlane[1];
         Vector1[2] = TestCoord[2] - CoordInPlane[2];

         double dist = fabs (DotProduct(Vector1, PlaneUnitNormalVector));

         return dist;
     }

     static inline void CoordProjectionOnPlane(double CoordOut[3], double CoordIn[3], double LocalCoordSystem[3][3], double IntersectionCoord[3])
     {
         double out_coord_local[3] = {0.0};
         double in_coord_local[3]  = {0.0};

          VectorGlobal2Local(LocalCoordSystem, CoordOut, out_coord_local);
          VectorGlobal2Local(LocalCoordSystem, CoordIn,  in_coord_local);

          double vector1[3] = {0.0};
          vector1[0] = out_coord_local[0];
          vector1[1] = out_coord_local[1];
          vector1[2] = in_coord_local [2];

          VectorLocal2Global(LocalCoordSystem, vector1, IntersectionCoord);

     }


     static inline void Compute2DimElementEdgeLocalSystem(double EdgeCoord1[3], double EdgeCoord2[3], double Centroid[3], double LocalCoordSystem[3][3])
     {
         double Vector1[3] = {0.0};
         double Vector2[3] = {0.0};
         double Vector3[3] = {0.0};
         double Normal [3] = {0.0};

         Vector1[0] = EdgeCoord2[0] - EdgeCoord1[0];
         Vector1[1] = EdgeCoord2[1] - EdgeCoord1[1];
         Vector1[2] = EdgeCoord2[2] - EdgeCoord1[2];

         Vector2[0] = Centroid[0] - EdgeCoord1[0];
         Vector2[1] = Centroid[1] - EdgeCoord1[1];
         Vector2[2] = Centroid[2] - EdgeCoord1[2];

         norm(Vector1);
         norm(Vector2);
         CrossProduct(Vector1, Vector2, Vector3);
         norm(Vector3);
         CrossProduct(Vector1, Vector3, Normal);
         norm(Normal);

         if (DotProduct(Vector2, Normal) > 0.0)
         {
             for (int ia = 0; ia < 3; ia++)
             {
                    LocalCoordSystem[0][ia] = -Vector1[ia];
                    LocalCoordSystem[1][ia] = -Vector3[ia];
                    LocalCoordSystem[2][ia] = -Normal [ia];
             }
         }
         else
         {
             for (int ia = 0; ia < 3; ia++)
             {
                    LocalCoordSystem[0][ia] = Vector1[ia];
                    LocalCoordSystem[1][ia] = Vector3[ia];
                    LocalCoordSystem[2][ia] = Normal [ia];
             }
         }
     }

     static inline void Compute3DimElementEdgeLocalSystem(double FaceCoord1[3], double FaceCoord2[3], double FaceCoord3[3], double Centroid[3], double LocalCoordSystem[3][3])
     {
         double Vector1[3] = {0.0};
         double Vector2[3] = {0.0};
         double Vector3[3] = {0.0};
         double Normal [3] = {0.0};

         Vector1[0] = FaceCoord2[0] - FaceCoord1[0];
         Vector1[1] = FaceCoord2[1] - FaceCoord1[1];
         Vector1[2] = FaceCoord2[2] - FaceCoord1[2];

         Vector2[0] = FaceCoord3[0] - FaceCoord2[0];
         Vector2[1] = FaceCoord3[1] - FaceCoord2[1];
         Vector2[2] = FaceCoord3[2] - FaceCoord2[2];

         norm(Vector1);
         norm(Vector2);
         CrossProduct(Vector1, Vector2, Normal);
         norm(Normal);

         CrossProduct(Normal, Vector1, Vector2);
         norm(Vector2);

         Vector3[0] = FaceCoord1[0] - Centroid[0];
         Vector3[1] = FaceCoord1[1] - Centroid[1];
         Vector3[2] = FaceCoord1[2] - Centroid[2];
         norm(Vector3);



         if (DotProduct(Vector3, Normal) > 0.0)
         {
             for (int ia = 0; ia < 3; ia++)
             {
                    LocalCoordSystem[0][ia] = Vector1[ia];
                    LocalCoordSystem[1][ia] = Vector2[ia];
                    LocalCoordSystem[2][ia] = Normal [ia];
             }
         }
         else
         {
             for (int ia = 0; ia < 3; ia++)
             {
                    LocalCoordSystem[0][ia] = -Vector1[ia];
                    LocalCoordSystem[1][ia] = -Vector2[ia];
                    LocalCoordSystem[2][ia] = -Normal [ia];
             }
         }
     }

     ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     ///////////////******Rotate a point over an arbitrary line though an arbitrary point******/////////////////////////////////////
     ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

     static inline void RotatePointAboutArbitraryLine(array_1d<double,3>& TargetPoint, const array_1d<double,3>& CentrePoint, const array_1d<double,3>& LineVector, const double RotationAngle)
     {

  
    
          

        const double O = RotationAngle;

        double x = TargetPoint[0], a = CentrePoint[0], u = LineVector[0];
        double y = TargetPoint[1], b = CentrePoint[1], v = LineVector[1];
        double z = TargetPoint[2], c = CentrePoint[2], w = LineVector[2];

        double L = u*u+v*v+w*w;


        /*
        KRATOS_WATCH(x)
        KRATOS_WATCH(y)
                KRATOS_WATCH(z)
         KRATOS_WATCH(u)
         KRATOS_WATCH(v)
                KRATOS_WATCH(w)
         KRATOS_WATCH(a)
                KRATOS_WATCH(b)
                KRATOS_WATCH(c)
         */

        if (L==0)
        {

           // KRATOS_WATCH("WARNING! Nul vector for the rotation!")

        }

        else

        {

        TargetPoint[0] = ((a*(v*v+w*w)-u*(b*v+c*w-u*x-v*y-w*z))*(1-cos(O))+L*x*cos(O)+sqrt(L)*(-c*w+b*w-w*y+v*z)*sin(O))*(1/L);
        TargetPoint[1] = ((b*(u*u+w*w)-v*(a*u+c*w-u*x-v*y-w*z))*(1-cos(O))+L*y*cos(O)+sqrt(L)*(c*u-a*w+w*x-u*z)*sin(O))*(1/L);
        TargetPoint[2] = ((c*(u*u+v*v)-w*(a*u+b*v-u*x-v*y-w*z))*(1-cos(O))+L*z*cos(O)+sqrt(L)*(-b*u+a*v-v*x+u*y)*sin(O))*(1/L);
 
        }

     }
     
     /////////****Quaternions****///////////////     

     static inline void Align(double normal[3], double quart_imag[3], double quart_angle, double quart_axis[3])
     {   
         double temp[3] = {0.0};
         
         temp[0] = normal[0] + quart_imag[0];
         temp[1] = normal[1] + quart_imag[1];
         temp[2] = normal[2] + quart_imag[2];
         
         norm(temp);
         
         double normal_norm = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
         
         quart_angle = DotProduct(normal, temp) / normal_norm;
         
         CrossProduct(normal, temp, quart_axis);
         
         quart_axis[0] = quart_axis[0] / normal_norm;
         quart_axis[1] = quart_axis[2] / normal_norm;
         quart_axis[2] = quart_axis[2] / normal_norm;
     }
     
     ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     ///////////////******EULER ANGLES from 2 vectors******/////////////////////////////////////////////////////////////////////////
     ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    /*static inline void CalculateEulerAngles(const array_1d<double,3>& OriginalVector_X, const array_1d<double,3>& OriginalVector_Z,
                const array_1d<double,3>& RotatedVector_X, const array_1d<double,3>& RotatedVector_Z, array_1d<double,3>& EulerAngles)
    {

        array_1d< double,3 > N = ZeroVector(3);



        CrossProduct( OriginalVector_Z, RotatedVector_Z, N);

        double return1 = DotProduct(N,OriginalVector_X);   //cos(Alpha)
        double return2 = DotProduct(OriginalVector_Z, RotatedVector_Z); //cos(Beta)
        double return3 = DotProduct(N,RotatedVector_X); //cos(Gamma)

        EulerAngles[0] = acos(return1);
        EulerAngles[1] = acos(return2);
        EulerAngles[2] = acos(return3);

    }*/
    static inline bool JudgeIfThisEdgeIsContactWithParticle(double EdgeCoord1[3], double EdgeCoord2[3], double Centroid[3], double Particle_Coord[3], double rad)
     {
         bool If_Conact = false;

         double LocalCoordSystem[3][3];
         Compute2DimElementEdgeLocalSystem(EdgeCoord1, EdgeCoord2, Centroid, LocalCoordSystem);

         double dist = DistancePointToPlane(EdgeCoord1, LocalCoordSystem[2], Particle_Coord);

         if(dist < rad)
         {
             double IntersectionCoord[3] = {0.0};

             CoordProjectionOnPlane(Particle_Coord, EdgeCoord1, LocalCoordSystem, IntersectionCoord);

             double Vector1[3] = {0.0};
             double Vector2[3] = {0.0};

             Vector1[0] = IntersectionCoord[0] - EdgeCoord1[0];
             Vector1[1] = IntersectionCoord[1] - EdgeCoord1[1];
             Vector1[2] = IntersectionCoord[2] - EdgeCoord1[2];

             Vector2[0] = IntersectionCoord[0] - EdgeCoord2[0];
             Vector2[1] = IntersectionCoord[1] - EdgeCoord2[1];
             Vector2[2] = IntersectionCoord[2] - EdgeCoord2[2];

             norm(Vector1);
             norm(Vector2);

             if( DotProduct(Vector1, Vector2) <= 0.0)
             {
                 If_Conact = true;
             }

         }

         return If_Conact;
     }

     static inline bool JudgeIfThisEdgeIsContactWithParticle(double EdgeCoord1[3], double EdgeCoord2[3], double Centroid[3], double Particle_Coord[3], double rad,
                                                             double LocalCoordSystem[3][3], double Weight[2], double &DistPToB)
     {
         bool If_Conact = false;

         Compute2DimElementEdgeLocalSystem(EdgeCoord1, EdgeCoord2, Centroid, LocalCoordSystem);

         DistPToB = DistancePointToPlane(EdgeCoord1, LocalCoordSystem[2], Particle_Coord);

         if(DistPToB < rad)
         {
             double IntersectionCoord[3];
             CoordProjectionOnPlane(Particle_Coord, EdgeCoord1, LocalCoordSystem, IntersectionCoord);

             double dist1 = DistanceOfTwoPoint(EdgeCoord1, EdgeCoord2);
             double dist2 = DistanceOfTwoPoint(EdgeCoord1, IntersectionCoord);
             double dist3 = DistanceOfTwoPoint(EdgeCoord2, IntersectionCoord);

             if( dist2 <= dist1 && dist3 <= dist1 )
             {
                 If_Conact = true;
                 Weight[0] = dist3 / dist1;
                 Weight[1] = 1.0 - Weight[0];
             }

             /*
             double Vector1[3] = {0.0};
             double Vector2[3] = {0.0};

             Vector1[0] = IntersectionCoord[0] - EdgeCoord1[0];
             Vector1[1] = IntersectionCoord[1] - EdgeCoord1[1];
             Vector1[2] = IntersectionCoord[2] - EdgeCoord1[2];

             Vector2[0] = IntersectionCoord[0] - EdgeCoord2[0];
             Vector2[1] = IntersectionCoord[1] - EdgeCoord2[1];
             Vector2[2] = IntersectionCoord[2] - EdgeCoord2[2];

             norm(Vector1);
             norm(Vector2);

             if( DotProduct(Vector1, Vector2) <= 0.0)
             {
                 If_Conact = true;
             }

              */

         }

         return If_Conact;
     }


     static inline bool JudgeIfThisFaceIsContactWithParticle(int FaceNodeTotal, double Coord[4][3], double Centroid[3], double Particle_Coord[3], double rad)
     {
         bool If_Conact = false;

         double LocalCoordSystem[3][3];
         Compute3DimElementEdgeLocalSystem(Coord[0], Coord[1], Coord[2], Centroid, LocalCoordSystem);

         double dist = DistancePointToPlane(Coord[0], LocalCoordSystem[2], Particle_Coord);

         if(dist < rad)
         {
             double IntersectionCoord[3];
             CoordProjectionOnPlane(Particle_Coord, Coord[0], LocalCoordSystem, IntersectionCoord);


             if(FaceNodeTotal == 3)
             {
                 double Weight[3] = {0.0};
                 TriAngleWeight(Coord[0], Coord[1], Coord[2], IntersectionCoord, Weight);

                 if( fabs(Weight[0] + Weight[1] + Weight[2] - 1.0) < 1.0e-3 )
                 {
                     If_Conact = true;
                 }

             }
             else if(FaceNodeTotal == 4)
             {
                 double FaceArea;
                 double TempFace[2];
                 TriAngleArea(Coord[0], Coord[1], Coord[2], TempFace[0]);
                 TriAngleArea(Coord[2], Coord[3], Coord[0], TempFace[1]);
                 FaceArea = TempFace[0] + TempFace[1];

                 double AreaComponent[4] = {0.0};
                 TriAngleArea(IntersectionCoord, Coord[0], Coord[1], AreaComponent[0]);
                 TriAngleArea(IntersectionCoord, Coord[1], Coord[2], AreaComponent[1]);
                 TriAngleArea(IntersectionCoord, Coord[2], Coord[3], AreaComponent[2]);
                 TriAngleArea(IntersectionCoord, Coord[3], Coord[0], AreaComponent[3]);

                 if(fabs( (AreaComponent[0] + AreaComponent[1] + AreaComponent[2] + AreaComponent[3] - FaceArea) / FaceArea) < 1.0e-3)
                 {
                     If_Conact = true;
                 }
             }
         }

         return If_Conact;
     }


     ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     ///////////////******The four Functions BELOW are used to calculate the weight coefficient for quadrilateral*******///////////
     ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    static inline void Coord_transform(double origin[3],double coordsystem[3][3],double Coord[3],double TransCoord[3])
    {
	TransCoord[0]=0.0;
	TransCoord[1]=0.0;
	TransCoord[2]=0.0;
	double vector[3];
	vector[0]=Coord[0] - origin[0];
	vector[1]=Coord[1] - origin[1];
	vector[2]=Coord[2] - origin[2];
	TransCoord[0]=DotProduct(coordsystem[0],vector);
	TransCoord[1]=DotProduct(coordsystem[1],vector);
	TransCoord[2]=DotProduct(coordsystem[2],vector);

    }
    static inline void N44(double xp,double yp,double xy[4][2],double& N1,double& N2,double& N3,double& N4)
    {
	double xc=(xy[0][0]+xy[1][0]+xy[2][0]+xy[3][0])/4.0;
        double yc=(xy[0][1]+xy[1][1]+xy[2][1]+xy[3][1])/4.0;

	double elelength_x=2.0*fabs(xy[0][0]-xc);
	double elelength_y=2.0*fabs(xy[0][1]-yc);

	double Eps,Ita;
	Eps=2.0*(xp-xc)/elelength_x;
	Ita=2.0*(yp-yc)/elelength_y;
	N1=0.25*(1+Eps*2*(xy[0][0]-xc)/elelength_x)*(1+Ita*2*(xy[0][1]-yc)/elelength_y);
	N2=0.25*(1+Eps*2*(xy[1][0]-xc)/elelength_x)*(1+Ita*2*(xy[1][1]-yc)/elelength_y);
	N3=0.25*(1+Eps*2*(xy[2][0]-xc)/elelength_x)*(1+Ita*2*(xy[2][1]-yc)/elelength_y);
	N4=0.25*(1+Eps*2*(xy[3][0]-xc)/elelength_x)*(1+Ita*2*(xy[3][1]-yc)/elelength_y);
    }

    //Cfeng: irregular quadrilateral transfer to regular quadrilateral
    static inline void gl_to_iso(double x0,double y0, double xy[4][2],double &x_exisp,double &y_etasp)
    {
	double exisp=0.0;
	double etasp=0.0;
	double  tolerance=1.0e-8;
	double  x,y,x1,x2,x3,x4,y1,y2,y3,y4,a1,a2,a3,a4,b1,b2,b3,b4,s1,t1,d1,g1,g2,g0;
	x1 = xy[0][0];
	x2 = xy[1][0];
	x3 = xy[2][0];
	x4 = xy[3][0];
        y1 = xy[0][1];
	y2 = xy[1][1];
	y3 = xy[2][1];
	y4 = xy[3][1];
	a1=1.0/4*(-x1+x2+x3-x4);
	a2=1.0/4*(x1-x2+x3-x4);
	a3=1.0/4*(-x1-x2+x3+x4);
	a4=1.0/4*(x1+x2+x3+x4);
	b1=1.0/4*(-y1+y2+y3-y4);
	b2=1.0/4*(y1-y2+y3-y4);
	b3=1.0/4*(-y1-y2+y3+y4);
	b4=1.0/4*(y1+y2+y3+y4);

	x = x0 - a4;
	y = y0 - b4;
	s1 = a2*b3 - a3*b2;
        t1 = b2*x - a2*y + a1*b3 - a3*b1;
        d1 = b1*x - a1*y;

        if (fabs(s1) < tolerance)
	{
            etasp = -d1/t1;
            exisp = (x-a3*etasp) / (a1+a2*etasp);
	}
	else
	{
            g1 = (-t1 + sqrt(t1*t1-4*s1*d1)) / (2*s1);
            g2 = (-t1 - sqrt(t1*t1-4*s1*d1)) / (2*s1);
            if (fabs(g1) < 1.0+tolerance)
            {
                g0 = (x-a3*g1) / (a1+a2*g1);
                if (fabs(g0) < 1.0+tolerance)
                {
                    etasp = g1;
                    exisp = g0;
                }
            }
            if(fabs(g2) < 1.0+tolerance)
            {
                g0 = (x-a3*g2) / (a1+a2*g2);
                if(fabs(g0) < 1.0+tolerance)
                {
                   etasp = g2;
                   exisp = g0;
                }
            }

	}
	x_exisp=exisp;
	y_etasp=etasp;
    }

    static inline void CalQuadWeightCoefficient(double Coord[4][3], double LocalCoordSystem[3][3], double IntersectionCoord[3], double Weight[4])
    {
        int j;

        double FaceCenter[3] = {0.0};
        for(j = 0; j < 4; j++)
        {
            FaceCenter[0] += Coord[j][0] * 0.25;
            FaceCenter[1] += Coord[j][1] * 0.25;
            FaceCenter[2] += Coord[j][2] * 0.25;
        }

	double TransCoord0[3],TransCoord1[3],TransCoord2[3],TransCoord3[3];
	double xy[4][2];
	double xy1[4][2]={{-1.0,-1.0},{1.0,-1.0},{1.0,1.0},{-1.0,1.0}};


	double TempLocalCoordSystem[3][3]={{0.0}, {0.0}, {0.0}};
	double vx[3]={1.0,0,0},vy[3]={0,1.0,0},vz[3]={0, 0, 1.0};

	if( DotProduct(LocalCoordSystem[2],vx)<0 || DotProduct(LocalCoordSystem[2],vy)<0 || DotProduct(LocalCoordSystem[2],vz)<0 )
	{
            for(j=0;j<3;j++)
            {
                TempLocalCoordSystem[0][j] =  LocalCoordSystem[0][j];
                TempLocalCoordSystem[1][j] =  LocalCoordSystem[1][j];
                TempLocalCoordSystem[2][j] = -LocalCoordSystem[2][j];
            }
	}
	else
	{
            for(j=0;j<3;j++)
            {
                TempLocalCoordSystem[0][j] = LocalCoordSystem[0][j];
                TempLocalCoordSystem[1][j] = LocalCoordSystem[1][j];
                TempLocalCoordSystem[2][j] = LocalCoordSystem[2][j];
            }
	}

	Coord_transform(FaceCenter, TempLocalCoordSystem, Coord[0], TransCoord0);
	Coord_transform(FaceCenter, TempLocalCoordSystem, Coord[1], TransCoord1);
	Coord_transform(FaceCenter, TempLocalCoordSystem, Coord[2], TransCoord2);
	Coord_transform(FaceCenter, TempLocalCoordSystem, Coord[3], TransCoord3);

	xy[0][0] = TransCoord0[0], xy[0][1] = TransCoord0[1];
	xy[1][0] = TransCoord1[0], xy[1][1] = TransCoord1[1];
	xy[2][0] = TransCoord2[0], xy[2][1] = TransCoord2[1];
	xy[3][0] = TransCoord3[0], xy[3][1] = TransCoord3[1];

	double in0=0.0, in1=0.0, in2=0.0, in3=0.0;
	double TransCoordp[3];
	double Coordp_iso[2];

	Coord_transform(FaceCenter, TempLocalCoordSystem, IntersectionCoord, TransCoordp);

	gl_to_iso(TransCoordp[0],TransCoordp[1],xy,Coordp_iso[0],Coordp_iso[1]);

	N44(Coordp_iso[0],Coordp_iso[1], xy1, in0, in1, in2, in3);

	Weight[0]=in0;
	Weight[1]=in1;
	Weight[2]=in2;
	Weight[3]=in3;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     ///////////////******The four Functions ABOVE are used to calculate the weight coefficient for quadrilateral*******///////////
     ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


     static inline bool JudgeIfThisFaceIsContactWithParticle(int FaceNodeTotal, double Coord[4][3], double Centroid[3], double Particle_Coord[3], double rad,
                                                            double LocalCoordSystem[3][3], double Weight[4], double &DistPToB)
     {
         bool If_Conact = false;

         Compute3DimElementEdgeLocalSystem(Coord[0], Coord[1], Coord[2], Centroid, LocalCoordSystem);

         DistPToB = DistancePointToPlane(Coord[0], LocalCoordSystem[2], Particle_Coord);

         if(DistPToB < rad)
         {
             double IntersectionCoord[3];
             CoordProjectionOnPlane(Particle_Coord, Coord[0], LocalCoordSystem, IntersectionCoord);


             if(FaceNodeTotal == 3)
             {
                 double TriWeight[3] = {0.0};
                 TriAngleWeight(Coord[0], Coord[1], Coord[2], IntersectionCoord, TriWeight);

                 if( fabs(TriWeight[0] + TriWeight[1] + TriWeight[2] - 1.0) < 1.0e-3 )
                 {
                     Weight[0] = TriWeight[0];
                     Weight[1] = TriWeight[1];
                     Weight[2] = TriWeight[2];

                     If_Conact = true;
                 }

             }
             else if(FaceNodeTotal == 4)
             {
                 double FaceArea;
                 double TempFace[2];
                 TriAngleArea(Coord[0], Coord[1], Coord[2], TempFace[0]);
                 TriAngleArea(Coord[2], Coord[3], Coord[0], TempFace[1]);
                 FaceArea = TempFace[0] + TempFace[1];

                 double AreaComponent[4] = {0.0};
                 TriAngleArea(IntersectionCoord, Coord[0], Coord[1], AreaComponent[0]);
                 TriAngleArea(IntersectionCoord, Coord[1], Coord[2], AreaComponent[1]);
                 TriAngleArea(IntersectionCoord, Coord[2], Coord[3], AreaComponent[2]);
                 TriAngleArea(IntersectionCoord, Coord[3], Coord[0], AreaComponent[3]);

                 if(fabs( (AreaComponent[0] + AreaComponent[1] + AreaComponent[2] + AreaComponent[3] - FaceArea) / FaceArea) < 1.0e-3)
                 {
                     If_Conact = true;

                     CalQuadWeightCoefficient(Coord, LocalCoordSystem, IntersectionCoord, Weight);
                 }
             }
         }

         return If_Conact;
     }


    }
}
#endif	/* _GEOMETRYFUNCTIONS_H */

