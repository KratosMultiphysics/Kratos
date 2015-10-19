/* 
 * File:   GeometryFunctions.h
 * Author: msantasusana, Chun Feng
 *
 * Created on 21 de mayo de 2012, 19:40
 */

#ifndef _GEOMETRYFUNCTIONS_H
#define	_GEOMETRYFUNCTIONS_H


#include <cmath>

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

    static inline void normalize(double Vector[3])
    {
            double distance = sqrt(Vector[0] * Vector[0] + Vector[1] * Vector[1] + Vector[2] * Vector[2]);
            
            double inv_distance = (distance != 0.0) ?  1.0 / distance : 0.00;
            Vector[0] = Vector[0] * inv_distance;
            Vector[1] = Vector[1] * inv_distance;
            Vector[2] = Vector[2] * inv_distance;
    }
    
    static inline void normalize( array_1d<double,3>& Vector, double& distance)
    {
            distance = sqrt(Vector[0] * Vector[0] + Vector[1] * Vector[1] + Vector[2] * Vector[2]);
            
            double inv_distance = (distance != 0.0) ?  1.0 / distance : 0.00;
            Vector[0] = Vector[0] * inv_distance;
            Vector[1] = Vector[1] * inv_distance;
            Vector[2] = Vector[2] * inv_distance;       
    }
    
    static inline void normalize( double Vector[3], double& distance)
    {
            distance = sqrt(Vector[0] * Vector[0] + Vector[1] * Vector[1] + Vector[2] * Vector[2]);
            
            double inv_distance = (distance != 0.0) ?  1.0 / distance : 0.00;
            Vector[0] = Vector[0] * inv_distance;
            Vector[1] = Vector[1] * inv_distance;
            Vector[2] = Vector[2] * inv_distance;       
    }
     
    static inline void normalize( array_1d<double,3>& Vector)
    {
            double distance = sqrt(Vector[0] * Vector[0] + Vector[1] * Vector[1] + Vector[2] * Vector[2]);
            
            double inv_distance = (distance != 0.0) ?  1.0 / distance : 0.00;
            Vector[0] = Vector[0] * inv_distance;
            Vector[1] = Vector[1] * inv_distance;
            Vector[2] = Vector[2] * inv_distance;
    }     
      
    static inline void module( const array_1d<double,3>& Vector, double& distance)
    {
            distance = sqrt(Vector[0] * Vector[0] + Vector[1] * Vector[1] + Vector[2] * Vector[2]);
    }
    
     static inline double module( const double Vector[3])
    {
            return sqrt(Vector[0] * Vector[0] + Vector[1] * Vector[1] + Vector[2] * Vector[2]);
    }
    
     static inline void module( const double Vector[3], double& distance)
    {
            distance = sqrt(Vector[0] * Vector[0] + Vector[1] * Vector[1] + Vector[2] * Vector[2]);
    }
    
    static inline double module( const array_1d<double,3>& Vector)
    {
            double distance = sqrt(Vector[0] * Vector[0] + Vector[1] * Vector[1] + Vector[2] * Vector[2]);
            return distance;
    }

    static inline void VectorGlobal2Local(const double LocalCoordSystem[3][3], const double GlobalVector[3], double LocalVector[3])
    {
        for (int i=0; i<3; i++) {
            LocalVector[i] = 0.0;
            for (int j=0; j<3; j++) {
                LocalVector[i]+=LocalCoordSystem[i][j]*GlobalVector[j];
            }
        }
    }
    static inline void VectorGlobal2Local(const double LocalCoordSystem[3][3], const array_1d<double, 3 > & GlobalVector, array_1d<double, 3 > & LocalVector)
    {
        for (int i=0; i<3; i++) {
            LocalVector[i] = 0.0;
            for (int j=0; j<3; j++) {
                LocalVector[i]+=LocalCoordSystem[i][j]*GlobalVector[j];
            }
        }
    }
    
    static inline void VectorGlobal2Local(const double LocalCoordSystem[3][3], const array_1d<double, 3 > & GlobalVector, double LocalVector[3])
    {
        for (int i=0; i<3; i++) {
            LocalVector[i] = 0.0;
            for (int j=0; j<3; j++) {
                LocalVector[i]+=LocalCoordSystem[i][j]*GlobalVector[j];
            }
        }
    }
        
    static inline void VectorLocal2Global(const double LocalCoordSystem[3][3], const double LocalVector[3], double GlobalVector[3])
    {
        for (int i=0; i<3; i++) {
            GlobalVector[i] = 0.0;
            for (int j=0; j<3; j++) {
                GlobalVector[i]+=LocalCoordSystem[j][i]*LocalVector[j];
            }
        }                
    }
    
    static inline void VectorLocal2Global(const double LocalCoordSystem[3][3], const array_1d<double, 3 > & LocalVector, array_1d<double, 3 > & GlobalVector)
    {
        for (int i=0; i<3; i++) {
            GlobalVector[i] = 0.0;
            for (int j=0; j<3; j++) {
                GlobalVector[i]+=LocalCoordSystem[j][i]*LocalVector[j];
            }
        }                
    }
    
    static inline void VectorLocal2Global(const double LocalCoordSystem[3][3], const array_1d<double, 3 > & LocalVector, double GlobalVector[3])
    {
        for (int i=0; i<3; i++) {
            GlobalVector[i] = 0.0;
            for (int j=0; j<3; j++) {
                GlobalVector[i]+=LocalCoordSystem[j][i]*LocalVector[j];
            }
        }                
    }
    

    static inline double DotProduct(double Vector1[3], double Vector2[3])
    {
        return Vector1[0] * Vector2[0] + Vector1[1] * Vector2[1] + Vector1[2] * Vector2[2];
    }

    static inline double DotProduct(const array_1d<double,3>& Vector1, const array_1d<double,3>& Vector2)
    {
        return Vector1[0] * Vector2[0] + Vector1[1] * Vector2[1] + Vector1[2] * Vector2[2];
    }

    static inline void CrossProduct(const double u[3], const double v[3], double ReturnVector[3])
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
    
    static inline void CrossProduct( const double u[3], const array_1d<double,3>& v, double ReturnVector[3])
    {
    	ReturnVector[0] = u[1]*v[2] - u[2]*v[1];
        ReturnVector[1] = v[0]*u[2] - u[0]*v[2];
        ReturnVector[2] = u[0]*v[1] - u[1]*v[0];
    }
    
    static inline void CrossProduct( const array_1d<double,3>& u, const double v[3], double ReturnVector[3])
    {
    	ReturnVector[0] = u[1]*v[2] - u[2]*v[1];
        ReturnVector[1] = v[0]*u[2] - u[0]*v[2];
        ReturnVector[2] = u[0]*v[1] - u[1]*v[0];
    }
    
    static inline void CrossProduct( const array_1d<double,3>& u, const double v[3], array_1d<double,3>& ReturnVector)
    {
    	ReturnVector[0] = u[1]*v[2] - u[2]*v[1];
        ReturnVector[1] = v[0]*u[2] - u[0]*v[2];
        ReturnVector[2] = u[0]*v[1] - u[1]*v[0];
    }
    
    static inline void CrossProduct( const array_1d<double,3>& u, const array_1d<double,3>& v,  double ReturnVector[3])
    {
    	ReturnVector[0] = u[1]*v[2] - u[2]*v[1];
        ReturnVector[1] = v[0]*u[2] - u[0]*v[2];
        ReturnVector[2] = u[0]*v[1] - u[1]*v[0];
    }

   //NOTE:: Modified by M. Santasusana Feb 2013 - simplification (the one proposed by F.Chun was for a more generalized case) 
    static inline void ComputeContactLocalCoordSystem(array_1d<double,3>& NormalDirection, const double& distance, double LocalCoordSystem[3][3])  //inline: modifies the LocalCoordSystem as it were a reference
    {
           
        double inv_distance = (distance != 0.0) ?  1.0 / distance : 0.00;
        NormalDirection[0] *= inv_distance;
        NormalDirection[1] *= inv_distance;
        NormalDirection[2] *= inv_distance;                      
       
      if(fabs(NormalDirection[0])>=0.577)
        {
            LocalCoordSystem[0][0]= - NormalDirection[1];
            LocalCoordSystem[0][1]= NormalDirection[0];
            LocalCoordSystem[0][2]= 0.0;
        }
        else if(fabs(NormalDirection[1])>=0.577)
        {
            LocalCoordSystem[0][0]= 0.0;
            LocalCoordSystem[0][1]= - NormalDirection[2];
            LocalCoordSystem[0][2]= NormalDirection[1];
        }        
        else
        {                   
            LocalCoordSystem[0][0]= NormalDirection[2];
            LocalCoordSystem[0][1]= 0.0;
            LocalCoordSystem[0][2]= - NormalDirection[0];
        }
       
        //normalize(Vector0);
        double distance0 = sqrt(LocalCoordSystem[0][0] * LocalCoordSystem[0][0] + LocalCoordSystem[0][1] * LocalCoordSystem[0][1] + LocalCoordSystem[0][2] * LocalCoordSystem[0][2]);
        double inv_distance0 = (distance0 != 0.0) ?  1.0 / distance0 : 0.00;
        LocalCoordSystem[0][0] = LocalCoordSystem[0][0] * inv_distance0;
        LocalCoordSystem[0][1] = LocalCoordSystem[0][1] * inv_distance0;
        LocalCoordSystem[0][2] = LocalCoordSystem[0][2] * inv_distance0;
        
        //CrossProduct(NormalDirection,Vector0,Vector1);
        LocalCoordSystem[1][0] = NormalDirection[1]*LocalCoordSystem[0][2] - NormalDirection[2]*LocalCoordSystem[0][1];
        LocalCoordSystem[1][1] = LocalCoordSystem[0][0]*NormalDirection[2] - NormalDirection[0]*LocalCoordSystem[0][2];
        LocalCoordSystem[1][2] = NormalDirection[0]*LocalCoordSystem[0][1] - NormalDirection[1]*LocalCoordSystem[0][0];

        //normalize(Vector1);
        
        LocalCoordSystem[2][0]=NormalDirection[0];
        LocalCoordSystem[2][1]=NormalDirection[1];
        LocalCoordSystem[2][2]=NormalDirection[2];        

    }
   
     static inline double DistanceOfTwoPoint(double coord1[3], double coord2[3])
     {
         double dx = coord1[0] - coord2[0];
         double dy = coord1[1] - coord2[1];
         double dz = coord1[2] - coord2[2];

         return sqrt(dx * dx + dy * dy + dz * dz);
     }
     
     static inline double DistanceOfTwoPointSquared(double coord1[3], double coord2[3])
    {
        double dx = coord1[0] - coord2[0];
        double dy = coord1[1] - coord2[1];
        double dz = coord1[2] - coord2[2];

        return (dx * dx + dy * dy + dz * dz);
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

     //TriAngle Weight, coord1,coord2,coord3,testcoord,weight
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
     
     //Quadrilatera Weight, coord1,coord2,coord3,testcoord,weight (Paper Zhang)
     static inline void QuadAngleWeight(double Coord1[3], double Coord2[3], double Coord3[3], double Coord4[3], double JudgeCoord[3], double Weight[4])
     {
         double area[4], s1, s2, s;
         TriAngleArea(Coord1, Coord2, JudgeCoord, area[0]);
         TriAngleArea(Coord2, Coord3, JudgeCoord, area[1]);
         TriAngleArea(Coord3, Coord4, JudgeCoord, area[2]);
         TriAngleArea(Coord4, Coord1, JudgeCoord, area[3]);
         
         TriAngleArea(Coord1, Coord2, Coord3, s1);//msimsi
         TriAngleArea(Coord1, Coord3, Coord4, s2);//msimsi
         
         s = s1 + s2;

         if( fabs(area[0] + area[1] + area[2] + area[3] - s) < 1.0e-15 ) //msimsi
         {
             double QuadNormArea = 1 / ((area[0] + area[2]) * (area[1] + area[3]));

             Weight[0] = (area[1] * area[2]) * QuadNormArea;
             Weight[1] = (area[2] * area[3]) * QuadNormArea;
             Weight[2] = (area[3] * area[0]) * QuadNormArea;
             Weight[3] = (area[0] * area[1]) * QuadNormArea;
         }
     }

     static inline double DistancePointToPlane(double CoordInPlane[3], double PlaneUnitNormalVector[3], double TestCoord[3])
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


     //MSI::NOT USED
     static inline void Compute2DimElementEdgeLocalSystem(double EdgeCoord1[3], double EdgeCoord2[3], double ParticleCoord[3], double LocalCoordSystem[3][3])
     {
         double Vector1[3] = {0.0};
         double Vector2[3] = {0.0};
         double Vector3[3] = {0.0};
         double Normal [3] = {0.0};

         Vector1[0] = EdgeCoord2[0] - EdgeCoord1[0];
         Vector1[1] = EdgeCoord2[1] - EdgeCoord1[1];
         Vector1[2] = EdgeCoord2[2] - EdgeCoord1[2];

         Vector2[0] = ParticleCoord[0] - EdgeCoord1[0];
         Vector2[1] = ParticleCoord[1] - EdgeCoord1[1];
         Vector2[2] = ParticleCoord[2] - EdgeCoord1[2];

         normalize(Vector1);
         normalize(Vector2);
         CrossProduct(Vector1, Vector2, Vector3);
         normalize(Vector3);
         CrossProduct(Vector1, Vector3, Normal);
         normalize(Normal);

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

    static inline void Compute3DimElementFaceLocalSystem(double FaceCoord1[3], double FaceCoord2[3], double FaceCoord3[3], double ParticleCoord[3], 
                                                         double LocalCoordSystem[3][3], double& normal_flag){
         
        //NOTE: this function is designed in a way that the normal always points the side where the centre of particle is found. Therefore should only be used in this way if the indentation is less than the radius value.
                //the fucntion returns a flag with the same value as the dotproduct of the normal of the triangle and the normal pointing to the particle.
         
         double Vector1[3] = {0.0};
         double Vector2[3] = {0.0};
         double Vector3[3] = {0.0};
         double Normal[3]  = {0.0};
         
         Vector1[0] = FaceCoord2[0] - FaceCoord1[0];
         Vector1[1] = FaceCoord2[1] - FaceCoord1[1];
         Vector1[2] = FaceCoord2[2] - FaceCoord1[2];

         Vector2[0] = FaceCoord3[0] - FaceCoord2[0];
         Vector2[1] = FaceCoord3[1] - FaceCoord2[1];
         Vector2[2] = FaceCoord3[2] - FaceCoord2[2];

         normalize(Vector1);
         CrossProduct(Vector1, Vector2, Normal);
         normalize(Normal);

         CrossProduct(Normal, Vector1, Vector2);
         normalize(Vector2);

         Vector3[0] = ParticleCoord[0] - FaceCoord1[0];
         Vector3[1] = ParticleCoord[1] - FaceCoord1[1];
         Vector3[2] = ParticleCoord[2] - FaceCoord1[2];
         normalize(Vector3);



         if (DotProduct(Vector3, Normal) > 0.0) 
         {
             for (int ia = 0; ia < 3; ia++)
             {
                    normal_flag             = 1.0;
                    LocalCoordSystem[0][ia] = Vector1[ia];
                    LocalCoordSystem[1][ia] = Vector2[ia];
                    LocalCoordSystem[2][ia] = Normal [ia];
             }
         }
         else
         {
             for (int ia = 0; ia < 3; ia++)
             {
                    normal_flag             = -1.0;
                    LocalCoordSystem[0][ia] = -Vector1[ia];
                    LocalCoordSystem[1][ia] = -Vector2[ia];
                    LocalCoordSystem[2][ia] = -Normal [ia];
             }
         }
     }//Compute3DimElementFaceLocalSystem

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

        if (L==0)
        {

           // KRATOS_WATCH("WARNING! Null vector for the rotation!")

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
         
         normalize(temp);
         
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

     ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     ///////////////******QuickDistanceForAKnownNeighbour******/////////////////////////////////////////////////////////////////////
     ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
/*    static inline void QuickDistanceForAKnownNeighbour(double Coord[4][3], double Particle_Coord[3],  double &DistPToB) //M.Santasusana
    {
     
          double Vector1[3] = {0.0};
          double Vector2[3] = {0.0};
          
          double Normal [3] = {0.0};
          
          Vector1[0] =  Coord[1][0] - Coord[0][0];
          Vector1[1] =  Coord[1][1] - Coord[0][1];
          Vector1[2] =  Coord[1][2] - Coord[0][2];

          Vector2[0] = Coord[2][0] - Coord[1][0];
          Vector2[1] = Coord[2][1] - Coord[1][1];
          Vector2[2] = Coord[2][2] - Coord[1][2];

          CrossProduct(Vector1, Vector2, Normal);
          
          normalize(Normal);

          DistPToB = DistancePointToPlane(Coord[0], Normal, Particle_Coord);


     } 
    
    static inline void QuickDistanceForAKnownEdgeNeighbour(double Coord[4][3], double Particle_Coord[3],  double &DistPToB) //M.Santasusana
    {
     
                 double Vector1[3] = {0.0};
                 double Vector2[3] = {0.0};
                 double Vector3[3] = {0.0};
                 double NormalV[3] = {0.0};
                 
                 Vector1[0] = Coord[1][0] - Coord[0][0];
                 Vector1[1] = Coord[1][1] - Coord[0][1];
                 Vector1[2] = Coord[1][2] - Coord[0][2];
                 //normalize(Vector1);
                 
                 Vector2[0] = Particle_Coord[0] - Coord[0][0];
                 Vector2[1] = Particle_Coord[1] - Coord[0][1];
                 Vector2[2] = Particle_Coord[2] - Coord[0][2];
                 //normalize(Vector2);
                 
                 CrossProduct(Vector1, Vector2, Vector3); 
                 //normalize(Vector3);
                 
                 CrossProduct(Vector3, Vector1, NormalV);
                 //normalize(NormalV);
         
                 if(DotProduct(NormalV,Vector2) < 0)
                 { 
                         NormalV[0] = -NormalV[0];
                         NormalV[1] = -NormalV[1]; 
                         NormalV[2] = -NormalV[2];
                 }
         
         normalize(NormalV);

        DistPToB = DistancePointToPlane(Coord[0], NormalV, Particle_Coord);
     }
     
    static inline void QuickDistanceForAKnownPointNeighbour(double Coord[4][3], double Particle_Coord[3],  double &DistPToB) //J.Irazábal
    {

        DistPToB = DistanceOfTwoPoint(Coord[0], Particle_Coord);

     }
*/

     ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     ///////////////******JudgeIfThisPointIsContactWithParticle******///////////////////////////////////////////////////////////////
     ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    static inline bool JudgeIfThisPointIsContactWithParticle(double PointCoord[3], double Particle_Coord[3], double rad, double LocalCoordSystem[3][3], double &DistPToB)
    {
        bool If_Contact = false;

        DistPToB = DistanceOfTwoPoint(PointCoord, Particle_Coord);

        if(DistPToB < rad)
        {
            If_Contact = true;

            //double Vector1[3] = {0.0};
            array_1d<double,3> Vector1;
            Vector1[0] = Particle_Coord[0] - PointCoord[0];
            Vector1[1] = Particle_Coord[1] - PointCoord[1];
            Vector1[2] = Particle_Coord[2] - PointCoord[2];
            normalize(Vector1);

            ComputeContactLocalCoordSystem(Vector1, 1.0, LocalCoordSystem); 
        }
        return If_Contact;

    }
    
    static inline bool JudgeIfThisPointIsContactWithParticle(double PointCoord[3], double Particle_Coord[3], double rad, double &dist)
    {
        bool If_Contact = false;
                
        dist = DistanceOfTwoPoint(PointCoord, Particle_Coord);
                
        if(dist < rad)
        {
            If_Contact = true;
                
        }
        
        return If_Contact;
    }
    
     ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     ///////////////******JudgeIfThisEdgeIsContactWithParticle******////////////////////////////////////////////////////////////////
     ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
    static inline bool JudgeIfThisEdgeIsContactWithParticle(double EdgeCoord1[3], double EdgeCoord2[3], double Particle_Coord[3], double rad, double &dist)
    {
         bool If_Contact = false;

         double LocalCoordSystem[3][3];
         
         double Vector1[3] = {0.0};
         double Vector2[3] = {0.0};
         double Vector3[3] = {0.0};
         double NormalV[3] = {0.0};
 
         Vector1[0] = EdgeCoord2[0] - EdgeCoord1[0];
         Vector1[1] = EdgeCoord2[1] - EdgeCoord1[1];
         Vector1[2] = EdgeCoord2[2] - EdgeCoord1[2];
         normalize(Vector1);

         Vector2[0] = Particle_Coord[0] - EdgeCoord1[0];
         Vector2[1] = Particle_Coord[1] - EdgeCoord1[1];
         Vector2[2] = Particle_Coord[2] - EdgeCoord1[2];
         normalize(Vector2);

         CrossProduct(Vector1, Vector2, Vector3); 
         normalize(Vector3);

         CrossProduct(Vector3, Vector1, NormalV);
         normalize(NormalV);

         LocalCoordSystem[0][0] = Vector3[0];
         LocalCoordSystem[0][1] = Vector3[1]; 
         LocalCoordSystem[0][2] = Vector3[2];
         
         LocalCoordSystem[1][0] = Vector1[0];
         LocalCoordSystem[1][1] = Vector1[1]; 
         LocalCoordSystem[1][2] = Vector1[2];

         LocalCoordSystem[2][0] = NormalV[0];
         LocalCoordSystem[2][1] = NormalV[1]; 
         LocalCoordSystem[2][2] = NormalV[2];

         if(DotProduct(NormalV,Vector2) < 0)
         {
             LocalCoordSystem[0][0] = -Vector3[0];
             LocalCoordSystem[0][1] = -Vector3[1]; 
             LocalCoordSystem[0][2] = -Vector3[2];
             
             LocalCoordSystem[1][0] = -Vector1[0];
             LocalCoordSystem[1][1] = -Vector1[1]; 
             LocalCoordSystem[1][2] = -Vector1[2];
             
             LocalCoordSystem[2][0] = -NormalV[0];
             LocalCoordSystem[2][1] = -NormalV[1]; 
             LocalCoordSystem[2][2] = -NormalV[2];
         }
         
         dist = DistancePointToPlane(EdgeCoord1, LocalCoordSystem[2], Particle_Coord);

         if(dist < rad)
         {
             double IntersectionCoord[3] = {0.0};

             CoordProjectionOnPlane(Particle_Coord, EdgeCoord1, LocalCoordSystem, IntersectionCoord);
             
             double dist1 = DistanceOfTwoPoint(EdgeCoord1, EdgeCoord2);
             double dist2 = DistanceOfTwoPoint(EdgeCoord1, IntersectionCoord);
             double dist3 = DistanceOfTwoPoint(EdgeCoord2, IntersectionCoord);

             if( fabs(((dist2 + dist3)/dist1) - 1.0) < 1.0e-15 )
             {
                 If_Contact = true;
             } 

//              Vector1[0] = IntersectionCoord[0] - EdgeCoord1[0];
//              Vector1[1] = IntersectionCoord[1] - EdgeCoord1[1];
//              Vector1[2] = IntersectionCoord[2] - EdgeCoord1[2];
// 
//              Vector2[0] = IntersectionCoord[0] - EdgeCoord2[0];
//              Vector2[1] = IntersectionCoord[1] - EdgeCoord2[1];
//              Vector2[2] = IntersectionCoord[2] - EdgeCoord2[2];
// 
//              normalize(Vector1);
//              normalize(Vector2);
// 
//              if( DotProduct(Vector1, Vector2) <= 0.0)
//              {
//                  If_Contact = true;
//              }

        }

        return If_Contact;
    }

     static inline bool JudgeIfThisEdgeIsContactWithParticle(double EdgeCoord1[3], double EdgeCoord2[3], double Particle_Coord[3], double rad,
                                                             double LocalCoordSystem[3][3], double Weight[2], double &DistPToB)
     {
         bool If_Contact = false;

         double Vector1[3] = {0.0};
         double Vector2[3] = {0.0};
         double Vector3[3] = {0.0};
         double NormalV[3] = {0.0};
 
         Vector1[0] = EdgeCoord2[0] - EdgeCoord1[0];
         Vector1[1] = EdgeCoord2[1] - EdgeCoord1[1];
         Vector1[2] = EdgeCoord2[2] - EdgeCoord1[2];
         normalize(Vector1);

         Vector2[0] = Particle_Coord[0] - EdgeCoord1[0];
         Vector2[1] = Particle_Coord[1] - EdgeCoord1[1];
         Vector2[2] = Particle_Coord[2] - EdgeCoord1[2];
         //normalize(Vector2);//useess

         CrossProduct(Vector1, Vector2, Vector3); 
         normalize(Vector3);

         CrossProduct(Vector3, Vector1, NormalV);
         normalize(NormalV);

         LocalCoordSystem[0][0] = Vector3[0];
         LocalCoordSystem[0][1] = Vector3[1]; 
         LocalCoordSystem[0][2] = Vector3[2];

         LocalCoordSystem[1][0] = Vector1[0];
         LocalCoordSystem[1][1] = Vector1[1]; 
         LocalCoordSystem[1][2] = Vector1[2];

         LocalCoordSystem[2][0] = NormalV[0];
         LocalCoordSystem[2][1] = NormalV[1]; 
         LocalCoordSystem[2][2] = NormalV[2];
 
         if(DotProduct(NormalV,Vector2) < 0)
         {
             LocalCoordSystem[0][0] = -Vector3[0];
             LocalCoordSystem[0][1] = -Vector3[1]; 
             LocalCoordSystem[0][2] = -Vector3[2];

             LocalCoordSystem[1][0] = -Vector1[0];
             LocalCoordSystem[1][1] = -Vector1[1]; 
             LocalCoordSystem[1][2] = -Vector1[2];

             LocalCoordSystem[2][0] = -NormalV[0];
             LocalCoordSystem[2][1] = -NormalV[1];
             LocalCoordSystem[2][2] = -NormalV[2];
         }

         DistPToB = DistancePointToPlane(EdgeCoord1, LocalCoordSystem[2], Particle_Coord);

         if(DistPToB < rad)
         {
             double IntersectionCoord[3];
             CoordProjectionOnPlane(Particle_Coord, EdgeCoord1, LocalCoordSystem, IntersectionCoord);

             double dist1 = DistanceOfTwoPoint(EdgeCoord1, EdgeCoord2);
             double dist2 = DistanceOfTwoPoint(EdgeCoord1, IntersectionCoord);
             double dist3 = DistanceOfTwoPoint(EdgeCoord2, IntersectionCoord);

             if( fabs(((dist2 + dist3)/dist1) - 1.0) < 1.0e-15 )
             {
                 If_Contact = true;
                 Weight[0] = dist3 / dist1;
                 Weight[1] = 1.0 - Weight[0];
             }
         }

         return If_Contact;
     }

   
    static inline void TauSqCalculation(double& tau_sq, double d0, double d1, double d2)
    {
        // Basically is Heron´s formula squared for two times the Area http://www.mathopenref.com/heronsformula.html
        tau_sq = 0.25*( d0 + d1 - d2)*( d0 - d1 + d2)*( -d0 + d1 + d2)*( d0 + d1 + d2);
        
    }
     
    static inline bool QuickClosestEdgeVertexDetermination(/*int& ContactType,*/ double Coord[4][3], double Particle_Coord[3], int size, double particle_radius) 
    {                  
        
        std::vector<double> dist_sq(size*2, 0.0); //vertices and edges of the face       
        std::vector<double> dist(size*2, 0.0);
        double base[4][3];
        double base_sq    = 0.0; 
        double tau_sq     = 0.0;
        double base_leng  = 0.0;
        double diag_back  = 0.0;
        double diag_front = 0.0;
        double particle_radius_sq = particle_radius*particle_radius;
        
        double NodeToP[4][3] = { {0.0},{0.0},{0.0},{0.0} };
        
        //VERTICES
        for (int n=0;n<size;n++){ //vertices
        
            
            for (int k= 0; k<3; k++){ //components
            
                NodeToP[n][k] =  Particle_Coord[k] - Coord[n][k];
                dist_sq[n]    += NodeToP[n][k]*NodeToP[n][k];
                
            }
            
            if ( dist_sq[n] <= particle_radius_sq )
            {
                //ContactType = 3;
                return true; 
            }
            
            dist[n] = sqrt(dist_sq[n]);
            
        }
       
        //EDGES
        for (int i=0;i<size;i++)
        {
           
            int j=(i+1)%size;
            base_sq = 0.0;
            for (int k=0;k<3;k++)
            {
                base[i][k] = (Coord[j][k] - Coord[i][k]); 
                base_sq += base[i][k]*base[i][k];
            
            }
                    
            base_leng = sqrt(base_sq);
            diag_back  = dist[i];
            diag_front = dist[j];
            
            double a2  =  dist_sq[i];
            double b2  =  base_sq;
            double c2  =  dist_sq[j];
            
            
            if ( (a2 <= b2 + c2) && (c2 <= a2 + b2) ){  //neglecting obtuse triangles.
                
                TauSqCalculation(tau_sq,base_leng,diag_back,diag_front);
                
                dist_sq[i+size] = tau_sq/base_sq; 
                
                if ( dist_sq[i+size] <= particle_radius_sq )
                {
                   
                    //ContactType=2;
                    return true;
                    
                }
            
            }
            
        }//for each edge starting at node i going to j
        
        return false;
     
    } //QuickClosestEdgeVertexDetermination

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////******JudgeIfThisFaceIsContactWithParticle******////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    static inline bool JudgeIfThisFaceIsContactWithParticle(int FaceNodeTotal, double Coord[4][3], double Particle_Coord[3], double rad,
                                                            double LocalCoordSystem[3][3], double Weight[4], double &DistPToB)
    {
 
        bool If_Contact = false;
      
        double dummy_normal_flag = 0.0;
        Compute3DimElementFaceLocalSystem(Coord[0], Coord[1], Coord[2], Particle_Coord, LocalCoordSystem, dummy_normal_flag);

        DistPToB = DistancePointToPlane(Coord[0], LocalCoordSystem[2], Particle_Coord);

        if (DistPToB < rad)
        {
            double IntersectionCoord[3];
            CoordProjectionOnPlane(Particle_Coord, Coord[0], LocalCoordSystem, IntersectionCoord);
           
            if (FaceNodeTotal == 3)
            {
                double TriWeight[3] = {0.0};
                TriAngleWeight(Coord[0], Coord[1], Coord[2], IntersectionCoord, TriWeight);
                      
                if (fabs(TriWeight[0] + TriWeight[1] + TriWeight[2] - 1.0) < 1.0e-15)
                {
                    Weight[0] = TriWeight[0];
                    Weight[1] = TriWeight[1];
                    Weight[2] = TriWeight[2];
 
                    If_Contact = true;
                }
            }
            
            if (FaceNodeTotal == 4)
            {
                double QuadWeight[4] = {0.0};
                QuadAngleWeight(Coord[0], Coord[1], Coord[2],  Coord[3], IntersectionCoord, QuadWeight);
                      
                if (fabs(QuadWeight[0] + QuadWeight[1] + QuadWeight[2] + QuadWeight[3] - 1.0) < 1.0e-15)
                {
                    Weight[0] = QuadWeight[0];
                    Weight[1] = QuadWeight[1];
                    Weight[2] = QuadWeight[2];
                    Weight[3] = QuadWeight[3];
               
                    If_Contact = true;
                }
            }
        }
 
        return If_Contact;
          
    } //JudgeIfThisFaceIsContactWithParticle;
     
    //M.Santasusana Februar 2015   
     
    static inline  bool SameSide(double Coord1[3], double Coord2[3], double ReferenceCoord[3], double JudgeCoord[3]){
        
        double cp1[3]  = {0.0};
        double cp2[3]  = {0.0};
        double b_a[3]  = {0.0};
        double p1_a[3] = {0.0};
        double p2_a[3] = {0.0};
        
        for (int i=0;i<3;i++){
            
            b_a[i]  = Coord2[i]-Coord1[i];
            p1_a[i] = JudgeCoord[i]-Coord1[i];
            p2_a[i] = ReferenceCoord[i]-Coord1[i];
            
        }
        
        GeometryFunctions::CrossProduct(b_a, p1_a, cp1);
        GeometryFunctions::CrossProduct(b_a, p2_a, cp2);
        
        if (GeometryFunctions::DotProduct(cp1, cp2) >= 0) return true;
            
        else return false;
        
    }//SameSide
    
    static inline  bool Fast_PointInsideTriangle(double Coord1[3], double Coord2[3], double Coord3[3], double JudgeCoord[3]){
        
        if( SameSide(Coord1, Coord2, Coord3, JudgeCoord) == false ) return false;
        if( SameSide(Coord2, Coord3, Coord1, JudgeCoord) == false ) return false;
        if( SameSide(Coord3, Coord1, Coord2, JudgeCoord) == false ) return false;
                               
        return true;    
        
    }
    
    static inline bool Fast_JudgeIfThisFaceIsContactWithParticle(int FaceNodeTotal, double Coord[4][3], double Particle_Coord[3], double rad,
                                                                double LocalCoordSystem[3][3], double &DistPToB)
    {
 
        bool If_Contact = false;

        double dummy_normal_flag = 0.0;
        Compute3DimElementFaceLocalSystem(Coord[0], Coord[1], Coord[2], Particle_Coord, LocalCoordSystem, dummy_normal_flag);

        DistPToB = DistancePointToPlane(Coord[0], LocalCoordSystem[2], Particle_Coord);
 
        if(DistPToB < rad )
        {

            double IntersectionCoord[3];
            CoordProjectionOnPlane(Particle_Coord, Coord[0], LocalCoordSystem, IntersectionCoord);
            
            If_Contact = Fast_PointInsideTriangle(Coord[0], Coord[1], Coord[2], IntersectionCoord);                  
 
            if( (FaceNodeTotal == 4) && (If_Contact==false)) //maybe not found inside first triangle, then check second triangular half of the quadrilateral
            {
            
                bool If_Contact2 = Fast_PointInsideTriangle(Coord[0], Coord[2], Coord[3], IntersectionCoord);
                
                if( (If_Contact == true) || (If_Contact2 == true) ) { If_Contact = true; }

            }
            
                        
        }
 
        return If_Contact;
        
    } //Fast_JudgeIfThisFaceIsContactWithParticle;
    
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

    /*static void CalQuadWeightCoefficient(double Coord[4][3], double LocalCoordSystem[3][3], double IntersectionCoord[3], double Weight[4])
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
        

        xy[0][0] = TransCoord0[0]; xy[0][1] = TransCoord0[1];
        xy[1][0] = TransCoord1[0]; xy[1][1] = TransCoord1[1];
        xy[2][0] = TransCoord2[0]; xy[2][1] = TransCoord2[1];
        xy[3][0] = TransCoord3[0]; xy[3][1] = TransCoord3[1];

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
        
        
        
    }*/

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     ///////////////******The four Functions ABOVE are used to calculate the weight coefficient for quadrilateral*******///////////
     ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    
     static inline void GetRotationMatrix(const array_1d<double, 3 > & EulerAngles, double rotation_matrix[3][3]){                    
         
            double cosA=cos(EulerAngles[0]);
            double sinA=sin(EulerAngles[0]);
            double cosB=cos(EulerAngles[1]);
            double sinB=sin(EulerAngles[1]);
            double cosC=cos(EulerAngles[2]);
            double sinC=sin(EulerAngles[2]);
         
            rotation_matrix[0][0] = cosC*cosA - cosB*sinA*sinC;
            rotation_matrix[0][1] = -sinC*cosA - cosB*sinA*cosC;
            rotation_matrix[0][2] = sinB*sinA;
            rotation_matrix[1][0] = cosC*sinA + cosB*cosA*sinC;
            rotation_matrix[1][1] = -sinC*sinA + cosB*cosA*cosC;
            rotation_matrix[1][2] = -sinB*cosA;
            rotation_matrix[2][0] = sinC*sinB;
            rotation_matrix[2][1] = cosC*sinB;
            rotation_matrix[2][2] = cosB;

            return;
      }
     
      static inline void GetEulerAngles(const double rotation_matrix[3][3], array_1d<double, 3 > & EulerAngles){
          
            if (rotation_matrix[2][2] < 1.0) {
                if (rotation_matrix[2][2] > -1.0) {
                    EulerAngles[0] = atan2(rotation_matrix[0][2], -rotation_matrix[1][2]);
                    EulerAngles[1] = acos(rotation_matrix[2][2]);
                    EulerAngles[2] = atan2(rotation_matrix[2][0], rotation_matrix[2][1]);
                } else // r22 = -1
                {
                    // Not a unique solution: thetaZ1 - thetaZ0 = atan2(-r01,r00)
                    EulerAngles[0] = -atan2(-rotation_matrix[0][1], rotation_matrix[0][0]);
                    EulerAngles[1] = KRATOS_M_PI;
                    EulerAngles[2] = 0;
                }
            } else // r22 = +1
            {
                // Not a unique solution: thetaZ1 + thetaZ0 = atan2(-r01,r00)
                EulerAngles[0] = atan2(-rotation_matrix[0][1], rotation_matrix[0][0]);
                EulerAngles[1] = 0;
                EulerAngles[2] = 0;
            }
            
            return;
      }
     
    
     ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     ///////////////****************************TRIANGLE - SPHERE INTERSECTION AREA CALCULATION**************************///////////
     ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
     static inline void AreaAndCentroidCircularSector(double C[3], double Radius, double P1[3], double P2[3], double Normal[3], double& Area, double CoMSC[3])
    {
        
        double a[3]           = {0.0};
        double c[3]           = {0.0};
        double bisection[3]   = {0.0};
        double norm_a         = 0.0;
        
        for (unsigned int index = 0;index<3;index++){        
        
            a[index] = P1[index]-C[index];
            c[index] = P2[index]-P1[index];
            
        }
        
        CrossProduct(Normal,c,bisection);
        normalize(bisection);       
        double dot_product = DotProduct(bisection,a);
        
        if(dot_product<0.0){
        
            for (unsigned int index = 0;index<3;index++){        
        
                bisection[index] = -bisection[index];
            }
            dot_product = -dot_product;
            
        }
        
        module(a,norm_a);
        
        double cos_alpha = dot_product/norm_a;
        double alpha = acos(cos_alpha);
        double sin_alpha = sin(alpha);
        
        Area = Radius*Radius*alpha;
        double dist = 0.666666666*(Radius*sin_alpha/alpha);
        for (unsigned int index = 0;index<3;index++){        
            CoMSC[index] = C[index]+dist*bisection[index];
        }
            
    }//AreaCircularSector

    static inline void AlternativeAreaCircularSegment(double Radius, double tol_Radius, double V0V1[3], double V0CC[3], double Normal[3], double& AreaSC, bool& flag)
    {
    
        double normal_outwards[3] = {0.0};
        flag = false;
        AreaSC = 0.0;
        
        CrossProduct(V0V1,Normal,normal_outwards);    
        normalize(normal_outwards);

        double dist   = DotProduct(normal_outwards,V0CC);
        double delta_circle  = Radius + dist; //dist can be positive or negative, depending on the side where the circle is
        
        if ( (delta_circle > tol_Radius) && ( delta_circle - 2*Radius < -tol_Radius ) ){//check for intersection
        
            
            flag = true;
            double b        = sqrt(delta_circle*(2*Radius-delta_circle));
            AreaSC   = 2*Radius*Radius*atan(delta_circle/b)-b*(Radius-delta_circle);
            
        }   
     
    }//AreaAndCentroidCircularSector1
    
    static inline void AreaAndCentroidCircularSegment(double Centre[3], double Radius, double tol_Radius, double V0[3], double V1[3], double Normal[3], double& AreaSegC, double CoMSegC[3], bool& flag)
    {
        double V0V1[3]            = {0.0};  double V0CC[3]      = {0.0}; 
        double a[3]               = {0.0};
        double normal_outwards[3] = {0.0};
        double Radius_SQ = 0.0;
        double distance_V0V1      = 0.0;
        double dist_CoM           = 0.0;
        flag = false;
        
        AreaSegC = 0.0;
        
        for (unsigned int index = 0; index<3; index++){        

            V0V1[index]     = V1[index] - V0[index];
            V0CC[index]     = Centre[index] - V0[index];
            
        }
       
        GeometryFunctions::CrossProduct(V0V1,Normal,normal_outwards); 
        GeometryFunctions::normalize(V0V1,distance_V0V1);
        
        double distV0 =  GeometryFunctions::DotProduct(V0CC,V0V1);
                   
        if( (distV0 > 0.0) && (distV0 < distance_V0V1)){
            
            GeometryFunctions::normalize(normal_outwards);
            double dist_normal   = GeometryFunctions::DotProduct(normal_outwards,V0CC);
            double delta_circle  = Radius + dist_normal; //dist can be positive or negative, depending on the side where the circle is
            
            if ( (delta_circle > tol_Radius) && ( delta_circle - 2*Radius < -tol_Radius ) ){//check for intersection
               
                Radius_SQ = Radius*Radius;
                double semi_dist = sqrt(Radius_SQ - dist_normal*dist_normal);
                flag = true;
                
                for (unsigned int index = 0;index<3;index++){        

                    a[index] = V0[index] + (distV0 - semi_dist)*V0V1[index] - Centre[index]; //Vector from Centre to first intersection point 
                    
                }
                
                double cos_alpha = GeometryFunctions::DotProduct(a,normal_outwards)/(GeometryFunctions::module(a)*GeometryFunctions::module(normal_outwards));
                double alpha = acos(cos_alpha);
                double sin_alpha = sin(alpha);
                
                AreaSegC = Radius_SQ*(alpha-sin_alpha*cos_alpha);

                if(abs(sin_alpha)<tol_Radius){dist_CoM=0.0;}
                else{ dist_CoM = 0.666666666*(Radius*sin_alpha*sin_alpha*sin_alpha/(alpha-sin_alpha*cos_alpha));}
                
                for (unsigned int index = 0;index<3;index++){        
                        CoMSegC[index] = Centre[index] + dist_CoM*normal_outwards[index];
                }
            } //if normal dist is okay
           
        } // if longitudinal dist is okay
                       
    }//AreaAndCentroidCircularSegment
    
    static inline void  AreaAndCentroidTriangle(double Coord1[3],double Coord2[3],double Coord3[3],double &area,double CoMTri[3])
    {
      
       TriAngleArea(Coord1,Coord2,Coord3,area);
       
       for (unsigned int index =0; index<3; index++)
       {
           CoMTri[index] = 0.3333333333*(Coord1[index]+Coord2[index]+Coord3[index]);
           
       }
       
    }//AreaAndCentroidTriangle
       
    
    
    } //namespace GeometryFunctions
    
} //namespace Kratos

#endif	/* _GEOMETRYFUNCTIONS_H */

