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
        int output;
        if(a < 0.0) output = -1;
        else if (a > 0.0) output = 1;
        else output = 0;
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
    static inline void ComputeContactLocalCoordSystem(array_1d<double,3> NormalDirection, const double& distance, double LocalCoordSystem[3][3])  //inline: modifies the LocalCoordSystem as it were a reference
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

     static inline double DistancePointToPlane(array_1d<double,3> CoordInPlane, double PlaneUnitNormalVector[3], double TestCoord[3])
     {
         double Vector1[3] = {0.0};

        for (unsigned int i = 0; i<3; i++)
        {
            Vector1[i] = TestCoord[i]- CoordInPlane[i];
        }

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

     static inline void CoordProjectionOnPlaneNew(double CoordOut[3], array_1d<double, 3> CoordIn, double LocalCoordSystem[3][3], double IntersectionCoord[3])
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

     static inline void Compute3DimElementFaceLocalSystem(array_1d <double,3> FaceCoord1, array_1d <double,3> FaceCoord2, array_1d <double,3> FaceCoord3, double ParticleCoord[3],
                                                          double LocalCoordSystem[3][3], double& normal_flag)
    {
        //NOTE: this function is designed in a way that the normal always points the side where the centre of particle is found. Therefore should only be used in this way if the indentation is less than the radius value.
        //the function returns a flag with the same value as the dotproduct of the normal of the triangle and the normal pointing to the particle.

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

    }
    
    //MSIMSI this one is being used only for distributed... adapt it
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


   
    static inline  bool InsideOutside(array_1d<double, 3> Coord1, array_1d<double, 3> Coord2, double JudgeCoord[3], double normal_out[3], double& area){

        double cp1[3]  = {0.0};
        double b_a[3]  = {0.0};
        double p1_a[3] = {0.0};

        for (int i=0;i<3;i++){

            b_a[i]  = Coord2[i]-Coord1[i];
            p1_a[i] = JudgeCoord[i]-Coord1[i];
        }

        GeometryFunctions::CrossProduct(b_a, p1_a, cp1);

        if (GeometryFunctions::DotProduct(cp1, normal_out) >= 0)
        {
            area = sqrt(cp1[0] * cp1[0] + cp1[1] * cp1[1] + cp1[2] * cp1[2]) / 2.0;
            return true;
        }
        else return false;

    }//PosCrossProd


    static inline void WeightsCalculation(std::vector<double> Area, std::vector<double> &Weight)
    {
      unsigned int facet_size = Area.size();
      if(facet_size == 3)
      {
        double total_area = Area[0]+Area[1]+Area[2];
        for (unsigned int i = 0; i< 3; i++)
        {
            Weight[i] = Area[(i+1)%facet_size]/total_area;
        }
      }
      else if (facet_size == 4)
      {
        double total_discriminant = Area[0]*Area[1]+Area[1]*Area[2]+Area[2]*Area[3]+Area[3]*Area[0]; //(Zhong et al 1993)
        for (unsigned int i = 0; i< 4; i++)
        {
            Weight[i] = (Area[(i+1)%facet_size]*Area[(i+2)%facet_size])/total_discriminant;
        }
      }
      else{KRATOS_WATCH("WEIGHTS FOR N-SIZE POLYGONAL FE TO BE IMPLEMENTED")}
    }//WeightsCalculation


    static inline bool FacetCheck(std::vector< array_1d <double,3> >Coord, double Particle_Coord[3], double rad,
                                  double LocalCoordSystem[3][3], double &DistPToB, std::vector<double> &Weight, unsigned int &current_edge_index) 
    {

        bool If_Contact = false;

        double normal_flag = 0.0;

        int facet_size = Coord.size();

        Compute3DimElementFaceLocalSystem(Coord[0], Coord[1], Coord[2], Particle_Coord, LocalCoordSystem, normal_flag);

        DistPToB = DistancePointToPlane(Coord[0], LocalCoordSystem[2], Particle_Coord);

        if(DistPToB < rad )
        {

            double IntersectionCoord[3];
            double PlaneUnitNormalVector[3];
            PlaneUnitNormalVector[0] = LocalCoordSystem[2][0]*normal_flag;
            PlaneUnitNormalVector[1] = LocalCoordSystem[2][1]*normal_flag;
            PlaneUnitNormalVector[2] = LocalCoordSystem[2][2]*normal_flag;

            //Projection of C onto plane
            for (unsigned int j = 0; j<3; j++)
            {
              IntersectionCoord[j] = Particle_Coord[j] - DistPToB*LocalCoordSystem[2][j]; //LocalCoordSystem[2] always points the particle

            }

            std::vector<double> Area;
            Area.resize(facet_size);

            for (int i = 0; i<facet_size; i++)
            {
              double this_area = 0.0;
              if( InsideOutside(Coord[i], Coord[(i+1)%facet_size], IntersectionCoord, PlaneUnitNormalVector, this_area) == false )
              {

                current_edge_index = i;
                return false;

              }
              else
              {
                Area[i] = this_area; //the area adjacent to vertex[ID] is assigned as Area[ID] so further treatment shall be done for the Weight calculation
              }

            }//for every vertex

            If_Contact = true;
            WeightsCalculation(Area,Weight);

      }//if DistPToB < rad

      return If_Contact;

    }

    static inline bool EdgeCheck( array_1d <double,3> Coord1, array_1d <double,3> Coord2, double Particle_Coord[3], double Radius,
                                  double LocalCoordSystem[3][3], double &DistEdgeToPlane, double& eta) 
    {
      double IntersectionCoordEdge[3];
      double normal_unit_vector[3];
      double edge_unit_vector[3];
      double module_edge_vector = 0.0;
      double particle_vector[3];

      for (unsigned int j = 0; j<3; j++)
      {
          edge_unit_vector[j] = Coord2[j] - Coord1[j];
          particle_vector[j]  = Particle_Coord[j] - Coord1[j];
      }

      normalize( edge_unit_vector, module_edge_vector);
      double projection_on_edge = DotProduct(particle_vector,edge_unit_vector);

      for (unsigned int j = 0; j<3; j++)
      {

        IntersectionCoordEdge[j] = Coord1[j] + projection_on_edge*edge_unit_vector[j];
        normal_unit_vector[j]   = Particle_Coord[j] - IntersectionCoordEdge[j];
      }

      normalize( normal_unit_vector, DistEdgeToPlane);

      if ( DistEdgeToPlane < Radius)
      {

        eta = projection_on_edge/module_edge_vector;

        if( (eta>=0.0) && (eta<=1.0) )
        {

          double dummy_length = 0.0;
          double auxiliar_unit_vector[3];
          CrossProduct(normal_unit_vector,edge_unit_vector,auxiliar_unit_vector);
          normalize( auxiliar_unit_vector, dummy_length);
          for (unsigned int j = 0; j<3; j++)
          {
            LocalCoordSystem[0][j] = edge_unit_vector[j];
            LocalCoordSystem[1][j] = auxiliar_unit_vector[j];
            LocalCoordSystem[2][j] = normal_unit_vector[j];
          }
          return true;
        }

      } //if distance to edge < radius

      return false;

    }//EdgeCheck

    static inline bool VertexCheck( array_1d <double,3> Coord, double Particle_Coord[3], double Radius, double LocalCoordSystem[3][3], double &DistEdgeToPlane)
    {

      double DistEdgeToPlaneSq = 0.0;
      array_1d<double,3> normal_v;
      for(unsigned int j = 0; j<3; j++)
      {
        normal_v[j] = Particle_Coord[j] - Coord[j];
        DistEdgeToPlaneSq += normal_v[j]*normal_v[j];
      }
      if(DistEdgeToPlaneSq <= Radius*Radius)
      {
        DistEdgeToPlane = sqrt(DistEdgeToPlaneSq);
        ComputeContactLocalCoordSystem(normal_v, DistEdgeToPlane, LocalCoordSystem);
        return true;
      }

      return false;
    }//VertexCheck

    
  
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

                if(fabs(sin_alpha)<tol_Radius){dist_CoM=0.0;}
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
