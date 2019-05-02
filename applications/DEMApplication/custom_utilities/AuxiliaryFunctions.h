/*
 * File:   AuxiliaryFunctions.h
 * Author: msantasusana
 *
 * Created on 21 de mayo de 2012, 19:40
 */

#ifndef _DEM_AUXILIARY_FUNCTIONS_H
#define	_DEM_AUXILIARY_FUNCTIONS_H

#include <cmath>
#include "containers/array_1d.h"
#include "../../../kratos/includes/serializer.h"
#include "includes/model_part.h"
#include "../DEM_application_variables.h"

namespace Kratos {

    namespace AuxiliaryFunctions {

	static inline void CalculateAlphaFactor3D(const int n_neighbours, const double external_sphere_area, const double total_equiv_area , double& alpha) {

	    double external_polyhedron_area = 0.0;

	    switch (n_neighbours) {

                case 4:
                    external_polyhedron_area = 3.30797*external_sphere_area;
		    break;

                case 5:
		    external_polyhedron_area = 2.60892*external_sphere_area;
		    break;

                case 6:
		    external_polyhedron_area = 1.90986*external_sphere_area;
		    break;

                case 7:
		    external_polyhedron_area = 1.78192*external_sphere_area;
		    break;

                case 8:
                    external_polyhedron_area = 1.65399*external_sphere_area;
		    break;

                case 9:
		    external_polyhedron_area = 1.57175*external_sphere_area;
		    break;

                case 10:
		    external_polyhedron_area = 1.48951*external_sphere_area;
		    break;

                case 11:
		    external_polyhedron_area = 1.40727*external_sphere_area;
		    break;

                case 12:
		    external_polyhedron_area = 1.32503*external_sphere_area;
		    break;

                case 13:
		    external_polyhedron_area = 1.31023*external_sphere_area;
		    break;

                case 14:
		    external_polyhedron_area = 1.29542*external_sphere_area;
		    break;

                case 15:
		    external_polyhedron_area = 1.28061*external_sphere_area;
		    break;

                case 16:
		    external_polyhedron_area = 1.26580*external_sphere_area;
		    break;

                case 17:
		    external_polyhedron_area = 1.25099*external_sphere_area;
		    break;

                case 18:
		    external_polyhedron_area = 1.23618*external_sphere_area;
		    break;

                case 19:
		    external_polyhedron_area = 1.22138*external_sphere_area;
		    break;

                case 20:
		    external_polyhedron_area = 1.20657*external_sphere_area;
		    break;

                default:
		    external_polyhedron_area = 1.15000*external_sphere_area;
		    break;

            }//switch (n_neighbours)

            alpha = external_polyhedron_area/total_equiv_area;

	}//CalculateAlphaFactor

	static inline void CalculateAlphaFactor2D(const int n_neighbours, const double external_circle_perimeter, const double total_equiv_perimeter , double& alpha) {

            double external_polygon_perimeter = 0.0;

            switch (n_neighbours) {

                case 3:
                    external_polygon_perimeter = 1.65399*external_circle_perimeter;
                    break;

                case 4:
                    external_polygon_perimeter = 1.27324*external_circle_perimeter;
                    break;

                case 5:
                    external_polygon_perimeter = 1.15633*external_circle_perimeter;
                    break;

                case 6:
                    external_polygon_perimeter = 1.10266*external_circle_perimeter;
                    break;

                case 7:
                    external_polygon_perimeter = 1.07303*external_circle_perimeter;
                    break;

                case 8:
                    external_polygon_perimeter = 1.05479*external_circle_perimeter;
                    break;

                case 9:
                    external_polygon_perimeter = 1.04270*external_circle_perimeter;
                    break;

                case 10:
                    external_polygon_perimeter = 1.03425*external_circle_perimeter;
                    break;

                case 11:
                    external_polygon_perimeter = 1.02811*external_circle_perimeter;
                    break;

                case 12:
                    external_polygon_perimeter = 1.02349*external_circle_perimeter;
                    break;

                case 13:
                    external_polygon_perimeter = 1.01993*external_circle_perimeter;
                    break;

                case 14:
                    external_polygon_perimeter = 1.01713*external_circle_perimeter;
                    break;

                default:
                    external_polygon_perimeter = 1.0*external_circle_perimeter;

                break;

            }//switch (n_neighbours)

            alpha = external_polygon_perimeter/total_equiv_perimeter;

        }//CalculateAlphaFactor

	static inline void SwitchCase(const int case_opt, bool& delta_OPTION, bool& continuum_simulation_OPTION) {

	    switch (case_opt) {

                case 0:
                    delta_OPTION = false;
	            continuum_simulation_OPTION = false;
	            break;

	        case 1:
                    delta_OPTION = true;
		    continuum_simulation_OPTION = false;
	            break;

	        case 2:
		    delta_OPTION = true;
		    continuum_simulation_OPTION = true;
		    break;

		case 3:
		    delta_OPTION = false;
		    continuum_simulation_OPTION = true;
	            break;

		default:
                    delta_OPTION = false;
	            continuum_simulation_OPTION = false;
	    }
	} //SwitchCase

	inline array_1d<double,3> LinearTimeIncreasingFunction(const array_1d<double,3>& external_total_applied_force, const double current_time, const double end_time)
	{
            array_1d<double,3> externally_applied_force_now = external_total_applied_force*current_time/end_time;
            return externally_applied_force_now;

        }// LinearTimeIncreasingFunction

        static inline void ComputeReactionOnTopAndBottomSpheres(ModelPart& r_model_part) {

            typedef ModelPart::ElementsContainerType ElementsArrayType;
            ElementsArrayType& pElements = r_model_part.GetCommunicator().LocalMesh().Elements();

            //double Y_coord    = 0.0;
            double Z_coord    = 0.0;
            double RX_bottom  = 0.0;
            double RY_bottom  = 0.0;
            double RZ_bottom  = 0.0;
            double RX_top     = 0.0;
            double RY_top     = 0.0;
            double RZ_top     = 0.0;
            double RX_average = 0.0;
            double RY_average = 0.0;
            double RZ_average = 0.0;
            double time       = 0.0;

            ElementsArrayType::iterator elem_iterator_begin = pElements.ptr_begin();
            ElementsArrayType::iterator elem_iterator_end = pElements.ptr_end();

            for (ElementsArrayType::iterator elem_iterator = elem_iterator_begin; elem_iterator != elem_iterator_end; ++elem_iterator) {
                //Node<3>& node = elem_iterator->GetGeometry()[0];
                array_1d<double, 3>& reaction_force = elem_iterator->GetGeometry()[0].FastGetSolutionStepValue(FORCE_REACTION);
                //Y_coord = elem_iterator->GetGeometry()[0].Coordinates()[1];
                Z_coord = elem_iterator->GetGeometry()[0].Coordinates()[2];
                //if (Y_coord < 0.15) { /////////////////////////////////// Take this correctly!!
                //if ((Z_coord < 0.3) || (Y_coord > 0.3)) {
                if (Z_coord < 0.25) {
                    RX_bottom += reaction_force[0];
                    RY_bottom += reaction_force[1];
                    RZ_bottom += reaction_force[2];
                } else {
                    RX_top += reaction_force[0];
                    RY_top += reaction_force[1];
                    RZ_top += reaction_force[2];
                }
            } //loop over particles

            time = r_model_part.GetProcessInfo()[TIME];
            RX_average = 0.5 * (RX_bottom - RX_top);
            RY_average = 0.5 * (RY_bottom - RY_top);
            RZ_average = 0.5 * (RZ_bottom - RZ_top);

            std::ofstream outputfileX("reaction_forces_X.txt", std::ios_base::out | std::ios_base::app);
            std::ofstream outputfileZ("reaction_forces_Z.txt", std::ios_base::out | std::ios_base::app);
            std::ofstream outputfileYZ("reaction_forces_YZ.txt", std::ios_base::out | std::ios_base::app);

            outputfileX  << time << " " << RX_bottom << " " << RX_top << " " << RX_average << "\n";
            outputfileZ  << time << " " << RZ_bottom << " " << RZ_top << " " << RZ_average << "\n";
            outputfileYZ << time << " " << RY_bottom << " " << RY_top << " " << RY_average << " " << RZ_bottom << " " << RZ_top << " " << RZ_average << "\n";

            outputfileX.close();
            outputfileZ.close();
            outputfileYZ.close();
        }

        static inline Vector EigenValuesDirectMethod(const Matrix& A) {
            // Given a real symmetric 3x3 matrix A, compute the eigenvalues
            const int dim= A.size1();
            Vector Result=ZeroVector(dim);

            const double p1 = A(0,1)*A(0,1) + A(0,2)*A(0,2) + A(1,2)*A(1,2);
            if (p1 == 0) {//A is diagonal.
                Result[0] = A(0,0);
                Result[1] = A(1,1);
                Result[2] = A(2,2);
                return Result;
            }

            const double q = 0.333333333333333333333333 * (A(0,0) + A(1,1) + A(2,2));
            const double p2 = (A(0,0) - q) * (A(0,0) - q) + (A(1,1) - q) * (A(1,1) - q) + (A(2,2) - q) * (A(2,2) - q) + 2.0 * p1;
            const double p = sqrt(0.16666666666666666666666667 * p2);

            Matrix B(3,3);
            const double inv_p = 1.0/p;

            // B = (1 / p) * (A - q * I)  where  I is the identity matrix
            B(0,0) = inv_p * (A(0,0) - q);
            B(1,1) = inv_p * (A(1,1) - q);
            B(2,2) = inv_p * (A(2,2) - q);
            B(0,1) = inv_p * A(0,1);
            B(1,0) = inv_p * A(1,0);
            B(0,2) = inv_p * A(0,2);
            B(2,0) = inv_p * A(2,0);
            B(1,2) = inv_p * A(1,2);
            B(2,1) = inv_p * A(2,1);

            //r = det(B) / 2
            double r = 0.5 * ( B(0,0)*B(1,1)*B(2,2) + B(0,1)*B(1,2)*B(2,0) + B(1,0)*B(2,1)*B(0,2) - B(2,0)*B(1,1)*B(0,2) - B(1,0)*B(0,1)*B(2,2) - B(0,0)*B(2,1)*B(1,2) );

            // In exact arithmetic for a symmetric matrix  -1 <= r <= 1
            // but computation error can leave it slightly outside this range.
            double phi = 0.0;
            if (r <= -1) { phi = 0.333333333333333333333333 * Globals::Pi; }
            else if (r >= 1) { phi = 0.0; }
            else { phi = 0.333333333333333333333333 * acos(r);}

            // the eigenvalues satisfy eig3 <= eig2 <= eig1 MAC: WATCH OUT!! This is not true, apparently!
            Result[0] = q + 2.0 * p * std::cos(phi);
            Result[2] = q + 2.0 * p * std::cos(phi + (0.6666666666666666666666*Globals::Pi));
            Result[1] = 3.0 * q - Result[0] - Result[2];     //% since trace(A) = eig1 + eig2 + eig3

            return Result;
        }

    }//namespace AuxiliaryFunctions

    namespace DemDebugFunctions {
        static inline void CheckIfNan(const array_1d<double,3>& vector, const std::string& sentence) {
            #ifdef KRATOS_DEBUG
            if(std::isnan(vector[0]) || std::isnan(vector[1]) || std::isnan(vector[2])){
                KRATOS_ERROR<<sentence;
            }
            #endif
        }

        static inline void CheckIfNan(const double vector[3], const std::string& sentence) {
            #ifdef KRATOS_DEBUG
            if(std::isnan(vector[0]) || std::isnan(vector[1]) || std::isnan(vector[2])){
                KRATOS_ERROR<<sentence;
            }
            #endif
        }

        static inline void CheckIfNan(const double& value, const std::string& sentence) {
            #ifdef KRATOS_DEBUG
            if(std::isnan(value)){
                KRATOS_ERROR<<sentence;
            }
            #endif
        }

    }

}//namespace Kratos

#endif	/* _KRATOSAUXILIARYFUNCTIONS_H */

