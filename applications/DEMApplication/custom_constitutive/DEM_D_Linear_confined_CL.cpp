
#include "DEM_D_Linear_confined_CL.h"
#include "custom_elements/spheric_particle.h"

namespace Kratos {

    DEMDiscontinuumConstitutiveLaw::Pointer DEM_D_Linear_confined::Clone() const {
        DEMDiscontinuumConstitutiveLaw::Pointer p_clone(new DEM_D_Linear_confined(*this));
        return p_clone;
    }

    void DEM_D_Linear_confined::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) {
        KRATOS_INFO("DEM") << "Assigning DEM_D_Linear_confined to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_DISCONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
    }

    double DEM_D_Linear_confined::CalculateNormalForce(SphericParticle* const element1, SphericParticle* const element2,
            const double indentation, double LocalCoordSystem[3][3]) {

        //const double my_young        = element1->GetYoung();
        //const double other_young     = element2->GetYoung();

        const double my_poisson      = element1->GetPoisson();
        const double other_poisson   = element2->GetPoisson();


            double equiv_poisson;
            if ((my_poisson + other_poisson) != 0.0) { equiv_poisson = 2.0 * my_poisson * other_poisson / (my_poisson + other_poisson); }
            else { equiv_poisson = 0.0; }


            //Get equivalent Radius
        const double my_radius       = element1->GetRadius();
        const double other_radius    = element2->GetRadius();
        const double radius_sum      = my_radius + other_radius;
        const double radius_sum_inv  = 1.0 / radius_sum;
        const double equiv_radius    = my_radius * other_radius * radius_sum_inv;


           // double calculated_radius  = sqrt(equiv_radius*indentation); // sosss

          // double  calculation_area = Globals::Pi * calculated_radius * calculated_radius;

        double  calculation_area = Globals::Pi * equiv_radius * indentation ;


       //double rmin = my_radius;
       //if (other_radius < my_radius) rmin = other_radius;
       //double  calculation_area = 3.14 * rmin*rmin;


        double normal_force = DEM_D_Linear_viscous_Coulomb::CalculateNormalForce(indentation);
           //double normal_force = 0.666666666666666666667 * mKn * indentation;
           //double normal_force = 0 * mKn * indentation;
                          // KRATOS_INFO("DEM") << "indentation  " << indentation << std::endl;

        double force[3];
        BoundedMatrix<double, 3, 3> average_stress_tensor = ZeroMatrix(3,3);

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                average_stress_tensor(i,j) = 0.5 * ((*(element1->mSymmStressTensor))(i,j) + (*(element2->mSymmStressTensor))(i,j));
                //average_stress_tensor(i,j) = 0.5 * ((*mSymmStressTensor)(i,j) + (*(element2->mSymmStressTensor))(i,j));


            //KRATOS_INFO("DEM") << average_stress_tensor(i,j)<< "\n ";
            //KRATOS_INFO("DEM") << average_stress_tensor(i,j) << " " ;

            }

                       //  KRATOS_INFO("DEM") << std::endl;


        }

        for (int i = 0; i < 3; i++) {

            force[i] = (average_stress_tensor)(i,0) * LocalCoordSystem[0][0] +
                       (average_stress_tensor)(i,1) * LocalCoordSystem[0][1] +
                       (average_stress_tensor)(i,2) * LocalCoordSystem[0][2]; // StressTensor*unitaryNormal0
        }

        double sigma_x = force[0] * LocalCoordSystem[0][0] +
                         force[1] * LocalCoordSystem[0][1] +
                         force[2] * LocalCoordSystem[0][2]; // projection to normal to obtain value of the normal stress

        for (int i = 0; i < 3; i++) {

            force[i] = (average_stress_tensor)(i,0) * LocalCoordSystem[1][0] +
                       (average_stress_tensor)(i,1) * LocalCoordSystem[1][1] +
                       (average_stress_tensor)(i,2) * LocalCoordSystem[1][2]; // StressTensor*unitaryNormal1
        }

        double sigma_y = force[0] * LocalCoordSystem[1][0] +
                         force[1] * LocalCoordSystem[1][1] +
                         force[2] * LocalCoordSystem[1][2]; // projection to normal to obtain value of the normal stress

        double poisson_force = calculation_area * equiv_poisson * (sigma_x + sigma_y);

        normal_force -= poisson_force;



      return normal_force;
    }


} // namespace Kratos
