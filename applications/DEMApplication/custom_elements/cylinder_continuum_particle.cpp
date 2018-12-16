//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//

#include "cylinder_continuum_particle.h"

namespace Kratos
{
     // using namespace GeometryFunctions;

      CylinderContinuumParticle::CylinderContinuumParticle()
      : SphericContinuumParticle(){/*mInitializedVariablesFlag = 0;*/}

      CylinderContinuumParticle::CylinderContinuumParticle(IndexType NewId, GeometryType::Pointer pGeometry)
      : SphericContinuumParticle(NewId, pGeometry){/*mInitializedVariablesFlag = 0;*/}

      CylinderContinuumParticle::CylinderContinuumParticle(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
      : SphericContinuumParticle(NewId, pGeometry, pProperties){/*mInitializedVariablesFlag = 0;*/}

      CylinderContinuumParticle::CylinderContinuumParticle(IndexType NewId, NodesArrayType const& ThisNodes)
      : SphericContinuumParticle(NewId, ThisNodes){/*mInitializedVariablesFlag = 0;*/}

      Element::Pointer CylinderContinuumParticle::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
      {
           return SphericContinuumParticle::Pointer(new CylinderContinuumParticle(NewId, GetGeometry().Create(ThisNodes), pProperties));

      }

      /// Destructor.
      CylinderContinuumParticle::~CylinderContinuumParticle(){} 
      
        void CylinderContinuumParticle::ContactAreaWeighting() //MISMI 10: POOYAN this could be done by calculating on the bars. not looking at the neighbous of my neighbours.
        { 

            double alpha = 1.0;
            double circle_perimeter = 2*Globals::Pi * GetRadius();       
            double total_equiv_perimeter = 0.0;
            unsigned int continuous_initial_neighbours_size = mContinuumInitialNeighborsSize;
            Vector& cont_ini_neigh_area = GetValue(NEIGHBOURS_CONTACT_AREAS);
        
            for (unsigned int i = 0; i < continuous_initial_neighbours_size; i++) {
                SphericParticle* ini_cont_neighbour_iterator = mNeighbourElements[i];
                double other_radius     = ini_cont_neighbour_iterator->GetInteractionRadius();
                double area = mContinuumConstitutiveLawArray[i]->CalculateContactArea(GetRadius(), other_radius, cont_ini_neigh_area); //This call fills the vector of areas only if the Constitutive Law wants.         
                total_equiv_perimeter += area;
            } //for every neighbour
      
            if (continuous_initial_neighbours_size >= 4) {
            
                if (!IsSkin()) {
                    AuxiliaryFunctions::CalculateAlphaFactor2D(continuous_initial_neighbours_size, circle_perimeter, total_equiv_perimeter, alpha); 
                    for (unsigned int i = 0; i < cont_ini_neigh_area.size(); i++) {
                        cont_ini_neigh_area[i] = alpha*cont_ini_neigh_area[i];                      
                    } //for every neighbour
                }
                else { //skin sphere                              
                    for (unsigned int i = 0; i < cont_ini_neigh_area.size(); i++) {
                        alpha            = 1.30*(1.10266)*(circle_perimeter/total_equiv_perimeter)*((double(continuous_initial_neighbours_size))/6); // 6 is mean coordination number.
                        cont_ini_neigh_area[i] = alpha*cont_ini_neigh_area[i];
                    }     //loop on cont neighs
                }           //skin particles.
            }               //if 3 neighbours or more.
        }                 //Contact Area Weighting
               
      double CylinderContinuumParticle::CalculateVolume() {
          return Globals::Pi * GetRadius() * GetRadius();
      }
      
      double CylinderContinuumParticle::CalculateMomentOfInertia() {
          return 0.5 * GetMass() * GetRadius() * GetRadius(); 
      }

      void CylinderContinuumParticle::AddContributionToRepresentativeVolume(const double distance,
                                                                         const double radius_sum,
                                                                         const double contact_area) {

        KRATOS_TRY

        double gap = distance - radius_sum;
        double real_distance = GetInteractionRadius() + 0.5 * gap;
        double& rRepresentative_Volume = this->GetGeometry()[0].FastGetSolutionStepValue(REPRESENTATIVE_VOLUME);
        rRepresentative_Volume += 0.5 * (real_distance * contact_area);

        KRATOS_CATCH("")
    }
//      void CylinderContinuumParticle::AddNeighbourContributionToStressTensor(double GlobalElasticContactForce[3],
//                                                                            array_1d<double,3> &other_to_me_vect,
//                                                                            const double &distance,
//                                                                            const double &radius_sum,
//                                                                            const double &calculation_area,
//                                                                            ParticleWeakIteratorType neighbour_iterator,
//                                                                            ProcessInfo& r_process_info,
//                                                                            double &rRepresentative_Volume){
//        double gap  = distance - radius_sum;
//        array_1d<double,3> normal_vector_on_contact =  -1 * other_to_me_vect; //outwards
//        double value = 0.0;
//        GeometryFunctions::normalize(normal_vector_on_contact,value); // Normalize to unitary module
//        array_1d<double,3> x_centroid      = (GetRadius() + 0.5*gap) * normal_vector_on_contact;
//        array_1d<double,3> surface_baricenter = x_centroid;
//        double result_product = GeometryFunctions::DotProduct(surface_baricenter,normal_vector_on_contact);
//        //Aproximation with error: surface_baricenter should be the baricenter of each surface, which can no be calculated because the surfaces are imaginary.
//        rRepresentative_Volume = rRepresentative_Volume + 0.5 * (result_product * calculation_area);
        
//        for (int i=0; i<3; i++){
//            for (int j=0; j<3; j++){
//                (*mStressTensor)(i,j) += (x_centroid[j]) * GlobalElasticContactForce[i]; //ref: Katalin Bagi 1995 Mean stress tensor
//            }
//        }
//      } //AddNeighbourContributionToStressTensor




      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************
//       void CylinderContinuumParticle::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info){}
//       void CylinderContinuumParticle::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output, const ProcessInfo& r_process_info){}
//       void CylinderContinuumParticle::Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& r_process_info){}
//       void CylinderContinuumParticle::Calculate(const Variable<Matrix >& rVariable, Matrix& Output, const ProcessInfo& r_process_info){}

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

}  // namespace Kratos.

