//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: Nelspon $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "spheric_particle.h"


namespace Kratos
{
      
      SphericParticle::SphericParticle( IndexType NewId, GeometryType::Pointer pGeometry) : DiscreteElement(NewId, pGeometry) {}
      
      SphericParticle::SphericParticle( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
      : DiscreteElement(NewId, pGeometry, pProperties)  
      {}
       
      SphericParticle::SphericParticle(IndexType NewId, NodesArrayType const& ThisNodes)
      : DiscreteElement(NewId, ThisNodes) 
      {}
           
      Element::Pointer SphericParticle::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const 
      {
         return DiscreteElement::Pointer(new SphericParticle(NewId, GetGeometry().Create( ThisNodes ), pProperties) );
      } 
      

      /// Destructor.
      SphericParticle::~SphericParticle(){}
      

      void SphericParticle::Initialize(){

          SetInitialContacts();








      }
      void SphericParticle::CalculateRightHandSide(VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo){}
      void SphericParticle::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo){}
      void SphericParticle::MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
      {
          
          double radius = GetGeometry()(0)->GetSolutionStepValue(RADIUS);
          double volume =   1.333333333333333*M_PI*radius*radius*radius;
          double density = GetGeometry()(0)->GetSolutionStepValue(PARTICLE_DENSITY);
          rMassMatrix.resize(1,1);
          rMassMatrix(0,0) = volume*density;
          
          KRATOS_WATCH("calul del massmatrix")
          KRATOS_WATCH(rMassMatrix)
          
          
      }
      


       void SphericParticle::SetInitialContacts()
        {
           /*
           Node<3>& my_node = GetGeometry()[0];
           //typedef WeakPointerVector<SphericParticle > ParticleWeakVectorType;  // per exemple per als initials neighbours teniem definida akest tipus pero.. on ho hauria de fer el el typedef, en el h del objecte oi=??
           ParticleWeakVectorType& initial_neighbours = my_node.GetValue(INITIAL_NEIGHBOUR_ELEMENTS);
           mInitialNeighbours = ParticleWeakVectorType();
           mInitialDelta = std::vector<double>();   ///M: com els defineixo?
           mContactFailureId = std::vector<int>();

        NeighbourElementsVectorType& neighbours_elements = my_node.GetValue(NEIGHBOUR_ELEMENTS);


        for(ParticleWeakIteratorType_ptr ineighbour = neighbours_elements.ptr_begin();
            ineighbour != neighbours_elements.ptr_end(); ineighbour++){

           initial_neighbours.push_back(*ineighbour);

           array_1d<double,3> other_to_me_vect = this->GetPosition() - ((*ineighbour).lock())->GetPosition();
           double distance                     = sqrt(other_to_me_vect[0] * other_to_me_vect[0] +
                                                 other_to_me_vect[1] * other_to_me_vect[1] +
                                                 other_to_me_vect[2] * other_to_me_vect[2]);

            double radius_sum                   = this->GetRadius() + ((*ineighbour).lock())->GetRadius();
            double initial_delta                = radius_sum - distance;


            mInitialDelta.push_back(initial_delta);
            mContactFailureId.push_back(1); //by default "generally detached" = 1

            } //end for: ParticleWeakIteratorType ineighbour
            * */
        } // SetInitialContacts

    void SphericParticle::AddContinuumContacts()
    {
        /*
        int continuum_group = mContinuumGroup;
         if (continuum_group == 0){
             return; //The particle neighbours'mContactFailureId are keep as 1 if the particle is not simulating the continuum.

         } //is it still continuum-simulating particle?

        unsigned int index = 0;
        for (ParticleWeakIteratorType iInitialNeighbours = mInitialNeighbours.begin(); iInitialNeighbours != mInitialNeighbours.end(); ++iInitialNeighbours)
        {

            int other_continuum_group = iInitialNeighbours->GetContinuumGroup();

            /* this loop will set only the 0 (contunuum simulating case) to the initial neighbours. The force calculator will change this
             * values depending of the type of failure as it is describre here:
             *
             *   mContactFailureId values:
             *      0 := Still a continuum simulating contact
             *      1 := General detachment (no initial continuum case: non continuum simulating particles or particles from diferent continuum group.)
             *      2 := Partially detached
             *      3 := tensile case
             *      4 := shear case
             *      5 :=von Misses.....M: define new cases...
             

            if(other_continuum_group != continuum_group ){
                mContactFailureId[index]=1;
              }
            else {     //same group
             mContactFailureId[index]=0;
             }
             index++;
          
        }//end for: ParticleWeakIteratorType iInitialNeighbours
      */
    }// AddContinuumContacts












      void SphericParticle::SphericParticle::DampMatrix(MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo){}
      
      
      
      void SphericParticle::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo){}
      void SphericParticle::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo){}
      void SphericParticle::FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo){}
      void SphericParticle::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo){}
      void SphericParticle::Calculate(const Variable<array_1d<double, 3 > >& rVariable, array_1d<double, 3 > & Output, const ProcessInfo& rCurrentProcessInfo){}     
      void SphericParticle::Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo){}
      void SphericParticle::Calculate(const Variable<Matrix >& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo){}
  
  
}  // namespace Kratos.



