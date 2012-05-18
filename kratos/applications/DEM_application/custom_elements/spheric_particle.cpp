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

        KRATOS_TRY

        mDimension = this->GetGeometry().WorkingSpaceDimension();

        double density      = GetGeometry()(0)->FastGetSolutionStepValue(PARTICLE_DENSITY);
        double radius       = GetGeometry()(0)->FastGetSolutionStepValue(RADIUS);
        double& mass        = GetGeometry()(0)->FastGetSolutionStepValue(NODAL_MASS);

        mContinuumGroup     = this->GetGeometry()[0].GetSolutionStepValue(PARTICLE_CONTINUUM);
        mFailureId          = !(mContinuumGroup);

        if(mDimension ==2)
        {
            mass     = M_PI * radius * radius * density;
            
            mRealMass = mass;

            mInertia = 0.25 * M_PI * radius * radius * radius  * radius ;

            mMomentOfInertia = 0.5 * radius * radius;


        }
        else
        {
            mass     = 4.0 / 3.0 * M_PI * radius * radius * radius * density;

            mRealMass = mass;

            mInertia = 0.25 * M_PI * radius * radius * radius  * radius ;

            mMomentOfInertia = 0.4 * radius * radius;
        }

    
        KRATOS_CATCH( "" )

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
        
      }


      void SphericParticle::SetInitialContacts(int case_opt) //vull ficar que sigui zero si no son veins cohesius.
      {

          // DEFINING THE REFERENCES TO THE MAIN PARAMETERS

          //Node<3>& my_node = GetGeometry()[0];
           //typedef WeakPointerVector<SphericContinuumParticle > ParticleWeakVectorType;  // per exemple per als initials neighbours teniem definida akest tipus pero.. on ho hauria de fer el el typedef, en el h del objecte oi=??


           //mContinuumGroup                                  = this->GetValue(PARTICLE_CONTINUUM);

         
           ParticleWeakVectorType& r_neighbours             = this->GetValue(NEIGHBOUR_ELEMENTS);
           ParticleWeakVectorType& r_initial_neighbours     = this->GetValue(INITIAL_NEIGHBOUR_ELEMENTS);
          
           //vector<double>& r_VectorContactFailureId            = this->GetValue(PARTICLE_CONTACT_FAILURE_ID); //M:temporarily double...
           //r_VectorContactFailureId.resize(r_neighbours.size());

           vector<double>& r_VectorContactInitialDelta       = this->GetValue(PARTICLE_CONTACT_INITIAL_DELTA);
           //r_VectorContactInitialDelta.resize(r_neighbours.size()); //the initial neighbours are the neighbours at the first step.

           unsigned int i=0;


           //SAVING THE INICIAL NEIGHBOURS, THE DELTAS AND THE FAILURE ID

           for(ParticleWeakIteratorType_ptr ineighbour = r_neighbours.ptr_begin();  //loop over the neighbours and store into a initial_neighbours vector.
            ineighbour != r_neighbours.ptr_end(); ineighbour++){

              if (this->Id() != ((*ineighbour).lock())->Id() ){

                   array_1d<double,3> other_to_me_vect = this->GetGeometry()(0)->Coordinates() - ((*ineighbour).lock())->GetGeometry()(0)->Coordinates();
                   double distance                     = sqrt(other_to_me_vect[0] * other_to_me_vect[0] +
                                                         other_to_me_vect[1] * other_to_me_vect[1] +
                                                         other_to_me_vect[2] * other_to_me_vect[2]);

                    double radius_sum                   = this->GetGeometry()(0)->GetSolutionStepValue(RADIUS) + ((*ineighbour).lock())->GetGeometry()(0)->GetSolutionStepValue(RADIUS);
                    double initial_delta                = radius_sum - distance;

                  //  KRATOS_WATCH(initial_delta)

                    int r_other_continuum_group = ((*ineighbour).lock())->GetGeometry()(0)->GetSolutionStepValue(PARTICLE_CONTINUUM);

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
                         */
 //KRATOS_WATCH(r_other_continuum_group)
   //      KRATOS_WATCH(mContinuumGroup)
                        if( (r_other_continuum_group == mContinuumGroup) || ( fabs(initial_delta)>1.0e-6 ) ) // R: TOLERANCIA PER DIR QUE TENEN UN IDENTACIÓ O SEPARACIÓ; HAURIA DANAR LLIGADA AMB LU DE LA BUSQUEDA. la busqueda es fa per intersecció, es a dir identació es detecta pero si estan separades un pel ja no es troba com a veí.
                        //THESE ARE THE CASES THAT NEED TO STORE THE INITIAL NEIGHBOURS
                        {


                            r_initial_neighbours.push_back(*ineighbour);
                            r_VectorContactInitialDelta[i]  =   initial_delta; //M:R:els guardo sempre, en forces i tal si és car accedir aquest valor podem demanar el flag pero no crec ke sigui car.

                            if (r_other_continuum_group == mContinuumGroup) {this->GetValue(PARTICLE_CONTACT_FAILURE_ID)[i]=0;  KRATOS_WATCH("SON DEL MATEIX GRUP")  }
                            else                                             this->GetValue(PARTICLE_CONTACT_FAILURE_ID)[i]=1; // MRMR:  //generally detached    //diferent group
                                                                                                                             //hi havia: r_VectorContactFailureId[i]=1; //generally detached    //diferent group
//KRATOS_WATCH(r_VectorContactFailureId[i])


                        } // FOR THE CASES THAT NEED STORING INITIAL NEIGHBOURS

                        else mFailureId=1;      //NO NEED to store the initial neighbour;
                        //M:also we can say that the particle is detached --> this can be a flag for the forces.

                        i++;
                    
               }//if I found myself.

            } //end for: ParticleWeakIteratorType ineighbour
      }//SET INITIAL CONTACTS.
 


        void SphericParticle::ComputeForcesGeneral(const ProcessInfo& rCurrentProcessInfo )

        {

            KRATOS_TRY
           //KRATOS_WATCH("HOOOLA")
        //    KRATOS_WATCH(" ")
         //   KRATOS_WATCH(" ")

            // ALIAS TO ELEMENT VALUES

          //  Node<3>& my_node = GetGeometry()[0];

            ParticleWeakVectorType& r_neighbours             = this->GetValue(NEIGHBOUR_ELEMENTS);

            //KRATOS_WATCH(r_neighbours.size())
            //KRATOS_WATCH(this->Id())

            //ParticleWeakVectorType& r_initial_neighbours     = my_node.GetValue(INITIAL_NEIGHBOUR_ELEMENTS);  M: no l'utilitzo

            vector<double>& r_VectorContactFailureId            = this->GetValue(PARTICLE_CONTACT_FAILURE_ID); //M:temporarily a vector of doubles... it should be vector of ints..
            vector<double>& r_VectorContactInitialDelta         = this->GetValue(PARTICLE_CONTACT_INITIAL_DELTA);


            // PROCESS INFO

            const array_1d<double,3>& gravity   = rCurrentProcessInfo[GRAVITY];

            double dt                           = rCurrentProcessInfo[DELTA_TIME];
            int damp_id                         = rCurrentProcessInfo[DAMP_TYPE];
            int type_id                         = rCurrentProcessInfo[FORCE_CALCULATION_TYPE];

            int case_OPTION                     = rCurrentProcessInfo[CASE_OPTION];
            bool delta_OPTION;
            bool continuum_simulation_OPTION;

                switch (case_OPTION) {
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
  
            // GETTING PARTICLE PROPERTIES

        //    int FailureId = mFailureId;
            int continuum_group     = mContinuumGroup;

            double LocalDampRatio   = this->GetGeometry()[0].GetSolutionStepValue(PARTICLE_LOCAL_DAMP_RATIO);
            double Tension          = this->GetGeometry()[0].GetSolutionStepValue(PARTICLE_TENSION);
            double Cohesion         = this->GetGeometry()[0].GetSolutionStepValue(PARTICLE_COHESION);
            double FriAngle         = this->GetGeometry()[0].GetSolutionStepValue(PARTICLE_FRICTION); //M: CFENG: es aixo el angle de friccio?
            double Friction         = tan( FriAngle / 180.0 * M_PI);

            double radius               = this->GetGeometry()[0].GetSolutionStepValue(RADIUS);
            double critic_damp_fraction = this->GetGeometry()[0].GetSolutionStepValue(PARTICLE_ZETA);
            double mass                 = mRealMass;

            double young                = this->GetGeometry()[0].GetSolutionStepValue(YOUNG_MODULUS);
            double poisson              = this->GetGeometry()[0].GetSolutionStepValue(POISSON_RATIO);

            array_1d<double,3>& force           = this->GetGeometry()[0].GetSolutionStepValue(FORCE);//total forces, we reset to 0. and we calculate again.
       //     KRATOS_WATCH(force)
            noalias(force)                      = ZeroVector(3);
        //    KRATOS_WATCH(force)
            array_1d<double,3> applied_force    = this->GetGeometry()[0].GetSolutionStepValue(APPLIED_FORCE);
        //    KRATOS_WATCH(applied_force)
            force  = mass*gravity + applied_force;
        //    KRATOS_WATCH(force)
            
            size_t iContactForce = 0;
            

            for(ParticleWeakIteratorType neighbour_iterator = r_neighbours.begin();
            neighbour_iterator != r_neighbours.end(); neighbour_iterator++)
            {
             // GETTING NEIGHBOUR PROPERTIES
//KRATOS_WATCH("HOOOLA")
                double other_radius                 = neighbour_iterator->GetGeometry()(0)->GetSolutionStepValue(RADIUS);
                double other_critic_damp_fraction   = neighbour_iterator->GetGeometry()(0)->GetSolutionStepValue(PARTICLE_ZETA);
                double equiv_visc_damp_ratio        = (critic_damp_fraction + other_critic_damp_fraction) / 2.0;   //M: is it correct to be a simple mean.
             // double other_mass                   = neighbour_iterator.mRealMass;
                double other_young                  = neighbour_iterator->GetGeometry()[0].GetSolutionStepValue(YOUNG_MODULUS);
                double other_poisson                = neighbour_iterator->GetGeometry()[0].GetSolutionStepValue(POISSON_RATIO);
                double other_tension                = neighbour_iterator->GetGeometry()[0].GetSolutionStepValue(PARTICLE_TENSION);
                double other_cohesion               = neighbour_iterator->GetGeometry()[0].GetSolutionStepValue(PARTICLE_COHESION);
   
             // double other_FriAngle               = neighbour_iterator->GetGeometry()(0)->GetSolutionStepValue(PARTICLE_FRICTION);

                // CONTINUUM SIMULATING PARAMETERS:

                double initial_delta = 0.0;
                double CTension = 0.0;
                double CCohesion = 0.0;

                ///WARNING: XAPUZAAAAA
                //double& mContactFailureId_double    = int(r_VectorContactFailureId[iContactForce]);
                //int mContactFailureId               = int(mContactFailureId_double);


                array_1d<double,3>& mContactForces  = this->GetValue(PARTICLE_CONTACT_FORCES)[iContactForce];

                KRATOS_WATCH("hshshshshs")

    KRATOS_WATCH(continuum_simulation_OPTION)
     KRATOS_WATCH(continuum_group)
          KRATOS_WATCH(this->GetValue(PARTICLE_CONTACT_FAILURE_ID)[iContactForce])
          KRATOS_WATCH((continuum_simulation_OPTION && (continuum_group!=0) && (this->GetValue(PARTICLE_CONTACT_FAILURE_ID)[iContactForce]==0)))

                if (continuum_simulation_OPTION && (continuum_group!=0) && (this->GetValue(PARTICLE_CONTACT_FAILURE_ID)[iContactForce]==0))
                {
                    KRATOS_WATCH("sadjksadksksk")
                    CTension    = 2* Tension * other_tension / (Tension + other_tension);
                    CCohesion   = 2* Cohesion * other_cohesion / (Cohesion + other_cohesion);             
                KRATOS_WATCH(CTension)
                        KRATOS_WATCH(CCohesion)

                }

                if( delta_OPTION && (iContactForce < r_VectorContactInitialDelta.size()) ) {

                    initial_delta = r_VectorContactInitialDelta[iContactForce];
                }

                // BASIC CALCULATIONS

                array_1d<double,3> other_to_me_vect = this->GetGeometry()(0)->Coordinates() - neighbour_iterator->GetGeometry()(0)->Coordinates();
                double distance                     = sqrt(other_to_me_vect[0] * other_to_me_vect[0] +
                                                      other_to_me_vect[1] * other_to_me_vect[1] +
                                                      other_to_me_vect[2] * other_to_me_vect[2]);

                double radius_sum                   = radius + other_radius;

                double indentation                  = radius_sum - distance - initial_delta; //M: Here, Initial_delta is expected to be positive if it is embeding and negative if it's separation.

               // KRATOS_WATCH(this->Id())
               // KRATOS_WATCH(neighbour_iterator->Id())
           

                double equiv_radius     = 2* radius * other_radius / (radius + other_radius);
                double equiv_area       = M_PI * equiv_radius * equiv_radius;
                double equiv_poisson    = 2* poisson * other_poisson / (poisson + other_poisson);
                double equiv_young      = 2 * young * other_young / (young + other_young);
            //   double equiv_young_star = young_star * other_young_star / (young_star + other_young_star);
                double kn               = M_PI * 0.5 * equiv_young * equiv_radius; //M: CANET FORMULA
               // KRATOS_WATCH(equiv_young)
                // KRATOS_WATCH(equiv_radius)
                double ks               = kn / (2.0 * (1.0 + equiv_poisson));

                // FORMING LOCAL CORDINATES

                double NormalDir[3]           = {0.0};
                double LocalCoordSystem[3][3] = {{0.0}, {0.0}, {0.0}};
                NormalDir[0] = other_to_me_vect[0];   // M. this way the compresion is positive.
                NormalDir[1] = other_to_me_vect[1];
                NormalDir[2] = other_to_me_vect[2];
                ComputeContactLocalCoordSystem(NormalDir, LocalCoordSystem);

                // VELOCITIES AND DISPLACEMENTS

                    array_1d<double, 3 > vel            = this->GetGeometry()(0)->GetSolutionStepValue(VELOCITY);
                    array_1d<double, 3 > other_vel      = neighbour_iterator->GetGeometry()(0)->GetSolutionStepValue(VELOCITY);

                    double DeltDisp[3] = {0.0};
                    double DeltVel [3] = {0.0};

                    DeltVel[0] = (vel[0] - other_vel[0]);
                    DeltVel[1] = (vel[1] - other_vel[1]);
                    DeltVel[2] = (vel[2] - other_vel[2]);

                    //DeltDisp in global cordinates

                    DeltDisp[0] = DeltVel[0] * dt;
                    DeltDisp[1] = DeltVel[1] * dt;
                    DeltDisp[2] = DeltVel[2] * dt;

                    double LocalDeltDisp[3] = {0.0};

                    double LocalContactForce[3]  = {0.0};
                    double GlobalContactForce[3] = {0.0};

                    GlobalContactForce[0] = mContactForces[0];   //M:aqui tenim guardades les del neighbour calculator.
                    GlobalContactForce[1] = mContactForces[1];
                    GlobalContactForce[2] = mContactForces[2];

                    VectorGlobal2Local(LocalCoordSystem, DeltDisp, LocalDeltDisp);
                    VectorGlobal2Local(LocalCoordSystem, GlobalContactForce, LocalContactForce);

             // FORCES

                    //M: aqui estem guardant la força total local, que es la suma de la total que teniem més l'increment, es el metode incremental.
                    //mcontact forces esta inicialitzat per mitjà de GetContactForces en el NeighboursCalculator.
                    //aqui em falta lu de pi i tal.... oi???
                    //M: these are arrays of vectors of 3 components.

                    KRATOS_WATCH(this->Id())
                    KRATOS_WATCH(neighbour_iterator->Id())
                    KRATOS_WATCH(this->GetValue(PARTICLE_CONTACT_FAILURE_ID)[iContactForce])
                    if ( (indentation > 0.0) || (this->GetValue(PARTICLE_CONTACT_FAILURE_ID)[iContactForce] == 0) )  // This conditions take in acount the fact that the particles must remember their initial delta's between initial neighbours.
                    //M: te una mica de tela entendre bé aquestes condicions pero crec que axi es compleixien tots el casos.
                    {
                        LocalContactForce[0] += - ks * LocalDeltDisp[0];  // 0: first tangential
                        LocalContactForce[1] += - ks * LocalDeltDisp[1];  // 1: second tangential
                        LocalContactForce[2] += - kn * LocalDeltDisp[2];  // 2: normal force
                              //  KRATOS_WATCH("primera")

                    }

                    //ABSOLUTE METHOD FOR NORMAL FORCE (Allows non-linearity)
                 //   KRATOS_WATCH(indentation)

                    if(type_id == 2)  //M: estem sobreeescribim els localContact que haviem fet, si es tipus 1 ja no ho reescribim.
                        //el feng ha pensat lu del type id per si volem incremental o absolut.   1--- incremental; 2 --- absolut i amb el cas hertzià
                    {

                        if(indentation > 0.0)
                        {
                            ///Cfeng:for compression,hertzian model
                            LocalContactForce[2] = kn * pow(indentation, 1.5);

                             //         KRATOS_WATCH("segona")
                       }
                        else
                        {
                      
                            LocalContactForce[2] = kn * indentation;
                        }
                    }

                  //  KRATOS_WATCH(LocalContactForce[2])


                    // per guardar els maxims i minims.


              // TENSION FAILURE
                        if (-LocalContactForce[2] > (CTension * equiv_area))  //M:si la tensio supera el limit es seteja tot a zero.
                        {
                            KRATOS_WATCH("HA FALLAT A TRACCIO")
                            LocalContactForce[0]  = 0.0;
                            LocalContactForce[1]  = 0.0;
                            LocalContactForce[2]  = 0.0;

                            this->GetValue(PARTICLE_CONTACT_FAILURE_ID)[iContactForce] = 3.0;  //tensile failure case.
                          //mContactFailureId        = 3;
                            //M: podria ficar que ara la cohesion = dinamic cohesion, not static... ??
                        }

                        // SHEAR FAILURE

                        else
                        {

                            double ShearForceMax = LocalContactForce[2] * Friction + CCohesion * equiv_area;  // MOHR COULOMB MODEL.
                            double ShearForceNow = sqrt(LocalContactForce[0] * LocalContactForce[0]
                                                 +      LocalContactForce[1] * LocalContactForce[1]);  //M: combinació de les dues direccions tangencials.

                            //Not normal contribution for the tensile case

                            if(LocalContactForce[2] < 0.0)
                            {
                                ShearForceMax = CCohesion * equiv_area;
                            }

                            //No cohesion or friction, no shear resistance

                            if(ShearForceMax == 0.0)
                            {
                                LocalContactForce[0] = 0.0;
                                LocalContactForce[1] = 0.0;
                            }

                            else if(ShearForceNow > ShearForceMax)     // the actual shear force is actualized to the maximum one, never over; also decomposed in the same directions.
                            {
                                KRATOS_WATCH("HA FALLAT A TALLANT")
                                LocalContactForce[0] = ShearForceMax / ShearForceNow * LocalContactForce[0];
                                LocalContactForce[1] = ShearForceMax / ShearForceNow * LocalContactForce[1];

                                //mContactFailureId = 4; // Shear failure case.
                               this->GetValue(PARTICLE_CONTACT_FAILURE_ID)[iContactForce] = 4.0;
                                //M: podria ficar que ara la cohesion = dinamic cohesion, not static... ??
                            }
                        }

                // VISCODAMPING (applyied locally)

                    if (damp_id == 2)
                    {
                        double visco_damping[3] = {0,0,0};
                        //the damping is never larger than the force.

                        //M: el damping tangencial dona petits problemes... cal realment un damping?

                        /*
                        if( abs(equiv_visc_damp_ratio * DeltVel[0]) > abs(LocalContactForce[0]) )   {visco_damping[0]= LocalContactForce[0]; }
                        else { visco_damping[0]= equiv_visc_damp_ratio * DeltVel[0]; }

                        if( abs(equiv_visc_damp_ratio * DeltVel[1]) > abs(LocalContactForce[1]) )   {visco_damping[1]= LocalContactForce[1]; }
                        else { visco_damping[1]= equiv_visc_damp_ratio * DeltVel[1]; }
                        */

                        if( abs(equiv_visc_damp_ratio * DeltVel[2]) > abs(LocalContactForce[2]) )   {visco_damping[2]= LocalContactForce[2]; }
                        else { visco_damping[2]= equiv_visc_damp_ratio * DeltVel[2]; }


                        LocalContactForce[0] = LocalContactForce[0] - visco_damping[0];
                        LocalContactForce[1] = LocalContactForce[1] - visco_damping[1];
                        LocalContactForce[2] = LocalContactForce[2] - visco_damping[2];

                    }

                // TRANSFORMING TO GLOBAL FORCES AND ADDING UP

                    VectorLocal2Global(LocalCoordSystem, LocalContactForce, GlobalContactForce);

                    force[0] += GlobalContactForce[0];
                    force[1] += GlobalContactForce[1];
                    force[2] += GlobalContactForce[2];

                // SAVING INTO THE LOCAL SYSTEM ARRAYS FOR NEXT STEPS

                    mContactForces[0] = GlobalContactForce[0];
                    mContactForces[1] = GlobalContactForce[1];
                    mContactForces[2] = GlobalContactForce[2];

                    iContactForce++;
           // KRATOS_WATCH("despres de calcular al forsa")
           // KRATOS_WATCH(force)
            }//for each neaighbour


//AFTER FORCES CALCULATIONS

            // DETECTING A COHESIVE PARTICLE THAT HAS BEEN COMPLETELY DETACHED.

            double num_of_detached_neighbours=0;

            for (unsigned int i=0; i< r_VectorContactFailureId.size(); i++ ) {
                if   ( int(r_VectorContactFailureId[i]) != 0)
                    {
                    num_of_detached_neighbours++;
                    }
            }

            // LOCAL DAMPING OPTION FOR THE UNBALANCED FORCES (IN GLOBAL CORDINATES).

            if (damp_id == 1)
            {
                for (int iDof = 0; iDof < 3; iDof++)
                {
                    if (this->GetGeometry()(0)->GetSolutionStepValue(VELOCITY)[iDof] > 0.0)
                    {
                        force[iDof] = force[iDof] - LocalDampRatio * fabs(force[iDof]);
                    }
                    else
                    {
                        force[iDof] = force[iDof] + LocalDampRatio * fabs(force[iDof]);
                    }
                }
            }

            // FOR PARTICLES WITH NO NEIGHBOURS

            if (iContactForce == 0)
            {
                this->GetGeometry()(0)->GetSolutionStepValue(PARTICLE_TENSION) = 0.0;
                this->GetGeometry()(0)->GetSolutionStepValue(PARTICLE_COHESION) = 0.0;

            }

            // SETTING THE CARACHTERISTIC FAILURE TYPE OF THE PARTICLE (THE DOMINANT ONE) FOR VISUALITZATION.

            /* ContactFailureId = 2 (partially detached) it's never considerer in a contact, but it can be a desired output on the visualitzation.
             *
             * If the particle has no neighbour, but it is not from ContinuumGroup = 0 its becouse it has lost all the neighbours,
             * we are not interested in its global failure type but in its neighbours failure type. See: //// DETECTING A COHESIVE PARTICLE THAT HAS BEEN COMPLETELY DETACHED. in AfterForceCalculation.
             *
             * If one of the failure types completely dominates over the others this is the FailureId choosen for the particle
             *
             * If its partially detached by some type(s) of failure but the dominant contact state is "still atached", so failureId =0, we chose FailureId = 2.
             *
             * If the particle is not simulating the continuum, the failure type is 1 (default set from the particle whose ContinuumGroup is 0); It won't be printed as output.
             *
             */


    ///MRMRM: MIQUEL ACTIVA AIXO AMB LA REFERENCIA A CONTACT FAILURE O EL NOM EXACTE IGUAL KE ABANS... HO TENS PROVIDONAL
/*
            if (mFailureId != 1)  // for mFailureId == 1 there's no failure to represent, the particle is not a continuum-simulating one or has been completelly detached already.
            {

                int tempType[5] = {0,0,0,0,0};
                for(size_t itype = 0; itype < r_VectorContactFailureId.size(); itype++)
                {
                    if( int(r_VectorContactFailureId[itype]) == 0)
                    {
                        tempType[0]++;
                    }
                    else if( int(r_VectorContactFailureId[itype]) == 1)
                    {
                        tempType[1]++;
                    }

                    // mContactFailureId == 2 intentionally skipped!

                    else if( int(r_VectorContactFailureId[itype]) == 3)
                    {
                        tempType[3]++;
                    }
                    else if( int(r_VectorContactFailureId[itype]) == 4)
                    {
                        tempType[4]++;
                    }

                }

                if ( tempType[0] == 0)
                {
                     mFailureId = 1.0;

                }   // no one neighbore is attached (but maybe still contacting).

                else if( (tempType[3] > tempType[4]) )
                {
                    mFailureId = 3.0;
                }

                else if( (tempType[4] > tempType[3]) )
                {
                    mFailureId = 4.0;
                }
                else if ( (tempType[4] > 0) || (tempType[3] > 0) )
                {

                   mFailureId = 2.0;  // Partially detached / mix case.
                }
                else
                {

                   mFailureId = 0.0;
                }

            }// if (mFailureId != 1)

            //M: si mFailureId = 1 already I think we dont need to do GetFailureId() = mFailureId; becouse the h declatation will be applyied for new objects of the class.

 */
            KRATOS_CATCH("")

        }//ComputeForcesGeneral


        void SphericParticle::SphericParticle::DampMatrix(MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo){}

        void SphericParticle::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo){

            ElementalDofList.resize( 0 );

            for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
            {
                ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
                ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );

                if ( GetGeometry().WorkingSpaceDimension() == 3 )
                {
                    ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );
                }
            }
        }

        void SphericParticle::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
        {

          int case_opt         = rCurrentProcessInfo[CASE_OPTION];

         
          int mSwitch          = rCurrentProcessInfo[DUMMY_SWITCH];
         
          if( (mSwitch==0) && (case_opt!=0) ){

 
                  SetInitialContacts(case_opt);  //si finalment nomes has de fer aixo el switch el pots ficar a la strategia i testalvies que i entrem cada cop a comprobar.

          }
          
        }
        void SphericParticle::FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo){}

        void SphericParticle::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo){

            if(rVariable==DUMMY_FORCES) {

                ComputeForcesGeneral(rCurrentProcessInfo); 
            }
        }

        void SphericParticle::Calculate(const Variable<array_1d<double, 3 > >& rVariable, array_1d<double, 3 > & Output, const ProcessInfo& rCurrentProcessInfo){}
        void SphericParticle::Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo){}
        void SphericParticle::Calculate(const Variable<Matrix >& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo){}
  

        void SphericParticle::norm(double Vector[3])
        {
            double distance = sqrt(Vector[0] * Vector[0] + Vector[1] * Vector[1] + Vector[2] * Vector[2]);
            Vector[0] = Vector[0] / distance;
            Vector[1] = Vector[1] / distance;
            Vector[2] = Vector[2] / distance;
        }

        void SphericParticle::VectorGlobal2Local(double LocalCoordSystem[3][3], double GlobalVector[3], double LocalVector[3])
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

        void SphericParticle::VectorLocal2Global(double LocalCoordSystem[3][3], double LocalVector[3], double GlobalVector[3])
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

        double SphericParticle::DotProduct(double Vector1[3], double Vector2[3])
        {
            return Vector1[0] * Vector2[0] + Vector1[1] * Vector2[1] + Vector1[2] * Vector2[2];
        }

        void SphericParticle::CrossProduct(double u[3], double v[3], double ReturnVector[3])
        {
            ReturnVector[0] = u[1]*v[2] - u[2]*v[1];
            ReturnVector[1] = v[0]*u[2] - u[0]*v[2];
            ReturnVector[2] = u[0]*v[1] - u[1]*v[0];
        }

        void SphericParticle::ComputeContactLocalCoordSystem(double NormalDirection[3], double LocalCoordSystem[3][3])
        {
            int ii;
            double Vector0[3] = {0.0},Vector1[3] = {0.0};

            norm(NormalDirection);

            double fix_coord[3]={0.0};

            //Ax+By+Cz=D
            double x0,y0,z0,D;

           // D=fix_coord[0]*NormalDirection[0] + fix_coord[1]*NormalDirection[1] +fix_coord[2]*NormalDirection[2];

            D = DotProduct(fix_coord, NormalDirection);

            if(fabs(NormalDirection[0])>=0.577)
            {
                y0=fix_coord[1]+1.0;
                z0=fix_coord[2];
                x0=( D-NormalDirection[1]*y0-NormalDirection[2]*z0 )/NormalDirection[0];
                Vector0[0]=x0-fix_coord[0];
                Vector0[1]=1.0;
                Vector0[2]=0.0;
            }
            else if(fabs(NormalDirection[1])>=0.577)
            {
                x0=fix_coord[0];
                z0=fix_coord[2]+1.0;
                y0=( D-NormalDirection[0]*x0-NormalDirection[2]*z0 )/NormalDirection[1];
                Vector0[0]=0.0;
                Vector0[1]=y0-fix_coord[1];
                Vector0[2]=1.0;
            }
            else
            {
                x0=fix_coord[0]+1.0;
                y0=fix_coord[1];
                z0=( D-NormalDirection[0]*x0-NormalDirection[1]*y0 )/NormalDirection[2];
                Vector0[0]=1.0;
                Vector0[1]=0.0;
                Vector0[2]=z0-fix_coord[2];
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

}  // namespace Kratos.


