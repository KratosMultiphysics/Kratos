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
#include "custom_utilities/GeometryFunctions.h"
#include "DEM_application.h"


namespace Kratos
{
     // using namespace GeometryFunctions;



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

        double & Inertia         = GetGeometry()(0)->FastGetSolutionStepValue(PARTICLE_INERTIA);
        double & MomentOfInertia = GetGeometry()(0)->FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA);


        mContinuumGroup     = this->GetGeometry()[0].GetSolutionStepValue(PARTICLE_CONTINUUM);
        mFailureId          = !(mContinuumGroup);
           
        if(mDimension ==2)
        {
            mass     = M_PI * radius * radius * density;

            mRealMass = mass;

            Inertia = 0.25 * M_PI * radius * radius * radius  * radius ;

            MomentOfInertia = 0.5 * radius * radius;

        }
        else
        {
            mass     = 4.0 / 3.0 * M_PI * radius * radius * radius * density;

            mRealMass = mass;

            Inertia = 0.25 * M_PI * radius * radius * radius  * radius ;

            MomentOfInertia = 0.4 * radius * radius;

        }


        KRATOS_CATCH( "" )

      }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************


      void SphericParticle::CalculateRightHandSide(VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo){

        ComputeParticleContactForce(rCurrentProcessInfo);


        if( (rCurrentProcessInfo[ROTATION_OPTION] != 0) && (rCurrentProcessInfo[ROTATION_SPRING_OPTION] != 0) )
        {
              ComputeParticleRotationSpring(rCurrentProcessInfo);
        }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

      }
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

           //vector<double>& r_VectorContactFailureId            = this->GetValue(PARTICLE_CONTACT_FAILURE_ID); //M:temporarily double...
           //r_VectorContactFailureId.resize(r_neighbours.size());

           //vector<double>& r_VectorContactInitialDelta       = this->GetValue(PARTICLE_CONTACT_INITIAL_DELTA);
           //r_VectorContactInitialDelta.resize(r_neighbours.size()); //the initial neighbours are the neighbours at the first step.

           ParticleWeakVectorType& r_neighbours             = this->GetValue(NEIGHBOUR_ELEMENTS);

           //this->GetValue(INITIAL_NEIGHBOUR_ELEMENTS).resize(r_neighbours.size());
           this->GetValue(PARTICLE_INITIAL_DELTA).resize(r_neighbours.size());

           ParticleWeakVectorType& r_initial_neighbours     = this->GetValue(INITIAL_NEIGHBOUR_ELEMENTS);

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

                        if( (r_other_continuum_group == mContinuumGroup) || ( fabs(initial_delta)>1.0e-6 ) ) // R: TOLERANCIA PER DIR QUE TENEN UN IDENTACIÓ O SEPARACIÓ; HAURIA DANAR LLIGADA AMB LU DE LA BUSQUEDA. la busqueda es fa per intersecció, es a dir identació es detecta pero si estan separades un pel ja no es troba com a veí.
                        //THESE ARE THE CASES THAT NEED TO STORE THE INITIAL NEIGHBOURS
                        {

                            r_initial_neighbours.push_back(*ineighbour);

                            this->GetValue(PARTICLE_INITIAL_DELTA)[i]  =   initial_delta; //M:R:els guardo sempre, en forces i tal si és car accedir aquest valor podem demanar el flag pero no crec ke sigui car.
                            this->GetValue(PARTICLE_CONTACT_DELTA)[i]  =   initial_delta;

                            if (r_other_continuum_group == mContinuumGroup && (mContinuumGroup != 0) ) {this->GetValue(PARTICLE_CONTACT_FAILURE_ID)[i]=0; }
                            else                                             this->GetValue(PARTICLE_CONTACT_FAILURE_ID)[i]=1; // MRMR:  //generally detached    //diferent group
                                                                                                                             //hi havia: r_VectorContactFailureId[i]=1; //generally detached    //diferent group



                        } // FOR THE CASES THAT NEED STORING INITIAL NEIGHBOURS

                        else mFailureId=1;      //NO NEED to store the initial neighbour;
                        //M:also we can say that the particle is detached --> this can be a flag for the forces.

                        i++;

               }//if I found myself.

            } //end for: ParticleWeakIteratorType ineighbour
      }//SET INITIAL CONTACTS.


        void SphericParticle::ComputeParticleContactForce(const ProcessInfo& rCurrentProcessInfo )

        {

            KRATOS_TRY

            ParticleWeakVectorType& r_neighbours             = this->GetValue(NEIGHBOUR_ELEMENTS);

          // vector<double>& r_VectorContactFailureId            = this->GetValue(PARTICLE_CONTACT_FAILURE_ID); //M:temporarily a vector of doubles... it should be vector of ints..
            vector<double>& r_VectorContactInitialDelta         = this->GetValue(PARTICLE_CONTACT_DELTA);

            // PROCESS INFO

            const array_1d<double,3>& gravity   = rCurrentProcessInfo[GRAVITY];

            double dt                           = rCurrentProcessInfo[DELTA_TIME];
            int damp_id                         = rCurrentProcessInfo[DAMP_TYPE];
            int type_id                         = rCurrentProcessInfo[FORCE_CALCULATION_TYPE];
            int rotation_OPTION                 = rCurrentProcessInfo[ROTATION_OPTION]; //M:  it's 1/0, should be a boolean

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

            double Tension          = this->GetGeometry()[0].GetSolutionStepValue(PARTICLE_TENSION);
            double Cohesion         = this->GetGeometry()[0].GetSolutionStepValue(PARTICLE_COHESION);
            double FriAngle         = this->GetGeometry()[0].GetSolutionStepValue(PARTICLE_FRICTION); //M: CFENG: es aixo el angle de friccio?
            double Friction         = tan( FriAngle / 180.0 * M_PI);

            double radius               = this->GetGeometry()[0].GetSolutionStepValue(RADIUS);
            double critic_damp_fraction = this->GetGeometry()[0].GetSolutionStepValue(VISCO_DAMP_COEFF);
            double mass                 = mRealMass;

            double young                = this->GetGeometry()[0].GetSolutionStepValue(YOUNG_MODULUS);
            double poisson              = this->GetGeometry()[0].GetSolutionStepValue(POISSON_RATIO);

            array_1d<double,3>& force           = this->GetGeometry()[0].GetSolutionStepValue(RHS);//total forces, we reset to 0. and we calculate again.
            //COMPROBAR QUE ES ZERO DEL INITIALIZE SOLUTION STEP

            array_1d<double,3> applied_force    = this->GetGeometry()[0].GetSolutionStepValue(APPLIED_FORCE); //Nelson: se llama force la que hauria d0haver aqui

            force  = mass*gravity + applied_force;
            //KRATOS_WATCH(mass)
            //KRATOS_WATCH(gravity)
            //KRATOS_WATCH(force)
            //KRATOS_WATCH(applied_force)

            array_1d<double, 3 > & mRota_Moment = this->GetGeometry()[0].GetSolutionStepValue(PARTICLE_MOMENT);


            size_t iContactForce = 0;


            for(ParticleWeakIteratorType neighbour_iterator = r_neighbours.begin();
            neighbour_iterator != r_neighbours.end(); neighbour_iterator++)
            {
             // GETTING NEIGHBOUR PROPERTIES


             
                double other_radius                 = neighbour_iterator->GetGeometry()(0)->GetSolutionStepValue(RADIUS);
                double other_critic_damp_fraction   = neighbour_iterator->GetGeometry()(0)->GetSolutionStepValue(VISCO_DAMP_COEFF);
                double equiv_visc_damp_ratio        = (critic_damp_fraction + other_critic_damp_fraction) / 2.0;   //M: is it correct to be a simple mean.
                // double other_mass                   = neighbour_iterator.mRealMass;
                double other_young                  = neighbour_iterator->GetGeometry()[0].GetSolutionStepValue(YOUNG_MODULUS);
                double other_poisson                = neighbour_iterator->GetGeometry()[0].GetSolutionStepValue(POISSON_RATIO);
                double other_tension                = neighbour_iterator->GetGeometry()[0].GetSolutionStepValue(PARTICLE_TENSION);
                double other_cohesion               = neighbour_iterator->GetGeometry()[0].GetSolutionStepValue(PARTICLE_COHESION);
                double other_FriAngle               = neighbour_iterator->GetGeometry()[0].GetSolutionStepValue(PARTICLE_FRICTION);

                // CONTINUUM SIMULATING PARAMETERS:

                double initial_delta = 0.0;
                double CTension = 0.0;
                double CCohesion = 0.0;

                ///WARNING: XAPUZAAAAA
                //double& mContactFailureId_double    = int(r_VectorContactFailureId[iContactForce]);
                //int mContactFailureId               = int(mContactFailureId_double);


                array_1d<double,3>& mContactForces  = this->GetValue(PARTICLE_CONTACT_FORCES)[iContactForce];



                if (continuum_simulation_OPTION && (continuum_group!=0) && (this->GetValue(PARTICLE_CONTACT_FAILURE_ID)[iContactForce]==0))
                {
                    // Test for Average, 120531, should be discussed
                    CTension  = (Tension + other_tension)   * 0.5;
                    CCohesion = (Cohesion + other_cohesion) * 0.5;
                    /*
                    CTensionUP    = 2* Tension * other_tension;
                    if(CTensionUP==0){

                        CTension=0;
                    }
                    else {

                    CTension    = 2* Tension * other_tension / (Tension + other_tension);

                    }
                    CCohesion   = 2* Cohesion * other_cohesion / (Cohesion + other_cohesion);
                     *
                    */
                }


                if( delta_OPTION && (iContactForce < r_VectorContactInitialDelta.size()) )
                {

                    initial_delta = r_VectorContactInitialDelta[iContactForce];
                }

                // BASIC CALCULATIONS

                array_1d<double,3> other_to_me_vect = this->GetGeometry()(0)->Coordinates() - neighbour_iterator->GetGeometry()(0)->Coordinates();
                double distance                     = sqrt(other_to_me_vect[0] * other_to_me_vect[0] +
                                                      other_to_me_vect[1] * other_to_me_vect[1] +
                                                      other_to_me_vect[2] * other_to_me_vect[2]);

                double radius_sum                   = radius + other_radius;

                double indentation                  = radius_sum - distance - initial_delta; //M: Here, Initial_delta is expected to be positive if it is embeding and negative if it's separation.

                double equiv_radius     = 2* radius * other_radius / (radius + other_radius);
                double equiv_area       = M_PI * equiv_radius * equiv_radius;
                double equiv_poisson    = 2* poisson * other_poisson / (poisson + other_poisson);
                double equiv_young      = 2 * young * other_young / (young + other_young);

                Friction                = tan( (FriAngle + other_FriAngle) * 0.5  / 180.0 * M_PI);

                double kn               = M_PI * 0.5 * equiv_young * equiv_radius; //M: CANET FORMULA
                double ks               = kn / (2.0 * (1.0 + equiv_poisson));
              
                //OÑATE. PROBETES.
                 //kn               = 197760000;
                 //ks               = 51593000;

                // FORMING LOCAL CORDINATES

                double NormalDir[3]           = {0.0};
                double LocalCoordSystem[3][3] = {{0.0}, {0.0}, {0.0}};
                NormalDir[0] = other_to_me_vect[0];   // M. this way the compresion is positive.
                NormalDir[1] = other_to_me_vect[1];
                NormalDir[2] = other_to_me_vect[2];
                GeometryFunctions::ComputeContactLocalCoordSystem(NormalDir, LocalCoordSystem);

                // VELOCITIES AND DISPLACEMENTS

                array_1d<double, 3 > vel            = this->GetGeometry()(0)->GetSolutionStepValue(VELOCITY);
                array_1d<double, 3 > other_vel      = neighbour_iterator->GetGeometry()(0)->GetSolutionStepValue(VELOCITY);

                array_1d<double, 3 > delta_displ            = this->GetGeometry()(0)->GetSolutionStepValue(DELTA_DISPLACEMENT);
                array_1d<double, 3 > other_delta_displ      = neighbour_iterator->GetGeometry()(0)->GetSolutionStepValue(DELTA_DISPLACEMENT);

                double DeltDisp[3] = {0.0};
                double RelVel [3] = {0.0};

                RelVel[0] = (vel[0] - other_vel[0]);
                RelVel[1] = (vel[1] - other_vel[1]);
                RelVel[2] = (vel[2] - other_vel[2]);

                //DeltDisp in global cordinates
  /*
               DeltDisp[0] = RelVel[0] * dt;
               DeltDisp[1] = RelVel[1] * dt;
               DeltDisp[2] = RelVel[2] * dt;
*/
               DeltDisp[0] = (delta_displ[0] - other_delta_displ[0]);
               DeltDisp[1] = (delta_displ[1] - other_delta_displ[1]);
               DeltDisp[2] = (delta_displ[2] - other_delta_displ[2]);
                
                    if ( rotation_OPTION == 1 )
                    {

                        double velA[3]   = {0.0};
                        double velB[3]   = {0.0};
                        double dRotaDisp[3] = {0.0};

                        array_1d<double, 3 > AngularVel       = this->GetGeometry()(0)->FastGetSolutionStepValue(ANGULAR_VELOCITY);
                        array_1d<double, 3 > Other_AngularVel = neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(ANGULAR_VELOCITY);

                        double Vel_Temp[3]       = {      AngularVel[0],       AngularVel[1],       AngularVel[2]};
                        double Other_Vel_Temp[3] = {Other_AngularVel[0], Other_AngularVel[1], Other_AngularVel[2]};
                        GeometryFunctions::CrossProduct(Vel_Temp,             LocalCoordSystem[2], velA);
                        GeometryFunctions::CrossProduct(Other_Vel_Temp, LocalCoordSystem[2], velB);

                        dRotaDisp[0] = -velA[0] * radius - velB[0] * other_radius;
                        dRotaDisp[1] = -velA[1] * radius - velB[1] * other_radius;
                        dRotaDisp[2] = -velA[2] * radius - velB[2] * other_radius;
                        //////contribution of the rotation vel
                        DeltDisp[0] += dRotaDisp[0] * dt;
                        DeltDisp[1] += dRotaDisp[1] * dt;
                        DeltDisp[2] += dRotaDisp[2] * dt;

                    }//if rotation_OPTION

                    double LocalDeltDisp[3] = {0.0};
                    double LocalContactForce[3]  = {0.0};
                    double GlobalContactForce[3] = {0.0};
                    //double GlobalContactForceOld[3] = {0.0};


                    GlobalContactForce[0] = mContactForces[0];   //M:aqui tenim guardades les del neighbour calculator.
                    GlobalContactForce[1] = mContactForces[1];
                    GlobalContactForce[2] = mContactForces[2];

                    GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, DeltDisp, LocalDeltDisp);
                    GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, GlobalContactForce, LocalContactForce);

             // FORCES

                    //M: aqui estem guardant la força total local, que es la suma de la total que teniem més l'increment, es el metode incremental.
                    //mcontact forces esta inicialitzat per mitjà de GetContactForces en el NeighboursCalculator.
                    //aqui em falta lu de pi i tal.... oi???
                    //M: these are arrays of vectors of 3 components.

                 if ( (indentation > 0.0) || (this->GetValue(PARTICLE_CONTACT_FAILURE_ID)[iContactForce] == 0) )  // for detached particles we enter only if the indentation is > 0.
                                                                                                                  // for attached particles we enter onlty if the particle is still attached.
                    {
                        LocalContactForce[0] += - ks * LocalDeltDisp[0];  // 0: first tangential
                        LocalContactForce[1] += - ks * LocalDeltDisp[1];  // 1: second tangential
                        LocalContactForce[2] += - kn * LocalDeltDisp[2];  // 2: normal force
                      
                        //ABSOLUTE METHOD FOR NORMAL FORCE (Allows non-linearity)

                        if(type_id == 2)  //  1--- incremental; 2 --- absolut i amb el cas hertzià
                        {
                            ///Cfeng:for compression,hertzian model
                            LocalContactForce[2] = kn * pow(indentation, 1.5);

                            if(indentation > 0.0) // usually theory considers the tension to be linear and compression non linear
                            {

                             LocalContactForce[2] = kn * indentation;

                            }

                        }
             
                    }

                    //ABSOLUTE METHOD FOR NORMAL FORCE (Allows non-linearity)

                    if(type_id == 2)  //M: estem sobreeescribim els localContact que haviem fet, si es tipus 1 ja no ho reescribim.
                        //  1--- incremental; 2 --- absolut i amb el cas hertzià
                    {

                        if(indentation > 0.0)
                        {
                            ///Cfeng:for compression,hertzian model
                            LocalContactForce[2] = kn * pow(indentation, 1.5);

                       }
                        else
                        {
                            LocalContactForce[2] = kn * indentation;
                        }
                    }

                  if ( (indentation < 0.0) )
                    {
                        LocalContactForce[0] = 0.0;
                        LocalContactForce[1] = 0.0;
                        LocalContactForce[2] = 0.0;
                    }


              // TENSION FAILURE

                        if (-LocalContactForce[2] > (CTension * equiv_area))  //M:si la tensio supera el limit es seteja tot a zero.
                      //OÑATE. PROBETES.
                       // if (-LocalContactForce[2] > (649.85))
                        {
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
                                //OÑATE //AMB FRICCIO
                            // double ShearForceMax = 5524.9 + 0.84*LocalContactForce[2];
                           double ShearForceMax = LocalContactForce[2] * Friction + CCohesion * equiv_area;  // MOHR COULOMB MODEL.
                            double ShearForceNow = sqrt(LocalContactForce[0] * LocalContactForce[0]
                                                 +      LocalContactForce[1] * LocalContactForce[1]);  //M: combinació de les dues direccions tangencials.

                            //Not normal contribution for the tensile case

                            if(LocalContactForce[2] < 0.0)
                            {
                                ShearForceMax = CCohesion * equiv_area;
                                //OÑATE
                                //ShearForceMax = 5524.9;

                            }

                            //No cohesion or friction, no shear resistance

                            if(ShearForceMax == 0.0)
                            {
                                LocalContactForce[0] = 0.0;
                                LocalContactForce[1] = 0.0;
                            }

                            else if(ShearForceNow > ShearForceMax)     // the actual shear force is actualized to the maximum one, never over; also decomposed in the same directions.
                            {
                                LocalContactForce[0] = ShearForceMax / ShearForceNow * LocalContactForce[0];
                                LocalContactForce[1] = ShearForceMax / ShearForceNow * LocalContactForce[1];

                                //mContactFailureId = 4; // Shear failure case.
                               this->GetValue(PARTICLE_CONTACT_FAILURE_ID)[iContactForce] = 4.0;
                                //M: podria ficar que ara la cohesion = dinamic cohesion, not static... ??
                            }
                        }

                // VISCODAMPING (applyied locally)

                        if (damp_id == 2 || damp_id == 3 )
                        {

                            double visco_damping[3] = {0,0,0};
                            //the damping is never larger than the force.

                            //M: el damping tangencial dona petits problemes... cal realment un damping?

                            /*
                            if( abs(equiv_visc_damp_ratio * RelVel[0]) > abs(LocalContactForce[0]) )   {visco_damping[0]= LocalContactForce[0]; }
                            else { visco_damping[0]= equiv_visc_damp_ratio * RelVel[0]; }

                            if( abs(equiv_visc_damp_ratio * RelVel[1]) > abs(LocalContactForce[1]) )   {visco_damping[1]= LocalContactForce[1]; }
                            else { visco_damping[1]= equiv_visc_damp_ratio * RelVel[1]; }
                            */

                            if( abs(equiv_visc_damp_ratio * RelVel[2]) > abs(LocalContactForce[2]) )   {visco_damping[2]= LocalContactForce[2]; }
                            else { visco_damping[2]= equiv_visc_damp_ratio * RelVel[2]; }


                            LocalContactForce[0] = LocalContactForce[0] - visco_damping[0];
                            LocalContactForce[1] = LocalContactForce[1] - visco_damping[1];
                            LocalContactForce[2] = LocalContactForce[2] - visco_damping[2];

                        }


                // TRANSFORMING TO GLOBAL FORCES AND ADDING UP

                    GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalContactForce, GlobalContactForce);

                    force[0] += GlobalContactForce[0];
                    force[1] += GlobalContactForce[1];
                    force[2] += GlobalContactForce[2];

                // SAVING INTO THE LOCAL SYSTEM ARRAYS FOR NEXT STEPS

                    mContactForces[0] = GlobalContactForce[0];
                    mContactForces[1] = GlobalContactForce[1];
                    mContactForces[2] = GlobalContactForce[2];

                    if ( rotation_OPTION == 1 )
                    {

                 ////global moment return back,for the ball self

                    double MA[3] = {0.0};
                    GeometryFunctions::CrossProduct(LocalCoordSystem[2], GlobalContactForce, MA);
                    mRota_Moment[0] -= MA[0] * radius;
                    mRota_Moment[1] -= MA[1] * radius;
                    mRota_Moment[2] -= MA[2] * radius;

                    }

                    iContactForce++;

            }//for each neaighbour

            KRATOS_CATCH("")

        }//ComputeParticleContactForce


      void SphericParticle::ApplyLocalForcesDamping(const ProcessInfo& rCurrentProcessInfo )

      {

          array_1d<double,3>& force           = this->GetGeometry()[0].GetSolutionStepValue(RHS);
          double LocalDampRatio   = this->GetGeometry()[0].GetSolutionStepValue(PARTICLE_LOCAL_DAMP_RATIO);

              // LOCAL DAMPING OPTION FOR THE UNBALANCED FORCES (IN GLOBAL CORDINATES).

                for (int iDof = 0; iDof < 3; iDof++)
                {
                    //if (this->GetGeometry()(0)->GetSolutionStepValue(VELOCITY)[iDof] > 0.0)
                    /*if (this->GetGeometry()(0)->GetSolutionStepValue(RHS)[iDof] > 0.0)
                    {
                        force[iDof] = force[iDof] - LocalDampRatio * fabs(force[iDof]);
                    }
                    else
                    {
                        force[iDof] = force[iDof] + LocalDampRatio * fabs(force[iDof]);
                    }*/
                    //KRATOS_WATCH(force[iDof])
                    force[iDof] = force[iDof] - LocalDampRatio * force[iDof];
                    //KRATOS_WATCH(force[iDof])
                }

      } //ApplyLocalForcesDamping

    void SphericParticle::ApplyLocalMomentsDamping(const ProcessInfo& rCurrentProcessInfo )

    {
       // int damp_id                         = rCurrentProcessInfo[DAMP_TYPE];  //M: revisable....en principi akest era per visco o local de forces no per moment.

        array_1d<double, 3 > & RotaMoment       = this->GetGeometry()[0].GetSolutionStepValue(PARTICLE_MOMENT);
        double LocalDampRatio                   = this->GetGeometry()[0].GetSolutionStepValue(PARTICLE_LOCAL_DAMP_RATIO);

        // LOCAL DAMPING OPTION FOR THE UNBALANCED FORCES (IN GLOBAL CORDINATES).


        for (int iDof = 0; iDof < 3; iDof++)
        {
            if (this->GetGeometry()(0)->GetSolutionStepValue(ANGULAR_VELOCITY)[iDof] > 0.0)
            {
                 RotaMoment[iDof] = RotaMoment[iDof] - LocalDampRatio * fabs(RotaMoment[iDof]);

            }
            else
            {
                 RotaMoment[iDof] = RotaMoment[iDof] + LocalDampRatio * fabs(RotaMoment[iDof]);
            }
        }


    } //ApplyLocalMomentsDamping



        void SphericParticle::ComputeParticleRotationSpring(const ProcessInfo& rCurrentProcessInfo)
        {

        double dt                           = rCurrentProcessInfo[DELTA_TIME]; //C.F.: neew


        double Tension        = this->GetGeometry()[0].GetSolutionStepValue(PARTICLE_TENSION);
        double Cohesion       = this->GetGeometry()[0].GetSolutionStepValue(PARTICLE_COHESION);
        double young          = this->GetGeometry()[0].GetSolutionStepValue(YOUNG_MODULUS);
        double poisson        = this->GetGeometry()[0].GetSolutionStepValue(POISSON_RATIO);
        double radius         = this->GetGeometry()[0].GetSolutionStepValue(RADIUS);
        double inertia        = this->GetGeometry()[0].GetSolutionStepValue(PARTICLE_INERTIA);

        array_1d<double, 3 > & mRota_Moment = GetGeometry()(0)->FastGetSolutionStepValue(PARTICLE_MOMENT);

        //Vector & mContactForces      = this->GetValue(PARTICLE_CONTACT_FORCE);
       //esta dintre del loop

        //Vector & mIfInitalContact    = this->GetValue(PARTICLE_IF_INITIAL_CONTACT);

        //C.F//ParticleWeakVector & rE      = this->GetValue(NEIGHBOUR_PARTICLE_ELEMENTS);
        ParticleWeakVectorType& rE             = this->GetValue(NEIGHBOUR_ELEMENTS);


//        Vector & mRotaSpringMoment       = this->GetValue(PARTICLE_ROTATE_SPRING_MOMENT);  //ESTA A DINTRE

        Vector & mRotaSpringFailureType  = this->GetValue(PARTICLE_ROTATE_SPRING_FAILURE_TYPE);

        size_t iContactForce = 0;

        for(ParticleWeakIteratorType ineighbour = rE.begin(); ineighbour != rE.end(); ineighbour++)
        {


            //C.F//if(mIfInitalContact[iContactForce] == 1 && mRotaSpringFailureType[iContactForce] == 0)

            {

       array_1d<double, 3 > & mRotaSpringMoment  = this->GetValue(PARTICLE_ROTATE_SPRING_MOMENT)[ iContactForce ];
//        Vector & mRotaSpringMoment
                double other_radius    = ineighbour->GetGeometry()(0)->FastGetSolutionStepValue(RADIUS);
                double other_young     = ineighbour->GetGeometry()[0].GetSolutionStepValue(YOUNG_MODULUS);
                double other_poisson   = ineighbour->GetGeometry()[0].GetSolutionStepValue(POISSON_RATIO);
                double other_tension   = ineighbour->GetGeometry()[0].GetSolutionStepValue(PARTICLE_TENSION);
                double other_cohesion  = ineighbour->GetGeometry()[0].GetSolutionStepValue(PARTICLE_COHESION);
                double other_inertia   = ineighbour->GetGeometry()(0)->FastGetSolutionStepValue(PARTICLE_INERTIA);

                Tension  = (Tension  + other_tension ) * 0.5;
                Cohesion = (Cohesion + other_cohesion) * 0.5;

                double equiv_radius     = (radius + other_radius) * 0.5 ;
                double equiv_area       = M_PI * equiv_radius * equiv_radius;
                double equiv_poisson    = (poisson + other_poisson) * 0.5 ;
                double equiv_young      = (young  + other_young)  * 0.5;

                double kn               = equiv_young * equiv_area / (2.0 * equiv_radius);
                double ks               = kn / (2.0 * (1.0 + equiv_poisson));

                //CF: NEEEEW!
                array_1d<double,3>& mContactForces = this->GetValue(PARTICLE_CONTACT_FORCES)[ iContactForce ];
                // BY MIKEL

                array_1d<double,3> other_to_me_vect = GetGeometry()(0)->Coordinates() - ineighbour->GetGeometry()(0)->Coordinates();

               /////Cfeng: Forming the Local Contact Coordinate system
                double NormalDir[3]           = {0.0};
                double LocalCoordSystem[3][3] = {{0.0}, {0.0}, {0.0}};
                NormalDir[0] = other_to_me_vect[0];
                NormalDir[1] = other_to_me_vect[1];
                NormalDir[2] = other_to_me_vect[2];
                GeometryFunctions::ComputeContactLocalCoordSystem(NormalDir, LocalCoordSystem);

                double LocalRotaSpringMoment[3]     = {0.0};
                double GlobalRotaSpringMoment[3]    = {0.0};
                double GlobalRotaSpringMomentOld[3] = {0.0};

                //////delt rotation, pay attention to the negative sign beforce the brackets,the negtive sign depends on the Force stored in the paticle
                array_1d<double, 3 > AngularVel       = GetGeometry()(0)->FastGetSolutionStepValue(ANGULAR_VELOCITY);
                array_1d<double, 3 > Other_AngularVel = ineighbour->GetGeometry()(0)->FastGetSolutionStepValue(ANGULAR_VELOCITY);
                double DeltRotaDisp[3] = {0.0};
                DeltRotaDisp[0] = -(AngularVel[0] - Other_AngularVel[0]) * dt;  //C.F. mtimestep
                DeltRotaDisp[1] = -(AngularVel[1] - Other_AngularVel[1]) * dt;
                DeltRotaDisp[2] = -(AngularVel[2] - Other_AngularVel[2]) * dt;

                double LocalDeltRotaDisp[3] = {0.0};
                GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, DeltRotaDisp, LocalDeltRotaDisp);




                GlobalRotaSpringMomentOld[0] = mRotaSpringMoment[ 0 ];

                GlobalRotaSpringMomentOld[1] = mRotaSpringMoment[ 1 ];

                GlobalRotaSpringMomentOld[2] = mRotaSpringMoment[ 2 ];

                GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, GlobalRotaSpringMomentOld, LocalRotaSpringMoment);


                double Inertia_I = (inertia + other_inertia) * 0.5;

                double Inertia_J = Inertia_I * 2.0;



                LocalRotaSpringMoment[0] +=  - Inertia_I * LocalDeltRotaDisp[0] * kn / equiv_area;

                LocalRotaSpringMoment[1] +=  - Inertia_I * LocalDeltRotaDisp[1] * kn / equiv_area;

                LocalRotaSpringMoment[2] +=  - Inertia_J * LocalDeltRotaDisp[2] * ks / equiv_area;



                ////Judge if the rotate spring is broken or not
                double GlobalContactForce[3]  = {0.0};
                double LocalContactForce [3]  = {0.0};

                GlobalContactForce[0] = mContactForces[ 0 ];
                GlobalContactForce[1] = mContactForces[ 1 ];
                GlobalContactForce[2] = mContactForces[ 2 ]; //GlobalContactForce[2] = mContactForces[3 * iContactForce  + 2 ];
                GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, GlobalContactForce, LocalContactForce);

                double ForceN  = LocalContactForce[2];
                double ForceS  = sqrt( LocalContactForce[0] * LocalContactForce[0] + LocalContactForce[1] * LocalContactForce[1]);
                double MomentS = sqrt(LocalRotaSpringMoment[0] * LocalRotaSpringMoment[0] + LocalRotaSpringMoment[1] * LocalRotaSpringMoment[1]);
                double MomentN = LocalRotaSpringMoment[2];

                //////bending stress and axial stress add together, use edge of the bar will failure first
                double TensiMax = -ForceN / equiv_area + MomentS        / Inertia_I * equiv_radius;
                double ShearMax = ForceS  / equiv_area + fabs(MomentN)  / Inertia_J * equiv_radius;

                if(TensiMax > Tension || ShearMax > Cohesion)
                {
                    mRotaSpringFailureType[iContactForce] = 1;


                    LocalRotaSpringMoment[0] = 0.0;
                    LocalRotaSpringMoment[1] = 0.0;
                    LocalRotaSpringMoment[2] = 0.0;
                }

                GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalRotaSpringMoment, GlobalRotaSpringMoment);



                mRotaSpringMoment[ 0 ] = GlobalRotaSpringMoment[0];
                mRotaSpringMoment[ 1 ] = GlobalRotaSpringMoment[1];
                mRotaSpringMoment[ 2 ] = GlobalRotaSpringMoment[2];

                ////feedback, contact moment----induce by rotation spring
                mRota_Moment[0] -= GlobalRotaSpringMoment[0];
                mRota_Moment[1] -= GlobalRotaSpringMoment[1];
                mRota_Moment[2] -= GlobalRotaSpringMoment[2];
            }

            iContactForce++;
            }


        }//ComputeParticleRotationSpring










        void SphericParticle::DampMatrix(MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo){}

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

          if( (mSwitch==0) && (case_opt!=0) )
          {
                SetInitialContacts(case_opt);  //si finalment nomes has de fer aixo el switch el pots ficar a la strategia i testalvies que i entrem cada cop a comprobar.
          }

          array_1d<double,3>& force           = this->GetGeometry()[0].GetSolutionStepValue(RHS);//total forces, we reset to 0. and we calculate again.
          noalias(force)                      = ZeroVector(3);



        }
        void SphericParticle::FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo){}

        void SphericParticle::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo)
        {

            if (rVariable == DELTA_TIME)
            {
                double E = this->GetGeometry()(0)->FastGetSolutionStepValue(YOUNG_MODULUS);
                double K = E * M_PI * this->GetGeometry()(0)->FastGetSolutionStepValue(RADIUS);
                Output = sqrt( mRealMass / K);

                if(rCurrentProcessInfo[ROTATION_OPTION] == 1)
                {
                    Output = Output * 0.5; //factor for critical time step when rotation is allowed.
                }
            } //CRITICAL DELTA CALCULATION

            if (rVariable == PARTICLE_LOCAL_DAMP_RATIO)
            {
                int damp_id             = rCurrentProcessInfo[DAMP_TYPE];
                int rotation_OPTION     = rCurrentProcessInfo[ROTATION_OPTION];

                if (damp_id == 1 || damp_id == 3 )
                {
                   ApplyLocalForcesDamping( rCurrentProcessInfo );

                   if ( rotation_OPTION != 0 )
                   {
                       ApplyLocalMomentsDamping( rCurrentProcessInfo );
                   }
                }
            } //DAMPING
        }


        void SphericParticle::Calculate(const Variable<array_1d<double, 3 > >& rVariable, array_1d<double, 3 > & Output, const ProcessInfo& rCurrentProcessInfo){}
        void SphericParticle::Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo){}
        void SphericParticle::Calculate(const Variable<Matrix >& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo){}

}  // namespace Kratos.