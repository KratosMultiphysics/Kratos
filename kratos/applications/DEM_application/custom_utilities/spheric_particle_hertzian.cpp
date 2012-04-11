//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: M.Santasusana $
//   Date:                $Date: 2012-3-12 19:44:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//

// System includes
#include <algorithm>
#include <iostream>
#include <ostream>
#include <math.h>//M: for adding pi number.
#include <vector>
#include <boost/weak_ptr.hpp> 
#include "utilities/timer.h"

// External includes 

// Project includes
#include "custom_utilities/spheric_particle_hertzian.h"
#include "spheric_particle_hertzian.h"


namespace Kratos{


    void SphericHertzianParticle::SetInitialContacts()
    {
        mInitialNeighbours = ParticleWeakVectorType();
        mInitialDelta = std::vector<double>();
        mContactFailureId = std::vector<int>();


      for(ParticleWeakIteratorType_ptr ineighbour = mNeighbours.ptr_begin();
        ineighbour != mNeighbours.ptr_end(); ineighbour++){

            mInitialNeighbours.push_back(*ineighbour);

           array_1d<double,3> other_to_me_vect = this->GetPosition() - ((*ineighbour).lock())->GetPosition();
           double distance                     = sqrt(other_to_me_vect[0] * other_to_me_vect[0] +
                                                 other_to_me_vect[1] * other_to_me_vect[1] +
                                                 other_to_me_vect[2] * other_to_me_vect[2]);

            double radius_sum                   = this->GetRadius() + ((*ineighbour).lock())->GetRadius();
            double initial_delta                = radius_sum - distance;

            
            mInitialDelta.push_back(initial_delta);
            mContactFailureId.push_back(1); //by default "generally detached" = 1
             
        } //end for: ParticleWeakIteratorType ineighbour
     } // SetInitialContacts

    void SphericHertzianParticle::AddContinuumContacts()
    {
       
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
             */

            if(other_continuum_group != continuum_group ){
                mContactFailureId[index]=1;
              }
            else {     //same group
             mContactFailureId[index]=0;
             }
             index++;

        }//end for: ParticleWeakIteratorType iInitialNeighbours
     
    }// AddContinuumContacts

    void SphericHertzianParticle::ComputeForcesGeneral( int type_id, int damp_id, double dt_input, array_1d<double,3>& gravity)
    {
        KRATOS_TRY

        bool continuum_simulation_OPTION = true;
        bool delta_OPTION = true;

        // GETTING PARTICLE PROPERTIES
               
  //    int FailureId = mFailureId;
        int continuum_group = mContinuumGroup;

        double dt             = dt_input;
        double LocalDampRatio = mLocalDampRatio;
        double Tension        = mTension;
        double Cohesion       = mCohesion;
        double FriAngle       = mFriction;
        double Friction       = tan( FriAngle / 180.0 * M_PI);

        double radius               = mRadius;
        double critic_damp_fraction = mZeta;
        double mass                 = mMass;
        double young                = mYoung;
  //    double young_star           = mYoungStar;
        double poisson              = mPoisson;

        array_1d<double,3>& force   = GetForce();
        noalias(force)              = ZeroVector(3);
        force                       = mass * gravity;

        size_t iContactForce = 0;

        for(ParticleWeakIteratorType neighbour_iterator = SphericHertzianParticle::mNeighbours.begin();
        neighbour_iterator != SphericHertzianParticle::mNeighbours.end(); neighbour_iterator++)
        {
           
        // GETTING NEIGHBOUR PROPERTIES

        double other_radius                 = neighbour_iterator->GetRadius();
        double other_critic_damp_fraction   = neighbour_iterator->GetZeta();
        double equiv_visc_damp_ratio        = (critic_damp_fraction + other_critic_damp_fraction) / 2.0;   //M: is it correct to be a simple mean.
     // double other_mass                   = neighbour_iterator->GetMass();
        double other_young                  = neighbour_iterator->GetYoung();
    //  double other_young_star             = neighbour_iterator->GetYoungStar();
        double other_poisson                = neighbour_iterator->GetPoisson();
        double other_tension                = neighbour_iterator->GetTension();
        double other_cohesion               = neighbour_iterator->GetCohesion();  

        // CONTINUUM SIMULATING PARAMETERS:
        
        double initial_delta = 0.0;
        double CTension = 0.0;
        double CCohesion = 0.0;
                        
        if (continuum_simulation_OPTION && continuum_group && mContactFailureId[iContactForce]==0)
        {
            CTension    = 2* Tension * other_tension / (Tension + other_tension);
            CCohesion   = 2* Cohesion * other_cohesion / (Cohesion + other_cohesion);
        }
       
        if(delta_OPTION) {
               
            neighbour_iterator->GetPointerToCenterNode()->Id();
             initial_delta = mContactInitialDelta[iContactForce];
     
        }
       
        // BASIC CALCULATIONS

        array_1d<double,3> other_to_me_vect = this->GetPosition() - neighbour_iterator->GetPosition();
        double distance                     = sqrt(other_to_me_vect[0] * other_to_me_vect[0] +
                                              other_to_me_vect[1] * other_to_me_vect[1] +
                                              other_to_me_vect[2] * other_to_me_vect[2]);
        
        double radius_sum                   = radius + other_radius;
        
        double indentation                  = radius_sum - distance - initial_delta; //M: Here, Initial_delta is expected to be positive if it is embeding and negative if it's separation.
        
        double equiv_radius     = 2* radius * other_radius / (radius + other_radius);
        double equiv_area       = M_PI * equiv_radius * equiv_radius;
        double equiv_poisson    = 2* poisson * other_poisson / (poisson + other_poisson);
        double equiv_young      = 2 * young * other_young / (young + other_young);
   //   double equiv_young_star = young_star * other_young_star / (young_star + other_young_star);
        double kn               = M_PI * 0.5 * equiv_young * equiv_radius; //M: CANET FORMULA
        double ks               = kn / (2.0 * (1.0 + equiv_poisson));
                
        // FORMING LOCAL CORDINATES

        double NormalDir[3]           = {0.0};
        double LocalCoordSystem[3][3] = {{0.0}, {0.0}, {0.0}};
        NormalDir[0] = other_to_me_vect[0];   // M. this way the compresion is positive.
        NormalDir[1] = other_to_me_vect[1];
        NormalDir[2] = other_to_me_vect[2];
        ComputeContactLocalCoordSystem(NormalDir, LocalCoordSystem);

        // VELOCITIES AND DISPLACEMENTS

            array_1d<double, 3 > vel = GetVelocity();
            array_1d<double, 3 > other_to_me_vel = neighbour_iterator->GetVelocity();

            double DeltDisp[3] = {0.0};
            double DeltVel [3] = {0.0};

            DeltVel[0] = (vel[0] - other_to_me_vel[0]);
            DeltVel[1] = (vel[1] - other_to_me_vel[1]);
            DeltVel[2] = (vel[2] - other_to_me_vel[2]);

            //DeltDisp in global cordinates

            DeltDisp[0] = DeltVel[0] * dt;
            DeltDisp[1] = DeltVel[1] * dt;
            DeltDisp[2] = DeltVel[2] * dt;

            double LocalDeltDisp[3] = {0.0};

            double LocalContactForce[3]  = {0.0};
            double GlobalContactForce[3] = {0.0};

            GlobalContactForce[0] = mContactForces[iContactForce][0];
            GlobalContactForce[1] = mContactForces[iContactForce][1];
            GlobalContactForce[2] = mContactForces[iContactForce][2];

            VectorGlobal2Local(LocalCoordSystem, DeltDisp, LocalDeltDisp);
            VectorGlobal2Local(LocalCoordSystem, GlobalContactForce, LocalContactForce);

     // FORCES
    
            //M: aqui estem guardant la força total local, que es la suma de la total que teniem més l'increment, es el metode incremental.
            //mcontact forces esta inicialitzat per mitjà de GetContactForces en el NeighboursCalculator.
            //aqui em falta lu de pi i tal.... oi???  
            //M: these are arrays of vectors of 3 components.          

            if ( (indentation > 0.0) || (mContactFailureId[iContactForce] == 0) )  // This conditions take in acount the fact that the particles must remember their initial delta's between initial neighbours.
            //M: te una mica de tela entendre bé aquestes condicions pero crec que axi es compleixien tots el casos.
            {
                LocalContactForce[0] += - ks * LocalDeltDisp[0];  // 0: first tangential
                LocalContactForce[1] += - ks * LocalDeltDisp[1];  // 1: second tangential
                LocalContactForce[2] += - kn * LocalDeltDisp[2];  // 2: normal force
            }
           
            //ABSOLUTE METHOD FOR NORMAL FORCE (Allows non-linearity)

            if(type_id == 2)  //M: estem sobreeescribim els localContact que haviem fet, si es tipus 1 ja no ho reescribim.
                //el feng ha pensat lu del type id per si volem incremental o absolut.   1--- incremental; 2 --- absolut i amb el cas hertzià
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

            // per guardar els maxims i minims.


            if(mTestVariable == 0){

               maxForce = 0;
               minForce = 0;
               mTestVariable = 1;

            }

            if (LocalContactForce[2]>maxForce){
            maxForce = LocalContactForce[2];
            }
            if(LocalContactForce[2]<minForce){
            minForce = LocalContactForce[2];
            }
            


                 // TENSION FAILURE
                if (-LocalContactForce[2] > (CTension * equiv_area))  //M:si la tensio supera el limit es seteja tot a zero.
                {
                    if( (mTestVariable==1) && (LocalContactForce[2]<-50) ){
                    mTestVariable=LocalContactForce[2];

                    
                    }
                    LocalContactForce[0]  = 0.0;
                    LocalContactForce[1]  = 0.0;
                    LocalContactForce[2]  = 0.0;

                    mContactFailureId[iContactForce] = 3;  //tensile failure case.
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
                        LocalContactForce[0] = ShearForceMax / ShearForceNow * LocalContactForce[0];
                        LocalContactForce[1] = ShearForceMax / ShearForceNow * LocalContactForce[1];

                        mContactFailureId[iContactForce] = 4; // Shear failure case.
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

            mContactForces[iContactForce][0] = GlobalContactForce[0];
            mContactForces[iContactForce][1] = GlobalContactForce[1];
            mContactForces[iContactForce][2] = GlobalContactForce[2];

            iContactForce++;

        }//for each neaighbour


      //AFTER FORCES CALCULATIONS

        // DETECTING A COHESIVE PARTICLE THAT HAS BEEN COMPLETELY DETACHED.

        double num_of_detached_neighbours=0;

        for (unsigned int i=0; i< mContactFailureId.size(); i++ ) {
            if   (mContactFailureId[i] != 0)
                {
                num_of_detached_neighbours++;
                }
        }

        // LOCAL DAMPING OPTION FOR THE UNBALANCED FORCES (IN GLOBAL CORDINATES).

        if (damp_id == 1)
        {
            for (int iDof = 0; iDof < 3; iDof++)
            {
                if (GetVelocity()[iDof] > 0.0)
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
            mTension = 0.0;
            mCohesion = 0.0;

            GetTension() = mTension;
            GetCohesion() = mCohesion;
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

        if (mFailureId != 1)  // for mFailureId == 1 there's no failure to represent, the particle is not a continuum-simulating one or has been completelly detached already.
        {

            int tempType[5] = {0,0,0,0,0};
            for(size_t itype = 0; itype < mContactFailureId.size(); itype++)
            {
                if(mContactFailureId[itype] == 0)
                {
                    tempType[0]++;
                }
                else if(mContactFailureId[itype] == 1)
                {
                    tempType[1]++;
                }

                // mContactFailureId[itype] == 2 intentionally skipped!

                else if(mContactFailureId[itype] == 3)
                {
                    tempType[3]++;
                }
                else if(mContactFailureId[itype] == 4)
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

        GetFailureId() = mFailureId;

        }// if (mFailureId != 1)

        //M: si mFailureId = 1 already I think we dont need to do GetFailureId() = mFailureId; becouse the h declatation will be applyied for new objects of the class.
     KRATOS_CATCH("")

    }


    void SphericHertzianParticle::norm(double Vector[3])
    {
        double distance = sqrt(Vector[0] * Vector[0] + Vector[1] * Vector[1] + Vector[2] * Vector[2]);
        Vector[0] = Vector[0] / distance;
        Vector[1] = Vector[1] / distance;
        Vector[2] = Vector[2] / distance;
    }

    void SphericHertzianParticle::VectorGlobal2Local(double LocalCoordSystem[3][3], double GlobalVector[3], double LocalVector[3])
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

    void SphericHertzianParticle::VectorLocal2Global(double LocalCoordSystem[3][3], double LocalVector[3], double GlobalVector[3])
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

    double SphericHertzianParticle::DotProduct(double Vector1[3], double Vector2[3])
    {
        return Vector1[0] * Vector2[0] + Vector1[1] * Vector2[1] + Vector1[2] * Vector2[2];
    }

    void SphericHertzianParticle::CrossProduct(double u[3], double v[3], double ReturnVector[3])
    {
    	ReturnVector[0] = u[1]*v[2] - u[2]*v[1];
	ReturnVector[1] = v[0]*u[2] - u[0]*v[2];
	ReturnVector[2] = u[0]*v[1] - u[1]*v[0];
    }

    void SphericHertzianParticle::ComputeContactLocalCoordSystem(double NormalDirection[3], double LocalCoordSystem[3][3])
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



