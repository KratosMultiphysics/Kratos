//
//   Project Name:        Kratos
//   Last Modified by:    $Author: G.Casas $
//   Date:                $Date: 2008-10-13 08:56:42 $
//   Revision:            $Revision: 1.5 $
//
//
//README::::look to the key word "VERSION" if you want to find all the points where you have to change something so that you can pass from a kdtree to a bin data search structure;

#if !defined(KRATOS_EXPLICIT_SOLVER_WITH_ROTATION)
#define  KRATOS_EXPLICIT_SOLVER_WITH_ROTATION

// /* External includes */

// System includes

// Project includes
#include "custom_utilities/circular_particle.h"
#include "custom_utilities/spheric_particle.h"
#include "custom_utilities/circular_particle_hertzian.h"
#include "custom_utilities/spheric_particle_hertzian.h"
#include "custom_utilities/spheric_rotating_particle.h"
#include "custom_utilities/spheric_rotating_particle_hertzian.h"
//Database includes

namespace Kratos {

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
 */
template< std::size_t TDim, class TParticle >
class Explicit_Solver_With_Rotation{
    ///@}

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}

public:

    ///@name Type Definitions
    typedef TParticle ParticleType;
    typedef typename TParticle::Pointer ParticlePointerType;
    typedef typename TParticle::ParticleWeakVectorType ParticleWeakVectorType;
    typedef typename TParticle::ParticleWeakIteratorType ParticleWeakIteratorType;
    typedef typename TParticle::DistanceVectorType DistanceVectorType;
    typedef typename TParticle::DistanceIteratorType DistanceIteratorType;
    typedef ModelPart::NodesContainerType::iterator PointPointerIterator;
    typedef std::vector<ParticleType> ParticleVectorType;
    typedef typename std::vector<ParticleType>::iterator ParticleIteratorType;
    typedef typename std::vector<typename ParticleType::Pointer> ParticlePointerVectorType;
    typedef typename std::vector<typename ParticleType::Pointer>::iterator ParticlePointerIteratorType;
    typedef std::vector<array_1d<double, 3 > > ComponentVectorType;
    typedef std::vector<array_1d<double, 3 > >::iterator ComponentIteratorType;

    ///@{

    /// Pointer definition of Explicit_Solver
    KRATOS_CLASS_POINTER_DEFINITION(Explicit_Solver_With_Rotation);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.

    Explicit_Solver_With_Rotation(int solver_id, double radius, double tol, ModelPart& model_part){
        mModelPart = &model_part;
        mRadiusSearch = radius;
        mProximityTol = tol;
        mMinRadius = radius;
        mMaxRadius = 0.0;
        mMinMass = 2.0 * model_part.NodesBegin()->GetSolutionStepValue(RADIUS);
        mNumberOfParticles = model_part.Nodes().size();
        mSolverId = solver_id;
        for (PointPointerIterator inode = model_part.NodesBegin(); inode != model_part.NodesEnd(); inode++){
            Node < 3 > ::Pointer center_pointer = *(inode.base());
            ParticlePointerType particle_pointer = ParticlePointerType(new ParticleType(tol, center_pointer));
            if (mMinRadius > particle_pointer->GetRadius()){mMinRadius = particle_pointer->GetRadius();}
            if (mMaxRadius < particle_pointer->GetRadius()){mMaxRadius = particle_pointer->GetRadius();}
            if (mMinMass > particle_pointer->GetMass()){mMinMass = particle_pointer->GetMass();}
            mListOfParticlePointers.push_back(particle_pointer);
            }
        }

    /// Destructor.

    virtual ~Explicit_Solver_With_Rotation(){}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void Search_Neighbours(){
        Neighbours_Calculator<TDim, ParticleType, ParticlePointerType, ParticleVectorType, ParticleWeakVectorType, ParticlePointerVectorType,
        ParticleWeakIteratorType, ParticleIteratorType, ParticlePointerIteratorType, DistanceVectorType, DistanceIteratorType>::
        Search_Neighbours(mListOfParticlePointers, *mModelPart, mRadiusSearch, mProximityTol);
    }

    void Evolve_Displacements(){
        KRATOS_TRY
        int size = mListOfParticlePointers.size();
        //#pragma omp parallel for
        for (int i = 0; i < size; i++){
            ParticlePointerIteratorType particle_it = mListOfParticlePointers.begin() + i;
            (*particle_it)->ComputeNewTangentialDisplacements();
            }
        KRATOS_CATCH("")
    }
    void Update_Contacting_Neighbours(){
        KRATOS_TRY
        int size = mListOfParticlePointers.size();
        //#pragma omp parallel for
        for (int i = 0; i < size; i++){
            ParticlePointerIteratorType particle_it = mListOfParticlePointers.begin() + i;
            (*particle_it)->UpdateContactsList();
            }
        KRATOS_CATCH("")
    }
    void Calculate_Forces(double delta_t, array_1d<double, 3 > gravity){
        KRATOS_TRY
        int size = mListOfParticlePointers.size();
        //#pragma omp parallel for
        for (int i = 0; i < size; i++){
            ParticlePointerIteratorType particle_it = mListOfParticlePointers.begin() + i;
            (*particle_it)->ComputeForcesOnCenterNode(delta_t, gravity);
            }
        KRATOS_CATCH("")
        }

    void Evolve_Motion(double delta_t, array_1d<double, 3 > gravity){
        int size = mListOfParticlePointers.size();
//****************************************************************************************************************************************************//
//                                                            F O R W A R D     E U L E R
//****************************************************************************************************************************************************//
        if (mSolverId == 1){
            double aux;
            for (ParticlePointerIteratorType particle_it = mListOfParticlePointers.begin();
                particle_it != mListOfParticlePointers.end(); ++particle_it){
                if (!(*particle_it)->GetPointerToCenterNode()->IsFixed(VELOCITY_X)){
                    array_1d<double, 3 > & vel = (*particle_it)->GetVelocity();
                    array_1d<double, 3 > & ang_vel = (*particle_it)->GetAngularVelocity();
                    array_1d<double, 3 > & displ = (*particle_it)->GetDisplacement();
                    array_1d<double, 3 > & rotation = (*particle_it)->GetRotation();
                    array_1d<double, 3 > & coor = (*particle_it)->GetPosition();
                    array_1d<double, 3 > & initial_coor = (*particle_it)->GetInitialPosition();
                    const array_1d<double, 3 > & force = (*particle_it)->GetForce();
                    const array_1d<double, 3 > & moment = (*particle_it)->GetMoment();
                    double aux_1 = delta_t / (*particle_it)->GetMass();
                    double aux_2 = delta_t / (*particle_it)->GetInertia();
                    noalias(vel) += aux_1 * force;
                    noalias(ang_vel) += aux_2 * moment;
                    //Evolution of position (u(n+1) = u(n) + v(n+0.5)*delta_t):
                    noalias(displ) += delta_t * vel;
                    noalias(coor) = initial_coor + displ;
                    noalias(rotation) += delta_t * ang_vel;
                    }
                }
                Evolve_Displacements();
            }//Forward Euler

//****************************************************************************************************************************************************//
//                                                             R U N G E    K U T T A    4
//****************************************************************************************************************************************************//
        else if(mSolverId == 2){
            double half_delta_t = 0.5 * delta_t;
            array_1d<double, 3 > aux_1, aux_2, vel_copy, displ_copy;
            ComponentVectorType vel_old, displ_old, kf, kf1, kf2, kf3, kf4, kv, kv1, kv2, kv3, kv4;
            ComponentVectorType ang_vel_old, rotation_old, lf, lf1, lf2, lf3, lf4, lv, lv1, lv2, lv3, lv4;
//***********************************       K1      ***************************************************************************************************//
            for (int i = 0; i < size; i++){
                ParticlePointerIteratorType particle_it = mListOfParticlePointers.begin() + i;
                array_1d<double, 3 > & vel = (*particle_it)->GetVelocity();
                array_1d<double, 3 > & ang_vel = (*particle_it)->GetAngularVelocity();
                array_1d<double, 3 > & displ = (*particle_it)->GetDisplacement();
                array_1d<double, 3 > & rotation = (*particle_it)->GetRotation();
                array_1d<double, 3 > & force = (*particle_it)->GetForce();
                array_1d<double, 3 > & moment = (*particle_it)->GetMoment();
                vel_old.push_back(vel);
                ang_vel_old.push_back(vel);
                displ_old.push_back(displ);
                rotation_old.push_back(displ);
                aux_1 = delta_t * force / (*particle_it)->GetMass();
                aux_2 = delta_t * moment / (*particle_it)->GetInertia();
                kv1.push_back(vel * delta_t);
                kf1.push_back(aux_1);
                lv1.push_back(ang_vel * delta_t);
                lf1.push_back(aux_2);
                if (!(*particle_it)->GetPointerToCenterNode()->IsFixed(VELOCITY_X)){
                    noalias(displ) = displ_old[i] + half_delta_t * vel;
                    noalias(vel) = vel_old[i] + 0.5 * aux_1;
                    noalias(rotation) = rotation_old[i] + half_delta_t * ang_vel;
                    noalias(ang_vel) = ang_vel_old[i] + 0.5 * aux_2;
                    }
                }
            Evolve_Displacements();
            Update_Contacting_Neighbours();
            Calculate_Forces(delta_t, gravity);
//***********************************       K2      ***************************************************************************************************//
            for (int i = 0; i < size; i++){
                ParticlePointerIteratorType particle_it = mListOfParticlePointers.begin() + i;
                    array_1d<double, 3 > & vel = (*particle_it)->GetVelocity();
                    array_1d<double, 3 > & ang_vel = (*particle_it)->GetAngularVelocity();
                    array_1d<double, 3 > & displ = (*particle_it)->GetDisplacement();
                    array_1d<double, 3 > & force = (*particle_it)->GetForce();
                    array_1d<double, 3 > & rotation = (*particle_it)->GetRotation();
                    array_1d<double, 3 > & moment = (*particle_it)->GetMoment();
                    aux_1 = delta_t * force / (*particle_it)->GetMass();
                    aux_2 = delta_t * moment / (*particle_it)->GetInertia();
                    kv2.push_back(vel * delta_t);
                    kf2.push_back(aux_1);
                    lv2.push_back(ang_vel * delta_t);
                    lf2.push_back(aux_2);
                    if (!(*particle_it)->GetPointerToCenterNode()->IsFixed(VELOCITY_X)){
                        noalias(displ) = displ_old[i] + half_delta_t * vel;
                        noalias(vel) = vel_old[i] + 0.5 * aux_1;
                        noalias(rotation) = rotation_old[i] + half_delta_t * ang_vel;
                        noalias(ang_vel) = ang_vel_old[i] + 0.5 * aux_2;
                        }
                }
            Evolve_Displacements();
            Update_Contacting_Neighbours();
            Calculate_Forces(delta_t, gravity);
//***********************************       K3      ***************************************************************************************************//
            for (int i = 0; i < size; i++){
                ParticlePointerIteratorType particle_it = mListOfParticlePointers.begin() + i;
                array_1d<double, 3 > & vel = (*particle_it)->GetVelocity();
                array_1d<double, 3 > & ang_vel = (*particle_it)->GetAngularVelocity();
                array_1d<double, 3 > & displ = (*particle_it)->GetDisplacement();
                array_1d<double, 3 > & force = (*particle_it)->GetForce();
                array_1d<double, 3 > & rotation = (*particle_it)->GetRotation();
                array_1d<double, 3 > & moment = (*particle_it)->GetMoment();
                aux_1 = delta_t * force / (*particle_it)->GetMass();
                aux_2 = delta_t * moment / (*particle_it)->GetInertia();
                kv3.push_back(vel * delta_t);
                kf3.push_back(aux_1);
                lv3.push_back(ang_vel * delta_t);
                lf3.push_back(aux_2);
                if (!(*particle_it)->GetPointerToCenterNode()->IsFixed(VELOCITY_X)){
                    noalias(displ) = displ_old[i] + delta_t * vel;
                    noalias(rotation) = rotation_old[i] + delta_t * ang_vel;
                    noalias(ang_vel) = ang_vel_old[i] + aux_2;
                        }
                }
            Evolve_Displacements();
            Update_Contacting_Neighbours();
            Calculate_Forces(delta_t, gravity);
//***********************************       K4      ***************************************************************************************************//
           for (int i = 0; i < size; i++){
                ParticlePointerIteratorType particle_it = mListOfParticlePointers.begin() + i;
                array_1d<double, 3 > & vel = (*particle_it)->GetVelocity();
                array_1d<double, 3 > & ang_vel = (*particle_it)->GetAngularVelocity();
                array_1d<double, 3 > & force = (*particle_it)->GetForce();
                array_1d<double, 3 > & moment = (*particle_it)->GetMoment();
                aux_1 = delta_t * force / (*particle_it)->GetMass();
                aux_2 = delta_t * moment / (*particle_it)->GetInertia();
                kv4.push_back(vel * delta_t);
                kf4.push_back(aux_1);
                lv4.push_back(ang_vel * delta_t);
                lf4.push_back(aux_2);
                }
//***********************************       MEAN K AND EVOLUTION      *********************************************************************************//

            for (int i = 0; i < size; i++){
                ParticlePointerIteratorType particle_it = mListOfParticlePointers.begin() + i;
                kf.push_back(kf1[i] + 2.0 * (kf2[i] + kf3[i]) + kf4[i]);
                kv.push_back(kv1[i] + 2.0 * (kv2[i] + kv3[i]) + kv4[i]);
                lf.push_back(lf1[i] + 2.0 * (lf2[i] + lf3[i]) + lf4[i]);
                lv.push_back(lv1[i] + 2.0 * (lv2[i] + lv3[i]) + lv4[i]);
                if (!(*particle_it)->GetPointerToCenterNode()->IsFixed(VELOCITY_X)){
                    array_1d<double, 3 > & vel = (*particle_it)->GetVelocity();
                    array_1d<double, 3 > & displ = (*particle_it)->GetDisplacement();
                    array_1d<double, 3 > & ang_vel = (*particle_it)->GetAngularVelocity();
                    array_1d<double, 3 > & rotation = (*particle_it)->GetRotation();
                    array_1d<double, 3 > & coor = (*particle_it)->GetPosition();
                    array_1d<double, 3 > & initial_coor = (*particle_it)->GetInitialPosition();
                    noalias(displ) = displ_old[i] + 0.166666666666666667 * kv[i];
                    noalias(coor) = initial_coor + displ;
                    noalias(vel) = vel_old[i] + 0.166666666666666667 * kf[i];
                    noalias(rotation) = rotation_old[i] + 0.166666666666666667 * lv[i];
                    noalias(ang_vel) = ang_vel_old[i] + 0.166666666666666667 * lf[i];
                    }
                }
            Evolve_Displacements();
    }//Runge Kutta 4

//****************************************************************************************************************************************************//
//                                                              M I D - P O I N T    R U L E
//****************************************************************************************************************************************************//
            else if (mSolverId == 3) {
                double half_delta_t = 0.5 * delta_t;
                array_1d<double, 3 > aux_1, vel_copy, displ_copy;
                array_1d<double, 3 > aux_2, ang_vel_copy, rotation_copy;
                ComponentVectorType vel_old, displ_new, kf, kv;
                ComponentVectorType ang_vel_old, rotation_new, lf, lv;
                //***********************************       K      ***************************************************************************************************//
                for (int i = 0; i < size; i++) {
                    ParticlePointerIteratorType particle_it = mListOfParticlePointers.begin() + i;
                    array_1d<double, 3 > & vel = (*particle_it)->GetVelocity();
                    array_1d<double, 3 > & ang_vel = (*particle_it)->GetAngularVelocity();
                    array_1d<double, 3 > & displ = (*particle_it)->GetDisplacement();
                    array_1d<double, 3 > & rotation = (*particle_it)->GetRotation();
                    array_1d<double, 3 > & force = (*particle_it)->GetForce();
                    array_1d<double, 3 > & moment = (*particle_it)->GetMoment();
                    aux_1 = half_delta_t * force / (*particle_it)->GetMass();
                    aux_2 = delta_t * moment / (*particle_it)->GetInertia();
                    vel_old.push_back(vel);
                    ang_vel_old.push_back(ang_vel);
                    displ_new.push_back(displ + delta_t * vel * (1 + half_delta_t));
                    rotation_new.push_back(rotation + delta_t * ang_vel * (1 + half_delta_t));
                    if (!(*particle_it)->GetPointerToCenterNode()->IsFixed(VELOCITY_X)) {
                        noalias(vel) += aux_1;
                        noalias(displ) += half_delta_t * vel;
                        noalias(ang_vel) += aux_2;
                        noalias(rotation) += half_delta_t * ang_vel;
                    }
                }
                Evolve_Displacements();
                Update_Contacting_Neighbours();
                Calculate_Forces(delta_t, gravity);
                //***********************************       EVOLUTION      ********************************************************************************************//
                for (int i = 0; i < size; i++) {
                    ParticlePointerIteratorType particle_it = mListOfParticlePointers.begin() + i;
                    if (!(*particle_it)->GetPointerToCenterNode()->IsFixed(VELOCITY_X)) {
                        array_1d<double, 3 > & vel = (*particle_it)->GetVelocity();
                        array_1d<double, 3 > & displ = (*particle_it)->GetDisplacement();
                        array_1d<double, 3 > & coor = (*particle_it)->GetPosition();
                        array_1d<double, 3 > & initial_coor = (*particle_it)->GetInitialPosition();
                        array_1d<double, 3 > & force = (*particle_it)->GetForce();
                        array_1d<double, 3 > & rotation = (*particle_it)->GetRotation();
                        array_1d<double, 3 > & ang_vel = (*particle_it)->GetAngularVelocity();
                        array_1d<double, 3 > & moment = (*particle_it)->GetMoment();
                        noalias(displ) = displ_new[i];
                        noalias(coor) = initial_coor + displ;
                        noalias(vel) = vel_old[i] + delta_t * force / (*particle_it)->GetMass();
                        noalias(rotation) = rotation_new[i];
                        noalias(ang_vel) = ang_vel_old[i] + delta_t * moment / (*particle_it)->GetInertia();
                    }
                }
                Evolve_Displacements();
            }// Mid-Point
        }//EVOLVE MOTION

        double Estimate_Time_Step_Circles(double alpha) { //NO ESTUDIADO CON CUIDADO, PERO TIENE QUE SER PARECIDO A SPHERES

            double max_nat_freq, max_nat_freq_2 = 0.0;
            double max_damping_ratio = 0.0;
            double dt;
            for (ParticlePointerIteratorType particle_it = mListOfParticlePointers.begin();
                    particle_it != mListOfParticlePointers.end(); ++particle_it) {
                double k = (*particle_it)->GetStiffness();
                double r = (*particle_it)->GetRadius();
                double m = (*particle_it)->GetMass();
                double new_freq_2 = k / (r * m);
                double new_damping_ratio = (*particle_it)->GetZeta();
                double new_critic_damp = new_damping_ratio * sqrt(new_freq_2);
                if (new_freq_2 > max_nat_freq_2) {
                    max_nat_freq_2 = new_freq_2;
                }
                if (new_damping_ratio > max_damping_ratio) {
                    max_damping_ratio = new_damping_ratio;
                }
            }
            max_nat_freq_2 = 0.25 * max_nat_freq_2; //natural frequency squared for a system of two equal particles with maximum natural frequences
            max_nat_freq = sqrt(max_nat_freq_2);
            double max_nat_freq_inv = 1 / max_nat_freq;
            double max_nat_freq_2_inv = 1 / max_nat_freq_2;
            if (mSolverId == 1) {
                dt = alpha * max_nat_freq_2_inv * (sqrt(1 + max_damping_ratio * max_damping_ratio) - max_damping_ratio);
            } else {
                dt = alpha * 2 * max_nat_freq_inv * (sqrt(1 + max_damping_ratio * max_damping_ratio) - max_damping_ratio);
            }
            return (dt);
        }//ESTIMATE TIME STEP CIRCLES

        double Estimate_Time_Step_Spheres(double alpha) {
            //EVITAR AVANZAR DEMASIADO EN UN PASO DE TIEMPO (PUEDE INTERESAR DESCOMENTAR...))
            //        array_1d<double, 3 > vel = (*mListOfParticlePointers.begin())->GetVelocity();
            //        array_1d<double, 3 > force = (*mListOfParticlePointers.begin())->GetForce();
            //        double vel_mod_2;
            //        double force_mod_2;
            //        double max_vel, max_acc, max_vel_2, max_force_2;
            //        double dt;
            //        for (ParticlePointerIteratorType particle_it = mListOfParticlePointers.begin();
            //            particle_it != mListOfParticlePointers.end(); ++particle_it) {
            //            vel = (*particle_it)->GetVelocity();
            //            force = (*particle_it)->GetForce();
            //            vel_mod_2 = vel(0) * vel(0) + vel(1) * vel(1);
            //            force_mod_2 = force(0) * force(0) + force(1) * force(1);
            //            if (max_vel_2 < vel_mod_2){max_vel_2 = vel_mod_2;}
            //            if (max_force_2 < force_mod_2){max_force_2 = force_mod_2;}
            //            }
            //            max_acc = 2 * sqrt(max_force_2) / mMinMass;
            //            max_vel = sqrt(max_vel_2);
            //            if (max_acc > max_vel){
            //                dt = (sqrt(max_vel_2 + max_acc * alpha * mMinRadius * alpha) - max_vel) / max_acc;
            //                }
            //            else{
            //                dt = mMinRadius * alpha / max_vel;
            //                }
            //            return(dt);
            double max_nat_freq, max_nat_freq_2 = 0.0;
            double max_damping_ratio = 0.0;
            double dt;
            double radius_quot = M_PI * (mMaxRadius + mMinRadius) / mMinRadius;
            for (ParticlePointerIteratorType particle_it = mListOfParticlePointers.begin();
                    particle_it != mListOfParticlePointers.end(); ++particle_it) {
                double k = (*particle_it)->GetNormalStiffness();
                double r = (*particle_it)->GetRadius();
                double m = (*particle_it)->GetMass();
                double new_freq_2 = k / (r * m);
                double new_damping_ratio = (*particle_it)->GetZeta();
                double new_critic_damp = new_damping_ratio * sqrt(new_freq_2);
                if (new_freq_2 > max_nat_freq_2) {
                    max_nat_freq_2 = new_freq_2;
                }
                if (new_damping_ratio > max_damping_ratio) {
                    max_damping_ratio = new_damping_ratio;
                }
            }
            max_nat_freq_2 = 0.25 * max_nat_freq_2; //natural frequency squared for a system of two equal particles with maximum natural frequences
            max_nat_freq = sqrt(max_nat_freq_2);
            double max_nat_freq_inv = 1 / max_nat_freq;
            double max_nat_freq_2_inv = 1 / max_nat_freq_2;
            if (mSolverId == 1) {
                dt = alpha * max_nat_freq_2_inv * (sqrt(1 + max_damping_ratio * max_damping_ratio) - max_damping_ratio);
            } else {
                dt = alpha * 2 * max_nat_freq_inv * (sqrt(1 + max_damping_ratio * max_damping_ratio) - max_damping_ratio);
            }
            return (dt);
        }//ESTIMATE TIME STEP SPHERES

        double Estimate_Time_Step_Hertzian_Circles(double alpha) { //NO ESTUDIADO CON CUIDADO, PERO TIENE QUE SER PARECIDO A SPHERES

            double max_freq, max_freq_2 = 0.0;
            double damp_ratio = 0.0;
            double dt;
            //        double radius_quot = pi * (mMaxRadius + mMinRadius) / mMinRadius;
            double radius_quot = 0.01;
            for (ParticlePointerIteratorType particle_it = mListOfParticlePointers.begin();
                    particle_it != mListOfParticlePointers.end(); ++particle_it) {
                double young = (*particle_it)->GetYoungStar() * radius_quot;
                double r = (*particle_it)->GetRadius();
                double m = (*particle_it)->GetMass();
                double new_freq_2 = young / (r * m);
                double new_damp_ratio = (*particle_it)->GetZeta();
                if (new_freq_2 > max_freq_2) {
                    max_freq_2 = new_freq_2;
                }
                if (new_damp_ratio > damp_ratio) {
                    damp_ratio = new_damp_ratio;
                }
            }
            max_freq = 0.5 * sqrt(max_freq_2); //natural frequency for a system of two equal particles
            double max_freq_inv = 1 / max_freq;
            dt = alpha * 2 * max_freq_inv * (sqrt(1 + damp_ratio * damp_ratio) - damp_ratio);
            return (dt);
        }//ESTIMATE TIME STEP HERTZIAN CIRCLES

        double Estimate_Time_Step_Hertzian_Spheres(double alpha) {

            double max_freq, max_freq_2 = 0.0;
            double damp_ratio = 0.0;
            double dt;
            //        double radius_quot = pi * (mMaxRadius + mMinRadius) / mMinRadius;
            double radius_quot = 0.01;
            for (ParticlePointerIteratorType particle_it = mListOfParticlePointers.begin();
                    particle_it != mListOfParticlePointers.end(); ++particle_it) {
                double young = (*particle_it)->GetYoungStar() * radius_quot;
                double r = (*particle_it)->GetRadius();
                double m = (*particle_it)->GetMass();
                double new_freq_2 = young / (r * m);
                double new_damp_ratio = (*particle_it)->GetZeta();
                if (new_freq_2 > max_freq_2) {
                    max_freq_2 = new_freq_2;
                }
                if (new_damp_ratio > damp_ratio) {
                    damp_ratio = new_damp_ratio;
                }
            }
            max_freq = 0.5 * sqrt(max_freq_2); //natural frequency for a system of two equal particles
            double max_freq_inv = 1 / max_freq;
            dt = alpha * 2 * max_freq_inv * (sqrt(1 + damp_ratio * damp_ratio) - damp_ratio);
            return (dt);
        }//ESTIMATE TIME STEP HERTZIAN SPHERES


        ///@}
        ///@name Access
        ///@{

        ParticlePointerVectorType& GetListOfParticlePointers() {
            return (mListOfParticlePointers);
        }

        ///@}
        ///@name Inquiry
        ///@{

        ///@}
        ///@name Input and output
        ///@{

        /// Turn back information as a stemplate<class T, std::size_t dim> tring.

        virtual std::string Info() const {
            return "";
        }

        /// Print information about this object.

        virtual void PrintInfo(std::ostream& rOStream) const {
        }

        /// Print object's data.

        virtual void PrintData(std::ostream& rOStream) const {
        }

        ///@}
        ///@name Friends
        ///@{

        ///@}

    protected:
        ///@name Protected static Member rVariables
        ///@{

        ///@}
        ///@name Protected member rVariables
        ///@{ template<class T, std::size_t dim>

        ///@}
        ///@name Protected Operators
        ///@{

        ///@}
        ///@name Protected Operations
        ///@{

        ///@}
        ///@name Protected  Access
        ///@{

        ///@}
        ///@name Protected Inquiry
        ///@{

        ///@}
        ///@name Protected LifeCycle
        ///@{

        ///@}

    private:
        ///@name Static Member rVariables
        ///@{

        ///@}
        ///@name Member rVariables
        ///@{

        int mNumberOfParticles;
        int mSolverId;
        double mRadiusSearch;
        double mProximityTol;
        double mMinRadius;
        double mMaxRadius;
        double mMinMass;
        ParticleVectorType mListOfParticles;
        ParticlePointerVectorType mListOfParticlePointers;
        ModelPart * mModelPart;

        ///@}
        ///@name Private Operators
        ///@{

        ///@}
        ///@name Private Operations
        ///@{

        inline void Clear(ModelPart::NodesContainerType::iterator node_it, int step_data_size) {
            unsigned int buffer_size = node_it->GetBufferSize();
            for (unsigned int step = 0; step < buffer_size; step++) {
                //getting the data of the solution step
                double* step_data = (node_it)->SolutionStepData().Data(step);

                //copying this data in the position of the vector we are interested in
                for (int j = 0; j < step_data_size; j++) {
                    step_data[j] = 0.0;
                }
            }
        }

        inline void ClearVariables(ModelPart::NodesContainerType::iterator node_it, Variable<array_1d<double, 3 > >& rVariable) {
            array_1d<double, 3 > & Aux_var = node_it->FastGetSolutionStepValue(rVariable, 0);
            noalias(Aux_var) = ZeroVector(3);
        }

        inline void ClearVariables(ModelPart::NodesContainerType::iterator node_it, Variable<double>& rVariable) {
            double& Aux_var = node_it->FastGetSolutionStepValue(rVariable, 0);
            Aux_var = 0.0;
        }

        ///@}
        ///@name Private  Access
        ///@{

        ///@}
        ///@name Private Inquiry
        ///@{

        ///@}
        ///@name Un accessible methods
        ///@{

        /// Assignment operator.
        Explicit_Solver_With_Rotation & operator=(Explicit_Solver_With_Rotation const& rOther);

        ///@}

    }; // Class Explicit_Solver_With_Rotation

    ///@}

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// output stream function
    // 	template<std::size_t TDim>
    // 	inline std::ostream& operator << (std::ostream& rOStream)
    // 	{
    // 		rThis.PrintInfo(rOStream);
    // 		rOStream << std::endl;
    // 		rThis.PrintData(rOStream);
    //
    // 		return rOStream;
    // 	}
    ///@}

} // namespace Kratos.
#endif // KRATOS_EXPLICIT_SOLVER_WITH_ROTATION  defined


