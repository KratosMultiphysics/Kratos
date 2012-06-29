/*
==============================================================================
KratosStructuralApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
 */

//
//   Project Name:        Kratos
//   Last modified by:    $Author: virginia $
//   Date:                $Date: 2009-01-23 14:39:59 $
//   Revision:            $Revision: 1.27 $
//
//


// System includes

// External includes


// Project includes
#include "custom_elements/DEM_FEM_Particle.h"
#include "includes/define.h"
#include "custom_utilities/GeometryFunctions.h"
#include "includes/constitutive_law.h"
#include "DEM_application.h"

//#include <omp.h>

namespace Kratos
{
    using namespace GeometryFunctions;

    DEM_FEM_Particle::DEM_FEM_Particle( IndexType NewId, GeometryType::Pointer pGeometry )
            : SphericParticle( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
    }

    //************************************************************************************
    //************************************************************************************

    DEM_FEM_Particle::DEM_FEM_Particle( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
            : SphericParticle( NewId, pGeometry, pProperties )
    {
        //         const unsigned int dim = GetGeometry().WorkingSpaceDimension();

    }

    Element::Pointer DEM_FEM_Particle::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
    {
        return Element::Pointer( new DEM_FEM_Particle( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
    }

    DEM_FEM_Particle::~DEM_FEM_Particle()
    {
    }


    void DEM_FEM_Particle::Initialize()
    {

        SphericParticle::Initialize();

        //Extra functions;

    }

     void DEM_FEM_Particle::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& CurrentProcessInfo )
    {
        //////Cfeng:Important calculate the particle contact force, add code here

        SphericParticle::CalculateRightHandSide(rRightHandSideVector, CurrentProcessInfo);


        ComputeParticleBlockContactForce(CurrentProcessInfo);

    }




    ////************************************************************************************
    ////************************************************************************************

    void DEM_FEM_Particle::InitializeSolutionStep( ProcessInfo& CurrentProcessInfo )
    {
        SphericParticle::InitializeSolutionStep(CurrentProcessInfo);

        //KRATOS_WATCH(mTimeStep);
    }


    ////************************************************************************************
    ////************************************************************************************

    void DEM_FEM_Particle::FinalizeSolutionStep( ProcessInfo& CurrentProcessInfo )
    {
  
    }

    void DEM_FEM_Particle::FindContactFaceOfBlockForParticle(ParticleWeakIterator rObj_2, int & RightFace, double LocalCoordSystem[3][3], double Coeff[4],double &DistPToB)
    {

        double Particle_Coord[3] = {0.0};
        Particle_Coord[0] = this->GetGeometry()(0)->Coordinates()[0];
        Particle_Coord[1] = this->GetGeometry()(0)->Coordinates()[1];
        Particle_Coord[2] = this->GetGeometry()(0)->Coordinates()[2];

        double rad = this->GetGeometry()(0)->FastGetSolutionStepValue(RADIUS);

        double Centroid[3] = {0.0};
        for(std::size_t inode = 0; inode < rObj_2->GetGeometry().size(); inode++)
        {
            Centroid[0] += rObj_2->GetGeometry()(inode)->Coordinates()[0] / (double)rObj_2->GetGeometry().size();
            Centroid[1] += rObj_2->GetGeometry()(inode)->Coordinates()[1] / (double)rObj_2->GetGeometry().size();
            Centroid[2] += rObj_2->GetGeometry()(inode)->Coordinates()[2] / (double)rObj_2->GetGeometry().size();
        }

         std::size_t dim = rObj_2->GetGeometry().WorkingSpaceDimension();
         if(dim == 2)
         {
             double Coord1[3]= {0.0};
             double Coord2[3]= {0.0};

            for (std::size_t iedge = 0; iedge < rObj_2->GetGeometry().Edges().size(); iedge++)
            {
                if( rObj_2->GetValue(IF_BOUNDARY_FACE)[iedge] == 1 )
                {
                   
                    Coord1[0] = rObj_2->GetGeometry().Edges()[iedge](0)->Coordinates()[0];
                    Coord1[1] = rObj_2->GetGeometry().Edges()[iedge](0)->Coordinates()[1];
                    Coord1[2] = rObj_2->GetGeometry().Edges()[iedge](0)->Coordinates()[2];

                    Coord2[0] = rObj_2->GetGeometry().Edges()[iedge](1)->Coordinates()[0];
                    Coord2[1] = rObj_2->GetGeometry().Edges()[iedge](1)->Coordinates()[1];
                    Coord2[2] = rObj_2->GetGeometry().Edges()[iedge](1)->Coordinates()[2];

                    double Weight[2];
                    bool If_PB_Contact = GeometryFunctions::JudgeIfThisEdgeIsContactWithParticle(Coord1, Coord2, Centroid, Particle_Coord, rad, LocalCoordSystem, Weight, DistPToB);

                    if(If_PB_Contact == true)
                    {
                        RightFace = iedge;

                        Coeff[0] = Weight[0];
                        Coeff[1] = Weight[1];

                        break;
                    }
                }

            }
         }
         else if(dim == 3)
         {
             double Coord[4][3] = { {0.0},{0.0},{0.0},{0.0} };
             for (std::size_t iface = 0; iface < rObj_2->GetGeometry().Faces().size(); iface++)
            {
                if( rObj_2->GetValue(IF_BOUNDARY_FACE)[iface] == 1 )
                {
                     ////Cfeng: Triangle
                     int FaceNodeTotal = 3;

                    Coord[0][0] = rObj_2->GetGeometry().Faces()[iface](0)->Coordinates()[0];
                    Coord[0][1] = rObj_2->GetGeometry().Faces()[iface](0)->Coordinates()[1];
                    Coord[0][2] = rObj_2->GetGeometry().Faces()[iface](0)->Coordinates()[2];

                    Coord[1][0] = rObj_2->GetGeometry().Faces()[iface](1)->Coordinates()[0];
                    Coord[1][1] = rObj_2->GetGeometry().Faces()[iface](1)->Coordinates()[1];
                    Coord[1][2] = rObj_2->GetGeometry().Faces()[iface](1)->Coordinates()[2];

                    Coord[2][0] = rObj_2->GetGeometry().Faces()[iface](2)->Coordinates()[0];
                    Coord[2][1] = rObj_2->GetGeometry().Faces()[iface](2)->Coordinates()[1];
                    Coord[2][2] = rObj_2->GetGeometry().Faces()[iface](2)->Coordinates()[2];

                    if(rObj_2->GetGeometry().Faces()[iface].size() == 4)
                    {
                        Coord[3][0] = rObj_2->GetGeometry().Faces()[iface](3)->Coordinates()[0];
                        Coord[3][1] = rObj_2->GetGeometry().Faces()[iface](3)->Coordinates()[1];
                        Coord[3][2] = rObj_2->GetGeometry().Faces()[iface](3)->Coordinates()[2];

                      ////Cfeng: Quadral
                      FaceNodeTotal = 4;
                    }

                    double Weight[4] = {0.0};

                    bool If_PB_Contact = GeometryFunctions::JudgeIfThisFaceIsContactWithParticle(FaceNodeTotal, Coord, Centroid, Particle_Coord, rad, LocalCoordSystem, Weight, DistPToB);

                    if(If_PB_Contact == true)
                    {
                        RightFace = iface;

                        Coeff[0] = Weight[0];
                        Coeff[1] = Weight[1];
                        Coeff[2] = Weight[2];
                        Coeff[3] = Weight[3];
                        break;
                    }
                }

            }

         }
    }



  



    void DEM_FEM_Particle::ComputeParticleBlockContactForce(ProcessInfo& CurrentProcessInfo)
    {
        double mTimeStep    = CurrentProcessInfo[DELTA_TIME];
        int rotation_OPTION = CurrentProcessInfo[ROTATION_OPTION]; 
        

        double Tension        = GetProperties()[PARTICLE_TENSION];
        double Cohesion       = GetProperties()[PARTICLE_COHESION];
        double FriAngle       = GetProperties()[PARTICLE_FRICTION];
        double young          = GetProperties()[YOUNG_MODULUS];
        double poisson        = GetProperties()[POISSON_RATIO];
        double Friction       = tan( FriAngle  / 180.0 * M_PI);
        double radius         = GetGeometry()(0)->FastGetSolutionStepValue(RADIUS);
        double area           = M_PI * radius * radius;
        double kn             = young * area / (2.0 * radius);
        double ks             = kn / (2.0 * (1.0 + poisson));


        array_1d<double,3>& force   = GetGeometry()(0)->FastGetSolutionStepValue(RHS);

        array_1d<double, 3 > & mRota_Moment = GetGeometry()(0)->FastGetSolutionStepValue(PARTICLE_MOMENT);
        
        Vector & mContactForces      = this->GetValue(PARTICLE_BLOCK_CONTACT_FORCE);
        Vector & mContactFailureType = this->GetValue(PARTICLE_BLOCK_CONTACT_FAILURE_ID);
        Vector & mIfInitalContact    = this->GetValue(PARTICLE_BLOCK_IF_INITIAL_CONTACT);
        ParticleWeakVector & rE      = this->GetValue(NEIGHBOUR_PARTICLE_BLOCK_ELEMENTS);


        size_t iContactForce = 0;

        for(ParticleWeakIterator ineighbour = rE.begin(); ineighbour != rE.end(); ineighbour++)
        {
            std::size_t dim = ineighbour->GetGeometry().WorkingSpaceDimension();

            ///Cfeng: contact failure type should be initialized zero
            mContactFailureType[iContactForce] = 0;

            // Not the initial contact any more, the tension and cohesion should be zero
            if(mIfInitalContact[iContactForce] == 0)
            {
                Tension  = 0.0;
                Cohesion = 0.0;
            }



            double LocalContactForce[3]  = {0.0};
            double GlobalContactForce[3] = {0.0};
            double GlobalContactForceOld[3] = {0.0};


            array_1d<double, 3 > vel = GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY);

            array_1d<double, 3 > other_to_me_vel;
            noalias(other_to_me_vel) = ZeroVector(3);


            int RightFace = -1;
            double LocalCoordSystem[3][3] = {{0.0}, {0.0}, {0.0}};
            double Coeff[4] = {0.0};
            double DistPToB = radius;

            FindContactFaceOfBlockForParticle(ineighbour, RightFace, LocalCoordSystem, Coeff, DistPToB);


            if(RightFace != -1)
            {
               if(dim == 2)
               {
                     other_to_me_vel += ineighbour->GetGeometry().Edges()[RightFace](0)->FastGetSolutionStepValue(VELOCITY) * Coeff[0];
                     other_to_me_vel += ineighbour->GetGeometry().Edges()[RightFace](1)->FastGetSolutionStepValue(VELOCITY) * Coeff[1];
               }
               else if( dim == 3)
               {
                   for(std::size_t ifnode = 0; ifnode < ineighbour->GetGeometry().Faces()[RightFace].size(); ifnode++)
                   {
                       other_to_me_vel += ineighbour->GetGeometry().Faces()[RightFace](ifnode)->FastGetSolutionStepValue(VELOCITY) * Coeff[ifnode];
                   }
               }
            }


            double DeltDisp[3] = {0.0};
            double DeltVel [3] = {0.0};

            DeltVel[0] = (vel[0] - other_to_me_vel[0]);
            DeltVel[1] = (vel[1] - other_to_me_vel[1]);
            DeltVel[2] = (vel[2] - other_to_me_vel[2]);

            // For translation movement delt displacement
            DeltDisp[0] = DeltVel[0] * mTimeStep;
            DeltDisp[1] = DeltVel[1] * mTimeStep;
            DeltDisp[2] = DeltVel[2] * mTimeStep;


            if ( rotation_OPTION == 1 )
            {
                double velA[3]   = {0.0};
                double dRotaDisp[3] = {0.0};

                array_1d<double, 3 > AngularVel= GetGeometry()(0)->FastGetSolutionStepValue(ANGULAR_VELOCITY);
                double Vel_Temp[3] = { AngularVel[0], AngularVel[1], AngularVel[2]};
                GeometryFunctions::CrossProduct(Vel_Temp, LocalCoordSystem[2], velA);

                dRotaDisp[0] = -velA[0] * radius;
                dRotaDisp[1] = -velA[1] * radius;
                dRotaDisp[2] = -velA[2] * radius;

                //////contribution of the rotation vel
                DeltDisp[0] += dRotaDisp[0] * mTimeStep;
                DeltDisp[1] += dRotaDisp[1] * mTimeStep;
                DeltDisp[2] += dRotaDisp[2] * mTimeStep;
            }


            double LocalDeltDisp[3] = {0.0};
            GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, DeltDisp, LocalDeltDisp);

            //////120323,for global storage

            GlobalContactForceOld[0] = mContactForces[3 * iContactForce  + 0 ];
            GlobalContactForceOld[1] = mContactForces[3 * iContactForce  + 1 ];
            GlobalContactForceOld[2] = mContactForces[3 * iContactForce  + 2 ];

            GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, GlobalContactForceOld, LocalContactForce);
            LocalContactForce[0] +=  - ks * LocalDeltDisp[0];
            LocalContactForce[1] +=  - ks * LocalDeltDisp[1];
            LocalContactForce[2] +=  - kn * LocalDeltDisp[2];


            /////Cfeng:120514 Absolutal method for normal force of the block and particle contact is adopted.
            /////But this method will lead to enery out of convergence, so command out now,120514, afternoon
            ////////LocalContactForce[2] = (radius - DistPToB) * kn;

            /// Cfeng: For Tension Failure
            if (-LocalContactForce[2] > (Tension * area))
            {
                LocalContactForce[0] = 0.0;
                LocalContactForce[1]  = 0.0;
                LocalContactForce[2]  = 0.0;

                mContactFailureType[iContactForce] = 1;
            }
            else
            {

                double ShearForceMax = LocalContactForce[2] * Friction + Cohesion * area;
                double ShearForceNow = sqrt(LocalContactForce[0] * LocalContactForce[0]
                                     +      LocalContactForce[1] * LocalContactForce[1]);


                //Cfeng: if tensile, only cohesion could be introduced
                if(LocalContactForce[2] < 0.0)
                {
                    ShearForceMax = Cohesion * area;
                }


                //Cfeng: for shear failure
                if(ShearForceMax == 0.0)
                {
                    LocalContactForce[0] = 0.0;
                    LocalContactForce[1] = 0.0;
                }
                else if(ShearForceNow > ShearForceMax)
                {
                    LocalContactForce[0] = ShearForceMax / ShearForceNow * LocalContactForce[0];
                    LocalContactForce[1] = ShearForceMax / ShearForceNow * LocalContactForce[1];

                    mContactFailureType[iContactForce] = 2;
                }
            }


            GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalContactForce, GlobalContactForce);

            mContactForces[3 * iContactForce + 0 ] = GlobalContactForce[0];
            mContactForces[3 * iContactForce + 1 ] = GlobalContactForce[1];
            mContactForces[3 * iContactForce + 2 ] = GlobalContactForce[2];

            force[0] += GlobalContactForce[0];
            force[1] += GlobalContactForce[1];
            force[2] += GlobalContactForce[2];

            if ( rotation_OPTION == 1 )
            {
                double MA[3] = {0.0};
                GeometryFunctions::CrossProduct(LocalCoordSystem[2], GlobalContactForce, MA);
                mRota_Moment[0] -= MA[0] * radius;
                mRota_Moment[1] -= MA[1] * radius;
                mRota_Moment[2] -= MA[2] * radius;
            }


            if(RightFace != -1)
            {
               if(dim == 2)
               {
                    ineighbour->GetGeometry().Edges()[RightFace](0)->FastGetSolutionStepValue(RHS)[0] -= GlobalContactForce[0] * Coeff[0];
                    ineighbour->GetGeometry().Edges()[RightFace](0)->FastGetSolutionStepValue(RHS)[1] -= GlobalContactForce[1] * Coeff[0];
                    ineighbour->GetGeometry().Edges()[RightFace](0)->FastGetSolutionStepValue(RHS)[2] -= GlobalContactForce[2] * Coeff[0];

                    ineighbour->GetGeometry().Edges()[RightFace](1)->FastGetSolutionStepValue(RHS)[0] -= GlobalContactForce[0] * Coeff[1];
                    ineighbour->GetGeometry().Edges()[RightFace](1)->FastGetSolutionStepValue(RHS)[1] -= GlobalContactForce[1] * Coeff[1];
                    ineighbour->GetGeometry().Edges()[RightFace](1)->FastGetSolutionStepValue(RHS)[2] -= GlobalContactForce[2] * Coeff[1];
               }
               else if( dim == 3)
               {
                   for(std::size_t ifnode = 0; ifnode < ineighbour->GetGeometry().Faces()[RightFace].size(); ifnode++)
                   {
                       ineighbour->GetGeometry().Faces()[RightFace](ifnode)->FastGetSolutionStepValue(RHS)[0] -= GlobalContactForce[0] * Coeff[ifnode];
                       ineighbour->GetGeometry().Faces()[RightFace](ifnode)->FastGetSolutionStepValue(RHS)[1] -= GlobalContactForce[1] * Coeff[ifnode];
                       ineighbour->GetGeometry().Faces()[RightFace](ifnode)->FastGetSolutionStepValue(RHS)[2] -= GlobalContactForce[2] * Coeff[ifnode];
                   }
               }
            }

            iContactForce++;

        }
    }

    void DEM_FEM_Particle::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo)
    {

        if (rVariable == DELTA_TIME)
        {
            double E = GetProperties()[YOUNG_MODULUS];
            double K = E * M_PI * this->GetGeometry()(0)->FastGetSolutionStepValue(RADIUS);
            Output = sqrt( mRealMass / K);

            if(rCurrentProcessInfo[ROTATION_OPTION] == 1)
            {
                Output = Output * 0.5;
            }
        }
    }

    //************************************************************************************
    //************************************************************************************
  
    void DEM_FEM_Particle::save( Serializer& rSerializer ) const
    {
//  std::cout << "Saving the Particle #" << Id() << std::endl;
        rSerializer.save( "Name", "Particle" );
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element );
    }

    void DEM_FEM_Particle::load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element );
    }



} // Namespace Kratos


