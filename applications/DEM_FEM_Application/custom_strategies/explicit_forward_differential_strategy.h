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
/* *********************************************************   
*          
*   Last Modified by:    $Author: Nelson $
*   Date:                $Date: 2009-09-18 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/


/// WARNING = Los desplazameitos son obtenidos en le paso
/// n+1; las velocidades  y aceleraciones son obtenidas para el paso n.

#if !defined(KRATOS_DEM_FEM_COUPLED_EXPLICIT_FORWARD_STRATEGY)
#define  KRATOS_DEM_FEM_COUPLED_EXPLICIT_FORWARD_STRATEGY


/* System includes */
#include <limits>
#include<iostream>
#include<iomanip>

/////////#define _OPENMP

/* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"
#include "includes/variables.h"
#include "containers/array_1d.h"
#include "particle_block_configure.h"
#include "spatial_containers/spatial_containers.h"
#include "dem_fem__application.h"
#include "custom_utilities/neighbours_calculator.h"
#include "custom_elements/discrete_element.h"
namespace Kratos
{
  
 template<
 class TSparseSpace,
 class TDenseSpace, 
 class TLinearSolver> 
 class DEM_FEM_Explicit_Forward_Differential_Strategy : public SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>
     {

	  public:

	  KRATOS_CLASS_POINTER_DEFINITION(DEM_FEM_Explicit_Forward_Differential_Strategy);

	  typedef SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver> BaseType;

	  typedef typename BaseType::TDataType TDataType;
	  
	  typedef TSparseSpace SparseSpaceType;

	  typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;
	  
	  typedef typename BaseType::TSchemeType TSchemeType;

	  typedef typename BaseType::DofsArrayType DofsArrayType;

	  typedef typename Element::DofsVectorType DofsVectorType;

	  typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

	  typedef typename BaseType::TSystemVectorType TSystemVectorType;

	  typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

	  typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

	  typedef ModelPart::NodesContainerType NodesArrayType;

	  typedef ModelPart::ElementsContainerType ElementsArrayType;

	  typedef ModelPart::ConditionsContainerType ConditionsArrayType;
	  
	  typedef ModelPart::ConditionsContainerType::ContainerType ConditionsContainerType;
      
	  typedef ConditionsContainerType::iterator                 ConditionsContainerIterator;

	  typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;

	  typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

          typedef ModelPart::PropertiesType PropertiesType;


          typedef Element::Pointer ParticlePointer;
          typedef typename std::vector<ParticlePointer> ParticlePointerVector;
          typedef typename std::vector<ParticlePointer>::iterator ParticlePointerIterator;

          typedef WeakPointerVector<Element > ParticleWeakVector;
          typedef WeakPointerVector<Element >::iterator ParticleWeakIterator;
	  
	  struct NFACE
          {
              ParticlePointerVector pBlock;
              std::vector< std::size_t > iface;
              std::vector<std::size_t>globnode;
          };


	  DEM_FEM_Explicit_Forward_Differential_Strategy(
	                ModelPart& model_part, 
			const int        dimension,
                        const int        damp_type,
			const double     damping_ratio,
                        const bool       virtual_mass,
                        const double     contact_stiffness_ratio,
                        const double     max_delta_time,
			const bool       CalculateReactions,
			const bool       ComputeFemFemContact,
			const bool       MoveMeshFlag,
			typename         TLinearSolver::Pointer pNewLinearSolver,
			typename         TSchemeType::Pointer pScheme,
			typename         TBuilderAndSolverType::Pointer pNewBuilderAndSolver
			)
			
	  : SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(model_part, MoveMeshFlag)
	      {
		std::cout<< "******************************************************"<< std::endl;
	        std::cout <<"   DYNAMIC SOLVER ANALYSIS FOR DEM_FEM COUPLED METHOD  "<< std::endl;
                std::cout <<"     TIME INTEGRATION METHOD  =  Forward Differences    "<< std::endl;
                std::cout <<"IMPLEMENTED BY = Chun feng Based on the Version of Nelson"<< std::endl;
                std::cout<< "******************************************************"<< std::endl;
               
		
	        mdimension                 = dimension; 
                mCalculateReactionsFlag    = CalculateReactions;                
                mElementsAreInitialized    = false;
		mConditionsAreInitialized  = false; 
		mInitializeWasPerformed    = false;
		mComputeTime               = false;

		mtimestep                  = max_delta_time;             

                mpLinearSolver       = pNewLinearSolver;
		mpBuilderAndSolver   = pNewBuilderAndSolver;
		mpScheme             = pScheme;		
		mpBuilderAndSolver->SetReshapeMatrixFlag(false);
                
  	
		mdamping_ratio       = damping_ratio;
		malpha_damp          = 0.00;
		mbeta_damp           = 0.00; 
		mpenalty_factor      = 1.0;


                mdamp_type               = damp_type;
                mvirtual_mass            = virtual_mass;
                mcontact_stiffness_ratio = contact_stiffness_ratio;
                mComputeFemFemContact    = ComputeFemFemContact;
                mnumber_step             = 0;

                mUnbalRatioCal           = 1;
                mIfHaveCalVirtualMass    = false;
                mIfParticleInitialSearch = true;
                mIfParticleBlockInitialSearch = true;
                mIfHaveInitialzeSolutionStep = false;
            
	      }

	  virtual ~DEM_FEM_Explicit_Forward_Differential_Strategy () {}
	           
double Solve()
      {

	KRATOS_TRY

        ModelPart& r_model_part              = BaseType::GetModelPart();

	ProcessInfo& CurrentProcessInfo      = r_model_part.GetProcessInfo();

        if(mInitializeWasPerformed == false)
        {
	  Initialize();
	}


	if(mComputeTime==false)
        {
	    ComputeCriticalTime();
	    mComputeTime = true;
	}


        ///Initialize solution step
        ///call this function once or in every time step should be ensured,this version is call this function each time step.
        if(mIfHaveInitialzeSolutionStep == false)
        {
            InitializeSolutionStep();
            mIfHaveInitialzeSolutionStep = true;
        }


        ////Cfeng: find the paticle-particle neighbours
        FindParticleNeighbours();



        ////Cfeng: find the block neighbours for particles
        FindParticleBlockNeighbours();



        ////////Cfeng: Calcualte the deformation force, and add them to the node
        Calculate_Element_Deformation_Force_Add_To_Node();


        ///Cfeng: cal damping forces
	ComputeDampingForces();


        ////////Cfeng:Calculate the movementment of the node;
	Calculate_Node_Movement_By_Newton_Law();


        ////Cfeng: Calculate Particle Rotation Evolvement
        if(CurrentProcessInfo[ROTATION_OPTION] == 1)
        {
           Calculate_Particle_Rotate_Evolvement();
        }

        ///Cfeng:Cal reaction
	if(mCalculateReactionsFlag)
        {
            CalculateReaction();
        }


	/// Actualizacion de los desplazamientos
 	if(BaseType::MoveMeshFlag() == true)
        {
            BaseType::MoveMesh();
        }



        ///Finalize solution step
	FinalizeSolutionStep();

        ///For Visualization, 120530, the dirived version, the number of neighbour is located on elements.
        TransferElementNeighboursToNodes();


        mnumber_step++;
        ////////////////////Cfeng: set the model time now
        r_model_part.CloneTimeStep(mtimestep * mnumber_step);


        if(mnumber_step % 100 == 0)
        {
            std::cout<< "Time_Step = "<<mnumber_step<<" Time_Now = " << mtimestep * mnumber_step << " The Unbalance Ratio = "<< mUnbalRatioCal  << std::endl<< std::endl;
        }

	return 0.00;
	KRATOS_CATCH("")
      }


void Initialize()
{

    KRATOS_TRY


   

    /// Inicializando los elemtos
    if(mElementsAreInitialized == false)
    {
        InitializeElements();

        ////Cfeng: Identify is the element is boundary element and if the face is boundary face.
        CreatAuxParticleBlockList();

        ////Cfeng: after connect to Miquel, Initialize neighbour search for particle120530
        FindParticleNeighbours();

        mElementsAreInitialized   = true;
    }

    /// Inicializando las condiciones
    if(mConditionsAreInitialized == false)
	InitializeConditions();

    if(mIfHaveCalVirtualMass == false)
        CalculateVirtualMass();


    mInitializeWasPerformed   = true;


    KRATOS_CATCH("")

}

void FindBoundaryFaceForFemElements()
{
    KRATOS_TRY

    ModelPart& r_model_part = BaseType::GetModelPart(); 

    ElementsArrayType& pElements     = r_model_part.Elements();
    typename ElementsArrayType::iterator it_begin=pElements.begin();
    typename ElementsArrayType::iterator it_end=pElements.end();

    /////Cfeng: Initialize the IF_BOUNDARY FLAG for element and edges (2d) or faces(3D)
    /////Cfeng: The initial value is the boundary =1 means boundary
    
    for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
    {     
       std::size_t dim = it->GetGeometry().WorkingSpaceDimension();
       if(dim == 2 && it->GetGeometry().size() > 1)   /////Cfeng: it->GetGeometry().size() > 1 FEM element == 1 particle
       {
           it->GetValue(IF_BOUNDARY_ELEMENT) = 1;

           std::size_t Nedge = it->GetGeometry().Edges().size();
           it->GetValue(IF_BOUNDARY_FACE).resize(Nedge);

           for (std::size_t iedge = 0; iedge < Nedge; iedge++)
           {
              it->GetValue(IF_BOUNDARY_FACE)[iedge] = 1;
           }

       }
       else if(dim == 3 && it->GetGeometry().size() > 1)
       {
           it->GetValue(IF_BOUNDARY_ELEMENT) = 1;

           std::size_t Nface = it->GetGeometry().Faces().size();
           it->GetValue(IF_BOUNDARY_FACE).resize(Nface);

           for (std::size_t iface = 0; iface < Nface; iface++)
           {
              it->GetValue(IF_BOUNDARY_FACE)[iface] = 1;
           }
       }
    }


    std::size_t TotalNodeNum = r_model_part.NumberOfNodes();

    std::vector< std::vector<std::size_t> > aiNodeFaces;
    aiNodeFaces.resize(TotalNodeNum);


    std::vector<NFACE> oface;

    std::size_t iface,inode,No;
    std::size_t i,j;
    bool if_insert = false;

    for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
    {
        std::size_t dim = it->GetGeometry().WorkingSpaceDimension();
        if(dim == 2 && it->GetGeometry().size() > 1)
        {
            for(iface = 0; iface < it->GetGeometry().Edges().size(); iface++)
            {
    
                std::vector<std::size_t> Gno;
                for(inode = 0; inode < it->GetGeometry().Edges()[iface].size(); inode++)
                {
                    No = it->GetGeometry().Edges()[iface](inode)->Id();

                    if_insert = false;

                    for(i = 0; i < Gno.size(); i++)
                    {
                        if(No < Gno[i])
                        {
                            Gno.insert(Gno.begin() + i, No);
                            if_insert = true;
                            break;
                        }
                    }

                    if(if_insert == false)
                    {
                        Gno.push_back(No);
                    }
                }

                std::vector<std::size_t> &rFaces = aiNodeFaces[Gno[0]];


                for(i = 0; i < rFaces.size(); i++)
                {
                    if(oface[rFaces[i]].globnode.size() == Gno.size())
                    {
                        std::size_t flag = 0;
                        for(j = 0; j < oface[rFaces[i]].globnode.size(); j++)
                        {
                            if(oface[rFaces[i]].globnode[j] == Gno[j])
                            {
                                flag++;
                            }
                        }
                        if(flag == Gno.size())
                        {
                            oface[rFaces[i]].pBlock.push_back( *(it.base()) );
                            oface[rFaces[i]].iface.push_back(iface);
                            break;
                        }

                    }
                }

                if(i ==  rFaces.size())
                {
                    rFaces.push_back(oface.size());

                    oface.resize(oface.size()+1);
                    for(j = 0 ; j < Gno.size(); j++)
                    {
                        oface[oface.size() - 1].globnode.push_back(Gno[j]);
                    }
                    oface[oface.size() - 1].pBlock.push_back(*(it.base()));
                    oface[oface.size() - 1].iface.push_back(iface);
                }
            }
        }
        else if(dim == 3 && it->GetGeometry().size() > 1)
        {
            for(iface = 0; iface < it->GetGeometry().Faces().size(); iface++)
            {

                std::vector<std::size_t> Gno;
                for(inode = 0; inode < it->GetGeometry().Faces()[iface].size(); inode++)
                {
                    No = it->GetGeometry().Faces()[iface](inode)->Id();

                    if_insert = false;

                    for(i = 0; i < Gno.size(); i++)
                    {
                        if(No < Gno[i])
                        {
                            Gno.insert(Gno.begin() + i, No);
                            if_insert = true;
                            break;
                        }
                    }

                    if(if_insert == false)
                    {
                        Gno.push_back(No);
                    }
                }

                std::vector<std::size_t> &rFaces = aiNodeFaces[Gno[0]];


                for(i = 0; i < rFaces.size(); i++)
                {
                    if(oface[rFaces[i]].globnode.size() == Gno.size())
                    {
                        std::size_t flag = 0;
                        for(j = 0; j < oface[rFaces[i]].globnode.size(); j++)
                        {
                            if(oface[rFaces[i]].globnode[j] == Gno[j])
                            {
                                flag++;
                            }
                        }
                        if(flag == Gno.size())
                        {
                            oface[rFaces[i]].pBlock.push_back(*(it.base()));
                            oface[rFaces[i]].iface.push_back(iface);
                            break;
                        }

                    }
                }

                if(i ==  rFaces.size())
                {
                    rFaces.push_back(oface.size());

                    oface.resize(oface.size()+1);
                    for(j = 0 ; j < Gno.size(); j++)
                    {
                        oface[oface.size() - 1].globnode.push_back(Gno[j]);
                    }
                    oface[oface.size() - 1].pBlock.push_back(*(it.base()));
                    oface[oface.size() - 1].iface.push_back(iface);
                }
            }
        }
    }
       
  
    for(i = 0 ; i< oface.size(); i++)
    {
        if( oface[i].pBlock.size() > 1)
        {
            ParticlePointer  b0 = oface[i].pBlock[0];
            ParticlePointer  b1 = oface[i].pBlock[1];

            int i0 = oface[i].iface[0];
            int i1 = oface[i].iface[1];
  
            b0->GetValue(IF_BOUNDARY_ELEMENT) = 0;
            b1->GetValue(IF_BOUNDARY_ELEMENT) = 0;
            
            b0->GetValue(IF_BOUNDARY_FACE)[i0] = 0;
            b1->GetValue(IF_BOUNDARY_FACE)[i1] = 0;                           
        }
    }


    for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
    {
       std::size_t dim = it->GetGeometry().WorkingSpaceDimension();
       if(dim == 2 && it->GetGeometry().size() > 1)   /////Cfeng: it->GetGeometry().size() > 1 FEM element == 1 particle
       {
           for (std::size_t iedge = 0; iedge < it->GetGeometry().Edges().size(); iedge++)
           {
               if(it->GetValue(IF_BOUNDARY_FACE)[iedge] == 1)
               {
                   it->GetValue(IF_BOUNDARY_ELEMENT) = 1;
                   break;
               }
           }

       }
       else if(dim == 3 && it->GetGeometry().size() > 1)
       {
           for (std::size_t iface = 0; iface < it->GetGeometry().Faces().size(); iface++)
           {
               if(it->GetValue(IF_BOUNDARY_FACE)[iface] == 1)
               {
                   it->GetValue(IF_BOUNDARY_ELEMENT) = 1;
                   break;
               }
           }
       }
    }


    KRATOS_CATCH("")
  
}


void CreatAuxParticleBlockList()
{
   KRATOS_TRY

   /////// Important, find the boundary face for particle block contact
   FindBoundaryFaceForFemElements();
  
   ModelPart& r_model_part = BaseType::GetModelPart();
   ElementsArrayType& pElements     = r_model_part.Elements();

   typename ElementsArrayType::iterator it_begin=pElements.begin();
   typename ElementsArrayType::iterator it_end=pElements.end();

   for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
   {
       if(it->GetGeometry().size() == 1)
       {
           aux_list_of_particles.push_back(   *(it.base())   );
           tempParticle.push_back         (   *(it.base())   );
       }
       else
       {         
           if(it->GetValue(IF_BOUNDARY_ELEMENT) == 1)
           {
               aux_list_of_blocks.push_back(   *(it.base())   );
           }
       }
   }

     KRATOS_CATCH("")
 
}



void ComputeDampingForces()
{
    if(mdamp_type == 1)
    {
        ComputeNonViscousDampingForces();
    }
    else if(mdamp_type == 2)
    {
        ComputeViscousDampingForces();
    }
}

void ComputeViscousDampingForces()
{
}


void ComputeNonViscousDampingForces()
{   
      KRATOS_TRY
      ModelPart& r_model_part          = BaseType::GetModelPart();
      NodesArrayType& pNodes           = r_model_part.Nodes(); 

	typename NodesArrayType::iterator it_begin=pNodes.ptr_begin();
	typename NodesArrayType::iterator it_end=pNodes.ptr_end();
	for (NodesArrayType::iterator it= it_begin; it!=it_end; ++it)
	 {
          
              ///Cfeng:120416,use &
            array_1d<double,3> &rhs        = it->FastGetSolutionStepValue(RHS);


	    const array_1d<double,3> Vel = it->FastGetSolutionStepValue(VELOCITY);


            for (int iDof = 0; iDof < 3; iDof++)
            {
                if (Vel[iDof] > 0.0)
                {
                    rhs[iDof] = rhs[iDof] - mdamping_ratio * fabs(rhs[iDof]);
                }
                else
                {
                    rhs[iDof] = rhs[iDof] + mdamping_ratio * fabs(rhs[iDof]);
                }
            }


	  }  
	
	KRATOS_CATCH("")
    }
 



void Calculate_Node_Movement_By_Newton_Law()
{

    KRATOS_TRY

    ModelPart& r_model_part          = BaseType::GetModelPart();
    ProcessInfo& CurrentProcessInfo  = r_model_part.GetProcessInfo();
    NodesArrayType& pNodes           = r_model_part.Nodes();


    #ifdef _OPENMP
    int number_of_threads = omp_get_max_threads();
    #else
    int number_of_threads = 1;
    #endif

    vector<unsigned int> node_partition;
    CreatePartition(number_of_threads, pNodes.size(), node_partition);

   // #pragma omp parallel for private(mid_pos_velocity, mid_neg_velocity)

    array_1d<double,3> UnbalForce;
    UnbalForce[0] = 0.0;
    UnbalForce[1] = 0.0;
    UnbalForce[2] = 0.0;

    array_1d<double,3> ReactionForce;
    ReactionForce[0] = 0.0;
    ReactionForce[1] = 0.0;
    ReactionForce[2] = 0.0;

    for(int k=0; k<number_of_threads; k++)
      {
	typename NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
	typename NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];

	for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)
	{
            array_1d<double,3>& acceleration = (i->FastGetSolutionStepValue(ACCELERATION));
            array_1d<double,3>& velocity     = (i->FastGetSolutionStepValue(VELOCITY));
	    array_1d<double,3>& displacement = i->FastGetSolutionStepValue(DISPLACEMENT);
	    array_1d<double,3>& nodeForce    = i->FastGetSolutionStepValue(RHS);

	    const double& nodal_mass         =  i->FastGetSolutionStepValue(NODAL_MASS);

            acceleration = nodeForce / nodal_mass;

            array_1d<double,3> velocity_old = velocity;
            array_1d<double,3> velocity_new = velocity_old + acceleration * mtimestep;




            ////Cfeng : fixed node should set acc-vel-dis zero
            if( (i->pGetDof(VELOCITY_X))->IsFixed() == true)
	    {
                acceleration[0] = 0.0;
                displacement[0] += velocity[0] * mtimestep;

                ReactionForce[0] += fabs (nodeForce[0]);
            }
            else
            {
               velocity    [0]  = (velocity_old[0] + velocity_new[0]) * 0.5;
               displacement[0] += velocity[0] * mtimestep;

               UnbalForce[0] += fabs (nodeForce[0]);
            }


            if( (i->pGetDof(VELOCITY_Y))->IsFixed() == true)
	    {
                acceleration[1] = 0.0;
                displacement[1] += velocity[1] * mtimestep;

                ReactionForce[1] += fabs (nodeForce[1]);
            }
            else
            {
               velocity    [1]  = (velocity_old[1] + velocity_new[1]) * 0.5;
               displacement[1] += velocity[1] * mtimestep;

               UnbalForce[1] += fabs (nodeForce[1]);
            }


            if( (i->pGetDof(VELOCITY_Z))->IsFixed() == true)
	    {
                acceleration[2] = 0.0;
                displacement[2] += velocity[2] * mtimestep;

                ReactionForce[2] += fabs (nodeForce[2]);
            }
            else
            {
                velocity    [2]  = (velocity_old[2] + velocity_new[2]) * 0.5;
                displacement[2] += velocity[2] * mtimestep;

                UnbalForce[2] += fabs (nodeForce[2]);
            }

        }
    }

    double UnBalForceRatio = sqrt(UnbalForce[0] * UnbalForce[0] + UnbalForce[1] * UnbalForce[1] + UnbalForce[2] * UnbalForce[2]);

    double Reaction = sqrt(ReactionForce[0] * ReactionForce[0] + ReactionForce[1] * ReactionForce[1] + ReactionForce[2] * ReactionForce[2]);


    mUnbalRatioCal = UnBalForceRatio / Reaction;

    if(Reaction < 1.0e-6)
    {
        mUnbalRatioCal = UnBalForceRatio;
    }


    CurrentProcessInfo[DEM_FEM_CONVERGENCE_RATIO] = mUnbalRatioCal;

    KRATOS_CATCH("")
}



void TransferElementNeighboursToNodes()
{
    KRATOS_TRY


      ModelPart& r_model_part          = BaseType::GetModelPart();
      ElementsArrayType& pElements     = r_model_part.Elements();


    typename ElementsArrayType::iterator it_begin=pElements.ptr_begin();
    typename ElementsArrayType::iterator it_end=pElements.ptr_end();
    for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
    {
        Element::GeometryType& geom = it->GetGeometry();
        for (unsigned int i = 0; i <geom.size(); i++)
        {
            geom(i)->FastGetSolutionStepValue(PARTICLE_NUMBER_OF_NEIGHBOURS) = it->GetValue(NUMBER_OF_NEIGHBOURS);
        }
    }


    KRATOS_CATCH("")
}




void Calculate_Particle_Rotate_Evolvement()
{
    KRATOS_TRY

    ModelPart& r_model_part          = BaseType::GetModelPart();
    ElementsArrayType& pElements     = r_model_part.Elements();

    typename ElementsArrayType::iterator it_begin=pElements.ptr_begin();
    typename ElementsArrayType::iterator it_end=pElements.ptr_end();
    for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
    {     
       Element::GeometryType& geom = it->GetGeometry();

       if (geom.size() == 1)
       {
            double PMass            = geom(0)->FastGetSolutionStepValue(NODAL_MASS);
            double PMomentOfInertia = geom(0)->FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA);

            array_1d<double, 3 > & AngularVel = geom(0)->FastGetSolutionStepValue(ANGULAR_VELOCITY);
            array_1d<double, 3 > & RotaMoment = geom(0)->FastGetSolutionStepValue(PARTICLE_MOMENT);
            array_1d<double, 3 > & Rota_Displace = geom(0)->FastGetSolutionStepValue(PARTICLE_ROTATION_ANGLE);
       
            bool If_Fix_Rotation[3] = {false, false, false};
            If_Fix_Rotation[0] = geom(0)->pGetDof(VELOCITY_X)->IsFixed();
            If_Fix_Rotation[1] = geom(0)->pGetDof(VELOCITY_Y)->IsFixed();
            If_Fix_Rotation[2] = geom(0)->pGetDof(VELOCITY_Z)->IsFixed();

            for(std::size_t i = 0 ; i < 3; i++)
            {
                if(If_Fix_Rotation[i] == false)
                {
                     double RotaAcc = 0.0;
                     if(AngularVel[i] > 0.0)
                     {
                         RotaAcc = (RotaMoment[i] - mdamping_ratio * fabs(RotaMoment[i])) / PMass / PMomentOfInertia;
                     }
                     else
                     {
                         RotaAcc = (RotaMoment[i] + mdamping_ratio * fabs(RotaMoment[i])) / PMass / PMomentOfInertia;
                     }

                     double RotaVelOld = AngularVel[i];
                     double RotaVelNew = RotaVelOld + RotaAcc * mtimestep;

                     AngularVel[i]  = 0.5 * (RotaVelOld + RotaVelNew);
                     Rota_Displace[i] += AngularVel[i] * mtimestep / M_PI * 180.0;
                }
                RotaMoment[i] = 0.0;
            }
       }     
    }

    KRATOS_CATCH("")
}

void CalculateVirtualMass()
{
    KRATOS_TRY

    if(mvirtual_mass == true)
    {   
      ModelPart& r_model_part          = BaseType::GetModelPart();
      ElementsArrayType& pElements     = r_model_part.Elements();

      ProcessInfo& CurrentProcessInfo  = r_model_part.GetProcessInfo();
      
      //NodesArrayType &pNode            = r_model_part.Nodes();
      typename NodesArrayType::iterator inode;
      for(inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); inode++)
      {
          inode->FastGetSolutionStepValue(NODAL_MASS) = 0.0;
      }

	typename ElementsArrayType::iterator it_begin=pElements.ptr_begin();
	typename ElementsArrayType::iterator it_end=pElements.ptr_end();
	for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
	{

           double Young  = it->GetProperties()[YOUNG_MODULUS];        
           double Length = it->GetGeometry().Length();
           double Volume = 0.0;
           double VirtualMass = 0.0;

           Element::GeometryType& geom = it->GetGeometry();

           if (geom.size() == 1)
           {
              VirtualMass = Young * M_PI * it->GetGeometry()(0)->FastGetSolutionStepValue(RADIUS);
              if(CurrentProcessInfo[ROTATION_OPTION] == 1)
              {
                  VirtualMass = VirtualMass * 2.5;
              }
           }
           else if (it->GetGeometry().Dimension() == 2 && geom.size() > 2)
           {
               Volume = it->GetGeometry().Area();

               VirtualMass = Young / (Length * Length) * Volume;
           } 
           else if (it->GetGeometry().Dimension() == 3 && geom.size() > 3 )
           {
               Volume = it->GetGeometry().Volume();

               VirtualMass = Young / (Length * Length) * Volume;
           }

	    
	    for (unsigned int i = 0; i <geom.size(); i++)
	     {
	        double& mass = geom(i)->FastGetSolutionStepValue(NODAL_MASS);
		mass  = mass + VirtualMass / (double)geom.size();
	     }
	}
    }

    mIfHaveCalVirtualMass = true;

    KRATOS_CATCH("")
}

//***************************************************************************
//***************************************************************************

void InitializeElements()
{
      KRATOS_TRY
      ModelPart& r_model_part          = BaseType::GetModelPart();
      ProcessInfo& CurrentProcessInfo  = r_model_part.GetProcessInfo();  
      ElementsArrayType& pElements     = r_model_part.Elements(); 

      Matrix MassMatrix;
      
      unsigned int index = 0;
      
      
	typename ElementsArrayType::iterator it_begin=pElements.ptr_begin();
	typename ElementsArrayType::iterator it_end=pElements.ptr_end();
	for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
	  {
	    Element::GeometryType& geom = it->GetGeometry(); // Nodos del elemento
	    (it)->Initialize(); 
	    (it)->GetValue(IS_INACTIVE) = false;
	    (it)->MassMatrix(MassMatrix, CurrentProcessInfo);
	    const unsigned int& dim   = geom.WorkingSpaceDimension();
         
            if(geom.size() > 1)
            {
                index = 0;
	        for (unsigned int i = 0; i <geom.size(); i++)
	        {
	            double& mass = geom(i)->FastGetSolutionStepValue(NODAL_MASS);
		    geom(i)->SetLock();
		    index = i*dim;
		    mass  = mass + MassMatrix(index,index);
		    geom(i)->UnSetLock();
	        }
            }
	 }
     
     KRATOS_CATCH("")
}

//***************************************************************************
//***************************************************************************

///WARNING = Falta colocar el contacto
void InitializeConditions()
{
  
       KRATOS_TRY
       
      ModelPart& r_model_part          = BaseType::GetModelPart();  
      //ProcessInfo& CurrentProcessInfo  = r_model_part.GetProcessInfo();   
      ConditionsArrayType& pConditions = r_model_part.Conditions();
      
      #ifdef _OPENMP
      int number_of_threads = omp_get_max_threads();
      #else
      int number_of_threads = 1;
      #endif

      vector<unsigned int> condition_partition;
      CreatePartition(number_of_threads, pConditions.size(), condition_partition);
      
      #pragma omp parallel for
      for(int k=0; k<number_of_threads; k++)
      {
	typename ConditionsArrayType::iterator it_begin=pConditions.ptr_begin()+condition_partition[k];
	typename ConditionsArrayType::iterator it_end=pConditions.ptr_begin()+condition_partition[k+1];
  
         for (ConditionsArrayType::iterator it= it_begin; it!=it_end; ++it)
          {
	    (it) -> Initialize();
	  }

      }
    
      
      mConditionsAreInitialized = true;
      
      KRATOS_CATCH("")
      
}
  
//***************************************************************************
//***************************************************************************

void InitializeSolutionStep()
{
    KRATOS_TRY
    ModelPart& r_model_part          = BaseType::GetModelPart();  
    ElementsArrayType& pElements     = r_model_part.Elements();
    ProcessInfo& CurrentProcessInfo  = r_model_part.GetProcessInfo();
    ConditionsArrayType& pConditions = r_model_part.Conditions();

    #ifdef _OPENMP
    int number_of_threads = omp_get_max_threads();
    #else
    int number_of_threads = 1;
    #endif

    vector<unsigned int> element_partition;
    vector<unsigned int> condition_partition;
    CreatePartition(number_of_threads, pElements.size(), element_partition);
    CreatePartition(number_of_threads, pConditions.size(), condition_partition);

    #pragma omp parallel for
    for(int k=0; k<number_of_threads; k++)
    {
      typename ElementsArrayType::iterator it_begin=pElements.ptr_begin()+element_partition[k];
      typename ElementsArrayType::iterator it_end=pElements.ptr_begin()+element_partition[k+1];

      for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
	  {
	    (it) -> InitializeSolutionStep(CurrentProcessInfo);
	  }
    }

    #pragma omp parallel for
    for(int k=0; k<number_of_threads; k++)
      {
	typename ConditionsArrayType::iterator it_begin=pConditions.ptr_begin()+condition_partition[k];
	typename ConditionsArrayType::iterator it_end=pConditions.ptr_begin()+condition_partition[k+1];
  
         for (ConditionsArrayType::iterator it= it_begin; it!=it_end; ++it)
          {
	    (it) -> InitializeSolutionStep(CurrentProcessInfo);
	  }

    }


    
    KRATOS_CATCH("")
}

//***************************************************************************
//***************************************************************************

void ComputeCriticalTime()
{
    KRATOS_TRY

    ModelPart& r_model_part          = BaseType::GetModelPart();
    ProcessInfo& CurrentProcessInfo  = r_model_part.GetProcessInfo();

    if(mvirtual_mass == true)
    {
        if(mtimestep > 0.9)
        {
            mtimestep = 0.9;
        }

        std::cout<<"******************Virtual Mass TimeStep is Used******************" <<std::endl;
    }
    else
    {          
          double TimeStepTemp = 0.0;
   
          ElementsArrayType& pElements     = r_model_part.Elements();

          typename ElementsArrayType::iterator it_begin = pElements.ptr_begin();
          typename ElementsArrayType::iterator it_end   = pElements.ptr_end();

          for(ElementsArrayType::iterator it = it_begin; it!= it_end; it++)
          {
              it->Calculate(DELTA_TIME, TimeStepTemp, CurrentProcessInfo);
	      
              if(mtimestep > TimeStepTemp)
              {
                  mtimestep = TimeStepTemp;
              }
          }

          std::cout<<"******************Real Mass TimeStep is Used******************" <<std::endl;
    }


    CurrentProcessInfo[DEM_DELTA_TIME]  = mtimestep;

    std::cout<<"******************Calculating TimeStep Is "<<mtimestep<<  "******************" <<std::endl;
  
    KRATOS_CATCH("")

}


//***************************************************************************
//***************************************************************************

void Calculate_Element_Deformation_Force_Add_To_Node()
{
      KRATOS_TRY
      
      /// Set to zero de RHS
      SetToZeroRHS();
      
      /// Compute the global external nodal force.
      Calculate_Conditions_RHS_and_Add();

      
      /// Compute the stress and body force of the element. ( No lineal analysis)
      Calculate_Elements_RHS_and_Add();

      
      KRATOS_CATCH("")

}

//***************************************************************************
//***************************************************************************

void SetToZeroRHS()
{
      KRATOS_TRY
      
      ModelPart& r_model_part  = BaseType::GetModelPart();
      NodesArrayType& pNodes   = r_model_part.Nodes(); 

      #ifdef _OPENMP
      int number_of_threads = omp_get_max_threads();
      #else
      int number_of_threads = 1;
      #endif

      vector<unsigned int> node_partition;
      CreatePartition(number_of_threads, pNodes.size(), node_partition);

      #pragma omp parallel for 
      for(int k=0; k<number_of_threads; k++)
	{
	  typename NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
	  typename NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];

          for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)      
	  {
	    array_1d<double,3>& normal    = (i->FastGetSolutionStepValue(NORMAL));
	    array_1d<double,3>& node_rhs  = (i->FastGetSolutionStepValue(RHS)); 
	    noalias(normal)               = ZeroVector(3);
	    noalias(node_rhs)             = ZeroVector(3);
	  }
	}
                             
     KRATOS_CATCH("")
}


//***************************************************************************
//***************************************************************************

void Calculate_Conditions_RHS_and_Add()
{
  KRATOS_TRY

      ModelPart& r_model_part          = BaseType::GetModelPart();   
      ProcessInfo& CurrentProcessInfo  = r_model_part.GetProcessInfo();   
      ConditionsArrayType& pConditions = r_model_part.Conditions();
      Vector rhs_cond;  
        
      #ifdef _OPENMP
      int number_of_threads = omp_get_max_threads();
      #else
      int number_of_threads = 1;
      #endif
      vector<unsigned int> condition_partition;
      CreatePartition(number_of_threads, pConditions.size(), condition_partition);
      unsigned int index; 
     #pragma omp parallel for private (index, rhs_cond)
      for(int k=0; k<number_of_threads; k++)
      {    
	typename ConditionsArrayType::iterator it_begin=pConditions.ptr_begin()+condition_partition[k];
	typename ConditionsArrayType::iterator it_end=pConditions.ptr_begin()+condition_partition[k+1];
  
        for (ConditionsArrayType::iterator it= it_begin; it!=it_end; ++it)
	{       
	  Condition::GeometryType& geom = it->GetGeometry();  
	  it->CalculateRightHandSide(rhs_cond,CurrentProcessInfo);    
	  const unsigned int& dim = geom.WorkingSpaceDimension();
          if(geom.size() > 1)
          {
             for (unsigned int i = 0; i <geom.size(); i++)
             {
	      index = i*dim;
	      array_1d<double,3>& node_rhs = geom(i)->FastGetSolutionStepValue(RHS);
 	      for(unsigned int kk=0; kk<dim; kk++)
	       {  geom(i)->SetLock();
	          node_rhs[kk] = node_rhs[kk] + rhs_cond[index+kk];
                  geom(i)->UnSetLock();
	       }
            } 
          }
        }
      }
 
	KRATOS_CATCH("")
}


//***************************************************************************
//***************************************************************************

        void Calculate_Elements_RHS_and_Add()
        {

            KRATOS_TRY
            ModelPart& r_model_part = BaseType::GetModelPart();
            ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
            ElementsArrayType& pElements = r_model_part.Elements();
          

            typename ElementsArrayType::iterator it_begin = pElements.ptr_begin() ;
            typename ElementsArrayType::iterator it_end   = pElements.ptr_end();
            for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it)
            {
                unsigned int index;
                Vector rhs_elem;

                Element::GeometryType& geom = it->GetGeometry();
                const unsigned int& dim = it->GetGeometry().WorkingSpaceDimension();
                it->CalculateRightHandSide(rhs_elem, CurrentProcessInfo);

                if (geom.size() > 1)
                {
                    for (unsigned int i = 0; i < geom.size(); i++)
                    {
                        index = i*dim;
                        array_1d<double, 3 > & node_rhs = geom(i)->FastGetSolutionStepValue(RHS);
                        for (unsigned int kk = 0; kk < dim; kk++)
                        {
                            geom(i)->SetLock();
                            node_rhs[kk] += /*node_rhs[kk]+*/ rhs_elem[index + kk];
                            geom(i)->UnSetLock();
                        }
                    }
                }
            }


            KRATOS_CATCH("")
        }





void FinalizeSolutionStep()
    {
    KRATOS_TRY
    ModelPart& r_model_part = BaseType::GetModelPart();  
    ElementsArrayType& pElements = r_model_part.Elements();
    ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
    ConditionsArrayType& pConditions = r_model_part.Conditions();

    #ifdef _OPENMP
    int number_of_threads = omp_get_max_threads();
    #else
    int number_of_threads = 1;
    #endif

    vector<unsigned int> element_partition;
    vector<unsigned int> condition_partition;
    CreatePartition(number_of_threads, pElements.size(), element_partition);
    CreatePartition(number_of_threads, pConditions.size(), condition_partition);
   

    #pragma omp parallel for
    for(int k=0; k<number_of_threads; k++)
    {
      typename ElementsArrayType::iterator it_begin=pElements.ptr_begin()+element_partition[k];
      typename ElementsArrayType::iterator it_end=pElements.ptr_begin()+element_partition[k+1];

      for (ElementsArrayType::iterator it= it_begin; it!=it_end; ++it)
	  {
	    (it) -> FinalizeSolutionStep(CurrentProcessInfo);
	  }
    }

    #pragma omp parallel for
    for(int k=0; k<number_of_threads; k++)
      {
	typename ConditionsArrayType::iterator it_begin=pConditions.ptr_begin()+condition_partition[k];
	typename ConditionsArrayType::iterator it_end=pConditions.ptr_begin()+condition_partition[k+1];
  
         for (ConditionsArrayType::iterator it= it_begin; it!=it_end; ++it)
          {
	    (it) -> FinalizeSolutionStep(CurrentProcessInfo);
	  }

    }
    
    KRATOS_CATCH("")

}



void CalculateReaction() 
{
    ModelPart& r_model_part = BaseType::GetModelPart();   
    NodesArrayType& pNodes  = r_model_part.Nodes(); 
    
    #ifdef _OPENMP
    int number_of_threads = omp_get_max_threads();
    #else
    int number_of_threads = 1;
    #endif
    vector<unsigned int> node_partition;
    CreatePartition(number_of_threads, pNodes.size(), node_partition);
    #pragma omp parallel for 
    for(int k=0; k<number_of_threads; k++)
      {
	typename NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
	typename NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];
	for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)
         //for(ModelPart::NodeIterator i = r_model_part.NodesBegin() ; i != r_model_part.NodesEnd() ; ++i)
         {
	    array_1d<double,3>& reaction                 = (i->FastGetSolutionStepValue(REACTION));
            const double& mass                           = (i->FastGetSolutionStepValue(NODAL_MASS));
            const array_1d<double,3>& rhs                = (i->FastGetSolutionStepValue(RHS));
            const array_1d<double,3>& acceleration       = (i->FastGetSolutionStepValue(ACCELERATION));
            const array_1d<double,3>& velocity           = (i->FastGetSolutionStepValue(VELOCITY));
            //array_1d<double,3> dif                 =  rhs; 

            if( (i->pGetDof(DISPLACEMENT_X))->IsFixed() == true)
 		{reaction[0] =  -rhs[0] + mass * acceleration[0] + mass * malpha_damp * velocity[0];}
                
            if( i->pGetDof(DISPLACEMENT_Y)->IsFixed() == true )
		{reaction[1] = -rhs[1] + mass * acceleration[1] + mass * malpha_damp * velocity[1];}
          
            if( i->HasDofFor(DISPLACEMENT_Z))
	        {
		if( i->pGetDof(DISPLACEMENT_Z)->IsFixed() == true )
		{reaction[2] = -rhs[2] + mass * acceleration[2] + mass * malpha_damp * velocity[2];}
                }
            }
      } 
 }


/*
  void FindParticleNeighbours()
   {
      KRATOS_TRY


     if(tempParticle.size() > 0)
     {
         ModelPart& r_model_part = BaseType::GetModelPart();
        typedef DiscreteElement   ParticleType;
        //typedef ElementsArrayType ParticleContainerType;

        typedef ModelPart::ElementsContainerType::ContainerType   ElementsContainerType;
        typedef ElementsContainerType ParticleContainerType;
        
        ProcessInfo& rCurrentProcessInfo = r_model_part.GetProcessInfo();
        ParticleContainerType& pElements = r_model_part.ElementsArray();

        if (mdimension == 2)
        {
             Neighbours_Calculator<2, ParticleType>::Search_Neighbours(pElements, rCurrentProcessInfo, mIfParticleInitialSearch);
        }
        else if (mdimension == 3)
        {
             Neighbours_Calculator<3, ParticleType>::Search_Neighbours(pElements,  rCurrentProcessInfo, mIfParticleInitialSearch);
        }
        mIfParticleInitialSearch = false;
     }
     

     KRATOS_CATCH("")
  }
*/


void FindParticleNeighbours()
   {
      KRATOS_TRY


     if(tempParticle.size() > 0)
     {
        ModelPart& r_model_part = BaseType::GetModelPart();
        ProcessInfo& rCurrentProcessInfo = r_model_part.GetProcessInfo();
        typedef Element  ParticleType;

        if (mdimension == 2)
        {
             Neighbours_Calculator<2, ParticleType>::Search_Neighbours(aux_list_of_particles, rCurrentProcessInfo, mIfParticleInitialSearch);
        }
        else if (mdimension == 3)
        {
             Neighbours_Calculator<3, ParticleType>::Search_Neighbours(aux_list_of_particles,  rCurrentProcessInfo, mIfParticleInitialSearch);
        }
        mIfParticleInitialSearch = false;
     }


     KRATOS_CATCH("")
  }


  void FindParticleBlockNeighbours()
   {
      KRATOS_TRY

     if(tempParticle.size() > 0 && aux_list_of_blocks.size() > 0)
     {

       unsigned int MaximumNumberOfResults = 100;
       ParticlePointerVector Results(MaximumNumberOfResults);

       typedef ParticleBlockConfigure < Element > ConfigureType;
       typedef BinsObjectStatic <ConfigureType> bins;
       bins particle_bin(aux_list_of_blocks.begin(), aux_list_of_blocks.end());



       //for(size_t ino = 0; ino < tempParticle.size(); ino++)
       for(ParticlePointerIterator ino = tempParticle.begin(); ino != tempParticle.end(); ino++)
       {
           ParticlePointer & rE = *ino;
         //ParticlePointer & rE = tempParticle[ino];

          typename ConfigureType::ResultIteratorType results_begin = Results.begin();
          
          int NumberofNeighbour = particle_bin.SearchObjects(*(ino.base()), results_begin, MaximumNumberOfResults) - 1;


          ParticleWeakVector   TempNeighbours;
          ParticleWeakVector & RealNeighbours = rE->GetValue(NEIGHBOUR_PARTICLE_BLOCK_ELEMENTS);
          TempNeighbours.swap(RealNeighbours);


          Vector TempContactForce;
          Vector & RealContactForce = rE->GetValue(PARTICLE_BLOCK_CONTACT_FORCE);
          TempContactForce.swap(RealContactForce);


          Vector TempContactFailureType;
          Vector & RealContactFailureType = rE->GetValue(PARTICLE_BLOCK_CONTACT_FAILURE_ID);
          TempContactFailureType.swap(RealContactFailureType);


          Vector TempIfInitialContact;
          Vector & RealIfInitialContact = rE->GetValue(PARTICLE_BLOCK_IF_INITIAL_CONTACT);
          TempIfInitialContact.swap(RealIfInitialContact);


          if(TempNeighbours.size() != TempContactFailureType.size())
          {
              KRATOS_WATCH("Neighbour size and Contact force size do not match!");
              KRATOS_WATCH(TempNeighbours.size());
              KRATOS_WATCH(TempContactFailureType.size());
              KRATOS_WATCH(TempContactForce.size());
              KRATOS_WATCH(TempIfInitialContact.size());

              KRATOS_WATCH(RealNeighbours.size());
              KRATOS_WATCH(RealContactFailureType.size());
              KRATOS_WATCH(RealContactForce.size());
              KRATOS_WATCH(RealIfInitialContact.size());
          }


          RealNeighbours.erase( RealNeighbours.begin(), RealNeighbours.end() );
          RealContactForce.clear();
          RealContactFailureType.clear();
          RealIfInitialContact.clear();

          int n_neighbours = NumberofNeighbour;
          int neighbour_counter = -1;

          // KRATOS_WATCH(NumberofNeighbour);

                for (ParticlePointerIterator neighbour_it = Results.begin(); neighbour_counter != n_neighbours; ++neighbour_it)
                {
                    //need not check, one is particle , the other one is block
                    //if (rE->GetGeometry()(0) != (*neighbour_it)->GetGeometry()(0))
                    {
                        
                        RealNeighbours.push_back(*neighbour_it);

                        size_t ForNo = RealContactForce.size();
                        RealContactForce.resize(ForNo + 3);
                        RealContactForce[ForNo + 0] = 0.0;
                        RealContactForce[ForNo + 1] = 0.0;
                        RealContactForce[ForNo + 2] = 0.0;

                        size_t FaiNo = RealContactFailureType.size();
                        RealContactFailureType.resize(FaiNo + 1);
                        RealContactFailureType[FaiNo] = 0;

                        RealIfInitialContact.resize(FaiNo + 1);
                        if(mIfParticleBlockInitialSearch == true)
                        {
                            RealIfInitialContact[FaiNo] = 1;
                        }
                        else
                        {
                            RealIfInitialContact[FaiNo] = 0;
                        }


                        int ContactForceNo = 0;
                        for (ParticleWeakIterator ineighbour = TempNeighbours.begin(); ineighbour != TempNeighbours.end(); ineighbour++)
                        {
                            if ((*neighbour_it)->GetGeometry()(0) == ineighbour->GetGeometry()(0))
                            {
                                RealContactForce[ForNo + 0] = TempContactForce[ContactForceNo * 3 + 0];
                                RealContactForce[ForNo + 1] = TempContactForce[ContactForceNo * 3 + 1];
                                RealContactForce[ForNo + 2] = TempContactForce[ContactForceNo * 3 + 2];

                                RealContactFailureType[FaiNo] = TempContactFailureType[ContactForceNo];

                                RealIfInitialContact[FaiNo]   = TempIfInitialContact[ContactForceNo];


                                break;
                            }

                            ContactForceNo++;
                        }
                    }
                    ++neighbour_counter;
                }


         
                int ContactForceNo = 0;
                for (ParticleWeakIterator ineighbour = TempNeighbours.begin(); ineighbour != TempNeighbours.end(); ineighbour++)
                {
                    if (TempContactFailureType[ContactForceNo] == 0)
                    {
                        
                        bool IfHaveSetUpRelation = false;
                        for(ParticleWeakIterator irightnei  = RealNeighbours.begin(); irightnei != RealNeighbours.end(); irightnei++)
                        {
                            if (irightnei->GetGeometry()(0) == ineighbour->GetGeometry()(0))
                            {
                                IfHaveSetUpRelation = true;
                                break;
                            }

                        }

                        if (IfHaveSetUpRelation == false)
                        {

                            RealNeighbours.push_back( TempNeighbours(ContactForceNo) );

                            size_t ForNo = RealContactForce.size();
                            RealContactForce.resize(ForNo + 3);
                            RealContactForce[ForNo + 0] = TempContactForce[ContactForceNo * 3 + 0];
                            RealContactForce[ForNo + 1] = TempContactForce[ContactForceNo * 3 + 1];
                            RealContactForce[ForNo + 2] = TempContactForce[ContactForceNo * 3 + 2];

                            size_t FaiNo = RealContactFailureType.size();
                            RealContactFailureType.resize(FaiNo + 1);
                            RealContactFailureType[FaiNo] = TempContactFailureType[ContactForceNo];

                            RealIfInitialContact.resize(FaiNo + 1);
                            RealIfInitialContact[FaiNo]   = TempIfInitialContact[ContactForceNo];

                            NumberofNeighbour++;
                        }
                    }

                    ContactForceNo++;
                }

            }

       }
     
       mIfParticleBlockInitialSearch = false;

       KRATOS_CATCH("")
    }



private:


unsigned int    mdimension;

bool   mElementsAreInitialized;
bool   mConditionsAreInitialized;
bool   mCalculateReactionsFlag;   
bool   mInitializeWasPerformed;
bool   mComputeTime;
double mdamping_ratio;
double malpha_damp;
double mbeta_damp; 
double mtimestep;  /// la suma de los delta time
double mpenalty_factor;

int    mdamp_type;
bool   mvirtual_mass;
double mcontact_stiffness_ratio;
bool   mComputeFemFemContact;
int    mnumber_step;

double mUnbalRatioCal;
bool   mIfHaveCalVirtualMass;
bool   mIfParticleInitialSearch;
bool   mIfParticleBlockInitialSearch;
bool   mIfHaveInitialzeSolutionStep;
ParticlePointerVector aux_list_of_particles;
ParticlePointerVector aux_list_of_blocks;
ParticlePointerVector tempParticle;


typename TBuilderAndSolverType::Pointer mpBuilderAndSolver;
typename TLinearSolver::Pointer mpLinearSolver;
typename TSchemeType::Pointer mpScheme;


//******************************************************************************************
//******************************************************************************************
inline void CreatePartition(unsigned int number_of_threads, const int number_of_rows, vector<unsigned int>& partitions)
    {
      partitions.resize(number_of_threads+1);
      int partition_size = number_of_rows / number_of_threads;
      partitions[0] = 0;
      partitions[number_of_threads] = number_of_rows;
      for(unsigned int i = 1; i<number_of_threads; i++)
      partitions[i] = partitions[i-1] + partition_size ;
  }


inline double Truncar_Delta_Time(double& num)
    {
      bool trunc = false;
      double num_trucado = num; 
      unsigned long int a     = 1;
      unsigned long int i     = 10;
      if(num!=0.00){
      while(trunc==false)
	{
	  num_trucado = num_trucado*i;
          a = a*i;
          if(num_trucado >= 1){
              num_trucado = static_cast<long unsigned int>(num_trucado);  
              num         = num_trucado/a;
              trunc       = true;
	   }
 	  }
	}
 
        return num;
 
    }

};

} /* namespace Kratos.*/
#endif /* KRATOS_RESIDUALBASED_CENTRAL_DIFERENCES_STRATEGY */


