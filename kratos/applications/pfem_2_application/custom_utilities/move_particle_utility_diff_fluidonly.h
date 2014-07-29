/*
==============================================================================
KratosIncompressibleFluidApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


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
//   Last Modified by:    $Author: pbecker $
//   Date:                $Date: 2011-09-21 12:30:32 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_MOVE_PART_UTILITY_FLUID_ONLY_DIFF2_INCLUDED )
#define  KRATOS_MOVE_PART_UTILITY_FLUID_ONLY_DIFF2_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <algorithm>

// External includes 


// Project includes
#include "includes/define.h"
#include "includes/node.h"

///
#include "includes/dof.h"
#include "includes/variables.h"
#include "containers/array_1d.h"
#include "containers/data_value_container.h"
#include "includes/mesh.h"
#include "utilities/math_utils.h"
#include "processes/node_erase_process.h" 
///

#include "utilities/geometry_utilities.h"

#include "includes/model_part.h"


#include "spatial_containers/spatial_containers.h"
#include "spatial_containers/bounding_box.h"
#include "spatial_containers/cell.h"
#include "spatial_containers/bins_dynamic_objects.h"

#include "utilities/spatial_containers_configure.h"

#include "geometries/line_2d_2.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/point.h"

#include "pfem_2_application.h"
#include "pfem_particle_fluidonly.h"

//#include "utilities/enrich_2d_2dofs.h"
#include "utilities/enrichment_utilities.h"
#include "utilities/openmp_utils.h"

#include "time.h"

//#include "processes/process.h"

namespace Kratos
{
	//this class is to be modified by the user to customize the interpolation process
	template< unsigned int TDim>
	class MoveParticleUtilityDiffFluidOnly 
	{
	public:
	
	    typedef SpatialContainersConfigure<TDim>     Configure;   
	    typedef typename Configure::PointType                      PointType; 
	    //typedef PointType::CoordinatesArrayType           CoordinatesArrayType;
        typedef typename Configure::ContainerType                  ContainerType;   
        //typedef Configure::PointerType                    PointerType;
        typedef typename Configure::IteratorType                   IteratorType; 
        typedef typename Configure::ResultContainerType            ResultContainerType;
	    //typedef Configure::ResultPointerType              ResultPointerType;
        typedef typename Configure::ResultIteratorType             ResultIteratorType; 
        typedef PointerVector< PFEM_Particle_Fluid, PFEM_Particle_Fluid*, std::vector<PFEM_Particle_Fluid*> > ParticlePointerVector;
        //typedef Configure::ContactPairType                ContactPairType;
        //typedef Configure::ContainerContactType           ContainerContactType; 
        //typedef Configure::IteratorContactType            IteratorContactType; 
        //typedef Configure::PointerContactType             PointerContactType; 
        //typedef Configure::PointerTypeIterator            PointerTypeIterator;

		KRATOS_CLASS_POINTER_DEFINITION(MoveParticleUtilityDiffFluidOnly);

		//template<unsigned int TDim>
		MoveParticleUtilityDiffFluidOnly(ModelPart& model_part, int maximum_number_of_particles)
			: mr_model_part(model_part) , mmaximum_number_of_particles(maximum_number_of_particles)  
		{
			KRATOS_WATCH("initializing moveparticle utility")
			//modd_timestep=false;
			mwrite_particle=0.0;
			mload_particle=0.0;
			mload_current_ParticlePointers=0.0;
			mload_old_ParticlePointers=0.0;
			mdo_nothing_per_particle=0.0;
			mintialized_transfer_tool=false;
			mcalculation_domain_complete_displacement=ZeroVector(3);
			mcalculation_domain_added_displacement=ZeroVector(3);
			
			//Condition const& rReferenceCondition = KratosComponents<Condition>::Get("Condition3D");      
			//KRATOS_WATCH("hola2")   
            //Properties::Pointer properties = mr_model_part.GetMesh().pGetProperties(1); 
            //KRATOS_WATCH("hola3aa")
            
            //storing water and air density and their inverses
            ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
			mDENSITY_AIR = CurrentProcessInfo[DENSITY_AIR];	
			mDENSITY_WATER = CurrentProcessInfo[DENSITY_WATER];	
			mINV_DENSITY_AIR = 1.0 / CurrentProcessInfo[DENSITY_AIR];	
			mINV_DENSITY_WATER = 1.0 / CurrentProcessInfo[DENSITY_WATER];	
			mDENSITY_AIRoverDENSITY_WATER = mDENSITY_AIR / mDENSITY_WATER;
			//mmaximum_number_of_particles = maximum_number_of_particles;
			
			
			//loop in elements to change their ID to their position in the array. Easier to get information later.
			//DO NOT PARALELIZE THIS! IT MUST BE SERIAL!!!!!!!!!!!!!!!!!!!!!!
			ModelPart::ElementsContainerType::iterator ielembegin = mr_model_part.ElementsBegin();
			for(unsigned int ii=0; ii<mr_model_part.Elements().size(); ii++)
			{
				ModelPart::ElementsContainerType::iterator ielem = ielembegin+ii;
				ielem->SetId(ii+1);
			}
			mlast_elem_id= (mr_model_part.ElementsEnd()-1)->Id();
            	
            // we look for the smallest edge. could be used as a weighting function when going lagrangian->eulerian instead of traditional shape functions(method currently used)
			ModelPart::NodesContainerType::iterator inodebegin = mr_model_part.NodesBegin();
			#pragma omp parallel for
			for(unsigned int ii=0; ii<mr_model_part.Nodes().size(); ii++)
			{
				ModelPart::NodesContainerType::iterator pnode = inodebegin+ii;
				array_1d<double,3> position_node;
				double distance=0.0;
				position_node = pnode->Coordinates();
				WeakPointerVector< Node<3> >& rneigh = pnode->GetValue(NEIGHBOUR_NODES);
				for( WeakPointerVector<Node<3> >::iterator inode = rneigh.begin(); inode!=rneigh.end(); inode++)
				{
					array_1d<double,3> position_difference;
					position_difference = inode->Coordinates() - position_node;
					double current_distance= sqrt(pow(position_difference[0],2)+pow(position_difference[1],2)+pow(position_difference[2],2));
					if (current_distance>distance)
						distance=current_distance;
				}
				pnode->FastGetSolutionStepValue(MEAN_SIZE)=distance;
			}
			
			//we also calculate the element mean size in the same way for the courant number
			//also we set the right size to the LHS vector
			#pragma omp parallel for
			for(unsigned int ii=0; ii<mr_model_part.Elements().size(); ii++)
			{
				ModelPart::ElementsContainerType::iterator ielem = ielembegin+ii;
				
				double mElemSize;
				array_1d<double,3> Edge(3,0.0);
				Edge = ielem->GetGeometry()[1].Coordinates() - ielem->GetGeometry()[0].Coordinates();
				mElemSize = Edge[0]*Edge[0];
				for (unsigned int d = 1; d < TDim; d++)
					mElemSize += Edge[d]*Edge[d];

				for (unsigned int i = 2; i < (TDim+1); i++)
					for(unsigned int j = 0; j < i; j++)
					{
						Edge = ielem->GetGeometry()[i].Coordinates() - ielem->GetGeometry()[j].Coordinates();
						double Length = Edge[0]*Edge[0];
						for (unsigned int d = 1; d < TDim; d++)
							Length += Edge[d]*Edge[d];
						if (Length < mElemSize) mElemSize = Length;
					}
				mElemSize = sqrt(mElemSize);
				
				ielem->GetValue(MEAN_SIZE) = mElemSize;
				if (TDim==3)
				{
					Vector & lhs_enrich = ielem->GetValue(ENRICH_LHS_ROW_3D);
					lhs_enrich.resize(4);
					lhs_enrich=ZeroVector(4);
				}
				//KRATOS_WATCH(mElemSize)
			}
			

            //matrix containing the position of the 4/15/45 particles that we will seed at the beggining
            boost::numeric::ublas::bounded_matrix<double, 5*(1+TDim), 3 > pos;
            boost::numeric::ublas::bounded_matrix<double, 5*(1+TDim), (1+TDim) > N;
            
            //int nnodes = mr_model_part.Nodes().size(); //at the beggining should be zero
            int particle_id=0;
			mnelems = mr_model_part.Elements().size();
			//mareas_vector.resize(nelems); UNUSED SO COMMENTED 
			KRATOS_WATCH("about to resize vectors")
			//setting the right size to the vector containing the particles assigned to each element			
			//particles vector. to be changed into pointers of particles
			mparticles_vector.resize(mnelems*mmaximum_number_of_particles);
			//ordered particles Ids, so that all particles that are inside an element are grouped togheter: first particles of element 1, then the ones inside element 2.. and so on
			//mparticles_ordered_pointers.resize(mnelems,mmaximum_number_of_particles*2);
			
			//and this vector contains the current number of particles that are in each element (currently zero)
			mnumber_of_particles_in_elems.resize(mnelems);
			mnumber_of_particles_in_elems=ZeroVector(mnelems);
			mnumber_of_particles_in_elems_aux.resize(mnelems);
			mpointers_to_particle_pointers_vectors.resize(mnelems);
			//new. water particles!
			mnumber_of_water_particles_in_elems.resize(mnelems);
			mnumber_of_water_particles_in_elems=ZeroVector(mnelems);

			int i_int=0; //careful! it's not the id, but the position inside the array!	
			KRATOS_WATCH("about to create particles")
			//now we seed: LOOP IN ELEMENTS
			//using loop index, can't paralelize! change lines : mparticles_in_elems_pointers((ii*mmaximum_number_of_particles)+mparticles_in_elems_integers(ii)) = pparticle; and the next one
			
			PFEM_Particle_Fluid& firstparticle =mparticles_vector[0];
			for(unsigned int ii=0; ii<mr_model_part.Elements().size(); ii++)
			{
				ModelPart::ElementsContainerType::iterator ielem = ielembegin+ii;
				//before doing anything we must reset the vector of nodes contained by each element
				///i->GetValue(NEIGHBOUR_NODES).resize(0); REPLACED BY THE vector: mparticles_in_elems
				//ielem->GetValue(FLUID_PARTICLE_POINTERS) = boost::shared_ptr< std::vector<PFEM_Particle_Fluid* > >( new std::vector<PFEM_Particle_Fluid* >() );	
				ParticlePointerVector&  particle_pointers =  (ielem->GetValue(FLUID_PARTICLE_POINTERS));
				mpointers_to_particle_pointers_vectors(ii) = &particle_pointers;
				//particle_pointers.resize(mmaximum_number_of_particles*2);
				for(unsigned int j=0; j<(mmaximum_number_of_particles*2); j++)
					particle_pointers.push_back(&firstparticle);
				int & number_of_particles = ielem->GetValue(NUMBER_OF_PARTICLES);
				int & number_of_water_particles = ielem->GetValue(NUMBER_OF_WATER_PARTICLES);
				
				Geometry< Node<3> >& geom = ielem->GetGeometry();
				//unsigned int elem_id = ielem->Id();
				//mareas_vector[i_int]=CalculateArea(geom); UNUSED SO COMMENTED 
				ComputeGaussPointPositions_initial(geom, pos, N); //we also have the standard (4), and 45
				//now we seed the particles in the current element
				for (unsigned int j = 0; j < pos.size1(); j++)
				{
					++particle_id;

					PFEM_Particle_Fluid& pparticle =mparticles_vector[particle_id-1];
					pparticle.X()=pos(j,0);
					pparticle.Y()=pos(j,1);
					pparticle.Z()=pos(j,2);
					
					pparticle.GetEraseFlag()=false;
					
					array_1d<double, 3 > & vel = pparticle.GetVelocity();
					//double & temp = pparticle.GetTemperature();
					double & distance= pparticle.GetDistance();
					//double & oxygen= pparticle.GetOxygen();
					//array_1d<double, 3 > & displacement= pparticle->FastGetSolutionStepValue(DISPLACEMENT);

					noalias(vel) = ZeroVector(3);
					//noalias(displacement) = ZeroVector(3);
					//temp = 0.0;
					distance=0.0;
					double pressure=0.0;
					//oxygen=0.0;
					//assigning values using shape functions (loop on nodes)	
                    for (unsigned int k = 0; k < (TDim+1); k++)
                    {
						noalias(vel) += N(j, k) * geom[k].FastGetSolutionStepValue(VELOCITY);
						//temp+= N(j, k) * geom[k].FastGetSolutionStepValue(TEMPERATURE);
						distance +=  N(j, k) * geom[k].FastGetSolutionStepValue(DISTANCE);
						pressure +=  N(j, k) * geom[k].FastGetSolutionStepValue(PRESSURE);

						//oxygen +=  N(j, k) * geom[k].FastGetSolutionStepValue(OXYGEN_FRACTION);
					}
					
					if(	ii % 100000 == 0)
						KRATOS_WATCH(particle_id);
					
					
						
					if (distance<=0.0)
					{
						//oxygen=0.0;
						distance=-1.0;

					}
					else 
					{
						//pparticle.GetEraseFlag()=true; //kill it!
						distance=1.0;

					}
					//if (distance>0.0)  pparticle.GetEraseFlag()=true;
					
                    //pparticle.GetElement() = Element::Pointer(*ielem.base());
					
					

						
					
					particle_pointers(j) = &pparticle;
					 number_of_particles++ ;
					if (distance<0)
						 number_of_water_particles++ ;
					

				}
				++i_int;
			}
			
			m_nparticles=particle_id; //we save the last particle created as the total number of particles we have. For the moment this is true.
			KRATOS_WATCH(m_nparticles);
			KRATOS_WATCH(mlast_elem_id);
		}
		

		~MoveParticleUtilityDiffFluidOnly()
		{}

		void MountBinDiff() //qué tipo de data es el object bin?
		{
			KRATOS_TRY

			//copy the elements to a new container, as the list will
			//be shuffled duringthe construction of the tree
			ContainerType& rElements           =  mr_model_part.ElementsArray();
	        IteratorType it_begin              =  rElements.begin();
	        IteratorType it_end                =  rElements.end();
	        //const int number_of_elem 		   =   rElements.size();

			typename BinsObjectDynamic<Configure>::Pointer paux = typename BinsObjectDynamic<Configure>::Pointer(new BinsObjectDynamic<Configure>(it_begin, it_end  ) );
			paux.swap(mpBinsObjectDynamic);
			//BinsObjectDynamic<Configure>  mpBinsObjectDynamic(it_begin, it_end ); 
				
			KRATOS_WATCH("inside MountBin")

			KRATOS_CATCH("")
		}
		
		
		//TOOL TO TRANSFER INFORMATION INITIALLY FROM ONE DOMAIN TO OTHER.
		void IntializeTransferTool(ModelPart::Pointer topographic_model_part, array_1d<double, 3 > initial_domains_offset, bool ovewrite_particle_data)
			//mtopographic_model_part(topographic_model_part)
		{
			KRATOS_TRY
			
			mintialized_transfer_tool=true;
			const unsigned int max_results = 1000;
			std::cout << "initializing transfer utility" << std::endl;
            ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
            mcalculation_domain_complete_displacement=initial_domains_offset;
            
            mtopographic_model_part_pointer =  topographic_model_part; //copying the pointer.
            
            //CONSTRUCTING BIN STRUCTURE
            ContainerType& rElements_topo           =  mtopographic_model_part_pointer->ElementsArray();
	        IteratorType it_begin_topo              =  rElements_topo.begin();
	        IteratorType it_end_topo                =  rElements_topo.end();
			typename BinsObjectDynamic<Configure>::Pointer paux = typename BinsObjectDynamic<Configure>::Pointer(new BinsObjectDynamic<Configure>(it_begin_topo, it_end_topo  ) );
			paux.swap(mpTopographicBinsObjectDynamic);
			

			std::cout << "Gathering Information From Topographic Domain for the first time" << std::endl;
			if(ovewrite_particle_data==false)
			{
				std::cout << "Not overwriting particle data (assuming correct initial conditions in calculation domain)" << std::endl;		  
			}
			else
			{
				std::cout << "Replacing particle information using the Topographic domain" << std::endl;		  
				const int offset = CurrentProcessInfo[PARTICLE_POINTERS_OFFSET]; //the array of pointers for each element has twice the required size so that we use a part in odd timesteps and the other in even ones.
				KRATOS_WATCH(offset)											//(flag managed only by MoveParticlesDiff
				ModelPart::ElementsContainerType::iterator ielembegin = mr_model_part.ElementsBegin();
				vector<unsigned int> element_partition;
				#ifdef _OPENMP
					int number_of_threads = omp_get_max_threads();
				#else
					int number_of_threads = 1;
				#endif
				OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Elements().size(), element_partition);
				
				//before doing anything we must reset the vector of nodes contained by each element (particles that are inside each element.
				#pragma omp parallel for
				for(int kkk=0; kkk<number_of_threads; kkk++)
				{
					ResultContainerType results(max_results);
					ResultIteratorType result_begin = results.begin();
					
					for(unsigned int ii=element_partition[kkk]; ii<element_partition[kkk+1]; ii++)
					{
						if (results.size()!=max_results)
							results.resize(max_results);
						//const int & elem_id = ielem->Id();
						ModelPart::ElementsContainerType::iterator ielem = ielembegin+ii;
						Element::Pointer pelement(*it_begin_topo.base());  //we have no idea in which element it might be from the topographic domain, so we just set it in the first element.
						
						//Geometry<Node<3> >& geom = ielem->GetGeometry(); 
						//array_1d<double,TDim+1> N;
						
						
						ParticlePointerVector&  element_particle_pointers =  (ielem->GetValue(FLUID_PARTICLE_POINTERS));
						int & number_of_particles_in_elem=ielem->GetValue(NUMBER_OF_PARTICLES);
						//std::cout << "elem " << ii << " with " << (unsigned int)number_of_particles_in_elem << " particles" << std::endl;
						
						for (unsigned int iii=0; iii<number_of_particles_in_elem ; iii++ )
						{
							//KRATOS_WATCH(iii)
							if (iii>mmaximum_number_of_particles) //it means we are out of our portion of the array, abort loop!
								break; 

							PFEM_Particle_Fluid & pparticle = element_particle_pointers[offset+iii];
							
							
							bool erase_flag= pparticle.GetEraseFlag();
							if (erase_flag==false)
							{
								OverwriteParticleDataUsingTopographicDomain(pparticle,pelement,mcalculation_domain_complete_displacement,result_begin, max_results); 
									
							}
							
							
						}
						
					}
				}
			}			
			KRATOS_CATCH("")

		}
		
		
		//TOOL TO TRANSFER INFORMATION FROM ONE DOMAIN TO OTHER when necessary. to be don
		void PreReseedUsingTopographicDomain(const int minimum_number_of_particles, array_1d<double, 3 > domains_added_displacement)
			//mtopographic_model_part(topographic_model_part)
		{
			KRATOS_TRY
			
			if(mintialized_transfer_tool==false)
				KRATOS_ERROR(std::logic_error, "TRANSFER TOOL NOT INITIALIZED!", "");
			const unsigned int max_results = 1000;
			std::cout << "executing transfer tool" << std::endl;
            ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
            mcalculation_domain_added_displacement = domains_added_displacement;
            mcalculation_domain_complete_displacement += domains_added_displacement;
             
            ContainerType& rElements_topo           =  mtopographic_model_part_pointer->ElementsArray();
	        IteratorType it_begin_topo              =  rElements_topo.begin();                        

			const int offset = CurrentProcessInfo[PARTICLE_POINTERS_OFFSET]; //the array of pointers for each element has twice the required size so that we use a part in odd timesteps and the other in even ones.
			KRATOS_WATCH(offset)											//(flag managed only by MoveParticlesDiff
			ModelPart::ElementsContainerType::iterator ielembegin = mr_model_part.ElementsBegin();
			vector<unsigned int> element_partition;
			#ifdef _OPENMP
				int number_of_threads = omp_get_max_threads();
			#else
				int number_of_threads = 1;
			#endif
			OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Elements().size(), element_partition);
			
			//before doing anything we must reset the vector of nodes contained by each element (particles that are inside each element.
			#pragma omp parallel for
			for(int kkk=0; kkk<number_of_threads; kkk++)
			{
				ResultContainerType results(max_results);
				ResultIteratorType result_begin = results.begin();
				
				Element::Pointer pelement(*it_begin_topo.base());  //we have no idea in which element it might be from the topographic domain, so we just set it in the first element.

				boost::numeric::ublas::bounded_matrix<double, (TDim+1), 3 > pos;
				boost::numeric::ublas::bounded_matrix<double, (TDim+1) , (TDim+1) > N;
				unsigned int freeparticle=0; //we start with the first position in the particles array
				
				for(unsigned int ii=element_partition[kkk]; ii<element_partition[kkk+1]; ii++)
				{
					if (results.size()!=max_results)
						results.resize(max_results);
					//const int & elem_id = ielem->Id();
					ModelPart::ElementsContainerType::iterator ielem = ielembegin+ii;
					
					ParticlePointerVector&  element_particle_pointers =  (ielem->GetValue(FLUID_PARTICLE_POINTERS));
					int & number_of_particles_in_elem=ielem->GetValue(NUMBER_OF_PARTICLES);
					if (number_of_particles_in_elem<(minimum_number_of_particles))// && (ielem->GetGeometry())[0].Y()<0.10 )
				    {
						//KRATOS_WATCH("elem with little particles")
						Geometry< Node<3> >& geom = ielem->GetGeometry();
						ComputeGaussPointPositionsForPreReseed(geom, pos, N);
						//double conductivity = ielem->GetProperties()[CONDUCTIVITY];
						//KRATOS_WATCH(conductivity);
						for (unsigned int j = 0; j < (pos.size1()); j++) //i am dropping the last one, the one in the middle of the element
						{
							bool keep_looking = true;
							while(keep_looking)
							{
								if (mparticles_vector[freeparticle].GetEraseFlag()==true)
								{
									#pragma omp critical
									{
										if (mparticles_vector[freeparticle].GetEraseFlag()==true)
										{
											mparticles_vector[freeparticle].GetEraseFlag()=false;
											keep_looking=false;
										}
									}
									if (keep_looking==false)
										break;
									/*	
									else if (freeparticle<(it_end_particle_model_part-1))
										freeparticle++;
									*/ 
									else
										freeparticle++;
										//break;
								}
								else
								{
									//if (freeparticle<(it_end_particle_model_part-1))
										freeparticle++;
									//else
										//break; //we finished the list and we couldnt find a free space
								}
							}

							PFEM_Particle_Fluid pparticle(pos(j,0),pos(j,1),pos(j,2));
							/*
							PFEM_Particle_Fluid & pparticle = mparticles_vector[freeparticle];
							pparticle.X() = pos(j,0);   
							pparticle.Y() = pos(j,1);    
							pparticle.Z() = pos(j,2);
							*/
							array_1d<double,TDim+1>aux2_N;
							bool is_found = CalculatePosition(geom,pos(j,0),pos(j,1),pos(j,2),aux2_N);
							if (is_found==false)
							{
								KRATOS_WATCH(aux2_N);
							}
							
							pparticle.GetEraseFlag()=false;
							OverwriteParticleDataUsingTopographicDomain(pparticle,pelement,mcalculation_domain_complete_displacement,result_begin, max_results); 
							 //and we copy it to the array:
							 mparticles_vector[freeparticle] =  pparticle;
							 element_particle_pointers(offset+number_of_particles_in_elem) = &mparticles_vector[freeparticle];
							number_of_particles_in_elem++;
							//KRATOS_WATCH(number_of_particles_in_elem);
							//KRATOS_WATCH(mparticles_vector[freeparticle])
							//KRATOS_WATCH(geom)
						
						  }
					  }
				  }
			}
			
			
			
			
			
			
			KRATOS_CATCH("")

		}
		
		void CalculateVelOverElemSize()
		{
			KRATOS_TRY
			
			//ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
			
			const double nodal_weight = 1.0/ (1.0 + double (TDim) );

			array_1d<double, 3 > vector_mean_velocity;
			
			ModelPart::ElementsContainerType::iterator ielembegin = mr_model_part.ElementsBegin();
			#pragma omp parallel for firstprivate(vector_mean_velocity)
			for(unsigned int ii=0; ii<mr_model_part.Elements().size(); ii++)
			{
				ModelPart::ElementsContainerType::iterator ielem = ielembegin+ii;
				Geometry<Node<3> >& geom = ielem->GetGeometry();

				vector_mean_velocity=ZeroVector(3);

				for (unsigned int i=0; i != (TDim+1) ; i++)
					vector_mean_velocity += geom[i].FastGetSolutionStepValue(VELOCITY);
				vector_mean_velocity *= nodal_weight;
				
				const double mean_velocity = sqrt ( pow(vector_mean_velocity[0],2) + pow(vector_mean_velocity[1],2) + pow(vector_mean_velocity[2],2) );
				ielem->GetValue(VELOCITY_OVER_ELEM_SIZE) = mean_velocity / ( ielem->GetValue(MEAN_SIZE) );
			}
			KRATOS_CATCH("")
		}
		
		

		//name self explained
		void ResetBoundaryConditions(bool fully_reset_nodes) 
		{
			KRATOS_TRY
			
			ModelPart::NodesContainerType::iterator inodebegin = mr_model_part.NodesBegin();
			if (fully_reset_nodes)
			{
				#pragma omp parallel for
				for(unsigned int ii=0; ii<mr_model_part.Nodes().size(); ii++)
				{
						ModelPart::NodesContainerType::iterator inode = inodebegin+ii;

						if (inode->IsFixed(VELOCITY_X))
						{
							inode->FastGetSolutionStepValue(VELOCITY_X)=inode->GetSolutionStepValue(VELOCITY_X,1);
						}
						if (inode->IsFixed(VELOCITY_Y))
						{
							inode->FastGetSolutionStepValue(VELOCITY_Y)=inode->GetSolutionStepValue(VELOCITY_Y,1);
						}
						if (TDim==3)
							if (inode->IsFixed(VELOCITY_Z))
							{
								inode->FastGetSolutionStepValue(VELOCITY_Z)=inode->GetSolutionStepValue(VELOCITY_Z,1);
							}
						
						if (inode->IsFixed(PRESSURE))
							inode->FastGetSolutionStepValue(PRESSURE)=inode->GetSolutionStepValue(PRESSURE,1);
						inode->GetSolutionStepValue(PRESSURE,1)=inode->FastGetSolutionStepValue(PRESSURE);		
				}
			}
			else  //for fractional step only!
			{
				#pragma omp parallel for
				for(unsigned int ii=0; ii<mr_model_part.Nodes().size(); ii++)
				{
						ModelPart::NodesContainerType::iterator inode = inodebegin+ii;

						const array_1d<double, 3 > original_velocity = inode->FastGetSolutionStepValue(VELOCITY);
						
						if (inode->IsFixed(VELOCITY_X)) // || inode->IsFixed(VELOCITY_Y) || inode->IsFixed(VELOCITY_Z))
						{
							const array_1d<double, 3 > & normal = inode->FastGetSolutionStepValue(NORMAL);
							const double normal_scalar_sq = normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2];
							const array_1d<double, 3 > normal_adimensionalized = normal / sqrt(normal_scalar_sq);
							array_1d<double, 3 > & velocity = inode->FastGetSolutionStepValue(VELOCITY);
							
							array_1d<double, 3 > normal_velocity;
							for (unsigned int j=0; j!=3; j++)
								normal_velocity[j] = fabs(normal_adimensionalized[j])*original_velocity[j];
							
							velocity[0] = normal_velocity[0]*inode->FastGetSolutionStepValue(NODAL_MASS)/inode->FastGetSolutionStepValue(NODAL_AREA);	
						}
						if (inode->IsFixed(VELOCITY_Y))
						{
							const array_1d<double, 3 > & normal = inode->FastGetSolutionStepValue(NORMAL);
							const double normal_scalar_sq = normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2];
							const array_1d<double, 3 > normal_adimensionalized = normal / sqrt(normal_scalar_sq);
							array_1d<double, 3 > & velocity = inode->FastGetSolutionStepValue(VELOCITY);
							
							array_1d<double, 3 > normal_velocity;
							for (unsigned int j=0; j!=3; j++)
								normal_velocity[j] = fabs(normal_adimensionalized[j])*original_velocity[j];
							
							velocity[1] = normal_velocity[1]*inode->FastGetSolutionStepValue(NODAL_MASS)/inode->FastGetSolutionStepValue(NODAL_AREA);	
							
						}
						if (TDim==3)
							if (inode->IsFixed(VELOCITY_Z))
							{
									const array_1d<double, 3 > & normal = inode->FastGetSolutionStepValue(NORMAL);
								const double normal_scalar_sq = normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2];
								const array_1d<double, 3 > normal_adimensionalized = normal / sqrt(normal_scalar_sq);
								array_1d<double, 3 > & velocity = inode->FastGetSolutionStepValue(VELOCITY);
								
								array_1d<double, 3 > normal_velocity;
								for (unsigned int j=0; j!=3; j++)
									normal_velocity[j] = fabs(normal_adimensionalized[j])*original_velocity[j];
								
								velocity[2] = normal_velocity[2]*inode->FastGetSolutionStepValue(NODAL_MASS)/inode->FastGetSolutionStepValue(NODAL_AREA);	
							}
						
						if (inode->IsFixed(PRESSURE))
							inode->FastGetSolutionStepValue(PRESSURE)=inode->GetSolutionStepValue(PRESSURE,1);
						inode->GetSolutionStepValue(PRESSURE,1)=inode->FastGetSolutionStepValue(PRESSURE);
						
				}
			}
			KRATOS_CATCH("")
		}
		
		void CalculateDeltaVelocity()//bool calculate_delta_pressure=false)
		{
			KRATOS_TRY
			//ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
			//const double one_over_delta_t = 1.0/CurrentProcessInfo[DELTA_TIME];
			
			//if we must do implicit correction, first we must subtract the RHS from the previoyus step
			ModelPart::NodesContainerType::iterator inodebegin = mr_model_part.NodesBegin();
			
			
				#pragma omp parallel for
				for(unsigned int ii=0; ii<mr_model_part.Nodes().size(); ii++)
				{
						ModelPart::NodesContainerType::iterator inode = inodebegin+ii;
						inode->FastGetSolutionStepValue(DELTA_VELOCITY) = inode->FastGetSolutionStepValue(VELOCITY) - inode->FastGetSolutionStepValue(MESH_VELOCITY) ;

				}

			KRATOS_CATCH("")
		}

		void FlagSplittedElementsAndTheirNodes() 
		{
			KRATOS_TRY
			ModelPart::NodesContainerType::iterator inodebegin = mr_model_part.NodesBegin();
			#pragma omp parallel for
			for(unsigned int ii=0; ii<mr_model_part.Nodes().size(); ii++)
			{
				ModelPart::NodesContainerType::iterator inode = inodebegin+ii;
				inode->FastGetSolutionStepValue(SPLIT_ELEMENT)=false;
			}
			
			ModelPart::ElementsContainerType::iterator ielembegin = mr_model_part.ElementsBegin();
			#pragma omp parallel for
			for(unsigned int ii=0; ii<mr_model_part.Elements().size(); ii++)
			{
				ModelPart::ElementsContainerType::iterator ielem = ielembegin+ii;
				Geometry<Node<3> >& geom = ielem->GetGeometry();
				bool is_split=false;
				double first_distance = geom[0].FastGetSolutionStepValue(DISTANCE);
				for (unsigned int i=1; i<(TDim+1);i++) //we do not check node zero
				{
					double & current_distance =  geom[i].FastGetSolutionStepValue(DISTANCE);
					if ((current_distance*first_distance)<0.0)
					{
						is_split=true;
						break;
					}
				}
				
				ielem->GetValue(SPLIT_ELEMENT)=is_split;
				
				if (is_split)
					for (unsigned int i=0; i<3;i++)
					{
						geom[i].SetLock();
						geom[i].FastGetSolutionStepValue(SPLIT_ELEMENT)=is_split;
						geom[i].UnSetLock();
					}
			}
			
			KRATOS_CATCH("")
		}
		
		
		

		//to move all the particles across the streamlines. heavy task!
		void MoveParticlesDiff(const bool viscosity_integrate, const bool add_gravity) //,const bool pressure_gradient_integrate) 
		{
			
			KRATOS_TRY
			
			ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
			
			const int offset = CurrentProcessInfo[PARTICLE_POINTERS_OFFSET]; //the array of pointers for each element has twice the required size so that we use a part in odd timesteps and the other in even ones.
																				//moveparticlesdiff reads from the pointers of one part (ie odd) and saves into the other part (ie even part)
																				//since it is the only function in the whole procedure that does this, it must use alternatively one part and the other.
			KRATOS_WATCH(offset)																		
																				
			bool even_timestep;
			if (offset!=0) even_timestep=false;
			else even_timestep=true;		
										
			const int post_offset = mmaximum_number_of_particles*int(even_timestep);	//and we also save the offset to know the location in which we will save the pointers after we've moved the particles																		
			KRATOS_WATCH(post_offset)	
			
			
			double delta_t = CurrentProcessInfo[DELTA_TIME];	
						
			array_1d<double,3> gravity= CurrentProcessInfo[GRAVITY];

			array_1d<double,TDim+1> N;
			const unsigned int max_results = 1000;
			
			const bool add_gravity_to_flying_particles = false; //!add_gravity; //when in most cases we do not want to add the gravity (in the fract vel), we still must do it for the "flying" particles, so for them we let the velocity be modified by the gravity.
			
			//double integration_distance= 2.0; 
			
			max_nsubsteps = 10;
			max_substep_dt=delta_t/double(max_nsubsteps);
			
			vector<unsigned int> element_partition;
			#ifdef _OPENMP
				int number_of_threads = omp_get_max_threads();
			#else
				int number_of_threads = 1;
			#endif
			OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Elements().size(), element_partition);

			ModelPart::ElementsContainerType::iterator ielembegin = mr_model_part.ElementsBegin();


			
			//before doing anything we must reset the vector of nodes contained by each element (particles that are inside each element.
			#pragma omp parallel for
			for(int kkk=0; kkk<number_of_threads; kkk++)
			{
				for(unsigned int ii=element_partition[kkk]; ii<element_partition[kkk+1]; ii++)
				{
					ModelPart::ElementsContainerType::iterator old_element = ielembegin+ii;
					
					int & number_of_particles = old_element->GetValue(NUMBER_OF_PARTICLES);
					
					mnumber_of_particles_in_elems_aux(ii)=number_of_particles;
					mnumber_of_particles_in_elems(ii)=0;
					mnumber_of_water_particles_in_elems(ii)=0;
					//we reset the local vectors for a faster access;
				}
			}
			KRATOS_WATCH("about to start transfering")
            //We move the particles across the fixed mesh and saving change data into them (using the function MoveParticle)
			//ModelPart::NodesContainerType::iterator iparticlebegin = mr_linea_model_part.NodesBegin();
			

			#pragma omp barrier
			
			#pragma omp parallel for
			for(int kkk=0; kkk<number_of_threads; kkk++)
			{
				
			  const array_1d<double,3> mesh_displacement = mcalculation_domain_added_displacement; //if it is a standard problem, displacements are zero and therefore nothing is added.
			  ResultContainerType results(max_results);
			
			  WeakPointerVector< Element > elements_in_trajectory;
			  elements_in_trajectory.resize(20);
				
			  for(unsigned int ielem=element_partition[kkk]; ielem<element_partition[kkk+1]; ielem++)
			  {
			//for(unsigned int ielem=0; ielem<mr_model_part.Elements().size(); ielem++)
			//{
				
				ModelPart::ElementsContainerType::iterator old_element = ielembegin+ielem;
				const int old_element_id = old_element->Id();

				ParticlePointerVector& old_element_particle_pointers = *mpointers_to_particle_pointers_vectors(old_element_id-1);
				
				if ( (results.size()) !=max_results)
						results.resize(max_results);
						
				
						
				unsigned int number_of_elements_in_trajectory=0; //excluding the origin one (current one, ielem)
				
				for(int ii=0; ii<(mnumber_of_particles_in_elems_aux(ielem)); ii++)
				{

					PFEM_Particle_Fluid & pparticle = old_element_particle_pointers[offset+ii];

					Element::Pointer pcurrent_element( *old_element.base() );
					ResultIteratorType result_begin = results.begin();
					bool & erase_flag=pparticle.GetEraseFlag();
					if (erase_flag==false){
						MoveParticle(viscosity_integrate,pparticle,pcurrent_element,elements_in_trajectory,number_of_elements_in_trajectory,result_begin,max_results, add_gravity_to_flying_particles,mesh_displacement); //saqué N de los argumentos, no lo necesito ya q empieza SIEMPRE en un nodo y no me importa donde termina

						const int current_element_id = pcurrent_element->Id();

						int & number_of_particles_in_current_elem = mnumber_of_particles_in_elems(current_element_id-1);
						int & number_of_water_particles_in_current_elem = mnumber_of_water_particles_in_elems(current_element_id-1);
						
						//original code with omp critical. too slow!
						if (number_of_particles_in_current_elem<mmaximum_number_of_particles && erase_flag==false) 
						{
							//if (old_element_id==current_element_id) //we can use the same vector of pointers that we gathered at the begining of the loop element
							//DO NOT USE, actually this means branching -> more expensive, it is better to leave it unbranched, for the worst case scenario.
							if (false)
							{
								#pragma omp critical
								{
									if (number_of_particles_in_current_elem<mmaximum_number_of_particles) // we cant go over this node, there's no room. otherwise we would be in the position of the first particle of the next element!!
									{	
										
										//mparticles_ordered_pointers((elem_id-1),(post_offset+number_of_particles_in_current_elem)) = &pparticle;
										
										old_element_particle_pointers(post_offset+number_of_particles_in_current_elem) = &pparticle;
										
										number_of_particles_in_current_elem++ ;
										if ((pparticle.GetDistance())<0.0)
											number_of_water_particles_in_current_elem++ ;
										
									}
									else 
										pparticle.GetEraseFlag()=true; //so we just delete it!
								}
							}
							else //a different element particle list must be used:
							{

								ParticlePointerVector& current_element_particle_pointers = *mpointers_to_particle_pointers_vectors(current_element_id-1);

								#pragma omp critical
								{
									if (number_of_particles_in_current_elem<mmaximum_number_of_particles) // we cant go over this node, there's no room. otherwise we would be in the position of the first particle of the next element!!
									{	

										current_element_particle_pointers(post_offset+number_of_particles_in_current_elem) = &pparticle;

										number_of_particles_in_current_elem++ ;
										if ((pparticle.GetDistance())<0.0)
											number_of_water_particles_in_current_elem++ ;
										if (number_of_particles_in_current_elem>mmaximum_number_of_particles)
											KRATOS_WATCH("MAL");
										
									}
									else 
										pparticle.GetEraseFlag()=true; //so we just delete it!
								}
							}
						}
						else 
							pparticle.GetEraseFlag()=true; //so we just delete it!
						

						
					}

					if (add_gravity)
						pparticle.GetVelocity() = pparticle.GetVelocity() + gravity*delta_t;
					

					
				}
			  }
			}
			
			//now we pass info from the local vector to the elements:
			#pragma omp parallel for
			for(unsigned int ii=0; ii<mr_model_part.Elements().size(); ii++)
			{
				ModelPart::ElementsContainerType::iterator old_element = ielembegin+ii;
				
				old_element->GetValue(NUMBER_OF_PARTICLES) = mnumber_of_particles_in_elems(ii);
				old_element->GetValue(NUMBER_OF_WATER_PARTICLES) = mnumber_of_water_particles_in_elems(ii);

			}
			
			
			//after having changed everything we change the status of the modd_timestep flag:
			CurrentProcessInfo[PARTICLE_POINTERS_OFFSET] = post_offset;; //

			KRATOS_CATCH("")
		}
		
		//void TransferLagrangianToEulerian(bool transfer_pressure, bool transfer_gradient_in_water_side)
		void TransferLagrangianToEulerian()
		{
			KRATOS_TRY

			ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
			//const double delta_t =CurrentProcessInfo[DELTA_TIME];
			const double threshold= 0.0/(double(TDim)+1.0);

				
			KRATOS_WATCH("About to Lagrangian->eulerian")
			
			const int offset = CurrentProcessInfo[PARTICLE_POINTERS_OFFSET]; //the array of pointers for each element has twice the required size so that we use a part in odd timesteps and the other in even ones.
			KRATOS_WATCH(offset)																	//(flag managed only by MoveParticlesDiff
			
			//we must project data from the particles (lagrangian)  into the eulerian mesh
			//ValuesVectorType eulerian_nodes_old_temperature;
			//int nnodes = mr_model_part.Nodes().size();
			//array_1d<double,(n_nodes)> eulerian_nodes_sumweights;
			
			//we save data from previous time step of the eulerian mesh in case we must reuse it later cos no particle was found around the nodes
			//though we could've use a bigger buffer, to be changed later!
			//after having saved data, we reset them to zero, this way it's easier to add the contribution of the surrounding particles.
			ModelPart::NodesContainerType::iterator inodebegin = mr_model_part.NodesBegin();
			#pragma omp parallel for
			for(unsigned int ii=0; ii<mr_model_part.Nodes().size(); ii++)
			{
				ModelPart::NodesContainerType::iterator inode = inodebegin+ii;
				inode->FastGetSolutionStepValue(DISTANCE)=0.0;
				inode->FastGetSolutionStepValue(MESH_VELOCITY)=ZeroVector(3); ///!!!!!!!!!!!!!!!!!!!!!!!!!!
				inode->FastGetSolutionStepValue(YP)=0.0;
				inode->FastGetSolutionStepValue(SOLID_YP)=1e-9;
	
			}
			
			//KRATOS_WATCH("ln 304")
			KRATOS_WATCH("About to Lagrangian->eulerian part 2")
			//adding contribution, loop on elements, since each element has stored the particles found inside of it
			ModelPart::ElementsContainerType::iterator ielembegin = mr_model_part.ElementsBegin();
			#pragma omp parallel for
			for(unsigned int ii=0; ii<mr_model_part.Elements().size(); ii++)
			{
				ModelPart::ElementsContainerType::iterator ielem = ielembegin+ii;
	
				array_1d<double,3*(TDim+1)> nodes_positions;
				array_1d<double,3*(TDim+1)> nodes_addedvel = ZeroVector(3*(TDim+1));
				
				array_1d<double,(TDim+1)> nodes_added_distance = ZeroVector((TDim+1));
				array_1d<double,(TDim+1)> nodes_addedweights = ZeroVector((TDim+1));
				array_1d<double,(TDim+1)> weighting_inverse_divisor;

				Geometry<Node<3> >& geom = ielem->GetGeometry();
				 
				for (int i=0 ; i!=(TDim+1) ; ++i) 
				{
					nodes_positions[i*3+0]=geom[i].X();
					nodes_positions[i*3+1]=geom[i].Y();
					nodes_positions[i*3+2]=geom[i].Z();
					weighting_inverse_divisor[i]=1.0/((geom[i].FastGetSolutionStepValue(MEAN_SIZE))*1.01); 
				}
				///KRATOS_WATCH(ielem->Id())
				///KRATOS_WATCH(ielem->GetValue(NEIGHBOUR_NODES).size());

				int & number_of_particles_in_elem= ielem->GetValue(NUMBER_OF_PARTICLES);
				ParticlePointerVector&  element_particle_pointers =  (ielem->GetValue(FLUID_PARTICLE_POINTERS));
				
				for (unsigned int iii=0; iii<number_of_particles_in_elem ; iii++ )
				{
					if (iii==mmaximum_number_of_particles) //it means we are out of our portion of the array, abort loop!
						break; 

					PFEM_Particle_Fluid & pparticle = element_particle_pointers[offset+iii];
					
					if (pparticle.GetEraseFlag()==false) 
					{
						
						array_1d<double,3> & position = pparticle.Coordinates();

						const array_1d<double,3>& velocity = pparticle.GetVelocity();

						const double& particle_distance = pparticle.GetDistance();  // -1 if water, +1 if air
					
						array_1d<double,TDim+1> N;
						bool is_found = CalculatePosition(geom,position[0],position[1],position[2],N);
						if (is_found==false) //something went wrong. if it was close enough to the edge we simply send it inside the element.
						{
							KRATOS_WATCH(N);
							for (int j=0 ; j!=(TDim+1); j++)
								if (N[j]<0.0 && N[j]> -1e-5)
									N[j]=1e-10;
						
						}
						
						for (int j=0 ; j!=(TDim+1); j++) //going through the 3/4 nodes of the element
						{
							double sq_dist = 0;
							//these lines for a weighting function based on the distance (or square distance) from the node insteadof the shape functions
							for (int k=0 ; k!=(TDim); k++) sq_dist += ((position[k] - nodes_positions[j*3+k])*(position[k] - nodes_positions[j*3+k]));
							double weight = (1.0 - (sqrt(sq_dist)*weighting_inverse_divisor[j] ) );
							
							weight=N(j);
							if (weight<threshold) weight=1e-10;
							if (weight<0.0) {KRATOS_WATCH(weight)}//;weight=0.0;KRATOS_WATCH(velocity);KRATOS_WATCH(N);KRATOS_WATCH(number_of_particles_in_elem);}//{KRATOS_WATCH(weight); KRATOS_WATCH(geom[j].Id()); KRATOS_WATCH(position);}
							else 
							{
								nodes_addedweights[j]+= weight;
								//nodes_addedtemp[j] += weight * particle_temp; 
								
								nodes_added_distance[j] += weight*particle_distance;
								
								//nodes_added_oxygen[j] += weight*particle_oxygen;

								for (int k=0 ; k!=(TDim); k++) //x,y,(z)
								{
									nodes_addedvel[j*3+k] += weight * velocity[k];
								}
								
							}//
						}
					}
				}
				
				for (int i=0 ; i!=(TDim+1) ; ++i) {
					geom[i].SetLock();
					geom[i].FastGetSolutionStepValue(DISTANCE) +=nodes_added_distance[i];
					geom[i].FastGetSolutionStepValue(MESH_VELOCITY_X) +=nodes_addedvel[3*i+0];
					geom[i].FastGetSolutionStepValue(MESH_VELOCITY_Y) +=nodes_addedvel[3*i+1];
					geom[i].FastGetSolutionStepValue(MESH_VELOCITY_Z) +=nodes_addedvel[3*i+2];  //we are updating info to the previous time step!!

					geom[i].FastGetSolutionStepValue(YP) +=nodes_addedweights[i];
					geom[i].UnSetLock();
				}
			}
			
			KRATOS_WATCH("ln 345")
			
			#pragma omp parallel for
			for(unsigned int ii=0; ii<mr_model_part.Nodes().size(); ii++)
			{
				ModelPart::NodesContainerType::iterator inode = inodebegin+ii;
				double sum_weights = inode->FastGetSolutionStepValue(YP);
				if (sum_weights>0.00001) 
				{
					//inode->FastGetSolutionStepValue(TEMPERATURE_OLD_IT)=(inode->FastGetSolutionStepValue(TEMPERATURE_OLD_IT))/sum_weights; //resetting the temperature
					double & dist = inode->FastGetSolutionStepValue(DISTANCE);
					 dist /=sum_weights; //resetting the density
					inode->FastGetSolutionStepValue(MESH_VELOCITY)=(inode->FastGetSolutionStepValue(MESH_VELOCITY))/sum_weights; //resetting the velocity
 					
				}
					
				else //this should never happen because other ways to recover the information have been executed before, but leaving it just in case..
				{
					inode->FastGetSolutionStepValue(DISTANCE)=inode->GetSolutionStepValue(DISTANCE,1); //resetting the temperature
					inode->FastGetSolutionStepValue(MESH_VELOCITY)=inode->GetSolutionStepValue(VELOCITY,1);

				}
				///finally, if there was an inlet that had a fixed position for the distance function, that has to remain unchanged:
				if (inode->IsFixed(DISTANCE))
					inode->FastGetSolutionStepValue(DISTANCE)=inode->GetSolutionStepValue(DISTANCE,1);
			}
			

			KRATOS_CATCH("")
		}
		
		

		void ReplaceParticlesVelocityAndDistance()
		{
			KRATOS_TRY
			
			//KRATOS_WATCH("ln 304")
			ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
			const int offset = CurrentProcessInfo[PARTICLE_POINTERS_OFFSET]; //the array of pointers for each element has twice the required size so that we use a part in odd timesteps and the other in even ones.
																				//(flag managed only by MoveParticlesDiff
			
			//adding contribution, loop on elements, since each element has stored the particles found inside of it
			ModelPart::ElementsContainerType::iterator ielembegin = mr_model_part.ElementsBegin();
			#pragma omp parallel for
			for(unsigned int ii=0; ii<mr_model_part.Elements().size(); ii++)
			{
				ModelPart::ElementsContainerType::iterator ielem = ielembegin+ii;
				array_1d<double,3*(TDim+1)> nodes_positions;
				array_1d<double,3*(TDim+1)> nodes_addedvel = ZeroVector(3*(TDim+1));
				//array_1d<double,(TDim+1)> nodes_addedtemp = ZeroVector((TDim+1));
				
				array_1d<double,(TDim+1)> nodes_added_distance = ZeroVector((TDim+1));
				//array_1d<double,(TDim+1)> nodes_added_is_water = ZeroVector((TDim+1));
				//array_1d<double,(TDim+1)> nodes_added_is_air = ZeroVector((TDim+1));
				array_1d<double,(TDim+1)> nodes_addedweights = ZeroVector((TDim+1));
				array_1d<double,(TDim+1)> weighting_inverse_divisor;
				//double mean_density=0.0;
				//array_1d<double,3*(TDim+1)> nodes_added_press_proj = ZeroVector(3*(TDim+1));
				//array_1d<double,(TDim+1)> nodes_added_pressure = ZeroVector((TDim+1));

				Geometry<Node<3> >& geom = ielem->GetGeometry();
				for (int i=0 ; i!=(TDim+1) ; ++i) {
					nodes_positions[i*3+0]=geom[i].X();
					nodes_positions[i*3+1]=geom[i].Y();
					nodes_positions[i*3+2]=geom[i].Z();
					weighting_inverse_divisor[i]=1.0/((geom[i].FastGetSolutionStepValue(MEAN_SIZE))*1.01); //ARGH!
				}
				///KRATOS_WATCH(ielem->Id())
				///KRATOS_WATCH(ielem->GetValue(NEIGHBOUR_NODES).size());
				//const int & elem_id = ielem->Id();

				ParticlePointerVector&  element_particle_pointers =  (ielem->GetValue(FLUID_PARTICLE_POINTERS));
				int & number_of_particles_in_elem=ielem->GetValue(NUMBER_OF_PARTICLES);


				for (unsigned int iii=0; iii<number_of_particles_in_elem ; iii++ )
				{
					if (iii>mmaximum_number_of_particles) //it means we are out of our portion of the array, abort loop!
						break; 

					PFEM_Particle_Fluid & pparticle = element_particle_pointers[offset+iii];
					
					if (pparticle.GetEraseFlag()==false) 
					{
						array_1d<double,3> & position = pparticle.Coordinates();
						//array_1d<double,3> velocity;
						//position[0] = iparticle->X()+iparticle->FastGetSolutionStepValue(DISPLACEMENT_X);
						//position[1] = iparticle->Y()+iparticle->FastGetSolutionStepValue(DISPLACEMENT_Y);
						//position[2] = iparticle->Z()+iparticle->FastGetSolutionStepValue(DISPLACEMENT_Z);

						array_1d<double,3>& velocity = pparticle.GetVelocity();
						//const array_1d<double,3>& press_proj = iparticle->FastGetSolutionStepValue(PRESS_PROJ);
						//const double& particle_pressure = iparticle->FastGetSolutionStepValue(PRESSURE);
						double& particle_distance = pparticle.GetDistance();  // -1 if water, +1 if air
						//velocity=ZeroVector(3);
						particle_distance=0.0;
						array_1d<double,TDim+1> N;
						CalculatePosition(geom,position[0],position[1],position[2],N);
						//KRATOS_WATCH(is_found)
						for (unsigned int i=0; i!=(TDim+1); i++) // 3/4 nodes of the element:
						{
							velocity+=N(i)*geom[i].FastGetSolutionStepValue(VELOCITY);
							particle_distance+=N(i)*geom[i].FastGetSolutionStepValue(DISTANCE);
						}
						//KRATOS_WATCH(distance)
					}
				}
			}
			
			
			
			KRATOS_CATCH("")
		}
		
		
		void AccelerateParticlesWithoutMovingUsingDeltaVelocity(const bool add_gravity, const bool update_stresses) 
		{
			KRATOS_TRY
			KRATOS_WATCH("accelerating particles using deltavelocity")
			ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
			double delta_t = CurrentProcessInfo[DELTA_TIME];	
			int fract_step_number = CurrentProcessInfo[FRACTIONAL_STEP];
			array_1d<double,3> & gravity= CurrentProcessInfo[GRAVITY];
			const bool update_stresses=true;
			const bool use_failure_criteria=true;
			
			//array_1d<double,TDim+1> N;
			//const int max_results = 1000;
			
			const int offset = CurrentProcessInfo[PARTICLE_POINTERS_OFFSET]; //the array of pointers for each element has twice the required size so that we use a part in odd timesteps and the other in even ones.
																				//(flag managed only by MoveParticlesDiff
			KRATOS_WATCH(offset)
			ModelPart::ElementsContainerType::iterator ielembegin = mr_model_part.ElementsBegin();
			
			
			vector<unsigned int> element_partition;
			#ifdef _OPENMP
				int number_of_threads = omp_get_max_threads();
			#else
				int number_of_threads = 1;
			#endif
			OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Elements().size(), element_partition);
			
			//before doing anything we must reset the vector of nodes contained by each element (particles that are inside each element.
			#pragma omp parallel for
			for(int kkk=0; kkk<number_of_threads; kkk++)
			{
				for(unsigned int ii=element_partition[kkk]; ii<element_partition[kkk+1]; ii++)
				{
					//const int & elem_id = ielem->Id();
					ModelPart::ElementsContainerType::iterator ielem = ielembegin+ii;
					Element::Pointer pelement(*ielem.base());
					Geometry<Node<3> >& geom = ielem->GetGeometry(); 
					
					ParticlePointerVector&  element_particle_pointers =  (ielem->GetValue(FLUID_PARTICLE_POINTERS));
					int & number_of_particles_in_elem=ielem->GetValue(NUMBER_OF_PARTICLES);
					//std::cout << "elem " << ii << " with " << (unsigned int)number_of_particles_in_elem << " particles" << std::endl;
					
					for (unsigned int iii=0; iii<number_of_particles_in_elem ; iii++ )
					{
						//KRATOS_WATCH(iii)
						if (iii>mmaximum_number_of_particles) //it means we are out of our portion of the array, abort loop!
							break; 

						PFEM_Particle_Fluid & pparticle = element_particle_pointers[offset+iii];
						
						
						bool erase_flag= pparticle.GetEraseFlag();
						if (erase_flag==false)
						{
							//Element::Pointer & pelement = pparticle.GetElement();
							AccelerateParticleUsingDeltaVelocity(pparticle,pelement,geom); //'lite' version, we pass by reference the geometry, so much cheaper
							if (fract_step_number==1 && add_gravity)
								pparticle.GetVelocity() = pparticle.GetVelocity() + gravity*delta_t;	
								
						}
						
						
					}
					
					
					Vector dummy_vector;
					if (update_stresses)
						ielem->Calculate(ELEMENT_MEAN_STRESS,dummy_vector,CurrentProcessInfo);
					
					//element_stresses = sum_particle_stresses/number_of_solid_particles;
					//KRATOS_WATCH(pparticle->FastGetSolutionStepValue(VELOCITY);

					
				}
			}
			KRATOS_CATCH("")
		}
		
		//**************************************************************************************************************
		//**************************************************************************************************************		

		
		template< class TDataType > void  AddUniqueWeakPointer
			(WeakPointerVector< TDataType >& v, const typename TDataType::WeakPointer candidate)
		{
			typename WeakPointerVector< TDataType >::iterator i = v.begin();
			typename WeakPointerVector< TDataType >::iterator endit = v.end();
			while ( i != endit && (i)->Id() != (candidate.lock())->Id())
			{
				i++;
			}
			if( i == endit )
			{
				v.push_back(candidate);
			}

		}    
		
		//**************************************************************************************************************
		//**************************************************************************************************************		
		
		void PreReseed(const bool viscosity_integrate, const bool add_gravity, int minimum_number_of_particles) 
		{
			KRATOS_TRY
			

			ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
			const int offset = CurrentProcessInfo[PARTICLE_POINTERS_OFFSET];
			double delta_t = CurrentProcessInfo[DELTA_TIME];	
			array_1d<double,3> & gravity= CurrentProcessInfo[GRAVITY];
			const int max_results = 1000;
			
			//tools for the paralelization
			//int last_id= (mr_linea_model_part.NodesEnd()-1)->Id();
			unsigned int number_of_threads = OpenMPUtils::GetNumThreads();
			//KRATOS_WATCH(number_of_threads);
			vector<unsigned int> elem_partition;
			int number_of_rows=mr_model_part.Elements().size();
			//KRATOS_WATCH(number_of_threads);
			//KRATOS_ERROR(std::logic_error, "Add  ----NODAL_H---- variable!!!!!! ERROR", "");
			elem_partition.resize(number_of_threads + 1);
			int elem_partition_size = number_of_rows / number_of_threads;
			elem_partition[0] = 0;
			elem_partition[number_of_threads] = number_of_rows;
			KRATOS_WATCH(elem_partition_size);
			for (unsigned int i = 1; i < number_of_threads; i++)
			elem_partition[i] = elem_partition[i - 1] + elem_partition_size;
			//typedef Node < 3 > PointType;
			//std::vector<PFEM_Particle_Fluid> aux;// aux;
			//aux.resize(number_of_threads);
			
			#pragma omp parallel firstprivate(elem_partition)
			{
				ResultContainerType results(max_results);
				int k = OpenMPUtils::ThisThread();
				ModelPart::ElementsContainerType::iterator it_begin = mr_model_part.ElementsBegin() +  elem_partition[k]; 
				ModelPart::ElementsContainerType::iterator it_end = mr_model_part.ElementsBegin() + elem_partition[k+1] ; 
				//ModelPart::NodesContainerType local_list=aux[k];
				//PointerVectorSet<PFEM_Particle_Fluid, IndexedObject> & list=aux[k];
				//KRATOS_WATCH(k);
			    boost::numeric::ublas::bounded_matrix<double, (TDim+1), 3 > pos;
				boost::numeric::ublas::bounded_matrix<double, (TDim+1) , (TDim+1) > N;
				unsigned int freeparticle=0; //we start with the first position in the particles array

				//int local_id=1;
				for (ModelPart::ElementsContainerType::iterator ielem = it_begin; ielem != it_end; ielem++)
				{
					results.resize(max_results);
					//const int & elem_id = ielem->Id();
					ParticlePointerVector&  element_particle_pointers =  (ielem->GetValue(FLUID_PARTICLE_POINTERS));
					int & number_of_particles_in_elem=ielem->GetValue(NUMBER_OF_PARTICLES);
					if (number_of_particles_in_elem<(minimum_number_of_particles))// && (ielem->GetGeometry())[0].Y()<0.10 )
				    {
						//KRATOS_WATCH("elem with little particles")
						Geometry< Node<3> >& geom = ielem->GetGeometry();
						ComputeGaussPointPositionsForPreReseed(geom, pos, N);
						//double conductivity = ielem->GetProperties()[CONDUCTIVITY];
						//KRATOS_WATCH(conductivity);
						for (unsigned int j = 0; j < (pos.size1()); j++) //i am dropping the last one, the one in the middle of the element
						{
							bool keep_looking = true;
							while(keep_looking)
							{
								if (mparticles_vector[freeparticle].GetEraseFlag()==true)
								{
									#pragma omp critical
									{
										if (mparticles_vector[freeparticle].GetEraseFlag()==true)
										{
											mparticles_vector[freeparticle].GetEraseFlag()=false;
											keep_looking=false;
										}
									}
									if (keep_looking==false)
										break;
									/*	
									else if (freeparticle<(it_end_particle_model_part-1))
										freeparticle++;
									*/ 
									else
										freeparticle++;
										//break;
								}
								else
								{
									//if (freeparticle<(it_end_particle_model_part-1))
										freeparticle++;
									//else
										//break; //we finished the list and we couldnt find a free space
								}
							}

							PFEM_Particle_Fluid pparticle(pos(j,0),pos(j,1),pos(j,2));
							/*
							PFEM_Particle_Fluid & pparticle = mparticles_vector[freeparticle];
							pparticle.X() = pos(j,0);   
							pparticle.Y() = pos(j,1);    
							pparticle.Z() = pos(j,2);
							*/
							array_1d<double,TDim+1>aux2_N;
							bool is_found = CalculatePosition(geom,pos(j,0),pos(j,1),pos(j,2),aux2_N);
							if (is_found==false)
							{
								KRATOS_WATCH(aux2_N);
							}
							
							//WeakPointerVector<Element > & pelement = (pparticle->GetValue(NEIGHBOUR_ELEMENTS));
							//pnode->GetSolutionStepValue(MEAN_SIZE)=1.0+0.01; //added a little to avoid problems when truncating with int(mean_size)
							pparticle.GetEraseFlag()=false;
							
						    ResultIteratorType result_begin = results.begin();
							//pparticle.GetElement()=(Element::Pointer( *ielem.base() ) );
							//Element::Pointer & pelement = pparticle.GetElement();
							Element::Pointer pelement( *ielem.base() );
							MoveParticle_inverse_way(viscosity_integrate, pparticle, pelement, result_begin, max_results);
							
							//pparticle.GetElement()=(Element::Pointer( *ielem.base() ) );
							
							//if (streamline_integrate)
							if (add_gravity)
								pparticle.GetVelocity() = pparticle.GetVelocity()+gravity*delta_t;///*double(nsubsteps);//only_integral; 

							
							 //and we copy it to the array:
							 mparticles_vector[freeparticle] =  pparticle;
							 
							 element_particle_pointers(offset+number_of_particles_in_elem) = &mparticles_vector[freeparticle];
							 pparticle.GetEraseFlag()=false;
							
							number_of_particles_in_elem++;
							//KRATOS_WATCH(number_of_particles_in_elem);
							//KRATOS_WATCH(mparticles_vector[freeparticle])
							//KRATOS_WATCH(geom)
						
						  }
					  }
				  }
			}



			
			KRATOS_CATCH("")
		}
		
		
		//**************************************************************************************************************
		//**************************************************************************************************************		
		
		
		void PostReseed(int minimum_number_of_particles, double mass_correction_factor ) //pooyan's way
		{
			KRATOS_TRY
			
			ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
			const int offset = CurrentProcessInfo[PARTICLE_POINTERS_OFFSET];
			
			if (mass_correction_factor>1.0) mass_correction_factor=1.0;
			if (mass_correction_factor<-1.0) mass_correction_factor=-1.0;
			//mass_correction_factor=0.0;
			
			//ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
			//const double delta_t = CurrentProcessInfo[DELTA_TIME];	
			//array_1d<double,3> & gravity= CurrentProcessInfo[GRAVITY];
			//const int max_results = 1000;
			
			const double threshold = mass_correction_factor*0.5;
			
			//TOOLS FOR THE PARALELIZATION
			//int last_id= (mr_linea_model_part.NodesEnd()-1)->Id();
			unsigned int number_of_threads = OpenMPUtils::GetNumThreads();
			//KRATOS_WATCH(number_of_threads);
			vector<unsigned int> elem_partition;
			int number_of_rows=mr_model_part.Elements().size();
			//KRATOS_WATCH(number_of_threads);
			//KRATOS_ERROR(std::logic_error, "Add  ----NODAL_H---- variable!!!!!! ERROR", "");
			elem_partition.resize(number_of_threads + 1);
			int elem_partition_size = number_of_rows / number_of_threads;
			elem_partition[0] = 0;
			elem_partition[number_of_threads] = number_of_rows;
			KRATOS_WATCH(elem_partition_size);
			for (unsigned int i = 1; i < number_of_threads; i++)
			elem_partition[i] = elem_partition[i - 1] + elem_partition_size;
			//typedef Node < 3 > PointType;
			//std::vector<ModelPart::NodesContainerType> aux;// aux;
			//aux.resize(number_of_threads);
			
			//ModelPart::NodesContainerType::iterator it_begin_particle_model_part = mr_linea_model_part.NodesBegin();
			//ModelPart::NodesContainerType::iterator it_end_particle_model_part = mr_linea_model_part.NodesEnd();

			#pragma omp parallel firstprivate(elem_partition) // firstprivate(results)//we will add the nodes in different parts of aux and later assemple everything toghether, remaming particles ids to get consecutive ids
			{
				unsigned int reused_particles=0;
				//unsigned int new_particles=0;
				
				//ModelPart::NodesContainerType::iterator freeparticle = it_begin_particle_model_part;
				unsigned int freeparticle = 0; //we start by the first position;
				//ResultContainerType results(max_results);
				//KRATOS_WATCH(elem_partition);
				
				int k = OpenMPUtils::ThisThread();
				//KRATOS_WATCH(k);
				ModelPart::ElementsContainerType::iterator it_begin = mr_model_part.ElementsBegin() +  elem_partition[k]; 
				ModelPart::ElementsContainerType::iterator it_end = mr_model_part.ElementsBegin() + elem_partition[k+1] ; 
				//ModelPart::NodesContainerType local_list=aux[k];
				//PointerVectorSet<Node<3>, IndexedObject> & list=aux[k];
				//KRATOS_WATCH(it_begin->Id());
				
				boost::numeric::ublas::bounded_matrix<double, (3+2*TDim), 3 > pos; //7 particles (2D) or 9 particles (3D)
				boost::numeric::ublas::bounded_matrix<double, (3+2*TDim), (TDim+1) > N;
				
				array_1d<double, 3 > vel_complete, vel_without_air_nodes;
				double sum_Ns_without_air_nodes;
				double mesh_distance;
				
				array_1d<double, (3+2*TDim) > distances;
				array_1d<int, (3+2*TDim) > positions;
				array_1d<bool, (3+2*TDim) > is_water_particle; //for both
				
				//bool has_water_node;
				//bool has_air_node;
				
				unsigned int number_of_reseeded_particles;
				//unsigned int number_of_water_reseeded_particles;

				//array_1d<double, 3 > nodes_distances;
				
				//int local_id=1;
				for (ModelPart::ElementsContainerType::iterator ielem = it_begin; ielem != it_end; ielem++)
				{
					//results.resize(max_results);
					
					int & number_of_particles_in_elem= ielem->GetValue(NUMBER_OF_PARTICLES);
					//int & number_of_water_particles_in_elem= ielem->GetValue(NUMBER_OF_WATER_PARTICLES);
					ParticlePointerVector&  element_particle_pointers =  (ielem->GetValue(FLUID_PARTICLE_POINTERS));

					
					
					//number_of_water_reseeded_particles=0;
					Geometry< Node<3> >& geom = ielem->GetGeometry();
					if ( (number_of_particles_in_elem<(minimum_number_of_particles)))// && (geom[0].Y()<0.10) ) || (number_of_water_particles_in_elem>2 && number_of_particles_in_elem<(minimum_number_of_particles) ) )
				    {
						
						//bool reseed_more=false;
						number_of_reseeded_particles=0;
						

						//reseed_more=true;
						number_of_reseeded_particles= 3+2*TDim;
						ComputeGaussPointPositionsForPostReseed(geom, pos, N);
						
						//if (number_of_water_particles_in_elem !=0 )
						//	number_of_water_reseeded_particles = double(number_of_reseeded_particles)*((mass_correction_factor+0.5) + double(number_of_water_particles_in_elem)) / double(number_of_particles_in_elem);
						

						
						distances = ZeroVector(3+2*TDim);
						
						//has_water_node=false;
						//has_air_node=false;
						/*
						for (unsigned int j = 0; j < (TDim+1); j++)
						{
							if ((geom[j].FastGetSolutionStepValue(DISTANCE))<0.0)
								has_water_node=true;
							else
								has_air_node=true;
						}
						*/
						for (unsigned int j = 0; j < number_of_reseeded_particles; j++) //first we order particles
						{
							positions[j]=j+1; //just creating a vector from 1 to 7 or whathever our lenght is (7 for 2d, 9 for 3d)
							for (unsigned int l = 0; l < (TDim+1); l++)
							{
								distances[j] +=  N(j, l) * geom[l].FastGetSolutionStepValue(DISTANCE);
							}
						}
						

							for (unsigned int j = 0; j < number_of_reseeded_particles ; j++) //first we order particles
							{
								if (distances[j]>threshold)
									is_water_particle[j]=false;
								else
									is_water_particle[j]=true;
							}

						bool fix_distance = false;
						unsigned int node_with_fixed_distance = 0;
						for (unsigned int j = 0; j < (TDim+1) ; j++) //we go over the 3/4 nodes:
						{
							if ((geom[j].IsFixed(DISTANCE)))
							{
								fix_distance = true;
								node_with_fixed_distance = j;
							}
						}
						// so now if the 3 were fixed, we assign the sign of the first node to all the particles:
						if (fix_distance)
						{
							bool is_water_for_all_particles=true;
							if ((geom[node_with_fixed_distance].FastGetSolutionStepValue(DISTANCE))>0.0)
								is_water_for_all_particles=false;
								
							for (unsigned int j = 0; j < number_of_reseeded_particles ; j++) //first we order particles
								is_water_particle[j]=is_water_for_all_particles;
						}
							
						
						for (unsigned int j = 0; j < number_of_reseeded_particles; j++) 
						{
							
							//Node < 3 > ::Pointer pnode = mr_linea_model_part.CreateNewNode(particle_id, pos(j,0), pos(j,1), pos(j,2));  //recordar que es el nueevo model part!!
							//now we have to find an empty space ( a particle that was about to be deleted) in the particles model part. once found. there will be our renewed particle:
							bool keep_looking = true;
							while(keep_looking)
							{
								if (mparticles_vector[freeparticle].GetEraseFlag()==true)
								{
									#pragma omp critical
									{
										if (mparticles_vector[freeparticle].GetEraseFlag()==true)
										{
											mparticles_vector[freeparticle].GetEraseFlag()=false;
											keep_looking=false;
										}
									}
									if (keep_looking==false)
										break;
									/*	
									else if (freeparticle<(it_end_particle_model_part-1))
										freeparticle++;
									*/ 
									else
										freeparticle++;
										//break;
								}
								else
								{
									//if (freeparticle<(it_end_particle_model_part-1))
										freeparticle++;
									//else
										//break; //we finished the list and we couldnt find a free space
								}
							}
							
							/*
							if (keep_looking) //it means we finished and we couldnt find a free space, so we need a completely new particle:
							{
								++local_id;
								pparticle = Node < 3 > ::Pointer (new Node < 3 >(local_id, pos(j,0), pos(j,1), pos(j,2)));  //recordar que es el nueevo model part!!
							}
							else //good luck!, no need to create!
							*/ 
							PFEM_Particle_Fluid pparticle(pos(j,0),pos(j,1),pos(j,2));
							

							//pnode->GetSolutionStepValue(CONDUCTIVITY)=conductivity;
							array_1d<double, 3 > & vel = pparticle.GetVelocity();
							//double & temp= pparticle.GetTemperature();
							double & distance= pparticle.GetDistance();
							//double & oxygen = pparticle.GetOxygen();
							//pnode->GetSolutionStepValue(DENSITY)=density;
							//array_1d<double, 3 > & disp = pnode->FastGetSolutionStepValue(DISPLACEMENT);
							//noalias(disp) = ZeroVector(3);
							
							array_1d<double,TDim+1>aux_N;
							bool is_found = CalculatePosition(geom,pos(j,0),pos(j,1),pos(j,2),aux_N);
							if (is_found==false)
							{
								KRATOS_WATCH(aux_N);
								KRATOS_WATCH(j)
								KRATOS_WATCH(ielem->Id())
							}
							
					
							noalias(vel_complete)=ZeroVector(3);
							noalias(vel_without_air_nodes)=ZeroVector(3);
							sum_Ns_without_air_nodes=0.0;
							
							noalias(vel) = ZeroVector(3);
							//temp = 0.0;
							distance=0.0;
							mesh_distance = 0.0;
							//oxygen = 0.0;
							
							double pressure=0.0;
							double solid_pressure=0.0;
							for (unsigned int l = 0; l < (TDim+1); l++)
							{
								//noalias(vel) += N(j, l) * geom[l].FastGetSolutionStepValue(VELOCITY);
								noalias(vel_complete) += N(j, l) * geom[l].FastGetSolutionStepValue(VELOCITY);
								pressure += N(j, l) * geom[l].FastGetSolutionStepValue(PRESSURE);
								solid_pressure += N(j, l) * geom[l].FastGetSolutionStepValue(SOLID_PRESSURE);
								//temp+= N(j, l) * geom[l].FastGetSolutionStepValue(TEMPERATURE);
								//oxygen +=  N(j, l) * geom[l].FastGetSolutionStepValue(OXYGEN_FRACTION);
								mesh_distance +=  N(j,l) * geom[l].FastGetSolutionStepValue(DISTANCE);
								if ((geom[l].FastGetSolutionStepValue(DISTANCE))<0.0)
								{
									sum_Ns_without_air_nodes+=N(j, l);
									noalias(vel_without_air_nodes) += N(j, l) * geom[l].FastGetSolutionStepValue(VELOCITY);
								}
							}
							
							///COMMENT TO GET A CONTINOUS DISTANCE FUNCTION FIELD
							if (is_water_particle[j])
							{
								distance=-1.0;
								//oxygen=0.0;
							}
							else
							{   
								if (mesh_distance<1.5)
									distance=1.0;
								else
									distance=2.0;
							}
								
							//if (distance<0.0 && sum_Ns_without_air_nodes>0.01)
							//	vel = vel_without_air_nodes / sum_Ns_without_air_nodes ;
							//else
								vel = vel_complete;
							
							//test:
							/*
							if (j==2 && add_2_water_node==true)
								distance=-1.0;
							if (j==3 && add_1_water_node==true)
								distance=-1.0;
							*/
							
							
							//pnode->GetValue(NEIGHBOUR_ELEMENTS)= Element::WeakPointer( *i.base() ) ;
							//pparticle->GetValue(NEIGHBOUR_ELEMENTS).resize(0);
							//pparticle.GetElement()=(Element::Pointer( *ielem.base() ) );
							//WeakPointerVector<Element > & pelement = (pparticle->GetValue(NEIGHBOUR_ELEMENTS));
							//pnode->GetSolutionStepValue(MEAN_SIZE)=1.0+0.01; //added a little to avoid problems when truncating with int(mean_size)
							pparticle.GetEraseFlag()=false;
							
							if (distance<=0.0)
							{
								distance=-1.0;
							}
							else //distance>1.5
							{
								distance=1.0;	
							}
							
							//pnode->GetSolutionStepValue(DISPLACEMENT_X)=0.0;
							//pnode->GetSolutionStepValue(DISPLACEMENT_Y)=0.0;
							//pnode->GetSolutionStepValue(DISPLACEMENT_Z)=0.0;
							//array_1d<double, 3 > & disp = pparticle->FastGetSolutionStepValue(DISPLACEMENT);
							//noalias(disp) = ZeroVector(3);

							//AddUniqueWeakPointer< Node<3> >(ielem->GetValue(NEIGHBOUR_NODES), pparticle);
							//if(distance<1.5)
							//{
								mparticles_vector[freeparticle]=pparticle;
								element_particle_pointers(offset+number_of_particles_in_elem) = &mparticles_vector[freeparticle];
								number_of_particles_in_elem++;
							//}
							//KRATOS_WATCH(number_of_particles_in_elem)
							if (keep_looking)
							{
								KRATOS_ERROR(std::logic_error, "FINISHED THE LIST AND COULDNT FIND A FREE CELL FOR THE NEW PARTICLE!", "");
							}
						    else
						    {
								reused_particles++;
								//KRATOS_WATCH(pparticle->FastGetSolutionStepValue(DISTANCE))
							}
						
						  }
					  }
				  }
				  //KRATOS_WATCH(reused_particles)
				  //KRATOS_WATCH(new_particles)
			}

			KRATOS_CATCH("")
		}
		
		
		void PostReseedOnlyInBoundingBox(int minimum_number_of_particles, double mass_correction_factor, const array_1d<double, 3 > bounding_box_lower_corner_full, const array_1d<double, 3 > bounding_box_upper_corner_full)
		{
			KRATOS_TRY
			
			ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
			const int offset = CurrentProcessInfo[PARTICLE_POINTERS_OFFSET];
			
			array_1d<double, TDim > bounding_box_lower_corner, bounding_box_upper_corner;
			for (unsigned int i = 0; i < TDim; i++)
			{
				bounding_box_lower_corner(i)=bounding_box_lower_corner_full(i);
				bounding_box_upper_corner(i)=bounding_box_upper_corner_full(i);
			}
			
			if (mass_correction_factor>1.0) mass_correction_factor=1.0;
			if (mass_correction_factor<-1.0) mass_correction_factor=-1.0;
			//mass_correction_factor=0.0;
			
			//ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
			//const double delta_t = CurrentProcessInfo[DELTA_TIME];	
			//array_1d<double,3> & gravity= CurrentProcessInfo[GRAVITY];
			//const int max_results = 1000;
			
			const double threshold = mass_correction_factor*0.5;
			
			//TOOLS FOR THE PARALELIZATION
			//int last_id= (mr_linea_model_part.NodesEnd()-1)->Id();
			unsigned int number_of_threads = OpenMPUtils::GetNumThreads();
			//KRATOS_WATCH(number_of_threads);
			vector<unsigned int> elem_partition;
			int number_of_rows=mr_model_part.Elements().size();
			//KRATOS_WATCH(number_of_threads);
			//KRATOS_ERROR(std::logic_error, "Add  ----NODAL_H---- variable!!!!!! ERROR", "");
			elem_partition.resize(number_of_threads + 1);
			int elem_partition_size = number_of_rows / number_of_threads;
			elem_partition[0] = 0;
			elem_partition[number_of_threads] = number_of_rows;
			KRATOS_WATCH(elem_partition_size);
			for (unsigned int i = 1; i < number_of_threads; i++)
			elem_partition[i] = elem_partition[i - 1] + elem_partition_size;
			//typedef Node < 3 > PointType;
			//std::vector<ModelPart::NodesContainerType> aux;// aux;
			//aux.resize(number_of_threads);
			
			//ModelPart::NodesContainerType::iterator it_begin_particle_model_part = mr_linea_model_part.NodesBegin();
			//ModelPart::NodesContainerType::iterator it_end_particle_model_part = mr_linea_model_part.NodesEnd();

			#pragma omp parallel firstprivate(elem_partition) // firstprivate(results)//we will add the nodes in different parts of aux and later assemple everything toghether, remaming particles ids to get consecutive ids
			{
				unsigned int reused_particles=0;
				//unsigned int new_particles=0;
				
				//ModelPart::NodesContainerType::iterator freeparticle = it_begin_particle_model_part;
				unsigned int freeparticle = 0; //we start by the first position;
				//ResultContainerType results(max_results);
				//KRATOS_WATCH(elem_partition);
				
				int k = OpenMPUtils::ThisThread();
				//KRATOS_WATCH(k);
				ModelPart::ElementsContainerType::iterator it_begin = mr_model_part.ElementsBegin() +  elem_partition[k]; 
				ModelPart::ElementsContainerType::iterator it_end = mr_model_part.ElementsBegin() + elem_partition[k+1] ; 
				//ModelPart::NodesContainerType local_list=aux[k];
				//PointerVectorSet<Node<3>, IndexedObject> & list=aux[k];
				//KRATOS_WATCH(it_begin->Id());
				
				boost::numeric::ublas::bounded_matrix<double, (3+2*TDim), 3 > pos; //7 particles (2D) or 9 particles (3D)
				boost::numeric::ublas::bounded_matrix<double, (3+2*TDim), (TDim+1) > N;
				
				array_1d<double, 3 > vel_complete, vel_without_air_nodes;
				double sum_Ns_without_air_nodes;
				double mesh_distance;
				
				array_1d<double, (3+2*TDim) > distances;
				array_1d<int, (3+2*TDim) > positions;
				array_1d<bool, (3+2*TDim) > is_water_particle; //for both
				
				//bool has_water_node;
				//bool has_air_node;
				
				unsigned int number_of_reseeded_particles;
				//unsigned int number_of_water_reseeded_particles;

				//array_1d<double, 3 > nodes_distances;
				
				//int local_id=1;
				for (ModelPart::ElementsContainerType::iterator ielem = it_begin; ielem != it_end; ielem++)
				{
					//results.resize(max_results);
					
					int & number_of_particles_in_elem= ielem->GetValue(NUMBER_OF_PARTICLES);
					//int & number_of_water_particles_in_elem= ielem->GetValue(NUMBER_OF_WATER_PARTICLES);

					
					
					//number_of_water_reseeded_particles=0;
					if ( (number_of_particles_in_elem<(minimum_number_of_particles)))// && (geom[0].Y()<0.10) ) || (number_of_water_particles_in_elem>2 && number_of_particles_in_elem<(minimum_number_of_particles) ) )
				    {
						Geometry< Node<3> >& geom = ielem->GetGeometry();
						bool elem_in_bounding_box=ChechIfElemIsInBoundingBox(geom,bounding_box_lower_corner,bounding_box_upper_corner);
						if(elem_in_bounding_box) //otherwise information might be coming from outside the domain, so we must use the ttopographic domain.
						{
							ParticlePointerVector&  element_particle_pointers =  (ielem->GetValue(FLUID_PARTICLE_POINTERS));
							
							//bool reseed_more=false;
							number_of_reseeded_particles=0;
							

							//reseed_more=true;
							number_of_reseeded_particles= 3+2*TDim;
							ComputeGaussPointPositionsForPostReseed(geom, pos, N);
							
							//if (number_of_water_particles_in_elem !=0 )
							//	number_of_water_reseeded_particles = double(number_of_reseeded_particles)*((mass_correction_factor+0.5) + double(number_of_water_particles_in_elem)) / double(number_of_particles_in_elem);
							

							
							distances = ZeroVector(3+2*TDim);
							
							//has_water_node=false;
							//has_air_node=false;
							/*
							for (unsigned int j = 0; j < (TDim+1); j++)
							{
								if ((geom[j].FastGetSolutionStepValue(DISTANCE))<0.0)
									has_water_node=true;
								else
									has_air_node=true;
							}
							*/
							for (unsigned int j = 0; j < number_of_reseeded_particles; j++) //first we order particles
							{
								positions[j]=j+1; //just creating a vector from 1 to 7 or whathever our lenght is (7 for 2d, 9 for 3d)
								for (unsigned int l = 0; l < (TDim+1); l++)
								{
									distances[j] +=  N(j, l) * geom[l].FastGetSolutionStepValue(DISTANCE);
								}
							}
							

								for (unsigned int j = 0; j < number_of_reseeded_particles ; j++) //first we order particles
								{
									if (distances[j]>threshold)
										is_water_particle[j]=false;
									else
										is_water_particle[j]=true;
								}

							bool fix_distance = false;
							unsigned int node_with_fixed_distance = 0;
							for (unsigned int j = 0; j < (TDim+1) ; j++) //we go over the 3/4 nodes:
							{
								if ((geom[j].IsFixed(DISTANCE)))
								{
									fix_distance = true;
									node_with_fixed_distance = j;
								}
							}
							// so now if the 3 were fixed, we assign the sign of the first node to all the particles:
							if (fix_distance)
							{
								bool is_water_for_all_particles=true;
								if ((geom[node_with_fixed_distance].FastGetSolutionStepValue(DISTANCE))>0.0)
									is_water_for_all_particles=false;
									
								for (unsigned int j = 0; j < number_of_reseeded_particles ; j++) //first we order particles
									is_water_particle[j]=is_water_for_all_particles;
							}
								
							
							for (unsigned int j = 0; j < number_of_reseeded_particles; j++) 
							{
								
								//Node < 3 > ::Pointer pnode = mr_linea_model_part.CreateNewNode(particle_id, pos(j,0), pos(j,1), pos(j,2));  //recordar que es el nueevo model part!!
								//now we have to find an empty space ( a particle that was about to be deleted) in the particles model part. once found. there will be our renewed particle:
								bool keep_looking = true;
								while(keep_looking)
								{
									if (mparticles_vector[freeparticle].GetEraseFlag()==true)
									{
										#pragma omp critical
										{
											if (mparticles_vector[freeparticle].GetEraseFlag()==true)
											{
												mparticles_vector[freeparticle].GetEraseFlag()=false;
												keep_looking=false;
											}
										}
										if (keep_looking==false)
											break;
										/*	
										else if (freeparticle<(it_end_particle_model_part-1))
											freeparticle++;
										*/ 
										else
											freeparticle++;
											//break;
									}
									else
									{
										//if (freeparticle<(it_end_particle_model_part-1))
											freeparticle++;
										//else
											//break; //we finished the list and we couldnt find a free space
									}
								}
								
								/*
								if (keep_looking) //it means we finished and we couldnt find a free space, so we need a completely new particle:
								{
									++local_id;
									pparticle = Node < 3 > ::Pointer (new Node < 3 >(local_id, pos(j,0), pos(j,1), pos(j,2)));  //recordar que es el nueevo model part!!
								}
								else //good luck!, no need to create!
								*/ 
								PFEM_Particle_Fluid pparticle(pos(j,0),pos(j,1),pos(j,2));
								

								//pnode->GetSolutionStepValue(CONDUCTIVITY)=conductivity;
								array_1d<double, 3 > & vel = pparticle.GetVelocity();
								//double & temp= pparticle.GetTemperature();
								double & distance= pparticle.GetDistance();
								//double & oxygen = pparticle.GetOxygen();
								//pnode->GetSolutionStepValue(DENSITY)=density;
								//array_1d<double, 3 > & disp = pnode->FastGetSolutionStepValue(DISPLACEMENT);
								//noalias(disp) = ZeroVector(3);
								
								array_1d<double,TDim+1>aux_N;
								bool is_found = CalculatePosition(geom,pos(j,0),pos(j,1),pos(j,2),aux_N);
								if (is_found==false)
								{
									KRATOS_WATCH(aux_N);
									KRATOS_WATCH(j)
									KRATOS_WATCH(ielem->Id())
								}
								
						
								noalias(vel_complete)=ZeroVector(3);
								noalias(vel_without_air_nodes)=ZeroVector(3);
								sum_Ns_without_air_nodes=0.0;
								
								noalias(vel) = ZeroVector(3);
								//temp = 0.0;
								distance=0.0;
								mesh_distance = 0.0;
								//oxygen = 0.0;
								
								double pressure=0.0;
								double solid_pressure=0.0;
								for (unsigned int l = 0; l < (TDim+1); l++)
								{
									//noalias(vel) += N(j, l) * geom[l].FastGetSolutionStepValue(VELOCITY);
									noalias(vel_complete) += N(j, l) * geom[l].FastGetSolutionStepValue(VELOCITY);
									pressure += N(j, l) * geom[l].FastGetSolutionStepValue(PRESSURE);
									solid_pressure += N(j, l) * geom[l].FastGetSolutionStepValue(SOLID_PRESSURE);
									//temp+= N(j, l) * geom[l].FastGetSolutionStepValue(TEMPERATURE);
									//oxygen +=  N(j, l) * geom[l].FastGetSolutionStepValue(OXYGEN_FRACTION);
									mesh_distance +=  N(j,l) * geom[l].FastGetSolutionStepValue(DISTANCE);
									if ((geom[l].FastGetSolutionStepValue(DISTANCE))<0.0)
									{
										sum_Ns_without_air_nodes+=N(j, l);
										noalias(vel_without_air_nodes) += N(j, l) * geom[l].FastGetSolutionStepValue(VELOCITY);
									}
								}
								
								///COMMENT TO GET A CONTINOUS DISTANCE FUNCTION FIELD
								if (is_water_particle[j])
								{
									distance=-1.0;
									//oxygen=0.0;
								}
								else
								{   
									if (mesh_distance<1.5)
										distance=1.0;
									else
										distance=2.0;
								}
									
								//if (distance<0.0 && sum_Ns_without_air_nodes>0.01)
								//	vel = vel_without_air_nodes / sum_Ns_without_air_nodes ;
								//else
									vel = vel_complete;
								
								//test:
								/*
								if (j==2 && add_2_water_node==true)
									distance=-1.0;
								if (j==3 && add_1_water_node==true)
									distance=-1.0;
								*/
								
								
								//pnode->GetValue(NEIGHBOUR_ELEMENTS)= Element::WeakPointer( *i.base() ) ;
								//pparticle->GetValue(NEIGHBOUR_ELEMENTS).resize(0);
								//pparticle.GetElement()=(Element::Pointer( *ielem.base() ) );
								//WeakPointerVector<Element > & pelement = (pparticle->GetValue(NEIGHBOUR_ELEMENTS));
								//pnode->GetSolutionStepValue(MEAN_SIZE)=1.0+0.01; //added a little to avoid problems when truncating with int(mean_size)
								pparticle.GetEraseFlag()=false;
								
								if (distance<=0.0)
								{
									distance=-1.0;
								}
								else //distance>1.5
								{
									distance=1.0;	
								}
								
								//pnode->GetSolutionStepValue(DISPLACEMENT_X)=0.0;
								//pnode->GetSolutionStepValue(DISPLACEMENT_Y)=0.0;
								//pnode->GetSolutionStepValue(DISPLACEMENT_Z)=0.0;
								//array_1d<double, 3 > & disp = pparticle->FastGetSolutionStepValue(DISPLACEMENT);
								//noalias(disp) = ZeroVector(3);

								//AddUniqueWeakPointer< Node<3> >(ielem->GetValue(NEIGHBOUR_NODES), pparticle);
								//if(distance<1.5)
								//{
									mparticles_vector[freeparticle]=pparticle;
									element_particle_pointers(offset+number_of_particles_in_elem) = &mparticles_vector[freeparticle];
									number_of_particles_in_elem++;
								//}
								//KRATOS_WATCH(number_of_particles_in_elem)
								if (keep_looking)
								{
									KRATOS_ERROR(std::logic_error, "FINISHED THE LIST AND COULDNT FIND A FREE CELL FOR THE NEW PARTICLE!", "");
								}
								else
								{
									reused_particles++;
									//KRATOS_WATCH(pparticle->FastGetSolutionStepValue(DISTANCE))
								}
							
							  }
						  }
					  }
				  }
				  //KRATOS_WATCH(reused_particles)
				  //KRATOS_WATCH(new_particles)
			}

			KRATOS_CATCH("")
		}
		
					

	protected:


	private:
	
	inline bool ChechIfElemIsInBoundingBox(Geometry< Node<3> >& geom, array_1d<double, TDim > bounding_box_lower_corner, array_1d<double, TDim > bounding_box_upper_corner)
	{
		array_1d<double,3> GP_position=ZeroVector(3);
		for(unsigned int j=0; j<(TDim+1); j++)
		{
			noalias(GP_position) += geom[j].Coordinates();
		}
		GP_position /= (double)(TDim+1);
		//bool is iniside=true;
		for(unsigned int j=0; j<(TDim); j++)
		{
			if(GP_position(j)<=bounding_box_lower_corner(j) || GP_position(j)>=bounding_box_upper_corner(j))
				return false;
				//is iniside=false;
		}
		//return is_inside;
		return true;
	}
	
	///this function moves a particle according to the "velocity" given
	///by "rVariable". The movement is performed in nsubsteps, during a total time
	///of Dt
	void MoveParticle(   const bool viscosity_integrate,
						 PFEM_Particle_Fluid & pparticle,
						 Element::Pointer & pelement,
						 WeakPointerVector< Element >& elements_in_trajectory,
						 unsigned int & number_of_elements_in_trajectory,						 
						 ResultIteratorType result_begin,
						 const unsigned int MaxNumberOfResults,
						 const bool & add_gravity_to_flying_particles,
						 const array_1d<double,3> mesh_displacement)
	{
		
		ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
		double delta_t = CurrentProcessInfo[DELTA_TIME];
		array_1d<double,3> & gravity= CurrentProcessInfo[GRAVITY];
		unsigned int nsubsteps;
		double substep_dt;
		
		
	    bool KEEP_INTEGRATING=false;
		bool is_found;
		//bool have_air_node;
		//bool have_water_node;
		
		array_1d<double,3> vel;
		array_1d<double,3> vel_without_other_phase_nodes=ZeroVector(3);
		array_1d<double,3> position;
		array_1d<double,3> mid_position;
		array_1d<double,TDim+1> N;
		
		//we start with the first position, then it will enter the loop.
		position = pparticle.Coordinates(); //initial coordinates
		
		//double weight;
		double g_value= 0.0;
		double g_value_without_other_phase_nodes= 0.0;
		double g_integral = 0.0;  //one for the left side, the other for the right side
		
		array_1d<double,3> viscosity_integral=ZeroVector(3);
		array_1d<double,3> viscosity_value=ZeroVector(3);
		array_1d<double,3> viscosity_value_without_other_phase_nodes=ZeroVector(3);
		/*
		array_1d<double,3> press_proj_integral=ZeroVector(3);
		array_1d<double,3> press_proj_value=ZeroVector(3);
		array_1d<double,3> press_proj_value_without_other_phase_nodes=ZeroVector(3);
		*/
		///*****
		const double particle_distance = pparticle.GetDistance();
		array_1d<double,3> particle_velocity = pparticle.GetVelocity();
		double distance=0.0;
		array_1d<double,3> last_useful_vel;
		double sum_Ns_without_other_phase_nodes;
		//double pressure=0.0;
		///*****
		bool flying_water_particle=true; //if a water particle does not find a water element in its whole path, then we add the gravity*dt
		double only_integral  = 0.0 ;
		///
		///vel=(pparticle)->FastGetSolutionStepValue(VELOCITY);
		///
		array_1d<double,3> initial_viscosity=ZeroVector(3);
		array_1d<double,3> final_viscosity=ZeroVector(3);
		
		
		is_found = FindNodeOnMesh(position, N ,pelement,result_begin,MaxNumberOfResults); //good, now we know where this point is:
		if(is_found == true)
		{
			KEEP_INTEGRATING=true;
			Geometry< Node<3> >& geom = pelement->GetGeometry();//the element we're in			
			vel=ZeroVector(3);
			vel_without_other_phase_nodes = ZeroVector(3);
			g_value=0.0;
			g_value_without_other_phase_nodes=0.0;
			viscosity_value=ZeroVector(3);
			viscosity_value_without_other_phase_nodes = ZeroVector(3);
			//press_proj_value=ZeroVector(3);
			//press_proj_value_without_other_phase_nodes = ZeroVector(3);
			
			//weight=0.0;
			sum_Ns_without_other_phase_nodes=0.0;
			//have_air_node = false;
			//have_water_node = false;
			distance=0.0;
			
			//if (particle_distance>0.0)
			if(true)
			{
				for(unsigned int j=0; j<(TDim+1); j++)
				{
					noalias(vel) += geom[j].FastGetSolutionStepValue(VELOCITY)*N[j]; 
					noalias(viscosity_value) += geom[j].GetSolutionStepValue(RHS,0)*N[j]; 
					//noalias(press_proj_value) += geom[j].GetSolutionStepValue(PRESS_PROJ_NO_RO,0)*N[j]; 
					g_value += geom[j].FastGetSolutionStepValue(G_VALUE)*N[j];
					//pressure += geom[j].FastGetSolutionStepValue(PRESSURE)*N[j];
				}
				flying_water_particle=false;
			}
			else //water particle, be careful!
			{
				for(unsigned int j=0; j<(TDim+1); j++)
				{
					if ((geom[j].FastGetSolutionStepValue(DISTANCE))<0.0) //ok. useful info!
					{
						sum_Ns_without_other_phase_nodes += N[j];
						noalias(vel_without_other_phase_nodes) += geom[j].FastGetSolutionStepValue(VELOCITY)*N[j]; 
						noalias(viscosity_value_without_other_phase_nodes) += geom[j].GetSolutionStepValue(RHS,0)*N[j]; 
						//noalias(press_proj_value_without_other_phase_nodes) += geom[j].GetSolutionStepValue(PRESS_PROJ_NO_RO,0)*N[j]; 
						g_value_without_other_phase_nodes += geom[j].FastGetSolutionStepValue(G_VALUE)*N[j];
						//have_water_node=true;
					}
					//else
						//have_air_node=true;
						
					noalias(vel) += geom[j].FastGetSolutionStepValue(VELOCITY)*N[j]; 
					noalias(viscosity_value) += geom[j].GetSolutionStepValue(RHS,0)*N[j]; 
					//noalias(press_proj_value) += geom[j].GetSolutionStepValue(PRESS_PROJ_NO_RO,0)*N[j]; 
					g_value += geom[j].FastGetSolutionStepValue(G_VALUE)*N[j];
					distance += geom[j].FastGetSolutionStepValue(DISTANCE)*N[j];
					//pressure += geom[j].FastGetSolutionStepValue(PRESSURE)*N[j];
						
						
				}
				
				//if (have_water_node)
				//if (distance<0.0)
				if (sum_Ns_without_other_phase_nodes>0.01)
				{
					vel  = vel_without_other_phase_nodes / sum_Ns_without_other_phase_nodes;
					viscosity_value = viscosity_value_without_other_phase_nodes / sum_Ns_without_other_phase_nodes;
					//press_proj_value = press_proj_value_without_other_phase_nodes / sum_Ns_without_other_phase_nodes;
					g_value = g_value_without_other_phase_nodes / sum_Ns_without_other_phase_nodes;
					flying_water_particle=false;
				}
				else
					vel = particle_velocity;
								
			}
			
			//we save pressure in case we want to convect it:
			//pparticle.GetPressure()=pressure; //no need to, done in accelerateparticleusingdeltavelocity
			
			//calculating substep to get +- courant(substep) = 1/4
			nsubsteps = 10.0 * (delta_t * pelement->GetValue(VELOCITY_OVER_ELEM_SIZE));
			if (nsubsteps<1)
				nsubsteps=1;
			substep_dt = delta_t / double(nsubsteps);
			
			g_integral+=g_value;//*double(nsubsteps);
			viscosity_integral=viscosity_value;//*double(nsubsteps);
			//press_proj_integral=press_proj_value;
			//accel_integral=ZeroVector(3);
			only_integral = 1.0;// weight;//*double(nsubsteps);

			///position += (vel+accel_value*0.5*substep_dt)*substep_dt;//weight;
			///vel += accel_value*substep_dt;
			position += vel*substep_dt;//weight;
			
			///*****
			last_useful_vel=vel;
			///*****
			//initial_accel=accel_value;
			//final_accel=accel_value;
			//initial_accel=ZeroVector(3);
			//final_accel=ZeroVector(3);

			//DONE THE FIRST LOCATION OF THE PARTICLE, NOW WE PROCEED TO STREAMLINE INTEGRATION USING THE MESH VELOCITY
			//////////////////////////////////////////////////////////////////////////////////////////////////////
			unsigned int check_from_element_number=0;
			
			
			for(unsigned int i=0; i<(nsubsteps-1); i++)// this is for the substeps n+1. in the first one we already knew the position of the particle.
			{ 
			  if (KEEP_INTEGRATING==true) 
			  {
				is_found = FindNodeOnMesh(position, N ,pelement,elements_in_trajectory,number_of_elements_in_trajectory,check_from_element_number,result_begin,MaxNumberOfResults); //good, now we know where this point is:
				if(is_found == true)
				{
					Geometry< Node<3> >& geom = pelement->GetGeometry();//the element we're in
					g_value=0.0;
					g_value_without_other_phase_nodes=0.0;
					//viscosity_value=ZeroVector(3);
					//viscosity_value_without_other_phase_nodes = ZeroVector(3);
					//press_proj_value=ZeroVector(3);
					viscosity_value_without_other_phase_nodes = ZeroVector(3);
					//press_proj_value_without_other_phase_nodes = ZeroVector(3);
					//weight=0.0;
					///*****
					distance=0.0;
					sum_Ns_without_other_phase_nodes=0.0;
					//have_air_node = false;
					//have_water_node = false;
					///*****
					
					//if (particle_distance>0.0)
					if(true)
					{
						double yp_total=0.0;
						vel_without_other_phase_nodes = ZeroVector(3);
						vel = ZeroVector(3);
						//press_proj_value=ZeroVector(3);
						viscosity_value=ZeroVector(3);
						
						for(unsigned int j=0; j<(TDim+1); j++)
						{
							if ((geom[j].FastGetSolutionStepValue(YP))>1e-12 && (geom[j].FastGetSolutionStepValue(DISTANCE)*particle_distance)>0.0)
							{
								noalias(vel_without_other_phase_nodes) += geom[j].FastGetSolutionStepValue(VELOCITY)*N[j] ;
								g_value += geom[j].FastGetSolutionStepValue(G_VALUE)*N[j];
								noalias(viscosity_value) += geom[j].FastGetSolutionStepValue(RHS)*N[j];  ///CHANGE THIS!!!!!!!!!, it should be ACCELeRATION!!  also lines above!
								//noalias(press_proj_value) += geom[j].FastGetSolutionStepValue(PRESS_PROJ_NO_RO)*N[j];  ///CHANGE THIS!!!!!!!!!, it should be ACCELeRATION!!  also lines above!
								///*****
								distance+=N[j]*geom[j].FastGetSolutionStepValue(DISTANCE);
								sum_Ns_without_other_phase_nodes += N[j];
							}
							noalias(vel) += geom[j].FastGetSolutionStepValue(VELOCITY)*N[j]; 
							///*****	
						}
						if (sum_Ns_without_other_phase_nodes>1e-12)
						//if(true)
						{
							//vel=vel_without_other_phase_nodes/sum_Ns_without_other_phase_nodes;
						}
						else
						{
							//vel = vel+gravity*substep_dt;
						}
						
					}
					else //water particle, be careful!
					{
						vel_without_other_phase_nodes = ZeroVector(3);
						
						for(unsigned int j=0; j<TDim+1; j++)
						{
							if ((geom[j].FastGetSolutionStepValue(DISTANCE))<0.0) //ok. useful info!
							{
								sum_Ns_without_other_phase_nodes += N[j];
								noalias(vel_without_other_phase_nodes) += geom[j].FastGetSolutionStepValue(VELOCITY)*N[j]; 
								noalias(viscosity_value_without_other_phase_nodes) += geom[j].GetSolutionStepValue(RHS,0)*N[j]; 
								//noalias(press_proj_value_without_other_phase_nodes) += geom[j].GetSolutionStepValue(PRESS_PROJ_NO_RO,0)*N[j]; 
								g_value_without_other_phase_nodes += geom[j].FastGetSolutionStepValue(G_VALUE)*N[j];
								//have_water_node=true;
							}
							//else
								//have_air_node=true;
							
							noalias(vel) += geom[j].FastGetSolutionStepValue(VELOCITY)*N[j]; 
							//noalias(viscosity_value) += geom[j].GetSolutionStepValue(RHS,0)*N[j]; 
							//noalias(press_proj_value) += geom[j].GetSolutionStepValue(PRESS_PROJ_NO_RO,0)*N[j]; 
							g_value += geom[j].FastGetSolutionStepValue(G_VALUE)*N[j];
							distance+=N[j]*geom[j].FastGetSolutionStepValue(DISTANCE);
						}
					
						//if (have_water_node)
						//if (distance<0.0)
						if (sum_Ns_without_other_phase_nodes>0.01)
						{
							vel  = vel_without_other_phase_nodes / sum_Ns_without_other_phase_nodes;
							viscosity_value = viscosity_value_without_other_phase_nodes / sum_Ns_without_other_phase_nodes;
							//press_proj_value = press_proj_value_without_other_phase_nodes / sum_Ns_without_other_phase_nodes;
							g_value = g_value_without_other_phase_nodes / sum_Ns_without_other_phase_nodes;
							flying_water_particle=false;
						}
						else
						{
							particle_velocity += substep_dt * gravity;
							
							//warning! no parabolic trajectory!
							//vel  = particle_velocity;
						}
					}
					
					//KRATOS_WATCH(vel)
					//KRATOS_WATCH(node_weight)
					//KRATOS_WATCH(weight)
					
					g_integral +=  g_value;
					viscosity_integral += (viscosity_value); //-initial_accel);
					//press_proj_integral += (press_proj_value); //-initial_accel);
					//final_accel=accel_value;
					only_integral += 1.0; //values saved for the current time step
					
					///position += (vel+accel_value*0.5*substep_dt)*substep_dt;//weight;
					///vel += accel_value*substep_dt;
					
					///*****
					//if (particle_distance<0.0 && (particle_distance*distance)<0.0) //meaning it is a water particle and we are entering air domain. WE CANNOT FOLLOW THIS STREAMLINE!
					//	vel=last_useful_vel;
					//else //good. we can do the usual thing. and we save this velocity in case in the next step we go out of the water domain
					//	last_useful_vel=vel;
					///*****
					position+=vel*substep_dt;//weight;
					


				  }
				  else 
				  {
					  KEEP_INTEGRATING=false;		
					  break;	
				  }
				}
				else
					break;
				

			}
		
		}
		//noalias(pparticle->Coordinates()) = position;
		
		
		//if there's a mesh velocity, we add it at the end in a single step:
		position-=mesh_displacement;
		
		
		if (KEEP_INTEGRATING==false) (pparticle.GetEraseFlag()=true);
		else is_found = FindNodeOnMesh(position, N ,pelement,result_begin,MaxNumberOfResults); //we must save the pointer of the last element that we're in (inside the pointervector pelement)
		
		//if (is_found==false) ( (pparticle)->FastGetSolutionStepValue(MEAN_SIZE) =-1.0);
		if (is_found==false) ( pparticle.GetEraseFlag()=true);

		 pparticle.Coordinates() = position;
		 //pparticle.GetElement() = pelement;
		
		//WARNING! TEMPEATURE PROBLEM IS OFF!
		///(pparticle)->FastGetSolutionStepValue(TEMPERATURE)= (pparticle)->FastGetSolutionStepValue(TEMPERATURE)+ g_integral/only_integral;
		
		
		if (add_gravity_to_flying_particles)
		{ 
			if (flying_water_particle==true)
				pparticle.GetVelocity() = pparticle.GetVelocity() + gravity*delta_t;
		}
		
		///WATCH THIS LINE!!!!! I'M DUMPING THE ORIGINAL DATA (and replacing it only based on fhe fixed mesh) INSTEAD OF ADDING A CONTRIBUTION!!! NOT THE WAY PFEM IS INTENDEED TO WORK!
		///FIX WAY TO SET BOUNDARY CONDITIONS!
		if (viscosity_integrate)
			pparticle.GetVelocity() = pparticle.GetVelocity() + (viscosity_integral  )*substep_dt; //- (final_accel-initial_accel)*1.0*double(nsubsteps)
			
			

		//if (pressure_gradient_integrate)
		//	pparticle.GetVelocity() = pparticle.GetVelocity() - (press_proj_integral  )*substep_dt; //- (final_accel-initial_accel)*1.0*double(nsubsteps)
		///(pparticle)->FastGetSolutionStepValue(VELOCITY)=vel;
		//(pparticle)->GetSolutionStepValue(VELOCITY,0)=accel_integral/only_integral; 
		///check 20 lines above, i'm using velocity instead of accelartion!!->	noalias(accel_value)
				

	}
	

	void AccelerateParticleUsingDeltaVelocity(
						 PFEM_Particle_Fluid & pparticle,
						 Element::Pointer & pelement,
						 Geometry< Node<3> >& geom)
	{
		//bool is_found;

		//array_1d<double,3> position;
		array_1d<double,TDim+1> N;
		
		ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
		double delta_t = CurrentProcessInfo[DELTA_TIME];
		array_1d<double,3> gravity = CurrentProcessInfo[GRAVITY];
		
		
		//we start with the first position, then it will enter the loop.
		array_1d<double,3> coords = pparticle.Coordinates();// + (pparticle)->FastGetSolutionStepValue(DISPLACEMENT); //initial coordinates
		double & particle_distance = pparticle.GetDistance();
	//	double & particle_temperature = pparticle.GetTemperature();
		//double temperature=0.0;
		//double temperature_without_air=0.0;
		//double temperature_without_water=0.0;
		//double delta_temperature= 0.0;
		//double delta_temperature_without_water= 0.0;
		//double delta_temperature_without_air= 0.0;
		double distance=0.0;
		//double pressure=0.0;
		array_1d<double,3> delta_velocity = ZeroVector(3);
		
		array_1d<double,3> delta_velocity_without_air = ZeroVector(3);		
		array_1d<double,3> delta_velocity_without_water = ZeroVector(3);		
			
		bool is_found = CalculatePosition(geom,coords[0],coords[1],coords[2],N);
		if(is_found == false)
		{
			KRATOS_WATCH(N)
			for (int j=0 ; j!=(TDim+1); j++)
								if (N[j]<0.0 )
									N[j]=1e-10;
			
			//KRATOS_ERROR(std::logic_error, "PARTICLE IN WRONG ELEMENT!", "");
		}
		
			//double distance=0.0;
			//g_value=0.0;
				
			//if (particle_distance>0.0) //no problem. air
			if(true)
			{
				double sum_Ns_without_water_nodes = 0.0;
				//bool have_air_node=false;
				for(unsigned int j=0; j<(TDim+1); j++)
				{
					if ((geom[j].FastGetSolutionStepValue(DISTANCE))>0.0)
					{
						noalias(delta_velocity_without_water) += geom[j].FastGetSolutionStepValue(DELTA_VELOCITY)*N[j];
						//delta_temperature_without_water += geom[j].FastGetSolutionStepValue(DELTA_TEMPERATURE)*N[j];
						sum_Ns_without_water_nodes += N[j];
						//temperature_without_water += geom[j].FastGetSolutionStepValue(TEMPERATURE)*N[j];
						//have_air_node=true;
						//g_value += geom[j].FastGetSolutionStepValue(G_VALUE)*N[j];
						
					}
					noalias(delta_velocity) += geom[j].FastGetSolutionStepValue(DELTA_VELOCITY)*N[j];
					distance += geom[j].FastGetSolutionStepValue(DISTANCE)*N[j];
					//delta_temperature += geom[j].FastGetSolutionStepValue(DELTA_TEMPERATURE)*N[j];
					//temperature += geom[j].FastGetSolutionStepValue(TEMPERATURE)*N[j];
					//pressure += geom[j].FastGetSolutionStepValue(PRESSURE)*N[j];
					
				}
				
				//if (distance>0.0)
				if (sum_Ns_without_water_nodes>0.01)
				//if(false)
				{
					//delta_velocity = delta_velocity_without_water/sum_Ns_without_water_nodes ;
					//delta_temperature = delta_temperature_without_water/sum_Ns_without_water_nodes ;
					//temperature = temperature_without_water/sum_Ns_without_water_nodes ;
				}
					
			}
			else //water particle
			{
				double sum_Ns_without_air_nodes = 0.0;
				//bool have_water_node=false;
				for(unsigned int j=0; j<(TDim+1); j++)
				{
					if ((geom[j].FastGetSolutionStepValue(DISTANCE))<0.0)
					{
						noalias(delta_velocity_without_air) += geom[j].FastGetSolutionStepValue(DELTA_VELOCITY)*N[j];
						sum_Ns_without_air_nodes += N[j];
						//delta_temperature_without_air += geom[j].FastGetSolutionStepValue(DELTA_TEMPERATURE)*N[j];
						//temperature_without_air += geom[j].FastGetSolutionStepValue(TEMPERATURE)*N[j];
						//have_water_node=true;
						//g_value += geom[j].FastGetSolutionStepValue(G_VALUE)*N[j];
						
					}
					noalias(delta_velocity) += geom[j].FastGetSolutionStepValue(DELTA_VELOCITY)*N[j];
					distance += geom[j].FastGetSolutionStepValue(DISTANCE)*N[j];
					//delta_temperature += geom[j].FastGetSolutionStepValue(DELTA_TEMPERATURE)*N[j];
					//temperature += geom[j].FastGetSolutionStepValue(TEMPERATURE)*N[j];
					//pressure += geom[j].FastGetSolutionStepValue(PRESSURE)*N[j];
				}
				
				
				//if (false)
				if (sum_Ns_without_air_nodes>0.01)
				{
					delta_velocity = delta_velocity_without_air/sum_Ns_without_air_nodes ;
					//delta_temperature = delta_temperature_without_air/sum_Ns_without_air_nodes ;
					//temperature = temperature_without_air/sum_Ns_without_air_nodes ;
				}
				else
				{
					if (mDENSITY_WATER>(10.0*mDENSITY_AIR))
					{
						delta_velocity=gravity*(1.0-mDENSITY_AIR/mDENSITY_WATER)*delta_t;
					}
				}

			}
			pparticle.GetVelocity() = pparticle.GetVelocity() + delta_velocity;//*substep_dt*double(nsubsteps)/particle_density;///*only_integral;//only_integral; 
			//pparticle.GetTemperature() = pparticle.GetTemperature() + delta_temperature;
			//pparticle.GetTemperature() = temperature;
			//pparticle.GetPressure() = pressure;
		


	}
		
	void MoveParticle_inverse_way(   
						 const bool viscosity_integrate,
						 PFEM_Particle_Fluid & pparticle,
						 Element::Pointer & pelement, //NOT A REFERENCE!! WE SHALL NOT OVERWRITE THE ELEMENT IT BELONGS TO!						 
						 ResultIteratorType result_begin,
						 const unsigned int MaxNumberOfResults)
	{
		
		ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
		double delta_t = CurrentProcessInfo[DELTA_TIME];
		unsigned int nsubsteps;
		double substep_dt;
		
		
	    bool KEEP_INTEGRATING=false;
		bool is_found;
		
		array_1d<double,3> vel;
		array_1d<double,3> position;
		array_1d<double,3> mid_position;
		array_1d<double,TDim+1> N;
		
		//array_1d<double,3> & press_proj = (pparticle)->FastGetSolutionStepValue(PRESS_PROJ);
		//double & pressure = (pparticle)->FastGetSolutionStepValue(PRESSURE);
		
		//we start with the first position, then it will enter the loop.
		position = pparticle.Coordinates(); // + (pparticle)->FastGetSolutionStepValue(DISPLACEMENT); //initial coordinates
		
		//double weight;
		double g_value= 0.0;
		double g_integral = 0.0;  //one for the left side, the other for the right side
		//double & temp = pparticle.GetTemperature();
		double & distance = pparticle.GetDistance();
		//double & oxygen = pparticle.GetOxygen();
		array_1d<double,3> viscosity_integral=ZeroVector(3);
		array_1d<double,3> viscosity_value=ZeroVector(3);
		
		//array_1d<double,3> press_proj_integral=ZeroVector(3);
		//array_1d<double,3> press_proj_value=ZeroVector(3);

		//array_1d<double,3> initial_accel=ZeroVector(3);
		//array_1d<double,3> final_accel=ZeroVector(3);
		
		double only_integral  = 0.0 ;
		//double pressure=0.0;
		//double gradient_discontinuity=0.0;
		double pressure=0.0;
		double solid_pressure=0.0;
		

		
		is_found = FindNodeOnMesh(position, N ,pelement,result_begin,MaxNumberOfResults); //good, now we know where this point is:
		if(is_found == true)
		{
			KEEP_INTEGRATING=true;
			Geometry< Node<3> >& geom = pelement->GetGeometry();//the element we're in			
			vel=ZeroVector(3);
			g_value=0.0;
			viscosity_value=ZeroVector(3);
			//oxygen=0.0;
			//weight=0.0;
			distance=0.0;
			//pressure=0.0;
			//press_proj_value=ZeroVector(3);
			
			for(unsigned int j=0; j<(TDim+1); j++)
			{
				distance +=  geom[j].FastGetSolutionStepValue(DISTANCE)*N(j);
				noalias(vel) += geom[j].FastGetSolutionStepValue(VELOCITY)*N[j]; 
				noalias(viscosity_value) += geom[j].GetSolutionStepValue(RHS,0)*N[j];  ///CHANGE THIS!!!!!!!!!, it should be ACCELeRATION!!  also lines below!
				//noalias(press_proj_value) += geom[j].GetSolutionStepValue(PRESS_PROJ_NO_RO,0)*N[j];  ///CHANGE THIS!!!!!!!!!, it should be ACCELeRATION!!  also lines below!
				g_value += geom[j].FastGetSolutionStepValue(G_VALUE)*N[j];
				solid_pressure += geom[j].FastGetSolutionStepValue(SOLID_PRESSURE)*N[j];
				pressure += geom[j].FastGetSolutionStepValue(PRESSURE)*N[j];
				//oxygen += geom[j].FastGetSolutionStepValue(OXYGEN_FRACTION)*N[j];
			}
			//calculating substep to get +- courant(substep) = 1/4
			nsubsteps = 10.0 * (delta_t * pelement->GetValue(VELOCITY_OVER_ELEM_SIZE));
			if (nsubsteps<1)
				nsubsteps=1;
			substep_dt = delta_t / double(nsubsteps);
			
			g_integral+=g_value;//*double(nsubsteps);
			viscosity_integral=viscosity_value;//*double(nsubsteps);
			//press_proj_integral=press_proj_value;//*double(nsubsteps);
			only_integral = 1.0;// weight;//*double(nsubsteps);
			position -= vel*substep_dt;//weight;
			//gradient_discontinuity=pelement->GetValue(GRADIENT_DISCONTINUITY);


			
			
			//initial_accel=accel_value;
			//final_accel=accel_value;
			
			for(unsigned int i=0; i<(nsubsteps-1); i++)// this is for the substeps n+1. in the first one we already knew the position of the particle.
			{ if (KEEP_INTEGRATING==true) {
				is_found = FindNodeOnMesh(position, N ,pelement,result_begin,MaxNumberOfResults); //good, now we know where this point is:
				if(is_found == true)
				{
					Geometry< Node<3> >& geom = pelement->GetGeometry();//the element we're in
			
					vel=ZeroVector(3);
					g_value=0.0;
					viscosity_value=ZeroVector(3);	
					//press_proj_value=ZeroVector(3);	
					//weight=0.0;
					//temp=0.0;
					distance=0.0;
					//oxygen=0.0;
					solid_pressure=0.0;
					pressure=0.0;
					
					for(unsigned int j=0; j<(TDim+1); j++)
					{
						//pressure += geom[j].FastGetSolutionStepValue(PRESSURE)*N[j];
						noalias(vel) += geom[j].FastGetSolutionStepValue(VELOCITY)*N[j] ;
						//temp += geom[j].FastGetSolutionStepValue(TEMPERATURE)*N[j];
						g_value += geom[j].FastGetSolutionStepValue(G_VALUE)*N[j];
						distance +=  geom[j].FastGetSolutionStepValue(DISTANCE)*N(j);
						noalias(viscosity_value) += geom[j].FastGetSolutionStepValue(RHS)*N[j];  ///CHANGE THIS!!!!!!!!!, it should be ACCELeRATION!!  also lines above!
						//oxygen += geom[j].FastGetSolutionStepValue(OXYGEN_FRACTION)*N[j];
						//noalias(press_proj_value) += geom[j].FastGetSolutionStepValue(PRESS_PROJ_NO_RO)*N[j];  ///CHANGE THIS!!!!!!!!!, it should be ACCELeRATION!!  also lines above!
						solid_pressure += geom[j].FastGetSolutionStepValue(SOLID_PRESSURE)*N(j);
						pressure += geom[j].FastGetSolutionStepValue(PRESSURE)*N(j);
					}

		
					g_integral +=  g_value ;
					viscosity_integral += viscosity_value;// - initial_accel;
					//press_proj_integral += press_proj_value;// - initial_accel;
					//final_accel=accel_value;
					only_integral += 1.0;//weight ; //values saved for the current time step
					position-=vel*substep_dt;//weight;
					//gradient_discontinuity=pelement->GetValue(GRADIENT_DISCONTINUITY);


				  }
				  else KEEP_INTEGRATING=false;			
				}
				

			}
			
			///COMMENT TO GET A A CONTINOUS DISTANCE FUNCTION FIELD!!!!!
			if(distance>0.0)
				distance=1.0;
			else
				distance=-1.0;
			

			
			pparticle.GetVelocity()=vel;
			
			if (viscosity_integrate)
				pparticle.GetVelocity() += (viscosity_integral  )*substep_dt; // - (final_accel-initial_accel)*1.0*double(nsubsteps)
			//if (pressure_gradient_integrate)
			//	pparticle.GetVelocity() -= (press_proj_integral  )*substep_dt;
				
			//pparticle.GetPressure() = pressure; 
			//pparticle.GetGradientDiscontinuity()=gradient_discontinuity;
		}
		else {KRATOS_WATCH(position); }		

	}
	
	
	
	
	
	
	
	void OverwriteParticleDataUsingTopographicDomain(
						 PFEM_Particle_Fluid & pparticle,
						 Element::Pointer & pelement,
						 array_1d<double,3> domains_offset,
						 ResultIteratorType result_begin,
						 const unsigned int MaxNumberOfResults)
	{
		array_1d<double,TDim+1> N;
		
		ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
		
		//we start with the first position, then it will enter the loop.
		array_1d<double,3> coords = pparticle.Coordinates()+domains_offset; 
		double & particle_distance = pparticle.GetDistance();
		bool is_found = FindNodeOnTopographicMesh(coords, N ,pelement,result_begin,MaxNumberOfResults); //good, now we know where this point is:
		
		if (is_found) //it is part of the solid topographic domain
		{
			particle_distance= -1.0;
			/*
			double distance=0.0;
			for(unsigned int j=0; j<(TDim+1); j++)
			{
				distance += geom[j].FastGetSolutionStepValue(DISTANCE)*N[j];
			}
			*/ 
		}
		else //it is outside the topographic domain, therefore it is air or whatever it means
		{
			particle_distance= 1.0;
		}
			

		pparticle.GetVelocity() = ZeroVector(3);
	}
	
	
	
	
	///this function should find the element into which a given node is located
	///and return a pointer to the element and the vector containing the 
	///shape functions that define the postion within the element
	///if "false" is devolved the element is not found
	bool FindNodeOnMesh( //int last_element,
						 array_1d<double,3>& position,
						 array_1d<double,TDim+1>& N,
						 //Element::Pointer pelement,
						 Element::Pointer & pelement,
						 ResultIteratorType result_begin,
						 const unsigned int MaxNumberOfResults)
	{
		typedef std::size_t SizeType; 
		
		const array_1d<double,3>& coords = position;
		 array_1d<double,TDim+1> aux_N;
	    //before using the bin to search for possible elements we check first the last element in which the particle was.
	    
		//ModelPart::ElementsContainerType::iterator i = mr_model_part.ElementsBegin()+last_element; 
		Geometry<Node<3> >& geom_default = pelement->GetGeometry(); //(*(i))->GetGeometry();
		bool is_found_1 = CalculatePosition(geom_default,coords[0],coords[1],coords[2],N);
		if(is_found_1 == true)
		{
			//pelement = (*(i));
			return true;
		}
		
		//KRATOS_WATCH("will look in another element")
		//KRATOS_WATCH(TDim)
		
		//to begin with we check the neighbour elements:
		WeakPointerVector< Element >& neighb_elems = pelement->GetValue(NEIGHBOUR_ELEMENTS);
		//the first we check is the one that has negative shape function, because it means it went outside in this direction:
		/*
		unsigned int checked_element=0;
		for (unsigned int i=0;i!=(TDim+1);i++)
		{
			if (N[i]<0.0)
			{
				checked_element=i;
				Geometry<Node<3> >& geom = neighb_elems[i].GetGeometry();
				bool is_found_2 = CalculatePosition(geom,coords[0],coords[1],coords[2],aux_N);
				if (is_found_2)
				{
					pelement=Element::Pointer(((neighb_elems(i))));
					N=aux_N;
					return true;
				}
				break;
			}
		}
		*/
		
		for (unsigned int i=0;i!=(neighb_elems.size());i++)
		{

				Geometry<Node<3> >& geom = neighb_elems[i].GetGeometry();
				bool is_found_2 = CalculatePosition(geom,coords[0],coords[1],coords[2],N);
				if (is_found_2)
				{
					pelement=Element::Pointer(((neighb_elems(i))));
					return true;
				}
		}
			
	    
		//ask to the container for the list of candidate elements
		SizeType results_found = mpBinsObjectDynamic->SearchObjectsInCell(coords, result_begin, MaxNumberOfResults );
		//KRATOS_WATCH(results_found)
				
		if(results_found>0){
		//loop over the candidate elements and check if the particle falls within
		for(SizeType i = 0; i< results_found; i++)
		{
			//std::cout<< "KIIIIIIIIIIIIII" << std::endl;
			//KRATOS_WATCH((*(result_begin+i))->Id());
			Geometry<Node<3> >& geom = (*(result_begin+i))->GetGeometry();
			
			
			//find local position
			bool is_found = CalculatePosition(geom,coords[0],coords[1],coords[2],N);
			
			//KRATOS_WATCH("ln243");
			//KRATOS_WATCH(N);
			
			if(is_found == true)
			{
				//pelement.clear();
				//pelement.push_back( Element::WeakPointer((*(result_begin+i).base())));
				pelement=Element::Pointer((*(result_begin+i).base()));
				return true;
			}
		}
	}
		
		//not found case
		return false;
	}
	
	
	// VERSION INCLUDING PREDEFINED ELEMENTS FOLLOWING A TRAJECTORY
		bool FindNodeOnMesh( //int last_element,
						 array_1d<double,3>& position,
						 array_1d<double,TDim+1>& N,
						 //Element::Pointer pelement,
						 Element::Pointer & pelement,
						 WeakPointerVector< Element >& elements_in_trajectory,
						 unsigned int & number_of_elements_in_trajectory,
						 unsigned int & check_from_element_number,
						 ResultIteratorType result_begin,
						 const unsigned int MaxNumberOfResults)
	{
		typedef std::size_t SizeType; 
		
		const array_1d<double,3>& coords = position;
		 array_1d<double,TDim+1> aux_N;
	    //before using the bin to search for possible elements we check first the last element in which the particle was.
	    
		//ModelPart::ElementsContainerType::iterator i = mr_model_part.ElementsBegin()+last_element; 
		Geometry<Node<3> >& geom_default = pelement->GetGeometry(); //(*(i))->GetGeometry();
		bool is_found_1 = CalculatePosition(geom_default,coords[0],coords[1],coords[2],N);
		if(is_found_1 == true)
		{
			//pelement = (*(i));
			return true;
		}
		
		//KRATOS_WATCH("will look in another element")
		//KRATOS_WATCH(TDim)
		
		//if it was not found in the first element, we can proceed to check in the following elements.
		for (unsigned int i=(check_from_element_number);i!=number_of_elements_in_trajectory;i++)
		{
			Geometry<Node<3> >& geom = elements_in_trajectory[i].GetGeometry();
			bool is_found_2 = CalculatePosition(geom,coords[0],coords[1],coords[2],aux_N);
			if (is_found_2)
			{
				pelement=Element::Pointer(((elements_in_trajectory(i))));
				N=aux_N;
				check_from_element_number = i+1 ; //now i element matches pelement, so to avoid cheching twice the same element we send the counter to the following element.
				return true;
			}
			
		}
		
		//now we check the neighbour elements:
		WeakPointerVector< Element >& neighb_elems = pelement->GetValue(NEIGHBOUR_ELEMENTS);
		//the first we check is the one that has negative shape function, because it means it went outside in this direction:
		/*
		unsigned int checked_element=0;
		for (unsigned int i=0;i!=(TDim+1);i++)
		{
			if (N[i]<0.0)
			{
				checked_element=i;
				Geometry<Node<3> >& geom = neighb_elems[i].GetGeometry();
				bool is_found_2 = CalculatePosition(geom,coords[0],coords[1],coords[2],aux_N);
				if (is_found_2)
				{
					pelement=Element::Pointer(((neighb_elems(i))));
					N=aux_N;
					return true;
				}
				break;
			}
		}
		*/
		for (unsigned int i=0;i!=(neighb_elems.size());i++)
		{

				Geometry<Node<3> >& geom = neighb_elems[i].GetGeometry();
				bool is_found_2 = CalculatePosition(geom,coords[0],coords[1],coords[2],N);
				if (is_found_2)
				{
					pelement=Element::Pointer(((neighb_elems(i))));
					if (number_of_elements_in_trajectory<20)
					{
						elements_in_trajectory(number_of_elements_in_trajectory)=pelement;
						number_of_elements_in_trajectory++;
						check_from_element_number = number_of_elements_in_trajectory;  //we do it after doing the ++ to the counter, so we woudlnt enter the loop that searches in the elements_in_trajectory list. we are the particle that is adding elements to the list
					}
					return true;
				}
		}
			
	    
		//ask to the container for the list of candidate elements
		SizeType results_found = mpBinsObjectDynamic->SearchObjectsInCell(coords, result_begin, MaxNumberOfResults );
		//KRATOS_WATCH(results_found)
				
		if(results_found>0){
		//loop over the candidate elements and check if the particle falls within
		for(SizeType i = 0; i< results_found; i++)
		{
			//std::cout<< "KIIIIIIIIIIIIII" << std::endl;
			//KRATOS_WATCH((*(result_begin+i))->Id());
			Geometry<Node<3> >& geom = (*(result_begin+i))->GetGeometry();
			
			
			//find local position
			bool is_found = CalculatePosition(geom,coords[0],coords[1],coords[2],N);
			
			//KRATOS_WATCH("ln243");
			//KRATOS_WATCH(N);
			
			if(is_found == true)
			{
				//pelement.clear();
				//pelement.push_back( Element::WeakPointer((*(result_begin+i).base())));
				pelement=Element::Pointer((*(result_begin+i).base()));
				if (number_of_elements_in_trajectory<20)
				{
				elements_in_trajectory(number_of_elements_in_trajectory)=pelement;
				number_of_elements_in_trajectory++;
				check_from_element_number = number_of_elements_in_trajectory;  //we do it after doing the ++ to the counter, so we woudlnt enter the loop that searches in the elements_in_trajectory list. we are the particle that is adding elements to the list
				}
				return true;
			}
		}
	}
		
		//not found case
		return false;
	}

	
///this function should find the element into which a given node is located
	///and return a pointer to the element and the vector containing the 
	///shape functions that define the postion within the element
	///if "false" is devolved the element is not found
	bool FindNodeOnTopographicMesh( //int last_element,
						 array_1d<double,3>& position,
						 array_1d<double,TDim+1>& N,
						 //Element::Pointer pelement,
						 Element::Pointer & pelement,
						 ResultIteratorType result_begin,
						 const unsigned int MaxNumberOfResults)
	{
		typedef std::size_t SizeType; 
		
		const array_1d<double,3>& coords = position;
		 array_1d<double,TDim+1> aux_N;
	    //before using the bin to search for possible elements we check first the last element in which the particle was.
	    
		//ModelPart::ElementsContainerType::iterator i = mr_model_part.ElementsBegin()+last_element; 
		Geometry<Node<3> >& geom_default = pelement->GetGeometry(); //(*(i))->GetGeometry();
		bool is_found_1 = CalculatePosition(geom_default,coords[0],coords[1],coords[2],N);
		if(is_found_1 == true)
		{
			//pelement = (*(i));
			return true;
		}
		
		//KRATOS_WATCH("will look in another element")
		//KRATOS_WATCH(TDim)
		
		//to begin with we check the neighbour elements:
		WeakPointerVector< Element >& neighb_elems = pelement->GetValue(NEIGHBOUR_ELEMENTS);
		//the first we check is the one that has negative shape function, because it means it went outside in this direction:
		/*
		unsigned int checked_element=0;
		for (unsigned int i=0;i!=(TDim+1);i++)
		{
			if (N[i]<0.0)
			{
				checked_element=i;
				Geometry<Node<3> >& geom = neighb_elems[i].GetGeometry();
				bool is_found_2 = CalculatePosition(geom,coords[0],coords[1],coords[2],aux_N);
				if (is_found_2)
				{
					pelement=Element::Pointer(((neighb_elems(i))));
					N=aux_N;
					return true;
				}
				break;
			}
		}
		*/
		
		for (unsigned int i=0;i!=(neighb_elems.size());i++)
		{

				Geometry<Node<3> >& geom = neighb_elems[i].GetGeometry();
				bool is_found_2 = CalculatePosition(geom,coords[0],coords[1],coords[2],N);
				if (is_found_2)
				{
					pelement=Element::Pointer(((neighb_elems(i))));
					return true;
				}
		}
			
	    
		//ask to the container for the list of candidate elements
		SizeType results_found = mpTopographicBinsObjectDynamic->SearchObjectsInCell(coords, result_begin, MaxNumberOfResults );
		//KRATOS_WATCH(results_found)
				
		if(results_found>0){
		//loop over the candidate elements and check if the particle falls within
		for(SizeType i = 0; i< results_found; i++)
		{
			//std::cout<< "KIIIIIIIIIIIIII" << std::endl;
			//KRATOS_WATCH((*(result_begin+i))->Id());
			Geometry<Node<3> >& geom = (*(result_begin+i))->GetGeometry();
			
			
			//find local position
			bool is_found = CalculatePosition(geom,coords[0],coords[1],coords[2],N);
			
			//KRATOS_WATCH("ln243");
			//KRATOS_WATCH(N);
			
			if(is_found == true)
			{
				//pelement.clear();
				//pelement.push_back( Element::WeakPointer((*(result_begin+i).base())));
				pelement=Element::Pointer((*(result_begin+i).base()));
				return true;
			}
		}
	}
		
		//not found case
		return false;
	}
		
	//***************************************
        //***************************************

        inline bool CalculatePosition(Geometry<Node < 3 > >&geom,
                const double xc, const double yc, const double zc,
                array_1d<double, 3 > & N
                )
        {
            double x0 = geom[0].X();
            double y0 = geom[0].Y();
            double x1 = geom[1].X();
            double y1 = geom[1].Y();
            double x2 = geom[2].X();
            double y2 = geom[2].Y();

            double area = CalculateVol(x0, y0, x1, y1, x2, y2);
            double inv_area = 0.0;
            if (area == 0.0)
            {
                KRATOS_ERROR(std::logic_error, "element with zero area found", "");
            } else
            {
                inv_area = 1.0 / area;
            }


            N[0] = CalculateVol(x1, y1, x2, y2, xc, yc) * inv_area;
            N[1] = CalculateVol(x2, y2, x0, y0, xc, yc) * inv_area;
            N[2] = CalculateVol(x0, y0, x1, y1, xc, yc) * inv_area;
			//KRATOS_WATCH(N);

            if (N[0] >= 0.0 && N[1] >= 0.0 && N[2] >= 0.0 && N[0] <= 1.0 && N[1] <= 1.0 && N[2] <= 1.0) //if the xc yc is inside the triangle return true
                return true;

            return false;
        }

	        //***************************************
        //***************************************

        inline bool CalculatePosition(Geometry<Node < 3 > >&geom,
                const double xc, const double yc, const double zc,
                array_1d<double, 4 > & N
                )
        {

            double x0 = geom[0].X();
            double y0 = geom[0].Y();
            double z0 = geom[0].Z();
            double x1 = geom[1].X();
            double y1 = geom[1].Y();
            double z1 = geom[1].Z();
            double x2 = geom[2].X();
            double y2 = geom[2].Y();
            double z2 = geom[2].Z();
            double x3 = geom[3].X();
            double y3 = geom[3].Y();
            double z3 = geom[3].Z();

            double vol = CalculateVol(x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3);

            double inv_vol = 0.0;
            if (vol < 0.0000000000001)
            {
                KRATOS_ERROR(std::logic_error, "element with zero vol found", "");
            } else
            {
                inv_vol = 1.0 / vol;
            }

            N[0] = CalculateVol(x1, y1, z1, x3, y3, z3, x2, y2, z2, xc, yc, zc) * inv_vol;
            N[1] = CalculateVol(x0, y0, z0, x1, y1, z1, x2, y2, z2, xc, yc, zc) * inv_vol;
            N[2] = CalculateVol(x3, y3, z3, x1, y1, z1, x0, y0, z0, xc, yc, zc) * inv_vol;
            N[3] = CalculateVol(x3, y3, z3, x0, y0, z0, x2, y2, z2, xc, yc, zc) * inv_vol;


            if (N[0] >= 0.0 && N[1] >= 0.0 && N[2] >= 0.0 && N[3] >= 0.0 &&
                    N[0] <= 1.0 && N[1] <= 1.0 && N[2] <= 1.0 && N[3] <= 1.0)
                //if the xc yc zc is inside the tetrahedron return true
                return true;

            return false;
        }

        inline double CalculateVol(const double x0, const double y0,
                const double x1, const double y1,
                const double x2, const double y2
                )
        {
            return 0.5 * ((x1 - x0)*(y2 - y0)- (y1 - y0)*(x2 - x0));
        }
        //***************************************
        //***************************************

        inline double CalculateVol(const double x0, const double y0, const double z0,
                const double x1, const double y1, const double z1,
                const double x2, const double y2, const double z2,
                const double x3, const double y3, const double z3
                )
        {
            double x10 = x1 - x0;
            double y10 = y1 - y0;
            double z10 = z1 - z0;

            double x20 = x2 - x0;
            double y20 = y2 - y0;
            double z20 = z2 - z0;

            double x30 = x3 - x0;
            double y30 = y3 - y0;
            double z30 = z3 - z0;

            double detJ = x10 * y20 * z30 - x10 * y30 * z20 + y10 * z20 * x30 - y10 * x20 * z30 + z10 * x20 * y30 - z10 * y20 * x30;
            return detJ * 0.1666666666666666666667;
        }
        
        
        
		void ComputeGaussPointPositions_4(Geometry< Node < 3 > >& geom, boost::numeric::ublas::bounded_matrix<double, 7, 3 > & pos,boost::numeric::ublas::bounded_matrix<double, 7, 3 > & N)
        {
            double one_third = 1.0 / 3.0;
            double one_sixt = 0.15; //1.0 / 6.0;
            double two_third = 0.7; //2.0 * one_third;

            N(0, 0) = one_sixt;
            N(0, 1) = one_sixt;
            N(0, 2) = two_third;
            N(1, 0) = two_third;
            N(1, 1) = one_sixt;
            N(1, 2) = one_sixt;
            N(2, 0) = one_sixt;
            N(2, 1) = two_third;
            N(2, 2) = one_sixt;
            N(3, 0) = one_third;
			N(3, 1) = one_third;
			N(3, 2) = one_third;
			
            //first
            pos(0, 0) = one_sixt * geom[0].X() + one_sixt * geom[1].X() + two_third * geom[2].X();
            pos(0, 1) = one_sixt * geom[0].Y() + one_sixt * geom[1].Y() + two_third * geom[2].Y();
            pos(0, 2) = one_sixt * geom[0].Z() + one_sixt * geom[1].Z() + two_third * geom[2].Z();
	    
            //second
            pos(1, 0) = two_third * geom[0].X() + one_sixt * geom[1].X() + one_sixt * geom[2].X();
            pos(1, 1) = two_third * geom[0].Y() + one_sixt * geom[1].Y() + one_sixt * geom[2].Y();
            pos(1, 2) = two_third * geom[0].Z() + one_sixt * geom[1].Z() + one_sixt * geom[2].Z();

            //third
            pos(2, 0) = one_sixt * geom[0].X() + two_third * geom[1].X() + one_sixt * geom[2].X();
            pos(2, 1) = one_sixt * geom[0].Y() + two_third * geom[1].Y() + one_sixt * geom[2].Y();
            pos(2, 2) = one_sixt * geom[0].Z() + two_third * geom[1].Z() + one_sixt * geom[2].Z();

            //fourth
			pos(3, 0) = one_third * geom[0].X() + one_third * geom[1].X() + one_third * geom[2].X();
            pos(3, 1) = one_third * geom[0].Y() + one_third * geom[1].Y() + one_third * geom[2].Y();
            pos(3, 2) = one_third * geom[0].Z() + one_third * geom[1].Z() + one_third * geom[2].Z();

        }
        
        
        void ComputeGaussPointPositionsForPostReseed(Geometry< Node < 3 > >& geom, boost::numeric::ublas::bounded_matrix<double, 7, 3 > & pos,boost::numeric::ublas::bounded_matrix<double, 7, 3 > & N) //2d
        {
            double one_third = 1.0 / 3.0;
            double one_eight = 0.12; //1.0 / 6.0;
            double three_quarters = 0.76; //2.0 * one_third;

            N(0, 0) = one_eight;
            N(0, 1) = one_eight;
            N(0, 2) = three_quarters;
            
            N(1, 0) = three_quarters;
            N(1, 1) = one_eight;
            N(1, 2) = one_eight;
            
            N(2, 0) = one_eight;
            N(2, 1) = three_quarters;
            N(2, 2) = one_eight;
            
            N(3, 0) = one_third;
			N(3, 1) = one_third;
			N(3, 2) = one_third;
			
			N(4, 0) = one_eight;
            N(4, 1) = 0.44;
            N(4, 2) = 0.44;
            
            N(5, 0) = 0.44;
            N(5, 1) = one_eight;
            N(5, 2) = 0.44;
            
            N(6, 0) = 0.44;
            N(6, 1) = 0.44;
            N(6, 2) = one_eight;
			
			
            //first
            pos(0, 0) = one_eight * geom[0].X() + one_eight * geom[1].X() + three_quarters * geom[2].X();
            pos(0, 1) = one_eight * geom[0].Y() + one_eight * geom[1].Y() + three_quarters * geom[2].Y();
            pos(0, 2) = one_eight * geom[0].Z() + one_eight * geom[1].Z() + three_quarters * geom[2].Z();
	    
            //second
            pos(1, 0) = three_quarters * geom[0].X() + one_eight * geom[1].X() + one_eight * geom[2].X();
            pos(1, 1) = three_quarters * geom[0].Y() + one_eight * geom[1].Y() + one_eight * geom[2].Y();
            pos(1, 2) = three_quarters * geom[0].Z() + one_eight * geom[1].Z() + one_eight * geom[2].Z();

            //third
            pos(2, 0) = one_eight * geom[0].X() + three_quarters * geom[1].X() + one_eight * geom[2].X();
            pos(2, 1) = one_eight * geom[0].Y() + three_quarters * geom[1].Y() + one_eight * geom[2].Y();
            pos(2, 2) = one_eight * geom[0].Z() + three_quarters * geom[1].Z() + one_eight * geom[2].Z();

            //fourth
			pos(3, 0) = one_third * geom[0].X() + one_third * geom[1].X() + one_third * geom[2].X();
            pos(3, 1) = one_third * geom[0].Y() + one_third * geom[1].Y() + one_third * geom[2].Y();
            pos(3, 2) = one_third * geom[0].Z() + one_third * geom[1].Z() + one_third * geom[2].Z();
            
            //fifth
            pos(4, 0) = one_eight * geom[0].X() + 0.44 * geom[1].X() + 0.44 * geom[2].X();
            pos(4, 1) = one_eight * geom[0].Y() + 0.44 * geom[1].Y() + 0.44 * geom[2].Y();
            pos(4, 2) = one_eight * geom[0].Z() + 0.44 * geom[1].Z() + 0.44 * geom[2].Z();
	    
            //sixth
            pos(5, 0) = 0.44 * geom[0].X() + one_eight * geom[1].X() + 0.44 * geom[2].X();
            pos(5, 1) = 0.44 * geom[0].Y() + one_eight * geom[1].Y() + 0.44 * geom[2].Y();
            pos(5, 2) = 0.44 * geom[0].Z() + one_eight * geom[1].Z() + 0.44 * geom[2].Z();

            //seventh
            pos(6, 0) = 0.44 * geom[0].X() + 0.44 * geom[1].X() + one_eight * geom[2].X();
            pos(6, 1) = 0.44 * geom[0].Y() + 0.44 * geom[1].Y() + one_eight * geom[2].Y();
            pos(6, 2) = 0.44 * geom[0].Z() + 0.44 * geom[1].Z() + one_eight * geom[2].Z();

            
            

        }
        
        void ComputeGaussPointPositionsForPostReseed(Geometry< Node < 3 > >& geom, boost::numeric::ublas::bounded_matrix<double, 9, 3 > & pos,boost::numeric::ublas::bounded_matrix<double, 9, 4 > & N) //3D
        {
            double one_quarter = 0.25;
            double small_fraction = 0.1; //1.0 / 6.0;
            double big_fraction = 0.7; //2.0 * one_third;
            double mid_fraction = 0.3; //2.0 * one_third;

            N(0, 0) = big_fraction;
            N(0, 1) = small_fraction;
            N(0, 2) = small_fraction;
            N(0, 3) = small_fraction;
            
            N(1, 0) = small_fraction;
            N(1, 1) = big_fraction;
            N(1, 2) = small_fraction;
            N(1, 3) = small_fraction;
            
            N(2, 0) = small_fraction;
            N(2, 1) = small_fraction;
            N(2, 2) = big_fraction;
            N(2, 3) = small_fraction;
            
            N(3, 0) = small_fraction;
            N(3, 1) = small_fraction;
            N(3, 2) = small_fraction;
            N(3, 3) = big_fraction;
			
			N(4, 0) = one_quarter;
			N(4, 1) = one_quarter;
			N(4, 2) = one_quarter;
			N(4, 3) = one_quarter;
			
			N(5, 0) = small_fraction;
            N(5, 1) = mid_fraction;
            N(5, 2) = mid_fraction;
            N(5, 3) = mid_fraction;
            
			N(6, 0) = mid_fraction;
            N(6, 1) = small_fraction;
            N(6, 2) = mid_fraction;
            N(6, 3) = mid_fraction;
            
			N(7, 0) = mid_fraction;
            N(7, 1) = mid_fraction;
            N(7, 2) = small_fraction;
            N(7, 3) = mid_fraction;
            
            N(8, 0) = mid_fraction;
            N(8, 1) = mid_fraction;
            N(8, 2) = mid_fraction;
            N(8, 3) = small_fraction;
			
			pos=ZeroMatrix(9,3);
            for (unsigned int i=0; i!=4; i++) //going through the 4 nodes
            {
				array_1d<double, 3 > & coordinates = geom[i].Coordinates();
				for (unsigned int j=0; j!=9; j++) //going through the 9 particles
				{
					for (unsigned int k=0; k!=3; k++) //x,y,z
						pos(j,k) += N(j,i) * coordinates[k];
				} 
			}
			
			
        }
        
        
        
        void ComputeGaussPointPositionsForPreReseed(Geometry< Node < 3 > >& geom, boost::numeric::ublas::bounded_matrix<double, 3, 3 > & pos,boost::numeric::ublas::bounded_matrix<double, 3, 3 > & N) //2D
        {

            N(0, 0) = 0.5;
            N(0, 1) = 0.25;
            N(0, 2) = 0.25;
            
            N(1, 0) = 0.25;
            N(1, 1) = 0.5;
            N(1, 2) = 0.25;
            
            N(2, 0) = 0.25;
            N(2, 1) = 0.25;
            N(2, 2) = 0.5;

            //first
            pos(0, 0) = 0.5 * geom[0].X() + 0.25 * geom[1].X() + 0.25 * geom[2].X();
            pos(0, 1) = 0.5 * geom[0].Y() + 0.25 * geom[1].Y() + 0.25 * geom[2].Y();
            pos(0, 2) = 0.5 * geom[0].Z() + 0.25 * geom[1].Z() + 0.25 * geom[2].Z();
	    
            //second
            pos(1, 0) = 0.25 * geom[0].X() + 0.5 * geom[1].X() + 0.25 * geom[2].X();
            pos(1, 1) = 0.25 * geom[0].Y() + 0.5 * geom[1].Y() + 0.25 * geom[2].Y();
            pos(1, 2) = 0.25 * geom[0].Z() + 0.5 * geom[1].Z() + 0.25 * geom[2].Z();

            //third
            pos(2, 0) = 0.25 * geom[0].X() + 0.25 * geom[1].X() + 0.5 * geom[2].X();
            pos(2, 1) = 0.25 * geom[0].Y() + 0.25 * geom[1].Y() + 0.5 * geom[2].Y();
            pos(2, 2) = 0.25 * geom[0].Z() + 0.25 * geom[1].Z() + 0.5 * geom[2].Z();

        }
        
        void ComputeGaussPointPositionsForPreReseed(Geometry< Node < 3 > >& geom, boost::numeric::ublas::bounded_matrix<double, 4, 3 > & pos,boost::numeric::ublas::bounded_matrix<double, 4, 4 > & N) //3D
        {
			//creating 4 particles, each will be closer to a node and equidistant to the other nodes
			

            N(0, 0) = 0.4;
            N(0, 1) = 0.2;
            N(0, 2) = 0.2;
            N(0, 3) = 0.2;
            
            N(1, 0) = 0.2;
            N(1, 1) = 0.4;
            N(1, 2) = 0.2;
            N(1, 3) = 0.2;
            
            N(2, 0) = 0.2;
            N(2, 1) = 0.2;
            N(2, 2) = 0.4;
            N(2, 3) = 0.2;
            
            N(3, 0) = 0.2;
            N(3, 1) = 0.2;
            N(3, 2) = 0.2;
            N(3, 3) = 0.4;

            pos=ZeroMatrix(4,3);
            for (unsigned int i=0; i!=4; i++) //going through the 4 nodes
            {
				array_1d<double, 3 > & coordinates = geom[i].Coordinates();
				for (unsigned int j=0; j!=4; j++) //going through the 4 particles
				{
					for (unsigned int k=0; k!=3; k++) //x,y,z
						pos(j,k) += N(j,i) * coordinates[k];
				} 
			}

        }
        

		
		void ComputeGaussPointPositions_45(Geometry< Node < 3 > >& geom, boost::numeric::ublas::bounded_matrix<double, 45, 3 > & pos,boost::numeric::ublas::bounded_matrix<double, 45, 3 > & N)
        {
			//std::cout << "NEW ELEMENT" << std::endl;
			unsigned int counter=0;
			for (unsigned int i=0; i!=9;i++)
			{
				for (unsigned int j=0; j!=(9-i);j++)
				{
					N(counter,0)=0.05+double(i)*0.1;
					N(counter,1)=0.05+double(j)*0.1;
					N(counter,2)=1.0 - ( N(counter,1)+ N(counter,0) ) ;
					pos(counter, 0) = N(counter,0) * geom[0].X() + N(counter,1) * geom[1].X() + N(counter,2) * geom[2].X();
					pos(counter, 1) = N(counter,0) * geom[0].Y() + N(counter,1) * geom[1].Y() + N(counter,2) * geom[2].Y();
					pos(counter, 2) = N(counter,0) * geom[0].Z() + N(counter,1) * geom[1].Z() + N(counter,2) * geom[2].Z();
					//std::cout << N(counter,0) << " " << N(counter,1) << " " << N(counter,2) << " " << std::endl;
					counter++;
					
				}
			}
					
        }
        
        void ComputeGaussPointPositions_initial(Geometry< Node < 3 > >& geom, boost::numeric::ublas::bounded_matrix<double, 15, 3 > & pos,boost::numeric::ublas::bounded_matrix<double, 15, 3 > & N) //2D
        {
			//std::cout << "NEW ELEMENT" << std::endl;
			unsigned int counter=0;
			for (unsigned int i=0; i!=5;i++)
			{
				for (unsigned int j=0; j!=(5-i);j++)
				{
					N(counter,0)=0.05+double(i)*0.2;
					N(counter,1)=0.05+double(j)*0.2;
					N(counter,2)=1.0 - ( N(counter,1)+ N(counter,0) ) ;
					pos(counter, 0) = N(counter,0) * geom[0].X() + N(counter,1) * geom[1].X() + N(counter,2) * geom[2].X();
					pos(counter, 1) = N(counter,0) * geom[0].Y() + N(counter,1) * geom[1].Y() + N(counter,2) * geom[2].Y();
					pos(counter, 2) = N(counter,0) * geom[0].Z() + N(counter,1) * geom[1].Z() + N(counter,2) * geom[2].Z();
					//std::cout << N(counter,0) << " " << N(counter,1) << " " << N(counter,2) << " " << std::endl;
					counter++;
					
				}
			}
					
        }
        
        void ComputeGaussPointPositions_initial(Geometry< Node < 3 > >& geom, boost::numeric::ublas::bounded_matrix<double, 20, 3 > & pos,boost::numeric::ublas::bounded_matrix<double, 20, 4 > & N) //3D
        {
			//std::cout << "NEW ELEMENT" << std::endl;
			//double total;
			double fraction_increment;
			unsigned int counter=0;
			for (unsigned int i=0; i!=4;i++) //going to build a particle "pyramid"(tetrahedra) by layers. the first layer will be made by a triangle of 4 base X 4 height. since it is a triangle, it means it will have 10 particles
			{
				//std::cout << "inside i" <<  i << std::endl;
				for (unsigned int j=0; j!=(4-i);j++)
				{
					//std::cout << "inside j" << j << std::endl;
					for (unsigned int k=0; k!=(4-i-j);k++)
					{
						//std::cout << "inside k" << k << std::endl;
						N(counter,0)= 0.27 * ( 0.175 + double(i) ) ; //this is our "surface" in which we will build each layer, so we must construct a triangle using what's left of the shape functions total (a total of 1)
						
						//total = 1.0 - N(counter,0); 
						fraction_increment = 0.27; // 
						
						N(counter,1)=fraction_increment * (0.175 + double(j));
						N(counter,2)=fraction_increment * (0.175 + double(k));
						N(counter,3)=1.0 - ( N(counter,0)+ N(counter,1) + N(counter,2) ) ;
						pos(counter, 0) = N(counter,0) * geom[0].X() + N(counter,1) * geom[1].X() + N(counter,2) * geom[2].X() + N(counter,3) * geom[3].X();
						pos(counter, 1) = N(counter,0) * geom[0].Y() + N(counter,1) * geom[1].Y() + N(counter,2) * geom[2].Y() + N(counter,3) * geom[3].Y();
						pos(counter, 2) = N(counter,0) * geom[0].Z() + N(counter,1) * geom[1].Z() + N(counter,2) * geom[2].Z() + N(counter,3) * geom[3].Z();
						//std::cout << N(counter,0) << " " << N(counter,1) << " " << N(counter,2) << " " << std::endl;
						counter++;
					}
					
				}
			}
					
        }
     
        
        // Bubble Sort Function for Descending Order 
		void BubbleSort(array_1d<double,7> &distances , array_1d<int,7 > &positions, unsigned int & arrange_number)
		{
			  int i, j;
			  bool flag = true;    // set flag to 1 to start first pass
			  double temp;             // holding variable
			  int temp_position;
			  int numLength = arrange_number; 
			  for(i = 1; (i <= numLength) && flag; i++)
			  {
				  flag = false;
				  for (j=0; j < (numLength -1); j++)
				 {
					   if (distances[j+1] < distances[j])      // descending order simply changes to >
					  { 
							temp = distances[j];             // swap elements
							distances[j] = distances[j+1];
							distances[j+1] = temp;
							
							temp_position = positions[j];  //swap positions
							positions[j] = positions[j+1];
							positions[j+1] = temp_position;
							
							flag = true;               // indicates that a swap occurred.
					   }
				  }
			  }
			  return;   //arrays are passed to functions by address; nothing is returned
		}
		
		
		
		void BubbleSort(array_1d<double,9> &distances , array_1d<int,9 > &positions, unsigned int & arrange_number)
		{
			  int i, j;
			  bool flag = true;    // set flag to 1 to start first pass
			  double temp;             // holding variable
			  int temp_position;
			  int numLength = arrange_number; 
			  for(i = 1; (i <= numLength) && flag; i++)
			  {
				  flag = false;
				  for (j=0; j < (numLength -1); j++)
				 {
					   if (distances[j+1] < distances[j])      // descending order simply changes to >
					  { 
							temp = distances[j];             // swap elements
							distances[j] = distances[j+1];
							distances[j+1] = temp;
							
							temp_position = positions[j];  //swap positions
							positions[j] = positions[j+1];
							positions[j+1] = temp_position;
							
							flag = true;               // indicates that a swap occurred.
					   }
				  }
			  }
			  return;   //arrays are passed to functions by address; nothing is returned
		}
		
	ModelPart& mr_model_part;
	ModelPart::Pointer mtopographic_model_part_pointer;
	array_1d<double, 3 > mcalculation_domain_complete_displacement;
	array_1d<double, 3 > mcalculation_domain_added_displacement;
	bool mintialized_transfer_tool;
	int m_nparticles;
	int mnelems;
	double mDENSITY_WATER;
	double mDENSITY_AIR;
	double mINV_DENSITY_WATER;
	double mINV_DENSITY_AIR;
	double mDENSITY_AIRoverDENSITY_WATER;
	//vector<double> mareas_vector; UNUSED SO COMMENTED 
	int max_nsubsteps;
	double max_substep_dt;
	unsigned int mmaximum_number_of_particles;
	//matrix< PFEM_Particle_Fluid* > mparticles_ordered_pointers; 
	std::vector< PFEM_Particle_Fluid  > mparticles_vector; //Point<3>
	unsigned int mlast_elem_id;
	bool modd_timestep;
	//ModelPart& mr_particle_model_part;
	
	vector<int> mnumber_of_particles_in_elems; 
	vector<int> mnumber_of_particles_in_elems_aux; 
	vector<int> mnumber_of_water_particles_in_elems; 
	vector<ParticlePointerVector*>  mpointers_to_particle_pointers_vectors;
	//std::vector< PointerVector< ParticlePointerVector >* > 
	double mwrite_particle;
	double mload_particle;
	double mload_current_ParticlePointers;
	double mload_old_ParticlePointers;
	double mdo_nothing_per_particle;
	//KRATOS_DEFINE_VARIABLE(PointerVector< ParticlePointerVector > , PARTICLE_POINTERS)
	
	//PointerVector< ParticlePointerVector >&  particle_pointers =  (ielem->GetValue(PARTICLE_POINTERS));
	
	
	
	typename BinsObjectDynamic<Configure>::Pointer  mpBinsObjectDynamic;
	typename BinsObjectDynamic<Configure>::Pointer  mpTopographicBinsObjectDynamic;

	};
	
}  // namespace Kratos.

#endif // KRATOS_MOVE_PART_UTILITY_DIFF2_INCLUDED  defined 


