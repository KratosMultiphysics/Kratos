#if !defined(KRATOS_CALCULATE_WATER_FRACTION_UTILITY_INCLUDED )
#define  KRATOS_CALCULATE_WATER_FRACTION_UTILITY_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <algorithm>

// Project includes
#include "includes/define.h"
#include "pfem_2_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/element.h"
#include "utilities/enrichment_utilities.h"


namespace Kratos
{
	template< unsigned int TDim>
	class CalculateWaterFraction
	{
	public:

		KRATOS_CLASS_POINTER_DEFINITION(CalculateWaterFraction);

		CalculateWaterFraction(ModelPart& model_part)
			: mr_model_part(model_part)              //mr_model_part is saved as private variable (declared at the end of the file)
		{
			KRATOS_TRY
			//std::cout << "Hello, I am the constructor of the Utility" << std::endl;
			KRATOS_CATCH("")
		}


		~CalculateWaterFraction()
		{}
		/*
		double Calculate() //water fraction
		{
			KRATOS_TRY
			//double area;                    //we create the needed variables
			double sum_area=0.0;
			double sum_water_area=0.0;
			//double one_third=1.0/3.0;
			ModelPart::NodesContainerType::iterator inodebegin = mr_model_part.NodesBegin();

			#pragma omp parallel for reduction(+:sum_area) reduction(+:sum_water_area)
			for(unsigned int ii=0; ii<mr_model_part.Nodes().size(); ii++)
			{
				ModelPart::NodesContainerType::iterator inode = inodebegin+ii;
				const double & nodal_area = inode->FastGetSolutionStepValue(NODAL_AREA); //resetting the temperature
				sum_area += nodal_area;
				if ((inode->FastGetSolutionStepValue(DISTANCE))<0.0)
					sum_water_area += nodal_area;
			}
			const double water_fraction = sum_water_area / sum_area;
			//std::cout << "Finished, the mean temperature is" << water_fraction << std::endl;   //we print the result
			return water_fraction;

			KRATOS_CATCH("")
		}
		*/
		double Calculate()
		{
			KRATOS_TRY


			double sum_areas=1.0e-100;
			//double sum_temperatures=0.0;
			//double nodal_weight=1.0/(1.0+double(TDim));

			ModelPart::ElementsContainerType::iterator ielembegin = mr_model_part.ElementsBegin();


			vector<unsigned int> element_partition;
			#ifdef _OPENMP
				int number_of_threads = omp_get_max_threads();
			#else
				int number_of_threads = 1;
			#endif
			OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Elements().size(), element_partition);

			#pragma omp parallel for reduction(+:sum_areas)
			for(int kkk=0; kkk<number_of_threads; kkk++)
			{
				double  thread_sum_areas=0.0;
				for(unsigned int ii=element_partition[kkk]; ii<element_partition[kkk+1]; ii++)
				{
					ModelPart::ElementsContainerType::iterator ielem = ielembegin+ii;
					double Area;
					BoundedMatrix<double, (TDim+1), TDim > DN_DX;
					array_1d<double, (TDim+1) > N;
					Geometry<Node<3> >& geom = ielem->GetGeometry();
					GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Area);
					//sum_areas+=area;
					int negative_nodes=0;
					int positive_nodes=0;
					for (unsigned int k = 0; k < (TDim+1); k++)
					{
						if(geom[k].FastGetSolutionStepValue(DISTANCE)<0.0)
							negative_nodes++;
						else
							positive_nodes++;
					}
					if (negative_nodes==(TDim+1))
						thread_sum_areas+=Area;

					else if (negative_nodes>0)
					{
						array_1d<double,(TDim+1)>  distances;
						for (unsigned int i = 0; i < (TDim+1); i++)
						{
							distances[i] = geom[i].FastGetSolutionStepValue(DISTANCE);
						}


						BoundedMatrix<double,3*(TDim-1), 2> Nenriched;
						array_1d<double,(3*(TDim-1))> volumes;
						BoundedMatrix<double,(TDim+1), TDim > coords;
						BoundedMatrix<double, 3*(TDim-1), (TDim+1) > Ngauss;
						array_1d<double,(3*(TDim-1))> signs;
						std::vector< Matrix > gauss_gradients(3*(TDim-1));
						//fill coordinates

						//unsigned int single_triangle_node = 1;
							for (unsigned int i = 0; i < (TDim+1); i++)
							{
								const array_1d<double, 3 > & xyz = geom[i].Coordinates();
								for (unsigned int j = 0; j < TDim; j++)
									coords(i,j)=xyz(j);
							}
						for (unsigned int i = 0; i < 3*(TDim-1); i++)
							gauss_gradients[i].resize(2, TDim, false);  //2 values of the 2 shape functions, and derivates in (xy) direction).
						unsigned int ndivisions = EnrichmentUtilities::CalculateEnrichedShapeFuncions(coords, DN_DX, distances, volumes, Ngauss, signs, gauss_gradients, Nenriched);
						for (unsigned int i=0;i!=ndivisions;i++)
							if (signs(i)<0.0) thread_sum_areas+=volumes(i);


					}
				}
				sum_areas = thread_sum_areas;

			}
			//const double mean_temperature = sum_temperatures / sum_areas;
			std::cout << "Finished, the water volume is " << sum_areas << std::endl;
			return sum_areas;

			KRATOS_CATCH("")
		}

 		double CalculateWaterHeight(double x_position)
		{
			KRATOS_TRY

			double all_threads_water_height=-100000000.0;
			const double tolerance=0.001;
			const double upper_limit=x_position+tolerance;
			const double lower_limit=x_position-tolerance;
			ModelPart::NodesContainerType::iterator inodebegin = mr_model_part.NodesBegin();

			vector<unsigned int> node_partition;
			#ifdef _OPENMP
				int number_of_threads = omp_get_max_threads();
			#else
				int number_of_threads = 1;
			#endif
			OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Nodes().size(), node_partition);

			#pragma omp parallel for
			for(int kkk=0; kkk<number_of_threads; kkk++)
			{
				double local_thread_water_height=-100000000.0;
				for(unsigned int ii=node_partition[kkk]; ii<node_partition[kkk+1]; ii++)
				{
					ModelPart::NodesContainerType::iterator inode = inodebegin+ii;
					double & distance = (inode->FastGetSolutionStepValue(DISTANCE));
					if ( distance <0.0)
					{
						if ((inode->X())<upper_limit && (inode->X())>lower_limit && (inode->Y())>local_thread_water_height)
							local_thread_water_height = (inode->Y());
					}
					//now we search for the node given certain criteria

				}


				if ( local_thread_water_height > all_threads_water_height )
				{
					#pragma omp critical
					{
						if ( local_thread_water_height > all_threads_water_height ) all_threads_water_height = local_thread_water_height;
					}
				}
			}

			return all_threads_water_height;

			KRATOS_CATCH("")
		}


		double CalculateMaxCourant()
		{
			KRATOS_TRY

			ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
			const double delta_t = CurrentProcessInfo[DELTA_TIME];

			double all_threads_max_courant = 0.0;
			ModelPart::ElementsContainerType::iterator ielembegin = mr_model_part.ElementsBegin();

			vector<unsigned int> element_partition;
			#ifdef _OPENMP
				int number_of_threads = omp_get_max_threads();
			#else
				int number_of_threads = 1;
			#endif
			OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Elements().size(), element_partition);

			#pragma omp parallel for
			for(int kkk=0; kkk<number_of_threads; kkk++)
			{
				double local_thread_max_courant = 0.0;
				for(unsigned int ii=element_partition[kkk]; ii<element_partition[kkk+1]; ii++)
				{
					ModelPart::ElementsContainerType::iterator ielem = ielembegin+ii;
					//double & distance = (inode->FastGetSolutionStepValue(DISTANCE));
					if ((ielem->GetValue(VELOCITY_OVER_ELEM_SIZE))>local_thread_max_courant)
						local_thread_max_courant = (ielem->GetValue(VELOCITY_OVER_ELEM_SIZE));
				}


				if ( local_thread_max_courant > all_threads_max_courant )
				{
					#pragma omp critical
					{
						if ( local_thread_max_courant > all_threads_max_courant ) all_threads_max_courant = local_thread_max_courant;
					}
				}
			}

			all_threads_max_courant *= delta_t * 1.414;

			return all_threads_max_courant;

			KRATOS_CATCH("")
		}

		double CalculateMaxCourantInNegativeElements()
		{
			KRATOS_TRY

			//using a nodal approach (faster!)
			ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
			const double delta_t = CurrentProcessInfo[DELTA_TIME];

			double all_threads_max_courant = 0.0;
			ModelPart::NodesContainerType::iterator inodebegin = mr_model_part.NodesBegin();

			vector<unsigned int> node_partition;
			#ifdef _OPENMP
				int number_of_threads = omp_get_max_threads();
			#else
				int number_of_threads = 1;
			#endif
			OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Nodes().size(), node_partition);

			#pragma omp parallel for
			for(int kkk=0; kkk<number_of_threads; kkk++)
			{
				double local_thread_max_courant = 0.0;
				for(unsigned int ii=node_partition[kkk]; ii<node_partition[kkk+1]; ii++)
				{
					ModelPart::NodesContainerType::iterator inode = inodebegin+ii;
					const double & distance = (inode->FastGetSolutionStepValue(DISTANCE));
					const double velocity = sqrt(pow(inode->FastGetSolutionStepValue(VELOCITY_X),2)+pow(inode->FastGetSolutionStepValue(VELOCITY_Y),2)+pow(inode->FastGetSolutionStepValue(VELOCITY_Z),2));
					const double nodal_courant =  (velocity*delta_t/inode->FastGetSolutionStepValue(MEAN_SIZE));
					if(nodal_courant>local_thread_max_courant && distance < 0.0) //only for negative nodes!
						local_thread_max_courant = nodal_courant;
				}


				if ( local_thread_max_courant > all_threads_max_courant )
				{
					#pragma omp critical
					{
						if ( local_thread_max_courant > all_threads_max_courant ) all_threads_max_courant = local_thread_max_courant;
					}
				}
			}

			//all_threads_max_courant *= delta_t * 1.414;

			return all_threads_max_courant;

			KRATOS_CATCH("")
		}

		double CalculateMeanCourant() //water fraction
		{
			KRATOS_TRY

			ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
			const double delta_t = CurrentProcessInfo[DELTA_TIME];

			//double area=0.0;                    //we create the needed variables
			//double  number_of_threads = double(OpenMPUtils::GetNumThreads());
			double sum_courant=0.0;

			ModelPart::ElementsContainerType::iterator ielembegin = mr_model_part.ElementsBegin();

			vector<unsigned int> element_partition;
			#ifdef _OPENMP
				int number_of_threads = omp_get_max_threads();
			#else
				int number_of_threads = 1;
			#endif
			OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Elements().size(), element_partition);

			#pragma omp parallel for reduction(+:sum_courant)
			for(int kkk=0; kkk<number_of_threads; kkk++)
			{
				double thread_sum_courant=0.0;
				for(unsigned int ii=element_partition[kkk]; ii<element_partition[kkk+1]; ii++)
				{
					ModelPart::ElementsContainerType::iterator ielem = ielembegin+ii;
					//const double & velocity_over_elem_size = (ielem->GetValue(VELOCITY_OVER_ELEM_SIZE));
					if ((ielem->GetValue(VELOCITY_OVER_ELEM_SIZE))>0.0)
						thread_sum_courant += ielem->GetValue(VELOCITY_OVER_ELEM_SIZE);
				}
				sum_courant += thread_sum_courant;
			}

			sum_courant *= delta_t * 1.414 / double(mr_model_part.Elements().size());

			return sum_courant;

			KRATOS_CATCH("")
		}

				//NOW ONLY VISCOUS. but since in the first step we cannot use the pressure we just add the viscoust forces. still, lines to use pressure can be uncommented
		double CalculateForce(int direction) //
		{
			KRATOS_TRY

			ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
			double viscosity = CurrentProcessInfo[VISCOSITY];
			//double delta_t = CurrentProcessInfo[DELTA_TIME];
			//array_1d<double,3> & gravity= CurrentProcessInfo[GRAVITY];


	        const array_1d<double,3> zero3 = ZeroVector(3);
			double nodal_weight = 1.0/ (double (TDim) );

			double force=0.0;

			ModelPart::ElementsContainerType::iterator ielembegin = mr_model_part.ElementsBegin();


			vector<unsigned int> element_partition;
			#ifdef _OPENMP
				int number_of_threads = omp_get_max_threads();
			#else
				int number_of_threads = 1;
			#endif
			OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Elements().size(), element_partition);

			#pragma omp parallel for reduction(+:force)
			for(int kkk=0; kkk<number_of_threads; kkk++)
			{
				double  thread_force=0.0;
				for(unsigned int ii=element_partition[kkk]; ii<element_partition[kkk+1]; ii++)
				{
					ModelPart::ElementsContainerType::iterator ielem = ielembegin+ii;
					if (ielem->Is(ACTIVE)) //elements can be inactive to add temporary walls. fractional velocity is integrated by parts so walls are seen as having zero velocity without doing anything
					{
						//double Area;
						Geometry<Node<3> >& geom = ielem->GetGeometry();

						array_1d<unsigned int, 4 > fixed_nodes; //unordered : i position in the array might not correspond to i node of the element
						array_1d<bool, 4 > is_node_fixed; //i position belongs to i node of the element
						unsigned int number_of_fixed_nodes=0;
						//bool boundary_element=false;
						BoundedMatrix<double, 4, 3 > velocities = ZeroMatrix(4, 3);

						for (unsigned int i = 0; i < (TDim+1); i++)
						{
							const array_1d<double, 3 > & velocity = geom[i].FastGetSolutionStepValue(VELOCITY);
							for (unsigned int j = 0; j < (TDim); j++)
								velocities(i,j) = velocity[j];

							if (TDim==2)
							{
								if (geom[i].IsFixed(FRACT_VEL_X) && geom[i].IsFixed(FRACT_VEL_Y))
								{
									fixed_nodes[number_of_fixed_nodes]=i;
									is_node_fixed[i]=true;
									number_of_fixed_nodes++;
								}
								else
									is_node_fixed[i]=false;
							}
							else // (TDim==3)
							{
								if (geom[i].IsFixed(FRACT_VEL_X) && geom[i].IsFixed(FRACT_VEL_Y) && geom[i].IsFixed(FRACT_VEL_Z) )
								{
									fixed_nodes[number_of_fixed_nodes]=i;
									number_of_fixed_nodes++;
									is_node_fixed[i]=true;
								}
								else
									is_node_fixed[i]=false;
							}
						}

						double plane_point_distance=1.0;
						double fixed_face_area_or_lenght=0.0;
						array_1d<double, 3 > boundary_stress;
						if (number_of_fixed_nodes==TDim) //it means we have cutted elements!
						{
							//boundary_element=true;
							array_1d<double, 3 > normal;
							unsigned int free_node=0;
							if (TDim==2)
							{
								fixed_face_area_or_lenght = fabs(sqrt(pow((geom[fixed_nodes[1]].Y()-geom[fixed_nodes[0]].Y()),2 ) + pow( (geom[fixed_nodes[1]].X()-geom[fixed_nodes[0]].X() ),2 ) ) );
								normal[0] = geom[fixed_nodes[1]].Y()-geom[fixed_nodes[0]].Y();
								normal[1] = - ( geom[fixed_nodes[1]].X()-geom[fixed_nodes[0]].X() );
								normal[2] = 0.0;
								normal /= sqrt(normal[0]*normal[0]+normal[1]*normal[1]);

								if (fixed_nodes[0]==0)
								{
									if (fixed_nodes[1]==1)
										free_node=2;
									else
										free_node=1;
								}
								else
									free_node=0;

								//the plane is composed by the unit normal and any of the points of fixed nodes. we will use fixed_nodes[0];
								plane_point_distance = inner_prod( (geom[free_node].Coordinates()-geom[fixed_nodes[0]].Coordinates()) , normal);
								//boundary_stress = geom[free_node].FastGetSolutionStepValue(VELOCITY)*viscosity/plane_point_distance;
								if (plane_point_distance<0.0)
								{
									plane_point_distance*=-1.0;
									normal *= -1.0;
								}
							}
							else //(TDim==3)
							{
								//the area is obtained from the crossproduct of the 2 vertices:
								MathUtils<double>::CrossProduct(normal, geom[fixed_nodes[1]].Coordinates() - geom[fixed_nodes[0]].Coordinates(), geom[fixed_nodes[2]].Coordinates() - geom[fixed_nodes[0]].Coordinates() );
								fixed_face_area_or_lenght = 0.5 * sqrt( pow(normal[0],2) + pow(normal[1],2) + pow(normal[2],2) );
								normal /= 2.0 * fixed_face_area_or_lenght;  //this way it is a unit vector. now we must find the distance from the plane generated by the triangles to the free node:

								//fixed_face_area_or_lenght = fabs(fixed_face_area_or_lenght);

								for (unsigned int j=0; j!=(TDim+1); j++)
								{
									if (is_node_fixed[j]==false)
									{
										free_node=j;
										break;
									}
								}

								//the plane is composed by the unit normal and any of the points of fixed nodes. we will use fixed_nodes[0];
								plane_point_distance = inner_prod( (geom[free_node].Coordinates()-geom[fixed_nodes[0]].Coordinates()) , normal);
								if (plane_point_distance<0.0)
									normal *= -1.0;
								{
									plane_point_distance*=-1.0;
									normal *= -1.0;
								}
								//boundary_stress = geom[free_node].FastGetSolutionStepValue(VELOCITY)*viscosity/plane_point_distance;
							}

							boundary_stress = - geom[free_node].FastGetSolutionStepValue(VELOCITY)*viscosity/(fabs(plane_point_distance));
							//KRATOS_WATCH(plane_point_distance)
							//KRATOS_WATCH(boundary_stress)
							//KRATOS_WATCH(fixed_face_area_or_lenght)

							//drag forces:
							thread_force += boundary_stress[direction]*fixed_face_area_or_lenght; // unit density! careful!
							//face_force+=fixed_face_area_or_lenght*normal[direction];
							//now pressure forces:

							for (unsigned int j=0; j!=(TDim); j++) // the 2 or 3 nodes that define the fixed face:
							{
								/*
								if ( (geom[fixed_nodes[j]].X())<5.0 )
									face_force += nodal_weight*(0.5*fixed_face_area_or_lenght*normal[direction]);
								else
									face_force -= nodal_weight*(0.5*fixed_face_area_or_lenght*normal[direction]);
								*/
								thread_force +=nodal_weight*(geom[fixed_nodes[j]].FastGetSolutionStepValue(PRESSURE))*fixed_face_area_or_lenght*normal[direction];

								//array_1d<double,3> & nodal_normal= (geom[fixed_nodes[j]].FastGetSolutionStepValue(NORMAL));
								//face_force += (geom[fixed_nodes[j]].FastGetSolutionStepValue(PRESSURE))*nodal_normal[direction];
							}




						}
					}

				}

				force+=thread_force;
			}

			return force;

			KRATOS_CATCH("")
		}



	protected:


	private:

		ModelPart& mr_model_part;

	};

}  // namespace Kratos.

#endif // KRATOS_CALCULATE__WATER_FRACTION_UTILITY_INCLUDED  defined
