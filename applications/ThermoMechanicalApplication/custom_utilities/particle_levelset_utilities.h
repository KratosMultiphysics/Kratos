/*
==============================================================================
KratosTestApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2010
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

The  above  copyright  notice  and  this permission  notice  sKRATOS_WATCH(disp);hall  be
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
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2007-03-06 10:30:31 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_PARTICLE_LEVELSET_UTILITIES_INCLUDED )
#define  KRATOS_PARTICLE_LEVELSET_UTILITIES_INCLUDED


// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes
#include "pybind11/pybind11.h"

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/deprecated_variables.h"
#include "includes/node.h"
#include "utilities/geometry_utilities.h"
#include "geometries/tetrahedra_3d_4.h"
#include "thermo_mechanical_application.h"
#include "spatial_containers/spatial_containers.h"
#include "utilities/timer.h"
#include "processes/node_erase_process.h"
#include "utilities/binbased_fast_point_locator.h"

#include "utilities/timer.h"

// #include <boost/random/linear_congruential.hpp>
// #include <boost/random/uniform_int.hpp>
// #include <boost/random/uniform_real.hpp>
// #include <boost/random/variate_generator.hpp>
// #include <boost/generator_iterator.hpp>
// #include <tr1/random>
#include <time.h>

#ifdef _OPENMP
#include "omp.h"
#endif



namespace Kratos
{

  template<std::size_t TDim> class ParticleLevelSetUtils
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ParticleLevelSetUtils<TDim>);

    //**********************************************************************************************
    //**********************************************************************************************
    //function to seed a list of new nodes

    void Seed(ModelPart& rEulerianModelPart, ModelPart& rLagrangianModelPart, const double max_seed_distance, const double min_edge_size)
    {
        KRATOS_TRY;

        rLagrangianModelPart.Nodes().clear();

	unsigned int ele_id = 1;
        for (ModelPart::ElementsContainerType::iterator el_it = rEulerianModelPart.ElementsBegin();
                el_it != rEulerianModelPart.ElementsEnd(); el_it++)
        {
	  el_it->SetId(ele_id);
	  ele_id++;
	}

      if(TDim==2){
        BoundedMatrix<double, 16, 3 > pos;
        BoundedMatrix<double, 16, 3 > N;

        CreateParticles2D(rEulerianModelPart,rLagrangianModelPart,pos,N,max_seed_distance,min_edge_size);

        }
      else
      {
// 	BoundedMatrix<double, 56, 3 > pos;
//         BoundedMatrix<double, 56, 4 > N;
//         CreateParticles3D(rEulerianModelPart,rLagrangianModelPart,pos,N,max_seed_distance,min_edge_size);

	BoundedMatrix<double, 10, 3 > pos;
        BoundedMatrix<double, 10, 4 > N;
        FewCreateParticles3D(rEulerianModelPart,rLagrangianModelPart,pos,N,max_seed_distance,min_edge_size);
      }



        for (ModelPart::NodesContainerType::iterator node_it = rLagrangianModelPart.NodesBegin();
                node_it != rLagrangianModelPart.NodesEnd(); node_it++)
        {
            node_it->FastGetSolutionStepValue(VELOCITY, 1) = node_it->FastGetSolutionStepValue(VELOCITY);
// 	    node_it->FastGetSolutionStepValue(DISTANCE, 1) = node_it->FastGetSolutionStepValue(DISTANCE);
        }


        KRATOS_CATCH("");
    }

    //**********************************************************************************************
    //**********************************************************************************************

    void StreamlineMove(const double dt, ModelPart& rEulerianModelPart, ModelPart& rLagrangianModelPart, BinBasedFastPointLocator<TDim>& node_locator)
    {
        KRATOS_TRY

        array_1d<double, 3 > veulerian;
        array_1d<double, TDim + 1 > N;
        const int max_results = 10000;
        typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);


        const int nparticles = rLagrangianModelPart.Nodes().size();
//KRATOS_WATCH("551")
        #pragma omp parallel for firstprivate(results,N,veulerian)
        for (int i = 0; i < nparticles; i++)
        {
            unsigned int substep = 0;
            unsigned int subdivisions = 1;
            double small_dt = dt;

            while(substep++ < subdivisions)
            {
                ModelPart::NodesContainerType::iterator iparticle = rLagrangianModelPart.NodesBegin() + i;
                (iparticle)->Set(TO_ERASE, true);

                Node < 3 > ::Pointer pparticle = *(iparticle.base());
                typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();

                Element::Pointer pelement;
                //      KRATOS_WATCH("561")

                bool is_found = node_locator.FindPointOnMesh(pparticle->Coordinates(), N, pelement, result_begin, max_results);
                //        KRATOS_WATCH("564")

                if (is_found == true)
                {
                    (pparticle)->GetValue(IS_VISITED) = 1;

                    Geometry< Node < 3 > >& geom = pelement->GetGeometry();

                    noalias(veulerian) = N[0] * geom[0].FastGetSolutionStepValue(VELOCITY);
                    for (unsigned int k = 1; k < geom.size(); k++)
                        noalias(veulerian) += N[k] * geom[k].FastGetSolutionStepValue(VELOCITY);

                    //compute adaptive subdivisions
                    if(substep == 1)
                    {
                        //compute h
                        double h = N[0] * geom[0].FastGetSolutionStepValue(NODAL_H);
                        for (unsigned int k = 1; k < geom.size(); k++)
                            h += N[k] * geom[k].FastGetSolutionStepValue(NODAL_H);

                        //compute number of subdivisions needed
                        const unsigned int min_subdivisions = 3;
                        const unsigned int max_subdivisions = 20;
                        double v = norm_2(veulerian);
                        subdivisions = double(floor(2*dt*v/h));
                        subdivisions = (subdivisions<min_subdivisions) ? min_subdivisions : (subdivisions>max_subdivisions) ? max_subdivisions : subdivisions;

                        //compute subdivisions time step
                        small_dt = dt / subdivisions;

//                         KRATOS_WATCH(subdivisions)

                    }

                    //move according to the streamline

                    array_1d<double, 3 > & disp = (iparticle)->FastGetSolutionStepValue(DISPLACEMENT);

                    noalias(disp) += small_dt*veulerian;

                    (pparticle)->Set(TO_ERASE, false);
                    //	KRATOS_WATCH("585")


                    //update position
                    noalias(iparticle->Coordinates()) = iparticle->GetInitialPosition();
                    noalias(iparticle->Coordinates()) += iparticle->FastGetSolutionStepValue(DISPLACEMENT);
                    (iparticle)->GetValue(IS_VISITED) = 0;
//KRATOS_WATCH("619")
                }
            }
        }

        KRATOS_CATCH("")
    }
    //**********************************************************************************************
    //**********************************************************************************************

    void ParticleLevelSetCorrection(ModelPart& rEulerianModelPart, ModelPart& rLagrangianModelPart, BinBasedFastPointLocator<TDim>& node_locator)
    {
        KRATOS_TRY
        //Initilize NAGATIVE_DISTANCE & POSETIVE_DISTANCE
        const int nnodes= rEulerianModelPart.Nodes().size();
        #pragma omp parallel for
        for (int jj = 0; jj < nnodes; jj++)
        {
          ModelPart::NodesContainerType::iterator node_itr = rEulerianModelPart.NodesBegin() + jj;

	  const double  nd_dist = node_itr->FastGetSolutionStepValue(DISTANCE);

	  node_itr->SetValue(POSETIVE_DISTANCE,nd_dist );
	  node_itr->SetValue(NAGATIVE_DISTANCE,nd_dist );

	}
        //loop over particles
	double particle_dist= 0.0;
        array_1d<double, TDim + 1 > N;
        const int max_results = 10000;
        typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);


        const int nparticles = rLagrangianModelPart.Nodes().size();
        #pragma omp parallel for firstprivate(results,N,particle_dist)
        for (int i = 0; i < nparticles; i++)
        {
          ModelPart::NodesContainerType::iterator particle_itr = rLagrangianModelPart.NodesBegin() + i;

	  Node < 3 > ::Pointer p_pointer = *(particle_itr.base());
	  typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();

	  Element::Pointer pelement;
	  bool is_found = node_locator.FindPointOnMesh(p_pointer->Coordinates(), N, pelement, result_begin, max_results);
	  if (is_found == true)
	  {
	    Geometry< Node < 3 > >& geom = pelement->GetGeometry();
            //interpolate the particle distance
	    particle_dist = N[0] * geom[0].FastGetSolutionStepValue(DISTANCE);
	    for (unsigned int k = 1; k < geom.size(); k++)
		particle_dist += N[k] * geom[k].FastGetSolutionStepValue(DISTANCE);

	    //check if correction is needed
	    const double p_sign = particle_itr->FastGetSolutionStepValue(IS_WATER);
	    const double p_radi = particle_itr->FastGetSolutionStepValue(PARTICLE_RADIUS);
	    if( particle_dist*p_sign < 0.0 && fabs(particle_dist) > p_radi)
	    {
	      double p_xx = particle_itr->X();
	      double p_yy = particle_itr->Y();
	      double p_zz = particle_itr->Z();
// 	      const Variable<double> posetive_negative_dist_var;

/*	      if( p_sign == -1.0 )
		posetive_negative_dist_var = NAGATIVE_DISTANCE;
	      else if( p_sign == 1.0 )
		posetive_negative_dist_var = POSETIVE_DISTANCE;	*/

	      for (unsigned int kk = 1; kk < geom.size(); kk++){
		p_xx -= geom[kk].X();
		p_yy -= geom[kk].Y();
		p_zz -= geom[kk].Z();

		double dd = p_xx*p_xx + p_yy*p_yy + p_zz*p_zz;
		dd = sqrt(dd);

		double dist_to_particle = p_sign * (p_radi - dd);

		 //correction due to particle distance and sign
		geom[kk].SetLock();
		if( p_sign == 1.0){
		  double& pos_distance = geom[kk].GetValue(POSETIVE_DISTANCE);
		  if ( dist_to_particle > pos_distance)
		       pos_distance = dist_to_particle;}
		else if( p_sign == -1.0){
		  double& neg_distance = geom[kk].GetValue(NAGATIVE_DISTANCE);
		  if ( dist_to_particle < neg_distance)
		       neg_distance = dist_to_particle;	}

		geom[kk].UnSetLock();

	      }
	    }
	  }
	}//end of loop over particles
	//final correction, choose between NAGATIVE_DISTANCE & POSETIVE_DISTANCE
//         const int nnodes= rEulerianModelPart.Nodes().size();
        #pragma omp parallel for
        for (int jj = 0; jj < nnodes; jj++)
        {
          ModelPart::NodesContainerType::iterator node_itr = rEulerianModelPart.NodesBegin() + jj;

	  double posetive = node_itr->GetValue(POSETIVE_DISTANCE);
	  double negative = node_itr->GetValue(NAGATIVE_DISTANCE);
	  double & nd_dist = node_itr->FastGetSolutionStepValue(DISTANCE);

	  if ( posetive != negative){
	    if( fabs(posetive) < fabs(negative) )
	      nd_dist = posetive;
	    else
	      nd_dist = negative;

	    node_itr->SetValue(POSETIVE_DISTANCE,nd_dist );
	    node_itr->SetValue(NAGATIVE_DISTANCE,nd_dist );
	  }
	}
        KRATOS_CATCH("")
    }
    //**********************************************************************************************
    //**********************************************************************************************
  void ResetParticleRadius(const double min_edge_length, ModelPart& rEulerianModelPart, ModelPart& rLagrangianModelPart, BinBasedFastPointLocator<TDim>& node_locator)
  {
    KRATOS_TRY;
    double particle_dist = 0.0;
    array_1d<double, TDim + 1 > N;
    const int max_results = 10000;
    typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);


    const int nparticles = rLagrangianModelPart.Nodes().size();
    #pragma omp parallel for firstprivate(results,N,particle_dist)
    for (int i = 0; i < nparticles; i++)
    {
	ModelPart::NodesContainerType::iterator particle_itr = rLagrangianModelPart.NodesBegin() + i;

	Node < 3 > ::Pointer p_pointer = *(particle_itr.base());
	typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();

	Element::Pointer pelement;
	bool is_found = node_locator.FindPointOnMesh(p_pointer->Coordinates(), N, pelement, result_begin, max_results);
	if (is_found == true)
	{
	  Geometry< Node < 3 > >& geom = pelement->GetGeometry();
	  //interpolate the particle distance
	  particle_dist = N[0] * geom[0].FastGetSolutionStepValue(DISTANCE);
	  for (unsigned int k = 1; k < geom.size(); k++)
	      particle_dist += N[k] * geom[k].FastGetSolutionStepValue(DISTANCE);

	  if( fabs(particle_dist) < 0.1*min_edge_length)
	    particle_itr->FastGetSolutionStepValue(PARTICLE_RADIUS) = 0.1*min_edge_length;
	  else if(fabs(particle_dist) > 0.5*min_edge_length)
	    particle_itr->FastGetSolutionStepValue(PARTICLE_RADIUS) = 0.5*min_edge_length;
	  else
	    particle_itr->FastGetSolutionStepValue(PARTICLE_RADIUS) = fabs(particle_dist);
	}
    }
    KRATOS_CATCH("")
  }
    //**********************************************************************************************
    //**********************************************************************************************
    void ParticleReseeding(ModelPart& rEulerianModelPart, ModelPart& rLagrangianModelPart, BinBasedFastPointLocator<TDim>& node_locator, const double max_seed_distance, const double min_edge_size)
    {
        KRATOS_TRY;
        //generate a tree with the position of the lagrangian nodes
//         typedef Node < 3 > PointType;
//         typedef Node < 3 > ::Pointer PointTypePointer;

        //unsigned int min_number_of_particles = 1;

        for (ModelPart::ElementsContainerType::iterator el_it = rEulerianModelPart.ElementsBegin();
                el_it != rEulerianModelPart.ElementsEnd(); el_it++)
        {
            el_it->SetValue(YOUNG_MODULUS,0.0);
        }
        for (ModelPart::NodesContainerType::iterator pparticle = rLagrangianModelPart.NodesBegin();
                pparticle != rLagrangianModelPart.NodesEnd(); pparticle++)
        {
            pparticle->Set(TO_ERASE,false);
            pparticle->SetValue(NL_ITERATION_NUMBER,(rEulerianModelPart.ElementsBegin())->Id());
            pparticle->SetValue(IS_ESCAPED,false);
            pparticle->SetValue(IS_VISITED,0);
        }

        //count particles that fall within an element
        array_1d<double, TDim + 1 > N;
        const int max_results = 10000;
        typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);
        const int nparticles = rLagrangianModelPart.Nodes().size();

        //count particles within an element
        #pragma omp parallel for firstprivate(results,N)
        for (int i = 0; i < nparticles; i++)
        {
            ModelPart::NodesContainerType::iterator iparticle = rLagrangianModelPart.NodesBegin() + i;

            Node < 3 > ::Pointer pparticle = *(iparticle.base());
            typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();

            Element::Pointer pelement;

            bool is_found = node_locator.FindPointOnMesh(pparticle->Coordinates(), N, pelement, result_begin, max_results);

            if (is_found == true)
            {
		const double particle_sign = iparticle->FastGetSolutionStepValue(IS_WATER);
	        Geometry< Node < 3 > >& geom = pelement->GetGeometry();
	        bool is_scaped = CheckIfEscaped(geom,N,particle_sign);
		iparticle->SetValue(IS_ESCAPED,is_scaped);

		if( CheckElemDist(geom,max_seed_distance) )// if it is inside the 3h band
		{
		  double& counter = pelement->GetValue(YOUNG_MODULUS);
		  #pragma omp atomic
		  counter += 1.0;

		  iparticle->SetValue(NL_ITERATION_NUMBER , pelement->Id());
		}
		else
		{
		  if( is_scaped == false) //delete if it is not an escaped particle
		    iparticle->Set(TO_ERASE,true);
		}
            }
        }

        //loop over close to the surface elements to ressed or delet particles
	if(TDim==2){
	    ReseedOrDelete2D(rEulerianModelPart, rLagrangianModelPart, max_seed_distance, min_edge_size);
	}
	else
	{
	    const int max_num_ptr = 16;//70;
	    const int num_ptr = 10;//56;
	    const int min_num_ptr = 6;//40;
	    MarkEraseExtraParticles3D(rEulerianModelPart, rLagrangianModelPart, max_seed_distance, min_edge_size, max_num_ptr, num_ptr);
	    ReseedPoorElements3D(rEulerianModelPart, rLagrangianModelPart, max_seed_distance, min_edge_size, min_num_ptr, num_ptr );

	    FewReseedPoorElements3D(rEulerianModelPart, rLagrangianModelPart, max_seed_distance, min_edge_size, min_num_ptr, num_ptr );
	}

        //perform the erase
        NodeEraseProcess(rLagrangianModelPart).Execute();

      KRATOS_CATCH("");
    }
    //**********************************************************************************************
    //**********************************************************************************************

    void VisualizationModelPart(ModelPart& rCompleteModelPart, ModelPart& rEulerianModelPart, ModelPart & rLagrangianModelPart)
    {
        KRATOS_TRY;

        rCompleteModelPart.Elements() = rEulerianModelPart.Elements();
        rCompleteModelPart.Nodes() = rEulerianModelPart.Nodes();


        unsigned int id;
        if(rEulerianModelPart.Nodes().size()!= 0)
            id = (rEulerianModelPart.Nodes().end() - 1)->Id() + 1;
        else
            id = 1;

        //preallocate the memory needed
        int tot_nodes = rEulerianModelPart.Nodes().size() + rLagrangianModelPart.Nodes().size();
        rCompleteModelPart.Nodes().reserve( tot_nodes );

        //note that here we renumber the nodes
        for (ModelPart::NodesContainerType::iterator node_it = rLagrangianModelPart.NodesBegin();
                node_it != rLagrangianModelPart.NodesEnd(); node_it++)
        {
            rCompleteModelPart.AddNode(*(node_it.base()));
            node_it->SetId(id++);
        }

        KRATOS_CATCH("");
    }

     //**********************************************************************************
    //**********************************************************************************
    void FindMaxMinEdgeSize(ModelPart& r_model_part, pybind11::list& maxmin)
    {
        KRATOS_TRY

        double max_edge = 0.0;
        double min_edge = 1000.0;

        for(ModelPart::ElementsContainerType::iterator it=r_model_part.ElementsBegin(); it!=r_model_part.ElementsEnd(); it++)
        {
            Geometry<Node<3> >&geom = it->GetGeometry();

            double loc_h_max = 0.0;
            double loc_h_min = 1000.0;

            for(unsigned int i=0; i<TDim+1; i++)
            {

                double xc = geom[i].X();
                double yc = geom[i].Y();
                double zc = geom[i].Z();
                for(unsigned int j=i+1; j<TDim+1; j++)
                {
                    double x = geom[j].X();
                    double y = geom[j].Y();
                    double z = geom[j].Z();
                    double l = (x - xc)*(x - xc);
                    l += (y - yc)*(y - yc);
                    l += (z - zc)*(z - zc);

                    if (l > loc_h_max) loc_h_max = l;
                    else if(l < loc_h_min) loc_h_min = l;
                }
            }

            loc_h_max = sqrt(loc_h_max);
            loc_h_min = sqrt(loc_h_min);

            if(loc_h_max > max_edge )  max_edge = loc_h_max;
            if(loc_h_min < min_edge )  min_edge = loc_h_min;

        }

//         r_model_part.GetCommunicator().MaxAll(h_max);


        maxmin.append(max_edge);
        maxmin.append(min_edge);

        KRATOS_CATCH("");
    }

private:

  void CreateParticles3D(ModelPart& rEulerianModelPart,
			 ModelPart& rLagrangianModelPart,
			 BoundedMatrix<double, 56, 3 > pos,
			 BoundedMatrix<double, 56, 4 > N,
			 const double max_seed_distance,
			 const double min_edge_size)
  {

         unsigned int id = (rEulerianModelPart.Nodes().end() - 1)->Id() + 1;

 	for (ModelPart::ElementsContainerType::iterator el_it = rEulerianModelPart.ElementsBegin();
                el_it != rEulerianModelPart.ElementsEnd(); el_it++)
        {
	 Geometry<Node < 3 > >& geom = el_it->GetGeometry();

	 if(CheckElemDist(geom,max_seed_distance))
	 {
            ComputeGaussPointPositions3D(geom, pos, N);
            for (unsigned int i = 0; i < pos.size1(); i++)
            {
                int node_id = id++;
                Node < 3 > ::Pointer pnode = rLagrangianModelPart.CreateNewNode(node_id, pos(i, 0), pos(i, 1), pos(i, 2));

                array_1d<double, 3 > & vel = pnode->FastGetSolutionStepValue(VELOCITY);
                noalias(vel) = ZeroVector(3);

// 		double & p_dist = pnode->FastGetSolutionStepValue(DISTANCE);
                double p_distance = 0.0;
                for (unsigned int j = 0; j < TDim + 1; j++){
                    noalias(vel) += N(i, j) * geom[j].FastGetSolutionStepValue(VELOCITY);
		    p_distance += N(i, j) * geom[j].FastGetSolutionStepValue(DISTANCE);
		}

		// Assign particle sign
		if(p_distance < 0.0)
		  pnode->FastGetSolutionStepValue(IS_WATER)=-1.0;
		else if(p_distance > 0.0)
		  pnode->FastGetSolutionStepValue(IS_WATER)= 1.0;

		pnode->Fix(IS_WATER);

		AssignParticleRadius(pnode,p_distance,min_edge_size);

            }
	  }
        }

  }

   void FewCreateParticles3D(ModelPart& rEulerianModelPart,
			 ModelPart& rLagrangianModelPart,
			 BoundedMatrix<double, 10, 3 > pos,
			 BoundedMatrix<double, 10, 4 > N,
			 const double max_seed_distance,
			 const double min_edge_size)
  {

         unsigned int id = (rEulerianModelPart.Nodes().end() - 1)->Id() + 1;

 	for (ModelPart::ElementsContainerType::iterator el_it = rEulerianModelPart.ElementsBegin();
                el_it != rEulerianModelPart.ElementsEnd(); el_it++)
        {
	 Geometry<Node < 3 > >& geom = el_it->GetGeometry();

	 if(CheckElemDist(geom,max_seed_distance))
	 {
            FewComputeGaussPointPositions3D(geom, pos, N);
            for (unsigned int i = 0; i < pos.size1(); i++)
            {
                int node_id = id++;
                Node < 3 > ::Pointer pnode = rLagrangianModelPart.CreateNewNode(node_id, pos(i, 0), pos(i, 1), pos(i, 2));

                array_1d<double, 3 > & vel = pnode->FastGetSolutionStepValue(VELOCITY);
                noalias(vel) = ZeroVector(3);

// 		double & p_dist = pnode->FastGetSolutionStepValue(DISTANCE);
                double p_distance = 0.0;
                for (unsigned int j = 0; j < TDim + 1; j++){
                    noalias(vel) += N(i, j) * geom[j].FastGetSolutionStepValue(VELOCITY);
		    p_distance += N(i, j) * geom[j].FastGetSolutionStepValue(DISTANCE);
		}

		// Assign particle sign
		if(p_distance < 0.0)
		  pnode->FastGetSolutionStepValue(IS_WATER)=-1.0;
		else if(p_distance > 0.0)
		  pnode->FastGetSolutionStepValue(IS_WATER)= 1.0;

		pnode->Fix(IS_WATER);

		AssignParticleRadius(pnode,p_distance,min_edge_size);

            }
	  }
        }

  }

   void CreateParticles2D(ModelPart& rEulerianModelPart,
			 ModelPart& rLagrangianModelPart,
			 BoundedMatrix<double, 16, 3 > pos,
			 BoundedMatrix<double, 16, 3 > N,
			 const double max_seed_distance,
			 const double min_edge_size)
  {
         unsigned int id = (rEulerianModelPart.Nodes().end() - 1)->Id() + 1;

 	for (ModelPart::ElementsContainerType::iterator el_it = rEulerianModelPart.ElementsBegin();
                el_it != rEulerianModelPart.ElementsEnd(); el_it++)
        {
	 Geometry<Node < 3 > >& geom = el_it->GetGeometry();

	 if(CheckElemDist(geom,max_seed_distance))
	 {
            ComputeGaussPointPositions2D(geom, pos, N);
            for (unsigned int i = 0; i < pos.size1(); i++)
            {
                int node_id = id++;
                Node < 3 > ::Pointer pnode = rLagrangianModelPart.CreateNewNode(node_id, pos(i, 0), pos(i, 1), pos(i, 2));

                array_1d<double, 3 > & vel = pnode->FastGetSolutionStepValue(VELOCITY);
                noalias(vel) = ZeroVector(3);

// 		double & p_dist = pnode->FastGetSolutionStepValue(DISTANCE);
                double p_distance = 0.0;
                for (unsigned int j = 0; j < TDim + 1; j++){
                    noalias(vel) += N(i, j) * geom[j].FastGetSolutionStepValue(VELOCITY);
		    p_distance += N(i, j) * geom[j].FastGetSolutionStepValue(DISTANCE);
		}

		// Assign particle sign
		if(p_distance < 0.0)
		  pnode->FastGetSolutionStepValue(IS_WATER)=-1.0;
		else if(p_distance > 0.0)
		  pnode->FastGetSolutionStepValue(IS_WATER)= 1.0;

		pnode->Fix(IS_WATER);

		AssignParticleRadius(pnode,p_distance,min_edge_size);

            }
	  }
        }

  }
    void ReseedOrDelete2D(ModelPart& rEulerianModelPart,
			  ModelPart& rLagrangianModelPart,
			  const double max_seed_distance,
			  const double min_edge_size)
    {
        int id;
        if (rLagrangianModelPart.Nodes().size() != 0)
            id = (rLagrangianModelPart.NodesEnd() - 1)->Id();
        else
            id = 1;
        const int nelements = rEulerianModelPart.Elements().size();
	const int nparticles = rLagrangianModelPart.Nodes().size();



        BoundedMatrix<double, 16, 3 > coord;
        BoundedMatrix<double, 16, 3 > NN;
//         #pragma omp parallel for firstprivate(NN,coord)
        for (int ne = 0; ne < nelements; ne++)
        {
          ModelPart::ElementsContainerType::iterator ielem = rEulerianModelPart.ElementsBegin() + ne;
	  Geometry<Node < 3 > >& geom = ielem->GetGeometry();
	  int n_ptr = int(ielem->GetValue(YOUNG_MODULUS));

	  if( n_ptr < 12 && CheckElemDist(geom,max_seed_distance) )//ressed in close to surface and poor element
	  {
	    //compute cooordinates
	    //RandomPariclePosition(geom, coord, NN);
	    ComputeGaussPointPositions2D(geom, coord, NN);
	    int aux_n_ptr = n_ptr;
	    int cnt = 0;
	   while( aux_n_ptr<16 ){
	        aux_n_ptr++;
		//COORDINATES
                int node_id = id++;
                Node < 3 > ::Pointer pnode = rLagrangianModelPart.CreateNewNode(node_id, coord(cnt,0), coord(cnt,1), coord(cnt,2));

                array_1d<double, 3 > & vel = pnode->FastGetSolutionStepValue(VELOCITY);
                noalias(vel) = ZeroVector(3);

// 		double & p_dist = pnode->FastGetSolutionStepValue(DISTANCE);
                double p_distance = 0.0;
                for (unsigned int j = 0; j < TDim + 1; j++){
                    noalias(vel) += NN(cnt,j) * geom[j].FastGetSolutionStepValue(VELOCITY);
		    p_distance += NN(cnt,j) * geom[j].FastGetSolutionStepValue(DISTANCE);
		}

		// Assign particle sign
		if(p_distance < 0.0)
		  pnode->FastGetSolutionStepValue(IS_WATER)=-1.0;
		else if(p_distance > 0.0)
		  pnode->FastGetSolutionStepValue(IS_WATER)= 1.0;

		pnode->Fix(IS_WATER);

		AssignParticleRadius(pnode,p_distance,min_edge_size);

		cnt++;
	    }
	  }
	  else if( n_ptr > 20 && CheckElemDist(geom,max_seed_distance) ){
	    const int ele_id = ielem->Id();
	    ModelPart::NodesContainerType element_particles;
            element_particles.reserve(64);
	    //save particle list
	    for (int kk = 0; kk < nparticles; kk++)
	    {
		ModelPart::NodesContainerType::iterator iparticle = rLagrangianModelPart.NodesBegin() + kk;

		const int ptr_nest = iparticle->GetValue(NL_ITERATION_NUMBER);
		if( ptr_nest==ele_id )
		{
		  iparticle->SetValue(SCALE, 0);
		  element_particles.push_back( *(iparticle.base()) );
		}
	    }

	    //loop to order based on the radius
            ModelPart::NodesContainerType::iterator ptr_begin = element_particles.begin();
	    unsigned int ptr_elem_size = element_particles.size();

	    for(unsigned int ii=0; ii < ptr_elem_size; ii++)
	      for(unsigned int jj=ii+1; jj < ptr_elem_size; jj++)
	      {
		double ii_radi = (ptr_begin + ii)->FastGetSolutionStepValue(PARTICLE_RADIUS);
		double jj_radi = (ptr_begin + jj)->FastGetSolutionStepValue(PARTICLE_RADIUS);

		(ii_radi>=jj_radi) ? (ptr_begin + ii)->GetValue(SCALE)+=1 : (ptr_begin + jj)->GetValue(SCALE)+=1;

	      }
	      //delete extra nodes
	       int aux_ptr_elem_size = int(ptr_elem_size);
	      while(aux_ptr_elem_size>16)
	      {
	        for(unsigned int ii=0; ii < ptr_elem_size; ii++){
		  bool swt = false;
		  for( int  kkk = ptr_elem_size; kkk>0; kkk-- )
		     if( (ptr_begin + ii)->GetValue(SCALE) == kkk && (ptr_begin + ii)->GetValue(IS_VISITED) == 0){
		        bool is_escaped = (ptr_begin + ii)->GetValue(IS_ESCAPED);
			if( is_escaped==false )
		            (ptr_begin + ii)->Set(TO_ERASE,true);//CHECK ESCASPED NODES
			(ptr_begin + ii)->SetValue(IS_VISITED,1);
			swt = true;
			break;
		     }
		     if(swt )
		       break;
		}
		aux_ptr_elem_size -= 1;
	      }
	  }


	}

    }

     void MarkEraseExtraParticles3D(ModelPart& rEulerianModelPart,
			  ModelPart& rLagrangianModelPart,
			  const double max_seed_distance,
			  const double min_edge_size,
		          const int max_num_particle,
			  const int num_particle)
    {
//        int id;
//        if (rLagrangianModelPart.Nodes().size() != 0)
//            id = (rLagrangianModelPart.NodesEnd() - 1)->Id();
//        else
//            id = 1;
        const int nelements = rEulerianModelPart.Elements().size();
	const int nparticles = rLagrangianModelPart.Nodes().size();

	std::vector< WeakPointerVector< Node< 3> > > particle_of_element(nelements);
// 	particle_of_element.reserve(nelements);

	std::vector< unsigned int > num_ptr_in_elem(nelements,0);
// 	num_ptr_in_elem.reserve(nelements);

	//loop on elements to resrve the size of particle in element list
        #pragma omp parallel for firstprivate(num_ptr_in_elem)
        for (int ne = 0; ne < nelements; ne++)
        {
          ModelPart::ElementsContainerType::iterator ielem = rEulerianModelPart.ElementsBegin() + ne;
	  int n_ptr = int(ielem->GetValue(YOUNG_MODULUS));
	  unsigned int ele_id = ielem->Id();

	  num_ptr_in_elem[ele_id-1] = n_ptr;

	  if(n_ptr > max_num_particle)
	    particle_of_element[ele_id-1].reserve(n_ptr);
	}

       //loop on particles to push_back particle related to full elements
	for (int kk = 0; kk < nparticles; kk++)
	{
	    ModelPart::NodesContainerType::iterator iparticle = rLagrangianModelPart.NodesBegin() + kk;

	    const int ptr_nest = iparticle->GetValue(NL_ITERATION_NUMBER);
        if( num_ptr_in_elem[ptr_nest-1] > static_cast<unsigned int>(max_num_particle) )
	      particle_of_element[ptr_nest-1].push_back( *(iparticle.base()) );
	}

	//loop over elements to reoreder the particle radius in over populated elements
        #pragma omp parallel for firstprivate(particle_of_element)
    for( int ii = 0; ii< static_cast<int>(particle_of_element.size()); ++ii)
	{
      if(particle_of_element[ii].size() > static_cast<unsigned int>(max_num_particle))
	  {
	    //sort
	    std::sort(particle_of_element[ii].ptr_begin(), particle_of_element[ii].ptr_end(), RadiusCompare() );

	    //delete extra nodes
	    WeakPointerVector< Node< 3> >::iterator ele_pt_ptr = particle_of_element[ii].begin();
	    const unsigned int this_ele_ptr = particle_of_element[ii].size();
	    int aux_ptr_elem_size = this_ele_ptr;

	      for( unsigned int ij = 0; (ij < this_ele_ptr && aux_ptr_elem_size > num_particle); ++ij)
	      {
		bool is_escaped = (ele_pt_ptr + ij)->GetValue(IS_ESCAPED);
		if( is_escaped==false ){
		   (ele_pt_ptr + ij)->Set(TO_ERASE,true);
		   aux_ptr_elem_size--;
		}
	      }

	   }

	 }

    }

struct RadiusCompare{
  template<class TRefrenceType>
bool operator()(const TRefrenceType  ptr_a, const TRefrenceType ptr_b)
    {
      double a_radi = ptr_a.lock()->FastGetSolutionStepValue(PARTICLE_RADIUS);
      double b_radi = ptr_b.lock()->FastGetSolutionStepValue(PARTICLE_RADIUS);
	return (a_radi > b_radi);
    }
};
     void ReseedPoorElements3D(ModelPart& rEulerianModelPart,
			  ModelPart& rLagrangianModelPart,
			  const double max_seed_distance,
			  const double min_edge_size,
			  const int min_num_particle,
			  const int num_particle)
    {
        int id;
        if (rLagrangianModelPart.Nodes().size() != 0)
            id = (rLagrangianModelPart.NodesEnd() - 1)->Id();
        else
            id = 1;
        const int nelements = rEulerianModelPart.Elements().size();
//	const int nparticles = rLagrangianModelPart.Nodes().size();



        BoundedMatrix<double, 56, 3 > coord;
        BoundedMatrix<double, 56, 4 > NN;
//         #pragma omp parallel for firstprivate(NN,coord)
        for (int ne = 0; ne < nelements; ne++)
        {
          ModelPart::ElementsContainerType::iterator ielem = rEulerianModelPart.ElementsBegin() + ne;
	  Geometry<Node < 3 > >& geom = ielem->GetGeometry();
	  int n_ptr = int(ielem->GetValue(YOUNG_MODULUS));

	  if( n_ptr < min_num_particle && CheckElemDist(geom,max_seed_distance) )//ressed in close to surface and poor element
	  {
	    //compute cooordinates
	    //RandomPariclePosition(geom, coord, NN);
	    ComputeGaussPointPositions3D(geom, coord, NN);
	    int aux_n_ptr = n_ptr;
	    int cnt = 0;
	   while( aux_n_ptr < num_particle ){
	        aux_n_ptr++;
		//COORDINATES
                int node_id = id++;
                Node < 3 > ::Pointer pnode = rLagrangianModelPart.CreateNewNode(node_id, coord(cnt,0), coord(cnt,1), coord(cnt,2));

                array_1d<double, 3 > & vel = pnode->FastGetSolutionStepValue(VELOCITY);
                noalias(vel) = ZeroVector(3);

// 		double & p_dist = pnode->FastGetSolutionStepValue(DISTANCE);
                double p_distance = 0.0;
                for (unsigned int j = 0; j < TDim + 1; j++){
                    noalias(vel) += NN(cnt,j) * geom[j].FastGetSolutionStepValue(VELOCITY);
		    p_distance += NN(cnt,j) * geom[j].FastGetSolutionStepValue(DISTANCE);
		}

		// Assign particle sign
		if(p_distance < 0.0)
		  pnode->FastGetSolutionStepValue(IS_WATER)=-1.0;
		else if(p_distance > 0.0)
		  pnode->FastGetSolutionStepValue(IS_WATER)= 1.0;

		pnode->Fix(IS_WATER);

		AssignParticleRadius(pnode,p_distance,min_edge_size);

		cnt++;
	    }
	  }
	}

    }

      void FewReseedPoorElements3D(ModelPart& rEulerianModelPart,
			  ModelPart& rLagrangianModelPart,
			  const double max_seed_distance,
			  const double min_edge_size,
			  const int min_num_particle,
			  const int num_particle)
    {
        int id;
        if (rLagrangianModelPart.Nodes().size() != 0)
            id = (rLagrangianModelPart.NodesEnd() - 1)->Id();
        else
            id = 1;
        const int nelements = rEulerianModelPart.Elements().size();
//	const int nparticles = rLagrangianModelPart.Nodes().size();



        BoundedMatrix<double, 10, 3 > coord;
        BoundedMatrix<double, 10, 4 > NN;
//         #pragma omp parallel for firstprivate(NN,coord)
        for (int ne = 0; ne < nelements; ne++)
        {
          ModelPart::ElementsContainerType::iterator ielem = rEulerianModelPart.ElementsBegin() + ne;
	  Geometry<Node < 3 > >& geom = ielem->GetGeometry();
	  int n_ptr = int(ielem->GetValue(YOUNG_MODULUS));

	  if( n_ptr < min_num_particle && CheckElemDist(geom,max_seed_distance) )//ressed in close to surface and poor element
	  {
	    //compute cooordinates
	    //RandomPariclePosition(geom, coord, NN);
	    FewComputeGaussPointPositions3D(geom, coord, NN);
	    int aux_n_ptr = n_ptr;
	    int cnt = 0;
	   while( aux_n_ptr < num_particle ){
	        aux_n_ptr++;
		//COORDINATES
                int node_id = id++;
                Node < 3 > ::Pointer pnode = rLagrangianModelPart.CreateNewNode(node_id, coord(cnt,0), coord(cnt,1), coord(cnt,2));

                array_1d<double, 3 > & vel = pnode->FastGetSolutionStepValue(VELOCITY);
                noalias(vel) = ZeroVector(3);

// 		double & p_dist = pnode->FastGetSolutionStepValue(DISTANCE);
                double p_distance = 0.0;
                for (unsigned int j = 0; j < TDim + 1; j++){
                    noalias(vel) += NN(cnt,j) * geom[j].FastGetSolutionStepValue(VELOCITY);
		    p_distance += NN(cnt,j) * geom[j].FastGetSolutionStepValue(DISTANCE);
		}

		// Assign particle sign
		if(p_distance < 0.0)
		  pnode->FastGetSolutionStepValue(IS_WATER)=-1.0;
		else if(p_distance > 0.0)
		  pnode->FastGetSolutionStepValue(IS_WATER)= 1.0;

		pnode->Fix(IS_WATER);

		AssignParticleRadius(pnode,p_distance,min_edge_size);

		cnt++;
	    }
	  }
	}

    }
//     void ReseedOrDelete3D(ModelPart& rEulerianModelPart,
// 			  ModelPart& rLagrangianModelPart,
// 			  const double max_seed_distance,
// 			  const double min_edge_size)
//     {
//         int id;
//         if (rLagrangianModelPart.Nodes().size() != 0)
//             id = (rLagrangianModelPart.NodesEnd() - 1)->Id();
//         else
//             id = 1;
//         const int nelements = rEulerianModelPart.Elements().size();
// 	const int nparticles = rLagrangianModelPart.Nodes().size();
//
//
//
//         BoundedMatrix<double, 56, 3 > coord;
//         BoundedMatrix<double, 56, 4 > NN;
// //         #pragma omp parallel for firstprivate(NN,coord)
//         for (int ne = 0; ne < nelements; ne++)
//         {
//           ModelPart::ElementsContainerType::iterator ielem = rEulerianModelPart.ElementsBegin() + ne;
// 	  Geometry<Node < 3 > >& geom = ielem->GetGeometry();
// 	  int n_ptr = int(ielem->GetValue(YOUNG_MODULUS));
//
// 	  if( n_ptr < 42 && CheckElemDist(geom,max_seed_distance) )//ressed in close to surface and poor element
// 	  {
// 	    //compute cooordinates
// 	    //RandomPariclePosition(geom, coord, NN);
// 	    ComputeGaussPointPositions3D(geom, coord, NN);
// 	    int aux_n_ptr = n_ptr;
// 	    int cnt = 0;
// 	   while( aux_n_ptr<56 ){
// 	        aux_n_ptr++;
// 		//COORDINATES
//                 int node_id = id++;
//                 Node < 3 > ::Pointer pnode = rLagrangianModelPart.CreateNewNode(node_id, coord(cnt,0), coord(cnt,1), coord(cnt,2));
//
//                 array_1d<double, 3 > & vel = pnode->FastGetSolutionStepValue(VELOCITY);
//                 noalias(vel) = ZeroVector(3);
//
// // 		double & p_dist = pnode->FastGetSolutionStepValue(DISTANCE);
//                 double p_distance = 0.0;
//                 for (unsigned int j = 0; j < TDim + 1; j++){
//                     noalias(vel) += NN(cnt,j) * geom[j].FastGetSolutionStepValue(VELOCITY);
// 		    p_distance += NN(cnt,j) * geom[j].FastGetSolutionStepValue(DISTANCE);
// 		}
//
// 		// Assign particle sign
// 		if(p_distance < 0.0)
// 		  pnode->FastGetSolutionStepValue(IS_WATER)=-1.0;
// 		else if(p_distance > 0.0)
// 		  pnode->FastGetSolutionStepValue(IS_WATER)= 1.0;
//
// 		pnode->Fix(IS_WATER);
//
// 		AssignParticleRadius(pnode,p_distance,min_edge_size);
//
// 		cnt++;
// 	    }
// 	  }
// 	  else if( n_ptr > 70 && CheckElemDist(geom,max_seed_distance) ){
// 	    const int ele_id = ielem->Id();
// 	    ModelPart::NodesContainerType element_particles;
//             element_particles.reserve(64);
// 	    //save particle list
// 	    for (int kk = 0; kk < nparticles; kk++)
// 	    {
// 		ModelPart::NodesContainerType::iterator iparticle = rLagrangianModelPart.NodesBegin() + kk;
//
// 		const int ptr_nest = iparticle->GetValue(NL_ITERATION_NUMBER);
// 		if( ptr_nest==ele_id )
// 		{
// 		  iparticle->SetValue(SCALE, 0);
// 		  element_particles.push_back( *(iparticle.base()) );
// 		}
// 	    }
//
// 	    //loop to order based on the radius
//             ModelPart::NodesContainerType::iterator ptr_begin = element_particles.begin();
// 	    unsigned int ptr_elem_size = element_particles.size();
//
// 	    for(unsigned int ii=0; ii < ptr_elem_size; ii++)
// 	      for(unsigned int jj=ii+1; jj < ptr_elem_size; jj++)
// 	      {
// 		double ii_radi = (ptr_begin + ii)->FastGetSolutionStepValue(PARTICLE_RADIUS);
// 		double jj_radi = (ptr_begin + jj)->FastGetSolutionStepValue(PARTICLE_RADIUS);
//
// 		(ii_radi>=jj_radi) ? (ptr_begin + ii)->GetValue(SCALE)+=1 : (ptr_begin + jj)->GetValue(SCALE)+=1;
//
// 	      }
// 	      //delete extra nodes
// 	       int aux_ptr_elem_size = int(ptr_elem_size);
// 	      while(aux_ptr_elem_size>56)
// 	      {
// 	        for(unsigned int ii=0; ii < ptr_elem_size; ii++){
// 		  bool swt = false;
// 		  for( int  kkk = ptr_elem_size; kkk>0; kkk-- )
// 		     if( (ptr_begin + ii)->GetValue(SCALE) == kkk && (ptr_begin + ii)->GetValue(IS_VISITED) == 0){
// 		        bool is_escaped = (ptr_begin + ii)->GetValue(IS_ESCAPED);
// 			if( is_escaped==false )
// 		            (ptr_begin + ii)->Set(TO_ERASE,true);//CHECK ESCASPED NODES
// 			(ptr_begin + ii)->SetValue(IS_VISITED,1);
// 			swt = true;
// 			break;
// 		     }
// 		     if(swt )
// 		       break;
// 		}
// 		aux_ptr_elem_size -= 1;
// 	      }
// 	  }
//
//
// 	}
//
//     }

    void ComputeGaussPointPositions2D(Geometry< Node < 3 > >& geom, BoundedMatrix<double, 16, 3 > & pos, BoundedMatrix<double, 16, 3 > & N)
    {
        //lower diagonal terms
        double ypos = 1.0 / 5.0;
        int pos_counter = 0;
        for (unsigned int i = 0; i < 4; i++)
        {
            double xpos = 1.0 / 8.0;
            for (unsigned int j = 0; j < (7-2*i); j++)
            {
                double N1 = xpos;
                double N2 = ypos;
                double N3 = 1.0 - xpos - ypos;

                pos(pos_counter, 0) = N1 * geom[0].X() + N2 * geom[1].X() + N3 * geom[2].X();
                pos(pos_counter, 1) = N1 * geom[0].Y() + N2 * geom[1].Y() + N3 * geom[2].Y();
                pos(pos_counter, 2) = N1 * geom[0].Z() + N2 * geom[1].Z() + N3 * geom[2].Z();

                N(pos_counter, 0) = N1;
                N(pos_counter, 1) = N2;
                N(pos_counter, 2) = N3;

                xpos += 1.0 / 8.0;
                pos_counter += 1;

            }
            ypos += 1.0 / 5.0;
        }

    }

    void ComputeGaussPointPositions3D(Geometry< Node < 3 > >& geom, BoundedMatrix<double, 56, 3 > & pos, BoundedMatrix<double, 56, 4 > & N)
    {
     int pos_counter = 0;
     const double one_seventh = 1.0/6.5;
     double close_point = 1.0/20;
     double zpos = close_point;

     for (unsigned int kk = 0; kk < 6; kk++)
      {
// 	double y_div = 1.0/(7 - kk);
        double ypos =  close_point;//one_seventh;// y_div * (1.0 - zpos);//one_seventh
        for (unsigned int i = 0; i < (6-kk); i++)
        {
// 	    double x_div =  1.0/(7 - kk);// -i
            double xpos = close_point;//one_seventh;//x_div* (1.0 - ypos) * (1.0 - zpos);//one_seventh
            for (unsigned int j = 0; j < (6-kk-i); j++)
            {
                double N1 = xpos;
                double N2 = ypos;
		double N3 = zpos;
                double N4 = 1.0 - xpos - ypos - zpos;

                pos(pos_counter, 0) = N1 * geom[0].X() + N2 * geom[1].X() + N3 * geom[2].X() + N4 * geom[3].X();
                pos(pos_counter, 1) = N1 * geom[0].Y() + N2 * geom[1].Y() + N3 * geom[2].Y() + N4 * geom[3].Y();
                pos(pos_counter, 2) = N1 * geom[0].Z() + N2 * geom[1].Z() + N3 * geom[2].Z() + N4 * geom[3].Z();

                N(pos_counter, 0) = N1;
                N(pos_counter, 1) = N2;
                N(pos_counter, 2) = N3;
                N(pos_counter, 3) = N4;

                xpos += one_seventh;//x_div * (1.0 - ypos) * (1.0 - zpos); //one_seventh
                pos_counter += 1;

            }
            ypos += one_seventh;//y_div * (1.0 - zpos);//one_seventh
        }
        zpos += one_seventh;
      }
    }

    void FewComputeGaussPointPositions3D(Geometry< Node < 3 > >& geom, BoundedMatrix<double, 10, 3 > & pos, BoundedMatrix<double, 10, 4 > & N)
    {
     int pos_counter = 0;
     const double one_third = 1.0/2.5;
     double close_point = 1.0/20;
     double zpos = close_point;

     for (unsigned int kk = 0; kk < 3; kk++)
      {
// 	double y_div = 1.0/(7 - kk);
        double ypos =  close_point;//one_seventh;// y_div * (1.0 - zpos);//one_seventh
        for (unsigned int i = 0; i < (3-kk); i++)
        {
// 	    double x_div =  1.0/(7 - kk);// -i
            double xpos = close_point;//one_seventh;//x_div* (1.0 - ypos) * (1.0 - zpos);//one_seventh
            for (unsigned int j = 0; j < (3-kk-i); j++)
            {
                double N1 = xpos;
                double N2 = ypos;
		double N3 = zpos;
                double N4 = 1.0 - xpos - ypos - zpos;

                pos(pos_counter, 0) = N1 * geom[0].X() + N2 * geom[1].X() + N3 * geom[2].X() + N4 * geom[3].X();
                pos(pos_counter, 1) = N1 * geom[0].Y() + N2 * geom[1].Y() + N3 * geom[2].Y() + N4 * geom[3].Y();
                pos(pos_counter, 2) = N1 * geom[0].Z() + N2 * geom[1].Z() + N3 * geom[2].Z() + N4 * geom[3].Z();

                N(pos_counter, 0) = N1;
                N(pos_counter, 1) = N2;
                N(pos_counter, 2) = N3;
                N(pos_counter, 3) = N4;

                xpos += one_third;//x_div * (1.0 - ypos) * (1.0 - zpos); //one_seventh
                pos_counter += 1;

            }
            ypos += one_third;//y_div * (1.0 - zpos);//one_seventh
        }
        zpos += one_third;
      }
    }


   void RandomPariclePosition(Geometry< Node < 3 > >& geom, BoundedMatrix<double, 16, 3 > & coord, BoundedMatrix<double, 16, 3 > & N_shape)
    {
      for(int ii=0;ii<16;ii++){
      double xi =  rand()* ( 1.0 / ( RAND_MAX + 1.0 ) );
      double etta =  (1.0 - xi) * ( rand()* ( 1.0 / ( RAND_MAX + 1.0 ) ) );
      double zetta = 1.0 - (xi + etta);

      coord(ii,0) = xi * geom[0].X() + etta * geom[1].X() + zetta * geom[2].X();
      coord(ii,1) = xi * geom[0].Y() + etta * geom[1].Y() + zetta * geom[2].Y();
      coord(ii,2) = xi * geom[0].Z() + etta * geom[1].Z() + zetta * geom[2].Z();

      N_shape(ii,0) = xi;
      N_shape(ii,1) = etta;
      N_shape(ii,1) = zetta;
      }
    }

  static int CheckElemDist(Geometry< Node < 3 > >& geom, const double max_dist)
  {
    for(unsigned int ii=0; ii < geom.size(); ++ii)
    {
     double nd_dist = geom[ii].FastGetSolutionStepValue(DISTANCE);
     if (fabs(nd_dist) < max_dist)
       return 1;
    }
    return 0;
  }

  bool CheckIfEscaped(Geometry< Node < 3 > >& geom, const array_1d<double, 3 > & N_shape,const double particle_sign)
  {
    double dist =  N_shape[0]*geom[0].FastGetSolutionStepValue(DISTANCE);
     for(unsigned int ii=1; ii < geom.size(); ++ii)
       dist += N_shape[ii]*geom[ii].FastGetSolutionStepValue(DISTANCE);

   if( dist*particle_sign < 0.0)
     return true;
   else
     return false;
  }

  void AssignParticleRadius(Node < 3 > ::Pointer nd_ptr, double& p_dist,const double min_edge_size)
  {
    if( fabs(p_dist) < 0.1*min_edge_size)
      nd_ptr->FastGetSolutionStepValue(PARTICLE_RADIUS) = 0.1*min_edge_size;
    else if(fabs(p_dist) > 0.5*min_edge_size)
       nd_ptr->FastGetSolutionStepValue(PARTICLE_RADIUS) = 0.5*min_edge_size;
    else
       nd_ptr->FastGetSolutionStepValue(PARTICLE_RADIUS) = fabs(p_dist);

  }

//   unsigned int time_seed()
//   {
//     time_t now = time ( 0 );
//     unsigned char *p = (unsigned char *)&now;
//     unsigned int seed = 0;
//     size_t i;
//
//     for ( i = 0; i < sizeof now; i++ )
//       seed = seed * ( UCHAR_MAX + 2U ) + p[i];
//
//    return seed;
//  }
};

}
#endif // KRATOS_LAGRANGIAN_PARTICLES_UTILITIES_INCLUDED  defined
