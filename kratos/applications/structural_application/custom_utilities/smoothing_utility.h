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
*   Last Modified by:    $Author: Nelson Lafontaine $
*   Date:                $Date: 01-27-2010$
*   Revision:            $Revision: 1.00   $
*
* ***********************************************************/

#if !defined(SMOOTHING_UTILITY)
#define SMOOTHING_UTILITY

#ifdef _OPENMP
#include <omp.h>
#endif

#include "boost/smart_ptr.hpp"
#include <boost/timer.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/lu.hpp>




/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/variables.h"
#include "containers/array_1d.h"
#include "processes/find_nodal_neighbours_process.h"
#include "processes/find_elements_neighbours_process.h"
#include "structural_application.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"
#include <cmath>
#include <algorithm>



namespace Kratos
{

    enum   polynomial_degree {Linear = 1, Cuadratic, Cubic, Cuartic};

    class Smoothing_Utility
    {

        public:


            typedef ModelPart::NodesContainerType NodesArrayType;
            typedef ModelPart::ElementsContainerType ElementsArrayType;
            typedef ModelPart::ConditionsContainerType ConditionsArrayType;
            typedef boost::numeric::ublas::vector<Matrix> Matrix_Order_Tensor;
            typedef boost::numeric::ublas::vector<Vector> Vector_Order_Tensor;


            Smoothing_Utility( ModelPart& model_part, int domain_size ) : mr_model_part( model_part )
            {
                mdomain_size = domain_size;
                minitialize_Setting_Variables = false;
            }

            ~Smoothing_Utility() {}



            inline void CreatePartition( unsigned int number_of_threads, const int number_of_rows, vector<unsigned int>& partitions )
            {
                partitions.resize( number_of_threads + 1 );
                int partition_size = number_of_rows / number_of_threads;
                partitions[0] = 0;
                partitions[number_of_threads] = number_of_rows;

                for ( unsigned int i = 1; i < number_of_threads; i++ )
                    partitions[i] = partitions[i-1] + partition_size ;
            }


           // it computes  the triburary volume or area of each nodes.
	   // only its valid for triangle and tetrahedra
           void CalculatetributaryFactor ( ModelPart& this_model_part,  const unsigned int& domain_size )
           {
	     
	        ElementsArrayType& pElements =  this_model_part.Elements();                         
                
                #ifdef _OPENMP
                int number_of_threads = omp_get_max_threads();
                #else
                int number_of_threads = 1;
                #endif

		vector<unsigned int> element_partition;
                CreatePartition( number_of_threads, pElements.size(), element_partition );
                
		switch(domain_size)
		{
		  case 2:
		  {
		    const double fact = 1.00/3.00; 
                    #pragma omp parallel for
                    for ( int k = 0; k < number_of_threads; k++ )
                     {

		      ElementsArrayType::iterator it_begin = pElements.ptr_begin() + element_partition[k];
		      ElementsArrayType::iterator it_end   = pElements.ptr_begin() + element_partition[k+1];

		      for ( ElementsArrayType::iterator it = it_begin; it != it_end; ++it )
		      {
			  Element::GeometryType& geom = it->GetGeometry(); 
			  for ( unsigned int i = 0; i < geom.size(); i++ )
			  {
			      geom[i].SetLock();
			      double&  Nodal_Area   = geom[i].GetValue(NODAL_AREA);
			      Nodal_Area           += fact  *  geom.Area();
			      geom[i].UnSetLock();
			}
	     
	              }
		   }
		   break;
		  }
		   ///********************************************************************************************
		  case 3: 
		  {           
		    
		    const double fact = 0.25;    
		    #pragma omp parallel for
                    for ( int k = 0; k < number_of_threads; k++ )
                     {

		      ElementsArrayType::iterator it_begin = pElements.ptr_begin() + element_partition[k];
		      ElementsArrayType::iterator it_end   = pElements.ptr_begin() + element_partition[k+1];

		      for ( ElementsArrayType::iterator it = it_begin; it != it_end; ++it )
		      {
			  Element::GeometryType& geom = it->GetGeometry();   
			  for ( unsigned int i = 0; i < geom.size(); i++ )
			  {
			      geom[i].SetLock();
			      double&  Nodal_Volume   = geom[i].GetValue(NODAL_VOLUME);
			      Nodal_Volume += fact  * geom.Volume(); 
			      geom[i].UnSetLock();
			}
	     
	              }
		   }
		   
		    break;
		  }
	      }
	   }
	   

            template<class TVariableType>
            void WeightedRecoveryGradients(const Variable<TVariableType>& rVariable, Variable<TVariableType>& rVariable_Smooth,
            ModelPart& this_model_part, const unsigned int& domain_size )
            {
              if( minitialize_Setting_Variables==false)
	      { 
		CalculatetributaryFactor(this_model_part, domain_size); 
		minitialize_Setting_Variables = true;
	      }
	      
	      switch(domain_size)
		{
		  case 2:  
		  { 
		    WeightedRecoveryGradients2D(rVariable, rVariable_Smooth, this_model_part);    
		    break;
		  }
		  case 3:
		    {
		     WeightedRecoveryGradients3D(rVariable, rVariable_Smooth, this_model_part);
		     break;
		    }
		}
	       
	    }
		
	
		// Only valid for triangle elements
		template<class TVariableType>
		void WeightedRecoveryGradients2D(const Variable<TVariableType>& rVariable, const Variable<TVariableType>& rVariable_Smooth, ModelPart& this_model_part)
		{
		  
		KRATOS_TRY
		
                ProcessInfo& CurrentProcessInfo    =  this_model_part.GetProcessInfo();
                ElementsArrayType& pElements       =  this_model_part.Elements();
                NodesArrayType& pNodes             =  this_model_part.Nodes();
                
                #ifdef _OPENMP
                int number_of_threads = omp_get_max_threads();
                #else
                int number_of_threads = 1;
                #endif

                vector<unsigned int> element_partition;
                CreatePartition( number_of_threads, pElements.size(), element_partition );
		vector<unsigned int> node_partition;
                CreatePartition( number_of_threads, pNodes.size(), node_partition );

                const double fact   = 1.00/3.00;
                std::vector<TVariableType>   Variable_Value;
                #pragma omp parallel for private(Variable_Value) 
                for ( int k = 0; k < number_of_threads; k++ )
                {
  
                    ElementsArrayType::iterator it_begin = pElements.ptr_begin() + element_partition[k];
                    ElementsArrayType::iterator it_end   = pElements.ptr_begin() + element_partition[k+1];

                    for ( ElementsArrayType::iterator it = it_begin; it != it_end; ++it )
                    {
                        Element::GeometryType& geom     = it->GetGeometry();
                        it->GetValueOnIntegrationPoints(rVariable, Variable_Value, CurrentProcessInfo );
                        for ( unsigned int i = 0; i < geom.size(); i++ )
                        {
                            geom[i].SetLock();
                            geom[i].GetValue(rVariable_Smooth) += fact * geom.Area() * Variable_Value[0];  
                            geom[i].UnSetLock();
                        }
                    }
                }
                

                #pragma omp parallel for
                for ( int k = 0; k < number_of_threads; k++ )
                {
                    NodesArrayType::iterator i_begin = pNodes.ptr_begin() + node_partition[k];
                    NodesArrayType::iterator i_end = pNodes.ptr_begin() + node_partition[k+1];

                    for (ModelPart::NodeIterator i = i_begin; i != i_end; ++i )
                    {
                        double& Area_Total           =  i->GetValue(NODAL_AREA );
			TVariableType& Nodal_Values  =  i->GetValue(rVariable_Smooth);
			Nodal_Values                 = (1.00 / Area_Total) * Nodal_Values;
		     }
		   
          	   }
                KRATOS_CATCH( "" )
                
            }

		// Only valid for tetrahedra elements
		// WARNING = NODAL_AREA SHOULD BE NODAL_VOLUME
		template<class TVariableType>
		void WeightedRecoveryGradients3D(const Variable<TVariableType>& rVariable, const Variable<TVariableType>& rVariable_Smooth, ModelPart& this_model_part)
		{

		
		KRATOS_TRY
		
                ProcessInfo& CurrentProcessInfo    =  this_model_part.GetProcessInfo();
                ElementsArrayType& pElements       =  this_model_part.Elements();
                NodesArrayType& pNodes             =  this_model_part.Nodes();
                
                #ifdef _OPENMP
                int number_of_threads = omp_get_max_threads();
                #else
                int number_of_threads = 1;
                #endif

                vector<unsigned int> element_partition;
                CreatePartition( number_of_threads, pElements.size(), element_partition );


                const double fact   = 0.25;
                std::vector<TVariableType>   Variable_Value;
                #pragma omp parallel for private(Variable_Value) 
                for ( int k = 0; k < number_of_threads; k++ )
                {
  
                    ElementsArrayType::iterator it_begin = pElements.ptr_begin() + element_partition[k];
                    ElementsArrayType::iterator it_end   = pElements.ptr_begin() + element_partition[k+1];

                    for ( ElementsArrayType::iterator it = it_begin; it != it_end; ++it )
                    {
                        Element::GeometryType& geom     = it->GetGeometry();
                        it->GetValueOnIntegrationPoints(rVariable, Variable_Value, CurrentProcessInfo );
                        for ( unsigned int i = 0; i < geom.size(); i++ )
                        {
                            geom[i].SetLock();
                            geom[i].GetValue(rVariable_Smooth) += fact * geom.Volume() * Variable_Value[0];         
                            geom[i].UnSetLock();
                        }
                    }
                }
                
                vector<unsigned int> node_partition;
                CreatePartition( number_of_threads, pNodes.size(), node_partition );

                #pragma omp parallel for
                for ( int k = 0; k < number_of_threads; k++ )
                {
                    NodesArrayType::iterator i_begin = pNodes.ptr_begin() + node_partition[k];
                    NodesArrayType::iterator i_end = pNodes.ptr_begin() + node_partition[k+1];

                    for ( ModelPart::NodeIterator i = i_begin; i != i_end; ++i )
                    {
                        double& Vol_Total            =  i->GetValue(NODAL_VOLUME);
			TVariableType& Nodal_Values  =  i->GetValue(rVariable_Smooth);
                        Nodal_Values                 = (1.00 / Vol_Total) * Nodal_Values;

		     }
          	   }
          	   
          	   
                KRATOS_CATCH( "" )
		  
            }

                
            template<class TVariableType>
            void InterpolatedRecoveryGradients(const Variable<TVariableType>& rVariable, Variable<TVariableType>& rVariable_Smooth,
            ModelPart& this_model_part, const unsigned int& domain_size)
            {

	        
                KRATOS_TRY
                NodesArrayType& pNodes             =  this_model_part.Nodes();
		ProcessInfo& CurrentProcessInfo    =  this_model_part.GetProcessInfo();

                #ifdef _OPENMP
                int number_of_threads = omp_get_max_threads();
                #else
                int number_of_threads = 1;
                #endif

                vector<unsigned int> node_partition;
                CreatePartition( number_of_threads, pNodes.size(), node_partition );

                // Variables Globales
                std::vector<TVariableType> Output_Values;
                array_1d<double, 3> Coord_Point = ZeroVector(3);
                array_1d<double, 3> Coord_Node  = ZeroVector(3);
                Vector Polynomial;

                unsigned int size_2 = 3;
                if ( domain_size == 3)
                    size_2 = 6;
               

                bool init = false;
                Matrix Value;
                Matrix_Order_Tensor Aux_b;
                Vector_Order_Tensor Aux_Poly;
                Vector_Order_Tensor a; 
                std::vector<int> work_array_elem;


                    for ( ModelPart::NodesContainerType::iterator i = this_model_part.NodesBegin();
                            i != this_model_part.NodesEnd(); i++ )
                    {
   
                        WeakPointerVector< Element > rneigh_el = i->GetValue( NEIGHBOUR_ELEMENTS ); // elementos vecinos al nodo
                        unsigned int num_of_elem_layer_one     = std::distance(rneigh_el.begin(), rneigh_el.end());

                        // Realizaremos una interpolacion cuadratica. En caso especial, donde el numero de elementos
                        // vecinos sean menor que 10 usaremos como ultima intancia una aproximacion lineal.

                        if ( num_of_elem_layer_one < 13 )
                        {
                            Compute_First_and_Second_Neighbour(i, work_array_elem );
                            num_of_elem_layer_one = work_array_elem.size();  
                            Aux_b.resize( num_of_elem_layer_one );
                            Aux_Poly.resize( num_of_elem_layer_one );

                            for ( unsigned int i_elem = 0; i_elem < work_array_elem.size(); i_elem ++ )
                            {

                                unsigned int elem = work_array_elem[i_elem];
                                Element::GeometryType& geom = ( this_model_part.Elements()[elem] ).GetGeometry();
                                Find_Coord_Gauss_Points( geom , Coord_Point );
                                CalculatePolynomialInterpolation( Polynomial, Coord_Point, num_of_elem_layer_one );
                                unsigned int size_Pol =  Polynomial.size();

                                if ( init == false )
                                {
                                    for ( unsigned k = 0; k < num_of_elem_layer_one; k++ )
                                    {
                                        Aux_Poly[k].resize( size_Pol );
                                        Aux_b[k].resize( 1, size_2 );
                                    }
                                }

                                (this_model_part.Elements()[elem]).GetValueOnIntegrationPoints(rVariable, Output_Values, CurrentProcessInfo ); 
                                noalias( Aux_b[i_elem] )    = Output_Values[0];
                                noalias( Aux_Poly[i_elem] ) = Polynomial;
                                init = true;
                            }

                        }

                        // Interpolacion Cuadratica
                        else
                        {
                            unsigned int counter = 0;
                            Aux_b.resize( num_of_elem_layer_one );
                            Aux_Poly.resize( num_of_elem_layer_one );

                            for ( WeakPointerVector<Element>::iterator it = rneigh_el.begin();
                                    it != rneigh_el.end(); it++ )
                            {

                                Element::GeometryType& geom = it->GetGeometry();
                                Find_Coord_Gauss_Points( geom , Coord_Point );
                                CalculatePolynomialInterpolation( Polynomial, Coord_Point, num_of_elem_layer_one );
                                unsigned int size_Pol =  Polynomial.size();

                                if ( init == false )
                                {
                                    for ( unsigned k = 0; k < num_of_elem_layer_one; k++ )
                                    {
                                        Aux_Poly[k].resize( size_Pol );
                                        Aux_b[k].resize( 1, size_2 );
                                    }
                                }

                                it->GetValueOnIntegrationPoints(rVariable, Output_Values, CurrentProcessInfo ); 
                                noalias( Aux_b[counter] )    = Output_Values[0];
                                noalias( Aux_Poly[counter] ) = Polynomial;
                                init = true;
                                counter++;
                            }
                        }

                        
                        Solve( Aux_b, Aux_Poly, a );
                        Coord_Node[0]  = this_model_part.Nodes()[i->Id()].X();
                        Coord_Node[1]  = this_model_part.Nodes()[i->Id()].Y();
                        Coord_Node[2]  = this_model_part.Nodes()[i->Id()].Z();

                        ///computing average
                        Compute_Interpolated_Sigma(Coord_Node, a, num_of_elem_layer_one, Value);
                        i->GetValue(rVariable_Smooth) = Value;
                        work_array_elem.clear();
                        init = false;
                     
                    }  
                                  
                KRATOS_CATCH( "" )
               
            }

            
            //WARNING = Solo para dos Dimesiones
            void CalculatePolynomialInterpolation( Vector& P, array_1d<double, 3>& Coord_Point, const unsigned int degree )
            {
                unsigned int dim = 0;
                unsigned int degree_aux = 0; // cantidad de puntos disponibles

                if ( degree >= 3 && degree < 15 )    // Linear
                {
                    degree_aux = 1;
                }
                else if ( degree >= 15 && degree < 40 ) // Cuadratic
                {
                    degree_aux = 2;
                }
                else if ( degree >= 40 ) // Cubic
                {
                    degree_aux = 3;
                }
                else if(degree>80)    // Cuartic
                  {
		    degree_aux=4;
		  }
                else
                {
                    std::cout << "Warning: No se puede realizar la interpolacion deseada" << std::endl;
                }

                mpolynomial_degree = static_cast<polynomial_degree>( degree_aux );

                switch ( mpolynomial_degree )
                {

                    case Linear:
                        {
                            dim  = 3;
                            P.resize( dim );
                            P[0] = 1.00;
                            P[1] = Coord_Point[0];
                            P[2] = Coord_Point[1];
                            break;
                        }

                    case Cuadratic:
                        {
                            dim  = 6;
                            P.resize( dim );
                            P[0] = 1.00;
                            P[1] = Coord_Point[0];
                            P[2] = Coord_Point[1];
                            P[3] = Coord_Point[0] * Coord_Point[1];
                            P[4] = Coord_Point[0] * Coord_Point[0];
                            P[5] = Coord_Point[1] * Coord_Point[1];
                            break;
                        }

                    case Cubic:
                        {
                            dim  = 10;
                            P.resize( dim );
                            P[0] = 1.00;
                            P[1] = Coord_Point[0];
                            P[2] = Coord_Point[1];
                            P[3] = Coord_Point[0] * Coord_Point[1];
                            P[4] = Coord_Point[0] * Coord_Point[0];
                            P[5] = Coord_Point[1] * Coord_Point[1];
                            P[6] = Coord_Point[0] * Coord_Point[0] * Coord_Point[1];
                            P[7] = Coord_Point[1] * Coord_Point[1] * Coord_Point[0];
                            P[8] = Coord_Point[0] * Coord_Point[0] * Coord_Point[0];
                            P[9] = Coord_Point[1] * Coord_Point[1] * Coord_Point[1];
                            break;
                        }

                    case Cuartic:
                        {
                            dim  = 15;
                            P.resize( dim );
                            P[0] = 1.00;
                            P[1] = Coord_Point[0];
                            P[2] = Coord_Point[1];
                            P[3] = Coord_Point[0] * Coord_Point[1];
                            P[4] = Coord_Point[0] * Coord_Point[0];
                            P[5] = Coord_Point[1] * Coord_Point[1];
                            P[6] = Coord_Point[0] * Coord_Point[0] * Coord_Point[1];
                            P[7] = Coord_Point[1] * Coord_Point[1] * Coord_Point[0];
                            P[8] = Coord_Point[0] * Coord_Point[0] * Coord_Point[0];
                            P[9] = Coord_Point[1] * Coord_Point[1] * Coord_Point[1];
                            P[10] = Coord_Point[0] * Coord_Point[0] * Coord_Point[1] * Coord_Point[1];
                            P[11] = Coord_Point[0] * Coord_Point[0] * Coord_Point[0] * Coord_Point[1];
                            P[12] = Coord_Point[1] * Coord_Point[1] * Coord_Point[1] * Coord_Point[0];
                            P[13] = Coord_Point[0] * Coord_Point[0] * Coord_Point[0] * Coord_Point[0];
                            P[14] = Coord_Point[1] * Coord_Point[1] * Coord_Point[1] * Coord_Point[1];
                            break;
                        }

                    default:
                        {
                            std::cout << "WARNING: CASE NOT VALID" << std::endl;
                        }
                }
            }


            ///WARNING = Valid for tethaedra and  triangle only
            void Find_Coord_Gauss_Points( Element::GeometryType& geom, array_1d<double, 3>&  Coord_Point )
            {
                double x    = 0.00;
                double y    = 0.00;
                double z    = 0.00;
                double fact = 1.00/3.00;

                if ( geom.size() == 4 )
                {
                    fact = 0.25;
                }

                Coord_Point = ZeroVector( 3 );
                for ( unsigned int i = 0; i < geom.size(); i++ )
                {

                    x = geom[i].X();
                    y = geom[i].Y();
                    z = geom[i].Z();

                    Coord_Point[0] += x;
                    Coord_Point[1] += y;
                    Coord_Point[2] += z;
                }

                noalias( Coord_Point ) = Coord_Point * fact;

            }



            void Compute_First_and_Second_Neighbour( ModelPart::NodesContainerType::iterator& this_node, std::vector<int>& work_array )
            {
                KRATOS_TRY

                work_array.reserve( 1000 );

                WeakPointerVector< Node<3> >& neighb_nodes    = this_node->GetValue( NEIGHBOUR_NODES );
                WeakPointerVector< Element >& neighb_elems    = this_node->GetValue( NEIGHBOUR_ELEMENTS );

                //filling the first neighbours list
                for ( WeakPointerVector<Element>::iterator i = neighb_elems.begin();
                        i != neighb_elems.end(); i++ )
                {
                    int index_j = i->Id();
                    work_array.push_back( index_j );
                }

                //adding the second neighbours
                for ( WeakPointerVector< Node<3> >::iterator i = neighb_nodes.begin();
                        i != neighb_nodes.end(); i++ )
                {
                    WeakPointerVector< Element > second_neighb_elems =  i->GetValue( NEIGHBOUR_ELEMENTS );

                    for ( WeakPointerVector<Element>::iterator  j = second_neighb_elems.begin();
                            j != second_neighb_elems.end(); j++ )
                    {
                        int second_neighb_index = j->Id();
                        work_array.push_back( second_neighb_index );
                    }
                }

                //sorting the indices and elminating the duplicates
                std::sort( work_array.begin(), work_array.end() );

                std::vector<int>::iterator new_end = std::unique( work_array.begin(), work_array.end() );

                unsigned int number_of_entries = new_end - work_array.begin();

                work_array.resize( number_of_entries, false );

                KRATOS_CATCH( "" )
            }


            void SettingNodalValues( ModelPart& this_model_part, const unsigned int& domain_size )
            {
                NodesArrayType& pNodes       =  this_model_part.Nodes();
                
		unsigned int size   = 3;
                if (domain_size == 3 )
                    size  = 6;
                
                const Matrix Nodal_Values   = ZeroMatrix(1, size );
                
                #ifdef _OPENMP
                int number_of_threads = omp_get_max_threads();
                #else
                int number_of_threads = 1;
                #endif
                
                vector<unsigned int> node_partition;
                CreatePartition( number_of_threads, pNodes.size(), node_partition );
                
		switch(domain_size)
		{
		  case 2:
		  {  
		  #pragma omp parallel for 
                  for ( int k = 0; k < number_of_threads; k++ )
                   {
                    NodesArrayType::iterator i_begin = pNodes.ptr_begin() + node_partition[k];
                    NodesArrayType::iterator i_end = pNodes.ptr_begin() + node_partition[k+1];

                    for ( ModelPart::NodeIterator i = i_begin; i != i_end; ++i )
                    {
                        {
			    i->GetValue( NODAL_AREA   )  = 0.00;
			    i->GetValue( NODAL_DAMAGE )  = 0.00;
                            i->GetValue( SPLIT_NODAL  )  = false;
                            i->GetValue( NODAL_STRESS )  = Nodal_Values;
                            i->GetValue( NODAL_STRAIN )  = Nodal_Values;    
                        }
                      }
		    }
		    break;
		  }
		  
		  case 3:
		  {
		  
		  #pragma omp parallel for 
                  for ( int k = 0; k < number_of_threads; k++ )
                   {
                    NodesArrayType::iterator i_begin = pNodes.ptr_begin() + node_partition[k];
                    NodesArrayType::iterator i_end = pNodes.ptr_begin() + node_partition[k+1];

                    for ( ModelPart::NodeIterator i = i_begin; i != i_end; ++i )
                    {
                        {
			    i->GetValue( NODAL_VOLUME)    = 0.00;
                            i->GetValue( NODAL_DAMAGE )   = 0.00;
                            i->GetValue( SPLIT_NODAL )    = false;
                            i->GetValue( NODAL_STRESS )   = Nodal_Values;
                            i->GetValue( NODAL_STRAIN )   = Nodal_Values;
                        }
                      }
		    }
		   break;
		  }

                }

            }

            
            void Solve( Matrix_Order_Tensor& Aux_b, Vector_Order_Tensor& Aux_Poly, Vector_Order_Tensor& Result )
            {

                Matrix A;
                Matrix A_inv;
                Vector b;
                Vector a;


                a     = ZeroVector( Aux_Poly[0].size() );
                b     = ZeroVector( Aux_Poly[0].size() );
                A     = ZeroMatrix( Aux_Poly[0].size(), Aux_Poly[0].size() );
                A_inv = ZeroMatrix( Aux_Poly[0].size(), Aux_Poly[0].size() );
                Result.resize( Aux_b[0].size2() );

                for ( unsigned int j = 0; j < Aux_Poly.size(); j++ )   //Numeros de elementos
                {
                    noalias( A ) = A + outer_prod( Aux_Poly[j], Aux_Poly[j] );
                }

                int singular = SD_MathUtils<double>::InvertMatrix( A, A_inv );

                if ( singular == 1 )
                {
                    KRATOS_WATCH( "MATRIX SINGULAR: MORE POINTS FOR EXTRAPOLATIONS" )
                }

                for ( unsigned i = 0; i < Aux_b[0].size2();i++ )  // Numeros de Sigma
                {
                    Result[i].resize( Aux_Poly[0].size() );

                    for ( unsigned int j = 0; j < Aux_Poly.size(); j++ )   //Numeros de elementos
                    {
                        noalias( b ) = b + Aux_Poly[j] * Aux_b[j]( 0, i );
                    }

                    noalias( a )   = prod( A_inv, b );

                    Result[i]    = a;

                    a  =  ZeroVector( Aux_Poly[0].size() );
                    b  =  ZeroVector( Aux_Poly[0].size() );
                }

            }

            int InvertMatrix( const Matrix& input, Matrix& inverse )
            {
                int singular = 0;

                using namespace boost::numeric::ublas;
                typedef permutation_matrix<std::size_t> pmatrix;
                Matrix A( input );
                pmatrix pm( A.size1() );
                singular = lu_factorize( A, pm );
                inverse.assign( identity_matrix<double>( A.size1() ) );
                lu_substitute( A, pm, inverse );
                return singular;
            }


            void Compute_Interpolated_Sigma(array_1d<double, 3>& Coord_Point,  Vector_Order_Tensor& a, const unsigned int degree, Matrix& Result )
            {
                unsigned int size = a.size();
                Vector Polynomial;
                Vector Aux;

                Result = ZeroMatrix( 1, size );
                CalculatePolynomialInterpolation( Polynomial, Coord_Point, degree );

                for ( unsigned int i = 0; i < size; i++ )
                {
                    Result( 0, i ) = inner_prod(Polynomial, a[i] );
                }
                
            }
             

	    
	    
            // Calcula el area tributaria de los nodos que han sido creados
            void Recompute_Values_For_New_Mesh( ModelPart& this_model_part,  const unsigned int domain_size )
            {
                
                KRATOS_TRY
                
                ElementsArrayType& pElements       =  this_model_part.Elements();

		#ifdef _OPENMP
		int number_of_threads = omp_get_max_threads();
		#else
		int number_of_threads = 1;
		#endif

                vector<unsigned int> element_partition;
                CreatePartition( number_of_threads, pElements.size(), element_partition );

              

                switch(domain_size)
		{
		  case 2:
		  {

		    double Area       = 0.00;
                    double fact       = 1.00/3.00;
                    #pragma omp parallel for private(Area) shared(fact)
                    for ( int k = 0; k < number_of_threads; k++ )
                      {

                        ElementsArrayType::iterator it_begin = pElements.ptr_begin() + element_partition[k];
                        ElementsArrayType::iterator it_end = pElements.ptr_begin() + element_partition[k+1];

                        for ( ElementsArrayType::iterator it = it_begin; it != it_end; ++it )
                        {
			  Element::GeometryType& geom = it->GetGeometry(); // Nodos del elemento
			  Area = geom.Area();

                        for ( unsigned int i = 0; i < geom.size(); i++ )
                        {
                            geom[i].SetLock();
                            double&  Nodal_Area   = geom[i].GetValue(NODAL_AREA );
                            Nodal_Area = Nodal_Area + fact * Area;
                            geom[i].GetValue(NODAL_AREA ) = Nodal_Area;
                            geom[i].UnSetLock();
                        }
                      }
                   }
                   break;
		  }
                   
		  case 3:
                  {
		   double volume     = 0.00;
                   double fact       = 0.25;
		   #pragma omp parallel for private(volume) shared(fact)
                   for ( int k = 0; k < number_of_threads; k++ )
                    {

                    ElementsArrayType::iterator it_begin = pElements.ptr_begin() + element_partition[k];
                    ElementsArrayType::iterator it_end = pElements.ptr_begin() + element_partition[k+1];

                    for ( ElementsArrayType::iterator it = it_begin; it != it_end; ++it )
                    {
                        Element::GeometryType& geom = it->GetGeometry(); // Nodos del elemento
                        volume = geom.Volume();    
                        for ( unsigned int i = 0; i < geom.size(); i++ )
                        {
                            geom[i].SetLock();
                            double&  Nodal_Volume   = geom[i].GetValue( NODAL_VOLUME );
                            Nodal_Volume +=  fact * volume;
                            geom[i].GetValue( NODAL_VOLUME ) = Nodal_Volume;
                            geom[i].UnSetLock();
                        }
                    }
                } 
                
	        break;
              }
 	 }
 	 
                 
                 KRATOS_CATCH( "" )
            }


            void Finalize()
            {
                minitialize_Setting_Variables = false;
            }


            bool minitialize_Setting_Variables;
            polynomial_degree mpolynomial_degree;
            ModelPart& mr_model_part;
            unsigned int mdomain_size;
    };
}

#endif

