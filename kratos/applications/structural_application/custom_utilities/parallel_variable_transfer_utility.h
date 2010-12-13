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
*   Last Modified by:    $Author: nelson $
*   Date:                $Date: 2008-12-09 15:23:36 $
*   Revision:            $Revision: 1.2 $
*
* ***********************************************************/

#if !defined(KRATOS_PARALLEL_VARIABLE_TRANSFER_UTILITY_INCLUDED )
#define  KRATOS_PARALLEL_VARIABLE_TRANSFER_UTILITY_INCLUDED
//System includes
//External includes
#include "boost/smart_ptr.hpp"

//Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "containers/array_1d.h"
#include "includes/element.h"
#include "integration/integration_point.h"
#include "geometries/geometry.h"
#include "linear_solvers/skyline_lu_factorization_solver.h"
#include "spaces/ublas_space.h"
#include "spaces/parallel_ublas_space.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/tetrahedra_3d_4.h"

namespace Kratos
{
    class ParallelVariableTransferUtility
    {
        public:
            typedef Dof<double> TDofType;
            typedef PointerVectorSet<TDofType, IndexedObject> DofsArrayType;
            typedef ModelPart::ElementsContainerType ElementsArrayType;
            typedef double* ContainerType;
            typedef Element::DofsVectorType DofsVectorType;
            typedef Geometry<Node<3> >::IntegrationPointsArrayType IntegrationPointsArrayType;
            typedef Geometry<Node<3> >::GeometryType GeometryType;
            typedef Geometry<Node<3> >::CoordinatesArrayType CoordinatesArrayType;
            typedef ParallelUblasSpace<double, CompressedMatrix, Vector> SpaceType;
            typedef UblasSpace<double, Matrix, Vector> ParallelLocalSpaceType;
            
			/** 
             * Constructor.
                         */
            ParallelVariableTransferUtility()
            {
                std::cout << "ParallelVariableTransferUtility created" << std::endl;
            }
			
			/** 
             * Destructor.
                         */
            virtual ~ParallelVariableTransferUtility()
            {}
			
			/**
             * Initializes elements of target model part.
             * @param rTarget new/target model part
             * KLUDGE: new model part instance is not automatically initialized 
                         */
            void InitializeModelPart( ModelPart& rTarget )
            {
                for( ModelPart::ElementIterator it = rTarget.ElementsBegin();
                     it!= rTarget.ElementsEnd(); it++ )
                {
                    (*it).Initialize();
                }
            }
			
			/**
             * Transfer of nodal solution step variables.
             * This Transfers all solution step variables from r_old_model_part
             * to r_new_model_part.
             * To cope with moved meshes, the source model_part is resetted to its
             * reference configuration temporarily!
             * @param r_old_model_part source model_part
             * @param r_new_model_part target model_part
             * TODO: find more elegant way to check existence of variables in each node
                         */
            void TransferNodalVariables(ModelPart& rSource, ModelPart& rTarget)
            {
				//reset source model part to reference configuration
                for( ModelPart::NodeIterator it = rSource.NodesBegin() ;
                     it != rSource.NodesEnd(); it++ )
                {
                    (*it).X() = (*it).X0();
                    (*it).Y() = (*it).Y0();
                    (*it).Z() = (*it).Z0();
                }
				//reset target model part to reference configuration
                for( ModelPart::NodeIterator it = rTarget.NodesBegin() ;
                     it != rTarget.NodesEnd(); it++ )
                {
                    (*it).X() = (*it).X0();
                    (*it).Y() = (*it).Y0();
                    (*it).Z() = (*it).Z0();
                }
				//time_target= time_source
                ProcessInfo SourceCurrentProcessInfo= rSource.GetProcessInfo();
                rTarget.CloneTimeStep(SourceCurrentProcessInfo[TIME]);

                ElementsArrayType& OldMeshElementsArray= rSource.Elements();
                Element correspondingElement;
//				FixDataValueContainer newNodalValues;
//				FixDataValueContainer oldNodalValues;
                Point<3>  localPoint;

                for(ModelPart::NodeIterator it = rTarget.NodesBegin() ; 
                    it != rTarget.NodesEnd() ; it++)
                {
                    if(FindPartnerElement(*(it), OldMeshElementsArray,
                       correspondingElement,localPoint))
                    {
						//TransferVariables from Old Mesh to new Node
                        if(it->HasDofFor(DISPLACEMENT_X)
                           || it->HasDofFor(DISPLACEMENT_Y)
                           || it->HasDofFor(DISPLACEMENT_Z))
                        {
                            noalias(it->GetSolutionStepValue(DISPLACEMENT_NULL))= 
                                    MappedValue(correspondingElement,
                                    localPoint,DISPLACEMENT_NULL );
                            noalias(it->GetSolutionStepValue(DISPLACEMENT_EINS))= 
                                    MappedValue(correspondingElement,
                                    localPoint,DISPLACEMENT_EINS );
                            noalias(it->GetSolutionStepValue(DISPLACEMENT_NULL_DT))= 
                                    MappedValue(correspondingElement,
                                    localPoint,DISPLACEMENT_NULL_DT );
                            noalias(it->GetSolutionStepValue(ACCELERATION_NULL))=
                                    MappedValue(correspondingElement,
                                    localPoint,ACCELERATION_NULL );
                            noalias(it->GetSolutionStepValue(DISPLACEMENT_OLD))= 
                                    MappedValue(correspondingElement,
                                    localPoint,DISPLACEMENT_OLD );
                        }
                        if(it->HasDofFor(WATER_PRESSURE))
                        {
                            it->GetSolutionStepValue(WATER_PRESSURE_NULL)= 
                                    MappedValuePressure(correspondingElement, localPoint,
                                    WATER_PRESSURE_NULL);
                            it->GetSolutionStepValue(WATER_PRESSURE_EINS)= 
                                    MappedValuePressure(correspondingElement, localPoint,
                                    WATER_PRESSURE_EINS);
                            it->GetSolutionStepValue(WATER_PRESSURE_NULL_DT)= 
                                    MappedValuePressure(correspondingElement, localPoint,
                                    WATER_PRESSURE_NULL_DT);
                            it->GetSolutionStepValue(WATER_PRESSURE_NULL_ACCELERATION)= 
                                    MappedValuePressure(correspondingElement, localPoint,
                                    WATER_PRESSURE_NULL_ACCELERATION);
                        }
                        if(it->HasDofFor(AIR_PRESSURE))
                        {
                            it->GetSolutionStepValue(AIR_PRESSURE_NULL)= 
                                    MappedValuePressure(correspondingElement, localPoint,
                                    AIR_PRESSURE_NULL);
                            it->GetSolutionStepValue(AIR_PRESSURE_EINS)= 
                                    MappedValuePressure(correspondingElement, localPoint,
                                    AIR_PRESSURE_EINS);
                            it->GetSolutionStepValue(AIR_PRESSURE_NULL_DT)= 
                                    MappedValuePressure(correspondingElement, localPoint,
                                    AIR_PRESSURE_NULL_DT);
                            it->GetSolutionStepValue(AIR_PRESSURE_NULL_ACCELERATION)= 
                                    MappedValuePressure(correspondingElement, localPoint,
                                    AIR_PRESSURE_NULL_ACCELERATION);
                        }
                        std::cout <<"VARIABLES TRANSFERRED" << std::endl;
                    }
                    else
                    {	
                        std::cout<<"###### NO PARTNER FOUND IN OLD MESH : TransferNodalVariables(...)#####"<<std::endl;
                    }
                }
				//restore source model_part
                for( ModelPart::NodeIterator it = rSource.NodesBegin() ;
                     it != rSource.NodesEnd(); it++ )
                {
                    (*it).X() = (*it).X0()+(*it).GetSolutionStepValue( DISPLACEMENT_X );
                    (*it).Y() = (*it).Y0()+(*it).GetSolutionStepValue( DISPLACEMENT_Y );
                    (*it).Z() = (*it).Z0()+(*it).GetSolutionStepValue( DISPLACEMENT_Z );
                }
				//restore target model_part
                for( ModelPart::NodeIterator it = rTarget.NodesBegin() ;
                     it != rTarget.NodesEnd(); it++ )
                {
                    (*it).X() = (*it).X0()+(*it).GetSolutionStepValue( DISPLACEMENT_X );
                    (*it).Y() = (*it).Y0()+(*it).GetSolutionStepValue( DISPLACEMENT_Y );
                    (*it).Z() = (*it).Z0()+(*it).GetSolutionStepValue( DISPLACEMENT_Z );
                }
            }
			
			/**
             * Transfer of INSITU_STRESS.
             * This transfers the in-situ stress from rSource to rTarget.
             * @param rSource the source model part
             * @param rTarget the target model part
                         */
            void TransferInSituStress( ModelPart& rSource, ModelPart& rTarget )
            {
				//reset original model part to reference configuration
                for( ModelPart::NodeIterator it = rSource.NodesBegin() ;
                     it != rSource.NodesEnd(); it++ )
                {
                    (*it).X() = (*it).X0();
                    (*it).Y() = (*it).Y0();
                    (*it).Z() = (*it).Z0();
                }
                for( ModelPart::NodeIterator it = rTarget.NodesBegin() ;
                     it != rTarget.NodesEnd(); it++ )
                {
                    (*it).X() = (*it).X0();
                    (*it).Y() = (*it).Y0();
                    (*it).Z() = (*it).Z0();
                }

                TransferVariablesToNodes(rSource, INSITU_STRESS);
// 				TransferVariablesBetweenMeshes(rSource, rTarget,INSITU_STRESS);

// 				TransferVariablesToGaussPoints(rTarget, INSITU_STRESS);
                TransferVariablesToGaussPoints(rSource, rTarget, INSITU_STRESS);
				//restore model_part
                for( ModelPart::NodeIterator it = rSource.NodesBegin() ;
                     it != rSource.NodesEnd(); it++ )
                {
                    (*it).X() = (*it).X0()+(*it).GetSolutionStepValue( DISPLACEMENT_X );
                    (*it).Y() = (*it).Y0()+(*it).GetSolutionStepValue( DISPLACEMENT_Y );
                    (*it).Z() = (*it).Z0()+(*it).GetSolutionStepValue( DISPLACEMENT_Z );
                }
                for( ModelPart::NodeIterator it = rTarget.NodesBegin() ;
                     it != rTarget.NodesEnd(); it++ )
                {
                    (*it).X() = (*it).X0()+(*it).GetSolutionStepValue( DISPLACEMENT_X );
                    (*it).Y() = (*it).Y0()+(*it).GetSolutionStepValue( DISPLACEMENT_Y );
                    (*it).Z() = (*it).Z0()+(*it).GetSolutionStepValue( DISPLACEMENT_Z );
                }
            }
			
			/**
             * Transfer of integration point variables.
             * This Transfers all variables on integration points from r_old_model_part
             * to r_new_model_part.
             * To cope with moved meshes, the source model_part is resetted to its
             * reference configuration temporarily!
             * @param r_old_model_part source model_part
             * @param r_new_model_part target model_part
             * TODO: find more elegant way to check existence of variables in each node
             * CAUTION: THIS MAY CREATE VARIABLES ON NODES THAT MIGHT CAUSE A SEGMENTATION
             * FAULT ON RUNTIME
                         */
            void TransferConstitutiveLawVariables(ModelPart& rSource, ModelPart& rTarget)
            {
				//reset source model part to reference configuration
                for( ModelPart::NodeIterator it = rSource.NodesBegin() ;
                     it != rSource.NodesEnd(); it++ )
                {
                    (*it).X() = (*it).X0();
                    (*it).Y() = (*it).Y0();
                    (*it).Z() = (*it).Z0();
                }
				//reset target model part to reference configuration
                for( ModelPart::NodeIterator it = rTarget.NodesBegin() ;
                     it != rTarget.NodesEnd(); it++ )
                {
                    (*it).X() = (*it).X0();
                    (*it).Y() = (*it).Y0();
                    (*it).Z() = (*it).Z0();
                }


                TransferVariablesToNodes(rSource, ELASTIC_LEFT_CAUCHY_GREEN_OLD);

// 				TransferVariablesBetweenMeshes(rSource, rTarget,ELASTIC_LEFT_CAUCHY_GREEN_OLD);
                // 
// 				TransferVariablesToGaussPoints(rTarget, ELASTIC_LEFT_CAUCHY_GREEN_OLD);

                TransferVariablesToGaussPoints( rSource, rTarget, ELASTIC_LEFT_CAUCHY_GREEN_OLD);

                for( ModelPart::NodeIterator it = rSource.NodesBegin() ;
                     it != rSource.NodesEnd(); it++ )
                {
                    (*it).X() = (*it).X0()+(*it).GetSolutionStepValue( DISPLACEMENT_X );
                    (*it).Y() = (*it).Y0()+(*it).GetSolutionStepValue( DISPLACEMENT_Y );
                    (*it).Z() = (*it).Z0()+(*it).GetSolutionStepValue( DISPLACEMENT_Z );
                }
// 				restore target model_part
                for( ModelPart::NodeIterator it = rTarget.NodesBegin() ;
                     it != rTarget.NodesEnd(); it++ )
                {
                    (*it).X() = (*it).X0()+(*it).GetSolutionStepValue( DISPLACEMENT_X );
                    (*it).Y() = (*it).Y0()+(*it).GetSolutionStepValue( DISPLACEMENT_Y );
                    (*it).Z() = (*it).Z0()+(*it).GetSolutionStepValue( DISPLACEMENT_Z );
                }
            }
	
			/**
             * Transfer of rThisVariable stored on nodes to integration point via
             * approximation by shape functions
             * @param model_part model_part on which the transfer should be done
             * @param rThisVariable Matrix-Variable which should be transferred
             * @see TransferVariablesToGaussPoints(ModelPart& model_part, 
            Variable<Kratos::Vector>& rThisVariable)
             * @see TransferVariablesToGaussPoints(ModelPart& model_part, 
            Variable<double>& rThisVariable)
                         */
            void TransferVariablesToGaussPoints(ModelPart& model_part, 
                    Variable<Kratos::Matrix>& rThisVariable)
            {

                ElementsArrayType& ElementsArray= model_part.Elements();

                int number_of_threads = omp_get_max_threads();
                vector<unsigned int> element_partition;
                CreatePartition(number_of_threads, ElementsArray.size(), element_partition);
                KRATOS_WATCH( number_of_threads );
                KRATOS_WATCH( element_partition );


                double start_prod = omp_get_wtime();
                
                #pragma omp parallel for 
                for(int k=0; k<number_of_threads; k++)
                {
                    ElementsArrayType::ptr_iterator it_begin=ElementsArray.ptr_begin()+element_partition[k];
                    ElementsArrayType::ptr_iterator it_end=ElementsArray.ptr_begin()+element_partition[k+1];
                    
                    for( ElementsArrayType::ptr_iterator it = it_begin; it != it_end; ++it )
                    {	
                        const IntegrationPointsArrayType& integration_points =
                                (*it)->GetGeometry().IntegrationPoints((*it)->GetIntegrationMethod());

                        std::vector<Matrix> ValuesOnIntPoint(integration_points.size());

                        const Matrix& Ncontainer = (*it)->GetGeometry().ShapeFunctionsValues((*it)->GetIntegrationMethod());

                        for( unsigned int PointNumber = 0; 
                             PointNumber<integration_points.size(); 
                             PointNumber++)
                        {
                            ValuesOnIntPoint[PointNumber].resize(3,3,false);

                            noalias(ValuesOnIntPoint[PointNumber])= ZeroMatrix(3,3);

                            for(unsigned int node= 0; node< (*it)->GetGeometry().size(); node++)
                            {

                                ValuesOnIntPoint[PointNumber]
                                        +=Ncontainer(PointNumber, node)*
                                        (*it)->GetGeometry()[node].GetSolutionStepValue(rThisVariable);
                            }
                        }

                        (*it)->SetValueOnIntegrationPoints( rThisVariable, ValuesOnIntPoint,
                        model_part.GetProcessInfo());
                    }
                }
                
                double stop_prod = omp_get_wtime();
                std::cout << "transfer time: " << stop_prod - start_prod << std::endl;
            }
			
			/**
             * Transfer of rThisVariable stored on nodes to integration point via
             * approximation by shape functions
             * @param model_part model_part on which the transfer should be done
             * @param rThisVariable Vector-Variable which should be transferred
             * @see TransferVariablesToGaussPoints(ModelPart& model_part, 
            Variable<Kratos::Matrix>& rThisVariable)
             * @see TransferVariablesToGaussPoints(ModelPart& model_part, 
            Variable<double>& rThisVariable)
                         */
            void TransferVariablesToGaussPoints(ModelPart& model_part, 
                    Variable<Kratos::Vector>& rThisVariable)
            {
                ElementsArrayType& ElementsArray= model_part.Elements();

                int number_of_threads = omp_get_max_threads();
                vector<unsigned int> element_partition;
                CreatePartition(number_of_threads, ElementsArray.size(), element_partition);
                KRATOS_WATCH( number_of_threads );
                KRATOS_WATCH( element_partition );


                double start_prod = omp_get_wtime();
                
                #pragma omp parallel for 
                for(int k=0; k<number_of_threads; k++)
                {
                    ElementsArrayType::ptr_iterator it_begin=ElementsArray.ptr_begin()+element_partition[k];
                    ElementsArrayType::ptr_iterator it_end=ElementsArray.ptr_begin()+element_partition[k+1];
                    
                    for( ElementsArrayType::ptr_iterator it = it_begin; it != it_end; ++it )
                    {
                        unsigned int NodesDispMin= 1;
                        unsigned int NodesDispMax= (*it)->GetGeometry().size();
					
                        const IntegrationPointsArrayType& integration_points =
                                (*it)->GetGeometry().IntegrationPoints((*it)->GetIntegrationMethod());

                        std::vector<Vector> ValuesOnIntPoint(integration_points.size());

                        const Matrix& Ncontainer = (*it)->GetGeometry().ShapeFunctionsValues((*it)->GetIntegrationMethod());

                        for( unsigned int PointNumber = 0; 
                             PointNumber<integration_points.size(); 
                             PointNumber++)
                        {
                            ValuesOnIntPoint[PointNumber].resize(6,false);

                            noalias(ValuesOnIntPoint[PointNumber])= ZeroVector(6);

                            for(unsigned int node= NodesDispMin-1; node< NodesDispMax; node++)
                            {
                                ValuesOnIntPoint[PointNumber]
                                        +=Ncontainer(PointNumber, node)*
                                        (*it)->GetGeometry()[node].GetSolutionStepValue(rThisVariable);
                            }
                        }

                        (*it)->SetValueOnIntegrationPoints( rThisVariable, ValuesOnIntPoint,
                          model_part.GetProcessInfo());
                    }
                }
                
                double stop_prod = omp_get_wtime();
                std::cout << "transfer time: " << stop_prod - start_prod << std::endl;
            }

			/**
             * Transfer of rThisVariable stored on nodes to integration point via
             * approximation by shape functions
             * @param model_part model_part on which the transfer should be done
             * @param rThisVariable double-Variable which should be transferred
             * @see TransferVariablesToGaussPoints(ModelPart& model_part, 
            Variable<Kratos::Matrix>& rThisVariable)
             * @see TransferVariablesToGaussPoints(ModelPart& model_part, 
            Variable<Kratos::Vector>& rThisVariable)
                         */
            void TransferVariablesToGaussPoints(ModelPart& model_part, 
                    Variable<double>& rThisVariable)
            {
                ElementsArrayType& ElementsArray= model_part.Elements();
				
                int number_of_threads = omp_get_max_threads();
                vector<unsigned int> element_partition;
                CreatePartition(number_of_threads, ElementsArray.size(), element_partition);
                KRATOS_WATCH( number_of_threads );
                KRATOS_WATCH( element_partition );


                double start_prod = omp_get_wtime();
                
                #pragma omp parallel for 
                for(int k=0; k<number_of_threads; k++)
                {
                    ElementsArrayType::ptr_iterator it_begin=ElementsArray.ptr_begin()+element_partition[k];
                    ElementsArrayType::ptr_iterator it_end=ElementsArray.ptr_begin()+element_partition[k+1];
                    
                    for( ElementsArrayType::ptr_iterator it = it_begin; it != it_end; ++it )
                    {
                        unsigned int NodesDispMin= 1;
                        unsigned int NodesDispMax= (*it)->GetGeometry().size();

                        const IntegrationPointsArrayType& integration_points =
                                (*it)->GetGeometry().IntegrationPoints((*it)->GetIntegrationMethod());
                        std::vector<double> ValuesOnIntPoint(integration_points.size());
					
                        const Matrix& Ncontainer = (*it)->GetGeometry().ShapeFunctionsValues((*it)->GetIntegrationMethod());
					
                        for( unsigned int PointNumber = 0; 
                             PointNumber<integration_points.size(); 
                             PointNumber++)
                        {
                            ValuesOnIntPoint[PointNumber]= 0.0;

                            for(unsigned int node= NodesDispMin-1; node< NodesDispMax; node++)
                            {
                                ValuesOnIntPoint[PointNumber]
                                        +=Ncontainer(PointNumber, node)*
                                        (*it)->GetGeometry()[node].GetSolutionStepValue(rThisVariable);
                            }
                        }
                        (*it)->SetValueOnIntegrationPoints( rThisVariable, ValuesOnIntPoint,
                          model_part.GetProcessInfo());
                    }
                }
                
                double stop_prod = omp_get_wtime();
                std::cout << "transfer time: " << stop_prod - start_prod << std::endl;
            }
			/**
             * Transfer of rThisVariable stored on nodes in source mesh to integration point of target 
             * mesh via approximation by shape functions
             * @param rSource 
             * @param rTarget 
             * @param rThisVariable Vector-Variable which should be transferred
             * @see TransferVariablesToGaussPoints(ModelPart& source_model_part, ModelPart&
            source_model_part, Variable<Kratos::Matrix>& rThisVariable)
             * @see TransferVariablesToGaussPoints(ModelPart& source_model_part, ModelPart&
            source_model_part, Variable<double>& rThisVariable)
                         */
            void TransferVariablesToGaussPoints(ModelPart& rSource, ModelPart& rTarget, 
                    Variable<Kratos::Vector>& rThisVariable)
            {
                ElementsArrayType& SourceMeshElementsArray= rSource.Elements();
                ElementsArrayType& TargetMeshElementsArray= rTarget.Elements();

                int number_of_threads = omp_get_max_threads();
                vector<unsigned int> element_partition;
                CreatePartition(number_of_threads, TargetMeshElementsArray.size(), element_partition);
                KRATOS_WATCH( number_of_threads );
                KRATOS_WATCH( element_partition );


                double start_prod = omp_get_wtime();
                
                #pragma omp parallel for 
                for(int k=0; k<number_of_threads; k++)
                {
                    ElementsArrayType::ptr_iterator it_begin=TargetMeshElementsArray.ptr_begin()+element_partition[k];
                    ElementsArrayType::ptr_iterator it_end=TargetMeshElementsArray.ptr_begin()+element_partition[k+1];
                    
                    for( ElementsArrayType::ptr_iterator it = it_begin; it != it_end; ++it )
                    {
                        const IntegrationPointsArrayType& integration_points 
                                = (*it)->GetGeometry().IntegrationPoints((*it)->GetIntegrationMethod());

                        std::vector<Vector> ValuesOnIntPoint(integration_points.size());

                        for(unsigned int point=0; point< integration_points.size(); point++)
                        {
                            Point<3> sourceLocalPoint;
                            Point<3> targetLocalPoint;
                            noalias(targetLocalPoint)= integration_points[point];
                            Point<3> targetGlobalPoint;
                            (*it)->GetGeometry().GlobalCoordinates(targetGlobalPoint,targetLocalPoint);
                            Element sourceElement;
                        //Calculate Value of rVariable(firstvalue, secondvalue) in OldMesh
                            ValuesOnIntPoint[point].resize(6,false);
                            if(FindPartnerElement(targetGlobalPoint, SourceMeshElementsArray,
                               sourceElement,sourceLocalPoint))
                            {
                                noalias(ValuesOnIntPoint[point])=
                                        ValueVectorInOldMesh(sourceElement, sourceLocalPoint, rThisVariable );
                            }
                        }
                        (*it)->SetValueOnIntegrationPoints( rThisVariable, ValuesOnIntPoint,
                          rTarget.GetProcessInfo());
                    }
                }
                
                double stop_prod = omp_get_wtime();
                std::cout << "transfer time: " << stop_prod - start_prod << std::endl;
            }
            
            /**
             * Transfer of rThisVariable stored on nodes in source mesh to integration point of target 
             * mesh via approximation by shape functions
             * @param rSource 
             * @param rTarget 
             * @param rThisVariable Matrix-Variable which should be transferred
             * @see TransferVariablesToGaussPoints(ModelPart& source_model_part, 
             *      ModelPart& source_model_part, Variable<Kratos::Vector>& rThisVariable)
             * @see TransferVariablesToGaussPoints(ModelPart& source_model_part, 
             *      ModelPart& source_model_part, Variable<double>& rThisVariable)
             */
           void TransferVariablesToGaussPoints(ModelPart& rSource, 
            	ModelPart& rTarget, 
                Variable<Matrix> const& rThisVariable)
            {
                ElementsArrayType const& SourceMeshElementsArray= rSource.Elements();
                ElementsArrayType& TargetMeshElementsArray= rTarget.Elements();
                
                int number_of_threads = omp_get_max_threads();
                vector<unsigned int> element_partition;
                CreatePartition(number_of_threads, TargetMeshElementsArray.size(), element_partition);
                KRATOS_WATCH( number_of_threads );
                KRATOS_WATCH( element_partition );

                double start_prod = omp_get_wtime();
                
#pragma omp parallel for 
                for(int k=0; k<number_of_threads; k++)
                {
                    ElementsArrayType::ptr_iterator it_begin=TargetMeshElementsArray.ptr_begin()+element_partition[k];
                    ElementsArrayType::ptr_iterator it_end=TargetMeshElementsArray.ptr_begin()+element_partition[k+1];
                    
                    for( ElementsArrayType::ptr_iterator it = it_begin; it != it_end; ++it )
                    {
                        const IntegrationPointsArrayType& integration_points 
                                = (*it)->GetGeometry().IntegrationPoints((*it)->GetIntegrationMethod());
                        std::vector<Matrix> ValuesOnIntPoint(integration_points.size());

                        for(unsigned int point=0; point< integration_points.size(); point++)
                        {
                            Point<3> sourceLocalPoint;
                            Point<3> targetLocalPoint;
                            noalias(targetLocalPoint)= integration_points[point];
                            Point<3> targetGlobalPoint;
                            (*it)->GetGeometry().GlobalCoordinates(targetGlobalPoint,targetLocalPoint);
                            Element sourceElement;
						//Calculate Value of rVariable(firstvalue, secondvalue) in OldMesh
                            if(FindPartnerElement(targetGlobalPoint, SourceMeshElementsArray,
                               sourceElement, sourceLocalPoint))
                            {
                                ValuesOnIntPoint[point].resize(3,3,false);
                                ValuesOnIntPoint[point] = ZeroMatrix(3,3);
				ValuesOnIntPoint[point] = ValueMatrixInOldMesh(sourceElement, sourceLocalPoint, rThisVariable );
                            }
                        }
                        (*it)->SetValueOnIntegrationPoints( rThisVariable, ValuesOnIntPoint,
                          rTarget.GetProcessInfo());
                    }
                }
                
                double stop_prod = omp_get_wtime();
                std::cout << "transfer time: " << stop_prod - start_prod << std::endl;
            }
           
/**
             * Transfer of rThisVariable stored on nodes in source mesh to integration point of target 
             * mesh via approximation by shape functions
             * @param rSource 
             * @param rTarget 
             * @param rThisVariable double-Variable which should be transferred
             * @see TransferVariablesToGaussPoints(ModelPart& source_model_part, ModelPart&
            source_model_part, Variable<Kratos::Matrix>& rThisVariable)
             * @see TransferVariablesToGaussPoints(ModelPart& source_model_part, ModelPart&
            source_model_part, Variable<Kratos::Vector>& rThisVariable)
 */
            void TransferVariablesToGaussPoints(ModelPart& rSource, ModelPart& rTarget, 
                    Variable<double>& rThisVariable)
            {
                ElementsArrayType& SourceMeshElementsArray= rSource.Elements();
                ElementsArrayType& TargetMeshElementsArray= rTarget.Elements();

                int number_of_threads = omp_get_max_threads();
                vector<unsigned int> element_partition;
                CreatePartition(number_of_threads, TargetMeshElementsArray.size(), element_partition);
                KRATOS_WATCH( number_of_threads );
                KRATOS_WATCH( element_partition );


                double start_prod = omp_get_wtime();
                
                #pragma omp parallel for 
                for(int k=0; k<number_of_threads; k++)
                {
                    ElementsArrayType::ptr_iterator it_begin=TargetMeshElementsArray.ptr_begin()+element_partition[k];
                    ElementsArrayType::ptr_iterator it_end=TargetMeshElementsArray.ptr_begin()+element_partition[k+1];
                    
                    for( ElementsArrayType::ptr_iterator it = it_begin; it != it_end; ++it )
                    {
                        const IntegrationPointsArrayType& integration_points 
                                = (*it)->GetGeometry().IntegrationPoints((*it)->GetIntegrationMethod());

                        std::vector<double> ValuesOnIntPoint(integration_points.size());

                        for(unsigned int point=0; point< integration_points.size(); point++)
                        {
                            Point<3> sourceLocalPoint;
                            Point<3> targetLocalPoint;
                            noalias(targetLocalPoint)= integration_points[point];
                            Point<3> targetGlobalPoint;
                            (*it)->GetGeometry().GlobalCoordinates(targetGlobalPoint,targetLocalPoint);
                            Element sourceElement;
						//Calculate Value of rVariable(firstvalue, secondvalue) in OldMesh
                            if(FindPartnerElement(targetGlobalPoint, SourceMeshElementsArray,
                               sourceElement,sourceLocalPoint))
                            {
                                ValuesOnIntPoint[point]=
                                        MappedValue(sourceElement, sourceLocalPoint, rThisVariable );
                            }
                        }
                        (*it)->SetValueOnIntegrationPoints( rThisVariable, ValuesOnIntPoint,
                          rTarget.GetProcessInfo());
                    }
                }
                
                double stop_prod = omp_get_wtime();
                std::cout << "transfer time: " << stop_prod - start_prod << std::endl;
            }
            
			/**
             * Transfer of rThisVariable defined on integration points to corresponding 
             * nodal values. The transformation is done in a form that ensures a minimization
             * of L_2-norm error (/sum{rThisVariable- f(x)) whereas 
             * f(x)= /sum{shape_func_i*rThisVariable_i}
             * @param model_part model_part on which the transfer should be done
             * @param rThisVariable Matrix-Variable which should be transferred
             * @see TransferVariablesToNodes(ModelPart& model_part, Variable<Kratos::Vector>& rThisVariable)
             * @see TransferVariablesToNodes(ModelPart& model_part, Variable<double>& rThisVariable)
             * @ref Jiao&Heath: "Common-refinement-based data transfer...", Int.
             * Journal for numer. meth. in eng. 61 (2004) 2402--2427
             * WARNING: this may cause segmentation faults as the respective variables
             * will be created on nodal level while they are originally intended to be 
             * stored on integration points!
                         */
            void TransferVariablesToNodes(ModelPart& model_part, Variable<Kratos::Matrix>& rThisVariable)
            {
                ElementsArrayType& ElementsArray= model_part.Elements();

				//loop over all master surfaces (global search)
                for(ModelPart::NodeIterator it = model_part.NodesBegin();
                    it != model_part.NodesEnd() ; it++)
                {
                    it->GetSolutionStepValue(rThisVariable)
                            = ZeroMatrix(3,3);
                }
				//SetUpEquationSystem
                SpaceType::MatrixType M(model_part.NumberOfNodes(),model_part.NumberOfNodes());
                SpaceType::VectorType g(model_part.NumberOfNodes());
                SpaceType::VectorType b(model_part.NumberOfNodes());
                noalias(M)= ZeroMatrix(model_part.NumberOfNodes(),model_part.NumberOfNodes());
                
                //create a partition of the element array
                int number_of_threads = omp_get_max_threads();
                vector<unsigned int> element_partition;
                CreatePartition(number_of_threads, ElementsArray.size(), element_partition);
                KRATOS_WATCH( number_of_threads );
                KRATOS_WATCH( element_partition );


                double start_prod = omp_get_wtime();
                
                #pragma omp parallel for 
                for(int k=0; k<number_of_threads; k++)
                {
                    ElementsArrayType::ptr_iterator it_begin=ElementsArray.ptr_begin()+element_partition[k];
                    ElementsArrayType::ptr_iterator it_end=ElementsArray.ptr_begin()+element_partition[k+1];
                    
                    for( ElementsArrayType::ptr_iterator it = it_begin; it != it_end; ++it )
                    {
                        const IntegrationPointsArrayType& integration_points 
                                = (*it)->GetGeometry().IntegrationPoints((*it)->GetIntegrationMethod());

                        GeometryType::JacobiansType J(integration_points.size());
                        J = (*it)->GetGeometry().Jacobian(J, (*it)->GetIntegrationMethod());
                        
                        const Matrix& Ncontainer = (*it)->GetGeometry().ShapeFunctionsValues((*it)->GetIntegrationMethod());
                        
                        Matrix InvJ(3,3);
                        double DetJ;
                        for(unsigned int point=0; point< integration_points.size(); point++)
                        {
                            MathUtils<double>::InvertMatrix(J[point],InvJ,DetJ);

                            double dV= DetJ*integration_points[point].Weight();

                            for(unsigned int prim=0; prim<(*it)->GetGeometry().size() ; prim++)
                            {
                                for(unsigned int sec=0; sec<(*it)->GetGeometry().size() ; sec++)
                                {
                                    M(((*it)->GetGeometry()[prim].Id()-1),
                                         ((*it)->GetGeometry()[sec].Id()-1))+=
                                            Ncontainer(point, prim)*Ncontainer(point, sec)*dV;
                                }
                            }
                        }
                    }

                }

//                 for( ElementsArrayType::ptr_iterator it = ElementsArray.ptr_begin(); 
//                      it != ElementsArray.ptr_end(); ++it )
//                 {
//                     const IntegrationPointsArrayType& integration_points 
//                             = (*it)->GetGeometry().IntegrationPoints((*it)->GetIntegrationMethod());
// 
//                     GeometryType::JacobiansType J(integration_points.size());
//                     J = (*it)->GetGeometry().Jacobian(J, (*it)->GetIntegrationMethod());
// 					
//                     const Matrix& Ncontainer = (*it)->GetGeometry().ShapeFunctionsValues((*it)->GetIntegrationMethod());
// // 					
//                     Matrix InvJ(3,3);
//                     double DetJ;							
//                     for(unsigned int point=0; point< integration_points.size(); point++)
//                     {
//                         MathUtils<double>::InvertMatrix(J[point],InvJ,DetJ);
// 
//                         double dV= DetJ*integration_points[point].Weight();
// 
//                         for(unsigned int prim=0; prim<(*it)->GetGeometry().size() ; prim++)
//                         {
//                             for(unsigned int sec=0; sec<(*it)->GetGeometry().size() ; sec++)
//                             {
//                                 M(((*it)->GetGeometry()[prim].Id()-1),
//                                      ((*it)->GetGeometry()[sec].Id()-1))+=
//                                         Ncontainer(point, prim)*Ncontainer(point, sec)*dV;
//                             }
//                         }
//                     }
//                 }
                
                double stop_prod = omp_get_wtime();
                std::cout << "assembling time: " << stop_prod - start_prod << std::endl;
                

                for(unsigned int firstvalue=0; firstvalue<3; firstvalue++)
                {
                    for(unsigned int secondvalue=0; secondvalue<3; secondvalue++)
                    {
                        noalias(g)= ZeroVector(model_part.NumberOfNodes());

                        noalias(b)= ZeroVector(model_part.NumberOfNodes());
						//Transfer of GaussianVariables to Nodal Variablias via L_2-Minimization 
						// see Jiao + Heath "Common-refinement-based data tranfer ..."
						// International Journal for numerical methods in engineering 61 (2004) 2402--2427 
						// for general description of L_2-Minimization
                        for( ElementsArrayType::ptr_iterator it = ElementsArray.ptr_begin(); 
                             it != ElementsArray.ptr_end();
                             ++it )
                        {
                            const IntegrationPointsArrayType& integration_points 
                                    = (*it)->GetGeometry().IntegrationPoints((*it)->GetIntegrationMethod());

                            GeometryType::JacobiansType J(integration_points.size());
                            J = (*it)->GetGeometry().Jacobian(J, (*it)->GetIntegrationMethod());
                            std::vector<Matrix> ValuesOnIntPoint(integration_points.size());
					
                            (*it)->GetValueOnIntegrationPoints(rThisVariable, ValuesOnIntPoint, model_part.GetProcessInfo());
							
                            const Matrix& Ncontainer = (*it)->GetGeometry().ShapeFunctionsValues((*it)->GetIntegrationMethod());
							
                            Matrix InvJ(3,3);
                            double DetJ;							
                            for(unsigned int point=0; point< integration_points.size(); point++)
                            {
                                MathUtils<double>::InvertMatrix(J[point],InvJ,DetJ);

                                double dV= DetJ*integration_points[point].Weight();

                                for(unsigned int prim=0 ; prim<(*it)->GetGeometry().size(); prim++)
                                {
                                    b(((*it)->GetGeometry()[prim].Id()-1))
                                            +=(ValuesOnIntPoint[point](firstvalue,secondvalue))
                                            *Ncontainer(point, prim)*dV;	
                                }
                            }
                        }
                        double start_solve = omp_get_wtime();
                        SkylineLUFactorizationSolver<SpaceType, SpaceType>().Solve(M, g, b);
                        double stop_solve = omp_get_wtime();
                        std::cout << "solving time: " << stop_solve - start_solve << std::endl;
                        for(ModelPart::NodeIterator it = model_part.NodesBegin() ; 
                            it != model_part.NodesEnd() ; it++)
                        {
                            it->GetSolutionStepValue(rThisVariable)(firstvalue,secondvalue)
                                    = g((it->Id()-1));
                        }
                    }//END firstvalue
                }//END secondvalue
            }
			
			/**
             * Transfer of rThisVariable defined on integration points to corresponding 
             * nodal values. The transformation is done in a form that ensures a minimization
             * of L_2-norm error (/sum{rThisVariable- f(x)) whereas 
             * f(x)= /sum{shape_func_i*rThisVariable_i}
             * @param model_part model_part on which the transfer should be done
             * @param rThisVariable Matrix-Variable which should be transferred
             * @see TransferVariablesToNodes(ModelPart& model_part, Variable<Kratos::Matrix>& rThisVariable)
             * @see TransferVariablesToNodes(ModelPart& model_part, Variable<double>& rThisVariable)
             * @ref Jiao&Heath: "Common-refinement-based data transfer...", Int.
             * Journal for numer. meth. in eng. 61 (2004) 2402--2427
             * WARNING: this may cause segmentation faults as the respective variables
             * will be created on nodal level while they are originally intended to be 
             * stored on integration points!
                         */
            void TransferVariablesToNodes(ModelPart& model_part, Variable<Kratos::Vector>& rThisVariable)
            {
                ElementsArrayType& ElementsArray= model_part.Elements();

				//loop over all master surfaces (global search)
                for(ModelPart::NodeIterator it = model_part.NodesBegin();
                    it != model_part.NodesEnd() ; it++)
                {
                    it->GetSolutionStepValue(rThisVariable)
                            = ZeroVector(6);
                }
				//SetUpEquationSystem
                SpaceType::MatrixType M(model_part.NumberOfNodes(),model_part.NumberOfNodes());
                SpaceType::VectorType g(model_part.NumberOfNodes());
                SpaceType::VectorType b(model_part.NumberOfNodes());
                noalias(M)= ZeroMatrix(model_part.NumberOfNodes(),model_part.NumberOfNodes());
/*
                for( ElementsArrayType::ptr_iterator it = ElementsArray.ptr_begin(); 
                it != ElementsArray.ptr_end(); ++it )
                {
									//SetUpEquationSystem
                SpaceType::MatrixType M((*it)->GetGeometry().size(),(*it)->GetGeometry().size());
                SpaceType::VectorType g((*it)->GetGeometry().size());
                SpaceType::VectorType b((*it)->GetGeometry().size());
                noalias(M)= ZeroMatrix((*it)->GetGeometry().size(),(*it)->GetGeometry().size());
					
                const IntegrationPointsArrayType& integration_points 
                = (*it)->GetGeometry().IntegrationPoints((*it)->GetIntegrationMethod());

                GeometryType::JacobiansType J(integration_points.size());
                J = (*it)->GetGeometry().Jacobian(J, (*it)->GetIntegrationMethod());
					
                std::vector<Vector> ValuesOnIntPoint(integration_points.size());
					
                (*it)->GetValueOnIntegrationPoints(rThisVariable, ValuesOnIntPoint, model_part.GetProcessInfo());
					
                const Matrix& Ncontainer = (*it)->GetGeometry().ShapeFunctionsValues((*it)->GetIntegrationMethod());
					
                Matrix InvJ(3,3);
                double DetJ;							
                for(unsigned int point=0; point< integration_points.size(); point++)
                {
                MathUtils<double>::InvertMatrix(J[point],InvJ,DetJ);

                double dV= DetJ*integration_points[point].Weight();

                for(unsigned int prim=0 ; prim<(*it)->GetGeometry().size(); prim++)
                {
                for(unsigned int sec=0 ; sec<(*it)->GetGeometry().size(); sec++)
                {
                M(prim,sec)+=
                Ncontainer(point, prim)*Ncontainer(point, sec)*dV;
            }
            }
            }
				
                for(unsigned int firstvalue=0; firstvalue<6; firstvalue++)
                {
                noalias(g)= ZeroVector((*it)->GetGeometry().size());

                noalias(b)= ZeroVector((*it)->GetGeometry().size());
						//Transfer of GaussianVariables to Nodal Variablias via L_2-Minimization 
						// see Jiao + Heath "Common-refinement-based data tranfer ..."
						// International Journal for numerical methods in engineering 61 (2004) 2402--2427 
						// for general description of L_2-Minimization

                Matrix InvJ(3,3);
                double DetJ;							
                for(unsigned int point=0; point< integration_points.size(); point++)
                {
                MathUtils<double>::InvertMatrix(J[point],InvJ,DetJ);

                double dV= DetJ*integration_points[point].Weight();

                for(unsigned int prim=0 ; prim<(*it)->GetGeometry().size(); prim++)
                {
                b(prim)
                +=(ValuesOnIntPoint[point](firstvalue))
                *Ncontainer(point, prim)*dV;	
            }
            }
					
                SkylineLUFactorizationSolver<SpaceType, SpaceType>().Solve(M, g, b);
						
                for(unsigned int prim=0 ; prim<(*it)->GetGeometry().size(); prim++)
                {
                (*it)->GetGeometry()[prim].GetSolutionStepValue(rThisVariable)(firstvalue)
                = g(prim);
            }
            }//END firstvalue
            }
				*/
				
                //create a partition of the element array
                int number_of_threads = omp_get_max_threads();
                vector<unsigned int> element_partition;
                CreatePartition(number_of_threads, ElementsArray.size(), element_partition);
                KRATOS_WATCH( number_of_threads );
                KRATOS_WATCH( element_partition );


                double start_prod = omp_get_wtime();
                
                #pragma omp parallel for 
                for(int k=0; k<number_of_threads; k++)
                {
                    ElementsArrayType::ptr_iterator it_begin=ElementsArray.ptr_begin()+element_partition[k];
                    ElementsArrayType::ptr_iterator it_end=ElementsArray.ptr_begin()+element_partition[k+1];
                    
                    for( ElementsArrayType::ptr_iterator it = it_begin; it != it_end; ++it )
                    {
                        const IntegrationPointsArrayType& integration_points 
                                = (*it)->GetGeometry().IntegrationPoints((*it)->GetIntegrationMethod());
                        GeometryType::JacobiansType J(integration_points.size());
                    
                        J = (*it)->GetGeometry().Jacobian(J, (*it)->GetIntegrationMethod());
                        const Matrix& Ncontainer =
                                (*it)->GetGeometry().ShapeFunctionsValues((*it)->GetIntegrationMethod());
                        
                        Matrix InvJ(3,3);
                        double DetJ;							
                        for(unsigned int point=0; point< integration_points.size(); point++)
                        {
                            MathUtils<double>::InvertMatrix(J[point],InvJ,DetJ);
                            double dV= DetJ*integration_points[point].Weight();
                            for(unsigned int prim=0 ; prim<(*it)->GetGeometry().size(); prim++)
                            {
                                for(unsigned int sec=0 ; sec<(*it)->GetGeometry().size(); sec++)
                                {
                                    M(((*it)->GetGeometry()[prim].Id()-1),
                                         ((*it)->GetGeometry()[sec].Id()-1))+=
                                            Ncontainer(point, prim)*Ncontainer(point, sec)*dV;
                                }
                            }
                        }
                    }
                }
                
                double stop_prod = omp_get_wtime();
                std::cout << "assembling time: " << stop_prod - start_prod << std::endl;
                
                for(unsigned int firstvalue=0; firstvalue<6; firstvalue++)
                {
                    noalias(g)= ZeroVector(model_part.NumberOfNodes());

                    noalias(b)= ZeroVector(model_part.NumberOfNodes());
					//Transfer of GaussianVariables to Nodal Variablias via L_2-Minimization 
					// see Jiao + Heath "Common-refinement-based data tranfer ..."
					// International Journal for numerical methods in engineering 61 (2004) 2402--2427 
					// for general description of L_2-Minimization
                    for( ElementsArrayType::ptr_iterator it = ElementsArray.ptr_begin(); 
                         it != ElementsArray.ptr_end();
                         ++it )
                    {
                        const IntegrationPointsArrayType& integration_points 
                                = (*it)->GetGeometry().IntegrationPoints( (*it)->GetIntegrationMethod());

                        GeometryType::JacobiansType J(integration_points.size());
                        J = (*it)->GetGeometry().Jacobian(J, (*it)->GetIntegrationMethod());
                        std::vector<Vector> ValuesOnIntPoint(integration_points.size());
					
                        (*it)->GetValueOnIntegrationPoints(rThisVariable, ValuesOnIntPoint, model_part.GetProcessInfo());
							
                        const Matrix& Ncontainer = (*it)->GetGeometry().ShapeFunctionsValues((*it)->GetIntegrationMethod());
							
                        Matrix InvJ(3,3);
                        double DetJ;							
                        for(unsigned int point=0; point< integration_points.size(); point++)
                        {
                            MathUtils<double>::InvertMatrix(J[point],InvJ,DetJ);

                            double dV= DetJ*integration_points[point].Weight();

                            for(unsigned int prim=0 ; prim<(*it)->GetGeometry().size(); prim++)
                            {
                                b(((*it)->GetGeometry()[prim].Id()-1))
                                        +=(ValuesOnIntPoint[point](firstvalue))
                                        *Ncontainer(point, prim)*dV;	
                            }
                        }
                    }
                    double start_solve = omp_get_wtime();
                    SkylineLUFactorizationSolver<SpaceType, SpaceType>().Solve(M, g, b);
                    double stop_solve = omp_get_wtime();
                    std::cout << "solving time: " << stop_solve - start_solve << std::endl;
                    for(ModelPart::NodeIterator it = model_part.NodesBegin() ; 
                        it != model_part.NodesEnd() ; it++)
                    {
                        it->GetSolutionStepValue(rThisVariable)(firstvalue)
                                = g((it->Id()-1));
                    }
                }//END firstvalue
            }

			/**
             * Transfer of rThisVariable defined on integration points to corresponding 
             * nodal values. The transformation is done in a form that ensures a minimization
             * of L_2-norm error (/sum{rThisVariable- f(x)) whereas 
             * f(x)= /sum{shape_func_i*rThisVariable_i}
             * @param model_part model_part on which the transfer should be done
             * @param rThisVariable Matrix-Variable which should be transferred
             * @see TransferVariablesToNodes(ModelPart& model_part, Variable<Kratos::Matrix>& rThisVariable)
             * @see TransferVariablesToNodes(ModelPart& model_part, Variable<Kratos::Vector>& rThisVariable)
             * @ref Jiao&Heath: "Common-refinement-based data transfer...", Int.
             * Journal for numer. meth. in eng. 61 (2004) 2402--2427
             * WARNING: this may cause segmentation faults as the respective variables
             * will be created on nodal level while they are originally intended to be 
             * stored on integration points!
                         */
            void TransferVariablesToNodes(ModelPart& model_part, Variable<double>& rThisVariable)
            {
                ElementsArrayType& ElementsArray= model_part.Elements();
				//loop over all master surfaces (global search)
                for(ModelPart::NodeIterator it = model_part.NodesBegin();
                    it != model_part.NodesEnd() ; it++)
                {
                    it->GetSolutionStepValue(rThisVariable)
                            = 0.0;
                }
					
				//SetUpEquationSystem	
                SpaceType::MatrixType M(model_part.NumberOfNodes(),model_part.NumberOfNodes());
                noalias(M)= ZeroMatrix(model_part.NumberOfNodes(),model_part.NumberOfNodes());
                SpaceType::VectorType g(model_part.NumberOfNodes());
                noalias(g)= ZeroVector(model_part.NumberOfNodes());
                SpaceType::VectorType b(model_part.NumberOfNodes());
                noalias(b)= ZeroVector(model_part.NumberOfNodes());
				//Transfer of GaussianVariables to Nodal Variablias via L_2-Minimization 
				// see Jiao + Heath "Common-refinement-based data tranfer ..."
				// International Journal for numerical methods in engineering 61 (2004) 2402--2427 
				// for general description of L_2-Minimization
                //create a partition of the element array
                int number_of_threads = omp_get_max_threads();
                vector<unsigned int> element_partition;
                CreatePartition(number_of_threads, ElementsArray.size(), element_partition);
                KRATOS_WATCH( number_of_threads );
                KRATOS_WATCH( element_partition );


                double start_prod = omp_get_wtime();
                
                #pragma omp parallel for 
                for(int k=0; k<number_of_threads; k++)
                {
                    ElementsArrayType::ptr_iterator it_begin=ElementsArray.ptr_begin()+element_partition[k];
                    ElementsArrayType::ptr_iterator it_end=ElementsArray.ptr_begin()+element_partition[k+1];
                    
                    for( ElementsArrayType::ptr_iterator it = it_begin; it != it_end; ++it )
                    {	
                        const IntegrationPointsArrayType& integration_points 
                                = (*it)->GetGeometry().IntegrationPoints((*it)->GetIntegrationMethod());
							
                        GeometryType::JacobiansType J(integration_points.size());
                        J = (*it)->GetGeometry().Jacobian(J, (*it)->GetIntegrationMethod());
						
                        std::vector<double> ValuesOnIntPoint(integration_points.size());
							
                        (*it)->GetValueOnIntegrationPoints(rThisVariable, ValuesOnIntPoint, model_part.GetProcessInfo());
							
                        const Matrix& Ncontainer = (*it)->GetGeometry().ShapeFunctionsValues((*it)->GetIntegrationMethod());

                        Matrix InvJ(3,3);
                        double DetJ;
							
                        for(unsigned int point=0; point< integration_points.size(); point++)
                        {
                            MathUtils<double>::InvertMatrix(J[point],InvJ,DetJ);

                            double dV= DetJ*integration_points[point].Weight();
                            for(unsigned int prim=0 ; prim<(*it)->GetGeometry().size(); prim++)
                            {
                                b(((*it)->GetGeometry()[prim].Id()-1))
                                        +=(ValuesOnIntPoint[point])
                                        *Ncontainer(point, prim)*dV;
                                for(unsigned int sec=0 ; sec<(*it)->GetGeometry().size(); sec++)
                                {
                                    M(((*it)->GetGeometry()[prim].Id()-1),
                                         ((*it)->GetGeometry()[sec].Id()-1))+=
                                            Ncontainer(point, prim)*Ncontainer(point, sec)*dV;
                                }
                            }
                        }
                    }
                }
                
                double stop_prod = omp_get_wtime();
                std::cout << "assembling time: " << stop_prod - start_prod << std::endl;
                
                double start_solve = omp_get_wtime();
                SkylineLUFactorizationSolver<SpaceType, SpaceType>().Solve(M, g, b);
                double stop_solve = omp_get_wtime();
                std::cout << "solving time: " << stop_solve - start_solve << std::endl;
                for(ModelPart::NodeIterator it = model_part.NodesBegin() ; 
                    it != model_part.NodesEnd() ; it++)
                {
                    it->GetSolutionStepValue(rThisVariable)
                            = g((it->Id()-1));
                }
            }
			/**
             * Transfer of rThisVariable stored on nodes form source mesh to target mesh.
             * The transformation is done in a way that ensures a minimization
             * of L_2-norm error (/sum{f_old(x)- f_new(x)) whereas 
             * f(x)_old/new= /sum{shape_func_i*rThisVariable_i}
             * @param rSource source model_part
             * @param rTarget target model_part
             * @param rThisVariable Matrix-Variable which should be transferred
             * @see TransferVariablesBetweenMeshes(ModelPart& rSource, ModelPart& rTarget,
            Variable<Kratos::Vector>& rThisVariable)
             * @see TransferVariablesBetweenMeshes(ModelPart& rSource, ModelPart& rTarget,
            Variable<double>& rThisVariable)
             * @ref Jiao&Heath: "Common-refinement-based data transfer...", Int.
             * Journal for numer. meth. in eng. 61 (2004) 2402--2427
                         */
            void TransferVariablesBetweenMeshes(ModelPart& rSource, ModelPart& rTarget,
                    Variable<Kratos::Matrix>& rThisVariable)
            {
                ElementsArrayType& SourceMeshElementsArray= rSource.Elements();

                ElementsArrayType& TargetMeshElementsArray= rTarget.Elements();
				//loop over all master surfaces (global search)
                for(ModelPart::NodeIterator it = rTarget.NodesBegin();
                    it != rTarget.NodesEnd() ; it++)
                {
                    it->GetSolutionStepValue(rThisVariable)
                            = ZeroMatrix(3,3);
                }
				//SetUpEquationSystem
                SpaceType::MatrixType M(rTarget.NumberOfNodes(),rTarget.NumberOfNodes());
                noalias(M)= ZeroMatrix(rTarget.NumberOfNodes(),rTarget.NumberOfNodes());
                SpaceType::VectorType g(rTarget.NumberOfNodes());
                SpaceType::VectorType b(rTarget.NumberOfNodes());
                for( ElementsArrayType::ptr_iterator it = TargetMeshElementsArray.ptr_begin(); 
                     it != TargetMeshElementsArray.ptr_end();
                     ++it )
                {
                    const IntegrationPointsArrayType& integration_points 
                            = (*it)->GetGeometry().IntegrationPoints((*it)->GetIntegrationMethod());

                    GeometryType::JacobiansType J(integration_points.size());
                    J = (*it)->GetGeometry().Jacobian(J, (*it)->GetIntegrationMethod());
							
                    const Matrix& Ncontainer = (*it)->GetGeometry().ShapeFunctionsValues((*it)->GetIntegrationMethod());

                    Matrix InvJ(3,3);
                    double DetJ;
                    for(unsigned int point=0; point< integration_points.size(); point++)
                    {
                        MathUtils<double>::InvertMatrix(J[point],InvJ,DetJ);

                        double dV= DetJ*integration_points[point].Weight();

                        for(unsigned int prim=0; prim<(*it)->GetGeometry().size() ; prim++)
                        {
                            for(unsigned int sec=0; sec<(*it)->GetGeometry().size() ; sec++)
                            {
                                M(((*it)->GetGeometry()[prim].Id()-1), ((*it)->GetGeometry()[sec].Id()-1))+=
                                        Ncontainer(point, prim)*Ncontainer(point, sec)*dV;
                            }
                        }
                    }
                }

                for(unsigned int firstvalue= 0; firstvalue< 3; firstvalue++)
                {
                    for(unsigned int secondvalue= 0; secondvalue< 3; secondvalue++)
                    {
                        noalias(b)= ZeroVector(rTarget.NumberOfNodes());
                        noalias(g)= ZeroVector(rTarget.NumberOfNodes());
						//Transfer of GaussianVariables to Nodal Variablias via L_2-Minimization 
						// see Jiao + Heath "Common-refinement-based data tranfer ..."
						// International Journal for numerical methods in engineering 61 (2004) 2402--2427 
						// for general description of L_2-Minimization
                        for( ElementsArrayType::ptr_iterator it = TargetMeshElementsArray.ptr_begin(); 
                             it != TargetMeshElementsArray.ptr_end();
                             ++it )
                        {
                            const IntegrationPointsArrayType& integration_points 
                                    = (*it)->GetGeometry().IntegrationPoints((*it)->GetIntegrationMethod());

                            GeometryType::JacobiansType J(integration_points.size());
                            J = (*it)->GetGeometry().Jacobian(J, (*it)->GetIntegrationMethod());
							
                            const Matrix& Ncontainer = (*it)->GetGeometry().ShapeFunctionsValues((*it)->GetIntegrationMethod());

                            Matrix InvJ(3,3);
                            double DetJ;

                            for(unsigned int point=0; point< integration_points.size(); point++)
                            {
                                MathUtils<double>::InvertMatrix(J[point],InvJ,DetJ);

                                Point<3> sourceLocalPoint;
                                Point<3> targetLocalPoint;
                                noalias(targetLocalPoint)= integration_points[point];
                                Point<3> targetGlobalPoint;
                                (*it)->GetGeometry().GlobalCoordinates(targetGlobalPoint,
                                  targetLocalPoint);
                                Element sourceElement;
                                double functionValue;
								//Calculate Value of rVariable(firstvalue, secondvalue) in OldMesh
                                if(FindPartnerElement(targetGlobalPoint, SourceMeshElementsArray,
                                   sourceElement,sourceLocalPoint))
                                {

                                    functionValue=
                                            ValueMatrixInOldMesh( sourceElement,sourceLocalPoint,rThisVariable, firstvalue, secondvalue );
                                }
                                else
                                {
                                    std::cout<<"###### NO PARTNER FOUND IN OLD MESH : TransferVariablesBetweenMeshes(...Matrix...)#####"<<std::endl;
                                    continue;
                                }

                                double dV= DetJ*integration_points[point].Weight();

                                for(unsigned int prim=0; prim<(*it)->GetGeometry().size(); prim++)
                                {
                                    b(((*it)->GetGeometry()[prim].Id()-1))
                                            +=functionValue
                                            *Ncontainer(point, prim)*dV;
                                }
                            }
                        }

                        SkylineLUFactorizationSolver<SpaceType, SpaceType>().Solve(M, g, b);
                        for(ModelPart::NodeIterator it = rTarget.NodesBegin() ; 
                            it != rTarget.NodesEnd() ; it++)
                        {
                            it->GetSolutionStepValue(rThisVariable)(firstvalue,secondvalue)
                                    = g((it->Id()-1));
                        }
                    }//END firstvalue
                }//END secondvalue
            }
			/**
             * Transfer of rThisVariable stored on nodes form source mesh to target mesh.
             * The transformation is done in a way that ensures a minimization
             * of L_2-norm error (/sum{f_old(x)- f_new(x)) whereas 
             * f(x)_old/new= /sum{shape_func_i*rThisVariable_i}
             * @param rSource source model_part
             * @param rTarget target model_part
             * @param rThisVariable Vector-Variable which should be transferred
             * @see TransferVariablesBetweenMeshes(ModelPart& rSource, ModelPart& rTarget,
            Variable<Kratos::Matrix>& rThisVariable)
             * @see TransferVariablesBetweenMeshes(ModelPart& rSource, ModelPart& rTargetw,
            Variable<double>& rThisVariable)
             * @ref Jiao&Heath: "Common-refinement-based data transfer...", Int.
             * Journal for numer. meth. in eng. 61 (2004) 2402--2427
                         */
            void TransferVariablesBetweenMeshes(ModelPart& rSource, ModelPart& rTarget,
                    Variable<Kratos::Vector>& rThisVariable)
            {
                ElementsArrayType& SourceMeshElementsArray= rSource.Elements();

                ElementsArrayType& TargetMeshElementsArray= rTarget.Elements();
				//loop over all master surfaces (global search)
                for(ModelPart::NodeIterator it = rTarget.NodesBegin();
                    it != rTarget.NodesEnd() ; it++)
                {
                    it->GetSolutionStepValue(rThisVariable)
                            = ZeroVector(6);
                }
				//SetUpEquationSystem
                SpaceType::MatrixType M(rTarget.NumberOfNodes(),rTarget.NumberOfNodes());
                noalias(M)= ZeroMatrix(rTarget.NumberOfNodes(),rTarget.NumberOfNodes());
                SpaceType::VectorType g(rTarget.NumberOfNodes());
                SpaceType::VectorType b(rTarget.NumberOfNodes());

                for( ElementsArrayType::ptr_iterator it = TargetMeshElementsArray.ptr_begin(); 
                     it != TargetMeshElementsArray.ptr_end();
                     ++it )
                {
                    const IntegrationPointsArrayType& integration_points 
                            = (*it)->GetGeometry().IntegrationPoints((*it)->GetIntegrationMethod());

                    GeometryType::JacobiansType J(integration_points.size());
                    J = (*it)->GetGeometry().Jacobian(J, (*it)->GetIntegrationMethod());
							
                    const Matrix& Ncontainer = (*it)->GetGeometry().ShapeFunctionsValues((*it)->GetIntegrationMethod());

                    Matrix InvJ(3,3);
                    double DetJ;

                    for(unsigned int point=0; point< integration_points.size(); point++)
                    {
                        MathUtils<double>::InvertMatrix(J[point],InvJ,DetJ);

                        double dV= DetJ*integration_points[point].Weight();

                        for(unsigned int prim=0 ;prim<(*it)->GetGeometry().size() ;prim++)
                        {
                            for(unsigned int sec=0 ;sec<(*it)->GetGeometry().size() ;sec++)
                            {
                                M(((*it)->GetGeometry()[prim].Id()-1), ((*it)->GetGeometry()[sec].Id()-1))+=
                                        Ncontainer(point, prim)*Ncontainer(point, sec)*dV;
                            }
                        }
                    }
                }
                for(unsigned int firstvalue= 0; firstvalue< 6; firstvalue++)
                {
                    noalias(b)= ZeroVector(rTarget.NumberOfNodes());
                    noalias(g)= ZeroVector(rTarget.NumberOfNodes());
					//Transfer of GaussianVariables to Nodal Variablias via L_2-Minimization 
					// see Jiao + Heath "Common-refinement-based data tranfer ..."
					// International Journal for numerical methods in engineering 61 (2004) 2402--2427 
					// for general description of L_2-Minimization
                    for( ElementsArrayType::ptr_iterator it = TargetMeshElementsArray.ptr_begin(); 
                         it != TargetMeshElementsArray.ptr_end();
                         ++it )
                    {
                        const IntegrationPointsArrayType& integration_points 
                                = (*it)->GetGeometry().IntegrationPoints((*it)->GetIntegrationMethod());

                        GeometryType::JacobiansType J(integration_points.size());
                        J = (*it)->GetGeometry().Jacobian(J, (*it)->GetIntegrationMethod());
							
                        const Matrix& Ncontainer = (*it)->GetGeometry().ShapeFunctionsValues((*it)->GetIntegrationMethod());

                        Matrix InvJ(3,3);
                        double DetJ;

                        for(unsigned int point=0; point< integration_points.size(); point++)
                        {
                            MathUtils<double>::InvertMatrix(J[point],InvJ,DetJ);

                            Point<3> sourceLocalPoint;
                            Point<3> targetLocalPoint;
                            noalias(targetLocalPoint)= integration_points[point];
                            Point<3> targetGlobalPoint;
                            (*it)->GetGeometry().GlobalCoordinates(targetGlobalPoint,
                              targetLocalPoint);
                            Element sourceElement;
                            double functionValue;
							//Calculate Value of rVariable(firstvalue, secondvalue) in OldMesh
                            if(FindPartnerElement(targetGlobalPoint, SourceMeshElementsArray,
                               sourceElement,sourceLocalPoint))
                            {
                                functionValue=
                                        ValueVectorInOldMesh( sourceElement,sourceLocalPoint,rThisVariable, firstvalue);
                            }
                            else
                            {
                                std::cout<<"###### NO PARTNER FOUND IN OLD MESH : TransferVariablesBetweenMeshes(...Vector...)#####"<<std::endl;
                                continue;
                            }

                            double dV= DetJ*integration_points[point].Weight();

                            for(unsigned int prim=0 ;prim<(*it)->GetGeometry().size() ;prim++)
                            {
                                b(((*it)->GetGeometry()[prim].Id()-1))
                                        +=functionValue
                                        *Ncontainer(point, prim)*dV;
                            }
                        }
                    }
                    SkylineLUFactorizationSolver<SpaceType, SpaceType>().Solve(M, g, b);
                    for(ModelPart::NodeIterator it = rTarget.NodesBegin() ; 
                        it != rTarget.NodesEnd() ; it++)
                    {
                        it->GetSolutionStepValue(rThisVariable)(firstvalue)
                                = g((it->Id()-1));
                    }
                }//END firstvalue
            }
			/**
             * Transfer of rThisVariable stored on nodes form source mesh to target mesh.
             * The transformation is done in a way that ensures a minimization
             * of L_2-norm error (/sum{f_old(x)- f_new(x)) whereas 
             * f(x)_old/new= /sum{shape_func_i*rThisVariable_i}
             * @param rSource source model_part
             * @param rTarget target model_part
             * @param rThisVariable double-Variable which should be transferred
             * @see TransferVariablesBetweenMeshes(ModelPart& rSource, ModelPart& rTarget,
            Variable<Kratos::Matrix>& rThisVariable)
             * @see TransferVariablesBetweenMeshes(ModelPart& rSource, ModelPart& rTarget,
            Variable<Kratos::Vector>& rThisVariable)
             * @ref Jiao&Heath: "Common-refinement-based data transfer...", Int.
             * Journal for numer. meth. in eng. 61 (2004) 2402--2427
                         */
            void TransferVariablesBetweenMeshes(ModelPart& rSource, ModelPart& rTarget,
                    Variable<double>& rThisVariable)
            {
                ElementsArrayType& SourceMeshElementsArray= rSource.Elements();

                ElementsArrayType& TargetMeshElementsArray= rTarget.Elements();
				//loop over all master surfaces (global search)
                for(ModelPart::NodeIterator it = rTarget.NodesBegin();
                    it != rTarget.NodesEnd() ; it++)
                {
                    it->GetSolutionStepValue(rThisVariable)
                            = 0.0;
                }
				//SetUpEquationSystem
                SpaceType::MatrixType M(rTarget.NumberOfNodes(),rTarget.NumberOfNodes());
                noalias(M)= ZeroMatrix(rTarget.NumberOfNodes(),rTarget.NumberOfNodes());
                SpaceType::VectorType g(rTarget.NumberOfNodes());
                noalias(g)= ZeroVector(rTarget.NumberOfNodes());
                SpaceType::VectorType b(rTarget.NumberOfNodes());
                noalias(b)= ZeroVector(rTarget.NumberOfNodes());
				//Transfer of GaussianVariables to Nodal Variablias via L_2-Minimization 
				// see Jiao + Heath "Common-refinement-based data tranfer ..."
				// International Journal for numerical methods in engineering 61 (2004) 2402--2427 
				// for general description of L_2-Minimization
                for( ElementsArrayType::ptr_iterator it = TargetMeshElementsArray.ptr_begin(); 
                     it != TargetMeshElementsArray.ptr_end();
                     ++it )
                {
                    const IntegrationPointsArrayType& integration_points 
                            = (*it)->GetGeometry().IntegrationPoints((*it)->GetIntegrationMethod());

                    GeometryType::JacobiansType J(integration_points.size());
                    J = (*it)->GetGeometry().Jacobian(J, (*it)->GetIntegrationMethod());
							
                    const Matrix& Ncontainer = (*it)->GetGeometry().ShapeFunctionsValues((*it)->GetIntegrationMethod());

                    Matrix InvJ(3,3);
                    double DetJ;

                    for(unsigned int point=0; point< integration_points.size(); point++)
                    {
                        MathUtils<double>::InvertMatrix(J[point],InvJ,DetJ);

                        Point<3> sourceLocalPoint;
                        Point<3> targetLocalPoint;
                        noalias(targetLocalPoint)= integration_points[point];
                        Point<3> targetGlobalPoint;
                        (*it)->GetGeometry().GlobalCoordinates(targetGlobalPoint,
                          targetLocalPoint);
                        Element sourceElement;
                        double functionValue;
						//Calculate Value of rVariable(firstvalue, secondvalue) in OldMesh
                        if(FindPartnerElement(targetGlobalPoint, SourceMeshElementsArray,
                           sourceElement,sourceLocalPoint))
                        {
                            functionValue=
                                    MappedValue( sourceElement,sourceLocalPoint,rThisVariable);
                        }
                        else
                        {
                            std::cout<<"###### NO PARTNER FOUND IN OLD MESH : TransferVariablesBetweenMeshes(...double...)#####"<<std::endl;
                            continue;
                        }

                        double dV= DetJ*integration_points[point].Weight();

                        for(unsigned int prim=0 ;prim<(*it)->GetGeometry().size() ;prim++)
                        {
                            b(((*it)->GetGeometry()[prim].Id()-1))
                                    +=functionValue*Ncontainer(point, prim)*dV;
                            for(unsigned int sec=0; sec<(*it)->GetGeometry().size(); sec++)
                            {
                                M(((*it)->GetGeometry()[prim].Id()-1), ((*it)->GetGeometry()[sec].Id()-1))+=
                                        Ncontainer(point, prim)*Ncontainer(point, sec)*dV;
                            }
                        }
                    }
                }
                SkylineLUFactorizationSolver<SpaceType, SpaceType>().Solve(M, g, b);
                for(ModelPart::NodeIterator it = rTarget.NodesBegin() ; 
                    it != rTarget.NodesEnd() ; it++)
                {
                    it->GetSolutionStepValue(rThisVariable)
                            = g((it->Id()-1));
                }
            }
			/**
             * Auxiliary function.
             * This one calculates the target value of given Variable by shape-function-based
             * interpolation of the nodal values from given source element to the given
             * target point that is assumed to lie within the source element
             * @return value of given variable in new point
             * @param oldElement corresponding element in source mesh
             * @param localPoint given target point to map the variable to
             * @param rThisVariable given variable to be transferred
             * @see ValueVectorInOldMesh(Element& oldElement, Point<3>&  localPoint, 
            const Variable<Kratos::Vector>& rThisVariable )
             * @see MappedValue( Element& sourceElement, Point<3>& targetPoint,
            const Variable<double>& rThisVariable)
                         */
            Matrix ValueMatrixInOldMesh(Element& oldElement, Point<3>&  localPoint, 
                                        const Variable<Kratos::Matrix>& rThisVariable )
            {
                Matrix newValue(3,3);
                noalias(newValue)= ZeroMatrix(3,3);
                Matrix temp(3,3);

                for(unsigned int i=0; i< oldElement.GetGeometry().size(); i++)
                {
                    noalias(temp)= oldElement.GetGeometry()[i].GetSolutionStepValue(rThisVariable);

                    for(unsigned int k=0; k<3; k++)
                        for(unsigned int l=0; l<3; l++)
                            newValue(k,l)
                                    +=oldElement.GetGeometry().ShapeFunctionValue( i, localPoint)
                                    *temp(k,l);
                }

                return newValue;
            }
			
			/**
             * Auxiliary function.
             * This one calculates the target value of given Variable by shape-function-based
             * interpolation of the nodal values from given source element to the given
             * target point that is assumed to lie within the source element
             * @return value of given variable in new point
             * @param oldElement corresponding element in source mesh
             * @param localPoint given target point to map the variable to
             * @param rThisVariable given variable to be transferred
             * @see ValueMatrixInOldMesh(Element& oldElement, Point<3>&  localPoint, 
            const Variable<Kratos::Matrix>& rThisVariable )
             * @see MappedValue( Element& sourceElement, Point<3>& targetPoint,
            const Variable<double>& rThisVariable)
                         */
            Vector ValueVectorInOldMesh(Element& oldElement, Point<3>&  localPoint, 
                                        const Variable<Kratos::Vector>& rThisVariable )
            {
                Vector newValue(6);
                noalias(newValue) = ZeroVector(6);
                Vector temp(6);	

                for(unsigned int i=0; i<oldElement.GetGeometry().size(); i++)
                {
                    noalias(temp)= oldElement.GetGeometry()[i].GetSolutionStepValue(rThisVariable);
                    for(unsigned int k=0; k<6; k++)
                        newValue(k) +=oldElement.GetGeometry().ShapeFunctionValue( i, localPoint)*temp(k);
                }
                return newValue;
            }

			/**
             * Auxiliary function.
             * This one calculates the target value of given Variable by shape-function-based
             * interpolation of the nodal values from given source element to the given
             * target point that is assumed to lie within the source element
             * @return value of given variable in new point
             * @param sourceElement corresponding element in source mesh
             * @param targetPoint given target point to map the variable to
             * @param rThisVariable given variable to be transferred
             * @see ValueMatrixInOldMesh(Element& oldElement, Point<3>&  localPoint, 
            const Variable<Kratos::Matrix>& rThisVariable )
             * @see ValueVectorInOldMesh(Element& oldElement, Point<3>&  localPoint, 
            const Variable<Kratos::Vector>& rThisVariable )
                         */
            double MappedValuePressure( Element& sourceElement, Point<3>& targetPoint,
                                        const Variable<double>& rThisVariable)
            {
                double newValue = 0.0;

                Geometry<Node<3> >::Pointer pPressureGeometry;

                if(sourceElement.GetGeometry().size()==20 || sourceElement.GetGeometry().size()==27)
                    pPressureGeometry= Geometry<Node<3> >::Pointer(new Hexahedra3D8 <Node<3> >(
                        sourceElement.GetGeometry()(0),sourceElement.GetGeometry()(1),
                        sourceElement.GetGeometry()(2),sourceElement.GetGeometry()(3),
                                sourceElement.GetGeometry()(4),sourceElement.GetGeometry()(5),
                                        sourceElement.GetGeometry()(6),sourceElement.GetGeometry()(7)));

                if(sourceElement.GetGeometry().size()==10 )
                    pPressureGeometry= Geometry<Node<3> >::Pointer(new Tetrahedra3D4 <Node<3> >(
                        sourceElement.GetGeometry()(0),sourceElement.GetGeometry()(1),
                        sourceElement.GetGeometry()(2),sourceElement.GetGeometry()(3)));

                for(unsigned int i= 0; i< pPressureGeometry->size(); i++)
                {
                    newValue+=
                            pPressureGeometry->ShapeFunctionValue( i, targetPoint)
                            *sourceElement.GetGeometry()[i].GetSolutionStepValue(rThisVariable);
                }
                return newValue;
            }

			/**
             * Auxiliary function.
             * This one calculates the target value of given Variable by shape-function-based
             * interpolation of the nodal values from given source element to the given
             * target point that is assumed to lie within the source element
             * @return value of given variable in new point
             * @param sourceElement corresponding element in source mesh
             * @param targetPoint given target point to map the variable to
             * @param rThisVariable given variable to be transferred
             * @see ValueMatrixInOldMesh(Element& oldElement, Point<3>&  localPoint, 
            const Variable<Kratos::Matrix>& rThisVariable )
             * @see ValueVectorInOldMesh(Element& oldElement, Point<3>&  localPoint, 
            const Variable<Kratos::Vector>& rThisVariable )
                         */
            double MappedValue( Element& sourceElement, Point<3>& targetPoint,
                                const Variable<double>& rThisVariable)
            {
                double newValue = 0.0;

                for(unsigned int i= 0; i< sourceElement.GetGeometry().size(); i++)
                {
                    newValue+=
                            sourceElement.GetGeometry().ShapeFunctionValue( i, targetPoint)
                            *sourceElement.GetGeometry()[i].GetSolutionStepValue(rThisVariable);
                }
                return newValue;
            }
			
			/**
             * Auxiliary function.
             * This one calculates the target value of given Variable by shape-function-based
             * interpolation of the nodal values from given source element to the given
             * target point that is assumed to lie within the source element
             * @return value of given variable in new point
             * @param sourceElement corresponding element in source mesh
             * @param targetPoint given target point to map the variable to
             * @param rThisVariable given variable to be transferred
                         */
            Vector MappedValue( Element& sourceElement, Point<3>& targetPoint,
                                const Variable<array_1d<double, 3 > >& rThisVariable)
            {
                Vector newValue = ZeroVector(3);
				
                for(unsigned int i=0; i<sourceElement.GetGeometry().size(); i++)
                {
                    newValue+=
                            sourceElement.GetGeometry().ShapeFunctionValue( i, targetPoint)
                            *sourceElement.GetGeometry()[i].GetSolutionStepValue(rThisVariable);
                }
                return newValue;
            }
            
            /**
             * calculates for a point given with the physical coords newNode
             * the element oldElement where it lays in and the natural coords
             * localPoint within this element
             * @return whether a corresponding element and natural coords could be found
             * @param newNode physical coordinates of given point
             * @param OldMeshElementsArray Array of elements wherein the search should be performed
             * @param oldElement corresponding element for newNode
             * @param rResult corresponding natural coords for newNode
             * TODO: find a faster method for outside search (hextree? etc.), maybe outside this
             * function by restriction of OldMeshElementsArray
             */
            bool FindPartnerElement( CoordinatesArrayType& newNode, 
                                     const ElementsArrayType& OldMeshElementsArray, 
                                     Element& oldElement, Point<3>&  rResult)
            {
                bool partner_found= false;
				//noalias(rResult)= ZeroVector(3);
                ElementsArrayType::Pointer OldElementsSet( new ElementsArrayType() );
                std::vector<double > OldMinDist;
                bool newMinDistFound= false;

                int counter = 0;
                do
                {
                    double minDist = 1.0e120;
                    newMinDistFound= false;
                    OldElementsSet->clear();
					//loop over all master surfaces (global search)
                    for( ElementsArrayType::ptr_const_iterator it = OldMeshElementsArray.ptr_begin(); 
                         it != OldMeshElementsArray.ptr_end(); ++it )
                    {
						//loop over all nodes in tested element
                        for( unsigned int n=0; n<(*it)->GetGeometry().size(); n++ )
                        {
                            double dist = ((*it)->GetGeometry().GetPoint(n).X0()-newNode[0])
                                        *((*it)->GetGeometry().GetPoint(n).X0()-newNode[0])
                                        +((*it)->GetGeometry().GetPoint(n).Y0()-newNode[1])
                                        *((*it)->GetGeometry().GetPoint(n).Y0()-newNode[1])
                                        +((*it)->GetGeometry().GetPoint(n).Z0()-newNode[2])
                                        *((*it)->GetGeometry().GetPoint(n).Z0()-newNode[2]);
                            if( fabs(dist-minDist) < 1e-7 )
                            {
                                OldElementsSet->push_back(*it);
                            }
                            else if( dist < minDist )
                            {
                                bool alreadyUsed= false;
                                for(unsigned int old_dist= 0; old_dist<OldMinDist.size(); old_dist++)  
                                {
                                    if(fabs(dist- OldMinDist[old_dist])< 1e-7 )
                                        alreadyUsed= true;
                                }
                                if(!alreadyUsed)
                                {
                                    OldElementsSet->clear();
                                    minDist = dist;
                                    OldElementsSet->push_back(*it);
                                    newMinDistFound= true;
                                }
                            }
                        }
                    }

                    OldMinDist.push_back(minDist);
//                     KRATOS_WATCH(OldElementsSet->size());

                    for( ElementsArrayType::ptr_const_iterator it = OldElementsSet->ptr_begin(); 
                         it != OldElementsSet->ptr_end(); ++it )
                    {
//                         std::cout << "checking elements list" << std::endl;
                        if( (*it)->GetGeometry().IsInside( newNode, rResult ) )
                        {
//                             std::cout << "isInside" << std::endl;
                            oldElement=**(it);
                            partner_found=true;
                            return partner_found;
                        }
                    }
                    std::cout << counter << std::endl;
                    counter++;
                    if( counter > 27 )
                        break;
                }while(newMinDistFound);

                if(!partner_found)
                    std::cout<<" !!!! NO PARTNER FOUND !!!! "<<std::endl;
                return partner_found;
            }
			
			//***************************************************************************
			//***************************************************************************
			/**
             * Auxiliary function.
             * This one calculates the target value of given Matrix-Variable at row firtsvalue
             * and column secondvalue by shape-function-based
             * interpolation of the nodal values from given source element to the given
             * target point that is assumed to lie within the source element
             * @return value of given variable in new point
             * @param sourceElement corresponding element in source mesh
             * @param targetPoint given target point to map the variable to
             * @param rThisVariable given variable to be transferred
             * @param firstvalue row index
             * @param secondvalue column index
             * @see ValueVectorInOldMesh(Element& oldElement, Point<3>&  localPoint, 
            const Variable<Kratos::Vector>& rThisVariable, unsigned int firstvalue)
                         */
            double ValueMatrixInOldMesh(Element& oldElement, Point<3>&  localPoint, 
                                        const Variable<Kratos::Matrix>& rThisVariable, unsigned int firstvalue, unsigned int secondvalue )
            {
                double newValue= 0.0;

                for(unsigned int i=0; i<oldElement.GetGeometry().size(); i++)
                {
                    newValue
                            +=oldElement.GetGeometry().ShapeFunctionValue( i, localPoint)
                            *oldElement.GetGeometry()[i].GetSolutionStepValue(rThisVariable)(firstvalue,secondvalue);
                }
                return newValue;
            }
			/**
             * Auxiliary function.
             * This one calculates the target value of given Vector-Variable at firtsvalue
             * by shape-function-based
             * interpolation of the nodal values from given source element to the given
             * target point that is assumed to lie within the source element
             * @return value of given variable in new point
             * @param sourceElement corresponding element in source mesh
             * @param targetPoint given target point to map the variable to
             * @param rThisVariable given variable to be transferred
             * @param firstvalue index
             * @see ValueVectorInOldMesh(Element& oldElement, Point<3>&  localPoint, 
            const Variable<Kratos::Vector>& rThisVariable, unsigned int firstvalue)
                         */
            double ValueVectorInOldMesh(Element& oldElement, Point<3>&  localPoint, 
                                        const Variable<Kratos::Vector>& rThisVariable, unsigned int firstvalue )
            {
                double newValue= 0.0;
				
                for(unsigned int i=0; i<oldElement.GetGeometry().size(); i++)
                {
                    newValue
                            +=oldElement.GetGeometry().ShapeFunctionValue( i, localPoint)
                            *oldElement.GetGeometry()[i].GetSolutionStepValue(rThisVariable)(firstvalue);
                }
                return newValue;
            }
            
            inline void CreatePartition(unsigned int number_of_threads,const int number_of_rows, vector<unsigned int>& partitions)
            {
                partitions.resize(number_of_threads+1);
                int partition_size = number_of_rows / number_of_threads;
                partitions[0] = 0;
                partitions[number_of_threads] = number_of_rows;
                for(unsigned int i = 1; i<number_of_threads; i++)
                    partitions[i] = partitions[i-1] + partition_size ;
            }

    };//Class Scheme
}//namespace Kratos.

#endif /* KRATOS_PARALLEL_VARIABLE_TRANSFER_UTILITY_INCLUDED  defined */
