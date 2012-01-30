/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).
Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:
Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.
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
//   Last Modified by:    $Author: jcotela $
//   Date:                $Date: 2010-10-09 10:34:00 $
//   Revision:            $Revision: 0.1 $
//
//
#if !defined(KRATOS_TWO_FLUID_TwoFluidVMS_H_INCLUDED )
#define  KRATOS_TWO_FLUID_TwoFluidVMS_H_INCLUDED
// System includes
#include <string>
#include <iostream>
// External includes
// Project includes
#include "containers/array_1d.h"
#include "includes/define.h"
#include "custom_elements/vms.h"
#include "includes/serializer.h"
#include "utilities/geometry_utilities.h"
#include "utilities/split_tetrahedra.h"
#include "utilities/enrichment_utilities.h"
// Application includes
#include "fluid_dynamics_application_variables.h"
#include "vms.h"
namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{
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
template< unsigned int TDim,
          unsigned int TNumNodes = TDim + 1 >
class TwoFluidVMS : public VMS<TDim, TNumNodes>
{
	public:
		///@name Type Definitions
		///@{
		/// Pointer definition of TwoFluidVMS
		KRATOS_CLASS_POINTER_DEFINITION(TwoFluidVMS);
		///base type: an IndexedObject that automatically has a unique number
		typedef IndexedObject BaseType;
		///Element from which it is derived
		typedef VMS<TDim, TNumNodes> ElementBaseType;
		///definition of node type (default is: Node<3>)
		typedef Node < 3 > NodeType;
		/**
		 * Properties are used to store any parameters
		 * related to the constitutive law
		 */
		typedef Properties PropertiesType;
		///definition of the geometry type with given NodeType
		typedef Geometry<NodeType> GeometryType;
		///definition of nodes container type, redefined from GeometryType
		typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
		typedef Vector VectorType;
		typedef typename ElementBaseType::MatrixType MatrixType;
		typedef std::size_t IndexType;
		typedef std::size_t SizeType;
		typedef std::vector<std::size_t> EquationIdVectorType;
		typedef std::vector< Dof<double>::Pointer > DofsVectorType;
		typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;
		typedef VectorMap<IndexType, DataValueContainer> SolutionStepsElementalDataContainerType;
		///@}
		///@name Life Cycle
		///@{
		//Constructors.
		/// Default constuctor.
		/**
		 * @param NewId Index number of the new element (optional)
		 */
		TwoFluidVMS(IndexType NewId = 0) :
			ElementBaseType(NewId)
		{
		}
		///Constructor using an array of nodes.
		/**
		 * @param NewId Index of the new element
		 * @param ThisNodes An array containing the nodes of the new element
		 */
		TwoFluidVMS(IndexType NewId, const NodesArrayType& ThisNodes) :
			ElementBaseType(NewId, ThisNodes)
		{
		}
		/// Constructor using a geometry object.
		/**
		 * @param NewId Index of the new element
		 * @param pGeometry Pointer to a geometry object
		 */
		TwoFluidVMS(IndexType NewId, GeometryType::Pointer pGeometry) :
			ElementBaseType(NewId, pGeometry)
		{
		}
		/// Constuctor using geometry and properties.
		/**
		 * @param NewId Index of the new element
		 * @param pGeometry Pointer to a geometry object
		 * @param pProperties Pointer to the element's properties
		 */
		TwoFluidVMS(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) :
			ElementBaseType(NewId, pGeometry, pProperties)
		{
		}
		/// Destructor.
		virtual ~TwoFluidVMS()
		{
		}
		///@}
		///@name Operators
		///@{
		///@}
		///@name Operations
		///@{
		/// Create a new element of this type
		/**
		 * Returns a pointer to a new TwoFluidVMS element, created using given input
		 * @param NewId: the ID of the new element
		 * @param ThisNodes: the nodes of the new element
		 * @param pProperties: the properties assigned to the new element
		 * @return a Pointer to the new element
		 */
		Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
		                        PropertiesType::Pointer pProperties) const
		{
			return Element::Pointer(new TwoFluidVMS(NewId, (this->GetGeometry()).Create(ThisNodes), pProperties));
		}
		/// Provides local contributions from body forces and projections to the RHS
		/**
		 * This is called during the assembly process and provides the RHS terms of the
		 * system that are either constant or computed explicitly (from the 'old'
		 * iteration variables). In this case this means the body force terms and the
		 * OSS projections, that are treated explicitly.
		 * @param rRightHandSideVector Will be filled with the elemental right hand side
		 * @param rCurrentProcessInfo ProcessInfo instance from the ModelPart. It is
		 * expected to contain values for OSS_SWITCH, DYNAMIC_TAU and DELTA_TIME
		 */
		virtual void CalculateRightHandSide(VectorType& rRightHandSideVector,
		                                    ProcessInfo& rCurrentProcessInfo)
		{
			const unsigned int LocalSize = (TDim + 1) * TNumNodes;
			// Check sizes and initialize
			if (rRightHandSideVector.size() != LocalSize)
				rRightHandSideVector.resize(LocalSize, false);
			noalias(rRightHandSideVector) = ZeroVector(LocalSize);
			// Calculate this element's geometric parameters
			double Area;
			array_1d<double, TNumNodes> N;
			boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim> DN_DX;
			GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);
			//get position of the cut surface
			Vector distances(TNumNodes);
			Matrix Nenriched(6, 1);
			Vector volumes(6);
			Matrix coords(TNumNodes, TDim);
			Matrix Ngauss(6, TNumNodes);
			Vector signs(6);
			std::vector< Matrix > gauss_gradients(6);
			//fill coordinates
			for (unsigned int i = 0; i < TNumNodes; i++)
			{
				const array_1d<double, 3 > & xyz = this->GetGeometry()[i].Coordinates();
				volumes[i] = 0.0;
				distances[i] = this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
				for (unsigned int j = 0; j < TDim; j++)
					coords(i, j) = xyz[j];
			}
			for (unsigned int i = 0; i < 6; i++)
				gauss_gradients[i].resize(1, TDim, false);
			unsigned int ndivisions = EnrichmentUtilities::CalculateTetrahedraEnrichedShapeFuncions(coords, DN_DX, distances, volumes, Ngauss, signs, gauss_gradients, Nenriched);
			//do integration
			for (unsigned int igauss = 0; igauss < ndivisions; igauss++)
			{
				//assigning the gauss data
				for (unsigned int k = 0; k < TNumNodes; k++)
					N[k] = Ngauss(igauss, k);
				double wGauss = volumes[igauss];
				// Calculate this element's fluid properties
				double Density;
				this->EvaluateInPoint(Density, DENSITY, N);
				// Calculate Momentum RHS contribution
				this->AddMomentumRHS(rRightHandSideVector, Density, N, wGauss);
				// For OSS: Add projection of residuals to RHS
				if (rCurrentProcessInfo[OSS_SWITCH] == 1)
				{
					array_1d<double, 3 > AdvVel;
					this->GetAdvectiveVel(AdvVel, N);
					double KinViscosity;
					this->EvaluateInPoint(KinViscosity, VISCOSITY, N);
					double Viscosity;
					this->GetEffectiveViscosity(Density, KinViscosity, N, DN_DX, Viscosity, rCurrentProcessInfo);
					// Calculate stabilization parameters
					double TauOne, TauTwo;
//                    if (ndivisions == 1)
					this->CalculateTau(TauOne, TauTwo, AdvVel, Area, Density, Viscosity, rCurrentProcessInfo);
//                else
//                {
//                    TauOne = 0.0;
//                    TauTwo = 0.0;
//                }
//                    this->CalculateTau(TauOne, TauTwo, AdvVel, Area, Viscosity, rCurrentProcessInfo);
					this->AddProjectionToRHS(rRightHandSideVector, AdvVel, Density, TauOne, TauTwo, N, DN_DX, wGauss, rCurrentProcessInfo[DELTA_TIME]);
				}
			}
		}
		/// Computes local contributions to the mass matrix
		/**
		 * Provides the local contributions to the mass matrix, which is defined here
		 * as the matrix associated to velocity derivatives. Note that the mass
		 * matrix implemented here is lumped.
		 * @param rMassMatrix Will be filled with the elemental mass matrix
		 * @param rCurrentProcessInfo the current process info instance
		 */
		virtual void MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
		{
			const unsigned int LocalSize = (TDim + 1) * TNumNodes;
			// Resize and set to zero
			if (rMassMatrix.size1() != LocalSize)
				rMassMatrix.resize(LocalSize, LocalSize, false);
			rMassMatrix = ZeroMatrix(LocalSize, LocalSize);
			// Get the element's geometric parameters
			double Area;
			array_1d<double, TNumNodes> N;
			boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim> DN_DX;
			GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);
			//get position of the cut surface
			Vector distances(TNumNodes);
			Matrix Nenriched(6, 1);
			Vector volumes(6);
			Matrix coords(TNumNodes, TDim);
			Matrix Ngauss(6, TNumNodes);
			Vector signs(6);
			std::vector< Matrix > gauss_gradients(6);
			//fill coordinates
			for (unsigned int i = 0; i < TNumNodes; i++)
			{
				const array_1d<double, 3 > & xyz = this->GetGeometry()[i].Coordinates();
				volumes[i] = 0.0;
				distances[i] = this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
				for (unsigned int j = 0; j < TDim; j++)
					coords(i, j) = xyz[j];
			}
			for (unsigned int i = 0; i < 6; i++)
				gauss_gradients[i].resize(1, TDim, false);
			unsigned int ndivisions = EnrichmentUtilities::CalculateTetrahedraEnrichedShapeFuncions(coords, DN_DX, distances, volumes, Ngauss, signs, gauss_gradients, Nenriched);
			//mass matrix
			for (unsigned int igauss = 0; igauss < ndivisions; igauss++)
			{
				//assigning the gauss data
				for (unsigned int k = 0; k < TNumNodes; k++)
					N[k] = Ngauss(igauss, k);
				double wGauss = volumes[igauss];
				// Calculate this element's fluid properties
				double Density;
				this->EvaluateInPoint(Density, DENSITY, N);
				// Consisten Mass Matrix
				this->AddConsistentMassMatrixContribution(rMassMatrix, N, Density, wGauss);
			}
			this->LumpMassMatrix(rMassMatrix);
//if (ndivisions > 1)
//{
//KRATOS_WATCH(this->Id());
//KRATOS_WATCH(rMassMatrix);
//
//for (unsigned int igauss = 0; igauss < ndivisions; igauss++)
//{
//    double dist = 0.0;
//    for(unsigned int k=0; k<4; k++)
//    {
//        dist += Ngauss(igauss, k)*this->GetGeometry()[k].FastGetSolutionStepValue(DISTANCE);
//    }
//    KRATOS_WATCH(dist);
//    KRATOS_WATCH(signs[igauss]);
//    if( signs[igauss] * dist < 0.0 )
//        KRATOS_ERROR(std::logic_error,"sign of partition does not coincide","")
//}
//
//}

			//stabilization terms
			for (unsigned int igauss = 0; igauss < ndivisions; igauss++)
			{
				//assigning the gauss data
				for (unsigned int k = 0; k < TNumNodes; k++)
					N[k] = Ngauss(igauss, k);
				double wGauss = volumes[igauss];
				// Calculate this element's fluid properties
				double Density;
				this->EvaluateInPoint(Density, DENSITY, N);
				/* For ASGS: add dynamic stabilization terms.
				 * These terms are not used in OSS, as they belong to the finite element
				 * space and cancel out with their projections.
				 */
				if (rCurrentProcessInfo[OSS_SWITCH] != 1)
				{
					double KinViscosity;
					this->EvaluateInPoint(KinViscosity, VISCOSITY, N);
					double Viscosity;
					this->GetEffectiveViscosity(Density, KinViscosity, N, DN_DX, Viscosity, rCurrentProcessInfo);
					// Get Advective velocity
					array_1d<double, 3 > AdvVel;
					this->GetAdvectiveVel(AdvVel, N);
					// Calculate stabilization parameters
					double TauOne, TauTwo;
					this->CalculateTau(TauOne, TauTwo, AdvVel, Area, Density, Viscosity, rCurrentProcessInfo);
//                    if (ndivisions == 1)
//                       this->CalculateTau(TauOne, TauTwo, AdvVel, Area, Viscosity, rCurrentProcessInfo);
//                    else
//                    {
//                        TauOne = 0.0;
//                        TauTwo = 0.0;
//                    }
					// Add dynamic stabilization terms ( all terms involving a delta(u) )
					this->AddMassStabTerms(rMassMatrix, Density, AdvVel, TauOne, N, DN_DX, wGauss);
				}
			}
		}
		/// Computes the local contribution associated to 'new' velocity and pressure values
		/**
		 * Provides local contributions to the system associated to the velocity and
		 * pressure terms (convection, diffusion, pressure gradient/velocity divergence
		 * and stabilization).
		 * @param rDampMatrix Will be filled with the velocity-proportional "damping" matrix
		 * @param rRightHandSideVector the elemental right hand side vector
		 * @param rCurrentProcessInfo the current process info instance
		 */
		virtual void CalculateLocalVelocityContribution(MatrixType& rDampMatrix,
		        VectorType& rRightHandSideVector,
		        ProcessInfo& rCurrentProcessInfo)
		{
			const unsigned int LocalSize = (TDim + 1) * TNumNodes;
			// Resize and set to zero the matrix
			// Note that we don't clean the RHS because it will already contain body force (and stabilization) contributions
			if (rDampMatrix.size1() != LocalSize)
				rDampMatrix.resize(LocalSize, LocalSize, false);
			noalias(rDampMatrix) = ZeroMatrix(LocalSize, LocalSize);
			// Get this element's geometric properties
			double Area;
			array_1d<double, TNumNodes> N;
			boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim> DN_DX;
			GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);
			//get position of the cut surface
			Vector distances(TNumNodes);
			Matrix Nenriched(6, 1);
			Vector volumes(6);
			Matrix coords(TNumNodes, TDim);
			Matrix Ngauss(6, TNumNodes);
			Vector signs(6);
			std::vector< Matrix > gauss_gradients(6);
			//fill coordinates
			for (unsigned int i = 0; i < TNumNodes; i++)
			{
				const array_1d<double, 3 > & xyz = this->GetGeometry()[i].Coordinates();
				volumes[i] = 0.0;
				distances[i] = this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
				for (unsigned int j = 0; j < TDim; j++)
					coords(i, j) = xyz[j];
			}
			for (unsigned int i = 0; i < 6; i++)
				gauss_gradients[i].resize(1, TDim, false);
			unsigned int ndivisions = EnrichmentUtilities::CalculateTetrahedraEnrichedShapeFuncions(coords, DN_DX, distances, volumes, Ngauss, signs, gauss_gradients, Nenriched);
			Vector enrichment_terms_vertical = ZeroVector(LocalSize);
			Vector enrichment_terms_horizontal = ZeroVector(LocalSize);
			double enrichment_diagonal = 0.0;
			double enriched_rhs = 0.0;
			array_1d<double,3> bf(3,0.0);
//             array_1d<double,3> Acceleration(3,0.0);
			double positive_volume = 0.0;
			double negative_volume = 0.0;
			//do integration
			for (unsigned int igauss = 0; igauss < ndivisions; igauss++)
			{
				//assigning the gauss data
				for (unsigned int k = 0; k < TNumNodes; k++)
					N[k] = Ngauss(igauss, k);
				double wGauss = volumes[igauss];
				if(signs[igauss] > 0)
					positive_volume += wGauss;
				else
					negative_volume += wGauss;
				// Calculate this element's fluid properties
				double Density, KinViscosity;
				this->EvaluateInPoint(Density, DENSITY, N);
				this->EvaluateInPoint(KinViscosity, VISCOSITY, N);
				this->EvaluateInPoint(bf,BODY_FORCE,N);
//                 this->EvaluateInPoint(Acceleration,ACCELERATION,N);
// if(Density != 1.0 && Density != 1000.0)
//   KRATOS_ERROR(std::logic_error,"taking the average of density on element ",this->Id());
				double Viscosity;
				this->GetEffectiveViscosity(Density, KinViscosity, N, DN_DX, Viscosity, rCurrentProcessInfo);
				// Get Advective velocity
				array_1d<double, 3 > AdvVel;
				this->GetAdvectiveVel(AdvVel, N);
				// Calculate stabilization parameters
				double TauOne, TauTwo;
//                if (ndivisions == 1)
				this->CalculateTau(TauOne, TauTwo, AdvVel, Area, Density, Viscosity, rCurrentProcessInfo);
//                else
//                {
//                    TauOne = rCurrentProcessInfo[DELTA_TIME];
//                    TauTwo = 0.0;
//                }
				double Dt =  rCurrentProcessInfo[DELTA_TIME];
				this->AddIntegrationPointVelocityContribution(rDampMatrix, rRightHandSideVector, Density, Viscosity, AdvVel, TauOne, TauTwo, N, DN_DX, wGauss);
				if (ndivisions > 1)
				{
//                    TauOne = rCurrentProcessInfo[DELTA_TIME];
					//compute enrichment terms contribution
					for (unsigned int inode = 0; inode < TNumNodes; inode++)
					{
						int base_index = (TDim + 1) * inode;
						array_1d<double,TNumNodes> AGradN(TNumNodes,0.0);
						this->GetConvectionOperator(AGradN,AdvVel,DN_DX);
						//momentum term
						for (unsigned int k = 0; k < TDim; k++)
						{
							double ConvTerm = wGauss * TauOne * gauss_gradients[igauss](0,k)* Density * AGradN[inode];
							enrichment_terms_vertical[base_index + k] += ConvTerm - wGauss * DN_DX(inode, k) * Nenriched(igauss, 0);
							enrichment_terms_horizontal[base_index + k] += ConvTerm + wGauss * DN_DX(inode, k) * Nenriched(igauss, 0);
//                             enrichment_terms_vertical[base_index + k] +=wGauss*N[inode]*gauss_gradients[igauss](0, k); //-= wGauss * DN_DX(inode, k) * Nenriched(igauss, 0);
//                            enrichment_terms_horizontal[base_index + k] -=Density*wGauss*N[inode]*gauss_gradients[igauss](0, k); //   += Density*wGauss * DN_DX(inode, k) * Nenriched(igauss, 0);
						}
						//pressure term
						for (unsigned int k = 0; k < TDim; k++)
						{
							double temp =  wGauss * TauOne* DN_DX(inode, k) * gauss_gradients[igauss](0, k);
							enrichment_terms_vertical[base_index + TDim] += temp;
							enrichment_terms_horizontal[base_index + TDim] += temp;
						}
						//add acceleration enrichment term
						const array_1d<double,3>& old_vnode = this->GetGeometry()[inode].FastGetSolutionStepValue(VELOCITY,1);
						for (unsigned int k = 0; k < TDim; k++)
						{
							double coeff = wGauss * TauOne *Density *  gauss_gradients[igauss](0,k)*Nenriched(igauss, 0) / Dt;
							enrichment_terms_horizontal[base_index + k] += coeff;
							enriched_rhs += coeff * old_vnode[k];
//                             enrichment_terms_vertical[base_index + k] +=wGauss*N[inode]*gauss_gradients[igauss](0, k); //-= wGauss * DN_DX(inode, k) * Nenriched(igauss, 0);
//                            enrichment_terms_horizontal[base_index + k] -=Density*wGauss*N[inode]*gauss_gradients[igauss](0, k); //   += Density*wGauss * DN_DX(inode, k) * Nenriched(igauss, 0);
						}
					}
					//compute diagonal enrichment term
					for (unsigned int k = 0; k < TDim; k++)
					{
						const Matrix& enriched_grad = gauss_gradients[igauss];
						enrichment_diagonal += wGauss  * TauOne * pow(enriched_grad(0, k), 2);
						enriched_rhs += wGauss * TauOne *Density * enriched_grad(0,k)*(bf[k]/*-Acceleration[k]*/);
					}
				}
			}
// if(this->Id() == 2059)
// {
//   KRATOS_WATCH("********************************************")
//   KRATOS_WATCH(distances)
//   KRATOS_WATCH(ndivisions)
//   KRATOS_WATCH(Ngauss);
//   KRATOS_WATCH(volumes);
//   for (unsigned int k = 0; k < ndivisions; k++)
//     KRATOS_WATCH(gauss_gradients[k])
//
// }
			/*	    for (unsigned int i = 0; i < (TDim+1)*TNumNodes; i++)
				    {
				      double b = rRightHandSideVector[i];
				      if ( b == rRightHandSideVector[i] + 10000000.0)
				      {
					for (unsigned int k = 0; k < TNumNodes; k++)
					{
			                    N[k] = Ngauss(0, k);
					    KRATOS_WATCH(this->GetGeometry()[k].FastGetSolutionStepValue(DISTANCE));
					}
			                double Density, KinViscosity;
			                this->EvaluateInPoint(Density, DENSITY, N);
			                this->EvaluateInPoint(KinViscosity, VISCOSITY, N);
			 		KRATOS_WATCH(N);
					KRATOS_WATCH(this->Id())
					KRATOS_WATCH(ndivisions);
					KRATOS_WATCH(Area);
					KRATOS_WATCH(enrichment_diagonal)
					KRATOS_WATCH(b)
					KRATOS_WATCH(rRightHandSideVector)
					KRATOS_WATCH(volumes)
					KRATOS_WATCH(rDampMatrix)
					KRATOS_WATCH(Ngauss)
					KRATOS_WATCH(Density)
					KRATOS_WATCH(DN_DX)
					KRATOS_WATCH(KinViscosity)
				      }
				    }	 */
			if (ndivisions > 1)
			{
//                KRATOS_WATCH(enrichment_terms_vertical);
//                KRATOS_WATCH(enrichment_terms_horizontal);
//                KRATOS_WATCH(enrichment_diagonal);
				/*		double aux_diag = 0.0;
						for (unsigned int inode = 0; inode < TNumNodes; inode++)
						  aux_diag += rDampMatrix((TDim + 1) * inode + TDim, (TDim + 1) * inode + TDim );*/
				//add to LHS enrichment contributions
				double inverse_diag_term = 1.0 / ( enrichment_diagonal);
// 		KRATOS_WATCH(inverse_diag_term);
// 		KRATOS_WATCH(enrichment_terms_vertical);
// 		KRATOS_WATCH(enrichment_terms_horizontal);
//
// 		for (unsigned int igauss = 0; igauss < ndivisions; igauss++)
// 		{
// 		    //assigning the gauss data
// 		    for (unsigned int k = 0; k < TNumNodes; k++)
// 			N[k] = Ngauss(igauss, k);
// 		    double wGauss = volumes[igauss];
//
// 		    KRATOS_WATCH(wGauss);
// 		}


/*		if(fabs(negative_volume) < Area*1e-4 || fabs(positive_volume) < Area*1e-4)
		{
 		  KRATOS_WATCH(this->Id());
			KRATOS_WATCH(negative_volume);
 		  KRATOS_WATCH(positive_volume);
 		}*/

// 		if(negative_volume > Area*1e-6 && negative_volume < Area*0.999999 ) //only add enrichment if the splitted volume is not too small
// 		{
				for (unsigned int i = 0; i < LocalSize; i++)
					for (unsigned int j = 0; j < LocalSize; j++)
						rDampMatrix(i, j) -= inverse_diag_term * enrichment_terms_vertical[i] * enrichment_terms_horizontal[j];
				//                VectorType U = ZeroVector(LocalSize);
				//                int LocalIndex = 0;
				rRightHandSideVector -= (inverse_diag_term*enriched_rhs )*enrichment_terms_vertical;
//  		}
//  		else
// 		{
// 		  KRATOS_WATCH(negative_volume);
// 		  KRATOS_WATCH(Area);
// 		  KRATOS_WATCH(positive_volume);
// 		}
			}
//             else
// 	    {
// 	       for(unsigned int i=1; i<4; i++)
// 	       {
//
// 		 if(this->GetGeometry()[i].FastGetSolutionStepValue(DENSITY) != this->GetGeometry()[0].FastGetSolutionStepValue(DENSITY))
// 		 {
// 		   KRATOS_WATCH(this->Id());
// 		   KRATOS_WATCH(i)
// 		   KRATOS_WATCH(this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE));
// 		   KRATOS_WATCH(this->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE));
//
// 		   KRATOS_WATCH(this->GetGeometry()[i].FastGetSolutionStepValue(DENSITY));
// 		   KRATOS_WATCH(this->GetGeometry()[0].FastGetSolutionStepValue(DENSITY));
// 		   KRATOS_ERROR(std::logic_error,"element not divided but with two distinct densities","");
// 		 }
// 	       }
// 	    }
			// Now calculate an additional contribution to the residual: r -= rDampMatrix * (u,p)
			VectorType U = ZeroVector(LocalSize);
			int LocalIndex = 0;
			for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
			{
				array_1d< double, 3 > & rVel = this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY);
				for (unsigned int d = 0; d < TDim; ++d) // Velocity Dofs
				{
					U[LocalIndex] = rVel[d];
					++LocalIndex;
				}
				U[LocalIndex] = this->GetGeometry()[iNode].FastGetSolutionStepValue(PRESSURE); // Pressure Dof
				++LocalIndex;
			}
			noalias(rRightHandSideVector) -= prod(rDampMatrix, U);
			this->SetValue(FLAG_VARIABLE,ndivisions);
//            if (ndivisions > 1)
//            {
//                KRATOS_WATCH(this->Id());
//                KRATOS_WATCH(this->GetGeometry());
//                for(unsigned int k=0; k<4; k++)
//                    std::cout << " " << this->GetGeometry()[k].FastGetSolutionStepValue(DISTANCE) ;
//                std::cout << std::endl;
//                KRATOS_WATCH(ndivisions)
//                KRATOS_WATCH(rRightHandSideVector);
//            }
			/*
			            //this is just a check! to be removed
			            Vector aaa(16,0.0);
			            Vector bbb(16,0.0);
			            if (ndivisions > 1)
			            {
			                KRATOS_WATCH(rRightHandSideVector);
			                double inverse_diag_term = 1.0 / enrichment_diagonal;
			                double pstar = inverse_diag_term * (enriched_rhs - inner_prod(enrichment_terms_horizontal, U));
			                //compute grad_p
			                array_1d<double, 4 > pressures;
			                array_1d<double, 3 > grad_p;
			                pressures[0] = this->GetGeometry()[0].FastGetSolutionStepValue(PRESSURE);
			                pressures[1] = this->GetGeometry()[1].FastGetSolutionStepValue(PRESSURE);
			                pressures[2] = this->GetGeometry()[2].FastGetSolutionStepValue(PRESSURE);
			                pressures[3] = this->GetGeometry()[3].FastGetSolutionStepValue(PRESSURE);
			                noalias(grad_p) = prod(trans(DN_DX), pressures);
			                if(this->Id()==2052)
			                {
			                    KRATOS_WATCH("2052")
			                    KRATOS_WATCH(pstar);
			                    KRATOS_WATCH(inner_prod(enrichment_terms_horizontal, U))
			                    KRATOS_WATCH(grad_p);
			                }
			                for (unsigned int igauss = 0; igauss < ndivisions; igauss++)
			                {
			                    //assigning the gauss data
			                    for (unsigned int k = 0; k < TNumNodes; k++)
			                        N[k] = Ngauss(igauss, k);
			                    double wGauss = volumes[igauss];
			                    // Calculate this element's fluid properties
			                    double Density, KinViscosity;
			                    this->EvaluateInPoint(Density, DENSITY, N);
			                    this->EvaluateInPoint(KinViscosity, VISCOSITY, N);
			                    this->EvaluateInPoint(bf, BODY_FORCE, N);
			                    double Viscosity;
			                    this->GetEffectiveViscosity(Density, KinViscosity, N, DN_DX, Viscosity, rCurrentProcessInfo);
			                    array_1d<double, 3 > AdvVel;
			                    this->GetAdvectiveVel(AdvVel, N);
			                    // Calculate stabilization parameters
			                    double TauOne, TauTwo;
			                    this->CalculateTau(TauOne, TauTwo, AdvVel, Area, Viscosity, rCurrentProcessInfo);
			                    for(unsigned int k=0; k<4; k++)
			                    {
			                        for(unsigned int t=0; t<3; t++)
			                        {
			                            aaa[k*(TDim+1)+t] += wGauss*N[k]*(Density*bf[t] - grad_p[t]);
			                            bbb[k*(TDim+1)+t] += wGauss*N[k]*(Density*bf[t]);
			                            aaa[k*(TDim+1)+3] += TauOne*wGauss*DN_DX(k,t)*( Density*bf[t] - grad_p[t]);
			                        }
			                    }
			                    if(this->Id() == 2052)
			                    {
			                        KRATOS_WATCH(igauss)
			                        KRATOS_WATCH(Density * bf);
			                        KRATOS_WATCH(pstar*row(gauss_gradients[igauss], 0))
			                    }
			                }
			                for(unsigned int k=0; k<4; k++)
			                    {
			                        for(unsigned int t=0; t<3; t++)
			                        {
			                             bbb[k*(TDim+1)+t] -= Area*0.25*(grad_p[t]);
			                        }
			                    }
			                KRATOS_WATCH(aaa );
			                KRATOS_WATCH(aaa-bbb );
			                KRATOS_WATCH(pstar*enrichment_terms_vertical);
			                KRATOS_WATCH(aaa - pstar*enrichment_terms_vertical);
			            }
			*/
		}
		/// Implementation of Calculate to compute the local OSS projections.
		/**
		 * If rVariable == ADVPROJ, This function computes the OSS projection
		 * terms using pressure and velocity values from the previous iteration. The
		 * projections are then added to the nodal variables ADVPROJ (Momentum residual)
		 * and DIVPROJ (Mass continuity residual). It is assumed that the scheme will
		 * divide the result by the assembled NODAL_AREA, which is equivalent to a
		 * nodal interpolation using a lumped mass matrix.
		 * @param rVariable Use ADVPROJ
		 * @param Output Will be overwritten with the elemental momentum error
		 * @param rCurrentProcessInfo Process info instance (unused)
		 */
		virtual void Calculate(const Variable<array_1d<double, 3 > >& rVariable,
		                       array_1d<double, 3 > & rOutput,
		                       const ProcessInfo& rCurrentProcessInfo)
		{
			if (rVariable == ADVPROJ) // Compute residual projections for OSS
			{
				// Get the element's geometric parameters
				double Area;
				array_1d<double, TNumNodes> N;
				boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim> DN_DX;
				GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);
				array_1d< double, 3 > ElementalMomRes(3, 0.0);
				double ElementalMassRes(0);
				//get position of the cut surface
				Vector distances(TNumNodes);
				Matrix Nenriched(6, 1);
				Vector volumes(6);
				Matrix coords(TNumNodes, TDim);
				Matrix Ngauss(6, TNumNodes);
				Vector signs(6);
				std::vector< Matrix > gauss_gradients(6);
				//fill coordinates
				for (unsigned int i = 0; i < TNumNodes; i++)
				{
					const array_1d<double, 3 > & xyz = this->GetGeometry()[i].Coordinates();
					volumes[i] = 0.0;
					distances[i] = this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
					for (unsigned int j = 0; j < TDim; j++)
						coords(i, j) = xyz[j];
				}
				for (unsigned int i = 0; i < 6; i++)
					gauss_gradients[i].resize(1, TDim, false);
				unsigned int ndivisions = EnrichmentUtilities::CalculateTetrahedraEnrichedShapeFuncions(coords, DN_DX, distances, volumes, Ngauss, signs, gauss_gradients, Nenriched);
				//do integration
				for (unsigned int igauss = 0; igauss < ndivisions; igauss++)
				{
					//assigning the gauss data
					for (unsigned int k = 0; k < TNumNodes; k++)
						N[k] = Ngauss(igauss, k);
					double wGauss = volumes[igauss];
					// Calculate this element's fluid properties
					double Density, KinViscosity;
					this->EvaluateInPoint(Density, DENSITY, N);
					this->EvaluateInPoint(KinViscosity, VISCOSITY, N);
					double Viscosity;
					this->GetEffectiveViscosity(Density, KinViscosity, N, DN_DX, Viscosity, rCurrentProcessInfo);
					// Get Advective velocity
					array_1d<double, 3 > AdvVel;
					this->GetAdvectiveVel(AdvVel, N);
					// Output containers
					ElementalMomRes = ZeroVector(3);
					ElementalMassRes = 0.0;
					this->AddProjectionResidualContribution(AdvVel, Density, ElementalMomRes, ElementalMassRes, N, DN_DX, wGauss);
					if (rCurrentProcessInfo[OSS_SWITCH] == 1)
					{
						// Carefully write results to nodal variables, to avoid parallelism problems
						for (unsigned int i = 0; i < TNumNodes; ++i)
						{
							this->GetGeometry()[i].SetLock(); // So it is safe to write in the node in OpenMP
							array_1d< double, 3 > & rAdvProj = this->GetGeometry()[i].FastGetSolutionStepValue(ADVPROJ);
							for (unsigned int d = 0; d < TDim; ++d)
								rAdvProj[d] += N[i] * ElementalMomRes[d];
							this->GetGeometry()[i].FastGetSolutionStepValue(DIVPROJ) += N[i] * ElementalMassRes;
							this->GetGeometry()[i].FastGetSolutionStepValue(NODAL_AREA) += wGauss * N[i];
							this->GetGeometry()[i].UnSetLock(); // Free the node for other threads
						}
					}
				}
				/// Return output
				rOutput = ElementalMomRes;
			}
			else if (rVariable == SUBSCALE_VELOCITY)
			{
				// Get the element's geometric parameters
				double Area;
				array_1d<double, TNumNodes> N;
				boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim> DN_DX;
				GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);
				array_1d< double, 3 > ElementalMomRes(3, 0.0);
				double ElementalMassRes(0);
				//get position of the cut surface
				Vector distances(TNumNodes);
				Matrix Nenriched(6, 1);
				Vector volumes(6);
				Matrix coords(TNumNodes, TDim);
				Matrix Ngauss(6, TNumNodes);
				Vector signs(6);
				std::vector< Matrix > gauss_gradients(6);
				//fill coordinates
				for (unsigned int i = 0; i < TNumNodes; i++)
				{
					const array_1d<double, 3 > & xyz = this->GetGeometry()[i].Coordinates();
					volumes[i] = 0.0;
					distances[i] = this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
					for (unsigned int j = 0; j < TDim; j++)
						coords(i, j) = xyz[j];
				}
				for (unsigned int i = 0; i < 6; i++)
					gauss_gradients[i].resize(1, TDim, false);
				unsigned int ndivisions = EnrichmentUtilities::CalculateTetrahedraEnrichedShapeFuncions(coords, DN_DX, distances, volumes, Ngauss, signs, gauss_gradients, Nenriched);
				//do integration
				for (unsigned int igauss = 0; igauss < ndivisions; igauss++)
				{
					//assigning the gauss data
					for (unsigned int k = 0; k < TNumNodes; k++)
						N[k] = Ngauss(igauss, k);
					double wGauss = volumes[igauss];
					// Calculate this element's fluid properties
					double Density, KinViscosity;
					this->EvaluateInPoint(Density, DENSITY, N);
					this->EvaluateInPoint(KinViscosity, VISCOSITY, N);
					double Viscosity;
					this->GetEffectiveViscosity(Density, KinViscosity, N, DN_DX, Viscosity, rCurrentProcessInfo);
					// Get Advective velocity
					array_1d<double, 3 > AdvVel;
					this->GetAdvectiveVel(AdvVel, N);
					// Output containers
					ElementalMomRes = ZeroVector(3);
					ElementalMassRes = 0.0;
					this->AddProjectionResidualContribution(AdvVel, Density, ElementalMomRes, ElementalMassRes, N, DN_DX, wGauss);
					if (rCurrentProcessInfo[OSS_SWITCH] == 1)
					{
						/* Projections of the elemental residual are computed with
						 * Newton-Raphson iterations of type M(lumped) dx = ElemRes - M(consistent) * x
						 */
						const double Weight = ElementBaseType::ConsistentMassCoef(wGauss); // Consistent mass matrix is Weigth * ( Ones(TNumNodes,TNumNodes) + Identity(TNumNodes,TNumNodes) )
						// Carefully write results to nodal variables, to avoid parallelism problems
						for (unsigned int i = 0; i < TNumNodes; ++i)
						{
							this->GetGeometry()[i].SetLock(); // So it is safe to write in the node in OpenMP
							// Add elemental residual to RHS
							array_1d< double, 3 > & rMomRHS = this->GetGeometry()[i].GetValue(ADVPROJ);
							double& rMassRHS = this->GetGeometry()[i].GetValue(DIVPROJ);
							for (unsigned int d = 0; d < TDim; ++d)
								rMomRHS[d] += N[i] * ElementalMomRes[d];
							rMassRHS += N[i] * ElementalMassRes;
							// Write nodal area
							this->GetGeometry()[i].FastGetSolutionStepValue(NODAL_AREA) += wGauss * N[i];
							// Substract M(consistent)*x(i-1) from RHS
							for (unsigned int j = 0; j < TNumNodes; ++j) // RHS -= Weigth * Ones(TNumNodes,TNumNodes) * x(i-1)
							{
								for (unsigned int d = 0; d < TDim; ++d)
									rMomRHS[d] -= Weight * this->GetGeometry()[j].FastGetSolutionStepValue(ADVPROJ)[d];
								rMassRHS -= Weight * this->GetGeometry()[j].FastGetSolutionStepValue(DIVPROJ);
							}
							for (unsigned int d = 0; d < TDim; ++d) // RHS -= Weigth * Identity(TNumNodes,TNumNodes) * x(i-1)
								rMomRHS[d] -= Weight * this->GetGeometry()[i].FastGetSolutionStepValue(ADVPROJ)[d];
							rMassRHS -= Weight * this->GetGeometry()[i].FastGetSolutionStepValue(DIVPROJ);
							this->GetGeometry()[i].UnSetLock(); // Free the node for other threads
						}
					}
				}
				/// Return output
				rOutput = ElementalMomRes;
			}
		}
		/// Checks the input and that all required Kratos variables have been registered.
		/**
		 * This function provides the place to perform checks on the completeness of the input.
		 * It is designed to be called only once (or anyway, not often) typically at the beginning
		 * of the calculations, so to verify that nothing is missing from the input
		 * or that no common error is found.
		 * @param rCurrentProcessInfo The ProcessInfo of the ModelPart that contains this element.
		 * @return 0 if no errors were found.
		 */
		virtual int Check(const ProcessInfo& rCurrentProcessInfo)
		{
			KRATOS_TRY
			// Perform basic element checks
			int ErrorCode = Kratos::Element::Check(rCurrentProcessInfo);
			if (ErrorCode != 0) return ErrorCode;
			// Check that all required variables have been registered
			if (DISTANCE.Key() == 0)
				KRATOS_ERROR(std::invalid_argument, "DISTANCE Key is 0. Check if the application was correctly registered.", "");
			if (VELOCITY.Key() == 0)
				KRATOS_ERROR(std::invalid_argument, "VELOCITY Key is 0. Check if the application was correctly registered.", "");
			if (MESH_VELOCITY.Key() == 0)
				KRATOS_ERROR(std::invalid_argument, "MESH_VELOCITY Key is 0. Check if the application was correctly registered.", "");
			if (ACCELERATION.Key() == 0)
				KRATOS_ERROR(std::invalid_argument, "ACCELERATION Key is 0. Check if the application was correctly registered.", "");
			if (PRESSURE.Key() == 0)
				KRATOS_ERROR(std::invalid_argument, "PRESSURE Key is 0. Check if the application was correctly registered.", "");
			if (DENSITY.Key() == 0)
				KRATOS_ERROR(std::invalid_argument, "DENSITY Key is 0. Check if the application was correctly registered.", "");
			if (VISCOSITY.Key() == 0)
				KRATOS_ERROR(std::invalid_argument, "VISCOSITY Key is 0. Check if the application was correctly registered.", "");
			if (OSS_SWITCH.Key() == 0)
				KRATOS_ERROR(std::invalid_argument, "OSS_SWITCH Key is 0. Check if the application was correctly registered.", "");
			if (DYNAMIC_TAU.Key() == 0)
				KRATOS_ERROR(std::invalid_argument, "DYNAMIC_TAU Key is 0. Check if the application was correctly registered.", "");
			if (DELTA_TIME.Key() == 0)
				KRATOS_ERROR(std::invalid_argument, "DELTA_TIME Key is 0. Check if the application was correctly registered.", "");
			if (ADVPROJ.Key() == 0)
				KRATOS_ERROR(std::invalid_argument, "ADVPROJ Key is 0. Check if the application was correctly registered.", "");
			if (DIVPROJ.Key() == 0)
				KRATOS_ERROR(std::invalid_argument, "DIVPROJ Key is 0. Check if the application was correctly registered.", "");
			if (NODAL_AREA.Key() == 0)
				KRATOS_ERROR(std::invalid_argument, "NODAL_AREA Key is 0. Check if the application was correctly registered.", "");
			if (C_SMAGORINSKY.Key() == 0)
				KRATOS_ERROR(std::invalid_argument, "C_SMAGORINSKY Key is 0. Check if the application was correctly registered.", "");
			if (ERROR_RATIO.Key() == 0)
				KRATOS_ERROR(std::invalid_argument, "ERROR_RATIO Key is 0. Check if the application was correctly registered.", "");
			// Additional variables, only required to print results:
			// SUBSCALE_VELOCITY, SUBSCALE_PRESSURE, TAUONE, TAUTWO, MU, VORTICITY.
			// Checks on nodes
			// Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
			for (unsigned int i = 0; i<this->GetGeometry().size(); ++i)
			{
				if (this->GetGeometry()[i].SolutionStepsDataHas(DISTANCE) == false)
					KRATOS_ERROR(std::invalid_argument, "missing DISTANCE variable on solution step data for node ", this->GetGeometry()[i].Id());
				if (this->GetGeometry()[i].SolutionStepsDataHas(VELOCITY) == false)
					KRATOS_ERROR(std::invalid_argument, "missing VELOCITY variable on solution step data for node ", this->GetGeometry()[i].Id());
				if (this->GetGeometry()[i].SolutionStepsDataHas(PRESSURE) == false)
					KRATOS_ERROR(std::invalid_argument, "missing PRESSURE variable on solution step data for node ", this->GetGeometry()[i].Id());
				if (this->GetGeometry()[i].SolutionStepsDataHas(MESH_VELOCITY) == false)
					KRATOS_ERROR(std::invalid_argument, "missing MESH_VELOCITY variable on solution step data for node ", this->GetGeometry()[i].Id());
				if (this->GetGeometry()[i].SolutionStepsDataHas(ACCELERATION) == false)
					KRATOS_ERROR(std::invalid_argument, "missing ACCELERATION variable on solution step data for node ", this->GetGeometry()[i].Id());
				if (this->GetGeometry()[i].HasDofFor(VELOCITY_X) == false ||
				        this->GetGeometry()[i].HasDofFor(VELOCITY_Y) == false ||
				        this->GetGeometry()[i].HasDofFor(VELOCITY_Z) == false)
					KRATOS_ERROR(std::invalid_argument, "missing VELOCITY component degree of freedom on node ", this->GetGeometry()[i].Id());
				if (this->GetGeometry()[i].HasDofFor(PRESSURE) == false)
					KRATOS_ERROR(std::invalid_argument, "missing PRESSURE component degree of freedom on node ", this->GetGeometry()[i].Id());
			}
			// Not checking OSS related variables NODAL_AREA, ADVPROJ, DIVPROJ, which are only required as SolutionStepData if OSS_SWITCH == 1
			// If this is a 2D problem, check that nodes are in XY plane
			if (this->GetGeometry().WorkingSpaceDimension() == 2)
			{
				for (unsigned int i = 0; i<this->GetGeometry().size(); ++i)
				{
					if (this->GetGeometry()[i].Z() != 0.0)
						KRATOS_ERROR(std::invalid_argument, "Node with non-zero Z coordinate found. Id: ", this->GetGeometry()[i].Id());
				}
			}
			return 0;
			KRATOS_CATCH("");
		}
		///@}
		///@name Access
		///@{
		///@}
		///@name Elemental Data
		///@{
		///@}
		///@name Inquiry
		///@{
		///@}
		///@name Input and output
		///@{
		/// Turn back information as a string.
		virtual std::string Info() const
		{
			std::stringstream buffer;
			buffer << "TwoFluidVMS #" << this->Id();
			return buffer.str();
		}
		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const
		{
			rOStream << "TwoFluidVMS" << TDim << "D";
		}
		//        /// Print object's data.
		//        virtual void PrintData(std::ostream& rOStream) const;
		///@}
		///@name Friends
		///@{
		///@}
	protected:
		///@name Protected static Member Variables
		///@{
		///@}
		///@name Protected member Variables
		///@{
		///@}
		///@name Protected Operators
		///@{
		///@}
		///@name Protected Operations
		///@{
		/// Add lumped mass matrix
		void LumpMassMatrix(MatrixType& rMassMatrix)
		{
			for (unsigned int i = 0; i < rMassMatrix.size1(); ++i)
			{
				double diag_factor = 0.0;
				for (unsigned int j = 0; j < rMassMatrix.size2(); ++j)
				{
					diag_factor += rMassMatrix(i, j);
					rMassMatrix(i, j) = 0.0;
				}
				rMassMatrix(i, i) = diag_factor;
			}
		}
		/// Add the weighted value of a variable at a point inside the element to a double
		/**
		 * Evaluate a scalar variable in the point where the form functions take the
		 * values given by rShapeFunc and add the result, weighted by Weight, to
		 * rResult. This is an auxiliary function used to compute values in integration
		 * points.
		 * @param rResult: The double where the value will be added to
		 * @param rVariable: The nodal variable to be read
		 * @param rShapeFunc: The values of the form functions in the point
		 * @param Step: The time Step (Defaults to 0 = Current)
		 * @param Weight: The variable will be weighted by this value before it is added to rResult
		 */
		virtual void AddPointContribution(double& rResult,
		                                  const Variable< double >& rVariable,
		                                  const array_1d< double, TNumNodes >& rShapeFunc,
		                                  const double Weight = 1.0)
		{
			double temp = 0.0;
			this->EvaluateInPoint(temp, rVariable, rShapeFunc);
			rResult += Weight*temp;
		}
		/// Write the value of a variable at a point inside the element to a double
		/**
		 * Evaluate a scalar variable in the point where the form functions take the
		 * values given by rShapeFunc and write the result to rResult.
		 * This is an auxiliary function used to compute values in integration points.
		 * @param rResult: The double where the value will be added to
		 * @param rVariable: The nodal variable to be read
		 * @param rShapeFunc: The values of the form functions in the point
		 * @param Step: The time Step (Defaults to 0 = Current)
		 */
		virtual void EvaluateInPoint(double& rResult,
		                             const Variable< double >& rVariable,
		                             const array_1d< double, TNumNodes >& rShapeFunc)
		{
			//compute sign of distance on gauss point
			double dist = 0.0;
			for (unsigned int i = 0; i < TNumNodes; i++)
				dist += rShapeFunc[i] * this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
//KRATOS_WATCH(dist)
//
//double test=0.0;
//for (unsigned int i = 0; i < TNumNodes; i++)
//    test+= rShapeFunc[i];
//if(test < 0.9999999)
//    KRATOS_ERROR(std::logic_error,"shape functions do not sum to 1","")
			double navg = 0.0;
			double value = 0.0;
			for (unsigned int i = 0; i < TNumNodes; i++)
			{
				if ( (dist * this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE)) > 0.0)
				{
					navg += 1.0;
					value += this->GetGeometry()[i].FastGetSolutionStepValue(rVariable);
//                    KRATOS_WATCH(this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE))
//                    KRATOS_WATCH(this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE))
//                    KRATOS_WATCH(value);
//                    KRATOS_WATCH(navg);
				}
			}
			if(navg != 0)
				value /= navg;
			else
				ElementBaseType::EvaluateInPoint(value,rVariable,rShapeFunc);
			rResult = value;
//            if(rVariable==DENSITY)
//                KRATOS_WATCH(rResult)
		}
		/// Add the weighted value of a variable at a point inside the element to a vector
		/**
		 * Evaluate a vector variable in the point where the form functions take the
		 * values given by rShapeFunc and add the result, weighted by Weight, to
		 * rResult. This is an auxiliary function used to compute values in integration
		 * points.
		 * @param rResult: The vector where the value will be added to
		 * @param rVariable: The nodal variable to be read
		 * @param rShapeFunc: The values of the form functions in the point
		 * @param Weight: The variable will be weighted by this value before it is added to rResult
		 */
		virtual void AddPointContribution(array_1d< double, 3 > & rResult,
		                                  const Variable< array_1d< double, 3 > >& rVariable,
		                                  const array_1d< double, TNumNodes>& rShapeFunc,
		                                  const double Weight = 1.0)
		{
			array_1d<double, 3 > temp(3, 0.0);
			this->EvaluateInPoint(temp, rVariable, rShapeFunc);
			rResult += Weight*temp;
		}
		/// Write the value of a variable at a point inside the element to a double
		/**
		 * Evaluate a scalar variable in the point where the form functions take the
		 * values given by rShapeFunc and write the result to rResult.
		 * This is an auxiliary function used to compute values in integration points.
		 * @param rResult: The double where the value will be added to
		 * @param rVariable: The nodal variable to be read
		 * @param rShapeFunc: The values of the form functions in the point
		 */
		virtual void EvaluateInPoint(array_1d< double, 3 > & rResult,
		                             const Variable< array_1d< double, 3 > >& rVariable,
		                             const array_1d< double, TNumNodes >& rShapeFunc)
		{
			//compute sign of distance on gauss point
			double dist = 0.0;
			for (unsigned int i = 0; i < TNumNodes; i++)
				dist += rShapeFunc[i] * this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
			double navg = 0.0;
			array_1d< double, 3 > value(3, 0.0);
			for (unsigned int i = 0; i < TNumNodes; i++)
			{
				if (dist * this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE) > 0.0)
				{
					navg += 1.0;
					value += this->GetGeometry()[i].FastGetSolutionStepValue(rVariable);
				}
			}
			if(navg != 0)
				value /= navg;
			else
				ElementBaseType::EvaluateInPoint(value,rVariable,rShapeFunc);
			rResult = value;
		}
		/// Return an estimate for the element size h, used to calculate the stabilization parameters
		/**
		 * Estimate the element size from its area or volume, required to calculate stabilization parameters.
		 * Note that its implementation is different for 2D or 3D elements.
		 * @see TwoFluidVMS2D, TwoFluidVMS3D for actual implementation
		 * @param Volume (in 3D) or Area (in 2D) of the element
		 * @return Element size h
		 */
		double ElementSize(const double);
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
		///@name Static Member Variables
		///@{
		///@}
		///@name Member Variables
		///@{
		///@}
		///@name Serialization
		///@{
		friend class Serializer;
		virtual void save(Serializer& rSerializer) const
		{
			KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ElementBaseType);
		}
		virtual void load(Serializer& rSerializer)
		{
			KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ElementBaseType);
		}
		///@}
		///@name Private Operators
		///@{
		///@}
		///@name Private Operations
		///@{
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
		TwoFluidVMS & operator=(TwoFluidVMS const& rOther);
		/// Copy constructor.
		TwoFluidVMS(TwoFluidVMS const& rOther);
		///@}
}; // Class TwoFluidVMS
///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
/// input stream function
template< unsigned int TDim,
          unsigned int TNumNodes >
inline std::istream & operator >>(std::istream& rIStream,
                                  TwoFluidVMS<TDim, TNumNodes>& rThis)
{
	return rIStream;
}
/// output stream function
template< unsigned int TDim,
          unsigned int TNumNodes >
inline std::ostream & operator <<(std::ostream& rOStream,
                                  const TwoFluidVMS<TDim, TNumNodes>& rThis)
{
	rThis.PrintInfo(rOStream);
	rOStream << std::endl;
	rThis.PrintData(rOStream);
	return rOStream;
}
///@}
///@} // Fluid Dynamics Application group
} // namespace Kratos.
#endif // KRATOS_TWO_FLUID_TwoFluidVMS_H_INCLUDED  defined
