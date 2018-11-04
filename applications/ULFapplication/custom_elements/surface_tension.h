//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main author:    Jordi Cotela
//  author     :    Alex Jarauta
//



#if !defined(KRATOS_SURFACE_TENSION_H_INCLUDED)
#define  KRATOS_SURFACE_TENSION_H_INCLUDED


// System includes
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>
#include <stdlib.h>
#include <iomanip>
#include <math.h>


// External includes
// #include "boost/smart_ptr.hpp"



// Project includes
#include <pybind11/pybind11.h>
#include "includes/define.h"
#include "includes/define_python.h"

#include "containers/array_1d.h"

#include "includes/checks.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/serializer.h"
#include "ULF_application_variables.h"
#include "includes/variables.h"
// #include "ULF_application_variables.h"
#include "includes/cfd_variables.h"
#include "utilities/geometry_utilities.h"
#include "includes/deprecated_variables.h"

namespace Kratos
{

///@addtogroup ULF
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

/// A stabilized element for the incompressible Navier-Stokes equations, utilizing lagrangian_Eulerian approach for droplet dynamics.
/**
 * This class is based on the VMS element in the fluid dynamics,
 * It is extended to solve for the droplet dynamics, and surface tensions.
 * for more information about the the governer equations and thier discretization, please find it in the below references:
 * Alex Jarauta, Pavel Ryzhakov, Marc Secanell, Prashant R Waghmare, and Jordi Pons-Prats. Numerical study of droplet dynamics in a polymer electrolyte fuel cell gas channel using an embedded eulerian-lagrangian approach. Journal of Power Sources, 323:201–212, 2016.
 */

template< unsigned int TDim,
          unsigned int TNumNodes = TDim + 1 >
class SurfaceTension : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SurfaceTension
    KRATOS_CLASS_POINTER_DEFINITION(SurfaceTension);

    ///base type: an IndexedObject that automatically has a unique number
    typedef IndexedObject BaseType;

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

    typedef Matrix MatrixType;

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    typedef std::vector<std::size_t> EquationIdVectorType;

    typedef std::vector< Dof<double>::Pointer > DofsVectorType;

    typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;

    typedef VectorMap<IndexType, DataValueContainer> SolutionStepsElementalDataContainerType;

     typedef array_1d<double, TNumNodes> ShapeFunctionsType;
     typedef BoundedMatrix<double, TNumNodes, TDim> ShapeFunctionDerivativesType;

    ///@}
    ///@name Life Cycle
    ///@{

    //Constructors.

    /// Default constuctor.
    /**
     * @param NewId Index number of the new element (optional)
     */
    SurfaceTension(IndexType NewId = 0) :
        Element(NewId)
    {}

    ///Constructor using an array of nodes.
    /**
     * @param NewId Index of the new element
     * @param ThisNodes An array containing the nodes of the new element
     */
    SurfaceTension(IndexType NewId, const NodesArrayType& ThisNodes) :
        Element(NewId, ThisNodes)
    {}

    /// Constructor using a geometry object.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     */
    SurfaceTension(IndexType NewId, GeometryType::Pointer pGeometry) :
        Element(NewId, pGeometry)
    {}

    /// Constuctor using geometry and properties.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     * @param pProperties Pointer to the element's properties
     */
    SurfaceTension(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) :
        Element(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    ~SurfaceTension() override
    {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /// Create a new element of this type
    /**
     * Returns a pointer to a new SurfaceTension element, created using given input
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
                            PropertiesType::Pointer pProperties) const override
    {

        return Kratos::make_shared< SurfaceTension<TDim, TNumNodes> >(NewId, GetGeometry().Create(ThisNodes), pProperties);

        //return Element::Pointer(new SurfaceTension(NewId, GetGeometry().Create(ThisNodes), pProperties));
    }

      Element::Pointer Create(IndexType NewId,
                             GeometryType::Pointer pGeom,
                             PropertiesType::Pointer pProperties) const override
      {
          return Kratos::make_shared< SurfaceTension<TDim, TNumNodes> >(NewId, pGeom, pProperties);
      }

    /// Provides local contributions from body forces and OSS projection terms
    /**
     * This is called during the assembly process and provides the terms of the
     * system that are either constant or computed explicitly (from the 'old'
     * iteration variables). In this case this means the body force terms and the
     * OSS projections, that are treated explicitly.
     * @param rLeftHandSideMatrix: the elemental left hand side matrix. Not used here, required for compatibility purposes only.
     * @param rRightHandSideVector: the elemental right hand side
     * @param rCurrentProcessInfo: the current process info
     */
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      ProcessInfo& rCurrentProcessInfo) override
    {
        const unsigned int LocalSize = (TDim + 1) * TNumNodes;

        // Check sizes and initialize matrix
        if (rLeftHandSideMatrix.size1() != LocalSize)
            rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);

        noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);

        // Calculate RHS
        this->CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
    }

    /// Returns a zero matrix of appropiate size (provided for compatibility with scheme)
    /**
     * @param rLeftHandSideMatrix Local matrix, will be filled with zeros
     * @param rCurrentProcessInfo Process info instance
     */
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                       ProcessInfo& rCurrentProcessInfo) override
    {
        const unsigned int LocalSize = (TDim + 1) * TNumNodes;

        if (rLeftHandSideMatrix.size1() != LocalSize)
            rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);

        noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);
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
    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                        ProcessInfo& rCurrentProcessInfo) override
    {
        const unsigned int LocalSize = (TDim + 1) * TNumNodes;

        // Check sizes and initialize
        if (rRightHandSideVector.size() != LocalSize)
            rRightHandSideVector.resize(LocalSize, false);

        noalias(rRightHandSideVector) = ZeroVector(LocalSize);

        // Calculate this element's geometric parameters
        double Area;
        array_1d<double, TNumNodes> N;
        BoundedMatrix<double, TNumNodes, TDim> DN_DX;
        GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

        // Calculate this element's fluid properties
        double Density;
        this->EvaluateInPoint(Density, DENSITY, N);

        // Calculate Momentum RHS contribution
        this->AddMomentumRHS(rRightHandSideVector, Density, N, Area);

        // For OSS: Add projection of residuals to RHS
        if (rCurrentProcessInfo[OSS_SWITCH] == 1)
        {
            array_1d<double, 3 > AdvVel;
            this->GetAdvectiveVel(AdvVel, N);

            double ElemSize = this->ElementSize(Area);
            double Viscosity = this->EffectiveViscosity(Density,N,DN_DX,ElemSize,rCurrentProcessInfo);

            // stabilization parameters
            double TauOne, TauTwo;
            this->CalculateTau(TauOne,TauTwo,AdvVel,ElemSize,Density,Viscosity,rCurrentProcessInfo);

            this->AddProjectionToRHS(rRightHandSideVector, AdvVel, Density, TauOne, TauTwo, N, DN_DX, Area,rCurrentProcessInfo[DELTA_TIME]);
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
    void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) override
    {
        const unsigned int LocalSize = (TDim + 1) * TNumNodes;

        // Resize and set to zero
        if (rMassMatrix.size1() != LocalSize)
            rMassMatrix.resize(LocalSize, LocalSize, false);

        rMassMatrix = ZeroMatrix(LocalSize, LocalSize);

        // Get the element's geometric parameters
        double Area;
        array_1d<double, TNumNodes> N;
        BoundedMatrix<double, TNumNodes, TDim> DN_DX;
        GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

        // Calculate this element's fluid properties
        double Density;
        this->EvaluateInPoint(Density, DENSITY, N);

        // Add 'classical' mass matrix (lumped)
        double Coeff = Density * Area / TNumNodes; //Optimize!
        this->CalculateLumpedMassMatrix(rMassMatrix, Coeff);

        /* For ASGS: add dynamic stabilization terms.
         These terms are not used in OSS, as they belong to the finite element
         space and cancel out with their projections.
         */
        if (rCurrentProcessInfo[OSS_SWITCH] != 1)
        {
            double ElemSize = this->ElementSize(Area);
            double Viscosity = this->EffectiveViscosity(Density,N,DN_DX,ElemSize,rCurrentProcessInfo);

            // Get Advective velocity
            array_1d<double, 3 > AdvVel;
            this->GetAdvectiveVel(AdvVel, N);

            // stabilization parameters
            double TauOne, TauTwo;
            this->CalculateTau(TauOne,TauTwo,AdvVel,ElemSize,Density,Viscosity,rCurrentProcessInfo);

            // Add dynamic stabilization terms ( all terms involving a delta(u) )
            this->AddMassStabTerms(rMassMatrix, Density, AdvVel, TauOne, N, DN_DX, Area);
        }
    }

    /// Computes the local contribution associated to 'new' velocity and pressure values
    /**
     * Provides local contributions to the system associated to the velocity and
     * pressure terms (convection, diffusion, pressure gradient/velocity divergence
     * and stabilization).
     * @param rDampingMatrix Will be filled with the velocity-proportional "damping" matrix
     * @param rRightHandSideVector the elemental right hand side vector
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateLocalVelocityContribution(MatrixType& rDampingMatrix,
            VectorType& rRightHandSideVector,
            ProcessInfo& rCurrentProcessInfo) override
    {
        const unsigned int LocalSize = (TDim + 1) * TNumNodes;

        // Resize and set to zero the matrix
        // Note that we don't clean the RHS because it will already contain body force (and stabilization) contributions
        if (rDampingMatrix.size1() != LocalSize)
            rDampingMatrix.resize(LocalSize, LocalSize, false);

        noalias(rDampingMatrix) = ZeroMatrix(LocalSize, LocalSize);

        // Get this element's geometric properties
        double Area;
        array_1d<double, TNumNodes> N;
        BoundedMatrix<double, TNumNodes, TDim> DN_DX;
        GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

        // Calculate this element's fluid properties
        double Density;
        this->EvaluateInPoint(Density, DENSITY, N);

        double ElemSize = this->ElementSize(Area);
        double Viscosity = this->EffectiveViscosity(Density,N,DN_DX,ElemSize,rCurrentProcessInfo);

        // Get Advective velocity
        array_1d<double, 3 > AdvVel;
        this->GetAdvectiveVel(AdvVel, N);

        // stabilization parameters
        double TauOne, TauTwo;
        this->CalculateTau(TauOne,TauTwo,AdvVel,ElemSize,Density,Viscosity,rCurrentProcessInfo);

        this->AddIntegrationPointVelocityContribution(rDampingMatrix, rRightHandSideVector, Density, Viscosity, AdvVel, TauOne, TauTwo, N, DN_DX, Area);



	// Surface tension contribution
	int k = 0;
	if(TDim < 3)
	{
	    array_1d<double,3> node_indx;
	    node_indx[0] = 0.0;
	    node_indx[1] = 0.0;
	    node_indx[2] = 0.0;
	    int node_indx_wrong = 5;
	    for(unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
	    {
	      if(this->GetGeometry()[iNode].FastGetSolutionStepValue(MEAN_CURVATURE_2D) != 0.0)
		k++;
	      if(this->GetGeometry()[iNode].FastGetSolutionStepValue(IS_BOUNDARY) < 0.1)
		node_indx_wrong = iNode;
	    }
	    if(k > 2 && node_indx_wrong != 5)
	    {
	      this->GetGeometry()[node_indx_wrong].FastGetSolutionStepValue(MEAN_CURVATURE_2D) = 0.0;
	    }
	    k = 0;
	    for(unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
	    {
	      if(this->GetGeometry()[iNode].FastGetSolutionStepValue(IS_FREE_SURFACE) != 0.0 || this->GetGeometry()[iNode].FastGetSolutionStepValue(TRIPLE_POINT) != 0.0)
	      {
		  node_indx[k] = iNode;
		  k++;
	      }
	    }
	    if(k > 1)
	      this->ApplySurfaceTensionContribution(rDampingMatrix, rRightHandSideVector, node_indx, k, rCurrentProcessInfo);
	}
	else
	{
	    array_1d<double,4> node_indx;
	    node_indx[0] = 0.0;
	    node_indx[1] = 0.0;
	    node_indx[2] = 0.0;
	    node_indx[3] = 0.0;
	    int node_indx_wrong = 5;
	    for(unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
	    {
	      if(this->GetGeometry()[iNode].FastGetSolutionStepValue(MEAN_CURVATURE_3D) != 0.0)
		k++;
	      if(this->GetGeometry()[iNode].FastGetSolutionStepValue(IS_BOUNDARY) < 0.1)
		node_indx_wrong = iNode;
	    }
	    if(k > 2 && node_indx_wrong != 5)
	    {
	      this->GetGeometry()[node_indx_wrong].FastGetSolutionStepValue(MEAN_CURVATURE_3D) = 0.0;
	    }
	    k = 0;
	    for(unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
	    {
	      if(this->GetGeometry()[iNode].FastGetSolutionStepValue(IS_FREE_SURFACE) > 0.999 || this->GetGeometry()[iNode].FastGetSolutionStepValue(TRIPLE_POINT) != 0.0)
	      {
		  node_indx[k] = iNode;
		  k++;
	      }
	    }
	    if(k >= 3)
	      this->ApplySurfaceTensionContribution3D(rDampingMatrix, rRightHandSideVector, node_indx, k, rCurrentProcessInfo);
	}


	// Viscous stress
	k = 0;
	for(unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
	{
	  if(this->GetGeometry()[iNode].FastGetSolutionStepValue(IS_WATER) == 0.0)
	  {
	    k++;
	  }
	}
	if(TDim < 3 && k > 2)
	    this->AddViscousStress2D();

        // Now calculate an additional contribution to the residual: r -= rDampingMatrix * (u,p)
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

        noalias(rRightHandSideVector) -= prod(rDampingMatrix, U);
    }

    void FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override
    {
    }

    /// Implementation of Calculate to compute an error estimate.
    /**
     * If rVariable == ERROR_RATIO, this function will provide an a posteriori
     * estimate of the norm of the subscale velocity, calculated as TauOne*||MomentumResidual||.
     * Note that the residual of the momentum equation is evaluated at the element center
     * and that the result has units of velocity (L/T).
     * The error estimate both saved as the elemental ERROR_RATIO variable and returned as rOutput.
     * If rVARIABLE == NODAL_AREA, the element's contribution to nodal area is added to its nodes.
     * @param rVariable Use ERROR_RATIO or NODAL_AREA
     * @param rOutput Returns the error estimate for ERROR_RATIO, unused for NODAL_AREA
     * @param rCurrentProcessInfo Process info instance (will be checked for OSS_SWITCH)
     * @see MarkForRefinement for a use of the error ratio
     */
   void Calculate(const Variable<double>& rVariable,
                           double& rOutput,
                           const ProcessInfo& rCurrentProcessInfo) override
    {
        if (rVariable == ERROR_RATIO)
        {
            rOutput = this->SubscaleErrorEstimate(rCurrentProcessInfo);
            this->SetValue(ERROR_RATIO,rOutput);
        }
        else if (rVariable == NODAL_AREA)
        {
            // Get the element's geometric parameters
            double Area;
            array_1d<double, TNumNodes> N;
            BoundedMatrix<double, TNumNodes, TDim> DN_DX;
            GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

            // Carefully write results to nodal variables, to avoid parallelism problems
            for (unsigned int i = 0; i < TNumNodes; ++i)
            {
                this->GetGeometry()[i].SetLock(); // So it is safe to write in the node in OpenMP
                this->GetGeometry()[i].FastGetSolutionStepValue(NODAL_AREA) += Area * N[i];
                this->GetGeometry()[i].UnSetLock(); // Free the node for other threads
            }
        }
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
   void Calculate(const Variable<array_1d<double, 3 > >& rVariable,
                           array_1d<double, 3 > & rOutput,
                           const ProcessInfo& rCurrentProcessInfo) override
    {
        if (rVariable == ADVPROJ) // Compute residual projections for OSS
        {
            // Get the element's geometric parameters
            double Area;
            array_1d<double, TNumNodes> N;
            BoundedMatrix<double, TNumNodes, TDim> DN_DX;
            GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

            // Calculate this element's fluid properties
            double Density;
            this->EvaluateInPoint(Density, DENSITY, N);

            // Get Advective velocity
            array_1d<double, 3 > AdvVel;
            this->GetAdvectiveVel(AdvVel, N);

            // Output containers
            array_1d< double, 3 > ElementalMomRes(3, 0.0);
            double ElementalMassRes(0);

            this->AddProjectionResidualContribution(AdvVel, Density, ElementalMomRes, ElementalMassRes, N, DN_DX, Area);

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
                    this->GetGeometry()[i].FastGetSolutionStepValue(NODAL_AREA) += Area * N[i];
                    this->GetGeometry()[i].UnSetLock(); // Free the node for other threads
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
            BoundedMatrix<double, TNumNodes, TDim> DN_DX;
            GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

            // Calculate this element's fluid properties
            double Density;
            this->EvaluateInPoint(Density, DENSITY, N);

            // Get Advective velocity
            array_1d<double, 3 > AdvVel;
            this->GetAdvectiveVel(AdvVel, N);

            // Output containers
            array_1d< double, 3 > ElementalMomRes(3,0.0);
            double ElementalMassRes(0.0);

            this->AddProjectionResidualContribution(AdvVel, Density, ElementalMomRes, ElementalMassRes, N, DN_DX, Area);

            if (rCurrentProcessInfo[OSS_SWITCH] == 1)
            {
                /* Projections of the elemental residual are computed with
                 * Newton-Raphson iterations of type M(lumped) dx = ElemRes - M(consistent) * x
                 */
                const double Weight = ConsistentMassCoef(Area); // Consistent mass matrix is Weigth * ( Ones(TNumNodes,TNumNodes) + Identity(TNumNodes,TNumNodes) )
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
                    this->GetGeometry()[i].FastGetSolutionStepValue(NODAL_AREA) += Area * N[i];

                    // Substract M(consistent)*x(i-1) from RHS
                    for(unsigned int j = 0; j < TNumNodes; ++j) // RHS -= Weigth * Ones(TNumNodes,TNumNodes) * x(i-1)
                    {
                        for(unsigned int d = 0; d < TDim; ++d)
                            rMomRHS[d] -= Weight * this->GetGeometry()[j].FastGetSolutionStepValue(ADVPROJ)[d];
                        rMassRHS -= Weight * this->GetGeometry()[j].FastGetSolutionStepValue(DIVPROJ);
                    }
                    for(unsigned int d = 0; d < TDim; ++d) // RHS -= Weigth * Identity(TNumNodes,TNumNodes) * x(i-1)
                        rMomRHS[d] -= Weight * this->GetGeometry()[i].FastGetSolutionStepValue(ADVPROJ)[d];
                    rMassRHS -= Weight * this->GetGeometry()[i].FastGetSolutionStepValue(DIVPROJ);

                    this->GetGeometry()[i].UnSetLock(); // Free the node for other threads
                }
            }

            /// Return output
            rOutput = ElementalMomRes;
        }
    }

    // The following methods have different implementations depending on TDim
    /// Provides the global indices for each one of this element's local rows
    /**
     * this determines the elemental equation ID vector for all elemental
     * DOFs
     * @param rResult A vector containing the global Id of each row
     * @param rCurrentProcessInfo the current process info object (unused)
     */
    void EquationIdVector(EquationIdVectorType& rResult,
                                  ProcessInfo& rCurrentProcessInfo) override;

    /// Returns a list of the element's Dofs
    /**
     * @param ElementalDofList the list of DOFs
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetDofList(DofsVectorType& rElementalDofList,
                            ProcessInfo& rCurrentProcessInfo) override;

    /// Returns VELOCITY_X, VELOCITY_Y, (VELOCITY_Z,) PRESSURE for each node
    /**
     * @param Values Vector of nodal unknowns
     * @param Step Get result from 'Step' steps back, 0 is current step. (Must be smaller than buffer size)
     */
    void GetFirstDerivativesVector(Vector& Values, int Step = 0) override;

    /// Returns ACCELERATION_X, ACCELERATION_Y, (ACCELERATION_Z,) 0 for each node
    /**
     * @param Values Vector of nodal second derivatives
     * @param Step Get result from 'Step' steps back, 0 is current step. (Must be smaller than buffer size)
     */
    void GetSecondDerivativesVector(Vector& Values, int Step = 0) override;

    /// Obtain an array_1d<double,3> elemental variable, evaluated on gauss points.
    /**
     * If the variable is VORTICITY, computes the vorticity (rotational of the velocity)
     * based on the current velocity values. Otherwise, it assumes that the input
     * variable is an elemental value and retrieves it. Implemented for a
     * single gauss point only.
     * @param rVariable Kratos vector variable to get
     * @param Output Will be filled with the values of the variable on integrartion points
     * @param rCurrentProcessInfo Process info instance
     */
    void GetValueOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
            std::vector<array_1d<double, 3 > >& rOutput,
            const ProcessInfo& rCurrentProcessInfo) override;

    /// Obtain a double elemental variable, evaluated on gauss points.
    /**
     * If the variable is TAUONE or TAUTWO, calculates the corresponding stabilization
     * parameter for the element, based on rCurrentProcessInfo's DELTA_TIME and
     * DYNAMIC_TAU. If the variable is MU, calculates the effective viscosity at the
     * element center due to Smagorinsky (in 'dynamic' units). Otherwise, it assumes
     * that the input variable is an elemental value and retrieves it.
     * Implemented for a single gauss point only.
     * @param rVariable Kratos vector variable to compute
     * @param Output Will be filled with the values of the variable on integrartion points
     * @param rCurrentProcessInfo Process info instance
     */
    void GetValueOnIntegrationPoints(const Variable<double>& rVariable,
            std::vector<double>& rValues,
            const ProcessInfo& rCurrentProcessInfo) override
    {
        if (rVariable == TAUONE || rVariable == TAUTWO || rVariable == MU || rVariable == TAU)
        {
            double Area;
            array_1d<double, TNumNodes> N;
            BoundedMatrix<double, TNumNodes, TDim> DN_DX;
            GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

            array_1d<double, 3 > AdvVel;
            this->GetAdvectiveVel(AdvVel, N);

            double Density;
            this->EvaluateInPoint(Density,DENSITY,N);

            double ElemSize = this->ElementSize(Area);
            double Viscosity = this->EffectiveViscosity(Density,N,DN_DX,ElemSize,rCurrentProcessInfo);

            // stabilization parameters
            double TauOne, TauTwo;
            this->CalculateTau(TauOne,TauTwo,AdvVel,ElemSize,Density,Viscosity,rCurrentProcessInfo);


            rValues.resize(1, false);
            if (rVariable == TAUONE)
            {
                rValues[0] = TauOne;
            }
            else if (rVariable == TAUTWO)
            {
                rValues[0] = TauTwo;
            }
            else if (rVariable == MU)
            {
                rValues[0] = Viscosity;
            }
            else if (rVariable == TAU)
            {
                double NormS = this->EquivalentStrainRate(DN_DX);
                rValues[0] = Viscosity*NormS;
            }
        }
        else if (rVariable == EQ_STRAIN_RATE)
        {
            double Area;
            array_1d<double, TNumNodes> N;
            BoundedMatrix<double, TNumNodes, TDim> DN_DX;
            GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

            rValues.resize(1, false);
            rValues[0] = this->EquivalentStrainRate(DN_DX);
        }
        else if(rVariable == SUBSCALE_PRESSURE)
        {
            double Area;
            array_1d<double, TNumNodes> N;
            BoundedMatrix<double, TNumNodes, TDim> DN_DX;
            GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

            array_1d<double, 3 > AdvVel;
            this->GetAdvectiveVel(AdvVel, N);

            double Density;
            this->EvaluateInPoint(Density,DENSITY,N);

            double ElemSize = this->ElementSize(Area);
            double Viscosity = this->EffectiveViscosity(Density,N,DN_DX,ElemSize,rCurrentProcessInfo);

            // stabilization parameters
            double TauOne, TauTwo;
            this->CalculateTau(TauOne,TauTwo,AdvVel,ElemSize,Density,Viscosity,rCurrentProcessInfo);

            double DivU = 0.0;
            for(unsigned int i=0; i < TNumNodes; i++)
            {
                for(unsigned int d = 0; d < TDim; d++)
                    DivU -= DN_DX(i,d) * this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[d];
            }

            rValues.resize(1, false);
            rValues[0] = TauTwo * DivU;
            if(rCurrentProcessInfo[OSS_SWITCH]==1)
            {
                double Proj = 0.0;
                for(unsigned int i=0; i < TNumNodes; i++)
                {
                    Proj += N[i]*this->GetGeometry()[i].FastGetSolutionStepValue(DIVPROJ);
                }
                rValues[0] -= TauTwo*Proj;
            }
        }
        else if (rVariable == NODAL_AREA && TDim == 3)
        {
            MatrixType J = ZeroMatrix(3,3);
            const array_1d<double,3>& X0 = this->GetGeometry()[0].Coordinates();
            const array_1d<double,3>& X1 = this->GetGeometry()[1].Coordinates();
            const array_1d<double,3>& X2 = this->GetGeometry()[2].Coordinates();
            const array_1d<double,3>& X3 = this->GetGeometry()[3].Coordinates();

            J(0,0) = X1[0]-X0[0];
            J(0,1) = X2[0]-X0[0];
            J(0,2) = X3[0]-X0[0];
            J(1,0) = X1[1]-X0[1];
            J(1,1) = X2[1]-X0[1];
            J(1,2) = X3[1]-X0[1];
            J(2,0) = X1[2]-X0[2];
            J(2,1) = X2[2]-X0[2];
            J(2,2) = X3[2]-X0[2];

            double DetJ = J(0,0)*( J(1,1)*J(2,2) - J(1,2)*J(2,1) ) + J(0,1)*( J(1,2)*J(2,0) - J(1,0)*J(2,2) ) + J(0,2)*( J(1,0)*J(2,1) - J(1,1)*J(2,0) );
            rValues.resize(1, false);
            rValues[0] = DetJ;
        }
        else if (rVariable == ERROR_RATIO)
        {
            rValues.resize(1,false);
            rValues[0] = this->SubscaleErrorEstimate(rCurrentProcessInfo);
        }
        else // Default behaviour (returns elemental data)
        {
            rValues.resize(1, false);
            /*
             The cast is done to avoid modification of the element's data. Data modification
             would happen if rVariable is not stored now (would initialize a pointer to &rVariable
             with associated value of 0.0). This is catastrophic if the variable referenced
             goes out of scope.
             */
            const SurfaceTension<TDim, TNumNodes>* const_this = static_cast<const SurfaceTension<TDim, TNumNodes>*> (this);
            rValues[0] = const_this->GetValue(rVariable);
        }
    }

    /// Empty implementation of unused CalculateOnIntegrationPoints overloads to avoid compilation warning
   void GetValueOnIntegrationPoints(const Variable<array_1d<double, 6 > >& rVariable,
            std::vector<array_1d<double, 6 > >& rValues,
            const ProcessInfo& rCurrentProcessInfo) override
    {}

    /// Empty implementation of unused CalculateOnIntegrationPoints overloads to avoid compilation warning
    void GetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
            std::vector<Vector>& rValues,
            const ProcessInfo& rCurrentProcessInfo) override
    {}

    /// Empty implementation of unused CalculateOnIntegrationPoints overloads to avoid compilation warning
    void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
            std::vector<Matrix>& rValues,
            const ProcessInfo& rCurrentProcessInfo) override
    {}

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Elemental Data
    ///@{

    /// Checks the input and that all required Kratos variables have been registered.
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo The ProcessInfo of the ModelPart that contains this element.
     * @return 0 if no errors were found.
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        // Perform basic element checks
        int ErrorCode = Kratos::Element::Check(rCurrentProcessInfo);
        if(ErrorCode != 0) return ErrorCode;

        // Check that all required variables have been registered
        KRATOS_CHECK_VARIABLE_KEY(VELOCITY);
        KRATOS_CHECK_VARIABLE_KEY(MESH_VELOCITY);
        KRATOS_CHECK_VARIABLE_KEY(ACCELERATION);
        KRATOS_CHECK_VARIABLE_KEY(PRESSURE);
        KRATOS_CHECK_VARIABLE_KEY(DENSITY);
        KRATOS_CHECK_VARIABLE_KEY(VISCOSITY);
        KRATOS_CHECK_VARIABLE_KEY(BODY_FORCE);
        KRATOS_CHECK_VARIABLE_KEY(OSS_SWITCH);
        KRATOS_CHECK_VARIABLE_KEY(DYNAMIC_TAU);
        KRATOS_CHECK_VARIABLE_KEY(DELTA_TIME);
        KRATOS_CHECK_VARIABLE_KEY(ADVPROJ);
        KRATOS_CHECK_VARIABLE_KEY(DIVPROJ);
        KRATOS_CHECK_VARIABLE_KEY(NODAL_AREA);
        KRATOS_CHECK_VARIABLE_KEY(C_SMAGORINSKY);
        KRATOS_CHECK_VARIABLE_KEY(ERROR_RATIO);
        // Additional variables, only required to print results:
        // SUBSCALE_VELOCITY, SUBSCALE_PRESSURE, TAUONE, TAUTWO, MU, VORTICITY.

        // Checks on nodes

        // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
        for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
        {
            Node<3> &rNode = this->GetGeometry()[i];
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY,rNode);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PRESSURE,rNode);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(MESH_VELOCITY,rNode);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ACCELERATION,rNode);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DENSITY,rNode);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VISCOSITY,rNode);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(BODY_FORCE,rNode);
            // Not checking OSS related variables NODAL_AREA, ADVPROJ, DIVPROJ, which are only required as SolutionStepData if OSS_SWITCH == 1

            KRATOS_CHECK_DOF_IN_NODE(VELOCITY_X,rNode);
            KRATOS_CHECK_DOF_IN_NODE(VELOCITY_Y,rNode);
            if (TDim == 3) KRATOS_CHECK_DOF_IN_NODE(VELOCITY_Z,rNode);
            KRATOS_CHECK_DOF_IN_NODE(PRESSURE,rNode);
        }
        // Not checking OSS related variables NODAL_AREA, ADVPROJ, DIVPROJ, which are only required as SolutionStepData if OSS_SWITCH == 1

        // If this is a 2D problem, check that nodes are in XY plane
        if (this->GetGeometry().WorkingSpaceDimension() == 2)
        {
            for (unsigned int i=0; i<this->GetGeometry().size(); ++i)
            {
                if (this->GetGeometry()[i].Z() != 0.0)
                    KRATOS_ERROR << "Node " << this->GetGeometry()[i].Id() << "has non-zero Z coordinate." << std::endl;            }
        }

        return 0;

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "SurfaceTension #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SurfaceTension" << TDim << "D";
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

    /// Calculate Stabilization parameters.
    /**
     * Calculates both tau parameters based on a given advective velocity.
     * Takes time step and dynamic coefficient from given ProcessInfo instance.
     * ProcessInfo variables DELTA_TIME and DYNAMIC_TAU will be used.
     * @param TauOne First stabilization parameter (momentum equation)
     * @param TauTwo Second stabilization parameter (mass equation)
     * @param rAdvVel advection velocity
     * @param ElemSize Characteristic element length
     * @param Density Density on integrartion point
     * @param Viscosity Dynamic viscosity (mu) on integrartion point
     * @param rCurrentProcessInfo Process info instance
     */
    virtual void CalculateTau(double& TauOne,
                              double& TauTwo,
                              const array_1d< double, 3 > & rAdvVel,
                              const double ElemSize,
                              const double Density,
                              const double Viscosity,
                              const ProcessInfo& rCurrentProcessInfo)
    {
        // Compute mean advective velocity norm
        double AdvVelNorm = 0.0;
        for (unsigned int d = 0; d < TDim; ++d)
            AdvVelNorm += rAdvVel[d] * rAdvVel[d];

        AdvVelNorm = sqrt(AdvVelNorm);

        double InvTau = Density * ( rCurrentProcessInfo[DYNAMIC_TAU] / rCurrentProcessInfo[DELTA_TIME] + 2.0*AdvVelNorm / ElemSize ) + 4.0*Viscosity/ (ElemSize * ElemSize);
        TauOne = 1.0 / InvTau;
        TauTwo = Viscosity + 0.5 * Density * ElemSize * AdvVelNorm;
    }

    /// Calculate momentum stabilization parameter (without time term).
    /**
     * Calculates the momentum tau parameter based on a given advective velocity.
     * The dynamic term is not taken into account. This function
     * is intended for error estimation only. In other cases use CalculateTau
     * @param TauOne First stabilization parameter (momentum equation)
     * @param rAdvVel advection velocity
     * @param ElemSize Characteristic element length
     * @param Density Density on integrartion point
     * @param Viscosity Dynamic viscosity (mu) on integrartion point
     */
    virtual void CalculateStaticTau(double& TauOne,
                                    const array_1d< double, 3 > & rAdvVel,
                                    const double ElemSize,
                                    const double Density,
                                    const double Viscosity)
    {
        // Compute mean advective velocity norm
        double AdvVelNorm = 0.0;
        for (unsigned int d = 0; d < TDim; ++d)
            AdvVelNorm += rAdvVel[d] * rAdvVel[d];

        AdvVelNorm = sqrt(AdvVelNorm);

        double InvTau = 2.0*Density*AdvVelNorm / ElemSize + 4.0*Viscosity/ (ElemSize * ElemSize);
        TauOne = 1.0 / InvTau;
    }

    /// Add the momentum equation contribution to the RHS (body forces)
    virtual void AddMomentumRHS(VectorType& F,
                                const double Density,
                                const array_1d<double, TNumNodes>& rShapeFunc,
                                const double Weight)
    {
        double Coef = Density * Weight;

        array_1d<double, 3 > BodyForce(3, 0.0);
        this->EvaluateInPoint(BodyForce, BODY_FORCE, rShapeFunc);

        // Add the results to the velocity components (Local Dofs are vx, vy, [vz,] p for each node)
        int LocalIndex = 0;

        for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
        {
            for (unsigned int d = 0; d < TDim; ++d)
            {
                F[LocalIndex++] += Coef * rShapeFunc[iNode] * BodyForce[d];
            }
            ++LocalIndex; // Skip pressure Dof
        }
    }

    /// Add OSS projection terms to the RHS
    virtual void AddProjectionToRHS(VectorType& RHS,
                                    const array_1d<double, 3 > & rAdvVel,
                                    const double Density,
                                    const double TauOne,
                                    const double TauTwo,
                                    const array_1d<double, TNumNodes>& rShapeFunc,
                                    const BoundedMatrix<double, TNumNodes, TDim>& rShapeDeriv,
                                    const double Weight,
                                    const double DeltaTime = 1.0)
    {
        const unsigned int BlockSize = TDim + 1;

        array_1d<double, TNumNodes> AGradN;
        this->GetConvectionOperator(AGradN, rAdvVel, rShapeDeriv); // Get a * grad(Ni)

        array_1d<double,3> MomProj(3,0.0);
        double DivProj = 0.0;
        this->EvaluateInPoint(MomProj,ADVPROJ,rShapeFunc);
        this->EvaluateInPoint(DivProj,DIVPROJ,rShapeFunc);

        MomProj *= TauOne;
        DivProj *= TauTwo;

        unsigned int FirstRow = 0;

        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            for (unsigned int d = 0; d < TDim; d++)
            {
                RHS[FirstRow+d] -= Weight * (Density * AGradN[i] * MomProj[d] + rShapeDeriv(i,d) * DivProj); // TauOne * ( a * Grad(v) ) * MomProjection + TauTwo * Div(v) * MassProjection
                RHS[FirstRow+TDim] -= Weight * rShapeDeriv(i,d) * MomProj[d]; // TauOne * Grad(q) * MomProjection
            }
            FirstRow += BlockSize;
        }
    }

    /// Add lumped mass matrix
    /**
     * Adds the lumped mass matrix to an elemental LHS matrix. Note that the time factor
     * (typically 1/(k*Dt) ) is added by the scheme outside the element.
     * @param rLHSMatrix: The local matrix where the result will be added
     * @param Mass: The weight assigned to each node (typically Density * Area / NumNodes or Density*Volume / NumNodes)
     */
    void CalculateLumpedMassMatrix(MatrixType& rLHSMatrix,
                                   const double Mass)
    {
        unsigned int DofIndex = 0;
        for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
        {
            for (unsigned int d = 0; d < TDim; ++d)
            {
                rLHSMatrix(DofIndex, DofIndex) += Mass;
                ++DofIndex;
            }
            ++DofIndex; // Skip pressure Dof
        }
    }

    void AddConsistentMassMatrixContribution(MatrixType& rLHSMatrix,
            const array_1d<double,TNumNodes>& rShapeFunc,
            const double Density,
            const double Weight)
    {
        const unsigned int BlockSize = TDim + 1;

        double Coef = Density * Weight;
        unsigned int FirstRow(0), FirstCol(0);
        double K; // Temporary results

        // Note: Dof order is (vx,vy,[vz,]p) for each node
        for (unsigned int i = 0; i < TNumNodes; ++i)
        {
            // Loop over columns
            for (unsigned int j = 0; j < TNumNodes; ++j)
            {
                // Delta(u) * TauOne * [ AdvVel * Grad(v) ] in velocity block
                K = Coef * rShapeFunc[i] * rShapeFunc[j];

                for (unsigned int d = 0; d < TDim; ++d) // iterate over dimensions for velocity Dofs in this node combination
                {
                    rLHSMatrix(FirstRow + d, FirstCol + d) += K;
                }
                // Update column index
                FirstCol += BlockSize;
            }
            // Update matrix indices
            FirstRow += BlockSize;
            FirstCol = 0;
        }
    }


    /// Add mass-like stabilization terms to LHS.
    /**
     * This function is only used in ASGS. For OSS, we avoid computing these
     * terms, as they shoud cancel out with the dynamic part of the projection
     * (which is not computed either)
     * @param rLHSMatrix Left hand side of the velocity-pressure system
     * @param Density Density on integration point
     * @param rAdvVel Advective velocity on integration point
     * @param TauOne Stabilization parameter for momentum equation
     * @param rShapeFunc Shape funcitions evaluated on integration point
     * @param rShapeDeriv Shape function derivatives evaluated on integration point
     * @param Weight Area (or volume) times integration point weight
     */
    void AddMassStabTerms(MatrixType& rLHSMatrix,
                          const double Density,
                          const array_1d<double, 3 > & rAdvVel,
                          const double TauOne,
                          const array_1d<double, TNumNodes>& rShapeFunc,
                          const BoundedMatrix<double, TNumNodes, TDim>& rShapeDeriv,
                          const double Weight)
    {
        const unsigned int BlockSize = TDim + 1;

        double Coef = Weight * TauOne;
        unsigned int FirstRow(0), FirstCol(0);
        double K; // Temporary results

        // If we want to use more than one Gauss point to integrate the convective term, this has to be evaluated once per integration point
        array_1d<double, TNumNodes> AGradN;
        this->GetConvectionOperator(AGradN, rAdvVel, rShapeDeriv); // Get a * grad(Ni)

        // Note: Dof order is (vx,vy,[vz,]p) for each node
        for (unsigned int i = 0; i < TNumNodes; ++i)
        {
            // Loop over columns
            for (unsigned int j = 0; j < TNumNodes; ++j)
            {
                // Delta(u) * TauOne * [ AdvVel * Grad(v) ] in velocity block
                K = Coef * Density * AGradN[i] * Density * rShapeFunc[j];

                for (unsigned int d = 0; d < TDim; ++d) // iterate over dimensions for velocity Dofs in this node combination
                {
                    rLHSMatrix(FirstRow + d, FirstCol + d) += K;
                    // Delta(u) * TauOne * Grad(q) in q * Div(u) block
                    rLHSMatrix(FirstRow + TDim, FirstCol + d) += Coef * Density * rShapeDeriv(i, d) * rShapeFunc[j];
                }
                // Update column index
                FirstCol += BlockSize;
            }
            // Update matrix indices
            FirstRow += BlockSize;
            FirstCol = 0;
        }
    }

    /// Add a the contribution from a single integration point to the velocity contribution
    void AddIntegrationPointVelocityContribution(MatrixType& rDampingMatrix,
            VectorType& rDampRHS,
            const double Density,
            const double Viscosity,
            const array_1d< double, 3 > & rAdvVel,
            const double TauOne,
            const double TauTwo,
            const array_1d< double, TNumNodes >& rShapeFunc,
            const BoundedMatrix<double, TNumNodes, TDim >& rShapeDeriv,
            const double Weight)
    {
        const unsigned int BlockSize = TDim + 1;

        // If we want to use more than one Gauss point to integrate the convective term, this has to be evaluated once per integration point
        array_1d<double, TNumNodes> AGradN;
        this->GetConvectionOperator(AGradN, rAdvVel, rShapeDeriv); // Get a * grad(Ni)

        // Build the local matrix and RHS
        unsigned int FirstRow(0), FirstCol(0); // position of the first term of the local matrix that corresponds to each node combination
        double K, G, PDivV, L, qF; // Temporary results

        array_1d<double,3> BodyForce(3,0.0);
        this->EvaluateInPoint(BodyForce,BODY_FORCE,rShapeFunc);
        BodyForce *= Density;

        for (unsigned int i = 0; i < TNumNodes; ++i) // iterate over rows
        {
            for (unsigned int j = 0; j < TNumNodes; ++j) // iterate over columns
            {
                // Calculate the part of the contributions that is constant for each node combination

                // Velocity block
                K = Density * rShapeFunc[i] * AGradN[j]; // Convective term: v * ( a * Grad(u) )
                //K = 0.5 * Density * (rShapeFunc[i] * AGradN[j] - AGradN[i] * rShapeFunc[j]); // Skew-symmetric convective term 1/2( v*grad(u)*u - grad(v) uu )
                K += TauOne * Density * AGradN[i] * Density * AGradN[j]; // Stabilization: (a * Grad(v)) * TauOne * (a * Grad(u))
                K *= Weight;

                // q-p stabilization block (reset result)
                L = 0;

                for (unsigned int m = 0; m < TDim; ++m) // iterate over v components (vx,vy[,vz])
                {
                    // Velocity block
                    //K += Weight * Viscosity * rShapeDeriv(i, m) * rShapeDeriv(j, m); // Diffusive term: Viscosity * Grad(v) * Grad(u)

                    // v * Grad(p) block
                    G = TauOne * Density * AGradN[i] * rShapeDeriv(j, m); // Stabilization: (a * Grad(v)) * TauOne * Grad(p)
                    PDivV = rShapeDeriv(i, m) * rShapeFunc[j]; // Div(v) * p

                    // Write v * Grad(p) component
                    rDampingMatrix(FirstRow + m, FirstCol + TDim) += Weight * (G - PDivV);
                    // Use symmetry to write the q * Div(u) component
                    rDampingMatrix(FirstCol + TDim, FirstRow + m) += Weight * (G + PDivV);

                    // q-p stabilization block
                    L += rShapeDeriv(i, m) * rShapeDeriv(j, m); // Stabilization: Grad(q) * TauOne * Grad(p)

                    for (unsigned int n = 0; n < TDim; ++n) // iterate over u components (ux,uy[,uz])
                    {
                        // Velocity block
                        rDampingMatrix(FirstRow + m, FirstCol + n) += Weight * TauTwo * rShapeDeriv(i, m) * rShapeDeriv(j, n); // Stabilization: Div(v) * TauTwo * Div(u)
                    }

                }

                // Write remaining terms to velocity block
                for (unsigned int d = 0; d < TDim; ++d)
                    rDampingMatrix(FirstRow + d, FirstCol + d) += K;

                // Write q-p stabilization block
                rDampingMatrix(FirstRow + TDim, FirstCol + TDim) += Weight * TauOne * L;


                // Update reference column index for next iteration
                FirstCol += BlockSize;
            }

            // Operate on RHS
            qF = 0.0;
            for (unsigned int d = 0; d < TDim; ++d)
            {
                rDampRHS[FirstRow + d] += Weight * TauOne * Density * AGradN[i] * BodyForce[d]; // ( a * Grad(v) ) * TauOne * (Density * BodyForce)
                qF += rShapeDeriv(i, d) * BodyForce[d];
            }
            rDampRHS[FirstRow + TDim] += Weight * TauOne * qF; // Grad(q) * TauOne * (Density * BodyForce)

            // Update reference indices
            FirstRow += BlockSize;
            FirstCol = 0;
        }

//            this->AddBTransCB(rDampingMatrix,rShapeDeriv,Viscosity*Weight);
        this->AddViscousTerm(rDampingMatrix,rShapeDeriv,Viscosity*Weight);
    }



    /// Add the surface tension term to the velocity contribution
    // ApplySurfaceTensionContribution(rDampingMatrix, rRightHandSideVector, node_indx, rCurrentProcessInfo);

    void ApplySurfaceTensionContribution(MatrixType& rDampingMatrix, VectorType& rRightHandSideVector,
            const array_1d< double, 3 >& node_indx, const int& k, const ProcessInfo& rCurrentProcessInfo)
    {
	const double gamma = rCurrentProcessInfo[SURFACE_TENSION_COEF]; //surface tension coefficient between air and water [N m-1]

	//adding dissipative_forces varialbes
        double zeta_dissapative_JM = rCurrentProcessInfo[DISSIPATIVE_FORCE_COEFF_JM];
        double zeta_dissapative_BM = rCurrentProcessInfo[DISSIPATIVE_FORCE_COEFF_BM];
        double zeta_dissapative_SM = rCurrentProcessInfo[DISSIPATIVE_FORCE_COEFF_SM];
// 	double gamma_sl = rCurrentProcessInfo[SOLID_LIQIUD_SURFTENS_COEFF];
// 	double gamma_sv = rCurrentProcessInfo[SOLID_AIR_SURFTENS_COEFF];

	double dt = rCurrentProcessInfo[DELTA_TIME];

	double theta_s = rCurrentProcessInfo[CONTACT_ANGLE_STATIC];
	double pi = 3.14159265359;
	theta_s = theta_s*pi/180.0;
	double sin_t = sin(theta_s);
	double cos_t = cos(theta_s);

	array_1d<double,2> m;

	//Clearing spurious curvatures:
	for(unsigned int i = 0; i < 3; i++){
	    if((this->GetGeometry()[i].FastGetSolutionStepValue(IS_BOUNDARY) == 0.0) || (this->GetGeometry()[i].FastGetSolutionStepValue(IS_FREE_SURFACE) == 0.0))
	        this->GetGeometry()[i].FastGetSolutionStepValue(MEAN_CURVATURE_2D) = 0.0;
	}

	//Flag counter to identify contact element:
	double flag_surf = 0.0;
	flag_surf += this->GetGeometry()[node_indx[0]].FastGetSolutionStepValue(IS_FREE_SURFACE);
	flag_surf += this->GetGeometry()[node_indx[1]].FastGetSolutionStepValue(IS_FREE_SURFACE);
	flag_surf += this->GetGeometry()[node_indx[2]].FastGetSolutionStepValue(IS_FREE_SURFACE);
	double flag_trip = 0.0;
	flag_trip += (this->GetGeometry()[node_indx[0]].FastGetSolutionStepValue(TRIPLE_POINT) != 0.0);
	flag_trip += (this->GetGeometry()[node_indx[1]].FastGetSolutionStepValue(TRIPLE_POINT) != 0.0);
	flag_trip += (this->GetGeometry()[node_indx[2]].FastGetSolutionStepValue(TRIPLE_POINT) != 0.0);
	double flag_struct = 0.0;
	flag_struct += (this->GetGeometry()[node_indx[0]].FastGetSolutionStepValue(IS_STRUCTURE) != 0.0);
	flag_struct += (this->GetGeometry()[node_indx[1]].FastGetSolutionStepValue(IS_STRUCTURE) != 0.0);
	flag_struct += (this->GetGeometry()[node_indx[2]].FastGetSolutionStepValue(IS_STRUCTURE) != 0.0);

	int ii = 5;
	int jj = 6;
	int kk = 7;
	double avg_curv = 0.0;

	//////////////////////////////////////////////////////////////////////////////////////////
	//Set the indexes as follows:
	// node "ii" -> triple point, if the element has one
	// node "jj" -> node with flag IS_FREE_SURFACE
	// node "kk" -> in elements with 3 nodes at the boundary:
	//		- if "ii" is TRIPLE_POINT, "kk" has flag IS_STRUCTURE (besides IS_FREE_SURFACE)
	//		- if there is no TRIPLE_POINT, "kk" is another IS_FREE_SURFACE node
	//////////////////////////////////////////////////////////////////////////////////////////
	if(k < 3) //General element with two nodes at the interface
	{
	    //Step to detect triple point. Node with index ii is TRIPLE_POINT, and node with index jj is IS_FREE_SURFACE
	    ii = node_indx[0];
	    jj = node_indx[1];

	    if(flag_trip > 0 && (this->GetGeometry()[node_indx[0]].FastGetSolutionStepValue(TRIPLE_POINT))*1000 == 0.0)
//             if(flag_trip > 0 && (this->GetGeometry()[node_indx[0]].FastGetSolutionStepValue(TRIPLE_POINT)) < 1e-15)
	    {
		  ii = node_indx[1];
		  jj = node_indx[0];
	    }
	}
	else //Element with three nodes at the free surface OR one at free surface, one triple point and one at the structure
	{
	  if(flag_trip == 0.0) //three nodes at interface
	  {
	    if(flag_struct == 0.0)
	    {
	      for(int i = 0; i < 3; i++)
	      {
		avg_curv += 0.333333333333*(this->GetGeometry()[i].FastGetSolutionStepValue(MEAN_CURVATURE_2D));
	      }
	      if(this->GetGeometry()[node_indx[0]].FastGetSolutionStepValue(MEAN_CURVATURE_2D) > avg_curv)
	      {
		ii = node_indx[1];
		jj = node_indx[2];
		kk = node_indx[0];
	      }
	      if(this->GetGeometry()[node_indx[1]].FastGetSolutionStepValue(MEAN_CURVATURE_2D) > avg_curv)
	      {
		ii = node_indx[0];
		jj = node_indx[2];
		kk = node_indx[1];
	      }
	      if(this->GetGeometry()[node_indx[2]].FastGetSolutionStepValue(MEAN_CURVATURE_2D) > avg_curv)
	      {
		ii = node_indx[0];
		jj = node_indx[1];
		kk = node_indx[2];
	      }
	    }
	    else //first time step, TRIPLE_POINT has not been set yet, but the element has a TRIPLE_POINT
	    {
	      if(this->GetGeometry()[node_indx[0]].FastGetSolutionStepValue(IS_FREE_SURFACE) != 0.0)
	      {
		jj = node_indx[0];
		if(this->GetGeometry()[node_indx[1]].FastGetSolutionStepValue(MEAN_CURVATURE_2D) > 1.0)
		{
		  ii = node_indx[1]; //TRIPLE_POINT
		  kk = node_indx[2]; //IS_STRUCTURE
		}
		else
		{
		  ii = node_indx[2]; //TRIPLE_POINT
		  kk = node_indx[1]; //IS_STRUCTURE
		}
	      }
	      if(this->GetGeometry()[node_indx[1]].FastGetSolutionStepValue(IS_FREE_SURFACE) != 0.0)
	      {
		jj = node_indx[1];
		if(this->GetGeometry()[node_indx[0]].FastGetSolutionStepValue(MEAN_CURVATURE_2D) > 1.0)
		{
		  ii = node_indx[0]; //TRIPLE_POINT
		  kk = node_indx[2]; //IS_STRUCTURE
		}
		else
		{
		  ii = node_indx[2]; //TRIPLE_POINT
		  kk = node_indx[0]; //IS_STRUCTURE
		}
	      }
	      if(this->GetGeometry()[node_indx[2]].FastGetSolutionStepValue(IS_FREE_SURFACE) != 0.0)
	      {
		jj = node_indx[2];
		if(this->GetGeometry()[node_indx[0]].FastGetSolutionStepValue(MEAN_CURVATURE_2D) > 1.0)
		{
		  ii = node_indx[0]; //TRIPLE_POINT
		  kk = node_indx[1]; //IS_STRUCTURE
		}
		else
		{
		  ii = node_indx[1]; //TRIPLE_POINT
		  kk = node_indx[0]; //IS_STRUCTURE
		}
	      }
	    }
	  }
	  else //Element has one node with TRIPLE_POINT
	  {
	    if(this->GetGeometry()[node_indx[0]].FastGetSolutionStepValue(TRIPLE_POINT) != 0.0)
	    {
	      ii = node_indx[0];
	      if(this->GetGeometry()[node_indx[1]].FastGetSolutionStepValue(IS_FREE_SURFACE) != 0.0)
	      {
		jj = node_indx[1];
		kk = node_indx[2]; //IS_STRUCTURE
	      }
	      else
	      {
		jj = node_indx[2];
		kk = node_indx[1]; //IS_STRUCTURE
	      }
	    }
	    if(this->GetGeometry()[node_indx[1]].FastGetSolutionStepValue(TRIPLE_POINT) != 0.0)
	    {
	      ii = node_indx[1];
	      if(this->GetGeometry()[node_indx[0]].FastGetSolutionStepValue(IS_FREE_SURFACE) != 0.0)
	      {
		jj = node_indx[0];
		kk = node_indx[2]; //IS_STRUCTURE
	      }
	      else
	      {
		jj = node_indx[2];
		kk = node_indx[0]; //IS_STRUCTURE
	      }
	    }
	    if(this->GetGeometry()[node_indx[2]].FastGetSolutionStepValue(TRIPLE_POINT) != 0.0)
	    {
	      ii = node_indx[2];
	      if(this->GetGeometry()[node_indx[0]].FastGetSolutionStepValue(IS_FREE_SURFACE) != 0.0)
	      {
		jj = node_indx[0];
		kk = node_indx[1];
	      }
	      else
	      {
		jj = node_indx[1];
		kk = node_indx[0];
	      }
	    }
	  }
	}

	array_1d<double,2> An1;
	array_1d<double,2> An2;
	An1[0] = this->GetGeometry()[ii].FastGetSolutionStepValue(NORMAL_X);
	An1[1] = this->GetGeometry()[ii].FastGetSolutionStepValue(NORMAL_Y);
	double norm1 = sqrt(An1[0]*An1[0] + An1[1]*An1[1]);
	An1 /= norm1;
	An2[0] = this->GetGeometry()[jj].FastGetSolutionStepValue(NORMAL_X);
	An2[1] = this->GetGeometry()[jj].FastGetSolutionStepValue(NORMAL_Y);
	double norm2 = sqrt(An2[0]*An2[0] + An2[1]*An2[1]);
	An2 /= norm2;

	double x1 = this->GetGeometry()[ii].X();
	double y1 = this->GetGeometry()[ii].Y();
	double x2 = this->GetGeometry()[jj].X();
	double y2 = this->GetGeometry()[jj].Y();
	//vector x12 is the vector pointing from node 1 [ii] to node 2 [jj]
	array_1d<double,2> x12;
	x12[0] = x2 - x1;
	x12[1] = y2 - y1;
	double dl = sqrt(x12[0]*x12[0] + x12[1]*x12[1]);
	x12 /= dl;

	double curv1 = this->GetGeometry()[ii].FastGetSolutionStepValue(MEAN_CURVATURE_2D);
	double curv2 = this->GetGeometry()[jj].FastGetSolutionStepValue(MEAN_CURVATURE_2D);

	array_1d<double,2> norm_eq;

	//elemental variables for the laplacian
	BoundedMatrix<double,3,3> msWorkMatrix;
        msWorkMatrix = ZeroMatrix(3,3);
	BoundedMatrix<double,3,2> msDN_Dx;
        msDN_Dx = ZeroMatrix(3,2);
	array_1d<double,3> msN; //dimension = number of nodes
        msN = ZeroVector(3);
	double Area;
	GeometryUtils::CalculateGeometryData(this->GetGeometry(),msDN_Dx,msN,Area);
	// 3 by 3 matrix that stores the laplacian
	msWorkMatrix = 0.5 * gamma * dl * prod(msDN_Dx,trans(msDN_Dx)) * dt;



	if(flag_trip == 0)
	{
	  if(this->GetGeometry()[ii].FastGetSolutionStepValue(IS_STRUCTURE) == 0.0)
	  {
	    //CSF:
	    rRightHandSideVector[3*ii]   -= 0.5*gamma*curv1*An1[0]*dl;
	    rRightHandSideVector[3*ii+1] -= 0.5*gamma*curv1*An1[1]*dl;

	    //Output force:
	    this->GetGeometry()[ii].FastGetSolutionStepValue(FORCE_X) = -gamma*curv1*An1[0]*dl;
	    this->GetGeometry()[ii].FastGetSolutionStepValue(FORCE_Y) = -gamma*curv1*An1[1]*dl;
	  }
	}
	else
	{
	    if (this->GetGeometry()[ii].FastGetSolutionStepValue(NORMAL_X) > 0.0)
	    {
		m[0] = -cos_t;
		m[1] = sin_t;
	    }
	    else
	    {
		m[0] = cos_t;
		m[1] = sin_t;
	    }
	    if(k < 3)
	    {

	      if (this->GetGeometry()[ii].FastGetSolutionStepValue(TRIPLE_POINT) != 0.0)
	      {
		  double coef = 1.0;
		  //MODEL 1 - contact angle condition with vector tangent to the surface:
		  rRightHandSideVector[3*ii] 	-= coef*gamma*(m[0]-x12[0]);
 		  rRightHandSideVector[3*ii+1]	-= coef*gamma*(m[1]-x12[1]);
 		  this->GetGeometry()[ii].FastGetSolutionStepValue(FORCE_X) = -coef*gamma*(m[0] - x12[0]);

//  		  this->GetGeometry()[ii].FastGetSolutionStepValue(FORCE_Y) = -coef*gamma*(m[1] - x12[1]);

		  //start of adding dissipative force: where v_clx here is the x_velocity at the contact line; and v_cly is the y_velocity at the contact line
                  double v_clx = this->GetGeometry()[ii].FastGetSolutionStepValue(VELOCITY_X);
                  double v_cly = this->GetGeometry()[ii].FastGetSolutionStepValue(VELOCITY_Y);
//
                  double mu;
                  mu  = this->GetGeometry()[ii].FastGetSolutionStepValue(VISCOSITY);

                  //capillary
		  double v_clx_abs = fabs (v_clx);
		  double v_cly_abs = fabs (v_cly) ;

                  double cap_x   =  mu *  v_clx_abs / gamma;
                  double cap_y   =  mu *  v_cly_abs / gamma;
                  // using Jiang's Model : gamma tanh(4.96 Ca^(0.702))
                  rRightHandSideVector[3*ii]	        -= zeta_dissapative_JM*gamma*tanh(4.96 * pow(cap_x,0.702));
                  rRightHandSideVector[3*ii+1]	        -= zeta_dissapative_JM*gamma*tanh(4.96 * pow(cap_y,0.702));
//
//                   // using Bracke's model : gamma 2.24 ca ^(0.54)
                  rRightHandSideVector[3*ii]	        -= zeta_dissapative_BM*gamma* 2.24 * pow(cap_x,0.54);
                  rRightHandSideVector[3*ii+1]	        -= zeta_dissapative_BM*gamma* 2.24 * pow(cap_y,0.54);
//
//                   // using Seeberg's model : gamm 2.24 ca ^(0.54) for Ca > 10^(-3), otherwise, 4.47 Ca^(0.42)
                  double cap = sqrt((cap_x * cap_x) + (cap_y * cap_y));
                  if (cap > 0.01)
                  {
                    rRightHandSideVector[3*ii]	        -= zeta_dissapative_SM*gamma* 2.24 * pow(cap_x,0.54);
                    rRightHandSideVector[3*ii+1]        -= zeta_dissapative_SM*gamma* 2.24 * pow(cap_y,0.54);
                  }
                  else
                  {
                    rRightHandSideVector[3*ii]	        -= zeta_dissapative_SM*gamma* 4.47 * pow(cap_x,0.42);
                    rRightHandSideVector[3*ii+1]        -= zeta_dissapative_SM*gamma* 4.47 * pow(cap_y,0.42);
                  }
                  // end of adding disppative force

	      }
	      else
	      {
		  this->GetGeometry()[ii].FastGetSolutionStepValue(VELOCITY_X) = 0.0;
	      }
	    }
	}

	//CSF:
	rRightHandSideVector[3*jj] 	-= 0.5*gamma*curv2*An2[0]*dl;
	rRightHandSideVector[3*jj+1] 	-= 0.5*gamma*curv2*An2[1]*dl;

	//Output force:
	this->GetGeometry()[jj].FastGetSolutionStepValue(FORCE_X) = -gamma*curv2*An2[0]*dl;
	this->GetGeometry()[jj].FastGetSolutionStepValue(FORCE_Y) = -gamma*curv2*An2[1]*dl;

	rDampingMatrix(3*ii,3*ii) += msWorkMatrix(ii,ii);
	rDampingMatrix(3*ii+1,3*ii+1) += msWorkMatrix(ii,ii);

	rDampingMatrix(3*ii,3*jj) += msWorkMatrix(ii,jj);
	rDampingMatrix(3*ii+1,3*jj+1) += msWorkMatrix(ii,jj);
//
	rDampingMatrix(3*jj,3*ii) += msWorkMatrix(jj,ii);
	rDampingMatrix(3*jj+1,3*ii+1) += msWorkMatrix(jj,ii);

	rDampingMatrix(3*jj,3*jj) += msWorkMatrix(jj,jj);
	rDampingMatrix(3*jj+1,3*jj+1) += msWorkMatrix(jj,jj);


	if(k > 2 && this->GetGeometry()[kk].FastGetSolutionStepValue(IS_STRUCTURE) == 0.0)
	{
	    array_1d<double,2> An3;
	    An3[0] = this->GetGeometry()[kk].FastGetSolutionStepValue(NORMAL_X);
	    An3[1] = this->GetGeometry()[kk].FastGetSolutionStepValue(NORMAL_Y);
	    double norm3 = sqrt(An3[0]*An3[0] + An3[1]*An3[1]);
	    An3 /= norm3;
	    double curv3 = this->GetGeometry()[kk].FastGetSolutionStepValue(MEAN_CURVATURE_2D);

	    double x3 = this->GetGeometry()[kk].X();
	    double y3 = this->GetGeometry()[kk].Y();
	    array_1d<double,2> x31;
	    array_1d<double,2> x32;
	    x31[0] = x1 - x3;
	    x31[1] = y1 - y3;
	    double dl1 = sqrt(x31[0]*x31[0] + x31[1]*x31[1]);
	    x32[0] = x2 - x3;
	    x32[1] = y2 - y3;
	    double dl2 = sqrt(x32[0]*x32[0] + x32[1]*x32[1]);
	    dl = dl1 + dl2;

	    //CSF:
	    rRightHandSideVector[3*kk] 		-= gamma*curv3*An3[0]*dl;
	    rRightHandSideVector[3*kk+1] 	-= gamma*curv3*An3[1]*dl;
	    msWorkMatrix = 1.0 * gamma * dl * prod(msDN_Dx,trans(msDN_Dx)) * dt;


//
	    rDampingMatrix(3*kk,3*kk) += msWorkMatrix(kk,kk);
	    rDampingMatrix(3*kk+1,3*kk+1) += msWorkMatrix(kk,kk);



	}

	//Clean spurious force values
	for(unsigned int i = 0; i < 3; i++)
	{
	  if(this->GetGeometry()[i].FastGetSolutionStepValue(IS_BOUNDARY) == 0.0)
	  {
	    this->GetGeometry()[i].FastGetSolutionStepValue(FORCE_X) = 0.0;
	    this->GetGeometry()[i].FastGetSolutionStepValue(FORCE_Y) = 0.0;
	    this->GetGeometry()[i].FastGetSolutionStepValue(FORCE_Z) = 0.0;
	  }
	}
    }


    void ApplySurfaceTensionContribution3D(MatrixType& rDampingMatrix, VectorType& rRightHandSideVector,
            const array_1d< double, 4 >& node_indx, const int& k, const ProcessInfo& rCurrentProcessInfo)
    {
	double dt = rCurrentProcessInfo[DELTA_TIME];
	double gamma = rCurrentProcessInfo[SURFACE_TENSION_COEF];
	double theta_s = rCurrentProcessInfo[CONTACT_ANGLE_STATIC];

	//Flag counter to identify contact element:
	double flag_surf = 0.0;
	flag_surf += this->GetGeometry()[node_indx[0]].FastGetSolutionStepValue(IS_FREE_SURFACE);
	flag_surf += this->GetGeometry()[node_indx[1]].FastGetSolutionStepValue(IS_FREE_SURFACE);
	flag_surf += this->GetGeometry()[node_indx[2]].FastGetSolutionStepValue(IS_FREE_SURFACE);
	flag_surf += this->GetGeometry()[node_indx[3]].FastGetSolutionStepValue(IS_FREE_SURFACE);
	double flag_trip = 0.0;
	flag_trip += (this->GetGeometry()[node_indx[0]].FastGetSolutionStepValue(TRIPLE_POINT) != 0.0);
	flag_trip += (this->GetGeometry()[node_indx[1]].FastGetSolutionStepValue(TRIPLE_POINT) != 0.0);
	flag_trip += (this->GetGeometry()[node_indx[2]].FastGetSolutionStepValue(TRIPLE_POINT) != 0.0);
	flag_trip += (this->GetGeometry()[node_indx[3]].FastGetSolutionStepValue(TRIPLE_POINT) != 0.0);
	double flag_struct = 0.0;
	flag_struct += (this->GetGeometry()[node_indx[0]].FastGetSolutionStepValue(IS_STRUCTURE) != 0.0);
	flag_struct += (this->GetGeometry()[node_indx[1]].FastGetSolutionStepValue(IS_STRUCTURE) != 0.0);
	flag_struct += (this->GetGeometry()[node_indx[2]].FastGetSolutionStepValue(IS_STRUCTURE) != 0.0);
	flag_struct += (this->GetGeometry()[node_indx[3]].FastGetSolutionStepValue(IS_STRUCTURE) != 0.0);

	int ii = 5;
	int jj = 6;
	int kk = 7;
	int ll = 8;
	double avg_curv = 0.0;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Set the indexes as follows:
	// node "ii" and "jj" -> triple point, if the element has (at least) two (those with one are not taken into account)
	// node "kk" -> node with flag either TRIPLE_POINT (corner element) or IS_FREE_SURFACE
	// node "ll" -> in elements with 4 nodes at the boundary:
	//		- if "ii" and "jj" are TRIPLE_POINT, "ll" has flag IS_STRUCTURE (besides IS_FREE_SURFACE)
	//		- if there is no TRIPLE_POINT, "ll" is another IS_FREE_SURFACE node
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	if(k == 3) //General element with three nodes at the interface
	{
	    //Step to detect triple point. Nodes with index ii and jj are TRIPLE_POINT, and node with index kk is IS_FREE_SURFACE
	    ii = node_indx[0];
	    jj = node_indx[1];
	    kk = node_indx[2];

	    if(flag_trip == 1)
	    {
	      if ((this->GetGeometry()[node_indx[1]].FastGetSolutionStepValue(TRIPLE_POINT))*1000 != 0.0)
//               if ((this->GetGeometry()[node_indx[1]].FastGetSolutionStepValue(TRIPLE_POINT)) < 1e-15)
	      {
		  ii = node_indx[1];
		  jj = node_indx[0];
	      }
	      if ((this->GetGeometry()[node_indx[2]].FastGetSolutionStepValue(TRIPLE_POINT))*1000 != 0.0)
//               if ((this->GetGeometry()[node_indx[2]].FastGetSolutionStepValue(TRIPLE_POINT)) < 1e-15 )
	      {
		  ii = node_indx[2];
		  kk = node_indx[0];
	      }
	    }
	    if(flag_trip > 1)
	    {
	      if ((this->GetGeometry()[node_indx[0]].FastGetSolutionStepValue(TRIPLE_POINT))*1000 == 0.0)
//               if ((this->GetGeometry()[node_indx[0]].FastGetSolutionStepValue(TRIPLE_POINT))  < 1e-15 )
	      {
		  ii = node_indx[1];
		  jj = node_indx[2];
		  kk = node_indx[0];
	      }
	      if ((this->GetGeometry()[node_indx[1]].FastGetSolutionStepValue(TRIPLE_POINT))*1000 == 0.0)
//               if ((this->GetGeometry()[node_indx[1]].FastGetSolutionStepValue(TRIPLE_POINT)) < 1e-15 )
	      {
		  jj = node_indx[2];
		  kk = node_indx[1];
	      }
	    }
	}
	else //Element with four nodes at the free surface OR one at free surface, two triple point and one at the structure OR one at free surface and three triple points
	{
	  ii = node_indx[0];
	  jj = node_indx[1];
	  kk = node_indx[2];
	  ll = node_indx[3];
	  if(flag_trip == 0.0) //four nodes at interface
//           if(flag_trip < 1e-15)
	  {
	    if(flag_struct == 0.0) //four nodes that are free surface
//             if(flag_struct < 1e-15)
	    {
	      for(int i = 0; i < 4; i++)
	      {
		avg_curv += 0.25*(this->GetGeometry()[i].FastGetSolutionStepValue(MEAN_CURVATURE_3D));
	      }
	      if(this->GetGeometry()[node_indx[0]].FastGetSolutionStepValue(MEAN_CURVATURE_3D) > avg_curv)
	      {
		ii = node_indx[0];
		jj = node_indx[1];
		kk = node_indx[2];
		ll = node_indx[3];
	      }
	      if(this->GetGeometry()[node_indx[1]].FastGetSolutionStepValue(MEAN_CURVATURE_3D) > avg_curv)
	      {
		ii = node_indx[1];
		jj = node_indx[0];
		kk = node_indx[2];
		ll = node_indx[3];
	      }
	      if(this->GetGeometry()[node_indx[2]].FastGetSolutionStepValue(MEAN_CURVATURE_3D) > avg_curv)
	      {
		ii = node_indx[2];
		jj = node_indx[0];
		kk = node_indx[1];
		ll = node_indx[3];
	      }
	      if(this->GetGeometry()[node_indx[3]].FastGetSolutionStepValue(MEAN_CURVATURE_3D) > avg_curv)
	      {
		ii = node_indx[3];
		jj = node_indx[0];
		kk = node_indx[1];
		ll = node_indx[2];
	      }
	    }
	    else //first time step, TRIPLE_POINT has not been set yet, but the element has three IS_STRUCTURE nodes
	    {
	      //OPTION 1
	      if(this->GetGeometry()[node_indx[0]].FastGetSolutionStepValue(IS_FREE_SURFACE) != 0.0)
	      {
		kk = node_indx[0];
		if(this->GetGeometry()[node_indx[1]].FastGetSolutionStepValue(MEAN_CURVATURE_3D) > 1.0)
		{
		  ii = node_indx[1]; //TRIPLE_POINT
		  if(this->GetGeometry()[node_indx[2]].FastGetSolutionStepValue(MEAN_CURVATURE_3D) > 1.0)
		  {
		    jj = node_indx[2]; //TRIPLE_POINT
		    ll = node_indx[3]; //IS_STRUCTURE
		  }
		  if(this->GetGeometry()[node_indx[3]].FastGetSolutionStepValue(MEAN_CURVATURE_3D) > 1.0)
		  {
		    jj = node_indx[3]; //TRIPLE_POINT
		    ll = node_indx[2]; //IS_STRUCTURE
		  }
		}
		if(this->GetGeometry()[node_indx[2]].FastGetSolutionStepValue(MEAN_CURVATURE_3D) > 1.0)
		{
		  ii = node_indx[2]; //TRIPLE_POINT
		  if(this->GetGeometry()[node_indx[1]].FastGetSolutionStepValue(MEAN_CURVATURE_3D) > 1.0)
		  {
		    jj = node_indx[1]; //TRIPLE_POINT
		    ll = node_indx[3]; //IS_STRUCTURE
		  }
		  if(this->GetGeometry()[node_indx[3]].FastGetSolutionStepValue(MEAN_CURVATURE_3D) > 1.0)
		  {
		    jj = node_indx[3]; //TRIPLE_POINT
		    ll = node_indx[1]; //IS_STRUCTURE
		  }
		}
		if(this->GetGeometry()[node_indx[3]].FastGetSolutionStepValue(MEAN_CURVATURE_3D) > 1.0)
		{
		  ii = node_indx[3]; //TRIPLE_POINT
		  if(this->GetGeometry()[node_indx[1]].FastGetSolutionStepValue(MEAN_CURVATURE_3D) > 1.0)
		  {
		    jj = node_indx[1]; //TRIPLE_POINT
		    ll = node_indx[2]; //IS_STRUCTURE
		  }
		  if(this->GetGeometry()[node_indx[2]].FastGetSolutionStepValue(MEAN_CURVATURE_3D) > 1.0)
		  {
		    jj = node_indx[2]; //TRIPLE_POINT
		    ll = node_indx[1]; //IS_STRUCTURE
		  }
		}
	      }

	      //OPTION 2
	      if(this->GetGeometry()[node_indx[1]].FastGetSolutionStepValue(IS_FREE_SURFACE) != 0.0)
	      {
		kk = node_indx[1];
		if(this->GetGeometry()[node_indx[0]].FastGetSolutionStepValue(MEAN_CURVATURE_3D) > 1.0)
		{
		  ii = node_indx[0]; //TRIPLE_POINT
		  if(this->GetGeometry()[node_indx[2]].FastGetSolutionStepValue(MEAN_CURVATURE_3D) > 1.0)
		  {
		    jj = node_indx[2]; //TRIPLE_POINT
		    ll = node_indx[3]; //IS_STRUCTURE
		  }
		  if(this->GetGeometry()[node_indx[3]].FastGetSolutionStepValue(MEAN_CURVATURE_3D) > 1.0)
		  {
		    jj = node_indx[3]; //TRIPLE_POINT
		    ll = node_indx[2]; //IS_STRUCTURE
		  }
		}
		if(this->GetGeometry()[node_indx[2]].FastGetSolutionStepValue(MEAN_CURVATURE_3D) > 1.0)
		{
		  ii = node_indx[2]; //TRIPLE_POINT
		  if(this->GetGeometry()[node_indx[0]].FastGetSolutionStepValue(MEAN_CURVATURE_3D) > 1.0)
		  {
		    jj = node_indx[0]; //TRIPLE_POINT
		    ll = node_indx[3]; //IS_STRUCTURE
		  }
		  if(this->GetGeometry()[node_indx[3]].FastGetSolutionStepValue(MEAN_CURVATURE_3D) > 1.0)
		  {
		    jj = node_indx[3]; //TRIPLE_POINT
		    ll = node_indx[0]; //IS_STRUCTURE
		  }
		}
		if(this->GetGeometry()[node_indx[3]].FastGetSolutionStepValue(MEAN_CURVATURE_3D) > 1.0)
		{
		  ii = node_indx[3]; //TRIPLE_POINT
		  if(this->GetGeometry()[node_indx[0]].FastGetSolutionStepValue(MEAN_CURVATURE_3D) > 1.0)
		  {
		    jj = node_indx[0]; //TRIPLE_POINT
		    ll = node_indx[2]; //IS_STRUCTURE
		  }
		  if(this->GetGeometry()[node_indx[2]].FastGetSolutionStepValue(MEAN_CURVATURE_3D) > 1.0)
		  {
		    jj = node_indx[2]; //TRIPLE_POINT
		    ll = node_indx[0]; //IS_STRUCTURE
		  }
		}
	      }

	      //OPTION 3
	      if(this->GetGeometry()[node_indx[2]].FastGetSolutionStepValue(IS_FREE_SURFACE) != 0.0)
	      {
		kk = node_indx[2];
		if(this->GetGeometry()[node_indx[1]].FastGetSolutionStepValue(MEAN_CURVATURE_3D) > 1.0)
		{
		  ii = node_indx[1]; //TRIPLE_POINT
		  if(this->GetGeometry()[node_indx[0]].FastGetSolutionStepValue(MEAN_CURVATURE_3D) > 1.0)
		  {
		    jj = node_indx[0]; //TRIPLE_POINT
		    ll = node_indx[3]; //IS_STRUCTURE
		  }
		  if(this->GetGeometry()[node_indx[3]].FastGetSolutionStepValue(MEAN_CURVATURE_3D) > 1.0)
		  {
		    jj = node_indx[3]; //TRIPLE_POINT
		    ll = node_indx[0]; //IS_STRUCTURE
		  }
		}
		if(this->GetGeometry()[node_indx[0]].FastGetSolutionStepValue(MEAN_CURVATURE_3D) > 1.0)
		{
		  ii = node_indx[0]; //TRIPLE_POINT
		  if(this->GetGeometry()[node_indx[1]].FastGetSolutionStepValue(MEAN_CURVATURE_3D) > 1.0)
		  {
		    jj = node_indx[1]; //TRIPLE_POINT
		    ll = node_indx[3]; //IS_STRUCTURE
		  }
		  if(this->GetGeometry()[node_indx[3]].FastGetSolutionStepValue(MEAN_CURVATURE_3D) > 1.0)
		  {
		    jj = node_indx[3]; //TRIPLE_POINT
		    ll = node_indx[1]; //IS_STRUCTURE
		  }
		}
		if(this->GetGeometry()[node_indx[3]].FastGetSolutionStepValue(MEAN_CURVATURE_3D) > 1.0)
		{
		  ii = node_indx[3]; //TRIPLE_POINT
		  if(this->GetGeometry()[node_indx[1]].FastGetSolutionStepValue(MEAN_CURVATURE_3D) > 1.0)
		  {
		    jj = node_indx[1]; //TRIPLE_POINT
		    ll = node_indx[0]; //IS_STRUCTURE
		  }
		  if(this->GetGeometry()[node_indx[0]].FastGetSolutionStepValue(MEAN_CURVATURE_3D) > 1.0)
		  {
		    jj = node_indx[0]; //TRIPLE_POINT
		    ll = node_indx[1]; //IS_STRUCTURE
		  }
		}
	      }

	      //OPTION 4
	      if(this->GetGeometry()[node_indx[3]].FastGetSolutionStepValue(IS_FREE_SURFACE) != 0.0)
	      {
		kk = node_indx[3];
		if(this->GetGeometry()[node_indx[1]].FastGetSolutionStepValue(MEAN_CURVATURE_3D) > 1.0)
		{
		  ii = node_indx[1]; //TRIPLE_POINT
		  if(this->GetGeometry()[node_indx[2]].FastGetSolutionStepValue(MEAN_CURVATURE_3D) > 1.0)
		  {
		    jj = node_indx[2]; //TRIPLE_POINT
		    ll = node_indx[0]; //IS_STRUCTURE
		  }
		  if(this->GetGeometry()[node_indx[0]].FastGetSolutionStepValue(MEAN_CURVATURE_3D) > 1.0)
		  {
		    jj = node_indx[0]; //TRIPLE_POINT
		    ll = node_indx[2]; //IS_STRUCTURE
		  }
		}
		if(this->GetGeometry()[node_indx[2]].FastGetSolutionStepValue(MEAN_CURVATURE_3D) > 1.0)
		{
		  ii = node_indx[2]; //TRIPLE_POINT
		  if(this->GetGeometry()[node_indx[1]].FastGetSolutionStepValue(MEAN_CURVATURE_3D) > 1.0)
		  {
		    jj = node_indx[1]; //TRIPLE_POINT
		    ll = node_indx[0]; //IS_STRUCTURE
		  }
		  if(this->GetGeometry()[node_indx[0]].FastGetSolutionStepValue(MEAN_CURVATURE_3D) > 1.0)
		  {
		    jj = node_indx[0]; //TRIPLE_POINT
		    ll = node_indx[1]; //IS_STRUCTURE
		  }
		}
		if(this->GetGeometry()[node_indx[0]].FastGetSolutionStepValue(MEAN_CURVATURE_3D) > 1.0)
		{
		  ii = node_indx[0]; //TRIPLE_POINT
		  if(this->GetGeometry()[node_indx[1]].FastGetSolutionStepValue(MEAN_CURVATURE_3D) > 1.0)
		  {
		    jj = node_indx[1]; //TRIPLE_POINT
		    ll = node_indx[2]; //IS_STRUCTURE
		  }
		  if(this->GetGeometry()[node_indx[2]].FastGetSolutionStepValue(MEAN_CURVATURE_3D) > 1.0)
		  {
		    jj = node_indx[2]; //TRIPLE_POINT
		    ll = node_indx[1]; //IS_STRUCTURE
		  }
		}
	      }

	    }
	  }
	  if (flag_trip == 2) //Element has two nodes with TRIPLE_POINT (those with one are not taken into accounts)
	  {
	    //OPTION 1
	    if(this->GetGeometry()[node_indx[0]].FastGetSolutionStepValue(TRIPLE_POINT) != 0.0)
	    {
	      ii = node_indx[0];
	      if(this->GetGeometry()[node_indx[1]].FastGetSolutionStepValue(TRIPLE_POINT) != 0.0)
	      {
		jj = node_indx[1];
		if(this->GetGeometry()[node_indx[2]].FastGetSolutionStepValue(IS_FREE_SURFACE) != 0.0)
		{
		  kk = node_indx[2]; //IS_FREE_SURFACE
		  ll = node_indx[3]; //IS_STRUCTURE
		}
		else
		{
		  kk = node_indx[3]; //IS_FREE_SURFACE
		  ll = node_indx[2]; //IS_STRUCTURE
		}
	      }
	      if(this->GetGeometry()[node_indx[2]].FastGetSolutionStepValue(TRIPLE_POINT) != 0.0)
	      {
		jj = node_indx[2];
		if(this->GetGeometry()[node_indx[1]].FastGetSolutionStepValue(IS_FREE_SURFACE) != 0.0)
		{
		  kk = node_indx[1]; //IS_FREE_SURFACE
		  ll = node_indx[3]; //IS_STRUCTURE
		}
		else
		{
		  kk = node_indx[3]; //IS_FREE_SURFACE
		  ll = node_indx[1]; //IS_STRUCTURE
		}
	      }
	      if(this->GetGeometry()[node_indx[3]].FastGetSolutionStepValue(TRIPLE_POINT) != 0.0)
	      {
		jj = node_indx[3];
		if(this->GetGeometry()[node_indx[1]].FastGetSolutionStepValue(IS_FREE_SURFACE) != 0.0)
		{
		  kk = node_indx[1]; //IS_FREE_SURFACE
		  ll = node_indx[2]; //IS_STRUCTURE
		}
		else
		{
		  kk = node_indx[2]; //IS_FREE_SURFACE
		  ll = node_indx[1]; //IS_STRUCTURE
		}
	      }
	    }
	    //OPTION 2
	    if(this->GetGeometry()[node_indx[1]].FastGetSolutionStepValue(TRIPLE_POINT) != 0.0)
	    {
	      ii = node_indx[1];
	      if(this->GetGeometry()[node_indx[2]].FastGetSolutionStepValue(TRIPLE_POINT) != 0.0)
	      {
		jj = node_indx[2];
		if(this->GetGeometry()[node_indx[0]].FastGetSolutionStepValue(IS_FREE_SURFACE) != 0.0)
		{
		  kk = node_indx[0]; //IS_FREE_SURFACE
		  ll = node_indx[3]; //IS_STRUCTURE
		}
		else
		{
		  kk = node_indx[3]; //IS_FREE_SURFACE
		  ll = node_indx[0]; //IS_STRUCTURE
		}
	      }
	      if(this->GetGeometry()[node_indx[3]].FastGetSolutionStepValue(TRIPLE_POINT) != 0.0)
	      {
		jj = node_indx[3];
		if(this->GetGeometry()[node_indx[0]].FastGetSolutionStepValue(IS_FREE_SURFACE) != 0.0)
		{
		  kk = node_indx[0]; //IS_FREE_SURFACE
		  ll = node_indx[2]; //IS_STRUCTURE
		}
		else
		{
		  kk = node_indx[2]; //IS_FREE_SURFACE
		  ll = node_indx[0]; //IS_STRUCTURE
		}
	      }
	    }
	    //OPTION 3
	    if(this->GetGeometry()[node_indx[2]].FastGetSolutionStepValue(TRIPLE_POINT) != 0.0)
	    {
	      ii = node_indx[2];
	      if(this->GetGeometry()[node_indx[3]].FastGetSolutionStepValue(TRIPLE_POINT) != 0.0)
	      {
		jj = node_indx[3];
		if(this->GetGeometry()[node_indx[0]].FastGetSolutionStepValue(IS_FREE_SURFACE) != 0.0)
		{
		  kk = node_indx[0]; //IS_FREE_SURFACE
		  ll = node_indx[1]; //IS_STRUCTURE
		}
		else
		{
		  kk = node_indx[1]; //IS_FREE_SURFACE
		  ll = node_indx[0]; //IS_STRUCTURE
		}
	      }
	    }
	  }
	  if (flag_trip == 3) //Element has three nodes with TRIPLE_POINT (those with one are not taken into account)
	  {
	    //OPTION 1
	    if(this->GetGeometry()[node_indx[0]].FastGetSolutionStepValue(TRIPLE_POINT) != 0.0)
	    {
	      ii = node_indx[0];
	      if(this->GetGeometry()[node_indx[1]].FastGetSolutionStepValue(TRIPLE_POINT) != 0.0)
	      {
		jj = node_indx[1];
		if(this->GetGeometry()[node_indx[2]].FastGetSolutionStepValue(TRIPLE_POINT) != 0.0)
		{
		  kk = node_indx[2]; //TRIPLE_POINT
		  ll = node_indx[3]; //IS_FREE_SURFACE
		}
		else
		{
		  kk = node_indx[3]; //TRIPLE_POINT
		  ll = node_indx[2]; //IS_FREE_SURFACE
		}
	      }
	      else //if second node is not triple point, the only option is that the rest are triple points and this one is free surface
	      {
		jj = node_indx[2];
		kk = node_indx[3];
		ll = node_indx[1];
	      }
	    }
	    else //if first node is not triple point, the only option is that the rest are triple points and this one is free surface
	    {
	      ii = node_indx[1];
	      jj = node_indx[2];
	      kk = node_indx[3];
	      ll = node_indx[0];
	    }
	  }
	}

	array_1d<double,3> An1 = this->GetGeometry()[ii].FastGetSolutionStepValue(NORMAL_GEOMETRIC);
	array_1d<double,3> An2 = this->GetGeometry()[jj].FastGetSolutionStepValue(NORMAL_GEOMETRIC);
	array_1d<double,3> An3 = this->GetGeometry()[kk].FastGetSolutionStepValue(NORMAL_GEOMETRIC);

	array_1d<double,3> m1 = - this->GetGeometry()[ii].FastGetSolutionStepValue(NORMAL_CONTACT_LINE) + this->GetGeometry()[ii].FastGetSolutionStepValue(NORMAL_CONTACT_LINE_EQUILIBRIUM);
	array_1d<double,3> m2 = - this->GetGeometry()[jj].FastGetSolutionStepValue(NORMAL_CONTACT_LINE) + this->GetGeometry()[jj].FastGetSolutionStepValue(NORMAL_CONTACT_LINE_EQUILIBRIUM);
	array_1d<double,3> m3 = - this->GetGeometry()[kk].FastGetSolutionStepValue(NORMAL_CONTACT_LINE) + this->GetGeometry()[kk].FastGetSolutionStepValue(NORMAL_CONTACT_LINE_EQUILIBRIUM);

	double fsign1 = this->GetGeometry()[ii].FastGetSolutionStepValue(CONTACT_ANGLE) - theta_s;
	fsign1 = 1.0;//fsign1/sqrt(fsign1*fsign1);
	double fsign2 = this->GetGeometry()[jj].FastGetSolutionStepValue(CONTACT_ANGLE) - theta_s;
	fsign2 = 1.0;//fsign2/sqrt(fsign2*fsign2);
	double fsign3 = this->GetGeometry()[kk].FastGetSolutionStepValue(CONTACT_ANGLE) - theta_s;
	fsign3 = 1.0;//fsign3/sqrt(fsign3*fsign3);

	//Check if there is a node with TRIPLE_POINT flag which shouldn't be
	if (flag_trip == 3)
	{
	    double temp = 0.0;
	    array_1d<double,3> vec_temp = ZeroVector(3);
	    for(unsigned int i_node = 0; i_node < 4; ++i_node)
	    {
		vec_temp = this->GetGeometry()[i_node].FastGetSolutionStepValue(NORMAL);
		NormalizeVec3D(vec_temp);
		temp = -(vec_temp[2]);
		if ( temp > 0.99 && this->GetGeometry()[i_node].FastGetSolutionStepValue(TRIPLE_POINT) != 0.0)
		    this->GetGeometry()[i_node].FastGetSolutionStepValue(TRIPLE_POINT) = 0.0;
	    }
	}

	double curv1 = this->GetGeometry()[ii].FastGetSolutionStepValue(MEAN_CURVATURE_3D);
	double curv2 = this->GetGeometry()[jj].FastGetSolutionStepValue(MEAN_CURVATURE_3D);
	double curv3 = this->GetGeometry()[kk].FastGetSolutionStepValue(MEAN_CURVATURE_3D);

	double nlen1 = this->GetGeometry()[ii].FastGetSolutionStepValue(NODAL_LENGTH);
	double nlen2 = this->GetGeometry()[jj].FastGetSolutionStepValue(NODAL_LENGTH);
	double nlen3 = this->GetGeometry()[kk].FastGetSolutionStepValue(NODAL_LENGTH);

	double x1 = this->GetGeometry()[ii].X();
	double y1 = this->GetGeometry()[ii].Y();
	double z1 = this->GetGeometry()[ii].Z();
	double x2 = this->GetGeometry()[jj].X();
	double y2 = this->GetGeometry()[jj].Y();
	double z2 = this->GetGeometry()[jj].Z();
	double x3 = this->GetGeometry()[kk].X();
	double y3 = this->GetGeometry()[kk].Y();
	double z3 = this->GetGeometry()[kk].Z();

	array_1d<double,3> r12 = Vector3D(x1,y1,z1,x2,y2,z2);
	array_1d<double,3> r13 = Vector3D(x1,y1,z1,x3,y3,z3);
	array_1d<double,3> r23 = Vector3D(x2,y2,z2,x3,y3,z3);
	array_1d<double,3> cprod = CrossProduct3D(r12,r13);
	double area_tr = 0.5*Norm3D(cprod);



	//elemental variables for the laplacian
	BoundedMatrix<double,4,4> msWorkMatrix = ZeroMatrix(4,4);
	BoundedMatrix<double,4,3> msDN_Dx = ZeroMatrix(4,3);
	array_1d<double,4> msN = ZeroVector(4); //dimension = number of nodes
	double Volume;
	GeometryUtils::CalculateGeometryData(this->GetGeometry(),msDN_Dx,msN,Volume);
	// 4 by 4 matrix that stores the laplacian
 	msWorkMatrix = 0.333333333333 * gamma * area_tr * prod(msDN_Dx,trans(msDN_Dx)) * dt;

	double coef_i = 0.333333333333; // 0.333333333333 | (1/num_neighs_i) | 1.0/(num_neighs_i-1)
	double coef_j = 0.333333333333; // 0.333333333333 | (1/num_neighs_j) | 1.0/(num_neighs_j-1)
	double coef_k = 0.333333333333; // 0.333333333333 | (1/num_neighs_k) | 1.0/(num_neighs_k-1)

        if(flag_trip == 0)
//         if(flag_trip < 1e-15)
	{
	  rRightHandSideVector[4*ii]   -= coef_i*gamma*curv1*An1[0]*area_tr;
	  rRightHandSideVector[4*ii+1] -= coef_i*gamma*curv1*An1[1]*area_tr;
	  rRightHandSideVector[4*ii+2] -= coef_i*gamma*curv1*An1[2]*area_tr;
	  this->GetGeometry()[ii].FastGetSolutionStepValue(FORCE_X) = -coef_i*gamma*curv1*An1[0]*area_tr;
	  this->GetGeometry()[ii].FastGetSolutionStepValue(FORCE_Y) = -coef_i*gamma*curv1*An1[1]*area_tr;
	  this->GetGeometry()[ii].FastGetSolutionStepValue(FORCE_Z) = -coef_i*gamma*curv1*An1[2]*area_tr;

	  rRightHandSideVector[4*jj]   -= coef_j*gamma*curv2*An2[0]*area_tr;
	  rRightHandSideVector[4*jj+1] -= coef_j*gamma*curv2*An2[1]*area_tr;
	  rRightHandSideVector[4*jj+2] -= coef_j*gamma*curv2*An2[2]*area_tr;
	  this->GetGeometry()[jj].FastGetSolutionStepValue(FORCE_X) = -coef_j*gamma*curv2*An2[0]*area_tr;
	  this->GetGeometry()[jj].FastGetSolutionStepValue(FORCE_Y) = -coef_j*gamma*curv2*An2[1]*area_tr;
	  this->GetGeometry()[jj].FastGetSolutionStepValue(FORCE_Z) = -coef_j*gamma*curv2*An2[2]*area_tr;

	  rRightHandSideVector[4*kk]   -= coef_k*gamma*curv3*An3[0]*area_tr;
	  rRightHandSideVector[4*kk+1] -= coef_k*gamma*curv3*An3[1]*area_tr;
	  rRightHandSideVector[4*kk+2] -= coef_k*gamma*curv3*An3[2]*area_tr;
	  this->GetGeometry()[kk].FastGetSolutionStepValue(FORCE_X) = -coef_k*gamma*curv3*An3[0]*area_tr;
	  this->GetGeometry()[kk].FastGetSolutionStepValue(FORCE_Y) = -coef_k*gamma*curv3*An3[1]*area_tr;
	  this->GetGeometry()[kk].FastGetSolutionStepValue(FORCE_Z) = -coef_k*gamma*curv3*An3[2]*area_tr;
	}

	double beta = 1.0;
	if(flag_trip >= 1)
	{
	  if(this->GetGeometry()[ii].FastGetSolutionStepValue(TRIPLE_POINT) != 0.0)
	  {
	    coef_i = 0.5;
	    rRightHandSideVector[4*ii]   -= coef_i*beta*fsign1*(gamma*m1[0]*nlen1);// + this->GetGeometry()[ii].FastGetSolutionStepValue(ADHESION_FORCE_X));
	    rRightHandSideVector[4*ii+1] -= coef_i*beta*fsign1*(gamma*m1[1]*nlen1);// + this->GetGeometry()[ii].FastGetSolutionStepValue(ADHESION_FORCE_Y));
	    rRightHandSideVector[4*ii+2] -= coef_i*beta*fsign1*(gamma*m1[2]*nlen1);// + this->GetGeometry()[ii].FastGetSolutionStepValue(ADHESION_FORCE_Z));
	    this->GetGeometry()[ii].FastGetSolutionStepValue(FORCE) = -coef_i*beta*fsign1*(gamma*m1*nlen1);
            // + this->GetGeometry()[ii].FastGetSolutionStepValue(ADHESION_FORCE_X));


	    rDampingMatrix(4*ii,4*ii)     += 0.5*gamma*dt*msN[ii]*msN[ii]*nlen1 - msWorkMatrix(ii,ii);
	    rDampingMatrix(4*ii+1,4*ii+1) += 0.5*gamma*dt*msN[ii]*msN[ii]*nlen1 - msWorkMatrix(ii,ii);
	    rDampingMatrix(4*ii+2,4*ii+2) += 0.5*gamma*dt*msN[ii]*msN[ii]*nlen1 - msWorkMatrix(ii,ii);
	    rDampingMatrix(4*ii,4*jj)     += 0.5*gamma*dt*msN[ii]*msN[jj]*0.5*(nlen1 + nlen2) - msWorkMatrix(ii,jj);
	    rDampingMatrix(4*ii+1,4*jj+1) += 0.5*gamma*dt*msN[ii]*msN[jj]*0.5*(nlen1 + nlen2) - msWorkMatrix(ii,jj);
	    rDampingMatrix(4*ii+2,4*jj+2) += 0.5*gamma*dt*msN[ii]*msN[jj]*0.5*(nlen1 + nlen2) - msWorkMatrix(ii,jj);
	  }
	  else
	  {
	    rRightHandSideVector[4*ii]   -= coef_i*beta*gamma*curv1*An1[0]*area_tr;
	    rRightHandSideVector[4*ii+1] -= coef_i*beta*gamma*curv1*An1[1]*area_tr;
	    rRightHandSideVector[4*ii+2] -= coef_i*beta*gamma*curv1*An1[2]*area_tr;
	    this->GetGeometry()[ii].FastGetSolutionStepValue(FORCE_X) = -coef_i*beta*gamma*curv1*An1[0]*area_tr;
	    this->GetGeometry()[ii].FastGetSolutionStepValue(FORCE_Y) = -coef_i*beta*gamma*curv1*An1[1]*area_tr;
	    this->GetGeometry()[ii].FastGetSolutionStepValue(FORCE_Z) = -coef_i*beta*gamma*curv1*An1[2]*area_tr;
	  }

// 	  beta = 1.0;
	  if(this->GetGeometry()[jj].FastGetSolutionStepValue(TRIPLE_POINT) != 0.0)
	  {
	    coef_j = 0.5;
	    rRightHandSideVector[4*jj]   -= coef_j*beta*fsign2*(gamma*m2[0]*nlen2);// + this->GetGeometry()[jj].FastGetSolutionStepValue(ADHESION_FORCE_X));
	    rRightHandSideVector[4*jj+1] -= coef_j*beta*fsign2*(gamma*m2[1]*nlen2);// + this->GetGeometry()[jj].FastGetSolutionStepValue(ADHESION_FORCE_Y));
	    rRightHandSideVector[4*jj+2] -= coef_j*beta*fsign2*(gamma*m2[2]*nlen2);// + this->GetGeometry()[jj].FastGetSolutionStepValue(ADHESION_FORCE_Z));
	    this->GetGeometry()[jj].FastGetSolutionStepValue(FORCE) = -coef_j*beta*fsign2*(gamma*m2*nlen2);// + this->GetGeometry()[jj].FastGetSolutionStepValue(ADHESION_FORCE_X));

	    rDampingMatrix(4*jj,4*jj)     += 0.5*gamma*dt*msN[jj]*msN[jj]*nlen2 - msWorkMatrix(jj,jj);
	    rDampingMatrix(4*jj+1,4*jj+1) += 0.5*gamma*dt*msN[jj]*msN[jj]*nlen2 - msWorkMatrix(jj,jj);
	    rDampingMatrix(4*jj+2,4*jj+2) += 0.5*gamma*dt*msN[jj]*msN[jj]*nlen2 - msWorkMatrix(jj,jj);
	    rDampingMatrix(4*jj,4*ii)     += 0.5*gamma*dt*msN[jj]*msN[ii]*0.5*(nlen1 + nlen2) - msWorkMatrix(jj,ii);
	    rDampingMatrix(4*jj+1,4*ii+1) += 0.5*gamma*dt*msN[jj]*msN[ii]*0.5*(nlen1 + nlen2) - msWorkMatrix(jj,ii);
	    rDampingMatrix(4*jj+2,4*ii+2) += 0.5*gamma*dt*msN[jj]*msN[ii]*0.5*(nlen1 + nlen2) - msWorkMatrix(jj,ii);
	  }
	  else
	  {
	    rRightHandSideVector[4*jj]   -= coef_j*beta*gamma*curv2*An2[0]*area_tr;
	    rRightHandSideVector[4*jj+1] -= coef_j*beta*gamma*curv2*An2[1]*area_tr;
	    rRightHandSideVector[4*jj+2] -= coef_j*beta*gamma*curv2*An2[2]*area_tr;
	    this->GetGeometry()[jj].FastGetSolutionStepValue(FORCE_X) = -coef_j*beta*gamma*curv2*An2[0]*area_tr;
	    this->GetGeometry()[jj].FastGetSolutionStepValue(FORCE_Y) = -coef_j*beta*gamma*curv2*An2[1]*area_tr;
	    this->GetGeometry()[jj].FastGetSolutionStepValue(FORCE_Z) = -coef_j*beta*gamma*curv2*An2[2]*area_tr;
	  }

// 	  beta = 1.0;
	  if(this->GetGeometry()[kk].FastGetSolutionStepValue(TRIPLE_POINT) != 0.0)
	  {
	    coef_k = 0.5;
	    rRightHandSideVector[4*kk]   -= coef_k*beta*fsign3*(gamma*m3[0]*nlen3);// + this->GetGeometry()[kk].FastGetSolutionStepValue(ADHESION_FORCE_X));
	    rRightHandSideVector[4*kk+1] -= coef_k*beta*fsign3*(gamma*m3[1]*nlen3);// + this->GetGeometry()[kk].FastGetSolutionStepValue(ADHESION_FORCE_Y));
	    rRightHandSideVector[4*kk+2] -= coef_k*beta*fsign3*(gamma*m3[2]*nlen3);// + this->GetGeometry()[kk].FastGetSolutionStepValue(ADHESION_FORCE_Z));
	    this->GetGeometry()[kk].FastGetSolutionStepValue(FORCE_X) = -coef_k*beta*fsign3*(gamma*m3[0]*nlen3);// + this->GetGeometry()[kk].FastGetSolutionStepValue(ADHESION_FORCE_X));
	    this->GetGeometry()[kk].FastGetSolutionStepValue(FORCE_Y) = -coef_k*beta*fsign3*(gamma*m3[1]*nlen3);// + this->GetGeometry()[kk].FastGetSolutionStepValue(ADHESION_FORCE_Y));
	    this->GetGeometry()[kk].FastGetSolutionStepValue(FORCE_Z) = -coef_k*beta*fsign3*(gamma*m3[2]*nlen3);// + this->GetGeometry()[kk].FastGetSolutionStepValue(ADHESION_FORCE_Z));
	  }
	  else
	  {
	    rRightHandSideVector[4*kk]   -= coef_k*beta*gamma*curv3*An3[0]*area_tr;
	    rRightHandSideVector[4*kk+1] -= coef_k*beta*gamma*curv3*An3[1]*area_tr;
	    rRightHandSideVector[4*kk+2] -= coef_k*beta*gamma*curv3*An3[2]*area_tr;
	    this->GetGeometry()[kk].FastGetSolutionStepValue(FORCE_X) = -coef_k*beta*gamma*curv3*An3[0]*area_tr;
	    this->GetGeometry()[kk].FastGetSolutionStepValue(FORCE_Y) = -coef_k*beta*gamma*curv3*An3[1]*area_tr;
	    this->GetGeometry()[kk].FastGetSolutionStepValue(FORCE_Z) = -coef_k*beta*gamma*curv3*An3[2]*area_tr;
	  }
	}

	//Implicit treatment of surface tension
	rDampingMatrix(4*ii,4*ii)     += msWorkMatrix(ii,ii);
	rDampingMatrix(4*ii+1,4*ii+1) += msWorkMatrix(ii,ii);
	rDampingMatrix(4*ii+2,4*ii+2) += msWorkMatrix(ii,ii);

	rDampingMatrix(4*ii,4*jj)     += msWorkMatrix(ii,jj);
	rDampingMatrix(4*ii+1,4*jj+1) += msWorkMatrix(ii,jj);
	rDampingMatrix(4*ii+2,4*jj+2) += msWorkMatrix(ii,jj);

	rDampingMatrix(4*ii,4*kk)     += msWorkMatrix(ii,kk);
	rDampingMatrix(4*ii+1,4*kk+1) += msWorkMatrix(ii,kk);
	rDampingMatrix(4*ii+2,4*kk+2) += msWorkMatrix(ii,kk);

	rDampingMatrix(4*jj,4*ii)     += msWorkMatrix(jj,ii);
	rDampingMatrix(4*jj+1,4*ii+1) += msWorkMatrix(jj,ii);
	rDampingMatrix(4*jj+2,4*ii+2) += msWorkMatrix(jj,ii);

	rDampingMatrix(4*kk,4*ii)     += msWorkMatrix(kk,ii);
	rDampingMatrix(4*kk+1,4*ii+1) += msWorkMatrix(kk,ii);
	rDampingMatrix(4*kk+2,4*ii+2) += msWorkMatrix(kk,ii);

	rDampingMatrix(4*jj,4*jj)     += msWorkMatrix(jj,jj);
	rDampingMatrix(4*jj+1,4*jj+1) += msWorkMatrix(jj,jj);
	rDampingMatrix(4*jj+2,4*jj+2) += msWorkMatrix(jj,jj);

	rDampingMatrix(4*jj,4*kk)     += msWorkMatrix(jj,kk);
	rDampingMatrix(4*jj+1,4*kk+1) += msWorkMatrix(jj,kk);
	rDampingMatrix(4*jj+2,4*kk+2) += msWorkMatrix(jj,kk);

	rDampingMatrix(4*kk,4*jj)     += msWorkMatrix(kk,jj);
	rDampingMatrix(4*kk+1,4*jj+1) += msWorkMatrix(kk,jj);
	rDampingMatrix(4*kk+2,4*jj+2) += msWorkMatrix(kk,jj);

	rDampingMatrix(4*kk,4*kk)     += msWorkMatrix(kk,kk);
	rDampingMatrix(4*kk+1,4*kk+1) += msWorkMatrix(kk,kk);
	rDampingMatrix(4*kk+2,4*kk+2) += msWorkMatrix(kk,kk);


	if(k > 3 && ll != 8)
	{
	    array_1d<double,3> An4 = this->GetGeometry()[ll].FastGetSolutionStepValue(NORMAL_GEOMETRIC);
	    array_1d<double,3> m4 = this->GetGeometry()[ll].FastGetSolutionStepValue(NORMAL_EQUILIBRIUM) - this->GetGeometry()[ll].FastGetSolutionStepValue(NORMAL_TRIPLE_POINT);

	    double fsign4 = this->GetGeometry()[ll].FastGetSolutionStepValue(CONTACT_ANGLE) - theta_s;
	    fsign4 = fsign4/abs(fsign4);

	    int num_neighs_l = 0;
	    WeakPointerVector< Node<3> >& neighb_l = this->GetGeometry()[ll].GetValue(NEIGHBOUR_NODES);
	    for (unsigned int i = 0; i < neighb_l.size(); i++)
	    {
	      if (neighb_l[i].FastGetSolutionStepValue(IS_BOUNDARY) != 0.0)
		num_neighs_l++;
	    }
	    double coef_l = 0.333333333333; // 0.333333333333 | (1/num_neighs_l) | 1.0/(num_neighs_l-1)

	    double curv4 = this->GetGeometry()[ll].FastGetSolutionStepValue(MEAN_CURVATURE_3D);
	    double nlen4 = this->GetGeometry()[ll].FastGetSolutionStepValue(NODAL_LENGTH);

// 	    beta = 1.0;
	    if(this->GetGeometry()[ll].FastGetSolutionStepValue(TRIPLE_POINT) != 0.0)
	    {
	      coef_l = 0.5;
	      rRightHandSideVector[4*ll]   -= coef_l*beta*fsign4*(gamma*m4[0]*nlen4);// + this->GetGeometry()[ll].FastGetSolutionStepValue(ADHESION_FORCE_X));
	      rRightHandSideVector[4*ll+1] -= coef_l*beta*fsign4*(gamma*m4[1]*nlen4);// + this->GetGeometry()[ll].FastGetSolutionStepValue(ADHESION_FORCE_Y));
	      rRightHandSideVector[4*ll+2] -= coef_l*beta*fsign4*(gamma*m4[2]*nlen4);// + this->GetGeometry()[ll].FastGetSolutionStepValue(ADHESION_FORCE_Z));
	      this->GetGeometry()[ll].FastGetSolutionStepValue(FORCE_X) = -coef_l*beta*fsign4*(gamma*m4[0]*nlen4);// + this->GetGeometry()[ll].FastGetSolutionStepValue(ADHESION_FORCE_X));
	      this->GetGeometry()[ll].FastGetSolutionStepValue(FORCE_Y) = -coef_l*beta*fsign4*(gamma*m4[1]*nlen4);// + this->GetGeometry()[ll].FastGetSolutionStepValue(ADHESION_FORCE_Y));
	      this->GetGeometry()[ll].FastGetSolutionStepValue(FORCE_Z) = -coef_l*beta*fsign4*(gamma*m4[2]*nlen4);// + this->GetGeometry()[ll].FastGetSolutionStepValue(ADHESION_FORCE_Z));
	    }
	    else
	    {
	      rRightHandSideVector[4*ll]   -= coef_l*beta*gamma*curv4*An4[0]*area_tr;
	      rRightHandSideVector[4*ll+1] -= coef_l*beta*gamma*curv4*An4[1]*area_tr;
	      rRightHandSideVector[4*ll+2] -= coef_l*beta*gamma*curv4*An4[2]*area_tr;
	      this->GetGeometry()[ll].FastGetSolutionStepValue(FORCE) = -coef_l*beta*gamma*curv4*An4*area_tr;
	    }
        }

    }

    /// Add the viscous stress to air domain in two dimensions
    // AddViscousStress2D();
    void AddViscousStress2D()
    {
	BoundedMatrix<double,3,2> DN_DX;
	array_1d<double,3> N;
	double Area;
	GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

	double x1 = this->GetGeometry()[0].X();
	double y1 = this->GetGeometry()[0].Y();
	double x2 = this->GetGeometry()[1].X();
	double y2 = this->GetGeometry()[1].Y();
	double x3 = this->GetGeometry()[2].X();
	double y3 = this->GetGeometry()[2].Y();

	double u1 = this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_X);
	double v1 = this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_Y);
	double u2 = this->GetGeometry()[1].FastGetSolutionStepValue(VELOCITY_X);
	double v2 = this->GetGeometry()[1].FastGetSolutionStepValue(VELOCITY_Y);
	double u3 = this->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY_X);
	double v3 = this->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY_Y);

	double x13 = x1 - x3;
	double x23 = x2 - x3;
	double y13 = y1 - y3;
	double y23 = y2 - y3;

	double du_dx = (0.5/Area)*(y23*(u1 - u3) - y13*(u2 - u3));
	double du_dy = (0.5/Area)*(-x23*(u1 - u3) + x13*(u2 - u3));
	double dv_dx = (0.5/Area)*(y23*(v1 - v3) - y13*(v2 - v3));
	double dv_dy = (0.5/Area)*(-x23*(v1 - v3) + x13*(y2 - y3));

	double mu,rho;
	for(unsigned int i = 0; i < TNumNodes; ++i)
	{
	  mu = this->GetGeometry()[i].FastGetSolutionStepValue(VISCOSITY);
	  rho = this->GetGeometry()[i].FastGetSolutionStepValue(DENSITY);
	  mu *= rho;
	  this->GetGeometry()[i].FastGetSolutionStepValue(VISCOUS_STRESSX_X) += mu * ( 2*du_dx );
	  this->GetGeometry()[i].FastGetSolutionStepValue(VISCOUS_STRESSY_X) += mu * ( du_dy + dv_dx );
	  this->GetGeometry()[i].FastGetSolutionStepValue(VISCOUS_STRESSX_Y) += mu * ( dv_dx + du_dy );
	  this->GetGeometry()[i].FastGetSolutionStepValue(VISCOUS_STRESSY_Y) += mu * ( 2*dv_dy );
	}

    }

    /// Assemble the contribution from an integration point to the element's residual.
    /** Note that the dynamic term is not included in the momentum equation.
     *  If OSS_SWITCH = 1, we don't take into account the 'dynamic' stabilization
     *  terms, as it they belong to the finite element space.
     */
    void AddProjectionResidualContribution(const array_1d< double, 3 > & rAdvVel,
                                           const double Density,
                                           array_1d< double, 3 > & rElementalMomRes,
                                           double& rElementalMassRes,
                                           const array_1d< double, TNumNodes >& rShapeFunc,
                                           const BoundedMatrix<double, TNumNodes, TDim >& rShapeDeriv,
                                           const double Weight)
    {
        // If we want to use more than one Gauss point to integrate the convective term, this has to be evaluated once per integration point
        array_1d<double, TNumNodes> AGradN;
        this->GetConvectionOperator(AGradN, rAdvVel, rShapeDeriv); // Get a * grad(Ni)

        // Compute contribution to Kij * Uj, with Kij = Ni * Residual(Nj); Uj = (v,p)Node_j (column vector)
        for (unsigned int i = 0; i < TNumNodes; ++i) // Iterate over element nodes
        {

            // Variable references
            const array_1d< double, 3 > & rVelocity = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
            const array_1d< double, 3 > & rBodyForce = this->GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE);
            const double& rPressure = this->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE);

            // Compute this node's contribution to the residual (evaluated at inegration point)
            for (unsigned int d = 0; d < TDim; ++d)
            {
                rElementalMomRes[d] += Weight * (Density * (rShapeFunc[i] * rBodyForce[d] - AGradN[i] * rVelocity[d]) - rShapeDeriv(i, d) * rPressure);
                rElementalMassRes -= Weight * rShapeDeriv(i, d) * rVelocity[d];
            }
        }
    }

    /// Assemble the contribution from an integration point to the element's residual.
    /**
     * ASGS version. Note that rElementalMomRes should be initialized before calling this.
     * @param rAdvVel Convection velocity (not including subscale)
     * @param Density Fluid density evaluated at integration point
     * @param rElementalMomRes Result
     * @param rShapeFunc Shape functions evaluated at integration point
     * @param rShapeDeriv Shape function derivatives evaluated at integration point
     * @param Weight Integration point weight (as a fraction of area or volume)
     */
    void ASGSMomResidual(const array_1d< double, 3 > & rAdvVel,
                         const double Density,
                         array_1d< double, 3 > & rElementalMomRes,
                         const array_1d< double, TNumNodes >& rShapeFunc,
                         const BoundedMatrix<double, TNumNodes, TDim >& rShapeDeriv,
                         const double Weight)
    {
        // If we want to use more than one Gauss point to integrate the convective term, this has to be evaluated once per integration point
        array_1d<double, TNumNodes> AGradN;
        this->GetConvectionOperator(AGradN, rAdvVel, rShapeDeriv); // Get a * grad(Ni)

        // Compute contribution to Kij * Uj, with Kij = Ni * Residual(Nj); Uj = (v,p)Node_j (column vector)
        for (unsigned int i = 0; i < TNumNodes; ++i) // Iterate over element nodes
        {

            // Variable references
            const array_1d< double, 3 > & rVelocity = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
            const array_1d< double, 3 > & rAcceleration = this->GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION);
            const array_1d< double, 3 > & rBodyForce = this->GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE);
            const double& rPressure = this->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE);

            // Compute this node's contribution to the residual (evaluated at inegration point)
            for (unsigned int d = 0; d < TDim; ++d)
            {
                rElementalMomRes[d] += Weight * (Density * (rShapeFunc[i] * (rBodyForce[d] - rAcceleration[d]) - AGradN[i] * rVelocity[d]) - rShapeDeriv(i, d) * rPressure);
            }
        }
    }

    /// Assemble the contribution from an integration point to the element's residual.
    /**
     * OSS version. Note that rElementalMomRes should be initialized before calling this.
     * @param rAdvVel Convection velocity (not including subscale)
     * @param Density Fluid density evaluated at integration point
     * @param rElementalMomRes Result
     * @param rShapeFunc Shape functions evaluated at integration point
     * @param rShapeDeriv Shape function derivatives evaluated at integration point
     * @param Weight Integration point weight (as a fraction of area or volume)
     */
    void OSSMomResidual(const array_1d< double, 3 > & rAdvVel,
                        const double Density,
                        array_1d< double, 3 > & rElementalMomRes,
                        const array_1d< double, TNumNodes >& rShapeFunc,
                        const BoundedMatrix<double, TNumNodes, TDim >& rShapeDeriv,
                        const double Weight)
    {
        // If we want to use more than one Gauss point to integrate the convective term, this has to be evaluated once per integration point
        array_1d<double, TNumNodes> AGradN;
        this->GetConvectionOperator(AGradN, rAdvVel, rShapeDeriv); // Get a * grad(Ni)

        // Compute contribution to Kij * Uj, with Kij = Ni * Residual(Nj); Uj = (v,p)Node_j (column vector)
        for (unsigned int i = 0; i < TNumNodes; ++i) // Iterate over element nodes
        {

            // Variable references
            const array_1d< double, 3 > & rVelocity = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
            const array_1d< double, 3 > & rBodyForce = this->GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE);
            const array_1d< double, 3 > & rProjection = this->GetGeometry()[i].FastGetSolutionStepValue(ADVPROJ);
            const double& rPressure = this->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE);

            // Compute this node's contribution to the residual (evaluated at inegration point)
            for (unsigned int d = 0; d < TDim; ++d)
            {
                rElementalMomRes[d] += Weight * (Density * (rShapeFunc[i] * rBodyForce[d] - AGradN[i] * rVelocity[d]) - rShapeDeriv(i, d) * rPressure);
                rElementalMomRes[d] -= Weight * rShapeFunc[i] * rProjection[d];
            }
        }
    }


    /**
     * @brief EffectiveViscosity Calculate the viscosity at given integration point, using Smagorinsky if enabled.
     *
     * The Smagorinsky model is used only if the C_SMAGORINSKY is defined on the elemental data container.
     *
     * @note: This function is redefined when using Non-Newtonian constitutive models. It is important to keep its
     * signature, otherwise non-Newtonian models will stop working.
     *
     * @param Density The fluid's density at the integration point.
     * @param rN Nodal shape functions evaluated at the integration points (area coordinates for the point).
     * @param rDN_DX Shape function derivatives at the integration point.
     * @param ElemSize Representative length of the element (used only for Smagorinsky).
     * @param rProcessInfo ProcessInfo instance passed from the ModelPart.
     * @return Effective viscosity, in dynamic units (Pa*s or equivalent).
     */
    virtual double EffectiveViscosity(double Density,
                                      const array_1d< double, TNumNodes > &rN,
                                      const BoundedMatrix<double, TNumNodes, TDim > &rDN_DX,
                                      double ElemSize,
                                      const ProcessInfo &rProcessInfo)
    {
        const double Csmag = (static_cast< const SurfaceTension<TDim> * >(this) )->GetValue(C_SMAGORINSKY);
        double Viscosity = 0.0;
        this->EvaluateInPoint(Viscosity,VISCOSITY,rN);

        if (Csmag > 0.0)
        {
            double StrainRate = this->EquivalentStrainRate(rDN_DX); // (2SijSij)^0.5
            double LengthScale = Csmag*ElemSize;
            LengthScale *= LengthScale; // square
            Viscosity += 2.0*LengthScale*StrainRate;
        }

        return Density*Viscosity;
    }



    /**
     * @brief EquivalentStrainRate Calculate the second invariant of the strain rate tensor GammaDot = (2SijSij)^0.5.
     *
     * @note Our implementation of non-Newtonian consitutive models such as Bingham relies on this funcition being
     * defined on all fluid elements.
     *
     * @param rDN_DX Shape function derivatives at the integration point.
     * @return GammaDot = (2SijSij)^0.5.
     */
    double EquivalentStrainRate(const BoundedMatrix<double, TNumNodes, TDim > &rDN_DX) const;


    /// Write the advective velocity evaluated at this point to an array
    /**
     * Writes the value of the advective velocity evaluated at a point inside
     * the element to an array_1d
     * @param rAdvVel: Output array
     * @param rShapeFunc: Shape functions evaluated at the point of interest
     */
    void GetAdvectiveVel(array_1d< double, 3 > & rAdvVel,
                                 const array_1d< double, TNumNodes >& rShapeFunc)
    {
        // Compute the weighted value of the advective velocity in the (Gauss) Point
        GeometryType& rGeom = this->GetGeometry();
        rAdvVel = rShapeFunc[0] * (rGeom[0].FastGetSolutionStepValue(VELOCITY) - rGeom[0].FastGetSolutionStepValue(MESH_VELOCITY));
        for (unsigned int iNode = 1; iNode < TNumNodes; ++iNode)
            rAdvVel += rShapeFunc[iNode] * (rGeom[iNode].FastGetSolutionStepValue(VELOCITY) - rGeom[iNode].FastGetSolutionStepValue(MESH_VELOCITY));
    }

    /// Write the advective velocity evaluated at this point to an array
    /**
     * Writes the value of the advective velocity evaluated at a point inside
     * the element to an array_1d
     * @param rAdvVel: Output array
     * @param rShapeFunc: Shape functions evaluated at the point of interest
     * @param Step: The time Step
     */
    virtual void GetAdvectiveVel(array_1d< double, 3 > & rAdvVel,
                                 const array_1d< double, TNumNodes >& rShapeFunc,
                                 const std::size_t Step)
    {
        // Compute the weighted value of the advective velocity in the (Gauss) Point
        GeometryType& rGeom = this->GetGeometry();
        rAdvVel = rShapeFunc[0] * (rGeom[0].FastGetSolutionStepValue(VELOCITY, Step) - rGeom[0].FastGetSolutionStepValue(MESH_VELOCITY, Step));
        for (unsigned int iNode = 1; iNode < TNumNodes; ++iNode)
            rAdvVel += rShapeFunc[iNode] * (rGeom[iNode].FastGetSolutionStepValue(VELOCITY, Step) - rGeom[iNode].FastGetSolutionStepValue(MESH_VELOCITY, Step));
    }

    /// Write the convective operator evaluated at this point (for each nodal funciton) to an array
    /**
     * Evaluate the convective operator for each node's shape function at an arbitrary point
     * @param rResult: Output vector
     * @param rVelocity: Velocity evaluated at the integration point
     * @param rShapeDeriv: Derivatives of shape functions evaluated at the integration point
     * @see GetAdvectiveVel provides rVelocity
     */
    void GetConvectionOperator(array_1d< double, TNumNodes >& rResult,
                               const array_1d< double, 3 > & rVelocity,
                               const BoundedMatrix<double, TNumNodes, TDim >& rShapeDeriv)
    {
        // Evaluate (and weight) the a * Grad(Ni) operator in the integration point, for each node i
        for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode) // Loop over nodes
        {
            // Initialize result
            rResult[iNode] = rVelocity[0] * rShapeDeriv(iNode, 0);
            for (unsigned int d = 1; d < TDim; ++d) // loop over components
                rResult[iNode] += rVelocity[d] * rShapeDeriv(iNode, d);
        }
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
        // Compute the weighted value of the nodal variable in the (Gauss) Point
        GeometryType& rGeom = this->GetGeometry();
        rResult = rShapeFunc[0] * rGeom[0].FastGetSolutionStepValue(rVariable);
        for (unsigned int iNode = 1; iNode < TNumNodes; ++iNode)
            rResult += rShapeFunc[iNode] * rGeom[iNode].FastGetSolutionStepValue(rVariable);
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
        // Compute the weighted value of the nodal variable in the (Gauss) Point
        GeometryType& rGeom = this->GetGeometry();
        rResult = rShapeFunc[0] * rGeom[0].FastGetSolutionStepValue(rVariable);
        for (unsigned int iNode = 1; iNode < TNumNodes; ++iNode)
            rResult += rShapeFunc[iNode] * rGeom[iNode].FastGetSolutionStepValue(rVariable);
    }

    /// Return an estimate for the element size h, used to calculate the stabilization parameters
    /**
     * Estimate the element size from its area or volume, required to calculate stabilization parameters.
     * Note that its implementation is different for 2D or 3D elements.
     * @see ST2D, ST3D for actual implementation
     * @param Volume (in 3D) or Area (in 2D) of the element
     * @return Element size h
     */
    double ElementSize(const double);


    /// Adds the contribution of the viscous term to the momentum equation.
    /**
     * The viscous term is written in stress-divergence (Cauchy) form.
     * @param rDampingMatrix Elemental Damping matrix
     * @param rShapeDeriv Elemental shape function derivatives
     * @param Weight Effective viscosity, in dynamic units, weighted by the integration point area
     */
    virtual void AddViscousTerm(MatrixType& rDampingMatrix,
                                const BoundedMatrix<double, TNumNodes, TDim >& rShapeDeriv,
                                const double Weight);

    /// Adds the contribution of the viscous term to the momentum equation (alternate).
    /**
     * This function is an alternate implementation of ST::AddViscousTerm.
     * This version works with ublas matrices, using the relationship between stress and
     * rate of strain given by ST::CalculateC. It is currently unused (as ST::AddViscousTerm
     * is a more efficient implementation of the Cauchy equation) but it is left here so derived
     * classes can use it to implement other constitutive equations.
     * @param rDampingMatrix Elemental Damping matrix
     * @param rShapeDeriv Elemental shape function derivatives
     * @param Weight Effective viscosity, in dynamic units, weighted by the integration point area
     */
    void AddBTransCB(MatrixType& rDampingMatrix,
                 BoundedMatrix<double, TNumNodes, TDim >& rShapeDeriv,
                     const double Weight)
    {
        BoundedMatrix<double, (TDim * TNumNodes)/2, TDim*TNumNodes > B;
        BoundedMatrix<double, (TDim * TNumNodes)/2, (TDim*TNumNodes)/2 > C;
        this->CalculateB(B,rShapeDeriv);
        this->CalculateC(C,Weight);

        const unsigned int BlockSize = TDim + 1;
        const unsigned int StrainSize = (TDim*TNumNodes)/2;

        vector<unsigned int> aux(TDim*TNumNodes);
        for(unsigned int i=0; i<TNumNodes; i++)
        {
            int base_index = TDim*i;
            int aux_index = BlockSize*i;
            for(unsigned int j=0; j<TDim; j++)
            {
                aux[base_index+j] = aux_index+j;
            }
        }

        for(unsigned int k=0; k< StrainSize; k++)
        {
            for(unsigned int l=0; l< StrainSize; l++)
            {
                const double Ckl = C(k,l);
                for (unsigned int i = 0; i < TDim*TNumNodes; ++i) // iterate over v components (vx,vy[,vz])
                {
                    const double Bki=B(k,i);
                    for (unsigned int j = 0; j < TDim*TNumNodes; ++j) // iterate over u components (ux,uy[,uz])
                    {
                        rDampingMatrix(aux[i],aux[j]) += Bki*Ckl*B(l,j);
                    }

                }

            }
        }
    }

    void ModulatedGradientDiffusion(MatrixType& rDampingMatrix,
            const BoundedMatrix<double, TNumNodes, TDim >& rDN_DX,
            const double Weight)
    {
        const GeometryType& rGeom = this->GetGeometry();

        // Velocity gradient
        MatrixType GradU = ZeroMatrix(TDim,TDim);
        for (unsigned int n = 0; n < TNumNodes; n++)
        {
            const array_1d<double,3>& rVel = this->GetGeometry()[n].FastGetSolutionStepValue(VELOCITY);
            for (unsigned int i = 0; i < TDim; i++)
                for (unsigned int j = 0; j < TDim; j++)
                    GradU(i,j) += rDN_DX(n,j)*rVel[i];
        }

        // Element lengths
        array_1d<double,3> Delta(3,0.0);
        Delta[0] = std::fabs(rGeom[TNumNodes-1].X()-rGeom[0].X());
        Delta[1] = std::fabs(rGeom[TNumNodes-1].Y()-rGeom[0].Y());
        Delta[2] = std::fabs(rGeom[TNumNodes-1].Z()-rGeom[0].Z());

        for (unsigned int n = 1; n < TNumNodes; n++)
        {
            double hx = std::fabs(rGeom[n].X()-rGeom[n-1].X());
            if (hx > Delta[0]) Delta[0] = hx;
            double hy = std::fabs(rGeom[n].Y()-rGeom[n-1].Y());
            if (hy > Delta[1]) Delta[1] = hy;
            double hz = std::fabs(rGeom[n].Z()-rGeom[n-1].Z());
            if (hz > Delta[2]) Delta[2] = hz;
        }

        double AvgDeltaSq = Delta[0];
        for (unsigned int d = 1; d < TDim; d++)
            AvgDeltaSq *= Delta[d];
        AvgDeltaSq = std::pow(AvgDeltaSq,2./TDim);

        Delta[0] = Delta[0]*Delta[0]/12.0;
        Delta[1] = Delta[1]*Delta[1]/12.0;
        Delta[2] = Delta[2]*Delta[2]/12.0;

        // Gij
        MatrixType G = ZeroMatrix(TDim,TDim);
        for (unsigned int i = 0; i < TDim; i++)
            for (unsigned int j = 0; j < TDim; j++)
                for (unsigned int d = 0; d < TDim; d++)
                    G(i,j) += Delta[d]*GradU(i,d)*GradU(j,d);

        // Gij:Sij
        double GijSij = 0.0;
        for (unsigned int i = 0; i < TDim; i++)
            for (unsigned int j = 0; j < TDim; j++)
                GijSij += 0.5*G(i,j)*( GradU(i,j) + GradU(j,i) );

        if (GijSij < 0.0) // Otherwise model term is clipped
        {
            // Gkk
            double Gkk = G(0,0);
            for (unsigned int d = 1; d < TDim; d++)
                Gkk += G(d,d);

            // C_epsilon
            const double Ce = 1.0;

            // ksgs
            double ksgs = -4*AvgDeltaSq*GijSij/(Ce*Ce*Gkk);

            // Assembly of model term
            unsigned int RowIndex = 0;
            unsigned int ColIndex = 0;

            for (unsigned int i = 0; i < TNumNodes; i++)
            {
                for (unsigned int j = 0; j < TNumNodes; j++)
                {
                    for (unsigned int d = 0; d < TDim; d++)
                    {
                        double Aux = rDN_DX(i,d) * Delta[0] * G(d,0)*rDN_DX(j,0);
                        for (unsigned int k = 1; k < TDim; k++)
                            Aux += rDN_DX(i,d) *Delta[k] * G(d,k)*rDN_DX(j,k);
                        rDampingMatrix(RowIndex+d,ColIndex+d) += Weight * 2.0*ksgs *  Aux;
                    }

                    ColIndex += TDim;
                }
                RowIndex += TDim;
                ColIndex = 0;
            }
        }

    }

    /// Calculate the strain rate matrix
    /**
     * Unused, left to support derived classes. @see SurfaceTension::AddBTransCB
     * @param rB Strain rate matrix
     * @param rShapeDeriv Nodal shape funcion derivatives
     */
    void CalculateB( BoundedMatrix<double, (TDim * TNumNodes) / 2, TDim * TNumNodes >& rB,
                     const BoundedMatrix<double, TNumNodes, TDim >& rShapeDeriv);

    /// Calculate a matrix that provides the stress given the strain rate
    /**
     * Unused, left to support derived classes. @see SurfaceTension::AddBTransCB.
     * Note that only non-zero terms are written, so the output matrix should be
     * initialized before calling this.
     * @param rC Matrix representation of the stress tensor (output)
     * @param Viscosity Effective viscosity, in dynamic units, weighted by the integration point area
     */
    virtual void CalculateC( BoundedMatrix<double, (TDim * TNumNodes) / 2, (TDim * TNumNodes) / 2 >& rC,
                             const double Viscosity);

    double ConsistentMassCoef(const double Area);


    double SubscaleErrorEstimate(const ProcessInfo& rProcessInfo)
    {
        // Get the element's geometric parameters
        double Area;
        array_1d<double, TNumNodes> N;
        BoundedMatrix<double, TNumNodes, TDim> DN_DX;
        GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

        // Calculate this element's fluid properties
        double Density;
        this->EvaluateInPoint(Density, DENSITY, N);

        double ElemSize = this->ElementSize(Area);
        double Viscosity = this->EffectiveViscosity(Density,N,DN_DX,ElemSize,rProcessInfo);

        // Get Advective velocity
        array_1d<double, 3 > AdvVel;
        this->GetAdvectiveVel(AdvVel, N);

        // Output container
        array_1d< double, 3 > ElementalMomRes(3, 0.0);

        // Calculate stabilization parameter. Note that to estimate the subscale velocity, the dynamic coefficient in TauOne is assumed zero.
        double TauOne;
        this->CalculateStaticTau(TauOne,AdvVel,ElemSize,Density,Viscosity);

        if ( rProcessInfo[OSS_SWITCH] != 1 ) // ASGS
        {
            this->ASGSMomResidual(AdvVel,Density,ElementalMomRes,N,DN_DX,1.0);
            ElementalMomRes *= TauOne;
        }
        else // OSS
        {
            this->OSSMomResidual(AdvVel,Density,ElementalMomRes,N,DN_DX,1.0);;
            ElementalMomRes *= TauOne;
        }

        // Error estimation ( ||U'|| / ||Uh_gauss|| ), taking ||U'|| = TauOne ||MomRes||
        double ErrorRatio(0.0);//, UNorm(0.0);
        //array_1d< double, 3 > UGauss(3, 0.0);
        //this->EvaluateInPoint(UGauss,VELOCITY,N);

        for (unsigned int i = 0; i < TDim; ++i)
        {
            ErrorRatio += ElementalMomRes[i] * ElementalMomRes[i];
            //UNorm += UGauss[i] * UGauss[i];
        }
        ErrorRatio = sqrt(ErrorRatio*Area);// / UNorm);
        //ErrorRatio /= Density;
        //this->SetValue(ERROR_RATIO, ErrorRatio);
        return ErrorRatio;
    }



    array_1d<double,2> Vector2D(const double x0, const double y0, const double x1, const double y1)
    {
      array_1d<double,2> r01;
      r01[0] = x1 - x0;
      r01[1] = y1 - y0;
      return (r01);
    }

    array_1d<double,3> Vector3D(const double x0, const double y0, const double z0, const double x1, const double y1, const double z1)
    {
      array_1d<double,3> r01;
      r01[0] = x1 - x0;
      r01[1] = y1 - y0;
      r01[2] = z1 - z0;
      return (r01);
    }

    array_1d<double,3> CrossProduct3D(const array_1d<double,3>& a, const array_1d<double,3>& b)
    {
      array_1d<double,3> c;
      c[0] = a[1]*b[2] - a[2]*b[1];
      c[1] = a[2]*b[0] - a[0]*b[2];
      c[2] = a[0]*b[1] - a[1]*b[0];
      return (c);
    }

    double Norm2D(const array_1d<double,2>& a)
    {
      return sqrt(a[0]*a[0] + a[1]*a[1]);
    }

    double Norm3D(const array_1d<double,3>& a)
    {
      return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
    }

    void NormalizeVec2D(array_1d<double,2>& input)
    {
      double norm = Norm2D(input);
      if (norm != 0.0)
      {
	input[0] /= norm;
	input[1] /= norm;
      }
    }

    void NormalizeVec3D(array_1d<double,3>& input)
    {
      double norm = Norm3D(input);
      if (norm != 0.0)
      {
	input[0] /= norm;
	input[1] /= norm;
	input[2] /= norm;
      }
    }

    double DotProduct2D(const array_1d<double,2>& a, const array_1d<double,2>& b)
    {
      return (a[0]*b[0] + a[1]*b[1]);
    }

    double DotProduct3D(const array_1d<double,3>& a, const array_1d<double,3>& b)
    {
      return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
    }

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




     void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
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
    SurfaceTension & operator=(SurfaceTension const& rOther);

    /// Copy constructor.
    SurfaceTension(SurfaceTension const& rOther);

    ///@}

}; // Class SurfaceTension

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< unsigned int TDim,
          unsigned int TNumNodes >
inline std::istream& operator >>(std::istream& rIStream,
                                 SurfaceTension<TDim, TNumNodes>& rThis)
{
    return rIStream;
}

/// output stream function
template< unsigned int TDim,
          unsigned int TNumNodes >
inline std::ostream& operator <<(std::ostream& rOStream,
                                 const SurfaceTension<TDim, TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} // ULF Application group

} // namespace Kratos.

#endif // KRATOS_ST_H_INCLUDED  defined
