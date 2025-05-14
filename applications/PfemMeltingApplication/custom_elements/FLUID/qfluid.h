//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//


#if !defined(KRATOS_QFLUID_H_INCLUDED )
#define  KRATOS_QFLUID_H_INCLUDED

// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "containers/array_1d.h"
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "includes/cfd_variables.h"
#include "utilities/geometry_utilities.h"

// Application includes
#include "pfem_melting_application_variables.h"

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

/// A stabilized element for the incompressible Navier-Stokes equations.
/**
 * This class implements a stabilized formulation based on the
 * Variational Multiscale framework. The the subscales can be modeled
 * using either Algebraic Subgird Scales (ASGS) or Orthogonal Subscales (OSS).
 * In the case of OSS, the projection terms are treated explicitly (computed
 * using the results of the previous iteration) and the subscales are not
 * tracked in time. The choice of subscale model is made based on the ProcessInfo
 * variable OSS_SWITCH (OSS if 1, ASGS otherwise).
 * This class implements both the 2D and 3D versions of the element.
 *
 * The ASGS implementation follows Ramon Codina, A stabilized finite element
 * method for generalized stationary incompressible flows, Computer Methods in
 * Applied Mechanics and Engineering. Vol. 190 (2001), 2681-2706.
 *
 * The OSS implementation corresponds to the case identified as explicit, quasi-
 * static orthogonal subscales in Ramon Codina, Stabilized finite element approximation
 * of transient incompressible flows using orthogonal subscales, Computer Methods
 * in Applied Mechanics and Engineering. Vol. 191 (2002), 4295-4321.
 *
 * In addition to the stabilization, this element implements the Smagorinsky
 * model of turbulence. This turbulent term is only activated if the elemental
 * value C_SMAGORINSKY is set to something other than zero.
 *
 * This class requires at least the following variables:\n
 * On each Node, as solution step variables VELOCITY, PRESSURE, ACCELERATION, MESH_VELOCITY, DENSITY, VISCOSITY.\n
 * On ProcessInfo OSS_SWITCH, DYNAMIC_TAU, DELTA_TIME.\n
 * If OSS is used, the nodes also require NODAL_AREA, ADVPROJ and DIVPROJ as solution step variables.\n
 * If Smagorinsky is used, C_SMAGORINSKY has to be defined on the elements.\n
 * Error estimation stores ERROR_RATIO on the elements.\n
 * Some additional variables can be used to print results on the element: SUBSCALE_VELOCITY, SUBSCALE_PRESSURE, TAUONE, TAUTWO, MU, VORTICITY.
 *
 * @see ResidualBasedEliminationBuilderAndSolver compatible monolithic solution strategy.
 * @see PressureSplittingBuilderAndSolver compatible segregated solution strategy.
 * @see TrilinosPressureSplittingBuilderAndSolver compatible mpi strategy.
 * @see DynamicSmagorinskyUtils to set the Smagorinsky parameter dynamically.
 * @see ResidualBasedPredictorCorrectorVelocityBossakScheme time scheme that can use
 * OSS stabilization.
 */
template< unsigned int TDim,
          unsigned int TNumNodes = TDim + 1 >
class KRATOS_API(PFEM_MELTING_APPLICATION) QFLUID : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of QFLUID
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(QFLUID);

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
    QFLUID(IndexType NewId = 0) :
        Element(NewId)
    {}

    ///Constructor using an array of nodes.
    /**
     * @param NewId Index of the new element
     * @param ThisNodes An array containing the nodes of the new element
     */
    QFLUID(IndexType NewId, const NodesArrayType& ThisNodes) :
        Element(NewId, ThisNodes)
    {}

    /// Constructor using a geometry object.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     */
    QFLUID(IndexType NewId, GeometryType::Pointer pGeometry) :
        Element(NewId, pGeometry)
    {}

    /// Constuctor using geometry and properties.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     * @param pProperties Pointer to the element's properties
     */
    QFLUID(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) :
        Element(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    ~QFLUID() override
    {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /// Create a new element of this type
    /**
     * Returns a pointer to a new QFLUID element, created using given input
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
                            PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive< QFLUID<TDim, TNumNodes> >(NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    Element::Pointer Create(IndexType NewId,
                           GeometryType::Pointer pGeom,
                           PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive< QFLUID<TDim, TNumNodes> >(NewId, pGeom, pProperties);
    }

    /// Provides local contributions from body forces and OSS projection terms
    /**
     * This is called during the assembly process and provides the terms of the
     * system that are either constant or computed explicitly (from the 'old'
     * iteration variables). In this case this means the body force terms and the
     * OSS projections, that are treated explicitly.
     * @param rLeftHandSideMatrix the elemental left hand side matrix. Not used here, required for compatibility purposes only.
     * @param rRightHandSideVector the elemental right hand side
     * @param rCurrentProcessInfo the current process info
     */
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      const ProcessInfo& rCurrentProcessInfo) override
    {
        const unsigned int LocalSize = TDim  * TNumNodes;

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
                                const ProcessInfo& rCurrentProcessInfo) override
    {
        const unsigned int LocalSize = TDim * TNumNodes;

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
                                const ProcessInfo& rCurrentProcessInfo) override
     {
        const unsigned int LocalSize = TDim * TNumNodes;
	//KRATOS_WATCH("ALLLLLLLLLLA")
        // Check sizes and initialize
        if (rRightHandSideVector.size() != LocalSize)
            rRightHandSideVector.resize(LocalSize, false);

        noalias(rRightHandSideVector) = ZeroVector(LocalSize);

	//KRATOS_THROW_ERROR(std::logic_error,"pressure calculation 3D not implemented","");

        /*if(rRightHandSideVector.size() != 6)
        	rRightHandSideVector.resize(6,false);*/

        noalias(rRightHandSideVector) = ZeroVector(LocalSize);
    	//noalias(rRightHandSideVector) = ZeroVector(6);
    	//KRATOS_THROW_ERROR(std::logic_error,"not dereeeeeeeeeeeeeeeeeeeeeeecha","");
/*    KRATOS_WATCH("SOLIDOOOOOOOOOOOOOOOOOOOOOOOOOO")
    KRATOS_WATCH("SOLIDOOOOOOOOOOOOOOOOOOOOOOOOOO")
    KRATOS_WATCH("SOLIDOOOOOOOOOOOOOOOOOOOOOOOOOO")*/

        // Calculate this element's geometric parameters
        double Area;
        array_1d<double, TNumNodes> N;
        BoundedMatrix<double, TNumNodes, TDim> DN_DX;
        GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);


        // Calculate this element's fluid properties
        double Density;
        this->EvaluateInPoint(Density, DENSITY, N);

        //KRATOS_WATCH("DENSITY")
        //KRATOS_WATCH(Density)

        // Calculate Momentum RHS contribution
        //this->AddMomentumRHS(rRightHandSideVector, Density, N, Area);


        //writing the body force
	    const array_1d<double,3>& body_force = 0.333333333*(this->GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE)+ this->GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE) + this->GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE));
	    //const array_1d<double,3>& body_force = GetProperties()[BODY_FORCE];

	    if(TDim==2){
	    for(unsigned int i = 0; i<TNumNodes; i++)
	    {
		rRightHandSideVector[i*2] = body_force[0]* Density * 0.3333333333333;
		rRightHandSideVector[i*2+1] = body_force[1] * Density * 0.3333333333333;

	    }}
	    else{
	    for(unsigned int i = 0; i<TNumNodes; i++)
	    {
	    const array_1d< double, 3 > & rBodyForce = this->GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE);
	    for (unsigned int d = 0; d < TDim; ++d)
		    {
		        rRightHandSideVector[d] += 0.25 * Density * rBodyForce[d];
		    }
            }

	    }

	    if(TDim==2)
	    {
	    //get the value of Cauchy stress at the Gauss point. It is given by:
		    const BoundedMatrix<double,2,2> & CauchyStress=this->GetValue(CAUCHY_STRESS_TENSOR);

		    //KRATOS_WATCH(CauchyStress)
		    /*
		    //dN1/dx*SigmaX + 0 + dN1/dy*TauXY
		    rRightHandSideVector[0] -= DN_DX(0,0)*CauchyStress(0,0) + DN_DX(0,1)*CauchyStress(0,1) ;
		    //0  +   dN1/dy*SigmaY + 0 + dN1/dx*TauXY
		    rRightHandSideVector[1] -= DN_DX(0,1)*CauchyStress(1,1) + DN_DX(0,0)*CauchyStress(0,1) ;
		    //dN2/dx*SigmaX + 0 + dN2/dy*TauXY
		    rRightHandSideVector[2] -= DN_DX(1,0)*CauchyStress(0,0) + DN_DX(1,1)*CauchyStress(0,1) ;
		    //0  +   dN2/dy*SigmaY + 0 + dN2/dx*TauXY
		    rRightHandSideVector[3] -= DN_DX(1,1)*CauchyStress(1,1) + DN_DX(1,0)*CauchyStress(0,1) ;
		    //dN3/dx*SigmaX + 0 + dN3/dy*TauXY
		    rRightHandSideVector[4] -= DN_DX(2,0)*CauchyStress(0,0) + DN_DX(2,1)*CauchyStress(0,1) ;
		    //0  +   dN3/dy*SigmaY + 0 + dN3/dx*TauXY
		    rRightHandSideVector[5] -= DN_DX(2,1)*CauchyStress(1,1) + DN_DX(2,0)*CauchyStress(0,1) ;

		    rRightHandSideVector*=Area;

		    */

		    double SXX= 0.333333333 * this->GetGeometry()[0].FastGetSolutionStepValue(HISTORICAL_SIGMA_XX) + 0.333333333 * this->GetGeometry()[1].FastGetSolutionStepValue(HISTORICAL_SIGMA_XX)+0.333333333 * this->GetGeometry()[2].FastGetSolutionStepValue(HISTORICAL_SIGMA_XX);
		    double SXY= 0.333333333 * this->GetGeometry()[0].FastGetSolutionStepValue(HISTORICAL_SIGMA_XY) + 0.333333333 * this->GetGeometry()[1].FastGetSolutionStepValue(HISTORICAL_SIGMA_XY)+0.333333333 * this->GetGeometry()[2].FastGetSolutionStepValue(HISTORICAL_SIGMA_XY);
		    //double SXZ= 0.333333333 * this->GetGeometry()[0].FastGetSolutionStepValue(HISTORICAL_SIGMA_XZ) + 0.333333333 * this->GetGeometry()[1].FastGetSolutionStepValue(HISTORICAL_SIGMA_XZ)+0.333333333 * this->GetGeometry()[2].FastGetSolutionStepValue(HISTORICAL_SIGMA_XZ);

		    double SYX= 0.333333333 * this->GetGeometry()[0].FastGetSolutionStepValue(HISTORICAL_SIGMA_YX) + 0.333333333 * this->GetGeometry()[1].FastGetSolutionStepValue(HISTORICAL_SIGMA_YX)+0.333333333 * this->GetGeometry()[2].FastGetSolutionStepValue(HISTORICAL_SIGMA_YX);
		    double SYY= 0.333333333 * this->GetGeometry()[0].FastGetSolutionStepValue(HISTORICAL_SIGMA_YY) + 0.333333333 * this->GetGeometry()[1].FastGetSolutionStepValue(HISTORICAL_SIGMA_YY)+0.333333333 * this->GetGeometry()[2].FastGetSolutionStepValue(HISTORICAL_SIGMA_YY);
		    //double SYZ= 0.333333333 * this->GetGeometry()[0].FastGetSolutionStepValue(HISTORICAL_SIGMA_YZ) + 0.333333333 * this->GetGeometry()[1].FastGetSolutionStepValue(HISTORICAL_SIGMA_YZ)+0.333333333 * this->GetGeometry()[2].FastGetSolutionStepValue(HISTORICAL_SIGMA_YZ);

		    //double SZX= 0.333333333 * this->GetGeometry()[0].FastGetSolutionStepValue(HISTORICAL_SIGMA_ZX) + 0.333333333 * this->GetGeometry()[1].FastGetSolutionStepValue(HISTORICAL_SIGMA_ZX)+0.333333333 * this->GetGeometry()[2].FastGetSolutionStepValue(HISTORICAL_SIGMA_ZX);
		    //double SZY= 0.333333333 * this->GetGeometry()[0].FastGetSolutionStepValue(HISTORICAL_SIGMA_ZY) + 0.333333333 * this->GetGeometry()[1].FastGetSolutionStepValue(HISTORICAL_SIGMA_ZY)+0.333333333 * this->GetGeometry()[2].FastGetSolutionStepValue(HISTORICAL_SIGMA_ZY);
		    //double SZZ= 0.333333333 * this->GetGeometry()[0].FastGetSolutionStepValue(HISTORICAL_SIGMA_ZZ) + 0.333333333 * this->GetGeometry()[1].FastGetSolutionStepValue(HISTORICAL_SIGMA_ZZ)+0.333333333 * this->GetGeometry()[2].FastGetSolutionStepValue(HISTORICAL_SIGMA_ZZ);


		    rRightHandSideVector[0] -= DN_DX(0,0)*SXX + DN_DX(0,1)*SXY;

		    //0  +   dN1/dy*SigmaY + 0 + dN1/dx*TauXY
		    rRightHandSideVector[1] -= DN_DX(0,0)*SYX + DN_DX(0,1)*SYY;

		    //dN2/dx*SigmaX + 0 + dN2/dy*TauXY
		    rRightHandSideVector[2] -= DN_DX(1,0)*SXX + DN_DX(1,1)*SXY;

		    //0  +   dN2/dy*SigmaY + 0 + dN2/dx*TauXY
		    rRightHandSideVector[3] -= DN_DX(1,0)*SYX + DN_DX(1,1)*SYY;

		    //dN3/dx*SigmaX + 0 + dN3/dy*TauXY
		    rRightHandSideVector[4] -= DN_DX(2,0)*SXX + DN_DX(2,1)*SXY;

		    //0  +   dN3/dy*SigmaY + 0 + dN3/dx*TauXY
		    rRightHandSideVector[5] -= DN_DX(2,0)*SYX + DN_DX(2,1)*SYY;

		    rRightHandSideVector*=Area;




    	    }
	else{
		    const BoundedMatrix<double,3,3> & CauchyStress=this->GetValue(CAUCHY_STRESS_TENSOR);
	            //KRATOS_WATCH(CauchyStress)
	            /*
		    //dN1/dx*SigmaX + 0 + dN1/dy*TauXY
		    rRightHandSideVector[0] -= DN_DX(0,0)*CauchyStress(0,0) + DN_DX(0,1)*CauchyStress(0,1) + DN_DX(0,2)*CauchyStress(0,2);
		    //0  +   dN1/dy*SigmaY + 0 + dN1/dx*TauXY
		    rRightHandSideVector[1] -= DN_DX(0,0)*CauchyStress(1,0) + DN_DX(0,1)*CauchyStress(1,1) + DN_DX(0,2)*CauchyStress(1,2);

		    rRightHandSideVector[2] -= DN_DX(0,0)*CauchyStress(2,0) + DN_DX(0,1)*CauchyStress(2,1) + DN_DX(0,2)*CauchyStress(2,2);


		    //dN2/dx*SigmaX + 0 + dN2/dy*TauXY
		    rRightHandSideVector[3] -= DN_DX(1,0)*CauchyStress(0,0) + DN_DX(1,1)*CauchyStress(0,1) + DN_DX(1,2)*CauchyStress(0,2);
		    //0  +   dN2/dy*SigmaY + 0 + dN2/dx*TauXY
		    rRightHandSideVector[4] -= DN_DX(1,0)*CauchyStress(1,0) + DN_DX(1,1)*CauchyStress(1,1) + DN_DX(1,2)*CauchyStress(1,2) ;

		    rRightHandSideVector[5] -= DN_DX(1,0)*CauchyStress(2,0) + DN_DX(1,1)*CauchyStress(2,1) + DN_DX(1,2)*CauchyStress(2,2) ;



		    //dN3/dx*SigmaX + 0 + dN3/dy*TauXY
		    rRightHandSideVector[6] -= DN_DX(2,0)*CauchyStress(0,0) + DN_DX(2,1)*CauchyStress(0,1) + DN_DX(2,2)*CauchyStress(0,2) ;
		    //0  +   dN3/dy*SigmaY + 0 + dN3/dx*TauXY
		    rRightHandSideVector[7] -= DN_DX(2,0)*CauchyStress(1,0) + DN_DX(2,1)*CauchyStress(1,1) + DN_DX(2,2)*CauchyStress(1,2) ;

		    rRightHandSideVector[8] -= DN_DX(2,0)*CauchyStress(2,0) + DN_DX(2,1)*CauchyStress(2,1) + DN_DX(2,2)*CauchyStress(2,2) ;

		     //dN3/dx*SigmaX + 0 + dN3/dy*TauXY
		    rRightHandSideVector[9] -= DN_DX(3,0)*CauchyStress(0,0) + DN_DX(3,1)*CauchyStress(0,1) + DN_DX(3,2)*CauchyStress(0,2) ;
		    //0  +   dN3/dy*SigmaY + 0 + dN3/dx*TauXY
		    rRightHandSideVector[10] -= DN_DX(3,0)*CauchyStress(1,0) + DN_DX(3,1)*CauchyStress(1,1) + DN_DX(3,2)*CauchyStress(1,2) ;

		    rRightHandSideVector[11] -= DN_DX(3,0)*CauchyStress(2,0) + DN_DX(3,1)*CauchyStress(2,1) + DN_DX(3,2)*CauchyStress(2,2) ;

		    */


		    //dN1/dx*SigmaX + 0 + dN1/dy*TauXY

		    double SXX= 0.25 * this->GetGeometry()[0].FastGetSolutionStepValue(HISTORICAL_SIGMA_XX) + 0.25 * this->GetGeometry()[1].FastGetSolutionStepValue(HISTORICAL_SIGMA_XX)+0.25 * this->GetGeometry()[2].FastGetSolutionStepValue(HISTORICAL_SIGMA_XX)+0.25 * this->GetGeometry()[3].FastGetSolutionStepValue(HISTORICAL_SIGMA_XX);
		    double SXY= 0.25 * this->GetGeometry()[0].FastGetSolutionStepValue(HISTORICAL_SIGMA_XY) + 0.25 * this->GetGeometry()[1].FastGetSolutionStepValue(HISTORICAL_SIGMA_XY)+0.25 * this->GetGeometry()[2].FastGetSolutionStepValue(HISTORICAL_SIGMA_XY)+0.25 * this->GetGeometry()[3].FastGetSolutionStepValue(HISTORICAL_SIGMA_XY);
		    double SXZ= 0.25 * this->GetGeometry()[0].FastGetSolutionStepValue(HISTORICAL_SIGMA_XZ) + 0.25 * this->GetGeometry()[1].FastGetSolutionStepValue(HISTORICAL_SIGMA_XZ)+0.25 * this->GetGeometry()[2].FastGetSolutionStepValue(HISTORICAL_SIGMA_XZ)+0.25 * this->GetGeometry()[3].FastGetSolutionStepValue(HISTORICAL_SIGMA_XZ);

		    double SYX= 0.25 * this->GetGeometry()[0].FastGetSolutionStepValue(HISTORICAL_SIGMA_YX) + 0.25 * this->GetGeometry()[1].FastGetSolutionStepValue(HISTORICAL_SIGMA_YX)+0.25 * this->GetGeometry()[2].FastGetSolutionStepValue(HISTORICAL_SIGMA_YX)+0.25 * this->GetGeometry()[3].FastGetSolutionStepValue(HISTORICAL_SIGMA_YX);
		    double SYY= 0.25 * this->GetGeometry()[0].FastGetSolutionStepValue(HISTORICAL_SIGMA_YY) + 0.25 * this->GetGeometry()[1].FastGetSolutionStepValue(HISTORICAL_SIGMA_YY)+0.25 * this->GetGeometry()[2].FastGetSolutionStepValue(HISTORICAL_SIGMA_YY)+0.25 * this->GetGeometry()[3].FastGetSolutionStepValue(HISTORICAL_SIGMA_YY);
		    double SYZ= 0.25 * this->GetGeometry()[0].FastGetSolutionStepValue(HISTORICAL_SIGMA_YZ) + 0.25 * this->GetGeometry()[1].FastGetSolutionStepValue(HISTORICAL_SIGMA_YZ)+0.25 * this->GetGeometry()[2].FastGetSolutionStepValue(HISTORICAL_SIGMA_YZ)+0.25 * this->GetGeometry()[3].FastGetSolutionStepValue(HISTORICAL_SIGMA_YZ);

		    double SZX= 0.25 * this->GetGeometry()[0].FastGetSolutionStepValue(HISTORICAL_SIGMA_ZX) + 0.25 * this->GetGeometry()[1].FastGetSolutionStepValue(HISTORICAL_SIGMA_ZX)+0.25 * this->GetGeometry()[2].FastGetSolutionStepValue(HISTORICAL_SIGMA_ZX)+0.25 * this->GetGeometry()[3].FastGetSolutionStepValue(HISTORICAL_SIGMA_ZX);
		    double SZY= 0.25 * this->GetGeometry()[0].FastGetSolutionStepValue(HISTORICAL_SIGMA_ZY) + 0.25 * this->GetGeometry()[1].FastGetSolutionStepValue(HISTORICAL_SIGMA_ZY)+0.25 * this->GetGeometry()[2].FastGetSolutionStepValue(HISTORICAL_SIGMA_ZY)+0.25 * this->GetGeometry()[3].FastGetSolutionStepValue(HISTORICAL_SIGMA_ZY);
		    double SZZ= 0.25 * this->GetGeometry()[0].FastGetSolutionStepValue(HISTORICAL_SIGMA_ZZ) + 0.25 * this->GetGeometry()[1].FastGetSolutionStepValue(HISTORICAL_SIGMA_ZZ)+0.25 * this->GetGeometry()[2].FastGetSolutionStepValue(HISTORICAL_SIGMA_ZZ)+0.25 * this->GetGeometry()[3].FastGetSolutionStepValue(HISTORICAL_SIGMA_ZZ);


		    //rRightHandSideVector[0] -= DN_DX(0,0)*this->GetGeometry()[0].FastGetSolutionStepValue(HISTORICAL_SIGMA_XX) + DN_DX(0,1)*this->GetGeometry()[0].FastGetSolutionStepValue(HISTORICAL_SIGMA_XY) + DN_DX(0,2)*this->GetGeometry()[0].FastGetSolutionStepValue(HISTORICAL_SIGMA_XZ);
		    rRightHandSideVector[0] -= DN_DX(0,0)*SXX + DN_DX(0,1)*SXY + DN_DX(0,2)*SXZ;

		    //0  +   dN1/dy*SigmaY + 0 + dN1/dx*TauXY
		    //rRightHandSideVector[1] -= DN_DX(0,0)*this->GetGeometry()[0].FastGetSolutionStepValue(HISTORICAL_SIGMA_YX) + DN_DX(0,1)*this->GetGeometry()[0].FastGetSolutionStepValue(HISTORICAL_SIGMA_YY) + DN_DX(0,2)*this->GetGeometry()[0].FastGetSolutionStepValue(HISTORICAL_SIGMA_YZ);
		    rRightHandSideVector[1] -= DN_DX(0,0)*SYX + DN_DX(0,1)*SYY + DN_DX(0,2)*SYZ;

		    //rRightHandSideVector[2] -= DN_DX(0,0)*this->GetGeometry()[0].FastGetSolutionStepValue(HISTORICAL_SIGMA_ZX) + DN_DX(0,1)*this->GetGeometry()[0].FastGetSolutionStepValue(HISTORICAL_SIGMA_ZY) + DN_DX(0,2)*this->GetGeometry()[0].FastGetSolutionStepValue(HISTORICAL_SIGMA_ZZ);
		    rRightHandSideVector[2] -= DN_DX(0,0)*SZX + DN_DX(0,1)*SZY + DN_DX(0,2)*SZZ;


		    //dN2/dx*SigmaX + 0 + dN2/dy*TauXY
		    //rRightHandSideVector[3] -= DN_DX(1,0)*this->GetGeometry()[1].FastGetSolutionStepValue(HISTORICAL_SIGMA_XX) + DN_DX(1,1)*this->GetGeometry()[1].FastGetSolutionStepValue(HISTORICAL_SIGMA_XY) + DN_DX(1,2)*this->GetGeometry()[1].FastGetSolutionStepValue(HISTORICAL_SIGMA_XZ);
		    rRightHandSideVector[3] -= DN_DX(1,0)*SXX + DN_DX(1,1)*SXY + DN_DX(1,2)*SXZ;

		    //0  +   dN2/dy*SigmaY + 0 + dN2/dx*TauXY
		    //rRightHandSideVector[4] -= DN_DX(1,0)*this->GetGeometry()[1].FastGetSolutionStepValue(HISTORICAL_SIGMA_YX) + DN_DX(1,1)*this->GetGeometry()[1].FastGetSolutionStepValue(HISTORICAL_SIGMA_YY) + DN_DX(1,2)*this->GetGeometry()[1].FastGetSolutionStepValue(HISTORICAL_SIGMA_YZ) ;
		    rRightHandSideVector[4] -= DN_DX(1,0)*SYX + DN_DX(1,1)*SYY + DN_DX(1,2)*SYZ;


		    //rRightHandSideVector[5] -= DN_DX(1,0)*this->GetGeometry()[1].FastGetSolutionStepValue(HISTORICAL_SIGMA_ZX) + DN_DX(1,1)*this->GetGeometry()[1].FastGetSolutionStepValue(HISTORICAL_SIGMA_ZY) + DN_DX(1,2)*this->GetGeometry()[1].FastGetSolutionStepValue(HISTORICAL_SIGMA_ZZ) ;
		    rRightHandSideVector[5] -= DN_DX(1,0)*SZX + DN_DX(1,1)*SZY + DN_DX(1,2)*SZZ;



		    //dN3/dx*SigmaX + 0 + dN3/dy*TauXY
		    //rRightHandSideVector[6] -= DN_DX(2,0)*this->GetGeometry()[2].FastGetSolutionStepValue(HISTORICAL_SIGMA_XX) + DN_DX(2,1)*this->GetGeometry()[2].FastGetSolutionStepValue(HISTORICAL_SIGMA_XY) + DN_DX(2,2)*this->GetGeometry()[2].FastGetSolutionStepValue(HISTORICAL_SIGMA_XZ) ;
		    rRightHandSideVector[6] -= DN_DX(2,0)*SXX + DN_DX(2,1)*SXY + DN_DX(2,2)*SXZ ;

		    //0  +   dN3/dy*SigmaY + 0 + dN3/dx*TauXY
		    //rRightHandSideVector[7] -= DN_DX(2,0)*this->GetGeometry()[2].FastGetSolutionStepValue(HISTORICAL_SIGMA_YX) + DN_DX(2,1)*this->GetGeometry()[2].FastGetSolutionStepValue(HISTORICAL_SIGMA_YY) + DN_DX(2,2)*this->GetGeometry()[2].FastGetSolutionStepValue(HISTORICAL_SIGMA_YZ) ;
		    rRightHandSideVector[7] -= DN_DX(2,0)*SYX + DN_DX(2,1)*SYY + DN_DX(2,2)*SYZ ;


		    //rRightHandSideVector[8] -= DN_DX(2,0)*this->GetGeometry()[2].FastGetSolutionStepValue(HISTORICAL_SIGMA_ZX) + DN_DX(2,1)*this->GetGeometry()[2].FastGetSolutionStepValue(HISTORICAL_SIGMA_ZY) + DN_DX(2,2)*this->GetGeometry()[2].FastGetSolutionStepValue(HISTORICAL_SIGMA_ZZ) ;
		    rRightHandSideVector[8] -= DN_DX(2,0)*SZX + DN_DX(2,1)*SZY + DN_DX(2,2)*SZZ;

		     //dN3/dx*SigmaX + 0 + dN3/dy*TauXY
		    //rRightHandSideVector[9] -= DN_DX(3,0)*this->GetGeometry()[3].FastGetSolutionStepValue(HISTORICAL_SIGMA_XX) + DN_DX(3,1)*this->GetGeometry()[3].FastGetSolutionStepValue(HISTORICAL_SIGMA_XY) + DN_DX(3,2)*this->GetGeometry()[3].FastGetSolutionStepValue(HISTORICAL_SIGMA_XZ) ;
		    rRightHandSideVector[9] -= DN_DX(3,0)*SXX + DN_DX(3,1)*SXY + DN_DX(3,2)*SXZ ;

		    //0  +   dN3/dy*SigmaY + 0 + dN3/dx*TauXY
		    //rRightHandSideVector[10] -= DN_DX(3,0)*this->GetGeometry()[3].FastGetSolutionStepValue(HISTORICAL_SIGMA_YX) + DN_DX(3,1)*this->GetGeometry()[3].FastGetSolutionStepValue(HISTORICAL_SIGMA_YY) + DN_DX(3,2)*this->GetGeometry()[3].FastGetSolutionStepValue(HISTORICAL_SIGMA_YZ) ;
		    rRightHandSideVector[10] -= DN_DX(3,0)*SYX + DN_DX(3,1)*SYY + DN_DX(3,2)*SYZ ;

		    //rRightHandSideVector[11] -= DN_DX(3,0)*this->GetGeometry()[3].FastGetSolutionStepValue(HISTORICAL_SIGMA_ZX) + DN_DX(3,1)*this->GetGeometry()[3].FastGetSolutionStepValue(HISTORICAL_SIGMA_ZY) + DN_DX(3,2)*this->GetGeometry()[3].FastGetSolutionStepValue(HISTORICAL_SIGMA_ZZ) ;
		    rRightHandSideVector[11] -= DN_DX(3,0)*SZX + DN_DX(3,1)*SZY + DN_DX(3,2)*SZZ ;

		    rRightHandSideVector*=Area;


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
    void CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo) override
{
    KRATOS_TRY
        const unsigned int LocalSize = TDim * TNumNodes;

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

/*    const double& density = 0.333333333*(this->GetGeometry()[0].FastGetSolutionStepValue(DENSITY)+
                                         this->GetGeometry()[1].FastGetSolutionStepValue(DENSITY) +
                                         this->GetGeometry()[2].FastGetSolutionStepValue(DENSITY));
    //lumped
    unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();
    unsigned int NumberOfNodes = this->GetGeometry().size();

    double mA0 = GeometryUtils::CalculateVolume2D(this->GetGeometry());

    if(rMassMatrix.size1() != 6)
        rMassMatrix.resize(6,6,false);

    noalias(rMassMatrix) = ZeroMatrix(6,6);

    double nodal_mass = mA0 * density * 0.333333333333333333;

    for(unsigned int i=0; i<NumberOfNodes; i++)
    {
        for(unsigned int j=0; j<dimension; j++)
        {
            unsigned int index = i*dimension + j;
            rMassMatrix(index,index) = nodal_mass;
        }
    }

    */
    //KRATOS_WATCH(rMassMatrix)

    KRATOS_CATCH("")
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
            const ProcessInfo& rCurrentProcessInfo) override
 {
    KRATOS_TRY
/*    if (rDampingMatrix.size1() != 6)
        rDampingMatrix.resize(6, 6, false);

    noalias(rDampingMatrix) = ZeroMatrix(6, 6);*/

    const unsigned int LocalSize = TDim * TNumNodes;

    // Resize and set to zero the matrix
    // Note that we don't clean the RHS because it will already contain body force (and stabilization) contributions
    if (rDampingMatrix.size1() != LocalSize)
        rDampingMatrix.resize(LocalSize, LocalSize, false);

    noalias(rDampingMatrix) = ZeroMatrix(LocalSize, LocalSize);


    if(TDim==2)
    {
	    //fill in the damping matrix
	    BoundedMatrix<double,3,6> msB;
	    BoundedMatrix<double,3,3> ms_constitutive_matrix;
	    BoundedMatrix<double,3,6> ms_temp;

	    BoundedMatrix<double,3,2> msDN_Dx;
	    array_1d<double,3> msN; //dimension = number of nodes

	    unsigned int NumberOfNodes = this->GetGeometry().size();
	    unsigned int dim = this->GetGeometry().WorkingSpaceDimension();

	    //getting data for the given geometry
	    double Area;
	    GeometryUtils::CalculateGeometryData(this->GetGeometry(), msDN_Dx, msN, Area);


	    double NU = this->GetProperties()[POISSON_RATIO];//  GetProperties()[POISSON_RATIO];
	    double E = this->GetProperties()[YOUNG_MODULUS];



	    double dt = rCurrentProcessInfo[DELTA_TIME];
	    //Lame constants. note that for the QFLUIDelastic solid the Lame constants must be multiplied by the timestep
	    //const
	    double MU=0.5*dt*E/(1.0+NU);
	    //const
	    double LAMBDA=NU*E*dt/((1.0+NU)*(1.0-2.0*NU));
	    //const
	    double KAPPA=LAMBDA+0.6666666*MU;

	    //SHEAR CONTRIBUTION TO THE "DAMPING" MATRIX
	    for (unsigned int i=0; i<NumberOfNodes; i++)
	    {
		unsigned int index = dim*i;
		msB(0,index+0)=msDN_Dx(i,0);
		msB(0,index+1)= 0.0;
		msB(1,index+0)=0.0;
		msB(1,index+1)= msDN_Dx(i,1);
		msB(2,index+0)= msDN_Dx(i,1);
		msB(2,index+1)= msDN_Dx(i,0);
	    }

	    //constitutive tensor
	    ms_constitutive_matrix(0,0) = (4.0/3.0)*MU;
	    ms_constitutive_matrix(0,1) = -2.0/3.0*MU;
	    ms_constitutive_matrix(0,2) = 0.0;
	    ms_constitutive_matrix(1,0) = -2.0/3.0*MU;
	    ms_constitutive_matrix(1,1) = 4.0/3.0*MU;
	    ms_constitutive_matrix(1,2) = 0.0;
	    ms_constitutive_matrix(2,0) = 0.0;
	    ms_constitutive_matrix(2,1) = 0.0;
	    ms_constitutive_matrix(2,2) = MU;

	    //calculating viscous contributions
	    ms_temp = prod( ms_constitutive_matrix , msB);
	    noalias(rDampingMatrix) = prod( trans(msB) , ms_temp);


	    //now we reuse the constitutive tensor
	    ms_constitutive_matrix(0,0) = KAPPA;
	    ms_constitutive_matrix(0,1) = KAPPA ;
	    ms_constitutive_matrix(0,2) = 0.0;
	    ms_constitutive_matrix(1,0) = KAPPA;
	    ms_constitutive_matrix(1,1) = KAPPA;
	    ms_constitutive_matrix(1,2) = 0.0;
	    ms_constitutive_matrix(2,0) = 0.0;
	    ms_constitutive_matrix(2,1) = 0.0;
	    ms_constitutive_matrix(2,2) = 0.0;

	    //calculating volumetric contribution
	    ms_temp = prod( ms_constitutive_matrix , msB);
	    //ms_temp*=dt;
	    rDampingMatrix+= prod( trans(msB) , ms_temp);

	    rDampingMatrix *= Area;

	    //Now calculate an additional contribution to the residual: r -= rDampingMatrix * (v)
	    array_1d< double, 6 > Vel;
	    Vel[0] = this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_X);
	    Vel[1] = this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_Y);
	    Vel[2] = this->GetGeometry()[1].FastGetSolutionStepValue(VELOCITY_X);
	    Vel[3] = this->GetGeometry()[1].FastGetSolutionStepValue(VELOCITY_Y);
	    Vel[4] = this->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY_X);
	    Vel[5] = this->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY_Y);

	    noalias(rRightHandSideVector) -= prod(rDampingMatrix, Vel);
    }
    else{

    	    //fill in the damping matrix
	    BoundedMatrix<double,6,12> msB = ZeroMatrix(6,12);
	    BoundedMatrix<double,6,6> ms_constitutive_matrix;
	    BoundedMatrix<double,6,12> ms_temp;

	    BoundedMatrix<double,4,3> msDN_Dx;
    	    array_1d<double,4> msN; //dimension = number of nodes

	    unsigned int NumberOfNodes = this->GetGeometry().size();
	    unsigned int dim = this->GetGeometry().WorkingSpaceDimension();

	    //getting data for the given geometry
	    double Area;
	    GeometryUtils::CalculateGeometryData(this->GetGeometry(), msDN_Dx, msN, Area);


	    double NU = this->GetProperties()[POISSON_RATIO];//  GetProperties()[POISSON_RATIO];
	    double E = this->GetProperties()[YOUNG_MODULUS];

	    /*KRATOS_WATCH(NU)
	    KRATOS_WATCH(E)*/



	    double dt = rCurrentProcessInfo[DELTA_TIME];
	    //Lame constants. note that for the QFLUIDelastic solid the Lame constants must be multiplied by the timestep
	    //const
	    double MU=0.5*dt*E/(1.0+NU);
	    //const
	    double LAMBDA=NU*E*dt/((1.0+NU)*(1.0-2.0*NU));
	    //const
	    double KAPPA=LAMBDA+0.6666666*MU;

	    //SHEAR CONTRIBUTION TO THE "DAMPING" MATRIX
	    for (unsigned int i=0; i<NumberOfNodes; i++)
	    {
		unsigned int start = dim*i;

		msB(0,start) =	msDN_Dx(i,0);
		msB(1,start+1)=	msDN_Dx(i,1);
		msB(2,start+2)= msDN_Dx(i,2);
		msB(3,start) =	msDN_Dx(i,1);
		msB(3,start+1) = msDN_Dx(i,0);
		msB(4,start) =	msDN_Dx(i,2);
		msB(4,start+2) = msDN_Dx(i,0);
		msB(5,start+1)= msDN_Dx(i,2);
		msB(5,start+2) = msDN_Dx(i,1);
	    }

    ms_constitutive_matrix(0,0) = (4.0/3.0)*MU;
    ms_constitutive_matrix(0,1) = -2.0/3.0*MU;
    ms_constitutive_matrix(0,2) = -2.0/3.0*MU;
    ms_constitutive_matrix(0,3) = 0.0;
    ms_constitutive_matrix(0,4) = 0.0;
    ms_constitutive_matrix(0,5) = 0.0;

    ms_constitutive_matrix(1,0) = -2.0/3.0*MU;
    ms_constitutive_matrix(1,1) = 4.0/3.0*MU;
    ms_constitutive_matrix(1,2) = -2.0/3.0*MU;
    ms_constitutive_matrix(1,3) = 0.0;
    ms_constitutive_matrix(1,4) = 0.0;
    ms_constitutive_matrix(1,5) = 0.0;

    ms_constitutive_matrix(2,0) = -2.0/3.0*MU;
    ms_constitutive_matrix(2,1) = -2.0/3.0*MU;
    ms_constitutive_matrix(2,2) = 4.0/3.0*MU;
    ms_constitutive_matrix(2,3) = 0.0;
    ms_constitutive_matrix(2,4) = 0.0;
    ms_constitutive_matrix(2,5) = 0.0;

    ms_constitutive_matrix(3,0) = 0.0;
    ms_constitutive_matrix(3,1) = 0.0;
    ms_constitutive_matrix(3,2) = 0.0;
    ms_constitutive_matrix(3,3) = MU;
    ms_constitutive_matrix(3,4) = 0.0;
    ms_constitutive_matrix(3,5) = 0.0;

    ms_constitutive_matrix(4,0) = 0.0;
    ms_constitutive_matrix(4,1) = 0.0;
    ms_constitutive_matrix(4,2) = 0.0;
    ms_constitutive_matrix(4,3) = 0.0;
    ms_constitutive_matrix(4,4) = MU;
    ms_constitutive_matrix(4,5) = 0.0;

    ms_constitutive_matrix(5,0) = 0.0;
    ms_constitutive_matrix(5,1) = 0.0;
    ms_constitutive_matrix(5,2) = 0.0;
    ms_constitutive_matrix(5,3) = 0.0;
    ms_constitutive_matrix(5,4) = 0.0;
    ms_constitutive_matrix(5,5) = MU;

	    //calculating volumetric contribution
	    ms_temp = prod( ms_constitutive_matrix , msB);
	    //ms_temp*=dt;
	    rDampingMatrix+= prod( trans(msB) , ms_temp);

	    //rDampingMatrix *= Area;
	    //constitutive tensor
	    ms_constitutive_matrix(0,0) = KAPPA;
	    ms_constitutive_matrix(0,1) = KAPPA ;
	    ms_constitutive_matrix(0,2) = KAPPA;
	    ms_constitutive_matrix(0,3)= 0.0;
	    ms_constitutive_matrix(0,4) = 0.0 ;
	    ms_constitutive_matrix(0,5) = 0.0;

	    ms_constitutive_matrix(1,0) = KAPPA;
	    ms_constitutive_matrix(1,1) = KAPPA;
	    ms_constitutive_matrix(1,2) = KAPPA;
	    ms_constitutive_matrix(1,3) = 0.0;
	    ms_constitutive_matrix(1,4) = 0.0;
	    ms_constitutive_matrix(1,5) = 0.0;

	    ms_constitutive_matrix(2,0) = KAPPA;
	    ms_constitutive_matrix(2,1) = KAPPA;
	    ms_constitutive_matrix(2,2) = KAPPA;
	    ms_constitutive_matrix(2,3) = 0.0;
	    ms_constitutive_matrix(2,4) = 0.0;
	    ms_constitutive_matrix(2,5) = 0.0;

	    ms_constitutive_matrix(3,0) = 0.0;
	    ms_constitutive_matrix(3,1) = 0.0;
	    ms_constitutive_matrix(3,2) = 0.0;
	    ms_constitutive_matrix(3,3) = 0.0;
	    ms_constitutive_matrix(3,4) = 0.0;
	    ms_constitutive_matrix(3,5) = 0.0;

	    ms_constitutive_matrix(4,0) = 0.0;
	    ms_constitutive_matrix(4,1) = 0.0;
	    ms_constitutive_matrix(4,2) = 0.0;
	    ms_constitutive_matrix(4,3) = 0.0;
	    ms_constitutive_matrix(4,4) = 0.0;
	    ms_constitutive_matrix(4,5) = 0.0;

	    ms_constitutive_matrix(5,0) = 0.0;
	    ms_constitutive_matrix(5,1) = 0.0;
	    ms_constitutive_matrix(5,2) = 0.0;
	    ms_constitutive_matrix(5,3) = 0.0;
	    ms_constitutive_matrix(5,4) = 0.0;
	    ms_constitutive_matrix(5,5) = 0.0;

   	    //calculating volumetric contribution
	    ms_temp = prod( ms_constitutive_matrix , msB);
	    //ms_temp*=dt;
	    rDampingMatrix+= prod( trans(msB) , ms_temp);

	    rDampingMatrix *= Area;

	    //Now calculate an additional contribution to the residual: r -= rDampingMatrix * (v)
	    /*array_1d< double, 12 > Vel;
	    Vel[0]  = this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_X);
	    Vel[1]  = this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_Y);
	    Vel[2]  = this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_Z);

	    Vel[3]  = this->GetGeometry()[1].FastGetSolutionStepValue(VELOCITY_X);
	    Vel[4]  = this->GetGeometry()[1].FastGetSolutionStepValue(VELOCITY_Y);
	    Vel[5]  = this->GetGeometry()[1].FastGetSolutionStepValue(VELOCITY_Z);

	    Vel[6]  = this->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY_X);
	    Vel[7]  = this->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY_Y);
	    Vel[8]  = this->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY_Z);

	    Vel[9]  = this->GetGeometry()[3].FastGetSolutionStepValue(VELOCITY_X);
	    Vel[10] = this->GetGeometry()[3].FastGetSolutionStepValue(VELOCITY_Y);
	    Vel[11] = this->GetGeometry()[3].FastGetSolutionStepValue(VELOCITY_Z);

	    noalias(rRightHandSideVector) -= prod(rDampingMatrix, Vel);*/

    }

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
        /*U[LocalIndex] = this->GetGeometry()[iNode].FastGetSolutionStepValue(PRESSURE); // Pressure Dof
        ++LocalIndex;*/
    }

    noalias(rRightHandSideVector) -= prod(rDampingMatrix, U);


    //KRATOS_WATCH(rDampingMatrix)
    //KRATOS_WATCH(rRightHandSideVector)
    KRATOS_CATCH("")
}
    void FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) override
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
        if (rVariable == NODAL_AREA)
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
                this->GetGeometry()[i].FastGetSolutionStepValue(NODAL_MASS) += Area * N[i];

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
/*    void Calculate(const Variable<array_1d<double, 3 > >& rVariable,
                           array_1d<double, 3 > & rOutput,
                           const ProcessInfo& rCurrentProcessInfo)
{


    virtual void Calculate(const Variable<Vector >& rVariable,
                           Vector& Output,
                           const ProcessInfo& rCurrentProcessInfo)
    {
    }

    virtual void Calculate(const Variable<Matrix >& rVariable,
                           Matrix& Output,
                           const ProcessInfo& rCurrentProcessInfo)

                           */
    virtual void Calculate(const Variable<Matrix >& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo) override
     {

  if(rVariable == CAUCHY_STRESS_TENSOR)
    {
    if(TDim==2)
    {
	    BoundedMatrix<double,2,2> CauchyStress=ZeroMatrix(2,2);
	    BoundedMatrix<double,2,2> HistoricalCauchyStress=ZeroMatrix(2,2);
	    //KRATOS_THROW_ERROR(std::logic_error,"not dereeeeeeeeeeeeeeeeeeeeeeecha","");
	    //noalias(rCauchyStress) = ZeroMatrix(2, 2);

	    const array_1d<double,3>& v0 = this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
	    const array_1d<double,3>& v1 = this->GetGeometry()[1].FastGetSolutionStepValue(VELOCITY);
	    const array_1d<double,3>& v2 = this->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY);

	    double Area;
	    BoundedMatrix<double,3,2> msDN_Dx;
	    array_1d<double,3> msN; //dimension = number of nodes
	    GeometryUtils::CalculateGeometryData(this->GetGeometry(), msDN_Dx, msN, Area);

	    double NU = this->GetProperties()[POISSON_RATIO];
	    double E = this->GetProperties()[YOUNG_MODULUS];

	    double dt = rCurrentProcessInfo[DELTA_TIME];
	    //Lame constants. note that for the QFLUIDelastic solid the Lame constants must be multiplied by the timestep
	    //const
	    double MU=0.5*dt*E/(1.0+NU);
	    //const
	    double LAMBDA=NU*E*dt/((1.0+NU)*(1.0-2.0*NU));
	    //const
	    double KAPPA=LAMBDA+0.6666666*MU;


	    CauchyStress(0,0)=2.0* (msDN_Dx(0,0)*v0[0]+msDN_Dx(1,0)*v1[0]+msDN_Dx(2,0)*v2[0]);
	    CauchyStress(0,1)=msDN_Dx(0,1)*v0[0]+msDN_Dx(1,1)*v1[0]+msDN_Dx(2,1)*v2[0] + msDN_Dx(0,0)*v0[1]+msDN_Dx(1,0)*v1[1]+msDN_Dx(2,0)*v2[1];

	    CauchyStress(1,0)=msDN_Dx(0,1)*v0[0]+msDN_Dx(1,1)*v1[0]+msDN_Dx(2,1)*v2[0] + msDN_Dx(0,0)*v0[1]+msDN_Dx(1,0)*v1[1]+msDN_Dx(2,0)*v2[1];
	    CauchyStress(1,1)=2.0*(msDN_Dx(0,1)*v0[1]+msDN_Dx(1,1)*v1[1]+msDN_Dx(2,1)*v2[1]);

	    CauchyStress*=MU;

	    //adding the volumetric part
	    double div_v = msDN_Dx(0,0)*v0[0] + msDN_Dx(0,1)*v0[1];
	    div_v+=       msDN_Dx(1,0)*v1[0] + msDN_Dx(1,1)*v1[1];
	    div_v+=	  msDN_Dx(2,0)*v2[0] + msDN_Dx(2,1)*v2[1];

	    CauchyStress(0,0)+=KAPPA*div_v;
	    CauchyStress(1,1)+=KAPPA*div_v;

	    this->SetValue(CAUCHY_STRESS_TENSOR, CauchyStress);
	}
    else{

	double Area;
        array_1d<double, TNumNodes> N;
        BoundedMatrix<double, TNumNodes, TDim> msDN_Dx;
        GeometryUtils::CalculateGeometryData(this->GetGeometry(), msDN_Dx, N, Area);

        BoundedMatrix<double,3,3> CauchyStress=ZeroMatrix(3,3);
	//BoundedMatrix<double,3,3> HistoricalCauchyStress=ZeroMatrix(3,3);
	/*KRATOS_WATCH("222222222222222222222222222222222")
	KRATOS_WATCH("222222222222222222222222222222222")
	KRATOS_WATCH("222222222222222222222222222222222")
	KRATOS_WATCH("222222222222222222222222222222222")	*/
	//KRATOS_THROW_ERROR(std::logic_error,"not dereeeeeeeeeeeeeeeeeeeeeeecha","");
	//noalias(rCauchyStress) = ZeroMatrix(2, 2);

	const array_1d<double,3>& v0 = this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
	const array_1d<double,3>& v1 = this->GetGeometry()[1].FastGetSolutionStepValue(VELOCITY);
	const array_1d<double,3>& v2 = this->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY);
        const array_1d<double,3>& v3 = this->GetGeometry()[3].FastGetSolutionStepValue(VELOCITY);

        double NU = this->GetProperties()[POISSON_RATIO];
	double E = this->GetProperties()[YOUNG_MODULUS];
        //NU = 3.000000e-01;
        //E = 3.00000e+05;

	double dt = rCurrentProcessInfo[DELTA_TIME];
	//Lame constants. note that for the QFLUIDelastic solid the Lame constants must be multiplied by the timestep
	//const
	double MU=0.5*dt*E/(1.0+NU);

	//const
	double LAMBDA=NU*E*dt/((1.0+NU)*(1.0-2.0*NU));
	//const
	double KAPPA=LAMBDA+0.6666666*MU;
        //KRATOS_WATCH(CauchyStress)
        CauchyStress(0,0)=2.0* (msDN_Dx(0,0)*v0[0]+msDN_Dx(1,0)*v1[0]+msDN_Dx(2,0)*v2[0]+msDN_Dx(3,0)*v3[0]);
	CauchyStress(0,1)=msDN_Dx(0,1)*v0[0]+msDN_Dx(1,1)*v1[0]+msDN_Dx(2,1)*v2[0]+msDN_Dx(3,1)*v3[0]+msDN_Dx(0,0)*v0[1]+msDN_Dx(1,0)*v1[1]+msDN_Dx(2,0)*v2[1]+msDN_Dx(3,0)*v3[1];
	CauchyStress(0,2)=msDN_Dx(0,2)*v0[0]+msDN_Dx(1,2)*v1[0]+msDN_Dx(2,2)*v2[0]+msDN_Dx(3,2)*v3[0]+msDN_Dx(0,0)*v0[2]+msDN_Dx(1,0)*v1[2]+msDN_Dx(2,0)*v2[2]+msDN_Dx(3,0)*v3[2];


	CauchyStress(1,0)=CauchyStress(0,1);//msDN_Dx(0,1)*v0[0]+msDN_Dx(1,1)*v1[0]+msDN_Dx(2,1)*v2[0]+msDN_Dx(3,1)*v3[0] + msDN_Dx(0,0)*v0[1]+msDN_Dx(1,0)*v1[1]+msDN_Dx(2,0)*v2[1]+msDN_Dx(3,0)*v3[1];
	CauchyStress(1,1)=2.0*(msDN_Dx(0,1)*v0[1]+msDN_Dx(1,1)*v1[1]+msDN_Dx(2,1)*v2[1]+msDN_Dx(3,1)*v3[1]);
	CauchyStress(1,2)=msDN_Dx(0,2)*v0[1]+msDN_Dx(1,2)*v1[1]+msDN_Dx(2,2)*v2[1] +msDN_Dx(3,2)*v3[1] + msDN_Dx(0,1)*v0[2]+msDN_Dx(1,1)*v1[2]+msDN_Dx(2,1)*v2[2]+msDN_Dx(3,1)*v3[2];


	CauchyStress(2,0)=CauchyStress(0,2);//msDN_Dx(0,2)*v0[0]+msDN_Dx(1,2)*v1[0]+msDN_Dx(2,2)*v2[0]+msDN_Dx(3,2)*v3[0] + msDN_Dx(0,0)*v0[2]+msDN_Dx(1,0)*v1[2]+msDN_Dx(2,0)*v2[2]+msDN_Dx(3,0)*v3[2];
	CauchyStress(2,1)=CauchyStress(1,2);//CauchyStress(1,2);//msDN_Dx(0,2)*v0[1]+msDN_Dx(1,2)*v1[1]+msDN_Dx(2,2)*v2[1] +msDN_Dx(3,2)*v3[1] + msDN_Dx(0,1)*v0[2]+msDN_Dx(1,1)*v1[2]+msDN_Dx(2,1)*v2[2]+msDN_Dx(2,1)*v3[2];
	CauchyStress(2,2)=2.0*(msDN_Dx(0,2)*v0[2]+msDN_Dx(1,2)*v1[2]+msDN_Dx(2,2)*v2[2]+msDN_Dx(3,2)*v3[2]);

       CauchyStress*=MU;

       //KRATOS_WATCH(CauchyStress)

       //KRATOS_WATCH(CauchyStress)
       //adding the volumetric part
       double div_v = msDN_Dx(0,0)*v0[0] + msDN_Dx(0,1)*v0[1] + msDN_Dx(0,2)*v0[2];
       div_v+= msDN_Dx(1,0)*v1[0] + msDN_Dx(1,1)*v1[1] + msDN_Dx(1,2)*v1[2];
       div_v+=	msDN_Dx(2,0)*v2[0] + msDN_Dx(2,1)*v2[1] + msDN_Dx(2,2)*v2[2];
       div_v+=	msDN_Dx(3,0)*v3[0] + msDN_Dx(3,1)*v3[1] + msDN_Dx(3,2)*v3[2];

       CauchyStress(0,0)+=KAPPA*div_v;
       CauchyStress(1,1)+=KAPPA*div_v;
       CauchyStress(2,2)+=KAPPA*div_v;
	    //KRATOS_WATCH(CauchyStress)
       //KRATOS_WATCH(HistoricalCauchyStress)
       //HistoricalCauchyStress=this->GetValue(CAUCHY_STRESS_TENSOR);
       //KRATOS_WATCH(HistoricalCauchyStress)
       //CauchyStress+=HistoricalCauchyStress;
       //KRATOS_WATCH(CauchyStress)

       this->SetValue(CAUCHY_STRESS_TENSOR, CauchyStress);
       //KRATOS_WATCH(CauchyStress)

    }





    }
    else
       KRATOS_ERROR << "Wrong variable. Calculate function of QFLUIDelastic element is meant to compute Cauchy stress only." << std::endl;


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
                          const ProcessInfo& rCurrentProcessInfo) const override;

    /// Returns a list of the element's Dofs
    /**
     * @param ElementalDofList the list of DOFs
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetDofList(DofsVectorType& rElementalDofList,
                    const ProcessInfo& rCurrentProcessInfo) const override;

    /// Returns VELOCITY_X, VELOCITY_Y, (VELOCITY_Z,) PRESSURE for each node
    /**
     * @param Values Vector of nodal unknowns
     * @param Step Get result from 'Step' steps back, 0 is current step. (Must be smaller than buffer size)
     */
    void GetFirstDerivativesVector(Vector& Values, int Step = 0) const override;

    /// Returns ACCELERATION_X, ACCELERATION_Y, (ACCELERATION_Z,) 0 for each node
    /**
     * @param Values Vector of nodal second derivatives
     * @param Step Get result from 'Step' steps back, 0 is current step. (Must be smaller than buffer size)
     */
    void GetSecondDerivativesVector(Vector& Values, int Step = 0) const override;

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
    void CalculateOnIntegrationPoints(
        const Variable<array_1d<double, 3 > >& rVariable,
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
    void CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo) override
    {}

    /// Empty implementation of unused CalculateOnIntegrationPoints overloads to avoid compilation warning
    void CalculateOnIntegrationPoints(
        const Variable<array_1d<double, 6 > >& rVariable,
        std::vector<array_1d<double, 6 > >& rValues,
        const ProcessInfo& rCurrentProcessInfo) override
    {}

    /// Empty implementation of unused CalculateOnIntegrationPoints overloads to avoid compilation warning
    void CalculateOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues,
        const ProcessInfo& rCurrentProcessInfo) override
    {}

    /// Empty implementation of unused CalculateOnIntegrationPoints overloads to avoid compilation warning
    void CalculateOnIntegrationPoints(
        const Variable<Matrix>& rVariable,
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
    int Check(const ProcessInfo& rCurrentProcessInfo) const override
    {
        KRATOS_TRY

        // Perform basic element checks
        int ErrorCode = Kratos::Element::Check(rCurrentProcessInfo);
        if(ErrorCode != 0) return ErrorCode;

        // Checks on nodes

        // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
        for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
        {
            const auto &rNode = this->GetGeometry()[i];
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
                    KRATOS_ERROR << "Node " << this->GetGeometry()[i].Id() << "has non-zero Z coordinate." << std::endl;
            }
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
        buffer << "QFLUID #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "QFLUID" << TDim << "D";
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

    /// Add the momentum equation contribution to the RHS (body forces)
    virtual void AddMomentumRHS(VectorType& F,
                                const double Density,
                                const array_1d<double, TNumNodes>& rShapeFunc,
                                const double Weight)
    {
        double Coef = Density * Weight;

        array_1d<double, 3 > BodyForce = ZeroVector(3);
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

    /// Add lumped mass matrix
    /**
     * Adds the lumped mass matrix to an elemental LHS matrix. Note that the time factor
     * (typically 1/(k*Dt) ) is added by the scheme outside the element.
     * @param rLHSMatrix The local matrix where the result will be added
     * @param Mass The weight assigned to each node (typically Density * Area / NumNodes or Density*Volume / NumNodes)
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
           // ++DofIndex; // Skip pressure Dof
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


    /// Add a the contribution from a single integration point to the velocity contribution

    /// Assemble the contribution from an integration point to the element's residual.
    /** Note that the dynamic term is not included in the momentum equation.
     *  If OSS_SWITCH = 1, we don't take into account the 'dynamic' stabilization
     *  terms, as it they belong to the finite element space.
     */


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
     * @param rAdvVel Output array
     * @param rShapeFunc Shape functions evaluated at the point of interest
     */
    virtual void GetAdvectiveVel(array_1d< double, 3 > & rAdvVel,
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
     * @param rAdvVel Output array
     * @param rShapeFunc Shape functions evaluated at the point of interest
     * @param Step The time Step
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
     * @param rResult Output vector
     * @param rVelocity Velocity evaluated at the integration point
     * @param rShapeDeriv Derivatives of shape functions evaluated at the integration point
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
     * @param rResult The double where the value will be added to
     * @param rVariable The nodal variable to be read
     * @param rShapeFunc The values of the form functions in the point
     * @param Step The time Step (Defaults to 0 = Current)
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
     * @param rResult The double where the value will be added to
     * @param rVariable The nodal variable to be read
     * @param rShapeFunc The values of the form functions in the point
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
     * @see VMS2D, VMS3D for actual implementation
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
    //virtual void AddViscousTerm(MatrixType& rDampingMatrix,
    //                            const BoundedMatrix<double, TNumNodes, TDim >& rShapeDeriv,
    //                            const double Weight);

    /// Adds the contribution of the viscous term to the momentum equation (alternate).
    /**
     * This function is an alternate implementation of VMS::AddViscousTerm.
     * This version works with ublas matrices, using the relationship between stress and
     * rate of strain given by VMS::CalculateC. It is currently unused (as VMS::AddViscousTerm
     * is a more efficient implementation of the Cauchy equation) but it is left here so derived
     * classes can use it to implement other constitutive equations.
     * @param rDampingMatrix Elemental Damping matrix
     * @param rShapeDeriv Elemental shape function derivatives
     * @param Weight Effective viscosity, in dynamic units, weighted by the integration point area
     */


    /// Calculate the strain rate matrix
    /**
     * Unused, left to support derived classes. @see VMS::AddBTransCB
     * @param rB Strain rate matrix
     * @param rShapeDeriv Nodal shape funcion derivatives
     */
    //void CalculateB( BoundedMatrix<double, (TDim * TNumNodes) / 2, TDim * TNumNodes >& rB,                    const BoundedMatrix<double, TNumNodes, TDim >& rShapeDeriv);

    /// Calculate a matrix that provides the stress given the strain rate
    /**
     * Unused, left to support derived classes. @see VMS::AddBTransCB.
     * Note that only non-zero terms are written, so the output matrix should be
     * initialized before calling this.
     * @param rC Matrix representation of the stress tensor (output)
     * @param Viscosity Effective viscosity, in dynamic units, weighted by the integration point area
     */
    //virtual void CalculateC( BoundedMatrix<double, (TDim * TNumNodes) / 2, (TDim * TNumNodes) / 2 >& rC,                             const double Viscosity);

    double ConsistentMassCoef(const double Area);


    /*double SubscaleErrorEstimate(const ProcessInfo& rProcessInfo)
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
        array_1d< double, 3 > ElementalMomRes = ZeroVector(3);

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
        //array_1d< double, 3 > UGauss = ZeroVector(3);
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
    }*/

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
    QFLUID & operator=(QFLUID const& rOther);

    /// Copy constructor.
    QFLUID(QFLUID const& rOther);

    ///@}

}; // Class QFLUID

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
                                 QFLUID<TDim, TNumNodes>& rThis)
{
    return rIStream;
}

/// output stream function
template< unsigned int TDim,
          unsigned int TNumNodes >
inline std::ostream& operator <<(std::ostream& rOStream,
                                 const QFLUID<TDim, TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} // Fluid Dynamics Application group

} // namespace Kratos.

#endif // KRATOS_VMS_H_INCLUDED  defined
