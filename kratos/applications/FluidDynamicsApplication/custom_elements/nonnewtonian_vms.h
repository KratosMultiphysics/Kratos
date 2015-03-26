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

#ifndef KRATOS_NONNEWTONIAN_VMS2D_H
#define	KRATOS_NONNEWTONIAN_VMS2D_H

// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "vms.h"

// Application includes
#include "fluid_dynamics_application_variables.h"

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

/// An element for the simulation of Bingham plastics.
/**
 * This class implements a stabilized formulation based on the
 * Variational Multiscale framework to simulate a non-Newtonian fluid using the
 * Bingham plasic model. The the subscales can be modeled
 * using either Algebraic Subgird Scales (ASGS) or Orthogonal Subscales (OSS).
 * In the case of OSS, the projection terms are treated explicitly (computed
 * using the results of the previous iteration) and the subscales are not
 * tracked in time. The choice of subscale model is made based on the ProcessInfo
 * variable OSS_SWITCH (OSS if 1, ASGS otherwise).
 * This class implements both the 2D and 3D versions of the element.
 *
 * Note that the terms in the stabilized momentum equation that depend on gradients
 * of viscosity are not taken into account, as in this model viscosity is constant
 * inside each element.
 *
 * This class requires at least the following variables:\n
 * On each Node, as solution step variables VELOCITY, PRESSURE, ACCELERATION, MESH_VELOCITY, DENSITY, VISCOSITY.\n
 * On ProcessInfo OSS_SWITCH, DYNAMIC_TAU, DELTA_TIME, YIELD_STRESS.\n
 * If OSS is used, the nodes also require NODAL_AREA, ADVPROJ and DIVPROJ as solution step variables.\n
 * To compute the effective viscosity, MU and TAU have to be defined on the elements.\n
 * Some additional variables can be used to print results on the element: TAUONE, TAUTWO, TAU, MU, VORTICITY.
 *
 * @see ResidualBasedEliminationBuilderAndSolver compatible monolithic solution strategy.
 * @see PressureSplittingBuilderAndSolver compatible segregated solution strategy.
 * @see TrilinosPressureSplittingBuilderAndSolver compatible mpi strategy.
 * @see ResidualBasedPredictorCorrectorVelocityBossakScheme time scheme that can use
 * OSS stabilization.
 */
template< unsigned int TDim,
          unsigned int TNumNodes = TDim + 1 >
class NonnewtonianVMS : public VMS<TDim,TNumNodes>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of VMS
    KRATOS_CLASS_POINTER_DEFINITION(NonnewtonianVMS);

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

    typedef VMS<TDim,TNumNodes> BaseElementType;

    ///@}
    ///@name Life Cycle
    ///@{

    //Constructors.

    /// Default constuctor.
    /**
     * @param NewId Index number of the new element (optional)
     */
    NonnewtonianVMS(IndexType NewId = 0) :
        BaseElementType(NewId)
    {}

    ///Constructor using an array of nodes.
    /**
     * @param NewId Index of the new element
     * @param ThisNodes An array containing the nodes of the new element
     */
    NonnewtonianVMS(IndexType NewId, const NodesArrayType& ThisNodes) :
        BaseElementType(NewId, ThisNodes)
    {}

    /// Constructor using a geometry object.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     */
    NonnewtonianVMS(IndexType NewId, GeometryType::Pointer pGeometry) :
        BaseElementType(NewId, pGeometry)
    {}

    /// Constuctor using geometry and properties.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     * @param pProperties Pointer to the element's properties
     */
    NonnewtonianVMS(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) :
        BaseElementType(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    virtual ~NonnewtonianVMS()
    {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /// Create a new element of this type
    /**
     * Returns a pointer to a new NonnewtonianVMS element, created using given input
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
                            PropertiesType::Pointer pProperties) const
    {
        return Element::Pointer(new NonnewtonianVMS(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
    }

    ///@}
    ///@name Access
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
        buffer << "NonnewtonianVMS #" << this->Id();
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "NonnewtonianVMS" << TDim << "D";
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


    /// Calculate the viscosity, as given by the Bingham model.
    /**
     * @param Density Fluid density evaluated at the integration point
     * @param MolecularViscosity Viscosity of the fluid, in kinematic units (m2/s)
     * @param rShapeFunc Elemental shape functions, evaluated on the integration point
     * @param rShapeDeriv Shape function derivatives, evaluated on the integration point
     * @param TotalViscosity Effective viscosity (output)
     * @param rCurrentProcessInfo ProcessInfo instance (Checked for YIELD_STRESS)
     */
    virtual double EffectiveViscosity(double Density,
                                      const array_1d< double, TNumNodes > &rN,
                                      const boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim > &rDN_DX,
                                      double ElemSize,
                                      const ProcessInfo &rProcessInfo)
    {
	KRATOS_TRY
	double mu = this->GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY);
	double m_coef = rProcessInfo[M];
	double app_mu = 0.0;
	double aux_1;
	double gamma_dot = this->EquivalentStrainRate(rDN_DX);
	
	unsigned int nodes_number = 3;
   	double yield = 0.0;
    	double friction_angle_tangent = 0.0;
	double cohesion = 0.0;

	for (unsigned int ii = 0; ii < nodes_number; ++ii)
    	{
		friction_angle_tangent += this->GetGeometry()[ii].FastGetSolutionStepValue(INTERNAL_FRICTION_ANGLE);
		cohesion +=  this->GetGeometry()[ii].FastGetSolutionStepValue(YIELD_STRESS);//this is the COHESION of Mohr Coulomb failure criteria!!!
		if(this->GetGeometry()[ii].FastGetSolutionStepValue(PRESSURE) > 0.0)  //nodes with zero solid pressure does not have  any yield. Otherwise negative yield values appears (if water pressure > 0)
		{
		    yield +=  this->GetGeometry()[ii].FastGetSolutionStepValue(PRESSURE);
// 		    if(GetGeometry()[ii].FastGetSolutionStepValue(WATER_PRESSURE) >= 0.0){
// 			water_pressure +=  GetGeometry()[ii].FastGetSolutionStepValue(WATER_PRESSURE);
// 		    }
		}
	}
	friction_angle_tangent /= nodes_number;
	yield /= nodes_number;
    	cohesion /= nodes_number;
	yield *= friction_angle_tangent;
    	yield += cohesion;

	////////EXPONENCIAL MODEL
	if (gamma_dot > 1e-8)
	{
		aux_1 = 1.0 - exp(-(m_coef * gamma_dot));
		app_mu = mu + (yield / gamma_dot) * aux_1;
		// 			gamma_dot_inv = 1.0/gamma_dot;
		if (app_mu < mu)
		{
		    KRATOS_ERROR(std::logic_error, "!!!!!!!!!!!  APPARENT VISCOSITY < VISCOSITY !!!!!!!!", this->Id());
		}
	}
	else
	{
		app_mu = mu + yield * m_coef ;
		// // 			gamma_dot_inv = 0.0;
	}

        return app_mu;
	KRATOS_CATCH("")
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

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseElementType);
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseElementType);
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
    NonnewtonianVMS & operator=(NonnewtonianVMS const& rOther);

    /// Copy constructor.
    NonnewtonianVMS(NonnewtonianVMS const& rOther);

    ///@}

}; // Class NonnewtonianVMS

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
                                 NonnewtonianVMS<TDim, TNumNodes>& rThis)
{
    return rIStream;
}

/// output stream function
template< unsigned int TDim,
          unsigned int TNumNodes >
inline std::ostream& operator <<(std::ostream& rOStream,
                                 const NonnewtonianVMS<TDim, TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} // Fluid Dynamics Application group

} // namespace Kratos.

#endif	/* KRATOS_BINGHAM_VMS_H */

