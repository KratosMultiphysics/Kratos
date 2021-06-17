//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#ifndef KRATOS_TWO_FLUIDS_NAVIER_STOKES_WALL_CONDITION_H
#define KRATOS_TWO_FLUIDS_NAVIER_STOKES_WALL_CONDITION_H

// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
// #include "includes/define.h"
// #include "includes/condition.h"
// #include "includes/model_part.h"
// #include "includes/serializer.h"
// #include "includes/process_info.h"

// Application includes
#include "custom_conditions/navier_stokes_wall_condition.h"

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

/// Implements a wall condition for the Navier-Stokes monolithic formulation.
/**
  It is intended to be used in combination with ASGS Navier-Stokes symbolic elements or their
  derived classes and the ResidualBasedIncrementalUpdateStaticSchemeSlip time scheme, which supports
  slip conditions.
  @see NavierStokes,EmbeddedNavierStokes,ResidualBasedIncrementalUpdateStaticSchemeSlip
 */
template< unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) TwoFluidNavierStokesWallCondition : public NavierStokesWallCondition<TDim, TNumNodes>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TwoFluidNavierStokesWallCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(TwoFluidNavierStokesWallCondition);

    typedef NavierStokesWallCondition<TDim, TNumNodes> BaseType;

    typedef typename BaseType::ConditionDataStruct ConditionDataStruct;

    typedef Node<3> NodeType;

    typedef Properties PropertiesType;

    typedef Geometry<NodeType> GeometryType;

    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

    typedef Vector VectorType;

    typedef Matrix MatrixType;

    typedef std::size_t IndexType;

    typedef std::vector<std::size_t> EquationIdVectorType;

    typedef std::vector< Dof<double>::Pointer > DofsVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /** Admits an Id as a parameter.
      @param NewId Index for the new         // Struct to pass around the data
        ElementDataStruct data;
        this->FillElementData(data, rCurrentProcessInfo);condition
      */
    TwoFluidNavierStokesWallCondition(IndexType NewId = 0)
        : NavierStokesWallCondition<TDim, TNumNodes>(NewId)
    {
    }

    /// Constructor using an array of nodes
    /**
     @param NewId Index of the new condition
     @param ThisNodes An array containing the nodes of the new condition
     */
    TwoFluidNavierStokesWallCondition(
        IndexType NewId,
        const NodesArrayType& ThisNodes)
        : NavierStokesWallCondition<TDim, TNumNodes>(NewId, ThisNodes)
    {
    }

    /// Constructor using Geometry
    /**
     @param NewId Index of the new condition
     @param pGeometry Pointer to a geometry object
     */
    TwoFluidNavierStokesWallCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : NavierStokesWallCondition<TDim, TNumNodes>(NewId, pGeometry)
    {
    }

    /// Constructor using Properties
    /**
     @param NewId Index of the new element
     @param pGeometry Pointer to a geometry object
     @param pProperties Pointer to the element's properties
     */
    TwoFluidNavierStokesWallCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : NavierStokesWallCondition<TDim, TNumNodes>(NewId, pGeometry, pProperties)
    {
    }

    /// Copy constructor.
    TwoFluidNavierStokesWallCondition(TwoFluidNavierStokesWallCondition const& rOther):
        NavierStokesWallCondition<TDim, TNumNodes>(rOther)
    {
    }

    /// Destructor.
    ~TwoFluidNavierStokesWallCondition() override {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator
    TwoFluidNavierStokesWallCondition & operator=(TwoFluidNavierStokesWallCondition const& rOther)
    {
        Condition::operator=(rOther);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /// Create a new TwoFluidNavierStokesWallCondition object.
    /**
      @param NewId Index of the new condition
      @param ThisNodes An array containing the nodes of the new condition
      @param pProperties Pointer to the element's properties
      */
    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<TwoFluidNavierStokesWallCondition>(NewId, BaseType::GetGeometry().Create(ThisNodes), pProperties);
    }

    /// Create a new TwoFluidNavierStokesWallCondition object.
    /**
      @param NewId Index of the new condition
      @param pGeom A pointer to the condition's geometry
      @param pProperties Pointer to the element's properties
      */
    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive< TwoFluidNavierStokesWallCondition >(NewId, pGeom, pProperties);
    }

    /**
     * Clones the selected element variables, creating a new one
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @return a Pointer to the new element
     */
    Condition::Pointer Clone(
        IndexType NewId,
        NodesArrayType const& rThisNodes) const override
    {
        Condition::Pointer pNewCondition = Create(NewId, BaseType::GetGeometry().Create( rThisNodes ), BaseType::pGetProperties() );

        pNewCondition->SetData(this->GetData());
        pNewCondition->SetFlags(this->GetFlags());

        return pNewCondition;
    }

    /**
     * @brief Condition check
     * Condition check. Derived from the base Navier-Stokes condition to check if the viscosity is a nodal variable
     * @param rCurrentProcessInfo Reference to the ProcessInfo container
     * @return int 0 if successful
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

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
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "TwoFluidNavierStokesWallCondition" << TDim << "D";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "TwoFluidNavierStokesWallCondition";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {}

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

    /**
     * @brief Computes the right-hand side of the Navier slip contribution as e.g. described in BEHR2004
     * The (Navier) slip length is read as a nodal variable.
     * If a smaller value is set, tangential velocities lead to a higher tangential traction.
     * Though only tangential velocities should appear, a tangetial projection is added.
     * (Reference BEHR2004: https://onlinelibrary.wiley.com/doi/abs/10.1002/fld.663)
     * @param rRightHandSideVector reference to the RHS vector
     * @param rDataStruct reference to a struct to hand over data
     */
    void ComputeGaussPointNavierSlipRHSContribution(
        array_1d<double,TNumNodes*(TDim+1)>& rRightHandSideVector,
        const ConditionDataStruct& rDataStruct) override;

    /**
     * @brief Computes the left-hand side of the Navier slip contribution as e.g. described in BEHR2004
     * The (Navier) slip length is read as a nodal variable.
     * If a smaller value is set, tangential velocities lead to a higher tangential traction.
     * Though only tangential velocities should appear, a tangetial projection is added.
     * (Reference BEHR2004: https://onlinelibrary.wiley.com/doi/abs/10.1002/fld.663)
     * @param rLeftHandSideMatrix reference to the LHS matrix
     * @param rDataStruct reference to a struct to hand over data
     */
    void ComputeGaussPointNavierSlipLHSContribution(
        BoundedMatrix<double,TNumNodes*(TDim+1),TNumNodes*(TDim+1)>& rLeftHandSideMatrix,
        const ConditionDataStruct& rDataStruct ) override;

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition );
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition );
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


    ///@}
}; // Class TwoFluidNavierStokesWallCondition

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
template< unsigned int TDim, unsigned int TNumNodes >
inline std::istream& operator >> (std::istream& rIStream, TwoFluidNavierStokesWallCondition<TDim,TNumNodes>& rThis)
{
    return rIStream;
}

/// output stream function
template< unsigned int TDim, unsigned int TNumNodes >
inline std::ostream& operator << (std::ostream& rOStream, const TwoFluidNavierStokesWallCondition<TDim,TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}
///@} addtogroup block
}  // namespace Kratos.

#endif // KRATOS_TWO_FLUIDS_NAVIER_STOKES_WALL_CONDITION_H