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

#ifndef KRATOS_NAVIER_STOKES_WALL_CONDITION_H
#define KRATOS_NAVIER_STOKES_WALL_CONDITION_H

// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/condition.h"
#include "includes/model_part.h"
#include "includes/serializer.h"
#include "includes/process_info.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "includes/deprecated_variables.h"
#include "includes/cfd_variables.h"
#include "custom_utilities/fluid_element_utilities.h"

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
template< unsigned int TDim, unsigned int TNumNodes = TDim >
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) NavierStokesWallCondition : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NavierStokesWallCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(NavierStokesWallCondition);

    struct ConditionDataStruct
    {
        double wGauss;                                  // Gauss point weight
        double charVel;                                 // Problem characteristic velocity (used in the outlet inflow prevention)
        double delta;                                   // Non-dimensional positive sufficiently small constant (used in the outlet inflow prevention)
        array_1d<double, 3> Normal;                     // Condition normal
        array_1d<double, TNumNodes> N;                  // Gauss point shape functions values
        Vector ViscousStress;                           // Viscous stresses that are retrieved from parent
    };

    typedef Node < 3 > NodeType;

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
    NavierStokesWallCondition(IndexType NewId = 0):Condition(NewId)
    {
    }

    /// Constructor using an array of nodes
    /**
     @param NewId Index of the new condition
     @param ThisNodes An array containing the nodes of the new condition
     */
    NavierStokesWallCondition(IndexType NewId, const NodesArrayType& ThisNodes):
        Condition(NewId,ThisNodes)
    {
    }

    /// Constructor using Geometry
    /**
     @param NewId Index of the new condition
     @param pGeometry Pointer to a geometry object
     */
    NavierStokesWallCondition(IndexType NewId, GeometryType::Pointer pGeometry):
        Condition(NewId,pGeometry)
    {
    }

    /// Constructor using Properties
    /**
     @param NewId Index of the new element
     @param pGeometry Pointer to a geometry object
     @param pProperties Pointer to the element's properties
     */
    NavierStokesWallCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):
        Condition(NewId,pGeometry,pProperties)
    {
    }

    /// Copy constructor.
    NavierStokesWallCondition(NavierStokesWallCondition const& rOther):
        Condition(rOther)
    {
    }

    /// Destructor.
    ~NavierStokesWallCondition() override {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator
    NavierStokesWallCondition & operator=(NavierStokesWallCondition const& rOther)
    {
        Condition::operator=(rOther);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /// Create a new NavierStokesWallCondition object.
    /**
      @param NewId Index of the new condition
      @param ThisNodes An array containing the nodes of the new condition
      @param pProperties Pointer to the element's properties
      */
    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<NavierStokesWallCondition>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    /// Create a new NavierStokesWallCondition object.
    /**
      @param NewId Index of the new condition
      @param pGeom A pointer to the condition's geometry
      @param pProperties Pointer to the element's properties
      */
    Condition::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive< NavierStokesWallCondition >(NewId, pGeom, pProperties);
    }

    /**
     * Clones the selected element variables, creating a new one
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @return a Pointer to the new element
     */
    Condition::Pointer Clone(IndexType NewId, NodesArrayType const& rThisNodes) const override
    {
        Condition::Pointer pNewCondition = Create(NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

        pNewCondition->SetData(this->GetData());
        pNewCondition->SetFlags(this->GetFlags());

        return pNewCondition;
    }


    /// Calculates the LHS and RHS condition contributions
    /**
     * Clones the selected element variables, creating a new one
     * @param rLeftHandSideMatrix reference to the LHS matrix
     * @param rRightHandSideVector reference to the RHS matrix
     * @param rCurrentProcessInfo reference to the ProcessInfo (unused)
     */
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      const ProcessInfo& rCurrentProcessInfo) override;


    /// Calculates the RHS condition contributions
    /**
     * Clones the selected element variables, creating a new one
     * @param rLeftHandSideMatrix reference to the LHS matrix
     * @param rCurrentProcessInfo reference to the ProcessInfo (unused)
     */
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                       const ProcessInfo& rCurrentProcessInfo) override;


    /// Calculates the RHS condition contributions
    /**
     * Clones the selected element variables, creating a new one
     * @param rRightHandSideVector reference to the RHS matrix
     * @param rCurrentProcessInfo reference to the ProcessInfo (unused)
     */
    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                        const ProcessInfo& rCurrentProcessInfo) override;



    /// Condition check
    /**
     * @param rCurrentProcessInfo reference to the ProcessInfo
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) const override;


    /// Provides the global indices for each one of this element's local rows.
    /** This determines the elemental equation ID vector for all elemental DOFs
     * @param rResult A vector containing the global Id of each row
     * @param rCurrentProcessInfo the current process info object (unused)
     */
    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const override;

    /// Returns a list of the element's Dofs
    /**
     * @param ElementalDofList the list of DOFs
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetDofList(DofsVectorType& ConditionDofList, const ProcessInfo& CurrentProcessInfo) const override;

    void Calculate(
        const Variable< array_1d<double,3> >& rVariable,
        array_1d<double,3>& Output,
        const ProcessInfo& rCurrentProcessInfo) override;

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
        buffer << "NavierStokesWallCondition" << TDim << "D";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "NavierStokesWallCondition";
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

    void CalculateNormal(array_1d<double,3>& An);

    void ComputeGaussPointLHSContribution(BoundedMatrix<double,TNumNodes*(TDim+1),TNumNodes*(TDim+1)>& lhs, const ConditionDataStruct& data);

    void ComputeGaussPointRHSContribution(array_1d<double,TNumNodes*(TDim+1)>& rhs, const ConditionDataStruct& data);

    void ComputeRHSNeumannContribution(array_1d<double,TNumNodes*(TDim+1)>& rhs, const ConditionDataStruct& data);

    void ComputeRHSOutletInflowContribution(array_1d<double,TNumNodes*(TDim+1)>& rhs, const ConditionDataStruct& data);

    /**
     * @brief Computes the right-hand side of the Navier slip contribution as e.g. described in BEHR2004
     * The (Navier) slip length is read as a nodal variable.
     * If a smaller value is set, tangential velocities lead to a higher tangential traction.
     * Though only tangential velocities should appear, a tangetial projection is added.
     * (Reference BEHR2004: https://onlinelibrary.wiley.com/doi/abs/10.1002/fld.663)
     * @param rRightHandSideVector reference to the RHS vector
     * @param rDataStruct reference to a struct to hand over data
     */
    virtual void ComputeGaussPointNavierSlipRHSContribution(
        array_1d<double, TNumNodes*(TDim+1)>& rRightHandSideVector,
        const ConditionDataStruct& rDataStruct );

    /**
     * @brief Computes the left-hand side of the Navier slip contribution as e.g. described in BEHR2004
     * The (Navier) slip length is read as a nodal variable.
     * If a smaller value is set, tangential velocities lead to a higher tangential traction.
     * Though only tangential velocities should appear, a tangetial projection is added.
     * (Reference BEHR2004: https://onlinelibrary.wiley.com/doi/abs/10.1002/fld.663)
     * @param rLeftHandSideMatrix reference to the LHS matrix
     * @param rDataStruct reference to a struct to hand over data
     */
    virtual void ComputeGaussPointNavierSlipLHSContribution(
        BoundedMatrix<double, TNumNodes*(TDim+1),TNumNodes*(TDim+1)>& rLeftHandSideMatrix,
        const ConditionDataStruct& rDataStruct);

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

    /**
     * @brief Computes the left-hand side contribution for the BEHR2004 slip condition
     * This specific implementation of the slip condition avoids spurious velocities
     * at points were the normal directions of the adjacent boundary geometries do not
     * coincide (Reference BEHR2004: https://onlinelibrary.wiley.com/doi/abs/10.1002/fld.663)
     * @param rLeftHandSideMatrix reference to the LHS matrix
     * @param rDataStruct reference to a struct to hand over data
     */
    void ComputeGaussPointBehrSlipLHSContribution(  BoundedMatrix<double,TNumNodes*(TDim+1),TNumNodes*(TDim+1)>& rLeftHandSideMatrix,
                                                    const ConditionDataStruct& rDataStruct );


    /**
     * @brief Computes the right-hand side contribution for the BEHR2004 slip condition
     * This specific implementation of the slip condition avoids spurious velocities
     * at points were the normal directions of the adjacent boundary geometries do not
     * coincide (Reference BEHR2004: https://onlinelibrary.wiley.com/doi/abs/10.1002/fld.663)
     * @param rRightHandSideVector reference to the RHS vector
     * @param rDataStruct reference to a struct to hand over data
     */
    void ComputeGaussPointBehrSlipRHSContribution(  array_1d<double,TNumNodes*(TDim+1)>& rRightHandSideVector,
                                                    const ConditionDataStruct& rDataStruct );

    /**
     * @brief Project the viscous stress
     * Provided a viscous stress tensor (in Voigt notation) and a unit normal vector,
     * this method calculates and returns the projection of the shear stress onto the normal
     * @param rViscousStress The viscous stress in Voigt notation
     * @param rNormal The unit normal vector to project onto
     * @param rProjectedViscousStress The projected viscous stress
     */
    void ProjectViscousStress(
        const Vector& rViscousStress,
        const array_1d<double,3> rNormal,
        array_1d<double,3> rProjectedViscousStress);

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

}; // Class NavierStokesWallCondition


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< unsigned int TDim, unsigned int TNumNodes >
inline std::istream& operator >> (std::istream& rIStream, NavierStokesWallCondition<TDim,TNumNodes>& rThis)
{
    return rIStream;
}

/// output stream function
template< unsigned int TDim, unsigned int TNumNodes >
inline std::ostream& operator << (std::ostream& rOStream, const NavierStokesWallCondition<TDim,TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block


}  // namespace Kratos.

#endif // KRATOS_NAVIER_STOKES_WALL_CONDITION_H