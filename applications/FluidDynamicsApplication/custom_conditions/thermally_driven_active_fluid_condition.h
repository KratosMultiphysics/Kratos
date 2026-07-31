//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolas Sibuet
//


#ifndef KRATOS_TEHRALLY_DRIVEN_ACTIVE_FLUID_CONDITION_H
#define KRATOS_TEHRALLY_DRIVEN_ACTIVE_FLUID_CONDITION_H

// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/process_info.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "includes/deprecated_variables.h"
#include "includes/cfd_variables.h"
#include "includes/kratos_flags.h"

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

/// Implements an empty condition for the ThermallyDrivenActiveFluidFormulation element.
/**
  It is intended to be used in combination with the ThermallyDrivenActiveFluidFormulation element.
  @see ThermallyDrivenActiveFluidFormulation
 */
template< unsigned int TDim>
class ThermallyDrivenActiveFluidCondition : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ThermallyDrivenActiveFluidCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(ThermallyDrivenActiveFluidCondition);

    static constexpr std::size_t VelocityNumNodes = TDim == 2 ? 6 : 10;

    static constexpr std::size_t PressureNumNodes = TDim == 2 ? 3 : 4;

    static constexpr std::size_t LocalSize = VelocityNumNodes*TDim*2 + PressureNumNodes*2 + VelocityNumNodes; // We have 2 velocity-like dofs (P2, vector), 2 pressure-like dofs (P1, scalar) and 1 P2 scalar dof.

    typedef Node NodeType;

    typedef Properties PropertiesType;

    typedef Geometry<NodeType> GeometryType;

    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

    typedef Vector VectorType;

    typedef Matrix MatrixType;

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    typedef std::vector<std::size_t> EquationIdVectorType;

    typedef std::vector< Dof<double>::Pointer > DofsVectorType;

    typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /** Admits an Id as a parameter.
      @param NewId Index for the new condition
      */
    ThermallyDrivenActiveFluidCondition(IndexType NewId = 0):
        Condition(NewId)
    {
    }

    /// Constructor using an array of nodes
    /**
     @param NewId Index of the new condition
     @param ThisNodes An array containing the nodes of the new condition
     */
    ThermallyDrivenActiveFluidCondition(IndexType NewId,
                            const NodesArrayType& ThisNodes):
        Condition(NewId,ThisNodes)
    {
    }

    /// Constructor using Geometry
    /**
     @param NewId Index of the new condition
     @param pGeometry Pointer to a geometry object
     */
    ThermallyDrivenActiveFluidCondition(IndexType NewId,
                            GeometryType::Pointer pGeometry):
        Condition(NewId,pGeometry)
    {
    }

    /// Constructor using Properties
    /**
     @param NewId Index of the new element
     @param pGeometry Pointer to a geometry object
     @param pProperties Pointer to the element's properties
     */
    ThermallyDrivenActiveFluidCondition(IndexType NewId,
                            GeometryType::Pointer pGeometry,
                            PropertiesType::Pointer pProperties):
        Condition(NewId,pGeometry,pProperties)
    {
    }

    /// Copy constructor.
    ThermallyDrivenActiveFluidCondition(ThermallyDrivenActiveFluidCondition const& rOther):
        Condition(rOther)
    {
    }

    /// Destructor.
    ~ThermallyDrivenActiveFluidCondition() override {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator
    ThermallyDrivenActiveFluidCondition & operator=(ThermallyDrivenActiveFluidCondition const& rOther)
    {
        Condition::operator=(rOther);

        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /// Create a new ThermallyDrivenActiveFluidCondition object.
    /**
      @param NewId Index of the new condition
      @param ThisNodes An array containing the nodes of the new condition
      @param pProperties Pointer to the element's properties
      */
    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<ThermallyDrivenActiveFluidCondition>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    }


    /// Create a new ThermallyDrivenActiveFluidCondition object.
    /**
      @param NewId Index of the new condition
      @param pGeom A pointer to the condition's geometry
      @param pProperties Pointer to the element's properties
      */
    Condition::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override {
        return Kratos::make_intrusive< ThermallyDrivenActiveFluidCondition >(NewId, pGeom, pProperties);
    }


    /// Return local contributions of the correct size, filled with zeros (for compatibility with time schemes).
    /** The actual local contributions are computed in the Damping functions
      @see CalculateLocalVelocityContribution
      */
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        if (rLeftHandSideMatrix.size1() != LocalSize)
            rLeftHandSideMatrix.resize(LocalSize,LocalSize);

        if (rRightHandSideVector.size() != LocalSize)
            rRightHandSideVector.resize(LocalSize);

        noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize,LocalSize);
        noalias(rRightHandSideVector) = ZeroVector(LocalSize);

        
        ApplyNeumannCondition(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
    }

    /// Return a matrix of the correct size, filled with zeros (for compatibility with time schemes).
    /** The actual local contributions are computed in the Damping functions
      @see DampingMatrix
      */
    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        if (rLeftHandSideMatrix.size1() != LocalSize)
            rLeftHandSideMatrix.resize(LocalSize,LocalSize);

        noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize,LocalSize);

        VectorType tmp;
        tmp.resize(LocalSize);
        ApplyNeumannCondition(rLeftHandSideMatrix,tmp,rCurrentProcessInfo);
    }

    /// Return local right hand side of the correct size, filled with zeros (for compatibility with time schemes).
    /** The actual local contributions are computed in the Damping functions
      @see CalculateLocalVelocityContribution
      */
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        if (rRightHandSideVector.size() != LocalSize)
            rRightHandSideVector.resize(LocalSize);

        noalias(rRightHandSideVector) = ZeroVector(LocalSize);

        MatrixType tmp;
        tmp.resize(LocalSize,LocalSize);
        ApplyNeumannCondition(tmp,rRightHandSideVector,rCurrentProcessInfo);
    }



    /// Check that all data required by this condition is available and reasonable
    int Check(const ProcessInfo& rCurrentProcessInfo) const override
    {
        KRATOS_TRY;

        int Check = Condition::Check(rCurrentProcessInfo); // Checks id > 0 and area > 0

        if (Check != 0)
        {
            return Check;
        }
        else
        {
                // Checks on nodes
                // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
                for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
                {
                    if(this->GetGeometry()[i].SolutionStepsDataHas(VELOCITY) == false)
                        KRATOS_THROW_ERROR(std::invalid_argument,"missing VELOCITY variable on solution step data for node ",this->GetGeometry()[i].Id());
                    if(this->GetGeometry()[i].SolutionStepsDataHas(VELOCITY_LAPLACIAN) == false)
                        KRATOS_THROW_ERROR(std::invalid_argument,"missing VELOCITY_LAPLACIAN variable on solution step data for node ",this->GetGeometry()[i].Id());
                    if(this->GetGeometry()[i].SolutionStepsDataHas(TEMPERATURE) == false)
                        KRATOS_THROW_ERROR(std::invalid_argument,"missing TEMPERATURE variable on solution step data for node ",this->GetGeometry()[i].Id());
                    if(this->GetGeometry()[i].SolutionStepsDataHas(PRESSURE) == false)
                        KRATOS_THROW_ERROR(std::invalid_argument,"missing PRESSURE variable on solution step data for node ",this->GetGeometry()[i].Id());
                    if(this->GetGeometry()[i].SolutionStepsDataHas(PRESSUREAUX) == false)
                        KRATOS_THROW_ERROR(std::invalid_argument,"missing PRESSUREAUX variable on solution step data for node ",this->GetGeometry()[i].Id());
                    if(this->GetGeometry()[i].HasDofFor(VELOCITY_X) == false ||
                            this->GetGeometry()[i].HasDofFor(VELOCITY_Y) == false ||
                            this->GetGeometry()[i].HasDofFor(VELOCITY_Z) == false)
                        KRATOS_THROW_ERROR(std::invalid_argument,"missing VELOCITY component degree of freedom on node ",this->GetGeometry()[i].Id());
                    if(this->GetGeometry()[i].HasDofFor(VELOCITY_LAPLACIAN_X) == false ||
                            this->GetGeometry()[i].HasDofFor(VELOCITY_LAPLACIAN_Y) == false ||
                            this->GetGeometry()[i].HasDofFor(VELOCITY_LAPLACIAN_Z) == false)
                        KRATOS_THROW_ERROR(std::invalid_argument,"missing VELOCITY_LAPLACIAN component degree of freedom on node ",this->GetGeometry()[i].Id());
                    if(this->GetGeometry()[i].HasDofFor(TEMPERATURE) == false)
                        KRATOS_THROW_ERROR(std::invalid_argument,"missing TEMPERATURE component degree of freedom on node ",this->GetGeometry()[i].Id());
                    if(this->GetGeometry()[i].HasDofFor(PRESSURE) == false)
                        KRATOS_THROW_ERROR(std::invalid_argument,"missing PRESSURE component degree of freedom on node ",this->GetGeometry()[i].Id());
                    if(this->GetGeometry()[i].HasDofFor(PRESSUREAUX) == false)
                            KRATOS_THROW_ERROR(std::invalid_argument,"missing PRESSUREAUX component degree of freedom on node ",this->GetGeometry()[i].Id());
                }

            return Check;
        }

        KRATOS_CATCH("");
    }


    /// Provides the global indices for each one of this element's local rows.
    /** This determines the elemental equation ID vector for all elemental DOFs
     * @param rResult A vector containing the global Id of each row
     * @param rCurrentProcessInfo the current process info object (unused)
     */
    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const override;


    /// Returns a list of the element's Dofs
    /**
     * @param ElementalDofList the list of DOFs
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetDofList(
        DofsVectorType& ConditionDofList,
        const ProcessInfo& CurrentProcessInfo) const override;


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
        buffer << "ThermallyDrivenActiveFluidCondition" << TDim << "D";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ThermallyDrivenActiveFluidCondition";
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


    void ApplyNeumannCondition(MatrixType &rLocalMatrix, VectorType &rLocalVector, const ProcessInfo &rCurrentProcessInfo);


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

}; // Class ThermallyDrivenActiveFluidCondition


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< unsigned int TDim>
inline std::istream& operator >> (std::istream& rIStream,
                                  ThermallyDrivenActiveFluidCondition<TDim>& rThis)
{
    return rIStream;
}

/// output stream function
template< unsigned int TDim>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ThermallyDrivenActiveFluidCondition<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block


}  // namespace Kratos.

#endif // KRATOS_TEHRALLY_DRIVEN_ACTIVE_FLUID_CONDITION_H
