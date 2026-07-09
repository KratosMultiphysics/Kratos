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


#ifndef KRATOS_FIRST_ORDER_STOKES_VARIABLE_VISCOSITY_CONDITION_H
#define KRATOS_FIRST_ORDER_STOKES_VARIABLE_VISCOSITY_CONDITION_H

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

/// Implements a wall condition for the monolithic stokes formulation - based on BDF2.
/**
  It is intended to be used with the Stokes2D and Stokes3D elements and a BDF2 time integration scheme
  @see Stokes3D
 */
template< unsigned int TDim, unsigned int TNumNodes = TDim >
class FirstOrderStokesVariableViscosityCondition : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FirstOrderStokesVariableViscosityCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(FirstOrderStokesVariableViscosityCondition);

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
    FirstOrderStokesVariableViscosityCondition(IndexType NewId = 0):
        Condition(NewId)
    {
    }

    /// Constructor using an array of nodes
    /**
     @param NewId Index of the new condition
     @param ThisNodes An array containing the nodes of the new condition
     */
    FirstOrderStokesVariableViscosityCondition(IndexType NewId,
                            const NodesArrayType& ThisNodes):
        Condition(NewId,ThisNodes)
    {
    }

    /// Constructor using Geometry
    /**
     @param NewId Index of the new condition
     @param pGeometry Pointer to a geometry object
     */
    FirstOrderStokesVariableViscosityCondition(IndexType NewId,
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
    FirstOrderStokesVariableViscosityCondition(IndexType NewId,
                            GeometryType::Pointer pGeometry,
                            PropertiesType::Pointer pProperties):
        Condition(NewId,pGeometry,pProperties)
    {
    }

    /// Copy constructor.
    FirstOrderStokesVariableViscosityCondition(FirstOrderStokesVariableViscosityCondition const& rOther):
        Condition(rOther)
    {
    }

    /// Destructor.
    ~FirstOrderStokesVariableViscosityCondition() override {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator
    FirstOrderStokesVariableViscosityCondition & operator=(FirstOrderStokesVariableViscosityCondition const& rOther)
    {
        Condition::operator=(rOther);

        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /// Create a new FirstOrderStokesVariableViscosityCondition object.
    /**
      @param NewId Index of the new condition
      @param ThisNodes An array containing the nodes of the new condition
      @param pProperties Pointer to the element's properties
      */
    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<FirstOrderStokesVariableViscosityCondition>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    }


    /// Create a new FirstOrderStokesVariableViscosityCondition object.
    /**
      @param NewId Index of the new condition
      @param pGeom A pointer to the condition's geometry
      @param pProperties Pointer to the element's properties
      */
    Condition::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override {
        return Kratos::make_intrusive< FirstOrderStokesVariableViscosityCondition >(NewId, pGeom, pProperties);
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
        const SizeType BlockSize = TDim + 1;
        const SizeType LocalSize = BlockSize * TNumNodes;

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
        const SizeType BlockSize = TDim + 1;
        const SizeType LocalSize = BlockSize * TNumNodes;

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
        const SizeType BlockSize = TDim + 1;
        const SizeType LocalSize = BlockSize * TNumNodes;

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
                    if(this->GetGeometry()[i].SolutionStepsDataHas(PRESSURE) == false)
                        KRATOS_THROW_ERROR(std::invalid_argument,"missing PRESSURE variable on solution step data for node ",this->GetGeometry()[i].Id());
                    if(this->GetGeometry()[i].SolutionStepsDataHas(DYNAMIC_VISCOSITY) == false)
                        KRATOS_THROW_ERROR(std::invalid_argument,"missing DYNAMIC_VISCOSITY variable on solution step data for node ",this->GetGeometry()[i].Id());
                    if(this->GetGeometry()[i].HasDofFor(VELOCITY_X) == false ||
                            this->GetGeometry()[i].HasDofFor(VELOCITY_Y) == false ||
                            this->GetGeometry()[i].HasDofFor(VELOCITY_Z) == false)
                        KRATOS_THROW_ERROR(std::invalid_argument,"missing VELOCITY component degree of freedom on node ",this->GetGeometry()[i].Id());
                    if(this->GetGeometry()[i].HasDofFor(PRESSURE) == false)
                        KRATOS_THROW_ERROR(std::invalid_argument,"missing PRESSURE component degree of freedom on node ",this->GetGeometry()[i].Id());
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
        buffer << "FirstOrderStokesVariableViscosityCondition" << TDim << "D";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "FirstOrderStokesVariableViscosityCondition";
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

}; // Class FirstOrderStokesVariableViscosityCondition


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< unsigned int TDim, unsigned int TNumNodes >
inline std::istream& operator >> (std::istream& rIStream,
                                  FirstOrderStokesVariableViscosityCondition<TDim,TNumNodes>& rThis)
{
    return rIStream;
}

/// output stream function
template< unsigned int TDim, unsigned int TNumNodes >
inline std::ostream& operator << (std::ostream& rOStream,
                                  const FirstOrderStokesVariableViscosityCondition<TDim,TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block


}  // namespace Kratos.

#endif // KRATOS_FIRST_ORDER_STOKES_VARIABLE_VISCOSITY_CONDITION_H
