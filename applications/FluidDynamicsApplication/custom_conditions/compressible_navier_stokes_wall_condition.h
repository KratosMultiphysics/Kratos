//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#ifndef KRATOS_COMPRESSIBLE_NAVIER_STOKES_WALL_CONDITION_H
#define KRATOS_COMPRESSIBLE_NAVIER_STOKES_WALL_CONDITION_H

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

/**
 * @brief Wall condition for compressible Navier-Stokes formulations
 * This class implements a wall condition for the compressible Navier-Stokes formulations
 * written in terms of the conservative variables (density, momentum and total energy)
 * In order to ensure the flux conservation, this condition integrates the fluxes over the domain boundaries
 * The fluxes are computed from the current solution values or, if a Neumann BC is imposed, they are
 * substituted by the corresponding flux on that part of the boundary.
 * @tparam TDim Dimension of the problem
 * @tparam TNumNodes Condition number of nodes
 */
template< unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize = TDim + 2>
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) CompressibleNavierStokesWallCondition : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CompressibleNavierStokesWallCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(CompressibleNavierStokesWallCondition);

    /**
     * @brief Structure to pass around the data
     * Auxiliary condition collecting the required condition data
     */
    struct ConditionDataStruct
    {
        // Geometrical properties
        double Area;
        array_1d<double, 3> Normal;
        array_1d<double, TNumNodes > N;
        BoundedMatrix<double, TNumNodes, TDim-1> DN_DX;

        // Physical properties
        double mu;
        double c_v;
        double gamma;
        double lambda;

        // Solution magnitudes
        // array_1d<double, TNumNodes> nu_sc_node;
        // array_1d<double, TNumNodes> alpha_sc_node;
        BoundedMatrix<double, TNumNodes, TBlockSize> U;
        
        // Other
        // bool ShockCapturing;
        array_1d<double, TNumNodes> beta_sc_nodes;
        array_1d<double, TNumNodes> lamb_sc_nodes;
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
      @param NewId Index for the new condition
      */
    CompressibleNavierStokesWallCondition(IndexType NewId = 0)
        : Condition(NewId)
    {
    }

    /// Constructor using an array of nodes
    /**
     @param NewId Index of the new condition
     @param ThisNodes An array containing the nodes of the new condition
     */
    CompressibleNavierStokesWallCondition(
        IndexType NewId,
        const NodesArrayType& ThisNodes)
        : Condition(NewId, ThisNodes)
    {
    }

    /// Constructor using Geometry
    /**
     @param NewId Index of the new condition
     @param pGeometry Pointer to a geometry object
     */
    CompressibleNavierStokesWallCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : Condition(NewId, pGeometry)
    {
    }

    /// Constructor using Properties
    /**
     @param NewId Index of the new element
     @param pGeometry Pointer to a geometry object
     @param pProperties Pointer to the element's properties
     */
    CompressibleNavierStokesWallCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : Condition(NewId, pGeometry, pProperties)
    {
    }

    /// Copy constructor.
    CompressibleNavierStokesWallCondition(CompressibleNavierStokesWallCondition const& rOther)
        : Condition(rOther)
    {
    }

    /// Destructor.
    ~CompressibleNavierStokesWallCondition() override
    {
    }


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator
    CompressibleNavierStokesWallCondition & operator=(CompressibleNavierStokesWallCondition const& rOther)
    {
        Condition::operator=(rOther);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /// Create a new CompressibleNavierStokesWallCondition object.
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
        return Kratos::make_intrusive<CompressibleNavierStokesWallCondition>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    /// Create a new CompressibleNavierStokesWallCondition object.
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
        return Kratos::make_intrusive<CompressibleNavierStokesWallCondition>(NewId, pGeom, pProperties);
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
        Condition::Pointer pNewCondition = Create(NewId, GetGeometry().Create(rThisNodes), pGetProperties());

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
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    /// Calculates the RHS condition contributions
    /**
     * Clones the selected element variables, creating a new one
     * @param rLeftHandSideMatrix reference to the LHS matrix
     * @param rCurrentProcessInfo reference to the ProcessInfo (unused)
     */
    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;


    /// Calculates the RHS condition contributions
    /**
     * Clones the selected element variables, creating a new one
     * @param rRightHandSideVector reference to the RHS matrix
     * @param rCurrentProcessInfo reference to the ProcessInfo (unused)
     */
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    /// Adds the condition explicit contribution
    /**
     * This is called during the assembling process in order
     * to calculate the elemental contribution in explicit calculation.
     * NodalData is modified Inside the function, so the
     * The "AddEXplicit" FUNCTIONS THE ONLY FUNCTIONS IN WHICH AN ELEMENT
     * IS ALLOWED TO WRITE ON ITS NODES.
     * the caller is expected to ensure thread safety hence
     * SET/UNSETLOCK MUST BE PERFORMED IN THE STRATEGY BEFORE CALLING THIS FUNCTION
      * @param rCurrentProcessInfo the current process info instance
     */
    void AddExplicitContribution(const ProcessInfo &rCurrentProcessInfo) override;

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
        buffer << "CompressibleNavierStokesWallCondition" << TDim << "D";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "CompressibleNavierStokesWallCondition";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {

    }

    ///@}
    ///@name Friends
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

    void FillConditionData(
        ConditionDataStruct& rData,
        const ProcessInfo& rCurrentProcessInfo);

    /**
     * @brief Internal CalculateRightHandSide() method
     * This auxiliary RHS calculated method is created to bypass the element API
     * In this way bounded vectors can be used in the explicit residual calculation
     * @param rRightHandSideBoundedVector Reference to the auxiliary RHS vector
     * @param rCurrentProcessInfo Reference to the current process info
     */
    void CalculateRightHandSideInternal(
        BoundedVector<double, TBlockSize * TNumNodes>& rRightHandSideBoundedVector,
        const ProcessInfo& rCurrentProcessInfo);

    void CalculateNormal(array_1d<double,3>& rNormal) const;

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
}; // Class CompressibleNavierStokesWallCondition

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
template< unsigned int TDim, unsigned int TNumNodes >
inline std::istream& operator >> (std::istream& rIStream, CompressibleNavierStokesWallCondition<TDim,TNumNodes>& rThis)
{
    return rIStream;
}

/// output stream function
template< unsigned int TDim, unsigned int TNumNodes >
inline std::ostream& operator << (std::ostream& rOStream, const CompressibleNavierStokesWallCondition<TDim,TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}
///@} addtogroup block
}  // namespace Kratos.

#endif // KRATOS_COMPRESSIBLE_NAVIER_STOKES_WALL_CONDITION_H