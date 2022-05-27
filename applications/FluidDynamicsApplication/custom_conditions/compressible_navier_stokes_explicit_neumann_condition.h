//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Eduard Gómez
//

#ifndef KRATOS_COMPRESSIBLE_NAVIER_STOKES_NEUMANN_CONDITION_H
#define KRATOS_COMPRESSIBLE_NAVIER_STOKES_NEUMANN_CONDITION_H

// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "geometries/geometry.h"
#include "includes/define.h"
#include "includes/condition.h"
#include "includes/model_part.h"
#include "includes/serializer.h"
#include "includes/process_info.h"
#include "includes/variables.h"
#include "utilities/element_size_calculator.h"
#include "utilities/atomic_utilities.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "includes/deprecated_variables.h"

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

template<unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) CompressibleNavierStokesExplicitNeumannCondition : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CompressibleNavierStokesExplicitNeumannCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(CompressibleNavierStokesExplicitNeumannCondition);

    using NodeType = Condition::NodeType;
    using PropertiesType = Condition::PropertiesType;
    using GeometryType = Condition::GeometryType;
    using NodesArrayType = Condition::NodesArrayType;
    using VectorType = Condition::VectorType;
    using MatrixType = Condition::MatrixType;
    using IndexType = Condition::IndexType;
    using EquationIdVectorType = Condition::EquationIdVectorType;
    using DofsVectorType = Condition::DofsVectorType;

    constexpr static unsigned int Dim = TDim;
    constexpr static unsigned int NumNodes = TNumNodes;
    constexpr static unsigned int BlockSize = Dim + 2;
    constexpr static unsigned int DofSize = NumNodes * BlockSize;

    struct FluxData
    {
        double density;
        array_1d<double, 3> momentum;
        double total_energy;
    };

    struct ConditionDataStruct
    {
        static constexpr SizeType NGauss = TNumNodes; // TODO: Not necessarily true
        std::array<FluxData, NGauss> fluxes;
        double volume;
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CompressibleNavierStokesExplicitNeumannCondition(IndexType NewId = 0):
        Condition(NewId)
    {
    }

    /// Constructor using an array of nodes
    /**
     @param NewId Index of the new condition
     @param ThisNodes An array containing the nodes of the new condition
     */
    CompressibleNavierStokesExplicitNeumannCondition(IndexType NewId, const NodesArrayType& ThisNodes):
        Condition(NewId, ThisNodes)
    {
    }

    /// Constructor using Geometry
    /**
     @param NewId Index of the new condition
     @param pGeometry Pointer to a geometry object
     */
    CompressibleNavierStokesExplicitNeumannCondition(IndexType NewId, GeometryType::Pointer pGeometry):
        Condition(NewId, pGeometry)
    {
    }

    /// Constructor using Properties
    /**
     @param NewId Index of the new element
     @param pGeometry Pointer to a geometry object
     @param pProperties Pointer to the element's properties
     */
    CompressibleNavierStokesExplicitNeumannCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):
        Condition(NewId, pGeometry, pProperties)
    {
    }

    /// Copy constructor.
    CompressibleNavierStokesExplicitNeumannCondition(CompressibleNavierStokesExplicitNeumannCondition const& rOther):
        Condition(rOther)
    {
    }

    /// Destructor.
    ~CompressibleNavierStokesExplicitNeumannCondition() override {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator
    CompressibleNavierStokesExplicitNeumannCondition & operator=(CompressibleNavierStokesExplicitNeumannCondition const& rOther)
    {
        Condition::operator=(rOther);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    GeometryData::IntegrationMethod GetIntegrationMethod() const override;

    /**
     * is called to initialize the condition
     * if the condition needs to perform any operation before any calculation is done
     * the condition variables will be initialized and set using this method
     */
    void Initialize(const ProcessInfo& rCurrentProcessInfo) override
    {
        // Fluxes are set to zero by default
        if(!this->Has(DENSITY_FLUX))         this->SetValue(DENSITY_FLUX, DENSITY_FLUX.Zero());
        if(!this->Has(MOMENTUM_FLUX))        this->SetValue(MOMENTUM_FLUX, MOMENTUM_FLUX.Zero());
        if(!this->Has(TOTAL_ENERGY_FLUX))    this->SetValue(TOTAL_ENERGY_FLUX, TOTAL_ENERGY_FLUX.Zero());
    }

    /// Create a new CompressibleNavierStokesExplicitNeumannCondition object.
    /**
      @param NewId Index of the new condition
      @param ThisNodes An array containing the nodes of the new condition
      @param pProperties Pointer to the element's properties
      */
    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<CompressibleNavierStokesExplicitNeumannCondition>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    /// Create a new CompressibleNavierStokesExplicitNeumannCondition object.
    /**
      @param NewId Index of the new condition
      @param pGeom A pointer to the condition's geometry
      @param pProperties Pointer to the element's properties
      */
    Condition::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive< CompressibleNavierStokesExplicitNeumannCondition >(NewId, pGeom, pProperties);
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
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        KRATOS_ERROR << "Calling the CalculateLocalSystem() method for the explicit compressible Navier-Stokes condition.";

        KRATOS_CATCH("")
    }


    /// Calculates the RHS condition contributions
    /**
     * Clones the selected element variables, creating a new one
     * @param rLeftHandSideMatrix reference to the LHS matrix
     * @param rCurrentProcessInfo reference to the ProcessInfo (unused)
     */
    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        KRATOS_ERROR << "Calling the CalculateLeftHandSide() method for the explicit compressible Navier-Stokes condition.";

        KRATOS_CATCH("")
    }


    /// Calculates the RHS condition contributions
    /**
     * Clones the selected element variables, creating a new one
     * @param rRightHandSideVector reference to the RHS matrix
     * @param rCurrentProcessInfo reference to the ProcessInfo (unused)
     */
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        KRATOS_ERROR << "Calling the CalculateRightHandSide() method for the explicit compressible Navier-Stokes condition.";

        KRATOS_CATCH("")
    }

    /**
     * This is called during the assembling process in order
     * to calculate the elemental contribution in explicit calculation.
     * NodalData is modified Inside the function, so the
     * The "AddExplicit" FUNCTIONS THE ONLY FUNCTIONS IN WHICH AN ELEMENT
     * IS ALLOWED TO WRITE ON ITS NODES.
     * the caller is expected to ensure thread safety hence
     * SET/UNSETLOCK MUST BE PERFORMED IN THE STRATEGY BEFORE CALLING THIS FUNCTION
      * @param rCurrentProcessInfo the current process info instance
     */
    void AddExplicitContribution(const ProcessInfo &rCurrentProcessInfo) override
    {
        // Calculate the explicit residual vector
        const auto rhs = CalculateRightHandSideInternal(rCurrentProcessInfo);

        // Add the residual contribution
        // Note that the reaction is indeed the formulation residual
        auto& r_geometry = GetGeometry();

        for (IndexType i_node = 0; i_node < NumNodes; ++i_node)
        {
            auto& r_node = r_geometry[i_node];

            IndexType i_dof = BlockSize * i_node;

            AtomicAdd(r_node.FastGetSolutionStepValue(REACTION_DENSITY), rhs[i_dof++]);

            for (IndexType d = 0; d < Dim; ++d)
            {
                AtomicAdd(r_node.FastGetSolutionStepValue(REACTION)[d], rhs[i_dof++]);
            }

            AtomicAdd(r_node.FastGetSolutionStepValue(REACTION_ENERGY),  rhs[i_dof]);
        }
    }



    /// Condition check
    /**
     * @param rCurrentProcessInfo reference to the ProcessInfo
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) const override
    {
        KRATOS_TRY

        // Perform basic element checks
        int error_code = Kratos::Condition::Check(rCurrentProcessInfo);
        if (error_code != 0) {
            return error_code;
        }

        // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].SolutionStepsDataHas(DENSITY)) << "Missing DENSITY variable on solution step data for node " << this->GetGeometry()[i].Id();
            KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].SolutionStepsDataHas(MOMENTUM)) << "Missing MOMENTUM variable on solution step data for node " << this->GetGeometry()[i].Id();
            KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].SolutionStepsDataHas(TOTAL_ENERGY)) << "Missing TOTAL_ENERGY variable on solution step data for node " << this->GetGeometry()[i].Id();

            // Activate as soon as we start using the explicit DOF based strategy
            KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].HasDofFor(DENSITY)) << "Missing DENSITY DOF in node ", this->GetGeometry()[i].Id();
            KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].HasDofFor(MOMENTUM_X) || this->GetGeometry()[i].HasDofFor(MOMENTUM_Y)) << "Missing MOMENTUM component DOF in node ", this->GetGeometry()[i].Id();
            if (TDim == 3) {
                KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].HasDofFor(MOMENTUM_Z)) << "Missing MOMENTUM component DOF in node ", this->GetGeometry()[i].Id();
            }
            KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].HasDofFor(TOTAL_ENERGY)) << "Missing TOTAL_ENERGY DOF in node ", this->GetGeometry()[i].Id();
        }

        return 0;

        KRATOS_CATCH("");
    }

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


    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    const Parameters GetSpecifications() const override
    {
        const Parameters specifications = Parameters(R"({
            "time_integration"           : ["explicit"],
            "framework"                  : "eulerian",
            "symmetric_lhs"              : false,
            "positive_definite_lhs"      : true,
            "output"                     : {
                "gauss_point"            : [],
                "nodal_historical"       : [],
                "nodal_non_historical"   : [],
                "entity"                 : []
            },
            "required_variables"         : ["DENSITY_FLUX","MOMENTUM_FLUX","TOTAL_ENERGY_FLUX"],
            "required_dofs"              : [],
            "flags_used"                 : [],
            "compatible_geometries"      : ["Line2D2"],
            "element_integrates_in_time" : true,
            "compatible_constitutive_laws": {
                "type"        : [],
                "dimension"   : [],
                "strain_size" : []
            },
            "required_polynomial_degree_of_geometry" : 1,
            "documentation"   :
                "This condition implements a compressible Navier-Stokes formulation written in conservative variables. This condition is compatible with both entropy-based and physics-based shock capturing techniques."
        })");

        std::vector<std::string> dofs(BlockSize, "");
        dofs[0] = "DENSITY";
        dofs[1] = "MOMENTUM_X";
        if(TDim >= 2) dofs[2] = "MOMENTUM_Y";
        if(TDim >= 3) dofs[3] = "MOMENTUM_Z";
        dofs[TDim+1] = "TOTAL_ENERGY";

        specifications["required_dofs"].SetStringArray(dofs);

        return specifications;
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "CompressibleNavierStokesExplicitNeumannCondition" << TDim << "D" << TNumNodes << "N";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "CompressibleNavierStokesExplicitNeumannCondition";
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

    static void ComputeGeometryData(
        const GeometryType& r_geometry,
        ConditionDataStruct& rData);


    /// TODO: This only works for simplex elements
    void ComputeGaussPointData(ConditionDataStruct& rData)
    {
        KRATOS_TRY

        const auto& flux_density = this->GetValue(DENSITY_FLUX);
        const auto& flux_momentum = this->GetValue(MOMENTUM_FLUX);
        const auto& flux_total_energy = this->GetValue(TOTAL_ENERGY_FLUX);

        for(IndexType g=0; g < ConditionDataStruct::NGauss; ++g)
        {
            // Fluxes
            rData.fluxes[g].density = flux_density;
            rData.fluxes[g].momentum = flux_momentum;
            rData.fluxes[g].total_energy = flux_total_energy;
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief Fill condition data
     * Auxiliary function to fill the element data structure
     * @param rData Reference to the element data structure to be filled
     * @param rCurrentProcessInfo Reference to the current process info
     */
    ConditionDataStruct ConditionData()
    {
        KRATOS_TRY

        ConditionDataStruct data;

        data.volume = GetGeometry().DomainSize();
        ComputeGaussPointData(data);

        return data;

        KRATOS_CATCH("")
    }

    /**
     * @brief Internal CalculateRightHandSide() method
     * This auxiliary RHS calculated method is created to bypass the element API
     * In this way bounded vectors can be used in the explicit residual calculation
     * @param rRightHandSideBoundedVector Reference to the auxiliary RHS vector
     * @param rCurrentProcessInfo Reference to the current process info
     */
    BoundedVector<double, BlockSize * TNumNodes> CalculateRightHandSideInternal(
        const ProcessInfo& rCurrentProcessInfo);


    ///@}
    ///@name Protected Operations
    ///@{


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

}; // Class CompressibleNavierStokesExplicitNeumannCondition


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< unsigned int TDim, unsigned int TNumNodes >
inline std::istream& operator >> (std::istream& rIStream, CompressibleNavierStokesExplicitNeumannCondition<TDim,TNumNodes>& rThis)
{
    return rIStream;
}

/// output stream function
template< unsigned int TDim, unsigned int TNumNodes >
inline std::ostream& operator << (std::ostream& rOStream, const CompressibleNavierStokesExplicitNeumannCondition<TDim,TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << '\n';
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block


}  // namespace Kratos.

#endif // KRATOS_COMPRESSIBLE_NAVIER_STOKES_NEUMANN_CONDITION_H