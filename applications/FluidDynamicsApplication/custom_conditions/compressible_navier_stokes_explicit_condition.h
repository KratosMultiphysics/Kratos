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

#ifndef KRATOS_COMPRESSIBLE_NAVIER_STOKES_CONDITION_H
#define KRATOS_COMPRESSIBLE_NAVIER_STOKES_CONDITION_H

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
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) CompressibleNavierStokesExplicitCondition : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CompressibleNavierStokesExplicitCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(CompressibleNavierStokesExplicitCondition);

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

    struct GradientData
    {
        array_1d<double, 3> density;
        BoundedMatrix<double, 3, 3> momentum;
        array_1d<double, 3> total_energy;
    };

    struct ConditionDataStruct
    {
        static constexpr SizeType NGauss = TNumNodes; // TODO: Not necessarily true

        BoundedMatrix<double, TNumNodes, BlockSize> U;
        BoundedMatrix<double, TNumNodes, BlockSize> dUdt;

        std::array<GradientData, NGauss> gradients;

        array_1d<double, TNumNodes> alpha_sc_nodes;
        array_1d<double, TNumNodes> mu_sc_nodes;
        array_1d<double, TNumNodes> beta_sc_nodes;
        array_1d<double, TNumNodes> lamb_sc_nodes;

        array_1d<double, TNumNodes> N;
        array_1d<double, 3> unit_normal {3, 0.0};


        double h;           // Element size
        double volume;      // In 2D: element area. In 3D: element volume
        double c_v;         // Heat capacity at constant volume
        double gamma;       // Heat capacity ratio
        double mu;          // Dynamic viscosity
        double lambda;      // Heat conductivity
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CompressibleNavierStokesExplicitCondition(IndexType NewId = 0):
        Condition(NewId)
    {
    }

    /// Constructor using an array of nodes
    /**
     @param NewId Index of the new condition
     @param ThisNodes An array containing the nodes of the new condition
     */
    CompressibleNavierStokesExplicitCondition(IndexType NewId, const NodesArrayType& ThisNodes):
        Condition(NewId, ThisNodes)
    {
    }

    /// Constructor using Geometry
    /**
     @param NewId Index of the new condition
     @param pGeometry Pointer to a geometry object
     */
    CompressibleNavierStokesExplicitCondition(IndexType NewId, GeometryType::Pointer pGeometry):
        Condition(NewId, pGeometry)
    {
    }

    /// Constructor using Properties
    /**
     @param NewId Index of the new element
     @param pGeometry Pointer to a geometry object
     @param pProperties Pointer to the element's properties
     */
    CompressibleNavierStokesExplicitCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):
        Condition(NewId, pGeometry, pProperties)
    {
    }

    /// Copy constructor.
    CompressibleNavierStokesExplicitCondition(CompressibleNavierStokesExplicitCondition const& rOther):
        Condition(rOther)
    {
    }

    /// Destructor.
    ~CompressibleNavierStokesExplicitCondition() override {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator
    CompressibleNavierStokesExplicitCondition & operator=(CompressibleNavierStokesExplicitCondition const& rOther)
    {
        Condition::operator=(rOther);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    GeometryData::IntegrationMethod GetIntegrationMethod() const override;

    /// Create a new CompressibleNavierStokesExplicitCondition object.
    /**
      @param NewId Index of the new condition
      @param ThisNodes An array containing the nodes of the new condition
      @param pProperties Pointer to the element's properties
      */
    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<CompressibleNavierStokesExplicitCondition>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    /// Create a new CompressibleNavierStokesExplicitCondition object.
    /**
      @param NewId Index of the new condition
      @param pGeom A pointer to the condition's geometry
      @param pProperties Pointer to the element's properties
      */
    Condition::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive< CompressibleNavierStokesExplicitCondition >(NewId, pGeom, pProperties);
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
            KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].SolutionStepsDataHas(BODY_FORCE)) << "Missing BODY_FORCE variable on solution step data for node " << this->GetGeometry()[i].Id();
            KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].SolutionStepsDataHas(HEAT_SOURCE)) << "Missing HEAT_SOURCE variable on solution step data for node " << this->GetGeometry()[i].Id();

            // Activate as soon as we start using the explicit DOF based strategy
            KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].HasDofFor(DENSITY)) << "Missing DENSITY DOF in node ", this->GetGeometry()[i].Id();
            KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].HasDofFor(MOMENTUM_X) || this->GetGeometry()[i].HasDofFor(MOMENTUM_Y)) << "Missing MOMENTUM component DOF in node ", this->GetGeometry()[i].Id();
            if (TDim == 3) {
                KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].HasDofFor(MOMENTUM_Z)) << "Missing MOMENTUM component DOF in node ", this->GetGeometry()[i].Id();
            }
            KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].HasDofFor(TOTAL_ENERGY)) << "Missing TOTAL_ENERGY DOF in node ", this->GetGeometry()[i].Id();
        }

        KRATOS_ERROR_IF_NOT(this->Has(NEIGHBOUR_ELEMENTS)) << "Missing NEIGHBOUR_ELEMENTS variable in condition #" << this->Id();
        KRATOS_ERROR_IF(this->GetValue(NEIGHBOUR_ELEMENTS).size() == 0) << "Empty NEIGHBOUR_ELEMENTS container in condition #" << this->Id();
        KRATOS_ERROR_IF(this->GetValue(NEIGHBOUR_ELEMENTS).size() > 1) << "NEIGHBOUR_ELEMENTS has more than one entry in condition #" << this->Id();

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

    void Calculate(
        const Variable< array_1d<double,3> >& rVariable,
        array_1d<double,3>& Output,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        // TODO
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

    const Parameters GetSpecifications() const override
    {
        const Parameters specifications = Parameters(R"({
            "time_integration"           : ["explicit"],
            "framework"                  : "eulerian",
            "symmetric_lhs"              : false,
            "positive_definite_lhs"      : true,
            "output"                     : {
                "gauss_point"            : ["SHOCK_SENSOR","SHEAR_SENSOR","THERMAL_SENSOR","ARTIFICIAL_CONDUCTIVITY","ARTIFICIAL_BULK_VISCOSITY","VELOCITY_DIVERGENCE"],
                "nodal_historical"       : ["DENSITY","MOMENTUM","TOTAL_ENERGY"],
                "nodal_non_historical"   : ["ARTIFICIAL_MASS_DIFFUSIVITY","ARTIFICIAL_DYNAMIC_VISCOSITY","ARTIFICIAL_BULK_VISCOSITY","ARTIFICIAL_CONDUCTIVITY","DENSITY_PROJECTION","MOMENTUM_PROJECTION","TOTAL_ENERGY_PROJECTION"],
                "entity"                 : []
            },
            "required_variables"         : ["DENSITY","MOMENTUM","TOTAL_ENERGY"],
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
        buffer << "CompressibleNavierStokesExplicitCondition" << TDim << "D" << TNumNodes << "N";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "CompressibleNavierStokesExplicitCondition";
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
    static void ComputeGaussPointData(
        Element& rParentElement,
        ConditionDataStruct& rData,
        const ProcessInfo& rCurrentProcessInfo)
    {
        std::vector<array_1d<double, 3>> grad_density;
        std::vector<array_1d<double, 3>> grad_total_energy;
        std::vector<Matrix> grad_momentum;

        rParentElement.CalculateOnIntegrationPoints(DENSITY_GRADIENT, grad_density, rCurrentProcessInfo);
        rParentElement.CalculateOnIntegrationPoints(MOMENTUM_GRADIENT, grad_momentum, rCurrentProcessInfo);
        rParentElement.CalculateOnIntegrationPoints(TOTAL_ENERGY_GRADIENT, grad_total_energy, rCurrentProcessInfo);

        for(IndexType g=0; g < ConditionDataStruct::NGauss; ++g)
        {
            rData.gradients[g].density = grad_density[0];
            rData.gradients[g].momentum = grad_momentum[0];
            rData.gradients[g].total_energy = grad_total_energy[0];
        }
    }

    /**
     * @brief Fill condition data
     * Auxiliary function to fill the element data structure
     * @param rData Reference to the element data structure to be filled
     * @param rCurrentProcessInfo Reference to the current process info
     */
    ConditionDataStruct ConditionData(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        ConditionDataStruct data;

        // Element adjacent to this condition
        auto& r_parent_element = *GetValue(NEIGHBOUR_ELEMENTS).begin();

        // Getting data for the given geometry
        const auto& r_geometry = GetGeometry();

        // Loads shape function info only if jacobian is uniform
        ComputeGeometryData(r_geometry, data);
        ComputeGaussPointData(r_parent_element, data, rCurrentProcessInfo);

        // Database access to all of the variables needed
        const Properties &r_properties = r_parent_element.GetProperties();
        data.mu = r_properties.GetValue(DYNAMIC_VISCOSITY);
        data.lambda = r_properties.GetValue(CONDUCTIVITY);
        data.c_v = r_properties.GetValue(SPECIFIC_HEAT); // TODO: WE SHOULD SPECIFY WHICH ONE --> CREATE SPECIFIC_HEAT_CONSTANT_VOLUME
        data.gamma = r_properties.GetValue(HEAT_CAPACITY_RATIO);

        // Magnitudes to calculate the time derivatives
        const double time_step = rCurrentProcessInfo[DELTA_TIME];
        const double theta = rCurrentProcessInfo[TIME_INTEGRATION_THETA];
        const double aux_theta = theta > 0 ? 1.0 / (theta * time_step) : 0.0;

        // Nodal data
        for(SizeType i=0; i<NumNodes; ++i)
        {
            // Shock capturing data
            const auto& r_node = r_geometry[i];
            data.alpha_sc_nodes(i) = r_node.GetValue(ARTIFICIAL_MASS_DIFFUSIVITY);
            data.mu_sc_nodes(i)    = r_node.GetValue(ARTIFICIAL_DYNAMIC_VISCOSITY);
            data.beta_sc_nodes(i)  = r_node.GetValue(ARTIFICIAL_BULK_VISCOSITY);
            data.lamb_sc_nodes(i)  = r_node.GetValue(ARTIFICIAL_CONDUCTIVITY);

            // Momentum data
            const array_1d<double, 3>& r_momentum = r_node.FastGetSolutionStepValue(MOMENTUM);
            const array_1d<double, 3>& r_momentum_old = r_node.FastGetSolutionStepValue(MOMENTUM, 1);
            const array_1d<double, 3> mom_inc = r_momentum - r_momentum_old;

            for (unsigned int d = 0; d < Dim; ++d) {
                data.U(i, 1+d)    = r_momentum[d];
                data.dUdt(i, 1+d) = aux_theta * mom_inc[d];
            }

            // Density data
            const double& r_rho = r_node.FastGetSolutionStepValue(DENSITY);
            const double& r_rho_old = r_node.FastGetSolutionStepValue(DENSITY, 1);
            data.U(i, 0)    = r_rho;
            data.dUdt(i, 0) = aux_theta * (r_rho - r_rho_old);

            // Total energy data
            const double& r_tot_ener = r_node.FastGetSolutionStepValue(TOTAL_ENERGY);
            const double& r_tot_ener_old = r_node.FastGetSolutionStepValue(TOTAL_ENERGY, 1);
            data.U(i, TDim+1)    = r_tot_ener;
            data.dUdt(i, TDim+1) = aux_theta * (r_tot_ener - r_tot_ener_old);
        }

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

}; // Class CompressibleNavierStokesExplicitCondition


///@}

///@name Type Definitions
///@{


/**
 *  This namespace is used to implement the compile-time choice of:
 *  1. Shape function calculator
 *  2. Element size calculator
 *
 *  [1] This is necessary because simplex geometries' shape functions are known
 *  at compile-time whereas non-simplex's ones are not (because they depend on
 *  the actual shape of the element).
 *
 *  [2] This is useful because GradientsElementSize is only implemented for
 *  simplex types.
 *
 *  This problem is solved using SFINAE
 */
namespace CompressibleNavierStokesExplicitConditionInternal
{
    template<unsigned int TDim, unsigned int TNumNodes>
    using ConditionDataStruct = typename CompressibleNavierStokesExplicitCondition<TDim, TNumNodes>::ConditionDataStruct;

    constexpr bool IsSimplex(const unsigned int Dimensions, const unsigned int NNodes)
    {
        return Dimensions == NNodes;
    }

    // Specialization for simplex geometries
    template<unsigned int TDim, unsigned int TNumNodes>
    inline typename std::enable_if<IsSimplex(TDim, TNumNodes), void>::type ComputeGeometryDataImpl(
        const Geometry<Node<3>> & rGeometry,
        ConditionDataStruct<TDim, TNumNodes>& rData)
    {
        BoundedMatrix<double, TNumNodes, 1> DN_DX; // unused
        GeometryUtils::CalculateGeometryData(rGeometry, DN_DX, rData.N, rData.volume);
        rData.unit_normal = rGeometry.UnitNormal(0);
    }

    /**
     * Specialization for geometries with non-uniform jacobian (non-simplex)
     * Shape functions cannot be obtained here. They will need to be calculated
     * during integration at each gauss point.
     */
    template<unsigned int TDim, unsigned int TNumNodes>
    inline typename std::enable_if<!IsSimplex(TDim, TNumNodes), void>::type ComputeGeometryDataImpl(
        const Geometry<Node<3>> & rGeometry,
        ConditionDataStruct<TDim, TNumNodes>& rData)
    {
        rData.volume = rGeometry.DomainSize();
        // Normal must be computed at each Gauss point
    }
}

template<unsigned int TDim, unsigned int TNumNodes>
void CompressibleNavierStokesExplicitCondition<TDim, TNumNodes>::ComputeGeometryData(
    const GeometryType& rGeometry,
    ConditionDataStruct& rData)
{
    CompressibleNavierStokesExplicitConditionInternal::ComputeGeometryDataImpl<TDim, TNumNodes>(rGeometry, rData);
}


///@}
///@name Input and output
///@{


/// input stream function
template< unsigned int TDim, unsigned int TNumNodes >
inline std::istream& operator >> (std::istream& rIStream, CompressibleNavierStokesExplicitCondition<TDim,TNumNodes>& rThis)
{
    return rIStream;
}

/// output stream function
template< unsigned int TDim, unsigned int TNumNodes >
inline std::ostream& operator << (std::ostream& rOStream, const CompressibleNavierStokesExplicitCondition<TDim,TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << '\n';
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block


}  // namespace Kratos.

#endif // KRATOS_COMPRESSIBLE_NAVIER_STOKES_CONDITION_H