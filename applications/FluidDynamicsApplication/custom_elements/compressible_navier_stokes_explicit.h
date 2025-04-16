
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#if !defined(KRATOS_COMPRESSIBLE_NAVIER_STOKES_EXPLICIT_H_INCLUDED)
#define KRATOS_COMPRESSIBLE_NAVIER_STOKES_EXPLICIT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/variables.h"
#include "includes/serializer.h"
#include "utilities/geometry_utilities.h"
#include "includes/cfd_variables.h"
#include "utilities/element_size_calculator.h"
#include "utilities/atomic_utilities.h"

// Application includes
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{

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
     * @brief Compressible Navier-Stokes explicit element
     * This element implements a compressible Navier-Stokes explicit formulation.
     * The formulation is written in conservative form so the element unknowns are
     * the DENSITY, MOMENTUM and TOTAL_ENERGY variables.
     * This element is intended to work with the Kratos explicit DOF based strategy.
     * Hence, the explicit residual is written in the corresponding REACTION variables.
     * @tparam TDim The space dimension (2 or 3)
     * @tparam TNumNodes The number of nodes
     */
    template <unsigned int TDim, unsigned int TNumNodes>
    class CompressibleNavierStokesExplicit : public Element
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Compile-time known quantities
        constexpr static unsigned int Dim = TDim;
        constexpr static unsigned int NumNodes = TNumNodes;
        constexpr static unsigned int BlockSize = Dim + 2;
        constexpr static unsigned int DofSize = NumNodes * BlockSize;

        /// Counted pointer of
        KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(CompressibleNavierStokesExplicit);

        struct ElementDataStruct
        {
            BoundedMatrix<double, TNumNodes, BlockSize> U;
            BoundedMatrix<double, TNumNodes, BlockSize> dUdt;
            BoundedMatrix<double, TNumNodes, BlockSize> ResProj;
            BoundedMatrix<double, TNumNodes, TDim> f_ext;
            array_1d<double, TNumNodes> m_ext;
            array_1d<double, TNumNodes> r_ext;
            array_1d<double, TNumNodes> nu_sc_node;
            array_1d<double, TNumNodes> alpha_sc_nodes;
            array_1d<double, TNumNodes> mu_sc_nodes;
            array_1d<double, TNumNodes> beta_sc_nodes;
            array_1d<double, TNumNodes> lamb_sc_nodes;

            array_1d<double, TNumNodes> N;
            BoundedMatrix<double, TNumNodes, TDim> DN_DX;

            double h;         // Element size
            double volume;    // In 2D: element area. In 3D: element volume
            double mu;        // Dynamic viscosity
            double nu;        // Kinematic viscosity
            double nu_sc;     // Kinematic viscosity (shock capturing)
            double lambda;    // Heat conductivity
            double lambda_sc; // Heat conductivity (shock capturing)
            double c_v;       // Heat capacity at constant volume
            double gamma;     // Heat capacity ratio
            double A_JWL;     // Parameter A in JWL EOS
            double B_JWL;     // Parameter B in JWL EOS
            double R1;        // Parameter R1 in JWL EOS
            double R2;        // Parameter R2 in JWL EOS
            double omega;     // Parameter omega in JWL EOS
            double rho_0;     // Parameter rho_0 in JWL EOS

            bool UseOSS;         // Use orthogonal subscales
            bool ShockCapturing; // Activate shock capturing
        };

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.
        CompressibleNavierStokesExplicit(
            IndexType NewId,
            GeometryType::Pointer pGeometry)
            : Element(NewId, pGeometry)
        {
        }

        CompressibleNavierStokesExplicit(
            IndexType NewId,
            GeometryType::Pointer pGeometry,
            PropertiesType::Pointer pProperties)
            : Element(NewId, pGeometry, pProperties)
        {
        }

        /// Destructor.
        ~CompressibleNavierStokesExplicit() override = default;

        ///@}
        ///@name Operators
        ///@{

        ///@}
        ///@name Operations
        ///@{

        Element::Pointer Create(
            IndexType NewId,
            NodesArrayType const &rThisNodes,
            PropertiesType::Pointer pProperties) const override
        {
            KRATOS_TRY
            return Kratos::make_intrusive<CompressibleNavierStokesExplicit<TDim, TNumNodes>>(NewId, this->GetGeometry().Create(rThisNodes), pProperties);
            KRATOS_CATCH("");
        }

        Element::Pointer Create(
            IndexType NewId,
            GeometryType::Pointer pGeom,
            PropertiesType::Pointer pProperties) const override
        {
            KRATOS_TRY
            return Kratos::make_intrusive<CompressibleNavierStokesExplicit<TDim, TNumNodes>>(NewId, pGeom, pProperties);
            KRATOS_CATCH("");
        }

        /**
         * This is called during the assembling process in order to
         * calculate all elemental contributions to the global system
         * matrix and the right hand side
         * Note that this is explicitly forbidden as this element is
         * conceived to only work with explicit time integration schemes
         * @param rLeftHandSideMatrix the elemental left hand side matrix
         * @param rRightHandSideVector the elemental right hand side
         * @param rCurrentProcessInfo the current process info instance
         */
        void CalculateLocalSystem(
            MatrixType &rLeftHandSideMatrix,
            VectorType &rRightHandSideVector,
            const ProcessInfo &rCurrentProcessInfo) override
        {
            KRATOS_TRY

            KRATOS_ERROR << "Calling the CalculateLocalSystem() method for the explicit compressible Navier-Stokes element.";

            KRATOS_CATCH("")
        }

        /**
         * This is called during the assembling process in order
         * to calculate the elemental right hand side vector only.
         * Note that this is explicitly forbidden as this element is
         * conceived to work with bounded arrays for the sake of efficiency.
         * A CalculateRightHandSideInternal() method is implemented instead.
         * @param rRightHandSideVector the elemental right hand side vector
         * @param rCurrentProcessInfo the current process info instance
         */
        void CalculateRightHandSide(
            VectorType &rRightHandSideVector,
            const ProcessInfo &rCurrentProcessInfo) override
        {
            KRATOS_TRY

            KRATOS_ERROR << "Calling the CalculateRightHandSide() method for the explicit compressible Navier-Stokes element. Call the CalculateRightHandSideInternal() instead.";

            KRATOS_CATCH("")
        }

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

        /**
         * This is called during the assembling process in order
         * to calculate the elemental mass matrix
         * @param rMassMatrix the elemental mass matrix
         * @param rCurrentProcessInfo the current process info instance
         */
        virtual void CalculateMassMatrix(
            MatrixType &rMassMatrix,
            const ProcessInfo &rCurrentProcessInfo) override;

        /**
         * @brief Calculate the lumped mass vector
         * This is called during the assembling process in order
         * to calculate the elemental lumped mass vector
         * @param rLumpedMassVector the elemental lumped mass vector
         * @param rCurrentProcessInfo the current process info instance
         */
        virtual void CalculateLumpedMassVector(
            VectorType &rLumpedMassVector,
            const ProcessInfo &rCurrentProcessInfo) const override;

        /**
         * This function provides the place to perform checks on the completeness of the input.
         * It is designed to be called only once (or anyway, not often) typically at the beginning
         * of the calculations, so to verify that nothing is missing from the input
         * or that no common error is found.
         * @param rCurrentProcessInfo The ProcessInfo of the ModelPart that contains this element.
         * @return 0 if no errors were found.
         */
        int Check(const ProcessInfo &rCurrentProcessInfo) const override;

        void Calculate(
            const Variable<double> &rVariable,
            double &Output,
            const ProcessInfo &rCurrentProcessInfo) override;

        void Calculate(
            const Variable<array_1d<double, 3>> &rVariable,
            array_1d<double, 3> &Output,
            const ProcessInfo &rCurrentProcessInfo) override;

        void Calculate(
            const Variable<Matrix> &rVariable,
            Matrix &Output,
            const ProcessInfo &rCurrentProcessInfo) override;

        void CalculateOnIntegrationPoints(
            const Variable<double> &rVariable,
            std::vector<double> &rOutput,
            const ProcessInfo &rCurrentProcessInfo) override;

        void CalculateOnIntegrationPoints(
            const Variable<array_1d<double, 3>> &rVariable,
            std::vector<array_1d<double, 3>> &rOutput,
            const ProcessInfo &rCurrentProcessInfo) override;

        ///@}
        ///@name Access
        ///@{

        ///@}
        ///@name Inquiry
        ///@{

        ///@}
        ///@name Input and output
        ///@{

        const Parameters GetSpecifications() const override;

        /// Turn back information as a string.
        std::string Info() const override
        {
            return "CompressibleNavierStokesExplicit #";
        }

        /// Print information about this object.
        void PrintInfo(std::ostream &rOStream) const override
        {
            rOStream << Info() << Id();
        }

        /// Print object's data.
        void PrintData(std::ostream &rOStream) const override
        {
            pGetGeometry()->PrintData(rOStream);
        }

    /**
     * @brief Calculate the density projection
     * Auxiliary method to calculate the density projections for the OSS.
     * Note that this method threadsafe adds the elemental RHS values of the L2 projection to the nodes.
     * The division by the lumped mass matrix values requires to be done at the strategy level.
     * @param rCurrentProcessInfo Reference to the current process info
     */
    void CalculateDensityProjection(const ProcessInfo& rCurrentProcessInfo);

        ///@}
    protected:
        ///@name Protected static member variables
        ///@{

        ///@}
        ///@name Protected member Variables
        ///@{

        /**
         * This determines the elemental equation ID vector for all elemental DOFs
         * @param rResult the elemental equation ID vector
         * @param rCurrentProcessInfo the current process info instance
         */
        void EquationIdVector(
            EquationIdVectorType &rResult,
            const ProcessInfo &rCurrentProcessInfo) const override;

        /**
         * Determines the elemental list of DOFs
         * @param ElementalDofList the list of DOFs
         * @param rCurrentProcessInfo the current process info instance
         */
        void GetDofList(
            DofsVectorType &ElementalDofList,
            const ProcessInfo &rCurrentProcessInfo) const override;

        ///@}
        ///@name Protected Operators
        ///@{

        CompressibleNavierStokesExplicit() = default;

        ///@}
        ///@name Protected Operations
        ///@{

        /**
         * @brief Fill element data
         * Auxiliary function to fill the element data structure
         * @param rData Reference to the element data structure to be filled
         * @param rCurrentProcessInfo Reference to the current process info
         */
        void FillElementData(
            ElementDataStruct &rData,
            const ProcessInfo &rCurrentProcessInfo);

        /**
         * @brief Internal CalculateRightHandSide() method
         * This auxiliary RHS calculated method is created to bypass the element API
         * In this way bounded vectors can be used in the explicit residual calculation
         * @param rRightHandSideBoundedVector Reference to the auxiliary RHS vector
         * @param rCurrentProcessInfo Reference to the current process info
         */
        void CalculateRightHandSideInternal(
            BoundedVector<double, BlockSize * TNumNodes> &rRightHandSideBoundedVector,
            const ProcessInfo &rCurrentProcessInfo);

        /**
         * @brief Calculate the momentum projection
         * Auxiliary method to calculate the momentum projections for the OSS.
         * Note that this method threadsafe adds the elemental RHS values of the L2 projection to the nodes.
         * The division by the lumped mass matrix values requires to be done at the strategy level.
         * @param rCurrentProcessInfo Reference to the current process info
         */
        void CalculateMomentumProjection(const ProcessInfo &rCurrentProcessInfo);

        /**
         * @brief Calculate the density projection
         * Auxiliary method to calculate the denstiy projections for the OSS.
         * Note that this method threadsafe adds the elemental RHS values of the L2 projection to the nodes.
         * The division by the lumped mass matrix values requires to be done at the strategy level.
         * @param rCurrentProcessInfo Reference to the current process info
         */
        void CalculateDensityProjection(const ProcessInfo &rCurrentProcessInfo);

        /**
         * @brief Calculate the total energy projection
         * Auxiliary method to calculate the total energy projections for the OSS.
         * Note that this method threadsafe adds the elemental RHS values of the L2 projection to the nodes.
         * The division by the lumped mass matrix values requires to be done at the strategy level.
         * @param rCurrentProcessInfo Reference to the current process info
         */
        void CalculateTotalEnergyProjection(const ProcessInfo &rCurrentProcessInfo);

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

        void save(Serializer &rSerializer) const override
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
        }

        void load(Serializer &rSerializer) override
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
        }

        ///@}
        ///@name Private Operations
        ///@{

        GeometryData::IntegrationMethod GetIntegrationMethod() const override;

        /**
         * @brief Calculate the midpoint velocity divergence
         * This method calculates the velocity divergence in the midpoint of the element
         * @return double Velocity divergence in the midpoint
         */
        double CalculateMidPointVelocityDivergence() const;

        /**
         * @brief Calculate the midpoint sound velocity
         * This method calculates the speed of sound velocity in the midpoint of the element
         * @return double Speed of sound velocity in the midpoint
         */
        double CalculateMidPointSoundVelocity() const;

        /**
         * @brief Calculate the midpoint density gradient
         * This method calculates the gradient of the density in the midpoint of the element
         * @return array_1d<double,3> Density gradient in the midpoint
         */
        array_1d<double, 3> CalculateMidPointDensityGradient() const;

        /**
         * @brief Calculate the midpoint temperature gradient
         * This method calculates the gradient of the temperature in the midpoint of the element
         * @return array_1d<double,3> Temperature gradient in the midpoint
         */
        array_1d<double, 3> CalculateMidPointTemperatureGradient() const;

        /**
         * @brief Calculate the midpoint velocity rotational
         * This method calculates the rotational of the velocity in the midpoint of the element
         * @return array_1d<double,3> Velocity rotational in the midpoint
         */
        array_1d<double, 3> CalculateMidPointVelocityRotational() const;

        /**
         * @brief Calculate the midpoint velocity gradient
         * This method calculates the gradient of the velocity in the midpoint of the element
         * @return BoundedMatrix<double, 3, 3> Velocity gradient in the midpoint
         */
        BoundedMatrix<double, 3, 3> CalculateMidPointVelocityGradient() const;

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
    };
    ///@}
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}

    /// Implementation of template-parameter independent methods

    template <unsigned int TDim, unsigned int TNumNodes>
    int CompressibleNavierStokesExplicit<TDim, TNumNodes>::Check(const ProcessInfo &rCurrentProcessInfo) const
    {
        KRATOS_TRY

        // Perform basic element checks
        int ErrorCode = Kratos::Element::Check(rCurrentProcessInfo);
        if (ErrorCode != 0)
        {
            return ErrorCode;
        }

        // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
        for (unsigned int i = 0; i < TNumNodes; ++i)
        {
            KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].SolutionStepsDataHas(DENSITY)) << "Missing DENSITY variable on solution step data for node " << this->GetGeometry()[i].Id();
            KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].SolutionStepsDataHas(MOMENTUM)) << "Missing MOMENTUM variable on solution step data for node " << this->GetGeometry()[i].Id();
            KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].SolutionStepsDataHas(TOTAL_ENERGY)) << "Missing TOTAL_ENERGY variable on solution step data for node " << this->GetGeometry()[i].Id();
            KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].SolutionStepsDataHas(BODY_FORCE)) << "Missing BODY_FORCE variable on solution step data for node " << this->GetGeometry()[i].Id();
            KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].SolutionStepsDataHas(HEAT_SOURCE)) << "Missing HEAT_SOURCE variable on solution step data for node " << this->GetGeometry()[i].Id();

            // Activate as soon as we start using the explicit DOF based strategy
            KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].HasDofFor(DENSITY)) << "Missing DENSITY DOF in node ", this->GetGeometry()[i].Id();
            KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].HasDofFor(MOMENTUM_X) || this->GetGeometry()[i].HasDofFor(MOMENTUM_Y)) << "Missing MOMENTUM component DOF in node ", this->GetGeometry()[i].Id();
            if constexpr (TDim == 3)
            {
                KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].HasDofFor(MOMENTUM_Z)) << "Missing MOMENTUM component DOF in node ", this->GetGeometry()[i].Id();
            }
            KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].HasDofFor(TOTAL_ENERGY)) << "Missing TOTAL_ENERGY DOF in node ", this->GetGeometry()[i].Id();
        }

        return 0;

        KRATOS_CATCH("");
    }

    template <unsigned int TDim, unsigned int TNumNodes>
    void CompressibleNavierStokesExplicit<TDim, TNumNodes>::Calculate(
        const Variable<double> &rVariable,
        double &Output,
        const ProcessInfo &rCurrentProcessInfo)
    {
        // Lumped projection terms
        if (rVariable == DENSITY_PROJECTION)
        {
            CalculateDensityProjection(rCurrentProcessInfo);
        }
        else if (rVariable == TOTAL_ENERGY_PROJECTION)
        {
            CalculateTotalEnergyProjection(rCurrentProcessInfo);
        }
        else if (rVariable == VELOCITY_DIVERGENCE)
        {
            Output = CalculateMidPointVelocityDivergence();
        }
        else if (rVariable == SOUND_VELOCITY)
        {
            Output = CalculateMidPointSoundVelocity();
        }
        else
        {
            KRATOS_ERROR << "Variable not implemented." << std::endl;
        }
    }

    template <unsigned int TDim, unsigned int TNumNodes>
    void CompressibleNavierStokesExplicit<TDim, TNumNodes>::Calculate(
        const Variable<array_1d<double, 3>> &rVariable,
        array_1d<double, 3> &Output,
        const ProcessInfo &rCurrentProcessInfo)
    {
        if (rVariable == DENSITY_GRADIENT)
        {
            Output = CalculateMidPointDensityGradient();
        }
        else if (rVariable == TEMPERATURE_GRADIENT)
        {
            Output = CalculateMidPointTemperatureGradient();
        }
        else if (rVariable == VELOCITY_ROTATIONAL)
        {
            Output = CalculateMidPointVelocityRotational();
        }
        else if (rVariable == MOMENTUM_PROJECTION)
        {
            CalculateMomentumProjection(rCurrentProcessInfo);
        }
        else
        {
            KRATOS_ERROR << "Variable not implemented." << std::endl;
        }
    }

    template <unsigned int TDim, unsigned int TNumNodes>
    void CompressibleNavierStokesExplicit<TDim, TNumNodes>::Calculate(
        const Variable<Matrix> &rVariable,
        Matrix &Output,
        const ProcessInfo &rCurrentProcessInfo)
    {
        if (rVariable == VELOCITY_GRADIENT)
        {
            Output = CalculateMidPointVelocityGradient();
        }
        else
        {
            KRATOS_ERROR << "Variable not implemented." << std::endl;
        }
    }

    template <unsigned int TDim, unsigned int TNumNodes>
    void CompressibleNavierStokesExplicit<TDim, TNumNodes>::CalculateOnIntegrationPoints(
        const Variable<double> &rVariable,
        std::vector<double> &rOutput,
        const ProcessInfo &rCurrentProcessInfo)
    {
        const auto &r_geometry = GetGeometry();
        const auto &r_integration_points = r_geometry.IntegrationPoints();
        if (rOutput.size() != r_integration_points.size())
        {
            rOutput.resize(r_integration_points.size());
        }

        if (rVariable == SHOCK_SENSOR)
        {
            const double sc = this->GetValue(SHOCK_SENSOR);
            for (unsigned int i_gauss = 0; i_gauss < r_integration_points.size(); ++i_gauss)
            {
                rOutput[i_gauss] = sc;
            }
        }
        else if (rVariable == SHEAR_SENSOR)
        {
            const double sc = this->GetValue(SHEAR_SENSOR);
            for (unsigned int i_gauss = 0; i_gauss < r_integration_points.size(); ++i_gauss)
            {
                rOutput[i_gauss] = sc;
            }
        }
        else if (rVariable == THERMAL_SENSOR)
        {
            const double sc = this->GetValue(THERMAL_SENSOR);
            for (unsigned int i_gauss = 0; i_gauss < r_integration_points.size(); ++i_gauss)
            {
                rOutput[i_gauss] = sc;
            }
        }
        else if (rVariable == ARTIFICIAL_CONDUCTIVITY)
        {
            const double k_star = this->GetValue(ARTIFICIAL_CONDUCTIVITY);
            for (unsigned int i_gauss = 0; i_gauss < r_integration_points.size(); ++i_gauss)
            {
                rOutput[i_gauss] = k_star;
            }
        }
        else if (rVariable == ARTIFICIAL_BULK_VISCOSITY)
        {
            const double beta_star = this->GetValue(ARTIFICIAL_BULK_VISCOSITY);
            for (unsigned int i_gauss = 0; i_gauss < r_integration_points.size(); ++i_gauss)
            {
                rOutput[i_gauss] = beta_star;
            }
        }
        else if (rVariable == VELOCITY_DIVERGENCE)
        {
            const double div_v = CalculateMidPointVelocityDivergence();
            for (unsigned int i_gauss = 0; i_gauss < r_integration_points.size(); ++i_gauss)
            {
                rOutput[i_gauss] = div_v;
            }
        }
        else
        {
            KRATOS_ERROR << "Variable not implemented." << std::endl;
        }
    }

    template <unsigned int TDim, unsigned int TNumNodes>
    void CompressibleNavierStokesExplicit<TDim, TNumNodes>::CalculateOnIntegrationPoints(
        const Variable<array_1d<double, 3>> &rVariable,
        std::vector<array_1d<double, 3>> &rOutput,
        const ProcessInfo &rCurrentProcessInfo)
    {
        const auto &r_geometry = GetGeometry();
        const auto &r_integration_points = r_geometry.IntegrationPoints();
        if (rOutput.size() != r_integration_points.size())
        {
            rOutput.resize(r_integration_points.size());
        }

        if (rVariable == DENSITY_GRADIENT)
        {
            const array_1d<double, 3> rho_grad = CalculateMidPointDensityGradient();
            for (unsigned int i_gauss = 0; i_gauss < r_integration_points.size(); ++i_gauss)
            {
                rOutput[i_gauss] = rho_grad;
            }
        }
        else if (rVariable == TEMPERATURE_GRADIENT)
        {
            const array_1d<double, 3> temp_grad = CalculateMidPointTemperatureGradient();

            for (unsigned int i_gauss = 0; i_gauss < r_integration_points.size(); ++i_gauss)
            {
                rOutput[i_gauss] = temp_grad;
            }
        }
        else if (rVariable == VELOCITY_ROTATIONAL)
        {
            const array_1d<double, 3> rot_v = CalculateMidPointVelocityRotational();
            for (unsigned int i_gauss = 0; i_gauss < r_integration_points.size(); ++i_gauss)
            {
                rOutput[i_gauss] = rot_v;
            }
        }
        else
        {
            KRATOS_ERROR << "Variable not implemented." << std::endl;
        }
    }

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
    namespace CompressibleNavierStokesExplicitInternal
    {
        template <unsigned int TDim, unsigned int TNumNodes>
        using ElementDataStruct = typename CompressibleNavierStokesExplicit<TDim, TNumNodes>::ElementDataStruct;

        constexpr bool IsSimplex(const unsigned int dimensions, const unsigned int nnodes)
        {
            return dimensions == nnodes - 1;
        }

        // Specialization for simplex geometries
        template <unsigned int TDim, unsigned int TNumNodes>
        static typename std::enable_if<IsSimplex(TDim, TNumNodes), void>::type ComputeGeometryData(
            const Geometry<Node> &rGeometry,
            ElementDataStruct<TDim, TNumNodes> &rData)
        {
            GeometryUtils::CalculateGeometryData(rGeometry, rData.DN_DX, rData.N, rData.volume);
            rData.h = ElementSizeCalculator<TDim, TNumNodes>::GradientsElementSize(rData.DN_DX);
        }

        /**
         * Specialization for geometries with non-uniform jacobian (non-simplex)
         * Shape functions cannot be obtained here. They will need to be calculated
         * during integration at each gauss point.
         */
        template <unsigned int TDim, unsigned int TNumNodes>
        static typename std::enable_if<!IsSimplex(TDim, TNumNodes), void>::type ComputeGeometryData(
            const Geometry<Node> &rGeometry,
            ElementDataStruct<TDim, TNumNodes> &rData)
        {
            rData.volume = rGeometry.DomainSize();
            rData.h = ElementSizeCalculator<TDim, TNumNodes>::AverageElementSize(rGeometry);
        }
    }

    template <unsigned int TDim, unsigned int TNumNodes>
    void CompressibleNavierStokesExplicit<TDim, TNumNodes>::FillElementData(
        ElementDataStruct &rData,
        const ProcessInfo &rCurrentProcessInfo)
    {
        // Getting data for the given geometry
        const auto &r_geometry = GetGeometry();

        // Loads shape function info only if jacobian is uniform
        CompressibleNavierStokesExplicitInternal::ComputeGeometryData<TDim, TNumNodes>(r_geometry, rData);

        // Database access to all of the variables needed
        Properties &r_properties = this->GetProperties();
        rData.mu = r_properties.GetValue(DYNAMIC_VISCOSITY);
        rData.lambda = r_properties.GetValue(CONDUCTIVITY);
        rData.c_v = r_properties.GetValue(SPECIFIC_HEAT); // TODO: WE SHOULD SPECIFY WHICH ONE --> CREATE SPECIFIC_HEAT_CONSTANT_VOLUME
        rData.gamma = r_properties.GetValue(HEAT_CAPACITY_RATIO);

        rData.UseOSS = rCurrentProcessInfo[OSS_SWITCH];
        rData.ShockCapturing = rCurrentProcessInfo[SHOCK_CAPTURING_SWITCH];

        // Magnitudes to calculate the time derivatives
        const double time_step = rCurrentProcessInfo[DELTA_TIME];
        const double theta = rCurrentProcessInfo[TIME_INTEGRATION_THETA];
        const double aux_theta = theta > 0 ? 1.0 / (theta * time_step) : 0.0;

        // Get nodal values
        if (rData.UseOSS)
        {
            for (unsigned int i = 0; i < TNumNodes; ++i)
            {
                const auto &r_node = r_geometry[i];
                // Vector data
                const array_1d<double, 3> &r_momentum = r_node.FastGetSolutionStepValue(MOMENTUM);
                const array_1d<double, 3> &r_momentum_old = r_node.FastGetSolutionStepValue(MOMENTUM, 1);
                const array_1d<double, 3> &r_momentum_projection = r_node.GetValue(MOMENTUM_PROJECTION);
                const array_1d<double, 3> mom_inc = r_momentum - r_momentum_old;
                const auto &r_body_force = r_node.FastGetSolutionStepValue(BODY_FORCE);
                for (unsigned int k = 0; k < TDim; ++k)
                {
                    rData.U(i, k + 1) = r_momentum[k];
                    rData.dUdt(i, k + 1) = aux_theta * mom_inc[k];
                    rData.ResProj(i, k + 1) = r_momentum_projection[k];
                    rData.f_ext(i, k) = r_body_force[k];
                }
                // Density data
                const double &r_rho = r_node.FastGetSolutionStepValue(DENSITY);
                const double &r_rho_old = r_node.FastGetSolutionStepValue(DENSITY, 1);
                const double rho_inc = r_rho - r_rho_old;
                rData.U(i, 0) = r_rho;
                rData.dUdt(i, 0) = aux_theta * rho_inc;
                rData.ResProj(i, 0) = r_node.GetValue(DENSITY_PROJECTION);
                // Total energy data
                const double &r_tot_ener = r_node.FastGetSolutionStepValue(TOTAL_ENERGY);
                const double &r_tot_ener_old = r_node.FastGetSolutionStepValue(TOTAL_ENERGY, 1);
                const double tot_ener_inc = r_tot_ener - r_tot_ener_old;
                rData.U(i, TDim + 1) = r_tot_ener;
                rData.dUdt(i, TDim + 1) = aux_theta * tot_ener_inc;
                rData.ResProj(i, TDim + 1) = r_node.GetValue(TOTAL_ENERGY_PROJECTION);
                // Source data
                rData.r_ext(i) = r_node.FastGetSolutionStepValue(HEAT_SOURCE);
                rData.m_ext(i) = r_node.FastGetSolutionStepValue(MASS_SOURCE);
                // Shock capturing data
                rData.alpha_sc_nodes(i) = r_node.GetValue(ARTIFICIAL_MASS_DIFFUSIVITY);
                rData.mu_sc_nodes(i) = r_node.GetValue(ARTIFICIAL_DYNAMIC_VISCOSITY);
                rData.beta_sc_nodes(i) = r_node.GetValue(ARTIFICIAL_BULK_VISCOSITY);
                rData.lamb_sc_nodes(i) = r_node.GetValue(ARTIFICIAL_CONDUCTIVITY);
            }
        }
        else
        {
            for (unsigned int i = 0; i < TNumNodes; ++i)
            {
                const auto &r_node = r_geometry[i];
                // Vector data
                const array_1d<double, 3> &r_momentum = r_node.FastGetSolutionStepValue(MOMENTUM);
                const array_1d<double, 3> &r_momentum_old = r_node.FastGetSolutionStepValue(MOMENTUM, 1);
                const array_1d<double, 3> mom_inc = r_momentum - r_momentum_old;
                const auto &r_body_force = r_node.FastGetSolutionStepValue(BODY_FORCE);
                for (unsigned int k = 0; k < TDim; ++k)
                {
                    rData.U(i, k + 1) = r_momentum[k];
                    rData.dUdt(i, k + 1) = aux_theta * mom_inc[k];
                    rData.f_ext(i, k) = r_body_force[k];
                }
                // Density data
                const double &r_rho = r_node.FastGetSolutionStepValue(DENSITY);
                const double &r_rho_old = r_node.FastGetSolutionStepValue(DENSITY, 1);
                rData.U(i, 0) = r_rho;
                rData.dUdt(i, 0) = aux_theta * (r_rho - r_rho_old);
                // Total energy data
                const double &r_tot_ener = r_node.FastGetSolutionStepValue(TOTAL_ENERGY);
                const double &r_tot_ener_old = r_node.FastGetSolutionStepValue(TOTAL_ENERGY, 1);
                rData.U(i, TDim + 1) = r_tot_ener;
                rData.dUdt(i, TDim + 1) = aux_theta * (r_tot_ener - r_tot_ener_old);
                // Source data
                rData.r_ext(i) = r_node.FastGetSolutionStepValue(HEAT_SOURCE);
                rData.m_ext(i) = r_node.FastGetSolutionStepValue(MASS_SOURCE);
                // Shock capturing data
                rData.alpha_sc_nodes(i) = r_node.GetValue(ARTIFICIAL_MASS_DIFFUSIVITY);
                rData.mu_sc_nodes(i) = r_node.GetValue(ARTIFICIAL_DYNAMIC_VISCOSITY);
                rData.beta_sc_nodes(i) = r_node.GetValue(ARTIFICIAL_BULK_VISCOSITY);
                rData.lamb_sc_nodes(i) = r_node.GetValue(ARTIFICIAL_CONDUCTIVITY);
            }
        }
    }

    template <unsigned int TDim, unsigned int TNumNodes>
    array_1d<double, 3> CompressibleNavierStokesExplicit<TDim, TNumNodes>::CalculateMidPointDensityGradient() const
    {
        // Get geometry data
        const auto &r_geom = GetGeometry();
        const unsigned int NumNodes = r_geom.PointsNumber();
        Geometry<Node>::ShapeFunctionsGradientsType dNdX_container;
        r_geom.ShapeFunctionsIntegrationPointsGradients(dNdX_container, GeometryData::IntegrationMethod::GI_GAUSS_1);
        const auto &r_dNdX = dNdX_container[0];

        // Calculate midpoint magnitudes
        array_1d<double, 3> midpoint_grad_rho = ZeroVector(3);
        for (unsigned int i_node = 0; i_node < NumNodes; ++i_node)
        {
            auto &r_node = r_geom[i_node];
            const auto node_dNdX = row(r_dNdX, i_node);
            const double &r_rho = r_node.FastGetSolutionStepValue(DENSITY);
            for (unsigned int d1 = 0; d1 < TDim; ++d1)
            {
                midpoint_grad_rho[d1] += node_dNdX(d1) * r_rho;
            }
        }

        return midpoint_grad_rho;
    }

    // TODO: IN HERE WE HAVE LINEARIZED THE DERIVATIVE... CHECK IF WE SHOULD PROPERLY COMPUTE IT
    template <unsigned int TDim, unsigned int TNumNodes>
    array_1d<double, 3> CompressibleNavierStokesExplicit<TDim, TNumNodes>::CalculateMidPointTemperatureGradient() const
    {
        // Get geometry data
        const auto &r_geom = GetGeometry();
        const unsigned int NumNodes = r_geom.PointsNumber();
        Geometry<Node>::ShapeFunctionsGradientsType dNdX_container;
        r_geom.ShapeFunctionsIntegrationPointsGradients(dNdX_container, GeometryData::IntegrationMethod::GI_GAUSS_1);
        const auto &r_dNdX = dNdX_container[0];

        // Calculate midpoint magnitudes
        const double c_v = GetProperties()[SPECIFIC_HEAT];
        array_1d<double, 3> midpoint_grad_temp = ZeroVector(3);
        for (unsigned int i_node = 0; i_node < NumNodes; ++i_node)
        {
            auto &r_node = r_geom[i_node];
            const auto node_dNdX = row(r_dNdX, i_node);
            const auto &r_mom = r_node.FastGetSolutionStepValue(MOMENTUM);
            const double &r_rho = r_node.FastGetSolutionStepValue(DENSITY);
            const double &r_tot_ener = r_node.FastGetSolutionStepValue(TOTAL_ENERGY);
            const array_1d<double, 3> vel = r_mom / r_rho;
            const double temp = (r_tot_ener / r_rho - 0.5 * inner_prod(vel, vel)) / c_v;
            for (unsigned int d1 = 0; d1 < TDim; ++d1)
            {
                midpoint_grad_temp[d1] += node_dNdX(d1) * temp;
            }
        }

        return midpoint_grad_temp;
    }

    template <unsigned int TDim, unsigned int TNumNodes>
    double CompressibleNavierStokesExplicit<TDim, TNumNodes>::CalculateMidPointSoundVelocity() const
    {
        // Get geometry data
        const auto &r_geom = GetGeometry();
        const unsigned int NumNodes = r_geom.PointsNumber();

        // Calculate midpoint magnitudes
        double midpoint_rho = 0.0;
        double midpoint_tot_ener = 0.0;
        array_1d<double, TDim> midpoint_mom = ZeroVector(TDim);
        for (unsigned int i_node = 0; i_node < NumNodes; ++i_node)
        {
            auto &r_node = r_geom[i_node];
            const auto &r_mom = r_node.FastGetSolutionStepValue(MOMENTUM);
            const double &r_rho = r_node.FastGetSolutionStepValue(DENSITY);
            const double &r_tot_ener = r_node.FastGetSolutionStepValue(TOTAL_ENERGY);
            midpoint_rho += r_rho;
            midpoint_tot_ener += r_tot_ener;
            for (unsigned int d1 = 0; d1 < TDim; ++d1)
            {
                midpoint_mom[d1] += r_mom(d1);
            }
        }
        midpoint_rho /= NumNodes;
        midpoint_mom /= NumNodes;
        midpoint_tot_ener /= NumNodes;

        // Calculate midpoint speed of sound
        const auto &r_prop = GetProperties();
        const double c_v = r_prop.GetValue(SPECIFIC_HEAT);
        const double gamma = r_prop.GetValue(HEAT_CAPACITY_RATIO);
        const double temp = (midpoint_tot_ener / midpoint_rho - inner_prod(midpoint_mom, midpoint_mom) / (2 * std::pow(midpoint_rho, 2))) / c_v;
        double midpoint_c = std::sqrt(gamma * (gamma - 1.0) * c_v * temp);
        return midpoint_c;
    }

    template <unsigned int TDim, unsigned int TNumNodes>
    double CompressibleNavierStokesExplicit<TDim, TNumNodes>::CalculateMidPointVelocityDivergence() const
    {
        // Get geometry data
        const auto &r_geom = GetGeometry();
        const unsigned int NumNodes = r_geom.PointsNumber();
        Geometry<Node>::ShapeFunctionsGradientsType dNdX_container;
        r_geom.ShapeFunctionsIntegrationPointsGradients(dNdX_container, GeometryData::IntegrationMethod::GI_GAUSS_1);
        const auto &r_dNdX = dNdX_container[0];

        // Calculate midpoint magnitudes
        double midpoint_rho = 0.0;
        double midpoint_div_mom = 0.0;
        array_1d<double, TDim> midpoint_mom = ZeroVector(TDim);
        array_1d<double, TDim> midpoint_grad_rho = ZeroVector(TDim);
        for (unsigned int i_node = 0; i_node < NumNodes; ++i_node)
        {
            auto &r_node = r_geom[i_node];
            const auto node_dNdX = row(r_dNdX, i_node);
            const auto &r_mom = r_node.FastGetSolutionStepValue(MOMENTUM);
            const double &r_rho = r_node.FastGetSolutionStepValue(DENSITY);
            midpoint_rho += r_rho;
            for (unsigned int d1 = 0; d1 < TDim; ++d1)
            {
                midpoint_mom[d1] += r_mom(d1);
                midpoint_div_mom += node_dNdX(d1) * r_mom(d1);
                midpoint_grad_rho[d1] += node_dNdX(d1) * r_rho;
            }
        }
        midpoint_rho /= NumNodes;
        midpoint_mom /= NumNodes;

        // Calculate velocity divergence
        // Note that the formulation is written in conservative variables. Hence we do div(mom/rho).
        double midpoint_div_v = (midpoint_rho * midpoint_div_mom - inner_prod(midpoint_mom, midpoint_grad_rho)) / std::pow(midpoint_rho, 2);
        return midpoint_div_v;
    }

    template <unsigned int TDim, unsigned int TNumNodes>
    void CompressibleNavierStokesExplicit<TDim, TNumNodes>::CalculateLumpedMassVector(
        VectorType &rLumpedMassVector,
        const ProcessInfo &rCurrentProcessInfo) const
    {
        // Initialize the lumped mass vector
        constexpr IndexType size = TNumNodes * BlockSize;
        if (rLumpedMassVector.size() != BlockSize)
        {
            rLumpedMassVector.resize(size, false);
        }

        // Fill the lumped mass vector
        const double nodal_mass = GetGeometry().DomainSize() / TNumNodes;
        std::fill(rLumpedMassVector.begin(), rLumpedMassVector.end(), nodal_mass);
    }

    template <unsigned int TDim, unsigned int TNumNodes>
    void CompressibleNavierStokesExplicit<TDim, TNumNodes>::AddExplicitContribution(const ProcessInfo &rCurrentProcessInfo)
    {
        // Calculate the explicit residual vector
        BoundedVector<double, DofSize> rhs;
        CalculateRightHandSideInternal(rhs, rCurrentProcessInfo);

        // Add the residual contribution
        // Note that the reaction is indeed the formulation residual
        auto &r_geometry = GetGeometry();

        for (IndexType i_node = 0; i_node < NumNodes; ++i_node)
        {
            auto &r_node = r_geometry[i_node];

            IndexType i_dof = BlockSize * i_node;

            AtomicAdd(r_node.FastGetSolutionStepValue(REACTION_DENSITY), rhs[i_dof++]);

            for (IndexType d = 0; d < Dim; ++d)
            {
                AtomicAdd(r_node.FastGetSolutionStepValue(REACTION)[d], rhs[i_dof++]);
            }

            AtomicAdd(r_node.FastGetSolutionStepValue(REACTION_ENERGY), rhs[i_dof]);
        }
    }

    template <unsigned int TDim, unsigned int TNumNodes>
    const Parameters CompressibleNavierStokesExplicit<TDim, TNumNodes>::GetSpecifications() const
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
        "required_variables"         : ["DENSITY","MOMENTUM","TOTAL_ENERGY","BODY_FORCE","HEAT_SOURCE"],
        "required_dofs"              : [],
        "flags_used"                 : [],
        "compatible_geometries"      : ["Triangle2D3","Quadrilateral2D4","Tetrahedra3D4"],
        "element_integrates_in_time" : true,
        "compatible_constitutive_laws": {
            "type"        : [],
            "dimension"   : [],
            "strain_size" : []
        },
        "required_polynomial_degree_of_geometry" : 1,
        "documentation"   :
            "This element implements a compressible Navier-Stokes formulation written in conservative variables. A Variational MultiScales (VMS) stabilization technique, both with Algebraic SubGrid Scales (ASGS) and Orthogonal Subgrid Scales (OSS), is used. This element is compatible with both entropy-based and physics-based shock capturing techniques."
    })");

        if constexpr (TDim == 2)
        {
            std::vector<std::string> dofs_2d({"DENSITY", "MOMENTUM_X", "MOMENTUM_Y", "TOTAL_ENERGY"});
            specifications["required_dofs"].SetStringArray(dofs_2d);
        }
        else
        {
            std::vector<std::string> dofs_3d({"DENSITY", "MOMENTUM_X", "MOMENTUM_Y", "MOMENTUM_Z", "TOTAL_ENERGY"});
            specifications["required_dofs"].SetStringArray(dofs_3d);
        }

        return specifications;
    }

} // namespace Kratos.

#endif // KRATOS_COMPRESSIBLE_NAVIER_STOKES_EXPLICIT_H_INCLUDED  defined
