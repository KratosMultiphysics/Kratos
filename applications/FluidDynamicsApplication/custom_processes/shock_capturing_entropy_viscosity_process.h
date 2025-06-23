//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Eduard Gómez
//

#if !defined(KRATOS_SHOCK_CAPTURING_ENTROPY_VISCOSITY_H_INCLUDED)
#define  KRATOS_SHOCK_CAPTURING_ENTROPY_VISCOSITY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "processes/process.h"

// Application includes


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
 * Shock Capturing Process for compressible navier stokes. Source:
 * 
 * 2011. Guermond, Pasquetti, Popov. Entropy viscosity method for nonlinear conservation laws
 * Journal of Computational Physics. Volume 230, Issue 11, 20 May 2011, Pages 4248-4267
 * 
 */
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) ShockCapturingEntropyViscosityProcess : public Process
{
public:

    ///@name Type Definitions
    ///@{

    typedef ModelPart::NodeType NodeType;

    struct InfNormData
    {
        double EntropyResidual;
        double Density;
        double TotalVelocity;
    };


    /**
     * @brief Small utility to compute the total derivative of a magnitude
     *
     * @param Value is the value of the magnitude at the current step
     * @param TimeDerivative is the value of the derivative of the magnitude in this step
     * @param Flux is the volumetric flowrate of the magnitude
     */
    struct TotalDerivativeUtil
    {
        Vector Value;
        Vector TimeDerivative;
        Matrix Flux;

        TotalDerivativeUtil()
            : TotalDerivativeUtil(0, 0)
        {}

        TotalDerivativeUtil(
            const std::size_t NumberOfDimensions,
            const std::size_t NumberOfNodes)
            : Value          {NumberOfNodes, 0.0},
              TimeDerivative {NumberOfNodes, 0.0},
              Flux           {NumberOfNodes, NumberOfDimensions, 0.0}
        {
        }

        void LoadNodalValues(
            const Variable<double>& rVariable,
            const NodeType& rNode,
            const std::size_t NodeIndex,
            const array_1d<double, 3>& Velocity,
            const double DeltaTime)
        {
            KRATOS_TRY

            const auto old_value = rNode.FastGetSolutionStepValue(rVariable, 1);
            Value[NodeIndex] = rNode.FastGetSolutionStepValue(rVariable);
            TimeDerivative[NodeIndex] = (Value[NodeIndex] - old_value) / DeltaTime;

            Flux(NodeIndex, 0) = Value[NodeIndex] * Velocity[0];
            Flux(NodeIndex, 1) = Value[NodeIndex] * Velocity[1];
            if(Flux.size2()==3){ // (if 3D)
                Flux(NodeIndex, 2) = Value[NodeIndex] * Velocity[2];
            }

            KRATOS_CATCH("")
        }

        double ComputeAtGaussPoint(
            const MatrixRow<const Matrix>& rShapeFunctions,
            const Matrix& rShapeFunctionsGradents) const
        {
            KRATOS_TRY

            const double divergence = Divergence(rShapeFunctionsGradents, Flux);
            const double time_derivative = inner_prod(TimeDerivative, rShapeFunctions);
            return time_derivative + divergence;

            KRATOS_CATCH("") // Catching possible error in divergence computation
        }

    private:
        /**
         * @brief Computes the divergence
         *
         * @param rShapeFunGradients The gradients of the shape functions [ndims x nnodes]
         * @param rNodalValues Values of the magnitude at the nodes. [ndims x nnodes]
         * @return Result
         */
        static double Divergence(const Matrix& rShapeFunGradients, const Matrix& rNodalValues);
    };


    /// Pointer definition of ShockCapturingEntropyViscosityProcess
    KRATOS_CLASS_POINTER_DEFINITION(ShockCapturingEntropyViscosityProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with model
    ShockCapturingEntropyViscosityProcess(
        Model& rModel,
        Parameters rParameters)
        : ShockCapturingEntropyViscosityProcess(rModel.GetModelPart(rParameters["model_part_name"].GetString()), rParameters) {};

    /// Constructor with model part
    ShockCapturingEntropyViscosityProcess(
        ModelPart& rModelPart,
        Parameters rParameters);

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void ExecuteInitializeSolutionStep() override;

    int Check() override;

    const Parameters GetDefaultParameters() const override;

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
        buffer << "ShockCapturingEntropyViscosityProcess";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override { rOStream << "ShockCapturingEntropyViscosityProcess"; }

    /// Print object's data.
    void PrintData(std::ostream&) const override {}

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

    ModelPart& mrModelPart;
    bool mIsInitialized = false;
    bool mComputeAreasEveryStep = false;
    double mEntropyConstant = 0.0;
    double mEnergyConstant = 0.0;
    double mArtificialMassDiffusivityPrandtl = 0.0;
    double mArtificialConductivityPrandtl = 0.0;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void ValidateAndAssignParameters(Parameters &rParameters);

    void UpdateNodalAreaProcess();

    /**
     * @brief Computes nodal entropies and initializes artificial variables to 0.
     *
     * @param WriteBufferIndex: The buffer index to write the result in. This is an
     * ugly hack to deal with the first time-step.
     */
    void ComputeNodalEntropies(const unsigned int WriteBufferIndex = 0);

    void ComputeArtificialMagnitudes();

    void DistributeVariablesToNodes(
        Element& rElement,
        const double ArtificialDynamicViscosity,
        const double ArtificialMassDiffusivity,
        const double ArtificialConductivity,
        const std::function<double(Geometry<Node>*)>& rGeometrySize) const;

    static double ComputeEntropy(
        const double Density, 
        const double Pressure, 
        const double Gamma)
    {
        return Density / (Gamma - 1.0) * std::log(Pressure / std::pow(Density, Gamma));
    }

    /**
     * @brief Computes the square of the element's shortest edge
     *
     * @param rElement
     * @return double
     */
    template <unsigned int TDim>
    static double MinimumEdgeLengthSquared(const Element &rElement);

    
    /**
     * @brief Computes the infinity norm of entropy and residual
     *
     * @param rElement: The element to compute them on
     * @param DeltaTime: The time step size
     * @return Tuple containing the inf norms of {entropy residual, density}
     */
    static InfNormData ComputeElementalInfNormData(
        const Element& rElement,
        const double DeltaTime,
        const double HeatCapacityRatio,
        const double SpecificHeatCV);

    /**
     * @brief Builds the TotalDerivativeUtil objects that will be used to compute inf norms
     *
     * @param rElement: The element to compute them on
     * @return Tuple containing {TotalDerivativeUtil for entroy, TotalDerivativeUtil for density, Vector with total velocities}
     */
    static std::tuple<TotalDerivativeUtil, TotalDerivativeUtil, Vector> BuildTotalDerivativeUtils(
        const Element& rElement,
        const double DeltaTime,
        const double HeatCapacityRatio,
        const double SpecificHeatCV);

    /**
     * @brief Computes entropy max value of residual, density and total velocity over all gauss points.
     *
     */
    static InfNormData ComputeInfNorms(
        const Geometry<NodeType>& rGeometry,
        const TotalDerivativeUtil& rEntropyTotalDerivative,
        const TotalDerivativeUtil& rDensityTotalDerivative,
        const Vector& rTotalVelocities);



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
    ShockCapturingEntropyViscosityProcess& operator=(ShockCapturingEntropyViscosityProcess const& rOther);

    /// Copy constructor.
    ShockCapturingEntropyViscosityProcess(ShockCapturingEntropyViscosityProcess const& rOther);

    ///@}

}; // Class ShockCapturingEntropyViscosityProcess

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const ShockCapturingEntropyViscosityProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    return rOStream;
}

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SHOCK_CAPTURING_ENTROPY_VISCOSITY_H_INCLUDED  defined
