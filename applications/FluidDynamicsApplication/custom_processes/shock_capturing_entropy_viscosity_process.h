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

class KRATOS_API(FLUID_DYNAMICS_APPLICATION) ShockCapturingEntropyViscosityProcess : public Process
{
public:

    ///@name Type Definitions
    ///@{

    typedef ModelPart::NodeType NodeType;

    struct InfNormData
    {
        double Entropy;
        double Density;
        double TotalVelocity;
    };


    class ElementAreaComputer
    {
    public:
        void SetDimension(unsigned int Dimension)
        {
            KRATOS_TRY

            switch(Dimension)
            {
                case 2: mComputeImpl = &ComputeGeometryVolume<2>;
                case 3: mComputeImpl = &ComputeGeometryVolume<3>;
                default: KRATOS_ERROR << "Only dimensions 1 and 2 supported" << std::endl;
            }

            KRATOS_CATCH("")
        }

        double operator()(const Element& rElement) const
        {
            return mComputeImpl(rElement.GetGeometry());
        }

        template<unsigned int TDim>
        static double ComputeGeometryVolume(const Geometry<NodeType>& rGeometry);

    private:
        auto mComputeImpl = &ComputeGeometryVolume<2>;
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

        VariableTotalDerivativeUtil(
            const std::size_t NumberOfDimensions,
            const std::size_t NumberOfNodes)
            : Value{NumberOfNodes, 0.0},
            TimeDerivative{NumberOfNodes, 0.0},
            Flux{NumberOfNodes, NumberOfDimensions, 0.0}
        {
        }

        LoadNodalValues(
            const Variable<double>& rVariable,
            const NodeType& rNode,
            const std::size_t NodeIndex,
            const array_1d<double, 3>& Velocity,
            const double DeltaTime)
        {
            KRATOS_TRY

            const auto old_value = rNode.FastGetSolutionStepValue(rVariable, 1);
            Value = rNode.FastGetSolutionStepValue(rVariable);
            TimeDerivative[NodeIndex] = (Value - old_value) / DeltaTime;
            column(Flux, i) = Value * Velocity;

            KRATOS_CATCH("")
        }

        ComputeAtGaussPoint(
            const Row<Matrix>& rShapeFunctions,
            const Matrix& rShapeFunctionsGradents) const
        {
            KRATOS_TRY

            const double divergence = Divergence(rShapeFunctionsGradents, EntropyRD.Flux);
            const double time_derivative = inner_prod(EntropyRD.TimeDerivative, rShapeFunctions);
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
        double Divergence(const Matrix& rShapeFunGradients, const Matrix& rNodalValues);
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
        Parameters rParameters)
        : Process()
        , mrModelPart(rModelPart);

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void ExecuteBeforeSolutionLoop() override;

    void ExecuteFinalizeSolutionStep() override;

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
    void PrintInfo(std::ostream &rOStream) const override { rOStream << "ShockCapturingEntropyViscosityProcess"; }

    /// Print object's data.
    void PrintData(std::ostream &rOStream) const override {}

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
    bool mComputeAreasEveryStep = false;
    double mTunableConstant = 0.0;
    double mTunableConstantMax = 0.0;
    double mArtificialBulkViscosityPrandtl = 0.0;
    double mArtificialConductivityPrandtl = 0.0;

    ElementVolumeComputer mElementVolumeComputer;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void ValidateAndAssignParameters(Parameters &rParameters);

    void ShockCapturingEntropyViscosityProcess::UpdateMeshDependentData();

    /**
     * @brief Computes nodal entropies and initializes artificial variables to 0
     */
    void ComputeNodalEntropies();


    void ComputeArtificialMagnitudes();
    
    static double ComputeEntropy(const double Density, const double Pressure, const double Gamma);
    
    /**
     * @brief Computes the square of the element's shortest edge
     * 
     * @param rElement 
     * @return double 
     */
    static double ComputeHSquared(const Element& rElement);

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
        const double HeatCapacityRatio);

    /**
     * @brief Buidls the TotalDerivativeUtil objects that will be used to compute inf norms
     * 
     * @param rElement: The element to compute them on
     * @return Tuple containing {TotalDerivativeUtil for entroy, TotalDerivativeUtil for density, Vector with total velocities}
     */
    static std::tuple<TotalDerivativeUtil, TotalDerivativeUtil, Vector> BuildTotalDerivativeUtils(
        const Element& rElement,
        const double DeltaTime,
        const double HeatCapacityRatio);
    
    /**
     * @brief Computes entropy max value of residual, density and total velocity over all gauss points.
     *  
     */
    static InfNormData ComputeInfNorms(
        const TotalDerivativeUtil& EntropyTotalDerivative,
        const TotalDerivativeUtil& DensityTotalDerivative,
        const Vector& TotalVelocities);



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
    const ShockCapturingEntropyViscosityProcess& rThis);

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SHOCK_CAPTURING_ENTROPY_VISCOSITY_H_INCLUDED  defined
