//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Ruben Zorrilla
//

#if !defined(KRATOS_SHOCK_CAPTURING_UTILITIES_H_INCLUDED)
#define  KRATOS_SHOCK_CAPTURING_UTILITIES_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"

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

class KRATOS_API(FLUID_DYNAMICS_APPLICATION) ShockCapturingProcess : public Process
{
public:

    ///@name Type Definitions
    ///@{

    /// The base process type
    typedef Process BaseType;

    /// Type for the metric calculation function
    typedef std::function<std::tuple<double, double, Matrix>(Geometry<Node<3>> &rGeometry)> ElementMetricFunctionType;

    /// Type for the 2D (linear triangle) TLS shock capturing container
    typedef std::tuple<double, Matrix, array_1d<double, 3>, array_1d<double, 3>, array_1d<double, 3>> ShockCapturingTLSType2D;

    /// Type for the 3D (linear tetrahedra) TLS shock capturing container
    typedef std::tuple<double, Matrix, array_1d<double, 3>, array_1d<double, 3>, array_1d<double, 3>> ShockCapturingTLSType3D;

    /// Pointer definition of ShockCapturingProcess
    KRATOS_CLASS_POINTER_DEFINITION(ShockCapturingProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with model
    ShockCapturingProcess(
        Model& rModel,
        Parameters& rParameters)
        : mrModelPart(rModel.GetModelPart(rParameters["model_part_name"].GetString()))
    {
        ValidateAndAssignParameters(rParameters);
    };

    /// Constructor with model part
    ShockCapturingProcess(
        ModelPart& rModelPart,
        Parameters& rParameters)
        : mrModelPart(rModelPart)
    {
        ValidateAndAssignParameters(rParameters);
    };

    /// Destructor.
    ~ShockCapturingProcess() {};

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void Execute() override;

    void ExecuteInitialize() override;

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
        buffer << "ShockCapturingProcess";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override { rOStream << "ShockCapturingProcess"; }

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

    bool mUpdateNodalArea;
    bool mShockSensor;
    bool mShearSensor;
    bool mThermalSensor;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void ValidateAndAssignParameters(Parameters &rParameters);

    /**
     * @brief Set the element metric function object
     * This method checks the first element type and sets the function to calculate
     * the metric tensor accordingly. Note that, for the sake of performance, it is
     * assumed that the element geometry type in the mesh is unique.
     * Such metric tensor calculation function returns a tuple that contains
     * 0 : The reference element size used to build the metric tensor
     * 1 : The infimum norm of the metric tensor
     * 2 : The metric tensor relative to the reference element size stored in 0
     * @return ElementMetricFunctionType Function to calculate the metric tensor
     */
    ElementMetricFunctionType SetElementMetricFunction();

    ElementMetricFunctionType SetShockCapturingTLSContainer();

    /**
     * @brief Physics-based shock capturing
     * This function calculates the artificial magnitudes using a physics-based shock capturing method.
     * References https://arc.aiaa.org/doi/abs/10.2514/6.2018-0062
     */
    void CalculatePhysicsBasedShockCapturing();

    double LimitingFunction(
        const double s,
        const double s_0,
        const double s_max,
        const double s_min = 0.0);

    double SmoothedLimitingFunction(
        const double s,
        const double s_0,
        const double s_max);

    // Smooth approximation of the max(s,0) function
    double SmoothedMaxFunction(const double s);

    // Smooth approximation of the min(s,0) function
    double SmoothedMinFunction(const double s);

    // https://es.wikipedia.org/wiki/Circunelipse_de_Steiner
    std::tuple<double, double, Matrix> CalculateTriangleMetricTensor(const Geometry<Node<3>> &rGeometry);

    // https://es.wikipedia.org/wiki/Circunelipse_de_Steiner --> 3D extension
    std::tuple<double, double, Matrix> CalculateTetrahedraMetricTensor(const Geometry<Node<3>> &rGeometry);

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
    ShockCapturingProcess& operator=(ShockCapturingProcess const& rOther);

    /// Copy constructor.
    ShockCapturingProcess(ShockCapturingProcess const& rOther);

    ///@}

}; // Class ShockCapturingProcess

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const ShockCapturingProcess& rThis);

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SHOCK_CAPTURING_UTILITIES_H_INCLUDED  defined
