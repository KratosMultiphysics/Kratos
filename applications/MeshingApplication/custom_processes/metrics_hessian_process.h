// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_HESSIAN_METRICS_PROCESS)
#define KRATOS_HESSIAN_METRICS_PROCESS

// Project includes
#include "meshing_application_variables.h"
#include "processes/process.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    /// The size type definition
    typedef std::size_t SizeType;

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
     * @struct AuxiliarHessianComputationVariables
     * @ingroup MeshingApplication
     * @brief This is an auxiliar struct to store remeshing variables
     * @author Vicente Mataix Ferrandiz
     */
    struct AuxiliarHessianComputationVariables
    {
        AuxiliarHessianComputationVariables(
            const double AnisotropicRatio,
            const double ElementMinSize,
            const double ElementMaxSize,
            const double NodalH,
            const bool EstimateInterpolationError,
            const double InterpolationError,
            const double MeshDependentConstant,
            const bool AnisotropyRemeshing,
            const bool EnforceAnisotropyRelativeVariable
        ) : mAnisotropicRatio(AnisotropicRatio),
            mElementMinSize(ElementMinSize),
            mElementMaxSize(ElementMaxSize),
            mNodalH(NodalH),
            mEstimateInterpolationError(EstimateInterpolationError),
            mInterpolationError(InterpolationError),
            mMeshDependentConstant(MeshDependentConstant),
            mAnisotropyRemeshing(AnisotropyRemeshing),
            mEnforceAnisotropyRelativeVariable(EnforceAnisotropyRelativeVariable)
        {

        };

        double mAnisotropicRatio;                /// The anisotropic ratio
        double mElementMinSize;                  /// The min size of element. This way we can impose as minimum as the previous size if we desire
        double mElementMaxSize;                  /// The maximal size of the elements. This way we can impose as maximum as the previous size if we desire
        double mNodalH;                          /// The size of the local node
        const bool mEstimateInterpolationError;        /// If the error of interpolation will be estimated
        const double mInterpolationError;              /// The error of interpolation allowed
        const double mMeshDependentConstant;           /// The mesh constant to remesh (depends of the element type)
        const bool mAnisotropyRemeshing;               /// If we consider anisotropic remeshing
        const bool mEnforceAnisotropyRelativeVariable; /// If we enforce a certain anisotropy relative toa  variable
    };

/**
 * @class ComputeHessianSolMetricProcess
 * @ingroup MeshingApplication
 * @brief This class is can be used to compute the metrics of the model part with an Hessian approach
 * @details References:
 *         [1] P.J. Frey, F. Alauzet; Anisotropic mesh adaptation for CFD computations; Comput. Methods Appl. Mech. Engrg. 194 (2005) 5068–5082
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(MESHING_APPLICATION) ComputeHessianSolMetricProcess
    : public Process
{
public:

    ///@name Type Definitions
    ///@{

    /// Node definition
    typedef Node <3>                                                                NodeType;

    /// Containers definitions
    typedef ModelPart::NodesContainerType                                     NodesArrayType;
    typedef ModelPart::ElementsContainerType                               ElementsArrayType;
    typedef ModelPart::ConditionsContainerType                           ConditionsArrayType;

    /// The index type definition
    typedef std::size_t                                                            IndexType;

    /// Pointer definition of ComputeHessianSolMetricProcess
    KRATOS_CLASS_POINTER_DEFINITION(ComputeHessianSolMetricProcess);

    ///@}
    ///@name  Enum's
    ///@{

    /**
     * @brief This enums allows to differentiate the interpolation types
     */
    enum class Interpolation {CONSTANT = 0, LINEAR = 1, EXPONENTIAL = 2};

    /**
     * @brief This enums allows to differentiate the normalization types
     */
    enum class Normalization {CONSTANT = 0, VALUE = 1, NORM_GRADIENT = 2};

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor

    /**
     * @brief This is the default constructor (pure parameters)
     * @param rThisModelPart The model part to be computed
     * @param ThisParameters The input parameters
     */
    ComputeHessianSolMetricProcess(
        ModelPart& rThisModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        );

    /**
     * @brief This is the default constructor
     * @param rThisModelPart The model part to be computed
     * @param rVariable The variable to compute
     * @param ThisParameters The input parameters
     */
    ComputeHessianSolMetricProcess(
        ModelPart& rThisModelPart,
        Variable<double>& rVariable,
        Parameters ThisParameters = Parameters(R"({})")
        );

    /// Destructor.
    ~ComputeHessianSolMetricProcess() override = default;

    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief We initialize the metrics of the MMG sol using the Hessian metric matrix approach
     */
    void Execute() override;

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
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
        return "ComputeHessianSolMetricProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ComputeHessianSolMetricProcess";
    }

    /// Print object"s data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

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
    ///@name Private static Member Variables
    ///@{

    ///@}
    ///@name Private member Variables
    ///@{

    ModelPart& mrModelPart;                                           /// The model part to compute

    bool mNonHistoricalVariable = false;                              /// If the variable is non-historical
    const Variable<double>* mpOriginVariable;                         /// The scalar variable list to compute
    const Variable<double>* mpRatioReferenceVariable;                 /// Variable used to compute the anisotropic ratio

    Parameters mThisParameters;                                       /// Here configurations are stored
    Interpolation mInterpolation;                                     /// The interpolation type

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This function is used to compute the Hessian Metric tensor
     * @details Note that when using the Hessian, more than one Metric can be defined simultaneously, so in consecuence we need to define the elipsoid which defines the volume of maximal intersection
     * @param Hessian The hessian tensor condensed already computed
     * @param rAuxiliarHessianComputationVariables Struct containing several variables
     */
    template<SizeType TDim>
    static array_1d<double, 3 * (TDim - 1)> ComputeHessianMetricTensor(
        const Vector& rHessian,
        const AuxiliarHessianComputationVariables& rAuxiliarHessianComputationVariables
        );

    /**
     * @brief This calculates the auxiliar hessian needed for the Metric
     */
    void CalculateAuxiliarHessian();

    /**
     * @brief This converts the interpolation string to an enum
     * @param Str The string that you want to convert in the equivalent enum
     * @return Interpolation: The equivalent enum (this requires less memmory and is eassier to compare than a std::string)
     */
    Interpolation ConvertInter(const std::string& Str)
    {
        if(Str == "Constant" || Str == "CONSTANT" || Str == "constant")
            return Interpolation::CONSTANT;
        else if(Str == "Linear" || Str == "LINEAR"  || Str == "linear")
            return Interpolation::LINEAR;
        else if(Str == "Exponential" || Str == "EXPONENTIAL"  || Str == "exponential")
            return Interpolation::EXPONENTIAL;
        else
            return Interpolation::LINEAR;
    }

    /**
     * @brief This converts the normalization string to an enum
     * @param Str The string that you want to convert in the equivalent enum
     * @return Normalization: The equivalent enum (this requires less memmory and is eassier to compare than a std::string)
     */
    Normalization ConvertNormalization(const std::string& Str)
    {
        if(Str == "Constant" || Str == "CONSTANT" || Str == "constant")
            return Normalization::CONSTANT;
        else if(Str == "Value" || Str == "VALUE" || Str == "value")
            return Normalization::VALUE;
        else if(Str == "Norm_Gradient" || Str == "NORM_GRADIENT" || Str == "norm_gradient")
            return Normalization::NORM_GRADIENT;
        else
            return Normalization::CONSTANT;
    }

    /**
     * @brief This calculates the anisotropic ratio
     * @param Distance Distance parameter
     * @param AnisotropicRatio The anisotropic ratio
     * @param BoundLayer The boundary layer limit
     * @param rInterpolation The type of interpolation
     */
    static double CalculateAnisotropicRatio(
        const double Distance,
        const double AnisotropicRatio,
        const double BoundLayer,
        const Interpolation rInterpolation
        );

    /**
     * @brief This method is the responsible to compute the metric of the problem
     */
    template<SizeType TDim>
    void CalculateMetric();

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    void InitializeVariables(Parameters ThisParameters);

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Private LifeCycle
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    ComputeHessianSolMetricProcess& operator=(ComputeHessianSolMetricProcess const& rOther);

    /// Copy constructor.
    //ComputeHessianSolMetricProcess(ComputeHessianSolMetricProcess const& rOther);

    ///@}
};// class ComputeHessianSolMetricProcess
///@}


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ComputeHessianSolMetricProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ComputeHessianSolMetricProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

};// namespace Kratos.
#endif /* KRATOS_HESSIAN_METRICS_PROCESS defined */
