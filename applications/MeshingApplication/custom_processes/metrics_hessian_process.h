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
#include "meshing_application.h"
#include "processes/process.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    /// The size type definition
    typedef std::size_t SizeType;

    /// Component variable type
    typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > >  ComponentType;

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
 * @class ComputeHessianSolMetricProcess
 * @ingroup MeshingApplication
 * @brief This class is can be used to compute the metrics of the model part with an Hessian approach
 * @details References:
 *         [1] F. Alauzet, Metric-based anisotropic mesh adaptation, CEA-EDF-INRIA schools: Numerical Analysis Summer School. CEA, Cadarache, France
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

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor

    /**
     * @brief This is the default constructor (double)
     * @param rThisModelPart The model part to be computed
     * @param rVariable The variable to compute
     * @param ThisParameters The input parameters
     */
    ComputeHessianSolMetricProcess(
        ModelPart& rThisModelPart,
        Variable<double>& rVariable,
        Parameters ThisParameters = Parameters(R"({})")
        );

    /**
     * @brief This is the default constructor (component)
     * @param rThisModelPart The model part to be computed
     * @param rVariable The variable to compute
     * @param ThisParameters The input parameters
     */
    ComputeHessianSolMetricProcess(
        ModelPart& rThisModelPart,
        ComponentType& rVariable,
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

    ModelPart& mThisModelPart;                                  /// The model part to compute
    std::vector<Variable<double>*> mrOriginVariableDoubleList;  /// The scalar variable list to compute
    std::vector<ComponentType*> mrOriginVariableComponentsList; /// The scalar variable list to compute (components)

    // TODO: Replace for Parameters
    std::string mRatioReferenceVariable = "DISTANCE";          /// Variable used to compute the anisotropic ratio
    double mMinSize;                                           /// The minimal size of the elements
    double mMaxSize;                                           /// The maximal size of the elements
    bool mEnforceCurrent;                                      /// With this we choose if we inforce the current nodal size (NODAL_H)
    bool mEstimateInterpError;                                 /// If the error of interpolation will be estimated
    double mInterpError;                                       /// The error of interpolation allowed
    double mMeshConstant;                                      /// The mesh constant to remesh (depends of the element type)
    double mAnisotropicRatio;                                  /// The minimal anisotropic ratio (0 < ratio < 1)
    double mBoundLayer;                                        /// The boundary layer limit distance
    Interpolation mInterpolation;                              /// The interpolation type

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
     * @param AnisotropicRatio The anisotropic ratio
     * @param ElementMinSize The min size of element
     * @param ElementMaxSize The maximal size of the elements
     * @param NodalH The size of the local node
     */
    template<SizeType TDim>
    array_1d<double, 3 * (TDim - 1)> ComputeHessianMetricTensor(
        const Vector& rHessian,
        const double AnisotropicRatio,
        const double ElementMinSize, // This way we can impose as minimum as the previous size if we desire
        const double ElementMaxSize, // This way we can impose as maximum as the previous size if we desire
        const double NodalH
        );

    /**
     * @brief This calculates the auxiliar hessian needed for the Metric
     */
    void CalculateAuxiliarHessian();

    /**
     * @brief This converts the interpolation string to an enum
     * @param Str The string that you want to comvert in the equivalent enum
     * @return Interpolation: The equivalent enum (this requires less memmory than a std::string)
     */
    Interpolation ConvertInter(const std::string& Str)
    {
        if(Str == "Constant" || Str == "CONSTANT")
            return Interpolation::CONSTANT;
        else if(Str == "Linear" || Str == "LINEAR")
            return Interpolation::LINEAR;
        else if(Str == "Exponential" || Str == "EXPONENTIAL")
            return Interpolation::EXPONENTIAL;
        else
            return Interpolation::LINEAR;
    }

    /**
     * @brief This calculates the anisotropic ratio
     * @param Distance Distance parameter
     * @param AnisotropicRatio The anisotropic ratio
     * @param BoundLayer The boundary layer limit
     * @param rInterpolation The type of interpolation
     */
    double CalculateAnisotropicRatio(
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
    Parameters GetDefaultParameters() const;

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
