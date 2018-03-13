// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____ 
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _ 
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
//

#if !defined(KRATOS_LEVELSET_METRICS_PROCESS)
#define KRATOS_LEVELSET_METRICS_PROCESS

// Project includes
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "includes/model_part.h"
#include "meshing_application.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    typedef ModelPart::NodesContainerType                        NodesArrayType;
    typedef ModelPart::ElementsContainerType                  ElementsArrayType;
    typedef ModelPart::ConditionsContainerType              ConditionsArrayType;
    typedef Node <3>                                                   NodeType;
    
///@}
///@name  Enum's
///@{
    
    #if !defined(INTERPOLATION_METRIC)
    #define INTERPOLATION_METRIC
        enum Interpolation {Constant = 0, Linear = 1, Exponential = 2};
    #endif
    
///@}
///@name  Functions
///@{
    
///@}
///@name Kratos Classes
///@{

//// This class is can be used to compute the metrics of the model part with a level set approach

template<unsigned int TDim>  
class ComputeLevelSetSolMetricProcess 
    : public Process
{
public:

    ///@name Type Definitions
    ///@{
    
    /// Pointer definition of ComputeLevelSetSolMetricProcess
    KRATOS_CLASS_POINTER_DEFINITION(ComputeLevelSetSolMetricProcess);
    
    ///@}
    ///@name Life Cycle
    ///@{
     
    // Constructor
    
    /**
     * This is the default constructor
     * @param rThisModelPart The model part to be computed
     * @param rVariableGradient The gradient variable
     * @param ThisParameters The input parameters
     */
    
    ComputeLevelSetSolMetricProcess(
        ModelPart& rThisModelPart,
        const Variable<array_1d<double,3>> rVariableGradient = DISTANCE_GRADIENT,
        Parameters ThisParameters = Parameters(R"({})")
        );
    
    /// Destructor.
    ~ComputeLevelSetSolMetricProcess() override = default; 
    
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
     * We initialize the metrics of the MMG sol using a level set approach
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
        return "ComputeLevelSetSolMetricProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ComputeLevelSetSolMetricProcess";
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

    ModelPart& mThisModelPart;                      // The model part to compute
    Variable<array_1d<double,3>> mVariableGradient; // The gradient variable
    double mMinSize;                                // The minimal size of the elements
    bool mEnforceCurrent;                           // With this we choose if we inforce the current nodal size (NODAL_H)
    double mAnisotropicRatio;                       // The minimal anisotropic ratio (0 < ratio < 1)
    double mBoundLayer;                             // The boundary layer limit Distance
    Interpolation mInterpolation;                   // The interpolation type
    
    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{
    
    /**
     * It calculates the tensor of the scalar, necessary to get the solution before remeshing
     * @param GradientValue The gradient of the scalar to remesh
     * @param Ratio The alpha parameter used to remesh
     * @param ElementSize The minimum size of the elements
     * @return The metric tensor
     */
        
    Vector ComputeLevelSetMetricTensor(
        const array_1d<double, 3>& GradientValue,
        const double Ratio,
        const double ElementSize
        );

    
    /**
     * This converts the interpolation string to an enum
     * @param str The string that you want to comvert in the equivalent enum
     * @return Interpolation: The equivalent enum (this requires less memmory than a std::string)
     */
        
    Interpolation ConvertInter(const std::string& str);
        
    /**
     * This calculates the anisotropic ratio
     * @param Distance Distance parameter
     * @param AnisotropicRatio The anisotropic ratio
     * @param BoundLayer The boundary layer limit
     * @param rInterpolation The type of interpolation
     */
    
    double CalculateAnisotropicRatio(
        const double Distance,
        const double AnisotropicRatio,
        const double BoundLayer,
        const Interpolation& rInterpolation
        );
    
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
    ComputeLevelSetSolMetricProcess& operator=(ComputeLevelSetSolMetricProcess const& rOther) 
    {
        return *this;
    };

    /// Copy constructor.
    //ComputeLevelSetSolMetricProcess(ComputeLevelSetSolMetricProcess const& rOther);

    ///@}
};// class ComputeLevelSetSolMetricProcess
///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
template<unsigned int TDim> 
inline std::istream& operator >> (std::istream& rIStream,
                                  ComputeLevelSetSolMetricProcess<TDim>& rThis);

/// output stream function
template<unsigned int TDim> 
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ComputeLevelSetSolMetricProcess<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

};// namespace Kratos.
#endif /* KRATOS_LEVELSET_METRICS_PROCESS defined */
