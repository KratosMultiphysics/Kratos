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
#include "utilities/math_utils.h"
#include "custom_utilities/metrics_math_utils.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"
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
     * @param rThisModelPart: The model part to be computed
     * @param ThisParameters: The input parameters
     */
    
    ComputeLevelSetSolMetricProcess(
        ModelPart& rThisModelPart,
        const Variable<array_1d<double,3>> rVariableGradient = DISTANCE_GRADIENT,
        Parameters ThisParameters = Parameters(R"({})")
        )
        :mThisModelPart(rThisModelPart),
        mVariableGradient(rVariableGradient),
        mMinSize(ThisParameters["minimal_size"].GetDouble()),
        mEnforceCurrent(ThisParameters["enforce_current"].GetBool())
    {   
        Parameters DefaultParameters = Parameters(R"(
        {
            "minimal_size"                         : 0.1, 
            "enforce_current"                      : true, 
            "anisotropy_remeshing"                 : true, 
            "anisotropy_parameters": 
            {
                "hmin_over_hmax_anisotropic_ratio"      : 1.0, 
                "boundary_layer_max_distance"           : 1.0, 
                "interpolation"                         : "Linear"
            }
        })" );
        ThisParameters.ValidateAndAssignDefaults(DefaultParameters);
        
        // In case we have isotropic remeshing (default values)
        if (ThisParameters["anisotropy_remeshing"].GetBool() == false)
        {
            mAnisRatio = DefaultParameters["anisotropy_parameters"]["hmin_over_hmax_anisotropic_ratio"].GetDouble();
            mBoundLayer = DefaultParameters["anisotropy_parameters"]["boundary_layer_max_distance"].GetDouble();
            mInterpolation = ConvertInter(DefaultParameters["anisotropy_parameters"]["interpolation"].GetString());
        }
        else
        {
            mAnisRatio = ThisParameters["anisotropy_parameters"]["hmin_over_hmax_anisotropic_ratio"].GetDouble();
            mBoundLayer = ThisParameters["anisotropy_parameters"]["boundary_layer_max_distance"].GetDouble();
            mInterpolation = ConvertInter(ThisParameters["anisotropy_parameters"]["interpolation"].GetString());
        }
    }
    
    /// Destructor.
    virtual ~ComputeLevelSetSolMetricProcess() {}
    
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
    
    virtual void Execute()
    {
        // Iterate in the nodes
        NodesArrayType& pNode = mThisModelPart.Nodes();
        int numNodes = pNode.end() - pNode.begin();
        
        #pragma omp parallel for 
        for(int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            
            if ( itNode->SolutionStepsDataHas( mVariableGradient ) == false )
            {
                KRATOS_ERROR << "Missing gradient variable on node " << itNode->Id() << std::endl;
            }
            
            const double Distance = itNode->FastGetSolutionStepValue(DISTANCE, 0);
            array_1d<double, 3> GradientValue = itNode->FastGetSolutionStepValue(mVariableGradient, 0);
            
            const double Ratio = CalculateAnisotropicRatio(Distance, mAnisRatio, mBoundLayer, mInterpolation);
            
            // For postprocess pourposes
            itNode->SetValue(ANISOTROPIC_RATIO, Ratio); 
            
            double ElementSize = mMinSize;
            const double NodalH = itNode->FastGetSolutionStepValue(NODAL_H, 0);
            if (((ElementSize > NodalH) && (mEnforceCurrent == true)) || (std::abs(Distance) > mBoundLayer))
            {
                ElementSize = NodalH;
            }
            
            const double Tolerance = 1.0e-12;
            const double NormGradientValue = norm_2(GradientValue);
            if (NormGradientValue > Tolerance)
            {
                GradientValue /= NormGradientValue;
            }
            
            // We compute the metric
            #ifdef KRATOS_DEBUG 
            if( itNode->Has(MMG_METRIC) == false) 
            {
                KRATOS_ERROR <<  " MMG_METRIC not defined for node " << itNode->Id();
            }
            #endif     
            Vector& Metric = itNode->GetValue(MMG_METRIC);
            
            #ifdef KRATOS_DEBUG 
            if(Metric.size() != TDim * 3 - 3) 
            {
                KRATOS_ERROR << "Wrong size of vector MMG_METRIC found for node " << itNode->Id() << " size is " << Metric.size() << " expected size was " << TDim * 3 - 3;
            }
            #endif
            
            const double NormMetric = norm_2(Metric);
            if (NormMetric > 0.0) // NOTE: This means we combine differents metrics, at the same time means that the metric should be reseted each time
            {
                const Vector OldMetric = itNode->GetValue(MMG_METRIC);
                const Vector NewMetric = ComputeLevelSetMetricTensor(GradientValue, Ratio, ElementSize);
                
                Metric = MetricsMathUtils<TDim>::IntersectMetrics(OldMetric, NewMetric);
            }
            else
            {
                Metric = ComputeLevelSetMetricTensor(GradientValue, Ratio, ElementSize);
            }
        }
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
    
    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "ComputeLevelSetSolMetricProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ComputeLevelSetSolMetricProcess";
    }

    /// Print object"s data.
    virtual void PrintData(std::ostream& rOStream) const
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
    double mAnisRatio;                              // The minimal anisotropic ratio (0 < ratio < 1)
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
     * @param GradientValue: The gradient of the scalar to remesh
     * @param Ratio: The alpha parameter used to remesh
     * @param ElementSize: The minimum size of the elements
     * @param node_id: The id of the node
     */
        
    Vector ComputeLevelSetMetricTensor(
        const array_1d<double, 3>& GradientValue,
        const double& Ratio,
        const double& ElementSize
    );

    
    /**
     * This converts the interpolation string to an enum
     * @param str: The string that you want to comvert in the equivalent enum
     * @return Interpolation: The equivalent enum (this requires less memmory than a std::string)
     */
        
    Interpolation ConvertInter(const std::string& str)
    {
        if(str == "Constant") 
        {
            return Constant;
        }
        else if(str == "Linear") 
        {
            return Linear;
        }
        else if(str == "Exponential") 
        {
            return Exponential;
        }
        else
        {
            return Linear;
        }
    }
        
    /**
     * This calculates the anisotropic ratio
     * @param Distance: Distance parameter
     */
    
    double CalculateAnisotropicRatio(
        const double& Distance,
        const double& rAnisRatio,
        const double& rBoundLayer,
        const Interpolation& rInterpolation
        )
    {
        const double Tolerance = 1.0e-12;
        double Ratio = 1.0; // NOTE: Isotropic mesh
        if (rAnisRatio < 1.0)
        {                           
            if (std::abs(Distance) <= rBoundLayer)
            {
                if (rInterpolation == Constant)
                {
                    Ratio = rAnisRatio;
                }
                else if (rInterpolation == Linear)
                {
                    Ratio = rAnisRatio + (std::abs(Distance)/rBoundLayer) * (1.0 - rAnisRatio);
                }
                else if (rInterpolation == Exponential)
                {
                    Ratio = - std::log(std::abs(Distance)/rBoundLayer) * rAnisRatio + Tolerance;
                    if (Ratio > 1.0)
                    {
                        Ratio = 1.0;
                    }
                }
            }
        }
        
        return Ratio;
    }
    
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
    ComputeLevelSetSolMetricProcess& operator=(ComputeLevelSetSolMetricProcess const& rOther) {};

    /// Copy constructor.
    //ComputeLevelSetSolMetricProcess(ComputeLevelSetSolMetricProcess const& rOther);

    ///@}
};// class ComputeLevelSetSolMetricProcess
///@}

///@name Explicit Specializations
///@{

    template<>  
    Vector ComputeLevelSetSolMetricProcess<2>::ComputeLevelSetMetricTensor(
        const array_1d<double, 3>& GradientValue,
        const double& Ratio,
        const double& ElementSize
    )
    {
        Vector Metric;
        Metric.resize(3, false);
        
        const double Coeff0 = 1.0/(ElementSize * ElementSize);
        const double Coeff1 = Coeff0/(Ratio * Ratio);
        
        const double v0v0 = GradientValue[0]*GradientValue[0];
        const double v0v1 = GradientValue[0]*GradientValue[1];
        const double v1v1 = GradientValue[1]*GradientValue[1];
        
        Metric[0] = Coeff0*(1.0 - v0v0) + Coeff1*v0v0;
        Metric[1] = Coeff0*(    - v0v1) + Coeff1*v0v1;  
        Metric[2] = Coeff0*(1.0 - v1v1) + Coeff1*v1v1;
        
        return Metric;
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    template<>  
    Vector ComputeLevelSetSolMetricProcess<3>::ComputeLevelSetMetricTensor(
        const array_1d<double, 3>& GradientValue,
        const double& Ratio,
        const double& ElementSize
    )
    {
        Vector Metric;
        Metric.resize(6, false);
        
        const double Coeff0 = 1.0/(ElementSize * ElementSize);
        const double Coeff1 = Coeff0/(Ratio * Ratio);
        
        const double v0v0 = GradientValue[0]*GradientValue[0];
        const double v0v1 = GradientValue[0]*GradientValue[1];
        const double v0v2 = GradientValue[0]*GradientValue[2];
        const double v1v1 = GradientValue[1]*GradientValue[1];
        const double v1v2 = GradientValue[1]*GradientValue[2];
        const double v2v2 = GradientValue[2]*GradientValue[2];
        
        Metric[0] = Coeff0*(1.0 - v0v0) + Coeff1*v0v0;
        Metric[1] = Coeff0*(    - v0v1) + Coeff1*v0v1; 
        Metric[2] = Coeff0*(    - v0v2) + Coeff1*v0v2; 
        Metric[3] = Coeff0*(1.0 - v1v1) + Coeff1*v1v1; 
        Metric[4] = Coeff0*(    - v1v2) + Coeff1*v1v2; 
        Metric[5] = Coeff0*(1.0 - v2v2) + Coeff1*v2v2;

        return Metric;
    }
    
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
