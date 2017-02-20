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
     * @param rVariableGradient: The gradient variable to compute
     * @param rMinSize: The min size of element
     * @param rAnisRatio: The anisotropic ratio
     * @param rBoundLayer: The boundary layer limit
     * @param rInterpolation: The interpolation type
     */
    
    ComputeLevelSetSolMetricProcess(
        ModelPart& rThisModelPart,
        const double rMinSize,
        const bool rEnforceCurrent = true,
        const Variable<array_1d<double,3>> rVariableGradient = DISTANCE_GRADIENT,
        const double rAnisRatio = 1.0,
        const double rBoundLayer =  1.0,
        const std::string rInterpolation = "Linear"
        )
        :mThisModelPart(rThisModelPart),
        mVariableGradient(rVariableGradient),
        mMinSize(rMinSize),
        mEnforceCurrent(rEnforceCurrent),
        mAnisRatio(rAnisRatio),
        mBoundLayer(rBoundLayer)
    {   
        mInterpolation = ConvertInter(rInterpolation);
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
            
            const double distance = itNode->FastGetSolutionStepValue(DISTANCE, 0);
            array_1d<double, 3> gradient_value = itNode->FastGetSolutionStepValue(mVariableGradient, 0);
            
            const double ratio = CalculateAnisotropicRatio(distance, mAnisRatio, mBoundLayer, mInterpolation);
            
            // For postprocess pourposes
            double& anisotropic_ratio = itNode->FastGetSolutionStepValue(ANISOTROPIC_RATIO, 0); 
            anisotropic_ratio = ratio;
            

            double element_size = mMinSize;
            const double nodal_h = itNode->FastGetSolutionStepValue(NODAL_H, 0);
            if (((element_size > nodal_h) && (mEnforceCurrent == true)) || (std::abs(distance) > mBoundLayer))
            {
                element_size = nodal_h;
            }
            
            const double tolerance = 1.0e-12;
            const double norm = norm_2(gradient_value);
            if (norm > tolerance)
            {
                gradient_value /= norm;
            }
            
            // We compute the metric
            #ifdef KRATOS_DEBUG 
            if( itNode->Has(MMG_METRIC) == false) 
            {
                KRATOS_ERROR <<  " MMG_METRIC not defined for node " << itNode->Id();
            }
            #endif     
            Vector& metric = itNode->GetValue(MMG_METRIC);
            
            #ifdef KRATOS_DEBUG 
            if(metric.size() != TDim * 3 - 3) 
            {
                KRATOS_ERROR << "Wrong size of vector MMG_METRIC found for node " << itNode->Id() << " size is " << metric.size() << " expected size was " << TDim * 3 - 3;
            }
            #endif
            
            const double normmetric = norm_2(metric);
            if (normmetric > 0.0) // NOTE: This means we combine differents metrics, at the same time means that the metric should be reseted each time
            {
                const Vector old_metric = itNode->GetValue(MMG_METRIC);
                const Vector new_metric = ComputeLevelSetMetricTensor(gradient_value, ratio, element_size);
                
                metric = MetricsMathUtils<TDim>::IntersectMetrics(old_metric, new_metric);
            }
            else
            {
                metric = ComputeLevelSetMetricTensor(gradient_value, ratio, element_size);
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

    /// Print object's data.
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
    double mBoundLayer;                             // The boundary layer limit distance
    Interpolation mInterpolation;                   // The interpolation type
    
    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{
    
    /**
     * It calculates the tensor of the scalar, necessary to get the solution before remeshing
     * @param gradient_value: The gradient of the scalar to remesh
     * @param ratio: The alpha parameter used to remesh
     * @param element_size: The minimum size of the elements
     * @param node_id: The id of the node
     */
        
    Vector ComputeLevelSetMetricTensor(
        const array_1d<double, 3>& gradient_value,
        const double& ratio,
        const double& element_size
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
     * @param distance: Distance parameter
     */
    
    double CalculateAnisotropicRatio(
        const double& distance,
        const double& rAnisRatio,
        const double& rBoundLayer,
        const Interpolation& rInterpolation
        )
    {
        const double tolerance = 1.0e-12;
        double ratio = 1.0; // NOTE: Isotropic mesh
        if (rAnisRatio < 1.0)
        {                           
            if (std::abs(distance) <= rBoundLayer)
            {
                if (rInterpolation == Constant)
                {
                    ratio = rAnisRatio;
                }
                else if (rInterpolation == Linear)
                {
                    ratio = rAnisRatio + (std::abs(distance)/rBoundLayer) * (1.0 - rAnisRatio);
                }
                else if (rInterpolation == Exponential)
                {
                    ratio = - std::log(std::abs(distance)/rBoundLayer) * rAnisRatio + tolerance;
                    if (ratio > 1.0)
                    {
                        ratio = 1.0;
                    }
                }
            }
        }
        
        return ratio;
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
        const array_1d<double, 3>& gradient_value,
        const double& ratio,
        const double& element_size
    )
    {
        Vector metric;
        metric.resize(3, false);
        
        const double coeff0 = 1.0/(element_size * element_size);
        const double coeff1 = coeff0/(ratio * ratio);
        
        const double v0v0 = gradient_value[0]*gradient_value[0];
        const double v0v1 = gradient_value[0]*gradient_value[1];
        const double v1v1 = gradient_value[1]*gradient_value[1];
        
        metric[0] = coeff0*(1.0 - v0v0) + coeff1*v0v0;
        metric[1] = coeff0*(    - v0v1) + coeff1*v0v1;  
        metric[2] = coeff0*(1.0 - v1v1) + coeff1*v1v1;
        
        return metric;
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    template<>  
    Vector ComputeLevelSetSolMetricProcess<3>::ComputeLevelSetMetricTensor(
        const array_1d<double, 3>& gradient_value,
        const double& ratio,
        const double& element_size
    )
    {
        Vector metric;
        metric.resize(6, false);
        
        const double coeff0 = 1.0/(element_size * element_size);
        const double coeff1 = coeff0/(ratio * ratio);
        
        const double v0v0 = gradient_value[0]*gradient_value[0];
        const double v0v1 = gradient_value[0]*gradient_value[1];
        const double v0v2 = gradient_value[0]*gradient_value[2];
        const double v1v1 = gradient_value[1]*gradient_value[1];
        const double v1v2 = gradient_value[1]*gradient_value[2];
        const double v2v2 = gradient_value[2]*gradient_value[2];
        
        metric[0] = coeff0*(1.0 - v0v0) + coeff1*v0v0;
        metric[1] = coeff0*(    - v0v1) + coeff1*v0v1; 
        metric[2] = coeff0*(    - v0v2) + coeff1*v0v2; 
        metric[3] = coeff0*(1.0 - v1v1) + coeff1*v1v1; 
        metric[4] = coeff0*(    - v1v2) + coeff1*v1v2; 
        metric[5] = coeff0*(1.0 - v2v2) + coeff1*v2v2;

        return metric;
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
