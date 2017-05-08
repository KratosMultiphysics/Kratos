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

#if !defined(KRATOS_ERROR_METRICS_PROCESS)
#define KRATOS_ERROR_METRICS_PROCESS

// Project includes
#include "utilities/math_utils.h"
#include "custom_utilities/metrics_math_utils.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "meshing_application.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    typedef ModelPart::NodesContainerType                                     NodesArrayType;
    typedef ModelPart::ElementsContainerType                               ElementsArrayType;
    typedef ModelPart::ConditionsContainerType                           ConditionsArrayType;
    typedef Node <3>                                                                NodeType;
    
///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{
    
///@}
///@name Kratos Classes
///@{

//// This class is can be used to compute the Metrics of the model part with an Hessian approach

template<unsigned int TDim>  
class ComputeErrorSolMetricProcess
    : public Process
{
public:

    ///@name Type Definitions
    ///@{
    
    /// Pointer definition of ComputeErrorSolMetricProcess
    KRATOS_CLASS_POINTER_DEFINITION(ComputeErrorSolMetricProcess);
    
    ///@}
    ///@name Life Cycle
    ///@{
     
    // Constructor
    
    /**
     * This is the default constructor
     * @param rThisModelPart: The model part to be computed
     * @param ThisParameters: The input parameters
     */
    
    ComputeErrorSolMetricProcess(
        ModelPart& rThisModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        )
        :mThisModelPart(rThisModelPart)
    {               
        Parameters DefaultParameters = Parameters(R"(
        {
            "minimal_size"                        : 0.1,
            "maximal_size"                        : 10.0, 
            "enforce_current"                     : true, 
            "error_strategy_parameters": 
            {
            }
        })" );
        ThisParameters.ValidateAndAssignDefaults(DefaultParameters);
         
        mMinSize = ThisParameters["minimal_size"].GetDouble();
        mMaxSize = ThisParameters["maximal_size"].GetDouble();
        mEnforceCurrent = ThisParameters["enforce_current"].GetBool();
    }
    
    /// Destructor.
    virtual ~ComputeErrorSolMetricProcess() {}
    
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
     * We initialize the Metrics of the MMG sol using the Hessian Metric matrix approach
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

            const double NodalH = itNode->FastGetSolutionStepValue(NODAL_H, 0);            
            
            double ElementMinSize = mMinSize;
            if ((ElementMinSize > NodalH) && (mEnforceCurrent == true))
            {
                ElementMinSize = NodalH;
            }
            double ElementMaxSize = mMaxSize;
            if ((ElementMaxSize > NodalH) && (mEnforceCurrent == true))
            {
                ElementMaxSize = NodalH;
            }
            
            // We compute the Metric
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
            
            const double NodalError = itNode->GetValue(NODAL_ERROR);
            const double NormMetric = norm_2(Metric);
            if (NormMetric > 0.0) // NOTE: This means we combine differents metrics, at the same time means that the metric should be reseted each time
            {
                const Vector OldMetric = itNode->GetValue(MMG_METRIC);
                const Vector NewMetric = ComputeErrorMetricTensor(NodalH, NodalError, ElementMinSize, ElementMaxSize);    
                
                Metric = MetricsMathUtils<TDim>::IntersectMetrics(OldMetric, NewMetric);
            }
            else
            {
                Metric = ComputeErrorMetricTensor(NodalH, NodalError, ElementMinSize, ElementMaxSize);    
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
        return "ComputeErrorSolMetricProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ComputeErrorSolMetricProcess";
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
    
    ModelPart& mThisModelPart;               // The model part to compute
    double mMinSize;                         // The minimal size of the elements
    double mMaxSize;                         // The maximal size of the elements
    bool mEnforceCurrent;                    // With this we choose if we inforce the current nodal size (NODAL_H)
    double mTolerance;                       // The error tolerance considered
    
    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * This function computes the metric tensor using as reference the nodal error
     * @param OldNodalH: The old size of the mesh in that node
     * @param NodalError: The weighted on the node  
     * @param ElementMinSize: This way we can impose as minimum as the previous size if we desire
     * @param ElementMaxSize: This way we can impose as maximum as the previous size if we desire
     */
        
    Vector ComputeErrorMetricTensor(
        const double& OldNodalH,
        const double& NodalError, 
        const double& ElementMinSize,
        const double& ElementMaxSize
        )
    {        
        // Calculating Metric parameters
        const double MinRatio = 1.0/(ElementMinSize * ElementMinSize);
//         const double MinRatio = 1.0/(mMinSize * mMinSize);
        const double MaxRatio = 1.0/(ElementMaxSize * ElementMaxSize);
//         const double MaxRatio = 1.0/(mMaxSize * mMaxSize);
        
        const double AverageNodalError = mThisModelPart.GetProcessInfo()[AVERAGE_NODAL_ERROR];
        
        boost::numeric::ublas::bounded_matrix<double, TDim, TDim> MetricMatrix = ZeroMatrix(TDim, TDim);
        
        const double NewNodalH = OldNodalH * AverageNodalError/NodalError;
        
        const double NewNodalRatio = NewNodalH > ElementMinSize ? MinRatio : NewNodalH < ElementMaxSize ? MaxRatio : 1.0/(NewNodalH * NewNodalH);
        
        // Right now just considering isotropic mesh
        for (unsigned int iDim = 0; iDim < TDim; iDim++)
        {
            MetricMatrix(iDim, iDim) = NewNodalRatio;
        }
        
        // Finally we transform to a vector
        const Vector Metric = MetricsMathUtils<TDim>::TensorToVector(MetricMatrix);
        
        return Metric;
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
    ComputeErrorSolMetricProcess& operator=(ComputeErrorSolMetricProcess const& rOther);

    /// Copy constructor.
    //ComputeErrorSolMetricProcess(ComputeErrorSolMetricProcess const& rOther);

    ///@}
};// class ComputeErrorSolMetricProcess
///@}


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
template<unsigned int TDim, class TVarType> 
inline std::istream& operator >> (std::istream& rIStream,
                                  ComputeErrorSolMetricProcess<TDim>& rThis);

/// output stream function
template<unsigned int TDim, class TVarType> 
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ComputeErrorSolMetricProcess<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

};// namespace Kratos.
#endif /* KRATOS_ERROR_METRICS_PROCESS defined */
