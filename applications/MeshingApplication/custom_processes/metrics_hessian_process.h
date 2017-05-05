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

#if !defined(KRATOS_HESSIAN_METRICS_PROCESS)
#define KRATOS_HESSIAN_METRICS_PROCESS

// Project includes
#include "utilities/math_utils.h"
#include "custom_utilities/metrics_math_utils.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"
#include "meshing_application.h"
#include "processes/compute_nodal_gradient_process.h" // TODO: Not prism or quadrilaterals implemented yet

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

//// This class is can be used to compute the metrics of the model part with an Hessian approach

template<unsigned int TDim, class TVarType>  
class ComputeHessianSolMetricProcess
    : public Process
{
public:

    ///@name Type Definitions
    ///@{
    
    /// Pointer definition of ComputeHessianSolMetricProcess
    KRATOS_CLASS_POINTER_DEFINITION(ComputeHessianSolMetricProcess);
    
    ///@}
    ///@name Life Cycle
    ///@{
     
    // Constructor
    
    /**
     * This is the default constructor
     * @param rThisModelPart: The model part to be computed
     * @param ThisParameters: The input parameters
     */
    
    ComputeHessianSolMetricProcess(
        ModelPart& rThisModelPart,
        TVarType& rVariable,
        Parameters ThisParameters = Parameters(R"({})")
        )
        :mThisModelPart(rThisModelPart),
        mVariable(rVariable)
    {               
        Parameters DefaultParameters = Parameters(R"(
        {
            "minimal_size"                        : 0.1,
            "maximal_size"                        : 10.0, 
            "enforce_current"                     : true, 
            "hessian_strategy_parameters": 
            { 
                "interpolation_error"                  : 1.0e-6, 
                "mesh_dependent_constant"              : 0.28125
            }, 
            "anisotropy_remeshing"                : true, 
            "anisotropy_parameters":
            {
                "hmin_over_hmax_anisotropic_ratio"     : 1.0, 
                "boundary_layer_max_distance"          : 1.0, 
                "interpolation"                        : "Linear"
            }
        })" );
        ThisParameters.ValidateAndAssignDefaults(DefaultParameters);
         
        mMinSize = ThisParameters["minimal_size"].GetDouble();
        mMaxSize = ThisParameters["maximal_size"].GetDouble();
        mEnforceCurrent = ThisParameters["enforce_current"].GetBool();
        
        // In case we have isotropic remeshing (default values)
        if (ThisParameters["anisotropy_remeshing"].GetBool() == false)
        {
            mInterpError = DefaultParameters["hessian_strategy_parameters"]["interpolation_error"].GetDouble();
            mMeshConstant = DefaultParameters["hessian_strategy_parameters"]["mesh_dependent_constant"].GetDouble();
            mAnisRatio = DefaultParameters["anisotropy_parameters"]["hmin_over_hmax_anisotropic_ratio"].GetDouble();
            mBoundLayer = DefaultParameters["anisotropy_parameters"]["boundary_layer_max_distance"].GetDouble();
            mInterpolation = ConvertInter(DefaultParameters["anisotropy_parameters"]["interpolation"].GetString());
        }
        else
        {
            mInterpError = ThisParameters["hessian_strategy_parameters"]["interpolation_error"].GetDouble();
            mMeshConstant = ThisParameters["hessian_strategy_parameters"]["mesh_dependent_constant"].GetDouble();
            mAnisRatio = ThisParameters["anisotropy_parameters"]["hmin_over_hmax_anisotropic_ratio"].GetDouble();
            mBoundLayer = ThisParameters["anisotropy_parameters"]["boundary_layer_max_distance"].GetDouble();
            mInterpolation = ConvertInter(ThisParameters["anisotropy_parameters"]["interpolation"].GetString());
        }
    }
    
    /// Destructor.
    virtual ~ComputeHessianSolMetricProcess() {}
    
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
     * We initialize the metrics of the MMG sol using the Hessian metric matrix approach
     */
    
    virtual void Execute()
    {
        // Iterate in the nodes
        NodesArrayType& pNode = mThisModelPart.Nodes();
        int numNodes = pNode.end() - pNode.begin();
        
        CalculateAuxiliarHessian();
        
        #pragma omp parallel for 
        for(int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            
            if ( itNode->SolutionStepsDataHas( mVariable ) == false )
            {
                KRATOS_ERROR << "Missing variable on node " << itNode->Id() << std::endl;
            }
            
            const double Distance = itNode->FastGetSolutionStepValue(DISTANCE, 0); // TODO: This should be changed for the varaible of interestin the future. This means that the value of the boundary value would be changed to a threshold value instead
            const Vector& Hessian = itNode->GetValue(AUXILIAR_HESSIAN);

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
            
            const double Ratio = CalculateAnisotropicRatio(Distance, mAnisRatio, mBoundLayer, mInterpolation);
            
            // For postprocess pourposes
            double& AnisotropicRatio = itNode->FastGetSolutionStepValue(ANISOTROPIC_RATIO, 0); 
            AnisotropicRatio = Ratio;
            
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
                const Vector NewMetric = ComputeHessianMetricTensor(Hessian, Ratio, ElementMinSize, ElementMaxSize);    
                
                Metric = MetricsMathUtils<TDim>::IntersectMetrics(OldMetric, NewMetric);
            }
            else
            {
                Metric = ComputeHessianMetricTensor(Hessian, Ratio, ElementMinSize, ElementMaxSize);    
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
        return "ComputeHessianSolMetricProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ComputeHessianSolMetricProcess";
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
    TVarType mVariable;            // The variable to calculate the hessian
    double mMinSize;                         // The minimal size of the elements
    double mMaxSize;                         // The maximal size of the elements
    bool mEnforceCurrent;                    // With this we choose if we inforce the current nodal size (NODAL_H)
    double mInterpError;                     // The error of interpolation allowed
    double mMeshConstant;                    // The mesh constant to remesh (depends of the element type)
    double mAnisRatio;                       // The minimal anisotropic ratio (0 < ratio < 1)
    double mBoundLayer;                      // The boundary layer limit distance
    Interpolation mInterpolation;            // The interpolation type
    
    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * This function is used to compute the Hessian Metric tensor, note that when using the Hessian, more than one Metric can be defined simultaneously, so in consecuence we need to define the elipsoid which defines the volume of maximal intersection
     * @param Hessian: The hessian tensor condensed already computed
     * @param AnisotropicRatio: The anisotropic ratio
     * @param ElementMinSize: The min size of element
     * @param ElementMaxSize: The maximal size of the elements
     */
        
    Vector ComputeHessianMetricTensor(
        const Vector& Hessian,
        const double& AnisotropicRatio,
        const double& ElementMinSize, // This way we can impose as minimum as the previous size if we desire
        const double& ElementMaxSize // This way we can impose as maximum as the previous size if we desire
        )
    {        
        // Calculating Metric parameters
        const double CEpsilon = mMeshConstant/mInterpError;
        const double MinRatio = 1.0/(ElementMinSize * ElementMinSize);
//         const double MinRatio = 1.0/(mMinSize * mMinSize);
        const double MaxRatio = 1.0/(ElementMaxSize * ElementMaxSize);
//         const double MaxRatio = 1.0/(mMaxSize * mMaxSize);
        
        typedef boost::numeric::ublas::bounded_matrix<double, TDim, TDim> TempType;
        
        // Declaring the eigen system
        boost::numeric::ublas::bounded_matrix<double, TDim, TDim> EigenVectorMatrix;
        boost::numeric::ublas::bounded_matrix<double, TDim, TDim> EigenValuesMatrix;

        // We first transform into a matrix
        const boost::numeric::ublas::bounded_matrix<double, TDim, TDim> HessianMatrix = MetricsMathUtils<TDim>::VectorToTensor(Hessian);
        
        MathUtils<double>::EigenSystem<TDim>(HessianMatrix, EigenVectorMatrix, EigenValuesMatrix, 1e-18, 20);
        
        // Recalculate the Metric eigen values
        for (unsigned int i = 0; i < TDim; i++)
        {
            EigenValuesMatrix(i, i) = MathUtils<double>::Min(MathUtils<double>::Max(CEpsilon * std::abs(EigenValuesMatrix(i, i)), MaxRatio), MinRatio);
        }
        
        // Considering anisotropic
        if (AnisotropicRatio < 1.0)
        {
            double EigenMax = EigenValuesMatrix(0, 0);
            double EigenMin = EigenValuesMatrix(1, 1);
            for (unsigned int i = 1; i < TDim - 1; i++)
            {
                EigenMax = MathUtils<double>::Max(EigenMax, EigenValuesMatrix(i, i));
                EigenMin = MathUtils<double>::Min(EigenMax, EigenValuesMatrix(i, i));
            }
            
            const double EigenRadius = std::abs(EigenMax - EigenMin) * (1.0 - AnisotropicRatio);
            const double RelativeEigenRadius = std::abs(EigenMax - EigenRadius);
            
            for (unsigned int i = 0; i < TDim; i++)
            {
                EigenValuesMatrix(i, i) = MathUtils<double>::Max(MathUtils<double>::Min(EigenValuesMatrix(i, i), EigenMax), RelativeEigenRadius);
            }
        }
        else // NOTE: For isotropic we should consider the maximum of the eigenvalues
        {
            double EigenMax = EigenValuesMatrix(0, 0);
            for (unsigned int i = 1; i < TDim - 1; i++)
            {
                EigenMax = MathUtils<double>::Max(EigenMax, EigenValuesMatrix(i, i));
            }
            for (unsigned int i = 0; i < TDim; i++)
            {
                EigenValuesMatrix(i, i) = EigenMax;
            }
            EigenVectorMatrix = IdentityMatrix(TDim, TDim);
        }
            
        // We compute the product
        const boost::numeric::ublas::bounded_matrix<double, TDim, TDim> Metric_matrix =  prod(trans(EigenVectorMatrix), prod<TempType>(EigenValuesMatrix, EigenVectorMatrix));
        
        // Finally we transform to a vector
        const Vector Metric = MetricsMathUtils<TDim>::TensorToVector(Metric_matrix);
        
        return Metric;
    }
    
    /**
     * This calculates the auxiliar hessian needed for the Metric
     * @param rThisModelPart: The original model part where we compute the hessian
     * @param rVariable: The variable to calculate the hessian
     */
    
    void CalculateAuxiliarHessian()
    {
        // Iterate in the nodes
        NodesArrayType& pNode = mThisModelPart.Nodes();
        auto numNodes = pNode.end() - pNode.begin();
        
//         #pragma omp parallel for // NOTE: Be careful with the parallel (MUST BE INITIALIZED TO BE THREAD SAFE)
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            
            Vector& hessian = itNode->GetValue(AUXILIAR_HESSIAN);  
            hessian = ZeroVector(3 * (TDim - 1));
        }
        
        // Compute auxiliar gradient
        ComputeNodalGradientProcess<TDim, TVarType> GradientProcess = ComputeNodalGradientProcess<TDim, TVarType>(mThisModelPart, mVariable, AUXILIAR_GRADIENT, NODAL_AREA);
        GradientProcess.Execute();
        
        // Iterate in the conditions
        ElementsArrayType& pElement = mThisModelPart.Elements();
        int numElements = pElement.end() - pElement.begin();
        
        #pragma omp parallel for
        for(int i = 0; i < numElements; i++) 
        {
            auto itElem = pElement.begin() + i;
            
            Element::GeometryType& geom = itElem->GetGeometry();

            double Volume;
            if (geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle2D3)
            {
                boost::numeric::ublas::bounded_matrix<double,3, 2> DN_DX;
                array_1d<double, 3> N;
    
                GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);
                
                boost::numeric::ublas::bounded_matrix<double,3, 2> values;
                for(unsigned int iNode = 0; iNode < 3; iNode++)
                {
                    const array_1d<double, 3> aux_grad = geom[iNode].FastGetSolutionStepValue(AUXILIAR_GRADIENT);
                    values(iNode, 0) = aux_grad[0];
                    values(iNode, 1) = aux_grad[1];
                }
                
                const boost::numeric::ublas::bounded_matrix<double,2, 2> Hessian = prod(trans(DN_DX), values); 
                const Vector HessianCond = MetricsMathUtils<2>::TensorToVector(Hessian);
                
                for(unsigned int iNode = 0; iNode < geom.size(); iNode++)
                {
                    for(unsigned int k = 0; k < 3; k++)
                    {
                        double& val = geom[iNode].GetValue(AUXILIAR_HESSIAN)[k];
                        
                        #pragma omp atomic
                        val += N[iNode] * Volume * HessianCond[k];
                    }
                }
            }
            else if (geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4)
            {
                boost::numeric::ublas::bounded_matrix<double,4,  3> DN_DX;
                array_1d<double, 4> N;
                
                GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);
                
                boost::numeric::ublas::bounded_matrix<double,4, 3> values;
                for(unsigned int iNode = 0; iNode < 4; iNode++)
                {
                    const array_1d<double, 3> aux_grad = geom[iNode].FastGetSolutionStepValue(AUXILIAR_GRADIENT);
                    values(iNode, 0) = aux_grad[0];
                    values(iNode, 1) = aux_grad[1];
                    values(iNode, 2) = aux_grad[2];
                }
                
                const boost::numeric::ublas::bounded_matrix<double, 3, 3> Hessian = prod(trans(DN_DX), values); 
                const Vector HessianCond = MetricsMathUtils<3>::TensorToVector(Hessian);
                
                for(unsigned int iNode = 0; iNode < geom.size(); iNode++)
                {
                    for(unsigned int k = 0; k < 6; k++)
                    {
                        double& val = geom[iNode].GetValue(AUXILIAR_HESSIAN)[k];
                        
                        #pragma omp atomic
                        val += N[iNode] * Volume * HessianCond[k];
                    }
                }
            }
            else
            {
                KRATOS_ERROR << "WARNING: YOU CAN USE JUST 2D TRIANGLES OR 3D TETRAEDRA RIGHT NOW IN THE GEOMETRY UTILS: " << geom.size() << std::endl;
            }
        }
            
        #pragma omp parallel for
        for(int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            itNode->GetValue(AUXILIAR_HESSIAN) /= itNode->FastGetSolutionStepValue(NODAL_AREA);
        }
    }
    
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
template<unsigned int TDim, class TVarType> 
inline std::istream& operator >> (std::istream& rIStream,
                                  ComputeHessianSolMetricProcess<TDim, TVarType>& rThis);

/// output stream function
template<unsigned int TDim, class TVarType> 
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ComputeHessianSolMetricProcess<TDim, TVarType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

};// namespace Kratos.
#endif /* KRATOS_HESSIAN_METRICS_PROCESS defined */
