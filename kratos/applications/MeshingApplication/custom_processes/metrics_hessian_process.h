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
     * @param rMinSize: The min size of element
     * @param rMaxSize: The maximal size of the elements
     * @param rInterpError: The interpolation error assumed
     * @param rMeshConstant: The constant that appears in papers an none knows where comes from, where...
     * @param rAnisRatio: The anisotropic ratio
     * @param rBoundLayer: The boundary layer limit
     * @param rInterpolation: The interpolation type
     * @param rVariable: The variable considered for remeshing
     */
    
    ComputeHessianSolMetricProcess(
        ModelPart& rThisModelPart,
        TVarType& rVariable,
        const double rMinSize,
        const double rMaxSize,
        const bool rEnforceCurrent = true,
        const double rInterpError = 1.0e-6,
        const double rMeshConstant = 0.28125, // TODO: This is for tetrahedron, look in literature for more
        const double rAnisRatio = 1.0,
        const double rBoundLayer =  1.0,
        const std::string rInterpolation = "Linear"
        )
        :mThisModelPart(rThisModelPart),
        mVariable(rVariable),
        mMinSize(rMinSize),
        mMaxSize(rMaxSize),
        mEnforceCurrent(rEnforceCurrent),
        mInterpError(rInterpError),
        mMeshConstant(rMeshConstant),
        mAnisRatio(rAnisRatio),
        mBoundLayer(rBoundLayer)
    {       
        mInterpolation = ConvertInter(rInterpolation);
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
                KRATOS_THROW_ERROR( std::invalid_argument, "Missing variable on node ", itNode->Id() )
            }
            
            const double distance = itNode->FastGetSolutionStepValue(DISTANCE, 0);
            const Vector& hessian = itNode->GetValue(AUXILIAR_HESSIAN);

            const double nodal_h = itNode->FastGetSolutionStepValue(NODAL_H, 0);            
            
            double element_min_size = mMinSize;
            if ((element_min_size > nodal_h) && (mEnforceCurrent == true))
            {
                element_min_size = nodal_h;
            }
            double element_max_size = mMaxSize;
            if ((element_max_size > nodal_h) && (mEnforceCurrent == true))
            {
                element_max_size = nodal_h;
            }
            
            const double ratio = CalculateAnisotropicRatio(distance, mAnisRatio, mBoundLayer, mInterpolation);
            
            // For postprocess pourposes
            double& anisotropic_ratio = itNode->FastGetSolutionStepValue(ANISOTROPIC_RATIO, 0); 
            anisotropic_ratio = ratio;
            
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
                const Vector new_metric = ComputeHessianMetricTensor(hessian, ratio, element_min_size, element_max_size);    
                
                metric = MetricsMathUtils<TDim>::IntersectMetrics(old_metric, new_metric);
            }
            else
            {
                metric = ComputeHessianMetricTensor(hessian, ratio, element_min_size, element_max_size);    
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
     * This function is used to compute the Hessian metric tensor, note that when using the Hessian, more than one metric can be defined simultaneously, so in consecuence we need to define the elipsoid which defines the volume of maximal intersection
     * @param hessian: The hessian tensor condensed already computed
     * @param ratio: The anisotropic ratio
     */
        
    Vector ComputeHessianMetricTensor(
        const Vector& hessian,
        const double& ratio,
        const double& element_min_size, // This way we can impose as minimum as the previous size if we desire
        const double& element_max_size // This way we can impose as maximum as the previous size if we desire
        )
    {        
        // Calculating metric parameters
        const double c_epsilon = mMeshConstant/mInterpError;
        const double min_ratio = 1.0/(element_min_size * element_min_size);
//         const double min_ratio = 1.0/(mMinSize * mMinSize);
        const double max_ratio = 1.0/(element_max_size * element_max_size);
//         const double max_ratio = 1.0/(mMaxSize * mMaxSize);
        
        typedef boost::numeric::ublas::bounded_matrix<double, TDim, TDim> temp_type;
        
        // Declaring the eigen system
        boost::numeric::ublas::bounded_matrix<double, TDim, TDim> eigen_vector_matrix;
        boost::numeric::ublas::bounded_matrix<double, TDim, TDim> eigen_values_matrix;

        // We first transform into a matrix
        const boost::numeric::ublas::bounded_matrix<double, TDim, TDim> hessian_matrix = MetricsMathUtils<TDim>::VectorToTensor(hessian);
        
        MetricsMathUtils<TDim>::EigenSystem(hessian_matrix, eigen_vector_matrix, eigen_values_matrix, 1e-18, 20);
        
        // Recalculate the metric eigen values
        for (unsigned int i = 0; i < TDim; i++)
        {
            eigen_values_matrix(i, i) = MathUtils<double>::Min(MathUtils<double>::Max(c_epsilon * std::abs(eigen_values_matrix(i, i)), max_ratio), min_ratio);
        }
        
        // Considering anisotropic
        if (ratio < 1.0)
        {
            double eigen_max = eigen_values_matrix(0, 0);
            double eigen_min = eigen_values_matrix(1, 1);
            for (unsigned int i = 1; i < TDim - 1; i++)
            {
                eigen_max = MathUtils<double>::Max(eigen_max, eigen_values_matrix(i, i));
                eigen_min = MathUtils<double>::Min(eigen_max, eigen_values_matrix(i, i));
            }
            
            const double eigen_radius = std::abs(eigen_max - eigen_min) * (1.0 - ratio);
            const double rel_eigen_radius = std::abs(eigen_max - eigen_radius);
            
            for (unsigned int i = 0; i < TDim; i++)
            {
                eigen_values_matrix(i, i) = MathUtils<double>::Max(MathUtils<double>::Min(eigen_values_matrix(i, i), eigen_max), rel_eigen_radius);
            }
        }
        else // NOTE: For isotropic we should consider the maximum of the eigenvalues
        {
            double eigen_max = eigen_values_matrix(0, 0);
            for (unsigned int i = 1; i < TDim - 1; i++)
            {
                eigen_max = MathUtils<double>::Max(eigen_max, eigen_values_matrix(i, i));
            }
            for (unsigned int i = 0; i < TDim; i++)
            {
                eigen_values_matrix(i, i) = eigen_max;
            }
            eigen_vector_matrix = IdentityMatrix(TDim, TDim);
        }
            
        // We compute the product
        const boost::numeric::ublas::bounded_matrix<double, TDim, TDim> metric_matrix =  prod(trans(eigen_vector_matrix), prod<temp_type>(eigen_values_matrix, eigen_vector_matrix));
        
        // Finally we transform to a vector
        const Vector metric = MetricsMathUtils<TDim>::TensorToVector(metric_matrix);
        
        return metric;
    }
    
    /**
     * This calculates the auxiliar hessian needed for the metric
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
                for(unsigned int i_node = 0; i_node < 3; i_node++)
                {
                    const array_1d<double, 3> aux_grad = geom[i_node].FastGetSolutionStepValue(AUXILIAR_GRADIENT);
                    values(i_node, 0) = aux_grad[0];
                    values(i_node, 1) = aux_grad[1];
                }
                
                const boost::numeric::ublas::bounded_matrix<double,2, 2> Hessian = prod(trans(DN_DX), values); 
                const Vector HessianCond = MetricsMathUtils<2>::TensorToVector(Hessian);
                
                for(unsigned int i_node = 0; i_node < geom.size(); i_node++)
                {
                    for(unsigned int k = 0; k < 3; k++)
                    {
                        double& val = geom[i_node].GetValue(AUXILIAR_HESSIAN)[k];
                        
                        #pragma omp atomic
                        val += N[i_node] * Volume * HessianCond[k];
                    }
                }
            }
            else if (geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4)
            {
                boost::numeric::ublas::bounded_matrix<double,4,  3> DN_DX;
                array_1d<double, 4> N;
                
                GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);
                
                boost::numeric::ublas::bounded_matrix<double,4, 3> values;
                for(unsigned int i_node = 0; i_node < 4; i_node++)
                {
                    const array_1d<double, 3> aux_grad = geom[i_node].FastGetSolutionStepValue(AUXILIAR_GRADIENT);
                    values(i_node, 0) = aux_grad[0];
                    values(i_node, 1) = aux_grad[1];
                    values(i_node, 2) = aux_grad[2];
                }
                
                const boost::numeric::ublas::bounded_matrix<double, 3, 3> Hessian = prod(trans(DN_DX), values); 
                const Vector HessianCond = MetricsMathUtils<3>::TensorToVector(Hessian);
                
                for(unsigned int i_node = 0; i_node < geom.size(); i_node++)
                {
                    for(unsigned int k = 0; k < 6; k++)
                    {
                        double& val = geom[i_node].GetValue(AUXILIAR_HESSIAN)[k];
                        
                        #pragma omp atomic
                        val += N[i_node] * Volume * HessianCond[k];
                    }
                }
            }
            else
            {
                KRATOS_THROW_ERROR( std::logic_error, "WARNING: YOU CAN USE JUST 2D TRIANGLES OR 3D TETRAEDRA RIGHT NOW IN THE GEOMETRY UTILS: ", geom.size() );
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
     * @param str: The string
     * @param Interpolation: The equivalent enum
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
