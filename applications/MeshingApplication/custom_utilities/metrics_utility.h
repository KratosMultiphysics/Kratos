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

#if !defined(KRATOS_METRICS_UTILITY)
#define KRATOS_METRICS_UTILITY

// Project includes
#include "utilities/math_utils.h"
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"
#include "processes/compute_nodal_gradient_process.h" // TODO: Not prism or quadrilaterals implemented yet

namespace Kratos
{
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

    /**
     * This function converts a reference to values that can be used in the computation of tensors in 2D
     */
    
    void TensorReferenceConversor2D(
        const unsigned int& ref,
        unsigned int& i,
        unsigned int& j
        )
    {
        if (ref == 0) 
        {
            i = 0;
            j = 0;
        }
        else if (ref == 1)
        {
            i = 0;
            j = 1;
        }
        else if (ref == 2) 
        {
            i = 1;
            j = 1;
        }
    }
    
    /**
     * This function converts a reference to values that can be used in the computation of tensors in 3D
     */
        
    void TensorReferenceConversor3D(
        const unsigned int& ref,
        unsigned int& i,
        unsigned int& j
        )
    {
        if (ref == 0) 
        {
            i = 0;
            j = 0;
        }
        else if (ref == 1)
        {
            i = 0;
            j = 1;
        }
        else if (ref == 2) 
        {
            i = 0;
            j = 2;
        }
        else if (ref == 3) 
        {
            i = 1;
            j = 1;
        }
        else if (ref == 4)
        {
            i = 1;
            j = 2;
        }
        else if (ref == 5) 
        {
            i = 2;
            j = 2;
        }
    }
    
    /**
     * Calculates the maximum of 3 values
     */  
    
    double max3(
        const double val1,
        const double val2,
        const double val3
        )
    {
        double max = val1;
        if (max < val2)
        {
            max = val2;
        }
        if (max < val3)
        {
            max = val3;
        }
        
        return max;
    }
    
    /**
     * Calculates the minimum of 3 values
     */  
    
    double min3(
        const double val1,
        const double val2,
        const double val3
        )
    {
        double min = val1;
        if (min > val2)
        {
            min = val2;
        }
        if (min > val3)
        {
            min = val3;
        }
        
        return min;
    }
    
    /**
     * Calculates the eigenvectors and eigenvalues of given symmetric 3x3 matrix
     * The eigenvectors and eigenvalues are calculated using the iterative Gauss-Seidel-method
     * @param A: The given symmetric matrix the eigenvectors are to be calculated.
     * @return eigen_vector_matrix: The result matrix (will be overwritten with the eigenvectors)
     * @return lambda: The result diagonal matrix with the eigenvalues
     * @param tolerance: The largest value considered to be zero
     * @param max_iterations: Maximum number of iterations
     */

    static inline bool EigenSystem3x3(
            const boost::numeric::ublas::bounded_matrix<double, 3, 3>& A,
            boost::numeric::ublas::bounded_matrix<double, 3, 3>& eigen_vector_matrix,
            boost::numeric::ublas::bounded_matrix<double, 3, 3>& eigen_values_matrix,
            const double tolerance = 1.0e-9,
            const unsigned int max_iterations = 10
            )
    {
        bool is_converged = false;
        boost::numeric::ublas::bounded_matrix<double, 3, 3> auxA;
        const boost::numeric::ublas::bounded_matrix<double, 3, 3> Identity = IdentityMatrix(3, 3);
        boost::numeric::ublas::bounded_matrix<double, 3, 3> V = Identity;
        boost::numeric::ublas::bounded_matrix<double, 3, 3> Rotation;

        for(unsigned int iterations = 0; iterations < max_iterations; iterations++)
        {
            is_converged = true;

            double a = 0.0;
            unsigned int index1 = 0;
            unsigned int index2 = 0;

            for(unsigned int i = 0; i < 3; i++)
            {
                for(unsigned int j = (i + 1); j < 3; j++)
                {
                    if((std::abs(A(i, j)) > a ) && (std::abs(A(i, j)) > tolerance))
                    {
                        a = std::abs(A(i,j));
                        index1 = i;
                        index2 = j;
                        is_converged = false;
                    }
                }
            }

            if(is_converged == true)
            {
                break;
            }

            // Calculation of Rotation angle
            const double gamma = (A(index2, index2)-A(index1, index1)) / (2.0 * A(index1, index2));
            double u = 1.0;

            if(std::abs(gamma) > tolerance && std::abs(gamma)< (1.0/tolerance))
            {
                u = gamma / std::abs(gamma) * 1.0 / (std::abs(gamma) + std::sqrt(1.0 + gamma * gamma));
            }
            else
            {
                if  (std::abs(gamma) >= (1.0/tolerance))
                {
                    u = 0.5 / gamma;
                }
            }

            const double c = 1.0 / (sqrt(1.0 + u * u));
            const double s = c * u;
            const double teta = s / (1.0 + c);

            // Rotation of the Matrix
            auxA = A;
            auxA(index2, index2) = A(index2,index2) + u * A(index1, index2);
            auxA(index1, index1) = A(index1,index1) - u * A(index1, index2);
            auxA(index1, index2) = 0.0;
            auxA(index2, index1) = 0.0;

            for(unsigned int i = 0; i < 3; i++)
            {
                if((i!= index1) && (i!= index2))
                {
                    auxA(index2, i) = A(index2, i) + s * (A(index1, i) - teta * A(index2, i));
                    auxA(i, index2) = A(index2, i) + s * (A(index1, i) - teta * A(index2, i));
                    auxA(index1, i) = A(index1, i) - s * (A(index2, i) + teta * A(index1, i));
                    auxA(i, index1) = A(index1, i) - s * (A(index2, i) + teta * A(index1, i));
                }
            }

            // Calculation of the eigenvectors V
            Rotation = Identity;
            Rotation(index2, index1) = -s;
            Rotation(index1, index2) =  s;
            Rotation(index1, index1) =  c;
            Rotation(index2, index2) =  c;

            eigen_vector_matrix = ZeroMatrix(3, 3);
            for(unsigned int i = 0; i < 3; i++)
            {
                for(unsigned int j = 0; j < 3; j++)
                {
                    for(unsigned int k = 0; k < 3; k++)
                    {
                        eigen_vector_matrix(j, i) += V(i, k) * Rotation(k, j);
                    }
                }
            }
        }

        if(!(is_converged))
        {
            std::cout<<" WARNING: Spectral decomposition not converged "<<std::endl;
        }

        for(unsigned int i = 0; i < 3; i++)
        {
            for(unsigned int j = 0; j < 3; j++)
            {
                if (i == j)
                {
                    eigen_values_matrix(i, i) = auxA(i, i);
                }
                else
                {
                    eigen_values_matrix(i, j) = 0.0;
                }
            }
        }

        return is_converged;
    }
    
///@}
///@name Kratos Classes
///@{

//// This class is can be used to compute

template<unsigned int TDim>  
class MetricsUtility
{
public:

    ///@name Type Definitions
    ///@{
    
    typedef ModelPart::NodesContainerType                        NodesArrayType;
    typedef ModelPart::ElementsContainerType                  ElementsArrayType;
    typedef ModelPart::ConditionsContainerType              ConditionsArrayType;
    typedef Node <3>                                                   NodeType;
    
    ///@}
    ///@name Life Cycle
    ///@{
     
    // Constructor
    
    /**
     * This is the default constructor
     * @param rThisModelPart: The model part to be computed
     * @param rMinSize: The min size of element
     * @param rAnisRatio: The anisotropic ratio
     * @param rBoundLayer: The boundary layer limit
     * @param rInterpolation: The interpolation type
     */
    
    MetricsUtility(
        const double rMinSize,
        const double rAnisRatio = 1.0,
        const double rBoundLayer =  1.0,
        const std::string rInterpolation = "Linear"
        )
        :mMinSize(rMinSize),
        mAnisRatio(rAnisRatio),
        mBoundLayer(rBoundLayer),
        mInterpolation(rInterpolation)
    {       
        std::cout << "Initializing metric utility" << std::endl;
    }
    
    /// Destructor.
    ~MetricsUtility() {}
    
    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
    
    /**
     * We initialize the metrics of the MMG sol using a level set approach
     * @param rVariableGradient: The gradient variable to compute
     */
    
    void ComputeLevelSetSolMetric(
        ModelPart& rThisModelPart,
        const Variable<array_1d<double,3>> rVariableGradient = DISTANCE_GRADIENT
        )
    {
        // Iterate in the nodes
        NodesArrayType& pNode = rThisModelPart.Nodes();
        auto numNodes = pNode.end() - pNode.begin();

//         #pragma omp parallel for 
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            
            const double distance = itNode->FastGetSolutionStepValue(DISTANCE, 0);
            array_1d<double, 3> gradient_value = itNode->FastGetSolutionStepValue(rVariableGradient, 0);
            
            double element_size = itNode->FastGetSolutionStepValue(NODAL_H, 0);
            if (element_size > mMinSize)
            {
                element_size = mMinSize;
            }
            
            const double ratio = CalculateAnisotropicRatio(distance);
            
            // For postprocess pourposes
            double& anisotropic_ratio = itNode->FastGetSolutionStepValue(ANISOTROPIC_RATIO, 0); 
            anisotropic_ratio = ratio;
            
            const double tolerance = 1.0e-12;
            const double norm = norm_2(gradient_value);
            if (norm > tolerance)
            {
                gradient_value /= norm;
            }
            
            // We compute the metric
            Vector& metric = itNode->GetValue(MMG_METRIC);
            metric = ComputeLevelSetMetricTensor(gradient_value, ratio, element_size);
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * We initialize the metrics of the MMG sol using the Hessian metric matrix approach
     * @param rMaxSize: The maximal size of the elements
     * @param rVariable: The variable considered for remeshing
     * @param rInterpError: The interpolation error assumed
     * @param rMeshConstant: The constant that appears in papers an none knows where comes from, where...
     */
    
    void ComputeHessianMetric(
        ModelPart& rThisModelPart,
        const double& rMaxSize,
        Variable<double> & rVariable,
        const double rInterpError = 1.0e-6,
        const double rMeshConstant = 0.28125 // TODO: This is for tetrahedron, look in literature for more
        )
    {
        // Iterate in the nodes
        NodesArrayType& pNode = rThisModelPart.Nodes();
        auto numNodes = pNode.end() - pNode.begin();

        CalculateAuxiliarHessian(rThisModelPart, rVariable);
        
//         #pragma omp parallel for 
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            
            const double distance = itNode->FastGetSolutionStepValue(DISTANCE, 0);
            const Vector& hessian = itNode->GetValue(AUXILIAR_HESSIAN);
            
            double element_size = itNode->FastGetSolutionStepValue(NODAL_H, 0);
            if (element_size > mMinSize)
            {
                element_size = mMinSize;
            }
            
            const double ratio = CalculateAnisotropicRatio(distance);
            
            // For postprocess pourposes
            double& anisotropic_ratio = itNode->FastGetSolutionStepValue(ANISOTROPIC_RATIO, 0); 
            anisotropic_ratio = ratio;
            
            // We compute the metric
            const array_1d<double, 3 * (1 + TDim)> metric = ComputeHessianMetricTensor(hessian, element_size, rMaxSize, rInterpError, rMeshConstant, ratio);
//             const array_1d<double, 3 * (1 + TDim)> metric = ComputeHessianMetricTensor(hessian, mMinSize, rMaxSize, rInterpError, rMeshConstant, ratio);         
        }
    }
    
    /**
     * The same, but with components
     */
    
    void ComputeHessianMetricComponents(
        ModelPart& rThisModelPart,
        const double& rMaxSize,
        Variable<array_1d<double,3>> & rVariable, // NOTE: In the case we consider some composed variable
        const double rInterpError = 1.0e-6,
        const double rMeshConstant = 0.28125 // TODO: This is for tetrahedron, look in literature for more
        )
    {//TODO: look cea2010_V2.pdf pag. 25
//         // Iterate in the nodes
//         NodesArrayType& pNode = rThisModelPart.Nodes();
//         auto numNodes = pNode.end() - pNode.begin();
// 
//         CalculateAuxiliarHessian(rThisModelPart, rVariable); // TODO: Solve this
//         
//         for(unsigned int i = 0; i < numNodes; i++) 
//         {
//             auto itNode = pNode.begin() + i;
//             
//             const Vector& hessian = itNode->GetValue(AUXILIAR_HESSIAN);
//             
//             // We compute the metric
//             const array_1d<double, 3 * (1 + TDim)> metric = ComputeHessianMetricTensor(hessian, mMinSize, rMaxSize, rInterpError, rMeshConstant, ratio);      
//         }
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
    
protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    double mMinSize;            // The minimal size of the elements
    double mAnisRatio;          // The minimal anisotropic ratio (0 < ratio < 1)
    double mBoundLayer;         // The boundary layer limit distance
    std::string mInterpolation; // The interpolation type
    
    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
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
    )
    {
        Vector metric;
        metric.resize(3 * (1 + TDim), false);
        
        const double coeff0 = 1.0/(element_size * element_size);
        
        if (TDim == 2) // 2D: The order of the metric is m11,m12,m22
        {
            const double coeff1 = coeff0/(ratio * ratio);
            
            const double v0v0 = gradient_value[0]*gradient_value[0];
            const double v0v1 = gradient_value[0]*gradient_value[1];
            const double v1v1 = gradient_value[1]*gradient_value[1];
            
            metric[0] = coeff0*(1.0 - v0v0) + coeff1*v0v0;
            metric[1] = coeff0*(    - v0v1) + coeff1*v0v1;  
            metric[2] = coeff0*(1.0 - v1v1) + coeff1*v1v1;

        }
        else // 3D: The order of the metric is m11,m12,m13,m22,m23,m33
        {
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
        }
        
        return metric;
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    /**
     * This function is used to compute the Hessian metric tensor, note that when using the Hessian, more than one metric can be defined simultaneously, so in consecuence we need to define the elipsoid which defines the volume of maximal intersection
     * @param hessian: The hessian tensor condensed already computed
     * @param mMinSize: The minimum size of the elements
     * @param rMaxSize: The maximal size of the elements
     * @param rInterpError: The interpolation error assumed
     * @param rMeshConstant: The constant that appears in papers an none knows where comes from, where...
     */
        
    Vector ComputeHessianMetricTensor(
        const Vector& hessian,
        const double& mMinSize,
        const double& rMaxSize,
        const double& rInterpError,
        const double& rMeshConstant,
        const double& ratio
        )
    {
        Vector metric;
        metric.resize(3 * (1 + TDim), false);
        
        // Calculating metric parameters
        const double c_epsilon = rMeshConstant/rInterpError;
        const double min_ratio = 1.0/(mMinSize * mMinSize);
        const double max_ratio = 1.0/(rMaxSize * rMaxSize);
        
        typedef boost::numeric::ublas::bounded_matrix<double, 2, 2> temp_type;
        
        if (TDim == 2)
        {
            const double tol = 1.0e-16;
            
            // Declaring the eigen system
            boost::numeric::ublas::bounded_matrix<double, 2, 2> eigen_vector_matrix;
            boost::numeric::ublas::bounded_matrix<double, 2, 2> eigen_values_matrix;
            
            // Taking the coefficients of the hessian
            const double h11 = hessian[0];
            const double h12 = hessian[1] + tol; // NOTE: To avoid problems with diagonal or almost diagonal matrices
            const double h22 = hessian[2];
            
            // The eigen values matrix 
            eigen_values_matrix(0, 0) = MathUtils<double>::Min(MathUtils<double>::Max(c_epsilon * 0.5 * std::abs(h11 + h22 - std::sqrt(h11 * h11 + 4.0 * h12 * h12 - 2.0 * h11 * h22 + h22 * h22)), max_ratio), min_ratio);
            eigen_values_matrix(0, 1) = 0.0; 
            eigen_values_matrix(1, 0) = 0.0; 
            eigen_values_matrix(1, 1) = MathUtils<double>::Min(MathUtils<double>::Max(c_epsilon * 0.5 * std::abs(h11 + h22 + std::sqrt(h11 * h11 + 4.0 * h12 * h12 - 2.0 * h11 * h22 + h22 * h22)), max_ratio), min_ratio);
            
            // The eigen vector matrix
            array_1d<double, 2> aux_array;
            // First eigen vector
            aux_array[0] = -(-h11 + h22 + std::sqrt(h11 * h11 + 4.0 * h12 * h12 - 2 * h11 * h22 + h22 * h22))/(2.0 * h12);
            aux_array[1] = 1.0;
            aux_array /= norm_2(aux_array);
            eigen_vector_matrix(0, 0) = aux_array[0];
            eigen_vector_matrix(0, 1) = aux_array[1];
            // Second eigen vector
            aux_array[0] = -(-h11 + h22 - std::sqrt(h11 * h11 + 4.0 * h12 * h12 - 2 * h11 * h22 + h22 * h22))/(2.0 * h12);
            aux_array[1] = 1.0;
            aux_array /= norm_2(aux_array);
            eigen_vector_matrix(1, 0) = aux_array[0];
            eigen_vector_matrix(1, 1) = aux_array[1];
            
            // Considering anisotropic
            if (ratio < 1.0)
            {
                const double eigen_max = MathUtils<double>::Max(eigen_values_matrix(0, 0), eigen_values_matrix(1, 1));
                const double eigen_min = MathUtils<double>::Min(eigen_values_matrix(0, 0), eigen_values_matrix(1, 1));
                const double eigen_radius = std::abs(eigen_max - eigen_min) * (1.0 - ratio);
                const double rel_eigen_radius = std::abs(eigen_max - eigen_radius);
                
                eigen_values_matrix(0, 0) = MathUtils<double>::Max(MathUtils<double>::Min(eigen_values_matrix(0, 0), eigen_max), rel_eigen_radius);
                eigen_values_matrix(1, 1) = MathUtils<double>::Max(MathUtils<double>::Min(eigen_values_matrix(1, 1), eigen_max), rel_eigen_radius);
            }
            
            // We compute the product
            const boost::numeric::ublas::bounded_matrix<double, 2, 2> metric_matrix =  prod(trans(eigen_vector_matrix), prod<temp_type>(eigen_values_matrix,eigen_vector_matrix));
            
            // Finally we transform to a vector
            for (unsigned int k = 0; k < 3; k++)
            {
                unsigned int index0 ,index1;
                TensorReferenceConversor2D(k, index0, index1);
                metric[k] = metric_matrix(index0, index1);
            }
        }
        else
        {
            // We first transform into a matrix
            boost::numeric::ublas::bounded_matrix<double, 3, 3> hessian_matrix;
            for(unsigned int k = 0; k < 6; k++)
            {
                unsigned int index0 ,index1;
                TensorReferenceConversor3D(k, index0, index1);
                hessian_matrix(index0, index1) = hessian[k];
            }
            
            // Declaring the eigen system
            boost::numeric::ublas::bounded_matrix<double, 3, 3> eigen_vector_matrix;
            boost::numeric::ublas::bounded_matrix<double, 3, 3> eigen_values_matrix;
            
            EigenSystem3x3(hessian_matrix, eigen_vector_matrix, eigen_values_matrix, 1e-18, 10);
            
            // Recalculate the metric eigen values
            for (unsigned int i = 0; i < 3; i++)
            {
                eigen_values_matrix(i, i) = MathUtils<double>::Min(MathUtils<double>::Max(c_epsilon * std::abs(eigen_values_matrix(i, i)), max_ratio), min_ratio);
            }
            
            // Considering anisotropic
            if (ratio < 1.0)
            {
                const double eigen_max = max3(eigen_values_matrix(0, 0), eigen_values_matrix(1, 1), eigen_values_matrix(2, 2));
                const double eigen_min = min3(eigen_values_matrix(0, 0), eigen_values_matrix(1, 1), eigen_values_matrix(2, 2));
                const double eigen_radius = std::abs(eigen_max - eigen_min) * (1.0 - ratio);
                const double rel_eigen_radius = std::abs(eigen_max - eigen_radius);
                
                eigen_values_matrix(0, 0) = MathUtils<double>::Max(MathUtils<double>::Min(eigen_values_matrix(0, 0), eigen_max), rel_eigen_radius);
                eigen_values_matrix(1, 1) = MathUtils<double>::Max(MathUtils<double>::Min(eigen_values_matrix(1, 1), eigen_max), rel_eigen_radius);
                eigen_values_matrix(1, 2) = MathUtils<double>::Max(MathUtils<double>::Min(eigen_values_matrix(2, 2), eigen_max), rel_eigen_radius);
            }
            
            // We compute the product
            const boost::numeric::ublas::bounded_matrix<double, 3, 3> metric_matrix =  prod(trans(eigen_vector_matrix), prod<temp_type>(eigen_values_matrix,eigen_vector_matrix));
            
            // Finally we transform to a vector
            for (unsigned int k = 0; k < 6; k++)
            {
                unsigned int index0 ,index1;
                TensorReferenceConversor3D(k, index0, index1);
                metric[k] = metric_matrix(index0, index1);
            }
        }
        
        return metric;
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * This calculates the auxiliar hessian needed for the metric
     * @param rThisModelPart: The original model part where we compute the hessian
     * @param rVariable: The variable to calculate the hessian
     */
    
    void CalculateAuxiliarHessian(
        ModelPart& rThisModelPart,
        Variable<double> & rVariable
        )
    {
        // Iterate in the nodes
        NodesArrayType& pNode = rThisModelPart.Nodes();
        auto numNodes = pNode.end() - pNode.begin();
        
//         #pragma omp parallel for // NOTE: Be careful with the parallel (MUST BE INITIALIZED TO BE THREAD SAFE)
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            
            Vector& hessian = itNode->GetValue(AUXILIAR_HESSIAN);  
            hessian = ZeroVector(3 * (TDim - 1));
        }
        
        // Compute auxiliar gradient
        ComputeNodalGradientProcess<TDim> GradientProcess = ComputeNodalGradientProcess<TDim>(rThisModelPart, rVariable, AUXILIAR_GRADIENT, NODAL_AREA);
        GradientProcess.Execute();
        
        // Iterate in the conditions
        ElementsArrayType& pElement = rThisModelPart.Elements();
        auto numElements = pElement.end() - pElement.begin();
        
        #pragma omp parallel for
        for(unsigned int i = 0; i < numElements; i++) 
        {
            auto itElem = pElement.begin() + i;
            
            Element::GeometryType& geom = itElem->GetGeometry();
            const unsigned int geom_size = geom.size();

            double Volume;
            if (TDim == 2)
            {
                if (geom_size == 3)
                {
                    boost::numeric::ublas::bounded_matrix<double,3, 2> DN_DX;
                    array_1d<double,3> N;
       
                    GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);
                    
                    boost::numeric::ublas::bounded_matrix<double,3, 2> values;
                    for(unsigned int i_node = 0; i_node < 3; i_node++)
                    {
                        const array_1d<double, 3> aux_grad = geom[i_node].FastGetSolutionStepValue(AUXILIAR_GRADIENT);
                        values(i_node, 0) = aux_grad[0];
                        values(i_node, 1) = aux_grad[1];
                    }
                    
                    const boost::numeric::ublas::bounded_matrix<double,2, 2> hessian = prod(trans(DN_DX), values); 
                    
                    for(unsigned int i_node = 0; i_node < geom_size; i_node++)
                    {
                        for(unsigned int k = 0; k < 3; k++)
                        {
                            unsigned int index0 ,index1;
                            TensorReferenceConversor2D(k, index0, index1);
                            
                            double& val = geom[i_node].FastGetSolutionStepValue(AUXILIAR_HESSIAN)[k];
                            
                            #pragma omp atomic
                            val += N[i_node] * Volume * hessian(index0, index1);
                        }
                    }
                }
                else
                {
                    KRATOS_THROW_ERROR( std::logic_error, "WARNING: YOU CAN USE JUST TRIANGLES RIGHT NOW IN THE GEOMETRY UTILS: ", geom_size );
                }
            }
            else
            {
                if (geom_size == 4)
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
                    
                    const boost::numeric::ublas::bounded_matrix<double,3, 3> hessian = prod(trans(DN_DX), values); 
                    
                    for(unsigned int i_node = 0; i_node < geom_size; i_node++)
                    {
                        for(unsigned int k = 0; k < 6; k++)
                        {
                            unsigned int index0 ,index1;
                            TensorReferenceConversor3D(k, index0, index1);

                            double& val = geom[i_node].FastGetSolutionStepValue(AUXILIAR_HESSIAN)[k];
                            
                            #pragma omp atomic
                            val += N[i_node] * Volume * hessian(index0, index1);
                        }
                    }
                }
                else
                {
                    KRATOS_THROW_ERROR( std::logic_error, "WARNING: YOU CAN USE JUST TETRAEDRA RIGHT NOW IN THE GEOMETRY UTILS: ", geom.size() );
                }
            }
        }
            
        #pragma omp parallel for
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            itNode->GetValue(AUXILIAR_HESSIAN) /= itNode->FastGetSolutionStepValue(NODAL_AREA);
            
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * This calculates the anisotropic ratio
     * @param distance: Distance parameter
     */
    
    double CalculateAnisotropicRatio(const double& distance)
    {
        const double tolerance = 1.0e-12;
        double ratio = 1.0; // NOTE: Isotropic mesh
        if (mAnisRatio < 1.0)
        {                           
            if (std::abs(distance) <= mBoundLayer)
            {
                if (mInterpolation.find("Constant") != std::string::npos)
                {
                    ratio = mAnisRatio;
                }
                else if (mInterpolation.find("Linear") != std::string::npos)
                {
                    ratio = mAnisRatio + (std::abs(distance)/mBoundLayer) * (1.0 - mAnisRatio);
                }
                else if (mInterpolation.find("Exponential") != std::string::npos)
                {
                    ratio = - std::log(std::abs(distance)/mBoundLayer) * mAnisRatio + tolerance;
                    if (ratio > 1.0)
                    {
                        ratio = 1.0;
                    }
                }
                else
                {
                    std::cout << "No interpolation defined, considering linear" << std:: endl;
                    ratio = 
                    mAnisRatio + (std::abs(distance)/mBoundLayer) * (1.0 - mAnisRatio);
                }
            }
        }
        
        return ratio;
    }
    
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
};// class MetricsUtility
///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

}// namespace Kratos.
#endif /* KRATOS_METRICS_UTILITY defined */
