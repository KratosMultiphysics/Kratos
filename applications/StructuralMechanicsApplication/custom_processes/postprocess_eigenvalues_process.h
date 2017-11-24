//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//

#if !defined(KRATOS_POSTPROCESS_EIGENVALUES_H_INCLUDED )
#define  KRATOS_POSTPROCESS_EIGENVALUES_H_INCLUDED


// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/kratos_parameters.h"

#include "custom_io/gid_eigen_io.h"


namespace Kratos {

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

/// Short class definition.
class PostprocessEigenvaluesProcess : public Process
{
  public:

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of PostprocessEigenvaluesProcess
    KRATOS_CLASS_POINTER_DEFINITION(PostprocessEigenvaluesProcess);

    typedef ModelPart::NodeIterator NodeIterator;
    typedef std::vector<std::vector<std::vector<std::vector<double>>>> AnimationResults;
    /* Structure of this AnimationResults
    AnimationStep
        EigenValues
            Nodes
                Coordinates
    */

    ///@}
    ///@name Life Cycle
    ///@{

    PostprocessEigenvaluesProcess(ModelPart &rModelPart,
                                  const int NumAnimationSteps) 
                                  : mrModelPart(rModelPart),
                                    mNumAnimationSteps(NumAnimationSteps) { }

    /// Destructor.
    ~PostprocessEigenvaluesProcess() = default;

    // Explicitly delete the other constructors
    PostprocessEigenvaluesProcess(const PostprocessEigenvaluesProcess&) = delete;
    PostprocessEigenvaluesProcess& operator=(const PostprocessEigenvaluesProcess&) = delete;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void ExecuteInitialize() override
    {
        const std::string result_file_name = mrModelPart.Name() + "_EigenResults";
        
        mpGidEigenIO = GidEigenIO::Pointer (new GidEigenIO( 
                            result_file_name,
                            GiD_PostBinary,//GiD_PostBinary,
                            MultiFileFlag::SingleFile,
                            WriteDeformedMeshFlag::WriteUndeformed,
                            WriteConditionsFlag::WriteConditions) );

        KRATOS_ERROR_IF_NOT(mpGidEigenIO) << "EigenIO could not be initialized!" << std::endl;
    }

    void ExecuteBeforeSolutionLoop() override
    {
        KRATOS_ERROR_IF_NOT(mpGidEigenIO) << " EigenIO is uninitialized!" << std::endl;

        mpGidEigenIO->InitializeMesh(0.0);
        mpGidEigenIO->WriteMesh(mrModelPart.GetMesh());
        mpGidEigenIO->WriteNodeMesh(mrModelPart.GetMesh());
        mpGidEigenIO->FinalizeMesh();
    }

    void ExecuteFinalize() override
    {
        KRATOS_ERROR_IF_NOT(mpGidEigenIO) << " EigenIO is uninitialized!" << std::endl;

        const int num_nodes = mrModelPart.NumberOfNodes();
        const auto& eigenvalue_vector = mrModelPart.GetProcessInfo()[EIGENVALUE_VECTOR];
        const int num_eigenvalues = eigenvalue_vector.size();

        AnimationResults animated_results(mNumAnimationSteps);
        for (int i=0; i<mNumAnimationSteps; ++i)
        {
            animated_results[i].resize(num_eigenvalues);
            for (int j=0; j<num_eigenvalues; ++j)
            {
                animated_results[i][j].resize(num_nodes);
                for (int k=0; k<num_nodes; ++k)
                {
                    animated_results[i][j][k].resize(3);
                }
            }
        }
        
        std::vector<std::string> labels(num_eigenvalues);
        std::vector<int> nodal_ids(num_nodes);

        for (int i=0; i<num_eigenvalues; ++i)
        {
            std::stringstream strstr;
            strstr << std::fixed << std::setprecision(3) << eigenvalue_vector[i];
            std::string eigenval = strstr.str();

            labels[i] = "EigenValue_" + eigenval;
        }

        NodeIterator it_begin = mrModelPart.NodesBegin();

        // #pragma omp parallel for // Don't modify, this is best suitable for this case
        for (int i=0; i<num_nodes; i++)
        {
            NodeIterator node_it = it_begin + i;

            nodal_ids[i] = node_it->Id();

            const auto& eigen_matrix = node_it->GetValue(EIGENVECTOR_MATRIX);
            // # NOTE: The DoF stored include the velocity and acceleration, solve this in the future

            for (int j=0; j < mNumAnimationSteps; ++j)
            {
                auto& step_results = animated_results[j];

                double angle = 2 * Globals::Pi * j / mNumAnimationSteps;

                for(int k=0; k<num_eigenvalues; ++k)
                {
                    auto& nodal_results = step_results[k][i];
                    nodal_results[0] = std::cos(angle) * eigen_matrix(k, 0);
                    nodal_results[1] = std::cos(angle) * eigen_matrix(k, 1);
                    nodal_results[2] = std::cos(angle) * eigen_matrix(k, 2);
                }
            }
        }

        for (int i=0; i < mNumAnimationSteps; ++i) // "Time"
        {
            mpGidEigenIO->WriteEigenResults(animated_results[i], nodal_ids, labels, i);
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
    virtual std::string Info() const {
        return "PostprocessEigenvaluesProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {
        rOStream << "PostprocessEigenvaluesProcess";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

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

    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{
    ModelPart& mrModelPart;
    int mNumAnimationSteps;
    GidEigenIO::Pointer mpGidEigenIO;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}

}; // Class PostprocessEigenvaluesProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.

#endif // KRATOS_POSTPROCESS_EIGENVALUES_H_INCLUDED  defined
