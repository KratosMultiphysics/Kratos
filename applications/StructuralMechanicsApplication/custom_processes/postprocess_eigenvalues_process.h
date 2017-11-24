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

/// Process to create the animated Eigenvectors
/** This process takes the results of an Eigenvalue Analysis and creates the 
 * animated Eigenvectors (Eigenmodes) for GiD using the GidEigenIO, which
 * is customized for this case
 * The Input should be the ComputingModelPart! (Otherwise nodal results migth be messed up)
 * It is of particular importance that all Nodes have the same Dofs!
 */
class PostprocessEigenvaluesProcess : public Process
{
  public:

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of PostprocessEigenvaluesProcess
    KRATOS_CLASS_POINTER_DEFINITION(PostprocessEigenvaluesProcess);

    typedef std::size_t SizeType;

    typedef ModelPart::NodeIterator NodeIterator;
    typedef ModelPart::NodeType::DofsContainerType DofsContainerType;
    
    typedef std::vector<std::vector<std::vector<std::vector<double>>>> AnimationResults;
    /* The main reason to do this is to access the nodal database only once and computing 
       the animated results
    Structure of this AnimationResults:
    AnimationStep
        EigenValues
            Nodes
                Coordinates
    */

    ///@}
    ///@name Life Cycle
    ///@{

    PostprocessEigenvaluesProcess(ModelPart &rModelPart,
                                  const SizeType NumAnimationSteps,
                                  const bool LabelEigenfreq) 
                                  : mrModelPart(rModelPart),
                                    mNumAnimationSteps(NumAnimationSteps),
                                    mLabelEigenfreq(LabelEigenfreq) { }

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

        const SizeType num_nodes = mrModelPart.NumberOfNodes();
        const auto& eigenvalue_vector = mrModelPart.GetProcessInfo()[EIGENVALUE_VECTOR];
        const SizeType num_eigenvalues = eigenvalue_vector.size();

        std::vector<std::string> labels(num_eigenvalues);
        std::vector<SizeType> nodal_ids(num_nodes);
        
        NodeIterator it_begin = mrModelPart.NodesBegin();
        
        const SizeType num_node_dofs_ref = it_begin->GetDofs().size();

        AnimationResults animated_results(mNumAnimationSteps);
        for (SizeType i=0; i<mNumAnimationSteps; ++i)
        {
            animated_results[i].resize(num_eigenvalues);
            for (SizeType j=0; j<num_eigenvalues; ++j)
            {
                animated_results[i][j].resize(num_nodes);
                for (SizeType k=0; k<num_nodes; ++k)
                    animated_results[i][j][k].resize(GetSizeOfVariable(mVariableName));
            }
        }
        // TODO pragma for
        for (SizeType i=0; i<num_nodes; i++)
        {
            NodeIterator node_it = it_begin + i;
            nodal_ids[i] = node_it->Id();
            const SizeType num_node_dofs = node_it->GetDofs().size();
            const auto& eigen_matrix = node_it->GetValue(EIGENVECTOR_MATRIX);
            const SizeType size_eigen_matrix = eigen_matrix.size2();
            
            // Validation
            KRATOS_ERROR_IF_NOT(num_node_dofs == num_node_dofs_ref) 
                << "Inconsistent number of dofs on node " << nodal_ids[i] << std::endl;
            KRATOS_ERROR_IF_NOT(num_node_dofs == size_eigen_matrix) 
                << "Number of results on node " << nodal_ids[i] << " is wrong" << std::endl;

            for (SizeType j=0; j < mNumAnimationSteps; ++j)
            {
                auto& step_results = animated_results[j];

                double angle = 2 * Globals::Pi * j / mNumAnimationSteps;

                for(SizeType k=0; k<num_eigenvalues; ++k)
                {                        
                    auto& nodal_results = step_results[k][i];
                    
                    for (SizeType m=0; m<num_node_dofs_ref; ++m)
                        nodal_results[m] = std::cos(angle) * eigen_matrix(k, m);
                }
            }
        }


        for (SizeType i=0; i<num_eigenvalues; ++i)
        {
            std::stringstream strstr;
            strstr << std::fixed << std::setprecision(3) << eigenvalue_vector[i];
            std::string eigenval = strstr.str();
            







            labels[i] = "EigenValue_" + eigenval;
        }

        for (SizeType i=0; i < mNumAnimationSteps; ++i) // "Time"
        {
            // for (SizeType j=0; j < num_eigenvalues; ++j)
                // mpGidEigenIO->WriteEigenResults(nodal_ids, labels[j], i, animated_results[i][j]);
                // mpGidEigenIO->WriteEigenResults(nodal_ids, labels[j], i, animated_results[i][j]);
        }
    }

    // std::cout << r_node_dofs[m].GetVariable().Name() << std::endl;

// SizeType m = 0;
// // for (auto itDof = std::begin(r_node_dofs); itDof != std::end(r_node_dofs); itDof++)
// for (const auto& dof : r_node_dofs)
// {
//     const std::string variable_name = dof.GetVariable().Name();
//     // std::cout << itDof->GetVariable().Name() << std::endl;
//     if( KratosComponents< Variable<double> >::Has( variable_name ) ) //case of double variable
//     {

//     }
//     else if( KratosComponents< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > >::Has(variable_name) ) //case of component variable
//     {
//         typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > component_type;
//         component_type var_component = KratosComponents< component_type >::Get(variable_name);
//         std::string source_var = var_component.GetSourceVariable().Name();
//     }
//     else
//     {
//         KRATOS_ERROR << std::endl;
//     }
//     ++m;
// }

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
    SizeType mNumAnimationSteps;
    bool mLabelEigenfreq;
    GidEigenIO::Pointer mpGidEigenIO;

    ///@}
    ///@name Private Operators
    ///@{

    struct DofResultLabel
    {
        bool mResultSize;
        std::string mResultLabel;
    };

    ///@}
    ///@name Private Operations
    ///@{

    SizeType GetSizeOfVariable(const std::string VariableName)
    {
        SizeType size_of_variable = 0;
        if( KratosComponents< Variable<double> >::Has( VariableName ) ) //case of double variable
        {
            size_of_variable = 1;
        }
        else if( KratosComponents< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > >::Has(variable_name) ) //case of component variable
        {
            size_of_variable = 3;
        }
        else
        {
            KRATOS_ERROR << "Unsupported Size of Variable " << VariableName << std::endl;
        }
        return size_of_variable;
    }

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
