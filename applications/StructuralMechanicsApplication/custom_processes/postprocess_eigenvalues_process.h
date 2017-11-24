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

    typedef ModelPart::NodeType::DofsContainerType DofsContainerType;
    
    ///@}
    ///@name Life Cycle
    ///@{

    PostprocessEigenvaluesProcess(ModelPart &rModelPart,
                                  Parameters OutputParameters) 
                                  : mrModelPart(rModelPart), 
                                    mOutputParameters(OutputParameters)
    {
        Parameters default_parameters(R"(
            {
                "animation_steps"   :  1,
                "use_eigenfrequency_in_label" : false,
                "eigen_results" : []
            }  )"
        );

        mOutputParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

        mNumAnimationSteps = mOutputParameters["animation_steps"].GetInt();
        mLabelEigenfreq = mOutputParameters["use_eigenfrequency_in_label"].GetBool();
    }

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
                            GiD_PostBinary,
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

        const auto& eigenvalue_vector = mrModelPart.GetProcessInfo()[EIGENVALUE_VECTOR];
        const SizeType num_eigenvalues = eigenvalue_vector.size();

        double angle = 0.0;
        std::string label = "";
        for (SizeType i=0; i < mNumAnimationSteps; ++i)
        {
            angle = 2 * Globals::Pi * i / mNumAnimationSteps;

            for(SizeType j=0; j<num_eigenvalues; ++j)
            {
                for(auto& node : mrModelPart.Nodes())
                {
                    // Copy the eigenvector to the Solutionstepvariable. Credit to Michael Andre
                    DofsContainerType& r_node_dofs = node.GetDofs();
                    Matrix& r_node_eigenvectors = node.GetValue(EIGENVECTOR_MATRIX);

                    KRATOS_ERROR_IF_NOT(r_node_dofs.size() == r_node_eigenvectors.size2()) 
                        << "Number of results on node " << node.Id() << " is wrong" << std::endl;

                    SizeType k = 0;
                    for (auto& r_dof : r_node_dofs)
                        r_dof.GetSolutionStepValue(0) = std::cos(angle) * r_node_eigenvectors(j,k++);
                }
                
                label = GetLabel(eigenvalue_vector[j]);

                const auto& requested_results = mOutputParameters["eigen_results"];

                for (const auto& var_name : requested_results)
                    mpGidEigenIO->WriteEigenResults(mrModelPart, DISPLACEMENT, label, i);
                    // mpGidEigenIO->WriteEigenResults(mrModelPart, GetVariable(var_name.GetString()), label, i);
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
    Parameters mOutputParameters;
    SizeType mNumAnimationSteps;
    bool mLabelEigenfreq;
    GidEigenIO::Pointer mpGidEigenIO;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{
    
    std::string GetLabel(double LabelNumber)
    {
        std::string label = "EigenValue_";

        // Compute the Eigenfrequency from the Eigenvalue (if the User specified it)
        if (mLabelEigenfreq)
        {
            label = "EigenFrequency_";
            LabelNumber = std::sqrt(LabelNumber) / (2 * Globals::Pi);
        }

        std::stringstream strstr;
        strstr << std::fixed << std::setprecision(3) << LabelNumber;
        return label + strstr.str();   
    }


    template<typename T>
    Variable<T> GetVariable(const std::string VariableName)
    {
        return DISPLACEMENT;
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
