//  KratosAdjointFluidApplication
//
//  License:		 BSD License
//					 license: AdjointFluidApplication/license.txt
//
//  Main authors:    Michael Andre, https://github.com/msandre
//

#if !defined(KRATOS_INPUT_PRIMAL_SOLUTION_PROCESS_H_INCLUDED )
#define KRATOS_INPUT_PRIMAL_SOLUTION_PROCESS_H_INCLUDED

// System includes
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <iomanip>

// External includes
#include "H5Cpp.h"

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "includes/model_part.h"

namespace Kratos
{

///@addtogroup AdjointFluidApplication
///@{

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

/// Process to input primal solution (e.g., velocity and pressure) needed for adjoint solution.
/**
 * Requires mesh topology and node numbering to remain constant.
 */
class InputPrimalSolutionProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(InputPrimalSolutionProcess);

    typedef ModelPart::NodeType NodeType;

    typedef hsize_t SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    InputPrimalSolutionProcess(ModelPart& rModelPart, Parameters& rParameters)
    : Process(), mrModelPart(rModelPart)
    {
        KRATOS_TRY

        Parameters DefaultParams(R"(
        {
            "model_part_name": "PLEASE_SPECIFY_MODEL_PART",
            "file_name": "PLEASE_SPECIFY_H5_FILE_NAME",
            "variable_list": ["VELOCITY", "ACCELERATION", "PRESSURE"]
        })");

        // try accessing parameters without defaults so that an error is thrown
        // if they don't exist
        rParameters["file_name"];

        rParameters.ValidateAndAssignDefaults(DefaultParams);

        // input filename
        std::stringstream FilenameStream;
        FilenameStream << rParameters["file_name"].GetString() << '_' << mrModelPart.GetCommunicator().MyPID() << ".h5";
        mFilename = FilenameStream.str();

        // nodal variable names to input
        mVariables.resize(rParameters["variable_list"].size());
        for (unsigned int i = 0; i < mVariables.size(); i++)
        {
            mVariables[i] = rParameters["variable_list"].GetArrayItem(i).GetString();
        }

        mNumNodes = 0;

        KRATOS_CATCH("")
    }

    /// Destructor.
    virtual ~InputPrimalSolutionProcess() {}

    ///@}
    ///@name Operators
    ///@{

    /// This operator simply calls the Execute method.
    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    virtual void Execute() {}

    virtual void ExecuteInitialize()
    {
        KRATOS_TRY

        // test for valid input parameters / model part
        for (unsigned int i = 0; i < mVariables.size(); i++)
        {
            if (KratosComponents< Variable<double> >::Has(mVariables[i]))
            {
                if (mrModelPart.GetNodalSolutionStepVariablesList().Has(KratosComponents< Variable<double> >::Get(mVariables[i])) == false)
                {
                    KRATOS_THROW_ERROR(std::runtime_error, "variable is not found in nodal solution steps variable list: ", mVariables[i])
                }
            }
            else if (KratosComponents< Variable<array_1d<double,3> > >::Has(mVariables[i]))
            {
                if (mrModelPart.GetNodalSolutionStepVariablesList().Has(KratosComponents< Variable<array_1d<double,3> > >::Get(mVariables[i])) == false)
                {
                    KRATOS_THROW_ERROR(std::runtime_error,"variable is not found in nodal solution steps variable list: ", mVariables[i])
                }
            }
            else
            {
                KRATOS_THROW_ERROR(std::invalid_argument, "variable type not supported: ", mVariables[i])
            }
        }

        // nodes should be ordered by node id before IO begins
        if (mrModelPart.Nodes().IsSorted() == false)
        {
            mrModelPart.Nodes().Sort();
        }

        // number of nodes must be constant for IO
        mNumNodes = static_cast<SizeType>(mrModelPart.Nodes().size());

        KRATOS_CATCH("")
    }

    virtual void ExecuteBeforeSolutionLoop()
    {
    }

    virtual void ExecuteInitializeSolutionStep()
    {
        KRATOS_TRY

        // open file to read data
        H5::H5File File(mFilename.c_str(), H5F_ACC_RDONLY);

        H5::DataSet NodeIdDataset = File.openDataSet("/NodalData/Id");

        hsize_t Dims[2];
        NodeIdDataset.getSpace().getSimpleExtentDims(Dims);

        if (mNumNodes != Dims[0])
            KRATOS_THROW_ERROR(std::runtime_error, "inconsistent dimension in file: ", mFilename)

        if (mNumNodes != static_cast<SizeType>(mrModelPart.Nodes().size()))
            KRATOS_THROW_ERROR(std::runtime_error, "detected change in number of nodes.", "")

        std::vector<unsigned int> NodeIdBuffer(mrModelPart.Nodes().size());
        NodeIdDataset.read(NodeIdBuffer.data(), H5::PredType::NATIVE_UINT);

        // check that the node ordering is consistent between the model part and file
        int Delta = 0;
        auto itNode = std::begin(NodeIdBuffer);
        std::for_each(std::begin(mrModelPart.Nodes()), std::end(mrModelPart.Nodes()),
           [&](NodeType& rNode)
           {
                Delta += std::abs(static_cast<int>(rNode.Id()) - static_cast<int>(*itNode));
                itNode++;
           });

        if (Delta != 0)
            KRATOS_THROW_ERROR(std::runtime_error, "detected mismatch of node ids in file: ", mFilename);

        // input time step data
        std::stringstream TimeStepPathStream;
        TimeStepPathStream << "/NodalData/" << std::fixed << std::setprecision(10) << mrModelPart.GetProcessInfo()[TIME];
        std::string TimeStepGroupName = TimeStepPathStream.str();

        for (unsigned int i = 0; i < mVariables.size(); i++)
        {
            std::string DatasetName = TimeStepGroupName + "/" + mVariables[i];
            H5::DataSet VariableDataset = File.openDataSet(DatasetName.c_str());
            if (KratosComponents< Variable<double> >::Has(mVariables[i]))
            {
                std::vector<double> VariableDataBuffer(mrModelPart.Nodes().size());
                VariableDataset.read(VariableDataBuffer.data(), H5::PredType::NATIVE_DOUBLE);
                auto itNode = std::begin(VariableDataBuffer);
                std::for_each(std::begin(mrModelPart.Nodes()), std::end(mrModelPart.Nodes()),
                   [&](NodeType& rNode)
                   {
                        rNode.FastGetSolutionStepValue(KratosComponents< Variable<double> >::Get(mVariables[i])) = *itNode;
                        itNode++;
                   });
            }
            else if (KratosComponents< Variable<array_1d<double,3> > >::Has(mVariables[i]))
            {
                std::vector<double> VariableDataBuffer(3 * mrModelPart.Nodes().size());
                VariableDataset.read(VariableDataBuffer.data(), H5::PredType::NATIVE_DOUBLE);
                unsigned int BlockBegin = 0;
                for (auto it = std::begin(mrModelPart.Nodes()); it != std::end(mrModelPart.Nodes()); it++)
                {
                    array_1d<double,3>& rData = it->FastGetSolutionStepValue(KratosComponents< Variable<array_1d<double,3> > >::Get(mVariables[i]));
                    for (unsigned int k = 0; k < 3; k++)
                        rData[k] = VariableDataBuffer[BlockBegin + k];
                    BlockBegin += 3;
                }
            }
        }

        KRATOS_CATCH("")
    }

    virtual void ExecuteFinalizeSolutionStep()
    {
    }

    virtual void ExecuteBeforeOutputStep()
    {
    }

    virtual void ExecuteAfterOutputStep()
    {
    }

    virtual void ExecuteFinalize()
    {
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
        return "InputPrimalSolutionProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << this->Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
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

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;
    std::string mFilename;
    SizeType mNumNodes;
    std::vector<std::string> mVariables;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
}; /* Class InputPrimalSolutionProcess */

///@}

///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  InputPrimalSolutionProcess& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const InputPrimalSolutionProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} // Adjoint Fluid Application group

}  /* namespace Kratos */

#endif /* KRATOS_INPUT_PRIMAL_SOLUTION_PROCESS_H_INCLUDED defined */
