//  KratosShapeOptimizationApplication
//
//  License:		 BSD License
//					 license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Long Chen
//                   
//


#if !defined(KRATOS_INPUT_ELEMENTAL_SENSITIVITY_PROCESS_H_INCLUDED )
#define KRATOS_INPUT_ELEMENTAL_SENSITIVITY_PROCESS_H_INCLUDED

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

///@addtogroup ShapeOptimizationApplication
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

/// Process to input primal solution (e.g., displacements) needed for adjoint solution.
/**
 * Requires mesh topology and node numbering to remain constant.
 */
class InputElementalSensitivitiyProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(InputElementalSensitivitiyProcess);

    typedef ModelPart::NodeType NodeType;
    typedef ModelPart::ElementType ElementType;

    typedef hsize_t SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    InputElementalSensitivitiyProcess(ModelPart& rModelPart, Parameters& rParameters)
    : Process(), mrModelPart(rModelPart)
    {
        KRATOS_TRY

        Parameters DefaultParams(R"(
        {
            "model_part_name": "PLEASE_SPECIFY_MODEL_PART",
            "file_name": "PLEASE_SPECIFY_H5_FILE_NAME",
            "variable_list": ["ADJIONT_DISPLACEMENT"]
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

        mNumElements = 0;

        // add nodal solution step variables


        KRATOS_CATCH("")
    }

    /// Destructor.
    ~InputElementalSensitivitiyProcess() override {}

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

    void Execute() override {}

    void ExecuteInitialize() override
    {
        KRATOS_TRY

            /*
        // test for valid input parameters / model part
        for (unsigned int i = 0; i < mVariables.size(); i++)
        {         
            
            KRATOS_WATCH(mrModelPart.GetNodalSolutionStepVariablesList());
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
                    KRATOS_WATCH(mrModelPart.GetNodalSolutionStepVariablesList());
                    KRATOS_THROW_ERROR(std::runtime_error,"variable is not found in nodal solution steps variable list: ", mVariables[i])
                }
            }
            else
            {
                KRATOS_THROW_ERROR(std::invalid_argument, "variable type not supported: ", mVariables[i])
            }
        }

        */
        // elements should be ordered by node id before IO begins

        if (mrModelPart.Elements().IsSorted() == false)
        {
            mrModelPart.Elements().Sort();
        }

        // number of nodes must be constant for IO
        mNumElements = static_cast<SizeType>(mrModelPart.Elements().size());

        KRATOS_CATCH("");

        std::cout << "LCHEN: ExecuteInitialize() finished" << std::endl;
    }

    void ExecuteBeforeSolutionLoop() override
    {
    }

    void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY;

        // open file to read data
        H5::H5File File(mFilename.c_str(), H5F_ACC_RDONLY);

        H5::DataSet ElementIdDataset = File.openDataSet("/ElementalData/Id");

        hsize_t Dims[2];
        ElementIdDataset.getSpace().getSimpleExtentDims(Dims);

        if (mNumElements != Dims[0])
            KRATOS_THROW_ERROR(std::runtime_error, "inconsistent dimension in file: ", mFilename);

        if (mNumElements != static_cast<SizeType> (mrModelPart.Elements().size()))
            KRATOS_THROW_ERROR(std::runtime_error, "detected change in number of elements.", "");

        std::vector<unsigned int> ElementIdBuffer(mrModelPart.Elements().size());
        ElementIdDataset.read(ElementIdBuffer.data(), H5::PredType::NATIVE_UINT);

        // check that the element ordering is consistent between the model part and file
        int Delta = 0;
        auto itElement = std::begin(ElementIdBuffer);
        std::for_each(std::begin(mrModelPart.Elements()), std::end(mrModelPart.Elements()),
            [&](ElementType& rElement)
        {
            Delta += std::abs(static_cast<int>(rElement.Id()) - static_cast<int>(*itElement));
            itElement++;
        });

        if (Delta != 0)
            KRATOS_THROW_ERROR(std::runtime_error, "detected mismatch of node ids in file: ", mFilename);

        // input time step data
        std::stringstream TimeStepPathStream;
        TimeStepPathStream << "/ElementalData/" << std::fixed << std::setprecision(10) << mrModelPart.GetProcessInfo()[TIME];
        std::string TimeStepGroupName = TimeStepPathStream.str();

        for (unsigned int i = 0; i < mVariables.size(); i++)
        {

            std::string DatasetName = TimeStepGroupName + "/" + mVariables[i];
            H5::DataSet VariableDataset = File.openDataSet(DatasetName.c_str());

            if (KratosComponents<Variable<double> >::Has(mVariables[i]))
            {
                KRATOS_WATCH(mVariables[i]);
                std::vector<double> VariableDataBuffer(mrModelPart.Elements().size());
                VariableDataset.read(VariableDataBuffer.data(), H5::PredType::NATIVE_DOUBLE);
                auto itElement = std::begin(VariableDataBuffer);

                std::for_each(std::begin(mrModelPart.Elements()), std::end(mrModelPart.Elements()),
                    [&](ElementType& rElement)
                {
                    std::vector<double> rValues;
                    rValues.resize(3);
                    rValues[0] = *itElement;
                    rValues[1] = *itElement;
                    rValues[2] = *itElement;
                    rElement.SetValue(KratosComponents<Variable<double> >::Get(mVariables[i]), *itElement);
                    //KRATOS_WATCH(KratosComponents<Variable<double> >::Get(mVariables[i]));

                    //rElement.SetValueOnIntegrationPoints(KratosComponents<Variable<double> >::Get(mVariables[i]), rValues,0);
                    //KRATOS_WATCH(rElement);
                    //rElement.SetValueOnIntegrationPoints
                    itElement++;
                });
            }
            else
            {
                KRATOS_THROW_ERROR(std::runtime_error, "this type of variable is not yet implemented!", mFilename);
            }
        }

        KRATOS_CATCH("");

        std::cout << "LCHEN: InputElementalSensitivitiyProcess finished!" << std::endl;
    }

    void ExecuteFinalizeSolutionStep() override
    {
    }

    void ExecuteBeforeOutputStep() override
    {
    }

    void ExecuteAfterOutputStep() override
    {
    }

    void ExecuteFinalize() override
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
    std::string Info() const override
    {
        return "InputElementalSensitivitiyProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << this->Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
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
    SizeType mNumElements;
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
}; /* Class InputElementalSensitivitiyProcess */

///@}

///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
    InputElementalSensitivitiyProcess& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const InputElementalSensitivitiyProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} // Shape Optimization Application group

}  /* namespace Kratos */

#endif /* KRATOS_INPUT_ELEMENTAL_SENSITIVITY_PROCESS_H_INCLUDED defined */
