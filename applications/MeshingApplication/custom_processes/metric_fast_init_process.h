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

#if !defined(KRATOS_METRIC_FAST_INIT_PROCESS)
#define KRATOS_METRIC_FAST_INIT_PROCESS

// System includes

// External includes

// Project includes
#include "processes/process.h"
#include "meshing_application.h"
#include "includes/model_part.h"

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

/// Short class definition.
// This process initializes the variables related with the ALM
/** Detail class definition.
*/
template<unsigned int TDim>
class MetricFastInit
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MetricFastInit
    KRATOS_CLASS_POINTER_DEFINITION(MetricFastInit);

    // General type definitions
    typedef Node<3>                                          NodeType;
    typedef Geometry<NodeType>                           GeometryType;
    typedef ModelPart::NodesContainerType              NodesArrayType;
    typedef ModelPart::ConditionsContainerType    ConditionsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MetricFastInit( ModelPart& rThisModelPart):mrThisModelPart(rThisModelPart)
    {
        KRATOS_TRY;

        KRATOS_CATCH("");
    }

    /// Destructor.
    virtual ~MetricFastInit()
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

    ///@}
    ///@name Friends
    ///@{

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

    virtual void Execute() override
    {
        KRATOS_TRY;

        constexpr unsigned int size = TDim == 2  ? 3: 6;

        const array_1d<double, size> zerovector(size, 0.0);

        // We iterate over the node
        NodesArrayType& pNodes = mrThisModelPart.Nodes();
        int numNodes = mrThisModelPart.NumberOfNodes();

        #pragma omp parallel for firstprivate(zerovector)
        for(int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNodes.begin() + i;

            // The metric
            itNode->SetValue(MMG_METRIC, zerovector);
        }

        KRATOS_CATCH("");
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
    virtual std::string Info() const override
    {
        return "MetricFastInit";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MetricFastInit";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
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
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{
    ModelPart& mrThisModelPart;

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

    /// Assignment operator.
    MetricFastInit& operator=(MetricFastInit const& rOther);

    /// Copy constructor.
    //MetricFastInit(MetricFastInit const& rOther);


    ///@}

}; // Class MetricFastInit

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
// inline std::istream& operator >> (std::istream& rIStream,
//                                   MetricFastInit& rThis);
//
// /// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
//                                   const MetricFastInit& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);
//
//     return rOStream;
// }

}
#endif /* KRATOS_METRIC_FAST_INIT_PROCESS defined */
