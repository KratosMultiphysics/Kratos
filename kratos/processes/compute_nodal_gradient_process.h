//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//  Colaborator:     Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_COMPUTE_GRADIENT_PROCESS_INCLUDED )
#define  KRATOS_COMPUTE_GRADIENT_PROCESS_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/enums.h"
#include "processes/process.h"
#include "includes/model_part.h"
#include "utilities/geometry_utilities.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > component_type;

///@}
///@name  Enum's
///@{
    
///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Compute Nodal Gradient process
/** This process computes the gradient of a certain variable stored in the nodes
*/

template< int TDim, class TVarType, HistoricalValues THist> 
class ComputeNodalGradientProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ComputeNodalGradientProcess
    KRATOS_CLASS_POINTER_DEFINITION(ComputeNodalGradientProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ComputeNodalGradientProcess(ModelPart& rModelPart
        , TVarType& rOriginVariable
        , Variable<array_1d<double,3> >& rGradientVariable
        , Variable<double>& rAreaVariable);

    /// Destructor.
    ~ComputeNodalGradientProcess() override
    {
    }


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

    void Execute() override
    {
        KRATOS_TRY;
        
        // Set to zero
        ClearGradient();
        
        bounded_matrix<double,TDim+1,  TDim> DN_DX;
        array_1d<double,TDim+1> N;
        double Volume;
        
        #pragma omp parallel for private(DN_DX,  N,  Volume)
        for(int i=0; i<static_cast<int>(mrModelPart.Elements().size()); i++)
        {
            auto it=mrModelPart.ElementsBegin()+i;
            Element::GeometryType& geom = it->GetGeometry();
            GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);
            
            array_1d<double, TDim+1> values;
            for(unsigned int i=0; i<TDim+1; i++)
            {
                values[i] = geom[i].FastGetSolutionStepValue(mrOriginVariable);
            }
            
            const array_1d<double,TDim> grad = prod(trans(DN_DX), values);
            
            for(unsigned int i=0; i<TDim+1; i++)
            {
                for(unsigned int k=0; k<TDim; k++)
                {
                    double& val = GetGradient(geom, i,k);
                    
                    #pragma omp atomic
                    val += N[i]*Volume*grad[k];
                }
                
                double& vol = geom[i].FastGetSolutionStepValue(mrAreaVariable);
                
                #pragma omp atomic
                vol += N[i]*Volume;
            }
        }
        
        PonderateGradient();

        KRATOS_CATCH("")
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
        return "ComputeNodalGradientProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ComputeNodalGradientProcess";
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
    TVarType& mrOriginVariable;
    Variable<array_1d<double,3> >& mrGradientVariable;
    Variable<double>& mrAreaVariable;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    // TODO: Try to use enable_if!!!
    
    /**
     * This clears the gradient
     */
    void ClearGradient();

    /**
     * This gets the gradient value
     * @param geom: The geometry of the element
     * @param i: The node index
     * @param k: The component index
     */
    double& GetGradient(
        Element::GeometryType& geom,
        unsigned int i, 
        unsigned int k
        );
    
    /**
     * This divides the gradient value by the nodal area
     */
    void PonderateGradient();

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
    ComputeNodalGradientProcess& operator=(ComputeNodalGradientProcess const& rOther);

    /// Copy constructor.
    //ComputeNodalGradientProcess(ComputeNodalGradientProcess const& rOther);


    ///@}

}; // Class ComputeNodalGradientProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
// inline std::istream& operator >> (std::istream& rIStream,
//                                   ComputeNodalGradientProcess& rThis);
// 
// /// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
//                                   const ComputeNodalGradientProcess& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);
// 
//     return rOStream;
// }
///@}

///@name Explicit Specializations
///@{

    template<>
    ComputeNodalGradientProcess<2, Variable<double>, Historical>::ComputeNodalGradientProcess(
        ModelPart& rModelPart, 
        Variable<double>& rOriginVariable, 
        Variable<array_1d<double,3> >& rGradientVariable, 
        Variable<double>& rAreaVariable)
        :mrModelPart(rModelPart), mrOriginVariable(rOriginVariable), mrGradientVariable(rGradientVariable), mrAreaVariable(rAreaVariable)
    {
        KRATOS_TRY
        
        if (rModelPart.GetNodalSolutionStepVariablesList().Has( rOriginVariable ) == false )
        {
            KRATOS_ERROR << "Missing variable " << rOriginVariable;
        }
        
        if (rModelPart.GetNodalSolutionStepVariablesList().Has( rGradientVariable ) == false )
        {
            KRATOS_ERROR << "Missing variable " << rGradientVariable;
        }
        
        if (rModelPart.GetNodalSolutionStepVariablesList().Has( rAreaVariable ) == false )
        {
            KRATOS_ERROR << "Missing variable " << rAreaVariable;
        }
        
        KRATOS_CATCH("")
    }
    
    template<>
    ComputeNodalGradientProcess<2, Variable<double>, NonHistorical>::ComputeNodalGradientProcess(
        ModelPart& rModelPart, 
        Variable<double>& rOriginVariable, 
        Variable<array_1d<double,3> >& rGradientVariable, 
        Variable<double>& rAreaVariable)
        :mrModelPart(rModelPart), mrOriginVariable(rOriginVariable), mrGradientVariable(rGradientVariable), mrAreaVariable(rAreaVariable)
    {
        KRATOS_TRY
        
        if (rModelPart.GetNodalSolutionStepVariablesList().Has( rOriginVariable ) == false )
        {
            KRATOS_ERROR << "Missing variable " << rOriginVariable;
        }
        
        if (rModelPart.GetNodalSolutionStepVariablesList().Has( rAreaVariable ) == false )
        {
            KRATOS_ERROR << "Missing variable " << rAreaVariable;
        }
        
        KRATOS_CATCH("")
    }
    
    template<>
    ComputeNodalGradientProcess<3, Variable<double>, Historical>::ComputeNodalGradientProcess(
        ModelPart& rModelPart, 
        Variable<double>& rOriginVariable, 
        Variable<array_1d<double,3> >& rGradientVariable, 
        Variable<double>& rAreaVariable)
        :mrModelPart(rModelPart), mrOriginVariable(rOriginVariable), mrGradientVariable(rGradientVariable), mrAreaVariable(rAreaVariable)
    {
        KRATOS_TRY
        
        if (rModelPart.GetNodalSolutionStepVariablesList().Has( rOriginVariable ) == false )
        {
            KRATOS_ERROR << "Missing variable " << rOriginVariable;
        }
        
        if (rModelPart.GetNodalSolutionStepVariablesList().Has( rGradientVariable ) == false )
        {
            KRATOS_ERROR << "Missing variable " << rGradientVariable;
        }
        
        if (rModelPart.GetNodalSolutionStepVariablesList().Has( rAreaVariable ) == false )
        {
            KRATOS_ERROR << "Missing variable " << rAreaVariable;
        }
        
        KRATOS_CATCH("")
    }
    
    template<>
    ComputeNodalGradientProcess<3, Variable<double>, NonHistorical>::ComputeNodalGradientProcess(
        ModelPart& rModelPart, 
        Variable<double>& rOriginVariable, 
        Variable<array_1d<double,3> >& rGradientVariable, 
        Variable<double>& rAreaVariable)
        :mrModelPart(rModelPart), mrOriginVariable(rOriginVariable), mrGradientVariable(rGradientVariable), mrAreaVariable(rAreaVariable)
    {
        KRATOS_TRY
        
        if (rModelPart.GetNodalSolutionStepVariablesList().Has( rOriginVariable ) == false )
        {
            KRATOS_ERROR << "Missing variable " << rOriginVariable;
        }
        
        if (rModelPart.GetNodalSolutionStepVariablesList().Has( rAreaVariable ) == false )
        {
            KRATOS_ERROR << "Missing variable " << rAreaVariable;
        }
        
        KRATOS_CATCH("")
    }
    
    template<>
    ComputeNodalGradientProcess<2, component_type, Historical>::ComputeNodalGradientProcess(
        ModelPart& rModelPart, 
        component_type& rOriginVariable, 
        Variable<array_1d<double,3> >& rGradientVariable, 
        Variable<double>& rAreaVariable)
        :mrModelPart(rModelPart), mrOriginVariable(rOriginVariable), mrGradientVariable(rGradientVariable), mrAreaVariable(rAreaVariable)
    {
        KRATOS_TRY
        
        if (rModelPart.GetNodalSolutionStepVariablesList().Has( rOriginVariable.GetSourceVariable() ) == false )
        {
            KRATOS_ERROR << "Missing variable " << rOriginVariable.GetSourceVariable();
        }
        
        if (rModelPart.GetNodalSolutionStepVariablesList().Has( rGradientVariable ) == false )
        {
            KRATOS_ERROR << "Missing variable " << rGradientVariable;
        }
        
        if (rModelPart.GetNodalSolutionStepVariablesList().Has( rAreaVariable ) == false )
        {
            KRATOS_ERROR << "Missing variable " << rAreaVariable;
        }
        
        KRATOS_CATCH("")
    }
    
    template<>
    ComputeNodalGradientProcess<2, component_type, NonHistorical>::ComputeNodalGradientProcess(
        ModelPart& rModelPart, 
        component_type& rOriginVariable, 
        Variable<array_1d<double,3> >& rGradientVariable, 
        Variable<double>& rAreaVariable)
        :mrModelPart(rModelPart), mrOriginVariable(rOriginVariable), mrGradientVariable(rGradientVariable), mrAreaVariable(rAreaVariable)
    {
        KRATOS_TRY
        
        if (rModelPart.GetNodalSolutionStepVariablesList().Has( rOriginVariable.GetSourceVariable() ) == false )
        {
            KRATOS_ERROR << "Missing variable " << rOriginVariable.GetSourceVariable();
        }
        
        if (rModelPart.GetNodalSolutionStepVariablesList().Has( rAreaVariable ) == false )
        {
            KRATOS_ERROR << "Missing variable " << rAreaVariable;
        }
        
        KRATOS_CATCH("")
    }
    
    template<>
    ComputeNodalGradientProcess<3, component_type, Historical>::ComputeNodalGradientProcess(
        ModelPart& rModelPart, 
        component_type& rOriginVariable, 
        Variable<array_1d<double,3> >& rGradientVariable, 
        Variable<double>& rAreaVariable)
        :mrModelPart(rModelPart), mrOriginVariable(rOriginVariable), mrGradientVariable(rGradientVariable), mrAreaVariable(rAreaVariable)
    {
        KRATOS_TRY
        
        if (rModelPart.GetNodalSolutionStepVariablesList().Has( rOriginVariable.GetSourceVariable() ) == false )
        {
            KRATOS_ERROR << "Missing variable " << rOriginVariable.GetSourceVariable();
        }
        
        if (rModelPart.GetNodalSolutionStepVariablesList().Has( rGradientVariable ) == false )
        {
            KRATOS_ERROR << "Missing variable " << rGradientVariable;
        }
        
        if (rModelPart.GetNodalSolutionStepVariablesList().Has( rAreaVariable ) == false )
        {
            KRATOS_ERROR << "Missing variable " << rAreaVariable;
        }
        
        KRATOS_CATCH("")
    }
    
    template<>
    ComputeNodalGradientProcess<3, component_type, NonHistorical>::ComputeNodalGradientProcess(
        ModelPart& rModelPart, 
        component_type& rOriginVariable, 
        Variable<array_1d<double,3> >& rGradientVariable, 
        Variable<double>& rAreaVariable)
        :mrModelPart(rModelPart), mrOriginVariable(rOriginVariable), mrGradientVariable(rGradientVariable), mrAreaVariable(rAreaVariable)
    {
        KRATOS_TRY
        
        if (rModelPart.GetNodalSolutionStepVariablesList().Has( rOriginVariable.GetSourceVariable() ) == false )
        {
            KRATOS_ERROR << "Missing variable " << rOriginVariable.GetSourceVariable();
        }
        
        if (rModelPart.GetNodalSolutionStepVariablesList().Has( rAreaVariable ) == false )
        {
            KRATOS_ERROR << "Missing variable " << rAreaVariable;
        }
        
        KRATOS_CATCH("")
    }

    template<>
    void ComputeNodalGradientProcess<2, Variable<double>, Historical>::ClearGradient()
    {
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(mrModelPart.Nodes().size()); i++)
        {
            auto it=mrModelPart.NodesBegin()+i;
            it->FastGetSolutionStepValue(mrAreaVariable) = 0.0;
            it->FastGetSolutionStepValue(mrGradientVariable).clear();
        }
    }
    
    template<>
    void ComputeNodalGradientProcess<3, Variable<double>, Historical>::ClearGradient()
    {
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(mrModelPart.Nodes().size()); i++)
        {
            auto it=mrModelPart.NodesBegin()+i;
            it->FastGetSolutionStepValue(mrAreaVariable) = 0.0;
            it->FastGetSolutionStepValue(mrGradientVariable).clear();
        }
    }
    
    template<>
    void ComputeNodalGradientProcess<2, component_type, Historical>::ClearGradient()
    {
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(mrModelPart.Nodes().size()); i++)
        {
            auto it=mrModelPart.NodesBegin()+i;
            it->FastGetSolutionStepValue(mrAreaVariable) = 0.0;
            it->FastGetSolutionStepValue(mrGradientVariable).clear();
        }
    }
    
    template<>
    void ComputeNodalGradientProcess<3, component_type, Historical>::ClearGradient()
    {
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(mrModelPart.Nodes().size()); i++)
        {
            auto it=mrModelPart.NodesBegin()+i;
            it->FastGetSolutionStepValue(mrAreaVariable) = 0.0;
            it->FastGetSolutionStepValue(mrGradientVariable).clear();
        }
    }

    template <>
    void ComputeNodalGradientProcess<2, Variable<double>, NonHistorical>::ClearGradient()
    {
        const array_1d<double, 3> AuxZeroVector = ZeroVector(3);
        
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(mrModelPart.Nodes().size()); i++)
        {
            auto it=mrModelPart.NodesBegin()+i;
            it->FastGetSolutionStepValue(mrAreaVariable) = 0.0;
            it->SetValue(mrGradientVariable, AuxZeroVector);
        }
    }

    template <>
    void ComputeNodalGradientProcess<3, Variable<double>, NonHistorical>::ClearGradient()
    {
        const array_1d<double, 3> AuxZeroVector = ZeroVector(3);
        
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(mrModelPart.Nodes().size()); i++)
        {
            auto it=mrModelPart.NodesBegin()+i;
            it->FastGetSolutionStepValue(mrAreaVariable) = 0.0;
            it->SetValue(mrGradientVariable, AuxZeroVector);
        }
    }

    template <>
    void ComputeNodalGradientProcess<2, component_type, NonHistorical>::ClearGradient()
    {
        const array_1d<double, 3> AuxZeroVector = ZeroVector(3);
        
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(mrModelPart.Nodes().size()); i++)
        {
            auto it=mrModelPart.NodesBegin()+i;
            it->FastGetSolutionStepValue(mrAreaVariable) = 0.0;
            it->SetValue(mrGradientVariable, AuxZeroVector);
        }
    }

    template <>
    void ComputeNodalGradientProcess<3, component_type, NonHistorical>::ClearGradient()
    {
        const array_1d<double, 3> AuxZeroVector = ZeroVector(3);
        
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(mrModelPart.Nodes().size()); i++)
        {
            auto it=mrModelPart.NodesBegin()+i;
            it->FastGetSolutionStepValue(mrAreaVariable) = 0.0;
            it->SetValue(mrGradientVariable, AuxZeroVector);
        }
    }

    template <>
    double& ComputeNodalGradientProcess<2, Variable<double>, Historical>::GetGradient(
        Element::GeometryType& geom,
        unsigned int i, 
        unsigned int k
        )
    {
        double& val = geom[i].FastGetSolutionStepValue(mrGradientVariable)[k];
        
        return val;
    }

    template <>
    double& ComputeNodalGradientProcess<3, Variable<double>, Historical>::GetGradient(
        Element::GeometryType& geom,
        unsigned int i, 
        unsigned int k
        )
    {
        double& val = geom[i].FastGetSolutionStepValue(mrGradientVariable)[k];
        
        return val;
    }
    
    template <>
    double& ComputeNodalGradientProcess<2, component_type, Historical>::GetGradient(
        Element::GeometryType& geom,
        unsigned int i, 
        unsigned int k
        )
    {
        double& val = geom[i].FastGetSolutionStepValue(mrGradientVariable)[k];
        
        return val;
    }

    template <>
    double& ComputeNodalGradientProcess<3, component_type, Historical>::GetGradient(
        Element::GeometryType& geom,
        unsigned int i, 
        unsigned int k
        )
    {
        double& val = geom[i].FastGetSolutionStepValue(mrGradientVariable)[k];
        
        return val;
    }

    template <>
    double& ComputeNodalGradientProcess<2, Variable<double>, NonHistorical>::GetGradient(
        Element::GeometryType& geom,
        unsigned int i, 
        unsigned int k
        )
    {
        double& val = geom[i].GetValue(mrGradientVariable)[k];
        
        return val;
    }
    
    template <>
    double& ComputeNodalGradientProcess<3, Variable<double>, NonHistorical>::GetGradient(
        Element::GeometryType& geom,
        unsigned int i, 
        unsigned int k
        )
    {
        double& val = geom[i].GetValue(mrGradientVariable)[k];
        
        return val;
    }

    template <>
    double& ComputeNodalGradientProcess<2, component_type, NonHistorical>::GetGradient(
        Element::GeometryType& geom,
        unsigned int i, 
        unsigned int k
        )
    {
        double& val = geom[i].GetValue(mrGradientVariable)[k];
        
        return val;
    }
    
    template <>
    double& ComputeNodalGradientProcess<3, component_type, NonHistorical>::GetGradient(
        Element::GeometryType& geom,
        unsigned int i, 
        unsigned int k
        )
    {
        double& val = geom[i].GetValue(mrGradientVariable)[k];
        
        return val;
    }

    template <>
    void ComputeNodalGradientProcess<2, Variable<double>, Historical>::PonderateGradient()
    {
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(mrModelPart.Nodes().size()); i++)
        {
            auto it=mrModelPart.NodesBegin()+i;
            it->FastGetSolutionStepValue(mrGradientVariable) /= it->FastGetSolutionStepValue(mrAreaVariable);
        }
    }
    
    template <>
    void ComputeNodalGradientProcess<3, Variable<double>, Historical>::PonderateGradient()
    {
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(mrModelPart.Nodes().size()); i++)
        {
            auto it=mrModelPart.NodesBegin()+i;
            it->FastGetSolutionStepValue(mrGradientVariable) /= it->FastGetSolutionStepValue(mrAreaVariable);
        }
    }

    template <>
    void ComputeNodalGradientProcess<2, component_type, Historical>::PonderateGradient()
    {
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(mrModelPart.Nodes().size()); i++)
        {
            auto it=mrModelPart.NodesBegin()+i;
            it->FastGetSolutionStepValue(mrGradientVariable) /= it->FastGetSolutionStepValue(mrAreaVariable);
        }
    }
    
    template <>
    void ComputeNodalGradientProcess<3, component_type, Historical>::PonderateGradient()
    {
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(mrModelPart.Nodes().size()); i++)
        {
            auto it=mrModelPart.NodesBegin()+i;
            it->FastGetSolutionStepValue(mrGradientVariable) /= it->FastGetSolutionStepValue(mrAreaVariable);
        }
    }

    template <>
    void ComputeNodalGradientProcess<2, Variable<double>, NonHistorical>::PonderateGradient()
    {
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(mrModelPart.Nodes().size()); i++)
        {
            auto it=mrModelPart.NodesBegin()+i;
            it->GetValue(mrGradientVariable) /= it->FastGetSolutionStepValue(mrAreaVariable);
        }
    }

    template <>
    void ComputeNodalGradientProcess<3, Variable<double>, NonHistorical>::PonderateGradient()
    {
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(mrModelPart.Nodes().size()); i++)
        {
            auto it=mrModelPart.NodesBegin()+i;
            it->GetValue(mrGradientVariable) /= it->FastGetSolutionStepValue(mrAreaVariable);
        }
    }
    
    template <>
    void ComputeNodalGradientProcess<2, component_type, NonHistorical>::PonderateGradient()
    {
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(mrModelPart.Nodes().size()); i++)
        {
            auto it=mrModelPart.NodesBegin()+i;
            it->GetValue(mrGradientVariable) /= it->FastGetSolutionStepValue(mrAreaVariable);
        }
    }

    template <>
    void ComputeNodalGradientProcess<3, component_type, NonHistorical>::PonderateGradient()
    {
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(mrModelPart.Nodes().size()); i++)
        {
            auto it=mrModelPart.NodesBegin()+i;
            it->GetValue(mrGradientVariable) /= it->FastGetSolutionStepValue(mrAreaVariable);
        }
    }
    
}  // namespace Kratos.

#endif // KRATOS_COMPUTE_GRADIENT_PROCESS_INCLUDED  defined 


