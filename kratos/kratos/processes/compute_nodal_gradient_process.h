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
//  Colaborator:     Vicente Mataix Ferrándiz
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
//erases the nodes marked as
/** Detail class definition.
*/

template< int TDim > // TODO: Consider the number of nodes, this way will be possible to compute for any geometry, not just triangles and tetrahedra
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
    ComputeNodalGradientProcess(ModelPart& model_part
        , Variable<double>& r_origin_variable
        , Variable<array_1d<double,3> >& r_gradient_variable
        , Variable<double>& r_area_variable)
        : mr_model_part(model_part), mr_origin_variable(r_origin_variable), mr_gradient_variable(r_gradient_variable), mr_area_variable(r_area_variable)
    {
        KRATOS_TRY
        
        if (model_part.GetNodalSolutionStepVariablesList().Has( r_origin_variable ) == false )
        {
            KRATOS_ERROR << "missing variable " << r_origin_variable;
        }
        
        
        if (model_part.GetNodalSolutionStepVariablesList().Has( r_gradient_variable ) == false )
        {
            KRATOS_ERROR << "missing variable " << r_gradient_variable;
        }
        
        if (model_part.GetNodalSolutionStepVariablesList().Has( r_area_variable ) == false )
        {
            KRATOS_ERROR << "missing variable " << r_area_variable;
        }
        
        KRATOS_CATCH("")
    }

    /// Destructor.
    virtual ~ComputeNodalGradientProcess()
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

    virtual void Execute()
    {
        KRATOS_TRY;
        
        //set to zero
        
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(mr_model_part.Nodes().size()); i++)
        {
            auto it=mr_model_part.NodesBegin()+i;
            it->FastGetSolutionStepValue(mr_area_variable) = 0.0;
            it->FastGetSolutionStepValue(mr_gradient_variable).clear();
        }
        
        boost::numeric::ublas::bounded_matrix<double,TDim+1,  TDim> DN_DX;
        array_1d<double,TDim+1> N;
        double Volume;
        
        #pragma omp parallel for private(DN_DX,  N,  Volume)
        for(int i=0; i<static_cast<int>(mr_model_part.Elements().size()); i++)
        {
            auto it=mr_model_part.ElementsBegin()+i;
            Element::GeometryType& geom = it->GetGeometry();
            GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);
            
            array_1d<double, TDim+1> values;
            for(unsigned int i=0; i<TDim+1; i++)
            {
                values[i] = geom[i].FastGetSolutionStepValue(mr_origin_variable);
            }
            
            const array_1d<double,TDim> grad = prod(trans(DN_DX), values);
            
            for(unsigned int i=0; i<TDim+1; i++)
            {
                for(unsigned int k=0; k<TDim; k++)
                {
                    double& val = geom[i].FastGetSolutionStepValue(mr_gradient_variable)[k];
                    
                    #pragma omp atomic
                    val += N[i]*Volume*grad[k];
                }
                
                double& vol = geom[i].FastGetSolutionStepValue(mr_area_variable);
                
                #pragma omp atomic
                vol += N[i]*Volume;
            }
        }
        
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(mr_model_part.Nodes().size()); i++)
        {
            auto it=mr_model_part.NodesBegin()+i;
            it->FastGetSolutionStepValue(mr_gradient_variable) /= it->FastGetSolutionStepValue(mr_area_variable);
        }

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
    virtual std::string Info() const
    {
        return "ComputeNodalGradientProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ComputeNodalGradientProcess";
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
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{
    ModelPart& mr_model_part;
    Variable<double>& mr_origin_variable;
    Variable<array_1d<double,3> >& mr_gradient_variable;
    Variable<double>& mr_area_variable;


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


}  // namespace Kratos.

#endif // KRATOS_COMPUTE_GRADIENT_PROCESS_INCLUDED  defined 


