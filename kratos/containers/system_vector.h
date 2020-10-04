//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//
#if !defined(KRATOS_SYSTEM_VECTOR_H_INCLUDED )
#define  KRATOS_SYSTEM_VECTOR_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "utilities/parallel_utilities.h"
#include "utilities/atomic_utilities.h"

namespace Kratos
{
///@addtogroup ApplicationNameApplication
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

/// Provides a SystemVector which implements FEM assemble capabilities
template<class TDataType=double, class TIndexType=std::size_t>
class SystemVector
{
public:
    ///@name Type Definitions
    ///@{
    typedef TIndexType IndexType;

    /// Pointer definition of SystemVector
    KRATOS_CLASS_POINTER_DEFINITION(SystemVector);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    template< class TGraphType >
    SystemVector(const TGraphType& rGraph){
        mData.resize(rGraph.Size(),false);
    }

    /// Destructor.
    virtual ~SystemVector(){}

    ///@}
    ///@name Operators
    ///@{
    void Clear()
    {
        mData.clear();
    }

    void SetValue(const TDataType value)
    {
        IndexPartition<IndexType>(mData.size()).for_each([&](IndexType i){
            mData[i] = value;
        });
    }

    IndexType Size() const
    {
        return mData.size();
    }

    TDataType& operator()(IndexType I){
        return mData[I];
    }


    const TDataType& operator()(IndexType I) const{
        return mData[I];
    }

    ///@}
    ///@name Operations
    ///@{
    void BeginAssemble(){} //the SMP version does nothing. This function is there to be implemented in the MPI case

    void FinalizeAssemble(){} //the SMP version does nothing. This function is there to be implemented in the MPI case

    template<class TVectorType, class TIndexVectorType >
    void Assemble(
        const TVectorType& rVectorInput,
        const TIndexVectorType& EquationId
    )
    {
        KRATOS_DEBUG_ERROR_IF(rVectorInput.size() != EquationId.size());

        for(unsigned int i=0; i<EquationId.size(); ++i){
            IndexType global_i = EquationId[i];
            KRATOS_DEBUG_ERROR_IF(global_i > mData.size());
            AtomicAdd(mData[global_i] , rVectorInput[i]);
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
std::stringstream buffer;
    buffer << "SystemVector" ;
    return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "SystemVector";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

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
    DenseVector<TDataType> mData;

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
    SystemVector& operator=(SystemVector const& rOther){}

    /// Copy constructor.
    SystemVector(SystemVector const& rOther){}

    ///@}

}; // Class SystemVector

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TDataType, class TIndexType>
inline std::istream& operator >> (std::istream& rIStream,
                SystemVector<TDataType,TIndexType>& rThis)
                {
                    return rIStream;
                }

/// output stream function
template<class TDataType, class TIndexType>
inline std::ostream& operator << (std::ostream& rOStream,
                const SystemVector<TDataType,TIndexType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SYSTEM_VECTOR_H_INCLUDED  defined


