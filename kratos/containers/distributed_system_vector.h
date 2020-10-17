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
#if !defined(KRATOS_DISTRIBUTED_SYSTEM_VECTOR_H_INCLUDED )
#define  KRATOS_DISTRIBUTED_SYSTEM_VECTOR_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "containers/distributed_system_vector.h"
#include "containers/distributed_sparse_graph.h"
#include "containers/distributed_numbering.h"
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

/// Provides a DistributedSystemVector which implements FEM assemble capabilities
template<class TDataType=double, class TIndexType=std::size_t>
class DistributedSystemVector
{
public:
    ///@name Type Definitions
    ///@{
    typedef TIndexType IndexType;
    typedef int MpiIndexType;


    /// Pointer definition of DistributedSystemVector
    KRATOS_CLASS_POINTER_DEFINITION(DistributedSystemVector);

    ///@}
    ///@name Life Cycle
    ///@{

    DistributedSystemVector(const DistributedSparseGraph<IndexType>& rGraph)
            :
            mrComm(rGraph.GetComm())
    {
        mpNumbering = Kratos::make_unique< DistributedNumbering<IndexType> >( rGraph.GetRowNumbering());

        mLocalData.resize(rGraph.LocalSize(),false);
        mNonLocalData.resize(mrComm.Size());

        //now add entries in nonlocal data;
        const auto& r_non_local_graphs = rGraph.GetNonLocalGraphs();

        IndexPartition<IndexType>(r_non_local_graphs.size()).for_each([&](IndexType cpu_id)
        {
            const auto& graph = r_non_local_graphs[cpu_id];

            for(auto item=graph.begin(); item!=graph.end(); ++item)
            {
                    IndexType row = item.GetRowIndex(); //note that this is a "remote local id"
                    mNonLocalData[cpu_id][row] = TDataType(); //first touching of nonlocaldata
                //}
            }
        });

        //compute communication colors
        const auto& nonlocal_graphs = rGraph.GetNonLocalGraphs();
        std::vector<int> send_list;
        for(unsigned int cpu_id = 0; cpu_id<nonlocal_graphs.size(); ++cpu_id)
                if( !nonlocal_graphs[cpu_id].IsEmpty())
                    send_list.push_back(cpu_id);

        mfem_assemble_colors = MPIColoringUtilities::ComputeCommunicationScheduling(send_list, GetComm());

        //fill the list of indices to be received
        for(auto color : mfem_assemble_colors)
        {
            if(color >= 0) //-1 would imply no communication
            {
                //compute the list of ids to send
                std::vector<IndexType> send_ids;
                for(const auto& item : mNonLocalData[color])
                    send_ids.push_back(item.first);
                mRecvIndicesByColor[color] = GetComm().SendRecv(send_ids, color, color); 
            }
        }

    }

    /// Copy constructor.
    explicit DistributedSystemVector(DistributedSystemVector const& rOther)
    :
    mrComm(rOther.mrComm)
    {
        mpNumbering = Kratos::make_unique<DistributedNumbering<IndexType>>(rOther.GetNumbering());

        //copying the data
        mLocalData.resize(rOther.LocalSize(),false); 
        IndexPartition<IndexType>(LocalSize()).for_each([&](IndexType i){
            (mLocalData)[i] += rOther[i];
        });

        mNonLocalData = rOther.mNonLocalData;
        mfem_assemble_colors = rOther.mfem_assemble_colors;
        mRecvIndicesByColor = rOther.mRecvIndicesByColor;
    }

    DistributedSystemVector(const DistributedNumbering<IndexType>& rNumbering)
            :
            mrComm(rNumbering.GetComm())
    {
        mpNumbering = Kratos::make_unique< DistributedNumbering<IndexType> >( rNumbering );

        mLocalData.resize(rNumbering.LocalSize(),false);
        mNonLocalData.resize(mrComm.Size());
    }

    /// Destructor.
    virtual ~DistributedSystemVector(){}

    ///@}
    ///@name Operators
    ///@{
    const DataCommunicator& GetComm(){
        return mrComm;
    }

    const DistributedNumbering<IndexType>& GetNumbering() const
    {
        return *mpNumbering;
    }

    void Clear()
    {
        mLocalData.clear();
    }

    void SetValue(const TDataType value)
    {
        IndexPartition<IndexType>(mLocalData.size()).for_each([&](IndexType i){
            mLocalData[i] = value;
        });
    }

    TDataType& operator()(IndexType I){
        return mLocalData[I];
    }


    const TDataType& operator()(IndexType I) const{
        return mLocalData[I];
    }

    TDataType& operator[](IndexType I){
        return mLocalData[I];
    }

    const TDataType& operator[](IndexType I) const{
        return mLocalData[I];
    }

    inline IndexType TotalSize() const{ //TODO discuss if this shall be called simply "Size"
        return mpNumbering->TotalSize(); 
    }

    inline IndexType LocalSize() const{
        return mpNumbering->LocalSize(); 
    }

    DenseVector<TDataType>& GetLocalData(){
        return mLocalData;
    }

    const DenseVector<TDataType>& GetLocalData() const{
        return mLocalData;
    }

    void Add(const double factor,
             const DistributedSystemVector& rOtherVector
            )
    {
        IndexPartition<IndexType>(LocalSize()).for_each([&](IndexType i){
            (mLocalData)[i] += factor*rOtherVector[i];
        });
    }

    /// Assignment operator.
    DistributedSystemVector& operator=(DistributedSystemVector const& rOtherVector){
        mLocalData.resize(rOtherVector.LocalSize(),false);
        IndexPartition<IndexType>(LocalSize()).for_each([&](IndexType i){
            (mLocalData)[i] = rOtherVector[i];
        });
    }


    void operator+=(const DistributedSystemVector& rOtherVector)
    {
        IndexPartition<IndexType>(LocalSize()).for_each([&](IndexType i){
            (mLocalData)[i] += rOtherVector[i];
        });
    }

    void operator-=(const DistributedSystemVector& rOtherVector)
    {
        IndexPartition<IndexType>(LocalSize()).for_each([&](IndexType i){
            (mLocalData)[i] -= rOtherVector[i];
        });
    }

    void operator*=(const TDataType& multiplier_factor)
    {
        IndexPartition<IndexType>(LocalSize()).for_each([&](IndexType i){
            (mLocalData)[i] *= multiplier_factor;
        });
    }

    void operator/=(const TDataType& divide_factor)
    {
        IndexPartition<IndexType>(LocalSize()).for_each([&](IndexType i){
            (mLocalData)[i] /= divide_factor;
        });
    }

    ///@}
    ///@name Operations
    ///@{
    void BeginAssemble(){
        //TODO set to zero nonlocal data prior to assembly
        IndexPartition<IndexType>(mNonLocalData.size()).for_each([&](IndexType cpu_id){
            for(auto & item : this->mNonLocalData[cpu_id])
            {
                item.second = 0.0;
            }
        });

    } 

    //communicate data to finalize the assembly
    void FinalizeAssemble(){

        std::vector<TDataType> send_buffer;
        std::vector<TDataType> recv_buffer;

        //sendrecv data
        for(auto color : mfem_assemble_colors)
        {
            if(color >= 0) //-1 would imply no communication
            {
                const auto& recv_indices = mRecvIndicesByColor[color];
                send_buffer.resize(mNonLocalData[color].size());
                recv_buffer.resize(mRecvIndicesByColor[color].size());

                IndexType counter = 0;
                for(const auto& it : mNonLocalData[color])
                {
                    send_buffer[counter++] = it.second;
                }

                // //NOTE: this can be made nonblocking 
                GetComm().SendRecv(send_buffer, color, 0, recv_buffer, color, 0); //TODO, we know all the sizes, we shall use that!

                for(IndexType i=0; i<recv_buffer.size(); ++i)
                {
                    IndexType local_i = recv_indices[i];

                    AtomicAdd( mLocalData[local_i] , recv_buffer[i]);
                }
            }
        }

    }

    template<class TVectorType, class TIndexVectorType >
    void Assemble(
        const TVectorType& rVectorInput,
        const TIndexVectorType& EquationId
    )
    {
        KRATOS_DEBUG_ERROR_IF(rVectorInput.size() != EquationId.size());

        for(unsigned int i=0; i<EquationId.size(); ++i){
            IndexType global_i = EquationId[i];

            if(GetNumbering().IsLocal(global_i))
            {
                IndexType local_i = GetNumbering().LocalId(global_i);
                AtomicAdd(mLocalData(local_i) , rVectorInput[i]);
            }
            else
            {
                auto owner_rank = GetNumbering().OwnerRank(global_i);
                IndexType local_i = GetNumbering().RemoteLocalId(global_i, owner_rank);
                auto it = (mNonLocalData[owner_rank].find( local_i ));
                KRATOS_DEBUG_ERROR_IF(it == mNonLocalData[owner_rank].end()) << "global_i = "<< global_i << " not in mNonLocalData" << std::endl;
                TDataType& value = (*it).second;
                AtomicAdd(value , rVectorInput[i]);
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
    virtual std::string Info() const
    {
std::stringstream buffer;
    buffer << "DistributedSystemVector" ;
    return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "DistributedSystemVector";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {
        std::cout << mLocalData << std::endl;
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
    const DataCommunicator& mrComm;
    typename DistributedNumbering<IndexType>::UniquePointer mpNumbering;

    DenseVector<TDataType> mLocalData; //contains the local data
    std::vector< std::unordered_map<IndexType, TDataType> > mNonLocalData;
    std::vector<int> mfem_assemble_colors; //coloring of communication

    std::unordered_map<IndexType, std::vector<IndexType> > mRecvIndicesByColor;

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

}; // Class DistributedSystemVector

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TDataType, class TIndexType>
inline std::istream& operator >> (std::istream& rIStream,
                DistributedSystemVector<TDataType,TIndexType>& rThis)
                {
                    return rIStream;
                }

/// output stream function
template<class TDataType, class TIndexType>
inline std::ostream& operator << (std::ostream& rOStream,
                const DistributedSystemVector<TDataType,TIndexType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_DISTRIBUTED_SYSTEM_VECTOR_H_INCLUDED  defined


