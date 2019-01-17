/**
 *
 * @file pastix/nompi.h
 *
 * PaStiX header to redefine all MPI keywords in order to allow compilation
 * without MPI.
 *
 * @copyright 2004-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2018-07-16
 *
 **/
#ifndef _pastix_nompi_h_
#define _pastix_nompi_h_

#define MPI_Datatype int
#define MPI_Op       int
#define MPI_Request  int
#define MPI_Aint     int
#define MPI_Comm     int
#define MPI_Fint     int

#define MPI_CHAR              1
#define MPI_UNSIGNED_CHAR     2
#define MPI_BYTE              3
#define MPI_SHORT             4
#define MPI_UNSIGNED_SHORT    5
#define MPI_INT               6
#define MPI_UNSIGNED          7
#define MPI_LONG              8
#define MPI_UNSIGNED_LONG     9
#define MPI_FLOAT             10
#define MPI_DOUBLE            11
#define MPI_LONG_DOUBLE       12
#define MPI_LONG_LONG_INT     13
#define MPI_LONG_LONG         13
#define MPI_PACKED            14
#define MPI_LB                15
#define MPI_UB                16
#define MPI_FLOAT_INT         17
#define MPI_DOUBLE_INT        18
#define MPI_LONG_INT          19
#define MPI_SHORT_INT         20
#define MPI_2INT              21
#define MPI_LONG_DOUBLE_INT   22
#define MPI_COMPLEX           23
#define MPI_DOUBLE_COMPLEX    24
#define MPI_LOGICAL           25
#define MPI_REAL              26
#define MPI_DOUBLE_PRECISION  27
#define MPI_INTEGER           28
#define MPI_2INTEGER          29
#define MPI_2COMPLEX          30
#define MPI_2DOUBLE_COMPLEX   31
#define MPI_2REAL             32
#define MPI_2DOUBLE_PRECISION 33
#define MPI_INTEGER4          34
#define MPI_INTEGER8          35
#define MPI_CHARACTER         1

#define MPI_SUCCESS           0
#define MPI_ERR_COUNT         -1
#define MPI_ERR_TYPE          -2

#define MPI_COMM_WORLD        0
#define MPI_COMM_SELF         1
#define MPI_ANY_SOURCE        2
#define MPI_ANY_TAG           3
#define MPI_REQUEST_NULL      4
#define MPI_UNDEFINED         5
#define MPI_SUM               8
#define MPI_MAX               9
#define MPI_BOTTOM            NULL
#define MPI_THREAD_MULTIPLE   11
#define MPI_THREAD_SINGLE     12
#define MPI_THREAD_FUNNELED   13
#define MPI_THREAD_SERIALIZED 14

#define MPI_MAX_PROCESSOR_NAME 1

typedef void (MPI_User_function) ( void *, void *, int *, MPI_Datatype * );

typedef struct MPI_Status{
    int MPI_SOURCE;
    int MPI_TAG;
    int MPI_ERROR;
} MPI_Status;

#define MPI_Get_processor_name(__name, __len) { (void)(__name); *(__len) = 0; }
#define MPI_Send(buf, count, datatype, dest, tag, comm)
#define MPI_Recv(buf, count, datatype, source, tag, comm, status)
#define MPI_Irecv(buf, count, datatype, source, tag, comm, request);
#define MPI_Isend(buf, count, datatype, dest, tag, comm, request)
#define MPI_Wait(request, status)
#define MPI_Waitany(count, array_of_requests, index, status);
#define MPI_Cancel(request)
#define MPI_Test(request, flag, status)
#define MPI_Testany(count, array_of_requests, index, flag, array_of_statuses)
#define MPI_Iprobe(source, tag, comm, flag, status)
#define MPI_Recv_init(buf, count, datatype, dest, tag, comm, request)
#define MPI_Start(request)
#define MPI_Startall(count, array_of_requests)
#define MPI_Type_contiguous(count, oldtype, newtype)
#define MPI_Type_struct(count, array_of_blocklengths, array_of_displacement, \
                        oldtype, newtype)
#define MPI_Address(location, newtype)
#define MPI_Type_commit(datatype)
#define MPI_Type_free(datatype)
#define MPI_Request_free(request)
#define MPI_Barrier(comm)
#define MPI_Op_create(function, commute, op)
#define MPI_Init(argc, argv)
#define MPI_Initialized(_init_) do { *(_init_) = 1; } while(0);
#define MPI_Finalize()
#define MPI_Comm_split(comm, color, id, new_comm)

#define MPI_Allgather(sendbuf, sendcount, sendtype, recvbuf,            \
                      recvcount, recvtype, comm)    recvbuf = sendbuf;

#define MPI_Get_count(status, datatype ,count) *(count) = 0;
#define MPI_Comm_size(comm, size) *(size)=1;
#define MPI_Comm_rank(comm, rank) *(rank)=0;


static inline int
pastix_nompi_copy( void *dst, const void *src, int count, MPI_Datatype datatype )
{
    if ( count < 0 ) {
        return MPI_ERR_COUNT;
    }

    if(  src != dst ) {
        switch (datatype) {
        case MPI_INT:
            memcpy(dst, src, count * sizeof(int));
            break;
        case MPI_LONG:
            memcpy(dst, src, count * sizeof(long));
            break;
        case MPI_INTEGER4:
            memcpy(dst, src, count * sizeof(int32_t));
            break;
        case MPI_INTEGER8:
            memcpy(dst, src, count * sizeof(int64_t));
            break;
        case MPI_FLOAT:
            memcpy(dst, src, count * sizeof(float));
            break;
        case MPI_DOUBLE:
            memcpy(dst, src, count * sizeof(double));
            break;
        case MPI_COMPLEX:
            memcpy(dst, src, count * 2 * sizeof(float));
            break;
        case MPI_DOUBLE_COMPLEX:
            memcpy(dst, src, count * 2 * sizeof(double));
            break;
        default:
            fprintf(stderr,"pastix_nompi_copy: Unknown MPI datatype\n");
            return MPI_ERR_TYPE;
        }
    }
    return MPI_SUCCESS;
}

static inline int
MPI_Gather( const void *sendbuf, int sendcount, MPI_Datatype sendtype,
            void *recvbuf, int recvcount, MPI_Datatype recvtype,
            int root, MPI_Comm comm )
{
    assert( sendcount == recvcount );
    assert( sendtype == recvtype );

    (void)recvcount;
    (void)recvtype;
    (void)root;
    (void)comm;
    return pastix_nompi_copy( recvbuf, sendbuf, sendcount, sendtype );
}

static inline int
MPI_Allreduce( const void *sendbuf, void *recvbuf, int count,
               MPI_Datatype datatype, MPI_Op op, MPI_Comm comm )
{
    (void)op;
    (void)comm;
    return pastix_nompi_copy( recvbuf, sendbuf, count, datatype );
}

static inline int
MPI_Alltoall( const void *sendbuf, int sendcount, MPI_Datatype sendtype,
              void *recvbuf, int recvcount, MPI_Datatype recvtype,
              MPI_Comm comm )
{
    assert( sendcount == recvcount );
    assert( sendtype == recvtype );

    (void)recvcount;
    (void)recvtype;
    (void)comm;
    return pastix_nompi_copy( recvbuf, sendbuf, sendcount, sendtype );
}

static inline int
MPI_Reduce( const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
            MPI_Op op, int root, MPI_Comm comm )
{
    (void)op;
    (void)root;
    (void)comm;
    return pastix_nompi_copy( recvbuf, sendbuf, count, datatype );
}

static inline int
MPI_Bcast( void *buffer, int count, MPI_Datatype datatype, int root,
           MPI_Comm comm )
{
    (void)buffer;
    (void)count;
    (void)datatype;
    (void)root;
    (void)comm;
    return MPI_SUCCESS;
}

#define MPI_Comm_f2c(comm) 0;
#define MPI_Init_thread(argc, argv,required,provided)

#endif /* _pastix_nompi_h_ */
