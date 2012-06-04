/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
mossaiby@yahoo.com
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

#include "includes/serializer.h"

#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"

class ParticleSpatialConfigure
{
public:

    static const std::size_t Dim = 3;
    static const std::size_t Dimension = 3;

    typedef std::size_t  SizeType;
    typedef std::size_t  IndexType;

    typedef Kratos::Point<Dim>                  PointType;
    typedef Kratos::Tetrahedra3D4<PointType>    ObjectType;

    typedef PointType*                          PtrPointType;
    typedef ObjectType*                         PtrObjectType;

    typedef PtrPointType*                       PointVector;
    typedef PtrPointType*                       PointIterator;

    typedef PtrObjectType*                      ObjectVector;
    typedef PtrObjectType*                      ObjectIterator;

    typedef double*                             DistanceVector;
    typedef double*                             DistanceIterator;

    typedef Kratos::SearchUtils::SquaredDistanceFunction<Dim,PointType> DistanceFunction;
    
    static inline void MPISave(std::vector<std::vector<PtrPointType> >& inputElements, std::string messages[]) 
    {
        int mpi_rank;
        int mpi_size;
      
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
      
        std::stringstream * serializer_buffer;
        
        char * recvBuffers[mpi_size];
        
        for(int i = 0; i < mpi_size; i++)
        {
            if(mpi_rank != i)
            {
                Kratos::Serializer particleSerializer;
                particleSerializer.save("nodes",inputElements[i]);
                  
                serializer_buffer = (std::stringstream *)particleSerializer.pGetBuffer();
                messages[i] = std::string(serializer_buffer->str());
            }
        }
    }

    static inline void MPILoad(std::vector<std::vector<PtrPointType> >& outputElements, std::string messages[]) 
    {
        int mpi_rank;
        int mpi_size;
      
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
        
        std::stringstream * serializer_buffer;
        
        for (int i = 0; i < mpi_size; i++)
        { 
            if (i != mpi_rank && messages[i].size())
            {
                Kratos::Serializer particleSerializer;
                serializer_buffer = (std::stringstream *)particleSerializer.pGetBuffer();
                serializer_buffer->write((char*)(messages[i].c_str()), messages[i].size());

                particleSerializer.load("nodes",outputElements[i]);
            }
        }
    }
};

