#include "hdf5_io.h"

#include <sstream>
#include "hdf5.h"

namespace Kratos
{

HDF5IO::HDF5IO(std::string FileName, Flags Options):
  IO(),
  mFileName(FileName)
{

}

HDF5IO::~HDF5IO()
{

}

void HDF5IO::WriteModelPart(ModelPart& rModelPart)
{
    hid_t FileId = H5Fcreate(mFileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    herr_t Status;

    // write nodes
    double *pNodeData = new double[rModelPart.NumberOfNodes()*3];
    int i = 0;
    for (ModelPart::NodeIterator iNode = rModelPart.NodesBegin(); iNode != rModelPart.NodesEnd(); ++iNode)
    {
        auto& Coordinates = iNode->Coordinates();
        pNodeData[i++] = Coordinates[0];
        pNodeData[i++] = Coordinates[1];
        pNodeData[i++] = Coordinates[2];
    }

    hsize_t NodeDimensions[2];
    NodeDimensions[0] = rModelPart.NumberOfNodes();
    NodeDimensions[1] = 3;

    hid_t NodeSpace = H5Screate_simple(2, NodeDimensions, NULL);
    hid_t NodeDataSet = H5Dcreate(FileId, "/Nodes", H5T_NATIVE_DOUBLE, NodeSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    Status = H5Dwrite(NodeDataSet, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pNodeData);
    Status = H5Dclose(NodeDataSet);
    Status = H5Sclose(NodeSpace);

    delete[] pNodeData;

    // Write elements
    int NumElements = rModelPart.NumberOfElements();
    int NodesInElement = rModelPart.ElementsBegin()->GetGeometry().PointsNumber(); //TODO: Do not assume constant sizes!
    int *pConnectivities = new int[NumElements*NodesInElement];
    i = 0;
    for (ModelPart::ElementIterator iElem = rModelPart.ElementsBegin(); iElem != rModelPart.ElementsEnd(); ++iElem)
    {
        Geometry< Node<3> > &rGeometry = iElem->GetGeometry();
        for (unsigned int j = 0; j < rGeometry.PointsNumber(); j++)
        {
            pConnectivities[i++] = rGeometry[j].Id() - 1; //TODO: This needs to go through a permutation in MPI!!!
        }
    }

    hsize_t ElementDimensions[2];
    ElementDimensions[0] = NumElements;
    ElementDimensions[1] = NodesInElement;

    hid_t ElementSpace = H5Screate_simple(2, ElementDimensions, NULL);
    hid_t ElementDataSet = H5Dcreate(FileId, "/Elements", H5T_NATIVE_INT, ElementSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    Status = H5Dwrite(ElementDataSet, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, pConnectivities);
    Status = H5Dclose(ElementDataSet);
    Status = H5Sclose(ElementSpace);

    delete[] pConnectivities;

    // Write results
    int ResultSize = rModelPart.NumberOfNodes();
    double *pResults = new double[ResultSize];
    i = 0;
    for (ModelPart::NodeIterator iNode = rModelPart.NodesBegin(); iNode != rModelPart.NodesEnd(); ++iNode)
    {
        auto& Value = iNode->FastGetSolutionStepValue(VELOCITY_X);
        pResults[i++] = Value;
    }

    hsize_t DataDimensions[1];
    DataDimensions[0] = rModelPart.NumberOfNodes();

    hid_t ResultDataSpace = H5Screate_simple(1, DataDimensions, NULL);
    hid_t ResultDataSet = H5Dcreate(FileId, "/VelocityX", H5T_NATIVE_DOUBLE, ResultDataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    Status = H5Dwrite(ResultDataSet, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pResults);
    Status = H5Dclose(ResultDataSet);
    Status = H5Sclose(ResultDataSpace);
    delete[] pResults;

    // Close file
    Status = H5Fclose(FileId);
}

/// Turn back information as a string.
std::string HDF5IO::Info() const
{
  std::stringstream msg;
  msg << "HDF5IO" << std::endl;
  return msg.str();
}

/// Print information about this object.
void HDF5IO::PrintInfo(std::ostream& rOStream) const
{
}

/// Print object's data.
void HDF5IO::PrintData(std::ostream& rOStream) const
{

}

}
