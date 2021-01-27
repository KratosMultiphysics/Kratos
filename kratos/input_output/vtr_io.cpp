//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   
//

// System includes
#include <fstream>

// External includes


// Project includes
#include "includes/define.h"
#include "input_output/vtr_io.h"


namespace Kratos
{

namespace{

    class VtkRectilinearGridPieceReader{
        std::vector<Internals::CartesianMeshColors>&  mMultiBlockMesh;
    public:
        VtkRectilinearGridPieceReader(std::vector<Internals::CartesianMeshColors>& MultiBlockMesh) : mMultiBlockMesh(MultiBlockMesh) {}

        void operator()(XmlIO& TheIO){
            auto& current_block_info = TheIO.GetCurrentBlockInfo();
            KRATOS_ERROR_IF_NOT(current_block_info.HasAttribute("Extent")) << " The RectilinearGrid Piece should provide the Extent attribute. i.e: WholeExtent=\"0 10 0 10 0 10\"" << std::endl;
            std::stringstream extent(current_block_info.GetAttribute("Extent"));

            int min_i, max_i, min_j, max_j, min_k, max_k;
            extent >> min_i >> max_i >> min_j >> max_j >> min_k >> max_k;

            KRATOS_ERROR_IF(min_i >= max_i) << "min extent=" << min_i << " should be smaller than max extent=" << max_i << std::endl;
            KRATOS_ERROR_IF(min_j >= max_j) << "min extent=" << min_j << " should be smaller than max extent=" << max_j << std::endl;
            KRATOS_ERROR_IF(min_k >= max_k) << "min extent=" << min_k << " should be smaller than max extent=" << max_k << std::endl;

            KRATOS_WATCH(min_i);           
            KRATOS_WATCH(max_i);           
            KRATOS_WATCH(min_j);           
            KRATOS_WATCH(max_j);           
            KRATOS_WATCH(min_k);           
            KRATOS_WATCH(max_k);     
    
            std::vector<double> x_coordinates(max_i - min_i + 1, -1.00);
            std::vector<double> y_coordinates(max_j - min_j + 1, -1.00);
            std::vector<double> z_coordinates(max_k - min_k + 1, -1.00);

            TheIO.SetBlockAction("Coordinates", [&](XmlIO& TheIO){
                TheIO.ReadBlock("DataArray", [&x_coordinates](XmlIO& TheIO){ // Here I am forcing to have DataArray block in the input with the given action.
                    TheIO.ReadBlockContent(x_coordinates);
                });

                TheIO.ReadBlock("DataArray", [&y_coordinates](XmlIO& TheIO){ // Here I am forcing to have DataArray block in the input with the given action.
                    TheIO.ReadBlockContent(y_coordinates);
                });

                TheIO.ReadBlock("DataArray", [&z_coordinates](XmlIO& TheIO){ // Here I am forcing to have DataArray block in the input with the given action.
                    TheIO.ReadBlockContent(z_coordinates);
                });

                this->mMultiBlockMesh.push_back(Internals::CartesianMeshColors());
                this->mMultiBlockMesh.back().SetCoordinates(x_coordinates, y_coordinates, z_coordinates);
            });

            TheIO.SetBlockAction("CellData", [&](XmlIO& TheIO){
                TheIO.SetBlockAction("DataArray", [&](XmlIO& TheIO){ 

                    auto& current_block_info = TheIO.GetCurrentBlockInfo();
                    KRATOS_ERROR_IF_NOT(current_block_info.HasAttribute("Name")) << " The cell data array should provide the Name attribute" << std::endl;

                    std::string const& data_name = current_block_info.GetAttribute("Name");
                    auto& mesh = this->mMultiBlockMesh.back();
                    std::size_t size = mesh.GetElementalColors().size();
                    std::vector<double> data(size);
                    TheIO.ReadBlockContent(data);

                    mesh.SetElementalData(data_name, data);

                });
                TheIO.Read();
            });

            TheIO.Read();

        }

    };

    class VtkRectilinearGridReader{
        std::vector<Internals::CartesianMeshColors>& mMultiBlockMesh;
        
        public:
        VtkRectilinearGridReader(std::vector<Internals::CartesianMeshColors>& MultiBlockMesh) : mMultiBlockMesh(MultiBlockMesh) {}

        void operator()(XmlIO& TheIO){
            auto& current_block_info = TheIO.GetCurrentBlockInfo();
            KRATOS_ERROR_IF_NOT(current_block_info.HasAttribute("WholeExtent")) << " The RectilinearGrid should provide the WholeExtent attribute. i.e: WholeExtent=\"0 10 0 10 0 10\"" << std::endl;
            std::stringstream extent(current_block_info.GetAttribute("WholeExtent"));
            int min_i, max_i, min_j, max_j, min_k, max_k;
            extent >> min_i >> max_i >> min_j >> max_j >> min_k >> max_k;
            // KRATOS_WATCH(min_i);           
            // KRATOS_WATCH(max_i);           
            // KRATOS_WATCH(min_j);           
            // KRATOS_WATCH(max_j);           
            // KRATOS_WATCH(min_k);           
            // KRATOS_WATCH(max_k);     


 
            TheIO.SetBlockAction("Piece", VtkRectilinearGridPieceReader(mMultiBlockMesh));

            TheIO.Read();
       }



    };

}

    /// Constructor with  filenames.
    VtrIO::VtrIO(std::string const& Filename)
        : mXmlIO(Filename)  {
        }

    VtrIO::VtrIO(std::iostream* pInputStream) : mXmlIO(pInputStream){}

    void VtrIO::Read(std::vector<Internals::CartesianMeshColors>& MultiBlockMesh){

        mXmlIO.SetBlockAction("VTKFile", [&MultiBlockMesh](XmlIO& TheIO){
            auto& current_block_info = TheIO.GetCurrentBlockInfo();
            KRATOS_ERROR_IF(current_block_info.Name() != "VTKFile") << " The file is not a valid VTK xml file" << std::endl;
            if(current_block_info.HasAttribute("type")){
                KRATOS_ERROR_IF(current_block_info.GetAttribute("type") != "RectilinearGrid") << " The type for vtr file should be RectilinearGrid" << std::endl;
            }
            if(current_block_info.HasAttribute("byte_order")){
                KRATOS_ERROR_IF(current_block_info.GetAttribute("byte_order") != "LittleEndian") << " Only LittleEndian byte order is supported" << std::endl;
            }

            TheIO.SetBlockAction("RectilinearGrid", VtkRectilinearGridReader(MultiBlockMesh));

            TheIO.Read();
        });

        mXmlIO.Read();
    }
 
    /// Turn back information as a string.
    std::string VtrIO::Info() const
    {
        return "VTR IO";
    }

    /// Print information about this object.
    void VtrIO::PrintInfo(std::ostream& rOStream) const{
        rOStream << Info();
    }

    /// Print object's data.
    void VtrIO::PrintData(std::ostream& rOStream) const{

    }



}  // namespace Kratos.


