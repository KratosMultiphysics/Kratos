#include <iostream>
#include <iomanip>
#include <sstream>

template <typename T>
class FileIO {
public:

    // Creator & destructor
    FileIO(){};
    ~FileIO(){};

    /**
     * Writes the mesh in raw format and wipes previous content in the file.
     * @grid:     Value of the grid in Local or Global coordinatr system
     * @X:        X-Size of the grid
     * @Y:        Y-Size of the grid
     * @Z:        Z-Size of the grid
     * @fileName: Name of the output file
     **/
    void WriteGridWipe(T * grid, 
        const uint &X, const uint &Y, const uint &Z,
        const char * fileName) {

      std::ofstream outputFile(fileName);

      outputFile << std::fixed;
      outputFile << std::setprecision(2);

      for(uint k = 0; k < Z + BW; k++) {
        for(uint j = 0; j < Y + BW; j++) {
          for(uint i = 0; i < X + BW; i++) {
            outputFile << grid[k*(Y+BW)*(X+BW)+j*(X+BW)+i] << " ";
          } 
          outputFile << std::endl;
        }
        outputFile << std::endl;
      }
      outputFile << "%" << std::endl;
      
    }

    /**
     * Writes the mesh in raw format without wiping previous content in the file.
     * @grid:     Value of the grid in Local or Global coordinatr system
     * @X:        X-Size of the grid
     * @Y:        Y-Size of the grid
     * @Z:        Z-Size of the grid
     * @fileName: Name of the output file
     **/
    void WriteGrid(T * grid, 
        const uint &X, const uint &Y, const uint &Z,
        const char * fileName) {

      std::ofstream outputFile(fileName,std::ofstream::app);

      outputFile << std::fixed;
      outputFile << std::setprecision(2);

      for(uint k = 0; k < Z + BW; k++) {
        for(uint j = 0; j < Y + BW; j++) {
          for(uint i = 0; i < X + BW; i++) {
            outputFile << grid[k*(Y+BW)*(X+BW)+j*(X+BW)+i] << " ";
          } 
          outputFile << std::endl;
        }
        outputFile << std::endl;
      }
      outputFile << "%" << std::endl;
      
    }

    /**
     * Writes the mesh in GiD format.
     * @grid:     Value of the grid in Local or Global coordinatr system
     * @X:        X-Size of the grid
     * @Y:        Y-Size of the grid
     * @Z:        Z-Size of the grid
     * @fileName: Name of the output file
     **/
    void WriteGidMesh(T * grid, 
        const uint &X, const uint &Y, const uint &Z,
        const char * fileName) {

      std::ofstream outputFile(fileName);

      outputFile << "MESH \"Grid\" dimension 3 ElemType Hexahedra Nnode 8" << std::endl;
      outputFile << "# color 96 96 96" << std::endl;
      outputFile << "Coordinates" << std::endl;
      outputFile << "# node number coordinate_x coordinate_y coordinate_z  " << std::endl;

      for(uint k = 0; k < Z + BW; k++) {
        for(uint j = 0; j < Y + BW; j++) {
          uint cell = k*(Z+BW)*(Y+BW)+j*(Y+BW)+BWP;
          for(uint i = 0; i < X + BW; i++) {
            outputFile << cell++ << "  " << i << "  " << j << "  " << k << std::endl;
          }
        }
      }

      outputFile << "end coordinates" << std::endl;
      outputFile << "Elements" << std::endl;
      outputFile << "# Element node_1 node_2 node_3 node_4 node_5 node_6 node_7 node_8" << std::endl;

      for(uint k = BWP; k < Z + BWP-1; k++) {
        for(uint j = BWP; j < Y + BWP-1; j++) {
          uint cell = k*(Z+BW)*(Y+BW)+j*(Y+BW)+BWP;
          for(uint i = BWP; i < X + BWP; i++) {
            outputFile << cell++ << " ";

            outputFile << cell                        << " " << cell+1                    << "  ";
            outputFile << cell+1+(Y+BW)               << " " << cell+(Y+BW)               << "  ";
            outputFile << cell+(Z+BW)*(Y+BW)          << " " << cell+1+(Z+BW)*(Y+BW)      << "  ";
            outputFile << cell+1+(Z+BW)*(Y+BW)+(Y+BW) << " " << cell+(Z+BW)*(Y+BW)+(Y+BW) << "  ";

            outputFile << std::endl;
          }
        }
      }

      outputFile << "end Elements" << std::endl;
    }

    /**
     * Writes the results in GiD format.
     * @grid:     Value of the grid in Local or Global coordinatr system
     * @X:        X-Size of the grid
     * @Y:        Y-Size of the grid
     * @Z:        Z-Size of the grid
     * @results:  Output file stream
     * @step:     Step of the result
     **/
    void WriteGidResults(T * grid, 
        const uint &X, const uint &Y, const uint &Z,
        std::ofstream& results, int step) {

      results << "Result \"Temperature\" \"Kratos\" " << step << " Scalar OnNodes" << std::endl;
      results << "Values" << std::endl;

      for(uint k = 0; k < Z + BW; k++) {
        for(uint j = 0; j < Y + BW; j++) {
          uint cell = k*(Z+BW)*(Y+BW)+j*(Y+BW)+BWP;
          for(uint i = 0; i < X + BW; i++) {
            results << cell << "  " << grid[cell-1] << std::endl; cell++;
          }
        }
      }

      results << "End Values" << std::endl;
    }
};