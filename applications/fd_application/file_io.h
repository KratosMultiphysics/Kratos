#include <iostream>
#include <iomanip>
#include <sstream>

template <typename T>
class FileIO {
private:

    std::stringstream name_mesh;
    std::stringstream name_post;

    std::ofstream * mesh_file;
    std::ofstream * post_file;

    double mdx;

    /**
     * Sets the name for the mesh file
     * @name:     name of the mesh file
     **/
    void SetMeshName(const char * mesh, const uint &N) {

      name_mesh << mesh << N << ".post.msh";
    }

    /**
     * Sets the name for the post file
     * @name:     name of the post file
     **/
    void setPostName(const char * post, const uint &N) {

      name_post << post << N << ".post.res";
    }

    /**
     * Initial formating for the post file
     **/
    void PreparePostFile() {
      (*post_file) << "GiD Post Results File 1.0" << std::endl << std::endl;
    }

public:

    // Creator & destructor
    FileIO(const char * name, const uint &N, const double &dx) {

      SetMeshName(name,N);
      setPostName(name,N);

      mdx = dx;

      mesh_file = new std::ofstream(name_mesh.str().c_str());
      post_file = new std::ofstream(name_post.str().c_str());

      PreparePostFile();
    };

    ~FileIO() {

      mesh_file->close();
      post_file->close();

      delete mesh_file;
      delete post_file;
    };

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
     **/
    void WriteGidMesh(T * grid, 
        const uint &X, const uint &Y, const uint &Z) {

      (*mesh_file) << "MESH \"Grid\" dimension 3 ElemType Hexahedra Nnode 8" << std::endl;
      (*mesh_file) << "# color 96 96 96" << std::endl;
      (*mesh_file) << "Coordinates" << std::endl;
      (*mesh_file) << "# node number coordinate_x coordinate_y coordinate_z  " << std::endl;

      for(uint k = 0; k < Z + BW; k++) {
        for(uint j = 0; j < Y + BW; j++) {
          uint cell = k*(Z+BW)*(Y+BW)+j*(Y+BW)+BWP;
          for(uint i = 0; i < X + BW; i++) {
            (*mesh_file) << cell++ << "  " << i * mdx << "  " << j * mdx << "  " << k * mdx << std::endl;
          }
        }
      }

      (*mesh_file) << "end coordinates" << std::endl;
      (*mesh_file) << "Elements" << std::endl;
      (*mesh_file) << "# Element node_1 node_2 node_3 node_4 node_5 node_6 node_7 node_8" << std::endl;

      for(uint k = BWP; k < Z + BWP-1; k++) {
        for(uint j = BWP; j < Y + BWP-1; j++) {
          uint cell = k*(Z+BW)*(Y+BW)+j*(Y+BW)+BWP;
          for(uint i = BWP; i < X + BWP; i++) {
            (*mesh_file) << cell++ << " ";

            (*mesh_file) << cell                        << " " << cell+1                    << "  ";
            (*mesh_file) << cell+1+(Y+BW)               << " " << cell+(Y+BW)               << "  ";
            (*mesh_file) << cell+(Z+BW)*(Y+BW)          << " " << cell+1+(Z+BW)*(Y+BW)      << "  ";
            (*mesh_file) << cell+1+(Z+BW)*(Y+BW)+(Y+BW) << " " << cell+(Z+BW)*(Y+BW)+(Y+BW) << "  ";

            (*mesh_file) << std::endl;
          }
        }
      }

      (*mesh_file) << "end Elements" << std::endl;
    }

    /**
     * Writes the results in GiD format.
     * @grid:     Value of the grid in Local or Global coordinatr system
     * @X:        X-Size of the grid
     * @Y:        Y-Size of the grid
     * @Z:        Z-Size of the grid
     * @step:     Step of the result
     **/
    void WriteGidResults(T * grid, 
        const uint &X, const uint &Y, const uint &Z, int step) {

      (*post_file) << "Result \"Temperature\" \"Kratos\" " << step << " Scalar OnNodes" << std::endl;
      (*post_file) << "Values" << std::endl;

      for(uint k = 0; k < Z + BW; k++) {
        for(uint j = 0; j < Y + BW; j++) {
          uint cell = k*(Z+BW)*(Y+BW)+j*(Y+BW)+BWP;
          for(uint i = 0; i < X + BW; i++) {
            (*post_file) << cell << "  " << grid[cell-1] << std::endl; cell++;
          }
        }
      }

      (*post_file) << "End Values" << std::endl;
    }
};