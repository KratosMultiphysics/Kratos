//          __        ____ ____
// ___  ___/  |_ __ _/_   /_   |
// \  \/ /\   __\  |  \   ||   |
//  \   /  |  | |  |  /   ||   |
//   \_/   |__| |____/|___||___|
//
//  License: BSD License ; see LICENSE
//

#ifndef VTU11_VTU11_HPP
#define VTU11_VTU11_HPP

#include "inc/alias.hpp"
#include "inc/writer.hpp"

namespace vtu11
{

struct Vtu11UnstructuredMesh
{
  std::vector<double>& points_;
  std::vector<size_t>& connectivity_;
  std::vector<size_t>& offsets_;
  std::vector<VtkCellType>& types_;

  std::vector<double>& points( ){ return points_; }
  std::vector<size_t>& connectivity( ){ return connectivity_; }
  std::vector<size_t>& offsets( ){ return offsets_; }
  std::vector<VtkCellType>& types( ){ return types_; }

  size_t numberOfPoints( ){ return points_.size( ) / 3; }
  size_t numberOfCells( ){ return types_.size( ); }
};

template<typename MeshGenerator, typename Writer = AsciiWriter>
void write( const std::string& filename,
            MeshGenerator& mesh,
            const std::vector<DataSet>& pointData,
            const std::vector<DataSet>& cellData,
            Writer writer = Writer( ) );

} // namespace vtu11

#include "impl/vtu11_impl.hpp"

#endif // VTU11_VTU11_HPP
