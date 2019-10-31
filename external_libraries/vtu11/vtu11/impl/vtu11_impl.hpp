#ifndef VTU11_VTU11_IMPL_HPP
#define VTU11_VTU11_IMPL_HPP

#include "inc/xml.hpp"
#include "inc/utilities.hpp"

#include <limits>

namespace vtu11
{

/*
 * todo
 * - check consistency of each point & cell data set
 */

namespace detail
{

template<typename Writer, typename DataType>
inline void addDataSet( Writer& writer,
                        std::ostream& output,
                        const std::vector<DataType>& data,
                        size_t numberOfComponents = 1,
                        const std::string& name = "" )
{
  StringStringMap attributes = { { "type", dataTypeString<DataType>( ) } };

  if( numberOfComponents > 1 )
  {
    attributes["NumberOfComponents"] = std::to_string( numberOfComponents );
  }

  if( name != "" )
  {
    attributes["Name"] = name;
  }

  writer.addDataAttributes( attributes );

  if( attributes["format"] != "appended" )
  {
    ScopedXmlTag dataArrayTag( output, "DataArray", attributes );

    writer.writeData( output, data );
  }
  else
  {
    writeEmptyTag( output, "DataArray", attributes );

    writer.writeData( output, data );
  }
}

} // namespace detail


template<typename MeshGenerator, typename Writer>
void write( const std::string& filename,
            MeshGenerator& mesh,
            const std::vector<DataSet>& pointData,
            const std::vector<DataSet>& cellData,
            Writer writer )
{
  std::ofstream output;

  output.open( filename );

  VTU11_CHECK( output.is_open( ), "Failed to open file \"" + filename + "\"" );

  output << "<?xml version=\"1.0\"?>\n";

  StringStringMap headerAttributes { { "byte_order",  endianness( )      },
                                     { "type"      ,  "UnstructuredGrid" },
                                     { "version"   ,  "0.1"              } };

  writer.addHeaderAttributes( headerAttributes );

  {
    ScopedXmlTag vtkFileTag( output, "VTKFile", headerAttributes );
    {
      ScopedXmlTag unstructuredGridFileTag( output, "UnstructuredGrid", { } );
      {
        ScopedXmlTag pieceTag( output, "Piece", { { "NumberOfPoints", std::to_string( mesh.numberOfPoints( ) ) },
                                                  { "NumberOfCells" , std::to_string( mesh.numberOfCells( ) )    } } );

        {
          ScopedXmlTag pointDataTag( output, "PointData", { } );

          for( const auto& dataSet : pointData )
          {
            detail::addDataSet( writer, output, std::get<2>( dataSet ), std::get<1>( dataSet ), std::get<0>( dataSet ) );
          }

        } // PointData

        {
          ScopedXmlTag cellDataTag( output, "CellData", { } );

          for( const auto& dataSet : cellData )
          {
            detail::addDataSet( writer, output, std::get<2>( dataSet ), std::get<1>( dataSet ), std::get<0>( dataSet ) );
          }

        } // CellData

        {
          ScopedXmlTag pointsTag( output, "Points", { } );

          detail::addDataSet( writer, output, mesh.points( ), 3 );

        } // Points

        {
          ScopedXmlTag pointsTag( output, "Cells", { } );

          detail::addDataSet( writer, output, mesh.connectivity( ), 1, "connectivity" );
          detail::addDataSet( writer, output, mesh.offsets( ), 1, "offsets" );
          detail::addDataSet( writer, output, mesh.types( ), 1, "types" );

        } // Cells

      } // Piece
    } // UnstructuredGrid

    auto appendedAttributes = writer.appendedAttributes( );

    if( !appendedAttributes.empty( ) )
    {
      ScopedXmlTag appendedDataTag( output, "AppendedData", appendedAttributes );

      output << "_";

      writer.writeAppended( output );

    } // AppendedData
  } // VTKFile

  output.close( );
}

} // namespace vtu11

#endif // VTU11_VTU11_IMPL_HPP
