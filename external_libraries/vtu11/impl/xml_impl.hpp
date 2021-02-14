//          __        ____ ____
// ___  ___/  |_ __ _/_   /_   |
// \  \/ /\   __\  |  \   ||   |
//  \   /  |  | |  |  /   ||   |
//   \_/   |__| |____/|___||___|
//
//  License: BSD License ; see LICENSE
//

#ifndef VTU11_XML_IMPL_HPP
#define VTU11_XML_IMPL_HPP

namespace vtu11
{
namespace detail
{

inline void writeTag( std::ostream& output,
                      const std::string& name,
                      const StringStringMap& attributes,
                      const std::string& tagEnd )
{
  output << "<" << name;

  for( const auto& attribute : attributes )
  {
    output << " " << attribute.first << "=\"" << attribute.second << "\"";
  }

  output << tagEnd << "\n";
}

} // namespace detail

inline ScopedXmlTag::ScopedXmlTag( std::ostream& output,
                                   const std::string& name,
                                   const StringStringMap& attributes ) :
   closeTag( [ &output, name ]( ){ output << "</" << name << ">\n"; } )
{
  detail::writeTag( output, name, attributes, ">" );
}


inline ScopedXmlTag::~ScopedXmlTag( )
{
  closeTag( );
}


inline void writeEmptyTag( std::ostream& output,
                           const std::string& name,
                           const StringStringMap& attributes )
{
  detail::writeTag( output, name, attributes, "/>" );
}

} // namespace vtu11

#endif // VTU11_XML_HPP
