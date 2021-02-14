//          __        ____ ____
// ___  ___/  |_ __ _/_   /_   |
// \  \/ /\   __\  |  \   ||   |
//  \   /  |  | |  |  /   ||   |
//   \_/   |__| |____/|___||___|
//
//  License: BSD License ; see LICENSE
//

#ifndef VTU11_XML_HPP
#define VTU11_XML_HPP

#include <functional>
#include "inc/alias.hpp"

namespace vtu11
{

class ScopedXmlTag final
{
public:
  ScopedXmlTag( std::ostream& output,
                const std::string& name,
                const StringStringMap& attributes );

  ~ScopedXmlTag( );

private:
  std::function<void()> closeTag;
};

void writeEmptyTag( std::ostream& output,
                    const std::string& name,
                    const StringStringMap& attributes );

} // namespace vtu11

#include "impl/xml_impl.hpp"

#endif // VTU11_XML_HPP
