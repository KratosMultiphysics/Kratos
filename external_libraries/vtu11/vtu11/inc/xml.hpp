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
