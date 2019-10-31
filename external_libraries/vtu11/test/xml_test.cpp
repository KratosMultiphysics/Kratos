#include "external/catch2/catch.hpp"
#include "inc/xml.hpp"

#include <sstream>

namespace vtu11
{

TEST_CASE("ScopedXmlTag_test")
{
  std::ostringstream output;

  {
    ScopedXmlTag tag1( output, "Test1", { { "attr1", "7" }, { "attr2", "nice" } } );
    {
      ScopedXmlTag tag2( output, "Test2", { { "attr3", "43.32" }, { "attr4", "[2, 3]" } } );

      output << "dadatata" << "\n";
    }
  }

  std::string expectedString =
    "<Test1 attr1=\"7\" attr2=\"nice\">\n"
        "<Test2 attr3=\"43.32\" attr4=\"[2, 3]\">\n"
            "dadatata\n"
        "</Test2>\n"
    "</Test1>\n";

  CHECK( output.str( ) == expectedString );
}

} // namespace vtu11

