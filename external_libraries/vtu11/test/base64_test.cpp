#include "external/catch2/catch.hpp"
#include "inc/base64.hpp"

#include <vector>
#include <string>

namespace vtu11
{

TEST_CASE("base64encode_test")
{
  std::string test1 = "hell";
  std::string test2 = "hello";
  std::string test3 = "hello1";
  std::string test4 = "Lorem ipsum dolor sit amet, consectetur adipiscing elit, "
                      "sed do eiusmod tempor incididunt ut labore et dolore magna aliqua";

  CHECK( base64Encode( test1.begin( ), test1.end( ) ) == "aGVsbA==" );
  CHECK( base64Encode( test2.begin( ), test2.end( ) ) == "aGVsbG8=" );
  CHECK( base64Encode( test3.begin( ), test3.end( ) ) == "aGVsbG8x" );

  std::string expected4 = "TG9yZW0gaXBzdW0gZG9sb3Igc2l0IGFtZXQsIGNvbnNlY3RldHVyIGFkaXBpc2NpbmcgZWxpdCwgc2VkIG"
                          "RvIGVpdXNtb2QgdGVtcG9yIGluY2lkaWR1bnQgdXQgbGFib3JlIGV0IGRvbG9yZSBtYWduYSBhbGlxdWE=";

  CHECK( base64Encode( test4.begin( ), test4.end( ) ) == expected4 );

  std::string empty;

  CHECK( base64Encode( empty.begin( ), empty.end( ) ) == "" );


  std::vector<double> test5 { -0.01, 0.05, 4.0 };
  std::vector<double> test6 { 32.23, 23.63, 33.52, 90.21 };
  std::vector<double> test7 { 23.046, 2.3, -457.129, 84762.3, -0.423 };
  std::vector<double> test8 { 1.0,-2.0,-3.0,4.0,5.0, 6.0 };

  CHECK( base64Encode( test5.begin( ), test5.end( ) ) == "exSuR+F6hL+amZmZmZmpPwAAAAAAABBA" );
  CHECK( base64Encode( test6.begin( ), test6.end( ) ) == "PQrXo3AdQEDhehSuR6E3QMP1KFyPwkBAPQrXo3CNVkA=" );
  CHECK( base64Encode( test7.begin( ), test7.end( ) ) == "sp3vp8YLN0BmZmZmZmYCQPLSTWIQknzAzczMzKSx9EDfT42XbhLbvw==" );
  CHECK( base64Encode( test8.begin( ), test8.end( ) ) == "AAAAAAAA8D8AAAAAAAAAwAAAAAAAAAjAAAAAAAAAEEAAAAAAAAAUQAAAAAAAABhA" );
}

} // namespace vtu11

