#if !defined(KRATOS_BFECC_UTILS)
#define KRATOS_BFECC_UTILS

#include <cstddef>
#include <vector>

enum BFECCFLAGS {
  FIXED_VELOCITY_X = 0x01,
  FIXED_VELOCITY_Y = 0x02,
  FIXED_VELOCITY_Z = 0x04,
  FIXED_PRESSURE   = 0x08,
  OUT_OF_BOUNDS    = 0x10
};

class BfeccUtils {
public:
  static void GlobalToLocal(double * coord, double f, const std::size_t &dim) {
    for(std::size_t d = 0; d < dim; d++)
      coord[d] *= f;
  }

  static void LocalToGlobal(double * coord, double f, const std::size_t &dim) {
    for(std::size_t d = 0; d < dim; d++)
      coord[d] /= f;
  }

  static std::size_t Index(
      const std::size_t & i, const std::size_t & j, const std::size_t & k,
      const std::vector<std::size_t> numCells, const std::vector<std::size_t> borderWidth ) {
    return k * (numCells[2] + borderWidth[2] * 2) * (numCells[1] + borderWidth[1] * 2) + j * (numCells[1] + borderWidth[1] * 2) + i;
  }

};

#endif // KRATOS_BFECC_UTILS defined
