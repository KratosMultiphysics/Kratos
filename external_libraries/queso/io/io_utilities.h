// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef IO_UTILTIES_H
#define IO_UTILTIES_H

//// STL includes
#include <fstream>
//// Project includes
#include "containers/element_container.hpp"
#include "containers/triangle_mesh.hpp"
#include "containers/boundary_integration_point.hpp"

namespace queso {

///@name QuESo Classes
///@{

/**
 * @class  IO
 * @author Manuel Messmer
 * @brief  Provides methods to parse data. Supports STL and VTK files.
*/
class IO{

public:
  ///@name Operations
  ///@{

  /// @brief Write TriangleMeshInterface to VTK-File
  /// @param rTriangleMesh
  /// @param Filename
  /// @param Binary
  /// @return bool
  static bool WriteMeshToVTK(const TriangleMeshInterface& rTriangleMesh,
                             const char* Filename,
                             const bool Binary);

  /// @brief Write TriangleMeshInterface to STL-File.
  /// @param rTriangleMesh
  /// @param Filename
  /// @param Binary
  /// @return bool
  static bool WriteMeshToSTL(const TriangleMeshInterface& rTriangleMesh,
                             const char* Filename,
                             const bool Binary);

  /// @brief Read TriangleMeshInterface from STL.
  /// @param rTriangleMesh
  /// @param Filename
  /// @return bool
  static bool ReadMeshFromSTL(TriangleMeshInterface& rTriangleMesh,
                              const char* Filename);

  /// @brief Write displacements to VTK-file. Append exisiting files, that contains vertices.
  /// @param rDisplacement
  /// @param Filename
  /// @param Binary
  /// @return bool
  static bool WriteDisplacementToVTK(const std::vector<Vector3d>& rDisplacement,
                                     const char* Filename,
                                     const bool Binary);
  /// @brief Write element container to VTK-file.
  /// @tparam TElementType
  /// @param rElementContainer
  /// @param Filename
  /// @param Binary
  /// @return bool
  template<typename TElementType>
  static bool WriteElementsToVTK(const ElementContainer<TElementType>& rElementContainer,
                                 const char* Filename,
                                 const bool Binary);

  /// @brief Write points to VTK. Interface for ElementContainer.
  /// @tparam TElementType
  /// @param rElementContainer
  /// @param Type
  /// @param Filename
  /// @param Binary
  /// @return bool
  template<typename TElementType>
  static bool WritePointsToVTK(const ElementContainer<TElementType>& rElementContainer,
                               const char* Type,
                               const char* Filename,
                               const bool Binary);

  /// @brief Write points to VTK.
  /// @tparam Type
  /// @param pPoints
  /// @param Filename
  /// @param Binary
  /// @return
  template<typename Type>
  static bool WritePointsToVTK(const std::vector<Type>& pPoints,
                               const char* Filename,
                               const bool Binary);

private:

  ///@}
  ///@name Private Operations
  ///@{

  template<typename T>
  static void SwapEnd(T& var)
  {
    char* varArray = reinterpret_cast<char*>(&var);
    for(long i = 0; i < static_cast<long>(sizeof(var)/2); i++)
      std::swap(varArray[sizeof(var) - 1 - i],varArray[i]);
  }

  template<typename T>
  static void WriteBinary(std::ofstream& stream, T& var){
    SwapEnd(var);
    stream.write(reinterpret_cast<char*>(&var), sizeof(T));
  }

  static bool ReadMeshFromSTL_Ascii(TriangleMeshInterface& rTriangleMesh,
                                    const char* Filename);
  static bool ReadMeshFromSTL_Binary(TriangleMeshInterface& rTriangleMesh,
                                    const char* Filename);

  static bool STLIsInASCIIFormat(const char* Filename);

  ///@}
}; // End class IO
///@} End QuESo Classes

} // End namespace queso

#endif // IO_UTILTIES_H