/*   ______     _____ _           ________
    / ____/___ / ___/(_)___ ___  /  _/ __ |
   / /   / __ \\__ \/ / __ `__ \ / // / / /
  / /___/ /_/ /__/ / / / / / / // // /_/ /
  \____/\____/____/_/_/ /_/ /_/___/\____/
  Kratos CoSimulationApplication

  License:         BSD License, see license.txt

  Main authors:    Philipp Bucher (https://github.com/philbucher)
*/

#ifndef CO_SIM_IO_C_MODEL_PART_INCLUDED
#define CO_SIM_IO_C_MODEL_PART_INCLUDED

#include "define_c.h"

typedef struct CoSimIO_Node
{
    void* PtrCppNode;
} CoSimIO_Node;

typedef struct CoSimIO_Element
{
    void* PtrCppElement;
} CoSimIO_Element;

typedef struct CoSimIO_ModelPart
{
    void* PtrCppModelPart;
} CoSimIO_ModelPart;

typedef enum
{
    CoSimIO_Hexahedra3D20,
    CoSimIO_Hexahedra3D27,
    CoSimIO_Hexahedra3D8,
    CoSimIO_Prism3D15,
    CoSimIO_Prism3D6,
    CoSimIO_Pyramid3D13,
    CoSimIO_Pyramid3D5,
    CoSimIO_Quadrilateral2D4,
    CoSimIO_Quadrilateral2D8,
    CoSimIO_Quadrilateral2D9,
    CoSimIO_Quadrilateral3D4,
    CoSimIO_Quadrilateral3D8,
    CoSimIO_Quadrilateral3D9,
    CoSimIO_Tetrahedra3D10,
    CoSimIO_Tetrahedra3D4,
    CoSimIO_Triangle2D3,
    CoSimIO_Triangle2D6,
    CoSimIO_Triangle3D3,
    CoSimIO_Triangle3D6,
    CoSimIO_Line2D2,
    CoSimIO_Line2D3,
    CoSimIO_Line3D2,
    CoSimIO_Line3D3,
    CoSimIO_Point2D,
    CoSimIO_Point3D
} CoSimIO_ElementType;


/* Node functions */
int CoSimIO_Node_Id(CoSimIO_Node I_Node);
double CoSimIO_Node_X(CoSimIO_Node I_Node);
double CoSimIO_Node_Y(CoSimIO_Node I_Node);
double CoSimIO_Node_Z(CoSimIO_Node I_Node);
double CoSimIO_Node_Coordinate(CoSimIO_Node I_Node, const int I_Index);


/* Element functions */
int CoSimIO_Element_Id(CoSimIO_Element I_Element);
CoSimIO_ElementType CoSimIO_Element_Type(CoSimIO_Element I_Element);
int CoSimIO_Element_NumberOfNodes(CoSimIO_Element I_Element);
CoSimIO_Node CoSimIO_Element_GetNodeByIndex(CoSimIO_Element I_Element, const int I_Index);

/* ModelPart functions */
CO_SIM_IO_NODISCARD CoSimIO_ModelPart CoSimIO_CreateModelPart(const char* I_Name);

int CoSimIO_FreeModelPart(CoSimIO_ModelPart I_ModelPart);

const char* CoSimIO_ModelPart_Name(CoSimIO_ModelPart I_ModelPart);
int CoSimIO_ModelPart_NumberOfNodes(CoSimIO_ModelPart I_ModelPart);
int CoSimIO_ModelPart_NumberOfLocalNodes(CoSimIO_ModelPart I_ModelPart);
int CoSimIO_ModelPart_NumberOfGhostNodes(CoSimIO_ModelPart I_ModelPart);
int CoSimIO_ModelPart_NumberOfElements(CoSimIO_ModelPart I_ModelPart);
CoSimIO_Node CoSimIO_ModelPart_GetNodeByIndex(CoSimIO_ModelPart I_ModelPart, const int I_Index);
CoSimIO_Node CoSimIO_ModelPart_GetLocalNodeByIndex(CoSimIO_ModelPart I_ModelPart, const int I_Index);
CoSimIO_Node CoSimIO_ModelPart_GetGhostNodeByIndex(CoSimIO_ModelPart I_ModelPart, const int I_Index);
CoSimIO_Node CoSimIO_ModelPart_GetNodeById(CoSimIO_ModelPart I_ModelPart, const int I_Id);
CoSimIO_Element CoSimIO_ModelPart_GetElementByIndex(CoSimIO_ModelPart I_ModelPart, const int I_Index);
CoSimIO_Element CoSimIO_ModelPart_GetElementById(CoSimIO_ModelPart I_ModelPart, const int I_Id);
void CoSimIO_ModelPart_Clear(CoSimIO_ModelPart I_ModelPart);

CoSimIO_Node CoSimIO_ModelPart_CreateNewNode(
    CoSimIO_ModelPart I_ModelPart,
    const int I_Id,
    const double I_X,
    const double I_Y,
    const double I_Z);

void CoSimIO_ModelPart_CreateNewNodes(
    CoSimIO_ModelPart I_ModelPart,
    const int I_NumberOfNodes,
    const int* I_Id,
    const double* I_X,
    const double* I_Y,
    const double* I_Z);

CoSimIO_Node CoSimIO_ModelPart_CreateNewGhostNode(
    CoSimIO_ModelPart I_ModelPart,
    const int I_Id,
    const double I_X,
    const double I_Y,
    const double I_Z,
    const int PartitionIndex);

void CoSimIO_ModelPart_CreateNewGhostNodes(
    CoSimIO_ModelPart I_ModelPart,
    const int I_NumberOfNodes,
    const int* I_Id,
    const double* I_X,
    const double* I_Y,
    const double* I_Z,
    const int* PartitionIndex);

CoSimIO_Element CoSimIO_ModelPart_CreateNewElement(
    CoSimIO_ModelPart I_ModelPart,
    const int I_Id,
    const CoSimIO_ElementType I_Type,
    const int* I_Connectivities);

void CoSimIO_ModelPart_CreateNewElements(
    CoSimIO_ModelPart I_ModelPart,
    const int I_NumberOfElements,
    const int* I_Id,
    const CoSimIO_ElementType* I_Type,
    const int I_NumberOfConnectivities,
    const int* I_Connectivities);


#endif /* CO_SIM_IO_C_MODEL_PART_INCLUDED */
