// Maps a Kratos element/condition/geometry type name to a VTK cell type id.
//
// The authoritative source is the Kratos core:
//   kratos/input_output/vtk_definition.cpp  (KratosGeometryType -> VTK id)
//   kratos/sources/kratos_application.cpp    (Element<dim>D<n>N -> geometry)
//
// Kratos names follow the convention `<Base><2D|3D><NodeCount>N`, e.g.
//   Element2D3N, SmallDisplacementElement3D4N, LineCondition2D2N.
// `Begin Geometries` blocks instead use bare geometry names like `Triangle2D3`
// or `Line2D2`. We handle both: first try the (dimension, nodeCount) decoding,
// then fall back to recognising the geometry-family word in the name.

/** VTK cell type ids (subset used by Kratos geometries). */
export const VtkCellType = {
  VERTEX: 1,
  LINE: 3,
  TRIANGLE: 5,
  QUAD: 9,
  TETRA: 10,
  HEXAHEDRON: 12,
  WEDGE: 13,
  PYRAMID: 14,
  QUADRATIC_EDGE: 21,
  QUADRATIC_TRIANGLE: 22,
  QUADRATIC_QUAD: 23,
  QUADRATIC_TETRA: 24,
  QUADRATIC_HEXAHEDRON: 25,
  QUADRATIC_WEDGE: 26,
  QUADRATIC_PYRAMID: 27,
  BIQUADRATIC_QUAD: 28,
  TRIQUADRATIC_HEXAHEDRON: 29,
} as const;

// (dimension, nodeCount) -> VTK cell type. Mirrors the core registration table.
const BY_DIM_NODES: Record<string, number> = {
  "2:1": VtkCellType.VERTEX,
  "2:2": VtkCellType.LINE,
  "2:3": VtkCellType.TRIANGLE,
  "2:4": VtkCellType.QUAD,
  "2:6": VtkCellType.QUADRATIC_TRIANGLE,
  "2:8": VtkCellType.QUADRATIC_QUAD,
  "2:9": VtkCellType.BIQUADRATIC_QUAD,
  "3:1": VtkCellType.VERTEX,
  "3:2": VtkCellType.LINE,
  "3:3": VtkCellType.TRIANGLE,
  "3:4": VtkCellType.TETRA,
  "3:5": VtkCellType.PYRAMID,
  "3:6": VtkCellType.WEDGE,
  "3:8": VtkCellType.HEXAHEDRON,
  "3:10": VtkCellType.QUADRATIC_TETRA,
  "3:13": VtkCellType.QUADRATIC_PYRAMID,
  "3:15": VtkCellType.QUADRATIC_WEDGE,
  "3:20": VtkCellType.QUADRATIC_HEXAHEDRON,
  "3:27": VtkCellType.TRIQUADRATIC_HEXAHEDRON,
};

// Geometry-family fallback: family word + nodeCount -> VTK cell type. Used for
// bare geometry names where the family is explicit (Triangle/Quadrilateral/...).
const BY_FAMILY_NODES: Record<string, number> = {
  "point:1": VtkCellType.VERTEX,
  "line:2": VtkCellType.LINE,
  "line:3": VtkCellType.QUADRATIC_EDGE,
  "triangle:3": VtkCellType.TRIANGLE,
  "triangle:6": VtkCellType.QUADRATIC_TRIANGLE,
  "quadrilateral:4": VtkCellType.QUAD,
  "quadrilateral:8": VtkCellType.QUADRATIC_QUAD,
  "quadrilateral:9": VtkCellType.BIQUADRATIC_QUAD,
  "tetrahedra:4": VtkCellType.TETRA,
  "tetrahedra:10": VtkCellType.QUADRATIC_TETRA,
  "hexahedra:8": VtkCellType.HEXAHEDRON,
  "hexahedra:20": VtkCellType.QUADRATIC_HEXAHEDRON,
  "hexahedra:27": VtkCellType.TRIQUADRATIC_HEXAHEDRON,
  "prism:6": VtkCellType.WEDGE,
  "prism:15": VtkCellType.QUADRATIC_WEDGE,
  "pyramid:5": VtkCellType.PYRAMID,
  "pyramid:13": VtkCellType.QUADRATIC_PYRAMID,
};

const FAMILY_WORDS = [
  "point",
  "line",
  "triangle",
  "quadrilateral",
  "tetrahedra",
  "hexahedra",
  "prism",
  "pyramid",
];

export interface DecodedName {
  dimension?: 2 | 3;
  nodeCount?: number;
  family?: string;
  vtkCellType?: number;
}

/**
 * Decode a Kratos type name into dimension / node count / family and resolve a
 * VTK cell type. Returns whatever could be determined; `vtkCellType` is
 * undefined when the name cannot be confidently mapped.
 */
export function decodeTypeName(rawName: string): DecodedName {
  const name = rawName.trim();
  const lower = name.toLowerCase();
  const result: DecodedName = {};

  const dimMatch = name.match(/([23])D/);
  if (dimMatch) {
    result.dimension = Number(dimMatch[1]) as 2 | 3;
  }

  const nodeMatch = name.match(/(\d+)N\b/) ?? name.match(/(\d+)$/);
  if (nodeMatch) {
    result.nodeCount = Number(nodeMatch[1]);
  }

  for (const fam of FAMILY_WORDS) {
    if (lower.includes(fam)) {
      result.family = fam;
      break;
    }
  }

  // Primary path: (dimension, nodeCount).
  if (result.dimension !== undefined && result.nodeCount !== undefined) {
    const hit = BY_DIM_NODES[`${result.dimension}:${result.nodeCount}`];
    if (hit !== undefined) {
      result.vtkCellType = hit;
      return result;
    }
  }

  // Fallback: (family, nodeCount).
  if (result.family !== undefined && result.nodeCount !== undefined) {
    const hit = BY_FAMILY_NODES[`${result.family}:${result.nodeCount}`];
    if (hit !== undefined) {
      result.vtkCellType = hit;
    }
  }

  return result;
}

/** Number of nodes VTK expects for a given cell type (linear connectivity). */
export function nodesForVtkCellType(vtkCellType: number): number | undefined {
  switch (vtkCellType) {
    case VtkCellType.VERTEX:
      return 1;
    case VtkCellType.LINE:
      return 2;
    case VtkCellType.TRIANGLE:
      return 3;
    case VtkCellType.QUAD:
    case VtkCellType.TETRA:
      return 4;
    case VtkCellType.PYRAMID:
      return 5;
    case VtkCellType.WEDGE:
      return 6;
    case VtkCellType.HEXAHEDRON:
      return 8;
    case VtkCellType.QUADRATIC_EDGE:
      return 3;
    case VtkCellType.QUADRATIC_TRIANGLE:
      return 6;
    case VtkCellType.QUADRATIC_QUAD:
    case VtkCellType.QUADRATIC_TETRA:
      return VtkCellType.QUADRATIC_QUAD === vtkCellType ? 8 : 10;
    case VtkCellType.BIQUADRATIC_QUAD:
      return 9;
    case VtkCellType.QUADRATIC_WEDGE:
      return 15;
    case VtkCellType.QUADRATIC_HEXAHEDRON:
      return 20;
    case VtkCellType.QUADRATIC_PYRAMID:
      return 13;
    case VtkCellType.TRIQUADRATIC_HEXAHEDRON:
      return 27;
    default:
      return undefined;
  }
}
