// Converts parsed MDPA entities into vtk.js PolyData for rendering.
//
// Surface cells (triangle/quad) become polygons; line cells become lines;
// points/unknown cells become vertices. Volume cells (tet/hex/wedge/pyramid)
// are reduced to their boundary surface: each cell contributes its faces, and
// faces shared by two cells cancel out, leaving only the outer skin. Quadratic
// elements are approximated by their corner nodes (good enough for preview).

import vtkPolyData from "@kitware/vtk.js/Common/DataModel/PolyData";
import { MdpaModel } from "../src/parser/types";
import { VtkCellType } from "../src/parser/geometryMap";

export interface PreparedNodes {
  /** Original node id -> contiguous base index into `coords`. */
  index: Map<number, number>;
  /** Flat x,y,z per node. */
  coords: Float32Array;
}

export function prepareNodes(model: MdpaModel): PreparedNodes {
  const index = new Map<number, number>();
  const coords = new Float32Array(model.nodes.length * 3);
  model.nodes.forEach((n, i) => {
    index.set(n.id, i);
    coords[i * 3] = n.x;
    coords[i * 3 + 1] = n.y;
    coords[i * 3 + 2] = n.z;
  });
  return { index, coords };
}

export interface Cell {
  cellType?: number;
  nodeIds: number[];
}

type Category = "point" | "line" | "surface" | "volume" | "unknown";
interface Topo {
  corners: number;
  category: Category;
  faces?: number[][];
}

const C = VtkCellType;

// Linear corner-node face templates (indices into the corner list).
const TET_FACES = [
  [0, 1, 2],
  [0, 3, 1],
  [0, 2, 3],
  [1, 3, 2],
];
const HEX_FACES = [
  [0, 1, 2, 3],
  [4, 7, 6, 5],
  [0, 4, 5, 1],
  [1, 5, 6, 2],
  [2, 6, 7, 3],
  [3, 7, 4, 0],
];
const WEDGE_FACES = [
  [0, 1, 2],
  [3, 5, 4],
  [0, 3, 4, 1],
  [1, 4, 5, 2],
  [2, 5, 3, 0],
];
const PYRAMID_FACES = [
  [0, 3, 2, 1],
  [0, 1, 4],
  [1, 2, 4],
  [2, 3, 4],
  [3, 0, 4],
];

function topo(cellType?: number): Topo {
  switch (cellType) {
    case C.VERTEX:
      return { corners: 1, category: "point" };
    case C.LINE:
    case C.QUADRATIC_EDGE:
      return { corners: 2, category: "line" };
    case C.TRIANGLE:
    case C.QUADRATIC_TRIANGLE:
      return { corners: 3, category: "surface" };
    case C.QUAD:
    case C.QUADRATIC_QUAD:
    case C.BIQUADRATIC_QUAD:
      return { corners: 4, category: "surface" };
    case C.TETRA:
    case C.QUADRATIC_TETRA:
      return { corners: 4, category: "volume", faces: TET_FACES };
    case C.HEXAHEDRON:
    case C.QUADRATIC_HEXAHEDRON:
    case C.TRIQUADRATIC_HEXAHEDRON:
      return { corners: 8, category: "volume", faces: HEX_FACES };
    case C.WEDGE:
    case C.QUADRATIC_WEDGE:
      return { corners: 6, category: "volume", faces: WEDGE_FACES };
    case C.PYRAMID:
    case C.QUADRATIC_PYRAMID:
      return { corners: 5, category: "volume", faces: PYRAMID_FACES };
    default:
      return { corners: 0, category: "unknown" };
  }
}

export interface BuiltMesh {
  polyData: ReturnType<typeof vtkPolyData.newInstance>;
}

/**
 * Build one PolyData from a list of cells, with its own compact point array so
 * that `actor.getBounds()` reflects only this layer (used for framing). Returns
 * null when the cells reference no valid geometry.
 */
export function buildPolyData(prep: PreparedNodes, cells: Cell[]): BuiltMesh | null {
  const localPoints: number[] = [];
  const localIndex = new Map<number, number>();
  const polys: number[] = [];
  const lines: number[] = [];
  const verts: number[] = [];

  const localOf = (id: number): number | undefined => {
    const cached = localIndex.get(id);
    if (cached !== undefined) {
      return cached;
    }
    const base = prep.index.get(id);
    if (base === undefined) {
      return undefined;
    }
    const li = localPoints.length / 3;
    localPoints.push(
      prep.coords[base * 3],
      prep.coords[base * 3 + 1],
      prep.coords[base * 3 + 2]
    );
    localIndex.set(id, li);
    return li;
  };

  // Boundary-face accumulation for volume cells (keyed by sorted global ids).
  const faceIds = new Map<string, number[]>();
  const faceCount = new Map<string, number>();

  for (const cell of cells) {
    const t = topo(cell.cellType);

    if (cell.cellType === undefined || t.category === "unknown") {
      for (const id of cell.nodeIds) {
        const li = localOf(id);
        if (li !== undefined) {
          verts.push(1, li);
        }
      }
      continue;
    }

    const corners = cell.nodeIds.slice(0, t.corners);
    if (
      corners.length < t.corners ||
      corners.some((id) => prep.index.get(id) === undefined)
    ) {
      continue; // incomplete or dangling connectivity
    }

    if (t.category === "point") {
      const li = localOf(corners[0]);
      if (li !== undefined) {
        verts.push(1, li);
      }
    } else if (t.category === "line") {
      const a = localOf(corners[0]);
      const b = localOf(corners[1]);
      if (a !== undefined && b !== undefined) {
        lines.push(2, a, b);
      }
    } else if (t.category === "surface") {
      const lis = corners.map(localOf) as number[];
      polys.push(lis.length, ...lis);
    } else if (t.category === "volume" && t.faces) {
      for (const face of t.faces) {
        const ids = face.map((fi) => corners[fi]);
        const key = [...ids].sort((p, q) => p - q).join(",");
        faceCount.set(key, (faceCount.get(key) ?? 0) + 1);
        if (!faceIds.has(key)) {
          faceIds.set(key, ids);
        }
      }
    }
  }

  for (const [key, ids] of faceIds) {
    if (faceCount.get(key) === 1) {
      const lis = ids.map(localOf) as number[];
      polys.push(lis.length, ...lis);
    }
  }

  if (polys.length === 0 && lines.length === 0 && verts.length === 0) {
    return null;
  }

  const polyData = vtkPolyData.newInstance();
  polyData.getPoints().setData(Float32Array.from(localPoints), 3);
  if (polys.length) {
    polyData.getPolys().setData(Uint32Array.from(polys));
  }
  if (lines.length) {
    polyData.getLines().setData(Uint32Array.from(lines));
  }
  if (verts.length) {
    polyData.getVerts().setData(Uint32Array.from(verts));
  }
  return { polyData };
}
