// Data model produced by the MDPA parser and consumed by the webview viewer.
// Kept free of any `vscode` or DOM dependency so it can be unit-tested under
// plain Node and reused on both the extension-host and webview sides.

/** A mesh node: integer id plus Cartesian coordinates. */
export interface MdpaNode {
  id: number;
  x: number;
  y: number;
  z: number;
}

/**
 * A single entity (element/condition/geometry) record: its id, optional
 * property id (geometries carry none), and the ordered connectivity node ids.
 */
export interface MdpaEntity {
  id: number;
  propertyId?: number;
  nodeIds: number[];
}

/** Kind of renderable entity block. */
export type EntityKind = "Elements" | "Conditions" | "Geometries";

/**
 * A `Begin Elements/Conditions/Geometries <Name>` block: all entities sharing
 * one Kratos type name, together with the resolved VTK cell type for rendering.
 */
export interface EntityBlock {
  kind: EntityKind;
  /** The Kratos type token, e.g. `Element2D3N`, `LineCondition2D2N`, `Triangle2D3`. */
  name: string;
  /** Resolved VTK cell type id, or undefined when the name could not be mapped. */
  vtkCellType?: number;
  entities: MdpaEntity[];
}

/**
 * A (possibly nested) SubModelPart. Only id references are stored here — the
 * actual geometry lives in the root entity blocks / node table. This is the
 * tree the outline renders and whose entries become toggleable layers.
 */
export interface SubModelPart {
  name: string;
  nodeIds: number[];
  elementIds: number[];
  conditionIds: number[];
  geometryIds: number[];
  constraintIds: number[];
  /** Hierarchical path from the root, e.g. `Inlets/Inlet1`. Unique per part. */
  path: string;
  children: SubModelPart[];
}

/** A non-rendered metadata block surfaced in the stats/outline panel. */
export interface MetaBlock {
  /** e.g. `ModelPartData`, `Properties 1`, `Table 1`, `NodalData DISPLACEMENT_X`. */
  label: string;
  /** Number of data lines inside the block. */
  lineCount: number;
}

/** A non-fatal problem encountered while parsing. */
export interface MdpaDiagnostic {
  line: number; // 1-based
  message: string;
}

/** The fully parsed model. */
export interface MdpaModel {
  nodes: MdpaNode[];
  blocks: EntityBlock[];
  subModelParts: SubModelPart[];
  meta: MetaBlock[];
  diagnostics: MdpaDiagnostic[];
  /** True when any node has a non-zero Z coordinate. */
  is3D: boolean;
  bounds: {
    min: [number, number, number];
    max: [number, number, number];
  };
}
