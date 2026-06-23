// Tolerant, line-oriented state-machine parser for Kratos `.mdpa` files.
//
// The format is whitespace-agnostic and block-structured: `Begin <Type> [args]`
// ... `End <Type>`. Blocks nest (SubModelPart within SubModelPart, Table within
// Properties). On any malformed line we record a diagnostic and keep going
// rather than throwing, so a partially-broken file still previews.
//
// See docs/pages/Kratos/For_Developers/IO/Input-data.md for the format.

import {
  EntityBlock,
  EntityKind,
  MdpaModel,
  MdpaNode,
  MetaBlock,
  SubModelPart,
} from "./types";
import { decodeTypeName } from "./geometryMap";

type SubListKey =
  | "nodeIds"
  | "elementIds"
  | "conditionIds"
  | "geometryIds"
  | "constraintIds";

interface Frame {
  type: string;
  entityBlock?: EntityBlock;
  subModelPart?: SubModelPart;
  listTarget?: { part: SubModelPart; key: SubListKey };
  meta?: MetaBlock;
}

const ENTITY_KINDS: Record<string, EntityKind> = {
  Elements: "Elements",
  Conditions: "Conditions",
  Geometries: "Geometries",
};

const SUBLIST_KEYS: Record<string, SubListKey> = {
  SubModelPartNodes: "nodeIds",
  SubModelPartElements: "elementIds",
  SubModelPartConditions: "conditionIds",
  SubModelPartGeometries: "geometryIds",
  SubModelPartConstraints: "constraintIds",
};

// Blocks treated as opaque metadata: surfaced in the stats panel, not rendered.
const META_TYPES = new Set([
  "ModelPartData",
  "Properties",
  "Table",
  "NodalData",
  "ElementalData",
  "ConditionalData",
  "Constraints",
  "Mesh",
  "MeshData",
  "MeshNodes",
  "MeshElements",
  "MeshConditions",
  "SubModelPartData",
  "SubModelPartTables",
  "SubModelPartProperties",
]);

/** Strip a `//` line comment (the format has no block comments). */
function stripComment(line: string): string {
  const idx = line.indexOf("//");
  return idx === -1 ? line : line.slice(0, idx);
}

export function parseMdpa(text: string): MdpaModel {
  const nodes: MdpaNode[] = [];
  const blocks: EntityBlock[] = [];
  const subModelParts: SubModelPart[] = [];
  const meta: MetaBlock[] = [];
  const diagnostics: { line: number; message: string }[] = [];

  // Reuse one EntityBlock per (kind, name) so repeated blocks form one layer.
  const blockIndex = new Map<string, EntityBlock>();
  const stack: Frame[] = [];

  const lines = text.split(/\r?\n/);

  const topSubModelPart = (): SubModelPart | undefined => {
    for (let i = stack.length - 1; i >= 0; i--) {
      if (stack[i].subModelPart) {
        return stack[i].subModelPart;
      }
    }
    return undefined;
  };

  for (let i = 0; i < lines.length; i++) {
    const lineNo = i + 1;
    const stripped = stripComment(lines[i]).trim();
    if (stripped.length === 0) {
      continue;
    }
    const tokens = stripped.split(/\s+/);
    const head = tokens[0];

    if (head === "Begin") {
      const blockType = tokens[1];
      const args = tokens.slice(2);
      if (!blockType) {
        diagnostics.push({ line: lineNo, message: "`Begin` without a block type." });
        stack.push({ type: "<unknown>" });
        continue;
      }

      if (blockType === "Nodes") {
        stack.push({ type: "Nodes" });
      } else if (ENTITY_KINDS[blockType]) {
        const kind = ENTITY_KINDS[blockType];
        const name = args[0] ?? "<unnamed>";
        const key = `${kind}::${name}`;
        let block = blockIndex.get(key);
        if (!block) {
          const decoded = decodeTypeName(name);
          block = {
            kind,
            name,
            vtkCellType: decoded.vtkCellType,
            entities: [],
          };
          blockIndex.set(key, block);
          blocks.push(block);
        }
        stack.push({ type: blockType, entityBlock: block });
      } else if (blockType === "SubModelPart") {
        const name = args[0] ?? "<unnamed>";
        const parent = topSubModelPart();
        const part: SubModelPart = {
          name,
          nodeIds: [],
          elementIds: [],
          conditionIds: [],
          geometryIds: [],
          constraintIds: [],
          path: parent ? `${parent.path}/${name}` : name,
          children: [],
        };
        if (parent) {
          parent.children.push(part);
        } else {
          subModelParts.push(part);
        }
        stack.push({ type: "SubModelPart", subModelPart: part });
      } else if (SUBLIST_KEYS[blockType]) {
        const part = topSubModelPart();
        if (!part) {
          diagnostics.push({
            line: lineNo,
            message: `${blockType} outside any SubModelPart.`,
          });
          stack.push({ type: blockType });
        } else {
          stack.push({
            type: blockType,
            listTarget: { part, key: SUBLIST_KEYS[blockType] },
          });
        }
      } else if (META_TYPES.has(blockType)) {
        const label = args.length ? `${blockType} ${args.join(" ")}` : blockType;
        const metaBlock: MetaBlock = { label, lineCount: 0 };
        meta.push(metaBlock);
        stack.push({ type: blockType, meta: metaBlock });
      } else {
        // Unknown block: keep nesting balanced, note it once.
        diagnostics.push({
          line: lineNo,
          message: `Unknown block type "${blockType}"; contents ignored.`,
        });
        stack.push({ type: blockType });
      }
      continue;
    }

    if (head === "End") {
      const endType = tokens[1];
      const frame = stack.pop();
      if (!frame) {
        diagnostics.push({ line: lineNo, message: `Stray "End ${endType ?? ""}".` });
      } else if (endType && frame.type !== endType && frame.type !== "<unknown>") {
        diagnostics.push({
          line: lineNo,
          message: `"End ${endType}" does not match open block "${frame.type}".`,
        });
      }
      continue;
    }

    // Data line: route by the current open block.
    const frame = stack[stack.length - 1];
    if (!frame) {
      diagnostics.push({ line: lineNo, message: "Data line outside any block." });
      continue;
    }

    if (frame.type === "Nodes") {
      if (tokens.length < 4) {
        diagnostics.push({ line: lineNo, message: "Node line needs id X Y Z." });
        continue;
      }
      nodes.push({
        id: parseInt(tokens[0], 10),
        x: Number(tokens[1]),
        y: Number(tokens[2]),
        z: Number(tokens[3]),
      });
    } else if (frame.entityBlock) {
      const id = parseInt(tokens[0], 10);
      if (frame.entityBlock.kind === "Geometries") {
        frame.entityBlock.entities.push({
          id,
          nodeIds: tokens.slice(1).map((t) => parseInt(t, 10)),
        });
      } else {
        frame.entityBlock.entities.push({
          id,
          propertyId: tokens.length > 1 ? parseInt(tokens[1], 10) : undefined,
          nodeIds: tokens.slice(2).map((t) => parseInt(t, 10)),
        });
      }
    } else if (frame.listTarget) {
      const id = parseInt(tokens[0], 10);
      if (!Number.isNaN(id)) {
        frame.listTarget.part[frame.listTarget.key].push(id);
      }
    } else if (frame.meta) {
      frame.meta.lineCount++;
    }
    // Frames for SubModelPart containers / unknown blocks ignore stray data.
  }

  if (stack.length > 0) {
    diagnostics.push({
      line: lines.length,
      message: `${stack.length} block(s) not closed by end of file.`,
    });
  }

  // Bounds + dimensionality.
  let is3D = false;
  const min: [number, number, number] = [Infinity, Infinity, Infinity];
  const max: [number, number, number] = [-Infinity, -Infinity, -Infinity];
  for (const n of nodes) {
    const c: [number, number, number] = [n.x, n.y, n.z];
    for (let d = 0; d < 3; d++) {
      if (c[d] < min[d]) min[d] = c[d];
      if (c[d] > max[d]) max[d] = c[d];
    }
    if (Math.abs(n.z) > 1e-12) {
      is3D = true;
    }
  }
  if (nodes.length === 0) {
    min[0] = min[1] = min[2] = 0;
    max[0] = max[1] = max[2] = 0;
  }

  return { nodes, blocks, subModelParts, meta, diagnostics, is3D, bounds: { min, max } };
}
