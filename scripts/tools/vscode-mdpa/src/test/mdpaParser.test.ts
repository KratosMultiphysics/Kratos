import { test } from "node:test";
import assert from "node:assert/strict";
import * as fs from "node:fs";
import * as path from "node:path";
import { parseMdpa } from "../parser/mdpaParser";
import { EntityBlock, SubModelPart } from "../parser/types";
import { VtkCellType } from "../parser/geometryMap";

// out/test -> out -> vscode-mdpa -> tools -> scripts -> repo root
const REPO_ROOT = path.resolve(__dirname, "../../../../../");

function fixture(rel: string): string {
  return fs.readFileSync(path.join(REPO_ROOT, rel), "utf8");
}

function findBlock(blocks: EntityBlock[], name: string): EntityBlock | undefined {
  return blocks.find((b) => b.name === name);
}

function findPart(parts: SubModelPart[], name: string): SubModelPart | undefined {
  for (const p of parts) {
    if (p.name === name) return p;
    const child = findPart(p.children, name);
    if (child) return child;
  }
  return undefined;
}

const READ_FIXTURE =
  "kratos/tests/auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_read.mdpa";

test("parses nodes, elements, conditions and geometries", () => {
  const model = parseMdpa(fixture(READ_FIXTURE));

  assert.equal(model.nodes.length, 6);

  const elements = findBlock(model.blocks, "Element2D3N");
  assert.ok(elements, "Element2D3N block present");
  assert.equal(elements!.entities.length, 4);
  assert.equal(elements!.vtkCellType, VtkCellType.TRIANGLE);
  assert.deepEqual(elements!.entities[0].nodeIds, [1, 2, 3]);
  assert.equal(elements!.entities[0].propertyId, 1);

  const conditions = findBlock(model.blocks, "LineCondition2D2N");
  assert.ok(conditions);
  assert.equal(conditions!.entities.length, 5);
  assert.equal(conditions!.vtkCellType, VtkCellType.LINE);

  const tris = findBlock(model.blocks, "Triangle2D3");
  const lines = findBlock(model.blocks, "Line2D2");
  assert.equal(tris!.entities.length, 4);
  assert.equal(lines!.entities.length, 5);
  // Geometries carry no property id.
  assert.equal(tris!.entities[0].propertyId, undefined);
  assert.deepEqual(tris!.entities[0].nodeIds, [1, 2, 3]);
});

test("builds the nested SubModelPart tree", () => {
  const model = parseMdpa(fixture(READ_FIXTURE));

  const inlets = findPart(model.subModelParts, "Inlets");
  assert.ok(inlets, "Inlets present at root");
  assert.equal(inlets!.path, "Inlets");
  assert.deepEqual(inlets!.nodeIds, [1, 2]);
  assert.deepEqual(inlets!.elementIds, [1]);
  assert.deepEqual(inlets!.conditionIds, [1, 1800]);
  assert.equal(inlets!.children.length, 2);

  const inlet1 = inlets!.children.find((c) => c.name === "Inlet1");
  assert.ok(inlet1);
  assert.equal(inlet1!.path, "Inlets/Inlet1");
  assert.deepEqual(inlet1!.nodeIds, [1, 3]);

  const inlet2 = inlets!.children.find((c) => c.name === "Inlet2");
  assert.ok(inlet2);
  assert.deepEqual(inlet2!.conditionIds, [1800, 1801]);

  const outlet = findPart(model.subModelParts, "Outlet");
  assert.ok(outlet, "Outlet present at root");
  assert.deepEqual(outlet!.conditionIds, [1948]);
});

test("detects 2D model and computes bounds", () => {
  const model = parseMdpa(fixture(READ_FIXTURE));
  assert.equal(model.is3D, false);
  assert.equal(model.bounds.min[2], 0);
  assert.equal(model.bounds.max[2], 0);
  assert.ok(model.bounds.max[0] >= model.bounds.min[0]);
});

test("produces no fatal diagnostics on the reference file", () => {
  const model = parseMdpa(fixture(READ_FIXTURE));
  const unclosed = model.diagnostics.filter((d) => /not closed|Stray|outside any/.test(d.message));
  assert.deepEqual(unclosed, []);
});

test("parses additional real fixtures without unbalanced blocks", () => {
  const files = [
    "applications/MetisApplication/tests/quads.mdpa",
    "applications/OptimizationApplication/tests/mdpas/shell.mdpa",
    "applications/MeshingApplication/tests/cube_with_5_faces.mdpa",
  ];
  for (const f of files) {
    const model = parseMdpa(fixture(f));
    assert.ok(model.nodes.length > 0, `${f}: has nodes`);
    assert.ok(model.blocks.length > 0, `${f}: has entity blocks`);
    const unbalanced = model.diagnostics.filter((d) => /not closed|Stray/.test(d.message));
    assert.deepEqual(unbalanced, [], `${f}: balanced Begin/End`);
  }
});

test("handles flag-form NodalData and irregular whitespace", () => {
  const text = [
    "Begin Nodes",
    "\t1\t  0.0   0.0\t0.0",
    "  2 1.0 0.0 0.0",
    "End Nodes",
    "Begin NodalData BOUNDARY",
    "1",
    "2",
    "End NodalData",
    "Begin NodalData DISPLACEMENT_X",
    "1 1 0.5",
    "End NodalData",
  ].join("\n");
  const model = parseMdpa(text);
  assert.equal(model.nodes.length, 2);
  assert.equal(model.nodes[0].id, 1);
  const boundary = model.meta.find((m) => m.label.includes("BOUNDARY"));
  assert.ok(boundary);
  assert.equal(boundary!.lineCount, 2);
});

test("performance: parses the large cube fixture quickly", () => {
  const text = fixture("applications/MetisApplication/tests/cube.mdpa");
  const start = Date.now();
  const model = parseMdpa(text);
  const elapsed = Date.now() - start;
  const totalEntities = model.blocks.reduce((sum, b) => sum + b.entities.length, 0);
  assert.ok(model.nodes.length > 100, "node count");
  assert.ok(totalEntities > 1000, "large entity count");
  assert.ok(elapsed < 2000, `parsed in ${elapsed}ms`);
});
