// Webview entry point: owns the VTK.js scene and the outline panel. Receives a
// parsed MdpaModel from the extension host, builds one renderable layer per
// entity block plus one toggleable overlay layer per SubModelPart, and wires
// the outline checkboxes (activate/deactivate layers) and labels (frame layer).

import "@kitware/vtk.js/Rendering/Profiles/Geometry";
import vtkGenericRenderWindow from "@kitware/vtk.js/Rendering/Misc/GenericRenderWindow";
import vtkActor from "@kitware/vtk.js/Rendering/Core/Actor";
import vtkMapper from "@kitware/vtk.js/Rendering/Core/Mapper";

import { MdpaModel, SubModelPart } from "../src/parser/types";
import { buildPolyData, Cell, prepareNodes, PreparedNodes } from "./meshBuilder";
import { OutlineNode, renderOutline } from "./outline";

declare function acquireVsCodeApi(): { postMessage(msg: unknown): void };
const vscode = acquireVsCodeApi();

type RGB = [number, number, number];

const PALETTE: RGB[] = [
  [0.26, 0.59, 0.98],
  [0.96, 0.62, 0.1],
  [0.2, 0.73, 0.4],
  [0.91, 0.3, 0.24],
  [0.61, 0.35, 0.71],
  [0.1, 0.74, 0.74],
  [0.85, 0.65, 0.13],
  [0.55, 0.55, 0.6],
];

interface Layer {
  id: string;
  actor: any;
  color: RGB;
  visible: boolean;
}

// --- DOM ----------------------------------------------------------------
const renderRoot = document.getElementById("render-root") as HTMLElement;
const viewport = document.getElementById("viewport") as HTMLElement;
const outlineEl = document.getElementById("outline") as HTMLElement;
const statsEl = document.getElementById("stats") as HTMLElement;

const labelsEl = document.createElement("div");
labelsEl.id = "labels";
viewport.appendChild(labelsEl);

const messageEl = document.createElement("div");
messageEl.id = "message";
viewport.appendChild(messageEl);

// --- VTK scene ----------------------------------------------------------
const grw: any = vtkGenericRenderWindow.newInstance({
  background: readThemeBackground(),
});
grw.setContainer(renderRoot);
const renderer: any = grw.getRenderer();
const renderWindow: any = grw.getRenderWindow();
const apiRW: any = grw.getApiSpecificRenderWindow
  ? grw.getApiSpecificRenderWindow()
  : grw.getOpenGLRenderWindow();

grw.resize();
new ResizeObserver(() => {
  grw.resize();
  if (showNodeIds) {
    requestLabelUpdate();
  }
}).observe(renderRoot);

// --- State --------------------------------------------------------------
let model: MdpaModel | undefined;
let prepared: PreparedNodes | undefined;
const layers = new Map<string, Layer>();
let actors: any[] = [];
let wireframe = false;
let showNodeIds = false;
const NODE_LABEL_LIMIT = 1000;

// --- Message handling ---------------------------------------------------
window.addEventListener("message", (event) => {
  const msg = event.data;
  switch (msg?.type) {
    case "model":
      model = msg.model as MdpaModel;
      buildScene();
      break;
    case "resetCamera":
      resetCamera();
      break;
    case "toggleNodeIds":
      setNodeIds(!showNodeIds);
      break;
    case "error":
      messageEl.textContent = `Parse error: ${msg.message}`;
      break;
  }
});

// --- Scene construction -------------------------------------------------
function clearScene(): void {
  for (const actor of actors) {
    renderer.removeActor(actor);
  }
  actors = [];
  layers.clear();
  labelsEl.textContent = "";
  messageEl.textContent = "";
}

function buildScene(): void {
  if (!model) {
    return;
  }
  clearScene();
  prepared = prepareNodes(model);

  // Resolve entity ids -> cells, per kind, for SubModelPart layers.
  const elementById = new Map<number, Cell>();
  const conditionById = new Map<number, Cell>();
  const geometryById = new Map<number, Cell>();
  for (const block of model.blocks) {
    const target =
      block.kind === "Elements"
        ? elementById
        : block.kind === "Conditions"
        ? conditionById
        : geometryById;
    for (const e of block.entities) {
      target.set(e.id, { cellType: block.vtkCellType, nodeIds: e.nodeIds });
    }
  }

  // One layer per entity block (visible by default).
  let colorIdx = 0;
  const blockNodes: OutlineNode[] = [];
  for (const block of model.blocks) {
    const color = PALETTE[colorIdx % PALETTE.length];
    colorIdx++;
    const cells: Cell[] = block.entities.map((e) => ({
      cellType: block.vtkCellType,
      nodeIds: e.nodeIds,
    }));
    const id = `block:${block.kind}:${block.name}`;
    const created = addLayer(id, cells, color, true);
    blockNodes.push({
      label: block.name + (block.vtkCellType === undefined ? " (?)" : ""),
      count: block.entities.length,
      layerId: created ? id : undefined,
      visible: true,
      color,
    });
  }

  // One overlay layer per SubModelPart (hidden by default), recursive.
  const partNodes: OutlineNode[] = model.subModelParts.map((p) =>
    buildPartLayer(p, elementById, conditionById, geometryById, () => {
      const c = PALETTE[colorIdx % PALETTE.length];
      colorIdx++;
      return c;
    })
  );

  const roots: OutlineNode[] = [];
  if (blockNodes.length) {
    roots.push({ label: "Mesh", section: true, children: blockNodes });
  }
  if (partNodes.length) {
    roots.push({ label: "SubModelParts", section: true, children: partNodes });
  }
  renderOutline(outlineEl, roots, {
    onToggle: (layerId, visible) => setLayerVisible(layerId, visible),
    onFocus: (layerId) => frameLayer(layerId),
  });

  renderStats();
  resetCamera();
}

function buildPartLayer(
  part: SubModelPart,
  elementById: Map<number, Cell>,
  conditionById: Map<number, Cell>,
  geometryById: Map<number, Cell>,
  nextColor: () => RGB
): OutlineNode {
  const cells: Cell[] = [];
  for (const eid of part.elementIds) {
    const c = elementById.get(eid);
    if (c) cells.push(c);
  }
  for (const cid of part.conditionIds) {
    const c = conditionById.get(cid);
    if (c) cells.push(c);
  }
  for (const gid of part.geometryIds) {
    const c = geometryById.get(gid);
    if (c) cells.push(c);
  }
  // Node-only group: render its nodes as a point cloud.
  if (cells.length === 0 && part.nodeIds.length) {
    for (const nid of part.nodeIds) {
      cells.push({ cellType: undefined, nodeIds: [nid] });
    }
  }

  const color = nextColor();
  const id = `smp:${part.path}`;
  const created = addLayer(id, cells, color, false);
  const total =
    part.elementIds.length +
    part.conditionIds.length +
    part.geometryIds.length +
    (part.elementIds.length + part.conditionIds.length + part.geometryIds.length === 0
      ? part.nodeIds.length
      : 0);

  return {
    label: part.name,
    count: total,
    layerId: created ? id : undefined,
    visible: false,
    color,
    children: part.children.map((child) =>
      buildPartLayer(child, elementById, conditionById, geometryById, nextColor)
    ),
  };
}

function addLayer(id: string, cells: Cell[], color: RGB, visible: boolean): boolean {
  if (!prepared) {
    return false;
  }
  const built = buildPolyData(prepared, cells);
  if (!built) {
    return false;
  }
  const mapper = vtkMapper.newInstance();
  mapper.setInputData(built.polyData);
  const actor = vtkActor.newInstance();
  actor.setMapper(mapper);
  const prop = actor.getProperty();
  prop.setColor(color[0], color[1], color[2]);
  prop.setEdgeVisibility(true);
  prop.setEdgeColor(color[0] * 0.5, color[1] * 0.5, color[2] * 0.5);
  prop.setPointSize(6);
  prop.setLineWidth(1.5);
  actor.setVisibility(visible);
  renderer.addActor(actor);
  actors.push(actor);
  layers.set(id, { id, actor, color, visible });
  return true;
}

// --- Interaction --------------------------------------------------------
function setLayerVisible(layerId: string, visible: boolean): void {
  const layer = layers.get(layerId);
  if (!layer) {
    return;
  }
  layer.visible = visible;
  layer.actor.setVisibility(visible);
  renderWindow.render();
}

function frameLayer(layerId: string): void {
  const layer = layers.get(layerId);
  if (!layer) {
    return;
  }
  const bounds = layer.actor.getBounds();
  if (bounds && bounds[0] <= bounds[1]) {
    renderer.resetCamera(bounds);
    renderWindow.render();
    if (showNodeIds) {
      requestLabelUpdate();
    }
  }
}

function resetCamera(): void {
  renderer.resetCamera();
  renderWindow.render();
  if (showNodeIds) {
    requestLabelUpdate();
  }
}

function setWireframe(on: boolean): void {
  wireframe = on;
  // Representation: 1 = wireframe, 2 = surface.
  for (const layer of layers.values()) {
    layer.actor.getProperty().setRepresentation(on ? 1 : 2);
  }
  renderWindow.render();
}

// --- Node id labels -----------------------------------------------------
let labelFrame: number | undefined;
function requestLabelUpdate(): void {
  if (labelFrame !== undefined) {
    return;
  }
  labelFrame = requestAnimationFrame(() => {
    labelFrame = undefined;
    updateNodeLabels();
  });
}

function setNodeIds(on: boolean): void {
  showNodeIds = on;
  const btn = document.querySelector('#toolbar button[data-action="nodeIds"]');
  btn?.classList.toggle("active", on);
  labelsEl.textContent = "";
  if (!on || !model) {
    messageEl.textContent = "";
    stopLabelLoop();
    return;
  }
  if (model.nodes.length > NODE_LABEL_LIMIT) {
    messageEl.textContent = `Node IDs hidden: ${model.nodes.length} nodes exceed the ${NODE_LABEL_LIMIT} label limit.`;
    showNodeIds = false;
    btn?.classList.remove("active");
    return;
  }
  // Pre-create one label element per node.
  for (const n of model.nodes) {
    const el = document.createElement("div");
    el.className = "node-label";
    el.textContent = String(n.id);
    el.dataset.x = String(n.x);
    el.dataset.y = String(n.y);
    el.dataset.z = String(n.z);
    labelsEl.appendChild(el);
  }
  startLabelLoop();
}

let labelLoop: number | undefined;
function startLabelLoop(): void {
  const tick = () => {
    updateNodeLabels();
    labelLoop = requestAnimationFrame(tick);
  };
  labelLoop = requestAnimationFrame(tick);
}
function stopLabelLoop(): void {
  if (labelLoop !== undefined) {
    cancelAnimationFrame(labelLoop);
    labelLoop = undefined;
  }
}

function updateNodeLabels(): void {
  if (!showNodeIds) {
    return;
  }
  const size = apiRW.getSize(); // device pixels [w, h]
  const dpr = window.devicePixelRatio || 1;
  const children = labelsEl.children;
  for (let i = 0; i < children.length; i++) {
    const el = children[i] as HTMLElement;
    const x = Number(el.dataset.x);
    const y = Number(el.dataset.y);
    const z = Number(el.dataset.z);
    const disp = apiRW.worldToDisplay(x, y, z, renderer); // origin bottom-left
    const left = disp[0] / dpr;
    const top = (size[1] - disp[1]) / dpr;
    el.style.left = `${left}px`;
    el.style.top = `${top}px`;
  }
}

// --- Stats panel --------------------------------------------------------
function renderStats(): void {
  if (!model) {
    return;
  }
  const count = (kind: string) =>
    model!.blocks
      .filter((b) => b.kind === kind)
      .reduce((s, b) => s + b.entities.length, 0);
  const unmapped = model.blocks.filter((b) => b.vtkCellType === undefined);
  const b = model.bounds;
  const fmt = (v: number) => (Number.isFinite(v) ? v.toPrecision(4) : "0");

  const rows: string[] = [
    row("Nodes", String(model.nodes.length)),
    row("Elements", String(count("Elements"))),
    row("Conditions", String(count("Conditions"))),
    row("Geometries", String(count("Geometries"))),
    row("SubModelParts", String(countParts(model.subModelParts))),
    row("Dimensionality", model.is3D ? "3D" : "2D"),
    row(
      "Bounds",
      `[${fmt(b.min[0])}, ${fmt(b.min[1])}, ${fmt(b.min[2])}] – [${fmt(
        b.max[0]
      )}, ${fmt(b.max[1])}, ${fmt(b.max[2])}]`
    ),
  ];
  if (unmapped.length) {
    rows.push(
      `<div class="stat-row warn"><span class="stat-key">Unmapped types</span><span>${unmapped
        .map((u) => u.name)
        .join(", ")}</span></div>`
    );
  }
  if (model.diagnostics.length) {
    rows.push(
      `<div class="stat-row warn"><span class="stat-key">Warnings</span><span>${model.diagnostics.length}</span></div>`
    );
  }
  statsEl.innerHTML = rows.join("");
}

function row(key: string, value: string): string {
  return `<div class="stat-row"><span class="stat-key">${key}</span><span>${value}</span></div>`;
}

function countParts(parts: SubModelPart[]): number {
  let n = parts.length;
  for (const p of parts) {
    n += countParts(p.children);
  }
  return n;
}

// --- Toolbar ------------------------------------------------------------
document.getElementById("toolbar")?.addEventListener("click", (e) => {
  const target = e.target as HTMLElement;
  const action = target.dataset.action;
  if (action === "reset") {
    resetCamera();
  } else if (action === "wireframe") {
    setWireframe(!wireframe);
    target.classList.toggle("active", wireframe);
  } else if (action === "nodeIds") {
    setNodeIds(!showNodeIds);
  }
});

// --- Helpers ------------------------------------------------------------
function readThemeBackground(): RGB {
  const css = getComputedStyle(document.body).backgroundColor;
  const m = css.match(/rgba?\(([^)]+)\)/);
  if (m) {
    const parts = m[1].split(",").map((s) => parseFloat(s.trim()));
    if (parts.length >= 3) {
      return [parts[0] / 255, parts[1] / 255, parts[2] / 255];
    }
  }
  return [0.12, 0.12, 0.14];
}

// Tell the host we are ready to receive the model.
vscode.postMessage({ type: "ready" });
