// Renders the ModelPart / SubModelPart outline as a tree of toggleable layers.
// Pure DOM, no framework. Each node with a `layerId` gets a visibility checkbox
// (the requested activate/deactivate of a layer) and a clickable label that
// frames that layer in the 3D view.

export interface OutlineNode {
  label: string;
  count?: number;
  /** Present when this row maps to a renderable layer/actor. */
  layerId?: string;
  visible?: boolean;
  /** RGB in 0..1 for the layer swatch. */
  color?: [number, number, number];
  /** True for non-toggleable section headers (e.g. "Mesh", "SubModelParts"). */
  section?: boolean;
  children?: OutlineNode[];
}

export interface OutlineHandlers {
  onToggle(layerId: string, visible: boolean): void;
  onFocus(layerId: string): void;
}

function rgbToCss(c: [number, number, number]): string {
  const to255 = (v: number) => Math.round(Math.max(0, Math.min(1, v)) * 255);
  return `rgb(${to255(c[0])}, ${to255(c[1])}, ${to255(c[2])})`;
}

export function renderOutline(
  container: HTMLElement,
  roots: OutlineNode[],
  handlers: OutlineHandlers
): void {
  container.textContent = "";
  for (const node of roots) {
    container.appendChild(buildNode(node, 0, handlers));
  }
}

function buildNode(
  node: OutlineNode,
  depth: number,
  handlers: OutlineHandlers
): HTMLElement {
  const wrapper = document.createElement("div");

  const row = document.createElement("div");
  row.className = node.section ? "outline-row outline-section" : "outline-row";
  row.style.paddingLeft = `${depth * 14}px`;

  if (node.layerId) {
    const checkbox = document.createElement("input");
    checkbox.type = "checkbox";
    checkbox.checked = node.visible ?? true;
    checkbox.title = "Show / hide layer";
    checkbox.addEventListener("change", () =>
      handlers.onToggle(node.layerId!, checkbox.checked)
    );
    row.appendChild(checkbox);

    if (node.color) {
      const swatch = document.createElement("span");
      swatch.className = "outline-swatch";
      swatch.style.background = rgbToCss(node.color);
      row.appendChild(swatch);
    }
  }

  const label = document.createElement("span");
  label.className = "outline-label";
  label.textContent = node.label;
  if (node.layerId) {
    label.title = "Click to frame in view";
    label.addEventListener("click", () => handlers.onFocus(node.layerId!));
  }
  row.appendChild(label);

  if (node.count !== undefined) {
    const count = document.createElement("span");
    count.className = "outline-count";
    count.textContent = `(${node.count})`;
    row.appendChild(count);
  }

  wrapper.appendChild(row);

  if (node.children) {
    for (const child of node.children) {
      wrapper.appendChild(buildNode(child, depth + 1, handlers));
    }
  }

  return wrapper;
}
