import * as vscode from "vscode";
import { parseMdpa } from "./parser/mdpaParser";

/**
 * A read-only custom editor that renders a `.mdpa` file as a 3D mesh preview
 * with a ModelPart/SubModelPart outline. The document text is parsed on the
 * extension-host side and the resulting model is posted to the webview, which
 * owns all rendering. Edits to the underlying text document re-parse and
 * re-post (debounced).
 */
export class MdpaEditorProvider implements vscode.CustomTextEditorProvider {
  public static readonly viewType = "kratos.mdpaPreview";

  /** The most recently focused preview, so window-level commands can reach it. */
  private activePanel: vscode.WebviewPanel | undefined;

  constructor(private readonly context: vscode.ExtensionContext) {}

  public postToActive(message: unknown): void {
    this.activePanel?.webview.postMessage(message);
  }

  public resolveCustomTextEditor(
    document: vscode.TextDocument,
    webviewPanel: vscode.WebviewPanel,
    _token: vscode.CancellationToken
  ): void {
    const mediaRoot = vscode.Uri.joinPath(this.context.extensionUri, "media");
    webviewPanel.webview.options = {
      enableScripts: true,
      localResourceRoots: [mediaRoot],
    };
    webviewPanel.webview.html = this.getHtml(webviewPanel.webview);

    this.activePanel = webviewPanel;

    const postModel = () => {
      try {
        const model = parseMdpa(document.getText());
        webviewPanel.webview.postMessage({
          type: "model",
          model,
          fileName: document.uri.path.split("/").pop() ?? "",
        });
      } catch (err) {
        webviewPanel.webview.postMessage({
          type: "error",
          message: err instanceof Error ? err.message : String(err),
        });
      }
    };

    // Debounce re-parsing on rapid edits to keep large files responsive.
    let debounce: ReturnType<typeof setTimeout> | undefined;
    const changeSub = vscode.workspace.onDidChangeTextDocument((e) => {
      if (e.document.uri.toString() !== document.uri.toString()) {
        return;
      }
      if (debounce) {
        clearTimeout(debounce);
      }
      debounce = setTimeout(postModel, 250);
    });

    const viewStateSub = webviewPanel.onDidChangeViewState((e) => {
      if (e.webviewPanel.active) {
        this.activePanel = e.webviewPanel;
      } else if (this.activePanel === e.webviewPanel) {
        this.activePanel = undefined;
      }
    });

    const msgSub = webviewPanel.webview.onDidReceiveMessage((msg) => {
      if (msg?.type === "ready") {
        postModel();
      }
    });

    webviewPanel.onDidDispose(() => {
      if (debounce) {
        clearTimeout(debounce);
      }
      changeSub.dispose();
      viewStateSub.dispose();
      msgSub.dispose();
      if (this.activePanel === webviewPanel) {
        this.activePanel = undefined;
      }
    });
  }

  private getHtml(webview: vscode.Webview): string {
    const mediaUri = (file: string) =>
      webview.asWebviewUri(
        vscode.Uri.joinPath(this.context.extensionUri, "media", file)
      );
    const scriptUri = mediaUri("webview.js");
    const styleUri = mediaUri("style.css");
    const nonce = getNonce();
    const csp = [
      `default-src 'none'`,
      `img-src ${webview.cspSource} https: data:`,
      `style-src ${webview.cspSource} 'unsafe-inline'`,
      `script-src 'nonce-${nonce}'`,
      `worker-src blob:`,
      `child-src blob:`,
    ].join("; ");

    return /* html */ `<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8" />
  <meta http-equiv="Content-Security-Policy" content="${csp}" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <link href="${styleUri}" rel="stylesheet" />
  <title>MDPA Preview</title>
</head>
<body>
  <div id="app">
    <aside id="sidebar">
      <div id="stats"></div>
      <div id="outline-header">Layers</div>
      <div id="outline"></div>
    </aside>
    <div id="viewport">
      <div id="toolbar">
        <button data-action="reset" title="Reset camera">Reset</button>
        <button data-action="wireframe" title="Toggle wireframe">Wireframe</button>
        <button data-action="nodeIds" title="Toggle node ids">Node IDs</button>
      </div>
      <div id="render-root"></div>
    </div>
  </div>
  <script nonce="${nonce}" src="${scriptUri}"></script>
</body>
</html>`;
  }
}

function getNonce(): string {
  const chars =
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
  let text = "";
  for (let i = 0; i < 32; i++) {
    text += chars.charAt(Math.floor(Math.random() * chars.length));
  }
  return text;
}
