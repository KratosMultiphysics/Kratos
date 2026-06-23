// Bundles the extension host (CommonJS for the Node VS Code runtime) and the
// webview script (ESM/IIFE for the browser-like webview) in one step.
// Never run `tsc` to emit the shipped code — esbuild owns that; `tsc` is only
// used for type-checking (`npm run typecheck`) and compiling tests.
const esbuild = require("esbuild");
const fs = require("fs");
const path = require("path");

const production = process.argv.includes("--production");
const watch = process.argv.includes("--watch");

// Copy the webview stylesheet (kept as source under webview/) into the shipped
// media/ folder, which is the webview's only allowed local resource root.
const copyStylePlugin = {
  name: "copy-style",
  setup(build) {
    build.onEnd(() => {
      const out = path.join(__dirname, "media");
      fs.mkdirSync(out, { recursive: true });
      fs.copyFileSync(
        path.join(__dirname, "webview", "style.css"),
        path.join(out, "style.css")
      );
    });
  },
};

/** @type {import('esbuild').BuildOptions} */
const extensionConfig = {
  entryPoints: ["src/extension.ts"],
  bundle: true,
  format: "cjs",
  platform: "node",
  target: "node18",
  outfile: "dist/extension.js",
  external: ["vscode"],
  sourcemap: !production,
  minify: production,
  logLevel: "info",
};

/** @type {import('esbuild').BuildOptions} */
const webviewConfig = {
  entryPoints: ["webview/main.ts"],
  bundle: true,
  format: "iife",
  platform: "browser",
  target: "es2021",
  outfile: "media/webview.js",
  sourcemap: !production,
  minify: production,
  logLevel: "info",
  // vtk.js expects a browser global and a NODE_ENV; esbuild's browser platform
  // does not provide these by default.
  define: {
    global: "globalThis",
    "process.env.NODE_ENV": production ? '"production"' : '"development"',
  },
  plugins: [copyStylePlugin],
};

async function main() {
  const contexts = await Promise.all([
    esbuild.context(extensionConfig),
    esbuild.context(webviewConfig),
  ]);

  if (watch) {
    await Promise.all(contexts.map((c) => c.watch()));
    console.log("[esbuild] watching...");
  } else {
    await Promise.all(contexts.map((c) => c.rebuild()));
    await Promise.all(contexts.map((c) => c.dispose()));
  }
}

main().catch((e) => {
  console.error(e);
  process.exit(1);
});
