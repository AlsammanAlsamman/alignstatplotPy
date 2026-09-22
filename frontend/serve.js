#!/usr/bin/env node
// Zero-dependency static file server for the alignstatplot web GUI frontend.
// The FastAPI backend (`alignstatplot web`, from the alignstatplot-py Python
// package) already serves these same files directly; this script exists so
// the frontend can also run standalone (e.g. pointed at a backend running
// on a different host) via `npx alignstatplot-web-gui`.
import http from "node:http";
import { readFile } from "node:fs/promises";
import path from "node:path";
import { fileURLToPath } from "node:url";

const ROOT = path.dirname(fileURLToPath(import.meta.url));
const PORT = process.env.PORT ? Number(process.env.PORT) : 5173;

const MIME = { ".html": "text/html", ".js": "text/javascript", ".css": "text/css" };

const server = http.createServer(async (req, res) => {
  const reqPath = req.url === "/" ? "/index.html" : req.url;
  const filePath = path.join(ROOT, path.normalize(reqPath).replace(/^(\.\.[/\\])+/, ""));
  try {
    const data = await readFile(filePath);
    res.writeHead(200, { "Content-Type": MIME[path.extname(filePath)] || "application/octet-stream" });
    res.end(data);
  } catch {
    res.writeHead(404);
    res.end("Not found");
  }
});

server.listen(PORT, () => {
  console.log(`alignstatplot web GUI frontend running at http://localhost:${PORT}`);
  console.log("Make sure the Python backend is running too: `alignstatplot web`");
});
