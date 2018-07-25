"""Stitch borders of a quad."""

import brainmesh as bm
from pathlib import Path

outdir = Path("results")
outdir.mkdir(exist_ok=True)

quad = bm.BrainSurface("data/full_border_quads.off")
print("Before stitching")
print(f"number of faces: {quad.num_faces()}")
print(f"number of edges: {quad.num_edges()}")
print(f"number of vertices: {quad.num_vertices()}")

quad.stitch_borders()

print("After stitching")
print(f"number of faces: {quad.num_faces()}")
print(f"number of edges: {quad.num_edges()}")
print(f"number of vertices: {quad.num_vertices()}")

quad.save(str(outdir/"stitched_quads.off"))
