# computational geometrical kernels

A c++ implementation of the following geometric kernels

- V-polyhedron to H-polyhedron
- H-polyhedron to V-polyhedron

- connectivity: V-polyhedron <--> H-polyhedron
  - which face indices constitute the vertex (H -> V)
  - which vertex indices constitute the face (V -> H)

- compute Chebyshev ball (inscribed ball) from H-polyhedron, possibly infeasible (non-intersection)

> Thanks to: [cddlib](https://github.com/cddlib/cddlib), [eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page)

