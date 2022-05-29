# Vega FEM

From the [website](https://viterbi-web.usc.edu/~jbarbic/vega):

> Vega is a computationally efficient and stable C/C++ physics library for three-dimensional deformable object simulation. It is designed to model large deformations, including geometric and material nonlinearities, and can also efficiently simulate linear systems. Vega contains about 145,000 lines of code, and is open-source and free. It is released under the [3-clause BSD license](http://en.wikipedia.org/wiki/BSD_licenses), which means that it can be used freely both in academic research and in commercial applications.

The fork in this repository consists mostly of providing CMake files for simpler and cross-platform build generation, fixing bugs and C-ism. Currently it is compiling and generating the binaries for Linux and Windows as proven by the [CI](https://github.com/FabienPean/VegaFEM/actions).

## Motivation

The Vega FEM library is relatively less known due to its official absence from principal git websites (GitHub, GitLab, BitBucket, etc). Yet, you can see it perform commercially in [Ziva VFX](https://docs.zivadynamics.com/vfx/third_party_licenses.html#vega-fem) from [Ziva Dynamics](https://zivadynamics.com/) and in [research](https://doi.org/10.1145/3197517.3201327).

The library provides many necessary tools for the development of physics-based animation with a permissive license. In particular, it contains a distance grid (level set) and tetrahedral mesh generators with little to no dependencies. Common alternatives are [TetGen](https://wias-berlin.de/software/index.jsp?id=TetGen), which is under a [copyleft license](https://tldrlegal.com/license/gnu-affero-general-public-license-v3-(agpl-3.0)), or [fTetWild](https://github.com/wildmeshing/fTetWild), which is a patchwork of many libraries difficult to build all together.

Former forks available on GitHub are relatively old and with no modification. This fork contains the most recent version accompanied with fixes and a build generator based on CMake.

## Quickstart

If you just need the utilities, the simplest is to download the binaries in the release section. Otherwise you can build the repository on your system assuming the following:

* The OS is Linux with Clang or Windows with Clang or MSVC
* [Intel oneAPI MKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html) is installed on the system, as well as the following dependencies: `intel-mkl opengl glui glew cgal openblas eigen3`  which can be obtained via [vcpkg](https://vcpkg.io/en/index.html)
* CMake 3.21

## License

The library itself is released under the BSD 3-clause. However, due to direct dependencies on the  [Polygon Mesh Processing](https://doc.cgal.org/latest/Polygon_mesh_processing/group__PkgPolygonMeshProcessingRef.html)  package of [CGAL](https://www.cgal.org/) of `libraries/libiglInterface`, the distributed libraries or binaries relying directly on it, or on `libraries/virtualTets` and `libraries/immersionMesher`, must be released under the GPL. 



