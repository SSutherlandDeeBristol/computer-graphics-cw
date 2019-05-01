28371 & 34517

Rasteriser:
- Full homogenous clipping
- Optimal triangulation algorithm
- Object loading from .obj files

Raytracer:
- Path Tracing
- Path Tracing reflection
- Photon Mapper
- Area Lights / Spotlight
- Reflection / Mirrors
- Refraction / Glass
- Soft Shadows
- Phong Shading
- Object loading from .obj files

HOW TO RUN CODE:
(assuming in top level directory of submission)

PATHTRACER:

cd pathtracer
make -B && ./Build/skeleton b n

where b = number of bounces
      n = number of samples at each bounce
      both are integers

PHOTON MAPPER:

cd raytracer
make -B && ./Build/skeleton n r

where n = number of photons to emit
      r = number of photons in radiance estimate

PHONG SHADING:

cd raytracer
make -B && ./Build/skeleton

RASTERISER:

cd rasteriser
make -B && ./Build/skeleton
