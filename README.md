#COMS30115 coursework by @Michael-JB and @SSutherlandDeeBristol

## Extensions

### Raytracer
* Phong Shading
* Path Tracer
    * Imperfect reflections
* Photon Mapper
    * Area lighting
    * Spot lighting
    * Reflections & mirrors
    * Refraction & glass
    * Soft shadows
* Spheres
* Object loading from .obj files
* Full camera movement with LookAt

### Rasteriser
* Full homogenous clipping
* Fan triangulation algorithm
* Optimal triangulation algorithm
* Object loading from .obj files
* Full camera movement with LookAt

### Raycaster
As well as extending the raytracer and rasteriser, we created a raycasting renderer. The ray- casting renderer is inspired by Wolfenstien 3D, building up a 3D scene from the intersection and distance data of rays cast in 2D space. The key features of our raycasting renderer are listed below.

* Minimap
* Fisheye Correction
* Distance Shading
* Floor and Sky
* Full camera movement with LookAt

## Gallery

### Photon Mapper

Balls Render               |  Balls Photon Map
:-------------------------:|:-------------------------:
![](https://github.com/SSutherlandDeeBristol/computer-graphics-cw/blob/master/submission/images/raytracer/photon%20mapper/photon3.png)  |  ![](https://github.com/SSutherlandDeeBristol/computer-graphics-cw/blob/master/submission/images/raytracer/photon%20mapper/photonmap3.png)

Bunny Render               |  Bunny Photon Map
:-------------------------:|:-------------------------:
![](https://github.com/SSutherlandDeeBristol/computer-graphics-cw/blob/master/submission/images/raytracer/photon%20mapper/photon2.png)  |  ![](https://github.com/SSutherlandDeeBristol/computer-graphics-cw/blob/master/submission/images/raytracer/photon%20mapper/photonmap2.png)

### Path Tracer

Tall Mirror Block Render |
:-----------------------:|
![](https://github.com/SSutherlandDeeBristol/computer-graphics-cw/blob/master/submission/images/raytracer/pathtracer/pathtracer1.png)

### Phong Shading

Cornell Box Render |
:-----------------------:|
![](https://github.com/SSutherlandDeeBristol/computer-graphics-cw/blob/master/submission/images/raytracer/phong%20shading/phong1.png)

### Rasteriser

A demonstration video can be viewed [here](https://www.youtube.com/watch?v=RgAZK1vxCeg).

Clipping to two planes               |  Optimal triangulation
:-------------------------:|:-------------------------:
![](https://github.com/SSutherlandDeeBristol/computer-graphics-cw/blob/master/submission/images/rasteriser/rasteriser1.png)  |  ![](https://github.com/SSutherlandDeeBristol/computer-graphics-cw/blob/master/submission/images/rasteriser/rasteriser2.png)

### Raycaster

A demonstration video can be viewed [here](https://www.youtube.com/watch?v=gmD1RMafxK0).

Snapshot 1 | Snapshot 2
:-------------------------:|:-------------------------:
![](https://github.com/SSutherlandDeeBristol/computer-graphics-cw/blob/master/submission/images/raycaster/raycaster1.png)  |  ![](https://github.com/SSutherlandDeeBristol/computer-graphics-cw/blob/master/submission/images/raycaster/raycaster2.png)

## How To Run Code

Before running anything, navigate to the top level of the submission directory

```bash
cd submission
```

### Photon Mapper

```bash
cd raytracer
make -B && ./Build/skeleton n r
```

where:
* **n** integer number of photons to emit
* **r** integer number of photons in radiance estimate

To start try:
* **n** = 20000
* **r** = 200

Increasing these values will give a nicer render but will significantly increase render time.
