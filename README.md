# FinalProject

CIS 700 - Procedural Graphics - Design Document - 4k Fractal Demoscene
Tabatha Hickman and Joseph Klinger

# Motivation
This project is motivated primarily by scenes containing fractal geometry created by other demo groups as well as by the various articles on IQ’s website and a general interest in using math to model complex and visually compelling geometry.

# Goals
The primary goals (generally in decreasing order of priority):
Create a scene with interesting fractal geometry
Implement dynamic camera movement to view our scene
Implement some interesting materials for the fractals in glsl
Make the project fit into a 4k executable

# Inspiration/Reference:
https://www.youtube.com/watch?v=LbWS79Cxykg
https://www.youtube.com/watch?v=TTpbP5BVtiA
https://www.youtube.com/watch?v=wiaSPOm_G9k

# Techniques:
3D Mandelbulb: http://mandelubber.blogspot.com/p/tutorials.html
http://www.iquilezles.org/www/articles/mandelbulb/mandelbulb.htm (as well as other fractal-related articles by IQ)
http://www.iquilezles.org/www/articles/ssao/ssao.htm ? screen space ao by IQ
http://www.iquilezles.org/www/articles/terrainmarching/terrainmarching.htm ? raymarching by IQ, will probably help with optimizations once Joe has finished the core engine (see timeline)
Various world/post processing effects discussed in class. This will be more solidified once Tabatha has finished designing the scene (see timeline)
Use some sort of program to compress the code into 4k (Crinkler?). This will be done after the project is mostly done. Will probably have to ask Rachel about this.



# Timeline:
By 4/17: Joe has raymarching engine done and preferably optimized. Tabatha plans out the actual scene and how to move the scene/camera for animation/viewing. Both learn about how to actually generate fractals/3d mandelbrot set/ julia sets in addition to the above. We would like to be viewing some sort of 3d mandelbulbs w/ Joe’s raymarcher by this date.

By 4/24: (Flexible week) Joe has implemented various shader effects to the scene (materials/shading, ao, some kind of background?). Tabatha might have helped with shaders and/or cleaned up details on the actual scene. We need to have started investigating methods for compressing the code into 4k by today (ask Rachel?)

By 5/1: Essentially finished. Depending on the state of the project, we will have spent the week adding effects, writing shaders, interactivity or otherwise. We may end up needing to spend time condensing code if our project doesn’t compress to 4k. Music will need to be added as well, which will like be some royalty-free music we find somewhere, we will not be composing our own.