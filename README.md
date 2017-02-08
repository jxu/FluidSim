FluidSim
========
Grid-based fluid simulation from approximations of Navier-Stokes equations, based on a paper by Stam. 
See source code for command line options. Use `--demo` or `--demo2` for demos. 

Requirements
------------
- Running on Windows: SDL2.dll from SDL2 in the same directory as executable file
- Running on Linux: libsdl2-2.0-0 package

Compiling on Linux
------------------
Install libsdl2-dev.

Compile with `g++ -O3 -Wall -o gui gui.cpp -lSDL2`
`
