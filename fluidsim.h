// Fluid simulation based on Stam's 2003 paper
// Graphics, UI implemented and some important routine changes

// To do: fix density boundary issues
#pragma once

#include <algorithm>

typedef std::vector<float> vfloat;
typedef std::vector<bool> vbool;

extern const int N, nsize;

inline int IX(int i, int j){return i + (N+2)*j;}

// bound (currently a box with solid walls)
void set_bnd(const int b, vfloat &x, vbool &bound)
{
    for (int i=1; i<=N; i++)
    {
        x[IX(0  ,i)] = b==1 ? -x[IX(1,i)] : x[IX(1,i)];
        x[IX(N+1,i)] = b==1 ? -x[IX(N,i)] : x[IX(N,i)];
        x[IX(i,  0)] = b==2 ? -x[IX(i,1)] : x[IX(i,1)];
        x[IX(i,N+1)] = b==2 ? -x[IX(i,N)] : x[IX(i,N)];
    }

    x[IX(0  ,0  )] = 0.5*(x[IX(1,0  )] + x[IX(0  ,1)]);
    x[IX(0  ,N+1)] = 0.5*(x[IX(1,N+1)] + x[IX(0  ,N)]);
    x[IX(N+1,0  )] = 0.5*(x[IX(N,0  )] + x[IX(N+1,1)]);
    x[IX(N+1,N+1)] = 0.5*(x[IX(N,N+1)] + x[IX(N+1,N)]);

    // Boundaries should be 2+ cells thick
    for (int i=1; i<=N; i++)
    {
        for (int j=1; j<=N; j++)
        {
            if (bound[IX(i,j)])
            {
                if (b==1)
                    x[IX(i,j)] = (bound[IX(i-1,j)] && bound[IX(i+1,j)]) ? 0 : - x[IX(i-1,j)] - x[IX(i+1,j)];
                else if (b==2)
                    x[IX(i,j)] = (bound[IX(i,j-1)] && bound[IX(i,j+1)]) ? 0 : - x[IX(i,j-1)] - x[IX(i,j+1)];

                else if (b==0)
                {
                    // Distribute density from bound to surrounding cells
                    int nearby_count = !bound[IX(i+1,j)] + !bound[IX(i-1,j)] + !bound[IX(i,j+1)] + !bound[IX(i,j-1)];
                    if (!nearby_count) x[IX(i,j)] = 0;
                    else
                        x[IX(i,j)] = ((bound[IX(i+1,j)] ? 0 : x[IX(i+1,j)]) +
                                      (bound[IX(i-1,j)] ? 0 : x[IX(i-1,j)]) +
                                      (bound[IX(i,j+1)] ? 0 : x[IX(i,j+1)]) +
                                      (bound[IX(i,j-1)] ? 0 : x[IX(i,j-1)])) / nearby_count;
                }
            }
        }
    }
}

inline void lin_solve(int b, vfloat &x, const vfloat &x0, float a, float c, vbool &bound)
{
    for (int k=0; k<20; k++)
    {
        for (int i=1; i<=N; i++)
        {
            for (int j=1; j<=N; j++)
                x[IX(i,j)] = (x0[IX(i,j)] +
                              a*(x[IX(i-1,j)]+x[IX(i+1,j)]+x[IX(i,j-1)]+x[IX(i,j+1)])) / c;
        }
        set_bnd (b, x, bound);
    }
}

// Add forces
void add_source(vfloat &x, const vfloat &s, float dt)
{
    for (int i=0; i<nsize; i++) x[i] += dt*s[i];
}

// Diffusion with Gauss-Seidel relaxation
void diffuse(int b, vfloat &x, const vfloat &x0, float diff, float dt, vbool &bound)
{
    float a = dt*diff*N*N;
    lin_solve(b, x, x0, a, 1+4*a+dt, bound); // Amazing fix due to Iwillnotexist Idonotexist
}

// Backwards advection
void advect(int b, vfloat &d, const vfloat &d0, const vfloat &u, const vfloat &v, float dt, vbool &bound)
{
    float dt0 = dt*N;
    for (int i=1; i<=N; i++)
    {
        for (int j=1; j<=N; j++)
        {
            float x = i - dt0*u[IX(i,j)];
            float y = j - dt0*v[IX(i,j)];
            if (x<0.5) x=0.5; if (x>N+0.5) x=N+0.5;
            int i0=(int)x; int i1=i0+1;
            if (y<0.5) y=0.5; if (y>N+0.5) y=N+0.5;
            int j0=(int)y; int j1=j0+1;

            float s1 = x-i0; float s0 = 1-s1; float t1 = y-j0; float t0 = 1-t1;
            d[IX(i,j)] = s0*(t0*d0[IX(i0,j0)] + t1*d0[IX(i0,j1)]) +
                         s1*(t0*d0[IX(i1,j0)] + t1*d0[IX(i1,j1)]);
        }
    }
    set_bnd(b, d, bound);
}

// Force velocity to be mass-conserving (Poisson equation black magic)
void project(vfloat &u, vfloat &v, vfloat &p, vfloat &div, vbool bound)
{
    float h = 1.0/N;
    for (int i=1; i<=N; i++)
    {
        for (int j=1; j<=N; j++)
        {
            div[IX(i,j)] = -0.5*h*(u[IX(i+1,j)] - u[IX(i-1,j)] +
                                   v[IX(i,j+1)] - v[IX(i,j-1)]);
            p[IX(i,j)] = 0;
        }
    }
    set_bnd(0, div, bound); set_bnd(0, p, bound);

    lin_solve(0, p, div, 1, 4, bound);

    for (int i=1; i<=N; i++)
    {
        for (int j=1; j<=N; j++)
        {
            u[IX(i,j)] -= 0.5*(p[IX(i+1,j)] - p[IX(i-1,j)])/h;
            v[IX(i,j)] -= 0.5*(p[IX(i,j+1)] - p[IX(i,j-1)])/h;
        }
    }
    set_bnd(1, u, bound); set_bnd(2, v, bound);
}

// Density solver
void dens_step(vfloat &x, vfloat &x0, vfloat &u, vfloat &v, float diff, float dt, vbool &bound)
{
    add_source(x, x0, dt);
    swap(x0, x); diffuse(0, x, x0, diff, dt, bound);
    swap(x0, x); advect(0, x, x0, u, v, dt, bound);
}

// Velocity solver: addition of forces, viscous diffusion, self-advection
void vel_step(vfloat &u, vfloat &v, vfloat &u0, vfloat &v0, float visc, float dt, vbool &bound)
{
    add_source(u, u0, dt); add_source(v, v0, dt);
    swap(u0, u); diffuse(1, u, u0, visc, dt, bound);
    swap(v0, v); diffuse(2, v, v0, visc, dt, bound);
    project(u, v, u0, v0, bound);
    swap(u0, u); swap(v0, v);
    advect(1, u, u0, u0, v0, dt, bound); advect(2, v, v0, u0, v0, dt, bound);
    project(u, v, u0, v0, bound);
}
