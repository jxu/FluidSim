// Fluid simulation based on Stam's 2003 paper
// Graphics, UI implemented and some important routine changes


#include <SDL2/SDL.h>
#include <stdio.h>
#include <iostream>

#include <algorithm> // swap


#define IX(i,j) ((i)+(N+2)*(j))

using namespace std;

// Constants
const int SCREEN_WIDTH = 640;
const int SCREEN_HEIGHT = 480;
const int N = 6; // Grid size

const float VISC = 0.5;
const float dt = 0.5;
const float DIFF = .1;

const int DISPLAY_MODE = 0; // Console or graphics
const bool CLEAR_CONSOLE = false;


// Code begins here
const int nsize = (N+2)*(N+2);

// Bounds (currently a box with solid walls)
void set_bnd(int N, int b, float *x)
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
}

// Add forces
void add_source(int N, float *x, float *s, float dt)
{
    for (int i=0; i<nsize; i++) x[i] += dt*s[i];
}

// Diffusion with Gauss-Seidel relaxation
void diffuse(int N, int b, float *x, float *x0, float diff, float dt)
{
    float a = dt*diff*N*N;
    float y[nsize] = {0};

    for (int k=0; k<20; k++)
    {
        for (int i=1; i<=N; i++)
        {
            for (int j=1; j<=N; j++)
            {
                y[IX(i,j)] = (x0[IX(i,j)] + a*(x[IX(i-1,j)] + x[IX(i+1,j)] +
                                               x[IX(i,j-1)] + x[IX(i,j+1)])) / (1+4*a);
            }
        }
        //std::copy(y, y+nsize, x);
        for (int i=0; i<nsize; i++)
            x[i] = y[i];

        set_bnd(N, b, x);
    }

}

// Backwards advection
void advect(int N, int b, float *d, float *d0, float *u, float *v, float dt)
{
    int i0, j0, i1, j1;
    float x, y, s0, t0, s1, t1, dt0;

    dt0 = dt*N;
    for (int i=1; i<=N; i++)
    {
        for (int j=1; j<=N; j++)
        {
            x = i - dt0*u[IX(i,j)];
            y = j - dt0*v[IX(i,j)];
            if (x<0.5) x=0.5; if (x>N+0.5) x=N+0.5;
            i0=(int)x; i1=i0+1;
            if (y<0.5) y=0.5; if (y>N+0.5) y=N+0.5;
            j0=(int)y; j1=j0+1;

            s1 = x-i0; s0 = 1-s1; t1 = y-j0; t0 = 1-t1;
            d[IX(i,j)] = s0*(t0*d0[IX(i0,j0)] + t1*d0[IX(i0,j1)]) +
                         s1*(t0*d0[IX(i1,j0)] + t1*d0[IX(i1,j1)]);
        }
    }
    set_bnd(N, b, d);
}

// Force velocity to be mass-conserving (Poisson equation black magic)
void project(int N, float *u, float *v, float *p, float *div)
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
    set_bnd(N, 0, div); set_bnd(N, 0, p);

    for (int k=0; k<20; k++)
    {
        for (int i=1; i<=N; i++)
        {
            for (int j=1; j<=N; j++)
            {
                p[IX(i,j)] = (div[IX(i,j)] + p[IX(i-1,j)] + p[IX(i+1,j)] +
                                             p[IX(i,j-1)] + p[IX(i,j+1)])/4;
            }
        }
        set_bnd(N, 0, p);
    }

    for (int i=1; i<=N; i++)
    {
        for (int j=1; j<=N; j++)
        {
            u[IX(i,j)] -= 0.5*(p[IX(i+1,j)] - p[IX(i-1,j)])/h;
            v[IX(i,j)] -= 0.5*(p[IX(i,j+1)] - p[IX(i,j-1)])/h;
        }
    }
    set_bnd(N, 1, u); set_bnd(N, 2, v);
}

// Density solver
void dens_step(int N, float *x, float *x0, float *u, float *v, float diff, float dt)
{
    add_source(N, x, x0, dt);
    swap(x0, x); diffuse(N, 0, x, x0, diff, dt);
    swap(x0, x); advect(N, 0, x, x0, u, v, dt);
}

// Velocity solver: addition of forces, viscous diffusion, self-advection
void vel_step(int N, float *u, float *v, float *u0, float *v0, float visc, float dt)
{
    add_source(N, u, u0, dt); add_source(N, v, v0, dt);
    swap(u0, u); diffuse(N, 1, u, u0, visc, dt);
    swap(v0, v); diffuse(N, 2, v, v0, visc, dt);
    project(N, u, v, u0, v0);
    swap(u0, u); swap(v0, v);
    advect(N, 1, u, u0, u0, v0, dt); advect(N, 2, v, v0, u0, v0, dt);
    project(N, u, v, u0, v0);
}


void console_draw(int N, float *x)
{
    for (int j=N+1; j>=0; j--)
    {
        for (int i=0; i<=N+1; i++)
            printf("%.3f\t", x[IX(i,j)]);
        cout << '\n';
    }
    cout << '\n';
}


void screen_test()
{
    SDL_Window* window = NULL;

    if ( SDL_Init( SDL_INIT_VIDEO ) < 0 )
    {
        printf( "SDL could not initialize! SDL_Error: %s\n", SDL_GetError() );
    }
    else
    {
        //Create window
        window = SDL_CreateWindow( "SDL Tutorial", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
                                  SCREEN_WIDTH, SCREEN_HEIGHT,
                                  SDL_WINDOW_SHOWN );
        if( window == NULL )
        {
            printf( "Window could not be created! SDL_Error: %s\n", SDL_GetError() );
        }
        else
        {
            SDL_Renderer* renderer = NULL;
            renderer = SDL_CreateRenderer(window, 0, SDL_RENDERER_ACCELERATED);

            SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);

            SDL_RenderClear(renderer);

            SDL_Rect r;
            r.x = 50;
            r.y = 50;
            r.w = 50;
            r.h = 50;


            SDL_SetRenderDrawColor(renderer, 0, 0, 255, 255);

            // Render rect
            SDL_RenderFillRect(renderer, &r);

            // Render the rect to the screen
            SDL_RenderPresent(renderer);

            // Wait for 5 sec
            SDL_Delay(5000);

        }
    }
    SDL_DestroyWindow(window);
    SDL_Quit();
}


int main(int argc, char* args[])
{
    float u[nsize] = {0}, v[nsize] = {0}, u_prev[nsize] = {0}, v_prev[nsize] = {0}; // Horizontal, vertical velocity
    float dens[nsize], dens_prev[nsize] = {0};
    float source[nsize] = {0};

    fill_n(dens, nsize, 0.0);
    //fill_n(dens_prev, nsize, 0.0);


    const int SIM_LEN = 10;

    for (int t=0; t<SIM_LEN; t++)
    {
        // Get from UI
        source[IX(1,1)] = 0.3;
        add_source(N, dens, source, dt);
        vel_step(N, u, v, u_prev, v_prev, VISC, dt);
        dens_step(N, dens, dens_prev, u, v, DIFF, dt);
        cout << "dens" << endl;
        console_draw(N, dens);
        //cout << "u, v" << endl;
        //console_draw(N, u_prev);
        //console_draw(N, v_prev);

    }


    return EXIT_SUCCESS;
}
