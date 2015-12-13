// Fluid simulation based on Stam's 2003 paper
// Graphics, UI implemented and some important routine changes

// Current issues: bounds
#include <SDL2/SDL.h>
#include <stdio.h>
#include <iostream>
#include <algorithm>

using namespace std;

typedef std::vector<float> vfloat;

// Constants
const int SCREEN_WIDTH = 800;
const int SCREEN_HEIGHT = 800;  // Should match SCREEN_WIDTH
const int N = 20;               // Grid size
const int SIM_LEN = 1000;
const int DELAY_LENGTH = 40;    // ms

const float VISC = .01;
const float dt = 0.1;
const float DIFF = 1;

const bool DISPLAY_CONSOLE = false; // Console or graphics
const bool DRAW_GRID = false; // implement later
const bool DRAW_VEL = true;


// Code begins here
const int nsize = (N+2)*(N+2);

inline int IX(int i, int j){return i + (N+2)*j;}

// Bounds (currently a box with solid walls)
void set_bnd(int N, const int b, vfloat &x, vector<bool> &bound)
{

    for (int i=1; i<=N; i++)
    {
        x[IX(0  ,i)] = b==1 ? -x[IX(1,i)] : 0.5*x[IX(1,i)];
        x[IX(N+1,i)] = b==1 ? -x[IX(N,i)] : 0.5*x[IX(N,i)];
        x[IX(i,  0)] = b==2 ? -x[IX(i,1)] : 0.5*x[IX(i,1)];
        x[IX(i,N+1)] = b==2 ? -x[IX(i,N)] : 0.5*x[IX(i,N)];
    }

    x[IX(0  ,0  )] = 0.5*(x[IX(1,0  )] + x[IX(0  ,1)]);
    x[IX(0  ,N+1)] = 0.5*(x[IX(1,N+1)] + x[IX(0  ,N)]);
    x[IX(N+1,0  )] = 0.5*(x[IX(N,0  )] + x[IX(N+1,1)]);
    x[IX(N+1,N+1)] = 0.5*(x[IX(N,N+1)] + x[IX(N+1,N)]);

    // Boundaries must be 2+ cells thick
    for (int i=1; i<=N; i++)
    {
        for (int j=1; j<=N; j++)
        {
            if (b==1 && bound[IX(i,j)])
            {
                x[IX(i,j)] = (bound[IX(i-1,j)] && bound[IX(i+1,j)]) ? 0 : - x[IX(i-1,j)] - x[IX(i+1,j)];
            }
            if (b==2 && bound[IX(i,j)])
            {
                x[IX(i,j)] = (bound[IX(i,j-1)] && bound[IX(i,j+1)]) ? 0 : - x[IX(i,j-1)] - x[IX(i,j+1)];
            }

        }
    }




}

inline void lin_solve(int N, int b, vfloat &x, const vfloat &x0, float a, float c, vector<bool> &bound)
{
    for (int k=0; k<20; k++)
    {
        for (int i=1; i<=N; i++)
        {
            for (int j=1; j<=N; j++)
            {
                x[IX(i,j)] = (x0[IX(i,j)] + a*(x[IX(i-1,j)]+x[IX(i+1,j)]+x[IX(i,j-1)]+x[IX(i,j+1)])) / c;
            }
        }
        set_bnd (N, b, x, bound);
    }
}

// Add forces
void add_source(vfloat &x, const vfloat &s, float dt)
{
    for (int i=0; i<nsize; i++) x[i] += dt*s[i];
}

// Diffusion with Gauss-Seidel relaxation
void diffuse(int N, int b, vfloat &x, const vfloat &x0, float diff, float dt, vector<bool> &bound)
{
    float a = dt*diff*N*N;
    lin_solve(N, b, x, x0, a, 1+4*a, bound);

}

// Backwards advection
void advect(int N, int b, vfloat &d, const vfloat &d0, const vfloat &u, const vfloat &v, float dt, vector<bool> &bound)
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
    set_bnd(N, b, d, bound);
}

// Force velocity to be mass-conserving (Poisson equation black magic)
void project(int N, vfloat &u, vfloat &v, vfloat &p, vfloat &div, vector<bool> bound)
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
    set_bnd(N, 0, div, bound); set_bnd(N, 0, p, bound);

    lin_solve(N, 0, p, div, 1, 4, bound);

    for (int i=1; i<=N; i++)
    {
        for (int j=1; j<=N; j++)
        {
            u[IX(i,j)] -= 0.5*(p[IX(i+1,j)] - p[IX(i-1,j)])/h;
            v[IX(i,j)] -= 0.5*(p[IX(i,j+1)] - p[IX(i,j-1)])/h;
        }
    }
    set_bnd(N, 1, u, bound); set_bnd(N, 2, v, bound);
}

// Density solver
void dens_step(int N, vfloat &x, vfloat &x0, vfloat &u, vfloat &v, float diff, float dt, vector<bool> &bound)
{
    add_source(x, x0, dt);
    swap(x0, x); diffuse(N, 0, x, x0, diff, dt, bound);
    swap(x0, x); advect(N, 0, x, x0, u, v, dt, bound);
}

// Velocity solver: addition of forces, viscous diffusion, self-advection
void vel_step(int N, vfloat &u, vfloat &v, vfloat &u0, vfloat &v0, float visc, float dt, vector<bool> &bound)
{
    add_source(u, u0, dt); add_source(v, v0, dt);
    swap(u0, u); diffuse(N, 1, u, u0, visc, dt, bound);
    swap(v0, v); diffuse(N, 2, v, v0, visc, dt, bound);
    project(N, u, v, u0, v0, bound);
    swap(u0, u); swap(v0, v);
    advect(N, 1, u, u0, u0, v0, dt, bound); advect(N, 2, v, v0, u0, v0, dt, bound);
    project(N, u, v, u0, v0, bound);
}


void console_draw(vfloat &x)
{
    for (int j=N+1; j>=0; j--)
    {
        for (int i=0; i<=N+1; i++)
            printf("%.3f\t", x[IX(i,j)]);
        cout << '\n';
    }
    cout << '\n';
}


void screen_draw(SDL_Renderer *renderer, vfloat &dens, vfloat &u, vfloat &v, vector<bool> &bound)
{
    const float r_size = (SCREEN_WIDTH / float(N+2));
    for (int i=0; i<=N+1; i++)
    {
        for (int j=0; j<=N+1; j++)
        {
            SDL_Rect r;
            r.w = r_size+1;
            r.h = r_size+1;
            r.x = round(i*r_size);
            r.y = round((N+1-j)*r_size);

            if (bound[IX(i,j)] == 0)
            {
                int color = min(int(dens[IX(i,j)] * 255.0), 255);
                SDL_SetRenderDrawColor(renderer, color, color, color, 0);
            }
            else
            {
                SDL_SetRenderDrawColor(renderer, 0, 100, 100, 0);
            }

            // Render rect
            SDL_RenderFillRect(renderer, &r);

            if (DRAW_VEL)
            {
                SDL_SetRenderDrawColor(renderer, 255, 0, 0, 0);

                int x1 = round((i+0.5)*r_size);
                int y1 = round((N+1-j+0.5)*r_size);
                int x2 = x1 + r_size*u[IX(i,j)];
                int y2 = y1 + r_size*v[IX(i,j)];
                SDL_RenderDrawLine(renderer, x1, y1, x2, y2);
            }
        }
    }

    // Render the rect to the screen
    SDL_RenderPresent(renderer);

    SDL_Delay(DELAY_LENGTH);
}


int main(int, char **)
{
    static vfloat u(nsize, 0), v(nsize, 0), u_prev(nsize, 0), v_prev(nsize, 0); // Horizontal, vertical velocity
    static vfloat dens(nsize, 0), dens_prev(nsize, 0);
    static vector<bool> bounds(nsize, 0);
    //fill_n(dens_prev, nsize, 0.0);

    // SDL initialize
    SDL_Window* window = NULL;

    if ( SDL_Init( SDL_INIT_VIDEO ) < 0 )
        printf( "SDL could not initialize! SDL_Error: %s\n", SDL_GetError() );

    window = SDL_CreateWindow( "SDL Window", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
                              SCREEN_WIDTH, SCREEN_HEIGHT,
                              SDL_WINDOW_SHOWN );
    if( window == NULL )
        printf( "Window could not be created! SDL_Error: %s\n", SDL_GetError() );


    SDL_Renderer* renderer = NULL;
    renderer = SDL_CreateRenderer(window, 0, SDL_RENDERER_ACCELERATED);

    SDL_SetRenderDrawColor(renderer, 255, 0, 255, 255); // Background color, should not see this
    SDL_RenderClear(renderer);

    for (int i=5; i<=7; i++)
    {
        for (int j=0; j<=10; j++)
            bounds[IX(i,j)] = 1;
    }

    for (int t=0; t<SIM_LEN; t++)
    {
        // Get from UI

        for (int j=1; j<=N/4.0; j++)
            u_prev[IX(1,j)] = 50.0;

        /*
        for (int i=N; i>3*N/4.0; i--)
            v_prev[IX(i,1)] = 20.0;

        for (int j=N; j>=3*N/4.0; j--)
            u_prev[IX(N,j)] = -20.0;
        */


        dens_prev[IX(3,3)] = (t<100) ? 100.0 : 0.0;


        vel_step(N, u, v, u_prev, v_prev, VISC, dt, bounds);
        dens_step(N, dens, dens_prev, u, v, DIFF, dt, bounds);

        if (DISPLAY_CONSOLE)
        {
            cout << "dens" << endl;
            console_draw(dens);
            cout << "u, v" << endl;
            console_draw(u);
            console_draw(v);
            cout << "dens_prev" << endl;
            console_draw(dens_prev);
        }

        screen_draw(renderer, dens, u, v, bounds);
    }
    SDL_Quit();
    return 0;
}
