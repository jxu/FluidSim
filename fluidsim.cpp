// Fluid simulation based on Stam's 2003 paper
// Graphics, UI implemented and some important routine changes

// Current issues: runaway values, bounds
#include <SDL2/SDL.h>
#include <stdio.h>
#include <iostream>
#include <algorithm>


using namespace std;

// Constants
const int SCREEN_WIDTH = 600;
const int SCREEN_HEIGHT = 600;  // Should match SCREEN_WIDTH
const int N = 20;               // Grid size
const int SIM_LEN = 1000;
const int DELAY_LENGTH = 40;    // ms

const float VISC = 0.01;
const float dt = 0.1;
const float DIFF = 0.1;

const bool DISPLAY_CONSOLE = false; // Console or graphics
const bool DRAW_GRID = false; // implement later


// Code begins here
const int nsize = (N+2)*(N+2);

inline int IX(int i, int j){return i + (N+2)*j;}

// Bounds (currently a box with solid walls)
void set_bnd(int N, int b, vector<float> &x)
{
    /*
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
    */

}

inline void lin_solve(int N, int b, vector<float> &x, const vector<float> &x0, float a, float c)
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
        set_bnd ( N, b, x );
    }
}

// Add forces
void add_source(vector<float> &x, const vector<float> &s, float dt)
{
    for (int i=0; i<nsize; i++) x[i] += dt*s[i];
}

// Diffusion with Gauss-Seidel relaxation
void diffuse(int N, int b, vector<float> &x, const vector<float> &x0, float diff, float dt)
{
    float a = dt*diff*N*N;
    lin_solve(N, b, x, x0, a, 1+4*a);

}

// Backwards advection
void advect(int N, int b, vector<float> &d, const vector<float> &d0, const vector<float> &u, const vector<float> &v, float dt)
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
    set_bnd(N, b, d);
}

// Force velocity to be mass-conserving (Poisson equation black magic)
void project(int N, vector<float> &u, vector<float> &v, vector<float> &p, vector<float> &div)
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

    lin_solve(N, 0, p, div, 1, 4);

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
void dens_step(int N, vector<float> &x, vector<float> &x0, vector<float> &u, vector<float> &v, float diff, float dt)
{
    add_source(x, x0, dt);
    swap(x0, x); diffuse(N, 0, x, x0, diff, dt);
    swap(x0, x); advect(N, 0, x, x0, u, v, dt);
}

// Velocity solver: addition of forces, viscous diffusion, self-advection
void vel_step(int N, vector<float> &u, vector<float> &v, vector<float> &u0, vector<float> &v0, float visc, float dt)
{
    add_source(u, u0, dt); add_source(v, v0, dt);
    swap(u0, u); diffuse(N, 1, u, u0, visc, dt);
    swap(v0, v); diffuse(N, 2, v, v0, visc, dt);
    project(N, u, v, u0, v0);
    swap(u0, u); swap(v0, v);
    advect(N, 1, u, u0, u0, v0, dt); advect(N, 2, v, v0, u0, v0, dt);
    project(N, u, v, u0, v0);
}


void console_draw(vector<float> &x)
{
    for (int j=N+1; j>=0; j--)
    {
        for (int i=0; i<=N+1; i++)
            printf("%.3f\t", x[IX(i,j)]);
        cout << '\n';
    }
    cout << '\n';
}


void screen_draw(SDL_Renderer *renderer, vector<float> &dens, vector<float> &u, vector<float> &v)
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

            int color = min(int(dens[IX(i,j)] * 255.0), 255);

            SDL_SetRenderDrawColor(renderer, color, color, color, 0);
            // Render rect
            SDL_RenderFillRect(renderer, &r);

            SDL_SetRenderDrawColor(renderer, 255, 0, 0, 0);

            int x1 = round((i+0.5)*r_size);
            int y1 = round((N+1-j+0.5)*r_size);
            int x2 = x1 + r_size*u[IX(i,j)];
            int y2 = y1 + r_size*v[IX(i,j)];
            SDL_RenderDrawLine(renderer, x1, y1, x2, y2);

        }
    }

    // Render the rect to the screen
    SDL_RenderPresent(renderer);

    SDL_Delay(DELAY_LENGTH);
}


int main(int, char **)
{
    static vector<float> u(nsize, 0), v(nsize, 0), u_prev(nsize, 0), v_prev(nsize, 0); // Horizontal, vertical velocity
    static vector<float> dens(nsize, 0), dens_prev(nsize, 0);

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


    for (int t=0; t<SIM_LEN; t++)
    {
        // Get from UI
        for (int j=1; j<=N/4.0; j++)
            u_prev[IX(1,j)] = 10.0;

        for (int i=N; i>3*N/4.0; i--)
            v_prev[IX(i,1)] = -10.0;



        dens_prev[IX(3,3)] = (t<100) ? 100.0 : 0.0;


        vel_step(N, u, v, u_prev, v_prev, VISC, dt);
        dens_step(N, dens, dens_prev, u, v, DIFF, dt);

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

        screen_draw(renderer, dens, u, v);
    }



    SDL_Quit();
    return 0;
}
