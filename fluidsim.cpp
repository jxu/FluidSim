// Fluid simulation based on Stam's 2003 paper
// Graphics, UI implemented and some important routine changes

// To do: fix density boundary issues
// Fix framerate (see DELAY_LENGTH below)
#include <SDL2/SDL.h>
#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <chrono>
#include <thread>

typedef std::vector<float> vfloat;

// Constants
const int SCREEN_WIDTH = 800;
const int SCREEN_HEIGHT = 800;  // Should match SCREEN_WIDTH
const int N = 50;               // Grid size
const int SIM_LEN = 3000;       // Based on actual framerate

// Locks framerate at ~64, see stackoverflow.com/q/23258650/3163618
const std::chrono::milliseconds DELAY_LENGTH(5);

const float VISC = 0.01;
const float dt = 0.005;
const float DIFF = 0.01;

const bool DISPLAY_CONSOLE = false; // Console or graphics
const bool DRAW_GRID = false; // implement later
const bool DRAW_VEL = true;

const float MOUSE_DENS = 100.0;


// Code begins here
const int nsize = (N+2)*(N+2);

inline int IX(int i, int j){return i + (N+2)*j;}

// Bounds (currently a box with solid walls)
void set_bnd(const int b, vfloat &x, std::vector<bool> &bound)
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
                    float spread = x[IX(i,j)] / nearby_count;
                    if (!bound[IX(i+1,j)]) x[IX(i+1,j)] += spread;
                    if (!bound[IX(i-1,j)]) x[IX(i-1,j)] += spread;
                    if (!bound[IX(i,j+1)]) x[IX(i,j+1)] += spread;
                    if (!bound[IX(i,j-1)]) x[IX(i,j-1)] += spread;

                    x[IX(i,j)] = 0;
                }
            }
        }
    }
}

inline void lin_solve(int b, vfloat &x, const vfloat &x0, float a, float c, std::vector<bool> &bound)
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
void diffuse(int b, vfloat &x, const vfloat &x0, float diff, float dt, std::vector<bool> &bound)
{
    float a = dt*diff*N*N;
    lin_solve(b, x, x0, a, 1+4*a+dt, bound); // Amazing fix due to Iwillnotexist Idonotexist

}

// Backwards advection
void advect(int b, vfloat &d, const vfloat &d0, const vfloat &u, const vfloat &v, float dt, std::vector<bool> &bound)
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
void project(vfloat &u, vfloat &v, vfloat &p, vfloat &div, std::vector<bool> bound)
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
void dens_step(vfloat &x, vfloat &x0, vfloat &u, vfloat &v, float diff, float dt, std::vector<bool> &bound)
{
    add_source(x, x0, dt);
    swap(x0, x); diffuse(0, x, x0, diff, dt, bound);
    swap(x0, x); advect(0, x, x0, u, v, dt, bound);
}

// Velocity solver: addition of forces, viscous diffusion, self-advection
void vel_step(vfloat &u, vfloat &v, vfloat &u0, vfloat &v0, float visc, float dt, std::vector<bool> &bound)
{
    add_source(u, u0, dt); add_source(v, v0, dt);
    swap(u0, u); diffuse(1, u, u0, visc, dt, bound);
    swap(v0, v); diffuse(2, v, v0, visc, dt, bound);
    project(u, v, u0, v0, bound);
    swap(u0, u); swap(v0, v);
    advect(1, u, u0, u0, v0, dt, bound); advect(2, v, v0, u0, v0, dt, bound);
    project(u, v, u0, v0, bound);
}


void console_write(vfloat &x)
{
    for (int j=N+1; j>=0; j--)
    {
        for (int i=0; i<=N+1; i++)
            printf("%.3f\t", x[IX(i,j)]);
        std::cout << '\n';
    }
    std::cout << '\n';
}


void screen_draw(SDL_Renderer *renderer, vfloat &dens, vfloat &u, vfloat &v, std::vector<bool> &bound)
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
                //if (dens[IX(i,j)] < 2.0)
                {
                    int color = std::min(int(dens[IX(i,j)] * 255), 255);
                    SDL_SetRenderDrawColor(renderer, color, color, color, 0);
                }
                //else // Negative density (error)
                //    SDL_SetRenderDrawColor(renderer, 255, 200, 255, 0);


            }
            else // Object boundary
                SDL_SetRenderDrawColor(renderer, 0, 100, 100, 0);

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
}

// Add density (or velocity) based on user input
void process_input(vfloat &dens_prev, vfloat &dens)
{
    int x, y;
    int *ptr_x = &x, *ptr_y = &y;

    float r_size = (SCREEN_WIDTH / float(N+2));

    SDL_PumpEvents();
    unsigned int mouse_state = SDL_GetMouseState(ptr_x, ptr_y);

    if (mouse_state & (SDL_BUTTON(SDL_BUTTON_LEFT) | SDL_BUTTON(SDL_BUTTON_RIGHT)))
    {
        int grid_x = round(x/r_size);
        int grid_y = N+2 - round(y/r_size);
        if (mouse_state & SDL_BUTTON(SDL_BUTTON_LEFT))
        {
            std::cout << "Left ";
            dens_prev[IX(grid_x,grid_y)] += MOUSE_DENS;
        }

        if (mouse_state & SDL_BUTTON(SDL_BUTTON_RIGHT))
        {
            std::cout << "Right ";
            dens[IX(grid_x,grid_y)] = 0.0f;
        }


        std::cout << "mouse: " << x << ' ' << y << '|' << grid_x << ' ' << grid_y << std::endl;
    }
}

int main(int, char **)
{
    static vfloat u(nsize, 0), v(nsize, 0), u_prev(nsize, 0), v_prev(nsize, 0); // Horizontal, vertical velocity
    static vfloat dens(nsize, 0), dens_prev(nsize, 0);
    static std::vector<bool> bounds(nsize, 0);
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

    //timeBeginPeriod(1); // Set period to 1ms
    std::chrono::time_point<std::chrono::system_clock> t_start, t_end;
    std::chrono::duration<double, std::milli> elapsed_ms;


    // Create boundary objects
    for (int i=15; i<=20; i++)
    {
        for (int j=20; j<=30; j++)
            bounds[IX(i,j)] = 1;
    }

    // Main loop
    for (int t=0; t<SIM_LEN; t++)
    {
        t_start = std::chrono::system_clock::now();

        process_input(dens_prev, dens);

        // Add some velocity
        for (int j=2*N/10.0; j<8*N/10.0; j++)
        {
            for (int i=0; i<10; i++)
                u_prev[IX(i,j)] = 200.0;
        }

        // Add some density
        for (int j=4*N/10.0; j<6*N/10.0;j++)
            dens_prev[IX(3,j)] = (t<100) ? 100.0 : 0.0;


        vel_step(u, v, u_prev, v_prev, VISC, dt, bounds);
        dens_step(dens, dens_prev, u, v, DIFF, dt, bounds);

        if (DISPLAY_CONSOLE)
        {
            std::cout << "dens" << std::endl;
            console_write(dens);
            std::cout << "u, v" << std::endl;
            console_write(u); console_write(v);
            std::cout << "dens_prev" << std::endl;
            console_write(dens_prev);
        }

        screen_draw(renderer, dens, u, v, bounds);

        t_end = std::chrono::system_clock::now();
        elapsed_ms = t_end - t_start;

        //if (elapsed_ms.count())
        //    std::cout << "ms: " << elapsed_ms.count() << '\n';

        std::this_thread::sleep_for(DELAY_LENGTH - elapsed_ms);
    }
    SDL_Quit();
    return 0;
}
