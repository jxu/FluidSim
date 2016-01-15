// To do: Fix framerate (see DELAY_LENGTH below)

#include "fluidsim.h"
#include <iostream>
#include <SDL2/SDL.h>
#include <stdio.h>
#include <chrono>
#include <thread>

// Constants
const int SCREEN_WIDTH = 800;
const int SCREEN_HEIGHT = 800;  // Should match SCREEN_WIDTH
const int N = 50;               // Grid size
const int SIM_LEN = -1;         // Based on actual framerate

// Locks framerate at ~64, see stackoverflow.com/q/23258650/3163618
const std::chrono::milliseconds DELAY_LENGTH(5);

const float VISC = 0.01;
const float dt = 0.003;
const float DIFF = 0.005;

const bool DISPLAY_CONSOLE = false; // Console or graphics
const bool DRAW_GRID = false; // implement later
const bool DRAW_VEL = true;

const float MOUSE_DENS = 100.0;
const float MOUSE_VEL = 300.0;

const int nsize = (N+2)*(N+2);


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

void screen_draw(SDL_Renderer *renderer, vfloat &dens, vfloat &u, vfloat &v, vbool &bound)
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
void process_input(vfloat &dens_prev, vfloat &dens, vfloat &u_prev, vfloat &v_prev)
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
        bool not_on_edge = (1<=grid_x && grid_x<=N && 1<=grid_y && grid_y<=N);

        if (mouse_state == SDL_BUTTON(SDL_BUTTON_LEFT))
        {
            std::cout << "Left ";
            dens_prev[IX(grid_x,grid_y)] += MOUSE_DENS;
        }

        if (mouse_state == SDL_BUTTON(SDL_BUTTON_RIGHT))
        {
            std::cout << "Right ";
            dens[IX(grid_x,grid_y)] = 0.0f;
            if (not_on_edge)
            {
                dens[IX(grid_x-1,grid_y)] = 0.0f;
                dens[IX(grid_x+1,grid_y)] = 0.0f;
                dens[IX(grid_x,grid_y+1)] = 0.0f;
                dens[IX(grid_x,grid_y-1)] = 0.0f;
            }
        }

        if (mouse_state == (SDL_BUTTON(SDL_BUTTON_LEFT) | SDL_BUTTON(SDL_BUTTON_RIGHT)))
        {
            std::cout << "Left+Right ";
            if (not_on_edge)
            {
                for (int i=grid_x-1; i<=grid_x+1; i++)
                {
                    v_prev[IX(i,grid_y+1)] -= MOUSE_VEL;
                    v_prev[IX(i,grid_y-1)] += MOUSE_VEL;
                }
                for (int j=grid_y-1; j<=grid_y+1; j++)
                {
                    u_prev[IX(grid_x+1,j)] += MOUSE_VEL;
                    u_prev[IX(grid_x-1,j)] -= MOUSE_VEL;
                }
            }
        }
        std::cout << "mouse: " << x << ' ' << y << '|' << grid_x << ' ' << grid_y << std::endl;
    }
}

int main(int, char **)
{
    static vfloat u(nsize, 0), v(nsize, 0), u_prev(nsize, 0), v_prev(nsize, 0); // Horizontal, vertical velocity
    static vfloat dens(nsize, 0), dens_prev(nsize, 0);
    static vbool bound(nsize, 0);
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

    SDL_Renderer* renderer = SDL_CreateRenderer(window, 0, SDL_RENDERER_ACCELERATED);

    SDL_SetRenderDrawColor(renderer, 255, 0, 255, 255); // Background color, should not see this
    SDL_RenderClear(renderer);

    //timeBeginPeriod(1); // Set period to 1ms
    std::chrono::time_point<std::chrono::system_clock> t_start, t_end;
    std::chrono::duration<double, std::milli> elapsed_ms;


    // Create boundary objects
    for (int i=15; i<=20; i++)
    {
        for (int j=24; j<=26; j++)
            bound[IX(i,j)] = 1;
    }

    // Main loop
    for (unsigned int t=0; t< unsigned(SIM_LEN); t++)
    {
        t_start = std::chrono::system_clock::now();

        process_input(dens_prev, dens, u_prev, v_prev);

        // Add some velocity
        for (int j=2*N/10.0; j<8*N/10.0; j++)
        {
            for (int i=0; i<=5; i++)
                u_prev[IX(i,j)] = 200.0;
        }

        for (int i=1*N/10.0; i<9*N/10.0; i++)
            u_prev[IX(i,10)] = 20.0;


        // Add some density
        for (int j=4*N/10.0; j<6*N/10.0;j++)
            dens_prev[IX(3,j)] = (t<100) ? 100.0 : 0.0;


        vel_step(u, v, u_prev, v_prev, VISC, dt, bound);
        dens_step(dens, dens_prev, u, v, DIFF, dt, bound);

        if (DISPLAY_CONSOLE)
        {
            std::cout << "dens" << std::endl;
            console_write(dens);
            std::cout << "u, v" << std::endl;
            console_write(u); console_write(v);
            std::cout << "dens_prev" << std::endl;
            console_write(dens_prev);
        }

        screen_draw(renderer, dens, u, v, bound);

        t_end = std::chrono::system_clock::now();
        elapsed_ms = t_end - t_start;

        //if (elapsed_ms.count())
        //    std::cout << "ms: " << elapsed_ms.count() << '\n';

        std::this_thread::sleep_for(DELAY_LENGTH - elapsed_ms);
    }
    SDL_Quit();
    return 0;
}