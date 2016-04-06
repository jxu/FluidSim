// GUI for Fluid Simulation
// Mouse and keyboard functions:
// Left or right click: add or remove density
// Left+right click: Add velocity
// Lshift+left or lshift+right: add or remove boundary cells

// To do: Fix framerate (see DELAY_LENGTH below)

#include "fluidsim.h"
#include <iostream>
#include <SDL2/SDL.h>
#include <stdio.h>
#include <chrono>
#include <thread>
#include <getopt.h>
#include <stdlib.h>

// Default values
int SCREEN_WIDTH = 800;
int SCREEN_HEIGHT;
int N = 50;                 // Grid size
int SIM_LEN = -1;           // Based on actual framerate

// Locks framerate at ~100, see stackoverflow.com/q/23258650/3163618
const std::chrono::milliseconds DELAY_LENGTH(10);

float VISC = 0.02;
float dt = 0.02;
float DIFF = 0.01;

// Flags
int DISPLAY_CONSOLE = false; // Console graphics
//const bool DRAW_GRID = false;
int DRAW_VEL = true;
int DEMO = 0;


float MOUSE_DENS = 100.0;
float MOUSE_VEL = 800.0;

int nsize;

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
    const float r_size = SCREEN_WIDTH / float(N+2);
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
                int color = std::min(int(dens[IX(i,j)] * 255), 255);
                SDL_SetRenderDrawColor(renderer, color, color, color, 0);
            }
            else // Object boundary
                SDL_SetRenderDrawColor(renderer, 0, 100, 200, 0);

            // Render rect
            SDL_RenderFillRect(renderer, &r);

            if (DRAW_VEL)
            {
                SDL_SetRenderDrawColor(renderer, 255, 0, 0, 0);

                int x1 = round((i+0.5)*r_size);
                int y1 = round((N+1-j+0.5)*r_size);
                int x2 = x1 + r_size*u[IX(i,j)];
                int y2 = y1 - r_size*v[IX(i,j)];
                SDL_RenderDrawLine(renderer, x1, y1, x2, y2);
            }
        }
    }

    // Render the rect to the screen
    SDL_RenderPresent(renderer);
}

// Add density (or velocity) based on user input
void process_input(vfloat &dens_prev, vfloat &dens, vfloat &u_prev, vfloat &v_prev, vbool &bound)
{
    int x, y;
    int *ptr_x = &x, *ptr_y = &y;

    float r_size = (SCREEN_WIDTH / float(N+2));

    SDL_PumpEvents();
    const unsigned int mouse_state = SDL_GetMouseState(ptr_x, ptr_y);
    const Uint8 *kb_state = SDL_GetKeyboardState(NULL);

    int grid_x = int(x/r_size);
    int grid_y = N+1 - int(y/r_size);

    if (mouse_state & (SDL_BUTTON(SDL_BUTTON_LEFT) | SDL_BUTTON(SDL_BUTTON_RIGHT)))
    {

        bool not_on_edge = (1<=grid_x && grid_x<=N && 1<=grid_y && grid_y<=N);

        if (kb_state[SDL_SCANCODE_LSHIFT])
        {
            if (mouse_state == SDL_BUTTON(SDL_BUTTON_LEFT))
                bound[IX(grid_x,grid_y)] = true;
            if (mouse_state == SDL_BUTTON(SDL_BUTTON_RIGHT))
                bound[IX(grid_x,grid_y)] = false;
        }
        else
        {
            if (mouse_state == SDL_BUTTON(SDL_BUTTON_LEFT))
            {
                //std::cout << "Left ";
                dens_prev[IX(grid_x,grid_y)] += MOUSE_DENS;
            }

            if (mouse_state == SDL_BUTTON(SDL_BUTTON_RIGHT))
            {
                //std::cout << "Right ";
                if (not_on_edge)
                {
                    dens[IX(grid_x,grid_y)]   = 0.0f;
                    dens[IX(grid_x-1,grid_y)] = 0.0f;
                    dens[IX(grid_x+1,grid_y)] = 0.0f;
                    dens[IX(grid_x,grid_y+1)] = 0.0f;
                    dens[IX(grid_x,grid_y-1)] = 0.0f;
                }
            }
        }
        //std::cout << "mouse: " << x << ' ' << y << '|' << grid_x << ' ' << grid_y << std::endl;
    }

    if (kb_state[SDL_SCANCODE_LEFT])    u_prev[IX(grid_x,grid_y)] -= MOUSE_VEL;
    if (kb_state[SDL_SCANCODE_RIGHT])   u_prev[IX(grid_x,grid_y)] += MOUSE_VEL;
    if (kb_state[SDL_SCANCODE_UP])      v_prev[IX(grid_x,grid_y)] += MOUSE_VEL;
    if (kb_state[SDL_SCANCODE_DOWN])    v_prev[IX(grid_x,grid_y)] -= MOUSE_VEL;
}

int main(int argc, char **argv)
{
    // Command line input options
    while (true)
    {
        static struct option long_options[] =
        {
            {"disp-console",    no_argument,        &DISPLAY_CONSOLE,   1},
            {"no-vel-draw",     no_argument,        &DRAW_VEL,          0},
            {"vel-draw",        no_argument,        &DRAW_VEL,          1},
            {"demo",            no_argument,        &DEMO,              1},
            {"demo2",           no_argument,        &DEMO,              2},
            {"screen-size",     required_argument,  0, 's'},
            {"grid",            required_argument,  0, 'N'},
            {"sim-len",         required_argument,  0, 'l'},
            {"visc",            required_argument,  0, 'v'},
            {"dt",              required_argument,  0, 't'},
            {"diff",            required_argument,  0, 'd'},
            {"mouse-dens",      required_argument,  0, 'm'},
            {"mouse-vel",       required_argument,  0, 'M'},
            {0, 0, 0, 0}
        };
        int option_index = 0;

        int c = getopt_long(argc, argv, "s:N:l:v:t:d:m:M:", long_options, &option_index);
        if (c == -1)
            break;

        switch (c)
        {
        case 0:
            if (long_options[option_index].flag != 0)
                break;
            printf ("option %s", long_options[option_index].name);
            if (optarg)
                printf (" with arg %s", optarg);
            printf ("\n");
            break;

        case 's':
            SCREEN_WIDTH = atoi(optarg);
            break;

        case 'N':
            N = atoi(optarg);
            break;

        case 'l':
            SIM_LEN = atoi(optarg);
            break;

        case 'v':
            VISC = atoi(optarg);
            break;

        case 't':
            dt = atoi(optarg);
            break;

        case 'd':
            DIFF = atoi(optarg);
            break;

        case 'm':
            MOUSE_DENS = atoi(optarg);
            break;

        case 'M':
            MOUSE_VEL = atoi(optarg);
            break;

        case '?':
            break;

        default:
            std::cout << "Flag case error!" << std::endl;
        }
    }

    nsize = (N+2)*(N+2);
    SCREEN_HEIGHT = SCREEN_WIDTH;

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


    if (DEMO==1)
    {
        if (N<50)
            std::cout << "N too small for demo!" << std::endl;
        // Create boundary objects
        for (int i=20; i<=25; i++)
        {
            for (int j=4; j<=30; j++)
                bound[IX(i,j)] = 1;
        }

        for (int i=20; i<=40; i++)
        {
            for (int j=4; j<=6; j++)
                bound[IX(i,j)] = 1;
        }
    }

    // Main loop
    for (unsigned int t=0; t < unsigned(SIM_LEN); t++)
    {
        t_start = std::chrono::system_clock::now();

        process_input(dens_prev, dens, u_prev, v_prev, bound);

        if (DEMO==1)
        {
             // Add some velocity
            for (int j=2*N/10.0; j<8*N/10.0; j++)
            {
                for (int i=0; i<=5; i++)
                    u_prev[IX(i,j)] = 200.0;
            }

            for (int i=20; i<=40; i++)
            {
                for (int j=1; j<4; j++)
                    u_prev[IX(i,j)] = 30.0;
            }

            for (int i=1*N/10.0; i<9*N/10.0; i++)
                u_prev[IX(i,10)] = 20.0;


            // Add some density
            for (int j=4*N/10.0; j<6*N/10.0;j++)
                dens_prev[IX(3,j)] = (t<400 && (t%40 < 20)) ? 100.0 : 0.0;
        }


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
