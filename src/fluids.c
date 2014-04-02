/*
 * TODO: make struct for simmulation containung all fields and size of the field
 */

#include <stdio.h>
#include <stdlib.h>
#include <SDL/SDL.h>
#include <signal.h>

// swap 2 variables or pointers or anything
#define SWAP(x, y) do { typeof(x) temp##x##y = x; x = y; y = temp##x##y; } while (0)

#define SIZE_X 512
#define SIZE_Y 512
#define TIMESTEP 0.2
#define DIFF 0.5

// do you whant me to end the programm?
volatile int stop = 0;

// callback to catch signals
void signal_callback(int i) {
        stop=1;
}

typedef struct vector {
    float x;
    float y;
} vector_t;

void add_float_field(float **field, float **sources, float dt) {
    int i,j;
    for(i=0; i<SIZE_X+2; i++) {
        for(j=0; i<SIZE_Y+2; j++) {
            field[i][j] += sources[i][j] * dt;
        }
    }
}

void add_vector_field(vector_t **field, vector_t **sources, float dt) {
    int i,j;
    for(i=0; i<SIZE_X+2; i++) {
        for(j=0; i<SIZE_Y+2; j++) {
            field[i][j].x += sources[i][j].x * dt;
            field[i][j].y += sources[i][j].y * dt;
        }
    }
}

void set_field_boundary(float **density) {

}

void set_vector_field_boundary(vector_t **velocity) {

}

void set_forcefield_boundary(vector_t **forcefield, int selector) {
    if(selector == 2) {
        //do both
    }
    else {
        //only do y part
    }

}

void apply_diffusion_1d(float **field, float **field_old, float diff, float dt) {
    int i, j, k;
    // TODO: diffusion rate?
    float a = dt*diff*SIZE_X*SIZE_Y;

    // TODO: explain Gauss-Seidel solving below
    for(k=0; k<20; k++) {
        for(i=1; i<=SIZE_X; i++) {
            for(j=1; j<=SIZE_Y; j++) {
                // calculate diffusion for every pixel
                field[i][j] = ( field_old[i][j] +
                               a * ( field[i-1][j] +
                                     field[i+1][j] +
                                     field[i][j-1] +
                                     field[i][j+1]
                                   )
                              ) / (1+4*a);
            }
        }
        // reset boundary after diffuse
        set_field_boundary(field);
    }
}


void apply_diffusion_2d(vector_t **vector_field, vector_t **vector_field_old, float diff, float dt) {
    int i, j, k;
    // TODO: diffusion rate?
    float a = dt*diff*SIZE_X*SIZE_Y;

    // TODO: explain Gauss-Seidel solving below
    for(k=0; k<20; k++) {
        for(i=1; i<=SIZE_X; i++) {
            for(j=1; j<=SIZE_Y; j++) {
                // calculate diffusion for every vector
                vector_field[i][j].x = ( vector_field_old[i][j].x +
                                        a * ( vector_field[i-1][j].x +
                                              vector_field[i+1][j].x +
                                              vector_field[i][j-1].x +
                                              vector_field[i][j+1].x
                                            )
                                       ) / (1+4*a);

                vector_field[i][j].y = ( vector_field_old[i][j].y +
                                        a * ( vector_field[i-1][j].y +
                                              vector_field[i+1][j].y +
                                              vector_field[i][j-1].y +
                                              vector_field[i][j+1].y
                                            )
                                       ) / (1+4*a);
            }
        }
        // reset boundary after diffuse
        set_vector_field_boundary(vector_field);
    }
}

void apply_advection_1d(float **field, float **field_old, vector_t **vector_field, float dt) {
    int i,j;

    float dt_x = dt*SIZE_X;
    float dt_y = dt*SIZE_Y;

    for(i=1; i<=SIZE_X; i++) {
        for(j=1; j<=SIZE_Y; j++) {

            float x = i - dt_x * vector_field[i][j].x;
            float y = j - dt_y * vector_field[i][j].y;

            // check if x coordinate is out of range
            if (x < 0.5) x = 0.5;
            if (x > SIZE_X + 0.5) x = SIZE_X + 0.5;

            int i0 = (int) x;
            int i1 = i0 + 1;

            // check if y coordinate is out of range
            if (y < 0.5) y = 0.5;
            if (y > SIZE_Y + 0.5) x = SIZE_Y + 0.5;

            int j0 = (int) y;
            int j1 = j0 +1;

            //TODO wtf?
            float s1 = x - i0;
            float s0 = 1 - s1;
            float t1 = y - j0;
            float t0 = 1 - t1;

            field[i][j] = s0 * (t0 * field_old[i0][j0] +
                                  t1 * field_old[i0][j1]) +
                            s1 * (t0 * field_old[i1][j0] +
                                  t1 * field_old[i1][j1]);

        }
    }
    set_field_boundary(field);
}


void apply_advection_2d(vector_t **vector_field, vector_t **vector_field_old, vector_t **in_vectors, float dt) {
    int i,j;

    float dt_x = dt*SIZE_X;
    float dt_y = dt*SIZE_Y;

    for(i=1; i<=SIZE_X; i++) {
        for(j=1; j<=SIZE_Y; j++) {

            float x = i - dt_x * in_vectors[i][j].x;
            float y = j - dt_y * in_vectors[i][j].y;

            // check if x coordinate is out of range
            if (x < 0.5) x = 0.5;
            if (x > SIZE_X + 0.5) x = SIZE_X + 0.5;

            int i0 = (int) x;
            int i1 = i0 + 1;

            // check if y coordinate is out of range
            if (y < 0.5) y = 0.5;
            if (y > SIZE_Y + 0.5) x = SIZE_Y + 0.5;

            int j0 = (int) y;
            int j1 = j0 +1;

            //TODO wtf?
            float s1 = x - i0;
            float s0 = 1 - s1;
            float t1 = y - j0;
            float t0 = 1 - t1;

            vector_field[i][j].x = s0 * (t0 * vector_field_old[i0][j0].x +
                                  t1 * vector_field_old[i0][j1].x) +
                            s1 * (t0 * vector_field_old[i1][j0].x +
                                  t1 * vector_field_old[i1][j1].x);

            vector_field[i][j].y = s0 * (t0 * vector_field_old[i0][j0].y +
                                  t1 * vector_field_old[i0][j1].y) +
                            s1 * (t0 * vector_field_old[i1][j0].y +
                                  t1 * vector_field_old[i1][j1].y);
        }
    }
    set_vector_field_boundary(vector_field);
}

void density_step(float **density, float **density_old, vector_t **velocity, float diff, float dt) {
    // apply old density field
    add_float_field(density, density_old, dt);

    // move result into old field
    SWAP(density_old, density);

    // apply diffusion to 2nd field with results from prev operation as oldfield
    apply_diffusion_1d(density, density_old, diff, dt);

    // move result again
    SWAP(density_old, density);

    // apply advection in same way diffusion is applied
    apply_advection_1d(density, density_old, velocity, dt);
}

void project(vector_t **velocity, vector_t **forcefield) {
    int i, j, k;

    // TODO: discribe this
    float h_x = 1.0/SIZE_X;
    float h_y = 1.0/SIZE_Y;

    // just to keep code readable
    vector_t **p = forcefield;
    vector_t **div = forcefield;


    for(i=1; i<=SIZE_X; i++) {
        for(j=1; j<=SIZE_Y; j++) {
            div[i][j].x = -0.5 * h_x * (velocity[i+1][j].x - velocity[i-1][j].x +
                                        velocity[i][j+1].y - velocity[i][j-1].y);
            p[i][j].y = 0;
        }
    }
    set_forcefield_boundary(forcefield,2);

    // TODO: explain Gauss-Seidel solving below
    for(k=0; k<20; k++) {
        for(i=1; i<=SIZE_X; i++) {
            for(j=1; j<=SIZE_Y; j++) {
                // y -> p             // x -> div
                p[i][j].y = (div[i][j].x + p[i-1][j].y + p[i+1][j].y +
                                           p[i][j-1].y + p[i][j+1].y) / 4;
            }
        }
        set_forcefield_boundary(p,1);
    }

    for(i=1; i<=SIZE_X; i++) {
        for(j=1; j<=SIZE_Y; j++) {
            velocity[i][j].x -= 0.5 * (p[i+1][j].y - p[i-1][j].y) / h_x;
            velocity[i][j].y -= 0.5 * (p[i][j+1].y - p[i][j-1].y) / h_y;
        }
    }
    set_vector_field_boundary(velocity);
}

float **alloc_float_field(unsigned int x_len, unsigned int y_len) {
    float **field = malloc(x_len * sizeof(float *));
    field[0] = malloc(x_len * y_len * sizeof(float));
    int i;
    for(i = 1; i < x_len; i++) {
        field[i] = field[0] + i * y_len;
    }
    return field;
}

vector_t **alloc_vector_field(unsigned int x_len, unsigned int y_len) {
vector_t **field = malloc(x_len * sizeof(vector_t *));
    field[0] = malloc(x_len * y_len * sizeof(vector_t));
    int i;
    for(i = 1; i < x_len; i++) {
        field[i] = field[0] + i * y_len;
    }
    return field;
}

void free_field(void **field) {
    free((void *) field[0]);
    free((void *) field);
}

void velocity_step(vector_t **velocity, vector_t **forcefield, float viscosity, float dt) {
    add_vector_field(velocity, forcefield, dt);

    SWAP(forcefield, velocity);

    apply_diffusion_2d(velocity, forcefield, viscosity, dt);

    project(velocity, forcefield);

    SWAP(forcefield, velocity);

    // in case of fail maybe look here
    apply_advection_2d(velocity, forcefield, forcefield, dt);

    project(velocity, forcefield);
}

// function that draws a single pixel to a surface
void draw_pixel(SDL_Surface *Surface, int x, int y, Uint16 color)
{
  Uint16 *Pixel;
  Pixel = (Uint16 *)Surface->pixels + y*Surface->pitch/2 + x*2;
  *Pixel = color;
}

int main(int argc, char *argv) {
    vector_t **velocity = alloc_vector_field(SIZE_X+2,SIZE_Y+2);
    vector_t **velocity_old = alloc_vector_field(SIZE_X+2, SIZE_Y+2);
    float **density = alloc_float_field(SIZE_X+2, SIZE_Y+2);
    float **density_old = alloc_float_field(SIZE_X+2, SIZE_Y+2);

    //register signal callback
    signal(SIGINT, signal_callback);    

    // surface that will hold the background
    SDL_Surface *screen = NULL;

    // initialize SDL video subsystem
    if(  SDL_Init(SDL_INIT_VIDEO) == -1 ) {
        return EXIT_FAILURE;
    }

    // Set up the screen with WINDOW_X*WINDOW_Y px size, 32 bit color and a
    // softwaresurface
    screen = SDL_SetVideoMode(SIZE_X, SIZE_Y, 32, SDL_SWSURFACE);
    if( screen == NULL ) {
        return EXIT_FAILURE;
    }

    // Set the window caption
    SDL_WM_SetCaption("ICE v0.0.1", NULL);

    // create a black surface which will hold the desity paterns
    SDL_Surface *fluid = SDL_CreateRGBSurface(SDL_SWSURFACE, SIZE_X, SIZE_Y, 32,0,0,0,0);

    while(!stop) {
        // -- set fields --

        // caluclate fluids
        //velocity_step();
        //density_step();

        // -- draw density --
        // Apply image to screen
        SDL_BlitSurface(fluid, NULL, screen, NULL);
        // Update Screen
        SDL_Flip(screen);
    }

    // Quit SDL
    SDL_Quit();

    return 0;
}
