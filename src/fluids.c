/*
 * TODO: make struct for simmulation containung all fields and size of the field
 */

#include <stdio.h>

// swap 2 variables or pointers or anything
#define SWAP(x, y) do { typeof(x) temp##x##y = x; x = y; y = temp##x##y; } while (0)

#define SIZE_X 512
#define SIZE_Y 512
#define TIMESTEP 0.2
#define DIFF 0.5


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

void set_density_boundary(float **density) {

}

void set_velocity_boundary(vector_t **velocity) {

}

void apply_diffusion(float **density, float **density_old, float diff, float dt) {
    int i, j, k;
    // TODO: diffusion rate?
    float a = dt*diff*SIZE_X*SIZE_Y;

    // TODO: explain Gauss-Seidel solving below
    for(k=0; k<20; k++) {
        for(i=1; i<=SIZE_X; i++) {
            for(j=1; j<=SIZE_Y; j++) {
                // calculate diffusion for every pixel
                density[i][j] = ( density_old[i][j] +
                                  a * ( density[i-1][j] +
                                        density[i+1][j] +
                                        density[i][j-1] +
                                        density[i][j+1] )
                                ) / (1+4*a);
            }
        }
        // reset boundary after diffuse
        set_density_boundary(density);
    }
}

void apply_advection(float **density, float **density_old, vector_t **velocity, float dt) {
    int i,j;

    float dt_x = dt*SIZE_X;
    float dt_y = dt*SIZE_Y;

    for(i=1; i<=SIZE_X; i++) {
        for(j=1; j<=SIZE_Y; j++) {

            float x = i - dt_x * velocity[i][j].x;
            float y = j - dt_y * velocity[i][j].y;

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

            density[i][j] = s0 * (t0 * density_old[i0][j0] +
                                  t1 * density_old[i0][j1]) +
                            s1 * (t0 * density_old[i1][j0] +
                                  t1 * density_old[i1][j1]);

        }
    }
    set_density_boundary(density);
}

void density_step(float **density, float **density_old, vector_t **velocity, float diff, float dt) {
    // apply old density field
    add_float_field(density, density_old, dt);

    // move result into old field
    SWAP(density_old, density);

    // apply diffusion to 2nd field with results from prev operation as oldfield
    apply_diffusion(density, density_old, diff, dt);

    // move result again
    SWAP(density_old, density);

    // apply advection in same way diffusion is applied
    apply_advection(density, density_old, velocity, dt);
}

void velocity_step(float **velocity, float **forcefield, float visc, float dt) {

}

int main(int argc, char *argv) {
    vector_t velocity[SIZE_X+2][SIZE_Y+2];
    vector_t velocity_old[SIZE_X+2][SIZE_Y+2];
    float density[SIZE_X+2][SIZE_Y+2];
    float density_old[SIZE_X+2][SIZE_Y+2];



    return 0;
}
