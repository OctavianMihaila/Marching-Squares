// Author: APD team, except where source was noted

#include "helpers.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>
#include <math.h>

#define CONTOUR_CONFIG_COUNT    16
#define FILENAME_MAX_SIZE       50
#define STEP                    8
#define SIGMA                   200
#define RESCALE_X               2048
#define RESCALE_Y               2048

#define CLAMP(v, min, max) if(v < min) { v = min; } else if(v > max) { v = max; }
#define min(a, b) ((a < b) ? (a) : (b))

typedef struct {
    int thread_id;
    int nr_threads;
    ppm_image image;
    ppm_image ***contour_map;
    ppm_image **scaled_image;
    unsigned char ***grid;
    int step_x;
    int step_y;
    int needs_rescale;
    pthread_barrier_t *rescale_barrier;
    pthread_barrier_t *sample_grid_barrier;
} arguments;


// Creates a map between the binary configuration (e.g. 0110_2) and the corresponding pixels
// that need to be set on the output image. An array is used for this map since the keys are
// binary numbers in 0-15. Contour images are located in the './contours' directory.
ppm_image **init_contour_map() {
    ppm_image **map = (ppm_image **)malloc(CONTOUR_CONFIG_COUNT * sizeof(ppm_image *));
    if (!map) {
        fprintf(stderr, "Unable to allocate memory 1\n");
        exit(1);
    }

    for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
        char filename[FILENAME_MAX_SIZE];
        sprintf(filename, "./contours/%d.ppm", i);
        map[i] = read_ppm(filename);
    }

    return map;
}

// Updates a particular section of an image with the corresponding contour pixels.
// Used to create the complete contour image.
void update_image(ppm_image *image, ppm_image *contour, int x, int y) {
    for (int i = 0; i < contour->x; i++) {
        for (int j = 0; j < contour->y; j++) {
            int contour_pixel_index = contour->x * i + j;
            int image_pixel_index = (x + i) * image->y + y + j;

            image->data[image_pixel_index].red = contour->data[contour_pixel_index].red;
            image->data[image_pixel_index].green = contour->data[contour_pixel_index].green;
            image->data[image_pixel_index].blue = contour->data[contour_pixel_index].blue;
        }
    }
}

// Corresponds to step 1 of the marching squares algorithm, which focuses on sampling the image.
// Builds a p x q grid of points with values which can be either 0 or 1, depending on how the
// pixel values compare to the `sigma` reference value. The points are taken at equal distances
// in the original image, based on the `step_x` and `step_y` arguments.
void sample_grid(ppm_image *image, unsigned char **grid, int thread_id, int nr_threads, int step_x, int step_y, unsigned char sigma) {
    int p = image->x / step_x;
    int q = image->y / step_y;

    int start = thread_id * (double)p / nr_threads;
    int end = min((thread_id + 1) * ((double)p / nr_threads), p);

    for (int i = start; i < end; i++) {
        for (int j = 0; j < q; j++) {
            ppm_pixel curr_pixel = image->data[i * step_x * image->y + j * step_y];

            unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

            if (curr_color > sigma) {
                grid[i][j] = 0;
            } else {
                grid[i][j] = 1;
            }
        }
    }

    grid[p][q] = 0;

    // last sample points have no neighbors below / to the right, so we use pixels on the
    // last row / column of the input image for them
    for (int i = start; i < end; i++) {
        ppm_pixel curr_pixel = image->data[i * step_x * image->y + image->x - 1];

        unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

        if (curr_color > sigma) {
            grid[i][q] = 0;
        } else {
            grid[i][q] = 1;
        }
    }

    for (int j = 0; j < q; j++) {
        ppm_pixel curr_pixel = image->data[(image->x - 1) * image->y + j * step_y];

        unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

        if (curr_color > sigma) {
            grid[p][j] = 0;
        } else {
            grid[p][j] = 1;
        }
    }

    return;
}


// Corresponds to step 2 of the marching squares algorithm, which focuses on identifying the
// type of contour which corresponds to each subgrid. It determines the binary value of each
// sample fragment of the original image and replaces the pixels in the original image with
// the pixels of the corresponding contour image accordingly.
void march(ppm_image *image, unsigned char **grid, ppm_image **contour_map, int thread_id, int nr_threads, int step_x, int step_y) {
    int p = image->x / step_x;
    int q = image->y / step_y;

    int start = thread_id * (double)p / nr_threads;
    int end = min((thread_id + 1) * ((double)p / nr_threads), p);

    for (int i = start; i < end; i++) {
        for (int j = 0; j < q; j++) {
            unsigned char k = 8 * grid[i][j] + 4 * grid[i][j + 1] + 2 * grid[i + 1][j + 1] + 1 * grid[i + 1][j];
            update_image(image, contour_map[k], i * step_x, j * step_y);
        }
    }
}

// Calls `free` method on the utilized resources.
void free_resources(ppm_image *image, ppm_image **contour_map, unsigned char **grid, int step_x) {
    for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
        free(contour_map[i]->data);
        free(contour_map[i]);
    }
    free(contour_map);

    for (int i = 0; i <= image->x / step_x; i++) {
        free(grid[i]);
    }
    free(grid);

    free(image->data);
    free(image);
}


void *rescale_image(ppm_image *image, ppm_image *scaled_image, int thread_id, int nr_threads) {
    uint8_t sample[3];

    int p = RESCALE_X;
    int start = thread_id * (double)p / nr_threads;
    int end = min((thread_id + 1) * ((double)p / nr_threads), p);

    for (int i = start; i < end; i++) {
        for (int j = 0; j < scaled_image->y; j++) {
            float u = (float)i / (float)(RESCALE_X - 1);
            float v = (float)j / (float)(RESCALE_Y - 1);
            
            // use bicubic interpolation for scaling
            sample_bicubic(image, u, v, sample);

            scaled_image->data[i * RESCALE_Y + j].red = sample[0];
            scaled_image->data[i * RESCALE_Y + j].green = sample[1];
            scaled_image->data[i * RESCALE_Y + j].blue = sample[2];
        }
    }

    return NULL;
}

void *routine_to_parallelize(void *arg) {
    arguments *data = (arguments *)arg;

    int thread_id = data->thread_id;
    int nr_threads = data->nr_threads;
    int step_x = data->step_x;
    int step_y = data->step_y;
    int needs_rescale = data->needs_rescale;

    unsigned char ***grid = data->grid;

    ppm_image *image = &(data->image);
    ppm_image **scaled_image = data->scaled_image;
    ppm_image ***contour_map = data->contour_map;

    pthread_barrier_t *rescale_barrier = data->rescale_barrier;
    pthread_barrier_t *sample_grid_barrier = data->sample_grid_barrier;
    
    // Rescale only when the image is too large (downscale).
    if (needs_rescale == 1) {
        rescale_image(image, *scaled_image, thread_id, nr_threads);
        pthread_barrier_wait(rescale_barrier);
    }

    // Perform step 1 of the marching squares algorithm
    sample_grid(*scaled_image, &(**grid), thread_id, nr_threads, step_x, step_y, SIGMA);
    pthread_barrier_wait(sample_grid_barrier);

    // Perform step 2 of the marching squares algorithm
    march(*scaled_image, &(**grid), &(**contour_map), thread_id, nr_threads, step_x, step_y);

    pthread_exit(NULL);
}

int main(int argc, char *argv[]) {
    if (argc < 4) {
        fprintf(stderr, "Usage: ./tema1 <in_file> <out_file> <P>\n");
        return 1;
    }

    ppm_image *image = read_ppm(argv[1]);
    
    int step_x = STEP;
    int step_y = STEP;
    int needs_rescale = 1;
    int nr_threads = atoi(argv[3]);

    pthread_t threads[nr_threads];
    pthread_barrier_t rescale_barrier;
    pthread_barrier_t sample_grid_barrier;
	
    pthread_barrier_init(&rescale_barrier, NULL, nr_threads);
	pthread_barrier_init(&sample_grid_barrier, NULL, nr_threads);

    // Initialize contour map
    ppm_image **contour_map = init_contour_map();
    ppm_image *scaled_image = NULL;

    // Check if the image needs to be scaled
    if (image->x <= RESCALE_X && image->y <= RESCALE_Y) {
        scaled_image = image;
        needs_rescale = 0;
    } else { // case when the image needs to be scaled
        scaled_image = malloc(sizeof(ppm_image));
        
        if (!scaled_image) {
            fprintf(stderr, "Memory allocation failed\n");
            exit(1);
        }

        scaled_image->x = RESCALE_X;
        scaled_image->y = RESCALE_Y;

        scaled_image->data = (ppm_pixel*)malloc(scaled_image->x * scaled_image->y * sizeof(ppm_pixel));
        if (!(scaled_image->data)) {
            fprintf(stderr, "Unable to allocate memory 4\n");
            exit(1);
        }
    }
   
    int p = scaled_image->x / step_x;
    int q = scaled_image->y / step_y;

    unsigned char **grid = (unsigned char **)malloc((p + 1) * sizeof(unsigned char*));
    if (!grid) {
        fprintf(stderr, "Unable to allocate memory 2\n");
        exit(1);
    }

    for (int i = 0; i <= p; i++) {
        grid[i] = (unsigned char *)malloc((q + 1) * sizeof(unsigned char));
        if (!grid[i]) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }
    }

    // Parallelize the marching squares algorithm
    for (int i = 0; i < nr_threads; i++) {
        arguments *arg = (arguments *)malloc(sizeof(arguments));

        if (!arg) {
            fprintf(stderr, "Memory allocation failed(arguments).\n");
            exit(1);
        }

        arg->thread_id = i;
        arg->nr_threads = nr_threads;
        arg->image = *image;
        arg->scaled_image = &scaled_image;
        arg->contour_map = &contour_map;
        arg->grid = &grid;
        arg->step_x = step_x;
        arg->step_y = step_y;
        arg->needs_rescale = needs_rescale;
        arg->rescale_barrier = &rescale_barrier;
        arg->sample_grid_barrier = &sample_grid_barrier;

        int r = pthread_create(&threads[i], NULL, routine_to_parallelize, arg);

        if (r) {
            fprintf(stderr, "Error creating thread %d\n", i);
            exit(1);
        }
    }

    for (int i = 0; i < nr_threads; i++) {
        int r = pthread_join(threads[i], NULL);

        if (r) {
            fprintf(stderr, "Error joining thread %d\n", i);
            exit(1);
        }
    }

    // Write output
    write_ppm(scaled_image, argv[2]);

    pthread_barrier_destroy(&rescale_barrier);
    pthread_barrier_destroy(&sample_grid_barrier);

    free_resources(scaled_image, contour_map, grid, step_x);

    return 0;
}
