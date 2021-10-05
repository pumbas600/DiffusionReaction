/*
 * A Reaction-Diffusion simulation based on https://www.karlsims.com/rd.html
 * Author: Josh Jeffers
 */

/*
 * MIT License
 *
 * Copyright (c) 2021 Josh Jeffers
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <malloc.h>
#include "LibBMP.h"

#define ITERATIONS 10000
#define IterationCapture 200
#define WIDTH 600
#define HEIGHT 600

// Simulation Settings
#define DIFFUSION_RATE_A 1.0F
#define DIFFUSION_RATE_B 0.5F
#define TIME_STEP 1.0F
#define ADJACENT_NEIGHBOUR_WEIGHT 0.2F
#define DIAGONAL_NEIGHBOUR_WEIGHT 0.05F
#define CENTRE_CIRCLE_RADIUS 20

typedef struct {
    float a, b;
} Cell;

static int iteration;
static const Colour colourA = { 0, 0, 0 };
static const Colour colourB = { 50, 230, 255 };

#define KILL_RATE_MAX 0.062F
#define KILL_RATE_MIN 0.062F

#define FEED_RATE_MAX 0.0545F
#define FEED_RATE_MIN 0.0545F

float getKillRate(int x, int y) {
    return KILL_RATE_MAX;
    //return ((float)x / WIDTH) * (KILL_RATE_MAX - KILL_RATE_MIN) + KILL_RATE_MIN;
}

float getFeedRate(int x, int y) {
    return FEED_RATE_MAX;
    //return ((float)y / HEIGHT) * (FEED_RATE_MAX - FEED_RATE_MIN) + FEED_RATE_MIN;
}

bool isWithinGrid(int x, int y) {
    return x >= 0 && x < WIDTH&& y >= 0 && y < HEIGHT;
}

int xyToIndex(int x, int y) {
    return x + y * WIDTH;
}

Cell calculateDifferenceBetweenCellAndNeighbours(Cell* currentGrid, Cell* cell, int xPos, int yPos) {
    Cell difference = { 0, 0 };
    for (int x = -1; x <= 1; x++) {
        for (int y = -1; y <= 1; y++) {
            if (x == 0 && y == 0) {
                difference.a -= cell->a;
                difference.b -= cell->b;
                continue;
            }

            int gridX = xPos + x;
            int gridY = yPos + y;

            if (isWithinGrid(gridX, gridY)) {
                Cell* gridCell = &currentGrid[xyToIndex(gridX, gridY)];
                float weight = (x == 0 || y == 0) ? ADJACENT_NEIGHBOUR_WEIGHT : DIAGONAL_NEIGHBOUR_WEIGHT;
                difference.a += weight * gridCell->a;
                difference.b += weight * gridCell->b;
            }
        }
    }
    return difference;
}

void calculateNewValue(Cell* currentGrid, Cell* newValue, int x, int y) {
    Cell* cell = &currentGrid[xyToIndex(x,y)];
    Cell difference = calculateDifferenceBetweenCellAndNeighbours(currentGrid, cell, x, y);
    float reaction = cell->a * cell->b * cell->b;
    float feedRate = getFeedRate(x, y);
    float killRate = getKillRate(x, y);

    newValue->a = cell->a + (DIFFUSION_RATE_A * difference.a - reaction + feedRate * (1 - cell->a)) * TIME_STEP;
    newValue->b = cell->b + (DIFFUSION_RATE_B * difference.b + reaction - (killRate + feedRate) * cell->b) * TIME_STEP;
}


void updateGrid(Cell* currentGrid, Cell* newGrid) {
    // Generate new values. Note this is stored in a temporary array as the new value depends on the old value of
    // surrounding cells, which may have already been updated otherwise.
    for (int x = 0; x < WIDTH; x++) {
        for (int y = 0; y < HEIGHT; y++) {
            calculateNewValue(currentGrid, &newGrid[xyToIndex(x,y)], x, y);
        }
    }
}

int square(int x) {
    return x * x;
}

bool withinDistanceOfCentre(int x, int y, int distance) {
    return square(x - WIDTH / 2) + square(y - HEIGHT / 2) < square(distance);
}


void seedGrid(Cell* currentGrid) {
    for (int x = 0; x < WIDTH; x++) {
        for (int y = 0; y < HEIGHT; y++) {
            int index = xyToIndex(x, y);
            // Check if this pixel is within the centre B circle
            if (withinDistanceOfCentre(x, y, CENTRE_CIRCLE_RADIUS)) {
                currentGrid[index].a = 0.0F;
                currentGrid[index].b = 1.0F;
            }
            else {
                currentGrid[index].a = 1.0F;
                currentGrid[index].b = 0.0F;
            }
        }
    }
}

Colour lerp(const Colour* a, const Colour* b, float t) {
    Colour colour = {
        a->red + (byte)((b->red - a->red) * t),
        a->green + (byte)((b->green - a->green) * t),
        a->blue + (byte)((b->blue - a->blue) * t)
    };
    return colour;
}

void saveSimulation(Cell* currentGrid, int iteration) {
    for (int x = 0; x < WIDTH; x++) {
        for (int y = 0; y < HEIGHT; y++) {
            Cell* cell = &currentGrid[xyToIndex(x, y)];

            float denominator = cell->b + cell->a;
            if (denominator != 0) {
                Colour pixelColour = lerp(&colourA, &colourB, cell->b / denominator);
                DrawPixel(y, x, pixelColour.red, pixelColour.green, pixelColour.blue);
            }
        }
    }
    char filename[50] = "Output";
    char iterationStr[30];

    sprintf(iterationStr, "%d", iteration);
    strcat(filename, iterationStr);
    strcat(filename, ".bmp");
    SaveBMPFile(filename, WIDTH, HEIGHT);
}

Cell* createGrid() {
    Cell* grid = (Cell*)malloc(WIDTH * HEIGHT * sizeof(Cell));
    return grid;
}

void startSimulation() {
    Cell* currentGrid = createGrid();;
    Cell* newGrid = createGrid();

    int width, height;
    LoadBMPFile("blank.bmp", &width, &height);
    seedGrid(currentGrid);

    for (iteration = 0; iteration < ITERATIONS; iteration++) {
        if (iteration % 2 == 0) {
            updateGrid(currentGrid, newGrid);
            if (iteration % IterationCapture == 0) {
                saveSimulation(currentGrid, iteration);
            }
        }
        else {
            updateGrid(newGrid, currentGrid);
            if (iteration % IterationCapture == 0) {
                saveSimulation(newGrid, iteration);
            }
        }
    }
    saveSimulation(ITERATIONS % 2 == 0 ? currentGrid : newGrid, ITERATIONS);
    //printf("%f %f\n", currentGrid[0][0].a, currentGrid[0][0].b);
}

int main()
{
    startSimulation();
    printf("Simulation finished!\n");
}
