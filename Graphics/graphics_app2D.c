//
// Created by jorda on 7/20/2024.
//

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include "graphics2d.h"
#include "graphics_app2D.h"

Canvas canvas;

void initGraphics() {
    // Initialize the canvas with a window
    initCanvas(&canvas, 800, 600, "2D Graphics Window");

    // Set the canvas window and viewport
    setCanvasWindow(&canvas, 0.0f, 800.0f, 0.0f, 600.0f);
    setCanvasViewport(&canvas, 0, 800, 0, 600);
}

void render() {
    // Clear the screen
    glClear(GL_COLOR_BUFFER_BIT);

    // Example drawing: draw a line from (100, 100) to (200, 200)
    canvasMoveTo(&canvas, 100.0f, 100.0f);
    canvasLineTo(&canvas, 200.0f, 200.0f);

    // Swap buffers
    glfwSwapBuffers(glfwGetCurrentContext());
}
