#include "scene.h"
#include "include/graphics2d.h"
#include <stdlib.h>

void initScene(Scene* scene) {
    if (!scene) return;
    clearFramebuffer(&scene->buffers[0]);
    clearFramebuffer(&scene->buffers[1]);
    scene->current = &scene->buffers[0];
    scene->next = &scene->buffers[1];
}

void drawScene(Scene* scene) {
    if (!scene || !scene->next) return;

    // Clear the next buffer
    clearFramebuffer(scene->next);

    // Use Graphics2D functions to draw shapes
//    drawCircle(scene->next, WIDTH / 2, HEIGHT / 2, 50);
//    drawLine(scene->next, 0, 0, WIDTH, HEIGHT);

    // Swap the buffers
    Framebuffer* temp = scene->current;
    scene->current = scene->next;
    scene->next = temp;
}

Framebuffer* getCurrentBuffer(Scene* scene) {
    return scene ? scene->current : NULL;
}
