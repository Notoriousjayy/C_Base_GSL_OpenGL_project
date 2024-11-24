//
// Created by jorda on 7/26/2024.
//

#ifndef SCENE_H
#define SCENE_H

#include "framebuffer.h"

typedef struct Scene {
    Framebuffer buffers[2];
    Framebuffer* current;
    Framebuffer* next;
} Scene;

void initScene(Scene* scene);
void drawScene(Scene* scene);
Framebuffer* getCurrentBuffer(Scene* scene);

#endif // SCENE_H
