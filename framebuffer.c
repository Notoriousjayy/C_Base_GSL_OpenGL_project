//
// Created by jorda on 7/26/2024.
//

#include "framebuffer.h"
#include <stdlib.h>
#include <string.h>
#include <GLFW/glfw3.h>
#include <stdio.h>

Framebuffer* createFramebuffer() {
    Framebuffer* fb = (Framebuffer*)malloc(sizeof(Framebuffer));
    if (fb) {
        clearFramebuffer(fb);
    }
    return fb;
}

void clearFramebuffer(Framebuffer* fb) {
    memset(fb->pixels, 0, WIDTH * HEIGHT); // Set all pixels to black (assuming 0 is black)
}

void drawPixel(Framebuffer* fb, int x, int y) {
    if (x >= 0 && x < WIDTH && y >= 0 && y < HEIGHT) {
        fb->pixels[(WIDTH * y) + x] = 0; // Assuming 0 represents BLACK
    }
}

const char* getPixels(const Framebuffer* fb) {
    return (const char*)fb->pixels;
}

void renderFramebuffer(const Framebuffer* fb) {
    glDrawPixels(WIDTH, HEIGHT, GL_LUMINANCE, GL_UNSIGNED_BYTE, fb->pixels);
    GLenum err = glGetError();
    if (err != GL_NO_ERROR) {
        printf("OpenGL error: %d\n", err);
    }
}

