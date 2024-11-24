//
// Created by jorda on 7/26/2024.
//

#ifndef FRAMEBUFFER_H
#define FRAMEBUFFER_H

#define WIDTH 1280
#define HEIGHT 720

typedef struct Framebuffer {
    unsigned char pixels[WIDTH * HEIGHT];
} Framebuffer;

Framebuffer* createFramebuffer();
void clearFramebuffer(Framebuffer* fb);
void drawPixel(Framebuffer* fb, int x, int y);
const char* getPixels(const Framebuffer* fb);
void renderFramebuffer(const Framebuffer* fb);

#endif // FRAMEBUFFER_H
