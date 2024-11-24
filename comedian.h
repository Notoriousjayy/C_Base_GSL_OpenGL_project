//
// Created by jorda on 7/26/2024.
//

#ifndef COMEDIAN_H
#define COMEDIAN_H

#include "actor.h"

typedef struct Comedian {
    Actor base;            // Base class
    Actor* facing;         // Pointer to another actor
} Comedian;

// Initializes a Comedian instance
void initComedian(Comedian* comedian, Actor* facing);

// Comedian-specific update function
void comedianUpdate(Actor* actor);

#endif // COMEDIAN_H
