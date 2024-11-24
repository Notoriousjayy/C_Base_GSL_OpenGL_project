//
// Created by jorda on 7/26/2024.
//

#ifndef STAGE_H
#define STAGE_H

#include "actor.h"

#define NUM_ACTORS 3

typedef struct Stage {
    Actor* actors[NUM_ACTORS];
} Stage;

// Initializes a Stage instance
void initStage(Stage* stage);

// Adds an Actor to the Stage at the specified index
void addActorToStage(Stage* stage, Actor* actor, int index);

// Updates all Actors on the Stage
void updateStage(Stage* stage);

#endif // STAGE_H
