//
// Created by jorda on 7/26/2024.
//

#include <stddef.h>
#include "stage.h"

// Initializes a Stage instance
void initStage(Stage* stage) {
    for (int i = 0; i < NUM_ACTORS; i++) {
        stage->actors[i] = NULL;
    }
}

// Adds an Actor to the Stage at the specified index
void addActorToStage(Stage* stage, Actor* actor, int index) {
    if (index >= 0 && index < NUM_ACTORS) {
        stage->actors[index] = actor;
    }
}

// Updates all Actors on the Stage
void updateStage(Stage* stage) {
    for (int i = 0; i < NUM_ACTORS; i++) {
        if (stage->actors[i] != NULL) {
            stage->actors[i]->update(stage->actors[i]);
            stage->actors[i]->swap(stage->actors[i]);
        }
    }
}
