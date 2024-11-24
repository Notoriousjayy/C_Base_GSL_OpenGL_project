//
// Created by jorda on 7/26/2024.
//

#include "actor.h"

// Initializes an Actor with the provided update function
void initActor(Actor* actor, void (*updateFunc)(Actor*)) {
    actor->currentSlapped = false;
    actor->nextSlapped = false;
    actor->update = updateFunc;
    actor->swap = swapActor;
    actor->slap = slapActor;
    actor->wasSlapped = wasActorSlapped;
}

// Swaps the buffer states
void swapActor(Actor* actor) {
    actor->currentSlapped = actor->nextSlapped;
    actor->nextSlapped = false;
}

// Sets the Actor as slapped for the next buffer
void slapActor(Actor* actor) {
    actor->nextSlapped = true;
}

// Checks if the Actor was slapped in the current buffer
bool wasActorSlapped(const Actor* actor) {
    return actor->currentSlapped;
}
