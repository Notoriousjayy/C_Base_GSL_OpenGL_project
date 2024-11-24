//
// Created by jorda on 7/26/2024.
//

#ifndef ACTOR_H
#define ACTOR_H

#include <stdbool.h>

// Define the Actor structure with function pointers for polymorphism
typedef struct Actor {
    bool currentSlapped;
    bool nextSlapped;
    void (*update)(struct Actor*);
    void (*swap)(struct Actor*);   // Function pointer for swapping states
    void (*slap)(struct Actor*);   // Function pointer for setting nextSlapped
    bool (*wasSlapped)(const struct Actor*); // Function pointer to check currentSlapped
} Actor;

// Initializes an Actor instance
void initActor(Actor* actor, void (*updateFunc)(Actor*));

// Swaps the buffer states
void swapActor(Actor* actor);

// Sets the Actor as slapped for the next buffer
void slapActor(Actor* actor);

// Checks if the Actor was slapped in the current buffer
bool wasActorSlapped(const Actor* actor);

#endif // ACTOR_H
