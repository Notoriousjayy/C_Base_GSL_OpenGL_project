//
// Created by jorda on 7/26/2024.
//

#include "comedian.h"

// Initializes a Comedian with a facing actor
void initComedian(Comedian* comedian, Actor* facing) {
    initActor(&comedian->base, comedianUpdate); // Use only update function
    comedian->facing = facing;
}

// Comedian-specific update function
void comedianUpdate(Actor* actor) {
    Comedian* comedian = (Comedian*)actor;  // Cast to Comedian
    if (wasActorSlapped((const Actor*)actor)) {  // Cast to const Actor* for wasSlapped
        slapActor(comedian->facing);  // Slap the facing actor
    }
}
