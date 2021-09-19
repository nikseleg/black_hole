#include "SDL2/SDL_main.h"
#include "SDL2/SDL.h"

#include "window.h"

#include <iostream>

void window::init() {

    SDL_Init(SDL_INIT_VIDEO);

    window = SDL_CreateWindow("Window", 
                                SDL_WINDOWPOS_UNDEFINED,
                                SDL_WINDOWPOS_UNDEFINED,
                                screen_width,
                                screen_height,
                                SDL_WINDOW_SHOWN);

    renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);

    return;
}


void window::run(SDL_Surface *image) {
    SDL_Texture *texture = SDL_CreateTextureFromSurface(renderer, image);

    SDL_Event e;
    while (1) {
        if (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT) {
                break;
            }
        }

        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);

        SDL_RenderCopy(renderer, texture, NULL, NULL);

        SDL_RenderPresent(renderer);
    }

    SDL_Quit();

    return;
}