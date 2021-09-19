#pragma once

#include "SDL2/SDL_main.h"
#include "SDL2/SDL.h"

#undef main

namespace window {
    static const unsigned int screen_width{1280};
    static const unsigned int screen_height{720};

    static SDL_Window *window;
    static SDL_Renderer *renderer;

    void init();
    void run(SDL_Surface *image);

}