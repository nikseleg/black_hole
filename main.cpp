
#include "window.h"
#include "calculator.h"

#include <iostream>
#include <vector>
#include <cmath>

int main(int argc, char const *argv[]) {
    window::init();

    SDL_Surface *image;
    image = SDL_CreateRGBSurface(0, window::screen_width, window::screen_height, 32, 0xff000000, 0x00ff0000, 0x0000ff00, 0x000000ff);

    calculator::image(image);

    window::run(image);

    return 0;
}