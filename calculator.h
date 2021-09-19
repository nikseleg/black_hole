#pragma once

#include "gmpxx.h"

#include "window.h"

#include <vector>
#include <mutex>

#define PI 3.14159

namespace calculator {
    std::vector<unsigned int> render(std::vector<size_t> x, std::vector<size_t> y, const size_t thread_number);
    void image(SDL_Surface *image);
    void print_progess(const unsigned int thread_number, const unsigned pixel);

    static std::mutex mut;
    
    static const mpf_class r_s("14852.32054");
    static const mpf_class c("299792458");

    static const mpf_class outer_boundary("297046.4108");
    static const mpf_class lambda("1000");
    static const unsigned long max_step{1000000};

    namespace camera {
        static const mpf_class plane_distance("148523.2054");    
        static const mpf_class plane_width("74261.6027");      
        static const mpf_class plane_height("41772.15151875");   
        static const mpf_class camera_distance("14852.32054");   
    }
}

