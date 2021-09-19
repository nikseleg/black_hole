#include "SDL2/SDL_main.h"
#include "SDL2/SDL.h"
#include "gmpxx.h"

#include "window.h"
#include "calculator.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <thread>
#include <future>
#include <mutex>

void calculator::print_progess(const unsigned int thread_number, const unsigned int pixel) {
    calculator::mut.lock();

    printf("Thread %u: %u\n", thread_number, pixel);

    calculator::mut.unlock();
}

std::vector<unsigned int> calculator::render(std::vector<size_t> x, std::vector<size_t> y, const size_t thread_number) {
    std::vector<unsigned int> color(x.size(), 0x00000000);

    std::vector<mpf_class> polar_pos(4), polar_vel(4), polar_acc(4);
    std::vector<mpf_class> cart_pos(4), cart_vel(4), cart_acc(4);
    mpf_class cart_pos_mag, cart_vel_mag, v_r0_mag, v_phi0_mag, w_mag, axis_mag;
    std::vector<mpf_class> v_r0(3), v_phi0(3), axis(3), w(3);

    double phi{0};

    mpf_class k_1("1.0"), k_2("1.0");

    unsigned long step;

    for (size_t i=0; i < x.size(); i++) {
        if (i % 1000 == 0) {
            calculator::print_progess(thread_number, i);
        } 

        cart_pos[0] = 0;                                                                                            // t
        cart_pos[1] = calculator::camera::plane_distance;                                                           // x
        cart_pos[2] = (x[i] - window::screen_width/2.0)/window::screen_width * calculator::camera::plane_width;     // y
        cart_pos[3] = -(y[i] - window::screen_height/2.0)/window::screen_height * calculator::camera::plane_height; // z

        cart_vel[0] = 0;
        cart_vel[1] = -calculator::camera::camera_distance;
        cart_vel[2] = cart_pos[2];
        cart_vel[3] = cart_pos[3];   

        cart_vel_mag = sqrt(cart_vel[1]*cart_vel[1] + cart_vel[2]*cart_vel[2] + cart_vel[3]*cart_vel[3]);
        
        cart_vel[0] = 0;
        cart_vel[1] = (cart_vel[1] / cart_vel_mag);
        cart_vel[2] = (cart_vel[2] / cart_vel_mag);
        cart_vel[3] = (cart_vel[3] / cart_vel_mag);


        v_r0[0] = cart_pos[1];
        v_r0[1] = cart_pos[2];
        v_r0[2] = cart_pos[3];  
    
        v_r0_mag = sqrt(v_r0[0]*v_r0[0] + v_r0[1]*v_r0[1] + v_r0[2]*v_r0[2]);
        
        v_r0[0] = v_r0[0] / v_r0_mag;
        v_r0[1] = v_r0[1] / v_r0_mag;
        v_r0[2] = v_r0[2] / v_r0_mag;

        axis[0] = 0;
        axis[1] = -v_r0[0]*v_r0[2];
        axis[2] = v_r0[0]*v_r0[1];  
        
        axis_mag = sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);

        if (axis_mag.get_d() != 0) {
            axis[0] = axis[0] / axis_mag;
            axis[1] = axis[1] / axis_mag;
            axis[2] = axis[2] / axis_mag;
        } else {
            axis[0] = 0;
            axis[1] = 0;
            axis[2] = 0;
        }

        w[0] = axis[1]*v_r0[2] - axis[2]*v_r0[1];
        w[1] = axis[2]*v_r0[0] - axis[0]*v_r0[2];
        w[2] = axis[0]*v_r0[1] - axis[1]*v_r0[0];

        v_phi0[0] = w[0];
        v_phi0[1] = w[1];
        v_phi0[2] = w[2]; 

        polar_pos[0] = 0;
        polar_pos[1] = sqrt(cart_pos[1]*cart_pos[1] + cart_pos[2]*cart_pos[2] + cart_pos[3]*cart_pos[3]);
        polar_pos[2] = PI/2;
        polar_pos[3] = 0;
        
        k_2 = polar_pos[1]*polar_pos[1] * ((v_phi0[0]*cart_vel[1] + v_phi0[1]*cart_vel[2] + v_phi0[2]*cart_vel[3])/ polar_pos[1]);  

        polar_vel[3] = k_2 / (polar_pos[1]*polar_pos[1]);
        polar_vel[2] = 0;
        polar_vel[1] = (v_r0[0]*cart_vel[1] + v_r0[1]*cart_vel[2] + v_r0[2]*cart_vel[3]);
        polar_vel[0] = sqrt(polar_pos[1]*polar_pos[1] * (polar_vel[3]*polar_vel[3] * polar_pos[1]*(polar_pos[1] - calculator::r_s) + polar_vel[1]*polar_vel[1])) / abs(polar_pos[1] - calculator::r_s);
        
        k_1 = polar_vel[0] - (calculator::r_s * polar_vel[0]) / polar_pos[1]; 

        step = 0;

        while (step < calculator::max_step) {
            
            polar_vel[0] = k_1 / (1 - calculator::r_s / polar_pos[1]);
            polar_vel[3] = k_2 / (polar_pos[1]*polar_pos[1]);

            polar_acc[1] = (polar_pos[1] - 1.5 * calculator::r_s) * (polar_vel[3]*polar_vel[3]);

            polar_vel[1] += polar_acc[1] * calculator::lambda;
            
            polar_pos[0] += polar_vel[0] * calculator::lambda;
            polar_pos[1] += polar_vel[1] * calculator::lambda;
            polar_pos[3] += polar_vel[3] * calculator::lambda;

            if (polar_pos[1].get_d() < calculator::r_s.get_d()) {
                   break;
            } else if (step > 0) {
                if (polar_pos[1].get_d() > calculator::outer_boundary.get_d()) {
                    break;
                }
            }
            if (step > calculator::max_step-2) {
                std::cout << "Step cap reached" << std::endl;
            }
            step++;
        }

        w[0] = axis[1]*cart_pos[3] - axis[2]*cart_pos[2];
        w[1] = axis[2]*cart_pos[1] - axis[0]*cart_pos[3];
        w[2] = axis[0]*cart_pos[2] - axis[1]*cart_pos[1];
        
        w_mag = sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);  
        
        cart_pos_mag = sqrt(cart_pos[1]*cart_pos[1] + cart_pos[2]*cart_pos[2] + cart_pos[3]*cart_pos[3]); 
        
        if (w_mag != 0) {
            cart_pos[0] = polar_pos[0];
            cart_pos[1] = cos(polar_pos[3].get_d())*cart_pos[1] + cart_pos_mag*(sin(polar_pos[3].get_d())/w_mag)*w[0];
            cart_pos[2] = cos(polar_pos[3].get_d())*cart_pos[2] + cart_pos_mag*(sin(polar_pos[3].get_d())/w_mag)*w[1];
            cart_pos[3] = cos(polar_pos[3].get_d())*cart_pos[3] + cart_pos_mag*(sin(polar_pos[3].get_d())/w_mag)*w[2];  
        } else {
            cart_pos[0] = polar_pos[0];
            cart_pos[1] = 0;
            cart_pos[2] = 0;
            cart_pos[3] = 0;
        }

        cart_pos_mag = sqrt(cart_pos[1]*cart_pos[1] + cart_pos[2]*cart_pos[2] + cart_pos[3]*cart_pos[3]);
        
        if (cart_pos_mag != 0) {
            cart_pos[1] = (cart_pos[1] / cart_pos_mag) * polar_pos[1];
            cart_pos[2] = (cart_pos[2] / cart_pos_mag) * polar_pos[1];
            cart_pos[3] = (cart_pos[3] / cart_pos_mag) * polar_pos[1];
        } else {
            cart_pos[1] = 0;
            cart_pos[2] = 0;
            cart_pos[3] = 0;
        }

        phi = atan2(cart_pos[2].get_d(),cart_pos[1].get_d()); 
        
        if (phi < 0) {
            phi = 2*PI + phi;
        }

        color[i] = 0x000000ff;

        if (polar_pos[1].get_d() < calculator::r_s.get_d()) {
            color[i] = 0x000000ff;
        } else {
            if (phi >= 0 && phi < 1 * ((2*PI)/6)) {
                color[i] += 0xff000000;
                color[i] += round(1.0/((2*PI)/6) * phi * 0xff) * 0x00010000;
            } else if (phi >= 1*((2*PI)/6) && phi < 2*((2*PI)/6)) {
                color[i] += 0x00ff0000;
                color[i] += round(-1.0/((2*PI)/6) * (phi - (2*PI)/6) * 0xff + 0xff) * 0x01000000;
            } else if (phi >= 2*((2*PI)/6) && phi < 3*((2*PI)/6)) {
                color[i] += 0x00ff0000;
                color[i] += round(1.0/((2*PI)/6) * (phi - 2*((2*PI)/6)) * 0xff) * 0x00000100;
            } else if (phi >= 3*((2*PI)/6) && phi < 4*((2*PI)/6)) {
                color[i] += 0x0000ff00;
                color[i] += round(-1.0/((2*PI)/6) * (phi - 3*((2*PI)/6)) * 0xff + 0xff) * 0x00010000;
            } else if (phi >= 4*((2*PI)/6) && phi < 5*((2*PI)/6)) {
                color[i] += 0x0000ff00;
                color[i] += round(1.0/((2*PI)/6) * (phi - 4*((2*PI)/6)) * 0xff) * 0x01000000;
            } else if (phi >= 5*((2*PI)/6) && phi < 6*((2*PI)/6)) {
                color[i] += 0xff000000;
                color[i] += round(-1.0/((2*PI)/6) * (phi - 5*((2*PI)/6)) * 0xff + 0xff) * 0x00000100;
            } else {
                color[i] += 0x00000000;
            }
        }
    }

    return color; 
}

void calculator::image(SDL_Surface *image) {
    SDL_SetSurfaceRLE(image, 1);
    SDL_LockSurface(image);

    unsigned int *pixel = (unsigned int*)(image->pixels);

    std::vector<std::vector<size_t>> x(8);
    std::vector<std::vector<size_t>> y(8);

    unsigned int mod;

    for (size_t yi=0; yi < window::screen_height; yi++) {
        for (size_t xi=0; xi < window::screen_width; xi++) {
            mod = (xi + yi * window::screen_width) % 8;
            x[mod].push_back(xi);
            y[mod].push_back(yi);
            
        }
    }

    
    std::vector<std::future<std::vector<unsigned int>>> color_future(8);

    for (size_t n=0; n < 8; n++) {
        color_future[n] = std::async(std::launch::async, calculator::render, x[n], y[n], n+1);
    }

    std::vector<std::vector<unsigned int>> color(8, std::vector<unsigned int>(115200, 0));

    for (size_t n=0; n < 8; n++) {
        color[n] = color_future[n].get();
    }
    

    //std::vector<std::vector<unsigned int>> color(1, std::vector<unsigned int>(921600, 0));

    //color[0] = calculator::render(x[0], y[0]);

    for (size_t n=0; n < 8; n++) {
        for (size_t i=0; i < color[n].size(); i++) {
            pixel[x[n][i] + y[n][i] * window::screen_width] = color[n][i];
        }
    }

    SDL_UnlockSurface(image);

    return;
}