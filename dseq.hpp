#pragma once

#include "dmth.hpp"

#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>



namespace dmth{



float solve(float (*f)(float,float*),float x0,float eps=0.001f,float y = 0.0f,float*params=nullptr){
    float xn,xn1;
    xn = x0;
    xn1 = xn - f(xn,params)/dfx(f,xn,eps,params);

    while (std::fabs(xn1-xn) > eps)
    { 
        std::cout << '.';
        xn = xn1;
        xn1 = xn - f(xn,params)/dfx(f,xn,eps,params);
    }
    return xn1;
}

/*/// @brief solve function with parameter support 
float solve(float (*f)(float,float*),float x0,float eps,float y = 0.0f,float * params){
    float df = eps+1.0f,xn,xn1;
    xn = x0;
    xn1 = f(xn,params)/dfx(f,xn,eps/10.0f,params);

    while (std::fabs(xn1-xn) > eps)
    { 
        xn = xn1;
        xn1 = xn - f(xn,params)/dfx(f,xn,eps/10.0f,params);
    }
    return xn1;
}
*/


}