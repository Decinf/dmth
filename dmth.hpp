#pragma once
#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>



namespace dmth{


float max(std::vector<float> * vec){
    float m = vec->at(0);
    for(uint i = 0; i < vec->size(); ++i) { if(vec->at(i) > m) m = vec->at(i); }
    return m;
}


float min(std::vector<float> * vec){
    float m = vec->at(0);
    for(uint i = 0; i < vec->size(); ++i) { if(vec->at(i) < m) m = vec->at(i); }
    return m;
}

float min(std::vector<float> * vec, float greater_then){
    float m = vec->at(0);
    for(uint i = 0; i < vec->size(); ++i) { if(vec->at(i) < m && vec->at(i)>greater_then) m = vec->at(i); }
    return m;
}


/// @brief find index by value
/// @return index (or -1 if not found)
int findi(std::vector<float> * vec, float value, float eps = 0.0001f){
    for(uint i = 0; i < vec->size(); ++i) if(fabs(vec->at(i)-value)<eps) return i;
    return -1;
}


/// @brief count number of approximately same elements
uint countsame(std::vector<float> * vec, float value, float eps = 0.00001f){
    uint ret = 0U;
    for(uint i = 0; i < vec->size(); ++i) ret+=(fabs(vec->at(i)-value)<eps);
    return ret;
}


std::vector<float> sub_vec(std::vector<float> * a, std::vector<float> * b){
    if(a->size() != b->size()) return std::vector<float>({});
    std::vector<float> ret;
    for(uint i=0;i<a->size();++i) ret.push_back(a->at(i)-b->at(i));
    return ret;
}




inline float dfx (float(*f)(float),float x,float eps){return (f(x+eps)-f(x))/eps;}
inline float dfx (float(*f)(float,float*),float x,float eps, float*params){return (f(x+eps,params)-f(x,params))/eps;}

inline float integral(float(*f)(float),float a, float b, float eps){
    float sum = 0.0f;
    for(float x = a; x<b; x+=eps) sum += ((f(x)) + 0.5f*dfx(f,x,eps))*eps;
    return sum;
}
inline float integral(float(*f)(float,float*),float a, float b, float eps, float*params){
    float sum = 0.0f;
    for(float x = a; x<=b; x+=eps) sum += ((f(x,params)) + 0.5f*dfx(f,x,eps,params))*eps;
    return sum;
}


class vec2
{
public:
    float x, y;
    vec2(float x, float y){this->x=x; this->y=y;}
    vec2(){this->x=0.0f; this->y=0.0f;}    
    
    vec2 operator +(vec2 const & v){return vec2(x+v.x,y+v.y);}
    vec2 operator -(vec2 const & v){return vec2(x-v.x,y-v.y);}
    vec2 operator *(float a){return vec2(x*a,y*a);}
    vec2 operator /(float a){return vec2(x/a,y/a);}
};


inline float minf(float a, float b){
    return a*(float)(a<=b) + b*(float)(b<a);
}

inline float maxf(float a, float b){
    return a*(float)(a>=b) + b*(float)(b>a);
}


}