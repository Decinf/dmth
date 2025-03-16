#pragma once

#include "dmth.hpp"
#include "dseq.hpp"

#include <iostream>
#include <queue>
#include <vector>
#include <math.h>
#include <algorithm>

namespace dmth{


struct interval
{
public:
    float min,max;
    uint n;
    interval(float min, float max, uint n){
        this->min = min;
        this->max = max;
        this->n = n;
    }
    inline float mid(){return (min+max)*0.5f;}
};



/// @brief returns float random [0;1]
#define randf ((float)std::rand()/(float)RAND_MAX)



#define str(Value) #Value 

/// @brief fills vector with random float values
/// @param vec 
/// @param min 
/// @param max 
/// @param seed leave 0 for no seed
void fillvec_rnd(std::vector<float> * vec, float min, float max, uint seed=0){
    std::srand(seed + !(bool)seed*time(NULL));
    float range = max-min;
    for(uint i = 0; i < vec->size(); ++i) { vec->at(i)=min+randf*range; }
}


void fillvec_func(std::vector<float> * vec, float (*func)()){
    for(uint i = 0; i < vec->size(); ++i) { vec->at(i)=func(); }
}

void pushvec_func(std::vector<float> * vec, uint len, float (*func)()){
    for(uint i = 0; i < len; ++i) { vec->push_back(func()); }
}


std::string ftostr(float f, uint8_t pop_symbols = 0){
            std::string ret = std::to_string(f);
            for(uint8_t i = 0; i < pop_symbols; ++i)ret.pop_back();
            return ret;
};
inline std::string vec_to_str(std::vector<float> * vec,bool format = false, char div = ' ',uint8_t pop_symbols = 0){
    std::string ret = "";
    for(uint i = 0; i < vec->size(); ++i) { if(i && format) ret+=','; ret+=div+ftostr(vec->at(i),pop_symbols); }
    ret += "";
    return ret;
}

inline std::string vec_to_html(std::vector<float> * vec, std::string head){
    std::string ret = "<table><tr><td>" + head + "</td><tr>";
    for(uint i = 0; i < vec->size(); ++i) { ret+="<tr><td>"+std::to_string(vec->at(i))+"</td></tr>"; }
    ret += "</table>";
    return ret;
}


inline float sum(std::vector<float> * x){
    float ret = 0.0f; 
    for(const auto&xi : *x) { ret+=xi; }
    return ret;
}

inline float sum(std::vector<interval> * x){
    float ret = 0.0f; 
    for(const interval&i : *x) { ret+=i.n*0.5f*(i.min+i.max); }
    return ret;
}


inline uint countn(std::vector<interval> * x){
    uint ret = 0.0f; 
    for(const interval&i : *x) { ret+=i.n; }
    return ret;
}

/// @brief expectance, m(x), E[x]
inline float m(std::vector<float> * x){
    return sum(x)/x->size();
}
/// @brief expectance, m(x), E[x]
inline float m(std::vector<interval> * x){
    float s = 0;
    float sumn = 0.0f;
    for(const interval&i : *x) { sumn+=i.n; }
    return sum(x)/sumn;
}

/// @brief start moment: v(x)=m(x^n)
inline float v(std::vector<float> * x, uint n){
    float ret = 0.0f; 
    for(const auto&xi : *x) { ret+=(float)std::pow(xi,n); }
    return ret/(x->size());
}
/// @brief start moment: v(x)=m(x^n)
inline float v(std::vector<interval> * x, uint n){
    float ret = 0.0f; 
    float sumn = 0.0f;
    for(const interval&xi : *x) {sumn+=xi.n; ret+=xi.n*(float)std::pow(0.5f*(xi.max+xi.min),n); }
    return ret/(sumn);
}


/// @brief central moment: u(x)=m(x-x^n)
inline float u(std::vector<float> * x, uint n){
    float ret = 0.0f; float mx = m(x);
    for(const auto&xi : *x) { ret+=(float)powf(xi-mx,n); }
    return ret/(x->size()-1);
}
/// @brief central moment: u(x)=m(x-x^n)
inline float u(std::vector<interval> * x, uint n){
    float ret = 0.0f; float mx = m(x);
    float sumn = 0.0f;
    for(const interval&xi : *x) {sumn+=xi.n; ret+=xi.n*(float)powf(0.5f*(xi.max+xi.min)-mx,n); }
    return ret/(sumn);
}

/// @brief central moment: u(x)=m(x-x^n)
inline float d(std::vector<float> * x){return u(x,2);}
/// @brief central moment: u(x)=m(x-x^n)
inline float d(std::vector<interval> * x){return u(x,2);}


/// @brief sqrt(d(x))
inline float sigma(std::vector<float> * x){return sqrtf(u(x,2));}
/// @brief sqrt(d(x))
inline float sigma(std::vector<interval> * x){return sqrtf(u(x,2));}

/// @brief excess(x)
inline float excess(std::vector<float> * x){
    float dis = d(x);
    return u(x,4)/(dis*dis) - 3.0f;
}
/// @brief excess(x)
inline float excess(std::vector<interval> * x){
    float dis = d(x);
    return u(x,4)/(dis*dis) - 3.0f;
}

/// @brief asym(x)
inline float asym(std::vector<float> * x){
    float dis = d(x);
    return u(x,3)/(dis*dis*dis);
}
/// @brief asym(x)
inline float asym(std::vector<interval> * x){
    float dis = d(x);
    return u(x,3)/(dis*dis*dis);
}


void sort(std::vector<float> * x){
    for(uint i = 0; i < x->size(); ++i){
        /*assume first elem is lowest*/
        float lowest = x->at(i);
        uint lowest_ind = i;
        /*try find lower elem*/
        for(uint j = i+1; j < x->size(); ++j){
            if(x->at(j)<lowest){
                lowest = x->at(j);
                lowest_ind = j;
            }
        }
        /*swap x[i] with x lowest*/
        if(lowest_ind != i){ x->at(lowest_ind) = x->at(i); x->at(i)=lowest; }

    }
}


std::vector<float> sorted_copy(std::vector<float> x){
    std::vector<float> ret = x;   
    for(uint i = 0; i < ret.size(); ++i){
        /*assume first elem is lowest*/
        float lowest = ret.at(i);
        uint lowest_ind = i;
        /*try find lower elem*/
        for(uint j = i+1; j < ret.size(); ++j){
            if(ret.at(j)<lowest){
                lowest = ret.at(j);
                lowest_ind = j;
            }
        }
        /*swap x[i] with x lowest*/
        if(lowest_ind != i){ ret.at(lowest_ind) = ret.at(i); ret.at(i)=lowest; }
    }
    return ret;
}


inline float med(std::vector<float> * x){
    std::vector<float> x_sort = sorted_copy(*x);
    uint len2floor = floor(x_sort.size()/2);
    if(x_sort.size() % 2 == 0){
        return 0.5f*(x_sort.at(len2floor) + x_sort.at(len2floor+1));
    }
    return  (float)(x_sort.size() % 2 == 0) * 0.5f * (x_sort.at(len2floor)+x_sort.at(len2floor+1)) + //find the everage between 2 elements
    (float)(x_sort.size() % 2 != 0) * x_sort.at(len2floor); // or just return single element
}


/// @brief difference between last and first elem
inline float ordered_range(std::vector<float> * x){
    return (x->at(x->size()-1)-x->at(0));
}

/// @brief Sturges' rule: n=1+log2(n) 
inline uint optimal_sturges(uint len){
    return 1+std::floor(std::log2f(len));
}


/* INTERVAL TABLE */
std::vector<interval> * get_intervals(std::vector<float> * x){
    std::vector<interval> * ret = new std::vector<interval>();
    float min_ = min(x), max_=max(x);
    
    uint n = optimal_sturges(x->size());
    float step = (max_ - min_)/(float)n;
    float first = min_; float last = min_+step;


    while (first <= max_)
    {
        if(last + step >= max_) last += 0.01f;
        interval iv(first,last,0);
        for(uint i = 0; i < x->size(); ++i){
            if((x->at(i) >= first && x->at(i) < last)) iv.n++;
        }
        if(iv.n > 0) ret->push_back(iv);
        first = last; last = first+step;
    }
    return ret;
}
inline std::string interval_to_str(interval i){
    return ""+std::to_string(i.min)+"\t"+std::to_string(i.max)+"\t"+std::to_string(i.n);
}

inline std::string interval_to_html(interval i, uint n, float N=1){
    return "<tr><td>"+std::to_string(n)+"</td><td>"+std::to_string(i.min)+"</td><td>"+std::to_string(i.max)+"</td><td>"+std::to_string(i.n)+"</td><td>"+std::to_string((float)i.n/N)+"</td><td>""</td></tr>";
}

std::string interval_table(std::vector<interval> * intervals){
    std::string ret = "<table><tr><td>#</td><td>min</td><td>max</td><td>n</td><td>p</td></tr>";
    int num = 0;
    float N=0;
    for(const auto & i : *intervals){
        N+= i.n;
    }

    uint sum = 0; float sump;
    for(const auto & i : *intervals){
        ++num;
        ret += interval_to_html(i,num,N);
        sum+=i.n;
        
    }
    ret += "</table>";
    return ret;
}




//____________________________________________________________


namespace normald{


inline float d(float x, float mx=0.0f, float sx=1.0f) {
    return (1.0f / (sx * std::sqrt(2.0f*M_PI))) * std::exp(-0.5f * std::pow((x-mx)/sx,2.0f));
}
/// @brief 
/// @param params[0] mx
/// @param params[1] sx 
inline float d(float x, float * params) {
    return (1.0f / (params[1] * std::sqrt(2.0f*M_PI))) * std::exp(-0.5f * std::pow((x-params[0])/params[1],2.0f));
}



float cd(float x,float mx=0.0f, float sx=1.0f, float eps = 0.001f){
    if(x>=mx){
        float sum = 0.5f;
        for(float xi = mx; xi < x;xi+=eps){
            sum+=(0.5f*(d(xi-eps,mx,sx)+d(xi+eps,mx,sx))*eps);
        }
        return sum;
    }
    else{
        return 1.0f - cd(2.0f*mx-x,mx,sx,eps);
    }
}

float cdint(float x0,float x1,float mx=0.0f, float sx=1.0f, float eps = 0.001f){
    float sum = 0.0f;
    //if(x >= mx+sx*5.0f) return 1.0f;
    for(float xi = x0; xi <= x1;xi+=eps){
        sum+=(0.5f*(d(xi-eps,mx,sx)+d(xi+eps,mx,sx))*eps);
    }
    return sum;
}


double cdr_approx(const double& p, double mx=0.0f, double sx=1.0f) {
    static double a[4] = {   2.50662823884,
                         -18.61500062529,
                          41.39119773534,
                         -25.44106049637};

    static double b[4] = {  -8.47351093090,
                          23.08336743743,
                         -21.06224101826,
                           3.13082909833};

    static double c[9] = {0.3374754822726147,
                        0.9761690190917186,
                        0.1607979714918209,
                        0.0276438810333863,
                        0.0038405729373609,
                        0.0003951896511919,
                        0.0000321767881768,
                        0.0000002888167364,
                        0.0000003960315187};

     if (p >= 0.5 && p <= 0.92) {
    double num = 0.0;
    double denom = 1.0;

    for (int i=0; i<4; i++) {
      num += a[i] * pow((p - 0.5), 2*i + 1);
      denom += b[i] * pow((p - 0.5), 2*i);
    }
    return num/denom;

  } else if (p > 0.92 && p < 1) {
    double num = 0.0;

    for (int i=0; i<9; i++) {
          num += c[i] * pow((log(-log(1-p))), i);
    }
    return mx+num*sx;

  } else {
        return -1.0*cdr_approx(1-p);
  }
}


/// @brief 
/// @param params[0] mx
/// @param params[1] sx 
float cd(float x, float * params, float eps = 0.001f){
    return cd(x,params[0],params[1],eps);
}

float cdr(float p, float mx=0.0f, float sx=1.0f, float eps = 0.0001f){
    if(p < 0.5f) return -cdr(1.0f-p); //return symmetrical value
    float x = mx, integral = 0.5f; // assume that if (p >= 0.5) then x >= 0

    while (true)
    {
        integral += 0.5f*eps*(d(x,mx,sx)+d(x+eps,mx,sx));
        x += eps;
        if(fabs(integral-p) < eps) return x;
    }   
}

}

//____________________________________________________________



namespace lognormald{


inline float d(float x, float mx=0.0f, float sx=1.0f) {
    if(x<=0.0f) return 0.0f;
    return (1.0f / (x*sx * std::sqrt(2.0f*M_PI))) * std::exp( std::pow((std::log(x)-mx),2.0f) / (2.0f*sx*sx) );
}


float cd(float x,float mx=0.0f, float sx=1.0f, float eps = 0.001f){
    if(x>=mx){
        float sum = 0.5f;
        for(float xi = mx; xi < x;xi+=eps){
            sum+=(0.5f*(d(xi-eps,mx,sx)+d(xi+eps,mx,sx))*eps);
        }
        return sum;
    }
    else{
        return 1.0f - cd(2.0f*mx-x,mx,sx,eps);
    }
}

float cdint(float x0,float x1,float mx=0.0f, float sx=1.0f, float eps = 0.001f){
    float sum = 0.0f;
    for(float xi = x0; xi <= x1;xi+=eps){
        sum+=(0.5f*(d(xi-eps,mx,sx)+d(xi+eps,mx,sx))*eps);
    }
    return sum;
}

}



namespace t{

    inline float d(float x, int v=1) {
        return expf(lgammaf((v + 1) / 2.0) 
                         - lgammaf(v / 2.0) 
                         - 0.5 * log(v * M_PI) 
                         - (v + 1) / 2.0 * log(1 + x * x / v));
    }

    float cd(float x, float v=1, float epsmul = 0.01f) {
        if(x > 0.0f){
            return 1.0f - cd(-x,v,epsmul);
        }
        int iter = 0;
        float sum = 0.0f;
        float xn = -1000.0f/v;
        float x1;

        float eps=4.0f*fabs(xn)*epsmul,tdf0,tdf1=d(xn+eps,v);
        while(xn<x){
            ++iter;
            tdf0=tdf1;
            tdf1 = d(xn+eps,v);
            
            sum += (tdf1+tdf0)*0.5f*eps;
            xn+=eps;
            eps = 4.0f*fabs(xn)*epsmul+epsmul;
        }
        std::cout << "\titer:" << iter << "\t";
        return sum;
    }

    inline float cdr_apr(float x, float n=1, float mx=0.0f, float sx=1.0f) {
        float z = normald::cdr_approx(x,mx,sx);
        return z*(1.0f + 
        (z*z+1.0f)/(4.0f*n) + 
        ((5.0f*z*z*z*z+16.0f*z*z+3.0f)/(96.0f*n*n)) + 
        ((3*std::pow(z,6)+19.0f*std::pow(z,4)+17.0f*std::pow(z,2)-15.0f)/(384.0f*n*n*n)));

        /*float z2 = z*z;
        float term1 = (z2 + 2.0f) / (4.0f * v);
        float term2 = (z2 + 2.0f) * (z2 + 3.0f) / (96.0f * v*v);
        
        float t_p = z * (1.0f + term1 + term2);
        
        return t_p;*/

        /*return (z*(1.0f+0.25f*(powf(z,2)+1.0f)/v)) /
                (1.0f-(powf(z,2)+2.0f)/v);*/
    }


    float cdr(float p, float v = 1.0f, float eps = 0.0001f){
        if(p < 0.5f) return -cdr(1.0f-p,v,eps); //return symmetrical value
        float x = 0.0f, integral = 0.5f; // assume that if (p >= 0.5) then x >= 0

        while (true)
        {
            integral += 0.5f*eps*(d(x,v)+d(x+eps,v));
            x += eps;
            if(p-integral < eps) return x;
            
        }   
    }

}



namespace chi2{

    inline float d(float x, float v=1.0f) {
        if(x<=0 || v < 1.0f) return -0.0f;
        float res = (v * 0.5f) * logf(0.5f) + (v * 0.5f - 1.0f) * logf(x) - (x * 0.5f) - lgammaf(v * 0.5f);
        return expf(res);
    }

    float cd(float x, float v=1, float epsmul = 0.01f) {
        if(x < 0.0f){
            return 0.0f;
        }
        int iter = 0;
        float sum = 0.0f;
        float xn = 0.0f;
        float x1;

        float eps=4.0f*epsmul,tdf0,tdf1=d(xn+eps,v);
        while(xn<x){
            ++iter;
            tdf0=tdf1;
            tdf1 = d(xn+eps,v);
            
            sum += (tdf1+tdf0)*0.5f*eps;
            xn+=eps;
            eps = 4.0f*fabs(tdf0)*epsmul+epsmul;
        }
        std::cout << "\titer:" << iter << "\t";
        return sum;
    }

    inline float cdr_approx(float p, float v=1, float mx=0.0f, float sx=1.0f) {
        float a = -1.37266f, b = 1.06807f, c = 2.13161f, d = -0.04589f;
        return powf(
            a + b*sqrtf(v) + (c+d*sqrtf(v)) * sqrtf(-log10f(p))
        ,2.0f);

    }

    float cdr(float p, float v = 1.0f, float _eps = 0.01f){
        float eps = _eps;
        //std::cout << "\nchi.cdr p=" << p << std::endl; 
        if(p <= 0.0f) return 0.0f;
        float x = 0.0f, integral = 0.0f, d0,d1; // assume that if (p >= 0.0) then x >= 0
        if(v>=100.0f && p>0.5f){
            integral=0.5f;
            x=v-0.66596;
        }
        int iter = 0; 
        while (true)
        {
            ++iter;
            d0 = d(x,v);
            d1 = d(x+eps,v);
            integral += 0.5f*eps*(d0+d1);
            x += eps;
            if(p-integral <= eps) { return x;}
            eps =_eps;
            //std::cout << " e="<<eps;
        }   
        
    }
}


namespace conf_intervals{
    /*this struct is being used for all types of confidence intervals*/
    struct ci
    {
        public:
        float value0;
        float value1;
        float eps;
        float quantile;
        float alpha;
        std::string desc;

        ci(float value0,float value1,float quantile,float beta,std::string desc){
            this->value0=value0;
            this->value1=value1;
            this->quantile=quantile;
            this->alpha=beta;
            this->desc = desc;
        };

        std::string to_string(){
            return desc + ": (" + std::to_string(value0) + ";" + std::to_string(value1) + ") quantile:" +  std::to_string(quantile) + " a=" + std::to_string(alpha);
        }
    };

    /// @brief  
    /// @param b beta 
    ci mx_rough_method(std::vector<float>* v,  float b){
        float ub,eps;
        uint n = v->size();
        float mx = m(v);
        float sx = sigma(v);
        ub = normald::cdr_approx(0.5f+b*0.5f);

        eps = ub * sx/sqrtf(n);

        return ci(mx-eps,mx+eps,ub,b,"математичне очікування [грубий]");
    }

    ci mx_exact_method(std::vector<float>* v, float b){
        float ub,eps;
        float mx = m(v);
        float sx = sigma(v);
        uint n = v->size();
        ub = t::cdr_apr(0.5f+b*0.5f,n-1);

        eps = ub * sx/sqrtf(n);

        return ci(mx-eps,mx+eps,ub,b,"математичне очікування [точний]");
    }

    enum ci_dx_mode {plain, normal};
    ci dx_rough_method(std::vector<float>* v, float b, ci_dx_mode mode){
        float ub,eps;
        float sx = sigma(v);
        uint n = v->size();
        std::string label = "дисперсія [грубий]";
        ub = normald::cdr(0.5f+b*0.5f,m(v),sx);

        float sDx;
        if(mode == ci_dx_mode::plain) {sDx = sqrtf((0.8f*n+1.2f)/(n*(n-1.0f)))*sx*sx; label+="[рівномірний]";}
        else if(mode == ci_dx_mode::normal) {sDx = sqrtf(2.0f/(n-1.0f))*sx*sx;label+="[нормльний]";}

        eps = ub * sDx;

        
        return ci(sx*sx-eps,sx*sx+eps,ub,b,label);
    }

    ci dx_exact_method(std::vector<float>* v, float b, float eps = 0.01f){
        float a = 1.0f - b;
        float sx = sigma(v);
        uint n = v->size();
        std::string label = "дисперсія [точний]";

        float chi0 = chi2::cdr(1-a*0.5f,n-1.0f,eps);
        float chi1 = chi2::cdr(a*0.5f,n-1.0f,eps);
        //std::cout << "\n\ta=" << (1.0f-a*0.5f) << "\t quantile:" << chi2::cdr((1.0f-a)*0.5f,n-1) << std::endl;
        //std::cout << "\ta=" << a*0.5f << "\t quantile:" << chi2::cdr(a*0.5f,optimal_sturges(n-1)) << std::endl;
        
        return ci(n*sx*sx/chi0,n*sx*sx/chi1,chi0,b,label);
    }
}



namespace pirson{

struct pirson_table
{
    std::vector<float> chi_criteria;
    std::vector<float> p_theoretical;
    float chi_sum = 0.0f;
    float chi_crit = 0.0f;
};

pirson_table normal_p(std::vector<interval>* inters, float alpha = 0.95f){
    pirson_table ret;
    float diff, mx = m(inters), sx = sigma(inters);
    float sum = 0.0f, np = 0.0f;
    float n = countn(inters);
    for(int i=0; i<inters->size(); ++i){
        interval inter = inters->at(i);
        float theor = normald::cdint(inter.min,inter.max,mx,sx) * n;
        diff = powf(theor-inter.n, 2.0f)/theor;
        ret.chi_criteria.push_back(diff);
        ret.p_theoretical.push_back(theor/n);
        sum += diff;
    }
    ret.chi_sum = sum;
    ret.chi_crit = chi2::cdr(alpha,inters->size()-1);
    return ret;
};
pirson_table lognormal_p(std::vector<interval>* inters, float alpha = 0.95f){
    pirson_table ret;
    float diff, mx = m(inters), sx = sigma(inters);
    float sum = 0.0f, np = 0.0f;
    float n = countn(inters);
    for(int i=0; i<inters->size(); ++i){
        interval inter = inters->at(i);
        float theor = lognormald::cdint(inter.min,inter.max,mx,sx) * n;
        diff = powf(theor-inter.n, 2.0f)/theor;
        ret.chi_criteria.push_back(diff);
        ret.p_theoretical.push_back(theor/n);
        sum += diff;
    }
    ret.chi_sum = sum;
    ret.chi_crit = chi2::cdr(alpha,inters->size()-1);
    return ret;
};
pirson_table plain_p(std::vector<interval>* inters, float alpha = 0.95f){
    pirson_table ret;
    float diff;
    float sum = 0.0f, np = 0.0f; float a = inters->at(0).min, b = inters->at(inters->size()-1).max;
    float n = countn(inters);
    float theor = (float)n/inters->size();
    for(int i=0; i<inters->size(); ++i){
        interval inter = inters->at(i);
        float p0 = (inter.min>=a)*(inter.min-a)/(b-a), p1 = (inter.max<=b)*(inter.max-a)/(b-a) + (inter.max>b);
        diff = powf(theor-inter.n, 2.0f)/theor;
        ret.chi_criteria.push_back(diff);
        ret.p_theoretical.push_back(theor/n);
        sum += diff;
    }
    ret.chi_sum = sum;
    ret.chi_crit = chi2::cdr(alpha,inters->size()-1);
    return ret;
};



inline std::string _pirson_interval_to_html(interval i, uint n, float N=1, float theoretical = 1.0f, float chi = 0.0f){
    return "<tr><td>"+std::to_string(n)+"</td><td>"+std::to_string(i.min)+"</td><td>"+std::to_string(i.max)+"</td><td>"+std::to_string((int)i.n)+"</td><td>"+std::to_string((float)i.n/N)+"</td><td>"+std::to_string(theoretical)+"</td><td>"+std::to_string(N*theoretical)+"</td><td>"+std::to_string(chi)+"</td><td>""</td></tr>";
}

std::string interval_table(std::vector<interval> * intervals, pirson_table pt){
    std::string ret = "<table><tr><td>#</td><td>min</td><td>max</td><td>mi</td><td>mi/N</td><td>pi</td><td>n*pi</td><td>chi2-criteria</td></tr>";
    int num = 0;
    float N=0;
    for(const auto & i : *intervals){
        N+= i.n;
    }

    uint sum = 0; float sump;
    for(const auto & i : *intervals){
        ++num;
        ret += _pirson_interval_to_html(i,num,N,pt.p_theoretical.at(num-1),pt.chi_criteria.at(num-1));
        sum+=i.n;
        
    }
    ret += "</table>";
    return ret;
}


}


namespace kolmohorov{

struct kolmohorov_table
{
    std::vector<float> diffs;
    std::vector<float> theoretical;
    std::vector<float> p_cumulative;
    float max_diff = 0.0f;
    float lambda = 0.0f;
    float lambda_crit = 0.0f;
};


kolmohorov_table normal_k(std::vector<interval>* inters, float alpha = 0.95f){
    kolmohorov_table ret;
    float diff, mx = m(inters), sx = sigma(inters), n=countn(inters);
    float sum = 0.0f;

    for(int i=0; i<inters->size(); ++i){
        interval inter = inters->at(i);
        sum+=inter.n/n;
        float theor = normald::cd(inter.max,mx,sx,0.001f);
        diff = fabs(sum-theor);
        ret.p_cumulative.push_back(sum); ret.diffs.push_back(diff); ret.theoretical.push_back(theor);
    }
    ret.max_diff = max(&ret.diffs);
    ret.lambda = ret.max_diff*sqrtf(n);
    return ret;
}


kolmohorov_table lognormal_k(std::vector<interval>* inters, float alpha = 0.95f){
    kolmohorov_table ret;
    float diff, mx = m(inters), sx = sigma(inters), n=countn(inters);
    float sum = 0.0f;

    for(int i=0; i<inters->size(); ++i){
        interval inter = inters->at(i);
        sum+=inter.n/n;
        float theor = lognormald::cd(inter.max,mx,sx,0.001f);
        diff = fabs(sum-theor);
        ret.p_cumulative.push_back(sum); ret.diffs.push_back(diff); ret.theoretical.push_back(theor);
    }
    ret.max_diff = max(&ret.diffs);
    ret.lambda = ret.max_diff*sqrtf(n);
    return ret;
}


kolmohorov_table plain_k(std::vector<interval>* inters, float alpha = 0.95f){
    kolmohorov_table ret;
    float diff, mx = m(inters), sx = sigma(inters), n=countn(inters);
    float sum = 0.0f;
    float a = inters->at(0).min, b = inters->at(inters->size()-1).max;

    for(int i=0; i<inters->size(); ++i){
        interval inter = inters->at(i);
        sum+=inter.n/n;
        float p0 = (inter.min>=a)*(inter.min-a)/(b-a), p1 = (inter.max<=b)*(inter.max-a)/(b-a) + (inter.max>b);
        float theor = (p1-p0);
        diff = fabs(sum-theor);
        ret.p_cumulative.push_back(sum); ret.diffs.push_back(diff); ret.theoretical.push_back(theor);
    }
    ret.max_diff = max(&ret.diffs);
    ret.lambda = ret.max_diff*sqrtf(n);
    return ret;
}


}


/// @brief only ca
/// @param vec use only with sorted values 
/// @return vector of lambda_i = |x_i-x_i-1|/sx
std::vector<float> irvin_criterium(std::vector<float> * vec){
    float sx = sigma(vec);
    std::vector<float> lambda;
    for(uint i = 1; i < vec->size(); ++i){
        lambda.push_back( fabs((vec->at(i)-vec->at(i-1))) / sx );
    }
    return lambda;
}



struct series_crit
{
public:
    std::vector<bool> signs; // + or -
    float me = 1.1f; // median(x)
    uint series_num = 1; //N if there at least one elem, there's already N=1 
    uint n1 = 0; //n1
    uint n2 = 0; //n2
    float z_criteria = 0.0f;

    series_crit(std::vector<float> * x){
        me = med(x);
        bool sign = x->at(0) > me;
        for(uint i = 0; i < x->size(); ++i){
            series_num += (uint)((x->at(i) > me) != sign); //current sign not equal to previous one! that's new serie
            sign = x->at(i) > me;
            signs.push_back(sign);
            n1 += (uint)(sign); n2 += (uint)(!sign);
        }
        float z_up = series_num - 2.0f*n1*n2/(float)(n1+n2) - 1.5f;
        float z_down = sqrtf( 2.0f*n1*n2*(2.0f*n1*n2-(n1+n2)) / (powf(n1+n2,2.0f)*(float)(n1+n2-1.0f)));
        z_criteria = z_up/z_down;
    }
};




std::vector<float> _rank_sorted(std::vector<float> * x){
    uint same_rank = 1;
    std::vector<float> ret;
    float prev_value = 0.0f;
    float expected_rank = 1.0f;
    for(uint i=0;i<x->size();++i){
        ret.push_back(expected_rank);
        if(x->at(i)==prev_value && i > 0){
            ++same_rank;
        }
        else if(same_rank>1 && i > 1) /*in that case, make average all previous ranks*/{
            float avg = 0.0f;
            for(int j=i-same_rank;j<=i-1;++j)avg+=ret.at(j);
            avg/=same_rank;
            for(int j=i-same_rank;j<=i-1;++j){ret.at(j) = avg; /*make all of them average*/}
            same_rank = 1; /*same elements ended, so drop value*/
        }
        ++expected_rank;/*simple increment of the rank every time, despite equal values*/
        prev_value = x->at(i);
    }
    return ret;
}



std::vector<float> rank_avg(std::vector<float> * x){
    std::vector<float> ret;
    std::vector<float> sorted_x = sorted_copy(*x);
    std::vector<float> ranks_sorted_x = _rank_sorted(&sorted_x);

    for(uint i = 0; i < x->size(); ++i){
        int ind = findi(&sorted_x,x->at(i));
        if(ind >= 0) ret.push_back(ranks_sorted_x.at(ind));
    }

    return ret;
}


struct wilcoxon{
    std::vector<float>*x; std::vector<float>*y; float a;
    std::vector<float> d;
    std::vector<float> ranks;
    float r0,r1,mt,st,u_crit,z;

public:
    wilcoxon(std::vector<float> * x, std::vector<float> * y, float a=0.05f){
        this->x = x; this->y = y; this->a = a; //input
        d = sub_vec(x,y);
        std::vector<float> module_d;
        for(uint i=0;i<d.size();++i) module_d.push_back(fabs(d.at(i)));
        ranks = rank_avg(&module_d);
        r0=0.0f, r1=0.0f;
        for(uint i=0;i<d.size();++i){
            r0 += (ranks.at(i)) * (d.at(i)<0);
            r1 += (ranks.at(i)) * (d.at(i)>0);
        }
        float n = x->size();
        mt = n*(n+1.0f)*0.25f, st = sqrtf(n*(n+1.0f)*(2.0f*n+1.0f)/24.0f);
        u_crit = normald::cdr_approx(1.0f-a*0.5f); 
        z = fabs((r1-mt)/st);
    }
};



struct wilcoxon_mann_u{
    std::vector<float>*x; std::vector<float>*y; float a; //input
    std::vector<float> xy_sorted;
    std::vector<float> xy_ranks = rank_avg(&xy_sorted);
    float r0,r1,n0,n1,w0,w1,w,u_crit,z;
public:
    wilcoxon_mann_u(std::vector<float> * x, std::vector<float> * y, float a=0.05f){
        this->x = x; this->y = y; this->a = a;
        xy_sorted = *x; 
        for(uint i=0; i<y->size();++i) xy_sorted.push_back(y->at(i));
        sort(&xy_sorted);
        xy_ranks = rank_avg(&xy_sorted);
        r0=0.0f,r1=0.0f;
        for(uint i=0; i<x->size();++i){int ind = findi(&xy_sorted,x->at(i)); if(ind>=0) r0+=xy_ranks.at(ind);}
        for(uint i=0; i<y->size();++i){int ind = findi(&xy_sorted,y->at(i)); if(ind>=0) r1+=xy_ranks.at(ind);}

        n0 = x->size(), n1 = y->size();
        w0 = n0*n1 + 0.5f*(n0*(n0+1.0f)) - r0; 
        w1 = n0*n1 + 0.5f*(n1*(n1+1.0f)) - r1;
        w = minf(w0,w1);
        u_crit = normald::cdr_approx(1.0f-a*0.5f);
        z = fabs((w-0.5f*n0*n1) / sqrtf(n0*n1*(n0+n1+1.0f)/12.0f));
        //return u0*(float)(u0<=u1) + u1*(float)(u1<u0);
    }

};




struct spearman
{
    std::vector<float>*x; std::vector<float>*y; float a; //input
    std::vector<float>r; std::vector<float> s;
    float p=0.0f,tr=0.0f,ts=0.0f,t,crit,square_d_sum=0.0f;
public:
    spearman(std::vector<float>*x, std::vector<float>*y, float a){
        this->x = x; this->y = y; this->a = a;
        r = rank_avg(x); s = rank_avg(y);
        for(uint i=0;i<r.size();++i){
            std::vector<float> used_values;
            if(findi(&used_values,r.at(i))==-1){
                uint trg_=countsame(&r,r.at(i));
                used_values.push_back(r.at(i));
                tr+=(powf(trg_,3.0f)-trg_)/12.0f;
            }

            square_d_sum += powf(r.at(i)-s.at(i),2.0f);
        }

        for(uint i=0;i<s.size();++i){
            std::vector<float> used_values;
            if(findi(&used_values,s.at(i))==-1){
                uint tsg_=countsame(&s,s.at(i));
                used_values.push_back(s.at(i));
                ts+=(powf(tsg_,3.0f)-tsg_)/12.0f;
            }
        }
        float n = s.size();
        p = 1.0f - 6.0f*square_d_sum/(n*(n*n-1.0f)-(tr+ts));
        t=p*sqrtf(n-2.0f)/sqrtf(1.0f-p*p);
        crit = t::cdr(1-a*0.5,n-2.0f);
    }
};




struct cendal
{
    std::vector<float>*x; std::vector<float>*y; float a; //input
    std::vector<float>r; std::vector<float> s;
    std::vector<float>r_sorted; std::vector<float> s_sorted;
    float p=0.0f,q=0.0f,pq=0.0f,tr=0.0f,ts=0.0f,t,crit,square_d_sum=0.0f,tau=0.0f;
public:
    cendal(std::vector<float>*x, std::vector<float>*y, float a){
        this->x = x; this->y = y; this->a = a;
        r = rank_avg(x); s = rank_avg(y);
        r_sorted = sorted_copy(r);
        for(uint i=0;i<r.size();++i){
            int ind = findi(&r,r_sorted.at(i));
            s_sorted.push_back(s.at(ind)); //depends on r position after sorting
        }
        for(uint i=0;i<s.size();++i){
            for(uint j=i+1;j<s.size();++j){
                p+=(s.at(i)>s.at(j));
                q+=(s.at(i)<s.at(j));
            }
        }
        pq = fabs(p-q);

        
        for(uint i=0;i<r.size();++i){
            std::vector<float> used_values;
            if(findi(&used_values,r.at(i))==-1){
                uint trg_=countsame(&r,r.at(i));
                used_values.push_back(r.at(i));
                tr+=trg_*(trg_-1.0f)/2.0f;
            }

            square_d_sum += powf(r.at(i)-s.at(i),2.0f);
        }

        for(uint i=0;i<s.size();++i){
            std::vector<float> used_values;
            if(findi(&used_values,s.at(i))==-1){
                uint tsg_=countsame(&s,s.at(i));
                used_values.push_back(s.at(i));
                ts+=tsg_*(tsg_-1.0f)/2.0f;
            }
        }

        float n = s.size();
        tau= pq / sqrtf( (0.5f*n*(n-1.0f)-tr)*(0.5f*n*(n-1.0f)-ts) );
        t = tau*sqrtf( (9.0f*(n-1.0f)) / (2.0f*(2.0f*n+5.0f)) );

        crit = t::cdr(1-a*0.5,n-2.0f);
    }
};




}