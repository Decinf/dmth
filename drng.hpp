#pragma once
#include <cmath>
#include <iostream>
#include <bitset>

namespace dmth{



float mid_sqr(long r0){
    long tmp = 1, n=0;
    while (tmp <= r0){tmp*=10;++n;} /*count n digits in r0*/
    //std::cout << "n=" << n << std::endl;
    long r0_ = r0;
    if(n%2 != 0) r0_*=10;
    long r1 = r0_*r0_;

    /*if n = 4*/
    /*example: 1234 -> 1522756 -> __3456__*/
    long div = pow(10,n/2);
    long mod = pow(10,n/2+n);

    return (r1%mod)/div;

}


uint lfsr(uint32_t x0, uint32_t x1){
    /*1010 << 0101*/
    uint32_t x00 = x0 & !(UINT32_MAX>>1);
    uint32_t x10 = x1 & !(UINT32_MAX>>1);
    x10^=x00;
    std::cout << "x0=" << std::bitset<sizeof(uint32_t)>(x0) << "\tx1=" << std::bitset<sizeof(uint32_t)>(x1) << "\tOR\txored first=" << std::bitset<sizeof(uint32_t)>(x00) << "\t=\t" << std::bitset<sizeof(uint32_t)>((x0 << 1)|x10) << std::endl;
    //ushort x10 = x1 & 0xFF;

    return ((x0 << 1)|x10);
}

int32_t _myran_seed = 54621;
int32_t myran(int32_t seed){
    uint32_t ret = ((1273+((seed+1)*574685)%RAND_MAX) >> 17) xor (((seed-1)*436347)%RAND_MAX << 13);
    _myran_seed = ret;
    return ret;
}

int32_t myran(){
    //++_myran_seed;
    uint32_t ret = ((1273+(_myran_seed*517466385)%RAND_MAX) >> 17) xor ((_myran_seed*436263647)%RAND_MAX << 13);
    _myran_seed = ret;
    return ret;
}

float myranf(){
    return ((float)myran()/(RAND_MAX))*0.5f + 0.5f;
}

float myranf(float min, float max){
    return myranf()*(max-min) + min;
}


}