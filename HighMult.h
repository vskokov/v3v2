#include <fstream>
#include <blitz/blitz.h>
#include <blitz/array.h>
#include <vector>
#include <complex>

#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"
#include "boost/random/normal_distribution.hpp"
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/math/special_functions/factorials.hpp>


using namespace std;
using namespace blitz;

typedef complex<double> cd;
typedef boost::mt19937 RNGType;

unsigned long mix(unsigned long a, unsigned long b, unsigned long c) {
    a=a-b;
    a=a-c;
    a=a^(c >> 13);
    b=b-c;
    b=b-a;
    b=b^(a << 8);
    c=c-a;
    c=c-b;
    c=c^(b >> 13);
    a=a-b;
    a=a-c;
    a=a^(c >> 12);
    b=b-c;
    b=b-a;
    b=b^(a << 16);
    c=c-a;
    c=c-b;
    c=c^(b >> 5);
    a=a-b;
    a=a-c;
    a=a^(c >> 3);
    b=b-c;
    b=b-a;
    b=b^(a << 10);
    c=c-a;
    c=c-b;
    c=c^(b >> 15);
    return c;
}

unsigned long seed = mix(clock(), time(NULL), getpid());

RNGType rng(seed);

#include <iostream>
using namespace std;


