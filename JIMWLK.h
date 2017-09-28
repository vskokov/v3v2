#include <fstream>
#include <fftw3.h>
#include <blitz/blitz.h>
#include <blitz/array.h>
#include "fftw_blitz.h"
#include <vector>
#include <complex>

#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"
#include "boost/random/normal_distribution.hpp"
#include  <boost/random/uniform_real_distribution.hpp>
#include <boost/math/special_functions/factorials.hpp>

using namespace std;
using namespace blitz;

typedef complex<double> cd;


typedef boost::mt19937 RNGType;

unsigned long mix(unsigned long a, unsigned long b, unsigned long c)
{
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
//unsigned long seed = 0.0;




RNGType rng(seed);
//RNGType rng(2311);




#include <iostream>
using namespace std;


template<class T, int N, int M>
inline TinyMatrix<T,M,N> transpose(const TinyMatrix<T,N,M>& x)
{
    TinyMatrix<T,M,N> m;
    for (int i=0; i < N; ++i)
        for (int j=0; j < M; ++j)
            m(i,j) = x(j,i);
    return m;
}




double trace (TinyMatrix <double,8,8>  matrix)
{
    double sum;
    sum=0;
    for (int i=0; i<8; i++) sum+=matrix(i,i);
    return sum;
}


cd trace (TinyMatrix <cd,3,3>  matrix)
{
    cd sum;
    sum=cd(0,0);
    for (int i=0; i<3; i++) sum+=matrix(i,i);
    return sum;
}

blitz::TinyMatrix<cd,3,3> Product(TinyMatrix<cd,3,3>   A, TinyMatrix<cd,3,3>   B)
{
    firstIndex i;
    secondIndex j;
    thirdIndex k;
    blitz::TinyMatrix<cd,3,3> C;
    C = sum(A(i,k) * B(k,j), k);
    return C;
}




blitz::TinyMatrix<double,8,8> Product(TinyMatrix<double,8,8>   A, TinyMatrix<double,8,8>   B)
{
    firstIndex i;
    secondIndex j;
    thirdIndex k;
    blitz::TinyMatrix<double,8,8> C;
    C = sum(A(i,k) * B(k,j), k);
    return C;
}


blitz::TinyMatrix<cd,3,3> matrix_exp(TinyMatrix<cd,3,3>&  M)
{
    blitz::TinyMatrix<cd,3,3>  sum;
    sum = cd(1.0,0.0),cd(0.0,0.0),cd(0.0,0.0),
    cd(0.0,0.0),cd(1.0,0.0),cd(0.0,0.0),
    cd(0.0,0.0),cd(0.0,0.0),cd(1.0,0.0);

    blitz::TinyMatrix<cd,3,3>  pro;
    pro=M;
    for(int i=1; i<31; i++)
    {
        sum = sum + pro / boost::math::factorial<double>(i);
        pro=Product(pro,M);
    }
    return sum;
}



blitz::TinyMatrix<cd,3,3> matrix_exp_evo(TinyMatrix<cd,3,3>&  M)
{
    blitz::TinyMatrix<cd,3,3>  sum;
    sum = cd(1.0,0.0),cd(0.0,0.0),cd(0.0,0.0),
    cd(0.0,0.0),cd(1.0,0.0),cd(0.0,0.0),
    cd(0.0,0.0),cd(0.0,0.0),cd(1.0,0.0);

    blitz::TinyMatrix<cd,3,3>  pro;
    pro=M;
    for(int i=1; i<25; i++)
    {
        sum = sum + pro / boost::math::factorial<double>(i);
        pro=Product(pro,M);
    }
    return sum;
}



#include "zheevv3.h"

blitz::TinyMatrix<cd,3,3> matrix_exp_an_antiH(TinyMatrix<cd,3,3>&  M)
{


    colorMat O;
    colorMat invO;
    blitz::TinyVector<double,3> w;

    zheevv3(M, O, w);

    colorMat out;
    colorMat D;
    D=
        exp(w(0)*cd(0.0,1.0)),cd(0.0,0.0),cd(0.0,0.0),
        cd(0.0,0.0),exp(w(1)*cd(0.0,1.0)),cd(0.0,0.0),
        cd(0.0,0.0),cd(0.0,0.0),exp(w(2)*cd(0.0,1.0));

    invO=blitz::conj(transpose(O));
    out = Product(O,Product(D,invO));

    return out;
}




inline double fac(int i)
{
    return boost::math::factorial<double>(i);
}



blitz::TinyMatrix<cd,3,3> dagger(TinyMatrix<cd,3,3>  M)
{
    blitz::TinyMatrix<cd,3,3> out;
    out = blitz::conj(transpose(M));
    return out;
}


blitz::TinyMatrix<double,8,8> dagger(TinyMatrix<double,8,8>  M)
{
    blitz::TinyMatrix<double,8,8> out;
    out = transpose(M);
    return out;
}





void FFTW(Array<complex<double>,2> &field, Array<complex<double>,2> &fftOfField)
{
    complex<double> *ptrField, *ptrFFT;
    fftw_plan plan;
    ptrField = field.data();
    ptrFFT = fftOfField.data();
    //plan = fftw_plan_dft_2d(field.rows(),field.cols(),reinterpret_cast<fftw_complex*>(ptrField),reinterpret_cast<fftw_complex*>(ptrFFT),FFTW_FORWARD,FFTW_ESTIMATE);
    plan = fftw_plan_dft_2d(field.rows(),field.cols(),reinterpret_cast<fftw_complex*>(ptrField),reinterpret_cast<fftw_complex*>(ptrFFT),FFTW_FORWARD,FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
}



void FFTW_b(Array<complex<double>,2> &field, Array<complex<double>,2> &fftOfField)
{
    complex<double> *ptrField, *ptrFFT;
    fftw_plan plan;
    ptrField = field.data();
    ptrFFT = fftOfField.data();
    plan = fftw_plan_dft_2d(field.rows(),field.cols(),reinterpret_cast<fftw_complex*>(ptrField),reinterpret_cast<fftw_complex*>(ptrFFT),FFTW_BACKWARD,FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
}









vector < TinyMatrix<cd,3,3>  > get_lambda()
{
    vector < TinyMatrix<cd,3,3>  > lambda(9);
    //for (int i =0; i<lambda.size(); i++) lambda.at(i).resize(3,3);

    lambda.at(0)=
        cd(0,0), cd(1,0), cd(0,0),
        cd(1,0), cd(0,0), cd(0,0),
        cd(0,0), cd(0,0), cd(0,0);

    lambda.at(1)=
        cd(0,0), cd(0,-1), cd(0,0),
        cd(0,1), cd(0,0), cd(0,0),
        cd(0,0), cd(0,0), cd(0,0);

    lambda.at(2)=
        cd(1,0), cd(0,0), cd(0,0),
        cd(0,0), cd(-1,0), cd(0,0),
        cd(0,0), cd(0,0), cd(0,0);

    lambda.at(3)=
        cd(0,0), cd(0,0), cd(1,0),
        cd(0,0), cd(0,0), cd(0,0),
        cd(1,0), cd(0,0), cd(0,0);

    lambda.at(4)=
        cd(0,0), cd(0,0), cd(0,-1),
        cd(0,0), cd(0,0), cd(0,0),
        cd(0,1), cd(0,0), cd(0,0);

    lambda.at(5)=
        cd(0,0), cd(0,0), cd(0,0),
        cd(0,0), cd(0,0), cd(1,0),
        cd(0,0), cd(1,0), cd(0,0);

    lambda.at(6)=
        cd(0,0), cd(0,0), cd(0,0),
        cd(0,0), cd(0,0), cd(0,-1),
        cd(0,0), cd(0,1), cd(0,0);

    lambda.at(7)=
        cd(1/sqrt(3),0), cd(0,0), cd(0,0),
        cd(0,0), cd(1/sqrt(3),0), cd(0,0),
        cd(0,0), cd(0,0), cd(-2/sqrt(3),0);


    lambda.at(8)=
        cd(1,0), cd(0,0), cd(0,0),
        cd(0,0), cd(1,0), cd(0,0),
        cd(0,0), cd(0,0), cd(1,0);

    lambda.at(8) *= cd(sqrt(2.0/3.0),0);

    for (int i=0; i<lambda.size(); i++) lambda.at(i)*=cd(0.5,0);

    return lambda;
}





vector < TinyMatrix<double,8,8>  > get_F()
{
    vector < TinyMatrix<double,8,8>  > F(8);

    vector < TinyMatrix<cd,3,3>  >  L = get_lambda();

    for (int i=0; i<8; i++) F.at(i)= 0.0;

    for (int a=0; a<8; a++)
        for (int b=0; b<8; b++)
            for (int c=0; c<8; c++)
            {
                blitz::TinyMatrix<cd,3,3> tmpab =  Product (L.at(a),L.at(b));
                blitz::TinyMatrix<cd,3,3> tmpba =  Product (L.at(b),L.at(a));
                blitz::TinyMatrix<cd,3,3> tmpd ;
                tmpd= tmpab-tmpba;
                F.at(a)(b,c) = real(complex<double>(0,-2.0)
                                    * trace( Product( tmpd, L.at(c)  ) ));

                //if(F.at(a)(b,c)!=0.0) cerr<< a << " " << b << " " << c << " " << F.at(a)(b,c) << endl;
            }



    return F;
}

