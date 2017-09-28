
#include "JIMWLK.h"
#include "string_conv.h"
#include <map>

using namespace std;
using namespace blitz;



//choice of parameters g^2 \mu = 1

const int size_x=512*4; // # steps in transverse direction
const int N_Y=100; // # steps in long direction

const double L_x=31;// step_x*size_x; // transverse extent
const double step_x=L_x/size_x; // transvers step size


const double step_x2 = step_x*step_x;
const int size_x2 = size_x*size_x;

double UV = step_x*step_x*1e-5;

string eventID;

boost::normal_distribution<> nd( 0.0, 1.0 / (step_x * sqrt(N_Y) ) );
boost::variate_generator< RNGType, boost::normal_distribution<> > Gaussian(rng, nd);


double g2mu = 0.25;

boost::normal_distribution<> nd_p( 0.0, g2mu / (step_x ) );
boost::variate_generator< RNGType, boost::normal_distribution<> > Gaussian_p(rng, nd_p);

typedef blitz::TinyVector<cd,8> AdComp;
typedef blitz::TinyMatrix<cd,3,3>  colorMat;
typedef blitz::TinyMatrix<double,8,8>  colorAdMat;
typedef Array<colorMat,2>  colorArr;
typedef Array<colorAdMat,2>  colorAdArr;

ofstream fileout;

/*
AdComp DecomposeAdMatrix(colorAdMat U)
{
    AdComp out;
    vector<colorAdMat> F=get_F();
    for(int i=0; i<8; i++)
    {
        double component = -1.0/3.0*trace(Product(F.at(i),U));
        out(i)=component;
    }
    return out;
}*/

double int_to_x(int i)
{
    return i*step_x;
}


double Qs(double Y)
{
	return 0.71; 
	//Saturation momentum on Y
    double a = 0.77190536;
    double b = 0.27141908;
    double c = -0.06149744;
    
	return a*exp(b*Y)+c;
}



int x_to_int(double x)
{
    int out = int(x/step_x+0.5);
    return out;
}

double x2(double x)
{
    return x*x;
}

double L_d_by_2pi = L_x/(2.0*M_PI);

void AdWilsonLine(colorAdArr& U, colorArr& V)
{
    vector < TinyMatrix<cd,3,3>  > lambda = get_lambda();
    for(int a=0; a<8; a++)
        for(int b=0; b<8; b++)
        {
            for(int i=0; i<size_x; i++)
                for(int j=0; j<size_x; j++)
                {
                    colorMat Vx, Vxdagger, UnderTr;
                    Vx = V(i,j);
                    Vxdagger = dagger(Vx);
                    UnderTr=Product(Product(lambda.at(a),Vx), Product(lambda.at(b),Vxdagger));
                    double Ux = 2.0*real(trace(UnderTr));
                    U(i,j)(a,b) = Ux;
                }
        }
}

void rho_generator(blitz::Array<double,2>& rho)
{
    for(int i=0; i<size_x; i++)
    {
        for(int j=0; j<size_x; j++)
        {
            rho(i,j) = Gaussian();
        }
    }
}


void rho_generator_p(blitz::Array<double,2>& rho)
{
    for(int i=0; i<size_x; i++)
    {
        for(int j=0; j<size_x; j++)
        {
            rho(i,j) = Gaussian_p();
        }
    }
}


blitz::Array<complex<double>,2> A_image (blitz::Array<complex<double>,2>&  rho_image)
{
    blitz::Array<complex<double>,2> out(size_x,size_x);
    out=0.0;
    for(int i=0; i<size_x; i++)
        for(int j=0; j<size_x; j++)
            out(i,j) = - 0.5*step_x*step_x * rho_image(i,j)
                       /
                       ( cos(2.0*M_PI*i/double(size_x)) + cos(2.0*M_PI*j/double(size_x)) - 2.0 );
    out(0,0)=cd(0.0,0.0); //removing zero mode
    return out;
}


int index_x_boundary(int in)
{
    //periodic boundaries
    int i = in;
    if(i<0) i = i + size_x;
    if(i>size_x-1) i = i - size_x;
    return i;
}


TinyMatrix<cd,3,3> fV(int i, int j, vector<Array<double,2> >&  A_a)
{
    vector < TinyMatrix<cd,3,3>  > lambda = get_lambda();

    blitz::TinyMatrix<cd,3,3>  V_out;

    V_out=cd(1.0,0.0),cd(0.0,0.0),cd(0.0,0.0),
    cd(0.0,0.0),cd(1.0,0.0),cd(0.0,0.0),
    cd(0.0,0.0),cd(0.0,0.0),cd(1.0,0.0);
    for(int ny=0; ny<N_Y; ny++)
    {

        blitz::TinyMatrix<cd,3,3>  in_exp;
        blitz::TinyMatrix<cd,3,3>  V;

        in_exp=cd(0.0,0.0);

        for(int a=0; a<8; a++)
            in_exp = in_exp + lambda.at(a) * A_a.at(ny*8+a)(i,j);
        in_exp=in_exp; // times i in the matrixix exponent
        V = matrix_exp_an_antiH(in_exp);
        V_out=Product(V,V_out);
    }

    return V_out;
}

blitz::Array<complex<double>,2> fft(blitz::Array<double,2>& rho)
{

    blitz::Array<complex<double>,2> image(size_x,size_x);

    blitz::Array<complex<double>,2> in(size_x,size_x);
    in = rho;
    FFTW(in,image);

    return image;
}

void IC_MV( colorArr& V )
{
    blitz::Array<double,2> rho_1(size_x, size_x);
    blitz::Array<double,2> rho_2(size_x, size_x);
    blitz::Array<double,2> rho_3(size_x, size_x);
    blitz::Array<double,2> rho_4(size_x, size_x);
    blitz::Array<double,2> rho_5(size_x, size_x);
    blitz::Array<double,2> rho_6(size_x, size_x);
    blitz::Array<double,2> rho_7(size_x, size_x);
    blitz::Array<double,2> rho_8(size_x, size_x);

    blitz::Array<double,2> rho_a[8]= {rho_1,rho_2,rho_3,rho_4,rho_5,rho_6,rho_7,rho_8};

    vector<blitz::Array<double,2> >  A_a(8);


    for(int i=0; i<size_x; i=i+1)
    {
        for(int j=0; j<size_x; j=j+1)
        {

            blitz::TinyMatrix<cd,3,3>  V_unit;
            V_unit=cd(1.0,0.0),cd(0.0,0.0),cd(0.0,0.0),
            cd(0.0,0.0),cd(1.0,0.0),cd(0.0,0.0),
            cd(0.0,0.0),cd(0.0,0.0),cd(1.0,0.0);
            V(i,j)=V_unit;
        }
    }

    for(int ny=0; ny<N_Y; ny++)
    {
        for(int i=0; i<8; i++) rho_generator(rho_a[i]);

        for(int i=0; i<8; i++)
        {
            blitz::Array<cd,2>  image(size_x, size_x);
            blitz::Array<cd,2>  A(size_x, size_x);
            blitz::Array<cd,2>  tmp(size_x, size_x);
            blitz::Array<double,2>  B(size_x, size_x);
            image=fft(rho_a[i]);
            tmp=A_image(image);
            FFTW_b(tmp, A);
            A=A /  double(size_x*size_x);
            B=real(A);
            A_a.at(i).resize(B.shape());
            A_a.at(i)=B;
        }

        for(int i=0; i<size_x; i=i+1)
        {
            for(int j=0; j<size_x; j=j+1)
            {
                vector < TinyMatrix<cd,3,3>  > lambda = get_lambda();


                blitz::TinyMatrix<cd,3,3>  in_exp;
                blitz::TinyMatrix<cd,3,3>  V_at_ij;

                in_exp=cd(0.0,0.0);

                for(int a=0; a<8; a++)
                    in_exp = in_exp + lambda.at(a) * A_a.at(a)(i,j);

                V_at_ij = matrix_exp_an_antiH(in_exp);
                V(i,j) = Product(V(i,j),V_at_ij);
            }
        }
    }
}

void IC_p(  vector<blitz::Array<double,2> >& A_a )
{

    blitz::Array<double,2> rho_1(size_x, size_x);
    blitz::Array<double,2> rho_2(size_x, size_x);
    blitz::Array<double,2> rho_3(size_x, size_x);
    blitz::Array<double,2> rho_4(size_x, size_x);
    blitz::Array<double,2> rho_5(size_x, size_x);
    blitz::Array<double,2> rho_6(size_x, size_x);
    blitz::Array<double,2> rho_7(size_x, size_x);
    blitz::Array<double,2> rho_8(size_x, size_x);

    blitz::Array<double,2> rho_a[8]= {rho_1,rho_2,rho_3,rho_4,rho_5,rho_6,rho_7,rho_8};

    for(int i=0; i<8; i++) rho_generator_p(rho_a[i]);

    for(int i=0; i<8; i++)
    {
        blitz::Array<cd,2>  image(size_x, size_x);
        blitz::Array<cd,2>  A(size_x, size_x);
        blitz::Array<cd,2>  tmp(size_x, size_x);
        blitz::Array<double,2>  B(size_x, size_x);
        image=fft(rho_a[i]);
        tmp=A_image(image);
        FFTW_b(tmp, A);
        A=A /  double(size_x*size_x);
        B=real(A);
        A_a.at(i).resize(B.shape());
        A_a.at(i)=B;
    }
}




#include <sys/stat.h>
#include "analyzer.h"


#include "nr3.h"
#include "interp_1d.h"
#include "interp_linear.h"
#include "interp_2d.h"
#include "fourier.h"


string dname;

complex<double> su3_group_element(colorMat V,  int a)
{
    vector < TinyMatrix<cd,3,3>  > lambda = get_lambda();
    colorMat P;
    P = Product(V, lambda.at(a));
    return 2.0*trace(P);
}

void output(double Y,colorArr& V_c)
{
    string fname;
    ofstream d_data;
    fname = dname+"/S_" + toString(Y) + ".dat";
    d_data.open(fname.c_str());



    blitz::Array<cd,2>  comp(size_x, size_x);
    blitz::Array<cd,2>  compFT(size_x, size_x);
    blitz::Array<cd,2>  sum(size_x, size_x);
    blitz::Array<cd,2>  sum_all(size_x, size_x);
    blitz::Array<cd,2>  sumqPerp(size_x, size_x);
    blitz::Array<cd,2>  pS(size_x, size_x);
    blitz::Array<cd,2>  pSPerp(size_x, size_x);
    sum=cd(0.0,0.0);
    sum_all=cd(0.0,0.0);
    sumqPerp=cd(0.0,0.0);
    comp=cd(0.0,0.0);

    int max_k_int = 100;
    double max_k = 15.0;
    double step_k = max_k/max_k_int;
    vector<double> distr(max_k_int);
    vector<double> normal(max_k_int);
    double Qs_Y = Qs(Y);

    for (int i=0; i<max_k; i++)
    {
        distr.at(i) = 0.0;
        normal.at(i) = 0.0;
    }

    for(int a=0; a<9; a++)
    {
        for(int i=0; i<size_x; i=i+1)
        {
            for(int j=0; j<size_x; j=j+1)
            {
                double x = int_to_x(i);
                double y = int_to_x(j);

                //double distance2 = ( x2(sin(x/L_d_by_2pi)) + x2(sin(y/L_d_by_2pi))  ) * x2(L_d_by_2pi );
                //if(distance2*Qs_Y*Qs_Y<1) comp(i,j) = su3_group_element(V_c(i,j), a);
                comp(i,j) = su3_group_element(V_c(i,j), a) ;//* exp(-0.5*distance2*Qs_Y*Qs_Y);
            }
        }

        FFTW(comp, compFT);

        for(int i=0; i<size_x; i=i+1)
        {
            double kx  = 2.0*M_PI*i/L_x;
            double kx_t  = 2.0/step_x*sin(kx*step_x/2.0);
            for(int j=0; j<size_x; j=j+1)
            {

                double ky  = 2.0*M_PI*j/L_x;
                double kx_t  = 1.0/step_x*sin(kx*step_x);
                double ky_t  = 1.0/step_x*sin(ky*step_x);
                double k2 = pow(2.0/step_x*sin(kx*step_x/2),2) +  pow(2.0/step_x*sin(ky*step_x/2),2);



                sum_all(i,j) += compFT(i,j)*conj(compFT(i,j));

                /*if(k2>pow(Qs_Y,2)) {
                    sum(i,j) += compFT(i,j)*conj(compFT(i,j));
                    sumqPerp(i,j) += compFT(i,j)*conj(compFT(i,j)) * k2 ;
                }*/

            }
        }
    }


    /*for (int ik=0;ik<max_k_int;ik++)
    {
    	fileout << Y << " " << ik*step_k << " "
    		    << distr.at(ik)/normal.at(ik) <<"\n" << flush;
    }*/

    /*
    	for(int i=0; i<size_x/4; i++)
        {
       double kx  = 2.0*M_PI*i/L_x;
            double kx_t  = 1.0/step_x*sin(kx*step_x);

            for(int j=0; j<size_x/4; j++)
    		{
    	            double ky  = 2.0*M_PI*j/L_x;
                double ky_t  = 1.0/step_x*sin(ky*step_x);
    			double k  = sqrt( kx_t*kx_t + ky_t*ky_t);
    			double S = real(sum(i,j))/6.0;
    			//cout << S << "\n";
    			d_data <<  k
    			  << " " << S << "\n" ;
    		}
    	}
    */


    FFTW_b(sum_all, pS);

	double dr = (10.1*step_x); 
	double N = int(L_x/dr); 

	vector<double> aS(N);
	vector<double> aN(N);

	for (int i=0; i<N; i++) 
	{
		aS[i]=0.0;
		aN[i]=0.0;
	}

    for(int i=0; i<size_x; i++)
    {
        for(int j=0; j<size_x; j++)
        {
            double r = sqrt(x2(int_to_x(i)) + x2(int_to_x(j)));
            int ir = int(r/dr);
			if(ir<N)
			{
				aS[ir] += 1.0-real(pS (i,j))*0.5/3.0/size_x2/size_x2 ;  
				aN[ir] += 1.0;
			}
            //fileout <<  sqrt(x2(int_to_x(i)) + x2(int_to_x(j)))  << " " << 1.0-real(pS (i,j))*0.5/3.0/size_x2/size_x2 << "\n" << flush;
        }
    }


    for(int ir=0; ir<N; ir++)
    {
        double r = (double(ir)+0.5)*dr;
        if (r<0.5*L_x) 
			fileout <<  r  << " " << aS[ir]/aN[ir] << "\n" << flush;
    }

    //d_data<< real(pS(0,0)) << " " <<  real(pSPerp(0,0)) << "\n" <<  flush;
    //cout << Y << " " << real(pS(0,0))*0.5/3.0/size_x2/size_x2 << " " <<  real(pSPerp(0,0)) *0.5/3.0/size_x2/size_x2<< " "<< Qs(Y) <<  "\n" <<    flush;
}


void Omega(colorArr& Omega_s, colorArr& Omega_a, colorArr V, vector<blitz::Array<double,2> > A_a)
{
    vector < TinyMatrix<cd,3,3>  > lambda = get_lambda();

    colorMat W_x, W_y;

    for(int i=0; i<size_x; i++)
        for(int j=0; j<size_x; j++)
        {


            Omega_s(i,j) = cd(0.0,0.0);
            Omega_a(i,j) = cd(0.0,0.0);


            for(int a=0; a<8; a++)
            {
                //x component
                double proton_part_x  =  0.5*(A_a.at(a)(index_x_boundary(i+1),j) - A_a.at(a)(index_x_boundary(i-1),j))/step_x;
                W_x =  0.5/step_x*(Product( Product( dagger(V(index_x_boundary(i+1),j)),lambda.at(a) ),  V(index_x_boundary(i+1),j) )
                                   -  Product( Product( dagger(V(index_x_boundary(i-1),j)),lambda.at(a) ),  V(index_x_boundary(i-1),j) ));

                //y component
                double proton_part_y  =  0.5*(A_a.at(a)(i,index_x_boundary(j+1)) - A_a.at(a)(i, index_x_boundary(j-1)))/step_x;
                W_y =  0.5/step_x*(Product( Product( dagger(V(i,index_x_boundary(j+1))),lambda.at(a) ),  V(i,index_x_boundary(j+1)) )
                                   -  Product( Product( dagger(V(i,index_x_boundary(j-1))),lambda.at(a) ),  V(i,index_x_boundary(j-1)) ));


                Omega_s(i,j) = Omega_s(i,j) + proton_part_x * W_x + proton_part_y * W_y;
                Omega_a(i,j) = Omega_a(i,j) + proton_part_x * W_y - proton_part_y * W_x;
            }
        }
}


double sign(double x)
{
    if(x<0) return -1.0;
    //if(x*x<UV*UV) return 0;
    return 1.0;
}



double SI(int ik, int jk, vector<blitz::Array<cd,2> > &Omega_s, vector<blitz::Array<cd,2> > & Omega_a)
{
    double kx  = 2.0*M_PI*ik/L_x;
    double ky  = 2.0*M_PI*jk/L_x;
    double kx_t  = 1.0/step_x*sin(kx*step_x);
    double ky_t  = 1.0/step_x*sin(ky*step_x);
    double k2 = pow(2.0/step_x*sin(kx*step_x/2),2) +  pow(2.0/step_x*sin(ky*step_x/2),2);

    cd sum = cd(0,0);

    for(int a=0; a<8; a++)
    {
        sum+=Omega_s.at(a)(ik,jk)*conj(Omega_s.at(a)(ik,jk)) + Omega_a.at(a)(ik,jk)*conj(Omega_a.at(a)(ik,jk));
    }


    return real(sum)/(k2+UV)/pow(2*M_PI,3)/pow(L_x,2);
}

//Forward multiply by a^2
//Backward divide by L^2


int four_bound(int n)
{
    if(n<0) return n+size_x;
    if(n>size_x-1) return n-size_x;
    return n;
}

double AsSI_old(int ik, int jk, vector<blitz::Array<cd,2> > &Omega_s, vector<blitz::Array<cd,2> > & Omega_a)
{

    cd Sum = cd(0.0,0.0);

    for(int i=0; i<size_x; i++)
        for(int j=0; j<size_x; j++)
        {
            {
                Sum+=
                    Omega_s.at(0)(i,j)*Omega_s.at(0)(four_bound(ik-i),four_bound(jk-j));

            }
        }


    //return real(Integrated(0,0))/pow(L_x,2);
    //return imag(Integrated(0,0))/pow(L_x,2);
    return imag(Sum) /pow(L_x,2);
}


double AsSI(int ik, int jk, vector<blitz::Array<cd,2> > &Omega_s, vector<blitz::Array<cd,2> > & Omega_a)
{
	//return 0;
    vector < TinyMatrix<double,8,8>  > F = get_F();
    blitz::Array<cd,2> UInt(size_x,size_x);
    blitz::Array<cd,2> Integrated(size_x,size_x);

    double kx  = 2.0*M_PI*ik/L_x;
    double ky  = 2.0*M_PI*jk/L_x;
    double kx_t  = 1.0/step_x*sin(kx*step_x);
    double ky_t  = 1.0/step_x*sin(ky*step_x);
    double k2 = pow(2.0/step_x*sin(kx*step_x/2),2) +  pow(2.0/step_x*sin(ky*step_x/2),2);

    for(int i=0; i<size_x; i++)
        for(int j=0; j<size_x; j++)
        {
            double qx = 2.0*M_PI*i/L_x;
            double qy = 2.0*M_PI*j/L_x;


            UInt(i,j) = cd(0.0,0.0);

            //if(((i!=0)||(j!=0))&&((index_x_boundary(-i+ik)!=0)||(index_x_boundary(-j+jk)!=0)))
            {

                double qx_t  = 1.0/step_x*sin(qx*step_x);
                double qy_t  = 1.0/step_x*sin(qy*step_x);

                double qxnkx = 1.0/step_x*sin(2.0*M_PI*four_bound(-i+ik)/L_x*step_x);
                double qynky = 1.0/step_x*sin(2.0*M_PI*four_bound(-j+jk)/L_x*step_x) ;

                double k_dot_k_minus_q = kx_t * qxnkx + ky_t * qynky;
                double q_dot_k_minus_q = qx_t * qxnkx + qy_t * qynky;

                double q2 = pow(2.0/step_x*sin(qx*step_x/2),2) +  pow(2.0/step_x*sin(qy*step_x/2),2);
                double kminusq2 = pow(2.0/step_x*sin((2.0*M_PI*four_bound(-i+ik)/L_x)*step_x/2),2)
                                  +  pow(2.0/step_x*sin((2.0*M_PI*four_bound(-j+jk)/L_x)*step_x/2),2) ;
                double kXq = qxnkx*qy_t- qynky *qx_t;


                cd sum = cd(0,0);

                //if((q2>=pow(2*M_PI/L_x,2))&&(kminusq2>=pow(2*M_PI/L_x,2)))
                //if(q2>=pow(2*M_PI/L_x,2))
                for(int a=0; a<8; a++)
                    for(int b=0; b<8; b++)
                        for(int c=0; c<8; c++)
                        {
                            sum+=
                                F.at(a)(b,c) *
                                sign(kXq) /( (q2+UV)*(kminusq2+UV)) *
                                (
                                    (k2 *  Omega_a.at(a)(i,j) * Omega_a.at(b)(four_bound(ik-i),four_bound(jk-j))
                                     -  q_dot_k_minus_q * ( Omega_a.at(a)(i,j) * Omega_a.at(b)(four_bound(ik-i),four_bound(jk-j))
                                                            + Omega_s.at(a)(i,j) * Omega_s.at(b)(four_bound(ik-i),four_bound(jk-j)) )
                                    ) * conj(Omega_a.at(c)(ik,jk))
                                    + 2.0 * k_dot_k_minus_q * Omega_a.at(a)(i,j) * Omega_s.at(b)(four_bound(ik-i),four_bound(jk-j)) * conj(Omega_a.at(c)(ik,jk))
                                );
                        }

                UInt(i,j) = sum;
            }
        }

    FFTW_b(UInt, Integrated);

    //return real(Integrated(0,0))/pow(L_x,2);
    //return imag(Integrated(0,0))/pow(L_x,2);
    return imag(Integrated(0,0))/(k2+UV)/pow(L_x,4)/pow(2*M_PI,3); //per unit length
}



void components(colorArr& Omega, vector<blitz::Array<cd,2> > & Omega_comp)
{
    for(int a=0; a<8; a++)
    {
        Omega_comp.at(a).resize(Omega.shape());

        for(int i=0; i<size_x; i++)
            for(int j=0; j<size_x; j++)
            {
                Omega_comp.at(a)(i,j) = su3_group_element(Omega(i,j), a);
            }
    }
}

int sign(int i)
{
    if(i<0) return -1;
    return 1;
}


int main(void)
{
    cin >> eventID;
    dname = "v3_data_test"+eventID;
    mkdir(dname.c_str(),S_IRWXU | S_IRWXG);

    string name = dname+"/MD_" + toString(eventID) + ".dat";
    fileout.open(name.c_str());

    //Target:
    colorArr V_c(size_x,size_x);
    IC_MV(V_c);
    cerr << "target done\n" << endl;

    //Projectile:
    //
    vector<blitz::Array<double,2> > A_a(8);
    IC_p(A_a);
    cerr << "projectile done\n" << endl;

    colorArr Omega_a(size_x,size_x), Omega_s(size_x,size_x);
    vector<blitz::Array<cd,2> > Omega_a_c(8), Omega_s_c(8);
    vector<blitz::Array<cd,2> > fftOmega_a_c(8), fftOmega_s_c(8);

    Omega(Omega_s, Omega_a, V_c, A_a);
    cerr << "omega done\n" << endl<< flush;

    components(Omega_s, Omega_s_c);
    components(Omega_a, Omega_a_c);

    for(int a=0; a<8; a++)
    {
        fftOmega_a_c.at(a).resize(Omega_a.shape());
        fftOmega_s_c.at(a).resize(Omega_a.shape());

        blitz::Array<cd,2> tmp1(size_x,size_x);
        blitz::Array<cd,2> tmp2(size_x,size_x);
        tmp1 = Omega_a_c.at(a);
        FFTW(tmp1,tmp2);
        fftOmega_a_c.at(a)=tmp2*pow(step_x,2);

        tmp1 = Omega_s_c.at(a);
        FFTW(tmp1,tmp2);
        fftOmega_s_c.at(a)=tmp2*pow(step_x,2);

        cerr << fftOmega_a_c.at(a)(0,0) << " " << fftOmega_s_c.at(a)(0,0) << "\n" << flush;
    }

    cerr << "components done\n" << endl << flush;

    output(0.0, V_c);

    cerr << "output done\n" << endl;

    double Qs2= pow(Qs(0),2);


    typedef std::map<std::pair<int, int>, double > Map;
    Map Kx, Ky, mapSI, mapASI;

    /*	mymap[std::make_pair(2, 4)] = 10;
    ...
    AMapT::iterator f = mymap.find(std::make_pair(3, 5));
    */

    double dK = 0.5;
	int Nphi =128; 
    double dphi = 2.0*M_PI/Nphi;
    for(double K= dK; K<5+dK; K+=dK)
    {        
		for(int iphi=0; iphi<Nphi; iphi+=1)
        {
			double phi = dphi*iphi; 
            double k_x = K*cos(phi);
            double k_y = K*sin(phi);

            double lat_k_x = 1.0/step_x*asin(k_x*step_x);
            double lat_k_y = 1.0/step_x*asin(k_y*step_x);

            int int_k_x = int(lat_k_x * L_x / (2.0*M_PI)) ;
            int int_k_y = int(lat_k_y * L_x / (2.0*M_PI)) ;

            Map::iterator f = mapSI.find(std::make_pair( four_bound(int_k_x), four_bound(int_k_y) ));
            if (f==mapSI.end())
            {
                double si = SI(four_bound(int_k_x), four_bound(int_k_y), fftOmega_s_c, fftOmega_a_c);
                double asi = AsSI(four_bound(int_k_x), four_bound(int_k_y), fftOmega_s_c, fftOmega_a_c);
                mapSI[std::make_pair( four_bound(int_k_x), four_bound(int_k_y))] = si;
                mapASI[std::make_pair( four_bound(int_k_x), four_bound(int_k_y))] = asi;
            }

            int s = sign(int_k_x);
            f = mapSI.find(std::make_pair(four_bound(int_k_x+s),four_bound(int_k_y)));
            if (f==mapSI.end())
            {
                double si = SI(four_bound(int_k_x+s), four_bound(int_k_y), fftOmega_s_c, fftOmega_a_c);
                double asi = AsSI(four_bound(int_k_x+s), four_bound(int_k_y), fftOmega_s_c, fftOmega_a_c);
                mapSI[std::make_pair( four_bound(int_k_x+s), four_bound(int_k_y))] = si;
                mapASI[std::make_pair( four_bound(int_k_x+s), four_bound(int_k_y))] = asi;
            }

            s = sign(int_k_y);
            f = mapSI.find(std::make_pair( four_bound(int_k_x) , four_bound(int_k_y+s) ));
            if (f==mapSI.end())
            {
                double si = SI(four_bound(int_k_x), four_bound(int_k_y+s), fftOmega_s_c, fftOmega_a_c);
                double asi = AsSI(four_bound(int_k_x), four_bound(int_k_y+s), fftOmega_s_c, fftOmega_a_c);
                mapSI[std::make_pair(four_bound(int_k_x), four_bound(int_k_y+s))] = si;
                mapASI[std::make_pair( four_bound(int_k_x), four_bound(int_k_y+s))] = asi;
            }

            int sx = sign(int_k_x);
            int sy = sign(int_k_y);
            f = mapSI.find(std::make_pair(four_bound(int_k_x+sx), four_bound(int_k_y+sy)));
            if (f==mapSI.end())
            {
                double si = SI(four_bound(int_k_x+sx), four_bound(int_k_y+sy), fftOmega_s_c, fftOmega_a_c);
                double asi = AsSI(four_bound(int_k_x+sx), four_bound(int_k_y+sy), fftOmega_s_c, fftOmega_a_c);
                mapSI[std::make_pair( four_bound(int_k_x+sx), four_bound(int_k_y+sy))] = si;
                mapASI[std::make_pair( four_bound(int_k_x+sx), four_bound(int_k_y+sy))] = asi;
            }
            cerr << K << " " << phi << endl << flush;
        }
    }

    for(double K=dK; K<5+dK; K+=dK)
    {
        for(int iphi=0; iphi<Nphi; iphi+=1)
        {
			double phi = dphi*iphi; 
            double k_x = K*cos(phi);
            double k_y = K*sin(phi);

            double lat_k_x = 1.0/step_x*asin(k_x*step_x);
            double lat_k_y = 1.0/step_x*asin(k_y*step_x);

            int int_k_x = int(lat_k_x * L_x / (2.0*M_PI)) ;
            int int_k_y = int(lat_k_y * L_x / (2.0*M_PI)) ;

            int sx = sign(int_k_x);
            int sy = sign(int_k_y);

            double k_x_m = 1.0/step_x*sin(2.0*M_PI*step_x/L_x*four_bound(int_k_x));
            double k_y_m = 1.0/step_x*sin(2.0*M_PI*step_x/L_x*four_bound(int_k_y));

            double k_x_m_shift = 1.0/step_x*sin(2.0*M_PI*step_x/L_x*four_bound(int_k_x+sx));
            double k_y_m_shift = 1.0/step_x*sin(2.0*M_PI*step_x/L_x*four_bound(int_k_y+sy));

            double t = (k_x - k_x_m)/(k_x_m_shift - k_x_m);
            double u = (k_y - k_y_m)/(k_y_m_shift - k_y_m);

			cerr << k_x << " " << k_x_m << " " << k_x_m_shift << "\n" << flush;

            double value_SI = (1-t)*(1-u)*mapSI[std::make_pair(four_bound(int_k_x), four_bound(int_k_y))]
                              + t*(1-u)*mapSI[std::make_pair( four_bound(int_k_x+sx), four_bound(int_k_y))]
                              + (1-t)*u*mapSI[std::make_pair(four_bound( int_k_x), four_bound(int_k_y+sy))]
                              + t*u*mapSI[std::make_pair( four_bound(int_k_x+sx), four_bound(int_k_y+sy))];

            double value_ASI = (1-t)*(1-u)*mapASI[std::make_pair(four_bound(int_k_x), four_bound(int_k_y))]
                               + t*(1-u)*mapASI[std::make_pair( four_bound(int_k_x+sx), four_bound(int_k_y))]
                               + (1-t)*u*mapASI[std::make_pair( four_bound(int_k_x), four_bound(int_k_y+sy))]
                               + t*u*mapASI[std::make_pair( four_bound(int_k_x+sx),four_bound( int_k_y+sy))];

            cout << K << " " << phi << " " << value_SI << " " << value_ASI << endl;
        }

    }

    return 1;

	int i = size_x/8;

        for(int j=0; j<size_x; j++)
        {
            double kx  = 2.0*M_PI*i/L_x;
            double ky  = 2.0*M_PI*j/L_x;
            double kx_t  = 1.0/step_x*sin(kx*step_x);
            double ky_t  = 1.0/step_x*sin(ky*step_x);
            double k2 = pow(2.0/step_x*sin(kx*step_x/2),2) +  pow(2.0/step_x*sin(ky*step_x/2),2);
 			if(k2<10*10){ 
            double si = SI(i, j, fftOmega_s_c, fftOmega_a_c);
            double asi = AsSI(i, j, fftOmega_s_c, fftOmega_a_c);


           cout << kx_t << " " << ky_t << " " << si << " " << asi <<  " ";
			
            kx  = 2.0*M_PI*(size_x-i)/L_x;
            ky  = 2.0*M_PI*(size_x-j)/L_x;
            kx_t  = 1.0/step_x*sin(kx*step_x);
            ky_t  = 1.0/step_x*sin(ky*step_x);
            
			si = SI(four_bound(size_x-i), four_bound(size_x-j), fftOmega_s_c, fftOmega_a_c);
            asi = AsSI(four_bound(size_x-i),four_bound(size_x-j), fftOmega_s_c, fftOmega_a_c);
           
			cout << si << " " << asi <<  endl;

        }
		}




	return 1;



    for(int i=0; i<size_x/2; i++)
        for(int j=0; j<size_x/2; j++)
        {
            double kx  = 2.0*M_PI*i/L_x;
            double ky  = 2.0*M_PI*j/L_x;
            double kx_t  = 1.0/step_x*sin(kx*step_x);
            double ky_t  = 1.0/step_x*sin(ky*step_x);
            double k2 = pow(2.0/step_x*sin(kx*step_x/2),2) +  pow(2.0/step_x*sin(ky*step_x/2),2);

            if((k2<50.0*Qs2)&&(k2>UV))
            {
                //cout <<  "# " << i << " " << j << " " << endl;
                double si = SI(i, j, fftOmega_s_c, fftOmega_a_c);
                double asi = AsSI(i, j, fftOmega_s_c, fftOmega_a_c);
                cout << kx << " " << ky << " " << si << " " << asi << endl;
            }
        }
    return 1;

}
