#include "JIMWLK.h"
#include "string_conv.h"
#include <map>
//#include <omp.h>


using namespace std;
using namespace blitz;

//choice of parameters g^2 \mu = 1

const int size_x=512; // # steps in transverse direction
const int N_Y=100; // # steps in long direction

const double L_x=32;// step_x*size_x; // transverse extent
const double step_x=L_x/size_x; // transvers step size

const double step_x2 = step_x*step_x;
const int size_x2 = size_x*size_x;

double UV = step_x*step_x*1e-5;

double L_d_by_2pi = L_x/(2.0*M_PI);
double step_k=2.0*M_PI/L_x; 

string eventID;

boost::normal_distribution<> nd( 0.0, 1.0 / (step_x * sqrt(N_Y) ) );
boost::variate_generator< RNGType, boost::normal_distribution<> > Gaussian(rng, nd);

double g2mu = 0.5;

boost::normal_distribution<> nd_p( 0.0, g2mu / (step_x ) );
boost::variate_generator< RNGType, boost::normal_distribution<> > Gaussian_p(rng, nd_p);

typedef blitz::TinyVector<cd,8> AdComp;
typedef blitz::TinyMatrix<cd,3,3>  colorMat;
typedef blitz::TinyMatrix<double,8,8>  colorAdMat;
typedef Array<colorMat,2>  colorArr;
typedef Array<colorAdMat,2>  colorAdArr;

ofstream fileout;

double A=1.0;
double B=1.0/A;


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
    return 0.82;
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

double int_to_k(int i)
{
    double k  = 2.0*M_PI*i/L_x;
    double k_t  = 1.0/step_x*sin(k*step_x);
    return k_t;
}

double int_to_k2(int i, int j)
{
	double kx  = 2.0*M_PI*i/L_x;
    double ky  = 2.0*M_PI*j/L_x;
    return  pow(2.0/step_x*sin(kx*step_x/2),2) +  pow(2.0/step_x*sin(ky*step_x/2),2);
}



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
    double xc = int_to_x(size_x/2);
    double yc = int_to_x(size_x/2);
    double R = 2.0;

    for(int i=0; i<size_x; i++)
    {
        for(int j=0; j<size_x; j++)
        {
            double x = int_to_x(i);
            double y = int_to_x(j);

            double ellipse = pow((x-xc)/(A*R),2)  + pow((y-yc)/(B*R),2) ;

            rho(i,j) = exp(-ellipse/2.0) * Gaussian_p();
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

TinyMatrix<cd,3,3> fV(int i, int j, const vector<Array<double,2> >&  A_a)
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

        cerr << "target slices" << ny << "\n";

        vector<blitz::Array<double,2> >  A_a(8);

        for(int i=0; i<8; i++)
            A_a.at(i).resize(size_x, size_x);

		for(int i=0; i<8; i++)
        {
			blitz::Array<double,2> rho_1(size_x, size_x);
			rho_generator(rho_1);

            blitz::Array<cd,2>  image(size_x, size_x);
            blitz::Array<cd,2>  A(size_x, size_x);
            blitz::Array<cd,2>  tmp(size_x, size_x);
            
			blitz::Array<double,2>  B(size_x, size_x);
            image=fft(rho_1);

            tmp=A_image(image);
            FFTW_b(tmp, A);
            A=A /  double(size_x*size_x);
            B=real(A);
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

complex<double> su3_group_element(const colorMat& V,  int a)
{
    vector < TinyMatrix<cd,3,3>  > lambda = get_lambda();
    colorMat P;
    P = Product(V, lambda.at(a));
    return 2.0*trace(P);
}

void output(double Y,const colorArr& V_c)
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


void Omega_old(colorArr& Omega_s, colorArr& Omega_a, const colorArr& V, const vector<blitz::Array<double,2> >& A_a)
{
    vector < TinyMatrix<cd,3,3>  > lambda = get_lambda();

    colorMat W_x, W_y;
    colorMat E_x, E_y;
    colorMat  Ua;
   	colorMat  Ud;
   	
	colorMat  diff;

    for(int i=0; i<size_x; i++)
        for(int j=0; j<size_x; j++)
        {

            Omega_s(i,j) = cd(0.0,0.0);
            Omega_a(i,j) = cd(0.0,0.0);


            for(int a=0; a<8; a++)
            {
            	Ud = dagger ( V(i,j) ); 
				Ua = V(i,j); 
            	Ua= Product( lambda.at(a), Ua);
            	Ua= Product( Ud , Ua);

				diff = 0.5*( V(index_x_boundary(i+1),j) -   V(index_x_boundary(i-1),j) )  / step_x;  
				E_x = Product(Ud, diff);
				diff = 0.5*( V(i, index_x_boundary(j+1)) -   V(i, index_x_boundary(j-1)) )  / step_x;
				E_y = Product(Ud, diff); 
            	//x component
                double proton_part_x  =  0.5*(A_a.at(a)(index_x_boundary(i+1),j) - A_a.at(a)(index_x_boundary(i-1),j))/step_x;
				W_x = Product(Ua, E_x) -  Product(E_x, Ua);  
               // W_x =  0.5/step_x*(Product( Product( dagger(V(index_x_boundary(i+1),j)),lambda.at(a) ),  V(index_x_boundary(i+1),j) )
                //                   -  Product( Product( dagger(V(index_x_boundary(i-1),j)),lambda.at(a) ),  V(index_x_boundary(i-1),j) ));

                //y component
                double proton_part_y  =  0.5*(A_a.at(a)(i,index_x_boundary(j+1)) - A_a.at(a)(i, index_x_boundary(j-1)))/step_x;
				W_y = Product(Ua, E_y) -  Product(E_y, Ua);  
                //W_y =  0.5/step_x*(Product( Product( dagger(V(i,index_x_boundary(j+1))),lambda.at(a) ),  V(i,index_x_boundary(j+1)) )
                 //                  -  Product( Product( dagger(V(i,index_x_boundary(j-1))),lambda.at(a) ),  V(i,index_x_boundary(j-1)) ));

                Omega_s(i,j) = Omega_s(i,j) + proton_part_x * W_x + proton_part_y * W_y;
                Omega_a(i,j) = Omega_a(i,j) + proton_part_x * W_y - proton_part_y * W_x;
            }
        }
}



void Omega(colorArr& Omega_s, colorArr& Omega_a, const colorArr& V, const vector<blitz::Array<double,2> >& A_a)
{
    vector < TinyMatrix<cd,3,3>  > lambda = get_lambda();

    colorMat W_x, W_y;
   
	colorMat  tmpPlus;
	colorMat  tmpMinus;

    for(int i=0; i<size_x; i++)
        for(int j=0; j<size_x; j++)
        {

            Omega_s(i,j) = cd(0.0,0.0);
            Omega_a(i,j) = cd(0.0,0.0);


            for(int a=0; a<8; a++)
            {

				tmpPlus = V(index_x_boundary(i+1),j); 			
				tmpPlus = Product( lambda.at(a), tmpPlus); 			
				tmpPlus = Product( dagger(V(index_x_boundary(i+1),j)), tmpPlus); 			
				
				tmpMinus = V(index_x_boundary(i-1),j); 			
				tmpMinus = Product( lambda.at(a), tmpMinus); 			
				tmpMinus = Product( dagger(V(index_x_boundary(i-1),j)), tmpMinus); 	

				W_x = 0.5*(tmpPlus - tmpMinus)/step_x;  

				tmpPlus = V(i,index_x_boundary(j+1)); 			
				tmpPlus = Product( lambda.at(a), tmpPlus); 			
				tmpPlus = Product( dagger(V(i,index_x_boundary(j+1))), tmpPlus); 			
				
				tmpMinus = V(i,index_x_boundary(j-1)); 			
				tmpMinus = Product( lambda.at(a), tmpMinus); 			
				tmpMinus = Product( dagger(V(i, index_x_boundary(j-1))), tmpMinus); 	

				W_y = 0.5*(tmpPlus - tmpMinus)/step_x;  


				/*W_x =  0.5/step_x*(Product( Product( dagger(V(index_x_boundary(i+1),j)),lambda.at(a) ),  V(index_x_boundary(i+1),j) )
                                   -  Product( Product( dagger(V(index_x_boundary(i-1),j)),lambda.at(a) ),  V(index_x_boundary(i-1),j) ));

                W_y =  0.5/step_x*(Product( Product( dagger(V(i,index_x_boundary(j+1))),lambda.at(a) ),  V(i,index_x_boundary(j+1)) )
                                   -  Product( Product( dagger(V(i,index_x_boundary(j-1))),lambda.at(a) ),  V(i,index_x_boundary(j-1)) ));
*/



                double proton_part_x  =  0.5*(A_a.at(a)(index_x_boundary(i+1),j) - A_a.at(a)(index_x_boundary(i-1),j))/step_x;
                double proton_part_y  =  0.5*(A_a.at(a)(i,index_x_boundary(j+1)) - A_a.at(a)(i, index_x_boundary(j-1)))/step_x;

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

int four_bound(int n)
{
    if(n<0) return n+size_x;
    if(n>size_x-1) return n-size_x;
    return n;
}


cd SI(const int ik, const int jk, const vector<blitz::Array<cd,2> > &Omega_s, const vector<blitz::Array<cd,2> > & Omega_a)
{
	//cout << ik << " " << jk << "\n"; 

	//double kx_t  = int_to_k(ik);
    //double ky_t  = int_to_k(jk);
    double k2 = int_to_k2(ik,jk);

    cd sum = cd(0,0);

    for(int a=0; a<8; a++)
    {
        sum+=Omega_s.at(a)(ik,jk)*Omega_s.at(a)(four_bound(size_x-ik),four_bound(size_x-jk)) + Omega_a.at(a)(ik,jk)*Omega_a.at(a)(four_bound(size_x-ik),four_bound(size_x-jk));
    }

    return (sum)/(k2+1e-7)/pow(2*M_PI,3)/pow(L_x,2);
    //return Omega_s.at(1)(ik,jk)/pow(2*M_PI,3)/pow(L_x,2);
}

//Forward multiply by a^2
//Backward divide by L^2





cd AsSI(const int ik, const int jk, const vector<blitz::Array<cd,2> > &Omega_s, const vector<blitz::Array<cd,2> > & Omega_a)
{
	
	vector < TinyMatrix<double,8,8>  > F = get_F();

    double kx_t  = int_to_k(ik);
    double ky_t  = int_to_k(jk);
    double k2 = 
		int_to_k2(ik,jk);
		//x2(kx_t) + x2(ky_t); //int_to_k2(ik,jk);

    cd Integrate = cd(0.0,0.0);
	for(int i=0; i<size_x; i++)
	{
     	//if( (i<size_x/4+1)||(i>3*size_x/4-1) )
		{
		for(int j=0; j<size_x; j++)
        {
     		//if( (j<size_x/4+1)||(j>3*size_x/4-1) )
            {

                int itemp = four_bound(ik-i);
				int jtemp = four_bound(jk-j);
				
				double qx_t  = int_to_k(i);
                double qy_t  = int_to_k(j);
                //double q2 = x2(qx_t) + x2(qy_t); //int_to_k2(i,j);
                double q2 = int_to_k2(i,j);

                double qxnkx =  int_to_k(itemp); 
                double qynky = 	int_to_k(jtemp);

                double k_dot_k_minus_q = kx_t * qxnkx + ky_t * qynky;
                double q_dot_k_minus_q = qx_t * qxnkx + qy_t * qynky;

                //double kminusq2 = x2(qxnkx) + x2(qynky) ; //int_to_k2(itemp,jtemp);
                double kminusq2 = int_to_k2(itemp,jtemp);

                double kXq = qxnkx*qy_t - qynky*qx_t;

                cd sum = cd(0,0);

                for(int a=0; a<8; a++)
                    for(int b=0; b<8; b++)
                        for(int c=0; c<8; c++)
                        {
							if(F.at(a)(b,c)*F.at(a)(b,c)>0){  
							sum+=
                                F.at(a)(b,c) *
                                sign(kXq) /( (q2+UV)*(kminusq2+UV)) *
                                (
                                    (k2 *  Omega_a.at(a)(i,j) * Omega_a.at(b)(itemp,jtemp)
                                     -  q_dot_k_minus_q * ( Omega_a.at(a)(i,j) * Omega_a.at(b)(itemp,jtemp)
                                                            + Omega_s.at(a)(i,j) * Omega_s.at(b)(itemp,jtemp) )
                                    ) * (Omega_a.at(c)(four_bound(size_x-ik),four_bound(size_x-jk)))
                                    + 2.0 * k_dot_k_minus_q * Omega_a.at(a)(i,j) * Omega_s.at(b)(itemp,jtemp) * (Omega_a.at(c)(four_bound(size_x-ik),four_bound(size_x-jk)))
									-
                                    (
									 (k2 *  Omega_a.at(a)(four_bound(size_x-i),four_bound(size_x-j)) * Omega_a.at(b)(four_bound(size_x-itemp),four_bound(size_x-jtemp))
                                     -  q_dot_k_minus_q * ( Omega_a.at(a)(four_bound(size_x-i),four_bound(size_x-j)) * Omega_a.at(b)(four_bound(size_x-itemp),four_bound(size_x-jtemp))
                                                            + Omega_s.at(a)(four_bound(size_x-i),four_bound(size_x-j)) * Omega_s.at(b)(four_bound(size_x-itemp),four_bound(size_x-jtemp)) )
                                    ) * (Omega_a.at(c)(ik,jk))
                                    + 2.0 * k_dot_k_minus_q * Omega_a.at(a)(four_bound(size_x-i),four_bound(size_x-j)) * Omega_s.at(b)(four_bound(size_x-itemp),four_bound(size_x-jtemp)) * (Omega_a.at(c)(ik,jk))
									)
                                );
							}
                        }
				Integrate += sum; 
            }
        }
		}
	}

    return cd(0.0,1.0)*(0.5*Integrate)/(k2)/pow(L_x,4)/pow(2*M_PI,3); //per unit length
}



void components(const colorArr& Omega, vector<blitz::Array<cd,2> > & Omega_comp)
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

class ParallelStream{
    std::ostringstream stdStream;
public:
    ParallelStream(){}
    template <class T>
    ParallelStream& operator<<(const T& inData){
        stdStream << inData;
        return *this;
    }
    std::string toString() const{
        return stdStream.str();
    }
};







int main(void)
{
	clock_t begin = clock();
	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

cin >> eventID;

	cin >> A;
	cin >> B;
    //dname = "v3_data_2018_"+eventID;
    //mkdir(dname.c_str(),S_IRWXU | S_IRWXG);

    //string name = dname+"/MD_" + eventID + ".dat";
    //fileout.open(name.c_str());

    //Target:
    colorArr V_c(size_x,size_x);
    begin = clock();
    IC_MV(V_c);
    end = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cerr << "target done\n" << "time " << elapsed_secs << flush  << endl;

    //Projectile:
    //
    vector<blitz::Array<double,2> > A_a(8);
    begin = clock();
    IC_p(A_a);
    end = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cerr << "projectile done\n" << "time " << elapsed_secs << flush  << endl;

    begin = clock();
    colorArr Omega_a(size_x,size_x), Omega_s(size_x,size_x);
    vector<blitz::Array<cd,2> > Omega_a_c(8), Omega_s_c(8);
    vector<blitz::Array<cd,2> > fftOmega_a_c(8), fftOmega_s_c(8);

    Omega(Omega_s, Omega_a, V_c, A_a);
    end = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    cerr << "omega done\n" << "time " << elapsed_secs << flush <<  endl;

    //output(0.0, V_c);

    V_c.free();
    for (int a=0; a<8; a++) A_a.at(a).free();
    cerr << "output done\n" << endl;

    components(Omega_s, Omega_s_c);
    components(Omega_a, Omega_a_c);

    Omega_s.free();
    Omega_a.free();
    cerr <<  "components done" << "\n";

    begin = clock();
    for(int a=0; a<8; a++)
    {
        cerr <<  "cycle " << a  << "\n";
        fftOmega_a_c.at(a).resize(size_x,size_x);
        fftOmega_s_c.at(a).resize(size_x,size_x);

        blitz::Array<cd,2> tmp1(size_x,size_x);
        blitz::Array<cd,2> tmp2(size_x,size_x);

        tmp1 = Omega_a_c.at(a);
        FFTW(tmp1,tmp2);
        fftOmega_a_c.at(a)=tmp2*pow(step_x,2);

        tmp1 = Omega_s_c.at(a);
        FFTW(tmp1,tmp2);
        fftOmega_s_c.at(a)=tmp2*pow(step_x,2);

        cerr << fftOmega_a_c.at(a)(0,0) << " " << fftOmega_s_c.at(a)(0,0) << "\n" << flush;

        Omega_s_c.at(a).free();
        Omega_a_c.at(a).free();

    }
    end = clock();

    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cerr << "components done\n" << "time " << elapsed_secs   << endl << flush;

    double Qs2= pow(Qs(0),2);

   	double dK = 2*step_k;
    double DK = dK*2;
    double Kmax = 10+DK;
	int NK = (Kmax-DK-dK)/DK;  
	int nk;



	vector<double> KV;


	for(double k=0.5; k<25;  k+=0.05)
	{
		KV.push_back(k); 
	}

	/*
	KV.push_back(0.5); 
	KV.push_back(1); 
	KV.push_back(1.5); 
	KV.push_back(2); 
	KV.push_back(2.5); 
	KV.push_back(3); 
	KV.push_back(3.5); 
	KV.push_back(4); 
	KV.push_back(4.5); 
	KV.push_back(5); 
	KV.push_back(5.5); 
	KV.push_back(6); 
	KV.push_back(6.5); 
	KV.push_back(7); 
	KV.push_back(8); 
	KV.push_back(9); 
	KV.push_back(10); 
	KV.push_back(11); 
	KV.push_back(12); 
	KV.push_back(13); 
	KV.push_back(14); 
	KV.push_back(15); */


  	for(nk=0; nk<KV.size(); nk++)
    {
		KV.at(nk) = KV.at(nk)*Qs(0);  
	}

	blitz::Array<cd,2> SIP(size_x,size_x);
	blitz::Array<cd,2> SIP_p(size_x,size_x);
	blitz::Array<cd,2> SIP_A(size_x,size_x);

	blitz::Array<cd,2> SIPx(size_x,size_x);

	for(int i1=0; i1<size_x; i1=i1+1)
        {
            for(int j1=0; j1<size_x; j1=j1+1)
            {
 				double k12 =  int_to_k2(i1,j1);
				SIP(i1,j1) = SI(i1, j1, fftOmega_s_c, fftOmega_a_c); //populate the matrix for future back FFT

				SIP_p(i1,j1) = SIP(i1,j1);
				SIP_A(i1,j1) = SIP(i1,j1);
				
				if(k12<g2mu*g2mu) 
				{
					SIP_p(i1,j1) = 0.0; 
				}

				if(k12<1.0) 
				{
					SIP_A(i1,j1) = 0.0; 
				}
			}
		}



    FFTW_b(SIP, SIPx);
	double Mult = real(SIPx(0,0)); 

    FFTW_b(SIP_p, SIPx);
	double Mult_p = real(SIPx(0,0)); 

	FFTW_b(SIP_A, SIPx);
	double Mult_A = real(SIPx(0,0)); 
    
	for(nk=0; nk<KV.size(); nk++)
    {
		cerr << nk << "\n" << flush;
		
		double K = KV.at(nk); 
		cerr << K << "\n" << flush;
		double v0 = 0.0; 
		cd v1 = cd(0.0,0.0);	
		cd v2 = cd(0.0,0.0);
		cd v3 = cd(0.0,0.0);
		cd v4 = cd(0.0,0.0);
		cd v5 = cd(0.0,0.0);

		
		
		
		int N=0;

        for(int i1=0; i1<size_x; i1=i1+1)
        {
			if( (i1<size_x/4+1)||(i1>3*size_x/4-1) ){
            double kx1_t  = int_to_k(i1);

            for(int j1=0; j1<size_x; j1=j1+1)
            {
				if( (j1<size_x/4+1)||(j1>3*size_x/4-1) )
				{
                double ky1_t  = int_to_k(j1);
                //double k12 =  x2(kx1_t) + x2(ky1_t); //   int_to_k2(i1,j1);
                double k12 =  int_to_k2(i1,j1);

				double k1 = sqrt(k12); 
                if( (k1-(K-0.5*dK) ) * ( k1-(K+0.5*dK) )  <0.0 )
                {
                    cd amp = SI(i1, j1, fftOmega_s_c, fftOmega_a_c); 
                    cd amp3 = cd(0.0); 
						//AsSI(i1, j1, fftOmega_s_c, fftOmega_a_c); 
                    double phi_1 = atan2(ky1_t,kx1_t);

					v0+=real(amp)  ; 
					v1+= amp3  * exp( cd(0.0,1.0) * phi_1) ; 
					v2+= amp *  exp( cd(0.0,2.0) * phi_1) ; 
					v3+= amp3 * exp( cd(0.0,3.0) * phi_1) ; 
					v5+= amp3 * exp( cd(0.0,5.0) * phi_1) ; 

					std::cout << (ParallelStream() 
     				<< K << " " << phi_1 << " " << real(k1*dK*amp) << " " << real(k1*dK*amp3)
					<< "\n").toString() << flush;
									
					N++;


                 }
				}

            }
			}

     	}

		cout << "# " << K<< " "   
			<< v0 << " " 
			<< real(v1)/v0 << " " << imag(v1)/v0 << " " 
			<< real(v2)/v0 << " " << imag(v2)/v0 << " " 
			<< real(v3)/v0 << " " << imag(v3)/v0 << " "  
			<< real(v5)/v0 << " " << imag(v5)/v0 << " " 
			<< Mult <<  " "  
			<< Mult_p <<  " "  
			<< Mult_A <<  " "  
			<<
			std::endl;  
		
    }
}

















int main_old(void)
{
	clock_t begin = clock();
	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

cin >> eventID;

	cin >> A;
	cin >> B;
    dname = "v3_data"+eventID;
    mkdir(dname.c_str(),S_IRWXU | S_IRWXG);

    string name = dname+"/MD_" + eventID + ".dat";
    fileout.open(name.c_str());

    //Target:
    colorArr V_c(size_x,size_x);
    begin = clock();
    IC_MV(V_c);
    end = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cerr << "target done\n" << "time " << elapsed_secs << flush  << endl;

    //Projectile:
    //
    vector<blitz::Array<double,2> > A_a(8);
    begin = clock();
    IC_p(A_a);
    end = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cerr << "projectile done\n" << "time " << elapsed_secs << flush  << endl;

    begin = clock();
    colorArr Omega_a(size_x,size_x), Omega_s(size_x,size_x);
    vector<blitz::Array<cd,2> > Omega_a_c(8), Omega_s_c(8);
    vector<blitz::Array<cd,2> > fftOmega_a_c(8), fftOmega_s_c(8);

    Omega(Omega_s, Omega_a, V_c, A_a);
    end = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    cerr << "omega done\n" << "time " << elapsed_secs << flush <<  endl;

    //output(0.0, V_c);

    V_c.free();
    for (int a=0; a<8; a++) A_a.at(a).free();
    cerr << "output done\n" << endl;

    components(Omega_s, Omega_s_c);
    components(Omega_a, Omega_a_c);

    Omega_s.free();
    Omega_a.free();
    cerr <<  "components done" << "\n";

    begin = clock();
    for(int a=0; a<8; a++)
    {
        cerr <<  "cycle " << a  << "\n";
        fftOmega_a_c.at(a).resize(size_x,size_x);
        fftOmega_s_c.at(a).resize(size_x,size_x);

        blitz::Array<cd,2> tmp1(size_x,size_x);
        blitz::Array<cd,2> tmp2(size_x,size_x);

        tmp1 = Omega_a_c.at(a);
        FFTW(tmp1,tmp2);
        fftOmega_a_c.at(a)=tmp2*pow(step_x,2);

        tmp1 = Omega_s_c.at(a);
        FFTW(tmp1,tmp2);
        fftOmega_s_c.at(a)=tmp2*pow(step_x,2);

        cerr << fftOmega_a_c.at(a)(0,0) << " " << fftOmega_s_c.at(a)(0,0) << "\n" << flush;

        Omega_s_c.at(a).free();
        Omega_a_c.at(a).free();

    }
    end = clock();

    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cerr << "components done\n" << "time " << elapsed_secs   << endl << flush;

    double Qs2= pow(Qs(0),2);

   	double dK = 2*step_k;
    double DK = dK*2;
    double Kmax = 10+DK;
	int NK = (Kmax-DK-dK)/DK;  
	int nk;
	//#pragma  omp parallel for schedule(static,5)
    for(nk=0; nk<NK; nk++)
    {
		cerr << nk << "\n" << flush;
		
		double K = nk*DK+DK+dK; 
		cerr << K << "\n" << flush;
		double v0 = 0.0; 
		cd v1 = cd(0.0,0.0);	
		cd v2 = cd(0.0,0.0);
		cd v3 = cd(0.0,0.0);
		cd v4 = cd(0.0,0.0);
		cd v5 = cd(0.0,0.0);

        int N=0;

        for(int i1=0; i1<size_x; i1=i1+1)
        {
			if( (i1<size_x/4+1)||(i1>3*size_x/4-1) ){
            double kx1_t  = int_to_k(i1);

            for(int j1=0; j1<size_x; j1=j1+1)
            {
				if( (j1<size_x/4+1)||(j1>3*size_x/4-1) )
				{
                double ky1_t  = int_to_k(j1);
                //double k12 =  x2(kx1_t) + x2(ky1_t); //   int_to_k2(i1,j1);
                double k12 =  int_to_k2(i1,j1);

				double k1 = sqrt(k12); 
                if( (k1-(K-0.5*dK) ) * ( k1-(K+0.5*dK) )  <0.0 )
                {
                    cd amp = SI(i1, j1, fftOmega_s_c, fftOmega_a_c); 
                    cd amp3 = AsSI(i1, j1, fftOmega_s_c, fftOmega_a_c); 
                    double phi_1 = atan2(ky1_t,kx1_t);

					std::cout << (ParallelStream() 
     				<< K << " " << phi_1 << " " << real(k1*dK*amp) << " " << real(k1*dK*amp3)
					<< "\n").toString() << flush;
									
					N++;


                 }
				}

            }
			}

     	}
		
    }
}

