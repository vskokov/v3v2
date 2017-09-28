#include "JIMWLK.h"
#include "GM.h"
#include "string_conv.h"

using namespace std;
using namespace blitz;



//choice of parameters g^2 \mu = 1

const int size_x=512; // # steps in transverse direction
const int N_Y=100; // # steps in long direction

const double step_x=0.05;
const double L_x=step_x*size_x; // transverse extent
//const double step_x=L_x/size_x; // transvers step size


const double step_x2 = step_x*step_x;
const int size_x2 = size_x*size_x;

const double dY_times_alpha=0.001;

const double prefactor = sqrt(dY_times_alpha)/M_PI;
double UV = step_x*step_x*1e-5;

string eventID;

boost::normal_distribution<> nd( 0.0, 1.0 / (step_x * sqrt(N_Y) ) );
boost::variate_generator< RNGType, boost::normal_distribution<> > Gaussian(rng, nd);


boost::normal_distribution<> nd_zeta( 0.0, 1.0 / step_x  );
boost::variate_generator< RNGType, boost::normal_distribution<> > Gaussian_zeta(rng, nd_zeta);

typedef blitz::TinyVector<cd,8> AdComp;
typedef blitz::TinyMatrix<cd,3,3>  colorMat;
typedef blitz::TinyMatrix<double,8,8>  colorAdMat;
typedef Array<colorMat,2>  colorArr;
typedef Array<colorAdMat,2>  colorAdArr;

ofstream fileout;


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
}

double int_to_x(int i)
{
    return i*step_x;
}


double Qs(double Y)
{
	double A0 = 0.269933;
	double A1 = 1.60715;
	double A2 = 0.430101;
	return A0*exp(A1*Y)+A2;
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

Array<cd,2> f_FT_Kx(void)
{
    Array<cd,2> Kx(size_x,size_x);
    Array<cd,2> FT_Kx(size_x,size_x);

    for(int i=0; i<size_x; i++)
    {
        for(int j=0; j<size_x; j++)
        {
            double x = int_to_x(i);
            double y = int_to_x(j);

            /*
            int nx_prime=0;
            int ny_prime=0;
            double distance_prime=x*x+y*y;
            for(int nx=-1;nx<2;nx++)
            for(int ny=-1;ny<2;ny++)
            {
            	double distance = (x+L_x*nx)*(x+L_x*nx)+(y+L_x*ny)*(y+L_x*ny);
            	if(distance<=distance_prime)
            	{
            		distance_prime=distance;
            		nx_prime=nx;
            		ny_prime=ny;
            	}
            }

            double x_p = x + L_x*nx_prime;
            double y_p = y + L_x*ny_prime;

            Kx(i,j) = x_p/(x_p*x_p+y_p*y_p+UV);
            */
            double denominator = ( x2(sin(0.5*x/L_d_by_2pi)) + x2(sin(0.5*y/L_d_by_2pi)) + 1e-20 ) * x2(2.0*L_d_by_2pi );
            Kx(i,j) = L_d_by_2pi * sin (x/L_d_by_2pi ) / denominator;
        }
    }

    FFTW(Kx,FT_Kx);
    return FT_Kx;
}

Array<cd,2> f_FT_Ky(void)
{
    Array<cd,2> Ky(size_x,size_x);
    Array<cd,2> FT_Ky(size_x,size_x);

    for(int i=0; i<size_x; i++)
    {
        for(int j=0; j<size_x; j++)
        {
            double x = int_to_x(i);
            double y = int_to_x(j);


            /*int nx_prime=0;
            int ny_prime=0;
            double distance_prime=x*x+y*y;
            for(int nx=-1;nx<2;nx++)
            for(int ny=-1;ny<2;ny++)
            {
            	double distance = (x+L_x*nx)*(x+L_x*nx)+(y+L_x*ny)*(y+L_x*ny);
            	if(distance<=distance_prime)
            	{
            		distance_prime=distance;
            		nx_prime=nx;
            		ny_prime=ny;
            	}
            }

            double x_p = x + L_x*nx_prime;
            double y_p = y + L_x*ny_prime;

            Ky(i,j) = y_p/(x_p*x_p+y_p*y_p+UV);
            */
            double denominator = ( x2(sin(0.5*x/L_d_by_2pi)) + x2(sin(0.5*y/L_d_by_2pi)) + 1e-20 ) * x2(2.0*L_d_by_2pi );
            Ky(i,j) = L_d_by_2pi * sin (y/L_d_by_2pi ) / denominator;

        }
    }

    FFTW(Ky,FT_Ky);
    return FT_Ky;
}




void generate_zeta0(blitz::Array<double,2>& zeta0 )
{

    for(int i=0; i<size_x; i++)
    {
        //cout << i << "\n" << flush;
        for(int j=0; j<size_x; j++)
        {
            zeta0(i,j) = Gaussian_zeta();
        }
    }

    // Periodic boundary conditions
    //for(int i=0; i<size_x; i++) zeta0(i,0)=zeta0(i,size_x-1);
    //for(int j=0; j<size_x; j++) zeta0(0,j)=zeta0(size_x-1,j);
}



void generate_zeta(colorArr& zeta)
{
    blitz::Array<double,2> zeta0[8];
    vector < TinyMatrix<cd,3,3>  > lambda = get_lambda();

    for(int a=0; a<8; a++)
    {
        zeta0[a].resize(size_x,size_x);
        generate_zeta0(zeta0[a]);
    }


    for(int i=0; i<size_x; i++)
    {
        for(int j=0; j<size_x; j++)
        {

            colorMat zeta_sum;
            zeta_sum=0;

            for(int a=0; a<8; a++)
            {
                zeta_sum += zeta0[a](i,j) * lambda.at(a);
            }

            zeta(i,j) = zeta_sum;
        }
    }
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




colorArr Integral1 (colorArr& V, colorArr& zeta_x, colorArr& zeta_y )
{
    colorArr Int(size_x,size_x);

    Array<cd,2> helper_out(size_x,size_x);

    Array<cd,2> helper_x_in(size_x,size_x);
    Array<cd,2> helper_x_out(size_x,size_x);

    Array<cd,2> helper_y_in(size_x,size_x);
    Array<cd,2> helper_y_out(size_x,size_x);



    colorArr V_zetax_Vdagger(size_x,size_x);
    colorArr V_zetay_Vdagger(size_x,size_x);



    Array<cd,2> FT_Kx(size_x,size_x);
    Array<cd,2> FT_Ky(size_x,size_x);

    FT_Kx = f_FT_Kx();
    FT_Ky = f_FT_Ky();

    for(int i=0; i<size_x; i++)
    {
        for(int j=0; j<size_x; j++)
        {
            colorMat Vdag = dagger(V(i,j));
            V_zetax_Vdagger(i,j)  =    Product(V(i,j) , Product( zeta_x(i,j) , Vdag) );
            V_zetay_Vdagger(i,j)  =    Product(V(i,j) , Product( zeta_y(i,j) , Vdag) );
        }
    }

    Array <cd,2> FT_K_V_z_Vd(size_x,size_x);

    for(int a=0; a<3; a++)
        for(int b=0; b<3; b++)
        {
            for(int i=0; i<size_x; i++)
            {
                for(int j=0; j<size_x; j++)
                {
                    helper_x_in(i,j) = 	V_zetax_Vdagger(i,j)(a,b);
                    helper_y_in(i,j) = 	V_zetay_Vdagger(i,j)(a,b);
                }
            }

            FFTW(helper_x_in,helper_x_out);
            FFTW(helper_y_in,helper_y_out);

            for(int i=0; i<size_x; i++)
            {
                for(int j=0; j<size_x; j++)
                {
                    //double k2 = k_x*k_x+k_y*k_y+1e-6;

                    FT_K_V_z_Vd(i,j) = (helper_x_out(i,j) * FT_Kx(i,j)+  helper_y_out(i,j) * FT_Ky(i,j));
                }
            }


            FFTW_b(FT_K_V_z_Vd,helper_out); //should be normalized by size_x^2

            for(int i=0; i<size_x; i++)
            {
                for(int j=0; j<size_x; j++)
                {
                    Int(i,j)(a,b) = (helper_out(i,j)) * step_x2/double(size_x2) ;
                }
            }


        }
    return Int;
}



colorArr Integral2 (colorArr& zeta_x, colorArr& zeta_y )
{
    colorArr Int(size_x,size_x);

    Array<cd,2> helper_out(size_x,size_x);

    Array<cd,2> helper_x_in(size_x,size_x);
    Array<cd,2> helper_x_out(size_x,size_x);

    Array<cd,2> helper_y_in(size_x,size_x);
    Array<cd,2> helper_y_out(size_x,size_x);


    Array<cd,2> FT_Kx(size_x,size_x);
    Array<cd,2> FT_Ky(size_x,size_x);

    FT_Kx = f_FT_Kx();
    FT_Ky = f_FT_Ky();


    Array <cd,2> FT_K_z(size_x,size_x);

    for(int a=0; a<3; a++)
        for(int b=0; b<3; b++)
        {
            for(int i=0; i<size_x; i++)
            {
                for(int j=0; j<size_x; j++)
                {
                    helper_x_in(i,j) = 	zeta_x(i,j)(a,b);
                    helper_y_in(i,j) = 	zeta_y(i,j)(a,b);
                }
            }

            FFTW(helper_x_in,helper_x_out);
            FFTW(helper_y_in,helper_y_out);

            for(int i=0; i<size_x; i++)
            {
                for(int j=0; j<size_x; j++)
                {
                    FT_K_z(i,j) = (helper_x_out(i,j) * FT_Kx(i,j) +  helper_y_out(i,j) * FT_Ky(i,j)) ;
                }
            }


            FFTW_b(FT_K_z,helper_out); //should be normalized by size_x^2

            for(int i=0; i<size_x; i++)
            {
                for(int j=0; j<size_x; j++)
                {
                    Int(i,j)(a,b) = (helper_out(i,j))  * step_x2/double(size_x2)  ;
                }
            }


        }
    return Int;
}




colorMat Evolution_kernel(int i, int j, colorArr& V, colorArr& Int1, colorArr& Int2)
{
    colorMat exp1;
    colorMat u_exp1;
    colorMat exp2;
    colorMat u_exp2;


    u_exp1 = -prefactor*Int1(i,j);
    u_exp2 = prefactor*Int2(i,j); //i's are in the matrix exponent!
    exp1 =  matrix_exp_an_antiH(u_exp1);
    exp2 =  matrix_exp_an_antiH(u_exp2);

    colorMat out;
    out = Product( Product( exp1,  V(i,j) ) , exp2   );
    return out;
}



void evo_step(colorArr &V_prev, colorArr &V_next)
{


    colorArr zeta_x(size_x,size_x);
    colorArr zeta_y(size_x,size_x);

    generate_zeta(zeta_x);
    generate_zeta(zeta_y);


    colorArr Int1(size_x,size_x);
    colorArr Int2(size_x,size_x);
    Int1 = Integral1(V_prev, zeta_x, zeta_y);
    Int2 = Integral2(zeta_x, zeta_y);

    for(int i=0; i<size_x; i++)
        for(int j=0; j<size_x; j++)
        {

            V_next(i,j) = Evolution_kernel(i,j,V_prev, Int1, Int2);
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

    // Periodic boundary conditions
    for(int i=0; i<size_x; i++) rho(i,0)=rho(i,size_x-1);
    for(int j=0; j<size_x; j++) rho(0,j)=rho(size_x-1,j);

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
    if(i<0) i = i + size_x - 1;
    if(i>size_x-1) i = i - (size_x-1);
    return i;
}


int index_y_boundary(int in)
{
    //periodic boundaries
    int i = in;
    if(i<0) i = i + size_x - 1;
    if(i>size_x-1) i = i - (size_x-1);
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

	for (int i=0;i<max_k;i++) 
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

				double distance2 = ( x2(sin(0.5*x/L_d_by_2pi)) + x2(sin(0.5*y/L_d_by_2pi))  ) * x2(2.0*L_d_by_2pi );
				//if(distance2*Qs_Y*Qs_Y<1) comp(i,j) = su3_group_element(V_c(i,j), a);
				comp(i,j) = su3_group_element(V_c(i,j), a) * exp(-0.5*distance2*Qs_Y*Qs_Y);
            }
        }
        FFTW(comp, compFT);
        for(int i=0; i<size_x; i=i+1)
        {        
			double kx  = 2.0*M_PI*i/L_x;
        	double kx_t  = 1.0/step_x*sin(kx*step_x);
            for(int j=0; j<size_x; j=j+1)
            {

            	double ky  = 2.0*M_PI*j/L_x;
            	double ky_t  = 1.0/step_x*sin(ky*step_x);
				double k2 = (kx_t*kx_t+ky_t*ky_t);

				sum_all(i,j) += compFT(i,j)*conj(compFT(i,j));

                if(k2>pow(Qs_Y,2))
				{
					sum(i,j) += compFT(i,j)*conj(compFT(i,j));
                	sumqPerp(i,j) += compFT(i,j)*conj(compFT(i,j)) * k2 ;
				}

			}
        }
    }


	//Flatten 
	for(int i=0; i<size_x; i=i+1)
    {        
		double kx  = 2.0*M_PI*i/L_x;
        double kx_t  = 1.0/step_x*sin(kx*step_x);
        for(int j=0; j<size_x; j=j+1)
        {

            double ky  = 2.0*M_PI*j/L_x;
            double ky_t  = 1.0/step_x*sin(ky*step_x);
			double k2 = (kx_t*kx_t+ky_t*ky_t);


			int ik =  int(sqrt(k2) / step_k + 0.5); 
			if (ik < max_k_int) 
			{
				distr.at(ik) += real(sum_all(i,j))*0.5/3.0/size_x2; 
				normal.at(ik) += 1.0; 
			}
			if( sqrt(k2) < 15*Qs(Y) ) fileout << Y << " " << kx_t/Qs(Y) << " " << ky_t/Qs(Y) << " " 
			    << real(sum_all(i,j))*0.5/3.0/size_x2 << " " << Qs(Y) <<"\n" << flush;

		}
		fileout << "\n"; 
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


    FFTW_b(sum, pS);
    FFTW_b(sumqPerp, pSPerp);

/*
    for(int i=0; i<size_x; i++)
    {
        for(int j=0; j<size_x; j++)
        {
           d_data <<  sqrt(x2(sin(int_to_x(i)*M_PI/L_x)) + x2(sin(int_to_x(j)*M_PI/L_x)))*L_x/M_PI  << " " << 1.0-real(pS (i,j))*0.5/3.0/size_x2/size_x2 << "\n" ;
        }
    }
	*/


    //d_data<< real(pS(0,0)) << " " <<  real(pSPerp(0,0)) << "\n" <<  flush;
    cout << Y << " " << real(pS(0,0))*0.5/3.0/size_x2/size_x2 << " " <<  real(pSPerp(0,0)) *0.5/3.0/size_x2/size_x2<< " "<< Qs(Y) <<  "\n" <<    flush;



}

void TMD(double Y, colorArr& V_c)
{

    Array<AdComp,2>  Ax_field(size_x,size_x);
    Array<AdComp,2>  Ay_field(size_x,size_x);

    Array<AdComp,2>  Akx_field(size_x,size_x);
    Array<AdComp,2>  Aky_field(size_x,size_x);

    Array<double,2>  TMD_xx_field(size_x,size_x);
    Array<double,2>  TMD_xy_field(size_x,size_x);
    Array<double,2>  TMD_yy_field(size_x,size_x);

    colorAdArr Um(size_x,size_x);
    AdWilsonLine(Um, V_c);
    for(int i=0; i<size_x; i=i+1)
    {
        for(int j=0; j<size_x; j=j+1)
        {

            colorAdMat Udagger;
            colorAdMat dUx, dUy;
            colorAdMat Ax, Ay;

            Udagger = dagger(Um(i,j));

            //Derivatives
            dUx = (Um(index_x_boundary(i+1), j) - Um(index_x_boundary(i-1) , j))/step_x*0.5;
            dUy = (Um(i, index_x_boundary(j+1)) - Um(i, index_x_boundary(j-1)))/step_x*0.5;

            Ax = Product(Udagger, dUx);
            Ay = Product(Udagger, dUy);

            AdComp components_x, components_y;
            components_x = DecomposeAdMatrix(Ax);
            components_y = DecomposeAdMatrix(Ay);
            Ax_field(i,j) = components_x;
            Ay_field(i,j) = components_y;
        }
    }

    //FFT

    for(int a=0; a<8; a++)
    {
        Array<cd,2> tmp(size_x, size_x);
        Array<cd,2> tmpk(size_x, size_x);

        //x component
        for(int i=0; i<size_x; i=i+1)
        {
            for(int j=0; j<size_x; j=j+1)
            {
                tmp(i,j) = Ax_field(i,j)(a);
            }
        }
        FFTW(tmp, tmpk);
        for(int i=0; i<size_x; i=i+1)
        {
            for(int j=0; j<size_x; j=j+1)
            {
                Akx_field(i,j)(a) = tmpk(i,j)*step_x2/(4.0*M_PI*M_PI);
            }
        }

        //y component
        for(int i=0; i<size_x; i=i+1)
        {
            for(int j=0; j<size_x; j=j+1)
            {
                tmp(i,j) = Ay_field(i,j)(a);
            }
        }
        FFTW(tmp, tmpk);
        for(int i=0; i<size_x; i=i+1)
        {
            for(int j=0; j<size_x; j=j+1)
            {
                Aky_field(i,j)(a) = tmpk(i,j)*step_x2/(4.0*M_PI*M_PI);
            }
        }
    }

    // xx componetnt
    for(int i=0; i<size_x; i=i+1)
    {
        for(int j=0; j<size_x; j=j+1)
        {
            // xx componetnt
            double Trace_sum =0;
            for(int a=0; a<8; a++)
            {
                Trace_sum+=real(Akx_field(i,j)(a)*conj(Akx_field(i,j)(a)));
            }
            TMD_xx_field(i,j) = -3.0*Trace_sum;

            // xy componetnt

            Trace_sum =0;
            for(int a=0; a<8; a++)
            {
                Trace_sum+=real(Akx_field(i,j)(a)*conj(Aky_field(i,j)(a)));
            }
            TMD_xy_field(i,j) = -3.0*Trace_sum;

            // yy componetnt

            Trace_sum =0;
            for(int a=0; a<8; a++)
            {
                Trace_sum+=real(Aky_field(i,j)(a)*conj(Aky_field(i,j)(a)));
            }
            TMD_yy_field(i,j) = -3.0*Trace_sum;
        }
    }

    ofstream t_data;
    string fTMDname = dname+"/TMD_k_" + toString(Y) + ".dat";
    t_data.open(fTMDname.c_str());

    for(int i=0; i<size_x; i=i+1)
    {
        double kx  = 2.0*M_PI*i/L_x;
        double kx_t  = 1.0/step_x*sin(kx*step_x);

        for(int j=0; j<size_x; j=j+1)
        {
            double ky  = 2.0*M_PI*j/L_x;
            double ky_t  = 1.0/step_x*sin(ky*step_x);

            double TMD_G =
                TMD_xx_field(i,j) +
                TMD_yy_field(i,j);
            double k2 = kx*kx+ky*ky;
            double k2_t = kx_t*kx_t+ky_t*ky_t;
            double TMD_H=TMD_G;
            if(k2>0.0) TMD_H+=
                    - 2.0*(
                        kx_t*kx_t/k2_t * TMD_xx_field(i,j)+
                        2.0*kx_t*ky_t/k2_t * TMD_xy_field(i,j)+
                        ky_t*ky_t/k2_t * TMD_yy_field(i,j)) ;
            t_data
                    << kx_t << " "
                    << ky_t << " "
                    << TMD_H << " "
                    << TMD_G << "\n";
        }
        t_data<<"\n";
    }

    t_data.close();
}



void binner_effective_action_FT(double Y, colorArr& V_c)
{

	ofstream fileoutFT;
    string nameFT = dname+"/eaFT_" + toString(Y) + ".dat";
    fileoutFT.open(nameFT.c_str());

    int NbinD = 150;
    int NbinK = 120;

    Array<double,2>  bin(NbinK+1, NbinD+1);
    bin=0.0;

    double dMax = 1.2;
    double dMin = -0.5;

	double Kcuttoff=15.0;
    double KMax = Kcuttoff;
    double KMin = 0;

    double stepD = (dMax-dMin)/NbinD;
    double stepK = (KMax-KMin)/NbinK;

	blitz::Array<cd,2>  comp(size_x, size_x);
    blitz::Array<cd,2>  compFT(size_x, size_x);
    blitz::Array<cd,2>  sum(size_x, size_x);
    blitz::Array<cd,2>  pS(size_x, size_x);
    sum=cd(0.0,0.0);

    for(int a=0; a<9; a++)
    {
        for(int i=0; i<size_x; i=i+1)
        {
            for(int j=0; j<size_x; j=j+1)
            {
                comp(i,j) = su3_group_element(V_c(i,j), a);
            }
        }
        FFTW(comp, compFT);
        for(int i=0; i<size_x; i=i+1)
        {
            for(int j=0; j<size_x; j=j+1)
            {
                sum(i,j) += compFT(i,j)*conj(compFT(i,j))*step_x2/(4.0*M_PI*M_PI)/double(size_x2);
            }
        }
    }

	for(int i=0; i<size_x/2; i++)
    {
   double kx  = 2.0*M_PI*i/L_x;
        double kx_t  = 1.0/step_x*sin(kx*step_x);

        for(int j=0; j<size_x/2; j++)
		{
	            double ky  = 2.0*M_PI*j/L_x;
            double ky_t  = 1.0/step_x*sin(ky*step_x);
			double k  = sqrt( kx_t*kx_t + ky_t*ky_t);
			double S = real(sum(i,j))/6.0;
			//cout << S << "\n";
			 if((k<KMax)&&(S<dMax)&&(S>dMin))
                            {
                                int iK = int((k - KMin)/stepK +0.5) ;
                                int iD = int((S - dMin)/stepD +0.5) ;
                                bin(iK,iD) += 1.0;
                            }

		}
	}


    //Normalization
    for(int ik = 0; ik < NbinK+1; ik++)
    {
        double sum = 0.0;
        for(int iD = 0; iD < NbinD+1; iD++)
        {
            sum = sum + bin(ik, iD);
        }

        for(int iD = 0; iD < NbinD+1; iD++)
        {
            bin(ik, iD) = bin(ik, iD)/sum;
        }
    }

    for(int ik=0; ik < NbinK+1; ik++)
    {
        for(int iD=0; iD < NbinD+1; iD++)
        {
            fileoutFT  << ik * stepK + KMin << " " << iD * stepD + dMin   << " "  << bin(ik,iD) << " "<< Y <<   "\n" <<flush;
        }
        fileoutFT << "\n";
    }
    fileoutFT.close();

}



void binner_effective_action_WW(double Y, colorArr& V_c)
{

    ofstream fileout;
    string name = dname+"/G_" + toString(Y) + ".dat";
    fileout.open(name.c_str());

    int NbinD = 150;
    int NbinR = 100;

    Array<double,2>  bin(NbinR+1, NbinD+1);
    Array<double,2>  bin1d(NbinR+1, NbinD+1);
    bin=0.0;
    bin1d=0.0;

    double dMax = 2.5;
    double dMin = -2.5;

	double rcuttoff=5.0;
    double rMax = rcuttoff/L_x;
    double rMin = 0;

    double stepD = (dMax-dMin)/NbinD;
    double stepR = (rMax-rMin)/NbinR;

    #pragma omp parallel
    #pragma omp for
    for(int ix=0; ix<size_x; ix++)
    {
        cout << "outter loop effective action " << ix << endl<<flush;
        for(int iy=0; iy<size_x; iy++)
        {
            if(fabs(int_to_x(ix)-int_to_x(iy))<rcuttoff)
            {
				for(int jx=0; jx<size_x; jx++)
                {
                    for(int jy=0; jy<size_x; jy++)
                    {

                        double r  =
                            sqrt(
                                pow(int_to_x(ix)-int_to_x(iy),2)
                                +
                                pow(int_to_x(jx)-int_to_x(jy),2)
                            );

                        if(r<rcuttoff)
                        {

							colorMat dVx_1, dVy_1; //component 1 for the points x and y  
							colorMat dVx_2, dVy_2; //component 2 for the points x and y  
								
							colorMat Ax_1, Ax_2; 
							colorMat Ay_1, Ay_2; 

                            colorMat  Vdaggerx, Vdaggery, summed;
                            
							Vdaggerx = dagger(V_c(ix,jx));
							Vdaggery = dagger(V_c(iy,jy));

            				//Derivatives
            				dVx_1 = (V_c(index_x_boundary(ix+1), jx) - V_c(index_x_boundary(ix-1) , jx))/step_x*0.5;
            				dVy_1 = (V_c(index_x_boundary(iy+1), jy) - V_c(index_x_boundary(iy-1) , jy))/step_x*0.5;
            				
							dVx_2 = (V_c(ix, index_x_boundary(jx+1)) - V_c(ix, index_x_boundary(jx-1)))/step_x*0.5;
							dVy_2 = (V_c(iy, index_x_boundary(jy+1)) - V_c(iy, index_x_boundary(jy-1)))/step_x*0.5;

            				Ax_1 = Product(Vdaggerx, dVx_1);
            				Ay_1 = Product(Vdaggery, dVy_1);

							Ax_2 = Product(Vdaggerx, dVx_2);
            				Ay_2 = Product(Vdaggery, dVy_2);

                            summed=Product(Ax_1, Ay_1) + Product(Ax_2,Ay_2);

                            double G = real(trace(summed))/3.0;

                            r = r / L_x;

                            if((r<rMax)&&(G<dMax)&&(G>dMin))
                            {
                                int ir = int((r - rMin)/stepR +0.5) ;
                                int iD = int((G - dMin)/stepD +0.5) ;
                                bin(ir,iD) += 1.0;
                            }
                        }
                    }
                }
            }
        }
    }

    //Normalization
    for(int ir = 0; ir < NbinR+1; ir++)
    {
        double sum = 0.0;
        for(int iD = 0; iD < NbinD+1; iD++)
        {
            sum = sum + bin(ir, iD);
        }

        for(int iD = 0; iD < NbinD+1; iD++)
        {
            bin(ir, iD) = bin(ir, iD)/sum;
        }
    }

    for(int ir=0; ir < NbinR+1; ir++)
    {
        for(int iD=0; iD < NbinD+1; iD++)
        {
            fileout  << ir * stepR + rMin << " " << iD * stepD + dMin   << " "  << bin(ir,iD) << " "<< Y <<   "\n" <<flush;
        }
        fileout << "\n";
    }
    fileout.close();
}





void binner_effective_action(double Y, colorArr& V_c)
{

    ofstream fileout;
    string name = dname+"/ea_" + toString(Y) + ".dat";
    fileout.open(name.c_str());

	ofstream fileoutFT;
    string nameFT = dname+"/eaFT_" + toString(Y) + ".dat";
    fileoutFT.open(nameFT.c_str());

    ofstream fileout1d;
    string name1d = dname+"/ea1d_" + toString(Y) + ".dat";
    fileout1d.open(name1d.c_str());

    int NbinD = 150;
    int NbinR = 100;

    Array<double,2>  bin(NbinR+1, NbinD+1);
    Array<double,2>  bin1d(NbinR+1, NbinD+1);
    bin=0.0;
    bin1d=0.0;

    double dMax = 1;
    double dMin = -0.5;

	double rcuttoff=5.0;
    double rMax = rcuttoff/L_x;
    double rMin = 0;

    double stepD = (dMax-dMin)/NbinD;
    double stepR = (rMax-rMin)/NbinR;

    #pragma omp parallel
    #pragma omp for
    for(int ix=0; ix<size_x; ix++)
    {
        cout << "outter loop effective action " << ix << endl<<flush;
        for(int iy=0; iy<size_x; iy++)
        {
            if(fabs(int_to_x(ix)-int_to_x(iy))<rcuttoff)
            {
				for(int jx=0; jx<size_x; jx++)
                {
                    for(int jy=0; jy<size_x; jy++)
                    {

                        double r  =
                            sqrt(
                                pow(int_to_x(ix)-int_to_x(iy),2)
                                +
                                pow(int_to_x(jx)-int_to_x(jy),2)
                            );

                        if(r<rcuttoff)
                        {
                            colorMat V;
                            colorMat  Vdagger;
                            V = V_c(ix,jx);
                            Vdagger = dagger(V_c(iy,jy));

                            colorMat  dipole;
                            dipole=Product(Vdagger, V);

                            double S = real(trace(dipole))/3.0;

                            r = r / L_x;

                            if((r<rMax)&&(S<dMax)&&(S>dMin))
                            {
                                int ir = int((r - rMin)/stepR +0.5) ;
                                int iD = int((S - dMin)/stepD +0.5) ;
                                bin(ir,iD) += 1.0;
                            }
                        }
                    }
                }
            }
        }
    }
	//1D
	int iR = 50; 
    double R  = int_to_x(iR) ;  

	 for(int i=0; i<size_x; i++)
    {
        for(int j=0; j<size_x; j++)
        {

						 colorMat V;
                          colorMat  Vdagger;
                            V = V_c(i,j);
                            Vdagger = dagger(V_c(index_x_boundary(i-iR),j));

                            colorMat  dipole;
                            dipole=Product(Vdagger, V);

                            double S = real(trace(dipole))/3.0;

                            if((S<dMax)&&(S>dMin))
                            {
                                int iD = int((S - dMin)/stepD +0.5) ;
                                bin1d(iD) += 1.0;
                            }
                        }
                    }



    //Normalization
    for(int ir = 0; ir < NbinR+1; ir++)
    {
        double sum = 0.0;
        for(int iD = 0; iD < NbinD+1; iD++)
        {
            sum = sum + bin(ir, iD);
        }

        for(int iD = 0; iD < NbinD+1; iD++)
        {
            bin(ir, iD) = bin(ir, iD)/sum;
        }
    }

 double sum1d = 0.0;
        for(int iD = 0; iD < NbinD+1; iD++)
        {
            sum1d = sum1d + bin1d(iD);
        }

        for(int iD = 0; iD < NbinD+1; iD++)
        {
            bin1d(iD) = bin1d(iD)/sum1d;
        }



    for(int ir=0; ir < NbinR+1; ir++)
    {
        for(int iD=0; iD < NbinD+1; iD++)
        {
            fileout  << ir * stepR + rMin << " " << iD * stepD + dMin   << " "  << bin(ir,iD) << " "<< Y <<   "\n" <<flush;
        }
    }
    fileout.close();


		for(int iD=0; iD < NbinD+1; iD++)
        {
            fileout1d<<  iD * stepD + dMin   << " "  << bin1d(iD) << " "<< Y << " " << R <<   "\n" <<flush;
        }
        fileout1d.close();

}



int main(void)
{
    cin >> eventID;
    dname = "TMD_data_"+eventID;
    mkdir(dname.c_str(),S_IRWXU | S_IRWXG);

    string name = "/efs/MD_" + toString(eventID) + ".dat";
    fileout.open(name.c_str());


    colorArr V_c(size_x,size_x);
    colorArr V_n(size_x,size_x);

    IC_MV(V_c);
    output(0, V_c);
    //TMD(0,V_c);
    //binner_effective_action (0,V_c);
    //binner_effective_action_WW (0,V_c);

    int i=0;
    for(double Y=0; Y<1.101; Y+=dY_times_alpha)
    {
        evo_step(V_c, V_n);
        V_c = V_n;
        //cout << Y+dY_times_alpha << "\n" << flush;
        i++;
        if(i==100)
        {
            output(Y+dY_times_alpha, V_c);
            //TMD(Y+dY_times_alpha, V_c);
           //binner_effective_action(Y+dY_times_alpha,V_c);
           //binner_effective_action_WW(Y+dY_times_alpha,V_c);
            //binner_effective_action_FT  (Y+dY_times_alpha,V_c);
            i=0;
        }
    }
}
