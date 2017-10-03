#include "JIMWLK.h"
#include "string_conv.h"
#include <map>

using namespace std;
using namespace blitz;



//choice of parameters g^2 \mu = 1

const int size_x=512*2; // # steps in transverse direction
const int N_Y=100; // # steps in long direction

const double L_x=32;// step_x*size_x; // transverse extent
const double step_x=L_x/size_x; // transvers step size


const double step_x2 = step_x*step_x;
const int size_x2 = size_x*size_x;

double UV = step_x*step_x*1e-5;

string eventID;

boost::normal_distribution<> nd( 0.0, 1.0 / (step_x * sqrt((double)N_Y) ) );
boost::variate_generator< RNGType, boost::normal_distribution<> > Gaussian(rng, nd);

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
        for(int b=0; b<8; b++) {
            for(int i=0; i<size_x; i++)
                for(int j=0; j<size_x; j++) {
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
    for(int i=0; i<size_x; i++) {
        for(int j=0; j<size_x; j++) {
            rho(i,j) = Gaussian();
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
    for(int ny=0; ny<N_Y; ny++) {

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
            blitz::TinyMatrix<cd,3,3>  V_unit;
            V_unit=cd(1.0,0.0),cd(0.0,0.0),cd(0.0,0.0),
            cd(0.0,0.0),cd(1.0,0.0),cd(0.0,0.0),
            cd(0.0,0.0),cd(0.0,0.0),cd(1.0,0.0);

			for(int i=0; i<size_x; i=i+1) {
        for(int j=0; j<size_x; j=j+1) {


            V(i,j)=V_unit;
        }
    }

    for(int ny=0; ny<N_Y; ny++) {

        cerr << "target slices" << ny << "\n";

        vector<blitz::Array<double,2> >  A_a(8);


        for(int i=0; i<8; i++)
            A_a.at(i).resize(size_x, size_x);

        for(int i=0; i<8; i++) {
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

        for(int i=0; i<size_x; i=i+1) {
            for(int j=0; j<size_x; j=j+1) {
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

    for (int i=0; i<max_k; i++) {
        distr.at(i) = 0.0;
        normal.at(i) = 0.0;
    }

    for(int a=0; a<9; a++) {
        for(int i=0; i<size_x; i=i+1) {
            for(int j=0; j<size_x; j=j+1) {
                double x = int_to_x(i);
                double y = int_to_x(j);

                //double distance2 = ( x2(sin(x/L_d_by_2pi)) + x2(sin(y/L_d_by_2pi))  ) * x2(L_d_by_2pi );
                //if(distance2*Qs_Y*Qs_Y<1) comp(i,j) = su3_group_element(V_c(i,j), a);
                comp(i,j) = su3_group_element(V_c(i,j), a) ;//* exp(-0.5*distance2*Qs_Y*Qs_Y);
            }
        }

        FFTW(comp, compFT);

        for(int i=0; i<size_x; i=i+1) {
            double kx  = 2.0*M_PI*i/L_x;
            double kx_t  = 2.0/step_x*sin(kx*step_x/2.0);
            for(int j=0; j<size_x; j=j+1) {

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

    for (int i=0; i<N; i++) {
        aS[i]=0.0;
        aN[i]=0.0;
    }

    for(int i=0; i<size_x; i++) {
        for(int j=0; j<size_x; j++) {
            double r = sqrt(x2(int_to_x(i)) + x2(int_to_x(j)));
            int ir = int(r/dr);
            if(ir<N) {
                aS[ir] += 1.0-real(pS (i,j))*0.5/3.0/size_x2/size_x2 ;
                aN[ir] += 1.0;
            }
            //fileout <<  sqrt(x2(int_to_x(i)) + x2(int_to_x(j)))  << " " << 1.0-real(pS (i,j))*0.5/3.0/size_x2/size_x2 << "\n" << flush;
        }
    }


    for(int ir=0; ir<N; ir++) {
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
        for(int j=0; j<size_x; j++) {


            Omega_s(i,j) = cd(0.0,0.0);
            Omega_a(i,j) = cd(0.0,0.0);


            for(int a=0; a<8; a++) {
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



double SI(const int ik, const int jk, const vector<blitz::Array<cd,2> > &Omega_s, const vector<blitz::Array<cd,2> > & Omega_a)
{
    //cout << ik << " " << jk << "\n";
    double kx  = 2.0*M_PI*ik/L_x;
    double ky  = 2.0*M_PI*jk/L_x;
    double kx_t  = 1.0/step_x*sin(kx*step_x);
    double ky_t  = 1.0/step_x*sin(ky*step_x);
    double k2 = pow(2.0/step_x*sin(kx*step_x/2),2) +  pow(2.0/step_x*sin(ky*step_x/2),2);

    cd sum = cd(0,0);

    for(int a=0; a<8; a++) {
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

void components(colorArr& Omega, vector<blitz::Array<cd,2> > & Omega_comp)
{
    for(int a=0; a<Omega_comp.size(); a++) {
        Omega_comp.at(a).resize(Omega.shape());

        for(int i=0; i<size_x; i++)
            for(int j=0; j<size_x; j++) {
                Omega_comp.at(a)(i,j) = su3_group_element(Omega(i,j), a);
            }
    }
}

int sign(int i)
{
    if(i<0) return -1;
    return 1;
}

class ParallelStream
{
    std::ostringstream stdStream;
public:
    ParallelStream() {}
    template <class T>
    ParallelStream& operator<<(const T& inData)
    {
        stdStream << inData;
        return *this;
    }
    std::string toString() const
    {
        return stdStream.str();
    }
};



cd Dipole(const int ik, const int jk, const vector<blitz::Array<cd,2> > &V_W_c)
{
	cd sum=cd(0,0);
	for(int a=0;a<V_W_c.size(); a++)
	{
		sum += conj( V_W_c.at(a)(ik,jk) ) * V_W_c.at(a)(ik,jk);  
	}
	sum *= 0.5; 
	return sum; 
}

cd Quadrupole(const int ik1, const int jk1, const int ik2, const int jk2, const vector<blitz::Array<cd,2> > &V_W_c)
{
	cd sum = cd(0,0);
  	vector < TinyMatrix<cd,3,3>  > lambda = get_lambda();
	TinyMatrix<cd,3,3> tmp;  
	for(int a=0;a<V_W_c.size(); a++)
	for(int b=0;b<V_W_c.size(); b++)
	for(int c=0;c<V_W_c.size(); c++)
	for(int d=0;d<V_W_c.size(); d++)
	{
		tmp=lambda.at(d);
		tmp=Product(lambda.at(c),tmp); 
		tmp=Product(lambda.at(b),tmp); 
		tmp=Product(lambda.at(a),tmp);

		sum += conj( V_W_c.at(a)(ik1,jk1) ) * V_W_c.at(b)(ik1,jk1) *
			   conj( V_W_c.at(c)(ik2,jk2) ) * V_W_c.at(d)(ik2,jk2) * trace(tmp) ;  
	}

	return sum; 
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



void outHBT(double A, double B, const colorArr& V_c)
{
	clock_t begin = clock();
    vector<blitz::Array<cd,2> > fftV_c_c(9);

    for(int a=0; a<9; a++) {
        cerr <<  "cycle " << a  << "\n";
        fftV_c_c.at(a).resize(size_x,size_x);


        blitz::Array<cd,2> tmp1(size_x,size_x);
        blitz::Array<cd,2> tmp2(size_x,size_x);


        for(int i=0; i<size_x; i=i+1) {
            for(int j=0; j<size_x; j=j+1) {
                double x = int_to_x(i);
                double y = int_to_x(j);
        	double x_c = int_to_x(size_x/2); 
        	double y_c = int_to_x(size_x/2); 
			double R = 5.0; 
                tmp1(i,j) = su3_group_element(V_c(i,j), a) 
					* exp(- pow((x-x_c)/R/A,2) - pow((y-y_c)/R/B,2)) ;
            }
        }

		
        cerr << tmp1(0,0) << "\n" << flush;

		FFTW(tmp1,tmp2);

        fftV_c_c.at(a)=tmp2*pow(step_x,2);

        cerr << fftV_c_c.at(a)(0,0) << "\n" << flush;

    }

    clock_t end = clock();

    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    cerr << "FFT components done\n" << "time " << elapsed_secs   << endl << flush;

	double dK = step_x*4;
	double DK = step_x*4;
	double Kmax = 5+DK;

	for(double K=DK;K<Kmax;K+=DK)
	{
	double v0_2 = 0.0;
	double v2_2 = 0.0; 
	double v3_2 = 0.0; 

	double vQ2_2 = 0.0; 
	double vQ3_2 = 0.0; 

	double vD2_2 = 0.0; 
	double vD3_2 = 0.0; 

	int N=0;

	//cerr<< "K=" <<K << endl;
    for(int i1=0; i1<size_x; i1=i1+1) {
      double kx1_t  = int_to_k(i1);
      for(int j1=0; j1<size_x; j1=j1+1) {
	
		  		
                double ky1_t  = int_to_k(j1);
                double k12 = int_to_k2(i1,j1); 
	
				if( ((sqrt(k12)-(K-0.5*dK)) * (sqrt(k12)-(K+0.5*dK)) <0.0 
	  ) 
		//	&& ((j1<size_x/4) || ((j1-3*size_x/4)*(j1-size_x/2)<0)) 
		//	&& ((i1<size_x/4) || ((i1-3*size_x/4)*(i1-size_x/2)<0)) 
			) 
	{
    for(int i2=0; i2<size_x; i2=i2+1) {
      double kx2_t  = int_to_k(i2);
      for(int j2=0; j2<size_x; j2=j2+1) {
                
                double ky2_t  = int_to_k(j2);
		  	double k22 = int_to_k2(i2,j2) ;
	
				if(( (sqrt(k22)-(K-0.5*dK)) * (sqrt(k22)-(K+0.5*dK)) <0.0
							  ) 
		//	&& ((j2<size_x/4) || ((j2-3*size_x/4)*(j2-size_x/2)<0)) 
		//	&& ((i2<size_x/4) || ((i2-3*size_x/4)*(i2-size_x/2)<0)) 
						) 
		{
			double phi_1 = atan2(ky1_t,kx1_t); 	
			double phi_2 = atan2(ky2_t,kx2_t); 

			double scalar_product = kx1_t*kx2_t+ky1_t*ky2_t; 
			double cross_product = kx1_t*ky2_t-ky1_t*kx2_t; 
			double delta_phi = atan2(cross_product,scalar_product); 

			//if(delta_phi<0) delta_phi=delta_phi+2*M_PI; 
			


			double D = dK*sqrt(k12*k22)*real(Dipole(i1,j1,fftV_c_c)*Dipole(i2,j2,fftV_c_c));  
			double Q = dK*sqrt(k12*k22)*real(Quadrupole(i1,j1,i2,j2,fftV_c_c)); 

			v0_2 += (D+Q);  
			v2_2 += (D+Q)*cos(2*delta_phi);  
			v3_2 += (D+Q)*cos(3*delta_phi); 

			vQ2_2 += Q*cos(2*delta_phi);  
			vQ3_2 += Q*cos(3*delta_phi);  

			vD2_2 += D*cos(2*delta_phi);  
			vD3_2 += D*cos(3*delta_phi);  

			N++;
			/*int phi = delta_phi * 100; 
			delta_phi = phi / 100.0;

			cout << K <<" "; 
			cout << delta_phi <<" "; 
			//cout << phi_1-phi_2 <<" "; 
			cout << D <<" "; 
			//cout << real(Dipole(i1,j1,fftV_c_c)) <<" "; 
			cout << Q << " "; 
			cout << phi_1 << " "; 
			cout << phi_2 << " "; 
			cout << kx1_t << " "; 
			cout << ky1_t << " "; 
			cout << kx2_t << " "; 
			cout << ky2_t << " "; 
			cout << endl; 
			*/
		}				

	  }}
	
	}
	
	}}
	
	if(v0_2>0) cout 
		<< A << " " << B << " " 
		<< K 
		<< " " << (v2_2/v0_2) 
		//<< " " << (v3_2/v0_2)
		<< " " << (vD2_2/v0_2) 
		//<< " " << (vD3_2/v0_2)
		<< " " << (vQ2_2/v0_2) 
		<< " " << (v0_2/N) 
		//<< " " << (vQ3_2/v0_2)
		<< endl; 
	}


}



int main(void)
{
    clock_t begin = clock();
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    cin >> eventID;

    //Target:
    begin = clock();
    colorArr V_c(size_x,size_x);
    IC_MV(V_c);
    end = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cerr << "target done\n" << "time " << elapsed_secs << flush  << endl;

	outHBT(1.0, 1.0, V_c);
    cerr << " one " << endl;
	outHBT(2.0, 0.5, V_c); 
    cerr << " two " << endl;
	outHBT(0.5, 2.0, V_c); 
    cerr << " three " << endl;
	outHBT(3.0, 3.0, V_c); 
    cerr << " four " << endl;

}
