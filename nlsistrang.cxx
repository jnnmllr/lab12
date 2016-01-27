#include <iostream>
#include <fstream>
#include <complex>
#include <sstream>
#include <cmath>
//-----------------------------------
using namespace std;
//-----------------------------------
typedef complex<double> cmplx;
//-----------------------------------
void init( cmplx* const psi0, const double eta, const double sigma, const double dx,
          const int Nx);

void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin);

void linstep(cmplx* const f1, cmplx* const f0,
          const double dt, const double dx, const int N);

void nonlinstep(cmplx* const psi1, cmplx* const psi0,
          const double dt, const int N);


//-----------------------------------
int main(){

	const int Nx = 4000;
	const double L = 800;
	const double xmin = 0;
	const double Tend = 150;
	const double dx = L / (Nx - 1);
	const double dt = dx  / 10;
	const int Na = 10;
	int Nk = int(Tend / Na / dt + 0.5);

	const double eta = 0.2;

	stringstream strm;

	cmplx* psi0 = new cmplx[Nx];
	cmplx* psi1 = new cmplx[Nx];
	cmplx* h;

	init(psi0, eta, dx, dt,Nx);
	


	writeToFile(psi0,"psistrang_0", dx,Nx,xmin);


	for (int i = 1; i <= Na; i++) {

		for (int j = 1; j <= Nk-1; j++) {	
		  linstep(psi1, psi0, (dt/2.0), dx, Nx);
		  nonlinstep(psi1, psi0, dt, Nx);
		  linstep(psi1, psi0, (dt/2.0), dx, Nx);
		        h = psi0;
			psi0 = psi1;
			psi1 = h;
		}
		strm.str("");
		strm << "psistrang_" << i;
		writeToFile(psi0,strm.str(), dx,Nx,xmin);
	}

	return 0;
}
//-----------------------------------
void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin)
{
	ofstream out(s.c_str());
	for(int i=0; i<Nx; i++){
		double x = xmin + i * dx;
		out << x << "\t" << norm(v[i]) << "\t" << v[i].real() << "\t" << v[i].imag() << endl;
	}
	out.close();
}
//-----------------------------------

void init( cmplx* const psi0, const double eta,  const double dx, const double dt,
          const int Nx)
{
	const double x0 = dx*Nx * 0.5;
	const double f = sqrt(2) * eta;
	for(int i=0;i<Nx; i++){
		double x = i*dx - x0;
		psi0[i] = 2*f/cosh(eta * x);
	}
}

void linstep(cmplx* const f1, cmplx* const f0,
          const double dt, const double dx, const int N)
{

  cmplx* d=new cmplx[N];
 cmplx j = cmplx(0.0,1.0);
 
  for(int i=0;i<N;i++) d[i] = 1.0 - 2.0*j*dt/(dx*dx);
			  const cmplx u = (1.0)*j*dt/(dx*dx);
			  const cmplx l = (1.0)*j*dt/(dx*dx);
  
      for(int i=1;i<N;i++){
      d[i] -= u*l/d[i-1];
      f0[i] -= f0[i-1]*l/d[i-1];
      }
      
  
      f1[N-1]= f0[N-1]/d[N-1];
      for(int i=N-2;i>=0;i--){
      f1[i]=(f0[i]-u*f1[i+1])/d[i];
      }


  delete[] d;
}

void nonlinstep(cmplx* const psi1, cmplx* const psi0,
          const double dt, const int N)
{
  cmplx i = cmplx(0.0,1.0);
  for(int j=0; j<N; j++){
    psi0[j] = psi1[j]*exp(-i*psi1[j]*conj(psi1[j])*dt);
  }
}
