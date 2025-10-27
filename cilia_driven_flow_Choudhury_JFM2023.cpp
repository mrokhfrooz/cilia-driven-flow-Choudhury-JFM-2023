
#include <cmath> 
#include <iostream>
#include <random>
#include <fstream>
#include <vector>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <unistd.h> 
#pragma warning(disable : 4996) 
using namespace std;



void Choudhury_JFM_newotonian(double b, double anji_x, double anji_y, double anji_t, double& anji_velocity_u, double& anji_velocity_v) {
    
	b= 6e-6;
	double pi = 3.14159265358979323846;
	double anji_freqency_Hz= 10; // unit= 1/sec taken from Choudhury Fig 1
	double anji_angular_frequency= 2*pi*anji_freqency_Hz; // \omega* in Anji paper
	double anji_velocity_scale=b*anji_angular_frequency;
  


	anji_x = anji_x/b; // made it dimensionsless
	anji_y = anji_y/b; // made it dimensionsless
	anji_t = anji_t* anji_angular_frequency ; //pi = especial case,  ; // make it dimensionsless
	
	anji_t= pi ;


	
		
    double anji_zeta= 1; // binary parameter, can be 0 or 1
	double anji_phi = 0e-6 /b  ; // the slip length, phi-star can be 10 micrometer for example
	double anji_a= (1.6e-6) / b;  // a* = 1.6 um, value taken from Anji paper
	double anji_box_width= 1*b; 
	
	
	double anji_k= ( (2*pi)*b ) / (anji_box_width); 
	double anji_k_prime= 2*anji_k;
	
	double anji_S_cont=  sinh(2*anji_k)      - (2*anji_k); // anji_S_cont is a big number
	double anji_S_prime= sinh(2*anji_k_prime)- (2*anji_k_prime);
	
	double anji_C_cont=  cosh(2*anji_k)       -1;
	double anji_C_prime= cosh(2*anji_k_prime) -1;
	
	
	double anji_coeff_A1= -anji_a/( anji_S_cont+ (2*anji_k*anji_phi*anji_C_cont)  );
	double anji_coeff_B1=  (anji_a*exp(-anji_k)*sinh(anji_k) )/( anji_k* (anji_S_cont+ (2*anji_k*anji_phi*anji_C_cont) )  );
	double anji_coeff_C1= - anji_coeff_A1;
	double anji_coeff_D1= anji_coeff_B1*exp(2*anji_k);	
	double anji_coeff_A2= 0.5*anji_a*anji_a*anji_k; 	
		
	double anji_coeff_A3= (anji_a*anji_a*anji_k_prime)/(4* (anji_S_prime+ (2*anji_k_prime*anji_phi*anji_C_prime) ) ); 
	double anji_coeff_B3= (-1*anji_a*anji_a* (1-exp(-2*anji_k_prime) ) ) / (8* (anji_S_prime+ (2*anji_k_prime*anji_phi*anji_C_prime) ) );
	double anji_coeff_C3= -anji_coeff_A3;
	double anji_coeff_D3= anji_coeff_B3*exp (2*anji_k_prime);
	




	
	
	for (int i=0;i<21;i++) { // discritization on x-axis
		for (int j=0;j<21;j++) { // discretization on y-axis 
			
			 double anji_x= i*0.6e-6 / (2*b); 
			 double anji_y= j*0.6e-6 / (2*b);

	double anji_analytical_veloc_v1 = ( ( ((anji_coeff_A1+anji_coeff_B1*anji_k*anji_y)*exp(anji_k*anji_y)) + ( (anji_coeff_C1+anji_coeff_D1*anji_k*anji_y)*exp(-anji_k*anji_y) )) *anji_k* sin(anji_k*anji_x+anji_t) ) ;
	double anji_analytical_veloc_v2 = ( anji_zeta*(((anji_coeff_A3+2*anji_coeff_B3*anji_k*anji_y)*exp(2*anji_k*anji_y)) + ((anji_coeff_C3+2*anji_coeff_D3*anji_k*anji_y)* exp(-2*anji_k*anji_y)))*(2*anji_k)*sin(2*(anji_k*anji_x+anji_t))  );
	double anji_analytical_veloc_v  = anji_analytical_veloc_v1+ anji_analytical_veloc_v2;
	
	double anji_analytical_veloc_u1 = (( (anji_k*exp(anji_k*anji_y)*(anji_coeff_A1+anji_coeff_B1+ (anji_coeff_B1*anji_k*anji_y) ) ) + ( anji_k*exp(-anji_k*anji_y)* (anji_coeff_D1- anji_coeff_C1-(anji_coeff_D1*anji_k*anji_y) ) ) )* (cos(anji_k*anji_x+anji_t)) );
	double anji_analytical_veloc_u2 =anji_zeta*anji_coeff_A2;
	double anji_analytical_veloc_u3 =anji_zeta* ( ( (anji_coeff_A3+ anji_coeff_B3+ 2*anji_coeff_B3*anji_k*anji_y)*2*anji_k*exp(2*anji_k*anji_y) ) + ( (anji_coeff_D3-anji_coeff_C3- (2*anji_coeff_D3*anji_k*anji_y) )*2*anji_k*exp(-2*anji_k*anji_y) ) ) * cos( 2*(anji_k*anji_x+anji_t) );
	double anji_analytical_veloc_u =anji_analytical_veloc_u1+anji_analytical_veloc_u2+anji_analytical_veloc_u3;
	
	cout << anji_x << ";" << anji_y << ";" << anji_analytical_veloc_u << ";" << anji_analytical_veloc_v << endl;
	

  
		}
	}

}







double dotProduct(vector<double> a, vector<double> b) {  // input must be 2 3-elements vectors
	double result_dot = 0;
	result_dot = (a[0] * b[0]) + (a[1] * b[1]) + (a[2] * b[2]);
	return result_dot;
}

vector<double> crossProduct(vector<double> a, vector<double> b) {
	vector<double> result(3);

	result[0] = a[1] * b[2] - a[2] * b[1];
	result[1] = a[2] * b[0] - a[0] * b[2];
	result[2] = a[0] * b[1] - a[1] * b[0];

	return result;
}











int main() {
	
	const int N = 3; 
	const int step = 10000;  // number of steps
	const int Saveeverystep = 1000;
	int verlet_restart = 25;  

	double shear_flow_rate_x_direction= 90*1e-6 ; // 
	int oscillation_direction=1;  
	int oscillation_reverse_step= 50*1000; 
	double b = 6000e-9; 
	
	double mucus_velocity_x_comp=0;
	double mucus_velocity_z_comp=0;
	double mucin_location_x=0;
	double mucin_location_z=0;
	double zaman=0;

	Choudhury_JFM_newotonian(b, mucin_location_x, mucin_location_z , zaman,  mucus_velocity_x_comp,  mucus_velocity_z_comp);
		
	

	return 0;

}


