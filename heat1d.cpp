# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cmath>

using namespace std;

# include "heat1d.hpp"

//Main function to calculate the evolution of heat

double *heat1d(double t, double dt, int x_num, double x[], double coeff, double x_cen,
     double *right(int x_num, double x[], double x_cen), void bc(int x_num, double h[]), double h[])
     	{
			double *f;
			double *h2;
			int j;

			// first we calculate the f, the source distribution

			f = right(x_num, x, x_cen);

		
			// Initialize the H2 matrix, the evolved heat vectors
			h2 = new double[x_num];
			h2[0] = 0.0;

			for (j=1; j< x_num-1; j++)
			 	{
			     h2[j] = h[j] + dt* f[j]+ coeff*(h[j-1]-2.0* h[j]+ h[j+1]);
				}	

			h2[x_num-1] = 0.0;

			bc(x_num,h2);

			delete [] f;
			return h2;

}


//Modified Coefficient for the Heat Equation

double coefficient(int x_num, double x_min, double x_max, int t_num, double t_min, double t_max,
	               double diff_const){
     
     double coeff;
     double dx;
     double dt;

     dx = (x_max-x_min)/ (double) (x_num-1);
     dt = (t_max-t_min)/ (double) (t_num-1);

     coeff = diff_const* dt/dx/dx;

     cout << "\n";
     cout << " the modified coefficient is coeff = " << coeff <<"\n";

     return coeff;
}


//Write an array to a file

void write_vec(string filename, int n, double x[]){
		int j;
		ofstream output;

		output.open(filename.c_str());

		if (!output)

		{
		  cerr << "\n";
		  cerr << "Could not open file!";
		  exit(1);

		} else {
		  for (j=0; j<n; j++){
		  output << " " << setw(24) << setprecision(16) << x[j] << "\n" ;
		  }

		  output.close();
		  return;

		}
	}


void write_matrix(string filename, int n, int m, double table[]) {
		int i;
		int j;

		ofstream output;

		output.open(filename.c_str());

		if (!output){
		  cerr << "\n";
		  cerr << "Could not open file!";
		  exit(1);
		} 

		else 
		  {
		  for (j=0; j<n; j++)
		  {
		  	for (i=0; i<m; i++)
		  	{
		  	output << ""<<setw(24)<<setprecision(16)<<table[i+j*m];
		    }	
            output << "\n" ;
		   }

		  output.close();
		  return;
		}

	}


//creates linearly spaced vectors, xmin+ (xmax-xmin)/(n-1)*i


double *linspace(int n, double x_min, double x_max){

		double *x;
		int     i;

		x = new double[n];

		if (n ==1)
		{
		  x[0]=(x_min + x_max)/2.0;

		}
		else
		{
		  for (i=0; i < n; i++)
		  {
		  	x[i]=((double)(n-1-i)*x_min+(double)(i)*x_max)/(double)(n-1);
		  }
		}

		return x;
		 
}

void timestamp(){
		# define TIME_SIZE 40

		static char time_buffer[TIME_SIZE];
		const struct std::tm *tm_ptr;

		size_t len;
		std::time_t now;

		now = std::time ( NULL );
		tm_ptr = std::localtime ( &now );

		len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

		std::cout << time_buffer << "\n";

		return;
		# undef TIME_SIZE
}







