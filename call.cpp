# include <iostream>
# include <cstdlib>
# include <iomanip>
# include <cmath>

using namespace std;

# include "heat1d.hpp"


int main();
void heat1d_main();
double *init_condition(int x_num);
void bound(int x_num, double h[]);
double *right(int x_num, double x[], double x_cen);



int main()

{

timestamp();

cout << "\n";
cout << " This calculates the heat transfer in 1D "<<"\n";
cout << " Written by Bailey C. Hsu with modification to the version "<<"\n";
cout << " created by John Burkadt. ";

heat1d_main();

cout << "\n";
cout << " Program ended"<<"\n";
timestamp();
return 0;

}

// heat1D main function, writes to the file, no return values
void heat1d_main()
{
  double diff_const;
  double t_max;
  double t_min;
  int t_num;
  double *t;
  double x_max;
  double x_min;
  int x_num;
  double *x;
  double dt;
  double coeff;
  double *h;
  double *hnew;
  double *hstep;
  int i;
  int j;
  double x_cen;
  

  cout << "\n";
  cout << "	 Heat flow numerical simulations:\n";
  cout << "  Compute an approximate solution to the time-dependent\n";
  cout << "  one dimensional heat equation:\n";
  cout << "\n";
  cout << "    dH/dt - K * d2H/dx2 = f(x,t)\n";
  cout << "\n";


  diff_const = 0.003;
  x_num = 11;
  x_min = 0.0;
  x_max = 2.0;
  x_cen = 0.50*(x_max-x_min);

  x = new double[x_num];
  x = linspace(x_num, x_min, x_max);
  

  t_num = 201;
  t_min = 0.0;
  t_max = 640.0;
   
  t = new double[t_num];
  t = linspace(t_num, t_min, t_max);
  dt = (t_max-t_min)/double(t_num-1);

  coeff = coefficient(x_num, x_min, x_max, t_num, t_min, t_max,
  	                  diff_const);

  // Setting initial conditions and 

  h = init_condition(x_num);
  bound(x_num, h);



  hnew = new double[x_num*t_num];

  // Setting the initial heat portfolio

 
  for (i = 0; i< x_num; i++)
  {
   hnew[i]=h[i];
  }

  
  for (j = 1; j< t_num; j++)
  {
   hstep = heat1d(t[j-1], dt, x_num, x, coeff, x_cen, right, bound, h);


    for (i = 0; i< x_num; i++)
    {
    hnew[i+j* x_num] = hstep[i];
    h[i] = hstep[i];
    }
    

    delete [] hstep;
	}

 

  write_matrix ( "temperature1.txt", t_num, x_num, hnew);
  write_vec ( "time1.txt", t_num, t );
  write_vec ( "position1.txt", x_num, x );

  delete [] h;
  delete [] hnew;
  delete [] t;
  delete [] x;

  return;
}



double *init_condition(int x_num)
{

		double *h;
		int j;

		h = new double[x_num];

		for (j = 0; j< x_num; j++)
		{
		   h[j] = 50.0;
		}
		return h;

}

void bound(int x_num, double h[])
{
		h[0] = 90.0;
		h[x_num-1] = 70.0;

		return;
}



double *right(int x_num, double x[], double x_cen)
{
		double *value;
		int j;

		value = new double[x_num];

		for (j=0; j< x_num; j++)
		{
		   value[j]= 0.2*exp(-pow(x[j]-x_cen, 2)*0.05);
		}
		 
		return value;
}
