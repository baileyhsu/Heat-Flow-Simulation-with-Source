double *heat1d(double t, double dt, int x_num, double x[], double coeff, double x_cen,
     double *right(int x_num, double x[], double x_cen), void bc(int x_num, double h[]), double h[]);

double coefficient(int x_num, double x_min, double x_max, int t_num, double t_min, double t_max,
	               double diff_const);


void write_vec(string filename, int n, double x[]);

void write_matrix(string filename, int n, int m, double table[]);

double *linspace(int n, double x_min, double x_max);

void timestamp();

