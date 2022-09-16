
#define sqr(x) ((x)*(x));
const double G = 6.67408 * pow(10,-11) // m^3 kg^-1 s^-2

int
rhs (double , double M[], const double y[], double lag[], void *params_ptr)
{
  int N = sizeof(M)/sizeof(M[0]);
  //double lag [4*N] = {0};

  for (int a = 0; a < N; a++)
    {
      double P = 0., Q = 0.;

      for (int b = 0; b < N; b++)
	{
	  if (a!=b)
	    {
	      double d_sqr = sqr(y[4*a]-y[4*b])+sqr(y[4*a+2]-y[4*b+2]); 
	      P += pow((y[4*a]-y[4*b])/d_sqr,3./2.);
	      Q += pow((y[4*a]-y[4*b])/d_sqr,3./2.);
	    }
	}
      lag[4*a] = (y[4*a+1]);
      lag[4*a+1] = -G*P;
      lag[4*a+2] = (y[4*a+3]);
      lag[4*a+3] = -G*Q;
    }

  return GSL_SUCCESS;
}
