/*
  --- simple program that prints out the Loop eigenstructures that were
  precomputed.

  Author : Jos Stam (jstam@aw.sgi.com), Alias|wavefront
*/

#include <stdio.h>
#include <stdlib.h>

#define IX(i,j,n) ((i)+(n)*(j))

typedef
struct
{
  double * val;
  double * vecI;
  double ** Phi;
} EVALSTRUCT;

/* --- this routine reads the eigen-data from the file structure and returns a
   pointer to a data structure containing it */

static EVALSTRUCT ** read_eval ( int * pNmax )
{
  EVALSTRUCT ** ev;
  FILE * f;
  int Nmax, i, N, K;

  if ( !(f=fopen("lpdata50NT.dat","rb")) ) return ( NULL );
  //printf("%d", sizeof(int));
  //printf("........\n");
  
  fread ( &Nmax, sizeof(int), 1, f );
  //Nmax = 6;
  ev = (EVALSTRUCT **) malloc ( (Nmax-2)*sizeof(EVALSTRUCT *) );
  
  for ( i=0 ; i<Nmax-2 ; i++ )
    {
      N = i+3;
      K = N+6;

      ev[i] = (EVALSTRUCT *) malloc ( sizeof(EVALSTRUCT) );
      ev[i]->val = (double *) malloc ( K*sizeof(double) );
      ev[i]->vecI = (double *) malloc ( K*K*sizeof(double) );
      ev[i]->Phi = (double **) malloc ( 3*sizeof(double) );
      ev[i]->Phi[0] = (double *) malloc ( K*12*sizeof(double) );
      ev[i]->Phi[1] = (double *) malloc ( K*12*sizeof(double) );
      ev[i]->Phi[2] = (double *) malloc ( K*12*sizeof(double) );

      fread ( ev[i]->val, sizeof(double), K, f );
      fread ( ev[i]->vecI, sizeof(double), K*K, f );
      fread ( ev[i]->Phi[0], sizeof(double), K*12, f );
      fread ( ev[i]->Phi[1], sizeof(double), K*12, f );
      fread ( ev[i]->Phi[2], sizeof(double), K*12, f );
    }

  fclose ( f );

  *pNmax = Nmax;

  return ( ev );
}

/* // the following two functions were copy-pasted from paper. iV is not a KxK matrix in the code */
/* void ProjectPoints ( point *Cp, point *C, int N ) { */

/*   for ( i=0 ; i<N+6 ; i++ ) { */
/*      Cp[i]=(0,0,0); */
/*      for ( j=0 ; j<N+6 ; j++ ) { */
/*        Cp[i] += eigen[N].iV[i][j] * C[j]; */
/*      } */
/*   } */
/* } */

/* void EvalSurf(point P, double v, double w, point* Cp, int N){ */
/*   // determine in which domaon Omega_k^n the parameter lies */

/*   n = floor(1-log2(v+w)); */
/*   pow2 = pow(2,n-1); */
/*   v *= pow2; w *= pow2; */
/*   if (v> 0.5){ */
/*     k=0; v=2*v-1; w=2*w; */
/*   } */
/*   else if (w>0.5){ */
/*     k=2; v=2*v; w=2*w-1; */
/*   } */
/*   else{ */
/*     k=1; v=1-2*v; w=1-2*w; */
/*   } */
/*   // Evaluate the surface */
/*   P  = {0,0,0}; */
/*   for (i=0; i<N+6; i++){ */
/*     P += pow(eigen[N].L[i],n-1) * EvalBasis(eigen[N].Phi[i][k],v,w) * Cp[i]; */
/*   } */
/* } */
/* --- just print out the eigen data */

static void print_eval ( EVALSTRUCT ** ev, int Nmax )
{
  int i, j, k, l, N, K;

  printf ( "Nmax = %d\n\n", Nmax );

  for ( i=0 ; i<Nmax-2 ; i++ )
    {
      N = i+3;
      K = N+6;

      printf ( "N=%d\n\n", N );

      printf ( "Eigenvalues:\n\n" );

      for ( j=0 ; j<K ; j++ ) printf ( "%6.3f ", ev[i]->val[j] );
      printf ( "\n\n" );

      printf ( "Inverse of Eigenvectors:\n\n" );

      for ( j=0 ; j<K ; j++ )
	{
	  for ( k=0 ; k<K ; k++ ) printf ( "%6.3f ", ev[i]->vecI[IX(j,k,K)] );
	  printf ( "\n" );
	}
      printf ( "\n\n" );

      printf ( "Coefficients of the Eigenbasis functions\n\n" );

      for ( k=0 ; k<3 ; k++ )
	{
	  printf ( "k=%d :\n\n", k );

	  for ( j=0 ; j<K ; j++ )
	    {
	      for ( l=0 ; l<12 ; l++ ) printf("%6.3f ",ev[i]->Phi[k][IX(j,l,K)] );
	      printf ( "\n" );
	    }
	  printf ( "\n\n" );
	}

      printf ( "\n\n" );

    }

  return;
}

int main ( int argc, char ** argv )
{
  EVALSTRUCT ** ev;
  int Nmax;

  ev = read_eval ( &Nmax );
  if ( !ev ) exit ( 1 );

  print_eval ( ev, Nmax );

  return(0);
}
