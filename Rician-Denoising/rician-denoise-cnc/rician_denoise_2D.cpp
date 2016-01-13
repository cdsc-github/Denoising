/*========================================================================
 *
 * RICIANDENOIESMX.C  Total variation minimization for Rician denoising
 *
 * u = riciandenoisemx(f,sigma,lambda,Tol) performs denoising on image f
 * with Rician noise with parameter sigma.  The denoised image u is found
 * as the minimizer of 
 *
 *         /                      / [ u^2 + f^2            u f    ]
 *    min  | |grad u| dx + lambda | [ --------- - log I0( ----- ) ] dx.
 *     u   /                      / [ 2 sigma^2          sigma^2  ]
 *
 * Parameter lambda >= 0 determines the strength of the denoising: smaller
 * lambda implies stronger denoising.  Tol specifies the stopping 
 * tolerance, the method stops when ||u^Iter - u^Iter-1||_inf < Tol.
 *
 *======================================================================*/
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <utility>

typedef std::pair<int, int> pair;

#include "rician_denoise_2D.h"

/* Method Parameters */
#define DT          5.0
#define EPSILON     1.0E-20
#define MAXITER     500
#define SIGMA       0.05
#define LAMBDA      0.065
#define TOL         2e-3
#define SIGMA2      0.0025
#define GAMMA       LAMBDA/SIGMA2

/* Macro functions */
#define SQR(x) ((x)*(x))

int compute_g::execute(const pair & t, rician_denoise_2D_context & c) const
{
	const int i = t.first;
	const int j = t.second;
	
	double value = 1.0 / sqrt( EPSILON + 
	               SQR( c.u.get(pair(i, j)) - c.u.get(pair(i, j-1)) ) + 
						SQR( c.u.get(pair(i, j)) - c.u.get(pair(i, j+1)) ) + 
						SQR( c.u.get(pair(i, j)) - c.u.get(pair(i-1, j)) ) + 
						SQR( c.u.get(pair(i, j)) - c.u.get(pair(i+1, j)) ) );	

	c.g.put(pair(i, j), value);
	c.ug_tag.put(pair(i, j));

	return CnC::CNC_Success;
}

int compute_r::execute(const pair & t, rician_denoise_2D_context & c) const
{
	const int i = t.first;
	const int j = t.second;

	double value = c.u.get(pair(i, j)) * c.f.get(pair(i, j)) / SIGMA2;

	value = ( value * (2.38944 + value * (0.950037 + value)) )
	/ ( 4.65314 + value * (2.57541 + value * (1.48937 + value)) );

	c.r.put(pair(i, j), value);

	return CnC::CNC_Success;
}

int compute_ug::execute(const pair & t, rician_denoise_2D_context & c) const
{
	const int i = t.first;
	const int j = t.second;

	double value = c.u.get(pair(i, j)) * c.g.get(pair(i, j));

	c.ug.put(pair(i, j), value);
	c.unext_tag.put(pair(i, j));

	return CnC::CNC_Success;
}

int compute_unext::execute(const pair & t, rician_denoise_2D_context & c) const
{
	const int i = t.first;
	const int j = t.second;

	double value = ( c.u.get(pair(i, j)) + 
	                 DT * ( c.ug.get(pair(i, j+1)) + c.ug.get(pair(i, j-1)) + 
						  c.ug.get(pair(i+1, j)) + c.ug.get(pair(i-1, j)) + 
						  GAMMA * c.f.get(pair(i, j)) * c.r.get(pair(i, j)) ) ) /
	                 ( 1.0 + DT * ( c.g.get(pair(i, j+1)) + c.g.get(pair(i, j-1)) + 
						    c.g.get(pair(i+1, j)) + c.g.get(pair(i-1, j)) + GAMMA ) );
	
	c.unext.put(pair(i, j), value);
	c.converged_tag.put(pair(i, j));

	return CnC::CNC_Success;
}

int compute_convergence::execute(const pair & t, rician_denoise_2D_context & c) const
{
	const int i = t.first;
	const int j = t.second;

	double diff = fabs( c.u.get(pair(i, j)) - c.unext.get(pair(i, j)) );

	c.difference.put(pair(i, j), diff);
	//printf("unext = %lf, diff = %lf, i = %d, j = %d\n", c.unext.get(pair(i, j)), diff, i, j);

	return CnC::CNC_Success;
}

void riciandenoise2D_Jacobi(double **u, double **f, int M, int N)
{
    int Converged;
    int m, n;
    int Iter;    
	 int i, j;

   /* Initialize u = f */
	for( i = 0; i < M; i++ )
		for( j = 0; j < N; j++ )
			u[i][j] = f[i][j];


	// Main gradient descent loop
	for(Iter = 1; Iter <= MAXITER; Iter++) {
		// Create an instance of the context class wich defines the graph
		rician_denoise_2D_context c;	
		
		// Give inputs from the environment to the graph
		Converged = 1;
		for( i = 0; i < M; i++ ) 
			for( j = 0; j < N; j++ ) {
				c.u.put(pair(i, j), u[i][j]);
				c.f.put(pair(i, j), f[i][j]); 
			}

		for( i = 1; i < M-1; i++ )
			for( j = 1; j < N-1; j++ )
				c.Iter_start.put(pair(i, j));

		for( i = 2; i < M-2; i++ )
			for( j = 2; j < N-2; j++ )
				c.unext_tag.put(pair(i, j));
		
		// Wait for all steps to finish
		c.wait();

		double temp;

		for( i = 2; i < M-2; i++ )
			for( j = 2; j < N-2; j++ ) {
				temp = c.difference.get(pair(i, j));
				//printf("diff(%d,%d):%lf\n", i, j, temp);
				
				if( temp > TOL ) {
					Converged = 0;
				}
				
				u[i][j] = c.unext.get(pair(i, j));
			}

		if(Converged)
			break;
		
		printf("Iter:%d\n", Iter);


	}

    /* Done, show exiting message */
    if(Converged)
        printf("Converged in %d iterations with tolerance %g.\n", Iter, TOL);
    else
        printf("Maximum iterations exceeded (MaxIter=%d).\n", MAXITER);

    return;
}

/* readImgArray2D: Reaa an input array from a file */
void readImgArray2D(char* in_file, double **img_array) {
	FILE* inputfile;
	int i, j, k;
	char tmp_str[100];
	int c;
	char tmp_c;


	inputfile = fopen(in_file, "r"); 

	if( inputfile == NULL ) {
		printf("Cannot open the input file!\n");
		exit(1);
	}	

	i = j = k = 0;

	while( !feof(inputfile) ) {
		c = fgetc(inputfile);

		if( c == '\n' ) {
			i++;
			j = 0;
		}
		else if(c == ' ') {
			tmp_str[k] = '\0';
			img_array[i][j] = atof(tmp_str);
			k = 0;
			j++;
		}
		else if(c == EOF) 
				break;
		else {
			tmp_c = (char) c;
			tmp_str[k]= tmp_c;
			k++;
		}
	}
	
	fclose(inputfile);
	
}

int main(int argc, char* argv[]) {
	if(argc < 4) {
		printf("riciandenoisemx M N inputfile\n");
		exit(0);
	}

	int M = atoi(argv[1]);
	int N = atoi(argv[2]);
	int i, j;
	double **img_array;
	double **processed_img_array;
	FILE* outputfile;	
	
	// Initialize image array
	img_array = (double**) calloc(M, sizeof(double*) );
	for( i = 0; i < M; i++ )
		img_array[i] = (double*) calloc(N, sizeof(double) );
	
	// Initialize the image array after processing
	processed_img_array = (double**) calloc(M, sizeof(double*) );
	for( i = 0; i < M; i++ )
		processed_img_array[i] = (double*) calloc(N, sizeof(double) );

	readImgArray2D(argv[3], img_array);
	riciandenoise2D_Jacobi(processed_img_array, img_array, M, N);

	outputfile = fopen("rician_denoise", "w");	

	for( i = 0; i < M; i++ ) {
		for( j = 0; j < N; j++ ) {
 			fprintf(outputfile, "%lf ", processed_img_array[i][j]);
		}
		fprintf(outputfile, "\n");
	}

	fclose(outputfile);
	
	return 1;
}

