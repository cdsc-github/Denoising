/* 3D Gaussian convolution using the method of Alvarez and Mazorra
 * Pascal Getreuer 2009
 */
#include<math.h>

#define GAUSSIAN_NUMSTEPS  3

/* GaussianBlur: Implements 3D Gaussian convolution with recursive (IIR)
 filters.  The parameter GAUSSIAN_NUMSTEPS determines the quality of the
 convolution.  More steps is a better approximation of the true Gaussian.*/
void GaussianBlur(double *u, const int Size[3], double Ksigma)
{
    const int PlaneStep = Size[0]*Size[1];
    double *uPtr, *uCopy, *uEnd;
    double lambda = (Ksigma*Ksigma)/(2.0*GAUSSIAN_NUMSTEPS);
    double nu = (1.0 + 2.0*lambda - sqrt(1.0 + 4.0*lambda))/(2.0*lambda);
    double BoundaryScale = 1.0/(1.0 - nu);
    double PostScale = pow(nu/lambda,3*GAUSSIAN_NUMSTEPS);
    int Step = GAUSSIAN_NUMSTEPS;
    int n[3];
    
    
    uEnd = u + PlaneStep*Size[2];
    
    do
    {
        for(n[2] = 0, uPtr = u; n[2] < Size[2]; ++n[2], uPtr+=PlaneStep)
        {
            uCopy = uPtr; 
            
            for(n[1] = 0; n[1] < Size[1]; ++n[1], uPtr+=Size[0])
            {
                /* Filter downwards */                
                uPtr[0] *= BoundaryScale;
                ++uPtr;
                
                for(n[0] = 1; n[0] < Size[0]; ++n[0], ++uPtr)
                {
                    uPtr[0] += nu*uPtr[-1];
                }
                
                /* Filter upwards */
                --uPtr;
                uPtr[0] *= BoundaryScale;
                --uPtr;
                
                for(n[0] = Size[0]-2; n[0] >= 0; --n[0], --uPtr)
                {
                    uPtr[0] += nu*uPtr[1];
                }
                
                ++uPtr;
            }
            
            uPtr = uCopy;
            
            /* Filter right */
            for(n[0] = 0; n[0] < Size[0]; ++n[0], ++uPtr)
            {
                uPtr[0] *= BoundaryScale;
            }
            
            for(n[1] = 1; n[1] < Size[1]; ++n[1])
            {
                for(n[0] = 0; n[0] < Size[0]; ++n[0], ++uPtr)
                {
                    uPtr[0] += nu*uPtr[-Size[0]];
                }
            }
            
            --uPtr;
            
            /* Filter left */
            for(n[0] = Size[0]-1; n[0] >= 0; --n[0], --uPtr)
            {
                uPtr[0] *= BoundaryScale;
            }
            
            for(n[1] = Size[1]-2; n[1] >= 0; --n[1])
            {
                for(n[0] = Size[0]-1; n[0] >= 0; --n[0], --uPtr)
                {
                    uPtr[0] += nu*uPtr[Size[0]];
                }
            }
            
            ++uPtr;
        }
        
        /* Filter out */
        n[0] = PlaneStep;
        uPtr = u;
            
        do
        {
            uPtr[0] *= BoundaryScale;
            ++uPtr;
        }while(--n[0]);
        
        for(n[2] = 1; n[2] < Size[2]; ++n[2])
        {
            n[0] = PlaneStep;
            
            do
            {
                uPtr[0] += nu*uPtr[-PlaneStep];
                ++uPtr;
            }while(--n[0]);
        }
        
        /* Filter in */
        n[0] = PlaneStep;
        
        do
        {
            --uPtr;
            uPtr[0] *= BoundaryScale;            
        }while(--n[0]);
        
        
        for(n[2] = Size[2]-2; n[2] >= 0; --n[2])
        {
            n[0] = PlaneStep;
            
            do
            {
                --uPtr;
                uPtr[0] += nu*uPtr[PlaneStep];
            }while(--n[0]);
        }
    }while(--Step);
    
    do
    {
        u[0] *= PostScale;
    }while(++u < uEnd);
}
