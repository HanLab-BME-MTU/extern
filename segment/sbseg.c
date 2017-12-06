/*				sbseg.c  by Tom Goldstein
 *   This code performs isotropic segmentation using the "Split Bregman" algorithm.
 * An image of dimensions mxn is denoised using the following command
 * 
 *   sbseg(f,edge,mu);
 * 
 * where:
 * 
 *   - "f" is the mxn noisy image to be segmented
 *   - "mu" is the weighting parameter for the fidelity term (mu should be about 1e-4 for images with pixels on the 0-255 scale)
 *   - "edge" is a 2D array of edge detector values at each pixel (e.g. the weights for TV regularizer).  For standard segmentation
 *         simply choose "edge = ones(size(f))"
 */

#include <math.h>
#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
typedef double num;


	
/*A method for performing an iteration of the method*/
double sbseg(double** f, double** edge,  double** u, double** x,double** y, double** bx, double** by,
			double *c1, double *c2, double mu, double lambda, int nIter, int width, int height);


/*  Methods for computing averages of each region*/

void updateC(double** f, double** u, double thresh, double *c1, double *c2, int rows, int cols);
void guessC(double** f, double** u, double *c1, double *c2, int rows, int cols);
	/******************Relaxation Methods*****************/
double gsUseg(double** f, double** u, double** x, double** y, double** bx, double** by, double c1, double c2, double mu, double lambda, int rows, int cols);
void gsSpaceseg(num** u, num** edge, num** x, num** y, num** bx, num** by, double lambda, int width, int height);
	 
	/************************Bregman***********************/
void bregmanXseg(num** x,num** u, num** bx, int width, int height);
void bregmanYseg(num** y,num** u, num** by, int width, int height);

/**********************Memory************/
	
num** newMatrixseg(int rows, int cols);
void deleteMatrixseg(num ** a);
void fillArray(num **a, double val, int rows, int cols);
double** get2dArray(const mxArray *mx, int isCopy);
double copyseg(double** source, double** dest, int rows, int cols);
	
/***********************The MEX Interface to rof_iso()***************/
	
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
		/* get the size of the image*/
	int rows = mxGetN(prhs[0]);
	int cols = mxGetM(prhs[0]);
	
		/* get the fidelity and convergence parameters*/
	double mu =  (double)(mxGetScalar(prhs[2]));
	double lambda = 0.5;
	/*double tol = (double)(mxGetScalar(prhs[3]));*/
	
	/* get the image, and declare memory to hold the auxillary variables*/
	double **f = get2dArray(prhs[0],0);
	double **edge = get2dArray(prhs[1],0);
	double **u = newMatrixseg(rows,cols);
	double **x = newMatrixseg(rows-1,cols);
	double **y = newMatrixseg(rows,cols-1);
	double **bx = newMatrixseg(rows-1,cols);
	double **by = newMatrixseg(rows,cols-1);
	
	double diff;
	int count=0;
	double c1,c2;
    double** uOld;
    double *outArray;
    double change;
    
    int i,j;
     /***********Check Conditions******/  
    if (nrhs != 3) {mexErrMsgTxt("Three input arguments required.");} 
    if (nlhs > 1){mexErrMsgTxt("Too many output arguments.");}
    if (!(mxIsDouble(prhs[0]))) {mexErrMsgTxt("Input array must be of type double.");}
    if(rows!=mxGetN(prhs[1]) || cols!=mxGetM(prhs[1])){mexErrMsgTxt("Input arrays must have same dimensions.");}
	
    
	/* segment the image*/

	uOld = newMatrixseg(rows,cols);
	guessC(f,u,&c1,&c2,rows,cols);
    
    
	do{
		change = sbseg(f,edge,u,x,y,bx,by,&c1,&c2,mu,lambda,5,rows,cols);
		diff = copyseg(u,uOld,rows,cols);
        count++;
     }while( (count<500 && change>.025) || count<5);

	/* copy denoised image to output vector*/
	plhs[0] = mxCreateDoubleMatrix(cols, rows, mxREAL); /*mxReal is our data-type*/
	outArray = mxGetPr(plhs[0]);

	for(i=0;i<rows;i++){
	    for(j=0;j<cols;j++){
	        outArray[(i*cols)+j] = u[i][j];
	    }
	}
   /* outArray[0] = (double)count;*/
	
	/* Free memory */
    deleteMatrixseg(u);
    deleteMatrixseg(uOld);
	deleteMatrixseg(x);
	deleteMatrixseg(y);
	deleteMatrixseg(bx);
	deleteMatrixseg(by);
	
    return;
}

double** get2dArray(const mxArray *mx, int isCopy){
	double* oned = mxGetPr(mx);
	int rowLen = mxGetN(mx);
	int colLen = mxGetM(mx);
	double** rval = (double**) malloc(rowLen*sizeof(double*));
    int r;
	if(isCopy){
		double *copy = (double*)malloc(rowLen*colLen*sizeof(double));
		int i, sent = rowLen*colLen;
		for(i=0;i<sent;i++)
			copy[i]=oned[i];
		oned=copy;
	}

	for(r=0;r<rowLen;r++)
		rval[r] = &oned[colLen*r];
	return rval;
}

void fillArray(num **a, double val, int rows, int cols){
	int r,c;
	for(r=0;r<rows;r++)
		for(c=0;c<cols;c++)
			a[r][c] = val;
	return;
}





/*                IMPLEMENTATION BELOW THIS LINE                         */

/******************Isotropic Segmentation**************/

double sbseg(double** f, double** edge,  double** u, double** x,double** y, double** bx, double** by,
			double *c1, double *c2, double mu, double lambda, int nIter, int width, int height){
	int j,n;
    double change=0;
	for(j=0;j<nIter;j++){	
        for(n=0;n<1;n++){
			change = gsUseg(f,u,x,y,bx,by,*c1,*c2,mu,lambda, width, height);
			gsSpaceseg(u,edge,x,y,bx,by,lambda, width, height);
            bregmanXseg(x, u, bx, width, height);
			bregmanYseg(y, u, by, width, height);
        }
			
			updateC(f, u, .5, c1, c2, width, height);
		}
    return change;
}

void updateC(double** f, double** u, double thresh, double *c1, double *c2, int rows, int cols){
	int r,c;
	double sum1=0, sum2=0;
	int n1=0,n2=0;
	double* ur,*fr;
	for(r=0;r<rows;r++){
		ur = u[r];
		fr = f[r];
		for(c=0;c<cols;c++){
			if(ur[c]>thresh){
				sum1+=fr[c];
				n1++;
			}else{
				sum2+=fr[c];
				n2++;
			}
		}
	}
	
	if(n1<1)n1=1;
	if(n2<1)n2=1;
	*c1 = sum1/n1;
	*c2 = sum2/n2;
	return;
}

void guessC(double** f, double** u, double *c1, double *c2, int rows, int cols){
		int r,c;
		double av=0,sum1=0, sum2=0;
		int n1=0,n2=0;
		double* fr;
		for(r=0;r<rows;r++){
			fr = f[r];
			for(c=0;c<cols;c++){
				av+=fr[c];
			}
		}
		av/=rows*cols;
		for(r=0;r<rows;r++){
			fr = f[r];
			for(c=0;c<cols;c++){
				if(fr[c]>av){
					sum1+=fr[c];
					n1++;
					u[r][c]=1;
				}else{
					sum2+=fr[c];
					n2++;
				}
			}
		}
		if(n1<1)n1=1;
		if(n2<1)n2=1;
		*c1 = sum1/n1;
		*c2 = sum2/n2;
		return;
	}

/****Relaxation operators****/

double gsUseg(double** f, double** u, double** x, double** y, double** bx, double** by, double c1, double c2, double mu, double lambda, int rows, int cols){
		int r,c;
		double g,sum,temp1,temp2;
		  /* optimization varables*/
		double *xr,*xrm1,*bxr,*bxrm1,*yr,*byr;
		double *fr,*ur,*urm1,*urp1;
		int cm1;
        double old, diff, change=0;
		for(r=1;r<rows-1;r++){
			xr = x[r];
			xrm1=x[r-1];
			bxr=bx[r];
			bxrm1=bx[r-1];
			yr=y[r];
			byr=by[r];
			fr = f[r];
			ur = u[r];
			urm1 = u[r-1];
			urp1 = u[r+1];
			for(c=1;c<cols-1;c++){
                old = ur[c];
				cm1 = c-1;
				temp1 = c1-fr[c];
				temp2 = c2-fr[c];
				g = temp1*temp1-temp2*temp2;
				sum = xrm1[c]-xr[c]-bxrm1[c]+bxr[c]+yr[cm1]-yr[c]-byr[cm1]+byr[c]-mu/lambda*g;
				ur[c] = 0.25*(urm1[c]+urp1[c]+ur[c+1]+ur[cm1]+sum);
				if(ur[c]>1) ur[c]=1;
				else if(ur[c]<0) ur[c]=0;
                diff = ur[c]-old; diff*=diff; if(diff>change) change=diff;
			}
		}
		r=0;
			for(c=1;c<cols-1;c++){
				temp1 = c1-f[r][c];
				temp2 = c2-f[r][c];
				g = temp1*temp1-temp2*temp2;
				sum = -x[r][c]+bx[r][c]+y[r][c-1]-y[r][c]-by[r][c-1]+by[r][c]-mu/lambda*g;
				u[r][c] = (u[r+1][c]+u[r][c+1]+u[r][c-1]+sum)/3.0;
				if(u[r][c]>1) u[r][c]=1;
				else if(u[r][c]<0) u[r][c]=0;
			}
		r = rows-1;
			for(c=1;c<cols-1;c++){
				temp1 = c1-f[r][c];
				temp2 = c2-f[r][c];
				g = temp1*temp1-temp2*temp2;
				sum = x[r-1][c]-bx[r-1][c]+y[r][c-1]-y[r][c]-by[r][c-1]+by[r][c]-mu/lambda*g;
				u[r][c] = (u[r-1][c]+u[r][c+1]+u[r][c-1]+sum)/3.0;
				if(u[r][c]>1) u[r][c]=1;
				else if(u[r][c]<0) u[r][c]=0;
			}
			
		c = 0;
			for(r=1;r<rows-1;r++){
						temp1 = c1-f[r][c];
						temp2 = c2-f[r][c];
						g = temp1*temp1-temp2*temp2;
						sum = x[r-1][c]-x[r][c]-bx[r-1][c]+bx[r][c]-y[r][c]+by[r][c]-mu/lambda*g;
						u[r][c] = 0.25*(u[r-1][c]+u[r+1][c]+u[r][c+1]+sum);
						if(u[r][c]>1) u[r][c]=1;
						else if(u[r][c]<0) u[r][c]=0;
					}

		c = cols-1;
			for(r=1;r<rows-1;r++){
					temp1 = c1-f[r][c];
					temp2 = c2-f[r][c];
					g = temp1*temp1-temp2*temp2;
					sum = x[r-1][c]-x[r][c]-bx[r-1][c]+bx[r][c]+y[r][c-1]-by[r][c-1]-mu/lambda*g;
					u[r][c] = (u[r-1][c]+u[r+1][c]+u[r][c-1]+sum)/3.0;
					if(u[r][c]>1) u[r][c]=1;
					else if(u[r][c]<0) u[r][c]=0;
				}
			
		r=c=0;{
				temp1 = c1-f[r][c];
				temp2 = c2-f[r][c];
				g = temp1*temp1-temp2*temp2;
				sum = -x[r][c]+bx[r][c]-y[r][c]+by[r][c]-mu/lambda*g;
				u[r][c] = (u[r+1][c]+u[r][c+1]+sum)/2.0;
				if(u[r][c]>1) u[r][c]=1;
				else if(u[r][c]<0) u[r][c]=0;
			}
		
		r=0;c=cols-1;{
				temp1 = c1-f[r][c];
				temp2 = c2-f[r][c];
				g = temp1*temp1-temp2*temp2;
				sum = -x[r][c]+bx[r][c]+y[r][c-1]-by[r][c-1]-mu/lambda*g;
				u[r][c] = (u[r+1][c]+u[r][c-1]+sum)/2.0;
				if(u[r][c]>1) u[r][c]=1;
				else if(u[r][c]<0) u[r][c]=0;
			}
		
		r=rows-1;c=0;{
				temp1 = c1-f[r][c];
				temp2 = c2-f[r][c];
				g = temp1*temp1-temp2*temp2;
				sum = x[r-1][c]-bx[r-1][c]-y[r][c]+by[r][c]-mu/lambda*g;
				u[r][c] = (u[r-1][c]+u[r][c+1]+sum)/2.0;
				if(u[r][c]>1) u[r][c]=1;
				else if(u[r][c]<0) u[r][c]=0;
			}
		
		r=rows-1;c=cols-1;{
				temp1 = c1-f[r][c];
				temp2 = c2-f[r][c];
				g = temp1*temp1-temp2*temp2;
				sum = x[r-1][c]-bx[r-1][c]+y[r][c-1]-by[r][c-1]-mu/lambda*g;
				u[r][c] = (u[r-1][c]+u[r][c-1]+sum)/2.0;
				if(u[r][c]>1) u[r][c]=1;
				else if(u[r][c]<0) u[r][c]=0;
			}
		return change;
}


void gsSpaceseg(num** u,num** edge, num** x, num** y, num** bx, num** by, double lambda, int width, int height){
	int w,h;
	num a,b,s;
	num flux = 1.0/lambda;
	num *uw,*uwp1,*bxw,*byw,*xw,*yw;
    num base;
	for(w=0;w<width-1;w++){	
		uw = u[w];uwp1=u[w+1];bxw=bx[w];byw=by[w];xw=x[w];yw=y[w];
		for(h=0;h<height-1;h++){
			flux = edge[w][h]/lambda;
			a =  uwp1[h]-uw[h]+bxw[h];
			b =  uw[h+1]-uw[h]+byw[h];
			s = a*a+b*b;
			if(s<flux*flux){xw[h]=0;yw[h]=0;continue;}
			s = sqrt(s);
			s=(s-flux)/s;
			xw[h] = s*a;
			yw[h] = s*b;
		}
	}		
	
	h = height-1;
	for(w=0;w<width-1;w++){	
			flux = edge[w][h]/lambda;
			base =  u[w+1][h]-u[w][h]+bx[w][h];
			if(base>flux) {x[w][h] = base-flux; continue;}
			if(base<-flux){x[w][h] = base+flux; continue;}
			x[w][h] = 0;
	}
	w = width-1;
	for(h=0;h<height-1;h++){	
		flux = edge[w][h]/lambda;
		base =  u[w][h+1]-u[w][h]+by[w][h];
		if(base>flux) {y[w][h] = base-flux; continue;}
		if(base<-flux){y[w][h] = base+flux; continue;}
		y[w][h] = 0;
	}
}



void bregmanXseg(num** x,num** u, num** bx, int width, int height){
		int w,h;
		double d;
		num* uwp1,*uw,*bxw,*xw;
		for(w=0;w<width-1;w++){
			uwp1=u[w+1];uw=u[w];bxw=bx[w];xw=x[w];
			for(h=0;h<height;h++){
				d = uwp1[h]-uw[h];
				bxw[h]+= d-xw[h];		
			}
		}
	}


void bregmanYseg(num** y,num** u, num** by, int width, int height){
		int w,h;
		double d;
		int hSent = height-1;
		num* uw,*byw,*yw;
		for(w=0;w<width;w++){
			uw=u[w];byw=by[w];yw=y[w];
			for(h=0;h<hSent;h++){
				d = uw[h+1]-uw[h];
				byw[h]+= d-yw[h];
			}
		}
	}
	
/************************memory****************/

double copyseg(double** source, double** dest, int rows, int cols){
	int r,c;
	double temp,sumDiff;
	for(r=0;r<rows;r++)
		for(c=0;c<cols;c++){
			temp = dest[r][c]- source[r][c];
			sumDiff +=temp*temp;
			
			dest[r][c]=source[r][c];
		
		}
	return sqrt(sumDiff/(rows*cols));
}



num** newMatrixseg(int rows, int cols){
		num* a = (num*) malloc(rows*cols*sizeof(num));
		num** rval = (num**) malloc(rows*sizeof(num*));
		int j,g;
		rval[0] = a;
		for(j=1;j<rows;j++)
			rval[j] = &a[j*cols];
		
		for(j=0;j<rows;j++)
			for(g=0;g<cols;g++)
				rval[j][g] = 0;
		return rval;
}
void deleteMatrixseg(num ** a){
	free(a[0]);
	free(a);
}
