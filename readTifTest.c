/*
 This C script is a test of using GDAL to read data from a GeoTIFF *.tif file
 Alexander Prescott
 December 22 2018
*/

#include<malloc.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<gdal.h>
#include<cpl_conv.h> /*for CPLMalloc() */

#define NR_END 1
#define FREE_ARG char*
#define PI 3.1415926


int Nx,Ny;
int **tifData;

void nrerror(error_text)
char error_text[];
/* Numerical Recipes standard error handler */
{
        void exit();

        fprintf(stderr,"Numerical Recipes run-time error...\n");
        fprintf(stderr,"%s\n",error_text);
        fprintf(stderr,"...now exiting to system...\n");
        exit(1);
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
        free((FREE_ARG) (v+nl-NR_END));
}

void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
   free((FREE_ARG) (v+nl-NR_END));
}   

void free_imatrix(int **m, long nrl,long nrh,long ncl,long nch)
/* free an int matrix allocated by imatrix() */
{
   free((FREE_ARG) (m[nrl]+ncl-1));
   free((FREE_ARG) (m+nrl-1));
}

void free_matrix(float **m, long nrl,long nrh,long ncl,long nch)
/* free a float matrix allocated by matrix() */
{
   free((FREE_ARG) (m[nrl]+ncl-1));
   free((FREE_ARG) (m+nrl-1));
} 

float *vector(nl,nh)
long nh,nl;
/* allocate a float vector with subscript range v[nl..nh] */
{
        float *v;

        v=(float *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(float)));
        if (!v) nrerror("allocation failure in vector()");
        return v-nl+NR_END;
}

int *ivector(nl,nh)
long nh,nl;
/* allocate an int vector with subscript range v[nl..nh] */
{
        int *v;

        v=(int *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(int)));
        if (!v) nrerror("allocation failure in ivector()");
        return v-nl+NR_END;
}

long *lvector(nl,nh)
long nh,nl;
/* allocate a long int vector with subscript range v[nl..nh] */
{
        long *v;

        v=(long *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(long)));
        if (!v) nrerror("allocation failure in ivector()");
        return v-nl+NR_END;
}

float **matrix(nrl,nrh,ncl,nch)
long nch,ncl,nrh,nrl;
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        float **m;

        /* allocate pointers to rows */
        m=(float **) malloc((unsigned int)((nrow+NR_END)*sizeof(float*)));
        if (!m) nrerror("allocation failure 1 in matrix()");
        m += NR_END;
        m -= nrl;

        /* allocate rows and set pointers to them */
        m[nrl]=(float *) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(float)));
        if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
        m[nrl] += NR_END;
        m[nrl] -= ncl;

        for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

        /* return pointer to array of pointers to rows */
        return m;
}

int **imatrix(nrl,nrh,ncl,nch)
long nch,ncl,nrh,nrl;
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        int **m;

        /* allocate pointers to rows */
        m=(int **) malloc((unsigned int)((nrow+NR_END)*sizeof(int*)));
        if (!m) nrerror("allocation failure 1 in matrix()");
        m += NR_END;
        m -= nrl;


        /* allocate rows and set pointers to them */
        m[nrl]=(int *) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(int)));
        if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
        m[nrl] += NR_END;
        m[nrl] -= ncl;

        for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

        /* return pointer to array of pointers to rows */
        return m;
}


void setupmatrices()
{    //int i,j;
	 
	 tifData=imatrix(1,Nx,1,Ny);

}

int main ()
{
	
	Nx = 43200;
	Ny = 20880;
	
	/* Prepare GDAL dataset */
	GDALDatasetH hDataset;
	GDALAllRegister();
	hDataset = GDALOpen("CHELSAcruts_prec_1_1901_V.1.0.tif",GA_ReadOnly);
	if (hDataset == NULL)
	{
		printf("DATASET LOAD FAILURE\n");
		return 1;
	}
	
	/* Do some Gdal Stuff */
	GDALDriverH   hDriver;
	double        adfGeoTransform[6];
	hDriver = GDALGetDatasetDriver( hDataset );
	printf( "Driver: %s/%s\n",
			GDALGetDriverShortName( hDriver ),
			GDALGetDriverLongName( hDriver ) );
	printf( "Size is %dx%dx%d\n",
			GDALGetRasterXSize( hDataset ),
			GDALGetRasterYSize( hDataset ),
			GDALGetRasterCount( hDataset ) );
	if( GDALGetProjectionRef( hDataset ) != NULL )
		printf( "Projection is `%s'\n", GDALGetProjectionRef( hDataset ) );
	if( GDALGetGeoTransform( hDataset, adfGeoTransform ) == CE_None )
	{
		printf( "Origin = (%.6f,%.6f)\n",
				adfGeoTransform[0], adfGeoTransform[3] );
		printf( "Pixel Size = (%.6f,%.6f)\n",
				adfGeoTransform[1], adfGeoTransform[5] );
	}
	
	/*
	for (j=1;j<=Ny;j++)
		if (tifData[j][j] >= 0)
			printf("%d \n",tifData[j][j]);
	*/
	
	return 0;
}