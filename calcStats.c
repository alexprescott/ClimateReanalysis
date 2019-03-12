/*
	This script is used in preparation of pulling CHELSA 20th century time series one at a time and calculating statistics on them using the Welford online single-pass algorithm.
	
	This version uses GDAL API and iterates on the downloaded CHELSA data sets. These files must be pre-downloaded using wget or something similar. The intention is to build a Docker container that compiles this script with the relevant library dependencies to allow for execution on the UA HPC system. Notably, the GDAL developer library libgdal-dev is required, linked upon compilation with the -lgdal flag.
	
	NOTE: wget stores files in /extra/alexprescott/prec20c; first line of bash script must be to copy all of the files from that directory to /tmp/, the temporary directory on the HPC execute node disk (maximum storage of ~800 GB). The /tmp/ directory is incrementally filled as users run their jobs and store files on it; the HPC manager automatically suspends jobs sent to a given node and clears the /tmp directory once it is filled past some threshold level. If the job fails and returns an error relating to not enough disk space, it is probably because /tmp was almost full and the requirements of my job exceeded the remaining available space. Note that I can not edit /tmp locations that other users have written to, even if their jobs have run to completion.
	
	Alexander Prescott
	March 12 2019
*/

int Nx,Ny;
int **Nprec;
float **precM1, **precM2;

#include<malloc.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<gdal.h>
#include<cpl_conv.h> /*for CPLMalloc() */

#define NR_END 1
#define FREE_ARG char*

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
{    int i,j;
	 
	 Nprec	=imatrix(1,Nx,1,Ny);
	 precM1	=matrix(1,Nx,1,Ny);
	 precM2	=matrix(1,Nx,1,Ny);
	 
	for (j=1;j<=Ny;j++)
	{
		for(i=1;i<=Nx;i++)
		{
			Nprec[i][j]=0;
			precM1[i][j]=0.0;
			precM2[i][j]=0.0;
		}
	}

}



int main ()
{
	FILE *fw0,*fw1,*fw2;
	int i,j,k,month,year;
	float del1,del2;
	char title[50];
	
	Nx = 43200;
	Ny = 20880;
	
	setupmatrices();
	
	short intlines[Nx];
	float fltlines[Nx];

	/* Initialize GDAL dataset */
	GDALDatasetH hDataset;
	GDALAllRegister();
	GDALRasterBandH hBand;

	/* Loop over all CHELSA 20th century reanalysis precipitation files starts here */
	for (k=1;k<=1392;k++)
	{
		/* print filename into title */
		month = (k-1)%12 + 1;
		year = (int) floor((k-1)/12) + 1901;
		sprintf(title,"/tmp/CHELSAcruts_prec_%d_%d_V.1.0.tif\0",month,year); // \0 is null string terminator
		
		/* Load dataset */
		hDataset = GDALOpen(title,GA_ReadOnly);
		if (hDataset == NULL)
		{
			printf("DATASET LOAD FAILURE ON ITERATION %d\n",k);
			return 1;
		}
		
		/* Choose relevant data band */
		hBand = GDALGetRasterBand( hDataset, 1 );
		
		// Load .tif data and update statistics
		for (j=1;j<=Ny;j++)
		{
			/* Read .tif data from hBand, with 0 x-offset and j-1 y-offset, read Nx entries into intlines, of type Int16 i.e. short integers; I believe this is parallelizable. */
			GDALRasterIO( hBand, GF_Read, 0, j-1, Nx, 1,intlines, Nx, 1, GDT_Int16,0, 0 );
			for (i=1;i<=Nx;i++)
			{
				if (intlines[i-1] >= 0) // if there is precip data, then do:
				{// Welford online algorithm
					Nprec[i][j]++;
					del1 = intlines[i-1] - precM1[i][j];
					precM1[i][j]+= del1/Nprec[i][j];
					del2 = intlines[i-1] - precM1[i][j];
					precM2[i][j]+= del1*del2;
				}
			}
		}
		
		if (k%40 == 20)
		{
			printf("Iteration: %d   Month: %d   Year: %d \n",k,month,year);
			fflush(stdout);
		}
		
	} /* end file load loop */
	
	printf("File loading and statistics done\n");
	
	// write statistics arrays to file
	fw0 = fopen("Chelsa20Cmean.flt","wb");
	fw1 = fopen("Chelsa20Ccv.flt","wb");
	fw2 = fopen("Chelsa20Ccount.bil","wb");
	for (j=1;j<=Ny;j++)
	{
		for (i=1;i<=Nx;i++)
			fltlines[i-1]=precM1[i][j];
		(void)fwrite(fltlines,sizeof(fltlines),1,fw0); // mean
		
		for (i=1;i<=Nx;i++)
			if (Nprec[i][j] > 2)
				fltlines[i-1]=sqrt(precM2[i][j]/(Nprec[i][j]-1))/precM1[i][j]; // coefficient of variation
		(void)fwrite(fltlines,sizeof(fltlines),1,fw1);
		
		for (i=1;i<=Nx;i++)
			intlines[i-1]=Nprec[i][j]; // number of records
		(void)fwrite(intlines,sizeof(intlines),1,fw2);
	}
	fclose(fw0);
	fclose(fw1);
	fclose(fw2);

	return 0;
}
