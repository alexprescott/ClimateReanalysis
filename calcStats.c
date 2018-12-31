/*
	This script is used in preparation of pulling CHELSA 20th century time series one at a time and calculating statistics on them using the Welford online single-pass algorithm.
	Alexander Prescott
	December 23 2018
*/

int Nx,Ny;
int **Nprec;
float **precM1, **precM2;

#include<malloc.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

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
	FILE *fr0,*fw0,*fw1,*fw2;
	int i,j,k,month,year;
	float del1,del2;
	char title[50],cmd1[100],cmd2[100],cmd3[50];
	
	Nx = 43200;
	Ny = 20880;
	
	setupmatrices();
	
	short intlines[Nx];
	float fltlines[Nx];

	// Loop over ftp server files starts here
	for (k=500;k<=500;k++)
	{
		// print filename
		month = (k-1)%12 + 1;
		year = (int) floor((k-1)/12) + 1901;
		sprintf(title,"CHELSAcruts_prec_%d_%d_V.1.0.tif",month,year); 

		// Run Bash commands from within this script
		// Note that strcpy and strcat here have not been controlled for null terminator
		// wget
		strcpy(cmd1,"wget -q https://www.wsl.ch/lud/chelsa/data/timeseries20c/prec/");
		strcat(cmd1,title);
		system(cmd1);
		printf("file loaded\n");
		
		// gdal_translate
		strcpy(cmd2,"gdal_translate -ot Int16 -of EHdr ");
		strcat(cmd2,title);
		strcat(cmd2," data.bil");
		system(cmd2);
		printf("file translated\n");

		// Load .bil data and update statistics
		fr0 = fopen("data.bil","rb");
		for (j=1;j<=Ny;j++)
		{
			(void) fread(intlines,sizeof(intlines),1,fr0);
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
		fclose(fr0);
		
		// rm
		strcpy(cmd3,"rm ");
		strcat(cmd3,title);
		strcat(cmd3," data*");
		system(cmd3);
		printf("files deleted\n");
		
		fflush(stdout);
		
	}
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
