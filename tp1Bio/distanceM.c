#include<stdio.h>
#include<unistd.h>
#include<stdlib.h>
#include<string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <limits.h>

#define size 5
struct Info{
	float min;
	int i_;
	int j_;
};


struct Info minMat(float * A, int n)
{
	int i,j;
	struct Info info;
	info.min = INT_MAX;

	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
			if (A[i*n+j]<info.min && A[i*n+j]!=0)
			{
				info.min=A[i*n+j];
				info.i_=i;
				info.j_=j;
			}

return info;
}

float * minV(float * v1, float * v2, int m)
{
	float * minvec;
	int i;
	minvec = malloc(m*sizeof(float));

	for (i=0;i<m;i++)
	{
	    if (v1[i]<=v2[i]) minvec[i] = v1[i];
		else minvec[i] = v2[i];
	}
return minvec;
}



int id(int i, int j, int sz)
{return i*sz+j;}

float rSize(char * inName)
{
	float sizeStat;
	struct stat stime;
	stat(inName, &stime);
	sizeStat=(float)stime.st_size;
	return sizeStat;
}

float d(char* fn0, char *fn1, int i, int j)
{
    FILE* fp;
	int retval;
	char str[100];
	float ta;
	char  zName[5][10];
    char * createZ[2];
    char * createZP[5];

	//if(strcmp(fn0,fn1)!=0)
	if(i!=j)
	{
        //Z(A0),Z(A1)

        sprintf(zName[0],"%d.zip",i);
		sprintf(str,"7z a -tzip %s %s -mx9",zName[0],fn0);
		printf("%d:%s\n",0,str);

        fp =fopen(zName[0],"r");
        if (fp==NULL)
            retval=system(str);
        else
            fclose(fp);

		if(retval==127)
		printf("Execution Failed\n");

        sprintf(zName[1],"%d.zip",j);

		sprintf(str,"7z a -tzip %s %s -mx9",zName[1],fn1);
		printf("%d:%s\n",1,str);

        fp =fopen(zName[1],"r");
        if (fp==NULL)
            retval=system(str);
        else
            fclose(fp);

		if(retval==127)
		printf("Execution Failed\n");
		//Z(A0&A0),Z(A1&A1)

		//A0
		sprintf(str,"cat %s >> temp.fna ; cat %s >> temp.fna",fn0,fn0);
		retval=system(str);
		 if(retval==127)
			printf("Execution Failed\n");

        sprintf(zName[2],"%d%d.zip",i,i);
		sprintf(str,"7z a -tzip %s temp.fna -mx9",zName[2]);
		//printf("%d:%s\n",i,str);

		 fp =fopen(zName[2],"r");
        if (fp==NULL)
            retval=system(str);
        else
            fclose(fp);

		if(retval==127)
			printf("Execution Failed\n");

		sprintf(str,"rm temp.fna");
		retval=system(str);
		 if(retval==127)
			printf("Execution Failed\n");
		//A1
		sprintf(str,"cat %s >> temp.fna ; cat %s >> temp.fna",fn1,fn1);
		retval=system(str);
		 if(retval==127)
			printf("Execution Failed\n");

        sprintf(zName[3],"%d%d.zip",j,j);
		sprintf(str,"7z a -tzip %s temp.fna -mx9",zName[3]);
		//printf("%d:%s\n",i,str);

		 fp =fopen(zName[3],"r");
        if (fp==NULL)
            retval=system(str);
        else
            fclose(fp);
		if(retval==127)
			printf("Execution Failed\n");

		sprintf(str,"rm temp.fna");
		retval=system(str);
		 if(retval==127)
			printf("Execution Failed\n");

        //Z(A0&A1)
        sprintf(zName[4],"%d%d.zip",i,j);
		sprintf(str,"7z a -tzip %s %s %s -mx9",zName[4],fn0,fn1);
		//printf("%d,%d,%s\n",i,j,str);
		fp =fopen(zName[4],"r");
        if (fp==NULL)
            retval=system(str);
        else
            fclose(fp);
		if(retval==127)
		    printf("Execution Failed\n");


        //sprintf(zName[0],"%d.zip",i);
        createZP[0]=zName[0];


        //sprintf(zName[1],"%d.zip",j);
        createZP[1]=zName[1];

        //sprintf(zName[2],"%d%d.zip",i,i);
        createZP[2]=zName[2];

        //sprintf(zName[3],"%d%d.zip",j,j);
        createZP[3]=zName[3];

        //sprintf(zName[4],"%d%d.zip",i,j);
        createZP[4]=zName[4];

        ta = rSize( createZP[4] ) / ( rSize(createZP[0]) + rSize( createZP[1])  )
            - rSize( createZP[2] ) / (4*rSize( createZP[0] ))
            - rSize( createZP[3] ) / (4*rSize(createZP[1]));


	}
	else
	{

        sprintf(zName[0],"%d.zip",i);

		sprintf(str,"7z a -tzip %s %s -mx9",zName[0],fn0);
		printf("%d:%s\n",0,str);
		fp =fopen(zName[0],"r");
        if (fp==NULL)
            retval=system(str);
        else
            fclose(fp);
		if(retval==127)
		printf("Execution Failed\n");


	//A0&A0
		sprintf(str,"cat %s >> temp.fna ; cat %s >> temp.fna",fn0,fn0);
		retval=system(str);
	 	if(retval==127)
	    	printf("Execution Failed\n");

        sprintf(zName[1],"%d%d.zip",i,i);
		sprintf(str,"7z a -tzip %s temp.fna -mx9",zName[1]);
		//printf("%d:%s\n",i,str);

		 fp =fopen(zName[1],"r");
        if (fp==NULL)
            retval=system(str);
        else
            fclose(fp);
		if(retval==127)
	    	printf("Execution Failed\n");

		sprintf(str,"rm temp.fna");
		retval=system(str);
	 	if(retval==127)
	    	printf("Execution Failed\n");



        //sprintf(zName[0],"%d.zip",i);
        //printf("file0=%s\n",zName[0]);
        createZ[0]=zName[0];


        //sprintf(zName[1],"%d%d.zip",i,i );
        //printf("file0=%s\n",zName[1]);
        createZ[1]=zName[1];

        printf("file0=%s\n",createZ[0]);
        printf("file1=%s\n",createZ[1]);

        printf("size0=%f\n",rSize(createZ[0]));
        printf("size1=%f\n",rSize(createZ[1]));

        ta = rSize( createZ[1] ) / ( rSize(createZ[0]) + rSize( createZ[0]) )
            - rSize( createZ[1] ) / (4*rSize(createZ[0]))
            - rSize( createZ[1] ) / (4*rSize(createZ[0]));

	}

	return ta;
}


int main(int argc, char** argv)
{
    int n;
    int m;
    int i,j,cpt=0;
    int retval;
    float A[size*size];
    char* fnaName[size];
    char str[100];
    float dab;

	for(i=0;i<size;i++){
        fnaName[i]=argv[i+1];
		printf(" Test size:%f\n", rSize(fnaName[i]) );
    }


	for(i=0;i<size;i++)
    for(j=0;j<size;j++)
    {
        dab=d(fnaName[i],fnaName[j],i,j);
        A[id(i,j,size)] = dab;
        printf("Distance: %f \n",dab);
    }
    //dab=d(fnaName[0],fnaName[1],0,1);

    printf("\n Test A:\n");
    for(i=0;i<size;i++){
    for(j=0;j<size;j++)
    {
        printf("%f ",A[id(i,j,size)]);
    }
    printf("\n");
    }



    n=size;
 	
 	m = n-1;
	while (m>0)
	{
		struct Info minA = minMat(A,n);

		printf("i_:%d,j_:%d \n",minA.i_, minA.j_);
		float B[m*m];

		for(i=0;i<minA.j_;i++)
		for(j=0;j<minA.j_;j++)
		B[i*m+j] = A[i*n+j];

		for(i=0;i<minA.j_;i++)
		for(j=minA.j_;j<m;j++)
		B[i*m+j] = A[i*n+(j+1)];

		for(i=minA.j_;i<m;i++)
		for(j=0;j<minA.j_;j++)
		B[i*m+j] = A[(i+1)*n+j];

		for(i=minA.j_;i<m;i++)
		for(j=minA.j_;j<m;j++)
		B[i*m+j] = A[(i+1)*n+(j+1)];


		for(i=0;i<m;i++)
		{
			for(j=0;j<m;j++)
				printf("B : %f", B[i*m+j]);
		 	printf("\n");
		}

		float * v;
		v = (float *)malloc(n*sizeof(float));
		float v1[n],v2[n];

		for(i=0;i<n;i++)
		{
			v1[i]=A[minA.i_*n+i];
			v2[i]=A[i*n+minA.j_];
		}

		for(i=0;i<n;i++)
		{
			printf("V1[%d] : %f  ",i,v1[i]);
		}

		printf("\n");

		for(i=0;i<n;i++)
		{
			printf("V2[%d] : %f  ",i,v2[i]);
		}

		printf("\n");
		v = minV(v1,v2,n);
		for(i=0;i<n;i++)
		{
			printf("V[%d] : %f  ",i,v[i]);
		}

		double v_[m];
		for (i=0;i<n;i++)
		{
			if(i<minA.j_)
			v_[i] = v[i];

			else if (i>=minA.j_)
			v_[i] = v[i+1];
		}

		for (i=0;i<m;i++)
			printf("v_ : %f  ",v_[i]);
		printf("\n");

		for (i=0;i<m;i++)
		{
			B[i*m+minA.i_] = v_[i];
			B[minA.i_*m+i] = v_[i];
		}

		printf("\n");
		for(i=0;i<m;i++)
		{
			for(j=0;j<m;j++)
			{
				printf(" %f ",B[i*m+j]);
			}
			printf("\n");
		}

		for (i=0;i<m;i++)
		for (j=0;j<m;j++)
		{
			A[i*m+j] = B[i*m+j];
		}

		n = m;
		m -= 1;
	}
    

	return 0;
}

