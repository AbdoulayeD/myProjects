#include<stdio.h>
#include<unistd.h>
#include<stdlib.h>
#include<string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#define size 2

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
	int retval;
	char str[100];
	float ta;
	char  zName[5][10];
    char * createZ[2];
    char * createZP[5];

	if(strcmp(fn0,fn1)!=0)
	{
        //Z(A0),Z(A1)
		sprintf(str,"7z a -tzip %d.zip %s -mx9",i,fn0);
		printf("%d:%s\n",0,str);
		retval=system(str);
		if(retval==127)
		printf("Execution Failed\n");

		sprintf(str,"7z a -tzip %d.zip %s -mx9",j,fn1);
		printf("%d:%s\n",1,str);
		retval=system(str);
		if(retval==127)
		printf("Execution Failed\n");
		//Z(A0&A0),Z(A1&A1)

		//A0
		sprintf(str,"cat %s >> temp.fna ; cat %s >> temp.fna",fn0,fn0);
		retval=system(str);
		 if(retval==127)
			printf("Execution Failed\n");

		sprintf(str,"7z a -tzip %d%d.zip temp.fna -mx9",i,i);
		//printf("%d:%s\n",i,str);

		retval=system(str);
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

		sprintf(str,"7z a -tzip %d%d.zip temp.fna -mx9",j,j);
		//printf("%d:%s\n",i,str);

		retval=system(str);
		if(retval==127)
			printf("Execution Failed\n");

		sprintf(str,"rm temp.fna");
		retval=system(str);
		 if(retval==127)
			printf("Execution Failed\n");

        //Z(A0&A1)
		sprintf(str,"7z a -tzip %d%d.zip %s %s -mx9",i,j,fn0,fn1);
		//printf("%d,%d,%s\n",i,j,str);
		retval=system(str);
		if(retval==127)
		    printf("Execution Failed\n");


        sprintf(zName[0],"%d.zip",i);
        createZP[0]=zName[0];


        sprintf(zName[1],"%d.zip",j);
        createZP[1]=zName[1];

        sprintf(zName[2],"%d%d.zip",i,i);
        createZP[2]=zName[2];

        sprintf(zName[3],"%d%d.zip",j,j);
        createZP[3]=zName[3];

        sprintf(zName[4],"%d%d.zip",i,j);
        createZP[4]=zName[4];

        ta = rSize( createZP[4] ) / ( rSize(createZP[0]) + rSize( createZP[1])  )
            - rSize( createZP[2] ) / (4*rSize( createZP[0] ))
            - rSize( createZP[3] ) / (4*rSize(createZP[1]));


	}
	else
	{
		sprintf(str,"7z a -tzip %d.zip %s -mx9",i,fn0);
		printf("%d:%s\n",0,str);
		retval=system(str);
		if(retval==127)
		printf("Execution Failed\n");


	//A0&A0
		sprintf(str,"cat %s >> temp.fna ; cat %s >> temp.fna",fn0,fn0);
		retval=system(str);
	 	if(retval==127)
	    	printf("Execution Failed\n");

		sprintf(str,"7z a -tzip %d%d.zip temp.fna -mx9",i,i);
		//printf("%d:%s\n",i,str);

		retval=system(str);
		if(retval==127)
	    	printf("Execution Failed\n");

		sprintf(str,"rm temp.fna");
		retval=system(str);
	 	if(retval==127)
	    	printf("Execution Failed\n");



        sprintf(zName[0],"%d.zip",i);
        //printf("file0=%s\n",zName[0]);
        createZ[0]=zName[0];


        sprintf(zName[1],"%d%d.zip",i,i );
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
    int sz;
    int rk;
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
        //if(i!=0j){
            dab=d(fnaName[i],fnaName[j],i,j);
            A[id(i,j,size)] = dab;
            printf("Distance: %f \n",dab);
        //}
    }


    //test 2 loop


    /*
    for(i=0;i<size;i++)
    {
        dab=d(fnaName[i],fnaName[i],i,i);
        A[id(i,i,size)] = dab;
        printf("Distance: %f \n",dab);
    }

    for(i=0;i<size-1;i++)
    {
        dab=d(fnaName[i],fnaName[i+1],i,i+1);

        A[id(i,i+1,size-1)] = dab;
        printf("Distance: %f \n",dab);
    }
    */



    printf("\n Test A:\n");
    for(i=0;i<size;i++){
    for(j=0;j<size;j++)
    {
        printf("%f ",A[id(i,j,size)]);
    }
    printf("\n");
    }
     //dab=d(fnaName[0],fnaName[1],0,1);
        //printf("Distance: %f \n",dab);


	//Code pour les 15
    /*

    for(i=0;i<15;i++)
    {
        sprintf(str,"7z a -tzip %s.zip %s -mx9",fnaName[i],fnaName[i]);
        cpt++;
        printf("%d:%s\n",i,str);
        retval=system(str);
        if(retval==127)
        printf("Execution Failed\n");
    }

	for(i=0;i<15;i++)
    {
        sprintf(str,"cat %s >> temp.fna ; cat %s >> temp.fna",fnaName[i],fnaName[i]);
        retval=system(str);
         if(retval==127)
            printf("Execution Failed\n");

        sprintf(str,"7z a -tzip %d_%d.zip temp.fna -mx9",i, i);
        //printf("%d:%s\n",i,str);

        retval=system(str);
        if(retval==127)
            printf("Execution Failed\n");

        sprintf(str,"rm temp.fna");
        retval=system(str);
         if(retval==127)
            printf("Execution Failed\n");
    }

    for(i=0;i<15;i++)
    for(j=0;j<15;j++)
    {
		if(i!=j)
		{
		    sprintf(str,"7z a -tzip %d_%d.zip %s %s -mx9",i,j,fnaName[i],fnaName[j]);
		    retval=system(str);
		    if(retval==127)
		        printf("Execution Failed\n");
		}
    }
   */


	return 0;
}
