//
//  main.cpp
//  Convolution
//
//  Created by DIOP Abdoulaye  on 20/10/2015.
//  Copyright © 2015 DIOP Abdoulaye . All rights reserved.
//

//typedef convolType int;
#include <iostream>
#include<cmath>
#include <omp.h>
#include "convolution.h"
#define rep 1
//#include <mpi.h>

int main(int argc, char ** argv)
{
    //Initialisation MPI
    MPI_Init(&argc,&argv);
    int rank;
    int nproc;
    int nthread;
    int  Q, R;
    std::string inputFilename;
    std::string outputFilename="output";
    std::string extention;
    std::string outputfile;
    double t0=0.0,t1=0.0,dt=0.0;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);

    // Arguments d'entrées
    if(argc<3 ){
        /*
        std::cerr<<"Not Enough input Arguments.You must specify:"
        <<std::endl<<"->The Filename"<<std::endl<<"->the extension (ppm or pgm)"
        <<std::endl<<"...."<<std::endl;
        return EXIT_FAILURE;
        */
        if (rank==0){
        std::cout<<"Default Mode:"
                 <<"Convolution pgm"
                 <<std::endl
                 <<"Image:lena.pgm"
                 <<std::endl;
        }
        extention=".pgm";

        outputfile =outputFilename+extention;
        inputFilename="./img/lena.pgm";
    }
    else
    {
        if ( static_cast<std::string>(argv[2]) == "ppm")
            extention=".ppm";
        if ( static_cast<std::string>(argv[2]) == "pgm")
            extention=".pgm";
      outputfile =outputFilename+extention;
        inputFilename=static_cast<std::string>(argv[1]);
    }
    //Chargement Image
        Image<int> img(inputFilename);
    // Préparation du fichier de sorties


       if(rank==0){
        std::cout<<"-->Input Image Carateristics:"<<std::endl;
        std::cout<<"Magic Number : "<<img.getmn()<<std::endl;
        std::cout<<"Image width : "<<img.getw()<<std::endl;
        std::cout<<"Image height : "<<img.geth()<<std::endl;
        std::cout<<"Maxval : "<<img.getmv()<<std::endl;
        std::cout<<std::endl;
    }
    struct infop infopip;
    infopip.rank  = rank;
    infopip.nproc = nproc;

    // Lancement de la Convolution
    if(nproc>1)
    {
        Q = ceil(img.geth()/nproc);
        R = (img.geth()%nproc);
        infopip.nloc = Q;
        infopip.ideb = infopip.rank * infopip.nloc;
        infopip.ifin = infopip.ideb + infopip.nloc;

        
        /*
        if( img.geth() % nproc == 0)
        {
            infopip.nloc = Q;
            infopip.ideb = infopip.rank * infopip.nloc;
            infopip.ifin = infopip.ideb + infopip.nloc;
        }
        
        else
        {
            if (infopip.rank < R) {
                infopip.nloc = Q+1;
                infopip.ideb = infopip.rank * infopip.nloc;
                infopip.ifin = infopip.ideb + infopip.nloc;

            } else {

                infopip.nloc = Q;
                infopip.ideb = R * (Q+1) + (infopip.rank-R) * Q;
                infopip.ifin = infopip.ideb + infopip.nloc;
            }

        }*/
        
        

        if(rank==0){
            std::cout<<"-->MPI Status:"<<std::endl;
            std::cout<<"Number of process :"<<nproc<<std::endl;
            std::cout<<std::endl;
        }
        Convolution<int> convol(img,sharpen<int>(2), infopip);

        MPI_Barrier(MPI_COMM_WORLD);

        for(int rc=0 ; rc < rep ; rc++)
        {   t0=MPI_Wtime();
            Convolution<float> convol(img,sharpen<float>(2), infopip);
            t1=MPI_Wtime();
            dt+=t1-t0;
        }
        if(infopip.rank==0){
                std::cout<<"Parallel Convolution Time Elapsed:"<<dt/rep<<std::endl;
                convol.save(outputfile,infopip);
                std::cout<<"Convolution and Saving Done!"<<std::endl;
            }
    }
    else
    {

        Convolution<int> convol(img,motionblur<int>());
        for(int rc = 0 ;  rc< rep ; rc++)
        {
            t0=MPI_Wtime();
            Convolution<float> convol(img,motionblur<float>());
                 t1=MPI_Wtime();
            dt+=t1-t0;
        }

        if(infopip.rank==0){
            std::cout<<"Sequential Convolution Time Elapsed:"<<dt/rep<<std::endl;
            convol.save(outputfile,infopip);
            std::cout<<"Convolution and Saving Done!"<<std::endl;
        }
    }

    MPI_Finalize();

    return 0;
}

