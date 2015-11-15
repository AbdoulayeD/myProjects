//
//  main.cpp
//  Convolution
//
//  Created by DIOP Abdoulaye  on 20/10/2015.
//  Copyright © 2015 DIOP Abdoulaye . All rights reserved.
//

//typedef convolType int;
#include <iostream>
#include "convolution.h"
//#include <mpi.h>

int main(int argc, char ** argv)
{
    //Initialisation
    MPI_Init(&argc,&argv);
    int nrank;
    int nproc;
    int  Q, R;
    MPI_Comm_rank(MPI_COMM_WORLD,&nrank);
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
    std::string outputFilename="output";
    std::string extention=".ppm";
    std::string outputfile =outputFilename+extention;
    //Chargement Image

    Image<int> img("lena.ppm");
       if(nrank==0){
        std::cout<<"Magic Number : "<<img.getmn()<<std::endl;
        std::cout<<"Image width : "<<img.getw()<<std::endl;
        std::cout<<"Image height : "<<img.geth()<<std::endl;
        std::cout<<"Maxval : "<<img.getmv()<<std::endl;
    }

     struct infop infopip;
    infopip.rank  = nrank;
    infopip.nproc = nproc;
    //Convolution

    if(nproc>1)
    {
        Q = img.geth()/nproc;
        R = img.geth()%nproc;

        if (infopip.rank < R) {

            infopip.nloc = Q+1;
            infopip.ideb = infopip.rank * (Q+1);
            infopip.ifin = infopip.ideb + infopip.nloc;

        } else {

            infopip.nloc = Q;
            infopip.ideb = R * (Q+1) + (infopip.rank - R) * Q;
            infopip.ifin = infopip.ideb + infopip.nloc;
        }
        if(infopip.rank==0)
            infopip.gatherdeb=infopip.ifin;

        Convolution<float> convol(img,sharpen<float>(2), infopip);
         if(infopip.rank==0)
            convol.save(outputfile,infopip);
    }
    else
    {

        Convolution<float> convol(img,motionblur<float>());
        convol.save(outputfile,infopip);
    }

    MPI_Finalize();

    return 0;
}

