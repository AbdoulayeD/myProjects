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

int main(int argc, char ** argv) {
    //Initialisation
    MPI_Init(&argc,&argv);
    int nrank;
    int nproc;
    int  Q, R;
    MPI_Comm_rank(MPI_COMM_WORLD,&nrank);
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
    std::string outputFilename="output";
    std::string extention=".pgm";
    std::string outputfile =outputFilename+extention;
    //std::cout<<std::endl<<"I am process number:"<<nrank<<std::endl;
    //Chargement Image
    Image<int> img("lena.ascii.pgm");
       if(nrank==0){
        std::cout<<"Image width : "<<img.getw()<<std::endl;
        std::cout<<"Image height : "<<img.geth()<<std::endl;
    }
    
    Q = img.geth()/nproc;
    R = img.geth()%nproc;
    struct infop infopip;
    infopip.rank  = nrank;
    infopip.nproc = nproc;
    
    /*infopip.nloc  = img.getw()/nproc;
    infopip.ideb  = nrank*infopip.nloc;
    infopip.ifin  = infopip.ideb+infopip.nloc;
    */
    if (infopip.rank < R) {
        
        infopip.nloc = Q+1;
        infopip.ideb =infopip.rank * (Q+1);
        infopip.ifin =infopip.ideb + infopip.nloc;
        
    } else {
        
        infopip.nloc = Q;
        infopip.ideb =R * (Q+1) + (infopip.rank - R) * Q;
        infopip.ifin =infopip.ideb + infopip.nloc;
    }


    //Convolution
    double fact=1.0;
    double bia=0.0;
    //Convolution<int> convol(img);
    //Convolution<float> convol(img,fact, bia);
    //Convolution<int> convol(img,edges<int>(9));
    //Convolution<float> convol(img,sharpen<float>(3),fact,bia);
    Convolution<float> convol(img,edges<float>(3),fact,bia,infopip);
    if(nrank == 0)
        convol.save(outputfile);
    
    MPI_Finalize();
    return 0;
}
