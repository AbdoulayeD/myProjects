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
#include <mpi.h>

struct infop{
    int ideb;
    int ifin;
    int nloc;
};

int main(int argc, char ** argv) {
    //Initialisation
    MPI_Init(&argc,&argv);
    int nrank;
    int nsize;
    MPI_Comm_rank(MPI_COMM_WORLD,&nrank);
    MPI_Comm_size(MPI_COMM_WORLD,&nsize);
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
    struct infop infopip;
    infopip.nloc=img.getw()/nsize;
    infopip.ideb=nrank*infopip.nloc;
    infopip.ifin=infopip.ideb+infopip.nloc;

    //Convolution
    double fact=1.0/3.0;
    double bia=0.0;
    //Convolution<int> convol(img);
    //Convolution<float> convol(img,fact, bia);
    //Convolution<int> convol(img,edges<int>(9));
    Convolution<float> convol(img,motionblur<float>(9),fact,bia);
    convol.save(outputfile);
    MPI_Finalize();
    return 0;
}
