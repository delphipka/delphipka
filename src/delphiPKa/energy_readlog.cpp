//
//  energy_readlog.cpp
//  DelPhiPKA
//
//  Created by Lin Wang on 10/22/14.
//  Copyright (c) 2014 Lin. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include "energy.h"


/*
 *  This funciton is used to read energy terms from text files: pairwise.txt and energy.txt generated by energy class.
 *  This is used in case that energy class is set to be skipped but the rest of program needs to run.
 *  Before the clustering class, call this function and the required energy terms are read into the data store class.
 *  Do not use this function if the energy class is set to be run, otherwise it will slow down the runtime.
 */

void CEnergy :: readFromFile()
{
    ionizableArr(); // This function is to genereate a 2d array to store charged residue and their atoms.
    
    
    vecIonRes.erase(vecIonRes.begin()); // vecIonRes[0][n] array is EMPTY. Therefore, remove from the vector
    
    int N = vecIonRes.size();
    int i=0, j=0, dim;
    dim = N * 2;
    EnergyPair.resize(dim);
    for(i=0;i<EnergyPair.size();i++) { EnergyPair[i].resize(dim); }
    
    
    ifstream polarRxnEnergy;
    polarRxnEnergy.open("energies.txt");
    
    string line, strtmp1, strtmp2;
    
    while(getline(polarRxnEnergy, line))
    {
        if(line.substr(0,8) == " RESNAME")  continue;
        
        strtmp1 = line.substr(17,10); remove_leadtail_whitespace_string(strtmp1);
        strtmp2 = line.substr(30,10); remove_leadtail_whitespace_string(strtmp2);
        
        if(line.substr(3,1) == "0") {
//            cout << stof(strtmp1) << "   " << stof(strtmp2) << endl;
            EnergyPolarNeu.push_back(stof(strtmp1) / 0.59);
            EnergyRxnNeu.push_back(stof(strtmp2) / 0.59);
        }
        
        if(line.substr(3,1) != "0") {
//            cout << stof(strtmp1) << "   " << stof(strtmp2) << endl;
            EnergyPolarCrg.push_back(stof(strtmp1) / 0.59);
            EnergyRxnCrg.push_back(stof(strtmp2) / 0.59);
        }
    }
    
    polarRxnEnergy.close();
    
    ifstream pairEnergy;
    pairEnergy.open("pairwise.txt");
    
    Energy ergtemp;
    vector<Energy> globalEnergyArr;
    globalEnergyArr.reserve(dim*dim);
    
    while(getline(pairEnergy, line))
    {
        if(line.substr(0,8) == " RESNAME")  continue;
        
        strtmp1 = line.substr(32,12); remove_leadtail_whitespace_string(strtmp1);
        ergtemp.resnum1 = stoi(line.substr(5,4));
        ergtemp.resnum2 = stoi(line.substr(22,4));
        ergtemp.fEnergy = stof(strtmp1);
        globalEnergyArr.push_back(ergtemp);
        
    }
    
    pairEnergy.close();
    
    for(i=0; i<vecIonRes.size(); i++) {
        for(j=0; j<vecIonRes.size(); j++) {
            
            EnergyPair[2*i  ][2*j  ].resnum1 = globalEnergyArr[ 4*N*i + j*4 + 0 ].resnum1;
            EnergyPair[2*i  ][2*j  ].resnum2 = globalEnergyArr[ 4*N*i + j*4 + 0 ].resnum2;
            EnergyPair[2*i  ][2*j  ].fEnergy = globalEnergyArr[ 4*N*i + j*4 + 0 ].fEnergy;
            
            EnergyPair[2*i][2*j + 1].resnum1 = globalEnergyArr[ 4*N*i + j*4 + 1 ].resnum1;
            EnergyPair[2*i][2*j + 1].resnum2 = globalEnergyArr[ 4*N*i + j*4 + 1 ].resnum2;
            EnergyPair[2*i][2*j + 1].fEnergy = globalEnergyArr[ 4*N*i + j*4 + 1 ].fEnergy;
            
            EnergyPair[2*i + 1][2*j].resnum1 = globalEnergyArr[ 4*N*i + j*4 + 2 ].resnum1;
            EnergyPair[2*i + 1][2*j].resnum2 = globalEnergyArr[ 4*N*i + j*4 + 2 ].resnum2;
            EnergyPair[2*i + 1][2*j].fEnergy = globalEnergyArr[ 4*N*i + j*4 + 2 ].fEnergy;
            
            EnergyPair[2*i + 1][2*j + 1].resnum1 = globalEnergyArr[ 4*N*i + j*4 + 3 ].resnum1;
            EnergyPair[2*i + 1][2*j + 1].resnum2 = globalEnergyArr[ 4*N*i + j*4 + 3 ].resnum2;
            EnergyPair[2*i + 1][2*j + 1].fEnergy = globalEnergyArr[ 4*N*i + j*4 + 3 ].fEnergy;
            
        }
    }
    
#ifdef MPI_PARALLEL
    int id = MPI::COMM_WORLD.Get_rank();
    if (0==id) {
        cout << "Done." << endl; cout << endl;
    }
#endif
    
#ifndef MPI_PARALLEL
    cout << "Done." << endl; cout << endl;
#endif
    
}

