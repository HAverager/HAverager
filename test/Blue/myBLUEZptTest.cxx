//-----------------------------------------------------------------------------
//
// BLUE: A ROOT class implementing the Best Linear Unbiased Estimate method.
// 
// Copyright (C) 2012-2014, Richard.Nisius@mpp.mpg.de
// All rights reserved
//
// This file is part of BLUE - Version 2.1.0.
//
// BLUE is free software: you can redistribute it and/or modify it under the 
// terms of the GNU Lesser General Public License as published by the Free 
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// For the licensing terms see the file COPYING or http://www.gnu.org/licenses.
//
//-----------------------------------------------------------------------------
#include "Blue.h"
#include <sstream>
#include <iostream>
#include <fstream>

using namespace std;

void myBLUEZptTest(){

	// The number of estimates, uncertainties and observables
	static const Int_t NumEst = 2*43;
	static const Int_t NumUnc = 2;
	static const Int_t NumObs = 43;

	// Index for which estimates determines which observable
	Int_t IWhichObs[NumEst];

	for(int i=0; i<NumEst; i++){
		IWhichObs[i] = i;
		if(i>=NumObs){
			IWhichObs[i] = i-NumObs;
		}
		printf("%i %i \n", i, IWhichObs[i]);
	}

  static const Int_t LenXEst = NumEst * (NumUnc+1);
  Double_t XEst[LenXEst];
  Double_t XCorr[NumEst*NumEst];
  Double_t BlueResult[LenXEst] = {0};
  Double_t BlueCov[NumObs*NumObs] = {0};

  //Reading text file
  string file1 = "Zpt/ZeePtBlue2.txt";
  string file2 = "Zpt/ZmmPtBlue2.txt";
  string file3 = "Zpt/ZPtBlueCor.txt";

  std::ifstream filestream(file1.c_str(), std::ios::in);
  int i=0;
  while(filestream.good()) {
	filestream >> XEst[i];
	i++;
  }

  i--;
  std::ifstream filestream2(file2.c_str(), std::ios::in);
  while(filestream2.good()) {
	filestream2 >> XEst[i];
	i++;
  }

  std::ifstream filestream3(file3.c_str(), std::ios::in);
  i=0;
  while(filestream3.good()) {
	filestream3 >> XCorr[i];
	i++;
  }

  for(int j=0; j<LenXEst; j++){
    if(!(j%(NumUnc+1))){
      printf("\n");
    }
     printf("%0.4e ", XEst[j]);
  }

  // Construct Object
  Blue *myBlue = new Blue(NumEst, NumUnc, NumObs, &IWhichObs[0]);
  myBlue->SetPrintLevel(1);
  myBlue->PrintStatus();

  // Fill estimates
  Int_t ind = 0;
  for(Int_t i = 0; i<NumEst; i++){
    myBlue->FillEst(i,&XEst[ind]);
    ind = ind + NumUnc + 1;
  }

	// Fill correlations
	myBlue->FillCor(0,0.0);
	myBlue->FillCor(1,XCorr);


	// Calculate
    myBlue->FixInp();
    myBlue->PrintEst(0);
    myBlue->PrintEst(1);
    myBlue->PrintRho();
    myBlue->Solve();
    myBlue->PrintResult();

	myBlue->GetResult(BlueResult);
	myBlue->GetCovRes(BlueCov);

	// Prepare output files
	ofstream myfile;
	myfile.open ("ZptBlueOut.txt");
  	for(Int_t k = 0; k<NumObs; k++){
  		for(Int_t i = 0; i<(NumUnc+1); i++){
			myfile << BlueResult[i + k*(NumUnc+1) ] << "  ";
		}
		myfile << "\n";
	}
	myfile.close();


	myfile.open ("ZptBlueCov.txt");
  	for(Int_t k = 0; k<NumObs; k++){
  		for(Int_t i = 0; i<NumObs; i++){
			myfile << BlueCov[i + k*(NumObs) ] << "  ";
		}
		myfile << "\n";
	}
	myfile.close();


  delete myBlue;
}


