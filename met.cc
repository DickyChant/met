#include <cmath>
#include <iostream>
#include <vector>

#include "TH1.h"
#include "TF1.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TRandom.h"

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>

#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_hermite.h>


using namespace Eigen;
using namespace std;

//constexpr int A=40;
constexpr int A=16;

constexpr double range = 0.17;

constexpr int tot = 1000;

constexpr double coeff_herm[3] = {0.1,0.1,0.1};
const double Norm = 0.090103;


struct particle{
	double x;
	double y;
	double z;
};

void init_particle(TH3D* PDF, particle* Ns, int size_A=tot){
	for (int i = 0 ; i < size_A ; i++){
		PDF->GetRandom3((Ns[i].x),(Ns[i].y),(Ns[i].z));
	}
}

void init_particle(TH1D* PDF, particle* Ns, int size_A=tot){
	for (int i = 0 ; i < size_A ; i++){
		Ns[i].x=PDF->GetRandom();
		Ns[i].y=PDF->GetRandom();
		Ns[i].y=PDF->GetRandom();
	}
}


double wave_func_sing(particle N){
	double Hx=0,Hy=0,Hz=0;
	for (int i = 0 ;i<3;i++){
		Hx+=gsl_sf_hermite_func(i,N.x)*coeff_herm[i];
		Hy+=gsl_sf_hermite_func(i,N.y)*coeff_herm[i];
		Hz+=gsl_sf_hermite_func(i,N.z)*coeff_herm[i];
	}
	return Hx*Hy*Hz/pow(Norm,3)*exp((-N.x*N.x-N.y*N.y-N.z*N.z)/2);
}

double corr_func(particle N1, particle N2){
	double L2 = (N1.x-N2.x)*(N1.x-N2.x)+(N1.y-N2.y)*(N1.y-N2.y)+(N1.z-N2.z)*(N1.z-N2.z);
	if (L2>(0.1*range*range)){
		return 1;
	}
	else{
		return 0;
	}
}

double total_prob(particle* Ns, int* N_sel, int size_A=A){
	double total = 1.;
	for ( int i = 0 ; i < size_A ; i++){
		total*=wave_func_sing(Ns[N_sel[i]]);
	}
	for ( int i = 0; i < size_A ; i++){
		for( int j = i ; j < size_A ; j++){
			total*=corr_func(Ns[N_sel[i]],Ns[N_sel[i]]);
		}
	}
	return total*total;
}

bool all_searched = false;

void recurrent(int* N_sel, int size=A, int pos = A -1 ,int total_number = tot ){
	if(N_sel[pos] == total_number - size + pos){
		if(pos == 0){
			for(int i = 0 ; i < size ; i++){
				N_sel[i] = 0;
			}
			all_searched = true;
		}
		else{
			recurrent(N_sel,size,pos-1,total_number);
		}
	}
	else{
		N_sel[pos] += 1;
		for (int i = 1 ; i < size - pos ; i++){
			N_sel[i+pos] = N_sel[pos]+i;
		}
	}
}


TH1D* generate_pdf(){
	TFile*f = new TFile("test.root","recreate");
	f->cd();
	TH1D* test_pdf = new TH1D("test_pdf","test_pdf",10000,0,range);
	for ( int i = 0 ; i < 10000 ; i++){
		test_pdf->SetBinContent(i+1,cos(0.7+i/10000.0)*10000);
	}
	test_pdf->Write();
	TCanvas c("test","test",800,600);
	test_pdf->Draw();
	c.SaveAs("test.pdf");
	f->Close();
	delete f;

	return test_pdf;
}


int main(){
	

	// initialize particles
	TH1D* test_pdf = generate_pdf();
	particle * Ns = new particle[tot];
	init_particle(test_pdf,Ns);
	
	// metropolis search
	// met init
	int N_sel[A]={0};//record selected
	int N_tmp[A]={0};//record current

	for (int i = 0 ; i < A ;i++){
		N_sel[i]=i;
		N_tmp[i]=N_sel[i];
	}

	all_searched = false;

	double p1=0,p2=0;
	TRandom* rdn = new TRandom();

	do{
		recurrent(N_sel);
		p1 = total_prob(Ns,N_sel)/total_prob(Ns,N_tmp);
		p2 = rdn->Rndm();
		if(p1>p2){
			for (int i = 0 ; i < A ;i++){
				N_tmp[i]=N_sel[i];
			}
		}
	}while(!all_searched);

	for(int i = 0 ; i < A ;i++){
		printf("(%f,%f,%f)",Ns[N_sel[i]].x,Ns[N_sel[i]].y,Ns[N_sel[i]].z);
	}
	return 0;
}