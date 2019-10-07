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
#include <gsl/gsl_sf_bessel.h>


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
//
constexpr int N_1 = 1E8;
constexpr int N_2 = 1E5;
constexpr int N_3 = 16;

constexpr double n0 = 0.17;

constexpr double cube_range = 838.;

class nuclei{
	private:
		double x;
		double y;
		double z;
		bool isproton;//true for proton, false for neutron
	public:
		nuclei(){
			this->x=0;
			this->y=0;
			this->z=0;
			this->isproton=false;
		}
		nuclei(double x0,double y0,double z0,bool isproton0){
			this->x=x0;
			this->y=y0;
			this->z=z0;
			this->isproton=isproton0;
		}
		~nuclei(){};
		
		double X(){
			return this->x;
		}
		double Y(){
			return this->y;
		}
		double Z(){
			return this->z;
		}
		double r2(){
			return this->x*this->x+this->y*this->y+this->z*this->z; 
		}
		double r(){
			return sqrt(this->r2());
		}
		double set_xyz(double x0,double y0,double z0){
			this->x=x0;
			this->y=y0;
			this->z=z0;
			return this->r2();
		}
		bool set_type(bool isproton0){
			this->isproton=isproton0;
			return this->isproton;
		}
		bool Type(){
			return this->isproton;
		}
};

double rel_dis2(nuclei n1,nuclei n2){
	double dx = n1.X()-n2.X();
	double dy = n1.Y()-n2.Y();
	double dz = n1.Z()-n2.Z();
	return dx*dx+dy*dy+dz*dz;
}

double rel_dis(nuclei n1,nuclei n2){
	return sqrt(rel_dis2(n1,n2));
}

void init_nucl(nuclei* ns, TRandom* Rdm ,int n_1 = N-1){
	for (int i = 0 ; i <n_1 ; i++){
		ns[i].set_xyz(cube_range*Rdm->Rndm(),cube_range*Rdm->Rndm(),cube_range*Rdm->Rndm());
	}
}

constexpr double kappa = 0.05;
constexpr double alpha = 0.05;

double gnn(nuclei n1,nuclei n2){
	double dis = rel_dis(n1,n2);
	if (dis < 0.9 ){
		return 1;
	}
	else{
		return (kappa-1+exp((0.9-alpha)/dis))/(kappa);
	}
}

constexpr double C = 0.05;

double contact(nuclei n1,nuclei n2){
	double dis = rel_dis(n1,n2);
	return dis;
}

double rho_0(nuclei n1,nuclei n2){
	double dis = rel_dis(n1,n2);
	return dis;
}

constexpr double kf = 300;

double exch(nuclei n1,nuclei n2, int Z = A/2){
	double dis = rel_dis(n1,n2);
	return Z/2/(Z-1)*pow((3*gsl_sf_bessel_j0(kf*dis)/kf/dis),2);
}

double phi_zero(nuclei n1,nuclei n2){
	double dis= rel_dis(n1,n2);
	return dis;
}

double corr(nuclei n1,nuclei n2){
	
	return 0;
}

double full_corr(nuclei* ns,int size =N_2){
	double total = 1.;
	for (int i = 0 ; i < size ; i++){
		for(int j = i;j<size;j++){
			total *= corr(ns[i],ns[j]);
		}
	}
	return total;
}

constexpr double para_ppnn[5] = {3.17,0.995,1.81,5.90,-9.87};
constexpr double para_pn[5] = {1.08,0.985,-0.432,-3.30,2.01};

double F_oH(nuclei n1,nuclei n2){
	const double* para;
	if (n1.Type() xor n2.Type()){
		para = para_pn;
	}
	else{
		para = para_ppnn;
	}
	double r = rel_dis(n1,n2);
	double dx = n1.X() - n2.X();
	double dy = n1.Y() - n2.Y();
	double dz = n1.Z() - n2.Z();

	double F = 1. - exp(-para[0]*r*r)*(para[1]+r*(dx*para[2]+dy*para[3]+dz*para[4]));
	return F;
}

constexpr double fermi_mom = 1.0688572E-4;

double exch_factor(nuclei n1,nuclei n2){
	double r = rel_dis(n1,n2);
	return 1.-pow(3*gsl_sf_bessel_j0(r*fermi_mom)/r/fermi_mom,2)/2.;
}

double F_classic(nuclei n1,nuclei n2){
	double F = F_oH(n1,n2);
	double pauli_factor = exch_factor(n1,n2);
	return F * pauli_factor;
}

double total_F_fermi(nuclei* ns, int* N_sel,int size = N_2){
	double total = 1.;
	for (int i = 0 ; i < size ; i++){
		for(int j = 0 ; j < size ; j++){
			total*= F_oH(ns[N_sel[i]],ns[N_sel[j]]);
		}
	}
	return total;
}

double total_F_classic(nuclei* ns, int* N_sel,int size = N_2){
	double total = 1.;
	for (int i = 0 ; i < size ; i++){
		for(int j = 0 ; j < size ; j++){
			total*= F_classic(ns[N_sel[i]],ns[N_sel[j]]);
		}
	}
	return total;
}

void met_search(nuclei* ns, int* N_sel ,int size = N_2){
	
}

int main(){
	TRandom* rdm = new TRandom();
	nuclei* Ns = new nuclei[N_1];
	init_nucl(Ns,rdm);

	

	return 0;
}