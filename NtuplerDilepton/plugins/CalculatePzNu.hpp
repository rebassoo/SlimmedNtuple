/**
 * \file CalculatePzNu.hpp
 * \author Orso Iorio
 * 
 * A slightly modified version of the routine to reconstruct the z-component of neutrino momentum
 * in single-top events.
 */

#pragma once

#include <TLorentzVector.h>
#include <TVector2.h>

#include <iostream>
#include <vector>
#include <complex>


//..................................................................................................
// Equation solver by Orso

template <class T>
std::vector< T > const EquationSolve(const T & a, const T & b,const T & c,const T & d){


	std::vector<T> result;


	std::complex<T> x1;
	std::complex<T> x2;
	std::complex<T> x3;

	if (a != 0) {

		T q = (3*a*c-b*b)/(9*a*a);
		T r = (9*a*b*c - 27*a*a*d - 2*b*b*b)/(54*a*a*a);
		T Delta = q*q*q + r*r;

		std::complex<T> s;
		std::complex<T> t;

		T rho=0;
		T theta=0;

		if( Delta<=0){
			rho = sqrt(-(q*q*q));

			theta = acos(r/rho);

			s = std::polar<T>(sqrt(-q),theta/3.0); 
			t = std::polar<T>(sqrt(-q),-theta/3.0); 
		}

		if(Delta>0){ 
			s = std::complex<T>(cbrt(r+sqrt(Delta)),0);
			t = std::complex<T>(cbrt(r-sqrt(Delta)),0);
		}

		std::complex<T> i(0,1.0); 


		x1 = s+t+std::complex<T>(-b/(3.0*a),0);
		x2 = (s+t)*std::complex<T>(-0.5,0)-std::complex<T>(b/(3.0*a),0)+(s-t)*i*std::complex<T>(sqrt(3)/2.0,0);
		x3 = (s+t)*std::complex<T>(-0.5,0)-std::complex<T>(b/(3.0*a),0)-(s-t)*i*std::complex<T>(sqrt(3)/2.0,0);

		if(fabs(x1.imag())<0.0001)result.push_back(x1.real());
		if(fabs(x2.imag())<0.0001)result.push_back(x2.real());
		if(fabs(x3.imag())<0.0001)result.push_back(x3.real());

		return result;
	}
	else{return result;}


	return result;
}


//..................................................................................................
/*
 * Modification of code by Orso.
 * WARNING! Use only pz component of the returned 4-momenta as the function adjusts the transverse
 * component.
 */


TLorentzVector Nu4Momentum(TLorentzVector const &Lepton, float const & metPt, float const & metPhi){

	// Set the flags from the constructor
	bool useNegativeDeltaSolutions_ = true;
	bool usePositiveDeltaSolutions_ = true;
	bool usePzMinusSolutions_ = false;
	bool usePzPlusSolutions_ = false;
	bool usePzAbsValMinimumSolutions_ = true;
	bool useMetForNegativeSolutions_ = false;
	bool usePxMinusSolutions_ = true;
	bool usePxPlusSolutions_ = true;

	double const mW = 80.38;

	//  double Wmt = sqrt(pow(Lepton.et()+MET.pt(),2) - pow(Lepton.px()+MET.px(),2) - pow(Lepton.py()+MET.py(),2) );

	float const metpx = metPt * cos(metPhi);
	float const metpy = metPt * sin(metPhi);

	double MisET2 = (metpx*metpx + metpy*metpy);
	double mu = (mW*mW)/2 + metpx*Lepton.Px() + metpy*Lepton.Py();
	double a  = (mu*Lepton.Pz())/(Lepton.Energy()*Lepton.Energy() - Lepton.Pz()*Lepton.Pz());
	double a2 = TMath::Power(a,2);
	double b  = (TMath::Power(Lepton.Energy(),2.)*(MisET2) - TMath::Power(mu,2.))/(TMath::Power(Lepton.Energy(),2) - TMath::Power(Lepton.Pz(),2));
	double pz1(0),pz2(0),pznu(0);
	//int nNuSol(0);

	TLorentzVector p4nu_rec;
	TLorentzVector p4W_rec;
	TLorentzVector p4b_rec;
	TLorentzVector p4Top_rec;
	TLorentzVector p4lep_rec;    

	p4lep_rec.SetPxPyPzE(Lepton.Px(),Lepton.Py(),Lepton.Pz(),Lepton.Energy());

	TLorentzVector p40_rec(0,0,0,0);

	if(a2-b > 0 ){
		if(!usePositiveDeltaSolutions_)
		{
			return p40_rec;
		}
		double root = sqrt(a2-b);
		pz1 = a + root;
		pz2 = a - root;
		//nNuSol = 2;     




		if(usePzPlusSolutions_)pznu = pz1;    
		if(usePzMinusSolutions_)pznu = pz2;
		if(usePzAbsValMinimumSolutions_){
			pznu = pz1;
			if(fabs(pz1)>fabs(pz2)) pznu = pz2;
		}


		double Enu = sqrt(MisET2 + pznu*pznu);

		p4nu_rec.SetPxPyPzE(metpx, metpy, pznu, Enu);

		return p4nu_rec;

	}
	else{

		if(!useNegativeDeltaSolutions_){
			return p40_rec;
		}
		//    double xprime = sqrt(mW;


		double ptlep = Lepton.Pt(),pxlep=Lepton.Px(),pylep=Lepton.Py();

		double EquationA = 1;
		double EquationB = -3*pylep*mW/(ptlep);
		double EquationC = mW*mW*(2*pylep*pylep)/(ptlep*ptlep)+mW*mW-4*pxlep*pxlep*pxlep*metpx/(ptlep*ptlep)-4*pxlep*pxlep*pylep*metpy/(ptlep*ptlep);
		double EquationD = 4*pxlep*pxlep*mW*metpy/(ptlep)-pylep*mW*mW*mW/ptlep;

		std::vector<long double> solutions = EquationSolve<long double>((long double)EquationA,(long double)EquationB,(long double)EquationC,(long double)EquationD);

		std::vector<long double> solutions2 = EquationSolve<long double>((long double)EquationA,-(long double)EquationB,(long double)EquationC,-(long double)EquationD);


		double deltaMin = 14000*14000;
		double zeroValue = -mW*mW/(4*pxlep); 
		double minPx=0;
		double minPy=0;

		//    std::cout<<"a "<<EquationA << " b " << EquationB  <<" c "<< EquationC <<" d "<< EquationD << std::endl; 

		if(usePxMinusSolutions_){
			for( int i =0; i< (int)solutions.size();++i){
				if(solutions[i]<0 ) continue;
				double p_x = (solutions[i]*solutions[i]-mW*mW)/(4*pxlep); 
				double p_y = ( mW*mW*pylep + 2*pxlep*pylep*p_x -mW*ptlep*solutions[i])/(2*pxlep*pxlep);
				double Delta2 = (p_x-metpx)*(p_x-metpx)+(p_y-metpy)*(p_y-metpy); 

				//      std::cout<<"intermediate solution1 met x "<<metpx << " min px " << p_x  <<" met y "<<metpy <<" min py "<< p_y << std::endl; 

				if(Delta2< deltaMin && Delta2 > 0){deltaMin = Delta2;
					minPx=p_x;
					minPy=p_y;}
				//     std::cout<<"solution1 met x "<<metpx << " min px " << minPx  <<" met y "<<metpy <<" min py "<< minPy << std::endl; 
			}
		} 

		if(usePxPlusSolutions_){
			for( int i =0; i< (int)solutions2.size();++i){
				if(solutions2[i]<0 ) continue;
				double p_x = (solutions2[i]*solutions2[i]-mW*mW)/(4*pxlep); 
				double p_y = ( mW*mW*pylep + 2*pxlep*pylep*p_x +mW*ptlep*solutions2[i])/(2*pxlep*pxlep);
				double Delta2 = (p_x-metpx)*(p_x-metpx)+(p_y-metpy)*(p_y-metpy); 
				//  std::cout<<"intermediate solution2 met x "<<metpx << " min px " << minPx  <<" met y "<<metpy <<" min py "<< minPy << std::endl; 
				if(Delta2< deltaMin && Delta2 > 0){deltaMin = Delta2;
					minPx=p_x;
					minPy=p_y;
				}
				//	std::cout<<"solution2 met x "<<metpx << " min px " << minPx  <<" met y "<<metpy <<" min py "<< minPy << std::endl; 
			}
		}

		double pyZeroValue= ( mW*mW*pxlep + 2*pxlep*pylep*zeroValue);
		double delta2ZeroValue= (zeroValue-metpx)*(zeroValue-metpx) + (pyZeroValue-metpy)*(pyZeroValue-metpy);

		if(deltaMin==14000*14000)return TLorentzVector();    
		//    else std::cout << " test " << std::endl;

		if(delta2ZeroValue < deltaMin){
			deltaMin = delta2ZeroValue;
			minPx=zeroValue;
			minPy=pyZeroValue;}

		//    std::cout<<" MtW2 from min py and min px "<< sqrt((minPy*minPy+minPx*minPx))*ptlep*2 -2*(pxlep*minPx + pylep*minPy)  <<std::endl;
		///    ////Y part   

		double mu_Minimum = (mW*mW)/2 + minPx*pxlep + minPy*pylep;
		double a_Minimum  = (mu_Minimum*Lepton.Pz())/(Lepton.Energy()*Lepton.Energy() - Lepton.Pz()*Lepton.Pz());
		pznu = a_Minimum;

		if(!useMetForNegativeSolutions_){
			double Enu = sqrt(minPx*minPx+minPy*minPy + pznu*pznu);
			p4nu_rec.SetPxPyPzE(minPx, minPy, pznu , Enu);
		}
		else{
			pznu = a;
			double Enu = sqrt(metpx*metpx+metpy*metpy + pznu*pznu);
			p4nu_rec.SetPxPyPzE(metpx, metpy, pznu , Enu);
		}
		return p4nu_rec;
	}
}
