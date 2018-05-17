#ifndef _proton_reconstruction_h_
#define _proton_reconstruction_h_

#include "TFile.h"
#include "TSpline.h"

#include <map>
#include <string>

#include "track_lite.h"

//----------------------------------------------------------------------------------------------------

struct ProtonData
{
	bool valid = false;

	double vtx_x = 0., vtx_y = 0.;	// m

	double th_x = 0., th_y = 0.;	// rad

	double xi = 0.;					// 1, positive when energy loss
	double xi_unc = 0.;
};

//----------------------------------------------------------------------------------------------------

TSpline3* PrepareOneFunction(const std::string &fn, const std::string &on)
{
	// load input graph
	TFile *f_in = TFile::Open(fn.c_str());
	if (f_in == NULL)
	{
		printf("ERROR: can't open file '%s'.\n", fn.c_str());
		exit(1);
	}

	TGraph *g_in = (TGraph *) f_in->Get(on.c_str());
	if (g_in == NULL)
	{
		printf("ERROR: can't load object '%s'.\n", on.c_str());
		exit(2);
	}

	// convert to positive xi
	TGraph *g_out = new TGraph();
	for (int i = 0; i < g_in->GetN(); i++)
	{
		double x, y;
		g_in->GetPoint(i, x, y);
		g_out->SetPoint(i, x, -y);
	}

	// build spline
	TSpline3 *s_out = new TSpline3("", g_out->GetX(), g_out->GetY(), g_out->GetN());

	// clean up
	delete f_in;

	return s_out;
}

//----------------------------------------------------------------------------------------------------

std::map<unsigned int, TSpline3 *> m_s_x_to_xi;

void InitReconstruction()
{
  

	m_s_x_to_xi[2] = PrepareOneFunction("xi_as_a_function_of_x_graph_b2.root", "XRPH_C6L5_B2");
	m_s_x_to_xi[3] = PrepareOneFunction("xi_as_a_function_of_x_graph_b2.root", "XRPH_D6L5_B2");

	m_s_x_to_xi[102] = PrepareOneFunction("xi_as_a_function_of_x_graph_b1.root", "XRPH_C6R5_B1");
	m_s_x_to_xi[103] = PrepareOneFunction("xi_as_a_function_of_x_graph_b1.root", "XRPH_D6R5_B1");

	m_s_x_to_xi[1980760064] = PrepareOneFunction("xi_as_a_function_of_x_graph_b2.root", "XRPH_C6L5_B2");
	m_s_x_to_xi[1981284352] = PrepareOneFunction("xi_as_a_function_of_x_graph_b2.root", "XRPH_D6L5_B2");

	m_s_x_to_xi[1997537280] = PrepareOneFunction("xi_as_a_function_of_x_graph_b1.root", "XRPH_C6R5_B1");
	m_s_x_to_xi[1998061568] = PrepareOneFunction("xi_as_a_function_of_x_graph_b1.root", "XRPH_D6R5_B1");


}

//----------------------------------------------------------------------------------------------------

void ReconstructProtonFromOneRP(const unsigned int rpId, const TrackData &track, ProtonData &result)
{
  //cout<<"I get in reconstruct proton from one RP"<<endl;
	// default result
	result.valid = false;

	// get x-to-xi spline
	const auto s_it = m_s_x_to_xi.find(rpId);
	//cout<<"I get here 0:"<<endl;
	//cout<<"rpId: "<<rpId<<endl;
	if (s_it == m_s_x_to_xi.end())
		return;
	//cout<<"I get here 1:"<<endl;
	const TSpline3 *s_x_to_xi = s_it->second;
	
	//cout<<"Track x: "<<track.x<<endl;

	// determine xi
	result.xi = s_x_to_xi->Eval(track.x*1E-3);		// spline expects x in m

	// determine uncertainty of xi
	const double si_x_alignment = 150E-6;			// in m, alignment uncertainty
	const double si_x_neglected_angle = 150E-6;		// in m, to (approximately) account for the neglected angular term in proton transport
	const double si_rel_D = 0.055;					// 1, relative uncertainty of dispersion

	const double si_x = sqrt( si_x_alignment*si_x_alignment + si_x_neglected_angle*si_x_neglected_angle );

	const double si_xi_from_x = s_x_to_xi->Eval(track.x*1E-3 + si_x) - result.xi;
	const double si_xi_from_D_x = si_rel_D * result.xi;

	result.xi_unc = sqrt( si_xi_from_x*si_xi_from_x + si_xi_from_D_x*si_xi_from_D_x );

	// reconstruction completed
	result.valid = true;
}

#endif
