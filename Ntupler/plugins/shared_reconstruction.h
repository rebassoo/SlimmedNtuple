#include "TFile.h"
#include "TSpline.h"

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

map<unsigned int, TSpline3 *> m_s_x_to_xi;

void InitReconstruction(const string &fin)
{
  //TFile *f_in = TFile::Open("/afs/cern.ch/work/j/jkaspar/analyses/ctpps/optics/optics.root");
  //TFile *f_in = TFile::Open("/home/users/rebassoo/work/2016_11_14_FinnNtupler/CMSSW_8_1_0_pre8/src/SlimmedNtuple/Ntupler/python/optics.root");
  TFile *f_in = TFile::Open(fin.c_str(),"r");
	if (!f_in)
	{
		printf("ERROR in InitReconstruction > Can't open input file.\n");
		exit(1);
	}

	m_s_x_to_xi[3]   = (TSpline3 *) f_in->Get("s_x_to_xi_L_1_F");
	m_s_x_to_xi[2]   = (TSpline3 *) f_in->Get("s_x_to_xi_L_1_N");
	m_s_x_to_xi[102] = (TSpline3 *) f_in->Get("s_x_to_xi_R_1_N");
	m_s_x_to_xi[103] = (TSpline3 *) f_in->Get("s_x_to_xi_R_1_F");
}

//----------------------------------------------------------------------------------------------------

void ReconstructFromOneRP(const map<unsigned int, TrackData> &tracks, unsigned int id, ProtonData &result)
{
	result.valid = false;

	auto it = tracks.find(id);
	if (it == tracks.end())
		return;
	
	const auto &tr = it->second;

	if (!tr.valid)
		return;

	TSpline3 *s_x_to_xi = m_s_x_to_xi[id];
	if (s_x_to_xi == NULL)
		return;

	result.valid = true;

	result.xi = s_x_to_xi->Eval(tr.x*1E-3);	// expects x in m

	double de_x = 0.150E-3;	// m
	double de_rel_D = 0.15;	// 1

	double de_xi = s_x_to_xi->Eval(tr.x*1E-3 + de_x) - result.xi;

	result.xi_unc = sqrt( pow(de_xi, 2.) + pow(de_rel_D * result.xi, 2.) );
}

//----------------------------------------------------------------------------------------------------

ProtonData ReconstructProton(const map<unsigned int, TrackData> &tracks, bool leftArm)
{
	ProtonData result;
	result.valid = false;

	//double D_1_F, D_1_N; // in mm
	unsigned int id_1_F, id_1_N;

	if (leftArm)
	{
		//D_1_N = 9.26E-02;
		//D_1_F = 9.22E-02;

		id_1_N = 2;
		id_1_F = 3;
	} else {
		//D_1_N = 5.81E-02;
		//D_1_F = 5.16E-02;

		id_1_N = 102;
		id_1_F = 103;
	}

	// try far
	ReconstructFromOneRP(tracks, id_1_F, result);

	// if far not successfull, try near
	if (!result.valid)
		ReconstructFromOneRP(tracks, id_1_N, result);

	return result;
}
