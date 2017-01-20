#include <map>

#include "DataFormats/CTPPSReco/interface/TotemRPLocalTrack.h"

using namespace std;
using namespace edm;

//----------------------------------------------------------------------------------------------------

struct TrackData
{
	bool valid = false;
	double x = 0.;
	double y = 0.;

	void operator= (const TotemRPLocalTrack &ft)
	{
		valid = ft.isValid();
		x = ft.getX0();
		y = ft.getY0();
	}
};

//----------------------------------------------------------------------------------------------------

struct TrackDataCollection : public map<unsigned int, TrackData>
{
};
