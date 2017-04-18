#ifndef _track_lite_h_
#define _track_lite_h_

#include <map>

#include "DataFormats/CTPPSReco/interface/TotemRPLocalTrack.h"

//----------------------------------------------------------------------------------------------------

struct TrackData
{
	bool valid = false;

	// track position, in mm
	double x = 0.;
	double y = 0.;

	// track angle, in rad
	double th_x = 0.;
	double th_y = 0.;

	void operator= (const TotemRPLocalTrack &ft)
	{
		valid = ft.isValid();

		x = ft.getX0();
		y = ft.getY0();

		th_x = ft.getTx();
		th_y = ft.getTy();
	}
};

//----------------------------------------------------------------------------------------------------

struct TrackDataCollection : public std::map<unsigned int, TrackData>
{
};

#endif
