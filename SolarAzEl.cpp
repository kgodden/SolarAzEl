
#include <time.h>
#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI (3.14159265358979323846264338327950288)
#endif /* M_PI */

// Programed by Darin C.Koblick 2 / 17 / 2009
//
//              Darin C.Koblick 4 / 16 / 2013 Vectorized for Speed
//                                         Allow for MATLAB Datevec input in
//                                         addition to a UTC string.
//                                         Cleaned up comments and code to
//                                         avoid warnings in MATLAB editor.
//
//				Kevin Godden 9/1/2020      Ported from Matlab to C++, tried to change as little as possible.
//                                         this is a non-vectorised port.
//
//--------------------------------------------------------------------------
//
// External Function Call Sequence :
//
// double lat = 52.975;
// double lon = -6.0494;
// double altitude = 0;
//
// double Az = 0.0;
// double El = 0.0;
// SolarAzEl(time(NULL), lat, lon, 0, &Az, &El);
// 
// printf("Azimuth: %f\n", Az);
// printf("Elevation: %f\n", El);
//
// Or to calculate Az & El for an arbitary UTC time:
//
//
// tm utc;
// tm_year is time since 1900
// utc.tm_year = y - 1900;
// Month is zero based, i.e. Jan is month 0
// utc.tm_mon = m - 1;
// utc.tm_mday = d;
// utc.tm_hour = 10;
// utc.tm_min = 16;
// utc.tm_sec = 00;
// utc.tm_isdst = 0;
// 
// Get UTC time_t val
// tim = timegm(&utc);	// or _mkgmtime() on windows
// 
// double altitude = 0;
// double Az = 0.0;
// double El = 0.0;
// 
// double lat = 52.975;
// double lon = -6.0494;
// 
// SolarAzEl(tim, lat, lon, 0, &Az, &El);
// 
// printf("Az: %f\n", Az);
// printf("El: %f\n", El);
// 
//
// Function Description :
//
// SolarAzEl will ingest a Universal Time, and specific site location on earth
// it will then output the solar Azimuth and Elevation angles relative to that
// site.
//
// Input Description :
//
// utc_time_point : time_t containing target time for sun position calculations.
//
// Lat : Site Latitude in degrees -90:90->S(-) N(+)
//
// Lon : Site Longitude in degrees -180:180 W(-) E(+)
//
// Alt : Altitude of the site above sea level(Km)
//
// Output Description :
//  Az    Azimuth location of the sun(deg)
//  El    Elevation location of the sun(deg)
//
//
// Source References :
// Solar Position obtained from :
// http ://stjarnhimlen.se/comp/tutorial.html#5
//

double julian_day(time_t utc_time_point);

void SolarAzEl(time_t utc_time_point, double Lat, double Lon, double Alt, double* Az, double* El) {
	double jd = julian_day(utc_time_point);

	double d = jd - 2451543.5;
	
	// Keplerian Elements for the Sun(geocentric)
	double w = 282.9404 + 4.70935e-5*d; // (longitude of perihelion degrees)
	// a = 1.000000; % (mean distance, a.u.)
	double e = 0.016709 - 1.151e-9*d; // (eccentricity)
	double M = fmod(356.0470 + 0.9856002585*d, 360.0); // (mean anomaly degrees)
		
	double L = w + M; // (Sun's mean longitude degrees)

	double oblecl = 23.4393 - 3.563e-7*d; // (Sun's obliquity of the ecliptic)

	// auxiliary angle
	double  E = M + (180 / M_PI)*e*sin(M*(M_PI / 180))*(1 + e*cos(M*(M_PI / 180)));

	// rectangular coordinates in the plane of the ecliptic(x axis toward perhilion)
	double x = cos(E*(M_PI / 180)) - e;
	double y = sin(E*(M_PI / 180))*sqrt(1 - pow(e, 2));

	// find the distance and true anomaly
	double r = sqrt(pow(x,2) + pow(y,2));
	double v = atan2(y, x)*(180 / M_PI);

	// find the longitude of the sun
	double lon = v + w;

	// compute the ecliptic rectangular coordinates
	double xeclip = r*cos(lon*(M_PI / 180));
	double yeclip = r*sin(lon*(M_PI / 180));
	double zeclip = 0.0;
	//rotate these coordinates to equitorial rectangular coordinates
	double xequat = xeclip;

	double yequat = yeclip*cos(oblecl*(M_PI / 180)) + zeclip * sin(oblecl*(M_PI / 180));

	double zequat = yeclip*sin(23.4406*(M_PI / 180)) + zeclip * cos(oblecl*(M_PI / 180));
	// convert equatorial rectangular coordinates to RA and Decl:
	r = sqrt(pow(xequat, 2) + pow(yequat, 2) + pow(zequat, 2)) - (Alt / 149598000); //roll up the altitude correction
	double RA = atan2(yequat, xequat)*(180 / M_PI);

	double delta = asin(zequat / r)*(180 / M_PI);
	
	// Following the RA DEC to Az Alt conversion sequence explained here :
	// http ://www.stargazing.net/kepler/altaz.html
	//	Find the J2000 value
	//	J2000 = jd - 2451545.0;
	//hourvec = datevec(UTC);
	//UTH = hourvec(:, 4) + hourvec(:, 5) / 60 + hourvec(:, 6) / 3600;

	// Get UTC representation of time / C++ Specific
	tm *ptm;
	ptm = gmtime(&utc_time_point);
	double UTH = (double)ptm->tm_hour + (double)ptm->tm_min / 60 + (double)ptm->tm_sec / 3600;

	// Calculate local siderial time
	double GMST0 = fmod(L + 180, 360.0) / 15;

	double SIDTIME = GMST0 + UTH + Lon / 15;
	
	// Replace RA with hour angle HA
	double HA = (SIDTIME*15 - RA);

	// convert to rectangular coordinate system
	x = cos(HA*(M_PI / 180))*cos(delta*(M_PI / 180));

	y = sin(HA*(M_PI / 180))*cos(delta*(M_PI / 180));
	double z = sin(delta*(M_PI / 180));

	// rotate this along an axis going east - west.
	double xhor = x*cos((90 - Lat)*(M_PI / 180)) - z*sin((90 - Lat)*(M_PI / 180));

	double yhor = y;
	double zhor = x*sin((90 - Lat)*(M_PI / 180)) + z*cos((90 - Lat)*(M_PI / 180));
	
	// Find the h and AZ
	*Az = atan2(yhor, xhor)*(180 / M_PI) + 180;
	*El = asin(zhor)*(180 / M_PI);
}

double julian_day(time_t utc_time_point) {

	// Extract UTC Time
	struct tm* tm = gmtime(&utc_time_point);

	double year = tm->tm_year + 1900;
	double month = tm->tm_mon + 1;
	double day = tm->tm_mday;
	double hour = tm->tm_hour;
	double min = tm->tm_min;
	double sec = tm->tm_sec;

	if (month <= 2) {
		year -= 1;
		month += 12;
	}

	double jd = floor(365.25*(year + 4716.0)) + floor(30.6001*(month + 1.0)) + 2.0 -
		floor(year / 100.0) + floor(floor(year / 100.0) / 4.0) + day - 1524.5 +
		(hour + min / 60 + sec / 3600) / 24;

	return jd;
}


/*
Copyright(c) 2010, Darin Koblick
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met :

*Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED.IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/
