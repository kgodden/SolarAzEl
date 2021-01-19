#pragma once

#include <time.h>

void SolarAzEl(time_t utc_time_point, double Lat, double Lon, double Alt, double* Az, double* El);
