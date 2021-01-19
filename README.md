# SolarAzEl
C++ code to estimate Solar Azimuth and Elevation given GPS position and time.

```cpp
double lat = 52.975;
double lon = -6.0494;
double altitude = 0;
 
double Az = 0.0;
double El = 0.0;
SolarAzEl(time(NULL), lat, lon, 0, &amp;Az, &amp;El);
 
printf("Azimuth: %f\n", Az);
printf("Elevation: %f\n", El);
```
