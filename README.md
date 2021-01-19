# SolarAzEl
C++ code to estimate Solar Azimuth and Elevation given GPS position and time without the need to access to any online/web services.

This can be handy if you're imaging outside and want to gauge the effect of chamging sunight based on i position and time of day!

This code was written in Mathlab by Darin C.Koblick, it is a reasonably straight port to C++, vectorisation is not supported.

The repo consists of 1 .cpp file and 1 .h file, they define a function called SolarAzEl().

Example of using SolarAzEl().

```cpp
double lat = 52.975;
double lon = -6.0494;
double altitude = 0;
 
double Az = 0.0;
double El = 0.0;
SolarAzEl(time(NULL), lat, lon, 0, &Az, &El);
 
printf("Azimuth: %f\n", Az);
printf("Elevation: %f\n", El);
```
