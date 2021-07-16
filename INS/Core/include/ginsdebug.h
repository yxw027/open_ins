#ifndef _GINS_DEBUG_
#define _GINS_DEBUG_
#include <stdint.h>
#include <stdio.h>
#include "lcstruct.h"

#ifdef __cplusplus
extern "C" {
#endif
extern int32_t gps_start_week;
int8_t setfilename(const char* inputname,const int32_t fileOutputType, const int32_t profiletype, const int32_t isKmlOutput, const int8_t useOdo);
int8_t ouputfile(const GnssInsSystem mGnssInsSystem, const int32_t  insOutputDataRate, const int32_t kmlOutputDateRate, const int32_t imuDataRate);
int8_t tracemis(const ImuData* pImuData, const int8_t* p1, const double*p2);
int8_t tracemis1(const double *a);
int8_t file_close();
int ouputgnssfile(const GnssData mGnssData, const int32_t gnssOutputDataRate, const int32_t kmlOutputDateRate);
int8_t tracerawimu(const ImuData* pImuData);
int8_t tracerawodo(const OdoData* pOdoData);

int8_t tracenativegnss(const double pnt[15]);

int8_t tracenativeodo(const double odo[15]);

int8_t tracenativeimu(const ImuData* pImuData);

int8_t tracenativertk(const GnssData* pgnssData);

#ifdef __cplusplus
}
#endif


#endif // !_GINS_DEBUG_
