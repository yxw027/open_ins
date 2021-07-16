#ifndef INS_INTERFACE_API
#define INS_INTERFACE_API
#include <stdint.h>
#include "datatype.h"
#include "loosecoupleset.h"
#ifdef __cplusplus
extern "C" {
#endif

#ifndef WIN32
#define ARM_MCU
#endif

#ifdef ARM_MCU
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wunused-const-variable="
#endif

int8_t insinitsystemfromcfg(const LCSetting lcSetting);
int8_t SetStartGNSSWeek(int32_t week);
int8_t IsStartGNSSWeek();
int8_t GetStartGNSSWeek(int32_t* week);
int8_t INSADDGNSSDATA(const GnssData mGnssData);
int8_t INSAddIMUData(const ImuData mImuData);
int8_t INSAddODOData(const OdoData odo_data);

int8_t stop();
int8_t GetKFStatus();
extern int8_t nhcflag;
extern int8_t GNSSOBSFLAG;
extern double lastgnsstime;
extern int KFStatus;
extern LCSetting mLCSetting;
/*--------------------------------------------------------------------*/
#ifdef __cplusplus
}
#endif
#endif
