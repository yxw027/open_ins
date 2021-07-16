#include <math.h>
#include "ins_interface_API.h"
#include "lcsystem.h"
#include "lcgnssupdate.h"
#include "lcinsupdate.h"
#include "loosecoupleset.h"
#include "ginsdebug.h"
#include "insoutmsg.h"

GnssData gnssData1 = {0};

int8_t nhcflag = 0;
int8_t GNSSOBSFLAG = 0;
double lastgnsstime = -1;
int8_t OBS_outpus = 0;
extern int KFStatus;// = 0; // 0¡êoempty 1¡êoUpdata 2¡êofeedback

#ifndef SECONDSOFWEEK
#define SECONDSOFWEEK (604800)
#endif // !SECONDSOFWEEK

#ifndef SCALEFACTOR_ODO 
#define SCALEFACTOR_ODO (0.97 * 10 * 0.640 * 3.1415926/2000)
#endif


int8_t insinitsystemfromcfg(const LCSetting lcSetting)
{
	initsystemfromcfg(lcSetting);
	return 1;
}

int8_t SetStartGNSSWeek(int32_t week)
{
	if (week > 1024 && week < 3072)
	{
		if (gps_start_week == -1)
		{
			gps_start_week = week;
		}
	}
	return 1;
}
int8_t IsStartGNSSWeek()
{
	if (gps_start_week != -1)
	{
		return 1;
	}
	else
	{
		return -1;
	}
}

int8_t GetStartGNSSWeek(int32_t* week)
{
	*week = gps_start_week;
	if (gps_start_week != -1)
	{
		return 1;
	}
	else
	{
		return -1;
	}
	
}
int8_t GetKFStatus()
{
	return KFStatus;
}
/*  ADD origin GNSS data and Updata process
	args:  GnssData *mGnssData              origin gnss data struct
	return:status ()

*/
int8_t INSADDGNSSDATA(const GnssData mGnssData)
{
	int8_t ret = -1;
	mGnssInsSystem.GNSSflag = 1;
	if (mGnssInsSystem.mImuData.timestamp + 0.2/mLCSetting.imuDataRate >= mGnssData.timestamp )
	{
		KFStatus = 1;
		ret = ADDGNSSDATA(mGnssData);
	
	    if (isEFKFinished)
	    {
	      	KFStatus = 2;
	    }
	    else
	    {
		    KFStatus = 0;
	     }
	    OBS_outpus = 1;
    }
	else
	{
		ret = 0;
		memcpy(&gnssData1,&mGnssData,sizeof(GnssData));
	}
    return ret;
}
int8_t INSAddIMUData(const ImuData mImuData)
{
	int8_t ret = -1;
	if (1 == mGnssInsSystem.GNSSflag
		&& 0 == KFStatus
		&& mGnssInsSystem.mImuData.timestamp + 0.2 / mLCSetting.imuDataRate >= gnssData1.timestamp)
	{
		KFStatus = 1;
		ADDGNSSDATA(gnssData1);

		if (isEFKFinished)
		{
			KFStatus = 2;
		}
		else
		{
			KFStatus = 0;
		}
		OBS_outpus = 1;
	}
	ret = AddIMUData(mImuData);
#ifndef DEBUGMESSAGE
	    int flag = OBS_outpus;
		int32_t outtime_week = mImuData.week;
		double  outtime_second = mImuData.timestamp;
		uint32_t outtime_second_m = (uint32_t)(outtime_second * 1000);
		if (mImuData.timestamp >= SECONDSOFWEEK)
		{
			int32_t df_week = floor(mImuData.timestamp / SECONDSOFWEEK);
			outtime_second = mImuData.timestamp -   df_week * SECONDSOFWEEK;
			outtime_second_m = (uint32_t)(outtime_second * 1000);
			outtime_week = mImuData.week + SECONDSOFWEEK;
		}

		int32_t  outtime_week_predict = 0;
		double  outtime_second_predict = 0.0;
		uint32_t outtime_second_predict_m = 0;

		if (1 == flag)
		{
			outtime_second_predict = mImuData.timestamped;
			outtime_week_predict = mImuData.week;
			outtime_second_predict_m = (uint32_t)(outtime_second_predict * 1000);
			if (mImuData.timestamped >= SECONDSOFWEEK)
			{
				int32_t df_week = floor(mImuData.timestamped / SECONDSOFWEEK);
				outtime_second_predict = mImuData.timestamped - df_week * SECONDSOFWEEK;
				outtime_second_predict_m = (uint32_t)(outtime_second_predict * 1000);
				outtime_week_predict = mImuData.week + df_week;
			}
		}

		if (1) //OUTPUTGNSSMSG
		{
			if (1 == OBS_outpus)
			{
				int32_t outtime_gnss_week = mGnssInsSystem.mGnssData.week;
				double  outtime_gnss_second = mGnssInsSystem.mGnssData.timestamp;
				if (mImuData.timestamp >= SECONDSOFWEEK)
				{
					int32_t df_week = floor(mGnssInsSystem.mGnssData.timestamp / SECONDSOFWEEK);
					outtime_gnss_second = mGnssInsSystem.mGnssData.timestamp - df_week * SECONDSOFWEEK;
					outtime_gnss_week = mGnssInsSystem.mGnssData.week + df_week;
				}
				uint32_t itow = (uint32_t)(outtime_gnss_second * 1000);
				writeGnssRawMsg(&OBS_outpus, outtime_gnss_week, itow, &mGnssInsSystem.mGnssData);
			}
		}
		if (1) //OUTPUTRAWMSG
		{
				writeRawImuMsg(flag, outtime_week, outtime_second_m, outtime_week_predict, outtime_second_predict_m, &mImuData);
		}
		if (1) //OUTPUTINSPVAX
		{
			writeINSPVAXMsg(outtime_week, outtime_second_m, &mGnssInsSystem);
		}
		if (1) //OUTPUTINSPVA
		{
			writeINSPVAMsg(outtime_week, outtime_second_m, &mGnssInsSystem);
		}
		if (1) //OUTPUTGGA
		{
			if (mGnssInsSystem.mlc_STATUS == 4 && fmod(outtime_second +0.001, 0.2) < 0.005)
			{
				writeGGAMsg(outtime_week, outtime_second, &mGnssInsSystem);
			}
		}
#endif
	return ret;
}

int8_t INSAddODOData(const OdoData odo_data)
{
	int8_t ret = -1;
#ifndef DEBUGMESSAGE
	int32_t outtime_week = odo_data.week;
	double  outtime_second = odo_data.timestamp;
	uint32_t outtime_second_m = (uint32_t)(outtime_second * 1000);
	if (odo_data.timestamp >= SECONDSOFWEEK)
	{
		int32_t df_week = floor(odo_data.timestamp / SECONDSOFWEEK);
		outtime_second = odo_data.timestamp - df_week * SECONDSOFWEEK;
		outtime_second_m = (uint32_t)(outtime_second * 1000);
		outtime_week = odo_data.week + SECONDSOFWEEK;
	}
	writeOdoDataMsg(&odo_data);
#endif // DEBUGMESSAGE
	double curent_vehicle_speed = 0.0;

		switch (odo_data.mode)
		{
		case 0:
		{
			ret = 1;
			mGnssInsSystem.Odomode = 1;
			curent_vehicle_speed = odo_data.vehicle_speed;
		}break;
		case 1:
		{
			if (odo_data.timestamp - mGnssInsSystem.mOdoData.timestamp > 0 && mGnssInsSystem.mOdoData.timestamp > 0)
			{
				double vehicle_speed;
				if (odo_data.timestamp - mGnssInsSystem.mOdoData.timestamp > 0.99 && odo_data.timestamp - mGnssInsSystem.mOdoData.timestamp < 1.01)
				{
					vehicle_speed = SCALEFACTOR_ODO * (odo_data.wheel_tick - mGnssInsSystem.mOdoData.wheel_tick);
				}
				else
				{
					vehicle_speed = SCALEFACTOR_ODO * (odo_data.wheel_tick - mGnssInsSystem.mOdoData.wheel_tick) * 0.1 / (odo_data.timestamp - mGnssInsSystem.mOdoData.timestamp);
				}
				 vehicle_speed = mGnssInsSystem.mOdoData.fwd == 1 ? vehicle_speed : -vehicle_speed;
				 if (fabs(vehicle_speed) < 0.01)
				 {
					 curent_vehicle_speed = 0;
				 }
				 else
				 {
					curent_vehicle_speed = vehicle_speed;
					// curent_vehicle_speed = 0.10 * mGnssInsSystem.mOdoData.vehicle_speed + 0.90 * vehicle_speed;
			      }	
				 ret = 1;
			}
			else
			{
				ret = 0;
				mGnssInsSystem.IsUseOdo = 0;
			}
		}break;
		default:
			break;
		}

	
#ifdef DEBUGMESSAGE
tracerawodo(&odo_data);
#endif // DEBUGMESSAGE
	mGnssInsSystem.mOdoData.week = odo_data.week;
	mGnssInsSystem.mOdoData.timestamp = odo_data.timestamp;
	mGnssInsSystem.mOdoData.vehicle_speed = curent_vehicle_speed;
	mGnssInsSystem.mOdoData.wheel_tick = odo_data.wheel_tick;
	mGnssInsSystem.mOdoData.fwd = odo_data.fwd;
	return ret;
}

int8_t stop()
{
	int8_t ret = -1;
#ifdef DEBUGMESSAGE
	file_close();
#endif // INSOUPUTFILE
	ret = 1;
	return ret;
}
