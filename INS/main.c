#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "configManger.h"
#include "util.h"
#include "datatype.h"
#include "ins_interface_API.h"

#include "include/ginsdebug.h"
#include "include/lcsystem.h"


extern LCSetting mLCSetting;
ImuData mImuData;
GnssData mGnssData;
OdoData mOdoData;

double LastIMUDataTime;
double LastGNSSDataTime = -1;
double LastGNSSDataSYSTime = -1;
double lastlat = 0;
double lastlon = 0;
double lastalt = 0;



#define DEBUGMESSAGE (1)


int8_t playLogFile(char *inputfilename)
{
	int8_t ret = -1;
	FILE *fdat = NULL;
	char buffer[150] = { 0 };
	char  *val[40];
	int32_t sentencesize;
	int8_t gnssPOSExistFlag = -1;
	int8_t gnssVelExistFlag = -1;

	double delytime = 0;

	fopen_s(&fdat, inputfilename, "r");
	if (!fdat)
	{
#ifdef DEBUGMESSAGE
		printf("Cant OPEN PROCFILE: %s\n", inputfilename);
#endif // DEBUGMESSAGE
		return ret;
	}

	while (!feof(fdat))
	{
		fgets(buffer, sizeof(buffer), fdat);
		if (strlen(buffer) <= 0) continue;	
		char *temp = strchr(buffer, '\n');
		if (temp != NULL) temp[0] = '\0';

		parse_fields(buffer, val,&sentencesize);
		if (strstr(val[0], "$GPODO") != NULL)
		{
			if (sentencesize == 7)
			{
				mOdoData.week = atoi(val[1]);
				mOdoData.timestamp = atof(val[2]);
				mOdoData.mode = atoi(val[3]);
				mOdoData.vehicle_speed = atof(val[4]);
				mOdoData.fwd = atoi(val[5]);
				mOdoData.wheel_tick = atoll(val[6]);
				if (mOdoData.week > 1024 && mOdoData.week < 3072)
				{
					if (IsStartGNSSWeek() != 1)
					{
						SetStartGNSSWeek(mOdoData.week);
					}
					INSAddODOData(mOdoData);
				}
				else
				{

				}
			}

		}
		if (strstr(val[0], "$GPIMU") != NULL)
		{
			if (sentencesize == 10)  //保证IMU数据完整性
			{
				memset(&mImuData, 0, sizeof(ImuData));
				mImuData.week = atoi(val[1]);
				mImuData.timestamp = atof(val[2]) + 0.0001;
				mImuData.timestamped = atof(val[3]) + 0.0001;
				mImuData.Accx = floor(atof(val[4])*10000000)* 9.780 / 9.806 /10000000.0;
				mImuData.Accy = floor(atof(val[5]) * 10000000)* 9.780 / 9.806 / 10000000.0;
				mImuData.Accz = floor(atof(val[6]) * 10000000)* 9.780 / 9.806 / 10000000.0;
				mImuData.Gyrox = atof(val[7]) * PI / 180;
				mImuData.Gyroy = atof(val[8]) * PI / 180;
				mImuData.Gyroz = atof(val[9]) * PI / 180;
				if (mImuData.week > 1024 && mImuData.week < 3072)
				{
					if (IsStartGNSSWeek() != 1)
					{
						SetStartGNSSWeek(mImuData.week);
					}
#ifdef DEBUGMESSAGE
					tracerawimu(&mImuData);
#endif
					INSAddIMUData(mImuData);
				}
				else
				{

				}
			}
			else
			{
#ifdef DEBUGMESSAGE
				printf("IMU Data: Format error ,gnss_week: %d gnss_time0fweek: %.4f\n", mImuData.week, mImuData.timestamped);
#endif // DEBUGMESSAGE
				continue;
			}
		}
		if (strstr(val[0], "$GPGNSS") != NULL)
		{
			if (sentencesize >= 10)
			{
				//if (fabs(atof(val[3])*PI / 180 - lastlat) < 0.000001
				//	&& fabs(atof(val[4])*PI / 180 - lastlon) < 0.000001
				//	&& fabs(atof(val[5]) - lastalt) < 0.000001
				//	&& val[3] != 0
				//	&& fabs(lastlat)> 0.1)
				//{
				//	continue;
				//}
				//else
				{
					mGnssData.week = atoi(val[1]);// atoi(val[1]);  mImuData.week
					mGnssData.timestamp = atof(val[2]) + 0.0001; //atof(val[2]); floor(mImuData.timestamp)
					mGnssData.timestampd = mImuData.timestamp;
					mGnssData.latitude = atof(val[3])*PI / 180;
					mGnssData.longitude = atof(val[4])*PI / 180;
					mGnssData.altitude = atof(val[5]);
					mGnssData.latitude_std = atof(val[6]);
					mGnssData.longitude_std = atof(val[7]);
					mGnssData.altitude_std = atof(val[8]);
					mGnssData.gpsFixType = atoi(val[9]);
					mGnssData.Mode = (int)atof(val[9]);
					if (mGnssData.Mode == 1)
					{
						mGnssData.latitude_std = 10;
						mGnssData.longitude_std = 10;
						mGnssData.altitude_std = 10;
					}
					if (fabs(mGnssData.latitude) > 0.1)
					{
						lastlat = mGnssData.latitude;
						lastlon = mGnssData.longitude;
						lastalt = mGnssData.altitude;
					}
					gnssPOSExistFlag = 1;
				}
			}
			else
			{
#ifdef DEBUGMESSAGE
				printf("GNSS Data: Format error  ,gnss_week: %d gnss_time0fweek: %.4f\n", mGnssData.week, mGnssData.timestamp);
#endif // DEBUGMESSAGE
				continue;
			}
		}
		if (strstr(val[0], "$GPVEL") != NULL)
		{
			if (sentencesize >= 6)
			{
				double timevel = atof(val[2]); //atof(val[2]); floor(mImuData.timestamp)
				if (fabs(timevel - mGnssData.timestamp) < 0.01)
				{
					mGnssData.heading = atof(val[4])*PI / 180;
					mGnssData.north_velocity = atof(val[3])*cos(mGnssData.heading);
					mGnssData.east_velocity = atof(val[3])*sin(mGnssData.heading);
					mGnssData.up_velocity = -atof(val[5]);
					if (gnssPOSExistFlag == 1)
						gnssVelExistFlag = 1;
				}
				else
				{
#ifdef DEBUGMESSAGE
					printf("GNSS VEL Data: Format error  ,gnss_week: %d gnss_time0fweek: %.4f\n", mGnssData.week, mGnssData.timestamp);
#endif // DEBUGMESSAGE
					continue;
				}
			}
			else
			{
#ifdef DEBUGMESSAGE
				printf("GNSS VEL Data: Format error  ,gnss_week: %d gnss_time0fweek: %.4f\n", mGnssData.week, mGnssData.timestamp);
#endif // DEBUGMESSAGE
				// "Gnss Data error,line" 
			}

				
			
		}
		if (mLCSetting.gnssDataType == 1)
		{
			if (gnssPOSExistFlag == 1 && fabs(mImuData.timestamp - floor(mImuData.timestamp)) > 0.0001 )
			{
				if (mLCSetting.gnssSensorType == 0)
				{
					if (mGnssData.Mode == 4)
					{
						if (mGnssData.latitude_std < 0.01)
							mGnssData.latitude_std = 0.01;
						if (mGnssData.longitude_std < 0.01)
							mGnssData.longitude_std = 0.01;
						if (mGnssData.altitude_std < 0.01)
							mGnssData.altitude_std = 0.01;
					}
					else if (mGnssData.Mode == 5)
					{
						if (mGnssData.latitude_std < 0.05)
							mGnssData.latitude_std = 0.1;
						if (mGnssData.longitude_std < 0.05)
							mGnssData.longitude_std = 0.1;
						if (mGnssData.altitude_std < 0.05)
							mGnssData.altitude_std = 0.1;

					}
				}

				if ((mGnssData.latitude_std > 0 && mGnssData.latitude_std < 100)
						 &&(mGnssData.longitude_std > 0 && mGnssData.longitude_std < 100)
						&&(mGnssData.altitude_std > 0 && mGnssData.altitude_std < 100)
						&&(mGnssData.longitude > -180 * PI / 180 && mGnssData.longitude < 180 * PI / 180)
						&(mGnssData.latitude > -90 * PI / 180 && mGnssData.latitude < 90 *PI /180)
						&&(mGnssData.altitude > -5000 && mGnssData.altitude < 5000))
					{
					    if (IsStartGNSSWeek() != 1)
					    {
							SetStartGNSSWeek(mGnssData.week);
					    }
					    
						INSADDGNSSDATA(mGnssData);
#ifdef DEBUGMESSAGE
						ouputgnssfile(mGnssData,mLCSetting.gnssOutputDataRate, mLCSetting.kmlOutputDateRate);
#endif
						LastGNSSDataTime = mGnssData.timestamp;
						gnssPOSExistFlag = 0;
					}
					else
					{
					    mGnssData.gpsFixType = 0;
					    mGnssData.Mode = 0;
					    INSADDGNSSDATA(mGnssData);

#ifdef DEBUGMESSAGE
					//printf("GNSS Data: Incorrect data ,gnss_week: %d gnss_time0fweek: %.4f\n", mGnssData.week, mGnssData.timestamp);
#endif // DEBUGMESSAGE
						continue;
					}

				}

		}
		else if (mLCSetting.gnssDataType == 3)
		{

			
//double tv = 0;
//if ((mGnssData.timestamp > 459593 + tv && mGnssData.timestamp < 459593 + tv + 60)
//	|| (mGnssData.timestamp > 461105 + tv && mGnssData.timestamp < 461105 + tv + 60)
//	|| (mGnssData.timestamp > 461934 + tv && mGnssData.timestamp < 461934 + tv + 60)
//	|| (mGnssData.timestamp > 463480 + tv && mGnssData.timestamp < 463480 + tv + 60)
//	)
	//if ((mGnssData.timestamp > 462317.400 + tv && mGnssData.timestamp < 462317.400 + tv + 60)
	//	|| (mGnssData.timestamp > 460838.200 + tv && mGnssData.timestamp < 460838.200 + tv + 60)
	//	|| (mGnssData.timestamp > 461334.000 + tv && mGnssData.timestamp < 461334.000 + tv + 60)
	//	|| (mGnssData.timestamp > 458465.200 + tv && mGnssData.timestamp < 458465.200 + tv + 60)
	//	)
		//if ((mGnssData.timestamp > 459986.0000 + tv && mGnssData.timestamp < 459986.0000 + tv + 60)
		//	|| (mGnssData.timestamp > 461140.0000 + tv && mGnssData.timestamp < 461140.0000 + tv + 60)
		//	//|| (mGnssData.timestamp > 461334.000 + tv && mGnssData.timestamp < 461334.000 + tv + 60)
		//	//|| (mGnssData.timestamp > 458465.200 + tv && mGnssData.timestamp < 458465.200 + tv + 60)
		//	)
	//		if ((mGnssData.timestamp > 462412.0000 + tv && mGnssData.timestamp < 462412.0000 + tv + 60)
	//			|| (mGnssData.timestamp > 463260.0000 + tv && mGnssData.timestamp < 463260.0000 + tv + 60)
	//			//|| (mGnssData.timestamp > 461334.000 + tv && mGnssData.timestamp < 461334.000 + tv + 60)
	//			//|| (mGnssData.timestamp > 458465.200 + tv && mGnssData.timestamp < 458465.200 + tv + 60)
	//			)
 //   {
	//gnssPOSExistFlag = 0;
	//gnssVelExistFlag = 0;
	//continue;
 //    }

			if (gnssPOSExistFlag == 1 && gnssVelExistFlag == 1 )//&& mGnssData.timestamp >287746.000)
			{
	

				double t_outstart[1] = {0.0};

				double dt_v = 0;
				double t_outtime = 0;
				for (int i = 0; i < 1; i++)
				{
					if (mGnssData.timestamp > t_outstart[i]+ dt_v && mGnssData.timestamp < t_outstart[i]+ dt_v + t_outtime)
					{
						gnssPOSExistFlag = 0;
						gnssVelExistFlag = 0;					
					}
				}
				if (gnssVelExistFlag == 0)
				{
					continue;
				}

				if ((mGnssData.latitude_std >= 0 && mGnssData.latitude_std < 100)
					&& (mGnssData.longitude_std >= 0 && mGnssData.longitude_std < 100)
					&& (mGnssData.altitude_std >= 0 && mGnssData.altitude_std < 100)
					&& (mGnssData.longitude > -180 * PI / 180 && mGnssData.longitude < 180 * PI / 180)
					&(mGnssData.latitude > -90 * PI / 180 && mGnssData.latitude < 90 * PI / 180)
					&& (mGnssData.altitude > -1000 && mGnssData.altitude < 1000))
				{
					if (IsStartGNSSWeek() != 1)
					{
						SetStartGNSSWeek(mGnssData.week);
					}
#ifdef DEBUGMESSAGE
					ouputgnssfile(mGnssData, mLCSetting.gnssOutputDataRate, mLCSetting.kmlOutputDateRate);
#endif
					INSADDGNSSDATA(mGnssData);
				}
				else
				{
					mGnssData.gpsFixType = 0;
					mGnssData.Mode = 0;
#ifdef DEBUGMESSAGE
					printf("GNSS Data: Incorrect data in (*-process) file, GNSS Week: %d gnss_time0fweek: %.4f\n", mGnssData.week, mGnssData.timestamp);
#endif // DEBUGMESSAGE
					INSADDGNSSDATA(mGnssData);
					LastGNSSDataTime = mGnssData.timestamp;
					LastGNSSDataSYSTime = mImuData.timestamp;
					gnssPOSExistFlag = 0;
					gnssVelExistFlag = 0;
					continue;
				}
				LastGNSSDataTime = mGnssData.timestamp;
				LastGNSSDataSYSTime = mImuData.timestamp;
				gnssPOSExistFlag = 0;
				gnssVelExistFlag = 0;
			}
		}

		//if (mImuData.timestamp - LastGNSSDataSYSTime >= 4.0 &&
		//	mImuData.timestamp - mGnssData.timestamp >= 0.5
		//	&&LastGNSSDataTime != -1
		//	&& mGnssInsSystem.AlignCTime > 0
		//	&& mImuData.timestamp - mGnssInsSystem.AlignCTime > 200)
		//{
		//	mGnssData.timestamp = mImuData.timestamp;
		//	mGnssData.Mode = 0;
		//	INSADDGNSSDATA(mGnssData);
		//	
		//}
	}

	if (fdat)fclose(fdat);
	return 1;
}

int8_t playGnssImuFile(char* gnssfile ,char* imufile)
{
	int8_t ret = -1;

	FILE *fGNSS = NULL;
	FILE *fIMU  = NULL;


	char buffer[150] = { 0 };
	char  *val[20];
	int32_t sentencesize;

	int8_t gnssPOSExistFlag = -1;

	fopen_s(&fGNSS, gnssfile, "r");
	fopen_s(&fIMU, imufile, "r");

	if (NULL == fGNSS)
	{
#ifdef DEBUGMESSAGE
		printf("Cant OPEN GNSSFILE: %s\n", gnssfile);
#endif // DEBUGMESSAGE
		return ret;
	}
	if (NULL == fIMU)
	{
#ifdef DEBUGMESSAGE
		printf("Cant OPEN IMUFILE: %s\n", imufile);
#endif // DEBUGMESSAGE
		return ret;
	}
	while (!feof(fIMU))
	{
		/*读取IMU数据*/
		fgets(buffer, sizeof(buffer), fIMU);
		if (strlen(buffer) <= 0) continue;
		char *temp = strchr(buffer, '\n');
		if (temp != NULL) temp[0] = '\0';
		parse_fields(buffer, val, &sentencesize);
		if (sentencesize == 9)  //保证IMU数据完整性
		{
			memset(&mImuData, 0, sizeof(ImuData));
			mImuData.week = atoi(val[0]);
			mImuData.timestamp = atof(val[1]);
			mImuData.timestamped = atof(val[2]);
			mImuData.Accx = atof(val[3])* 9.780 / 9.806;
			mImuData.Accy = atof(val[4])* 9.780 / 9.806;
			mImuData.Accz = atof(val[5])* 9.780 /9.806;
			mImuData.Gyrox = atof(val[6]) * PI / 180;
			mImuData.Gyroy = atof(val[7]) * PI / 180;
			mImuData.Gyroz = atof(val[8]) * PI / 180;
			if (mImuData.week > 1024 && mImuData.week < 3072)
			{
				if (IsStartGNSSWeek() != 1)
				{
					SetStartGNSSWeek(mImuData.week);
				}
#ifdef DEBUGMESSAGE
				tracerawimu(&mImuData);
#endif
				INSAddIMUData(mImuData);
			}
			else
			{
				continue;
			}
		}
		else
		{
			continue;
		}
		while (!feof(fGNSS))
		{
			if (mGnssData.timestamp < mImuData.timestamp - 0.0001)
			{
				if (1 == gnssPOSExistFlag)
				{
					if ((mGnssData.latitude_std > 0 && mGnssData.latitude_std < 100)
						&& (mGnssData.longitude_std > 0 && mGnssData.longitude_std < 100)
						&& (mGnssData.altitude_std > 0 && mGnssData.altitude_std < 100)
						&& (mGnssData.longitude > -180 * PI / 180 && mGnssData.longitude < 180 * PI / 180)
						&& (mGnssData.latitude > -90 * PI / 180 && mGnssData.latitude < 90 * PI / 180)
						&& (mGnssData.altitude > -1000 && mGnssData.altitude < 1000))
					{
						if (IsStartGNSSWeek() != 1)
						{
							SetStartGNSSWeek(mGnssData.week);
						}
						mGnssData.timestampd = mImuData.timestamp;

						INSADDGNSSDATA(mGnssData);
#ifdef DEBUGMESSAGE
						ouputgnssfile(mGnssData, mLCSetting.gnssOutputDataRate, mLCSetting.kmlOutputDateRate);
						tracenativertk(&mGnssData);
#endif
					}
					else
					{
						mGnssData.gpsFixType = 0;
						mGnssData.Mode = 0;
						INSADDGNSSDATA(mGnssData);
#ifdef DEBUGMESSAGE
						printf("GNSS Data: Incorrect data in (GNSS + IMU) files, gnss_week: %d gnss_time0fweek: %.4f\n", mGnssData.week, mGnssData.timestamp);
#endif // DEBUGMESSAGE
					}
					gnssPOSExistFlag = 0;
					LastGNSSDataTime = mGnssData.timestamp;
					LastGNSSDataSYSTime = mImuData.timestamp;

				}
				fgets(buffer, sizeof(buffer), fGNSS);
				if (strlen(buffer) <= 0) continue;
				char *temp = strchr(buffer, '\n');
				if (temp != NULL) temp[0] = '\0';
				parse_fields(buffer, val, &sentencesize);
				if (sentencesize >= 9)
				{
					mGnssData.week = atoi(val[0]);
					mGnssData.timestamp = atof(val[1]);
					mGnssData.latitude = atof(val[2])*PI / 180;
					mGnssData.longitude = atof(val[3])*PI / 180;
					mGnssData.altitude = atof(val[4]);
					mGnssData.latitude_std = atof(val[5]);
					mGnssData.longitude_std = atof(val[6]);
					mGnssData.altitude_std = atof(val[7]);
					mGnssData.gpsFixType = atoi(val[8]);
					mGnssData.Mode = (int)atof(val[8]);
					mGnssData.north_velocity = atof(val[9]);
					mGnssData.east_velocity = atof(val[10]);
					mGnssData.up_velocity = atof(val[11]);
					gnssPOSExistFlag = 1;
				}
				else
				{
#ifdef DEBUGMESSAGE
					printf("GNSS Data: Format error  ,gnss_week: %d gnss_time0fweek: %.4f\n", mGnssData.week, mGnssData.timestamp);
#endif // DEBUGMESSAGE
					continue;
				}
			}
			else
			{
				break;
			}
		}
	}

	if (fGNSS)
	{
		fclose(fGNSS);
	}
	if (fIMU)
	{
		fclose(fIMU);
	}
	return 1;
}

int8_t process(const int8_t procType)
{
	switch (procType)
	{
	case 0:   // Simulate real-time data flow
	{
		playLogFile(mLCSetting.procfileNme);
	}break;
	case 1:   //Read GNSS data and IMU data respectively
	{
		playGnssImuFile(mLCSetting.gnssfileNme,mLCSetting.insfileName);
	}
	default:
		break;
	}

	return 1;
}

int8_t readconfig(const char* procfname, LCSetting*  mLCSetting)
{
	int8_t ret = -1;
	FILE *fdat = NULL;
	char buffer[1024] = { 0 };
	char  *val[40];
	int32_t sentencesize;
	fopen_s(&fdat, procfname, "r");
	if (!fdat)
	{
#ifdef DEBUGMESSAGE
		printf("Cant OPEN PROCFILE: %s\n", procfname);
#endif // DEBUGMESSAGE
		return ret;
	}

	while (!feof(fdat))
	{
		fgets(buffer, sizeof(buffer), fdat);
		if (strlen(buffer) <= 0) continue;
		char *temp = strchr(buffer, '\n');
		if (temp != NULL) temp[0] = '\0';

		parse_fields(buffer, val, &sentencesize);
		if (strstr(val[0], "$CONFIG") != NULL)
		{
			mLCSetting->procType = 0;
			memcpy(mLCSetting->procfileNme ,procfname, strlen(procfname) * sizeof(char));
			//通过模式判断是否使用 NHC  0：车载   1：机载   2：船载
			mLCSetting->isUseNHC = atoi(val[2]) == 0 ? 1 : 0;
			mLCSetting->initialAttitudeMode = atoi(val[3]);
			mLCSetting->isUseOdo = atoi(val[4]);

			mLCSetting->imuSensorType = atoi(val[5]);
			mLCSetting->imuDataRate = atoi(val[6]);
			mLCSetting->gnssSensorType = atoi(val[7]);
			mLCSetting->gnssDataType = atoi(val[8]);
			mLCSetting->gnssDataRate = atoi(val[9]);

			mLCSetting->priLeverArm[0] = atof(val[10]);
			mLCSetting->priLeverArm[1] = atof(val[11]);
			mLCSetting->priLeverArm[2] = atof(val[12]);

			mLCSetting->rotationRBV[0] = atof(val[13]);
			mLCSetting->rotationRBV[1] = atof(val[14]);
			mLCSetting->rotationRBV[2] = atof(val[15]);

			mLCSetting->attitueRPH[0]  = atof(val[16]);
			mLCSetting->attitueRPH[1]  = atof(val[17]);
			mLCSetting->attitueRPH[2]  = atof(val[18]);

			mLCSetting->userLeverArm[0] = atof(val[19]);
			mLCSetting->userLeverArm[1] = atof(val[20]);
			mLCSetting->userLeverArm[2] = atof(val[21]);

			if (fabs(mLCSetting->userLeverArm[0]) < 0.01
				&& fabs(mLCSetting->userLeverArm[0]) < 0.01
				&& fabs(mLCSetting->userLeverArm[0]) < 0.01)
			{
				mLCSetting->isUserOutPut == 0;
			}

			mLCSetting->odoLeverArm[0] = atof(val[22]);
			mLCSetting->odoLeverArm[1] = atof(val[23]);
			mLCSetting->odoLeverArm[2] = atof(val[24]);

			mLCSetting->isUseMisAlignment = 0;

			/*其他配置，但无需对外*/
			mLCSetting->isUseGNSSVel = 1;
			mLCSetting->odoScale = 1;
			mLCSetting->isUseZUPT = 1;
			mLCSetting->isUseZUPTLOCK = 1;
			mLCSetting->isUseHeadOnline = 1;
			mLCSetting->isUseExpQC = 1;
			/*后处理 输出频率*/
			mLCSetting->fileOutputType = 7;
			mLCSetting->profiletype = 255;
			mLCSetting->isKmlOutput = 1;
			mLCSetting->kmlOutputDateRate = 1;
			mLCSetting->gnssOutputDataRate = 1;
			mLCSetting->insOutputDataRate = 10;
			ret = 1;
			break;
		}
	}
	if (fdat)
	{
		fclose(fdat);
	}
	return ret;

}


int main(int argc, char **argv)
{
	if (3 == argc)
	{
		readconfigfromfile(argv[1], &mLCSetting);
		// process message
#ifdef DEBUGMESSAGE
		printf("OpenIns processing %s\n", mLCSetting.procfileNme);
		setfilename(mLCSetting.procfileNme, mLCSetting.fileOutputType, 
			mLCSetting.profiletype, mLCSetting.isKmlOutput, mLCSetting.isUseOdo);
		if (mLCSetting.profiletype & INSDEBUGMESSAGE)
		{
		}
#endif
	}
	else if (2 == argc)
	{
		readconfig(argv[1], &mLCSetting);
#ifdef DEBUGMESSAGE
		printf("OpenIns processing %s\n", mLCSetting.procfileNme);
		setfilename(mLCSetting.procfileNme, mLCSetting.fileOutputType, 
			mLCSetting.profiletype, mLCSetting.isKmlOutput, mLCSetting.isUseOdo);
		if (mLCSetting.profiletype & INSDEBUGMESSAGE)
		{
		}
#endif
	}
	 insinitsystemfromcfg(mLCSetting);
	 mLCSetting.isusefliter = 0;
	 process(mLCSetting.procType);

	 stop();

	 return 0;
}




