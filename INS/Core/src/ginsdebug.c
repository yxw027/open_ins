#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdio.h> 
#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include "insoutmsg.h"
#include "datatype.h"
#include "lcstruct.h"
#include "earth.h"
//#include "log.h"
//#include "com.h"
#include "ginsdebug.h"
#include "loosecoupleset.h"
#include "lcsystem.h"

#ifndef PI
#define PI 3.14159265358979
#endif
#define HEADKML1 "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
#define HEADKML2 "<kml xmlns=\"http://www.opengis.net/kml/2.2\">"
#define MARKICON "http://maps.google.com/mapfiles/kml/shapes/track.png"


#define LINUX 1
#define R2D   (180/PI)
#define SIZP     0.4            /* mark size of rover positions */
#define SIZR     0.4            /* mark size of reference position */
#define TINT     30.0           /* time label interval (sec) */
//#define SYS_OS   0

FILE *output = NULL;
FILE *output_12 = NULL;
FILE *output_debug = NULL;
FILE *output_kml = NULL;
FILE* output_GGA = NULL; /* output GGA for INS solution, always use 4 for temp */
FILE* output_GGA_GPS = NULL; /* output GGA for GPS input */
FILE* output_GP_GNSS = NULL; /* output GNSS for GPS input */
FILE* output_GP_KML = NULL;
FILE* output_CSV_GINS_12 = NULL; /* output GGA for GPS input */
FILE* output_CSV_GINS = NULL; /* output GGA for GPS input */
FILE* output_rtkcsv_GINS = NULL; /* output GGA for GPS input */
FILE* output_rtkcsv_GINS_rtk = NULL; /* output GGA for GPS input */
FILE* output_rawimu = NULL; /* output Rawimu */
FILE* output_corrimu = NULL; /* output corrimu */
FILE* output_rawodo = NULL; /* output Rawodo */
FILE* output_proc = NULL; /* output Rawodo */
FILE *output_mis = NULL;



extern int32_t gps_start_week;
typedef struct {        /* time struct */
	time_t time;        /* time (s) expressed by standard time_t */
	double sec;         /* fraction of second under 1 s */
} gtime_t;
#define MAXLEAPS    64                  /* max number of leap seconds table */
const static double gpst0[] = {1980,1,6,0,0,0 }; /* gps time reference */
static double leaps[MAXLEAPS + 1][7] = { /* leap seconds (y,m,d,h,m,s,utc-gpst) */
	{2017,1,1,0,0,0,-18},
	{2015,7,1,0,0,0,-17},
	{2012,7,1,0,0,0,-16},
	{2009,1,1,0,0,0,-15},
	{2006,1,1,0,0,0,-14},
	{1999,1,1,0,0,0,-13},
	{1997,7,1,0,0,0,-12},
	{1996,1,1,0,0,0,-11},
	{1994,7,1,0,0,0,-10},
	{1993,7,1,0,0,0, -9},
	{1992,7,1,0,0,0, -8},
	{1991,1,1,0,0,0, -7},
	{1990,1,1,0,0,0, -6},
	{1988,1,1,0,0,0, -5},
	{1985,7,1,0,0,0, -4},
	{1983,7,1,0,0,0, -3},
	{1982,7,1,0,0,0, -2},
	{1981,7,1,0,0,0, -1},
	{0}
};
static  gtime_t timeadd(gtime_t t, double sec)
{
	double tt;

	t.sec += sec; tt = floor(t.sec); t.time += (int)tt; t.sec -= tt;
	return t;
}
static void time2epoch(gtime_t t, double *ep)
{
	const int mday[] = { /* # of days in a month */
		31,28,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31,
		31,29,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31
	};
	int days, sec, mon, day;

	/* leap year if year%4==0 in 1901-2099 */
	days = (int)(t.time / 86400);
	sec = (int)(t.time - (time_t)days * 86400);
	for (day = days % 1461, mon = 0; mon < 48; mon++) {
		if (day >= mday[mon]) day -= mday[mon]; else break;
	}
	ep[0] = 1970 + days / 1461 * 4 + mon / 12; ep[1] = mon % 12 + 1; ep[2] = day + 1;
	ep[3] = sec / 3600; ep[4] = sec % 3600 / 60; ep[5] = sec % 60 + t.sec;
}
static gtime_t epoch2time(const double *ep)
{
	const int doy[] = { 1,32,60,91,121,152,182,213,244,274,305,335 };
	gtime_t time = { 0 };
	int days, sec, year = (int)ep[0], mon = (int)ep[1], day = (int)ep[2];

	if (year < 1970 || 2099 < year || mon < 1 || 12 < mon) return time;

	/* leap year if year%4==0 in 1901-2099 */
	days = (year - 1970) * 365 + (year - 1969) / 4 + doy[mon - 1] + day - 2 + (year % 4 == 0 && mon >= 3 ? 1 : 0);
	sec = (int)floor(ep[5]);
	time.time = (time_t)days * 86400 + (int)ep[3] * 3600 + (int)ep[4] * 60 + sec;
	time.sec = ep[5] - sec;
	return time;
}
static gtime_t gpst2time(int week, double sec)
 {
	 gtime_t t = epoch2time(gpst0);

	 if (sec < -1E9 || 1E9 < sec) sec = 0.0;
	 t.time += 86400 * 7 * week + (int)sec;
	 t.sec = sec - (int)sec;
	 return t;
 }
static gtime_t gpst2utc(gtime_t t)
 {
	 gtime_t tu;

	 ////for (i = 0; leaps[i][0] > 0; i++) {
		//// tu = timeadd(t, leaps[i][6]);
		//// if (timediff(tu, epoch2time(leaps[i])) >= 0.0) return tu;
	 ////}
	 tu = timeadd(t, leaps[0][6]);
	 return tu;
 }
static void  DegtoDMS(const double mdeg, int* Deg, int* Min, double* Sec)
 {
	 double AM;
	 *Deg = (int)mdeg;
	 AM = (mdeg - *Deg) * 60.0;
	 *Min = (int)AM;
	 *Sec = (AM - *Min) * 60.0;
 }
static void deg2dms(double deg, double* dms)
  {
	  double sign = deg < 0.0 ? -1.0 : 1.0, a = fabs(deg);
	  dms[0] = floor(a); a = (a - dms[0]) * 60.0;
	  dms[1] = floor(a); a = (a - dms[1]) * 60.0;
	  dms[2] = a; dms[0] *= sign;
  }

  /* output solution in the form of nmea GGA sentence --------------------------*/
 static int outnmea_gga(unsigned char* buff, double time, int type, double* blh, int ns, double dop, double age)
  {
	  double h, ep[6], dms1[3], dms2[3];
	  char* p = (char*)buff, * q, sum;

	  if (type != 1 && type != 4 && type != 5) {
#ifdef SYS_OS 
		  p += snprintf(p, 255, "$GPGGA,,,,,,,,,,,,,,");
#else
		  p += sprintf_s(p, 255, "$GPGGA,,,,,,,,,,,,,,");
#endif
		  for (q = (char*)buff + 1, sum = 0; *q; q++) sum ^= *q;
#ifdef SYS_OS 
		  p += snprintf(p, 255,"*%02X%c%c", sum, 0x0D, 0x0A);
#else
		  p += sprintf_s(p, 255,"*%02X%c%c", sum, 0x0D, 0x0A);
#endif
		  return p - (char*)buff;
	  }
	  time -= 18.0;
	  ep[2] = floor(time / (24 * 3600));
	  time -= ep[2] * 24 * 3600.0;
	  ep[3] = floor(time / 3600);
	  time -= ep[3] * 3600;
	  ep[4] = floor(time / 60);
	  time -= ep[4] * 60;
	  ep[5] = time;
	  h = 0.0;
	  deg2dms(fabs(blh[0]) * 180.0 / PI, dms1);
	  deg2dms(fabs(blh[1]) * 180.0 / PI, dms2);
#ifdef SYS_OS 
	  p += snprintf(p, 255, "$GPGGA,%02.0f%02.0f%05.2f,%02.0f%010.7f,%s,%03.0f%010.7f,%s,%d,%02d,%.1f,%.3f,M,%.3f,M,%.1f,",
		  ep[3], ep[4], ep[5], dms1[0], dms1[1] + dms1[2] / 60.0, blh[0] >= 0 ? "N" : "S",
		  dms2[0], dms2[1] + dms2[2] / 60.0, blh[1] >= 0 ? "E" : "W", type,
		  ns, dop, blh[2] - h, h, age);
#else
	  p += sprintf_s(p, 255, "$GPGGA,%02.0f%02.0f%05.2f,%02.0f%010.7f,%s,%03.0f%010.7f,%s,%d,%02d,%.1f,%.3f,M,%.3f,M,%.1f,",
		  ep[3], ep[4], ep[5], dms1[0], dms1[1] + dms1[2] / 60.0, blh[0] >= 0 ? "N" : "S",
		  dms2[0], dms2[1] + dms2[2] / 60.0, blh[1] >= 0 ? "E" : "W", type,
		  ns, dop, blh[2] - h, h, age);
#endif
	  for (q = (char*)buff + 1, sum = 0; *q; q++) sum ^= *q; /* check-sum */
#ifdef SYS_OS 
	  p += snprintf(p, 255,"*%02X%c%c", sum, 0x0D, 0x0A);
#else
	  p += sprintf_s(p, 255,"*%02X%c%c", sum, 0x0D, 0x0A);
#endif
	  return p - (char*)buff;
  }

int8_t setfilename(const char* inputname, 
	const int32_t fileOutputType, const int32_t profiletype, 
	const int32_t isKmlOutput, const int8_t useOdo)
{
	int8_t ret = -1;
	char filname[255];
#ifdef SYS_OS 
	strncpy(filname, inputname, 255);
#else
	strcpy_s(filname,255, inputname);
#endif
	char fileName[255] = { 0 }, dirname[255] = { 0 };
	char outfilename[255] = { 0 };
#ifdef SYS_OS 
	strncpy(dirname, filname, 255);
#else
	strcpy_s(dirname, 255, filname);
#endif
	char* result1 = strrchr(dirname, '\\');
	if (result1 == NULL) result1 = strrchr(dirname, '/');
	if (result1 != NULL)
	{
#ifdef SYS_OS 
		strncpy(fileName, result1 + 1, 255);
#else
		strcpy_s(fileName, 255, result1 + 1);
#endif
		result1[0] = '/';
		result1[1] = '\0';
	}
	else
	{
#ifdef SYS_OS 
		strncpy(fileName, dirname, 255);
#else
		strcpy_s(fileName, 255, dirname);
#endif
		dirname[0] = '\0';
	}

	result1 = strrchr(fileName, '.');
	if (result1 != NULL) result1[0] = '\0';

	if (fileOutputType & GNSSOUTPUT)
	{
#ifdef SYS_OS 
		snprintf(outfilename, 255, "%s%s_gps.nmea", dirname, fileName);
		output_GGA_GPS = fopen(outfilename, "w");
#else
		sprintf_s(outfilename, 255, "%s%s_gps.nmea", dirname, fileName);
		fopen_s(&output_GGA_GPS, outfilename, "w");
#endif

		if (isKmlOutput == 1)
		{
#ifdef SYS_OS 
			snprintf(outfilename, 255, "%s%s_gnss.kml", dirname, fileName);
			output_GP_KML = fopen(outfilename, "w");
#else           
			sprintf_s(outfilename, 255, "%s%s_gnss.kml", dirname, fileName);
			fopen_s(&output_GP_KML, outfilename, "w");
#endif
			char *color[] = {
			"ffffffff","ff0000ff","ffff00ff","50FF78F0","ff00ff00","ff00aaff"
			};  //B-G-R 白色 绿色 浅黄  红色  黄色 青色
			fprintf(output_GP_KML, "%s\n%s\n", HEADKML1, HEADKML2);
			fprintf(output_GP_KML, "<Document>\n");
			for (int i = 0; i < 6; i++) {
				fprintf(output_GP_KML, "<Style id=\"P%d\">\n", i);
				fprintf(output_GP_KML, "  <IconStyle>\n");
				fprintf(output_GP_KML, "    <color>%s</color>\n", color[i]);
				fprintf(output_GP_KML, "    <scale>%.1f</scale>\n", i == 0 ? SIZR : SIZP);
				fprintf(output_GP_KML, "    <Icon><href>%s</href></Icon>\n", MARKICON);
				fprintf(output_GP_KML, "  </IconStyle>\n");
				fprintf(output_GP_KML, "</Style>\n");
			}
		}

	}

	if (fileOutputType & INSOUTPUT)
	{
		if (profiletype & GGA)
		{
#ifdef SYS_OS 
			snprintf(outfilename, 255, "%s%s_ins.nmea", dirname, fileName);
			output_GGA = fopen(outfilename, "w");
#else		
			if (useOdo) {
				sprintf_s(outfilename, 255, "%s%s_ins_odo.nmea", dirname, fileName);
			}
			else {
				sprintf_s(outfilename, 255, "%s%s_ins.nmea", dirname, fileName);
			}
			fopen_s(&output_GGA, outfilename, "w");
#endif
		}

		/*if (profiletype & GGA)
		{
#ifdef SYS_OS 
			snprintf(outfilename, 255, "%s%s_ins.nmea", dirname, fileName);
			output_GGA = fopen(outfilename, "w");
#else
			sprintf_s(outfilename, 255, "%s%s_ins.nmea", dirname, fileName);
			fopen_s(&output_GGA, outfilename, "w");
#endif
		}*/

		if (profiletype & INSDEBUG1)
		{
#ifdef SYS_OS 
			snprintf(outfilename, 255, "%s%s_result.txt", dirname, fileName);
			output = fopen(outfilename, "w");
#else
			if (useOdo) {
				sprintf_s(outfilename, 255, "%s%s_result_odo.txt", dirname, fileName);
			}
			else {
				sprintf_s(outfilename, 255, "%s%s_result.txt", dirname, fileName);
			}
			fopen_s(&output, outfilename, "w");
#endif
		}

		if (profiletype & INSDEBUG2)
		{
#ifdef SYS_OS 
			snprintf(outfilename, 255, "%s%s_deb.txt", dirname, fileName);
			output_debug = fopen(outfilename, "w");
#else
			if (useOdo) {
				sprintf_s(outfilename, 255, "%s%s_deb_odo.txt", dirname, fileName);
			}
			else {
				sprintf_s(outfilename, 255, "%s%s_deb.txt", dirname, fileName);
			}
			fopen_s(&output_debug, outfilename, "w");
#endif
		}

		if (profiletype & CORRIMU)
		{
#ifdef SYS_OS 
			snprintf(outfilename, 255, "%s%s_corrimu.txt", dirname, fileName);
			output = fopen(outfilename, "w");
#else
			sprintf_s(outfilename, 255, "%s%s_corrimu.txt", dirname, fileName);
			fopen_s(&output_corrimu, outfilename, "w");
#endif
		}

		if (profiletype & TESLAGNSSCSV)
		{
#ifdef SYS_OS 
			snprintf(outfilename, 255, "%s%s_tesla_rover.csv", dirname, fileName);
			output_rtkcsv_GINS = fopen(outfilename, "w");
#else
			sprintf_s(outfilename, 255, "%s%s_tesla_rover.csv", dirname, fileName);
			fopen_s(&output_rtkcsv_GINS, outfilename, "w");
#endif
			fprintf(output_rtkcsv_GINS, "GPS_week,week second in ms,Error_N,Error_E,Error_U,Error_2D,Error_3D,Corrections Valid,Latitude,Longitude,Height_ell_m,HDOP,PDOP,Num_Sattelites,SDEast_m,SDNorth_m,SDHeight_m,Accuracy\n");
		}

		if (profiletype & TESLAREFCSV)
		{
#ifdef SYS_OS 
			snprintf(outfilename, 255, "%s%s_tesla_ref.csv", dirname, fileName);
			output_CSV_GINS = fopen(outfilename, "w");
#else
			sprintf_s(outfilename, 255, "%s%s_tesla_ref.csv", dirname, fileName);
			fopen_s(&output_CSV_GINS, outfilename, "w");
#endif
			fprintf(output_CSV_GINS, "      Week,  GPSTime,          Roll,         Pitch,       Heading,  VX-ECEF,  VY-ECEF,  VZ-ECEF,    VEast,   VNorth,      VUp, AngRateX, AngRateY, AngRateZ,AccBdyX,AccBdyY,AccBdyZ,      X-ECEF,      Y-ECEF,      Z-ECEF,        RollSD,       PitchSD,        HdngSD,      SDEast,     SDNorth,    SDHeight,    SD-VE,    SD-VN,    SD-VH,      Latitude,     Longitude,       H-Ell,NS,  HDOP");
			fprintf(output_CSV_GINS, "   (weeks),    (sec),         (deg),         (deg),         (deg),    (m/s),    (m/s),    (m/s),    (m/s),    (m/s),    (m/s),  (deg/s),  (deg/s),  (deg/s),(m/s^2),(m/s^2),(m/s^2),         (m),         (m),         (m),         (deg),         (deg),         (deg),         (m),         (m),         (m),    (m/s),    (m/s),    (m/s),         (deg),         (deg),         (m),  , (dop)");
		}

		if (isKmlOutput == 1)
		{
#ifdef SYS_OS 
			snprintf(outfilename, 255, "%s%s_ins.kml", dirname, fileName);
			output_kml = fopen(outfilename, "w");
#else
			if (useOdo) {
				sprintf_s(outfilename, 255, "%s%s_ins_odo.kml", dirname, fileName);
			}
			else {
				sprintf_s(outfilename, 255, "%s%s_ins.kml", dirname, fileName);
			}
			fopen_s(&output_kml, outfilename, "w");
#endif
			char *color[] = {
				// white ,cyan, purple, red, green, yellow
			"ffffffff","50FF78F0","ffff00ff","ff0000ff","ff00ff00","ff00aaff"
			};  //B-G-R 白色 绿色 浅黄  红色  黄色 青色
			fprintf(output_kml, "%s\n%s\n", HEADKML1, HEADKML2);
			fprintf(output_kml, "<Document>\n");
			for (int i = 0; i < 6; i++) {
				fprintf(output_kml, "<Style id=\"P%d\">\n", i);
				fprintf(output_kml, "  <IconStyle>\n");
				fprintf(output_kml, "    <color>%s</color>\n", color[i]);
				fprintf(output_kml, "    <scale>%.1f</scale>\n", i == 0 ? SIZR : SIZP);
				fprintf(output_kml, "    <Icon><href>%s</href></Icon>\n", MARKICON);
				fprintf(output_kml, "  </IconStyle>\n");
				fprintf(output_kml, "</Style>\n");
			}
		}
	}

	if (fileOutputType & RAWMEAS)
	{
#ifdef SYS_OS 
		snprintf(outfilename, 255, "%s%s_rawimu.txt", dirname, fileName);
		output_rawimu = fopen(outfilename, "w");
#else
		sprintf_s(outfilename, 255, "%s%s_rawimu.txt", dirname, fileName);
		fopen_s(&output_rawimu, outfilename, "w");
#endif

#ifdef SYS_OS 
		snprintf(outfilename, 255, "%s%s_rawodo.txt", dirname, fileName);
		output_rawodo = fopen(outfilename, "w");
#else
		sprintf_s(outfilename, 255, "%s%s_rawodo.txt", dirname, fileName);
		fopen_s(&output_rawodo, outfilename, "w");
#endif

#ifdef SYS_OS 
		snprintf(outfilename, 255, "%s%s_proc.txt", dirname, fileName);
		output_proc = fopen(outfilename, "w");
#else
		sprintf_s(outfilename, 255, "%s%s_proc.txt", dirname, fileName);
		fopen_s(&output_proc, outfilename, "w");
#endif


	}

	if (fileOutputType & MISALIGMENT)
	{
#ifdef SYS_OS 
		snprintf(outfilename, 255, "%s%s_mis.txt", dirname, fileName);
		output_mis = fopen(outfilename, "w");
#else
		sprintf_s(outfilename, 255, "%s%s_mis.txt", dirname, fileName);
		fopen_s(&output_mis, outfilename, "w");
#endif
	}


	return 1;
}
int ouputgnssfile(const GnssData mGnssData, const int32_t gnssOutputDataRate, const int32_t kmlOutputDateRate)
{
	if (output_GGA_GPS)
	{
		unsigned char ggaBuffer[400] = { 0 };
		double blh[3] = { mGnssData.latitude,mGnssData.longitude, mGnssData.altitude };
		int len = outnmea_gga(ggaBuffer, mGnssData.timestamp ,mGnssData.Mode, blh, 10, 1.0, 1.0);
		fprintf(output_GGA_GPS, "%s", ggaBuffer);
	}
	if (output_GP_KML && fmod(mGnssData.timestamp + 0.01, 1.0/kmlOutputDateRate) < 0.05)
	{
		double ep[6] = { 0.0 }, alt = 0.0;
		char str[256] = "";
		double timestamp = mGnssData.timestamp;
		gtime_t gpstime = gpst2time(gps_start_week, timestamp);
		gtime_t utctime = gpst2utc(gpstime);
		time2epoch(utctime, ep);
			fprintf(output_GP_KML, "<Placemark>\n");
			int style = 1;
#ifdef SYS_OS 
			snprintf(str, 255, "%04.0f-%02.0f-%02.0fT%02.0f:%02.0f:%06.3fZ",
				ep[0], ep[1], ep[2], ep[3], ep[4], ep[5]);
#else
			sprintf_s(str, 255, "%04.0f-%02.0f-%02.0fT%02.0f:%02.0f:%06.3fZ",
				ep[0], ep[1], ep[2], ep[3], ep[4], ep[5]);
#endif
			fprintf(output_GP_KML, "<TimeStamp><when>%s</when></TimeStamp>\n", str);
			fprintf(output_GP_KML, "<description><![CDATA[\n");
			fprintf(output_GP_KML, "<TABLE border=\"1\" width=\"100%\" Align=\"center\">\n");
			fprintf(output_GP_KML, "<TR ALIGN=RIGHT>\n");
			fprintf(output_GP_KML, "<TR ALIGN = RIGHT><TD ALIGN = LEFT>Time:</TD><TD>");
			fprintf(output_GP_KML, "%4d</TD><TD>", gps_start_week);
			fprintf(output_GP_KML, "%11.4f</TD><TD>", mGnssData.timestamp);
			fprintf(output_GP_KML, "%02.0f:%02.0f:%06.3f</TD><TD>", ep[3], ep[4], ep[5]);
			fprintf(output_GP_KML, "%2d/%2d/%3d</TD></TR>\n", (int)ep[0], (int)ep[1], (int)ep[2]);
			fprintf(output_GP_KML, "<TR ALIGN=RIGHT><TD ALIGN=LEFT>Position:</TD><TD>");
			fprintf(output_GP_KML, "%11.7f</TD><TD>", mGnssData.latitude * 180 / PI);
			fprintf(output_GP_KML, "%11.7f</TD><TD>", mGnssData.longitude * 180 / PI);
			fprintf(output_GP_KML, "%8.4f</TD>", mGnssData.altitude);
			fprintf(output_GP_KML, "<TD>(DMS,m)</TD></TR>\n");
			fprintf(output_GP_KML, "<TR ALIGN=RIGHT><TD ALIGN=LEFT>STD(N,E,D):</TD><TD>");
			fprintf(output_GP_KML, "%8.4f</TD><TD>", mGnssData.latitude_std);
			fprintf(output_GP_KML, "%8.4f</TD><TD>", mGnssData.longitude_std);
			fprintf(output_GP_KML, "%8.4f</TD>", mGnssData.altitude_std);
			fprintf(output_GP_KML, "<TD>(m/s)</TD></TR>\n");
			fprintf(output_GP_KML, "<TR ALIGN=RIGHT><TD ALIGN=LEFT>Att(r,p,y):</TD><TD>");
			fprintf(output_GP_KML, "%d</TD><TD>", mGnssData.Mode);
			fprintf(output_GP_KML, "%8.4f</TD><TD>", -100);
			fprintf(output_GP_KML, "%8.4f</TD>", -100);
			fprintf(output_GP_KML, "<TD>(deg)</TD></TR>\n");
			fprintf(output_GP_KML, "</TABLE>\n");
			fprintf(output_GP_KML, "]]></description>\n");

			fprintf(output_GP_KML, "<styleUrl>#P%d</styleUrl>\n", style);
			fprintf(output_GP_KML, "<Style>\n");
			fprintf(output_GP_KML, "<IconStyle>\n");
			fprintf(output_GP_KML, "<heading>%f</heading>\n", -100);
			fprintf(output_GP_KML, "</IconStyle>\n");
			fprintf(output_GP_KML, "</Style>\n");


			fprintf(output_GP_KML, "<Point>\n");
			fprintf(output_GP_KML, "<coordinates>%13.9f,%12.9f,%5.3f</coordinates>\n", mGnssData.longitude * R2D,
				mGnssData.latitude * R2D, mGnssData.altitude);
			fprintf(output_GP_KML, "</Point>\n");
			fprintf(output_GP_KML, "</Placemark>\n");
		
	}
	return 1;
}

static ins_status[] = {"INS_INACTIVE", "INS_ALIGNING", "INS_HIGH_VARIANCE", "INS_SOLUTION_GOOD", "INS_SOLUTION_FREE", "INS_ALIGNMENT_COMPLETE"};
static ins_postype[] = { "INS_NONE", "INS_PSRSP", "INS_PSRDIFF", "INS_PROPOGATED", "INS_RTKFIXED", "INS_RTKFLOAT"};
int8_t ouputfile(const GnssInsSystem mGnssInsSystem,const int32_t insOutputDataRate, const int32_t kmlOutputDateRate,const int32_t imuDataRate)
{
	int8_t ret = -1;
	double ep[6] = {0.0}, alt = 0.0;
	char str[256] = "";
	double timestamp = mGnssInsSystem.mImuData.timestamp;
	gtime_t gpstime = gpst2time(gps_start_week, timestamp);
	gtime_t utctime = gpst2utc(gpstime);
	time2epoch(utctime, ep);
	double P_ins[16][16] = { 0.0 };
	for (int i = 0; i < mGnssInsSystem.mKalmanStruct.n; i++)
	{
		P_ins[i][i] = mGnssInsSystem.mKalmanStruct.P[i * mGnssInsSystem.mKalmanStruct.n + i];
	}
	if (output_GGA && fmod(mGnssInsSystem.mImuData.timestamp + 0.5 / imuDataRate, 1.0) < 1.0 / imuDataRate)
	{
		unsigned char ggaBuffer[400] = { 0 };
		double blh[3] = { mGnssInsSystem.outNav.lat, mGnssInsSystem.outNav.lon, mGnssInsSystem.outNav.height };
		int len = outnmea_gga(ggaBuffer, timestamp, 4, blh, 10, 1.0, 1.0);
		fprintf(output_GGA, "%s", ggaBuffer);
	}
	if (output && fmod(mGnssInsSystem.mImuData.timestamp + 0.5 / imuDataRate, 1.0/insOutputDataRate) < 1.0 / imuDataRate)
	{
		fprintf(output, "%4d,%14.10f,%14.10f,%14.10f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%d,%d,%d,%d,%d,%10.5f\n",
			mGnssInsSystem.mImuData.week,mGnssInsSystem.mImuData.timestamp, mGnssInsSystem.outNav.lat * 180 / PI, mGnssInsSystem.outNav.lon * 180 / PI, mGnssInsSystem.outNav.height,
			mGnssInsSystem.outNav.vn, mGnssInsSystem.outNav.ve, mGnssInsSystem.outNav.vd,
			mGnssInsSystem.outNav.roll * 180 / PI, mGnssInsSystem.outNav.pitch * 180 / PI, mGnssInsSystem.outNav.heading * 180 / PI,
			mGnssInsSystem.outNav.sensorbias.bias_gyro_x * 180 / PI * 3600, mGnssInsSystem.outNav.sensorbias.bias_gyro_y * 180 / PI * 3600, mGnssInsSystem.outNav.sensorbias.bias_gyro_z * 180 / PI * 3600,
			mGnssInsSystem.outNav.sensorbias.bias_acc_x * 1.0e5, mGnssInsSystem.outNav.sensorbias.bias_acc_y * 1.0e5, mGnssInsSystem.outNav.sensorbias.bias_acc_z * 1.0e5,
			mGnssInsSystem.outPerMeasUpdata,
			(int)mGnssInsSystem.CurIsUseZupt,
			mGnssInsSystem.ins_status,
			mGnssInsSystem.ins_positin_type,
			mGnssInsSystem.isUseOdo,
			mGnssInsSystem.mNav.Odo_scale
		);
	}
	if (output_corrimu)// && fmod(mGnssInsSystem.mImuData.timestamp + 0.5 / imuDataRate, 1.0 / insOutputDataRate) < 1.0 / imuDataRate)
	{
		fprintf(output_corrimu, "%4d,%10.4f,%10.4f,%14.4f,%14.4f,%14.4f,%14.4f,%14.4f,%14.4f\n",
			mGnssInsSystem.mImuData.week, mGnssInsSystem.mImuData.timestamp,
			mGnssInsSystem.outNav.a_b[0], mGnssInsSystem.outNav.a_b[1], mGnssInsSystem.outNav.a_b[2],
			mGnssInsSystem.outNav.w_b[0]* 180/PI, mGnssInsSystem.outNav.w_b[1] * 180 / PI, mGnssInsSystem.outNav.w_b[2] * 180 / PI
			);
	}
	if (output_debug && fmod(mGnssInsSystem.mImuData.timestamp + 0.5 / imuDataRate, 1.0/insOutputDataRate) < 1.0 / imuDataRate)
	{
		fprintf(output_debug, "%4d, %14.10f,%14.10f,%14.10f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f\n",
			mGnssInsSystem.mImuData.week,mGnssInsSystem.mImuData.timestamp,
			sqrt(P_ins[0][0]), sqrt(P_ins[1][1]), sqrt(P_ins[2][2]),
			sqrt(P_ins[3][3]), sqrt(P_ins[4][4]), sqrt(P_ins[5][5]),
			sqrt(P_ins[6][6]) * 180 / PI, sqrt(P_ins[7][7]) * 180 / PI, sqrt(P_ins[8][8]) * 180 / PI,
			sqrt(P_ins[9][9]) * 180 / PI * 3600, sqrt(P_ins[10][10]) * 180 / PI * 3600, sqrt(P_ins[11][11]) * 180 / PI * 3600,
			sqrt(P_ins[12][12]) * 1.0e5, sqrt(P_ins[13][13]) * 1.0e5, sqrt(P_ins[14][14]) * 1.0e5, sqrt(P_ins[15][15])
		);
	}
	if (output_rtkcsv_GINS && fmod(mGnssInsSystem.mImuData.timestamp + 0.5 / imuDataRate, 1.0/insOutputDataRate) < 1.0 / imuDataRate)
	{
		fprintf(output_rtkcsv_GINS, "%d,%14.4f,0,0,0,0,0,1,%14.10f,%14.10f,%10.5f,0,0,0,%10.5f,%10.5f,%10.5f,%d\n",
			mGnssInsSystem.mImuData.week, mGnssInsSystem.mImuData.timestamp * 1000,
			mGnssInsSystem.outNav.lat * 180 / PI, mGnssInsSystem.outNav.lon * 180 / PI, mGnssInsSystem.outNav.height,
			sqrt(P_ins[1][1]), sqrt(P_ins[0][0]), sqrt(P_ins[2][2]), mGnssInsSystem.ins_positin_type);
	}
	if (output_CSV_GINS && fmod(mGnssInsSystem.mImuData.timestamp + 0.5 / imuDataRate, 1.0/insOutputDataRate) < 1.0 / imuDataRate)
	{
		fprintf(output_CSV_GINS, "%d,%14.4f,%10.5f,%10.5f,%10.5f,0,0,0,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,0,0,0,%10.5f,%10.5f,%10.5f, %10.5f, %10.5f, %10.5f, %10.5f, %10.5f, %10.5f, %14.10f, %14.10f, %10.5f ,0,0\n",
			mGnssInsSystem.mImuData.week, mGnssInsSystem.mImuData.timestamp,
			mGnssInsSystem.outNav.roll * 180 / PI, mGnssInsSystem.outNav.pitch * 180 / PI, mGnssInsSystem.outNav.heading * 180 / PI,

			mGnssInsSystem.outNav.ve, mGnssInsSystem.outNav.vn, -mGnssInsSystem.outNav.vd,
			mGnssInsSystem.outNav.wnb_n[0] * 180 / PI, mGnssInsSystem.outNav.wnb_n[1] * 180 / PI, mGnssInsSystem.outNav.wnb_n[2] * 180 / PI,
			mGnssInsSystem.outNav.a_n[0], mGnssInsSystem.outNav.a_n[1], mGnssInsSystem.outNav.a_n[2],

			sqrt(P_ins[6][6]) * 180 / PI, sqrt(P_ins[7][7]) * 180 / PI, sqrt(P_ins[8][8]) * 180 / PI,
			sqrt(P_ins[1][1]), sqrt(P_ins[0][0]), sqrt(P_ins[2][2]),
			sqrt(P_ins[4][4]), sqrt(P_ins[3][3]), sqrt(P_ins[5][5]),
			mGnssInsSystem.outNav.lat * 180 / PI, mGnssInsSystem.outNav.lon * 180 / PI, mGnssInsSystem.outNav.height
		);
	}
	if (output_kml && fmod(mGnssInsSystem.mImuData.timestamp + 0.5 / imuDataRate, 1.0/kmlOutputDateRate) < 1.0 / imuDataRate)
	{
			fprintf(output_kml, "<Placemark>\n");
			int style = 1;
			switch (mGnssInsSystem.ins_positin_type)  /*///			
               INS_NONE,
				INS_PSRSP,
				INS_PSRDIFF,
				INS_RTKFLOAT,
				INS_RTKFIXED*/
			{
				//B-G-R # white, cyan, red, purple, light-yellow, green, yellow
			case 0:
				//ros_msg->status = "INS_INACTIVE";
				style = 0;
				break;
			case 1:
				//ros_msg->status = "INS_ALIGNING";
				style = 1;
				break;
			case 2:
				//ros_msg->status = "INS_HIGH_VARIANCE";
				style = 2;
				break;
			case 3:
				//ros_msg->status = "INS_SOLUTION_GOOD";
				style = 3;
				break;
			case 4:
				//ros_msg->status = "INS_SOLUTION_FREE";
				style = 4;
				break;
			case 5:
				//ros_msg->status = "INS_SOLUTION_FREE";
				style = 5;
				break;
			default:
				break;
	}
#ifdef SYS_OS 
			snprintf(str,255, "%04.0f-%02.0f-%02.0fT%02.0f:%02.0f:%06.3fZ",
				ep[0], ep[1], ep[2], ep[3], ep[4], ep[5]);
#else
			sprintf_s(str,255, "%04.0f-%02.0f-%02.0fT%02.0f:%02.0f:%06.3fZ",
				ep[0], ep[1], ep[2], ep[3], ep[4], ep[5]);
#endif
			fprintf(output_kml, "<TimeStamp><when>%s</when></TimeStamp>\n", str);
			fprintf(output_kml, "<description><![CDATA[\n");
			fprintf(output_kml, "<TABLE border=\"1\" width=\"100%\" Align=\"center\">\n");
			fprintf(output_kml, "<TR ALIGN=RIGHT>\n");
			fprintf(output_kml, "<TR ALIGN = RIGHT><TD ALIGN = LEFT>Time:</TD><TD>");
			fprintf(output_kml, "%4d</TD><TD>", gps_start_week);
			fprintf(output_kml, "%11.4f</TD><TD>", mGnssInsSystem.mImuData.timestamp);
			fprintf(output_kml, "%02.0f:%02.0f:%06.3f</TD><TD>", ep[3], ep[4], ep[5]);
			fprintf(output_kml, "%2d/%2d/%3d</TD></TR>\n", (int)ep[0], (int)ep[1], (int)ep[2]);
			fprintf(output_kml, "<TR ALIGN=RIGHT><TD ALIGN=LEFT>Position:</TD><TD>");
			fprintf(output_kml, "%11.7f</TD><TD>", mGnssInsSystem.outNav.lat * 180 / PI);
			fprintf(output_kml, "%11.7f</TD><TD>", mGnssInsSystem.outNav.lon * 180 / PI);
			fprintf(output_kml, "%8.4f</TD>", mGnssInsSystem.outNav.height);
			fprintf(output_kml, "<TD>(DMS,m)</TD></TR>\n");
			fprintf(output_kml, "<TR ALIGN=RIGHT><TD ALIGN=LEFT>Vel(N,E,D):</TD><TD>");
			fprintf(output_kml, "%8.4f</TD><TD>", mGnssInsSystem.outNav.vn);
			fprintf(output_kml, "%8.4f</TD><TD>", mGnssInsSystem.outNav.ve);
			fprintf(output_kml, "%8.4f</TD>", mGnssInsSystem.outNav.vd);
			fprintf(output_kml, "<TD>(m/s)</TD></TR>\n");
			fprintf(output_kml, "<TR ALIGN=RIGHT><TD ALIGN=LEFT>Att(r,p,y):</TD><TD>");
			fprintf(output_kml, "%8.4f</TD><TD>", mGnssInsSystem.outNav.roll * 180 / PI);
			fprintf(output_kml, "%8.4f</TD><TD>", mGnssInsSystem.outNav.pitch * 180 / PI);
			fprintf(output_kml, "%8.4f</TD>", mGnssInsSystem.outNav.heading * 180 / PI);
			fprintf(output_kml, "<TD>(deg)</TD></TR>\n");
			fprintf(output_kml, "<TR ALIGN=RIGHT><TD ALIGN=LEFT>Misc Info:</TD><TD>");
			fprintf(output_kml, "%s</TD><TD>", ins_status[mGnssInsSystem.ins_status]);
			fprintf(output_kml, "%s</TD><TD>", ins_postype[style]);
			//fprintf(output_kml, "%d</TD>", 0);
			fprintf(output_kml, "<TD></TD></TR>\n");
			fprintf(output_kml, "</TABLE>\n");
			fprintf(output_kml, "]]></description>\n");

			fprintf(output_kml, "<styleUrl>#P%d</styleUrl>\n", style);
			fprintf(output_kml, "<Style>\n");
			fprintf(output_kml, "<IconStyle>\n");
			fprintf(output_kml, "<heading>%f</heading>\n", mGnssInsSystem.outNav.heading * 180 / PI);
			fprintf(output_kml, "</IconStyle>\n");
			fprintf(output_kml, "</Style>\n");


			fprintf(output_kml, "<Point>\n");
			fprintf(output_kml, "<coordinates>%13.9f,%12.9f,%5.3f</coordinates>\n", mGnssInsSystem.outNav.lon * R2D,
				mGnssInsSystem.outNav.lat * R2D, mGnssInsSystem.outNav.height);
			fprintf(output_kml, "</Point>\n");
			fprintf(output_kml, "</Placemark>\n");
	}
	return 1;
}

int8_t tracemis(const ImuData* pImuData, const int8_t* p1,const double*p2)
{
	//fprintf(output_mis, "%4d,%10.4f,%10.4f,%d,%d,%d,%14.4f,%14.4f,%14.4f\n",
	//	pImuData->week, pImuData->timestamp, pImuData->timestamped,
	//	(int)(p1[0]), (int)(p1[1]), (int)(p1[2]),
	//	p2[0] * 180 / PI, p2[1] * 180 / PI, p2[2] * 180 / PI);
}

int8_t tracemis1(const double *a)
{
	if (output_mis)
	{
		fprintf(output_mis, "%10.4f,%10.4f,%10.4f,%10.4f,%10.4f,%10.4f,%10.4f,%10.4f,%10.4f,%10.4f%10.4f\n",
			a[0], a[1] * 180 / PI, a[2] * 180 / PI,
			a[3], a[4], a[5],
			a[6] * 180 / PI, a[7] * 180 / PI, a[8] * 180 / PI, a[9] * 180 / PI
			, a[10]);
	}

}

int8_t tracerawimu(const ImuData* pImuData)
{
	if (output_rawimu)
	fprintf(output_rawimu, "%4d,%10.4f,%10.4f,%14.4f,%14.4f,%14.4f,%14.4f,%14.4f,%14.4f\n",
		pImuData->week,pImuData->timestamp,pImuData->timestamped,
		pImuData->Accx,pImuData->Accy,pImuData->Accz,
		pImuData->Gyrox * 180 /PI, pImuData->Gyroy * 180 / PI, pImuData->Gyroz * 180 / PI
	);


}
int8_t tracerawodo(const OdoData* pOdoData)
{
	if(output_rawodo)
	fprintf(output_rawodo, "%4d,%10.4f,%10.4f,%14.4f\n",
		pOdoData->week, pOdoData->timestamp, pOdoData->timestampd,
		pOdoData->vehicle_speed
	);
}

int8_t tracenativegnss(const double pnt[15])
{
	if(output_proc)
	fprintf(output_proc, "$GPSPP,%4.0f,%10.4f,%10.4f,%14.4f\n",
		pnt[0], pnt[1], pnt[2], pnt[3]
	);
	return 1;
}

int8_t tracenativeodo(const double odo[15])
{
	if(output_proc)
	fprintf(output_proc, "$GPODO,%4d,%10.4f,%10.4f,%14.4f\n",
		odo[0], odo[1], odo[2], odo[3]
	);
	return 1;
}

int8_t tracenativeimu(const ImuData* pImuData)
{
	if(output_proc)
	fprintf(output_proc, "GPIMU,%4d,%10.4f,%10.4f,%14.4f,%14.4f,%14.4f,%14.4f,%14.4f,%14.4f\n",
		pImuData->week, pImuData->timestamp, pImuData->timestamped,
		pImuData->Accx, pImuData->Accy, pImuData->Accz,
		pImuData->Gyrox * 180 / PI, pImuData->Gyroy * 180 / PI, pImuData->Gyroz * 180 / PI
	);
	return 1;
}

int8_t tracenativertk(const GnssData* pgnssData)
{
	if(output_proc)
	fprintf(output_proc, "GPRTK,%4d,%10.4f,%10.4f,%14.4f,%14.4f,%14.4f\n",
		pgnssData->week, pgnssData->timestamp, pgnssData->timestampd,
		pgnssData->latitude, pgnssData->longitude, pgnssData->altitude
	);
	return 1;
}



int8_t file_close()
{
	int8_t ret = -1;
	if (output)
	{
		fclose(output);
	}
	if (output_12) fclose(output_12);

	if (output_debug) fclose(output_debug);
	if (output_GGA) fclose(output_GGA);
	if (output_GGA_GPS) fclose(output_GGA_GPS);
	if (output_CSV_GINS_12) fclose(output_CSV_GINS_12);
	if (output_CSV_GINS) fclose(output_CSV_GINS);
	if (output_rtkcsv_GINS) fclose(output_rtkcsv_GINS);
	if (output_rtkcsv_GINS_rtk) fclose(output_rtkcsv_GINS_rtk);

	if (output_rawimu) fclose(output_rawimu);
	if (output_rawodo) fclose(output_rawodo);

	if (output_proc) fclose(output_proc);
	if (output_corrimu) fclose(output_corrimu);

	
	if (output_mis) fclose(output_mis);


	if (output_kml)
	{
		fprintf(output_kml, "</Document>\n");
		fprintf(output_kml, "</kml>\n");
		if (output_kml)
		{
			fclose(output_kml);
		}
	}
	if(output_GP_KML)
	{
		fprintf(output_GP_KML, "</Document>\n");
		fprintf(output_GP_KML, "</kml>\n");
		if (output_GP_KML)
		{
			fclose(output_GP_KML);
		}
	}
	return 1;
}
