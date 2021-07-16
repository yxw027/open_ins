#ifndef _GNSS_INS_H_
#define _GNSS_INS_H_

#include "stdint.h"
#include <stdio.h> 

#ifdef __cplusplus
extern "C" {
#endif

/*-------------------------------------------------------------*/
/* define constants */
#ifndef PI
#define PI (3.1415926535897932384626433832795)
#endif

/*-------------------------------------------------------------*/
/* define math functions */

uint8_t MatrixAdd(const double *matrix_a, const double *matrix_b, const int matrix_a_row, const int matrix_a_colume, double *matrix_result);
uint8_t MatrixSub(const double *matrix_a, const double *matrix_b, const int matrix_a_row, const int matrix_a_colume, double *matrix_result);
uint8_t MatrixMutiply(const double *matrix_a, const double *matrix_b, const int matrix_a_row, const int matrix_a_colume, const int matrix_b_column, double *matrix_result);
uint8_t MatrixTranspose(const double *matrix_a, const int matrix_row, const int matrix_column, double *matrix_result);
uint8_t MatrixCholosky(const double *matrix, const int matrix_row, double *lower_tri_matrix);
//uint8_t MatrixInverse(const double *lower_tri_matrix, const int matrix_row, double *matrix_result);
uint8_t CrossProduct(double a[3], double b[3], double c[3]);
uint8_t GetSkewSymmetricMatrixOfVector(const double  mvector[3], double *M);
uint8_t quatprod(const double p[4], const double q[4], double r[4]);
uint8_t rvec2quat(double rot_vect[3], double q[4]);
uint8_t quat2pos(double q[4], double *r);
uint8_t norm_quat(double q[4]);
uint8_t quatinv(const double p[4], double q[4]);
uint8_t quat2rvec(double q[4], double rot_vec[3]);
uint8_t MatrixInverse(uint8_t n,double* a); //Matrix Inverse

/*-------------------------------------------------------------*/
/* define struct for coordinate computation */
typedef  struct EarthParameter
{
	double a ;       //Ellipsoid long axis
	double b ;       //Ellipsoid short axis
	double f;       //Ellipsoidal oblate 
	double e;       //first Eccentricity of Elliopsoid 
	double e2;
	//double ep;
	//double ep2;     //second Eccentricity of Elliopsoid 
	double wie;     //rotational angular velocity of the earths  
	double GM;      //geocentric gravitational constant 
} EarthParameter;

typedef struct GravityParameter
{
	double g0;
	double g1;
	double g2;
	double g3;
	double g4;
	double g5;
} GravityParameter;


uint8_t UpdateMN(const double *BLH, double *M, double *N);
uint8_t UpdateW(const double *BLH, const double M, const double N, const double *vn, double *wnie, double *wnen);

uint8_t UpdateGravity(const double *BLH, double* ng);

//Earth-Centered, Earth-Fixed
typedef struct ECEF
{
	double x;
	double y;
	double z;
} ECEF;
//Geographic coordinates 
typedef struct Geo
{
	double lat;
	double lon;
	double height;
} Geo;

int8_t llh2ecef(const Geo* mGeo, ECEF* mECEF);
int8_t ecef2llh(const ECEF* mECEF,Geo* mGeo);


/*-------------------------------------------------------------*/
/* define basic imu, gps and odometer's data struct */

typedef struct  ImuData
{
	double timestamp;  //GPS time 
	double timestamped; //received time 
	double Gyrox;	// angle velocity along x axis  rad/s
	double Gyroy;	// angle velocity along y axis  rad/s
	double Gyroz;	// angle velocity along z axis  rad/s

	double Accx;	// line acceleration along x axis  m/s^2
	double Accy;	// line acceleration along y axis  m/s^2
	double Accz;	// line acceleration along z axis  m/s^2
}ImuData;

typedef struct  GnssData
{
	double flag;
	double timestamp;  //GPS time 
	double timestampd; //received time

	double longitude;	// degrees
	double latitude;	// degrees
	double altitude;	// ellipsoidal height - WGS84 (m)

	double north_velocity;	// m/s
	double east_velocity;	// m/s
	double up_velocity;		// m/s

	double longitude_std;	// longitude standard deviation
	double latitude_std;	// latitude standard deviation
	double altitude_std;	// altitude standard deviation

	double north_velocity_std;	// velocity standard deviation
	double east_velocity_std;
	double up_velocity_std;

	double length; //Dual antenna baseline length
	double pitch;
	double heading; //Dual antenna

	double pitch_std;
	double heading_std; //Dual antenna

	uint16_t Num_track; //number of stars visible
	uint16_t Num_sol; //number of stars visible

	uint16_t average_snr;
	double HDOP;
	double PDOP;
	uint16_t Mode;//0:NGNSS  1:spp, 2:PSR,  4:fixed,5:float
}GnssData;

typedef struct OdoData
{
	double timestamp;   //GNSS time
	double timestampd;  //systemtime
	double vehicle_speed;  //m/s
}OdoData;

/*-------------------------------------------------------------*/
typedef struct  INS
{
	double timestamp;
	double timestamped;
	double gyro_x;      // rad/s
	double gyro_y;
	double gyro_z;
	double acc_x;       // m/s
	double acc_y;
	double acc_z;
}INS;
typedef struct Par
{
	double Rm;
	double Rn;
	double w_ie[3];  //earth rotation rate in n-frame;
	double w_en[3];  //trans rate in n-frame
	double g[3];
	double f_n[3];   //specific force 
	double f_b[3];   //specific force 
	double w_b[3];   //B系下陀螺观测量
}Par;

typedef struct SensorOut
{
	double f_n[3];
	double f_b[3];
	double w_n[3];
	double w_b[3];
}SensorOut;
typedef struct Sensorbias
{
	double bias_acc_x;
	double bias_acc_y;
	double bias_acc_z;
	double bias_gyro_x;
	double bias_gyro_y;
	double bias_gyro_z;
	double Odo_scale;
}Sensorbias;
typedef struct Nav
{
	double lat;
	double lon;
	double height;
	double vn;
	double ve;
	double vd;
	double roll;
	double pitch;
	double heading;

	double dv_n[3];
	double wnb_n[3];   //载体角速度，在n系下b相对于n系的角速度，rad/s
	double a_n[3];
	double c_bn[3][3];
	double q_bn[4];
	double q_ne[4];
	double att[3];

	SensorOut mSensorOut;
	Sensorbias sensorbias;
	double Odo_scale;
}Nav;
enum MeasUpdataType
{
	MeasNone = 0,
	GNSSPoitionUpdate = 1 << 0,
	GNSSVelocityUpdate = 1 << 1,
	ZuptUpdate = 1 << 2,
	NHCUpdate = 1 << 3,
	OdoUpdate = 1 << 4,
	GNSSYawUpdate = 1 << 5,
};
typedef struct PVA
{
	double longitude;
	double latitude;
	double altitude;

	double northVelocity;
	double eastVelocity;
	double downVelocity;

	double roll;
	double pitch;
	double yaw;

	double g_bias[3];
	double a_bias[3];
	double Odo_scale;

	double xAcc;
	double yAcc;
	double zAcc;

	double rollRate;
	double pitchRate;
	double yawRate;

	double wnb_n[3]; //载体角速度，在n系下b系相对于n系的角速度，rad/s
	double a_n[3]; //载体加速度，n系下b系相对于n系的加速度，m/s^2
}PVA;
enum LC_STATUS
{
	INSDATA_NREADY = 0,
	INSDATA_READY,
	ALIGNMENT_ING,
	ALIGNMENT_COMPLETE,
	INS_FUSING
};
typedef struct ImuSensor
{
	uint16_t IMUDataRate;

	double Gyr_bias_std[3];//deg/h
	double Acc_bias_std[3];//mGal
	double Gyr_scale_std[3];//ppm
	double Acc_scale_std[3];//ppm
	double Gyr_noise_arw;//deg/sqrt(h)
	double Acc_noise_vrw;//m/s/sqrt(h)
	double Gyr_bias_CorTime;//h
	double Acc_bias_CorTime;//h
	double Gyr_scale_CorTime;//h
	double Acc_scale_CorTime;//h

	double bg_model[3];
	double ba_model[3];

	double q_bg[3];
	double q_ba[3];
}ImuSensor;
typedef struct GnssInsLC
{
	uint16_t ProcType;
	double Leverarm[3];
	double OdoLeverarm[3];
	double INSRotationRvb[3];
	ImuSensor mIMUSensor;
	double OdoRate;
	uint8_t   IsCar;
	double CarLeverarm[3];
	uint8_t   IsUseZUPT;
	uint8_t   IsUseNHC;
	uint8_t   IsUseOdo;
	uint16_t   IMUDataRate;
	uint8_t   InitialAttitudeMode;//0:motion Alignment  1:Dual antenna
	uint8_t    GnssDataType;

	double g_bias[3];
	double g_scale[3];
	double a_bias[3];
	double a_scale[3];
	double Odo_scale;
}GnssInsLC;
typedef struct UpdataStruct
{
	double LCTime;
	PVA mPVA;
	double C_bn[3][3];
	double w_b[3];
	double X[16];
	double P[16][16];
	double PHI[16][16];
	double Q[16][16];
	int8_t IsUseZUPT;
}UpdataStruct;
typedef struct KalmanStruct
{
	int16_t n;
	double P[16][16];
	double X[16];
	double H[3][16];
	double Z[3];
	double R[3][3];
}KalmanStruct;
typedef struct GnssInsSystem
{
	int8_t flag;
	double prevIMUTime;
	ImuData  mImuData;
	ImuData  mpreImuData;

	INS  mInsData;
	INS  mpreInsData;

	//Sensorbias mSensorbias;
	Nav mNav;
	Nav mpreNav;

	Par mPar;

	GnssData premGnssData;
	GnssData mGnssData;
	OdoData mOdoData;
	double OdoV;

	PVA mCarPVA;   //目标位置
	PVA mPVA;
	PVA mprePVA;
	ImuSensor mIMUSensor;
	uint8_t  GNSSflag;
	int64_t Sow;  //整秒GNSS时间

	int16_t IsUseNHC;
	int16_t IsUseOdo;

	int16_t IsUseNHCflag;
	double lastNHCtime;


	enum LC_STATUS mlc_STATUS;
	//double x[16];
	//double prex[16];
	//double phi[16][16];
	 double Q[16];     //system Q
	//double P[16][16];
	//double preP[16][16];
	KalmanStruct mKalmanStruct;

	double bg_model[3];
	double ba_model[3];
	double AlignCTime;

	int8_t InitNavi ;
	int8_t InsIsTimeOut ;
	int8_t KFIsConverage ;
	double TimeOutZuptTime ;
	double NaviStartTime ;
	double TimeOutStartTime;
	double LastBdsRejectedTime ;
	double LastBdsUpdateTime ;
	double LastInsUpdateTime ;
	int8_t GyroBiasIsConverage;
	double InsZuptTimeLen;

	int8_t NHCIsTimeOut;  //标志车载约束等待时间是否结束，暂定0.5s  ， -1表示正常没有等待到时，-2表示时间已经使用

	PVA lastPVA;

	int8_t CurIsUseZupt;
	enum MeasUpdataType mMeasUpdataType;
	int16_t PerMeasUpdata;
	int8_t IsIMUBiasFeedBack;

	double lastGNSSLCTime ;
	double lastGNSSUseTime ;
	double firstGNSSUseTime ;

	//0~DELY NEXT 1 LAST 0	
	double nextLCTime;	
	double lastLCTime;


	double lastGNSSLCtime;
	double lastGNSSLCUseTime;

	double lastNHCLCTime;
	double lastNHCUseTime;
	double nextNHCLCTime;

	double lastZuptTime ;

	//安装角
	double C_InstallationAngle[3][3]; // Rotation matrix from  vehicle frame to IMU frame C_vb

}GnssInsSystem;

/*-------------------------------------------------------------*/
/* INS KF */
int8_t KF_feedback(Par* mPar, double x[16], Nav* mNav, int IsBiasFeedBack, GnssInsSystem *mGnssInsSystem);

int8_t DataChangeLC(const double dt, const ImuData* mImuData, INS *cins);

int8_t compensate(const double dt, const Sensorbias* bias, INS* ins);



int8_t INS_MECH(const INS* ins_pre, const INS* ins_cur, const Nav* nav, Nav* nav1, Par* mPar);

/***************************************************************************************
Function: KF_predict_16PHI
Description: ;Calculate PHI(transition matrix)
Input :dt  time increment
	   nav  Navigation information
Output:PHI  transition matrix
Return :
Others:
********************************************************************************************/
int8_t KF_predict_16PHI(const double dt, const Nav *mNav, const Par *mPar, const ImuSensor *mImuSensor, double PHI[16][16]);
/*******************************************************************************************
Function: KF_predict
Description: ;
Input :x    state vect,old
	   P    Covariance,old
	   PHI  transition matrix
	   Q    the continuous-time system noise
	   dt   time increment
Output:x1  state vect new
	   P1  Covariance new
Return :
Others:
*************************************************************************************************/
int8_t  KF_predict(const double PHI[16][16], const double Q[16], const double dt, double x[16], double P[16][16], UpdataStruct* mUpdataStruct);



/*-------------------------------------------------------------*/
/* orientation related functions */
	//euler angles to direction cosine matrix  
void euler2dcm(const double eular[3], double dc[3][3]);
	//euler angles to quaternions  
void  euler2quat(const double eular[3], double q[4]);
	//direction cosine matrix to euler angles
void dcm2euler(const double dc[3][3], double eular[3]);
	//direction cosine matrix to quaternions
void dcm2quat(const double dc[3][3], double q[4]);
	//quaternions to euler angles
void quat2euler();
	//quaternions to direction cosine matrix 
void quat2dcm(const double q[4], double dc[3][3]);

	//referance frames and transformation
	//Compute the DCM (C_ne) representing the attitude of the navigation frame from the latitude and longitude.
void pos2dcm(double lat, double lon, double C_ne[3][3]);
	//Compute the quaternion (q_ne) representing the attitude of the navigation frame.
void pos2quat(double lat, double lon, double q[4]);
	// Position error to rotation vector conversion
	// void dpos2rvec(double delta[3], double rv[3])
	//{
	//	rv[0] = delta[2] * cos(delta[0]);
	//	rv[1] = -delta[1];
	//	rv[2] = -delta[2] * sin(delta[0]);
	//}
	// Position error to rotation vector conversion
void dpos2rvec(double lat, double delta_lat, double delta_lon, double rv[3]);


/*-------------------------------------------------------------*/

typedef struct LCSetting
{
	int16_t ProcType;
	char InsfileName[255], GnssfileNme[255], OdofileName[255];
	int16_t mIMUSenorType;  //IMU type
	int16_t Mode ;       //0-MOTION 1-GIVEN  2 STATIC
	double Leverarm[3];  //Offset from the IMU center of navigation to the phase center of the primary GNSS antenna.
	double OdoLeverarm[3];//Offset from the IMU center of navigation to the  center of the Odometer.
	double INSRotationRvb[3];//Rotation from the vehicle frame to the IMU body frame.
	int8_t IsUser;
	double UserLeverarm[3];//Offset from the IMU center of navigation to the  center of the output user.
	int8_t IsOutputKML;

	uint8_t   GnssDataType;
	uint8_t   InitialAttitudeMode;

	uint8_t   IsUseZUPT;
	uint8_t   IsUseNHC;
	uint8_t   IsUseOdo;
	int16_t   IMUDataRate;
	int16_t   ODODataRate;

	//SystemInit by Given attitue
	double Attitue_RPH[3]; // deg
	//init imu sensor zero offset
	double a_bias[3];   //mGal
	double g_bias[3];   //deg/h
	double a_scale[3];   //ppm
	double g_scale[3];   //ppm
	double Odo_scale;

}LCSetting;


int8_t ADDGNSSDATA(const GnssData msg);
int8_t AddIMUData(const ImuData mImudata, FILE* fSOL, FILE* fLOG);

int8_t readconfigfromfile(const char* cfgfname, LCSetting*  mLCSetting);

int8_t ouputfile(GnssInsSystem GnssInsSystem, FILE *output, FILE *output_debug);


int8_t initsystemfromcfg(GnssInsSystem* pGnssInsSystem);

int8_t initsystemfromGNSS(GnssInsSystem* pGnssInsSystem);

int8_t initsystemSoftreset();

int8_t SetLCtime(double lctime);

/* zupt */
int8_t SetZuptDetectData(const ImuData mImudata, int16_t imudataRata);
int8_t GetZuptVal();

int8_t process_adi16485_imu(const char* fname);

#ifdef __cplusplus
}
#endif

#endif