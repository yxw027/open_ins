#ifndef _INS_OUTMESGS_
#define _INS_OUTMESGS_
#include <stdint.h>
#include <stdio.h>
#include<string.h>
#include "datatype.h"
#include "lcstruct.h"

/*
typedef enum SolutionStatus
{
    SOL_COMPUTED,           //!< solution computed
    INSUFFICIENT_OBS,       //!< insufficient observations
    NO_CONVERGENCE,         //!< noconvergence
    SINGULARITY,            //!< singularity at parameters matrix
    COV_TRACE,              //!< covariance trace exceeds maximum (trace>1000m)
    TEST_DIST,              //!< test distance exceeded (max of 3 rejections if distance > 10km)
    COLD_START,             //!< not yet converged from cold start
    V_H_LIMIT,              //!< height or velocity limits exceeded
    VARIANCE,               //!< variance exceeds limits
    RESIDUALS,              //!< residuals are too large
    DELTA_POS,              //!< delta position is too large
    NEGATIVE_VAR,           //!< negative variance
    INTEGRITY_WARNING = 13, //!< large residuals make position unreliable
    INS_INACTIVE,           //!< ins has not started yet
    INS_ALIGNING,           //!< ins doing its coarse alignment
    INS_BAD,                //!< ins position is bad
    IMU_UNPLUGGED,          //!< no imu detected
    PENDING,                //!< when a fix position command is entered, the receiver computes its own position and determines if the fixed position is valid
    INVALID_FIX,            //!< the fixed position entered using the fix position command is not valid
    UNAUTHORIZED
} SolutionStatus;*/
/*
typedef enum PositionType
{
    NONE = 0,
    FIXEDPOS = 1,
    FIXEDHEIGHT = 2,
    Reserved = 3,
    FLOATCONV = 4,
    WIDELANE = 5,
    NARROWLANE = 6,
    DOPPLER_VELOCITY = 8,
    SINGLE = 16,
    PSRDIFF = 17,
    WAAS = 18,
    PROPOGATED = 19,
    OMNISTAR = 20,
    L1_FLOAT = 32,
    IONOFREE_FLOAT = 33,
    NARROW_FLOAT = 34,
    L1_INT = 48,
    WIDE_INT = 49,
    NARROW_INT = 50,
    RTK_DIRECT_INS = 51,
    INS1 = 52,
    INS_PSRSP = 53,
    INS_PSRDIFF = 54,
    INS_RTKFLOAT = 55,
    INS_RTKFIXED = 56,
    OMNISTAR_HP = 64,
    OMNISTAR_XP = 65,
    CDGPS = 66,
} PositionType;*/
/*
typedef enum DatumID
{
    ADIND = 1,
    ARC50,
    ARC60,
    AGD66,
    AGD84,
    BUKIT,
    ASTRO,
    CHATM,
    CARTH,
    CAPE,
    DJAKA,
    EGYPT,
    ED50,
    ED79,
    GUNSG,
    GEO49,
    GRB36,
    GUAM,
    HAWAII,
    KAUAI,
    MAUI,
    OAHU,
    HERAT,
    HJORS,
    HONGK,
    HUTZU,
    INDIA,
    IRE65,
    KERTA,
    KANDA,
    LIBER,
    LUZON,
    MINDA,
    MERCH,
    NAHR,
    NAD83,
    CANADA,
    ALASKA,
    NAD27,
    CARIBB,
    MEXICO,
    CAMER,
    MINNA,
    OMAN,
    PUERTO,
    QORNO,
    ROME,
    CHUA,
    SAM56,
    SAM69,
    CAMPO,
    SACOR,
    YACAR,
    TANAN,
    TIMBA,
    TOKYO,
    TRIST,
    VITI,
    WAK60,
    WGS72,
    WGS84,
    ZANDE,
    USER,
    CSRS,
    ADIM,
    ARSM,
    ENW,
    HTN,
    INDB,
    INDI,
    IRL,
    LUZA,
    LUZB,
    NAHC,
    NASP,
    OGBM,
    OHAA,
    OHAB,
    OHAC,
    OHAD,
    OHIA,
    OHIB,
    OHIC,
    OHID,
    TIL,
    TOYM
} DatumID;
*/
typedef enum MessageFormat //!< Bits 5-6 of MessageType struct
{
    BINARY = 0b00,
    ASCII = 0b01,
    ABREVIATED_ASCII = 0b10,
    NMEA = 0b11,
} MessageFormat;
typedef enum ResponseBit //!< Last bit (7) of MessageType struct
{
    ORIGINAL_MESSAGE = 0b0,
    RESPONSE_MESSAGE = 0b1,
} ResponseBit;
#pragma pack(1)
typedef struct MessageType
{
    unsigned reserved : 5;
    MessageFormat format : 2;
    ResponseBit response : 1;
}MessageType;
// __attribute__ ((__aligned__(1)))
typedef struct Oem4BinaryHeader
{
    uint8_t sync1;            //!< start of packet first byte (0xAA)
    uint8_t sync2;            //!< start of packet second byte (0x44)
    uint8_t sync3;            //!< start of packet third  byte (0x12)
    uint8_t header_length;    //!< Length of the header in bytes ( From start of packet )
    uint16_t message_id;      //!< Message ID number
    MessageType message_type; //!< Message type - binary, ascii, nmea, etc...
    uint8_t port_address;     //!< Address of the data port the log was received on
    uint16_t message_length;  //!< Message length (Not including header or CRC)
    uint16_t sequence;        //!< Counts down from N-1 to 0 for multiple related logs
    uint8_t idle;             //!< Time the processor was idle in last sec between logs with same ID
    uint8_t time_status;      //!< Indicates the quality of the GPS time
    uint16_t gps_week;        //!< GPS Week number
    uint32_t gps_millisecs;   //!< Milliseconds into week
    uint32_t status;          //!< Receiver status word
    uint16_t Reserved;        //!< Reserved for internal use
    uint16_t version;         //!< Receiver software build number (0-65535)
} Oem4BinaryHeader;
typedef struct Position
{
    Oem4BinaryHeader header;                  //!< Message header
    uint32_t solution_status;                 //!< Solution status
    uint32_t position_type;                   //!< Position type
    double latitude;                          //!< latitude (deg)
    double longitude;                         //!< longitude (deg)
    double height;                            //!< height above mean sea level (m)
    float undulation;                         //!< Undulation - the relationship between the geoid and the ellipsoid (m)
    uint32_t datum_id;                        //!< datum id number
    float latitude_standard_deviation;        //!< latitude standard deviation (m)
    float longitude_standard_deviation;       //!< longitude standard deviation (m)
    float height_standard_deviation;          //!< height standard deviation (m)
    int8_t base_station_id[4];                //!< base station id
    float differential_age;                   //!< differential position age (sec)
    float solution_age;                       //!< solution age (sec)
    uint8_t number_of_satellites;             //!< number of satellites tracked
    uint8_t number_of_satellites_in_solution; //!< number of satellites used in solution
    uint8_t num_gps_plus_glonass_l1;          //!< number of GPS plus GLONASS L1 satellites used in solution
    uint8_t num_gps_plus_glonass_l2;          //!< number of GPS plus GLONASS L2 satellites used in solution
    uint8_t reserved;                         //!< reserved
    uint8_t extended_solution_status;         //!< extended solution status - OEMV and greater only
    uint8_t reserved2;                        //!< reserved
    uint8_t signals_used_mask;                //!< signals used mask - OEMV and greater only
    uint8_t crc[4];                           //!< 32-bit cyclic redundancy check (CRC)
} Position;
/*!
 * Velocity Message Structure
 * This log contains the best available velocity
 * information computed by the receiver. In addition,
 * it reports a velocity status indicator, which is
 * useful in indicating whether or not the corresponding
 * data is valid. The velocity measurements sometimes
 * have a latency associated with them. The time of validity
 * is the time tag in the log minus the latency value.
 *
 * This structure represents the format of the following messages:
 *  - BESTVEL
 *  - RTKVEL
 *  - PSRVEL
 */
typedef struct Velocity
{
	Oem4BinaryHeader header;			//!< Message header
	uint32_t solution_status;		//!< Solution status
	uint32_t position_type;			//!< Position type
	float latency;						//!< measure of the latency of the velocity time tag in seconds
	float age;							//!< differential age in seconds
	double horizontal_speed;			//!< horizontal speed in m/s
	double track_over_ground;			//!< direction of travel in degrees
	double vertical_speed; 				//!< vertical speed in m/s
	float reserved;
	int8_t crc[4];
}Velocity;
/*!
 * INSPVA Message Structure
 * This log allows INS position, velocity and
 * attitude to be collected in one log, instead
 * of using three separate logs.
 */
 typedef struct InsPositionVelocityAttitude
{
	Oem4BinaryHeader header;	//!< Message header
	uint32_t gps_week;			//!< GPS week number
	double gps_millisecs;		//!< Milliseconds into GPS week
	double latitude;			//!< latitude - WGS84 (deg)
	double longitude;			//!< longitude - WGS84 (deg)
	double height;				//!< Ellipsoidal height - WGS84 (m)
	double north_velocity;		//!< velocity in a northerly direction (m/s)
	double east_velocity;		//!< velocity in an easterly direction (m/s)
	double up_velocity;			//!< velocity in an up direction
	double roll;				//!< right handed rotation around y-axis (degrees)
	double pitch;				//!< right handed rotation aruond x-axis (degrees)
	double azimuth;				//!< right handed rotation around z-axis (degrees)
	int32_t status;			//!< status of the INS system
	int8_t crc[4];
}InsPositionVelocityAttitude;
typedef struct INSPVAX
{
	Oem4BinaryHeader header;	//!< Message header
	int32_t ins_status;         //!< Solution status
	int32_t pos_type;           //!< Position type
	double latitude;			//!< latitude - WGS84 (deg)
	double longitude;			//!< longitude - WGS84 (deg)
	double height;				//!< Height above mean sea level  - WGS84 (m)
	float undulation;          //!< Undulation (m)
	double north_velocity;		//!< velocity in a northerly direction (m/s)
	double east_velocity;		//!< velocity in an easterly direction (m/s)
	double up_velocity;			//!< velocity in an up direction
	double roll;				//!< right handed rotation around y-axis (degrees)
	double pitch;				//!< right handed rotation aruond x-axis (degrees)
	double azimuth;				//!< right handed rotation around z-axis (degrees)
	float latitude_std;
	float longitude_std;
	float altitude_std;
	float north_velocity_std;
	float east_velocity_std;
	float up_velocity_std;
	float roll_std;
	float pitch_std;
	float azimuth_std;
	int32_t Ext_sol_stat;			//!< Extended solution status
	int16_t time_since_update;      //!< Elapsed time since the last ZUPT or positionupdate (seconds)
	int8_t crc[4];
}INSPVAX;
typedef struct ImuStatus
{
    unsigned counter : 4;                    //!< 4 byte counter
    unsigned imu_test : 1;                   //!< IMU test: Passed=0, Failed=1
    unsigned z_axis_path_length_control : 1; //!< Z-axis path length control: Good=0, Reset=1
    unsigned y_axis_path_length_control : 1; //!< Y-axis path length control: Good=0, Reset=1
    unsigned x_axis_path_length_control : 1; //!< X-axis path length control: Good=0, Reset=1
    unsigned accelerometer_temperature : 8;  //!< Accelerometer temperature
    unsigned software_version : 8;           //!< IMU software version number
    unsigned reserved : 3;                   //!< Reserved
    unsigned gyro_test : 1;                  //!< Gyro tests: Passed=0, Failed=1
    unsigned accel_test : 1;                 //!< Accelerometer tests: Passed=0, Failed=1
    unsigned other_tests : 1;                //!< Other tests: Passed=0, Failed=1
    unsigned memory_tests : 1;               //!< Memory tests: Passed=0, Failed=1
    unsigned processor_tests : 1;            //!< Processor tests: Passed=0, Failed=1
} ImuStatus;
typedef struct RawImu
{
    Oem4BinaryHeader header;  //!< Message header
    uint32_t gps_week;        //!< GPS week number
    double gps_millisecs;     //!< Milliseconds into GPS week
    ImuStatus imuStatus;      //!< Status of the IMU
    float z_acceleration;     //!< change in velocity along z axis in scaled m/s
    float y_acceleration_neg; //!< -change in velocity along y axis in scaled m/s
    float x_acceleration;     //!< change in velocity along x axis in scaled m/s
    float z_gyro_rate;        //!< change in angle around z axis in radians
    float y_gyro_rate_neg;    //!< -(change in angle around y axis) in radians
    float x_gyro_rate;        //!< change in angle around x axis in radians
    int8_t crc[4];
} RawImu;
typedef struct OdoSpeed
{
    Oem4BinaryHeader header;  //!< Message header
    int week;
    double gps_millisecs;     //!< Milliseconds into GPS week
    char mode;
    double speed;
    char fwd;
    uint64_t wheel_tick;
    int8_t crc[4];
} OdoSpeed;
#pragma pack()



unsigned long CalculateBlockCRC32(unsigned long ulCount, unsigned char *ucBuffer); /* Number of bytes in the data block  Data block */

int writeRawImuMsg(const int gps_update,const int32_t week1,const uint32_t sec1, const int32_t week2, const uint32_t sec2,const ImuData* p_ImuData);

int writePositionMsg(const int week, const uint32_t gps_millisecs,const GnssData* p_GnssData);

int writeVelocityMsg(const int week, const uint32_t gps_millisecs,const GnssData* p_GnssData);

int writeGnssRawMsg(int8_t* gnss_obs_flag, const int week, const uint32_t gps_millisecs, const GnssData* p_GnssData);

int writeOdoDataMsg(const OdoData* p_OdoData);

int writeGGAMsg(int week, double time, const GnssInsSystem* p_gnssInsSystem);

int writeINSPVAXMsg(int32_t week, uint32_t itow, const GnssInsSystem* p_gnssInsSystem);

int writeINSPVAMsg(int32_t week, uint32_t itow, const GnssInsSystem* p_gnssInsSystem);


extern RawImu imuStr;
extern Position positionStr ;
extern Velocity velocityStr;
extern InsPositionVelocityAttitude inspvastr;
extern INSPVAX inspvaxstr;
extern char ggaBuff[120];
extern char pashrBuff[120];
extern char rmcBuff[200];
extern char vtgBuff[120];

#endif // !_INS_OUTMESGS_
