#include "gnss_ins.h"
#include <math.h>
#include <string.h>
#include <stdlib.h> 

static UpdataStruct mUpdataStruct;
static GnssInsSystem mGnssInsSystem;
static int8_t IsFinished;

static LCSetting mLCSetting;

/* coordinate computations */
const GravityParameter normg = { 9.7803267715 ,0.0052790414 ,0.0000232718,-0.000003087691089,0.00000000439773 , 0.000000000000721 };
const  EarthParameter WGS84 = { 6378137.0, 6356752.3142, 0.0033528106643315515,0.081819190837555025,0.0066943799893122479 ,  7.2922115147e-5,398600441800000.00 };

uint8_t UpdateMN(const double *BLH, double *M, double *N)
{
	double sinB = sin(*BLH);
	double temp = 1 - WGS84.e2 * sinB * sinB;
	double sqrttemp = sqrt(temp);
	*M = WGS84.a * (1 - WGS84.e2) / (sqrttemp*temp);
	*N = WGS84.a / sqrttemp;
	return 1;
};
uint8_t UpdateW(const double *BLH, const double M, const double N, const double *vn, double *wnie, double *wnen)
{
	double cosB = cos(*BLH);
	double sinB = sin(*BLH);
	double tanB = tan(*BLH);
	wnie[0] = WGS84.wie * cosB;
	wnie[1] = 0.0;
	wnie[2] = -WGS84.wie * sinB;
	wnen[0] = vn[1] / (N + BLH[2]);
	wnen[1] = -vn[0] / (M + BLH[2]);
	wnen[2] = -vn[1] * tanB / (N + BLH[2]);
	return 1;
};

uint8_t UpdateGravity(const double *BLH, double* ng)
{
	double sinB = sin(*BLH);
	double s2 = sinB * sinB;
	double s4 = s2 * s2;
	*ng = normg.g0 * (1.0f + normg.g1 * s2 + normg.g2 * s4) + (normg.g3 + normg.g4 * s2) * BLH[2] + normg.g5 * BLH[2] * BLH[2];
	return 1;
}

int8_t llh2ecef(const Geo* mGeo, ECEF* mECEF)
{
	double M, N;
	UpdateMN(&(mGeo->lat), &M, &N);
	double c_lat, s_lat, c_lon, s_lon, Rn_h;
	c_lat = cos(mGeo->lat);
	s_lat = sin(mGeo->lat);
	c_lon = cos(mGeo->lon);
	s_lon = sin(mGeo->lon);
	Rn_h = N + mGeo->height;
	mECEF->x = Rn_h * c_lat*c_lon;
	mECEF->y = Rn_h * c_lat*s_lon;
	mECEF->z = (N*(1 - WGS84.e2) + mGeo->height)*s_lat;
	return 1;
}
int8_t ecef2llh(const ECEF* mECEF, Geo* mGeo)
{
	double M, N;
	mGeo->lon = atan2(mECEF->y, mECEF->x);
	double p = sqrt(mECEF->x*mECEF->x + mECEF->y * mECEF->y);
	mGeo->lat = 0.0f;
	UpdateMN(&(mGeo->lat), &M, &N);
	mGeo->height = p / cos(mGeo->lat) - N;
	mGeo->lat = atan2(mECEF->z, p * (1 - WGS84.e2 * N / (N + mGeo->height)));
	for (int i = 0; i < 4; i++)
	{
		UpdateMN(&(mGeo->lat), &M, &N);
		mGeo->height = p / cos(mGeo->lat) - N;
		mGeo->lat = atan2(mECEF->x, p * (1 - WGS84.e2 * N / (N + mGeo->height)));
	}
	return 1;
}
/* matrix functions */


	uint8_t MatrixAdd(const double *matrix_a, const double *matrix_b, const int matrix_a_row, const int matrix_a_colume, double *matrix_result)
	{
		for (int i = 0; i < matrix_a_row*matrix_a_colume; i++)
		{
			*(matrix_result + i) = *(matrix_a + i) + *(matrix_b + i);
		}
		return 1;
	}
	uint8_t MatrixSub(const double *matrix_a, const double *matrix_b, const int matrix_a_row, const int matrix_a_colume, double *matrix_result)
	{
		for (int i = 0; i < matrix_a_row*matrix_a_colume; i++)
		{
			*(matrix_result + i) = *(matrix_a + i) - *(matrix_b + i);
		}
		return 1;

	}
	uint8_t MatrixMutiply(const double *matrix_a, const double *matrix_b, const int matrix_a_row, const int matrix_a_column, const int matrix_b_column, double *matrix_result)
	{
		double sum = 0;
		double median = 0;
		for (int i = 0; i < matrix_a_row; i++)
		{
			for (int k = 0; k < matrix_b_column; k++)
			{
				for (int j = 0; j < matrix_a_column; j++)
				{
					median = matrix_a[matrix_a_column*i + j] * matrix_b[matrix_b_column*j + k];
					sum = sum + median;
				}
				matrix_result[matrix_b_column*i + k] = sum;
				sum = 0;
			}
		}
		return 1;

	}
	uint8_t MatrixTranspose(const double *matrix_a, const int matrix_a_row, const int matrix_a_column, double *matrix_result)
	{
		for (int i = 0; i < matrix_a_column; i++)
		{
			for (int j = 0; j < matrix_a_row; j++)
			{
				matrix_result[matrix_a_row*i + j] = matrix_a[matrix_a_column*j + i];
			}
		}
		return 1;

	}
	/*************************************************
	Function: MatrixCholosky
	Description:
	Input:
	Output:lower_tri_matrix
	Return:
	Others:
	*************************************************/
	uint8_t MatrixCholosky(const double *matrix, const int matrix_row, double *lower_tri_matrix)
	{
		int i, j, k, l, u;
		double sum = 0.0;
		if (matrix[0] <= 0.0)
		{
#ifdef DEBUG
			printf(" The matrix can't be decomposed with Cholosky decomposition method.\n");
#endif
		}
		lower_tri_matrix[0] = sqrt(matrix[0]);

		for (i = 1; i < matrix_row; i++)
			lower_tri_matrix[i*matrix_row] = matrix[i*matrix_row] / lower_tri_matrix[0];


		for (k = 1; k < matrix_row; k++)
		{
			l = k * matrix_row + k;
			for (j = 0; j < k; j++)
				sum = sum + lower_tri_matrix[k*matrix_row + j] * lower_tri_matrix[k*matrix_row + j];

			if ((lower_tri_matrix[l] - sum) <= 0.0)
			{
#ifdef DEBUG
				printf(" The matrix can't be decomposed with Cholosky decomposition method.\n");
#endif
			}
			lower_tri_matrix[l] = sqrt(matrix[l] - sum);
			sum = 0.0;
			for (i = k + 1; i < matrix_row; i++)
			{
				u = i * matrix_row + k;
				for (j = 0; j < k; j++)
					sum = sum + lower_tri_matrix[i*matrix_row + j] * lower_tri_matrix[k*matrix_row + j];
				lower_tri_matrix[u] = (matrix[u] - sum) / lower_tri_matrix[l];
				sum = 0.0;
			}
		}

		for (i = 0; i < matrix_row - 1; i++)
		{
			for (j = i + 1; j < matrix_row; j++)
				lower_tri_matrix[i*matrix_row + j] = 0.0;
		}
		return 1;
	}
	/*************************************************
	Function: MatrixInverse
	Description:
	Input:
	Output:matrix_result
	Return:
	Others:
	*************************************************/
//	uint8_t MatrixInverse(const double *lower_tri_matrix, const int matrix_row, double *matrix_result)
//	{
//		int ret = 0;
//		int i, j, k, l, u;
//		double sum;
//		for (i = 0; i < matrix_row; i++)
//		{
//			l = i * matrix_row + i;
//			if (lower_tri_matrix[l] <= 0)
//			{
//#ifdef DEBUG
//				printf(" The matrix does not exist inverse matrix.\n");
//				return -1;
//#endif
//			}
//			matrix_result[l] = 1.0 / lower_tri_matrix[l];
//		}
//
//		for (i = 1; i < matrix_row; i++)
//		{
//			sum = 0.0;
//			for (j = 0; j < i; j++)
//			{
//				for (k = j; k < i; k++)
//				{
//					l = i * matrix_row + k;
//					u = k * matrix_row + j;
//					sum = sum + lower_tri_matrix[l] * matrix_result[u];
//				}
//				matrix_result[i*matrix_row + j] = -matrix_result[i*matrix_row + i] * sum;
//				sum = 0.0;
//			}
//		}
//
//		for (i = 0; i < matrix_row - 1; i++)
//		{
//			for (j = i + 1; j < matrix_row; j++)
//			{
//				l = i * matrix_row + j;
//				matrix_result[l] = 0.0;
//			}
//		}
//		return 1;
//	}
	static unsigned short l, u;
	static char is[25], js[25];
	uint8_t MatrixInverse(uint8_t n, double* a)
	{
		int i, j, k;
		unsigned short v;
		double d, p;
		for (k = 0; k <= n - 1; k++)
		{
			d = 0.0;
			for (i = k; i <= n - 1; i++)
				for (j = k; j <= n - 1; j++)
				{
					l = i * n + j; p = fabs(a[l]);
					if (p > d) { d = p; is[k] = i; js[k] = j; }
				}
			if (is[k] != k)
				for (j = 0; j <= n - 1; j++)
				{
					u = k * n + j; v = is[k] * n + j;
					p = a[u]; a[u] = a[v]; a[v] = p;
				}
			if (js[k] != k)
				for (i = 0; i <= n - 1; i++)
				{
					u = i * n + k; v = i * n + js[k];
					p = a[u]; a[u] = a[v]; a[v] = p;
				}
			l = k * n + k;
			a[l] = 1.0 / a[l];
			for (j = 0; j <= n - 1; j++)
				if (j != k)
				{
					u = k * n + j; a[u] = a[u] * a[l];
				}
			for (i = 0; i <= n - 1; i++)
				if (i != k)
					for (j = 0; j <= n - 1; j++)
						if (j != k)
						{
							u = i * n + j;
							a[u] = a[u] - a[i*n + k] * a[k*n + j];
						}
			for (i = 0; i <= n - 1; i++)
				if (i != k)
				{
					u = i * n + k; a[u] = -a[u] * a[l];
				}
		}
		for (k = n - 1; k >= 0; k--)
		{
			if (js[k] != k)
				for (j = 0; j <= n - 1; j++)
				{
					u = k * n + j; v = js[k] * n + j;
					p = a[u]; a[u] = a[v]; a[v] = p;
				}
			if (is[k] != k)
				for (i = 0; i <= n - 1; i++)
				{
					u = i * n + k; v = i * n + is[k];
					p = a[u]; a[u] = a[v]; a[v] = p;
				}
		}
		return 1;

	}
	// Cross product of two vectors(c = a x b)
	uint8_t CrossProduct(double a[3], double b[3], double c[3])
	{
		c[0] = a[1] * b[2] - b[1] * a[2];
		c[1] = b[0] * a[2] - a[0] * b[2];
		c[2] = a[0] * b[1] - b[0] * a[1];
		return 1;

	}
	uint8_t GetSkewSymmetricMatrixOfVector(const double  mvector[3], double *M)
	{
		M[0] = 0; M[1] = -mvector[2]; M[2] = mvector[1];
		M[3] = mvector[2]; M[4] = 0; M[5] = -mvector[0];
		M[6] = -mvector[1]; M[7] = mvector[0]; M[8] = 0;
		return 1;

	}
	uint8_t quatprod(const double p[4], const double q[4], double r[4])
	{
		r[0] = p[0] * q[0] - p[1] * q[1] - p[2] * q[2] - p[3] * q[3];
		r[1] = p[1] * q[0] + p[0] * q[1] - p[3] * q[2] + p[2] * q[3];
		r[2] = p[2] * q[0] + p[3] * q[1] + p[0] * q[2] - p[1] * q[3];
		r[3] = p[3] * q[0] - p[2] * q[1] + p[1] * q[2] + p[0] * q[3];
		return 1;

	}
	uint8_t quat2rvec(double q[4], double rot_vec[3])
	{
		int i;
		if (fabs(q[0]) > 1.0E-15)
		{
			double atan_05zata = atan2(sqrt(q[1] * q[1] + q[2] * q[2] + q[3] * q[3]), q[0]);
			double atan2_05zata = atan_05zata * atan_05zata;
			double atan4_05zata = atan2_05zata * atan2_05zata;
			double atan6_05zata = atan4_05zata * atan2_05zata;
			double f = 0.5*(1 - atan2_05zata / 6 + atan4_05zata / 120 - atan6_05zata / 5040);
			for (i = 0; i < 3; ++i)
			{
				rot_vec[i] = q[i + 1] / f;
			}
		}
		else
		{
			for (i = 0; i < 3; ++i)
			{
				rot_vec[i] = q[i + 1] * PI;
			}
		}
		return 1;

	}

	uint8_t rvec2quat(double rot_vec[3], double q[4])
	{
		double mag2, c, s;
		mag2 = rot_vec[0] * rot_vec[0] + rot_vec[1] * rot_vec[1] + rot_vec[2] * rot_vec[2];

		if (mag2 < PI*PI)
		{
			mag2 = 0.25*mag2;

			c = 1.0 - mag2 / 2.0*(1.0 - mag2 / 12.0*(1.0 - mag2 / 30.0));
			s = 1.0 - mag2 / 6.0*(1.0 - mag2 / 20.0*(1.0 - mag2 / 42.0));

			if (c < 0)
			{
				q[0] = -c;
				q[1] = -0.5*s*rot_vec[0];
				q[2] = -0.5*s*rot_vec[1];
				q[3] = -0.5*s*rot_vec[2];
			}
			else
			{
				q[0] = c;
				q[1] = 0.5*s*rot_vec[0];
				q[2] = 0.5*s*rot_vec[1];
				q[3] = 0.5*s*rot_vec[2];
			}
		}
		else
		{
			c = sqrt(mag2);
			s = sin(c / 2);
			mag2 = s / c;

			q[0] = cos(c / 2);
			q[1] = rot_vec[0] * mag2;
			q[2] = rot_vec[1] * mag2;
			q[3] = rot_vec[2] * mag2;
			if (q[0] < 0)
			{
				q[0] = -q[0];
				q[1] = -q[1];
				q[2] = -q[2];
				q[3] = -q[3];
			}
		}
		return 1;

	}
	uint8_t quat2pos(double q[4], double *r)
	{
		r[0] = -2 * atan(q[2] / q[0]) - PI / 2;
		r[1] = 2 * atan2(q[3], q[0]);
		return 1;

	}
	uint8_t norm_quat(double q[4])
	{
		double e = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3] - 1) / 2;

		q[0] = (1 - e)*q[0];
		q[1] = (1 - e)*q[1];
		q[2] = (1 - e)*q[2];
		q[3] = (1 - e)*q[3];
		return 1;

	}
	uint8_t quatinv(const double p[4], double q[4])
	{
		double a = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2] + p[3] * p[3]);
		if (a < 1.0E-15)
		{
			q[0] = 1.0; q[1] = 0.0; q[2] = 0.0; q[3] = 0.0;
			return -1;
		}
		q[0] = p[0] / a;
		q[1] = -p[1] / a;
		q[2] = -p[2] / a;
		q[3] = -p[3] / a;
		return 1;

	}



/* orientation functions */

//euler angles to direction cosine matrix  
void euler2dcm(const double eular[3], double dc[3][3])
{
	double roll = eular[0];
	double  pitch = eular[1];
	double  heading = eular[2];
	double  cr, cp, ch, sr, sp, sh;
	cr = cos(roll); cp = cos(pitch); ch = cos(heading);
	sr = sin(roll); sp = sin(pitch); sh = sin(heading);

	dc[0][0] = cp * ch;
	dc[0][1] = -cr * sh + sr * sp*ch;
	dc[0][2] = sr * sh + cr * sp*ch;

	dc[1][0] = cp * sh;
	dc[1][1] = cr * ch + sr * sp*sh;
	dc[1][2] = -sr * ch + cr * sp * sh;

	dc[2][0] = -sp;
	dc[2][1] = sr * cp;
	dc[2][2] = cr * cp;
}
//euler angles to quaternions  
void  euler2quat(const double eular[3], double q[4])
{
	double c_r = cos(eular[0] / 2);
	double c_p = cos(eular[1] / 2);
	double c_h = cos(eular[2] / 2);
	double s_r = sin(eular[0] / 2);
	double s_p = sin(eular[1] / 2);
	double s_h = sin(eular[2] / 2);
	q[0] = c_r * c_p * c_h + s_r * s_p * s_h;
	q[1] = s_r * c_p * c_h - c_r * s_p * s_h;
	q[2] = c_r * s_p * c_h + s_r * c_p * s_h;
	q[3] = c_r * c_p * s_h - s_r * s_p * c_h;
}
//direction cosine matrix to euler angles
void dcm2euler(const double dc[3][3], double eular[3])
{
	double roll, pitch, heading;
	pitch = atan(-dc[2][0] / sqrt(dc[2][1] * dc[2][1] + dc[2][2] * dc[2][2]));
	if (dc[2][0] <= -0.999)
	{
		roll = atan2(dc[2][1], dc[2][2]);
		heading = atan2((dc[1][2] - dc[0][1]), (dc[0][2] + dc[1][1]));
	}
	else if (dc[2][0] >= 0.999)
	{
		roll = atan2(dc[2][1], dc[2][2]);
		heading = PI + atan2((dc[1][2] + dc[0][1]), (dc[0][2] - dc[1][1]));
	}
	else
	{
		roll = atan2(dc[2][1], dc[2][2]);
		heading = atan2(dc[1][0], dc[0][0]);
	}
	if (heading < 0)// heading 0-360
	{
		heading = 2 * PI + heading;
	}
	eular[0] = roll;
	eular[1] = pitch;
	eular[2] = heading;
}
//direction cosine matrix to quaternions
void dcm2quat(const double dc[3][3], double q[4])
{

}
//quaternions to euler angles
void quat2euler()
{

}
//quaternions to direction cosine matrix 
void quat2dcm(const double q[4], double dc[3][3])
{
	dc[0][0] = q[0] * q[0] + q[1] * q[1] - q[2] * q[2] - q[3] * q[3];
	dc[0][1] = 2 * (q[1] * q[2] - q[0] * q[3]);
	dc[0][2] = 2 * (q[1] * q[3] + q[0] * q[2]);

	dc[1][0] = 2 * (q[1] * q[2] + q[0] * q[3]);
	dc[1][1] = q[0] * q[0] - q[1] * q[1] + q[2] * q[2] - q[3] * q[3];
	dc[1][2] = 2 * (q[2] * q[3] - q[0] * q[1]);

	dc[2][0] = 2 * (q[1] * q[3] - q[0] * q[2]);
	dc[2][1] = 2 * (q[2] * q[3] + q[0] * q[1]);
	dc[2][2] = q[0] * q[0] - q[1] * q[1] - q[2] * q[2] + q[3] * q[3];
}

//referance frames and transformation
//Compute the DCM (C_ne) representing the attitude of the navigation frame from the latitude and longitude.
void pos2dcm(double lat, double lon, double C_ne[3][3])
{
	double sinB = sin(lat);
	double sinL = sin(lon);
	double cosB = cos(lat);
	double cosL = cos(lon);

	C_ne[0][0] = -sinB * cosL;
	C_ne[0][1] = -sinL;
	C_ne[0][2] = -cosB * cosL;

	C_ne[1][0] = -sinB * sinL;
	C_ne[1][1] = cosL;
	C_ne[1][2] = -cosB * sinL;

	C_ne[2][0] = cosB;
	C_ne[2][1] = 0;
	C_ne[2][2] = -sinB;
}
//Compute the quaternion (q_ne) representing the attitude of the navigation frame.
void pos2quat(double lat, double lon, double q[4])
{
	double cos_L = cos(lon / 2);
	double sin_L = sin(lon / 2);
	double cos_B = cos(-PI / 4 - lat / 2);
	double sin_B = sin(-PI / 4 - lat / 2);
	q[0] = cos_B * cos_L;
	q[1] = -sin_B * sin_L;
	q[2] = sin_B * cos_L;
	q[3] = cos_B * sin_L;
}
// Position error to rotation vector conversion
// void dpos2rvec(double delta[3], double rv[3])
//{
//	rv[0] = delta[2] * cos(delta[0]);
//	rv[1] = -delta[1];
//	rv[2] = -delta[2] * sin(delta[0]);
//}
// Position error to rotation vector conversion
void dpos2rvec(double lat, double delta_lat, double delta_lon, double rv[3])
{
	rv[0] = delta_lon * cos(lat);
	rv[1] = -delta_lat;
	rv[2] = -delta_lon * sin(lat);
}

/*--------------------------------------------------------*/
/* ZUPT */
static double ZuptdetectData[50][6];
static int16_t DateLength = 0;
static double lastZuptDetTime = -999.99;

static int16_t ZuptValue1 = 0;
static int16_t lastZuptValue1 = 0;
static int16_t ZuptValue2 = 0;
int8_t SetZuptDetectData(const ImuData mImudata, int16_t imudataRata)
{
	int32_t CheckDataLegth = imudataRata;
	if (DateLength == CheckDataLegth)
	{
		for (int16_t i = 0; i < DateLength - 1; i++)
		{
			for (int16_t j = 0; j < 6; j++)
			{
				ZuptdetectData[i][j] = ZuptdetectData[i + 1][j];
			}
		}
		ZuptdetectData[DateLength - 1][0] = mImudata.Gyrox*imudataRata;
		ZuptdetectData[DateLength - 1][1] = mImudata.Gyroy*imudataRata;
		ZuptdetectData[DateLength - 1][2] = mImudata.Gyroz*imudataRata;
		ZuptdetectData[DateLength - 1][3] = mImudata.Accx*imudataRata;
		ZuptdetectData[DateLength - 1][4] = mImudata.Accy*imudataRata;
		ZuptdetectData[DateLength - 1][5] = mImudata.Accz*imudataRata;
	}
	if (DateLength < CheckDataLegth)
	{
		ZuptdetectData[DateLength][0] = mImudata.Gyrox*imudataRata;
		ZuptdetectData[DateLength][1] = mImudata.Gyroy*imudataRata;
		ZuptdetectData[DateLength][2] = mImudata.Gyroz*imudataRata;
		ZuptdetectData[DateLength][3] = mImudata.Accx*imudataRata;
		ZuptdetectData[DateLength][4] = mImudata.Accy*imudataRata;
		ZuptdetectData[DateLength][5] = mImudata.Accz*imudataRata;
		DateLength += 1;
		return 1;
	}

	double Max[6] = { 0.0 };
	double Mean[6] = { 0.0 };
	double Var[6] = { 0.0 };
	double Std[6] = { 0.0 };
	double dt = mImudata.timestamp - floor(mImudata.timestamp);
	double Dt = 1.0/ imudataRata;
	if ((dt < Dt || (dt > 0.5 - Dt && dt < 0.5 + Dt)
		&& mImudata.timestamp > lastZuptDetTime + 0.45))
	{
		ZuptValue1 = 0;
		for (int16_t i = 0; i < DateLength; i++)
		{
			for (int16_t j = 0; j < 6; j++)
			{
				Max[j] += ZuptdetectData[i][j];
			}
		}
		for (int16_t j = 0; j < 6; j++)
		{
			Mean[j] = Max[j] / DateLength;
		}
		for (int16_t i = 0; i < DateLength; i++)
		{
			for (int16_t j = 0; j < 6; j++)
			{
				Var[j] += (ZuptdetectData[i][j] - Mean[j])*(ZuptdetectData[i][j] - Mean[j]);
			}
		}
		for (int16_t j = 0; j < 6; j++)
		{
			Std[j] = sqrt(Var[j] / DateLength);
		}
		if (Std[0] < 0.0006 && Std[1] < 0.0006 && Std[2] < 0.0006
			&&Std[3] < 0.03 && Std[4] < 0.03)  //&& Std[5] < 0.6
		{
			ZuptValue1 = 1;
		}
		if (ZuptValue1 && lastZuptValue1)
		{
			ZuptValue2 += 1;
		}
		if (ZuptValue2 > 5)
		{
			ZuptValue2 = 5;
		}
		if (!ZuptValue1 && lastZuptValue1)
		{
			ZuptValue2 = (ZuptValue2 - 4 > 0) ? (ZuptValue2 - 4) : 0;
		}
		if (!ZuptValue1 && !lastZuptValue1)
		{
			ZuptValue2 = 0;
		}

		lastZuptValue1 = ZuptValue1;
		lastZuptDetTime = mImudata.timestamp;
		return 1;
		//DateLength = 0;
	}
	return 1;
}

int8_t GetZuptVal()
{
	int8_t ret = 0;
	if (ZuptValue2)
	{
		ret = 1;
	}
	return ret;
}

/* transition matrix */

#define SMD(i) ((i)*((i)+1)/2)
#define SMI(i, j)	((i) > (j) ? (SMD(i) + (j)) : (SMD(j) + (i)))
#define NX_INS 16          //15 + 1
#define IsOdo 1

//PHI*P
void PHI_P(const double* PHI, const double* P, double* PHIP)
{
	unsigned char i, j;
	const double* phi;
	double* c;
	memset(PHIP, 0, sizeof(double) * NX_INS * NX_INS);

	phi = PHI;
	c = PHIP;
	for (j = 0; j < 3; ++j)
	{
		for (i = 0; i < NX_INS; ++i)
		{
			c[i] = phi[0] * P[i] + phi[1] * P[i + NX_INS] + phi[2] * P[i + 2 * NX_INS] + phi[3 + j] * P[i + (3 + j) * NX_INS];
		}
		phi += NX_INS;
		c += NX_INS;
	}

	for (j = 0; j < 3; ++j)
	{
		for (i = 0; i < NX_INS; ++i)
		{
			c[i] = phi[0] * P[i] + phi[1] * P[i + NX_INS] + phi[2] * P[i + 2 * NX_INS] + phi[3] * P[i + 3 * NX_INS] + phi[4] * P[i + 4 * NX_INS] + phi[5] * P[i + 5 * NX_INS]
				+ phi[6] * P[i + 6 * NX_INS] + phi[7] * P[i + 7 * NX_INS] + phi[8] * P[i + 8 * NX_INS]
				+ phi[12] * P[i + 12 * NX_INS] + phi[13] * P[i + 13 * NX_INS] + phi[14] * P[i + 14 * NX_INS];
		}
		phi += NX_INS;
		c += NX_INS;
	}

	for (j = 0; j < 3; ++j)
	{
		for (i = 0; i < NX_INS; ++i)
		{
			c[i] = phi[6] * P[i + 6 * NX_INS] + phi[7] * P[i + 7 * NX_INS] + phi[8] * P[i + 8 * NX_INS]
				+ phi[9] * P[i + 9 * NX_INS] + phi[10] * P[i + 10 * NX_INS] + phi[11] * P[i + 11 * NX_INS];
		}
		phi += NX_INS;
		c += NX_INS;
	}

	for (j = 0; j < 6; ++j)
	{
		for (i = 0; i < NX_INS; ++i)
		{
			c[i] = phi[9 + j] * P[i + (9 + j) * NX_INS];
		}
		phi += NX_INS;
		c += NX_INS;
	}

	if (IsOdo)
	{
		for (i = 0; i < NX_INS; ++i)
		{
			c[i] = P[i + (NX_INS - 1) * NX_INS];
		}
	}
}

//PHI*Q
void PHI_Q(const double* PHI, const double* Q, double* PHIQ)
{
	unsigned char j;
	const double* phi;
	double* c;

	memset(PHIQ, 0, sizeof(double) * NX_INS * NX_INS);
	phi = PHI;
	c = PHIQ;
	for (j = 0; j < 3; ++j)
	{
		c[0] = phi[0] * Q[0];
		c[1] = phi[1] * Q[1];
		c[2] = phi[2] * Q[2];
		c[3 + j] = phi[3 + j] * Q[3 + j];
		phi += NX_INS;
		c += NX_INS;
	}

	for (j = 0; j < 3; ++j)
	{
		c[0] = phi[0] * Q[0];
		c[1] = phi[1] * Q[1];
		c[2] = phi[2] * Q[2];

		c[3] = phi[3] * Q[3];
		c[4] = phi[4] * Q[4];
		c[5] = phi[5] * Q[5];

		c[6] = phi[6] * Q[6];
		c[7] = phi[7] * Q[7];
		c[8] = phi[8] * Q[8];

		c[12] = phi[12] * Q[12];
		c[13] = phi[13] * Q[12];
		c[14] = phi[4] * Q[14];

		phi += NX_INS;
		c += NX_INS;
	}

	for (j = 0; j < 3; ++j)
	{
		c[6] = phi[6] * Q[6];
		c[7] = phi[7] * Q[7];
		c[8] = phi[8] * Q[8];

		c[9] = phi[9] * Q[9];
		c[10] = phi[10] * Q[10];
		c[11] = phi[11] * Q[11];

		phi += NX_INS;
		c += NX_INS;
	}

	for (j = 0; j < 6; ++j)
	{
		c[9 + j] = phi[j + 9] * Q[j + 9];
		phi += NX_INS;
		c += NX_INS;
	}

	if (IsOdo)
	{
		c[15] = Q[21];
	}
}

//PHIP*PHI_T
void PHIP_PHIT(const double* PHIP, const double* PHI, double* PHIPPHIT)
{
	unsigned char i, j;
	const double* phi, * phip;
	double* c;
	phi = PHI;
	c = PHIPPHIT;
	phip = PHIP;

	for (i = 0; i < NX_INS; ++i)
	{
		phi = PHI;
		for (j = 0; j < 3; ++j)
		{
			c[j] = phip[0] * phi[0] + phip[1] * phi[1] + phip[2] * phi[2] + phip[3 + j] * phi[3 + j];
			phi += NX_INS;
		}
		for (j = 3; j < 6; ++j)
		{
			c[j] = phip[0] * phi[0] + phip[1] * phi[1] + phip[2] * phi[2]
				+ phip[3] * phi[3] + phip[4] * phi[4] + phip[5] * phi[5]
				+ phip[6] * phi[6] + phip[7] * phi[7] + phip[8] * phi[8]
				+ phip[12] * phi[12] + phip[13] * phi[13] + phip[14] * phi[14];
			phi += NX_INS;
		}
		for (j = 6; j < 9; ++j)
		{
			c[j] = phip[9] * phi[9] + phip[10] * phi[10] + phip[11] * phi[11]
				+ phip[6] * phi[6] + phip[7] * phi[7] + phip[8] * phi[8];
			phi += NX_INS;
		}

		for (j = 9; j < 15; ++j)
		{
			c[j] = phip[j] * phi[j];
			phi += NX_INS;
		}

		if (IsOdo)
		{
			c[15] = phip[15];
		}

		phip += NX_INS;
		c += NX_INS;
	}
}

//0.5*(PHI*Q+Q*PHIT)*dt
void PHIQ_QPHIT(double* PHIQ, double dt, double* c)
{
	unsigned char i, j;
	dt = dt * 0.5;
	for (i = 0; i < NX_INS; ++i)
	{
		for (j = 0; j < NX_INS; ++j)
		{
			c[i + j * NX_INS] = dt * (PHIQ[i + j * NX_INS] + PHIQ[j + i * NX_INS]);
		}
	}
}

/* KF functions */
int8_t KF_feedback(Par* mPar, double x[16], Nav* mNav, int IsBiasFeedBack, GnssInsSystem *mGnssInsSystem)
{
	int ret = -1;
	int i = 0;
	//Î»ÖÃ¸üÐÂ
	double d_lat = x[0] / (mNav->height + mPar->Rm);
	double d_lon = x[1] / (mPar->Rn + mNav->height) / cos(mNav->lat);
	double d_theta[3], d_ftheta[3];
	dpos2rvec(mNav->lat, d_lat, d_lon, d_theta);
	for (i = 0; i < 3; i++)
	{
		d_ftheta[i] = -d_theta[i];
	}
	double qn[4];
	rvec2quat(d_ftheta, qn);
	double q_ne[4];
	quatprod(mNav->q_ne, qn, q_ne);
	double r_n[3];
	quat2pos(q_ne, r_n);
	mNav->lat = r_n[0];
	mNav->lon = r_n[1];
	mNav->height = mNav->height + x[2];
	pos2quat(mNav->lat, mNav->lon, mNav->q_ne);


	//ËÙ¶È¸üÐÂ
	double eye[3][3] = { 1,0,0,0,1,0,0,0,1 };
	double C_cn[3][3], C1[3][3];
	GetSkewSymmetricMatrixOfVector(d_theta, *C1);
	MatrixAdd(*eye, *C1, 3, 3, *C_cn);

	double v1[3], v[3];

	v1[0] = mNav->vn - x[3];
	v1[1] = mNav->ve - x[4];
	v1[2] = mNav->vd - x[5];

	MatrixMutiply(*C_cn, v1, 3, 3, 1, v);
	mNav->vn = v[0];
	mNav->ve = v[1];
	mNav->vd = v[2];

	//×ËÌ¬¸üÐÂ
	double phi_ang[3];
	double  phi[3] = { x[6],x[7],x[8] };
	MatrixAdd(phi, d_theta, 3, 1, phi_ang);
	double qnp[4];
	rvec2quat(phi_ang, qnp);
	double qbn[4];
	quatprod(qnp, mNav->q_bn, qbn);
	memcpy(mNav->q_bn, qbn, 4 * sizeof(double));

	quat2dcm(mNav->q_bn, mNav->c_bn);
	dcm2euler(mNav->c_bn, mNav->att);
	mNav->roll = mNav->att[0];
	mNav->pitch = mNav->att[1];
	mNav->heading = mNav->att[2];



	//ÁãÆ«Ð£Õý
	if (IsBiasFeedBack)
	{
		mNav->sensorbias.bias_gyro_x += x[9];
		mNav->sensorbias.bias_gyro_y += x[10];
		mNav->sensorbias.bias_gyro_z += x[11];
		mNav->sensorbias.bias_acc_x += x[12];
		mNav->sensorbias.bias_acc_y += x[13];
		mNav->sensorbias.bias_acc_z += x[14];
	}


	//±ÈÀýÒò×ÓÎó²îÐ£Õý
	if (!(mGnssInsSystem->PerMeasUpdata&ZuptUpdate)
		&& (mGnssInsSystem->mInsData.timestamp - mGnssInsSystem->lastGNSSUseTime < 1.2)
		&& fabs(mGnssInsSystem->OdoV) > 9.0
		&&fabs(mGnssInsSystem->mNav.wnb_n[2]) < 0.03)
	{
		mNav->sensorbias.Odo_scale -= x[15];
	}

	memset(x, 0, 16 * sizeof(double));

	return 1;
}


int8_t DataChangeLC(const double dt, const ImuData* mImuData, INS *cins)
{
	int ret = -1;
	cins->timestamp = mImuData->timestamp;
	cins->timestamped = mImuData->timestamped;
	//cins->acc_x = mImuData.Gyrox / dt;
	//cins->acc_y = mImuData.Gyroy / dt;
	//cins->acc_z = mImuData.Gyroz / dt;
	//cins->gyro_x = mImuData.Accx / dt;
	//cins->gyro_y = mImuData.Accy / dt;
	//cins->gyro_z = mImuData.Accz / dt;
	cins->gyro_x = mImuData->Gyrox;
	cins->gyro_y = mImuData->Gyroy;
	cins->gyro_z = mImuData->Gyroz;
	cins->acc_x = mImuData->Accx;
	cins->acc_y = mImuData->Accy;
	cins->acc_z = mImuData->Accz;
	return 1;
}
int8_t compensate(const double dt, const Sensorbias* bias, INS* ins)
{
	ins->acc_x -= bias->bias_acc_x;
	ins->acc_y -= bias->bias_acc_y;
	ins->acc_z -= bias->bias_acc_z;
	ins->gyro_x -= bias->bias_gyro_x;
	ins->gyro_y -= bias->bias_gyro_y;
	ins->gyro_z -= bias->bias_gyro_z;
	return 1;
}



int8_t INS_MECH(const INS* ins_pre, const INS* ins_cur, const Nav* nav, Nav* nav1, Par* mPar)
{
	unsigned char i;
	double temp1[3] = { 0.0 };
	double temp2[3] = { 0.0 };
	double temp3[3] = { 0.0 };
	double q1[4] = { 0.0 };
	double q2[4] = { 0.0 };
	double q3[4] = { 0.0 };
	double mid_r[3] = { 0.0 };
	double mid_v[3] = { 0.0 };
	double Cn[3][3] = { 0.0 };
	double C1[3][3] = { 0.0 };
	double C2[3][3] = { 0.0 };
	double C3[3][3] = { 0.0 };

	double dv_fb[3] = { 0.0 };
	double dv_fn[3] = { 0.0 };
	double dv_n[3] = { 0.0 };

	//Calculat earthhparents  use navparements
	double dt = ins_cur->timestamp - ins_pre->timestamp;
	//Date prepare
	double dvel_b_prev[3] = { ins_pre->acc_x*dt, ins_pre->acc_y * dt,ins_pre->acc_z * dt };
	double dtheta_b_prev[3] = { ins_pre->gyro_x * dt, ins_pre->gyro_y * dt, ins_pre->gyro_z * dt };
	double dvel_b_cur[3] = { ins_cur->acc_x * dt, ins_cur->acc_y * dt, ins_cur->acc_z * dt };
	double dtheta_b_cur[3] = { ins_cur->gyro_x * dt,ins_cur->gyro_y * dt,ins_cur->gyro_z * dt };

	double coning[3] = { 0.0 };
	CrossProduct(dtheta_b_prev, dtheta_b_cur, temp1);
	for (int i = 0; i < 3; i++)
	{
		coning[i] = 1.0 / 12.0 * temp1[i];
	}
	double sculling[3] = { 0.0 };
	CrossProduct(dtheta_b_cur, dvel_b_cur, temp1);
	CrossProduct(dtheta_b_prev, dvel_b_cur, temp2);
	CrossProduct(dvel_b_prev, dtheta_b_cur, temp3);
	for (int i = 0; i < 3; i++)
	{
		sculling[i] = 0.5*temp1[i] + (temp2[i] + temp3[i]) / 12;
		dv_fb[i] = dvel_b_cur[i] + sculling[i];
	}
	//V*************************Velcocity Updata***************************************
	//position extrapolation
	//earth and angular rate updateing
	double pos_pre[3] = { nav->lat, nav->lon ,nav->height };
	double v_pre[3] = { nav->vn, nav->ve, nav->vd };
	double att_pre[3] = { nav->roll, nav->pitch, nav->heading };

	double M = 0.0, N = 0.0, wnie[3] = { 0.0 }, wnen[3] = { 0.0 }, ng = 0.0;

	UpdateMN(pos_pre, &M, &N);
	UpdateW(pos_pre, M, N,  v_pre,wnie, wnen);
	// navigation frame rotation vector
	//zetahalf t-1µ½t-0.5 µÄnÏµÐý×ªÊ¸Á¿
	double zeta[3] = { 0.0 }, zetahalf[3] = { 0.0 };
	//double cn[3][3] = { 0.0 }, cn1[3][3] = { 0.0 };
	for (i = 0; i < 3; i++)
	{
		zetahalf[i] = 0.5*(wnie[i] + wnen[i]) * dt;
	}
	rvec2quat(zetahalf, q1);  // q_n
	temp3[0] = 0.0;
	temp3[1] = 0.0;
	temp3[2] = -WGS84.wie*dt / 2;
	rvec2quat(temp3, q2);//q2£ºq_e
	quatprod(nav->q_ne, q1, q3);
	quatprod(q2, q3, q1);//q1£ºÖÐ¼äÊ±¿ÌµÄq_ne
	quat2pos(q1, mid_r);
	mid_r[2] = pos_pre[2] - 0.5 * v_pre[2] * dt;

	//extrapolate velocity
	mid_v[0] = nav->vn + 0.5 * nav->dv_n[0];
	mid_v[1] = nav->ve + 0.5 * nav->dv_n[1];
	mid_v[2] = nav->vd + 0.5 * nav->dv_n[2];

	//¼ÆËãn_k-1µ½n_kµÄÐý×ªÊ¸Á¿
	UpdateMN(pos_pre, &M, &N);
	UpdateW(mid_r,M,N, mid_v, wnie, wnen);

	for (i = 0; i < 3; i++)
	{
		zetahalf[i] = 0.5*(wnie[i] + wnen[i]) * dt;
	}
	double eye[3][3] = { 1,0,0,0,1,0,0,0,1 };
	GetSkewSymmetricMatrixOfVector(zetahalf, *C1);
	MatrixSub(*eye, *C1, 3, 3, *Cn);
	MatrixMutiply(*Cn, *nav->c_bn, 3, 3, 3, *C2);
	MatrixMutiply(*C2, dv_fb, 3, 3, 1, dv_fn);
	//ÖØÁ¦¡¢¸çÊÏ¼ÓËÙ¶ÈµÈÒýÆðµÄËÙ¶ÈÔöÁ¿
	UpdateGravity(mid_r, &ng);
	double gn[3] = { 0.0,0.0,ng };

	for (i = 0; i < 3; i++)
	{
		temp1[i] = 2 * wnie[i] + wnen[i];
	}
	CrossProduct(temp1, mid_v, temp2);
	for (i = 0; i < 3; i++)
	{
		dv_n[i] = dv_fn[i] + (gn[i] - temp2[i])*dt;
		nav1->dv_n[i] = dv_n[i];
		nav1->a_n[i] = dv_n[i] / dt;
	}
	nav1->vn = nav->vn + dv_n[0];
	nav1->ve = nav->ve + dv_n[1];
	nav1->vd = nav->vd + dv_n[2];
	double v_cur[3] = { nav1->vn, nav1->ve, nav1->vd };

	for (i = 0; i < 3; i++)
	{
		nav1->mSensorOut.f_n[i] = dv_n[i] / dt;
	}

	//**************************Postion Updata**************************
	for (i = 0; i < 3; i++)
	{
		mid_v[i] = 0.5 * (v_pre[i] + v_cur[i]);
	}
	UpdateW(mid_r,  M, N, mid_v, wnie, wnen);
	for (i = 0; i < 3; i++)
	{
		zeta[i] = (wnie[i] + wnen[i]) * dt;
	}
	rvec2quat(zeta, q1);  // q_n
	temp3[0] = 0.0;
	temp3[1] = 0.0;
	temp3[2] = -WGS84.wie*dt;

	rvec2quat(temp3, q2);//q2£ºq_e
	quatprod(nav->q_ne, q1, q3);
	quatprod(q2, q3, nav1->q_ne);
	norm_quat(nav1->q_ne);
	double cur_r[3] = { 0.0 };
	quat2pos(nav1->q_ne, cur_r);  //´Ë´¦mid_r±íÊ¾kÊ±¿ÌÎ»ÖÃ

	nav1->lat = cur_r[0];
	nav1->lon = cur_r[1];
	nav1->height = nav->height - mid_v[2] * dt;
	//*************************Attitude Updata*****************************
		// body frame rotation

		//Î»ÖÃÄÚ²å
	quatinv(nav->q_ne, q3);
	quatprod(q3, nav1->q_ne, q2); //q2 ÖÐ¼äÊ±¿Ìqne
	quat2rvec(q2, temp3);
	for (int i = 0; i < 3; ++i)
	{
		temp3[i] /= 2;
	}
	rvec2quat(temp3, q3);
	quatprod(nav->q_ne, q3, q2);//µÃµ½q£¨ÖÐ¼äÊ±¿ÌµÄq_ne£©

	quat2pos(q2, mid_r);
	mid_r[2] = (nav->height + nav1->height) / 2;

	//Calculate rotational and sculling motion
	double rotational[3] = { 0.0 };
	CrossProduct(dtheta_b_prev, dtheta_b_cur, temp1);
	for (int i = 0; i < 3; i++)
	{
		rotational[i] = temp1[i] / 12;
	}
	for (i = 0; i < 3; ++i)
	{
		temp2[i] = dtheta_b_cur[i] + rotational[i];
	}
	rvec2quat(temp2, q2); //q2£ºq_b

	//ÔÙ´Î¼ÆËãnÏµk-1Ê±¿Ìµ½kÊ±¿ÌµÄÐý×ªÊ¸Á¿
	UpdateW(mid_r, M, N,  mid_v, wnie, wnen);
	for (i = 0; i < 3; i++)
	{
		zeta[i] = (wnie[i] + wnen[i]) * dt;
	}
	for (i = 0; i < 3; i++)
	{
		temp3[i] = -zeta[i];
	}
	rvec2quat(temp3, q1);  //q1 qn

	quatprod(nav->q_bn, q2, q3);
	quatprod(q1, q3, nav1->q_bn);
	norm_quat(nav1->q_bn);
	quat2dcm(nav1->q_bn, nav1->c_bn);
	dcm2euler(nav1->c_bn, nav1->att);

	nav1->roll = nav1->att[0];
	nav1->pitch = nav1->att[1];
	nav1->heading = nav1->att[2];



	memcpy(mPar->w_ie, wnie, 3 * sizeof(double));
	memcpy(mPar->w_en, wnen, 3 * sizeof(double));
	for (i = 0; i < 3; i++)
	{
		mPar->f_n[i] = dv_fn[i] / dt;
		mPar->f_b[i] = dvel_b_cur[i] / dt;
		mPar->w_b[i] = dtheta_b_cur[i] / dt;
		temp2[i] = (wnie[i] + wnen[i]);
	}
	MatrixMutiply(*(nav1->c_bn), mPar->w_b, 3, 3, 1, temp1);
	for (i = 0; i < 3; i++)
	{
		nav1->wnb_n[i] = temp1[i] - temp2[i];
	}

	mPar->Rm = M;
	mPar->Rn = N;

	UpdateGravity(mid_r, &ng);
	mPar->g[0] = 0.0;
	mPar->g[1] = 0.0;
	mPar->g[2] = ng;
	return 1;
}

/***************************************************************************************
Function: KF_predict_16PHI
Description: ;Calculate PHI(transition matrix)
Input :dt  time increment
	   nav  Navigation information
Output:PHI  transition matrix
Return :
Others:
********************************************************************************************/
int8_t KF_predict_16PHI(const double dt, const Nav *mNav, const Par *mPar, const ImuSensor *mImuSensor, double PHI[16][16])
{
	double M = mPar->Rm;
	double N = mPar->Rn;
	double wie[3], wen[3];
	memcpy(wie, mPar->w_ie, 3 * sizeof(double));
	memcpy(wen, mPar->w_en, 3 * sizeof(double));
	double pos[3] = { mNav->lat, mNav->lon, mNav->height };
	double g[3];
	memcpy(g, mPar->g, 3 * sizeof(double));

	double  f_b[3], f_n[3];
	memcpy(f_b, mPar->f_b, 3 * sizeof(double));
	memcpy(f_n, mPar->f_n, 3 * sizeof(double));
	double cbn[3][3];
	memcpy(cbn, mNav->c_bn, 9 * sizeof(double));
	//PHI = 16x16 transition matrix
	double eye[3][3] = { 1,0,0,0,1,0,0,0,1 };
	memset(PHI, 0, 256 * sizeof(double));
	//Postion error dynamics
	//pos to pos
	double fwnen[3], Mwnen[3][3];
	for (int i = 0; i < 3; i++)
	{
		fwnen[i] = -wen[i];
	}
	GetSkewSymmetricMatrixOfVector(fwnen, *Mwnen);
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			PHI[i][j] = eye[i][j] + Mwnen[i][j] * dt;
		}
	}
	//pos to vel
	for (int i = 0; i < 3; i++)
	{
		for (int j = 3; j < 6; j++)
		{
			PHI[i][j] = eye[i][j - 3] * dt;
		}
	}
	//Velocity error dynamics
	//vel to pose
	double R = sqrt(M * N);
	PHI[3][0] = -g[2] / (R + pos[2])*dt;
	PHI[4][1] = PHI[3][0];
	PHI[5][2] = -2 * PHI[3][0];
	//vel to vel
	double w1[3];
	for (int i = 0; i < 3; i++)
	{
		w1[i] = -2 * wie[i] - wen[i];
	}
	double mw1[3][3];
	GetSkewSymmetricMatrixOfVector(w1, *mw1);
	for (int i = 3; i < 6; i++)
	{
		for (int j = 3; j < 6; j++)
		{
			PHI[i][j] = eye[i - 3][j - 3] + mw1[i - 3][j - 3] * dt;
		}
	}
	//vel to att
	double mfn[3][3];
	GetSkewSymmetricMatrixOfVector(f_n, *mfn);
	for (int i = 3; i < 6; i++)
	{
		for (int j = 6; j < 9; j++)
		{
			PHI[i][j] = mfn[i - 3][j - 6] * dt;
		}
	}
	//vel to acc bias
	for (int i = 3; i < 6; i++)
	{
		for (int j = 12; j < 15; j++)
		{
			PHI[i][j] = cbn[i - 3][j - 12] * dt;
		}
	}
	//attitude error dynamic
	double fwnin[3];
	for (int i = 0; i < 3; i++)
	{
		fwnin[i] = -wie[i] - wen[i];
	}
	double mfwnin[3][3];
	GetSkewSymmetricMatrixOfVector(fwnin, *mfwnin);
	for (int i = 6; i < 9; i++)
	{
		for (int j = 6; j < 9; j++)
		{
			PHI[i][j] = eye[i - 6][j - 6] + mfwnin[i - 6][j - 6] * dt;
		}
	}
	//att to gyro bias
	for (int i = 6; i < 9; i++)
	{
		for (int j = 9; j < 12; j++)
		{
			PHI[i][j] = -cbn[i - 6][j - 9] * dt;
		}
	}

	//gyro bias
	PHI[9][9] = mImuSensor->bg_model[1];
	PHI[10][10] = mImuSensor->bg_model[1];
	PHI[11][11] = mImuSensor->bg_model[1];
	//acc bias
	PHI[12][12] = mImuSensor->ba_model[1];
	PHI[13][13] = mImuSensor->ba_model[1];
	PHI[14][14] = mImuSensor->ba_model[1];
	//odo
	PHI[15][15] = 1.0;
	return 1;
}
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
int8_t  KF_predict(const double PHI[16][16], const double Q[16], const double dt, double x[16], double P[16][16], UpdataStruct* mUpdataStruct)
{
	//MatrixMutiply(*PHI, x, 16, 16, 1, x1);
	double Qd[16][16];
	double M1[16][16],M2[16][16];

	PHI_Q(*PHI, Q, *M1);//PHI*Q
	PHIQ_QPHIT(*M1, dt, *Qd);


	//ÀÛ¼ÆQ£¬ÓÃÓÚ×´Ì¬×ªÒÆ
	PHI_P(*PHI, *(mUpdataStruct->Q), *M1); //PHI*Q
	PHIP_PHIT(*M1, *PHI, *M2);  //PHI*Q*PHI_T
	MatrixAdd(*M2, *Qd, 16, 16, *(mUpdataStruct->Q));

	//ÀÛ¼ÆPHI£¬ÓÃÓÚ×´Ì¬×ªÒÆ
	PHI_P(*PHI, *(mUpdataStruct->PHI), *M1);
	memcpy(mUpdataStruct->PHI, M1, 16 * 16 * sizeof(double));

	//¸üÐÂP
	PHI_P(*PHI, *P, *M1);
	PHIP_PHIT(*M1, *PHI, *M2);
	MatrixAdd(*M2, *Qd, 16, 16, *P);
	return 1;
}

/*--------------------------------------------------------*/

static int8_t MotionAlign(const GnssData* premGnssData, const GnssData* mGnssData, uint8_t GnssType, Nav *mNav)
{
	int8_t retal = 0;
	double Vel_NED[3] = { 0 };
	double vehchile_speed = 0.0;
	if (GnssType == 0)
	{
		double prepos[3] = { premGnssData->latitude, premGnssData->longitude ,premGnssData->altitude };
		double pos[3] = { mGnssData->latitude, mGnssData->longitude ,mGnssData->altitude };

		double M, N;
		UpdateMN(pos, &M, &N);
		double dt = mGnssData->timestamp - premGnssData->timestamp;
		if (premGnssData->timestamp < 0) return 0;
		if (dt < 2.0 && dt > 0.9)
		{
			Vel_NED[0] = (M + pos[2])*(pos[0] - prepos[0]) / dt;
			Vel_NED[1] = (N + pos[2])*(pos[1] - prepos[1])*cos(pos[0]) / dt;
			Vel_NED[2] = -(pos[2] - prepos[2]) / dt;
			vehchile_speed = sqrt(Vel_NED[0] * Vel_NED[0] + Vel_NED[1] * Vel_NED[1]);
		}
	}
	if (GnssType == 1)
	{
		Vel_NED[0] = mGnssData->north_velocity;
		Vel_NED[1] = mGnssData->east_velocity;
		Vel_NED[2] = -mGnssData->up_velocity;
		vehchile_speed = sqrt(Vel_NED[0] * Vel_NED[0] + Vel_NED[1] * Vel_NED[1]);
	}
	if (vehchile_speed > 0.5) //
	{
		mNav->lon = mGnssData->longitude;
		mNav->lat = mGnssData->latitude;
		mNav->height = mGnssData->altitude;
		mNav->vn = Vel_NED[0];
		mNav->ve = Vel_NED[1];
		mNav->vd = Vel_NED[2];
		mNav->roll = 0.0;
		mNav->pitch = 0.0;
		double heading = atan2(mNav->ve, mNav->vn);
		if (heading < 0)
		{
			heading += 2 * PI;
		}
		if (heading > 2 * PI)
		{
			heading -= 2 * PI;
		}
		mNav->heading = heading;


		retal = 1;
	}
	return retal;
}
static int8_t DualantennaAlign(const GnssData* premGnssData, const GnssData* mGnssData, uint8_t GnssType, Nav *mNav)
{
	int8_t retal = 0;
	double Vel_NED[3];
	double vehchile_speed = 0.0;
	
	if (GnssType == 0)
	{
		double prepos[3] = { premGnssData->latitude, premGnssData->longitude ,premGnssData->altitude };
		double pos[3] = { mGnssData->latitude, mGnssData->longitude ,mGnssData->altitude };

		double M, N;
		UpdateMN(pos, &M, &N);
		double dt = (mGnssData->timestamp - premGnssData->timestamp);
		if (dt < 1.1 && dt > 0.9)
		{
			Vel_NED[0] = (M + pos[2])*(pos[0] - prepos[0]) / dt;
			Vel_NED[1] = (N + pos[2])*(pos[1] - prepos[1])*cos(pos[0]) / dt;
			Vel_NED[2] = -(pos[2] - prepos[2]) / dt;
			vehchile_speed = sqrt(Vel_NED[0] * Vel_NED[0] + Vel_NED[1] * Vel_NED[1]);
		}
	}
	if (GnssType == 1)
	{
		Vel_NED[0] = mGnssData->north_velocity;
		Vel_NED[0] = mGnssData->east_velocity;
		Vel_NED[0] = -mGnssData->up_velocity;
		vehchile_speed = sqrt(Vel_NED[0] * Vel_NED[0] + Vel_NED[1] * Vel_NED[1]);
	}
	{
		mNav->lon = mGnssData->longitude;
		mNav->lat = mGnssData->latitude;
		mNav->height = mGnssData->altitude;
		mNav->vn = 0;// Vel_NED[0];
		mNav->ve = 0; //Vel_NED[1];
		mNav->vd = 0;// Vel_NED[2];
		mNav->roll = 0.0;
		mNav->pitch = 0.0;
		double heading = 3.14;//mGnssData->heading;
		if (heading < 0)
		{
			heading += 2 * PI;
		}
		if (heading >= 2 * PI)
		{
			heading -= 2 * PI;
		}
		mNav->heading = heading;


		retal = 1;
	}
	return retal;
}


static int8_t SetMeasZ(int16_t mMeasUpdataType, GnssData mGnssData)
{
	int ret = -1;
	if (mMeasUpdataType & GNSSPoitionUpdate)
	{
		memset(mGnssInsSystem.mKalmanStruct.Z, 0, sizeof(mGnssInsSystem.mKalmanStruct.Z));

		ECEF r_gps_E = { 0 };
		Geo r_gps = { 0 };

		r_gps.lat = mGnssData.latitude;
		r_gps.lon = mGnssData.longitude;
		r_gps.height = mGnssData.altitude;
		llh2ecef(&r_gps, &r_gps_E);
		double r_gps_e[3] = { r_gps_E.x, r_gps_E.y, r_gps_E.z };


		Geo r_ins = { 0 };
		ECEF r_ins_E = { 0 };
		r_ins.lat = mUpdataStruct.mPVA.latitude;
		r_ins.lon = mUpdataStruct.mPVA.longitude;
		r_ins.height = mUpdataStruct.mPVA.altitude;
		llh2ecef(&r_ins, &r_ins_E);
		double r_ins_e[3] = { r_ins_E.x, r_ins_E.y, r_ins_E.z };

		double Cne[3][3] = { 0.0 };
		double Cen[3][3] = { 0.0 };
		pos2dcm(r_gps.lat, r_gps.lon, Cne);
		MatrixTranspose(*Cne, 3, 3, *Cen);


		double leverarm_b[3] = { 0.0 };
		MatrixMutiply(*mGnssInsSystem.C_InstallationAngle, mLCSetting.Leverarm, 3, 3, 1, leverarm_b);

		double  larm_rn[3] = { 0.0 };
		MatrixMutiply(*mUpdataStruct.C_bn, leverarm_b, 3, 3, 1, larm_rn);

		double diff_r_e[3] = { 0.0 }, diff_r_n[3] = { 0.0 };
		MatrixSub(r_ins_e, r_gps_e, 3, 1, diff_r_e);
		MatrixMutiply(*Cen, diff_r_e, 3, 3, 1, diff_r_n);
		MatrixAdd(diff_r_n, larm_rn, 3, 1, mGnssInsSystem.mKalmanStruct.Z);

		double Mlarm_rn[3][3] = { 0.0 };
		GetSkewSymmetricMatrixOfVector(larm_rn, *Mlarm_rn);

		memset(mGnssInsSystem.mKalmanStruct.H, 0, sizeof(mGnssInsSystem.mKalmanStruct.H));
		mGnssInsSystem.mKalmanStruct.H[0][0] = 1.0;
		mGnssInsSystem.mKalmanStruct.H[1][1] = 1.0;
		mGnssInsSystem.mKalmanStruct.H[2][2] = 1.0;

		memset(mGnssInsSystem.mKalmanStruct.R, 0, sizeof(mGnssInsSystem.mKalmanStruct.R));
		mGnssInsSystem.mKalmanStruct.R[0][0] = mGnssData.latitude_std * mGnssData.latitude_std;
		mGnssInsSystem.mKalmanStruct.R[1][1] = mGnssData.longitude_std * mGnssData.longitude_std;
		mGnssInsSystem.mKalmanStruct.R[2][2] = mGnssData.altitude_std * mGnssData.altitude_std;
	}

	if (mMeasUpdataType&ZuptUpdate)
	{
		memset(mGnssInsSystem.mKalmanStruct.Z, 0, sizeof(mGnssInsSystem.mKalmanStruct.Z));
		mGnssInsSystem.mKalmanStruct.Z[0] = mUpdataStruct.mPVA.northVelocity;
		mGnssInsSystem.mKalmanStruct.Z[1] = mUpdataStruct.mPVA.eastVelocity;
		mGnssInsSystem.mKalmanStruct.Z[2] = mUpdataStruct.mPVA.downVelocity;

		memset(mGnssInsSystem.mKalmanStruct.H, 0, sizeof(mGnssInsSystem.mKalmanStruct.H));
		mGnssInsSystem.mKalmanStruct.H[0][3] = 1;
		mGnssInsSystem.mKalmanStruct.H[1][4] = 1;
		mGnssInsSystem.mKalmanStruct.H[2][5] = 1;

		memset(mGnssInsSystem.mKalmanStruct.R, 0, sizeof(mGnssInsSystem.mKalmanStruct.R));
		mGnssInsSystem.mKalmanStruct.R[0][0] = 0.1 * 0.1;
		mGnssInsSystem.mKalmanStruct.R[1][1] = 0.1 * 0.1;
		mGnssInsSystem.mKalmanStruct.R[2][2] = 0.1 * 0.1;
	}

	if (mMeasUpdataType&NHCUpdate)
	{
		double Mwnb_n[3][3] = { 0.0 }, Mbv[3] = { 0.0 };
		GetSkewSymmetricMatrixOfVector(mUpdataStruct.mPVA.wnb_n, *Mwnb_n);
		

		MatrixMutiply(*Mwnb_n, mLCSetting.OdoLeverarm, 3, 3, 1, Mbv);   
		double Vel_NED[3] = { mUpdataStruct.mPVA.northVelocity, mUpdataStruct.mPVA.eastVelocity,  mUpdataStruct.mPVA.downVelocity };
		double INS_Nv_Odo[3] = { 0.0 };
		MatrixAdd(Vel_NED, Mbv, 3, 1, INS_Nv_Odo); 

		double win_n[3] = { 0.0 }, Mwin[3][3] = { 0.0 };
		double pos_pre[3] = { mUpdataStruct.mPVA.latitude, mUpdataStruct.mPVA.longitude,mUpdataStruct.mPVA.altitude };
		double M = 0.0, N = 0.0, wnie[3] = { 0.0 }, wnen[3] = { 0.0 };
		UpdateMN(pos_pre, &M, &N);
		UpdateW(pos_pre, M, N, Vel_NED, wnie, wnen);


		MatrixAdd(wnie, wnen, 3, 1, win_n);
		GetSkewSymmetricMatrixOfVector(win_n, *Mwin);

		double C_bv[3][3] = { 0.0 }, C_vb[3][3] = { 0.0 }, C_vn[3][3] = { 0.0 }, C_nv[3][3] = { 0.0 };
		euler2dcm(mLCSetting.INSRotationRvb, C_bv);
		MatrixTranspose(*C_bv, 3, 3, *C_vb);
		MatrixMutiply(*mUpdataStruct.C_bn, *C_vb, 3, 3, 3, *C_vn);
		MatrixTranspose(*C_vn, 3, 3, *C_nv);

		double INS_Vv_Odo[3] = { 0.0 };
		MatrixMutiply(*C_nv, INS_Nv_Odo, 3, 3, 1, INS_Vv_Odo);  //V??³??????

		memset(mGnssInsSystem.mKalmanStruct.Z, 0, sizeof(mGnssInsSystem.mKalmanStruct.Z));
		mGnssInsSystem.mKalmanStruct.Z[0] = INS_Vv_Odo[1];
		mGnssInsSystem.mKalmanStruct.Z[1] = INS_Vv_Odo[2];


		//
		double ML[3][3] = { 0.0 }, CbnML[3][3] = { 0.0 }, MwinCbnML[3][3] = { 0.0 }, MLWb[3] = { 0.0 }, C_bnMLWb[3] = { 0.0 }, MC_bnMLWb[3][3] = { 0.0 };
		GetSkewSymmetricMatrixOfVector(mLCSetting.OdoLeverarm, *ML);
		MatrixMutiply(*mUpdataStruct.C_bn, *ML, 3, 3, 3, *CbnML);

		MatrixMutiply(*Mwin, *CbnML, 3, 3, 3, *MwinCbnML);

		CrossProduct(mLCSetting.OdoLeverarm, mUpdataStruct.w_b, MLWb);
		MatrixMutiply(*mUpdataStruct.C_bn, MLWb, 3, 3, 1, C_bnMLWb);
		GetSkewSymmetricMatrixOfVector(C_bnMLWb, *MC_bnMLWb);

		double C1[3][3] = { 0.0 };
		MatrixAdd(*MwinCbnML, *MC_bnMLWb, 3, 3, *C1);

		double A_psi[3][3] = { 0.0 }, A_g[3][3] = { 0.0 };
		MatrixMutiply(*C_nv, *C1, 3, 3, 3, *A_psi);
		MatrixMutiply(*C_nv, *CbnML, 3, 3, 3, *A_g);

		memset(mGnssInsSystem.mKalmanStruct.H, 0, sizeof(mGnssInsSystem.mKalmanStruct.H));
		for (int16_t i = 1; i < 3; i++)
		{
			for (int16_t j = 3; j < 6; j++)
			{
				mGnssInsSystem.mKalmanStruct.H[i - 1][j] = C_nv[i][j - 3];
			}
		}
		for (int16_t i = 1; i < 3; i++)
		{
			for (int16_t j = 6; j < 9; j++)
			{
				mGnssInsSystem.mKalmanStruct.H[i - 1][j] = A_psi[i][j - 6];
			}
		}
		for (int16_t i = 1; i < 3; i++)
		{
			for (int16_t j = 9; j < 12; j++)
			{
				mGnssInsSystem.mKalmanStruct.H[i - 1][j] = -A_g[i][j - 9];
			}
		}

		memset(mGnssInsSystem.mKalmanStruct.R, 0, sizeof(mGnssInsSystem.mKalmanStruct.R));
		mGnssInsSystem.mKalmanStruct.R[0][0] = 0.15 * 0.15;
		mGnssInsSystem.mKalmanStruct.R[1][1] = 0.15 * 0.15;

	}

	return ret;

}

static int8_t KalmanUpate()
{
	int ret = -1;
	int i;
	double Ht[16][3] = { 0.0 };
	double PHt[16][3] = { 0.0 };
	double HPHt[3][3] = { 0.0 };
	double RHPHt[3][3] = { 0.0 };

	MatrixTranspose(*mGnssInsSystem.mKalmanStruct.H, 3, 16, *Ht);
	MatrixMutiply(*mUpdataStruct.P, *Ht, 16, 16, 3, *PHt);
	MatrixMutiply(*mGnssInsSystem.mKalmanStruct.H, *PHt, 3, 16, 3, *HPHt);
	MatrixAdd(*mGnssInsSystem.mKalmanStruct.R, *HPHt, 3, 3, *RHPHt);


	double U1[3][3], U[3][3], Uit[3][3];
	int16_t ndf = MatrixCholosky(*RHPHt, 3, *U1);


	if (1 == ndf)
	{
		double Ui[3][3], UU[3][3], K[16][3];
		//MatrixTranspose(*U1, 3, 3, *U);
		//MatrixInverse(*U, 3, *Ui);
		//MatrixTranspose(*Ui, 3, 3, *Uit);

		//MatrixMutiply(*Ui, *Uit, 3, 3, 3, *UU);
		memcpy(*UU, *RHPHt, 3 * 3 * sizeof(double));
		MatrixInverse(3,*UU);
		MatrixMutiply(*PHt, *UU, 16, 3, 3, *K);


		double HX[3], inno[3], dx[16];
		MatrixMutiply(*mGnssInsSystem.mKalmanStruct.H, mUpdataStruct.X, 3, 16, 1, HX);
		MatrixSub(mGnssInsSystem.mKalmanStruct.Z, HX, 3, 1, inno);
		MatrixMutiply(*K, inno, 16, 3, 1, dx);
		MatrixAdd(mUpdataStruct.X, dx, 16, 1, mUpdataStruct.X);

		//????P
		double KH[16][16], IKH[16][16], IKHt[16][16], eye[16][16], IKHP[16][16];
		memset(eye, 0, 16 * 16 * sizeof(double));
		for (int16_t i = 0; i < 16; i++)
		{
			eye[i][i] = 1;
		}
		MatrixMutiply(*K, *mGnssInsSystem.mKalmanStruct.H, 16, 3, 16, *KH);
		MatrixSub(*eye, *KH, 16, 16, *IKH);
		MatrixMutiply(*IKH, *mUpdataStruct.P, 16, 16, 16, *IKHP);
		MatrixTranspose(*IKH, 16, 16, *IKHt);

		double IKHPIKHt[16][16];
		MatrixMutiply(*IKHP, *IKHt, 16, 16, 16, *IKHPIKHt);

		double Kt[3][16], KR[16][3], KRkt[16][16];
		MatrixTranspose(*K, 16, 3, *Kt);
		MatrixMutiply(*K, *mGnssInsSystem.mKalmanStruct.R, 16, 3, 3, *KR);
		MatrixMutiply(*KR, *Kt, 16, 3, 16, *KRkt);

		MatrixAdd(*IKHPIKHt, *KRkt, 16, 16, *mUpdataStruct.P);
	}
	else
	{


	}
	return 1;
}

static int8_t NhcUpdate()
{
	int ret = -1;
	int i;
	double H[2][16] = { 0.0 };
	double Ht[16][2] = { 0.0 };
	double PHt[16][2] = { 0.0 };
	double HPHt[2][2] = { 0.0 };
	double Z[2][1] = { 0.0 };
	double R[2][2] = { 0.0 };
	double RHPHt[2][2] = { 0.0 };

	memcpy(*H, *mGnssInsSystem.mKalmanStruct.H, 16 * 2 * sizeof(double));
	R[0][0] = mGnssInsSystem.mKalmanStruct.R[0][0];
	R[1][1] = mGnssInsSystem.mKalmanStruct.R[1][1];


	MatrixTranspose(*H, 2, 16, *Ht);
	MatrixMutiply(*mUpdataStruct.P, *Ht, 16, 16, 2, *PHt);
	MatrixMutiply(*H, *PHt, 2, 16, 2, *HPHt);
	MatrixAdd(*R, *HPHt, 2, 2, *RHPHt);

	double U1[2][2], U[2][2], Uit[2][2];
	int16_t ndf = MatrixCholosky(*RHPHt, 2, *U1);
	if (1 == ndf)
	{
		double Ui[2][2], UU[2][2], K[16][2];
		//MatrixTranspose(*U1, 3, 3, *U);
		//MatrixInverse(*U, 3, *Ui);
		//MatrixTranspose(*Ui, 3, 3, *Uit);

		//MatrixMutiply(*Ui, *Uit, 3, 3, 3, *UU);
		memcpy(*UU, *RHPHt, 2 * 2 * sizeof(double));
		MatrixInverse(2,*UU);
		MatrixMutiply(*PHt, *UU, 16, 2, 2, *K);

		double HX[2], inno[2], dx[16];
		MatrixMutiply(*H, mUpdataStruct.X, 2, 16, 1, HX);
		MatrixSub(mGnssInsSystem.mKalmanStruct.Z, HX, 2, 1, inno);
		MatrixMutiply(*K, inno, 16, 2, 1, dx);
		MatrixAdd(mUpdataStruct.X, dx, 16, 1, mUpdataStruct.X);

		double KH[16][16], IKH[16][16], IKHt[16][16], eye[16][16], IKHP[16][16];
		memset(eye, 0, 16 * 16 * sizeof(double));
		for (int16_t i = 0; i < 16; i++)
		{
			eye[i][i] = 1;
		}
		MatrixMutiply(*K, *H, 16, 2, 16, *KH);
		MatrixSub(*eye, *KH, 16, 16, *IKH);
		MatrixMutiply(*IKH, *mUpdataStruct.P, 16, 16, 16, *IKHP);
		MatrixTranspose(*IKH, 16, 16, *IKHt);

		double IKHPIKHt[16][16];
		MatrixMutiply(*IKHP, *IKHt, 16, 16, 16, *IKHPIKHt);

		double Kt[2][16], KR[16][2], KRkt[16][16];
		MatrixTranspose(*K, 16, 2, *Kt);
		MatrixMutiply(*K, *R, 16, 2, 2, *KR);
		MatrixMutiply(*KR, *Kt, 16, 2, 16, *KRkt);

		MatrixAdd(*IKHPIKHt, *KRkt, 16, 16, *mUpdataStruct.P);
		ret = 1;
	}
	else
	{

	}
	return ret;
}

int8_t ObsUpdata(double systemtime, GnssData mGnssData)
{
	int ret = -1;
	if (mUpdataStruct.LCTime != mGnssData.timestamp)
	{
		return -1;
	}

	//mGnssInsSystem.CurIsUseZupt = GetZuptVal();
	if (!mGnssData.flag)mGnssInsSystem.GNSSflag = 0;
	//if (mUpdataStruct.IsUseZUPT)  // mGnssInsSystem.CurIsUseZupt
	//{
	//	SetMeasZ(4, mGnssInsSystem.mGnssData);
	//	KalmanUpate();
	//	mGnssInsSystem.PerMeasUpdata = mGnssInsSystem.PerMeasUpdata | mGnssInsSystem.mMeasUpdataType;
	//	mGnssInsSystem.lastZuptTime = mGnssData.timestamp;
	//	mGnssInsSystem.lastNHCUseTime = mGnssData.timestamp;
	//	if (mGnssInsSystem.mGnssData.Mode != 4)
	//	{
	//		mGnssInsSystem.mGnssData.Mode = 0;
	//		mGnssInsSystem.GNSSflag = 0;
	//	}
	//	IsFinished = 1;
	//}


	if (mGnssInsSystem.GNSSflag == 1 && fabs(mGnssInsSystem.mGnssData.timestampd - mUpdataStruct.LCTime) < 4.0)
	{
		SetMeasZ(1, mGnssInsSystem.mGnssData);
		//Quality Control
		//if (mGnssInsSystem.mGnssData.Mode )  //!= 4
		//{
		//	if (mUpdataStruct.P[0][0] < 9 * mGnssInsSystem.mKalmanStruct.R[0][0] &&
		//		mUpdataStruct.P[1][1] < 9 * mGnssInsSystem.mKalmanStruct.R[1][1] &&
		//		mUpdataStruct.P[2][2] < 9 * mGnssInsSystem.mKalmanStruct.R[2][2])
		//	{
		//		if (fabs(mGnssInsSystem.mKalmanStruct.Z[0]) > 6 * sqrt(mUpdataStruct.P[0][0] + mGnssInsSystem.mKalmanStruct.R[0][0])
		//			|| fabs(mGnssInsSystem.mKalmanStruct.Z[1]) > 6 * sqrt(mUpdataStruct.P[1][1] + mGnssInsSystem.mKalmanStruct.R[1][1])
		//			|| fabs(mGnssInsSystem.mKalmanStruct.Z[2]) > 6 * sqrt(mUpdataStruct.P[2][2] + mGnssInsSystem.mKalmanStruct.R[2][2]))
		//		{
		//			mGnssInsSystem.GNSSflag = 0;
		//		}
		//		for (int16_t i = 0; i < 3; i++)
		//		{
		//			if (fabs(mGnssInsSystem.mKalmanStruct.Z[i]) > 2.5 * sqrt(mUpdataStruct.P[i][i] + mGnssInsSystem.mKalmanStruct.R[i][i]))
		//			{
		//				mGnssInsSystem.mKalmanStruct.R[i][i] *= fabs(mGnssInsSystem.mKalmanStruct.Z[i]) / sqrt(mUpdataStruct.P[i][i] + mGnssInsSystem.mKalmanStruct.R[i][i]);
		//			}
		//		}
		//	}
		//}
		if (mGnssInsSystem.GNSSflag != 0)
		{
			KalmanUpate();
			mGnssInsSystem.PerMeasUpdata = mGnssInsSystem.PerMeasUpdata | mGnssInsSystem.mMeasUpdataType;
			mGnssInsSystem.lastGNSSUseTime = mGnssInsSystem.lastLCTime;
			mGnssInsSystem.GNSSflag = 0;
			if (mGnssInsSystem.firstGNSSUseTime < 0.5)
			{
				mGnssInsSystem.firstGNSSUseTime = mGnssInsSystem.lastGNSSUseTime;
			}
		}
	}
	//if (mLCSetting.IsUseNHC)  //mLCSetting.IsUseNHC
	//{
	//	if (fabs(mGnssInsSystem.mGnssData.timestampd - mGnssInsSystem.lastLCTime) < 0.9
	//		&& !mGnssInsSystem.CurIsUseZupt)
	//	{
	//		int16_t mesatype = mLCSetting.IsUseOdo ? 16 : 8;
	//		SetMeasZ(mesatype, mGnssInsSystem.mGnssData);
	//		if (mLCSetting.IsUseOdo)
	//		{
	//			KalmanUpate();
	//		}
	//		else
	//		{
	//			NhcUpdate();
	//		}
	//		mGnssInsSystem.PerMeasUpdata = mGnssInsSystem.PerMeasUpdata | mesatype;
	//		mGnssInsSystem.lastNHCLCTime = mGnssInsSystem.mGnssData.timestampd;
	//		mGnssInsSystem.lastNHCUseTime = mGnssInsSystem.mGnssData.timestampd;
	//	}
	//}
	IsFinished = 1;
	return 1;

}

int8_t  ADDGNSSDATA(const GnssData mGnssData)
{
	int8_t ret = -1;

	memcpy(&(mGnssInsSystem.mGnssData), &(mGnssData), sizeof(mGnssData));
	mGnssInsSystem.GNSSflag = 1;

	switch (mGnssInsSystem.mlc_STATUS)
	{
	case INSDATA_NREADY:
	{
	}break;
	case INSDATA_READY:
	{
		switch (mLCSetting.InitialAttitudeMode)
		{
		case 0:
		{
			int8_t AlignMentflag = MotionAlign(&mGnssInsSystem.premGnssData, &mGnssInsSystem.mGnssData, mLCSetting.GnssDataType, &mGnssInsSystem.mNav);
			if (AlignMentflag)
			{
				mGnssInsSystem.mlc_STATUS = ALIGNMENT_COMPLETE;
				mGnssInsSystem.AlignCTime = mGnssData.timestamp;

				//double DeltaT = mGnssInsSystem.mImuData.timestamp - mGnssInsSystem.mpreImuData.timestamp;
				//DataChangeLC(DeltaT, &mGnssInsSystem.mImuData, &mGnssInsSystem.mpreInsData);

			    initsystemfromGNSS(&mGnssInsSystem);
                initsystemSoftreset();

				mGnssInsSystem.mlc_STATUS = INS_FUSING;
			}
		}break;
		case 1:
		{
			if (fabs(mGnssInsSystem.mGnssData.timestamp - mGnssInsSystem.premGnssData.timestamp) < 1.5)
			{
				int8_t AlignMentflag = DualantennaAlign(&mGnssInsSystem.premGnssData, &mGnssInsSystem.mGnssData, mLCSetting.GnssDataType, &mGnssInsSystem.mNav);

				if (AlignMentflag)
				{
					mGnssInsSystem.mlc_STATUS = ALIGNMENT_COMPLETE;
					mGnssInsSystem.AlignCTime = mGnssData.timestamp;

					//double DeltaT = mGnssInsSystem.mImuData.timestamp - mGnssInsSystem.mpreImuData.timestamp;
					//DataChangeLC(DeltaT, mGnssInsSystem.mImuData, &mGnssInsSystem.mpreInsData);
					initsystemfromGNSS(&mGnssInsSystem);
					initsystemSoftreset();
					mGnssInsSystem.mlc_STATUS = INS_FUSING;
				}
			}


		}break;
		default:
			return ret;
		}
	}break;
	case ALIGNMENT_COMPLETE:
	{
		mGnssInsSystem.mlc_STATUS = INS_FUSING;
	}break;
	case INS_FUSING:
	{
		ObsUpdata(mGnssData.timestampd, mGnssData);
	}break;
	default:
		return ret;
	}
	memcpy(&(mGnssInsSystem.premGnssData), &(mGnssData), sizeof(mGnssData));

	ret = 1;
	return  ret;
}

/*--------------------------------------------------------*/
static int8_t IsNeedLC(double lastinstime, double nowinstime, double NEXTLCTIME)
{
	int8_t ret = 0;
	if (lastinstime < NEXTLCTIME && (nowinstime+0.00001) >= NEXTLCTIME)
	{
		ret = NEXTLCTIME - lastinstime > nowinstime - NEXTLCTIME ? 1 : 2;
		if (!(ret == 1 && fabs(nowinstime - NEXTLCTIME) < 0.0001 || ret == 2 && fabs(NEXTLCTIME - lastinstime) < 0.0001))
		{
			ret = 3;
		}
	}
	return ret;
};

static int8_t GetLCData(Nav* mNav, double LCTime)
{
	int8_t ret = 0;
	mUpdataStruct.LCTime = LCTime;
	//mUpdataStruct.IsUseZUPT = GetZuptVal();
	mGnssInsSystem.lastLCTime = LCTime;
	//mGnssInsSystem.nextLCTime = LCTime + 1;
	mUpdataStruct.mPVA.latitude = mNav->lat;
	mUpdataStruct.mPVA.longitude = mNav->lon;
	mUpdataStruct.mPVA.altitude = mNav->height;
	mUpdataStruct.mPVA.northVelocity = mNav->vn;
	mUpdataStruct.mPVA.eastVelocity = mNav->ve;
	mUpdataStruct.mPVA.downVelocity = mNav->vd;
	mUpdataStruct.mPVA.roll = mNav->roll;
	mUpdataStruct.mPVA.pitch = mNav->pitch;
	mUpdataStruct.mPVA.yaw = mNav->heading;

	memcpy(mUpdataStruct.mPVA.wnb_n, mNav->wnb_n, 3 * sizeof(double));
	memcpy(mUpdataStruct.mPVA.a_n, mNav->a_n, 3 * sizeof(double));


	memcpy(mUpdataStruct.P, mGnssInsSystem.mKalmanStruct.P, 16 * 16 * sizeof(double));
	memcpy(mUpdataStruct.w_b, mGnssInsSystem.mPar.w_b, 3 * sizeof(double));
	memcpy(mUpdataStruct.C_bn, mGnssInsSystem.mNav.c_bn, 3 * 3 * sizeof(double));
	memset(*mUpdataStruct.Q, 0, 16 * 16 * sizeof(double));
	memset(*mUpdataStruct.PHI, 0, 16 * 16 * sizeof(double));

	memset(mUpdataStruct.X, 0, 9 * sizeof(double)); //???X????

	if (1)
	{
		memset(mUpdataStruct.X + 9, 0, 7 * sizeof(double)); //?????t?????????
	}

	for (int8_t i = 0; i < 16; i++)
	{
		for (int8_t j = 0; j < 16; j++)
		{
			if (i == j)
			{
				mUpdataStruct.PHI[i][j] = 1;
				continue;
			}
		}
	}
	mGnssInsSystem.flag = 1;
	return 1;
}

static int8_t ConvertIMUData(const ImuData *imudata, ImuData *SySimudata)
{
	int8_t ret = 1;
	double dt = imudata->timestamp - mGnssInsSystem.mpreImuData.timestamp;
	if (mGnssInsSystem.mpreImuData.timestamp < 0.01)
	{
		dt = 0.01;
	}
	if (dt < 0.0001)
	{
		ret = -1;
		return ret;
	}
	ret = 1;

	SySimudata->timestamp = imudata->timestamp;
	SySimudata->timestamped = imudata->timestamped;
	SySimudata->Accx = imudata->Accx;
	SySimudata->Accy = imudata->Accy;
	SySimudata->Accz = imudata->Accz;
	SySimudata->Gyrox = imudata->Gyrox;
	SySimudata->Gyroy = imudata->Gyroy;
	SySimudata->Gyroz = imudata->Gyroz;

	return ret;
}

static int8_t MechanizationINT(const ImuData* imudata)
{
	double DeltaT = 0.0;
	DeltaT = imudata->timestamp - mGnssInsSystem.mpreImuData.timestamp;

	if (DeltaT < 0.00001)
	{
		return -1;                                     //repeated data
	}
	//SetZuptDetectData(*imudata, 50);
	DataChangeLC(DeltaT, imudata, &mGnssInsSystem.mInsData);
	compensate(DeltaT, &mGnssInsSystem.mNav.sensorbias, &mGnssInsSystem.mInsData);

	if (IsFinished)  //???????
	{
		MatrixMutiply(*mUpdataStruct.PHI, mUpdataStruct.X, 16, 16, 1, mGnssInsSystem.mKalmanStruct.X); //PHI*X 
		//memcpy(mGnssInsSystem.mKalmanStruct.X, mUpdataStruct.X,16*sizeof(double));
		double PHIP[16][16], PHIt[16][16];
		MatrixMutiply(*mUpdataStruct.PHI, *mUpdataStruct.P, 16, 16, 16, *PHIP);
		MatrixTranspose(*mUpdataStruct.PHI, 16, 16, *PHIt);
		MatrixMutiply(*PHIP, *PHIt, 16, 16, 16, *(mGnssInsSystem.mKalmanStruct.P));
		MatrixAdd(*(mGnssInsSystem.mKalmanStruct.P), *mUpdataStruct.Q, 16, 16, *(mGnssInsSystem.mKalmanStruct.P));

		int16_t feedback = 1;
		if (imudata->timestamp - mGnssInsSystem.AlignCTime < 3)feedback = 0;
		KF_feedback(&(mGnssInsSystem.mPar), mGnssInsSystem.mKalmanStruct.X, &(mGnssInsSystem.mNav), feedback, &mGnssInsSystem);
		IsFinished = 0;
	}

	memcpy(&(mGnssInsSystem.mpreNav), &(mGnssInsSystem.mNav), sizeof(Nav));

	int8_t type = IsNeedLC(mGnssInsSystem.mpreInsData.timestamp, mGnssInsSystem.mInsData.timestamp, mGnssInsSystem.nextLCTime);

	double phi[16][16];

	if (type == 1)
	{
		INS_MECH(&mGnssInsSystem.mpreInsData, &mGnssInsSystem.mInsData, &mGnssInsSystem.mpreNav, &mGnssInsSystem.mNav, &mGnssInsSystem.mPar);

		KF_predict_16PHI(DeltaT, &mGnssInsSystem.mNav, &mGnssInsSystem.mPar, &mGnssInsSystem.mIMUSensor, phi);
		KF_predict(phi, mGnssInsSystem.Q, DeltaT, mGnssInsSystem.mKalmanStruct.X, mGnssInsSystem.mKalmanStruct.P , &mUpdataStruct);

		GetLCData(&mGnssInsSystem.mNav, mGnssInsSystem.nextLCTime);
	}
	else if (type == 2)
	{
		GetLCData(&mGnssInsSystem.mNav, mGnssInsSystem.nextLCTime);
		INS_MECH(&mGnssInsSystem.mpreInsData, &mGnssInsSystem.mInsData, &mGnssInsSystem.mpreNav, &mGnssInsSystem.mNav, &mGnssInsSystem.mPar);

		double PHI[16][16] = { 0.0 };
		KF_predict_16PHI(DeltaT, &mGnssInsSystem.mNav, &mGnssInsSystem.mPar, &mGnssInsSystem.mIMUSensor, phi);
		KF_predict(phi, mGnssInsSystem.Q, DeltaT, mGnssInsSystem.mKalmanStruct.X, mGnssInsSystem.mKalmanStruct.P, &mUpdataStruct);
	}
	else if (type == 3)
	{
		INS_MECH(&mGnssInsSystem.mpreInsData, &mGnssInsSystem.mInsData, &mGnssInsSystem.mpreNav, &mGnssInsSystem.mNav, &mGnssInsSystem.mPar);
		KF_predict_16PHI(DeltaT, &mGnssInsSystem.mNav, &mGnssInsSystem.mPar, &mGnssInsSystem.mIMUSensor, phi);
		KF_predict(phi, mGnssInsSystem.Q, DeltaT, mGnssInsSystem.mKalmanStruct.X, mGnssInsSystem.mKalmanStruct.P, &mUpdataStruct);

		double dt = imudata->timestamp - mGnssInsSystem.nextLCTime;
		dt = (imudata->timestamp - mGnssInsSystem.nextLCTime) / DeltaT;
		mGnssInsSystem.mpreNav.lat = (1 - dt)*mGnssInsSystem.mNav.lat + mGnssInsSystem.mpreNav.lat *dt;
		mGnssInsSystem.mpreNav.lon = (1 - dt)*mGnssInsSystem.mNav.lon + mGnssInsSystem.mpreNav.lon *dt;
		mGnssInsSystem.mpreNav.height = (1 - dt)*mGnssInsSystem.mNav.height + mGnssInsSystem.mpreNav.height *dt;

		GetLCData(&mGnssInsSystem.mpreNav, mGnssInsSystem.nextLCTime);
	}
	else
	{
		INS_MECH(&mGnssInsSystem.mpreInsData, &mGnssInsSystem.mInsData, &mGnssInsSystem.mpreNav, &mGnssInsSystem.mNav, &mGnssInsSystem.mPar);
		KF_predict_16PHI(DeltaT, &mGnssInsSystem.mNav, &mGnssInsSystem.mPar, &mGnssInsSystem.mIMUSensor, phi);
		KF_predict(phi, mGnssInsSystem.Q, DeltaT, mGnssInsSystem.mKalmanStruct.X, mGnssInsSystem.mKalmanStruct.P, &mUpdataStruct);
	}

	memcpy(&mGnssInsSystem.mpreInsData, &mGnssInsSystem.mInsData, sizeof(INS));

	return 1;
}

static int8_t firstfusionimu = 1;
int8_t AddIMUData(const ImuData mImudata, FILE *fSOL, FILE *fLOG)
{
	int8_t retal = -1;
	retal = ConvertIMUData(&mImudata, &mGnssInsSystem.mImuData);
	if (retal == -1)
	{
		return retal;
	}
	switch (mGnssInsSystem.mlc_STATUS)
	{
	case INSDATA_NREADY:
	{
		int8_t IsImuNormal = 1;// CheckImuDataBase(mGnssInsSystem.mImuData);
		if (IsImuNormal)
		{
			mGnssInsSystem.mlc_STATUS = INSDATA_READY;
			//std::cout << "IMU Data is Ready :" << setw(15) << mGnssInsSystem.mImuData.timestamp << "(s)" << std::endl;
		}
	}break;
	case INSDATA_READY:
	{
	}break;
	case ALIGNMENT_COMPLETE:
	{

	}break;
	case INS_FUSING:
	{

		if (firstfusionimu)
		{
			double DeltaT = mGnssInsSystem.mImuData.timestamp - mGnssInsSystem.mpreImuData.timestamp;
			DataChangeLC(DeltaT, &mGnssInsSystem.mpreImuData, &mGnssInsSystem.mpreInsData);
			firstfusionimu = 0;
		}
		if (mImudata.timestamped > 1)
		{
			SetLCtime(mGnssInsSystem.mImuData.timestamped);
		}

	    MechanizationINT(&mGnssInsSystem.mImuData);

		mGnssInsSystem.mMeasUpdataType = MeasNone;
		mGnssInsSystem.PerMeasUpdata = 0;


		ouputfile(mGnssInsSystem, fSOL, fLOG);

	}break;
	default:
		return  retal;
	}
	memcpy(&(mGnssInsSystem.mpreImuData), &(mGnssInsSystem.mImuData), sizeof(mGnssInsSystem.mImuData));
	retal = 1;

	return  retal;
}

/*--------------------------------------------------------*/
int8_t initsystemfromcfg(GnssInsSystem* pGnssInsSystem)
{
	int ret = -1;
	pGnssInsSystem->mNav.sensorbias.bias_acc_x = 0.0;//mLCSetting.a_bias[0];
	pGnssInsSystem->mNav.sensorbias.bias_acc_y = 0.0;// mLCSetting.a_bias[1];
	pGnssInsSystem->mNav.sensorbias.bias_acc_z = 0.0;// mLCSetting.a_bias[2];
	pGnssInsSystem->mNav.sensorbias.bias_gyro_x = 0.0 * PI / 180;// mLCSetting.g_bias[0] * PI / 180;
	pGnssInsSystem->mNav.sensorbias.bias_gyro_y = 0.0 * PI / 180;// mLCSetting.g_bias[1] * PI / 180;
	pGnssInsSystem->mNav.sensorbias.bias_gyro_z = 0.0;// mLCSetting.g_bias[2] * PI / 180;
	pGnssInsSystem->mNav.sensorbias.Odo_scale = 1.0;// mLCSetting.Odo_scale;


	for (uint8_t i = 0; i < 3; ++i)
	{
		pGnssInsSystem->mIMUSensor.Gyr_bias_std[i] = 2 * PI / 180 / 3600;     //
		pGnssInsSystem->mIMUSensor.Acc_bias_std[i] = 20 * 1.0e-5;
	}
	pGnssInsSystem->mIMUSensor.Gyr_noise_arw = 0.2 * (PI / 180 / 60);
	pGnssInsSystem->mIMUSensor.Acc_noise_vrw = 0.04 / 60;
	pGnssInsSystem->mIMUSensor.Gyr_bias_CorTime = 0.5 * 3600;
	pGnssInsSystem->mIMUSensor.Acc_bias_CorTime = 0.5 * 3600;
	pGnssInsSystem->mIMUSensor.bg_model[0] = exp(-1.0 / pGnssInsSystem->mIMUSensor.Gyr_bias_CorTime / mLCSetting.IMUDataRate);
	pGnssInsSystem->mIMUSensor.bg_model[1] = exp(-1.0 / pGnssInsSystem->mIMUSensor.Gyr_bias_CorTime / mLCSetting.IMUDataRate);
	pGnssInsSystem->mIMUSensor.bg_model[2] = exp(-1.0 / pGnssInsSystem->mIMUSensor.Gyr_bias_CorTime / mLCSetting.IMUDataRate);
	pGnssInsSystem->mIMUSensor.ba_model[0] = exp(-1.0 / pGnssInsSystem->mIMUSensor.Acc_bias_CorTime / mLCSetting.IMUDataRate);
	pGnssInsSystem->mIMUSensor.ba_model[1] = exp(-1.0 / pGnssInsSystem->mIMUSensor.Acc_bias_CorTime / mLCSetting.IMUDataRate);
	pGnssInsSystem->mIMUSensor.ba_model[2] = exp(-1.0 / pGnssInsSystem->mIMUSensor.Acc_bias_CorTime / mLCSetting.IMUDataRate);

	memset(&pGnssInsSystem->Q, 0, 16 * sizeof(double));
	pGnssInsSystem->Q[0] = 0.0;
	pGnssInsSystem->Q[1] = pGnssInsSystem->Q[0];
	pGnssInsSystem->Q[2] = pGnssInsSystem->Q[0];

	pGnssInsSystem->Q[3] = pGnssInsSystem->mIMUSensor.Acc_noise_vrw * pGnssInsSystem->mIMUSensor.Acc_noise_vrw;
	pGnssInsSystem->Q[4] = pGnssInsSystem->Q[3];
	pGnssInsSystem->Q[5] = pGnssInsSystem->Q[3];

	pGnssInsSystem->Q[6] = pGnssInsSystem->mIMUSensor.Gyr_noise_arw * pGnssInsSystem->mIMUSensor.Gyr_noise_arw;
	pGnssInsSystem->Q[7] = pGnssInsSystem->Q[6];
	pGnssInsSystem->Q[7] = pGnssInsSystem->Q[6];

	pGnssInsSystem->Q[9] = 2 * pGnssInsSystem->mIMUSensor.Gyr_bias_std[0] * pGnssInsSystem->mIMUSensor.Gyr_bias_std[0] / pGnssInsSystem->mIMUSensor.Gyr_bias_CorTime;
	pGnssInsSystem->Q[10] = 2 * pGnssInsSystem->mIMUSensor.Gyr_bias_std[1] * pGnssInsSystem->mIMUSensor.Gyr_bias_std[1] / pGnssInsSystem->mIMUSensor.Gyr_bias_CorTime;
	pGnssInsSystem->Q[11] = 2 * pGnssInsSystem->mIMUSensor.Gyr_bias_std[2] * pGnssInsSystem->mIMUSensor.Gyr_bias_std[2] / pGnssInsSystem->mIMUSensor.Gyr_bias_CorTime;

	pGnssInsSystem->Q[12] = 2 * pGnssInsSystem->mIMUSensor.Acc_bias_std[0] * pGnssInsSystem->mIMUSensor.Acc_bias_std[0] / pGnssInsSystem->mIMUSensor.Acc_bias_CorTime;
	pGnssInsSystem->Q[13] = 2 * pGnssInsSystem->mIMUSensor.Acc_bias_std[1] * pGnssInsSystem->mIMUSensor.Acc_bias_std[1] / pGnssInsSystem->mIMUSensor.Acc_bias_CorTime;
	pGnssInsSystem->Q[14] = 2 * pGnssInsSystem->mIMUSensor.Acc_bias_std[2] * pGnssInsSystem->mIMUSensor.Acc_bias_std[2] / pGnssInsSystem->mIMUSensor.Acc_bias_CorTime;

	pGnssInsSystem->Q[15] = 0.0001;

	pGnssInsSystem->premGnssData.timestamp = -999;



	return 1;

}

void acc2att(double* accel, double* roll, double* pitch)
{
	double accelMag = sqrt(accel[0] * accel[0] + accel[1] * accel[1] + accel[2] * accel[2]);
	double invAccelMag = 1.0f / accelMag;
	double unitGravityVector[3] = { 0.0 };
	unitGravityVector[0] = -accel[0] * invAccelMag;
	unitGravityVector[1] = -accel[1] * invAccelMag;
	unitGravityVector[2] = -accel[2] * invAccelMag;
	*roll = (atan2(unitGravityVector[1], unitGravityVector[2]));
	*pitch = (-asin(unitGravityVector[0]));
	return;
}

int8_t initsystemfromGNSS(GnssInsSystem* pGnssInsSystem)
{
	int8_t ret = -1;
	double att[3] = { 0.0 };
	//pGnssInsSystem->nextLCTime = floor(pGnssInsSystem->mGnssData.timestamp) + 1;
	pGnssInsSystem->GNSSflag = 0;

	double fxyz[3] = { pGnssInsSystem->mImuData.Accx, pGnssInsSystem->mImuData.Accy, pGnssInsSystem->mImuData.Accz };

	acc2att(fxyz, att + 0, att + 1);

	mLCSetting.INSRotationRvb[0] = floor(att[0] / PI + 0.5) * PI;
	mLCSetting.INSRotationRvb[1] = floor(att[1] / PI + 0.5) * PI;

	double C_bv[3][3] = { 0.0 };
	euler2dcm(mLCSetting.INSRotationRvb, C_bv);
	MatrixTranspose(*C_bv, 3, 3, *pGnssInsSystem->C_InstallationAngle);

	double eular[3] = { pGnssInsSystem->mNav.roll, pGnssInsSystem->mNav.pitch,pGnssInsSystem->mNav.heading };
	double C_vn[3][3] = { 0.0 };
	euler2dcm(eular, C_vn);
	MatrixMutiply(*C_vn, *C_bv, 3, 3, 3, *pGnssInsSystem->mNav.c_bn);
	dcm2euler(pGnssInsSystem->mNav.c_bn, eular);
	pGnssInsSystem->mNav.roll = eular[0];
	pGnssInsSystem->mNav.pitch = eular[1];
	pGnssInsSystem->mNav.heading = eular[2];
	euler2quat(eular, pGnssInsSystem->mNav.q_bn);

	double leverarm_b[3] = { 0.0 };
	MatrixMutiply(*pGnssInsSystem->C_InstallationAngle, mLCSetting.Leverarm, 3, 3, 1, leverarm_b);

	double leverarm_n[3] = { 0.0 };
	MatrixMutiply(*pGnssInsSystem->mNav.c_bn, leverarm_b, 3, 3, 1, leverarm_n);

	double pos[3] = { pGnssInsSystem->mGnssData.latitude, pGnssInsSystem->mGnssData.longitude, pGnssInsSystem->mGnssData.altitude };
	double M = 0.0, N = 0.0;
	UpdateMN(pos, &M, &N);

	double d_leverarm[3] = { 0.0 };
	d_leverarm[0] = leverarm_n[0] / (M + pos[2]);
	d_leverarm[1] = leverarm_n[1] / ((N + pos[2]) * cos(pos[0]));
	d_leverarm[2] = -leverarm_n[2];

	MatrixSub(pos, d_leverarm, 3, 1, pos);

	pGnssInsSystem->mNav.lat = pos[0] + (pGnssInsSystem->mGnssData.timestampd - pGnssInsSystem->mGnssData.timestamp) * pGnssInsSystem->mNav.vn / (M + pos[2]);
	pGnssInsSystem->mNav.lon = pos[1] + (pGnssInsSystem->mGnssData.timestampd - pGnssInsSystem->mGnssData.timestamp) * pGnssInsSystem->mNav.ve / ((N + pos[2]) * cos(pos[0]));
	pGnssInsSystem->mNav.height = pos[2];

	pos2quat(pos[0], pos[1], pGnssInsSystem->mNav.q_ne);

	memset(&pGnssInsSystem->mKalmanStruct.P, 0, 16 * 16 * sizeof(double));
	double roll_std = 5 * PI / 180;
	double pitch_std = 5 * PI / 180;
	double heading_std = 10 * PI / 180;

	int scale = 10 * 10;
	pGnssInsSystem->mKalmanStruct.P[0][0] = 10;//0 * pGnssInsSystem->mGnssData.latitude_std * pGnssInsSystem->mGnssData.latitude_std;
	pGnssInsSystem->mKalmanStruct.P[1][1] = 10;// 0 * pGnssInsSystem->mGnssData.longitude_std * pGnssInsSystem->mGnssData.longitude_std;
	pGnssInsSystem->mKalmanStruct.P[2][2] = 10;// 0 * pGnssInsSystem->mGnssData.altitude_std * pGnssInsSystem->mGnssData.altitude_std;

	pGnssInsSystem->mKalmanStruct.P[3][3] = 1;// scale * 0.01*0.01;//pGnssInsSystem->mGnssData.north_velocity_std * pGnssInsSystem->mGnssData.north_velocity_std;
	pGnssInsSystem->mKalmanStruct.P[4][4] = 1;//scale * 0.01*0.01;//pGnssInsSystem->mGnssData.east_velocity_std * pGnssInsSystem->mGnssData.east_velocity_std;
	pGnssInsSystem->mKalmanStruct.P[5][5] = 1;//scale * 0.01*0.01;//pGnssInsSystem->mGnssData.up_velocity_std * pGnssInsSystem->mGnssData.up_velocity_std;

	pGnssInsSystem->mKalmanStruct.P[6][6] = roll_std * roll_std;
	pGnssInsSystem->mKalmanStruct.P[7][7] = pitch_std * pitch_std;
	pGnssInsSystem->mKalmanStruct.P[8][8] = heading_std * heading_std;

	pGnssInsSystem->mKalmanStruct.P[9][9] = (0.02 * (PI / 180) / 3600) * (0.02 * (PI / 180) / 3600);//10000*pGnssInsSystem->mIMUSensor.Gyr_bias_std[0] * pGnssInsSystem->mIMUSensor.Gyr_bias_std[0];
	pGnssInsSystem->mKalmanStruct.P[10][10] = (0.02 * (PI / 180) / 3600) * (0.02 * (PI / 180) / 3600);//10000*pGnssInsSystem->mIMUSensor.Gyr_bias_std[1] * pGnssInsSystem->mIMUSensor.Gyr_bias_std[1];
	pGnssInsSystem->mKalmanStruct.P[11][11] = (0.02 * (PI / 180) / 3600) * (0.02 * (PI / 180) / 3600);//10000*pGnssInsSystem->mIMUSensor.Gyr_bias_std[2] * pGnssInsSystem->mIMUSensor.Gyr_bias_std[2];

	pGnssInsSystem->mKalmanStruct.P[12][12] = (0.02) * (0.02);//10000*pGnssInsSystem->mIMUSensor.Acc_bias_std[0] * pGnssInsSystem->mIMUSensor.Acc_bias_std[0];
	pGnssInsSystem->mKalmanStruct.P[13][13] = (0.02) * (0.02);//10000*pGnssInsSystem->mIMUSensor.Acc_bias_std[1] * pGnssInsSystem->mIMUSensor.Acc_bias_std[1];
	pGnssInsSystem->mKalmanStruct.P[14][14] = (0.02) * (0.02); //10000*pGnssInsSystem->mIMUSensor.Acc_bias_std[2] * pGnssInsSystem->mIMUSensor.Acc_bias_std[2];


	pGnssInsSystem->mKalmanStruct.P[15][15] = 0.1 * 0.1;

	return 1;

}
int8_t initsystemSoftreset()
{
	int ret = 0;
	double DeltaT = 0;
	mGnssInsSystem.CurIsUseZupt = 0;

	if (mGnssInsSystem.InitNavi == -1)
	{
		mGnssInsSystem.InitNavi = -1;
		mGnssInsSystem.InsIsTimeOut = -1;
		mGnssInsSystem.KFIsConverage = -1;
		mGnssInsSystem.TimeOutZuptTime = 0.0;
		mGnssInsSystem.NaviStartTime = mGnssInsSystem.mGnssData.timestampd;
		mGnssInsSystem.TimeOutStartTime = -999.999;
		mGnssInsSystem.LastBdsRejectedTime = -999.999;
		mGnssInsSystem.LastBdsUpdateTime = -999.999;
		mGnssInsSystem.LastInsUpdateTime = -999.999;
		mGnssInsSystem.GyroBiasIsConverage = -1;
		mGnssInsSystem.InsZuptTimeLen = -999.999;
	}
	return 1;
}

int8_t SetLCtime(double lctime)
{
	mGnssInsSystem.nextLCTime = lctime;
	return 1;
}

/*--------------------------------------------------------*/

int8_t readconfigfromfile(const char* cfgfname, LCSetting*  mLCSetting)
{
	int8_t ret = -1;
	FILE* fcfg = NULL;
	size_t len;
	uint8_t* data = NULL;

    fopen_s(&fcfg,cfgfname, "r");
	if (fcfg == NULL)
	{
		return ret;
	}
	/* get the length */
	fseek(fcfg, 0, SEEK_END);
	len = ftell(fcfg);
	fseek(fcfg, 0, SEEK_SET);

	data = (uint8_t*)malloc(len + 1);
	memset(data, 0, len + 1);
	fread(data, 1, len, fcfg);
	data[len] = '\0';
	fclose(fcfg);

	// parse data
	do {
		#if 0
			/* replace JSON decoder: TODO */
		cJSON *root = NULL;

		root = cJSON_Parse(data);
		if (!root)
		{
			break;
		}
		mLCSetting->ProcType = cJSON_GetObjectItem(root, "PROCESS")->valueint;
		char* filename = cJSON_GetObjectItem(root, "InsfileName")->valuestring;
		memcpy(mLCSetting->InsfileName, filename, strlen(filename) * sizeof (char));

		mLCSetting->mIMUSenorType = cJSON_GetObjectItem(root, "IMUSenorType")->valueint;
		mLCSetting->IMUDataRate = cJSON_GetObjectItem(root, "IMUDataRate")->valueint;
		mLCSetting->ODODataRate = cJSON_GetObjectItem(root, "ODODataRate")->valueint;


		cJSON * leverarmArrayItem = cJSON_GetObjectItem(root, "Leverarm");
		for (int i = 0; i < 3; i++)
		{
			mLCSetting->Leverarm[i] = cJSON_GetArrayItem(leverarmArrayItem, i)->valuedouble;
		}

		cJSON * odoLevararmArrayItem = cJSON_GetObjectItem(root, "OdoLeverarm");
		for (int i = 0; i < 3; i++)
		{
			mLCSetting->OdoLeverarm[i] = cJSON_GetArrayItem(odoLevararmArrayItem, i)->valuedouble;
		}

		cJSON * INSRotationRvbArrayItem = cJSON_GetObjectItem(root, "INSRotationRvb");
		for (int i = 0; i < 3; i++)
		{
			mLCSetting->INSRotationRvb[i] = 0.017453292519943* cJSON_GetArrayItem(INSRotationRvbArrayItem, i)->valuedouble;
		}

		mLCSetting->IsUser = cJSON_GetObjectItem(root, "IsUser")->valueint;

		cJSON * UserLeverarmArrayItem = cJSON_GetObjectItem(root, "UserLeverarm");
		for (int i = 0; i < 3; i++)
		{
			mLCSetting->UserLeverarm[i] = cJSON_GetArrayItem(UserLeverarmArrayItem, i)->valuedouble;
		}

		mLCSetting->IsOutputKML = cJSON_GetObjectItem(root, "IsOutputKML")->valueint;

		mLCSetting->GnssDataType = cJSON_GetObjectItem(root, "GnssDataType")->valueint;
		mLCSetting->InitialAttitudeMode = cJSON_GetObjectItem(root, "InitialAttitudeMode")->valueint;


		mLCSetting->IsUseNHC = cJSON_GetObjectItem(root, "IsUseNHC")->valueint;
		mLCSetting->IsUseOdo = cJSON_GetObjectItem(root, "IsUseOdo")->valueint;
		mLCSetting->IsUseZUPT = cJSON_GetObjectItem(root, "IsUseZUPT")->valueint;



		cJSON * Attitue_RPHArrayItem = cJSON_GetObjectItem(root, "Attitue_RPH");
		for (int i = 0; i < 3; i++)
		{
			mLCSetting->Attitue_RPH[i] = cJSON_GetArrayItem(Attitue_RPHArrayItem, i)->valuedouble;
		}

		cJSON * a_biasArrayItem = cJSON_GetObjectItem(root, "a_bias");
		for (int i = 0; i < 3; i++)
		{
			mLCSetting->a_bias[i] = cJSON_GetArrayItem(a_biasArrayItem, i)->valuedouble;
		}

		cJSON * g_biasArrayItem = cJSON_GetObjectItem(root, "g_bias");
		for (int i = 0; i < 3; i++)
		{
			mLCSetting->g_bias[i] = cJSON_GetArrayItem(g_biasArrayItem, i)->valuedouble;
		}

		cJSON * a_scaleArrayItem = cJSON_GetObjectItem(root, "a_scale");
		for (int i = 0; i < 3; i++)
		{
			mLCSetting->a_scale[i] = cJSON_GetArrayItem(a_scaleArrayItem, i)->valuedouble;
		}

		cJSON * g_scaleArrayItem = cJSON_GetObjectItem(root, "g_scale");
		for (int i = 0; i < 3; i++)
		{
			mLCSetting->g_scale[i] = cJSON_GetArrayItem(g_scaleArrayItem, i)->valuedouble;
		}
		mLCSetting->Odo_scale = cJSON_GetObjectItem(root, "Odo_scale")->valuedouble;

		ret = 1;
		cJSON_Delete(root);
		#endif
	} while (0);

	free(data);
	data = NULL;
	return ret;
		 
}

#define MAXFIELD 100

int8_t  parse_fields(char* const buffer, char** val, int* size)
{
	int8_t ret = -1;
	char* p, * q;
	int n = 0;

	/* parse fields */
	for (p = buffer; *p && n < MAXFIELD; p = q + 1) {
		if ((q = strchr(p, ',')) || (q = strchr(p, '*'))) {
			val[n++] = p;
			*q = '\0';
		}
		else
		{
			val[n++] = p;
			break;
		}

	}
	*size = n;
	return 1;

}



int8_t playLogFile(const char *cfgfname)
{
	int8_t ret = -1;
	FILE *fdat = NULL;
	char buffer[150] = { 0 };
	char  *val[20];
	int32_t sentencesize;
	int8_t gnsflag = -1;
	double delytime = 0;

	LCSetting mLCSetting = { 0 };

	readconfigfromfile(cfgfname, &mLCSetting);

	initsystemfromcfg(&mGnssInsSystem);

	fopen_s(&fdat,mLCSetting.InsfileName, "r");
	if (!fdat)
	{
		return ret;
	}

	ImuData mImuData = { 0 };
	GnssData mGnssData = { 0 };

	double LastIMUDataTime = 0.0;

	while (!feof(fdat))
	{
		fgets(buffer, sizeof(buffer), fdat);
		if (strlen(buffer) <= 0) continue;	
		char *temp = strchr(buffer, '\n');
		if (temp != NULL) temp[0] = '\0';

		parse_fields(buffer, val,&sentencesize);

		if (strstr(val[0], "$GPIMU") != NULL)
		{
			if (sentencesize == 9)
			{
				memset(&mImuData, 0, sizeof(ImuData));
				mImuData.timestamp = atof(val[2]);  
				mImuData.timestamped = atof(val[1]);

			    double dt = 0.01;
				if (fabs(LastIMUDataTime) > 1)
				{
					dt = mImuData.timestamp - LastIMUDataTime;
				}

				mImuData.Accx = atof(val[3])* 0.02;
				mImuData.Accy = atof(val[4])* 0.02;
				mImuData.Accz = atof(val[5])* 0.02;
				mImuData.Gyrox = atof(val[6])* 0.02 * PI / 180;
				mImuData.Gyroy = atof(val[7])* 0.02* PI / 180;
				mImuData.Gyroz = atof(val[8])* 0.02 * PI / 180;
				if (fabs(mImuData.Accx) > 6 * 0.02
					|| fabs(mImuData.Accy) > 6 * 0.02
					|| fabs(mImuData.Accz + 9.8 * 0.02) > 3 * 0.02
					|| fabs(mImuData.Gyrox) > 5 * 0.02  * PI / 180
					|| fabs(mImuData.Gyroy) > 5 * 0.02  * PI / 180
					)
				{
					// "Imudata error" 
					continue;
				}
				AddIMUData(mImuData, NULL, NULL);

			}
			else
			{
				//"IMU Data error,line"
				continue;
			}

		}
		if (strstr(val[0], "$GPGNSS") != NULL)
		{
			if (sentencesize == 10)
			{
				mGnssData.timestampd = mImuData.timestamp;//atof(val[1]);  
				mGnssData.timestamp = atof(val[2]);  

				mGnssData.longitude = atof(val[4])*PI / 180;
				mGnssData.latitude = atof(val[3])*PI / 180;
				mGnssData.altitude = atof(val[5]);
				mGnssData.longitude_std = atof(val[7]);
				mGnssData.latitude_std = atof(val[6]);
				mGnssData.altitude_std = atof(val[8]);
				mGnssData.Mode = atof(val[9]);
				//if (mGnssData.Mode == 5)
				//{
				//	mGnssData.longitude_std *= 2;
				//	mGnssData.latitude *= 2;
				//	mGnssData.altitude_std *= 2;

				//}
				mGnssData.Mode = 4;
				//	atoi(val[10]);
				//mGnssData.heading = atoi(val[14]) * PI / 180;
				//if (mGnssData.timestamp > 1557846910 && mGnssData.timestamp < 1557846970)
				//{
				//	mGnssData.flag = 0;
				//	ADDGNSSDATA(mGnssData);
				//}
				//else
				//{

				//	gnsflag = 1;
				//}
				gnsflag = 1;

			}
			else
			{
				// "Gnss Data error,line" 
			}
			if (gnsflag)
			{
				mGnssData.flag = 1;
				//if (mImuData.timestamp - mGnssData.timestamp > delytime)  //ADD delytime 
				//{
					ADDGNSSDATA(mGnssData);
					gnsflag = 0;
				/*}*/
			}
		}


		}

	if (fdat)fclose(fdat);

	return 1;

}

int8_t process_adi16485_imu(const char* fname)
{
	int8_t ret = -1;
	FILE* fdat = fopen(fname, "r");
	FILE* fSOL = fopen("sol.txt", "w");
	FILE* fLOG = fopen("log.txt", "w");
	char buffer[1024] = { 0 };
	char* val[20];

	mLCSetting.GnssDataType = 1;
	mLCSetting.IMUDataRate = 100;

	initsystemfromcfg(&mGnssInsSystem);

	ImuData mImuData = { 0 };
	GnssData mGnssData = { 0 };

	double data[50] = { 0.0 };
	double data_[50] = { 0.0 };

	double LastIMUDataTime = 0.0;
	unsigned long numofepoch = 0;

	while (fdat!=NULL&&!feof(fdat))
	{
		fgets(buffer, sizeof(buffer), fdat);
		if (strlen(buffer) <= 0) continue;
		char* temp = strchr(buffer, '\n');
		if (temp != NULL) temp[0] = '\0';

		sscanf(buffer, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf"
			, data + 0 /* time */
			, data + 1, data + 2, data + 3 /* fxyz */
			, data + 4, data + 5, data + 6 /* wxyz */
			, data + 7 /* time GPS */
			, data + 8, data + 9, data + 10 /* lat, lon, ht */
			, data + 11, data + 12, data + 13 /* vn, ve, vd */
			, data + 14, data + 15, data + 16 /* rms N, rms E, rms U */
			, data + 17, data + 18, data + 19 /* rms Vn, rms Ve, rms Vu */
			, data + 20, data + 21 /* pos type, vel type */
			);

		++numofepoch;
		if (numofepoch == 1)
		{
			memcpy(data_, data, sizeof(data));
			continue;
		}
		double dt = data[0] - data_[0];

		memset(&mImuData, 0, sizeof(ImuData));
		memset(&mGnssData, 0, sizeof(GnssData));

		mImuData.timestamp = data[0];
		mImuData.timestamped = data[0];

		mImuData.Accx = data[1];
		mImuData.Accy = data[2];
		mImuData.Accz = data[3];
		mImuData.Gyrox = data[4];
		mImuData.Gyroy = data[5];
		mImuData.Gyroz = data[6];

		AddIMUData(mImuData, fSOL, fLOG);

		double dt_gps = data[7] - data[0];
		if (fabs(dt_gps) < 0.005)
		{
			/* new GPS data */
			mGnssData.timestampd = data[0];  
			mGnssData.timestamp = data[7];

			mGnssData.latitude = data[8];
			mGnssData.longitude = data[9];
			mGnssData.altitude = data[10];
			mGnssData.north_velocity = data[11];
			mGnssData.east_velocity = data[12];
			mGnssData.up_velocity =-data[13];
			mGnssData.latitude_std = data[14];
			mGnssData.longitude_std = data[15];
			mGnssData.altitude_std = data[16];
			mGnssData.north_velocity_std = data[17];
			mGnssData.east_velocity_std = data[18];
			mGnssData.up_velocity_std = data[19];

			mGnssData.Mode = (data[20]==50)?(4):(5);

			mGnssData.flag = 1;

			ADDGNSSDATA(mGnssData);

		}
		else
		{
			dt_gps = data[0] - data_[7];
			if (fabs(dt_gps) < 0.005)
			{
				/* new GPS data */
				mGnssData.timestampd = data[0];
				mGnssData.timestamp = data_[7];

				mGnssData.latitude = data_[8];
				mGnssData.longitude = data_[9];
				mGnssData.altitude = data_[10];
				mGnssData.north_velocity = data_[11];
				mGnssData.east_velocity = data_[12];
				mGnssData.up_velocity = -data_[13];
				mGnssData.latitude_std = data_[14];
				mGnssData.longitude_std = data_[15];
				mGnssData.altitude_std = data_[16];
				mGnssData.north_velocity_std = data_[17];
				mGnssData.east_velocity_std = data_[18];
				mGnssData.up_velocity_std = data_[19];

				mGnssData.Mode = (data_[20] == 50) ? (4) : (5);

				mGnssData.flag = 1;

				ADDGNSSDATA(mGnssData);
			}
		}

		memcpy(data_, data, sizeof(data));

	}

	if (fdat)fclose(fdat);
	if (fSOL)fclose(fSOL);
	if (fLOG)fclose(fLOG);

	return 1;

}



int8_t ouputfile(GnssInsSystem mGnssInsSystem, FILE *output, FILE *output_debug)
{
	int8_t ret = -1;

	if (output)
	{
		fprintf(output, "%14.10f,%14.10f,%14.10f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%d,%d,%d,%10.5f\n",
			mGnssInsSystem.mInsData.timestamp, mGnssInsSystem.mNav.lat * 180 / PI, mGnssInsSystem.mNav.lon * 180 / PI, mGnssInsSystem.mNav.height,
			mGnssInsSystem.mNav.vn, mGnssInsSystem.mNav.ve, mGnssInsSystem.mNav.vd,
			mGnssInsSystem.mNav.roll * 180 / PI, mGnssInsSystem.mNav.pitch * 180 / PI, mGnssInsSystem.mNav.heading * 180 / PI,
			mGnssInsSystem.mNav.sensorbias.bias_gyro_x * 180 / PI * 3600, mGnssInsSystem.mNav.sensorbias.bias_gyro_y * 180 / PI * 3600, mGnssInsSystem.mNav.sensorbias.bias_gyro_z * 180 / PI * 3600,
			mGnssInsSystem.mNav.sensorbias.bias_acc_x * 1.0e5, mGnssInsSystem.mNav.sensorbias.bias_acc_y * 1.0e5, mGnssInsSystem.mNav.sensorbias.bias_acc_z * 1.0e5,
			mGnssInsSystem.mNav.sensorbias.Odo_scale,
			mGnssInsSystem.PerMeasUpdata,
			(int)mLCSetting.InitialAttitudeMode,
			(int)mGnssInsSystem.CurIsUseZupt,
			mGnssInsSystem.OdoV
		);
	}
	if (output_debug)
	{
		fprintf(output_debug, "%14.10f,%14.10f,%14.10f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f\n",
			mGnssInsSystem.mInsData.timestamp,
			sqrt(mGnssInsSystem.mKalmanStruct.P[0][0]), sqrt(mGnssInsSystem.mKalmanStruct.P[1][1]), sqrt(mGnssInsSystem.mKalmanStruct.P[2][2]),
			sqrt(mGnssInsSystem.mKalmanStruct.P[3][3]), sqrt(mGnssInsSystem.mKalmanStruct.P[4][4]), sqrt(mGnssInsSystem.mKalmanStruct.P[5][5]),
			sqrt(mGnssInsSystem.mKalmanStruct.P[6][6]) * 180 / PI, sqrt(mGnssInsSystem.mKalmanStruct.P[7][7]) * 180 / PI, sqrt(mGnssInsSystem.mKalmanStruct.P[8][8]) * 180 / PI,
			sqrt(mGnssInsSystem.mKalmanStruct.P[9][9]) * 180 / PI * 3600, sqrt(mGnssInsSystem.mKalmanStruct.P[10][10]) * 180 / PI * 3600, sqrt(mGnssInsSystem.mKalmanStruct.P[11][11]) * 180 / PI * 3600,
			sqrt(mGnssInsSystem.mKalmanStruct.P[12][12]) * 1.0e5, sqrt(mGnssInsSystem.mKalmanStruct.P[13][13]) * 1.0e5, sqrt(mGnssInsSystem.mKalmanStruct.P[14][14]) * 1.0e5,
			sqrt(mGnssInsSystem.mKalmanStruct.P[15][15])
		);
	}
	return 1;
}
