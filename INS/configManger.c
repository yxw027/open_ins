#include <string.h>
#include <stdlib.h>
#include "cjson/cJSON.h"
#include "configManger.h"

#ifndef DEG2RAD 
#define DEG2RAD (0.017453292519943)
#endif
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
		cJSON *root = NULL;

		root = cJSON_Parse(data);
		if (!root)
		{
			break;
		}
		mLCSetting->procType = cJSON_GetObjectItem(root, "procType")->valueint;
		char* procfileNme = cJSON_GetObjectItem(root, "procfileNme")->valuestring;
		memcpy(mLCSetting->procfileNme, procfileNme, strlen(procfileNme) * sizeof (char));

		char* gnssfileNme = cJSON_GetObjectItem(root, "gnssfileNme")->valuestring;
		memcpy(mLCSetting->gnssfileNme, gnssfileNme, strlen(gnssfileNme) * sizeof(char));


		char* insfileName = cJSON_GetObjectItem(root, "insfileName")->valuestring;
		memcpy(mLCSetting->insfileName, insfileName, strlen(insfileName) * sizeof(char));

		mLCSetting->initialAttitudeMode = cJSON_GetObjectItem(root, "initialAttitudeMode")->valueint;
		cJSON * attitueRPHArrayItem = cJSON_GetObjectItem(root, "attitueRPH");
		for (int i = 0; i < 3; i++)
		{
			mLCSetting->attitueRPH[i] = cJSON_GetArrayItem(attitueRPHArrayItem, i)->valuedouble;
		}

		mLCSetting->imuSensorType = cJSON_GetObjectItem(root, "imuSensorType")->valueint;
		mLCSetting->imuDataRate = cJSON_GetObjectItem(root, "imuDataRate")->valueint;

		mLCSetting->gnssSensorType = cJSON_GetObjectItem(root, "gnssSensorType")->valueint;
		mLCSetting->gnssDataType = cJSON_GetObjectItem(root, "gnssDataType")->valueint;
		mLCSetting->gnssDataRate = cJSON_GetObjectItem(root, "gnssDataRate")->valueint;
		mLCSetting->odoDataRate = cJSON_GetObjectItem(root, "odoDataRate")->valueint;


		cJSON * priLeverArmArrayItem = cJSON_GetObjectItem(root, "priLeverArm");
		for (int i = 0; i < 3; i++)
		{
			mLCSetting->priLeverArm[i] = cJSON_GetArrayItem(priLeverArmArrayItem, i)->valuedouble;
		}
		mLCSetting->isUseDuaAnt = cJSON_GetObjectItem(root, "isUseDuaAnt")->valueint;
		cJSON * secLeverArmArrayItem = cJSON_GetObjectItem(root, "secLeverArm");
		for (int i = 0; i < 3; i++)
		{
			mLCSetting->secLeverArm[i] = cJSON_GetArrayItem(secLeverArmArrayItem, i)->valuedouble;
		}

		cJSON * rotationRBVArrayItem = cJSON_GetObjectItem(root, "rotationRBV");
		for (int i = 0; i < 3; i++)
		{
			mLCSetting->rotationRBV[i] = DEG2RAD * cJSON_GetArrayItem(rotationRBVArrayItem, i)->valuedouble;
		}

		mLCSetting->isUseMisAlignment = cJSON_GetObjectItem(root, "isUseMisAlignment")->valueint;
		cJSON * misAlignmentAiaxArrayItem = cJSON_GetObjectItem(root, "MisAlignmentAiax");
		for (int i = 0; i < 3; i++)
		{
			mLCSetting->misAlignment[i] = cJSON_GetArrayItem(misAlignmentAiaxArrayItem, i)->valueint;
		}
		cJSON * misAlignmentArrayItem = cJSON_GetObjectItem(root, "MisAlignment");
		for (int i = 0; i < 3; i++)
		{
			mLCSetting->misAlignment[i] = DEG2RAD * cJSON_GetArrayItem(misAlignmentArrayItem, i)->valuedouble;
		}
		mLCSetting->isOnlineMisAlignmentEst = cJSON_GetObjectItem(root, "isOnlineMisAlignmentEst")->valueint;


		mLCSetting->isUseNHC = cJSON_GetObjectItem(root, "isUseNHC")->valueint;
		cJSON * odoLeverarmArrayItem = cJSON_GetObjectItem(root, "odoLeverarm");
		for (int i = 0; i < 3; i++)
		{
			mLCSetting->odoLeverArm[i] = cJSON_GetArrayItem(odoLeverarmArrayItem, i)->valuedouble;
		}
		mLCSetting->isUseGNSSVel = cJSON_GetObjectItem(root, "isUseGNSSVel")->valueint;
		mLCSetting->isUseOdo = cJSON_GetObjectItem(root, "isUseOdo")->valueint;
		mLCSetting->odoScale = cJSON_GetObjectItem(root, "odoScale")->valuedouble;

		mLCSetting->isUseZUPT = cJSON_GetObjectItem(root, "isUseZUPT")->valueint;
		mLCSetting->isUseZUPTLOCK = cJSON_GetObjectItem(root, "isUseZUPTLOCK")->valueint;


		cJSON * a_biasArrayItem = cJSON_GetObjectItem(root, "accBias");
		for (int i = 0; i < 3; i++)
		{
			mLCSetting->accBias[i] = cJSON_GetArrayItem(a_biasArrayItem, i)->valuedouble;
		}

		cJSON * g_biasArrayItem = cJSON_GetObjectItem(root, "gyroBias");
		for (int i = 0; i < 3; i++)
		{
			mLCSetting->gyroBias[i] = cJSON_GetArrayItem(g_biasArrayItem, i)->valuedouble;
		}

		cJSON * a_scaleArrayItem = cJSON_GetObjectItem(root, "accScale");
		for (int i = 0; i < 3; i++)
		{
			mLCSetting->accScale[i] = cJSON_GetArrayItem(a_scaleArrayItem, i)->valuedouble;
		}

		cJSON * g_scaleArrayItem = cJSON_GetObjectItem(root, "gyroScale");
		for (int i = 0; i < 3; i++)
		{
			mLCSetting->gyroScale[i] = cJSON_GetArrayItem(g_scaleArrayItem, i)->valuedouble;
		}

		mLCSetting->isUserOutPut = cJSON_GetObjectItem(root, "isUserOutPut")->valueint;
		cJSON * userLeverArmArrayItem = cJSON_GetObjectItem(root, "userLeverArm");
		for (int i = 0; i < 3; i++)
		{
			mLCSetting->userLeverArm[i] = cJSON_GetArrayItem(userLeverArmArrayItem, i)->valuedouble;
		}

		mLCSetting->fileOutputType = cJSON_GetObjectItem(root, "fileOutputType")->valueint;
		mLCSetting->profiletype = cJSON_GetObjectItem(root, "profiletype")->valueint;
		mLCSetting->isKmlOutput = cJSON_GetObjectItem(root, "isKmlOutput")->valueint;
		mLCSetting->kmlOutputDateRate = cJSON_GetObjectItem(root, "kmlOutputDateRate")->valueint;
		mLCSetting->gnssOutputDataRate = cJSON_GetObjectItem(root, "gnssOutputDataRate")->valueint;
		mLCSetting->insOutputDataRate = cJSON_GetObjectItem(root, "insOutputDataRate")->valueint;
		ret = 1;
	    cJSON_Delete(root);
	} while (0);

	free(data);
	data = NULL;
	return ret;
		 
}
