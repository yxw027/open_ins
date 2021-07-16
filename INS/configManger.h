#ifndef CONFIG_MANAGER_H
#define CONFIG_MANAGER_H
#include <stdio.h>
#include <stdint.h>
#include "loosecoupleset.h"


int8_t readconfigfromfile(const char* cfgfname, LCSetting*  mLCSetting);


#endif