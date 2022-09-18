/*==============================================================================
	HEADER: getparam.h			[General_libs]
==============================================================================*/


#ifndef _getparam_h
#define _getparam_h

void InitParam(string *, string *);

string GetParam(string);

int GetiParam(string);
bool GetbParam(string);
double GetdParam(string);

#define getargv0()    (GetParam("argv0"))
#define getversion()  (GetParam("Version"))

int GetParamStat(string);

#define DEFPARAM        001                     
#define REQPARAM        002                     
#define ARGPARAM        004                     
#define INPARAM         010                     
#define OUTPARAM        020                     
#define HIDPARAM        040                     

#endif  // ! _getparam_h

