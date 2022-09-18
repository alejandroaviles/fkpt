/*==============================================================================
	MODULE: getparam.c			[General_libs]
Note: modified version of the zeno version
==============================================================================*/

#include "stdinc.h"
#include "getparam.h"
#include <string.h>

typedef struct {
    string name;
	string name_alias;
    string value;
    string comment;
    int flags;
} param;

local int CountDefaults(string *);
local void SetProgram(param *, string);
local void CopyDefaults(param *, string *);
local void CheckHelp(param *, string);
local void PrintHeader(string, string);
local void PrintItem(string, string, string);
local void SetArguments(param *, string *);
local void ReqArguments(param *);
local param *FindParam(string, param *);
local string ParName(string);

local string ParValue(string);                  

local param *paramvec = NULL;                   

local string progname = NULL;                   

void InitParam(string *argv, string *defv)
{
    int nparam;
    param *pvec;

    progname = argv[0];                         
    nparam = 1 + CountDefaults(defv);
    pvec = (param *) allocate(sizeof(param) * (nparam + 1));
    SetProgram(pvec, argv[0]);                  
    CopyDefaults(pvec, defv);                   
    CheckHelp(pvec, argv[1]);
    SetArguments(pvec, argv);
    ReqArguments(pvec);
    paramvec = pvec;
}

local int CountDefaults(string *defv)
{
    int ndefault;
    string *dp;

    ndefault = 0;
    for (dp = defv; *dp != NULL; dp++)          
        if (**dp != ';')                        
            ndefault++;                         
    return (ndefault);
}

local void SetProgram(param *pvec, string argv0)
{
    pvec->name = "argv0";
    pvec->name_alias = NULL;
    pvec->value = argv0;
    pvec->comment = NULL;
    pvec->flags = ARGPARAM;
}

local void CopyDefaults(param *pvec, string *defv)
{
    param *pp;
    string *dp, name, value;

    pp = pvec;                                  
    for (dp = defv; *dp != NULL; dp++)          
        if (**dp == ';') {                      
            if (pp->comment != NULL)            
                error("copydefaults: cant append comments\n");
            pp->comment = strdup(*dp + 1);      
        } else {
			if (**dp == ':') {
				pp->name_alias = strdup(*dp + 1);
			} else {
				pp++;
				name = ParName(*dp);
				value = ParValue(*dp);
				if (name == NULL || value == NULL)  
					error("copydefaults: bad parameter %s\n", *dp);
				pp->name = strdup(name);            
				pp->value = strdup(value);          
				pp->comment = NULL;                 
				pp->flags = DEFPARAM;               
				if (streq(pp->value, "???"))        
					pp->flags |= REQPARAM;          
				if (**dp == '<')                    
					pp->flags |= INPARAM;           
				else if (**dp == '>')               
					pp->flags |= OUTPARAM;          
				if (name[0] == '.')                 
					pp->flags |= HIDPARAM;
			}
        }
    pp++;                                       
    pp->name = NULL;                            
}

local void CheckHelp(param *pvec, string argv1)
{
    param *pp;
    char buf[128];

    if (argv1 != NULL && streq(argv1, "-clue")) {
                                               
        printf("%s", pvec->value);
        for (pp = pvec+1; pp->name != NULL; pp++)
            printf(" %s=%s", pp->name, pp->value);
        printf("\n");
        exit(0);
    }
    if (argv1 != NULL && streq(argv1, "-help")) {
        PrintHeader(pvec->value, pvec->comment);
        for (pp = pvec+1; pp->name != NULL; pp++) {
            sprintf(buf, "  %s=%s", pp->name, pp->value);
            PrintItem(buf, pp->comment, pp->name_alias);
        }
        exit(0);
    }
}

local void PrintHeader(string item, string comment)
{
    char buf[128];

    if (comment == NULL)
        printf("%s\n", item);
    else
        if (strlen(item) < 32) {
            sprintf(buf,"\n%s\n%s\n",item,comment);
            printf("%s\n", buf);
        } else
            printf("%s\n\t  %s\n", item, comment);
}

local void PrintItem(string item, string comment, string alias)
{
    if (comment == NULL)
        printf("%s\n", item);
    else
        if (strlen(item) < 32)
			if (alias == NULL)
				printf("%-32s  %s\n", item, comment);
			else
				printf("%-32s  %s [a: %s]\n", item, comment, alias);
        else
			if (alias == NULL)
				printf("%s\n\t\t\t\t  %s\n", item, comment);
			else
				printf("%s\n\t\t\t\t  %s [a: %s]\n", item, comment, alias);
}

local void SetArguments(param *pvec, string *argv)
{
    bool scanpos;
    param *pp;
    string *ap, name;

    scanpos = TRUE;
    pp = pvec;
    for (ap = argv + 1; *ap != NULL; ap++) { 
        name = ParName(*ap);                    
        scanpos = scanpos && (name == NULL);    
        if (scanpos) {                          
            pp++;                               
            if (pp->name == NULL)               
                error("%s: too many arguments\n", progname);
            pp->value = strdup(*ap);            
        } else {                                
            if (name == NULL)                   
                error("%s: nameless arg %s\n", progname, *ap);
            pp = FindParam(name, pvec);         
            if (pp == NULL)                     
                error("%s: parameter %s unknown\n", progname, name);
            if (pp->flags & ARGPARAM)           
                error("%s: parameter or alias %s duplicated\n", progname, name);
            pp->value = strdup(ParValue(*ap));  
        }
        pp->flags = (pp->flags & ~DEFPARAM) | ARGPARAM;
                                                
    }
}

local void ReqArguments(param *pvec)
{
    bool needarg;
    param *pp;

    needarg = FALSE;                            
    for (pp = pvec+1; pp->name != NULL; pp++)   
        if ((pp->flags & REQPARAM) &&           
              (pp->flags & DEFPARAM))           
            needarg = TRUE;                     
    if (needarg) {                              
        eprintf("Usage: %s", progname);
        for (pp = pvec+1; pp->name != NULL; pp++)
            if ((pp->flags & REQPARAM))
                eprintf(" %s=???", pp->name);
        error("%s: required arguments missing\n", progname);
    }
}

string GetParam(string name)
{
    param *par;

    if (paramvec == NULL) {
        if (streq(name, "argv0") && progname != NULL)
            return (progname);
        else
            error("getparam: called before initparam\n");
	}
    par = FindParam(name, paramvec);
    if (par == NULL)
        error("getparam in %s: parameter %s unknown\n", progname, name);
    return (par->value);
}

int GetParamStat(string name)
{
    param *par;

    par = FindParam(name, paramvec);
    return (par != NULL ? par->flags : 0);
}

int GetiParam(string name)
{
    return (atoi(GetParam(name)));              
}

double GetdParam(string name)
{
    return (atof(GetParam(name)));              
}

bool GetbParam(string name)
{
    char *val;

    val = GetParam(name);                       
    if (strchr("tTyY1", *val) != NULL)          
        return (TRUE);
    if (strchr("fFnN0", *val) != NULL)          
        return (FALSE);
    error("getbparam: %s=%s not bool\n", name, val);
    return (FALSE);                             
}

local param *FindParam(string name, param *pvec)
{
    param *pp;

    for (pp = pvec; pp->name != NULL; pp++) {
		if (pp->name_alias == NULL) {
			if (streq(name, pp->name))
				return (pp);
		}
		else
			if (streq(name, pp->name)) {
				return (pp);
			}
			else
				if (streq(name, pp->name_alias)) {
					return (pp);
				}
	}
    return (NULL);
}


// Removing the warning:
// warning: cast from pointer to integer of different size [-Wpointer-to-int-cast]
local string ParName(string arg)
{
    char *ap, *ep;
    static char namebuf[64];
    
    ap = (char *) arg;
    if (*ap == '<' || *ap == '>')
        ap++;
    strncpy(namebuf, ap, 63);
    namebuf[63] = '\0';
    ep = strchr(namebuf, '=');
    if (ep == NULL)
        return (NULL);
    *ep = '\0';
    return (namebuf);
}

local string ParValue(string arg)
{
    char *ep;

    ep = strchr(arg, '=');
    if (ep == NULL)
        return (NULL);
    return (ep + 1);
}
