// csrc/pyext_link_shims.c
#ifndef FKPT_PYEXT
#define FKPT_PYEXT 1
#endif

#define main fkpt_main
#include "../src/main.c"          // pulls in globaldefs.h etc.

#include "../libs/stdinc.h"
#include <string.h>
#include <stdlib.h>

/* ---- minimal headline parser ---- */
static const char* _find_in_headlines(const char *key) {
    static char buf[256];
    const char *hs[4] = { gd.headline0, gd.headline1, gd.headline2, gd.headline3 };
    size_t klen = key ? strlen(key) : 0;
    if (!klen) return NULL;

    for (int i = 0; i < 4; ++i) {
        const char *h = hs[i];
        if (!h || !*h) continue;
        const char *p = h;
        while ((p = strstr(p, key))) {
            if ((p == h || p[-1] == ' ') && p[klen] == '=') {
                p += klen + 1;
                size_t j = 0;
                while (p[j] && p[j] != ' ' && j < sizeof(buf)-1) buf[j++] = p[j];
                buf[j] = '\0';
                return buf;
            }
            ++p;
        }
    }
    return NULL;
}

void   InitParam(string *h0, string *h1){(void)h0;(void)h1;}
string GetParam(string name){
    const char *v = _find_in_headlines(name);
    size_t n = v ? strlen(v) : 0;
    char *out = (char*)malloc(n + 1);
    if (!out) { out = (char*)malloc(1); if (out) out[0]='\0'; return out; }
    if (n) memcpy(out, v, n);
    out[n] = '\0';
    return out;
}
int    GetParamStat(string name){ return _find_in_headlines(name) ? ARGPARAM : 0; }
int    GetiParam(string name){ const char *v = _find_in_headlines(name); return v ? atoi(v) : 0; }
bool   GetbParam(string name){
    const char *v = _find_in_headlines(name);
    if (!v || !*v) return FALSE;
    if (strchr("tTyY1", *v)) return TRUE;
    if (strchr("fFnN0", *v)) return FALSE;
    return FALSE;
}
double GetdParam(string name){ const char *v = _find_in_headlines(name); return v ? atof(v) : 0.0; }

/* ---- keep IO stubs if libs/inout.c is excluded ---- */
void inout_InputData    (string f,int c1,int c2,int *n){(void)f;(void)c1;(void)c2; if(n)*n=0;}
void inout_InputData_1c (string f,int c2,int *n){(void)f;(void)c2; if(n)*n=0;}
void inout_InputData_3c (string f,int c1,int c2,int c3,int *n){(void)f;(void)c1;(void)c2;(void)c3; if(n)*n=0;}
void inout_InputData_4c (string f,int c1,int c2,int c3,int c4,int *n){(void)f;(void)c1;(void)c2;(void)c3;(void)c4; if(n)*n=0;}