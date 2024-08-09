#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#define STACK_BULGE1  1     /* stacking energies for bulges of size 1 */
#define NEW_NINIO     1     /* new asymetry penalty */
#define MAXSECTORS    500   /* dimension for a backtrack array */
#define LOCALITY      0.    /* locality parameter for base-pairs */
#define UNIT 100
#define MINPSCORE -2 * UNIT
#define NONE -10000         /* score for forbidden pairs */
#define PUBLIC
#define PRIVATE static
#define MAXALPHA 20       /* maximal length of alphabet */

static short  alias[MAXALPHA + 1];
static int    pair[MAXALPHA + 1][MAXALPHA + 1];
/* rtype[pair[i][j]]:=pair[j][i] */
static int    rtype[8] = {
  0, 2, 1, 4, 3, 6, 5, 7
};
#define NBASES 8
static const char Law_and_Order[]         = "_ACGUTXKI";
static int        BP_pair[NBASES][NBASES] =
  /* _  A  C  G  U  X  K  I */
{ { 0, 0, 0, 0, 0, 0, 0, 0 },
  { 0, 0, 0, 0, 5, 0, 0, 5 },
  { 0, 0, 0, 1, 0, 0, 0, 0 },
  { 0, 0, 2, 0, 3, 0, 0, 0 },
  { 0, 6, 0, 4, 0, 0, 0, 6 },
  { 0, 0, 0, 0, 0, 0, 2, 0 },
  { 0, 0, 0, 0, 0, 1, 0, 0 },
  { 0, 6, 0, 0, 5, 0, 0, 0 } };

/** The gas constant */
#define GASCONST 1.98717  /* in [cal/K] */
/** 0 deg Celsius in Kelvin */
#define K0  273.15
/** Infinity as used in minimization routines */
#define INF 10000000 /* (INT_MAX/10) */

#define EMAX (INF/10)
/** forbidden */
#define FORBIDDEN 9999
/** bonus contribution */
#define BONUS 10000
/** The number of distinguishable base pairs */
#define NBPAIRS 7
/** The minimum loop length */
#define TURN 3
/** The maximum loop length */
#define MAXLOOP 30
PUBLIC double Tmeasure = 37+K0;

PUBLIC double lxc37=107.856;
PUBLIC int ML_intern37=-90;
PUBLIC int ML_interndH=-220;
PUBLIC int ML_closing37=930;
PUBLIC int ML_closingdH=3000;
PUBLIC int ML_BASE37=0;
PUBLIC int ML_BASEdH=0;
PUBLIC int MAX_NINIO=300;
PUBLIC int ninio37=60;
PUBLIC int niniodH=320;
PUBLIC int TerminalAU37=50;
PUBLIC int TerminalAUdH=370;
PUBLIC int DuplexInit37=410;
PUBLIC int DuplexInitdH=360;
PUBLIC int TripleC37=100;
PUBLIC int TripleCdH=1860;
PUBLIC int MultipleCA37=30;
PUBLIC int MultipleCAdH=340;
PUBLIC int MultipleCB37=160;
PUBLIC int MultipleCBdH=760;

PUBLIC int GQuadAlpha37 = -1800;
PUBLIC int GQuadAlphadH = -11934;
PUBLIC int GQuadBeta37 = 1200;
PUBLIC int GQuadBetadH = 0;
PUBLIC int GQuadLayerMismatch37   = 300;
PUBLIC int GQuadLayerMismatchH    = 0;
PUBLIC int GQuadLayerMismatchMax  = 1;

#define saltT md->temperature+K0

#define VRNA_OPTION_WINDOW (1 << 4)
#define VRNA_OPTION_EVAL_ONLY (1 << 3)

#define WITH_PTYPE 1L
#define WITH_PTYPE_COMPAT 2L
#define VRNA_OPTION_EVAL_ONLY (1 << 3)

PUBLIC int stack37[NBPAIRS+1][NBPAIRS+1] =
{{   INF,   INF,   INF,   INF,   INF,   INF,   INF,   INF}
,{   INF,  -240,  -330,  -210,  -140,  -210,  -210,  -140}
,{   INF,  -330,  -340,  -250,  -150,  -220,  -240,  -150}
,{   INF,  -210,  -250,   130,   -50,  -140,  -130,   130}
,{   INF,  -140,  -150,   -50,    30,   -60,  -100,    30}
,{   INF,  -210,  -220,  -140,   -60,  -110,   -90,   -60}
,{   INF,  -210,  -240,  -130,  -100,   -90,  -130,   -90}
,{   INF,  -140,  -150,   130,    30,   -60,   -90,   130}};
PUBLIC int stackdH[NBPAIRS+1][NBPAIRS+1] =
{{   INF,   INF,   INF,   INF,   INF,   INF,   INF,   INF}
,{   INF, -1060, -1340, -1210,  -560, -1050, -1040,  -560}
,{   INF, -1340, -1490, -1260,  -830, -1140, -1240,  -830}
,{   INF, -1210, -1260, -1460, -1350,  -880, -1280,  -880}
,{   INF,  -560,  -830, -1350,  -930,  -320,  -700,  -320}
,{   INF, -1050, -1140,  -880,  -320,  -940,  -680,  -320}
,{   INF, -1040, -1240, -1280,  -700,  -680,  -770,  -680}
,{   INF,  -560,  -830,  -880,  -320,  -320,  -680,  -320}};

PUBLIC int hairpin37[31] = {   INF,   INF,   INF,   540,   560,   570,   540,   600,   550,   640,   650,   660,   670,   680,   690,   690,   700,   710,   710,   720,   720,   730,   730,   740,   740,   750,   750,   750,   760,   760,   770};
PUBLIC int hairpindH[31] = {   INF,   INF,   INF,   130,   480,   360,  -290,   130,  -290,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500};
PUBLIC int bulge37[31] = {   INF,   380,   280,   320,   360,   400,   440,   460,   470,   480,   490,   500,   510,   520,   530,   540,   540,   550,   550,   560,   570,   570,   580,   580,   580,   590,   590,   600,   600,   600,   610};
PUBLIC int bulgedH[31] = {   INF,  1060,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710};
PUBLIC int internal_loop37[31] = {   INF,   INF,   100,   100,   110,   200,   200,   210,   230,   240,   250,   260,   270,   280,   290,   290,   300,   310,   310,   320,   330,   330,   340,   340,   350,   350,   350,   360,   360,   370,   370};
PUBLIC int internal_loopdH[31] = {   INF,   INF,   -720,   -720,  -720,  -680,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130};

PUBLIC int mismatchI37[NBPAIRS+1][5][5] =
{{{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 }
,{{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,   -80,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,  -100,     0,  -100,     0}
 ,{     0,     0,     0,     0,   -60}
 }
,{{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,   -80,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,  -100,     0,  -100,     0}
 ,{     0,     0,     0,     0,   -60}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,   -10,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,   -30,    70,   -30,    70}
 ,{    70,    70,    70,    70,    10}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,   -10,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,   -30,    70,   -30,    70}
 ,{    70,    70,    70,    70,    10}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,   -10,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,   -30,    70,   -30,    70}
 ,{    70,    70,    70,    70,    10}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,   -10,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,   -30,    70,   -30,    70}
 ,{    70,    70,    70,    70,    10}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,   -10,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,   -30,    70,   -30,    70}
 ,{    70,    70,    70,    70,    10}
 }};
PUBLIC int mismatchIdH[NBPAIRS+1][5][5] =
{{{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 }
,{{   280,     0,     0,   280,     0}
 ,{     0,     0,     0,  -340,     0}
 ,{     0,     0,     0,     0,     0}
 ,{   280,  -760,     0,   280,     0}
 ,{     0,     0,     0,     0,  -580}
 }
,{{   280,     0,     0,   280,     0}
 ,{     0,     0,     0,  -340,     0}
 ,{     0,     0,     0,     0,     0}
 ,{   280,  -760,     0,   280,     0}
 ,{     0,     0,     0,     0,  -580}
 }
,{{   790,   500,   500,   790,   500}
 ,{   500,   500,   500,   170,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   790,  -260,   500,   790,   500}
 ,{   500,   500,   500,   500,   -80}
 }
,{{   790,   500,   500,   790,   500}
 ,{   500,   500,   500,   170,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   790,  -260,   500,   790,   500}
 ,{   500,   500,   500,   500,   -80}
 }
,{{   790,   500,   500,   790,   500}
 ,{   500,   500,   500,   170,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   790,  -260,   500,   790,   500}
 ,{   500,   500,   500,   500,   -80}
 }
,{{   790,   500,   500,   790,   500}
 ,{   500,   500,   500,   170,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   790,  -260,   500,   790,   500}
 ,{   500,   500,   500,   500,   -80}
 }
,{{   790,   500,   500,   790,   500}
 ,{   500,   500,   500,   170,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   790,  -260,   500,   790,   500}
 ,{   500,   500,   500,   500,   -80}
 }};

PUBLIC int mismatchH37[NBPAIRS+1][5][5] =
{{{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 }
,{{   -80,  -100,  -110,  -100,   -80}
 ,{  -140,  -150,  -150,  -140,  -150}
 ,{   -80,  -100,  -110,  -100,   -80}
 ,{  -150,  -230,  -150,  -240,  -150}
 ,{  -100,  -100,  -140,  -100,  -210}
 }
,{{   -50,  -110,   -70,  -110,   -50}
 ,{  -110,  -110,  -150,  -130,  -150}
 ,{   -50,  -110,   -70,  -110,   -50}
 ,{  -150,  -250,  -150,  -220,  -150}
 ,{  -100,  -110,  -100,  -110,  -160}
 }
,{{    20,    20,   -20,   -10,   -20}
 ,{    20,    20,   -50,   -30,   -50}
 ,{   -10,   -10,   -20,   -10,   -20}
 ,{   -50,  -100,   -50,  -110,   -50}
 ,{   -10,   -10,   -30,   -10,  -100}
 }
,{{     0,   -20,   -10,   -20,     0}
 ,{   -30,   -50,   -30,   -60,   -30}
 ,{     0,   -20,   -10,   -20,     0}
 ,{   -30,   -90,   -30,  -110,   -30}
 ,{   -10,   -20,   -10,   -20,   -90}
 }
,{{   -10,   -10,   -20,   -10,   -20}
 ,{   -30,   -30,   -50,   -30,   -50}
 ,{   -10,   -10,   -20,   -10,   -20}
 ,{   -50,  -120,   -50,  -110,   -50}
 ,{   -10,   -10,   -30,   -10,  -120}
 }
,{{     0,   -20,   -10,   -20,     0}
 ,{   -30,   -50,   -30,   -50,   -30}
 ,{     0,   -20,   -10,   -20,     0}
 ,{   -30,  -150,   -30,  -150,   -30}
 ,{   -10,   -20,   -10,   -20,   -90}
 }
,{{    20,    20,   -10,   -10,     0}
 ,{    20,    20,   -30,   -30,   -30}
 ,{     0,   -10,   -10,   -10,     0}
 ,{   -30,   -90,   -30,  -110,   -30}
 ,{   -10,   -10,   -10,   -10,   -90}
 }};
PUBLIC int mismatchHdH[NBPAIRS+1][5][5] =
{{{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 }
,{{   560,  -570,   560,  -560,  -270}
 ,{  -560,  -910,  -560,  -560,  -560}
 ,{  -270,  -570,  -340,  -570,  -270}
 ,{   560, -1400,   560,  -920,  -560}
 ,{  -530,  -570,  -530,  -570, -1440}
 }
,{{    50,  -520,    50,  -560,  -400}
 ,{  -400,  -520,  -400,  -560,  -400}
 ,{    50,  -720,    50,  -720,  -420}
 ,{  -400, -1290,  -400,  -620,  -400}
 ,{   -30,  -720,   -30,  -720, -1080}
 }
,{{   970,   140,   970,   140,   570}
 ,{   570,    30,   570,    20,   570}
 ,{   970,   140,   970,   140,   340}
 ,{   570,  -270,   570,    20,   570}
 ,{   830,   140,   830,   140,   -50}
 }
,{{   230,   100,   230,   220,   190}
 ,{  -110,  -110,  -260,  -520,  -260}
 ,{   190,   -60,  -140,   -60,   190}
 ,{   220,   100,  -260,   220,  -260}
 ,{   230,   -60,   230,   -60,   -70}
 }
,{{   970,   140,   970,   140,   570}
 ,{   570,   -20,   570,    20,   570}
 ,{   970,   140,   970,   140,   340}
 ,{   570,  -520,   570,    20,   570}
 ,{   830,   140,   830,   140,  -380}
 }
,{{   230,   -30,   230,   -60,   190}
 ,{   -30,   -30,  -260,  -520,  -260}
 ,{   190,   -60,  -140,   -60,   190}
 ,{  -260,  -590,  -260,  -520,  -260}
 ,{   230,   -60,   230,   -60,   -70}
 }
,{{   970,   140,   970,   220,   570}
 ,{   570,    30,   570,    20,   570}
 ,{   970,   140,   970,   140,   340}
 ,{   570,   100,   570,   220,   570}
 ,{   830,   140,   830,   140,   -50}
 }};

/* mismatch_multi */
PUBLIC int mismatchM37[NBPAIRS+1][5][5] =
{{ /* NP.. */
  {   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 },
 { /* CG.. */
  {   -50,  -110,   -50,  -140,   -70}
 ,{  -110,  -110,  -110,  -160,  -110}
 ,{   -70,  -150,   -70,  -150,  -100}
 ,{  -110,  -130,  -110,  -140,  -110}
 ,{   -50,  -150,   -50,  -150,   -70}
 },
 { /* GC.. */
  {   -80,  -140,   -80,  -140,  -100}
 ,{  -100,  -150,  -100,  -140,  -100}
 ,{  -110,  -150,  -110,  -150,  -140}
 ,{  -100,  -140,  -100,  -160,  -100}
 ,{   -80,  -150,   -80,  -150,  -120}
 },
 { /* GU.. */
  {   -50,   -80,   -50,   -50,   -50}
 ,{   -50,  -100,   -70,   -50,   -70}
 ,{   -60,   -80,   -60,   -80,   -60}
 ,{   -70,  -110,   -70,   -80,   -70}
 ,{   -50,   -80,   -50,   -80,   -50}
 },
 { /* UG.. */
  {   -30,   -30,   -60,   -60,   -60}
 ,{   -30,   -30,   -60,   -60,   -60}
 ,{   -70,  -100,   -70,  -100,   -80}
 ,{   -60,   -80,   -60,   -80,   -60}
 ,{   -60,  -100,   -70,  -100,   -60}
 },
 { /* AU.. */
  {   -50,   -80,   -50,   -80,   -50}
 ,{   -70,  -100,   -70,  -110,   -70}
 ,{   -60,   -80,   -60,   -80,   -60}
 ,{   -70,  -110,   -70,  -120,   -70}
 ,{   -50,   -80,   -50,   -80,   -50}
 },
 { /* UA.. */
  {   -60,   -80,   -60,   -80,   -60}
 ,{   -60,   -80,   -60,   -80,   -60}
 ,{   -70,  -100,   -70,  -100,   -80}
 ,{   -60,   -80,   -60,   -80,   -60}
 ,{   -70,  -100,   -70,  -100,   -80}
 },
 { /* NN.. */
  {   -30,   -30,   -50,   -50,   -50}
 ,{   -30,   -30,   -60,   -50,   -60}
 ,{   -60,   -80,   -60,   -80,   -60}
 ,{   -60,   -80,   -60,   -80,   -60}
 ,{   -50,   -80,   -50,   -80,   -50}
 }};

/* mismatch_multi_enthalpies */
PUBLIC int mismatchMdH[NBPAIRS+1][5][5] =
{{ /* NP.. */
  {   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 },
 { /* CG.. */
  {    50,  -400,    50,  -400,   -30}
 ,{  -520,  -520,  -720,  -710,  -720}
 ,{    50,  -400,    50,  -400,   -30}
 ,{  -560,  -560,  -720,  -620,  -720}
 ,{  -400,  -400,  -420,  -400,  -500}
 },
 { /* GC.. */
  {  -270,  -560,  -270,  -560,  -530}
 ,{  -570,  -910,  -570,  -820,  -570}
 ,{  -340,  -560,  -340,  -560,  -530}
 ,{  -560,  -560,  -570,  -920,  -570}
 ,{  -270,  -560,  -270,  -560,  -860}
 },
 { /* GU.. */
  {   310,  -480,  -180,   310,   140}
 ,{   310,  -480,  -430,   310,  -430}
 ,{  -140,  -630,  -510,  -630,  -140}
 ,{  -150,  -890,  -430,  -150,  -430}
 ,{   140,  -630,  -180,  -630,   140}
 },
 { /* UG.. */
  {   600,   200,   600,   200,   460}
 ,{   -60,  -340,  -230,   -60,  -230}
 ,{   600,   200,   600,   200,   460}
 ,{  -230,  -350,  -230,  -350,  -230}
 ,{   200,   200,   -30,   200,   160}
 },
 { /* AU.. */
  {   140,  -400,  -180,  -380,   140}
 ,{  -380,  -400,  -430,  -380,  -430}
 ,{  -140,  -630,  -510,  -630,  -140}
 ,{  -430,  -890,  -430,  -890,  -430}
 ,{   140,  -630,  -180,  -630,   140}
 },
 { /* UA.. */
  {   600,   200,   600,   200,   460}
 ,{  -230,  -390,  -230,  -310,  -230}
 ,{   600,   200,   600,   200,   460}
 ,{  -230,  -350,  -230,  -350,  -230}
 ,{   200,   200,   -30,   200,  -170}
 },
 { /* NN.. */
  {   600,   200,   600,   310,   460}
 ,{   310,  -340,  -230,   310,  -230}
 ,{   600,   200,   600,   200,   460}
 ,{  -150,  -350,  -230,  -150,  -230}
 ,{   200,   200,   -30,   200,   160}
 }};

PUBLIC int mismatch1nI37[NBPAIRS+1][5][5] =
{{{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 }
,{{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 }
,{{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 }};
PUBLIC int mismatch1nIdH[NBPAIRS+1][5][5] =
{{{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 }
,{{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 }
,{{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 }
,{{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 }
,{{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 }
,{{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 }
,{{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 }
,{{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 }};

PUBLIC int mismatch23I37[NBPAIRS+1][5][5] =
{{{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 }
,{{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,   -50,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,  -110,     0,   -70,     0}
 ,{     0,     0,     0,     0,   -30}
 }
,{{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,  -120,     0,   -70,     0}
 ,{     0,     0,     0,     0,   -30}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,   -40,    70,     0,    70}
 ,{    70,    70,    70,    70,    40}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    20,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,   -40,    70,     0,    70}
 ,{    70,    70,    70,    70,    40}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,   -40,    70,     0,    70}
 ,{    70,    70,    70,    70,    40}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    20,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,   -40,    70,     0,    70}
 ,{    70,    70,    70,    70,    40}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,   -40,    70,     0,    70}
 ,{    70,    70,    70,    70,    40}
 }};
PUBLIC int mismatch23IdH[NBPAIRS+1][5][5] =
{{{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 }
,{{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,  -570,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,  -860,     0,  -900,     0}
 ,{     0,     0,     0,     0,  -640}
 }
,{{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0, -1090,     0,  -900,     0}
 ,{     0,     0,     0,     0,  -640}
 }
,{{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,  -580,   500,  -400,   500}
 ,{   500,   500,   500,   500,  -140}
 }
,{{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   -60,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,  -360,   500,  -400,   500}
 ,{   500,   500,   500,   500,  -140}
 }
,{{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,  -580,   500,  -400,   500}
 ,{   500,   500,   500,   500,  -140}
 }
,{{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   -60,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,  -360,   500,  -400,   500}
 ,{   500,   500,   500,   500,  -140}
 }
,{{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,  -360,   500,  -400,   500}
 ,{   500,   500,   500,   500,  -140}
 }};

/* mismatch_exterior */
PUBLIC int mismatchExt37[NBPAIRS+1][5][5] =
{{ /* NP.. */
  {   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 },
 { /* CG.. */
  {   -50,  -110,   -50,  -140,   -70}
 ,{  -110,  -110,  -110,  -160,  -110}
 ,{   -70,  -150,   -70,  -150,  -100}
 ,{  -110,  -130,  -110,  -140,  -110}
 ,{   -50,  -150,   -50,  -150,   -70}
 },
 { /* GC.. */
  {   -80,  -140,   -80,  -140,  -100}
 ,{  -100,  -150,  -100,  -140,  -100}
 ,{  -110,  -150,  -110,  -150,  -140}
 ,{  -100,  -140,  -100,  -160,  -100}
 ,{   -80,  -150,   -80,  -150,  -120}
 },
 { /* GU.. */
  {   -50,   -80,   -50,   -50,   -50}
 ,{   -50,  -100,   -70,   -50,   -70}
 ,{   -60,   -80,   -60,   -80,   -60}
 ,{   -70,  -110,   -70,   -80,   -70}
 ,{   -50,   -80,   -50,   -80,   -50}
 },
 { /* UG.. */
  {   -30,   -30,   -60,   -60,   -60}
 ,{   -30,   -30,   -60,   -60,   -60}
 ,{   -70,  -100,   -70,  -100,   -80}
 ,{   -60,   -80,   -60,   -80,   -60}
 ,{   -60,  -100,   -70,  -100,   -60}
 },
 { /* AU.. */
  {   -50,   -80,   -50,   -80,   -50}
 ,{   -70,  -100,   -70,  -110,   -70}
 ,{   -60,   -80,   -60,   -80,   -60}
 ,{   -70,  -110,   -70,  -120,   -70}
 ,{   -50,   -80,   -50,   -80,   -50}
 },
 { /* UA.. */
  {   -60,   -80,   -60,   -80,   -60}
 ,{   -60,   -80,   -60,   -80,   -60}
 ,{   -70,  -100,   -70,  -100,   -80}
 ,{   -60,   -80,   -60,   -80,   -60}
 ,{   -70,  -100,   -70,  -100,   -80}
 },
 { /* NN.. */
  {   -30,   -30,   -50,   -50,   -50}
 ,{   -30,   -30,   -60,   -50,   -60}
 ,{   -60,   -80,   -60,   -80,   -60}
 ,{   -60,   -80,   -60,   -80,   -60}
 ,{   -50,   -80,   -50,   -80,   -50}
 }};

/* mismatch_exterior_enthalpies */
PUBLIC int mismatchExtdH[NBPAIRS+1][5][5] =
{{ /* NP.. */
  {   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 },
 { /* CG.. */
  {    50,  -400,    50,  -400,   -30}
 ,{  -520,  -520,  -720,  -710,  -720}
 ,{    50,  -400,    50,  -400,   -30}
 ,{  -560,  -560,  -720,  -620,  -720}
 ,{  -400,  -400,  -420,  -400,  -500}
 },
 { /* GC.. */
  {  -270,  -560,  -270,  -560,  -530}
 ,{  -570,  -910,  -570,  -820,  -570}
 ,{  -340,  -560,  -340,  -560,  -530}
 ,{  -560,  -560,  -570,  -920,  -570}
 ,{  -270,  -560,  -270,  -560,  -860}
 },
 { /* GU.. */
  {   310,  -480,  -180,   310,   140}
 ,{   310,  -480,  -430,   310,  -430}
 ,{  -140,  -630,  -510,  -630,  -140}
 ,{  -150,  -890,  -430,  -150,  -430}
 ,{   140,  -630,  -180,  -630,   140}
 },
 { /* UG.. */
  {   600,   200,   600,   200,   460}
 ,{   -60,  -340,  -230,   -60,  -230}
 ,{   600,   200,   600,   200,   460}
 ,{  -230,  -350,  -230,  -350,  -230}
 ,{   200,   200,   -30,   200,   160}
 },
 { /* AU.. */
  {   140,  -400,  -180,  -380,   140}
 ,{  -380,  -400,  -430,  -380,  -430}
 ,{  -140,  -630,  -510,  -630,  -140}
 ,{  -430,  -890,  -430,  -890,  -430}
 ,{   140,  -630,  -180,  -630,   140}
 },
 { /* UA.. */
  {   600,   200,   600,   200,   460}
 ,{  -230,  -390,  -230,  -310,  -230}
 ,{   600,   200,   600,   200,   460}
 ,{  -230,  -350,  -230,  -350,  -230}
 ,{   200,   200,   -30,   200,  -170}
 },
 { /* NN.. */
  {   600,   200,   600,   310,   460}
 ,{   310,  -340,  -230,   310,  -230}
 ,{   600,   200,   600,   200,   460}
 ,{  -150,  -350,  -230,  -150,  -230}
 ,{   200,   200,   -30,   200,   160}
 }};

/* dangle5 */
PUBLIC int dangle5_37[NBPAIRS+1][5] =
{ /*           N      A      C      G      U */
/* NP */ {   INF,   INF,   INF,   INF,   INF},
/* CG */ {   -10,   -50,   -30,   -20,   -10},
/* GC */ {    -0,   -20,   -30,    -0,    -0},
/* GU */ {   -20,   -30,   -30,   -40,   -20},
/* UG */ {   -10,   -30,   -10,   -20,   -20},
/* AU */ {   -20,   -30,   -30,   -40,   -20},
/* UA */ {   -10,   -30,   -10,   -20,   -20},
/* NN */ {    -0,   -20,   -10,    -0,    -0}
};

/* dangle3 */
PUBLIC int dangle3_37[NBPAIRS+1][5] =
{ /*           N      A      C      G      U */
/* NP */ {   INF,   INF,   INF,   INF,   INF},
/* CG */ {   -40,  -110,   -40,  -130,   -60},
/* GC */ {   -80,  -170,   -80,  -170,  -120},
/* GU */ {   -10,   -70,   -10,   -70,   -10},
/* UG */ {   -50,   -80,   -50,   -80,   -60},
/* AU */ {   -10,   -70,   -10,   -70,   -10},
/* UA */ {   -50,   -80,   -50,   -80,   -60},
/* NN */ {   -10,   -70,   -10,   -70,   -10}
};

/* dangle5_enthalpies */
PUBLIC int dangle5_dH[NBPAIRS+1][5] =
{ /*           N      A      C      G      U */
/* NP */ {   INF,   INF,   INF,   INF,   INF},
/* CG */ {   330,  -240,   330,    80,  -140},
/* GC */ {    70,  -160,    70,  -460,   -40},
/* GU */ {   310,   160,   220,    70,   310},
/* UG */ {   690,   -50,   690,    60,    60},
/* AU */ {   310,   160,   220,    70,   310},
/* UA */ {   690,   -50,   690,    60,    60},
/* NN */ {   690,   160,   690,    80,   310}
};

/* dangle3_enthalpies */
PUBLIC int dangle3_dH[NBPAIRS+1][5] =
{ /*           N      A      C      G      U */
/* NP */ {   INF,   INF,   INF,   INF,   INF},
/* CG */ {  -280,  -740,  -280,  -640,  -360},
/* GC */ {  -410,  -900,  -410,  -860,  -750},
/* GU */ {   -70,  -570,   -70,  -580,  -220},
/* UG */ {   -90,  -490,   -90,  -550,  -230},
/* AU */ {   -70,  -570,   -70,  -580,  -220},
/* UA */ {   -90,  -490,   -90,  -550,  -230},
/* NN */ {   -70,  -490,   -70,  -550,  -220}
};

PUBLIC char Triloops[241] =
  "CAACG "
  "GUUAC "
;
PUBLIC int Triloop37[40] = {   680,   690};
PUBLIC int TriloopdH[40] = {  2370,  1080};

PUBLIC char Tetraloops[281] =
  "CAACGG "
  "CCAAGG "
  "CCACGG "
  "CCCAGG "
  "CCGAGG "
  "CCGCGG "
  "CCUAGG "
  "CCUCGG "
  "CUAAGG "
  "CUACGG "
  "CUCAGG "
  "CUCCGG "
  "CUGCGG "
  "CUUAGG "
  "CUUCGG "
  "CUUUGG "
;
PUBLIC int Tetraloop37[40] = {   550,   330,   370,   340,   350,   360,   370,   250,   360,   280,   370,   270,   280,   350,   370,   370};
PUBLIC int TetraloopdH[40] = {   690, -1030,  -330,  -890,  -660,  -750,  -350, -1390,  -760, -1070,  -660, -1290, -1070,  -620, -1530,  -680};

PUBLIC char Hexaloops[361] =
  "ACAGUACU "
  "ACAGUGAU "
  "ACAGUGCU "
  "ACAGUGUU "
;
PUBLIC int Hexaloop37[40] = {   280,   360,   290,   180};
PUBLIC int HexaloopdH[40] = { -1680, -1140, -1280, -1540};

#define   VRNA_GQUAD_MAX_STACK_SIZE     7
#define   VRNA_GQUAD_MIN_STACK_SIZE     2
#define   VRNA_GQUAD_MAX_LINKER_LENGTH  15
#define   VRNA_GQUAD_MIN_LINKER_LENGTH  1
#define   VRNA_GQUAD_MIN_BOX_SIZE       ((4 * VRNA_GQUAD_MIN_STACK_SIZE) + \
                                         (3 * VRNA_GQUAD_MIN_LINKER_LENGTH))
#define   VRNA_GQUAD_MAX_BOX_SIZE       ((4 * VRNA_GQUAD_MAX_STACK_SIZE) + \
                                         (3 * VRNA_GQUAD_MAX_LINKER_LENGTH))

#define VRNA_MODEL_DEFAULT_TEMPERATURE    37.0

#define VRNA_MODEL_DEFAULT_TEMPERATURE    37.0

/**
 *  @brief  Default scaling factor for partition function computations
 *
 *  @see  #vrna_exp_param_t.pf_scale, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_PF_SCALE       -1

/**
 *  @brief  Default scaling factor for absolute thermodynamic temperature in Boltzmann factors
 *
 *  @see    #vrna_exp_param_t.alpha, #vrna_md_t.betaScale, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_BETA_SCALE     1.

/**
 *  @brief  Default dangling end model
 *
 *  @see  #vrna_md_t.dangles, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_DANGLES        2

/**
 *  @brief  Default model behavior for lookup of special tri-, tetra-, and hexa-loops
 *
 *  @see    #vrna_md_t.special_hp, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_SPECIAL_HP     1

/**
 *  @brief  Default model behavior for so-called 'lonely pairs'
 *
 *  @see    #vrna_md_t.noLP, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_NO_LP          0

/**
 *  @brief  Default model behavior for G-U base pairs
 *
 *  @see    #vrna_md_t.noGU, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_NO_GU          0

/**
 *  @brief  Default model behavior for G-U base pairs closing a loop
 *
 *  @see    #vrna_md_t.noGUclosure, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_NO_GU_CLOSURE  0

/**
 *  @brief  Default model behavior to treat a molecule as a circular RNA (DNA)
 *  @see    #vrna_md_t.circ, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_CIRC           0

/**
 *  @brief  Default model behavior regarding the treatment of G-Quadruplexes
 *
 *  @see    #vrna_md_t.gquad, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_GQUAD          0

/**
 *  @brief  Default behavior of the model regarding unique multi-branch loop decomposition
 *
 *  @see    #vrna_md_t.uniq_ML, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_UNIQ_ML        0

/**
 *  @brief  Default model behavior on which energy set to use
 *
 *  @see    #vrna_md_t.energy_set, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_ENERGY_SET     0

/**
 *  @brief  Default model behavior with regards to backtracking of structures
 *  @see    #vrna_md_t.backtrack, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_BACKTRACK      1

/**
 *  @brief  Default model behavior on what type of backtracking to perform
 *
 *  @see    #vrna_md_t.backtrack_type, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_BACKTRACK_TYPE 'F'

/**
 *  @brief  Default model behavior with regards to computing base pair probabilities
 *
 *  @see    #vrna_md_t.compute_bpp, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_COMPUTE_BPP    1

/**
 *  @brief  Default model behavior for the allowed maximum base pair span
 *
 *  @see    #vrna_md_t.max_bp_span, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_MAX_BP_SPAN    -1

/**
 *  @brief  Default model behavior for the sliding window approach
 *
 *  @see    #vrna_md_t.window_size, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_WINDOW_SIZE    -1

/**
 *  @brief  Default model behavior on how to evaluate the energy contribution of multi-branch loops
 *
 *  @see    #vrna_md_t.logML, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_LOG_ML         0

/**
 *  @brief  Default model behavior for consensus structure energy evaluation
 *  @see    #vrna_md_t.oldAliEn, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_ALI_OLD_EN     0

/**
 *  @brief  Default model behavior for consensus structure co-variance contribution assessment
 *
 *  @see    #vrna_md_t.ribo, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_ALI_RIBO       0

/**
 *  @brief  Default model behavior for weighting the co-variance score in consensus structure prediction
 *
 *  @see    #vrna_md_t.cv_fact, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_ALI_CV_FACT    1.

/** @brief  Default model behavior for weighting the nucleotide conservation? in consensus structure prediction
 *
 *  @see    #vrna_md_t.nc_fact, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_ALI_NC_FACT    1.


#define VRNA_MODEL_DEFAULT_PF_SMOOTH      1

/**
 *  @brief  Default model salt concentration (M)
 */
#define VRNA_MODEL_DEFAULT_SALT           1.021


/**
 *  @brief  Default model lower bound of multiloop size for salt correction fiting
 */
#define VRNA_MODEL_DEFAULT_SALT_MLLOWER    6


/**
 *  @brief  Default model upper bound of multiloop size for salt correction fiting
 */
#define VRNA_MODEL_DEFAULT_SALT_MLUPPER    24


/**
 *  @brief  Default model value to turn off user-provided salt correction for duplex initializtion
 */
#define VRNA_MODEL_DEFAULT_SALT_DPXINIT       99999

#define VRNA_MODEL_SALT_DPXINIT_FACT_RNA      -45.324
#define VRNA_MODEL_SALT_DPXINIT_FACT_DNA      -58.389


#define VRNA_MODEL_DEFAULT_SALT_DPXINIT_FACT  VRNA_MODEL_SALT_DPXINIT_FACT_RNA

/* Geometric parameters for RNA and DNA */

#define VRNA_MODEL_HELICAL_RISE_RNA   2.8
#define VRNA_MODEL_HELICAL_RISE_DNA   3.4
/**
 *  @brief  Default helical rise
 */
#define VRNA_MODEL_DEFAULT_HELICAL_RISE   VRNA_MODEL_HELICAL_RISE_RNA

#define VRNA_MODEL_BACKBONE_LENGTH_RNA   6.0
#define VRNA_MODEL_BACKBONE_LENGTH_DNA   6.76
/**
 *  @brief  Default backbone length
 */
#define VRNA_MODEL_DEFAULT_BACKBONE_LENGTH   VRNA_MODEL_BACKBONE_LENGTH_RNA


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

#ifndef MAXALPHA
/**
 *  @brief Maximal length of alphabet
 */
#define MAXALPHA              20

#endif

#endif

#define   BP_REV_DEFAULT        { 0, 2, 1, 4, 3, 6, 5, 7 }

#define   BP_ALIAS_DEFAULT      { 0, 1, 2, 3, 4, 3, 2, 0 }

#define   BP_ENCODING_DEFAULT \
  /*  _  A  C  G  U  X  K  I */ \
  { { 0, 0, 0, 0, 0, 0, 0, 0 }, \
    { 0, 0, 0, 0, 5, 0, 0, 5 }, \
    { 0, 0, 0, 1, 0, 0, 0, 0 }, \
    { 0, 0, 2, 0, 3, 0, 0, 0 }, \
    { 0, 6, 0, 4, 0, 0, 0, 6 }, \
    { 0, 0, 0, 0, 0, 0, 2, 0 }, \
    { 0, 0, 0, 0, 0, 1, 0, 0 }, \
    { 0, 6, 0, 0, 5, 0, 0, 0 } }

#define   DM_DEFAULT \
  { { 0, 0, 0, 0, 0, 0, 0 }, /* hamming distance between pairs */ \
    { 0, 0, 2, 2, 1, 2, 2 } /* CG */, \
    { 0, 2, 0, 1, 2, 2, 2 } /* GC */, \
    { 0, 2, 1, 0, 2, 1, 2 } /* GU */, \
    { 0, 1, 2, 2, 0, 2, 1 } /* UG */, \
    { 0, 2, 2, 1, 2, 0, 2 } /* AU */, \
    { 0, 2, 2, 2, 1, 2, 0 } /* UA */ }

typedef double FLT_OR_DBL;


double          temperature     = VRNA_MODEL_DEFAULT_TEMPERATURE;
double          pf_scale        = VRNA_MODEL_DEFAULT_PF_SCALE;
int             dangles         = VRNA_MODEL_DEFAULT_DANGLES;
int             tetra_loop      = VRNA_MODEL_DEFAULT_SPECIAL_HP;
int             noLonelyPairs   = VRNA_MODEL_DEFAULT_NO_LP;
int             noGU            = 1;
int             no_closingGU    = VRNA_MODEL_DEFAULT_NO_GU_CLOSURE;
int             circ            = VRNA_MODEL_DEFAULT_CIRC;
int             gquad           = VRNA_MODEL_DEFAULT_GQUAD;
int             uniq_ML         = VRNA_MODEL_DEFAULT_UNIQ_ML;
int             energy_set      = VRNA_MODEL_DEFAULT_ENERGY_SET;
int             do_backtrack    = VRNA_MODEL_DEFAULT_COMPUTE_BPP;
char            backtrack_type  = VRNA_MODEL_DEFAULT_BACKTRACK_TYPE;
char            *nonstandards   = NULL;
int             max_bp_span     = VRNA_MODEL_DEFAULT_MAX_BP_SPAN;
int             oldAliEn        = VRNA_MODEL_DEFAULT_ALI_OLD_EN;
int             ribo            = VRNA_MODEL_DEFAULT_ALI_RIBO;
double          cv_fact         = VRNA_MODEL_DEFAULT_ALI_CV_FACT;
double          nc_fact         = VRNA_MODEL_DEFAULT_ALI_NC_FACT;
int             logML           = VRNA_MODEL_DEFAULT_LOG_ML;

/* below are some more deprecated global symbols we need to get rid off */

int             james_rule        = 1;    /* interior loops of size 2 get energy 0.8Kcal and
                                           * no mismatches (no longer used) */
char            *RibosumFile      = NULL; /* TODO: compile ribosums into program
                                           * Warning: this variable will vanish */
int             csv               = 0;    /*generate comma seperated output*/
// vrna_bp_stack_t *base_pair        = NULL;
FLT_OR_DBL      *pr               = NULL; /* base pairing prob. matrix */
int             *iindx            = NULL; /* pr[i,j] -> pr[iindx[i]-j] */
int             fold_constrained  = 0;    /* fold with constraints */

double          salt             = VRNA_MODEL_DEFAULT_SALT;
int             saltDPXInit      = VRNA_MODEL_DEFAULT_SALT_DPXINIT;
float           saltDPXInitFact  = VRNA_MODEL_DEFAULT_SALT_DPXINIT_FACT;
float           helical_rise     = VRNA_MODEL_DEFAULT_HELICAL_RISE;
float           backbone_length  = VRNA_MODEL_DEFAULT_BACKBONE_LENGTH;

/*
 #################################
 # GLOBAL VARIABLES              #
 #################################
 */

/*
 #################################
 # PRIVATE VARIABLES             #
 #################################
 */

struct vrna_md_s {
  double  temperature;                      /**<  @brief  The temperature used to scale the thermodynamic parameters */
  double  betaScale;                        /**<  @brief  A scaling factor for the thermodynamic temperature of the Boltzmann factors */
  int     pf_smooth;                        /**<  @brief  A flat specifying whether energies in Boltzmann factors need to be smoothed */
  int     dangles;                          /**<  @brief  Specifies the dangle model used in any energy evaluation (0,1,2 or 3)
                                             *
                                             *    If set to 0 no stabilizing energies are assigned to bases adjacent to
                                             *    helices in free ends and multiloops (so called dangling ends). Normally
                                             *    (dangles = 1) dangling end energies are assigned only to unpaired
                                             *    bases and a base cannot participate simultaneously in two dangling ends. In
                                             *    the partition function algorithm vrna_pf() these checks are neglected.
                                             *    To provide comparability between free energy minimization and partition function
                                             *    algorithms, the default setting is 2.
                                             *    This treatment of dangling ends gives more favorable energies to helices
                                             *    directly adjacent to one another, which can be beneficial since such
                                             *    helices often do engage in stabilizing interactions through co-axial
                                             *    stacking.\n
                                             *    If set to 3 co-axial stacking is explicitly included for
                                             *    adjacent helices in multiloops. The option affects only mfe folding
                                             *    and energy evaluation (vrna_mfe() and vrna_eval_structure()), as
                                             *    well as suboptimal folding (vrna_subopt()) via re-evaluation of energies.
                                             *    Co-axial stacking with one intervening mismatch is not considered so far.
                                             *    Note, that some function do not implement all dangle model but only a subset of
                                             *    (0,1,2,3). In particular, partition function algorithms can only handle
                                             *    0 and 2. Read the documentation of the particular recurrences or
                                             *    energy evaluation function for information about the provided dangle
                                             *    model.
                                             */
  int     special_hp;                       /**<  @brief  Include special hairpin contributions for tri, tetra and hexaloops */
  int     noLP;                             /**<  @brief  Only consider canonical structures, i.e. no 'lonely' base pairs */
  int     noGU;                             /**<  @brief  Do not allow GU pairs */
  int     noGUclosure;                      /**<  @brief  Do not allow loops to be closed by GU pair */
  int     logML;                            /**<  @brief  Use logarithmic scaling for multiloops */
  int     circ;                             /**<  @brief  Assume RNA to be circular instead of linear */
  int     gquad;                            /**<  @brief  Include G-quadruplexes in structure prediction */
  int     uniq_ML;                          /**<  @brief  Flag to ensure unique multi-branch loop decomposition during folding */
  int     energy_set;                       /**<  @brief  Specifies the energy set that defines set of compatible base pairs */
  int     backtrack;                        /**<  @brief  Specifies whether or not secondary structures should be backtraced */
  char    backtrack_type;                   /**<  @brief  Specifies in which matrix to backtrack */
  int     compute_bpp;                      /**<  @brief  Specifies whether or not backward recursions for base pair probability (bpp) computation will be performed */
  char    nonstandards[64];                 /**<  @brief  contains allowed non standard bases */
  int     max_bp_span;                      /**<  @brief  maximum allowed base pair span */

  int     min_loop_size;                    /**<  @brief  Minimum size of hairpin loops
                                             *
                                             *    The default value for this field is #TURN, however, it may
                                             *    be 0 in cofolding context.
                                             */
  int     window_size;                      /**<  @brief  Size of the sliding window for locally optimal structure prediction */
  int     oldAliEn;                         /**<  @brief  Use old alifold energy model */
  int     ribo;                             /**<  @brief  Use ribosum scoring table in alifold energy model */
  double  cv_fact;                          /**<  @brief  Co-variance scaling factor for consensus structure prediction */
  double  nc_fact;                          /**<  @brief  Scaling factor to weight co-variance contributions of non-canonical pairs */
  double  sfact;                            /**<  @brief  Scaling factor for partition function scaling */
  int     rtype[8];                         /**<  @brief  Reverse base pair type array */
  short   alias[MAXALPHA + 1];              /**<  @brief  alias of an integer nucleotide representation */
  int     pair[MAXALPHA + 1][MAXALPHA + 1]; /**<  @brief  Integer representation of a base pair */
  float   pair_dist[7][7];                  /**<  @brief  Base pair dissimilarity, a.k.a. distance matrix */
  double  salt;                             /**<  @brief  Salt (monovalent) concentration (M) in buffer */
  int     saltMLLower;                      /**<  @brief  Lower bound of multiloop size to use in loop salt correction linear fitting */
  int     saltMLUpper;                      /**<  @brief  Upper bound of multiloop size to use in loop salt correction linear fitting */
  int     saltDPXInit;                      /**<  @brief  User-provided salt correction for duplex initialization (in dcal/mol).
                                             *    If set to 99999 the default salt correction is used.
                                             *    If set to 0 there is no salt correction for duplex initialization.
                                             */
  float   saltDPXInitFact;                  /**<  @brief  */
  float   helical_rise;                     /**<  @brief  */
  float   backbone_length;                  /**<  @brief  */
};



typedef struct vrna_md_s vrna_md_t;


struct vrna_param_s {
  int       id;
  int       stack[NBPAIRS + 1][NBPAIRS + 1];
  int       hairpin[31];
  int       bulge[MAXLOOP + 1];
  int       internal_loop[MAXLOOP + 1];
  int       mismatchExt[NBPAIRS + 1][5][5];
  int       mismatchI[NBPAIRS + 1][5][5];
  int       mismatch1nI[NBPAIRS + 1][5][5];
  int       mismatch23I[NBPAIRS + 1][5][5];
  int       mismatchH[NBPAIRS + 1][5][5];
  int       mismatchM[NBPAIRS + 1][5][5];
  int       dangle5[NBPAIRS + 1][5];
  int       dangle3[NBPAIRS + 1][5];
  int       int11[NBPAIRS + 1][NBPAIRS + 1][5][5];
  int       int21[NBPAIRS + 1][NBPAIRS + 1][5][5][5];
  int       int22[NBPAIRS + 1][NBPAIRS + 1][5][5][5][5];
  int       ninio[5];
  double    lxc;
  int       MLbase;
  int       MLintern[NBPAIRS + 1];
  int       MLclosing;
  int       TerminalAU;
  int       DuplexInit;
  int       Tetraloop_E[200];
  char      Tetraloops[1401];
  int       Triloop_E[40];
  char      Triloops[241];
  int       Hexaloop_E[40];
  char      Hexaloops[1801];
  int       TripleC;
  int       MultipleCA;
  int       MultipleCB;
  int       gquad[VRNA_GQUAD_MAX_STACK_SIZE + 1][3 * VRNA_GQUAD_MAX_LINKER_LENGTH + 1];
  int       gquadLayerMismatch;
  int       gquadLayerMismatchMax;

  double    temperature;      /**<  @brief  Temperature used for loop contribution scaling */

  vrna_md_t model_details;    /**<  @brief  Model details to be used in the recursions */
  char      param_file[256];  /**<  @brief  The filename the parameters were derived from, or empty string if they represent the default */
  int       SaltStack;
  int       SaltLoop[MAXLOOP + 2];
  double    SaltLoopDbl[MAXLOOP + 2];
  int       SaltMLbase;
  int       SaltMLintern;
  int       SaltMLclosing;
  int       SaltDPXInit;
};


typedef struct  vrna_param_s vrna_param_t;
PRIVATE vrna_md_t defaults = {
  VRNA_MODEL_DEFAULT_TEMPERATURE,
  1.,
  VRNA_MODEL_DEFAULT_PF_SMOOTH,
  VRNA_MODEL_DEFAULT_DANGLES,
  VRNA_MODEL_DEFAULT_SPECIAL_HP,
  VRNA_MODEL_DEFAULT_NO_LP,
  VRNA_MODEL_DEFAULT_NO_GU,
  VRNA_MODEL_DEFAULT_NO_GU_CLOSURE,
  VRNA_MODEL_DEFAULT_LOG_ML,
  VRNA_MODEL_DEFAULT_CIRC,
  VRNA_MODEL_DEFAULT_GQUAD,
  VRNA_MODEL_DEFAULT_UNIQ_ML,
  VRNA_MODEL_DEFAULT_ENERGY_SET,
  VRNA_MODEL_DEFAULT_BACKTRACK,
  VRNA_MODEL_DEFAULT_BACKTRACK_TYPE,
  VRNA_MODEL_DEFAULT_COMPUTE_BPP,
  { 0 },
  VRNA_MODEL_DEFAULT_MAX_BP_SPAN,
  TURN,
  VRNA_MODEL_DEFAULT_WINDOW_SIZE,
  VRNA_MODEL_DEFAULT_ALI_OLD_EN,
  VRNA_MODEL_DEFAULT_ALI_RIBO,
  VRNA_MODEL_DEFAULT_ALI_CV_FACT,
  VRNA_MODEL_DEFAULT_ALI_NC_FACT,
  1.07,
  BP_REV_DEFAULT,
  BP_ALIAS_DEFAULT,
  BP_ENCODING_DEFAULT,
  DM_DEFAULT,
  VRNA_MODEL_DEFAULT_SALT,
  VRNA_MODEL_DEFAULT_SALT_MLLOWER,
  VRNA_MODEL_DEFAULT_SALT_MLUPPER,
  VRNA_MODEL_DEFAULT_SALT_DPXINIT,
  VRNA_MODEL_DEFAULT_SALT_DPXINIT_FACT,
  VRNA_MODEL_DEFAULT_HELICAL_RISE,
  VRNA_MODEL_DEFAULT_BACKBONE_LENGTH
};

PRIVATE vrna_param_t  *P = NULL;
PRIVATE int           **c = NULL;     /* energy array, given that i-j pair */
PRIVATE short         *S1 = NULL, *SS1 = NULL, *S2 = NULL, *SS2 = NULL;
PRIVATE int           n1, n2;         /* sequence lengths */

typedef enum {
  VRNA_FC_TYPE_SINGLE,      /**< Type is suitable for single, and hybridizing sequences */
  VRNA_FC_TYPE_COMPARATIVE  /**< Type is suitable for sequence alignments (consensus structure prediction) */
} vrna_fc_type_e;

struct vrna_elem_prob_s {
  int   i;    /**<  @brief  Start position (usually 5' nucleotide that starts the element, e.g. base pair) */
  int   j;    /**<  @brief  End position (usually 3' nucleotide that ends the element, e.g. base pair) */
  float p;    /**<  @brief  Probability of the element */
  int   type; /**<  @brief  Type of the element */
};

typedef struct vrna_elem_prob_s vrna_ep_t;


struct vrna_cstr_s {
  char          *string;
  size_t        size;
  FILE          *output;
  unsigned char istty;
};

typedef struct vrna_cstr_s *vrna_cstr_t;
struct output_stream {
  vrna_cstr_t data;
  vrna_cstr_t err;
};


#define VRNA_OPTION_DEFAULT 0U
#define VRNA_OPTION_HYBRID (1 << 2)

struct id_data {
  char      *name;
  int       auto_id;
  char      *prefix;
  char      *delimiter;
  int       digits;
  long int  number;
};
typedef struct id_data *dataset_id;


typedef void (*vrna_stream_output_f)(void        *auxdata,
                                            unsigned int i,
                                            void         *data);

struct vrna_ordered_stream_s {
  unsigned int                start;      /* first element index in queue, i.e. start of queue */
  unsigned int                end;        /* last element index in queue */
  unsigned int                size;       /* available memory size for 'data' and 'provided' */
  unsigned int                shift;      /* pointer offset for 'data' and 'provided' */

  vrna_stream_output_f output;    /* callback to execute if consecutive elements from head are available */
  void                        **data;     /* actual data passed to the callback */
  unsigned char               *provided;  /* for simplicity we use unsigned char instead of single bits per element */
  void                        *auxdata;   /* auxiliary data passed to the callback */
#if VRNA_WITH_PTHREADS
  pthread_mutex_t             mtx;        /* semaphore to provide concurrent access */
#endif
};

typedef struct vrna_ordered_stream_s *vrna_ostream_t;

#define MAX_ALPHABET (6)
#define MAX_PAIRS (NBPAIRS + 1 + 25)

struct vrna_sc_mod_param_s {
  unsigned int  available;

  char          *name;
  char          one_letter_code;
  char          unmodified;
  char          fallback;
  char          pairing_partners[7];
  unsigned int  pairing_partners_encoding[7];
  unsigned int  unmodified_encoding;
  unsigned int  fallback_encoding;

  size_t        num_ptypes;
  size_t        ptypes[MAX_ALPHABET][MAX_ALPHABET];

  int           stack_dG[MAX_PAIRS][MAX_ALPHABET][MAX_ALPHABET];
  int           stack_dH[MAX_PAIRS][MAX_ALPHABET][MAX_ALPHABET];

  int           dangle5_dG[MAX_PAIRS][MAX_ALPHABET];
  int           dangle5_dH[MAX_PAIRS][MAX_ALPHABET];
  int           dangle3_dG[MAX_PAIRS][MAX_ALPHABET];
  int           dangle3_dH[MAX_PAIRS][MAX_ALPHABET];

  int           mismatch_dG[MAX_PAIRS][MAX_ALPHABET][MAX_ALPHABET];
  int           mismatch_dH[MAX_PAIRS][MAX_ALPHABET][MAX_ALPHABET];

  int           terminal_dG[MAX_PAIRS];
  int           terminal_dH[MAX_PAIRS];
};

typedef struct vrna_sc_mod_param_s *vrna_sc_mod_param_t;
struct options {
  int             filename_full;
  char            *filename_delim;
  int             pf;
  int             doT;
  int             doC;
  int             noPS;
  int             noconv;
  int             centroid;
  int             MEA;
  double          MEAgamma;
  double          bppmThreshold;
  int             verbose;
  vrna_md_t       md;
  // vrna_cmd_t      commands;

  dataset_id      id_control;

  char            *concentration_file;

  char            *constraint_file;
  int             constraint_batch;
  int             constraint_enforce;
  int             constraint_canonical;

  int             shape;
  char            *shape_file;
  char            *shape_method;
  char            *shape_conversion;

  vrna_sc_mod_param_t *mod_params;

  int             csv_output;
  int             csv_header;
  char            csv_output_delim;

  int             jobs;
  int             keep_order;
  unsigned int    next_record_number;
  vrna_ostream_t  output_queue;
};


struct record_data {
  unsigned int    number;
  char            *id;
  char            *sequence;
  char            *SEQ_ID;
  char            *input_filename;
  int             multiline_input;
  struct options  *options;
  int             tty;
};

#ifdef VRNA_WARN_DEPRECATED
# if defined(__clang__)
#  define DEPRECATED(func, msg) func __attribute__ ((deprecated("", msg)))
# elif defined(__GNUC__)
#  define DEPRECATED(func, msg) func __attribute__ ((deprecated(msg)))
# else
#  define DEPRECATED(func, msg) func
# endif
#else
# define DEPRECATED(func, msg) func
#endif

typedef enum {
  VRNA_SEQ_UNKNOWN = 0,   /**< @brief Nucleotide sequence represents an Unkown type */
  VRNA_SEQ_RNA = 1,       /**< @brief Nucleotide sequence represents an RNA type */
  VRNA_SEQ_DNA = 2        /**< @brief Nucleotide sequence represents a DNA type */
} vrna_seq_type_e;

struct vrna_sequence_s {
  vrna_seq_type_e type;       /**< @brief The type of sequence */
  char            *name;
  char            *string;    /**< @brief The string representation of the sequence */
  short           *encoding;  /**< @brief The integer representation of the sequence */
  short           *encoding5;
  short           *encoding3;
  unsigned int    length;     /**< @brief The length of the sequence */
};


typedef struct vrna_sequence_s vrna_seq_t;

typedef enum {
  VRNA_HC_DEFAULT = 0,  /**<  @brief  Default Hard Constraints */
  VRNA_HC_WINDOW  = 1  /**<  @brief  Hard Constraints suitable for local structure prediction using
                     *    window approach.
                     *    @see    vrna_mfe_window(), vrna_mfe_window_zscore(), pfl_fold()
                     */
} vrna_hc_type_e;

typedef unsigned char (*vrna_hc_eval_f)(int           i,
                                        int           j,
                                        int           k,
                                        int           l,
                                        unsigned char d,
                                        void          *data);

typedef void (*vrna_auxdata_free_f)(void *data);

struct hc_nuc {
  int           direction;
  unsigned char context;
  unsigned char nonspec;
};

struct hc_basepair {
  size_t        list_size;
  size_t        list_mem;
  unsigned int  *j;
  unsigned int  *strand_j;
  unsigned char *context;
};

struct vrna_hc_depot_s {
  unsigned int            strands;
  size_t                  *up_size;
  struct hc_nuc           **up;
  size_t                  *bp_size;
  struct hc_basepair      **bp;
};


typedef struct vrna_hc_depot_s vrna_hc_depot_t;

struct vrna_hc_s {
  vrna_hc_type_e  type;
  unsigned int    n;

  unsigned char   state;

#ifndef VRNA_DISABLE_C11_FEATURES
  /* C11 support for unnamed unions/structs */
  union {
    struct {
#endif
  unsigned char *mx;
#ifndef VRNA_DISABLE_C11_FEATURES
};
struct {
#endif
  unsigned char **matrix_local;
#ifndef VRNA_DISABLE_C11_FEATURES
};
};
#endif

  int                 *up_ext;    /**<  @brief  A linear array that holds the number of allowed
                                   *            unpaired nucleotides in an exterior loop
                                   */
  int                 *up_hp;     /**<  @brief  A linear array that holds the number of allowed
                                   *            unpaired nucleotides in a hairpin loop
                                   */
  int                 *up_int;    /**<  @brief  A linear array that holds the number of allowed
                                   *            unpaired nucleotides in an interior loop
                                   */
  int                 *up_ml;     /**<  @brief  A linear array that holds the number of allowed
                                   *            unpaired nucleotides in a multi branched loop
                                   */

  vrna_hc_eval_f      f;          /**<  @brief  A function pointer that returns whether or
                                   *            not a certain decomposition may be evaluated
                                   */

  void                *data;      /**<  @brief  A pointer to some structure where the user
                                   *            may store necessary data to evaluate its
                                   *            generic hard constraint function
                                   */

  vrna_auxdata_free_f free_data;  /**<  @brief  A pointer to a function to free memory
                                   *            occupied by auxiliary data
                                   *
                                   *    The function this pointer is pointing to will be
                                   *    called upon destruction of the #vrna_hc_s, and
                                   *    provided with the vrna_hc_s.data pointer that
                                   *    may hold auxiliary data. Hence, to avoid leaking
                                   *    memory, the user may use this pointer to free
                                   *    memory occupied by auxiliary data.
                                   */

  vrna_hc_depot_t *depot;
};


typedef struct  vrna_hc_s vrna_hc_t;


typedef enum {
  VRNA_MX_DEFAULT = 0,  /**<  @brief  Default DP matrices */
  VRNA_MX_WINDOW = 1,   /**<  @brief  DP matrices suitable for local structure prediction using
                     *    window approach.
                     *    @see    vrna_mfe_window(), vrna_mfe_window_zscore(), pfl_fold()
                     */
  VRNA_MX_2DFOLD = 2    /**<  @brief  DP matrices suitable for distance class partitioned structure prediction
                     *    @see  vrna_mfe_TwoD(), vrna_pf_TwoD()
                     */
} vrna_mx_type_e;

struct vrna_mx_mfe_s {
  /** @name Common fields for MFE matrices
   *  @{
   */
  const vrna_mx_type_e  type;     /**< Type of the DP matrices */
  unsigned int          length;   /**<  @brief  Length of the sequence, therefore an indicator of the size of the DP matrices */
  unsigned int          strands;  /**< Number of strands */
  /**
   *  @}
   */

#ifndef VRNA_DISABLE_C11_FEATURES
  /* C11 support for unnamed unions/structs */
  union {
    struct {
#endif
  /** @name Default DP matrices
   *  @note These data fields are available if
   *        @code vrna_mx_mfe_t.type == VRNA_MX_DEFAULT @endcode
   * @{
   */
  int *c;           /**<  @brief  Energy array, given that i-j pair */
  int *f5;          /**<  @brief  Energy of 5' end */
  int *f3;          /**<  @brief  Energy of 3' end */
  int **fms5;       /**<  @brief  Energy for connected interstrand configurations */
  int **fms3;       /**<  @brief  nergy for connected interstrand configurations */
  int *fML;         /**<  @brief  Multi-loop auxiliary energy array */
  int *fM1;         /**<  @brief  Second ML array, only for unique multibrnach loop decomposition */
  int *fM2;         /**<  @brief  Energy for a multibranch loop region with exactly two stems, extending to 3' end */
  int *ggg;         /**<  @brief  Energies of g-quadruplexes */
  int Fc;           /**<  @brief  Minimum Free Energy of entire circular RNA */
  int FcH;          /**<  @brief  Minimum Free Energy of hairpin loop cases in circular RNA */
  int FcI;          /**<  @brief  Minimum Free Energy of internal loop cases in circular RNA */
  int FcM;          /**<  @brief  Minimum Free Energy of multibranch loop cases in circular RNA */
  /**
   * @}
   */

#ifndef VRNA_DISABLE_C11_FEATURES
  /* C11 support for unnamed unions/structs */
};
struct {
#endif
  /** @name Local Folding DP matrices using window approach
   *  @note These data fields are available if
   *        @code vrna_mx_mfe_t.type == VRNA_MX_WINDOW @endcode
   * @{
   */
  int **c_local;            /**<  @brief  Energy array, given that i-j pair */
  int *f3_local;            /**<  @brief  Energy of 5' end */
  int **fML_local;          /**<  @brief  Multi-loop auxiliary energy array */
  int **ggg_local;          /**<  @brief  Energies of g-quadruplexes */
  /**
   * @}
   */
#ifndef VRNA_DISABLE_C11_FEATURES
  /* C11 support for unnamed unions/structs */
};
struct {
#endif

  /** @name Distance Class DP matrices
   *  @note These data fields are available if
   *        @code vrna_mx_mfe_t.type == VRNA_MX_2DFOLD @endcode
   * @{
   */
  int           ***E_F5;
  int           **l_min_F5;
  int           **l_max_F5;
  int           *k_min_F5;
  int           *k_max_F5;

  int           ***E_F3;
  int           **l_min_F3;
  int           **l_max_F3;
  int           *k_min_F3;
  int           *k_max_F3;

  int           ***E_C;
  int           **l_min_C;
  int           **l_max_C;
  int           *k_min_C;
  int           *k_max_C;

  int           ***E_M;
  int           **l_min_M;
  int           **l_max_M;
  int           *k_min_M;
  int           *k_max_M;

  int           ***E_M1;
  int           **l_min_M1;
  int           **l_max_M1;
  int           *k_min_M1;
  int           *k_max_M1;

  int           ***E_M2;
  int           **l_min_M2;
  int           **l_max_M2;
  int           *k_min_M2;
  int           *k_max_M2;

  int           **E_Fc;
  int           *l_min_Fc;
  int           *l_max_Fc;
  int           k_min_Fc;
  int           k_max_Fc;

  int           **E_FcH;
  int           *l_min_FcH;
  int           *l_max_FcH;
  int           k_min_FcH;
  int           k_max_FcH;

  int           **E_FcI;
  int           *l_min_FcI;
  int           *l_max_FcI;
  int           k_min_FcI;
  int           k_max_FcI;

  int           **E_FcM;
  int           *l_min_FcM;
  int           *l_max_FcM;
  int           k_min_FcM;
  int           k_max_FcM;

  /* auxilary arrays for remaining set of coarse graining (k,l) > (k_max, l_max) */
  int           *E_F5_rem;
  int           *E_F3_rem;
  int           *E_C_rem;
  int           *E_M_rem;
  int           *E_M1_rem;
  int           *E_M2_rem;

  int           E_Fc_rem;
  int           E_FcH_rem;
  int           E_FcI_rem;
  int           E_FcM_rem;

#ifdef COUNT_STATES
  unsigned long ***N_F5;
  unsigned long ***N_C;
  unsigned long ***N_M;
  unsigned long ***N_M1;
#endif

  /**
   * @}
   */

#ifndef VRNA_DISABLE_C11_FEATURES
  /* C11 support for unnamed unions/structs */
};
};
#endif
};



typedef struct  vrna_mx_mfe_s vrna_mx_mfe_t;


struct vrna_mx_pf_s {
  /** @name Common fields for DP matrices
   *  @{
   */
  const vrna_mx_type_e  type;       /**< Type of the DP matrices */
  unsigned int          length;     /**< Size of the DP matrices (i.e. sequence length) */
  FLT_OR_DBL            *scale;     /**< Boltzmann factor scaling */
  FLT_OR_DBL            *expMLbase; /**< Boltzmann factors for unpaired bases in multibranch loop */

  /**
   *  @}
   */

#ifndef VRNA_DISABLE_C11_FEATURES
  /* C11 support for unnamed unions/structs */
  union {
    struct {
#endif

  /** @name Default PF matrices
   *  @note These data fields are available if
   *        @code vrna_mx_pf_t.type == VRNA_MX_DEFAULT @endcode
   *  @{
   */
  FLT_OR_DBL *q;
  FLT_OR_DBL *qb;
  FLT_OR_DBL *qm;
  FLT_OR_DBL *qm1;
  FLT_OR_DBL *probs;
  FLT_OR_DBL *q1k;
  FLT_OR_DBL *qln;
  FLT_OR_DBL *G;

  FLT_OR_DBL qo;
  FLT_OR_DBL *qm2;
  FLT_OR_DBL qho;
  FLT_OR_DBL qio;
  FLT_OR_DBL qmo;

  /**
   *  @}
   */

#ifndef VRNA_DISABLE_C11_FEATURES
  /* C11 support for unnamed unions/structs */
};
struct {
#endif

  /** @name Local Folding DP matrices using window approach
   *  @note These data fields are available if
   *        @code vrna_mx_mfe_t.type == VRNA_MX_WINDOW @endcode
   * @{
   */
  FLT_OR_DBL **q_local;
  FLT_OR_DBL **qb_local;
  FLT_OR_DBL **qm_local;
  FLT_OR_DBL **pR;
  FLT_OR_DBL **qm2_local;
  FLT_OR_DBL **QI5;
  FLT_OR_DBL **q2l;
  FLT_OR_DBL **qmb;
  FLT_OR_DBL **G_local;
  /**
   *  @}
   */

#ifndef VRNA_DISABLE_C11_FEATURES
  /* C11 support for unnamed unions/structs */
};
struct {
#endif

  /** @name Distance Class DP matrices
   *  @note These data fields are available if
   *        @code vrna_mx_pf_t.type == VRNA_MX_2DFOLD @endcode
   *  @{
   */
  FLT_OR_DBL ***Q;
  int **l_min_Q;
  int **l_max_Q;
  int *k_min_Q;
  int *k_max_Q;


  FLT_OR_DBL ***Q_B;
  int **l_min_Q_B;
  int **l_max_Q_B;
  int *k_min_Q_B;
  int *k_max_Q_B;

  FLT_OR_DBL ***Q_M;
  int **l_min_Q_M;
  int **l_max_Q_M;
  int *k_min_Q_M;
  int *k_max_Q_M;

  FLT_OR_DBL ***Q_M1;
  int **l_min_Q_M1;
  int **l_max_Q_M1;
  int *k_min_Q_M1;
  int *k_max_Q_M1;

  FLT_OR_DBL ***Q_M2;
  int **l_min_Q_M2;
  int **l_max_Q_M2;
  int *k_min_Q_M2;
  int *k_max_Q_M2;

  FLT_OR_DBL **Q_c;
  int *l_min_Q_c;
  int *l_max_Q_c;
  int k_min_Q_c;
  int k_max_Q_c;

  FLT_OR_DBL **Q_cH;
  int *l_min_Q_cH;
  int *l_max_Q_cH;
  int k_min_Q_cH;
  int k_max_Q_cH;

  FLT_OR_DBL **Q_cI;
  int *l_min_Q_cI;
  int *l_max_Q_cI;
  int k_min_Q_cI;
  int k_max_Q_cI;

  FLT_OR_DBL **Q_cM;
  int *l_min_Q_cM;
  int *l_max_Q_cM;
  int k_min_Q_cM;
  int k_max_Q_cM;

  /* auxilary arrays for remaining set of coarse graining (k,l) > (k_max, l_max) */
  FLT_OR_DBL *Q_rem;
  FLT_OR_DBL *Q_B_rem;
  FLT_OR_DBL *Q_M_rem;
  FLT_OR_DBL *Q_M1_rem;
  FLT_OR_DBL *Q_M2_rem;

  FLT_OR_DBL Q_c_rem;
  FLT_OR_DBL Q_cH_rem;
  FLT_OR_DBL Q_cI_rem;
  FLT_OR_DBL Q_cM_rem;
  /**
   *  @}
   */

#ifndef VRNA_DISABLE_C11_FEATURES
  /* C11 support for unnamed unions/structs */
};
};
#endif
};

typedef struct  vrna_mx_pf_s vrna_mx_pf_t;


typedef void (*vrna_recursion_status_f)(unsigned char status,
                                              void          *data);



struct vrna_structured_domains_s {
  char __placeholder; /* dummy placeholder to not leave this struct empty */
};
typedef struct vrna_structured_domains_s vrna_sd_t;





struct vrna_exp_param_s {
  int     id;   /**<  @brief  An identifier for the data structure
                 *    @deprecated This attribute will be removed in version 3
                 */
  double  expstack[NBPAIRS + 1][NBPAIRS + 1];
  double  exphairpin[31];
  double  expbulge[MAXLOOP + 1];
  double  expinternal[MAXLOOP + 1];
  double  expmismatchExt[NBPAIRS + 1][5][5];
  double  expmismatchI[NBPAIRS + 1][5][5];
  double  expmismatch23I[NBPAIRS + 1][5][5];
  double  expmismatch1nI[NBPAIRS + 1][5][5];
  double  expmismatchH[NBPAIRS + 1][5][5];
  double  expmismatchM[NBPAIRS + 1][5][5];
  double  expdangle5[NBPAIRS + 1][5];
  double  expdangle3[NBPAIRS + 1][5];
  double  expint11[NBPAIRS + 1][NBPAIRS + 1][5][5];
  double  expint21[NBPAIRS + 1][NBPAIRS + 1][5][5][5];
  double  expint22[NBPAIRS + 1][NBPAIRS + 1][5][5][5][5];
  double  expninio[5][MAXLOOP + 1];
  double  lxc;
  double  expMLbase;
  double  expMLintern[NBPAIRS + 1];
  double  expMLclosing;
  double  expTermAU;
  double  expDuplexInit;
  double  exptetra[40];
  double  exptri[40];
  double  exphex[40];
  char    Tetraloops[1401];
  double  expTriloop[40];
  char    Triloops[241];
  char    Hexaloops[1801];
  double  expTripleC;
  double  expMultipleCA;
  double  expMultipleCB;
  double  expgquad[VRNA_GQUAD_MAX_STACK_SIZE + 1][3 * VRNA_GQUAD_MAX_LINKER_LENGTH + 1];
  double  expgquadLayerMismatch;
  int     gquadLayerMismatchMax;

  double  kT;
  double  pf_scale;           /**<  @brief    Scaling factor to avoid over-/underflows */

  double  temperature;        /**<  @brief    Temperature used for loop contribution scaling */
  double  alpha;              /**<  @brief    Scaling factor for the thermodynamic temperature
                               *    @details  This allows for temperature scaling in Boltzmann
                               *              factors independently from the energy contributions.
                               *              The resulting Boltzmann factors are then computed by
                               *              @f$ e^{-E/(\alpha \cdot K \cdot T)} @f$
                               */

  vrna_md_t model_details;    /**<  @brief  Model details to be used in the recursions */
  char      param_file[256];  /**<  @brief  The filename the parameters were derived from, or empty string if they represent the default */

  double    expSaltStack;
  double    expSaltLoop[MAXLOOP + 2];
  double    SaltLoopDbl[MAXLOOP + 2];
  int       SaltMLbase;
  int       SaltMLintern;
  int       SaltMLclosing;
  int       SaltDPXInit;
};


typedef struct  vrna_exp_param_s vrna_exp_param_t;



typedef void (*vrna_ud_production_f)(vrna_fold_compound_t *fc,
                                           void                 *data);
typedef void (*vrna_ud_exp_production_f)(vrna_fold_compound_t *fc,
                                               void                 *data);

typedef int (*vrna_ud_f)(vrna_fold_compound_t  *fc,
                                      int                   i,
                                      int                   j,
                                      unsigned int          loop_type,
                                      void                  *data);

typedef FLT_OR_DBL (*vrna_ud_exp_f)(vrna_fold_compound_t *fc,
                                                 int                  i,
                                                 int                  j,
                                                 unsigned int         loop_type,
                                                 void                 *data);

typedef void (*vrna_ud_add_probs_f)(vrna_fold_compound_t  *fc,
                                          int                   i,
                                          int                   j,
                                          unsigned int          loop_type,
                                          FLT_OR_DBL            exp_energy,
                                          void                  *data);


typedef FLT_OR_DBL (*vrna_ud_get_probs_f)(vrna_fold_compound_t  *fc,
                                                int                   i,
                                                int                   j,
                                                unsigned int          loop_type,
                                                int                   motif,
                                                void                  *data);


struct vrna_unstructured_domain_s {
  /*
   **********************************
   * Keep track of all motifs added
   **********************************
   */
  int           uniq_motif_count;                   /**<  @brief The unique number of motifs of different lengths */
  unsigned int  *uniq_motif_size;                   /**<  @brief An array storing a unique list of motif lengths */

  int           motif_count;                        /**<  @brief Total number of distinguished motifs */
  char          **motif;                            /**<  @brief Motif sequences */
  char          **motif_name;                       /**<  @brief Motif identifier/name */
  unsigned int  *motif_size;                        /**<  @brief Motif lengths */
  double        *motif_en;                          /**<  @brief Ligand binding free energy contribution */
  unsigned int  *motif_type;                        /**<  @brief Type of motif, i.e. loop type the ligand binds to */

  /*
   **********************************
   * Grammar extension for ligand
   * binding
   **********************************
   */
  vrna_ud_production_f     prod_cb;        /**<  @brief Callback to ligand binding production rule, i.e. create/fill DP free energy matrices
                                                   *    @details This callback will be executed right before the actual secondary structure decompositions,
                                                   *    and, therefore, any implementation must not interleave with the regular DP matrices.
                                                   */
  vrna_ud_exp_production_f exp_prod_cb;    /**<  @brief Callback to ligand binding production rule, i.e. create/fill DP partition function matrices */
  vrna_ud_f         energy_cb;      /**<  @brief Callback to evaluate free energy of ligand binding to a particular unpaired stretch */
  vrna_ud_exp_f     exp_energy_cb;  /**<  @brief Callback to evaluate Boltzmann factor of ligand binding to a particular unpaired stretch */
  void                            *data;          /**<  @brief Auxiliary data structure passed to energy evaluation callbacks */
  vrna_auxdata_free_f      free_data;      /**<  @brief Callback to free auxiliary data structure */
  vrna_ud_add_probs_f      probs_add;      /**<  @brief Callback to store/add outside partition function */
  vrna_ud_get_probs_f      probs_get;      /**<  @brief Callback to retrieve outside partition function */
};


typedef struct vrna_unstructured_domain_s vrna_ud_t;


struct vrna_fc_s {
  /**
   *  @name Common data fields
   *  @{
   */
  const vrna_fc_type_e type;        /**<  @brief  The type of the #vrna_fold_compound_t.
                                     * @details Currently possible values are #VRNA_FC_TYPE_SINGLE, and #VRNA_FC_TYPE_COMPARATIVE
                                     * @warning Do not edit this attribute, it will be automagically set by
                                     *      the corresponding get() methods for the #vrna_fold_compound_t.
                                     *      The value specified in this attribute dictates the set of other
                                     *      attributes to use within this data structure.
                                     */
  unsigned int      length;         /**<  @brief  The length of the sequence (or sequence alignment) */
#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
  DEPRECATED(int cutpoint,
             "Use strand_* members instead");
                                    /**<  @brief  The position of the (cofold) cutpoint within the provided sequence.
                                     * If there is no cutpoint, this field will be set to -1
                                     */
#endif
  unsigned int      *strand_number; /**<  @brief  The strand number a particular nucleotide is associated with */
  unsigned int      *strand_order;  /**<  @brief  The strand order, i.e. permutation of current concatenated sequence */
  unsigned int      *strand_order_uniq; /**<  @brief  The strand order array where identical sequences have the same ID */
  unsigned int      *strand_start;  /**<  @brief  The start position of a particular strand within the current concatenated sequence */
  unsigned int      *strand_end;    /**<  @brief  The end (last) position of a particular strand within the current concatenated sequence */

  unsigned int      strands;        /**<  @brief  Number of interacting strands */
  vrna_seq_t        *nucleotides;   /**<  @brief  Set of nucleotide sequences */

  vrna_hc_t         *hc;            /**<  @brief  The hard constraints data structure used for structure prediction */

  vrna_mx_mfe_t     *matrices;      /**<  @brief  The MFE DP matrices */
  vrna_mx_pf_t      *exp_matrices;  /**<  @brief  The PF DP matrices  */

  vrna_param_t      *params;        /**<  @brief  The precomputed free energy contributions for each type of loop */
  vrna_exp_param_t  *exp_params;    /**<  @brief  The precomputed free energy contributions as Boltzmann factors  */

  int               *iindx;         /**<  @brief  DP matrix accessor  */
  int               *jindx;         /**<  @brief  DP matrix accessor  */

  /**
   *  @}
   *
   *  @name User-defined data fields
   *  @{
   */
  vrna_recursion_status_f   stat_cb;       /**<  @brief  Recursion status callback (usually called just before, and
                                            *            after recursive computations in the library
                                            *    @see    vrna_recursion_status_f(), vrna_fold_compound_add_callback()
                                            */

  void                      *auxdata;      /**<  @brief  A pointer to auxiliary, user-defined data
                                            *    @see vrna_fold_compound_add_auxdata(), #vrna_fold_compound_t.free_auxdata
                                            */

  vrna_auxdata_free_f       free_auxdata;  /**<  @brief A callback to free auxiliary user data whenever the fold_compound itself is free'd
                                            *    @see  #vrna_fold_compound_t.auxdata, vrna_auxdata_free_f()
                                            */

  /**
   *  @}
   *
   *  @name Secondary Structure Decomposition (grammar) related data fields
   *  @{
   */

  /* data structure to adjust additional structural domains, such as G-quadruplexes */
  vrna_sd_t     *domains_struc;             /**<  @brief  Additional structured domains */

  /* data structure to adjust additional contributions to unpaired stretches, e.g. due to protein binding */
  vrna_ud_t     *domains_up;                /**<  @brief  Additional unstructured domains */

  /* auxiliary (user-defined) extension to the folding grammar */
  vrna_gr_aux_t *aux_grammar;               /**<  @brief  Additional decomposition grammar rules */

  /**
   *  @}
   */

#ifndef VRNA_DISABLE_C11_FEATURES
  /* C11 support for unnamed unions/structs */
  union {
    struct {
#endif

  /**
   *  @name Data fields available for single/hybrid structure prediction
   *  @{
   */
      char *sequence;                   /**<  @brief  The input sequence string
                                         *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_SINGLE @endverbatim
                                         */
      short *sequence_encoding;         /**<  @brief  Numerical encoding of the input sequence
                                         *    @see    vrna_sequence_encode()
                                         *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_SINGLE @endverbatim
                                         */
      short *encoding5;
      short *encoding3;

      short *sequence_encoding2;


      char *ptype;                      /**<  @brief  Pair type array
                                         *
                                         *    Contains the numerical encoding of the pair type for each pair (i,j) used
                                         *    in MFE, Partition function and Evaluation computations.
                                         *    @note This array is always indexed via jindx, in contrast to previously
                                         *    different indexing between mfe and pf variants!
                                         *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_SINGLE @endverbatim
                                         *    @see    vrna_idx_col_wise(), vrna_ptypes()
                                         */
      char *ptype_pf_compat;            /**<  @brief  ptype array indexed via iindx
                                         *    @deprecated  This attribute will vanish in the future!
                                         *    It's meant for backward compatibility only!
                                         *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_SINGLE @endverbatim
                                         */
      vrna_sc_t *sc;                    /**<  @brief  The soft constraints for usage in structure prediction and evaluation
                                         *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_SINGLE @endverbatim
                                         */

  /**
   *  @}
   */

#ifndef VRNA_DISABLE_C11_FEATURES
  /* C11 support for unnamed unions/structs */
};
struct {
#endif

  /**
   *  @name Data fields for consensus structure prediction
   *  @{
   */
      char          **sequences;        /**<  @brief  The aligned sequences
                                         *    @note   The end of the alignment is indicated by a NULL pointer in the second dimension
                                         *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
                                         */
      unsigned int  n_seq;              /**<  @brief  The number of sequences in the alignment
                                         *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
                                         */
      char          *cons_seq;          /**<  @brief  The consensus sequence of the aligned sequences
                                         *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
                                         */
      short         *S_cons;            /**<  @brief  Numerical encoding of the consensus sequence
                                         *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
                                         */
      short         **S;                /**<  @brief  Numerical encoding of the sequences in the alignment
                                         *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
                                         */
      short         **S5;               /**<  @brief    S5[s][i] holds next base 5' of i in sequence s
                                         *    @warning  Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
                                         */
      short         **S3;               /**<  @brief    Sl[s][i] holds next base 3' of i in sequence s
                                         *    @warning  Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
                                         */
  char          **Ss;
  unsigned int  **a2s;
      int           *pscore;              /**<  @brief  Precomputed array of pair types expressed as pairing scores
                                           *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
                                           */
      int           **pscore_local;       /**<  @brief  Precomputed array of pair types expressed as pairing scores
                                           *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
                                           */
      short         *pscore_pf_compat;    /**<  @brief  Precomputed array of pair types expressed as pairing scores indexed via iindx
                                           *    @deprecated  This attribute will vanish in the future!
                                           *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
                                           */
      vrna_sc_t     **scs;                /**<  @brief  A set of soft constraints (for each sequence in the alignment)
                                           *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
                                           */
  int           oldAliEn;

  /**
   *  @}
   */
#ifndef VRNA_DISABLE_C11_FEATURES
};
};
#endif

  /**
   *  @name Additional data fields for Distance Class Partitioning
   *
   *  These data fields are typically populated with meaningful data only if used in the context of Distance Class Partitioning
   *  @{
   */
  unsigned int  maxD1;            /**<  @brief  Maximum allowed base pair distance to first reference */
  unsigned int  maxD2;            /**<  @brief  Maximum allowed base pair distance to second reference */
  short         *reference_pt1;   /**<  @brief  A pairtable of the first reference structure */
  short         *reference_pt2;   /**<  @brief  A pairtable of the second reference structure */

  unsigned int  *referenceBPs1;   /**<  @brief  Matrix containing number of basepairs of reference structure1 in interval [i,j] */
  unsigned int  *referenceBPs2;   /**<  @brief  Matrix containing number of basepairs of reference structure2 in interval [i,j] */
  unsigned int  *bpdist;          /**<  @brief  Matrix containing base pair distance of reference structure 1 and 2 on interval [i,j] */

  unsigned int  *mm1;             /**<  @brief  Maximum matching matrix, reference struct 1 disallowed */
  unsigned int  *mm2;             /**<  @brief  Maximum matching matrix, reference struct 2 disallowed */

  /**
   *  @name Additional data fields for local folding
   *
   *  These data fields are typically populated with meaningful data only if used in the context of local folding
   *  @{
   */
  int   window_size;              /**<  @brief  window size for local folding sliding window approach */
  char  **ptype_local;            /**<  @brief  Pair type array (for local folding) */

};



typedef struct vrna_fc_s vrna_fold_compound_t;





typedef void (*vrna_grammar_cond_f)(vrna_fold_compound_t *fc,
                                     unsigned char        stage,
                                     void                 *data);

typedef int (*vrna_grammar_rule_f)(vrna_fold_compound_t  *fc,
                                    int                   i,
                                    int                   j,
                                    void                  *data);


typedef void (*vrna_grammar_rule_f_aux)(vrna_fold_compound_t  *fc,
                                    int                   i,
                                    int                   j,
                                    void                  *data);


typedef FLT_OR_DBL (*vrna_grammar_rule_f_exp)(vrna_fold_compound_t *fc,
                                               int                  i,
                                               int                  j,
                                               void                 *data);


typedef void (*vrna_grammar_rule_f_aux_exp)(vrna_fold_compound_t *fc,
                                               int                  i,
                                               int                  j,
                                               void                 *data);


typedef void (*vrna_grammar_data_free_f)(void *data);



struct vrna_gr_aux_s {
  vrna_grammar_cond_f       cb_proc; /**< @brief A callback for pre- and post-processing of auxiliary grammar rules */

  vrna_grammar_rule_f       cb_aux_f;
  vrna_grammar_rule_f       cb_aux_c;
  vrna_grammar_rule_f       cb_aux_m;
  vrna_grammar_rule_f       cb_aux_m1;
  vrna_grammar_rule_f_aux       cb_aux;

  vrna_grammar_rule_f_exp   cb_aux_exp_f;
  vrna_grammar_rule_f_exp   cb_aux_exp_c;
  vrna_grammar_rule_f_exp     cb_aux_exp_m;
  vrna_grammar_rule_f_exp     cb_aux_exp_m1;
  vrna_grammar_rule_f_aux_exp   cb_aux_exp;

  void                        *data;
  vrna_grammar_data_free_f  free_data;
};


typedef struct vrna_gr_aux_s vrna_gr_aux_t;


typedef enum {
  VRNA_SC_DEFAULT = 0,  /**<  @brief  Default Soft Constraints */
  VRNA_SC_WINDOW = 1   /**<  @brief  Soft Constraints suitable for local structure prediction using
                     *    window approach.
                     *    @see    vrna_mfe_window(), vrna_mfe_window_zscore(), pfl_fold()
                     */
} vrna_sc_type_e;

typedef struct {
  unsigned int  interval_start;
  unsigned int  interval_end;
  int           e;
} vrna_sc_bp_storage_t;

typedef int (*vrna_sc_f)(int            i,
                         int            j,
                         int            k,
                         int            l,
                         unsigned char  d,
                         void           *data);

struct vrna_basepair_s {
  int i;
  int j;
};

typedef struct vrna_basepair_s vrna_basepair_t;

typedef vrna_basepair_t *(*vrna_sc_bt_f)(int            i,
                                         int            j,
                                         int            k,
                                         int            l,
                                         unsigned char  d,
                                         void           *data);

typedef FLT_OR_DBL (*vrna_sc_exp_f)(int           i,
                                    int           j,
                                    int           k,
                                    int           l,
                                    unsigned char d,
                                    void          *data);


typedef int (*vrna_auxdata_prepare_f)(vrna_fold_compound_t  *fc,
                                      void                  *data,
                                      unsigned int          event,
                                      void                  *event_data);




struct vrna_sc_s {
  const vrna_sc_type_e  type;
  unsigned int          n;

  unsigned char         state;

  int                   **energy_up;      /**<  @brief Energy contribution for stretches of unpaired nucleotides */
  FLT_OR_DBL            **exp_energy_up;  /**<  @brief Boltzmann Factors of the energy contributions for unpaired sequence stretches */

  int                   *up_storage;      /**<  @brief  Storage container for energy contributions per unpaired nucleotide */
  vrna_sc_bp_storage_t  **bp_storage;     /**<  @brief  Storage container for energy contributions per base pair */

#ifndef VRNA_DISABLE_C11_FEATURES
  /* C11 support for unnamed unions/structs */
  union {
    struct {
#endif
  int *energy_bp;                               /**<  @brief Energy contribution for base pairs */
  FLT_OR_DBL *exp_energy_bp;                    /**<  @brief Boltzmann Factors of the energy contribution for base pairs */
#ifndef VRNA_DISABLE_C11_FEATURES
  /* C11 support for unnamed unions/structs */
};
struct {
#endif
  int         **energy_bp_local;                        /**<  @brief Energy contribution for base pairs (sliding window approach) */
  FLT_OR_DBL  **exp_energy_bp_local;                    /**<  @brief Boltzmann Factors of the energy contribution for base pairs (sliding window approach) */
#ifndef VRNA_DISABLE_C11_FEATURES
  /* C11 support for unnamed unions/structs */
};
};
#endif

  int           *energy_stack;                    /**<  @brief Pseudo Energy contribution per base pair involved in a stack */
  FLT_OR_DBL    *exp_energy_stack;                /**<  @brief Boltzmann weighted pseudo energy contribution per nucleotide involved in a stack */

  /* generic soft contraints below */
  vrna_sc_f     f;            /**<  @brief  A function pointer used for pseudo
                               *            energy contribution in MFE calculations
                               *    @see    vrna_sc_add_f()
                               */

  vrna_sc_bt_f  bt;           /**<  @brief  A function pointer used to obtain backtraced
                               *            base pairs in loop regions that were altered
                               *            by soft constrained pseudo energy contributions
                               *    @see    vrna_sc_add_bt()
                               */

  vrna_sc_exp_f exp_f;        /**<  @brief  A function pointer used for pseudo energy
                               *            contribution boltzmann factors in PF
                               *            calculations
                               *    @see    vrna_sc_add_exp_f()
                               */

  void                *data;  /**<  @brief  A pointer to the data object provided for
                               *            for pseudo energy contribution functions of the
                               *            generic soft constraints feature
                               */

  vrna_auxdata_prepare_f  prepare_data;
  vrna_auxdata_free_f     free_data;
};


typedef struct  vrna_sc_s vrna_sc_t;




struct vrna_dimer_pf_s {
  /* free energies for: */
  double  F0AB; /**< @brief Null model without DuplexInit */
  double  FAB;  /**< @brief all states with DuplexInit correction */
  double  FcAB; /**< @brief true hybrid states only */
  double  FA;   /**< @brief monomer A */
  double  FB;   /**< @brief monomer B */
};

typedef struct vrna_dimer_pf_s vrna_dimer_pf_t;

struct vrna_dimer_conc_s {
  double  Ac_start;   /**< @brief start concentration A */
  double  Bc_start;   /**< @brief start concentration B */
  double  ABc;        /**< @brief End concentration AB */
  double  AAc;
  double  BBc;
  double  Ac;
  double  Bc;
};



typedef struct vrna_dimer_conc_s vrna_dimer_conc_t;


#define VRNA_OPTION_PF (1 << 1)



PUBLIC int
vrna_gr_set_aux_exp_c(vrna_fold_compound_t      *fc,
                      vrna_grammar_rule_f_exp cb);

typedef struct {
  FLT_OR_DBL  *prm_l;
  FLT_OR_DBL  *prm_l1;
  FLT_OR_DBL  *prml;

  int         ud_max_size;
  FLT_OR_DBL  **pmlu;
  FLT_OR_DBL  *prm_MLbu;
} helper_arrays;



typedef struct {
  struct hc_ext_def_dat     hc_dat_ext;
  vrna_hc_eval_f hc_eval_ext;

  struct hc_hp_def_dat      hc_dat_hp;
  vrna_hc_eval_f hc_eval_hp;

  struct hc_int_def_dat     hc_dat_int;
  eval_hc                   hc_eval_int;

  struct hc_mb_def_dat      hc_dat_mb;
  vrna_hc_eval_f hc_eval_mb;

  struct sc_ext_exp_dat     sc_wrapper_ext;
  struct sc_hp_exp_dat      sc_wrapper_hp;
  struct sc_int_exp_dat     sc_wrapper_int;
  struct sc_mb_exp_dat      sc_wrapper_mb;
} constraints_helper;

PRIVATE helper_arrays *
get_ml_helper_arrays(vrna_fold_compound_t *fc);

PRIVATE constraints_helper *
get_constraints_helper(vrna_fold_compound_t *fc);

PRIVATE void
compute_bpp_internal(vrna_fold_compound_t *fc,
                     int                  l,
                     vrna_ep_t            **bp_correction,
                     int                  *corr_cnt,
                     int                  *corr_size,
                     FLT_OR_DBL           *Qmax,
                     int                  *ov,
                     constraints_helper   *constraints);


PRIVATE void
compute_bpp_multibranch(vrna_fold_compound_t  *fc,
                        int                   l,
                        helper_arrays         *ml_helpers,
                        FLT_OR_DBL            *Qmax,
                        int                   *ov,
                        constraints_helper    *constraints);

PRIVATE void
compute_bpp_external(vrna_fold_compound_t *fc,
                     constraints_helper   *constraints);



PRIVATE void
multistrand_update_Y5(vrna_fold_compound_t  *fc,
                      int                   l,
                      FLT_OR_DBL            *Y5,
                      FLT_OR_DBL            **Y5p,
                      constraints_helper    *constraints);



PRIVATE void
multistrand_update_Y3(vrna_fold_compound_t  *fc,
                      int                   l,
                      FLT_OR_DBL            **Y3,
                      FLT_OR_DBL            **Y3p,
                      constraints_helper    *constraints);



PRIVATE void
multistrand_contrib(vrna_fold_compound_t  *fc,
                    int                   l,
                    FLT_OR_DBL            *Y5,
                    FLT_OR_DBL            **Y3,
                    constraints_helper    *constraints,
                    FLT_OR_DBL            *Qmax,
                    int                   *ov);




PUBLIC char *
vrna_db_from_probs(const FLT_OR_DBL *p,
                   unsigned int     length);



PRIVATE void
free_ml_helper_arrays(helper_arrays *ml_helpers);


PRIVATE void
free_constraints_helper(constraints_helper *helper);



PUBLIC vrna_dimer_pf_t
vrna_pf_dimer(vrna_fold_compound_t  *fc,
              char                  *structure);