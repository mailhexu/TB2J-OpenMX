import cffi

ffi = cffi.FFI()
ffi.set_source(
    "TB2J_OpenMX.cmod._scfout_parser", r"""

/**********************************************************************
  read_scfout.c:

  read_scfout.c is a subroutine to read a binary file,
  filename.scfout.
 
  Log of read_scfout.c:

  02/Jul./2003  Released by T. Ozaki
  22/Aug./2011  Released by H. Kotaka
  17/Sep./2019  Modified by N. Yamaguchi

 ***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAX_LINE_SIZE 256
#define fp_bsize         1048576     /* buffer size for setvbuf */
#define SCFOUT_VERSION 3

void free_scfout();

int T_NumOrbs; /* total number of orbs*/
int *MP ;         /* Full system version of MP */
double ****HR;
double ****HR_imag;
double ***SR;




/*******************************************************
  int atomnun;
  the number of total atoms
 *******************************************************/
int atomnum;   

/*******************************************************
  int Catomnun;
  the number of atoms in the central region
 *******************************************************/
int Catomnum;

/*******************************************************
  int Latomnun;
  the number of atoms in the left lead
 *******************************************************/
int Latomnum;

/*******************************************************
  int Ratomnun;
  the number of atoms in the left lead
 *******************************************************/
int Ratomnum;   

/*******************************************************
  int SpinP_switch;
0: non-spin polarized 
1: spin polarized
3: non-collinear calculation
 *******************************************************/
int SpinP_switch;

/*******************************************************
  int TCpyCell;
  the total number of periodic cells
 *******************************************************/
int TCpyCell;

/*******************************************************
  int Solver;
  method for solving eigenvalue problem
 *******************************************************/
int Solver;

/*******************************************************
  double ChemP;
  chemical potential
 *******************************************************/
double ChemP;

/*******************************************************
  int Valence_Electrons;
  total number of valence electrons
 *******************************************************/
int Valence_Electrons;

/*******************************************************
  double Total_SpinS;
  total value of Spin (2*Total_SpinS = muB)
 *******************************************************/
double Total_SpinS;

/*******************************************************
  double E_Temp;
  electronic temperature
 *******************************************************/
double E_Temp;

/*******************************************************
  int *Total_NumOrbs; 
  the number of atomic orbitals in each atom
size: Total_NumOrbs[atomnum+1]
 *******************************************************/
int *Total_NumOrbs;

/*******************************************************
  int *FNAN; 
  the number of first neighboring atoms of each atom
size: FNAN[atomnum+1]
 *******************************************************/
int *FNAN;

/*******************************************************
  int **natn; 
  grobal index of neighboring atoms of an atom ct_AN
size: natn[atomnum+1][FNAN[ct_AN]+1]
 *******************************************************/
int **natn;

/*******************************************************
  int **ncn; 
  grobal index for cell of neighboring atoms of
  an atom ct_AN
size: ncn[atomnum+1][FNAN[ct_AN]+1]
 *******************************************************/
int **ncn;

/*******************************************************
  double **atv;
  x,y,and z-components of translation vector of  
  periodically copied cells
size: atv[TCpyCell+1][4];
 *******************************************************/
double **atv;

/*******************************************************
  int **atv_ijk;
  i,j,and k number of periodically copied cells
size: atv_ijk[TCpyCell+1][4];
 *******************************************************/
int **atv_ijk;

/*******************************************************
  double tv[4][4];
  unit cell vectors in Bohr
 *******************************************************/
double tv[4][4];

/*******************************************************
  double rtv[4][4]:
  reciprocal unit cell vectors in Bohr^{-1}

note:
tv_i \dot rtv_j = 2PI * Kronecker's delta_{ij}
 *******************************************************/
double rtv[4][4];

/*******************************************************
  double Gxyz[atomnum+1][60];
  atomic coordinates in Bohr
 *******************************************************/
double **Gxyz;

/*******************************************************
  double *****Hks;
  Kohn-Sham matrix elements of basis orbitals
size: Hks[SpinP_switch+1]
[atomnum+1]
[FNAN[ct_AN]+1]
[Total_NumOrbs[ct_AN]]
[Total_NumOrbs[h_AN]] 
 *******************************************************/
double *****Hks;

/*******************************************************
  double *****iHks;
  imaginary Kohn-Sham matrix elements of basis orbitals
  for alpha-alpha, beta-beta, and alpha-beta spin matrices
  of which contributions come from spin-orbit coupling 
  and Hubbard U effective potential.
size: iHks[3]
[atomnum+1]
[FNAN[ct_AN]+1]
[Total_NumOrbs[ct_AN]]
[Total_NumOrbs[h_AN]] 
 *******************************************************/
double *****iHks;

/*******************************************************
  double ****OLP;
  overlap matrix
size: OLP[atomnum+1]
[FNAN[ct_AN]+1]
[Total_NumOrbs[ct_AN]]
[Total_NumOrbs[h_AN]]
 *******************************************************/
double ****OLP;

/*******************************************************
  double *****DM;
  density matrix
size: DM[SpinP_switch+1]
[atomnum+1]
[FNAN[ct_AN]+1]
[Total_NumOrbs[ct_AN]]
[Total_NumOrbs[h_AN]]
 *******************************************************/
double *****DM;

/*******************************************************
  double *****iDM;
  density matrix
size: iDM[2]
[atomnum+1]
[FNAN[ct_AN]+1]
[Total_NumOrbs[ct_AN]]
[Total_NumOrbs[h_AN]]
 *******************************************************/
double *****iDM;

/*******************************************************
  double dipole_moment_core[4];
 *******************************************************/
double dipole_moment_core[4];

/*******************************************************
  int version; (added by N. Yamaguchi)
 *******************************************************/
int version;

/*******************************************************
  int order_max; (added by N. Yamaguchi for HWC)
 *******************************************************/
int order_max;

/*******************************************************
  double *cc_vec; (added by N. Yamaguchi for HWC)
 *******************************************************/
double *cc_vec;

/*******************************************************
  double ******OLPpo; (added by N. Yamaguchi for HWC)
 *******************************************************/
double ******OLPpo;

/*******************************************************
  double *****OLPmo; 
 *******************************************************/
double *****OLPmo;

/*******************************************************
  double dipole_moment_background[4];
 *******************************************************/
double dipole_moment_background[4];




/* Added by N. Yamaguchi ***/
/* ***
 Note: Use "FREAD" instead of "fread" to avoid the mismatch of endianness. 
 Every "fread" in "read_scfout.c" was replaced with "FREAD".
 *** */
#define LATEST_VERSION 3
#define FREAD(POINTER, SIZE, NUMBER, FILE_POINTER)\
  do {\
    fread(POINTER, SIZE, NUMBER, FILE_POINTER);\
    if (conversionSwitch){\
      int dATA;\
      for (dATA=0; dATA<NUMBER; dATA++){\
	char *out=(char*)(POINTER+dATA);\
	int bYTE;\
	for (bYTE=0; bYTE<SIZE/2; bYTE++){\
	  char tmp=out[bYTE];\
	  out[bYTE]=out[SIZE-bYTE-1];\
	  out[SIZE-bYTE-1]=tmp;\
	}\
      }\
    }\
  } while (0)
/* ***/

static void Input( FILE *fp );



void read_scfout(char *argv[])
{
  static FILE *fp;

  char buf[fp_bsize];          /* setvbuf */


  /* Added by N. Yamaguchi ***/
  /* Modified by N. Yamaguchi ***/

  if ((fp = fopen(argv[1],"r")) != NULL){

    setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */

      printf("\nRead the scfout file (%s)\n",argv[1]);fflush(stdout);

    Input(fp);
    fclose(fp);
  }

  else {
    printf("Failure of reading the scfout file (%s).\n",argv[1]);fflush(stdout);
  }

}


void Input( FILE *fp )
{
  static int Gc_AN,ct_AN,h_AN,i,j,can,Gh_AN;
  static int wan1,wan2,TNO1,TNO2,spin,Rn,num_lines;
  static int k,q_AN,Gq_AN;
  static int i_vec[20],*p_vec;
  static double d_vec[20];
  static char makeinp[100];
  static char strg[MAX_LINE_SIZE];
  int direction, order;
  FILE *fp_makeinp;
  char buf[fp_bsize];          /* setvbuf */

  /****************************************************
    atomnum
    spinP_switch
    version (added by N. Yamaguchi)
   ****************************************************/

  fread(i_vec,sizeof(int),6,fp);

  /* Disabled by N. Yamaguchi
   * atomnum      = i_vec[0];
   * SpinP_switch = i_vec[1];
   * Catomnum =     i_vec[2];
   * Latomnum =     i_vec[3];
   * Ratomnum =     i_vec[4];
   * TCpyCell =     i_vec[5];
   */

  /* Added by N. Yamaguchi ***/
  int conversionSwitch;
  if (i_vec[1]==0 && i_vec[1]<0 || i_vec[1]>(LATEST_VERSION)*4+3){
    conversionSwitch=1;
    int i;
    for (i=0; i<6; i++){
      int value=*(i_vec+i);
      char *in=(char*)&value;
      char *out=(char*)(i_vec+i);
      int j;
      for (j=0; j<sizeof(int); j++){
	out[j]=in[sizeof(int)-j-1];
      }
    }
    if (i_vec[1]==0 && i_vec[1]<0 || i_vec[1]>(LATEST_VERSION)*4+3){
      puts("Error: Mismatch of the endianness");fflush(stdout);
    }
  } else {
    conversionSwitch=0;
  }
  /* ***/

  atomnum      = i_vec[0];

  /* Disabled by N. Yamaguchi ***
  SpinP_switch = i_vec[1];
  * ***/

  /* Added by N. Yamaguchi ***/
  SpinP_switch=i_vec[1]%4;
  version=i_vec[1]/4;


  char* openmxVersion;

  if (version==0){
    openmxVersion="3.7, 3.8 or an older distribution";
  } else if (version==1){
    openmxVersion="3.7.x (for development of HWC)";
  } else if (version==2){
    openmxVersion="3.7.x (for development of HWF)";
  } else if (version==3){
    openmxVersion="3.9";
  }

  if (version!=SCFOUT_VERSION){
      printf("The file format of the SCFOUT file:  %d\n", version);
      printf("The vesion is not supported by the current read_scfout\n");
    exit(0);
  }

      //printf("SpinP_switch %d", SpinP_switch);

    puts("***");

    printf("The file format of the SCFOUT file:  %d\n", version);
    //puts("And it supports the following functions:");
    //puts("- jx");
    //puts("- polB");
    //puts("- kSpin");
    //puts("- Z2FH");
    //puts("- calB");
    puts("***");

  /* ***/

  Catomnum =     i_vec[2];
  Latomnum =     i_vec[3];
  Ratomnum =     i_vec[4];
  TCpyCell =     i_vec[5];

  /****************************************************
    order_max (added by N. Yamaguchi for HWC)
   ****************************************************/

  FREAD(i_vec, sizeof(int), 1, fp);
  order_max=i_vec[0];

  /****************************************************
    allocation of arrays:

    double atv[TCpyCell+1][4];
  ****************************************************/

  atv = (double**)malloc(sizeof(double*)*(TCpyCell+1));
  for (Rn=0; Rn<=TCpyCell; Rn++){
    atv[Rn] = (double*)malloc(sizeof(double)*4);
  }

  /****************************************************
    read atv[TCpyCell+1][4];
   ****************************************************/

  for (Rn=0; Rn<=TCpyCell; Rn++){
    FREAD(atv[Rn],sizeof(double),4,fp);
  }

  /****************************************************
    allocation of arrays:

    int atv_ijk[TCpyCell+1][4];
   ****************************************************/

  atv_ijk = (int**)malloc(sizeof(int*)*(TCpyCell+1));
  for (Rn=0; Rn<=TCpyCell; Rn++){
    atv_ijk[Rn] = (int*)malloc(sizeof(int)*4);
  }

  /****************************************************
    read atv_ijk[TCpyCell+1][4];
   ****************************************************/

  for (Rn=0; Rn<=TCpyCell; Rn++){
    FREAD(atv_ijk[Rn],sizeof(int),4,fp);
  }

  /****************************************************
    allocation of arrays:

    int Total_NumOrbs[atomnum+1];
    int FNAN[atomnum+1];
   ****************************************************/

  Total_NumOrbs = (int*)malloc(sizeof(int)*(atomnum+1));
  FNAN = (int*)malloc(sizeof(int)*(atomnum+1));

  /****************************************************
    the number of orbitals in each atom
   ****************************************************/

  p_vec = (int*)malloc(sizeof(int)*atomnum);
  FREAD(p_vec,sizeof(int),atomnum,fp);
  Total_NumOrbs[0] = 1;
  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    Total_NumOrbs[ct_AN] = p_vec[ct_AN-1];
  }
  free(p_vec);

  /****************************************************
    FNAN[]:
    the number of first nearest neighbouring atoms
   ****************************************************/

  p_vec = (int*)malloc(sizeof(int)*atomnum);
  FREAD(p_vec,sizeof(int),atomnum,fp);
  FNAN[0] = 0;
  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    FNAN[ct_AN] = p_vec[ct_AN-1];
  }
  free(p_vec);

  /****************************************************
    allocation of arrays:

    int natn[atomnum+1][FNAN[ct_AN]+1];
    int ncn[atomnum+1][FNAN[ct_AN]+1];
   ****************************************************/

  natn = (int**)malloc(sizeof(int*)*(atomnum+1));
  for (ct_AN=0; ct_AN<=atomnum; ct_AN++){
    natn[ct_AN] = (int*)malloc(sizeof(int)*(FNAN[ct_AN]+1));
  }

  ncn = (int**)malloc(sizeof(int*)*(atomnum+1));
  for (ct_AN=0; ct_AN<=atomnum; ct_AN++){
    ncn[ct_AN] = (int*)malloc(sizeof(int)*(FNAN[ct_AN]+1));
  }

  /****************************************************
    natn[][]:
    grobal index of neighboring atoms of an atom ct_AN
   ****************************************************/

  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    FREAD(natn[ct_AN],sizeof(int),FNAN[ct_AN]+1,fp);
  }

  /****************************************************
    ncn[][]:
    grobal index for cell of neighboring atoms
    of an atom ct_AN
   ****************************************************/

  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    FREAD(ncn[ct_AN],sizeof(int),FNAN[ct_AN]+1,fp);
  }

  /****************************************************
    tv[4][4]:
    unit cell vectors in Bohr
   ****************************************************/

  FREAD(tv[1],sizeof(double),4,fp);
  FREAD(tv[2],sizeof(double),4,fp);
  FREAD(tv[3],sizeof(double),4,fp);

  /****************************************************
    rtv[4][4]:
    unit cell vectors in Bohr
   ****************************************************/

  FREAD(rtv[1],sizeof(double),4,fp);
  FREAD(rtv[2],sizeof(double),4,fp);
  FREAD(rtv[3],sizeof(double),4,fp);

  /****************************************************
    Gxyz[][1-3]:
    atomic coordinates in Bohr
   ****************************************************/

  Gxyz = (double**)malloc(sizeof(double*)*(atomnum+1));
  for (i=0; i<(atomnum+1); i++){
    Gxyz[i] = (double*)malloc(sizeof(double)*60);
  }

  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    FREAD(Gxyz[ct_AN],sizeof(double),4,fp);
  }

  /****************************************************
    allocation of arrays:

    Kohn-Sham Hamiltonian

    dooble Hks[SpinP_switch+1]
    [atomnum+1]
    [FNAN[ct_AN]+1]
    [Total_NumOrbs[ct_AN]]
    [Total_NumOrbs[h_AN]];

    Overlap matrix

    dooble OLP[atomnum+1]
    [FNAN[ct_AN]+1]
    [Total_NumOrbs[ct_AN]]
    [Total_NumOrbs[h_AN]]; 

    Overlap matrix with position operator x, y, z

    double ******OLPpo;
    [3]
    [1]
    [atomnum+1]
    [FNAN[ct_AN]+1]
    [Total_NumOrbs[ct_AN]]
    [Total_NumOrbs[h_AN]]; 

    Overlap matrix with momentum operator px, py, pz

    double *****OLPmo;
    [3]
    [atomnum+1]
    [FNAN[ct_AN]+1]
    [Total_NumOrbs[ct_AN]]
    [Total_NumOrbs[h_AN]]; 

    Density matrix

    dooble DM[SpinP_switch+1]
    [atomnum+1]
    [FNAN[ct_AN]+1]
    [Total_NumOrbs[ct_AN]]
    [Total_NumOrbs[h_AN]];
   ****************************************************/

  Hks = (double*****)malloc(sizeof(double****)*(SpinP_switch+1));
  for (spin=0; spin<=SpinP_switch; spin++){

    Hks[spin] = (double****)malloc(sizeof(double***)*(atomnum+1));
    for (ct_AN=0; ct_AN<=atomnum; ct_AN++){
      TNO1 = Total_NumOrbs[ct_AN];
      Hks[spin][ct_AN] = (double***)malloc(sizeof(double**)*(FNAN[ct_AN]+1));
      for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++){
	Hks[spin][ct_AN][h_AN] = (double**)malloc(sizeof(double*)*TNO1);

	if (ct_AN==0){
	  TNO2 = 1;
	}
	else{
	  Gh_AN = natn[ct_AN][h_AN];
	  TNO2 = Total_NumOrbs[Gh_AN];
	}
	for (i=0; i<TNO1; i++){
	  Hks[spin][ct_AN][h_AN][i] = (double*)malloc(sizeof(double)*TNO2);
	}
      }
    }
  }

  //printf("Hks col alloc OK!\n");

  /* Added by N. Yamaguchi ***/
  if (SpinP_switch==3){
    /* ***/

    iHks = (double*****)malloc(sizeof(double****)*3);
    for (spin=0; spin<3; spin++){

      iHks[spin] = (double****)malloc(sizeof(double***)*(atomnum+1));
      for (ct_AN=0; ct_AN<=atomnum; ct_AN++){
	TNO1 = Total_NumOrbs[ct_AN];
	iHks[spin][ct_AN] = (double***)malloc(sizeof(double**)*(FNAN[ct_AN]+1));
	for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++){
	  iHks[spin][ct_AN][h_AN] = (double**)malloc(sizeof(double*)*TNO1);

	  if (ct_AN==0){
	    TNO2 = 1;
	  }
	  else{
	    Gh_AN = natn[ct_AN][h_AN];
	    TNO2 = Total_NumOrbs[Gh_AN];
	  }
	  for (i=0; i<TNO1; i++){
	    iHks[spin][ct_AN][h_AN][i] = (double*)malloc(sizeof(double)*TNO2);
	    for (j=0; j<TNO2; j++) iHks[spin][ct_AN][h_AN][i][j] = 0.0;
	  }
	}
      }
    }

  //printf("Hks non-col alloc OK!\n");

    /* Added by N. Yamaguchi ***/
  }
  /* ***/

  OLP = (double****)malloc(sizeof(double***)*(atomnum+1));
  for (ct_AN=0; ct_AN<=atomnum; ct_AN++){
    TNO1 = Total_NumOrbs[ct_AN];
    OLP[ct_AN] = (double***)malloc(sizeof(double**)*(FNAN[ct_AN]+1));
    for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++){
      OLP[ct_AN][h_AN] = (double**)malloc(sizeof(double*)*TNO1);

      if (ct_AN==0){
	TNO2 = 1;
      }
      else{
	Gh_AN = natn[ct_AN][h_AN];
	TNO2 = Total_NumOrbs[Gh_AN];
      }
      for (i=0; i<TNO1; i++){
	OLP[ct_AN][h_AN][i] = (double*)malloc(sizeof(double)*TNO2);
      }
    }
  }

  //printf("Overlap matrix read OK!\n");

  /* Added by N. Yamaguchi ***/

  OLPpo=(double******)malloc(sizeof(double*****)*3);
  for (direction=0; direction<3; direction++){
    OLPpo[direction]=(double*****)malloc(sizeof(double****)*order_max);
    for (order=0; order<order_max; order++){
      OLPpo[direction][order]=(double****)malloc(sizeof(double***)*(atomnum+1));
      for (ct_AN=0; ct_AN<=atomnum; ct_AN++){
	TNO1 = Total_NumOrbs[ct_AN];
	OLPpo[direction][order][ct_AN] = (double***)malloc(sizeof(double**)*(FNAN[ct_AN]+1));
	for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++){
	  OLPpo[direction][order][ct_AN][h_AN] = (double**)malloc(sizeof(double*)*TNO1);

	  if (ct_AN==0){
	    TNO2 = 1;
	  }
	  else{
	    Gh_AN = natn[ct_AN][h_AN];
	    TNO2 = Total_NumOrbs[Gh_AN];
	  }
	  for (i=0; i<TNO1; i++){
	    OLPpo[direction][order][ct_AN][h_AN][i] = (double*)malloc(sizeof(double)*TNO2);
	  }
	}
      }
    }
  }

  OLPmo=(double*****)malloc(sizeof(double****)*3);
  for (direction=0; direction<3; direction++){
    OLPmo[direction]=(double****)malloc(sizeof(double***)*(atomnum+1));
    for (ct_AN=0; ct_AN<=atomnum; ct_AN++){
      TNO1 = Total_NumOrbs[ct_AN];
      OLPmo[direction][ct_AN] = (double***)malloc(sizeof(double**)*(FNAN[ct_AN]+1));
      for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++){
	OLPmo[direction][ct_AN][h_AN] = (double**)malloc(sizeof(double*)*TNO1);

	if (ct_AN==0){
	  TNO2 = 1;
	}
	else{
	  Gh_AN = natn[ct_AN][h_AN];
	  TNO2 = Total_NumOrbs[Gh_AN];
	}
	for (i=0; i<TNO1; i++){
	  OLPmo[direction][ct_AN][h_AN][i] = (double*)malloc(sizeof(double)*TNO2);
	}
      }
    }
  }

  DM = (double*****)malloc(sizeof(double****)*(SpinP_switch+1));
  for (spin=0; spin<=SpinP_switch; spin++){

    DM[spin] = (double****)malloc(sizeof(double***)*(atomnum+1));
    for (ct_AN=0; ct_AN<=atomnum; ct_AN++){
      TNO1 = Total_NumOrbs[ct_AN];
      DM[spin][ct_AN] = (double***)malloc(sizeof(double**)*(FNAN[ct_AN]+1));
      for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++){
	DM[spin][ct_AN][h_AN] = (double**)malloc(sizeof(double*)*TNO1);

	if (ct_AN==0){ 
	  TNO2 = 1;
	}
	else{ 
	  Gh_AN = natn[ct_AN][h_AN];
	  TNO2 = Total_NumOrbs[Gh_AN];
	}
	for (i=0; i<TNO1; i++){
	  DM[spin][ct_AN][h_AN][i] = (double*)malloc(sizeof(double)*TNO2);
	}
      }
    }
  }

  iDM = (double*****)malloc(sizeof(double****)*2);
  for (spin=0; spin<2; spin++){

    iDM[spin] = (double****)malloc(sizeof(double***)*(atomnum+1));
    for (ct_AN=0; ct_AN<=atomnum; ct_AN++){
      TNO1 = Total_NumOrbs[ct_AN];
      iDM[spin][ct_AN] = (double***)malloc(sizeof(double**)*(FNAN[ct_AN]+1));
      for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++){
	iDM[spin][ct_AN][h_AN] = (double**)malloc(sizeof(double*)*TNO1);

	if (ct_AN==0){ 
	  TNO2 = 1;
	}
	else{ 
	  Gh_AN = natn[ct_AN][h_AN];
	  TNO2 = Total_NumOrbs[Gh_AN];
	}
	for (i=0; i<TNO1; i++){
	  iDM[spin][ct_AN][h_AN][i] = (double*)malloc(sizeof(double)*TNO2);
	}
      }
    }
  }

  /****************************************************
    Hamiltonian matrix
   ****************************************************/

  for (spin=0; spin<=SpinP_switch; spin++){
    for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      TNO1 = Total_NumOrbs[ct_AN];
      for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++){
	Gh_AN = natn[ct_AN][h_AN];
	TNO2 = Total_NumOrbs[Gh_AN];
	for (i=0; i<TNO1; i++){
	  FREAD(Hks[spin][ct_AN][h_AN][i],sizeof(double),TNO2,fp);
	}
      }
    }
  }

  /*********************************************************
   iHks:
   imaginary Kohn-Sham matrix elements of basis orbitals
   for alpha-alpha, beta-beta, and alpha-beta spin matrices
   of which contributions come from spin-orbit coupling 
   and Hubbard U effective potential.
  **********************************************************/

  if (SpinP_switch==3){
    for (spin=0; spin<3; spin++){
      for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
	TNO1 = Total_NumOrbs[ct_AN];
	for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++){
	  Gh_AN = natn[ct_AN][h_AN];
	  TNO2 = Total_NumOrbs[Gh_AN];
	  for (i=0; i<TNO1; i++){
	    FREAD(iHks[spin][ct_AN][h_AN][i],sizeof(double),TNO2,fp);
	  }
	}
      }
    }
  }

  /****************************************************
      Overlap matrix
   ****************************************************/

  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    TNO1 = Total_NumOrbs[ct_AN];
    for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++){
      Gh_AN = natn[ct_AN][h_AN];
      TNO2 = Total_NumOrbs[Gh_AN];
      for (i=0; i<TNO1; i++){
	FREAD(OLP[ct_AN][h_AN][i],sizeof(double),TNO2,fp);
      }
    }
  }

  /****************************************************
    Overlap matrix with position operator
    (added by N. Yamaguchi for HWC)
   ****************************************************/

  /* Added by N. Yamaguchi ***/

  for (direction=0; direction<3; direction++){
    for (order=0; order<order_max; order++){
      for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
	TNO1 = Total_NumOrbs[ct_AN];
	for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++){
	  Gh_AN = natn[ct_AN][h_AN];
	  TNO2 = Total_NumOrbs[Gh_AN];
	  for (i=0; i<TNO1; i++){
	    FREAD(OLPpo[direction][order][ct_AN][h_AN][i],sizeof(double),TNO2,fp);
	  }
	}
      }
    }
  }

  /****************************************************
    Overlap matrix with momentum operator
   ****************************************************/

  for (direction=0; direction<3; direction++){
    for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      TNO1 = Total_NumOrbs[ct_AN];
      for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++){
	Gh_AN = natn[ct_AN][h_AN];
	TNO2 = Total_NumOrbs[Gh_AN];
	for (i=0; i<TNO1; i++){
	  FREAD(OLPmo[direction][ct_AN][h_AN][i],sizeof(double),TNO2,fp);
	}
      }
    }
  }

  /****************************************************
             Density matrix: DM and iDM
   ****************************************************/

  for (spin=0; spin<=SpinP_switch; spin++){
    for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      TNO1 = Total_NumOrbs[ct_AN];
      for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++){
	Gh_AN = natn[ct_AN][h_AN];
	TNO2 = Total_NumOrbs[Gh_AN];
	for (i=0; i<TNO1; i++){
	  FREAD(DM[spin][ct_AN][h_AN][i],sizeof(double),TNO2,fp);
	}
      }
    }
  }

  for (spin=0; spin<2; spin++){
    for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
      TNO1 = Total_NumOrbs[ct_AN];
      for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++){
	Gh_AN = natn[ct_AN][h_AN];
	TNO2 = Total_NumOrbs[Gh_AN];
	for (i=0; i<TNO1; i++){
	  FREAD(iDM[spin][ct_AN][h_AN][i],sizeof(double),TNO2,fp);
	}
      }
    }
  }

  /****************************************************
    Solver
   ****************************************************/

  FREAD(i_vec,sizeof(int),1,fp);
  Solver = i_vec[0];

  /****************************************************
    ChemP
    Temp
   ****************************************************/

  FREAD(d_vec,sizeof(double),10,fp);
  ChemP  = d_vec[0];
  E_Temp = d_vec[1];
  dipole_moment_core[1] = d_vec[2]; 
  dipole_moment_core[2] = d_vec[3]; 
  dipole_moment_core[3] = d_vec[4]; 
  dipole_moment_background[1] = d_vec[5]; 
  dipole_moment_background[2] = d_vec[6]; 
  dipole_moment_background[3] = d_vec[7]; 
  Valence_Electrons = d_vec[8]; 
  Total_SpinS = d_vec[9];

  /****************************************************
    input file
   ****************************************************/

  FREAD(i_vec, sizeof(int), 1, fp);
  num_lines = i_vec[0];

  sprintf(makeinp,"temporal_12345.input");

  if ((fp_makeinp = fopen(makeinp,"w")) != NULL){

    setvbuf(fp_makeinp,buf,_IOFBF,fp_bsize);  /* setvbuf */
    for (i=1; i<=num_lines; i++){
      FREAD(strg, sizeof(char), MAX_LINE_SIZE, fp);
      fprintf(fp_makeinp,"%s",strg);
    }
    fclose(fp_makeinp);
  }
}


void free_scfout()
{
  int Rn,ct_AN,i,h_AN,spin,TNO1,order,direction;

  for (Rn=0; Rn<=TCpyCell; Rn++){
    free(atv[Rn]);
  }
  free(atv);

  for (Rn=0; Rn<=TCpyCell; Rn++){
    free(atv_ijk[Rn]);
  }
  free(atv_ijk);

  for (ct_AN=0; ct_AN<=atomnum; ct_AN++){
    free(natn[ct_AN]);
  }
  free(natn);

  for (ct_AN=0; ct_AN<=atomnum; ct_AN++){
    free(ncn[ct_AN]);
  }
  free(ncn);

  for (i=0; i<(atomnum+1); i++){
    free(Gxyz[i]);
  }
  free(Gxyz);

  for (spin=0; spin<=SpinP_switch; spin++){

    for (ct_AN=0; ct_AN<=atomnum; ct_AN++){
      TNO1 = Total_NumOrbs[ct_AN];
      for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++){
	for (i=0; i<TNO1; i++){
	  free(Hks[spin][ct_AN][h_AN][i]);
	}
        free(Hks[spin][ct_AN][h_AN]);
      }
      free(Hks[spin][ct_AN]);
    }
    free(Hks[spin]);
  }
  free(Hks);

  if (SpinP_switch==3){

    for (spin=0; spin<3; spin++){

      for (ct_AN=0; ct_AN<=atomnum; ct_AN++){
	TNO1 = Total_NumOrbs[ct_AN];
	for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++){
	  for (i=0; i<TNO1; i++){
	    free(iHks[spin][ct_AN][h_AN][i]);
	  }
          free(iHks[spin][ct_AN][h_AN]);
	}
        free(iHks[spin][ct_AN]);
      }
      free(iHks[spin]);
    }
    free(iHks);
  }

  for (ct_AN=0; ct_AN<=atomnum; ct_AN++){
    TNO1 = Total_NumOrbs[ct_AN];
    for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++){
      for (i=0; i<TNO1; i++){
	free(OLP[ct_AN][h_AN][i]);
      }
      free(OLP[ct_AN][h_AN]);
    }
    free(OLP[ct_AN]);
  }
  free(OLP);

  for (direction=0; direction<3; direction++){
    for (order=0; order<order_max; order++){
      for (ct_AN=0; ct_AN<=atomnum; ct_AN++){
	TNO1 = Total_NumOrbs[ct_AN];
	for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++){
	  for (i=0; i<TNO1; i++){
	    free(OLPpo[direction][order][ct_AN][h_AN][i]);
	  }
          free(OLPpo[direction][order][ct_AN][h_AN]);
	}
        free(OLPpo[direction][order][ct_AN]);
      }
      free(OLPpo[direction][order]);
    }
    free(OLPpo[direction]);
  }
  free(OLPpo);

  for (direction=0; direction<3; direction++){
    for (ct_AN=0; ct_AN<=atomnum; ct_AN++){
      TNO1 = Total_NumOrbs[ct_AN];
      for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++){
	for (i=0; i<TNO1; i++){
	  free(OLPmo[direction][ct_AN][h_AN][i]);
	}
        free(OLPmo[direction][ct_AN][h_AN]);
      }
      free(OLPmo[direction][ct_AN]);
    }
    free(OLPmo[direction]);
  }
  free(OLPmo);

  for (spin=0; spin<=SpinP_switch; spin++){
    for (ct_AN=0; ct_AN<=atomnum; ct_AN++){
      TNO1 = Total_NumOrbs[ct_AN];
      for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++){
	for (i=0; i<TNO1; i++){
	  free(DM[spin][ct_AN][h_AN][i]);
	}
        free(DM[spin][ct_AN][h_AN]);
      }
      free(DM[spin][ct_AN]);
    }
    free(DM[spin]);
  }
  free(DM);

  for (spin=0; spin<2; spin++){
    for (ct_AN=0; ct_AN<=atomnum; ct_AN++){
      TNO1 = Total_NumOrbs[ct_AN];
      for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++){
	for (i=0; i<TNO1; i++){
	  free(iDM[spin][ct_AN][h_AN][i]);
	}
        free(iDM[spin][ct_AN][h_AN]);
      }
      free(iDM[spin][ct_AN]);
    }
    free(iDM[spin]);
  }
  free(iDM);

  free(Total_NumOrbs);
  free(FNAN);
}

void prepare_HSR(){
  static int i,j,l,n, II, JJ;
  static int ct_AN,h_AN,Gh_AN;
  static int spin,Rn, iR;
  int end_switch;    /* switch in the evalusion loop of J */


  MP = (int*)malloc(sizeof(int)*(atomnum+1));
  T_NumOrbs = 0;
  for (i=0; i<atomnum; i++){
    ct_AN = i+1;
    MP[i+1] = T_NumOrbs+1;
    T_NumOrbs = T_NumOrbs + Total_NumOrbs[ct_AN];
  }


  // allocate HR
  //[iR, spin,  iorb, jorb]
  HR = (double****)malloc(sizeof(double***)*(TCpyCell+1));
  for (iR=0; iR<TCpyCell+1; iR++){
    HR[iR]= (double***)malloc(sizeof(double**)*(SpinP_switch+1));
    for (spin=0; spin<=SpinP_switch; spin++){
      HR[iR][spin] = (double**)malloc(sizeof(double*)*(T_NumOrbs));
      for (i=0; i<T_NumOrbs; i++){
        HR[iR][spin][i] = (double*)malloc(sizeof(double)*(T_NumOrbs));
        for (j=0; j<T_NumOrbs; j++){
          HR[iR][spin][i][j]=0.0;
        }
      }
    }
  }

if(SpinP_switch==3) {
  HR_imag = (double****)malloc(sizeof(double***)*(TCpyCell+1));
  for (iR=0; iR<TCpyCell+1; iR++){
    HR_imag[iR]= (double***)malloc(sizeof(double**)*3);
    for (spin=0; spin<=2; spin++){
      HR_imag[iR][spin] = (double**)malloc(sizeof(double*)*(T_NumOrbs));
      for (i=0; i<T_NumOrbs; i++){
        HR_imag[iR][spin][i] = (double*)malloc(sizeof(double)*(T_NumOrbs));
        for (j=0; j<T_NumOrbs; j++){
          HR_imag[iR][spin][i][j]=0.0;
        }
      }
    }
  }
}


  // allocate SR
  //[iR, spin,  iorb, jorb]
  SR = (double***)malloc(sizeof(double***)*(TCpyCell+1));
  for (iR=0; iR<TCpyCell+1; iR++){
      SR[iR] = (double**)malloc(sizeof(double*)*(T_NumOrbs+1));
      for (i=0; i<=T_NumOrbs; i++){
        SR[iR][i] = (double*)malloc(sizeof(double)*(T_NumOrbs+1));
        for (j=0; j<=T_NumOrbs; j++){
          SR[iR][i][j]=0;

      }
    }
  }



// fill HR and SR
for (spin=0; spin<=SpinP_switch; spin++){
  for (II=0; II<atomnum; II++){
    for (JJ=0; JJ<atomnum; JJ++){
      ct_AN = II+1 ; //index of atom I
      // loop over NN list
      for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++){ 
        // global index of NN atom J
        Gh_AN = natn[ct_AN][h_AN]; // why not h_AN+1?
        // or if JJ+1 in NN list
        if (Gh_AN == JJ+1){
          iR = ncn[ct_AN][h_AN];
          //if(iR==0 || iR==342){
          //   printf("iR: %d\n", iR);
          //}
          //printf("Norb for atom %d: %d\n", ct_AN, Total_NumOrbs[ct_AN]);
          for (i=0; i<Total_NumOrbs[ct_AN]; i++){
            //printf("iorb: %d\n", i+MP[II+1]);
            for (j=0; j<Total_NumOrbs[Gh_AN]; j++){
              HR[iR][spin][i+MP[II+1]-1][j+MP[JJ+1]-1]
                = Hks[spin][ct_AN][h_AN][i][j];
              SR[iR][i+MP[II+1]-1][j+MP[JJ+1]-1]
                = OLP[ct_AN][h_AN][i][j];
            }
          }
        }
      }
    }
  }
 }



// fill HR_imag
if (SpinP_switch==3){
for (spin=0; spin<3; spin++){
  for (II=0; II<atomnum; II++){
    for (JJ=0; JJ<atomnum; JJ++){
      ct_AN = II+1 ; //index of atom I
      // loop over NN list
      for (h_AN=0; h_AN<=FNAN[ct_AN]; h_AN++){ 
        // global index of NN atom J
        Gh_AN = natn[ct_AN][h_AN]; // why not h_AN+1?
        // or if JJ+1 in NN list
        if (Gh_AN == JJ+1){
          iR = ncn[ct_AN][h_AN];
          for (i=0; i<Total_NumOrbs[ct_AN]; i++){
            for (j=0; j<Total_NumOrbs[Gh_AN]; j++){
              HR_imag[iR][spin][i+MP[II+1]-1][j+MP[JJ+1]-1]
                = iHks[spin][ct_AN][h_AN][i][j];
            }
          }
        }
      }
    }
  }
 }
 }


  //printf("Hks nc read OK!\n");

}




void free_HSR(){
  int spin, i, iR;

  for (iR=0; iR<TCpyCell+1; iR++){
    for (spin=0; spin<=SpinP_switch; spin++){
      for (i=0; i<T_NumOrbs; i++){
        free(HR[iR][spin][i]);
      }
      free(HR[iR][spin]);
    }
    free(HR[iR]);
  }
  free(HR);



  if(SpinP_switch==3) {
  for (iR=0; iR<TCpyCell+1; iR++){
    for (spin=0; spin<=2; spin++){
      for (i=0; i<T_NumOrbs; i++){
        free(HR_imag[iR][spin][i]);
      }
      free(HR_imag[iR][spin]);
    }
    free(HR_imag[iR]);
  }
  free(HR_imag);
  }



  for (iR=0; iR<TCpyCell+1; iR++){
    for (i=0; i<T_NumOrbs; i++){
      free(SR[iR][i]);
    }
    free(SR[iR]);
  }
  free(SR);

  free(MP);
}


""")

ffi.cdef("""
        extern int atomnum;
        extern int Catomnum;
        extern int Ratomnum;
        extern int SpinP_switch;
        extern int TCpyCell;
        extern int Solver;
        extern double ChemP;
        extern int Valence_Electrons;
        extern double Total_SpinS;
        extern double E_Temp;
        extern int *Total_NumOrbs;
        extern int *FNAN;
        extern int **natn;
        extern int **ncn;
        extern double **atv;
        extern int **atv_ijk;
        extern double tv[4][4];
        extern double rtv[4][4];
        extern double **Gxyz;
        extern double *****Hks;
        extern double *****iHks;
        extern double ****OLP;
        extern double *****DM;
        extern double dipole_moment_core[4];
        extern int version;
        extern int order_max;
        extern double *cc_vec;
        extern double ******OLPpo;
        extern double *****OLPmo;
        extern double dipole_moment_background[4];
        extern int T_NumOrbs; /* total number of orbs*/
        extern int *MP ;         /* Full system version of MP */
        extern double ****HR;
        extern double ****HR_imag;
        extern double ***SR;

        void read_scfout(char *argv[]);
        void free_scfout();
        void prepare_HSR();
        void free_HSR();

        """)
if __name__ == '__main__':
    ffi.compile(verbose=True)
