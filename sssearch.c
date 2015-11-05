/*************************************************************************

   Program:    sssearch
   File:       ssearch.c
   
   Version:    V1.1
   Date:       04.07.96
   Function:   Search for potential disulphide bonding sites
   
   Copyright:  (c) SciTech Software 1993-6
   Author:     Dr. Andrew C. R. Martin
   Address:    SciTech Software
               23, Stag Leys,
               Ashtead,
               Surrey,
               KT21 2TD.
   Phone:      +44 (0)1372 275775
   EMail:      andrew@stagleys.demon.co.uk
               
**************************************************************************

   This program is not in the public domain, but it may be freely copied
   and distributed for no charge providing this header is included.
   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If someone
   else breaks this code, I don't want to be blamed for code that does not
   work! The code may not be sold commercially without prior permission 
   from the author, although it may be given away free with commercial 
   products, providing it is made clear that this program is free and that 
   the source code is provided with the program.

**************************************************************************

   Description:
   ============
   Performs a simple search for potential disulphide bond introduction
   sites by looking at the distances between CBs.

   Obviously, this rules out any residues which are currently glycines,
   but this is probably a good idea anyway.

   Optionally sites involving prolines may also be eliminated.

**************************************************************************

   Usage:
   ======

**************************************************************************

   Notes:
   ======
   
**************************************************************************

   Revision History:
   =================
   V1.0  12.10.93 Original
   V1.1  04.07.96 Moved ; one to right in output and added stdlib.h

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "bioplib/SysDefs.h"
#include "bioplib/MathType.h"
#include "bioplib/pdb.h"
#include "bioplib/macros.h"
#include "bioplib/parse.h"

/************************************************************************/
/* Defines
*/
#define MAXSTRLEN   160
#define MAXSTRPARAM 1
#define MAXNUMPARAM 2

#define KEY_CBDIST  0
#define KEY_CADIST  1
#define KEY_GO      2
#define KEY_QUIT    3
#define KEY_HELP    4
#define KEY_DOPRO   5
#define KEY_NOPRO   6
#define NCOMM       7

/* Default parameters for search ranges
*/
#define CAMIN       3.9
#define CAMAX       8.3
#define CBMIN       2.8
#define CBMAX       4.6

/************************************************************************/
/* Globals
*/
KeyWd gKeyWords[NCOMM];
char  *gStrParam[MAXSTRPARAM];
REAL  gNumParam[MAXNUMPARAM],
      gCADist[2],
      gCBDist[2];
BOOL  gDoProlines = FALSE,
      gRanking    = TRUE;

/*************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL CheckCmdLine(int argc, char **argv, char *filename);
BOOL SetupParser(void);
PDB *OpenAndReadPDB(char *file, FILE *Msgfp);
BOOL DoParseLoop(PDB *pdb);
void InitDefaults(void);
void ShowHelp(void);
void Usage(BOOL ShowUse);
void DoCalculations(PDB *pdb);
void DoRanking(void);

/*************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for disulphide searching

   12.10.93 Original   By: ACRM
*/
int main(int argc, char **argv)
{
   char filename[160];
   PDB  *pdb = NULL;

   /* Initialise defaults for parameters                                 */
   InitDefaults();

   /* Check command line                                                 */
   if(CheckCmdLine(argc, argv, filename))
   {
      if((pdb = OpenAndReadPDB(filename,stderr)) != NULL)
      {
         /* Initialise command parser                                    */
         if(SetupParser())
         {
            if(DoParseLoop(pdb) && gRanking)
               DoRanking();
         }
      }
   }
}

/*************************************************************************/
/*>BOOL CheckCmdLine(int argc, char **argv, char *filename)
   --------------------------------------------------------
   Check the command line outputting the filename and returning a flag
   to say if all was OK. Calls Usage() if not.

   12.10.93 Original   By: ACRM
*/
BOOL CheckCmdLine(int argc, char **argv, char *filename)
{
   argc--; argv++;

   while(argc > 1)
   {
      if(argv[0][0] == '-')
      {
         /* Process command line switches                                */
      }
      else
      {
         fprintf(stderr,"Error in command line at %s\n",argv[0]);
         return(FALSE);
      }
      argc--; argv++;
   }

   if(argc != 1) 
   {
      Usage(TRUE);
      return(FALSE);
   }
   strcpy(filename,argv[0]);
   return(TRUE);
}

/************************************************************************/
/*>BOOL SetupParser(void)
   ----------------------
   Set up the command parser

   10.10.93 Original    By: ACRM
*/
BOOL SetupParser(void)
{
   int i;

   /* Initialise returned string array                                  */
   for(i=0; i<MAXSTRPARAM; i++)
   {
      if((gStrParam[i] = (char *)malloc(MAXSTRLEN * sizeof(char))) == NULL)
      {
         return(FALSE);
      }
   }

   /* Create keywords                                                   */
   MAKEKEY(gKeyWords[KEY_CBDIST], "CBDIST", NUMBER, 2);
   MAKEKEY(gKeyWords[KEY_CADIST], "CADIST", NUMBER, 2);
   MAKEKEY(gKeyWords[KEY_GO],     "GO",     STRING, 0);
   MAKEKEY(gKeyWords[KEY_QUIT],   "QUIT",   STRING, 0);
   MAKEKEY(gKeyWords[KEY_HELP],   "HELP",   STRING, 0);
   MAKEKEY(gKeyWords[KEY_DOPRO],  "DOPRO",  STRING, 0);
   MAKEKEY(gKeyWords[KEY_NOPRO],  "NOPRO",  STRING, 0);

   /* Check on MAKEKEYs                                                 */
   for(i=0; i<NCOMM; i++)
      if(gKeyWords[i].name == NULL) return(FALSE);

   return(TRUE);
}

/************************************************************************/
/*>PDB *OpenAndReadPDB(char *file, FILE *Msgfp)
   --------------------------------------------
   Opens a PDB file specified by name and reads it with ReadPDB. Error
   messages are sent to the specified file. If NULL, no messages will
   be issued. Returns a pointer to the PDB linked list or NULL if failed.

   10.10.93 Original    By: ACRM
   11.10.93 Added legal NULL Msgfp
*/
PDB *OpenAndReadPDB(char *file, FILE *Msgfp)
{
   FILE *fp = NULL;
   PDB  *pdb = NULL;
   int  natoms;

   if((fp=fopen(file,"r"))==NULL)
   {
      if(Msgfp!=NULL) 
         fprintf(Msgfp,"Unable to open file: %s\n",file);
   }
   else
   {
      pdb = ReadPDB(fp, &natoms);
      if(pdb == NULL || natoms == 0)
      {
         if(Msgfp!=NULL) 
            fprintf(Msgfp,"No atoms read from file: %s\n",file);
         if(pdb!=NULL) FREELIST(pdb,PDB);
         pdb = NULL;
      }
      fclose(fp);
   }
   return(pdb);
}

/************************************************************************/
/*>void DoParseLoop(PDB *pdb)
   --------------------------
   The main command parser loop.

   10.10.93 Original    By: ACRM
*/
BOOL DoParseLoop(PDB *pdb)
{
   char buffer[160];
   int key;

   fprintf(stderr,"SSSearch> ");

   while(fgets(buffer,160,stdin))
   {
      TERMINATE(buffer);
      key = parse(buffer,NCOMM,gKeyWords,gNumParam,gStrParam);
      switch(key)
      {
      case PARSE_ERRC:
         fprintf(stderr,"Error in command: %s\n",buffer);
         break;
      case PARSE_ERRP:
         fprintf(stderr,"Error in parameters: %s\n",buffer);
         break;
      case KEY_CBDIST:
         gCBDist[0] = gNumParam[0];
         gCBDist[1] = gNumParam[1];
         break;
      case KEY_CADIST:
         gCADist[0] = gNumParam[0];
         gCADist[1] = gNumParam[1];
         break;
      case KEY_HELP:
         ShowHelp();
         break;
      case KEY_DOPRO:
         gDoProlines = TRUE;
         break;
      case KEY_NOPRO:
         gDoProlines = FALSE;
         break;
      case KEY_GO:
         DoCalculations(pdb);
         return(TRUE);
      case KEY_QUIT:
         return(FALSE);
      }
      fprintf(stderr,"SSSearch> ");
   }
}

/*************************************************************************/
/*>void InitDefaults(void)
   -----------------------
   Initialise default values of variables

   12.10.93 Original   By: ACRM
*/
void InitDefaults(void)
{
   gCADist[0]  = CAMIN;
   gCADist[1]  = CAMAX;
   gCBDist[0]  = CBMIN;
   gCBDist[1]  = CBMAX;
   gDoProlines = FALSE;
}

/*************************************************************************/
/*>void ShowHelp(void)
   -------------------
   Display help message

   12.10.93 Original   By: ACRM
*/
void ShowHelp(void)
{
   Usage(FALSE);

   fprintf(stderr,"CADIST <min> <max>   Specify the CA distance range\n");
   fprintf(stderr,"CBDIST <min> <max>   Specify the CB distance range\n");
   fprintf(stderr,"DOPRO                Include prolines in the \
results\n");
   fprintf(stderr,"NOPRO                Exclude prolines from the \
results\n");
   fprintf(stderr,"GO                   Start calculations and end \
program\n");
   fprintf(stderr,"QUIT                 End program without \
calculations\n\n");

   fprintf(stderr,"Current parameters are:\n");
   fprintf(stderr,"C-alpha distance range: %lf %lf\n",
           gCADist[0],gCADist[1]);
   fprintf(stderr,"C-beta  distance range: %lf %lf\n",
           gCBDist[0],gCBDist[1]);
   fprintf(stderr,"Prolines will%sbe included in the results\n\n",
          ((gDoProlines)?" ":" NOT ")); 
}

/*************************************************************************/
/*>void Usage(BOOL ShowUse)
   ------------------------
   Display usage message or start of help message

   12.10.93 Original   By: ACRM
*/
void Usage(BOOL ShowUse)
{
   fprintf(stderr,"SSSearch V1.1 (c) 1993-6, Dr. A.C.R. Martin, SciTech \
Software, DKfz\n");
   fprintf(stderr,"This program is freely distributable providing no \
profit is made in so doing\n");
   fprintf(stderr,"SSSearch searches for potential disulphide bond sites \
by looking for\n");
   fprintf(stderr,"residues with CB distances and CA distances in \
specified ranges.\n\n");
   
   if(ShowUse)
   {
      fprintf(stderr,"Usage: sssearch [>search.out] <file.pdb>\n\n");
   }
}

/*************************************************************************/
/*>void DoCalculations(PDB *pdb)
   -----------------------------
   Do the actual calculations on a PDB linked list

   12.10.93 Original   By: ACRM
   04.07.96 Moved ; one space to right in output
*/
void DoCalculations(PDB *pdb)
{
   PDB  *start1 = NULL,
        *end1   = NULL,
        *start2 = NULL,
        *end2   = NULL,
        *p      = NULL,
        *KeyCA  = NULL,
        *KeyCB  = NULL,
        *CA     = NULL,
        *CB     = NULL;
   REAL DistCA,
        DistCB;
   FILE *tempfp = NULL;

   if((tempfp=fopen("sssearch.tmp","w"))==NULL)
   {
      tempfp = stdout;
      fprintf(stderr,"Unable to open temp file; ranking not possible\n");
      gRanking = FALSE;
   }

   /* Step through PDB linked list a residue at a time                   */
   for(start1=pdb; start1!=NULL; start1 = end1)
   {
      /* Find next residue                                               */
      end1 = FindEndPDB(start1);

      /* Step through this residue looking for CA and CB                 */
      KeyCA = KeyCB = NULL;
      for(p=start1; p!=end1; NEXT(p))
      {
         if(!strncmp(p->atnam,"CA  ",4)) KeyCA = p;
         if(!strncmp(p->atnam,"CB  ",4)) KeyCB = p;
      }

      /* If we found both, then search all following residues            */
      if(KeyCA != NULL && KeyCB != NULL)
      {
         for(start2=end1; start2!=NULL; start2 = end2)
         {
            /* Find next residue                                         */
            end2 = FindEndPDB(start2);

            /* Step through this residue looking for CA and CB           */
            CA = CB = NULL;
            for(p=start2; p!=end2; NEXT(p))
            {
               if(!strncmp(p->atnam,"CA  ",4)) CA = p;
               if(!strncmp(p->atnam,"CB  ",4)) CB = p;
            }

            if(CA != NULL && CB != NULL)
            {
               /* If not PRO, or gDoProlines flag set, calc dists        */
               if(gDoProlines || 
                  (strncmp(KeyCA->resnam,"PRO ",4) && 
                   strncmp(CA->resnam,"PRO ",4)))
               {
                  DistCA = DIST(KeyCA,CA);
                  DistCB = DIST(KeyCB,CB);

                  /* If in ranges print message                          */
                  if(DistCA >= gCADist[0] &&
                     DistCA <= gCADist[1] &&
                     DistCB >= gCBDist[0] &&
                     DistCB <= gCBDist[1])
                  {
                     fprintf(tempfp,"%s %c%5d%c with %s %c%5d%c ; \
CA %5.2lf, CB %5.2lf\n",
KeyCA->resnam, KeyCA->chain[0], KeyCA->resnum, KeyCA->insert[0],
   CA->resnam,    CA->chain[0],    CA->resnum,    CA->insert[0],
   DistCA, DistCB);
                  }
               }
            }
         }
      }
   }
   if(gRanking) fclose(tempfp);
}

/*************************************************************************/
/*>void DoRanking(void)
   --------------------
   Sorts the output of the calculations using the UNIX sort utility which
   outputs to stdout. Then removes the temporary file which was used to
   store the unsorted results

   12.10.93 Original   By: ACRM
*/
void DoRanking(void)
{
   system("sort +11n +9n sssearch.tmp");
   remove("sssearch.tmp");
}

