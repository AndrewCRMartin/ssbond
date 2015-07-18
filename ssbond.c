/*************************************************************************

   Program:    ssbond
   File:       ssbond.c
   
   Version:    V2.0
   Date:       08.07.96
   Function:   Do conf search for potential disulphides
   
   Copyright:  (c) Dr. Andrew C. R. Martin 1995
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   EMail:      andrew@bioinf.org.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If someone
   else breaks this code, I don't want to be blamed for code that does not
   work! 

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0  25.10.88 Original FORTRAN version
   V1.1  08.10.93 Various fixes while at DKfz
   V2.0  08.07.96 Rewritten in C to support chain names and inserts

*************************************************************************/
/* Includes
*/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "bioplib/pdb.h"
#include "bioplib/general.h"
#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/macros.h"
#include "bioplib/angle.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF    160

#define DEF_BOND   ((REAL)1.813) /* CB-SG bond length                   */
#define DEF_ANGLE  ((REAL)114.0) /* CA-CB-SG angle (Sigma = 20degrees   */
#define DEF_DIDEAL ((REAL)2.03)  /* SG-SG bond length                   */
#define DEF_STEP   ((REAL)5.0)   /* Rotation step (degrees)             */
#define DEF_NBEST  10            /* Number of confs to keep             */

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  REAL *bond, REAL *angle, REAL *dideal,
                  REAL *step, int *nbest, char *rspec1, char *rspec2);
void Usage(void);
BOOL DoSearch(FILE *out, PDB *pdb, REAL bond, REAL angle, REAL dIdeal, 
              REAL step, int nbest, char *rspec1, char *rspec2);
BOOL GetResSpecs(char *rspec1, char *rspec2);
BOOL FindCoords(PDB *pdb, char *rspec, VEC3F *N, VEC3F *CA, VEC3F *CB);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for disulphide conformational search

   08.07.96 Original   By: ACRM
   09.07.96 Now loops if residues not on command line
*/
int main(int argc, char **argv)
{
   char infile[MAXBUFF],
        outfile[MAXBUFF],
        rspec1[16],
        rspec2[16];
   FILE *in    = stdin,
        *out   = stdout;
   REAL bond,                               /* CB-SG                    */
        angle,                              /* CA-CB-SG                 */
        dideal,                             /* SG1-SG2                  */
        step;
   int  nbest,          
        natoms;
   PDB  *pdb;
   BOOL DoLoop = TRUE;
   
   if(ParseCmdLine(argc, argv, infile, outfile, &bond, &angle, &dideal,
                   &step, &nbest, rspec1, rspec2))
   {
      /* If residue specs given on command line, we don't loop          */
      if(rspec1[0] != '\0' && rspec2[0] != '\0')
         DoLoop = FALSE;

      if(OpenStdFiles(infile, outfile, &in, &out))
      {
         if((pdb=ReadPDB(in, &natoms))!=NULL)
         {
            if(DoLoop)
            {
               while(GetResSpecs(rspec1, rspec2))
               {
                  DoSearch(out, pdb, bond, angle, dideal, step, nbest, 
                           rspec1, rspec2);
               }
            }
            else
            {
               DoSearch(out, pdb, bond, angle, dideal, step, nbest, 
                        rspec1, rspec2);
            }
         }
         else
         {
            fprintf(stderr,"Unable to read input PDB file\n");
            return(1);
         }
      }
   }
   else
   {
      Usage();
   }
   return(0);
}


/************************************************************************/
/*>BOOL DoSearch(FILE *out, PDB *pdb, REAL bond, REAL angle, REAL dIdeal, 
                 REAL step, int nbest, char *rspec1, char *rspec2)
   -----------------------------------------------------------------------
   Does the real work of conformational search on a residue pair.
   Allocates temporary arrays, finds the coordinates, does the search
   sorts and displays the results and frees the temporary arrays.

   09.07.96 Original   By: ACRM
*/
BOOL DoSearch(FILE *out, PDB *pdb, REAL bond, REAL angle, REAL dIdeal, 
              REAL step, int nbest, char *rspec1, char *rspec2)
{
   VEC3F N1,   N2,
         CA1,  CA2,
         CB1,  CB2,
         SG1,  SG2;
   REAL  *diffArray,
         *ang1Array,
         *ang2Array,
         ang1, ang2, d2,
         chi3, dist,
         dIdealSq  = dIdeal * dIdeal;
   int   ArraySize,
         *Indx,
         Count,
         i;
   BOOL  RetVal = TRUE;
   
   /* Calculate the array size. This is the number of steps squared     */
   ArraySize = 1+(int)(2.0*PI/step);
   ArraySize *= ArraySize;

   /* Allocate storage for the differences from ideal, the torsions used
      and the index array
   */
   if((diffArray=(REAL *)malloc(ArraySize*sizeof(REAL)))==NULL)
   {
      return(FALSE);
   }
   if((ang1Array=(REAL *)malloc(ArraySize*sizeof(REAL)))==NULL)
   {
      RetVal = FALSE;
      goto Cleanup;
   }
   if((ang2Array=(REAL *)malloc(ArraySize*sizeof(REAL)))==NULL)
   {
      RetVal = FALSE;
      goto Cleanup;
   }
   if((Indx=(int *)malloc(ArraySize*sizeof(int)))==NULL)
   {
      RetVal = FALSE;
      goto Cleanup;
   }
   
   /* Find the coordinates for the atoms in the two residues            */
   if(!FindCoords(pdb, rspec1, &N1, &CA1, &CB1))
   {
      fprintf(stderr,"Unable to find N, CA and CB for %s\n",rspec1);
      RetVal = FALSE;
      goto Cleanup;
   }
   if(!FindCoords(pdb, rspec2, &N2, &CA2, &CB2))
   {
      fprintf(stderr,"Unable to find N, CA and CB for %s\n",rspec2);
      RetVal = FALSE;
      goto Cleanup;
   }

   /* Do the conformational search, storing the differences and angles  */
   Count = 0;
   for(ang1=0.0; ang1<(2.0*PI-(step/2.0)); ang1+=step)
   {
      TorToCoor(N1, CA1, CB1, bond, angle, ang1, &SG1);
      
      for(ang2=0.0; ang2<(2.0*PI-(step/2.0)); ang2+=step)
      {
         TorToCoor(N2, CA2, CB2, bond, angle, ang2, &SG2);
         d2 = DISTSQ((&SG1),(&SG2));
         ang1Array[Count] = ang1;
         ang2Array[Count] = ang2;
         diffArray[Count] = ABS((d2-dIdealSq));
         Count++;
      }
   }
   
   /* Index sort the diffArray                                          */
   IndexReal(diffArray, Indx, Count);
   
   /* Print the best orientations                                       */
   fprintf(out,"\nBest orientations are (CA-CB-SG = %.2f):\n",
           (REAL)180.0*angle/PI);
   fprintf(out,"Optimum bond length:  %f\n",dIdeal);
   fprintf(out,"Residue %5s                      Residue %5s\n",
           rspec1, rspec2);
   fprintf(out,"Torsion           SG Coords        Torsion           \
SG Coords\n");
   fprintf(out,"N-CA-CB-SG    X       Y       Z    N-CA-CB-SG    \
X       Y       Z       DIST    CHI3\n");
   
   for(i=0; i<nbest; i++)
   {
      TorToCoor(N1, CA1, CB1, bond, angle, ang1Array[Indx[i]], &SG1);
      TorToCoor(N2, CA2, CB2, bond, angle, ang2Array[Indx[i]], &SG2);
   
      chi3=phi(CB1.x, CB1.y, CB1.z,
               SG1.x, SG1.y, SG1.z, 
               SG2.x, SG2.y, SG2.z,
               CB2.x, CB2.y, CB2.z);

      dist=DIST((&SG1), (&SG2));
      
      fprintf(out,"%6.1f    %8.3f%8.3f%8.3f %6.1f    %8.3f%8.3f%8.3f\
%8.3f   %6.1f\n",
              (REAL)180.0*ang1Array[Indx[i]]/PI,
              SG1.x, SG1.y, SG1.z, 
              (REAL)180.0*ang2Array[Indx[i]]/PI,
              SG2.x, SG2.y, SG2.z, 
              dist, 
              (REAL)180.0*chi3/PI);
   }

Cleanup:
   if(Indx!=NULL)      free(Indx);
   if(ang2Array!=NULL) free(ang2Array);
   if(ang1Array!=NULL) free(ang1Array);
   if(diffArray!=NULL) free(diffArray);

   return(TRUE);
}

/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                     REAL *bond, REAL *angle, REAL *dideal,
                     REAL *step, int *nbest, char *rspec1, char *rspec2)
   ---------------------------------------------------------------------
   Parse the command line

   09.07.96 Original   By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  REAL *bond, REAL *angle, REAL *dideal,
                  REAL *step, int *nbest, char *rspec1, char *rspec2)
{
   argc--;
   argv++;

   infile[0] = outfile[0] = '\0';
   *bond     = DEF_BOND;
   *angle    = DEF_ANGLE * PI / 180.0;
   *dideal   = DEF_DIDEAL;
   *step     = DEF_STEP * PI / 180.0;
   *nbest    = DEF_NBEST;
   rspec1[0] = '\0';
   rspec2[0] = '\0';
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'b':
            argc--;
            argv++;
            if(argc<0) return(FALSE);
            if(!sscanf(argv[0],"%lf",bond))
               return(FALSE);
            break;
         case 'a':
            argc--;
            argv++;
            if(argc<0) return(FALSE);
            if(!sscanf(argv[0],"%lf",angle))
               return(FALSE);
            *angle *= PI/180.0;
            break;
         case 'd':
            argc--;
            argv++;
            if(argc<0) return(FALSE);
            if(!sscanf(argv[0],"%lf",dideal))
               return(FALSE);
            break;
         case 's':
            argc--;
            argv++;
            if(argc<0) return(FALSE);
            if(!sscanf(argv[0],"%lf",step))
               return(FALSE);
            *step *= PI/180.0;
            break;
         case 'n':
            argc--;
            argv++;
            if(argc<0) return(FALSE);
            if(!sscanf(argv[0],"%d",nbest))
               return(FALSE);
            break;
         case 'r':
            argc--;
            argv++;
            if(argc<0) return(FALSE);
            strcpy(rspec1, argv[0]);
            argc--;
            argv++;
            if(argc<0) return(FALSE);
            strcpy(rspec2, argv[0]);
            break;
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are only 1 or 2 arguments left             */
         if(argc > 2)
            return(FALSE);
         
         /* Copy the first to infile                                    */
         strcpy(infile, argv[0]);
         
         /* If there's another, copy it to outfile                      */
         argc--;
         argv++;
         if(argc)
            strcpy(outfile, argv[0]);
            
         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   return(TRUE);
}


/************************************************************************/
/*>void Usage(void)
   ----------------
   Prints a usage message

   09.07.96 Original   By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"\nSSBond V2.0 (c) 1989-1996 Dr. Andrew C.R. Martin \
LMB Oxford, DKfz, UCL\n");

   fprintf(stderr,"\nUsage: ssbond [-b bondlen] [-a angle] [-d dideal] \
[-s step] [-n nbest] \n");
   fprintf(stderr,"              [-r rspec1 rspec2] [in.pdb \
[out.txt]]\n");
   fprintf(stderr,"       -b Specify the CB-SG bond length [%.2f]\n", 
           DEF_BOND);
   fprintf(stderr,"       -a Specify the CA-CB-SG angle [%.2f]\n",
           DEF_ANGLE);
   fprintf(stderr,"       -d Specify the ideal SG-SG bond length \
[%.2f]\n", DEF_DIDEAL);
   fprintf(stderr,"       -s Specify the rotation step (degrees) \
[%.2f]\n", DEF_STEP);
   fprintf(stderr,"       -n Specify the number of conformations to \
keep [%d]\n", DEF_NBEST);
   fprintf(stderr,"       -r Specify the two residues to search\n");

   fprintf(stderr,"\nIf input and output files are not specified \
stdin/stdout will be used.\n");

   fprintf(stderr,"\nSSBond substitutes Cys residues at two sites and \
performs a conformational \n");
   fprintf(stderr,"search by spinning the CA-CB torsion angles in an \
attempt to optimise the\n");
   fprintf(stderr,"SG-SG separation.\n");

   fprintf(stderr,"\nResidues to be searched are specified as [c]NNN[i] \
where [c] is an optional\n");
   fprintf(stderr,"chain name, NNN is a residue number and [i] is an \
optional insertion code.\n");

   fprintf(stderr,"\nIf the residue pair is not specified on the command \
line, the program will\n");
   fprintf(stderr,"prompt for a residue pair to be entered. This may be \
repeated until QUIT \n");
   fprintf(stderr,"or EXIT is entered to end the program.\n\n");
}


/************************************************************************/
/*>BOOL GetResSpecs(char *rspec1, char *rspec2)
   --------------------------------------------
   Gets residue specifications from stdin. Prints a prompt if stdin
   is a tty

   09.07.96 Original   By: ACRM
*/
BOOL GetResSpecs(char *rspec1, char *rspec2)
{
   char buffer[MAXBUFF];
   
   PROMPT(stdin,"Enter residue specifications as [c]NNN[i] [c]NNN[i]: ");

   if(!fgets(buffer,MAXBUFF,stdin))
      return(FALSE);
   
   if(sscanf(buffer,"%s %s",rspec1,rspec2) != 2)
      return(FALSE);
   if(!upstrncmp(rspec1,"QUIT",4) || !upstrncmp(rspec1,"EXIT",4))
      return(FALSE);
   
   return(TRUE);
}


/************************************************************************/
/*>BOOL FindCoords(PDB *pdb, char *rspec, VEC3F *N, VEC3F *CA, VEC3F *CB)
   ----------------------------------------------------------------------
   Finds the coordinates for the atoms required

   09.07.96 Original   By: ACRM
*/
BOOL FindCoords(PDB *pdb, char *rspec, VEC3F *N, VEC3F *CA, VEC3F *CB)
{
   PDB *start,
       *stop,
       *p;
   
   /* Set NULL coordinates for the 3 atoms                              */
   N->x  = N->y  = N->z  = (REAL)9999.000;
   CA->x = CA->y = CA->z = (REAL)9999.000;
   CB->x = CB->y = CB->z = (REAL)9999.000;
   
   /* Find the boundaries of the specified residue                      */
   if((start = FindResidueSpec(pdb, rspec)) == NULL)
      return(FALSE);
   stop = FindNextResidue(start);
   
   /* Find the atoms in that residue                                    */
   for(p=start; p!=stop; NEXT(p))
   {
      if(!strncmp(p->atnam,"N   ",4))
      {
         N->x = p->x;
         N->y = p->y;
         N->z = p->z;
      }
      else if(!strncmp(p->atnam,"CA  ",4))
      {
         CA->x = p->x;
         CA->y = p->y;
         CA->z = p->z;
      }
      else if(!strncmp(p->atnam,"CB  ",4))
      {
         CB->x = p->x;
         CB->y = p->y;
         CB->z = p->z;
      }
   }

   if((N->x > 9998.0) || (CA->x > 9998.0) || (CB->x > 9998.0))
      return(FALSE);
   
   return(TRUE);
}

