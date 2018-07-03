/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

// Syntax: widomCC_<machine> P in.lammps
//         P = # of procs to run LAMMPS on
//             must be <= # of procs the driver code itself runs on
//         in.lammps = LAMMPS input script
// See README for compilation instructions

#include <string>
#include <iostream>
#include <fstream>
#include <limits>

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "mpi.h"

#include "lammps.h"         // these are LAMMPS include files
#include "domain.h"
#include "input.h"
#include "atom.h"
#include "force.h"
#include "pair.h"
#include "pair_lj_cut.h"
#include "library.h"
#include "math_const.h"

#include "geometry.h"
#include "output_widom.h"

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
};

using namespace LAMMPS_NS;

int main(int narg, char **arg)
{
  // setup MPI and various communicators
  // driver runs on all procs in MPI_COMM_WORLD
  // comm_lammps only has 1st P procs (could be all or any subset)

  MPI_Init(&narg,&arg);

  if (narg != 3) {
    printf("Syntax: mpirun -n P widomCC P in.lammps\n");
    exit(1);
  }

  int me,nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD,&me);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

  int nprocs_lammps = atoi(arg[1]);
  if (nprocs_lammps > nprocs) {
    if (me == 0)
      printf("ERROR: LAMMPS cannot use more procs than available\n");
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  int lammps;
  if (me < nprocs_lammps) lammps = 1;
  else lammps = MPI_UNDEFINED;
  MPI_Comm comm_lammps;
  MPI_Comm_split(MPI_COMM_WORLD,lammps,0,&comm_lammps);
  
  // open LAMMPS input script

  FILE *fp;
  if (me == 0) {
    fp = fopen(arg[2],"r");
    if (fp == NULL) {
      printf("ERROR: Could not open LAMMPS input script\n");
      MPI_Abort(MPI_COMM_WORLD,1);
    }
  }

  // run the input script thru LAMMPS one line at a time until end-of-file
  // driver proc 0 reads a line, Bcasts it to all procs
  // (could just send it to proc 0 of comm_lammps and let it Bcast)
  // all LAMMPS procs call input->one() on the line
  
  LAMMPS *lmp;
  if (lammps == 1) lmp = new LAMMPS(0,NULL,comm_lammps);

  int n;
  char line[1024];
  while (1) {
    if (me == 0) {
      if (fgets(line,1024,fp) == NULL) n = 0;
      else n = strlen(line) + 1;
      if (n == 0) fclose(fp);
    }
    MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
    if (n == 0) break;
    MPI_Bcast(line,n,MPI_CHAR,0,MPI_COMM_WORLD);
    if (lammps == 1) lmp->input->one(line);
  }

  // run 10 more steps
  // get coords from LAMMPS
  // change coords of 1st atom
  // put coords back into LAMMPS
  // run a single step with changed coords

  if (lammps == 1) {

    lmp->input->one("run 10");

    double *bbox = new double[3];

    double xlo = *( (double*) lammps_extract_global(lmp,"boxxlo") );
    double ylo = *( (double*) lammps_extract_global(lmp,"boxylo") );
    double zlo = *( (double*) lammps_extract_global(lmp,"boxzlo") );
    double xhi = *( (double*) lammps_extract_global(lmp,"boxxhi") );
    double yhi = *( (double*) lammps_extract_global(lmp,"boxyhi") );
    double zhi = *( (double*) lammps_extract_global(lmp,"boxzhi") );

    bbox[0] = xhi - xlo;
    bbox[1] = yhi - ylo;
    bbox[2] = zhi - zlo;

    if (me == 0) {
      std::cout << "x: " << xlo << " " << xhi << std::endl;
      std::cout << "y: " << ylo << " " << yhi << std::endl;
      std::cout << "z: " << zlo << " " << zhi << std::endl;
    }

    double x, y, z;
    x = fRand( xlo, xhi );
    y = fRand( ylo, yhi );
    z = fRand( zlo, zhi );

    // Add one atom
    char add_command[216];
    sprintf(add_command, "create_atoms 1 single %f %f %f", x, y, z);
    lmp->input->one(add_command);
    int natoms = static_cast<int> (lmp->atom->natoms);
    int test_particle = natoms - 1;

    
    // Use read-in xyz for atoms for 1st step
    int *type = new int[natoms];
    double *r = new double[3*natoms];


    std::ifstream file;
    file.open( "dump.all.lammpstrj" );
    if( !file.is_open() ) {
      printf("Unable to open file.\n");
      exit(1);
    }

    bool xyz = false;
    bool lammps = true;
    if( xyz && lammps ) {
      std::cout << "ERROR: Pick either xyz or lammps, not both." << std::endl; 
      exit(1);
    }
    
    std::string element;
    int id, atype;
    double rsq, f, eij;
    double wtest, wtest_sq;
    double expE_kbT;
    double beta_mu, stdev_beta_mu;

    double T = 1.0;
    double kb; 
    kb = static_cast<double> (lmp->force->boltz);
    double beta = 1.0 / ( kb * T );

    int nsamples = 0;

    std::ofstream outfile( "out.file" );
    OutputWidom output;
    if (me == 0) {
      output.header(outfile);
      output.header(std::cout);
    }

    // Gather types
    type = ( (int*) lammps_extract_atom(lmp,"type") );

    // Add tail corrections to pair style (Lennard-Jones only)
    bool tail_flag = true;
    double *etail_ij;
    double V, rho, etail;
    double eps = 1.0;
    double sig = 1.0;
    double rc = 3.5;
    double irc = 1.0 / rc;
    double irc3 = irc*irc*irc;
    double irc6 = irc3*irc3;
    double irc9 = irc6*irc3;

    etail_ij = new double[natoms];     
    if( tail_flag ) {


      for( int i = 0; i < natoms; i++ ) {

        // Get the number of each particle type
        double count[2];
        count[0] = count[1] = 0.0;
        for (int k = 0; k < natoms; k++) {
          if (type[k] == type[i])        count[0] += 1.0; //# of particle type i
          if (type[k] == type[test_particle]) count[1] += 1.0; //# of inserted particle type
        }
       
        double tempcut; 
        tempcut = lmp->force->pair->init_one( type[i], type[test_particle] );
        etail_ij[i]  = lmp->force->pair->etail_ij;
        etail_ij[i] /= count[1];
        V = (lmp->domain->xprd * lmp->domain->yprd * lmp->domain->zprd);
        rho = (natoms - 1) / V; //Subtract test particle in rho
        etail_ij[i] /= V;
        etail = (8.0/3.0)*M_PI*rho*eps*pow(sig, 3.0);
        etail *= ( (1.0/3.0)*irc9 - irc3 );

      }
    }

    // initialize random seed
    srand(3987);

    // Loop over snapshots of trajectory file
    while( !file.eof() ) {

      // gather atom coordinates for this timestep
      lammps_gather_atoms(lmp,"x",1,3,r);

      // XYZ file
      if( xyz ) {

        // skip 1st 2 lines of XYZ
        for( int i = 0; i < 2; i++ )
          file.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );

        // read in xyz 
        for( int i = 0; i < (natoms-1); i++ )  
          file >> element >> r[3*i + 0] >> r[3*i + 1] >> r[3*i + 2];
      
        file.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
      }

      if( lammps ) {
        // skip 1st 9 lines of XYZ
        for( int i = 0; i < 9; i++ )
          file.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );

        // read in xyz 
        for( int i = 0; i < (natoms-1); i++ )  
          file >> id >> atype >> r[3*i + 0] >> r[3*i + 1] >> r[3*i + 2];
      
        file.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
      }

      // initialize energy sum for insertion
      int ninserts = 10000;
      wtest    = 0.0;
      wtest_sq = 0.0;

      // Peform test particle insertion loop
      double cut   = lmp->force->pair->cutforce;
      double cutsq = cut * cut;
      for( int insert = 0; insert < ninserts; insert++ ) {

        // use random coordinates for insertion atom 
        r[3*(test_particle) + 0] = fRand( xlo, xhi );
        r[3*(test_particle) + 1] = fRand( ylo, yhi );
        r[3*(test_particle) + 2] = fRand( zlo, zhi );

        // scatter atoms so we can compute energy w/ new coordinates       
        lammps_scatter_atoms(lmp,"x",1,3,r);

        // Compute PE for insertion atom
        double energy = 0.0;
        for( int i = 0; i < natoms-1; i++ ) {
          rsq     = getDistanceSq( i, test_particle, r, bbox );
        
          if( rsq <= cutsq ) {
            eij = lmp->force->pair->single( i, test_particle, type[i], type[test_particle], rsq, 1.0, 1.0, f);
            energy += eij;
          }
        }

        // collect energy sum for insertion
        expE_kbT  = exp( -1.0 * beta * energy);
        wtest    += expE_kbT;
        wtest_sq += expE_kbT * expE_kbT;

        // collect atoms back up to re-insert last atom
        lammps_gather_atoms(lmp,"x",1,3,r);
      }

      // collect chemical potential & stats
      beta_mu = -1.0 * log( wtest / (double) ninserts );
      stdev_beta_mu = ( wtest_sq - wtest )  / (double) ninserts;
      stdev_beta_mu = sqrt( stdev_beta_mu ) / wtest;

      if( tail_flag )
        beta_mu += 2.0*etail;
        //beta_mu += 2.0*etail_ij[test_particle];

      // scatter atoms back out so we can gather them back up in the loop
      lammps_scatter_atoms(lmp,"x",1,3,r);

      nsamples++;

      if (me == 0) {
        output.add_line( outfile, nsamples, beta, beta_mu, stdev_beta_mu );
        output.add_line( std::cout, nsamples, beta, beta_mu, stdev_beta_mu );
      }
    }


    delete [] r;
    delete [] bbox;
    if( tail_flag )
      delete [] etail_ij;
    file.close();

  }

  if (lammps == 1) delete lmp;

  // close down MPI

  MPI_Finalize();
}

