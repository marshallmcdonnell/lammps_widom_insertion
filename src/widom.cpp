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
#include <sstream>
#include <limits>
#include <vector>

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "mpi.h"

#include "lammps.h"         // these are LAMMPS include files
#include "domain.h"
#include "input.h"
#include "compute.h"
#include "modify.h"
#include "neighbor.h"
#include "atom.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "library.h"
#include "math_const.h"
#include "comm.h"
#include "update.h"
#include "read_dump.h"
#include "error.h"

#include "geometry.h"
#include "output_widom.h"

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
};

using namespace LAMMPS_NS;

void init(int narg, char **arg, int& lammps, LAMMPS *&lmp) {
  // setup MPI and various communicators
  // driver runs on all procs in MPI_COMM_WORLD
  // comm_lammps only has 1st P procs (could be all or any subset)

  MPI_Init(&narg,&arg);

  if (narg != 7) {
    printf("Syntax: mpirun -n P widomCC P in.lammps trajectory_file ninserts seed temp\n");
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
}


double energy_full(void *ptr)
{
  LAMMPS *lmp = (LAMMPS *) ptr;

  if (lmp->domain->triclinic) lmp->domain->x2lamda(lmp->atom->nlocal);
  lmp->domain->pbc();
  lmp->comm->exchange();
  lmp->atom->nghost = 0;
  lmp->comm->borders();
  if (lmp->domain->triclinic) lmp->domain->lamda2x(lmp->atom->nlocal+lmp->atom->nghost);
  if (lmp->modify->n_pre_neighbor) lmp->modify->pre_neighbor();
  lmp->neighbor->build(1); // 1 is topology flag and the default throughout lammps
  int eflag = 1;
  int vflag = 0;
  // clear forces so they don't accumulate over multiple
  // calls within fix gcmc timestep, e.g. for fix shake

  size_t nbytes = sizeof(double) * (lmp->atom->nlocal + lmp->atom->nghost);
  if (nbytes) memset(&lmp->atom->f[0][0],0,3*nbytes);

  if (lmp->modify->n_pre_force) lmp->modify->pre_force(vflag);

  if (lmp->force->pair) lmp->force->pair->compute(eflag,vflag);

  if (lmp->atom->molecular) {
    if (lmp->force->bond) lmp->force->bond->compute(eflag,vflag);
    if (lmp->force->angle) lmp->force->angle->compute(eflag,vflag);
    if (lmp->force->dihedral) lmp->force->dihedral->compute(eflag,vflag);
    if (lmp->force->improper) lmp->force->improper->compute(eflag,vflag);
  }

  if (lmp->force->kspace) lmp->force->kspace->compute(eflag,vflag);

  // unlike Verlet, not performing a reverse_comm() or forces here
  // b/c GCMC does not care about forces
  // don't think it will mess up energy due to any post_force() fixes

  if (lmp->modify->n_post_force) lmp->modify->post_force(vflag);
  if (lmp->modify->n_end_of_step) lmp->modify->end_of_step();

  // NOTE: all fixes with THERMO_ENERGY mask set and which
  //   operate at pre_force() or post_force() or end_of_step()
  //   and which user has enable via fix_modify thermo yes,
  //   will contribute to total MC energy via pe->compute_scalar()

  // Setup pe compute (c_pe)
  Compute *c_pe;
  char *id_pe = (char *) "thermo_pe";
  int ipe = lmp->modify->find_compute(id_pe);
  c_pe = lmp->modify->compute[ipe];

  lmp->update->eflag_global = lmp->update->ntimestep;
  double total_energy = c_pe->compute_scalar();

  return total_energy;
}


int main(int narg, char **arg)
{
  // Initialization of lammps
  int lammps;
  LAMMPS *lmp;

  init(narg, arg, lammps, lmp);

  if (lammps != 1) {
    printf("ERROR: LAMMPS failed to initialize\n");
    exit(1);
  }

  int me;
  MPI_Comm_rank(MPI_COMM_WORLD,&me);

  lmp->input->one("run 0");

  double xlo = *( (double*) lammps_extract_global(lmp,"boxxlo") );
  double ylo = *( (double*) lammps_extract_global(lmp,"boxylo") );
  double zlo = *( (double*) lammps_extract_global(lmp,"boxzlo") );
  double xhi = *( (double*) lammps_extract_global(lmp,"boxxhi") );
  double yhi = *( (double*) lammps_extract_global(lmp,"boxyhi") );
  double zhi = *( (double*) lammps_extract_global(lmp,"boxzhi") );

  double coord[3];
  coord[0] = fRand( xlo, xhi );
  coord[1] = fRand( ylo, yhi );
  coord[2] = fRand( zlo, zhi );

  // Get initial number of atoms
  int natoms_old = static_cast<int> (lmp->atom->natoms);

  // Setup add atom command
  char add_command[216];
  sprintf(add_command, "create_atoms 1 single %f %f %f units box", 
          coord[0], coord[1], coord[2] );

  // Setup delete atom command
  char delete_command[216];
  sprintf(delete_command, "delete_atoms group test_particle");

  // Add atom
  lmp->input->one(add_command);
  int natoms_new = static_cast<int> (lmp->atom->natoms);

  // Get original N system group and test particle N+1 system group
  char group_original_command[216];
  sprintf(group_original_command, "group original_system id %d:%d", 1, natoms_old);
  lmp->input->one(group_original_command);
  char group_test_particle_command[216];
  sprintf(group_test_particle_command, "group test_particle id %d:%d", natoms_old+1, natoms_new);
  lmp->input->one(group_test_particle_command);  

  // Delete atom to reset before loop
  lmp->input->one(delete_command); 

  double energy;
  energy = energy_full(lmp);

  // Gather xyz and check gather + scatter
  double *r = new double[3*natoms_old];
  lammps_gather_atoms(lmp,"x",1,3,r);
  lammps_scatter_atoms(lmp,"x",1,3,r);

  // Open trajectory file 
  std::ifstream file;
  file.open( arg[3] );
  if( !file.is_open() ) {
    printf("Unable to open trajectory file.\n");
    exit(1);
  }

  // Initialize variables
  std::string element;
  int id, atype;
  double rsq, f, eij;
  double wtest, wtest_sq, wtest_sum, wtest_sq_sum, var_wtest;
  double expE_kbT;
  double beta_mu, stdev_beta_mu;
  int nsamples = 0;


  // Setup beta term
  double T = atof(arg[6]);
  double kb; 
  double beta = 1.0 / ( lmp->force->boltz * T );

  // Setup output file
  std::ofstream outfile( "out.file" );
  OutputWidom output;
  if (me == 0) {
    output.header(outfile);
    output.header(std::cout);
  }

  // initialize random seed
  srand(atoi(arg[5]));

  // Setup reader
  ReadDump *rd = new ReadDump(lmp);
  int timestep = 0;
  int read_dump_narg = 13;
  char **read_dump_args = new char*[read_dump_narg];
  read_dump_args[0] = (char *) "dump.all.lammpstrj";
  read_dump_args[1] = 0;
  read_dump_args[2] = (char *) "x";
  read_dump_args[3] = (char *) "y";
  read_dump_args[4] = (char *) "z";
  read_dump_args[5] = (char *) "scaled";
  read_dump_args[6] = (char *) "yes";
  read_dump_args[7] = (char *) "box";
  read_dump_args[8] = (char *) "yes";
  read_dump_args[9] = (char *) "replace";
  read_dump_args[10] = (char *) "yes";
  read_dump_args[11] = (char *) "format";
  read_dump_args[12] = (char *) "native";
  rd->store_files(1,&read_dump_args[0]);
  int nremain = read_dump_narg - 2;
  if (nremain)
    nremain = rd->fields_and_keywords(nremain,&read_dump_args[read_dump_narg-nremain]);
  else nremain = rd->fields_and_keywords(0,NULL);
  if (nremain) rd->setup_reader(nremain,&read_dump_args[read_dump_narg-nremain]);
  else rd->setup_reader(0,NULL);

  // Turn off screen
  lmp->screen = NULL;

  // Loop over snapshots of trajectory file
  int nstep = 0;
  while( 1 ) {

    // Read in nextsnapshot
    bigint ntimestep = rd->seek(nstep,0); // 0 => find first step >= nstep
    nstep = ntimestep; // if we didn't find exactly nstep, this sets it to what was found
    if (ntimestep < 0)
      break;
    rd->header(1);
    lmp->update->reset_timestep(nstep);
    rd->atoms();
      
    // Post-read, Pre-energy setup
    if (lmp->force->kspace) lmp->force->kspace->qsum_qsq();
    if (lmp->force->pair->tail_flag) lmp->force->pair->reinit();

    // get initial pe for N system to get diff from N+1 system
    double energy_before;
    energy_before = energy_full(lmp);

    // initialize energy sum for insertion
    int ninserts = atoi(arg[4]);
    wtest    = 0.0;
    wtest_sq = 0.0;
      
    // Peform test particle insertion loop
    double energy_after, delta_energy;
    for( int insert = 0; insert < ninserts; insert++ ) {

      // get random point for COM of test particle
      coord[0] = fRand(xlo, xhi);
      coord[1] = fRand(ylo, yhi);
      coord[2] = fRand(zlo, zhi);

      // Add test particle
      sprintf(add_command, "create_atoms 1 single %f %f %f units box",
              coord[0], coord[1], coord[2] );
      lmp->input->one(add_command);

      // Post-insertion, Pre-energy setup
      if (lmp->domain->triclinic) lmp->domain->x2lamda(lmp->atom->nlocal);
      lmp->comm->borders();
      if (lmp->domain->triclinic) lmp->domain->lamda2x(lmp->atom->nlocal+lmp->atom->nghost);
      if (lmp->force->kspace) lmp->force->kspace->qsum_qsq();
      if (lmp->force->pair->tail_flag) lmp->force->pair->reinit();

      // Get energy diff from N -> N+1 system
      energy_after = energy_full(lmp);
      delta_energy = energy_after - energy_before;

      // collect energy sum for insertion
      expE_kbT  = exp( -1.0 * beta * delta_energy);
      wtest    += expE_kbT;
      wtest_sq += expE_kbT * expE_kbT;

      // Delete one atom again
      lmp->input->one(group_test_particle_command);   
      lmp->input->one(delete_command);

      // Post-deletion breakdown
      if (lmp->force->kspace) lmp->force->kspace->qsum_qsq();
      if (lmp->force->pair->tail_flag) lmp->force->pair->reinit();

    }

    // collect chemical potential & stats
    beta_mu = -1.0 * log( wtest / (double) ninserts );

    wtest_sum = wtest / (double) ninserts;
    wtest_sq_sum = wtest_sq / (double) ninserts;
    var_wtest = wtest_sq_sum - wtest_sum * wtest_sum;
    stdev_beta_mu = sqrt( var_wtest ) / wtest;

    nsamples++;

    if (me == 0) {
      output.add_line( outfile, nsamples, beta, beta_mu, stdev_beta_mu );
      output.add_line( std::cout, nsamples, beta, beta_mu, stdev_beta_mu );
    }

    // Increment nstep 
    nstep++;

  }


  delete [] r;

  file.close();

  if (lammps == 1) delete lmp;

  // close down MPI

  MPI_Finalize();
}

