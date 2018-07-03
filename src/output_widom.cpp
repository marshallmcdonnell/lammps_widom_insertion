/*_._._._._._._._._._._._._._._._._._._._._

* File Name : output_widom.cpp

* Purpose :

* Creation Date : 14-10-2015

* Last Modified : Thu 22 Oct 2015 02:09:29 PM EDT

* Created by :

_._._._._._._._._._._._._._._._._._._._._._ */

#include "output_widom.h"

OutputWidom::OutputWidom(){ };

void OutputWidom::header( std::ostream& stream )
{
//  stream << "Timestep   " << "β*μ       " << "σ(βμ)       " << "μ " << std::endl;
  stream << "-------------------------------------------" << std::endl;
  stream << "Timestep   " << "BetaMu " << "STD(BetaMu)   " << "Mu " << std::endl;
  return;
}

void OutputWidom::add_line( std::ostream& stream, int step, double beta, double beta_mu, double stdev_beta_mu )
{
  stream << step << "        " << beta_mu << " " << stdev_beta_mu << " " << beta_mu / beta << std::endl;
  return;
}


