#ifndef IC_DEFECT_HEADER
#define IC_DEFECT_HEADER

#include "parser.hpp"
#include "global_defect.hpp"
#include "metadata.hpp"
#include "background.hpp"

using namespace std;
using namespace LATfield2;

void generateIC_defects(cosmology & cosmo, defects_metadata & defects_sim, const double fourpiG, DefectBase *defects, GlobalDefect  & defects_, double z_ic_defect,double z_in ,double dx, string h5filename)
{

	defects = &defects_; 

	defects_.generate_init_cond();

	COUT << "Starting the prevolution of defects!" << endl;

	double endz = z_in;
	int i =0;
	double tmp = 1 / (z_ic_defect + 1);
	double sum = 0;
	double sum_;
	int N = 10;


	ofstream phifile;
	if(parallel.rank() == 0)
	{
		phifile.open (h5filename + "/average_rho_phi_defect.txt",ios::trunc);
		phifile << "#i" << " " << "#z" << " " << "#a" << " " << "#adotovera" << " " << "#average phi" << " " << "#average rho" << endl;
		phifile.close();
	}


	do 
	{
		double vari = Hconf(tmp, fourpiG, cosmo);
		double DT = 0.001;
		defects->update_phi(&DT);
		rungekutta4bg(tmp, fourpiG, cosmo, DT/2.0);
		defects->update_pi(&DT,&tmp,&vari);
		rungekutta4bg(tmp, fourpiG, cosmo, DT/2.0);

		defects->compute_Tuv_defect(tmp);
		

		double averagephi = defects_.averagephi();
		double averagerho = defects_.averagerhodefect();
		i++;

		if(parallel.rank() == 0)
		{
			phifile.open (h5filename + "average_rho_phi_defect.txt",std::ios_base::app);
			phifile << i << " " << (1/tmp) - 1 << " " << tmp << " " << vari << " " << averagephi << " " << averagerho << endl;
			phifile.close();
		}

		sum += averagephi;
		if(i%N == 0)
		{
			sum = sum/N;
			sum_ = i;
			COUT << sum << " " ;
			if(tmp <= 1/(endz + 1))
			{
				sum = 0;
			}
		}

	}
	while (tmp <= 1/(endz + 1) );

	sum = sum /(i-sum_);
	double z = (1-tmp)/tmp;
	COUT << endl << endl << "The total steps for defects prevolution is: " << i << " and the redshift at the end of defect prevolution now is: " << z << endl;

	if(sum < 0.988 || sum > 1.02)
	{
		COUT << endl << COLORTEXT_RED << " error " << COLORTEXT_RESET << ": The defect has not reached stable condition !!!" << endl;
#ifdef LATFIELD2_HPP
		parallel.abortForce();
#endif
	}
	else if(sum >= 0.988 || sum <= 1.02)
	{
		defects_.write_Tuv_defect(h5filename, i);
		COUT << endl << COLORTEXT_BLUE << "The defect has reached stable condition in prevolution " << COLORTEXT_RESET << endl;
	}


}

#endif
