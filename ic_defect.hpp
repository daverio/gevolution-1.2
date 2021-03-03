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

	COUT << "Starting the prevolution of defects!" << endl;
	defects = &defects_; 
	double endz = 0;
	int i =0;
	double temp = 1 / (z_ic_defect + 1);
	double sum = 0;
	double sum_;
	int N = 10;


	ofstream phifile;
	if(parallel.rank() == 0)
	{
		phifile.open ("average_rho_phi_defect.txt",ios::trunc);
		phifile << "#i" << " " << "#z" << " " << "#a" << " " << "#adotovera" << " " << "#average phi" << " " << "#average rho" << endl;
		phifile.close();
	}


	do 
	{
		i++;
		double vari = Hconf(temp, fourpiG, cosmo);
		double DT = 0.001;
		defects->update_phi(&DT);
		rungekutta4bg(temp, fourpiG, cosmo, DT/2.0);
		defects->update_pi(&DT,&temp,&vari);
		rungekutta4bg(temp, fourpiG, cosmo, DT/2.0);
		double averagephi = defects_.averagephi();
		double averagerho = defects_.averagerhodefect(temp);
		

		if(parallel.rank() == 0)
		{
			phifile.open ("average_rho_phi_defect.txt",std::ios_base::app);
			phifile << i << " " << (1/temp) - 1 << " " << temp << " " << vari << " " << averagephi << " " << averagerho << endl;
			phifile.close();
		}

		sum += averagephi;
		if(i%N == 0)
		{
			sum = sum/N;
			sum_ = i;
			if(temp <= 1/(endz + 1))
			{
				sum = 0;
			}
		}
		
		if(i%1000 == 0)
		{ 
		    defects->compute_Tuv_defect(temp, h5filename, i);
		}
		
	}
	while (temp <= 1/(endz + 1) );

	sum = sum /(i-sum_);
	double z = (1-temp)/temp;
	COUT << endl << endl << "The total steps for defects prevolution is: " << i << " and the redshift at the end of defect prevolution now is: " << z << endl;

	if(sum < 0.988 || sum > 1.02)
	{
		defects->compute_Tuv_defect(temp, h5filename, i);
		COUT << endl << COLORTEXT_RED << " /!\ error " << COLORTEXT_RESET << ": The defect has not reached stable condition !!!" << endl;
#ifdef LATFIELD2_HPP
		parallel.abortForce();
#endif
	}
	else if(sum >= 0.988 || sum <= 1.02)
	{
		defects->compute_Tuv_defect(temp, h5filename, i);
		COUT << endl << COLORTEXT_BLUE << "The defect has reached stable condition in prevolution " << COLORTEXT_RESET << endl;
	}


}

#endif
