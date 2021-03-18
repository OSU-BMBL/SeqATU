#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <vector> 
#include <iomanip>  
#include "QuadProg++.hh"
#include "SeqATUinput.h"
#include "quadprog_out.h" 

using namespace QuadProgPP;

int main (int argc, char *const argv[]) 
{
	Matrix<double> G, CE, CI;
	Vector<double> g0, ce0, ci0, x;
	int n, m, p, s, t, i, j, genenumber;
	double sum = 0.0;
	char ch;

	std::vector<string> vec_input_1;
	std::vector<string> vec_input_2;
	std::vector<string> vec_input_3;

	std::string temp;
	std::string temp_1;
	std::ifstream ifs;
	std::istringstream istr;
	std::ofstream ofs(argv[4]);


	std::vector<string> vec_gene;
	ATU_Out A;
	std::vector<std::string> vec_x;
	std::string temp_0924;
	std::stringstream ss_0924; 


//************SeqATU--preprocessing*********
	ifs.open(argv[1]);
	while(getline(ifs,temp))
	{
		vec_input_1.push_back(temp);
	}
	ifs.close();

	ifs.open(argv[2]);
	while(getline(ifs,temp))
	{
		vec_input_2.push_back(temp);
	}
	ifs.close();

	ifs.open(argv[3]);
	while(getline(ifs,temp))
	{
		vec_input_3.push_back(temp);
	}
	ifs.close();

	double cutoff_1;
	double cutoff_2;
	cutoff_1=atof(argv[5]);
	cutoff_2=atof(argv[6]);



	std::vector<std::vector<string> > Vec_seqatu_input;
	SeqATU_Input SeqATU;
	Vec_seqatu_input=SeqATU.ATUinput(vec_input_1,vec_input_2,vec_input_3,cutoff_1,cutoff_2);
//************SeqATU--preprocessing*********

	std::cout<<"Preprocessing was completed successfully!"<<"\n";
//************SeqATU--convex quadratic programming*********
	for(int e=0;e<Vec_seqatu_input.size();e++)
	{
		vec_gene=Vec_seqatu_input[e];

		
		if(Vec_seqatu_input[e].size()==2)
		{
			ofs<<Vec_seqatu_input[e][0]<<"\t"<<Vec_seqatu_input[e][1]<<"\n";
		}
		if(Vec_seqatu_input[e].size()>2)
		{
			ofs<<Vec_seqatu_input[e][0];
			for(int z=1;z<Vec_seqatu_input[e].size();z++)
			{
				ofs<<"\t"<<Vec_seqatu_input[e][z];
			}
			ofs<<"\n";
		}
		genenumber=0;
		genenumber=Vec_seqatu_input[e].size();
		genenumber=genenumber-1;
		s=0.5*genenumber*(genenumber+1);
		t=2*genenumber-1;
		n=s+t;
		m=t;
		e=e+1;
		
		G.resize(n, n);
		{
			i=0;j=0;
			for(int z=0;z<Vec_seqatu_input[e].size();z++)
			{
				temp_1=Vec_seqatu_input[e][z];
				if(i!=n)
				{
					G[i][j]=atof(temp_1.c_str());
					j++;
					if(j==n)
					{
						i++;
						j=0;
					}
					
				}
			}	
		}
		
			
		
		e=e+1;
		g0.resize(n);
		{
			i=0;
			for(int z=0;z<Vec_seqatu_input[e].size();z++)
			{
				temp_1=Vec_seqatu_input[e][z];
				if(i!=n)
				{
					g0[i]=atof(temp_1.c_str());
					i++;
				}
			}		
			
		}


		e=e+1;
		CE.resize(n, m);
		{
			i=0;j=0;
			for(int z=0;z<Vec_seqatu_input[e].size();z++)
			{
				temp_1=Vec_seqatu_input[e][z];
				if(i!=n)
				{
					CE[i][j]=atof(temp_1.c_str());
					j++;
					if(j==m)
					{
						i++;
						j=0;
					}
				}
			}
			
		} 
		
		e=e+1;
		ce0.resize(m);
		{
			i=0;
			for(int z=0;z<Vec_seqatu_input[e].size();z++)
			{
				temp_1=Vec_seqatu_input[e][z];
				if(i!=m)
				{
					ce0[i]=atof(temp_1.c_str());
					i++;
				}
			}
		}

		e=e+1;
		temp=Vec_seqatu_input[e][0];
		p=atoi(temp.c_str());
		e=e+1;
		CI.resize(n, p);
		{
			i=0;j=0;
			for(int z=0;z<Vec_seqatu_input[e].size();z++)
			{
				temp_1=Vec_seqatu_input[e][z];
				if(i!=n)
				{
					CI[i][j]=atof(temp_1.c_str());
					j++;
					if(j==p)
					{
						i++;
						j=0;
					}
				}
			}
		}
		
	
		e=e+1;
		ci0.resize(p);
		{
			i=0;
			for(int z=0;z<Vec_seqatu_input[e].size();z++)
			{
				temp_1=Vec_seqatu_input[e][z];
				if(i!=p)
				{
					ci0[i]=atof(temp_1.c_str());
					i++;
				}
			}
		}

		x.resize(n);
    	solve_quadprog(G, g0, CE, ce0, CI, ci0, x);	
//************SeqATU--convex quadratic programming*********		




//************SeqATU--output of ATUs*********
		ss_0924.clear();
		if(x[0]>0.0000001)
		{
			ss_0924<<std::setprecision(15)<<x[0];
		}
		
		if(x[0]<=0.0000001)
		{
			ss_0924<<"0.000";
			x[0]=0.000;
		}
		for(int k=1;k<s;k++)
		{
			if(x[k]>0.0000001)
			{
				ss_0924<<","<<std::setprecision(15)<<x[k];
			}
		
			if(x[k]<=0.0000001)
			{
				ss_0924<<","<<"0.000";
				x[k]=0.000;
			}
		}
		for(int k=s;k<s+t;k++)
		{
			ss_0924<<","<<std::setprecision(15)<<x[k];
		}
		ss_0924>>temp_0924;
		ss_0924.clear();
		vec_x.clear();
		vec_x=A.ATU_main(vec_gene,temp_0924);
		for(int k=0;k<vec_x.size();k++)
		{
			ofs<<vec_x[k]<<"\t";
		}
		ofs<<"\n";
//************SeqATU--output of ATUs*********	
	
	
	}
	ifs.close();
	std::cout<<"SeqATU was completed successfully!"<<"\n";	
	vec_input_1.clear();
	vec_input_2.clear();
	vec_input_3.clear();
	vec_gene.clear();
	vec_x.clear();
}	


	
