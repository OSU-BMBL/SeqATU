#ifndef ATUOUT_H
#define	ATUOUT_H
#include<iostream>
#include<vector>
#include<string>
#include<sstream>
#include<fstream>
#include<stdlib.h>

using namespace std;


class ATU_Out
{
	public:
	vector<string> vec_atu;
	ATU_Out();
	~ATU_Out();
	vector<string> ATU_main(vector<string>&vec_G,string&temp_X);
};

ATU_Out::ATU_Out()
{
 return;
}
ATU_Out::~ATU_Out()
{
 return; 
}


vector<string> ATU_Out::ATU_main(vector<string>&vec_G,string&temp_X)
{

	vec_atu.clear();
	std::string temp_0924;
	std::stringstream ss_0924;

	string temp;
	string temp_1;
	istringstream istr;

	vector<int> vec_flag_1;
	vector<int> vec_flag_2;
	vector<string> vec_genename;
	vector<float> vec_out;
	int genenum;
	int count;

	vec_genename.clear();
	vec_out.clear();
	vec_flag_1.clear();
	vec_flag_2.clear();
	vec_genename=vec_G;
	vec_genename.pop_back();//remove forward or reverse
	genenum=vec_genename.size();
	
	temp=temp_X;
	for(int i=0;i<temp.size();i++)
	{
		if(temp[i]==',')
		{
			temp[i]=' ';
		}
	}
	istr.str(temp);
	while(istr>>temp_1)
	{
		vec_out.push_back(atof(temp_1.c_str()));
	}
	istr.clear();
	
	
	count=2*genenum-1;
	for(int i=1;i<=count;i++)
	{
		vec_out.pop_back();
	}//left are TU expression
	
	for(int i=1;i<=genenum;i++)
	{
		for(int j=1;j<=genenum-i+1;j++)
		{
			vec_flag_1.push_back(i);
		}
	}
	for(int i=1;i<=genenum;i++)
	{
		for(int j=1;j<=genenum-i+1;j++)
		{
			vec_flag_2.push_back(i+j-1);
		}
	}

	ss_0924.clear();
	for(int i=0;i<vec_out.size();i++)
	{
		if(vec_out[i]!=0)
		{
			for(int j=vec_flag_1[i];j<=vec_flag_2[i];j++)
			{
				if(j==vec_flag_1[i])	{ss_0924<<vec_genename[j-1];}
				if(j!=vec_flag_1[i]) {ss_0924<<"|"<<vec_genename[j-1];}
				
			}
			ss_0924>>temp_0924;
			ss_0924.clear();
			vec_atu.push_back(temp_0924);	
		}

	}

	vec_flag_1.clear();
	vec_flag_2.clear();
	vec_genename.clear();
	vec_out.clear();
	return vec_atu;
}

#endif
