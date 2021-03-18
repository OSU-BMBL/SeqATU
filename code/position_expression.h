#ifndef POSITIONEXPRESSION_H
#define POSITIONEXPRESSION_H
#include<iostream>
#include<vector>
#include<string>
#include<fstream>
#include<sstream>
#include<stdlib.h>
#include<map>
#include<algorithm>

using namespace std;
class P_Expression
{

public:

P_Expression();
~P_Expression();

vector<int> eachbase_readcount(int&F,int&L,vector<int>&Read_l,vector<int>&Read_r);
vector<double> eachgene_readcount(int&F,int&L,vector<int>&Vec_lo,vector<int>&Vec_bp);

};

P_Expression::P_Expression()
{
 return;
}
P_Expression::~P_Expression()
{
 return; 
}


//function1--->返回vec_bp
vector<int> P_Expression::eachbase_readcount(int&F,int&L,vector<int>&Read_l,vector<int>&Read_r)
{
int First=F;
int Last=L;
int leng=Last-First+1;
vector<int> Read_left=Read_l;
vector<int> Read_right=Read_r;
vector<int> vec_bp;
int start;
int end;


	for(int i=0;i<leng;i++)
	{
		vec_bp.push_back(0);
	}

	for(int i=0;i<Read_left.size();i++)
	{
		start=Read_left[i];
//cout<<"start:"<<start<<"\n";
		end=Read_right[i];
		if((First<=start)&&(start<=Last)&&(First<=end)&&(end<=Last))
		{
			for(int j=start-First;j<end-First+1;j++)
			{
				vec_bp[j]++;
			}
		}

		if((start<First)&&(First<=end)&&(end<=Last))
		{
			for(int j=0;j<end-First+1;j++)
			{
				vec_bp[j]++;
			}
		}		
		if((First<=start)&&(start<=Last)&&(Last<end))
		{
			for(int j=start-First;j<vec_bp.size();j++) 
			{
				vec_bp[j]++;
			}
		}
		if((start<First)&&(Last<end))
		{
			for(int j=0;j<vec_bp.size();j++)
			{
				vec_bp[j]++;
			}
		}
	}	
	return vec_bp;
}

//function2--->return vec_result
vector<double> P_Expression::eachgene_readcount(int&F,int&L,vector<int>&Vec_lo,vector<int>&Vec_bp)
{
int First=F;
int Last=L;	
vector<int> vec_bp=Vec_bp;
vector<int> Vec_location=Vec_lo;
vector<int> vec_bp_part;
vector<double> vec_result;


double sum_add=0;
double average_add;

vector<double> vec_time;//储存暂时的碱基表达量，用来求min 


	for(int i=0;i<Vec_location.size();i++)//vec_result的0元素无意义，输出时不会输出
	{
		vec_result.push_back(-1);
	}
	
	//int leng=Last-First+1;
	for(int i=First-1;i<Last;i++)
	{
		vec_bp_part.push_back(vec_bp[i]);
	}
	
	for(int k=0;k<vec_bp_part.size();k++)
	{
		sum_add=sum_add+vec_bp_part[k];
	}
	average_add=sum_add/vec_bp_part.size();

		
	int tag;
	for(int i=0;i<vec_bp_part.size();i++)
	{	
		tag=i+First;
		for(int j=0;j<Vec_location.size();j++)
		{
			if((Vec_location[j]<=tag)&&(tag<Vec_location[j+1]))
			{
				vec_result[j+1]=vec_result[j+1]+vec_bp_part[i];
				
			}	
		}
	}

	for(int i=1;i<vec_result.size();i++)
	{
		if(Vec_location[i]>Vec_location[i-1])
		{
			vec_result[i]=vec_result[i]/(Vec_location[i]-Vec_location[i-1]);
		}
	}

	for(int i=1;i<vec_result.size();i++)
	{
		if(Vec_location[i]<=Vec_location[i-1])
		{
			vec_result[i]=(vec_result[i-1]+vec_result[i+1])/2;
		}

	}	
	vec_result.push_back(average_add);
	vec_bp_part.clear();
	return vec_result;
}


#endif


