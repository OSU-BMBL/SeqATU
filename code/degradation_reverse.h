#ifndef REVERSE_H
#define REVERSE_H
#include<iostream>
#include<vector>
#include<string>
#include<math.h>
#include<algorithm>
using namespace std;

class Dec_R
{

public:
vector<int> vec_part;//u
vector<float> vec_part_x;//x
vector<float> vec_out;//out
vector<float> vec_out_new;
vector<int> vec_location;//main u
vector<float> vec_fd_1;
vector<float> vec_fd_2;
vector<float> vec_main_out;
vector<int>::iterator iter;

Dec_R();
~Dec_R();

void dec_main(vector<int>&vec_A);
vector<float> dec_small(vector<int>&vec_B);

};

Dec_R::Dec_R()
{
 return;
}
Dec_R::~Dec_R()
{
 return; 
}

void Dec_R::dec_main(vector<int>&vec_A)
{
	vec_location=vec_A;
	iter=max_element(vec_location.begin(),vec_location.end());
if(*iter==vec_location[vec_location.size()-1])
{
		while(vec_location.size()!=0)
		{
			vec_main_out=dec_small(vec_location);
			for(int i=0;i<vec_main_out.size();i++)
			{
				if(i%2==0){vec_fd_1.push_back(vec_main_out[i]);}
				if(i%2!=0){vec_fd_2.push_back(vec_main_out[i]);}
			}
			vec_location.pop_back();
			vec_location.pop_back();
			vec_main_out.clear();
		}
		vec_location.clear();
}
	vec_out.clear();
	vec_out_new.clear();
	vec_location.clear();
	vec_main_out.clear();

}


vector<float> Dec_R::dec_small(vector<int>&vec_B)
{
vec_part=vec_B;//u
//���洴��matlab�ж�Ӧ��x
vec_out.clear();
vec_out_new.clear();
int k=vec_part.size();
int length=vec_part[k-1]-vec_part[0]+1;

float initial=vec_part[0];//x1
float finish=vec_part[vec_part.size()-1];//xmax

for(int i=0;i<length;i++)
{
	vec_part_x.push_back(initial+i);
}

for(int i=0;i<vec_part_x.size();i++)
{
	vec_part_x[i]=(vec_part_x[i]-initial)/(finish-initial)*1000;
}



/*
//x������ϣ����洴��d
float mean;
float sum=0;
float va;//variance
float std;//standard deviation
for(int i=0;i<vec_part_x.size();i++)
{
	sum=sum+vec_part_x[i];
}
mean=sum/vec_part_x.size();
sum=0;
for(int i=0;i<vec_part_x.size();i++)
{
	sum=sum+pow(vec_part_x[i]-mean,2);
}
va=sum/vec_part_x.size();//����
std=sqrt(va);//��׼��
for(int i=0;i<vec_part_x.size();i++)
{
	vec_part_x[i]=(vec_part_x[i]-mean)/std;
}
//���¹�����xΪd

*/


int width;
int tag=0;
float count=0;
for(int i=1;i<vec_part.size();i++)
{
	width=vec_part[i]-vec_part[i-1]+1;
	if(width>0)
	{
		for(int j=tag;j<tag+width;j++)
		{
			count=count+0.25563225*exp(0.00128*vec_part_x[j]);
		}
		vec_out.push_back(count/width);
	}
	if(width<=0)
	{
		vec_out.push_back(0);//��϶Ϊ0 or <0�Ĳ���,����Ϊ0������ȡǰ���ƽ��
	}
	tag=tag+width-1;
count=0;

}

for(int i=vec_out.size()-1;i>=0;i--)//���������
{
	vec_out_new.push_back(vec_out[i]);
}
for(int i=0;i<vec_out_new.size();i++)
{
	if(vec_out_new[i]==0)//width <=0 ��ȡƽ��
	{
		vec_out_new[i]=(vec_out_new[i-1]+vec_out_new[i+1])/2;
	}
}
	
	


//vec_out�������
vec_part.clear();
vec_part_x.clear();
vec_out.clear();

return vec_out_new;
}


#endif




