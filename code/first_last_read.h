#ifndef FIRSTLAST_H
#define FIRSTLAST_H
#include<iostream>
#include<vector>
#include<string>
#include<fstream>
#include<sstream>
#include<stdlib.h>

using namespace std;

class First_Last_Read
{
public:
vector<vector<int> > Vec_first_last;//output 0708
vector<string> vec_samread;//input 0708
string temp;
string temp_1;
string temp_2;
string temp_3;
string temp_4;
string temp_5;
string temp_6;
istringstream istr;
vector<int> read_start;
vector<int> read_end;
vector<char> s1_letter;
vector<int> s1_number;
int start;
int end;
int start_new;
int end_new;
int temp_2_flag;
int min;
int max;
int option;


First_Last_Read();
~First_Last_Read();

int Two_type_1(int m);
vector<vector<int> > firstlast_readmain(vector<string>&vec_samread_old);

};

First_Last_Read::First_Last_Read()
{
 return;
}
First_Last_Read::~First_Last_Read()
{
 return; 
}

vector<vector<int> > First_Last_Read::firstlast_readmain(vector<string>&vec_samread_old)
{
vec_samread=vec_samread_old;
int t_1=0;

for(int k=0;k<vec_samread.size();k++)
{
	temp=vec_samread[k];
	if(temp[0]!='@')
	{
		istr.str(temp);
		istr>>temp_1>>temp_2>>temp_3>>temp_4>>temp_5>>temp_6;
		istr.clear();
		if(temp_6!="*")
		{
			start=atoi(temp_4.c_str());
			temp_2_flag=Two_type_1(atoi(temp_2.c_str()));
			//read_start.push_back(start);
			
			end=start;
			for(int i=0;i<temp_6.size();i++)
			{
				if((temp_6[i]>='A')&&(temp_6[i]<='Z'))
				{
					s1_letter.push_back(temp_6[i]);
					s1_number.push_back(atoi(temp_6.substr(t_1,i-t_1).c_str()));
					t_1=i+1;
				}
			}
			for(int j=0;j<s1_letter.size();j++)
			{
			if((s1_letter[j]=='M')||(s1_letter[j]=='D'))
				{
					end=end+s1_number[j];
				}
			}
			end=end-1;
			s1_letter.clear();
			s1_number.clear();
			t_1=0;

			if(temp_2_flag==1)
			{
				temp=vec_samread[k+1];
				k=k+1;
				istr.str(temp);
				istr>>temp_1>>temp_2>>temp_3>>temp_4>>temp_5>>temp_6;
				istr.clear();
				if(temp_6!="*")
				{
					start_new=atoi(temp_4.c_str());
					temp_2_flag=Two_type_1(atoi(temp_2.c_str()));			
					end_new=start_new;	
					for(int i=0;i<temp_6.size();i++)
					{

						if((temp_6[i]>='A')&&(temp_6[i]<='Z'))
						{
							s1_letter.push_back(temp_6[i]);
							s1_number.push_back(atoi(temp_6.substr(t_1,i-t_1).c_str()));
							t_1=i+1;
						}
					}
					for(int j=0;j<s1_letter.size();j++)
					{
						if((s1_letter[j]=='M')||(s1_letter[j]=='D'))
						{
							end_new=end_new+s1_number[j];
						}
					}
					end_new=end_new-1;
					s1_letter.clear();
					s1_number.clear();
					t_1=0;
				}
			
			
				if(start>start_new)
				{
					start=start_new;
				}
				if(end_new>end)
				{
					end=end_new;
				}
		
		
			}
			
			read_start.push_back(start);
			read_end.push_back(end);	
		}
	}
}
Vec_first_last.push_back(read_start);
Vec_first_last.push_back(read_end);

read_start.clear();
read_end.clear();

return Vec_first_last;
}

int First_Last_Read::Two_type_1(int m)
{  
int M=m;
int u;
int w;
int r_2;
  vector<int> two_type_2;
  while(M!=0)
      {  u=M%2;
         two_type_2.push_back(u);
         M=(int)(M/2);
      }
if(two_type_2.size()<11)
   {
      w=11-two_type_2.size();
     for(int i=1;i<=w;i++)
        two_type_2.push_back(0);
   }
 r_2=two_type_2[1];
 two_type_2.clear();
 return r_2;
}


#endif



