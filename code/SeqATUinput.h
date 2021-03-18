#ifndef SEQATUINPUT_H
#define	SEQATUINPUT_H
#include<iostream>
#include<vector>
#include<string>
#include"degradation_reverse.h"
#include"degradation_forward.h"
#include"first_last_read.h"
#include"position_expression.h"
#include"motification_of_maximalATUcluster.h"
#include<fstream>
#include<sstream>
#include<stdlib.h>
#include<map>
#include<algorithm>
#include<ctime>
using namespace std;

class SeqATU_Input
{
public:
vector<vector<string> > Vec_seqatu_input;

SeqATU_Input();
~SeqATU_Input();

vector<vector<string> > ATUinput(vector<string>&vec_1,vector<string>&vec_2,
									vector<string>&vec_3,double&Cut_1,double&Cut_2);

};

SeqATU_Input::SeqATU_Input()
{
 return;
}
SeqATU_Input::~SeqATU_Input()
{
 return; 
}


vector<vector<string> > SeqATU_Input::ATUinput(vector<string>&vec_1,vector<string>&vec_2,
									vector<string>&vec_3,double&Cut_1,double&Cut_2)
{

vector<string> vec_input_1;
vec_input_1=vec_1;
vector<string> vec_input_2;
vec_input_2=vec_2;
vector<string> vec_input_3;
vec_input_3=vec_3;
double cutoff_1=Cut_1;
double cutoff_2=Cut_2;

Dec_R R;
Dec_F F;

string temp;
string temp_1;
istringstream istr;

string temp_2;
string temp_3;
string temp_4;
string temp_5;
string temp_6;
string temp_7;
string temp_8;
string temp_9;
string temp_gene;
string str;
vector<string> vec_hash;
map<string,vector<string> > hash;
map<string,vector<string> >::iterator iter;
stringstream ss;

vector<int> vec_location_reverse;
vec_location_reverse.push_back(0);
vector<int> vec_location_forward;
vec_location_forward.push_back(0);
map<int,string> hash_reverse;
map<string,int> hash_tag_reverse;
map<string,int>::iterator iter_2;
map<int,string> hash_forward;
map<string,int> hash_tag_forward;
int tag_reverse=1;
int tag_forward=1;

map<string,string> hash_1;
map<string,string>::iterator iter_0227;//add 0227
string temp_Dbxref;
map<string,int> hash_forward_c;//gene-->tag
map<string,int>::iterator iter_c;
map<string,int> hash_reverse_c;
int tag_1=0;
int tag_2=0;

//opposite gene
int flag_gene=0;
map<string,int> hash_flag;
map<int,string> hash_strand;

//location
vector<string> vec_samread;
vector<vector<int> > Vec_firstlast_read;


for(int i=0;i<vec_input_1.size();i++)
{
	temp=vec_input_1[i];
	if(temp[0]!='@')
	{
		vec_samread.push_back(temp);
	}
}



First_Last_Read FLR;
Vec_firstlast_read=FLR.firstlast_readmain(vec_samread);

for(int e=0;e<vec_input_2.size();e++)
{
	temp=vec_input_2[e];
	if(temp[0]!='#')
	{
		istr.str(temp);
		istr>>temp_1>>temp_2>>temp_3>>temp_4>>temp_5>>temp_6>>temp_7>>temp_8>>temp_9;
		istr.clear();
		if((temp_3=="gene")||(temp_3=="pseudogene"))//find gene row
		{
			flag_gene++;
			vec_hash.push_back(temp_4);
			vec_hash.push_back(temp_5);
			if(temp_7=="+")
			{
				vec_hash.push_back("forward");
			}
			if(temp_7=="-")
			{
				vec_hash.push_back("reverse");
			}
			for(int i=0;i<temp_9.size();i++)
			{
				if(temp_9[i]==';')
				{
					temp_9[i]=' ';
				}
				if(temp_9[i]=='=')
				{
					temp_9[i]=' ';
				}

			}

			istr.str(temp_9);
			while(istr>>temp_1)
			{
				//if(temp_1=="gene")
				if(temp_1=="Name")     //change 0227
				{
					istr>>temp_gene;
					hash[temp_gene]=vec_hash;
					vec_hash.clear();
					break;
				}
				if(temp_1=="Dbxref")                  //maximal ATU cluster
				{
					istr>>temp_Dbxref;
					//hash_1[temp_Dbxref]=temp_gene;    // !!!!!20210227 best to test
				}
			}
			istr.clear();
			hash_1[temp_Dbxref]=temp_gene;            //maximal ATU cluster
			
			hash_flag[temp_gene]=flag_gene;           //opposite gene insert
	
			ss << temp_4.size();
			ss >> str;
			ss.clear();
			str.push_back('|');
			str.push_back(temp_4[0]);
			if(temp_7=="-")
			{
				iter_2=hash_tag_reverse.find(str);
				if(iter_2==hash_tag_reverse.end())
				{
					hash_tag_reverse[str]=tag_reverse;
				}
				hash_reverse[tag_reverse]=temp_gene;
				tag_reverse=tag_reverse+2;
				vec_location_reverse.push_back(atoi(temp_4.c_str()));
				vec_location_reverse.push_back(atoi(temp_5.c_str()));
				
				iter_c=hash_reverse_c.find(temp_gene);
				if(iter_c==hash_reverse_c.end())
				{
					tag_2++;
					hash_reverse_c[temp_gene]=tag_2;
				}
				if(iter_c!=hash_reverse_c.end())
				{
				}
				
				hash_strand[flag_gene]="reverse";
				
			}
			if(temp_7=="+")
			{
				iter_2=hash_tag_forward.find(str);
				if(iter_2==hash_tag_forward.end())
				{
					hash_tag_forward[str]=tag_forward;
				}
				hash_forward[tag_forward]=temp_gene;
				tag_forward=tag_forward+2;
				vec_location_forward.push_back(atoi(temp_4.c_str()));
				vec_location_forward.push_back(atoi(temp_5.c_str()));
				
				iter_c=hash_forward_c.find(temp_gene);
				if(iter_c==hash_forward_c.end())
				{
					tag_1++;
					hash_forward_c[temp_gene]=tag_1;
				}
				if(iter_c!=hash_forward_c.end())
				{
				}
				
				hash_strand[flag_gene]="forward";
			}
		}
	}
}




//*********genenitc and intergenic region expression*****************
vector<int>::iterator iter_max;
int min=0;
int max=0;
vector<int> vec_bp;
P_Expression PE;


vector<string> vec_cluster;
vector<string> vec_cluster_new;
vector<vector<string> > Vec_Cluster;
vector<vector<string> > Vec_Cluster_new;
vector<string> vec_hash_new;
vector<int> vec_flag;
vector<int> vec_location;
vector<int>::iterator iter_1;
int count=0;
int tag_include;//indicate if cluster have case which genes containing
double Average_add;
vector<double> result;
vector<vector<double> > Vec_expression;//expression
vector<double> vec_expression;
vector<vector<int> > Vec_genelocation;//location
vector<int> vec_genelocation;
vector<double> vec_average;//average
vector<int> vec_length;


iter_max=max_element(Vec_firstlast_read[1].begin(),Vec_firstlast_read[1].end());// change
min=1;


vec_bp=PE.eachbase_readcount(min,*iter_max,Vec_firstlast_read[0],Vec_firstlast_read[1]);



//*******************maximal ATU cluster*******************
vector<vector<string> > Vec_inital_cluster;//initial cluster 
vector<string> vec_inicluster;
int flag_cluster;
vector<string> vec_temp;
vector<vector<string> > Vec_forward_c;
vector<vector<string> > Vec_reverse_c;
int location_c;
for(int i=0;i<=hash_forward_c.size();i++)
{
	vec_inicluster.clear();
	vec_inicluster.push_back("0");
	Vec_forward_c.push_back(vec_inicluster);
}
	
for(int i=0;i<=hash_reverse_c.size();i++)
{
	vec_inicluster.clear();
	vec_inicluster.push_back("0");
	Vec_reverse_c.push_back(vec_inicluster);
}



for(int e=0;e<vec_input_3.size();e++)
{
	temp=vec_input_3[e];
	vec_temp.clear();
	flag_cluster=0;
	istr.str(temp);
	while(istr>>temp_1)
	{
		flag_cluster++;
		if(flag_cluster==4)//strand
		{
			temp_4=temp_1;
		}
		if(flag_cluster>=6)//maximal ATU cluster
		{	
			vec_temp.push_back(temp_1);
		}
	}
	istr.clear();
	
	vec_inicluster.clear();	
	//for(int i=0;i<vec_temp.size();i++)
	//{
		//vec_inicluster.push_back(hash_1[vec_temp[i]]);
	//}
	//add 0227
	for(int i=0;i<vec_temp.size();i++)
	{
		iter_0227=hash_1.find(vec_temp[i]);
		if(iter_0227!=hash_1.end())
		{
			vec_inicluster.push_back(hash_1[vec_temp[i]]);
		}
		if(iter_0227==hash_1.end())
		{
			vec_inicluster.push_back(vec_temp[i]);
		}
	}
	//add 0227
	
	if(temp_4=="+")
	{
		vec_inicluster.push_back("forward");
	}
	if(temp_4=="-")
	{
		vec_inicluster.push_back("reverse");
	}
	Vec_inital_cluster.push_back(vec_inicluster);
	vec_inicluster.clear();
}

for(int i=0;i<Vec_inital_cluster.size();i++)
{
	vec_inicluster.clear();
	vec_inicluster=Vec_inital_cluster[i];
	if(vec_inicluster[vec_inicluster.size()-1]=="forward")
	{
		iter_c=hash_forward_c.find(vec_inicluster[0]);
		if(iter_c!=hash_forward_c.end())
		{
			location_c=iter_c->second;
			Vec_forward_c[location_c]=vec_inicluster;
		}
	}
	if(vec_inicluster[vec_inicluster.size()-1]=="reverse")
	{
		iter_c=hash_reverse_c.find(vec_inicluster[0]);
		if(iter_c!=hash_reverse_c.end())
		{
			location_c=iter_c->second;
			Vec_reverse_c[location_c]=vec_inicluster;
		}
	}
}
//update Vec_inital_cluster;
Vec_inital_cluster.clear();

for(int i=1;i<Vec_forward_c.size();i++)
{
	vec_inicluster.clear();
	vec_inicluster=Vec_forward_c[i];
	if(vec_inicluster.size()>1)
	{
		Vec_inital_cluster.push_back(vec_inicluster);
	}
}
for(int i=1;i<Vec_reverse_c.size();i++)
{
	vec_inicluster.clear();
	vec_inicluster=Vec_reverse_c[i];
	if(vec_inicluster.size()>1)
	{
		Vec_inital_cluster.push_back(vec_inicluster);
	}
}	
//new Vec_inital_cluster	




//***************modification of maximal ATU cluster*********
MATUcluster MATU;
vector<vector<string> > Vec_inital_cluster_new;
Vec_inital_cluster_new=MATU.maxATU_main(hash,Vec_inital_cluster,cutoff_1,cutoff_2,hash_flag,hash_strand,vec_bp);

Vec_Cluster=Vec_inital_cluster_new;


int qiqi=0;
//vec_cluster,the name of these genes are stored
for(int k=0;k<Vec_Cluster.size();k++)
{	
	vec_cluster.clear();
	vec_genelocation.clear();
	vec_cluster=Vec_Cluster[k];
	for(int i=0;i<vec_cluster.size();i++)// each gene group
	{
		iter=hash.find(vec_cluster[i]);
		if(iter!=hash.end())
		{
			vec_hash_new=iter->second;
			if(vec_hash_new[2]=="forward")
			{
				vec_flag.push_back(1);
			}
			if(vec_hash_new[2]=="reverse")
			{
				vec_flag.push_back(2);
			}
			vec_location.push_back(atoi(vec_hash_new[0].c_str()));
			vec_location.push_back(atoi(vec_hash_new[1].c_str()));
			vec_hash_new.clear();
		}
	}
	// get the position of the genes : vec_location
	for(int i=0;i<vec_flag.size();i++)
	{
		count=count+vec_flag[i];
	}
	iter_1=min_element(vec_location.begin(),vec_location.end());
	min=*iter_1;
	iter_1=max_element(vec_location.begin(),vec_location.end());
	max=*iter_1;

	//remove genes containing
	tag_include=0;
	for(int i=2;i<vec_location.size();i=i+2)
	{
		if((vec_location[i-2]<=vec_location[i])&&(vec_location[i]<=vec_location[i-1]))
		{
			if((vec_location[i-2]<=vec_location[i+1])&&(vec_location[i+1]<=vec_location[i-1]))
			{
				tag_include=1;
			}
		}
	}
	if(tag_include==0)
	{
		result=PE.eachgene_readcount(min,max,vec_location,vec_bp);
		Average_add=result[result.size()-1];
		result.pop_back();
//		vec_average.push_back(Average_add); 

		for(int i=1;i<vec_location.size();i++)
		{
			vec_length.push_back(vec_location[i]-vec_location[i-1]);
		}

		if(count==2*vec_flag.size())//reverse
		{
			vec_average.push_back(Average_add);
			vec_cluster_new.clear();
			vec_genelocation.clear();
			vec_expression.clear();
			for(int i=vec_cluster.size()-2;i>=0;i--)
			{
				vec_cluster_new.push_back(vec_cluster[i]);
			}
			vec_cluster_new.push_back(vec_cluster[vec_cluster.size()-1]);
			
			Vec_Cluster_new.push_back(vec_cluster_new);//new cluster

			vec_cluster_new.clear();
		
			for(int i=0;i<vec_location.size();i++)
			{
				vec_genelocation.push_back(vec_location[i]);
			}
			Vec_genelocation.push_back(vec_genelocation);//gene location
			vec_genelocation.clear();
		
			for(int i=result.size()-1;i>0;i=i-2)
			{
				vec_expression.push_back(result[i]);
			}
			for(int i=result.size()-2;i>0;i=i-2)
			{
				vec_expression.push_back(result[i]);
			}
			Vec_expression.push_back(vec_expression);//gene and intergenic expression
			vec_expression.clear();
		}
		
		
		
		if(count==vec_flag.size())//forward
		{
			vec_average.push_back(Average_add);
			vec_cluster_new.clear();
			vec_genelocation.clear();
			vec_expression.clear();
			Vec_Cluster_new.push_back(vec_cluster);
			
			for(int i=0;i<vec_location.size();i++)
			{
				vec_genelocation.push_back(vec_location[i]);
			}			
			Vec_genelocation.push_back(vec_genelocation);//gene location
			vec_genelocation.clear();
		
			for(int i=1;i<result.size();i=i+2)
			{
				vec_expression.push_back(result[i]);
			}
			for(int i=2;i<result.size();i=i+2)
			{
				vec_expression.push_back(result[i]);
			}
			Vec_expression.push_back(vec_expression);//gene and intergenic expression
			vec_expression.clear();
		}	
	}
	vec_cluster_new.clear();
	vec_cluster.clear();
	vec_genelocation.clear();
	vec_expression.clear();
	vec_length.clear();	
	vec_location.clear();
	vec_flag.clear();
	vec_cluster.clear();
	result.clear();
	count=0;
}

Vec_Cluster.clear();  



//******genenitc and intergenic(gap) degradation*************************
vector<vector<float> > Vec_gene_degradation;
vector<vector<float> > Vec_gap_degradation;
vec_genelocation.clear();
vec_cluster_new.clear();
for(int i=0;i<Vec_genelocation.size();i++)
{
	vec_genelocation.clear();
	vec_genelocation=Vec_genelocation[i];
	vec_cluster_new.clear();
	vec_cluster_new=Vec_Cluster_new[i];
	
	F.vec_fd_1.clear();
	F.vec_fd_2.clear();
	R.vec_fd_1.clear();
	R.vec_fd_2.clear();
	
	if(vec_cluster_new[vec_cluster_new.size()-1]=="forward")
	{
		F.dec_main(vec_genelocation);
		Vec_gene_degradation.push_back(F.vec_fd_1);
		Vec_gap_degradation.push_back(F.vec_fd_2);
		vec_genelocation.clear();
	}

	if(vec_cluster_new[vec_cluster_new.size()-1]=="reverse")
	{
		R.dec_main(vec_genelocation);
		Vec_gene_degradation.push_back(R.vec_fd_1);
		Vec_gap_degradation.push_back(R.vec_fd_2);
		vec_genelocation.clear();
	}
}





//************convec input************************
int gene_number;
int s;
int t;
int a;
int b;
int flag;
int c;
int col;
int p;
float float_2;
vector<string> vec;
vector<vector<string> > Vec_ce;
vector<vector<string> > Vec_ci;
vector<vector<string> > Vec_A;
vector<vector<string> > Vec_B;
vector<vector<string> > Vec_u;
vector<vector<string> > Vec_v;


vector<string> vec_time_0924;
for(int x=0;x<Vec_Cluster_new.size();x++)
{

	
	Vec_seqatu_input.push_back(Vec_Cluster_new[x]);
	vec_time_0924.clear();
	gene_number=Vec_Cluster_new[x].size()-1;
	s=0.5*gene_number*(gene_number+1);
	t=2*gene_number-1;
	for(int i=0;i<s;i++)
	{
		for(int j=0;j<s+t;j++)
		{
			if(i==j)
			{
				vec_time_0924.push_back("0.00000001");
			}
			if(i!=j)
			{
				vec_time_0924.push_back("0.00");
			}
		}
	}

	
	for(int i=s;i<s+t;i++)
	{
		for(int j=0;j<s+t;j++)
		{
			if(i==j)
			{
				vec_time_0924.push_back("2.00");
			}
			if(i!=j)
			{
				vec_time_0924.push_back("0.00");
			}
		}
	}//G
	Vec_seqatu_input.push_back(vec_time_0924);
	vec_time_0924.clear();
	string temp_0924;
    stringstream ss_0924;  
	

	
	for(int i=0;i<s;i++)
	{
		ss_0924<<vec_average[x];
		ss_0924>>temp_0924;
		ss_0924.clear();
		vec_time_0924.push_back(temp_0924);
	}
	for(int i=s;i<s+t;i++)
	{
		vec_time_0924.push_back("0.00");
	}
	Vec_seqatu_input.push_back(vec_time_0924);
	vec_time_0924.clear();
	
	a=1;
	b=0;
	flag=gene_number;
	vec.clear();
	Vec_u.clear();
	vec.push_back("0.00");
	Vec_u.push_back(vec);
	for(int j=0;j<Vec_gene_degradation[x].size();j++)
	{
		ss.clear();
		ss<<Vec_gene_degradation[x][j];
		ss>>temp_1;
		ss.clear();
		vec.push_back(temp_1);
		b++;
		if(b==flag)
		{
			Vec_u.push_back(vec);
			vec.clear();
			vec.push_back("0.00");
			for(int i=0;i<a;i++)
			{
				vec.push_back("0.00");
			}
			a++;
			flag=a*gene_number-0.5*a*(a-1);
		}
	}//u


	
	vec.clear();
	Vec_A.clear();
	for(int i=0;i<=gene_number;i++)
	{
		for(int j=0;j<=s;j++)
		{
			vec.push_back("0.00");
		}
		Vec_A.push_back(vec);
		vec.clear();
	}//initial A

	for(int i=1;i<=gene_number;i++)
	{
		c=1;
		for(int j=1;j<=gene_number;j++)
		{
			col=gene_number-j+1;
			for(int k=1;k<=col;k++)
			{
				if(k<=i-j)
				{
					c++;
				}
				if(k>i-j)
				{
					Vec_A[i][c]=Vec_u[j][i];
					c++;
				}
			}
		}
	}//A


	
	a=1;
	b=0;
	flag=gene_number;
	vec.clear();
	Vec_v.clear();
	vec.push_back("0.00");
	Vec_v.push_back(vec);
	vec.push_back("0.00");
	b++;
	for(int j=0;j<Vec_gap_degradation[x].size();j++)
	{
	    ss.clear();
		ss<<Vec_gap_degradation[x][j];
		ss>>temp_1;
		ss.clear();
		vec.push_back(temp_1);
		b++;
		if(b==flag)
		{
			Vec_v.push_back(vec);
			vec.clear();
			vec.push_back("0.00");
			for(int i=0;i<a;i++)
			{
				vec.push_back("0.00");
			}
			a++;
			flag=a*(gene_number)-0.5*a*(a-1);
			vec.push_back("0.00");
			b++;
		}
	}
		
	vec.clear();
	for(int i=0;i<gene_number+1;i++)
	{
		vec.push_back("0.00");
	}
	Vec_v.push_back(vec);
	vec.clear();
	//v

	
	vec.clear();
	Vec_B.clear();
	for(int i=0;i<=gene_number;i++)
	{
		for(int j=0;j<=s;j++)
		{
			vec.push_back("0.00");
		}
		Vec_B.push_back(vec);
		vec.clear();
	}//initial B

	for(int i=1;i<=gene_number;i++)
	{
		c=1;
		for(int j=1;j<=gene_number;j++)
		{
			col=gene_number-j+1;
			for(int k=1;k<=col;k++)
			{
				if(k<=i-j)
				{
					c++;
				}
				if(k>i-j)
				{
					Vec_B[i][c]=Vec_v[j][i];
					c++;
				}
			}
		}
	}//B



	vec.clear();
	for(int i=1;i<=t;i++)
	{
		if(i<=gene_number)
		{
			for(int j=1;j<=s;j++)
			{
				vec.push_back(Vec_A[i][j]);
			}
		}
		if(i>gene_number)
		{
			for(int j=1;j<=s;j++)
			{
				vec.push_back(Vec_B[i-gene_number+1][j]);
			}
		}
		for(int j=1;j<=t;j++)
		{
			if(i==j)
			{
				vec.push_back("-1.00");
			}
			if(i!=j)
			{
				vec.push_back("0.00");
			}
		}
		Vec_ce.push_back(vec);
		vec.clear();
	}//Ce

	vec_time_0924.clear();
	for(int i=0;i<s+t;i++)
	{
		for(int j=0;j<t;j++)
		{
			vec_time_0924.push_back(Vec_ce[j][i]);
		}
	}
	Vec_seqatu_input.push_back(vec_time_0924);
	vec_time_0924.clear();
	Vec_ce.clear();
	
	double expression_0924;
	vec_time_0924.clear();
	for(int j=0;j<Vec_expression[x].size();j++)
	{
		expression_0924=-1*Vec_expression[x][j];
		ss_0924.clear();
		ss_0924<<expression_0924;
		ss_0924>>temp_0924;
		vec_time_0924.push_back(temp_0924);
	}
	Vec_seqatu_input.push_back(vec_time_0924);
	vec_time_0924.clear();
	//Ge0

	vec.clear();
	p=0;
	for(int i=0;i<s;i++)
	{
		for(int j=0;j<s+t;j++)
		{
			if(i==j)
			{
				vec.push_back("1.00");
			}
			if(i!=j)
			{
				vec.push_back("0.00");
			}
		}
		Vec_ci.push_back(vec);
		vec.clear();
		p++;
	}

	
	vec_time_0924.clear();
	ss_0924.clear();
	ss_0924<<p;
	ss_0924>>temp_0924;
	vec_time_0924.push_back(temp_0924);
	Vec_seqatu_input.push_back(vec_time_0924);
	vec_time_0924.clear();
	
	for(int i=0;i<s+t;i++)
	{
		for(int j=0;j<p;j++)
		{
			vec_time_0924.push_back(Vec_ci[j][i]);
		}
	}
	Vec_seqatu_input.push_back(vec_time_0924);
	vec_time_0924.clear();
	Vec_ci.clear();
	
	for(int i=0;i<s;i++)
	{
		vec_time_0924.push_back("0.00");
	}
	Vec_seqatu_input.push_back(vec_time_0924);
	vec_time_0924.clear();

}

return Vec_seqatu_input;
}

#endif