#ifndef Maximal_ATUcluster_motification_H
#define Maximal_ATUcluster_motification_H
#include<iostream>
#include<vector>
#include<string>
#include<sstream>
#include<fstream>
#include<stdlib.h>
#include<map>
#include<algorithm>

using namespace std;

class MATUcluster
{
public:
map<string,vector<string> > hash;
vector<vector<string> > Vec_inital_cluster;
//vector<string> vec_crossgeneread;
map<string,int> hash_flag_gff;
map<int,string> hash_strand;
vector<int> vec_bp;
double cutoff_1;
double cutoff_2;
istringstream istr;


MATUcluster();
~MATUcluster();

vector<vector<string> > maxATU_main(map<string,vector<string> >&Hash,vector<vector<string> >&Vec_Inital_Cluster,
									double&cutoff_1,double&cutoff_2,
									map<string,int>&Hash_flag_gff,map<int,string>&Hash_strand,
									vector<int>&vec_BP);

vector<vector<string> > cluster_combine(vector<vector<string> >&Vec_Inital_Cluster);
vector<vector<string> > oppositegene_split(vector<vector<string> >&Vec_Cluster_1);
vector<vector<string> > average_select(vector<vector<string> >&Vec_Cluster_2);


};


MATUcluster::MATUcluster()
{
 return;
}

MATUcluster::~MATUcluster()
{
 return; 
}

vector<vector<string> > MATUcluster::maxATU_main(map<string,vector<string> >&Hash,
												vector<vector<string> >&Vec_Inital_Cluster,
												double&Cutoff_1,double&Cutoff_2,
												map<string,int>&Hash_flag_gff,map<int,string>&Hash_strand,
												vector<int>&vec_BP)
{
	hash=Hash;
	Vec_inital_cluster=Vec_Inital_Cluster;
	//vec_crossgeneread=vec_Crossgeneread;
	cutoff_1=Cutoff_1;
	cutoff_2=Cutoff_2;
	hash_flag_gff=Hash_flag_gff;
	hash_strand=Hash_strand;
	vec_bp=vec_BP;
	vector<vector<string> > Vec_cluster_1;//combine之后的cluster
	Vec_cluster_1=cluster_combine(Vec_inital_cluster);//input Vec_1
//cout<<"1"<<"\n";
	vector<vector<string> > Vec_cluster_2;
	Vec_cluster_2=oppositegene_split(Vec_cluster_1);
//cout<<"2"<<"\n";
	vector<vector<string> > Vec_cluster_3;
	Vec_cluster_3=average_select(Vec_cluster_2);
//cout<<"3"<<"\n";
	return Vec_cluster_3;

}

vector<vector<string> > MATUcluster::cluster_combine(vector<vector<string> >&Vec_Inital_Cluster)
{
	vector<vector<string> > Vec_input_1=Vec_Inital_Cluster;
	vector<vector<string> > Vec_output_1;
	vector<vector<string> > Vec_forward;//
	vector<vector<string> > Vec_reverse;//
	stringstream ss;
	string temp;
	string temp_1;
	string temp_2;
	string temp_1_1;
	string temp_2_1;
	vector<string> vec_1;
	vector<string> vec_2;
	int left;
	int right;
	string left_gene;
	string right_gene;
	string time;
	vector<vector<string> > Vec_cluster;//
	int tag=0;//indicate 是否没有满足阈值
	vector<string> vec_time;//中间vector
	double A=cutoff_1;
	
	

	for(int i=0;i<Vec_input_1.size();i++)
	{
		vec_time.clear();
		vec_time=Vec_input_1[i];
		if(vec_time[vec_time.size()-1]=="forward")
		{
			Vec_forward.push_back(vec_time);
		}
		if(vec_time[vec_time.size()-1]=="reverse")
		{
			Vec_reverse.push_back(vec_time);
		}
	}
//cout<<"Vec_reverse.size():"<<Vec_reverse.size()<<"\n";
//cout<<"Vec_forward.size():"<<Vec_forward.size()<<"\n";	
/* 	map<string,string> hash_gene;
	map<string,string>::iterator iter_gene;
	map<string,int> hash_flag;
	map<string,int>::iterator iter_flag;
	
	for(int i=0;i<vec_crossgeneread.size();i++)
	{
		temp=vec_crossgeneread[i];
		for(int i =0;i<temp.size();i++)
		{
			if(temp[i]=='|')
			{
				temp[i]=' ';
			}
		}
		istr.str(temp);
		istr>>temp_1;
		while(istr>>temp_2)
		{
//cout<<"temp_1:"<<temp_1<<"\t"<<temp_2<<"\n";
			if(temp_1!=temp_2)
			{
				iter_gene=hash_gene.find(temp_1);
				if(iter_gene==hash_gene.end())
				{
					hash_gene[temp_1]=temp_2;
				}
				if(iter_gene!=hash_gene.end())
				{
					if(temp_2!=iter_gene->second)
					{
						//cout<<"error:"<<temp<<"\n";
					}
				}
				
				iter_flag=hash_flag.find(temp_1);
				if(iter_flag==hash_flag.end())
				{
					hash_flag[temp_1]=1;
				}
				if(iter_flag!=hash_flag.end())
				{
					hash_flag[temp_1]++;
				}
				temp_1=temp_2;	
			}
		}
		istr.clear();
	}//for */

//cout<<"hash_flag.size():"<<hash_flag.size()<<"\n";	
//cout<<"hash_gene.size():"<<hash_gene.size()<<"\n";	
	
//cout<<"forward"<<"\n";
	
	//forward
	vec_1=Vec_forward[0];
	Vec_cluster.push_back(vec_1);
	for(int i=1;i<Vec_forward.size();i++)
	{
		tag=0;
		vec_2=Vec_forward[i];
		
		left_gene=vec_1[vec_1.size()-2];//gap的左端基因
//cout<<"left_gene:"<<left_gene<<"\n";		
		vec_time.clear();
		vec_time=hash[left_gene];
		time=vec_time[1];
		left=atoi(time.c_str());
//cout<<"left:"<<left<<"\n";	
		
		right_gene=vec_2[0];
//cout<<"right_gene:"<<right_gene<<"\n";		
		vec_time.clear();
		vec_time=hash[right_gene];
		time=vec_time[0];
		right=atoi(time.c_str());
//cout<<"right:"<<right<<"\n";		
		vec_time.clear();
		//add
		/* iter_gene=hash_gene.find(left_gene);
		if(iter_gene!=hash_gene.end())
		{
			if(right_gene==iter_gene->second)//有cross-gene支持
			{
				tag=1;
			}
		} */
		//add
		if(right-left>=A)
		{
			if(Vec_cluster.size()==1)
			{
				vec_time.clear();
				vec_time=Vec_cluster[0];
				Vec_output_1.push_back(vec_time);
				Vec_cluster.clear();
			}
			if(Vec_cluster.size()>1)
			{
				vec_time.clear();
				for(int i=0;i<Vec_cluster.size();i++)
				{
					for(int j=0;j<Vec_cluster[i].size()-1;j++)
					{
						vec_time.push_back(Vec_cluster[i][j]);
					}
				}
				vec_time.push_back("forward");
				Vec_output_1.push_back(vec_time);
				vec_time.clear();
				Vec_cluster.clear();
			}
			Vec_cluster.push_back(vec_2);	
		}
		if(right-left<A)
		{
			Vec_cluster.push_back(vec_2);
		}
		vec_1.clear();
		vec_1=vec_2;
		vec_2.clear();
	}
	
	vec_time.clear();
	for(int i=0;i<Vec_cluster.size();i++)
	{
		for(int j=0;j<Vec_cluster[i].size()-1;j++)
		{
			vec_time.push_back(Vec_cluster[i][j]);
		}
	}
	vec_time.push_back("forward");
	Vec_output_1.push_back(vec_time);
	vec_time.clear();
	Vec_cluster.clear();
	
	
//cout<<"reverse"<<"\n";	
	
	
	//reverse
	vec_1=Vec_reverse[0];
	Vec_cluster.push_back(vec_1);
	
	for(int i=1;i<Vec_reverse.size();i++)
	{
		tag=0;
		vec_2=Vec_reverse[i];
		
		left_gene=vec_1[vec_1.size()-2];//gap的左端基因
		vec_time.clear();
		vec_time=hash[left_gene];
		time=vec_time[1];
		left=atoi(time.c_str());
		//left=hash_right[left_gene];
		
		right_gene=vec_2[0];
		vec_time.clear();
		vec_time=hash[right_gene];
		time=vec_time[0];
		right=atoi(time.c_str());
		//right=hash_left[right_gene];
		
		//add
		/* iter_gene=hash_gene.find(left_gene);
		if(iter_gene!=hash_gene.end())
		{
			if(right_gene==iter_gene->second)//有cross-gene支持
			{
				tag=1;
			}
		} */
		//add
		
		//if(right-left>=A)
		if(right-left>=A)
		{
			if(Vec_cluster.size()==1)
			{
				vec_time.clear();
				vec_time=Vec_cluster[0];
				Vec_output_1.push_back(vec_time);
				Vec_cluster.clear();
			}
			if(Vec_cluster.size()>0)
			{
				vec_time.clear();
				for(int i=0;i<Vec_cluster.size();i++)
				{
					for(int j=0;j<Vec_cluster[i].size()-1;j++)
					{
						vec_time.push_back(Vec_cluster[i][j]);
					}
				}
				vec_time.push_back("reverse");
				Vec_output_1.push_back(vec_time);
				vec_time.clear();
				Vec_cluster.clear();
			}
			Vec_cluster.push_back(vec_2);
			
		}

		//if(right-left<A)
		if(right-left<A)
		{
			Vec_cluster.push_back(vec_2);
		}
	
		vec_1.clear();
		vec_1=vec_2;
		vec_2.clear();
	}
	
	
	vec_time.clear();
	for(int i=0;i<Vec_cluster.size();i++)
	{
		for(int j=0;j<Vec_cluster[i].size()-1;j++)
		{
			vec_time.push_back(Vec_cluster[i][j]);
		}
	}
	vec_time.push_back("reverse");
	Vec_output_1.push_back(vec_time);
	vec_time.clear();
	Vec_cluster.clear();
	
	
	
	vec_1.clear();
	vec_2.clear();
	Vec_cluster.clear();
	Vec_forward.clear();
	Vec_reverse.clear();
	/* hash_gene.clear(); */
	/* hash_flag.clear(); */
	
	return Vec_output_1;
}

vector<vector<string> > MATUcluster::oppositegene_split(vector<vector<string> >&Vec_Cluster_1)
{
	vector<vector<string> > Vec_input_2=Vec_Cluster_1;
	vector<vector<string> > Vec_output_2;
	vector<string> vec_time;//中间向量
	string temp;
	string temp_1;
	string temp_2;
	map<string,int>::iterator iter_f;
	map<int,string>::iterator iter_s;
	vector<string> vec_temp;
	vector<int> vec_flag;
	string strand;
	vector<int> vec_strand;
	vector<int> vec_transcript;
	vector<int> vec_tag;//记录cluster需要拆分的位置
	int choose=0;
	
/* 	map<string,string> hash_gene;
	map<string,string>::iterator iter_gene;
	map<string,int> hash_flag_1;
	map<string,int>::iterator iter_flag;
	int tag; */
	
/* //***cross-gene read	
	for(int i=0;i<vec_crossgeneread.size();i++)
	{
		temp=vec_crossgeneread[i];
		for(int i =0;i<temp.size();i++)
		{
			if(temp[i]=='|')
			{
				temp[i]=' ';
			}
		}
		istr.str(temp);
		istr>>temp_1;
		while(istr>>temp_2)
		{
//cout<<"temp_1:"<<temp_1<<"\t"<<temp_2<<"\n";
			if(temp_1!=temp_2)
			{
				iter_gene=hash_gene.find(temp_1);
				if(iter_gene==hash_gene.end())
				{
					hash_gene[temp_1]=temp_2;
				}
				if(iter_gene!=hash_gene.end())
				{
					if(temp_2!=iter_gene->second)
					{
	//cout<<"errortemp_2:"<<temp_2<<"\t"<<iter_gene->second<<"\n";
						cout<<"error:"<<temp<<"\n";
					}
				}
				
				iter_flag=hash_flag_1.find(temp_1);
				if(iter_flag==hash_flag_1.end())
				{
					hash_flag_1[temp_1]=1;
				}
				if(iter_flag!=hash_flag_1.end())
				{
					hash_flag_1[temp_1]++;
				}
				temp_1=temp_2;	
			}
		}
		istr.clear();
	}//for */
	
	
//Maximal_ATUcluster update			
	for(int i=0;i<Vec_input_2.size();i++)
	{
		vec_temp.clear();
		vec_flag.clear();
		vec_tag.clear();
		vec_tag.push_back(0);
		vec_transcript.clear();
		//tag=0;
		for(int j=0;j<Vec_input_2[i].size();j++)
		{
			temp_1=Vec_input_2[i][j];
			iter_f=hash_flag_gff.find(temp_1);
			if(iter_f!=hash_flag_gff.end())
			{
				vec_temp.push_back(temp_1);
				vec_flag.push_back(iter_f->second);
			}
			if(iter_f==hash_flag_gff.end())
			{
				strand=temp_1;
			}
		}

		if(vec_flag.size()==1)//单基因cluster
		{
			if(choose==0)
			{
				vec_time.clear();
				vec_time.push_back(vec_temp[0]);
				vec_time.push_back(strand);
				Vec_output_2.push_back(vec_time);
				vec_time.clear();
			}
		}
		if(vec_flag.size()>1)
		{
//cout<<"right"<<"\n";
			for(int i=1;i<vec_flag.size();i++)
			{
				if(vec_flag[i]-vec_flag[i-1]!=1)//如果间隔大于1
				{                                                   //在这插
//cout<<vec_temp[i]<<"\n";
					//add
					//tag=0;
					/* iter_gene=hash_gene.find(vec_temp[i-1]);
					if(iter_gene!=hash_gene.end())
					{
						if(vec_temp[i]==iter_gene->second)//有cross-gene支持
						{
							tag=1;
						}
					} */
					strand=hash_strand[vec_flag[i-1]];
					vec_strand.clear();
					//if(tag==0)//没有支持才可断
					//{	
						//记录cluster的链方向
						for(int j=vec_flag[i-1]+1;j<vec_flag[i];j++)
						{
							iter_s=hash_strand.find(j);
							if(strand!=iter_s->second)//cluster的方向和它插入基因的方向不同
							{
								vec_strand.push_back(1);
							}	
						}
						if(vec_strand.size()>0)
						{
							vec_transcript.push_back(vec_strand.size());//在i基因和i-1之间插有反向基因
							vec_tag.push_back(i);
						}
					//}
				}
			}
			//vec_tag记录了该cluster在哪有插入反向基因情况
			//vec_transcript的size为有几处有反向基因转录情况，元素为在此处的反向基因个数
			vec_tag.push_back(vec_temp.size());
			if(vec_transcript.size()!=0)
			{
				for(int i=0;i<vec_tag.size()-1;i++)
				{
					vec_time.clear();
					for(int j=vec_tag[i];j<vec_tag[i+1];j++)
					{
						vec_time.push_back(vec_temp[j]);
					}
					vec_time.push_back(strand);
					Vec_output_2.push_back(vec_time);
					vec_time.clear();
				}
			
			}
				
			if(vec_transcript.size()==0)//不含反向基因插入情况的
			{
				vec_time.clear();
				for(int i=0;i<vec_temp.size();i++)
				{
					vec_time.push_back(vec_temp[i]);
				}
				vec_time.push_back(strand);
				Vec_output_2.push_back(vec_time);
				vec_time.clear();
			}
		}	
	}//while
	
	vec_flag.clear();
	vec_strand.clear();
	vec_tag.clear();
	vec_temp.clear();
	vec_transcript.clear();
	/* hash_gene.clear(); */
	return Vec_output_2;
}
	
vector<vector<string> > MATUcluster::average_select(vector<vector<string> >&Vec_Cluster_2)
{
	vector<vector<string> > Vec_input_3=Vec_Cluster_2;
	vector<vector<string> > Vec_output_3;
	vector<string> vec_time;//中间向量
	vector<int> vec_location;
	vector<int>::iterator iter_1;
	int max;
	int min;
	vector<int> vec_bp_part;
	double sum_add;
	double average_add;
	
	
	for(int i=0;i<Vec_input_3.size();i++)
	{
		vec_location.clear();
//cout<<"Vec_input_3[i].size():"<<Vec_input_3[i].size()<<"\n";
		for(int j=0;j<Vec_input_3[i].size()-1;j++)
		{
			vec_time.clear();
			vec_time=hash[Vec_input_3[i][j]];
			vec_location.push_back(atoi(vec_time[0].c_str()));
			vec_location.push_back(atoi(vec_time[1].c_str()));
			vec_time.clear();
		}
//cout<<"vec_time.size():"<<vec_time.size()<<"\n";
		iter_1=min_element(vec_location.begin(),vec_location.end());
		min=*iter_1;
		iter_1=max_element(vec_location.begin(),vec_location.end());
		max=*iter_1;
//cout<<"max:"<<max<<"\n";
//cout<<"min:"<<min<<"\n";
		
		vec_bp_part.clear();
		for(int k=min-1;k<max;k++)
		{
			vec_bp_part.push_back(vec_bp[k]);
		}
		sum_add=0;
		for(int k=0;k<vec_bp_part.size();k++)
		{
			sum_add=sum_add+vec_bp_part[k];
		}	
		average_add=sum_add/vec_bp_part.size();
		if(average_add>cutoff_2)
		{
			Vec_output_3.push_back(Vec_input_3[i]);
		}
	}
	
	return Vec_output_3;
}
		
		
			
	
	
#endif	
