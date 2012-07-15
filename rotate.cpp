//rotate.cpp
//by Jeff Sung
//Last updated 2011-10-20 at 11:26
//rotate.cpp takes a PDB, reads in selections of atom numbers to rotate by a
//specific angle IN DEGREES, and then rotates those atoms, saving in a new pdb.
//The code below contains a simplistic (read: lazy) method of parsing PDBs.
//
//Syntax for file 'input'
//Odd Numbered lines: list of atom numbers to be moved (indexing starts from 1)
//Even Numbered lines: axis of rotation ('x', 'y', or 'z') & angle of rotation
//Will rotate arbitrary numbers of atom selections
//
//Notes: no error checking for existence or syntax of file 'input'
//       input and output are PDB files; no error checking for correct syntax
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <math.h>

using namespace std;

struct atom{
	int serial;
	int resSeq;
	double x, y, z;
	string before;
	string after;
	
	friend ostream &operator<<(ostream &os, const atom &at);
	friend bool operator==(atom &a1, atom &a2);
};

struct residue{
	vector <atom> atoms;
	double theta; //angle
	char axis; // possible values: 'x', 'y', or 'z'
};

double str2double(const string &str);
int str2int(const string & str);
atom readatom(const string & str);
string trim(const string & str);
vector <atom> split(const string & str);


double str2double(const string &str){
	stringstream ss(str);
	double d;
	if (!(ss>>d)) d=0;
	return d;
}
int str2int(const string &str){
	stringstream ss(str);
	int i;
	if (!(ss>>i)) i=0;
	return i;
}
atom readatom(const string &str){
	atom at;
	at.before=str.substr(0,30);
	at.after=str.substr(54,26);
	at.serial=str2int(str.substr(6,5));
	at.resSeq=str2int(str.substr(22,4));
	at.x=str2double(str.substr(30,8));
	at.y=str2double(str.substr(38,8));
	at.z=str2double(str.substr(46,8));
	return at;
}

string trim(const string & strorig){
	return strorig;
}

ostream &operator<<(ostream &os, const atom &at){
	os<<at.before;
	os<<setiosflags(ios_base::right|ios_base::fixed);
	os<<setw(8)<<setprecision(3)<<at.x;
	os<<setw(8)<<setprecision(3)<<at.y;
	os<<setw(8)<<setprecision(3)<<at.z;
	os<<at.after;
	return os;
}

bool operator==(atom &a1, atom &a2){
	return a1.serial==a2.serial;
}

vector <atom> split(const string & str){
	vector <atom> ret(1);
	atom at;
	stringstream ss(str);
	string s2;
	int i=0;
	for (; ss>>s2;){
		at.serial=str2int(s2);
		ret[i++]=at;
		if (i==ret.size()) ret.resize(2*i);
	}
	ret.resize(i);
	return ret;
}

int main(int argc,char *argv[]){
	if (argc < 3){
		cout<<"Syntax: "<<argv[0]<<" pdb_in pdb_out"<<endl;
	}
	else{
		int which=0, i, j=0;
		atom curat;
		residue curesid;
		string line, s2;
		double cur_x, cur_y, cur_z;

		//Read in the atoms that need to be moved and the requisite vectors
		ifstream infile("input");
		vector <residue> movelist(1);
		while (getline(infile,line)){
			movelist[which].atoms=split(line);
			infile>>movelist[which].axis>>movelist[which].theta;
			movelist[which].theta*=3.141592/180.0;
//			cout<<movelist[which].axis<<" "<<movelist[which].theta<<endl;	
			which++;
			infile.ignore(100,'\n');
			if (movelist.size()==which) movelist.resize(2*which);

		}	
		movelist.resize(which);
		infile.close();
		
		//Read in the given PDB, rotate any necessary atoms, and then write the resultant PDB
		infile.open(argv[1]);
		ofstream outfile(argv[2]);
		while (getline(infile,line)){
			if (line.find("ATOM")==string::npos&&line.find("HETATM")==string::npos){
				outfile<<line<<endl;
				continue;
			}
			curat=readatom(line);
			i=0;
			for (i=0; i<movelist.size(); i++){
				for (j=0; j<movelist[i].atoms.size(); j++){
					if (movelist[i].atoms[j]==curat){
							
						if ( movelist[i].axis=='x'){
							cur_x=curat.x;
							cur_y=curat.y*cos(movelist[i].theta)-curat.z*sin(movelist[i].theta);
							cur_z=-curat.y*sin(movelist[i].theta)+curat.z*cos(movelist[i].theta);
							curat.x=cur_x;
							curat.y=cur_y;
							curat.z=cur_z;
						}
						else if (movelist[i].axis=='y'){
							cur_x=curat.x*cos(movelist[i].theta)+curat.z*sin(movelist[i].theta);
							cur_y=curat.y;
							cur_z=-curat.x*sin(movelist[i].theta)+curat.z*cos(movelist[i].theta);
							curat.x=cur_x;
							curat.y=cur_y;
							curat.z=cur_z;
						}
						else if (movelist[i].axis=='z'){
							cur_x=curat.x*cos(movelist[i].theta)-curat.y*sin(movelist[i].theta);
							cur_y=curat.x*sin(movelist[i].theta)+curat.y*cos(movelist[i].theta);
							cur_z=curat.z;
							curat.x=cur_x;
							curat.y=cur_y;
							curat.z=cur_z;
						}

					}
				}
			}
			//Write
			outfile<<curat<<endl;	
		}
		infile.close();
		outfile.close();

	}
	return 0;
}
