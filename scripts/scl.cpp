//run for 50kb            10 100 3
#include <sstream>
#include <cstdlib>
#include <iostream>
#include <string>
#include <map>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <time.h>
#include <sys/time.h> 
#include <unistd.h>
#include <iomanip>
#include <string.h>
//Debug
//#define USEDEBUG

#ifdef USEDEBUG
#define Debug(x) std::cout << x
#else
#define Debug(x) 
#endif 

//Debug2
//#define USEDEBUG2

#ifdef USEDEBUG2
#define Debug2(x) std::cout << x
#else
#define Debug2(x)
#endif



using namespace std;

struct node{
	double x;
        double y;
      	double z;
     	struct node *next;
     	struct node *prev;
}*start;

class Double_list{
	struct node{
        	double x;
        	double y;
        	double z;
        	struct node *next;
        	struct node *prev;
	}*start;

	public:
		int num_nodes;
	public:
		Double_list();
		void append_node(double x, double y, double z);
		double get_x_value(double position);
		double get_y_value(double position); //position is starting from 1, end at sequence length, so the first node has 1
		double get_z_value(double position);
		void set_value(double pos, double x, double y, double z);	//position starting from 1
		void clear();	//erase everything, free memory
};

Double_list::Double_list(){
	this -> num_nodes = 0;
	start = NULL;
}

void Double_list::append_node(double x, double y, double z){
	struct node *temp, *current;
	temp = new(struct node);
	temp -> x = x;
	temp -> y = y;
	temp -> z = z;
	temp -> next = NULL;
	if(start == NULL){	//empty list
		temp -> prev = NULL;
		start = temp;
	}
	else{
		current = start;
		while(current -> next != NULL){
			current = current -> next;
		}
		current -> next = temp;
		temp -> prev = current;
	}
	this -> num_nodes++;
}

double Double_list::get_x_value(double pos){
	if(start == NULL){	//that means it is empty
		cout << "Fatal Error: Cannot get the value of node, list is empty!\n";
		exit(1);
	}
	if(pos <= 0){
		cout << "Position of the list starting from 1. Invalid position number!\n";
		exit(1);
	}
	struct node *current;
	current = start;
	int i = 1;
	while(i < pos){
		current = current -> next;
                if(current == NULL){
                        cout << "This is beyond the number of node in the list! Cannot get value!\n";
                        exit(1);
                }
		i++;
	}
	return current -> x;
}

double Double_list::get_y_value(double pos){
        if(start == NULL){      //that means it is empty
                cout << "Cannot get the value of node, list is empty!\n";
                exit(1);
        }
        if(pos <= 0){
                cout << "Position of the list starting from 1. Invalid position number!\n";
	          exit(1);
        }
        struct node *current;
        current = start;
        int i = 1;
        while(i < pos){
                current = current -> next;
                if(current == NULL){
                        cout << "This is beyond the number of node in the list! Cannot get value!\n";
                        exit(1);
                }
		i++;
        }
        return current -> y;
}

double Double_list::get_z_value(double pos){
        if(start == NULL){      //that means it is empty
                cout << "Cannot get the value of node, list is empty!\n";
                exit(1);
        }
        if(pos <= 0){
                cout << "Position of the list starting from 1. Invalid position number!\n";
                exit(1);
        }
        struct node *current;
        current = start;
        int i = 1;
        while(i < pos){
                current = current -> next;
		if(current == NULL){
			cout << "This is beyond the number of node in the list! Cannot get value!\n";
			exit(1);
		}
		i++;
        }
        return current -> z;
}

void Double_list::set_value(double pos, double x, double y, double z){
        if(start == NULL){      //that means it is empty
                cout << "Cannot set the value of node, list is empty!\n";
                exit(1);
        }
        else{
                struct node *current;
                current = start;
                int i = 1;
                while(i < pos){
                        current = current -> next;
                        if(current == NULL){
                                cout << "This is beyond the number of node in the list! Cannot set value!\n";
                                exit(1);
                        }
                        i++;
                }
                current -> x = x;
		current -> y = y;
		current -> z = z;
        }
}

void Double_list::clear(){
	struct node *pDel = start;
	while(pDel != NULL){
		start = start -> next;
		delete pDel;
		pDel = start;
	}
	start = NULL;
}


//The following is a class that performs the simulation
//
//
class Simulation{
		int size_cube, length_seq,length_seq2;
		double delta0;
		double theta1;
		double mu1;
		double mu2;
		double rho;
		double beta;	
		double phi;
		double d0;
		double tau;
		double e1;
		//double open;
		double d1;
		double min_dist, max_dist;	//The min and max distance allowed, usually they are [2, sqrt(10)];
	public:
		double temperature;
		vector <vector <double> > real_dist;
		vector <vector <double> > delta; 
		Double_list list;	
		Simulation(int size, int length, double min, double max, vector < vector <double> > real_dist,double temper, vector < vector <double> > delta, double r, double the_1, double m_1, double m_2, double tau_1,double b,double del_0,double d_0,double p_1,/*double open_c,*/double d_1,int length2, vector < vector <double> > real_dist2, vector < vector <double> > delta2,double e1);
		//Simulation(int size, int length, double min, double max,vector < vector <double> > real_dist,double temper, vector <vector <double> > delta, double a, double b, double t_0, double d_0, double d_1);
		bool rand_initialize();	//initialize the sequence in the cube, randomly
		vector<node> find_neighbours_init(Double_list list);	//Find all the possible positions legal. This is for the initialization part,
									//which means the legal positions should be [2, sqrt root of 10] away from all the previous nodes
									//that have existent in the list, which were saved in the list
		void print_pdb_format();
		double calculate_cost1(double rand1,double zero,vector< vector <double> > delta,vector < vector <double> > real_dist,double length_seq);
		double calculate_cost2(double rand1);
		bool simulate_once(int lstart,int open);	//This function randomly select one node, move it, and make decision of acceptance or reject
		void print_list_to_file(char *file_name);	//Print final structure into file
		double En;
			
		Double_list list2;//a new list used to save the expand list
		bool expand_list(double res);//expand_list fuction
		double expand_length;
		vector <vector <double> > delta2;
		vector <vector <double> > real_dist2;
		vector<node> around_pos(int zero,int rand1,int length_seq,double d_1);//to find the candidates surrounding the randomly number
};

Simulation::Simulation(int size, int length, double min, double max, vector < vector <double> > real_dist,double temper, vector < vector <double> > delta, double r, double the_1, double m_1, double m_2, double tau_1,double b,double del_0,double d_0,double p_1,double d_1,int length2,vector < vector <double> > real_dist2, vector < vector <double> > delta2,double e1){
	this -> size_cube = size;
	this -> length_seq = length;
	this -> length_seq2 = length2;
	this -> min_dist = min;
	this -> max_dist = max;
	this -> real_dist = real_dist;
	this -> real_dist2 = real_dist2;
	this -> delta2 = delta2; 
	this -> delta = delta;
	this -> temperature = temper;
	this -> rho = r;
	this -> beta = b;
	this -> delta0 = del_0;
	this -> d0 = d_0;
	this -> mu1 = m_1;
	this -> mu2 = m_2;
	this -> tau = tau_1;
	this -> phi = p_1;
	this -> theta1 = the_1;
	this -> d1 = d_1;
	this -> e1 = e1;
}

//For the first node, randomly put it in the cube. For the others call function find_neigbbours_init
//Return true if initilization is successful, i.e., all the nodes have been put in tube, with
//enough space

bool Simulation::rand_initialize(){
	//debug3 cout << "Create a new list of nodes\n";
	Double_list list;
	this -> list = list;
	//cout << "Initialization: randomly put sequence in the cube:\n\n";
	for(int i = 0; i < length_seq; i++){
		//cout<<"	Node number: " << i + 1 << "\n";//@@@
		if(i == 0){	//This is the first node in the list or sequence
			timeval t1;
		        gettimeofday(&t1, NULL);
        		srand(t1.tv_usec * t1.tv_sec);
			//srand (time(NULL));
			int rand1 = rand() % size_cube; //The range of this random number is from 0 to size_cube - 1
                        int rand2 = rand() % size_cube;
                        int rand3 = rand() % size_cube;
			//cout<<"		node put at coordinates: " << rand1 << " " << rand2 << " " << rand3 << "\n";//@@@
			this -> list.append_node(rand1, rand2, rand3);
			////this -> list.append_node(0, 0, 0);	//This doesn't help. trying to maximize the radius of gyration(dont need this, bug sloved, BUG: the new node doesn't have
								//to be smaller than sqrt(10), but just > 2 with each of the nodes already in the list)
								//In this way, every random model will be starting from the origin. This may solve the problem 
								//that there are no enought space as a random first node may leave space not enough to have 
								//available candidate positions satisfying [2, sqrt(10)]
			double n = this -> list.get_x_value(1);
			//cout<<this -> list.get_x_value(1) << this -> list.get_y_value(1) << this -> list.get_z_value(1);//@@@
			//p
			//sleep(1);//@@@	
		}
		else{
			vector<node> candidates = find_neighbours_init(this -> list);
			if(candidates.size() == 0){
				Debug("The size of candidate nodes are 0! initialization failed!There are no enough space. Erase everything and restart\n");
				cout<<"The size of candidate nodes are 0! initialization failed! There are no enough space. Erase everything and restar"<<endl;
				this -> list.clear();
				//sleep(5);
				return false;
			}
			this -> list.append_node(candidates[0].x, candidates[0].y, candidates[0].z);		
			//cout<<candidates[0].x << " " << candidates[0].y << " " << candidates[0].z << " added to list\n";//@@@
			
		}
	}
	 //cout << "Initialization done! Initial nodes are:\n";
	//print_pdb_format();
	return true;
}


void Simulation::print_pdb_format(){
        for(int i = 0; i < length_seq; i++){
		//cout << "ATOM         N   MET A 108                                                   \n";
                 cout << "ATOM  " << right << setw(5) << "    " << " " << setw(4) << "CA" << " " << setw(3) << "MET" << " " << "A" << "        " << setw(8) <<this -> list.get_x_value(i + 1) << setw(8) << this -> list.get_y_value(i + 1) << setw(8) << this ->list.get_z_value(i + 1) << "\n";
		//cout << "ATOM         C   MET A 108                                                   \n";
		//cout << "ATOM         O   MET A 108                                                   \n";
		//cout << "ATOM         CB  MET A 108                                                   \n"; 
		//cout << "ATOM         CG  MET A 108                                                   \n";  
		//cout << "ATOM         SD  MET A 108                                                   \n";
		//cout << "ATOM         CE  MET A 108                                                   \n";
		//cout << "ATOM         H   MET A 108                                                   \n";
        }	
}



void Simulation::print_list_to_file(char *file_name){//use to print the first list
	ofstream file;
	file.open(file_name);
	if(file.is_open()){
		//cout<<"length of the list is" << length_seq << "\n";//@@@
		for(int i = 0; i < length_seq; i++){
			file << setw(8) << this -> list.get_x_value(i + 1) << setw(8) << this -> list.get_y_value(i + 1) << setw(8) << this -> list.get_z_value(i + 1) << "\n"; 
		}
		file.close();
	}
	else{
		cout
 << "Unable to open the output file! Print the coordinates onto screen!";
                for(int i = 0; i < length_seq; i++){
                        cout << setw(8) << this -> list.get_x_value(i + 1) << setw(8) << this -> list.get_y_value(i + 1) << setw(8) << this -> list.get_z_value(i + 1);
                }

	}
}

//Input is a linked list containing all the previously existent nodes. This function at first
//find all the nodes around the last node in the list, i.e., [x-sqrt(10), y - sqrt(10), z - sqrt(10)], x, y, z
//are the coordinates for the last node in the list, and iterate each one of them to 
//see whether it is [2, sqrt(10)] away from the existent nodes. At last, it returns all the legal candidate positions 
vector<node> Simulation::find_neighbours_init(Double_list list){
	vector<node> possible_pos;
	if(list.num_nodes == 0){
		cout << "The doubly linked list is empty, cannot find legal positions for next node!";
		exit(1);
	}
	//This newly added nodes must be within [2, sqrt(10)] of the last node in the list
	//so the new node must locate in [x - max - 1, x + max + 1], +-1 is because rounding, for
	//safe reason, because C truncate double when convert to int, so include more possibilities
	//So get all of the possible positions based on this criteria first, then graduately
	//eliminate unqualified ones
	double x_bottom_bound = ( (double)list.get_x_value(list.num_nodes) ) - (this -> max_dist) - 1;
	double x_upper_bound = ( (double)list.get_x_value(list.num_nodes) ) + this -> max_dist + 1;
	for(int x = (int)x_bottom_bound; x <= (int)x_upper_bound; x++){
		if(x < 0 || x >= this->size_cube){	//beyond the boundery of the cube, cube coordinates starting from 0
			continue;
		}
//		cout<<"x is " << x << "\n";//@@@
		double y_bottom_bound = ( (double)list.get_y_value(list.num_nodes) ) - this -> max_dist - 1;
		double y_upper_bound = ( (double)list.get_y_value(list.num_nodes) ) + this -> max_dist + 1;
		for(int y = (int)y_bottom_bound; y <= (int)y_upper_bound; y++){
			if(y < 0 || y >= this -> size_cube){
				continue;
			}
//			cout<<"y is " << y << "\n";//@@@
			double z_bottom_bound = ( (double)list.get_z_value(list.num_nodes) ) - this -> max_dist - 1;
			double z_upper_bound = ( (double)list.get_z_value(list.num_nodes) ) + this -> max_dist + 1;
			for(int z = (int)z_bottom_bound; z <= (int)z_upper_bound; z++){
				if(z < 0 || z >= this -> size_cube){
					continue;
				}
//				cout<<"z is " << z << "\n";//@@@
				struct node *temp;
			        temp = new(struct node);
			        temp -> x = x;
			        temp -> y = y;
			        temp -> z = z;
				//cout<<"		" << x << " " << y << " " << z << " " << "added as candidate positions\n";//@@@
				possible_pos.push_back(*temp);
				delete temp;
			}
		}
	}
//	cout<<"		number of candidates now: " << possible_pos.size() << "\n";//@@@

	//These possible positions should not contain any ones that are without [2, sqrt(10)] towards
	//to any of the pre-existed nodes saved in list. Now get ride of these ndoes.
	//cout<<"		Checking distance between candidate positions and each node in list\n";//@@@
	vector<node> possible_pos_2;	//Define another vector saving the legal positions
	bool random_find=false;
	int count=0;
	while(random_find==false){
		count++;
		if(count==possible_pos.size()){
			break;
		}
		timeval t1;
                gettimeofday(&t1, NULL);
                srand(t1.tv_usec * t1.tv_sec);
		
		int rand1= rand() % possible_pos.size()-1;	
		int x_p = possible_pos[rand1].x;
		int y_p = possible_pos[rand1].y;
		int z_p = possible_pos[rand1].z;
		//cout<<"Checking candidate position with each nodes already in the list..." << m << "\n";//@@@
		//cout<<x_p << " " << y_p << " " << z_p << "\n"; //@@@
		bool removed = false;
		for(int i = 1; i <= list.num_nodes; i++){	//List position starting from 1
			//cout<<"		checking the " << i << " node in the list already\n";	//@@@
			double x_e = list.get_x_value(i);
			double y_e = list.get_y_value(i);
			double z_e = list.get_z_value(i);
			//cout<<"		calculating the dist to node " << x_e << " " << y_e << " " << z_e << "\n";//@@@
			double dist = sqrt( (x_p - x_e) * (x_p - x_e) + (y_p - y_e) * (y_p - y_e) + (z_p - z_e) * (z_p - z_e) );
			//cout<<"		distance = " << dist << "\n";//@@@
			if(i < list.num_nodes){
				if(dist < this -> min_dist || dist == sqrt(8)){		//|| dist > this -> max_dist is removed, because it doesn't need to be smaller than sqrt(10) with all the other nodes, but only 
										//the last node in the list
					removed = true;
					//cout<<"		removed\n";//@@@
					break;
				}
			}
			if(i == list.num_nodes){	//This is judging with the last node in the list, this one can use dist > this -> max_dist
				if(dist < this -> min_dist || dist == sqrt(8) || dist > this -> max_dist){
					removed = true;
					//cout<<"         removed\n";//@@@
					break;
				}
			}
		}
		if(removed == false){
			possible_pos_2.push_back(possible_pos[rand1]);
			random_find=true;
		}
	}
	if(possible_pos_2.size() < 0){
		cout << "Fatal Error: the number of possible positions cannot be a negative number. Something is very wrong!\n";
		exit(1);
	}
	//cout<<"		size of candidates kept: " << possible_pos_2.size() << "\n";//@@@changed debug
	return possible_pos_2;
}


//Read from a file that contains the correlation of Hi-C, return a vector of vector

vector< vector<double> > read_correlation(double res,char* file){//read two position bead pairs to generate a hic matrix
        vector< vector<double> > vector2;
        vector< vector<double> > vector3;
        vector<double> row;
        string line;
        double temp;
        ifstream myfile;
        myfile.open(file);

        if (! myfile.is_open()) {
                cout<<"File :< " << file << " >open error"<< endl;
                exit(EXIT_FAILURE);
        }
        int big = 0;
        while (getline(myfile, line)){//save two columns position in array
                istringstream istr(line);
                while (istr >> temp) {
                        row.push_back(temp);
                        if(big<=temp){
                                big=temp;//find the biggest position int the hic contact, set it as the nodes number of the chromosome
                        }
                        //cout<<temp<<" ";
                }
                //cout<<endl;
                vector2.push_back(row);
                if(vector2.back().empty()){
                        cout<<"File :< " << file << " >has an empty row"<< endl;
                        exit(EXIT_FAILURE);
                }
                row.clear();
                istr.clear();
                line.clear();
        }
        //cout<<"biggest is "<<big<<endl;
        myfile.close();
        res = res*1000000;
        int big_num = big/res+1;
        //cout<<"big bead is "<<big_num<<endl;
        for(int i = 0; i<big_num;i++){//generatre a null matrix
                vector <double> t1;
                for(int j =0; j<big_num;j++){
                        int temp_num = 0;
                        t1.push_back(temp_num);
                        //cout<<temp_num<<" ";
                }
                vector3.push_back(t1);
                //cout<<endl;
        }
        for(int i=0;i<vector2.size();i++){//give the value[i,j] in matrix 1 as a bead pair
                //cout<<" i is "<<i<<endl;
                int i0 = vector2[i][0]/res+1;
                int i1 = vector2[i][1]/res+1;
                //cout<<i0<<" "<<i1<<endl;
                if(i0!=i1){//dont consider the situation bead pair are two same position
                        vector3[i0-1][i1-1]=1;
                        vector3[i1-1][i0-1]=1;
                       // cout<<vector3[i0-1][i1-1]<<endl;
                }
        }
        return vector3;

}


//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//This function calculate the delta matrix
vector< vector<double> > calculate_t1(double length_seq, vector < vector <double> > real_dist,double mu2,double d0){//to generate the delta matrix, d0 is useless here
        vector< vector<double> > t1;
        vector< vector<double> > t2;
        vector< vector<double> > t3;
        
        for(int i=0; i<length_seq;i++){//get all the real hic contact bead pairs
                for(int j=0; j<length_seq;j++){
                        if(real_dist[i][j]==1){
                                vector <double> row;
                                row.push_back(i);
                                row.push_back(j);
                                t2.push_back(row);
                        }
                }
        } 
        for(int i =0; i<length_seq; i++){//build half a marix to calculate the delta which not be 1, than calculate all the beads pairs equals to 1 in the matrix
                vector<double> row;
                for(int j = 0; j<length_seq; j++){
                    double temp1=0;
                    if(j>=i){
                        if(real_dist[i][j]==0 && i != j){
                                for(int k =0; k< t2.size(); k++){
                                        if(t2[k][0]==i && t2[k][1]==j){
                                                temp1=1;
                                        }else{
                                                double i1 = t2[k][0];
                                                double j1 = t2[k][1];
                                                double temp2= (i1-i)*(i1-i)+(j1-j)*(j1-j);
                                                temp2 = temp2/mu2;
                                                temp2 = exp(-temp2);
                                                temp1 = temp1 + temp2;
                                        }
                                }
                                if(temp1>1){
                                        temp1 = 1;
                                }
                        }else{
                                temp1=0;
                        }
                    }else{

                    }
                    row.push_back(temp1);
                }
                t3.push_back(row);
        }
	
        for(int i=0;i<length_seq;i++){//put the half matrix to a full one
                vector <double> row;
                for(int j=0;j<length_seq;j++){
                        double temp1=0;
                        if(i<=j){
                                temp1=t3[i][j];
                        }
                        if(i>j){
                                temp1=t3[j][i];
                        }
			row.push_back(temp1);
                }
                t1.push_back(row);
        }
        return t1;
}

//this function is used to calcualte loss_function                    
double Simulation::calculate_cost1(double rand1,double zero,vector< vector <double> > delta,vector< vector <double> > real_dist,double length_seq){
	double cost1 = 0.0;
	double cost1_1 = 0.0;
	double cost1_2 = 0.0;
	double cost1_3 = 0.0;
	for(int i = 0; i < length_seq; i++){
		if(rand1 != i + 1){//
	//		cout<<"rand1 is "<<rand1<<";delta0 is "<<delta0<<"; theta1 is "<<theta1<<"; mu1 is "<<mu1<<"; beta is "<<beta<<"; tau is  "<<tau<<"; phi is "<<phi<<"; rho is "<<rho<<";"<<endl;
			double x_i,y_i,z_i,x_j,y_j,z_j;
			if(zero==1){
				x_i = list.get_x_value(i + 1);
        	                y_i = list.get_y_value(i + 1);
                	        z_i = list.get_z_value(i + 1);
                        	x_j = list.get_x_value(rand1);
                        	y_j = list.get_y_value(rand1);
                        	z_j = list.get_z_value(rand1);
			}
			if(zero==2){
	                        x_i = list2.get_x_value(i + 1);
        	                y_i = list2.get_y_value(i + 1);
                	        z_i = list2.get_z_value(i + 1);
                        	x_j = list2.get_x_value(rand1);
                        	y_j = list2.get_y_value(rand1);
                        	z_j = list2.get_z_value(rand1);
			}
                        double dist = sqrt( (x_i - x_j) * (x_i - x_j) + (y_i - y_j) * (y_i - y_j) + (z_i - z_j) * (z_i - z_j) );
                        if(real_dist[rand1-1][i]==1 || delta[rand1-1][i]==1){
				cost1_1 = cost1_1 + ((dist-delta0)*(dist-delta0))/(delta0*delta0);
			}
			if(real_dist[rand1-1][i]==0 && delta[rand1-1][i]<1 && delta[rand1-1][i]>theta1){
                                double t1 = delta[rand1-1][i];
                                t1 = pow(t1,1.0/3);
                                t1 = delta0/t1;
                                cost1_2=cost1_2 + beta*(1- exp(-(dist-t1)*(dist-t1)/mu1));
                        }
                        if(real_dist[rand1-1][i]==0 && delta[rand1-1][i]<=theta1){
                                double t2 = 0;
                                t2 = pow(theta1,1.0/3);
                                t2 = delta0/t2;
                                cost1_3=cost1_3 + tau*(1- 1/(1+exp(-(dist-(t2-rho)))/phi)) ;
			}
		}
	}
	cost1=cost1_1+cost1_2+cost1_3;
	//cout<<"cost is "<<cost1<<endl;
	return cost1;
}
double Simulation::calculate_cost2(double rand1){
        double cost1 = 0.0;
        double cost1_1 = 0.0;
        double cost1_2 = 0.0;
        double cost1_3 = 0.0;
        for(int i = 0; i < length_seq; i++){
                
		for(int j = i+1; j<length_seq;j++){
                       // cout<<"rand1 is "<<rand1<<";delta0 is "<<delta0<<"; theta1 is "<<theta1<<"; mu1 is "<<mu1<<"; beta is "<<beta<<"; tau is  "<<tau<<"; phi is "<<phi<<"; rho is "<<rho<<";"<<endl;
                        double x_i = list.get_x_value(i + 1);
                        double y_i = list.get_y_value(i + 1);
                        double z_i = list.get_z_value(i + 1);
                        double x_j = list.get_x_value(j + 1);
                        double y_j = list.get_y_value(j + 1);
                        double z_j = list.get_z_value(j + 1);
                        double dist = sqrt( (x_i - x_j) * (x_i - x_j) + (y_i - y_j) * (y_i - y_j) + (z_i - z_j) * (z_i - z_j) );
                        if(real_dist[i][j]==1 || delta[i][j]==1){
                                cost1_1 = cost1_1 + ((dist-delta0)*(dist-delta0))/(delta0*delta0);
                        }
                        if(real_dist[i][j]==0 && delta[i][j]<1 && delta[i][j]>theta1){
                                double t1 = delta[i][j];
                                t1 = pow(t1,1.0/3);
                                t1 = delta0/t1;
                                cost1_2=cost1_2 + beta*(1- exp(-(dist-t1)*(dist-t1)/mu1));
                        }
                        if(real_dist[i][j]==0 && delta[i][j]<=theta1){
                                double t2 = 0;
                                t2 = pow(theta1,1.0/3);
                                t2 = delta0/t2;
                                cost1_3=cost1_3 + tau*(1- 1/(1+exp(-(dist-(t2-rho)))/phi)) ;
                        }
                }
        }
        cost1=cost1_1+cost1_2+cost1_3;
	//cout<<"cost is "<<cost1<<endl;
        return cost1;
}
vector<node> Simulation::around_pos(int zero,int rand1,int length_seq,double d_1){
	vector<node> candidates_pos_1;
	double x,y,z,xmin,xmax,ymin,ymax,zmin,zmax;
	if(zero==1){
		x = list.get_x_value(rand1);
		y = list.get_y_value(rand1);
		z = list.get_z_value(rand1);
	}
	if(zero==2){
		x = list2.get_x_value(rand1);
                y = list2.get_y_value(rand1);
                z = list2.get_z_value(rand1);
	}
        for(double i = x -1; i <= x + 1;i= i+1){
		if(i<0||i>=length_seq*5){continue;}
                for(double j = y - 1; j <= y + 1;j= j+1){
                        if(j<0||j>=length_seq*5){continue;}
			for(double k = z - 1; k <= z + 1;k= k+1){
				if(k<0||k>=length_seq*5){continue;}
                                if(x == i && y == j && z == k){
                                        continue;
                                }
                                else{
                                        struct node *temp;
                                        temp = new(struct node);
                                        temp -> x = i;
                                        temp -> y = j;
                                        temp -> z = k;
                                        if(rand1 > 1 && rand1 < length_seq ){
						if(zero==1){
					                xmin = this->list.get_x_value(rand1-1);
					                xmax = this->list.get_x_value(rand1+1);
					                ymin = this->list.get_y_value(rand1-1);
					                ymax = this->list.get_y_value(rand1+1);
					                zmin = this->list.get_z_value(rand1-1);
					                zmax = this->list.get_z_value(rand1+1);
						}
						if(zero==2){
							xmin = this->list2.get_x_value(rand1-1);
		                			xmax = this->list2.get_x_value(rand1+1);
                					ymin = this->list2.get_y_value(rand1-1);
        					        ymax = this->list2.get_y_value(rand1+1);
              						zmin = this->list2.get_z_value(rand1-1);
                					zmax = this->list2.get_z_value(rand1+1);
						}                                                  
                                                double dist1 = sqrt((i-xmin)*(i-xmin)+(j-ymin)*(j-ymin)+(k-zmin)*(k-zmin));
                                                double dist2 = sqrt((i-xmax)*(i-xmax)+(j-ymax)*(j-ymax)+(k-zmax)*(k-zmax));                                                
                                               
						if(dist1>(d_1/1) || dist2>(d_1/1) ){
                                                        continue;
                                                }
                                        }
                                        if(rand1  == 1){
						if(zero==1){
                                                        xmax = this->list.get_x_value(rand1+1);
                                                        ymax = this->list.get_y_value(rand1+1);
                                                        zmax = this->list.get_z_value(rand1+1);
						}if(zero==2){
                                                        xmax = this->list2.get_x_value(rand1+1);
                                                        ymax = this->list2.get_y_value(rand1+1);
                                                        zmax = this->list2.get_z_value(rand1+1);
						}
						double dist2 = sqrt((i-xmax)*(i-xmax)+(j-ymax)*(j-ymax)+(k-zmax)*(k-zmax));
                                                if(dist2 >(d_1/1)){
                                                        continue;
                                                }
                                        }
                                        if(rand1==length_seq){
						if(zero==2){
                                                        xmin = this->list2.get_x_value(rand1-1);
                                                        ymin = this->list2.get_y_value(rand1-1);
                                                        zmin = this->list2.get_z_value(rand1-1);
						}
						if(zero==1){
                                                        xmin = this->list.get_x_value(rand1-1);
                                                        ymin = this->list.get_y_value(rand1-1);
                                                        zmin = this->list.get_z_value(rand1-1);
						}
                                                double dist1 = sqrt((i-xmin)*(i-xmin)+(j-ymin)*(j-ymin)+(k-zmin)*(k-zmin));
						if(dist1 > (d_1/1)){
                                                        continue;
                                                }
                                        }
                                        candidates_pos_1.push_back(*temp);
                                        delete temp;
                                }
                        }
                }
        }
	return candidates_pos_1;
}
             

bool Simulation::simulate_once(int lstart,int open){
	
	timeval t1;
	gettimeofday(&t1, NULL);
	srand(t1.tv_usec * t1.tv_sec);
	int rand1 =0;
	if(open==2){
		rand1 = (rand() % this -> list2.num_nodes)+1;//random number is from 0 to list2.num_nodes - 1; we should plus them 1 to use for getting value from the list, the list idex from 1 to n;
	}
	if(open==0){
		rand1 = (rand() % this -> list.num_nodes) + 1; //The range of this random number is from 0+1 to num_nodes - 1 + 1
	}
	//cout <<" randomly choose the node at position: " << rand1 << " to move...\n";	
	//Get the number of possible candidate positions
	//put the 27 - 1 possible positions into a vector of nodes
	
	vector<node> candidates_pos;       
    if(open==2){
	candidates_pos = around_pos(2,rand1,length_seq2,d1);
	if(candidates_pos.size()==0){	
		return false;
	}
    }
    if(open==0){
	candidates_pos = around_pos(1,rand1,length_seq,d1);                                
        if(candidates_pos.size()==0){
                return false;
        }

    }
	//cout<<"now have "<<num_cand<<" candidates\n";//!!!
	//select the movement in one of the candidate positions, randomly
        gettimeofday(&t1, NULL);
        srand(t1.tv_usec * t1.tv_sec);
	int rand2 = rand() % candidates_pos.size();        //The range is from 0 to num_cand - 1
	//debug3 cout<<"the chosed node is "<<candidates_pos[rand2].x<<" "<<candidates_pos[rand2].y<<" "<<candidates_pos[rand2].z<<"\n";	
////////////////////////////////////////////////////////////////////////////////////////////////////////////
	double cost1=0.0;
	double cost2=0.0;
	if(open==0&&lstart!=1){
		cost1 = calculate_cost1(rand1,1,delta,real_dist,length_seq);
	}
	if(lstart==1&&open==0){
		cost1 = calculate_cost2(rand1);
	}
	if(open==2){
		cost1 = calculate_cost1(rand1,2,delta2,real_dist2,length_seq2);
	}
//	cout<<"cost_before is "<<cost1<<endl;
	double oldx,oldy,oldz;
	if(open==2){	
		oldx = this -> list2.get_x_value(rand1);
		oldy = this -> list2.get_y_value(rand1);
       	 	oldz = this -> list2.get_z_value(rand1);	
		//cout<<" selected node after is "<<candidates_pos[rand2].x<<" "<<candidates_pos[rand2].y<<" "<<candidates_pos[rand2].z<<endl;
		list2.set_value(rand1, candidates_pos[rand2].x, candidates_pos[rand2].y, candidates_pos[rand2].z);
	}
	if(open==0){
		oldx = this -> list.get_x_value(rand1);
        	oldy = this -> list.get_y_value(rand1);
        	oldz = this -> list.get_z_value(rand1);
		list.set_value(rand1, candidates_pos[rand2].x, candidates_pos[rand2].y, candidates_pos[rand2].z);
	}
	//double cost2 = calculate_cost1(rand1);	
	if(open==0&&lstart!=1){
                cost2 = calculate_cost1(rand1,1,delta,real_dist,length_seq);
        }
        if(open==0&&lstart==1){
                cost2 = calculate_cost2(rand1);
        }
	if(open==2){
		cost2 = calculate_cost1(rand1,2,delta2,real_dist2,length_seq2);
	}
	//cout<<"cost_after is "<<cost2<<endl;
	double Z;        
 	Z = cost2 - cost1;
	//cout<<"cost difference is "<<Z<<endl;
	//cout<<"temperature is "<<temperature<<endl;
	if(Z <= 0){
		//if(open==2){	
		//	cout << "cost difference less than 0, accept !!!!!!!!!!!!!!!!!!!!"<<endl;
		//}
		En = cost2;
		return true;
	} 
	else{
		double rand3=0.0;
                gettimeofday(&t1, NULL);
                srand(t1.tv_usec * t1.tv_sec);
		//cout<<"cost difference is bigger than 1"<<endl;
		//cout<<"temperature is "<<temperature<<endl;
        	//cout<<"exp(-Z/temperature) is "<<exp(-Z/temperature)<<endl;
                rand3 = rand()%(999999+1)/(float)(999999+1);
		//cout<<"rand3 is "<<rand3<<endl;
		if(rand3 < exp( -Z/temperature) ){
			//if(open==2){
			//	cout<< "random number is "<<rand3<<", less  than key,  accept 22222222222"<<endl;
			//}
			En = cost2;
			return true;
		}else{
			if(open==2){
				list2.set_value(rand1, oldx, oldy, oldz);
			}
			if(open==0){
				list.set_value(rand1,oldx,oldy,oldz);
			}
                	En = cost1;
                	//if(open==2){
			//	cout <<"random number is "<<rand3<<" bigger than key, Rejected!!!!\n";
			//}
			return false;
		}	
	}
	
}
bool Simulation::expand_list(double res){
	Double_list list2;
	//cout<<"get into function"<<endl;
	this -> list2 = list2;
	res = res*1000000;
	//cout<<"res is "<<res<<endl;
	int bad_num=0;
	int g = (int)(0.5*1000000)/res;
	int g1=g;
	//cout<<"g is "<<g<<"; should be 10"<<endl;
	double ratio = e1*g;
	for(int i = 1; i<=length_seq;i++){
		if(i<length_seq){
			if(i==length_seq-1){
				g1=length_seq2-(i-1)*g-1;
			}
			double add_pos=0;
			int num=0;//the numbers already pushed between two beads 
			double r = ratio;
			double x_2=list.get_x_value(i+1)*r;
			double y_2=list.get_y_value(i+1)*r;
			double z_2=list.get_z_value(i+1)*r;
			double x_1=list.get_x_value(i)*r;
			double y_1=list.get_y_value(i)*r;
			double z_1=list.get_z_value(i)*r;
			this->list2.append_node(x_1,y_1,z_1);
			//cout<<"x2 y2 z2 is "<<x_2<<" "<<y_2<<" "<<z_2<<endl;
 //                       cout<<" x1 y1 z1 is "<<x_1<<" "<<y_1<<" "<<z_1<<endl;
//			cout<<(i-1)*10+1<<" node is done!!"<<endl;
		    for(int l=(i-1)*g+1;l<(i-1)*g+g1;l++){
			//cout<<l<<"node is "<<endl;
			num++;
			//cout<<" this is "<<num<<"node in "<<g<<endl;
			//cout<<this->max_dist<<" is the max"<<endl;
			//cout<<"now has "<<cout1<<" nodes in the list2"<<endl;	
			vector <node> possible_pos;
			double x_bottom_bound = ( (double)this->list2.get_x_value(l) ) - sqrt(10) - 1;
		        double x_upper_bound = ( (double)this->list2.get_x_value(l) ) + sqrt(10) + 1;
        		for(int x = (int)x_bottom_bound; x <= (int)x_upper_bound; x++){
				if(x<0){continue;}
				double y_bottom_bound = ( (double)this->list2.get_y_value(l) ) - sqrt(10) - 1;
                		double y_upper_bound = ( (double)this->list2.get_y_value(l) ) + sqrt(10) + 1;
                		for(int y = (int)y_bottom_bound; y <= (int)y_upper_bound; y++){
					if(y<0){continue;}
					double z_bottom_bound = ( (double)this->list2.get_z_value(l) ) - sqrt(10)  - 1;
                        		double z_upper_bound = ( (double)this->list2.get_z_value(l) ) + sqrt(10) + 1;
					for(int z = (int)z_bottom_bound; z <= (int)z_upper_bound; z++){
						if(z<0){continue;}
						//cout<<"x y z is"<<x<<" "<<y<<" "<<z<<endl;
						//cout<<"x1 y1 z1 is "<<x_1<<" "<<y_1<<" "<<z_1<<endl;
						//cout<<"x2 "
						struct node *temp;
                                		temp = new(struct node);
                                		temp -> x = x;
                                		temp -> y = y;
                                		temp -> z = z;
						double dist = sqrt( (x_2-x_1)*(x_2-x_1)+(y_2-y_1)*(y_2-y_1)+(z_2-z_1)*(z_2-z_1));
						double dist2 = sqrt( (x_2-x)*(x_2-x)+(y_2-y)*(y_2-y)+(z_2-z)*(z_2-z));
						//cout<<"new dist is "<<dist2<<"; old dist is "<<dist<<endl;
						if(dist2<dist){
							possible_pos.push_back(*temp);
						}
		                                delete temp;
					}
				}
			}
			//cout<<"size of possible_pos is"<<possible_pos.size()<<endl;
			//this->list2.append_node(x_1,y_1,z_1);
			if(possible_pos.size()==0){
				//cout<<"no ones dist less than before"<<endl;
				add_pos=1;
			}
			vector<node> possible_pos_2;
			bool random_find=false;
			int count=0;
        		while(random_find==false){
				count++;
				double r = ratio;
				if(add_pos==1){
                                        struct node *temp;
                                        temp = new(struct node);
                                        double xa = this->list2.get_x_value(l);
                                        double ya = this->list2.get_y_value(l);
                                        double za = this->list2.get_z_value(l);
                                        temp -> x = xa+(x_2 - xa)/(g1-num+1);
                                        temp -> y = ya+(y_2 - ya)/(g1-num+1);
                                        temp -> z = za+(z_2 - za)/(g1-num+1);
                                        possible_pos_2.push_back(*temp);
          //                              cout<<"push "<<possible_pos_2[0].x<<" "<<possible_pos_2[0].y<<" "<<possible_pos_2[0].z<<endl;
                                        delete temp;
                                        this->list2.append_node(possible_pos_2[0].x,possible_pos_2[0].y,possible_pos_2[0].z);
                                        bad_num++;
                                        //cout<<l+1<<" nodes is a bad one, bad nodes are "<<bad_num<<endl;
					add_pos=0;
                                        break;
				}
                		timeval t1;
                		gettimeofday(&t1, NULL);
                		srand(t1.tv_usec * t1.tv_sec);

                		int rand1= rand() % possible_pos.size();
                		int x_p = possible_pos[rand1].x;
                		int y_p = possible_pos[rand1].y;
                		int z_p = possible_pos[rand1].z;
	                	bool removed = false;
				//cout<<"here is good"<<endl;
                		for(int i = 1; i <= list.num_nodes; i++){
					double x_e = list.get_x_value(i)*r;
                        		double y_e = list.get_y_value(i)*r;
                        		double z_e = list.get_z_value(i)*r;
					double dist1 = sqrt( (x_p - x_e) * (x_p - x_e) + (y_p - y_e) * (y_p - y_e) + (z_p - z_e) * (z_p - z_e) );
					if(dist1<2||dist1==sqrt(8) ){
						removed = true;
						break;
					}
				}
				for(int i = 1; i <= this->list2.num_nodes; i++){
                                    if(i!=this->list2.num_nodes){ 	
				        double x_e = this->list2.get_x_value(i);
                                        double y_e = this->list2.get_y_value(i);
                                        double z_e = this->list2.get_z_value(i);
                                        double dist1 = sqrt( (x_p - x_e) * (x_p - x_e) + (y_p - y_e) * (y_p - y_e) + (z_p - z_e) * (z_p - z_e) );
                                        if(dist1<2||dist1==sqrt(8) ){
                                                removed = true;
                                                break;
                                        }
				    }
				    if(i==this->list2.num_nodes){
                                        double x_e = this->list2.get_x_value(i);
                                        double y_e = this->list2.get_y_value(i);
                                        double z_e = this->list2.get_z_value(i);
                                        double dist1 = sqrt( (x_p - x_e) * (x_p - x_e) + (y_p - y_e) * (y_p - y_e) + (z_p - z_e) * (z_p - z_e) );
                                        if(dist1<2||dist1==sqrt(8)||dist1>sqrt(10) ){
                                                removed = true;
                                                break;
                                        }
                                    }
                                }
				if(removed==false){
					possible_pos_2.push_back(possible_pos[rand1]);
					//cout<<"x1 y1 z1 is"<<possible_pos_2[0].x<<" "<<possible_pos_2[0].y<<" "<<possible_pos_2[0].z<<endl;
					this->list2.append_node(possible_pos_2[0].x,possible_pos_2[0].y,possible_pos_2[0].z);
//					cout<<"x~y~z~ is "<<this->list2.get_x_value(l+1)<<" "<<this->list2.get_y_value(l+1)<<" "<<this->list2.get_z_value(l+1)<<endl;
//					cout<<l+1<<" node is done !!!!!!!!!!!!!!!"<<endl;
					random_find=true;
				}
				if(count==possible_pos.size()*10&&removed==true){
					struct node *temp;
                                	temp = new(struct node);
	                                double xa = this->list2.get_x_value(l);
        	                        double ya = this->list2.get_y_value(l);
                	                double za = this->list2.get_z_value(l);
                        	        temp -> x = xa+(x_2 - xa)/(g1-num+1);
                      		        temp -> y = ya+(y_2 - ya)/(g1-num+1);
	                                temp -> z = za+(z_2 - za)/(g1-num+1);
	                                possible_pos_2.push_back(*temp);
//	                                cout<<"push "<<possible_pos_2[0].x<<" "<<possible_pos_2[0].y<<" "<<possible_pos_2[0].z<<endl;
        	                        delete temp;
                	                this->list2.append_node(possible_pos_2[0].x,possible_pos_2[0].y,possible_pos_2[0].z);
                        	        bad_num++;
                               		break;
				}
			}
		    }	
		}if(i==length_seq){
			double x=this->list.get_x_value(i);
                        double y=this->list.get_y_value(i);
                        double z=this->list.get_z_value(i);
                        //this->list2.append_node(x,y,z);
			double xx_1=x*ratio;
                        double yy_1=y*ratio;
                        double zz_1=z*ratio;
//                        cout<<length_seq2<<"th done x* y* z* "<<xx_1<<" "<<yy_1<<" "<<zz_1<<endl;  
                        this->list2.append_node(xx_1,yy_1,zz_1);
			//cout<<"x* y* z* "<<this->list2.get_x_value(cou)/10<<" "<<this->list2.get_y_value(cou)/10<<" "<<this->list2.get_z_value(cou)/10<<endl;
		}	
	}
        cout<<bad_num<<" out of "<<length_seq2<<"segments cannot be expanded using the preferred approach due to no enough space for newly-added high-resolution beads. Usually this is OK; and you do not need to anything. But if you want, you can slightly increase the expansion index by -e parameter. The current expansion index is "<<e1<<"."<<endl;	
	return true;
}


int main(int argc, char* argv[]){

	int seq_length=0;
	double temperature = 10.0;
	string dist_fi = "none";
	string dist_fi2="none";
	string out_f = "none1";
	string out_f2 = "none2";
	string out_f3 = "none3";
	double rho = 1.0;
	double phi = 0.1;
	double theta1=0.7;
	double beta = 1.0;
	double delta0 = 8.0;
	double mu1 = 20.0;
	double mu2 = 2.0;
	double tau=1.0;
	int d0 = 5;
	int open = 0;
	double d1 = 8.0;
	int seq_length2=0;
	double res = 0.5;
	double temperature2=0.1;
	double e1=0.12;	
	for(int i = 1; i < argc - 1; i++){
                if(strcmp(argv[i],"-i")==0){
                        dist_fi = argv[i+1];
                }else if(strcmp(argv[i],"-i2")==0){
			dist_fi2 = argv[i+1];
		}else if(strcmp(argv[i],"-o")==0){
                        out_f = argv[i+1];
                }else if(strcmp(argv[i],"-o2")==0){
			out_f2 = argv[i+1];
		}else if(strcmp(argv[i],"-o3")==0){
			out_f3 = argv[i+1];
		}else if(strcmp(argv[i],"-l")==0){
                        string length (argv[i+1]);
                        seq_length = atoi(length.c_str());
                }else if(strcmp(argv[i],"-t")==0){
                        string temper (argv[i+1]);
                        temperature = atof(temper.c_str());
                }else if(strcmp(argv[i],"-rho")==0){
			string r (argv[i+1]);
			rho = atof(r.c_str()); 
		}else if(strcmp(argv[i],"-theta1")==0){
                        string r (argv[i+1]);
                        theta1 = atof(r.c_str()); 
                }else if(strcmp(argv[i],"-delta0")==0){
                        string r (argv[i+1]);
                        delta0 = atof(r.c_str()); 
                }else if(strcmp(argv[i],"-mu1")==0){
                        string r (argv[i+1]);
                        mu1 = atof(r.c_str()); 
                }else if(strcmp(argv[i],"-beta")==0){
                        string r (argv[i+1]);
                        beta = atof(r.c_str()); 
                }else if(strcmp(argv[i],"-tau")==0){
                        string r (argv[i+1]);
                        tau = atof(r.c_str()); 
                }else if(strcmp(argv[i],"-phi")==0){
                        string r (argv[i+1]);
                        phi = atof(r.c_str()); 
                }else if(strcmp(argv[i],"-times")==0){//d0 now is the argument to modified the simulation times
                        string r (argv[i+1]);
                        d0 = atoi(r.c_str()); 
                }else if(strcmp(argv[i],"-mu2")==0){
                        string r (argv[i+1]);
                        mu2 = atof(r.c_str()); 
                }else if(strcmp(argv[i],"-cost")==0){
			string r (argv[i+1]);
			open = atoi(r.c_str());
		}else if(strcmp(argv[i],"-d1")==0){
			string r (argv[i+1]);
			d1 = atof(r.c_str());
		}else if(strcmp(argv[i],"-res")==0){
			string r (argv[i+1]);
			res = atof(r.c_str());
		}else if(strcmp(argv[i],"-e")==0){
			string r (argv[i+1]);
			e1 = atof(r.c_str());
		}else if(strcmp(argv[i],"-t2")==0){//for the expand list simulation start from
			string r (argv[i+1]);
			temperature2 = atof(r.c_str());
		}	
        }
	char *dist_file=const_cast<char*>(dist_fi.c_str());
	if(dist_fi=="none"){
                cout<<"you should insert a hic contacts file"<<endl;
                exit(1);
        }
	vector < vector <double> > real_dist;//to save the hic matrix
	vector < vector <double> > real_dist2;//to save the high resolution hic matrix
	vector < vector <double> > delta;//to save the delta matrix
	vector < vector <double> > delta2;//to save the high resolution delta matrix
	if(res>=0.5){
		real_dist = read_correlation(res,dist_file);
		seq_length=real_dist.size();	
		cout<<"Target resolution: "<<res<<" Mb. A structure with "<<seq_length<<" beads will be generated."<<endl;
		cout<<"Reading Hi-C data and generating the delta matrix."<<endl;
		delta = calculate_t1(seq_length,real_dist,mu2,d0);//d0 is useless
		cout<<"Generating "<<res*1000<<" kb structure"<<endl;
	}
	if(res<0.5){
		real_dist = read_correlation(0.5,dist_file);
		seq_length=real_dist.size();
		real_dist2 = read_correlation(res,dist_file);
		seq_length2=real_dist2.size();
		cout<<"Target resolution: "<<res<<" Mb. A structure at 500 kb resolution with "<<seq_length<<" beads will be generated first and then extended to "<<seq_length2<<" beads."<<endl;
		//cout<<"number of beads for matrix 1 is "<<seq_length<<"; 2 is "<<seq_length2<<endl;
		cout<<"Reading Hi-C data and generating the delta matrix."<<endl;
		delta = calculate_t1(seq_length,real_dist,mu2,d0);//d0 is useless
		delta2 = calculate_t1(seq_length2,real_dist2,mu2,d0);
		cout<<"Generating 500 kb structure."<<endl;
	}
	int size_cube = seq_length*5;//set the size cube as number of beads times 5 to let the structure moves sufficiently
        int sim_times = seq_length*100; 
	char *out_file=const_cast<char*>(out_f.c_str());
	char *out_file2=const_cast<char*>(out_f2.c_str());
	char *out_file3=const_cast<char*>(out_f3.c_str());

       
/* 	//print out the delta matrix for large one
      	string s= "delta_matrix_out";//store the delata matrix
        char *p= const_cast<char*>(s.c_str());
        ofstream file1;
        file1.open(p);
	if(!file1.is_open()){
                cout << "Unable to open the output the delta matrix! ";
                exit(1);
        }
	for(int i = 0; i< seq_length2;i++){
		for(int j=0; j< seq_length2;j++){
			file1<<delta[i][j]<<" ";
		}
		file1<<endl;
	}
	file1.close();
*/
	Simulation simulation(size_cube, seq_length, 2, sqrt(10), real_dist, temperature, delta,rho,theta1,mu1,mu2,tau,beta,delta0,d0,phi,d1,seq_length2,real_dist2,delta2,e1);
	cout<<"Initialization of the 3D structure started."<<endl;	
	bool initial_status = false;
	while(initial_status == false){
		initial_status = simulation.rand_initialize();
	}
	int n;
	vector<double> tem;
	vector<double> energy;
	vector<int> num;
	n = seq_length;
	double time;
	int failed_times = 0;
	cout<<"Cooling started."<<endl;
	for(int i =  1; i <= 1000000; i++ ){
		int counter = 0;        		
		int accept_counter = 0;
		cout<<"Temperature: "<<temperature<<endl;
		while(counter < sim_times){
			counter++;
			if(simulation.simulate_once(open,0) == true){
				accept_counter ++;	
			}
			energy.push_back(simulation.En);
			num.push_back(i);
			if(accept_counter== seq_length*10){
				failed_times = 0;
				break;
			}
			if(counter == sim_times - 1){
				failed_times++;
			}
		}
	//	cout<<"accept_counter is "<<accept_counter<<endl;
		tem.push_back(temperature);
		if(open==1){
			cout<<"Temperature: "<<temperature<<"; Cost: "<<simulation.En<<endl;
		}
		temperature = temperature * 0.9;
		simulation.temperature = temperature;
		if(failed_times == 3){
			break;
		}	
	}
	if(res>=0.5){	                        			
		simulation.print_list_to_file(out_file);
		cout<<res*1000<<" kb structure has been generated"<<endl;
	}
	if(res<0.5){
		//simulation.print_list_to_file(out_file2);
		cout<<"A structure at 500 kb resolution has been generated. Expanding the beads for a high-resolution structure."<<endl;
		bool expand_status = false;
        	if(expand_status == false){
                	expand_status == simulation.expand_list(res);
        	}
		cout<<"Fine-tuning started."<<endl;
        	temperature=temperature2;         
        	simulation.temperature=temperature;
                int counter=0;
                cout<<"Temperature:"<<temperature<<endl;
                while(counter<seq_length2*d0){
                	counter++;
                        if(simulation.simulate_once(2,2)==true){
                        }
                }
		cout<<"Fine-tuning is done. High-resolution structure has been successfully generated."<<endl;
        	ofstream file11;
        	file11.open(out_file);
       		if(!file11.is_open()){
                	cout << "Unable to open the output the cost file ";
                	exit(1);
        	}
		for(int i = 0 ; i <simulation.list2.num_nodes; i = i+1){
                	file11<<i<<"        "<<simulation.list2.get_x_value(i+1)/1<<" "<<simulation.list2.get_y_value(i+1)/1<<" "<<simulation.list2.get_z_value(i+1)/1<<endl;
		}
        	file11.close();
	}
	return 0;

}

