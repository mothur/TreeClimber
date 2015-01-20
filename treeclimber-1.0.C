//  Parts of functions used to read in a tree were inspired from Bret Larget's BAMBE source code, which
//  can be found at http://www.mathcs.duq.edu/larget/bambe.html

#include <iostream>

using namespace std;


#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <map>

#define MAX_LINE		513
#define SKIPLINE(f,c)	{while((c=f.get())!=EOF && ((c) != '\n')){}}

/**************************************************************************************************/

void usageError(char *name)
{
	cerr << "\n*******************************************************************************\n";
	cerr << "Usage: " << name << " -t <treefile> -n <name file> [-i Iterations (<1000>)]" << endl;
	cerr << "\nOptions:\n";
	cerr << " YOU MUST EITHER SET BOTH -t AND -n OR -r!!!\n";
	cerr << " -t: Phylip or NEXUS file containing list of newick-formatted trees\n";
	cerr << " -n: Two column file where first column contains sequence names and second\n";
	cerr << "     contains a label for that sequences group.\n";
	cerr << " -r: Generate and score random trees without an input tree.\n";
	cerr << " -i: Number of iterations (default = 1000)\n";
	cerr << "*******************************************************************************\n\n";
	exit(1);
}

/**************************************************************************************************/

int readSpecialChar(istream& f, char c, const char *name)
{
	char d;

	while(isspace(d=f.get()))		{;}
	if(d == EOF){
		cerr << "Error: Input file ends prematurely, expecting a " << name << "\n";
		exit(1);
	}
	if(d != c){
		cerr << "Error: Expected " << name << " in input file.  Found " << d << ".\n";
		exit(1);
	}
	if(d == ')' && f.peek() == '\n'){
		while(isspace(d=f.get()))		{;}
		f.putback(d);
	}	
	

	return d;
}

/**************************************************************************************************/

int readNodeChar(istream& f)
{
	char d;
	while(isspace(d=f.get()))		{;}
	if(d == EOF){
		cerr << "Error: Input file ends prematurely, expecting a left parenthesis\n";
		exit(1);
	}
	return d;
}

/**************************************************************************************************/

double readBranchLength(istream& f)
{
	double b;
	
	if(!(f >> b)){
		cerr << "Error: Missing branch length in input tree.\n";
		exit(1);
	}
	return b;
}

/**************************************************************************************************/

int binarySearchString(vector<string> a, string value)
{
   	int midval;
	int left = 0;
	int right = a.size()-1;
	while(left <= right){
         midval = (left+right)/2;
         if(a[midval] == value){
             return midval;
		 }
         else if(value < a[midval]){
			 right = midval-1;
		 }
         else if(value > a[midval]){
			 left  = midval+1;
		 }
	}
	return -1;
}

/**************************************************************************************************/

template<typename T>
string tostring(const T&x){
	ostringstream out;
	out << x;
	return out.str();
}

/**************************************************************************************************/

class Node {
	public:
		Node(int i, string name, string group) :
			number(i),
			taxaname(name),
			lchild(NULL),
			rchild(NULL),
			parent(NULL),
			penalty(0)
		{	parsimony.push_back(group);}
		Node *lchild, *rchild;
		Node *parent;
		int number;
		string taxaname;
		double blength;
		vector<string> parsimony;
		int penalty;
		void setParsimony(Node *, Node *);
		void setChildren(Node *l, Node *r)	{	lchild = l; rchild = r;		}
		void setParent(Node *p)				{	parent = p;					}
		void setBranchLength(double b)		{	blength = b;				}
		int getNumber()						{ 	return number; 				}
		void setNodes(vector<Node*>& nodes) { 	nodes[number] = this; 		}
		void printnode(void);
};

class Tree {
	public:
		Tree(const vector<pair<string, string> >&, istream&, int&);
 		~Tree(){
			root->setNodes(nodes);
			for(int i=0;i<nodes.size();i++)
				delete nodes[i]; 
		}
	protected:
		int numleaves;
		int numnodes;
		int numtaxa;
		vector<string> names;
		void output_inputs(void);
		void read_tree(istream&);
		Node* readNewickInt(istream&, int&, vector<Node*>);
		int getParsimonyScore() {	return root->penalty;	}
		vector<Node*> nodes;
		Node *root;
		void printtree(void);
};

class RandomTree {
	public:
		RandomTree(const vector<pair<string, string> >&, int&);
 		~RandomTree(){
			root->setNodes(nodes);
			for(int i=0;i<nodes.size();i++)
				delete nodes[i]; 
		}
	protected:
		int numleaves;
		int numnodes;
		int numtaxa;
		vector<string> names;
		vector<Node*> nodes;
		Node *root;
		void createRandomTree(void);
		void printtree(void);
		int getParsimonyScore() {	return root->penalty;	}
};

class Output {
	public:
		Output(map<int,int>, map<int,int>, int, int, string);
	protected:
		map<int,int> userdata;
		map<int,int> randdata;
		int no_iters;
		int no_trees;
		double rand_cum;
		double user_cum;
		void printoutput(string);
		void printprobs(int, int, int, double&, double&, ofstream&);
		ofstream f;
};

/**************************************************************************************************/

Tree::Tree(const vector<pair<string,string> >& names_groups, istream& filename, int& score)
{
	numtaxa = names_groups.size();
	numleaves = numtaxa;
	numnodes = 2*numleaves-1;
	nodes.resize(numnodes);

	for(int i=0;i<numleaves;i++){
		names.push_back(names_groups[i].first);
		nodes[i] = new Node(i, names[i], names_groups[i].second);
	}
	for(int i=numleaves;i<numnodes;i++){
		nodes[i] = new Node(i, "", "");
	}
//	output_inputs();
	read_tree(filename);
//	printtree();
	score = getParsimonyScore();
}

/**************************************************************************************************/

void Tree::output_inputs(void)
{
	cout << "Names of taxa:\n";
	for(int i=0;i<numtaxa;i++){		cout << i << "\t" << names[i] << endl;	}
	cout << "\nNumber of leaves: " << numleaves << endl;
	cout << "Numbmer of nodes: " << numnodes << endl;
	cout << "Size of nodes vector: " << nodes.size() << endl << endl;
}

/**************************************************************************************************/

void Tree::read_tree(istream& tfile)
{
	int n = 0;
	Node *lc, *rc, *zc;
	int rooted = 0;

	int ch = tfile.peek();						
	if(ch == '('){
		n = numleaves;
		lc = readNewickInt(tfile,n,nodes);
		if(tfile.peek()==','){							
			readSpecialChar(tfile,',',"comma");
		}												
		else if((ch=tfile.peek())==';' || ch=='['){		
			rooted = 1;									
		}												
		if(rooted != 1){								
			rc = readNewickInt(tfile,n,nodes);
			if(tfile.peek() == ')'){					
				readSpecialChar(tfile,')',"right parenthesis");
			}											
		}												
	}
	else{
		tfile.putback(ch);
		char name[MAX_LINE];
		tfile.get(name, MAX_LINE,'\n');
		SKIPLINE(tfile,ch);
		n = binarySearchString(names,name);
		if(n!=0){
			cerr << "Internal error: The only taxon is not taxon 0.\n";
			exit(1);
		}
		lc = rc = NULL;
	}
	while((ch=tfile.get())!=';'){;}						

	if(rooted == 1){									
		root = (Node*)(nodes[numnodes-1]);				
	}													
	else{												
		root = (Node*)(nodes[n]);
		root -> setChildren(lc,rc);
		root -> setParsimony(lc,rc);
		root -> setBranchLength(0);
		root -> setParent(NULL);
		if(lc!=NULL){		lc -> setParent(root);		}
		if(rc!=NULL){		rc -> setParent(root);		}
	}													
}
	
/**************************************************************************************************/

Node* Tree::readNewickInt(istream& f, int& n, vector<Node*> nds)
{
	int c = readNodeChar(f);
	if(c == '('){
		Node* lc = readNewickInt(f, n, nds);
		readSpecialChar(f,',',"comma");
		Node* rc = readNewickInt(f, n, nds);

		if(f.peek()==')'){											
			readSpecialChar(f,')',"right parenthesis");					
		}															
		if(f.peek() == ':'){										
			readSpecialChar(f,':',"colon");								
			if(n >= numnodes){										
				cerr << "Error: Too many nodes in input tree\n";	
				exit(1);											
			}														
			nds[n] -> setBranchLength(readBranchLength(f));
		}															
		else{nds[n] -> setBranchLength(0.1);}						

		nds[n] -> setChildren(lc,rc);
		nds[n] -> setParsimony(lc,rc);
		lc -> setParent(nds[n]);
		rc -> setParent(nds[n]);
		return nds[n++];
	}
	else{
		f.putback(c);
		string name = "";
		char d=f.get();
		while(d != ':' && d != ',' && d!=')' && d!='\n'){						
			name += d;
			d=f.get();
		}
		int blen = 0;
		if(d == ':')	{		blen = 1;			}				
		f.putback(d);
		
		int n1 = binarySearchString(names, name);
		
		if(n1 == -1){cerr << "Name: " << name << " not found\n"; exit(1);}
	
		nds[n1] -> setChildren(NULL,NULL);

		if(blen == 1){															
			f.get();															
			nds[n] -> setBranchLength(readBranchLength(f));
		}																		

		while((c=f.get())!=0 && (c != ':' && c != ',' && c!=')') )		{;}		
		f.putback(c);															

		return nds[n1];
	}
}

/**************************************************************************************************/

void Tree::printtree(void)
{
//	cout << root->getNumber() << " " << numleaves << endl;
	cout << "taxon node lchild rchild parent parscr set\n";
	root->setNodes(nodes);
	for(int i=0;i<numnodes;i++){
		nodes[i]->printnode();
	}
}

/**************************************************************************************************/

void Node::printnode(void)
{
    cout	<< setw(4) << number+1 	<< ' ' << setw(5) << number << ' '
			<< setw(6) << (lchild?lchild->number:-1) << ' ' 
			<< setw(6) << (rchild?rchild->number:-1) << ' ' 
			<< setw(6) << (parent?parent->number:-1) << ' '
			<< setw(6) << penalty << ' ';
	cout << "[";
	if(parsimony.size()>0){
		for(int i=0;i<parsimony.size()-1;i++){
			cout << parsimony[i] << ",";
		}
		cout << parsimony[parsimony.size()-1] << "]\n";
	}
	else{
		cout << " ]\n";
	}
}

/**************************************************************************************************/

void Node::setParsimony(Node *lc, Node *rc)
{
	vector<string> merge;
	vector<string> v1 = lc->parsimony;
	vector<string> v2 = rc->parsimony;

	for(int i=0;i<v1.size();i++){
		if(v1[i] != "xxx"){
			merge.push_back(v1[i]);
		}
	}
	int size_v1=merge.size();
	
	for(int i=0;i<v2.size();i++){
		if(v2[i] != "xxx"){
			merge.push_back(v2[i]);
		}
	}
	int size_v2=merge.size()-size_v1;
	
	penalty = lc->penalty + rc->penalty;
	parsimony.pop_back();
	sort(merge.begin(),merge.end());
	
	if(size_v1>0 && size_v2>0){
	
		if(merge.size() > 0){
			for(int i=0;i<merge.size()-1;i++){
				if(merge[i] == merge[i+1]){
					parsimony.push_back(merge[i]);
					i++;
				}
			}
			if(parsimony.size() == 0){
				parsimony = merge;
				penalty++;
			}
		}
	}
	else if(size_v1==0 && size_v2==0){
		parsimony.push_back("xxx");
	}
	else if(size_v1==0 || size_v2==0){
		parsimony = merge;
	}
	
}	

/**************************************************************************************************/

RandomTree::RandomTree(const vector<pair<string,string> >& names_groups, int& score)
{
	numtaxa = names_groups.size();
	numleaves = numtaxa;
	numnodes = 2*numleaves-1;
	nodes.resize(numnodes);

	for(int i=0;i<numleaves;i++){
		names.push_back(names_groups[i].first);
		nodes[i] = new Node(i, names[i], names_groups[i].second);
	}
	for(int i=numleaves;i<numnodes;i++){
		nodes[i] = new Node(i, "", "");
	}
	
	createRandomTree();
//	printtree();

	score = getParsimonyScore();
	if(score == 8){
		cout << score << endl;
		printtree();
	}
}

/**************************************************************************************************/

void RandomTree::createRandomTree(void)
{
	for(int i=numleaves;i<numnodes;i++){
		int escape =0;
		int rnd_index1, rnd_index2;
		while(escape == 0){
			rnd_index1 = (int)(((double)rand() / (double) RAND_MAX)*i);
			if(nodes[rnd_index1]->parent == NULL){escape = 1;}
		}

		escape = 0;
		while(escape == 0){
			rnd_index2 = (int)(((double)rand() / (double) RAND_MAX)*i);
			if(rnd_index2 != rnd_index1 && nodes[rnd_index2]->parent == NULL){
				escape = 1;
			}		
		}

		nodes[i]->setChildren(nodes[rnd_index1],nodes[rnd_index2]);
		nodes[i]->setParent(NULL);
		nodes[rnd_index1]->setParent(nodes[i]);
		nodes[rnd_index2]->setParent(nodes[i]);
		nodes[i]->setParsimony(nodes[rnd_index1],nodes[rnd_index2]);
	}
	root = nodes[numnodes-1];	
}

/**************************************************************************************************/

void RandomTree::printtree(void)
{
//	cout << root->getNumber() << " " << numleaves << endl;
	cout << "taxon node lchild rchild parent parscr set\n";
	root->setNodes(nodes);
	for(int i=0;i<numnodes;i++){
		nodes[i]->printnode();
	}
}

/**************************************************************************************************/

void getNamesNGroups(string fname, vector< pair<string,string> >& namesngroups)
{
	ifstream f(fname.c_str());
	if(!f){
		cerr << "Error: We can't open the file containing the names and group assignments\n";
		exit(1);
	}

	while(!f.eof()){
		string name, group;
		f >> name;
		f >> group;
		namesngroups.push_back(make_pair(name, group));
	}
	if(namesngroups[namesngroups.size()-1].first == ""){
		namesngroups.pop_back();
	}
	sort(namesngroups.begin(),namesngroups.end());
}

/**************************************************************************************************/

Output::Output(map<int,int> user, map<int,int> rand, int trees, int iter, string outfile)
{
	userdata = user;
	randdata = rand;
	no_iters = iter;
	no_trees = trees;	

	printoutput(outfile);	
}

/**************************************************************************************************/

void Output::printoutput(string outfile)
{
	ofstream f(outfile.c_str(), ios::trunc);
	f.setf(ios::fixed, ios::floatfield);
	f.setf(ios::showpoint);
	
	double rand_cumprob = 0.0;
	double user_cumprob = 0.0;

	int minpars = userdata.begin()->first;
	if(minpars > randdata.begin()->first){
		minpars = randdata.begin()->first;
	}

	int maxpars = userdata.rbegin()->first;
	if(maxpars < randdata.rbegin()->first){
		maxpars = randdata.rbegin()->first;
	}

	cout << "Score\tUserFreq\tUserCumul\tRandFreq\tRandCumul\n";	
	f << "Score\tUserFreq\tUserCumul\tRandFreq\tRandCumul\n";	
	for(int i=minpars;i<=maxpars;i++){
		if(userdata[i] != 0 || randdata[i] != 0){
			printprobs(i, userdata[i], randdata[i], user_cumprob, rand_cumprob, f);
		}
	}	
}

/**************************************************************************************************/

void Output::printprobs(int index, int count1, int count2, double& cum1, double& cum2, ofstream& f)
{
	cum1 += (double)(count1)/(double)no_trees;
	cum2 += (double)(count2)/(double)no_iters;
	
	cout << setw(5) << index << "\t" << (double)(count1)/(double)no_trees << "\t" << cum1 << "\t";
	cout <<	(double)(count2)/(double)no_iters << "\t" << cum2 << endl;

	f << setw(5) << index << "\t" << (double)(count1)/(double)no_trees << "\t" << cum1 << "\t";
	f <<	(double)(count2)/(double)no_iters << "\t" << cum2 << endl;
}

/**************************************************************************************************/

void nexus(ifstream& f, vector< pair<string,string> >& namesgroups)
{
	string holder = "";
	int comment = 0;
	while(holder != "translate" && holder != "Translate"){	
	//	cout << holder << endl;
		if(holder == "[" || holder == "[!"){
			comment = 1;
		}
		if(holder == "]"){
			comment = 0;
		}
		f >> holder; 
		if(holder == "tree" && comment != 1){return;}
	}

	vector< pair<string,string> > translation;
	for(int i=0;i<namesgroups.size();i++){
		string number, name;
		f >> number;
		f >> name;
		name.erase(name.end()-1);
		translation.push_back(make_pair(name,number));
	}
	sort(translation.begin(),translation.end());

	for(int i=0;i<namesgroups.size();i++){
//		cout << i << " " << namesgroups[i].first << " " << namesgroups[i].second << " ";
//		cout << 			translation[i].first << " " << translation[i].second << endl;
		namesgroups[i].first = translation[i].second;
	}
	sort(namesgroups.begin(),namesgroups.end());
}

/**************************************************************************************************/

void getrandomdata(vector<pair<string,string> >& names_groups){
	
	int numGroups;
	int size;
	int index = 0;
	
	cout << endl;	
	cout << "Enter the number of groups:  ";
	cin >> numGroups;

	cout << "Enter the number of sequences in each group:  ";
	for(int i=0;i<numGroups;i++){
		cin >> size;
		for(int j=0;j<size;j++){
			names_groups.push_back(make_pair(tostring(index),tostring(i)));
		}
		index++;
	}
	cout << endl;
}

/**************************************************************************************************/
		
int main(int argc, char *argv[])
{
//	srand(54321);
	srand( (unsigned)time( NULL ) );

	cout.setf(ios::fixed, ios::floatfield);
	cout.setf(ios::showpoint);
	cerr.setf(ios::fixed, ios::floatfield);
	cerr.setf(ios::showpoint);
	
	string treefile = "random";
	string namesngroupfile = "";
	int iter = 1000;
	int random = 0;
	
	char **p;
	if(argc>1){	
		for(p=argv+1;p<argv+argc;p++){
			if(strcmp(*p, "-t")==0){
				if(++p>=argv+argc)		{		usageError(argv[0]);		}
				istringstream f(*p);
				if(!(f >> treefile))		{		usageError(argv[0]);		}
			}
			else if(strcmp(*p, "-n")==0){
				if(++p>=argv+argc) 		{		usageError(argv[0]);		}
				istringstream f(*p);
				if(!(f >> namesngroupfile))  	{		usageError(argv[0]);		}	
			}
			else if(strcmp(*p, "-i")==0){
				if(++p>=argv+argc) 		{		usageError(argv[0]);		}
				istringstream f(*p);
				if(!(f >> iter)) 		{		usageError(argv[0]);		}
			}
			else if(strcmp(*p,"-r")==0) {
				if(p>=argv+argc)		{		usageError(argv[0]);		}
			  random=1;
		  }
		}
	}
	else										{		usageError(argv[0]);		}

	if(random != 1 && (treefile == "" || namesngroupfile == ""))	{
		cerr << "Either the treefile or the file with the names and group id's for each\n";
		cerr << "sequence was not set.\n";
		exit(1);
	}

	vector< pair<string,string> > names_groups;
	map<int,int> score_freq;
	int count_trees;	

	if(random == 1){
		getrandomdata(names_groups);
		score_freq[0]++;
		count_trees = 1;
	}

	cout << "|---------------------------------------------------------------|\n";
	cout << "| % Progress      |   0  10  20  30  40  50  60  70  80  90 100 |\n";
	cout << "|---------------------------------------------------------------|\n";

	if(random != 1){
		getNamesNGroups(namesngroupfile, names_groups);
	
		ifstream f(treefile.c_str());
		if(!f){
			cerr << "Error: We can't open the input tree file\n";
			exit(1);
		}
	
		char c = f.get();
//		char nexus_format = 0;
		if(c == '#'){
			nexus(f, names_groups);
//			nexus_format = 1;
		}
		else{
			f.putback(c);
		}
		
		cout << "| User trees read |****";
		
		count_trees = 0;
		while(1){
			char d;
			while((d=f.get())!='('){
				if(d == '['){
					char e;
					while((e=f.get())!=']'){;}
				}
			}
			f.putback(d);
			cout.flush();
	
			int user_score;
			count_trees++;
			if((count_trees % 100) == 0){	cout << "****";		}
			
			Tree usertree(names_groups, f, user_score);
			score_freq[user_score]++;
//			usertree;
		
			while(isspace(d=f.get()))		{;}
			if(d!='(' && d!='t' && d!='['){ break;}
			f.putback(d);
		}
		cout << " " << count_trees << " trees read\n";
	}
	
	cout << "| Random trees    |";
	map<int,int> random_freq;
	int every = int((double)iter * 0.1);
	for(int i=0;i<iter;i++){

		int rand_score = 0;
		RandomTree rndtree(names_groups, rand_score);
		random_freq[rand_score]++;
//		rndtree;
		if(i % every == 0){			//	If iter is under 10, program with throw a
			cout << "****";			//	floating point exception
			cout.flush();
		}
	}
	cout << "**** |\n";
	cout << "|---------------------------------------------------------------|\n\n\n";

	cout.flush();
	string outfile = treefile + ".pars";	
	Output data(score_freq,random_freq,count_trees,iter,outfile);

	return 0;
}
