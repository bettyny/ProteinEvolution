using boost::mt19937;
using boost::uniform_01;

struct Indiv_struct {
	double ind_score; // this AA fitness score
	char AA;  //what the current AA is
	std::string codon; //what the current codon is
};


std::pair<std::string,std::string> codontoAA(std::string codon, int pos);
std::pair<std::string,std::string>  Attempt_Mutate(char x);

long double random01(mt19937 generator)
{
    static uniform_01<mt19937> dist(generator);
    return dist();
}

std::vector<double> Get_OBS(){
	std::string s;
	std::vector<double> obs;
	std::ifstream myReadFile;
	myReadFile.open("Fitness.txt");
 	while (!myReadFile.eof()) {
		getline(myReadFile,s);
		std::string prob;
		for(int i=0; i<s.length();i++){	
			if(isdigit(s[i]) || s[i]=='.' ||  s[i]=='-'){
				prob=prob+s[i];
			}
		}
		//keep track of sum of the elements to matrix can be normalized
		std::transform(s.begin(), s.end(),s.begin(), ::toupper);
		s.erase(std::remove(s.begin(), s.end(), ' '),s.end());
		obs.push_back(atof(prob.c_str()));
	}
	myReadFile.close();		
	return obs;
}

extern std::vector<double> obs1(Get_OBS());

std::string Get_Eve_Seq(){
	std::string s;
	std::vector<double> obs;
	std::ifstream myReadFile;
	myReadFile.open("testseq");
	std::string h; //string contain human seq
	  while (!myReadFile.eof()) {
		std::getline(myReadFile,s);
		int H_flag=0;//set to 1 once done reading human seq;
		if(s[0]=='>'){
			H_flag++;
			std::getline(myReadFile,s);
			}
			if(H_flag==1){
			  h=h+s;
			}
			  if(H_flag==2){//if done reading human seq break
					break;
				}
			}
	myReadFile.close();		
	return h;
}

void Output(std::string evo, std::string EVE, int file_num){
       // std::cout << "attemptin to output" << std::endl;
	std::string name =boost::lexical_cast<std::string>(file_num)+".fasta";
	std::string file_name_out="/home/ateufel/Indy_Sim/all_align_fast/"+name;
	std::ofstream myOutFile;
	myOutFile.open(file_name_out.c_str());
	myOutFile << ">Org_H "<< std::endl;
	myOutFile << EVE << std::endl;
	myOutFile << ">New_H_score"<< std::endl;
	myOutFile << evo << std::endl;
	myOutFile.close();
}



std::string Get_Codon(char x){
	std::string AA=boost::lexical_cast<std::string>(x);
	int pos = rand() % 3;
	//std::cout << "attempt mutate: " << pos << std::endl;
	if(AA.compare("F") == 0){
		int index = rand() % 2;
		std::string f[]={"UUU","UUC"};
		return (f[index]);
	}
	if(AA.compare("L") == 0){
		int index = rand() % 6;
		std::string l[]={"UUA","UUG","CUU","CUC","CUA","CUG"};
		return( l[index]);
	}
	if(AA.compare("I") == 0){
		int index = rand() % 3;
		std::string i[]={"AUU","AUC","AUA"};
		return( i[index]);
	}
	if(AA.compare("M") == 0){
		return( "AUG");
	}
	if(AA.compare("V") == 0){
		int index = rand() % 4;
		std::string v[]={"GUU","GUC","GUA", "GUG"};
		return (v[index]);
	}
	if(AA.compare("S") == 0){
		int index = rand() % 6;
		std::string s[]={"UCU","UCC","UCA", "UCG","AGU","AGC"};
		return (s[index]);
	}
	if(AA.compare("P") == 0){
		int index = rand() % 4;
		std::string p[]={"CCU","CCC","CCA", "CCG"};
		return (p[index]);
	}
	if(AA.compare("T") == 0){
		int index = rand() % 4;
		std::string t[]={"ACU","ACC","ACA", "ACG"};
		return (t[index]);
	}
	if(AA.compare("A") == 0){
		int index = rand() % 4;
		std::string a[]={"GCU","GCC","GCA", "GCG"};
		return (a[index]);
	}
	if(AA.compare("Y") == 0){
		int index = rand() % 2;
		std::string y[]={"UAU","UAC"};
		return (y[index]);
	}
	if(AA.compare("H") == 0){
		int index = rand() % 2;
		std::string h[]={"CAU","CAC"};
		return (h[index]);
	}
	if(AA.compare("Q") == 0){
		int index = rand() % 2;
		std::string q[]={"CAA","CAG"};
		return (q[index]);
	}
	if(AA.compare("N") == 0){
		int index = rand() % 2;
		std::string n[]={"AAU","AAC"};
		return (n[index]);
	}
	if(AA.compare("K") == 0){
		int index = rand() % 2;
		std::string k[]={"AAA","AAG"};
		return (k[index]);
	}
	if(AA.compare("D") == 0){
		int index = rand() % 2;
		std::string d[]={"GAU","GAC"};
		return (d[index]);
	}
	if(AA.compare("E") == 0){
		int index = rand() % 2;
		std::string e[]={"GAA","GAG"};
		return (e[index]);
	}
	if(AA.compare("C") == 0){
		int index = rand() % 2;
		std::string c[]={"UGU","UGC"};
		return (c[index]);
	}
	if(AA.compare("W") == 0){
		return ("UGG");
	}
	if(AA.compare("R") == 0){
		int index = rand() % 6;
		std::string r[]={"CGU","CGC","CGA","CGG","AGA","AGG"};
		return (r[index]);
	}
	if(AA.compare("G") == 0){
		int index = rand() % 4;
		std::string g[]={"GGU","GGC","GGA","GGG"};
		return (g[index]);
	}
	if(AA.compare("-") == 0){
		return ("---");
	}
	if(AA.compare("X") == 0){ //if AA is unknown just introduce a gap
		return ("---");
	}
	std::cout << "broken: "<< AA << "assci: " << (int) x << " char: " << x <<  std::endl;
	return ("---");
}


std::pair<std::string, std::string> Attempt_Mutate(char x){
	std::string AA=boost::lexical_cast<std::string>(x);
	int pos = rand() % 3;
	//std::cout << "attempt mutate: " << pos << std::endl;
	if(AA.compare("F") == 0){
		int index = rand() % 2;
		std::string f[]={"UUU","UUC"};
		return codontoAA(f[index],pos);
	}
	if(AA.compare("L") == 0){
		int index = rand() % 6;
		std::string l[]={"UUA","UUG","CUU","CUC","CUA","CUG"};
		return codontoAA(l[index],pos);
	}
	if(AA.compare("I") == 0){
		int index = rand() % 3;
		std::string i[]={"AUU","AUC","AUA"};
		return codontoAA(i[index],pos);
	}
	if(AA.compare("M") == 0){
		return codontoAA("AUG",pos);
	}
	if(AA.compare("V") == 0){
		int index = rand() % 4;
		std::string v[]={"GUU","GUC","GUA", "GUG"};
		return codontoAA(v[index],pos);
	}
	if(AA.compare("S") == 0){
		int index = rand() % 6;
		std::string s[]={"UCU","UCC","UCA", "UCG","AGU","AGC"};
		return codontoAA(s[index],pos);
	}
	if(AA.compare("P") == 0){
		int index = rand() % 4;
		std::string p[]={"CCU","CCC","CCA", "CCG"};
		return codontoAA(p[index],pos);
	}
	if(AA.compare("T") == 0){
		int index = rand() % 4;
		std::string t[]={"ACU","ACC","ACA", "ACG"};
		return codontoAA(t[index],pos);
	}
	if(AA.compare("A") == 0){
		int index = rand() % 4;
		std::string a[]={"GCU","GCC","GCA", "GCG"};
		return codontoAA(a[index],pos);
	}
	if(AA.compare("Y") == 0){
		int index = rand() % 2;
		std::string y[]={"UAU","UAC"};
		return codontoAA(y[index],pos);
	}
	if(AA.compare("H") == 0){
		int index = rand() % 2;
		std::string h[]={"CAU","CAC"};
		return codontoAA(h[index],pos);
	}
	if(AA.compare("Q") == 0){
		int index = rand() % 2;
		std::string q[]={"CAA","CAG"};
		return codontoAA(q[index],pos);
	}
	if(AA.compare("N") == 0){
		int index = rand() % 2;
		std::string n[]={"AAU","AAC"};
		return codontoAA(n[index],pos);
	}
	if(AA.compare("K") == 0){
		int index = rand() % 2;
		std::string k[]={"AAA","AAG"};
		return codontoAA(k[index],pos);
	}
	if(AA.compare("D") == 0){
		int index = rand() % 2;
		std::string d[]={"GAU","GAC"};
		return codontoAA(d[index],pos);
	}
	if(AA.compare("E") == 0){
		int index = rand() % 2;
		std::string e[]={"GAA","GAG"};
		return codontoAA(e[index],pos);
	}
	if(AA.compare("C") == 0){
		int index = rand() % 2;
		std::string c[]={"UGU","UGC"};
		return codontoAA(c[index],pos);
	}
	if(AA.compare("W") == 0){
		return codontoAA("UGG",pos);
	}
	if(AA.compare("R") == 0){
		int index = rand() % 6;
		std::string r[]={"CGU","CGC","CGA","CGG","AGA","AGG"};
		return codontoAA(r[index],pos);
	}
	if(AA.compare("G") == 0){
		int index = rand() % 4;
		std::string g[]={"GGU","GGC","GGA","GGG"};
		return codontoAA(g[index],pos);
	}
	if(AA.compare("-") == 0){
		return codontoAA("---",0);
	}
	if(AA.compare("X") == 0){ //if AA is unknown just introduce a gap
		return codontoAA("---",0);
	}
	std::cout << "broken: "<< AA << "assci: " << (int) x << " char: " << x <<  std::endl;
	return codontoAA("---",0);
}

std::pair<std::string,std::string> codontoAA(std::string codon, int pos){
	//std::cout << "codon: " << codon << std::endl;
	std::string temp=codon; //transitions happen twice as often as trasnversions
	codon=codon[pos];
	int t = rand() % 6;
	if(t > 1 ){ //transition
		if(codon.compare("A") == 0) //a to g , c to u
			temp[pos]='G';
		if(codon.compare("G") == 0)
			temp[pos]='A';
		if(codon.compare("C") == 0)
			temp[pos]='U';
		if(codon.compare("U") == 0)
			temp[pos]='C';
	}
	if(t == 0) { //transversion
		if(codon.compare("A") == 0) //a to g , c to u
			temp[pos]='U';
		if(codon.compare("G") == 0)
			temp[pos]='C';
		if(codon.compare("C") == 0)
			temp[pos]='G';
		if(codon.compare("U") == 0)
			temp[pos]='A';
	}
	if(t == 1){
		if(codon.compare("A") == 0) //a to g , c to u
			temp[pos]='C';
		if(codon.compare("G") == 0)
			temp[pos]='U';
		if(codon.compare("C") == 0)
			temp[pos]='A';
		if(codon.compare("U") == 0)
			temp[pos]='G';
	}
	//std::cout << "temp: " << temp << std::endl;	
codon=temp;

if(codon.compare("UUU") == 0 ){
	std::pair<std::string,std::string> x("F","UUU");
	return x;
	 }
if(codon.compare( "UUC") == 0 ){
	std::pair<std::string,std::string> x("F","UUC");
	return x;
	 }
if(codon.compare( "UUA") == 0 ){
	std::pair<std::string,std::string> x("L","UUA");
	return x;
	 }
if(codon.compare( "UUG") == 0 ){
	std::pair<std::string,std::string> x("L","UUG");
	return x;
	 }
if(codon.compare( "CUU") == 0 ){
	std::pair<std::string,std::string> x("L","CUU");
	return x;
	 }
if(codon.compare( "CUC") == 0 ){
	std::pair<std::string,std::string> x("L","CUC");
	return x;
	 }
if(codon.compare( "CUA") == 0 ){
	std::pair<std::string,std::string> x("L","CUA");
	return x;
	 }
if(codon.compare( "CUG") == 0 ){
	std::pair<std::string,std::string> x("L","CUG");
	return x;
	 }
if(codon.compare( "AUU") == 0 ){
	std::pair<std::string,std::string> x("I","AUU");
	return x;
	 }
if(codon.compare( "AUC") == 0 ){
	std::pair<std::string,std::string> x("I","AUC");
	return x;
	 }
if(codon.compare( "AUA") == 0 ){
	std::pair<std::string,std::string> x("I","AUA");
	return x;
	 }
if(codon.compare( "AUG") == 0 ){
	std::pair<std::string,std::string> x("M","AUG");
	return x;
	 }
if(codon.compare( "GUU") == 0 ){
	std::pair<std::string,std::string> x("V","GUU");
	return x;
	 }
if(codon.compare( "GUC") == 0 ){
	std::pair<std::string,std::string> x("C","GUC");
	return x;
	 }
if(codon.compare( "GUA") == 0 ){
	std::pair<std::string,std::string> x("V","GUA");
	return x;
	 }
if(codon.compare( "GUG") == 0 ){
	std::pair<std::string,std::string> x("V","GUG");
	return x;
	 }
if(codon.compare( "UCU") == 0 ){
	std::pair<std::string,std::string> x("S","UCU");
	return x;
	 }
if(codon.compare( "UCC") == 0 ){
	std::pair<std::string,std::string> x("S","UCC");
	return x;
	 }
if(codon.compare( "UCA") == 0 ){
	std::pair<std::string,std::string> x("S","UCA");
	return x;
	 }
if(codon.compare( "UCG") == 0 ){
	std::pair<std::string,std::string> x("S","UCG");
	return x;
	 }
if(codon.compare( "CCU") == 0 ){
	std::pair<std::string,std::string> x("P","CCU");
	return x;
	 }
if(codon.compare( "CCC") == 0 ){
	std::pair<std::string,std::string> x("P","CCC");
	return x;
	 }
if(codon.compare( "CCA") == 0 ){
	std::pair<std::string,std::string> x("P","CCA");
	return x;
	 }
if(codon.compare( "CCG") == 0 ){
	std::pair<std::string,std::string> x("P","CCG");
	return x;
	 }
if(codon.compare( "ACU") == 0 ){
	std::pair<std::string,std::string> x("T","ACU");
	return x;
	 }
if(codon.compare( "ACC") == 0 ){
	std::pair<std::string,std::string> x("T","ACC");
	return x;
	 }
if(codon.compare( "ACA") == 0 ){
	std::pair<std::string,std::string> x("T","ACA");
	return x;
	 }
if(codon.compare( "ACG") == 0 ){
	std::pair<std::string,std::string> x("T","ACG");
	return x;
	 }
if(codon.compare( "GCU") == 0 ){
	std::pair<std::string,std::string> x("A","GCU");
	return x;
	 }
if(codon.compare( "GCC") == 0 ){
	std::pair<std::string,std::string> x("A","GCC");
	return x;
	 }
if(codon.compare( "GCA") == 0 ){
	std::pair<std::string,std::string> x("A","GCA");
	return x;
	 }
if(codon.compare( "GCG") == 0 ){
	std::pair<std::string,std::string> x("A","GCG");
	return x;
	 }
if(codon.compare( "UAU") == 0 ){
	std::pair<std::string,std::string> x("Y","UAU");
	return x;
	 }
if(codon.compare( "UAC") == 0 ){
	std::pair<std::string,std::string> x("Y","UAC");
	return x;
	 }
if(codon.compare( "CAU") == 0 ){
	std::pair<std::string,std::string> x("H","CAU");
	return x;
	 }
if(codon.compare( "CAC") == 0 ){
	std::pair<std::string,std::string> x("H","CAC");
	return x;
	 }
if(codon.compare( "CAA") == 0 ){
	std::pair<std::string,std::string> x("Q","CAA");
	return x;
	 }
if(codon.compare( "CAG") == 0 ){
	std::pair<std::string,std::string> x("Q","CAG");
	return x;
	 }
if(codon.compare( "AAU") == 0 ){
	std::pair<std::string,std::string> x("N","AAU");
	return x;
	 }
if(codon.compare( "AAC") == 0 ){
	std::pair<std::string,std::string> x("N","AAC");
	return x;
	 }
if(codon.compare( "AAA") == 0 ){
	std::pair<std::string,std::string> x("K","AAA");
	return x;
	 }
if(codon.compare( "AAG") == 0 ){
	std::pair<std::string,std::string> x("K","AAG");
	return x;
	 }
if(codon.compare( "GAU") == 0 ){
	std::pair<std::string,std::string> x("D","GAU");
	return x;
	 }
if(codon.compare( "GAC") == 0 ){
	std::pair<std::string,std::string> x("D","GAC");
	return x;
	 }
if(codon.compare( "GAA") == 0 ){
	std::pair<std::string,std::string> x("E","GAA");
	return x;
	 }
if(codon.compare( "GAG") == 0 ){
	std::pair<std::string,std::string> x("E","GAG");
	return x;
	 }
if(codon.compare( "UGU") == 0 ){
	std::pair<std::string,std::string> x("C","UGU");
	return x;
	 }
if(codon.compare( "UGC") == 0 ){
	std::pair<std::string,std::string> x("C","UGC");
	return x;
	 }
if(codon.compare( "UGG") == 0 ){
	std::pair<std::string,std::string> x("W","UGG");
	return x;
	 }
if(codon.compare( "CGU") == 0 ){
	std::pair<std::string,std::string> x("R","CGU");
	return x;
	 }
if(codon.compare( "CGC") == 0 ){
	std::pair<std::string,std::string> x("R","CGC");
	return x;
	 }
if(codon.compare( "CGA") == 0 ){
	std::pair<std::string,std::string> x("R","CGA");
	return x;
	 }
if(codon.compare( "CGG") == 0 ){
	std::pair<std::string,std::string> x("R","CGG");
	return x;
	 }
if(codon.compare( "AGU") == 0 ){
	std::pair<std::string,std::string> x("S","AGU");
	return x;
	 }
if(codon.compare( "AGC") == 0 ){
	std::pair<std::string,std::string> x("S","AGC");
	return x;
	 }
if(codon.compare( "AGA") == 0 ){
	std::pair<std::string,std::string> x("R","AGA");
	return x;
	 }
if(codon.compare( "AGG") == 0 ){
	std::pair<std::string,std::string> x("R","AGG");
	return x;
	 }
if(codon.compare( "GGU") == 0 ){
	std::pair<std::string,std::string> x("G","GGU");
	return x;
	 }
if(codon.compare( "GGC") == 0 ){
	std::pair<std::string,std::string> x("G","GGC");
	return x;
	 }
if(codon.compare( "GGA") == 0 ){
	std::pair<std::string,std::string> x("G","GGA");
	return x;
	 }
if(codon.compare( "GGG") == 0 ){
	std::pair<std::string,std::string> x("G","GGG");
	return x;
	 }
//all of the stop codons!
if(codon.compare( "UAA") == 0 ){
	std::pair<std::string,std::string> x("x","UAA");
	return x;
	 }
if(codon.compare( "UAG") == 0 ){
	std::pair<std::string,std::string> x("x","UAG");
	return x;
	 }
if(codon.compare( "UGA") == 0 ){
	std::pair<std::string,std::string> x("x","UGA");
	return x;
	 }

 std::cout << "mutation didn not make sense returning - : " << codon << std::endl;
 std::pair<std::string,std::string> x("-","---");
 return x;
} //end function 



double Get_Score(char aa){
	aa=tolower(aa);//convert to lower case
	char AA[] = {'a', 's', 't', 'c', 'v', 'l', 'i', 'm', 'p', 'f', 'y', 'w', 'd', 'e', 'n', 'q', 'h', 'k', 'r', 'g'};
	std::vector<char> AAs(AA,AA+20);
	std::vector<char>::iterator it;
	it=find(AAs.begin(),AAs.end(),aa); //find what the amino acids position is
	size_t i=it-AAs.begin();
	std::vector<double> Fixation(obs1); 
	return(Fixation[i]); //return the score
}
