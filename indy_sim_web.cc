#include <iostream>
#include <fstream>
#include <vector>
#include <valarray>
#include <cmath>
#include <boost/tuple/tuple.hpp>
#include <boost/concept_check.hpp>
#include <utility>
#include <boost/signals.hpp>
#include <numeric>
#include <algorithm>
#include "math.h"
#include <string>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <dirent.h> 
#include <cstdio>
#include <ctime>
#include <boost/math/distributions/binomial.hpp>
#include <iostream>
#include <boost/random.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include "indy_sim_web.h"
//most of these are unnessary. 

#define pop_size 100
//g++ indy_sim_web.cc  

int
main (int argc, char *argv[])
{ 
  srand((unsigned)time(0)); //set random generators up
  mt19937 generator(time(0));
  std::string Eve=Get_Eve_Seq();
  std::cout << "EVE: " << Eve << std::endl;
  std::vector<std::string> Final(Eve.length(),"-");
 // #pragma omp parallel for 
for(int s=0; s<Eve.length(); s++){
	char AA=Eve[s];
	std::cout << AA << std::endl;
	if(AA=='-'){ //if there is a gap then just add the gap
		//Final=Final+'-';
		Final[s]="-";
	}

	if(AA!='-'){ //make sure its not a gap   
		Indiv_struct I;
		I.AA=AA;
		I.codon = Get_Codon(AA);
		I.ind_score=1; //initall the score is 1. the change between 1+(to-from) happens when there is a mutation
		std::vector<Indiv_struct> Pop;	
		Pop.assign(pop_size, I);
		
	int flag=0; //mutation flag  
	
for(int time=0; time<200000; time++){ //assuming 6 million years of evolution and a generation time of 30 years 200,000 generations will occur
	clock_t start = clock();
	for(int i=0; i<pop_size; i++){
	//std::cout << "Indv: " << i << std::endl;
		// figure out if AA mutates
		const long double p = 1e-6; // mutation rate from Estimate of the Mutation Rate Per Nucleotide in Humans by Nachman and Crowell 1e-9
		long double k = random01(generator);
  		if(k < p){ //if value is smaller then a mutation has happened
			flag=1;
			int mut_pos = rand() % 3; //choose a random position in the codon to mutate
			std::pair<std::string,std::string> mut = codontoAA(Pop[i].codon,mut_pos); //return AA and the codon
			if(mut.first[0]=='x'){ //a mutation to a stop codon has occured, this is lethal so remove it from the population
				Pop[i].ind_score=0;
				break;
			}
			
			if(Pop[i].AA!=mut.first[0]){ // the mutation has changed the amino acid
				double fit=1+( Get_Score(mut.first[0]) - Get_Score(Pop[i].AA)); //calculate what the indivuals new score is going to be
				Pop[i].AA=mut.first[0];//get the new AA
				Pop[i].codon=mut.second;//get the codon
				if(Pop[i].AA == AA){ //if the mutation is back to wild type, note wild type is defined as whatever is currently there, 2nd lvl effects are not considered, this is a problem to be fixed later.
					Pop[i].ind_score=1;
				}
				else{ //if the mutation isnt to wild type
					Pop[i].ind_score=fit;//update the score
				}
				break;
			}
			if(Pop[i].AA==mut.first[0]){ //the mutation has not change the amino acid, but may have change the codon
				Pop[i].codon=mut.second;//get the codon
				break;
			}
			
  		}//end mutation
	}//end i

	if(flag!=0){
		double s_tot=0;
	
		for(int i=0; i < pop_size; i++){
			s_tot+=Pop[i].ind_score; //figure out total so weighted selection can happen
		}
		
		std::vector<double> prob1;
		for(int i=0; i < pop_size; i++){// normalize weight scores
			if(Pop[i].ind_score > 0){ 
			prob1.push_back(Pop[i].ind_score/s_tot); 
			}
			if(Pop[i].ind_score == 0){
			prob1.push_back(0); //if a mutation has happened to a stop codon, this is leathal and these indivuals can not be picked to be part of the next generats
			} 
		}
		// set up weighted distribution to select indexes
		//http://www.boost.org/doc/libs/1_43_0/doc/html/boost_random/tutorial.html
		std::vector<double> prob;
		std::partial_sum(&prob1[0], &prob1[0] + pop_size, std::back_inserter(prob));
		boost::uniform_real<> dist(0, prob.back());
		boost::variate_generator<boost::mt19937&, boost::uniform_real<> > indexs(generator, dist);
	
		std::vector<Indiv_struct> Pop_temp;
		//Pop_temp.reserve(pop_size); //i dont rememeber why i had this here..it might have something to do with fast mem manganment. 
		
		//fill up temp with weighted choosen indexs
		for(int i=0; i < pop_size; i++) {
		  Pop_temp.push_back(Pop[std::lower_bound(prob.begin(), prob.end(), indexs()) - prob.begin()]);
		}
		//dont copy swap. its way faster
		std::swap(Pop,Pop_temp);
	}//end mutation flag

	}//end time

	//Final=Final+Pop[1].AA; //increase the alignment by 1
	
	Final[s]=Pop[1].AA;
	std::cout << "Final:" << Pop[1].AA<< std::endl;
	std::cout << "score: " << Pop[1].ind_score << std::endl;
	}//end single AA
}//end sequence


	//Output(Final,Eve,atoi(argv[1]));
	for(int i=0; i < Final.size(); i++){
	  std::cout << Final[i];
	}
	std::cout << std::endl;

}//end sim
