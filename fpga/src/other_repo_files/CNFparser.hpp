//
//  CNFparser.hpp
//  SpMM MAX-SAT with zero padding
//

#ifndef CNFparser_hpp
#define CNFparser_hpp
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <stdlib.h>
#include <fstream>
#include <unistd.h>
#include <iostream>
#include <assert.h>
#include <stdio.h>
using namespace std;
typedef struct ParsedVector {
    int num_var;
    int num_clauses;
    int k;
    std::vector<int> problem_vector;
    std::vector<int> index_vector;
    std::vector<int> cnt_vector;
    
    std::vector<int> T_problem_vector;
    std::vector<int> T_cnt_vector;
     
    std::vector<int> weight_vector;
}ParsedVector;

ParsedVector CNFparser(std::string filepath, bool index);
int* zero_padding(std::vector<int> dense_array, std::vector<int> ind_ptr, int padded_coldim, int padded_rowdim);
int solution_checker(const float* s, ParsedVector parsed_result,int num_trial);
#endif /* CNFparser_hpp */
