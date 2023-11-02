//
//  CNFparser.cpp
//  SpMM MAX-SAT with zero padding
//
#include "CNFparser.hpp"

ParsedVector CNFparser(std::string filepath, bool index){
    ParsedVector parsed_result;
    
    std::ifstream ProblemFile;
    ProblemFile.open(filepath);
    if (ProblemFile.fail())
    {
       std::cout << "cnf file not found!"
                 << "\n";
       return parsed_result;
    }
    
    int cnt_lit;
    bool weighted = false;
    parsed_result.k = 0;
    parsed_result.cnt_vector.push_back(0);
    parsed_result.T_cnt_vector.push_back(0);
    
    
    
    
    // parser
    for ( std::string line; std::getline( ProblemFile, line ); ) {
          cnt_lit = 0;
          if ( line[0] == 'c') {
             continue;
          }
        
        if ( line[0] == 'p' ) {
            std::string config = line.substr(line.find("cnf")+4);
            std::stringstream iss( config );
            iss >> parsed_result.num_var;
            iss >> parsed_result.num_clauses;
            // check if weighted
            if (line.find("p wcnf") != std::string::npos){
                weighted = true;
            }
            continue;
        }
        
        // store the weight in the weight vector
        if(weighted){
            std::stringstream iss( line );
            int weight;
            iss >> weight;
            parsed_result.weight_vector.push_back(weight);
            line = line.substr(1);
        }else{
            parsed_result.weight_vector.push_back(1); // for unweighted instance
        }
          

          std::stringstream iss( line );
          int number;
          while ( iss >> number ) {
              if ( number == 0 ) {
                  if (cnt_lit > parsed_result.k){
                      parsed_result.k = cnt_lit;
                  }
                   parsed_result.cnt_vector.push_back(cnt_lit+parsed_result.cnt_vector.back());
                   break; // end of line
              }
              
              // add index option to choose between two encodings
              if(index){
                      parsed_result.problem_vector.push_back(number/abs(number));
                      parsed_result.index_vector.push_back(abs(number));
              }
              else{
                  parsed_result.problem_vector.push_back(number);
              }
              cnt_lit++;
          }
    }
    ProblemFile.close();
    
    
    // Generate CSC transposed matrix with CSR matrix
    //transposed matrix 2D vector
    std::vector<std::vector<int>> T_Cmi(parsed_result.num_var);
    for(int j = 0; j < parsed_result.num_clauses; j ++){
        for(int i = parsed_result.cnt_vector[j]; i < parsed_result.cnt_vector[j+1]; i ++ ){
            int rind = parsed_result.problem_vector[i];
            T_Cmi[abs(rind)-1].push_back((j+1) * rind/abs(rind));
        }
    }
    
    for(int i = 0; i < parsed_result.num_var; i ++){
        // populate CSC transposed matrix
        for (std::vector<int>::iterator it = T_Cmi[i].begin() ; it != T_Cmi[i].end(); ++it){
            parsed_result.T_problem_vector.push_back(*it);
        }
        
        // populate row cnt vector for transposed matrix
        parsed_result.T_cnt_vector.push_back(int(T_Cmi[i].size())+parsed_result.T_cnt_vector.back());
    }
    
     
    
    return parsed_result;
}






int* zero_padding(std::vector<int> dense_array, std::vector<int> ind_ptr, int padded_coldim, int padded_rowdim){
    int* padded_array = (int *) malloc(padded_coldim * padded_rowdim * sizeof(int));

    for (int i = 0; i < padded_rowdim; i ++){
        for(int j = 0; j < padded_coldim; j ++){
            if(ind_ptr[i+1]-ind_ptr[i] < j+1){
                padded_array[i*padded_coldim+j] = 0;
            }
            else{
                padded_array[i*padded_coldim+j] = dense_array[ind_ptr[i]+j];
                
            }
        }
    }
    
    return padded_array;
    
}



int solution_checker(const float* s, ParsedVector parsed_result,int num_trial){

    int maxsat_COST= 0; // weighted sum of unsatisified clauses

    int cost = 0;
    int max_runid = 0;
    
    float s_check[parsed_result.num_var* num_trial];
    std::memcpy(s_check, s, sizeof(float) * num_trial* parsed_result.num_var);

    // Solution encoding
    for (int n = 0; n < num_trial; n ++){
        for(int i = 0; i < parsed_result.num_var; i ++) {
            if (s[i*num_trial+n] > 0)
                s_check[i*num_trial+n] = 1;
            if (s[i*num_trial+n] < 0)
                s_check[i*num_trial+n] = -1;
            if (s[i*num_trial+n] == 0)
                return max_runid;
        }
    }
    
    // Sparse matrix multiplication checking satisfaction condition
    for(int n = 0; n < num_trial; n ++){
        cost = 0;
        
        // for single trial
        for (int i = 0; i <  parsed_result.num_clauses; i ++){
            float res = 0;
            for (int j = 0; j <  parsed_result.k; j ++){
                if ( parsed_result.problem_vector[i* parsed_result.k+j] != 0){
                    if( parsed_result.problem_vector[j* parsed_result.k+i] > 0){
                        res += 1 * s_check[( parsed_result.problem_vector[j* parsed_result.k + i]-1)*num_trial+n];
                    }
                    if(parsed_result.problem_vector[j*parsed_result.k+i] < 0){
                        res += (-1) * s_check[(-parsed_result.problem_vector[j* parsed_result.k + i]-1)*num_trial+n];
                    }
                        
                }
            }
           
            // check whether the clause is satisfied
            // if not add the cost
            if (res > -(parsed_result.cnt_vector[i+1] - parsed_result.cnt_vector[i])){
                continue;
            }
            else{
                cost += parsed_result.weight_vector[i];
            }
        }
        
        if (n==0)
        {
            maxsat_COST = cost;
            max_runid = n;
        }
        
        // find the num_maxsat, min_cost across trials
        
        if (cost < maxsat_COST){
            maxsat_COST = cost;
            max_runid = n;
        }
  
    }
    
    return maxsat_COST;
}
