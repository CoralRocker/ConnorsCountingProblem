#include <iostream>
#include <cmath>
#include <string>
#include <array>
#include <vector>
#include <iomanip>
#include <future>
#include <chrono>

using namespace std;

typedef vector<int> adjstr;

chrono::milliseconds halfsec(500);

/**
 * @brief Returns the factors of n as a vector, not including 1
 * 
 * @param n number to factorize
 * @return vector<int> factors, not including 1
 */
vector<int> factors(const int n){
    adjstr fact;
    for( int i = (n/2); i > 1; i-- ){
        if( (n/i)*i == n )
            fact.push_back(i);
    }
    return fact;
}

/**
 * @brief Print the n x n matrix
 * 
 * @param matrix square input matrix
 * @param n size of side of matrix
 */
void printmat(const int* matrix, const int n){
    for( int i = 0; i < n; i++ ){
        cout << "[ ";
        for( int j = 0; j < n; j++ ){
            cout << setw(3) << matrix[i*n + j] <<" ";
        }
        cout << "]\n";
    }
    cout << endl;
}

/**
 * @brief Returns an array of the adjacent positions in the matrix.
 * 
 * @param matrix Square input matrix, flattened to vector
 * @param n size of matrix side
 * @param pos which point in the matrix to get the adjacency for. 
 * @return array<int*, 4> pointers to all adjacent items in the matrix. Nullptr if a position is out of bounds
 */
array<const int*, 4> adjacent(const int* matrix, int n, int pos){
    array<const int*, 4> adjacents;
    adjacents.fill(nullptr);

    if(  pos > 0 && ((pos - 1) % n) < (pos % n)){
        adjacents[3] = &matrix[pos-1];
    }
    if( ((pos+1)%n) > (pos%n) ){
        adjacents[1] = &matrix[pos+1];
    }
    if( (pos+n) < (n*n) ){
        adjacents[2] = &matrix[pos+n];
    }
    if( (pos-n) >= 0 ){
        adjacents[0] = &matrix[pos-n];
    }

    return adjacents;
}


/**
 * @brief Adds an item to a vector, returning a copy
 * 
 * @tparam T Type held in vector
 * @param vec vector to add to
 * @param val value to push_back
 * @return vector<T> New vector produced
 */
template<class T>
vector<T> operator+(const vector<T> vec, const T val){
    vector<T> v(vec);
    v.push_back(val);
    return v;
}

/**
 * @brief Get a position vector from an adjacency array
 * 
 * @param matrix square matrix, flattened to vector
 * @param adj size of sides of matrix
 * @return adjstr a vector<int> type with all valid indexes in it
 */
adjstr posFromPtr(const int* matrix, const array<const int*, 4> adj){
    adjstr pos;
    for(const int* ptr : adj ){
        if( ptr ){
            pos.push_back(ptr - matrix);
        }
    }
    return pos;
}

/**
 * @brief Checks that two adjacency strings are equal.
 * Note: Adjacency is the same forward and backward
 * 
 * @param vec1 
 * @param vec2 
 * @return true vec1 and vec2 are equal
 * @return false vec1 and vec2 aren't equal
 */
bool equal(const adjstr vec1, const adjstr vec2){
    if( vec1.size() != vec2.size() ) return false;
    bool same = true, rsame = true;
    for( unsigned i = 0; i < vec1.size(); i++ ){
        if( vec1.at(i) != vec2.at(i) ) same = false;
        if( vec1.at(i) != *(vec2.rbegin()+i) ) rsame = false;
    }

    return same || rsame;
}

/**
 * @brief Checks if all elements of vec1 are different that all elements of vec2
 * 
 * @param vec1 
 * @param vec2 
 * @return true all elements are distinct
 * @return false at least 2 elements are the same
 */
bool disjoint(const adjstr vec1, const adjstr vec2){
    for( int i : vec1 ){
        for( int j : vec2 ){
            if( i == j )
                return false;
        }
    }
    return true;
}

/**
 * @brief Check if a value is in an adjacency string
 * 
 * @param vec adjacency string
 * @param test value to test
 * @return true value was found 
 * @return false value wasn't found
 */
bool in(const adjstr vec, int test){
    for( int t : vec ){
        if( t == test )
            return true;
    }
    return false;
}

/**
 * @brief Recursively create an adjacency string of a given depth. Adds the string to mainvec when it is at the right size
 * 
 * @param matrix matrix to get info from
 * @param n size of matrix side
 * @param mainvec vector of all adjacency strings
 * @param adj current adjacency string to add to
 * @param depth current depth in adjacency string
 */
void recurseAddString(const int* matrix, int n, vector<adjstr> &mainvec, adjstr adj, int depth){
    if( depth == 0 ){
        mainvec.push_back(adj);
        return;
    }
    auto moves = posFromPtr(matrix, adjacent(matrix, n, adj.back()));
    for( int move : moves ){
        if( !in(adj, move) )
            recurseAddString(matrix, n, mainvec, adj+move, depth-1);
    }
}

/**
 * @brief Get the all the adjacency strings of a given size in the matrix for a given position
 * 
 * @param matrix square matrix
 * @param n size of side in matrix
 * @param pos position to start from
 * @param size size of adjacency string
 * @return vector<adjstr> all the adjacency strings
 */
vector<adjstr> getAdjStrings(const int* matrix, int n, int pos, int size){
    vector<adjstr> strings;
    recurseAddString(matrix, n, strings, {pos}, size-1);
    return strings;
}

/**
 * @brief Prints out an adjacency string
 * 
 * @param matrix matrix 
 * @param adj adjacency string
 */
void printAdj(adjstr adj){
    cout << "[ ";
    for( int pos : adj ){
        cout << pos << " ";
    }
    cout << "]" << endl;
}

/**
 * @brief Merges two adjacency string vectors such that all adjstrs in the result are unique
 * 
 * @param rhs 
 * @param lhs 
 * @return vector<adjstr> 
 */
vector<adjstr> merge(const vector<adjstr> &rhs, const vector<adjstr> &lhs){
    vector<adjstr> result(lhs);
    result.reserve(rhs.size() + lhs.size());

    for( auto &adjstring : rhs ){
        bool found = false;
        for( auto &cmpstring : lhs ){
            found = equal(adjstring, cmpstring);
        }
        if( !found ){
            result.push_back(adjstring);
        }
    }
    return result;
}

/**
 * @brief Get all possible adjacency strings of given size for the matrix
 * 
 * @param matrix square matrix, flattened to a vector
 * @param n size of one side of matrix
 * @param size size of adjacency strings
 * @return vector<adjstr> 
 */
vector<adjstr> getAllAdjStrings(const int* matrix, int n, int size){
    vector<adjstr> result;
    for( int i = 0; i < n*n; i++ ){
        result = merge(result, getAdjStrings(matrix, n, i, size));
    }
    return result;
}

/**
 * @brief Get setsleft disjoint adjstrs and add them to a potential solution. If setsleft is 0, a solution is added to the solution set.
 * Note: recursive
 * 
 * @param solnset full set of all the solutions to the problem
 * @param allstrings All adjacency strings
 * @param sol potential solution being tested
 * @param setsleft number of sets left to be added to solution
 */
void getDisjointSets(vector<vector<adjstr>> &solnset, const vector<adjstr> &allstrings, vector<adjstr> sol, int setsleft){
    if( setsleft == 0 ){
        solnset.push_back(sol);
        return;
    }
    for( const adjstr &adj : allstrings ){
        bool dis4all = true;
        for( adjstr &tmp : sol ){
            if( !disjoint(adj, tmp) ){ 
                dis4all = false;
                break;
            }
        }
        if( dis4all ){
            getDisjointSets(solnset, allstrings, sol + adj, setsleft - 1);
        }
    }
}

/**
 * @brief Get all possible solutions to the adjacency string problem of a given size
 * 
 * @param allstrings all adjacency strings for the matrix
 * @param solsize number of adjacency strings in a solution
 * @return vector<vector<adjstr>> 
 */
vector<vector<adjstr>> getSolutions(const vector<adjstr> allstrings, const int solsize){
    vector<vector<adjstr>> result;
    for( const adjstr &adj : allstrings ){
        vector<adjstr> sol = {adj};
        getDisjointSets(result, allstrings, sol, solsize-1);
    }
    return result;
}

/**
 * @brief Get all possible solutions to the adjacency problem for the given starting strings
 * 
 * @param allstrings all adjacency strings for the matrix
 * @param taskset adjacency strings to test for
 * @param solsize number of strings in solution
 * @return vector<vector<adjstr>> 
 */
vector<vector<adjstr>> getSolutionsInTask(const vector<adjstr> allstrings, const vector<adjstr> taskset, const int solsize){
    vector<vector<adjstr>> result;
    for( const adjstr &adj : taskset ){
        vector<adjstr> sol = {adj};
        getDisjointSets(result, allstrings, sol, solsize-1);
    }
    return result;
}

/**
 * @brief Prints a solution vector
 * 
 * @param sol solution vector of adjacency strings
 */
void printsol(const vector<adjstr> sol){
    cout << "{ " << endl;
    for( adjstr adj : sol ){
        cout << "\t[ ";
        for( int i : adj ){
            cout << i << " ";
        }
        cout << "]" << endl;
    }
    cout << "}" << endl;
}

/**
 * @brief Returns a unique integer hash of an adjacency string.
 * 
 * @param adj String to hash. Must be with the smallest end point at the start.
 * @return unsigned long long hashed string
 */
unsigned long long hashAdjstr(const adjstr &adj){

    int g = 31;
    unsigned long long hash = 0;
    for( int pos : adj ){
        hash = g*hash + pos;
    }

    return hash;
}

/**
 * @brief Returns a hash of all the adjacency strings in a solution set, normalized and summed.
 * 
 * @param sol 
 * @return unsigned long long 
 */
unsigned long long hashSol(const vector<adjstr> &sol){
    unsigned long long hash = 0;
    for( const adjstr &adj : sol ){
        adjstr tmp;
        if( adj.back() < adj.front() ){
            tmp.insert(tmp.begin(), adj.crbegin(), adj.crend());
        }else{
            tmp.insert(tmp.begin(), adj.cbegin(), adj.cend());
        }
        hash += hashAdjstr(tmp);
    }
    return hash;
}

/**
 * @brief Reduces the set of all possible solutions to the set of all unique solutions
 * 
 * @param solns Set of all possible solutions
 * @return vector<vector<adjstr>> 
 */
vector<vector<adjstr>> uniqueSolns(const vector<vector<adjstr>> &solns){
    vector<vector<adjstr>> result;
    vector<unsigned long long> result_checksums;
    vector<unsigned long long> checksums(solns.size());
    for( unsigned i = 0; i < solns.size(); i++ ){
        checksums.at(i) = hashSol(solns.at(i));

    }

    for( auto it = solns.begin(); it < solns.end(); it++ ){
        unsigned long long hash = checksums.at(it - solns.begin());
        bool in = false;
        for( unsigned i = 0; i < result.size(); i++ ){
            if( hash == result_checksums.at(i) ){
                in = true;
                break;
            }
        }
        if( !in ){
            result.push_back(*it);
            result_checksums.push_back(hash);
        }
    }
    return result;
}

/**
 * @brief Splits the given set of adjacency strings into numThreads new vectors
 * 
 * @param allstrings the set of all adjacency strings
 * @param numThreads number of vectors to split into
 * @return vector<vector<adjstr>> 
 */
vector<vector<adjstr>> splitTask(const vector<adjstr> &allstrings, int numThreads){
    vector<vector<adjstr>> result(numThreads);
    int counter = 0;
    for(const adjstr &adj : allstrings ){
        result[counter%numThreads].push_back(adj);
        counter++;
    }
    return result;
}

/**
 * @brief Get the solutions using multithreading
 * 
 * @param allstrings 
 * @param numThreads 
 * @param solsize 
 * @return vector<vector<adjstr>> 
 */
vector<vector<adjstr>> getSolutionsMT(const vector<adjstr> allstrings, const int numThreads, const int solsize){
    auto tasksets = splitTask(allstrings, numThreads);
    vector<future<vector<vector<adjstr>>>> results(numThreads);
    vector<vector<adjstr>> resultant;

    for( int t = 0; t < numThreads; t++ ){
        results.at(t) = async(getSolutionsInTask, allstrings, tasksets.at(t), solsize);
    }

    for( auto &res : results ){
        while( res.wait_for(halfsec) == future_status::timeout );
        auto soln = res.get();
        resultant.insert(resultant.end(), soln.begin(), soln.end());
    }

    return resultant;
}



int main(int argc, char** argv){
    int n = 4;
    int num_threads = 16;
    if( argc > 1 )
        n = stoi(argv[1]);
    
    int* matrix = new int[n*n];
    for( int i = 1; i <= n*n; i++ ){ matrix[i-1] = i; }

    printmat(matrix, n);

    cout << "Factors of " << n*n << ": ";
    auto fact = factors(n*n);
    for( int i : fact ){
        cout << i << " ";
    }
    cout << endl;

    for( int factor : fact ){
        cout << "!!! FOR STRINGS OF SIZE " << factor << " !!!" << endl;
        auto all = getAllAdjStrings(matrix, n, factor); 
        cout << "There are " << all.size() << " adjacency strings" << endl;        
        auto sols = getSolutionsMT(all, num_threads, (n*n)/factor);
        cout << "Solutions: " << sols.size() << endl;
        auto unq = uniqueSolns(sols);
        cout << "Unique solutions: " << unq.size() << endl;
    }

    delete[] matrix;
}