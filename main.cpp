#include <iostream>
#include <cmath>
#include <string>
#include <array>
#include <vector>
#include <iomanip>

using namespace std;

typedef vector<int> adjstr;


vector<int> factors(int n){
    adjstr fact;
    for( int i = 2; i <= (n/2); i++ ){
        if( (n/i)*i == n )
            fact.push_back(i);
    }
    return fact;
}

void printmat(int* matrix, int n){
    for( int i = 0; i < n; i++ ){
        cout << "[ ";
        for( int j = 0; j < n; j++ ){
            cout << setw(3) << matrix[i*n + j] <<" ";
        }
        cout << "]\n";
    }
    cout << endl;
}

array<int*, 4> adjacent(int* matrix, int n, int pos){
    array<int*, 4> adjacents;
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
template<class T>
vector<T> operator+(const vector<T> vec, const T val){
    vector<T> v(vec);
    v.push_back(val);
    return v;
}

adjstr posFromPtr(int* matrix, array<int*, 4> adj){
    adjstr pos;
    for( int* ptr : adj ){
        if( ptr ){
            pos.push_back(ptr - matrix);
        }
    }
    return pos;
}

bool equal(adjstr vec1, adjstr vec2){
    bool same = true, rsame = true;
    for( int i = 0; i < vec1.size(); i++ ){
        if( vec1.at(i) != vec2.at(i) ) same = false;
        if( vec1.at(i) != *(vec2.rbegin()+i) ) rsame = false;
    }

    return same || rsame;
}

bool disjoint(adjstr vec1, adjstr vec2){
    for( int i : vec1 ){
        for( int j : vec2 ){
            if( i == j )
                return false;
        }
    }
    return true;
}

bool in(adjstr vec, int test){
    for( int t : vec ){
        if( t == test )
            return true;
    }
    return false;
}

void recurseAddString(int* matrix, int n, vector<adjstr> &mainvec, adjstr adj, int depth){
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

vector<adjstr> getAdjStrings(int* matrix, int n, int pos, int size){
    vector<adjstr> strings;
    recurseAddString(matrix, n, strings, {pos}, size-1);
    return strings;
}


void printAdj(int* matrix, adjstr adj){
    cout << "[ ";
    for( int pos : adj ){
        cout << pos << " ";
    }
    cout << "]" << endl;
}

vector<adjstr> merge(const vector<adjstr> &rhs, const vector<adjstr> &lhs){
    vector<adjstr> result(lhs);
    result.reserve(rhs.size() + lhs.size());

    for( auto &adjstring : rhs ){
        bool found = false;
        for( auto &cmpstring : lhs ){
            if( adjstring.size() != cmpstring.size() ) continue;

            bool same = true, rsame = true; // forward and reverse sameness
            for( int i = 0; i < adjstring.size(); i++ ){
                if( adjstring.at(i) != cmpstring.at(i) ) same = false;
                if( adjstring.at(i) != *(cmpstring.rbegin()+i) ) rsame = false;
            }
            if( same || rsame ){
                found = true;
            }

        }
        if( !found ){
            result.push_back(adjstring);
        }
    }
    return result;
}

vector<adjstr> getAllAdjStrings(int* matrix, int n, int size){
    vector<adjstr> result;
    for( int i = 0; i < n*n; i++ ){
        result = merge(result, getAdjStrings(matrix, n, i, size));
    }
    return result;
}

void getDisjointSets(vector<vector<adjstr>> &solnset, vector<adjstr> &allstrings, vector<adjstr> sol, int setsleft){
    if( setsleft == 0 ){
        solnset.push_back(sol);
        return;
    }
    for( adjstr &adj : allstrings ){
        bool dis4all = true;
        for( adjstr &tmp : sol ){
            if( !disjoint(adj, tmp) ) 
                dis4all = false;
        }
        if( dis4all ){
            getDisjointSets(solnset, allstrings, sol + adj, setsleft - 1);
        }
    }
}

vector<vector<adjstr>> getSolutions(vector<adjstr> allstrings, int solsize){
    vector<vector<adjstr>> result;
    for( adjstr &adj : allstrings ){
        vector<adjstr> sol = {adj};
        getDisjointSets(result, allstrings, sol, solsize-1);
    }
    return result;
}

void printsol(vector<adjstr> sol){
    cout << "{ " << endl;
    for( adjstr &adj : sol ){
        cout << "\t[ ";
        for( int i : adj ){
            cout << i << " ";
        }
        cout << "]" << endl;
    }
    cout << "}" << endl;
}

vector<vector<adjstr>> uniqueSolns(const vector<vector<adjstr>> &solns){
    vector<vector<adjstr>> result;
    for( auto it = solns.begin(); it < solns.end(); it++ ){
        bool unique = true;
        for( auto other = it+1; other < solns.end(); other++ ){
            bool same = true;
            for( const adjstr &adj : *it ){
                bool diff = true;
                for( const adjstr &cmp : *other ){
                    if( equal(adj, cmp) ){ 
                        diff = false;
                        break;
                    }
                }
                if( diff ){
                    same = false;
                    break;
                }
            }
            if( same ){ 
                unique = false;
                break;
            }
        }
        if( unique ){
            result.push_back(*it);
        }
    }
    return result;
}

int main(int argc, char** argv){
    int n = 4;
    if( argc > 1 )
        n = stoi(argv[1]);
    
    int* matrix = new int[n*n];
    for( int i = 1; i <= n*n; i++ ){ matrix[i-1] = i; }

    printmat(matrix, n);

    cout << "Factors of " << n*n << ": ";
    for( int i : factors(n*n)){
        cout << i << " ";
    }
    cout << endl;

    int factor = 2;
    cout << "All adjacent strings of size "<< factor <<" in matrix: " << endl;
    auto sols = getSolutions(getAllAdjStrings(matrix, n, factor), (n*n)/factor);
    auto unq = uniqueSolns(sols);
    cout << "Solution set size: " << sols.size() << endl;
    cout << "Unique set size: " << unq.size() << endl;

    //for( auto sol : unq ){
    //    printsol(sol);
    //}

    delete[] matrix;
}