#include <bitset>

using namespace std;

// rewrite following function without using bitset 
vector<int> t_hop(int n_, int np_) {
    Initialize();
    vector<int> temp;
    vector<int> hop;
    int count;
    // replace following line with something equivalent
    bitset<np_> foo(n_);
    foo.to_string();
    for (int i = 0; i < N_s; ++i) {
        temp.clear();
        count = 0;
        for (int j = 0; j < N; ++j) {
            temp.push_back(foo.test(i+j*N_s));
            count += int(foo.test(i+j*N_s));
        }
        if (count > 0 && count < N) {
            bitset<np_> state_(foo);
            vector<int> temp_(temp);
            sort(temp.begin(),temp.end());
            do {
                if (temp != temp_) {
                    for (int k = 0; k < N; ++k) {
                        state_.set(i+k*N_s, temp[k]);
                    }
                    hop.push_back((int) (state_.to_ulong()));
                }
            } while ( next_permutation(temp.begin(), temp.end()) );
        }
    }
    return hop;
}