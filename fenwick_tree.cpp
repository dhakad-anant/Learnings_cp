#include <bits/stdc++.h>
using namespace std;

#define MAXSIZE 10000

/*--------------- BINARY INDEXED TREE. //BEGINS -----------*/
vector<int>BIT;

/**
 * @brief This is a utility function which calculate 
 * the sum of array[0 ... idx] using BIT data structure.
 * 
 * @param idx Index tree node.
 * @return int 
 */
int Sum(int idx){
    int _sum = 0;

    /* Traversing the ancestors of idx node. */
    while(idx > 0){
        _sum += BIT[idx]; // Add current element in count.

        idx -= (idx & -idx); // Move to 'sum' parent node.
    }
    // Return count of smaller element on right side.
    return _sum;
}

/**
 * @brief This function updates a node of Binary Indexed Tree(BIT) 
 * at given index 'idx'. Given value 'change' is added to BIT[idx] node 
 * and its ancestors.
 * 
 * @param idx Index tree node.
 * @param change Effective difference that needs to be updated.
 * @return void 
 */
void Update(int idx, int change){

    /* Traverse over all ancestors of BIT[idx] and add 'change' */
    while(idx < BIT.size()){    
        BIT[idx] += change; // Add 'change' in current node.

        idx += (idx & -idx); // Move to 'update' parent node.
    }
}
/*--------------- BINARY INDEXED TREE. //ENDS   -----------*/

vector<int> tree(4*MAXSIZE + 1);

int sum(int i){
    int ans = 0;
    while(i > 0){
        ans += tree[i];
        i -= (i & -i);
    }
    return ans;
}

void update(int i, int k){
    while(i < tree.size()){
        tree[i] += k;
        i += (i & -i); // subtract right side bit.
    }
}

// Building fenwick in O(N);
void buildTree(vector<int>&arr){
    int n = arr.size();
    for(int i=0; i<n; i++)tree[i] = arr[i];

    for(int i=1; i< tree.size(); i++){
        int p = i + (i & -i); // parent node
        if(p[i] < tree.size()){
            tree[p] = tree[p] + tree[i];
        } 
    }
    return tree;
}

int main(){
    // ios_base::sync_with_stdio(0);cin.tie(0);cout.tie(0);
    // freopen("input.txt","r",stdin);
    // freopen("output.txt","w",stdout);
    solve();
    return 0;
}
