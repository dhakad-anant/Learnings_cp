/**
 * @file ExploreAllpath_UndirGraph.cpp
 * @author Anant Dhakad
 * @brief 
 * @version 0.1
 * @date 2021-12-26
 * @link https://www.geeksforgeeks.org/shortest-distance-between-given-nodes-in-a-bidirectional-weighted-graph-by-removing-any-k-edges/
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include <bits/stdc++.h>
using namespace std;

#define F first
#define S second

// To store given graph.
vector<vector<pair<int,int>>> graph; 

// To keep track of visited nodes.
vector<bool> vis; 

// To store all the edges for every path traversed during the DFS call.
vector<vector<int>> all_paths; 

// n - number of nodes in the graph.
int n;

// e - number of edges in the graph.
int e;

// k - number of edges whose weight can be reduced to zero.
int k;

/**
 * @brief DFS to find all possible paths from source to destination nodes.
 * 
 * @param src Source Node. 
 * @param dest Destination Node.
 * @param temp_path vector to store edge weights of current path
 */
void DFS_findallpath(int src, int dest, vector<int> &temp_path){

    // Reached to destination that means we found 
    // a path.
    if(src == dest){
        all_paths.push_back(temp_path);
        return;
    }

    // Marking current node as visited.
    vis[src] = true;

    // Iterate over all connections of SRC node.
    for(auto p : graph[src]){
        int v = p.F, w = p.S;

        // If this node is node visited earlier.
        if(!vis[v]){
            // Push the weight of src-->v edge in current path.
            temp_path.push_back(w);

            // Calling DFS recursively.
            DFS_findallpath(v, dest, temp_path);

            // Pop the last added edge weight so as to 
            // explore all paths.
            temp_path.pop_back();
        }
    }

    // Marking current node as Unvisited for exploring 
    // more possible paths.
    vis[src] = false;
}

/**
 * @brief Get the Minimum Distance from source to destination
 * after reducing weights of atmost K edges to zero.
 * 
 * @return int 
 */
int getMinimumDistance(){
    // To store the minimum distance.
    int minDist = INT_MAX;

    // If all_paths is empty, 
    // i.e. no path exists btw src and dest
    if(all_paths.size() == 0){
        return -1;
    }

    /**
     * @brief Iterate over all stored path to 
     * find the minimum distance path.
     */
    for(auto path: all_paths){
        
        // If some path has less than equal to K edges
        // we can reduce them all to Zero.
        if(path.size() <= k){
            return 0;
        }

        // Sorting the edge weights in decreasing order.
        sort(path.rbegin(), path.rend());

        int totalWeight=0, Ksum=0;
        for(int i=0; i<path.size(); i++){
            totalWeight += path[i];
            if(i < k)Ksum += path[i];
        }

        // Update the minimum weighted path.
        minDist = min(totalWeight-Ksum, minDist);
    } 
    return minDist;
}

void solve(){
    cin >> n >> e >> k;
    graph.resize(n);
    vis.resize(n, false);

    for(int i=0; i<e; i++){
        int u,v,w; cin >> u >> v >> w;
        graph[u].push_back({v, w});
        graph[v].push_back({u, w});
    }

    /* Source and destination nodes. */
    int A, B; cin >> A >> B;

    // For storing the edge of a particular path.
    vector<int> temp_path;
    DFS_findallpath(A, B, temp_path);

    cout << "Minimum distance btw "<< A << " & " << B << " is: " << getMinimumDistance() << endl;

    graph.clear();
    vis.clear();
    all_paths.clear();
}

int main(){
    int tt;cin >> tt;
    while(tt--){
        solve();
    }
    return 0;
}
