//shortest path session by manas sir problem 6
//OR problem 1 of https://codeforces.com/blog/entry/45897

#include <bits/stdc++.h>
using namespace std;

#define pii pair<int,int>
#define INF INT_MAX
#define pb push_back
#define rep(i,s,e) for(int i=s; i<e; i++)
#define F first
#define S second

const int mxN = 1e3;
const int mxM = 1e5;
const int mxC = 1e3;    

int n,m,c,a,b;
vector<pii>adj[mxN];
int dist[mxN][mxC];
int cost[mxN];

struct State{
    int v;//vertex
    int c;//currentFuel
    int d;//distance from src node
    State(int a,int b,int _c){
        v = a,c = b,d = _c;
    }
    bool operator>(const State& s) const {return d > s.d;}
};


int dijkstra(){
    rep(i,0,n){
        rep(j,0,C+1)dist[i][j] = INF;
    }
    priority_queue<State, vector<State>, greater<State>>pq;
    pq.push(State(a,C,0));
    State cur;
    while(!pq.empty()){
        cur = pq.top();
        if(dist[cur.v][cur.c] < cur.d)continue;
        dist[cur.v][cur.c] = cur.d;

        for(auto i : adj[cur.v]){
            int v = i.F;int w = i.S;
            for(int buy=C-cur.c; buy>=0 && curr.c+buy-w>=0; buy--){
                if(dist[v][curr.c+buy-w] > s.d + buy*cost[curr.v]){
                    dist[v][curr.c+buy-w] = s.d + buy*cost[curr.v];
                    pq.push(State(v,curr.c+buy-w,dist[v][curr.c+buy-w]));
                }
            }
        }
    }
    int ans = INF;
    rep(i,0,C+1)ans = min(ans, dist[b][i]);
    return ans;
}


int main(){
    cin>>n>>m>>c>>a>>b;
    a--;b--;
    int u,v,w;
    rep(i,0,m){
        cin>>u>>v>>w;
        --u;--v;
        adj[u].pb({v,w});
        adj[v].pb({u,w});
    }

    int ans = dijkstra();
    cout << (ans==INT_MAX ? -1 : ans) << endl;
}