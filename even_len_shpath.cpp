//OR problem 2 of https://codeforces.com/blog/entry/45897

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

int n,m,a,b;
vector<pii>adj[mxN];
int dist[mxN][mxC];

struct State{
    int v;//vertex
    int p;//parity of currNode;
    int d;//distance from src node;
    State(int a,int b,int _c){
        v = a,p = b,d = _c;
    }
    bool operator>(const State& s) const {return d > s.d;}
};


int dijkstra(){
    rep(i,0,n)dist[i][0] = dist[i][1] = INF;

    priority_queue<State, vector<State>, greater<State>>pq;
    pq.push(State(a,0,0));
    dist[a][0]=0;
    State cur;
    while(!pq.empty()){
        cur = pq.top();pq.pop();
        if(dist[cur.v][cur.p] < cur.d)continue;

        for(auto i : adj[cur.v]){
            int v = i.F;int w = i.S;
            if(dist[v][cur.p^1] > dist[cur.v][cur.p] + w){
                dist[v][cur.p^1] = dist[cur.v][cur.p] + w;
                pq.push(State(v,cur.p^1,dist[v][cur.p^1]));
            }
        }
        
    }
    return dist[b][0];
}


int main(){
    cin>>n>>m>>a>>b;
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