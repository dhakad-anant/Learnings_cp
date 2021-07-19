#include <bits/stdc++.h>
using namespace std;
//codechef
#define rep(i,s,e) for(int i=s ; i < e ; i++)
#define rrep(i,s,e) for(int i=s ; i > e ; i--)
#define srep(i,s,e,j) for(int i=s ; i < e ; i+=j)
#define ll long long
#define vi vector<int>
#define vpi vector<pii>
#define maxpqi priority_queue<int>
#define minpqi priority_queue <int, vector<int>, greater<int> >
#define pii pair<int,int> 
#define F first
#define S second
#define mk make_pair
#define pb push_back
#define endl '\n'
#define br printf("\n")
#define inf 1000000
#define imap map<int,int>
#define uimap unordered_map<int,int>
#define gcd(a,b) __gcd(a,b)
#define max3(a,b,c) max(max((a),(b)),(c))
#define max4(a,b,c,d) max(max((a),(b)),max((c),(d)))
#define min3(a,b,c) min(min((a),(b)),(c))
#define min4(a,b,c,d) min(min((a),(b)),min((c),(d)))
#define mod 1000000007
#define clr(x) memset(x,0,sizeof(x))
#define fill(x,y) memset(x,y,sizeof(x))
// x << y    ====   x*(2^y)
// x >> y    ====   x/(2^y)

void solve(){
    int tt;
    cin >> tt;
    while(tt--){
        int n,q;
        cin >> n >> q;

        int a[n+1]; clr(a);
        int l=0,r=0;
        vpi queries;
        rep(i,0,q){
            cin >> l >> r;
            queries.pb(mk(l,r));
            a[l]++;
            if(r+1 <= n)a[r+1]--;
        }
        rep(i,2,n+1)a[i] += a[i-1];

        rep(i,0,q){
            l = queries[i].F;
            r = queries[i].S;
            if(r+1 <= n)a[r+1] -= (r - l + 1);
        }
        rep(i,2,n+1)a[i] += a[i-1];
        
        rep(i,1,n+1)cout << a[i] << " ";
        cout << endl;
    }
}

int main()
{
    // ios_base::sync_with_stdio(false);
    // cin.tie(NULL);
    // cout.tie(NULL);
    solve();
    return 0;
}
