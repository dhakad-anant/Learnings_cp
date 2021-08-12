
// https://codeforces.com/contest/1375/problem/D

#include <bits/stdc++.h>
using namespace std;

#define fastio ios_base::sync_with_stdio(0);cin.tie(0);cout.tie(0);
#define rep(i,s,e) for(int i=s ; i < e ; i++)
#define rrep(i,s,e) for(int i=s ; i > e ; i--)
#define srep(i,s,e,j) for(int i=s ; i < e ; i+=j)
#define tr(i,x) for(auto i : x)
#define vinp(a) for(int i=0 ; i<a.size() ; i++)cin>>a[i]
#define ainp(a,n) for(int i=0; i<n; i++)cin>>a[i]
#define int long long
#define vi vector<int>
#define vs vector<string>
#define vb vector<bool>
#define vpi vector<pii>
#define maxpqi priority_queue<int>
#define minpqi priority_queue <int, vector<int>, greater<int> >
#define pii pair<int,int> 
#define F first
#define S second
#define mk make_pair
#define pb push_back
#define pf push_front
#define endl '\n'
#define gcd(a,b) __gcd(a,b)
#define clr(x) memset(x,0,sizeof(x))
#define lb lower_bound
#define ub upper_bound
#define npos string::npos
#define all(x) x.begin(),x.end()
#define sayyes cout << "YES" << endl
#define sayno cout << "NO" << endl

void solve(){
    int tt = 1;
    cin >> tt;
    while(tt--){
        int n; cin >> n;
        int a[n];ainp(a, n);
        set<int>mex;
        // unordered_map<int, int>f;
        int f[n+1]; fill(f, f+n+1, 0);  
        rep(i, 0, n+1)mex.insert(i);

        rep(i, 0, n){
            f[a[i]]++;
            if(mex.find(a[i]) != mex.end()){
                mex.erase(a[i]);
            }
        }
        bool sorted = true;
        rep(i, 1, n){
            if(a[i-1] > a[i]){
                sorted = false;
                break;
            }
        }
        vi ans;
        int cnt = 0;
        int curmex = *mex.begin();
        while(cnt < 2*n && !sorted){
            if(curmex == n){
                rep(i, 0, n){
                    if(a[i] != i){
                        ans.pb(i);
                        if(--f[a[i]] == 0)
                            mex.insert(a[i]);

                        a[i] = curmex;
                        f[curmex]++;
                        mex.erase(curmex);        
                        break;
                    }
                }
            }else{
                ans.pb(curmex);
                if(--f[a[curmex]] == 0)
                    mex.insert(a[curmex]);

                a[curmex] = curmex;
                f[curmex]++;
                mex.erase(curmex);
            }
            sorted = true;
            curmex = *mex.begin();
            rep(i, 1, n){
                if(a[i-1] > a[i]){
                    sorted = false;
                    break;
                }
            }
            cnt++;
        }
        cout << ans.size() << endl;
        tr(i, ans)cout << i+1 << " ";
        cout<<endl;
    }   
}

int32_t main()
{
    fastio
    // freopen("input.txt","r",stdin);
    // freopen("output.txt","w",stdout);
    solve();
    return 0;
}
