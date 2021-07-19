#include <bits/stdc++.h>
using namespace std;

// https://www.geeksforgeeks.org/ordered-set-gnu-c-pbds/
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
#define fill(x,y) memset(x,y,sizeof(x))
#define lb lower_bound
#define ub upper_bound
#define npos string::npos
#define all(x) x.begin(),x.end()

const int mxn = 20;
int mat[mxn][mxn];
int n;
int dr[] = { 2, 1, -1, -2, -2, -1, 1, 2 };
int dc[] = { 1, 2, 2, 1, -1, -2, -2, -1 };
int movesleft;

bool path(int x,int y){
    mat[x][y] = n*n - movesleft;
    movesleft--;
    if(movesleft == 0){
        return true;
    }
    int r,c;
    rep(i,0,8){
        r = x + dr[i];
        c = y + dc[i];
        if(0 <= r && r < n && 0 <= c && c < n){
            if(mat[r][c] == -1 && path(r,c))return true;
        }
    }
    mat[x][y] = -1;
    movesleft++;
    return false;
}

void solve(){
    int tt = 1;
    // cin >> tt;
    while(tt--){
        cin >> n;
        movesleft = n*n;
        memset(mat,-1,sizeof(mat));
        bool ans = false;
        rep(i,0,n){
            rep(j,0,n){
                if(path(i,j)){
                    ans = true;
                    break;
                }
            }
            if(ans)break;
        }
        set<int>st;
        rep(i,0,n){
            rep(j,0,n){
                cout<<mat[i][j] << " ";
                st.insert(mat[i][j]);
            }
            cout<<endl;
        }
        if(st.size() == n*n){
            cout << "VALID" << endl;
        }else{
            cout << "INVALID" << endl;
        }
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
