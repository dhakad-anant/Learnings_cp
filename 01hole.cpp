#include <bits/stdc++.h>
using namespace std;

#define fastio ios_base::sync_with_stdio(0);cin.tie(0);cout.tie(0);
#define rep(i,s,e) for(int i=s ; i < e ; i++)
#define rrep(i,s,e) for(int i=s ; i > e ; i--)
#define srep(i,s,e,j) for(int i=s ; i < e ; i+=j)
#define tr(i,x) for(auto i : x)
#define inp(a) for(int i=0 ; i<a.size() ; i++)cin>>a[i]
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
#define lb lower_bound
#define ub upper_bound
#define npos string::npos
#define all(x) x.begin(),x.end()

void solve(){
    int n;cin>>n;
    vi a(n,0);inp(a);
    vi le(n,0),ri(n,0);
    le[0]=a[0];
    rep(i,1,n)le[i] = max(le[i-1]+a[i], a[i]);

    ri[n-1]=a[n-1];
    rrep(i,n-2,-1)ri[i] = max(ri[i+1]+a[i], a[i]);

    int fans=a[0];
    int ans=a[0];
    rep(i,1,n){
        ans = max(ans+a[i], a[i]);
        fans = max(ans, fans);
    }
    // cout << "fans " << fans << endl;
    if(n==1){
        cout << fans << endl;
        return;
    }
    rep(i,0,n){
        if(i==0){
            fans = max(fans, ri[i+1]);
        }else if(i==n-1){
            fans = max(fans, le[i-1]);
        }else{
            fans = max(fans, le[i-1]+ri[i+1]);
        }
    }
    cout << fans << endl;    
}

int32_t main()
{
    fastio
    // ios::sync_with_stdio(false);
    // freopen("input.txt","r",stdin);
    // freopen("output.txt","w",stdout);
    solve();
    return 0;
}

/*
// Sort the array in descending order 
    sort(arr, arr + n, greater<int>()); 
*/