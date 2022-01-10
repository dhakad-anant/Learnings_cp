#include <bits/stdc++.h>
using namespace std;

// https://www.geeksforgeeks.org/ordered-set-gnu-c-pbds/
#define rep(i,s,e) for(int i=s ; i < e ; i++)
#define rrep(i,s,e) for(int i=s ; i > e ; i--)
#define srep(i,s,e,j) for(int i=s ; i < e ; i+=j)
#define tr(i,x) for(auto i : x)
#define int long long
// #define ll long long
#define vi vector<int>
#define vll vector<ll>
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
#define sumv(a) accumulate(a.begin(),a.end(),0)
#define suma(a,n) accumulate(a,a+n,0)
#define maxv(a) max_element(a.begin(),a.end()) // gives pointer to max value 
#define maxa(a,n) max_element(a,a+n) // gives pointer to max value
#define max3(a,b,c) max(max((a),(b)),(c))
#define max4(a,b,c,d) max(max((a),(b)),max((c),(d)))
#define min3(a,b,c) min(min((a),(b)),(c))
#define min4(a,b,c,d) min(min((a),(b)),min((c),(d)))
#define mod 1000000007
#define clr(x) memset(x,0,sizeof(x))
#define fill(x,y) memset(x,y,sizeof(x))
#define lb lower_bound
#define ub upper_bound
#define all(x) x.begin(),x.end()
#define descen greater<int>()
#define si1(a) scanf("%d",&a);
#define si2(a,b) scanf("%d %d",&a,&b);
#define pi1(a) printf("%d ",a);
#define pi2(a,b) printf("%d %d ",a,b)
#define pl1(a) printf("%lld ",a)
#define pt printf

//bitset<n>(x).to_string() -- return string // here n is bit length.
// x << y    ====   x*(2^y)
// x >> y    ====   x/(2^y)
// x&1 --- 1 then x is odd

// Driver function to sort the vector elements 
// by second element of pairs 
// -- IN ASCENDING ORDER --
bool sortbysec(const pair<int,int> &a, const pair<int,int> &b) { 
    return (a.second < b.second); 
}

// Returns value of Binomial Coefficient C(n, k)
int binomialCoeff(int n, int k){
    int res = 1;
    // Since C(n, k) = C(n, n-k)
    if (k > n - k)
        k = n - k;
    // Calculate value of [n*(n-1)*---*(n-k+1)] /
    // [k*(k-1)*---*1]
    for (int i = 0; i < k; ++i) {
        res *= (n - i);
        res /= (i + 1);
    }
    return res;
}

const int N1 = 1e7+7;
int primefactcnt[N1]; // if primefactcnt[i]==0  ==> i is a prime number.
void cntprimefact(){
    fill(primefactcnt,0);
    
    for(int p=2; p<N1 ; p++){
        if(primefactcnt[p] != 0)continue;
        else primefactcnt[p] = 1;
        // allow if p is a prime number
        for(int i=2*p; i<N1 ; i+=p)primefactcnt[i]++;
    }
}

// ----------------- string matching algorithm kmp starts here -----------------
vector<int>prefix_function(string s){
    int n = s.size();
    vector<int>p(n, 0);
    int j;
    for(int i=1; i<n; i++){
        j = p[i-1];
        while(j > 0 && s[i] != s[j])j = p[j-1];

        if(s[i] == s[j])j++;
        p[i] = j;
    }
    return p;
}
int KMPSearch(string txt, string pat){ // return number of times pat is found in txt;
    vector<int>lps = prefix_function(pat);
    int n = txt.length();
    int m = pat.length();
    int cnt = 0;
    int i=0, j=0;
    while(i < n){
        if(pat[j] == txt[i]){
            i++;j++;
        }
        if(j == m){
            cnt++;
            j = lps[j-1];
            continue;
        }
        if(j>0 && pat[j] != txt[i]){
            j = lps[j-1];
        }else if(pat[j] != txt[i]){
            i++;
        }
    }
    return cnt;
}
// ----------------- string matching algorithm kmp ends here -----------------


// ----------------- Generating all subsets using recursion. .//BEGINS -----------------
vector<int>subset;
void genAllSubsets(int a[], int idx, int n){
    if(idx == n){
        // Print current subset.
        for(auto x : subset)cout << x << " ";
        cout<<endl;
    }else{
        // Include a[idx]
        subset.push_back(a[idx]);
        genAllSubsets(a, idx+1, n);

        // Exclude a[idx]
        subset.pop_back();
        genAllSubsets(a, idx+1, n);
    }
}
// ----------------- Generating all subsets using recursion. .//ENDS -----------------



bool binarySearch(vector<int>a,int x){
    int n = a.size();
    int lo = 0;
    int hi = n-1;
    int mid;
    while(lo <= h){
        mid = lo + (hi - lo)/2;
        if(a[mid] == x)return true;
        else if(a[mid] < x){
            lo = mid+1;
        }else if(a[mid] > x){
            hi = mid-1;
        }
    }
    return false;
}
//Longest common subsequence
int LCS(string s1,string s2){
    int n=s1.length(),m=s2.length();
    if(n > m){
        swap(n,m);
        string tmp = s1;
        s1 = s2;
        s2 = tmp;
    }
    int dp[n+1][m+1];
    clr(dp);
    rep(i,1,n+1){
        rep(j,1,m+1){
            if(s1[i-1] == s2[j-1]){
                dp[i][j] = 1 + dp[i-1][j-1];
            }else{
                dp[i][j] = max3(dp[i-1][j-1], dp[i][j-1], dp[i-1][j]);
            }
        }
    }
    return dp[n][m];
}
//LONGEST INCREASING SUBSEQUENCE
//nlogn approach.
int LIS(vi arr){
    int n = arr.size();
    vi dp(n+1, INT_MAX);
    rep(i,0,n){
        *lower_bound(all(dp) , arr[i]) = arr[i];
    }
    rep(i,0,n+1){
        if(dp[i] == INT_MAX)return i;
    }
}
int Froot(int parent[],int i){
    while( i != parent[i] ){
        parent[i] = parent[parent[i]];
        i = parent[i];
    }
    return i;
}
bool U(int sz[],int parent[],int x,int y){
    int rx = Froot(parent,x);
    int ry = Froot(parent,y);
    if(rx == ry)return false;
    if(sz[rx] > sz[ry]){
        sz[rx] += sz[ry];
        parent[ry] = rx;
    }else{
        sz[ry] += sz[rx];
        parent[rx] = ry;
    }
    return true;
}

bool isPrime(int n) { 
    // Corner cases 
    if (n <= 1) return false; 
    if (n <= 3) return true; 
  
    // This is checked so that we can skip 
    // middle five numbers in below loop 
    if (n % 2 == 0 || n % 3 == 0) return false; 
  
    for (int i = 5; i * i <= n; i = i + 6){
        if (n % i == 0 || n % (i + 2) == 0) return false;
    } 
    return true; 
} 

void Factorization(int n,imap & pfactorfreq){  
    // Print the number of 2s that divide n  
    int c = 0;
    while (n % 2 == 0){ 
        c++;
        n = n/2;  
    }  
    if(c > 0)pfactorfreq[2] = c;
    // n must be odd at this point. So we can skip  
    // one element (Note i = i +2)  
    for (int i = 3; i <= sqrt(n); i = i + 2)  
    {  
        // While i divides n, print i and divide n  
        c = 0;
        while (n % i == 0){  
            c++;
            n = n/i;  
        }  
        if(c > 0)pfactorfreq[i] = c;
    }  
    // This condition is to handle the case when n  
    // is a prime number greater than 2  
    if (n > 2){
        pfactorfreq[n] += 1;
    }
}
const int N = 1e5;
bool isprime[N];
void seive(){
    fill(isprime,true); isprime[0]=isprime[1]=false;

    for(int p=2 ; p*p<=N ; p++){
        if(isprime[p]){
            for(int i=p*p ; i<=N ; i+=p)isprime[i]=false;
        }
    }
}
// inv(a) = a^{mod-2} // useful formula.
int mod_pow(int x,int y,int p){ // (x^y)%p
    int res = 1;
    x = x%p;
    if(x == 0 && y != 0)return 0;
    while(y > 0){
        // if y is odd multiply res with x.
        if( y&1 )
            res = (res*x)%p;
        y = y>>1; // y = y/2;
        x = (x*x)%p;
    }   
    return res;
}   


//list of NEXT GREATER ELEMENT IN O(N)
vi next_greater_element(vi a){
    int n = a.size();
    vi ans(n);
    stack<pii>stk;
    stk.push({a[0],0});
    int nxt=0;
    rep(i,1,n){
        nxt = a[i];
        while(!stk.empty() && stk.top().F <= nxt){
            ans[stk.top().S] = i;
            stk.pop();
        }
        stk.push({nxt,i});
    }
    while(!stk.empty()){
        ans[stk.top().S] = -1;
        stk.pop();
    }
    return ans;
}

vi optimized_next_greater_element(vi a){ // here a is 1-index based (a1, a2, a3...)
    int n = a.size();
    vi r(n);
    rrep(i, n-1, 0){
        r[i] = i;
        while(r[i] < n-1 && a[i] >= a[r[i]+1]){
            r[i] = r[r[i]+1];
        }
    }
    return r;
}

vi nearest_smaller_num_onLeft(vi a){
    int n = a.size();
    vi ans(n);
    stack<int>st;
    rep(i,0,n){
        while(!st.empty() && st.top() >= a[i])st.pop();
        if(st.empty()){
            ans[i] = -1;
        }else{
            ans[i] = st.top();
        }
        st.push(a[i]);
    }
    return ans;
}

int SumOfXorOfAllPairs(int a[],int n,int numbits){
    int sum = 0;
    int cnt[2];
    rep(i,0,numbits){
        clr(cnt);
        rep(j,0,n)cnt[(a[i]>>i)&1]++;
        sum += (1LL<<i)%mod*cnt[0]%mod*cnt[1];
        sum = sum%mod;
    }
    return sum;
}

// Dijkstra with complexity of O(ElogV)
vi dijkstra(vpi adj[],int src,int n){
	set<pii>finalized;
	vi dist(n,INT_MAX);
	dist[src] = 0;
	finalized.insert({0,src});
	pii curr;
	int u,v,w;
	while(!finalized.empty()){
		curr = *finalized.begin();
		finalized.erase(finalized.begin());

		u = curr.S;
		tr(i, adj[u]){
			v = i.F;
			w = i.S;
			if(dist[v] > dist[u] + w){

				if(dist[v] != INT_MAX)finalized.erase(finalized.find({dist[v],v}));

				dist[v] = dist[u] + w;
				finalized.insert({dist[v], v});
			}
		}
	}
	return dist;
}

/* REMOVING CONSECUTIVE DUPLICATES ELEMENT FROM THE LIST.
'''
The idea is to keep track of two indexes, index of current character in str and 
index of next distinct character in str.
Whenever we see same character, we only increment current character index. 
We see different character, we increment index of distinct character.
'''
        int j = 0; // lastdistinct index

        loop(i,1,n){
            if(a[j] != a[i]){
                j++;
                a[j] = a[i];
            }
        }
        int lastINDEX = j; // USEFUL LIST IS TILL lastINDEX.

        -- THIS IS OPTIONAL -- 
        j++;
        for(; j < n ; j++){
            a[j] = 0;
        }

''' -- stl way -- '''
        vector<int> a(n);
		for (auto &it : a) cin >> it;  // take input.
		n = unique(a.begin(), a.end()) - a.begin(); // pointer to new end of list.
		a.resize(n);
*/


/*  ////////////FINDING LARGEST FACTOR OF G LESS THAN N.
for(int i=1 ; i*i <= G ; i++){
    if(G%i == 0){
        int t1 = i,t2 = G/i;
        if(n >= t1)ans = max(ans,t1);
        if(n >= t2)ans = max(ans,t2);
    }
}
*/

// printf("%.10f\n",ans/2.0);

//------------VERY LARGE INPUT 
// __int128

// a&(a-1)  == 0 // if a is a power of 2.

// ------------LINKED LIST USEFUL 
// Custom list node.
template <typename T>
class ListNode {
public:
	T value;
	ListNode<T>* next;
	ListNode(T value) {
		this->value = value;
		this->next = nullptr;
	}
};
ListNode<int> *reverseList(ListNode<int> *head){
    ListNode<int> *prev, *cur;
    prev = nullptr;
    cur = head;

    while (cur)
    {
        ListNode<int> *next = cur->next;
        cur->next = prev;
        prev = cur;
        cur = next;
    }

    return prev;
}
ListNode<int>* findMid(ListNode<int> *head){
    ListNode<int> *slow, *fast;
    slow = head;
    fast = head;

    while (fast and fast->next){
        slow = slow->next;
        fast = fast->next->next;
    } 
    return slow;
}

// ------------LINKED LIST USEFUL end --------------------------------


// ----------------- TRIE IMPLE USING 2D-ARRAY -----------------
const int MXN = 1e5;
const int ALPH = 26;
class Trie(){
    int adj[MXN][ALPH];
    bool hotnode[MXN];
    int nodecnt;

    Trie(){
        memset(adj,-1,sizeof(adj));
        memset(hotnode,false,sizeof(hotnode));
        nodecnt = 0;
    }
    void insert(const string s){
        int n = s.length();
        int cur = 0;
        for(int i=0; i<n; i++){
            if(adj[cur][s[i]-'a'] == -1)adj[cur][s[i]-'a'] = ++nodecnt;
            cur = adj[cur][s[i]-'a'];
        }   
        hotnode[cur] = true;
    }
    bool ExactSearch(const string s){
        int n = s.length();
        int cur = 0;
        for(int i=0; i<n; i++){
            if(adj[cur][s[i]-'a'] == -1)return false;
            cur = adj[cur][s[i]-'a'];
        }
        return hotnode[cur];
    }
};
// ----------------- TRIE IMPLE USING 2D-ARRAY END HERE -----------------


// ----------------- modular arithmetic and much more -----------------
const int MOD = 1e9 + 7;
struct mod_int {
    int val;
 
    mod_int(long long v = 0) {
        if (v < 0)
            v = v % MOD + MOD;
 
        if (v >= MOD)
            v %= MOD;
 
        val = v;
    }
 
    static int mod_inv(int a, int m = MOD) {
        int g = m, r = a, x = 0, y = 1;
 
        while (r != 0) {
            int q = g / r;
            g %= r; swap(g, r);
            x -= q * y; swap(x, y);
        }
        return x < 0 ? x + m : x;
    }
 
    explicit operator int() const {
        return val;
    }
 
    mod_int& operator+=(const mod_int &other) {
        val += other.val;
        if (val >= MOD) val -= MOD;
        return *this;
    }
 
    mod_int& operator-=(const mod_int &other) {
        val -= other.val;
        if (val < 0) val += MOD;
        return *this;
    }
 
    static unsigned fast_mod(uint64_t x, unsigned m = MOD) {
           #if !defined(_WIN32) || defined(_WIN64)
                return x % m;
           #endif
           unsigned x_high = x >> 32, x_low = (unsigned) x;
           unsigned quot, rem;
           asm("divl %4\n"
            : "=a" (quot), "=d" (rem)
            : "d" (x_high), "a" (x_low), "r" (m));
           return rem;
    }
 
    mod_int& operator*=(const mod_int &other) {
        val = fast_mod((uint64_t) val * other.val);
        return *this;
    }
 
    mod_int& operator/=(const mod_int &other) {
        return *this *= other.inv();
    }
 
    friend mod_int operator+(const mod_int &a, const mod_int &b) { return mod_int(a) += b; }
    friend mod_int operator-(const mod_int &a, const mod_int &b) { return mod_int(a) -= b; }
    friend mod_int operator*(const mod_int &a, const mod_int &b) { return mod_int(a) *= b; }
    friend mod_int operator/(const mod_int &a, const mod_int &b) { return mod_int(a) /= b; }
 
    mod_int& operator++() {
        val = val == MOD - 1 ? 0 : val + 1;
        return *this;
    }
 
    mod_int& operator--() {
        val = val == 0 ? MOD - 1 : val - 1;
        return *this;
    }
 
    mod_int operator++(int32_t) { mod_int before = *this; ++*this; return before; }
    mod_int operator--(int32_t) { mod_int before = *this; --*this; return before; }
 
    mod_int operator-() const {
        return val == 0 ? 0 : MOD - val;
    }
 
    bool operator==(const mod_int &other) const { return val == other.val; }
    bool operator!=(const mod_int &other) const { return val != other.val; }
 
    mod_int inv() const {
        return mod_inv(val);
    }
 
    mod_int pow(long long p) const {
        assert(p >= 0);
        mod_int a = *this, result = 1;
 
        while (p > 0) {
            if (p & 1)
                result *= a;
 
            a *= a;
            p >>= 1;
        }
 
        return result;
    }
 
    friend ostream& operator<<(ostream &stream, const mod_int &m) {
        return stream << m.val;
    }
    friend istream& operator >> (istream &stream, mod_int &m) {
        return stream>>m.val;   
    }
};
// ----------------- modular arithmetic and much more end here -----------------


//----------------------------------------------------------------floyd warshall
void floyd(vector<vector<int>>&dp){
    // i->i is zero
    // i->j if there is no such edge then dp[i, j] = INT_MAX;
    int n = adj.size();
    
    //selecting the intermediate set as {0, 1, 2, ..k}
    for(int k=0; k<n; k++){ 
        //selecting the source vertex i 
        for(int i=0; i<n; i++){
            //selecting the source vertex j
            for(int j=0; j<n; j++){
                if(dp[i][k] != INT_MAX && dp[k][j] != INT_MAX){
                    dp[i][j] = min(dp[i][j], dp[i][k]+dp[k][j]);
                }
            }
        }
    }
}
//------------------------------------------------------------floyd warshall end

bool IsCyclePresent(int src){
    vis[src] = 1; // Grey mean still processing
    for(auto v : adj[src]){
        if(vis[v] == 1)return true; // backedge
        if(vis[v] == 0 && IsCyclePresent(v))return true; // cycle in subtree
    }
    vis[src] = 2; // Black means subtree processed
    return false;
}

//------------------------------- BINARY LIFTING // nodes are zero indexed; -------------------------------
const int mxN = 1e5;
const int LOG = ceil(1.0*log2(mxN));
vector<vector<int>>up(mxN, vector<int>(LOG));
vector<int>depth(mxN);
void PreComputationBinaryLifting(int u, int p){
    up[u][0] = p;
    depth[u] = depth[p] + 1;
    for(int j=1; j<LOG; j++){
        up[u][j] = up[ up[u][j-1] ][j-1];
    }
    for(auto v : g[u]){
        if(v == p)continue;
        PreComputationBinaryLifting(v, u);
    }
}

int LCA(int u, int v){
    if(depth[u] < depth[v])swap(u, v);
    // 1) get the same depth 
    int k = depth[u] - depth[v];
    for(int i=LOG-1; i>=0; i--){
        if(k & (1 << i))
            u = up[u][i];
        // if(depth[u] - depth[v] >= (1 << i))
        //     u = up[u][i];
    }
    // 2) if v was ancestor of u then now u==v
    if(u == v)return v;

    // 3) move both u and v with powers of two
    for(int i=LOG-1; i>=0; i--){
        if(up[u][i] != up[v][i]){
            u = up[u][i];
            v = up[v][i];
        }
    }
    return up[u][0];
}
//------------------------------- BINARY LIFTING ends -------------------------------

//------------------------------- SEGMENT TREE TEMPLATES ------------------------------- 
class LazySegTree{
public:
    /* NOTATIONS
        si - Index of current node in segment tree 
             Initial value is passed as 0.
        [ss, se] - Starting and ending indexes of array elements 
                   covered under this node of segment tree.
                   Initial values passed as 0 and n-1.
        [qs, qe] - Starting and ending indexes of query.
        [us, ue] - Starting and ending indexes of update query.
    */
    int n;
    vector<int> tree;
    vector<int> lazy;
    int lazy_init = 0;

    //initializing segment tree
    LazySegTree(int x){
        n = x;
        tree.resize(4*n + 5, 0);
        lazy.resize(4*n + 5, lazy_init);
    }
    
    void build(int si, int ss, int se, vector<int>&a){
        if(ss == se){
            tree[si] = a[ss];
            return;
        }
        int left = 2*si;
        int right= left+1;
        int mid = ss + (se-ss)/2;
        build(left, ss, mid, a);
        build(right, mid+1, se, a);
        // tree[si] = tree[left] + tree[right];
        tree[si] = min(tree[left], tree[right]);
    }

    void push_down(int si, int ss, int se){
        if(lazy[si] != lazy_init){
            tree[si] += (se-ss+1)*lazy[si]; // in case of sum segment tree
            // tree[si] += lazy[si];
            if(ss != se){
                lazy[2*si] += lazy[si];
                lazy[2*si + 1] += lazy[si];
            }
            lazy[si] = lazy_init;
        }
    }

    void update(int us, int ue, int val){
        _update(1, 0, n-1, us, ue, val);
    }
    void _update(int si, int ss, int se, int us, int ue, int val){
        // first make any pending updates.
        push_down(si, ss, se);

        //no intersection
        if(ss > se || ss > ue || se < us)return;

        //complete overlap
        if(us <= ss && se <= ue){
            lazy[si] += val;
            push_down(si, ss, se);
            return;
        }
        //some overlap
        int left = 2*si;
        int right= left+1;
        int mid = ss + (se-ss)/2;
        _update(left, ss, mid, us, ue, val);
        _update(right, mid+1, se, us, ue, val);
        // tree[si] = tree[left] + tree[right];
        tree[si] = min(tree[left], tree[right]);
    }

    int query(int qs, int qe){
        return _q(1, 0, n-1, qs, qe);
    }
    int _q(int si, int ss, int se, int qs, int qe){
        // first make any pending updates.
        push_down(si, ss, se);

        //no intersection
        if(ss > se || ss > qe || se < qs)return 0;

        //complete overlap
        if(qs <= ss && se <= qe){
            return tree[si];
        }
        //some overlap
        int left = 2*si;
        int right= left+1;
        int mid = ss + (se-ss)/2;
        int q1 = _q(left, ss, mid, qs, qe);
        int q2 = _q(right, mid+1, se, qs, qe);
        // return q1 + q2;
        return min(q1, q2);
    }
};

class SegTree{
public:
    /* NOTATIONS
        si - Index of current node in segment tree 
             Initial value is passed as 0.
        [ss, se] - Starting and ending indexes of array elements 
                   covered under this node of segment tree.
                   Initial values passed as 0 and n-1.
        [qs, qe] - Starting and ending indexes of query.
        [us, ue] - Starting and ending indexes of update query.
    */
    int n;
    vector<int> tree;

    //initializing segment tree
    SegTree(int x){
        n = x;
        tree.resize(4*n + 5, 0);
    }

    void build(int si, int ss, int se, vector<int>&a){
        if(ss == se){
            tree[si] = a[ss];
            return;
        }
        int left = 2*si;
        int right= left+1;
        int mid = ss + (se-ss)/2;
        build(left, ss, mid, a);
        build(right, mid+1, se, a);
        // tree[si] = tree[left] + tree[right];
        tree[si] = min(tree[left], tree[right]);
    }

    void update(int us, int ue, int val){
        _update(1, 0, n-1, us, ue, val);
    }
    void _update(int si, int ss, int se, int us, int ue, int val){
        //no intersection
        if(ss > se || ss > ue || se < us)return;

        //complete overlap
        if(us <= ss && se <= ue){
            tree[si] = val;
            return;
        }
        //some overlap
        int left = 2*si;
        int right= left+1;
        int mid = ss + (se-ss)/2;
        _update(left, ss, mid, us, ue, val);
        _update(right, mid+1, se, us, ue, val);
        // tree[si] = tree[left] + tree[right];
        tree[si] = min(tree[left], tree[right]);
    }

    int query(int qs, int qe){
        return _q(1, 0, n-1, qs, qe);
    }
    int _q(int si, int ss, int se, int qs, int qe){
        //no intersection
        if(ss > se || ss > qe || se < qs)return 0;

        //complete overlap
        if(qs <= ss && se <= qe){
            return tree[si];
        }
        //some overlap
        int left = 2*si;
        int right= left+1;
        int mid = ss + (se-ss)/2;
        int q1 = _q(left, ss, mid, qs, qe);
        int q2 = _q(right, mid+1, se, qs, qe);
        // return q1 + q2;
        return min(q1, q2);
    }
};
//------------------------------- SEGMENT TREE TEMPLATES ENDS------------------------------- 

class DSU{
public:
    // 1 indexed;
    int n;
    vector<int>p;
    // vector<int>sz;
    DSU(int _n){
        n = _n;
        p.resize(n+1, -1);
        // p.resize(n+1);iota(all(p), 0);
        // sz.resize(n+1, 1);
    }
    int leader(int i){
        if(p[i] < 0)return i;
        // if(p[i] == i)return i;
        return p[i] = leader(p[i]);
    }
    void merge(int pu, int pv){
        if(p[pu] > p[pv])swap(pu, pv);
        // pu is a bigger component;
        p[pu] += p[pv];
        p[pv] = pu;

        // if(sz[pu] < sz[pv])swap(pu, pv);
        // // pu is a bigger component;
        // sz[pu] += sz[pv];
        // p[pv] = pu;
    }
};

int main(){
    int n; cin >> n;
    imap pf;
    primeFactors(n,pf);
    tr(i,pf)cout << i.F << " ";
    cout << endl;
    return 0;
}


//nCk % MOD
//https://ideone.com/sctrLM

/* FOR MEX QUERIES
 https://codeforces.com/blog/entry/81287?#comment-677837
*/