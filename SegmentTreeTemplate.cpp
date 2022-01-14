/**
 * @author Anant Dhakad
 * @brief You can use this template of segment tree in coding contest.
 * @version 0.1
 * @date 2022-01-14
 * 
 */
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
        int left = si << 1;
        int right = si << 1 | 1;
        int mid = ss + ((se-ss) >> 1);
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
        int left = si << 1;
        int right = si << 1 | 1;
        int mid = ss + ((se-ss) >> 1);
        _update(left, ss, mid, us, ue, val);
        _update(right, mid+1, se, us, ue, val);
        // tree[si] = tree[left] + tree[right];
        tree[si] = min(tree[left], tree[right]);
    }

    int query(int qs, int qe){
        return _query(1, 0, n-1, qs, qe);
    }

    int _query(int si, int ss, int se, int qs, int qe){
        //no intersection
        if(ss > se || ss > qe || se < qs)return 0;

        //complete overlap
        if(qs <= ss && se <= qe){
            return tree[si];
        }
        //some overlap
        int left = si << 1;
        int right = si << 1 | 1;
        int mid = ss + ((se-ss) >> 1);
        int q1 = _query(left, ss, mid, qs, qe);
        int q2 = _query(right, mid+1, se, qs, qe);
        // return q1 + q2;
        return min(q1, q2);
    }
};