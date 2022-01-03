/**
 * @file CntLisSegtree.cpp
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2022-01-03
 * @link https://www.geeksforgeeks.org/length-of-longest-increasing-subsequences-lis-using-segment-tree/
 * @copyright Copyright (c) 2022
 */
#include <bits/stdc++.h>
using namespace std;

#define MAXSIZE 100000
#define pii pair<int, int>
#define F first
#define S second

/*-------------------- Auxillary function & objects for segment tree. BEGINS ----------------------------------------*/

// For storing the segment tree. (The root of the segment tree contains the length of the LIS and it's count).
vector<pii> tree;

/**
 * @brief Function for updating the segment tree. 
 * 
 * @param start left end of the Range/Interval which treeIdx represents.
 * @param end right end of the Range/Interval which treeIdx represents.
 * @param treeIdx Segment tree index in the array(Which represents the tree).
 * @param updateIdx Index to be updated.
 * @param length Length of LIS.
 * @param count Count of this LIS length.
 */
void update(int start, int end, int treeIdx, int updateIdx, int length, int count){

    /* If both intervals overlaps completely */
    if(start == end && start == updateIdx){
        /* Update the current tree node's length and count */
        tree[treeIdx].F = max(tree[treeIdx].F, length);
        tree[treeIdx].S = count;
        return;
    }

    /* If desired Index is completely outside the current Interval */
    if(updateIdx < start || end < updateIdx){
        /* Then ignore this subtree(do nothing) and return */
        return;
    }

    /* If intervals overlaps partially */
    int mid = (start + end) >> 1;

    /* Updating left subtree */
    update(start, mid, 2*treeIdx, updateIdx, length, count);
    /* Updating right subtree */
    update(mid+1, end, 2*treeIdx+1, updateIdx, length, count);

    /* If length of left & right child are equal */
    if(tree[2*treeIdx].F == tree[2*treeIdx+1].F){
        /* Then take the length same as left child, and update the count
        as the sum of left and right child */
        tree[treeIdx].F = tree[2*treeIdx].F;
        tree[treeIdx].S = tree[2*treeIdx].S + tree[2*treeIdx+1].S;
    }
    /* If length of left child is greater right child*/
    else if(tree[2*treeIdx].F > tree[2*treeIdx+1].F){
        /* Then take the length of Left child as it is greater. Count is also
        taken same as of the left child */
        tree[treeIdx] = tree[2*treeIdx];
    }
    /* If length of left child is less than right child */
    else if(tree[2*treeIdx].F > tree[2*treeIdx+1].F){
        /* Then take the length of Left child as it is greater. Count is also
        taken same as of the left child */
        tree[treeIdx] = tree[2*treeIdx+1];
    }
}

/**
 * @brief Function for finding the LIS and it's count in a given range.
 * 
 * @param start left end of the Range/Interval which treeIdx node represents.
 * @param end right end of the Range/Interval which treeIdx node represents.
 * @param treeIdx Index of current node in segment tree which represents the Interval [start, end].
 * @param queryStart left end of the queried Interval.
 * @param queryEnd right end of the queried Interval.
 * @return pair<int, int> 
 */
pii query(int start, int end, int treeIdx, int queryStart, int queryEnd){
    /* If both queried intervals is completely inside the CURRENT INTERVAL */
    if(queryStart <= start && end <= queryEnd){
        // Return length and count of current LIS 
        return tree[treeIdx];
    }

    /* If the intervals are completely DISJOINT */
    if(end < queryStart || queryEnd < start){
        return {INT_MIN, 0};
    }

    /* If both overlap partially */
    int mid = (start + end) >> 2;

    pii left = query(start, mid, 2*treeIdx, queryStart, queryEnd);
    pii right = query(mid+1, end, 2*treeIdx+1, queryStart, queryEnd);

    /* If length of left child is greater right child*/
    if(left.F > right.F)
        return left;

    /* If length of left child is less than right child */
    if(right.F > left.F)
        return right;
    
    /* If lenght of both children is equal then return count as 
    sum of left count and right count. */
    return {left.F, left.S + right.S};
}

/*-------------------- Auxillary function & objects for segment tree. ENDS ----------------------------------------*/


/* Comparator function for sorting array of pairs 
ACSENDING order of the 1st element and thereafter
in DESCENDING order of the 2nd element. */
bool customComparator(pii &a, pii &b){
    if(a.F == b.F)
        return a.S > b.S;
    return a.F < b.F;
}

/**
 * @brief Main function to solve the question
 * (finds the count of LIS in the given input array).
 * 
 * @param arr Input array.
 * @param n Length of Input array.
 * @return int 
 */
int solve(int arr[], int n){
    

    vector<pii> v(n);
    // v[i].F stores arr[i] element
    // v[i].S stores index of arr[i] in the original array.
    for(int i=0; i<n; i++)v[i].F = arr[i], v[i].S = i;

    /* Sort array of pairs in ASCENDING ORDER of elements value. */
    sort(v.begin(), v.end(), customComparator);

    for(int i=0; i<n; i++){
        
        int updateIdx = v[i].S;     

        if(updateIdx == 0){
            update(0, n-1, 1, updateIdx, 1, 1);
            continue;
        }

        // Find the LIS and it's count over the INTERVAL [0, updateIdx-1].
        pii tmp = query(0, n-1, 1, 0, updateIdx-1);

        // Updating the segment tree.
        update(0, n-1, 1, updateIdx, tmp.F+1, max(1, tmp.S));
    }

    pii finalAns = query(0, n-1, 1, 0, n-1);

    // return count of LIS of the array.
    return finalAns.S;
}


// Driver Function (It takes input and prints the final answer.)
int main(){
    int tt; cin >> tt;
    while(tt--){
        int n; cin >> n;
        int arr[n];

        tree.resize(4 * MAXSIZE + 1, {0, 0});

        for(int i=0; i<n; i++)cin >> arr[i];

        cout << solve(arr, n) << endl;
    }
    return 0;
}
