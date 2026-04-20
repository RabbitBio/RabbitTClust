#ifndef UNIONFIND_H
#define UNIONFIND_H


class UnionFind {
    int *parent, *ranks, _size;
    int _capacity;  // number of elements allocated in parent / ranks
public:
    UnionFind(){
        parent = nullptr; ranks = nullptr; _size = 0; _capacity = 0;
    }
    UnionFind(int size){
        parent = new int[size]; ranks = new int[size];
        for(int element = 0 ; element < size ; element++){
            parent[element] = element , ranks[element] = 0 ;
        }
        _size = size;
        _capacity = size;
    }
    void resize(int size){
        parent = new int[size]; ranks = new int[size];
        for(int element = 0 ; element < size ; element++){
            parent[element] = element , ranks[element] = 0 ;
        }
        _size = size;
        _capacity = size;
    }
    // Grow the structure by k new singleton elements. Useful for incremental
    // append where we do not know up front how many new clusters will be created.
    // Each new element starts in its own set (so _size grows by k too).
    void extend(int k){
        if(k <= 0) return;
        int old_cap = _capacity;
        int new_cap = _capacity + k;
        int* new_parent = new int[new_cap];
        int* new_ranks = new int[new_cap];
        for(int i = 0; i < old_cap; i++){
            new_parent[i] = parent[i];
            new_ranks[i] = ranks[i];
        }
        for(int i = old_cap; i < new_cap; i++){
            new_parent[i] = i;
            new_ranks[i] = 0;
        }
        delete[] parent; delete[] ranks;
        parent = new_parent;
        ranks = new_ranks;
        _capacity = new_cap;
        _size += k;
    }
    int capacity() const { return _capacity; }
    int find(int element){
        if(parent[element] == element){
            return element;
        }
        else{
            return parent[element] = find(parent[element]);          // Path Compression algorithm
        }
    }
    bool connected(int x,int y){
        if(find(x) == find(y)){
            return true;
        }
        else{
            return false;
        }
    }
    void merge(int x,int y){
        x = find(x);
        y = find(y);
        if(x != y){                                                   // Union by Rank algorithm
            if(ranks[x] > ranks[y]){
                parent[y] = x;
            }
            else if(ranks[x] < ranks[y]){
                parent[x] = y;
            }
            else{
                parent[x] = y; ranks[y] ++ ;
            }
            _size--;
        }
    }
    void clear(){
        delete [] parent; delete [] ranks;
    }
    int size(){
        return _size;
    }
};




#endif
