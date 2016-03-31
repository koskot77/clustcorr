#ifndef unionfind_h
#define unionfind_h

#include <list>
#include <map>

class UnionFind {
private:
    std::map<int,int>              node2cluster;
    std::map<int, std::list<int> > cluster2nodes;

public:
    int joinClusters(int cluster1, int cluster2){
        if( cluster1 == cluster2 ) return cluster1;
        std::list<int> &nodes1 = cluster2nodes[cluster1];
        std::list<int> &nodes2 = cluster2nodes[cluster2];
        if( nodes1.size()==0 || nodes2.size()==0 ) return 0;
        int newCluster = 0;
        if( nodes1.size() < nodes2.size() ){
            newCluster = cluster2;
            for(std::list<int>::const_iterator n = nodes1.begin(); n != nodes1.end(); n++)
                node2cluster[*n] = newCluster;
            nodes2.insert(nodes2.end(),nodes1.begin(),nodes1.end());
            cluster2nodes.erase(cluster1);
        } else {
            newCluster = cluster1;
            for(std::list<int>::const_iterator n = nodes2.begin(); n != nodes2.end(); n++)
                node2cluster[*n] = newCluster;
            nodes1.insert(nodes1.end(),nodes2.begin(),nodes2.end());
            cluster2nodes.erase(cluster2);
        }
        return 0;
    }

    int findCluster(int node) { return node2cluster[node]; }

    int nClusters(void) const { return cluster2nodes.size(); }

    const std::map<int, std::list<int> >& clusters(void) const { return cluster2nodes; }

    UnionFind(int maxNodes){
        for(int i=1; i<=maxNodes; i++){
             node2cluster [i] = i;
             cluster2nodes[i].push_back(i);
        }
    }
};

#endif
