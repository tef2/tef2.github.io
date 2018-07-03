/*
 This file is part of TEF.

    TEF is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    TEF is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with TEF.  If not, see <http://www.gnu.org/licenses/>.


(C) 2018 Stratmann D, Pathmanathan JS, Postic G, Rey J, Chomilier J
Contact: dirk.stratmann@sorbonne-universite.fr
*/


#include <map>
#include <iostream>
#include <fstream>
#include <cmath>
#include "fastVector.h"
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <deque>

//#define LOG_FILE

using namespace std;

// overlapping
unsigned MIN = 2;

// maximum gap between two TEF, starting from the first possible TEF
unsigned MAX = 100;

// maximum CA-CA distance of TEF-ends
double dist_ca = 10.0;
double avg_dist_ca;

double tef_weight = 0.0;

double gap_weight = 1.0;

double dist_weight = 1.0;

ofstream logfile;

///////////////////////////////////////////////////////////////////////////////////
// from Robert Sedgewick's book:
// Algorithms in C++, Part 5, Graph Algorithms, Thrid Edition


class WeightedEdge
{
    public:
        WeightedEdge(){
            init(-999, -999, -999);
        }
        WeightedEdge(int a_, int b_, double weight_) {
            init(a_, b_, weight_);
        }
        void init (int a_, int b_, double weight_) {
            a = a_; b = b_; weight = weight_;
        }

        int v() const { return a; }
        int w() const { return b; }
        double wt() const { return weight; }
        ostream& show(ostream& output){
            output << v() << "-" << w() << " " << wt() << endl;
            return output;
        }
        /*bool from() { return true; }
        int other(int a_) {
            if (from(a_)){
                return w();
            } else{
                return v();
            }
        }*/
    protected:
        int a, b;
        double weight;
};

class TEFedge : public WeightedEdge
{
    public:
        TEFedge(){
            init(-999,-999,-999,-999,-999,-999);
        }
        TEFedge(int a_, int b_){
            init(a_, b_, -999,-999,-999,-999);
        }
        TEFedge(int a_, int b_, double weight_, int gap_, double dist_, int numTEF_) {
            init(a_, b_, weight_, gap_, dist_, numTEF_);
        }
        void init (int a_, int b_, double weight_, int gap_, double dist_, int numTEF_) {
            a = a_; b = b_; weight = weight_;
            gap = gap_; dist = dist_; numTEF = numTEF_;
        }
        void setWeight (double weight_, int gap_, double dist_, int numTEF_){
            weight = weight_;
            gap = gap_; dist = dist_; numTEF = numTEF_;
        }

        ostream& show(ostream& output){
            output << v() << "-" << w() << " " << wt() << " gap: " << gap << " dist: " << dist << " numTEF: " << numTEF << endl;
            return output;
        }

    private:
        int gap;
        double dist;
        int numTEF;
};

template <class Graph, class Edge> void showGraph(const Graph &G, ostream& output)
{
    for (int s = 0; s < G.V(); s++) {
        output.width(2);
        output << s << ":" << endl;
        typename Graph::adjIterator A(G, s);
        for (Edge* e = A.beg(); !A.end(); e = A.nxt()){
            output.width(2);
            e->show(output);
        }
        output << endl;
    }
}

template <class Edge> class DenseGRAPH
{
    int Vcnt, Ecnt;
    bool digraph;
    vector <vector <Edge *> > adj;
public:
    DenseGRAPH(int V, bool digraph = false) :
        Vcnt(V), Ecnt(0), digraph(digraph), adj(V) {
        for (int i = 0; i < V; i++)
            adj[i].assign(V, (Edge*)0);
    }
    int V() const {
        return Vcnt;
    }
    int E() const {
        return Ecnt;
    }
    bool directed() const {
        return digraph;
    }
    void insert(Edge *e) {
        int v = e->v(), w = e->w();
        if (adj[v][w] == NULL) Ecnt++;
        adj[v][w] = e;
        if (!digraph) adj[w][v] = e;
    }
    void remove(Edge *e) {
        int v = e->v(), w = e->w();
        if (adj[v][w] != NULL) Ecnt--;
        adj[v][w] = NULL;
        if (!digraph) adj[w][v] = NULL;
    }
    Edge* edge(int v, int w) const {
        return adj[v][w];
    }
    class adjIterator;
    friend class adjIterator;
};

template <class Edge>
class DenseGRAPH<Edge>::adjIterator
{
    int i, v;
    const DenseGRAPH<Edge> &G;
public:
    adjIterator(const DenseGRAPH<Edge> &G, int v) :
        i(0), v(v), G(G) { }
    Edge *beg() {
        i = -1;
        return nxt();
    }
    Edge *nxt() {
        for (i++; i < G.V(); i++)
            if (G.edge(v, i)) return G.adj[v][i];
        return 0;
    }
    bool end() const {
        return i >= G.V();
    }
};

template <class Dag, class Edge> class dagTS
{
    bool reverse;
    vector<int> pre, post, postI;
    const Dag &D;
    int tcnt, cnt;
    // gives a reverse topological sort
    void tsR_reverse(int v) {
        pre[v] = cnt++;
        typename Dag::adjIterator A(D, v);
        //for (int t = A.beg(); !A.end(); t = A.nxt())
        for (Edge* e = A.beg(); !A.end(); e = A.nxt())
            if (pre[e->w()] == -1) tsR_reverse(e->w());
        post[v] = tcnt;
        postI[tcnt++] = v;
    }
    // gives a topological sort (not reverse)
    void tsR(int v)
    {
      pre[v] = cnt++;
      for (int w = 0; w < D.V(); w++)
        if (D.edge(w, v))
          if (pre[w] == -1) tsR(w);
      post[v] = tcnt; postI[tcnt++] = v;
    }

public:
    dagTS(const Dag &D_, bool reverse_) : reverse(reverse_), pre(D_.V(), -1), post(D_.V(), -1), postI(D_.V(), -1), D(D_), tcnt(0), cnt(0) {
        for (int v = 0; v < D.V(); v++)
            if (pre[v] == -1){
                if (reverse) tsR_reverse(v);
                else tsR(v);
            }
    }
    int operator[](int v) const {
        return postI[v];
    }
    int relabel(int v) const {
        return post[v];
    }
};


template <class Item>
class QUEUE
{
private:
    struct node {
        Item item;
        node* next;
        node(Item x) {
            item = x;
            next = 0;
        }
    };
    typedef node *link;
    link head, tail;
public:
    QUEUE() {
        head = 0;
    }
    int empty() const {
        return head == 0;
    }
    void put(Item x) {
        link t = tail;
        tail = new node(x);
        if (head == 0)
            head = tail;
        else t->next = tail;
    }
    Item get() {
        Item v = head->item;
        link t = head->next;
        delete head;
        head = t;
        return v;
    }
};

template <class Dag> class dagTS_queue
{
    const Dag &D;
    vector<int> in, ts, tsI;
public:
    dagTS_queue(const Dag &D) : D(D),
        in(D.V(), 0), ts(D.V(), -1), tsI(D.V(), -1) {
        QUEUE<int> Q;
        for (int v = 0; v < D.V(); v++) {
            typename Dag::adjIterator A(D, v);
            for (int t = A.beg(); !A.end(); t = A.nxt())
                in[t]++;
        }
        for (int v = 0; v < D.V(); v++)
            if (in[v] == 0) Q.put(v);
        for (int j = 0; !Q.empty(); j++) {
            ts[j] = Q.get();
            tsI[ts[j]] = j;
            typename Dag::adjIterator A(D, ts[j]);
            for (int t = A.beg(); !A.end(); t = A.nxt())
                if (--in[t] == 0) Q.put(t);
        }
    }
    int operator[](int v) const {
        return ts[v];
    }
    int relabel(int v) const {
        return tsI[v];
    }
};


template <class Graph, class Edge> class SPTdag_queue
{
    const Graph &G;
    vector<double> wt;
    vector<Edge *> spt;
    vector<int> in, ts, tsI;
    deque<int> solution;

    void generateSolution(){
        int v = G.V() - 1;
        solution.push_front(v);
        Edge* edge = pathR(v);
        while(edge){
            v = edge->v();
            edge = pathR(v);
            solution.push_front(v);
        }
    }

public:
    SPTdag_queue(const Graph &G) : G(G), wt(G.V(), 999999), spt(G.V(), (Edge*)0),
        in(G.V(), 0), ts(G.V(), -1), tsI(G.V(), -1) {
        int w;
        wt[0] = 0.0;
        QUEUE<int> Q;
        for (int v = 0; v < G.V(); v++) {
            typename Graph::adjIterator A(G, v);
            for (Edge* e = A.beg(); !A.end(); e = A.nxt())
                in[e->w()]++;
        }
        /*logfile << "Sources: " << endl;
        for (int v = 0; v < G.V(); v++)
            if (in[v] == 0){
                logfile << v << endl;
                Q.put(v);
            }*/
        Q.put(0); // only the node 0 (=start of the amino acid sequence) should be a source
        // other sources may appear with low MAX_GAP values, but will be ignored
        /*for (int v = 0; v < G.V(); v++)
            if (in[v] == 0){
                typename Graph::adjIterator A(G, v);
                if (v==0) Q.put(v);
                else if (A.beg()!=0) logfile << "ERROR found another source v= " << v << endl;
            }*/

        //logfile << "Queue-search:" << endl;
        for (int j = 0; !Q.empty(); j++) {
            ts[j] = Q.get();
            int v = ts[j];
//            logfile << endl;
//            logfile << "v=" << v << " wt[v]=" << wt[v] << endl;
            tsI[v] = j;
            typename Graph::adjIterator A(G, v);
            for (Edge* e = A.beg(); !A.end(); e = A.nxt()){
                w = e->w();
                // DEBUG
                //double wt_w = wt[w];
                //double wt_v = wt[v];
                //double ewt = e->wt();
                //logfile << "  w=" << w << " wt[w]=" << wt_w << " e->wt()=" << ewt << endl;
                // END DEBUG
                if (wt[w] > wt[v] + e->wt()) {
                    wt[w] = wt[v] + e->wt();
                    spt[w] = e;
                }

                if (--in[w] == 0){
                    Q.put(w);
                    //logfile << "  put w on queue" << endl;
                }
            }
        }
        generateSolution();
    }
    Edge *pathR(int v) const {
        return spt[v];
    }
    double dist(int v) const {
        return wt[v];
    }

    const deque<int>& pathFromSource() const {
        return solution;
    }
};


template <class Graph, class Edge> class SPTdag
{
    const Graph &G;
    vector<double> wt;
    vector<Edge *> spt;
    deque<int> solution;

    void generateSolution(){
        int v = G.V() - 1;
        solution.push_front(v);
        Edge* edge = pathR(v);
        while(edge){
            v = edge->v();
            edge = pathR(v);
            solution.push_front(v);
        }
    }

public:
    SPTdag(const Graph &G) : G(G), wt(G.V(), 999999), spt(G.V(), (Edge*)0) {
        int j, w;
        dagTS<Graph, Edge> ts(G, false);
        //wt[wt.size()-1] = 0.0;
        wt[0] = 0.0;
        for (int v = ts[j = 0]; j < G.V(); v = ts[++j]) {
            typename Graph::adjIterator A(G, v);
            for (Edge* e = A.beg(); !A.end(); e = A.nxt()){
                w = e->w();
                // DEBUG
                //double wt_w = wt[w];
                //double wt_v = wt[v];
                //double ewt = e->wt();
                // END DEBUG
                if (wt[w] > wt[v] + e->wt()) {
                    wt[w] = wt[v] + e->wt();
                    spt[w] = e;
                }
            }
        }

        generateSolution();
    }
    Edge *pathR(int v) const {
        return spt[v];
    }
    double dist(int v) const {
        return wt[v];
    }

    const deque<int>& pathFromSource() const {
        return solution;
    }

};


template <class Graph, class Edge> class LPTdag
{
    const Graph &G;
    vector<double> wt;
    vector<Edge *> lpt;
public:
    LPTdag(const Graph &G) : G(G), wt(G.V(), 0), lpt(G.V(), (Edge*)0) {
        int j, w;
        dagTS<Graph, Edge> ts(G, true);
        for (int v = ts[j = 0]; j < G.V(); v = ts[++j]) {
            typename Graph::adjIterator A(G, v);
            for (Edge* e = A.beg(); !A.end(); e = A.nxt()){
                w = e->w();
                // DEBUG
                double wt_w = wt[w];
                double wt_v = wt[v];
                double ewt = e->wt();
                // END DEBUG
                if (wt[w] < wt[v] + e->wt()) {
                    wt[w] = wt[v] + e->wt();
                    lpt[w] = e;
                }
            }
        }
    }
    Edge *pathR(int v) const {
        return lpt[v];
    }
    double dist(int v) const {
        return wt[v];
    }
};


// END of graph-algorithms from Robert Sedgewick's book:
// Algorithms in C++, Part 5, Graph Algorithms, Thrid Edition
////////////////////////////////////////////////////////////////




struct Solution_struct {
    fastDeque<unsigned> solution;
    double score;
};

class TEF_list
{
public:
    TEF_list(char *filname);

    int firstResID, lastResID;
    fastVector<int> tef_beg;
    fastVector<int> tef_end;
    fastVector<double> tef_cad; //CA-CA distance between the ends of a TEF
    map<int, fastDeque<unsigned> > tef_index;
    unsigned numTEF;
};

TEF_list::TEF_list(char *filename)
{
    int beg,end;
    double cad;
    fastDeque<int> t_beg, t_end;
    fastDeque<double> t_cad;

    ifstream input(filename);
    input >> firstResID;
    input >> lastResID;
    while (! input.eof() ) {
        input >> beg;
        input >> end;
        input >> cad;
        //cout << beg << " " << end << " " << cad << endl;
        if (input.good()) {
            t_beg.push_back(beg);
            t_end.push_back(end);
            t_cad.push_back(cad);
        }
    }
    tef_beg = t_beg;
    tef_end = t_end;
    tef_cad = t_cad;

    numTEF = tef_beg.size();

    for (unsigned i=0; i < numTEF; i++) {
        tef_index[tef_beg[i]].push_back(i);
    }
}


class Search
{
public:
    Search(TEF_list& tef_list) : tef_graph(tef_list.tef_beg.size()+2, true), tef_edges(100000) {
        firstResID = tef_list.firstResID;
        lastResID = tef_list.lastResID;
        Nres = lastResID - firstResID + 1;
        firstTEFres.resize(Nres);
        tef_beg = tef_list.tef_beg;
        tef_end = tef_list.tef_end;
        tef_cad = tef_list.tef_cad;
        tef_index = tef_list.tef_index;
    }
    void doGraphSearch();
    void testGraph();
    void writeSolution(const char *filename);
    fastDeque<unsigned> getBestSolution() const{
        return best_solution;
    }
    void clearBestSolution(){
        best_solution.clear();
    }

private:
    void buildGraph();
    double getSingleCost(int oldEnd, int newStart, double dist, unsigned numTEF);
    double getSingleCost(int oldEnd, int newStart, double dist, unsigned numTEF, TEFedge& edge);

    fastVector<int> tef_beg;
    fastVector<int> tef_end;
    fastVector<double> tef_cad; //CA-CA distance between the ends of a TEF
    map<int, fastDeque<unsigned> > tef_index;
    fastVector<int> firstTEFres;


    int lastResID;
    int firstResID;
    int Nres;
    fastDeque<unsigned> best_solution;
    double best_score;

#ifdef LOG_FILE
    DenseGRAPH<TEFedge> tef_graph;
    fastDeque<TEFedge> tef_edges;
#else
    DenseGRAPH<WeightedEdge> tef_graph;
    fastDeque<WeightedEdge> tef_edges;
#endif

};


void Search::testGraph()
{
    DenseGRAPH<WeightedEdge> graph(10, true);
    fastDeque<WeightedEdge> edges;

    edges.push_back(WeightedEdge(0,1,0.41));
    edges.push_back(WeightedEdge(0,7,0.41));
    edges.push_back(WeightedEdge(0,9,0.41));
    edges.push_back(WeightedEdge(1,2,0.51));
    edges.push_back(WeightedEdge(7,8,0.32));
    edges.push_back(WeightedEdge(7,3,0.32));
    edges.push_back(WeightedEdge(6,8,0.21));
    edges.push_back(WeightedEdge(6,3,0.21));
    edges.push_back(WeightedEdge(8,2,0.32));
    edges.push_back(WeightedEdge(9,4,0.29));
    edges.push_back(WeightedEdge(9,6,0.29));
    edges.push_back(WeightedEdge(2,10,0.50));
    edges.push_back(WeightedEdge(3,10,0.36));
    edges.push_back(WeightedEdge(4,10,0.38));
    edges.push_back(WeightedEdge(5,10,0.45));
}


void Search::buildGraph()
{
    vector<int> indegree(tef_graph.V(), 0);

    unsigned numTEF = tef_beg.size();
    for (int i=0; i<(int)firstTEFres.size(); i++){
        firstTEFres[i] = 999999;
        for (map<int, fastDeque<unsigned> >::iterator it = tef_index.begin(); it != tef_index.end(); it++) {
            if (it->first > firstResID + i - (int)MIN) {
                if (it->first < firstTEFres[i]) firstTEFres[i] = it->first;
            }
        }
    }
    int minBegin = lastResID;
    for (int i=0; i<(int)tef_beg.size(); i++){
        if (tef_beg[i]<minBegin) minBegin = tef_beg[i];
    }
    int maxEnd = firstResID;
    for (int i=0; i<(int)tef_end.size(); i++){
        if (tef_end[i]>maxEnd) maxEnd = tef_end[i];
    }
    for (int i=0; i<(int)tef_beg.size(); i++) {
        if (tef_beg[i] < minBegin + (int)MAX){
#ifdef LOG_FILE
            TEFedge edge(0, i+1);
            getSingleCost(firstResID-1, tef_beg[i], tef_cad[i], 1, edge);
            tef_edges.push_back(edge);
#else
            double weight = getSingleCost(firstResID-1, tef_beg[i], tef_cad[i], 1);
            tef_edges.push_back(WeightedEdge(0, i+1, weight));
#endif
            indegree[i+1]++;
        }
        int firstTEF = firstTEFres[tef_end[i]-firstResID];
        if (firstTEF == 999999){ // no more TEF after this one
#ifdef LOG_FILE
            TEFedge edge(i+1, numTEF+1);
            getSingleCost(tef_end[i], lastResID+1, avg_dist_ca, 0, edge);
            tef_edges.push_back(edge);
#else
            double weight = getSingleCost(tef_end[i], lastResID+1, avg_dist_ca, 0);
            tef_edges.push_back(WeightedEdge(i+1, numTEF+1, weight));
#endif
            indegree[numTEF+1]++;
        }
        int tefEnd = tef_end[i];
        int startRes = -999;
        for (map<int, fastDeque<unsigned> >::iterator it = tef_index.begin(); it != tef_index.end(); it++) {
            if (startRes == -999 && it->first > tefEnd - (int)MIN) {
                startRes = it->first;
            }
            if (startRes > -999) {
                if (it->first <= startRes + (int)MAX) { // to change into (tef_end + MAX) for startRes < tef_end ?? NO, as there can be large gaps without any TEF.
                    unsigned numElem = it->second.size();
                    for (unsigned j=0; j < numElem; j++) {
                        unsigned tefID = it->second[j];
                        if (i != (int)tefID){
#ifdef LOG_FILE
                            TEFedge edge(i+1, tefID+1);
                            getSingleCost(tefEnd, tef_beg[tefID], tef_cad[tefID], 1, edge);
                            tef_edges.push_back(edge);
#else
                            double weight = getSingleCost(tefEnd, tef_beg[tefID], tef_cad[tefID], 1);
                            tef_edges.push_back(WeightedEdge(i+1, tefID+1, weight));
#endif
                            indegree[tefID+1]++;
                        }
                    }
                } else {
                    break;
                }
            }
        }
    }

    //logfile << "buildGraph():" << endl;
    for (int i=0; i<(int)tef_edges.size(); i++){
        WeightedEdge& edge = tef_edges[i];
        if (edge.v() == 0 || indegree[edge.v()] > 0){
            tef_graph.insert(tef_edges.getPointer(i));
        } else{
            indegree[edge.w()]--;
            //logfile << "skip edge from source v=" << edge.v() << endl;
        }
    }

#ifdef LOG_FILE
    //showGraph<DenseGRAPH<TEFedge>, TEFedge>(tef_graph, logfile);
    //showGraph<DenseGRAPH<WeightedEdge>, WeightedEdge>(tef_graph, logfile);
#endif

}


double Search::getSingleCost(int oldEnd, int newStart, double dist, unsigned numTEF)
{
    int gap = abs(newStart - oldEnd - 1);
    dist -= avg_dist_ca;
    double cost = gap_weight * gap + dist_weight * dist - 10 * tef_weight * numTEF;
    return cost;
}

double Search::getSingleCost(int oldEnd, int newStart, double dist, unsigned numTEF, TEFedge& edge)
{
    int gap = abs(newStart - oldEnd - 1);
    dist -= avg_dist_ca;
    double cost = gap_weight * gap + dist_weight * dist - 10 * tef_weight * numTEF;
    edge.setWeight(cost, gap, dist, numTEF);
    return cost;
}



void Search::doGraphSearch()
{
#ifdef LOG_FILE
    clock_t t1, t2;
    logfile << "buildGraph()..." << endl;
    t1 = clock();
#endif

    buildGraph();

#ifdef LOG_FILE
    t2 = clock();
    logfile << "clocks elapsed = " << t2 - t1 << endl;
    logfile << "search shortest Path, first method..." << endl;

    t1 = clock();
    for (int i=0; i<1000; i++){
        SPTdag<DenseGRAPH<TEFedge>, TEFedge> shortestPath(tef_graph);
    }
    t2 = clock();
    logfile << "clocks elapsed = " << t2 - t1 << endl;
#endif

#ifdef LOG_FILE
    SPTdag<DenseGRAPH<TEFedge>, TEFedge> shortestPath(tef_graph);
#else
    SPTdag<DenseGRAPH<WeightedEdge>, WeightedEdge> shortestPath(tef_graph);
#endif

#ifdef LOG_FILE
    logfile << "search shortest Path, second method..." << endl;
    t1 = clock();
    for (int i=0; i<1000; i++){
        SPTdag_queue<DenseGRAPH<TEFedge>, TEFedge> shortestPath_queue(tef_graph);
    }
    t2 = clock();
    logfile << "clocks elapsed = " << t2 - t1 << endl; //Both methods take about the same time (the second one a bit longer)
#endif

#ifdef LOG_FILE
    SPTdag_queue<DenseGRAPH<TEFedge>, TEFedge> shortestPath_queue(tef_graph);
#else
    SPTdag_queue<DenseGRAPH<WeightedEdge>, WeightedEdge> shortestPath_queue(tef_graph);
#endif

#ifdef LOG_FILE
    logfile << endl << "doGraphSearch() result:" << endl;
    logfile << "first method:" << endl;
#endif
/*    int v = tef_graph.V()-1;
    TEFedge* edge = shortestPath.pathR(v);
    while(edge){
        double dist = shortestPath.dist(v);
        logfile << "w=" << v << " ";
        v = edge->v();
        logfile << "v=" << v << " dist=" << dist << endl;
        edge = shortestPath.pathR(v);
    }*/
#ifdef LOG_FILE
    const deque<int>& solution = shortestPath.pathFromSource();
    for (unsigned i=0; i<solution.size(); i++){
        logfile << solution[i] << " dist=" << shortestPath.dist(solution[i]) << endl;
    }

    logfile << "second method:" << endl;
#endif
/*    v = tef_graph.V()-1;
    edge = shortestPath_queue.pathR(v);
    while(edge){
        double dist = shortestPath_queue.dist(v);
        logfile << "w=" << v << " ";
        v = edge->v();
        logfile << "v=" << v << " dist=" << dist << endl;
        edge = shortestPath_queue.pathR(v);
    }*/
    const deque<int>& solutionTwo = shortestPath_queue.pathFromSource();
#ifdef LOG_FILE
    for (unsigned i=0; i<solutionTwo.size(); i++){
        logfile << solutionTwo[i] << " dist=" << shortestPath_queue.dist(solutionTwo[i]) << endl;
    }
#endif

    if (fabs(shortestPath.dist(tef_graph.V()-1) - shortestPath_queue.dist(tef_graph.V()-1)) > 1e-7){
#ifdef LOG_FILE
        logfile << "ERROR: the two graph-based methods yields different solutions" << endl;
#endif
        cout << "ERROR: the two graph-based methods yields different solutions" << endl;
    } else if (solutionTwo.size() > 0){
        best_score = shortestPath_queue.dist(tef_graph.V()-1);
        best_solution.clear();
        for (unsigned i=1; i<solutionTwo.size()-1; i++){
            best_solution.push_back(solutionTwo[i]-1);
        }
    } else{
#ifdef LOG_FILE
        logfile << "ERROR: no solution found!" << endl;
#endif
        cout << "ERROR: no solution found!" << endl;
    }
}


void Search::writeSolution(const char *filename)
{
    ofstream output(filename);
    output << best_score << endl;
    for (unsigned i=0; i < best_solution.size(); i++) {
        output << best_solution[i] << endl;
    }
}


int main (int argc, char *argv[])
{
    if (argc != 9) {
        cout << "usage: tef_search inputFile outputFile MIN MAX dist_ca tef_weight gap_weight dist_weight" << endl;
        return -1;
    }
#ifdef LOG_FILE
    logfile.open("tef_search.log");
#endif
    // chevauchement
    MIN = atoi(argv[3]);
    //printf("MIN = %u \n", MIN);
    // gap lors de la selection du prochain tef
    MAX = atoi(argv[4]);
    //printf("MAX = %u \n", MAX);
    // distance max entre deux carbones alpha
    dist_ca = atof(argv[5]);
    avg_dist_ca = dist_ca / 2.0 + (dist_ca - dist_ca / 2.0) / 2.0;
    //cout << "avg_dist_ca = " << avg_dist_ca << endl;
    // poids d'un tef pour le calcul du score
    tef_weight = atof(argv[6]);
    //printf("tef weight = %f \n", tef_weight);
    gap_weight = atof(argv[7]);
    //printf("gap weight = %f \n", gap_weight);
    dist_weight = atof(argv[8]);
    //printf("distance weight = %f \n", dist_weight);
    TEF_list tef_list(argv[1]);
    Search s(tef_list);
    s.doGraphSearch();
    s.writeSolution(argv[2]);
    int returnCode = 0;
#ifdef LOG_FILE
    logfile.close();
#endif
    return returnCode;
}


