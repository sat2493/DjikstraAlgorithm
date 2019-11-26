#include <algorithm>
#include <iostream>
#include <utility>
#include <vector>
#include <stack>
#include <string>
#include <fstream>
#include <sstream>
#include "index_min_pq.h"

class Edge {
 public:
  Edge();
  Edge(unsigned int d, double w): dest(d), weight(w) {}
  ~Edge();
  unsigned int get_dest() { return dest; }
  double get_weight() { return weight; }
 private:
  unsigned int dest;
  double weight;
};

Edge::Edge() {}
Edge::~Edge() {}

class Vertex {
 public:
  Vertex();
  ~Vertex() {
    // delete Edge pointers in Edges vector
    for (std::vector<Edge*>::iterator it = Edges.begin();
         it != Edges.end(); it ++) {
      delete (*it);
    }
    Edges.clear();
  }
  void addEdge(Edge* e) { Edges.push_back(e); }
  Edge* getEdge() { return Edges[0]; }
  std::vector<Edge*> getEdges() { return Edges; }
 private:
  std::vector<Edge*> Edges;
};

Vertex::Vertex() {}

class Graph {
 public:
  Graph();
  ~Graph() {
    // delete Vertex pointers in Vertices vector
    for (std::vector<std::pair<unsigned int, Vertex*>*>
         ::iterator it = Vertices.begin();
         it != Vertices.end(); it++) {
      delete (*it);
    }
    Vertices.clear();
  }
  void addVertex(std::pair<unsigned int, Vertex*>* coordinates) {
    for (std::vector<std::pair<unsigned int, Vertex*>*>
         ::iterator it = Vertices.begin();
         it != Vertices.end(); it++) {
      // if Vertex already has been inserted
      if ((*it)->first == coordinates->first) {
        // insert another Edge
        (*it)->second->addEdge(coordinates->second->getEdge());
        // delete unnecessary memory held by ptr, coordinates
        delete coordinates;
        return;
      }
    }
    Vertices.push_back(coordinates);
  }

  std::vector<Edge*> Adj(unsigned int src) {
    std::vector<Edge*> adj;
    for (std::vector<std::pair<unsigned int, Vertex*>*>
         ::iterator it = Vertices.begin();
         it != Vertices.end(); it++) {
      if ((*it)->first == src) {
        adj = (*it)->second->getEdges();
        break;
      }
    }
    return adj;
  }

 private:
  std::vector<std::pair<unsigned int, Vertex*>*> Vertices;
};

void Dijkstra(Graph g, unsigned int src, unsigned int dest, unsigned int nv) {
  IndexMinPQ<double> q(nv);

  std::vector<double> dist;
  std::vector<int> prev;
  unsigned int i = 0;
  // initialize dist values to infinity
  // SOURCE for setting to infinity:
  // https://stackoverflow.com/questions/8690567/setting-an-int-to-infinity-
  // in-c
  // initialize prev values to undefined, -1
  while (i < nv) {
    dist.push_back(std::numeric_limits<double>::infinity());
    prev.push_back(-1);
    i++;
  }

  // Check if the number is in range of the vertices for user input
  if (src >= nv) {
    std::cerr << "Error: invalid source vertex number " << src << std::endl;
    exit(1);
  }

  if (dest >= nv) {
    std::cerr << "Error: invalid dest vertex number " << dest << std::endl;
    exit(1);
  }

  dist[src] = 0;

  q.Push(dist[src], src);

  unsigned int u;
  while (q.Size() != 0) {
    u = q.Top();
    if (u == dest)  break;
    std::vector<Edge*> adj_neighbors = g.Adj(u);
    for (std::vector<Edge*>::iterator it = adj_neighbors.begin();
         it != adj_neighbors.end(); it++) {
      unsigned int v = (*it)->get_dest();
      double alt = dist[u] + (*it)->get_weight();
      if (alt < dist[v]) {
        dist[v] = alt;
        prev[v] = u;
        if (q.Contains(v)) {
          q.ChangeKey(dist[v], v);
        } else {
          q.Push(dist[v], v);
        }
      }
    }
    q.Pop();
  }

  std::stack<unsigned int> nodesTraveled;
  u = dest;
  if (dist[u] != std::numeric_limits<double>::infinity()) {
    while (u != -1) {
      nodesTraveled.push(u);
      u = prev[u];
    }
  }

  // PRINTING AREA
  std::cout << src << " to " << dest << ": ";
  while (nodesTraveled.size()) {
    std::cout << nodesTraveled.top() << " ";
    if (nodesTraveled.size() != 1) {
      std::cout << "=> ";
    }
      nodesTraveled.pop();
  }
  if (dist[dest] == std::numeric_limits<double>::infinity()) {
    std::cout << "no path" << std::endl;
  } else {
    std::cout << "(" << dist[dest] << ")" << std::endl;
  }
}

Graph::Graph() {}

bool IntCheck(std::string input_str2) {
  std::istringstream cs2(input_str2);
  int complex_int;
  if (cs2 >> complex_int) {
    return true;
  } else {
    return false;
  }
}

int main(int argc, char* argv[]) {
  if (argc == 1) {
    std::cerr << "Usage: " << argv[0] << " <graph.dat> src dst" << std::endl;
    exit(1);
  }

  // invalid argv[2] or arv[3]

  std::ifstream graphFile;
  graphFile.open(argv[1]);
  if (!graphFile.good()) {
    std::cerr << "Error: cannot open file " << argv[1] << std::endl;
    exit(1);
  }

  std::string line;
  getline(graphFile, line);
  unsigned int nv = 0;

  int bool_checker = IntCheck(line);
  if (bool_checker) {
    nv = std::stoi(line);
  } else {
    std::cerr << "Error: invalid graph size" << std::endl;
    exit(1);
  }

  // build Graph
  Graph g;

  while (getline(graphFile, line)) {
    std::istringstream is(line);
    unsigned int src;
    unsigned int dest;
    double weight;
    while (is >> src >> dest >> weight) {
      // Check if the number is in the range of number of vertices from
      // file input
      if (src >= nv) {
        std::cerr << "Invalid source vertex number " << src << std::endl;
        exit(1);
      }

      if (dest >= nv) {
        std::cerr << "Invalid dest vertex number " << dest << std::endl;
        exit(1);
      }
      if (weight < 0) {
        std::cerr << "Invalid weight " << weight << std::endl;
        exit(1);
      }
      Vertex* v = new Vertex;
      Edge* e = new Edge(dest, weight);
      v->addEdge(e);

      std::pair<unsigned int, Vertex*>* coordinate =
      new std::pair<unsigned int, Vertex*>(src, v);
      g.addVertex(coordinate);
    }
  }
  unsigned int src = std::stoi(argv[2]);
  unsigned int dest = std::stoi(argv[3]);

  Dijkstra(g, src, dest, nv);

  return 0;
}
