//
// Algoritmos e Estruturas de Dados --- 2023/2024
//
// Joaquim Madeira, Joao Manuel Rodrigues - June 2021, Nov 2023
//
// Graph - Using a list of adjacency lists representation
//

#include "Graph.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "SortedList.h"
#include "instrumentation.h"

#define ADJVERTEX InstrCount[1]
struct _Vertex
{
  unsigned int id;
  unsigned int inDegree;
  unsigned int outDegree;
  List *edgesList;
};

struct _Edge
{
  unsigned int adjVertex;
  double weight;
};

struct _GraphHeader
{
  int isDigraph;
  int isComplete;
  int isWeighted;
  unsigned int numVertices;
  unsigned int numEdges;
  List *verticesList;
};

// The comparator for the VERTICES LIST

int graphVerticesComparator(const void *p1, const void *p2)
{
  unsigned int v1 = ((struct _Vertex *)p1)->id;
  unsigned int v2 = ((struct _Vertex *)p2)->id;
  int d = v1 - v2;
  return (d > 0) - (d < 0);
}

// The comparator for the EDGES LISTS

int graphEdgesComparator(const void *p1, const void *p2)
{
  unsigned int v1 = ((struct _Edge *)p1)->adjVertex;
  unsigned int v2 = ((struct _Edge *)p2)->adjVertex;
  int d = v1 - v2;
  return (d > 0) - (d < 0);
}

Graph *GraphCreate(unsigned int numVertices, int isDigraph, int isWeighted)
{
  Graph *g = (Graph *)malloc(sizeof(struct _GraphHeader));
  if (g == NULL)
    abort();

  g->isDigraph = isDigraph;
  g->isComplete = 0;
  g->isWeighted = isWeighted;

  g->numVertices = numVertices;
  g->numEdges = 0;

  g->verticesList = ListCreate(graphVerticesComparator);

  for (unsigned int i = 0; i < numVertices; i++)
  {
    struct _Vertex *v = (struct _Vertex *)malloc(sizeof(struct _Vertex));
    if (v == NULL)
      abort();

    v->id = i;
    v->inDegree = 0;
    v->outDegree = 0;

    v->edgesList = ListCreate(graphEdgesComparator);

    ListInsert(g->verticesList, v);
  }

  assert(g->numVertices == ListGetSize(g->verticesList));

  return g;
}

Graph *GraphCreateComplete(unsigned int numVertices, int isDigraph)
{
  Graph *g = GraphCreate(numVertices, isDigraph, 0);

  g->isComplete = 1;

  List *vertices = g->verticesList;
  ListMoveToHead(vertices);
  unsigned int i = 0;
  for (; i < g->numVertices; ListMoveToNext(vertices), i++)
  {
    struct _Vertex *v = ListGetCurrentItem(vertices);
    List *edges = v->edgesList;
    for (unsigned int j = 0; j < g->numVertices; j++)
    {
      if (i == j)
      {
        continue;
      }
      struct _Edge *new = (struct _Edge *)malloc(sizeof(struct _Edge));
      if (new == NULL)
        abort();
      new->adjVertex = j;
      new->weight = 1;

      ListInsert(edges, new);
    }
    if (g->isDigraph)
    {
      v->inDegree = g->numVertices - 1;
      v->outDegree = g->numVertices - 1;
    }
    else
    {
      v->outDegree = g->numVertices - 1;
    }
  }
  if (g->isDigraph)
  {
    g->numEdges = numVertices * (numVertices - 1);
  }
  else
  {
    g->numEdges = numVertices * (numVertices - 1) / 2;
  }

  return g;
}

void GraphDestroy(Graph **p)
{
  assert(*p != NULL);
  Graph *g = *p;

  List *vertices = g->verticesList;
  if (ListIsEmpty(vertices) == 0)
  {
    ListMoveToHead(vertices);
    unsigned int i = 0;
    for (; i < g->numVertices; ListMoveToNext(vertices), i++)
    {
      struct _Vertex *v = ListGetCurrentItem(vertices);

      List *edges = v->edgesList;
      if (ListIsEmpty(edges) == 0)
      {
        unsigned int i = 0;
        ListMoveToHead(edges);
        for (; i < ListGetSize(edges); ListMoveToNext(edges), i++)
        {
          struct _Edge *e = ListGetCurrentItem(edges);
          free(e);
        }
      }
      ListDestroy(&(v->edgesList));
      free(v);
    }
  }

  ListDestroy(&(g->verticesList));
  free(g);

  *p = NULL;
}

Graph *GraphCopy(const Graph *g)
{
  assert(g != NULL);

  // TO BE COMPLETED !!
  Graph *copy = GraphCreate(g->numVertices, g->isDigraph, g->isWeighted);
  if(copy == NULL){
    return NULL;
  }

  List *vertices = g->verticesList; // lista de vertices do grafo original
  ListMoveToHead(vertices);
  ListMoveToHead(copy->verticesList);

  for (unsigned int i = 0; i < g->numVertices; i++)
  {
    struct _Vertex *originalVertex = ListGetCurrentItem(vertices);                 // vertice atual da lista de vertices do grafo original
    struct _Vertex *copyVertex = ListGetCurrentItem(copy->verticesList); // criar vertice para o grafo copia
    
    copyVertex->id = originalVertex->id;
    copyVertex->inDegree = originalVertex->inDegree; // copiar informacao do vertice original para o vertice copia
    copyVertex->outDegree = originalVertex->outDegree;

    List *originalEdges = originalVertex->edgesList; // lista de arestas do vertice original
    ListMoveToHead(originalEdges);

    for (unsigned int j = 0; j < originalVertex->outDegree; j++)
    {
      struct _Edge *originalEdge = ListGetCurrentItem(originalEdges);        // aresta atual da lista de arestas do vertice original
      struct _Edge *copyEdge = (struct _Edge *)malloc(sizeof(struct _Edge)); // criar aresta para o vertice copia
      if (copyEdge == NULL)
      {
        free(copyEdge);
        GraphDestroy(&copy);
        abort();
      }
      copyEdge->adjVertex = originalEdge->adjVertex; // copiar informacao da aresta original para a aresta copia
      ADJVERTEX++;
      copyEdge->weight = originalEdge->weight;

      ListInsert(copyVertex->edgesList, copyEdge); // inserir aresta na lista de arestas do vertice copia
      copy->numEdges++;                              
      ListMoveToNext(originalEdges);               // mover para a proxima aresta da lista de arestas do vertice original
    }

    ListMoveToNext(vertices);                   // mover para o proximo vertice da lista de vertices do grafo original
    ListMoveToNext(copy->verticesList);  
  }

  return copy;
}

static int _addEdge(Graph *g, unsigned int v, unsigned int w, double weight);

Graph *GraphFromFile(FILE *f)
{
  assert(f != NULL);
  // TO BE COMPLETED !!

  unsigned int numVertices, numEdges, weighted, directed; // variaveis para guardar informacao do grafo

  if(fscanf(f, "%u %u %u %u", &directed, &weighted, &numVertices, &numEdges) != 4) // ler informacao do grafo
  {
    fprintf(stderr, "Error: invalid file format.\n");
    return NULL;
  }

  Graph *g = GraphCreate(numVertices, directed, weighted); // criar grafo

  for(unsigned int i = 0; i < numEdges; i++){ // ler arestas do grafo
    unsigned int vi, vf;
    double weight = 1.0;
    if(weighted){
      if(fscanf(f, "%u %u %lf", &vi, &vf, &weight) != 3){ // ler aresta com peso
        fprintf(stderr, "Error: invalid file format.\n");
        GraphDestroy(&g);
        return NULL;
      }
    }
    else{
      if(fscanf(f, "%u %u", &vi, &vf) != 2){ // ler aresta sem peso
        fprintf(stderr, "Error: invalid file format.\n");
        GraphDestroy(&g);
        return NULL;
      }
    }

    if(vi == vf){ // verificar se a aresta e um lacete
      continue;
    }
    
    int validEdge = 1;
    unsigned int *adjacent = GraphGetAdjacentsTo(g, vi);
    for(unsigned int vertex = 1; vertex <= adjacent[0]; vertex++){
      if(adjacent[vertex] == vf){
        validEdge = 0;
        break;
      }
    }
    free(adjacent);
    if(!validEdge){
      continue;
    }
    _addEdge(g,vi,vf,weight); // adicionar aresta ao grafo
  }

  return g;
}

// Graph

int GraphIsDigraph(const Graph *g) { return g->isDigraph; }

int GraphIsComplete(const Graph *g) { return g->isComplete; }

int GraphIsWeighted(const Graph *g) { return g->isWeighted; }

unsigned int GraphGetNumVertices(const Graph *g) { return g->numVertices; }

unsigned int GraphGetNumEdges(const Graph *g) { return g->numEdges; }

//
// For a graph
//
double GraphGetAverageDegree(const Graph *g)
{
  assert(g->isDigraph == 0);
  return 2.0 * (double)g->numEdges / (double)g->numVertices;
}

static unsigned int _GetMaxDegree(const Graph *g)
{
  List *vertices = g->verticesList;
  if (ListIsEmpty(vertices))
    return 0;

  unsigned int maxDegree = 0;
  ListMoveToHead(vertices);
  unsigned int i = 0;
  for (; i < g->numVertices; ListMoveToNext(vertices), i++)
  {
    struct _Vertex *v = ListGetCurrentItem(vertices);
    if (v->outDegree > maxDegree)
    {
      maxDegree = v->outDegree;
    }
  }
  return maxDegree;
}

//
// For a graph
//
unsigned int GraphGetMaxDegree(const Graph *g)
{
  assert(g->isDigraph == 0);
  return _GetMaxDegree(g);
}

//
// For a digraph
//
unsigned int GraphGetMaxOutDegree(const Graph *g)
{
  assert(g->isDigraph == 1);
  return _GetMaxDegree(g);
}

// Vertices

//
// returns an array of size (outDegree + 1)
// element 0, stores the number of adjacent vertices
// and is followed by indices of the adjacent vertices
//
unsigned int *GraphGetAdjacentsTo(const Graph *g, unsigned int v)
{
  assert(v < g->numVertices);

  // Node in the list of vertices
  List *vertices = g->verticesList;
  ListMove(vertices, v);
  struct _Vertex *vPointer = ListGetCurrentItem(vertices);
  unsigned int numAdjVertices = vPointer->outDegree;

  unsigned int *adjacent =
      (unsigned int *)calloc(1 + numAdjVertices, sizeof(unsigned int));

  if (numAdjVertices > 0)
  {
    adjacent[0] = numAdjVertices;
    List *adjList = vPointer->edgesList;
    ListMoveToHead(adjList);
    for (unsigned int i = 0; i < numAdjVertices; ListMoveToNext(adjList), i++)
    {
      struct _Edge *ePointer = ListGetCurrentItem(adjList);
      adjacent[i + 1] = ePointer->adjVertex;
    }
  }

  return adjacent;
}

//
// returns an array of size (outDegree + 1)
// element 0, stores the number of adjacent vertices
// and is followed by the distances to the adjacent vertices
//
double *GraphGetDistancesToAdjacents(const Graph *g, unsigned int v)
{
  assert(v < g->numVertices);

  // Node in the list of vertices
  List *vertices = g->verticesList;
  ListMove(vertices, v);
  struct _Vertex *vPointer = ListGetCurrentItem(vertices);
  unsigned int numAdjVertices = vPointer->outDegree;

  double *distance = (double *)calloc(1 + numAdjVertices, sizeof(double));

  if (numAdjVertices > 0)
  {
    distance[0] = numAdjVertices;
    List *adjList = vPointer->edgesList;
    ListMoveToHead(adjList);
    for (unsigned int i = 0; i < numAdjVertices; ListMoveToNext(adjList), i++)
    {
      struct _Edge *ePointer = ListGetCurrentItem(adjList);
      distance[i + 1] = ePointer->weight;
    }
  }

  return distance;
}

//
// For a graph
//
unsigned int GraphGetVertexDegree(Graph *g, unsigned int v)
{
  assert(g->isDigraph == 0);
  assert(v < g->numVertices);

  ListMove(g->verticesList, v);
  struct _Vertex *p = ListGetCurrentItem(g->verticesList);

  return p->outDegree;
}

//
// For a digraph
//
unsigned int GraphGetVertexOutDegree(Graph *g, unsigned int v)
{
  assert(g->isDigraph == 1);
  assert(v < g->numVertices);

  ListMove(g->verticesList, v);
  struct _Vertex *p = ListGetCurrentItem(g->verticesList);

  return p->outDegree;
}

//
// For a digraph
//
unsigned int GraphGetVertexInDegree(Graph *g, unsigned int v)
{
  assert(g->isDigraph == 1);
  assert(v < g->numVertices);

  ListMove(g->verticesList, v);
  struct _Vertex *p = ListGetCurrentItem(g->verticesList);

  return p->inDegree;
}

// Edges

static int _addEdge(Graph *g, unsigned int v, unsigned int w, double weight)
{
  struct _Edge *edge = (struct _Edge *)malloc(sizeof(struct _Edge));
  edge->adjVertex = w;
  edge->weight = weight;

  ListMove(g->verticesList, v);
  struct _Vertex *vertex = ListGetCurrentItem(g->verticesList);
  int result = ListInsert(vertex->edgesList, edge);

  if (result == -1)
  {
    return 0;
  }
  else
  {
    g->numEdges++;
    vertex->outDegree++;

    ListMove(g->verticesList, w);
    struct _Vertex *destVertex = ListGetCurrentItem(g->verticesList);
    destVertex->inDegree++;
  }

  if (g->isDigraph == 0)
  {
    // Bidirectional edge
    struct _Edge *edge = (struct _Edge *)malloc(sizeof(struct _Edge));
    edge->adjVertex = v;
    edge->weight = weight;

    ListMove(g->verticesList, w);
    struct _Vertex *vertex = ListGetCurrentItem(g->verticesList);
    result = ListInsert(vertex->edgesList, edge);

    if (result == -1)
    {
      return 0;
    }
    else
    {
      // g->numEdges++; // Do not count the same edge twice on a undirected
      // graph !!
      vertex->outDegree++;
    }
  }

  return 1;
}

int GraphAddEdge(Graph *g, unsigned int v, unsigned int w)
{
  assert(g->isWeighted == 0);
  assert(v != w);
  assert(v < g->numVertices);
  assert(w < g->numVertices);

  return _addEdge(g, v, w, 1.0);
}

int GraphAddWeightedEdge(Graph *g, unsigned int v, unsigned int w,double weight)
{
  assert(g->isWeighted == 1);
  assert(v != w);
  assert(v < g->numVertices);
  assert(w < g->numVertices);

  return _addEdge(g, v, w, weight);
}

int GraphRemoveEdge(Graph *g, unsigned int v, unsigned int w)
{
  assert(g != NULL);
  // TO BE COMPLETED !!

  ListMove(g->verticesList, v); // Mover para o vertice de origem da aresta
  struct _Vertex *vertex = ListGetCurrentItem(g->verticesList); // Vertice de origem da aresta

  if (vertex == NULL) {
    return 0;
  }

  ListMoveToHead(vertex->edgesList); 
  while (ListGetCurrentItem(vertex->edgesList) != NULL) 
  {
    struct _Edge *edge = ListGetCurrentItem(vertex->edgesList);

    if (edge->adjVertex == w)
    {
      ListRemoveCurrent(vertex->edgesList); // Remover a aresta da lista
      free(edge); // Libertar a memoria alocada para a aresta 
      g->numEdges--;   // Decrementar o numero de arestas do grafo 
      vertex->outDegree--;  // Decrementar o outDegree do vertice 
      break;    // Sair do ciclo porque a aresta foi encontrada
    }

    ListMoveToNext(vertex->edgesList); // Mover para a proxima aresta da lista
  }

  ListMove(g->verticesList, w); // Mover para o vertice de destino da aresta
  struct _Vertex *destVertex = ListGetCurrentItem(g->verticesList); // Vertice de destino da aresta

  if (destVertex == NULL) {
    return 0;
  }

  if (g->isDigraph == 0)
  {
    ListMoveToHead(destVertex->edgesList);
    while (ListGetCurrentItem(destVertex->edgesList) != NULL)
    {
      struct _Edge *edge = ListGetCurrentItem(destVertex->edgesList);

      if (edge->adjVertex == v)
      {
        ListRemoveCurrent(destVertex->edgesList); 
        free(edge);                       
        destVertex->outDegree--;        
        break;                             
      }

      ListMoveToNext(destVertex->edgesList); 
    }
  }
  else{
    destVertex->inDegree--;
  }

  return 1;
}

// CHECKING

int GraphCheckInvariants(const Graph *g)
{
  assert(g != NULL);
  // TO BE COMPLETED !!
  unsigned int EdgesCount = 0, listSize, vertixEdgeListSize, i, j;

  listSize = ListGetSize(g->verticesList);

  if(listSize != g->numVertices){
    return 0;
  }

  ListMoveToHead(g->verticesList);
  for(i=0; i<listSize; i++){
    struct _Vertex *Vertix = ListGetCurrentItem(g->verticesList);
    if(Vertix->id != i){
      return 0;
    }
    vertixEdgeListSize = ListGetSize(Vertix->edgesList);
    if(vertixEdgeListSize != Vertix->outDegree){
      return 0;
    }
    EdgesCount += vertixEdgeListSize;
    for(j=0; j<vertixEdgeListSize; j++){
      struct _Edge *Edge = ListGetCurrentItem(Vertix->edgesList);
      if(Edge->adjVertex == Vertix->id){
        return 0;
      }
      ListMoveToNext(Vertix->edgesList);
    }
    ListMoveToNext(g->verticesList);    
  }

  if(g->isDigraph == 0){
    if(EdgesCount != g->numEdges){
      return 0;
    }
  }
  else{
    if(EdgesCount != g->numEdges*2){
      return 0;
    }
  }

  return 1;
}

// DISPLAYING on the console

void GraphDisplay(const Graph *g)
{
  printf("---\n");
  if (g->isWeighted)
  {
    printf("Weighted ");
  }
  if (g->isComplete)
  {
    printf("COMPLETE ");
  }
  if (g->isDigraph)
  {
    printf("Digraph\n");
    printf("Max Out-Degree = %d\n", GraphGetMaxOutDegree(g));
  }
  else
  {
    printf("Graph\n");
    printf("Max Degree = %d\n", GraphGetMaxDegree(g));
  }
  printf("Vertices = %2d | Edges = %2d\n", g->numVertices, g->numEdges);

  List *vertices = g->verticesList;
  ListMoveToHead(vertices);
  unsigned int i = 0;
  for (; i < g->numVertices; ListMoveToNext(vertices), i++)
  {
    printf("%2d ->", i);
    struct _Vertex *v = ListGetCurrentItem(vertices);
    if (ListIsEmpty(v->edgesList))
    {
      printf("\n");
    }
    else
    {
      List *edges = v->edgesList;
      unsigned int i = 0;
      ListMoveToHead(edges);
      for (; i < ListGetSize(edges); ListMoveToNext(edges), i++)
      {
        struct _Edge *e = ListGetCurrentItem(edges);
        if (g->isWeighted)
        {
          printf("   %2d(%4.2f)", e->adjVertex, e->weight);
        }
        else
        {
          printf("   %2d", e->adjVertex);
        }
      }
      printf("\n");
    }
  }
  printf("---\n");
}

void GraphListAdjacents(const Graph *g, unsigned int v)
{
  printf("---\n");

  unsigned int *array = GraphGetAdjacentsTo(g, v);

  printf("Vertex %d has %d adjacent vertices -> ", v, array[0]);

  for (unsigned int i = 1; i <= array[0]; i++)
  {
    printf("%d ", array[i]);
  }

  printf("\n");

  free(array);

  printf("---\n");
}
