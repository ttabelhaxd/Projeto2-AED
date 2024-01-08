//
// Algoritmos e Estruturas de Dados --- 2023/2024
//
// Topological Sorting
//

#include "GraphTopologicalSorting.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "Graph.h"
#include "IntegersQueue.h"
#include "instrumentation.h"

#define VERTEXACCESS InstrCount[0]
#define ADJVERTEX InstrCount[1] 
#define QUEUEACCESS InstrCount[2]
#define NUMITERATIONS InstrCount[3]


struct _GraphTopoSort {
  int* marked;                     // Aux array
  unsigned int* numIncomingEdges;  // Aux array
  unsigned int* vertexSequence;    // The result
  int validResult;                 // 0 or 1
  unsigned int numVertices;        // From the graph
  Graph* graph;
};

// AUXILIARY FUNCTION
// Allocate memory for the struct
// And for its array fields
// Initialize all struct fields
//
static GraphTopoSort* _create(Graph* g) {
  assert(g != NULL);

  // TO BE COMPLETED
  
  GraphTopoSort* p = (GraphTopoSort*)malloc(sizeof(GraphTopoSort));
  if(p == NULL) {
    return NULL;
  }
  p->numVertices = GraphGetNumVertices(g);
  p->marked = (int*)calloc(sizeof(int), p->numVertices);
  if(p->marked == NULL) {
    free(p);
    return NULL;
  }
  p->numIncomingEdges = (unsigned int*)calloc(sizeof(unsigned int), p->numVertices);
  if(p->numIncomingEdges == NULL) {
    free(p->marked);
    free(p);
    return NULL;
  }
  p->vertexSequence = (unsigned int*)calloc(sizeof(unsigned int), p->numVertices);
  if(p->vertexSequence == NULL) {
    free(p->marked);
    free(p->numIncomingEdges);
    free(p);
    return NULL;
  }
  p->validResult = 0;
  p->graph = g;

  return p;
}

//
// Computing the topological sorting, if any, using the 1st algorithm:
// 1 - Create a copy of the graph
// 2 - Successively identify vertices without incoming edges and remove their
//     outgoing edges
// Check if a valid sorting was computed and set the isValid field
// For instance, by checking if the number of elements in the vertexSequence is
// the number of graph vertices
//

//A complexidade deste algoritmo é O(V^2 + E)
GraphTopoSort* GraphTopoSortComputeV1(Graph* g) {
  assert(g != NULL && GraphIsDigraph(g) == 1);

  // Create and initialize the struct

  GraphTopoSort* topoSort = _create(g);

  // Build the topological sorting

  // TO BE COMPLETED
  unsigned int vCounter = 0;
  Graph* copy = GraphCopy(g); //O(V^2)
  if(copy == NULL) {
    GraphTopoSortDestroy(&topoSort);
    return NULL;
  }

  unsigned int vertex = 0; 
  while(vertex < topoSort->numVertices) { //O(E)
    NUMITERATIONS++;
    VERTEXACCESS++;
    if(GraphGetVertexInDegree(copy, vertex) == 0 && !topoSort->marked[vertex]) {
      topoSort->marked[vertex] = 1;
      topoSort->vertexSequence[vCounter] = vertex;
      vCounter++;
      unsigned int *adjVertex = GraphGetAdjacentsTo(copy,vertex);
      for(unsigned int i = 1; i <= adjVertex[0]; i++){
          NUMITERATIONS++;
          unsigned int v = adjVertex[i];
          ADJVERTEX++;
          GraphRemoveEdge(copy, vertex, v);
      }
      free(adjVertex);
      vertex = 0;
    }
    else {
      vertex++;
    }

  }
  
  if(vCounter == topoSort->numVertices) {
    topoSort->validResult = 1;
  }
  GraphDestroy(&copy);

  return topoSort;
}

//
// Computing the topological sorting, if any, using the 2nd algorithm
// Check if a valid sorting was computed and set the isValid field
// For instance, by checking if the number of elements in the vertexSequence is
// the number of graph vertices
//

//A complexidade deste algoritmo é O(V^2 + E)
GraphTopoSort* GraphTopoSortComputeV2(Graph* g) {
  assert(g != NULL && GraphIsDigraph(g) == 1);

  // Create and initialize the struct

  GraphTopoSort* topoSort = _create(g);

  // Build the topological sorting

  // TO BE COMPLETED
  unsigned int vCounter = 0;
  unsigned int vertex;

  for(vertex = 0; vertex < topoSort->numVertices; vertex++) { //O(E)
    topoSort->numIncomingEdges[vertex] = GraphGetVertexInDegree(g, vertex);
    NUMITERATIONS++;
  }

  vertex = 0;

  while(vertex < topoSort->numVertices) { //O(V)
    VERTEXACCESS++;
    NUMITERATIONS++;
    if(topoSort->numIncomingEdges[vertex] == 0 && !topoSort->marked[vertex]) {
      topoSort->marked[vertex] = 1;
      topoSort->vertexSequence[vCounter] = vertex;
      vCounter++;
      unsigned int *adjVertex = GraphGetAdjacentsTo(g,vertex);
      for(unsigned int i = 1; i <= adjVertex[0]; i++){ //O(V)
          NUMITERATIONS++;
          unsigned int v = adjVertex[i];
          ADJVERTEX++;
          topoSort->numIncomingEdges[v]--;
      }
      free(adjVertex);
      vertex = 0;
    }
    else {
      vertex++;
    }

  }
  
  if(vCounter == topoSort->numVertices) {
    topoSort->validResult = 1;
  }

  return topoSort;
}

//
// Computing the topological sorting, if any, using the 3rd algorithm
// Check if a valid sorting was computed and set the isValid field
// For instance, by checking if the number of elements in the vertexSequence is
// the number of graph vertices
//

//A complexidade deste algoritmo é O(V^2 + E)
GraphTopoSort* GraphTopoSortComputeV3(Graph* g) {
  assert(g != NULL && GraphIsDigraph(g) == 1);

  // Create and initialize the struct

  GraphTopoSort* topoSort = _create(g);

  // Build the topological sorting

  // TO BE COMPLETED
  unsigned int vCounter = 0;
  unsigned int vertex;

  for(vertex = 0; vertex < topoSort->numVertices; vertex++) { //O(E)
    topoSort->numIncomingEdges[vertex] = GraphGetVertexInDegree(g, vertex);
    NUMITERATIONS++;
  }

  Queue *queue = QueueCreate(topoSort->numVertices);
  if(queue == NULL) {
    GraphTopoSortDestroy(&topoSort);
    return NULL;
  }

  for(vertex = 0; vertex < topoSort->numVertices; vertex++) { //O(V)
    NUMITERATIONS++;
    if(topoSort->numIncomingEdges[vertex] == 0) {
      QueueEnqueue(queue, vertex);
      QUEUEACCESS++;
      
    }
  }

  while(!QueueIsEmpty(queue)) { //O(V)
    vertex = QueueDequeue(queue);
    QUEUEACCESS++;
    topoSort->vertexSequence[vCounter] = vertex;
    vCounter++;
    VERTEXACCESS++;
    NUMITERATIONS++;
    unsigned int *adjVertex = GraphGetAdjacentsTo(g,vertex);
    for(unsigned int i = 1; i <= adjVertex[0]; i++){ //O(V)
        unsigned int v = adjVertex[i];
        ADJVERTEX++;
        NUMITERATIONS++;
        topoSort->numIncomingEdges[v]--;
        
        if(topoSort->numIncomingEdges[v] == 0) {
          QueueEnqueue(queue, v);
          QUEUEACCESS++;
        }
    }
    free(adjVertex);
  }

  QueueDestroy(&queue);

  if(vCounter == topoSort->numVertices) {
    topoSort->validResult = 1;
  }

  return topoSort;
}

void GraphTopoSortDestroy(GraphTopoSort** p) {
  assert(*p != NULL);

  GraphTopoSort* aux = *p;

  free(aux->marked);
  free(aux->numIncomingEdges);
  free(aux->vertexSequence);

  free(*p);
  *p = NULL;
}

//
// A valid sorting was computed?
//
int GraphTopoSortIsValid(const GraphTopoSort* p) { return p->validResult; }

//
// Getting an array containing the topological sequence of vertex indices
// Or NULL, if no sequence could be computed
// MEMORY IS ALLOCATED FOR THE RESULTING ARRAY
//
unsigned int* GraphTopoSortGetSequence(const GraphTopoSort* p) {
  assert(p != NULL);
  // TO BE COMPLETED
  if(!p->validResult) {
    return NULL;
  }
  unsigned int* vertexSequence = (unsigned int*)malloc(sizeof(unsigned int) * p->numVertices);
  if(vertexSequence == NULL) {
    return NULL;
  }
  for(unsigned int vertex = 0; vertex < p->numVertices; vertex++) {
    vertexSequence[vertex] = p->vertexSequence[vertex];
  }
  return vertexSequence;
}

// DISPLAYING on the console

//
// The toplogical sequence of vertex indices, if it could be computed
//
void GraphTopoSortDisplaySequence(const GraphTopoSort* p) {
  assert(p != NULL);

  if (p->validResult == 0) {
    printf(" *** The topological sorting could not be computed!! *** \n");
    return;
  }

  printf("Topological Sorting - Vertex indices:\n");
  for (unsigned int i = 0; i < GraphGetNumVertices(p->graph); i++) {
    printf("%d ", p->vertexSequence[i]);
  }
  printf("\n");
}

//
// The toplogical sequence of vertex indices, if it could be computed
// Followed by the digraph displayed using the adjecency lists
// Adjacency lists are presented in topologic sorted order
//
void GraphTopoSortDisplay(const GraphTopoSort* p) {
  assert(p != NULL);

  // The topological order
  GraphTopoSortDisplaySequence(p);

  if (p->validResult == 0) {
    return;
  }

  // The Digraph
  for (unsigned int i = 0; i < GraphGetNumVertices(p->graph); i++) {
    GraphListAdjacents(p->graph, p->vertexSequence[i]);
  }
  printf("\n");
}
