#ifndef GRAPH_H
#define GRAPH_H

#include <stdio.h>
#include <stdlib.h>

#ifdef SYSV
#include <string.h>
#else
#include <strings.h>
#endif

#undef min  // Avoid conflicts with other definitions

/* Utility union to hold different types of pointers and values */
typedef union {
    struct vertex_struct *V;
    struct arc_struct *A;
    struct graph_struct *G;
    char *S;
    long I;
} util;

/* Vertex structure representing a node in the graph */
typedef struct vertex_struct {
    struct arc_struct *arcs;  // Pointer to the list of arcs
    char *name;               // Name of the vertex
    util u, v, w, x, y, z;    // Utility fields for additional data
} Vertex;

/* Arc structure representing an edge in the graph */
typedef struct arc_struct {
    struct vertex_struct *tip; // Pointer to the destination vertex
    struct arc_struct *next;   // Pointer to the next arc in the list
    long len;                  // Length of the arc
    util a, b;                 // Utility fields for additional data
} Arc;

/* Macro to initialize an area */
#define init_area(s) *s = NULL

/* Structure for managing memory areas */
struct area_pointers {
    char *first;                // Pointer to the first allocated memory block
    struct area_pointers *next; // Pointer to the next area
};

/* Typedef for an array of area pointers */
typedef struct area_pointers *Area[1];

/* Constants */
#define ID_FIELD_SIZE 161

/* Graph structure representing the entire graph */
typedef struct graph_struct {
    Vertex *vertices;          // Pointer to the array of vertices
    long n;                    // Number of vertices
    long m;                    // Number of arcs
    char id[ID_FIELD_SIZE];    // Identifier for the graph
    char util_types[15];       // Types for utility fields
    Area data;                 // Data area
    Area aux_data;             // Auxiliary data area
    util uu, vv, ww, xx, yy, zz; // Utility fields for additional data
} Graph;

/* Typedef for size type */
typedef unsigned long siz_t;

/* External variables */
extern long verbose;
extern long panic_code;
extern long trouble_code;
extern siz_t edge_trick;

/* Error codes */
#define alloc_fault       (-1)  // a previous memory request failed
#define no_room            1    // the current memory request failed
#define early_data_fault  10    // error detected at beginning of .dat file
#define late_data_fault   11    //  error detected at end of .dat file
#define syntax_error      20    // error detected while reading .dat file
#define bad_specs         30    // parameter out of range or otherwise disallowed
#define very_bad_specs    40    // parameter far out or otherwise stupid
#define missing_operand   50    // graph parameter is null
#define invalid_operand   60    // graph parameter doesn't obey assumptions
#define impossible        90    // this can't happen

/* Function declarations */
extern char *graph_alloc(long n, Area s);
extern void graph_free(Area s);

extern long extra_n;
extern char null_string[];

#define new_graph nugraph
#define new_arc nuarc
#define new_edge nuedge

extern Graph *new_graph(long n);
extern void new_arc(Vertex *u, Vertex *v, long len);
extern Arc *virgin_arc();
extern void new_edge(Vertex *u, Vertex *v, long len);
extern char *save_string(char *s);
extern void switch_to_graph(Graph *g);
extern void recycle(Graph *g);

/* Hashing function declarations */
extern void hash_in(Vertex *v);
extern Vertex *hash_out(char *s);
extern void hash_setup(Graph *g);
extern Vertex *hash_lookup(char *s, Graph *g);
extern void make_double_compound_id(Graph *g, char *s1, Graph *gg, char *s2, Graph *ggg, char *s3);
extern void make_compound_id(Graph *g, char *s1, Graph *gg, char *s2);

/* Macros for graph operations */
#define n_1  uu.I
#define mark_bipartite(g, n1) do { \
    g->n_1 = n1; \
    g->util_types[8] = 'I'; \
} while (0)

#define typed_alloc(n, t, s) (t*) graph_alloc((long) ((n) * sizeof(t)), s)

#endif // GRAPH_H