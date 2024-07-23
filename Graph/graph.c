//
// Created by jorda on 7/17/2024.
//

#ifndef GRAPH_H
#define GRAPH_H

#ifdef SYSV
#include <string.h>
#else
#include <strings.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Macros for memory allocation and graph operations */
#define typed_alloc(n, t, s) (t*) graph_alloc((long) ((n) * sizeof(t)), s)

#define n_1 uu.I
#define arcs_per_block 102
#define new_graph nugraph
#define new_arc nuarc
#define new_edge nuedge
#define string_block_size 1016

/* Macros for hashing operations */
#define hash_link u.V
#define hash_head v.V
#define HASH_MULT 314159
#define HASH_PRIME 516595003

/* Utility union to hold different types of pointers and values */
typedef union {
    struct vertex_struct *V;    // pointer to Vertex
    struct arc_struct *A;       // pointer to Arc
    struct graph_struct *G;     // pointer to Graph
    char *S;                    // pointer to string
    long I;                     // integer
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
typedef unsigned long siz_t;    // basic machine address, as signless integer

/* External variables */
extern long verbose;
extern long panic_code;
extern long trouble_code;
extern siz_t edge_trick;

/* Error codes */
#define alloc_fault       (-1)
#define no_room            1
#define early_data_fault  10
#define late_data_fault   11
#define syntax_error      20
#define bad_specs         30
#define very_bad_specs    40
#define missing_operand   50
#define invalid_operand   60
#define impossible        90

/* Function declarations */
extern char *graph_alloc(long n, Area s);
extern void graph_free(Area s);

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

/* Static variables */
static Arc *next_arc;                   // the next Arc available for allocation
static Arc *bad_arc;                    // but if next_arc = bad_arc, that Arc isn't there
static char *next_string;               // the next byte available for storing strings
static char *bad_string;                // but if next_string = bad_string, don't byte
static Arc dummy_arc[2];                // an Arc record to point to in an emergency
static Graph dummy_graph;               // a Graph record that's normally unused
static Graph *cur_graph = &dummy_graph; // the Graph most recently created

/* Default values */
long verbose = 0;
long panic_code = 0;    // nonzero  if "verbose" output is desired
long trouble_code = 0;  // set nonzero if graph generator returns null pointer
long extra_n = 4;
char null_string[1] = { '\0' };
siz_t edge_trick = sizeof(Arc) - (sizeof(Arc) & (sizeof(Arc) - 1)); // least significant bit in sizeof(Arc)

/* Memory allocation function */
char* graph_alloc(long n, Area s) {
    if (n <= 0 || n > 0xffff00 - 2 * sizeof(char*)) {
        fprintf(stderr, "Error: Illegal memory allocation request.\n");
        trouble_code |= 2;
        return NULL;
    }

    long m = sizeof(char*);
    n = ((n + m - 1) / m) * m;

    char *loc = (char*) calloc((unsigned)((n + 2 * m + 255) / 256), 256);
    if (loc) {
        struct area_pointers *t = (struct area_pointers*)(loc + n);
        t->first = loc;
        t->next = *s;
        *s = t;
    } else {
        fprintf(stderr, "Error: Memory allocation failed.\n");
        trouble_code |= 1;
        return NULL;
    }
    return loc;
}

/* Memory graph_free function */
void graph_free(Area s) {
    while (*s) {
        struct area_pointers *t = (*s)->next;
        free((*s)->first);
        *s = t;
    }
}

/* Function to create a new graph */
Graph* new_graph(long n) {
    if (n <= 0) {
        fprintf(stderr, "Error: Number of vertices must be positive.\n");
        return NULL;
    }

    cur_graph = (Graph*) calloc(1, sizeof(Graph));
    if (!cur_graph) {
        fprintf(stderr, "Error: Memory allocation for graph failed.\n");
        return NULL;
    }

    cur_graph->vertices = typed_alloc(n + extra_n, Vertex, cur_graph->data);
    if (!cur_graph->vertices) {
        fprintf(stderr, "Error: Memory allocation for vertices failed.\n");
        free(cur_graph);
        cur_graph = NULL;
        return NULL;
    }

    Vertex *p;
    cur_graph->n = n;
    for (p = cur_graph->vertices + n + extra_n - 1; p >= cur_graph->vertices; p--) {
        p->name = null_string;
    }

    snprintf(cur_graph->id, ID_FIELD_SIZE, "new_graph(%ld)", n);
    strcpy(cur_graph->util_types, "ZZZZZZZZZZZZZZ");

    next_arc = bad_arc = NULL;
    next_string = bad_string = NULL;
    trouble_code = 0;

    return cur_graph;
}

/* Function to create a compound ID */
void make_compound_id(Graph *g, char *s1, Graph *gg, char *s2) {
    if (!g || !s1 || !gg || !s2) {
        fprintf(stderr, "Error: Invalid parameters for make_compound_id.\n");
        return;
    }

    int avail = ID_FIELD_SIZE - strlen(s1) - strlen(s2);
    char tmp[ID_FIELD_SIZE];
    strncpy(tmp, gg->id, ID_FIELD_SIZE);

    if (strlen(tmp) < avail) {
        snprintf(g->id, ID_FIELD_SIZE, "%s%s%s", s1, tmp, s2);
    } else {
        snprintf(g->id, ID_FIELD_SIZE, "%s%.*s...)%s", s1, avail - 5, tmp, s2);
    }
}

/* Function to create a double compound ID */
void make_double_compound_id(Graph *g, char *s1, Graph *gg, char *s2, Graph *ggg, char *s3) {
    if (!g || !s1 || !gg || !s2 || !ggg || !s3) {
        fprintf(stderr, "Error: Invalid parameters for make_double_compound_id.\n");
        return;
    }

    int avail = ID_FIELD_SIZE - strlen(s1) - strlen(s2) - strlen(s3);
    if (strlen(gg->id) + strlen(ggg->id) < avail) {
        snprintf(g->id, ID_FIELD_SIZE, "%s%s%s%s%s", s1, gg->id, s2, ggg->id, s3);
    } else {
        snprintf(g->id, ID_FIELD_SIZE, "%s%.*s...)%s%.*s...)%s", s1, avail / 2 - 5, gg->id, s2, (avail - 9) / 2, ggg->id, s3);
    }
}

/* Function to get a new arc */
Arc* virgin_arc() {
    if (!cur_graph) {
        fprintf(stderr, "Error: No current graph set.\n");
        return NULL;
    }

    if (next_arc == bad_arc) {
        next_arc = typed_alloc(arcs_per_block, Arc, cur_graph->data);
        if (!next_arc) {
            next_arc = dummy_arc;
        } else {
            bad_arc = next_arc + arcs_per_block;
        }
    }

    Arc *cur_arc = next_arc;
    next_arc++;

    return cur_arc;
}

/* Function to create a new arc */
void new_arc(Vertex *u, Vertex *v, long len) {
    if (!u || !v) {
        fprintf(stderr, "Error: Invalid vertices provided.\n");
        return;
    }

    Arc *cur_arc = virgin_arc();
    if (!cur_arc) {
        fprintf(stderr, "Error: Failed to create a new arc.\n");
        return;
    }

    cur_arc->tip = v;
    cur_arc->next = u->arcs;
    cur_arc->len = len;
    u->arcs = cur_arc;
    cur_graph->m++;
}

/* Function to create a new edge */
void new_edge(Vertex *u, Vertex *v, long len) {
    if (!u || !v) {
        fprintf(stderr, "Error: Invalid vertices provided.\n");
        return;
    }

    Arc *cur_arc = virgin_arc();
    if (!cur_arc || cur_arc == dummy_arc) {
        fprintf(stderr, "Error: Failed to create a new edge.\n");
        return;
    }

    next_arc++;
    if (u < v) {
        cur_arc->tip = v;
        cur_arc->next = u->arcs;
        (cur_arc + 1)->tip = u;
        (cur_arc + 1)->next = v->arcs;
        u->arcs = cur_arc;
        v->arcs = cur_arc + 1;
    } else {
        (cur_arc + 1)->tip = v;
        (cur_arc + 1)->next = u->arcs;
        u->arcs = cur_arc + 1;
        cur_arc->tip = u;
        cur_arc->next = v->arcs;
        v->arcs = cur_arc;
    }

    cur_arc->len = (cur_arc + 1)->len = len;
    cur_graph->m += 2;
}

/* Function to save a string */
char* save_string(char *s) {
    if (!s) {
        fprintf(stderr, "Error: Null string provided.\n");
        return null_string;
    }

    char *p = s;
    long len;
    while (*p++) ;
    len = p - s;

    p = next_string;
    if (p + len > bad_string) {
        long size = string_block_size;
        if (len > size) {
            size = len;
        }

        p = graph_alloc(size, cur_graph->data);
        if (!p) {
            fprintf(stderr, "Error: Memory allocation for string failed.\n");
            return null_string;
        }

        bad_string = p + size;
    }

    while (*s) {
        *p++ = *s++;
    }
    *p++ = '\0';

    next_string = p;
    return p - len;
}

/* Function to switch to a different graph */
void switch_to_graph(Graph *g) {
    if (cur_graph) {
        cur_graph->ww.A = next_arc;
        cur_graph->xx.A = bad_arc;
        cur_graph->yy.S = next_string;
        cur_graph->zz.S = bad_string;
    }

    cur_graph = (g ? g : &dummy_graph);
    next_arc = cur_graph->ww.A;
    bad_arc = cur_graph->xx.A;
    next_string = cur_graph->yy.S;
    bad_string = cur_graph->zz.S;

    cur_graph->ww.A = NULL;
    cur_graph->xx.A = NULL;
    cur_graph->yy.S = NULL;
    cur_graph->zz.S = NULL;
}

/* Function to recycle a graph */
void recycle(Graph *g) {
    if (g) {
        graph_free(g->data);
        graph_free(g->aux_data);
        free(g);
    }
}

/* Hashing function to insert a vertex */
void hash_in(Vertex *v) {
    if (!v || !v->name) {
        fprintf(stderr, "Error: Invalid vertex or vertex name provided.\n");
        return;
    }

    char *t = v->name;
    long h = 0;

    while (*t) {
        h += (h ^ (h >> 1)) + HASH_MULT * (unsigned char) *t++;
        while (h >= HASH_PRIME) {
            h -= HASH_PRIME;
        }
    }

    Vertex *u = cur_graph->vertices + (h % cur_graph->n);
    v->hash_link = u->hash_head;
    u->hash_head = v;
}

/* Hashing function to look up a vertex by name */
Vertex* hash_out(char *s) {
    if (!s) {
        fprintf(stderr, "Error: Null string provided.\n");
        return NULL;
    }

    char *t = s;
    long h = 0;

    while (*t) {
        h += (h ^ (h >> 1)) + HASH_MULT * (unsigned char) *t++;
        while (h >= HASH_PRIME) {
            h -= HASH_PRIME;
        }
    }

    Vertex *u = cur_graph->vertices + (h % cur_graph->n);
    for (u = u->hash_head; u; u = u->hash_link) {
        if (strcmp(s, u->name) == 0) {
            return u;
        }
    }

    return NULL;
}

/* Function to set up the hash table for a graph */
void hash_setup(Graph *g) {
    if (!g || g->n <= 0) {
        fprintf(stderr, "Error: Invalid graph or number of vertices.\n");
        return;
    }

    Graph *save_cur_graph = cur_graph;
    cur_graph = g;

    for (Vertex *v = g->vertices; v < g->vertices + g->n; v++) {
        v->hash_head = NULL;
    }

    for (Vertex *v = g->vertices; v < g->vertices + g->n; v++) {
        hash_in(v);
    }

    g->util_types[0] = g->util_types[1] = 'V';
    cur_graph = save_cur_graph;
}

/* Function to look up a vertex by name in a specified graph */
Vertex* hash_lookup(char *s, Graph *g) {
    if (!g || g->n <= 0 || !s) {
        fprintf(stderr, "Error: Invalid graph, number of vertices, or string provided.\n");
        return NULL;
    }

    Graph *save_cur_graph = cur_graph;
    cur_graph = g;
    Vertex *v = hash_out(s);
    cur_graph = save_cur_graph;

    return v;
}

#endif // GRAPH_H
