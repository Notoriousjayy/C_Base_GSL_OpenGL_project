#include <string.h>
#include "../include/graph.h"
#include "../include/graph_basic.h"

/* Macros for error handling and buffer sizes */
#define panic(c) \
{ \
    panic_code = c; \
    graph_free(working_storage); \
    trouble_code = 0; \
    return NULL; \
}

#define BUF_SIZE 4096
#define MAX_D 91
#define MAX_NNN 1000000000.0
#define UL_BITS 8 * sizeof(unsigned long)
#define vert_offset(v, delta) ((Vertex*) (((siz_t) v) + delta))

/* Macros for temporary variables */
#define tmp u.V
#define tlen z.A
#define mult v.I
#define minlen w.I
#define map z.V
#define ind z.I
#define IND_GRAPH 1000000000
#define subst y.G

/* Static variables */
static Area working_storage;
static char buffer[BUF_SIZE];
static long nn[MAX_D + 1];
static long wr[MAX_D + 1];
static long del[MAX_D + 1];
static long sig[MAX_D + 2];
static long xx[MAX_D + 1], yy[MAX_D + 1];
static char* short_imap = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_^~&@,;.:?!%#$+-*/|<=>()[]{}`'";

/* Placeholder function declarations */
void nuedge(Vertex *u, Vertex *v, long len);
void nuarc(Vertex *u, Vertex *v, long len);
Graph* nugraph(long n);

/* Function implementations */

Graph* board(long n1, long n2, long n3, long n4, long piece, long wrap, long directed) {
    Graph* g;
    register long i, j, k;
    register long d;
    register Vertex* v;
    register long s;
    long n;
    long p;
    long l;

    if (piece == 0) piece = 1;
    if (n1 <= 0) {
        n1 = n2 = 8;
        n3 = 0;
    }
    nn[1] = n1;
    if (n2 <= 0) {
        k = 2;
        d = -n2;
        n3 = n4 = 0;
    } else {
        nn[2] = n2;
        if (n3 <= 0) {
            k = 3;
            d = -n3;
            n4 = 0;
        } else {
            nn[3] = n3;
            if (n4 <= 0) {
                k = 4;
                d = -n4;
            } else {
                nn[4] = n4;
                d = 4;
                goto done;
            }
        }
    }
    if (d == 0) {
        d = k - 1;
        goto done;
    }

    if (d > MAX_D) panic(bad_specs);
    for (j = 1; k <= d; j++, k++) nn[k] = nn[j];

    done:

    {
        float nnn;
        for (n = 1, nnn = 1.0, j = 1; j <= d; j++) {
            nnn *= (float)nn[j];
            if (nnn > MAX_NNN) panic(very_bad_specs);
            n *= nn[j];
        }
        g = new_graph(n);
        if (g == NULL)
        panic(no_room);
        sprintf(g->id, "board(%ld,%ld,%ld,%ld,%ld,%ld,%d)",
                n1, n2, n3, n4, piece, wrap, directed ? 1 : 0);
        strcpy(g->util_types, "ZZZIIIZZZZZZZZ");

        {
            register char* q;
            nn[0] = xx[0] = xx[1] = xx[2] = xx[3] = 0;
            for (k = 4; k <= d; k++) xx[k] = 0;
            for (v = g->vertices;; v++) {
                q = buffer;
                for (k = 1; k <= d; k++) {
                    sprintf(q, ".%ld", xx[k]);
                    while (*q) q++;
                }
                v->name = save_string(&buffer[1]);
                v->x.I = xx[1];
                v->y.I = xx[2];
                v->z.I = xx[3];
                for (k = d; xx[k] + 1 == nn[k]; k--) xx[k] = 0;
                if (k == 0) break;
                xx[k]++;
            }
        }
    }

    {
        register long w = wrap;
        for (k = 1; k <= d; k++, w >>= 1) {
            wr[k] = w & 1;
            del[k] = sig[k] = 0;
        }
        sig[0] = del[0] = sig[d + 1] = 0;
    }
    p = piece;
    if (p < 0) p = -p;
    while (1) {
        for (k = d; sig[k] + (del[k] + 1) * (del[k] + 1) > p; k--) del[k] = 0;
        if (k == 0) break;
        del[k]++;
        sig[k + 1] = sig[k] + del[k] * del[k];
        for (k++; k <= d; k++) sig[k + 1] = sig[k];
        if (sig[d + 1] < p) continue;

        while (1) {
            for (k = 1; k <= d; k++) xx[k] = 0;
            for (v = g->vertices;; v++) {
                for (k = 1; k <= d; k++) yy[k] = xx[k] + del[k];
                for (l = 1;; l++) {
                    for (k = 1; k <= d; k++) {
                        if (yy[k] < 0) {
                            if (!wr[k]) goto no_more;
                            do yy[k] += nn[k]; while (yy[k] < 0);
                        } else if (yy[k] >= nn[k]) {
                            if (!wr[k]) goto no_more;
                            do yy[k] -= nn[k]; while (yy[k] >= nn[k]);
                        }
                    }

                    if (piece < 0) {
                        for (k = 1; k <= d; k++) if (yy[k] != xx[k]) goto unequal;
                        goto no_more;
                        unequal:;
                    }

                    for (k = 2, j = yy[1]; k <= d; k++) j = nn[k] * j + yy[k];
                    if (directed) new_arc(v, g->vertices + j, l);
                    else new_edge(v, g->vertices + j, l);
                    if (piece > 0) goto no_more;
                    for (k = 1; k <= d; k++) yy[k] += del[k];
                }
                no_more:;
                for (k = d; xx[k] + 1 == nn[k]; k--) xx[k] = 0;
                if (k == 0) break;
                xx[k]++;
            }

            for (k = d; del[k] <= 0; k--) del[k] = -del[k];
            if (sig[k] == 0) break;
            del[k] = -del[k];
        }
    }

    if (trouble_code) {
        recycle(g);
        panic(alloc_fault);
    }
    return g;
}

Graph* simplex(unsigned long n, long n0, long n1, long n2, long n3, long n4, long directed) {
    Graph* g;
    register long i, j, k;
    register long d;
    register Vertex* v;
    register long s;

    if (n0 == 0) n0 = -2;
    if (n0 < 0) {
        k = 2;
        nn[0] = n;
        d = -n0;
        n1 = n2 = n3 = n4 = 0;
    } else {
        if (n0 > n) n0 = n;
        nn[0] = n0;
        if (n1 <= 0) {
            k = 2;
            d = -n1;
            n2 = n3 = n4 = 0;
        } else {
            if (n1 > n) n1 = n;
            nn[1] = n1;
            if (n2 <= 0) {
                k = 3;
                d = -n2;
                n3 = n4 = 0;
            } else {
                if (n2 > n) n2 = n;
                nn[2] = n2;
                if (n3 <= 0) {
                    k = 4;
                    d = -n3;
                    n4 = 0;
                } else {
                    if (n3 > n) n3 = n;
                    nn[3] = n3;
                    if (n4 <= 0) {
                        k = 5;
                        d = -n4;
                    } else {
                        if (n4 > n) n4 = n;
                        nn[4] = n4;
                        d = 4;
                        goto done;
                    }
                }
            }
        }
        if (d == 0) {
            d = k - 2;
            goto done;
        }
        nn[k - 1] = nn[0];
    }

    if (d > MAX_D) panic(bad_specs);
    for (j = 1; k <= d; j++, k++) nn[k] = nn[j];

    done:

    {
        long nverts;
        register long* coef = typed_alloc(n + 1, long, working_storage);
        if (trouble_code) panic(no_room + 1);
        for (k = 0; k <= nn[0]; k++) coef[k] = 1;

        for (j = 1; j <= d; j++) {
            for (k = n, i = n - nn[j] - 1; i >= 0; k--, i--) coef[k] -= coef[i];
            s = 1;
            for (k = 1; k <= n; k++) {
                s += coef[k];
                if (s > 1000000000) panic(very_bad_specs);
                coef[k] = s;
            }
        }

        nverts = coef[n];
        graph_free(working_storage);
        g = new_graph(nverts);
        if (g == NULL)
        panic(no_room);
    }
    sprintf(g->id, "simplex(%lu,%ld,%ld,%ld,%ld,%ld,%d)",
            n, n0, n1, n2, n3, n4, directed ? 1 : 0);
    strcpy(g->util_types, "VVZIIIZZZZZZZZ");

    v = g->vertices;
    yy[d + 1] = 0;
    sig[0] = n;
    for (k = d; k >= 0; k--) yy[k] = yy[k + 1] + nn[k];
    if (yy[0] >= n) {
        k = 0;
        xx[0] = (yy[1] >= n ? 0 : n - yy[1]);
        while (1) {
            for (s = sig[k] - xx[k], k++; k <= d; s -= xx[k], k++) {
                sig[k] = s;
                if (s <= yy[k + 1]) xx[k] = 0;
                else xx[k] = s - yy[k + 1];
            }
            if (s != 0) panic(impossible + 1);

            {
                register char* p = buffer;
                for (k = 0; k <= d; k++) {
                    sprintf(p, ".%ld", xx[k]);
                    while (*p) p++;
                }
                v->name = save_string(&buffer[1]);
                v->x.I = xx[0];
                v->y.I = xx[1];
                v->z.I = xx[2];
            }
            hash_in(v);

            for (j = 0; j < d; j++)
                if (xx[j]) {
                    register Vertex* u;
                    xx[j]--;
                    for (k = j + 1; k <= d; k++)
                        if (xx[k] < nn[k]) {
                            register char* p = buffer;
                            xx[k]++;
                            for (i = 0; i <= d; i++) {
                                sprintf(p, ".%ld", xx[i]);
                                while (*p) p++;
                            }
                            u = hash_out(&buffer[1]);
                            if (u == NULL) panic(impossible + 2);
                            if (directed) new_arc(u, v, 1L);
                            else new_edge(u, v, 1L);
                            xx[k]--;
                        }
                    xx[j]++;
                }
            v++;
            for (k = d - 1;; k--) {
                if (xx[k] < sig[k] && xx[k] < nn[k]) break;
                if (k == 0) goto last;
            }
            xx[k]++;
        }
    }
    last:
    if (v != g->vertices + g->n)
    panic(impossible);
    if (trouble_code) {
        recycle(g);
        panic(alloc_fault);
    }
    return g;
}

Graph* subsets(unsigned long n, long n0, long n1, long n2, long n3, long n4, unsigned long size_bits, long directed) {
    Graph* g;
    register long i, j, k;
    register long d;
    register Vertex* v;
    register long s;

    if (n0 == 0) n0 = -2;
    if (n0 < 0) {
        k = 2;
        nn[0] = n;
        d = -n0;
        n1 = n2 = n3 = n4 = 0;
    } else {
        if (n0 > n) n0 = n;
        nn[0] = n0;
        if (n1 <= 0) {
            k = 2;
            d = -n1;
            n2 = n3 = n4 = 0;
        } else {
            if (n1 > n) n1 = n;
            nn[1] = n1;
            if (n2 <= 0) {
                k = 3;
                d = -n2;
                n3 = n4 = 0;
            } else {
                if (n2 > n) n2 = n;
                nn[2] = n2;
                if (n3 <= 0) {
                    k = 4;
                    d = -n3;
                    n4 = 0;
                } else {
                    if (n3 > n) n3 = n;
                    nn[3] = n3;
                    if (n4 <= 0) {
                        k = 5;
                        d = -n4;
                    } else {
                        if (n4 > n) n4 = n;
                        nn[4] = n4;
                        d = 4;
                        goto done;
                    }
                }
            }
        }
        if (d == 0) {
            d = k - 2;
            goto done;
        }
        nn[k - 1] = nn[0];
    }

    if (d > MAX_D) panic(bad_specs);
    for (j = 1; k <= d; j++, k++) nn[k] = nn[j];

    done:

    {
        long nverts;
        register long* coef = typed_alloc(n + 1, long, working_storage);
        if (trouble_code) panic(no_room + 1);
        for (k = 0; k <= nn[0]; k++) coef[k] = 1;

        for (j = 1; j <= d; j++) {
            for (k = n, i = n - nn[j] - 1; i >= 0; k--, i--) coef[k] -= coef[i];
            s = 1;
            for (k = 1; k <= n; k++) {
                s += coef[k];
                if (s > 1000000000) panic(very_bad_specs);
                coef[k] = s;
            }
        }

        nverts = coef[n];
        graph_free(working_storage);
        g = new_graph(nverts);
        if (g == NULL)
        panic(no_room);
    }
    sprintf(g->id, "subsets(%lu,%ld,%ld,%ld,%ld,%ld,0x%lx,%d)",
            n, n0, n1, n2, n3, n4, size_bits, directed ? 1 : 0);
    strcpy(g->util_types, "ZZZIIIZZZZZZZZ");

    v = g->vertices;
    yy[d + 1] = 0;
    sig[0] = n;
    for (k = d; k >= 0; k--) yy[k] = yy[k + 1] + nn[k];
    if (yy[0] >= n) {
        k = 0;
        xx[0] = (yy[1] >= n ? 0 : n - yy[1]);
        while (1) {
            for (s = sig[k] - xx[k], k++; k <= d; s -= xx[k], k++) {
                sig[k] = s;
                if (s <= yy[k + 1]) xx[k] = 0;
                else xx[k] = s - yy[k + 1];
            }
            if (s != 0) panic(impossible + 1);

            {
                register char* p = buffer;
                for (k = 0; k <= d; k++) {
                    sprintf(p, ".%ld", xx[k]);
                    while (*p) p++;
                }
                v->name = save_string(&buffer[1]);
                v->x.I = xx[0];
                v->y.I = xx[1];
                v->z.I = xx[2];
            }

            {
                register Vertex* u;
                for (u = g->vertices; u <= v; u++) {
                    register char* p = u->name;
                    long ss = 0;
                    for (j = 0; j <= d; j++, p++) {
                        for (s = (*p++) - '0'; *p >= '0'; p++) s = 10 * s + *p - '0';
                        if (xx[j] < s) ss += xx[j];
                        else ss += s;
                    }
                    if (ss < UL_BITS && (size_bits & (((unsigned long)1) << ss))) {
                        if (directed) new_arc(u, v, 1L);
                        else new_edge(u, v, 1L);
                    }
                }
            }
            v++;
            for (k = d - 1;; k--) {
                if (xx[k] < sig[k] && xx[k] < nn[k]) break;
                if (k == 0) goto last;
            }
            xx[k]++;
        }
    }
    last:
    if (v != g->vertices + g->n)
    panic(impossible);
    if (trouble_code) {
        recycle(g);
        panic(alloc_fault);
    }
    return g;
}

Graph* perms(long n0, long n1, long n2, long n3, long n4, unsigned long max_inv, long directed) {
    Graph* g;
    register long i, j, k;
    register long d;
    register Vertex* v;
    register long s;
    register long n;

    if (n0 == 0) {
        n0 = 1;
        n1 = 0;
    } else if (n0 < 0) {
        n1 = n0;
        n0 = 1;
    }
    n = BUF_SIZE;

    if (n0 == 0) n0 = -2;
    if (n0 < 0) {
        k = 2;
        nn[0] = n;
        d = -n0;
        n1 = n2 = n3 = n4 = 0;
    } else {
        if (n0 > n) n0 = n;
        nn[0] = n0;
        if (n1 <= 0) {
            k = 2;
            d = -n1;
            n2 = n3 = n4 = 0;
        } else {
            if (n1 > n) n1 = n;
            nn[1] = n1;
            if (n2 <= 0) {
                k = 3;
                d = -n2;
                n3 = n4 = 0;
            } else {
                if (n2 > n) n2 = n;
                nn[2] = n2;
                if (n3 <= 0) {
                    k = 4;
                    d = -n3;
                    n4 = 0;
                } else {
                    if (n3 > n) n3 = n;
                    nn[3] = n3;
                    if (n4 <= 0) {
                        k = 5;
                        d = -n4;
                    } else {
                        if (n4 > n) n4 = n;
                        nn[4] = n4;
                        d = 4;
                        goto done;
                    }
                }
            }
        }
        if (d == 0) {
            d = k - 2;
            goto done;
        }
        nn[k - 1] = nn[0];
    }

    if (d > MAX_D) panic(bad_specs);
    for (j = 1; k <= d; j++, k++) nn[k] = nn[j];

    done:

    {
        register long ss;
        for (k = 0, s = ss = 0; k <= d; ss += s * nn[k], s += nn[k], k++)
            if (nn[k] >= BUF_SIZE) panic(bad_specs);

        if (s >= BUF_SIZE) panic(bad_specs + 1);
        n = s;
        if (max_inv == 0 || max_inv > ss) max_inv = ss;
    }

    {
        long nverts;
        register long* coef = typed_alloc(max_inv + 1, long, working_storage);
        if (trouble_code) panic(no_room + 1);
        coef[0] = 1;
        for (j = 1, s = nn[0]; j <= d; s += nn[j], j++)
            for (k = 1; k <= nn[j]; k++) {
                register long ii;
                for (i = max_inv, ii = i - k - s; ii >= 0; ii--, i--) coef[i] -= coef[ii];
                for (i = k, ii = 0; i <= max_inv; i++, ii++) {
                    coef[i] += coef[ii];
                    if (coef[i] > 1000000000) panic(very_bad_specs + 1);
                }
            }
        for (k = 1, nverts = 1; k <= max_inv; k++) {
            nverts += coef[k];
            if (nverts > 1000000000) panic(very_bad_specs);
        }
        graph_free(working_storage);
        g = new_graph(nverts);
        if (g == NULL)
        panic(no_room);
        sprintf(g->id, "perms(%ld,%ld,%ld,%ld,%ld,%lu,%d)",
                n0, n1, n2, n3, n4, max_inv, directed ? 1 : 0);
        strcpy(g->util_types, "VVZZZZZZZZZZZZ");
    }

    {
        register long* xtab, * ytab, * ztab;
        long m = 0;

        xtab = typed_alloc(3 * n + 3, long, working_storage);
        if (trouble_code) {
            recycle(g);
            panic(no_room + 2);
        }
        ytab = xtab + (n + 1);
        ztab = ytab + (n + 1);
        for (j = 0, k = 1, s = nn[0];; k++) {
            xtab[k] = ztab[k] = j;
            if (k == s) {
                if (++j > d) break;
                else s += nn[j];
            }
        }

        v = g->vertices;
        while (1) {
            {
                register char* p;
                register long* q;
                for (p = &buffer[n - 1], q = &xtab[n]; q > xtab; p--, q--) *p = short_imap[*q];
                v->name = save_string(buffer);
                hash_in(v);
            }

            for (j = 1; j < n; j++)
                if (xtab[j] > xtab[j + 1]) {
                    register Vertex* u;
                    buffer[j - 1] = short_imap[xtab[j + 1]];
                    buffer[j] = short_imap[xtab[j]];
                    u = hash_out(buffer);
                    if (u == NULL) panic(impossible + 2);
                    if (directed) new_arc(u, v, 1L);
                    else new_edge(u, v, 1L);
                    buffer[j - 1] = short_imap[xtab[j]];
                    buffer[j] = short_imap[xtab[j + 1]];
                }
            v++;
            for (k = n; k; k--) {
                if (m < max_inv && ytab[k] < k - 1)
                    if (ytab[k] < ytab[k - 1] || ztab[k] > ztab[k - 1]) goto move;
                if (ytab[k]) {
                    for (j = k - ytab[k]; j < k; j++) xtab[j] = xtab[j + 1];
                    m -= ytab[k];
                    ytab[k] = 0;
                    xtab[k] = ztab[k];
                }
            }
            goto last;
            move:
            j = k - ytab[k];
            xtab[j] = xtab[j - 1];
            xtab[j - 1] = ztab[k];
            ytab[k]++;
            m++;
        }
        last:
        if (v != g->vertices + g->n)
        panic(impossible);
        graph_free(working_storage);
    }

    if (trouble_code) {
        recycle(g);
        panic(alloc_fault);
    }
    return g;
}

Graph* parts(unsigned long n, unsigned long max_parts, unsigned long max_size, long directed) {
    Graph* g;
    register long i, j, k;
    register long d;
    register Vertex* v;
    register long s;

    if (max_parts == 0 || max_parts > n) max_parts = n;
    if (max_size == 0 || max_size > n) max_size = n;
    if (max_parts > MAX_D) panic(bad_specs);

    {
        long nverts;
        register long* coef = typed_alloc(n + 1, long, working_storage);
        if (trouble_code) panic(no_room + 1);
        coef[0] = 1;
        for (k = 1; k <= max_parts; k++) {
            for (j = n, i = n - k - max_size; i >= 0; i--, j--) coef[j] -= coef[i];
            for (j = k, i = 0; j <= n; i++, j++) {
                coef[j] += coef[i];
                if (coef[j] > 1000000000) panic(very_bad_specs);
            }
        }
        nverts = coef[n];
        graph_free(working_storage);
        g = new_graph(nverts);
        if (g == NULL)
        panic(no_room);
        sprintf(g->id, "parts(%lu,%lu,%lu,%d)",
                n, max_parts, max_size, directed ? 1 : 0);
        strcpy(g->util_types, "VVZZZZZZZZZZZZ");
    }

    v = g->vertices;
    xx[0] = max_size;
    sig[1] = n;
    for (k = max_parts, s = 1; k > 0; k--, s++) yy[k] = s;
    if (max_size * max_parts >= n) {
        k = 1;
        xx[1] = (n - 1) / max_parts + 1;
        while (1) {
            for (s = sig[k] - xx[k], k++; s; k++) {
                sig[k] = s;
                xx[k] = (s - 1) / yy[k] + 1;
                s -= xx[k];
            }
            d = k - 1;

            {
                register char* p = buffer;
                for (k = 1; k <= d; k++) {
                    sprintf(p, "+%ld", xx[k]);
                    while (*p) p++;
                }
                v->name = save_string(&buffer[1]);
                hash_in(v);
            }

            if (d < max_parts) {
                xx[d + 1] = 0;
                for (j = 1; j <= d; j++) {
                    if (xx[j] != xx[j + 1]) {
                        long a, b;
                        for (b = xx[j] / 2, a = xx[j] - b; b; a++, b--) {
                            register Vertex* u;
                            register char* p = buffer;
                            for (k = j + 1; xx[k] > a; k++) nn[k - 1] = xx[k];
                            nn[k - 1] = a;
                            for (; xx[k] > b; k++) nn[k] = xx[k];
                            nn[k] = b;
                            for (; k <= d; k++) nn[k + 1] = xx[k];
                            for (k = 1; k <= d + 1; k++) {
                                sprintf(p, "+%ld", nn[k]);
                                while (*p) p++;
                            }
                            u = hash_out(&buffer[1]);
                            if (u == NULL) panic(impossible + 2);
                            if (directed) new_arc(v, u, 1L);
                            else new_edge(v, u, 1L);
                        }
                    }
                    nn[j] = xx[j];
                }
            }
            v++;
            for (k = d - 1;; k--) {
                if (xx[k] < sig[k] && xx[k] < xx[k - 1]) break;
                if (k == 1) goto last;
            }
            xx[k]++;
        }
    }
    last:
    if (v != g->vertices + g->n)
    panic(impossible);
    if (trouble_code) {
        recycle(g);
        panic(alloc_fault);
    }
    return g;
}

Graph* binary(unsigned long n, unsigned long max_height, long directed) {
    Graph* g;
    register long i, j, k;
    register long d;
    register Vertex* v;
    register long s;

    if (2 * n + 2 > BUF_SIZE) panic(bad_specs);
    if (max_height == 0 || max_height > n) max_height = n;
    if (max_height > 30) panic(very_bad_specs);

    {
        long nverts;
        if (n >= 20 && max_height >= 6) {
            register float ss;
            d = (1L << max_height) - 1 - n;
            if (d > 8) panic(bad_specs + 1);
            if (d < 0) nverts = 0;
            else {
                nn[0] = nn[1] = 1;
                for (k = 2; k <= d; k++) nn[k] = 0;
                for (j = 2; j <= max_height; j++) {
                    for (k = d; k; k--) {
                        for (ss = 0.0, i = k; i >= 0; i--) ss += ((float)nn[i]) * ((float)nn[k - i]);
                        if (ss > MAX_NNN) panic(very_bad_specs + 1);
                        for (s = 0, i = k; i >= 0; i--) s += nn[i] * nn[k - i];
                        nn[k] = s;
                    }
                    i = (1L << j) - 1;
                    if (i <= d) nn[i]++;
                }
                nverts = nn[d];
            }
        } else {
            nn[0] = nn[1] = 1;
            for (k = 2; k <= n; k++) nn[k] = 0;
            for (j = 2; j <= max_height; j++)
                for (k = n - 1; k; k--) {
                    for (s = 0, i = k; i >= 0; i--) s += nn[i] * nn[k - i];
                    nn[k + 1] = s;
                }
            nverts = nn[n];
        }
        g = new_graph(nverts);
        if (g == NULL)
        panic(no_room);
        sprintf(g->id, "binary(%lu,%lu,%d)",
                n, max_height, directed ? 1 : 0);
        strcpy(g->util_types, "VVZZZZZZZZZZZZ");
    }

    {
        register long* xtab, * ytab, * ltab, * stab;

        xtab = typed_alloc(8 * n + 4, long, working_storage);
        if (trouble_code) {
            recycle(g);
            panic(no_room + 2);
        }
        d = n + n;
        ytab = xtab + (d + 1);
        ltab = ytab + (d + 1);
        stab = ltab + (d + 1);
        ltab[0] = 1L << max_height;
        stab[0] = n;

        v = g->vertices;
        if (ltab[0] > n) {
            k = 0;
            xtab[0] = n ? 1 : 0;
            while (1) {
                for (j = k + 1; j <= d; j++) {
                    if (xtab[j - 1]) {
                        ltab[j] = ltab[j - 1] >> 1;
                        ytab[j] = ytab[j - 1] + ltab[j];
                        stab[j] = stab[j - 1];
                    } else {
                        ytab[j] = ytab[j - 1] & (ytab[j - 1] - 1);
                        ltab[j] = ytab[j - 1] - ytab[j];
                        stab[j] = stab[j - 1] - 1;
                    }
                    if (stab[j] <= ytab[j]) xtab[j] = 0;
                    else xtab[j] = 1;
                }

                {
                    register char* p = buffer;
                    for (k = 0; k <= d; k++, p++) *p = (xtab[k] ? '.' : 'x');
                    v->name = save_string(buffer);
                    hash_in(v);
                }

                for (j = 0; j < d; j++)
                    if (xtab[j] == 1 && xtab[j + 1] == 1) {
                        for (i = j + 1, s = 0; s >= 0; s += (xtab[i + 1] << 1) - 1, i++) xtab[i] = xtab[i + 1];
                        xtab[i] = 1;
                        {
                            register char* p = buffer;
                            register Vertex* u;
                            for (k = 0; k <= d; k++, p++) *p = (xtab[k] ? '.' : 'x');
                            u = hash_out(buffer);
                            if (u) {
                                if (directed) new_arc(v, u, 1L);
                                else new_edge(v, u, 1L);
                            }
                        }
                        for (i--; i > j; i--) xtab[i + 1] = xtab[i];
                        xtab[i + 1] = 1;
                    }
                v++;
                for (k = d - 1;; k--) {
                    if (k <= 0) goto last;
                    if (xtab[k]) break;
                }
                for (k--;; k--) {
                    if (xtab[k] == 0 && ltab[k] > 1) break;
                    if (k == 0) goto last;
                }
                xtab[k]++;
            }
        }
    }
    last:
    if (v != g->vertices + g->n)
    panic(impossible);
    graph_free(working_storage);

    if (trouble_code) {
        recycle(g);
        panic(alloc_fault);
    }
    return g;
}

Graph* complement(Graph* g, long copy, long self, long directed) {
    Graph* new_g;
    register long i, j, k;
    register long d;
    register Vertex* v;
    register long s;
    register long n;
    register Vertex* u;
    register siz_t delta;

    if (g == NULL) panic(missing_operand);

    n = g->n;
    new_g = new_graph(n);
    if (new_g == NULL)
    panic(no_room);
    delta = ((siz_t)(new_g->vertices)) - ((siz_t)(g->vertices));
    for (u = new_g->vertices, v = g->vertices; v < g->vertices + n; u++, v++)
        u->name = save_string(v->name);

    sprintf(buffer, ",%d,%d,%d)", copy ? 1 : 0, self ? 1 : 0, directed ? 1 : 0);
    make_compound_id(new_g, "complement(", g, buffer);

    for (v = g->vertices; v < g->vertices + n; v++) {
        u = vert_offset(v, delta);

        {
            register Arc* a;
            for (a = v->arcs; a; a = a->next) vert_offset(a->tip, delta)->tmp = u;
        }
        if (directed) {
            for (Vertex* vv = new_g->vertices; vv < new_g->vertices + n; vv++)
                if ((vv->tmp == u && copy) || (vv->tmp != u && !copy))
                    if (vv != u || self) new_arc(u, vv, 1L);
        } else {
            for (Vertex* vv = (self ? u : u + 1); vv < new_g->vertices + n; vv++)
                if ((vv->tmp == u && copy) || (vv->tmp != u && !copy))
                    new_edge(u, vv, 1L);
        }
    }
    for (v = new_g->vertices; v < new_g->vertices + n; v++) v->tmp = NULL;

    if (trouble_code) {
        recycle(new_g);
        panic(alloc_fault);
    }
    return new_g;
}

Graph* gunion(Graph* g, Graph* gg, long multi, long directed) {
    Graph* new_g;
    register long i, j, k;
    register long d;
    register Vertex* v;
    register long s;
    register long n;
    register Vertex* u;
    register siz_t delta, ddelta;

    if (g == NULL || gg == NULL) panic(missing_operand);

    n = g->n;
    new_g = new_graph(n);
    if (new_g == NULL)
    panic(no_room);
    delta = ((siz_t)(new_g->vertices)) - ((siz_t)(g->vertices));
    for (u = new_g->vertices, v = g->vertices; v < g->vertices + n; u++, v++)
        u->name = save_string(v->name);

    sprintf(buffer, ",%d,%d)", multi ? 1 : 0, directed ? 1 : 0);
    make_double_compound_id(new_g, "gunion(", g, ",", gg, buffer);
    ddelta = ((siz_t)(new_g->vertices)) - ((siz_t)(gg->vertices));

    for (v = g->vertices; v < g->vertices + n; v++) {
        register Arc* a;
        register Vertex* vv = vert_offset(v, delta);

        for (a = v->arcs; a; a = a->next) {
            u = vert_offset(a->tip, delta);
            if (directed) {
                if (multi || u->tmp != vv) new_arc(vv, u, a->len);
                else {
                    register Arc* b = u->tlen;
                    if (a->len < b->len) b->len = a->len;
                }
                u->tmp = vv;
                u->tlen = vv->arcs;
            } else if (u >= vv) {
                if (multi || u->tmp != vv) new_edge(vv, u, a->len);
                else {
                    register Arc* b = u->tlen;
                    if (a->len < b->len) b->len = (b + 1)->len = a->len;
                }
                u->tmp = vv;
                u->tlen = vv->arcs;
                if (u == vv && a->next == a + 1) a++;
            }
        }
        if (vv < gg->vertices + gg->n) for (a = vv->arcs; a; a = a->next) {
                u = vert_offset(a->tip, ddelta);
                if (u < new_g->vertices + n)
                    if (directed) {
                        if (multi || u->tmp != vv) new_arc(vv, u, a->len);
                        else {
                            register Arc* b = u->tlen;
                            if (a->len < b->len) b->len = a->len;
                        }
                        u->tmp = vv;
                        u->tlen = vv->arcs;
                    } else if (u >= vv) {
                        if (multi || u->tmp != vv) new_edge(vv, u, a->len);
                        else {
                            register Arc* b = u->tlen;
                            if (a->len < b->len) b->len = (b + 1)->len = a->len;
                        }
                        u->tmp = vv;
                        u->tlen = vv->arcs;
                        if (u == vv && a->next == a + 1) a++;
                    }
            }
    }
    for (v = new_g->vertices; v < new_g->vertices + n; v++)
        v->tmp = NULL, v->tlen = NULL;

    if (trouble_code) {
        recycle(new_g);
        panic(alloc_fault);
    }
    return new_g;
}

Graph* intersection(Graph* g, Graph* gg, long multi, long directed) {
    Graph* new_g;
    register long i, j, k;
    register long d;
    register Vertex* v;
    register long s;
    register long n;
    register Vertex* u;
    register siz_t delta, ddelta;

    if (g == NULL || gg == NULL) panic(missing_operand);

    n = g->n;
    new_g = new_graph(n);
    if (new_g == NULL)
    panic(no_room);
    delta = ((siz_t)(new_g->vertices)) - ((siz_t)(g->vertices));
    for (u = new_g->vertices, v = g->vertices; v < g->vertices + n; u++, v++)
        u->name = save_string(v->name);

    sprintf(buffer, ",%d,%d)", multi ? 1 : 0, directed ? 1 : 0);
    make_double_compound_id(new_g, "intersection(", g, ",", gg, buffer);
    ddelta = ((siz_t)(new_g->vertices)) - ((siz_t)(gg->vertices));

    for (v = g->vertices; v < g->vertices + n; v++) {
        register Arc* a;
        register Vertex* vv = vert_offset(v, delta);

        register Vertex* vvv = vert_offset(vv, -ddelta);

        if (vvv >= gg->vertices + gg->n) continue;

        for (a = v->arcs; a; a = a->next) {
            u = vert_offset(a->tip, delta);
            if (u->tmp == vv) {
                u->mult++;
                if (a->len < u->minlen) u->minlen = a->len;
            } else {
                u->tmp = vv;
                u->mult = 0;
                u->minlen = a->len;
            }
            if (u == vv && !directed && a->next == a + 1) a++;
        }

        for (a = vvv->arcs; a; a = a->next) {
            u = vert_offset(a->tip, ddelta);
            if (u >= new_g->vertices + n) continue;
            if (u->tmp == vv) {
                long l = u->minlen;
                if (a->len > l) l = a->len;
                if (u->mult < 0) {
                    register Arc* b = u->tlen;
                    if (l < b->len) {
                        b->len = l;
                        if (!directed) (b + 1)->len = l;
                    }
                } else {
                    if (directed) new_arc(vv, u, l);
                    else {
                        if (vv <= u) new_edge(vv, u, l);
                        if (vv == u && a->next == a + 1) a++;
                    }
                    if (!multi) {
                        u->tlen = vv->arcs;
                        u->mult = -1;
                    } else if (u->mult == 0) u->tmp = NULL;
                    else u->mult--;
                }
            }
        }
    }

    for (v = new_g->vertices; v < new_g->vertices + n; v++) {
        v->tmp = NULL;
        v->tlen = NULL;
        v->mult = 0;
        v->minlen = 0;
    }

    if (trouble_code) {
        recycle(new_g);
        panic(alloc_fault);
    }
    return new_g;
}

Graph* lines(Graph* g, long directed) {
    Graph* new_g;
    register long i, j, k;
    register long d;
    register Vertex* v;
    register long s;
    register long m;
    register Vertex* u;

    if (g == NULL) panic(missing_operand);

    m = (directed ? g->m : (g->m) / 2);
    new_g = new_graph(m);
    if (new_g == NULL)
    panic(no_room);
    make_compound_id(new_g, "lines(", g, directed ? ",1)" : ",0)");
    u = new_g->vertices;
    for (v = g->vertices + g->n - 1; v >= g->vertices; v--) {
        register Arc* a;
        register long mapped = 0;
        for (a = v->arcs; a; a = a->next) {
            register Vertex* vv = a->tip;
            if (!directed) {
                if (vv < v) continue;
                if (vv >= g->vertices + g->n) goto near_panic;
            }

            u->u.V = v;
            u->v.V = vv;
            u->w.A = a;
            if (!directed) {
                if (u >= new_g->vertices + m || (a + 1)->tip != v) goto near_panic;
                if (v == vv && a->next == a + 1) a++;
                else (a + 1)->tip = u;
            }
            sprintf(buffer, "%.*s-%c%.*s", (BUF_SIZE - 3) / 2, v->name,
                    directed ? '>' : '-', BUF_SIZE / 2 - 1, vv->name);
            u->name = save_string(buffer);

            if (!mapped) {
                u->map = v->map;
                v->map = u;
                mapped = 1;
            }
            u++;
        }
    }
    if (u != new_g->vertices + m) goto near_panic;

    if (directed) {
        for (u = new_g->vertices; u < new_g->vertices + m; u++) {
            v = u->v.V;
            if (v->arcs) {
                v = v->map;
                do {
                    new_arc(u, v, 1L);
                    v++;
                } while (v->u.V == u->v.V);
            }
        }
    } else {
        for (u = new_g->vertices; u < new_g->vertices + m; u++) {
            register Vertex* vv;
            register Arc* a;
            register long mapped = 0;
            v = u->u.V;
            for (vv = v->map; vv < u; vv++) new_edge(u, vv, 1L);
            v = u->v.V;
            for (a = v->arcs; a; a = a->next) {
                vv = a->tip;
                if (vv < u && vv >= new_g->vertices) new_edge(u, vv, 1L);
                else if (vv >= v && vv < g->vertices + g->n) mapped = 1;
            }
            if (mapped && v > u->u.V)
                for (vv = v->map; vv->u.V == v; vv++) new_edge(u, vv, 1L);
        }
    }

    for (u = new_g->vertices, v = NULL; u < new_g->vertices + m; u++) {
        if (u->u.V != v) {
            v = u->u.V;
            v->map = u->map;
            u->map = NULL;
        }
        if (!directed) ((u->w.A) + 1)->tip = v;
    }

    if (trouble_code) {
        recycle(new_g);
        panic(alloc_fault);
    }
    return new_g;
    near_panic:
    m = u - new_g->vertices;
    for (u = new_g->vertices, v = NULL; u < new_g->vertices + m; u++) {
        if (u->u.V != v) {
            v = u->u.V;
            v->map = u->map;
            u->map = NULL;
        }
        if (!directed) ((u->w.A) + 1)->tip = v;
    }
    recycle(new_g);
    panic(invalid_operand);
}

Graph* product(Graph* g, Graph* gg, long type, long directed) {
    Graph* new_g;
    register long i, j, k, d, s, n;
    register Vertex* v;
    register Vertex* u;
    register Vertex* vv;

    if (g == NULL || gg == NULL) panic(missing_operand);

    float test_product = ((float)(g->n)) * ((float)(gg->n));
    if (test_product > MAX_NNN) panic(very_bad_specs);
    n = (g->n) * (gg->n);
    new_g = new_graph(n);
    if (new_g == NULL) panic(no_room);
    for (u = new_g->vertices, v = g->vertices, vv = gg->vertices; u < new_g->vertices + n; u++) {
        sprintf(buffer, "%.*s,%.*s", BUF_SIZE / 2 - 1, v->name, (BUF_SIZE - 1) / 2, vv->name);
        u->name = save_string(buffer);
        if (++vv == gg->vertices + gg->n) vv = gg->vertices, v++;
    }
    sprintf(buffer, ",%d,%d)", (type ? 2 : 0) - (int)(type & 1), directed ? 1 : 0);
    make_double_compound_id(new_g, "product(", g, ",", gg, buffer);

    if ((type & 1) == 0) {
        register Vertex* uu, * uuu;
        register Arc* a;
        register siz_t delta;
        delta = ((siz_t)(new_g->vertices)) - ((siz_t)(gg->vertices));
        for (u = gg->vertices; u < gg->vertices + gg->n; u++)
            for (a = u->arcs; a; a = a->next) {
                v = a->tip;
                if (!directed) {
                    if (u > v) continue;
                    if (u == v && a->next == a + 1) a++;
                }
                for (uu = vert_offset(u, delta), vv = vert_offset(v, delta); uu < new_g->vertices + n; uu += gg->n, vv += gg->n)
                    if (directed) new_arc(uu, vv, a->len);
                    else new_edge(uu, vv, a->len);
            }

        for (u = g->vertices, uu = new_g->vertices; uu < new_g->vertices + n; u++, uu += gg->n)
            for (a = u->arcs; a; a = a->next) {
                v = a->tip;
                if (!directed) {
                    if (u > v) continue;
                    if (u == v && a->next == a + 1) a++;
                }
                vv = new_g->vertices + ((gg->n) * (v - g->vertices));
                for (uuu = uu; uuu < uu + gg->n; uuu++, vv++)
                    if (directed) new_arc(uuu, vv, a->len);
                    else new_edge(uuu, vv, a->len);
            }
    }

    if (type) {
        Vertex* uu;
        Arc* a;
        siz_t delta0 = ((siz_t)(new_g->vertices)) - ((siz_t)(gg->vertices));
        siz_t del = (gg->n) * sizeof(Vertex);
        register siz_t delta, ddelta;
        for (uu = g->vertices, delta = delta0; uu < g->vertices + g->n; uu++, delta += del)
            for (a = uu->arcs; a; a = a->next) {
                vv = a->tip;
                if (!directed) {
                    if (uu > vv) continue;
                    if (uu == vv && a->next == a + 1) a++;
                }
                ddelta = delta0 + del * (vv - g->vertices);
                for (u = gg->vertices; u < gg->vertices + gg->n; u++) {
                    register Arc* aa;
                    for (aa = u->arcs; aa; aa = aa->next) {
                        long length = a->len;
                        if (length > aa->len) length = aa->len;
                        v = aa->tip;
                        if (directed)
                            new_arc(vert_offset(u, delta), vert_offset(v, ddelta), length);
                        else new_edge(vert_offset(u, delta), vert_offset(v, ddelta), length);
                    }
                }
            }
    }

    if (trouble_code) {
        recycle(new_g);
        panic(alloc_fault);
    }
    return new_g;
}

Graph* induced(Graph* g, char* description, long self, long multi, long directed) {
    Graph* new_g;
    register long i, j, k, d, s, n, nn;
    register Vertex* v;
    register Vertex* u;

    if (g == NULL) panic(missing_operand);

    for (v = g->vertices; v < g->vertices + g->n; v++)
        if (v->ind > 0) {
            if (n > IND_GRAPH) panic(very_bad_specs);
            if (v->ind >= IND_GRAPH) {
                if (v->subst == NULL) panic(missing_operand + 1);
                n += v->subst->n;
            } else n += v->ind;
        } else if (v->ind < -nn) nn = -(v->ind);
    if (n > IND_GRAPH || nn > IND_GRAPH) panic(very_bad_specs + 1);
    n += nn;

    new_g = new_graph(n);
    if (new_g == NULL) panic(no_room);

    for (k = 1, u = new_g->vertices; k <= nn; k++, u++) {
        u->mult = -k;
        sprintf(buffer, "%ld", -k);
        u->name = save_string(buffer);
    }
    for (v = g->vertices; v < g->vertices + g->n; v++)
        if ((k = v->ind) < 0) v->map = (new_g->vertices) - (k + 1);
        else if (k > 0) {
            u->mult = k;
            v->map = u;
            if (k <= 2) {
                u->name = save_string(v->name);
                u++;
                if (k == 2) {
                    sprintf(buffer, "%s'", v->name);
                    u->name = save_string(buffer);
                    u++;
                }
            } else if (k >= IND_GRAPH) {
                register Graph* gg = v->subst;
                register Vertex* vv = gg->vertices;
                register Arc* a;
                siz_t delta = ((siz_t)u) - ((siz_t)vv);
                for (j = 0; j < v->subst->n; j++, u++, vv++) {
                    sprintf(buffer, "%.*s:%.*s", BUF_SIZE / 2 - 1, v->name, (BUF_SIZE - 1) / 2, vv->name);
                    u->name = save_string(buffer);
                    for (a = vv->arcs; a; a = a->next) {
                        register Vertex* vvv = a->tip;
                        Vertex* uu = vert_offset(vvv, delta);
                        if (vvv == vv && !self) continue;
                        if (uu->tmp == u && !multi) {
                            register Arc* b = uu->tlen;
                            if (a->len < b->len) {
                                b->len = a->len;
                                if (!directed)(b + 1)->len = a->len;
                            }
                            continue;
                        }
                        if (!directed) {
                            if (vvv < vv) continue;
                            if (vvv == vv && a->next == a + 1) a++;
                            new_edge(u, uu, a->len);
                        } else new_arc(u, uu, a->len);
                        uu->tmp = u;
                        uu->tlen = ((directed || u <= uu) ? u->arcs : uu->arcs);
                    }
                }
            } else for (j = 0; j < k; j++, u++) {
                    sprintf(buffer, "%.*s:%ld", BUF_SIZE - 12, v->name, j);
                    u->name = save_string(buffer);
                }
        }

    sprintf(buffer, ",%s,%d,%d,%d)", description ? description : null_string, self ? 1 : 0, multi ? 1 : 0, directed ? 1 : 0);
    make_compound_id(new_g, "induced(", g, buffer);

    for (v = g->vertices; v < g->vertices + g->n; v++) {
        u = v->map;
        if (u) {
            register Arc* a;
            register Vertex* uu, * vv;
            k = u->mult;
            if (k < 0) k = 1;
            else if (k >= IND_GRAPH) k = v->subst->n;
            for (; k; k--, u++) {
                if (!multi)
                    for (a = u->arcs; a; a = a->next) {
                        a->tip->tmp = u;
                        if (directed || a->tip > u || a->next == a + 1) a->tip->tlen = a;
                        else a->tip->tlen = a + 1;
                    }
                for (a = v->arcs; a; a = a->next) {
                    vv = a->tip;
                    uu = vv->map;
                    if (uu == NULL) continue;
                    j = uu->mult;
                    if (j < 0) j = 1;
                    else if (j >= IND_GRAPH) j = vv->subst->n;
                    if (!directed) {
                        if (vv < v) continue;
                        if (vv == v) {
                            if (a->next == a + 1) a++;
                            j = k, uu = u;
                        }
                    }
                    for (; j; j--, uu++) {
                        if (u == uu && !self) continue;
                        if (uu->tmp == u && !multi) {
                            register Arc* b = uu->tlen;
                            if (a->len < b->len) {
                                b->len = a->len;
                                if (!directed)(b + 1)->len = a->len;
                            }
                            continue;
                        }
                        if (directed) new_arc(u, uu, a->len);
                        else new_edge(u, uu, a->len);
                        uu->tmp = u;
                        uu->tlen = ((directed || u <= uu) ? u->arcs : uu->arcs);
                    }
                }
            }
        }
    }

    for (v = g->vertices; v < g->vertices + g->n; v++)
        if (v->map) v->ind = v->map->mult;
    for (v = new_g->vertices; v < new_g->vertices + n; v++)
        v->u.I = v->v.I = v->z.I = 0;

    if (trouble_code) {
        recycle(new_g);
        panic(alloc_fault);
    }
    return new_graph;
}

Graph* bi_complete(unsigned long n1, unsigned long n2, long directed) {
    Graph* g = board(2L, 0L, 0L, 0L, 1L, 0L, directed);
    if (g) {
        g->vertices->ind = n1;
        (g->vertices + 1)->ind = n2;
        g = induced(g, NULL, 0L, 0L, directed);
        if (g) {
            sprintf(g->id, "bi_complete(%lu,%lu,%d)", n1, n2, directed ? 1 : 0);
            mark_bipartite(g, n1);
        }
    }
    return g;
}

Graph* wheel(unsigned long n, unsigned long n1, long directed) {
    Graph* g = board(2L, 0L, 0L, 0L, 1L, 0L, directed);

    if (g) {
        g->vertices->ind = n1;
        (g->vertices + 1)->ind = IND_GRAPH;
        (g->vertices + 1)->subst = board(n, 0L, 0L, 0L, 1L, 1L, directed);

        g = induced(g, NULL, 0L, 0L, directed);
        if (g) {
            sprintf(g->id, "wheel(%lu,%lu,%d)", n, n1, directed ? 1 : 0);
        }
    }
    return g;
}
