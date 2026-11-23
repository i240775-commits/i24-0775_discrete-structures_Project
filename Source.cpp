#define _CRT_SECURE_NO_WARNINGS   // Visual Studio: allow fopen/fscanf/sprintf (must be the first line)
// CS1005 Discrete Structures - Fall 2026 - Network motifs (minimal libraries)
// R = 2877  (775 + 1050 + 1052)  ->  pair 27 = (1,10), 45 random networks, 8*edges attempts, seed 2877
// Only uses <stdio.h>, <stdlib.h>, <math.h>. No STL: plain arrays and our own random generator.
// Build: g++ -O2 motifs_min.cpp -o motifs
// Run:   ./motifs          (reads email-Eu-core.txt and email-Eu-core-department-labels.txt
//                           from the current folder; then: python3 plot.py)
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAXN 1010   // node ids are 0..1004
#define MAXE 26000  // at most 25571 edges
#define MAXG 60     // at most 20+30 = 50 random networks

char full[MAXN][MAXN];  // cleaned full graph: full[u][v]=1 if u->v
int dept[MAXN];         // department of each node
char g[MAXN][MAXN];     // directed graph currently being processed
char und[MAXN][MAXN];   // undirected version (for neighbours)
int nb[MAXN][MAXN];     // neighbour lists
int deg[MAXN];

// ---------------------------------------------------------------- random numbers
// splitmix64: a small generator written by us, seeded with R
unsigned long long state;
unsigned long long rnd() {
    state += 0x9E3779B97F4A7C15ULL;
    unsigned long long z = state;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

// ---------------------------------------------------------------- T3
// bit k is 1 if edge k is present: 0:0->1  1:0->2  2:1->0  3:1->2  4:2->0  5:2->1
int EP[6][2] = { {0, 1}, {0, 2}, {1, 0}, {1, 2}, {2, 0}, {2, 1} };
int P[6][3] = { {0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0} };
int CODES[13] = { 3, 6, 7, 10, 11, 15, 21, 23, 25, 27, 30, 31, 63 };

// returns class 1..13 (0 if not connected). Uses the matrix g.
int triadClass(int a, int b, int c) {
    int nodes[3] = { a, b, c };
    int best = 64;
    for (int p = 0; p < 6; p++) {              // all 6 orderings
        int code = 0;
        for (int k = 0; k < 6; k++)
            if (g[nodes[P[p][EP[k][0]]]][nodes[P[p][EP[k][1]]]]) code |= (1 << k);
        if (code < best) best = code;          // canonical code = smallest
    }
    for (int k = 0; k < 13; k++)
        if (CODES[k] == best) return k + 1;
    return 0;
}

// ---------------------------------------------------------------- T4
// Count each connected triad exactly once.
// Triad {v,a,b} is seen from every node adjacent to the other two.
// If a,b are not adjacent, only v is such a node -> count it.
// If a,b are adjacent (triangle in the undirected sense), all three are -> count only from the smallest.
void census(int* eu, int* ev, int m, long long* cnt) {
    for (int i = 0; i < m; i++) g[eu[i]][ev[i]] = 1;
    for (int i = 0; i < m; i++) {
        int u = eu[i], v = ev[i];
        if (!und[u][v]) {
            und[u][v] = und[v][u] = 1;
            nb[u][deg[u]++] = v;
            nb[v][deg[v]++] = u;
        }
    }
    for (int k = 0; k < 13; k++) cnt[k] = 0;
    for (int v = 0; v < MAXN; v++)
        for (int i = 0; i < deg[v]; i++)
            for (int j = i + 1; j < deg[v]; j++) {
                int a = nb[v][i], b = nb[v][j];
                if (!und[a][b] || (v < a && v < b)) cnt[triadClass(v, a, b) - 1]++;
            }
    for (int i = 0; i < m; i++) {              // clean up for next call
        int u = eu[i], v = ev[i];
        g[u][v] = 0;
        und[u][v] = und[v][u] = 0;
        deg[u] = deg[v] = 0;
    }
}

// ---------------------------------------------------------------- T5
// double-edge swaps; result written to ru, rv
void randomize(int* eu, int* ev, int m, long long attempts, int* ru, int* rv) {
    for (int i = 0; i < m; i++) { ru[i] = eu[i]; rv[i] = ev[i]; g[ru[i]][rv[i]] = 1; }
    for (long long t = 0; t < attempts; t++) {
        int i = rnd() % m, j = rnd() % m;
        int a = ru[i], b = rv[i], c = ru[j], d = rv[j];
        if (a == d || c == b || g[a][d] || g[c][b]) continue;   // reject
        g[a][b] = 0; g[c][d] = 0;
        g[a][d] = 1; g[c][b] = 1;
        rv[i] = d;           // edge i becomes (a,d)
        ru[j] = c; rv[j] = b; // edge j becomes (c,b)
    }
    for (int i = 0; i < m; i++) g[ru[i]][rv[i]] = 0;
}

void degrees(int* eu, int* ev, int m, int* out, int* in) {
    for (int i = 0; i < MAXN; i++) out[i] = in[i] = 0;
    for (int i = 0; i < m; i++) { out[eu[i]]++; in[ev[i]]++; }
}

// ---------------------------------------------------------------- T6
int ru[MAXE], rv[MAXE];
long long samples[MAXG][13];

void getZ(int* eu, int* ev, int m, int nGraphs, long long attempts, long long* real, double* Z) {
    static int o0[MAXN], i0[MAXN], o1[MAXN], i1[MAXN];
    census(eu, ev, m, real);
    degrees(eu, ev, m, o0, i0);
    for (int s = 0; s < nGraphs; s++) {
        randomize(eu, ev, m, attempts, ru, rv);
        degrees(ru, rv, m, o1, i1);
        for (int x = 0; x < MAXN; x++)
            if (o0[x] != o1[x] || i0[x] != i1[x]) { printf("ERROR: degrees changed\n"); exit(1); }
        census(ru, rv, m, samples[s]);
    }
    for (int k = 0; k < 13; k++) {
        double mean = 0, var = 0;
        for (int s = 0; s < nGraphs; s++) mean += samples[s][k];
        mean /= nGraphs;
        for (int s = 0; s < nGraphs; s++) var += (samples[s][k] - mean) * (samples[s][k] - mean);
        double sd = sqrt(var / nGraphs);       // population std
        Z[k] = (sd == 0) ? 0 : (real[k] - mean) / sd;
    }
}

// ---------------------------------------------------------------- main (T1, T2)
int eu[2][MAXE], ev[2][MAXE];

int main(int argc, char** argv) {
    long long R = 2877;                          // our group's R (optional override: motifs <R> [data_dir])
    if (argc > 1) R = atoll(argv[1]);
    const char* dir = argc > 2 ? argv[2] : ".";  // data files are in the same folder as the program
    char path[300];

    // load edges; skip self-loops; the matrix removes duplicates
    sprintf(path, "%s/email-Eu-core.txt", dir);
    FILE* f = fopen(path, "r");
    if (!f) { printf("cannot open %s\n", path); return 1; }
    int u, v;
    while (fscanf(f, "%d %d", &u, &v) == 2) if (u != v) full[u][v] = 1;
    fclose(f);

    for (int i = 0; i < MAXN; i++) dept[i] = -1;
    sprintf(path, "%s/email-Eu-core-department-labels.txt", dir);
    f = fopen(path, "r");
    if (!f) { printf("cannot open %s\n", path); return 1; }
    while (fscanf(f, "%d %d", &u, &v) == 2) dept[u] = v;
    fclose(f);

    // internal edges per department, then top 20 (ties -> smaller id)
    int internal[100] = { 0 }, used[100] = { 0 }, top[20];
    long long total = 0;
    for (u = 0; u < MAXN; u++)
        for (v = 0; v < MAXN; v++)
            if (full[u][v]) { total++; if (dept[u] >= 0 && dept[u] == dept[v]) internal[dept[u]]++; }
    for (int r = 0; r < 20; r++) {
        int best = -1;
        for (int d = 0; d < 100; d++)
            if (!used[d] && (best == -1 || internal[d] > internal[best])) best = d;
        used[best] = 1;
        top[r] = best;
    }

    // pair number R % 190, lexicographic
    int idx = (int)(R % 190), pi = 0, pj = 0, c = 0;
    for (int i = 0; i < 20; i++)
        for (int j = i + 1; j < 20; j++) { if (c == idx) { pi = i; pj = j; } c++; }
    int D[2] = { top[pi], top[pj] };
    int nGraphs = 20 + (int)(R % 31);
    long long mult = 5 + R % 6;

    printf("Cleaned graph: %lld edges\n", total);
    printf("R=%lld pair#=%d -> top20[%d]=%d and top20[%d]=%d\n", R, idx, pi, D[0], pj, D[1]);
    printf("networks=%d attempts=%lld*edges seed=%lld\n", nGraphs, mult, R);

    state = (unsigned long long)R;
    long long cnt[2][13];
    double Z[2][13];
    for (int s = 0; s < 2; s++) {
        int m = 0;                              // induced subgraph of department D[s]
        for (u = 0; u < MAXN; u++)
            for (v = 0; v < MAXN; v++)
                if (full[u][v] && dept[u] == D[s] && dept[v] == D[s]) { eu[s][m] = u; ev[s][m] = v; m++; }
        printf("Dept %c (%d): %d edges\n", 'A' + s, D[s], m);
        getZ(eu[s], ev[s], m, nGraphs, mult * m, cnt[s], Z[s]);
    }

    f = fopen("results.csv", "w");
    fprintf(f, "R,%lld\ndeptA,%d\ndeptB,%d\nclass,count_A,Z_A,count_B,Z_B\n", R, D[0], D[1]);
    printf("\nClass | N_A    Z_A | N_B    Z_B\n");
    for (int k = 0; k < 13; k++) {
        printf("%5d | %lld  %.2f | %lld  %.2f\n", k + 1, cnt[0][k], Z[0][k], cnt[1][k], Z[1][k]);
        fprintf(f, "%d,%lld,%.4f,%lld,%.4f\n", k + 1, cnt[0][k], Z[0][k], cnt[1][k], Z[1][k]);
    }
    fclose(f);
    return 0;
}