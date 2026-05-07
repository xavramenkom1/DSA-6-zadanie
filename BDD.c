#include "bdd.h"
#include <ctype.h>
#include <math.h>
#include <time.h>

/* ═══════════════════════════════════════════════════════════════════════════
   SECTION 1 – DNF PARSER & EVALUATOR
   ═══════════════════════════════════════════════════════════════════════════ */

/*
 * Parse the set of variable names that appear in a DNF expression.
 * Fills vars[] with unique variable indices (A=0,B=1,...,Z=25).
 * Returns number of unique variables found.
 */
static int parse_vars(const char *expr, int *vars_out) {
    int seen[MAX_VARS] = {0};
    int count = 0;
    for (int i = 0; expr[i]; i++) {
        if (isupper(expr[i])) {
            int idx = expr[i] - 'A';
            if (!seen[idx]) {
                seen[idx] = 1;
                vars_out[count++] = idx;
            }
        }
    }
    return count;
}

/*
 * Evaluate one minterm (product term) of a DNF.
 * Returns 1 if the term is satisfied, 0 otherwise.
 * term_end is set to the index after the term.
 */
static int eval_term(const char *expr, int start, const int *vals, int *term_end) {
    int i = start;
    int result = 1;
    while (expr[i] && expr[i] != '+') {
        if (!isupper(expr[i])) { i++; continue; }
        int var = expr[i] - 'A';
        int neg = (expr[i+1] == '\'');
        int v = vals[var];
        if (neg) v = !v;
        result &= v;
        i += neg ? 2 : 1;
    }
    *term_end = i;
    return result;
}

/*
 * Evaluate a full DNF expression for variable assignment vals[].
 * vals[i] = value of variable i (A=0, B=1, ...).
 */
int eval_dnf(const char *expr, const int *vals, int num_vars) {
    (void)num_vars;
    int i = 0;
    int len = (int)strlen(expr);
    while (i < len) {
        int end;
        int t = eval_term(expr, i, vals, &end);
        if (t) return 1;
        i = end;
        if (expr[i] == '+') i++;
    }
    return 0;
}

/* ═══════════════════════════════════════════════════════════════════════════
   SECTION 2 – UNIQUE TABLE (hash table for ROBDD reduction)
   ═══════════════════════════════════════════════════════════════════════════
   Key: (var, low_id, high_id)  →  Value: BDDNode*
   This ensures structural sharing and canonical form.
   ═══════════════════════════════════════════════════════════════════════════ */

#define HTAB_SIZE (1 << 20)   /* must be power of 2 */
#define HTAB_MASK (HTAB_SIZE - 1)

typedef struct HTEntry {
    int var, low_id, high_id;
    BDDNode *node;
    struct HTEntry *next;
} HTEntry;

typedef struct {
    HTEntry *buckets[HTAB_SIZE];
    HTEntry *pool;
    int pool_used;
    int pool_cap;
} UniqueTable;

static UniqueTable *ht_new(void) {
    UniqueTable *t = calloc(1, sizeof(UniqueTable));
    t->pool_cap = 1 << 20;
    t->pool = malloc(t->pool_cap * sizeof(HTEntry));
    t->pool_used = 0;
    return t;
}

static void ht_free(UniqueTable *t) {
    free(t->pool);
    free(t);
}

static unsigned int ht_hash(int var, int low_id, int high_id) {
    unsigned int h = (unsigned int)(var * 2654435761u)
                   ^ (unsigned int)(low_id * 2246822519u)
                   ^ (unsigned int)(high_id * 3266489917u);
    return h & HTAB_MASK;
}

static BDDNode *ht_lookup(UniqueTable *t, int var, int low_id, int high_id) {
    unsigned int h = ht_hash(var, low_id, high_id);
    for (HTEntry *e = t->buckets[h]; e; e = e->next) {
        if (e->var == var && e->low_id == low_id && e->high_id == high_id)
            return e->node;
    }
    return NULL;
}

static void ht_insert(UniqueTable *t, int var, int low_id, int high_id, BDDNode *node) {
    unsigned int h = ht_hash(var, low_id, high_id);
    HTEntry *e = &t->pool[t->pool_used++];
    e->var = var; e->low_id = low_id; e->high_id = high_id; e->node = node;
    e->next = t->buckets[h];
    t->buckets[h] = e;
}

/* ═══════════════════════════════════════════════════════════════════════════
   SECTION 3 – BDD BUILDER (Shannon expansion on DNF, with reduction)
   ═══════════════════════════════════════════════════════════════════════════
   We build the ROBDD top-down using the Shannon expansion:
     f(x1,...,xn) = (NOT xi) * f|xi=0  +  xi * f|xi=1
   Terminal cases:
     - If the restricted DNF is identically 0 → return terminal-0
     - If the restricted DNF is identically 1 → return terminal-1
   Reduction:
     - We use a unique table so that identical sub-BDDs are shared.
     - If low==high we return low directly (elimination rule).
   ═══════════════════════════════════════════════════════════════════════════ */

/* Builder context (per BDD_create call) */
typedef struct {
    UniqueTable *utable;
    BDDNode     *terminal0;
    BDDNode     *terminal1;
    /* node memory pool */
    BDDNode     *pool;
    int          pool_used;
    int          pool_cap;
    int          next_id;
    /* variable ordering: level -> variable index */
    int          order[MAX_VARS];
    int          num_vars;
} BuildCtx;

static BDDNode *alloc_node(BuildCtx *ctx) {
    if (ctx->pool_used >= ctx->pool_cap) {
        ctx->pool_cap *= 2;
        ctx->pool = realloc(ctx->pool, ctx->pool_cap * sizeof(BDDNode));
    }
    return &ctx->pool[ctx->pool_used++];
}

static BuildCtx *ctx_new(const int *order, int num_vars) {
    BuildCtx *ctx = malloc(sizeof(BuildCtx));
    ctx->utable   = ht_new();
    ctx->pool_cap = 1 << 16;
    ctx->pool     = malloc(ctx->pool_cap * sizeof(BDDNode));
    ctx->pool_used = 0;
    ctx->next_id   = 0;
    ctx->num_vars  = num_vars;
    memcpy(ctx->order, order, num_vars * sizeof(int));

    /* terminal nodes */
    ctx->terminal0 = alloc_node(ctx);
    ctx->terminal0->var = -1; ctx->terminal0->val = 0;
    ctx->terminal0->low = ctx->terminal0->high = NULL;
    ctx->terminal0->id  = ctx->next_id++;

    ctx->terminal1 = alloc_node(ctx);
    ctx->terminal1->var = -1; ctx->terminal1->val = 1;
    ctx->terminal1->low = ctx->terminal1->high = NULL;
    ctx->terminal1->id  = ctx->next_id++;

    return ctx;
}

/*
 * Restrict DNF: fix variable `var` to value `val`.
 * Returns a newly allocated string (caller must free).
 */
static char *restrict_dnf(const char *expr, int var, int val) {
    char letter = 'A' + var;
    /* Worst case: same length as original */
    int len = (int)strlen(expr);
    char *out = malloc(len + 2);
    int oi = 0;
    int i = 0;

    while (expr[i]) {
        /* Collect one term */
        int term_start = oi;
        int term_ok = 1;   /* becomes 0 if term is always 0 */
        int term_trivial = 1; /* all literals in term are 1 → term is 1 */

        while (expr[i] && expr[i] != '+') {
            if (!isupper(expr[i])) { i++; continue; }
            char c = expr[i];
            int neg = (expr[i+1] == '\'');
            if (c == letter) {
                /* literal involves our variable */
                int lit_val = neg ? !val : val;
                if (!lit_val) { term_ok = 0; }
                /* either way: skip this literal */
                i += neg ? 2 : 1;
            } else {
                /* keep literal */
                out[oi++] = c;
                if (neg) out[oi++] = '\'';
                term_trivial = 0;
                i += neg ? 2 : 1;
            }
        }

        if (term_ok) {
            if (term_trivial) {
                /* term = 1 → whole function = 1 */
                free(out);
                char *one = malloc(2); one[0]='1'; one[1]='\0';
                return one;
            }
            out[oi++] = '+';
        } else {
            /* term vanishes → rewind output */
            oi = term_start;
        }

        if (expr[i] == '+') i++;
    }

    /* Remove trailing '+' */
    if (oi > 0 && out[oi-1] == '+') oi--;
    if (oi == 0) { out[0]='0'; out[1]='\0'; return out; }
    out[oi] = '\0';
    return out;
}

/*
 * Check if a DNF string is trivially 0 or 1 (after restriction).
 */
static int dnf_is_zero(const char *expr) {
    for (int i = 0; expr[i]; i++)
        if (isupper(expr[i])) return 0;
    return (expr[0] == '0' || expr[0] == '\0');
}

static int dnf_is_one(const char *expr) {
    return (expr[0] == '1');
}

/*
 * get_or_create – look up node in unique table; create if missing.
 * This enforces the two ROBDD reduction rules:
 *   1. Elimination: if low == high, return low (skip this variable).
 *   2. Merging: if (var,low,high) already exists, return existing node.
 */
static BDDNode *get_or_create(BuildCtx *ctx, int var, BDDNode *low, BDDNode *high) {
    /* Rule 1: elimination */
    if (low == high) return low;

    /* Rule 2: merging via unique table */
    BDDNode *existing = ht_lookup(ctx->utable, var, low->id, high->id);
    if (existing) return existing;

    BDDNode *n = alloc_node(ctx);
    n->var  = var;
    n->val  = -1;
    n->low  = low;
    n->high = high;
    n->id   = ctx->next_id++;

    ht_insert(ctx->utable, var, low->id, high->id, n);
    return n;
}

/*
 * Recursive BDD builder.
 * level: current level in the ordering (0 = top).
 * expr:  current (restricted) DNF expression.
 */
static BDDNode *build_rec(BuildCtx *ctx, int level, const char *expr) {
    /* Terminal cases */
    if (dnf_is_zero(expr)) return ctx->terminal0;
    if (dnf_is_one(expr))  return ctx->terminal1;
    if (level >= ctx->num_vars) {
        /* No more variables to branch on; evaluate */
        return ctx->terminal0; /* should not normally reach here */
    }

    int var = ctx->order[level];

    char *expr0 = restrict_dnf(expr, var, 0);
    char *expr1 = restrict_dnf(expr, var, 1);

    BDDNode *low  = build_rec(ctx, level + 1, expr0);
    BDDNode *high = build_rec(ctx, level + 1, expr1);

    free(expr0);
    free(expr1);

    return get_or_create(ctx, var, low, high);
}

/* ═══════════════════════════════════════════════════════════════════════════
   SECTION 4 – COUNT NODES (exclude terminals)
   ═══════════════════════════════════════════════════════════════════════════ */

static void count_nodes_rec(BDDNode *n, int *visited, int *count) {
    if (!n || n->var == -1) return;
    if (visited[n->id]) return;
    visited[n->id] = 1;
    (*count)++;
    count_nodes_rec(n->low,  visited, count);
    count_nodes_rec(n->high, visited, count);
}

static int count_nodes(BDDNode *root, int max_id) {
    int *visited = calloc(max_id + 1, sizeof(int));
    int count = 0;
    count_nodes_rec(root, visited, &count);
    free(visited);
    return count;
}

/* ═══════════════════════════════════════════════════════════════════════════
   SECTION 5 – PUBLIC API
   ═══════════════════════════════════════════════════════════════════════════ */

BDD *BDD_create(const char *bfunkcia, const char *poradie) {
    if (!bfunkcia || !poradie) return NULL;

    /* Parse variable ordering from the order string (e.g. "ABCD") */
    int order[MAX_VARS];
    int num_vars = 0;
    for (int i = 0; poradie[i] && num_vars < MAX_VARS; i++) {
        if (isupper(poradie[i]))
            order[num_vars++] = poradie[i] - 'A';
    }
    if (num_vars == 0) return NULL;

    BuildCtx *ctx = ctx_new(order, num_vars);

    BDDNode *root = build_rec(ctx, 0, bfunkcia);

    int node_count = count_nodes(root, ctx->next_id);

    BDD *bdd = malloc(sizeof(BDD));
    bdd->num_vars  = num_vars;
    bdd->size      = node_count;
    bdd->root      = root;
    memcpy(bdd->order, order, num_vars * sizeof(int));

    /* Transfer pool ownership */
    bdd->pool      = ctx->pool;
    bdd->pool_size = ctx->pool_used;
    bdd->pool_cap  = ctx->pool_cap;

    ht_free(ctx->utable);
    free(ctx);
    return bdd;
}

BDD *BDD_create_with_best_order(const char *bfunkcia) {
    if (!bfunkcia) return NULL;

    /* Discover variables */
    int vars[MAX_VARS];
    int num_vars = parse_vars(bfunkcia, vars);
    if (num_vars == 0) return NULL;

    /* We need at least N unique orderings tried.
       Strategy: cyclic shifts (N orderings) + random shuffles up to max(N, 10) total */
    int max_tries = num_vars < 10 ? 10 : num_vars;
    if (num_vars <= 8) {
        /* For small N try all N! permutations */
        int fact = 1;
        for (int i = 1; i <= num_vars; i++) fact *= i;
        if (fact < max_tries) max_tries = fact;
    }

    /* Helper: build order string from permutation of vars[] */
    char order_str[MAX_VARS + 2];

    BDD *best = NULL;

    /* ── Generate orderings ─────────────────────────────────────────────── */

    /* Seen-set: store hashes of tried permutations to avoid duplicates */
    unsigned long long seen[2048] = {0};
    int seen_count = 0;
    int seen_cap = 2048;

    /* perm[]: current permutation of indices into vars[] */
    int perm[MAX_VARS];
    for (int i = 0; i < num_vars; i++) perm[i] = i;

    srand((unsigned)time(NULL));

    int tries = 0;

    /* Generate all permutations for small N, random for large N */
    /* We use Heap's algorithm for small N */

    if (num_vars <= 8) {
        /* Heap's algorithm */
        int c[MAX_VARS] = {0};
        /* Try initial permutation */
        for (int i = 0; i < num_vars; i++) order_str[i] = 'A' + vars[perm[i]];
        order_str[num_vars] = '\0';
        BDD *b = BDD_create(bfunkcia, order_str);
        if (!best || (b && b->size < best->size)) { if(best) BDD_free(best); best = b; }
        else if(b) BDD_free(b);
        tries++;

        int k = 0;
        while (k < num_vars && tries < max_tries) {
            if (c[k] < k) {
                if (k % 2 == 0) { int t = perm[0]; perm[0]=perm[k]; perm[k]=t; }
                else            { int t = perm[c[k]]; perm[c[k]]=perm[k]; perm[k]=t; }

                for (int i = 0; i < num_vars; i++) order_str[i] = 'A' + vars[perm[i]];
                order_str[num_vars] = '\0';

                BDD *b2 = BDD_create(bfunkcia, order_str);
                if (!best || (b2 && b2->size < best->size)) { if(best) BDD_free(best); best = b2; }
                else if(b2) BDD_free(b2);
                tries++;
                c[k]++;
                k = 0;
            } else {
                c[k] = 0;
                k++;
            }
        }
    } else {
        /* For large N: cyclic shifts + random shuffles */
        for (int shift = 0; shift < num_vars && tries < max_tries; shift++) {
            for (int i = 0; i < num_vars; i++)
                order_str[i] = 'A' + vars[(i + shift) % num_vars];
            order_str[num_vars] = '\0';

            BDD *b = BDD_create(bfunkcia, order_str);
            if (!best || (b && b->size < best->size)) { if(best) BDD_free(best); best = b; }
            else if(b) BDD_free(b);
            tries++;
        }

        /* Fill remaining with random shuffles */
        while (tries < max_tries) {
            for (int i = 0; i < num_vars; i++) perm[i] = i;
            for (int i = num_vars - 1; i > 0; i--) {
                int j = rand() % (i + 1);
                int t = perm[i]; perm[i] = perm[j]; perm[j] = t;
            }
            for (int i = 0; i < num_vars; i++) order_str[i] = 'A' + vars[perm[i]];
            order_str[num_vars] = '\0';

            /* Simple duplicate check via string hash */
            unsigned long long h = 0;
            for (int i = 0; order_str[i]; i++) h = h * 31 + order_str[i];
            int dup = 0;
            for (int i = 0; i < seen_count && i < seen_cap; i++)
                if (seen[i] == h) { dup = 1; break; }
            if (dup) continue;
            if (seen_count < seen_cap) seen[seen_count++] = h;

            BDD *b = BDD_create(bfunkcia, order_str);
            if (!best || (b && b->size < best->size)) { if(best) BDD_free(best); best = b; }
            else if(b) BDD_free(b);
            tries++;
        }
    }

    return best;
}

char BDD_use(BDD *bdd, const char *vstupy) {
    if (!bdd || !vstupy) return (char)-1;
    if ((int)strlen(vstupy) < bdd->num_vars) return (char)-1;

    /* Build variable value map: var_index -> value */
    int vals[MAX_VARS];
    for (int i = 0; i < bdd->num_vars; i++) {
        char c = vstupy[i];
        if (c != '0' && c != '1') return (char)-1;
        vals[bdd->order[i]] = c - '0';
    }

    BDDNode *cur = bdd->root;
    while (cur && cur->var != -1) {
        int v = vals[cur->var];
        cur = v ? cur->high : cur->low;
    }
    if (!cur) return (char)-1;
    return (char)('0' + cur->val);
}

void BDD_free(BDD *bdd) {
    if (!bdd) return;
    free(bdd->pool);
    free(bdd);
}