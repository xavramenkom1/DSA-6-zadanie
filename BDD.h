#ifndef BDD_H
#define BDD_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ─── DNF representation ─────────────────────────────────────────────────────
   Input format:  AB'+CD   means (A AND NOT-B) OR (C AND D)
   Variables:     uppercase letters A-Z
   Negation:      letter followed by '
   Terms (AND):   concatenated letters, e.g. AB'C
   Sum (OR):      terms separated by +
   ─────────────────────────────────────────────────────────────────────────── */

#define MAX_VARS  26
#define MAX_NODES 1000000

/* ── BDD node ──────────────────────────────────────────────────────────────── */
typedef struct BDDNode {
    int var;               /* variable index (0..N-1); -1 for terminal nodes  */
    int val;               /* terminal value: 0 or 1 (valid only if var == -1)*/
    struct BDDNode *low;   /* child when var = 0                              */
    struct BDDNode *high;  /* child when var = 1                              */
    int id;                /* unique node id (for reduction hash table)       */
} BDDNode;

/* ── BDD structure ─────────────────────────────────────────────────────────── */
typedef struct {
    int     num_vars;      /* number of variables in the Boolean function     */
    int     size;          /* number of nodes (excluding terminal nodes)      */
    BDDNode *root;         /* pointer to root node                            */
    int     order[MAX_VARS]; /* order[i] = variable index at BDD level i     */
    /* memory pool */
    BDDNode *pool;
    int      pool_size;
    int      pool_cap;
} BDD;

/* ─── Public API ──────────────────────────────────────────────────────────── */

/*
 * BDD_create – build a reduced BDD from a DNF expression string.
 *
 * @param bfunkcia   DNF expression, e.g. "AB'+CD"
 * @param poradie    variable ordering string, e.g. "ABCD" means A first, etc.
 * @return           pointer to the newly created BDD, or NULL on error
 */
BDD *BDD_create(const char *bfunkcia, const char *poradie);

/*
 * BDD_create_with_best_order – try multiple variable orderings and return the
 * BDD with the fewest nodes.  At least N unique orderings are tried.
 *
 * @param bfunkcia   DNF expression string
 * @return           pointer to the best BDD found
 */
BDD *BDD_create_with_best_order(const char *bfunkcia);

/*
 * BDD_use – evaluate the BDD for a concrete input assignment.
 *
 * @param bdd     pointer to the BDD
 * @param vstupy  assignment string, index i = value of variable order[i]
 *                e.g. "1010" sets order[0]=1, order[1]=0, order[2]=1, order[3]=0
 * @return        '0' or '1'; -1 on error
 */
char BDD_use(BDD *bdd, const char *vstupy);

/*
 * BDD_free – release all memory associated with a BDD.
 */
void BDD_free(BDD *bdd);

/*
 * eval_dnf – evaluate a DNF expression for a given variable assignment.
 * vars[i] = value (0 or 1) of variable i  (A=0, B=1, …)
 *
 * @return 0 or 1; -1 on parse error
 */
int eval_dnf(const char *expr, const int *vars, int num_vars);

#endif /* BDD_H */