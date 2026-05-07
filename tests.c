#include "bdd.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

/* ═══════════════════════════════════════════════════════════════════════════
   HELPER: random DNF generator
   ═══════════════════════════════════════════════════════════════════════════ */
static void random_dnf(int num_vars, int num_terms, char *buf)
{
    int pos = 0;
    for (int t = 0; t < num_terms; t++)
    {
        if (t > 0)
            buf[pos++] = '+';
        int any = 0;
        for (int v = 0; v < num_vars; v++)
        {
            int r = rand() % 3; /* 0: skip, 1: positive, 2: negated */
            if (r == 0)
                continue;
            any = 1;
            buf[pos++] = 'A' + v;
            if (r == 2)
                buf[pos++] = '\'';
        }
        if (!any)
        {
            buf[pos++] = 'A';
        }
    }
    buf[pos] = '\0';
}

static void default_order(int num_vars, char *buf)
{
    for (int i = 0; i < num_vars; i++)
        buf[i] = 'A' + i;
    buf[num_vars] = '\0';
}

static long long full_bdd_size(int n)
{
    return (1LL << n) - 1;
}

/*
 * Build the input string for BDD_use given:
 *   vals[v] = value of variable v (A=0,B=1,...)
 *   bdd->order[i] = which variable is at position i
 * So vstupy[i] = vals[ bdd->order[i] ]
 */
static void make_input_str(const int *vals, BDD *bdd, char *buf)
{
    for (int i = 0; i < bdd->num_vars; i++)
        buf[i] = '0' + vals[bdd->order[i]];
    buf[bdd->num_vars] = '\0';
}

/* ═══════════════════════════════════════════════════════════════════════════
   TEST 1 – Basic correctness (example from the assignment)
   ═══════════════════════════════════════════════════════════════════════════ */
static void test_basic(void)
{
    printf("===========================================\n");
    printf("TEST 1: Basic correctness  (f = AB+C)\n");
    printf("===========================================\n");

    BDD *bdd = BDD_create("AB+C", "ABC");
    if (!bdd)
    {
        printf("FAIL: BDD_create returned NULL\n");
        return;
    }
    printf("BDD size (internal nodes): %d\n", bdd->size);

    /* Truth table for f(A,B,C) = AB + C */
    const char *inputs[] = {"000", "001", "010", "011", "100", "101", "110", "111"};
    const char expected[] = {'0', '1', '0', '1', '0', '1', '1', '1'};

    int errors = 0;
    for (int i = 0; i < 8; i++)
    {
        char res = BDD_use(bdd, inputs[i]);
        if (res != expected[i])
        {
            printf("  FAIL: input=%s  got=%c  expected=%c\n",
                   inputs[i], res, expected[i]);
            errors++;
        }
    }
    if (errors == 0)
        printf("  All 8 cases PASSED.\n");
    BDD_free(bdd);
    printf("\n");
}

/* ═══════════════════════════════════════════════════════════════════════════
   TEST 2 – Exhaustive correctness
   ═══════════════════════════════════════════════════════════════════════════ */
static void test_exhaustive(int num_vars, int num_functions)
{
    printf("===========================================\n");
    printf("TEST 2: Exhaustive  N=%d  functions=%d\n", num_vars, num_functions);
    printf("===========================================\n");

    long long total_inputs = 1LL << num_vars;
    int total_errors = 0;
    double total_reduction = 0.0;
    double total_size = 0.0;
    double time_create = 0.0;
    double time_use = 0.0;
    long long full_size = full_bdd_size(num_vars);

    char dnf[8192];
    char order[MAX_VARS + 1];
    char input_str[MAX_VARS + 1];

    default_order(num_vars, order);

    for (int fn = 0; fn < num_functions; fn++)
    {
        int num_terms = 1 + rand() % (num_vars + 2);
        random_dnf(num_vars, num_terms, dnf);

        clock_t t0 = clock();
        BDD *bdd = BDD_create(dnf, order);
        clock_t t1 = clock();
        time_create += (double)(t1 - t0) / CLOCKS_PER_SEC;

        if (!bdd)
        {
            printf("  FAIL: BDD_create NULL for: %s\n", dnf);
            continue;
        }

        total_size += bdd->size;
        double reduction = (full_size > 0)
                               ? 100.0 * (1.0 - (double)bdd->size / (double)full_size)
                               : 0.0;
        total_reduction += reduction;

        int vals[MAX_VARS];
        clock_t t2 = clock();
        for (long long mask = 0; mask < total_inputs; mask++)
        {
            for (int v = 0; v < num_vars; v++)
                vals[v] = (mask >> (num_vars - 1 - v)) & 1;

            make_input_str(vals, bdd, input_str);
            char got = BDD_use(bdd, input_str);
            int expected = eval_dnf(dnf, vals, num_vars);

            if (got != '0' + expected)
            {
                if (total_errors < 5)
                    printf("  FAIL fn=%d: DNF=%s  input=%s  got=%c  expected=%c\n",
                           fn, dnf, input_str, got, '0' + expected);
                total_errors++;
            }
        }
        clock_t t3 = clock();
        time_use += (double)(t3 - t2) / CLOCKS_PER_SEC;

        BDD_free(bdd);
    }

    printf("  Errors:              %d / %lld\n",
           total_errors, (long long)num_functions * total_inputs);
    printf("  Avg BDD size:        %.1f nodes  (full tree: %lld)\n",
           total_size / num_functions, full_size);
    printf("  Avg reduction:       %.1f%%\n", total_reduction / num_functions);
    printf("  Total BDD_create:    %.4f s  (avg %.6f s)\n",
           time_create, time_create / num_functions);
    printf("  Total BDD_use:       %.4f s  (avg per call %.3f us)\n",
           time_use, 1e6 * time_use / ((double)num_functions * (double)total_inputs));
    if (total_errors == 0)
        printf("  ALL CORRECT\n");
    printf("\n");
}

/* ═══════════════════════════════════════════════════════════════════════════
   TEST 3 – BDD_create vs BDD_create_with_best_order
   ═══════════════════════════════════════════════════════════════════════════ */
static void test_best_order(int num_vars, int num_functions)
{
    printf("===========================================\n");
    printf("TEST 3: Best-order benefit  N=%d  functions=%d\n", num_vars, num_functions);
    printf("===========================================\n");

    char dnf[8192];
    char order[MAX_VARS + 1];
    default_order(num_vars, order);

    double sum_create = 0, sum_best = 0;
    double time_create = 0, time_best = 0;
    int errors = 0;
    long long total_inputs = 1LL << num_vars;
    char input_str[MAX_VARS + 1];
    int vals[MAX_VARS];

    for (int fn = 0; fn < num_functions; fn++)
    {
        int num_terms = 1 + rand() % (num_vars + 2);
        random_dnf(num_vars, num_terms, dnf);

        clock_t t0 = clock();
        BDD *b1 = BDD_create(dnf, order);
        clock_t t1 = clock();
        time_create += (double)(t1 - t0) / CLOCKS_PER_SEC;

        clock_t t2 = clock();
        BDD *b2 = BDD_create_with_best_order(dnf);
        clock_t t3 = clock();
        time_best += (double)(t3 - t2) / CLOCKS_PER_SEC;

        if (!b1 || !b2)
        {
            if (b1)
                BDD_free(b1);
            if (b2)
                BDD_free(b2);
            continue;
        }

        sum_create += b1->size;
        sum_best += b2->size;

        /* Verify b2 is also correct for all inputs */
        for (long long mask = 0; mask < total_inputs; mask++)
        {
            for (int v = 0; v < num_vars; v++)
                vals[v] = (mask >> (num_vars - 1 - v)) & 1;

            /* Build input string matching b2->order */
            make_input_str(vals, b2, input_str);
            char got = BDD_use(b2, input_str);
            int expected = eval_dnf(dnf, vals, num_vars);
            if (got != '0' + expected)
                errors++;
        }

        BDD_free(b1);
        BDD_free(b2);
    }

    double avg_create = sum_create / num_functions;
    double avg_best = sum_best / num_functions;
    double extra_red = (avg_create > 0)
                           ? 100.0 * (avg_create - avg_best) / avg_create
                           : 0.0;

    printf("  Avg nodes (default order): %.1f\n", avg_create);
    printf("  Avg nodes (best order):    %.1f\n", avg_best);
    printf("  Extra reduction via order: %.1f%%\n", extra_red);
    printf("  Time BDD_create (total):   %.4f s\n", time_create);
    printf("  Time best_order (total):   %.4f s\n", time_best);
    printf("  Correctness errors:        %d\n", errors);
    if (errors == 0)
        printf("  ALL CORRECT\n");
    printf("\n");
}

/* ═══════════════════════════════════════════════════════════════════════════
   TEST 4 – Reduction stats table per N
   ═══════════════════════════════════════════════════════════════════════════ */
static void test_reduction_stats(void)
{
    printf("===========================================\n");
    printf("TEST 4: Reduction stats by N (100 functions each)\n");
    printf("  N | Full nodes | Avg BDD nodes | Reduction%%\n");
    printf("----+------------+---------------+-----------\n");

    char dnf[8192];
    char order[MAX_VARS + 1];

    for (int n = 2; n <= 13; n++)
    {
        default_order(n, order);
        long long full = full_bdd_size(n);
        double sum_size = 0;
        int fn_count = 100;

        for (int fn = 0; fn < fn_count; fn++)
        {
            int num_terms = 1 + rand() % (n + 2);
            random_dnf(n, num_terms, dnf);
            BDD *b = BDD_create(dnf, order);
            if (b)
            {
                sum_size += b->size;
                BDD_free(b);
            }
        }

        double avg = sum_size / fn_count;
        double red = 100.0 * (1.0 - avg / (double)full);
        printf(" %2d | %10lld | %13.1f | %8.1f%%\n", n, full, avg, red);
    }
    printf("\n");
}

/* ═══════════════════════════════════════════════════════════════════════════
   MAIN
   ═══════════════════════════════════════════════════════════════════════════ */
int main(void)
{
    srand((unsigned)time(NULL));

    printf("\n");
    printf("*=================================================*\n");
    printf("|       BDD Implementation -- Test Suite          |\n");
    printf("*=================================================*\n\n");

    test_basic();

    for (int n = 2; n <= 13; n++)
        test_exhaustive(n, 100);

    for (int n = 3; n <= 8; n++)
        test_best_order(n, 100);

    test_reduction_stats();

    printf("Done.\n");
    return 0;
}