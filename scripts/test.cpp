#include <pbc.h>
#include <pbc_test.h>
#include <gmp.h>
#include <bits/stdc++.h>
#include <chrono>

using namespace std;

#define NC "\e[0m"
#define RED "\e[0;31m"
#define GRN "\e[0;32m"

typedef struct
{
    field_t Fq, Fq2, Eq;
    int exp2, exp1;
    int sign1;
} *a_pairing_data_ptr;

typedef struct
{
    field_ptr field;
    element_t a, b;
    mpz_ptr cofac;
    element_t gen_no_cofac;
    element_t gen;
    mpz_ptr quotient_cmp;
} *curve_data_ptr;

void Test(int argc, char **argv, int rp)
{
    cout << "Calculating test (avg: " << rp << ")" << endl;
    chrono::system_clock::time_point start, end;
    double time;
    pairing_t pairing;
    element_t P, eP, r;
    pbc_demo_pairing_init(pairing, argc, argv);
    if (!pairing_is_symmetric(pairing))
        pbc_die(" Error: pairing must be symmetric");
    element_init_Zr(r, pairing);
    element_random(r);
    element_init_G1(P, pairing);
    element_set(P, ((curve_data_ptr)((a_pairing_data_ptr)pairing->data)->Eq->data)->gen);
    element_init_GT(eP, pairing);
    element_pairing(eP, P, P);

    cout << " Testing |Cg| pow-cost" << endl;
    start = chrono::system_clock::now();
    for(int i = 0; i < rp; i++){
        element_pow_zn(P, P, r);
    }
    end = chrono::system_clock::now();
    time = static_cast<double>(chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0);
    printf("  time %lf[ms]\n\n", time / rp);

    cout << " Testing |Cgt| pow-cost" << endl;
    start = chrono::system_clock::now();
    for(int i = 0; i < rp; i++){
        element_pow_zn(eP, eP, r);
    }
    end = chrono::system_clock::now();
    time = static_cast<double>(chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0);
    printf("  time %lf[ms]\n\n", time / rp);

    cout << " Testing |Cpair| pow-cost" << endl;
    start = chrono::system_clock::now();
    for(int i = 0; i < rp; i++){
        element_pairing(eP, P, P);
    }
    end = chrono::system_clock::now();
    time = static_cast<double>(chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0);
    printf("  time %lf[ms]\n", time / rp);

    element_clear(P);
    element_clear(eP);
    cout << "Finishing test" << endl;
}

int main(int argc, char **argv)
{
    // Setup (attribute = 0 ~ 13)
    Test(argc, argv, 100);

}
