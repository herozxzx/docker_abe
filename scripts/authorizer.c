#include <pbc.h>
#include <pbc_test.h>
#include <stdio.h>
#include <math.h>
#include <gmp.h>
#include <string.h>
#include "message_handle.h"

typedef struct
{
    field_t Fq, Fq2, Eq;
    int exp2, exp1;
    int sign1;
} * a_pairing_data_ptr;

typedef struct
{
    field_ptr field; // The field where the curve is defined.
    element_t a, b;  // The curve is E: Y^2 = X^3 + a X + b.
    // cofac == NULL means we're using the whole group of points.
    // otherwise we're working in the subgroup of order #E / cofac,
    // where #E is the number of points in E.
    mpz_ptr cofac;
    // A generator of E.
    element_t gen_no_cofac;
    // A generator of the subgroup.
    element_t gen;
    // A non-NULL quotient_cmp means we are working with the quotient group of
    // order #E / quotient_cmp, and the points are actually coset
    // representatives. Thus for a comparison, we must multiply by quotient_cmp
    // before comparing.
    mpz_ptr quotient_cmp;
} * curve_data_ptr;

typedef struct
{
    element_t sk1;
    element_t sk2;
    element_t sk3;
} secrets;

typedef struct
{
    int i;
    int users[1000];
    element_t hi;
    element_t pk;
    secrets sk;
} auth;

typedef struct
{
    element_t c1;
    element_t c2;
    element_t c3;
    element_t c4;
} cipher;

pairing_t pairing;
element_t P, Q, R, S;

void Setup(int argc, char **argv)
{
    printf("Setup Starting\n");
    pbc_demo_pairing_init(pairing, argc, argv);
    if (!pairing_is_symmetric(pairing))
        pbc_die("pairing must be symmetric");
    element_init_G1(P, pairing);
    element_init_G1(Q, pairing);
    element_init_G1(R, pairing);
    element_init_G1(S, pairing);

    element_t q, r, s;
    element_init_Zr(q, pairing);
    element_init_Zr(r, pairing);
    element_init_Zr(s, pairing);
    element_random(q);
    element_random(r);
    element_random(s);

    element_set(P, ((curve_data_ptr)((a_pairing_data_ptr)pairing->data)->Eq->data)->gen);
    element_pow_zn(Q, P, q);
    element_pow_zn(R, P, q);
    element_pow_zn(S, P, q);

    element_clear(q);
    element_clear(r);
    element_clear(s);
    printf("Setup Finished\n");
}

auth CreateAuthority(int i)
{
    auth aa;
    element_t x, y;
    element_init_Zr(x, pairing);
    element_init_Zr(y, pairing);
    element_init_Zr(aa.hi, pairing);
    element_random(x);
    element_random(y);
    element_random(aa.hi);

    element_init_G1(aa.pk, pairing);
    element_pow_zn(aa.pk, P, aa.hi);

    element_t tmp;
    element_init_G1(aa.sk.sk1, pairing);
    element_init_G1(aa.sk.sk2, pairing);
    element_init_G1(aa.sk.sk3, pairing);
    element_init_G1(tmp, pairing);
    // sk1
    element_pow_zn(tmp, R, y);
    element_mul(tmp, Q, tmp);
    element_pow_zn(tmp, tmp, aa.hi);
    element_pow_zn(aa.sk.sk1, S, x);
    element_mul(aa.sk.sk1, tmp, aa.sk.sk1);
    // sk2
    element_pow_zn(aa.sk.sk2, P, x);
    // sk3
    element_pow_zn(aa.sk.sk3, P, y);
    element_pow_zn(aa.sk.sk3, aa.sk.sk3, aa.hi);

    aa.i = i;
    for (int u = 0; u < 1000; u++)
        aa.users[u] = 0;

    element_clear(x);
    element_clear(y);
    element_clear(tmp);

    return aa;
}

element_t *RequestAttributePK()
{
    element_t tmp;

    element_t *pk = &tmp;
    return pk;
}

secrets RequestAttributeSK()
{
    secrets sk;

    return sk;
}

void Encrypt(cipher *ciphertexts, char *raw_m, int condc, int *condn, int (*conds)[100], element_t **pk)
{
    printf("\nStarting Encyption\n Message: ");
    char message_dec[2048], message[2048];
    mpz_t message_mpz;
    element_t message_ele;
    element_init_GT(message_ele, pairing);
    mpz_init(message_mpz);
    puts(raw_m);
    messageToValue(raw_m, message_mpz, message_dec);
    strcpy(message, "[");
    strcat(message, message_dec);
    strcat(message, ",0]");
    element_set_str(message_ele, message, 10);
    element_printf("enc_m = %B\n", message_ele);

    for (int l = 0; l < condc; l++)
    {
        int auth = conds[l][0];
        element_t *pkl = pk[conds[l][0]];
        printf("Auth = %d (%d, 0)\n", auth, l);
        // element_printf("PK = %B\n", *pkl);
        for (int m = 1; m < condn[l]; m++)
        {
            int auth = conds[l][m];
            element_t *curr_pk = pk[auth];
            printf("Auth = %d (%d, %d)\n", auth, l, m);
            // element_printf("PK = %B\n", *curr_pk);
            element_mul(*pkl, *pkl, *curr_pk);
        }
        
        element_t pair1, r, tmp;
        element_init_GT(ciphertexts[l].c1, pairing);
        element_init_G1(ciphertexts[l].c2, pairing);
        element_init_G1(ciphertexts[l].c3, pairing);
        element_init_G1(ciphertexts[l].c4, pairing);
        element_init_G1(tmp, pairing);
        element_init_GT(pair1, pairing);
        element_init_Zr(r, pairing);

        element_random(r);
        element_pairing(pair1, Q, *pkl);
        element_pow_zn(pair1, pair1, r);
        element_mul(ciphertexts[l].c1, message_ele, pair1);
        element_pow_zn(ciphertexts[l].c2, P, r);
        element_pow_zn(ciphertexts[l].c3, R, r);
        element_pow_zn(ciphertexts[l].c4, S, r);
    }
    printf("Finished Encyption\n");
}

void Decrypt(char *dec_m, cipher ciphertext, int condn, int *cond, secrets *sk)
{
    printf("\nStarting Decryption\n");
    char message_dec[2048], message[2048];
    mpz_t message_mpz;
    mpz_init(message_mpz);
    int auth = cond[0];
    secrets sks = sk[auth];
    printf("Auth = %d (%d)\n", auth, 0);
    // element_printf("SK1 = %B\n", sks.sk1);
    // element_printf("SK2 = %B\n", sks.sk2);
    // element_printf("SK3 = %B\n", sks.sk3);
    for (int m = 1; m < condn; m++)
    {
        int auth = cond[m];
        printf("Auth = %d (%d)\n", auth, m);
        // element_printf("SK = %B\n", sk[auth].sk1);
        element_mul(sks.sk1, sks.sk1, sk[auth].sk1);
        element_mul(sks.sk2, sks.sk2, sk[auth].sk2);
        element_mul(sks.sk3, sks.sk3, sk[auth].sk3);
    }    
    // element_printf("SK1 = %B\n", sks.sk1);
    // element_printf("SK2 = %B\n", sks.sk2);
    // element_printf("SK3 = %B\n", sks.sk3);
    // element_printf("PK1 = %B\n", ciphertext.c1);
    // element_printf("PK2 = %B\n", ciphertext.c2);
    // element_printf("PK3 = %B\n", ciphertext.c3);
    // element_printf("PK4 = %B\n", ciphertext.c4);
    element_t pair1, pair2, pair3, m;
    element_init_GT(pair1, pairing);
    element_init_GT(pair2, pairing);
    element_init_GT(pair3, pairing);
    element_init_GT(m, pairing);
    element_pairing(pair1, ciphertext.c3, sks.sk3);
    element_pairing(pair2, ciphertext.c4, sks.sk2);
    element_pairing(pair3, ciphertext.c2, sks.sk1);
    element_mul(m, ciphertext.c1, pair1);
    element_mul(m, m, pair2);
    element_div(m, m, pair3);
    element_printf("dec_m = %B\n", m);
    element_to_mpz(message_mpz, m);
    valueToMessage(dec_m, message_mpz);

    printf("Decrypt successfully. The message is:\n");
    puts(dec_m);
}

int main(int argc, char **argv)
{
    Setup(argc, argv);
    auth aa1 = CreateAuthority(1);
    auth aa2 = CreateAuthority(2);
    auth aa3 = CreateAuthority(3);

    // Encrypt (condition = 1 or (2 and 3))
    char raw_message[2048] = "Test message\0";
    int condc = 2;
    int condn[2] = {1, 2};
    int conds[2][100];
    conds[0][0] = 1;
    conds[1][0] = 2; conds[1][1] = 3;
    int max = 3;
    element_t *pk[max + 1];
    for (int k = 0; k < max + 1; k++)
        pk[k] = (element_t *)malloc(sizeof(element_t));
    pk[1] = &aa1.pk; // attr 1
    pk[2] = &aa2.pk; // attr 2
    pk[3] = &aa3.pk; // attr 3
    cipher ct[max + 1];
    Encrypt(ct, raw_message, condc, condn, conds, pk);

    // Decrypt user (id = u1, auths = 1)
    int user1_a[1] = {1};
    int condc1 = 1;
    int cond1[1] = {1};
    max = 1;
    secrets sk1[max + 1];
    sk1[user1_a[0]] = aa1.sk;
    char dec_message1[2048];
    Decrypt(dec_message1, ct[0], condc1, cond1, sk1);

    // Decrypt user (id = u2, auths = 2, 3)
    int user2_a[2] = {2, 3};
    int condc2 = 2;
    int cond2[2] = {2, 3};
    max = 3;
    secrets sk2[max + 1];
    sk2[user2_a[0]] = aa2.sk;
    sk2[user2_a[1]] = aa3.sk;
    char dec_message2[2048];
    Decrypt(dec_message2, ct[1], condc2, cond2, sk2);
}
