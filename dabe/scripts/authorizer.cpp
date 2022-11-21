#include "message_handle.h"
#include <pbc.h>
#include <pbc_test.h>
#include <gmp.h>
#include <bits/stdc++.h>
#include <chrono>
#include <thread>

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
} *curve_data_ptr;

typedef struct
{
    element_t sk1;
    element_t sk2;
    element_t sk3;
} secrets;

typedef struct
{
    int attr;
    vector<int> users;
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

typedef struct
{
    int condc;   // number of conditions
    int *condn;  // number of attributes in each condition
    int **conds; // conditions
    int max;     // max attribute
} condition;

pairing_t pairing;
element_t P, Q, R, S;
vector<element_t> T;

void Setup(int argc, char **argv)
{
    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();
    cout << "Starting Setup" << endl;
    pbc_demo_pairing_init(pairing, argc, argv);
    if (!pairing_is_symmetric(pairing))
        pbc_die(" Error: pairing must be symmetric");
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
    end = chrono::system_clock::now();
    double time = static_cast<double>(chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0);
    printf("Finished Setup (time %lf[ms])\n", time);
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

    // pk
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

    aa.attr = i;

    element_clear(x);
    element_clear(y);
    element_clear(tmp);

    return aa;
}
/*
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
*/
bool Check_Message(vector<char[2048]> &dec_m_tmp, char *dec_m)
{
    bool res = false;
    for (char *m : dec_m_tmp)
    {
        try
        {
            string m_head(m, 8);
            string m_str(m, 8, string(m).size() - 8);
            if (m_head == "message:")
            {
                strcpy(dec_m, m_str.c_str());
                res = true;
            }
        }
        catch (...)
        {
        }
    }
    return res;
}

void Encrypt(vector<cipher> &ciphertexts, char *raw_m, vector<vector<int>> &w, map<int, element_t *> &pks)
{
    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();
    cout << "\nStarting Encyption\n Message = " << string(raw_m, 8, string(raw_m).size() - 8) << endl;
    char message_dec[2048], message[2048];
    mpz_t message_mpz;
    element_t message_ele;
    element_init_GT(message_ele, pairing);
    mpz_init(message_mpz);
    messageToValue(raw_m, message_mpz, message_dec);
    strcpy(message, "[");
    strcat(message, message_dec);
    strcat(message, ",0]");
    element_set_str(message_ele, message, 10);
    // element_printf(" enc_m = %B\n", message_ele);
    cout << " W = ";
    int count1 = 0;
    for (vector<int> &cond : w)
    {
        if (count1 != 0)
            cout << " or ";
        cout << "(";
        element_t pkl;
        element_init_G1(pkl, pairing);
        int count2 = 0;
        for (int &attr : cond)
        {
            if (count2 == 0)
            {
                cout << attr;
                element_set(pkl, *pks[attr]);
            }
            else
            {
                cout << " and " << attr;
                element_t *curr_pk = pks[attr];
                element_mul(pkl, pkl, *curr_pk);
            }
            // element_printf("PK = %B\n", *curr_pk);
            count2++;
        }
        cout << ")";
        element_t pair1, r, tmp;
        element_init_GT(ciphertexts[count1].c1, pairing);
        element_init_G1(ciphertexts[count1].c2, pairing);
        element_init_G1(ciphertexts[count1].c3, pairing);
        element_init_G1(ciphertexts[count1].c4, pairing);
        element_init_G1(tmp, pairing);
        element_init_GT(pair1, pairing);
        element_init_Zr(r, pairing);

        element_random(r);
        element_pairing(pair1, Q, pkl);
        element_pow_zn(pair1, pair1, r);
        element_mul(ciphertexts[count1].c1, message_ele, pair1);
        element_pow_zn(ciphertexts[count1].c2, P, r);
        element_pow_zn(ciphertexts[count1].c3, R, r);
        element_pow_zn(ciphertexts[count1].c4, S, r);

        element_clear(pkl);
        element_clear(pair1);
        element_clear(r);
        element_clear(tmp);
        count1++;
    }
    end = chrono::system_clock::now();
    double time = static_cast<double>(chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0);
    printf("\nFinished Encyption (time %lf[ms])\n", time);
}

void Decrypt(cipher &ciphertext, char *dec_m, vector<int> &cond, map<int, secrets *> &sks, bool out)
{
    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();
    if (out)
        cout << "\nStarting Decryption" << endl;
    char message_dec[2048];
    mpz_t message_mpz;
    mpz_init(message_mpz);
    if (out)
        cout << " Auth = ";
    secrets skl;
    element_init_G1(skl.sk1, pairing);
    element_init_G1(skl.sk2, pairing);
    element_init_G1(skl.sk3, pairing);
    int count = 0;
    for (int &attr : cond)
    {
        if (count == 0)
        {
            if (out)
                cout << attr;
            secrets *curr_sk = sks[attr]; // deep copy
            element_set(skl.sk1, curr_sk->sk1);
            element_set(skl.sk2, curr_sk->sk2);
            element_set(skl.sk3, curr_sk->sk3);
        }
        else
        {
            if (out)
                cout << " and " << attr;
            secrets *curr_sk = sks[attr]; // shallow copy
            element_mul(skl.sk1, skl.sk1, curr_sk->sk1);
            element_mul(skl.sk2, skl.sk2, curr_sk->sk2);
            element_mul(skl.sk3, skl.sk3, curr_sk->sk3);
        }
        count++;
    }
    if (out)
        cout << endl;
    element_t pair1, pair2, pair3, m;
    element_init_GT(pair1, pairing);
    element_init_GT(pair2, pairing);
    element_init_GT(pair3, pairing);
    element_init_GT(m, pairing);

    element_pairing(pair1, ciphertext.c3, skl.sk3);
    element_pairing(pair2, ciphertext.c4, skl.sk2);
    element_pairing(pair3, ciphertext.c2, skl.sk1);
    element_mul(m, ciphertext.c1, pair1);
    element_mul(m, m, pair2);
    element_div(m, m, pair3);
    // if(out) element_printf(" dec_m = %B\n", m);
    element_to_mpz(message_mpz, m);
    valueToMessage(dec_m, message_mpz);

    element_clear(skl.sk1);
    element_clear(skl.sk2);
    element_clear(skl.sk3);
    element_clear(pair1);
    element_clear(pair2);
    element_clear(pair3);
    element_clear(m);

    end = chrono::system_clock::now();
    double time = static_cast<double>(chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0);
    if (out)
    {
        vector<char[2048]> tmp_v(1);
        strcpy(tmp_v[0], dec_m);
        if (Check_Message(tmp_v, dec_m))
            printf(" Message = %s\n", dec_m);
        else
            printf(" Failed to decrypt message\n");
        printf("Finished Decryption (time %lf[ms])\n", time);
    }
}

void P_Decrypt(vector<cipher> &ciphertext, char *dec_m, vector<int> &cond, map<int, secrets *> &sks)
{
    cout << "\nStarting Decryption (Parallel)" << endl;
    cout << " Auth = " << cond[0];
    for (int i = 1; i < cond.size(); i++)
        cout << " and " << cond[i];
    cout << endl;

    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();
    const int thread_n = ciphertext.size();
    vector<char[2048]> dec_m_tmp(thread_n);
    vector<thread> threads;
    bool res = false;

    for (int i = 0; i < thread_n; ++i)
        threads.push_back(thread(Decrypt, ref(ciphertext[i]), dec_m_tmp[i], ref(cond), ref(sks), false));
    cout << " " << threads.size() << " threads started." << endl;
    for (thread &t : threads)
        t.join();
    if (Check_Message(dec_m_tmp, dec_m))
        cout << " Message = " << dec_m << endl;
    else
        cout << " Failed to decrypt message." << endl;
    cout << "Finished Decryption ";
    end = chrono::system_clock::now();
    double time = static_cast<double>(chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0);
    printf("(time %lf[ms])\n", time);
}

int main(int argc, char **argv)
{
    // Setup (attribute = 0 ~ 13)
    // 0 | 2 | 6 | 10
    // 1 | 3 | 7 | 11
    //   | 4 | 8 | 12
    //   | 5 | 9 | 13
    Setup(argc, argv);
    map<int, auth> aas;
    for (int attr = 0; attr <= 13; attr++)
        aas.insert(make_pair(attr, CreateAuthority(attr)));

    // Encrypt (condition = 0 and (2 or 3) and (6 or 7) and (10 or 11))
    vector<int> w1 = {0, 2, 6, 10}; // 0 and 2 and 6 and 10
    vector<int> w2 = {0, 2, 6, 11}; // 0 and 2 and 6 and 11
    vector<int> w3 = {0, 2, 7, 10}; // 0 and 2 and 7 and 10
    vector<int> w4 = {0, 2, 7, 11}; // 0 and 2 and 7 and 11
    vector<int> w5 = {0, 3, 6, 10}; // 0 and 3 and 6 and 10
    vector<int> w6 = {0, 3, 6, 11}; // 0 and 3 and 6 and 11
    vector<int> w7 = {0, 3, 7, 10}; // 0 and 3 and 7 and 10
    vector<int> w8 = {0, 3, 7, 11}; // 0 and 3 and 7 and 11
    vector<vector<int>> conds = {w1, w2, w3, w4, w5, w6, w7, w8};
    vector<cipher> ct(conds.size()); // ciphertext
    map<int, element_t *> pks;       // PK dict
    string message("hello world!!"); // message
    string head("message:");
    char raw_message[2048];
    strcpy(raw_message, (head + message).c_str());
    for (vector<int> &cond : conds)
        for (int &attr : cond)
            if (pks.find(attr) == pks.end())
                pks.insert(make_pair(attr, &aas[attr].pk));
    Encrypt(ct, raw_message, conds, pks);

    // Decrypt user1(auths = 0 and 2 and 6 and 10)
    vector<int> user1_a = {0, 2, 6, 10};
    char dec_message1[2048];
    map<int, secrets *> sks;
    for (int &attr : user1_a)
        sks.insert(make_pair(attr, &aas[attr].sk));
    P_Decrypt(ct, dec_message1, user1_a, sks);
    Decrypt(ct[0], dec_message1, user1_a, sks, true);

    // Decrypt user2(auths = 1 and 2 and 6 and 10)
    vector<int> user2_a = {1, 2, 6, 10};
    char dec_message2[2048];
    sks.clear();
    for (int &attr : user2_a)
        sks.insert(make_pair(attr, &aas[attr].sk));
    P_Decrypt(ct, dec_message2, user2_a, sks);

    // Decrypt user3(auths = 0 and 3 and 7 and 11)
    vector<int> user3_a = {0, 3, 7, 11};
    char dec_message3[2048];
    sks.clear();
    for (int &attr : user3_a)
        sks.insert(make_pair(attr, &aas[attr].sk));
    P_Decrypt(ct, dec_message3, user3_a, sks);

    // Decrypt user4(auths = 0 and 3 and 7 and 12)
    vector<int> user4_a = {0, 3, 7, 12};
    char dec_message4[2048];
    sks.clear();
    for (int &attr : user4_a)
        sks.insert(make_pair(attr, &aas[attr].sk));
    P_Decrypt(ct, dec_message4, user4_a, sks);
}
