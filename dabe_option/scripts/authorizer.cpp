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
#define MAX_SIZE 1000

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
    int attr_c;
    element_t hi;
    element_t pk;
    secrets sk;
    vector<int> users;
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
    element_t c1;
    element_t c2;
    element_t c3;
    element_t c4;
    vector<vector<element_t *>> c5;
} cipher2;

// Public information
pairing_t pairing;
element_t P, Q, R, S;
map<int, int> attr_dict;
map<int, element_t *> T;
int non_opt_w = 1;
int opt_w = 3;
vector<int> non_opt_h{2};
vector<int> opt_h{10, 5, 10};

void print_2d(vector<vector<int>> &list, int split)
{
    int int_w = 2;
    string sp(int_w, ' ');
    string ind = "          ";
    int h = list.size();
    int w = list[0].size();
    cout << ind;
    for (int i = 0; i < w; i++)
    {
        if (i == split)
            cout << " |";
        cout << " [" << setw(int_w) << i << " ]";
    }
    cout << endl;
    for (int i = 0; i < h; i++)
    {
        cout << ind;
        for (int j = 0; j < w; j++)
        {
            if (j == split)
                cout << "||";
            if (list[i][j] == -1)
                cout << "|  " << sp << " ";
            else
                cout << "|  " << setw(int_w) << list[i][j] << " ";
        }
        cout << "|" << endl;
    }
}

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

    int attr_w, attr_h;
    attr_w = non_opt_w + opt_w;
    attr_h = max(*max_element(non_opt_h.begin(), non_opt_h.end()), *max_element(opt_h.begin(), opt_h.end()));
    vector<vector<int>> attr_l(attr_h, vector<int>(attr_w, -1));
    int attr_c = 0;
    for (int i = 0; i < non_opt_w; i++)
    {
        for (int j = 0; j < non_opt_h[i]; j++)
        {
            attr_l[j][i] = attr_c;
            attr_dict.insert(make_pair(attr_c, i));
            attr_c++;
        }
    }
    element_t Ts[MAX_SIZE];
    for (int i = 0; i < opt_w; i++)
    {
        for (int j = 0; j < opt_h[i]; j++)
        {
            element_t tmp_z;
            element_init_G1(Ts[attr_c], pairing);
            element_init_Zr(tmp_z, pairing);
            element_random(tmp_z);
            element_pow_zn(Ts[attr_c], P, tmp_z);
            element_clear(tmp_z);
            //element_printf("T %B\n", Ts[attr_c]);
            T.insert(make_pair(attr_c, &Ts[attr_c]));
            int opt_attr = i + non_opt_w;
            attr_l[j][opt_attr] = attr_c;
            attr_dict.insert(make_pair(attr_c, opt_attr));
            attr_c++;
        }
    }
    cout << " Non-option | option attribute:" << endl;
    print_2d(attr_l, non_opt_w);
    end = chrono::system_clock::now();
    double time = static_cast<double>(chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0);
    printf("Finished Setup (time %lf[ms])\n", time);
}

auth CreateAuthority(int i)
{
    auth aa;
    aa.attr_c = i;
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
    if (T.find(i) != T.end())
    {
        element_t tmp2, tmp3;
        element_init_G1(tmp2, pairing);
        element_init_G1(tmp3, pairing);
        element_set(tmp3, *T[i]);
        element_pow_zn(tmp2, tmp3, y);
        element_mul(tmp, tmp2, tmp);
        element_clear(tmp2);
        element_clear(tmp3);
    }
    element_pow_zn(tmp, tmp, aa.hi);
    element_pow_zn(aa.sk.sk1, S, x);
    element_mul(aa.sk.sk1, tmp, aa.sk.sk1);
    // sk2
    element_pow_zn(aa.sk.sk2, P, x);
    // sk3
    element_pow_zn(aa.sk.sk3, P, y);
    element_pow_zn(aa.sk.sk3, aa.sk.sk3, aa.hi);

    element_clear(x);
    element_clear(y);
    element_clear(tmp);

    return aa;
}

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
void Encrypt(cipher2 &ciphertext, char *raw_m, vector<int> &w, map<int, element_t *> &pks)
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
    element_t tmp_rTs[MAX_SIZE];
    vector<vector<element_t *>> tmp_c5(opt_w);
    element_t pkl, r;
    element_init_G1(pkl, pairing);
    element_init_Zr(r, pairing);
    element_random(r);
    int count = 0;
    string non_opt_out = "", opt_out = "";
    for (int &attr_c : w)
    {
        element_t *curr_pk = pks[attr_c];
        if (count == 0)
            element_set(pkl, *curr_pk);
        else
            element_mul(pkl, pkl, *curr_pk);
        if (attr_dict[attr_c] < non_opt_w)
        {
        }
        else
        {
            int opt_attr = attr_dict[attr_c] - non_opt_w;
            cout << attr_c << " " << opt_attr << endl;
            element_init_G1(tmp_rTs[attr_c], pairing);
            element_set(tmp_rTs[attr_c], *T[attr_c]);
            element_pow_zn(tmp_rTs[attr_c], tmp_rTs[attr_c], r);
            tmp_c5[opt_attr].push_back(&tmp_rTs[attr_c]);
        }
        count++;
    }

    cout << " W = " << endl;
    element_t pair1, tmp;
    element_init_GT(ciphertext.c1, pairing);
    element_init_G1(ciphertext.c2, pairing);
    element_init_G1(ciphertext.c3, pairing);
    element_init_G1(ciphertext.c4, pairing);
    element_init_G1(tmp, pairing);
    element_init_GT(pair1, pairing);

    element_pairing(pair1, Q, pkl);
    element_pow_zn(pair1, pair1, r);
    element_mul(ciphertext.c1, message_ele, pair1);
    element_pow_zn(ciphertext.c2, P, r);
    element_pow_zn(ciphertext.c3, R, r);
    element_pow_zn(ciphertext.c4, S, r);
    ciphertext.c5 = tmp_c5;
    element_clear(pkl);
    element_clear(pair1);
    element_clear(r);
    element_clear(tmp);

    end = chrono::system_clock::now();
    double time = static_cast<double>(chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0);
    printf("\nFinished Encyption (time %lf[ms])\n", time);
}

void Decrypt(cipher2 &ciphertext, char *dec_m, vector<int> &cond, vector<int> &cond_c, map<int, secrets *> &sks, bool out)
{
    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();
    if (out)
        cout << "\nStarting Decryption" << endl;
    if (cond_c.size() != opt_w)
    {
        cout << "\nERROR(decrypt): invalid opt_attrs size" << endl;
        return;
    }
    char message_dec[2048];
    mpz_t message_mpz;
    mpz_init(message_mpz);
    if (out)
        cout << " Auth = ";
    vector<int> opt_attrs;
    secrets skl;
    element_init_G1(skl.sk1, pairing);
    element_init_G1(skl.sk2, pairing);
    element_init_G1(skl.sk3, pairing);
    int count = 0;
    for (int &attr_c : cond)
    {
        if (non_opt_w <= attr_dict[attr_c])
            opt_attrs.push_back(attr_c);
        if (count == 0)
        {
            if (out)
                cout << attr_c;
            secrets *curr_sk = sks[attr_c]; // deep copy
            element_set(skl.sk1, curr_sk->sk1);
            element_set(skl.sk2, curr_sk->sk2);
            element_set(skl.sk3, curr_sk->sk3);
        }
        else
        {
            if (out)
                cout << " and " << attr_c;
            secrets *curr_sk = sks[attr_c]; // shallow copy
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
    for (int &attr_c : opt_attrs)
    {
        int x, y;
        x = attr_dict[attr_c] - non_opt_w;
        y = cond_c[x];
        cout << x << " " << y << endl;
        element_printf("%B\n", *ciphertext.c5[x][y]);
        element_t pair_tmp;
        element_init_GT(pair_tmp, pairing);
        element_pairing(pair_tmp, *ciphertext.c5[x][y], sks[attr_c]->sk3);
        element_mul(m, m, pair_tmp);
        element_clear(pair_tmp);
    }
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
        printf(" Message = %s\n", dec_m);
        if (Check_Message(tmp_v, dec_m))
            printf(" Message = %s\n", dec_m);
        else
            printf(" Failed to decrypt message\n");
        printf("Finished Decryption (time %lf[ms])\n", time);
    }
}

/*
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
*/

int main(int argc, char **argv)
{
    // Setup (attribute = 0 ~ 13)
    // 0 | 2 | 6 | 10
    // 1 | 3 | 7 | 11
    //   | 4 | 8 | 12
    //   | 5 | 9 | 13
    Setup(argc, argv);
    map<int, auth> aas;
    for (int attr_c = 0; attr_c <= 26; attr_c++)
        aas.insert(make_pair(attr_c, CreateAuthority(attr_c)));

    cipher2 ct;
    vector<int> w{0, 2, 12, 17};
    map<int, element_t *> pks;
    string message("hello world!!"); // message
    string head("message:");
    char raw_message[2048];
    strcpy(raw_message, (head + message).c_str());
    for (int &attr_c : w)
        pks.insert(make_pair(attr_c, &aas[attr_c].pk));
    for (int attr_c = 2; attr_c <= 26; attr_c++)
        element_printf("T %d %p %B\n", attr_c, T[attr_c], *T[attr_c]);
    Encrypt(ct, raw_message, w, pks);
    
    for(int i = 0; i < ct.c5.size(); i++)
        for(int j = 0; j < ct.c5[i].size(); j++)
            element_printf("ct (%d, %d) : %B\n", i, j, *ct.c5[i][j]);

    vector<int> user{0, 2, 12, 17};
    vector<int> cc{0,0,0};
    char dec_message1[2048];
    map<int, secrets *> sks;
    for (int &attr_c : user)
        sks.insert(make_pair(attr_c, &aas[attr_c].sk));
    Decrypt(ct, dec_message1, user, cc, sks, true);
    /*
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
        for (int &attr_c : cond)
            if (pks.find(attr_c) == pks.end())
                pks.insert(make_pair(attr_c, &aas[attr_c].pk));
    Encrypt(ct, raw_message, conds, pks);

    // Decrypt user1(auths = 0 and 2 and 7 and 10)
    vector<int> user1_a = {0, 2, 7, 10};
    char dec_message1[2048];
    map<int, secrets *> sks;
    for (int &attr_c : user1_a)
        sks.insert(make_pair(attr_c, &aas[attr_c].sk));
    P_Decrypt(ct, dec_message1, user1_a, sks);

    // Decrypt user2(auths = 1 and 2 and 7 and 10)
    vector<int> user2_a = {1, 2, 7, 10};
    char dec_message2[2048];
    sks.clear();
    for (int &attr_c : user2_a)
        sks.insert(make_pair(attr_c, &aas[attr_c].sk));
    P_Decrypt(ct, dec_message2, user2_a, sks);

    // Decrypt user3(auths = 0 and 3 and 7 and 11)
    vector<int> user3_a = {0, 3, 7, 11};
    char dec_message3[2048];
    sks.clear();
    for (int &attr_c : user3_a)
        sks.insert(make_pair(attr_c, &aas[attr_c].sk));
    P_Decrypt(ct, dec_message3, user3_a, sks);
    Decrypt(ct[7], dec_message3, user3_a, sks, true);

    // Decrypt user4(auths = 0 and 3 and 7 and 12)
    vector<int> user4_a = {0, 3, 7, 12};
    char dec_message4[2048];
    sks.clear();
    for (int &attr_c : user4_a)
        sks.insert(make_pair(attr_c, &aas[attr_c].sk));
    P_Decrypt(ct, dec_message4, user4_a, sks);
    Decrypt(ct[7], dec_message4, user4_a, sks, true);
    */
}
