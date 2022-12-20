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
    int attr;
    element_t hi;
    element_t pk;
    map<int, secrets *> sks;
    vector<int> users;
} auth2;

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
map<int, element_t *> T;    // key = attribute child, value = T (G1)
map<int, int> attr_dict;    // key = attribute child, value = attribute
vector<int> non_opt_h{2};   // non-option attributes size (only one value)
vector<int> opt_h{4, 4, 4}; // option attributes size

void print_2d(vector<vector<int>> &list, int split, int rows)
{
    int int_w = 2;
    string sp(int_w, ' ');
    string ind = "          ";
    int h = rows;
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
            if (list[i][j] == -1 || j >= list[i].size())
                cout << "|  " << sp << " ";
            else
                cout << "|  " << setw(int_w) << list[i][j] << " ";
        }
        cout << "|" << endl;
    }
}

int Get_cipher_size(cipher &ciphertext)
{
    int tmp, sum = 0;
    element_t tmp1, tmp2;
    element_init_GT(tmp1, pairing);
    element_init_G1(tmp2, pairing);
    element_set(tmp1, ciphertext.c1);
    tmp = pairing_length_in_bytes_GT(pairing);
    ;
    sum += tmp;
    // unsigned char *data1 = (unsigned char*)pbc_malloc(tmp);
    // element_to_bytes_compressed(data1, tmp1);

    element_set(tmp2, ciphertext.c2);
    tmp = pairing_length_in_bytes_compressed_G1(pairing);
    ;
    sum += tmp;
    // unsigned char *data2 = (unsigned char*)pbc_malloc(tmp);
    // element_to_bytes_compressed(data2, tmp2);

    element_set(tmp2, ciphertext.c3);
    tmp = pairing_length_in_bytes_compressed_G1(pairing);
    ;
    sum += tmp;
    // unsigned char *data3 = (unsigned char*)pbc_malloc(tmp);
    // element_to_bytes_compressed(data3, tmp2);

    element_set(tmp2, ciphertext.c4);
    tmp = pairing_length_in_bytes_compressed_G1(pairing);
    ;
    sum += tmp;
    // unsigned char *data4 = (unsigned char*)pbc_malloc(tmp);
    // element_to_bytes_compressed(data4, tmp2);

    element_clear(tmp1);
    element_clear(tmp2);

    return sum;
}

int Get_cipher2_size(cipher2 &ciphertext)
{
    int tmp, sum = 0;
    element_t tmp1, tmp2;
    element_init_GT(tmp1, pairing);
    element_init_G1(tmp2, pairing);
    element_set(tmp1, ciphertext.c1);
    tmp = pairing_length_in_bytes_GT(pairing);
    sum += tmp;
    unsigned char *data1 = (unsigned char*)pbc_malloc(tmp);
    element_to_bytes(data1, tmp1);

    element_set(tmp2, ciphertext.c2);
    tmp = pairing_length_in_bytes_compressed_G1(pairing);
    sum += tmp;
    unsigned char *data2 = (unsigned char*)pbc_malloc(tmp);
    element_to_bytes_compressed(data2, tmp2);

    element_set(tmp2, ciphertext.c3);
    tmp = pairing_length_in_bytes_compressed_G1(pairing);
    sum += tmp;
    unsigned char *data3 = (unsigned char*)pbc_malloc(tmp);
    element_to_bytes_compressed(data3, tmp2);

    element_set(tmp2, ciphertext.c4);
    tmp = pairing_length_in_bytes_compressed_G1(pairing);
    sum += tmp;
    unsigned char *data4 = (unsigned char*)pbc_malloc(tmp);
    element_to_bytes_compressed(data4, tmp2);
    vector<unsigned char *> data5s;
    for (vector<element_t *> &c5_l : ciphertext.c5)
        for (element_t *c5 : c5_l)
        {
            tmp = pairing_length_in_bytes_compressed_G1(pairing);
            ;
            sum += tmp;
            unsigned char *data5 = (unsigned char*)pbc_malloc(tmp);
            element_to_bytes_compressed(data5, tmp2);
            data5s.push_back(data5);
        }

    element_clear(tmp1);
    element_clear(tmp2);

    return sum;
}

bool Check_Message(vector<string> &dec_m_tmp, char *dec_m)
{
    bool res = false;
    for (string &m : dec_m_tmp)
    {
        try
        {
            if (m.substr(0, 8) == "message:")
            {
                strcpy(dec_m, m.substr(8).c_str());
                res = true;
            }
        }
        catch (...)
        {
        }
    }
    return res;
}

void mul_pair_dec(cipher2 &ciphertext, map<int, secrets *> &sks, element_t &m, vector<int> &cond_c, vector<int> &opt_attrs)
{
    int non_opt_w = non_opt_h.size();
    for (int &attr_c : opt_attrs)
    {
        int x, y;
        x = attr_dict[attr_c] - non_opt_w;
        y = cond_c[x];
        // cout << attr_c << " " << x << " " << y << endl;
        // element_printf("%B\n", *ciphertext.c5[x][y]);
        element_t pair_tmp;
        element_init_GT(pair_tmp, pairing);
        element_pairing(pair_tmp, *ciphertext.c5[x][y], sks[attr_c]->sk3);
        element_mul(m, m, pair_tmp);
        element_clear(pair_tmp);
    }
}

void Setup(int argc, char **argv)
{
    int non_opt_w = non_opt_h.size();
    int opt_w = opt_h.size();
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
            // cout << attr_c << endl;
            T.insert(make_pair(attr_c, &Ts[attr_c]));
            int opt_attr = i + non_opt_w;
            attr_l[j][opt_attr] = attr_c;
            attr_dict.insert(make_pair(attr_c, opt_attr));
            attr_c++;
        }
    }
    cout << " Non-option(0) | Option attributes :" << endl;
    print_2d(attr_l, non_opt_w, attr_l[0].size());
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

    aa.attr_c = i;

    element_clear(x);
    element_clear(y);
    element_clear(tmp);

    return aa;
}

auth2 CreateAuthority2(int i)
{
    int non_opt_w = non_opt_h.size();
    int opt_w = opt_h.size();
    auth2 aa;
    aa.attr = i;
    vector<int> attr_cs;
    for (auto it = attr_dict.begin(); it != attr_dict.end(); it++)
        if (it->second == i)
            attr_cs.push_back(it->first);

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

    secrets tmp_sk[MAX_SIZE];
    for (int &attr_c : attr_cs)
    {
        element_t tmp;
        element_init_G1(tmp_sk[attr_c].sk1, pairing);
        element_init_G1(tmp_sk[attr_c].sk2, pairing);
        element_init_G1(tmp_sk[attr_c].sk3, pairing);
        element_init_G1(tmp, pairing);
        // sk1
        element_pow_zn(tmp, R, y);
        element_mul(tmp, Q, tmp);
        if (T.find(attr_c) != T.end())
        {
            element_t tmp2, tmp3;
            element_init_G1(tmp2, pairing);
            element_init_G1(tmp3, pairing);
            element_set(tmp3, *T[attr_c]);
            element_pow_zn(tmp2, tmp3, y);
            element_mul(tmp, tmp2, tmp);
            element_clear(tmp2);
            element_clear(tmp3);
        }
        else
        {
            element_clear(tmp);
            break;
        }
        element_pow_zn(tmp, tmp, aa.hi);
        element_pow_zn(tmp_sk[attr_c].sk1, S, x);
        element_mul(tmp_sk[attr_c].sk1, tmp, tmp_sk[attr_c].sk1);
        // sk2
        element_pow_zn(tmp_sk[attr_c].sk2, P, x);
        // sk3
        element_pow_zn(tmp_sk[attr_c].sk3, P, y);
        element_pow_zn(tmp_sk[attr_c].sk3, tmp_sk[attr_c].sk3, aa.hi);
        aa.sks.insert(make_pair(attr_c, &tmp_sk[attr_c]));
        // cout << attr_c << " " << &tmp_sk[attr_c] << endl;
        element_clear(tmp);
    }
    element_clear(x);
    element_clear(y);

    return aa;
}

void Encrypt(cipher2 &ciphertext, char *raw_m, vector<int> &w, vector<element_t *> &pks)
{
    int non_opt_w = non_opt_h.size();
    int opt_w = opt_h.size();
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
    int attr_w, attr_h;
    attr_w = non_opt_w + opt_w;
    attr_h = max(*max_element(non_opt_h.begin(), non_opt_h.end()), *max_element(opt_h.begin(), opt_h.end()));
    vector<vector<int>> attr_l(attr_h, vector<int>(attr_w, -1));
    vector<int> max(attr_w, 0);
    for (element_t *curr_pk : pks)
    {
        if (count == 0)
            element_set(pkl, *curr_pk);
        else
            element_mul(pkl, pkl, *curr_pk);
        count++;
    }
    for (int &attr_c : w)
    {
        // cout << max[attr_dict[attr_c]] << " " << attr_dict[attr_c] << " " << attr_c << endl;
        attr_l[max[attr_dict[attr_c]]][attr_dict[attr_c]] = attr_c;
        max[attr_dict[attr_c]]++;
        if (attr_dict[attr_c] != 0)
        {
            int opt_attr = attr_dict[attr_c] - non_opt_w;
            // cout << attr_c << " " << opt_attr;
            element_init_G1(tmp_rTs[attr_c], pairing);
            element_set(tmp_rTs[attr_c], *T[attr_c]);
            element_pow_zn(tmp_rTs[attr_c], tmp_rTs[attr_c], r);
            // element_printf(" %B\n", tmp_rTs[attr_c]);
            tmp_c5[opt_attr].push_back(&tmp_rTs[attr_c]);
        }
    }

    cout << " Non-option(0) | Option attributes (W):" << endl;
    print_2d(attr_l, non_opt_w, *max_element(max.begin(), max.end()));
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
    printf("Finished Encyption (time %lf[ms])\n", time);
}

void Decrypt(cipher2 &ciphertext, char *dec_m, vector<int> &cond, vector<int> &cond_c, map<int, secrets *> &sks, bool out)
{
    int non_opt_w = non_opt_h.size();
    int opt_w = opt_h.size();
    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();
    if (out)
        cout << "\nStarting Decryption" << endl;
    if (cond_c.size() != opt_w)
    {
        cout << "\nERROR(decrypt): invalid opt_attrs size" << endl;
        return;
    }
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
    mul_pair_dec(ciphertext, sks, m, cond_c, opt_attrs);

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
        vector<string> tmp_v{string(dec_m)};
        if (Check_Message(tmp_v, dec_m))
            printf(" Message = %s\n", dec_m);
        else
            printf(" Failed to decrypt message\n");
        printf("Finished Decryption (time %lf[ms])\n", time);
    }
}

void P_Decrypt(cipher2 &ciphertext, char *dec_m, vector<int> &cond, map<int, secrets *> &sks)
{
    int non_opt_w = non_opt_h.size();
    int opt_w = opt_h.size();
    cout << "\nStarting Decryption (Parallel)" << endl;
    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();    

    mpz_t message_mpz;
    mpz_init(message_mpz);
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
            cout << attr_c;
            secrets *curr_sk = sks[attr_c]; // deep copy
            element_set(skl.sk1, curr_sk->sk1);
            element_set(skl.sk2, curr_sk->sk2);
            element_set(skl.sk3, curr_sk->sk3);
        }
        else
        {
            cout << " and " << attr_c;
            secrets *curr_sk = sks[attr_c]; // shallow copy
            element_mul(skl.sk1, skl.sk1, curr_sk->sk1);
            element_mul(skl.sk2, skl.sk2, curr_sk->sk2);
            element_mul(skl.sk3, skl.sk3, curr_sk->sk3);
        }
        count++;
    }
    cout << endl;

    int len = ciphertext.c5.size();
    vector<int> ps(len, 0);
    vector<vector<int>> c_indexs;
    bool stop = true;
    while (stop)
    {
        vector<int> c_index;
        for (int i = 0; i < len; i++)
            c_index.push_back(ps[i]);
        c_indexs.push_back(c_index);
        for (int i = 0; i < len; i++)
        {
            if (++ps[i] < ciphertext.c5[i].size())
                break;
            else
            {
                if (i >= len - 1)
                    stop = false;
                ps[i] = 0;
            }
        }
    }
    int thread_n = c_indexs.size();
    vector<thread> threads;
    vector<string> dec_m_tmp(thread_n);
    element_t m_tmp[thread_n];
    for (int i = 0; i < thread_n; i++)
    {
        element_init_GT(m_tmp[i], pairing);
        element_set1(m_tmp[i]);
        threads.push_back(thread(mul_pair_dec, ref(ciphertext), ref(sks), ref(m_tmp[i]), ref(c_indexs[i]), ref(opt_attrs)));
    }
    cout << " " << threads.size() << " threads started." << endl;

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

    for (int i = 0; i < thread_n; i++)
    {
        threads[i].join();
        element_mul(m_tmp[i], m_tmp[i], m);
        element_to_mpz(message_mpz, m_tmp[i]);
        valueToMessage(dec_m, message_mpz);
        element_clear(m_tmp[i]);
        dec_m_tmp[i] = string(dec_m);
    }
    element_clear(skl.sk1);
    element_clear(skl.sk2);
    element_clear(skl.sk3);
    element_clear(pair1);
    element_clear(pair2);
    element_clear(pair3);
    element_clear(m);
    if (Check_Message(dec_m_tmp, dec_m))
        cout << " Message = " << dec_m << endl;
    else
        cout << " Failed to decrypt message." << endl;
    cout << "Finished Decryption ";
    end = chrono::system_clock::now();
    double time = static_cast<double>(chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0);
    printf("(time %lf[ms])\n", time);
}

void Get_SKs(vector<int> &user, map<int, auth> &non_opt_aas, map<int, auth2> &opt_aas, map<int, secrets *> &sks)
{
    for (int &attr_c : user)
    {
        int attr = attr_dict[attr_c];
        if (attr == 0)
            sks.insert(make_pair(attr_c, &non_opt_aas[attr_c].sk));
        else
            sks.insert(make_pair(attr_c, opt_aas[attr].sks[attr_c]));
    }
}

int main(int argc, char **argv)
{
    // Setup
    Setup(argc, argv);
    map<int, auth> non_opt_aas;
    map<int, auth2> opt_aas;
    for (auto it = attr_dict.begin(); it != attr_dict.end(); it++)
    {
        int attr = it->second;
        int attr_c = it->first;
        if (attr == 0)
        {
            if (non_opt_aas.find(attr_c) == non_opt_aas.end())
                non_opt_aas.insert(make_pair(attr_c, CreateAuthority(attr_c)));
        }
        else
        {
            if (opt_aas.find(attr) == opt_aas.end())
                opt_aas.insert(make_pair(attr, CreateAuthority2(attr)));
        }
    }

    // Encrypt
    cipher2 ct;                           // ciphertext
    vector<element_t *> pks;              // PKs
    vector<int> w{0, 2, 3, 6, 7, 10, 11}; // condition
    string message("hello world!!");      // message
    string head("message:");
    char raw_message[2048];
    strcpy(raw_message, (head + message).c_str());
    for (int &attr_c : w)
        if (attr_dict[attr_c] == 0)
            pks.push_back(&non_opt_aas[attr_c].pk);
    for (int i = non_opt_h.size(); i < opt_h.size() + non_opt_h.size(); i++)
        pks.push_back(&opt_aas[i].pk);
    Encrypt(ct, raw_message, w, pks);

    // Decrypt user1
    vector<int> user1{0, 2, 6, 10};
    char dec_message1[2048];
    map<int, secrets *> sks;
    Get_SKs(user1, non_opt_aas, opt_aas, sks);
    P_Decrypt(ct, dec_message1, user1, sks);
    vector<int> cc{0, 0, 0};
    Decrypt(ct, dec_message1, user1, cc, sks, true);

    // Decrypt user2
    vector<int> user2{1, 2, 6, 10};
    char dec_message2[2048];
    sks.clear();
    Get_SKs(user2, non_opt_aas, opt_aas, sks);
    P_Decrypt(ct, dec_message2, user2, sks);

    // Decrypt user3
    vector<int> user3{0, 3, 7, 11};
    char dec_message3[2048];
    sks.clear();
    Get_SKs(user3, non_opt_aas, opt_aas, sks);
    P_Decrypt(ct, dec_message3, user3, sks);

    // Decrypt user4
    vector<int> user4{0, 3, 7, 12};
    char dec_message4[2048];
    sks.clear();
    Get_SKs(user4, non_opt_aas, opt_aas, sks);
    P_Decrypt(ct, dec_message4, user4, sks);

    cout << "\nCipher text size = " << Get_cipher2_size(ct) << " bytes" << endl;
    cout << "(G1 points are compressed)" << endl;
    /*
    // debug
    for (auto it = non_opt_aas.begin(); it != non_opt_aas.end(); it++)
        cout << " " << it-> first;
    cout << endl;
    for (auto it = opt_aas.begin(); it != opt_aas.end(); it++)
        cout << " " << it-> first;
    cout << endl;
    for (int attr_c = 2; attr_c <= 26; attr_c++)
        element_printf("T %d %p %B\n", attr_c, T[attr_c], *T[attr_c]);
    */
}
