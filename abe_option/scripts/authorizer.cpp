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
    map<int, element_t *> sk3;
} secrets;

typedef struct
{
    element_t w;
    vector<int> users;
} ca;

typedef struct
{
    element_t c1;
    element_t c2;
    element_t c3;
    vector<vector<element_t *>> c4;
} cipher;

// Public information
pairing_t pairing;
element_t P, ePPw;
map<int, int> attr_dict;
map<int, element_t *> Q, R;
int non_opt = 2;
vector<int> opt_h{4, 4, 4};
ca center;

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
    sum += tmp;
    unsigned char *data1 = (unsigned char *)pbc_malloc(tmp);
    element_to_bytes(data1, tmp1);

    element_set(tmp2, ciphertext.c2);
    tmp = pairing_length_in_bytes_compressed_G1(pairing);
    sum += tmp;
    unsigned char *data2 = (unsigned char *)pbc_malloc(tmp);
    element_to_bytes_compressed(data2, tmp2);

    element_set(tmp2, ciphertext.c3);
    tmp = pairing_length_in_bytes_compressed_G1(pairing);
    sum += tmp;
    unsigned char *data3 = (unsigned char *)pbc_malloc(tmp);
    element_to_bytes_compressed(data3, tmp2);

    vector<unsigned char *> data4s;
    for (vector<element_t *> &c4_l : ciphertext.c4)
        for (element_t *c4 : c4_l)
        {
            tmp = pairing_length_in_bytes_compressed_G1(pairing);
            sum += tmp;
            unsigned char *data4 = (unsigned char *)pbc_malloc(tmp);
            element_to_bytes_compressed(data4, tmp2);
            data4s.push_back(data4);
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

void mul_pair_dec(cipher &ciphertext, secrets &sk, element_t &m, vector<int> &cond_c, vector<int> &opt_attrs)
{
    for (int &attr_c : opt_attrs)
    {
        int x, y;
        x = attr_dict[attr_c] - 1;
        y = cond_c[x];
        // cout << attr_c << " " << x << " " << y << endl;
        // element_printf("%B\n", *ciphertext.c4[x][y]);
        element_t pair_tmp;
        element_init_GT(pair_tmp, pairing);
        element_pairing(pair_tmp, *ciphertext.c4[x][y], *sk.sk3[attr_c]);
        element_mul(m, m, pair_tmp);
        element_clear(pair_tmp);
    }
}

void Setup(int argc, char **argv)
{
    int opt_w = opt_h.size();
    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();
    cout << "Starting Setup" << endl;
    pbc_demo_pairing_init(pairing, argc, argv);
    if (!pairing_is_symmetric(pairing))
        pbc_die(" Error: pairing must be symmetric");

    element_init_G1(P, pairing);
    element_init_GT(ePPw, pairing);

    element_t tmp;
    element_init_G1(tmp, pairing);
    element_init_Zr(center.w, pairing);
    element_random(center.w);

    element_set(P, ((curve_data_ptr)((a_pairing_data_ptr)pairing->data)->Eq->data)->gen);
    element_pow_zn(tmp, P, center.w);
    element_pairing(ePPw, P, tmp);

    element_clear(tmp);

    int attr_c = 0;
    int attr_w, attr_h;
    attr_w = opt_w + 1;
    attr_h = max(non_opt, *max_element(opt_h.begin(), opt_h.end()));
    vector<vector<int>> attr_l(attr_h, vector<int>(attr_w, -1));
    static element_t Qs[MAX_SIZE];
    for (int i = 0; i < non_opt; i++)
    {
        element_t tmp_z;
        element_init_G1(Qs[attr_c], pairing);
        element_init_Zr(tmp_z, pairing);
        element_random(tmp_z);
        element_pow_zn(Qs[attr_c], P, tmp_z);
        element_clear(tmp_z);
        Q.insert(make_pair(attr_c, &Qs[attr_c]));
        attr_l[i][0] = attr_c;
        attr_dict.insert(make_pair(attr_c, 0));
        attr_c++;
    }

    static element_t Rs[MAX_SIZE];
    for (int i = 0; i < opt_w; i++)
    {
        for (int j = 0; j < opt_h[i]; j++)
        {
            element_t tmp_z;
            element_init_G1(Rs[attr_c], pairing);
            element_init_Zr(tmp_z, pairing);
            element_random(tmp_z);
            element_pow_zn(Rs[attr_c], P, tmp_z);
            element_clear(tmp_z);
            R.insert(make_pair(attr_c, &Rs[attr_c]));
            // element_printf("R %d %p %B\n", attr_c, R[attr_c], *R[attr_c]);
            attr_l[j][i + 1] = attr_c;
            attr_dict.insert(make_pair(attr_c, i + 1));
            attr_c++;
        }
    }
    cout << " Non-option(0) | Option attributes:" << endl;
    print_2d(attr_l, 1, attr_l[0].size());
    end = chrono::system_clock::now();
    double time = static_cast<double>(chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0);
    printf("Finished Setup (time %lf[ms])\n", time);
}

void KeyGen(vector<int> &attrs, secrets &sk)
{
    element_t u, k, s_sum;
    static element_t xP[MAX_SIZE];
    element_init_Zr(u, pairing);
    element_init_Zr(k, pairing);
    element_init_Zr(s_sum, pairing);
    element_random(u);
    element_random(k);
    element_set0(s_sum);

    // sk2
    element_init_G1(sk.sk2, pairing);
    element_pow_zn(sk.sk2, P, u);

    // sk1 & sk3
    element_t tmp;
    element_init_G1(tmp, pairing);
    element_set1(tmp);
    for (int &attr_c : attrs)
    {
        int attr = attr_dict[attr_c];
        // cout << attr_c << " " << attr << endl;
        if (attr == 0)
        {
            element_t tmp2;
            element_init_G1(tmp2, pairing);
            element_pow_zn(tmp2, *Q[attr_c], u);
            element_mul(tmp, tmp, tmp2);

            element_clear(tmp2);
        }
        else
        {
            element_t x, s;
            element_init_Zr(x, pairing);
            element_init_Zr(s, pairing);
            element_random(x);
            element_random(s);
            element_add(s_sum, s_sum, s);

            element_t tmp2, tmp3;
            element_init_G1(tmp2, pairing);
            element_init_G1(tmp3, pairing);
            element_pow_zn(tmp2, P, s);
            element_pow_zn(tmp3, *R[attr_c], x);
            element_mul(tmp, tmp, tmp2);
            element_mul(tmp, tmp, tmp3);
            

            element_init_G1(xP[attr_c], pairing);
            element_pow_zn(xP[attr_c], P, x);
            sk.sk3.insert(make_pair(attr_c, &xP[attr_c]));

            element_clear(x);
            element_clear(s);
            element_clear(tmp2);
            element_clear(tmp3);
        }
    }
    element_t tmp2, tmp3, negk;
    element_init_G1(tmp2, pairing);
    element_init_G1(tmp3, pairing);
    element_init_Zr(negk, pairing);
    element_sub(s_sum, center.w, s_sum);
    element_add(s_sum, s_sum, k);
    element_pow_zn(tmp2, P, s_sum);

    element_neg(negk, k);
    element_pow_zn(tmp3, P, negk);

    element_mul(tmp, tmp, tmp2);
    element_mul(tmp, tmp, tmp3);

    element_init_G1(sk.sk1, pairing);
    element_set(sk.sk1, tmp);

    element_clear(u);
    element_clear(k);
    element_clear(tmp);
    element_clear(tmp2);
    element_clear(tmp3);
}

void Encrypt(cipher &ciphertext, string &raw_message, vector<int> &attrs)
{
    int opt_w = opt_h.size();
    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();
    cout << "\nStarting Encyption\n Message = " << raw_message << endl;
    string head("message:");
    char raw_m[2048];
    strcpy(raw_m, (head + raw_message).c_str());
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

    element_t r;
    element_init_Zr(r, pairing);
    element_random(r);

    element_init_GT(ciphertext.c1, pairing);
    element_pow_zn(ciphertext.c1, ePPw, r);
    element_mul(ciphertext.c1, ciphertext.c1, message_ele);

    element_init_G1(ciphertext.c2, pairing);
    element_pow_zn(ciphertext.c2, P, r);

    element_init_G1(ciphertext.c3, pairing);
    vector<vector<element_t *>> c4(opt_w);
    ciphertext.c4 = c4;
    static element_t rR[MAX_SIZE];
    int attr_w, attr_h;
    attr_w = opt_w + 1;
    attr_h = max(non_opt, *max_element(opt_h.begin(), opt_h.end()));
    vector<vector<int>> attr_l(attr_h, vector<int>(attr_w, -1));
    vector<int> attr_i(attr_w, 0);
    for (int &attr_c : attrs)
    {
        int attr = attr_dict[attr_c];
        if (attr == 0)
        {
            element_t tmp;
            element_init_G1(tmp, pairing);
            element_pow_zn(tmp, *Q[attr_c], r);
            element_mul(ciphertext.c3, ciphertext.c3, tmp);
        }
        else
        {
            element_init_G1(rR[attr_c], pairing);
            element_pow_zn(rR[attr_c], *R[attr_c], r);
            ciphertext.c4[attr - 1].push_back(&rR[attr_c]);
        }
        attr_l[attr_i[attr]][attr] = attr_c;
        attr_i[attr]++;
    }
    cout << " Non-option(0) | Option attributes (W):" << endl;
    print_2d(attr_l, 1, *max_element(attr_i.begin(), attr_i.end()));
    end = chrono::system_clock::now();
    double time = static_cast<double>(chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0);
    printf("Finished Encyption (time %lf[ms])\n", time);
}

void Decrypt(cipher &ciphertext, vector<int> &cond_c, char *dec_m, vector<int> &attrs, secrets &sk)
{
    int opt_w = opt_h.size();
    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();
    cout << "\nStarting Decryption" << endl;
    if (cond_c.size() != opt_w)
    {
        cout << "\nERROR(decrypt): invalid opt_attrs size" << endl;
        return;
    }
    mpz_t message_mpz;
    mpz_init(message_mpz);
    cout << " Auth = ";
    vector<int> opt_attrs;
    int count = 0;
    for (int &attr_c : attrs)
    {
        if (1 <= attr_dict[attr_c])
            opt_attrs.push_back(attr_c);
        if (count == 0)
            cout << attr_c;
        else
            cout << " and " << attr_c;
        count++;
    }
    cout << endl;
    
    element_t pair1, pair2, pair3, m;
    element_init_GT(pair1, pairing);
    element_init_GT(pair2, pairing);
    element_init_GT(m, pairing);

    element_pairing(pair1, ciphertext.c3, sk.sk2);
    element_pairing(pair2, ciphertext.c2, sk.sk1);
    element_mul(m, ciphertext.c1, pair1);
    element_div(m, m, pair2);
    mul_pair_dec(ciphertext, sk, m, cond_c, opt_attrs);

    // if(out) element_printf(" dec_m = %B\n", m);
    element_to_mpz(message_mpz, m);
    valueToMessage(dec_m, message_mpz);

    element_clear(pair1);
    element_clear(pair2);
    element_clear(m);

    end = chrono::system_clock::now();
    double time = static_cast<double>(chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0);
    vector<string> tmp_v{string(dec_m)};
    if (Check_Message(tmp_v, dec_m))
        printf(" Message = %s\n", dec_m);
    else
        printf(" Failed to decrypt message\n");
    printf("Finished Decryption (time %lf[ms])\n", time);
}

void P_Decrypt(cipher &ciphertext, char *dec_m, vector<int> &attrs, secrets &sk)
{
    int opt_w = opt_h.size();
    cout << "\nStarting Decryption (Parallel)" << endl;
    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();

    mpz_t message_mpz;
    mpz_init(message_mpz);
    cout << " Auth = ";
    vector<int> opt_attrs;
    int count = 0;
    for (int &attr_c : attrs)
    {
        if (1 <= attr_dict[attr_c])
            opt_attrs.push_back(attr_c);
        if (count == 0)
            cout << attr_c;
        else
            cout << " and " << attr_c;
        count++;
    }
    cout << endl;

    int len = ciphertext.c4.size();
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
            if (++ps[i] < ciphertext.c4[i].size())
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
        threads.push_back(thread(mul_pair_dec, ref(ciphertext), ref(sk), ref(m_tmp[i]), ref(c_indexs[i]), ref(opt_attrs)));
    }
    cout << " " << threads.size() << " threads started." << endl;
    
    element_t pair1, pair2, pair3, m;
    element_init_GT(pair1, pairing);
    element_init_GT(pair2, pairing);
    element_init_GT(m, pairing);

    element_pairing(pair1, ciphertext.c3, sk.sk2);
    element_pairing(pair2, ciphertext.c2, sk.sk1);
    element_mul(m, ciphertext.c1, pair1);
    element_div(m, m, pair2);

    for (int i = 0; i < thread_n; i++)
    {
        threads[i].join();
        element_mul(m_tmp[i], m_tmp[i], m);
        element_to_mpz(message_mpz, m_tmp[i]);
        valueToMessage(dec_m, message_mpz);
        element_clear(m_tmp[i]);
        dec_m_tmp[i] = string(dec_m);
    }

    element_clear(pair1);
    element_clear(pair2);
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

int main(int argc, char **argv)
{
    // Setup
    Setup(argc, argv);

    // Encrypt
    cipher ct;                            // ciphertext
    vector<int> w{0, 2, 3, 6, 7, 10, 11}; // condition
    string message("hello world!!");      // message
    Encrypt(ct, message, w);

    // Decrypt user1
    vector<int> user1{0, 2, 6, 10};
    secrets user1_sk;
    KeyGen(user1, user1_sk); // Key Genaration
    char dec_message1[2048];
    P_Decrypt(ct, dec_message1, user1, user1_sk);
    vector<int> cond_c1{0, 0, 0};
    Decrypt(ct, cond_c1, dec_message1, user1, user1_sk);

    // Decrypt user2
    vector<int> user2{1, 2, 6, 10};
    secrets user2_sk;
    KeyGen(user2, user2_sk); // Key Genaration
    char dec_message2[2048];
    P_Decrypt(ct, dec_message2, user2, user2_sk);
    
    // Decrypt user3
    vector<int> user3{0, 3, 7, 11};
    secrets user3_sk;
    KeyGen(user3, user3_sk); // Key Genaration
    char dec_message3[2048];
    P_Decrypt(ct, dec_message3, user3, user3_sk);

    // Decrypt user4
    vector<int> user4{0, 3, 7, 12};
    secrets user4_sk;
    KeyGen(user4, user4_sk); // Key Genaration
    char dec_message4[2048];
    P_Decrypt(ct, dec_message4, user4, user4_sk);

    cout << "\nCipher text size = " << Get_cipher_size(ct) << " bytes" << endl;
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
