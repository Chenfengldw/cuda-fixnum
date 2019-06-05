#include <cstdio>
#include <cstring>
#include <cassert>
#include <assert.h>
#include <iostream>
#include <fstream>

#include "fixnum/warp_fixnum.cu"
#include "array/fixnum_array.h"
#include "functions/modexp.cu"
#include "functions/multi_modexp.cu"
#include "modnum/modnum_monty_redc.cu"
#include "modnum/modnum_monty_cios.cu"
#include "functions/paillier_encrypt.cu"
#include "functions/paillier_decrypt.cu"

using namespace std;
using namespace cuFIXNUM;


#include "functions/quorem_preinv.cu"


template< typename fixnum >
struct add {
    __device__ void operator()(fixnum &r, fixnum a, fixnum b) {
        fixnum s;
        fixnum::add(s, a, b);
        r = s;
    }
};


template< typename fixnum >
struct mod_square {
    __device__ void operator()(fixnum &r, fixnum n, fixnum a, fixnum b) {
        fixnum n2,s;        
        fixnum::sqr_lo(n2,n);
        quorem_preinv<fixnum> mod_n2(n2);        
        fixnum c_hi,c_lo;
        fixnum::mul_wide(c_hi,c_lo, a, b);
        mod_n2(s, c_hi, c_lo);
        r = s;
    }
};


template< typename fixnum >
struct mul_lo {
    __device__ void operator()(fixnum &r, fixnum a, fixnum b) {
        fixnum s;
        fixnum::mul_lo(s, a, b);
        r = s;
    }
};

template< typename fixnum >
struct mul_wide {
    __device__ void operator()(fixnum &r, fixnum a) {
        fixnum rr, ss;
        fixnum::mul_wide(ss, rr, a, a);
        r = ss;
    }
};


template< typename fixnum >
struct pencrypt {
    __device__ void operator()(fixnum &ct,fixnum n, fixnum m, fixnum r) {        
        fixnum kk;         
        paillier_encrypt<fixnum> enc(n);
        enc(kk, m, r);        
        ct = kk;
    };
};



template< typename fixnum >
struct pdecrypt {
    __device__ void operator()(fixnum &z, fixnum & res,fixnum ct, fixnum p, fixnum q, fixnum m, fixnum r) {
        if (fixnum::cmp(p, q) == 0
              || fixnum::cmp(r, p) == 0
              || fixnum::cmp(r, q) == 0) {
            z = fixnum::zero();
            return;
        }
        paillier_decrypt<fixnum> dec(p, q);
        fixnum zz;
        dec(zz, fixnum::zero(), ct);
        res = zz;
        // z = (z != m)
        z = fixnum::digit( !! fixnum::cmp(zz, m));
    };
};


template< typename fixnum >
struct pdecrypt_no_check {
    __device__ void operator()(fixnum & res,fixnum ct, fixnum p, fixnum q) {        
        paillier_decrypt<fixnum> dec(p, q);
        fixnum zz;
        dec(zz, fixnum::zero(), ct);
        res = zz;
    };
};


template< typename fixnum >
struct encrypt_decrypt {
	__device__ void operator()(fixnum &r, fixnum in1, fixnum in2, fixnum in3, fixnum in4) {

        fixnum n;
        fixnum::mul_lo(n, in1, in2);
        paillier_encrypt<fixnum> enc(n);
        paillier_decrypt<fixnum> dec(in1, in2);
        enc(r, in3, in4);
        fixnum rr;
		dec(rr, fixnum::zero(), r);
        r = rr;
	}
};

template< typename fixnum >
struct cmp {
    __device__ void operator()(fixnum x, fixnum y){
        int result = fixnum::cmp(x, y);
        //assert(result == 0);
    }
};


// fn_bytes is how many bytes a single element has, and is used for key size here.
// word_fixnum type size, which is also word_fixnum::BYTES is used for digit size, content to be encrytped.
// nelts is number of elements.
template<int fn_bytes, int data_bytes, typename word_fixnum>
void my_bench_paillier(int nelts, const char* filename)  {

    typedef warp_fixnum<fn_bytes, word_fixnum> fixnum;
    typedef fixnum_array<fixnum> fixnum_array;

    uint8_t *input1, *input2, *input3, *input4;
    fixnum_array *cipher,*cipher2,*cipher3, *output, *in1, *in2, *in3, *in4, *p_arr,*q_arr,*n_arr, *res,*res2;
    
    int half_fn_bytes = fn_bytes / 4;

    ifstream data(filename, ios_base::binary);
    uint8_t* p = new uint8_t[half_fn_bytes];
    uint8_t* q = new uint8_t[half_fn_bytes];

    data.seekg(0U, ios_base::beg); // Move the input position indicator to the beginning of the file for reading
    data.read(reinterpret_cast<char*>(p), half_fn_bytes); 
    data.read(reinterpret_cast<char*>(q), half_fn_bytes);
    data.close();



    /**fixnum_array* pkey;
    pkey = fixnum_array::create(p, half_fn_bytes * 1, half_fn_bytes);
    cout << pkey << endl;

    fixnum_array* qkey;
    qkey = fixnum_array::create(q, half_fn_bytes * 1, half_fn_bytes);
    cout << qkey << endl;**/


    // % 256 because one byte at most represent 256 numbers.

    input1= new uint8_t[data_bytes * nelts];
    for (int i = 0; i < data_bytes * nelts; ++i)
        input1[i] = (i * 17 + 11) % 256;
    
    input2 = new uint8_t[half_fn_bytes * nelts];
    for (int i = 0; i < half_fn_bytes * nelts; ++i)
        input2[i] = (i * 13 + 11) % 256;

    input3 = new uint8_t[data_bytes * nelts];
    for (int i = 0; i < data_bytes * nelts; ++i)
        input3[i] = (i * 23 + 11) % 256;
    
    input4 = new uint8_t[half_fn_bytes * nelts];
    for (int i = 0; i < half_fn_bytes * nelts; ++i)
        input4[i] = (i * 19 + 11) % 256;

    p_arr = fixnum_array::create(p, half_fn_bytes, half_fn_bytes);
    q_arr = fixnum_array::create(q, half_fn_bytes, half_fn_bytes);   
    
    my_paillier_decrypt<fixnum> my_decrypt(p_arr,q_arr);

    p_arr = p_arr->repeat(nelts);
    q_arr = q_arr->repeat(nelts);
    n_arr = fixnum_array::create(nelts);

    in1 = fixnum_array::create(input1, data_bytes * nelts, data_bytes);
    in2 = fixnum_array::create(input2, half_fn_bytes * nelts, half_fn_bytes);
    in3 = fixnum_array::create(input3, data_bytes * nelts, data_bytes);
    in4 = fixnum_array::create(input4, half_fn_bytes * nelts, half_fn_bytes);

    output = fixnum_array::create(nelts);    

    cipher = fixnum_array::create(nelts);            
    res = fixnum_array::create(nelts);
    
    cipher2 = fixnum_array::create(nelts);       
    res2 = fixnum_array::create(nelts);

    cipher3 = fixnum_array::create(nelts);


    if (nelts == 0) {
        puts(" -*-  nelts == 0; skipping...  -*-");
        return;
    }

    fixnum_array::template map<mul_lo>(n_arr,p_arr,q_arr);    

    clock_t c = clock();
    fixnum_array::template map<pencrypt>(cipher,n_arr,in1,in2);    
    c = clock() - c;
    
    double secinv = (double)CLOCKS_PER_SEC / c;
    double total_MiB = fixnum::BYTES * (double)nelts / (1 << 20);
    
    printf("enc    %4d   %3d    %6.1f   %7.3f  %12.1f\n",
           fixnum::BITS, fixnum::digit::BITS, total_MiB,
           1/secinv, nelts * 1e-3 * secinv);


    c = clock();
    fixnum_array::template map<pencrypt>(cipher2,n_arr,in3,in4);    
    c = clock() - c;
    
    secinv = (double)CLOCKS_PER_SEC / c;
    total_MiB = fixnum::BYTES * (double)nelts / (1 << 20);

    printf("enc2   %4d   %3d    %6.1f   %7.3f  %12.1f\n",
           fixnum::BITS, fixnum::digit::BITS, total_MiB,
           1/secinv, nelts * 1e-3 * secinv);


    c = clock();
    fixnum_array::template map<mod_square>(output,n_arr,cipher,cipher2);    
    c = clock() - c;
    
    secinv = (double)CLOCKS_PER_SEC / c;
    total_MiB = fixnum::BYTES * (double)nelts / (1 << 20);

    printf("add    %4d   %3d    %6.1f   %7.3f  %12.1f\n",
           fixnum::BITS, fixnum::digit::BITS, total_MiB,
           1/secinv, nelts * 1e-3 * secinv);
    
    

    c = clock();    
    fixnum_array::template map<my_paillier_decrypt>(my_decrypt,res,output);
    c = clock() - c;

    secinv = (double)CLOCKS_PER_SEC / c;
    total_MiB = fixnum::BYTES * (double)nelts / (1 << 20);

    printf("dec    %4d   %3d    %6.1f   %7.3f  %12.1f\n",
           fixnum::BITS, fixnum::digit::BITS, total_MiB,
           1/secinv, nelts * 1e-3 * secinv);


    fixnum_array::template map<add>(res2,in1,in3);
    fixnum_array::template map<pencrypt>(cipher3,n_arr,res2,in2);    

    /*cout << "p: " << p_arr << endl;
    cout << "q: " << q_arr << endl;    
    cout << "input1: " << in1 <<endl;
    cout << "input2: " << in3 << endl;   
    cout << "encrypted 1: " << cipher <<endl;  
    cout << "encrypted 2: " << cipher2 <<endl;  
    cout << "add result: " << output <<endl;
    cout << "add result2: " << cipher3 <<endl;
    cout << "decrypted: " << res <<endl;
    cout << "for check: " << res2 <<endl;*/
    

    //bool *in_out_same = new bool[nelts];    
    //fixnum_array::template map<cmp>(in3, res);

    
    cout<<endl;

    delete p_arr;delete q_arr;
    delete in1; delete in2; delete in3; delete in4;
    delete res;delete res2;
    delete cipher;delete cipher2;
    delete input1; delete input2;delete input3; delete input4;
    delete output;

}

void my_bench_func(const char *fn_name, int nelts) {
    printf("Function: %s, #elts: %d\n", fn_name, nelts);
    printf("func  fixnum digit  total data   time       Kops/s\n");
    printf("       bits  bits     (MiB)    (seconds)\n"); 

    cout<<endl;
    /**my_bench_paillier<16, 4, u32_fixnum>(nelts,"./generated_keys/key_128");
    my_bench_paillier<32, 4, u32_fixnum>(nelts,"./generated_keys/key_256");
    my_bench_paillier<64, 4, u32_fixnum>(nelts,"./generated_keys/key_512");
    my_bench_paillier<128,4, u32_fixnum>(nelts,"./generated_keys/key_1024");**/
    my_bench_paillier<32, 4, u64_fixnum>(nelts,"./generated_keys/key_128");
    my_bench_paillier<64, 4, u64_fixnum>(nelts,"./generated_keys/key_256");
    my_bench_paillier<128, 4, u64_fixnum>(nelts,"./generated_keys/key_512");
    my_bench_paillier<256, 4, u64_fixnum>(nelts,"./generated_keys/key_1024");
    puts("");
}


int main(int argc, char *argv[]) {
    long m = 1;
    if (argc > 1)
        m = atol(argv[1]);
    m = std::max(m, 1L);    
    
    //my_bench_mul_lo<32, u32_fixnum, mul_lo>(1);
    my_bench_func("encrypt_decrypt", m);
    puts("");

    return 0;
}
