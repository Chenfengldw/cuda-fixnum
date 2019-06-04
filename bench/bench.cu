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
    __device__ void operator()(fixnum &ct, fixnum n, fixnum m, fixnum r) {
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
        assert(result == 0);
    }
};


// fn_bytes is the size of fixnum, which should be able to contain n^2 which is 4 times the size of p and q.
// data_bytes is used for data size, content to be encrytped.
// nelts is number of elements.
template<int fn_bytes, int data_bytes, typename word_fixnum>
void my_bench_paillier(int nelts, const char* filename)  {

    typedef warp_fixnum<fn_bytes, word_fixnum> fixnum;
    typedef fixnum_array<fixnum> fixnum_array;

    uint8_t *input1, *input2, *input3, *input4;
    fixnum_array *cipher, *in1, *in2, *in3, *in4, *n, *res;
    
    // half_fn_bytes is the actual size of p in byte
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
    input1 = new uint8_t[half_fn_bytes * nelts];
    input2 = new uint8_t[half_fn_bytes * nelts];

    for (int i = 0; i < nelts; ++i){
        memcpy(input1+i*half_fn_bytes,p,half_fn_bytes);
        memcpy(input2+i*half_fn_bytes,q,half_fn_bytes);   
    }

    input3 = new uint8_t[data_bytes * nelts];
    for (int i = 0; i < data_bytes * nelts; ++i)
        input3[i] = (i * 23 + 11) % 256;
    
    input4 = new uint8_t[half_fn_bytes * nelts];
    for (int i = 0; i < half_fn_bytes * nelts; ++i)
        input4[i] = (i * 19 + 11) % 256;

    in1 = fixnum_array::create(input1, half_fn_bytes * nelts, half_fn_bytes);
    in2 = fixnum_array::create(input2, half_fn_bytes * nelts, half_fn_bytes);
    in3 = fixnum_array::create(input3, data_bytes * nelts, data_bytes);
    in4 = fixnum_array::create(input4, half_fn_bytes * nelts, half_fn_bytes);
    cipher = fixnum_array::create(nelts);        
    //output = fixnum_array::create(nelts);
    n = fixnum_array::create(nelts);
    res = fixnum_array::create(nelts);


    

    if (nelts == 0) {
        puts(" -*-  nelts == 0; skipping...  -*-");
        return;
    }


    fixnum_array::template map<mul_lo>(n, in1, in2);    

    //cout << "txt before encrypt" << in3 << endl;
    clock_t c = clock();
    fixnum_array::template map<pencrypt>(cipher, n, in3, in4);    
    c = clock() - c;
    
    double secinv = (double)CLOCKS_PER_SEC / c;
    double total_MiB = fixnum::BYTES * (double)nelts / (1 << 20);
    
    printf("enc    %4d   %3d    %6.1f   %7.3f  %12.3f\n",
           half_fn_bytes*16, data_bytes*8, total_MiB,
           1/secinv, nelts * 1e-3 * secinv);

  
    c = clock();
    fixnum_array::template map<pdecrypt_no_check>(res, cipher, in1, in2);
    c = clock() - c;

    
    secinv = (double)CLOCKS_PER_SEC / c;
    total_MiB = fixnum::BYTES * (double)nelts / (1 << 20);

    printf("dec    %4d   %3d    %6.1f   %7.3f  %12.3f\n",
           half_fn_bytes*16, data_bytes*8, total_MiB,
           1/secinv, nelts * 1e-3 * secinv);
    
    //cout << "txt before encrypt" << in3 << endl;
    //cout << "txt after encrypt" << cipher << endl;
    //cout << "txt after decrypt" << res << endl;

    /*cout << "p: " << in1 << endl;
    cout << "q: " << in2 << endl;
    cout << "n: " << n << endl;
    cout << "input: " << in3 <<endl;
    cout << "rand: " << in4 << endl;   
    cout << "encrypted: " << cipher <<endl;  
    cout << "decrypted: " << res <<endl;
    cout << "check decrypted==input: " << output <<endl;*/

    //bool *in_out_same = new bool[nelts];    
    //fixnum_array::template map<cmp>(in3, res);

    
    cout<<endl;


    delete in1; delete in2; delete in3; delete in4;
    delete input1; delete input2; delete input3; delete input4;
    delete cipher;delete n;delete res;

}

void my_bench_func(const char *fn_name, int nelts) {
    printf("Function: %s, #elts: %d\n", fn_name, nelts);
    printf("func  key_size  data_size  total data   time       Kops/s\n");
    printf("       bits       bits      (MiB)     (seconds)\n"); 

    cout<<endl;

    //here key_128 means the file size is 128, which contains p and q, so p q size is half of 128

    my_bench_paillier<8, 4, u64_fixnum>(nelts,"./generated_keys/key_32");
    my_bench_paillier<16, 4, u64_fixnum>(nelts,"./generated_keys/key_64");
    my_bench_paillier<32, 4, u64_fixnum>(nelts,"./generated_keys/key_128");
    my_bench_paillier<64, 4, u64_fixnum>(nelts,"./generated_keys/key_256");
    my_bench_paillier<128, 4, u64_fixnum>(nelts,"./generated_keys/key_512");
    my_bench_paillier<256,4, u64_fixnum>(nelts,"./generated_keys/key_1024");
    //my_bench_paillier<512,4, u32_fixnum>(nelts,"./generated_keys/key_2048"); 
    //my_bench_paillier<1024,4, u32_fixnum>(nelts,"./generated_keys/key_4096");
    my_bench_paillier<8, 8, u64_fixnum>(nelts,"./generated_keys/key_32");
    my_bench_paillier<16, 8, u64_fixnum>(nelts,"./generated_keys/key_64");
    my_bench_paillier<32, 8, u64_fixnum>(nelts,"./generated_keys/key_128");
    my_bench_paillier<64, 8, u64_fixnum>(nelts,"./generated_keys/key_256");
    my_bench_paillier<128, 8, u64_fixnum>(nelts,"./generated_keys/key_512");
    my_bench_paillier<256, 8, u64_fixnum>(nelts,"./generated_keys/key_1024");
    //my_bench_paillier<512, 8, u64_fixnum>(nelts,"./generated_keys/key_2048");
    //my_bench_paillier<1024, 8, u64_fixnum>(nelts,"./generated_keys/key_4096");
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
