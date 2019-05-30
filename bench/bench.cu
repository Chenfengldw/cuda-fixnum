#include <cstdio>
#include <cstring>
#include <cassert>
#include <assert.h>

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



template< typename fixnum >
struct mul_lo {
    __device__ void operator()(fixnum &r, fixnum a) {
        fixnum s;
        fixnum::mul_lo(s, a, a);
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
struct sqr_wide {
    __device__ void operator()(fixnum &r, fixnum a) {
        fixnum rr, ss;
        fixnum::sqr_wide(ss, rr, a);
        r = ss;
    }
};

template< typename fixnum >
struct encrypt {
	__device__ void operator()(fixnum &ctxt, fixnum p, fixnum q, fixnum m, fixnum r) {

        fixnum n;
        fixnum::mul_lo(n, p, q);
		paillier_encrypt<fixnum> enc(n);
		enc(ctxt, m, r);
	}
};

template< typename fixnum >
struct decrypt {
	__device__ void operator()(fixnum &r, fixnum p, fixnum q, fixnum hi, fixnum lo) {
		paillier_decrypt<fixnum> dec(p, q);
		fixnum rr;
		dec(rr, fixnum::zero(), r);
		r = rr;
	}
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


template< typename modnum >
struct my_modexp {
    typedef typename modnum::fixnum fixnum;

    __device__ void operator()(fixnum &z, fixnum x) {
        modexp<modnum> me(x, x);
        fixnum zz;
        me(zz, x);
        z = zz;
    };
};

template< typename modnum >
struct my_multi_modexp {
    typedef typename modnum::fixnum fixnum;

    __device__ void operator()(fixnum &z, fixnum x) {
        multi_modexp<modnum> mme(x);
        fixnum zz;
        mme(zz, x, x);
        z = zz;
    };
};

template< int fn_bytes, typename word_fixnum, template <typename> class Func >
void bench(int nelts) {
    typedef warp_fixnum<fn_bytes, word_fixnum> fixnum;
    typedef fixnum_array<fixnum> fixnum_array;

    if (nelts == 0) {
        puts(" -*-  nelts == 0; skipping...  -*-");
        return;
    }

    uint8_t *input = new uint8_t[fn_bytes * nelts];
    for (int i = 0; i < fn_bytes * nelts; ++i)
        input[i] = (i * 17 + 11) % 256;

    fixnum_array *res, *in;
    in = fixnum_array::create(input, fn_bytes * nelts, fn_bytes);
    res = fixnum_array::create(nelts);


    clock_t c = clock();
    fixnum_array::template map<Func>(res, in);
    c = clock() - c;

    double secinv = (double)CLOCKS_PER_SEC / c;
    double total_MiB = fixnum::BYTES * (double)nelts / (1 << 20);
    printf(" %4d   %3d    %6.1f   %7.3f  %12.1f\n",
           fixnum::BITS, fixnum::digit::BITS, total_MiB,
           1/secinv, nelts * 1e-3 * secinv);

    delete in;
    delete res;
    delete[] input;
}

template< template <typename> class Func >
void bench_func(const char *fn_name, int nelts) {
    printf("Function: %s, #elts: %de3\n", fn_name, (int)(nelts * 1e-3));
    printf("fixnum digit  total data   time       Kops/s\n");
    printf(" bits  bits     (MiB)    (seconds)\n");
    bench<4, u32_fixnum, Func>(nelts);
    bench<8, u32_fixnum, Func>(nelts);
    bench<16, u32_fixnum, Func>(nelts);
    bench<32, u32_fixnum, Func>(nelts);
    bench<64, u32_fixnum, Func>(nelts);
    bench<128, u32_fixnum, Func>(nelts);
    puts("");

    bench<8, u64_fixnum, Func>(nelts);
    bench<16, u64_fixnum, Func>(nelts);
    bench<32, u64_fixnum, Func>(nelts);
    bench<64, u64_fixnum, Func>(nelts);
    bench<128, u64_fixnum, Func>(nelts);
    bench<256, u64_fixnum, Func>(nelts);
    puts("");
}

template< int fn_bytes, typename word_fixnum, template <typename> class Func >
void my_bench(int nelts)  {

    typedef warp_fixnum<fn_bytes, word_fixnum> fixnum;
    typedef fixnum_array<fixnum> fixnum_array;
    
    uint8_t *input1, *input2, *input3, *input4;
    fixnum_array *res, *output, *in1, *in2, *in3, *in4;


    input1 = new uint8_t[fn_bytes * nelts];
    for (int i = 0; i < fn_bytes * nelts; ++i)
        input1[i] = (i * 17 + 11) % 256;
    input2 = new uint8_t[fn_bytes * nelts];
    for (int i = 0; i < fn_bytes * nelts; ++i)
        input2[i] = (i * 13 + 11) % 256;
    input3 = new uint8_t[fn_bytes * nelts];
    for (int i = 0; i < fn_bytes * nelts; ++i)
        input3[i] = (i * 23 + 11) % 256;
    input4 = new uint8_t[fn_bytes * nelts];
    for (int i = 0; i < fn_bytes * nelts; ++i)
        input4[i] = (i * 19 + 11) % 256;

    in1 = fixnum_array::create(input1, fn_bytes * nelts, fn_bytes);
    in2 = fixnum_array::create(input2, fn_bytes * nelts, fn_bytes);
    in3 = fixnum_array::create(input3, fn_bytes * nelts, fn_bytes);
    in4 = fixnum_array::create(input4, fn_bytes * nelts, fn_bytes);
    res = fixnum_array::create(nelts);
    output = fixnum_array::create(nelts);

   

    if (nelts == 0) {
        puts(" -*-  nelts == 0; skipping...  -*-");
        return;
    }

    clock_t c = clock();
    //fixnum_array::template map<Func>(res, in1, in2, in3, in4);
    fixnum_array::template map<encrypt>(res, in1, in2, in3, in4);
    c = clock() - c;

    double secinv = (double)CLOCKS_PER_SEC / c;
    double total_MiB = fixnum::BYTES * (double)nelts / (1 << 20);
    printf("enc   %4d   %3d    %6.1f   %7.3f  %12.1f\n",
           fixnum::BITS, fixnum::digit::BITS, total_MiB,
           1/secinv, nelts * 1e-3 * secinv);

    
  
    c = clock();
    //fixnum_array::template map<Func>(res, in1, in2, in3, in4);
    fixnum_array::template map<decrypt>(output, in1, in2, in3, in4);
    c = clock() - c;

    secinv = (double)CLOCKS_PER_SEC / c;
    total_MiB = fixnum::BYTES * (double)nelts / (1 << 20);
    printf("dec   %4d   %3d    %6.1f   %7.3f  %12.1f\n",
           fixnum::BITS, fixnum::digit::BITS, total_MiB,
           1/secinv, nelts * 1e-3 * secinv);


    //bool *in_out_same = new bool[nelts];    
    fixnum_array::template map<cmp>(res, output);
    cerr << cudaGetErrorString(cudaPeekAtLastError()) << endl;


    
}

template< template <typename> class Func >
void my_bench_func(const char *fn_name, int nelts) {

    printf("Function: %s, #elts: %de3\n", fn_name, (int)(nelts * 1e-3));
    printf("func  fixnum digit  total data   time       Kops/s\n");
    printf("       bits  bits     (MiB)    (seconds)\n");
    my_bench<4, u32_fixnum, Func>(nelts);
    my_bench<8, u32_fixnum, Func>(nelts);
    my_bench<16, u32_fixnum, Func>(nelts);
    my_bench<32, u32_fixnum, Func>(nelts);
    my_bench<64, u32_fixnum, Func>(nelts);
    my_bench<128, u32_fixnum, Func>(nelts);
    puts("");

    my_bench<8, u64_fixnum, Func>(nelts);
    my_bench<16, u64_fixnum, Func>(nelts);
    my_bench<32, u64_fixnum, Func>(nelts);
    my_bench<64, u64_fixnum, Func>(nelts);
    my_bench<128, u64_fixnum, Func>(nelts);
    my_bench<256, u64_fixnum, Func>(nelts);
    puts("");
}


template< typename fixnum >
using modexp_redc = my_modexp< modnum_monty_redc<fixnum> >;

template< typename fixnum >
using modexp_cios = my_modexp< modnum_monty_cios<fixnum> >;

template< typename fixnum >
using multi_modexp_redc = my_multi_modexp< modnum_monty_redc<fixnum> >;

template< typename fixnum >
using multi_modexp_cios = my_multi_modexp< modnum_monty_cios<fixnum> >;


int main(int argc, char *argv[]) {
    long m = 1;
    if (argc > 1)
        m = atol(argv[1]);
    m = std::max(m, 1000L);
    my_bench_func<encrypt_decrypt>("encrypt_decrypt", m);
    puts("");


    // delete in1, in2, in3, in4;
    // delete res;
    // delete[] input1, input2, input3, input4;

    return 0;
}
