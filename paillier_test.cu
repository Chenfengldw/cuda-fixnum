#include <cstdlib>
#include <cassert>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <time.h>
#include <stdlib.h>

#include "fixnum/warp_fixnum.cu"
#include "functions/paillier_encrypt.cu"
#include "functions/paillier_decrypt.cu"
#include "fixnum/warp_fixnum.cu"
#include "array/fixnum_array.h"

using namespace std;
using namespace cuFIXNUM;

#define KEY_BITS 1024
#define KEY_BYTES ((KEY_BITS)/8)
#define KEY_DIGITS ((KEY_BYTES)/4)

template< typename fixnum >
fixnum* createFixnum() {
	uint8_t* bytes = (uint8_t*)malloc(fixnum::BYTES);
	for (int i = 0; i < fixnum::BYTES; i++) {
		bytes[i] = rand() % 256;
	}
	fixnum* res;
	cudaMallocManaged((void**)&res, fixnum::BYTES);
	fixnum::from_bytes(reinterpret_cast<uint8_t*>(res), bytes, fixnum::BYTES);
	return res;
}

template< typename fixnum >
string fixnumToString(fixnum* fn) {
	ostringstream ss;
	uint8_t* bytes = reinterpret_cast<uint8_t*>(fn);

	for (int i = fixnum::BYTES - 1; i >= 0; i--) {
		ss << hex << setw(2) << setfill('0') << (int)bytes[i];
	}
	return ss.str();
}

template< typename fixnum >
struct Add {
	__device__ void operator()(fixnum& res, fixnum a, fixnum b) {
		fixnum s;
		fixnum::add(s, a, b);
		res = s;
	}
};

template< typename fixnum >
struct Decrypt {
	//paillier_decrypt<fixnum> _dec;
	//__device__ Decrypt(fixnum p, fixnum q) : _dec(p, q) {}

	__device__ void operator()(fixnum &ptxt, fixnum ctxt_hi, fixnum ctxt_lo, fixnum p, fixnum q) {
		paillier_decrypt<fixnum> _dec(p, q);
		fixnum tmp;
		_dec(tmp, ctxt_hi, ctxt_lo);
		ptxt = tmp;
	}
};

// Encrypt added at 5.28 dw

template< typename fixnum >
struct Encrypt {

	__device__ void operator()(fixnum &ctxt, fixnum m, fixnum r, fixnum n) {
		paillier_encrypt<fixnum> _enc(n);
		_enc(ctxt, m, r);
	}
};



// template< typename fixnum >
// struct DecryptFunc {
// 	paillier_decrypt<fixnum> _dec;
// 	__device__ Decrypt(fixnum p, fixnum q) : _dec(p, q) {}

// 	__device__ void operator()(fixnum &ptxt, fixnum ctxt_hi, fixnum ctxt_lo) {
// 		fixnum tmp;
// 		_dec(tmp, ctxt_hi, ctxt_lo);
// 		ptxt = tmp;
// 	}
// };

template< typename fixnum >
__global__ void dispatch(fixnum* res, fixnum* a, fixnum* b) {
	int bid = blockIdx.x * blockDim.x;
	int tid = threadIdx.x;
	int idx = bid + tid;

	Add<fixnum> add_fn;

	if (idx < KEY_DIGITS) {
		add_fn(res[idx], a[idx], b[idx]);
	} else {
		return;
	}
}

template< typename fixnum >
__global__ void dispatch_func(fixnum* res, fixnum* hi, fixnum* lo, fixnum* p, fixnum* q) {
	int bid = blockIdx.x * blockDim.x;
	int tid = threadIdx.x;
	int idx = bid + tid;

	//Func<fixnum> fn;
	Decrypt<fixnum> dec;

	if (idx < KEY_DIGITS) {
		dec(res[idx], hi[idx], lo[idx], p[idx], q[idx]);
	} else {
		return;
	}
}


// add by dw 5.29
template< typename fixnum >
__global__ void dispatch_func(fixnum* sum, fixnum* p, fixnum* q, int testTime) {
	int bid = blockIdx.x * blockDim.x;
	int tid = threadIdx.x;
	int idx = bid + tid;

	Encrypt<fixnum> enc;
	Decrypt<fixnum> dec;
	fixnum n;

	if (idx < KEY_DIGITS) {
		for(int i=0; i<testTime; ++i){
			fixnum::mul_digit(n, p[idx], q[idx]);
			enc(sum[idx], p[idx], q[idx], n);
			dec(sum[idx], p[idx], q[idx], p[idx], q[idx]);
		}


	} else {
		return;
	}
}

// template< typename fixnum >
// __global__ void dispatch_func(fixnum* res, fixnum* hi, fixnum* lo, DecryptFunc& dec) {
// 	int bid = blockIdx.x * blockDim.x;
// 	int tid = threadIdx.x;
// 	int idx = bid + tid;

// 	if (idx < KEY_DIGITS) {
// 		dec(res[idx], hi[idx], lo[idx]);
// 	} else {
// 		return;
// 	}
// }


//dw to test encrypt and decrypt


template< typename fixnum >
__global__ void parallel_tests_func(int testTime, double& timeCost, fixnum* p, fixnum* q, fixnum* sum){


	int bid = blockIdx.x * blockDim.x;
 	int tid = threadIdx.x;
	int idx = bid + tid;

	//double* timeCost = new double[testTime];
	
	fixnum n;
	fixnum::mul_digit(n, p[idx], q[idx]);
	paillier_encrypt<fixnum> enc(n);
	paillier_decrypt<fixnum> dec(p[idx], q[idx]);


	
	// cudaStream_t stream;
	// cudaStreamCreate(&stream);
	// cudaStreamSynchronize(stream);

	//the following enc and dec run on __device__

	clock_t encStart, decEnd;

	encStart = clock();
	for(int i = 0; i < testTime; i++){
		enc(sum[idx], p[idx], q[idx]);
		dec(sum[idx], p[idx], q[idx]);

	}
	
	// cudaStreamSynchronize(stream);
	// cudaStreamDestroy(stream);

	decEnd = clock();

	timeCost = (double)(decEnd-encStart)/CLOCKS_PER_SEC;
	return;

}


template< typename fixnum >
__global__ void single_tests_func(double& timeCost, fixnum* p, fixnum* q, fixnum* sum){

	int bid = blockIdx.x * blockDim.x;
 	int tid = threadIdx.x;
	int idx = bid + tid;

	//timeCost = new double[testTime];
	
	fixnum n;
	fixnum::mul_digit(n, p[idx], q[idx]);
	paillier_encrypt<fixnum> enc(n);
	paillier_decrypt<fixnum> dec(p[idx], q[idx]);

	//the following enc and dec run on __device__

	clock_t encStart, decEnd;

	encStart = clock();
	enc(sum[idx], p[idx], q[idx]);
	dec(sum[idx], p[idx], q[idx]);
	decEnd = clock();

	timeCost = (double)(decEnd-encStart)/CLOCKS_PER_SEC;

	return ;

}



// int main(){

// 	typedef warp_fixnum<KEY_BYTES, u32_fixnum> fixnum;  
// 	typedef fixnum_array<fixnum> farray;

// 	int testTime = 100;
// 	double runTogetherTime;

// 	fixnum* p = createFixnum<fixnum>();
// 	fixnum* q = createFixnum<fixnum>();
// 	fixnum* sum = createFixnum<fixnum>();

//  	cudaStream_t stream;
// 	cudaStreamCreate(&stream);
// 	//cudaStreamSynchronize(stream);
// 	//parallel_tests_func<fixnum><<<1, 32, 0, stream>>>(testTime, runTogetherTime, p, q, sum);
// 	//cudaStreamSynchronize(stream);
// 	//cout << "run 100 together cost " << runTogetherTime << "ms;" << endl;

// 	double runSingleTime;

// 	cudaStreamSynchronize(stream);
// 	single_tests_func<fixnum><<<1, 32, 0, stream>>>(runSingleTime, p, q, sum);
// 	cudaStreamSynchronize(stream);

// 	//cout << "*********run 100 seperately ***********" << endl;
// 	cout << runSingleTime << "ms;" << endl;
	

// 	return 0;
// }




int main(int arg, char* args[]) {
	typedef warp_fixnum<KEY_BYTES, u32_fixnum> fixnum;    // 64 bits
	typedef fixnum_array<fixnum> farray;

	fixnum* p = createFixnum<fixnum>();
	fixnum* q = createFixnum<fixnum>();
	fixnum* sum = createFixnum<fixnum>();

	//cout << fixnumToString(p) << endl;
	//cout << fixnumToString(q) << endl;

	int testTime;
	if(arg < 2){
		testTime = 1;
	}else {
		testTime = strtol(args[1], NULL, 10);
	}
	

	cudaStream_t stream;
	cudaStreamCreate(&stream);
	cudaStreamSynchronize(stream);

	dispatch<fixnum><<<1, 32, 0, stream>>>(sum, p, q);

	cudaStreamSynchronize(stream);
	cudaStreamDestroy(stream);

	//cout << fixnumToString(sum) << endl;

	//Decrypt<fixnum> dec(*p, *q);
	cudaStreamCreate(&stream);

	dispatch_func<fixnum><<<1, 32, 0, stream>>>(sum, p, q, p, q);

	//add by dw 5.28
	//cout << "Before enc: " << fixnumToString(sum) << endl;

	clock_t startEncTime, endEncTime, endDecTime;

	cudaStreamSynchronize(stream);
	startEncTime = clock();
	dispatch_func<fixnum><<<1, 32, 0, stream>>>(sum, p, q, testTime);
	cudaStreamSynchronize(stream);
	endEncTime = clock();
	//cout << "After dec: " << fixnumToString(sum) << endl;
	cout << "total time use: " << (double)1000*(endEncTime - startEncTime)/CLOCKS_PER_SEC << "ms" << endl;

	cudaStreamDestroy(stream);

	cerr << cudaGetErrorString(cudaPeekAtLastError()) << endl;

	return 0;
}
