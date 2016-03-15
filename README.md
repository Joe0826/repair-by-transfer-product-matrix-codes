# repair-by-transfer-product-matrix-codes

This repository is associated with the USENIX FAST 2015 paper ["Having Your Cake and Eating It Too: Jointly Optimal Erasure Codes for I/O, Storage, and Network-bandwidth", K. V. Rashmi, Preetum Nakkiran, Jingyan Wang, Nihar B. Shah, and Kannan Ramchandran](https://www.usenix.org/conference/fast15/technical-sessions/presentation/rashmi).

Most of this implementation is done by two UC Berkeley undergraduate students Preetum Nakkiran and Jingyan Wang.

The work in the above mentioned USENIX FAST 2015 paper builds on top of a powerful class of distributed storage codes called the "Product-Matrix" codes  proposed in the IEEE Transactions on Information Theory 2011 paper ["Optimal Exact-Regenerating Codes for Distributed Storage at the MSR and MBR Points via a Product-Matrix Construction", K. V. Rashmi, Nihar B. Shah, P. Vijay Kumar](http://eecs.berkeley.edu/~rashmikv/papers/product_matrix_codes.pdf). 

This repository also contains an implementation of the vanilla Product-Matrix codes. Ankush Gupta and Diivanand Ramalingam, two other UC Berkeley undergraduate students, also contributed to the implementation of the vanilla Product-Matrix codes as a part of their [EE121](https://inst.eecs.berkeley.edu/~ee121/fa13/) course project.


##Instructions for compiling:
1. This project needs gf-complete and jerasure libraries which are available at http://jerasure.org/groups/jerasure. First download and install these two libraries in this "parent" directory following the instructions as specified by gf-complete and jerasure respectively. After this step, you should have two directories, gf-complete-kmg and jerasure-kmg, as siblings to the pmc directory and both installed. 

2. cd pmc and make

3. cd pmc/filecoder and make


## Structure
The core product-matrix routines are in [pmc](./pmc) directory.

File-based encoder/collector/repair examples are in [pmc/filecoder](./pmc/filecoder). Documentation and usage examples for filecoder are in the corresponding [README](./pmc/filecoder/README.md) file.


## Extras
[pmc/filecoder](./pmc/filecoder) also contains a toy distributed storage application: thin python wrappers for filecoder that enable storage, collection, and repair over the network.

Some other simple local examples are in [pmc/tests](./pmc/tests) directory.

## More details: 
Uses Jerasure (http://web.eecs.utk.edu/~plank/plank/papers/CS-08-627.html) library for efficient Galois field and matrix operations. 

Implements the MSR Product-Matrix codes for parameters d > 2k-2 using the unified form proposed in the IEEE Transactions on Information Theory 2015 paper ["An Unified Form of Exact-MSR Codes via
Product-Matrix Framework", Lin, Chung](http://www.citi.sinica.edu.tw/papers/whc/3668-F.pdf).

