# Generalized-Weak-incNGs-Consistency
Generalized weak-incNGs consistency


This is the implementation of the Generalized-Weak-incNGs-Consistency (GWIC) on the incNGs global constraint proposed in the paper "Jimmy H.M. Lee and Zichen Zhu. Filtering Nogoods Lazily in Dynamic Symmetry Breaking During Search, Proceedings of the 24th International Joint Conference on Artificial Intelligence (IJCAI 2015), pages 339-345, Buenos Aires, Argentina, July, 2015"

Please use Gecode Solver 4.2.0 to run these files.

Put GWIC folder into gecode folder and efpa_lresbds_GWIC.cpp file into example folder.

To run the EFPA problem (e.g. (5 3 3 4)) using one of the symmetry breaking methods NONE/DoubleLex/SnakeLex/LDSB/LReSBDS:

./efpa_lresbds_GWIC -search none 5 3 3 4

./efpa_lresbds_GWIC -search double 5 3 3 4

./efpa_lresbds_GWIC -search snake 5 3 3 4

./efpa_lresbds_GWIC -search ldsb 5 3 3 4

./efpa_lresbds_GWIC -search lresbds 5 3 3 4
