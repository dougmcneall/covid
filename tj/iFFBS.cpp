#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::imat iFFBS(
    arma::ivec DI_prime,
    arma::ivec DH_prime,
    arma::ivec t,
    int Npop,
    int Nposs,
    int Niter,
    double pH,
    double pHD,
    double pI1,
    double pI1D,
    double pI1H,
    double pI2,
    double pP,
    double pE,
    double pEP,
    double pA,
    double beta,
    double betaA,
    double pini,
    int print = 0
) {
    // initialise counters and temporary vars
    int i, j, k, l, iter;
    double norm, u;
    
    // states: S, E, A, RA, P, I1, I2, RI, DI, H, RH, DH
    //         0, 1, 2, 3,  4, 5,  6,  7,  8,  9, 10, 11

    // initialise states so they are consistent
    // with the observed data
    arma::imat Xind;
    Xind.zeros(Nposs, t.n_elem);
    k = 0;
    for(j = (t.n_elem - 1); j >= 0; j--) {
        for(i = 0; i < DI_prime(j); i++) {
            Xind(k, j) = 8;
            Xind(k, j - 1) = 6;
            Xind(k, j - 2) = 5;
            Xind(k, j - 3) = 4;
            Xind(k, j - 4) = 1;
            k++;
        }
    }
    for(j = (t.n_elem - 1); j >= 0; j--) {
        for(i = 0; i < DH_prime(j); i++) {
            Xind(k, j) = 11;
            Xind(k, j - 1) = 9;
            Xind(k, j - 2) = 5;
            Xind(k, j - 3) = 4;
            Xind(k, j - 4) = 1;
            k++;
        }
    }
    for(i = 0; i < Nposs; i++) {
        for(j = 1; j < t.n_elem; j++) {
            if(Xind(i, j) != Xind(i, j - 1) && Xind(i, j) == 0) {
                Xind(i, j) = Xind(i, j - 1);
            }
        }
    }
    
    // set up observed counts
    arma::ivec DI(DI_prime.n_elem);
    DI(0) = DI_prime(0);
    for(i = 1; i < DI_prime.n_elem; i++) DI(i) = DI(i - 1) + DI_prime(i);
    arma::ivec DH(DH_prime.n_elem);
    DH(0) = DH_prime(0);
    for(i = 1; i < DH_prime.n_elem; i++) DH(i) = DH(i - 1) + DH_prime(i);
    
    // print outputs if required
    if(print == 1) {
        Rprintf("DI:\n");
        for(int c = 0; c < t.n_elem; c++) {
            Rprintf("%d ", DI(c));
        }
        Rprintf("\nDH:\n");
        for(int c = 0; c < t.n_elem; c++) {
            Rprintf("%d ", DH(c));
        }
        Rprintf("\n");
        
        Rprintf("Xind:\n");
        for(int r = 0; r < Nposs; r++) {
            for(int c = 0; c < t.n_elem; c++) {
                Rprintf("%d ", Xind(r, c));
            }
            Rprintf("\n");
        }
    }
    
    // initialise number in states
    // states: S, E, A, RA, P, I1, I2, RI, DI, H, RH, DH
    //         0, 1, 2, 3,  4, 5,  6,  7,  8,  9, 10, 11
    arma::imat nXsum;
    nXsum.zeros(t.n_elem, 12);
    arma::imat nXcum;
    nXcum.zeros(t.n_elem, 12);
    for(i = 0; i < Nposs; i++) {
        for(j = 0; j < t.n_elem; j++) {
            nXsum(j, Xind(i, j))++;
        }
    }
    
    // print outputs if required
    if(print == 1) {
        Rprintf("nXsum:\n");
        for(int r = 0; r < t.n_elem; r++) {
            for(int c = 0; c < 12; c++) {
                Rprintf("%d ", nXsum(r, c));
            }
            Rprintf("\n");
        }
    }
    
    for(i = 0; i < Nposs; i++) {
        if(Xind(i, 0) > 0) {
            nXcum(0, Xind(i, 0))++;
        }
    }
    for(j = 1; j < t.n_elem; j++) {
        for(i = 0; i < 12; i++) {
            nXcum(j, i) = nXcum(j - 1, i);
        }
        for(i = 0; i < Nposs; i++) {
            if(Xind(i, j - 1) != Xind(i, j)) {
                nXcum(j, Xind(i, j))++;
            }
        }
    }
    
    // print outputs if required
    if(print == 1) {
        Rprintf("nXcum:\n");
        for(int r = 0; r < t.n_elem; r++) {
            for(int c = 0; c < 12; c++) {
                Rprintf("%d ", nXcum(r, c));
            }
            Rprintf("\n");
        }
    }
    
    // set up output vector
    arma::imat output;
    output.zeros(Niter, t.n_elem * 12);
    
    // initialise transition probabilities
    arma::mat pXsum;
    pXsum.zeros(12, 12);
    
    // S to S
    pXsum(0, 0) = exp(-beta * (nXsum(0, 4) + nXsum(0, 5) + nXsum(0, 6) + betaA * nXsum(0, 2)) / ((double) Npop));
    // S to E
    pXsum(0, 1) = (1.0 - pXsum(0, 0)) + pini - (1.0 - pXsum(0, 0)) * pini;
    pXsum(0, 0) = 1.0 - pXsum(0, 1);
    
    // E to E
    pXsum(1, 1) = 1.0 - pE;
    // E to A
    pXsum(1, 2) = (1.0 - pEP) * pE;
    // E to P
    pXsum(1, 4) = pEP * pE;
    
    // A to A
    pXsum(2, 2) = 1.0 - pA;
    // A to RA
    pXsum(2, 3) = pA;
    
    // RA to RA
    pXsum(3, 3) = 1.0;
    
    // P to P
    pXsum(4, 4) = 1.0 - pP;
    // P to I1
    pXsum(4, 5) = pP;
    
    // I1 to I1
    pXsum(5, 5) = 1.0 - pI1;
    // I1 to I2
    pXsum(5, 6) = (1.0 - pI1H - pI1D) * pI1;
    // I1 to DI
    pXsum(5, 8) = pI1D * pI1;
    // I1 to H
    pXsum(5, 9) = pI1H * pI1;
    
    // I2 to I2
    pXsum(6, 6) = 1.0 - pI2;
    // I2 to RI
    pXsum(6, 7) = pI2;
    
    // RI to RI
    pXsum(7, 7) = 1.0;
    
    // DI to DI
    pXsum(8, 8) = 1.0;
    
    // H to H
    pXsum(9, 9) = 1.0 - pH;
    // H to RH
    pXsum(9, 10) = (1.0 - pHD) * pH;
    // H to DH
    pXsum(9, 11) = pHD * pH;
    
    // RH to RH
    pXsum(10, 10) = 1.0;
    
    // DH to DH
    pXsum(11, 11) = 1.0; 
    
    // print outputs if required
    if(print == 1) {
        Rprintf("pXsum:\n");
        for(int r = 0; r < 12; r++) {
            for(int c = 0; c < 12; c++) {
                Rprintf("%f ", pXsum(r, c));
            }
            Rprintf("\n");
        }
    }
    
    // initialise filtering probabilities
    arma::cube pXind;
    pXind.zeros(Nposs, t.n_elem, 12);
    
    // initialise backwards sampling probabilities
    arma::vec pXback;
    pXback.zeros(12);
    
    // loop over iterations
    for(iter = 0; iter < Niter; iter++) {
    
        // check for user interrupt
        R_CheckUserInterrupt();
    
        // loop over individuals
        for(i = 0; i < Nposs; i++) {
            
            // set initial state probability
            for(k = 0; k < 12; k++) pXind(i, 0, k) = 0.0;
            pXind(i, 0, 1) = pini;
            pXind(i, 0, 0) = 1.0 - pXind(i, 0, 1);
            
            // remove individual i from states
            // and cumulative states to get correct
            // powers in subsequent updates
            nXsum(0, Xind(i, 0))--;
            for(k = 0; k < t.n_elem; k++) nXcum(k, Xind(i, 0))--;
            for(k = 1; k < t.n_elem; k++) {
                if(Xind(i, k - 1) != Xind(i, k)) {
                    for(j = k; j < t.n_elem; j++) nXcum(j, Xind(i, k))--;
                }   
                nXsum(k, Xind(i, k))--;
            }
            
            // print outputs if required
            if(print == 1) {
                Rprintf("Xind -1:\n");
                for(int r = 0; r < t.n_elem; r++) Rprintf("%d ", Xind(i, r));
                Rprintf("\n");
                Rprintf("nXsum:\n");
                for(int r = 0; r < t.n_elem; r++) {
                    for(int c = 0; c < 12; c++) {
                        Rprintf("%d ", nXsum(r, c));
                    }
                    Rprintf("\n");
                }
                
                Rprintf("nXcum:\n");
                for(int r = 0; r < t.n_elem; r++) {
                    for(int c = 0; c < 12; c++) {
                        Rprintf("%d ", nXcum(r, c));
                    }
                    Rprintf("\n");
                }
            }
                    
            // update transmission probabilities with individual i removed     
            // S to S
            pXsum(0, 0) = exp(-beta * (nXsum(0, 4) + nXsum(0, 5) + nXsum(0, 6) + betaA * nXsum(0, 2)) / ((double) Npop));        
            // S to E
            pXsum(0, 1) = (1.0 - pXsum(0, 0)) + pini - (1.0 - pXsum(0, 0)) * pini;
            pXsum(0, 0) = 1.0 - pXsum(0, 1);
        
            // loop over possible states
            norm = 0.0;
            for(k = 0; k < 12; k++) {
                
                if(pXind(i, 0, k) > 0.0) {
                    
                    // set observation process
                    if((nXsum(0, 8) + (k == 8 ? 1:0)) == DI(0) && (nXsum(0, 11) + (k == 11 ? 1:0)) == DH(0)) {
                        // if infective state then update transmission probs
                        if(k == 4 || k == 5 || k == 6) {
                            // update transmission probabilities
                            // S to S
                            pXsum(0, 0) = exp(-beta * (1 + nXsum(0, 4) + nXsum(0, 5) + nXsum(0, 6) + betaA * nXsum(0, 2)) / ((double) Npop));        
                            // S to E
                            pXsum(0, 1) = (1.0 - pXsum(0, 0)) + pini - (1.0 - pXsum(0, 0)) * pini;
                            pXsum(0, 0) = 1.0 - pXsum(0, 1);
                        }
                        if(k == 2) {
                            // update transmission probabilities
                            // S to S
                            pXsum(0, 0) = exp(-beta * (nXsum(0, 4) + nXsum(0, 5) + nXsum(0, 6) + betaA * (1 + nXsum(0, 2))) / ((double) Npop));        
                            // S to E
                            pXsum(0, 1) = (1.0 - pXsum(0, 0)) + pini - (1.0 - pXsum(0, 0)) * pini;
                            pXsum(0, 0) = 1.0 - pXsum(0, 1);
                        }
                        
                        // now calculate other individual's transition probabilities
                            
                        // states: S, E, A, RA, P, I1, I2, RI, DI, H, RH, DH
                        //         0, 1, 2, 3,  4, 5,  6,  7,  8,  9, 10, 11
                        
                        // S to S
                        pXind(i, 0, k) *= pow(pXsum(0, 0), (double) nXsum(1, 0));
                        // S to E
                        pXind(i, 0, k) *= pow(pXsum(0, 1), (double) nXcum(1, 1) - nXcum(0, 1));
                        
                        // E to E
                        pXind(i, 0, k) *= pow(pXsum(1, 1), (double) nXsum(0, 1) - (nXcum(1, 2) - nXcum(0, 2)) - (nXcum(1, 4) - nXcum(0, 4))); 
                        // E to A
                        pXind(i, 0, k) *= pow(pXsum(1, 2), (double) nXcum(1, 2) - nXcum(0, 2));
                        // E to P
                        pXind(i, 0, k) *= pow(pXsum(1, 4), (double) nXcum(1, 4) - nXcum(0, 4));
                        
                        // A to A
                        pXind(i, 0, k) *= pow(pXsum(2, 2), (double) nXsum(0, 2) - (nXcum(1, 3) - nXcum(0, 3)));
                        // A to RA
                        pXind(i, 0, k) *= pow(pXsum(2, 3), (double) nXcum(1, 3) - nXcum(0, 3));
                        
                        // RA to RA
                        pXind(i, 0, k) *= pow(pXsum(3, 3), (double) nXsum(0, 3));
                        
                        // P to P
                        pXind(i, 0, k) *= pow(pXsum(4, 4), (double) nXsum(0, 4) - (nXcum(1, 5) - nXcum(0, 5)));
                        // P to I1
                        pXind(i, 0, k) *= pow(pXsum(4, 5), (double) nXcum(1, 5) - nXcum(0, 5));
                        
                        // I1 to I1
                        pXind(i, 0, k) *= pow(pXsum(5, 5), (double) nXsum(0, 5) - (nXcum(1, 6) - nXcum(0, 6)) - (nXcum(1, 8) - nXcum(0, 8)) - (nXcum(1, 9) - nXcum(0, 9)));
                        // I1 to I2
                        pXind(i, 0, k) *= pow(pXsum(5, 6), (double) nXcum(1, 6) - nXcum(0, 6));
                        // I1 to DI
                        pXind(i, 0, k) *= pow(pXsum(5, 8), (double) nXcum(1, 8) - nXcum(0, 8));
                        // I1 to H
                        pXind(i, 0, k) *= pow(pXsum(5, 9), (double) nXcum(1, 9) - nXcum(0, 9));
                        
                        // I2 to I2
                        pXind(i, 0, k) *= pow(pXsum(6, 6), (double) nXsum(0, 6) - (nXcum(1, 7) - nXcum(0, 7)));
                        // I2 to RI
                        pXind(i, 0, k) *= pow(pXsum(6, 7), (double) nXcum(1, 7) - nXcum(0, 7));
                        
                        // RI to RI
                        pXind(i, 0, k) *= pow(pXsum(7, 7), (double) nXsum(0, 7));
                        
                        // H to H
                        pXind(i, 0, k) *= pow(pXsum(9, 9), (double) nXsum(0, 9) - (nXcum(1, 10) - nXcum(0, 10)) - (nXcum(1, 11) - nXcum(0, 11)));
                        // H to RH
                        pXind(i, 0, k) *= pow(pXsum(9, 10), (double) nXcum(1, 10) - nXcum(0, 10));
                        // H to DH
                        pXind(i, 0, k) *= pow(pXsum(9, 11), (double) nXcum(1, 11) - nXcum(0, 11));
                        
                        // RH to RH
                        pXind(i, 0, k) *= pow(pXsum(10, 10), (double) nXsum(0, 10));
                        
                        // DH to DH
                        pXind(i, 0, k) *= pow(pXsum(11, 11), (double) nXsum(0, 11));
                    } else {
                        pXind(i, 0, k) = 0.0;
                    }
                }
                norm += pXind(i, 0, k);
            }
            
            // normalise probabilities
            for(k = 0; k < 12; k++) pXind(i, 0, k) /= norm;
        
            // print outputs if required
            if(print == 1) {
                Rprintf("pXind(%d, 0): ", i);
                for(int c = 0; c < 12; c++) {
                    Rprintf("%f ", pXind(i, 0, c));
                }
                Rprintf("\n");
            }
            
            // loop over other time points
            for(j = 1; j < (t.n_elem - 1); j++) {
            
                // update transmission probabilities
                // S to S
                pXsum(0, 0) = exp(-beta * (nXsum(j, 4) + nXsum(j, 5) + nXsum(j, 6) + betaA * (nXsum(j, 2))) / ((double) Npop));        
                // S to E
                pXsum(0, 1) = (1.0 - pXsum(0, 0)) + pini - (1.0 - pXsum(0, 0)) * pini;
                pXsum(0, 0) = 1.0 - pXsum(0, 1);
    
//                // print outputs if required
//                if(print == 1) {
//                    Rprintf("pXsum:\n");
//                    for(int r = 0; r < 12; r++) {
//                        for(int c = 0; c < 12; c++) {
//                            Rprintf("%f ", pXsum(r, c));
//                        }
//                        Rprintf("\n");
//                    }
//                }
            
                // loop over possible states to move to
                for(l = 0; l < 12; l++) {
                    // loop over possible states to move from
                    pXind(i, j, l) = 0.0;
                    for(k = 0; k < 12; k++) {
                        pXind(i, j, l) += pXind(i, j - 1, k) * pXsum(k, l);
                    }
                }
                
                // loop over possible states to move to
                norm = 0.0;
                for(k = 0; k < 12; k++) {

                    if(pXind(i, j, k) > 0.0) {
                        
                        // set observation process
                        if((nXsum(j, 8) + (k == 8 ? 1:0)) == DI(j) && (nXsum(j, 11) + (k == 11 ? 1:0)) == DH(j)) {
                            // if infective state then update transmission probs
                            if(k == 4 || k == 5 || k == 6) {
                                // update transmission probabilities
                                // S to S
                                pXsum(0, 0) = exp(-beta * (1 + nXsum(j, 4) + nXsum(j, 5) + nXsum(j, 6) + betaA * nXsum(j, 2)) / ((double) Npop));        
                                // S to E
                                pXsum(0, 1) = (1.0 - pXsum(0, 0)) + pini - (1.0 - pXsum(0, 0)) * pini;
                                pXsum(0, 0) = 1.0 - pXsum(0, 1);
                            }
                            if(k == 2) {
                                // update transmission probabilities
                                // S to S
                                pXsum(0, 0) = exp(-beta * (nXsum(j, 4) + nXsum(j, 5) + nXsum(j, 6) + betaA * (1 + nXsum(j, 2))) / ((double) Npop));        
                                // S to E
                                pXsum(0, 1) = (1.0 - pXsum(0, 0)) + pini - (1.0 - pXsum(0, 0)) * pini;
                                pXsum(0, 0) = 1.0 - pXsum(0, 1);
                            }
                            
                            // now calculate other individual's transition probabilities
                                
                            // states: S, E, A, RA, P, I1, I2, RI, DI, H, RH, DH
                            //         0, 1, 2, 3,  4, 5,  6,  7,  8,  9, 10, 11
                            
                            // S to S
                            pXind(i, j, k) *= pow(pXsum(0, 0), (double) nXsum(j, 0));
                            // S to E
                            pXind(i, j, k) *= pow(pXsum(0, 1), (double) nXcum(j, 1) - nXcum(j - 1, 1));
                            
                            // E to E
                            pXind(i, j, k) *= pow(pXsum(1, 1), (double) nXsum(j, 1) - (nXcum(j, 2) - nXcum(j - 1, 2)) - (nXcum(j, 4) - nXcum(j - 1, 4))); 
                            // E to A
                            pXind(i, j, k) *= pow(pXsum(1, 2), (double) nXcum(j, 2) - nXcum(j - 1, 2));
                            // E to P
                            pXind(i, j, k) *= pow(pXsum(1, 4), (double) nXcum(j, 4) - nXcum(j - 1, 4));
                            
                            // A to A
                            pXind(i, j, k) *= pow(pXsum(2, 2), (double) nXsum(j, 2) - (nXcum(j, 3) - nXcum(j - 1, 3)));
                            // A to RA
                            pXind(i, j, k) *= pow(pXsum(2, 3), (double) nXcum(j, 3) - nXcum(j - 1, 3));
                            
                            // RA to RA
                            pXind(i, j, k) *= pow(pXsum(3, 3), (double) nXsum(j, 3));
                            
                            // P to P
                            pXind(i, j, k) *= pow(pXsum(4, 4), (double) nXsum(j, 4) - (nXcum(j, 5) - nXcum(j - 1, 5)));
                            // P to I1
                            pXind(i, j, k) *= pow(pXsum(4, 5), (double) nXcum(j, 5) - nXcum(j - 1, 5));
                            
                            // I1 to I1
                            pXind(i, j, k) *= pow(pXsum(5, 5), (double) nXsum(j, 5) - (nXcum(j, 6) - nXcum(j - 1, 6)) - (nXcum(j, 8) - nXcum(j - 1, 8)) - (nXcum(j, 9) - nXcum(j - 1, 9)));
                            // I1 to I2
                            pXind(i, j, k) *= pow(pXsum(5, 6), (double) nXcum(j, 6) - nXcum(j - 1, 6));
                            // I1 to DI
                            pXind(i, j, k) *= pow(pXsum(5, 8), (double) nXcum(j, 8) - nXcum(j - 1, 8));
                            // I1 to H
                            pXind(i, j, k) *= pow(pXsum(5, 9), (double) nXcum(j, 9) - nXcum(j - 1, 9));
                            
                            // I2 to I2
                            pXind(i, j, k) *= pow(pXsum(6, 6), (double) nXsum(j, 6) - (nXcum(j, 7) - nXcum(j - 1, 7)));
                            // I2 to RI
                            pXind(i, j, k) *= pow(pXsum(6, 7), (double) nXcum(j, 7) - nXcum(j - 1, 7));
                            
                            // RI to RI
                            pXind(i, j, k) *= pow(pXsum(7, 7), (double) nXsum(j, 7));
                            
                            // H to H
                            pXind(i, j, k) *= pow(pXsum(9, 9), (double) nXsum(j, 9) - (nXcum(j, 10) - nXcum(j - 1, 10)) - (nXcum(j, 11) - nXcum(j - 1, 11)));
                            // H to RH
                            pXind(i, j, k) *= pow(pXsum(9, 10), (double) nXcum(j, 10) - nXcum(j - 1, 10));
                            // H to DH
                            pXind(i, j, k) *= pow(pXsum(9, 11), (double) nXcum(j, 11) - nXcum(j - 1, 11));
                            
                            // RH to RH
                            pXind(i, j, k) *= pow(pXsum(10, 10), (double) nXsum(j, 10));
                            
                            // DH to DH
                            pXind(i, j, k) *= pow(pXsum(11, 11), (double) nXsum(j, 11));
                        } else {
                            pXind(i, j, k) = 0.0;
                        }
                    }
                    norm += pXind(i, j, k);
                }
                
                // normalise probabilities
                for(k = 0; k < 12; k++) pXind(i, j, k) /= norm;
        
                // print outputs if required
                if(print == 1) {
                    Rprintf("pXind(%d, %d): ", i, j);
                    for(int c = 0; c < 12; c++) {
                        Rprintf("%f ", pXind(i, j, c));
                    }
                    Rprintf("\n");
                }
            }
            
            // update final time point
            
            // update transmission probabilities
            // S to S
            pXsum(0, 0) = exp(-beta * (nXsum(j, 4) + nXsum(j, 5) + nXsum(j, 6) + betaA * nXsum(j, 2)) / ((double) Npop));        
            // S to E
            pXsum(0, 1) = (1.0 - pXsum(0, 0)) + pini - (1.0 - pXsum(0, 0)) * pini;
            pXsum(0, 0) = 1.0 - pXsum(0, 1);
        
            // loop over possible states to move to
            for(l = 0; l < 12; l++) {
                // loop over possible states to move from
                pXind(i, j, l) = 0.0;
                for(k = 0; k < 12; k++) {
                    pXind(i, j, l) += pXind(i, j - 1, k) * pXsum(k, l);
                } 
            }
            
            // loop over possible states to move to
            norm = 0.0;
            for(k = 0; k < 12; k++) {
                if(pXind(i, j, k) > 0.0) {          
                    // set observation process
                    if((nXsum(j, 8) + (k == 8 ? 1:0)) != DI(j) || (nXsum(j, 11) + (k == 11 ? 1:0)) != DH(j)) {
                        pXind(i, j, k) = 0.0;
                    }
                }
                norm += pXind(i, j, k);
            }   
            // normalise probabilities
            for(k = 0; k < 12; k++) pXind(i, j, k) /= norm;
        
            // print outputs if required
            if(print == 1) {
                Rprintf("pXind(%d, %d): ", i, j);
                for(int c = 0; c < 12; c++) {
                    Rprintf("%f ", pXind(i, j, c));
                }
                Rprintf("\n");
            }
            
            // backwards sampling step
            u = R::runif(0.0, 1.0);
            k = 0;
            norm = pXind(i, t.n_elem - 1, k);
            while(u > norm) {
                k++;
                norm += pXind(i, t.n_elem - 1, k);
            }
            // set sample for individual i
            Xind(i, t.n_elem - 1) = k;
        
            // print outputs if required
            if(print == 1) {
                Rprintf("\npXback(%d): ", t.n_elem - 1);
                for(int c = 0; c < 12; c++) {
                    Rprintf("%f ", pXind(i, t.n_elem - 1, c));
                }
                Rprintf("\n");
            }
            
            // run through rest of time points
            for(j = (t.n_elem - 2); j >= 0; j--) {
            
                // update transmission probabilities
                // S to S
                pXsum(0, 0) = exp(-beta * (nXsum(j, 4) + nXsum(j, 5) + nXsum(j, 6) + betaA * nXsum(j, 2)) / ((double) Npop));        
                // S to E
                pXsum(0, 1) = (1.0 - pXsum(0, 0)) + pini - (1.0 - pXsum(0, 0)) * pini;
                pXsum(0, 0) = 1.0 - pXsum(0, 1);
            
                // calculate event probability
                norm = 0.0;
                for(k = 0; k < 12; k++) {
                    norm += pXsum(k, Xind(i, j + 1)) * pXind(i, j, k);
                }
                for(l = 0; l < 12; l++) {
                    pXback(l) = pXsum(l, Xind(i, j + 1)) * pXind(i, j, l);
                    pXback(l) /= norm;
                }
                
                // print outputs if required
                if(print == 1) {
                    Rprintf("pXback(%d): ", j);
                    for(int c = 0; c < 12; c++) {
                        Rprintf("%f ", pXback(c));
                    }
                    Rprintf("\n");
                }
                
                // sample state
                u = R::runif(0.0, 1.0);
                k = 0;
                norm = pXback(k);
                while(u > norm) {
                    k++;
                    norm += pXback(k);
                }
                // set sample for individual i
                Xind(i, j) = k;
            }
            
            // print outputs if required
            if(print == 1) {
                Rprintf("Xind(%d):\n", i);
                for(j = 0; j < t.n_elem; j++) Rprintf("%d ", Xind(i, j));
                Rprintf("\n");
            }
            
            // add individual i from states
            // and cumulative states to get correct
            // powers in subsequent updates
            nXsum(0, Xind(i, 0))++;
            for(k = 0; k < t.n_elem; k++) nXcum(k, Xind(i, 0))++;
            for(k = 1; k < t.n_elem; k++) {
                if(Xind(i, k - 1) != Xind(i, k)) {
                    for(j = k; j < t.n_elem; j++) nXcum(j, Xind(i, k))++;
                }   
                nXsum(k, Xind(i, k))++;
            }
            
            // print outputs if required
            if(print == 1) {
                Rprintf("nXsum:\n");
                for(int r = 0; r < t.n_elem; r++) {
                    for(int c = 0; c < 12; c++) {
                        Rprintf("%d ", nXsum(r, c));
                    }
                    Rprintf("\n");
                }
                
                Rprintf("nXcum:\n");
                for(int r = 0; r < t.n_elem; r++) {
                    for(int c = 0; c < 12; c++) {
                        Rprintf("%d ", nXcum(r, c));
                    }
                    Rprintf("\n");
                }
            }
                    
            // update transmission probabilities with individual i removed     
            // S to S
            pXsum(0, 0) = exp(-beta * (nXsum(0, 4) + nXsum(0, 5) + nXsum(0, 6) + betaA * nXsum(0, 2)) / ((double) Npop));        
            // S to E
            pXsum(0, 1) = (1.0 - pXsum(0, 0)) + pini - (1.0 - pXsum(0, 0)) * pini;
            pXsum(0, 0) = 1.0 - pXsum(0, 1);
        }
        
        // save output
        for(int r = 0; r < t.n_elem; r++) {
            for(int c = 0; c < 12; c++) {
                output(iter, r * 12 + c) = nXcum(r, c);
            }
        }
        if(iter % 100 == 0) Rprintf("i = %d\n", iter);
    }
            
    // return alphas
    return output;
}
