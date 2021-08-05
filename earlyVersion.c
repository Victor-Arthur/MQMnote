#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
/*Struct ket
*k indicates the coefficient.
*dimension indicates the dimensions of the bases.
*base indicates the base vector by using binary method.
*example:
*3(10)=11(2), with dimension of 2 and k of 0.25, this ket is like:
*1/4|1,1>
*/
struct ket {
    double k = 1;
    unsigned long base = 0;
    unsigned long dimension = 1;
};
/*
*Construct a ket.
*/
void getKet(double k, unsigned long base, unsigned long n, ket *r) {
    r->k = k;
    r->base = base;
    r->dimension = n;
}
/*
*Construct a copy from a known ket.
*/
void copy(ket *s1, ket *s2) {
    s2->k = s1->k;
    s2->base = s1->base;
    s2->dimension = s1->dimension;
}
/*
*X component of Hamiltonian, indicates Sx1Sx2.
*/
void spinx(ket *s, unsigned long n1, unsigned long n2) {
    s->k = (s->k) / 4;
    unsigned long n = s->base;
    n = (1 << n2) ^ ((1 << n1) ^ n);
    s->base = n;
}
/*
*The coefficient after Sy twice will still be a real number.
*/
void spiny(ket *s, unsigned long n1, unsigned long n2) {
    s->k = -(s->k) / 4;
    unsigned long n = s->base;
    if (((1 << n1) & n) != 0) {
        s->k = -s->k;
    }
    if (((1 << n2) & n) != 0) {
        s->k = -s->k;
    }
    n = (1 << n2) ^ ((1 << n1) ^ n);
    s->base = n;
}

void spinz(ket *s, unsigned long n1, unsigned long n2) {
    s->k = (s->k) / 4;
    unsigned long n = s->base;
    if (((1 << n1) & n) != 0) {
        s->k = -s->k;
    }
    if (((1 << n2) & n) != 0) {
        s->k = -s->k;
    }
    s->base = n;
}

double innerProduct(ket *s1, ket *s2) {
    if (s1->base != s2->base) {
        return 0;
    }
    return s1->k * s2->k;
}

double hamiltonian(ket *s1, ket *s2) {
    unsigned long n = s1->dimension;
    double result = 0;
    for (int i = 0; i < n - 1; i++) {
        int j = i + 1;
        ket x;
        ket y;
        ket z;
        copy(s2, &x);
        copy(s2, &y);
        copy(s2, &z);
        spinx(&x, i, j);
        spiny(&y, i, j);
        spinz(&z, i, j);
        result = result + innerProduct(s1, &x);
        result = result + innerProduct(s1, &y);
        result = result + innerProduct(s1, &z);
    }
    return result;
}

size_t getCombination(size_t n, size_t k){
    if(k > n|k < 0)
        return 0;
    size_t l = k<(n-k)?k:(n-k);
    size_t up = 1, lo = 1;
    for (int i = 0 ; i < l ; i++){
        up = up * (n-i);
        lo = lo * (1+i);
    }
    return up/lo;
}

size_t getLocation(size_t n, size_t k){
    size_t c = 0;
    for (unsigned long i = 0; i < n; i++){
        unsigned long  r = 1 << i;
        if ((k&r)>0)
            c++;
    }
    return c;
}
/*
*Construct the M0 matrix.
*/
gsl_matrix* getSimplifiedHamiltonian(long n){
    gsl_vector_ulong *v = gsl_vector_ulong_alloc(1 << (n-1));
    unsigned long c = 0;
    unsigned long i = (1 << ((n+1)/2))-1;

    while(c < getCombination(n,(n+1)/2)){
        if (getLocation(n, i)==(n+1)/2){
            gsl_vector_ulong_set(v,c,i);
            c++;
        }
        i++;
    }
    static gsl_matrix *r = gsl_matrix_calloc(c,c);
    for (i = 0; i < c; i++) {
        for (unsigned long j = 0; j < c; j++) {
            ket s1;
            ket s2;
            getKet(1, gsl_vector_ulong_get(v,i), n, &s1);
            getKet(1, gsl_vector_ulong_get(v,j), n, &s2);
            double a = hamiltonian(&s1, &s2);
            gsl_matrix_set(r,i,j,a);
        }
    }
    gsl_vector_ulong_free(v);
    return r;
}