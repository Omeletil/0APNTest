#include <stdio.h>
#include <gmp.h>

/* Generate representatives from all cyclotomic cosets in dimension n */
int main(int argc, char ** argv) {
  if (argc != 2) {
    printf("Usage: %s DIMENSION\n", argv[0]);
    return 1;
  }

  size_t n;
  sscanf(argv[1], "%lu", &n);
  
  if (n > 63) {
    printf("Error: DIMENSION >= 64\n");
    return 2;
  }
  
  size_t e0 = 1 + 0 * 2;
  size_t step = 1 * 2;

  size_t m = (1L << n) - 1;
  size_t emax = (1L << (n - 1)) - 1;
  
  mpz_t gmpM;
  mpz_init_set_ui(gmpM, m);

  mpz_t gmpR;
  mpz_init(gmpR);

  /* Go through all odd exponents */
  if (n & 1) {
    mpz_t gmpRI;
    mpz_init(gmpRI);
  
    for(size_t e = e0; e <= emax; e += step) {
      /* Check GCD */
      mpz_set_ui(gmpR, e);

      size_t gcd = mpz_gcd_ui((mpz_ptr) 0, gmpR, m);
      if(gcd != 1) {
        continue;
      }
      
      mpz_invert(gmpRI, gmpR, gmpM);
      
      size_t r = e;
      size_t ri = mpz_get_ui(gmpRI);

      /* Go through the cyclotomic coset and search for a smaller representative */
      for(size_t k = 0; k < n; ++k) {
        r = (r << 1) % m;
        ri = (ri << 1) % m;
        if(r < e || ri < e) {
          /* We have found something smaller; abort */
          goto next_e_odd;
        }
      }
      
      printf("%lu ", e);
      
      next_e_odd: continue;
    }
  } else {
    for(size_t e = e0; e <= emax; e += step) {
      /* Check GCD */
      mpz_set_ui(gmpR, e);

      size_t gcd = mpz_gcd_ui((mpz_ptr) 0, gmpR, m);
      if(gcd != 3) {
        continue;
      }
      
      size_t r = e;

      /* Go through the cyclotomic coset and search for a smaller representative */
      for(size_t k = 0; k < n; ++k) {
        r = (r << 1) % m;
        if(r < e) {
          /* We have found something smaller; abort */
          goto next_e_even;
        }
      }
      
      printf("%lu ", e);
      
      next_e_even: continue;
    }
  }

  printf("\n");
  return 0;
}
