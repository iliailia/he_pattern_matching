/* Copyright (C) 2019 IBM Corp.
 * This program is Licensed under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 *   http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License. See accompanying LICENSE file.
 */
#include <iostream>
#include <time.h>
#include <random>

#include <helib/helib.h>
#include <helib/debugging.h>

using namespace std;
using namespace NTL;
using namespace helib;

double intlog(unsigned int base, unsigned int input)
{
  return floor(log2(input)/log2(base));
}

void convertToBaset(vector<long>& decomp, unsigned int input, unsigned int base, int nslots)
{
  decomp.clear();
  decomp.resize(nslots,0);
  int power = static_cast<int>(intlog(base, input)) + 1;
  if (power > nslots)
  {
    cout << "Input character is too big to be converted" << endl;
    exit(1);
  }
  unsigned int rest = input;
  unsigned int coeff;

  int i = 0;
  while(i < power)
    {
      coeff = rest % base;
      decomp[i] = static_cast<long>(coeff);
      rest = (rest - coeff) / base;
      i++;
    }
}


void rotate_and_add(helib::Ctxt& x, int nb_slots, long pat_len, const helib::EncryptedArray& ea)
{
  if (nb_slots == 0)
  {
    std::vector<long> ptxt1(nb_slots,1);
    ea.encrypt(x,x.getPubKey(),ptxt1);
  }
  helib::Ctxt y(x.getPubKey());
  std::vector<long> ptxt0(nb_slots,0);
  ea.encrypt(y,x.getPubKey(),ptxt0);
  int rot_x = 1;
  //rot_y = -2^K, where K is floor(log_2(pat_len))
  int rot_y = pat_len;
  helib::Ctxt temp(x.getPubKey());
  //std::vector<long> decryptedx(nb_slots);
  //std::vector<long> decryptedy(nb_slots);
  int tmp_len = pat_len;
  while (tmp_len > 1)
    {
      if (tmp_len % 2 == 0)
        {
          temp = x;
          ea.rotate(temp, -rot_x);
          x += temp;
          rot_x *=2;
          tmp_len = tmp_len/2;
        }
      else //if (nb_slots % 2 == 1)
        {
          rot_y -= rot_x;
          temp = x;
          ea.rotate(temp, -rot_y);
          y += temp;
          temp = x;
          ea.rotate(temp, -rot_x);
          x += temp;
          rot_x *= 2;
          tmp_len = (tmp_len - 1)/2;
        }
    }
  x +=y;
}


void prob_eq_circuit(Ctxt& ctxt_res, Ctxt& ctxtx, helib::Ctxt& ctxty, int nslots, long pat_len, int p, int ord_p, bool wildcard, SecKey& secret_key, const helib::EncryptedArray& ea)
{
  FHE_NTIMER_START(EqualityCircuit);
  //first compute the difference slotwise
  ctxt_res = ctxtx;
  ctxt_res -= ctxty;
  //multiply by a random polynomial
  Ptxt<BGV> poly_r(ea.getContext());
  poly_r.random();

  ctxt_res.multByConstant(poly_r);

  // additional multiplication when using wildcards
  if(wildcard)
    ctxt_res.multiplyBy(ctxty);

  vector<ZZX> decrypted(nslots);
  ea.decrypt(ctxt_res, secret_key, decrypted);

  //rotate and add
  rotate_and_add(ctxt_res, nslots, pat_len, ea);

  mapTo01(ea, ctxt_res);

  std::vector<long> ptxtones(nslots,1);
  NTL::ZZX poly_ones;
  ea.encode(poly_ones,ptxtones);
  ctxt_res.negate();
  ctxt_res.addConstant(poly_ones,1);
  FHE_NTIMER_STOP(EqualityCircuit);
}

// main call examples
// ./pattern_matching p m q pattern_len experiment_runs wildcard_bool
// ./pattern_matching 7 21177 290 5 100 0
// ./pattern_matching 7 21177 320 7 100 1
// ./pattern_matching 17 18913 330 3 90 0
// ./pattern_matching 17 18913 360 3 90 1

int main(int argc, char *argv[]) {
  // initialize the random generator
  random_device rd;
  mt19937 eng(rd());
  uniform_int_distribution<unsigned int> distr_u;
  uniform_int_distribution<int> distr_i; 

  // Wildcard. If true, then the wildcard character is encoded by 0.
  bool wildcard;
  if (atoi(argv[6]) == 0)
  {
    wildcard = false;
  }
  else if (atoi(argv[6]) == 1)
  {
    wildcard = true;
  }
  else
  {
    cout << "Wildcard parameter should be either 0 or 1" << endl;
    return 1;
  }

  // Plaintext prime modulus
  unsigned int p = atol(argv[1]);
  // Cyclotomic polynomial - defines phi(m)
  unsigned long m = atol(argv[2]);
  // Hensel lifting (default = 1)
  unsigned long r = 1;
  // Number of ciphertext prime bits in the modulus chain = depth of the computation
  unsigned long nb_primes = atol(argv[3]);
  // Number of columns of Key-Switching matix (default = 2 or 3)
  unsigned long c = 3;
  std::cout << "Initialising context object..." << std::endl;
  // Intialise context
  Context context(m, p, r);
  context.scale = 6;
  // Modify the context, adding primes to the modulus chain
  cout  << "Building modulus chain..." << endl;
  buildModChain(context, nb_primes, c);

  // Print the context
  context.zMStar.printout();
  cout << endl;

  //determine the order of p in (Z/mZ)*
  long ord_p = context.zMStar.getOrdP();

  // Print the security level
  cout << "Q size: " << context.logOfProduct(context.ctxtPrimes)/log(2.0) << endl;
  cout << "Q*P size: " << context.logOfProduct(context.fullPrimes())/log(2.0) << endl;
  cout << "Security: " << context.securityLevel() << endl;

  // Secret key management
  cout << "Creating secret key..." << endl;
  // Create a secret key associated with the context
  SecKey secret_key(context);
  // Generate the secret key
  secret_key.GenSecKey();
  cout << "Generating key-switching matrices..." << endl;
  // Compute key-switching matrices that we need
  add1DMatrices(secret_key);
  addFrbMatrices(secret_key);

  // Public key management
  // Set the secret key (upcast: SecKey is a subclass of PubKey)
  const PubKey& public_key = secret_key;

  // Get the EncryptedArray of the context
  const EncryptedArray& ea = *(context.ea);

  // Get the number of slot (phi(m))
  long nslots = ea.size();
  cout << "Number of slots: " << nslots << endl;
  cout << "Extension degree of a slot:  " << ord_p << endl;

  //pattern length
  long pat_len = atol(argv[4]);
  //pattern copies in one ciphertext
  long pat_copies = nslots/pat_len;

  //timers
  setTimersOn();

  //repeat experiments several times
  int runs = atoi(argv[5]);
  long min_capacity = 1000;
  long capacity;
  for (int run = 0; run < runs; run++)
  {
    vector<ZZX> expected_result(nslots);
    vector<ZZX> decrypted(nslots);

    // Create the plaintext polynomials for the text and for the pattern
    vector<ZZX> pol_txt(nslots);
    vector<ZZX> pol_pat(nslots);
    
    // Generates copies of a random pattern of length len with characters in the set {0,..., 2^32-1} and encode each character into F_(p^{ord_p})
    unsigned int input_pat_coef;
    ZZX pol_pat_slot;

    for (int i = 0; i < pat_len; i++)
    {
      if (wildcard && (i > 0))
      {
        // generate a wildcard character with probability 0.33
        input_pat_coef = distr_u(eng) * static_cast<unsigned int>(distr_i(eng) % 3);
      }
      else
      {
        input_pat_coef = distr_u(eng);
      }
      vector<long> decomp_char;
      convertToBaset(decomp_char, input_pat_coef, p, ord_p);
      for (int j = 0; j < ord_p; j++)
      {
        SetCoeff(pol_pat_slot, j, decomp_char[j]);
      }
      for (int j = 0; j < pat_copies; j++)
      {
        pol_pat[i + j * pat_len] = pol_pat_slot;
      }
    }
    
    /*
    cout << "Input pattern: " << endl;
    for (int i =0; i< nslots; i++)
    {
        printZZX(cout, pol_pat[i], ord_p);
        cout << "Is zero?" << IsZero(pol_pat[i]);
        cout << endl;
    }
    */

    // text generation
    // indicates whether a snippet of the text will be a copy of the pattern
    int is_equal;
    ZZX pol_txt_slot;
    unsigned int input_txt_coef;
    int iChar = 0;

    while (iChar < nslots)
    {
      //randomly decide whether the text substring is equal to the pattern
      if (iChar <= nslots - pat_len)
      {
        is_equal = distr_i(eng) % 20;
        is_equal = (is_equal > 0) ? 0 : 1;
      }
      else
      {
        is_equal = 0;
      }

      expected_result[iChar] = ZZX(INIT_MONO, 0, is_equal);

      if (is_equal)
      {
        for (int i = 0; i < pat_len; i++)
        {
          if (IsZero(pol_pat[i]) && wildcard)
          {
            input_txt_coef = distr_u(eng)|1u;
            vector<long> decomp_char;
            convertToBaset(decomp_char, input_txt_coef,p, ord_p);
            for (int j = 0; j < ord_p; j++)
            {
              SetCoeff(pol_txt_slot, j, decomp_char[j]);
            }
            pol_txt[iChar + i] = pol_txt_slot;
          }
          else
            pol_txt[iChar + i] = pol_pat[i];
        }

        iChar += pat_len;
      }
      else
      {
        if(wildcard)
        {
          input_txt_coef = distr_u(eng)|1u;
        }
        else
        {
          input_txt_coef = distr_u(eng);
        }
        vector<long> decomp_char;
        convertToBaset(decomp_char, input_txt_coef,p, ord_p);
        for (int i = 0; i < ord_p; i++)
        {
          SetCoeff(pol_txt_slot, i, decomp_char[i]);
        }
        if (pol_txt_slot != pol_pat[0])
        {
          pol_txt[iChar] = pol_txt_slot;
          iChar++;
        }  
      }
    }
    
    /*
    cout << "Input: " << endl;
    for (int i =0; i< nslots; i++)
    {
        printZZX(cout, pol_pat[i], ord_p);
        cout << '\t';
        printZZX(cout, pol_txt[i], ord_p);
        cout << endl;
    }
    
    cout << "Expected results:" << endl;
    for(int i=0; i < expected_result.size(); i++)
    {
      printZZX(cout, expected_result[i], ord_p);
    }
    cout << endl;
    */

    Ctxt ctxt_pat(public_key);
    Ctxt ctxt_txt(public_key);
    ea.encrypt(ctxt_pat, public_key, pol_pat);
    ea.encrypt(ctxt_txt, public_key, pol_txt);

    //results
    Ctxt ctxt_res(ZeroCtxtLike, ctxt_txt);    
    //compute the equality
    printf("Run %d started\n", run);

    FHE_NTIMER_START(PatternMatching);

    for (int iRot = 0; iRot < pat_len; iRot++)
    {
      //printf("Rotation %d\n", iRot);
      Ctxt ctxt_tmp(public_key);
      Ctxt ctxt_pat_rot(ctxt_pat);

      if (iRot > 0)
      {
        ea.rotate(ctxt_pat_rot, iRot);
      }

      prob_eq_circuit(ctxt_tmp, ctxt_txt, ctxt_pat_rot, nslots, pat_len, p, ord_p, wildcard, secret_key, ea);

      // select correct slots and zeroized the rest
      vector<long> sel_vec(nslots,0);
      //index of the first slot containing the equality function result
      int iSlot = iRot;
      while(iSlot < nslots)
      {
        if (iSlot <= nslots - pat_len)
        {
          sel_vec[iSlot] = 1;
          iSlot += pat_len;
        }
        else
        {
          sel_vec[iSlot] = 0;
          iSlot++;
        }
      }
      //cout << "Selector slots: " << helib::vecToStr(sel_vec) << endl;
      ZZX sel_poly;
      ea.encode(sel_poly, sel_vec);
      ctxt_tmp.multByConstant(sel_poly);

      ctxt_res += ctxt_tmp;
    }
    printNamedTimer(cout, "EqualityCircuit");
    FHE_NTIMER_STOP(PatternMatching);
    printNamedTimer(cout, "PatternMatching");

    // remove the line below if it gives bizarre results 
    ctxt_res.cleanUp();
    capacity = ctxt_res.bitCapacity();
    cout << "Final capacity: " << capacity << endl;
    if (capacity < min_capacity)
      min_capacity = capacity;
    cout << "Min. capacity: " << min_capacity << endl;
    cout << "Final size: " << ctxt_res.logOfPrimeSet()/log(2.0) << endl;
    ea.decrypt(ctxt_res, secret_key, decrypted);

    for(int i = 0; i < nslots; i++)
    {
      
      if (decrypted[i] != expected_result[i])
      {
        printf("Slot %d: ", i);
        printZZX(cout, decrypted[i], ord_p);
        cout << endl;
        cout << "Failure" << endl;
        return 1;
      }
    }
  }

  return 0;
}
