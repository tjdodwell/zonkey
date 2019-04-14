//
//  KLFunctions.h
//
//
//  Created by Tim Dodwell on 01/12/2015.
//
//

#ifndef KLFunctions_hh
#define KLFunctions_hh

// ------------------------------------------------------------------------------------
// Define K-L Expansion
// ------------------------------------------------------------------------------------


template<class X>
void inline evaluate1DeigenValues(double ell, X& lambda, X& freq){
  //  returns eigenvalues for KL expansion.
    double c = 1 / ell;
    for (int i = 0; i < lambda.size(); i++){
        lambda[i] = 2.0 * c /(freq[i] * freq[i] + c * c);
    }
}

template<class X,class Y>
void bubbleSort(X& index, Y& lam)
{
    const int N = lam.size();

    for (int i = 0; i < N; i++){ index[i] = i; }

    int sw = 1; double tmpreal;
    int tmpint;

    while (sw == 1)
    {
        sw = 0;
        for(int i = 0; i < N-1; i++)
        {
            if (lam[i] < lam[i+1])
            {
                tmpreal = lam[i+1];
                tmpint = index[i+1];
                lam[i+1] = lam[i];
                index[i+1] = index[i];
                lam[i] = tmpreal;
                index[i] = tmpint;
                sw = 1;
            }
        }
    }
}

void inline construct3DeigenValues(std::vector<double>& lam1D, std::vector<double>& lambda, std::vector<std::vector<int>>& new_index){

    const int dim = 3;

    const int N = lam1D.size(); // Number of 1D eigenvalues

    

    std::vector<std::vector<int>> index(dim);

    for (int i = 0; i < dim; i++){
      index[i].resize(std::pow(N,dim));
      new_index[i].resize(std::pow(N,dim));
    }

    std::cout << "lam1D = " << lam1D.size() << std::endl;

    lambda.resize(std::pow(N,dim));
    std::vector<int> ind(std::pow(N,dim));

    int counter = 0;
    for (int i = 0; i < N; i++){
      for (int j = 0; j < N; j++){
          for (int k = 0; k < N; k++){
            lambda[counter] = lam1D[i] * lam1D[j] * lam1D[k];
            ind[counter] = counter;
            index[0][counter] = i;  index[1][counter] = j;  index[2][counter] = k;
            counter += 1;
            }
        }
    }

    bubbleSort(ind,lambda);

    for (int i = 0; i < lambda.size(); i++){
        new_index[0][i] = index[0][ind[i]];
        new_index[1][i] = index[1][ind[i]];
        new_index[2][i] = index[2][ind[i]];
    }

}

double res(double x, double ell){
  double c = 1.0 / ell;
  double g = std::tan(x) - 2.0 * c * x / (x * x - c * c);
  return g;
}


double findRoot(double ell, double a, double b){
  double fa, fb, x, fx, error, m;
  error = 1.0;
  while (error > 1e-6){
    fa = res(a,ell);
    fb = res(b,ell);
    m = (fb - fa) / (b - a);
    x = a - fa / m;
    fx = res(x,ell);
    if (((fa < 0) & (fx < 0)) | ((fa > 0) & (fx > 0))) { a = x; }
    else { b = x; }
    error = std::abs(fx);
  }
  return x;
}

void rootFinder(int M, double ell, std::vector<double>& answer){
double c = 1.0 / ell;
std::vector<double> freq(M+2);
// For all intervals

int m = -1;
for (int i = 0; i < M + 1; i++){

//  std::cout << "Root i = " << i << std::endl;

  double w_min = (i - 0.4999) * M_PI;
  double w_max = (i + 0.4999) * M_PI;

//  std::cout << "w_min = " << w_min << std::endl;
  //std::cout << "w_max = " << w_max << std::endl;
  if ((w_min <= c) && (w_max >= c)){

    // If not first interval look for solution near left boundary
    if (w_min > 0.0){
      m += 1;
      freq[m] = findRoot(ell,w_min,0.5*(c+w_min));
    }
    // Always look for solution near right boundary
    m += 1;
    freq[m] = findRoot(ell,0.5*(c + w_max),w_max);
  }
  else{
    m += 1;
    freq[m] = findRoot(ell,w_min,w_max);
  }

}

for (int i = 0; i < M; i++){
  answer[i] = freq[i+1];
}



}



#endif /* KLFunctions_h */
