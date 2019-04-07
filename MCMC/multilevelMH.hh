


#include <iostream>
#include <cstdio>
#include <ctime>


namespace Zonkey {
  namespace MCMC {


template<class ForwardModel, class Hierarchy, class CoarseProposal>
class MultiLevelMetropolisHastings{

  public:

    MultiLevelMetropolisHastings(int L_, ForwardModel& F_, CoarseProposal& coarseProposal_):
      L(L_),
      F(F_),
      coarseProposal(coarseProposal_){

        hierarchy.resize(L + 1);

        // Setup a hierarchy of Multilevel chains
        for (int i = 0; i < L + 1; i++){
          hierarchy[i].setLevel(i);
        }

        isBurnt_in.resize(L + 1);

        // Setup initial Samples fro each level
        NStar.resize(L + 1);
        NStar[0] = 1000;
        for (int i = 1; i < L + 1; i++){
          NStar[i] = std::max(NStar[i-1]/2,10);
        }

        //

        payOff.resize(L + 1);

    }

    void ApplyGreedy(){
      int ell = 0;
      // Coarse level first
      runCoarse(ell,Nstar[0]);
      Neff[0] = hierarchy[0].getEss(0)
      ACT[0] = Neff[0] / Nstar[0]; // Estimated Auto Correlation Time (ACT) for Coarse Chain
      mean[0] = hierarchy[0].getMean();
      variance[0] = hierarchy[0].getVar();

      samplingError[0] = variance * variance / Neff[0];


      double totalVariance = 0.0;


      for (int l = 0; l < L + 1; l++){

          // Compute some initial samples

          std::clock_t start;

          start = std::clock();

          computeSamples(l,Nstar[l]);

          Neff[l] = hierarchy[l].getEss(l)
          ACT[l] = Neff[l] / Nstar[l]; // Estimated Auto Correlation Time (ACT) for Coarse Chain

          mean[l] = hierarchy[l].getMean();
          Var[l]  = hierarchy[l].getVar();

          samplingError[l] = Var * Var / Neff[l];

          time_per_sample[l] = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

          payOff[l] = samplingError[l] / time_per_sample[l];

          totalVariance = 0.0;
          for (int i = 0; i < l + 1; i++){
            totalVariance += samplingError[i];
          }

          while (totalVariance > tolerance){

            // Find level less than l with greatest payoff

            theLevel = findMaxPayOff(l);

            // Double Samples on Level "theLevel"

            start = std::clock();

            int newSamples = hierarchy[theLevel].size(theLevel);

            computeSamples(theLevel,newSamples);

            time_per_sample[l] = (std::clock() - start ) / ((double) CLOCKS_PER_SEC); 







          }







      } // For each level





    }

    std::vector<int> howManySample(int level, int numSamples){
      // Returns how many samples require across all multichains
      std::vector<int> N(level + 1);
      N[level] = numSamples;
      for (int j = 1; j < level + 1; j++){
        N[level - j] = numSamples;
        for (int i = 1; i < j + 1; i++){
          N[level - j] *= ACT[level - i];
        }
      }
      return N;
    }

    void burnin(int level, int N){

      std::vector<int> numSamples = howManySamples(level,N);

      // Burning hierarchy[level]
      for (int i = 0; i < level + 1; i++){

          if(i == 0){ // This is i = 0, just a normal chain
            runCoarse(level,numSamples[0]);
          }
          else{
            run(i,level,numSamples[i]);
          }
      }

      isBurnt_in[level] = true;

       // Need to record where burn is.


    } // end burnin


    void inline runCoarse(int l, int numSamples){
      if(hierarchy[l].size(0) < 1){ // If this is the first sample
        Eigen::VectorXd theta_fP = F.samplePrior(0); // Sample prior
        Link firstPoint(theta_fP);
        hierarchy[l].addLink(firstPoint,0);
        F.apply(0,hierarchy[l].back(0));
        numSamples -= 1;
      }
      for (int i = 0; i < numSamples; i++){ // For each sample
        auto theta_p = coarseProposal.apply(hierarchy[l].back(0)); // Make a proposal
        F.apply(0,theta_p);  // Apply forward Model
        // Accept / Reject Step
        bool accept = coarseProposal.acceptReject(hierarchy[l].back(0),  theta_p);
        if(accept){ theta_p.setAccepted(1); hierarchy[l].addLink(theta_p,0);
        }
        else {  hierarchy[l].addLink(hierarchy[l].back(0),0);  } // accept / reject
      } // end for each sample
    } // end run coarseU

    void inline run(int l, int level, int numSamples){

      if(hierarchy[level].size(l) < 1){
        Eigen::VectorXd thetaC = hierarchy[level].getTheta(l-1,0);
        // Need to propose extra modes
        Link firstPoint(theta_p);
        hierarchy[level].addLink(firstPoint,l);
        F.apply(l,hierarchy[level].back(l));
        numSamples -= 1;
      }
      for (int i = 0; i < numSamples; i++){ // For each sample

        auto theta_p_C = hierarchy[level].getTheta(l-1,??);




      }






    }











  private:

    std::vector<Zonkey::MCMC:MultiChain> hierarchy;

    std::vector<bool> isBurnt_in;

    std::vector<int> burnin;

    std::vector<int> ACT;
    std::vector<double> mean, variance, payOff;

    ForwardModel& F_;

    CoarseProposal& coarseProposal;

};

}
}
