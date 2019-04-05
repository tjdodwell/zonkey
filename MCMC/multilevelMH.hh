





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

    }

    std::vector<int> howManySample(int level, int numSamples){
      // Initiating getting a new proposal on level
      std::vector<int> N(level + 1);
      N[level] = numSamples;
      for (int i = 0; i < level; i++){
        N[i] = numSamples;
        for (int j = 0; j < level + 1; j++){
          N[i] *= ess[level-j];
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

    std::vector<int> ess;

    ForwardModel& F_;

    CoarseProposal& coarseProposal;

};

}
}
