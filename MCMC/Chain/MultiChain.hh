#ifndef ZONKEY_MCMC_MULTI_CHAIN_HH
#define ZONKEY_MCMC_MULTI_CHAIN_HH



namespace Zonkey {
  namespace MCMC {

template<class Chain>
class MultiChain{


  public:

    typedef typename Chain::theLink Link;

    MultiChain(){  }

    void setLevel(int numLevels){ multiChain.resize(numLevels); }

    Chain getLevel(int l){ return multiChain[l]; }

    theLink& back(int l){ return multiChain[l].back(); }

    void addLink(theLink& newLink, int level){  multiChain[level].push_back(newLink);   }

    Eigen::VectorXd getTheta(int level, int index){
      auto thisChain = multiChain[level]
      return thisChain[index].getTheta();
    }

    int size(int level){ return multichain[level].size() }

  private:

    std::vector<Chain> multiChain;

};

}
}

#endif // end ZONKEY_MCMC_MULTI_CHAIN_HH
