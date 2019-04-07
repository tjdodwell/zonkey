#ifndef ZONKEY_MCMC_MULTI_CHAIN_HH
#define ZONKEY_MCMC_MULTI_CHAIN_HH



namespace Zonkey {
  namespace MCMC {

template<class Link,class Chain>
class MultiChain{


  public:

    MultiChain(){  }

    void setLevel(int numLevels){ multiChain.resize(numLevels); }

    Chain getLevel(int l){ return multiChain[l]; }

    Link& back(int l){ return multiChain[l].back(); }

    void addLink(Link& newLink, int level){  multiChain[level].push_back(newLink);   }

    Eigen::VectorXd getTheta(int level, int index){
      auto thisChain = multiChain[level];
      return thisChain[index].getTheta();
    }

    int size(int level){ return multiChain[level].size(); }

  private:

    std::vector<Chain> multiChain;

};

}
}

#endif // end ZONKEY_MCMC_MULTI_CHAIN_HH
