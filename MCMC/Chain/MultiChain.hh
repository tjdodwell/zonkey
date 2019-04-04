#ifndef ZONKEY_MCMC_MULTI_CHAIN_HH
#define ZONKEY_MCMC_MULTI_CHAIN_HH


template<class Chain>
class MultilevelChain{


  public:

    typedef typename Chain::theLink Link;

    MultilevelChain(int numLevels){

      multiChain.resize(numLevels);

    }

    Chain getLevel(int l){ return multiChain[l]; }

    theLink& back(int l){ return multiChain[l].back(); }

  private:

    std::vector<Chain> multiChain;





};


#endif // end ZONKEY_MCMC_MULTI_CHAIN_HH
