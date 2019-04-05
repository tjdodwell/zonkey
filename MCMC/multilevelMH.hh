


std::vector<Zonkey::MCMC::MultiChain> hierarchy(L + 1);


namespace Zonkey {
  namespace MCMC {


template<class Hierarchy>
class MultiLevelMetropolisHastings{

  public:

    MultiLevelMetropolisHastings(Hierarchy& hierarchy_):
      hierarchy(hierarchy_){

        // Setup a hierarchy of Multilevel chains
        for (int i = 0; i < L + 1; i++){
          hierarchy[i].setLevel(i);
        }

    }

    void burnin(int level, int est_ACT, int factor = 10){






      
    }







  private:

    Hierarhcy hierarchy;

};

}
}
