#ifndef ZONKEY_MCMC_SINGLE_CHAIN_HH
#define ZONKEY_MCMC_SINGLE_CHAIN_HH

#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>

#include <algorithm>

// basic file operations
#include <iostream>
#include <fstream>
#include <string>


#include "getEss.hh"

using namespace Eigen;
using namespace std;

#include <vector>

namespace Zonkey {
  namespace MCMC {

  template<class Link>
  class SingleChain{

  public:

      SingleChain(){ }

      SingleChain(std::string fileName){

        readChainFromFile(fileName);

        std::cout << "Chain built from file = " << fileName << std::endl;


      }


      void addLink(Link& newLink, int accept = 0){
        theChain.push_back(newLink);
        theChain.back().setAccepted(accept);
      }

      Link& back(){return theChain.back(); }

      Link& operator[] (const int index){
        return theChain[index];
      }

      void inline resize(int N){  theChain.resize(N); }

      int inline size(){ return theChain.size(); }

      Eigen::VectorXd inline mean(bool param = false){ 
        int numParam = theChain[0].size(param); // Number of Parameters
        Eigen::VectorXd mu(numParam); 
        for (int j = 0; j < numParam; j++){
            mu(j) = 0.0;
            for (int i = burninLength; i < this->size(); i++){
              if(param == false){
                mu(j) += theChain[i].getQ(j);
              }
              else{
                mu(j) += theChain[i].getTheta(j);
              }
            }
            mu(j) /= (double) (this->size() - burninLength);
        }
        return mu;
      }

      Eigen::VectorXd inline variance(bool param = false){
        int numParam = theChain[0].size(param);
        Eigen::VectorXd var(numParam);
        auto mu  = this->mean(param);
        for (int j = 0; j < numParam; j++){
            var(j) = 0.0;
            for (int i = burninLength; i < this->size(); i++){
              if(param == false){
                var(j) += std::pow(theChain[i].getQ(j) - mu,2);
              }
              else{
                var(j) += std::pow(theChain[i].getTheta(j) - mu,2);
              }
            }
            var(j) /= (double) (this->size() - burninLength - 1);
        }
        return var;
      }

      Eigen::VectorXd inline samplingError(bool param = true){
        Eigen::VectorXd es(0);
        auto var = this->variance(param);
        auto Neff = this->EffectiveSampleSizes(param);
        int numParam = var.size(param);
        es.resize(numParam);
        for (int j = 0; j < numParam; j++){
              es(j) = var(j) / Neff(j);
        }
        return es;
      }

     

      Eigen::VectorXd EffectiveSampleSizes(bool param = false){

        int num = 10000000;

        int numSamplesUsed = std::min(this->size(),num);

        if(param == false){

          Eigen::VectorXd ESS(theChain[0].getNumQoI());
          for (int j = 0; j < theChain[0].getNumQoI(); j++){
            std::vector<double> vals(numSamplesUsed);
            for (int i = this->size() - numSamplesUsed; i < this->size(); i++){

              auto tmp = theChain[i].getQ();
              vals[i] = tmp[j];
            }
            ESS(j) = getESS(vals);
          }

          return ESS;

        }
        else{

        int numParam = theChain[0].size(); // Number of Parameters
        Eigen::VectorXd ESS(numParam);

        for (int j = 0; j < numParam; j++){
          std::vector<double> vals(numSamplesUsed);
          for (int i = this->size() - numSamplesUsed; i < this->size(); i++){
            vals[i] = theChain[i].getTheta(j);
          }
          ESS(j) = getESS(vals);
        }

        return ESS;

        }

      }

      Eigen::VectorXd getESS_All(){ return this->EffectiveSampleSizes(); }

      double getMaxESS(bool param = false){
        Eigen::VectorXd ESS = this->EffectiveSampleSizes(param);
        double max = ESS.maxCoeff();
        return max;
      }

      double getMinESS(bool param = false){
        Eigen::VectorXd ESS = this->EffectiveSampleSizes(param);
        return ESS.minCoeff();
      }

      double acceptRatio(int lastNSamples = -1){
        if (lastNSamples < 0){lastNSamples = this->size();}

        int numAccept = 0;
         for (int i = this->size() - lastNSamples; i < this->size(); i++){
            numAccept += theChain[i].getAccept();
         }
         double ratio = (double) numAccept / lastNSamples;
         return ratio;
      }

      void setBurninLength(int N){burninLength = N;}

      void setbestObserved(Link& best_){ bestObserved = best_;}

      Link getbestObserved(){ return bestObserved; }

      void readChainFromFile(std::string fileName){

        ifstream myfile(fileName);

        std::string line;

        std::getline(myfile,line); // First line is level

        chainLevel = std::stoi( line );

        std::cout << "This is a level = " << chainLevel << " chain" << std::endl;

        std::getline(myfile,line); 
        int Stochastic_DIM = std::stoi( line );

        std::cout << "Stochastic_DIM = " << line << std::endl;

        std::getline(myfile,line); 
        int numQoI = std::stoi( line );

        std::getline(myfile,line);
        int N = std::stoi(line);

        std::getline(myfile,line);
        int burninLength = std::stoi(line);

        for (int i = burninLength; i < N; i++){  // for each sample

          // Initiate Link

          std::getline(myfile,line);
          std::cout << "Reading --> " << line << std::endl;

          std::getline(myfile,line); // Accept / Reject

          int acceptSample = std::stoi(line);

          std::getline(myfile,line); // LogLikelihood

          double like = std::stod(line);

          std::getline(myfile,line); // LogPrior

          double prior = 0.0; // std::stod(line);


          Eigen::VectorXd QoI_vals(numQoI);
          for (int j = 0; j < numQoI; j++){
            std::getline(myfile,line); // QoI j
            QoI_vals[j] = std::stod(line);
          }


          Eigen::VectorXd vals(N);

          for (int j = 0; j < Stochastic_DIM; j++){
            std::getline(myfile,line);
            vals(j) = std::stod(line);
          }

          Link newLink(vals);

          newLink.setlogPi0(prior);
          newLink.setlogPhi(like);

          addLink(newLink,acceptSample);
      
        }



        myfile.close();

      }

      void writeQoI(std::string filename, int level){

        ofstream myfile;
        myfile.open(filename+".txt");

        for (int i = burninLength; i < theChain.size(); i++){
          for (int j = 0; j < theChain[0].getNumQoI(); j++){
            myfile << theChain[i].getQ() << "\t";
          }
          myfile  << "\n";
        }

      }

      void write2File(std::string fileName, int level){

        ofstream myfile;
        myfile.open (fileName+".txt");
       
        myfile << level << "\n";
        myfile << theChain[0].getSDIM() << "\n";
        myfile << theChain[0].getNumQoI() << "\n";
        myfile << theChain.size() << "\n";
        myfile << burninLength << "\n";



        //  Write all samples to file
        for (int i = burninLength; i < theChain.size(); i++){

          int cnt = i - burninLength;
          myfile << "#Sample " << cnt << "\n"; 
          myfile << theChain[cnt].getAccept() << "\n";
          myfile << theChain[cnt].getlogPhi() << "\n";
          myfile << theChain[cnt].getlogPi0() << "\n";
          myfile << theChain[cnt].getQ() << "\n";
          myfile << theChain[cnt].getTheta() << "\n";
        }

         myfile.close();



      }



  private:

    std::vector<Link> theChain;

    int chainLevel;

    int burninLength;

    Link bestObserved;

    int bestObserved_id;


  };

}
}
#endif /* chain_h */
