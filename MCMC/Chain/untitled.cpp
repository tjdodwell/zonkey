



// SubSample Batch


int subSamplingRate = 100;

int batch = 5000; // Not this is hard coded at the moment in write file

int totalFile = (coarseSamples / batch) - 1;

CHAIN coarseChain;

std::string fileName = "chainDataOutput";

int count = 0; // Global Count of Samples independent of file.

for (int i = 0; i < totalFile; i++){

	ifstream myfile(fileName + to_string(i) + ".txt");

	for (int j = 0; j < 5; j++){
		std::getline(myfile,line);
	}


	for (int j = 0; j < batch; j++){

		Link newLink readSample(myfile,line);

		if (count % subSamplingRate == 0){

			coarseChain.addLink(newLink);
			
		}
	}

	myfile.close();

}















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

        for (int i = 0; i < N; i++){  // for each sample

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


          if (i >= burninLength){
            addLink(newLink,acceptSample);
          }
      
        }



        myfile.close();

      }


