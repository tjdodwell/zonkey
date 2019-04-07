

class ParallelStats
{

public:

	double Bias, VoE, expectedValue;
	int greedyLevel;

	ParallelStats(){};

	ParallelStats(double Bias_, double VoE_, int greedyLevel_, double expVal_) :
		Bias(Bias_),
		VoE(VoE_),
		greedyLevel(greedyLevel_),
		expectedValue(expVal_){
		}

	void inline setEQ(int l, double value){	EQs.push_back(value);	}

	double inline getEQ(int l){	return EQs[l];	}

	void inline setVQ(int l, double value){ VQs.push_back(value);	}

	double inline getVQ(int l){ return VQs[l];	}

private:

	std::vector<double> EQs, VQs;



};



template<class Samples>
ParallelStats inline parallel_Stats_Update(Samples& Q,int rank, int L, int nproc, double alpha_rate){

	std::vector<double> localEQ(L + 1), localEQ2(L + 1);

	Q.getEQs(localEQ);  Q.getEQ2s(localEQ2);


	double * EQRoot = NULL;
    double * EQ2Root = NULL;

    if (rank == 0){
        EQRoot = new double[nproc];
        EQ2Root = new double[nproc];
    }

    double Bias, Var;

    std::vector<double> PayOff(L + 1), VarEstim(L + 1);

    double sum = 0.0;

    std::vector<double> levelEQ(L+1), levelVQ(L+1);

    for (int l = 0; l < L + 1; l++){ // For each level

          // Communicate EQ and EQ2 from each processor to rank 0 processor
          MPI_Gather(&localEQ[l],1,MPI_DOUBLE,EQRoot,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
          MPI_Gather(&localEQ2[l],1,MPI_DOUBLE,EQ2Root,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

          if (rank == 0){

          	// * First compute EQ and E[Q^2] - from those collected from each processor
            double EQ_all = 0.0;
            double EQ2_all = 0.0;

            for (int i = 0; i < nproc; i++){
            	EQ_all += EQRoot[i];
            	EQ2_all += EQ2Root[i];
            }
            EQ_all /= nproc;
            EQ2_all /= nproc;

            sum += EQ_all;

            levelEQ[l] = EQ_all;

            Var = EQ2_all - EQ_all * EQ_all;

            levelVQ[l] = Var;

            if(l == L){	Bias = std::abs(EQ_all / std::pow(4.0,alpha_rate));  }

            VarEstim[l] =  Var / (Q.numSamples(l) * nproc); // Variance of Estimate

            PayOff[l] = VarEstim[l] / Q.ExpectedCost(l);

          }

    }

    double VoE = 0.0;

    int greedyLevel = 0;

    if (rank == 0){

    	VoE = std::accumulate(VarEstim.begin(), VarEstim.end(),0.0);

    	// Find level with maximum payoff

        std::vector<double>::iterator maxPayOffLevel = std::max_element(PayOff.begin(),PayOff.end());

        greedyLevel = std::distance(PayOff.begin(),maxPayOffLevel);



    }


    MPI_Bcast(&Bias,1,MPI_DOUBLE,0,MPI_COMM_WORLD); // BroadCast Bias
    MPI_Bcast(&VoE,1,MPI_DOUBLE,0,MPI_COMM_WORLD); // BroadCast VoE
    MPI_Bcast(&greedyLevel,1,MPI_INT,0,MPI_COMM_WORLD);

    ParallelStats myStats(Bias,VoE,greedyLevel,sum);

    for (int i = 0; i < L+1;i++){
    	myStats.setEQ(i,levelEQ[i]);
    	myStats.setVQ(i,levelVQ[i]);
    }

    return myStats;

}
