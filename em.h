#ifndef EM_H
#define EM_H

namespace qi
{
    class cEM
    {
    private:
        int _nVec;     // Number of the vectors to be clustered
        int _nDim;     // Number of dimensions
        int _nClu;     // Number of clusters
        
        double **_ppVector;   // Vectors to be clustered
        double **_ppGamma;        // Probability of the ith vector belongs to the jth pattern
        double **_ppMiu;      // The mean vector of the ith pattern
        double *_pPi;         // The probability a priori of the ith pattern
        double ***_pppSigma;  // The covariance matrix of the ith pattern
        
        double _LL;           // For storing the log likelihood value, calculated in the E-Step
                              // Objective: Maximize _LL
        
    public:

        //Initialization
        cEM(int nVec, int nDim, int nClu);
        
        //Release memory block
        virtual ~cEM();

        //Clustering
        void Cluster(double **ppVector);
        
        //Get the result
        void GetResult(double ** &ppGamma, double ** &ppMiu, double *** &pppSigma);
        
    private:

        //Find clustering centers
        void InitCluster();
        
        //E-Step
        void Expectation();

        //M-Step
        void Maximization();

        //Solve determinant of matrix ppMat
        double Determinant(double **ppMat);

        //Find the inverse of matrix ppMat
        double ** Inverse(double **ppMat);
        
        //Find the PDF of Multivariate Gaussian distribution
        double Gaussian(double *pVec, double *pMiu, double **ppSigma);
    };
}

#endif
