#include "iostream.h"
#include "string.h"
#include "memory.h"
#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "time.h"
#include "em.h"

namespace qi
{
    #define PI 3.1415926
    #define CONST_E 0.000000000001
    
    cEM::cEM(int nVec, int nDim, int nClu)
    {
        _nVec = nVec;
        _nDim = nDim;
        _nClu = nClu;

        //_ppGamma[i][j] : Probability of the ith vector belongs to the jth cluster, initialized to 0
        _ppGamma = new double *[_nVec];
        for (int i = 0; i < _nVec; i++)
        {
            _ppGamma[i] = new double[_nClu];
            memset(_ppGamma[i], 0, sizeof(double) * _nClu);
        }
        
        //_ppMiu[i]: The mean vector of the ith cluster
        _ppMiu = new double *[_nClu];
        for (i = 0; i < _nClu; i++)
        {
            _ppMiu[i] = new double[_nDim];
        }
        
        //_pppSigma[i]: The covariance matrix of the ith cluster, initialized to unit matirx
        _pppSigma = new double **[_nClu];
        for (i = 0; i < _nClu; i++)
        {
            _pppSigma[i] = new double *[_nDim];
            for (int j = 0; j < _nDim; j++)
            {
                _pppSigma[i][j] = new double[_nDim];
                memset(_pppSigma[i][j], 0, sizeof(double) * _nDim);
                _pppSigma[i][j][j] = 1;
            }
        }
        
        //_pPi[i]: The probability a priori of the ith cluster, initialized to 1/_nClu
        _pPi = new double [_nClu];
        for (i = 0; i < _nClu; i++) _pPi[i] = 1.0 / _nClu;
    }

    cEM::~cEM()
    {
        for (int i = 0; i < _nVec; i++) delete []_ppGamma[i];
        delete []_ppGamma;

        for (i = 0; i < _nClu; i++) delete[]_ppMiu[i];
        delete []_ppMiu;
        
        for (i = 0; i < _nClu; i++)
        {
            for (int j = 0; j < _nDim; j++) delete[]_pppSigma[i][j];
            delete []_pppSigma[i];
        }
        delete []_pppSigma;
        delete []_pPi;
    }

    double **cEM::Inverse(double **ppMatSrc)
    {
        // Copy the original matrix
        double **ppMat = new double *[_nDim];
        for (int i = 0; i < _nDim; i++)
        {
            ppMat[i] = new double [_nDim];
            memcpy(ppMat[i], ppMatSrc[i], sizeof(double) * _nDim);
        }
        
        // Create a unit matrix
        double **ppI = new double *[_nDim];
        for (i = 0; i < _nDim; i++)
        {
            ppI[i] = new double [_nDim];
            memset(ppI[i], 0, sizeof(double) * _nDim);
            ppI[i][i] = 1;
        }
        
        // Matrix inverse
	    // Gaussian elimination: Forward Elimination
        for (i = 0; i < _nDim; i++)
        {
            double Tmp = ppMat[i][i];
            for (int j = 0; j < _nDim; j++)
            {
                ppMat[i][j] /= Tmp;
                ppI[i][j] /= Tmp;
            }
            
            for (j = i + 1; j < _nDim; j++)
            {
                double Tmp = -ppMat[j][i];
                for (int k = 0; k < _nDim; k++)
                {
                    ppMat[j][k] += Tmp * ppMat[i][k];
                    ppI[j][k] += Tmp * ppI[i][k];
                }
            }
		}
        // Gaussian elimination: Back Substitution
        for (i = _nDim - 1; i >= 0; i--)
            for (int j = i - 1; j >= 0; j--)
            {
                double Tmp = -ppMat[j][i];
                for (int k = 0; k < _nDim; k++)
                {
                    ppMat[j][k] += ppMat[i][k] * Tmp;
                    ppI[j][k] += ppI[i][k] * Tmp;
                }
            }

		// Release memory block
        for (i = 0; i < _nDim; i++) delete []ppMat[i];
        delete []ppMat;
        
        return ppI;
    }
    
    double cEM::Determinant(double **ppMatSrc)
    {
		// Copy the original matrix
        double **ppMat = new double *[_nDim];
        for (int i = 0; i < _nDim; i++)
        {
            ppMat[i] = new double [_nDim];
            memcpy(ppMat[i], ppMatSrc[i], sizeof(double) * _nDim);
        }
        
        // Gaussian elimination to the diagonal matrix
        for (i = 0; i < _nDim; i++)
        {
            for (int j = i + 1; j < _nDim; j++)
            {
                double Tmp = -ppMat[j][i] / ppMat[i][i];
                for (int k = 0; k < _nDim; k++) ppMat[j][k] += Tmp * ppMat[i][k];
            }
        }
        
		// Calculate the multiplication of elements on the diagonal
        double Rs = 1;
        for (i = 0; i < _nDim; i++) Rs *= ppMat[i][i];
            
		// Release memory block
        for (i = 0; i < _nDim; i++) delete []ppMat[i];
        delete []ppMat;
            
        return Rs;
    }
    
    double cEM::Gaussian(double *pVec, double *pMiu, double **ppSigma)
    {
        // Find the PDF of Multivariate Gaussian distribution
        // N(x|Miu,Sigma) = 1/((2PI)^(nDim/2))*(1/(abs(Sigma))^0.5)*exp(-1/2*(x-Miu)'Sigma^(-1)*(x-Miu))  
        double *pTmp1 = new double[_nDim];
        double *pTmp2 = new double[_nDim];
        
        for (int i = 0; i < _nDim; i++) pTmp1[i] = pVec[i] - pMiu[i];
        memset(pTmp2, 0, sizeof(double) * _nDim);
        
        double **ppInvSigma = Inverse(ppSigma);

        for (i = 0; i < _nDim; i++)
            for (int j = 0; j < _nDim; j++)
                pTmp2[i] += pTmp1[j] * ppInvSigma[j][i];
                
        double Tmp1 = 0;
        for (i = 0; i < _nDim; i++)
            Tmp1 += pTmp2[i] * pTmp1[i];
        Tmp1 /= -2;
        Tmp1 = exp(Tmp1);
        
        double Sigma = Determinant(ppSigma);
        double Tmp2 = 1 / sqrt (pow(2 * PI, _nDim) * Sigma) * Tmp1;
        
        for (i = 0; i < _nDim; i++) delete []ppInvSigma[i];
        delete []ppInvSigma;
        delete []pTmp1;
        delete []pTmp2;

        if (Tmp2 < CONST_E) Tmp2 = CONST_E;
        return Tmp2;
    }

    void cEM::Expectation()
    {
        //E-Step
        _LL = 0;
        for (int i = 0; i < _nVec; i++)
        {
            double Sum = 0;// For calculating the log likelihood
            for (int j = 0; j < _nClu; j++)
            {
                double Tmp = 0;
                for (int l = 0; l < _nClu; l++)
                {
                    double G = Gaussian(_ppVector[i], _ppMiu[l], _pppSigma[l]);
					//Gaussian posterior probability   
                    Tmp += _pPi[l] * G;
                }
                double G = Gaussian(_ppVector[i], _ppMiu[j], _pppSigma[j]);
                _ppGamma[i][j] = _pPi[j] * G / Tmp; 
				// numerator: pi(k) * N(xi | Miu(k), Sigma(k)) 
				// denominator: Sumj ( pi(j) * N(xi | Miu(j), Sigma(j)) )
                Sum += _pPi[j] * G;
            }
            _LL += log(Sum);
        }
    }
    
    void cEM::Maximization()
    {
        // M-Step
        for (int j = 0; j < _nClu; j++)
        {
            // The sum of probability that the jth Gaussian generated each vector
			double Tmp1 = 0;
            for (int i = 0; i < _nVec; i++) Tmp1 += _ppGamma[i][j];
            
            // Updata the mean vector matrix _ppMiu[j]
			// through MLE, i.e. let the derivative equals to zero
            memset(_ppMiu[j], 0, sizeof(double) * _nDim);
            for (i = 0; i < _nVec; i++)
                for (int k = 0; k < _nDim; k++)
                    _ppMiu[j][k] += _ppGamma[i][j] * _ppVector[i][k];
            for (int k = 0; k < _nDim; k++) _ppMiu[j][k] /= Tmp1;

            // Update the covariance matrix _pppSigma[j]
            for (i = 0; i < _nDim; i++)
                memset(_pppSigma[j][i], 0, sizeof(double) * _nDim);
            for (i = 0; i < _nVec; i++)
            {
                double *pTmp = new double[_nDim];
                for (int a = 0; a < _nDim; a++) pTmp[a] = _ppVector[i][a] - _ppMiu[j][a];

                // Only consider the individual variances
                for (a = 0; a < _nDim; a++)
                    _pppSigma[j][a][a] += _ppGamma[i][j] * pTmp[a] * pTmp[a];
                delete[]pTmp;
            }
            for (int a = 0; a < _nDim; a++)
            {
                _pppSigma[j][a][a] = _pppSigma[j][a][a] / Tmp1;
                _pppSigma[j][a][a] += CONST_E;
            }
                    
            // Update the weight of each Gaussian
            _pPi[j] = Tmp1 / _nVec;
        }
    }
    
    void cEM::GetResult(double ** &ppGamma, double ** &ppMiu, double *** &pppSigma)
    {
        ppGamma = _ppGamma;
        ppMiu = _ppMiu;
        pppSigma = _pppSigma;
    }
    
    void cEM::InitCluster()
	// _ppMiu[i][j] is randomly initialized with a value between pMin[j] and pMax[j]
	// where pMin[j] and pMax[j] are the minimum and maximum value of the jth feature 
    {
        double *pMax = new double[_nDim];
        double *pMin = new double[_nDim];
        
        memcpy(pMax, _ppVector[0], sizeof(double) * _nDim);
        memcpy(pMin, _ppVector[0], sizeof(double) * _nDim);
        
        for (int i = 1; i < _nVec; i++)
            for (int j = 0; j < _nDim; j++)
            {
                if (_ppVector[i][j] < pMin[j]) pMin[j] = _ppVector[i][j];
                if (_ppVector[i][j] > pMax[j]) pMax[j] = _ppVector[i][j];
            }

        srand((unsigned)time(NULL));
        rand();
        
        for (i = 0; i < _nClu; i++)
        {
            for (int j = 0; j < _nDim; j++)
            {
                _ppMiu[i][j] = (double)rand() / RAND_MAX * (pMax[j] - pMin[j]) + pMin[j];
            }
        }
 
        delete []pMax;
        delete []pMin;
    }
        
    void cEM::Cluster(double **ppVector)
    {
        _ppVector = ppVector;
        
        // Initialize the cluster centers
        InitCluster();

        // EM algorithm
        
        Expectation();
        Maximization();
            
        int nCount = 0;
        double LL_bf; // Log likelihood before
        do
        {
            LL_bf = _LL;
            nCount++;
            Expectation();
            Maximization();
        } while(fabs(_LL - LL_bf) > CONST_E); // Check for convergence
    }
}
