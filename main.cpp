#include "stdio.h"
#include "iostream.h"
#include "fstream.h"
#include "string.h"
#include "algorithm"
#include "em.h"

using namespace qi;

int main()
{	
	int DataNum;
	for (DataNum = 1; DataNum <= 8; DataNum++)
	{
		int nVec, nDim, nClu;
		char DataText[20];
		char ResultText[20];
		char PlotText[20];

		// Feature data readin
		sprintf(DataText, "./InputData/data%d.txt",DataNum);
		FILE *fp = fopen(DataText, "r");
		fscanf(fp, "%d%d%d", &nVec, &nDim, &nClu);
    
		double **a = new double*[nVec];
		for (int i = 0; i < nVec; i++)
		{
			a[i] = new double[nDim];
			for (int j = 0; j < nDim; j++) fscanf(fp, "%lf", &a[i][j]);
		}
		fclose(fp);

		// Result data output
		sprintf(ResultText, "./OutputData/Result%d.txt",DataNum);
		sprintf(PlotText, "./OutputData/Plot%d.txt",DataNum);

		FILE *fpr= fopen(ResultText, "w");
		fprintf(fpr, "\nNumber of vectors: %d",nVec);
		fprintf(fpr, "\nNumber of dimensions: %d",nDim);
		fprintf(fpr, "\nNumber of cluster: %d\n",nClu);
		
		// EM algorithm for clustering
		cEM em(nVec, nDim, nClu);
		em.Cluster(a);	    
		double **ppGamma, **ppMiu, ***pppSigma;
		em.GetResult(ppGamma, ppMiu, pppSigma);

		// For each cluster, find the players belongs to it
		// For the purpose of comparing 8 results, the clusters are ordered according to # of players

		int k=0; // For output elements in each cluster
		int n=0; // For output format
		int ** Result; // Result matrix, [nClu, nEle]
		int * EleNum; // # of elements in each cluster
		int * TempEleNum; // For sorting clusters
		int * ClusterIndex; // Cluster index, elements number descent

		EleNum = new int [nClu];
		TempEleNum = new int [nClu];
		Result = new int *[nClu];
		ClusterIndex = new int [nClu];

		for (i = 0; i <nClu; i++)
		{
			Result[i] = new int [nVec];
			memset(Result[i], 0, sizeof(int) * nVec);
		}
		for (int j = 0; j < nClu; j++)
		{
			k = 0;
			for (i = 0; i < nVec; i++)
			{
				if (ppGamma[i][j]>0.95) // 95% confidence				
				{
					Result[j][k]=i+1;
					k++; // Counting # of players belong to the jth cluster
				}
			}
			EleNum[j] = k;
		}

		for (i=0;i<nClu;i++)
			TempEleNum[i] = EleNum[i];
		
		std::sort(TempEleNum, TempEleNum + nClu); // Ordering the cluster: # of players ascent	
		std::reverse(TempEleNum, TempEleNum + nClu); // Ordering the cluster: # of players descent

		for (i=0 ; i<nClu; i++)
		{
			for (j=0; j<nClu; j++)
			{
				if (TempEleNum[i] == EleNum[j])
				{
					ClusterIndex[i] = j; // Record the sorted order of clusters
					break;
				}
			}
		}	
		
		// Miu matrix for each cluster
		for (i = 0; i < nClu; i++)
		{
			fprintf(fpr, "Cluster %d's mean Matrix\n\n", i+1);
			for (int j = 0; j < nDim; j++)			
				fprintf(fpr, "%e\t\t",ppMiu[ClusterIndex[i]][j]);		
			fprintf(fpr,"\n\n");
		}

		// Sigma matrix for each cluster
		for (i = 0; i < nClu; i++)
		{
			fprintf(fpr, "Cluster %d's covariance Matrix\n\n", i+1);
			for (int j = 0; j < nDim; j++)
			{
				for (int a = 0; a < nDim; a++)
					fprintf(fpr,"%e\t\t",pppSigma[ClusterIndex[i]][j][a]);		
				fprintf(fpr,"\n");
			}
			fprintf(fpr,"\n");
		}
		
		// Clustering result for 300 players 
		for (i=0 ; i<nClu ; i++)
		{
			n=0;
			fprintf (fpr, "\nCluster %d\n",i+1);
			for (j=0 ; j<nVec ;j++)
			{
				if (Result[ClusterIndex[i]][j]!=0)
				{
					if (n%10 == 0)fprintf(fpr, "\n");
					fprintf (fpr, "%d\t",Result[ClusterIndex[i]][j]);
					n++;
				}			
			}
			fprintf (fpr, "\n");
		}
		fclose (fpr);
		
		// For matlab plotting
		FILE *fpp= fopen(PlotText, "w");
		for (i = 0; i < nClu; i++)
		{
			for (int j = 0; j < nVec; j++) 
			{
				if (Result[ClusterIndex[i]][j]!=0)
				{
					fprintf (fpp, "%d\t",Result[ClusterIndex[i]][j]);
				}	
			}
			fprintf(fpp,"\n");
		}
		fclose(fpp);
		for (i = 0; i < nVec; i++) delete []a[i];
		delete[]a;
		
	}
	printf("\nClustering Finished!\n");
	return 0;
}
