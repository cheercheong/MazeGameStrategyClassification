# MazeGameStrategyClassification
Implementation of Expectation–maximization (EM) algorithm in C++ to classify the strategies applied by 300 plays when playing a maze game

//**************************************************************************//

Files field1.txt through field8.txt each contain the “topology” of ﬁeld 1 through 8. That is, each
ﬁle, a text ﬁle, is a 512-by-512 matrix, where the matrix entry is the diﬃculty level of the ﬁeld at the
corresponding location.

Files playOnField1.txt to playOnField8.txt each contain the trajectories of the 300 players when
they played on respectively ﬁeld 1 through 8. The trajectory of a player is the row-column coordinate of
the player at each time instant. 

//**************************************************************************//

Matlab part: DataPreprocessing

Src: data_preprocessing.m , for data preprocessing and feature extraction, Observation.m , for checking the result

InputData: 8 field map text files and 8 corresponding playonfield text files.

OutputData: 8 feature data files, have been put into the StrategyClustering_cpp folder.

Warning: running data_preprocessing takes a lot of time, you can set the k value range at the beginning from 1 to 1. 

//**************************************************************************//

C++ part: StrategyClustering

Src: Source code folder, including main.cpp, em.cpp, em.h

InputData: 8 data text files generated by matlab, containing the feature extracted

OutputData: 16 data text files, each field has two data text files, for example, for field1, the resulting files
	     are Result1.txt and Plot1.txt, the first one is for checking the clustering result, and the second one
	     is for matlab plotting.

Open the StrategyClustering.dsw, compile and run.

//**************************************************************************//

