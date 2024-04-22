The HeteroAnaysis program will automattically handle all of the histogram fitting (and fit plots) for your data. 
Add this entire folder to your path and run HeteroAnalysis. 
It will prompt you to select the root directory, and will handle all processing from there. 
If you want to run data individually, you can call the functons from the command line. 
weightedHistogram will create your histogram, and freqCurve will fit it with a GMM. To plot this outcome, use freqCurvePlot. 
To adjust the criterion you wish to use, open and edit the freqCurve function. Multiple options are available in the comments.